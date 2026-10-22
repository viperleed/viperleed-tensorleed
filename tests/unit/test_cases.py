"""Tests of case loading, run execution and completeness checks of harness.cases."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import json
import tomllib

import pytest

from harness import cases
from harness.cases import (CASES_DIR, CaseError, expected_outputs, load_cases,
                           load_cases_from, tolerance_for)
from harness.runner import ProgramFailed

CASE_TOML = """
program = "refcalc"
description = "synthetic"
covers = ["x"]
[[runs]]
name = "r1"
stdin = "inputs/FIN"
param = "inputs/PARAM"
outputs = ["fd.out", "T_*"]
"""


@pytest.fixture
def case_root(tmp_path):
    """Return a root holding one synthetic refcalc case with one run."""
    c = tmp_path / "refcalc" / "c1"
    (c / "inputs").mkdir(parents=True)
    (c / "expected" / "r1").mkdir(parents=True)
    (c / "case.toml").write_text(CASE_TOML)
    for f in ("inputs/FIN", "inputs/PARAM", "expected/r1/fd.out",
              "expected/r1/T_2", "expected/r1/T_1"):
        (c / f).write_text("1.0")
    return tmp_path


def test_load_and_resolve(case_root):
    """Loading a case resolves its run paths and lists expected outputs in sorted order."""
    (case,) = load_cases_from(case_root, "refcalc")
    assert case.name == "c1" and case.runs[0].stdin.name == "FIN"
    assert case.runs[0].param.name == "PARAM" and case.runs[0].files == {}
    assert [p.name for p in expected_outputs(case, case.runs[0])] == ["T_1", "T_2", "fd.out"]


@pytest.mark.parametrize("make_dir", [False, True], ids=["missing", "empty"])
def test_missing_or_empty_kind_dir_raises(tmp_path, make_dir):
    """A missing or empty case directory is an error, not zero tests."""
    if make_dir:
        (tmp_path / "benchmark").mkdir()
    with pytest.raises(CaseError, match="no cases"):
        load_cases_from(tmp_path, "benchmark")


def test_empty_glob_is_error(case_root):
    """An output glob that matches no expected file raises instead of comparing nothing."""
    (case_root / "refcalc/c1/expected/r1/T_1").unlink()
    (case_root / "refcalc/c1/expected/r1/T_2").unlink()
    (case,) = load_cases_from(case_root, "refcalc")
    with pytest.raises(FileNotFoundError):
        expected_outputs(case, case.runs[0])


def test_tolerance_lookup():
    """tolerance_for finds a known name and raises KeyError for an unknown one."""
    tol, ignore, skip = tolerance_for("T_3")
    assert tol.rtol > 0
    with pytest.raises(KeyError):
        tolerance_for("unknown.txt")


def test_execute_run_copies_files_and_runs(tmp_path, stub_script, stub_build):
    """execute_run feeds stdin, copies extra files in and records the run info."""
    stub = stub_script(tmp_path / "stub.py", """
        import sys, pathlib
        pathlib.Path('fd.out').write_text(sys.stdin.read() + pathlib.Path('WEXPEL').read_text())
    """)
    stub_build(stub)
    (tmp_path / "in").write_text("A")
    (tmp_path / "w").write_text("B")
    case = cases.Case("c", tmp_path, "stub", [])
    run = cases.Run("r", tmp_path / "in", None, {"WEXPEL": tmp_path / "w"}, ["fd.out"])
    res = cases.execute_run(case, run, tmp_path, None, tmp_path / "wd")
    assert (tmp_path / "wd" / "fd.out").read_text() == "AB" and res.returncode == 0
    info = json.loads((tmp_path / "wd" / ".runinfo.json").read_text())
    assert info["returncode"] == 0 and "max_rss_kb" in info


@pytest.mark.parametrize("listed", [True, False], ids=["listed", "not_listed"])
def test_execute_run_stores_stdout_if_listed(tmp_path, stub_script, stub_build, listed):
    """The standard output becomes the file 'stdout' exactly when run.outputs lists it."""
    stub_build(stub_script(tmp_path / "stub.py", "print('*BEAM NO.  1')"))
    (tmp_path / "in").write_text("")
    case = cases.Case("c", tmp_path, "stub", [])
    run = cases.Run("r", tmp_path / "in", None, {}, ["stdout"] if listed else [])
    cases.execute_run(case, run, tmp_path, None, tmp_path / "wd")
    stored = tmp_path / "wd" / "stdout"
    assert stored.is_file() == listed
    if listed:
        assert stored.read_text() == "*BEAM NO.  1\n"


def test_execute_run_failing_program_raises(tmp_path, stub_script, stub_build):
    """A program that exits non-zero makes execute_run raise ProgramFailed."""
    stub_build(stub_script(tmp_path / "stub.py", "import sys; sys.exit(2)"))
    (tmp_path / "in").write_text("A")
    case = cases.Case("c", tmp_path, "stub", [])
    run = cases.Run("r", tmp_path / "in", None, {}, ["fd.out"])
    with pytest.raises(ProgramFailed):
        cases.execute_run(case, run, tmp_path, None, tmp_path / "wd")


def test_tolerance_lookup_reads_optional_atol(mocker):
    """An absolute tolerance given in the configuration reaches the Tolerance."""
    mocker.patch.object(cases, "load_tolerances", return_value={"numeric": {
        "ROUT": {"rtol": 0.0, "atol_frac": 0.0, "atol": 1.5e-4}}})
    tol, _, _ = cases.tolerance_for("ROUT")
    assert tol.atol == 1.5e-4


FIN_2E = "title\n   138.00   141.01     3.00   EI,EF,DE\n rest\n"
FD_HEAD = "t\n    1\n    2   1.00000   0.00000  1\n"


def test_refcalc_complete_when_all_energies_present(tmp_path):
    """A fd.out holding every energy of the EI,EF,DE range counts as complete."""
    (tmp_path / "FIN").write_text(FIN_2E)
    (tmp_path / "fd.out").write_text(FD_HEAD + " 138.00 0.0001   0.1E-02\n\n 141.00 0.0001   0.2E-02\n")
    cases.check_complete("refcalc", tmp_path / "FIN", tmp_path)


def test_refcalc_incomplete_when_program_stopped_early(tmp_path):
    """A fd.out lacking an energy, or missing altogether, counts as incomplete.

    Legacy error paths print a message and STOP with exit code 0.
    """
    (tmp_path / "FIN").write_text(FIN_2E)
    (tmp_path / "fd.out").write_text(FD_HEAD + " 138.00 0.0001   0.1E-02\n")
    with pytest.raises(ProgramFailed, match="incomplete"):
        cases.check_complete("refcalc", tmp_path / "FIN", tmp_path)
    (tmp_path / "fd.out").unlink()
    with pytest.raises(ProgramFailed, match="incomplete"):
        cases.check_complete("refcalc", tmp_path / "FIN", tmp_path)


def test_rfactor_incomplete_without_summary(tmp_path):
    """A ROUT without the average R-factor summary line counts as incomplete."""
    (tmp_path / "ROUT").write_text("AVERAGE   0  0.0001  3.00\n")
    with pytest.raises(ProgramFailed, match="incomplete"):
        cases.check_complete("rfactor", None, tmp_path)
    (tmp_path / "ROUT").write_text("WITH INNER POTENTIAL SHIFT    3.00 EV, AVERAGE R-FACTOR =  0.2054\n")
    cases.check_complete("rfactor", None, tmp_path)


def test_cases_match_the_manifest():
    """The loaded (case, run) pairs equal the manifest.

    A deleted case or one with no runs must not silently shrink the parametrization.
    """
    manifest = tomllib.loads((CASES_DIR / "MANIFEST.toml").read_text())
    for kind in ("refcalc", "rfactor"):
        loaded = [f"{c.name}/{r.name}" for c in load_cases(kind) for r in c.runs]
        assert loaded, f"no {kind} cases"
        assert sorted(loaded) == sorted(manifest[kind])


def test_case_without_runs_or_case_toml_raises(tmp_path):
    """A case with no runs, or a directory without case.toml, is an error."""
    (tmp_path / "refcalc" / "norun").mkdir(parents=True)
    (tmp_path / "refcalc" / "norun" / "case.toml").write_text('program = "refcalc"\nruns = []\n')
    with pytest.raises(CaseError, match="no runs"):
        load_cases_from(tmp_path, "refcalc")
    (tmp_path / "rfactor" / "nodef").mkdir(parents=True)
    with pytest.raises(CaseError, match="without case.toml"):
        load_cases_from(tmp_path, "rfactor")
