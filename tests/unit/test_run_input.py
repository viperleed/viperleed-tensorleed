"""Unit tests of the input a ref-calc run is fed: converted FIN or original, by tree."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest

from harness import cases, perf
from harness.fin_convert import converted_values, reads_converted_fin
from harness.param_size import explicit_values, read_param

CASE, RUN = next((c, r) for c, r in cases.iter_runs(("refcalc",), only=("ir100_2x1_o",))
                 if r.name == "E050_L06")
BENCHMARK_CASE = next(c for c in cases.load_cases_from(cases.CASES_DIR, "benchmark")
                      if c.program == "refcalc")
# a stub program that stores its standard input in the work directory
KEEP_STDIN = """
    import sys
    open("seen", "w").write(sys.stdin.read())
"""


def _tree(root, ported):
    """A TL tree under root; a ported one has the new reader."""
    lib = root / "TensErLEED" / "lib"
    lib.mkdir(parents=True)
    if ported:
        (lib / "leed_input.f90").write_text("module leed_input\nend module leed_input\n")
    return root


@pytest.fixture
def no_completeness_check(mocker):
    """The stub programs write no fd.out."""
    mocker.patch.object(cases, "check_complete")


def test_tree_with_the_ported_reader_reads_the_converted_fin(tmp_path):
    """Only a tree that has lib/leed_input.f90 reads the converted FIN."""
    assert reads_converted_fin(_tree(tmp_path / "ported", True))
    assert not reads_converted_fin(_tree(tmp_path / "legacy", False))


def test_run_stdin_is_the_original_for_a_legacy_tree(tmp_path):
    """A tree without the ported reader is fed the run's FIN as it is."""
    legacy = _tree(tmp_path / "legacy", False)
    assert cases.run_stdin(CASE, RUN, legacy, tmp_path) == RUN.stdin


def test_run_stdin_is_a_converted_copy_for_a_ported_tree(tmp_path):
    """A tree with the ported reader is fed a converted copy, with the values of the run's PARAM."""
    ported = _tree(tmp_path / "ported", True)
    path = cases.run_stdin(CASE, RUN, ported, tmp_path)
    assert path != RUN.stdin and path.parent == tmp_path
    assert converted_values(path.read_text()) == explicit_values(read_param(RUN.param))
    assert RUN.stdin.read_text() != path.read_text()


def test_run_stdin_of_rfactor_is_never_converted(tmp_path):
    """The rfactor input has no converted form."""
    case = cases.Case("c", tmp_path, "rfactor", [])
    assert cases.run_stdin(case, RUN, _tree(tmp_path / "ported", True), tmp_path) == RUN.stdin


def test_perf_runs_give_each_tree_its_own_input(tmp_path, stub_script, stub_build, mocker,
                                                no_completeness_check):
    """compare_performance feeds the baseline tree the original FIN, the candidate a converted one."""
    stub_build(stub_script(tmp_path / "keep.py", KEEP_STDIN))
    seen = {}
    real = cases.run_stdin

    def spy(case, run, tl_root, workdir):
        """Record the text the real function chooses for the tree."""
        path = real(case, run, tl_root, workdir)
        seen.setdefault(Path(tl_root).name, set()).add(Path(path).read_text())
        return path
    mocker.patch.object(cases, "run_stdin", spy)
    perf.compare_performance(CASE, RUN, _tree(tmp_path / "cand", True),
                             _tree(tmp_path / "base", False), reps=1,
                             workroot=tmp_path / "w")
    assert seen["base"] == {RUN.stdin.read_text()}
    (cand_text,) = seen["cand"]
    assert converted_values(cand_text) == explicit_values(read_param(RUN.param))


def test_benchmark_candidate_gets_a_converted_fin(tmp_path, stub_script, stub_build,
                                                  no_completeness_check):
    """execute_run on a ported tree feeds a benchmark run its converted FIN."""
    run = BENCHMARK_CASE.runs[0]
    stub_build(stub_script(tmp_path / "keep.py", KEEP_STDIN))
    cases.execute_run(BENCHMARK_CASE, run, _tree(tmp_path / "cand", True), None, tmp_path / "wd")
    seen = (tmp_path / "wd" / "seen").read_text()
    assert converted_values(seen) == explicit_values(read_param(run.param))
    cases.execute_run(BENCHMARK_CASE, run, _tree(tmp_path / "base", False), None, tmp_path / "wb")
    assert (tmp_path / "wb" / "seen").read_text() == run.stdin.read_text()
