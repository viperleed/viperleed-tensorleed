"""Baseline-tree and coverage guards of the generators (tests/generate).

The viperleed generators are imported with a stub fortranformat module
and a dummy viperleed_examples location, so these tests run in tl-tests.
"""
import importlib
import subprocess
import sys
import types
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "generate"))
import baseline_tree  # noqa: E402
from baseline_tree import GeneratorError, check_baseline_tree  # noqa: E402


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
                           *args], check=True, capture_output=True, text=True).stdout.strip()


@pytest.fixture
def repo(tmp_path):
    """A git tree with a committed TensErLEED/; never the real baseline commit."""
    tree = tmp_path / "tl"
    (tree / "TensErLEED").mkdir(parents=True)
    (tree / "TensErLEED" / "a.f").write_text("      END\n")
    _git(tree, "init", "-q")
    _git(tree, "add", ".")
    _git(tree, "commit", "-qm", "init")
    return tree


@pytest.fixture
def at_baseline(repo, monkeypatch):
    monkeypatch.setattr(baseline_tree, "BASELINE_COMMIT", _git(repo, "rev-parse", "HEAD"))
    return repo


def _dirty(tree):
    (tree / "TensErLEED" / "a.f").write_text("      STOP\n      END\n")


def _generator(name, monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "fortranformat", types.ModuleType("fortranformat"))
    monkeypatch.setenv("TL_TEST_VIPERLEED_EXAMPLES", str(tmp_path))
    monkeypatch.setenv("TL_TEST_BASELINE", str(tmp_path))
    return importlib.import_module(name)


def _must_not_run(what):
    def fail(*args, **kwargs):
        raise AssertionError(f"{what} ran")
    return fail


# ------------------------------------------------------------ the checker ---

def test_wrong_commit_refused(repo):
    with pytest.raises(GeneratorError, match="not the baseline"):
        check_baseline_tree(repo)


def test_dirty_tree_refused(at_baseline):
    _dirty(at_baseline)
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(at_baseline)


def test_dirty_tree_refused_with_other_commit_allowed(repo):
    _dirty(repo)
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(repo, allow_other_commit=True)


def test_untracked_file_refused(at_baseline):
    (at_baseline / "TensErLEED" / "b.f").write_text("      END\n")
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(at_baseline)


def test_ignored_file_refused(at_baseline):
    (at_baseline / ".git" / "info" / "exclude").write_text("*.f90\n")
    (at_baseline / "TensErLEED" / "b.f90").write_text("end\n")
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(at_baseline)


def test_ignored_subdirectory_of_baseline_refused(at_baseline):
    """An ignored copy of TensErLEED/ inside a checkout at the baseline
    must not pass as that checkout."""
    sub = at_baseline / "foo"
    (sub / "TensErLEED").mkdir(parents=True)
    (sub / ".gitignore").write_text("*\n")
    (sub / "TensErLEED" / "a.f").write_text("      STOP\n      END\n")
    with pytest.raises(GeneratorError, match="not the top level"):
        check_baseline_tree(sub)


def test_clean_baseline_accepted(at_baseline):
    assert check_baseline_tree(at_baseline) == _git(at_baseline, "rev-parse", "HEAD")


def test_expected_outputs_come_from_the_pinned_baseline():
    """cases/BASELINE.txt names the commit the generators are pinned to."""
    info = (Path(__file__).parents[1] / "cases" / "BASELINE.txt").read_text()
    (commit,) = [l.split(":", 1)[1].strip() for l in info.splitlines()
                 if l.startswith("commit:")]
    assert commit.startswith(baseline_tree.BASELINE_COMMIT), (commit, baseline_tree.BASELINE_COMMIT)


# -------------------------------------------------------------- generators ---

def test_make_expected_refuses_dirty_tree(repo, monkeypatch):
    import make_expected
    monkeypatch.setattr(make_expected, "generate", _must_not_run("generate"))
    monkeypatch.setattr(make_expected, "write_baseline_info", _must_not_run("write_baseline_info"))
    _dirty(repo)
    with pytest.raises(SystemExit, match="local changes"):
        make_expected.main(["--tl-root", str(repo), "--allow-non-baseline"])


def test_baseline_info_has_no_absolute_path(tmp_path, monkeypatch):
    import make_expected
    monkeypatch.setattr(make_expected, "CASES_DIR", tmp_path)
    make_expected.write_baseline_info("3419e83abc", False)
    text = (tmp_path / "BASELINE.txt").read_text()
    assert "tree: --tl-root" in text and str(Path.home()) not in text


def test_failing_coverage_writes_no_case(tmp_path, monkeypatch):
    gen = _generator("make_refcalc_cases", monkeypatch, tmp_path)

    def prepare(sysdef, calc):
        work = calc / "work"
        work.mkdir(parents=True)
        (work / "refcalc-FIN").write_text("fin\n")
        (work / "refcalc-PARAM").write_text("### PARAM file for LMAX = 5 ###\nparam\n")
        (work / "refcalc-x.log").write_text("".join(
            f"### STARTING LOG FOR refcalc-part_{e:.2f}eV ###\nLMAX = 5\n"
            for e in (50, 52, 54, 56)))

    monkeypatch.setattr(gen, "prepare_calc_dir", prepare)
    monkeypatch.setattr(gen, "run_viperleed", lambda *a, **k: None)
    monkeypatch.setattr(gen, "set_run", lambda *a: None)
    monkeypatch.setattr(gen, "check_coverage", lambda fin, checks: {c: False for c in checks})
    monkeypatch.setattr(gen, "_viperleed_exe", lambda: "viperleed")
    monkeypatch.setattr(gen.tier2, "viperleed_version", lambda cmd: "0")
    monkeypatch.setattr(gen, "write_case", _must_not_run("write_case"))
    sysdef = types.SimpleNamespace(name="x", must_cover=["nsub_gt_1"],
                                   parameters={"THEO_ENERGIES": "50 56 2"})
    results, names, _ = gen.generate(sysdef, tmp_path / "work", tmp_path)
    assert results == {"nsub_gt_1": False} and not names


@pytest.mark.parametrize("module, args, work", [
    ("make_tier2_expected", ["--calibrate"], "calibrate"),
    ("make_rfactor_cases", ["--verify-only"], "verify"),
])
def test_check_only_modes_refuse_wrong_commit(module, args, work, repo, tmp_path, monkeypatch):
    gen = _generator(module, monkeypatch, tmp_path)
    monkeypatch.setattr(sys.modules["make_refcalc_cases"], "location", lambda key: repo)
    monkeypatch.setattr(gen, work, _must_not_run(work))
    with pytest.raises(GeneratorError, match="not the baseline"):
        gen.main([*args, "--workroot", str(tmp_path)])


def test_hand_cases_match_their_derivation():
    """The committed hand-made cases are what make_hand_cases.py derives
    from the inputs of their source cases."""
    import make_hand_cases
    assert make_hand_cases.main(["--check"]) == 0
