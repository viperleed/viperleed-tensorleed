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

from generate import baseline_tree, make_expected
from generate.baseline_tree import GeneratorError, check_baseline_tree


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
def at_baseline(repo, mocker):
    """Return the repo with BASELINE_COMMIT patched to its HEAD."""
    mocker.patch.object(baseline_tree, "BASELINE_COMMIT", _git(repo, "rev-parse", "HEAD"))
    return repo


def _dirty(tree):
    (tree / "TensErLEED" / "a.f").write_text("      STOP\n      END\n")


def _generator(name, monkeypatch, mocker, tmp_path):
    mocker.patch.dict(sys.modules, {"fortranformat": types.ModuleType("fortranformat")})
    monkeypatch.setenv("TL_TEST_VIPERLEED_EXAMPLES", str(tmp_path))
    monkeypatch.setenv("TL_TEST_BASELINE", str(tmp_path))
    return importlib.import_module(f"generate.{name}")


# ------------------------------------------------------------ the checker ---

def test_wrong_commit_refused(repo):
    """A tree not at the baseline commit is refused."""
    with pytest.raises(GeneratorError, match="not the baseline"):
        check_baseline_tree(repo)


def test_dirty_tree_refused(at_baseline):
    """A modified tracked file is refused."""
    _dirty(at_baseline)
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(at_baseline)


def test_dirty_tree_refused_with_other_commit_allowed(repo):
    """Local changes are refused even when another commit is allowed."""
    _dirty(repo)
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(repo, allow_other_commit=True)


def test_untracked_file_refused(at_baseline):
    """An untracked file under TensErLEED/ is refused."""
    (at_baseline / "TensErLEED" / "b.f").write_text("      END\n")
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(at_baseline)


def test_ignored_file_refused(at_baseline):
    """A git-ignored file under TensErLEED/ is refused."""
    (at_baseline / ".git" / "info" / "exclude").write_text("*.f90\n")
    (at_baseline / "TensErLEED" / "b.f90").write_text("end\n")
    with pytest.raises(GeneratorError, match="local changes"):
        check_baseline_tree(at_baseline)


def test_ignored_subdirectory_of_baseline_refused(at_baseline):
    """An ignored copy of TensErLEED/ inside a checkout at the baseline
    must not pass as that checkout.
    """
    sub = at_baseline / "foo"
    (sub / "TensErLEED").mkdir(parents=True)
    (sub / ".gitignore").write_text("*\n")
    (sub / "TensErLEED" / "a.f").write_text("      STOP\n      END\n")
    with pytest.raises(GeneratorError, match="not the top level"):
        check_baseline_tree(sub)


def test_clean_baseline_accepted(at_baseline):
    """A clean tree at the baseline commit is accepted and its commit returned."""
    assert check_baseline_tree(at_baseline) == _git(at_baseline, "rev-parse", "HEAD")


def test_expected_outputs_come_from_the_pinned_baseline():
    """cases/BASELINE.txt names the commit the generators are pinned to."""
    info = (Path(__file__).parents[1] / "cases" / "BASELINE.txt").read_text()
    (commit,) = [l.split(":", 1)[1].strip() for l in info.splitlines()
                 if l.startswith("commit:")]
    assert commit.startswith(baseline_tree.BASELINE_COMMIT), (commit, baseline_tree.BASELINE_COMMIT)


# -------------------------------------------------------------- generators ---

def test_make_expected_refuses_dirty_tree(repo, mocker):
    """make_expected stops on a dirty tree before generating or writing anything."""
    generate = mocker.patch.object(make_expected, "generate")
    write_info = mocker.patch.object(make_expected, "write_baseline_info")
    _dirty(repo)
    with pytest.raises(SystemExit, match="local changes"):
        make_expected.main(["--tl-root", str(repo), "--allow-non-baseline"])
    generate.assert_not_called()
    write_info.assert_not_called()


def test_baseline_info_has_no_absolute_path(tmp_path, mocker):
    """BASELINE.txt records no absolute path of the machine."""
    mocker.patch.object(make_expected, "CASES_DIR", tmp_path)
    make_expected.write_baseline_info("3419e83abc", False)
    text = (tmp_path / "BASELINE.txt").read_text()
    assert "tree: --tl-root" in text and str(Path.home()) not in text


def test_failing_coverage_writes_no_case(tmp_path, monkeypatch, mocker):
    """A system that fails its coverage checks produces no case."""
    gen = _generator("make_refcalc_cases", monkeypatch, mocker, tmp_path)

    def prepare(sysdef, calc):
        """Write the files a viperleed run would leave in the calculation directory."""
        work = calc / "work"
        work.mkdir(parents=True)
        (work / "refcalc-FIN").write_text("fin\n")
        (work / "refcalc-PARAM").write_text("### PARAM file for LMAX = 5 ###\nparam\n")
        (work / "refcalc-x.log").write_text("".join(
            f"### STARTING LOG FOR refcalc-part_{e:.2f}eV ###\nLMAX = 5\n"
            for e in (50, 52, 54, 56)))

    mocker.patch.object(gen, "prepare_calc_dir", side_effect=prepare)
    mocker.patch.object(gen, "run_viperleed")
    mocker.patch.object(gen, "set_run")
    mocker.patch.object(gen, "check_coverage",
                        side_effect=lambda fin, checks: {c: False for c in checks})
    mocker.patch.object(gen, "_viperleed_exe", return_value="viperleed")
    mocker.patch.object(gen.tier2, "viperleed_version", return_value="0")
    write_case = mocker.patch.object(gen, "write_case")
    sysdef = types.SimpleNamespace(name="x", must_cover=["nsub_gt_1"],
                                   parameters={"THEO_ENERGIES": "50 56 2"})
    results, names, _ = gen.generate(sysdef, tmp_path / "work", tmp_path)
    assert results == {"nsub_gt_1": False} and not names
    write_case.assert_not_called()


@pytest.mark.parametrize("module, args, work", [
    ("make_tier2_expected", ["--calibrate"], "calibrate"),
    ("make_rfactor_cases", ["--verify-only"], "verify"),
])
def test_check_only_modes_refuse_wrong_commit(module, args, work, repo, tmp_path,
                                              monkeypatch, mocker):
    """The check-only generator modes refuse a tree at the wrong commit before any work."""
    gen = _generator(module, monkeypatch, mocker, tmp_path)
    mocker.patch.object(sys.modules["generate.make_refcalc_cases"], "location",
                        return_value=repo)
    action = mocker.patch.object(gen, work)
    with pytest.raises(GeneratorError, match="not the baseline"):
        gen.main([*args, "--workroot", str(tmp_path)])
    action.assert_not_called()
