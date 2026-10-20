"""Tests of the viperleed result parsing, comparison and version pin of harness.tier2."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import re
import tomllib
import zipfile

import pytest

from harness import tier2
from harness.cases import CaseError, load_tolerances
from unit.gitrepo import git

ROUT = """tl-test
AVERAGE                 0  0.0001    3.00   45.00  216.05    479.50  0.2054  <---

BEST GEOMETRY ENCOUNTERED, NO.    1, CHARACTERIZED BY VALUE  0.0001
WITH INNER POTENTIAL SHIFT    3.00 EV, AVERAGE R-FACTOR =  0.2054
"""
CSV = """ Displacement (z) [Å],      R
              -0.0200, 0.1706
               0.0000, 0.2054
               0.0200, 0.2760
"""


def _workdir(tmp_path, csv=CSV, rout=ROUT):
    out = tmp_path / "OUT"
    out.mkdir(parents=True)
    (out / "R_OUT_refcalc_R=0.2054").write_text(rout)
    with zipfile.ZipFile(out / "Errors.zip", "w") as z:
        z.writestr("Errors_geo_atoms#1_z.csv", csv)
    return tmp_path


def test_parse_results(tmp_path):
    """R-factor, inner-potential shift and error curves are read from a work directory."""
    res = tier2.parse_results(_workdir(tmp_path))
    assert res["refcalc"] == {"r": 0.2054, "v0_shift": 3.0}
    curve = res["errors"]["Errors_geo_atoms#1_z.csv"]
    assert curve == {"Displacement (z) [Å]": [-0.02, 0.0, 0.02], "R": [0.1706, 0.2054, 0.276]}


def test_toml_roundtrip_and_compare(tmp_path):
    """Results survive a TOML round trip and compare equal to themselves."""
    res = tier2.parse_results(_workdir(tmp_path))
    back = tier2.from_toml(tomllib.loads(tier2.to_toml(res, "comment")))
    assert back == res
    assert tier2.compare_results(res, back, 1.5e-4, 1e-6) == []


def test_compare_detects_changed_r_and_missing_curve(tmp_path):
    """A changed R value and a missing error curve are both reported."""
    exp = tier2.parse_results(_workdir(tmp_path / "e"))
    act = tier2.parse_results(_workdir(tmp_path / "a", csv=CSV.replace("0.2760", "0.2763")))
    fails = tier2.compare_results(act, exp, 1.5e-4, 1e-6)
    assert len(fails) == 1 and "R[2]" in fails[0]
    act["errors"] = {}
    assert any("error curves differ" in f for f in tier2.compare_results(act, exp, 1.5e-4, 1e-6))


def test_missing_outputs_fail(tmp_path):
    """A work directory without viperleed outputs raises ViperleedFailed."""
    with pytest.raises(tier2.ViperleedFailed):
        tier2.parse_results(tmp_path)


def test_set_parameters_replaces_and_appends(tmp_path):
    """set_parameters replaces existing parameters in place and appends new ones."""
    p = tmp_path / "PARAMETERS"
    p.write_text("RUN = 1 5\nN_CORES = 2\n")
    tier2.set_parameters(p, {"N_CORES": "4", "FORTRAN_COMP": "'gfortran -O0'"})
    assert p.read_text() == "RUN = 1 5\nN_CORES = 4\nFORTRAN_COMP = 'gfortran -O0'\n"


@pytest.mark.parametrize("side", ["actual", "expected"])
def test_compare_fails_on_nan(tmp_path, side):
    """A NaN in either the actual or the expected curve fails the comparison."""
    nan_csv = CSV.replace("0.2760", "nan")
    exp = tier2.parse_results(_workdir(tmp_path / "e", csv=nan_csv if side == "expected" else CSV))
    act = tier2.parse_results(_workdir(tmp_path / "a", csv=nan_csv if side == "actual" else CSV))
    fails = tier2.compare_results(act, exp, 1.5e-4, 1e-6)
    assert len(fails) == 1 and "R[2]" in fails[0]


@pytest.mark.parametrize("row", ["0.0200", "0.0200, 0.2760, 0.1"], ids=["short", "long"])
def test_row_length_mismatch_fails(row):
    """A CSV row with another number of columns than the header raises."""
    with pytest.raises(tier2.ViperleedFailed, match="header has 2"):
        tier2.parse_error_csv(CSV.replace("0.0200, 0.2760", row))


BANNER = "ViPErLEED (Vienna Package for Erlangen LEED) v{}"


@pytest.mark.parametrize("suffix", ["0", ".dev0", "+g1"])
def test_viperleed_version_mismatch_stops(mocker, suffix):
    """A viperleed other than the pinned version, even one that starts with it,
    stops with a message naming both versions."""
    pinned = load_tolerances()["tier2"]["viperleed_version"]
    for installed in ("0.13.0", pinned + suffix):
        version = mocker.patch.object(tier2, "viperleed_version",
                                      return_value=BANNER.format(installed))
        with pytest.raises(tier2.ViperleedFailed,
                           match=rf"{re.escape(pinned)}.*runs {re.escape(installed)}$"):
            tier2.check_viperleed_version(["viperleed"], pinned)
        version.assert_called_once_with(["viperleed"])


def test_viperleed_version_match_passes(mocker):
    """The pinned viperleed version passes."""
    pinned = load_tolerances()["tier2"]["viperleed_version"]
    version = mocker.patch.object(tier2, "viperleed_version", return_value=BANNER.format(pinned))
    tier2.check_viperleed_version(["viperleed"], pinned)
    version.assert_called_once_with(["viperleed"])


def test_system_without_expected_toml_raises(tmp_path):
    """A system directory without expected.toml raises instead of dropping out of the tests."""
    for name in ("a", "b"):
        (tmp_path / name).mkdir()
    (tmp_path / "a" / "expected.toml").write_text("")
    with pytest.raises(CaseError, match=r"without expected.toml.*\['b'\]"):
        tier2.load_systems(tmp_path)


def test_empty_systems_dir_raises(tmp_path):
    """A directory without systems raises, ignoring non-system entries such as __pycache__."""
    (tmp_path / "__pycache__").mkdir()
    with pytest.raises(CaseError, match="no systems"):
        tier2.load_systems(tmp_path)


def test_load_systems_sorted(tmp_path):
    """Systems are returned sorted by name."""
    for name in ("b", "a"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "expected.toml").write_text("")
    assert [p.name for p in tier2.load_systems(tmp_path)] == ["a", "b"]


def _tl_tree(path, tag=None, modified=False):
    path.mkdir(parents=True)
    (path / "f90").mkdir()
    (path / "f90" / "a.f90").write_text("original\n")
    git(path, "init", "-q")
    git(path, "add", "f90")
    git(path, "commit", "-q", "-m", "c")
    if tag:
        git(path, "tag", tag)
    if modified:
        (path / "f90" / "a.f90").write_text("changed\n")
    return path


@pytest.mark.parametrize("kwargs, released", [
    ({"tag": "v2.1.0"}, True),
    ({"tag": "v2.1.0", "modified": True}, False),
    ({"tag": "v2.1.0-rc1"}, False),
    ({"tag": "nightly"}, False),
    ({}, False),
], ids=["released", "modified", "prerelease-tag", "other-tag", "untagged"])
def test_is_released_tl_tree(tmp_path, kwargs, released):
    """Only an unmodified checkout at a vX.Y.Z tag counts as a released TL tree."""
    assert tier2.is_released_tl_tree(_tl_tree(tmp_path / "tl", **kwargs)) is released


def test_tree_without_git_is_not_released(tmp_path, monkeypatch):
    """A directory that is not in any git checkout is not a released tree."""
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    (tmp_path / "tl").mkdir()
    assert tier2.is_released_tl_tree(tmp_path / "tl") is False


def test_folder_inside_released_tree_is_not_released(tmp_path):
    """A folder without git of its own inside a released checkout is not released."""
    tl = _tl_tree(tmp_path / "tl", tag="v2.1.0")
    (tl / "nested").mkdir()
    assert tier2.is_released_tl_tree(tl / "nested") is False


def test_subdirectory_of_released_tree_is_not_released(tmp_path):
    """A tracked subdirectory of a released checkout is not a released tree."""
    tl = _tl_tree(tmp_path / "tl", tag="v2.1.0")
    assert tier2.is_released_tl_tree(tl / "f90") is False


def test_worktree_at_release_tag_is_released(tmp_path):
    """A git worktree (with a .git file) at a release tag is a released tree."""
    tl = _tl_tree(tmp_path / "tl", tag="v2.1.0")
    git(tl, "worktree", "add", "-q", "--detach", str(tmp_path / "wt"), "v2.1.0")
    assert tier2.is_released_tl_tree(tmp_path / "wt") is True


def _point_git_env_at(monkeypatch, repo):
    monkeypatch.setenv("GIT_DIR", str(repo / ".git"))
    monkeypatch.setenv("GIT_WORK_TREE", str(repo))


def test_git_env_of_released_tree_does_not_release_other_folder(tmp_path, monkeypatch):
    """GIT_DIR/GIT_WORK_TREE naming a released tree do not make another folder released."""
    _point_git_env_at(monkeypatch, _tl_tree(tmp_path / "tl", tag="v2.1.0"))
    (tmp_path / "plain").mkdir()
    assert tier2.is_released_tl_tree(tmp_path / "plain") is False


def test_git_env_of_other_repo_does_not_hide_release(tmp_path, monkeypatch):
    """GIT_DIR/GIT_WORK_TREE naming an untagged repo do not hide a released tree."""
    tl = _tl_tree(tmp_path / "tl", tag="v2.1.0")
    _point_git_env_at(monkeypatch, _tl_tree(tmp_path / "other"))
    assert tier2.is_released_tl_tree(tl) is True


def _system(tmp_path):
    system = tmp_path / "system"
    system.mkdir()
    (system / "PARAMETERS").write_text("RUN = 1 5\nTL_IGNORE_CHECKSUM = True\n")
    return system


@pytest.mark.parametrize("tag, expected", [("v2.1.0", "False"), (None, "True")],
                         ids=["released", "untagged"])
def test_checksum_check_follows_release_status(mocker, tmp_path, tag, expected):
    """viperleed's checksum check stays on for a released tl_root and is off otherwise.

    The frozen PARAMETERS of the system say TL_IGNORE_CHECKSUM = True either way."""
    tl = _tl_tree(tmp_path / "tl", tag=tag)
    mocker.patch.object(tier2, "parse_results")
    tier2.run_viperleed_calc(_system(tmp_path), tl, ["true"], tmp_path / "work")
    text = (tmp_path / "work" / "PARAMETERS").read_text()
    assert f"TL_IGNORE_CHECKSUM = {expected}\n" in text
    assert text.count("TL_IGNORE_CHECKSUM") == 1


def test_caller_cannot_set_checksum_parameter(tmp_path):
    """A TL_IGNORE_CHECKSUM in the caller's parameters raises ValueError."""
    with pytest.raises(ValueError, match="TL_IGNORE_CHECKSUM"):
        tier2.run_viperleed_calc(_system(tmp_path), tmp_path, ["true"], tmp_path / "work",
                                 parameters={"TL_IGNORE_CHECKSUM": "True"})
