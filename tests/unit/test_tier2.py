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
