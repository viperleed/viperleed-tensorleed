"""Tier 2: reference calculation + error curves through viperleed.

pytest tests -m tier2 --tl-root <TL tree> [--viperleed-cmd "<cmd>"]

Each system in tests/tier2/systems/ is run with 'viperleed calc -t <tl_root>'
(RUN = 1 5: reference calculation, then delta amplitudes, superpos and
rfactor for the error curves) and compared with expected.toml.
"""
import tomllib
from pathlib import Path

import pytest

from harness import tier2
from harness.cases import load_tolerances

SYSTEMS = tier2.load_systems(Path(__file__).parent / "systems")


@pytest.fixture(scope="session")
def checked_viperleed_cmd(viperleed_cmd):
    """viperleed_cmd, after checking it runs the pinned viperleed version."""
    try:
        tier2.check_viperleed_version(viperleed_cmd,
                                      load_tolerances()["tier2"]["viperleed_version"])
    except tier2.ViperleedFailed as err:
        pytest.fail(str(err), pytrace=False)
    return viperleed_cmd


@pytest.mark.tier2
@pytest.mark.parametrize("system", SYSTEMS, ids=[p.name for p in SYSTEMS])
def test_viperleed_refcalc_and_error_curves(system, tl_root, checked_viperleed_cmd,
                                            tmp_path):
    tol = load_tolerances()["tier2"]
    expected = tier2.from_toml(tomllib.loads((system / "expected.toml").read_text()))
    actual = tier2.run_viperleed_calc(system, tl_root, checked_viperleed_cmd, tmp_path)
    fails = tier2.compare_results(actual, expected, tol["r_atol"], tol["other_atol"])
    assert not fails, "\n".join(fails[:30])
