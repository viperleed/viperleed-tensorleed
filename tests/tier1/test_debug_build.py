"""Debug build (bounds checks, FP traps): programs must terminate cleanly.

Known baseline failures are listed in tests/known_debug_failures.toml: those
runs must fail with exactly the known error, and must run clean with all
checks except bounds. Outputs are compared with relaxed tolerances
(rtol >= [debug].min_rtol in tests/tolerances.toml) as a sanity
check only.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import dataclasses
import tomllib

import pytest

from harness.build import DEBUG, DEBUG_NO_BOUNDS
from harness.cases import (TESTS_DIR, compare_outputs, execute_run, expected_outputs,
                           format_failures, iter_runs, load_tolerances)
from harness.runner import ProgramFailed

KNOWN = {(f["case"], f["run"]): f for f in tomllib.loads(
    (TESTS_DIR / "known_debug_failures.toml").read_text())["failure"]}


def _params():
    for case, run in iter_runs():
        yield pytest.param(case, run, id=f"{case.name}/{run.name}")


@pytest.mark.tier1
def test_known_failures_exist_as_runs():
    """Every entry of the known debug failures names an existing case and run."""
    runs = {(c.name, r.name) for c, r in iter_runs()}
    assert set(KNOWN) <= runs


def _compare_relaxed(case, run, workdir):
    min_rtol = load_tolerances()["debug"]["min_rtol"]
    return format_failures(compare_outputs(
        workdir, expected_outputs(case, run),
        adjust=lambda tol: dataclasses.replace(tol, rtol=max(tol.rtol, min_rtol))))


@pytest.mark.tier1
@pytest.mark.parametrize("min_rtol,ok", [(1e-3, True), (1e-4, False)])
def test_relaxed_rtol_from_config(min_rtol, ok, mocker, tmp_path):
    """The debug comparison uses the configured minimum rtol (a 5e-4 deviation passes at 1e-3)."""
    # ROUT: rtol 0, atol 1.5e-4; a relative deviation of 5e-4 passes only with min_rtol 1e-3
    (tmp_path / "exp").mkdir()
    (tmp_path / "exp" / "ROUT").write_text("R = 1.0000\n")
    (tmp_path / "ROUT").write_text("R = 1.0005\n")
    mocker.patch.dict(globals(), {
        "expected_outputs": lambda c, r: [tmp_path / "exp" / "ROUT"],
        "load_tolerances": lambda: {"debug": {"min_rtol": min_rtol}}})
    assert (not _compare_relaxed(None, None, tmp_path)) == ok


@pytest.mark.debug
@pytest.mark.parametrize("case,run", list(_params()))
def test_debug_build_runs_clean(case, run, tl_root, tmp_path):
    """The debug build (bounds checks, FP traps) terminates cleanly and gives sane outputs.

    A run with a known failure must fail with exactly the recorded error and then run clean
    without bounds checks. Failing means a new runtime error, or a known one that vanished
    or changed (update the known failures only for a deliberate fix)."""
    known = KNOWN.get((case.name, run.name))
    if known:
        with pytest.raises(ProgramFailed, match=known["match"]):
            execute_run(case, run, tl_root, DEBUG, tmp_path / "debug")
        workdir = tmp_path / "nobounds"
        execute_run(case, run, tl_root, DEBUG_NO_BOUNDS, workdir)
    else:
        workdir = tmp_path
        execute_run(case, run, tl_root, DEBUG, workdir)
    failures = _compare_relaxed(case, run, workdir)
    assert not failures, "\n".join(failures)
