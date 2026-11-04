"""Debug build (bounds checks, FP traps): programs must terminate cleanly.

Known baseline failures are listed in tests/known_debug_failures.toml: those
runs must fail with exactly the known error, and must run clean with all
checks except bounds. Runs with debug = false in their case.toml are not run
(undefined behaviour shared with the legacy code; README, 'Rules'). Outputs are compared with relaxed tolerances
(rtol >= [debug].min_rtol in tests/tolerances.toml) as a sanity
check only. Before the first debug run, a small Fortran probe built with the
same flags must stop at an out-of-bounds index and at each trapped
floating-point exception (invalid, zero, overflow).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import dataclasses
import signal
import tomllib

import pytest

from harness.build import DEBUG, DEBUG_NO_BOUNDS
from harness.cases import (TESTS_DIR, compare_outputs, execute_run, expected_outputs,
                           format_failures, iter_runs, load_tolerances)
from harness.runner import ProgramFailed, trap_probe

KNOWN = {(f["case"], f["run"]): f for f in tomllib.loads(
    (TESTS_DIR / "known_debug_failures.toml").read_text()).get("failure", [])}


# fault of the trap probe: (exit code, text in stderr) when the debug build stops it
_TRAPS = {"bounds": (2, "Fortran runtime error: Index '5' of dimension 1 of array 'a' "
                        "above upper bound of 4"),
          "invalid": (-signal.SIGFPE, "Program received signal SIGFPE"),
          "zero": (-signal.SIGFPE, "Program received signal SIGFPE"),
          "overflow": (-signal.SIGFPE, "Program received signal SIGFPE")}


def _untrapped(workdir, fault):
    """None if the probe built with DEBUG stops at fault as expected, else what it did."""
    res = trap_probe(workdir, DEBUG, fault)
    code, text = _TRAPS[fault]
    if res.returncode == code and text in res.stderr:
        return None
    return (f"{fault}: exit code {res.returncode}, stdout {res.stdout.strip()!r}, "
            f"stderr tail {res.stderr[-300:]!r}")


@pytest.fixture(scope="module")
def debug_traps(tmp_path_factory):
    """Check that the DEBUG flags stop the trap probe at every fault, before any debug run."""
    workdir = tmp_path_factory.mktemp("trap_probe")
    failures = [f for f in (_untrapped(workdir, fault) for fault in _TRAPS) if f]
    assert not failures, "debug flags do not trap:\n" + "\n".join(failures)


def _params():
    for case, run in iter_runs():
        if run.debug:
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
@pytest.mark.parametrize("fault", list(_TRAPS))
def test_debug_flags_stop_trap_probe(fault, tmp_path):
    """The trap probe built with the DEBUG flags stops at the fault (bounds, invalid, zero, overflow).

    Failing means DEBUG lost a bounds check or a floating-point trap, and the
    debug runs below (which error) would pass without that check."""
    assert _untrapped(tmp_path, fault) is None


@pytest.mark.debug
@pytest.mark.parametrize("case,run", list(_params()))
def test_debug_build_runs_clean(case, run, tl_root, debug_traps, tmp_path):
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
