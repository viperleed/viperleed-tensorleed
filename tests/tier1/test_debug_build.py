"""Debug build (bounds checks, FP traps): programs must terminate cleanly.

Known baseline failures are listed in tests/known_debug_failures.toml: those
runs must fail with exactly the known error, and must run clean with all
checks except bounds. Runs with debug = false in their case.toml are not run
(undefined behaviour shared with the legacy code; README, 'Rules'). Outputs are compared with relaxed tolerances
(rtol >= [debug].min_rtol in tests/tolerances.toml) as a sanity
check only. Before the first debug run, a small Fortran probe built with the
same flags must stop at an out-of-bounds index and at each trapped
floating-point exception (invalid, zero, overflow), and in Makefile mode
every object of the tree's Makefile build of ref-calc and rfactor must record
the debug flags, without a later override, as its compile options.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import dataclasses
import re
import shutil
import signal
import subprocess
import tomllib

import pytest

from harness.build import DEBUG, DEBUG_NO_BOUNDS, uses_makefile
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


# option words that decide the checks and optimisation of an object
_CHECK_OPTION = re.compile(r"-(fcheck=|ffpe-trap=|fno-|O|g)")


def _checked_options(words):
    return sorted(w for w in words if _CHECK_OPTION.match(w))


def _producer(obj):
    """The compiler and its effective options as recorded in obj (DW_AT_producer), or None."""
    info = subprocess.run(["readelf", "--debug-dump=info", str(obj)],
                          capture_output=True, text=True).stdout
    match = re.search(r"DW_AT_producer\s*:.*?:\s*(GNU Fortran.*)", info)
    return match.group(1) if match else None


def _makefile_without_debug_flags(tl_root, workdir):
    """What is wrong with the objects of a Makefile build of ref-calc and rfactor.

    Both programs are built with DEBUG's flags; every object must record in
    its DW_AT_producer exactly DEBUG's checking, optimisation and debug-info
    options, so a recipe that drops or overrides them (later -fcheck=no-...,
    -fno-..., -ffpe-trap=, another -O) is reported. Also reports a program
    that builds no object.
    """
    param = next(run.param for case, run in iter_runs() if case.program == "refcalc")
    tree = workdir / "tree" / "TensErLEED"
    shutil.copytree(tl_root / "TensErLEED", tree, symlinks=True)
    want = _checked_options(DEBUG.compile.split())
    bad = []
    for program in ("refcalc", "rfactor"):
        builddir = workdir / program
        proc = subprocess.run(
            ["make", "-B", "-C", str(tree), program, f"BUILDDIR={builddir}",
             "FC=gfortran", f"FFLAGS={DEBUG.compile}", f"LDLIBS={DEBUG.link}",
             f"PARAM={param}"],
            capture_output=True, text=True)
        if proc.returncode != 0:
            bad.append(f"{program}: make failed:\n{proc.stderr[-2000:]}")
            continue
        objects = sorted(builddir.rglob("*.o"))
        if not objects:
            bad.append(f"{program}: make builds no object")
        for obj in objects:
            producer = _producer(obj)
            if producer is None:
                bad.append(f"{program}: {obj.name} has no DW_AT_producer (compiled without -g?)")
            elif _checked_options(producer.split()) != want:
                bad.append(f"{program}: {obj.name} was compiled with: {producer}")
    return bad


@pytest.fixture(scope="module")
def makefile_debug_problems(tl_root, tmp_path_factory):
    """The problems of the Makefile build with DEBUG's flags (built once)."""
    return _makefile_without_debug_flags(tl_root, tmp_path_factory.mktemp("makefile_debug"))


@pytest.fixture(scope="module")
def debug_traps(tl_root, tmp_path_factory, makefile_debug_problems):
    """Check, before any debug run, that the DEBUG flags stop the trap probe at
    every fault and that the Makefile build applies them."""
    workdir = tmp_path_factory.mktemp("trap_probe")
    failures = [f for f in (_untrapped(workdir, fault) for fault in _TRAPS) if f]
    assert not failures, "debug flags do not trap:\n" + "\n".join(failures)
    if uses_makefile(tl_root):
        assert not makefile_debug_problems, ("the Makefile build drops the debug flags:\n"
                                             + "\n".join(makefile_debug_problems))


@pytest.mark.tier1
@pytest.mark.debug
def test_makefile_applies_debug_flags(tl_root, makefile_debug_problems):
    """Every object of the Makefile build of ref-calc and rfactor was compiled
    with the DEBUG flags and no later override (read from its DW_AT_producer).

    Failing means a recipe drops or overrides $(FFLAGS), and the debug builds
    of that program would run without bounds checks or FP traps."""
    if not uses_makefile(tl_root):
        pytest.skip("legacy tree: the harness compiles with the flags itself")
    assert not makefile_debug_problems, "\n".join(makefile_debug_problems)


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


@pytest.mark.tier1
@pytest.mark.debug
@pytest.mark.parametrize("fault", list(_TRAPS))
def test_debug_flags_stop_trap_probe(fault, tmp_path):
    """The trap probe built with the DEBUG flags stops at the fault (bounds, invalid, zero, overflow).

    Failing means DEBUG lost a bounds check or a floating-point trap, and the
    debug runs below, and the tier-1 tests built with DEBUG, would pass
    without that check."""
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
