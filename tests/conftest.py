"""pytest options and fixtures shared by all test tiers."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import pathlib
import sys

import pytest

from harness.build import DEBUG, RELEASE
from harness.paths import location, viperleed_command

TESTS_DIR = pathlib.Path(__file__).parent


def pytest_addoption(parser):
    """Add the options selecting the tree under test, the baseline, the build and viperleed."""
    group = parser.getgroup("tensorleed")
    group.addoption("--tl-root", default=str(TESTS_DIR.parent),
                    help="TensErLEED tree under test (default: this repository)")
    group.addoption("--baseline-tl-root", default=None,
                    help="baseline TensErLEED tree for performance A/B tests and the "
                         "build unit tests (default: TL_TEST_BASELINE or "
                         "baseline_tl_root in tests/local.toml)")
    group.addoption("--build", choices=("release", "debug"), default="release",
                    help="build flags for tier-1 correctness tests")
    group.addoption("--viperleed-cmd", default=None,
                    help="command that runs viperleed (tier 2; default: "
                         "TL_TEST_VIPERLEED_CMD, viperleed_cmd in tests/local.toml, "
                         "or viperleed on PATH)")


def pytest_configure(config):
    """Under pytest-xdist, distribute by xdist_group and refuse runs that must be serial.

    Distributing by group keeps the tests that share a session-cached run
    (tier1/conftest.py: run_case) on the same worker; without groups,
    loadgroup distributes like load.
    """
    if getattr(config.option, "dist", "no") == "load":
        config.option.dist = "loadgroup"
    # workers re-parse the command line and only see "--dist loadgroup" if it
    # was given explicitly; the controller passes the switch above on
    if getattr(config, "workerinput", {}).get("tl_loadgroup"):
        config.option.loadgroup = True
    _refuse_serial_only_under_xdist(config)


# Timing tests never run under pytest -n (CPU times need a quiet
# machine; parallel workers compete for cores and memory bandwidth).
SERIAL_ONLY = ("perf", "benchmark")
SERIAL_ONLY_MSG = ("-m {}: performance tests must run serially and alone on the "
                   "machine; run them without -n (pytest-xdist)")


# Tier 2 takes at most TIER2_MAX_WORKERS workers (each viperleed run
# uses N_CORES = 4 of the reference machine's 8 cores).
TIER2_MAX_WORKERS = 2
TIER2_MSG = ("-m tier2 with {} workers: tier 2 takes at most -n {} (each viperleed "
             "run uses 4 cores)")


def _markers_selected(markexpr, markers):
    """Those of markers whose tests the -m expression can select."""
    if not markexpr.strip():
        return list(markers)
    try:
        from _pytest.mark.expression import Expression
        expr = Expression.compile(markexpr)
        return [m for m in markers
                if expr.evaluate(lambda name, /, **kwargs: name == m)]
    except Exception:  # pytest internals changed: the worker check decides
        return []


def _serial_only_selected(markexpr):
    """Markers of SERIAL_ONLY whose tests the -m expression can select."""
    return _markers_selected(markexpr, SERIAL_ONLY)


def _tier2_selected(markexpr):
    """Whether the -m expression can select tier-2 tests."""
    return bool(_markers_selected(markexpr, ("tier2",)))


def _n_workers(tx):
    """Number of workers the pytest-xdist --tx specs start ('3*popen' is 3)."""
    n = 0
    for spec in tx:
        count, star, _ = spec.partition("*")
        n += int(count) if star and count.isdigit() else 1
    return n


def _refuse_serial_only_under_xdist(config):
    """pytest-xdist controller: refuse before any worker starts."""
    if getattr(config.option, "dist", "no") == "no" or hasattr(config, "workerinput"):
        return
    markexpr = config.option.markexpr or ""
    timed = _serial_only_selected(markexpr)
    if timed:
        raise pytest.UsageError(SERIAL_ONLY_MSG.format("/".join(timed)))
    n_workers = _n_workers(getattr(config.option, "tx", None) or [])
    if n_workers > TIER2_MAX_WORKERS and _tier2_selected(markexpr):
        raise pytest.UsageError(TIER2_MSG.format(n_workers, TIER2_MAX_WORKERS))


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    """Fail a timing or tier-2 test that reaches a worker it must not run on.

    Exact backstop for the controller-side check, before any fixture runs a timing.
    """
    if hasattr(item.config, "workerinput"):
        timed = [m for m in SERIAL_ONLY if item.get_closest_marker(m)]
        if timed:
            pytest.fail(SERIAL_ONLY_MSG.format("/".join(timed)), pytrace=False)
        n_workers = item.config.workerinput.get("workercount", 0)
        if n_workers > TIER2_MAX_WORKERS and item.get_closest_marker("tier2"):
            pytest.fail(TIER2_MSG.format(n_workers, TIER2_MAX_WORKERS), pytrace=False)


@pytest.hookimpl(optionalhook=True)
def pytest_configure_node(node):
    """pytest-xdist controller: tell each worker whether loadgroup is on."""
    node.workerinput["tl_loadgroup"] = node.config.getvalue("dist") == "loadgroup"


def pytest_collection_modifyitems(config, items):
    """Mark the harness unit tests as unit and tier1, so every tier-1 run includes them."""
    unit_dir = TESTS_DIR / "unit"
    for item in items:
        if unit_dir in pathlib.Path(item.path).parents:
            item.add_marker(pytest.mark.unit)
            item.add_marker(pytest.mark.tier1)


@pytest.fixture(scope="session")
def tl_root(request):
    """The TensErLEED tree under test (--tl-root)."""
    return pathlib.Path(request.config.getoption("--tl-root")).resolve()


@pytest.fixture(scope="session")
def baseline_tl_root(request):
    """The baseline TensErLEED tree; raise if it has no TensErLEED/ directory."""
    root = request.config.getoption("--baseline-tl-root") or location("baseline_tl_root")
    root = pathlib.Path(root).resolve()
    if not (root / "TensErLEED").is_dir():
        raise FileNotFoundError(f"baseline TL tree {root} has no TensErLEED/")
    return root


@pytest.fixture(scope="session")
def build_flags(request):
    """The BuildFlags of the --build choice for tier-1 correctness tests."""
    return {"release": RELEASE, "debug": DEBUG}[request.config.getoption("--build")]


@pytest.fixture(scope="session")
def viperleed_cmd(request):
    """The command (list of words) that runs viperleed for tier 2."""
    return viperleed_command(request.config.getoption("--viperleed-cmd"))


def pytest_terminal_summary(terminalreporter):
    """Print the performance table if performance tests ran."""
    mod = sys.modules.get("tier1.test_performance") or sys.modules.get("test_performance")
    rows = getattr(mod, "REPORT", None)
    if not rows:
        return
    tr = terminalreporter
    tr.section("performance (candidate vs baseline)")
    tr.write_line(f"{'run':36s} {'inner':>5s} {'cpu cand':>9s} {'cpu base':>9s} {'ratio':>6s} "
                  f"{'rss cand':>9s} {'rss base':>9s} {'ratio':>6s}")
    for name, r in rows:
        tr.write_line(f"{name:36s} {r.inner:5d} {r.cand_cpu_median:9.3f} {r.base_cpu_median:9.3f} "
                      f"{r.cpu_ratio:6.3f} {r.cand_rss_max:9d} {r.base_rss_max:9d} "
                      f"{r.rss_ratio:6.3f}")
