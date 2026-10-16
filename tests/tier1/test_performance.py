"""Performance A/B tests: candidate (--tl-root) vs baseline (--baseline-tl-root).

Marker 'perf' (opt-in): pytest tests -m perf --tl-root <candidate>
Thresholds and repetitions: [perf] in tests/tolerances.toml.
Short runs are repeated within one sample until the sample takes at least
sample_cpu_s user CPU time. Peak RSS is tested per run. CPU time is tested
per case: the sum of the per-run medians of the case.
Per-run CPU comparisons of these millisecond-scale runs are too noisy (A/A
deviations up to 12 % with 0.3 s samples), and a sum per program dilutes
regressions in one code path (a +30 % slowdown of the composite-layer case
passes it). The benchmark tests CPU per run for large cases.
"""
import pytest

from harness import perf
from harness.cases import load_cases, load_tolerances

from .conftest import case_run_params

PARAMS = case_run_params("refcalc") + case_run_params("rfactor")
REPORT = []  # rows for the terminal summary (see tests/conftest.py)


@pytest.fixture(scope="session")
def thresholds():
    """The [perf] table of the tolerances file."""
    return load_tolerances()["perf"]


@pytest.fixture(scope="session")
def perf_results(tl_root, baseline_tl_root, thresholds, tmp_path_factory):
    """{(case, run): PerfResult} for all runs, computed once per session."""
    workroot = tmp_path_factory.mktemp("perf")
    results = {}
    for kind in ("refcalc", "rfactor"):
        for case in load_cases(kind):
            # the CPU gate sums the runs of a case: each case's samples
            # together take at least case_sample_cpu_s
            sample_s = max(thresholds["sample_cpu_s"],
                           thresholds["case_sample_cpu_s"] / len(case.runs))
            for run in case.runs:
                res = perf.compare_performance(case, run, tl_root, baseline_tl_root,
                                               thresholds["repetitions"], workroot, sample_s)
                results[(case.name, run.name)] = res
                REPORT.append((f"{case.name}/{run.name}", res))
    return results


@pytest.mark.perf
@pytest.mark.parametrize("case,run", PARAMS)
def test_peak_memory(case, run, perf_results, thresholds):
    """The candidate peak RSS of a run is at most the allowed increase above the baseline."""
    res = perf_results[(case.name, run.name)]
    assert perf.rss_ok(res, thresholds), (
        f"peak RSS {res.cand_rss_max} kB vs baseline {res.base_rss_max} kB "
        f"(allowed +{perf.rss_allowed_increase_kb(res.base_rss_max, thresholds):.0f} kB)")


CASES = [pytest.param(c.name, id=c.name)
         for kind in ("refcalc", "rfactor") for c in load_cases(kind)]


@pytest.mark.perf
@pytest.mark.parametrize("case", CASES)
def test_cpu_time_per_case(case, perf_results, thresholds):
    """The summed median user CPU time of a case is within the baseline threshold."""
    rows = [r for (c, _), r in perf_results.items() if c == case]
    assert rows, f"no runs of {case}"
    agg = perf.PerfResult(0, 0, sum(r.cand_cpu_median for r in rows),
                          sum(r.base_cpu_median for r in rows), 0, 0, 0)
    agg.cpu_ratio = agg.cand_cpu_median / agg.base_cpu_median
    REPORT.append((f"SUM {case} ({len(rows)} runs)", agg))
    assert perf.cpu_ok(agg, thresholds), (
        f"{case}: summed median user CPU {agg.cand_cpu_median:.3f} s vs baseline "
        f"{agg.base_cpu_median:.3f} s (ratio {agg.cpu_ratio:.3f})")
