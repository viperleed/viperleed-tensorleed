"""Large-system benchmark (opt-in): pytest tests -m benchmark --tl-root <candidate>

Cases: tests/cases/benchmark/ (inputs only). For each run, the candidate
(--tl-root) and the baseline (--baseline-tl-root) are compared
- A/B in user CPU time and peak RSS (perf thresholds, [perf].benchmark_repetitions
  interleaved repetitions each), and
- in their outputs (candidate vs a fresh baseline run, tier-1 tolerances).
"""
import pytest

from harness import perf
from harness.build import RELEASE
from harness.cases import (CASES_DIR, compare_outputs, execute_run, file_set_failure,
                           format_failures, load_cases_from, load_tolerances)

PARAMS = [pytest.param(c, r, id=f"{c.name}/{r.name}")
          for c in load_cases_from(CASES_DIR, "benchmark") for r in c.runs]


@pytest.mark.benchmark
@pytest.mark.parametrize("case,run", PARAMS)
def test_benchmark_outputs_match_fresh_baseline(case, run, tl_root, baseline_tl_root, tmp_path):
    """Candidate outputs of a benchmark run match those of a fresh baseline run.

    Failing means numerical results or the set of output files changed on a large system."""
    base, cand = tmp_path / "base", tmp_path / "cand"
    execute_run(case, run, baseline_tl_root, RELEASE, base)
    execute_run(case, run, tl_root, RELEASE, cand)
    files = sorted({p.name for pat in run.outputs for p in base.glob(pat) if p.is_file()})
    assert files, "baseline produced no outputs"
    failures = format_failures(compare_outputs(cand, [base / name for name in files]))
    problem = file_set_failure(cand, run.outputs, files)
    if problem:
        failures.append(f"output files differ: {problem}")
    assert not failures, "\n".join(failures)


@pytest.mark.benchmark
@pytest.mark.parametrize("case,run", PARAMS)
def test_benchmark_performance(case, run, tl_root, baseline_tl_root, tmp_path):
    """Candidate CPU time and peak RSS of a benchmark run are within the baseline thresholds."""
    th = load_tolerances()["perf"]
    res = perf.compare_performance(case, run, tl_root, baseline_tl_root,
                                   th["benchmark_repetitions"], tmp_path)
    print(f"\n{case.name}/{run.name}: CPU {res.cand_cpu_median:.2f} s vs {res.base_cpu_median:.2f} s "
          f"(ratio {res.cpu_ratio:.3f}); RSS {res.cand_rss_max} kB vs {res.base_rss_max} kB")
    assert perf.cpu_ok(res, th), f"CPU ratio {res.cpu_ratio:.3f}"
    assert perf.rss_ok(res, th), f"RSS {res.cand_rss_max} vs {res.base_rss_max} kB"


@pytest.mark.tier1
@pytest.mark.parametrize("reps", [2, 4])
def test_benchmark_repetitions_from_config(reps, mocker):
    """compare_performance is called with the benchmark_repetitions of the configuration."""
    compare = mocker.patch.object(perf, "compare_performance",
                                  return_value=perf.PerfResult(1, 1, 1, 1, 1, 1))
    th = {**load_tolerances()["perf"], "benchmark_repetitions": reps}
    mocker.patch.dict(globals(), {"load_tolerances": lambda: {"perf": th}})
    case, run = PARAMS[0].values
    test_benchmark_performance(case, run, "cand", "base", "w")
    compare.assert_called_once_with(case, run, "cand", "base", reps, "w")
