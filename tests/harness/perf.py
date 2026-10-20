"""A/B performance comparison (user CPU time, peak RSS) of two TL trees.

Thresholds come from [perf] in tests/tolerances.toml:
CPU: candidate median user CPU <= (1 + cpu_rel_threshold) * baseline median;
RSS: candidate max RSS - baseline max RSS <= max(rss_rel_threshold * baseline,
rss_abs_floor_kb).
Sample sizing: sizing_runs, min_sizing_cpu_s, max_inner (same table).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import math
import shutil
import statistics
from dataclasses import dataclass
from pathlib import Path

from . import cases
from .build import RELEASE


@dataclass
class PerfResult:
    """Candidate-to-baseline CPU and RSS ratios with the underlying medians and maxima."""
    cpu_ratio: float
    rss_ratio: float
    cand_cpu_median: float
    base_cpu_median: float
    cand_rss_max: int
    base_rss_max: int
    inner: int = 1  # executions per sample


def compare_performance(case, run, candidate_root, baseline_root, reps, workroot,
                        sample_cpu_s=0.0):
    """Run baseline and candidate interleaved (B, C, B, C, ...) reps times each.

    One sample = `inner` executions in fresh run directories under workroot,
    with `inner` chosen from the median of sizing_runs baseline executions so that
    a sample takes at least sample_cpu_s user CPU time (short runs are otherwise
    dominated by timer resolution and start-up noise; [perf] in tolerances.toml).
    Returns the medians of the per-sample user CPU times and the maxima of peak
    RSS.
    """
    workroot = Path(workroot)
    exes = {"base": cases.build_program(baseline_root, case.program, run.param, RELEASE),
            "cand": cases.build_program(candidate_root, case.program, run.param, RELEASE)}
    counter = iter(range(10**9))

    def once(label):
        """Execute the program of label ("base" or "cand") once in a fresh run directory."""
        workdir = workroot / f"{case.name}-{run.name}-{label}-{next(counter)}"
        shutil.rmtree(workdir, ignore_errors=True)
        res = cases.run_built(case, run, exes[label], workdir)
        shutil.rmtree(workdir, ignore_errors=True)
        return res

    cfg = cases.load_tolerances()["perf"]
    first = [once("base") for _ in range(cfg["sizing_runs"])]
    inner = 1
    if sample_cpu_s > 0:
        typical = max(statistics.median(r.user_cpu_s for r in first), cfg["min_sizing_cpu_s"])
        inner = min(cfg["max_inner"], max(1, math.ceil(sample_cpu_s / typical)))
    cpu = {"base": [], "cand": []}
    rss = {"base": [r.max_rss_kb for r in first], "cand": []}
    for _ in range(reps):
        for label in ("base", "cand"):
            total = 0.0
            for _ in range(inner):
                res = once(label)
                total += res.user_cpu_s
                rss[label].append(res.max_rss_kb)
            cpu[label].append(total)
    base_cpu, cand_cpu = statistics.median(cpu["base"]), statistics.median(cpu["cand"])
    base_rss, cand_rss = max(rss["base"]), max(rss["cand"])
    return PerfResult(cand_cpu / base_cpu if base_cpu > 0 else float("inf"),
                      cand_rss / base_rss, cand_cpu, base_cpu, cand_rss, base_rss, inner)


def cpu_ok(res, thresholds):
    """True if the candidate CPU time is within the relative threshold of the baseline."""
    return res.cand_cpu_median <= (1 + thresholds["cpu_rel_threshold"]) * res.base_cpu_median


def rss_allowed_increase_kb(base_rss_kb, thresholds):
    """Largest allowed peak-RSS increase: the relative threshold or the absolute floor."""
    return max(thresholds["rss_rel_threshold"] * base_rss_kb, thresholds["rss_abs_floor_kb"])


def rss_ok(res, thresholds):
    """True if the candidate peak RSS exceeds the baseline by no more than the allowed increase."""
    return res.cand_rss_max - res.base_rss_max <= rss_allowed_increase_kb(res.base_rss_max,
                                                                          thresholds)
