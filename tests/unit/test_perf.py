import sys, textwrap
from pathlib import Path

import pytest

from harness import cases, perf

PERF = cases.load_tolerances()["perf"]


def _stub(path, cpu_s, mb, log=None):
    path.write_text("#!" + sys.executable + "\n" + textwrap.dedent(f"""
        import resource
        if {str(log)!r} != "None":
            with open({str(log)!r}, "a") as f:
                f.write({path.stem!r} + "\\n")
        b = bytearray({mb}*1024*1024); b[::4096] = b'x'*len(b[::4096])
        while resource.getrusage(resource.RUSAGE_SELF).ru_utime < {cpu_s}:
            pass
    """))
    path.chmod(0o755)
    return path


def test_compare_performance_detects_slower_and_bigger_candidate(tmp_path, monkeypatch):
    log = tmp_path / "order.log"
    base = _stub(tmp_path / "base.py", 0.1, 10, log)
    cand = _stub(tmp_path / "cand.py", 0.4, 60, log)

    def fake_build(tl_root, program, param, flags):
        return cand if Path(tl_root).name == "cand" else base
    monkeypatch.setattr(cases, "build_program", fake_build)
    (tmp_path / "in").write_text("")
    case = cases.Case("c", tmp_path, "stub", [])
    run = cases.Run("r", tmp_path / "in", None, {}, [])
    res = perf.compare_performance(case, run, tmp_path / "cand", tmp_path / "base", reps=3,
                                   workroot=tmp_path / "w")
    assert res.cpu_ratio > 2.0
    assert res.cand_rss_max - res.base_rss_max > 25 * 1024
    assert res.base_cpu_median >= 0.1
    # sizing_runs baseline executions to size the samples, then interleaved
    assert log.read_text().split() == ["base"] * PERF["sizing_runs"] + ["base", "cand"] * 3


def test_rss_and_cpu_checks():
    th = PERF
    cpu = th["cpu_rel_threshold"]
    ok = perf.PerfResult(1 + cpu / 2, 1.0, 1 + cpu / 2, 1.0, 100_000, 100_000)
    assert perf.cpu_ok(ok, th) and perf.rss_ok(ok, th)
    slow = perf.PerfResult(1 + 2 * cpu, 1.0, 1 + 2 * cpu, 1.0, 100_000, 100_000)
    assert not perf.cpu_ok(slow, th)
    # large program: the relative threshold applies; small program: the absolute floor
    big = 100 * th["rss_abs_floor_kb"] / th["rss_rel_threshold"]
    small = th["rss_abs_floor_kb"] / th["rss_rel_threshold"] / 100
    for base, allowed in ((big, th["rss_rel_threshold"] * big), (small, th["rss_abs_floor_kb"])):
        assert perf.rss_allowed_increase_kb(base, th) == allowed
        assert perf.rss_ok(perf.PerfResult(1, 1, 1, 1, base + 0.9 * allowed, base), th)
        assert not perf.rss_ok(perf.PerfResult(1, 1, 1, 1, base + 1.1 * allowed, base), th)


def test_short_runs_are_repeated_until_a_sample_is_long_enough(tmp_path, monkeypatch):
    base = _stub(tmp_path / "base.py", 0.02, 10)
    monkeypatch.setattr(cases, "build_program", lambda *a, **k: base)
    (tmp_path / "in").write_text("")
    case = cases.Case("c", tmp_path, "stub", [])
    run = cases.Run("r", tmp_path / "in", None, {}, [])
    res = perf.compare_performance(case, run, tmp_path / "a", tmp_path / "b", reps=2,
                                   workroot=tmp_path / "w", sample_cpu_s=0.2)
    assert res.inner >= 2
    # inner is sized from an estimate: samples reach about sample_cpu_s
    assert res.base_cpu_median >= 0.16 and res.cand_cpu_median >= 0.16


def _sized_inner(tmp_path, monkeypatch, run_cpu_s, sample_cpu_s, **cfg):
    """inner chosen by compare_performance with [perf] keys overridden by cfg."""
    log = tmp_path / "order.log"
    base = _stub(tmp_path / "base.py", run_cpu_s, 1, log)
    monkeypatch.setattr(cases, "build_program", lambda *a, **k: base)
    monkeypatch.setattr(cases, "load_tolerances", lambda: {"perf": {**PERF, **cfg}})
    (tmp_path / "in").write_text("")
    case = cases.Case("c", tmp_path, "stub", [])
    run = cases.Run("r", tmp_path / "in", None, {}, [])
    res = perf.compare_performance(case, run, tmp_path / "a", tmp_path / "b", reps=1,
                                   workroot=tmp_path / "w", sample_cpu_s=sample_cpu_s)
    return res.inner, len(log.read_text().split())


@pytest.mark.parametrize("sizing_runs", [1, 2])
def test_sizing_runs_from_config(tmp_path, monkeypatch, sizing_runs):
    inner, n_exec = _sized_inner(tmp_path, monkeypatch, 0.0, 0.0, sizing_runs=sizing_runs)
    assert (inner, n_exec) == (1, sizing_runs + 2)


@pytest.mark.parametrize("max_inner", [2, 3])
def test_max_inner_from_config(tmp_path, monkeypatch, max_inner):
    inner, _ = _sized_inner(tmp_path, monkeypatch, 0.0, 100.0, sizing_runs=1,
                            max_inner=max_inner)
    assert inner == max_inner


@pytest.mark.parametrize("floor,expected", [(0.5, 2), (0.25, 4)])
def test_min_sizing_cpu_s_from_config(tmp_path, monkeypatch, floor, expected):
    # a run far below the floor: inner = ceil(sample_cpu_s / floor)
    inner, _ = _sized_inner(tmp_path, monkeypatch, 0.0, 1.0, sizing_runs=1,
                            min_sizing_cpu_s=floor)
    assert inner == expected
