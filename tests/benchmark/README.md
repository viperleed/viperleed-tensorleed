# Large-system benchmark

`pytest tests -m benchmark --tl-root=<candidate> [--baseline-tl-root=<baseline>]` (opt-in, about 6 min for A/A: 336–375 s wall).

**ref-calc:** **Cu(111)-Te** (superstructure `5 0, 1 2`, viperleed example), case `tests/cases/benchmark/bench_cu111_te` (inputs only), generated with `python -m generate.make_refcalc_cases --benchmark` (from `tests/`) in an environment with viperleed installed (`THEO_ENERGIES = 100 400 10`, `LMAX = 8-14`; the 10 eV step matters because viperleed derives the phase-shift energy grid from it). Two single-energy runs at viperleed's LMAX for that energy:

| run | NT/LMAX | baseline user CPU | baseline peak RSS |
|---|---|---|---|
| `E100_L08` | 100 eV, LMAX 8 | 1.07–1.41 s | 113 MB |
| `E400_L14` | 400 eV, LMAX 14 | 25.4–30.2 s | 256 MB |

(8 cores / 15 GB, gfortran 15.2 `-O2`, system BLAS/LAPACK; release build of `3419e83` (ref-calc sources identical to v2.1.0). Measured 2026-10-07: CPU of `E100_L08` and the RSS values from one A/A run, CPU of `E400_L14` from 10 back-to-back runs, wall times from both.)

Choice: Cu(111)-Te is the benchmark system because, with the same settings, the 400 eV run (LMAX 10) of Fe₂O₃(012)-(1×1) has a smaller peak RSS and a shorter CPU time; Cu(111)-Te still stays < 4 GB and < 10 min CPU. Fe₂O₃(012) is defined in `systems.py` (`BENCHMARK_CANDIDATES`).

**rfactor:** case `tests/cases/benchmark/bench_ag100_rfactor`, of which only `case.toml` is committed: its inputs (~8 MB) are generated at test time, once per session, into the pytest basetemp by `build_benchmark_case` in `generate/rfactor_inputs.py` (standard library only, about 0.5 s), deterministically from the committed `ag100_rfactor` inputs (no viperleed run), and the case is loaded from there (`GENERATED` in `test_benchmark.py`). `python -m generate.make_rfactor_cases --benchmark` (from `tests/`) generates the case into a temporary directory, writes its `case.toml` and verifies the runs with the baseline. The case: 100 experimental and 120 theory beams (modulated copies of the four Ag(100) beams), 50 geometries, Pendry R factor. One run, `nb100_t120_ns50`: peak RSS 13.7 MB, user CPU 0.6 s (release build of the refactored tree, one run, 2026-10-23). The case exists for the memory gate: at the ~4 MB peak RSS of the tier-1 rfactor runs, the 1 MB floor lets an increase of ~25 % pass, here at most ~7.5 %.

Tests per run: (1) candidate outputs vs a fresh baseline run with the tier-1 tolerances (no stored expected outputs); (2) A/B CPU and peak RSS with the `[perf]` thresholds, `benchmark_repetitions` (3) interleaved repetitions.

A/A (baseline vs baseline, 2026-10-07, four runs): CPU ratio `E100_L08` 0.759–1.041, `E400_L14` 1.010–1.123 (one run failed at 1.123); RSS ratio ≤ 1.0005.

**The CPU gate is not reliable on the reference machine** (WSL2 on a hybrid CPU: the same executable's user CPU time is bimodal, about 25.5 s or 29.5 s for `E400_L14`), so it can fail an A/A run and pass a real regression. The RSS gate is not affected.
