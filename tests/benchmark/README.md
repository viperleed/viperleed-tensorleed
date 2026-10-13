# Large-system benchmark

`pytest tests -m benchmark --tl-root <candidate> [--baseline-tl-root <baseline>]` (opt-in, about 6 min for A/A: 336–375 s wall).

System: **Cu(111)-Te** (superstructure `5 0, 1 2`, viperleed example), case `tests/cases/benchmark/bench_cu111_te` (inputs only), generated with `python tests/generate/make_refcalc_cases.py --benchmark` in an environment with viperleed installed (`THEO_ENERGIES = 100 400 10`, `LMAX = 8-14`; the 10 eV step matters because viperleed derives the phase-shift energy grid from it). Two single-energy runs at viperleed's LMAX for that energy:

| run | NT/LMAX | baseline user CPU | baseline peak RSS |
|---|---|---|---|
| `E100_L08` | 100 eV, LMAX 8 | 1.07–1.41 s | 113 MB |
| `E400_L14` | 400 eV, LMAX 14 | 25.4–30.2 s | 256 MB |

(8 cores / 15 GB, gfortran 15.2 `-O2`, system BLAS/LAPACK; release build of `3419e83` (ref-calc sources identical to v2.1.0). Measured 2026-10-07: CPU of `E100_L08` and the RSS values from one A/A run, CPU of `E400_L14` from 10 back-to-back runs, wall times from both.)

Choice: Cu(111)-Te is the benchmark system because, with the same settings, the 400 eV run (LMAX 10) of Fe₂O₃(012)-(1×1) has a smaller peak RSS and a shorter CPU time; Cu(111)-Te still stays < 4 GB and < 10 min CPU. Fe₂O₃(012) is defined in `systems.py` (`BENCHMARK_CANDIDATES`).

Tests per run: (1) candidate outputs vs a fresh baseline run with the tier-1 tolerances (no stored expected outputs); (2) A/B CPU and peak RSS with the `[perf]` thresholds, `benchmark_repetitions` (3) interleaved repetitions.

A/A (baseline vs baseline, 2026-10-07, four runs): CPU ratio `E100_L08` 0.759–1.041, `E400_L14` 1.010–1.123 (one run failed at 1.123); RSS ratio ≤ 1.0005.

**The CPU gate is not reliable on the reference machine** (WSL2 on a hybrid CPU: the same executable's user CPU time is bimodal, about 25.5 s or 29.5 s for `E400_L14`), so it can fail an A/A run and pass a real regression. The RSS gate is not affected.
