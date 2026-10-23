# TensErLEED test suite

Regression tests for the TensErLEED programs `ref-calc` and `rfactor`, built for the F77 → modern Fortran refactor. They check correctness (outputs compared with a frozen baseline within calibrated tolerances), clean termination under a debug build, CPU time and peak memory (A/B against the baseline), and end-to-end behaviour through viperleed (R factors and error curves).

Baseline (ground truth): TensErLEED v2.1.0 (`22e3998`), which includes the ADREF2T, YSMOOTH and DER fixes (PRs #5, #7, #9); details in `cases/BASELINE.txt`.

## Quick start

```bash
conda env create -f tests/environment.yml          # creates env tl-tests: python 3.12, numpy, pytest, pytest-xdist, pytest-mock
# requirements: gfortran, gcc, BLAS/LAPACK (-llapack -lblas)

conda run -n tl-tests pytest tests -m tier1 --tl-root=<TL tree>    # unit + tier 1 (~20 s)
conda run -n tl-tests pytest tests -m debug --tl-root=<TL tree>    # debug build (~14 s)
conda run -n tl-tests pytest tests -m perf  --tl-root=<candidate> --baseline-tl-root=<baseline>   # A/B (~3 min)
conda run -n tl-tests pytest tests -m tier2 --tl-root=<TL tree>    # through viperleed (use -n 2, ~2.5 min)
conda run -n tl-tests pytest tests -m benchmark --tl-root=<candidate>   # large system (~6 min)
```

A "TL tree" is a directory containing `TensErLEED/src` and `TensErLEED/lib` (a checkout of this repository). `--tl-root` defaults to this repository. Machine-specific locations (baseline tree, and for the generators viperleed's example systems and a work directory) come from environment variables or `tests/local.toml` (git-ignored; copy `tests/local.toml.example`): `TL_TEST_BASELINE` / `baseline_tl_root` is the default of `--baseline-tl-root` and is needed by the perf and benchmark A/B tests and the build unit tests; a missing setting is a clear error naming both (`harness/paths.py`). Plain `pytest tests` runs everything except perf, debug, benchmark and tier2. Tier 2 additionally needs viperleed in the version pinned by `viperleed_version` in `tests/tolerances.toml` `[tier2]` (by default the command `viperleed` on `PATH`; `--viperleed-cmd`, `TL_TEST_VIPERLEED_CMD` or `viperleed_cmd` in `tests/local.toml` change it, see [Tier 2](#tier-2-viperleed)); it stops before any run on a mismatch.

All program runs use `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`.

## Parallel runs (`-n`, pytest-xdist)

`-n <N>` runs one suite in N worker processes (pytest-xdist); without `-n` it runs serially.

```bash
conda run -n tl-tests pytest tests -m tier1 -n 6 --tl-root=<TL tree>
conda run -n tl-tests pytest tests -m debug -n 6 --tl-root=<TL tree>
conda run -n tl-tests pytest tests -m tier2 -n 2 --tl-root=<TL tree>   # ~2.5 min
```

- **Limits** (for the reference machine, 8 cores and 15 GB): tier 1 / debug / unit up to `-n 6`; tier 2 at most `-n 2` (each viperleed run uses `N_CORES = 4`), and nothing else beside it. Independent suite runs started side by side count against the same budget.
- **Tier 2 refuses more than 2 workers.** With `-n 3` or more, a `-m` expression that can select tier 2 (`tier2`, `not tier1`, `-m ""`, …) is a usage error before any worker starts (exit code 4); as a backstop, a tier-2 test fails in setup on a worker before viperleed starts.
- **perf and benchmark run serially, alone on the machine,** because CPU times need a quiet machine. With `-n`, a `-m` expression that can select them (`perf`, `benchmark`, `not tier1`, `-m ""`, …) is a usage error before any worker starts (exit code 4); as a backstop, such a test fails in setup on a worker before it measures anything. Run them without `-n`, with `uptime` load < 1.
- **pytest internal API.** The early checks above evaluate the `-m` expression with `_pytest.mark.expression.Expression`, a pytest internal. If a pytest upgrade changes or removes it, the early check passes silently and the worker backstop decides: workers start, and every perf/benchmark test (any `-n`) or tier-2 test (`-n` > 2) fails in setup before it measures or starts viperleed, so timings and tier 2 never run beyond their limits. `tests/unit/test_xdist.py` then fails, which shows the change.
- **One run per (case, run).** The tests of a tier-1 run share its session-cached execution (`run_case`), and each worker has its own session; each (case, run) is therefore an `xdist_group`, and `-n` (default `--dist load`) is switched to `--dist loadgroup`, so its tests go to the same worker. Under `-n` the test ids carry the group as a suffix (`…[ag100/E045_L06]@ag100/E045_L06`).
- **Builds.** Workers share `tests/.build-cache`; a per-entry lock (`tests/.build-cache.locks/`, git-ignored) makes one worker build a program while the others wait for it, also across separate suite runs.
- **Temp dirs.** Each worker uses `<basetemp>/popen-gw<i>/`. Parallel suite runs side by side need their own `--basetemp`.

Checked 2026-10-06: `-m tier1`, `-m debug` (each from an empty build cache and with a warm one) and `-m tier2` give the same pass/fail/skip sets serially and with `-n`.

## Layout

```
tests/
  conftest.py, pytest.ini        options (--tl-root, --baseline-tl-root, --build, --viperleed-cmd), markers,
                                 pytest -n setup (loadgroup, perf/benchmark refused, tier 2 at most -n 2)
  tolerances.toml                numerical tolerances, tier-2 tolerances, perf thresholds
  local.toml.example             machine-specific paths (copy to local.toml, git-ignored)
  known_debug_failures.toml      baseline runs that fail the debug build (xfail strict)
  harness/                       tokens.py (Fortran output tokenizer), compare.py, tensor.py (T_* layout),
                                 runner.py (+ rusage_launcher.c), build.py, cases.py, paths.py, perf.py, tier2.py
  unit/                          tests of the harness (marked unit + tier1)
  tier1/                         test_refcalc.py, test_rfactor.py, test_array_temporaries.py,
                                 test_debug_build.py, test_performance.py, test_rfacsb_routines.py
                                 (+ rfacsb_routines/: driver.f90, values.json)
  tier2/                         test_viperleed.py, systems/<name>/ (frozen viperleed inputs + expected.toml)
  benchmark/                     test_benchmark.py, README.md
  cases/refcalc|rfactor|benchmark/<case>/   case.toml, inputs/, expected/<run>/ (benchmark: no expected/;
                                 bench_ag100_rfactor: inputs/ generated at test time, see benchmark/README.md)
  cases/MANIFEST.toml            all (case, run) pairs; tier 1 fails if the loaded cases differ
  generate/                      generators (run in an environment with viperleed installed unless noted; not run by pytest,
                                 except rfactor_inputs.py, which the benchmark tests use to generate inputs)
```

## Rules: tolerances, reference data, known debug failures

The suite checks that the refactored code reproduces the baseline. These rules keep the suite itself from drifting.

- **Tolerances** (`tolerances.toml` only, never inline in a test; derivation under [Tolerances](#tolerances)). They are calibrated on the baseline's own floating-point noise: the same sources built with gfortran 15.2 at `-O0`, `-O1`, `-O3`, `-O2 -march=native` and other variants, and linked against other BLAS/LAPACK builds, compared with the `-O2` reference. Each tolerance is set just above that noise, so that real changes are detected; loose tolerances would hide regressions. A tolerance may be loosened only with evidence that the deviation is pure floating-point reordering, for example that it disappears when both sides are built with `-O0 -ffp-contract=off`, or that a specific reordered sum explains it. Reordering during the refactor is expected, but it must be shown, not assumed. The PR that loosens a tolerance states the evidence. Exempt: the routine check (`tier1/test_rfacsb_routines.py`) requires the ported routines to agree bit for bit with the legacy routines on the same inputs, so it compares exactly and has no entry in `tolerances.toml`.
- **Reference data** (frozen inputs and `expected/` outputs under `cases/`, `cases/BASELINE.txt`, `tier2/systems/*/expected.toml`) is the ground truth. It changes only (a) for a deliberate bug fix or behaviour change, or (b) when a tolerance is shown to be too tight for legitimate floating-point reordering (see above). Inputs are frozen: only the expected outputs are regenerated (procedure under [Checking a TensErLEED fix](#checking-a-tenserleed-fix-pull-requests)); the regeneration commit contains only the new data, and the PR that changes reference data states why.
- **Known debug failures** (`known_debug_failures.toml`, see [Debug build](#debug-build)). A run may stay in the suite although the baseline debug build fails on it only if the failure is a bounds-check artefact of a declaration: a dummy array is declared smaller than the indices used, while the actual array is large enough, so no memory outside it is touched. The run must finish cleanly with every other check (`-fcheck=all,no-bounds -ffpe-trap=invalid,zero,overflow`), and its outputs from that build must equal the release outputs, i.e. the input shows no undefined behaviour. Inputs that trigger real out-of-bounds accesses or other undefined behaviour (e.g. more than one gap per beam) are not tested for exact values; they are documented under [Known baseline issues](#known-baseline-issues). To add an entry, give `case`, `run`, a `match` regex for the debug failure (on the message with the stderr tail, naming the source file and the array), a `reason` with the routine, the declaration and `file:line`, and a `spec_ref` pointing to the description in this README, and describe the issue under Known baseline issues. When a change fixes a listed failure, its test fails with "DID NOT RAISE"; that change removes the entry.
- **Debug opt-out** (`debug = false` in a run of `case.toml`, see [Case format](#case-format)). A tier-1 run whose input drives undefined behaviour that the baseline code shares (an out-of-bounds read, 0/0) may keep its release-build comparison and opt out of `-m debug`, where it is not collected. A comment above the flag names the legacy behaviour (routine, `file:line` or array), and the issue is described under [Known baseline issues](#known-baseline-issues). Used only where the run guards a path no other run reaches; a fix that removes the undefined behaviour removes the flag.

## Checking a TensErLEED fix (pull requests)

For a PR that changes `TensErLEED/`, on a clean clone (Python ≥ 3.11 with numpy, pytest, pytest-xdist, pytest-mock, e.g. the `tl-tests` env; gfortran, BLAS/LAPACK):

```bash
# baseline tree: the pin before the fix; create it before re-pinning (or set baseline_tl_root, see local.toml.example)
git worktree add --detach ../tl-baseline "$(sed -n 's/^commit: //p' tests/cases/BASELINE.txt)"
pytest tests -m tier1 -n 6 --baseline-tl-root=../tl-baseline   # unit + tier 1
pytest tests -m debug -n 6                                      # debug build
pytest tests -m tier2 -n 2      # if viperleed's results can change (needs viperleed, see Tier 2)
pytest tests -m perf --baseline-tl-root=../tl-baseline         # if performance can change, machine quiet
```

`--tl-root` defaults to the clone itself. The perf A/B runs against this pre-fix tree, so create it before the PR re-pins (below): afterwards `cases/BASELINE.txt` names the fix, and an A/B against that tree compares the fix with itself. On a PR that has already re-pinned, the pre-fix pin is the one at the merge base with master: `git show "$(git merge-base HEAD origin/master)":tests/cases/BASELINE.txt`.

If the fix changes outputs on purpose (see [Rules](#rules-tolerances-reference-data-known-debug-failures)), commit the fix and pin it as the new baseline (`BASELINE_COMMIT` in `generate/baseline_tree.py`; the generators refuse other commits). The inputs are frozen, so regenerate only the expected outputs, on the committed inputs (commands under [Regenerating](#regenerating-cases-and-expected-outputs)): `python -m generate.make_expected`, and for tier 2 `python -m generate.make_tier2_expected --skip-freeze` in an environment with viperleed installed, with a `tests/local.toml` (`baseline_tl_root` here a checkout of the fix with `eeasisss` built, not the pre-fix tree used for perf; `generate_workroot`). Then `git status` and `git diff --stat` must show only `generate/baseline_tree.py`, `cases/BASELINE.txt`, `expected/` files of the affected runs and `tier2/systems/*/expected.toml`; the PR states why each output changed. Any other change (inputs, `tolerances.toml`, `known_debug_failures.toml`) needs its own justification under the Rules. A rebase gives the fix a new hash: pin and regenerate again. Merge the PR with a merge commit, not squashed: CI checks out the pinned commit from `cases/BASELINE.txt`, and a squashed fix commit is no longer reachable once the PR branch is deleted.

CI (`.github/workflows/tests.yml`) runs unit + tier 1 and debug on pushes to `master` and `refactor/f90`, on pull requests and when started manually (Ubuntu 26.04, gfortran 15, reference BLAS/LAPACK). Tier 2, perf and benchmark do not run in CI; run them locally where they apply and report the result in the PR.

## Build contract

`harness/build.py` builds a program from a TL tree and caches the executable in `tests/.build-cache/<program>-<hash>/` (git-ignored). An executable from another tree state or other flags is never reused. The hash covers:
- the build mode and the flags;
- the full `gfortran --version` text and the resolved `libblas`/`liblapack` (`build.toolchain()`);
- the PARAM file, if used;
- the source files with their paths relative to `TensErLEED/` — every file under `TensErLEED/` in Makefile mode, the program's sources in legacy mode.

- **Makefile mode:** used whenever `TensErLEED/Makefile` exists, so it applies during the transition too, when modules and PARAM coexist. The harness copies `TensErLEED/` to a fresh directory and runs `make -B -C <copy> <refcalc|rfactor> BUILDDIR=<dir> FC=gfortran FFLAGS="<flags>" LDLIBS="<link flags>" [PARAM=<absolute path>]`. `-B` rebuilds everything: the copy keeps mtimes, so objects that an in-tree `make` left in `TensErLEED/` (other flags) would otherwise be relinked. The tree must build `<dir>/<refcalc|rfactor>` and must not depend on files outside `TensErLEED/`. `PARAM=` is passed (and hashed) as long as a main source of the program (`ref-calc.f`, `lib.tleed.f`, resp. `rfactor.f`, `rfacsb.f`) contains `INCLUDE 'PARAM'`. The Makefile must use it as the PARAM include for that build. The refactored tree has to provide this Makefile. Its rfactor link step uses no libraries and ignores `LDLIBS` (rfactor calls no BLAS or LAPACK routine). `tests/unit/test_build.py` checks the contract with a toy Makefile.
- **Legacy mode:** no `TensErLEED/Makefile`; the program's main sources include PARAM (the baseline). The harness compiles as viperleed does, with the run's PARAM file: `ref-calc.f` + `lib.tleed.f` (+ `GLOBAL`), resp. `rfacsb.f` + `rfactor.f`. It runs `gfortran <flags> -c …` and links with `-llapack -lpthread -lblas -fbacktrace`.

Source layout in the refactored tree: `refcalc` is still built from the fixed-form `src/ref-calc.f` + `lib/lib.tleed.f` with PARAM. `rfactor` is free form (`-std=f2018`) and needs no PARAM: `lib/rfacsb.f90` holds the modules in dependency order (kinds, array growth, output, interpolation, derivatives, R-factor terms, domain averaging, input of control parameters, experiment and theory, their preparation, the inner potential scan) and `src/rfactor.f90` the main program; all arrays are sized from the input. `make rfactor` builds it alone (`make rfactor PARAM=…` is accepted but the file is not used); viperleed compiles the same two files with `-c` each and links them. The header of `lib/rfacsb.f90` maps the routine and variable names of the fixed-form `rfactor.f`/`rfacsb.f` to the modules, procedures and derived-type components that hold them now.

Concurrent processes (`pytest -n`, suite runs side by side) build each cache entry once: the build holds a lock file `tests/.build-cache.locks/<program>-<hash>.lock` (git-ignored, never removed), and processes waiting for it use the finished executable.

Out of tree: `make` does not track flags, so objects left in the tree by a release build would be silently relinked into a "debug" executable; a fresh copy per build prevents this. In both modes the build fails if the hashed sources change during the build. A failed build keeps only its log (`.build-cache/failed-<program>-<hash>.log`); temporary build directories of processes that no longer run are removed.

Release flags: `-O2`. Debug flags: `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall` (warnings are logged, not errors; + link flags incl. `-fbacktrace`).

A run fails if the exit code is non-zero; stdout/stderr contain `Fortran runtime error`, `Program received signal`, `Error termination` or a `STOP <n≠0>` line; or the output is incomplete (ref-calc: `fd.out` lacks an energy of the EI,EF,DE range; rfactor: `ROUT` lacks the `AVERAGE R-FACTOR =` summary). The last check is needed because many legacy error paths print a message and `STOP` with exit code 0. Comparisons fail on missing files, a different number of lines or of tokens on a line, a number in a different format, an overflow field (`*`, `**…`; except the two `*` of an rfactor beam label `*( 1  0)*` and a `*` in column 1 before a letter or blank, which marks a comment line of rfactor's standard output) where the expected file has a number or an overflow field of other text or end column, and any differing integer or text token. An overflow field matches a textually identical one ending in the same column.

Routine check (`tier1/test_rfacsb_routines.py`): single routines of `lib/rfacsb.f90` (`der`, `smooth`, `ysmooth` of module `derivatives`) are called by a small driver (`tier1/rfacsb_routines/driver.f90`, built with gfortran at the release and the debug flags) on curves of 1, 2, 3 and a few more points, and their results must equal, bit for bit, those of the legacy routines on the same inputs (`tier1/rfacsb_routines/values.json`, reals as hex bit patterns). The values come from `generate/make_rfacsb_routine_values.py`, which builds `generate/rfacsb_routines_legacy_driver.f` against the baseline's `lib/rfacsb.f` at `-O2` and `-O0` and requires both to agree.

Static check (`tier1/test_array_temporaries.py`): each program is built with `-O2 -Warray-temporaries`; any warning fails (large arrays must not get compiler-generated temporaries; the baseline has none, and array temporaries are invisible to the peak-RSS gate unless they coincide with the peak).

## Case format

`cases/<refcalc|rfactor>/<case>/case.toml`:

```toml
program = "refcalc"
description = "..."
covers = ["bravais", "tensors"]          # informative
[[runs]]
name = "E045_L06"                        # expected outputs in expected/E045_L06/
stdin = "inputs/FIN_E045_L06"            # fed to the program on stdin
param = "inputs/PARAM_L06"               # legacy build mode only
files = { WEXPEL = "inputs/WEXPEL_x" }   # optional, copied into the run directory
outputs = ["fd.out", "amp.out", "T_*"]   # globs; resolved in expected/<run>/
debug = false                            # optional; not run in the debug build (Rules)
```

The output name `stdout` stands for the program's standard output, which the harness stores under that name in the run directory (`harness/cases.py`, `run_built`); every rfactor run lists it.

Every expected output file is compared line by line and token by token (same number of lines, same tokens per line, each number in the same format — F, E, ES or D, number of decimals, exponent width, column where the field ends; viperleed reads these files by fixed format) with the tolerance of the first matching pattern in `tolerances.toml`; the set of produced output files must equal the expected set.

## Regenerating cases and expected outputs

Expected outputs are the ground truth; when they may change is under [Rules](#rules-tolerances-reference-data-known-debug-failures), how under [Checking a TensErLEED fix](#checking-a-tenserleed-fix-pull-requests). The inputs of existing cases (FINs, PHASESHIFTS, tier-2 system inputs) are frozen: the viperleed generators below only create new cases, and the suite does not follow changes in viperleed's inputs (those are for viperleed's own tests). The commands with `python` run in an environment with viperleed installed, those with `conda run -n tl-tests` in the suite's environment.

Run the generators as modules from `tests/` (`cd tests`):

```bash
python -m generate.make_refcalc_cases [--only CASE] [--check-only]                          # FINs, PARAMs
python -m generate.make_rfactor_cases [--only CASE] [--benchmark] [--verify-only]           # WEXPEL, spectra (--benchmark: case.toml only)
conda run -n tl-tests  python -m generate.make_expected --tl-root=<baseline>                # expected/
conda run -n tl-tests  python -m generate.calibrate_tolerances --tl-root=<baseline> [--write]
python -m generate.make_tier2_expected [--only NAME] [--calibrate]
python -m generate.make_refcalc_cases --benchmark
conda run -n tl-tests  python -m generate.make_hand_cases [--only CASE] [--check]   # hand-made cases from frozen inputs
conda run -n tl-tests  python -m generate.make_rfacsb_routine_values --tl-root=<baseline>  # tier1/rfacsb_routines/values.json
```

The viperleed generators keep their work directories in `generate_workroot` and run viperleed with `-t <baseline tree>` (`baseline_tl_root`, see Quick start) and its own `eeasisss`. The tree must be at the baseline commit with an unmodified `TensErLEED/`, and `eeasisss` must have been built there with `make gcc GCC_FFLAGS="-O3 -fno-finite-math-only -std=legacy"` (gfortran ≥ 14 rejects a legacy FORMAT in eeasisss at run time otherwise); the generators check this. viperleed compiles in its own work directories and does not modify the tree. Regenerating the existing cases this way reproduces the committed inputs except for the time stamp in the phase-shift header line of the FIN, and except for `ag100_mix` (ref-calc case and tier-2 system): viperleed's phase shifts for its mixed Ag/Pd site vary from run to run. They write the fixed title `tl-test <case>` into all inputs, so outputs contain no directory names or timestamps and no output line needs to be skipped. `make_expected.py` refuses a tree whose HEAD is not the baseline unless `--allow-non-baseline` is given, and always a modified `TensErLEED/`; `BASELINE.txt` records the commit, not the path of the tree.

## Test systems

| case | code paths (checked in the FIN by `make_refcalc_cases.py`) |
|---|---|
| `ag100` | Bravais layers (NSUB = 1 everywhere), energy-dependent LMAX 6–9, tensor output for 5 layers, layer doubling |
| `ag100_oblique` | `BEAM_INCIDENCE = THETA 12, PHI 33` (non-zero THETA FI) |
| `ag100_bulk2` | `N_BULK_LAYERS = 2`, `BULK_REPEAT = [0 0 4.07294]`: top bulk layer type ≠ bottom bulk layer type, with different sites (lower bulk layer `Ag_bulkb`, vibration amplitude 0.10 Å instead of 0.069 Å), so the two layers differ numerically (`DBLGAS`) |
| `ag100_mix` | `ELEMENT_MIX Ag = Ag Pd`, surface site 0.7 Ag / 0.3 Pd (two non-zero concentrations, ATA) |
| `ir100_2x1_o` | composite layers (NSUB > 1), 2 beam sets, 2 elements, LMAX 6–8 |
| `ir100_2x1_o_thick` | Ir(100)-(2×1)-O with other layer cuts: O alone, then one layer of two Ir planes (4 sublayers with pairwise equal interplanar vectors: `NGEQ` copy in `SRTLAY`/`GHMAT`), and a bulk layer of two Ir planes (composite layer with bulk periodicity under a (2×1) overlayer: `ref-calc.f` `LAY = 2` branch, `RTINV_SIMPLE` substrate lattice sums) |
| `ir100_2x1_o_unsorted` | `ir100_2x1_o` with the sublayers of the composite top layer listed as Ir, Ir, O instead of O, Ir, Ir (out of order along +x: `SRTLAY` re-sorting, `NORD` mapping back to the input order for the tensor files); derived by `make_hand_cases.py`, 3 runs: E118 and E110-122 as listed, whose outputs equal those of `ir100_2x1_o` except the title line, and E118_three_heights with the Ir N=2 sublayer moved deeper, to x = 1.2004 (x grows into the crystal), and the sublayers listed as Ir N=2, O, Ir N=3 (three distinct heights, deepest first), which has its own outputs |
| `ag100_psformat2` | `ag100` with the phase shifts in the alternative format (`PSFORMAT = 2`, `READIN`: one line per element and energy, energies in eV divided by `HARTREE`; header `EM = 0` and `C1`–`C4` of `ag100`, so the inner potential is the same); derived by `make_hand_cases.py`, 2 runs (E045, E141). The energies are rounded to 0.001 eV, so `fd.out` equals that of `ag100` and `amp.out` is within its tolerance, while the `T_*` files differ by up to 2e-5 of their maximum (beyond tolerance) |
| `cu111` | Cu(111)-(1×1), hexagonal cell, ABC stacking, LMAX 5–7, 50–140 eV: the bulk repeat vector has a lateral shift whose mirror image is the other hollow, so in-plane sign errors (e.g. in `DBLGAS`) are not equivalent as on fcc(100). Defined in `generate/systems.py` (no viperleed example), beams from an `IVBEAMS` file |

Each ref-calc case made by `make_refcalc_cases.py` has 4 runs (the hand-made ones 2 or 3): lowest energy at the lowest LMAX, a middle energy at a middle LMAX, the highest energy at the highest LMAX (each as viperleed assigns LMAX per energy), and 4 consecutive energies at one LMAX.

rfactor cases (from the `ag100` data; every run also compares its standard output, `stdout`): `ag100_rfactor` — `pendry`, `r2`, `zj`, `rs` (all WR = 0 → R_S), `mixed` (R2 + Pendry), `smooth` (`ISMOTH = 2`), `averaging` (theory and experimental beams 3+4 averaged, groups 2 and 3 fractional), `v0_wide` (V0 range −6…+10); `ag100_rfactor_errors/mns` (7-geometry superpos spectrum of an error calculation); `ag100_rfactor_gaps` — `explicit_gap` (`GAP = 1`, one gap in beam 2, 100–120 eV) and `auto_gap` (beam 3 without points in 150–165 eV, detected by `EXPAV`). `ag100_rfactor_paths` (paths not reached by the other cases) — `theo_past_exp` (experimental data end at 190 eV, theory at 210 eV: the `COMNEI` branch for theory beyond the experiment), `permuted_weighted` (`IBP` = 2 1 3 0, so experimental beam 1 is compared with theory beam 2 and `MITTEL` differs between the two; theory beams 3+4 averaged with `SYM` codes 2 and 10, i.e. `RINTAV` weights 2 and 8), `r2_gap` and `zj_gap` (the explicit gap with R2 resp. ZJ weighting). The `pendry` run reproduces viperleed's own R factor for that system (0.2054 at V0 shift +3.00). `ag100_rfactor_short` — `pendry` and `rs` with the experimental beams cut to 2 (beam 1, from 120 eV), 23 (beam 2, from 85.5 eV), 22 (beam 3, from 150 eV) and 15 points (beam 4, from 170 eV): one curve for each `DER` branch (NE = 2, NE = 23 at the long/short threshold, the short formula). `ag100_rfactor_hand` (hand-made from the `ag100_rfactor` inputs by `generate/make_hand_cases.py`) — `pseudo_exp` (`EOT = 1`: the experiment is read by `READT` as a pseudo-experiment in theory-spectrum layout, geometry 1 of the `ag100_rfactor_errors` spectrum, layer 1 displaced by −0.06 Å; beams 3 and 4 start with zero intensities; R = 0.1316 at V0 shift +1.50) and `norm_max` (`NORM = 2`: `exp.column`/`theo.column` normalised to the highest intensity of each beam). `ag100_rfactor_params` (input parameters, one run each, hand-made from the `ag100_rfactor` pendry inputs by `generate/make_hand_cases.py`; the 2001 CPC paper describes these parameters only for the R-factor input of `search`, Table 7 and p. 415) — `theory_beam_unmatched` (experimental beam 3 removed, so theory beam 3 has no experimental counterpart: `KAV = 0` for it, `MITTEL` 1 1 1 0, `MNBED = 3`, as viperleed writes it; CPC: `KAV(k) = j` puts theory beam k into group j, whose beams are averaged; `KAV = 0` is not described), `weights_all` (`WR` = 0.1, 0.2, 0.3, 0.4 for R2, R_ZJ, R_P, R_S; not in CPC, where `search` uses either R_P or R2), `smooth_gap` (`explicit_gap` of `ag100_rfactor_gaps` with `ISMOTH = 2`; CPC: the experimental data are smoothed `ISMOTH` times with a three-point formula; gaps are not described), `all_fractional` (`MITTEL = 2` for every beam, so `NINT = 0`; CPC: `MITTEL` groups the beams into two sets for R-factor averaging, e.g. integer and fractional), `ipr2` (`IPR = 2` with `ISMOTH = 1`, standard output with the interpolated and smoothed curves and the Y functions; CPC: `IPR` = 0, 1, 2 sets the amount of standard output, for debugging) and `gaps3_all` (one explicit gap in each of experimental beams 1, 2 and 4, `MNGAP = 3`, all four weights: the gap offsets of later beams). `ag100_rfactor_beams` (hand-made from the `ag100_rfactor` pendry inputs: scaled copies of the 4 beams) — `nb25` and `nb30` (25 and 30 experimental and theory beams: `KAV`/`NBEA`/`IBP`/`MITTEL` fill one full 25I4 record, followed by a blank line, resp. two records).

Generator runs use `N_CORES = 2` (with 1 core viperleed uses a single LMAX).

Total size of the committed files in `tests/cases` + `tests/tier2`: 12.0 MB (all of `tests/`: 12.2 MB).

## Runtimes (8 cores, 15 GB, WSL2)

Measured 2026-10-07 (8 cores, WSL2, warm build cache, wall clock including `conda run`): unit + tier 1 ~20 s serial, ~11 s with `-n 3`; debug ~14 s serial, ~8 s with `-n 3`; tier 2 `-n 2` ~2.5 min; perf A/A ~3 min; benchmark A/A ~6 min. An empty build cache adds the builds (several per tree and flag set).

## Known baseline issues

Baseline = TensErLEED v2.1.0 (`22e3998`, with the ADREF2T, YSMOOTH and DER fixes). None of these are fixed by the test suite; they are listed so that the refactor preserves or deliberately fixes them. File names and line numbers refer to the baseline sources (fixed-form `rfactor.f`, `rfacsb.f`); each entry says what the refactored tree does.

### ref-calc

- **`OUTXIST`, debug build (bounds check).** `OUTXIST` declares `SPQF(2,NT)` but indexes it up to `KNT` (`lib.tleed.f:3003`). Ten tier-1 runs (`ir100_2x1_o` 3, `ir100_2x1_o_thick` 3, `ir100_2x1_o_unsorted` 3, `cu111` 1) fail `-fcheck=bounds` there. Declaration artefact; the runs are clean without bounds checks.

### rfactor

- **Gap input, debug build (bounds check).** Every input with `GAP = 1` fails `-fcheck=bounds` in `GAPIN` (`rfacsb.f:598`): `BM(NGAP)` is declared with the *output* argument `NGAP`, which is 0 on entry. The actual array (`GBM(MNGAP)`) is large enough, so this is a declaration artefact, not an out-of-bounds access in memory. Refactored tree: gone; `gapin` sizes the gap arrays from the number of gaps read.
- **Auto-detected gap, debug build (bounds check).** Every input with a gap detected by `EXPAV` fails `-fcheck=bounds` in `GRID` (`rfacsb.f:631`): `A(NGP)` is declared with `NGP`, which is undefined in the caller (`rfactor.f:603`). Only the declaration uses `NGP`; the actual array `A(MNGP)` is large enough. Refactored tree: gone; the caller defines the working grid size `ngp` with which `grid` declares its arrays.
- **Zanazzi-Jona, debug build (bounds check).** `WR(2) > 0` fails `-fcheck=bounds` in `SUM` (`rfacsb.f:1216`, called from `VARSUM`, `rfacsb.f:1294`): `SUM` declares `Y(1,1,NGP)` but receives `YY(10*NGP)` with up to `10*(N-1)+1` points. Declaration artefact; memory access is within `YY`. Refactored tree: gone; `SUM` is the private `trapezoid` in module `rfactor_terms` and takes the actual number of points.
- **Gaps, all builds.** In the per-gap loops of `rfactor.f`, the gap integrals (`SG`, `SG2`, `SGZJ`, `SGYP2`/`SGYS2`) are overwritten instead of summed, which is wrong for more than one gap per beam; `EXPAV` writes auto-detected gaps into arrays of size `MNGAP` (1 as viperleed writes it) without a bounds check. Refactored tree: the overwriting is kept (bug fixes are out of scope); the gap arrays grow to the number of gaps, so more gaps no longer write out of bounds; several gaps in one beam stay untested.
- **Gaps with an experimental data step ≠ `EINCR`, all builds.** With a gap, `EE` stays on the raw data grid, so `COMNEI` gets wrong indices: garbage overlaps and R factors, out-of-bounds reads. Refactored tree: same behaviour (no input check: non-uniform data with gaps is valid and works); not tested.
- **Gap outside the common interval, all builds.** For a gap entirely outside the energy interval common to theory and experiment, the theory gap integral reads one element before the array (`rfactor.f:1078`, `AT(IBT,0)`). The integral then runs over a single point (N = 1), so the value read does not change the result. Refactored tree: the same read (`varsum`, `a1(0)`). Tier-1 run `ag100_rfactor_params/autogap_edge_r2_v0std` covers it (a detected gap at the upper end of the data); it opts out of the debug build, which stops there (see [Rules](#rules-tolerances-reference-data-known-debug-failures)).
- **No overlap of theory and experiment, all builds.** An input in which no beam overlaps its theory beam divides 0 by 0 in the average over the beams (`average_over_beams` in the refactored tree, likewise in the baseline); the release build writes `*******`, the debug build stops with SIGFPE. Tier-1 run `ag100_rfactor_params/no_overlap` covers it; `ROUT` and stdout hold the overflow fields, which the comparison matches textually; the run opts out of the debug build.
- **Skipped geometries (`NSSK`), all builds; dropped: rfactor rejects any `NSSK` ≠ 0 with an error.** `RDPL` reads the value of each of the NS geometries into `SP(J)`, J = 1 … NS, but declares the automatic array `SP(NSS)` with the number of geometries left after skipping (`rfacsb.f:787`, written at `rfacsb.f:815`). Every input that skips a geometry (NSS < NS) writes past the end of `SP`; the release build crashed in `free` (segmentation fault) on the `ag100_rfactor_errors` spectrum with `NSSK = 2,4,7`. Refactored tree: the error stops before `RDPL`; `tier1/test_rfactor_rejected_inputs.py` tests it.

Gap test inputs (`cases/rfactor/ag100_rfactor_gaps`, the gap runs of `ag100_rfactor_paths` and `ag100_rfactor_params`) have at most one gap per beam (one in total, except `gaps3_all`), so they avoid the gap bugs that depend on the number of gaps (see above). They run clean with `-O0 -g -fcheck=all,no-bounds -ffpe-trap=invalid,zero,overflow`, and their outputs are identical between that build and the release build (checked for `ag100_rfactor_params/smooth_gap` and `weights_all` on 2026-10-22).

## Tolerances

`tests/tolerances.toml` holds one tolerance per output-file pattern. A float `a` matches the expected `b` if

    |a-b| <= rtol*|b| + atol_frac*S + atol + ulps*ulp(b)

where `ulp(b)` is the unit in the last printed digit of `b` (e.g. 1e-7 for `0.16631E-02`) and `S` is the largest |float| of the expected file (`scale = "file"`) or of `b`'s line (`scale = "line"`). Integers and text must match exactly; line and token counts and number formats must match. The exponent letter must match (E vs D is a difference), and a mantissa with a non-zero integer digit (`-9.79E-02`, as `ES` or `1PE` write it) differs from `-0.979E-01` (a zero value counts as plain E). The exponent width counts letter, sign and digits, so `0.1234-100` (Ew.d for |exponent| > 99) has the same format as `0.1234E-99`; integer digits are not compared (`-.1234` = `-0.1234`).

Calibrated with `tests/generate/calibrate_tolerances.py` (2026-10-05, gfortran 15.2, system BLAS/LAPACK): the baseline is built with `-O0`, `-O1`, `-O3`, `-O2 -march=native` and `-O2 -fno-fast-math`, all tier-1 runs are executed, and the outputs are compared with the expected outputs (`-O2`). ifort is not available on the calibration machine.

Findings:
- `-O0`, `-O1`, `-O3` and `-O2 -fno-fast-math` give outputs **bit-identical** to `-O2`.
- `-O2 -march=native` differs. The cause is only FMA contraction: `-O2 -march=native -ffp-contract=off` is bit-identical to `-O2`, and `-O2 -mfma` reproduces all deviations. Largest deviations: `fd.out` 6 units in the last printed digit; `*.column` 1 unit; `amp.out` 1.1e-6 absolute (relative up to 6e-5 for amplitudes ~1e-2); `T_*` 1.5e-4 absolute (1.9e-5 of the file maximum; single-precision noise, including near-zero components that change sign); `ROUT`/`ROUTSHORT` none.

Rule: for each pattern, candidate tolerances that cover the observed deviations are built, and the one with the lowest *insensitivity* is chosen. Insensitivity is the share of non-zero expected values for which a 1 % change would pass.
- **Print-level noise** (≤ 10 ulps): `ulps = max(2, 3 × observed)`. Rounding flips are quantised, so a factor of 3 is used instead of 10.
- **File-scale absolute part:** `atol_frac = 10 × max |a-b|/S_file`, plus `ulps = 2` and `rtol = 10 ×` the relative deviation not covered by these (floor 1e-6).
- **Line-scale absolute part:** the same with `S_line`.
- **`ROUT`/`ROUTSHORT`:** fixed at `atol = 1.5e-4`, so one unit in the 4th decimal of an R factor passes and two fail. All other numbers in these files (V0 shifts, energies, overlaps; 2 decimals) must match exactly.
- **rfactor `stdout`** (the per-beam R_2/R_ZJ/R_P/R_S per V0, the normalisation constants, CIFE/CIFT, the echoed input): fixed at `ulps = 1.5`, one unit in the last printed digit of each value passes, two fail, as for `ROUT`. Noise (`check_noise_variants.py`, rfactor cases): `-O0`, `-O1`, `-O3`, reference BLAS and the OpenBLAS kernels are bit-identical; FMA (`-march=native`, `-mfma`) and reassociation flip values by 1 unit, except in `ag100_rfactor_short` (up to 30 units with FMA, 32 with reassociation, on its short beams), where `ROUT` deviates by the same amount. The output contains no paths, timings or dates. Insensitivity 0 %.

Result of that calibration (insensitivity in brackets): `fd.out` 18 ulps (8.4 %), `*.column` 3 ulps (0.3 %), `amp.out` file scale (13.5 %), `T_*` file scale, atol_frac 1.9e-4 (19.9 %), `ROUT` (17 %).

**Structural rules.** The file-scale rules for `T_*` and `amp.out` are set by the noisiest values of the file and check the rest 100–10 000× more loosely than their noise. These patterns use structural rules instead; each is capped at the file-scale rule (`file_cap`), so no value is checked more loosely than under it:
- **`T_*`** (`scale = "tensor"`): `harness/tensor.py` parses the OUTXIST/OUTAMP layout (header line, L1, t-matrix CAF, XIST, then per exit beam NEXIT, k-vector line, ALM, ending with `-1`; any deviation fails the comparison). The absolute part is `sections[s]` × max |value| of the value's WRITE block, with hdr 0, kvec 0, tmat 1.2e-6, xist 1.4e-4, alm 3e-4 (10× the largest noise of any variant, FMA, relative to the block maximum).
- **`amp.out`** (`scale = "line"`): 2e-4 × the largest amplitude of the line.
- **Exact columns:** F-format floats with ≤ 5 decimals in `fd.out` and `amp.out` (energy, AP, beam list) and the energy columns of `*.column` (`exact_columns = ["odd"]`); all are bit-identical in every noise variant.

Noise variants (`generate/check_noise_variants.py`, all tier-1 runs): `-O0`, `-O1`, `-O3`, `-O2 -march=native`, `-O2 -mfma`, reassociation (`-O3 -fassociative-math -fno-signed-zeros -fno-trapping-math -funroll-loops`), reference BLAS/LAPACK, OpenBLAS with `OPENBLAS_CORETYPE=Sandybridge` and `=Prescott`. Every run passes every variant, except `ag100_rfactor_short` under FMA (`-O2 -march=native`, `-O2 -mfma`) and reassociation: its `ROUT` and `stdout` differ by up to ~32 units in the last digit on its short beams. This exception is deliberate. The `ROUT`/`stdout` tolerances stay at one unit, because a looser one would nearly pass a 1 % change of an R-factor, and `ag100_rfactor_short` is the only run that detects a precision change in a single expression: one interpolation expression computed in double precision, or the end correction of the trapezoid sum regrouped, fails its two runs and passes all others. A build-flag change that adds FMA or reassociation is therefore checked on every run except this case. The machine's OpenBLAS kernel changes thousands of `T_*` values, all within the tolerances. `calibrate_tolerances.py` keeps the structural rules unchanged.

Insensitivity with these rules (non-zero values; exact values count as sensitive): `T_*` 12.9 % (file-scale rule: 19.9 %), `amp.out` amplitudes 4.8 % (13.6 %), `fd.out` 0 % (without exact columns: 8.4 %), `*.column` 0.3 %, `ROUT` 17 %. Detected by them and passed by the file-scale rules: a 1e-4 error of one t-matrix element, a 1 % error of all XIST amplitudes or 3 % of one, an energy 0.1 eV off in `fd.out`. The remaining `T_*` insensitivity is ALM components far below their block maximum (single-precision noise).

When a tolerance may be loosened: [Rules](#rules-tolerances-reference-data-known-debug-failures).

## Known coverage gaps

Line coverage of all tier-1 runs on the baseline (`python -m generate.coverage_report --missed 2`, gcov, 2026-10-06): `ref-calc.f` 96.3 %, `lib.tleed.f` 90.5 %, `rfactor.f` 90.3 %, `rfacsb.f` 83.6 % of executable lines. Never executed (percentages measured at `3419e83`; line numbers refer to v2.1.0):

- **Dead code** (cannot be reached by any input): `GHMAT` 2191-2200 and `SRTLAY` 4618-4640, the reuse of interplanar propagators of an *earlier* call (`NGOL`), because `NEW` is always 1 (`ref-calc.f:997`); hence a deliberately injected bug in that copy cannot be detected. `ZSU`, `DEBWAL`/`Debwal` (unused).
- **Not produced by viperleed inputs:** unformatted tensor output in `OPENOUT`, `OUTAMP`, `OUTXIST` and `DELIMIT` (`IFORM = 0`; viperleed always writes `IFORM = 1`). The other paths of this kind in the measurement (`SRTLAY` 4561-4578, `READIN` `PSFORMAT = 2`, `READT`, `COLUMN` with `NORM = 2`) are now covered by hand-made cases (below; for `SRTLAY` only in part).
- **Numerical branches not hit by the test systems:** pivot row swaps in `ZGE` (5473-5475) and `CXMTXT` (933-958); single-beam path of `SUBRAS_BLAS` (4741-4747); zero vibration in `PSTEMP`; `GHD` 2026-2027, `TSCATF` 5287-5288, `sphrm` 4437-4440; `ref-calc.f` 738-742 (alternative energy-dependent inner potential) and 1145-1147.
- **rfactor:** a gap fully outside the range (665-670) or at the ends of the common interval (937-969), no overlap of theory and experiment (981-985); `RINTAV` weights 3, 6, 12 (988-992); `DER` for NE < 23 (`rfacsb.f` 361-377; N = 1 branch added by the DER fix, not measured); `XNTRP2`; `YPEND` negative intensity; `EXPAV` 566-568, `GAPIN` 600-605, `GRID` 636-638, `INTPOL` 692-694. The dropped input paths (`IBP = 0` among the first NBE entries, `NBEA = 0`, `LIMFIL` ≠ 1, `NSSK` ≠ 0, `V0RR` ≠ 0) stop with an error, tested in `tier1/test_rfactor_rejected_inputs.py`.
- **Refactored rfactor** (`lib/rfacsb.f90`, gcov over all tier-1 rfactor runs, 92 % of 1405 lines): `write_no_overlap` and the `beam_rfactor` branches for a gap at the start of or outside the common interval (inputs not in tier 1, see Known baseline issues); `write_too_few_points`; `write_suspect_intensity` (stop on an averaged theory intensity > 1); the `prepare_experiment` branch for a gap that ends before the first non-zero intensity; `write_smoothed` on the gap path (`IPR = 2` with smoothing and a gap); `rintav` weights 3, 4, 6 and 12; `xntrp2`; the early returns of `grid` and `intpol` for a curve that is zero up to its last point and of `grid_points` for fewer than two points; `der` for one point; `smooth` with non-uniform steps; `binsrx` below the first point; the error stops (tested in `tier1/test_rfactor_rejected_inputs.py`, not in this measurement). Gaps with a data step ≠ `EINCR` are not tested (see Known baseline issues).
- **Error exits** (`WRITE … STOP` on invalid input or failed convergence): most of the remaining short blocks, e.g. in `READGEO`, `READIN`, `SUBRAS_BLAS` 4764-4768, `rfactor.f` 452-453 and 769-772.

Covered by cases not in the measurement above: composite bulk layer (`ref-calc.f` 1058-1079, `RTINV_SIMPLE` substrate branch) and the `NGEQ` copy in `GHMAT` (`ir100_2x1_o_thick`); `COMNEI` for theory beyond the experiment, permuted `IBP` with `MITTEL` differing, `RINTAV` weights 2 and 8, R2 and ZJ inside gaps (`ag100_rfactor_paths`).

Covered by hand-made cases (`generate/make_hand_cases.py`, not in the measurement above; checked with `python -m generate.coverage_report --lines`): `SRTLAY` 4561-4578 (`ir100_2x1_o_unsorted`: sublayers out of order along +x, re-sorted with one move when the two Ir sublayers at equal height are listed before O, and with two moves when three distinct heights are listed deepest, topmost, middle; the order of sublayers at equal height, i.e. whether the sort is stable (`.GE.` vs `.GT.`), does not show in the outputs and is not tested), `READIN` `PSFORMAT = 2` 3674-3703 except its two error stops (`ag100_psformat2`), `READT` except the beam-count error stop (its `ELSE IF (NEE(K).EQ.0)` guard is never false: the pseudo-experiment has zeros only before a beam's first non-zero point, as real theory output does, so zeros after the onset are untested) and `COLUMN` with `NORM = 2` (`ag100_rfactor_hand`). Covered by `ag100_rfactor_params`, each with a mutant inside the branch that fails only its run: `KAV = 0` (`RINTAV` leaves the beam out of every group), all four `WR` weights in the average R factor, smoothing with gaps (`rfactor.f` 560-591), `NINT = 0` (no integer-beam averages) and the `IPR = 2` prints. Also in `ag100_rfactor_params`, without mutants of their own: a detected gap outside the common interval for some inner-potential shifts (`autogap_edge_r2_v0std`) and no overlap of theory and experiment (`no_overlap`); see Known baseline issues.

Covered by `ag100_rfactor_short` (not in the measurement above): every `DER` branch except NE = 1, i.e. NE = 2 (`rfacsb.f` 374-377), 3 ≤ NE ≤ 22 (361-368) and NE = 23 at the threshold (345). The NE = 1 branch (370-372) cannot be reached through rfactor: a 1-point beam fails the bounds check in `COMNEI` (`rfacsb.f:317`) first.

## Debug build

`pytest tests -m debug` builds with `-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall -fbacktrace` and requires every run except those with `debug = false` (see [Rules](#rules-tolerances-reference-data-known-debug-failures)) to terminate cleanly (outputs compared with `rtol >= min_rtol` from `[debug]` in `tolerances.toml` (1e-3) as a sanity check).

Known baseline failures are listed in `tests/known_debug_failures.toml`. For each one the test requires:
1. the debug build to fail with exactly the known error (regex `match` on the message with the stderr tail). When the refactor fixes it, the test fails ("DID NOT RAISE") and the list has to be updated deliberately;
2. the same run to finish cleanly with all checks except bounds (`-fcheck=all,no-bounds -ffpe-trap=…`). So new FP traps or other check failures in these paths are still caught: an injected division by zero in `GAPIN` fails this check (a plain `xfail` would pass it).

Refactored tree: rfactor has no known debug failures (`SUM`, `GAPIN` and `GRID` artefacts gone, see Known baseline issues), so `known_debug_failures.toml` lists only the 10 `OUTXIST` entries; the rfactor runs that the baseline debug build stopped on (ZJ, all `WR`, the gap runs) finish cleanly in it.

Status of the 60 tier-1 runs of the debug build (all 62 except the 2 with `debug = false`): refactored tree 50 clean, 10 known failures (`OUTXIST`). Baseline 42 clean, 18 known failures: 10 × `OUTXIST` SPQF bounds (`ir100_2x1_o` 3, `ir100_2x1_o_thick` 3, `ir100_2x1_o_unsorted` 3, `cu111` 1), 2 × `SUM` (`ag100_rfactor/zj`, `ag100_rfactor_params/weights_all`), 5 × `GAPIN` (`explicit_gap`, `r2_gap`, `zj_gap`, `smooth_gap`, `gaps3_all`) and 1 × `GRID` (`auto_gap`), all bounds-only declaration artefacts (see Known baseline issues); all 18 run clean without bounds checks, the 8 rfactor runs with outputs identical to the release build. The 10 `SPQF` entries stay until a change corrects the declaration (see [Rules](#rules-tolerances-reference-data-known-debug-failures)).

## Performance (A/B)

`pytest tests -m perf --tl-root=<candidate> [--baseline-tl-root=<baseline>]` (~3 min) builds both trees (release) and runs every tier-1 run in `repetitions` (5) samples each, interleaved (baseline, candidate, baseline, …), in fresh directories. A sample repeats a short run `inner` times, so that it takes at least `sample_cpu_s` (0.1 s) of user CPU and the samples of all runs of a case together take at least `case_sample_cpu_s` (0.6 s). `inner` is sized from the median of `sizing_runs` (3) baseline executions (floored at `min_sizing_cpu_s`, 1 ms) and capped at `max_inner` (500). All thresholds are relative to the baseline measured in the same session, so the tests do not depend on the machine's speed; they do depend on the machine being quiet. Performance results are authoritative only on the reference machine (the 8-core WSL2 workstation the thresholds were calibrated on); on other machines and in CI they are advisory.

Thresholds (`[perf]` in `tolerances.toml`; memory should essentially not increase, CPU time is noisier):
- **Peak RSS** per run: increase ≤ max(2 %, 1 MB).
- **CPU time** per case: sum of the per-run median sample times of the case ≤ 1.10 × baseline.

Per-run CPU comparisons are too noisy: these runs take milliseconds, and A/A deviations reach 12 % per run even with 0.3 s samples. A sum per program dilutes slowdowns of one code path: a +30 % slowdown of the composite-layer runs only gives a ref-calc sum of 1.084 and passes it; per case it fails (`ir100_2x1_o` 1.24). The benchmark (`-m benchmark`) tests CPU per run for two large ref-calc runs and one large rfactor run. A table is printed in the terminal summary.

Peak RSS is measured through `harness/rusage_launcher.c` (built with `gcc` on first use, cached per source and `gcc --version` text): on Linux `ru_maxrss` survives `execve`, so a program started directly from Python reports at least Python's own peak RSS (~35 MB here; ref-calc needs 8–13 MB, rfactor ~4 MB).

Observed A/A noise (baseline vs baseline, 3 sessions on 8 cores / 15 GB, WSL2, 2026-10-06): CPU per case within ±4 % (0.972–1.038); RSS within ±140 kB (far below the 1 MB floor). With only `sample_cpu_s` (no per-case minimum), the single-run case `ag100_rfactor_errors` (0.08 s per sample) reached 1.16 in A/A.

Sensitivity (manual, not committed): a 2 MB static array touched once fails the RSS test for all 20 ref-calc runs. Doubling the energy-independent `CELMG` or the per-energy `FMAT` call changes CPU by < 2 % (both are cheap) and is not detected.

## Tier 2 (viperleed)

`pytest tests -m tier2 --tl-root=<TL tree> [--viperleed-cmd "<cmd>"]` (default command: `viperleed` on `PATH`; to run another one, e.g. `conda run --no-capture-output -n <env> viperleed`, set it with `--viperleed-cmd`, the environment variable `TL_TEST_VIPERLEED_CMD` or `viperleed_cmd` in `tests/local.toml`, which take precedence in this order; the viperleed version is pinned in `tests/tolerances.toml`, `[tier2] viperleed_version`; the checksum rule is below). For each system in `tests/tier2/systems/` it runs `viperleed calc -t <tl_root>` with `RUN = 1 5` (reference calculation + error calculation using delta, superpos and rfactor of the same tree) and compares the refcalc R factor, the best V0 shift and every point of every `Errors_*.csv` in `OUT/Errors.zip` with `expected.toml`.

**Checksum check.** The frozen `PARAMETERS` of the systems contain `TL_IGNORE_CHECKSUM = True`; the harness overrides it in the work directory for each run. If `--tl-root` is a released TL version, viperleed's checksum check stays active (`TL_IGNORE_CHECKSUM = False`), so tier 2 also verifies that the tree matches the released sources. For any other tree, `TL_IGNORE_CHECKSUM = True` is set, because viperleed rejects every changed file and every unknown TL version. "Released" means: `--tl-root` is the top level of its own git checkout (a git worktree counts), whose `HEAD` carries a tag `vMAJOR.MINOR.PATCH` and has no modified tracked files (`harness/tier2.py`, `is_released_tl_tree`). A tree without a git checkout of its own counts as not released, including a folder inside another checkout. The tags must have been fetched: a clone made with `--no-tags` or a shallow clone without tags counts as not released (`git fetch --tags` fixes that).

| system | error curves |
|---|---|
| `ag100` | z of layer 1, z of layer 2 (±0.06 Å, step 0.02); vibration `Ag_surf` (±0.02 Å, step 0.01) |
| `ag100_mix` | occupation of `Ag_surf`: Ag 0.5→0.9, Pd 0.5→0.1 (step 0.1) |
| `ir100_2x1_o` | z of `O_ads`, z of `Ir_surf` (±0.06 Å, step 0.02); vibration `O_ads` (±0.02 Å, step 0.01) |

Inputs are frozen after viperleed's initialization (annotated POSCAR, PHASESHIFTS, IVBEAMS, VIBROCC), so the tree under test needs no eeasisss. Regenerate the expected values with `python -m generate.make_tier2_expected --skip-freeze` in an environment with viperleed installed (see [Rules](#rules-tolerances-reference-data-known-debug-failures)). The expected values, first generated at `3419e83`, are unchanged at v2.1.0: no tier-2 curve has N ≤ 22 (smallest N = 129), so the `DER` fix does not act on them.

Tolerances (`[tier2]` in `tolerances.toml`): R values `r_atol = 1.5e-4` (one unit in the 4th decimal), displacements/V0 exact. Calibration (`python -m generate.make_tier2_expected --calibrate`): `FORTRAN_COMP = 'gfortran -O0'` identical to `-O2` for all systems; `-O2 -march=native` changes one R value of `ir100_2x1_o` by 1e-4.

Runtime: ~2.5 min for the three systems with `-n 2` (N_CORES = 4).

Negative check: with the ADREF2T fix reverted, all three systems fail in their error curves (`ag100`: both geometry curves and the vibration curve; `ag100_mix`: occupation; `ir100_2x1_o`: all curves); refcalc R factors are unaffected, as expected.
