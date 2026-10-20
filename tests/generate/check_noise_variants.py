"""Check that the configured tolerances pass all legitimate noise variants.

Builds the baseline with every variant in VARIANTS (compiler flags, BLAS
library, OpenBLAS kernel; all of them "the same program"), runs every tier-1
run, and compares the outputs with the expected outputs using the tolerances
of tests/tolerances.toml, exactly as tier 1 does. Every comparison must pass.

Also reports, per output-file pattern, the insensitivity of the configured
tolerances: the share of non-zero expected floats (as in calibrate_tolerances.py)
for which a 1 % change would pass.

Outputs are kept in --out/<variant>/<case>-<run>/ (default: a temporary
directory); --compare-only reuses them (e.g. after a tolerance change).

Usage (env tl-tests):
    python -m generate.check_noise_variants --tl-root=<baseline> [--out DIR]
        [--compare-only] [--variants NAME ...] [--cases CASE ...]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import argparse
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

from harness.build import BuildFlags
from harness.cases import (compare_outputs, execute_run, expected_outputs, iter_runs,
                           load_tolerances, pattern_for, tolerance_for)
from harness.compare import scaled_atols
from harness.tokens import tokenize

_LIB = "/usr/lib/x86_64-linux-gnu"
REFERENCE_BLAS_LINK = (f"{_LIB}/lapack/liblapack.so.3 {_LIB}/blas/libblas.so.3 "
                       f"-Wl,-rpath,{_LIB}/blas -Wl,-rpath,{_LIB}/lapack -lpthread -fbacktrace")
REASSOC = "-O3 -fassociative-math -fno-signed-zeros -fno-trapping-math -funroll-loops"

# name: (BuildFlags, extra environment)
VARIANTS = {
    "O0": (BuildFlags("O0", "-O0"), {}),
    "O1": (BuildFlags("O1", "-O1"), {}),
    "O3": (BuildFlags("O3", "-O3"), {}),
    "native": (BuildFlags("O2-native", "-O2 -march=native"), {}),
    "fma": (BuildFlags("O2-fma", "-O2 -mfma"), {}),
    "reassoc": (BuildFlags("reassoc", REASSOC), {}),
    "refblas": (BuildFlags("O2-refblas", "-O2", REFERENCE_BLAS_LINK), {}),
    "obl_sandy": (BuildFlags("release", "-O2"), {"OPENBLAS_CORETYPE": "Sandybridge"}),
    "obl_prescott": (BuildFlags("release", "-O2"), {"OPENBLAS_CORETYPE": "Prescott"}),
}
INSENSITIVE_REL = 0.01



def produce(tl_root, out, variants, cases):
    """Run all tier-1 runs with every variant (flags and environment) into out/<variant>/."""
    for name in variants:
        flags, env = VARIANTS[name]
        saved = {k: os.environ.get(k) for k in env}
        os.environ.update(env)
        try:
            for case, run in iter_runs(only=cases):
                wd = out / name / f"{case.name}-{run.name}"
                wd.mkdir(parents=True, exist_ok=True)
                try:
                    execute_run(case, run, tl_root, flags, wd)
                except Exception as err:  # noqa: BLE001
                    print(f"  {name} {case.name}/{run.name}: RUN FAILED "
                          f"{str(err).splitlines()[0]}")
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        print(f"produced {name}", flush=True)


def compare(out, variants, cases):
    """Number of failing runs over all variants."""
    n_fail_total = 0
    for name in variants:
        n_fail = n_runs = 0
        for case, run in iter_runs(only=cases):
            n_runs += 1
            wd = out / name / f"{case.name}-{run.name}"
            bad = [f"{name}: {msgs[0][:100]}" for name, msgs in
                   compare_outputs(wd, expected_outputs(case, run)).items()]
            if bad:
                n_fail += 1
                print(f"  FAIL {name} {case.name}/{run.name}: {bad[:3]}")
        print(f"== {name}: {n_fail}/{n_runs} runs fail")
        n_fail_total += n_fail
    return n_fail_total


def insensitive_counts(toks, tol):
    """(n non-zero floats in toks, n of them for which a 1 % change passes tol)."""
    n = n_ins = 0
    for t, scaled, exact in zip(toks, scaled_atols(toks, tol), tol.exact_mask(toks)):
        if t.kind == "float" and t.value != 0:
            n += 1
            n_ins += (not exact
                      and tol.allowed_abs(t.value, t.ulp, scaled)
                      >= INSENSITIVE_REL * abs(t.value))
    return n, n_ins


def expected_insensitive_counts(exp, tol):
    """insensitive_counts of the expected output file exp (a path), read with
    the ignore and skip settings configured for its name."""
    _, ignore, skip = tolerance_for(exp.name)
    return insensitive_counts(tokenize(exp.read_text(errors="replace"), ignore, skip), tol)


def insensitivity(cases):
    """{pattern: (n non-zero floats, n insensitive to a 1 % change)}."""
    patterns = list(load_tolerances()["numeric"])
    n, n_ins = Counter(), Counter()
    for case, run in iter_runs(only=cases):
        for exp in expected_outputs(case, run):
            pattern = pattern_for(exp.name, patterns)
            counts = expected_insensitive_counts(exp, tolerance_for(exp.name)[0])
            n[pattern] += counts[0]
            n_ins[pattern] += counts[1]
    return {p: (n[p], n_ins[p]) for p in n}


def main(argv=None):
    """Produce (unless --compare-only) and compare the variant outputs; report insensitivity."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tl-root", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--compare-only", action="store_true")
    parser.add_argument("--variants", nargs="*", default=list(VARIANTS))
    parser.add_argument("--cases", nargs="*", default=[])
    args = parser.parse_args(argv)
    with tempfile.TemporaryDirectory() as tmp:
        out = args.out or Path(tmp)
        if not args.compare_only:
            if args.tl_root is None:
                parser.error("--tl-root is required unless --compare-only")
            produce(args.tl_root, out, args.variants, args.cases)
        n_fail = compare(out, args.variants, args.cases)
    print("\ninsensitivity (share of non-zero values for which a 1 % change passes):")
    for pattern, (n, n_ins) in insensitivity(args.cases).items():
        print(f"  {pattern:10s} n = {n:7d}  insensitive {100 * n_ins / n:5.1f} %")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
