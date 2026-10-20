"""Calibrate the numerical tolerances of tests/tolerances.toml.

Builds the baseline with compiler-flag variants that must count as "the same
program" (-O0, -O1, -O3, -O2 -march=native, -O2 -fno-fast-math), runs all
tier-1 runs with each, and compares the outputs with the stored expected
outputs (release build, -O2). For every float it records the deviation |a-b|,
the expected value's unit in the last printed digit (ulp) and the largest
|float| of its line.

Rule per output-file pattern (see tests/README.md, 'Tolerances'):
- Print-level noise (every deviation <= MAX_PRINT_ULPS ulps): tolerance in
  units of each value's own last printed digit,
  ulps = max(ULPS_FLOOR, ULPS_SAFETY * max observed ulps); rtol = atol_frac = 0.
  The tolerance then scales with each value, not with the largest number in
  the file (in fd.out that is the energy, which made a per-file absolute
  tolerance blind to intensity changes).
- Larger deviations (single-precision cancellation in small tensor
  components): an absolute part relative to the largest |value| of the same
  line, atol_frac = SAFETY * max(|a-b| / line max), scale = "line", plus
  ulps = ULPS_FLOOR; rtol = SAFETY * relative deviation not covered by that,
  floor 1e-6.
- ROUT/ROUTSHORT (R factors with 4 decimals): atol = 1.5e-4, i.e. one unit in
  the last digit of an R factor passes, two fail; V0 shifts, energies and
  overlaps (2 decimals) must match exactly.

For each pattern the script also reports the insensitivity of the result:
the share of non-zero expected floats whose tolerance allows a change of 1 %
or more (such changes go undetected).

ifort/ifx is not available on the calibration machine; skipped.

Usage (env tl-tests):
    python -m generate.calibrate_tolerances --tl-root=<baseline> [--write]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import argparse
import datetime
import sys
import tempfile
from pathlib import Path

from harness.build import BuildFlags
from harness.cases import (TOLERANCES_FILE, execute_run, expected_outputs, pattern_for,
                           iter_runs, load_tolerances)
from harness.compare import Tolerance, float_scales
from harness.tokens import tokenize
from .check_noise_variants import expected_insensitive_counts

VARIANTS = [BuildFlags("O0", "-O0"), BuildFlags("O1", "-O1"), BuildFlags("O3", "-O3"),
            BuildFlags("O2-native", "-O2 -march=native"),
            BuildFlags("O2-no-fast-math", "-O2 -fno-fast-math")]
SAFETY = 10.0           # line-scale rule
FLOOR = 1e-6            # rtol floor (single precision)
MAX_PRINT_ULPS = 10.0   # deviations up to this many ulps count as print-level noise
ULPS_SAFETY = 3.0       # print-level rule: rounding flips are quantised
ULPS_FLOOR = 2.0
LAST_DIGIT_PATTERNS = {"ROUT": 1e-4, "ROUTSHORT": 1e-4}  # pattern: last printed unit


def pairs(actual_text, expected_text):
    """[(dev, |b|, ulp(b), line max, file max)] for all floats, or (None, problem)."""
    act, exp = tokenize(actual_text), tokenize(expected_text)
    if len(act) != len(exp):
        return None, f"token count {len(act)} vs {len(exp)}"
    line_max = float_scales(exp, "line")
    file_max = float_scales(exp, "file").default
    out = []
    for a, e in zip(act, exp):
        if a.kind != e.kind or (e.kind != "float" and a.value != e.value):
            return None, f"line {e.line}: {a.value!r} vs {e.value!r}"
        if e.kind == "float":
            out.append((abs(a.value - e.value), abs(e.value), e.ulp, line_max[e.line], file_max))
    return out, None


def collect(tl_root, variants):
    """({pattern: {variant: [pairs]}}, problems)."""
    patterns = list(load_tolerances()["numeric"])
    data, problems = {}, []
    for case, run in iter_runs():
        for flags in variants:
            with tempfile.TemporaryDirectory() as tmp:
                try:
                    execute_run(case, run, tl_root, flags, Path(tmp))
                except Exception as err:  # noqa: BLE001
                    problems.append(f"{case.name}/{run.name} [{flags.name}]: "
                                    f"run failed: {str(err).splitlines()[0]}")
                    continue
                for exp in expected_outputs(case, run):
                    act = Path(tmp) / exp.name
                    where = f"{case.name}/{run.name}/{exp.name} [{flags.name}]"
                    if not act.is_file():
                        problems.append(f"{where}: missing")
                        continue
                    p, prob = pairs(act.read_text(errors="replace"),
                                    exp.read_text(errors="replace"))
                    if prob:
                        problems.append(f"{where}: {prob}")
                        continue
                    (data.setdefault(pattern_for(exp.name, patterns), {})
                         .setdefault(flags.name, []).extend(p))
    return data, problems


def _scaled_candidate(all_pairs, scale):
    """Absolute part relative to the line or file maximum, rtol for the rest."""
    idx = 3 if scale == "line" else 4
    atol_frac = SAFETY * max((p[0] / p[idx] for p in all_pairs if p[idx] > 0), default=0.0)
    beyond = max((max(0.0, d - atol_frac * p[idx] - ULPS_FLOOR * u) / b
                  for p in all_pairs for d, b, u in [p[:3]] if b > 0), default=0.0)
    return (Tolerance(max(SAFETY * beyond, FLOOR), atol_frac, ulps=ULPS_FLOOR, scale=scale),
            f"{scale}-scale absolute part, rel. beyond {beyond:.1e}")


def propose_one(pattern, all_pairs):
    """[(Tolerance, note)] candidates that each cover the observed noise."""
    if pattern in LAST_DIGIT_PATTERNS:
        unit = LAST_DIGIT_PATTERNS[pattern]
        return [(Tolerance(0.0, 0.0, atol=1.5 * unit), "one unit in the last digit of R")]
    max_ulps = max((p[0] / p[2] for p in all_pairs if p[2] > 0), default=0.0)
    cands = [_scaled_candidate(all_pairs, "file"), _scaled_candidate(all_pairs, "line")]
    if max_ulps <= MAX_PRINT_ULPS:
        cands.append((Tolerance(0.0, 0.0, ulps=max(ULPS_FLOOR, ULPS_SAFETY * max_ulps)),
                      f"print-level noise, max {max_ulps:.1f} ulps"))
    return cands


def insensitivity(tol):
    """Share of non-zero expected floats (all expected files) allowing a 1 % change
    (check_noise_variants.insensitive_counts, the metric tier 1 is reported with)."""
    n = n_bad = 0
    for case, run in iter_runs():
        for exp in expected_outputs(case, run):
            t = tol.get(pattern_for(exp.name, list(tol)))
            if not t:
                continue
            counts = expected_insensitive_counts(exp, t)
            n += counts[0]
            n_bad += counts[1]
    return n_bad / n if n else 0.0


def kept_patterns():
    """Patterns with a structural rule (scale "tensor", file_cap or exact F
    columns): not recalibrated by this script, kept as they are;
    check them with check_noise_variants.py."""
    return {p for p, t in load_tolerances()["numeric"].items()
            if t.get("scale") == "tensor"
            or {"file_cap", "exact_f_decimals", "exact_columns"} & set(t)}


def _toml_block(text, pattern):
    """The [numeric."pattern"] table of text, verbatim (up to the next table)."""
    start = text.index(f'[numeric."{pattern}"]')
    end = text.find("\n[", start)
    return text[start:end if end >= 0 else len(text)].rstrip("\n")


def propose(data):
    """Per pattern, the candidate with the lowest insensitivity."""
    out = {}
    for pattern, by_variant in data.items():
        if pattern in kept_patterns():
            continue
        all_pairs = [p for ps in by_variant.values() for p in ps]
        best = None
        for tol, note in propose_one(pattern, all_pairs):
            ins = insensitivity({pattern: tol})
            print(f"  candidate {pattern:10s} {note:45s} insensitive {100 * ins:5.1f} %")
            if best is None or ins < best["insensitive"]:
                best = {"tol": tol, "note": note, "insensitive": ins}
        ok = all(p[0] <= best["tol"].allowed(p[1], p[2], p[3] if best["tol"].scale == "line"
                                             else p[4]) for p in all_pairs)
        out[pattern] = {**best, "ok": ok}
    return out


def render_toml(proposed, variants):
    """Text of a new tolerances.toml: proposed and kept patterns, then the old tier2 part."""
    old = TOLERANCES_FILE.read_text()
    rest = old[old.index("[tier2]"):]
    lines = [
        "# Numerical tolerances per output file pattern: a float a matches the expected b if",
        "#   |a-b| <= rtol*|b| + atol_frac*S + atol + ulps*ulp(b)",
        "# where ulp(b) is b's unit in the last printed digit and S is the largest |float| of the",
        "# file (scale = \"file\", default) or of b's line (scale = \"line\");",
        "# scale = \"tensor\": per WRITE block of a T_* file (see that table; kept as is",
        "# by calibrate_tolerances.py, checked by check_noise_variants.py).",
        "# Patterns are matched with fnmatch against the output file name; the first match wins.",
        "# Optional: atol, ulps, scale, skip_lines (1-based), ignore_line_patterns (regexes).",
        "# skip_lines and ignore_line_patterns hide differences: only entries on SKIP_ALLOWLIST in",
        "# tests/unit/test_tolerance_config.py may use them (see the comment there for the rule).",
        "# Calibrated by tests/generate/calibrate_tolerances.py (see README, 'Tolerances'):",
        f"#   date {datetime.date.today().isoformat()}, baseline -O2 vs variants "
        + ", ".join(f"'{v.compile}'" for v in variants) + " (gfortran; ifort not available).",
        "# 'insensitive' = share of non-zero expected values for which a 1 % change passes.",
        "# Loosen only with evidence of pure floating-point reordering (README, 'Rules').",
    ]
    order = ["fd.out", "amp.out", "T_*", "ROUT", "ROUTSHORT", "*.column"]
    for pattern in order + [p for p in proposed if p not in order]:
        if pattern in kept_patterns():
            lines += ["", _toml_block(old, pattern)]
            continue
        p = proposed[pattern]
        t = p["tol"]
        lines += ["", f'[numeric."{pattern}"]',
                  f"# {p['note']}; insensitive {100 * p['insensitive']:.1f} %",
                  f"rtol = {t.rtol:.1e}", f"atol_frac = {t.atol_frac:.1e}"]
        if t.atol:
            lines.append(f"atol = {t.atol:.1e}")
        if t.ulps:
            lines.append(f"ulps = {t.ulps:g}")
        if t.scale != "file":
            lines.append(f'scale = "{t.scale}"')
    return "\n".join(lines) + "\n\n" + rest


def main(argv=None):
    """Print the observed deviations and the proposed tolerances; with --write, store them."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tl-root", type=Path, required=True)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    data, problems = collect(args.tl_root, VARIANTS)
    print(f"{'pattern':12s} {'variant':16s} {'max abs':>10s} {'max ulps':>10s} {'max /line':>10s}")
    for pattern, by_variant in data.items():
        for name, ps in by_variant.items():
            print(f"{pattern:12s} {name:16s} {max(d for d, *_ in ps):10.2e} "
                  f"{max((p[0] / p[2] for p in ps if p[2] > 0), default=0):10.1f} "
                  f"{max((p[0] / p[3] for p in ps if p[3] > 0), default=0):10.2e}")
    proposed = propose(data)
    print("\nproposed:")
    for pattern, p in proposed.items():
        t = p["tol"]
        print(f"  {pattern:10s} rtol {t.rtol:.1e} atol_frac {t.atol_frac:.1e} atol {t.atol:.1e} "
              f"ulps {t.ulps:g} scale {t.scale}  [{p['note']}; insensitive "
              f"{100 * p['insensitive']:.1f} %]{'' if p['ok'] else '   <-- SUSPICIOUS'}")
    for prob in problems:
        print("PROBLEM:", prob)
    if args.write:
        if problems or not all(p["ok"] for p in proposed.values()):
            print("not writing: problems or suspicious deviations (see above)")
            return 1
        TOLERANCES_FILE.write_text(render_toml(proposed, VARIANTS))
        print(f"wrote {TOLERANCES_FILE}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
