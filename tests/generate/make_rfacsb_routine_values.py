"""Generate the expected values of the rfacsb routine tests from the baseline.

Builds generate/rfacsb_routines_legacy_driver.f against the baseline's
legacy lib/rfacsb.f (gfortran -w, -O2 and -O0), runs it on every input set
defined in CASES below, checks that both builds give the same bits, and
writes the inputs (stdin of the driver) and the results to
tier1/rfacsb_routines/values.json, which tier1/test_rfacsb_routines.py
compares with tier1/rfacsb_routines/driver.f90 built against the tree
under test. Reals are 32-bit and stored as the 8 hex digits of their bit
pattern, so the comparison is exact.

Usage (env tl-tests):
    python -m generate.make_rfacsb_routine_values --tl-root=<baseline tree> \
        [--allow-non-baseline]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-22'
__license__ = 'GPLv3+'

import math
import subprocess
from pathlib import Path

from generate.routine_values import hex32 as _hex
from generate.routine_values import generate, main, run_driver

TESTS_DIR = Path(__file__).parents[1]
LEGACY_DRIVER = Path(__file__).with_name("rfacsb_routines_legacy_driver.f")
VALUES = TESTS_DIR / "tier1" / "rfacsb_routines" / "values.json"


def _curve(n, offset=0.0):
    """A smooth, non-polynomial test curve of n points."""
    return [1.0 + 0.8 * math.sin(0.7 * (i + offset)) + 0.05 * (i + offset)**2
            for i in range(1, n + 1)]


def _der(y, h=0.5):
    return ["DER", str(len(y)), _hex(h), *map(_hex, y)]


def _smooth(ismoth, e, a):
    return ["SMOOTH", str(ismoth), str(len(e)), *map(_hex, e), *map(_hex, a)]


def _ysmooth(a, ap, vi=4.5, eincr=0.5):
    return ["YSMOOTH", str(len(a)), _hex(vi), _hex(eincr), *map(_hex, a), *map(_hex, ap)]


# (id, what it covers, stdin lines of the driver)
CASES = [
    ("der_n1", "DER, one point: derivative zero", _der([3.7])),
    ("der_n2", "DER, two points: two-point difference at both", _der([1.3, 2.9])),
    ("der_n3", "DER, three points: central difference and 3-point ends",
     _der([1.0, 2.5, 1.7])),
    ("der_short", "DER, 8 points (3 <= NE <= 22): central differences, 3-point ends",
     _der(_curve(8))),
    ("der_end_order", "DER, 5 points: the last three values make any reordering of the "
     "3-point end formula for Y1(NE) round differently",
     _der([1.2, 0.9, 0.569, 2.384, 1.845])),
    ("der_long", "DER, 25 points (NE >= 23): 7-point central, 4-point ends",
     _der(_curve(25))),
    ("smooth_n1", "SMOOTH, one point: unchanged", _smooth(1, [10.0], [2.3])),
    ("smooth_n2", "SMOOTH, two points: unchanged", _smooth(1, [10.0, 10.5], [2.3, 1.1])),
    ("smooth_n3", "SMOOTH, three points on equal intervals",
     _smooth(1, [10.0, 10.5, 11.0], [1.0, 2.5, 1.7])),
    ("smooth_short", "SMOOTH twice, 8 points on equal intervals",
     _smooth(2, [10.0 + 0.5 * i for i in range(8)], _curve(8))),
    ("smooth_nonuniform", "SMOOTH, 6 points on unequal intervals",
     _smooth(1, [10.0, 10.5, 11.2, 11.7, 12.5, 13.0], _curve(6, offset=0.3))),
    ("ysmooth_n1", "YSMOOTH, one point: no second derivative",
     _ysmooth([2.3], [0.4])),
    ("ysmooth_n2", "YSMOOTH, two points: no second derivative",
     _ysmooth([2.3, 1.1], [0.4, -0.7])),
    ("ysmooth_n3", "YSMOOTH, three points: one-sided second derivative at both ends",
     _ysmooth([1.0, 0.2, 1.5], [-1.6, 0.5, 2.6])),
    ("ysmooth_short", "YSMOOTH, 8 points: zero and tiny intensities, both curvatures",
     _ysmooth([0.0, 1e-12, 0.6, 2.0, 2.4, 1.9, 1.2, 1.3],
              [0.0, 1.2, 2.0, 1.8, -0.1, -1.2, -0.6, 0.2])),
]


def _build(tl, workdir, flags):
    """Compile the legacy driver with the baseline rfacsb.f; return the executable."""
    exe = workdir / f"legacy{flags}"
    subprocess.run(["gfortran", flags, "-w", "-o", str(exe), str(LEGACY_DRIVER),
                    str(tl / "TensErLEED" / "lib" / "rfacsb.f")],
                   check=True, cwd=workdir)
    return exe


def _run(exe, stdin):
    """Run the driver on the stdin lines; return its output lines (lowercase hex)."""
    return [line.strip().lower() for line in run_driver(exe, stdin).splitlines()]


def _generate(argv=None):
    """Generate tier1/rfacsb_routines/values.json from the baseline."""
    generate(argv, __doc__.splitlines()[0], _build, _run, CASES, VALUES,
             "lib/rfacsb.f")


if __name__ == "__main__":
    main(_generate)
