"""Generate the expected values of the ref-calc routine tests from the baseline.

Builds generate/refcalc_routines_legacy_driver.f against the baseline's
legacy lib/lib.tleed.f (gfortran -w, -O2 and -O0), runs it on every input
set defined in CASES below, checks that both builds give the same bits,
and writes the inputs (stdin of the driver) and the summaries of the
results to tier1/refcalc_routines/values.json, which
tier1/test_refcalc_routines.py compares with tier1/refcalc_routines/driver.f90
built against the tree under test. Covered: CAAA, CELMG, LXGENT and CPPP
(arrays CAA, CLM, FAC, LX, LXI, LT, LXM, PPP) for every LMAX from 2 to 18,
and SLIND (V, VL, JJS) for the lattices in LATTICES. Arrays are summarized
by their count, a SHA-256 of their bits and their first values
(harness/routine_blocks.py).

Result of the -O2 versus -O0 comparison (gfortran 15.2): every array of
every case is bit-identical in the two builds; the generator stops with an
error otherwise.

Usage (env tl-tests):
    python -m generate.make_refcalc_routine_values --tl-root=<baseline tree> \
        [--allow-non-baseline]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import subprocess
from pathlib import Path

import numpy as np

from generate.routine_values import hex32 as _hex
from generate.routine_values import generate, main, run_driver
from harness.routine_blocks import summarize

TESTS_DIR = Path(__file__).parents[1]
LEGACY_DRIVER = Path(__file__).with_name("refcalc_routines_legacy_driver.f")
VALUES = TESTS_DIR / "tier1" / "refcalc_routines" / "values.json"
LINK = ["-llapack", "-lblas"]
LMAX_RANGE = range(2, 19)  # the range of the DATA tables NCA and NLMS

_F32 = np.float32
_BOHR = _F32(0.529177)
_PI = _F32(3.14159265)
_SQ = (2.88, 0.0), (0.0, 2.88)
_IR = (2.7116, 0.0), (0.0, 2.7116)
_HEX = (2.556, 0.0), (-1.278, 2.2136)

# (id, what, substrate vectors, superlattice vectors, IDEG, NL1, NL2); vectors
# in Angstrom as in FIN (F9.4), the substrate vectors as the (1x1) cell
LATTICES = [
    ("ag100", "p(1x1) square (Ag(100) of tier 1)", _SQ, _SQ, 4, 1, 1),
    ("cu111", "p(1x1) hexagonal, IDEG 3 (Cu(111) of tier 1)", _HEX, _HEX, 3, 1, 1),
    ("ir100_2x1", "p(2x1) rectangular, IDEG 2 (Ir(100)-O of tier 1)",
     _IR, ((5.4232, 0.0), (0.0, 2.7116)), 2, 2, 1),
    ("square_p2x2", "p(2x2) square, IDEG 4",
     _SQ, ((5.76, 0.0), (0.0, 5.76)), 4, 2, 2),
    ("square_c2x2", "c(2x2) square (non-rectangular superlattice), IDEG 4",
     _SQ, ((2.88, 2.88), (2.88, -2.88)), 4, 2, 1),
    ("hex_p1x1_ideg6", "p(1x1) hexagonal, IDEG 6", _HEX, _HEX, 6, 1, 1),
    ("hex_p2x2_ideg3", "p(2x2) hexagonal, IDEG 3",
     _HEX, ((5.112, 0.0), (-2.556, 4.4272)), 3, 2, 2),
    ("hex_p2x2_ideg6", "p(2x2) hexagonal, IDEG 6",
     _HEX, ((5.112, 0.0), (-2.556, 4.4272)), 6, 2, 2),
    ("hex_r3_ideg3", "hexagonal (sqrt3 x sqrt3) R30, IDEG 3",
     _HEX, ((3.834, 2.2136), (-3.834, 2.2136)), 3, 3, 1),
    ("hex_r3_ideg6", "hexagonal (sqrt3 x sqrt3) R30, IDEG 6",
     _HEX, ((3.834, 2.2136), (-3.834, 2.2136)), 6, 3, 1),
]


def _reciprocal(a):
    """Reciprocal vectors of the cell a (rows: first and second basis vector) as
    read_input of ref-calc computes them, in 32-bit arithmetic."""
    area = abs(a[0][0] * a[1][1] - a[1][0] * a[0][1])
    atv = _F32(2) * _PI / area
    return [(a[1][1] * atv, -a[1][0] * atv), (-a[0][1] * atv, a[0][0] * atv)]


def _sublattice_stdin(substrate, superlattice, ideg, nl1, nl2):
    """The stdin lines of an SLIND case: the cells as read_input of ref-calc
    leaves them (Angstrom to bohr), and the reciprocal superlattice."""
    a, b = ([[_F32(x) / _BOHR for x in vec] for vec in cell]
            for cell in (substrate, superlattice))
    rb = _reciprocal(b)
    values = [*a[0], *a[1], *b[0], *b[1], *rb[0], *rb[1]]
    return ["SLIND", f"{ideg} {nl1} {nl2}", *map(_hex, values)]


# (id, what it covers, stdin lines of the driver)
CASES = [
    (f"angular_lmax{lmax:02d}", f"CAAA, CELMG, LXGENT and CPPP at LMAX = {lmax}",
     ["ANGULAR", str(lmax)])
    for lmax in LMAX_RANGE
] + [
    (f"slind_{ident}", f"SLIND, {what}, NL1 x NL2 = {nl1} x {nl2}",
     _sublattice_stdin(substrate, superlattice, ideg, nl1, nl2))
    for ident, what, substrate, superlattice, ideg, nl1, nl2 in LATTICES
]


def _build(tl, workdir, flags):
    """Compile the legacy driver with the baseline lib.tleed.f; return the executable."""
    lib = workdir / f"lib{flags}.o"
    exe = workdir / f"legacy{flags}"
    subprocess.run(["gfortran", flags, "-w", "-I", str(tl / "TensErLEED" / "src"),
                    "-c", "-o", str(lib), str(tl / "TensErLEED" / "lib" / "lib.tleed.f")],
                   check=True, cwd=workdir)
    subprocess.run(["gfortran", flags, "-w", "-o", str(exe), str(LEGACY_DRIVER),
                    str(lib), *LINK], check=True, cwd=workdir)
    return exe


def _run(exe, stdin):
    """Run the driver on the stdin lines; return the summary of its output."""
    return summarize(run_driver(exe, stdin))


def _generate(argv=None):
    """Generate tier1/refcalc_routines/values.json from the baseline."""
    generate(argv, __doc__.splitlines()[0], _build, _run, CASES, VALUES,
             "lib/lib.tleed.f")


if __name__ == "__main__":
    main(_generate)
