"""Generate the expected values of the ref-calc routine tests from the baseline.

Builds generate/refcalc_routines_legacy_driver.f against the baseline's
legacy lib/lib.tleed.f (gfortran -w, -O2 and -O0), runs it on every input
set defined in CASES below, checks that both builds give the same bits,
and writes the inputs (stdin of the driver) and the summaries of the
results to tier1/refcalc_routines/values.json, which
tier1/test_refcalc_routines.py compares with tier1/refcalc_routines/driver.f90
built against the tree under test. Covered: CAAA, CELMG, LXGENT and CPPP
(arrays CAA, CLM, FAC, LX, LXI, LT, LXM, PPP) for every LMAX from 2 to 18,
and SLIND (V, VL, JJS) for the lattices in LATTICES, one of them with an
IDEG inconsistent with the lattices, so that the no-match result
JJS = NL + 1 is frozen. Of the energy loop: TSCATF (with PSTEMP and
BESSEL) and the concentration-weighted sum of its results over the elements
of a site (a vibrating site whose high-L phase shifts PSTEMP leaves zero
after an early exit, no vibration, no early exit, extrapolation above the
table of phase shifts, a mixed site); BEAMS over energy sweeps for normal and
oblique incidence, with one and with several beam sets, where beams stop
propagating as the energy falls; SPHRM_MOD on a grid of angles for every
LMAX from 2 to 18; BEAMS_YLM for each beam set of a beam list (nonzero
offsets, rows of CYLM beyond the set untouched) and for a window of the
list, at normal and oblique incidence; FMAT for the SLIND lattices and every
IDEG that fits them, at two energies and k-parallel values; and BEAMS, FMAT and
RSMF_SIMPLE (with XM and MFOLD_SIMPLE) chained as in ref-calc, for one and
for several beam sets (RA, TA, AMULT, CYLM, XEVST, XODST); the composite
layer: SRTLAY with GHMAT (GHD, SH, SB, GHSC), TAUMAT (XMT, CXMTXT) and
RTINV_SIMPLE (THMAT, TAUY_SIMPLE, MFOLT_SIMPLE; R/T top and bottom, TSTORE)
on sublayer stacks of Ir(100)-(2x1)-O (three sublayers; coplanar sublayers
in arrays for three; four sublayers with equal interplanar vectors, and with
vectors that differ only in z;
sublayers out of order along the normal; t-matrices scaled up so that the
elimination of CXMTXT interchanges rows), at normal and oblique incidence,
for LMAX < MLMAX, for one call over all beams and for one call per beam set
(offsets NA, NS > 0), as overlayer and as bulk layer; and FMAT with a
superlattice vector 0.15 % longer than the substrate vector (first and second). Phase shifts
and beam lists are those of the tier-1 cases in tests/cases. Arrays are
summarized by their count, a SHA-256 of their bits and their first values
(harness/routine_blocks.py).

The results of the cases whose ids start with a key of LIVE_ARRAYS are
live for the arrays it lists: of the routines driven here only
RSMF_SIMPLE (RA, TA) and RTINV_SIMPLE (R_TOP, T_TOP, R_BOT, T_BOT, TSTORE)
call BLAS/LAPACK (CGETRF, CGETRS), and these results can differ in the last
bit between OpenBLAS and the reference BLAS (all other arrays, of these and
of all other cases, are bit-identical under both; TAUMAT inverts with its
own elimination, CXMTXT). values.json holds the other arrays of these cases
and lists the live ones under "live"; the tier-1 test builds this driver
from the baseline tree (build_legacy_driver) and compares them with its
output, run with the same BLAS/LAPACK as the driver under test.

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

import json
import subprocess
from pathlib import Path

import numpy as np

from generate.baseline_tree import GeneratorError
from generate.routine_values import hex32 as _hex
from generate.routine_values import generate, main, run_driver
from harness.routine_blocks import summarize

TESTS_DIR = Path(__file__).parents[1]
LEGACY_DRIVER = Path(__file__).with_name("refcalc_routines_legacy_driver.f")
VALUES = TESTS_DIR / "tier1" / "refcalc_routines" / "values.json"
LINK = ["-llapack", "-lblas"]
# cases with legacy results that depend on the BLAS/LAPACK the driver runs
# with, and those results
LIVE_ARRAYS = {
    "rsmf_": ("RA", "TA"),
    "rtinv_": ("R_TOP", "T_TOP", "R_BOT", "T_BOT", "TSTORE"),
    }
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
    ("square_p2x2_ideg3", "p(2x2) square with IDEG 3, inconsistent with the "
     "lattices: the 120-degree rotations of the sublattices other than the "
     "first match none, so JJS = NL + 1 for them",
     _SQ, ((5.76, 0.0), (0.0, 5.76)), 3, 2, 2),
]


def _reciprocal(a):
    """Reciprocal vectors of the cell a (rows: first and second basis vector) as
    read_input of ref-calc computes them, in 32-bit arithmetic."""
    area = abs(a[0][0] * a[1][1] - a[1][0] * a[0][1])
    atv = _F32(2) * _PI / area
    return [(a[1][1] * atv, -a[1][0] * atv), (-a[0][1] * atv, a[0][0] * atv)]


def _bohr_cell(cell):
    """The cell (rows: basis vectors, in Angstrom) in bohr, as ref-calc has it."""
    return [[_F32(x) / _BOHR for x in vec] for vec in cell]


def _lattice_lines(substrate, superlattice, ideg, nl1, nl2):
    """The stdin lines of a lattice: the cells as read_input of ref-calc
    leaves them (Angstrom to bohr), and the reciprocal superlattice."""
    a, b = _bohr_cell(substrate), _bohr_cell(superlattice)
    rb = _reciprocal(b)
    values = [*a[0], *a[1], *b[0], *b[1], *rb[0], *rb[1]]
    return [f"{ideg} {nl1} {nl2}", *map(_hex, values)]


CASES_DIR = TESTS_DIR / "cases" / "refcalc"
_HARTREE = _F32(27.211396)
# the energy offset inside the crystal (eV) and the imaginary part of the
# energy (V0i of the tier-1 inputs), and the beam-selection criterion that
# ref-calc derives from TST = 0.001 and an interlayer distance of 2.0365 A
_INNER_EV = 10.0
_VPI = -_F32(4.5) / _HARTREE
_TST = _F32((np.log(_F32(0.001)) / (_F32(2.0365) / _BOHR))**2)
# up and down again: NPUC is not reset between the energies, so a beam to be
# output that propagates on the way up and not on the way down must read zero
_SWEEP_EV = (10.0, 20.0, 35.0, 60.0, 100.0, 150.0, 250.0,
             150.0, 100.0, 60.0, 35.0, 20.0, 10.0)


def _read_fin(case, run):
    """The pieces of the input of a tier-1 run that the routine cases use:
    cells, beam sets (index pairs), theta and phi (degrees), and the phase
    shift table of the first LMAX + 1 values of the elements in a format-1
    file: (energies in Hartree, values[energy][element][l])."""
    lines = (CASES_DIR / case / "inputs" / f"FIN_{run}").read_text().splitlines()
    cells = [[float(x) for x in lines[i].split()[:2]] for i in range(2, 6)]
    i = 7
    sets = []
    while "TST" not in lines[i]:
        n = int(lines[i])
        sets.append([tuple(float(x) for x in ln.split()[:2])
                     for ln in lines[i + 1:i + 1 + n]])
        i += n + 1
    theta, phi = (float(x) for x in lines[i + 2].split()[:2])
    lmax = int(lines[i + 5].split()[0])
    j = next(k for k, ln in enumerate(lines) if "PSFORMAT" in ln) + 2
    nel = int(lines[j].split()[0])
    j += 1
    energies, values = [], []
    while not lines[j].startswith("---"):
        energies.append(_F32(lines[j][:7]))
        values.append([[_F32(x) for x in
                        (lines[j + 1 + 2 * k][i7:i7 + 7] for i7 in range(0, 70, 7))][:lmax + 1]
                       for k in range(nel)])
        j += 1 + 2 * nel
    return {"cells": cells, "sets": sets, "theta": theta, "phi": phi,
            "lmax": lmax, "es": energies, "phases": values}


def _crystal_energy(ev):
    """Energy in the crystal (Hartree) for the vacuum energy ev (eV)."""
    return _F32(ev) / _HARTREE + _F32(_INNER_EV) / _HARTREE


def _k_parallel(ev, theta, phi):
    """Parallel components of the incident wave vector (bohr^-1)."""
    ksq = _F32(2) * (_F32(ev) / _HARTREE)
    k = np.sqrt(ksq) * np.sin(_F32(np.radians(theta)))
    return k * np.cos(_F32(np.radians(phi))), k * np.sin(_F32(np.radians(phi)))


def _beam_lines(fin):
    """knbs, set sizes, then the momenta (substrate reciprocal vectors) and
    index pairs of all beams as the driver's BEAMS and RSMF inputs have them."""
    rec = _reciprocal(_bohr_cell(fin["cells"][:2]))
    beams = [b for s in fin["sets"] for b in s]
    spq = [_F32(i1) * rec[0][c] + _F32(i2) * rec[1][c] for i1, i2 in beams for c in (0, 1)]
    spqf = [_F32(x) for b in beams for x in b]
    return [str(len(fin["sets"])), " ".join(str(len(s)) for s in fin["sets"])], spq, spqf


def _beams_stdin(fin, theta, phi, energies):
    head, spq, spqf = _beam_lines(fin)
    lines = ["BEAMS", head[0], head[1], "4", "2 6 10 14", _hex(_TST),
             *map(_hex, spq), *map(_hex, spqf), str(len(energies))]
    for ev in energies:
        ak2, ak3 = _k_parallel(ev, theta, phi)
        lines += map(_hex, [_crystal_energy(ev), ak2, ak3])
    return lines


def _tscatf_stdin(table, lmax, nel, energy, site):
    """TSCATF input: the phase shift table cut to lmax and nel elements, the
    energy (Hartree), and the site as [(iel, concentration, vibration)] or
    [(iel, concentration, vibration perpendicular, vibration parallel)]."""
    es = table["es"]
    phss = [[[table["phases"][ie][iel][l] for l in range(lmax + 1)]
             for iel in range(nel)] for ie in range(len(es))]
    flat = [phss[ie][iel][l] for l in range(lmax + 1) for iel in range(nel)
            for ie in range(len(es))]
    lines = ["TSCATF", f"{lmax} {nel} {len(es)}", *map(_hex, es), *map(_hex, flat),
             _hex(energy), str(len(site))]
    for iel, conc, vib, *parallel in site:
        lines += [str(iel), _hex(conc), _hex(vib), _hex(*parallel or [vib])]
    return lines


def _sphrm_stdin():
    """SPHRM input: real and complex direction cosines on a grid of angles."""
    points = []
    for theta in (0.0, 15.0, 35.0, 60.0, 80.0, 89.5):
        for phi in ((0.0,) if theta == 0.0 else (0.0, 40.0, 135.0, 250.0)):
            t, f = _F32(np.radians(theta)), _F32(np.radians(phi))
            points.append(((np.cos(t), 0.0), (np.sin(t) if theta else _F32(0.0), 0.0),
                           (np.cos(f), np.sin(f))))
    for sinh, phi in ((0.7, 30.0), (1.5, 100.0), (0.1, 200.0)):  # evanescent beams
        f = _F32(np.radians(phi))
        ct = (0.0, _F32(sinh))
        points.append((ct, (_F32(np.sqrt(1 + sinh**2)), 0.0), (np.cos(f), np.sin(f))))
    values = [v for k in range(3) for p in points for v in p[k]]
    return ["SPHRM", str(len(points)), *map(_hex, values)]


def _bylm_stdin(fin, lmax, ev, angles):
    """BYLM input: the beam list of fin at the energy ev (eV)."""
    ak2, ak3 = _k_parallel(ev, *angles)
    head, spq, spqf = _beam_lines(fin)
    return ["BYLM", *map(_hex, [_crystal_energy(ev), _VPI, ak2, ak3]), str(lmax),
            *head, _hex(_TST), *map(_hex, spq), *map(_hex, spqf)]


def _fmat_stdin(substrate, superlattice, ideg, nl1, nl2, lmax,
                incidences=((50.0, (0.0, 0.0)), (120.0, (15.0, 25.0)))):
    """FMAT input: the incidences (energy in eV, angles) in turn; by default
    normal incidence at 50 eV and oblique incidence at 120 eV."""
    lines = ["FMAT", *_lattice_lines(substrate, superlattice, ideg, nl1, nl2),
             str(lmax), str(len(incidences))]
    for ev, angles in incidences:
        ak2, ak3 = _k_parallel(ev, *angles)
        lines += map(_hex, [_crystal_energy(ev), _VPI, ak2, ak3])
    return lines


def _rsmf_stdin(fin, ideg, nl1, nl2, lmax, ev, angles, mode):
    """RSMF input: the lattices and beam list of fin at the energy ev (eV);
    the atomic t-matrix of its first phase shifts, as sin(d) exp(i d)."""
    cells = _bohr_cell(fin["cells"])
    area = lambda cell: abs(cell[0][0] * cell[1][1] - cell[1][0] * cell[0][1])
    tv = area(cells[2:] if mode == "O" else cells[:2])
    ak2, ak3 = _k_parallel(ev, *angles)
    delta = [fin["phases"][0][0][l] if l < len(fin["phases"][0][0]) else _F32(0)
             for l in range(lmax + 1)]
    tmat = [x for d in delta for x in
            ((np.sin(d) * np.cos(d)), (np.sin(d) * np.sin(d)))]  # sin(d) exp(i d)
    head, spq, spqf = _beam_lines(fin)
    return ["RSMF", *_lattice_lines(fin["cells"][:2], fin["cells"][2:], ideg, nl1, nl2),
            str(lmax),
            *map(_hex, [_crystal_energy(ev), _VPI, ak2, ak3, tv]), mode,
            *map(_hex, tmat), *head, _hex(_TST), *map(_hex, spq), *map(_hex, spqf)]


def _read_layers(case, run):
    """The layer types of a tier-1 input: [{"lay": LAY, "subs": [(site
    type, x, y, z)]}] with the coordinates in Angstrom (x along the surface
    normal)."""
    lines = (CASES_DIR / case / "inputs" / f"FIN_{run}").read_text().splitlines()
    i = next(k for k, ln in enumerate(lines) if "NLTYPE" in ln)
    layers = []
    for _ in range(int(lines[i].split()[0])):
        i += 1
        while lines[i].startswith("-"):
            i += 1
        nsub = int(lines[i + 1].split()[0])
        subs = [(int(w[0]), *(float(x) for x in w[1:4]))
                for w in (ln.split() for ln in lines[i + 2:i + 2 + nsub])]
        layers.append({"lay": int(lines[i].split()[0]), "subs": subs})
        i += 1 + nsub
    return layers


def _sublayer(layer, isub):
    """Sublayer isub (0-based) of a layer type, zeros beyond its last one."""
    return layer["subs"][isub] if isub < len(layer["subs"]) else (0, 0.0, 0.0, 0.0)


def _composite_stdin(name, fin, layers, setup):
    """The stdin of a composite-layer case (GHMAT, TAUMAT or RTINV): the
    layer type setup["layer"] (1-based) of layers on the lattice of fin, as
    ref-calc sets it up: the t-matrices of its site types contracted to
    those in use (TSF, LPS), the sublayer positions of all layer types
    (SUBPOS(layer type, sublayer, 3), in bohr), the incidence, and the beam
    list of fin. In mode O one call covers all beams, in mode B one call
    per beam set."""
    lmax, mlmax, mnsub = setup["lmax"], setup["mlmax"], setup["mnsub"]
    layer = layers[setup["layer"] - 1]
    lay = layer["lay"]
    cells = _bohr_cell(fin["cells"])
    area = lambda cell: abs(cell[0][0] * cell[1][1] - cell[1][0] * cell[0][1])
    tv = area(cells[2:] if lay == 1 else cells[:2])
    ak2, ak3 = _k_parallel(setup["ev"], *setup["angles"])
    subpos = [_F32(_sublayer(lay_, isub)[1 + k]) / _BOHR for k in range(3)
              for isub in range(mnsub) for lay_ in layers]
    sites = {}
    lps = []
    tsf = [[(_F32(0), _F32(0))] * (mlmax + 1) for _ in range(mnsub)]
    for site, *_ in layer["subs"]:
        if site not in sites:
            sites[site] = len(sites) + 1
            for l in range(lmax + 1):
                d = fin["phases"][0][site - 1][l]
                tsf[sites[site] - 1][l] = (np.sin(d) * np.cos(d), np.sin(d) * np.sin(d))
                tsf[sites[site] - 1][l] = tuple(
                    _F32(setup.get("tsf_scale", 1)) * x for x in tsf[sites[site] - 1][l])
        lps.append(sites[site])
    tsf_flat = [x for l in range(mlmax + 1) for row in tsf for x in row[l]]
    head, spq, spqf = _beam_lines(fin)
    return [name, *_lattice_lines(fin["cells"][:2], fin["cells"][2:], setup["ideg"],
                                  setup["nl1"], setup["nl2"]),
            str(lmax), f"{mlmax} {mnsub} {len(layer['subs'])} {lay} {len(layers)} "
            f"{setup['layer']}", setup["mode"],
            *map(_hex, [_crystal_energy(setup["ev"]), _VPI, ak2, ak3, tv, _F32(0.001)]),
            *map(_hex, subpos), str(len(sites)), " ".join(map(str, lps)),
            *map(_hex, tsf_flat), *head, _hex(_TST), *map(_hex, spq), *map(_hex, spqf)]


_AG = _read_fin("ag100_mix", "E141_L08")
_AG_OBL = _read_fin("ag100_oblique", "E141_L08")
_IR = _read_fin("ir100_2x1_o", "E118_L07")
_CU = _read_fin("cu111", "E140_L07")
_TABLE = {"es": _AG["es"][::8], "phases": _AG["phases"][::8]}
_E_TABLE = _TABLE["es"][3] * _F32(1.1)  # between two tabulated energies
_E_EARLY = _TABLE["es"][0] + (_TABLE["es"][1] - _TABLE["es"][0]) * _F32(0.3)
_E_ABOVE = _TABLE["es"][-1] * _F32(1.3)  # above the last one: extrapolation
_VIB = _F32(0.18)  # about 0.1 Angstrom in bohr
if not _TABLE["es"][0] < _E_TABLE < _TABLE["es"][-1] < _E_ABOVE:
    raise GeneratorError("the energies of the TSCATF cases do not fit the table")

# (id, what, number of L values, element count, energy, site)
_TSCATF_CASES = [
    ("vibrating", "a vibrating Ag site; PSTEMP exits early, the DEL of L > 6 stay zero",
     8, 2, _E_EARLY, [(1, _F32(1), _VIB)]),
    ("no_vibration", "VIB = 0 (the small-FALFE branch of PSTEMP)",
     8, 2, _E_TABLE, [(1, _F32(1), _F32(0))]),
    ("no_early_exit", "a vibrating Pd site at LMAX = 2: PSTEMP never exits early",
     2, 2, _E_TABLE, [(2, _F32(1), _VIB)]),
    ("extrapolation", "an energy above the table of phase shifts",
     8, 2, _E_ABOVE, [(1, _F32(1), _VIB)]),
    ("anisotropic", "different vibration amplitudes perpendicular and parallel "
     "to the surface (averaged to DR)",
     8, 2, _E_TABLE, [(1, _F32(1), _VIB, _F32(1.7) * _VIB)]),
    ("ata_mix", "a mixed site: Ag 0.7, Pd 0.3 with other vibration amplitudes, "
     "and an element of concentration 0 that is skipped",
     8, 2, _E_TABLE, [(1, _F32(0.7), _VIB), (2, _F32(0.3), _F32(1.2) * _VIB),
                      (1, _F32(0), _VIB)]),
]
# (id, what, input data, angles of incidence)
_BEAM_CASES = [
    ("ag100_normal", "one beam set, normal incidence", _AG, (0.0, 0.0)),
    ("ag100_oblique", "one beam set, oblique incidence", _AG_OBL,
     (_AG_OBL["theta"], _AG_OBL["phi"])),
    ("ir100_sets_normal", "p(2x1) superstructure, two beam sets, normal incidence",
     _IR, (0.0, 0.0)),
    ("ir100_sets_oblique", "p(2x1) superstructure, two beam sets, oblique incidence",
     _IR, (20.0, 10.0)),
]
# (id, what, input data, LMAX, energy (eV), angles)
_BYLM_CASES = [
    ("ir100_sets_normal", "two beam sets (offset > 0), normal incidence",
     _IR, 7, 100.0, (0.0, 0.0)),
    ("ir100_sets_oblique", "two beam sets (offset > 0), oblique incidence",
     _IR, 7, 100.0, (20.0, 10.0)),
    ("ir100_sets_low_energy", "two beam sets, oblique incidence, at an energy "
     "where beams are evanescent", _IR, 4, 20.0, (30.0, 70.0)),
    ("cu111_oblique", "one beam set, hexagonal lattice, oblique incidence",
     _CU, 5, 80.0, (15.0, 20.0)),
    ("ag100_theta0p05", "one beam set, almost normal incidence",
     _AG_OBL, 6, 45.0, (0.05, 30.0)),
]
# FMAT: the lattices of SLIND with their IDEG, 180-degree axes (IDEG 2) of
# lattices of a higher symmetry, and 20 x 1 and 1 x 20 superlattices (IDEG 2)
# whose sublattice vectors fall at 0.05, ..., 0.95 of a superlattice vector
# (the folding of the sublattice vectors into the superlattice cell moves
# them only for fractions beyond 0.999 or below -0.001 and 0.001), also
# along the opposite direction of the vectors
_SQ_20 = 20 * 2.88
_FMAT_LATTICES = [row for row in LATTICES if row[0] != "square_p2x2_ideg3"] + [
    (f"{ident}_ideg2", f"{what}, as IDEG 2", sub, sup, 2, nl1, nl2)
    for ident, what, sub, sup, _, nl1, nl2 in LATTICES
    if ident in ("ag100", "cu111", "square_p2x2", "hex_p2x2_ideg3")
] + [
    ("square_20x1", "p(20x1) square, sublattices along the first vector", _SQ,
     ((_SQ_20, 0.0), (0.0, 2.88)), 2, 20, 1),
    ("square_20x1_reversed", "p(20x1) square, vectors against the substrate "
     "vectors", _SQ, ((-_SQ_20, 0.0), (0.0, -2.88)), 2, 20, 1),
    ("square_1x20", "p(1x20) square, sublattices along the second vector", _SQ,
     ((2.88, 0.0), (0.0, _SQ_20)), 2, 1, 20),
]
# FMAT at small parallel wave vector (the thresholds on its length in the
# lattice sums and in the beam matrices), and with a superlattice vector 0.08 %
# longer than the substrate vector: the second sublattice vector, equal to the
# substrate vector, is at 0.9992 of it, beyond the 0.999 of the folding into
# the superlattice cell (IDEG 1: the rotation matches the vector with itself)
_SMALL_THETA = (0.05, 1.0)
_FMAT_SMALL_K = [
    (f"ag100_theta{str(theta).replace('.', 'p')}",
     f"p(1x1) square at {theta:g} degrees", _SQ, _SQ, 4, 1, 1,
     ((45.0, (theta, 30.0)), (120.0, (theta, 200.0))))
    for theta in _SMALL_THETA
] + [
    ("square_fold", "square lattice with a first superlattice vector 0.08 % "
     "longer than the substrate vector", _SQ, ((2.8823, 0.0), (0.0, 2.88)), 1, 2, 1,
     ((50.0, (0.0, 0.0)), (120.0, (15.0, 25.0)))),
]
# (id, what, input data, IDEG, NL1, NL2, LMAX, energy (eV), angles, mode)
_RSMF_CASES = [
    *[(f"ag100_theta{str(theta).replace('.', 'p')}_overlayer",
       f"one beam set, {theta:g} degrees, overlayer-type layer",
       _AG_OBL, 4, 1, 1, 6, 45.0, (theta, 30.0), "O") for theta in _SMALL_THETA],
    ("ag100_normal_overlayer", "one beam set, normal incidence, overlayer-type layer",
     _AG, 4, 1, 1, 6, 100.0, (0.0, 0.0), "O"),
    ("ag100_oblique_overlayer", "one beam set, oblique incidence, overlayer-type layer",
     _AG_OBL, 4, 1, 1, 6, 100.0, (12.0, 33.0), "O"),
    ("ir100_sets_normal_substrate", "two beam sets (offset > 0), normal incidence, "
     "substrate-type layer of two sublattices", _IR, 2, 2, 1, 7, 100.0, (0.0, 0.0), "B"),
    ("ir100_sets_oblique_substrate", "two beam sets (offset > 0), oblique incidence, "
     "substrate-type layer of two sublattices", _IR, 2, 2, 1, 7, 100.0, (20.0, 10.0), "B"),
    ("ir100_sets_overlayer", "two beam sets in one call, overlayer-type layer",
     _IR, 2, 2, 1, 7, 100.0, (0.0, 0.0), "O"),
    ("cu111_oblique_substrate", "hexagonal lattice, IDEG 3, oblique incidence, "
     "substrate-type layer", _CU, 3, 1, 1, 5, 80.0, (15.0, 20.0), "B"),
    # 2D superstructure: four sublattices (one at 45 degrees) and beam sets
    # whose first beams are fractional, so that the sublattice phase of
    # every set has both components (bk2, bk3) and sums four terms
    ("square_p2x2_sets_substrate", "p(2x2) square, IDEG 4, four beam sets (first "
     "beams (1/2, 1/2) etc.), substrate-type layer of four sublattices",
     {**_AG_OBL, "cells": [*_AG_OBL["cells"][:2], [5.76, 0.0], [0.0, 5.76]],
      "sets": [[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)],
               [(0.5, 0.0), (-0.5, 0.0)], [(0.0, 0.5), (0.0, -0.5)],
               [(0.5, 0.5), (-0.5, 0.5), (0.5, -0.5), (-0.5, -0.5)]]},
     4, 2, 2, 6, 100.0, (12.0, 33.0), "B"),
]


_IR_LAYERS = _read_layers("ir100_2x1_o", "E118_L07")
_UNSORTED_LAYERS = _read_layers("ir100_2x1_o_unsorted", "E118_L07")
_THICK_LAYERS = _read_layers("ir100_2x1_o_thick", "E118_L07")
# a composite bulk layer: two sublayers of the substrate 2.7116 Angstrom
# square cell, 1.9174 Angstrom apart along the normal, shifted by half a cell
_BULK_LAYERS = [{"lay": 2, "subs": [(3, 0.0, 0.0, 0.0), (3, 1.9174, 1.3558, 1.3558)]}]
# the thick layer with its last sublayer moved to (1.9174, 0, 2.7116): the
# interplanar vectors to it and to the one before differ only in z
_Z_ONLY_LAYERS = [
    *_THICK_LAYERS[:1],
    {**_THICK_LAYERS[1], "subs": [*_THICK_LAYERS[1]["subs"][:3],
                                  (3, 1.9174, 0.0, 2.7116)]},
    *_THICK_LAYERS[2:],
]
# a composite bulk layer of two coplanar sublayers (direct lattice sum only)
_BULK_COPLANAR_LAYERS = [{"lay": 2, "subs": [(3, 0.0, 0.0, 0.0), (3, 0.0, 1.3558, 1.3558)]}]
_IR_LATTICE = {"ideg": 2, "nl1": 2, "nl2": 1}
# (id, what, input data, layer types, setup): setup holds the layer type, LMAX
# and MLMAX (the PARAM size: it sets L2M = 2 * MLMAX + 1), MNSUB (the size of
# the sublayer arrays), the energy and angles of incidence, the mode (one call
# for all beams, O, or one per beam set, B) and the lattice
_COMPOSITE_SETUPS = [
    ("ir100_top_normal", "the three sublayers of the top layer of Ir(100)-(2x1)-O "
     "(O and two Ir at different heights), overlayer-type, normal incidence",
     _IR, _IR_LAYERS, dict(layer=1, lmax=7, mlmax=7, mnsub=3, ev=100.0, angles=(0.0, 0.0),
                           mode="O", **_IR_LATTICE)),
    ("ir100_top_oblique", "the same layer, oblique incidence",
     _IR, _IR_LAYERS, dict(layer=1, lmax=7, mlmax=7, mnsub=3, ev=118.0, angles=(20.0, 10.0),
                           mode="O", **_IR_LATTICE)),
    ("ir100_unsorted_oblique", "the sublayers of that layer listed out of order along "
     "the normal, oblique incidence, LMAX 6",
     _IR, _UNSORTED_LAYERS, dict(layer=1, lmax=6, mlmax=6, mnsub=3, ev=80.0,
                                 angles=(15.0, 40.0), mode="O", **_IR_LATTICE)),
    ("ir100_coplanar", "two coplanar Ir sublayers (direct lattice sum only) in arrays "
     "for three (NLAY < MNSUB)",
     _IR, _IR_LAYERS, dict(layer=3, lmax=6, mlmax=6, mnsub=3, ev=100.0, angles=(0.0, 0.0),
                           mode="O", **_IR_LATTICE)),
    ("ir100_thick_equal_vectors", "four Ir sublayers with equal interplanar vectors "
     "(NGEQ), oblique incidence",
     _IR, _THICK_LAYERS, dict(layer=2, lmax=5, mlmax=5, mnsub=4, ev=100.0,
                              angles=(10.0, 20.0), mode="O", **_IR_LATTICE)),
    ("ir100_bulk_coplanar", "a bulk-type layer of two coplanar sublayers (direct "
     "lattice sum only), oblique incidence, in arrays for three sublayers",
     _IR, _BULK_COPLANAR_LAYERS, dict(layer=1, lmax=6, mlmax=6, mnsub=3, ev=100.0,
                                      angles=(20.0, 10.0), mode="O", **_IR_LATTICE)),
    ("ir100_thick_z_only_vectors", "the same four sublayers with the last one moved so "
     "that some interplanar vectors differ only in the third component, oblique "
     "incidence",
     _IR, _Z_ONLY_LAYERS, dict(layer=2, lmax=5, mlmax=5, mnsub=4, ev=100.0,
                               angles=(10.0, 20.0), mode="O", **_IR_LATTICE)),
    ("ir100_strong_scatterers", "the top layer with t-matrices 20 times as large "
     "(not physical: the pivot search of the elimination in CXMTXT interchanges rows), "
     "oblique incidence",
     _IR, _IR_LAYERS, dict(layer=1, lmax=5, mlmax=5, mnsub=3, ev=100.0, angles=(20.0, 10.0),
                           mode="O", tsf_scale=20, **_IR_LATTICE)),
    ("ir100_lmax5_mlmax7", "the top layer at LMAX 5 with arrays for MLMAX 7 "
     "(LMAX < MLMAX), oblique incidence",
     _IR, _IR_LAYERS, dict(layer=1, lmax=5, mlmax=7, mnsub=3, ev=90.0, angles=(12.0, 33.0),
                           mode="O", **_IR_LATTICE)),
    ("ir100_overlayer_sets", "the top layer, one call per beam set (offset > 0), "
     "normal incidence",
     _IR, _IR_LAYERS, dict(layer=1, lmax=5, mlmax=5, mnsub=3, ev=100.0, angles=(0.0, 0.0),
                           mode="B", **_IR_LATTICE)),
    ("ir100_bulk_sets", "a bulk-type layer of two sublayers, one call per beam set "
     "(offset > 0), oblique incidence, in arrays for three sublayers",
     _IR, _BULK_LAYERS, dict(layer=1, lmax=6, mlmax=6, mnsub=3, ev=100.0,
                             angles=(20.0, 10.0), mode="B", **_IR_LATTICE)),
]
_COMPOSITE_ROUTINES = (
    ("GHMAT", "SRTLAY and GHMAT (GHD, SH, SB and GHSC)"),
    ("TAUMAT", "TAUMAT (XMT and CXMTXT)"),
    ("RTINV", "RTINV_SIMPLE (TAUMAT, SRTLAY, GHMAT, THMAT, TAUY_SIMPLE, MFOLT_SIMPLE)"),
)
# FMAT with a superlattice vector 0.15 % longer than the substrate vector: the
# second sublattice vector, equal to the substrate vector, is at 0.9985 of it,
# in (0.998, 0.999], for the first and for the second superlattice vector
_FMAT_FOLD_015 = [
    ("square_fold_015", "square lattice with a first superlattice vector 0.15 % longer "
     "than the substrate vector", _SQ, ((2.8843, 0.0), (0.0, 2.88)), 1, 2, 1),
    ("square_fold_015_second", "square lattice with a second superlattice vector 0.15 % "
     "longer than the substrate vector", _SQ, ((2.88, 0.0), (0.0, 2.8843)), 1, 1, 2),
]
_FOLD_015_INCIDENCES = ((50.0, (0.0, 0.0)), (120.0, (15.0, 25.0)))

# (id, what it covers, stdin lines of the driver)
CASES = [
    (f"angular_lmax{lmax:02d}", f"CAAA, CELMG, LXGENT and CPPP at LMAX = {lmax}",
     ["ANGULAR", str(lmax)])
    for lmax in LMAX_RANGE
] + [
    (f"slind_{ident}", f"SLIND, {what}, NL1 x NL2 = {nl1} x {nl2}",
     ["SLIND", *_lattice_lines(substrate, superlattice, ideg, nl1, nl2)])
    for ident, what, substrate, superlattice, ideg, nl1, nl2 in LATTICES
] + [
    (f"tscatf_{ident}", f"TSCATF, {what}",
     _tscatf_stdin(_TABLE, lmax, nel, energy, site))
    for ident, what, lmax, nel, energy, site in _TSCATF_CASES
] + [
    (f"beams_{ident}", f"BEAMS, {what}, at {len(_SWEEP_EV)} energies, "
     f"{_SWEEP_EV[0]:g} eV up to {max(_SWEEP_EV):g} eV and down again, NPUC not reset",
     _beams_stdin(fin, *angles, _SWEEP_EV))
    for ident, what, fin, angles in _BEAM_CASES
] + [
    ("sphrm_grid", "SPHRM_MOD on a grid of real and complex angles, LMAX = 2 to 18",
     _sphrm_stdin()),
] + [
    (f"bylm_{ident}", f"BEAMS_YLM, {what}, LMAX = {lmax}, {energy:g} eV",
     _bylm_stdin(fin, lmax, energy, angles))
    for ident, what, fin, lmax, energy, angles in _BYLM_CASES
] + [
    (f"fmat_{ident}", f"FMAT, {what}, NL1 x NL2 = {nl1} x {nl2}, IDEG {ideg}, "
     "LMAX = 4, normal incidence at 50 eV and oblique at 120 eV",
     _fmat_stdin(substrate, superlattice, ideg, nl1, nl2, 4))
    for ident, what, substrate, superlattice, ideg, nl1, nl2 in _FMAT_LATTICES
] + [
    (f"fmat_{ident}", f"FMAT, {what}, NL1 x NL2 = {nl1} x {nl2}, IDEG {ideg}, "
     "LMAX = 4, " + " and ".join(f"{ev:g} eV, {angles[0]:g} degrees"
                                 for ev, angles in incidences),
     _fmat_stdin(substrate, superlattice, ideg, nl1, nl2, 4, incidences))
    for ident, what, substrate, superlattice, ideg, nl1, nl2, incidences
    in _FMAT_SMALL_K
] + [
    (f"rsmf_{ident}", f"BEAMS, FMAT and RSMF_SIMPLE, {what}, LMAX = {lmax}, "
     f"{energy:g} eV", _rsmf_stdin(fin, ideg, nl1, nl2, lmax, energy, angles, mode))
    for ident, what, fin, ideg, nl1, nl2, lmax, energy, angles, mode in _RSMF_CASES
] + [
    (f"{name.lower()}_{ident}", f"{routines}, {what}, LMAX = {setup['lmax']}, "
     f"MLMAX = {setup['mlmax']}, mode {setup['mode']}, {setup['ev']:g} eV, "
     f"{setup['angles'][0]:g} degrees", _composite_stdin(name, fin, layers, setup))
    for name, routines in _COMPOSITE_ROUTINES
    for ident, what, fin, layers, setup in _COMPOSITE_SETUPS
] + [
    (f"fmat_{ident}", f"FMAT, {what}, NL1 x NL2 = {nl1} x {nl2}, IDEG {ideg}, LMAX = 4, "
     "50 eV, 0 degrees and 120 eV, 15 degrees",
     _fmat_stdin(substrate, superlattice, ideg, nl1, nl2, 4, _FOLD_015_INCIDENCES))
    for ident, what, substrate, superlattice, ideg, nl1, nl2 in _FMAT_FOLD_015
]

# (case id, array, index of the first value in the leading ones, count): values
# that must be zero (PSTEMP leaves the DEL of the high L at zero) or nonzero
_ZERO_DEL = ("tscatf_vibrating", "DEL_01", 14, 2)
_NONZERO_DEL = ("tscatf_no_early_exit", "DEL_01", 0, 6)


def build_legacy_driver(tl, workdir, flags):
    """Compile the legacy driver with the lib.tleed.f of the tree tl into
    workdir; return the executable."""
    lib = workdir / f"lib{flags}.o"
    exe = workdir / f"legacy{flags}"
    subprocess.run(["gfortran", flags, "-w", "-I", str(tl / "TensErLEED" / "src"),
                    "-c", "-o", str(lib), str(tl / "TensErLEED" / "lib" / "lib.tleed.f")],
                   check=True, cwd=workdir)
    subprocess.run(["gfortran", flags, "-w", "-o", str(exe), str(LEGACY_DRIVER),
                    str(lib), *LINK], check=True, cwd=workdir)
    return exe


def run_legacy_driver(exe, stdin):
    """Run the driver on the stdin lines; return the summary of its output."""
    return summarize(run_driver(exe, stdin))


def _check_composite_cases():
    """Raise GeneratorError unless the GHMAT cases show what they are for:
    equal interplanar vectors (NGEQ), interplanar vectors that differ only in
    the third component and are not equal (all NGEQ zero), sublayers sorted by
    SRTLAY (NORD), and coplanar sublayers (the x component of the interplanar
    vector is zero)."""
    cases = {case["id"]: case["expected"] for case in json.loads(VALUES.read_text())["cases"]}
    heads = lambda case_id, array: cases[f"ghmat_ir100_{case_id}"][array]["head"]
    checks = [
        ("thick_equal_vectors", any(n != "0" for n in heads("thick_equal_vectors", "NGEQ"))),
        ("thick_z_only_vectors", all(n == "0" for n in heads("thick_z_only_vectors", "NGEQ"))),
        ("unsorted_oblique", heads("unsorted_oblique", "NORD") != ["1", "2", "3"]),
        ("coplanar", heads("coplanar", "DRL")[0] == "00000000"),
        ]
    failed = [case_id for case_id, ok in checks if not ok]
    if failed:
        raise GeneratorError(f"the cases {failed} do not show what they are for")


def _check_pstemp_exit():
    """Raise GeneratorError unless the TSCATF cases exercise the early exit
    of PSTEMP (the DEL of the highest L stay zero) and its absence."""
    cases = {case["id"]: case["expected"] for case in json.loads(VALUES.read_text())["cases"]}
    for (case_id, array, first, count), zero in ((_ZERO_DEL, True), (_NONZERO_DEL, False)):
        head = cases[case_id][array]["head"][first:first + count]
        zeros = [bits == "00000000" for bits in head]
        if not (all(zeros) if zero else not any(zeros)):
            raise GeneratorError(f"{case_id}: {array} {head} does not show "
                                 f"{'an early exit' if zero else 'no early exit'}")


def _generate(argv=None):
    """Generate tier1/refcalc_routines/values.json from the baseline."""
    generate(argv, __doc__.splitlines()[0], build_legacy_driver, run_legacy_driver,
             CASES, VALUES, "lib/lib.tleed.f",
             live=lambda case_id: next((arrays for prefix, arrays in LIVE_ARRAYS.items()
                                        if case_id.startswith(prefix)), ()))
    _check_pstemp_exit()
    _check_composite_cases()


if __name__ == "__main__":
    main(_generate)
