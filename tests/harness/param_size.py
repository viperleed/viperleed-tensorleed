"""Array sizes of a ref-calc run (the PARAM values) derived from its input.

The legacy ref-calc takes its array sizes from a PARAM file compiled into
the executable. Each size is one of three kinds (SIZES, the column 'kind'):
  strict:   READIN/READGEO (lib/lib.tleed.f) stop unless the value equals
            the count in the input;
  bound:    the value only needs to be at least the count in the input;
  explicit: the value reaches ref-calc only through PARAM and cannot be
            derived from the input (IDEG, NL1/NL2, NPUN, NPSI); it is taken
            from the run's PARAM.
derive_sizes reads the stdin input of ref-calc ("FIN", which also holds the
phase shifts) with the explicit values and returns every size, bounds at
their minimum. param_text writes a PARAM file as viperleed's writePARAM
does (viperleed 0.14.2, calc/files/iorefcalc.py:234-336), including the
derived quantities of its section 4.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import math
import re
from pathlib import Path
from typing import NamedTuple

STRICT, BOUND, EXPLICIT = "strict", "bound", "explicit"


class Size(NamedTuple):
    """One PARAM value: its name, kind and how it follows from the input."""
    name: str
    kind: str
    source: str


SIZES = (
    Size("MIDEG", EXPLICIT, "PARAM: degree of the symmetry used in SLIND/FMAT; not in FIN"),
    Size("MNL1", EXPLICIT, "PARAM: superlattice repetitions along a1; not in FIN"),
    Size("MNL2", EXPLICIT, "PARAM: superlattice repetitions along a2; not in FIN"),
    Size("MKNBS", STRICT, "FIN beam list: number of beam sets; READIN stops unless "
                          "it equals MNL1*MNL2, so MNL1*MNL2 sets are read"),
    Size("MKNT", STRICT, "FIN beam list: sum of the beam counts of the sets; "
                         "READIN stops unless equal"),
    Size("MNPUN", EXPLICIT, "PARAM: number of output beams (MNT0 = MNPUN); FIN "
                            "lists them without a count"),
    Size("MNPSI", EXPLICIT, "PARAM: number of phase-shift energies; FIN gives "
                            "no count"),
    Size("MNEL", STRICT, "FIN phase-shift block: NEL on its first line; READIN "
                         "stops unless equal"),
    Size("MLMAX", STRICT, "FIN: LMAX line; READIN stops unless LMAX = MLMAX"),
    Size("MNLMO", BOUND, "FIN: NCA(LMAX), the size of CAA used (ref-calc.f DATA "
                         "NCA); no check"),
    Size("MNLM", BOUND, "FIN: NLMS(LMAX), the size of CLM used (ref-calc.f DATA "
                        "NLMS); no check"),
    Size("MNSITE", STRICT, "FIN geometry: NSITE; READGEO stops unless equal"),
    Size("MNLTYPE", STRICT, "FIN geometry: NLTYPE; READGEO stops unless equal"),
    Size("MNBRAV", STRICT, "FIN geometry: layer types with one subplane; READGEO "
                           "stops unless equal"),
    Size("MNSUB", BOUND, "FIN geometry: largest number of subplanes of a layer "
                         "type; READGEO stops only if a layer type has more"),
    Size("MNSTACK", STRICT, "FIN geometry: NSTACK; READGEO stops unless equal"),
)
KIND = {s.name: s.kind for s in SIZES}

# Sizes of CAA (NCA) and CLM (NLMS) for LMAX = 1..18, as in viperleed's
# writePARAM (iorefcalc.py:247-250); ref-calc.f's DATA NLMS has 0 for LMAX = 1.
MNLMO_BY_LMAX = (1, 70, 264, 759, 1820, 3836, 7344, 13053, 21868, 34914, 53560,
                 79443, 114492, 160952, 221408, 298809, 396492, 518206)
MNLM_BY_LMAX = (1, 76, 284, 809, 1925, 4032, 7680, 13593, 22693, 36124, 55276,
                81809, 117677, 165152, 226848, 305745, 405213, 529036)

_ASSIGNMENT = re.compile(r"\b(M[A-Z0-9]*)\s*=\s*(\d+)\s*[,)]")


def read_param(path):
    """{name: value} of the integer constants assigned in a PARAM file."""
    text = Path(path).read_text()
    return {name: int(value) for name, value in _ASSIGNMENT.findall(text)}


def explicit_values(param):
    """The explicit values of a {name: value} PARAM dict."""
    return {s.name: param[s.name] for s in SIZES if s.kind == EXPLICIT}


def _int(line, width):
    """Integer in the first width columns of line, read as Fortran Iw (blank is 0)."""
    field = line[:width].strip()
    return int(field) if field else 0


class _Records:
    """The lines of a FIN, read one record at a time."""

    def __init__(self, text):
        self.lines = text.splitlines()
        self.pos = 0

    def next(self):
        if self.pos >= len(self.lines):
            raise ValueError("FIN ends early")
        self.pos += 1
        return self.lines[self.pos - 1]

    def skip(self, n):
        for _ in range(n):
            self.next()


def _read_record(rec, names, width, explicit):
    """Read a record of Iwidth integers into explicit under names (a converted FIN's)."""
    line = rec.next()
    for i, name in enumerate(names):
        explicit[name] = _int(line[i*width:], width)


def derive_sizes(fin_text, explicit):
    """{name: value} for every Size, from a FIN and the explicit values.

    The FIN is read in the order and with the record layout of READIN and
    READGEO; explicit ({name: value} of the EXPLICIT sizes) supplies the
    record counts the FIN does not give (beam sets, output beams, phase-shift
    energies). With explicit=None the FIN is a converted one (TL >= 3.0.0
    layout, see fin_convert.py) and holds these values itself. Bounds are at
    their minimum. ValueError if the FIN does not have the layout expected,
    or for LMAX < 2, where no size is correct (ref-calc writes past CAA at
    LMAX = 1).
    """
    from_fin = explicit is None
    explicit = {} if from_fin else dict(explicit)
    rec = _Records(fin_text)
    rec.skip(7)  # title, EI/EF/DE, ARA1, ARA2, ARB1, ARB2, VPI/ASE/WORKFN
    if from_fin:
        _read_record(rec, ("MIDEG", "MNL1", "MNL2"), 5, explicit)
    sizes = dict(explicit)
    sizes["MKNBS"] = explicit["MNL1"] * explicit["MNL2"]
    knt = 0
    for _ in range(sizes["MKNBS"]):
        knb = _int(rec.next(), 5)
        rec.skip(knb)
        knt += knb
    sizes["MKNT"] = knt
    tst = rec.next()
    if "." not in tst[:11]:
        raise ValueError(f"FIN: expected TST after {sizes['MKNBS']} beam sets, got {tst!r}")
    if from_fin:
        _read_record(rec, ("MNPUN",), 5, explicit)
        sizes["MNPUN"] = explicit["MNPUN"]
    rec.skip(math.ceil(explicit["MNPUN"] / 15))  # NPU, format 15I5
    rec.skip(3)  # THETA/FI, EPS, LITER
    lmax = _int(rec.next(), 3)
    if lmax < 2:
        # legacy bug (README, known issues): CAAA writes CAA(2) past NCA(1) = 1
        raise ValueError(f"FIN: LMAX = {lmax}; ref-calc needs LMAX >= 2: at LMAX = 1 "
                         "it writes past CAA for any MNLMO, so no PARAM size is correct")
    sizes["MLMAX"] = lmax
    sizes["MNLMO"] = MNLMO_BY_LMAX[lmax - 1]
    sizes["MNLM"] = MNLM_BY_LMAX[lmax - 1]
    psfmt = _int(rec.next(), 3)
    rec.skip(1)  # IFORM
    if from_fin:
        _read_record(rec, ("MNPSI",), 5, explicit)
        sizes["MNPSI"] = explicit["MNPSI"]
    nel = _int(rec.next(), 3)
    sizes["MNEL"] = nel
    npsi = explicit["MNPSI"]
    if psfmt == 1:
        # per energy: ES, then per element LMAX+1 values in 10F7.4, and
        # one more record skipped if LMAX < 10
        per_element = math.ceil((lmax + 1) / 10) + (1 if lmax < 10 else 0)
        rec.skip(npsi * (1 + nel * per_element))
    elif psfmt == 2:
        rec.skip(2)  # PAR0/PAR1, PSL1 and coefficients
        rec.skip(nel * (2 + npsi))  # per element: two headers, one record per energy
    else:
        raise ValueError(f"FIN: phase-shift format {psfmt} is neither 1 nor 2")

    rec.skip(3)
    nsite = _int(rec.next(), 3)
    sizes["MNSITE"] = nsite
    rec.skip(nsite * (1 + nel))
    rec.skip(3)
    nltype = _int(rec.next(), 3)
    sizes["MNLTYPE"] = nltype
    nsub = []
    for _ in range(nltype):
        rec.skip(2)  # header, LATT
        nsub.append(_int(rec.next(), 3))
        rec.skip(nsub[-1])
    sizes["MNBRAV"] = nsub.count(1)
    sizes["MNSUB"] = max(nsub)
    rec.skip(3)
    rec.skip(5)  # TSLAB and four records (ASA, TOPLAYB, BOTLAYB, ASBULK or skipped)
    rec.skip(3)
    nstack = _int(rec.next(), 3)
    sizes["MNSTACK"] = nstack
    for _ in range(nstack):
        ltype = _int(rec.next(), 3)
        if _int(rec.next(), 3) == 1:  # TENS: one file name per subplane
            rec.skip(nsub[ltype - 1])
    return sizes


def param_text(sizes):
    """PARAM file text for sizes ({name: value} of all SIZES), as writePARAM writes it.

    Layout and the section-4 formulas follow viperleed's writePARAM
    (iorefcalc.py:253-336): MNBRAV2 is 1 without Bravais layers, and the
    composite-layer sizes (MNCOMP ... MLMNI) are 1 if MNSUB is 1.
    """
    s = sizes
    composite = s["MNSUB"] > 1

    def comp(formula):
        return formula if composite else "1 "

    return (
        "C  Dimension statements for Tensor LEED reference calculation, \n"
        "C  version v1.2\n\n"
        "C  1. lattice symmetry\n\n"
        f"      PARAMETER (MIDEG={s['MIDEG']},MNL1={s['MNL1']},MNL2={s['MNL2']})\n"
        "      PARAMETER (MNL = MNL1*MNL2)\n"
        "\nC  2. General calculational quantities\n\n"
        f"      PARAMETER (MKNBS = {s['MKNBS']})\n"
        f"      PARAMETER (MKNT =  {s['MKNT']})\n"
        f"      PARAMETER (MNPUN = {s['MNPUN']}, MNT0 = {s['MNPUN']})\n"
        f"      PARAMETER (MNPSI = {s['MNPSI']}, MNEL = {s['MNEL']})\n"
        f"      PARAMETER (MLMAX = {s['MLMAX']})\n"
        f"      PARAMETER (MNLMO = {s['MNLMO']}, MNLM = {s['MNLM']})\n"
        "\nC  3. Parameters for (3D) geometry within (2D) unit mesh\n\n"
        f"      PARAMETER (MNSITE  = {s['MNSITE']})\n"
        f"      PARAMETER (MNLTYPE = {s['MNLTYPE']})\n"
        f"      PARAMETER (MNBRAV  = {s['MNBRAV']})\n"
        f"      PARAMETER (MNSUB   = {s['MNSUB']})\n"
        f"      PARAMETER (MNSTACK = {s['MNSTACK']})\n"
        "\nC  4. some derived quantities that must be treated explicitly "
        "(dummy dimensions for\n"
        "C     special cases necessary\n\n"
        "      PARAMETER (MLMAX1=MLMAX+1)\n"
        "      PARAMETER (MLMMAX = MLMAX1*MLMAX1)\n\n"
        f"      PARAMETER (MNBRAV2 = {'MNBRAV' if s['MNBRAV'] > 0 else '1'})\n\n"
        f"      PARAMETER (MNCOMP= {comp('MNLTYPE-MNBRAV')})\n"
        f"      PARAMETER (MLMT  = {comp('MNSUB*MLMMAX')})\n"
        f"      PARAMETER (MNSUB2= {comp('MNSUB * (MNSUB-1)/2')})\n"
        f"      PARAMETER (MLMG  = {comp('MNSUB2*MLMMAX*2')})\n"
        f"      PARAMETER (MLMN  = {comp('MNSUB * MLMMAX')})\n"
        f"      PARAMETER (MLM2N = {comp('2*MLMN')})\n"
        f"      PARAMETER (MLMNI = {comp('MNSUB*MLMMAX')})\n"
    )


def size_mismatches(derived, param):
    """Messages for each strict size not equal to, and bound size above, its PARAM value."""
    problems = []
    for s in SIZES:
        if s.kind == STRICT and derived[s.name] != param[s.name]:
            problems.append(f"{s.name}: PARAM {param[s.name]}, input {derived[s.name]}")
        elif s.kind == BOUND and derived[s.name] > param[s.name]:
            problems.append(f"{s.name}: PARAM {param[s.name]} < input {derived[s.name]}")
    return problems


def oversized(derived, param):
    """{name: (PARAM value, derived value)} of the bound sizes PARAM sets above the minimum."""
    return {s.name: (param[s.name], derived[s.name]) for s in SIZES
            if s.kind == BOUND and param[s.name] > derived[s.name]}


def run_sizes(run):
    """(PARAM values, derived values) of a ref-calc case run."""
    param = read_param(run.param)
    return param, derive_sizes(Path(run.stdin).read_text(), explicit_values(param))
