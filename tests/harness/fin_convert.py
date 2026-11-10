"""Conversion of a legacy ref-calc FIN to the TL >= 3.0.0 layout.

The legacy ref-calc takes IDEG, NL1/NL2, NPUN and NPSI from a PARAM file
compiled into the executable. From TL 3.0.0 on, FIN holds them in three
extra records, each next to the data it sizes:
  3I5  IDEG, NL1, NL2  directly after the VPI/ASE/WORKFN record;
  I5   NPUN            directly before the NPU record(s);
  I5   NPSI            directly after the IFORM record.
Text after the numbers is a comment, as elsewhere in FIN. Nothing else
changes: line 2 does not move and no added comment contains '   LMAX', so
viperleed's edit_fin_energy_lmax (viperleed 0.14.x,
calc/sections/refcalc.py:436-455), which edits line 2 and splits at the
first '   LMAX', still works.
convert_fin adds the records to a legacy FIN, strip_fin removes them,
converted_values reads them. reads_converted_fin tells from a TL tree
whether its ref-calc reads the converted layout, and converted_input_file
writes the FIN a run of such a tree is fed.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

import math
from pathlib import Path

from .param_size import _int, _Records, read_param

_COMMENT_COLUMN = 45  # as the LMAX, PSFORMAT and IFORM lines
_RECORDS = (  # (name of the PARAM values, comment)
    (("MIDEG", "MNL1", "MNL2"), "IDEG NL1 NL2 (TL >= 3.0.0)"),
    (("MNPUN",), "NPUN (TL >= 3.0.0)"),
    (("MNPSI",), "NPSI (TL >= 3.0.0)"),
)


def _positions(fin_text, param_values=None):
    """Line indices of the three records, and the NL1*NL2 used to find them.

    FIN is read in READIN's order. For a legacy FIN (param_values given)
    the indices are those the records would be inserted at; for a converted
    one (param_values None) those the records are at.
    """
    converted = param_values is None
    rec = _Records(fin_text)
    rec.skip(7)  # title, EI/EF/DE, ARA1, ARA2, ARB1, ARB2, VPI/ASE/WORKFN
    first = rec.pos
    if converted:
        line = rec.next()
        nsets = _int(line[5:], 5) * _int(line[10:], 5)
    else:
        nsets = param_values["MNL1"] * param_values["MNL2"]
    for _ in range(nsets):
        rec.skip(_int(rec.next(), 5))
    rec.next()  # TST
    second = rec.pos
    if converted:
        npun = _int(rec.next(), 5)
    else:
        npun = param_values["MNPUN"]
    rec.skip(math.ceil(npun / 15) + 3 + 3)  # NPU, THETA/FI, EPS, LITER, LMAX, PSFMT, IFORM
    return first, second, rec.pos


def _record(values, comment, ending):
    return "".join(f"{v:5d}" for v in values).ljust(_COMMENT_COLUMN) + comment + ending


def convert_fin(fin_text, param_values):
    """The TL >= 3.0.0 FIN for a legacy FIN and the PARAM values of its run.

    param_values is a {name: value} dict with at least MIDEG, MNL1, MNL2,
    MNPUN and MNPSI. Line endings and every other byte are preserved; the
    new records use the line ending of the line they follow.
    """
    lines = fin_text.splitlines(keepends=True)
    for index, (names, comment) in zip(reversed(_positions(fin_text, param_values)),
                                       reversed(_RECORDS)):
        previous = lines[index - 1]
        ending = previous[len(previous.rstrip("\r\n")):] or "\n"
        lines.insert(index, _record([param_values[n] for n in names], comment, ending))
    return "".join(lines)


def strip_fin(fin_text):
    """The legacy FIN of a converted one: the inverse of convert_fin."""
    lines = fin_text.splitlines(keepends=True)
    for index in reversed(_positions(fin_text)):
        del lines[index]
    return "".join(lines)


def converted_values(fin_text):
    """{name: value} of the PARAM values held by the records of a converted FIN."""
    lines = fin_text.splitlines()
    values = {}
    for index, (names, _) in zip(_positions(fin_text), _RECORDS):
        for i, name in enumerate(names):
            values[name] = _int(lines[index][5*i:], 5)
    return values


def reads_converted_fin(tl_root):
    """True if the ref-calc of the tree reads the TL >= 3.0.0 FIN.

    The ported reader (lib/leed_input.f90) reads only that layout; the
    legacy ref-calc, which takes the values from PARAM, reads only the
    original.
    """
    return (Path(tl_root) / "TensErLEED" / "lib" / "leed_input.f90").is_file()


def converted_input_file(run, workdir):
    """Write the converted FIN of a run (values from its PARAM) to workdir/.stdin; return the path."""
    path = Path(workdir) / ".stdin"
    path.write_text(convert_fin(Path(run.stdin).read_text(), read_param(run.param)))
    return path
