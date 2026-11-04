"""The ref-calc design document covers every legacy routine and COMMON block.

TensErLEED/doc/refcalc-design.md maps each program unit (PROGRAM, SUBROUTINE,
FUNCTION, ENTRY) and each COMMON block of src/ref-calc.f and lib/lib.tleed.f
to its place in the free-form port. These tests take the names from the
Fortran sources, so a routine or block the document does not mention fails.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import re

import pytest

LEGACY_SOURCES = ("src/ref-calc.f", "lib/lib.tleed.f")
DESIGN_DOC = "doc/refcalc-design.md"

# Fixed form: a comment has C, c, * or ! in column 1; statements may start
# with blanks or a tab.
_UNIT = re.compile(r"^[ \t]+(?:[\w*]+[ \t]+)*(?:PROGRAM|SUBROUTINE|FUNCTION|ENTRY)"
                   r"[ \t]+(\w+)", re.I | re.M)
_COMMON = re.compile(r"^[ \t]+COMMON[ \t]*(/[ \t]*\w*[ \t]*/)?", re.I | re.M)


def _legacy_text(tl_root):
    """The text of the legacy ref-calc sources, without comment lines."""
    tl = tl_root / "TensErLEED"
    lines = []
    for source in LEGACY_SOURCES:
        lines += [line for line in (tl / source).read_text(errors="replace").splitlines()
                  if line[:1] not in "Cc*!"]
    return "\n".join(lines)


def _routine_names(tl_root):
    """The names of all program units and entries of the legacy sources."""
    return sorted({name.upper() for name in _UNIT.findall(_legacy_text(tl_root))})


def _common_names(tl_root):
    """The COMMON blocks of the legacy sources: '/NAME/', or 'blank COMMON'."""
    names = set()
    for block in _COMMON.findall(_legacy_text(tl_root)):
        name = re.sub(r"[ \t]", "", block).upper()
        names.add(name if name not in ("", "//") else "blank COMMON")
    return sorted(names)


@pytest.fixture(name="design_doc")
def fixture_design_doc(tl_root):
    """The text of the design document."""
    return (tl_root / "TensErLEED" / DESIGN_DOC).read_text()


def test_legacy_names_found(tl_root):
    """The source scan finds the main program, a function and blank COMMON."""
    assert {"LEED", "FACT", "TSCATF"} <= set(_routine_names(tl_root))
    assert {"blank COMMON", "/SL/"} <= set(_common_names(tl_root))


def test_every_routine_in_design_doc(tl_root, design_doc):
    """Every routine of ref-calc.f and lib.tleed.f appears in the document as a word."""
    missing = [name for name in _routine_names(tl_root)
               if not re.search(rf"\b{name}\b", design_doc)]
    assert not missing, f"{DESIGN_DOC} does not mention {missing}"


def test_every_common_block_in_design_doc(tl_root, design_doc):
    """Every COMMON block of ref-calc.f and lib.tleed.f appears in the document."""
    missing = [name for name in _common_names(tl_root) if name not in design_doc]
    assert not missing, f"{DESIGN_DOC} does not mention {missing}"
