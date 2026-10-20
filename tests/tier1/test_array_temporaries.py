"""Static check: no compiler-generated array temporaries.

Array temporaries are invisible to the peak-RSS gate unless they coincide
with the peak (an injected ~10 MB temporary of the array REPSTO passes
perf and benchmark). The baseline has none, so every program
is built with -Warray-temporaries and any such warning fails. A probe tree
with a known temporary (array_temporaries_probe/) checks that the same build
and regex do report one.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import re
from pathlib import Path

import pytest

from harness.build import BuildFlags, build_program
from harness.cases import load_cases

WARN_FLAGS = BuildFlags("warn-array-temporaries", "-O2 -Warray-temporaries")
_TEMPORARY = re.compile(r"Warning: Creating array temporary")
_PROBE = Path(__file__).with_name("array_temporaries_probe")


def _temporaries(tl_root, program, param):
    """Build program with WARN_FLAGS; return the reported temporaries with context."""
    exe = build_program(tl_root, program, param, WARN_FLAGS)
    lines = (exe.parent / "build.log").read_text(errors="replace").splitlines()
    return sorted({" | ".join(lines[max(0, i - 4):i + 1]) for i, line in enumerate(lines)
                   if _TEMPORARY.search(line)})


@pytest.mark.tier1
@pytest.mark.parametrize("program", ["refcalc", "rfactor"])
def test_no_array_temporaries(program, tl_root):
    """Building a program with -Warray-temporaries reports none.

    A hit means the code makes the compiler copy an array into a hidden temporary, which
    the peak-memory test alone may not notice."""
    case = load_cases(program)[0]
    hits = _temporaries(tl_root, program, case.runs[0].param)
    assert not hits, f"{len(hits)} array temporaries:\n" + "\n".join(hits[:10])


@pytest.mark.tier1
def test_probe_temporary_is_reported():
    """The check can see a temporary: the probe tree with a known one is reported.

    Failing means the build flags or the warning regex no longer detect temporaries,
    so test_no_array_temporaries would pass vacuously."""
    hits = _temporaries(_PROBE, "rfactor", _PROBE / "PARAM")
    assert any("A(1:N:2)" in h for h in hits), hits
