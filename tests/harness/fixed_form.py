"""Routines of fixed-form Fortran source text (lib.tleed.f, ref-calc.f)."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import re

# a SUBROUTINE or FUNCTION statement (not a comment line), and a bare END
_ROUTINE_START = re.compile(r"^[ \t]+(?:[\w*]+[ \t]+)*(?:SUBROUTINE|FUNCTION)[ \t]+(\w+)",
                            re.I | re.M)
_END = re.compile(r"^[ \t]+END[ \t]*$", re.I | re.M)


def routine(text, name):
    """The text of routine name in fixed-form source text, from its
    SUBROUTINE/FUNCTION statement to the first bare END after it.

    ValueError unless name is defined exactly once and ends before the next
    routine starts."""
    starts = [m for m in _ROUTINE_START.finditer(text) if m[1].upper() == name]
    if len(starts) != 1:
        raise ValueError(f"{len(starts)} definitions of {name}")
    end = _END.search(text, starts[0].end())
    if not end:
        raise ValueError(f"no END after {name}")
    body = text[starts[0].start():end.end()]
    if len(_ROUTINE_START.findall(body)) != 1:
        raise ValueError(f"{name} runs into the next routine")
    return body
