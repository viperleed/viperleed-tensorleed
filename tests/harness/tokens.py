"""Tokenize text output of Fortran programs into numbers and text.

Fortran fixed-width output frequently concatenates numbers without
separators, omits the leading zero (-.1653E+00), and writes three-digit
exponents without the 'E' (0.123-100). With an 'E', the exponent always
has two digits (0.142878E+020.61E-06 is two numbers; see _FLOAT_E). A field filled with '*' means the
value did not fit the format; this is always an error, also a lone '*'
(a field of width 1). The only exception are the two '*' of an rfactor
beam label such as '*( 1  0)*', which are text.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import re
from typing import Literal, NamedTuple


class FortranOverflowError(ValueError):
    """A field was filled with asterisks: the value did not fit its format."""


class Token(NamedTuple):
    """One number or text item of Fortran output; numbers carry their print format."""
    kind: Literal["float", "int", "text"]
    value: object
    line: int
    ulp: float = 0.0  # floats: unit in the last printed digit
    fixed: bool = False  # floats: printed without exponent (F format)
    fmt: str = ""  # numbers: format as compared (number_format)
    end: int = 0  # numbers: column of the last character (1-based)


def _ulp(s: str) -> float:
    """Unit in the last printed digit of a float literal."""
    m = re.fullmatch(r"[-+]?\d*\.(\d*)(?:[EeDd]?([-+]?\d+))?", s)
    decimals = len(m.group(1))
    exponent = int(m.group(2)) if m.group(2) else 0
    return 10.0 ** (exponent - decimals)


# Ew.d writes 2 exponent digits, so 0.142878E+020.61E-06 is two numbers; a
# 3-digit exponent (1.0E-100) is read only if no digit or '.' follows, a
# 1-digit one (1.5E+5, list-directed style) only if no digit follows.
_FLOAT_E = r"[-+]?(?:\d+\.\d*|\.\d+)[EeDd][-+]?(?:\d{3}(?![\d.])|\d\d|\d(?!\d))"
# E field with |exponent| > 99 written without 'E': 0.1234-100. Not if a '.'
# follows: 0.1234-123.45 is two F fields.
_FLOAT_3EXP = r"[-+]?0?\.\d+[-+]\d{3}(?!\.)"
_FLOAT_PLAIN = r"[-+]?(?:\d+\.\d*|\.\d+)"
_INT = r"[-+]?\d+"
_TOKEN_RE = re.compile(
    rf"(?P<fe>{_FLOAT_E})|(?P<f3>{_FLOAT_3EXP})|(?P<fp>{_FLOAT_PLAIN})"
    rf"|(?P<int>{_INT})|(?P<stars>\*+)|(?P<text>[^\s\d.+\-*]+|[.+\-])"
)


def number_format(s: str, kind: str) -> str:
    """'I'; 'F.<decimals>'; or '<letter>.<decimals>x<exponent width>' for floats.

    The letter is E or D; the letter-less -100 that Ew.d writes for
    |exponent| > 99 counts as E. The exponent width counts letter, sign and
    digits: E+02 and -100 are both 4 wide (the same edit descriptor), E+002
    is 5. The number of integer digits is not part of the format, so -.1234
    and -0.1234 have the same one. A mantissa with a non-zero integer digit
    (-9.79E-02, as ES or 1PE write it) adds S to the letter; a zero value is
    ambiguous and stays E.
    """
    if kind == "int":
        return "I"
    m = re.fullmatch(r"[-+]?(\d*)\.(\d*)(.*)", s)
    if kind == "fp":
        return f"F.{len(m.group(2))}"
    exp = m.group(3)
    letter = "D" if exp[:1] in ("D", "d") else "E"
    if m.group(1).strip("0"):
        letter += "S"
    return f"{letter}.{len(m.group(2))}x{len(exp)}"


# rfactor beam label (ROUT): *( 1  0)*
_BEAM_LABEL_RE = re.compile(r"\*\([^()*]*\)\*")


def _to_float(s: str) -> float:
    s = s.replace("D", "E").replace("d", "E")
    if "E" not in s and "e" not in s:
        m = re.fullmatch(r"([-+]?0?\.\d+)([-+]\d{3})", s)
        if m:
            return float(f"{m.group(1)}E{m.group(2)}")
    return float(s)


def kept_lines(text, ignore_line_patterns=(), skip_lines=()):
    """[(lineno, line)] of the lines not ignored (see tokenize)."""
    ignore = [re.compile(p) for p in ignore_line_patterns]
    skip = set(skip_lines)
    return [(lineno, line) for lineno, line in enumerate(text.splitlines(), start=1)
            if lineno not in skip and not any(p.search(line) for p in ignore)]


def tokenize(text, ignore_line_patterns=(), skip_lines=()):
    """Split text into Token(kind, value, line, ...); skip ignored lines.

    ignore_line_patterns: regexes; a line matching any is skipped.
    skip_lines: 1-based line numbers to skip.
    Raises FortranOverflowError on '*' overflow fields.
    """
    tokens = []
    for lineno, line in kept_lines(text, ignore_line_patterns, skip_lines):
        label_stars = {i for b in _BEAM_LABEL_RE.finditer(line) for i in (b.start(), b.end() - 1)}
        for m in _TOKEN_RE.finditer(line):
            kind = m.lastgroup
            s = m.group()
            if kind == "stars" and s == "*" and m.start() in label_stars:
                kind = "text"
            if kind == "stars":
                raise FortranOverflowError(f"line {lineno}: overflow field {s!r}")
            if kind in ("fe", "f3", "fp"):
                tokens.append(Token("float", _to_float(s), lineno, _ulp(s), kind == "fp",
                                    number_format(s, kind), m.end()))
            elif kind == "int":
                tokens.append(Token("int", int(s), lineno, fmt="I", end=m.end()))
            else:
                tokens.append(Token("text", s, lineno))
    return tokens
