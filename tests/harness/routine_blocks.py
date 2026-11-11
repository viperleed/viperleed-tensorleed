"""Summaries of the array blocks that the ref-calc routine-check drivers print.

A block is a header line "name kind count" followed by the values: kind r4
and r8 as the hex digits of the bit patterns (8 and 16 per value), kind i4
in decimal; 8 values per line (4 for r8). The summary of a block is its
count, the SHA-256 of its values (lower-case hex tokens or decimals, joined
by single spaces) and its first HEAD values.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import hashlib

HEAD = 16
_HEX_WIDTH = {"r4": 8, "r8": 16}
_PER_LINE = {"r4": 8, "r8": 4, "i4": 8}


def summarize(stdout):
    """Return {name: {"kind", "n", "sha256", "head"}} for the blocks in stdout."""
    summary = {}
    lines = iter(stdout.splitlines())
    for header in lines:
        name, kind, count = header.split()
        count = int(count)
        if kind in _HEX_WIDTH:
            width = _HEX_WIDTH[kind]
            body = "".join(line.strip() for line in _body(lines, count, kind))
            values = [body[i:i + width].lower() for i in range(0, len(body), width)]
        else:
            values = [token for line in _body(lines, count, kind)
                      for token in line.split()]
        if len(values) != count:
            raise ValueError(f"block {name}: {len(values)} values, expected {count}")
        summary[name] = {"kind": kind, "n": count,
                         "sha256": hashlib.sha256(" ".join(values).encode()).hexdigest(),
                         "head": values[:HEAD]}
    return summary


def _body(lines, count, kind):
    """The lines of one block of count values of the given kind."""
    for _ in range(-(-count // _PER_LINE[kind])):
        yield next(lines)
