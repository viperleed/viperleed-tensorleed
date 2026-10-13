"""Structural parser for formatted TensErLEED tensor files (T_*).

Layout per energy, as written by OUTXIST and OUTAMP (lib.tleed.f):
    E, VPIS, VPIO, VV                  hdr   one line, 4 floats (5E16.10)
    L1                                       one line, 1 int (I5)
    CAF(1:L1)                          tmat  2*L1 floats, 5 per line
    XIST(1:NT0)                        xist  2*NT0 floats, 5 per line
    for each exit beam:
      NEXIT                                  one line, 1 int >= 0
      PQ1, PQ2, AK21, AK31             kvec  one line, 4 floats (5E12.6)
      ALM(1:L1**2)                     alm   2*L1**2 floats, 5 per line
    -1                                       one line, end of energy
Every WRITE starts a new line, so each block starts on its own line. Any
deviation from this layout raises TensorLayoutError: a changed layout must
fail the comparison, never shift tokens into the wrong section.
"""
from itertools import groupby

_PER_LINE = 5


class TensorLayoutError(ValueError):
    pass


def label_tensor(tokens):
    """[(section, block) for float tokens, None for int tokens].

    tokens: harness.tokens.tokenize() of a whole T_* file. block numbers
    the WRITE blocks of the file consecutively.
    """
    lines = [(ln, [i for i, _ in grp])
             for ln, grp in groupby(enumerate(tokens), key=lambda it: it[1].line)]
    labels = [None] * len(tokens)
    pos = 0
    block = -1

    def fail(msg):
        line = lines[pos][0] if pos < len(lines) else "EOF"
        raise TensorLayoutError(f"tensor layout, line {line}: {msg}")

    def kinds(idx):
        return {tokens[i].kind for i in idx}

    def single_int():
        nonlocal pos
        if pos >= len(lines):
            fail("expected an integer line")
        _, idx = lines[pos]
        if len(idx) != 1 or tokens[idx[0]].kind != "int":
            fail("expected a line with one integer")
        pos += 1
        return tokens[idx[0]].value

    def floats(section, n=None):
        """Consume one WRITE block of floats: exactly n, or up to the next int line."""
        nonlocal pos, block
        block += 1
        count = 0
        while pos < len(lines) and (n is None or count < n):
            _, idx = lines[pos]
            if kinds(idx) != {"float"}:
                if n is None and kinds(idx) == {"int"}:
                    break
                fail(f"expected only floats in section {section}")
            if len(idx) > _PER_LINE or (count and prev_len != _PER_LINE):
                fail(f"section {section}: not {_PER_LINE} values per line")
            for i in idx:
                labels[i] = (section, block)
            count += len(idx)
            prev_len = len(idx)
            pos += 1
        if n is not None and count != n:
            fail(f"section {section}: {count} values, expected {n}")
        if count == 0 or count % 2:
            fail(f"section {section}: {count} values, expected an even number > 0")
        return count

    def one_line(section, n):
        nonlocal pos, block
        if pos >= len(lines):
            fail(f"expected section {section}")
        _, idx = lines[pos]
        if kinds(idx) != {"float"} or len(idx) != n:
            fail(f"section {section}: expected one line of {n} floats")
        block += 1
        for i in idx:
            labels[i] = (section, block)
        pos += 1

    n_energies = 0
    while pos < len(lines):
        one_line("hdr", 4)
        l1 = single_int()
        if l1 < 1:
            fail(f"L1 = {l1}")
        floats("tmat", 2 * l1)
        floats("xist")
        while (nexit := single_int()) != -1:
            if nexit < 0:
                fail(f"exit beam number {nexit}")
            one_line("kvec", 4)
            floats("alm", 2 * l1 ** 2)
        n_energies += 1
    if not n_energies:
        raise TensorLayoutError("tensor layout: empty file")
    return labels


def block_scales(tokens, labels):
    """{block: max |float| of the block}."""
    scales = {}
    for t, lab in zip(tokens, labels):
        if lab:
            scales[lab[1]] = max(scales.get(lab[1], 0.0), abs(t.value))
    return scales
