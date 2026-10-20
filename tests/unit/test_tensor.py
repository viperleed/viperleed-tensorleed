"""Structural T_* parser and the section-aware T_* tolerance."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

from collections import Counter

import pytest

from harness.cases import CASES_DIR, tolerance_for
from harness.compare import Tolerance, compare_texts, scaled_atols
from harness.tensor import TensorLayoutError, label_tensor
from harness.tokens import tokenize

# ag100, 1 energy, L1 = 7: 5 exit beams
T_FILE = CASES_DIR / "refcalc" / "ag100" / "expected" / "E045_L06" / "T_1"
# ir100_2x1_o, 4 energies (110-122 eV)
T_MULTI = CASES_DIR / "refcalc" / "ir100_2x1_o" / "expected" / "E110-122_L07" / "T_1"


def _labels(path):
    toks = tokenize(path.read_text())
    return toks, label_tensor(toks)


def test_parser_labels_a_real_tensor_file():
    """Floats of a real T_* file get their section label, integers none.

    The section sizes follow from L1 and the number of exit beams.
    """
    toks, labels = _labels(T_FILE)
    sections = Counter(lab[0] for lab in labels if lab)
    l1 = toks[4].value  # line 2
    n_exit = sum(1 for t, lab in zip(toks, labels) if lab is None) - 2  # minus L1 and -1
    assert l1 == 7
    assert sections["hdr"] == 4
    assert sections["tmat"] == 2 * l1
    assert sections["kvec"] == 4 * n_exit
    assert sections["alm"] == 2 * l1 ** 2 * n_exit
    assert sections["xist"] > 0 and sections["xist"] % 2 == 0
    # every float is labelled, every int is not
    assert all((lab is None) == (t.kind == "int") for t, lab in zip(toks, labels))


def test_parser_handles_several_energies():
    """A file with several energies is labelled energy by energy."""
    toks, labels = _labels(T_MULTI)
    assert sum(1 for t in toks if t.kind == "int" and t.value == -1) == 4
    assert Counter(lab[0] for lab in labels if lab)["hdr"] == 4 * 4


@pytest.mark.parametrize("edit", [
    lambda lines: lines[:3] + lines[4:],                 # one CAF line missing
    lambda lines: [lines[0] + "0.1000000000E+01"] + lines[1:],  # 5 header values
    lambda lines: lines[:-1],                            # terminating -1 missing
    lambda lines: lines[:1] + ["    0"] + lines[2:],     # L1 = 0
    lambda lines: [ln.replace("    0", "   -2") for ln in lines],  # bad exit number
])
def test_layout_mismatch_raises(edit):
    """Each listed corruption of the T_* layout raises TensorLayoutError."""
    lines = T_FILE.read_text().splitlines()
    with pytest.raises(TensorLayoutError):
        label_tensor(tokenize("\n".join(edit(lines))))


def test_layout_mismatch_fails_comparison():
    """A moved line break fails the comparison instead of shifting the block labels.

    The tokens stay the same, so only the layout check can detect it.
    """
    text = T_FILE.read_text()
    lines = text.splitlines()
    moved = "\n".join(lines[:2] + [lines[2] + lines[3]] + lines[4:]) + "\n"
    tol, *_ = tolerance_for("T_1")
    assert len(tokenize(moved)) == len(tokenize(text))
    res = compare_texts(moved, text, tol)
    assert not res.ok and "layout" in res.first_failures[0]


def _perturb(text, section, factor, index=0):
    """Scale the index-th float of the first block of section by factor."""
    toks = tokenize(text)
    labels = label_tensor(toks)
    target = [t for t, lab in zip(toks, labels) if lab and lab[0] == section][index]
    lines = text.splitlines()
    line = lines[target.line - 1]
    field = _fortran_e(target.value)
    assert field in line, (field, line)
    lines[target.line - 1] = line.replace(field, _fortran_e(target.value * factor), 1)
    return "\n".join(lines) + "\n"


def _fortran_e(x, d=10):
    """x in Fortran E16.10 form (0.1234567890E+01, -.1234567890E+01)."""
    if x == 0:
        return f"0.{'0' * d}E+00"
    s = f"{abs(x):.{d - 1}E}"
    mant, exp = s.split("E")
    digits = mant.replace(".", "")
    e = int(exp) + 1
    return f"{'-' if x < 0 else '0'}.{digits}E{'-' if e < 0 else '+'}{abs(e):02d}"


# the file-scale rule (the cap of the T_* rule): one absolute part for the
# whole file (1.9e-4 * file max)
FILE_RULE = Tolerance(rtol=1e-6, atol_frac=1.9e-4, ulps=2)


def test_perturbed_xist_value_is_caught():
    """A 1 % change of a mid-sized XIST component fails only the section-aware rule.

    The file-scale rule would accept it.
    """
    text = T_FILE.read_text()
    tol, *_ = tolerance_for("T_1")
    perturbed = _perturb(text, "xist", 1.01, 2)
    assert compare_texts(text, text, tol).ok
    assert compare_texts(perturbed, text, FILE_RULE).ok
    assert not compare_texts(perturbed, text, tol).ok


def test_perturbed_tmat_value_is_caught():
    """A 1e-4 change of the dominant t-matrix element fails only the section-aware rule.

    The file-scale rule would accept it.
    """
    text = T_FILE.read_text()
    tol, *_ = tolerance_for("T_1")
    perturbed = _perturb(text, "tmat", 1 + 1e-4, 1)
    assert compare_texts(perturbed, text, FILE_RULE).ok
    assert not compare_texts(perturbed, text, tol).ok


def test_block_rule_never_looser_than_file_rule():
    """The section-aware tolerance never exceeds the file-scale rule in any block.

    No existing tolerance is loosened.
    """
    tol, *_ = tolerance_for("T_1")
    toks = tokenize(T_FILE.read_text())
    file_rule = scaled_atols(toks, Tolerance(tol.rtol, tol.file_cap, ulps=tol.ulps))
    assert all(n <= f for n, f in zip(scaled_atols(toks, tol), file_rule))
