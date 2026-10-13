"""Compare Fortran text outputs token by token within tolerances."""
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .tensor import TensorLayoutError, block_scales, label_tensor
from .tokens import FortranOverflowError, kept_lines, tokenize


@dataclass(frozen=True)
class Tolerance:
    rtol: float
    atol_frac: float  # absolute tolerance as a fraction of max |expected float|
    atol: float = 0.0  # absolute tolerance
    ulps: float = 0.0  # tolerance in units of the expected value's last printed digit
    scale: str = "file"  # atol_frac relative to max |float| of the "file", the "line",
    #                      or "tensor": per WRITE block of a T_* file (see sections)
    sections: tuple = ()  # scale "tensor": ((section, fraction of the block max), ...)
    file_cap: float | None = None  # scaled part <= file_cap * max |float| of the file
    exact_f_decimals: int = -1  # F-format floats with <= this many decimals: exact
    exact_columns: tuple = ()  # 1-based token positions in a line that must be exact;
    #                            "odd"/"even" for all odd/even positions

    def exact_mask(self, tokens):
        """[True for each float token that must match exactly]."""
        mask = []
        col, line = 0, None
        for t in tokens:
            col, line = (col + 1 if t.line == line else 1), t.line
            mask.append(t.kind == "float" and (
                (t.fixed and round(-math.log10(t.ulp)) <= self.exact_f_decimals)
                or col in self.exact_columns
                or ("odd" in self.exact_columns and col % 2 == 1)
                or ("even" in self.exact_columns and col % 2 == 0)))
        return mask

    def allowed(self, expected, ulp, scale):
        """Allowed |a-b| for expected value b (ulp: its last-digit unit)."""
        return self.allowed_abs(expected, ulp, self.atol_frac * scale)

    def allowed_abs(self, expected, ulp, scaled_atol):
        """Allowed |a-b|, with the scaled absolute part already computed."""
        return self.rtol * abs(expected) + scaled_atol + self.atol + self.ulps * ulp


@dataclass
class CompareResult:
    ok: bool
    n_compared: int = 0
    first_failures: list = field(default_factory=list)


_MAX_REPORTED = 10


def compare_texts(actual, expected, tol, ignore_line_patterns=(), skip_lines=()):
    """Compare two texts line by line and token by token, including the number format.

    The texts must have the same number of lines and the same number of
    tokens on each line (after ignored and skipped lines), and each number
    the same format (tokens.number_format) ending in the same column.
    Floats match if |a-b| <= tol.allowed(b, ulp(b), scale), i.e.
    rtol*|b| + atol_frac*max|expected float| (of the file or of b's line)
    + atol + ulps*ulp(b), where ulp(b) is b's unit in the last printed digit;
    ints and text must be identical. A line or token count mismatch, an empty
    expected text, or an overflow field ('***') is a failure.
    """
    res = CompareResult(ok=True)
    try:
        act = tokenize(actual, ignore_line_patterns, skip_lines)
        exp = tokenize(expected, ignore_line_patterns, skip_lines)
    except FortranOverflowError as err:
        res.ok = False
        res.first_failures.append(f"overflow field: {err}")
        return res
    if not exp:
        res.ok = False
        res.first_failures.append("expected output contains no tokens")
        return res
    try:
        scaled = scaled_atols(exp, tol)
        if tol.scale == "tensor" and label_tensor(act) != label_tensor(exp):
            raise TensorLayoutError("tensor layout differs from the expected file")
    except TensorLayoutError as err:
        res.ok = False
        res.first_failures.append(str(err))
        return res
    structure = _line_structure_difference(
        act, kept_lines(actual, ignore_line_patterns, skip_lines),
        exp, kept_lines(expected, ignore_line_patterns, skip_lines))
    if structure:
        res.ok = False
        res.first_failures.append(structure)
        return res
    exact = tol.exact_mask(exp)
    for i, (a, e) in enumerate(zip(act, exp)):
        res.n_compared += 1
        bad = None
        if a.kind != e.kind:
            bad = f"kind {a.kind} vs {e.kind} ({a.value!r} vs {e.value!r})"
        elif e.kind == "float":
            dev = abs(a.value - e.value)
            allowed = 0.0 if exact[i] else tol.allowed_abs(e.value, e.ulp, scaled[i])
            if dev > allowed:
                bad = (f"{a.value!r} vs {e.value!r} "
                       f"(dev {dev:.3e}, allowed {allowed:.3e})")
        elif a.value != e.value:
            bad = f"{a.value!r} vs {e.value!r}"
        if e.kind != "text" and (a.fmt, a.end) != (e.fmt, e.end):
            fmt = (f"format {a.fmt} ending in column {a.end} vs "
                   f"{e.fmt} ending in column {e.end} ({a.value!r} vs {e.value!r})")
            bad = f"{bad}; {fmt}" if bad else fmt
        if bad:
            res.ok = False
            if len(res.first_failures) < _MAX_REPORTED:
                res.first_failures.append(f"line {e.line}: {bad}")
    return res


def _line_structure_difference(act, act_lines, exp, exp_lines):
    """Message naming the first line whose token count differs, or ''.

    With different line counts, that is the first kept line pair whose token
    counts differ, else the first line beyond the shorter file (numbered as
    in the expected file).
    """
    act_counts = Counter(t.line for t in act)
    exp_counts = Counter(t.line for t in exp)
    first = next(((la, le) for (la, _), (le, _) in zip(act_lines, exp_lines)
                  if act_counts[la] != exp_counts[le]), None)
    if len(act_lines) != len(exp_lines):
        if first is None:
            n = len(act_lines)
            le = exp_lines[n][0] if n < len(exp_lines) else exp_lines[-1][0] + 1
        else:
            le = first[1]
        return (f"line {le}: line count differs: actual {len(act_lines)}, "
                f"expected {len(exp_lines)}")
    if first:
        la, le = first
        return (f"line {le}: token count differs: actual {act_counts[la]} "
                f"(actual line {la}), expected {exp_counts[le]}")
    return ""


def scaled_atols(tokens, tol):
    """Scaled absolute part of the tolerance per token (0 for non-floats).

    scale "file"/"line": atol_frac * max |float| of the file or line; floats
    compared exactly (tol.exact_mask: e.g. energies in fd.out) do not count
    for the maximum.
    scale "tensor" (T_* files, parsed by harness.tensor): sections[section]
    * max |float| of the value's WRITE block. Raises TensorLayoutError if
    the file does not have the tensor layout.
    With file_cap, the result is at most file_cap * max |float| of the file
    (all floats): a finer rule capped by a file-scale rule never allows more
    than that rule.
    """
    if tol.scale == "tensor":
        af = dict(tol.sections)
        labels = label_tensor(tokens)
        blocks = block_scales(tokens, labels)
        scaled = [af[lab[0]] * blocks[lab[1]] if lab else 0.0 for lab in labels]
    else:
        scales = float_scales([t for t, x in zip(tokens, tol.exact_mask(tokens)) if not x],
                              tol.scale)
        scaled = [tol.atol_frac * scales[t.line] if t.kind == "float" else 0.0
                  for t in tokens]
    if tol.file_cap is not None:
        cap = tol.file_cap * float_scales(tokens, "file").default
        scaled = [min(x, cap) for x in scaled]
    return scaled


def float_scales(tokens, mode):
    """{line: max |float|} of the line ("line") or of all tokens ("file")."""
    per_line = {}
    for t in tokens:
        if t.kind == "float":
            per_line[t.line] = max(per_line.get(t.line, 0.0), abs(t.value))
    if mode == "line":
        return _Default(per_line)
    if mode != "file":
        raise ValueError(f"unknown tolerance scale {mode!r}")
    return _Default({}, max(per_line.values(), default=0.0))


class _Default(dict):
    def __init__(self, data, default=0.0):
        super().__init__(data)
        self.default = default

    def __missing__(self, key):
        return self.default


def compare_files(actual: Path, expected: Path, tol, ignore_line_patterns=(),
                  skip_lines=()):
    """compare_texts on two files; a missing actual file is a failure."""
    if not Path(actual).is_file():
        return CompareResult(ok=False,
                             first_failures=[f"actual output missing: {actual}"])
    return compare_texts(Path(actual).read_text(errors="replace"),
                         Path(expected).read_text(errors="replace"),
                         tol, ignore_line_patterns, skip_lines)
