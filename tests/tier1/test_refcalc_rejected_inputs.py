"""Tier 1: ref-calc rejects the inputs that its readers check.

The inputs are the converted FINs (TL >= 3.0.0 layout) of the cases; the
tests skip for a tree whose ref-calc does not read that layout.
Each test corrupts one value of the input of ir100_2x1_o/E050_L06 (or of
ag100_psformat2/E045_L06 for the checks of phase-shift format 2), or one size
in its PARAM, so that exactly one check of the reader fires,
runs refcalc and checks that it stops with a non-zero exit code and the
check's message on stderr, and that the message is not on stdout.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

import pytest

from harness.build import DEBUG, RELEASE, build_program
from harness.cases import RUN_TIMEOUT_S, iter_runs
from harness.fin_convert import convert_fin, reads_converted_fin
from harness.param_size import read_param
from harness.runner import run_program

CASE, RUN = next((c, r) for c, r in iter_runs(("refcalc",), only=("ir100_2x1_o",))
                 if r.name == "E050_L06")
CASE_MIX, RUN_MIX = next((c, r) for c, r in iter_runs(("refcalc",), only=("ag100_mix",))
                         if r.name == "E045_L06")
CASE_PS2, RUN_PS2 = next((c, r) for c, r in iter_runs(("refcalc",), only=("ag100_psformat2",))
                         if r.name == "E045_L06")


@pytest.fixture(autouse=True)
def _tree_reads_converted_fin(tl_root):
    if not reads_converted_fin(tl_root):
        pytest.skip("the ref-calc of this tree does not read the TL >= 3.0.0 FIN")


def _converted(run=RUN):
    """The converted FIN text of a run."""
    return convert_fin(run.stdin.read_text(), read_param(run.param))


def _fin(marker, value, occurrence=0, offset=0, run=RUN):
    """Converted FIN text with the integer (first three columns) of a line set to value.

    The line is offset lines after the occurrence-th (from 0) one that
    contains marker."""
    lines = _converted(run).splitlines(keepends=True)
    hits = [i for i, line in enumerate(lines) if marker in line]
    assert len(hits) > occurrence, f"{marker!r} not found {occurrence + 1} times in FIN"
    i = hits[occurrence] + offset
    lines[i] = f"{value:3d}" + lines[i][3:]
    return "".join(lines)


def _record(marker, value, column=0):
    """Converted FIN text with one integer (5 columns) of the record commented marker set to value."""
    lines = _converted().splitlines(keepends=True)
    (i,) = [i for i, line in enumerate(lines) if marker in line]
    start = 5 * column
    lines[i] = lines[i][:start] + f"{value:5d}" + lines[i][start + 5:]
    return "".join(lines)


def _line(marker, text):
    """Converted FIN text with the line containing marker replaced by text."""
    lines = _converted().splitlines(keepends=True)
    (i,) = [i for i, line in enumerate(lines) if marker in line]
    lines[i] = text + "\n"
    return "".join(lines)


def _param(old, new):
    """PARAM text with the one occurrence of old replaced by new."""
    text = RUN.param.read_text()
    assert text.count(old) == 1, f"{old!r} not found once in PARAM"
    return text.replace(old, new)


# (id, FIN text, PARAM text or None, message that must be on stderr)
REJECTED = [
    pytest.param(lambda: _fin("NSITE:", 4), None,
                 "NSITE .ne. MNSITE", id="nsite"),
    pytest.param(lambda: _fin("NLTYPE:", 8), None,
                 "NLTYPE .ne. MNLTYPE", id="nltype"),
    pytest.param(lambda: _fin("number of Bravais sublayers in layer 1", 4), None,
                 "Number of subplanes in layer type 1 exceeds dimension MNSUB",
                 id="nsub-above-mnsub"),
    pytest.param(lambda: _fin("number of Bravais sublayers in layer 2", 0), None,
                 "Number of subplanes in layer type 2 less than 1",
                 id="nsub-zero"),
    # NSITE = 2 differs from NEL = 4 here, so the bound cannot be mistaken
    pytest.param(lambda: _fin("Atom N=2 (Ag)", 3, run=RUN_MIX),
                 lambda: RUN_MIX.param.read_text(),
                 "Site type 3 in sublayer 1 of layer 2: No such site defined",
                 id="site-type-above-nsite"),
    pytest.param(lambda: _fin("Atom N=2 (Ir)", 0), None,
                 "Site type 0 in sublayer 2 of layer 1: No such site defined",
                 id="site-type-zero"),
    pytest.param(_converted,
                 lambda: _param("MNBRAV  = 1", "MNBRAV  = 2"),
                 "Number of Bravais layers in actual input does not match "
                 "dimension MNBRAV", id="nbrav"),
    pytest.param(lambda: _fin("top layer of bulk unit", 0), None,
                 "Top bulk layer is illegal layer type", id="toplayb-zero"),
    pytest.param(lambda: _fin("top layer of bulk unit", 8), None,
                 "Top bulk layer is illegal layer type", id="toplayb-above-nltype"),
    pytest.param(lambda: _fin("bottom layer of bulk unit", 0), None,
                 "Bottom bulk layer is illegal layer type", id="botlayb-zero"),
    pytest.param(lambda: _fin("bottom layer of bulk unit", 8), None,
                 "Bottom bulk layer is illegal layer type", id="botlayb-above-nltype"),
    pytest.param(lambda: _fin("TSLAB =", 2), None,
                 "TSLAB should be 0 or 1", id="tslab"),
    pytest.param(lambda: _fin("NSTACK:", 5), None,
                 "NSTACK .ne. MNSTACK", id="nstack"),
    pytest.param(lambda: _fin("layer 6: layer type 6", 0), None,
                 "Stacked layer number 6 has illegal layer type", id="ltype-zero"),
    pytest.param(lambda: _fin("layer 6: layer type 6", 8), None,
                 "Stacked layer number 6 has illegal layer type", id="ltype-above-nltype"),
    pytest.param(lambda: _fin("LAY = 1: layer type no. 6", 2), None,
                 "Tensor LEED is currently not possible for a substrate type layer",
                 id="tensor-substrate-layer-with-superlattice"),
    pytest.param(lambda: _fin("TENSOR_OUTPUT", 2), None,
                 "Tensor LEED flag for layer 6 not 0 or 1", id="tensor-flag"),
    pytest.param(_converted,
                 lambda: _param("MKNBS = 2", "MKNBS = 3"),
                 "Number of beam sets, MKNBS, differs from number of "
                 "sublattices", id="nl-neq-knbs"),
    pytest.param(_converted,
                 lambda: _param("MKNT =  85", "MKNT =  86"),
                 "Number of beams in full dynamic beam list, 85 is "
                 "inconsistent with dimension MKNT = 86", id="knt"),
    pytest.param(lambda: _record("IDEG NL1 NL2", 3, 0), None,
                 "IDEG = 3 in the input differs from MIDEG = 2 in PARAM", id="ideg"),
    pytest.param(lambda: _record("IDEG NL1 NL2", 4, 1), None,
                 "NL1 = 4 in the input differs from MNL1 = 2 in PARAM", id="nl1"),
    pytest.param(lambda: _record("IDEG NL1 NL2", 2, 2), None,
                 "NL2 = 2 in the input differs from MNL2 = 1 in PARAM", id="nl2"),
    pytest.param(lambda: _record("NPUN (TL", 11), None,
                 "NPUN = 11 in the input differs from MNPUN = 12 in PARAM", id="npun"),
    pytest.param(lambda: _record("NPSI (TL", 43), None,
                 "NPSI = 43 in the input differs from MNPSI = 44 in PARAM", id="npsi"),
    # A record with a value missing is not completed from the next line
    pytest.param(lambda: _line("IDEG NL1 NL2", "    2    2"), None,
                 "The record of IDEG, NL1, NL2 needs 3 integer(s)", id="ideg-record-short"),
    pytest.param(lambda: _line("IDEG NL1 NL2", "    2    2   IDEG NL1, NL2 missing"), None,
                 "The record of IDEG, NL1, NL2 needs 3 integer(s)",
                 id="ideg-record-short-comment"),
    pytest.param(lambda: _line("NPUN (TL", ""), None,
                 "The record of NPUN needs 1 integer(s)", id="npun-record-short"),
    pytest.param(lambda: _line("NPUN (TL", "     NPUN missing"), None,
                 "The record of NPUN needs 1 integer(s)", id="npun-record-short-comment"),
    pytest.param(lambda: _line("NPSI (TL", ""), None,
                 "The record of NPSI needs 1 integer(s)", id="npsi-record-short"),
    pytest.param(lambda: _line("NPSI (TL", "     NPSI missing"), None,
                 "The record of NPSI needs 1 integer(s)", id="npsi-record-short-comment"),
    pytest.param(lambda: _fin("LMAX", 5), None,
                 "Please enforce LMAX in input = MLMAX in PARAM", id="lmax"),
    pytest.param(lambda: _fin("PSFORMAT", 3), None,
                 "Phaseshift format incorrectly assigned", id="psformat"),
    pytest.param(lambda: _fin("IFORM", 4, offset=2), None,
                 "NEL from phaseshift file not consistent with MNEL",
                 id="nel-psformat1"),
    ]


@pytest.mark.tier1
@pytest.mark.parametrize("fin,param,message", REJECTED)
def test_refcalc_rejects_invalid_input(fin, param, message, tl_root, build_flags, tmp_path):
    """refcalc stops with a non-zero exit code and the check's message on stderr.

    Failing means an invalid input is no longer rejected, or is rejected
    with another message or on another channel."""
    _assert_rejected(RUN.param, fin(), param, message, tl_root, build_flags, tmp_path)


@pytest.mark.tier1
@pytest.mark.parametrize("flags", [RELEASE, DEBUG], ids=lambda f: f.name)
def test_refcalc_rejects_ideg_that_fits_no_sublattice(flags, tl_root, tmp_path):
    """IDEG = 3 (in FIN and PARAM) on the p(2x1) lattice of Ir(100)-O is rejected.

    The 120-degree rotation maps the second sublattice onto none, where the
    legacy code indexed past FLMS in FMAT (a trap in debug builds, changed
    output in release builds). Both builds are run because only the debug
    build of the legacy code stops at all. Failing means the check of
    SLIND is gone, or reports on another channel."""
    _assert_rejected(RUN.param, _record("IDEG NL1 NL2", 3, 0),
                     lambda: _param("MIDEG=2", "MIDEG=3"),
                     "IDEG = 3 does not fit", tl_root, flags, tmp_path)


@pytest.mark.tier1
def test_refcalc_rejects_wrong_number_of_elements_in_psformat2(tl_root, build_flags, tmp_path):
    """The check of the number of elements of a phase-shift file of format 2 fires.

    Failing means the check is gone, or it reports on another channel."""
    _assert_rejected(RUN_PS2.param, _fin("NEL", 3, run=RUN_PS2), None,
                     "NEL from phaseshift file not consistent with MNEL",
                     tl_root, build_flags, tmp_path)


def _assert_rejected(param_file, fin, param, message, tl_root, build_flags, tmp_path):
    """Run refcalc on fin (and the PARAM text param(), if any, else param_file).

    Assert a non-zero exit code and, as first text on stderr, the message,
    which must not be on stdout."""
    if param:
        param_file = tmp_path / "PARAM"
        param_file.write_text(param())
    exe = build_program(tl_root, CASE.program, param_file, build_flags)
    (tmp_path / "FIN").write_text(fin)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    res = run_program(exe, run_dir, tmp_path / "FIN", timeout=RUN_TIMEOUT_S)
    assert res.returncode != 0, f"exit code 0; stderr:\n{res.stderr}"
    assert message in res.stderr, res.stderr
    first_line = res.stderr.splitlines()[0]
    assert first_line.startswith("ref-calc: ") and message in first_line, res.stderr
    assert message not in res.stdout


# The values of a record written in another style than five-column fields
RECORD_STYLES = {
    # three-column fields, as the neighbouring LMAX, PSFORMAT and IFORM
    # records (`  4  1  1`)
    "columns": lambda values: "".join(f"{v:3d}" for v in values) + "   ",
    # left-aligned, separated by one blank (`4 1 1 `)
    "left": lambda values: " ".join(str(v) for v in values) + " ",
    # separated by commas (`4,1,1 `)
    "comma": lambda values: ",".join(str(v) for v in values) + " ",
    }


def _reformatted(fin, marker, new_values, style, comment):
    """FIN text with the numbers of the record commented marker rewritten.

    The values are written in style (a key of RECORD_STYLES) and followed
    by comment."""
    lines = fin.splitlines(keepends=True)
    (i,) = [i for i, line in enumerate(lines) if marker in line]
    lines[i] = RECORD_STYLES[style](new_values) + f"{comment}\n"
    return "".join(lines)


@pytest.mark.tier1
@pytest.mark.parametrize("style", RECORD_STYLES)
@pytest.mark.parametrize("comment", ("", "IDEG NL1 NL2, with / a slash"))
def test_refcalc_reads_new_records_list_directed(style, comment, tl_root, build_flags,
                                                 tmp_path):
    """The IDEG/NL1/NL2, NPUN and NPSI records are read for any field width.

    The records written in narrow fields, left-aligned or separated by
    commas, with or without text with separators after them, give the same
    run as the records in five-column fields. Failing means a record is
    read with a fixed format again."""
    param = read_param(RUN.param)
    base = _converted()
    fin = base
    for marker, names in (("IDEG NL1 NL2", ("MIDEG", "MNL1", "MNL2")),
                          ("NPUN (TL", ("MNPUN",)),
                          ("NPSI (TL", ("MNPSI",))):
        fin = _reformatted(fin, marker, [param[n] for n in names], style, comment)
    assert fin != base
    exe = build_program(tl_root, CASE.program, RUN.param, build_flags)
    results = []
    for name, text in (("fixed", base), ("free", fin)):
        run_dir = tmp_path / name
        run_dir.mkdir()
        (run_dir / "FIN").write_text(text)
        results.append(run_program(exe, run_dir, run_dir / "FIN", timeout=RUN_TIMEOUT_S))
    fixed, free = results
    assert fixed.returncode == 0, fixed.stderr
    assert free.returncode == 0, free.stderr
    assert free.stdout == fixed.stdout


@pytest.mark.tier1
def test_refcalc_does_not_check_bravais_subplane_position(tl_root, build_flags, tmp_path):
    """A Bravais layer with a subplane position other than (0,0,0) is not rejected.

    The legacy check for it tests the subplane counter after its loop, where
    it is at least 2, so it never fires; the port keeps the dead check (the
    behaviour is preserved, not fixed). Failing means the check became
    effective."""
    fin = _converted().replace("3   0.0000   0.0000   0.0000  Atom N=14",
                               "3   0.0000   0.0000   0.1000  Atom N=14")
    assert fin != _converted()
    (tmp_path / "FIN").write_text(fin)
    exe = build_program(tl_root, CASE.program, RUN.param, build_flags)
    res = run_program(exe, tmp_path, tmp_path / "FIN", timeout=RUN_TIMEOUT_S)
    assert res.returncode == 0, res.stderr
    assert "Bravais layer" not in res.stderr
