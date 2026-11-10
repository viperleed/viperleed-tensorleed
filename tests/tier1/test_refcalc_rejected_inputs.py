"""Tier 1: ref-calc rejects the inputs that its readers check.

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

from harness.build import build_program
from harness.cases import RUN_TIMEOUT_S, iter_runs
from harness.runner import run_program

CASE, RUN = next((c, r) for c, r in iter_runs(("refcalc",), only=("ir100_2x1_o",))
                 if r.name == "E050_L06")
CASE_PS2, RUN_PS2 = next((c, r) for c, r in iter_runs(("refcalc",), only=("ag100_psformat2",))
                         if r.name == "E045_L06")


def _fin(marker, value, occurrence=0, offset=0, run=RUN):
    """FIN text with the integer (first three columns) of a line set to value.

    The line is offset lines after the occurrence-th (from 0) one that
    contains marker."""
    lines = run.stdin.read_text().splitlines(keepends=True)
    hits = [i for i, line in enumerate(lines) if marker in line]
    assert len(hits) > occurrence, f"{marker!r} not found {occurrence + 1} times in FIN"
    i = hits[occurrence] + offset
    lines[i] = f"{value:3d}" + lines[i][3:]
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
    pytest.param(lambda: _fin("Atom N=2 (Ir)", 4), None,
                 "Site type 4 in sublayer 2 of layer 1: No such site defined",
                 id="site-type-above-nsite"),
    pytest.param(lambda: _fin("Atom N=2 (Ir)", 0), None,
                 "Site type 0 in sublayer 2 of layer 1: No such site defined",
                 id="site-type-zero"),
    pytest.param(lambda: RUN.stdin.read_text(),
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
    pytest.param(lambda: RUN.stdin.read_text(),
                 lambda: _param("MKNBS = 2", "MKNBS = 3"),
                 "Number of beam sets, MKNBS, differs from number of "
                 "sublattices", id="nl-neq-knbs"),
    pytest.param(lambda: RUN.stdin.read_text(),
                 lambda: _param("MKNT =  85", "MKNT =  86"),
                 "Number of beams in full dynamic beam list, 85 is "
                 "inconsistent with dimension MKNT = 86", id="knt"),
    pytest.param(lambda: _fin("LMAX", 5), None,
                 "Please enforce LMAX in input = MLMAX in PARAM", id="lmax"),
    pytest.param(lambda: _fin("PSFORMAT", 3), None,
                 "Phaseshift format incorrectly assigned", id="psformat"),
    pytest.param(lambda: _fin("IFORM", 4, offset=1), None,
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
def test_refcalc_rejects_wrong_number_of_elements_in_psformat2(tl_root, build_flags, tmp_path):
    """The check of the number of elements of a phase-shift file of format 2 fires.

    Failing means the check is gone, or it reports on another channel."""
    _assert_rejected(RUN_PS2.param, _fin("NEL", 3, run=RUN_PS2), None,
                     "NEL from phaseshift file not consistent with MNEL",
                     tl_root, build_flags, tmp_path)


def _assert_rejected(param_file, fin, param, message, tl_root, build_flags, tmp_path):
    """Run refcalc on fin (and the PARAM text param(), if any, else param_file).

    Assert a non-zero exit code and message on stderr but not on stdout."""
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
    assert message not in res.stdout


@pytest.mark.tier1
def test_refcalc_does_not_check_bravais_subplane_position(tl_root, build_flags, tmp_path):
    """A Bravais layer with a subplane position other than (0,0,0) is not rejected.

    The legacy check for it tests the subplane counter after its loop, where
    it is at least 2, so it never fires; the port keeps the dead check (the
    behaviour is preserved, not fixed). Failing means the check became
    effective."""
    fin = RUN.stdin.read_text().replace("3   0.0000   0.0000   0.0000  Atom N=14",
                                        "3   0.0000   0.0000   0.1000  Atom N=14")
    assert fin != RUN.stdin.read_text()
    (tmp_path / "FIN").write_text(fin)
    exe = build_program(tl_root, CASE.program, RUN.param, build_flags)
    res = run_program(exe, tmp_path, tmp_path / "FIN", timeout=RUN_TIMEOUT_S)
    assert res.returncode == 0, res.stderr
    assert "Bravais layer" not in res.stderr
