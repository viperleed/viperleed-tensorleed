"""Tier 1: rfactor rejects the unsupported input paths with an error.

Each test edits the WEXPEL of the ag100_rfactor/pendry run so that it
selects one unsupported path, runs rfactor and checks that it stops with a
non-zero exit code and a message on stderr that names the parameter.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-22'
__license__ = 'GPLv3+'

import dataclasses

import pytest

from harness.build import build_program
from harness.cases import (RUN_TIMEOUT_S, compare_outputs, expected_outputs,
                           format_failures, iter_runs, run_built)
from harness.runner import run_program

CASE, RUN = next((c, r) for c, r in iter_runs(("rfactor",), only=("ag100_rfactor",))
                 if r.name == "pendry")


def _wexpel_lines():
    return CASE.path.joinpath(RUN.files["WEXPEL"]).read_text().splitlines(keepends=True)


def _replace(old, new):
    """WEXPEL text with the one line equal to old (stripped) replaced by new (or removed)."""
    lines = _wexpel_lines()
    idx = [i for i, line in enumerate(lines) if line.strip() == old]
    assert len(idx) == 1, f"{old!r} not found once in WEXPEL"
    lines[idx[0]] = new + "\n" if new else ""
    return "".join(lines)


def _zero_first_entry(block_line):
    """WEXPEL text with the first entry of the block_line-th integer line after NL1 set to 0.

    The integer lines after NL1 are KAV (0), NBEA (1), IBP (2), MITTEL (3).
    """
    lines = _wexpel_lines()
    first = next(i for i, line in enumerate(lines) if line.strip() == "&END") + 1
    line = lines[first + block_line]
    lines[first + block_line] = "   0" + line[4:]
    return "".join(lines)


REJECTED = [
    pytest.param(lambda: _zero_first_entry(2), "IBP", id="ibp-zero"),
    pytest.param(lambda: _zero_first_entry(1), "NBEA", id="nbea-zero"),
    pytest.param(lambda: _replace("LIMFIL=      1,", " LIMFIL=      2,"), "LIMFIL",
                 id="limfil-two"),
    pytest.param(lambda: _replace("NSSK=    0,", " NSSK=    1,"), "NSSK",
                 id="nssk-nonzero"),
    pytest.param(lambda: _replace("V0RR=      0.0,", " V0RR=      1.0,"), "V0RR",
                 id="v0rr-nonzero"),
    ]


@pytest.fixture(scope="module")
def rfactor_exe(tl_root, build_flags):
    return build_program(tl_root, CASE.program, RUN.param, build_flags)


@pytest.mark.tier1
@pytest.mark.parametrize("wexpel,name", REJECTED)
def test_rfactor_rejects_unsupported_input(wexpel, name, rfactor_exe, tmp_path):
    """rfactor stops with a non-zero exit code and names the parameter on stderr.

    Failing means an unsupported input path is no longer rejected."""
    (tmp_path / "WEXPEL").write_text(wexpel())
    res = run_program(rfactor_exe, tmp_path, RUN.stdin, timeout=RUN_TIMEOUT_S)
    assert res.returncode != 0, f"exit code 0; stderr:\n{res.stderr}"
    assert f"rfactor: {name}" in res.stderr and "not supported" in res.stderr, res.stderr


@pytest.mark.tier1
def test_rfactor_missing_v0rr_counts_as_zero(rfactor_exe, tmp_path):
    """A WEXPEL without V0RR gives the same results as one with V0RR = 0.

    Failing means a missing V0RR is no longer treated as 0."""
    src = tmp_path / "WEXPEL_no_v0rr"
    src.write_text(_replace("V0RR=      0.0,", ""))
    run = dataclasses.replace(RUN, files={"WEXPEL": src})
    run_built(CASE, run, rfactor_exe, tmp_path / "wd")
    failures = format_failures(compare_outputs(tmp_path / "wd", expected_outputs(CASE, RUN)))
    assert not failures, "\n".join(failures)
