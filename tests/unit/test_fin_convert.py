"""Unit tests of harness/fin_convert.py over the FIN of every tier-1 ref-calc run."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-10'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest

from harness.cases import iter_runs
from harness.fin_convert import convert_fin, converted_values, strip_fin
from harness.param_size import derive_sizes, explicit_values, read_param

RUNS = [pytest.param(run, id=f"{case.name}-{run.name}")
        for case, run in iter_runs(("refcalc",))]
LMAX_MARK = "   LMAX"  # the text edit_fin_energy_lmax splits at


def _inputs(run):
    """(legacy FIN text, PARAM {name: value}) of a run."""
    return Path(run.stdin).read_text(), read_param(run.param)


@pytest.mark.parametrize("run", RUNS)
def test_removing_the_records_gives_back_the_original(run):
    """strip_fin undoes convert_fin byte for byte."""
    fin, param = _inputs(run)
    assert strip_fin(convert_fin(fin, param)) == fin


@pytest.mark.parametrize("run", RUNS)
def test_converted_values_equal_param(run):
    """The three records hold the run's MIDEG, MNL1, MNL2, MNPUN and MNPSI."""
    fin, param = _inputs(run)
    assert converted_values(convert_fin(fin, param)) == explicit_values(param)


@pytest.mark.parametrize("run", RUNS)
def test_records_are_at_the_positions_of_the_format(run):
    """3I5 follows the VPI line, I5 NPUN precedes NPU, I5 NPSI follows IFORM."""
    fin, param = _inputs(run)
    lines = convert_fin(fin, param).splitlines()
    assert "VPI" in lines[6] or "WORKFN" in lines[6]
    assert [int(lines[7][i:i+5]) for i in (0, 5, 10)] == [
        param["MIDEG"], param["MNL1"], param["MNL2"]]
    tst = next(i for i, line in enumerate(lines) if line.rstrip().endswith("TST"))
    assert int(lines[tst + 1][:5]) == param["MNPUN"]
    iform = next(i for i, line in enumerate(lines) if "IFORM" in line)
    assert int(lines[iform + 1][:5]) == param["MNPSI"]
    assert int(lines[iform + 2][:3]) >= 1  # the NEL record, as before


@pytest.mark.parametrize("run", RUNS)
def test_energy_and_lmax_lines_are_where_viperleed_edits_them(run):
    """Line 2 stays at index 1; the first '   LMAX' line is unchanged and none is added.

    edit_fin_energy_lmax (viperleed 0.14.x, calc/sections/refcalc.py:436-455)
    replaces the second line, and splits the rest at the first '   LMAX',
    dropping the text of that line before the mark.
    """
    fin, param = _inputs(run)
    old, new = fin.splitlines(), convert_fin(fin, param).splitlines()
    assert new[:2] == old[:2]
    assert next(line for line in new if LMAX_MARK in line) == (
        next(line for line in old if LMAX_MARK in line))
    assert sum(LMAX_MARK in line for line in new) == sum(LMAX_MARK in line for line in old)


@pytest.mark.parametrize("run", RUNS)
def test_sizes_from_converted_fin_equal_sizes_from_param(run):
    """derive_sizes(explicit=None) on the converted FIN gives the PARAM-based sizes."""
    fin, param = _inputs(run)
    assert derive_sizes(convert_fin(fin, param), None) == (
        derive_sizes(fin, explicit_values(param)))


def test_line_endings_of_the_new_records_follow_the_file():
    """A FIN with CRLF endings gets CRLF endings on the added records."""
    run = RUNS[0].values[0]
    fin, param = _inputs(run)
    crlf = fin.replace("\n", "\r\n")
    converted = convert_fin(crlf, param)
    assert "\n" not in converted.replace("\r\n", "")
    assert strip_fin(converted) == crlf
