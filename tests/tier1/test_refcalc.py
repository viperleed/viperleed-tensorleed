"""Tier 1: refcalc outputs compared with the baseline (tests/cases/refcalc)."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import pytest

from harness.cases import (compare_outputs, expected_outputs, file_set_failure,
                           format_failures)

from .conftest import case_run_params

PARAMS = case_run_params("refcalc")


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_refcalc_outputs_match_baseline(case, run, run_case):
    """Every expected output of a refcalc run matches the candidate output within tolerances.

    Failing means the numerical results (or their print format) changed."""
    workdir = run_case(case, run)
    failures = format_failures(compare_outputs(workdir, expected_outputs(case, run)))
    assert not failures, "\n".join(failures)


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_refcalc_output_file_set_matches_baseline(case, run, run_case):
    """A refcalc run produces exactly the expected set of output files.

    Failing means a file went missing or a new one appeared."""
    workdir = run_case(case, run)
    names = [p.name for p in expected_outputs(case, run)]
    problem = file_set_failure(workdir, run.outputs, names)
    assert not problem, problem
