"""Tier 1: rfactor outputs compared with the baseline (tests/cases/rfactor)."""
import pytest

from harness.cases import (compare_outputs, expected_outputs, file_set_failure,
                           format_failures)

from .conftest import case_run_params

PARAMS = case_run_params("rfactor")


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_rfactor_outputs_match_baseline(case, run, run_case):
    """Every expected output of an rfactor run matches the candidate output within tolerances.

    Failing means the numerical results (or their print format) changed."""
    workdir = run_case(case, run)
    failures = format_failures(compare_outputs(workdir, expected_outputs(case, run)))
    assert not failures, "\n".join(failures)


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_rfactor_output_file_set_matches_baseline(case, run, run_case):
    """An rfactor run produces exactly the expected set of output files.

    Failing means a file went missing or a new one appeared."""
    workdir = run_case(case, run)
    names = [p.name for p in expected_outputs(case, run)]
    problem = file_set_failure(workdir, run.outputs, names)
    assert not problem, problem
