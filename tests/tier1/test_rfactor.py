"""Tier 1: rfactor outputs compared with the baseline (tests/cases/rfactor)."""
import pytest

from harness.cases import (compare_outputs, expected_outputs, file_set_failure,
                           format_failures)

from .conftest import case_run_params

PARAMS = case_run_params("rfactor")


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_rfactor_outputs_match_baseline(case, run, run_case):
    workdir = run_case(case, run)
    failures = format_failures(compare_outputs(workdir, expected_outputs(case, run)))
    assert not failures, "\n".join(failures)


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_rfactor_output_file_set_matches_baseline(case, run, run_case):
    workdir = run_case(case, run)
    names = [p.name for p in expected_outputs(case, run)]
    problem = file_set_failure(workdir, run.outputs, names)
    assert not problem, problem
