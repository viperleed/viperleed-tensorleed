"""Tier 1: rfactor outputs compared with the baseline (tests/cases/rfactor)."""
import pytest

from harness.cases import expected_outputs, tolerance_for
from harness.compare import compare_files

from .conftest import case_run_params

PARAMS = case_run_params("rfactor")


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_rfactor_outputs_match_baseline(case, run, run_case):
    workdir = run_case(case, run)
    failures = []
    for exp in expected_outputs(case, run):
        tol, ignore, skip = tolerance_for(exp.name)
        res = compare_files(workdir / exp.name, exp, tol, ignore, skip)
        if not res.ok:
            failures.append(f"{exp.name}: " + "; ".join(res.first_failures))
    assert not failures, "\n".join(failures)


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_rfactor_output_file_set_matches_baseline(case, run, run_case):
    workdir = run_case(case, run)
    produced = {p.name for pattern in run.outputs for p in workdir.glob(pattern) if p.is_file()}
    expected = {p.name for p in expected_outputs(case, run)}
    assert produced == expected, (f"missing: {sorted(expected - produced)}, "
                                  f"extra: {sorted(produced - expected)}")
