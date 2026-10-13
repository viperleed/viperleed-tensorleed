"""Shared helpers for tier-1 tests: run each (case, run) once per session."""
import pytest

from harness.cases import execute_run, load_cases


def case_run_params(kind):
    # xdist_group: with pytest -n, all tests of one (case, run) go to the
    # same worker, so run_case executes the run once (see tests/conftest.py)
    return [pytest.param(case, run, id=f"{case.name}/{run.name}",
                         marks=pytest.mark.xdist_group(f"{case.name}/{run.name}"))
            for case in load_cases(kind) for run in case.runs]


@pytest.fixture(scope="session")
def run_case(tl_root, build_flags, tmp_path_factory):
    """run_case(case, run) -> work dir with the outputs (cached per session)."""
    done = {}

    def _run(case, run):
        key = (case.path, run.name)
        if key not in done:
            workdir = tmp_path_factory.mktemp(f"{case.name}-{run.name}")
            try:
                execute_run(case, run, tl_root, build_flags, workdir)
            except Exception as err:  # cache failures too, re-raise every time
                done[key] = err
            else:
                done[key] = workdir
        if isinstance(done[key], Exception):
            raise done[key]
        return done[key]
    return _run
