"""Tier 1: the PARAM of every refcalc run follows from its input (harness/param_size.py)."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest

from harness.build import build_program, needs_param
from harness.cases import run_built
from harness.param_size import param_text, run_sizes, size_mismatches

from .conftest import case_run_params

PARAMS = case_run_params("refcalc")


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_refcalc_param_sizes_follow_from_input(case, run):
    """Each strict PARAM size equals, and each bound size is at least, the value derived from the FIN.

    Failing means the derivation in harness/param_size.py misreads the
    input, or the case's PARAM does not fit its FIN."""
    param, derived = run_sizes(run)
    problems = size_mismatches(derived, param)
    assert not problems, "\n".join(problems)


def _output_files(workdir, patterns):
    return {p.name: p for pattern in patterns
            for p in Path(workdir).glob(pattern) if p.is_file()}


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_refcalc_minimal_param_gives_identical_outputs(case, run, run_case, tl_root,
                                                       build_flags, tmp_path):
    """refcalc built with the derived minimal PARAM writes byte-identical outputs to the case PARAM build.

    The reference is the same tree built with the case's PARAM (run_case),
    so only the sizes differ. If the minimal PARAM is byte-identical to the
    case's, both builds are the same cache entry and nothing is run.
    Failing means an output depends on an array size beyond what the input
    needs."""
    if not needs_param(tl_root, case.program):
        pytest.skip("refcalc of this tree takes no PARAM")
    param, derived = run_sizes(run)
    minimal = param_text(derived)
    if minimal == Path(run.param).read_text():
        return
    minimal_param = tmp_path / Path(run.param).name
    minimal_param.write_text(minimal)
    exe = build_program(tl_root, case.program, minimal_param, build_flags)
    workdir = tmp_path / "run"
    run_built(case, run, exe, workdir)
    reference = _output_files(run_case(case, run), run.outputs)
    actual = _output_files(workdir, run.outputs)
    assert sorted(actual) == sorted(reference)
    differ = [name for name in sorted(reference)
              if actual[name].read_bytes() != reference[name].read_bytes()]
    assert not differ, f"outputs differ from the case PARAM build: {differ}"
