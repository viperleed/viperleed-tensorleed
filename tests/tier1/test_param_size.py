"""Tier 1: the PARAM of every refcalc run follows from its input (harness/param_size.py)."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import re
import shutil
from pathlib import Path

import pytest

from harness.build import build_program, needs_param
from harness.cases import iter_runs, run_built
from harness.fixed_form import routine
from harness.param_size import (BOUND, EXPLICIT, KIND, STRICT, param_text, read_param,
                                run_sizes, size_mismatches)

from .conftest import case_run_params

PARAMS = case_run_params("refcalc")

# Oversized PARAM variants, oversized_param/<case>/<run>/PARAM: the run's
# PARAM with bound sizes above what the input needs (header of each file).
OVERSIZED = Path(__file__).with_name("oversized_param")
OVERSIZED_PARAMS = [
    pytest.param(case, run, OVERSIZED / case.name / run.name / "PARAM",
                 id=f"{case.name}/{run.name}",
                 marks=pytest.mark.xdist_group(f"{case.name}/{run.name}"))
    for case, run in iter_runs(("refcalc",))
    if (OVERSIZED / case.name / run.name / "PARAM").is_file()
    ]


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", PARAMS)
def test_refcalc_param_sizes_follow_from_input(case, run):
    """Each strict PARAM size equals, and each bound size is at least, the value derived from the FIN.

    Failing means the derivation in harness/param_size.py misreads the
    input, or the case's PARAM does not fit its FIN."""
    param, derived = run_sizes(run)
    problems = size_mismatches(derived, param)
    assert not problems, "\n".join(problems)


# IF (a .ne. b) THEN, and its block up to END IF or ELSE
_CHECK = re.compile(r"^[ \t]+IF[ \t]*\(([^\n]+?)\.ne\.([^\n]+?)\)[ \t]*THEN(.*?)"
                    r"^[ \t]+(?:END[ \t]*IF|ELSE)", re.I | re.M | re.S)
# NAME = MNAME in ref-calc.f: a local copy of a PARAM value
_COPY = re.compile(r"^[ \t]+(\w+)[ \t]*=[ \t]*(M\w+)[ \t]*$", re.M)
_NAME = re.compile(r"[A-Za-z]\w*")


def _param_definitions(text):
    """{name: expression} of the PARAMETER statements in PARAM text."""
    defs = {}
    for body in re.findall(r"PARAMETER\s*\((.*)\)", text, re.I):
        for item in body.split(","):
            name, _, expr = item.partition("=")
            defs[name.strip().upper()] = expr.upper()
    return defs


def _equality_checks(tl_root):
    """The sizes (SIZES names) compared in an IF (a.ne.b) THEN ... STOP of
    READIN or READGEO."""
    tl = tl_root / "TensErLEED"
    lib = (tl / "lib" / "lib.tleed.f").read_text()
    copies = {local.upper(): param.upper()
              for local, param in _COPY.findall((tl / "src" / "ref-calc.f").read_text())}
    defs = _param_definitions(Path(next(iter_runs(("refcalc",)))[1].param).read_text())

    def sizes(name):
        """The SIZES names a PARAM name is defined from."""
        if name in KIND:
            return {name}
        return {s for n in _NAME.findall(defs.get(name, "")) for s in sizes(n.upper())}

    found = set()
    for name in ("READIN", "READGEO"):
        for lhs, rhs, block in _CHECK.findall(routine(lib, name)):
            if not re.search(r"^[ \t]+STOP\b", block, re.I | re.M):
                continue
            for word in _NAME.findall(lhs + " " + rhs):
                if word.upper() in copies:
                    found |= sizes(copies[word.upper()])
    return found


@pytest.mark.tier1
def test_strict_sizes_are_the_ones_ref_calc_checks(tl_root):
    """The strict sizes are exactly those READIN/READGEO compare with .ne. and stop on.

    The checks are read from lib.tleed.f (IF (a.ne.b) THEN ... STOP), their
    variables traced to PARAM through ref-calc.f's copies (NSITE = MNSITE)
    and the PARAMETER definitions (MLMAX1 = MLMAX+1). Explicit sizes in the
    same checks (MNL1*MNL2 against MKNBS) do not count, and no bound size
    is among them. Failing means SIZES classifies a size differently from
    the source."""
    explicit = {name for name, kind in KIND.items() if kind == EXPLICIT}
    assert _equality_checks(tl_root) - explicit == {
        name for name, kind in KIND.items() if kind == STRICT}


def _output_files(workdir, patterns):
    return {p.name: p for pattern in patterns
            for p in Path(workdir).glob(pattern) if p.is_file()}


def _assert_same_outputs(workdir, reference_dir, patterns):
    """Assert workdir has the same output files as reference_dir, byte for byte."""
    reference = _output_files(reference_dir, patterns)
    actual = _output_files(workdir, patterns)
    assert sorted(actual) == sorted(reference)
    differ = [name for name in sorted(reference)
              if actual[name].read_bytes() != reference[name].read_bytes()]
    assert not differ, f"outputs differ from the case PARAM build: {differ}"


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
    _assert_same_outputs(workdir, run_case(case, run), run.outputs)


@pytest.mark.tier1
def test_oversized_param_variants_exist():
    """Both oversized PARAM variants under oversized_param/ belong to a refcalc run.

    Failing means a variant was removed or its case or run renamed, so the
    oversized-PARAM tests below would silently run fewer builds."""
    found = sorted(p.relative_to(OVERSIZED).as_posix()
                   for p in OVERSIZED.glob("*/*/PARAM"))
    assert found == ["ag100/E045_L06/PARAM", "ir100_2x1_o/E050_L06/PARAM"]
    assert len(OVERSIZED_PARAMS) == len(found)


@pytest.mark.tier1
@pytest.mark.parametrize("case,run,variant", OVERSIZED_PARAMS)
def test_oversized_param_exceeds_case_param_in_bounds_only(case, run, variant):
    """The variant differs from the run's PARAM only in bound sizes, each set larger.

    Failing means the variant no longer oversizes anything (the comparison
    of outputs would be vacuous) or changes a strict or explicit size."""
    case_param = read_param(run.param)
    param = read_param(variant)
    changed = {name for name in KIND if param[name] != case_param[name]}
    assert changed, "the variant has the same sizes as the run's PARAM"
    assert all(KIND[name] == BOUND for name in changed), sorted(changed)
    assert all(param[name] > case_param[name] for name in changed), sorted(changed)


@pytest.mark.tier1
@pytest.mark.parametrize("case,run,variant", OVERSIZED_PARAMS)
def test_oversized_param_gives_identical_outputs(case, run, variant, run_case, tl_root,
                                                 build_flags, tmp_path):
    """refcalc built with an oversized PARAM writes byte-identical outputs to the case PARAM build.

    viperleed writes bound sizes (MNSUB of bulk layers) larger than the input
    needs. The reference is the same tree built with the run's PARAM
    (run_case); the variant must be a separate build (cache entry).
    Failing means an output depends on an array size beyond what the input
    needs."""
    if not needs_param(tl_root, case.program):
        pytest.skip("refcalc of this tree takes no PARAM")
    reference_dir = run_case(case, run)
    # under the name of the run's PARAM, so that only the content can make
    # the cache key differ
    param = tmp_path / Path(run.param).name
    shutil.copy(variant, param)
    exe = build_program(tl_root, case.program, param, build_flags)
    case_exe = build_program(tl_root, case.program, run.param, build_flags)
    assert exe.parent != case_exe.parent, "the variant was not built on its own"
    workdir = tmp_path / "run"
    run_built(case, run, exe, workdir)
    _assert_same_outputs(workdir, reference_dir, run.outputs)
