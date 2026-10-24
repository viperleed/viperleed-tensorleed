"""Tier 1: runs with poison_heap = true, repeated with a poisoned heap.

glibc fills new memory with a non-zero pattern and serves every allocation
from the heap (harness.runner.POISONED_HEAP_ENV). A read of an array element
the program never wrote then changes the outputs, whereas in a normal run
it mostly reads fresh, zeroed pages (README, 'Rules').
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-24'
__license__ = 'GPLv3+'

import pytest

from harness.cases import compare_outputs, execute_run, expected_outputs, format_failures, iter_runs
from harness.runner import POISONED_HEAP_ENV


def _params():
    for case, run in iter_runs():
        if run.poison_heap:
            yield pytest.param(case, run, id=f"{case.name}/{run.name}")


@pytest.mark.tier1
@pytest.mark.parametrize("case,run", list(_params()))
def test_outputs_match_baseline_with_poisoned_heap(case, run, tl_root, build_flags, tmp_path):
    """Every expected output matches the baseline also when new memory is not zero.

    Failing means the program reads memory it never wrote."""
    execute_run(case, run, tl_root, build_flags, tmp_path, env=POISONED_HEAP_ENV)
    failures = format_failures(compare_outputs(tmp_path, expected_outputs(case, run)))
    assert not failures, "\n".join(failures)
