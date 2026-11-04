"""Print the PARAM sizes of ref-calc and the runs whose PARAM is larger than needed.

The first table lists every PARAM value with its kind (strict, bound,
explicit; see harness/param_size.py) and how it follows from the input. The
second lists, for every tier-1 refcalc run, the bound values its PARAM sets
above the minimum derived from the FIN. Needs no viperleed.

Usage (env tl-tests, from tests/):
    python -m generate.param_size_table
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

from harness.cases import iter_runs
from harness.param_size import SIZES, oversized, run_sizes


def main():
    print("| value | kind | derived from |")
    print("|---|---|---|")
    for s in SIZES:
        print(f"| `{s.name}` | {s.kind} | {s.source} |")
    print()
    print("| run | value | PARAM | derived |")
    print("|---|---|---|---|")
    n_runs = n_over = 0
    for case, run in iter_runs(("refcalc",)):
        n_runs += 1
        param, derived = run_sizes(run)
        over = oversized(derived, param)
        n_over += bool(over)
        for name, (in_param, minimum) in over.items():
            print(f"| {case.name}/{run.name} | `{name}` | {in_param} | {minimum} |")
    print(f"\n{n_over} of {n_runs} refcalc runs have an oversized PARAM.")


if __name__ == "__main__":
    main()
