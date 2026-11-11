"""Unit tests of harness/param_size.py on edited case inputs, and of its
LMAX size tables against the DATA tables of ref-calc.f and the routine-check
values."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import json
import re
from pathlib import Path

import pytest

from harness.cases import iter_runs
from harness.param_size import (MNLM_BY_LMAX, MNLMO_BY_LMAX, derive_sizes,
                                explicit_values, read_param)

ROUTINE_VALUES = Path(__file__).parents[1] / "tier1" / "refcalc_routines" / "values.json"


def _fin_and_explicit(lmax):
    """The FIN of ag100/E045_L06 with its LMAX line set to lmax, and the run's explicit sizes."""
    run = next(r for c, r in iter_runs(("refcalc",))
               if c.name == "ag100" and r.name == "E045_L06")
    lines = Path(run.stdin).read_text().splitlines(keepends=True)
    i = next(k for k, line in enumerate(lines) if line.rstrip().endswith("LMAX"))
    lines[i] = f"{lmax:3d}" + lines[i][3:]
    return "".join(lines), explicit_values(read_param(run.param))


@pytest.mark.parametrize("lmax", [0, 1])
def test_derive_sizes_rejects_lmax_below_2(lmax):
    """LMAX 0 and 1 raise: at LMAX 1 ref-calc writes past CAA for any size."""
    with pytest.raises(ValueError, match=f"LMAX = {lmax}; ref-calc needs LMAX >= 2"):
        derive_sizes(*_fin_and_explicit(lmax))


def test_derive_sizes_at_lmax_2():
    """LMAX 2, the smallest accepted, derives the table values of LMAX 2."""
    sizes = derive_sizes(*_fin_and_explicit(2))
    assert (sizes["MLMAX"], sizes["MNLMO"], sizes["MNLM"]) == (2, 70, 76)


def _data_table(source, name):
    """The values of the fixed-form statement DATA name(1),...,name(n)/v1,...,vn/."""
    lines = source.splitlines()
    start = next(i for i, line in enumerate(lines)
                 if re.match(rf"\s+DATA\s+{name}\(1\)", line, re.I))
    statement = lines[start]
    for line in lines[start + 1:]:
        if len(line) < 6 or line[0] in "cC*!" or line[5] in " 0":
            break
        statement += line[6:]
    values = statement.split("/")[1]
    return tuple(int(v) for v in values.split(","))


def _routine_value_counts(array):
    """{LMAX: element count of array} from the routine-check values."""
    cases = json.loads(ROUTINE_VALUES.read_text())["cases"]
    return {int(case["stdin"][1]): case["expected"][array]["n"]
            for case in cases if case["stdin"][0] == "ANGULAR"}


def test_lmax_tables_match_ref_calc_data_and_routine_values(tl_root):
    """param_size's MNLMO/MNLM per LMAX equal ref-calc's DATA NCA/NLMS, which
    cut the CAA/CLM sections passed to the setup routines, and the CAA/CLM
    lengths of the legacy routine-check values (LMAX 2 to 18). At LMAX 1
    (rejected, see test_derive_sizes_rejects_lmax_below_2) NLMS is 0 while
    MNLM is 1, the least array size."""
    source = (tl_root / "TensErLEED" / "src" / "ref-calc.f").read_text()
    nca, nlms = _data_table(source, "NCA"), _data_table(source, "NLMS")
    assert len(nca) == len(nlms) == len(MNLMO_BY_LMAX) == len(MNLM_BY_LMAX) == 18
    assert nca == MNLMO_BY_LMAX
    assert nlms[1:] == MNLM_BY_LMAX[1:]
    assert (nlms[0], MNLM_BY_LMAX[0]) == (0, 1)
    caa, clm = _routine_value_counts("CAA"), _routine_value_counts("CLM")
    assert sorted(caa) == sorted(clm) == list(range(2, 19))
    assert all(caa[lmax] == nca[lmax - 1] for lmax in caa)
    assert all(clm[lmax] == nlms[lmax - 1] for lmax in clm)
