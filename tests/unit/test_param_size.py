"""Unit tests of harness/param_size.py on edited case inputs."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

from pathlib import Path

import pytest

from harness.cases import iter_runs
from harness.param_size import derive_sizes, explicit_values, read_param


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
