"""The kind parameters of ref-calc and rfactor are the same.

ref-calc takes them from lib/kinds.f90, rfactor from its own copy of module
kinds in lib/rfacsb.f90 (viperleed compiles rfactor from rfacsb* and
rfactor* only). Both copies must define the module identically.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import re

_KINDS_MODULE = re.compile(r"^module kinds\n.*?^end module kinds$", re.M | re.S)


def _kinds_module(path):
    """The text of module kinds (from 'module kinds' to 'end module kinds') in path."""
    found = _KINDS_MODULE.findall(path.read_text())
    assert len(found) == 1, f"{path}: {len(found)} definitions of module kinds"
    return found[0]


def test_kinds_modules_identical(tl_root):
    """module kinds in lib/kinds.f90 and lib/rfacsb.f90 are identical."""
    lib = tl_root / "TensErLEED" / "lib"
    assert _kinds_module(lib / "kinds.f90") == _kinds_module(lib / "rfacsb.f90")
