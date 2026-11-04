"""Tier 1: delta still builds from the tree as viperleed builds it.

delta is not refactored: it links the legacy library lib/lib.tleed.f, which
viperleed compiles from whatever file matches lib/lib.tleed*. The ref-calc
refactor moves routines out of that file; the routines delta calls must stay
there, unchanged. Checked here: delta compiles and links from the tree's
sources, picked and compiled as viperleed does (with a fixed PARAM,
delta/PARAM), TSCATF, PSTEMP, BESSEL and CPPP in lib/lib.tleed.f are
textually identical to those of the baseline commit, and so is src/GLOBAL
(the constants HARTREE and BOHR, which viperleed copies for delta).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import re
import shutil
import subprocess
from pathlib import Path

import pytest

from generate.baseline_tree import BASELINE_COMMIT
from harness.build import RELEASE

DATA = Path(__file__).with_name("delta")
SHARED_ROUTINES = ("TSCATF", "PSTEMP", "BESSEL", "CPPP")

# a SUBROUTINE or FUNCTION statement (not a comment line), and a bare END
_ROUTINE_START = re.compile(r"^[ \t]+(?:[\w*]+[ \t]+)*(?:SUBROUTINE|FUNCTION)[ \t]+(\w+)",
                            re.I | re.M)
_END = re.compile(r"^[ \t]+END[ \t]*$", re.I | re.M)


def _viperleed_delta_sources(tl):
    """The files viperleed compiles delta from, in its compile order, and GLOBAL.

    As DeltaCompileTask.get_source_files in viperleed 0.14.2
    (calc/sections/deltas.py:125-135): the first match of src/delta*,
    lib/lib.tleed*, lib/lib.delta*. A glob with several matches fails here:
    which one viperleed takes would depend on the directory order."""
    picked = []
    for subdir, pattern in (("src", "delta*"), ("lib", "lib.tleed*"),
                            ("lib", "lib.delta*")):
        found = sorted((tl / subdir).glob(pattern))
        assert len(found) == 1, f"{tl / subdir}/{pattern} matches {found}, need one file"
        picked.append(found[0])
    return picked, tl / "src" / "GLOBAL"


@pytest.mark.tier1
def test_delta_links_with_viperleed_file_set(tl_root, tmp_path):
    """delta compiles and links from the tree as viperleed 0.14.2 builds it.

    Mirrors compile_delta (calc/sections/deltas.py:313-379): PARAM, GLOBAL
    and the three sources copied to one directory, each source compiled
    separately to main.o, lib.tleed.o, lib.delta.o, then linked; the gfortran
    flags of Rparams.getFortranComp (calc/classes/rparams/rparams.py:714-716).
    Failing at the link step (undefined reference) means lib.tleed.f lost a
    routine delta calls."""
    sources, global_file = _viperleed_delta_sources(tl_root / "TensErLEED")
    for path in (*sources, global_file):
        shutil.copy(path, tmp_path / path.name)
    shutil.copy(DATA / "PARAM", tmp_path / "PARAM")
    objects = ("main.o", "lib.tleed.o", "lib.delta.o")
    commands = [f"gfortran {RELEASE.compile} -o {obj} -c {src.name} {RELEASE.link}"
                for src, obj in zip(sources, objects)]
    commands.append(f"gfortran {RELEASE.compile} -o delta {' '.join(objects)} {RELEASE.link}")
    for cmd in commands:
        proc = subprocess.run(cmd.split(), cwd=tmp_path, capture_output=True, text=True)
        assert proc.returncode == 0, f"{cmd}\n{proc.stdout}{proc.stderr[-3000:]}"
    assert (tmp_path / "delta").is_file()


def _routine(text, name):
    """The text of routine name in fixed-form source text, from its
    SUBROUTINE/FUNCTION statement to the first bare END after it."""
    starts = [m for m in _ROUTINE_START.finditer(text) if m[1].upper() == name]
    assert len(starts) == 1, f"{len(starts)} definitions of {name}"
    end = _END.search(text, starts[0].end())
    assert end, f"no END after {name}"
    body = text[starts[0].start():end.end()]
    assert len(_ROUTINE_START.findall(body)) == 1, f"{name} runs into the next routine"
    return body


def _baseline_text(path):
    """The text of path (relative to the repository root) at the baseline commit."""
    return subprocess.run(["git", "-C", str(Path(__file__).parent), "show",
                           f"{BASELINE_COMMIT}:{path}"],
                          check=True, capture_output=True, text=True).stdout


@pytest.mark.tier1
@pytest.mark.parametrize("name", SHARED_ROUTINES)
def test_delta_routine_in_lib_tleed_frozen(tl_root, name):
    """The routine in lib/lib.tleed.f is textually identical to the baseline's.

    delta links these routines from lib.tleed.f; they stay as in the
    baseline commit until delta itself is refactored."""
    path = "TensErLEED/lib/lib.tleed.f"
    current = (tl_root / path).read_text()
    assert _routine(current, name) == _routine(_baseline_text(path), name)


@pytest.mark.tier1
def test_global_frozen(tl_root):
    """src/GLOBAL is identical to the baseline's.

    viperleed copies it for delta, whose numerics use HARTREE and BOHR; no
    ref-calc output depends on them, so only this test guards them."""
    path = "TensErLEED/src/GLOBAL"
    assert (tl_root / path).read_text() == _baseline_text(path)
