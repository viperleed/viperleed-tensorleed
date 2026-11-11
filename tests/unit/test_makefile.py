"""Tests of the tree's own TensErLEED/Makefile (not of the harness build).

Free-form sources must conform to Fortran 2018: the Makefile compiles every
.f90 source it lists for a program (print-<PROGRAM>_SOURCES) with -std=f2018,
and they compile cleanly with it.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-04'
__license__ = 'GPLv3+'

import shlex
import shutil
import subprocess

from harness.build import program_sources
from harness.cases import iter_runs

PROGRAMS = ("refcalc", "rfactor")


def test_free_form_sources_compile_with_std_f2018(tl_root, tmp_path):
    """A build of both programs compiles each free-form source with -std=f2018
    and succeeds, so a non-standard construct in any of them fails it."""
    tree = tmp_path / "TensErLEED"
    shutil.copytree(tl_root / "TensErLEED", tree, symlinks=True)
    free_form = [str(p.relative_to(tl_root / "TensErLEED"))
                 for prog in PROGRAMS for p in program_sources(tl_root, prog)
                 if p.suffix == ".f90"]
    assert free_form
    param = next(run.param for case, run in iter_runs() if case.program == "refcalc")
    proc = subprocess.run(["make", "-B", "-C", str(tree), *PROGRAMS,
                           f"BUILDDIR={tmp_path / 'build'}", "FFLAGS=-O0",
                           f"PARAM={param}"],
                          capture_output=True, text=True)
    errors = [line for line in proc.stderr.splitlines() if "Error" in line]
    assert proc.returncode == 0, "\n".join(errors) or proc.stderr[-3000:]
    compile_lines = [shlex.split(line) for line in proc.stdout.splitlines()
                     if " -c " in line]
    for src in free_form:
        lines = [words for words in compile_lines if words[-1] == src]
        assert len(lines) == 1, f"{src} is not compiled exactly once"
        assert "-std=f2018" in lines[0], f"{src}: {shlex.join(lines[0])}"


def test_incremental_make_rebuilds_objects_using_a_changed_module(tl_root, tmp_path):
    """After a module of ref-calc changes, make without -B rebuilds the
    driver, which uses it, and not lib.tleed.f, which does not."""
    tree = tmp_path / "TensErLEED"
    shutil.copytree(tl_root / "TensErLEED", tree, symlinks=True)
    param = next(run.param for case, run in iter_runs() if case.program == "refcalc")
    make = ["make", "-C", str(tree), "refcalc", f"BUILDDIR={tmp_path / 'build'}",
            "FFLAGS=-O0", f"PARAM={param}"]
    first = subprocess.run(make, capture_output=True, text=True)
    assert first.returncode == 0, first.stderr[-3000:]
    kinds = tree / "lib" / "kinds.f90"
    kinds.write_text(kinds.read_text().replace(
        "end module", "integer, parameter :: new_kind = 1\nend module", 1))
    second = subprocess.run(make, capture_output=True, text=True)
    assert second.returncode == 0, second.stderr[-3000:]
    compiled = [shlex.split(line)[-1] for line in second.stdout.splitlines()
                if " -c " in line]
    assert "src/ref-calc.f" in compiled
    assert "lib/lib.tleed.f" not in compiled
