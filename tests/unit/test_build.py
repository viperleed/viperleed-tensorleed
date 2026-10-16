"""Tests of the build contract of harness/build.py, on toy and baseline trees.

The contract (tests/README.md, "Build contract") covered here:

- Out-of-tree builds: in Makefile mode make runs on a fresh copy of
  TensErLEED/ with BUILDDIR pointing at a scratch directory, so objects left
  in the tree by an in-tree make are never reused and the tree is not written.
- Cache key: the executable is cached per program and hash; the hash covers
  the sources (every file under TensErLEED/ in Makefile mode, the program's
  sources in legacy mode), the PARAM file, the flags, the make invocation and
  the toolchain (full compiler version, BLAS/LAPACK). Anything that changes
  the executable changes the key.
- Makefile precedence: a TensErLEED/Makefile is used whenever it exists, and
  PARAM is passed to it while a main source still INCLUDEs 'PARAM'; without
  a Makefile the legacy compile needs a PARAM file.
- What a build must not touch: the source tree (a build that changes it
  fails), the cache (a failed build leaves only its log, stale temporary
  directories of dead processes are removed), and concurrent builds of one
  entry run once.

Tests using the baseline tree need the baseline location (see harness/paths.py).
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from harness import build
from harness.build import (BuildError, DEBUG, RELEASE, build_program,
                           needs_param, source_files)

RFACTOR_PARAM = ("      PARAMETER (MNBED = 7, MNBTD = 7)\n"
                 "      PARAMETER (MNET = 86, MNGP = 1422)\n"
                 "      PARAMETER (MNS = 1)\n"
                 "      PARAMETER (MNGAP = 1)\n")

# A Makefile that builds a shell script into BUILDDIR, copying PARAM there
# if it is passed.
TOY_MAKEFILE = """\
$(BUILDDIR)/rfactor: lib/rfacsb.f src/rfactor.f
\tprintf '#!/bin/sh\\necho built FFLAGS=$(FFLAGS)\\n' > $@
\tchmod +x $@
\tif [ -n "$(PARAM)" ]; then cp $(PARAM) $(BUILDDIR)/PARAM.used; fi
rfactor: $(BUILDDIR)/rfactor
"""

# A Makefile that writes objects into the tree and does not track flags.
# Built in the tree, a debug build after a release build only relinks the
# release object and does not trap.
TOY_INTREE_MAKEFILE = """\
refcalc: src/prog.o
\t$(FC) $(FFLAGS) -o $(BUILDDIR)/refcalc src/prog.o $(LDLIBS)
src/prog.o: src/prog.f90
\t$(FC) $(FFLAGS) -c -o src/prog.o src/prog.f90
"""
TOY_OOB_PROGRAM = """\
program refcalc
  implicit none
  integer :: a(3), i
  i = 4
  a = 0
  a(i) = 1          ! out of bounds: caught only with -fcheck=bounds
  print *, 'done', a(1)
end program
"""

# Script run by a child process: build through a given cache directory.
_CHILD = """\
import sys
from pathlib import Path
import harness.build as build
build.CACHE = Path(sys.argv[1])
print(build.build_program(Path(sys.argv[2]), "rfactor", None, build.RELEASE))
"""


@pytest.fixture
def tree(tmp_path, baseline_tl_root):
    """Return a writable copy of the baseline TL tree (legacy mode)."""
    t = tmp_path / "tl"
    shutil.copytree(baseline_tl_root / "TensErLEED", t / "TensErLEED")
    return t


@pytest.fixture
def param(tmp_path):
    """Return the path of a valid rfactor PARAM file."""
    p = tmp_path / "PARAM"
    p.write_text(RFACTOR_PARAM)
    return p


@pytest.fixture
def bad_source(tree, tmp_path):
    """Return (tree, PARAM) where the rfactor library source does not compile."""
    (tree / "TensErLEED/lib/rfacsb.f").write_text("      THIS IS NOT FORTRAN\n")
    bad_param = tmp_path / "PARAM"
    bad_param.write_text("      PARAMETER (MNBED = 1)\n")
    return tree, bad_param


@pytest.fixture
def make_tree(tmp_path):
    """Return a toy TL tree with a Makefile (see TOY_MAKEFILE)."""
    t = tmp_path / "mk"
    (t / "TensErLEED" / "lib").mkdir(parents=True)
    (t / "TensErLEED" / "src").mkdir()
    (t / "TensErLEED" / "lib" / "rfacsb.f").write_text("! module stand-in\n")
    (t / "TensErLEED" / "src" / "rfactor.f").write_text("! no PARAM include\n")
    (t / "TensErLEED" / "Makefile").write_text(TOY_MAKEFILE)
    return t


@pytest.fixture
def add_recipe(make_tree):
    """Return a function appending a shell line to the toy rfactor recipe."""
    def add(line):
        makefile = make_tree / "TensErLEED" / "Makefile"
        makefile.write_text(makefile.read_text().replace(
            "\tchmod +x $@\n", f"\tchmod +x $@\n\t{line}\n"))
    return add


@pytest.fixture
def toy_intree(tmp_path):
    """Return a toy TL tree whose Makefile would build in the tree."""
    t = tmp_path / "toy"
    (t / "TensErLEED" / "src").mkdir(parents=True)
    (t / "TensErLEED" / "src" / "prog.f90").write_text(TOY_OOB_PROGRAM)
    (t / "TensErLEED" / "Makefile").write_text(TOY_INTREE_MAKEFILE)
    return t


def _tree_state(root):
    """Return the sorted (relative path, content) pairs of all files under root."""
    return sorted((str(p.relative_to(root)), p.read_bytes())
                  for p in root.rglob("*") if p.is_file())


def test_baseline_needs_param(baseline_tl_root):
    """The baseline sources INCLUDE 'PARAM', so both programs need a PARAM file."""
    assert needs_param(baseline_tl_root, "refcalc") and needs_param(baseline_tl_root, "rfactor")


def test_source_files_refcalc(baseline_tl_root):
    """Legacy mode hashes exactly the refcalc sources and GLOBAL, nothing else."""
    names = sorted(p.name for p in source_files(baseline_tl_root, "refcalc"))
    assert names == ["GLOBAL", "lib.tleed.f", "ref-calc.f"]


def test_cache_key_changes_with_source(tree, param):
    """An unchanged tree reuses the cached executable; a changed source rebuilds."""
    exe1 = build_program(tree, "rfactor", param, RELEASE)
    assert exe1.is_file()
    assert build_program(tree, "rfactor", param, RELEASE) == exe1        # cached
    src = tree / "TensErLEED/lib/rfacsb.f"
    src.write_text(src.read_text() + "\n")
    assert build_program(tree, "rfactor", param, RELEASE) != exe1        # rebuilt


def test_cache_key_changes_with_param_and_flags(tree, param):
    """Other flags or another PARAM file give another executable."""
    exe1 = build_program(tree, "rfactor", param, RELEASE)
    assert build_program(tree, "rfactor", param, DEBUG) != exe1
    param.write_text(RFACTOR_PARAM.replace("MNS = 1", "MNS = 2"))
    assert build_program(tree, "rfactor", param, RELEASE) != exe1


def test_missing_param_in_legacy_mode_raises(tree):
    """Legacy mode without a PARAM file is a BuildError naming PARAM."""
    with pytest.raises(BuildError, match="PARAM"):
        build_program(tree, "rfactor", None, RELEASE)


def test_compile_error_raises(bad_source):
    """A source that does not compile is reported as 'build failed'."""
    tree_, param_ = bad_source
    with pytest.raises(BuildError, match="build failed"):
        build_program(tree_, "rfactor", param_, RELEASE)


def test_makefile_mode_builds_into_builddir(make_tree):
    """Makefile mode passes the flags and puts the executable into BUILDDIR.

    Without a main source that includes PARAM, no PARAM is passed to make.
    """
    exe = build_program(make_tree, "rfactor", None, RELEASE)
    assert exe.is_file() and "FFLAGS=-O2" in exe.read_text()
    assert not (exe.parent / "PARAM.used").exists()


def test_makefile_mode_hash_covers_all_files_under_tensorleed(make_tree):
    """In Makefile mode any new file under TensErLEED/ changes the cache key."""
    exe1 = build_program(make_tree, "rfactor", None, RELEASE)
    (make_tree / "TensErLEED" / "modules").mkdir()
    (make_tree / "TensErLEED" / "modules" / "kinds.f90").write_text("module kinds\nend\n")
    assert build_program(make_tree, "rfactor", None, RELEASE) != exe1


def test_makefile_wins_and_gets_param_while_sources_include_it(make_tree, param):
    """The Makefile is used even if the main source still INCLUDEs PARAM.

    This is the transition state (modules and Makefile, PARAM not yet gone):
    make must receive the run's PARAM file.
    """
    (make_tree / "TensErLEED" / "src" / "rfactor.f").write_text("      INCLUDE 'PARAM'\n")
    exe = build_program(make_tree, "rfactor", param, RELEASE)
    assert (exe.parent / "PARAM.used").read_text() == param.read_text()


def test_makefile_mode_builds_out_of_tree(toy_intree):
    """Building leaves the tree untouched, and a debug build is not a relinked release.

    With an in-tree make the release object left behind would be reused and
    the debug executable would not trap on the out-of-bounds access.
    """
    before = _tree_state(toy_intree)
    release = build_program(toy_intree, "refcalc", None, RELEASE)
    debug = build_program(toy_intree, "refcalc", None, DEBUG)
    assert _tree_state(toy_intree) == before          # no objects in the tree
    assert subprocess.run([release], capture_output=True).returncode == 0
    proc = subprocess.run([debug], capture_output=True, text=True)
    assert proc.returncode != 0 and "Fortran runtime error" in proc.stderr


def test_debug_build_ignores_release_object_in_tree(toy_intree):
    """A release object already in the tree (newer than its source) is not reused.

    The copy keeps mtimes, so a plain make in the copy would relink the
    object and the debug build would not trap.
    """
    src = toy_intree / "TensErLEED" / "src"
    subprocess.run(["gfortran", "-O2", "-c", "-o", "prog.o", "prog.f90"],
                   cwd=src, check=True)
    debug = build_program(toy_intree, "refcalc", None, DEBUG)
    proc = subprocess.run([debug], capture_output=True, text=True)
    assert proc.returncode != 0 and "Fortran runtime error" in proc.stderr


def test_build_that_modifies_the_tree_fails(make_tree, add_recipe, monkeypatch):
    """A recipe writing into the original tree makes the build fail.

    The target is an absolute path passed through the environment, which make
    imports as variables.
    """
    target = make_tree / "TensErLEED" / "src" / "generated.f90"
    add_recipe('if [ -n "$(TOUCH)" ]; then echo x > $(TOUCH); fi')
    monkeypatch.setenv("TOUCH", str(target))
    with pytest.raises(BuildError, match="changed during the build"):
        build_program(make_tree, "rfactor", None, RELEASE)


def test_cache_key_covers_compiler_and_blas(make_tree, mocker):
    """Another compiler version or BLAS library gives another executable."""
    exe1 = build_program(make_tree, "rfactor", None, RELEASE)
    assert "GNU Fortran" in build.toolchain() and "blas" in build.toolchain()
    mocker.patch.object(build, "toolchain", return_value="GNU Fortran 99|/other/libblas.so|x")
    assert build_program(make_tree, "rfactor", None, RELEASE) != exe1


def test_failed_build_leaves_only_its_log(bad_source, tmp_path, mocker):
    """A failed build leaves one failed-<program>-*.log in the cache, no directory."""
    tree_, param_ = bad_source
    mocker.patch.object(build, "CACHE", tmp_path / "cache")
    with pytest.raises(BuildError, match="build log: .*failed-rfactor-"):
        build_program(tree_, "rfactor", param_, RELEASE)
    left = [p.name for p in (tmp_path / "cache").iterdir()]
    assert len(left) == 1 and left[0].startswith("failed-rfactor-") and left[0].endswith(".log")


def test_stale_tmp_dirs_of_dead_processes_are_removed(make_tree, tmp_path, mocker):
    """A build removes temporary directories of dead processes, not of live ones."""
    cache = tmp_path / "cache"
    mocker.patch.object(build, "CACHE", cache)
    dead = cache / "rfactor-0123456789abcdef.tmp999999999"   # no such pid
    alive = cache / f"rfactor-fedcba9876543210.tmp{os.getpid()}"   # running process
    dead.mkdir(parents=True)
    alive.mkdir()
    build_program(make_tree, "rfactor", None, RELEASE)
    assert not dead.exists() and alive.exists()


def test_concurrent_processes_build_once(make_tree, add_recipe, tmp_path, monkeypatch):
    """Several processes asking for one program on a cold cache run make once.

    This is the situation under pytest -n; the others wait for the lock and
    use the finished executable.
    """
    count = tmp_path / "count"
    add_recipe("echo x >> $(COUNT); sleep 1")
    monkeypatch.setenv("COUNT", str(count))
    tests_dir = Path(__file__).parents[1]
    procs = [subprocess.Popen([sys.executable, "-c", _CHILD,
                               str(tmp_path / "cache"), str(make_tree)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                              env={**os.environ, "PYTHONPATH": str(tests_dir)})
             for _ in range(4)]
    results = [p.communicate(timeout=60) for p in procs]
    assert all(p.returncode == 0 for p in procs), [err for _, err in results]
    exes = {out.strip() for out, _ in results}
    assert len(exes) == 1 and Path(exes.pop()).is_file()
    assert count.read_text() == "x\n"            # built once, the others waited


def test_cache_key_covers_make_invocation(make_tree, mocker):
    """Changing how make is called gives another executable."""
    exe1 = build_program(make_tree, "rfactor", None, RELEASE)
    mocker.patch.object(build, "MAKE", "make -B -j1")
    assert build_program(make_tree, "rfactor", None, RELEASE) != exe1


def test_toolchain_keeps_full_compiler_version(mocker):
    """toolchain() hashes the whole compiler version text, not just its first line."""
    version = "GNU Fortran 15.2\nline 2\n"

    def fake_run(cmd, **kw):
        """Return the current version text for --version, a library path otherwise."""
        out = version if cmd[1:] == ("--version",) else "/usr/lib/x.so"
        return subprocess.CompletedProcess(cmd, 0, stdout=out)
    mocker.patch.object(build.subprocess, "run", side_effect=fake_run)
    mocker.patch.object(build, "toolchain", build.toolchain.__wrapped__)
    key1 = build.toolchain()
    version = "GNU Fortran 15.2\nline 2 changed\n"
    assert build.toolchain() != key1
