from pathlib import Path
import shutil, subprocess, pytest
from harness.build import needs_param, source_files, build_program, RELEASE, DEBUG, BuildError

RFACTOR_PARAM = ("      PARAMETER (MNBED = 7, MNBTD = 7)\n      PARAMETER (MNET = 86, MNGP = 1422)\n"
                 "      PARAMETER (MNS = 1)\n      PARAMETER (MNGAP = 1)\n")

@pytest.fixture
def tree(tmp_path, baseline_tl_root):
    t = tmp_path / "tl"
    shutil.copytree(baseline_tl_root / "TensErLEED", t / "TensErLEED")
    return t

@pytest.fixture
def param(tmp_path):
    p = tmp_path / "PARAM"
    p.write_text(RFACTOR_PARAM)
    return p

def test_baseline_needs_param(baseline_tl_root):
    assert needs_param(baseline_tl_root, "refcalc") and needs_param(baseline_tl_root, "rfactor")

def test_source_files_refcalc(baseline_tl_root):
    names = sorted(p.name for p in source_files(baseline_tl_root, "refcalc"))
    assert names == ["GLOBAL", "lib.tleed.f", "ref-calc.f"]

def test_cache_key_changes_with_source(tree, param):
    exe1 = build_program(tree, "rfactor", param, RELEASE)
    assert exe1.is_file()
    assert build_program(tree, "rfactor", param, RELEASE) == exe1        # cached
    src = tree / "TensErLEED/lib/rfacsb.f"; src.write_text(src.read_text() + "\n")
    assert build_program(tree, "rfactor", param, RELEASE) != exe1        # rebuilt

def test_cache_key_changes_with_param_and_flags(tree, param):
    exe1 = build_program(tree, "rfactor", param, RELEASE)
    assert build_program(tree, "rfactor", param, DEBUG) != exe1
    param.write_text(RFACTOR_PARAM.replace("MNS = 1", "MNS = 2"))
    assert build_program(tree, "rfactor", param, RELEASE) != exe1

def test_missing_param_in_legacy_mode_raises(tree):
    with pytest.raises(BuildError, match="PARAM"):
        build_program(tree, "rfactor", None, RELEASE)

def test_compile_error_raises(tree, tmp_path):
    (tree / "TensErLEED/lib/rfacsb.f").write_text("      THIS IS NOT FORTRAN\n")
    param = tmp_path / "PARAM"; param.write_text("      PARAMETER (MNBED = 1)\n")
    with pytest.raises(BuildError, match="build failed"):
        build_program(tree, "rfactor", param, RELEASE)


TOY_MAKEFILE = """\
$(BUILDDIR)/rfactor: lib/rfacsb.f src/rfactor.f
\tprintf '#!/bin/sh\\necho built FFLAGS=$(FFLAGS)\\n' > $@
\tchmod +x $@
\tif [ -n "$(PARAM)" ]; then cp $(PARAM) $(BUILDDIR)/PARAM.used; fi
rfactor: $(BUILDDIR)/rfactor
"""


@pytest.fixture
def make_tree(tmp_path):
    t = tmp_path / "mk"
    (t / "TensErLEED" / "lib").mkdir(parents=True)
    (t / "TensErLEED" / "src").mkdir()
    (t / "TensErLEED" / "lib" / "rfacsb.f").write_text("! module stand-in\n")
    (t / "TensErLEED" / "src" / "rfactor.f").write_text("! no PARAM include\n")
    (t / "TensErLEED" / "Makefile").write_text(TOY_MAKEFILE)
    return t


def test_makefile_mode_builds_into_builddir(make_tree):
    exe = build_program(make_tree, "rfactor", None, RELEASE)
    assert exe.is_file() and "FFLAGS=-O2" in exe.read_text()
    assert not (exe.parent / "PARAM.used").exists()


def test_makefile_mode_hash_covers_all_files_under_tensorleed(make_tree):
    exe1 = build_program(make_tree, "rfactor", None, RELEASE)
    (make_tree / "TensErLEED" / "modules").mkdir()
    (make_tree / "TensErLEED" / "modules" / "kinds.f90").write_text("module kinds\nend\n")
    assert build_program(make_tree, "rfactor", None, RELEASE) != exe1


def test_makefile_wins_and_gets_param_while_sources_include_it(make_tree, param):
    # transition: modules + Makefile, main source still includes PARAM
    (make_tree / "TensErLEED" / "src" / "rfactor.f").write_text("      INCLUDE 'PARAM'\n")
    exe = build_program(make_tree, "rfactor", param, RELEASE)
    assert (exe.parent / "PARAM.used").read_text() == param.read_text()


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


@pytest.fixture
def toy_intree(tmp_path):
    t = tmp_path / "toy"
    (t / "TensErLEED" / "src").mkdir(parents=True)
    (t / "TensErLEED" / "src" / "prog.f90").write_text(TOY_OOB_PROGRAM)
    (t / "TensErLEED" / "Makefile").write_text(TOY_INTREE_MAKEFILE)
    return t


def _tree_state(root):
    return sorted((str(p.relative_to(root)), p.read_bytes())
                  for p in root.rglob("*") if p.is_file())


def test_makefile_mode_builds_out_of_tree(toy_intree):
    before = _tree_state(toy_intree)
    release = build_program(toy_intree, "refcalc", None, RELEASE)
    debug = build_program(toy_intree, "refcalc", None, DEBUG)
    assert _tree_state(toy_intree) == before          # no objects in the tree
    assert subprocess.run([release], capture_output=True).returncode == 0
    proc = subprocess.run([debug], capture_output=True, text=True)
    assert proc.returncode != 0 and "Fortran runtime error" in proc.stderr


# A release object left in the tree by an in-tree make, newer than its
# source. The copy keeps mtimes, so a plain make in the copy would relink it
# and the debug build would not trap.
def test_debug_build_ignores_release_object_in_tree(toy_intree):
    src = toy_intree / "TensErLEED" / "src"
    subprocess.run(["gfortran", "-O2", "-c", "-o", "prog.o", "prog.f90"],
                   cwd=src, check=True)
    debug = build_program(toy_intree, "refcalc", None, DEBUG)
    proc = subprocess.run([debug], capture_output=True, text=True)
    assert proc.returncode != 0 and "Fortran runtime error" in proc.stderr


def test_build_that_modifies_the_tree_fails(make_tree, monkeypatch):
    # a recipe writing into the original tree (absolute path, here passed
    # through the environment, which make imports as variables)
    target = make_tree / "TensErLEED" / "src" / "generated.f90"
    (make_tree / "TensErLEED" / "Makefile").write_text(
        TOY_MAKEFILE.replace("\tchmod +x $@\n",
                             "\tchmod +x $@\n\tif [ -n \"$(TOUCH)\" ]; then echo x > $(TOUCH); fi\n"))
    monkeypatch.setenv("TOUCH", str(target))
    with pytest.raises(BuildError, match="changed during the build"):
        build_program(make_tree, "rfactor", None, RELEASE)


def test_cache_key_covers_compiler_and_blas(make_tree, monkeypatch):
    import harness.build as build
    exe1 = build_program(make_tree, "rfactor", None, RELEASE)
    assert "GNU Fortran" in build.toolchain() and "blas" in build.toolchain()
    monkeypatch.setattr(build, "toolchain", lambda: "GNU Fortran 99|/other/libblas.so|x")
    assert build_program(make_tree, "rfactor", None, RELEASE) != exe1


def test_failed_build_leaves_only_its_log(tree, tmp_path, monkeypatch):
    import harness.build as build
    monkeypatch.setattr(build, "CACHE", tmp_path / "cache")
    (tree / "TensErLEED/lib/rfacsb.f").write_text("      THIS IS NOT FORTRAN\n")
    param = tmp_path / "PARAM"; param.write_text("      PARAMETER (MNBED = 1)\n")
    with pytest.raises(BuildError, match="build log: .*failed-rfactor-"):
        build_program(tree, "rfactor", param, RELEASE)
    left = [p.name for p in (tmp_path / "cache").iterdir()]
    assert len(left) == 1 and left[0].startswith("failed-rfactor-") and left[0].endswith(".log")


def test_stale_tmp_dirs_of_dead_processes_are_removed(make_tree, tmp_path, monkeypatch):
    import os
    import harness.build as build
    cache = tmp_path / "cache"
    monkeypatch.setattr(build, "CACHE", cache)
    dead = cache / "rfactor-0123456789abcdef.tmp999999999"   # no such pid
    alive = cache / f"rfactor-fedcba9876543210.tmp{os.getpid()}"   # running process
    dead.mkdir(parents=True); alive.mkdir()
    build_program(make_tree, "rfactor", None, RELEASE)
    assert not dead.exists() and alive.exists()


_CHILD = """\
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import harness.build as build
build.CACHE = Path(sys.argv[2])
print(build.build_program(Path(sys.argv[3]), "rfactor", None, build.RELEASE))
"""


def test_concurrent_processes_build_once(make_tree, tmp_path, monkeypatch):
    # pytest -n: several workers ask for the same program on a cold cache
    import sys
    count = tmp_path / "count"
    (make_tree / "TensErLEED" / "Makefile").write_text(
        TOY_MAKEFILE.replace("\tchmod +x $@\n",
                             "\tchmod +x $@\n\techo x >> $(COUNT); sleep 1\n"))
    monkeypatch.setenv("COUNT", str(count))
    tests_dir = Path(__file__).parents[1]
    procs = [subprocess.Popen([sys.executable, "-c", _CHILD, str(tests_dir),
                               str(tmp_path / "cache"), str(make_tree)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
             for _ in range(4)]
    results = [p.communicate(timeout=60) for p in procs]
    assert all(p.returncode == 0 for p in procs), [err for _, err in results]
    exes = {out.strip() for out, _ in results}
    assert len(exes) == 1 and Path(exes.pop()).is_file()
    assert count.read_text() == "x\n"            # built once, the others waited


def test_cache_key_covers_make_invocation(make_tree, monkeypatch):
    import harness.build as build
    exe1 = build_program(make_tree, "rfactor", None, RELEASE)
    monkeypatch.setattr(build, "MAKE", "make -B -j1")
    assert build_program(make_tree, "rfactor", None, RELEASE) != exe1


def test_toolchain_keeps_full_compiler_version(monkeypatch):
    import harness.build as build
    version = "GNU Fortran 15.2\nline 2\n"

    def fake_run(cmd, **kw):
        out = version if cmd[1:] == ("--version",) else "/usr/lib/x.so"
        return subprocess.CompletedProcess(cmd, 0, stdout=out)
    monkeypatch.setattr(build.subprocess, "run", fake_run)
    monkeypatch.setattr(build, "toolchain", build.toolchain.__wrapped__)
    key1 = build.toolchain()
    version = "GNU Fortran 15.2\nline 2 changed\n"
    assert build.toolchain() != key1
