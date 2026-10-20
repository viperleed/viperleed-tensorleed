"""Build ref-calc / rfactor from a TensErLEED source tree.

See tests/README.md, section 'Build contract':
- Makefile mode (TensErLEED/Makefile exists): copy TensErLEED/ to a fresh
  directory and run make -B -C <copy> <program> BUILDDIR=<dir> FC=gfortran
  FFLAGS=... LDLIBS=... [PARAM=<file>] there, with PARAM passed while the
  program's main sources still INCLUDE 'PARAM' (during the refactor).
- Legacy mode (no TensErLEED/Makefile): compile the program's sources with
  the case's PARAM the same way viperleed does.
Executables are cached in tests/.build-cache/<program>-<hash>/, where the
hash covers the build mode, the flags, the full compiler version text and the resolved
BLAS/LAPACK libraries (toolchain()), the PARAM file if used, and every
file under TensErLEED/ (Makefile mode) or the program's source files
(legacy mode), including their paths relative to the tree. The build fails
if that hash changes during the build (the build must not modify the tree).
Concurrent processes (pytest -n, parallel suite runs) build each cache entry
once: the build of an entry holds a lock file in tests/.build-cache.locks/,
and processes waiting for it use the finished executable.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import contextlib
import fcntl
import functools
import hashlib
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

CACHE = Path(__file__).parents[1] / ".build-cache"
_PARAM_INCLUDE = re.compile(r"^\s*include\s+['\"]PARAM['\"]", re.I | re.M)


class BuildError(RuntimeError):
    """The build failed, or the source tree changed while it was built."""


@dataclass(frozen=True)
class BuildFlags:
    """Compiler and linker flags of one build variant; part of the cache key."""
    name: str
    compile: str
    link: str = "-llapack -lpthread -lblas -fbacktrace"


RELEASE = BuildFlags("release", "-O2")
DEBUG = BuildFlags("debug", "-O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -Wall")
# how Makefile mode calls make; part of the cache key, so entries built
# with another invocation are not reused
MAKE = "make -B"
# all debug checks except array bounds: for runs whose bounds failures are
# known declaration artefacts (tests/known_debug_failures.toml)
DEBUG_NO_BOUNDS = BuildFlags("debug-nobounds",
                             "-O0 -g -fcheck=all,no-bounds -ffpe-trap=invalid,zero,overflow -Wall")

# program: (sources in compile order [(subdir, file)], copied-only files)
_LEGACY = {
    "refcalc": ([("src", "ref-calc.f"), ("lib", "lib.tleed.f")], [("src", "GLOBAL")]),
    "rfactor": ([("lib", "rfacsb.f"), ("src", "rfactor.f")], []),
}


def uses_makefile(tl_root):
    """True if the tree has a TensErLEED/Makefile (Makefile mode), else legacy mode."""
    return (Path(tl_root) / "TensErLEED" / "Makefile").is_file()


def needs_param(tl_root, program):
    """True if a main source of program (ref-calc.f, lib.tleed.f, resp.
    rfactor.f, rfacsb.f) INCLUDEs 'PARAM'."""
    tl = Path(tl_root) / "TensErLEED"
    comp, _ = _LEGACY[program]
    return any((tl / d / f).is_file()
               and _PARAM_INCLUDE.search((tl / d / f).read_text(errors="replace"))
               for d, f in comp)


def source_files(tl_root, program):
    """All files whose content determines the executable of program."""
    tl = Path(tl_root) / "TensErLEED"
    if uses_makefile(tl_root):
        return sorted(p for p in tl.rglob("*") if p.is_file())
    comp, extra = _LEGACY[program]
    return [tl / d / f for d, f in comp + extra]


@functools.lru_cache(maxsize=None)
def toolchain():
    """Full compiler version text and the resolved BLAS/LAPACK libraries (part of the
    cache key: switching the compiler or the BLAS alternative rebuilds)."""
    def run(*cmd):
        """Run a command and return its stdout, or "missing" if it cannot be started."""
        try:
            return subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
        except OSError:
            return "missing"
    libs = [os.path.realpath(run("gfortran", f"-print-file-name=lib{lib}.so"))
            for lib in ("blas", "lapack")]
    return "|".join([run("gfortran", "--version")] + libs)


def _hash(paths, root, extra=""):
    h = hashlib.sha256(extra.encode())
    for p in paths:
        p = Path(p)
        name = p.relative_to(root) if p.is_relative_to(root) else p.name
        h.update(str(name).encode() + b"\0")
        h.update(p.read_bytes())
    return h.hexdigest()[:16]


def _run(cmd, cwd, log):
    proc = subprocess.run(cmd, cwd=cwd, shell=True, capture_output=True, text=True)
    log.write(f"$ {cmd}\n{proc.stdout}{proc.stderr}\n")
    if proc.returncode != 0:
        raise BuildError(f"build failed: {cmd}\n{proc.stderr[-3000:]}")


def _build_in(tmp, tl_root, program, param, flags, use_make, use_param):
    """Build program into tmp (legacy compile or out-of-tree make)."""
    with open(tmp / "build.log", "w") as log:
        if not use_make:
            if not use_param:
                raise BuildError(f"{tl_root} has neither TensErLEED/Makefile nor PARAM "
                                 f"includes in the {program} sources")
            comp, extra = _LEGACY[program]
            tl = tl_root / "TensErLEED"
            for d, f in comp + extra:
                shutil.copy(tl / d / f, tmp / f)
            shutil.copy(param, tmp / "PARAM")
            objs = []
            for _, f in comp:
                obj = Path(f).stem.replace(".", "_") + ".o"
                _run(f"gfortran {flags.compile} -o {obj} -c {f} {flags.link}", tmp, log)
                objs.append(obj)
            _run(f"gfortran {flags.compile} -o {program} {' '.join(objs)} {flags.link}",
                 tmp, log)
        else:
            # out of tree: make runs on a fresh copy, so the tree is never
            # written to; -B because the copy keeps mtimes, and objects that
            # an in-tree make left there (other flags) would look up to date
            copy = tmp / "tree" / "TensErLEED"
            shutil.copytree(tl_root / "TensErLEED", copy, symlinks=True)
            param_arg = f" PARAM={Path(param).resolve()}" if use_param else ""
            _run(f"{MAKE} -C {copy} {program} BUILDDIR={tmp} "
                 f"FC=gfortran FFLAGS='{flags.compile}' LDLIBS='{flags.link}'{param_arg}",
                 tmp, log)
            shutil.rmtree(tmp / "tree")


def _remove_stale_tmp_dirs():
    """Remove temporary build dirs left by processes that no longer run."""
    for d in CACHE.glob("*.tmp*"):
        pid = d.suffix[len(".tmp"):]
        if pid.isdigit() and not Path(f"/proc/{pid}").exists():
            shutil.rmtree(d, ignore_errors=True)


@contextlib.contextmanager
def _locked(name):
    """Hold an exclusive lock on <CACHE>.locks/<name>.lock (blocks until free).

    The lock files live beside the cache, not in it, and are never removed
    (removing a lock file another process waits on would break the lock)."""
    locks = CACHE.with_name(CACHE.name + ".locks")
    locks.mkdir(parents=True, exist_ok=True)
    with open(locks / f"{name}.lock", "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def build_program(tl_root, program, param, flags):
    """Build program from tl_root (cached); return the executable's path."""
    tl_root = Path(tl_root)
    use_make = uses_makefile(tl_root)
    use_param = needs_param(tl_root, program)
    if use_param and (param is None or not Path(param).is_file()):
        raise BuildError(f"{program} includes PARAM but the run has no PARAM file")
    mode = MAKE if use_make else "legacy"

    def tree_key():
        """Cache key: hash of the flags, the mode, the toolchain and the input files."""
        key_files = source_files(tl_root, program) + ([Path(param)] if use_param else [])
        return _hash(key_files, tl_root / "TensErLEED",
                     f"{program}|{flags}|{mode}|{'param' if use_param else 'noparam'}"
                     f"|{toolchain()}")
    key = tree_key()
    out = CACHE / f"{program}-{key}"
    exe = out / program
    if exe.is_file():
        return exe
    with _locked(f"{program}-{key}"):
        if exe.is_file():  # built by another process while we waited
            return exe
        return _build_entry(tl_root, program, param, flags, use_make, use_param,
                            key, tree_key)


def _build_entry(tl_root, program, param, flags, use_make, use_param, key, tree_key):
    """Build program into CACHE/<program>-<key>/ (caller holds its lock)."""
    out = CACHE / f"{program}-{key}"
    exe = out / program
    _remove_stale_tmp_dirs()
    tmp = CACHE / f"{program}-{key}.tmp{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        _build_in(tmp, tl_root, program, param, flags, use_make, use_param)
        if tree_key() != key:
            raise BuildError(f"the source tree {tl_root / 'TensErLEED'} (or PARAM) "
                             "changed during the build")
        if not (tmp / program).is_file():
            raise BuildError("no executable produced")
    except BaseException as err:
        # keep only the log of a failed build, never a half-built directory
        failed_log = CACHE / f"failed-{program}-{key}.log"
        if (tmp / "build.log").is_file():
            shutil.copy(tmp / "build.log", failed_log)
        shutil.rmtree(tmp, ignore_errors=True)
        if isinstance(err, BuildError):
            raise BuildError(f"{err}\n(build log: {failed_log})") from err
        raise
    tmp.rename(out)
    return exe
