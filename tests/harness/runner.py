"""Run a program and measure user CPU time and peak resident memory.

Programs are started through a small C launcher (rusage_launcher.c, built on
first use into tests/.build-cache) that reports the program's rusage: on
Linux ru_maxrss survives execve, so a program forked directly from Python
would report at least Python's own peak RSS.
"""
import hashlib
import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

_LAUNCHER_SRC = Path(__file__).with_name("rusage_launcher.c")
_CACHE = Path(__file__).parents[1] / ".build-cache"

_FAILURE_MARKERS = ("Fortran runtime error", "Program received signal",
                    "Error termination")
_NONZERO_STOP = re.compile(r"^\s*STOP\s+0*[1-9]\d*\s*$", re.M)
_THREAD_ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
               "MKL_NUM_THREADS": "1"}


class ProgramFailed(RuntimeError):
    pass


@dataclass
class RunResult:
    returncode: int
    stdout: str
    stderr: str
    user_cpu_s: float
    sys_cpu_s: float
    max_rss_kb: int
    wall_s: float

    def check(self):
        """Raise ProgramFailed unless the program terminated cleanly."""
        text = self.stdout + self.stderr
        if (self.returncode != 0
                or any(m in text for m in _FAILURE_MARKERS)
                or _NONZERO_STOP.search(text)):
            raise ProgramFailed(
                f"exit code {self.returncode}; stderr tail:\n{self.stderr[-2000:]}"
                f"\nstdout tail:\n{self.stdout[-1000:]}")
        return self


def _launcher():
    """Build the launcher once per source and C compiler version; return its path."""
    gcc = subprocess.run(["gcc", "--version"], capture_output=True, text=True,
                         check=True).stdout
    key = hashlib.sha256(_LAUNCHER_SRC.read_bytes() + gcc.encode()).hexdigest()[:16]
    exe = _CACHE / f"rusage_launcher-{key}"
    if not exe.is_file():
        _CACHE.mkdir(parents=True, exist_ok=True)
        tmp = exe.with_suffix(f".tmp{os.getpid()}")
        subprocess.run(["gcc", "-O2", "-o", str(tmp), str(_LAUNCHER_SRC)], check=True)
        tmp.rename(exe)
    return exe


def run_program(exe, cwd, stdin_path, timeout):
    """Run exe in cwd with stdin from stdin_path (or none).

    stdout/stderr are kept in cwd/.stdout and cwd/.stderr. CPU times and
    peak RSS (KiB) are taken from the child's rusage (os.wait4).
    """
    env = dict(os.environ, **_THREAD_ENV)
    cwd = Path(cwd).resolve()  # the launcher runs inside cwd
    out_path, err_path = cwd / ".stdout", cwd / ".stderr"
    usage_path = cwd / ".rusage"
    usage_path.unlink(missing_ok=True)
    launcher = _launcher()
    stdin = open(stdin_path, "rb") if stdin_path else subprocess.DEVNULL
    try:
        with open(out_path, "wb") as out, open(err_path, "wb") as err:
            t0 = time.perf_counter()
            proc = subprocess.Popen(
                [str(launcher), str(usage_path), str(Path(exe).resolve())],
                cwd=cwd, stdin=stdin, stdout=out, stderr=err, env=env)
            try:
                proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()  # the program gets SIGKILL via PR_SET_PDEATHSIG
                proc.wait()
                raise ProgramFailed(f"timeout after {timeout} s: {exe}") from None
            wall = time.perf_counter() - t0
    finally:
        if stdin_path:
            stdin.close()
    if proc.returncode != 0 or not usage_path.is_file():
        raise ProgramFailed(f"rusage launcher failed (exit {proc.returncode}): "
                            f"{err_path.read_text(errors='replace')[-2000:]}")
    status, user, system, maxrss = usage_path.read_text().split()
    return RunResult(os.waitstatus_to_exitcode(int(status)),
                     out_path.read_text(errors="replace"),
                     err_path.read_text(errors="replace"),
                     float(user), float(system), int(maxrss), wall)
