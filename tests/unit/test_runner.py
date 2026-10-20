"""Tests of harness.runner: running programs, failure detection and measurement."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import resource
import subprocess
from pathlib import Path

import pytest

from harness import runner
from harness.runner import ProgramFailed, run_program


def test_success_and_stdin(tmp_path, stub_script):
    """A program gets stdin, and its output, peak memory and CPU time are recorded."""
    exe = stub_script(tmp_path / "prog.py", "import sys; print(sys.stdin.read().upper())")
    (tmp_path / "in").write_text("abc")
    r = run_program(exe, tmp_path, tmp_path / "in", timeout=30)
    r.check()
    assert "ABC" in r.stdout and r.max_rss_kb > 0 and r.user_cpu_s >= 0


def test_nonzero_exit_fails(tmp_path, stub_script):
    """A non-zero exit code fails the check."""
    exe = stub_script(tmp_path / "prog.py", "import sys; sys.exit(3)")
    with pytest.raises(ProgramFailed):
        run_program(exe, tmp_path, None, timeout=30).check()


def test_fortran_runtime_error_text_fails(tmp_path, stub_script):
    """A "Fortran runtime error" message fails the check even with exit code 0."""
    exe = stub_script(tmp_path / "prog.py", "import sys; "
                      "sys.stderr.write('Fortran runtime error: Index 3 out of bounds')")
    with pytest.raises(ProgramFailed):
        run_program(exe, tmp_path, None, timeout=30).check()


def test_stop_with_nonzero_code_text_fails(tmp_path, stub_script):
    """A "STOP <n>" line with n != 0 fails the check."""
    exe = stub_script(tmp_path / "prog.py", "import sys; sys.stderr.write('STOP 2\\n')")
    with pytest.raises(ProgramFailed):
        run_program(exe, tmp_path, None, timeout=30).check()


def test_stop_zero_is_ok(tmp_path, stub_script):
    """A "STOP 0" line passes the check."""
    exe = stub_script(tmp_path / "prog.py", "import sys; sys.stderr.write('STOP 0\\n')")
    run_program(exe, tmp_path, None, timeout=30).check()


def test_timeout_fails(tmp_path, stub_script):
    """A program running past the timeout fails with "timeout"."""
    exe = stub_script(tmp_path / "prog.py", "import time; time.sleep(30)")
    with pytest.raises(ProgramFailed, match="timeout"):
        run_program(exe, tmp_path, None, timeout=0.5)


def test_thread_env_is_set(tmp_path, stub_script):
    """The program runs with single-threaded BLAS and OpenMP environment."""
    exe = stub_script(tmp_path / "prog.py", "import os; print(os.environ['OPENBLAS_NUM_THREADS'])")
    assert run_program(exe, tmp_path, None, timeout=30).stdout.strip() == "1"


def test_memory_is_measured(tmp_path, stub_script):
    """The peak resident memory of the program is measured."""
    exe = stub_script(tmp_path / "prog.py",
                      "b = bytearray(200*1024*1024); b[::4096] = b'x'*len(b[::4096])")
    r = run_program(exe, tmp_path, None, timeout=60)
    assert r.max_rss_kb > 150 * 1024


def test_peak_rss_excludes_the_parent_python_process(tmp_path):
    """Peak memory does not include the parent.

    On Linux ru_maxrss survives execve: a child forked from this (large) Python process would
    report its RSS for a tiny program.
    """
    ballast = bytearray(100 * 1024 * 1024)
    ballast[::4096] = b"x" * len(ballast[::4096])
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 100 * 1024
    r = run_program(Path("/bin/true"), tmp_path, None, timeout=30)
    assert 0 < r.max_rss_kb < 10 * 1024


def test_killed_by_signal_gives_negative_returncode(tmp_path, stub_script):
    """A program killed by a signal has a negative return code and fails the check."""
    exe = stub_script(tmp_path / "prog.py",
                      "import os, signal; os.kill(os.getpid(), signal.SIGSEGV)")
    r = run_program(exe, tmp_path, None, timeout=30)
    assert r.returncode == -11
    with pytest.raises(ProgramFailed):
        r.check()


def test_relative_cwd_and_paths(tmp_path, stub_script, monkeypatch):
    """Relative executable, working directory and stdin paths are resolved correctly."""
    exe = stub_script(tmp_path / "prog.py", "import sys; print(sys.stdin.read())")
    (tmp_path / "wd").mkdir()
    (tmp_path / "in").write_text("xyz")
    monkeypatch.chdir(tmp_path)
    r = run_program(Path("prog.py"), Path("wd"), Path("in"), timeout=30).check()
    assert "xyz" in r.stdout


def test_launcher_rebuilt_for_other_c_compiler_version(tmp_path, mocker):
    """The launcher is rebuilt exactly when the C compiler version text changes."""
    version = "gcc 15.2\nline 2\n"

    def fake_run(cmd, **kw):
        """Fake subprocess.run: report the compiler version or create the compile output."""
        if cmd[1:] == ["--version"]:
            return subprocess.CompletedProcess(cmd, 0, stdout=version)
        Path(cmd[cmd.index("-o") + 1]).touch()
        return subprocess.CompletedProcess(cmd, 0)
    mocker.patch.object(runner.subprocess, "run", side_effect=fake_run)
    mocker.patch.object(runner, "_CACHE", tmp_path)
    exe1 = runner._launcher()
    assert runner._launcher() == exe1
    version = "gcc 15.2\nline 2 changed\n"
    assert runner._launcher() != exe1
