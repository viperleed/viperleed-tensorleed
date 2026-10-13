"""The baseline commit and the check that a TL tree is that baseline.

Shared by all generators; imports only the standard library (usable from
the tl-tests environment and from an environment with viperleed installed).
"""
import subprocess
from pathlib import Path

BASELINE_COMMIT = "22e3998"


class GeneratorError(RuntimeError):
    pass


def check_baseline_tree(tl, allow_other_commit=False):
    """Return the HEAD commit of tl after checking that tl is the top level
    of a git checkout, at BASELINE_COMMIT (unless allow_other_commit), and
    that TensErLEED/ has no local changes, untracked or git-ignored files
    included (always: the commit must describe the code that ran)."""
    def git(*args):
        return subprocess.run(["git", "-C", str(tl), *args], check=True,
                              capture_output=True, text=True).stdout.strip()
    top = git("rev-parse", "--show-toplevel")
    if Path(top).resolve() != Path(tl).resolve():
        raise GeneratorError(f"{tl} is not the top level of its git checkout ({top})")
    head = git("rev-parse", "HEAD")
    if not head.startswith(BASELINE_COMMIT) and not allow_other_commit:
        raise GeneratorError(f"{tl} is at {head[:7]}, not the baseline {BASELINE_COMMIT}")
    dirty = git("status", "--porcelain", "--ignored", "TensErLEED")
    if dirty:
        raise GeneratorError(f"{tl}/TensErLEED has local changes:\n{dirty}")
    return head
