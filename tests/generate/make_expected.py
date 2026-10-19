"""Generate the expected outputs of all tier-1 cases from the baseline.

Runs every run of every case (tests/cases/<kind>/<case>/case.toml) with the
release build of --tl-root and copies the outputs listed in run.outputs to
expected/<run>/. Writes tests/cases/BASELINE.txt (commit, compiler, BLAS/
LAPACK, flags, date).

Regenerating expected outputs changes the ground truth of the test suite:
only allowed for a deliberate bug fix/behaviour change or a tolerance shown
to be too tight for legitimate floating-point reordering; the PR that
changes reference data states why (tests/README.md, 'Rules').

Usage (env tl-tests):
    python -m generate.make_expected --tl-root=<baseline tree> \
        [--only CASE ...] [--allow-non-baseline]
"""
import argparse
import datetime
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from harness.build import RELEASE
from harness.cases import CASES_DIR, execute_run, iter_runs

from .baseline_tree import BASELINE_COMMIT, GeneratorError, check_baseline_tree


def _cmd(args):
    return subprocess.run(args, capture_output=True, text=True, shell=isinstance(args, str)).stdout.strip()


def write_baseline_info(head, allowed_non_baseline):
    """Write cases/BASELINE.txt: the commit, toolchain, flags and date of the expected outputs."""
    lines = [
        f"commit: {head}",
        "tree: --tl-root, TensErLEED/ unmodified at this commit",
        f"compiler: {_cmd(['gfortran', '--version']).splitlines()[0]}",
        "blas/lapack:",
        *("  " + l.strip() for l in _cmd("ldconfig -p | grep -E 'libblas|liblapack'").splitlines()),
        f"compile flags: {RELEASE.compile}",
        f"link flags: {RELEASE.link}",
        f"date: {datetime.date.today().isoformat()}",
    ]
    if allowed_non_baseline:
        lines.append(f"WARNING: generated with --allow-non-baseline (expected {BASELINE_COMMIT})")
    (CASES_DIR / "BASELINE.txt").write_text("\n".join(lines) + "\n")


def generate(case, run, tl_root):
    """Run one run with the release build, store its outputs as expected; return their names."""
    out_dir = case.path / "expected" / run.name
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        execute_run(case, run, tl_root, RELEASE, work)
        files = set()
        for pattern in run.outputs:
            matches = [p for p in work.glob(pattern) if p.is_file()]
            if not matches:
                raise RuntimeError(f"{case.name}/{run.name}: no output matches {pattern!r}")
            files.update(matches)
        if out_dir.exists():
            shutil.rmtree(out_dir)
        out_dir.mkdir(parents=True)
        for f in sorted(files):
            shutil.copy(f, out_dir / f.name)
    return sorted(f.name for f in files)


def main(argv=None):
    """Regenerate the expected outputs of the selected tier-1 runs from the baseline tree."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tl-root", type=Path, required=True)
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--allow-non-baseline", action="store_true")
    args = parser.parse_args(argv)
    try:
        head = check_baseline_tree(args.tl_root, args.allow_non_baseline)
    except GeneratorError as err:
        sys.exit(f"{err}; --allow-non-baseline overrides only the commit check")
    for case, run in iter_runs(only=args.only):
        names = generate(case, run, args.tl_root)
        print(f"{case.program}/{case.name}/{run.name}: {len(names)} files ({', '.join(names[:6])}"
              f"{', ...' if len(names) > 6 else ''})")
    write_baseline_info(head, not head.startswith(BASELINE_COMMIT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
