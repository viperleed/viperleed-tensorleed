"""Line coverage (gcov) of the tier-1 runs on a TL tree.

Builds ref-calc / rfactor with -O0 --coverage once per run, as the tier-1
harness does (through TensErLEED/Makefile if the tree has one, else the
legacy compile; with the run's PARAM while the sources need it), executes
every tier-1 run, and sums the gcov line counts over all runs. Reports
per routine how many executable lines ran, and lists the never-executed
lines of chosen ranges (--lines FILE:A-B). A gcov failure, or a source with
routines but no gcov output, stops the report.

The list of known gaps in tests/README.md ("Known coverage gaps") comes from this script.

Usage (env tl-tests):
    python -m generate.coverage_report --tl-root=<baseline> --work DIR
        [--cases CASE ...] [--reuse] [--lines lib.tleed.f:2176-2200 ...]
        [--never-only] [--missed N]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from harness.build import (BuildFlags, build_in, needs_param, program_sources,
                           uses_makefile)
from harness.cases import iter_runs
from harness.runner import run_program

COVERAGE = BuildFlags("coverage", "-O0 --coverage", "-llapack -lblas")
_ROUTINE = re.compile(r"^\s{6}[^!]*\b(SUBROUTINE|FUNCTION|PROGRAM)\s+(\w+)", re.I)
_GCOV_LINE = re.compile(r"\s*([^:]+):\s*(\d+):")
# a program unit with executable code; a module without one (constants
# only, such as kinds.f90) has no gcov output
_CODE_UNIT = re.compile(r"^[^!\n]*\b(SUBROUTINE|FUNCTION|PROGRAM)\s+\w+", re.I | re.M)


def run_with_coverage(tl_root, case, run, workdir):
    """Build case.program with --coverage in workdir, run it, return {file: {line: count}}."""
    workdir = workdir.resolve()
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    program = case.program
    build_in(workdir, Path(tl_root).resolve(), program, run.param, COVERAGE,
             uses_makefile(tl_root), needs_param(tl_root, program), keep_tree=True)
    for name, src in run.files.items():
        shutil.copy(src, workdir / name)
    res = run_program(workdir / program, workdir, run.stdin, timeout=900)
    if res.returncode != 0:
        print(f"  WARNING: {case.name}/{run.name} exit code {res.returncode}")
    # gcov resolves the source paths from the directory the build compiled in
    compiled_in = workdir / "tree" / "TensErLEED"
    if not compiled_in.is_dir():
        compiled_in = workdir
    for gcno in workdir.rglob("*.gcno"):
        proc = subprocess.run(["gcov", "-o", str(gcno.parent), str(gcno)],
                              cwd=compiled_in, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"gcov failed on {gcno} (exit code {proc.returncode}):\n"
                               f"{proc.stderr[-2000:]}")
    counts = {}
    for src in program_sources(tl_root, program):
        f = src.name
        gcov = compiled_in / f"{f}.gcov"
        if not gcov.is_file():
            if _CODE_UNIT.search(src.read_text(errors="replace")):
                raise RuntimeError(f"gcov wrote no {gcov.name} for {src}, which has code")
            continue  # no executable lines: a module of constants
        cnt = counts.setdefault(f, {})
        for line in gcov.read_text(errors="replace").splitlines():
            m = _GCOV_LINE.match(line)
            if not m or m.group(1).strip() == "-" or int(m.group(2)) == 0:
                continue
            c = m.group(1).strip()
            cnt[int(m.group(2))] = 0 if c[0] in "#=" else int(c.rstrip("*"))
    return counts


def aggregate(tl_root, work, cases):
    """Execution counts per source file and line, summed over the given runs of the cases."""
    agg = {}
    for case, run in iter_runs(only=cases):
        counts = run_with_coverage(tl_root, case, run, work / f"{case.name}-{run.name}")
        for f, cnt in counts.items():
            tot = agg.setdefault(f, {})
            for ln, c in cnt.items():
                tot[ln] = tot.get(ln, 0) + c
        print(f"{case.name}/{run.name}", flush=True)
    return agg


def source_path(tl_root, fname):
    """Path of the Fortran source file fname of either program in the tree tl_root."""
    return next(p for prog in ("refcalc", "rfactor")
                for p in program_sources(tl_root, prog) if p.name == fname)


def routines(lines):
    """[(first line, last line, name)] of the routines of a fixed-form source."""
    starts = [(i, m.group(2)) for i, l in enumerate(lines, 1)
              if not re.match(r"^[cC*!]", l) and (m := _ROUTINE.match(l))]
    starts.append((len(lines) + 1, None))
    return [(s, e - 1, name) for (s, name), (e, _) in zip(starts, starts[1:])]


def report(tl_root, agg, never_only):
    """Print, per source file and routine, how many executable lines ran."""
    for f, cnt in agg.items():
        lines = source_path(tl_root, f).read_text(errors="replace").splitlines()
        executable = len(cnt)
        hit = sum(1 for c in cnt.values() if c > 0)
        print(f"==== {f}: {hit}/{executable} executable lines ran "
              f"({100 * hit / executable:.1f} %)")
        for s, e, name in routines(lines):
            ex = [ln for ln in range(s, e + 1) if ln in cnt]
            n_hit = sum(1 for ln in ex if cnt[ln] > 0)
            if ex and (not never_only or n_hit == 0):
                print(f"  {name:24s} {s:5d}-{e:5d} ran {n_hit:4d}/{len(ex):4d}"
                      f"{'  NEVER' if n_hit == 0 else ''}")


def missed_ranges(tl_root, agg, min_lines):
    """Contiguous blocks of never-executed lines (>= min_lines), per routine."""
    for f, cnt in agg.items():
        lines = source_path(tl_root, f).read_text(errors="replace").splitlines()
        print(f"==== {f}: never-executed blocks of >= {min_lines} executable lines")
        for s, e, name in routines(lines):
            block = []
            for ln in [ln for ln in range(s, e + 1) if ln in cnt] + [None]:
                if ln is not None and cnt[ln] == 0:
                    block.append(ln)
                    continue
                if len(block) >= min_lines:
                    print(f"  {name:20s} {block[0]:5d}-{block[-1]:5d} ({len(block)} lines)")
                block = []


def show_lines(tl_root, agg, spec):
    """Print the source lines of spec ("<file>:<first>-<last>") with their execution counts."""
    fname, rng = spec.split(":")
    a, b = (int(x) for x in rng.split("-"))
    lines = source_path(tl_root, fname).read_text(errors="replace").splitlines()
    cnt = agg.get(fname, {})
    print(f"---- {spec}")
    for ln in range(a, b + 1):
        c = cnt.get(ln)
        mark = "     -" if c is None else ("#####" if c == 0 else f"{c:6d}")
        print(f"{mark} {ln:5d}: {lines[ln - 1]}")


def main(argv=None):
    """Collect (or reuse) the coverage counts of the tier-1 runs and print the reports."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tl-root", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--cases", nargs="*", default=[])
    parser.add_argument("--reuse", action="store_true",
                        help="reuse <work>/coverage.json instead of running")
    parser.add_argument("--lines", nargs="*", default=[])
    parser.add_argument("--never-only", action="store_true",
                        help="list only routines that never ran")
    parser.add_argument("--missed", type=int, metavar="N", default=0,
                        help="list never-executed blocks of at least N executable lines")
    args = parser.parse_args(argv)
    store = args.work / "coverage.json"
    if args.reuse:
        agg = {f: {int(k): v for k, v in d.items()}
               for f, d in json.loads(store.read_text()).items()}
    else:
        agg = aggregate(args.tl_root, args.work, args.cases)
        store.write_text(json.dumps(agg))
    report(args.tl_root, agg, args.never_only)
    if args.missed:
        missed_ranges(args.tl_root, agg, args.missed)
    for spec in args.lines:
        show_lines(args.tl_root, agg, spec)
    return 0


if __name__ == "__main__":
    sys.exit(main())
