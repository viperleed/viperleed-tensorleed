"""Line coverage (gcov) of the tier-1 runs on a TL tree (legacy build mode).

Builds ref-calc / rfactor with -O0 --coverage once per run (with the run's
PARAM), executes every tier-1 run, and sums the gcov line counts over all
runs. Reports per routine how many executable lines ran, and lists the
never-executed lines of chosen ranges (--lines FILE:A-B).

The list of known gaps in tests/README.md ("Known coverage gaps") comes from this script.

Usage (env tl-tests):
    python -m generate.coverage_report --tl-root <baseline> --work DIR
        [--cases CASE ...] [--reuse] [--lines lib.tleed.f:2176-2200 ...]
        [--never-only] [--missed N]
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from harness.cases import load_cases
from harness.runner import run_program

SOURCES = {"refcalc": [("src", "ref-calc.f"), ("lib", "lib.tleed.f")],
           "rfactor": [("lib", "rfacsb.f"), ("src", "rfactor.f")]}
EXTRA = {"refcalc": [("src", "GLOBAL")], "rfactor": []}
_ROUTINE = re.compile(r"^\s{6}[^!]*\b(SUBROUTINE|FUNCTION|PROGRAM)\s+(\w+)", re.I)
_GCOV_LINE = re.compile(r"\s*([^:]+):\s*(\d+):")


def run_with_coverage(tl_root, case, run, workdir):
    """Build case.program with --coverage in workdir, run it, return {file: {line: count}}."""
    tl = Path(tl_root) / "TensErLEED"
    if workdir.exists():
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)
    for d, f in SOURCES[case.program] + EXTRA[case.program]:
        shutil.copy(tl / d / f, workdir / f)
    shutil.copy(run.param, workdir / "PARAM")
    objs = []
    for _, f in SOURCES[case.program]:
        obj = f.replace(".", "_") + ".o"
        subprocess.run(["gfortran", "-O0", "--coverage", "-c", f, "-o", obj],
                       cwd=workdir, check=True, capture_output=True)
        objs.append(obj)
    subprocess.run(["gfortran", "--coverage", "-o", "prog", *objs, "-llapack", "-lblas"],
                   cwd=workdir, check=True, capture_output=True)
    for name, src in run.files.items():
        shutil.copy(src, workdir / name)
    res = run_program(workdir / "prog", workdir, run.stdin, timeout=900)
    if res.returncode != 0:
        print(f"  WARNING: {case.name}/{run.name} exit code {res.returncode}")
    counts = {}
    for _, f in SOURCES[case.program]:
        obj = f.replace(".", "_") + ".o"
        subprocess.run(["gcov", "-o", obj, f], cwd=workdir, capture_output=True)
        cnt = counts.setdefault(f, {})
        for line in (workdir / f"{f}.gcov").read_text(errors="replace").splitlines():
            m = _GCOV_LINE.match(line)
            if not m or m.group(1).strip() == "-" or int(m.group(2)) == 0:
                continue
            c = m.group(1).strip()
            cnt[int(m.group(2))] = 0 if c[0] in "#=" else int(c.rstrip("*"))
    return counts


def aggregate(tl_root, work, cases):
    agg = {}
    for kind in ("refcalc", "rfactor"):
        for case in load_cases(kind):
            if cases and case.name not in cases:
                continue
            for run in case.runs:
                counts = run_with_coverage(tl_root, case, run, work / f"{case.name}-{run.name}")
                for f, cnt in counts.items():
                    tot = agg.setdefault(f, {})
                    for ln, c in cnt.items():
                        tot[ln] = tot.get(ln, 0) + c
                print(f"{case.name}/{run.name}", flush=True)
    return agg


def source_path(tl_root, fname):
    sub = next(d for prog in SOURCES.values() for d, f in prog if f == fname)
    return Path(tl_root) / "TensErLEED" / sub / fname


def routines(lines):
    """[(first line, last line, name)] of the routines of a fixed-form source."""
    starts = [(i, m.group(2)) for i, l in enumerate(lines, 1)
              if not re.match(r"^[cC*!]", l) and (m := _ROUTINE.match(l))]
    starts.append((len(lines) + 1, None))
    return [(s, e - 1, name) for (s, name), (e, _) in zip(starts, starts[1:])]


def report(tl_root, agg, never_only):
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
