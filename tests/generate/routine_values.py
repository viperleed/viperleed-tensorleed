"""The common part of the routine-check generators.

make_rfacsb_routine_values.py and make_refcalc_routine_values.py each build
a legacy driver against the baseline at every flag set of FLAG_SETS, run it
on their input sets, require all builds to give the same results, and write
the inputs and results to a values.json that their tier-1 test reads.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-11-11'
__license__ = 'GPLv3+'

import argparse
import json
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

from generate.baseline_tree import GeneratorError, check_baseline_tree

FLAG_SETS = ("-O2", "-O0")


def hex32(x):
    """The 8 hex digits (big-endian, as Fortran's Z8 format) of x rounded to
    a 32-bit real."""
    return struct.pack(">f", x).hex()


def run_driver(exe, stdin):
    """Run the driver exe on the stdin lines; return its stdout."""
    return subprocess.run([str(exe)], input="\n".join(stdin) + "\n", text=True,
                          capture_output=True, check=True).stdout


def results_per_flag_set(build, run, cases, workdir):
    """Return [(case_id, what, stdin, result)] for cases [(case_id, what,
    stdin)], after building the driver for every flag set with
    build(workdir, flags) -> executable and running it with
    run(exe, stdin) -> result. Raises GeneratorError if the results of any
    case differ between the flag sets."""
    exes = [build(workdir, flags) for flags in FLAG_SETS]
    rows = []
    for case_id, what, stdin in cases:
        results = [run(exe, stdin) for exe in exes]
        differing = [flags for flags, result in zip(FLAG_SETS, results)
                     if result != results[0]]
        if differing:
            raise GeneratorError(f"{case_id}: the legacy results of {differing} "
                                 f"differ from those of {FLAG_SETS[0]}")
        rows.append((case_id, what, stdin, results[0]))
    return rows


def generate(argv, description, build, run, cases, values, legacy_source):
    """Parse the command line, check the baseline tree, and write values.

    build(tl, workdir, flags) builds the legacy driver against the baseline
    tree tl; run(exe, stdin) returns the result of one case; legacy_source
    names the legacy file the driver links, for the "source" entry."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--tl-root", required=True, type=Path,
                        help="baseline TensErLEED tree")
    parser.add_argument("--allow-non-baseline", action="store_true",
                        help="accept a tree at another commit than the baseline")
    args = parser.parse_args(argv)
    tl = args.tl_root.resolve()
    commit = check_baseline_tree(tl, allow_other_commit=args.allow_non_baseline)
    with tempfile.TemporaryDirectory() as tmp:
        rows = results_per_flag_set(lambda workdir, flags: build(tl, workdir, flags),
                                    run, cases, Path(tmp))
    compiler = subprocess.run(["gfortran", "--version"], capture_output=True, text=True,
                              check=True).stdout.splitlines()[0]
    values.write_text(json.dumps({
        "source": f"legacy {legacy_source} at {commit[:7]}, {compiler}, "
                  f"flags {' and '.join(FLAG_SETS)} (identical)",
        "cases": [{"id": case_id, "covers": what, "stdin": stdin, "expected": result}
                  for case_id, what, stdin, result in rows]}, indent=1) + "\n")
    print(f"wrote {len(rows)} cases to {values}")


def main(generate_values):
    """Run generate_values() and turn a GeneratorError into an error exit."""
    try:
        generate_values()
    except GeneratorError as err:
        sys.exit(f"error: {err}")
