"""Generate frozen rfactor test cases.

Base inputs (rfactor-WEXPEL, rfactor-PARAM, refcalc-fd.out) come from the
viperleed work dir of the 'ag100' refcalc generator run
(<workroot>/ag100/calc/work, see make_refcalc_cases.py). The multi-geometry
case (MNS > 1) comes from an ag100 error calculation (RUN = 1 5, geometry
displacements of layer 1) run here in <workroot>/ag100_rfactor_errors.

Cases and runs:
  ag100_rfactor:        pendry, r2, zj, rs, mixed, smooth, averaging, v0_wide
  ag100_rfactor_errors: mns
  ag100_rfactor_gaps:   explicit_gap, auto_gap
  ag100_rfactor_paths:  theo_past_exp, permuted_weighted, r2_gap, zj_gap
  ag100_rfactor_short:  pendry, rs (short experimental beams, DER for NE < 24)

With --benchmark, the benchmark case bench_ag100_rfactor is generated instead,
into a temporary directory (many beams and geometries made from the committed
ag100_rfactor inputs, no viperleed needed; see rfactor_inputs.make_many_beams),
and its case.toml is written to tests/cases/benchmark/bench_ag100_rfactor. Its
inputs are not committed: the benchmark tests generate them at test time.

The title lines of all frozen WEXPEL files and theory spectra are replaced
by 'tl-test <case>'. Every run is executed with the baseline rfactor
(legacy build: release, debug, debug without bounds checks) via the test
harness; the release run must succeed and write an R factor; debug failures
are reported. Gap inputs must run clean in the debug build without bounds
checks (see DEBUG_NO_BOUNDS; known gap bugs: tests/README.md, 'Known baseline issues').

Usage (in an environment with viperleed installed):
    python -m generate.make_rfactor_cases \
        [--workroot DIR] [--only CASE ...] [--benchmark] [--verify-only]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import argparse
import re
import shutil
import sys
import tempfile
from pathlib import Path

from .case_files import set_namelist_value
from .rfactor_inputs import (BENCHMARK_DIR, RFACTOR_DIR as CASES_DIR, build_benchmark_case,
                             join_exp_data, set_beam_arrays, split_exp_data,
                             write_rfactor_case)
from .make_refcalc_cases import (GeneratorError,
                                 ensure_generator_tl, run_viperleed, set_run)
from harness.build import DEBUG, RELEASE, BuildFlags, build_program
from harness.paths import location
from harness.runner import run_program

# Debug build without bounds checks: the baseline's gap path fails the bounds
# checks for any gap input because of dummy-argument declarations
# (GAPIN: BM(NGAP) with NGAP = 0 on entry; GRID: A(NGP) with NGP undefined),
# while the actual arrays are large enough. Gap inputs must run clean with
# all other checks and FP traps (see tests/README.md, known baseline issues).
DEBUG_NO_BOUNDS = BuildFlags(
    "debug-nobounds", "-O0 -g -fcheck=all,no-bounds -ffpe-trap=invalid,zero,overflow")
ERROR_DISPLACEMENTS = "= GEO_DELTA\nAg L(1) z = -0.06 0.06 0.02\n"

# explicit gap: beam number (1-based exp beam), energy range (eV)
EXPLICIT_GAP = (2, 100.0, 120.0)
# auto gap: remove the experimental points of this beam in this range
AUTO_GAP = (3, 150.0, 165.0)
# ag100_rfactor_paths/theo_past_exp: experimental data end here (theory: 210 eV)
EXP_END = 190.0
# ag100_rfactor_paths/permuted_weighted: SYM codes of the theory beams
# (weights in RINTAV: 1 -> 1, 2 -> 2, 10 -> 8)
PERMUTED_SYM = [1, 1, 2, 10]
# ag100_rfactor_short: (beam, first energy (eV), number of points) kept of the
# experimental data. Each beam takes another DER branch (rfacsb.f): NE = 2,
# NE = 23 (the long formula at its threshold), NE = 22 and 15 (the short
# formula with its one-sided end points). Each window lies where the
# branches differ most for the baseline data.
SHORT_BEAMS = [(1, 120.0, 2), (2, 85.5, 23), (3, 150.0, 22), (4, 170.0, 15)]


# ------------------------------------------------------------ WEXPEL edits ---

def remove_exp_points(text, beam, emin, emax):
    """Remove the points of experimental beam (1-based) with energy in [emin, emax]."""
    head, fmt, beams, tail = split_exp_data(text)
    name, rest, points = beams[beam - 1]
    kept = [p for p in points if not emin <= p[0] <= emax]
    if len(kept) == len(points):
        raise GeneratorError(f"no exp points of beam {beam} in [{emin}, {emax}]")
    beams[beam - 1] = (name, rest, kept)
    return join_exp_data(head, fmt, beams, tail)


def keep_exp_points(text, beam, e_first, n):
    """Keep only n experimental points of beam (1-based), from e_first on."""
    head, fmt, beams, tail = split_exp_data(text)
    name, rest, points = beams[beam - 1]
    kept = [p for p in points if p[0] >= e_first][:n]
    if len(kept) != n:
        raise GeneratorError(f"beam {beam} has fewer than {n} exp points from {e_first} eV")
    beams[beam - 1] = (name, rest, kept)
    return join_exp_data(head, fmt, beams, tail)


def add_explicit_gap(text, beam, e1, e2):
    """GAP = 1 and the GAPIN block (I3 / I3,2F7.2) after the exp data."""
    text = set_namelist_value(text, "GAP", "     1")
    head, fmt, beams, tail = split_exp_data(text)
    gap = ["  1", f"{beam:3d}{e1:7.2f}{e2:7.2f}"]
    return join_exp_data(head, fmt, beams, gap + tail)


def truncate_exp(text, emax):
    """Drop the experimental points above emax (all beams)."""
    head, fmt, beams, tail = split_exp_data(text)
    beams = [(name, rest, [p for p in points if p[0] <= emax])
             for name, rest, points in beams]
    return join_exp_data(head, fmt, beams, tail)


def set_spectrum_sym(text, syms):
    """Set the SYM code of each beam of a theory spectrum (beam header lines,
    read by rfactor with FORMAT(I5,2F10.5,I3): SYM in columns 26-28)."""
    lines = text.split("\n")
    if len(syms) != int(lines[1]):
        raise GeneratorError(f"{len(syms)} SYM codes for {lines[1].strip()} beams")
    for i, sym in enumerate(syms):
        line = lines[2 + i]
        if not re.fullmatch(r"[ \d]{5}[ \d.-]{20}[ \d]{3}", line):
            raise GeneratorError(f"unexpected spectrum beam line {line!r}")
        lines[2 + i] = f"{line[:25]}{sym:3d}"
    return "\n".join(lines)


# ------------------------------------------------------------------ cases ---

def _wr(text, value):
    return set_namelist_value(text, "WR", value)


RUN_EDITS = {
    "pendry": lambda t: _wr(t, "0.,0.,1.,"),
    "r2": lambda t: _wr(t, "1.,0.,0.,"),
    "zj": lambda t: _wr(t, "0.,1.,0.,"),
    "rs": lambda t: _wr(t, "0.,0.,0.,"),
    "mixed": lambda t: _wr(t, "0.5,0.,0.5,"),
    "smooth": lambda t: set_namelist_value(t, "ISMOTH", "     2"),
    # theory beams 3 and 4 averaged into group 3, exp beams 3 and 4 likewise;
    # groups 2 and 3 labelled fractional (MITTEL = 2)
    "averaging": lambda t: set_beam_arrays(t, [1, 2, 3, 3], [1, 2, 3, 3],
                                           [1, 2, 3, 0], [1, 2, 2, 0]),
    "v0_wide": lambda t: set_namelist_value(
        set_namelist_value(t, "V01", "    -6.00"), "V02", "    10.00"),
}


def run_error_calc(workroot, tl):
    """ag100 error calculation (RUN = 1 5) -> its work dir."""
    src = workroot / "ag100" / "calc"
    calc = workroot / "ag100_rfactor_errors" / "calc"
    if calc.exists():
        shutil.rmtree(calc)
    calc.mkdir(parents=True)
    for f in ("POSCAR", "PARAMETERS", "VIBROCC", "PHASESHIFTS", "IVBEAMS",
              "EXPBEAMS.csv"):
        shutil.copy(src / f, calc / f)
    set_run(calc / "PARAMETERS", "1 5")
    (calc / "DISPLACEMENTS").write_text(ERROR_DISPLACEMENTS)
    run_viperleed(calc, tl, "-k", "-w", "work")
    return calc / "work"


def build_cases(workroot, tl, only):
    """Write the rfactor cases (all, or those named in only) and return their directories."""
    base = workroot / "ag100" / "calc" / "work"
    if not (base / "rfactor-WEXPEL").is_file():
        raise GeneratorError(f"{base} missing; run python -m generate.make_refcalc_cases"
                             " --only ag100 first (from tests/)")
    wexpel = (base / "rfactor-WEXPEL").read_text()
    param = (base / "rfactor-PARAM").read_text()
    spectrum = (base / "refcalc-fd.out").read_text()
    cases = []
    if not only or "ag100_rfactor" in only:
        runs = [(r, edit(wexpel), param, spectrum) for r, edit in RUN_EDITS.items()]
        cases.append(write_rfactor_case(
            "ag100_rfactor", "Ag(100) rfactor: R-factor types, smoothing, averaging, V0 range",
            ["pendry", "r2", "zj", "rs", "mixed-weights", "smoothing",
             "beam-averaging", "fractional-beams", "v0-range"], runs))
    if not only or "ag100_rfactor_errors" in only:
        err = run_error_calc(workroot, tl)
        cases.append(write_rfactor_case(
            "ag100_rfactor_errors", "Ag(100) rfactor on a multi-geometry superpos spectrum",
            ["mns"], [("mns", (err / "rfactor-WEXPEL").read_text(),
                       (err / "rfactor-PARAM").read_text(),
                       (err / "superpos-spec.out").read_text())]))
    if not only or "ag100_rfactor_paths" in only:
        # code paths not reached by the other cases
        cases.append(write_rfactor_case(
            "ag100_rfactor_paths",
            "Ag(100) rfactor: theory beyond the experiment, permuted beams with weights, "
            "R2 and ZJ with a gap",
            ["theory-past-experiment", "permuted-ibp", "averaging-weights", "r2-gap",
             "zj-gap"],
            [("theo_past_exp", truncate_exp(wexpel, EXP_END), param, spectrum),
             # exp beam 1 <-> theory beam 2 (integer vs fractional label),
             # theory beams 3 and 4 averaged with weights 2 and 8
             ("permuted_weighted",
              set_beam_arrays(wexpel, [1, 2, 3, 3], [1, 2, 3, 3], [2, 1, 3, 0],
                              [1, 2, 2, 0]),
              param, set_spectrum_sym(spectrum, PERMUTED_SYM)),
             ("r2_gap", _wr(add_explicit_gap(wexpel, *EXPLICIT_GAP), "1.,0.,0.,"),
              param, spectrum),
             ("zj_gap", _wr(add_explicit_gap(wexpel, *EXPLICIT_GAP), "0.,1.,0.,"),
              param, spectrum)]))
    if not only or "ag100_rfactor_short" in only:
        short = wexpel
        for cut in SHORT_BEAMS:
            short = keep_exp_points(short, *cut)
        cases.append(write_rfactor_case(
            "ag100_rfactor_short",
            "Ag(100) rfactor on short experimental beams (2, 23, 22 and 15 points)",
            ["der-ne2", "der-ne23", "der-short"],
            [(r, RUN_EDITS[r](short), param, spectrum) for r in ("pendry", "rs")]))
    if not only or "ag100_rfactor_gaps" in only:
        cases.append(write_rfactor_case(
            "ag100_rfactor_gaps", "Ag(100) rfactor with gaps in the experimental data",
            ["gap-explicit", "gap-auto-detected"],
            [("explicit_gap", add_explicit_gap(wexpel, *EXPLICIT_GAP), param, spectrum),
             ("auto_gap", remove_exp_points(wexpel, *AUTO_GAP), param, spectrum)]))
    return cases


# ----------------------------------------------------------- verification ---

def verify(case_dir, tl):
    """Run every run with the baseline release and debug builds."""
    import tomllib
    data = tomllib.loads((case_dir / "case.toml").read_text())
    report = []
    for run in data["runs"]:
        status = {}
        for flags in (RELEASE, DEBUG, DEBUG_NO_BOUNDS):
            exe = build_program(tl, "rfactor", case_dir / run["param"], flags)
            with tempfile.TemporaryDirectory() as tmp:
                shutil.copy(case_dir / run["files"]["WEXPEL"], Path(tmp) / "WEXPEL")
                res = run_program(exe, Path(tmp), case_dir / run["stdin"], timeout=600)
                rout = Path(tmp) / "ROUT"
                ok = True
                try:
                    res.check()
                except Exception as err:  # noqa: BLE001 - report all failures
                    ok, msg = False, str(err).splitlines()[0:4]
                if ok and not (rout.is_file() and "AVERAGE R-FACTOR =" in rout.read_text()):
                    ok, msg = False, ["no 'AVERAGE R-FACTOR =' in ROUT"]
                status[flags.name] = "OK" if ok else f"FAILED {' | '.join(msg)}"
                if flags is RELEASE and ok:
                    status["R"] = best_r(rout.read_text())
        report.append((run["name"], status))
    return report


def best_r(rout):
    """(V0 shift, R) of the best geometry from the ROUT summary."""
    m = re.search(r"INNER POTENTIAL SHIFT\s+(-?[\d.]+) EV, AVERAGE R-FACTOR =\s*([\d.]+)",
                  rout)
    return (float(m.group(1)), float(m.group(2))) if m else None


def check_pendry_against_viperleed(workroot, report):
    """pendry run == viperleed's own refcalc R factor."""
    base = workroot / "ag100" / "calc" / "work"
    (rfile,) = base.glob("R_OUT_refcalc_R=*")
    want = best_r(rfile.read_text())
    got = dict(report)["pendry"].get("R")
    print(f"pendry sanity check: harness {got}, viperleed {want} ({rfile.name})")
    return got == want


def main(argv=None):
    """Build the rfactor cases (or take the committed ones) and verify them with the baseline."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--workroot", type=Path, default=None,
                        help="default: generate_workroot (tests/harness/paths.py)")
    parser.add_argument("--only", nargs="+",
                        choices=["ag100_rfactor", "ag100_rfactor_errors", "ag100_rfactor_gaps",
                                 "ag100_rfactor_paths", "ag100_rfactor_short"])
    parser.add_argument("--benchmark", action="store_true",
                        help="the benchmark case instead (inputs generated into a temporary"
                             " directory; with --verify-only, case.toml is not written)")
    parser.add_argument("--verify-only", action="store_true",
                        help="only build/run the baseline on the committed cases")
    args = parser.parse_args(argv)
    args.workroot = args.workroot or location("generate_workroot")
    tl = ensure_generator_tl()
    with tempfile.TemporaryDirectory() as tmp:
        if args.benchmark:
            case_dirs = [build_benchmark_case(Path(tmp))]
            if not args.verify_only:
                shutil.copy(case_dirs[0] / "case.toml",
                            BENCHMARK_DIR / case_dirs[0].name / "case.toml")
        elif args.verify_only:
            case_dirs = [CASES_DIR / c for c in (args.only or sorted(
                p.name for p in CASES_DIR.iterdir() if p.is_dir()))]
        else:
            case_dirs = build_cases(args.workroot, tl, args.only)
        failed = False
        for case_dir in case_dirs:
            report = verify(case_dir, tl)
            for run, status in report:
                print(f"{case_dir.name}/{run}: "
                      + "; ".join(f"{k} {v}" for k, v in status.items()))
                failed |= status["release"] != "OK"
                if "gap" in run:
                    failed |= status["debug-nobounds"] != "OK"
            if case_dir.name == "ag100_rfactor":
                failed |= not check_pendry_against_viperleed(args.workroot, report)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
