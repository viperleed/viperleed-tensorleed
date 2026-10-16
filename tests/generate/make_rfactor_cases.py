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

The title lines of all frozen WEXPEL files and theory spectra are replaced
by 'tl-test <case>'. Every run is executed with the baseline rfactor
(legacy build: release, debug, debug without bounds checks) via the test
harness; the release run must succeed and write an R factor; debug failures
are reported. Gap inputs must run clean in the debug build without bounds
checks (see DEBUG_NO_BOUNDS; known gap bugs: tests/README.md, 'Known baseline issues').

Usage (in an environment with viperleed installed):
    python -m generate.make_rfactor_cases \
        [--workroot DIR] [--only CASE ...]
"""
import argparse
import re
import shutil
import sys
import tempfile
from pathlib import Path

import fortranformat as ff

from .make_refcalc_cases import (GeneratorError,
                                ensure_generator_tl, run_viperleed, set_run)
from harness.build import DEBUG, RELEASE, BuildFlags, build_program
from harness.paths import location
from harness.runner import run_program

TESTS_DIR = Path(__file__).parents[1]
CASES_DIR = TESTS_DIR / "cases" / "rfactor"
OUTPUTS = ["ROUT", "ROUTSHORT", "exp.column", "theo.column"]
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

def set_namelist_value(text, key, value):
    """Set KEY= value, in a namelist (viperleed layout: one key per line)."""
    if key == "WR":
        new, n = re.subn(r"(?m)^(\s*WR=\s*).*$", rf"\g<1>{value}", text)
    else:
        new, n = re.subn(rf"(?m)^(\s*{key}=\s*)[^,\n]*,", rf"\g<1>{value},", text)
    if n != 1:
        raise GeneratorError(f"expected one '{key}=' line in WEXPEL, found {n}")
    return new


def set_title(text, title):
    """Replace the title line (first line after the &NL3 ... &END block)."""
    lines = text.split("\n")
    i_nl3 = next(i for i, l in enumerate(lines) if l.strip() == "&NL3")
    i_end = next(i for i in range(i_nl3, len(lines)) if lines[i].strip() == "&END")
    lines[i_end + 1] = f" {title}"
    return "\n".join(lines)


def format_int_line(values):
    """25I4 lines as written by viperleed's writeWEXPEL."""
    out = ff.FortranRecordWriter("25I4").write(values) + "\n"
    if len(values) % 25 == 0:
        out += "\n"
    return out


def set_beam_arrays(text, kav, nbea, ibp, mittel):
    """Replace the KAV/NBEA/IBP/MITTEL lines between &NL1 ... &END and &NL2."""
    lines = text.split("\n")
    i_end = next(i for i, l in enumerate(lines) if l.strip() == "&END")
    i_nl2 = next(i for i, l in enumerate(lines) if l.strip() == "&NL2")
    block = "".join(format_int_line(v) for v in (kav, nbea, ibp, mittel))
    return "\n".join(lines[:i_end + 1] + block.rstrip("\n").split("\n") + lines[i_nl2:])


def split_exp_data(text):
    """WEXPEL -> (head up to and incl. format line, [beams], tail).

    Each beam is (name_line, fac_line_rest, [(E, I), ...]).
    """
    lines = text.split("\n")
    i_fmt = next(i for i, l in enumerate(lines) if re.fullmatch(r"\s*\(.*E.*\)\s*", l))
    head = lines[:i_fmt + 1]
    fmt = lines[i_fmt].strip()
    reader = ff.FortranRecordReader(fmt)
    beams, i = [], i_fmt + 1
    while i < len(lines) and lines[i].strip():
        name = lines[i]
        n, rest = int(lines[i + 1][:4]), lines[i + 1][4:]
        i += 2
        values = []
        while len(values) < 2 * n:
            values += [v for v in reader.read(lines[i]) if v is not None]
            i += 1
        values = values[:2 * n]
        beams.append((name, rest, list(zip(values[0::2], values[1::2]))))
    return head, fmt, beams, lines[i:]


def join_exp_data(head, fmt, beams, tail):
    """Inverse of split_exp_data: assemble the experimental-data part from its pieces."""
    writer = ff.FortranRecordWriter(fmt)
    out = list(head)
    for name, rest, points in beams:
        out.append(name)
        out.append(f"{len(points):4d}{rest}")
        flat = [x for p in points for x in p]
        out += writer.write(flat).split("\n")
    return "\n".join(out + tail)


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


def set_spectrum_title(text, title):
    """Replace the first (title) line of a theoretical spectrum file."""
    _, rest = text.split("\n", maxsplit=1)
    return title + "\n" + rest


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


def _write_case(name, description, covers, runs):
    """runs: [(run_name, wexpel_text, param_text, spectrum_text)]."""
    case_dir = CASES_DIR / name
    inputs = case_dir / "inputs"
    if inputs.exists():
        shutil.rmtree(inputs)
    inputs.mkdir(parents=True)
    title = f"tl-test {name}"
    toml = ["# Generated by tests/generate/make_rfactor_cases.py",
            'program = "rfactor"', f'description = "{description}"',
            "covers = [" + ", ".join(f'"{c}"' for c in covers) + "]", ""]
    written = {}
    for run, wexpel, param, spectrum in runs:
        (inputs / f"WEXPEL_{run}").write_text(set_title(wexpel, title))
        # share identical PARAM / spectrum files between runs
        for kind, content in (("PARAM", param),
                              ("spectrum", set_spectrum_title(spectrum, title))):
            key = (kind, content)
            if key not in written:
                fname = kind if not any(k[0] == kind for k in written) else f"{kind}_{run}"
                (inputs / fname).write_text(content)
                written[key] = fname
        toml += ["[[runs]]", f'name = "{run}"',
                 f'stdin = "inputs/{written[("spectrum", set_spectrum_title(spectrum, title))]}"',
                 f'param = "inputs/{written[("PARAM", param)]}"',
                 f'files = {{ WEXPEL = "inputs/WEXPEL_{run}" }}',
                 "outputs = [" + ", ".join(f'"{o}"' for o in OUTPUTS) + "]", ""]
    (case_dir / "case.toml").write_text("\n".join(toml))
    return case_dir


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
        raise GeneratorError(f"{base} missing; run make_refcalc_cases.py --only ag100 first")
    wexpel = (base / "rfactor-WEXPEL").read_text()
    param = (base / "rfactor-PARAM").read_text()
    spectrum = (base / "refcalc-fd.out").read_text()
    cases = []
    if not only or "ag100_rfactor" in only:
        runs = [(r, edit(wexpel), param, spectrum) for r, edit in RUN_EDITS.items()]
        cases.append(_write_case(
            "ag100_rfactor", "Ag(100) rfactor: R-factor types, smoothing, averaging, V0 range",
            ["pendry", "r2", "zj", "rs", "mixed-weights", "smoothing",
             "beam-averaging", "fractional-beams", "v0-range"], runs))
    if not only or "ag100_rfactor_errors" in only:
        err = run_error_calc(workroot, tl)
        cases.append(_write_case(
            "ag100_rfactor_errors", "Ag(100) rfactor on a multi-geometry superpos spectrum",
            ["mns"], [("mns", (err / "rfactor-WEXPEL").read_text(),
                       (err / "rfactor-PARAM").read_text(),
                       (err / "superpos-spec.out").read_text())]))
    if not only or "ag100_rfactor_paths" in only:
        # code paths not reached by the other cases
        cases.append(_write_case(
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
        cases.append(_write_case(
            "ag100_rfactor_short",
            "Ag(100) rfactor on short experimental beams (2, 23, 22 and 15 points)",
            ["der-ne2", "der-ne23", "der-short"],
            [(r, RUN_EDITS[r](short), param, spectrum) for r in ("pendry", "rs")]))
    if not only or "ag100_rfactor_gaps" in only:
        cases.append(_write_case(
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
    parser.add_argument("--verify-only", action="store_true",
                        help="only build/run the baseline on the committed cases")
    args = parser.parse_args(argv)
    args.workroot = args.workroot or location("generate_workroot")
    tl = ensure_generator_tl()
    if args.verify_only:
        case_dirs = [CASES_DIR / c for c in (args.only or sorted(
            p.name for p in CASES_DIR.iterdir() if p.is_dir()))]
    else:
        case_dirs = build_cases(args.workroot, tl, args.only)
    failed = False
    for case_dir in case_dirs:
        report = verify(case_dir, tl)
        for run, status in report:
            print(f"{case_dir.name}/{run}: " + "; ".join(f"{k} {v}" for k, v in status.items()))
            failed |= status["release"] != "OK"
            if "gap" in run:
                failed |= status["debug-nobounds"] != "OK"
        if case_dir.name == "ag100_rfactor":
            failed |= not check_pendry_against_viperleed(args.workroot, report)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
