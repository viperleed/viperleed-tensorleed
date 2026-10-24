"""Derive hand-made test cases from the frozen inputs of existing cases.

Some input paths of TensErLEED are never produced by viperleed (README,
'Known coverage gaps'). Each case here is a deterministic text edit of the
committed inputs of an existing case of the same system, so the physics
stays the same and only the input form differs:

- ir100_2x1_o_unsorted: the sublayers of the composite layer type 1 of
  ir100_2x1_o (O, Ir, Ir) are listed as Ir, Ir, O, i.e. not in ascending
  order along the surface normal (+x). SRTLAY has to re-sort them. The
  site type and the tensor file name of each sublayer move with it, so the
  structure and the tensor files are those of the source case: the
  expected outputs of E118_L07 and E110-122_L07 equal those of the source
  case apart from the title line.
  Run E118_L07_three_heights also moves the Ir N=2 sublayer deeper
  (DEEPER_IR_X; x grows into the crystal) and lists the sublayers as
  Ir N=2, O, Ir N=3: three distinct heights, deepest (largest x) first,
  so that SRTLAY has to move two sublayers and pick the smallest x of
  several candidates. Its structure
  differs from the source case, so its outputs are its own.
- ag100_psformat2: the phase shifts of ag100 in the alternative format
  read by READIN for PSFORMAT = 2 ('Rundgren_v1.7' in the FIN). Format 1 gives,
  for each energy (Hartree, F7.4), the phase shifts of each element on two
  lines of 10F7.4; format 2 gives, for each element, two header lines and
  one line per energy with the energy in eV (F8.3) and up to 20 phase
  shifts (F8.4). The phase shifts are copied digit for digit; the energies
  are converted with the HARTREE constant of TensErLEED (27.211396 eV) and
  rounded to 0.001 eV, so after READIN divides them by HARTREE they differ
  from the format-1 values by up to 2e-5 Hartree. fd.out is identical to
  that of ag100 apart from the title line; amp.out and the T_* files
  differ slightly (own expected outputs).
  The format-2 header holds the inner-potential parameters (format 202:
  PSL1, EM, C1..C8). ref-calc uses EM to choose the formula for the inner
  potential: EM = 0 selects V = WORKFN - max(C1 + C2/sqrt(E + WORKFN + C3),
  C4), the formula format-1 runs use (EM is never set there), so EM = 0
  and C1..C4 are those of the format-1 header. C5..C8 only enter for EM > 0
  and are 0. PSL1 is read but never used; it is set to the number of phase
  shifts per energy. The title line (PAR1) is only printed.
- ag100_rfactor_hand: the pendry run of ag100_rfactor, edited in WEXPEL.
  pseudo_exp sets EOT = 1, so rfactor reads the experiment with READT as a
  pseudo-experiment in the layout of a theory spectrum: one geometry of the
  ag100_rfactor_errors spectrum (PSEUDO_EXP_GEOMETRY) replaces the
  experimental data after the title. The (2 0) and (2 1) beams start with
  zero intensities (not yet emerging), which READT skips. norm_max sets
  NORM = 2, so exp.column and theo.column are normalised to the highest
  intensity of each beam instead of the integral.
- ag100_rfactor_params: rfactor input parameters, one run each, from the
  ag100_rfactor pendry run unless noted. The 2001 CPC paper (Table 7 and
  p. 415) describes them only for the R-factor input of search (rf.info).
  theory_beam_unmatched: experimental beam 3 is removed, so theory beam 3
  has no experimental counterpart; as viperleed writes it, KAV = 0 for that
  beam, MITTEL = 0 after the experimental beams and MNBED = 3. CPC: KAV(k) = j
  puts theory beam k into group j, whose beams are averaged; KAV = 0 is not
  described. weights_all: WR = 0.1, 0.2, 0.3, 0.4 weights R2, R_ZJ, R_P and
  R_S in the average R factor; not in CPC, where search uses either R_P or
  R2. smooth_gap: the explicit_gap run of ag100_rfactor_gaps (one gap) with
  ISMOTH = 2; CPC: ISMOTH smooths the experimental data that many times
  with a three-point formula; gaps are not described. all_fractional:
  MITTEL = 2 for every beam, so there are no integer beams (NINT = 0); CPC:
  MITTEL groups the beams into two sets for R-factor averaging, e.g.
  integer and fractional. ipr2: IPR = 2 with ISMOTH = 1, so the smoothed
  curves are printed as well; CPC: IPR (0, 1, 2) sets the amount of
  standard output, for debugging. gaps3_all: GAPS3, one explicit gap in
  each of experimental beams 1, 2 and 4, with all four weights of
  weights_all; MNGAP = 3. Each beam has at most one gap: the gap integrals
  of several gaps in one beam are overwritten in legacy (README, known
  baseline issues). autogap_edge_r2_v0std: the experimental points of beam
  4 from 205 to 212 eV are removed (EDGE_HOLE), with R2 only; rfactor
  detects the hole as a gap, which lies outside the energy interval common
  to theory and experiment for some inner-potential shifts.
  gap_below_theory_r2: an explicit gap LOW_GAP just above the first data
  point of beam 1, V0 from -4 to +4 eV, R2 only; for the larger shifts the
  gap lies partly before the start of the theory. gap_edge_hi_allw: an
  explicit gap HIGH_GAP near the upper end of beam 1, V0 from -6 to +10 eV,
  all four weights of weights_all; for the larger shifts the gap lies
  partly beyond the end of the theory. no_overlap: all experimental
  points below NO_OVERLAP_EXP_MIN are removed and EMAX is raised, so no beam
  overlaps the theory (up to 216 eV); legacy ends normally and prints
  asterisks for the average R factors in ROUT and stdout.
  gap_below_theory_r2 and no_overlap are not run in the debug build
  (README, 'Rules'): their inputs drive undefined behaviour that the legacy
  code shares, named in case.toml. autogap_edge_r2_v0std,
  gap_below_theory_r2 and gap_edge_hi_allw are also run with a poisoned
  heap (README, 'Rules').
- ag100_rfactor_beams: the ag100_rfactor pendry run with N_BEAMS
  experimental and theory beams instead of 4, each a copy of one of the 4
  beams with its intensities scaled by a smooth factor that differs between
  beams. KAV, NBEA, IBP and MITTEL are read as records of 25 entries: nb25
  fills one record exactly (followed by a blank line, as viperleed writes
  it), nb30 needs two.

Writes tests/cases/<program>/<case>/{case.toml, inputs/}; the expected
outputs come from make_expected.py. With --check, the cases are derived
into a temporary directory and compared with the committed ones.

Usage (env tl-tests):
    python -m generate.make_hand_cases [--only CASE ...] [--check]
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-14'
__license__ = 'GPLv3+'

import argparse
import math
import re
import sys
import tempfile
import tomllib
from pathlib import Path

from .case_files import set_first_line_title, set_namelist_value, set_wexpel_title, write_case

CASES_DIR = Path(__file__).parents[1] / "cases"
HEADER = "# Generated by tests/generate/make_hand_cases.py"
# ag100_rfactor_hand/pseudo_exp: geometry of the ag100_rfactor_errors spectrum
# used as the pseudo-experiment (1: layer 1 displaced by -0.06 A; the
# ag100_rfactor spectrum is geometry 4, no displacement)
PSEUDO_EXP_GEOMETRY = 1
HARTREE = 27.211396   # eV, as in TensErLEED's GLOBAL
# ir100_2x1_o_unsorted/E118_L07_three_heights: x (A, along the surface
# normal, growing into the crystal) of the Ir N=2 sublayer, moved deeper
# from 1.1504 so that the three sublayers of the composite layer have
# distinct heights
DEEPER_IR_X = 1.2004
# ag100_rfactor_params/gaps3_all: explicit gaps, (experimental beam, first and
# last energy in eV), one in each of three beams
GAPS3 = ((1, 70.0, 75.0), (2, 120.0, 126.0), (4, 190.0, 195.0))
# ag100_rfactor_params/autogap_edge_r2_v0std: experimental points removed,
# (experimental beam, first and last energy in eV)
EDGE_HOLE = (4, 205.0, 212.0)
# ag100_rfactor_params/gap_below_theory_r2: explicit gap (experimental beam,
# first and last energy in eV) just above the first data point of beam 1
# (48 eV; the theory starts at 45 eV), and the inner-potential range (eV)
LOW_GAP = (1, 48.5, 52.5)
LOW_GAP_V0 = (-4.0, 4.0)
# ag100_rfactor_params/gap_edge_hi_allw: explicit gap (experimental beam,
# first and last energy in eV) near the upper end of beam 1 (the theory ends
# at 216 eV), and the inner-potential range (eV)
HIGH_GAP = (1, 208.0, 214.0)
HIGH_GAP_V0 = (-6.0, 10.0)
# ag100_rfactor_params/no_overlap: experimental points below this energy (eV)
# are removed from every beam, and EMAX is set to NO_OVERLAP_EMAX
NO_OVERLAP_EXP_MIN = 220.0
NO_OVERLAP_EMAX = 300.0
# ag100_rfactor_beams: number of experimental and theory beams per run
N_BEAMS = (25, 30)


def _source_runs(program, source, run_names):
    """The runs of an existing case and the texts of their input files by path."""
    case_dir = CASES_DIR / program / source
    runs = {r["name"]: r for r in tomllib.loads((case_dir / "case.toml").read_text())["runs"]}
    runs = [runs[name] for name in run_names]
    paths = {p for run in runs
             for p in (run["stdin"], run["param"], *run.get("files", {}).values())}
    return runs, {path: (case_dir / path).read_text() for path in sorted(paths)}


def permute_sublayers(fin, ltype, order):
    """Reorder the sublayers of layer type `ltype` in a ref-calc FIN.

    order[i] is the (1-based) source sublayer that becomes sublayer i + 1.
    Permutes the sublayer lines of the layer-type block (site type and
    position) and the tensor file names of every stacked layer of that type.
    """
    lines = fin.split("\n")
    n = len(order)

    def permute(first):
        """Reorder the n lines starting at index `first` according to `order`."""
        block = lines[first:first + n]
        lines[first:first + n] = [block[k - 1] for k in order]

    start = lines.index(f"-   layer type {ltype} ---")
    if not re.match(rf"\s*{n}\s+number of Bravais sublayers in layer {ltype}\b",
                    lines[start + 2]):
        raise ValueError(f"layer type {ltype} does not have {n} sublayers")
    permute(start + 3)
    stacked = [i for i, l in enumerate(lines)
               if re.search(rf": layer type {ltype}, interlayer vector below", l)]
    if not stacked:
        raise ValueError(f"layer type {ltype} is not stacked")
    for i in stacked:
        if lines[i + 1].split()[0] == "1":       # tensor output: one file name per sublayer
            names = lines[i + 2:i + 2 + n]
            if not all(f"sublayer {k}" in l for k, l in enumerate(names, 1)):
                raise ValueError(f"tensor file names after line {i + 1} not found")
            # (file name with its padding, comment)
            fields = [re.match(r"(\S+\s+)(.*)", l).groups() for l in names]
            for j, k in enumerate(order):
                name, comment = fields[k - 1][0], fields[j][1]
                width = max(len(name), len(fields[j][0]))
                lines[i + 2 + j] = name.strip().ljust(width) + comment
    return "\n".join(lines)


def set_sublayer_x(fin, ltype, sublayer, x):
    """Set the x position (along the surface normal) of one sublayer of layer
    type `ltype` in a ref-calc FIN, in the fixed columns of the sublayer line
    (site type I3, then x, y, z as F9.4)."""
    lines = fin.split("\n")
    i = lines.index(f"-   layer type {ltype} ---") + 2 + sublayer
    lines[i] = lines[i][:3] + f"{x:9.4f}" + lines[i][12:]
    return "\n".join(lines)


def phase_shifts_to_format2(fin):
    """Rewrite the phase-shift block of a ref-calc FIN from PSFORMAT 1 to 2.

    The title of the phase shifts (PAR1, read as A8, A70 after a blank PAR0)
    names the title of the FIN.
    """
    lines = fin.split("\n")
    ifmt = next(i for i, l in enumerate(lines) if l.split()[1:2] == ["PSFORMAT"])
    if lines[ifmt][:3] != "  1":
        raise ValueError("phase shifts are not in format 1")
    head = ifmt + 2
    nel = int(lines[head][:3])
    c4, c1, c2, c3 = (float(lines[head][4 + 8 * k:12 + 8 * k]) for k in range(4))
    end = next(i for i in range(head, len(lines)) if lines[i].startswith("---"))
    block = lines[head + 1:end]
    per_energy = 1 + 2 * nel
    if len(block) % per_energy:
        raise ValueError("phase-shift block is not NPSI x (energy + NEL x 2 lines)")
    energies, shifts = [], [[] for _ in range(nel)]
    for k in range(0, len(block), per_energy):
        energies.append(float(block[k]))
        for iel in range(nel):
            text = block[k + 1 + 2 * iel] + block[k + 2 + 2 * iel]
            shifts[iel].append([text[j:j + 7] for j in range(0, len(text), 7)])
    nps = len(shifts[0][0])
    if nps > 20 or any(len(s) != nps for el in shifts for s in el):
        raise ValueError("unexpected number of phase shifts per energy")
    new = [f"{nel:3d}" + " " * 42 + "NEL",
           f"{'':8}{f'phase shifts of {lines[0].strip()}, Rundgren_v1.7 layout':<70}",
           f"{'':3}{nps:3d}{'':2}{0.0:6.2f}{'':1}" + "".join(f"{c:7.2f}" for c in (c1, c2, c3, c4))
           + f"{'':2}{0.0:7.2f}" + "".join(f"{0.0:8.4f}" for _ in range(3))]
    for iel in range(nel):
        new += [f"element {iel + 1}",
                "  E (eV)" + "".join(f"{f'l={l}':>8}" for l in range(nps))]
        new += [f"{e * HARTREE:8.3f}" + "".join(f"{float(p):8.4f}" for p in ps)
                for e, ps in zip(energies, shifts[iel])]
    lines[ifmt] = "  2" + lines[ifmt][3:]
    lines[head:end] = new
    return "\n".join(lines)


def ag100_psformat2():
    """Spec of ag100_psformat2: the phase shifts of ag100 in PSFORMAT 2."""
    case, source = "ag100_psformat2", "ag100"
    runs, inputs = _source_runs("refcalc", source, ["E045_L06", "E141_L08"])
    for run in runs:
        inputs[run["stdin"]] = phase_shifts_to_format2(
            set_first_line_title(inputs[run["stdin"]], f"tl-test {case}"))
    return {"program": "refcalc", "name": case,
            "description": "Ag(100) as ag100, phase shifts in the alternative format "
                           "(PSFORMAT = 2, energies in eV)",
            "covers": ["bravais", "psformat-2", "tensors"],
            "runs": runs, "inputs": inputs}


def ir100_2x1_o_unsorted():
    """Spec of ir100_2x1_o_unsorted: ir100_2x1_o with out-of-order sublayers."""
    case, source = "ir100_2x1_o_unsorted", "ir100_2x1_o"
    runs, inputs = _source_runs("refcalc", source, ["E118_L07", "E110-122_L07"])
    for run in runs:
        inputs[run["stdin"]] = permute_sublayers(
            set_first_line_title(inputs[run["stdin"]], f"tl-test {case}"),
            ltype=1, order=[2, 3, 1])
    (source_run,), source_inputs = _source_runs("refcalc", source, ["E118_L07"])
    three = {**source_run, "name": "E118_L07_three_heights",
             "stdin": "inputs/FIN_E118_L07_three_heights"}
    fin = set_first_line_title(source_inputs[source_run["stdin"]], f"tl-test {case}")
    inputs[three["stdin"]] = permute_sublayers(
        set_sublayer_x(fin, ltype=1, sublayer=2, x=DEEPER_IR_X), ltype=1, order=[2, 1, 3])
    runs.append(three)
    return {"program": "refcalc", "name": case,
            "description": "Ir(100)-(2x1)-O as ir100_2x1_o, sublayers of the composite "
                           "top layer listed out of order along the surface normal",
            "covers": ["composite-layer", "unsorted-sublayers", "tensors"],
            "runs": runs, "inputs": inputs}


def wexpel_head(wexpel):
    """The lines of an rfactor WEXPEL up to and including the title (the line
    after the &NL3 block), without the experimental data after it."""
    lines = wexpel.split("\n")
    i_title = lines.index(" &END", lines.index(" &NL3")) + 1
    return "\n".join(lines[:i_title + 1]) + "\n"


def pseudo_experiment(spectrum):
    """One geometry of a theory spectrum, in the layout READT reads.

    The theory spectrum (title, number of beams, beam lines, then for each
    energy one record per geometry: energy, geometry number, intensities) is
    the layout of READT without its title line. Only the records of
    PSEUDO_EXP_GEOMETRY are kept.
    """
    lines = spectrum.split("\n")
    nbeams = int(lines[1])
    out = lines[1:2 + nbeams]
    for line in lines[2 + nbeams:]:
        if line.strip() and round(float(line[7:14]) * 1e4) == PSEUDO_EXP_GEOMETRY:
            out += [line, ""]
    return "\n".join(out) + "\n"


def ag100_rfactor_hand():
    """Spec of ag100_rfactor_hand: a pseudo-experiment and max-normalised plots."""
    case, source = "ag100_rfactor_hand", "ag100_rfactor"
    (pendry,), inputs = _source_runs("rfactor", source, ["pendry"])
    inputs[pendry["stdin"]] = set_first_line_title(inputs[pendry["stdin"]], f"tl-test {case}")
    wexpel = set_wexpel_title(inputs.pop(pendry["files"]["WEXPEL"]), f"tl-test {case}")
    (mns,), mns_inputs = _source_runs("rfactor", "ag100_rfactor_errors", ["mns"])
    pseudo = {**pendry, "name": "pseudo_exp", "files": {"WEXPEL": "inputs/WEXPEL_pseudo_exp"}}
    inputs["inputs/WEXPEL_pseudo_exp"] = (wexpel_head(set_namelist_value(wexpel, "EOT", 1))
                                          + pseudo_experiment(mns_inputs[mns["stdin"]]))
    norm_max = {**pendry, "name": "norm_max", "files": {"WEXPEL": "inputs/WEXPEL_norm_max"}}
    inputs["inputs/WEXPEL_norm_max"] = set_namelist_value(wexpel, "NORM", 2)
    return {"program": "rfactor", "name": case,
            "description": "Ag(100) rfactor as ag100_rfactor pendry: a pseudo-experiment "
                           "(another geometry of the error calculation), plot "
                           "normalisation to the highest intensity",
            "covers": ["pseudo-experiment", "plot-norm-max"],
            "runs": [pseudo, norm_max], "inputs": inputs}


def set_beam_lines(wexpel, kav, nbea, ibp, mittel):
    """Replace the beam-index lines between &NL1 and &NL2 of an rfactor WEXPEL
    (KAV, NBEA, IBP, MITTEL; format 25I4, 25 entries per line). As viperleed
    writes them, a list whose last line is full is followed by a blank line."""
    lines = wexpel.split("\n")
    first, nl2 = lines.index(" &END") + 1, lines.index("&NL2")
    block = []
    for row in (kav, nbea, ibp, mittel):
        block += ["".join(f"{k:4d}" for k in row[i:i + 25]) for i in range(0, len(row), 25)]
        if len(row) % 25 == 0:
            block.append("")
    lines[first:nl2] = block
    return "\n".join(lines)


def drop_exp_beam(wexpel, beam):
    """Remove experimental beam `beam` (1-based) from an rfactor WEXPEL: its
    label line, the NEE/FAC line and its data, up to the next label."""
    lines = wexpel.split("\n")
    labels = [i for i, l in enumerate(lines) if re.fullmatch(r"\*\(.*\)\*\s*", l)]
    if beam >= len(labels):
        raise ValueError("only beams followed by another beam can be dropped")
    del lines[labels[beam - 1]:labels[beam]]
    return "\n".join(lines)


def drop_exp_points(wexpel, beam, emin, emax):
    """Remove the points with energies in [emin, emax] from experimental beam
    `beam` (1-based) of an rfactor WEXPEL. The remaining (energy, intensity)
    fields (E14.5) are kept digit for digit, 12 per line, and the number of
    points (I4) on the NEE/FAC line is updated."""
    lines = wexpel.split("\n")
    labels = [i for i, l in enumerate(lines) if re.fullmatch(r"\*\(.*\)\*\s*", l)]
    start = labels[beam - 1] + 1
    end = labels[beam] if beam < len(labels) else lines.index("", start)
    fields = [l[k:k + 14] for l in lines[start + 1:end] for k in range(0, len(l), 14)]
    kept = [f for k in range(0, len(fields), 2) if not emin <= float(fields[k]) <= emax
            for f in fields[k:k + 2]]
    lines[start:end] = ([f"{len(kept) // 2:4d}{lines[start][4:]}"]
                        + ["".join(kept[k:k + 12]) for k in range(0, len(kept), 12)])
    return "\n".join(lines)


def set_param(param, name, value):
    """Set one PARAMETER value (NAME = value) in a PARAM file."""
    new, n = re.subn(rf"\b{name} = \d+", f"{name} = {value}", param)
    if n != 1:
        raise ValueError(f"expected one '{name} =' in PARAM, found {n}")
    return new


def add_gaps(wexpel, gaps):
    """GAP = 1 and the gap list read by GAPIN (I3, then I3,2F7.2 per gap) after
    the experimental data, before the trailing blank line."""
    lines = set_namelist_value(wexpel, "GAP", 1).split("\n")
    if lines[-2:] != ["", ""]:
        raise ValueError("WEXPEL must end with a blank line")
    block = [f"{len(gaps):3d}"] + [f"{beam:3d}{e1:7.2f}{e2:7.2f}" for beam, e1, e2 in gaps]
    return "\n".join(lines[:-2] + block + lines[-2:])


def ag100_rfactor_params():
    """Spec of ag100_rfactor_params: rfactor input parameters, one run each."""
    case = "ag100_rfactor_params"
    (pendry,), inputs = _source_runs("rfactor", "ag100_rfactor", ["pendry"])
    title = f"tl-test {case}"
    inputs[pendry["stdin"]] = set_first_line_title(inputs[pendry["stdin"]], title)
    wexpel = set_wexpel_title(inputs.pop(pendry["files"]["WEXPEL"]), title)
    runs = []

    def add(name, text, **source):
        path = f"inputs/WEXPEL_{name}"
        inputs[path] = text
        runs.append({**pendry, **source, "name": name, "files": {"WEXPEL": path}})

    # experimental beam 3, (2 0), removed: theory beam 3 has no experimental
    # counterpart, and viperleed writes KAV = 0 for it, MITTEL = 0 after the
    # experimental beams, and PARAM with MNBED = 3
    inputs["inputs/PARAM_unmatched"] = set_param(inputs[pendry["param"]], "MNBED", 3)
    add("theory_beam_unmatched",
        set_beam_lines(drop_exp_beam(wexpel, 3), kav=[1, 2, 0, 3], nbea=[1, 2, 3],
                       ibp=[1, 2, 3], mittel=[1, 1, 1, 0]),
        param="inputs/PARAM_unmatched")
    add("weights_all", set_namelist_value(wexpel, "WR", "0.1,0.2,0.3,0.4,"))
    (gap,), gap_inputs = _source_runs("rfactor", "ag100_rfactor_gaps", ["explicit_gap"])
    add("smooth_gap", set_namelist_value(
        set_wexpel_title(gap_inputs[gap["files"]["WEXPEL"]], title), "ISMOTH", 2))
    add("all_fractional", set_beam_lines(wexpel, kav=[1, 2, 3, 4], nbea=[1, 2, 3, 4],
                                         ibp=[1, 2, 3, 4], mittel=[2, 2, 2, 2]))
    add("ipr2", set_namelist_value(set_namelist_value(wexpel, "IPR", 2), "ISMOTH", 1))
    # one gap in each of three beams, so that the gap offsets of the later
    # beams are those of the earlier ones; MNGAP = number of gaps
    inputs["inputs/PARAM_gaps3"] = set_param(inputs[pendry["param"]], "MNGAP", len(GAPS3))
    add("gaps3_all", add_gaps(set_namelist_value(wexpel, "WR", "0.1,0.2,0.3,0.4,"), GAPS3),
        param="inputs/PARAM_gaps3")
    add("autogap_edge_r2_v0std",
        drop_exp_points(set_namelist_value(wexpel, "WR", "1.,0.,0.,"), *EDGE_HOLE),
        poison_heap=True)
    low_gap = set_namelist_value(wexpel, "WR", "1.,0.,0.,")
    for key, v0 in zip(("V01", "V02"), LOW_GAP_V0):
        low_gap = set_namelist_value(low_gap, key, f"{v0:.2f}")
    add("gap_below_theory_r2", add_gaps(low_gap, [LOW_GAP]),
        no_debug=["The gap integrals of experiment and theory read theory intensities",
                  "at grid indices below 1 (out of bounds) for the larger inner-potential",
                  "shifts, as legacy rfactor.f does with AT(IBT,NG1+NEET), where the",
                  "legacy build reads zeros; these reads count as zero. The one data",
                  "point before the gap is also interpolated as a two-point curve, as",
                  "legacy XNTRP2 does (reads YS(2), out of bounds)."],
        poison_heap=True)
    high_gap = set_namelist_value(wexpel, "WR", "0.1,0.2,0.3,0.4,")
    for key, v0 in zip(("V01", "V02"), HIGH_GAP_V0):
        high_gap = set_namelist_value(high_gap, key, f"{v0:.2f}")
    add("gap_edge_hi_allw", add_gaps(high_gap, [HIGH_GAP]), poison_heap=True)
    no_overlap = set_namelist_value(wexpel, "EMAX", f"{NO_OVERLAP_EMAX:.2f}")
    for beam in range(1, 5):
        no_overlap = drop_exp_points(no_overlap, beam, 0.0, NO_OVERLAP_EXP_MIN - 0.01)
    add("no_overlap", no_overlap,
        no_debug=["No beam overlaps the theory: the averages over the beams divide",
                  "0 by 0, as in legacy rfactor.f (the legacy debug build stops with",
                  "SIGFPE there)."])
    return {"program": "rfactor", "name": case,
            "description": "Ag(100) rfactor input parameters, one run each: a theory beam "
                           "without experimental counterpart (KAV = 0), all four R-factor "
                           "weights, smoothing with a gap, only fractional beams "
                           "(MITTEL = 2), IPR = 2 with smoothing, gaps in three beams, a "
                           "detected gap at the upper end of the data, a gap below the start "
                           "of the theory, a gap near the end of the theory with all weights, no "
                           "overlap with the theory; CPC 2001 Table 7 "
                           "covers KAV, ISMOTH, MITTEL and IPR, not KAV = 0, WR or gaps",
            "covers": ["kav-zero", "weights-all", "smoothing-gap", "fractional-only",
                       "ipr-2", "gaps-several-beams", "gap-outside-common-interval",
                       "gap-below-theory", "gap-near-theory-end", "no-overlap"],
            "runs": runs, "inputs": inputs}


def fortran_e14_5(x):
    """x in Fortran format E14.5 (0.ddddd E+xx)."""
    if x == 0:
        return "   0.00000E+00"
    mantissa, exponent = f"{x:.4E}".split("E")
    digits = mantissa.lstrip("-").replace(".", "")
    return f"{'-' if x < 0 else ''}0.{digits}E{int(exponent) + 1:+03d}".rjust(14)


def beam_scale(k, energy):
    """Intensity factor of the copies of a beam: smooth in energy, different for
    each beam k (0-based), so that every copy has its own R factor."""
    return 1.0 + 0.25 * math.sin(0.037 * energy * (k + 1) + k)


def copy_exp_beams(wexpel, n):
    """Experimental data of an rfactor WEXPEL with n beams: beam k (0-based) is
    beam k mod 4 of the source with its intensities scaled by beam_scale."""
    lines = wexpel.split("\n")
    i_data = lines.index(" (12E14.5)") + 1
    source, i = [], i_data
    while i < len(lines) and lines[i].strip():
        npts, rest = int(lines[i + 1][:4]), lines[i + 1][4:]
        values, i = [], i + 2
        while len(values) < 2 * npts:
            values += [float(v) for v in lines[i].split()]
            i += 1
        source.append((rest, values[0::2], values[1::2]))
    data = []
    for k in range(n):
        rest, energies, intensities = source[k % len(source)]
        data += [f"*({k:3d}  0)*", f"{len(energies):4d}{rest}"]
        pairs = [v for e, a in zip(energies, intensities) for v in (e, a * beam_scale(k, e))]
        data += ["".join(fortran_e14_5(v) for v in pairs[j:j + 12])
                 for j in range(0, len(pairs), 12)]
    return "\n".join(lines[:i_data] + data + lines[i:])


def copy_theory_beams(spectrum, n):
    """A theory spectrum of 4 beams with n beams, as copy_exp_beams: beam k is
    beam k mod 4 scaled by beam_scale. Each record is the energy and geometry
    (F7.2, F7.4) with the first 4 intensities, then 5 intensities per line."""
    lines = spectrum.split("\n")
    if int(lines[1]) != 4:
        raise ValueError("the source spectrum must have 4 beams")
    out = [lines[0], f"{n:5d}"]
    out += [f"{k + 1:5d}{1.0 + k:10.5f}{0.0:10.5f}{1:3d}" for k in range(n)]
    for line in lines[6:]:
        if not line.strip():
            continue
        energy = float(line[:7])
        source = [float(line[14 + 14 * j:28 + 14 * j]) for j in range(4)]
        values = [fortran_e14_5(source[k % 4] * beam_scale(k, energy)) for k in range(n)]
        out.append(line[:14] + "".join(values[:4]))
        out += ["".join(values[j:j + 5]) for j in range(4, n, 5)]
        if (n - 4) % 5 == 0:
            out.append("")
    return "\n".join(out) + "\n"


def ag100_rfactor_beams():
    """Spec of ag100_rfactor_beams: more experimental and theory beams than one
    25I4 record of KAV/NBEA/IBP/MITTEL holds."""
    case = "ag100_rfactor_beams"
    (pendry,), inputs = _source_runs("rfactor", "ag100_rfactor", ["pendry"])
    title = f"tl-test {case}"
    spectrum = set_first_line_title(inputs.pop(pendry["stdin"]), title)
    wexpel = set_wexpel_title(inputs.pop(pendry["files"]["WEXPEL"]), title)
    param = inputs.pop(pendry["param"])
    runs = []
    for n in N_BEAMS:
        name = f"nb{n}"
        inputs[f"inputs/spectrum_{name}"] = copy_theory_beams(spectrum, n)
        inputs[f"inputs/PARAM_{name}"] = set_param(set_param(param, "MNBED", n), "MNBTD", n)
        beams = list(range(1, n + 1))
        inputs[f"inputs/WEXPEL_{name}"] = set_beam_lines(
            copy_exp_beams(wexpel, n), kav=beams, nbea=beams, ibp=beams,
            mittel=[1 + k % 2 for k in range(n)])
        runs.append({**pendry, "name": name, "stdin": f"inputs/spectrum_{name}",
                     "param": f"inputs/PARAM_{name}",
                     "files": {"WEXPEL": f"inputs/WEXPEL_{name}"}})
    return {"program": "rfactor", "name": case,
            "description": "Ag(100) rfactor as ag100_rfactor pendry with 25 and 30 "
                           "experimental and theory beams (scaled copies of the 4 beams): "
                           "beam lists of one full 25I4 record and of two records",
            "covers": ["beams-25", "beams-over-25"],
            "runs": runs, "inputs": inputs}


CASES = {f.__name__: f for f in (ir100_2x1_o_unsorted, ag100_psformat2, ag100_rfactor_hand,
                                 ag100_rfactor_params, ag100_rfactor_beams)}


def _write(spec, cases_dir):
    write_case(cases_dir / spec["program"] / spec["name"], HEADER, spec["program"],
               spec["description"], spec["covers"], spec["runs"], spec["inputs"])


def _differences(spec, derived_root):
    rel = Path(spec["program"]) / spec["name"]
    derived, committed = derived_root / rel, CASES_DIR / rel
    files = sorted({p.relative_to(derived) for p in derived.rglob("*") if p.is_file()}
                   | {p.relative_to(committed) for p in (committed / "inputs").rglob("*")
                      if p.is_file()} | {Path("case.toml")})
    return [str(rel / f) for f in files
            if not (derived / f).is_file() or not (committed / f).is_file()
            or (derived / f).read_bytes() != (committed / f).read_bytes()]


def main(argv=None):
    """Write the hand-made cases, or with --check compare them with the committed ones."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", nargs="+", choices=sorted(CASES))
    parser.add_argument("--check", action="store_true",
                        help="derive into a temporary directory and compare with the "
                             "committed cases")
    args = parser.parse_args(argv)
    specs = [CASES[name]() for name in (args.only or CASES)]
    if not args.check:
        for spec in specs:
            _write(spec, CASES_DIR)
            print(f"{spec['program']}/{spec['name']}: {len(spec['runs'])} runs")
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        differ = []
        for spec in specs:
            _write(spec, Path(tmp))
            differ += _differences(spec, Path(tmp))
    if differ:
        print("differs from the derived inputs:\n  " + "\n  ".join(differ))
        return 1
    print(f"{len(specs)} cases match their derivation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
