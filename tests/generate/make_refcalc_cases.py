"""Generate frozen ref-calc test cases with viperleed.

For each system in systems.SYSTEMS:
1. unpack the example inputs into <workroot>/<case>/calc, apply the
   PARAMETERS overrides, run viperleed with RUN = 0 (initialization: annotated
   POSCAR, PHASESHIFTS, IVBEAMS, VIBROCC), then RUN = 1 with -k -w work;
2. take the complete refcalc-FIN, the per-LMAX PARAM blocks (refcalc-PARAM)
   and viperleed's energy -> LMAX assignment (refcalc-*.log) from the work dir;
3. check that the FIN exhibits the code paths listed in must_cover; if a
   check fails, nothing is written;
4. write tests/cases/refcalc/<case>/{case.toml, inputs/FIN_<run>,
   inputs/PARAM_L<n>}: three single-energy runs (lowest energy at the lowest
   LMAX, a middle energy at a middle LMAX, the highest energy at the highest
   LMAX) and one run over 4 consecutive energies at one LMAX. The per-run FIN
   is edited exactly as viperleed does it (edit_fin_energy_lmax) and its
   title line is replaced by 'tl-test <case>'.

viperleed runs use the baseline tree itself (configured baseline_tl_root,
see tests/harness/paths.py), which must be at BASELINE_COMMIT with an
unmodified TensErLEED/ and contain eeasisss (built there with
make gcc GCC_FFLAGS="-O3 -fno-finite-math-only -std=legacy"; see
ensure_generator_tl). viperleed compiles in its own work directories and
does not modify the tree.

Usage (in an environment with viperleed installed):
    python tests/generate/make_refcalc_cases.py \
        [--only CASE ...] [--workroot DIR] [--check-only] [--benchmark]
"""
import argparse
import math
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import fortranformat as ff

sys.path.insert(0, str(Path(__file__).parent))
import case_files  # noqa: E402
from baseline_tree import BASELINE_COMMIT, GeneratorError, check_baseline_tree  # noqa: E402
from systems import BENCHMARKS, SYSTEMS  # noqa: E402

from harness import tier2  # noqa: E402
from harness.paths import location  # noqa: E402

TESTS_DIR = Path(__file__).parents[1]
CASES_DIR = TESTS_DIR / "cases" / "refcalc"
BENCHMARK_DIR = TESTS_DIR / "cases" / "benchmark"
N_MULTI = 4  # energies in the multi-energy run


# ---------------------------------------------------------------- TL tree ---

def ensure_generator_tl():
    """Return the baseline tree (configured baseline_tl_root), after checking
    that it is at BASELINE_COMMIT, TensErLEED/ is unmodified and eeasisss
    exists."""
    tl = location("baseline_tl_root")
    check_baseline_tree(tl)
    if not (tl / "eeasisss").is_file():
        # gfortran >= 14 rejects a legacy FORMAT in eeasisss at run time
        # unless built with -std=legacy
        raise GeneratorError(f"{tl}/eeasisss missing; build it there with "
                             "make gcc GCC_FFLAGS=\"-O3 -fno-finite-math-only -std=legacy\"")
    return tl


# --------------------------------------------------------------- viperleed ---

def _viperleed_exe():
    exe = Path(sys.executable).parent / "viperleed"
    if not exe.is_file():
        raise GeneratorError(f"viperleed not found next to {sys.executable}; "
                             "run this script in an environment with viperleed installed")
    return exe


def run_viperleed(calc_dir, tl, *args, allow_section_errors=()):
    log = calc_dir / f"viperleed-{'-'.join(a.strip('-') for a in args) or 'run'}.out"
    with open(log, "w") as out:
        proc = subprocess.run([str(_viperleed_exe()), "calc", "-t", str(tl), *args],
                              cwd=calc_dir, stdout=out, stderr=subprocess.STDOUT)
    text = log.read_text(errors="replace")
    failed_sections = set(re.findall(r"# Error in section (.+?)\s*$", text, re.M))
    if failed_sections and failed_sections <= set(allow_section_errors):
        return  # tolerated (benchmark: only the refcalc inputs are needed)
    if proc.returncode != 0 or "### ERROR ###" in text:
        raise GeneratorError(f"viperleed failed in {calc_dir} (see {log}):\n"
                             + text[-3000:])


def _param_key(line):
    s = line.split("!")[0]
    if "=" not in s:
        return None
    return " ".join(s.split("=")[0].split()).upper()


def edit_parameters(text, overrides):
    """Comment out lines whose key is overridden, append the overrides."""
    keys = {" ".join(k.split()).upper() for k in overrides}
    lines = []
    for line in text.splitlines():
        if _param_key(line) in keys:
            line = "! " + line + "   ! overridden by tl-test generator"
        lines.append(line)
    lines.append("")
    lines.append("! ---- tl-test generator overrides ----")
    lines += [f"{k} = {v}" for k, v in overrides.items() if v is not None]
    return "\n".join(lines) + "\n"


def set_run(params_file, run):
    text = params_file.read_text()
    text, n = re.subn(r"(?m)^RUN\s*=.*$", f"RUN = {run}", text)
    if n != 1:
        raise GeneratorError(f"expected exactly one RUN line in {params_file}")
    params_file.write_text(text)


def prepare_calc_dir(sysdef, calc_dir):
    if calc_dir.exists():
        shutil.rmtree(calc_dir)
    calc_dir.mkdir(parents=True)
    if sysdef.source_zip is None:
        for name, content in sysdef.base_files.items():
            (calc_dir / name).write_text(content)
    else:
        with zipfile.ZipFile(sysdef.source_zip) as z:
            for m in z.namelist():
                name = Path(m)
                if (m.startswith("__MACOSX") or m.endswith("/")
                        or name.name.startswith(".")):
                    continue
                (calc_dir / name.name).write_bytes(z.read(m))
    overrides = {**sysdef.parameters, "RUN": "0",
                 "TL_IGNORE_CHECKSUM": "True",
                 # N_CORES > 1: viperleed then compiles one executable per
                 # LMAX (energy-dependent LMAX); with 1 core only LMAX max
                 "N_CORES": "2"}
    params = calc_dir / "PARAMETERS"
    params.write_text(edit_parameters(params.read_text(), overrides))
    for name, content in sysdef.extra_files.items():
        (calc_dir / name).write_text(content)


# ------------------------------------------------------- FIN / PARAM / log ---

def edit_fin_energy_lmax(fin, energies, step, lmax, title):
    """viperleed's edit_fin_energy_lmax (TL >= 1.7), plus fixed title.

    Single energy: EI = E, EF = E + 0.01, DE = 1 (as viperleed). Several
    energies: EI = first, EF = last + 0.01, DE = step.
    """
    _, _, rest = fin.split("\n", maxsplit=2)
    eformatter = ff.FortranRecordWriter("3F9.2")
    if len(energies) == 1:
        evals = [energies[0], energies[0] + 0.01, 1.0]
    else:
        evals = [energies[0], energies[-1] + 0.01, step]
    energy = eformatter.write(evals).ljust(30) + "EI,EF,DE"
    before_lmax, after_lmax = rest.split("   LMAX", maxsplit=1)
    before_lmax, *_ = before_lmax.rsplit("\n", maxsplit=1)
    after_lmax = str(lmax).rjust(3).ljust(45) + "LMAX" + after_lmax
    return "\n".join((title, energy, before_lmax, after_lmax))


def split_param(text):
    """refcalc-PARAM -> {lmax: PARAM text}."""
    blocks = re.split(r"(?m)^### PARAM file for LMAX = (\d+) ###\s*$", text)
    out = {}
    for lmax, body in zip(blocks[1::2], blocks[2::2]):
        out[int(lmax)] = body.strip("\n") + "\n"
    if not out:
        raise GeneratorError("no PARAM blocks found in refcalc-PARAM")
    return out


def energy_lmax_map(work):
    """{energy: lmax} from viperleed's refcalc log in the work dir."""
    emap = {}
    for log in work.glob("refcalc-*.log"):
        parts = re.split(r"### STARTING LOG FOR refcalc-part_([\d.]+)eV ###",
                         log.read_text(errors="replace"))
        for en, body in zip(parts[1::2], parts[2::2]):
            m = re.search(r"LMAX\s*=\s*(\d+)", body)
            if not m:
                raise GeneratorError(f"no LMAX in log part for {en} eV")
            emap[float(en)] = int(m.group(1))
    if not emap:
        raise GeneratorError(f"no per-energy refcalc log found in {work}")
    return dict(sorted(emap.items()))


def select_runs(emap):
    """Three single-energy runs + one 4-energy run (see module docstring)."""
    lvals = sorted(set(emap.values()))
    lmin, lmid, lmax = lvals[0], lvals[len(lvals) // 2], lvals[-1]
    by_l = {l: [e for e, v in emap.items() if v == l] for l in lvals}
    mid_es = by_l[lmid]
    singles = [(min(by_l[lmin]), lmin), (mid_es[len(mid_es) // 2], lmid),
               (max(by_l[lmax]), lmax)]
    singles = list(dict.fromkeys(singles))  # if only one LMAX value
    energies = list(emap)
    windows = [energies[i:i + N_MULTI] for i in range(len(energies) - N_MULTI + 1)
               if len({emap[e] for e in energies[i:i + N_MULTI]}) == 1]
    if not windows:
        raise GeneratorError("no 4 consecutive energies with the same LMAX")
    preferred = [w for w in windows if emap[w[0]] == lmid] or windows
    multi = preferred[len(preferred) // 2]
    return singles, (multi, emap[multi[0]])


# ---------------------------------------------------------- coverage checks ---

def _fin_lines(fin, marker):
    return [l for l in fin.splitlines() if marker in l]


def _ints(lines):
    return [int(l.split()[0]) for l in lines]


def _sites(fin):
    """[[(element, occ), ...] per site] from the site block."""
    sites, cur = [], None
    for line in fin.splitlines():
        if re.match(r"^-\s+site type \d+", line):
            cur = []
            sites.append(cur)
        m = re.search(r"Occ & VibAmp for (\w+) in", line)
        if m and cur is not None:
            cur.append((m.group(1), float(line.split()[0])))
    return sites


def _layer_sites(fin):
    """{layer type: [site type of each Bravais sublayer]}."""
    layers, cur = {}, None
    for line in fin.splitlines():
        m = re.match(r"^-\s+layer type (\d+)", line)
        if m:
            cur = layers.setdefault(int(m.group(1)), [])
        elif cur is not None and re.search(r"Atom N=\d+", line):
            cur.append(int(line.split()[0]))
        elif "define bulk stacking" in line:
            cur = None
    return layers


def _n_beam_sets(fin):
    """Number of beam-set headers (lone integer lines) before the TST line."""
    lines = fin.splitlines()
    tst = next(i for i, l in enumerate(lines) if "TST" in l)
    return sum(1 for l in lines[6:tst] if re.fullmatch(r"\s*\d+\s*", l))


def check_coverage(fin, checks):
    nsub = _ints(_fin_lines(fin, "number of Bravais sublayers"))
    results = {}
    for check in checks:
        if check == "nsub_all_1":
            ok = bool(nsub) and all(n == 1 for n in nsub)
        elif check == "nsub_gt_1":
            ok = any(n > 1 for n in nsub)
        elif check == "nsub_ge_4":
            ok = any(n >= 4 for n in nsub)
        elif check == "bulk_composite":
            lay = _ints(_fin_lines(fin, "LAY = "))
            ok = any(l == 2 and n > 1 for l, n in zip(lay, nsub))
        elif check.startswith("tensor_layers="):
            want = int(check.split("=")[1])
            flags = _ints(_fin_lines(fin, "Tensor output is required"))
            ok = sum(flags) == want
        elif check == "theta_nonzero":
            (line,) = _fin_lines(fin, "THETA FI")
            ok = float(line.split()[0]) != 0.0
        elif check == "bulk_two_layers":
            (top,) = _fin_lines(fin, "top layer of bulk unit")
            (bot,) = _fin_lines(fin, "bottom layer of bulk unit")
            ok = top.split()[0] != bot.split()[0]
        elif check == "bulk_layers_distinct_sites":
            (top,) = _fin_lines(fin, "top layer of bulk unit")
            (bot,) = _fin_lines(fin, "bottom layer of bulk unit")
            sites = _layer_sites(fin)
            ok = sites[int(top.split()[0])] != sites[int(bot.split()[0])]
        elif check == "two_elements":
            ok = len({el for site in _sites(fin) for el, _ in site}) >= 2
        elif check == "mixed_site":
            ok = any(sum(occ > 0 for _, occ in site) >= 2 for site in _sites(fin))
        elif check == "hexagonal":
            (a1,) = _fin_lines(fin, "ARA1")
            (a2,) = _fin_lines(fin, "ARA2")
            v1, v2 = ([float(x) for x in l.split()[:2]] for l in (a1, a2))
            cos = (v1[0] * v2[0] + v1[1] * v2[1]) / (math.hypot(*v1) * math.hypot(*v2))
            ok = abs(abs(cos) - 0.5) < 1e-4
        elif check == "beam_sets_gt_1":
            ok = _n_beam_sets(fin) > 1
        else:
            raise GeneratorError(f"unknown coverage check {check!r}")
        results[check] = ok
    return results


# -------------------------------------------------------------- case output ---

def _run_name(energies, lmax):
    if len(energies) == 1:
        return f"E{energies[0]:03.0f}_L{lmax:02d}"
    return f"E{energies[0]:03.0f}-{energies[-1]:03.0f}_L{lmax:02d}"


def write_case(sysdef, fin, params, runs, step, provenance, cases_dir=CASES_DIR):
    """runs: [([energies], lmax)]."""
    title = f"tl-test {sysdef.name}"
    inputs, toml_runs = {}, []
    for energies, lmax in runs:
        name = _run_name(energies, lmax)
        fin_path, param_path = f"inputs/FIN_{name}", f"inputs/PARAM_L{lmax:02d}"
        inputs[fin_path] = edit_fin_energy_lmax(fin, energies, step, lmax, title)
        inputs[param_path] = params[lmax]
        toml_runs.append({"name": name, "stdin": fin_path, "param": param_path,
                          "outputs": ["fd.out", "amp.out", "T_*"]})
    case_files.write_case(cases_dir / sysdef.name,
                          f"# Generated by tests/generate/make_refcalc_cases.py; {provenance}",
                          "refcalc", sysdef.description, sysdef.covers, toml_runs, inputs)
    return [run["name"] for run in toml_runs]


def _theo_step(sysdef):
    return float(sysdef.parameters["THEO_ENERGIES"].split()[2])


def generate(sysdef, workroot, tl, benchmark=False):
    calc = workroot / sysdef.name / "calc"
    prepare_calc_dir(sysdef, calc)
    run_viperleed(calc, tl)                      # RUN = 0: initialization
    set_run(calc / "PARAMETERS", "1")
    # RUN = 1: reference calc. Benchmarks use a sparse energy grid on which
    # viperleed's R-factor step fails; only the refcalc inputs are needed.
    run_viperleed(calc, tl, "-k", "-w", "work",
                  allow_section_errors=("R-FACTOR CALCULATION",) if benchmark else ())
    work = calc / "work"
    fin = (work / "refcalc-FIN").read_text()
    results = check_coverage(fin, sysdef.must_cover)
    if not all(results.values()):
        return results, [], None
    params = split_param((work / "refcalc-PARAM").read_text())
    emap = energy_lmax_map(work)
    if benchmark:  # lowest and highest energy, each at viperleed's LMAX
        e_lo, e_hi = min(emap), max(emap)
        runs = [([e_lo], emap[e_lo]), ([e_hi], emap[e_hi])]
    else:
        singles, multi = select_runs(emap)
        runs = [([e], l) for e, l in singles] + [multi]
    version = tier2.viperleed_version([str(_viperleed_exe())])
    names = write_case(sysdef, fin, params, runs, _theo_step(sysdef),
                       f"viperleed {version}; TL {BASELINE_COMMIT}",
                       BENCHMARK_DIR if benchmark else CASES_DIR)
    return results, names, emap


def check_committed(names):
    failed = False
    for name in names:
        fins = sorted((CASES_DIR / name / "inputs").glob("FIN_*"))
        if not fins:
            print(f"{name}: no FIN files")
            failed = True
        for fin in fins:
            for check, ok in check_coverage(fin.read_text(),
                                            SYSTEMS[name].must_cover).items():
                print(f"{name}/{fin.name} coverage {check}: {'OK' if ok else 'FAILED'}")
                failed |= not ok
    return 1 if failed else 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", nargs="+", choices=sorted(SYSTEMS) + sorted(BENCHMARKS))
    parser.add_argument("--benchmark", action="store_true",
                        help="generate the benchmark cases (tests/cases/benchmark) instead")
    parser.add_argument("--workroot", type=Path, default=None,
                        help="default: generate_workroot (tests/harness/paths.py)")
    parser.add_argument("--check-only", action="store_true",
                        help="only re-run the coverage checks on the committed FINs")
    args = parser.parse_args(argv)
    if args.check_only:
        return check_committed(args.only or SYSTEMS)
    args.workroot = args.workroot or location("generate_workroot")
    tl = ensure_generator_tl()
    failed = False
    systems = BENCHMARKS if args.benchmark else SYSTEMS
    for name in args.only or systems:
        results, runs, emap = generate(systems[name], args.workroot, tl, args.benchmark)
        if runs:
            lm = sorted(set(emap.values()))
            print(f"{name}: runs {', '.join(runs)}; LMAX values used by viperleed {lm}")
        else:
            print(f"{name}: coverage check failed, case not written")
        for check, ok in results.items():
            print(f"  coverage {check}: {'OK' if ok else 'FAILED'}")
            failed |= not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
