"""Tier 2: run viperleed calc on a frozen system and read R factors/error curves.

A system directory (tests/tier2/systems/<name>/) contains viperleed inputs
(PARAMETERS with RUN = 1 5, POSCAR, VIBROCC, PHASESHIFTS, IVBEAMS,
EXPBEAMS.csv, DISPLACEMENTS) and expected.toml with the baseline results.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

from .cases import case_dirs

INPUT_FILES = ("PARAMETERS", "POSCAR", "VIBROCC", "PHASESHIFTS", "IVBEAMS",
               "EXPBEAMS.csv", "DISPLACEMENTS")
TIMEOUT_S = 3600
RELEASE_TAG = re.compile(r"v\d+\.\d+\.\d+")


class ViperleedFailed(RuntimeError):
    """A viperleed run failed, or its results are missing or malformed."""


def load_systems(systems_dir):
    """System directories in systems_dir, sorted by name.

    CaseError if there are none or one has no expected.toml (a broken system
    must fail, not disappear from the parametrization).
    """
    return case_dirs(systems_dir, "expected.toml", "system")


def set_parameters(params_file, values):
    """Set 'KEY = value' lines in PARAMETERS (replace existing, else append)."""
    text = Path(params_file).read_text()
    for key, value in values.items():
        line = f"{key} = {value}"
        pattern = rf"(?m)^{re.escape(key)}\s*=.*$"
        if re.search(pattern, text):
            text = re.sub(pattern, line.replace("\\", "\\\\"), text)
        else:
            text = text.rstrip("\n") + f"\n{line}\n"
    Path(params_file).write_text(text)


def viperleed_version(viperleed_cmd):
    """Version banner printed by 'viperleed calc --version'."""
    return subprocess.run([*viperleed_cmd, "calc", "--version"],
                          capture_output=True, text=True).stdout.strip()


def check_viperleed_version(viperleed_cmd, pinned):
    """Raise ViperleedFailed unless viperleed_cmd runs viperleed version pinned."""
    banner = viperleed_version(viperleed_cmd)
    m = re.search(r"\bv(\d\S*)$", banner)
    if not m or m.group(1) != pinned:
        installed = m.group(1) if m else repr(banner)
        raise ViperleedFailed(f"tier 2 needs viperleed {pinned} (tests/tolerances.toml "
                              f"[tier2] viperleed_version); "
                              f"'{' '.join(viperleed_cmd)}' runs {installed}")


def is_released_tl_tree(tl_root):
    """True if tl_root is a git checkout of a released TL version, unmodified.

    Released: tl_root is the top level of its own git checkout (a worktree
    counts), HEAD carries a tag 'vMAJOR.MINOR.PATCH' and no tracked file is
    modified. Only such a tree can pass viperleed's checksum check; anything
    else cannot: branch tips, modified trees, trees without tags (e.g. a
    shallow or --no-tags clone) and trees without a git checkout of their own,
    including a folder inside another checkout. Repository variables such as
    GIT_DIR in the environment are ignored.
    """
    local_vars = subprocess.run(["git", "rev-parse", "--local-env-vars"],
                                capture_output=True, text=True).stdout.split()
    env = {k: v for k, v in os.environ.items() if k not in local_vars}

    def git(*args):
        return subprocess.run(["git", "-C", str(tl_root), *args],
                              capture_output=True, text=True, env=env)
    top = git("rev-parse", "--show-toplevel")
    if top.returncode or Path(top.stdout.strip()).resolve() != Path(tl_root).resolve():
        return False
    tags = git("tag", "--points-at", "HEAD")
    status = git("status", "--porcelain", "--untracked-files=no")
    if tags.returncode or status.returncode:
        return False
    released = any(RELEASE_TAG.fullmatch(t) for t in tags.stdout.split())
    return released and not status.stdout.strip()


def run_viperleed_calc(system_dir, tl_root, viperleed_cmd, workdir, parameters=None):
    """Copy the inputs to workdir, run viperleed calc, return parse_results().

    TL_IGNORE_CHECKSUM is set from is_released_tl_tree(tl_root), overriding
    the system's PARAMETERS.
    parameters: optional {KEY: value} set in PARAMETERS (e.g. FORTRAN_COMP);
    ValueError if it sets TL_IGNORE_CHECKSUM.
    Raises ViperleedFailed (with the log tail) on a non-zero exit code or an
    '### ERROR ###' in the output.
    """
    parameters = dict(parameters or {})
    if "TL_IGNORE_CHECKSUM" in parameters:
        raise ValueError("TL_IGNORE_CHECKSUM is set by the harness, not by the caller")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    for name in INPUT_FILES:
        src = Path(system_dir) / name
        if src.is_file():
            shutil.copy(src, workdir / name)
    ignore_checksum = not is_released_tl_tree(tl_root)
    set_parameters(workdir / "PARAMETERS",
                   {"TL_IGNORE_CHECKSUM": str(ignore_checksum), **parameters})
    log = workdir / "viperleed.out"
    with open(log, "w") as out:
        proc = subprocess.run([*viperleed_cmd, "calc", "-t", str(tl_root)], cwd=workdir,
                              stdout=out, stderr=subprocess.STDOUT, timeout=TIMEOUT_S)
    text = log.read_text(errors="replace")
    if proc.returncode != 0 or "### ERROR ###" in text:
        raise ViperleedFailed(f"viperleed failed (exit {proc.returncode}) in {workdir}; "
                              f"log tail:\n{text[-4000:]}")
    return parse_results(workdir)


def parse_rout(text):
    """(V0 shift, R) of the best result in an R_OUT/ROUT file."""
    m = re.search(r"INNER POTENTIAL SHIFT\s+(-?[\d.]+) EV, AVERAGE R-FACTOR =\s*([\d.]+)", text)
    if not m:
        raise ViperleedFailed("no 'AVERAGE R-FACTOR' summary in R_OUT file")
    return float(m.group(1)), float(m.group(2))


def parse_error_csv(text):
    """{column name: [values]} of an Errors_*.csv file."""
    lines = [l for l in text.splitlines() if l.strip()]
    header = [h.strip() for h in lines[0].split(",")]
    cols = {h: [] for h in header}
    for n, line in enumerate(lines[1:], start=2):
        fields = line.split(",")
        if len(fields) != len(header):
            raise ViperleedFailed(f"Errors csv line {n} has {len(fields)} fields, "
                                  f"header has {len(header)}")
        for h, v in zip(header, fields):
            cols[h].append(float(v))
    return cols


def parse_results(workdir):
    """Read the best R factor, V0 shift and the error curves from the OUT of a calc run."""
    out = Path(workdir) / "OUT"
    routs = sorted(out.glob("R_OUT_refcalc_R=*"))
    if len(routs) != 1:
        raise ViperleedFailed(f"expected one R_OUT_refcalc file in {out}, found {routs}")
    v0, r = parse_rout(routs[0].read_text(errors="replace"))
    errors_zip = out / "Errors.zip"
    if not errors_zip.is_file():
        raise ViperleedFailed(f"no Errors.zip in {out}")
    with zipfile.ZipFile(errors_zip) as z:
        errors = {n: parse_error_csv(z.read(n).decode()) for n in sorted(z.namelist())
                  if n.endswith(".csv")}
    if not errors:
        raise ViperleedFailed("Errors.zip contains no error curves")
    return {"refcalc": {"r": r, "v0_shift": v0}, "errors": errors}


def to_toml(results, comment=""):
    """Serialize results (from parse_results) to expected.toml text, after the comment lines."""
    lines = [f"# {l}" for l in comment.splitlines()]
    lines += ["[refcalc]", f"r = {results['refcalc']['r']!r}",
              f"v0_shift = {results['refcalc']['v0_shift']!r}", ""]
    for name, cols in results["errors"].items():
        lines += ["[[errors]]", f'file = "{name}"']
        lines += ["columns = [" + ", ".join(f'"{c}"' for c in cols) + "]"]
        for i, (col, vals) in enumerate(cols.items()):
            lines.append(f"col{i} = [" + ", ".join(repr(v) for v in vals) + "]")
        lines.append("")
    return "\n".join(lines)


def from_toml(data):
    """Inverse of to_toml: results from the parsed expected.toml."""
    errors = {}
    for e in data["errors"]:
        errors[e["file"]] = {c: e[f"col{i}"] for i, c in enumerate(e["columns"])}
    return {"refcalc": dict(data["refcalc"]), "errors": errors}


def compare_results(actual, expected, r_atol, other_atol):
    """List of differences; R values within r_atol, other numbers within other_atol.

    NaN on either side is a difference.
    """
    fails = []
    for key, atol in (("r", r_atol), ("v0_shift", other_atol)):
        a, e = actual["refcalc"][key], expected["refcalc"][key]
        if not abs(a - e) <= atol:
            fails.append(f"refcalc {key}: {a} vs {e}")
    if set(actual["errors"]) != set(expected["errors"]):
        fails.append(f"error curves differ: {sorted(actual['errors'])} vs "
                     f"{sorted(expected['errors'])}")
    for name in sorted(set(actual["errors"]) & set(expected["errors"])):
        acols, ecols = actual["errors"][name], expected["errors"][name]
        if list(acols) != list(ecols):
            fails.append(f"{name}: columns {list(acols)} vs {list(ecols)}")
            continue
        for col in ecols:
            atol = r_atol if col == "R" else other_atol
            if len(acols[col]) != len(ecols[col]):
                fails.append(f"{name}/{col}: {len(acols[col])} vs {len(ecols[col])} points")
                continue
            for i, (a, e) in enumerate(zip(acols[col], ecols[col])):
                if not abs(a - e) <= atol:
                    fails.append(f"{name}/{col}[{i}]: {a} vs {e}")
    return fails
