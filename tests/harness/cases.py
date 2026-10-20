"""Test-case definitions (cases/<kind>/<case>/case.toml) and tolerances.

case.toml:
    program = "refcalc" | "rfactor"
    description = "..."       # informative only, not loaded
    covers = ["code path", ...]  # informative only, not loaded
    [[runs]]
    name = "..."              # expected outputs live in expected/<name>/
    stdin = "inputs/..."      # file fed to the program on stdin
    param = "inputs/PARAM"    # PARAM file, used only in legacy build mode
    files = { WEXPEL = "inputs/..." }  # optional; copied into the run dir
    outputs = ["fd.out", "T_*"]        # globs, resolved in expected/<name>/
All paths are relative to the case directory.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import fnmatch
import json
import re
import shutil
import tomllib
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path

from .build import build_program
from .compare import Tolerance, compare_files
from .runner import ProgramFailed, run_program

RUN_TIMEOUT_S = 900
KINDS = ("refcalc", "rfactor")  # case kinds run in tier 1

TESTS_DIR = Path(__file__).parents[1]
CASES_DIR = TESTS_DIR / "cases"
TOLERANCES_FILE = TESTS_DIR / "tolerances.toml"


@dataclass
class Run:
    """One execution of a program within a case: its inputs and the expected output globs."""
    name: str
    stdin: Path
    param: Path | None
    files: dict = field(default_factory=dict)  # name in run dir -> source Path
    outputs: list = field(default_factory=list)


@dataclass
class Case:
    """A test case: a program and the runs made with it."""
    name: str
    path: Path
    program: str
    runs: list


class CaseError(ValueError):
    """A case or system directory is missing or malformed."""


def _load_case(case_dir):
    data = tomllib.loads((case_dir / "case.toml").read_text())
    if not data.get("runs"):
        raise CaseError(f"{case_dir / 'case.toml'} defines no runs")
    runs = []
    for r in data["runs"]:
        runs.append(Run(
            name=r["name"],
            stdin=case_dir / r["stdin"],
            param=case_dir / r["param"] if r.get("param") else None,
            files={k: case_dir / v for k, v in r.get("files", {}).items()},
            outputs=list(r["outputs"]),
        ))
    return Case(case_dir.name, case_dir, data["program"], runs)


def load_cases_from(root, kind):
    """All cases in root/kind/*/case.toml, sorted by name.

    CaseError for a missing or empty root/kind, a case directory without
    case.toml or a case without runs (a broken case must fail, not disappear
    from the parametrization).
    """
    return [_load_case(d) for d in case_dirs(Path(root) / kind, "case.toml", "case")]


def case_dirs(parent, marker, what):
    """Directories in parent, sorted by name, each of which must hold marker.

    CaseError for a missing or empty parent or a directory without marker;
    what ("case", "system") names the directories in the message.
    """
    parent = Path(parent)
    dirs = sorted(d for d in parent.iterdir()
                  if d.is_dir() and not d.name.startswith((".", "__"))) if parent.is_dir() else []
    if not dirs:
        raise CaseError(f"no {what}s in {parent}")
    missing = [d.name for d in dirs if not (d / marker).is_file()]
    if missing:
        raise CaseError(f"{what} directories without {marker} in {parent}: {missing}")
    return dirs


def load_cases(kind):
    """All cases of a kind in the committed cases directory (see load_cases_from)."""
    return load_cases_from(CASES_DIR, kind)


def iter_runs(kinds=KINDS, only=None):
    """Yield (case, run) for every run of every case of the given kinds.

    Cases come in load_cases order (kinds as given, then sorted by name),
    runs in the order of their case.toml. If only is given and not empty,
    cases whose name is not in it are skipped.
    """
    for kind in kinds:
        for case in load_cases(kind):
            if only and case.name not in only:
                continue
            for run in case.runs:
                yield case, run


def expected_outputs(case, run):
    """Expected output files of run (globs resolved, sorted by name)."""
    exp_dir = case.path / "expected" / run.name
    found = set()
    for pattern in run.outputs:
        matches = [p for p in exp_dir.glob(pattern) if p.is_file()]
        if not matches:
            raise FileNotFoundError(
                f"{case.name}/{run.name}: no expected output matches {pattern!r} in {exp_dir}")
        found.update(matches)
    return sorted(found, key=lambda p: p.name)


@lru_cache(maxsize=None)
def load_tolerances(path=TOLERANCES_FILE):
    """The parsed tolerances file (cached per path)."""
    return tomllib.loads(Path(path).read_text())


def pattern_for(name, patterns):
    """First of patterns matching name (fnmatch, case-sensitive), or None."""
    return next((p for p in patterns if fnmatch.fnmatchcase(name, p)), None)


def tolerance_for(filename):
    """(Tolerance, ignore_line_patterns, skip_lines) for an output file name.

    The first pattern in [numeric] matching filename (fnmatch) is used
    (keys rtol, atol_frac; optional atol, ulps, scale, file_cap,
    exact_f_decimals, exact_columns, and for scale "tensor" the table
    sections = {section: fraction of the block max});
    KeyError if none matches.
    """
    numeric = load_tolerances()["numeric"]
    pattern = pattern_for(filename, numeric)
    if pattern is None:
        raise KeyError(f"no tolerance configured for output file {filename!r}")
    t = numeric[pattern]
    sections = tuple((k, float(v)) for k, v in t.get("sections", {}).items())
    return (Tolerance(float(t["rtol"]), float(t["atol_frac"]),
                      float(t.get("atol", 0.0)), float(t.get("ulps", 0.0)),
                      t.get("scale", "file"), sections,
                      float(t["file_cap"]) if "file_cap" in t else None,
                      int(t.get("exact_f_decimals", -1)),
                      tuple(t.get("exact_columns", ()))),
            list(t.get("ignore_line_patterns", [])),
            list(t.get("skip_lines", [])))


def compare_outputs(actual_dir, expected_files, adjust=None):
    """{file name: failure messages} for each expected file that does not match.

    Each file in expected_files (paths) is compared with the file of the same
    name in actual_dir, using the tolerance configured for that name
    (tolerance_for); a missing actual file is a failure. adjust, if given,
    maps each Tolerance to the one actually used. An empty dict means all
    files match.
    """
    failures = {}
    for exp in expected_files:
        tol, ignore, skip = tolerance_for(exp.name)
        if adjust is not None:
            tol = adjust(tol)
        res = compare_files(Path(actual_dir) / exp.name, exp, tol, ignore, skip)
        if not res.ok:
            failures[exp.name] = res.first_failures
    return failures


def format_failures(failures):
    """Lines 'name: message; message' for a compare_outputs result."""
    return [f"{name}: " + "; ".join(msgs) for name, msgs in failures.items()]


def file_set_failure(actual_dir, outputs, expected_names):
    """Message naming missing and extra files, or '' if the sets agree.

    The files produced are those in actual_dir matching a glob in outputs;
    expected_names are the file names that must be produced.
    """
    produced = {p.name for pattern in outputs
                for p in Path(actual_dir).glob(pattern) if p.is_file()}
    expected = set(expected_names)
    if produced == expected:
        return ""
    return (f"missing: {sorted(expected - produced)}, "
            f"extra: {sorted(produced - expected)}")


def execute_run(case, run, tl_root, flags, workdir):
    """Build case.program from tl_root, run it in workdir, return RunResult.

    Copies run.files into workdir (under their keys), feeds run.stdin,
    writes the RunResult (without stdout/stderr) to workdir/.runinfo.json,
    and raises ProgramFailed if the program did not terminate cleanly.
    """
    exe = build_program(tl_root, case.program, run.param, flags)
    return run_built(case, run, exe, workdir)


def run_built(case, run, exe, workdir):
    """Run an already built executable; see execute_run."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    for name, src in run.files.items():
        shutil.copy(src, workdir / name)
    res = run_program(exe, workdir, run.stdin, timeout=RUN_TIMEOUT_S)
    info = {k: v for k, v in asdict(res).items() if k not in ("stdout", "stderr")}
    info.update(case=case.name, run=run.name, exe=str(exe))
    (workdir / ".runinfo.json").write_text(json.dumps(info, indent=1))
    res.check()
    check_complete(case.program, run.stdin, workdir)
    return res


def check_complete(program, stdin, workdir):
    """Raise ProgramFailed unless the program produced complete results.

    Many legacy error paths print a message and STOP with exit code 0 (e.g.
    'TOO LOW ENERGY FOR AVAILABLE PHASE SHIFTS'). refcalc: fd.out must
    contain every energy of the EI,EF,DE range on line 2 of the FIN.
    rfactor: ROUT must contain the best-fit summary line.
    """
    workdir = Path(workdir)
    if program == "refcalc":
        ei, ef, de = (float(x) for x in Path(stdin).read_text().splitlines()[1].split()[:3])
        n_expected = int((ef - ei) / de + 1e-4) + 1
        fd = workdir / "fd.out"
        found = []
        if fd.is_file():
            found = re.findall(r"(?m)^\s*(\d+\.\d\d) \d\.\d{4}\s", fd.read_text(errors="replace"))
        if len(found) != n_expected:
            raise ProgramFailed(f"incomplete output: fd.out has {len(found)} energies, "
                                f"expected {n_expected} (EI={ei}, EF={ef}, DE={de}); "
                                f"stdout tail:\n{(workdir / '.stdout').read_text(errors='replace')[-1000:] if (workdir / '.stdout').is_file() else ''}")
    elif program == "rfactor":
        rout = workdir / "ROUT"
        if not (rout.is_file() and "AVERAGE R-FACTOR =" in rout.read_text(errors="replace")):
            raise ProgramFailed("incomplete output: ROUT has no 'AVERAGE R-FACTOR =' summary")
