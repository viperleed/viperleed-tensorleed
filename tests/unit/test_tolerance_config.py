"""Checks of tests/tolerances.toml against the stored expected outputs."""
import fnmatch

from harness.cases import CASES_DIR, load_tolerances
from harness.tokens import tokenize


def _expected_files():
    return [p for p in CASES_DIR.glob("*/*/expected/*/*") if p.is_file()]


def test_every_expected_output_has_a_tolerance():
    patterns = list(load_tolerances()["numeric"])
    missing = {p.name for p in _expected_files()
               if not any(fnmatch.fnmatchcase(p.name, pat) for pat in patterns)}
    assert not missing


def test_no_lines_skipped_because_titles_are_fixed():
    # The generators write the fixed title 'tl-test <case>' into all inputs,
    # so outputs contain no directory names or timestamps and no line has
    # to be skipped. Skipping lines would hide differences.
    cfg = load_tolerances()["numeric"]
    assert {k: v.get("skip_lines", []) for k, v in cfg.items() if v.get("skip_lines")} == {}
    assert {k: v.get("ignore_line_patterns", []) for k, v in cfg.items()
            if v.get("ignore_line_patterns")} == {}


def test_titles_in_outputs_are_the_fixed_case_title():
    for path in _expected_files():
        case = path.parents[2].name
        for line in path.read_text(errors="replace").splitlines():
            if "tl-test" in line:
                assert line.strip() in (f"tl-test {case}", f"*tl-test {case}"), (path, line)


def test_skipped_lines_contain_no_numbers():
    # Guard for future configs that do skip lines.
    for pattern, cfg in load_tolerances()["numeric"].items():
        skip = set(cfg.get("skip_lines", []))
        for path in _expected_files():
            if skip and fnmatch.fnmatchcase(path.name, pattern):
                lines = path.read_text(errors="replace").splitlines()
                text = "\n".join(lines[i - 1] for i in sorted(skip) if i <= len(lines))
                assert not [t for t in tokenize(text) if t.kind == "float"], (path, text)


def test_perf_and_debug_keys():
    perf = load_tolerances()["perf"]
    assert type(perf["repetitions"]) is int and perf["repetitions"] >= 3
    assert type(perf["benchmark_repetitions"]) is int and perf["benchmark_repetitions"] >= 3
    assert type(perf["sizing_runs"]) is int and 1 <= perf["sizing_runs"] <= 10
    assert type(perf["max_inner"]) is int and perf["max_inner"] >= 1
    assert isinstance(perf["min_sizing_cpu_s"], float) and 0 < perf["min_sizing_cpu_s"] < 0.1
    for key in ("sample_cpu_s", "case_sample_cpu_s"):
        assert isinstance(perf[key], float) and 0 < perf[key] < 10
    min_rtol = load_tolerances()["debug"]["min_rtol"]
    assert isinstance(min_rtol, float) and 0 < min_rtol <= 1e-2


def test_tier2_viperleed_version_key():
    version = load_tolerances()["tier2"]["viperleed_version"]
    assert isinstance(version, str) and version.strip()
