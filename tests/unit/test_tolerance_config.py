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


# Entries of tolerances.toml[numeric] that may use skip_lines or
# ignore_line_patterns. Skipping lines hides differences, so an entry is
# added here only together with a decisions entry that carries the evidence
# that the skipped lines hold no physics (the same rule as for loosening a
# tolerance), and the lines must still contain no numbers, see
# test_skipped_lines_contain_no_numbers. Empty: the generators write the
# fixed title 'tl-test <case>' into all inputs, so outputs contain no
# directory names or timestamps and no line has to be skipped.
SKIP_ALLOWLIST = frozenset()

_SKIP_KEYS = ("skip_lines", "ignore_line_patterns")


def _unlisted_skips(numeric_cfg, allowlist):
    """Return the entries that skip lines without being on the allowlist."""
    return sorted(k for k, v in numeric_cfg.items()
                  if any(v.get(key) for key in _SKIP_KEYS) and k not in allowlist)


def test_lines_skipped_only_for_allowlisted_entries():
    """Fail if tolerances.toml skips lines in an entry not in SKIP_ALLOWLIST."""
    assert _unlisted_skips(load_tolerances()["numeric"], SKIP_ALLOWLIST) == []


def test_skip_allowlist_detects_unlisted_entries():
    """A skip on an unlisted entry is reported; listing the entry allows it."""
    cfg = {"a.out": {"skip_lines": [1]}, "b.out": {"ignore_line_patterns": ["x"]},
           "c.out": {"rtol": 1e-6}, "d.out": {"skip_lines": []}}
    assert _unlisted_skips(cfg, frozenset()) == ["a.out", "b.out"]
    assert _unlisted_skips(cfg, frozenset({"a.out"})) == ["b.out"]


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
