"""Checks of tests/tolerances.toml against the stored expected outputs."""
import fnmatch

from harness.cases import CASES_DIR, load_tolerances
from harness.tokens import kept_lines, tokenize


def _expected_files():
    return [p for p in CASES_DIR.glob("*/*/expected/*/*") if p.is_file()]


def test_every_expected_output_has_a_tolerance():
    """Every stored expected file name is matched by a tolerance entry."""
    patterns = list(load_tolerances()["numeric"])
    missing = {p.name for p in _expected_files()
               if not any(fnmatch.fnmatchcase(p.name, pat) for pat in patterns)}
    assert not missing


# Entries of tolerances.toml[numeric] that may use skip_lines or
# ignore_line_patterns. Skipping lines hides differences, so an entry is
# added here only with evidence that the skipped lines hold no physics, and
# the PR that adds it states that evidence (the same rule as for loosening a
# tolerance, README 'Rules'). The skipped or ignored lines must still contain
# no numbers, see test_skipped_lines_contain_no_numbers. Empty: the generators
# write the fixed title 'tl-test <case>' into all inputs, so outputs contain
# no directory names or timestamps and no line has to be skipped.
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
    """The only title line in expected outputs is the fixed title of the case."""
    for path in _expected_files():
        case = path.parents[2].name
        for line in path.read_text(errors="replace").splitlines():
            if "tl-test" in line:
                assert line.strip() in (f"tl-test {case}", f"*tl-test {case}"), (path, line)


def _dropped_text(text, ignore_line_patterns, skip_lines):
    """The lines of text that comparison leaves out, joined by newlines."""
    kept = {lineno for lineno, _ in kept_lines(text, ignore_line_patterns, skip_lines)}
    return "\n".join(line for lineno, line in enumerate(text.splitlines(), start=1)
                     if lineno not in kept)


def test_skipped_lines_contain_no_numbers():
    """Lines that tolerances.toml skips or ignores hold no floating-point numbers.

    Covers both skip_lines and ignore_line_patterns, for entries on
    SKIP_ALLOWLIST as for any other.
    """
    for pattern, cfg in load_tolerances()["numeric"].items():
        ignore, skip = cfg.get("ignore_line_patterns", []), cfg.get("skip_lines", [])
        if not (ignore or skip):
            continue
        for path in _expected_files():
            if fnmatch.fnmatchcase(path.name, pattern):
                text = _dropped_text(path.read_text(errors="replace"), ignore, skip)
                assert not [t for t in tokenize(text) if t.kind == "float"], (path, text)


def test_dropped_text_has_skipped_and_ignored_lines():
    """_dropped_text returns the lines dropped by skip_lines and by a pattern."""
    text = "title\n1.5 2.5\n# note 3.5\n4.5\n"
    assert _dropped_text(text, [r"^#"], [1]) == "title\n# note 3.5"


def test_perf_and_debug_keys():
    """The [perf] and [debug] keys exist with plausible types and ranges."""
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
    """The [tier2] section pins a non-empty viperleed version."""
    version = load_tolerances()["tier2"]["viperleed_version"]
    assert isinstance(version, str) and version.strip()
