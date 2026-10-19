"""Unit tests for the shared output-comparison helpers of harness.cases."""
import dataclasses

import pytest

from harness import cases
from harness.compare import Tolerance

TOL = Tolerance(rtol=1e-5, atol_frac=0.0)


@pytest.fixture(autouse=True)
def fixed_tolerance(mocker):
    """Use the same tolerance for every file name."""
    mocker.patch.object(cases, "tolerance_for", return_value=(TOL, [], []))


@pytest.fixture
def dirs(tmp_path):
    """(actual dir, expected dir), both created."""
    act, exp = tmp_path / "act", tmp_path / "exp"
    act.mkdir()
    exp.mkdir()
    return act, exp


@pytest.mark.parametrize("patterns,expected", [
    (["T_*", "*"], "T_*"),
    (["*", "T_*"], "*"),
    (["fd.out"], None),
])
def test_pattern_for_returns_first_match(patterns, expected):
    """The first matching pattern is returned, None if there is none."""
    assert cases.pattern_for("T_1", patterns) == expected


def test_compare_outputs_accepts_matching_files(dirs):
    """Files within tolerance give no failures."""
    act, exp = dirs
    (exp / "a").write_text("1.0000\n")
    (act / "a").write_text("1.0000\n")
    assert cases.compare_outputs(act, [exp / "a"]) == {}


def test_compare_outputs_reports_failing_comparison(dirs):
    """Each file deviating beyond tolerance is reported under its name.

    Failing files come first and last with a matching file between, so a
    failure is not lost when a later file matches or fails.
    """
    act, exp = dirs
    for name, actual in (("bad1", "2.0000"), ("good", "1.0000"), ("bad2", "3.0000")):
        (exp / name).write_text("1.0000\n")
        (act / name).write_text(actual + "\n")
    failures = cases.compare_outputs(act, [exp / n for n in ("bad1", "good", "bad2")])
    assert set(failures) == {"bad1", "bad2"}
    assert "2.0" in failures["bad1"][0]
    assert "3.0" in failures["bad2"][0]


def test_compare_outputs_reports_missing_actual_file(dirs):
    """An expected file that was not produced is a failure."""
    act, exp = dirs
    (exp / "a").write_text("1.0\n")
    assert "a" in cases.compare_outputs(act, [exp / "a"])


def test_compare_outputs_applies_tolerance_adjustment(dirs):
    """The adjust callable maps the configured tolerance to the one used."""
    act, exp = dirs
    (exp / "a").write_text("1.0000\n")
    (act / "a").write_text("1.0005\n")
    assert cases.compare_outputs(act, [exp / "a"])
    relaxed = cases.compare_outputs(
        act, [exp / "a"], adjust=lambda t: dataclasses.replace(t, rtol=1e-3))
    assert relaxed == {}


def test_format_failures_joins_messages_per_file():
    """Each failing file gives one line with its messages joined."""
    assert cases.format_failures({"a": ["x", "y"], "b": ["z"]}) == ["a: x; y", "b: z"]


def test_file_set_failure_empty_when_sets_agree(tmp_path):
    """Produced files matching the globs equal to the expected names are fine."""
    (tmp_path / "fd.out").touch()
    (tmp_path / "T_1").touch()
    assert cases.file_set_failure(tmp_path, ["fd.out", "T_*"], ["T_1", "fd.out"]) == ""


def test_file_set_failure_names_missing_file(tmp_path):
    """An expected file that is not produced is named as missing."""
    (tmp_path / "fd.out").touch()
    msg = cases.file_set_failure(tmp_path, ["fd.out", "T_*"], ["fd.out", "T_1"])
    assert "missing: ['T_1']" in msg
    assert "extra: []" in msg


def test_file_set_failure_names_extra_file(tmp_path):
    """A produced file matching a glob but not expected is named as extra."""
    (tmp_path / "T_1").touch()
    (tmp_path / "T_2").touch()
    (tmp_path / "unrelated").touch()
    msg = cases.file_set_failure(tmp_path, ["T_*"], ["T_1"])
    assert "missing: []" in msg
    assert "extra: ['T_2']" in msg
