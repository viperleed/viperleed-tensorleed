"""Shared iteration over cases and runs (harness.cases.iter_runs) and the
per-file insensitivity count shared by the tolerance generators."""
from types import SimpleNamespace

from generate import check_noise_variants
from harness import cases


def _case(name, *run_names):
    return SimpleNamespace(name=name, runs=[SimpleNamespace(name=r) for r in run_names])


def _patch_cases(mocker, by_kind):
    return mocker.patch.object(cases, "load_cases", side_effect=lambda kind: by_kind[kind])


def test_iter_runs_yields_every_run_in_kind_order(mocker):
    """Kinds in the order given, cases and runs in their loaded order."""
    a, b = _case("a", "a1", "a2"), _case("b", "b1")
    _patch_cases(mocker, {"refcalc": [a], "rfactor": [b]})
    assert [(c.name, r.name) for c, r in cases.iter_runs()] == [
        ("a", "a1"), ("a", "a2"), ("b", "b1")]
    assert [(c.name, r.name) for c, r in cases.iter_runs(("rfactor", "refcalc"))] == [
        ("b", "b1"), ("a", "a1"), ("a", "a2")]


def test_iter_runs_restricts_to_kinds_and_names(mocker):
    """A kind not asked for is not loaded; only= filters by case name."""
    a, b = _case("a", "a1"), _case("b", "b1")
    load = _patch_cases(mocker, {"refcalc": [a, b], "rfactor": []})
    assert [c.name for c, _ in cases.iter_runs(("refcalc",), only=["b"])] == ["b"]
    load.assert_called_once_with("refcalc")
    assert len(list(cases.iter_runs(("refcalc",), only=[]))) == 2


def test_expected_insensitive_counts_uses_configured_ignore_and_skip(tmp_path, mocker):
    """The file is tokenized with the ignore/skip settings of its name."""
    exp = tmp_path / "fd.out"
    exp.write_text("header 9.9\n1.0 2.0\n")
    mocker.patch.object(check_noise_variants, "tolerance_for",
                        return_value=(None, ["header"], []))
    tokenize = mocker.spy(check_noise_variants, "tokenize")
    count = mocker.patch.object(check_noise_variants, "insensitive_counts",
                                return_value=(2, 1))
    tol = object()
    assert check_noise_variants.expected_insensitive_counts(exp, tol) == (2, 1)
    tokenize.assert_called_once_with("header 9.9\n1.0 2.0\n", ["header"], [])
    assert count.call_args.args == (tokenize.spy_return, tol)
