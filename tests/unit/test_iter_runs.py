"""Tests of the shared iteration over cases and runs (harness.cases.iter_runs)."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-16'
__license__ = 'GPLv3+'

from types import SimpleNamespace

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
