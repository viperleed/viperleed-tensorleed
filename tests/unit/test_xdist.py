"""pytest -n (pytest-xdist): which -m expressions are refused (tests/conftest.py)."""
import pytest

from conftest import _n_workers, _serial_only_selected, _tier2_selected


@pytest.mark.parametrize("markexpr,refused", [
    ("perf", ["perf"]),
    ("benchmark", ["benchmark"]),
    ("tier1 or perf", ["perf"]),
    ("not tier1", ["perf", "benchmark"]),
    ("", ["perf", "benchmark"]),                      # -m "" selects everything
    ("tier1", []),
    ("debug", []),
    ("tier2", []),
    ("not perf and not debug and not benchmark and not tier2", []),   # addopts
])
def test_serial_only_selected(markexpr, refused):
    assert _serial_only_selected(markexpr) == refused


@pytest.mark.parametrize("markexpr,selected", [
    ("tier2", True),
    ("tier1 or tier2", True),
    ("not tier1", True),
    ("", True),                                       # -m "" selects everything
    ("tier1", False),
    ("debug", False),
    ("not perf and not debug and not benchmark and not tier2", False),   # addopts
])
def test_tier2_selected(markexpr, selected):
    assert _tier2_selected(markexpr) is selected


@pytest.mark.parametrize("tx,n", [
    ([], 0),
    (["popen"] * 3, 3),
    (["2*popen"], 2),
    (["3*popen//python=python3", "popen"], 4),
])
def test_n_workers(tx, n):
    assert _n_workers(tx) == n
