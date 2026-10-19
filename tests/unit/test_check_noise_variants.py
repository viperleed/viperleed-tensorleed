"""Tests of the per-file insensitivity count of generate.check_noise_variants.

The count is shared with calibrate_tolerances.
"""
from generate import check_noise_variants


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
