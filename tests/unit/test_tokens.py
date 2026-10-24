"""Tests of the Fortran output tokenizer of harness.tokens."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2026 ViPErLEED developers'
__created__ = '2026-10-13'
__license__ = 'GPLv3+'

import pytest

from harness.tokens import FortranOverflowError, tokenize


def _vals(text, **kw):
    """Return the (kind, value) pairs of the tokens of text."""
    return [(t.kind, t.value) for t in tokenize(text, **kw)]


@pytest.mark.parametrize("text,expected", [
    pytest.param("0.2069128752E+01-.1653718948E+00",
                 [("float", 2.069128752), ("float", -0.1653718948)],
                 id="concatenated_floats_without_leading_zero"),
    pytest.param("0.123456-100 0.5E+00", [("float", 0.123456e-100), ("float", 0.5)],
                 id="three_digit_exponent_without_E"),
    pytest.param("   45.00 0.0001   7  -3",
                 [("float", 45.0), ("float", 0.0001), ("int", 7), ("int", -3)],
                 id="plain_decimals_and_ints"),
    pytest.param("( 1| 0) BEST GEOMETRY",
                 [("text", "("), ("int", 1), ("text", "|"), ("int", 0), ("text", ")"),
                  ("text", "BEST"), ("text", "GEOMETRY")],
                 id="text_tokens_kept"),
    # rfactor ROUT / WEXPEL beam labels: *( 1  0)*
    pytest.param("*( 1  0)*  1  0.0001",
                 [("text", "*"), ("text", "("), ("int", 1), ("int", 0), ("text", ")"),
                  ("text", "*"), ("int", 1), ("float", 0.0001)],
                 id="single_star_in_beam_label_is_text"),
    # rfactor stdout comment lines: *BEAM NO.  1 ..., * <title>
    pytest.param("*BEAM NO.  1", [("text", "*"), ("text", "BEAM"), ("text", "NO"),
                                  ("text", "."), ("int", 1)],
                 id="star_before_letter_in_column_1_is_text"),
    pytest.param("* tl-test", [("text", "*"), ("text", "tl"), ("text", "-"),
                               ("text", "test")],
                 id="star_before_blank_in_column_1_is_text"),
    # T_* files: 0.142878E+02 followed directly by 0.610153E-06
    pytest.param("0.142878E+020.610153E-06", [("float", 14.2878), ("float", 0.610153e-06)],
                 id="exponent_has_two_digits_when_next_number_starts_with_digit"),
    pytest.param("0.4231456220-1000.4308617413E+00",
                 [("float", 0.4231456220e-100), ("float", 0.4308617413)],
                 id="three_digit_exponent_followed_by_positive_number"),
    pytest.param("-.1234-100-.5678E+00", [("float", -0.1234e-100), ("float", -0.5678)],
                 id="three_digit_exponent_without_leading_zero_concatenated"),
    pytest.param("12.345-123.456", [("float", 12.345), ("float", -123.456)],
                 id="f_fields_are_not_read_as_three_digit_exponents"),
    # not 0.1234E-123 and .45, which would let a change of the integer part pass
    pytest.param("0.1234-123.45", [("float", 0.1234), ("float", -123.45)],
                 id="f_field_followed_by_f_field_with_three_integer_digits"),
    pytest.param("1.5E+5", [("float", 1.5e5)], id="one_digit_e_exponent"),
    pytest.param("1.0E-100 2.0", [("float", 1.0e-100), ("float", 2.0)],
                 id="three_digit_e_exponent_followed"),
    pytest.param("1.0E-100", [("float", 1.0e-100)], id="three_digit_e_exponent"),
])
def test_tokens(text, expected):
    """The tokenizer splits Fortran-formatted text into the listed kinds and values."""
    assert _vals(text) == expected


@pytest.mark.parametrize("text", [
    pytest.param("0.1E+01********0.2E+01", id="overflow_field"),
    pytest.param("  1 ** 2", id="two_star_overflow"),
    # a 1-wide field that overflowed, e.g. I4 or F4.1 written as '   *'
    pytest.param("  12   *  0.5000", id="lone_star_in_numeric_column"),
    pytest.param("0.1E+01*", id="lone_star_next_to_number"),
    pytest.param("*1  0.5000", id="lone_star_in_column_1_before_digit"),
    pytest.param("** BEAM", id="two_stars_in_column_1"),
    pytest.param(" *BEAM", id="star_before_letter_not_in_column_1"),
])
def test_overflow_raises(text):
    """Overflow fields (asterisks) raise FortranOverflowError."""
    with pytest.raises(FortranOverflowError):
        tokenize(text)


def test_overflow_mode_returns_overflow_tokens():
    """With overflow=True, an overflow field is a token with its text and end column."""
    toks = tokenize(" 1 *******  0.5000 *( 1  0)*", overflow=True)
    assert [(t.kind, t.value) for t in toks if t.kind != "text"] == [
        ("int", 1), ("overflow", "*******"), ("float", 0.5), ("int", 1), ("int", 0)]
    assert next(t for t in toks if t.kind == "overflow").end == 10


def test_ignored_lines_are_skipped_and_line_numbers_kept():
    """Ignored lines give no tokens, and later tokens keep their line numbers."""
    toks = tokenize("Ag100_kept 261005\n 1.5\n", ignore_line_patterns=[r"^Ag100"])
    assert [(t.kind, t.value, t.line) for t in toks] == [("float", 1.5, 2)]


def test_skip_lines_by_number():
    """Lines skipped by number give no tokens."""
    toks = tokenize("title 12 3.0\n 1.5\n", skip_lines=[1])
    assert [(t.kind, t.value, t.line) for t in toks] == [("float", 1.5, 2)]


def test_float_tokens_carry_unit_in_last_digit():
    """Each float token knows the value of one unit in its last digit."""
    toks = tokenize("0.16631E-02 -.1653718948E+00 45.00 0.123456-100")
    assert [t.ulp for t in toks] == pytest.approx([1e-7, 1e-10, 1e-2, 1e-106])
