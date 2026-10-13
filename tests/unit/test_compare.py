import pytest
from pathlib import Path
from harness.compare import Tolerance, compare_texts, compare_files

TOL = Tolerance(rtol=1e-5, atol_frac=1e-7)
EXACT = Tolerance(rtol=0.0, atol_frac=0.0)  # format cases: values unchanged


@pytest.mark.parametrize("actual,expected,tol,ok", [
    pytest.param("1.000001E+00", "1.000000E+00", TOL, True,
                 id="small_float_deviation_within_tolerance"),
    pytest.param("3", "4", TOL, False, id="int_mismatch_fails"),
    pytest.param("AVERAGE", "AV.-INT", TOL, False, id="text_mismatch_fails"),
    pytest.param("", "1.0", TOL, False, id="empty_actual_fails"),
    pytest.param(" -.1234  -.5E+00\n", "-0.1234 -0.5E+00\n", EXACT, True,
                 id="dropped_leading_zero_of_negative_number_is_same_format"),
    # Ew.d writes 0.1234E-99 but 0.1234-100: the same edit descriptor
    pytest.param(" 0.1234-100\n", " 0.1234E-99\n", Tolerance(0.0, 0.0, atol=1e-99), True,
                 id="three_digit_exponent_without_e_is_same_format_as_ew_d"),
    pytest.param(" 1.0\n\n 2.0\n", " 1.0\n 2.0\n", EXACT, False,
                 id="blank_lines_count_for_line_structure"),
    pytest.param(" -.1234E-01\n", "-0.1234E-01\n", EXACT, True,
                 id="dropped_leading_zero_of_e_mantissa_is_same_format"),
    # ES writes zero as 0.0000E+00, like E: ambiguous, counts as E
    pytest.param(" 0.0000E+00 -0.0000E+00\n", " 0.0000E+00  -.0000E+00\n", EXACT, True,
                 id="zero_value_is_plain_e_format"),
])
def test_compare_texts_ok(actual, expected, tol, ok):
    assert compare_texts(actual, expected, tol).ok is ok


def test_identical_passes():
    r = compare_texts("1.0E+00 2 x", "1.0E+00 2 x", TOL)
    assert r.ok and r.n_compared == 3

def test_large_float_deviation_fails():
    r = compare_texts("1.001E+00", "1.000E+00", TOL)
    assert not r.ok and r.first_failures

def test_atol_scales_with_largest_expected_value():
    # 1e-3 absolute deviation on a tiny value is fine if the file contains 1e5
    # scale = 1e5: atol_frac 1e-9 -> atol 1e-4 < 1e-3 (fail); atol_frac 1e-7 -> atol 1e-2 (pass)
    assert not compare_texts("0.1E+06 0.1E-02", "0.1E+06 0.0E+00", Tolerance(0.0, 1e-9)).ok
    assert compare_texts("0.1E+06 0.1E-02", "0.1E+06 0.0E+00", Tolerance(0.0, 1e-7)).ok

def test_fewer_tokens_fails():
    r = compare_texts("1.0 2.0", "1.0 2.0 3.0", TOL)
    assert not r.ok and "token count" in r.first_failures[0]

def test_missing_file_fails(tmp_path: Path):
    exp = tmp_path / "e"; exp.write_text("1.0")
    r = compare_files(tmp_path / "missing", exp, TOL)
    assert not r.ok and "missing" in r.first_failures[0]

def test_ignored_title_line():
    assert compare_texts("dirA 1\n 1.0", "dirB 2\n 1.0", TOL, ignore_line_patterns=[r"^dir"]).ok

def test_overflow_in_actual_fails():
    r = compare_texts("0.1E+01********", "0.1E+01 0.2E+01", TOL)
    assert not r.ok and "overflow" in r.first_failures[0]

def test_absolute_atol_one_unit_in_last_digit():
    # ROUT prints R factors with 4 decimals: 1 unit passes, 2 units fail,
    # independent of the largest value in the file
    tol = Tolerance(rtol=0.0, atol_frac=0.0, atol=1.5e-4)
    assert compare_texts("479.50 0.2055", "479.50 0.2054", tol).ok
    assert not compare_texts("479.50 0.2056", "479.50 0.2054", tol).ok
    assert not compare_texts("479.51 0.2054", "479.50 0.2054", tol).ok

def test_ulps_tolerance_is_per_value():
    tol = Tolerance(rtol=0.0, atol_frac=0.0, ulps=2)
    # energy 210 does not loosen the check of a small intensity
    assert compare_texts(" 210.00 0.16633E-02", " 210.00 0.16631E-02", tol).ok
    assert not compare_texts(" 210.00 0.16634E-02", " 210.00 0.16631E-02", tol).ok


def test_line_scale_uses_largest_value_of_the_line():
    tol = Tolerance(rtol=0.0, atol_frac=1e-3, scale="line")
    exp = "100.00 1.00\n1.00 1.00\n"
    assert compare_texts("100.00 1.05\n1.00 1.00\n", exp, tol).ok        # 0.05 <= 1e-3*100
    assert not compare_texts("100.00 1.00\n1.00 1.05\n", exp, tol).ok    # 0.05 > 1e-3*1


def test_configured_fd_out_tolerance_detects_one_percent_intensity_change():
    # a per-file scale (max = energy) would pass +19 %
    from harness.cases import CASES_DIR, tolerance_for
    path = CASES_DIR / "refcalc/ag100/expected/E210_L09/fd.out"
    text = path.read_text()
    assert "0.67269E-03" in text
    tol, ignore, skip = tolerance_for("fd.out")
    assert compare_texts(text, text, tol, ignore, skip).ok
    changed = text.replace("0.67269E-03", "0.67942E-03")   # +1 %
    assert not compare_texts(changed, text, tol, ignore, skip).ok


def test_file_cap_limits_the_scaled_part():
    # the line rule alone (1.0 * line max) would allow 1.0; the cap 1e-4 * 100 allows 0.01
    tol = Tolerance(rtol=0.0, atol_frac=1.0, scale="line", file_cap=1e-4)
    exp = "100.000\n  1.000\n"
    assert compare_texts("100.000\n  1.005\n", exp, tol).ok
    assert not compare_texts("100.000\n  1.500\n", exp, tol).ok


def test_configured_amp_out_tolerance_scales_with_the_line_amplitudes():
    # the file max of amp.out is the energy, so under a file-scale rule a 1 %
    # change of a 2e-3 amplitude passes (2.3e-7 * 210 = 4.8e-5 > 2e-5)
    from harness.cases import tolerance_for
    tol, ignore, skip = tolerance_for("amp.out")
    exp = " 210.00 0.0001  0.1000000000E-01  0.2000000000E-02\n"
    file_scale = Tolerance(rtol=1e-6, atol_frac=2.3e-7, ulps=2)
    changed = exp.replace("0.2000000000E-02", "0.2020000000E-02")
    assert compare_texts(changed, exp, file_scale).ok
    assert not compare_texts(changed, exp, tol, ignore, skip).ok
    noise = exp.replace("0.1000000000E-01", "0.1000001000E-01")  # 1e-6 of the line max
    assert compare_texts(noise, exp, tol, ignore, skip).ok


def test_low_precision_f_columns_are_exact():
    # ulps alone would give fd.out energies 18 ulps (0.18 eV)
    from harness.cases import tolerance_for
    tol, ignore, skip = tolerance_for("fd.out")
    exp = "    2   1.00000   0.00000  1\n 138.00 0.0001   0.15603E-02\n"
    assert compare_texts(exp, exp, tol, ignore, skip).ok
    for changed in (exp.replace("138.00", "138.01"), exp.replace("0.0001", "0.0002"),
                    exp.replace("1.00000", "1.00001")):
        assert not compare_texts(changed, exp, tol, ignore, skip).ok
    # the E-format intensity keeps its print-level tolerance
    assert compare_texts(exp.replace("0.15603E-02", "0.15604E-02"), exp, tol, ignore, skip).ok


def test_column_file_energies_are_exact():
    from harness.cases import tolerance_for
    tol, ignore, skip = tolerance_for("theo.column")
    exp = "   48.00  145.55   63.00  999.99\n"
    assert compare_texts("   48.00  145.56   63.00  999.98\n", exp, tol, ignore, skip).ok
    assert not compare_texts("   48.01  145.55   63.00  999.99\n", exp, tol, ignore, skip).ok
    assert not compare_texts("   48.00  145.55   63.01  999.99\n", exp, tol, ignore, skip).ok


def test_exact_columns_by_position():
    tol = Tolerance(rtol=0.1, atol_frac=0.0, exact_columns=(2,))
    assert compare_texts("1.05 2.00 3.05", "1.00 2.00 3.00", tol).ok
    assert not compare_texts("1.00 2.05 3.00", "1.00 2.00 3.00", tol).ok


@pytest.mark.parametrize("actual,expected,what", [
    ("  1.0  2.0\n  3.0\n", "  1.0  2.0  3.0\n", "line count"),            # line split in two
    ("1.0\n2.0 3.0\n", "1.0 2.0\n3.0\n", "token count"),                   # value moved to next line
    (" 0.1234E+01\n", "     1.2340\n", "format"),                          # F -> E
    ("  1.23400\n", "   1.2340\n", "format"),                              # one decimal more
    ("  1.2340\n", "   1.2340\n", "format"),                               # field end column shifted
    ("0.1234E+001\n", " 0.1234E+01\n", "format"),                          # exponent digits 2 -> 3
    ("    12\n", "     12\n", "format"),                                   # integer field end column
    (" 1.2340D-01\n", " 1.2340E-01\n", "format"),                          # exponent letter E -> D
    ("0.1234E+01 0.5000E+00", "  1.234\n  0.50000", "line count"),         # same tokens, other lines
    (" -9.7988896070E-02\n", " -0.9798889607E-01\n", "format"),            # E -> ES
    ("1.23450E-01\n", "0.12345E+00\n", "format"),                          # E -> ES
], ids=["split", "moved", "F_to_E", "decimals", "end_column", "exp_digits", "int_end",
        "exp_letter", "same_tokens_other_lines", "E_to_ES", "E_to_ES_short"])
def test_format_and_line_structure_changes_fail(actual, expected, what):
    r = compare_texts(actual, expected, EXACT)
    assert not r.ok
    assert what in r.first_failures[0] and "line" in r.first_failures[0]


def test_line_count_failure_names_first_differing_line():
    r = compare_texts("x\n  1.0  2.0\n  3.0\n", "x\n  1.0  2.0  3.0\n", EXACT)
    msg = r.first_failures[0]
    assert msg.startswith("line 2:") and "line count" in msg


def test_line_count_failure_names_first_extra_line():
    r = compare_texts("x\n  1.0\n  2.0\n", "x\n  1.0\n", EXACT)
    msg = r.first_failures[0]
    assert msg.startswith("line 3:") and "line count" in msg


def test_format_failure_names_line_and_difference():
    r = compare_texts("x\n 0.1234E+01\n", "x\n     1.2340\n", EXACT)
    msg = r.first_failures[0]
    assert msg.startswith("line 2:") and "E" in msg and "F" in msg
