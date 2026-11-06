"""M6a double-entry tooling: pass diff and per-row arithmetic check (fast)."""

import pytest

from fork.data import transcription as tr


def _row(page, line, day, t, **kw):
    r = dict.fromkeys(tr.COLUMNS, "")
    r.update(page=page, line=line, day=day, mean_time=t, **kw)
    return r


def test_diff_aligns_rows_so_a_skipped_row_does_not_cascade():
    a = [_row("129", "1", "23", "9 41 25", ra_obs="55 49 19.7"),
         _row("129", "2", "2", "9 46 47"),
         _row("129", "3", "4", "12 43 35", lon_c_minus_o="-65.4")]  # fmt: skip
    b = [
        _row("129", "1", "23", "9 41 25", ra_obs="55 49 19.1"),  # one digit differs
        _row("129", "2", "4", "12 43 35", lon_c_minus_o="-65.4"),
    ]  # pass 2 skipped a row
    d = tr.diff_passes(a, b)
    assert [(x.field, x.pass1, x.pass2) for x in d] == [
        ("ra_obs", "55 49 19.7", "55 49 19.1"),
        ("row", "2|9 46 47", ""),
    ]
    assert d[1].line1 == "2" and d[1].line2 is None
    assert tr.diff_passes(a, a) == []


def test_arithmetic_check_on_printed_rows():
    """Rows as printed on p. 129 (1690 Dec 23; 1781 Sep 25) and p. 133 (1811 Feb 17,
    southern Dec): printed C − O matches computed seconds − observed seconds."""
    r1690 = _row("129", "1", "23", "9 41 25", ra_obs="55 49 19.7", ra_calc_sec="12.8",
                 ra_c_minus_o="-66.9", dec_obs="+19 35 14.4", dec_calc_sec="54.1",
                 dec_c_minus_o="-20.3")  # fmt: skip
    c = tr.arithmetic_check(r1690)
    assert c["ra"] == pytest.approx(0.0, abs=1e-6) and c["dec"] == pytest.approx(0.0, abs=1e-6)
    r1781 = _row("129", "15", "25", "18 0 15", ra_obs="6 12 26.63", ra_calc_sec="25.39",
                 ra_c_minus_o="-1.24", dec_obs="» » »", dec_calc_sec="»")  # fmt: skip
    c = tr.arithmetic_check(r1781)
    assert c["ra"] == pytest.approx(0.0, abs=1e-6) and c["dec"] is None
    south = _row("133", "1", "17", "17 16 56", dec_obs="-17 6 13.3", dec_calc_sec="1.2",
                 dec_c_minus_o="+12.1")  # fmt: skip
    assert tr.arithmetic_check(south)["dec"] == pytest.approx(0.0, abs=1e-6)
    bad = {**r1690, "ra_c_minus_o": "-66.4"}  # a misread digit
    assert tr.arithmetic_check(bad)["ra"] == pytest.approx(0.5, abs=1e-6)
