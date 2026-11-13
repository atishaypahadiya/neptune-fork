"""Centesimal to sexagesimal conversion for Bouvard's 1821 tables."""

import pytest
from hypothesis import given
from hypothesis import strategies as st

from fork.data import centesimal as c


def test_the_circle_and_the_right_angle():
    assert c.grades_to_deg(400) == pytest.approx(360)
    assert c.grades_to_deg(100) == pytest.approx(90)
    assert c.centesimal_dms_to_deg(100, 0, 0) == pytest.approx(90)


def test_units_of_each_order():
    assert c.grades_to_deg(1) == pytest.approx(0.9)  # 54 arcminutes
    assert c.centesimal_dms_to_deg(0, 1, 0) * 3600 == pytest.approx(32.4)  # one minute
    assert c.centesimal_dms_to_deg(0, 0, 1) * 3600 == pytest.approx(0.324)  # one second
    assert c.centesimal_sec_to_arcsec(1) == pytest.approx(0.324)
    assert c.centesimal_sec_to_arcsec(1_000_000) == pytest.approx(90 * 3600)


def test_a_second_is_not_an_arcsecond():
    # The trap: reading a centesimal value as sexagesimal overstates it by 1 / 0.324.
    assert c.centesimal_sec_to_arcsec(100.0) == pytest.approx(32.4)
    assert c.centesimal_sec_to_arcsec(100.0) != pytest.approx(100.0, rel=0.5)


def test_printed_forms_agree():
    # 185° 79' 43" centesimal is 185.7943 centesimal degrees.
    assert c.centesimal_dms_to_deg(185, 79, 43) == pytest.approx(c.grades_to_deg(185.7943))
    assert c.centesimal_dms_to_deg(185, 79, 43) == pytest.approx(167.21487)


def test_precession_constant_in_both_divisions():
    # Bouvard's uniform precession, 154".63 centesimal a year, is 50".10 sexagesimal.
    assert c.centesimal_sec_to_arcsec(154.63) == pytest.approx(50.10012)


def test_decimal_time():
    assert c.decimal_time_to_day(5, 0, 0) == pytest.approx(0.5)  # noon
    assert c.decimal_time_to_day(0, 0, 1) * 86400 == pytest.approx(0.864)
    assert c.decimal_time_to_day(9, 99, 99) < 1


@pytest.mark.parametrize("bad", [(1, 100, 0), (1, 0, 100), (-1, 0, 0), (1, -1, 0)])
def test_rejects_what_is_not_centesimal(bad):
    with pytest.raises(ValueError):
        c.centesimal_dms_to_deg(*bad)


def test_rejects_sexagesimal_hours():
    with pytest.raises(ValueError):
        c.decimal_time_to_day(12, 0, 0)


@given(st.integers(0, 399), st.integers(0, 99), st.integers(0, 99))
def test_round_trip_through_arcseconds(d, m, s):
    deg = c.centesimal_dms_to_deg(d, m, s)
    total_centesimal_seconds = d * 10000 + m * 100 + s
    assert deg * 3600 == pytest.approx(
        c.centesimal_sec_to_arcsec(total_centesimal_seconds), abs=1e-6
    )
    assert 0 <= deg < 360
