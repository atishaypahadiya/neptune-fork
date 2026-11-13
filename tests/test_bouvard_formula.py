"""Bouvard 1821, the printed formula for Uranus: the two entries and unit sanity checks.
No place is computed from it and nothing is fitted."""

import csv
import math
from pathlib import Path

import pytest

from fork.data.centesimal import centesimal_sec_to_arcsec, grades_to_deg

D = Path(__file__).resolve().parents[1] / "data" / "historical" / "adams_challis" / "bouvard_1821"


def read(name: str) -> list[dict[str, str]]:
    with (D / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


T = read("bouvard_1821_uranus_formula.csv")


def val(group: str, label: str) -> str:
    (row,) = [r for r in T if r["group"] == group and r["label"] == label]
    return row["value_as_printed"]


def test_two_entries_agree_in_every_cell():
    p1 = read("passes/bouvard_1821_formula_pass1.csv")
    p2 = read("passes/bouvard_1821_formula_pass2.csv")
    assert len(p1) == len(p2) == len(T) == 104
    assert (
        [r["value"] for r in p1] == [r["value"] for r in p2] == [r["value_as_printed"] for r in T]
    )


def test_saturn_mass_is_kept_as_printed():
    assert val("masses", "m′, denominator") == "5512"
    (row,) = [r for r in T if r["n"] == "002"]
    assert "misprint candidate, open" in row["note"]


def test_precession_is_centesimal():
    p = float(val("V″, first line", "coefficient of t after φ″, seconds"))
    assert centesimal_sec_to_arcsec(p) == pytest.approx(50.1, abs=0.05)  # arcsec a year


def test_mean_motion_is_centesimal_degrees_a_year():
    n = float(val("mean longitudes", "ε″ + n″t, coefficient of t"))
    assert 360.0 / grades_to_deg(n) == pytest.approx(84.0, abs=0.1)  # Uranus's period, years
    nj = float(val("mean longitudes", "ε + nt, coefficient of t, degrees"))
    assert 360.0 / grades_to_deg(nj) == pytest.approx(11.86, abs=0.01)


def test_equation_of_the_centre_matches_the_eccentricity_loosely():
    # First term 2e − e³/4 radians. The series and the eccentricity of p. ii differ by
    # about 0.16 %; why is not established (docs/notes/P1B_bouvard_1821.md).
    e = float(val("eccentricities", "e″"))
    c1 = float(val("V″ group I", "line 1: constant in the bracket"))
    rad = math.radians(centesimal_sec_to_arcsec(c1) / 3600)
    assert rad == pytest.approx(2 * e - e**3 / 4, rel=5e-3)


def test_perihelion_and_node_are_in_range():
    for lab in ("ϖ″, constant", "θ″, constant"):
        g = float(
            val("perihelia and nodes (printed lines only, not the handwriting at the right)", lab)
        )
        assert 0 <= grades_to_deg(g) < 360
    peri = "perihelia and nodes (printed lines only, not the handwriting at the right)"
    assert grades_to_deg(float(val(peri, "ϖ″, constant"))) == pytest.approx(167.5, abs=0.1)
