"""Step 3d, validation: the transcriptions for V1 and V2, and V1's outcome
(docs/notes/P1B_bouvard_attribution.md, § 5 and the results below its lines)."""

import csv
from pathlib import Path

import pytest

from fork.export import p1b_bouvard as P
from fork.historical import bouvard_1821 as B

D = Path(__file__).resolve().parents[1] / "data" / "historical" / "adams_challis"
PASSES = D / "bouvard_1821" / "passes"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


@pytest.mark.parametrize(
    ("stem", "n"),
    [
        ("bouvard_1821_table_rows", 764),
        ("bouvard_1821_great_inequalities", 24),
        ("adams_memoir_s7_terms", 15),
        ("bouvard_1821_table_rows_v1p", 609),
    ],
)
def test_two_entries_agree_in_every_cell(stem, n):
    a, b = read(PASSES / f"{stem}_pass1.csv"), read(PASSES / f"{stem}_pass2.csv")
    key = lambda r: (r["table"], r["row"], r["column"])  # noqa: E731
    va, vb = {key(r): r["value"] for r in a}, {key(r): r["value"] for r in b}
    assert len(va) == len(vb) == n
    assert va == vb


def test_reconciled_rows_differ_from_the_reading_only_where_logged():
    rows = read(D / "bouvard_1821" / "bouvard_1821_uranus_table_rows.csv")
    changed = {
        (r["table"], r["row"], r["column"]): (r["value_as_printed"], r["value_used"])
        for r in rows
        if r["value_as_printed"] != r["value_used"]
    }
    assert changed == {
        ("X", "324", "eq"): ("4.29.83,8", "4.28.83,8"),
        ("X", "327", "eq"): ("4.39.52,6", "4.38.52,6"),
        ("X", "25", "eq"): ("2.35.14,3", "2.33.14,3"),
        ("XII", "2600", "eq"): ("3.41,6", "3.31,6"),
    }
    assert all(r["note"] for r in rows if r["value_as_printed"] != r["value_used"])


def test_table_arguments_in_parts_of_the_circle():
    assert B.centesimal_seconds("111.58.62,6") == pytest.approx(1_115_862.6)
    assert B.centesimal_seconds("81.02.48") == pytest.approx(810_248.0)
    assert B.centesimal_seconds("1.21,3") == pytest.approx(121.3)
    assert B.centesimal_seconds("26,3") == pytest.approx(26.3)
    assert B.ARG_PARTS["III"] == 10_000 and B.ARG_PARTS["II"] == 1_000


def test_table_ix_is_the_great_inequality_of_p_iii():
    for inst in B.INSTANTS:
        diff = B.formula(inst.t)["great_inequality"] - B.table_longitude(inst)["great_inequality"]
        assert abs(diff) < 0.2  # centesimal seconds; the table is printed to 0.1


def test_s7_terms():
    terms = B.s7_terms()
    assert [t[0] for t in terms] == pytest.approx([1.918, 1.085, 32.00, -8.35, -1.49])
    assert [t[1] for t in terms] == [(0, 1, -2), (1, 0, -1), (2, -6, 3), (2, -6, 2), (2, -6, 4)]
    assert terms[0][2] == pytest.approx(-(13 + 1.5 / 60))


def test_rounding_budgets_are_those_of_the_plan():
    b = {i.name: B.table_longitude(i) for i in B.INSTANTS}
    assert [b[n]["entries"] for n in ("1783", "1798", "1813", "1828", "1846", "1715")] == [
        17, 19, 19, 19, 19, 20,
    ]  # fmt: skip
    assert b["1783"]["budget_arcsec"] == pytest.approx(0.2754)
    assert b["1846"]["budget_arcsec"] == pytest.approx(0.3078)


def test_v1_fails_as_reported():
    r = P.v1()
    got = [r["instants"][n]["formula_minus_tables_arcsec"] for n in P.V1_NAMES]
    assert got == pytest.approx([-2.93, 0.26, 3.40, 2.88, -0.38], abs=0.01)
    assert r["outcome"] == "fail"
    assert r["over_budget"] == ["1783", "1813", "1828", "1846"]


def test_the_difference_is_the_centre_coefficient_and_the_mean_motion():
    r = P.v1()
    rows = [r["instants"][n] for n in P.V1_NAMES]
    # Table X answers to 59417″.54, ten seconds under the printed 59427″.54
    assert [x["with_centre_from_tables_cs"] for x in rows] == pytest.approx(
        [0.15, 1.12, 1.55, 1.04, 3.27], abs=0.01
    )
    # Table I's mean longitude falls behind p. ii by a tenth of a second a year from 1780
    for x in rows[1:]:
        assert x["mean_longitude_drift_cs_per_year_since_1780"] == pytest.approx(0.105, abs=0.002)
    left = [abs(x["with_centre_from_tables_less_mean_and_precession_cs"]) for x in rows]
    assert max(left) == pytest.approx(1.97, abs=0.01)
    assert all(
        v <= x["budget_cs"] + x["argument_allowance_cs"] for v, x in zip(left, rows, strict=True)
    )
    assert r["eccentricity"]["of_table_x_59417.54"] == pytest.approx(0.046679, abs=2e-6)


def test_the_sign_before_group_iii_is_minus():
    assert P.v1()["sign_before_group_iii"]["largest_difference_arcsec_with_plus"] > 80.0


def test_corrected_formula_at_the_v1_instants_is_the_post_hoc_removal():
    r = P.v1p()["not_part_of_the_rule"]
    got = [r[n]["corrected_minus_tables_arcsec"] for n in (*P.V1_NAMES, "1715")]
    assert got == pytest.approx([0.05, -0.06, -0.23, -0.57, -0.04, 0.73], abs=0.01)


def test_v1_prime_as_reported():
    r = P.v1p()
    rows = r["instants"]
    names = (*P.MODERN_V1P, *P.ANCIENT_V1P)
    assert [rows[n]["corrected_minus_tables_arcsec"] for n in names] == pytest.approx(
        [-0.50, 0.00, 0.02, -0.12, -1.05, 0.65, -0.48], abs=0.01
    )
    assert [rows[n]["entries"] for n in names] == [18, 19, 19, 19, 20, 20, 19]
    assert all(
        abs(rows[n]["corrected_minus_tables_arcsec"]) <= rows[n]["limit_arcsec"] for n in names
    )
    # one instant passes only with a reading replaced by what the printed differences give
    assert rows["1789"]["corrected_minus_tables_as_read_arcsec"] == pytest.approx(-2.40, abs=0.01)
    assert r["outcomes"] == {n: "stop and ask" if n == "1789" else "pass" for n in names}
    assert r["consequence"] == "stop and ask"
    # the mean motion and the equation of the centre agree at every fresh instant
    for n in names:
        assert abs(rows[n]["parts_cs"]["mean_longitude"]) < 0.1
        assert abs(rows[n]["parts_cs"]["equation_of_centre"]) < 0.4


def test_v1_prime_ruling_all_pass_with_the_adopted_readings():
    r = P.v1p()["ruled_2026_11_13"]
    assert set(r["outcomes"].values()) == {"pass"}


def test_misprint_log_matches_the_rows_used():
    log = read(D / "bouvard_1821" / "bouvard_1821_misprints.csv")
    rows = {
        (r["table"], r["row"]): r
        for r in read(D / "bouvard_1821" / "bouvard_1821_uranus_table_rows.csv")
        if r["column"] in ("eq", "red")
    }
    for entry in log:
        if not entry["place"].startswith("Table"):
            continue
        table, row = (
            entry["place"].split(",")[0].split()[1],
            entry["place"].split("row ")[1].split(",")[0],
        )
        assert rows[table, row]["value_as_printed"] == entry["reads"]
        assert rows[table, row]["value_used"] == entry["used"]


def test_v2_s8_times_the_factor_is_s9():
    r = P.v2()
    assert r["s8_to_s9_worst_arcsec"]["modern"] <= P.TOL_S9_MODERN
    assert r["s8_to_s9_worst_arcsec"]["ancient"] <= P.TOL_S9_ANCIENT
    assert r["s8_to_s9_pass"]


def test_v2_s7_terms_at_his_epochs():
    r = P.v2()
    lo, hi = r["s7_range_arcsec"]
    assert lo == pytest.approx(-29.91, abs=0.02) and hi == pytest.approx(-9.89, abs=0.02)
    # the 32″ term is nearly a constant and a trend over his 150 years
    assert r["each_term_after_constant_and_trend_max_abs_arcsec"][2] < 1.0
    assert sorted(r["s5_shift_arcsec"]) == [1690, 1789]
    assert abs(r["s5_shift_arcsec"][1690]) < 0.03
    assert r["s5_shift_arcsec"][1789] == pytest.approx(1.91, abs=0.01)


def test_c_bouvard_is_the_corrected_formula_plus_s7():
    inst = B.INSTANTS_V1P[2]
    want = B.corrected_formula(inst)["longitude"] * 0.9e-4 + B.s7_correction(inst.t) / 3600.0
    assert B.c_bouvard_deg(inst.t) == pytest.approx(want % 360.0, abs=1e-6)


def test_v3_gates_and_frame_tests_as_recorded():
    """V3 needs the ephemeris at his 30 epochs, so the recorded run is checked, not rerun."""
    import json

    v3 = json.loads((P.OUT_JSON).read_text(encoding="utf-8"))["V3"]
    for hyp, gates in P.GATES.items():
        for name, must in gates.items():
            assert v3["gates"][hyp][name]["gives"] == pytest.approx(must, abs=P.TOL_GATE_DEG)
    for method in ("d1c", "chain"):
        assert abs(v3["frame_test_offset"][method]["change_deg"]) <= P.TOL_FRAME_OFFSET_DEG
    assert v3["frame_test_curve"]["max_abs_difference_arcsec"] <= P.TOL_FRAME_CURVE_ARCSEC
    assert v3["propagation_check"]["d_arcsec"] == pytest.approx(1.05)


def test_attribution_algebra_and_classes():
    a = P.attribute(1.6, -3.2, -3.8, -9.5)
    assert a["R"] + a["Ob"] == pytest.approx(a["A"] - a["D"])
    assert a["R"] - a["Ob"] == pytest.approx(a["C"] - a["B"])
    assert a["I"] == pytest.approx(0.5 * ((1.6 + 3.2) - (-3.8 + 9.5)))
    assert [P.effect_class(v) for v in (6.0, 5.99, 3.01, 3.0, -3.0, -3.01)] == [
        "carries", "partial", "partial", "does not carry", "does not carry", "works against",
    ]  # fmt: skip


@pytest.mark.parametrize(
    ("cells", "row", "override"),
    [
        ((2.0, -8.0, 1.0, -9.0), "The reference carries the shift.", False),
        ((2.0, 1.0, -8.0, -9.0), "The observations carry the shift.", False),
        ((2.0, -3.5, -3.5, -9.0), "Both, in the shares found.", False),
        ((2.0, -12.0, 6.0, -9.0), "The reference carries more than the whole shift", False),
        ((2.0, -9.0, -9.0, -9.0), "Both, in the shares found.", True),
    ],
)
def test_reading_rows_on_constructed_effects(cells, row, override):
    a = P.attribute(*cells)
    assert a["row"].startswith(row)
    assert a["interaction_override"] is override


def test_cells_as_recorded():
    """The cells take minutes and need the ephemeris cache: the recorded run is checked."""
    import json

    c = json.loads(P.OUT_JSON.read_text(encoding="utf-8"))["cells"]
    h1 = c["cells"]["H1"]
    for cell, key in (("A", "A"), ("D", "D")):
        assert h1[cell]["d1c"]["offset"] == pytest.approx(P.GATES["H1"][f"{key}_d1c"], abs=1e-3)
        assert h1[cell]["chain"]["offset"] == pytest.approx(P.GATES["H1"][f"{key}_chain"], abs=1e-3)
    r = c["reading"]
    assert (r["B"], r["C"]) == pytest.approx((-3.170, -3.801), abs=2e-3)
    assert (r["R"], r["Ob"], r["I"]) == pytest.approx((5.264, 5.895, -0.447), abs=2e-3)
    assert r["R"] + r["Ob"] == pytest.approx(r["A"] - r["D"])
    assert (r["R_class"], r["Ob_class"]) == ("partial", "partial")
    assert len(r["overrides_and_remarks"]) == 2
    chain = c["attribution"]["H1"]["chain"]
    assert (chain["R_class"], chain["Ob_class"]) == ("carries", "does not carry")
    assert set(c["sensitivities"]) >= {"S2 eccentricity of p. ii", "S3 without § 7"}
    assert len(c["D"]["years"]) == 30


def test_s8_as_recorded():
    import json

    s = json.loads(P.OUT_JSON.read_text(encoding="utf-8"))["cells"]["sensitivities"]
    b = s["S8b the 1780 place on Table I's branch"]
    c = s["S8c the 1780 place: D interpolated from its neighbours"]
    assert b["delta_arcsec"]["1780"] == pytest.approx(3.402, abs=1e-3)
    assert c["delta_arcsec"]["1780"] == pytest.approx(2.797, abs=2e-3)
    assert b["H1"]["d1c"]["change_B"] == pytest.approx(4.81, abs=0.01)
    assert c["H1"]["d1c"]["change_B"] == pytest.approx(3.90, abs=0.01)


@pytest.mark.slow
def test_cells_b_and_c_recomputed_from_the_transcribed_inputs():
    """Cells B and C solved again from the transcribed tables, the memoir and the model,
    without the recorded JSON, against the recorded offsets."""
    import json

    got = P.cells(sensitivities=False)["cells"]
    want = json.loads(P.OUT_JSON.read_text(encoding="utf-8"))["cells"]["cells"]
    for hyp in ("H1", "H2"):
        for cell in ("B", "C"):
            for method in ("d1c", "chain"):
                assert got[hyp][cell][method]["offset"] == pytest.approx(
                    want[hyp][cell][method]["offset"], abs=1e-3
                ), (hyp, cell, method)
    assert got["H1"]["B"]["d1c"]["offset"] == pytest.approx(-3.170, abs=2e-3)
    assert got["H1"]["C"]["d1c"]["offset"] == pytest.approx(-3.801, abs=2e-3)
