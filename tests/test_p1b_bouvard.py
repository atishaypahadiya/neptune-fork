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
