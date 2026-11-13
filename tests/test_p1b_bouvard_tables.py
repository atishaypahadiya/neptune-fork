"""Step 3e: Bouvard's tables at Adams's 30 epochs (docs/notes/P1B_bouvard_tables.md)."""

import csv
from pathlib import Path

from fork.export import p1b_bouvard_tables as T
from fork.historical import bouvard_1821 as B

D = Path(__file__).resolve().parents[1] / "data" / "historical" / "adams_challis"
PASSES = D / "bouvard_1821" / "passes"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def test_inputs_have_the_hashes_of_the_plan():
    got = T.check_inputs(require_cache=False)
    assert set(got) == set(T.INPUT_HASHES)


def test_the_two_entries_agree_and_the_rows_are_their_reading():
    key = lambda r: (r["table"], r["row"], r["column"])  # noqa: E731
    a = {key(r): r["value"] for r in read(PASSES / "bouvard_1821_table_rows_route_b_pass1.csv")}
    b = {key(r): r["value"] for r in read(PASSES / "bouvard_1821_table_rows_route_b_pass2.csv")}
    assert len(a) == len(b) == 1855
    assert {k for k in a if a[k] != b[k]} == {("I", "1801", "VIII")}
    printed = B.route_b_rows(used=False)
    assert printed == {k: (v if v == b[k] else "150") for k, v in a.items()}
    # the two rows files do not overlap, so step 3d's instants read the rows they read before
    assert not set(B.table_rows()) & set(printed)
    assert len(B.all_table_rows()) == len(B.table_rows()) + len(printed)


def test_misprint_log_is_the_check_of_the_printed_differences():
    on_disk = read(B.ROUTE_B_MISPRINTS)
    assert on_disk == T.misprint_rows()
    assert [r["place"] for r in on_disk][:5] == [
        "Table IX, between rows 1790 and 1800, Différences",
        "Table X, between rows 84 and 85, Différences premières",
        "Table X, between rows 142 and 143, Différences premières",
        "Table XIII, row 720, Équation",
        "Table XIII, between rows 40 and 50, Différences",
    ]
    assert all(r["touches"] for r in on_disk)


def test_no_reading_of_the_new_rows_is_contradicted_on_both_sides():
    check = T.difference_check()
    readings = [f for f in check["findings"] if f["kind"] == "reading"]
    assert [(f["table"], f["row"], f["new"]) for f in readings] == [("XIII", 720, False)]
    # step 3d's reading stands, as step 3d ruled it
    assert B.route_b_rulings() == {}
    assert B.all_table_rows()["XIII", "720", "eq"] == "2,6"
    assert B.all_table_rows() == {**B.table_rows(), **B.route_b_rows(used=False)}
    # every new reading with a printed difference was checked on at least one side
    assert all(c["new_none"] == 0 for c in check["counts"].values())


# The tables' longitude (centesimal seconds) at the six instants of V1 and the seven of V1′,
# with the readings used and as printed, from the table path as it stood before
# table_longitude_at was added (commit df3d312).
RECORDED = {
    "1783": (1060041.2637539774, 1060032.62909073),
    "1798": (1839233.6160782645, 1839233.6160782645),
    "1813": (2621018.3559883614, 2621018.3559883614),
    "1828": (3336453.321454607, 3336453.321454607),
    "1846": (136265.68524779548, 136265.68524779548),
    "1715": (1879563.6619253624, 1879563.6619253624),
    "1789": (1365289.4925139016, 1365295.3415908809),
    "1804": (2157360.8093644893, 2157358.609934781),
    "1822": (3059152.151606344, 3059152.151606344),
    "1837": (3738788.38990761, 3738788.38990761),
    "1690": (662583.5428621583, 662583.5525372885),
    "1750": (3603644.5673999167, 3603644.5673999167),
    "1769": (402625.40601711813, 402625.40601711813),
}


def _t_of_midnight(y: int, m: int, d: int, frac: float = 0.0) -> float:
    import datetime

    return ((datetime.date(y, m, d) - datetime.date(1800, 1, 1)).days + frac) / 365.25


def test_v1_and_v1p_instants_unchanged():
    for inst in (*B.INSTANTS, *B.INSTANTS_V1P):
        used, printed = RECORDED[inst.name]
        assert abs(B.table_longitude(inst)["longitude"] - used) < 1e-6
        assert abs(B.table_longitude(inst, used=False)["longitude"] - printed) < 1e-6
        # the extended path, reading both rows files, at the same Paris midnight
        t = _t_of_midnight(*inst.date)
        at = B.table_longitude_at(t)
        assert abs(at["longitude"] - used) < 1e-6, inst.name
        assert abs(B.table_longitude_at(t, used=False)["longitude"] - printed) < 1e-6
        assert at["entries"] == B.table_longitude(inst)["entries"]
        assert (
            abs(at["argument_allowance_cs"] - B.table_longitude(inst)["argument_allowance_cs"])
            < 1e-12
        )
        assert (at["time_of_day"]["h"], at["time_of_day"]["m"], at["time_of_day"]["s"]) == (0, 0, 0)


def test_rounding_budget_on_constructed_dates():
    half = 0.05 * B.CENTESIMAL_SEC_IN_ARCSEC  # 0″.0162
    cases = [
        ((1783, 1, 1, 0.0), 17),  # Table I alone
        ((1798, 3, 7, 0.0), 19),  # a month and a day
        ((1798, 3, 1, 0.0), 18),  # the first of a month: no Table IV
        ((1798, 1, 7, 0.0), 18),  # January: no Table III
        ((1715, 3, 4, 0.0), 20),  # Table II as well
        ((1798, 3, 7, 0.5), 20),  # 5h 00m 00s: Table V only
        ((1798, 3, 7, 0.01233), 21),  # 0h 12m 33s: minutes and seconds of Table VI
        ((1798, 3, 7, 0.12445), 22),  # 1h 24m 45s
        ((1715, 3, 3, 0.62444605), 23),  # the 1715 place: 6h 24m 45s, a century, a month, a day
        ((1801, 3, 20, 0.00029), 20),  # 0h 00m 29s: the seconds alone
    ]
    for (y, m, d, frac), n in cases:
        ti = B.table_instant(_t_of_midnight(y, m, d, frac))
        assert B.rounding_entries(ti) == n, (y, m, d, frac)
    # the path counts the same entries where its rows were entered
    for (y, m, d, frac), n in (cases[0], cases[1], cases[4], cases[8], cases[9]):
        r = B.table_longitude_at(_t_of_midnight(y, m, d, frac))
        assert r["entries"] == n, (y, m, d, frac)
        assert abs(r["budget_arcsec"] - n * half) < 1e-12
    assert abs(17 * half - 0.275) < 0.001 and abs(23 * half - 0.373) < 0.001


def test_time_of_day_is_added_to_the_mean_longitude_only():
    t0 = _t_of_midnight(1798, 3, 7)
    a, b = B.table_longitude_at(t0), B.table_longitude_at(_t_of_midnight(1798, 3, 7, 0.12445))
    assert (b["time_of_day"]["h"], b["time_of_day"]["m"], b["time_of_day"]["s"]) == (1, 24, 45)
    rows = B.all_table_rows()
    motion = (
        sum(B.centesimal_seconds(rows[k]) for k in (("V", "1", "lon"), ("VI", "24", "lon")))
        + B.centesimal_seconds(rows["VI", "45", "lon"]) / 100
    )
    assert abs(b["mean_longitude"] - a["mean_longitude"] - motion) < 1e-9
    assert b["args"] == a["args"]


def test_a_missing_row_stops_the_epoch():
    import pytest

    with pytest.raises(B.MissingRow):
        B.table_longitude_at(_t_of_midnight(1798, 3, 7, 0.0222))  # VI row 22 was not entered


def test_the_1780_place_on_both_branches():
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A

    t = B.t_of_jd(R.place_epochs(A.load_memoir())[1780])
    a = B.table_longitude_at(t, branch="century")
    b = B.table_longitude_at(t, branch="own")
    assert a["epoch_rows"] == [("I", "1879", 1), ("II", "-100", 1), ("III", "Décembre common", 1),
                               ("IV", "18", 1)]  # fmt: skip
    assert b["epoch_rows"] == [("I", "1780", 1), ("IV", "14", -1)]
    assert (b["time_of_day"]["h"], b["time_of_day"]["m"], b["time_of_day"]["s"]) == (9, 12, 1)
    assert (a["entries"], b["entries"]) == (22, 21)


def test_the_rule_of_section_2_on_constructed_readings():
    # read 3,0 between 1,0 and 3,0 with differences 1,0 and 1,0: both sides give 2,0
    r = T.rule_of_section_2(30, 10, 10, 30, 10)
    assert r["contradicted_on_both_sides"] and r["replace_with"] == 20
    # the two sides give different values: contradicted, left as read
    r = T.rule_of_section_2(30, 10, 10, 50, 10)
    assert r["contradicted_on_both_sides"] and r["replace_with"] is None
    # one side only: not replaced
    r = T.rule_of_section_2(30, 10, 20, 50, 10)
    assert (r["below"], r["above"], r["replace_with"]) == (True, False, None)
    # a side whose row was not entered is not checked
    r = T.rule_of_section_2(30, None, 10, 40, 10)
    assert (r["below"], r["above"], r["contradicted_on_both_sides"]) == (None, True, False)
    # the parser of the log reads a replaced reading back
    row = {"place": "Table XII, row 2600, Équation", "status": B.REPLACED + ": both",
           "used": "3.31,6"}  # fmt: skip
    table, rw, name = (p.strip() for p in row["place"].split(","))
    assert (table.split()[1], rw.split()[1], name) == ("XII", "2600", B.COLUMN_NAMES["eq"])


def _state(d1c_cells, chain_cells):
    from fork.export import p1b_bouvard as P

    return T.reading_state(P.attribute(*d1c_cells), P.attribute(*chain_cells))


def test_same_reading_on_constructed_effects():
    # A and D fixed at the gates' values; B and C chosen to give the effects wanted
    both = (1.6, -3.2, -3.8, -9.5)  # R +5.3, Ob +5.9: partial, partial ("Both")
    obs = (1.6, 1.0, -8.0, -9.5)  # R +0.6, Ob +10.6: the observations carry it
    ref = (1.6, -8.5, 1.0, -9.5)  # R +10.3, Ob +0.6: the reference carries it
    # the same rows by both methods, in both: the same reading
    r = T.same_reading(_state(obs, obs), _state(obs, (1.0, 0.5, -8.5, -9.9)))
    assert r["same"] and r["differs_in"] == []
    # D1(c) in another row: condition 1
    r = T.same_reading(_state(obs, obs), _state(both, both))
    assert not r["same"] and r["differs_in"] == [
        "1: D1(c)'s row",
        "2: his chain's row (the override applies in neither)",
    ]
    # the chain alone moves to another row, D1(c) stays: condition 2 only
    r = T.same_reading(_state(obs, ref), _state(obs, obs))
    assert r["differs_in"] == ["2: depends on the method"]
    # "depends on the method" in one and not the other: condition 2
    r = T.same_reading(_state(both, both), _state(both, ref))
    assert r["differs_in"] == ["2: depends on the method"]
    # depends on the method in both: the same, whatever the chain's row
    r = T.same_reading(_state(both, ref), _state(both, obs))
    assert r["same"]
    # the interaction override in one only: condition 3
    big_i = (1.6, -9.0, -9.0, -9.5)  # R +5.55, Ob +5.55 (partial, partial), I +5.05
    r = T.same_reading(_state(big_i, big_i), _state(both, both))
    assert r["differs_in"] == ["3: the interaction override"]


def test_r1_rule_on_constructed_differences():
    half = 0.3 / T.CS  # budget 0″.3 in centesimal seconds
    allow = 1.0 / T.CS
    assert T.r1_outcome(0.2 / T.CS, half, allow) == T.OUTCOMES[0]
    assert T.r1_outcome(-1.2 / T.CS, half, allow) == T.OUTCOMES[1]
    assert T.r1_outcome(1.5 / T.CS, half, allow) == T.OUTCOMES[2]
    assert T.r1_outcome(-2.1 / T.CS, half, allow) == T.OUTCOMES[3]


def test_the_two_row_apparatus():
    d = B.all_table_rows()
    # Table VII: the printed entries at whole seconds and years, nothing for no second
    # difference, and between rows and columns
    assert abs(B._table_vii(d, 5.0, 5.0) - 0.6) < 1e-12
    assert abs(B._table_vii(d, 4.0, 1.0) - 0.2) < 1e-12
    assert abs(B._table_vii(d, 4.0, 9.0) - 0.2) < 1e-12  # the column of 1 and 9 years
    assert B._table_vii(d, 0.0, 5.0) == 0.0 and B._table_vii(d, 3.0, 0.0) == 0.0
    assert abs(B._table_vii(d, 4.5, 5.0) - 0.55) < 1e-12
    # Table IX: the sign printed at the head of a run is carried down it
    s = B._signed_ix_differences(d)
    assert (s[1690], s[1730], s[1740], s[1790], s[1860]) == (19.4, 1.0, -3.6, -25.4, -42.1)
    # Table VIII: a printed date, and a date between two
    assert abs(B._table_viii(d, (1782, 12, 31)) - 1.0) < 1e-12
    assert abs(B._table_viii(d, (1834, 8, 13)) - (0.60 + 0.03 * 5 / 10)) < 1e-12
    # Table X between two rows at a whole degree is the row itself
    assert abs(B._centre_two_row(d, 84.0) - B._centre_row(d, 84)[0]) < 1e-9


def _recorded() -> dict:
    import json

    return json.loads(T.OUT_JSON.read_text(encoding="utf-8"))


def test_r1_as_recorded():
    """R1 is table lookups: recomputed and set against the recorded run."""
    got, rec = T.r1(), _recorded()["R1"]
    for y, row in rec["epochs"].items():
        assert (
            abs(
                got["epochs"][y]["formula_minus_tables_arcsec"] - row["formula_minus_tables_arcsec"]
            )
            < 1e-9
        )
        assert got["epochs"][y]["outcome"] == row["outcome"]
    assert len(got["epochs"]) == 29 and "1780" not in got["epochs"]
    counts = [got["all_29"]["counts"][o] for o in T.OUTCOMES]
    assert counts == [18, 11, 0, 0]
    assert [got["fresh_17"]["counts"][o] for o in T.OUTCOMES] == [12, 5, 0, 0]
    assert [got["near_v1_or_v1p_12"]["counts"][o] for o in T.OUTCOMES] == [6, 6, 0, 0]
    assert got["all_29"]["largest"]["epoch"] == "1690"
    assert abs(got["all_29"]["largest"]["arcsec"] + 1.073) < 0.001
    assert got["fresh_17"]["largest"]["epoch"] == "1819"
    assert got["all_29"]["outside"] == [] and got["stop_and_ask"] == []


def test_gates_and_r2_as_recorded():
    from fork.export import p1b_bouvard as P

    j = _recorded()
    assert j["gates_pass"] and len(j["gates"]) == 8
    assert j["step_3d_columns_reproduced_max_abs_arcsec"] < 1e-9
    r = j["R2"]["reading"]
    assert abs(r["R"] - 3.16) < 0.005 and abs(r["Ob"] - 7.99) < 0.005 and abs(r["I"] + 1.04) < 0.005
    assert abs(r["R"] + r["Ob"] - (r["A"] - r["D"])) < 1e-9
    assert (r["R_class"], r["Ob_class"]) == ("partial", "carries")
    assert r["row"] == P.READING_ROWS["partial", "carries"]
    assert r["chain_row"] == r["row"]
    assert not r["state"]["depends_on_the_method"] and not r["state"]["interaction_override"]


def test_r3_as_recorded():
    j = _recorded()["R3"]
    primary, century, own = j.values()
    assert not primary["same"] and primary["differs_in"][0] == "1: D1(c)'s row"
    assert primary["formula_state"]["d1c_classes"] == ["does not carry", "carries"]
    assert not century["same"] and century["differs_in"] == ["1: D1(c)'s row"]
    assert own["same"]
    # the rule applied again to the recorded states gives the recorded outcome
    for v in j.values():
        assert T.same_reading(v["tables_state"], v["formula_state"])["same"] == v["same"]


def test_s_midnight_not_run_for_want_of_rows():
    import numpy as np

    ep = T.epochs()
    s = T.s_midnight_status({"years": list(ep), "t": np.array(list(ep.values()))})
    assert not s["run"]
    assert sorted(s["missing"]) == ["1764", "1792", "1810", "1819", "1840"]
    assert s == _recorded()["sensitivities"]["S-midnight"]
