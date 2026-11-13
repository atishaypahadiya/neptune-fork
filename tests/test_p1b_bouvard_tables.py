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
