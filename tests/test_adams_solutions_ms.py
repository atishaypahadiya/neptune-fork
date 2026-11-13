"""Manuscript check of Adams's solutions (St John's College, MS W16): the two entries, the
sheets' own arithmetic, and the comparison with the committed tables. No image is read."""

import csv
from pathlib import Path

import pytest

D = Path(__file__).resolve().parents[1] / "data" / "historical" / "adams_challis"


def read(name: str) -> list[dict[str, str]]:
    with (D / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


P1, P2 = read("passes_ms/stjohns_w16_pass1.csv"), read("passes_ms/stjohns_w16_pass2.csv")
REC = read("adams_solutions_ms_reconciliation.csv")


def block(rows, image, name):
    return {r["label"]: r["value"] for r in rows if r["image"] == image and r["block"] == name}


def test_entries_differ_only_where_the_reconciliation_says():
    assert len(P1) == len(P2) == 200
    diff = {(a["image"], a["block"], a["label"]) for a, b in zip(P1, P2, strict=True)
            if a["value"] != b["value"]}  # fmt: skip
    logged = {(r["image"], r["block"], r["label"]) for r in REC if r["pass1"] != r["pass2"]}
    assert diff == logged and len(diff) == 4


def test_the_sheets_own_sums():
    b = block(P1, "16870", "first approximation")
    six = sum(float(b[f"left column term {i}"]) for i in range(1, 6))
    six += float(b["left column term 6 (below these, above the rule)"])
    last = float(b["below the rule, line 3 (underlined)"])
    assert six == pytest.approx(last, abs=5e-5)
    assert float(b["below the rule, line 1"]) + float(b["below the rule, line 2"]) == pytest.approx(
        last, abs=5e-5
    )
    for name, n_left in (("theta = -45 sums", 5), ("theta = -43 sums", 6)):
        b = block(P1, "16868", name)
        left = sum(float(b[f"left list {i}"]) for i in range(1, n_left + 1))
        assert left == pytest.approx(float(b["left list total"]), abs=5e-4)
        right = sum(float(b[f"right list {i}"]) for i in range(1, 7))
        assert right == pytest.approx(float(b["right list total"]), abs=5e-4)
        assert float(b["left list total"]) + float(b["right list total"]) == pytest.approx(
            float(b["whole sum"]), abs=5e-4
        )
    b43 = block(P2, "16871", "column 1846")  # the second entry's −23.27 is the one that sums
    assert sum(float(b43[f"term {i}"]) for i in range(1, 11)) == pytest.approx(
        float(b43["total"]), abs=5e-3
    )
    b = block(P1, "16871", "column 1843")
    assert sum(float(b[f"term {i}"]) for i in range(1, 11)) == pytest.approx(
        float(b["total"]), abs=5e-3
    )
    for year in ("1834", "1840", "1846"):
        b = block(P1, "16869", f"column {year}")
        terms = sum(float(b[f"term {i}"]) for i in range(1, 12))
        assert terms == pytest.approx(
            float(b["sum of terms"]), abs=1.5e-6
        )  # 1840, 1846: one unit off


def test_theta_of_august_1846_is_where_the_two_sums_cross_zero():
    s45 = float(block(P1, "16868", "theta = -45 sums")["whole sum"])
    s43 = float(block(P1, "16868", "theta = -43 sums")["whole sum"])
    theta = -43.0 + 2.0 * (-s43) / (s43 - s45)
    assert theta == pytest.approx(-(42 + 51 / 60), abs=0.01)


def test_manuscript_values_against_the_committed_tables():
    ms = read("adams_solutions_ms.csv")
    assert [r["comparison"] for r in ms].count("agrees") == 11
    assert not [r for r in ms if r["comparison"] == "contradicts"]
    s4 = next(r for r in read("adams_solutions.csv") if r["solution_id"] == "S4")
    assert "0.1428" in s4["eccentricity"] and "0.000173" in s4["mass_sun"]
    assert "321° 40′" in s4["stated_heliocentric_longitude"]
    a16e = next(r for r in read("adams_solutions_statements.csv") if r["statement_id"] == "A16e")
    for r in ms:
        if r["image"] == "16869":
            assert r["value_ms"].lstrip("+").replace(".", "·") in a16e["other_values"]


def test_the_authors_spot_check_covers_the_reconciled_entries():
    import collections
    import hashlib
    import json

    f = D.parents[2] / "docs" / "review" / "stjohns_spotcheck_human_results.json"
    sc = json.loads(f.read_text(encoding="utf-8"))
    p1 = D / "passes_ms" / "stjohns_w16_pass1.csv"
    assert sc["pass1_sha256"] == hashlib.sha256(p1.read_bytes()).hexdigest()
    assert sc["reviewer"] == "human" and sc["seed"] == 20261112
    cells = sc["cells"]
    assert len(cells) == 20 and all(c["matches"] is True and c["reading"] == "" for c in cells)
    # The kit shows the reconciled table: the first entry, with the reconciliation applied.
    first = {(r["image"], r["block"], r["label"]): r["value"] for r in P1}
    first |= {
        (r["image"], r["block"], r["label"]): r["resolved"]
        for r in REC
        if r["block"] != "letter head"
    }
    keys = [(c["image"], c["block"], c["label"]) for c in cells]
    assert len(set(keys)) == 20
    assert all(first[k] == c["value"] for k, c in zip(keys, cells, strict=True))
    how = collections.Counter(c["reason"] for c in cells)
    assert how == {"random": 9, "one per image": 4, "flagged doubtful": 4, "entries differed": 3}
    # Every difference in a value is among the checked cells, and none is left doubtful.
    differed = {k for k, c in zip(keys, cells, strict=True) if c["reason"] == "entries differed"}
    values = {(r["image"], r["block"], r["label"]) for r in REC if r["pass1"] != r["pass2"]}
    assert differed == {k for k in values if k[1] != "letter head"}
    assert all(r["resolved"] != "doubtful" for r in REC)
    assert all(
        first[k] == r["resolved"]
        for r in REC
        if (k := (r["image"], r["block"], r["label"])) in differed
    )
    assert 1 - 0.025 ** (1 / 9) == pytest.approx(0.336, abs=0.001)
