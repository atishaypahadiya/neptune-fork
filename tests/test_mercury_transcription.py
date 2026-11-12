"""Newcomb 1882, pp. 451–452 (data/historical/mercury/): the two entries agree, and the
checks the table itself allows. Nothing is fitted to these data."""

import csv
from collections import defaultdict
from pathlib import Path

M = Path(__file__).resolve().parents[1] / "data" / "historical" / "mercury"
FIELDS = ("table", "row", "contact", "year", "month", "day", "obs_h", "obs_m", "obs_s",
          "k_sign", "k_coeff", "weight")  # fmt: skip


def read(name: str) -> list[dict[str, str]]:
    with open(M / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


FINAL = read("newcomb_1882_contacts.csv")


def seconds(r: dict[str, str]) -> float:
    """Seconds from the noon that begins the row's astronomical day."""
    return int(r["obs_h"]) * 3600 + int(r["obs_m"]) * 60 + float(r["obs_s"])


def test_the_two_entries_agree_and_the_final_table_is_theirs():
    p1 = read("passes/newcomb_1882_contacts_pass1.csv")
    p2 = read("passes/newcomb_1882_contacts_pass2.csv")
    assert len(p1) == len(p2) == len(FINAL) == 58
    for a, b, c in zip(p1, p2, FINAL, strict=True):
        for k in FIELDS:
            assert a[k] == b[k] == c[k], (a["table"], a["row"], k)


def test_rows_per_table_and_ranges():
    n = defaultdict(int)
    for r in FINAL:
        n[r["table"]] += 1
        assert 0 <= int(r["obs_h"]) < 24 and 0 <= int(r["obs_m"]) < 60
        assert 0 <= float(r["obs_s"]) < 60
        assert 0 < float(r["weight"]) <= 6.0  # p. 451: "6 is the maximum weight"
        assert r["month"] == ("Nov" if r["table"] in ("I", "III") else "May")
        assert (r["contact"] in ("II", "III")) == (r["table"] in ("I", "II"))
    assert dict(n) == {"I": 25, "II": 12, "III": 15, "IV": 6}


def test_one_k_term_per_transit_across_the_four_tables():
    k = defaultdict(set)
    for r in FINAL:
        k[r["year"]].add((r["k_sign"] or "-", r["k_coeff"]) if r["k_coeff"] != "0" else "0")
    assert all(len(v) == 1 for v in k.values()), {y: v for y, v in k.items() if len(v) > 1}


def test_contacts_are_in_order_within_a_transit():
    """Within a year: the first interior contact precedes the second by 1 to 8 hours, and
    the exterior contact follows the last interior one by 1 to 8 minutes. Days as printed are
    not used (two rows are in the misprint log); times are taken modulo a day."""
    interior = defaultdict(list)
    for r in FINAL:
        if r["table"] in ("I", "II"):
            interior[r["year"]].append(seconds(r))
    for year, t in interior.items():
        if len(t) == 2:
            assert 3600 <= (t[1] - t[0]) % 86400 <= 8 * 3600, year
    for r in FINAL:
        if r["table"] in ("III", "IV") and r["year"] in interior:
            gap = (seconds(r) - interior[r["year"]][-1]) % 86400
            assert 60 <= gap <= 480, (r["year"], gap)


def test_logged_misprints_point_at_printed_values():
    for m in read("newcomb_1882_misprints.csv"):
        row = next(r for r in FINAL if (r["table"], r["row"]) == (m["table"], m["row"]))
        assert row[m["field"]] == m["printed"]
        assert m["action"] == "log_only" and m["corrected"] == ""
