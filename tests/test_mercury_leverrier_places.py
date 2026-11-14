"""Le Verrier 1859, Annales V, pp. 51-61, and Newcomb 1882, p. 450 (data/historical/mercury/):
the two entries agree, the final files are theirs, and the checks the table allows. Nothing
is fitted to these data, and the table's error columns are not in the repository."""

import csv
from collections import Counter
from pathlib import Path

M = Path(__file__).resolve().parents[1] / "data" / "historical" / "mercury"
VALUES = ("page", "row", "group", "year", "year_printed", "month", "month_printed", "day",
          "time_h", "time_m", "time_s", "lon_d", "lon_m", "lon_s", "lat_sign", "lat_d", "lat_m",
          "lat_s")  # fmt: skip
ROWS_PER_PAGE = {51: 17, 52: 39, 53: 39, 54: 39, 55: 41, 56: 41, 57: 38, 58: 39, 59: 40,
                 60: 42, 61: 22}  # fmt: skip
MONTHS = {"Janv.", "Fév.", "Mars", "Avril", "Mai", "Juin", "Juill.", "Août", "Sept.", "Oct.",
          "Nov.", "Déc."}  # fmt: skip


def read(name: str) -> list[dict[str, str]]:
    with open(M / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


FINAL = read("leverrier_1859_paris_places.csv")


def test_the_two_entries_agree_and_the_final_table_is_theirs():
    p1 = read("passes/leverrier_1859_paris_places_pass1.csv")
    p2 = read("passes/leverrier_1859_paris_places_pass2.csv")
    assert len(p1) == len(p2) == len(FINAL) == 397
    for a, b, c in zip(p1, p2, FINAL, strict=True):
        for k in VALUES:
            assert a[k].strip() == b[k].strip() == c[k], (a["page"], a["row"], k)


def test_no_error_column_is_present():
    assert not any("err" in k.lower() or "excess" in k.lower() for k in FINAL[0])


def test_rows_per_page_and_ranges():
    assert Counter(int(r["page"]) for r in FINAL) == ROWS_PER_PAGE
    years = [int(r["year"]) for r in FINAL]
    assert years == sorted(years)
    assert set(years) <= set(range(1801, 1829)) | set(range(1836, 1843))
    assert sum(y <= 1828 for y in years) == 240 and sum(y >= 1836 for y in years) == 157
    for r in FINAL:
        assert r["month"] in MONTHS and 1 <= int(r["day"]) <= 31
        assert 0 <= int(r["time_h"]) < 24 and 0 <= int(r["time_m"]) < 60
        if r["lon_d"]:
            assert 0 <= int(r["lon_d"]) < 360 and 0 <= int(r["lon_m"]) < 60
            assert 0 <= float(r["lon_s"] or 0) < 60
        if r["lat_d"]:
            assert 0 <= int(r["lat_d"]) <= 7 and 0 <= float(r["lat_s"] or 0) < 60
            assert r["lat_sign"] in ("", "-")


def test_group_numbers_run_1_to_195_with_one_logged_exception():
    g = [int(r["group"]) for r in FINAL]
    assert g[0] == 1 and g[-1] == 195
    drops = [(r["page"], r["row"]) for r, a, b in zip(FINAL[1:], g, g[1:], strict=False) if b < a]
    assert drops == [("58", "9")]
    assert sorted(set(g)) == [x for x in range(1, 196) if x != 124]  # 24 printed for 124


def test_blank_cells_are_the_guillemets():
    for r in FINAL:
        blank = [k for k in ("lon_d", "lat_d") if not r[k]]
        if blank:
            assert "guillemet" in r["flags"], (r["page"], r["row"])
    assert sum("guillemet" in r["flags"] for r in FINAL) == 9


def test_misprint_log_names_printed_values_that_are_kept():
    log = read("leverrier_1859_paris_places_misprints.csv")
    assert {(m["page"], m["row"]) for m in log} == {("58", "9"), ("57", "15;16")}
    row = {(r["page"], r["row"]): r for r in FINAL}
    assert row[("58", "9")]["group"] == "24"
    assert (row[("57", "15")]["lon_d"], row[("57", "16")]["lon_d"]) == ("163", "154")
    assert all(m["action"] == "log_only" and m["corrected"] == "" for m in log)


def test_p450_final_text_is_both_entries_after_reconciliation():
    final = (M / "newcomb_1882_p450.txt").read_text(encoding="utf-8")
    final.encode("ascii")
    assert "[?]" not in final
    for name in ("passes/newcomb_1882_p450_pass1.txt", "passes/newcomb_1882_p450_pass2.txt"):
        body = (M / name).read_text(encoding="utf-8").split("\nNOTES:")[0]
        body = body.replace("mo[?]tion", "motion").replace("18s[?]", "18s")
        assert body.replace("1677--1881", "1677-1881").rstrip("\n") == final.rstrip("\n")
    assert "- 18s + 60s T^2" in final and "kDelta_t" in final


def test_time_increases_and_daily_steps_are_possible_but_for_the_logged_pair():
    """Within 6 days, Mercury moves under 2.3 degrees a day in geocentric longitude and
    under 0.45 degrees a day in latitude. The only row breaking this is in the misprint log."""
    import datetime as dt

    months = {m: i for i, m in enumerate(MONTHS_IN_ORDER, 1)}
    prev, bad = None, []
    for r in FINAL:
        if not (r["lon_d"] and r["lat_d"]):
            prev = None
            continue
        t = dt.datetime(int(r["year"]), months[r["month"]], int(r["day"]))
        t += dt.timedelta(hours=int(r["time_h"]), minutes=int(r["time_m"]))
        lon = int(r["lon_d"]) + int(r["lon_m"]) / 60 + float(r["lon_s"] or 0) / 3600
        lat = int(r["lat_d"]) + int(r["lat_m"]) / 60 + float(r["lat_s"] or 0) / 3600
        lat = -lat if r["lat_sign"] == "-" else lat
        if prev:
            days = (t - prev[0]).total_seconds() / 86400
            assert days > 0, (r["page"], r["row"])
            step = (lon - prev[1] + 180) % 360 - 180
            if days < 6 and (
                abs(step) > 2.3 * days + 0.1 or abs(lat - prev[2]) > 0.45 * days + 0.05
            ):
                bad.append((r["page"], r["row"]))
        prev = (t, lon, lat)
    assert bad == [("57", "16")]


MONTHS_IN_ORDER = ("Janv.", "Fév.", "Mars", "Avril", "Mai", "Juin", "Juill.", "Août", "Sept.",
                   "Oct.", "Nov.", "Déc.")  # fmt: skip


def test_the_spot_check_by_a_person_matches_the_final_table():
    """20 cells checked against the page images (docs/review/, 13 November 2026): all match."""
    import hashlib
    import json

    path = M.parents[2] / "docs" / "review" / "leverrier_places_spotcheck_human_results.json"
    sc = json.loads(path.read_text(encoding="utf-8"))
    table = M / "leverrier_1859_paris_places.csv"
    assert sc["places_sha256"] == hashlib.sha256(table.read_bytes()).hexdigest()
    assert sc["reviewer"] == "human" and sc["seed"] == 20261114
    cells = sc["cells"]
    assert len(cells) == 20 and all(c["matches"] is True and c["reading"] == "" for c in cells)
    rows = {(r["page"], r["row"]): r for r in FINAL}
    keys = [(c["page"], c["row"], c["field"]) for c in cells]
    assert len(set(keys)) == 20
    assert all(rows[k[:2]][k[2]] == c["value"] for k, c in zip(keys, cells, strict=True))
    assert Counter(c["reason"] for c in cells) == {
        "random": 10, "marked doubtful": 6, "misprint log": 3, "guillemet cell": 1
    }  # fmt: skip
