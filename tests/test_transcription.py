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
    signed = _row("136", "11", "13", "12 18 19", dec_obs="-2 6 8.3", dec_calc_sec="-43.7",
                  dec_c_minus_o="+24.6")  # fmt: skip  # sign printed on the seconds
    assert tr.arithmetic_check(signed)["dec"] == pytest.approx(0.0, abs=1e-6)
    bad = {**r1690, "ra_c_minus_o": "-66.4"}  # a misread digit
    assert tr.arithmetic_check(bad)["ra"] == pytest.approx(0.5, abs=1e-6)


def _res(page, line, field, p1, p2, resolved):
    return {(page, line, field): dict(page=page, line=line, field=field, pass1=p1, pass2=p2,
                                      resolved=resolved, witness="w", reason="r")}  # fmt: skip


def test_reconcile_takes_logged_resolutions_and_refuses_gaps():
    a = [_row("129", "1", "23", "9 41 25", ra_obs="55 49 19.7", notes="x"),
         _row("129", "2", "2", "9 46 47", year="1780")]  # fmt: skip
    b = [_row("129", "1", "23", "9 41 25", ra_obs="55 49 19.7"),
         _row("129", "2", "2", "9 46 47", year="1789", notes="y")]  # fmt: skip
    log = _res("129", "2", "year", "1780", "1789", "1789")
    out = tr.reconcile(a, b, log)
    assert [r["year"] for r in out] == ["", "1789"]
    assert [r["notes"] for r in out] == ["p1: x", "p2: y"]
    with pytest.raises(ValueError, match="without a resolution"):
        tr.reconcile(a, b, {})
    with pytest.raises(ValueError, match="match no disagreement"):
        tr.reconcile(a, b, {**log, **_res("129", "1", "day", "1", "2", "1")})
    with pytest.raises(ValueError, match="log says"):
        tr.reconcile(a, b, _res("129", "2", "year", "1780", "1788", "1789"))


def test_committed_table77_is_the_reconciliation_of_the_committed_passes():
    """The reconciled file on disk equals reconcile(pass1, pass2, log): no hand edits."""
    a, b = (tr.read_pass(p) for p in tr.PASS_FILES)
    rows = tr.reconcile(a, b, tr.read_resolutions(tr.RESOLUTION_FILE))
    assert tr.read_pass(tr.RECONCILED_FILE) == rows
    assert len(rows) == 279
    assert sum(int(r["year"] or 0) < 1781 for r in rows if r["year"]) == 10  # ancient year groups
    assert [r["year"] for r in rows].index("1781") == 19  # 19 ancient rows, 260 modern


def test_resolution_log_and_misprint_table_parse_strictly():
    """Every log row has exactly its columns (a printed decimal comma must be quoted), and every
    misprint entry names the table's printed value and carries a reason."""
    table = tr.read_pass(tr.RECONCILED_FILE)
    log = tr.read_resolutions(tr.RESOLUTION_FILE)
    assert all(r["witness"].startswith("cdt_1849_additions.pdf p.") for r in log.values())
    m = tr.read_misprints(tr.MISPRINT_FILE, table)
    assert {(x["page"], x["line"], x["field"], x["action"]) for x in m} == {
        ("129", "10", "dec_c_minus_o", "log_only"),
        ("129", "12", "mean_time", "log_only"),
        ("129", "19", "dec_obs", "apply"),
        ("132", "5", "mean_time", "apply"),
        ("132", "21", "mean_time", "apply"),
        ("135", "3", "dec_obs", "apply"),
        ("135", "3", "ra_calc_sec", "log_only"),
    }
    assert all(x["uses_ephemeris"] == "no" for x in m)


def _digit_positions(cell: str) -> list[tuple[int, str]]:
    """(index, part) for each digit of a cell: part is top/min (leading sexagesimal parts) or
    sec_tens/sec_units/frac (the last part: seconds, or a seconds-only / C − O value)."""
    out, off = [], len(cell) - len(cell.lstrip("+-"))
    parts = cell.lstrip("+-").split()
    for pi, p in enumerate(parts):
        whole = p.partition(".")[0]
        for j, ch in enumerate(p):
            if ch.isdigit():
                if pi < len(parts) - 1:
                    part = "top" if pi == 0 else "min"
                elif j >= len(whole):
                    part = "frac"
                else:
                    part = "sec_tens" if len(whole) == 2 and j == 0 else "sec_units"
                out.append((off + j, part))
        off += len(p) + 1
    return out


def test_arithmetic_check_reach_on_every_single_digit_misread():
    """Pins the release README's statement: in every checkable row the printed C − O closes
    exactly; of all single-digit substitutions in the digits the check reads (seconds and
    decimals of the observed RA/Dec, the computed seconds, the C − O), it detects 98.5 %, and
    every miss is a tens-of-seconds digit off by exactly 6 (a change of 60 that the unprinted
    computed minute absorbs)."""
    rows = tr.read_pass(tr.RECONCILED_FILE)
    misprint = ("129", "10", "dec")  # 1764 Jan 15, in the misprint table
    fields = ("ra_obs", "ra_calc_sec", "ra_c_minus_o", "dec_obs", "dec_calc_sec", "dec_c_minus_o")
    n_checked = {"ra": 0, "dec": 0}
    total = caught = 0
    misses = set()
    for r in rows:
        base = tr.arithmetic_check(r)
        for coord, v in base.items():
            if v is not None:
                n_checked[coord] += 1
                if (r["page"], r["line"], coord) != misprint:
                    assert abs(v) < 1e-6, (r["page"], r["line"], coord, v)
        for f in fields:
            coord = f.split("_")[0]
            if base[coord] is None or (r["page"], r["line"], coord) == misprint:
                continue
            cell = r[f]
            for i, part in _digit_positions(cell):
                if f.endswith("_obs") and part in ("top", "min"):
                    continue  # degrees/hours/minutes: not read by the check, by design
                for d in "0123456789":
                    if d == cell[i]:
                        continue
                    v = tr.arithmetic_check({**r, f: cell[:i] + d + cell[i + 1 :]})[coord]
                    total += 1
                    if abs(v) > 1e-6:
                        caught += 1
                    else:
                        misses.add((part, abs(int(d) - int(cell[i]))))
    assert n_checked == {"ra": 279, "dec": 250}
    assert misses == {("sec_tens", 6)}
    assert caught / total == pytest.approx(0.985, abs=0.0005)
