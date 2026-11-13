"""Phase 1b, step 3e: the four cells of step 3d with Bouvard's tables themselves (route B).

The plan is docs/notes/P1B_bouvard_tables.md, frozen before any table longitude was formed
at the 30 epochs. Its order (§ 8): the input hashes; the readings of the new rows against
their printed differences, logged before any longitude is formed; the table path at the
epochs; R1 (tables against the corrected formula); the gates; R2 (the four cells); R3 (the
comparison with step 3d); the sensitivities of § 6.

    uv run python -m fork.export.p1b_bouvard_tables --misprints   # the misprint log
    uv run python -m fork.export.p1b_bouvard_tables   # docs/notes/P1B_bouvard_tables.json
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

from fork.historical import bouvard_1821 as B

ROOT = Path(__file__).resolve().parents[3]
OUT_JSON = ROOT / "docs" / "notes" / "P1B_bouvard_tables.json"
_BV = "data/historical/adams_challis/bouvard_1821"

# § 1 of the plan: the inputs and their SHA-256
INPUT_HASHES = {
    "data/cache/adams_challis/bouvard/bouvard_1821_SfZ4Eiw9yR0C.pdf":
        "ebd3d9fbe1d44a9cc23eeb23225e15f44b226ca29b20119df0bd0c4ea4f6d1f6",
    f"{_BV}/bouvard_1821_uranus_table_rows.csv":
        "5f8919c5477e3bf959828697f3a0db876bf0ad376a44b2c7a512a6f55e4272c5",
    f"{_BV}/bouvard_1821_misprints.csv":
        "459eb0914691df29e9acaae265b2d85db94413ee91da0c196424fec0ae9c0a89",
    f"{_BV}/bouvard_1821_uranus_formula.csv":
        "42ceb111a5d0aece67cd43a7d57e5f0a46aaac889c15aa7b35a160ef2c7be1d3",
    f"{_BV}/bouvard_1821_great_inequalities.csv":
        "01b494d839ac2fa6079da5d18d3a08fc59302ebfdad8c723918e3c180747f25b",
    "data/historical/adams_challis/adams_memoir_normal_places.csv":
        "52aa0a9cb1459ab7181991e7441f22c303e87a8fbf4166564a8e6d1481f12963",
    "data/historical/adams_challis/adams_memoir_s7_terms.csv":
        "53339c28a8e9f9ce223fb5d92a2699f94720d36ca14298eb2ce5485cbf1888e6",
    "docs/notes/P1B_bouvard_attribution.json":
        "3a8f7ea23449375110fe213d590575cfd0d08f196570b86043e2a08f851f51f3",
    f"{_BV}/bouvard_1821_uranus_table_rows_route_b.csv":
        "6b2efe05cf1f7844755bbc1e132ce441d204a2a1779c1ba8cfff1ee4d0cb07c1",
    f"{_BV}/bouvard_1821_route_b_reconciliation.csv":
        "3fd4b42239d77cdcee976a2a572024febd011e444c30d82533c24956ccf5d098",
    f"{_BV}/passes/bouvard_1821_table_rows_route_b_pass1.csv":
        "540984f284ba7d1e2c5578c68843edf54fffebfa026bc908415e689940b1be76",
    f"{_BV}/passes/bouvard_1821_table_rows_route_b_pass2.csv":
        "faea3e0e47f3b38adf7f6ea54b0ba6192ed76fcb4f96507eea5c8285544225e8",
}  # fmt: skip
# the scan is cached, not committed: it is checked when present (the runner) and skipped by
# the tests, which run where the cache is absent
CACHED_ONLY = ("data/cache/adams_challis/bouvard/bouvard_1821_SfZ4Eiw9yR0C.pdf",)


class InputChanged(RuntimeError):
    """An input of § 1 does not have its recorded hash: the runner refuses to run."""


def check_inputs(require_cache: bool = True) -> dict[str, str]:
    """§ 1: every input's SHA-256 against the plan. Raises ``InputChanged`` on a difference."""
    got = {}
    for rel, want in INPUT_HASHES.items():
        path = ROOT / rel
        if rel in CACHED_ONLY and not require_cache and not path.exists():
            got[rel] = "not present (cache)"
            continue
        h = hashlib.sha256(path.read_bytes()).hexdigest()
        if h != want:
            raise InputChanged(f"{rel}: {h} differs from the plan's {want}")
        got[rel] = h
    return got


# ------------------------------------------------------------------ the printed differences

# table -> (column of the reading, column of its first difference, step of the argument,
#           period of the argument or None)
DIFFERENCED = {
    "IX": ("eq", "diff", 10, None),
    "X": ("eq", "d1", 1, 400),
    "XI": ("eq", "diff", 10, 1000),
    "XII": ("eq", "diff", 100, 10_000),
    "XIII": ("eq", "diff", 10, 1000),
    "XIV": ("eq", "diff", 10, 1000),
}


def _tenths(cs: float) -> int:
    return round(cs * 10)


def _fmt(tenths: int, like: str) -> str:
    """Tenths of a centesimal second written as the table prints them (m.s,t or s,t)."""
    v = abs(tenths)
    sec, t = divmod(v, 10)
    m, s = divmod(sec, 100)
    sign = "-" if tenths < 0 and like.startswith("-") else ""
    if m or like.count(".") >= 1:
        return f"{sign}{m}.{s:02d},{t}"
    return f"{sign}{s},{t}"


def _fmt_centre(tenths: int, row: int) -> str:
    """A Table X equation (tenths of a centesimal second, constant subtracted as printed)
    written as the line prints it."""
    v = tenths  # a degree is 100,000 tenths of a centesimal second
    if v < 0:
        v += 40_000_000 if row <= 200 else 1_000_000  # 39x printed in full, or its last figure
    deg, rest = divmod(v, 100_000)
    m, rest = divmod(rest, 1000)
    s, t = divmod(rest, 10)
    return f"{deg}.{m:02d}.{s:02d},{t}"


def _label(d: dict, table: str, row: int, col: str) -> int | None:
    """The label under which a line was entered: the line of argument 0 and the line of the
    full circle are one line (Table X 0 and 400; XI 0 and 1000; XII 0 and 10000)."""
    period = DIFFERENCED[table][3]
    cands = [row]
    if period is not None:
        r = row % period
        cands += [r] + ([period] if r == 0 else [])
    for lab in cands:
        if d.get((table, str(lab), col)):
            return lab
    return None


def _reading(d: dict, table: str, row: int) -> int | None:
    lab = _label(d, table, row, "eq")
    if lab is None:
        return None
    if table == "X":
        return _tenths(B._centre_row(d, lab)[0])
    return _tenths(B.centesimal_seconds(d[table, str(lab), "eq"]))


def _difference(d: dict, table: str, row: int, col: str) -> int | None:
    lab = _label(d, table, row, col)
    if lab is None:
        return None
    return _tenths(B.centesimal_seconds(d[table, str(lab), col].lstrip("+-")))


def rule_of_section_2(
    e: int, e_lo: int | None, d_lo: int | None, e_hi: int | None, d_hi: int | None
) -> dict:
    """§ 2 (step 3d's first amendment) on one reading, all in tenths of a centesimal second:
    is it contradicted by the printed differences on both sides, and if so do both give the
    same value? Magnitudes only (a sign is printed on the first line of a run alone)."""
    below = None if e_lo is None or d_lo is None else abs(e - e_lo) == d_lo
    above = None if e_hi is None or d_hi is None else abs(e_hi - e) == d_hi
    both = below is False and above is False
    common = sorted({e_lo + d_lo, e_lo - d_lo} & {e_hi + d_hi, e_hi - d_hi}) if both else []
    return {
        "below": below,
        "above": above,
        "contradicted_on_both_sides": both,
        "by_the_differences": common,
        "replace_with": common[0] if both and len(common) == 1 else None,
    }


def _rows_of(d: dict, table: str, col: str) -> list[int]:
    return sorted({int(r) for (t, r, c) in d if t == table and c == col})


def difference_check() -> dict:
    """Every reading with a printed difference beside it, against those differences
    (magnitudes: the tables print a sign only on the first line of a run). A reading both
    of whose differences contradict it is replaced only when both give the same value (the
    rule of § 2). Printed differences that the readings at both their ends contradict, and
    Table X's second differences that its first differences contradict, are reported too.
    No longitude is formed."""
    new = B.route_b_rows(used=False)
    d = {**B.table_rows(used=True), **new}  # neighbours as step 3d uses them
    counts: dict[str, dict[str, int]] = {}
    findings: list[dict] = []
    for table, (ec, dc, step, _) in DIFFERENCED.items():
        rows = _rows_of(d, table, ec)
        c = counts.setdefault(table, {"new_both": 0, "new_one": 0, "new_none": 0,
                                      "old_checked": 0})  # fmt: skip

        for k in rows:
            e = _reading(d, table, k)
            lo, hi = k - step, k + step
            e_lo, e_hi = _reading(d, table, lo), _reading(d, table, hi)
            d_lo, d_hi = _difference(d, table, lo, dc), _difference(d, table, k, dc)
            rule = rule_of_section_2(e, e_lo, d_lo, e_hi, d_hi)
            is_new = (table, str(k), ec) in new
            n_sides = (rule["below"] is not None) + (rule["above"] is not None)
            if is_new:
                c[["new_none", "new_one", "new_both"][n_sides]] += 1
            elif n_sides:
                c["old_checked"] += 1
            if rule["contradicted_on_both_sides"]:
                findings.append({
                    "kind": "reading",
                    "table": table, "row": k, "column": ec, "new": is_new,
                    "reads": d[table, str(k), ec],
                    "neighbours": (lo, e_lo, hi, e_hi), "differences": (d_lo, d_hi),
                    "by_the_differences": rule["by_the_differences"],
                    "same_value": rule["replace_with"] is not None,
                })  # fmt: skip
        # printed first differences against the readings at both their ends
        for k in _rows_of(d, table, dc):
            hi = k + step
            e0, e1, dv = _reading(d, table, k), _reading(d, table, hi), _difference(d, table, k, dc)
            if e0 is None or e1 is None or dv is None or abs(e1 - e0) == dv:
                continue
            findings.append({
                "kind": "first difference",
                "table": table, "row": k, "column": dc,
                "new": (table, str(k), dc) in new,
                "reads": d[table, str(k), dc], "rows": (k, hi), "readings": (e0, e1),
                "by_the_readings": abs(e1 - e0),
            })  # fmt: skip
    # Table X: the second difference on a line against the first differences on either side
    signed = {}
    for k in _rows_of(d, "X", "d1"):
        mag = _difference(d, "X", k, "d1")
        if mag is None:
            continue
        e0, e1 = _reading(d, "X", k), _reading(d, "X", k + 1)
        if e0 is not None and e1 is not None and e1 != e0:
            sign = 1 if e1 > e0 else -1
        else:
            sign = -1 if d["X", str(k), "d1"].startswith("-") else 0
        signed[k] = (sign, mag)
    for k in _rows_of(d, "X", "d2"):
        v = _difference(d, "X", k, "d2")
        if v is None or k not in signed or k - 1 not in signed:
            continue
        (s0, m0), (s1, m1) = signed[k - 1], signed[k]
        if 0 in (s0, s1):  # a sign the readings cannot give: compare magnitudes of the change
            s0 = s1 = s0 or s1 or 1
        by_first = abs(s1 * m1 - s0 * m0)
        if by_first != v:
            findings.append({
                "kind": "second difference",
                "table": "X", "row": k, "column": "d2",
                "new": ("X", str(k), "d2") in new,
                "reads": d["X", str(k), "d2"], "rows": (k - 1, k),
                "first_differences": (d["X", str(k - 1), "d1"], d["X", str(k), "d1"]),
                "by_the_first_differences": by_first,
            })  # fmt: skip
    # a first difference that disagrees only because a reading at one of its ends is
    # contradicted on both sides is reported with that reading, not again
    explained = {(f["table"], f["row"]) for f in findings if f["kind"] == "reading"}
    findings = [
        f
        for f in findings
        if f["kind"] != "first difference"
        or not {(f["table"], f["rows"][0]), (f["table"], f["rows"][1])} & explained
    ]
    return {"counts": counts, "findings": findings}


# Which epochs a logged cell reaches. Located from the entry plan's list of rows by epoch
# (rows found from the corrected formula's arguments, as the entry plan located them); no
# table longitude is formed for this. The primary path interpolates Tables IX, X and XII in
# four rows of readings and the other tables in two, and uses no printed difference; only
# S-two-row (§ 6) uses the differences, Table IX's and Table X's.
TOUCHES = {
    ("IX", 1790, "diff"): "not in the primary path; S-two-row at 1792, 1795, 1798 "
    "(the decade 1790 to 1800)",
    ("X", 84, "d1"): "not in the primary path; S-two-row at 1816 (anomaly 84.7)",
    ("X", 142, "d1"): "not in the primary path; S-two-row at 1828 (anomaly 142.6)",
    ("XIII", 40, "diff"): "none: no path uses Table XIII's printed differences",
    ("XIII", 720, "eq"): "the primary path at 1771 and 1813 (step 3d's reading, as used there)",
    ("X", 4, "d2"): "none: no epoch interpolates between anomalies 3 and 5",
    ("X", 83, "d2"): "none: no epoch interpolates between anomalies 82 and 84",
    ("X", 101, "d2"): "none: no epoch interpolates between anomalies 100 and 102",
    ("X", 141, "d2"): "none: no epoch interpolates between anomalies 140 and 142",
    ("X", 143, "d2"): "S-two-row at 1828 (anomaly 142.6)",
    ("X", 189, "d2"): "none: no epoch interpolates between anomalies 188 and 190",
}


def _printed(d: dict, table: str, row: int, col: str) -> str:
    lab = _label(d, table, row, col)
    return d[table, str(lab), col] if lab is not None else "not entered"


def _page(table: str, row: int, col: str) -> str:
    for path in (B.ROUTE_B_ROWS, B.BOUVARD / "bouvard_1821_uranus_table_rows.csv"):
        with path.open(encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if (r["table"], r["row"], r["column"]) == (table, str(row), col):
                    return r["page"]
    return ""


def misprint_rows(check: dict | None = None) -> list[dict[str, str]]:
    """The rows of ``bouvard_1821_misprints_route_b.csv``, in the columns of step 3d's log."""
    check = check or difference_check()
    d = {**B.table_rows(used=True), **B.route_b_rows(used=False)}
    out = []
    for f in check["findings"]:
        t, k, col = f["table"], f["row"], f["column"]
        step = DIFFERENCED[t][2]
        name = B.COLUMN_NAMES[col]
        fmt = (lambda v, r: _fmt_centre(v, r)) if t == "X" else (lambda v, r: _fmt(v, ""))
        if f["kind"] == "reading":
            lo, _, hi, _ = f["neighbours"]
            by = f["by_the_differences"]
            evidence = (
                f"printed differences above and below ({_printed(d, t, lo, DIFFERENCED[t][1])} "
                f"and {_printed(d, t, k, DIFFERENCED[t][1])}); with the readings of {lo} "
                f"({_printed(d, t, lo, 'eq')}) and {hi} ({_printed(d, t, hi, 'eq')}) "
                + ("each gives the same value" if f["same_value"] else "they give different values")
            )
            if f["new"]:
                status = (
                    f"{B.REPLACED}: both printed differences give the same value (§ 2)"
                    if f["same_value"]
                    else "contradicted on both sides; the two differences give different "
                    "values; left as read (§ 2)"
                )
                used = fmt(by[0], k) if f["same_value"] else f["reads"]
            else:
                status = (
                    "a reading of step 3d's rows, in its misprint log; its ruling stands "
                    "(left as read; § 9, choice 5); not replaced here"
                )
                used = f["reads"]
            out.append({
                "place": f"Table {t}, row {k}, {name}", "page": _page(t, k, col),
                "reads": f["reads"], "evidence": evidence,
                "value_by_the_evidence": " or ".join(fmt(v, k) for v in by) or "none",
                "status": status, "used": used, "touches": TOUCHES.get((t, k, col), ""),
            })  # fmt: skip
        elif f["kind"] == "first difference":
            a, b = f["rows"]
            ends = []
            for end, other, other_d in ((a, a - step, a - step), (b, b + step, b)):
                e_end, e_other = _reading(d, t, end), _reading(d, t, other)
                dv = _difference(d, t, other_d, DIFFERENCED[t][1])
                if e_other is None or dv is None:
                    ends.append(f"{end} cannot be checked on its other side")
                elif abs(e_end - e_other) == dv:
                    ends.append(
                        f"{end} agrees with its other printed difference "
                        f"({_printed(d, t, other_d, DIFFERENCED[t][1])})"
                    )
                else:
                    ends.append(f"{end} disagrees with its other printed difference too")
            evidence = (
                f"the readings of {a} ({_printed(d, t, a, 'eq')}) and {b} "
                f"({_printed(d, t, b, 'eq')}) differ by {_fmt(f['by_the_readings'], f['reads'])}; "
                + "; ".join(ends)
            )
            if t == "X":
                evidence += (
                    f"; the printed second differences on {a} ({_printed(d, t, a, 'd2')}) and "
                    f"{b} ({_printed(d, t, b, 'd2')})"
                )
            out.append({
                "place": f"Table {t}, between rows {a} and {b}, {name}",
                "page": _page(t, k, col), "reads": f["reads"], "evidence": evidence,
                "value_by_the_evidence": _fmt(f["by_the_readings"], f["reads"]),
                "status": "a printed difference the readings at both its ends contradict; no "
                "reading is contradicted on both sides, so nothing is replaced (§ 2)",
                "used": "the readings as read; the difference as printed",
                "touches": TOUCHES.get((t, k, col), ""),
            })  # fmt: skip
        else:
            p, q = f["first_differences"]
            out.append({
                "place": f"Table {t}, row {k}, {name}", "page": _page(t, k, col),
                "reads": f["reads"],
                "evidence": f"the printed first differences on either side ({p} above, {q} "
                f"below) differ by {_fmt(f['by_the_first_differences'], '')}",
                "value_by_the_evidence": _fmt(f["by_the_first_differences"], ""),
                "status": "a printed second difference the first differences contradict; not a "
                "reading of the rule of § 2; nothing is replaced",
                "used": "as printed",
                "touches": TOUCHES.get((t, k, col), ""),
            })  # fmt: skip
    return out


def write_misprint_log(path: Path | None = None) -> Path:
    """Step 2 of the order: the log, written before any longitude is formed."""
    path = path or B.ROUTE_B_MISPRINTS
    rows = misprint_rows()
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=B.MISPRINT_COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path


if __name__ == "__main__":
    check_inputs()
    if "--misprints" in sys.argv:
        print(write_misprint_log())
