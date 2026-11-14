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

import numpy as np

from fork.historical import bouvard_1821 as B

CS = B.CENTESIMAL_SEC_IN_ARCSEC  # centesimal seconds to arcseconds

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


# ------------------------------------------------------------------ R1

FAIL_ARCSEC = 2.0  # § 3: over this, the entries of the epoch's rows are checked again
OUTCOMES = (
    "inside the rounding",
    "inside the rounding of the arguments",
    "outside: the formula does not stand for the tables at that epoch",
    "outside, over 2″: the entry of every row the epoch uses is to be checked again",
)


def r1_outcome(diff_cs: float, budget_cs: float, allowance_cs: float) -> str:
    """§ 3's rule (that of V1′), on the absolute difference, formula minus tables."""
    a = abs(diff_cs)
    if a <= budget_cs:
        return OUTCOMES[0]
    if a <= budget_cs + allowance_cs:
        return OUTCOMES[1]
    if a * CS <= FAIL_ARCSEC:
        return OUTCOMES[2]
    return OUTCOMES[3]


def _wrap(x: float) -> float:
    return (x + B.CIRCLE_CS / 2) % B.CIRCLE_CS - B.CIRCLE_CS / 2


def epochs() -> dict[int, float]:
    """The 30 epochs of step 3d as t (Julian years from the Paris midnight of 1800)."""
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A

    return {y: B.t_of_jd(jd) for y, jd in R.place_epochs(A.load_memoir()).items()}


def compared_before(t: float) -> str | None:
    """The instant of V1 or V1′ within a day of t, if any (§ 3: those epochs share most of
    their rows with it, and three corrections of the formula were found from five of them)."""
    for inst in (*B.INSTANTS, *B.INSTANTS_V1P):
        if abs(t - inst.t) * 365.25 < 1.0:
            return f"{inst.name} ({'V1' if inst in B.INSTANTS else 'V1′'})"
    return None


def r1_row(t: float, branch: str = "century") -> dict:
    tab = B.table_longitude_at(t, branch=branch)
    printed = B.table_longitude_at(t, used=False, branch=branch)
    f = B.corrected_formula_at(t, branch=branch)
    diff = _wrap(f["longitude"] - tab["longitude"])
    diff_printed = _wrap(f["longitude"] - printed["longitude"])
    out = r1_outcome(diff, tab["budget_cs"], tab["argument_allowance_cs"])
    out_printed = r1_outcome(diff_printed, printed["budget_cs"], printed["argument_allowance_cs"])
    return {
        "t_years_from_1800": t,
        "paris_date": "{:04d}-{:02d}-{:02d}".format(*tab["date"]),
        "decimal_time": "{h}h {m:02d}m {s:02d}s".format(**tab["time_of_day"]),
        "epoch_rows": [f"{'−' if s < 0 else ''}{tb} {r}" for tb, r, s in tab["epoch_rows"]],
        "anomaly_grades": tab["anomaly_grades"],
        "entries": tab["entries"],
        "budget_arcsec": tab["budget_arcsec"],
        "argument_allowance_arcsec": tab["argument_allowance_cs"] * CS,
        "limit_arcsec": (tab["budget_cs"] + tab["argument_allowance_cs"]) * CS,
        "tables_cs": tab["longitude"],
        "formula_cs": f["longitude"],
        "formula_minus_tables_arcsec": diff * CS,
        "outcome": out,
        "as_printed": {
            "formula_minus_tables_arcsec": diff_printed * CS,
            "outcome": out_printed,
            "differs": abs(diff_printed - diff) > 1e-9,
        },
        "parts_cs": {
            "mean_longitude": _wrap(f["mean_longitude"] - tab["mean_longitude"]),
            "great_inequality": f["great_inequality"] - tab["great_inequality"],
            "equation_of_centre": f["centre"] - tab["centre"] - tab["secular"],
            "perturbations": f["perturbations"] - tab["perturbations"],
            "reduction": f["reduction"] - tab["reduction"],
        },
        "argument_formula_minus_table_parts": {
            a: (f[f"arg_{a}"] - tab["args"][a] + B.ARG_PARTS[a] / 2) % B.ARG_PARTS[a]
            - B.ARG_PARTS[a] / 2
            for a in B.ARGS
        },
        "compared_before": compared_before(t),
        "rows_read": len(tab["rows_read"]),
        "reads_a_route_b_ruling": bool(set(tab["rows_read"]) & set(B.route_b_rulings())),
    }


def r1() -> dict:
    ep = epochs()
    rows = {str(y): r1_row(t) for y, t in ep.items() if y != 1780}
    beside = {
        "1780 (S8(a), the century branch)": r1_row(ep[1780], "century"),
        "1780 (S8(b), Table I's own branch)": r1_row(ep[1780], "own"),
    }
    fresh = [y for y, r in rows.items() if r["compared_before"] is None]
    seen = [y for y, r in rows.items() if r["compared_before"] is not None]

    def summary(names: list[str]) -> dict:
        diffs = {y: rows[y]["formula_minus_tables_arcsec"] for y in names}
        worst = max(diffs, key=lambda y: abs(diffs[y]))
        return {
            "n": len(names),
            "counts": {o: sum(rows[y]["outcome"] == o for y in names) for o in OUTCOMES},
            "largest": {"epoch": worst, "arcsec": diffs[worst]},
            "outside": [y for y in names if rows[y]["outcome"] in OUTCOMES[2:]],
        }

    return {
        "epochs": rows,
        "beside_1780": beside,
        "all_29": summary(list(rows)),
        "fresh_17": summary(fresh),
        "near_v1_or_v1p_12": summary(seen),
        # § 2: an epoch that reads a reading replaced under its rule, whose outcome differs
        # with the reading as printed (no reading of the new rows was replaced)
        "stop_and_ask": [
            y
            for y, r in rows.items()
            if r["reads_a_route_b_ruling"] and r["as_printed"]["outcome"] != r["outcome"]
        ],
        "over_2_arcsec": [y for y, r in rows.items() if r["outcome"] == OUTCOMES[3]],
    }


# ------------------------------------------------------------------ R2, the cells

GATE_TOL_DEG = 0.001
BOUNDARIES = (6.0, 3.0, -3.0)


def context() -> dict:
    """Step 3d's inputs to the cells (``fork.export.p1b_bouvard.cells``), unchanged: the
    memoir, the clean truth and reference, the conversions to mean longitude."""
    from fork.export import p1b_adams_repro as R
    from fork.export import p1b_bouvard as P
    from fork.historical import adams_1846 as A

    m = A.load_memoir()
    syn = R.synthetic_normal_places(m)
    years = syn["years"]
    jd = syn["jd"]
    clean_true = np.array([syn["oc_true_arcsec"][y] for y in years])
    clean_mean = np.array([syn["oc_mean_arcsec"][y] for y in years])
    v2r = P.v2()["s8_to_s9"]
    lon_c = np.degrees(syn["lon_true_rad"]) - clean_true / 3600.0
    return {
        "m": m,
        "ctx": {
            "truth": R.perturber_lon_1846(syn["truth_theta"], syn["ctx"]),
            "eps_c": syn["c_fit"]["eps_1810_ecl_j2000_deg"],
        },
        "years": years,
        "t": np.array([B.t_of_jd(x) for x in jd]),
        "clean_true": clean_true,
        "clean_mean": clean_mean,
        "factor_clean": clean_mean / clean_true,
        "factor_his": np.array([v2r[y]["factor"] for y in years]),
        "s8": np.array([p["oc_true"] for p in m.places]),
        "s9": {p["year"]: p["oc_mean"] for p in m.places},
        "lon_c_date": np.array(
            [R.lon_j2000_to_of_date(v % 360.0, x) for v, x in zip(lon_c, jd, strict=True)]
        ),
    }


def d_of(cx: dict, lon_b_deg) -> np.ndarray:
    """D = C − C_clean in arcseconds, of date, at the 30 places (step 3d's formula)."""
    return ((np.asarray(lon_b_deg) - cx["lon_c_date"] + 180.0) % 360.0 - 180.0) * 3600.0


def columns(cx: dict, d) -> tuple[dict, dict]:
    """Cells B and C's columns from D, as step 3d forms them."""
    years = cx["years"]
    col_b = dict(zip(years, ((cx["s8"] + d) * cx["factor_clean"]).tolist(), strict=True))
    col_c = dict(zip(years, ((cx["clean_true"] - d) * cx["factor_his"]).tolist(), strict=True))
    return col_b, col_c


def c_tables_deg(t: float, **kw) -> float:
    """C_tables: the tables' longitude of date (degrees of 360) plus the five terms of § 7."""
    return (
        B.table_longitude_at(t, **kw)["longitude"] * 0.9e-4 + B.s7_correction(t) / 3600.0
    ) % 360.0


def d_tables(cx: dict, place_1780: str = "interpolated", **kw):
    """D from the tables at the 30 places. The 1780 place: "interpolated" (the cubic through
    1769, 1771, 1783, 1786, § 1), "century" (S-a) or "own" (S-b), at its own date."""
    from fork.export import p1b_bouvard as P

    years = cx["years"]
    i80 = years.index(1780)
    lon = [
        c_tables_deg(t, branch=("own" if place_1780 == "own" else "century"), **kw)
        if y == 1780
        else c_tables_deg(t, **kw)
        for y, t in zip(years, cx["t"], strict=True)
    ]
    d = d_of(cx, lon)
    if place_1780 == "interpolated":
        d[i80] = P._d_1780_interpolated(years, cx["t"], d)
    return d


def solve_bc(cx: dict, d, log=None) -> dict:
    from fork.export import p1b_bouvard as P

    col_b, col_c = columns(cx, d)
    out = {}
    for hyp in (1, 2):
        out[f"H{hyp}"] = {
            "B": P.solve_cell(cx["m"], hyp, col_b, "clean", cx["ctx"]),
            "C": P.solve_cell(cx["m"], hyp, col_c, "bouvard", cx["ctx"]),
        }
        if log:
            log(f"  B and C solved, Hypothesis {'I' * hyp}")
    return out


def near_boundary(v: float) -> list[float]:
    return [b for b in BOUNDARIES if abs(v - b) < 1.0]


def reading(att: dict) -> dict:
    """§ 4's reading rule (step 3d's § 6) on the attributions of both methods and both
    hypotheses, as step 3d's runner states it, with the remarks of § 4."""
    h1, h2 = att["H1"], att["H2"]
    r = {k: h1["d1c"][k] for k in ("A", "B", "C", "D", "shift", "R", "Ob", "I", "R_class",
                                    "Ob_class", "row")}  # fmt: skip
    notes = []
    if h1["d1c"]["interaction_override"]:
        notes.append("interaction of 3° or more: no single attribution is stated")
    if (h1["chain"]["R_class"], h1["chain"]["Ob_class"]) != (
        h1["d1c"]["R_class"],
        h1["d1c"]["Ob_class"],
    ):
        notes.append("D1(c) and his chain fall in different rows: depends on the method")
    if (h2["d1c"]["R_class"], h2["d1c"]["Ob_class"]) != (
        h1["d1c"]["R_class"],
        h1["d1c"]["Ob_class"],
    ):
        notes.append("Hypothesis II falls in different classes")
    for k in ("R", "Ob"):
        for b in near_boundary(h1["d1c"][k]):
            notes.append(f"{k} ({h1['d1c'][k]:+.2f}°) lies within 1° of the boundary {b:+.0f}°")
    r["overrides_and_remarks"] = notes
    r["chain_row"] = h1["chain"]["row"]
    r["state"] = reading_state(h1["d1c"], h1["chain"])
    return r


def reading_state(d1c: dict, chain: dict) -> dict:
    """What § 5's three conditions compare, on Hypothesis I."""
    return {
        "d1c_classes": [d1c["R_class"], d1c["Ob_class"]],
        "d1c_row": d1c["row"],
        "chain_classes": [chain["R_class"], chain["Ob_class"]],
        "depends_on_the_method": (chain["R_class"], chain["Ob_class"])
        != (d1c["R_class"], d1c["Ob_class"]),
        "interaction_override": abs(d1c["I"]) >= 3.0,
    }


def same_reading(tables: dict, formula: dict) -> dict:
    """§ 5, fixed in advance: the reading with the tables is the same as with the formula
    when (1) D1(c)'s pair of classes falls in the same row; (2) "depends on the method"
    applies in both or neither (and where it does not, his chain falls in that same row in
    both); (3) the interaction override applies in both or neither. Hypothesis I."""
    c1 = tables["d1c_classes"] == formula["d1c_classes"]
    c2_override = tables["depends_on_the_method"] == formula["depends_on_the_method"]
    c2_chain = (
        tables["depends_on_the_method"]
        or formula["depends_on_the_method"]
        or tables["chain_classes"] == formula["chain_classes"]
    )
    c3 = tables["interaction_override"] == formula["interaction_override"]
    differs = []
    if not c1:
        differs.append("1: D1(c)'s row")
    if not c2_override:
        differs.append("2: depends on the method")
    elif not c2_chain:
        differs.append("2: his chain's row (the override applies in neither)")
    if not c3:
        differs.append("3: the interaction override")
    return {
        "same": not differs,
        "conditions": {
            "1: D1(c)'s row": c1,
            "2: depends on the method": c2_override,
            "2: his chain's row where the override applies in neither": c2_chain,
            "3: the interaction override": c3,
        },
        "differs_in": differs,
    }


def attribution(cells_: dict) -> dict:
    from fork.export import p1b_bouvard as P

    return {
        h: {meth: P.attribute(*(cells_[h][k][meth]["offset"] for k in "ABCD"))
            for meth in ("d1c", "chain")}
        for h in ("H1", "H2")
    }  # fmt: skip


def _formula_side() -> dict:
    """Step 3d's recorded numbers, the formula side of § 5's comparisons."""
    import json

    from fork.export import p1b_bouvard as P

    c = json.loads(P.OUT_JSON.read_text(encoding="utf-8"))["cells"]
    sides = {}
    for key, name in (("c", "primary: S8(c), 1780 interpolated"),
                      ("a", "beside: the registered run, 1780 on the century branch"),
                      ("b", "beside: S8(b), 1780 on Table I's own branch")):  # fmt: skip
        v = c["S8"]["variants"][key]
        att = {}
        for h in ("H1", "H2"):
            att[h] = {}
            for meth in ("d1c", "chain"):
                a = P.attribute(
                    c["cells"][h]["A"][meth]["offset"], v[h][meth]["B"], v[h][meth]["C"],
                    c["cells"][h]["D"][meth]["offset"],
                )  # fmt: skip
                att[h][meth] = a
        sides[key] = {"name": name, "attribution": att}
    return sides


def compare(tab_att: dict, form_att: dict) -> dict:
    t1, f1 = tab_att["H1"], form_att["H1"]
    st, sf = reading_state(t1["d1c"], t1["chain"]), reading_state(f1["d1c"], f1["chain"])
    beside = {}
    for h in ("H1", "H2"):
        for meth in ("d1c", "chain"):
            ta, fa = tab_att[h][meth], form_att[h][meth]
            beside[f"{h} {meth}"] = {
                "change_B": ta["B"] - fa["B"], "change_C": ta["C"] - fa["C"],
                "change_R": ta["R"] - fa["R"], "change_Ob": ta["Ob"] - fa["Ob"],
                "tables": {k: ta[k] for k in ("R", "Ob", "I", "R_class", "Ob_class")},
                "formula": {k: fa[k] for k in ("R", "Ob", "I", "R_class", "Ob_class")},
            }  # fmt: skip
    h2t, h2f = tab_att["H2"]["d1c"], form_att["H2"]["d1c"]
    return {
        "tables_state": st,
        "formula_state": sf,
        **same_reading(st, sf),
        "beside": beside,
        "hypothesis_II_changes_class": (h2t["R_class"], h2t["Ob_class"])
        != (h2f["R_class"], h2f["Ob_class"]),
    }


def run(log=print) -> dict:
    """The whole order of § 8 after the misprint log: R1, the gates, R2, R3, § 6."""
    import json

    from fork.export import p1b_bouvard as P

    check_inputs()
    out: dict = {"plan": "docs/notes/P1B_bouvard_tables.md", "inputs": "hashes as in § 1"}
    out["difference_check"] = {
        "counts": difference_check()["counts"],
        "logged": len(misprint_rows()),
        "replaced": B.route_b_rulings(),
    }
    log("R1")
    out["R1"] = r1()
    cx = context()
    # step 3d's columns from step 3d's C_Bouvard through this runner's code
    rec = json.loads(P.OUT_JSON.read_text(encoding="utf-8"))["cells"]["D"]
    d3d = d_of(cx, [B.c_bouvard_deg(x) for x in cx["t"]])
    cb, cc = columns(cx, d3d)
    out["step_3d_columns_reproduced_max_abs_arcsec"] = max(
        max(abs(cb[y] - rec["columns_arcsec"]["B"][str(y)]) for y in cx["years"]),
        max(abs(cc[y] - rec["columns_arcsec"]["C"][str(y)]) for y in cx["years"]),
    )
    log("gates: cells A and D")
    col_d = dict(zip(cx["years"], cx["clean_mean"].tolist(), strict=True))
    ad = {
        f"H{h}": {
            "A": P.solve_cell(cx["m"], h, cx["s9"], "bouvard", cx["ctx"]),
            "D": P.solve_cell(cx["m"], h, col_d, "clean", cx["ctx"]),
        }
        for h in (1, 2)
    }
    gates = {}
    for h, g in P.GATES.items():
        for name, must in g.items():
            cell, meth = name.split("_")
            got = ad[h][cell]["d1c" if meth == "d1c" else "chain"]["offset"]
            gates[f"{h} {name}"] = {
                "must_give": must, "gives": got,
                "pass": got is not None and abs(got - must) <= GATE_TOL_DEG,
            }  # fmt: skip
    out["gates"] = gates
    out["gates_pass"] = all(v["pass"] for v in gates.values())
    if not out["gates_pass"]:
        log("a gate failed: the work stops (§ 4)")
        return out
    log("R2: cells B and C, the 1780 place interpolated")
    d0 = d_tables(cx)
    out["D"] = {
        "years": cx["years"],
        "t_years_from_1800": cx["t"].tolist(),
        "D_tables_arcsec": d0.tolist(),
        "D_formula_arcsec_step_3d": d3d.tolist(),
        "tables_minus_formula_arcsec": (d0 - d3d).tolist(),
    }
    bc = solve_bc(cx, d0, log)
    cells_ = {h: {**ad[h], **bc[h]} for h in ("H1", "H2")}
    att = attribution(cells_)
    out["R2"] = {"cells": cells_, "attribution": att, "reading": reading(att)}
    form = _formula_side()
    variants = {"interpolated": (att, form["c"]["attribution"])}
    out["sensitivities"] = {}
    for name, kw in (("S-a", {"place_1780": "century"}), ("S-b", {"place_1780": "own"}),
                     ("S-printed", {"used": False}),
                     ("S-two-row", {"interpolation": "two-row"})):  # fmt: skip
        log(name)
        d = d_tables(cx, **kw)
        bc = solve_bc(cx, d, log)
        c = {h: {**ad[h], **bc[h]} for h in ("H1", "H2")}
        a = attribution(c)
        out["sensitivities"][name] = {
            "delta_D_arcsec": {
                str(y): float(v)
                for y, v in zip(cx["years"], d - d0, strict=True)
                if abs(v) > 0.005
            },
            "change": {
                h: {meth: {"B": c[h]["B"][meth]["offset"] - cells_[h]["B"][meth]["offset"],
                           "C": c[h]["C"][meth]["offset"] - cells_[h]["C"][meth]["offset"]}
                    for meth in ("d1c", "chain")}
                for h in ("H1", "H2")
            },
            "cells": {h: {k: {meth: c[h][k][meth]["offset"] for meth in ("d1c", "chain")}
                          for k in "BC"} for h in ("H1", "H2")},
            "attribution": a,
            "reading": reading(a),
        }  # fmt: skip
        if name == "S-a":
            variants["century"] = (a, form["a"]["attribution"])
        if name == "S-b":
            variants["own"] = (a, form["b"]["attribution"])
    out["sensitivities"]["S-midnight"] = s_midnight_status(cx)
    out["R3"] = {
        "primary (1780 interpolated; S8(c))": compare(*variants["interpolated"]),
        "beside (1780 on the century branch; S-a against the registered run)": compare(
            *variants["century"]
        ),
        "beside (1780 on Table I's own branch; S-b against S8(b))": compare(*variants["own"]),
    }
    return out


def s_midnight_status(cx: dict) -> dict:
    """S-midnight needs the rows of each epoch's nearest Paris midnight. Where one was not
    entered, § 2 stops at that epoch and nothing is filled in; the variant is not run."""
    missing = {}
    for y, t in zip(cx["years"], cx["t"], strict=True):
        if y == 1780:
            continue  # interpolated from its neighbours in the primary treatment
        try:
            B.table_longitude_at(float(t), midnight=True)
        except B.MissingRow as e:
            missing[str(y)] = {
                "midnight": "{:04d}-{:02d}-{:02d}".format(
                    *B.paris_date_and_time(float(t), midnight=True)["date"]
                ),
                "row_not_entered": list(e.args[0]),
            }
    return {
        "run": not missing,
        "status": "not run: rows not entered at these epochs' nearest midnights" if missing else "",
        "missing": missing,
    }


def main() -> None:
    import json
    import time

    from fork.export import p1b_bouvard as P

    t0 = time.time()
    progress = ROOT / "data" / "cache" / "p1b_bouvard_tables_progress.log"

    def log(msg: str) -> None:
        line = f"[{time.time() - t0:7.1f} s] {msg}"
        print(line, flush=True)
        try:
            with progress.open("a", encoding="utf-8") as f:
                f.write(line + "\n")
        except OSError:
            pass

    out = run(log)
    OUT_JSON.write_text(
        json.dumps(out, default=P._plain, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    log(f"written {OUT_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    check_inputs()
    if "--misprints" in sys.argv:
        print(write_misprint_log())
    else:
        main()
