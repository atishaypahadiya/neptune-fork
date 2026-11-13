"""Phase 1b, step 3d: the validation stages of docs/notes/P1B_bouvard_attribution.md.

V1 (Bouvard's printed formula against his tables) failed its pre-registered threshold. V1′
(the amendment of 13 November 2026) tests the corrected formula at seven further instants.
Nothing here runs a cell.

    uv run python -m fork.export.p1b_bouvard   # docs/notes/P1B_bouvard_attribution.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from fork.data.centesimal import CENTESIMAL_SEC_IN_ARCSEC as CS
from fork.historical import bouvard_1821 as B

ROOT = Path(__file__).resolve().parents[3]
OUT_JSON = ROOT / "docs" / "notes" / "P1B_bouvard_attribution.json"
FAIL_ARCSEC = 2.0  # § 5: over this at any instant, the formula is not the tables
CENTRE_FROM_TABLES = B.CENTRE_FROM_TABLES  # Table X's leading coefficient (printed 59427″.54)
V1_NAMES = ("1783", "1798", "1813", "1828", "1846")
RAD_IN_CS = 4_000_000.0 / (2.0 * math.pi)


def _wrap(x: float) -> float:
    return (x + B.CIRCLE_CS / 2) % B.CIRCLE_CS - B.CIRCLE_CS / 2


def eccentricity_of(lead_cs: float) -> float:
    """The eccentricity whose equation of the centre has this leading coefficient
    (2e − e³/4 + 5e⁵/96, in centesimal seconds)."""
    e = lead_cs / RAD_IN_CS / 2
    for _ in range(20):
        e = (lead_cs / RAD_IN_CS + e**3 / 4 - 5 * e**5 / 96) / 2
    return e


def v1_row(inst: B.Instant) -> dict:
    t = B.table_longitude(inst)
    as_read = B.table_longitude(inst, used=False)
    f = B.formula(inst.t)
    f2 = B.formula(inst.t, centre_coefficient=CENTRE_FROM_TABLES)
    d_mean = _wrap(f["mean_longitude"] - t["mean_longitude"])
    diff = _wrap(f["longitude"] - t["longitude"])
    diff2 = _wrap(f2["longitude"] - t["longitude"])
    return {
        "instant": inst.name,
        "date_paris_midnight": "{:04d}-{:02d}-{:02d}".format(*inst.date),
        "t_years_from_1800": inst.t,
        "anomaly_grades": t["anomaly_grades"],
        "entries": t["entries"],
        "budget_cs": t["budget_cs"],
        "budget_arcsec": t["budget_arcsec"],
        "argument_allowance_cs": t["argument_allowance_cs"],
        "formula_minus_tables_cs": diff,
        "formula_minus_tables_arcsec": diff * CS,
        "formula_minus_tables_as_read_cs": _wrap(f["longitude"] - as_read["longitude"]),
        "parts_cs": {
            "mean_longitude": d_mean,
            "precession_inequality": -t["precession_inequality"],
            "great_inequality": f["great_inequality"] - t["great_inequality"],
            "equation_of_centre": f["centre"] - t["centre"] - t["secular"],
            "perturbations": f["perturbations"] - t["perturbations"],
            "reduction": f["reduction"] - t["reduction"],
        },
        "mean_longitude_drift_cs_per_year_since_1780": d_mean / (inst.year - 1780.0),
        "with_centre_from_tables_cs": diff2,
        "with_centre_from_tables_less_mean_and_precession_cs": diff2
        - d_mean
        + t["precession_inequality"],
        "argument_formula_minus_table_parts": {
            a: (f[f"arg_{a}"] - t["args"][a] + B.ARG_PARTS[a] / 2) % B.ARG_PARTS[a]
            - B.ARG_PARTS[a] / 2
            for a in B.ARGS
        },
    }


def v1() -> dict:
    rows = {i.name: v1_row(i) for i in B.INSTANTS}
    five = [rows[n] for n in V1_NAMES]
    worst = max(abs(r["formula_minus_tables_arcsec"]) for r in five)
    over_budget = [r["instant"] for r in five if abs(r["formula_minus_tables_cs"]) > r["budget_cs"]]
    outcome = "fail" if worst > FAIL_ARCSEC else "stop and ask" if over_budget else "pass"
    yrs = np.array([float(n) for n in V1_NAMES]) - 1810.0
    design = np.vstack([np.ones(5), yrs]).T
    detrended = {}
    for key in ("formula_minus_tables_cs", "with_centre_from_tables_cs"):
        y = np.array([r[key] for r in five])
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        detrended[key] = {
            "constant_cs": float(coef[0]),
            "trend_cs_per_year": float(coef[1]),
            "residual_cs": dict(zip(V1_NAMES, (y - design @ coef).tolist(), strict=True)),
        }
    return {
        "instants": rows,
        "largest_difference_arcsec": worst,
        "over_budget": over_budget,
        "outcome": outcome,
        "constant_and_trend_removed": detrended,
        "sum_of_added_constants_cs": B.centre_constant(),
        "eccentricity": {
            "p_ii": 0.0466108,
            "of_printed_coefficient_59427.54": eccentricity_of(59427.54),
            "of_table_x_59417.54": eccentricity_of(CENTRE_FROM_TABLES),
            "adams": 0.046679,
        },
        "sign_before_group_iii": {
            "used": "-",
            "largest_difference_arcsec_with_plus": max(
                abs(
                    _wrap(
                        B.formula(i.t, sign_group_iii=1)["longitude"]
                        - B.table_longitude(i)["longitude"]
                    )
                )
                * CS
                for i in B.INSTANTS[:5]
            ),
        },
    }


def v1p_row(inst: B.Instant) -> dict:
    """One instant of V1′: the corrected formula against the tables, with the readings the
    printed differences contradict taken both ways."""
    t = B.table_longitude(inst)
    as_read = B.table_longitude(inst, used=False)
    f = B.corrected_formula(inst)
    diff = _wrap(f["longitude"] - t["longitude"])
    diff_read = _wrap(f["longitude"] - as_read["longitude"])
    limit = t["budget_cs"] + t["argument_allowance_cs"]
    ok, ok_read = abs(diff) <= limit, abs(diff_read) <= limit
    return {
        "instant": inst.name,
        "date_paris_midnight": "{:04d}-{:02d}-{:02d}".format(*inst.date),
        "t_years_from_1800": inst.t,
        "anomaly_grades": t["anomaly_grades"],
        "entries": t["entries"],
        "budget_arcsec": t["budget_arcsec"],
        "argument_allowance_arcsec": t["argument_allowance_cs"] * CS,
        "limit_arcsec": limit * CS,
        "corrected_minus_tables_arcsec": diff * CS,
        "corrected_minus_tables_as_read_arcsec": diff_read * CS,
        "within_budget_alone": abs(diff) <= t["budget_cs"],
        "outcome": "pass"
        if ok and ok_read
        else "fail"
        if not ok and not ok_read
        else "stop and ask",
        "printed_formula_minus_tables_arcsec": _wrap(
            B.formula(inst.t)["longitude"] - t["longitude"]
        )
        * CS,
        "parts_cs": {
            "mean_longitude": _wrap(f["mean_longitude"] - t["mean_longitude"]),
            "great_inequality": f["great_inequality"] - t["great_inequality"],
            "equation_of_centre": f["centre"] - t["centre"] - t["secular"],
            "perturbations": f["perturbations"] - t["perturbations"],
            "reduction": f["reduction"] - t["reduction"],
        },
        "argument_formula_minus_table_parts": {
            a: (f[f"arg_{a}"] - t["args"][a] + B.ARG_PARTS[a] / 2) % B.ARG_PARTS[a]
            - B.ARG_PARTS[a] / 2
            for a in B.ARGS
        },
    }


MODERN_V1P = ("1789", "1804", "1822", "1837")
ANCIENT_V1P = ("1690", "1750", "1769")


def v1p() -> dict:
    rows = {i.name: v1p_row(i) for i in B.INSTANTS_V1P}
    outcomes = {n: r["outcome"] for n, r in rows.items()}
    if any(o == "stop and ask" for o in outcomes.values()):
        consequence = "stop and ask"
    elif any(outcomes[n] == "fail" for n in MODERN_V1P):
        consequence = "route B for all 30 epochs"
    elif any(outcomes[n] == "fail" for n in ANCIENT_V1P):
        consequence = "route B for the ancient epochs; the corrected formula for the modern"
    else:
        consequence = "the corrected formula is C_Bouvard"
    return {
        "instants": rows,
        "outcomes": outcomes,
        "consequence": consequence,
        "not_part_of_the_rule": {i.name: v1p_row(i) for i in B.INSTANTS},
    }


def main() -> None:
    payload = {
        "plan": "docs/notes/P1B_bouvard_attribution.md",
        "V1": v1(),
        "V1_prime": v1p(),
        "V2": "the five terms of § 7 are transcribed (two entries, 15 cells, no disagreement); "
        "its checks are not run: the work stops after V1′ for a report",
        "V3": "not run",
        "cells": "not run",
    }
    OUT_JSON.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    r = payload["V1"]
    print("V1′:", payload["V1_prime"]["outcomes"], "->", payload["V1_prime"]["consequence"])
    print("V1:", r["outcome"], f"largest {r['largest_difference_arcsec']:.2f} arcsec")


if __name__ == "__main__":
    main()
