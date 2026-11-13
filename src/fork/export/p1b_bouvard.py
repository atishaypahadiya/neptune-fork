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


def _plain(v):
    return bool(v) if isinstance(v, np.bool_) else float(v)


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
        "ruled_2026_11_13": {
            "table_xii_row_2600": "3.31,6",
            "outcomes": dict.fromkeys(outcomes, "pass")
            if all(
                abs(r["corrected_minus_tables_arcsec"]) <= r["limit_arcsec"] for r in rows.values()
            )
            else None,
            "consequence": "the corrected formula is C_Bouvard for all 30 epochs",
        },
        "not_part_of_the_rule": {i.name: v1p_row(i) for i in B.INSTANTS},
    }


# ------------------------------------------------------------------ V2

TOL_S9_MODERN = 0.02  # arcseconds, his printing precision
TOL_S9_ANCIENT = 0.1


def _detrended_max(t: np.ndarray, y: np.ndarray) -> float:
    design = np.vstack([np.ones(len(t)), t]).T
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return float(np.abs(y - design @ coef).max())


def v2() -> dict:
    """The § 7 terms at his 30 epochs, and his § 8 column times r²/(ab) against his § 9."""
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A

    m = A.load_memoir()
    e = m.const(1, 22, "e")
    beta = A.setup_from_print(m, 1).beta
    n = m.angle(1, 13, "n")
    conv = {}
    for p in m.places:
        mean_anom = math.radians(beta + n * float(p["t_h1"]))
        ecc_anom = mean_anom
        for _ in range(60):
            ecc_anom = mean_anom + e * math.sin(ecc_anom)
        factor = (1 - e * math.cos(ecc_anom)) ** 2 / math.sqrt(1 - e * e)
        conv[p["year"]] = {
            "kind": p["kind"],
            "factor": factor,
            "s8": p["oc_true"],
            "s9": p["oc_mean"],
            "s8_times_factor_minus_s9": p["oc_true"] * factor - p["oc_mean"],
        }
    worst = {
        k: max(abs(v["s8_times_factor_minus_s9"]) for v in conv.values() if v["kind"] == k)
        for k in ("modern", "ancient")
    }
    epochs = R.place_epochs(m)
    t = np.array([B.t_of_jd(jd) for jd in epochs.values()])
    s7 = np.array([B.s7_correction(x) for x in t])
    s7_great = np.array([B.s7_correction(x, great_inequalities=True) for x in t])
    lon_date = []
    for x in t:  # the same terms with the mean longitudes counted from the equinox of date
        lon = [v + 154.63 * 0.324 * x / 3600.0 for v in B.mean_longitudes_deg(x)]
        lon_date.append(
            sum(
                c * math.sin(math.radians(sum(k * v for k, v in zip(mu, lon, strict=True)) + k0))
                for c, mu, k0 in B.s7_terms()
            )
        )
    lon_date = np.array(lon_date)
    return {
        "terms": [
            {
                "coefficient_arcsec": c,
                "multipliers_jupiter_saturn_uranus": list(mu),
                "constant_deg": k,
            }
            for c, mu, k in B.s7_terms()
        ],
        "s7_at_the_30_epochs_arcsec": dict(zip(epochs, s7.tolist(), strict=True)),
        "s7_range_arcsec": [float(s7.min()), float(s7.max())],
        "s7_after_constant_and_trend_max_abs_arcsec": _detrended_max(t, s7),
        "each_term_after_constant_and_trend_max_abs_arcsec": [
            _detrended_max(t, np.array([B.s7_correction(x, terms=(k,)) for x in t]))
            for k in range(5)
        ],
        "ambiguity_great_inequalities_after_constant_and_trend_max_abs_arcsec": _detrended_max(
            t, s7_great - s7
        ),
        "ambiguity_equinox_of_date_after_constant_and_trend_max_abs_arcsec": _detrended_max(
            t, lon_date - s7
        ),
        "s8_to_s9": conv,
        "s8_to_s9_worst_arcsec": worst,
        "s8_to_s9_pass": worst["modern"] <= TOL_S9_MODERN and worst["ancient"] <= TOL_S9_ANCIENT,
        "s5_shift_arcsec": {
            y: B.misprint_row_2600_shift(x)
            for y, x in zip(epochs, t, strict=True)
            if abs(B.misprint_row_2600_shift(x)) > 0.005
        },
    }


# ------------------------------------------------------------------ V3

GATES = {  # § 0 and § 5 of the plan, degrees
    "H1": {"A_d1c": 1.647, "A_chain": 1.021, "D_d1c": -9.512, "D_chain": -9.923},
    "H2": {"A_d1c": 4.964, "A_chain": 4.447, "D_d1c": -7.621, "D_chain": -8.127},
}
TOL_GATE_DEG = 0.001
TOL_FRAME_OFFSET_DEG = 0.05
TOL_FRAME_CURVE_ARCSEC = 0.05
PROPAGATION_D_ARCSEC = 1.05  # the amendment of 13 November 2026, ruling 5


def _solve_column(m, hyp: int, col: dict, eps_c: float, lon_truth: float) -> dict:
    """A cell's column by both methods: the offset of the solution's true longitude at
    1846.762 from the truth, in the data's frame (as steps 3b and 3c)."""
    import copy

    from fork.export import p1b_adams_noise as N
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A

    s = A.setup_from_print(m, hyp)
    rows = A.place_rows(m, s, oc=col)
    w = A.weights_for(m, hyp, "c", rows)
    ls = A.solve_least_squares(rows, w, s, A.printed_solution(m, hyp))
    x = np.array([ls["x"][k] for k in A.PARAMS])
    lon = N.longitudes_of(x, eps_c, s.n_prime_deg)
    ms = copy.deepcopy(m)
    for p in ms.places:
        p["oc_mean"] = col[p["year"]]
    ch = N.chain_solve(ms, hyp, rows, w, s)
    out = {"d1c": R._wrap(lon["true_lon_1846_ecl_j2000"] - lon_truth), "chain": None}
    if ch["converged"]:
        lon_ch = N.longitudes_of(ch["x"], eps_c, s.n_prime_deg)
        out["chain"] = R._wrap(lon_ch["true_lon_1846_ecl_j2000"] - lon_truth)
    return out


def v3() -> dict:
    """The regression gates (cells A and D by both methods, both hypotheses) and the two
    frame tests. Cells B and C are not touched."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export import p1b_adams_noise as N
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A
    from fork.model.forward import predict_states_ecl_j2000
    from fork.model.observation import I_URANUS

    m = A.load_memoir()
    syn = R.synthetic_normal_places(m)
    truth = R.perturber_lon_1846(syn["truth_theta"], syn["ctx"])
    col = {int(k): float(v) for k, v in syn["oc_mean_arcsec"].items()}
    eps_c = syn["c_fit"]["eps_1810_ecl_j2000_deg"]
    got, gates = {}, {}
    for hyp in (1, 2):
        real = N._real_data_offsets(m, hyp, truth["of_date"])["O"]
        clean = _solve_column(m, hyp, col, eps_c, truth["ecl_j2000"])
        got[f"H{hyp}"] = {
            "A_d1c": real["d1c"],
            "A_chain": real["C2"],
            "D_d1c": clean["d1c"],
            "D_chain": clean["chain"],
        }
        gates[f"H{hyp}"] = {
            k: {
                "must_give": GATES[f"H{hyp}"][k],
                "gives": v,
                "pass": v is not None and abs(v - GATES[f"H{hyp}"][k]) <= TOL_GATE_DEG,
            }
            for k, v in got[f"H{hyp}"].items()
        }
    # frame test 1: a constant and a trend added to cell D's column
    shifted = {
        y: v + 10.0 + 0.2 * (A.EPOCH_YEAR + A.UNIT_YR * float(p["t_h1"]) - A.EPOCH_YEAR)
        for (y, v), p in zip(col.items(), m.places, strict=True)
    }
    moved = _solve_column(m, 1, shifted, eps_c, truth["ecl_j2000"])
    frame1 = {
        k: {
            "cell_D": got["H1"][f"D_{k}"],
            "with_10_arcsec_and_0.2_arcsec_a_year": moved[k],
            "change_deg": moved[k] - got["H1"][f"D_{k}"],
            "pass": abs(moved[k] - got["H1"][f"D_{k}"]) <= TOL_FRAME_OFFSET_DEG,
        }
        for k in ("d1c", "chain")
    }
    # frame test 2: the ephemeris's Uranus minus C_clean, in ecl_j2000 and of date
    jd = syn["jd"]
    pos_c, _ = predict_states_ecl_j2000(syn["c_fit"]["theta"], jd, syn["ctx"])
    pos_c = pos_c[:, I_URANUS]
    pos_e = about("uranus", jd)[:, :3]

    def lon_j(p):
        return np.degrees(np.arctan2(p[:, 1], p[:, 0]))

    def lon_d(p):
        return np.array([R.lon_of_date_deg(q, t) for q, t in zip(p, jd, strict=True)])

    def wrap(d):
        return ((d + 180.0) % 360.0 - 180.0) * 3600.0

    years = (jd - R.JD_1810) / 365.25
    design = np.vstack([np.ones(len(jd)), years]).T
    curves = {}
    for name, f in (("ecl_j2000", lon_j), ("of_date", lon_d)):
        y = wrap(f(pos_e) - f(pos_c))
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        curves[name] = y - design @ coef
    diff = float(np.abs(curves["ecl_j2000"] - curves["of_date"]).max())
    return {
        "gates": gates,
        "gates_pass": all(v["pass"] for g in gates.values() for v in g.values()),
        "frame_test_offset": frame1,
        "frame_test_curve": {
            "curve_range_arcsec": [
                float(curves["ecl_j2000"].min()),
                float(curves["ecl_j2000"].max()),
            ],
            "max_abs_difference_arcsec": diff,
            "pass": diff <= TOL_FRAME_CURVE_ARCSEC,
        },
        "propagation_check": {
            "d_arcsec": PROPAGATION_D_ARCSEC,
            "errors": [
                "d sin(M), M Uranus's mean anomaly",
                "d cos(M)",
                "d sin(III), III the argument φ′ − 2φ″",
                "d cos(III)",
            ],
            "status": "needs cells B and C; computed with the cells",
        },
    }


def main() -> None:
    payload = {
        "plan": "docs/notes/P1B_bouvard_attribution.md",
        "V1": v1(),
        "V1_prime": v1p(),
        "V2": v2(),
        "V3": v3(),
        "cells": "not run",
    }
    OUT_JSON.write_text(
        json.dumps(payload, default=_plain, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("V2: § 8 to § 9", payload["V2"]["s8_to_s9_worst_arcsec"])
    print("V3: gates", payload["V3"]["gates_pass"], payload["V3"]["frame_test_offset"])
    print("V3: frame curve", payload["V3"]["frame_test_curve"])
    r = payload["V1"]
    print("V1′:", payload["V1_prime"]["outcomes"], "->", payload["V1_prime"]["consequence"])
    print("V1:", r["outcome"], f"largest {r['largest_difference_arcsec']:.2f} arcsec")


if __name__ == "__main__":
    main()
