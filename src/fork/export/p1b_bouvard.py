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


# ------------------------------------------------------------------ the cells

SHIFT_CLASSES = ((6.0, "carries"), (3.0, "partial"), (-3.0, "does not carry"))
INTERACTION_OVERRIDE_DEG = 3.0
E_P_II = 0.0466108


def effect_class(v: float) -> str:
    """§ 6: ≥ +6 carries; over +3 partial; −3 to +3 does not carry; under −3 works against."""
    if v >= 6.0:
        return "carries"
    if v > 3.0:
        return "partial"
    if v >= -3.0:
        return "does not carry"
    return "works against"


READING_ROWS = {
    ("carries", "does not carry"): "The reference carries the shift.",
    ("does not carry", "carries"): "The observations carry the shift.",
    ("carries", "partial"): "Mostly the reference, with a share from the observations.",
    ("partial", "carries"): "Mostly the observations, with a share from the reference.",
    ("partial", "partial"): "Both, in the shares found.",
    ("carries", "works against"): "The reference carries more than the whole shift and the "
    "observations take some back.",
    ("works against", "carries"): "The observations carry more than the whole shift and the "
    "reference takes some back.",
}


def attribute(a: float, b: float, c: float, d: float) -> dict:
    """§ 4 and § 6 on four offsets."""
    r = 0.5 * ((a - b) + (c - d))
    ob = 0.5 * ((a - c) + (b - d))
    inter = 0.5 * ((a - b) - (c - d))
    cls = (effect_class(r), effect_class(ob))
    out = {
        "A": a, "B": b, "C": c, "D": d,
        "shift": a - d,
        "R": r,
        "Ob": ob,
        "I": inter,
        "R_class": cls[0],
        "Ob_class": cls[1],
        "row": READING_ROWS.get(cls, "no row of § 6 fits these classes"),
        "reference_effect_at_each_level": {"real observations": a - b, "clean observations": c - d},
        "observations_effect_at_each_level": {"Bouvard": a - c, "clean reference": b - d},
    }  # fmt: skip
    out["interaction_override"] = abs(inter) >= INTERACTION_OVERRIDE_DEG
    return out


def _of_date_lon(x: np.ndarray, s) -> float:
    """A solution's true longitude at 1846.762 with his ε and precession (cells A and C)."""
    from fork.historical import adams_1846 as A

    eps_prime = (s.eps - float(x[4])) % 360.0
    mean = (eps_prime + 12.0 * s.n_prime_deg + s.prec) % 360.0
    varpi = (s.eps - float(x[7]) + s.prec) % 360.0
    return A.true_longitude(mean, float(x[6]) / 20.0, varpi)


def solve_cell(m, hyp: int, col: dict, frame: str, ctx: dict) -> dict:
    """One column by both methods. ``frame`` is "bouvard" (his ε, of date) or "clean"
    (the C_clean orbit's ε, ecl_j2000), the frame of the cell's C."""
    import copy

    from fork.export import p1b_adams_noise as N
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A

    s = A.setup_from_print(m, hyp)
    rows = A.place_rows(m, s, oc=col)
    w = A.weights_for(m, hyp, "c", rows)
    ls = A.solve_least_squares(rows, w, s, A.printed_solution(m, hyp))
    x = np.array([ls["x"][k] for k in A.PARAMS])
    ms = copy.deepcopy(m)
    for p in ms.places:
        p["oc_mean"] = col[p["year"]]
    ch = N.chain_solve(ms, hyp, rows, w, s)

    def offset(v):
        if frame == "clean":
            lon = N.longitudes_of(v, ctx["eps_c"], s.n_prime_deg)["true_lon_1846_ecl_j2000"]
            return R._wrap(lon - ctx["truth"]["ecl_j2000"])
        return R._wrap(_of_date_lon(v, s) - ctx["truth"]["of_date"])

    return {
        "d1c": {
            "offset": offset(x),
            "theta": float(x[4]),
            "m_prime": float(x[5]),
            "e_prime_units": float(x[6]),
            "chi2": ls["chi2"],
            "rms_arcsec": ls["rms_arcsec"],
            "n_distinct_minima": ls["n_distinct_minima"],
        },
        "chain": {
            "offset": offset(ch["x"]) if ch["converged"] else None,
            "converged": ch["converged"],
            "n_positive_mass_roots": ch["n_positive"],
            "theta": float(ch["x"][4]) if ch["converged"] else None,
            "m_prime": float(ch["x"][5]) if ch["converged"] else None,
            "chi2_c_weights": ch.get("chi2"),
        },
    }


S8_VARIANTS = {
    "a": None,  # the century branch: the cells as run
    "b": "S8b the 1780 place on Table I's branch",
    "c": "S8c the 1780 place: D interpolated from its neighbours",
}


def s8_summary(c: dict) -> dict:
    """S8 (defined after the cells were seen): the attribution each variant of the 1780
    place would give, whether the two methods fall in one row of § 6, and the ranges of R
    and Ob over variants (b) and (c). Arithmetic on the cells and sensitivities only."""
    out: dict = {"variants": {}}
    r_all, ob_all = [], []
    for v, name in S8_VARIANTS.items():
        entry: dict = {}
        for hyp in ("H1", "H2"):
            base = c["cells"][hyp]
            rows = {}
            for meth in ("d1c", "chain"):
                a, d = base["A"][meth]["offset"], base["D"][meth]["offset"]
                if name is None:
                    b, cc = base["B"][meth]["offset"], base["C"][meth]["offset"]
                else:
                    s = c["sensitivities"][name][hyp][meth]
                    b, cc = s["B"], s["C"]
                att = attribute(a, b, cc, d)
                rows[meth] = {k: att[k] for k in ("B", "C", "R", "Ob", "I", "R_class",
                                                  "Ob_class", "row")}  # fmt: skip
                if name is not None:
                    r_all.append(att["R"])
                    ob_all.append(att["Ob"])
            rows["methods_in_one_row"] = rows["d1c"]["row"] == rows["chain"]["row"]
            rows["methods_differ_in_R_deg"] = abs(rows["d1c"]["R"] - rows["chain"]["R"])
            entry[hyp] = rows
        if name is not None:
            delta = c["sensitivities"][name]["delta_arcsec"]
            entry["delta_1780_arcsec"] = delta[1780] if 1780 in delta else delta["1780"]
        out["variants"][v] = entry
    out["b_and_c"] = {
        "R_min": min(r_all), "R_max": max(r_all),
        "Ob_min": min(ob_all), "Ob_max": max(ob_all),
        "n_rows": len(r_all),
    }  # fmt: skip
    out["label"] = "post hoc: S8 was defined after the cells were seen"
    return out


def _d_1780_interpolated(years: list, t: np.ndarray, d0: np.ndarray) -> float:
    """S8 (c): the cubic through D at 1769, 1771, 1783 and 1786, at the 1780 place's epoch."""
    idx = [years.index(y) for y in (1769, 1771, 1783, 1786)]
    x = t[years.index(1780)]
    total = 0.0
    for i in idx:
        w = 1.0
        for j in idx:
            if j != i:
                w *= (x - t[j]) / (t[i] - t[j])
        total += w * d0[i]
    return float(total)


def cells(sensitivities: bool = True) -> dict:
    """The four cells for both methods and both hypotheses, the attribution, the reading
    and the sensitivities (the plan's §§ 2 to 7 and its amendments)."""
    from fork.export import p1b_adams_repro as R
    from fork.historical import adams_1846 as A

    m = A.load_memoir()
    syn = R.synthetic_normal_places(m)
    ctx = {
        "truth": R.perturber_lon_1846(syn["truth_theta"], syn["ctx"]),
        "eps_c": syn["c_fit"]["eps_1810_ecl_j2000_deg"],
    }
    years = syn["years"]
    jd = syn["jd"]
    t = np.array([B.t_of_jd(x) for x in jd])
    clean_true = np.array([syn["oc_true_arcsec"][y] for y in years])
    clean_mean = np.array([syn["oc_mean_arcsec"][y] for y in years])
    factor_clean = clean_mean / clean_true  # r²/(ab) of the C_clean ellipse
    v2r = v2()["s8_to_s9"]
    factor_his = np.array([v2r[y]["factor"] for y in years])
    s8 = np.array([p["oc_true"] for p in m.places])
    s9 = {p["year"]: p["oc_mean"] for p in m.places}
    lon_c = np.degrees(syn["lon_true_rad"]) - clean_true / 3600.0
    lon_c_date = np.array(
        [R.lon_j2000_to_of_date(v % 360.0, x) for v, x in zip(lon_c, jd, strict=True)]
    )
    lon_b = np.array([B.c_bouvard_deg(x) for x in t])
    d0 = ((lon_b - lon_c_date + 180.0) % 360.0 - 180.0) * 3600.0  # D, arcseconds
    use = np.array([y != 1690 for y in years])
    design = np.vstack([np.ones(use.sum()), t[use]]).T
    coef, *_ = np.linalg.lstsq(design, d0[use], rcond=None)
    d_detr = d0 - (coef[0] + coef[1] * t)

    def columns(delta: np.ndarray) -> tuple[dict, dict]:
        d = d0 + delta
        col_b = dict(zip(years, ((s8 + d) * factor_clean).tolist(), strict=True))
        col_c = dict(zip(years, ((clean_true - d) * factor_his).tolist(), strict=True))
        return col_b, col_c

    zero = np.zeros(len(years))
    col_b, col_c = columns(zero)
    col_d = dict(zip(years, clean_mean.tolist(), strict=True))
    out: dict = {
        "D": {
            "years": years,
            "t_years_from_1800": t.tolist(),
            "c_bouvard_of_date_deg": lon_b.tolist(),
            "c_clean_of_date_deg": lon_c_date.tolist(),
            "D_arcsec": d0.tolist(),
            "constant_arcsec": float(coef[0]),
            "trend_arcsec_per_year": float(coef[1]),
            "D_less_constant_and_trend_arcsec": d_detr.tolist(),
            "columns_arcsec": {"A": s9, "B": col_b, "C": col_c, "D": col_d},
        },
        "cells": {},
        "attribution": {},
    }
    base = {}
    for hyp in (1, 2):
        c = {
            "A": solve_cell(m, hyp, s9, "bouvard", ctx),
            "B": solve_cell(m, hyp, col_b, "clean", ctx),
            "C": solve_cell(m, hyp, col_c, "bouvard", ctx),
            "D": solve_cell(m, hyp, col_d, "clean", ctx),
        }
        out["cells"][f"H{hyp}"] = c
        base[hyp] = c
        out["attribution"][f"H{hyp}"] = {
            meth: attribute(*(c[k][meth]["offset"] for k in "ABCD"))
            if all(c[k][meth]["offset"] is not None for k in "ABCD")
            else None
            for meth in ("d1c", "chain")
        }
    h1 = out["attribution"]["H1"]
    reading = dict(h1["d1c"])
    notes = []
    if h1["d1c"]["interaction_override"]:
        notes.append("interaction of 3° or more: no single attribution is stated")
    if h1["chain"] is None:
        notes.append("his chain gave no solution in a cell")
    elif (h1["chain"]["R_class"], h1["chain"]["Ob_class"]) != (
        h1["d1c"]["R_class"],
        h1["d1c"]["Ob_class"],
    ):
        notes.append("D1(c) and his chain fall in different rows: depends on the method")
    h2 = out["attribution"]["H2"]["d1c"]
    if h2 and (h2["R_class"], h2["Ob_class"]) != (h1["d1c"]["R_class"], h1["d1c"]["Ob_class"]):
        notes.append("Hypothesis II falls in different classes")
    reading["overrides_and_remarks"] = notes
    out["reading"] = reading

    if not sensitivities:
        return out

    # sensitivities: a change δ of C_Bouvard (arcseconds, true longitude) at the 30 epochs
    form = [B.corrected_formula_at(x) for x in t]
    anomaly = np.array([f["anomaly_grades"] for f in form]) * math.pi / 200.0
    arg3 = np.array([f["arg_III"] for f in form]) * 2 * math.pi / B.ARG_PARTS["III"]
    e_tables = eccentricity_of(CENTRE_FROM_TABLES)
    d = PROPAGATION_D_ARCSEC
    deltas = {
        "S2 eccentricity of p. ii": np.array(
            [B.eccentricity_shift(x, e_tables, E_P_II) for x in t]
        ),
        "S3 without § 7": -np.array([B.s7_correction(x) for x in t]),
        "S5 Table XII row 2600 as printed": np.array([B.misprint_row_2600_shift(x) for x in t]),
        "S5b Table X row 25 as printed": np.array([B.misprint_row_25_shift(x) for x in t]),
        "S6 § 7 with the great inequalities": np.array(
            [B.s7_correction(x, great_inequalities=True) - B.s7_correction(x) for x in t]
        ),
        "S7 § 7 from the equinox of date": np.array(
            [B.s7_correction_of_date(x) - B.s7_correction(x) for x in t]
        ),
        "propagation d sin M": d * np.sin(anomaly),
        "propagation d cos M": d * np.cos(anomaly),
        "propagation d sin III": d * np.sin(arg3),
        "propagation d cos III": d * np.cos(arg3),
        # post hoc, not in the plan: a frame rate in true longitude
        "post hoc: 0.1 arcsec a year in true longitude": 0.1 * (t - t[use].mean()),
        # S8 (defined after the cells were seen): the 1780 place falls 0.04 year before
        # 1780, so the first amendment's rule gives it the century row (variant a, as run);
        # (b) puts it on Table I's own branch; (c) interpolates D from its neighbours
        "S8c the 1780 place: D interpolated from its neighbours": np.array(
            [
                _d_1780_interpolated(years, t, d0) - d0[years.index(1780)] if y == 1780 else 0.0
                for y in years
            ]
        ),
        "S8b the 1780 place on Table I's branch": np.array(
            [B.MEAN_SLIP_PER_YEAR * 100.0 * CS if y == 1780 else 0.0 for y in years]
        ),
    }
    sens = {}
    for name, delta in deltas.items():
        cb, cc = columns(delta)
        entry = {
            "delta_arcsec": {
                y: float(v) for y, v in zip(years, delta, strict=True) if abs(v) > 0.005
            }
            if name.startswith(("S5", "S8"))
            else {"min": float(delta.min()), "max": float(delta.max())},
        }
        for hyp in (1, 2):
            sb = solve_cell(m, hyp, cb, "clean", ctx)
            sc = solve_cell(m, hyp, cc, "bouvard", ctx)
            entry[f"H{hyp}"] = {
                meth: {
                    "B": sb[meth]["offset"],
                    "C": sc[meth]["offset"],
                    "change_B": None
                    if sb[meth]["offset"] is None or base[hyp]["B"][meth]["offset"] is None
                    else sb[meth]["offset"] - base[hyp]["B"][meth]["offset"],
                    "change_C": None
                    if sc[meth]["offset"] is None or base[hyp]["C"][meth]["offset"] is None
                    else sc[meth]["offset"] - base[hyp]["C"][meth]["offset"],
                }
                for meth in ("d1c", "chain")
            }
        sens[name] = entry
    out["sensitivities"] = sens
    out["S8"] = s8_summary(out)
    out["S1"] = "not run: no primary source for Bouvard's mass of Saturn"
    return out


def figure(c: dict) -> None:
    """D at the 29 places with its constant and trend removed."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = c["D"]
    yr = 1800.0 + np.array(d["t_years_from_1800"])
    keep = np.array([y != 1690 for y in d["years"]])
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    ax.axhline(0.0, color="0.75", lw=0.8)
    ax.plot(yr[keep], np.array(d["D_less_constant_and_trend_arcsec"])[keep], "o-",
            color="#1f4e79", ms=4, lw=1.2)  # fmt: skip
    ax.set_xlabel("year")
    ax.set_ylabel("arcseconds")
    ax.set_title("Bouvard's reference minus the clean one, constant and trend removed")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_JSON.with_suffix(".png"), dpi=160)
    plt.close(fig)


def main() -> None:
    payload = {
        "plan": "docs/notes/P1B_bouvard_attribution.md",
        "V1": v1(),
        "V1_prime": v1p(),
        "V2": v2(),
        "V3": v3(),
        "cells": cells(),
    }
    OUT_JSON.write_text(
        json.dumps(payload, default=_plain, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    figure(payload["cells"])
    print(
        "reading:",
        {k: payload["cells"]["reading"][k] for k in ("A", "B", "C", "D", "R", "Ob", "I", "row")},
    )
    print("V2: § 8 to § 9", payload["V2"]["s8_to_s9_worst_arcsec"])
    print("V3: gates", payload["V3"]["gates_pass"], payload["V3"]["frame_test_offset"])
    print("V3: frame curve", payload["V3"]["frame_test_curve"])
    r = payload["V1"]
    print("V1′:", payload["V1_prime"]["outcomes"], "->", payload["V1_prime"]["consequence"])
    print("V1:", r["outcome"], f"largest {r['largest_difference_arcsec']:.2f} arcsec")


if __name__ == "__main__":
    main()
