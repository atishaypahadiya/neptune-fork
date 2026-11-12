"""Phase 1b step 3b: the Adams reproduction (docs/notes/P1B_adams_reproduction.md, pre-registered).

Stages B and C live in ``fork.historical.adams_1846``; this runner adds stage D and the report:
D1 his equations and model as one weighted least-squares problem (three weightings), D2 the α
profile 0.40–0.60 with the equations rebuilt at each α, D3 one noise-free truth fitted three ways
((i) his layout and method, (ii) his layout and the project's model, (iii) set D's layout through
step 3's H1(G) nested job), D4 the placements in step 3's fixed-a posteriors at 1846.762, then
the readings of § 5 in order and the figure. Results: docs/notes/P1B_adams_reproduction.{json,png};
intermediate stages are cached under data/cache/p1b_adams_repro/ (gitignored).
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.export import m6c
from fork.export import p1b_adams as pa
from fork.export.to_json import A1_EPOCH_JD_TDB
from fork.historical import adams_1846 as A

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "p1b_adams_repro"
PROGRESS_LOG = CACHE / "progress.log"
DRIVER_LOG = CACHE / "driver.log"
CHEAP_JSON = CACHE / "cheap.json"
D3III_JSON = CACHE / "d3iii.json"
NOTES = ROOT / "docs" / "notes"
OUT_JSON = NOTES / "P1B_adams_reproduction.json"
OUT_PNG = NOTES / "P1B_adams_reproduction.png"
MEMOIR_DIR = ROOT / "data" / "historical" / "adams_challis"
MEMOIR_FILES = ("normal_places", "series", "equations", "solution", "statements")

HELD_A_AU = 38.4  # D3: every fit at Hypothesis I's distance (α = 0.5)
SEED_NAMESPACE = (2026, 11, 3, 7)  # D3(iii): [2026, 11, 3, 7, m, k + 1], m = ERRORS.index("g")
SUMMARY_RNG = pa.SUMMARY_RNG  # D4 reuses step 3's draws; D3(iii) takes index 7 in it
D3III_RNG_INDEX = 7
NLIVE = pa.NLIVE_DEFAULT
WORKERS_PER_CELL = pa.WORKERS_HELD
BOUVARD_BASE_YEARS = (1781, 1821)  # D3(i): "C" is fitted to this part only
SIGMA_SYN_ARCSEC = 1.0  # D3(i), (ii): uniform weights on the noise-free longitudes
NEAR_DEG = 3.0  # § 5: "near", "reproduces the shortfall", "past" thresholds
BROOKES_SIN_DEG = (29, 30, 31, 32)
BROOKES_TABLE_X_TRUE_LON = {29: 323.8, 30: 327.1, 31: 330.4, 32: 333.2}  # λ′ at discovery
BROOKES_MEAN_LON_1846 = {29: 323.1, 30: 324.3, 31: 325.5, 32: 326.7}  # his conclusion, l′
STEP3_HELD_VARIANTS = ("h1_g", "h2_g", "h1_t", "h2_t")
PLACED = ("adams_printed", "brookes_1970", "C2", "D1a", "D1b", "D1c")


# ------------------------------------------------------------------ epochs and frames


def jd_of_year(y: float) -> float:
    """Julian-year convention: 1846.762 is JD 2395574.82 (1846 October 5.3 TDB; S5 and S6 were
    carried at 1846-10-06.0 in step 3, 0.7 day later, 0.004° of Neptune's motion)."""
    return 2451545.0 + (y - 2000.0) * 365.25


JD_1810 = jd_of_year(A.EPOCH_YEAR)
JD_1846 = jd_of_year(A.REPORT_YEAR)


def lon_of_date_deg(pos_ecl_j2000: np.ndarray, jd: float) -> float:
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date

    p = precession_matrix_ecl_j2000_to_ecl_of_date(jd) @ np.asarray(pos_ecl_j2000)
    return math.degrees(math.atan2(p[1], p[0])) % 360.0


def lon_j2000_to_of_date(lon_deg: float, jd: float) -> float:
    """A longitude in the ecliptic of J2000 (zero latitude) carried to the ecliptic of date."""
    lo = math.radians(lon_deg)
    return lon_of_date_deg(np.array([math.cos(lo), math.sin(lo), 0.0]), jd)


def _wrap(d: float) -> float:
    return (d + 180.0) % 360.0 - 180.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def memoir_hashes() -> dict[str, str]:
    return {f: sha256(MEMOIR_DIR / f"adams_memoir_{f}.csv") for f in MEMOIR_FILES}


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with DRIVER_LOG.open("a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def _json_default(o):
    if isinstance(o, A.Solution):
        return o.as_dict()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def _dump(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=1, ensure_ascii=False, default=_json_default) + "\n",
        encoding="utf-8",
    )


# ------------------------------------------------------------------ R2, R3: stages B and C


def stage_b(m: A.Memoir) -> dict:
    out = {"b1": A.b1_table(m), "b2": A.b2_table(m), "b3": {}, "b4": {}, "b5": {}}
    for hyp in (1, 2):
        b3 = A.b3_table(m, hyp)
        out["b3"][f"H{hyp}"] = {
            "rows": b3,
            "coef_max": max(abs(r["delta"]) for r in b3 if r["term"] != "lhs"),
            "lhs_max": max(abs(r["delta"]) for r in b3 if r["term"] == "lhs"),
            "outside_0_0005": [r for r in b3 if r["term"] != "lhs" and abs(r["delta"]) > 0.0005],
        }
        b4 = A.b4_table(m, hyp)
        out["b4"][f"H{hyp}"] = {
            "rows": b4,
            "coef_max": max(abs(r["delta"]) for r in b4 if r["term"] not in ("lhs", "de", "dn")),
            "outside_0_0005": [
                r for r in b4 if r["term"] not in ("lhs", "de", "dn") and abs(r["delta"]) > 0.0005
            ],
        }
        out["b5"][f"H{hyp}"] = {str(k): v for k, v in A.b5_table(m, hyp).items()}
    out["b1_max_log10"] = max(abs(r["delta"]) for r in out["b1"])
    out["b2_max_arcsec"] = max(
        abs(r["delta"]) for r in out["b2"] if not r["item"].startswith("§28")
    )
    return out


def stage_c(m: A.Memoir) -> dict:
    c = A.c_table(m)
    for rec in c.values():
        rec["solution"] = rec["solution"]  # already a dict
        rec["true_lon_1846_of_date"] = rec["true_lon_1846"]
    return c


def locate(m: A.Memoir) -> dict:
    """R3's second part: where C1 and C2 part (the sensitivity table of both hypotheses)."""
    return {f"H{hyp}": A.sensitivity_table(m, hyp) for hyp in (1, 2)}


# ------------------------------------------------------------------ R4: D1


def _d1_record(m: A.Memoir, hyp: int, scheme: str, s: A.Setup, rows, printed) -> dict:
    w = A.weights_for(m, hyp, scheme, rows)
    r = A.solve_least_squares(rows, w, s, printed)
    sol = r["solution"]
    ch = A.chain_with_print_fallback(
        m, hyp, (A.run_hypothesis_1 if hyp == 1 else A.run_hypothesis_2)(m, "C0")
    )
    his_x = A.x_of(ch.solution, A.element_corrections(ch, m, ch.solution))
    rec = {
        "scheme": scheme,
        "weights": {
            "modern": float(w[[k == "modern" for _, k, _, _ in rows].index(True)]),
            "ancient": float(w[[k == "ancient" for _, k, _, _ in rows].index(True)]),
        },
        "x": r["x"],
        "sigma": r["sigma"],
        "sigma_derived": {
            "mean_lon_1846_deg": r["sigma"]["theta"],
            "varpi_prime_deg": r["sigma"]["beta_prime"],
            "eccentricity": r["sigma"]["e_prime_units"] * A.ECC_UNIT,
            "m_prime_rel": (r["sigma"]["m_prime"] / sol.m_prime if sol.m_prime else None),
        },
        "solution": sol.as_dict(),
        "true_lon_1846": r["true_lon_1846"],
        "chi2": r["chi2"],
        "dof": r["dof"],
        "n_rows": r["n_rows"],
        "rms_arcsec": r["rms_arcsec"],
        "his_solution_chi2_same_weights": A.chi2_of(rows, w, s, his_x),
        "vs_print": A.tier_of(sol, A.printed_target(m, hyp)),
        "vs_brookes": A.tier_of(sol, A.BROOKES_1970[hyp]),
        "mean_lon_minus_printed_deg": _wrap(sol.mean_lon_1846 - printed["mean_lon_1846"]),
        "true_lon_minus_printed_deg": _wrap(
            r["true_lon_1846"]
            - A.true_longitude(
                printed["mean_lon_1846"], printed["eccentricity"], printed["varpi_prime_1846"]
            )
        ),
        "n_distinct_minima": r["n_distinct_minima"],
        "residuals": r["residuals"],
    }
    return rec


def d1(m: A.Memoir) -> dict:
    out = {}
    for hyp in (1, 2):
        s = A.setup_from_print(m, hyp)
        rows = A.place_rows(m, s)
        printed = A.printed_solution(m, hyp)
        out[f"H{hyp}"] = {
            "alpha": s.alpha,
            "printed_residual_rms_arcsec": A.residual_rms(m, hyp),
            "printed": printed
            | {
                "true_lon_1846": A.true_longitude(
                    printed["mean_lon_1846"], printed["eccentricity"], printed["varpi_prime_1846"]
                )
            },
            "schemes": {sch: _d1_record(m, hyp, sch, s, rows, printed) for sch in "abc"},
        }
        c = out[f"H{hyp}"]["schemes"]
        out[f"H{hyp}"]["a_minus_c_true_lon_deg"] = _wrap(
            c["a"]["true_lon_1846"] - c["c"]["true_lon_1846"]
        )
        out[f"H{hyp}"]["b_minus_c_true_lon_deg"] = _wrap(
            c["b"]["true_lon_1846"] - c["c"]["true_lon_1846"]
        )
    return out


# ------------------------------------------------------------------ R5: D2


def _profile_point(m: A.Memoir, alpha: float, rows_t, w: np.ndarray, printed, drop: bool) -> dict:
    s = A.setup_at_alpha(m, alpha, 1, drop_second_order=drop)
    rows = [(y, k, t, A.condition_equation(t, c, s.th)) for y, k, t, c in rows_t]
    r = A.solve_least_squares(rows, w, s, printed)
    sol = r["solution"]
    return {
        "alpha": alpha,
        "a_au": 19.18 / alpha,
        "n_prime_over_n": s.th.n_prime_over_n,
        "flags": A.alpha_flags(alpha),
        "chi2": r["chi2"],
        "rms_arcsec": r["rms_arcsec"],
        "theta": sol.theta,
        "mean_lon_1846": sol.mean_lon_1846,
        "eccentricity": sol.eccentricity,
        "m_prime": sol.m_prime,
        "mass_msun": sol.mass,
        "varpi_prime_1846": sol.varpi_prime_1846,
        "true_lon_1846": r["true_lon_1846"],
        "sigma_theta": r["sigma"]["theta"],
        "c_p3e": s.th.c["c_p3e"],
        "c_sin2": s.th.c["c_sin2"],
        "n_distinct_minima": r["n_distinct_minima"],
    }


def d2(m: A.Memoir) -> dict:
    """The α profile with D1(c)'s weights (Hypothesis I's § 52 residuals), the equations rebuilt
    at each α from the fixed dates; the same with the two second-order constants dropped."""
    printed = A.printed_solution(m, 1)
    rows_t = [
        (p["year"], p["kind"], float(p["t_h1"]), p["oc_mean"])
        for p in m.places
        if p["year"] != 1690
    ]
    w = A.weights_for(m, 1, "c", [(y, k, t, None) for y, k, t, _ in rows_t])
    grid = A.alpha_grid()
    prof = [_profile_point(m, a, rows_t, w, printed, False) for a in grid]
    prof_drop = [_profile_point(m, a, rows_t, w, printed, True) for a in grid]
    points = {
        str(d): _profile_point(m, math.sin(math.radians(d)), rows_t, w, printed, False)
        for d in BROOKES_SIN_DEG
    }
    for d, rec in points.items():
        rec["brookes_table_x_true_lon"] = BROOKES_TABLE_X_TRUE_LON[int(d)]
        rec["brookes_mean_lon_1846"] = BROOKES_MEAN_LON_1846[int(d)]
        rec["true_lon_minus_brookes"] = _wrap(
            rec["true_lon_1846"] - BROOKES_TABLE_X_TRUE_LON[int(d)]
        )
        rec["mean_lon_minus_brookes"] = _wrap(rec["mean_lon_1846"] - BROOKES_MEAN_LON_1846[int(d)])

    def minimum(points_):
        clean = [p for p in points_ if not p["flags"]]
        pm = min(points_, key=lambda p: p["chi2"])
        pc = min(clean, key=lambda p: p["chi2"])
        return {
            "alpha": pm["alpha"],
            "chi2": pm["chi2"],
            "flags": pm["flags"],
            "alpha_outside_bands": pc["alpha"],
            "chi2_outside_bands": pc["chi2"],
        }

    def at(points_, alpha):
        return next(p for p in points_ if abs(p["alpha"] - alpha) < 1e-9)

    mn, mnd = minimum(prof), minimum(prof_drop)
    return {
        "grid": grid,
        "weights": {"modern": float(w[0]), "ancient": float(w[-1])},
        "profile": prof,
        "profile_second_order_dropped": prof_drop,
        "minimum": mn,
        "minimum_second_order_dropped": mnd,
        "delta_chi2_at_0_50": at(prof, 0.5)["chi2"] - mn["chi2_outside_bands"],
        "delta_chi2_at_0_515": at(prof, 0.515)["chi2"] - mn["chi2_outside_bands"]
        if any(abs(p["alpha"] - 0.515) < 1e-9 for p in prof)
        else None,
        "poles": A.POLES,
        "bands": {
            "truncation_limited": "alpha >= 0.58 (plan)",
            "pole_limited_3_1": "|alpha - 0.4807| <= 0.015 (found when the grid ran: the "
            "i = 3 eccentric terms carry the divisor n - 3n', zero at alpha = 3^(-2/3))",
        },
        "second_order_constants": "carried as the printed Hypothesis I values scaled by "
        "c_p3e(alpha)/c_p3e(0.5); the dropped-terms profile is the sensitivity to that choice",
        "points_at_sin": points,
    }


# ------------------------------------------------------------------ R6: D3


def truth():
    """The project's model with Uranus and the perturber from the cached Horizons states at the
    integrator's epoch: the 13-entry θ (inclined perturber), the context and Uranus's (i, Ω)."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.model import bodies
    from fork.model.forward import (
        ForwardContext,
        inclined_perturber_theta_from_state,
        uranus_inc_raan_at_epoch,
        uranus_theta_from_state,
    )

    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    th_u = uranus_theta_from_state(about("uranus", [A1_EPOCH_JD_TDB])[0])
    th_p = inclined_perturber_theta_from_state(
        about("neptune", [A1_EPOCH_JD_TDB])[0], bodies.NEPTUNE.mass_msun, inc_raan
    )
    return np.concatenate([th_u, th_p]), ctx, inc_raan


def perturber_lon_1846(theta: np.ndarray, ctx) -> dict:
    """The perturber's heliocentric true longitude at 1846.762 in both frames, and its distance."""
    from fork.model.forward import predict_states_ecl_j2000
    from fork.model.hidden_mass import I_PERTURBER

    pos, _ = predict_states_ecl_j2000(theta, [JD_1846], ctx)
    p = pos[0, I_PERTURBER]
    return {
        "ecl_j2000": math.degrees(math.atan2(p[1], p[0])) % 360.0,
        "of_date": lon_of_date_deg(p, JD_1846),
        "r_au": float(np.linalg.norm(p)),
    }


def place_epochs(m: A.Memoir, hyp: int = 1) -> dict[int, float]:
    return {
        p["year"]: jd_of_year(A.EPOCH_YEAR + A.UNIT_YR * float(p[f"t_h{hyp}"])) for p in m.places
    }


def uranus_mean_longitude_at(theta_u: np.ndarray, jd: float, ctx) -> float:
    """Uranus's mean longitude (ecl_j2000) at ``jd`` from a Uranus-only orbit: its integrated
    true longitude converted with the orbit's e and ϖ (Bouvard's ε is the tables' mean longitude
    at the epoch)."""
    from fork.model.forward import predict_uranus_lon

    nu = float(predict_uranus_lon(theta_u, [jd], ctx)[0]) - theta_u[4]
    e = theta_u[1]
    E = 2.0 * math.atan2(math.sqrt(1 - e) * math.sin(nu / 2), math.sqrt(1 + e) * math.cos(nu / 2))
    M = E - e * math.sin(E)
    return math.degrees(M + theta_u[4]) % 360.0


def synthetic_normal_places(m: A.Memoir) -> dict:
    """D3(i)'s data: Uranus's heliocentric true longitude at his 30 epochs from the truth; "C"
    the same model without the perturber fitted to 1781–1821; O − C converted by r²/(ab) to
    mean-longitude corrections."""
    from fork.data.synthetic import LongitudeObservations
    from fork.inference.fit import fit_uranus_elements
    from fork.model.forward import predict_states_ecl_j2000, predict_uranus_lon
    from fork.model.observation import I_URANUS

    th, ctx, inc_raan = truth()
    epochs = place_epochs(m)
    years = list(epochs)
    jds = np.array([epochs[y] for y in years])
    lon_true = predict_uranus_lon(th, jds, ctx)
    base = np.array([BOUVARD_BASE_YEARS[0] <= y <= BOUVARD_BASE_YEARS[1] for y in years])
    sig = math.radians(SIGMA_SYN_ARCSEC / 3600.0)
    obs = LongitudeObservations(
        jds[base],
        lon_true[base],
        np.full(int(base.sum()), sig),
        np.array(["synthetic"] * int(base.sum())),
    )
    fit_c = fit_uranus_elements(obs, ctx, inc_raan)
    lon_c = predict_uranus_lon(fit_c.theta, jds, ctx)
    pos, _ = predict_states_ecl_j2000(fit_c.theta, jds, ctx)
    r = np.linalg.norm(pos[:, I_URANUS], axis=1)
    a, e = fit_c.theta[0], fit_c.theta[1]
    oc_true = ((lon_true - lon_c + math.pi) % (2 * math.pi) - math.pi) * A.ARCSEC_PER_RAD
    oc_mean = oc_true * r**2 / (a * a * math.sqrt(1 - e * e))
    eps_c = uranus_mean_longitude_at(fit_c.theta, JD_1810, ctx)
    return {
        "years": years,
        "jd": jds,
        "lon_true_rad": lon_true,
        "oc_true_arcsec": dict(zip(years, oc_true.tolist(), strict=True)),
        "oc_mean_arcsec": dict(zip(years, oc_mean.tolist(), strict=True)),
        "c_fit": {
            "theta": fit_c.theta,
            "chi2": fit_c.chi2,
            "n": int(base.sum()),
            "rms_arcsec": math.degrees(math.sqrt(np.mean(fit_c.residual_rad**2))) * 3600,
            "a_au": float(a),
            "e": float(e),
            "varpi_deg": math.degrees(fit_c.theta[4]) % 360,
            "eps_1810_ecl_j2000_deg": eps_c,
            "beta_deg": (eps_c - math.degrees(fit_c.theta[4])) % 360.0,
        },
        "truth_theta": th,
        "ctx": ctx,
        "inc_raan": inc_raan,
    }


def _his_model_longitudes(
    theta: float, e: float, beta_prime: float, eps_c: float, n_prime_deg: float
) -> dict:
    """His solution's perturber longitude at 1846.762 in the data's frame (ecl_j2000: the mean
    longitude ε′ + 12 n′ without his 30′ of precession) and of date."""
    eps_prime = (eps_c - theta) % 360.0
    mean = (eps_prime + 12.0 * n_prime_deg) % 360.0
    varpi = (eps_c - beta_prime) % 360.0
    true = A.true_longitude(mean, e, varpi)
    return {
        "eps_prime_1810": eps_prime,
        "mean_lon_1846_ecl_j2000": mean,
        "varpi_prime_ecl_j2000": varpi,
        "true_lon_1846_ecl_j2000": true,
        "true_lon_1846_of_date": lon_j2000_to_of_date(true, JD_1846),
    }


def d3_i(m: A.Memoir, syn: dict) -> dict:
    """(i) his layout, his method: D1(c) on the synthetic normal places, and his own procedure
    (the C2 chain with the § 9 column replaced)."""
    s = A.setup_from_print(m, 1)
    rows = A.place_rows(m, s, oc=syn["oc_mean_arcsec"])
    w = A.weights_for(m, 1, "c", rows)
    r = A.solve_least_squares(rows, w, s, A.printed_solution(m, 1))
    sol = r["solution"]
    eps_c = syn["c_fit"]["eps_1810_ecl_j2000_deg"]
    out = {
        "d1c": {
            "x": r["x"],
            "sigma": r["sigma"],
            "solution": sol.as_dict(),
            "chi2": r["chi2"],
            "rms_arcsec": r["rms_arcsec"],
            "n_distinct_minima": r["n_distinct_minima"],
            "longitudes": _his_model_longitudes(
                sol.theta, sol.eccentricity, sol.beta_prime, eps_c, s.n_prime_deg
            ),
        }
    }
    # his own procedure on the same numbers
    ms = copy.deepcopy(m)
    for p in ms.places:
        p["oc_mean"] = syn["oc_mean_arcsec"][p["year"]]
    ch = A.run_hypothesis_1(ms, "C2")
    roots = []
    for root in ch.roots:
        if root["m_prime"] <= 0:
            continue
        th_, _S, _e, beta, k, r_cos, r_sin, _R, eps, n_prime, prec = A._setup(ms, 1, "print")
        cs = A.complete_solution(
            1,
            root["theta"],
            root["P"],
            root["Q"],
            th_.c,
            ch.steps[(22, "1")],
            k,
            r_cos,
            r_sin,
            beta,
            eps,
            n_prime,
            prec,
        )
        corr = A.element_corrections(ch, ms, cs)
        x = A.x_of(cs, corr)
        roots.append(
            {
                "theta": cs.theta,
                "m_prime": cs.m_prime,
                "eccentricity": cs.eccentricity,
                "beta_prime": cs.beta_prime,
                "chi2_c_weights": A.chi2_of(rows, w, s, x),
                "solution": cs,
            }
        )
    out["his_procedure"] = {"n_roots": len(ch.roots), "positive_mass_roots": len(roots)}
    if roots:
        best = min(roots, key=lambda q: q["chi2_c_weights"])
        cs = best["solution"]
        out["his_procedure"] |= {
            "solution": cs.as_dict(),
            "chi2_c_weights": best["chi2_c_weights"],
            "roots": [{k: v for k, v in q.items() if k != "solution"} for q in roots],
            "longitudes": _his_model_longitudes(
                cs.theta, cs.eccentricity, cs.beta_prime, eps_c, s.n_prime_deg
            ),
            "picked_by": "smallest chi2 under D1(c)'s weights among the positive-mass roots",
        }
    return out


def d3_ii(m: A.Memoir, syn: dict) -> dict:
    """(ii) his layout, the project's model: the 29 noise-free heliocentric longitudes fitted by
    ``fit_perturber_multistart`` with a held at 38.4 AU and the headline's λ₀ starts."""
    from fork.data.synthetic import LongitudeObservations
    from fork.inference.fit import fit_perturber_multistart

    keep = np.array([y != 1690 for y in syn["years"]])
    sig = math.radians(SIGMA_SYN_ARCSEC / 3600.0)
    obs = LongitudeObservations(
        syn["jd"][keep],
        syn["lon_true_rad"][keep],
        np.full(int(keep.sum()), sig),
        np.array(["synthetic"] * int(keep.sum())),
    )
    res = fit_perturber_multistart(
        obs,
        syn["ctx"],
        syn["inc_raan"],
        a_fixed_au=HELD_A_AU,
        n_workers=min(8, os.cpu_count() or 1),
    )
    best = res.best
    th = best.theta
    lon = perturber_lon_1846(th, syn["ctx"])
    return {
        "n_obs": int(keep.sum()),
        "sigma_arcsec": SIGMA_SYN_ARCSEC,
        "status_counts": res.status_counts(),
        "minima": [
            {
                "chi2": mn.best.chi2,
                "n_starts": mn.n_starts,
                "mass_earth": mn.best.mass_earth,
                "e": float(mn.best.theta[8]),
                "lam0_deg": math.degrees(mn.best.theta[10]) % 360,
                "true_lon_1846": perturber_lon_1846(mn.best.theta, syn["ctx"]),
            }
            for mn in res.minima
        ],
        "best": {
            "theta": th,
            "chi2": best.chi2,
            "dof": best.dof,
            "mass_earth": best.mass_earth,
            "e": float(th[8]),
            "varpi_deg": math.degrees(th[9]) % 360,
            "rms_arcsec": math.degrees(math.sqrt(np.mean(best.residual_rad**2))) * 3600,
            "longitudes": lon,
        },
        "uranus_only_chi2": res.uranus_only.chi2,
    }


def d3_cheap(m: A.Memoir) -> dict:
    syn = synthetic_normal_places(m)
    T = perturber_lon_1846(syn["truth_theta"], syn["ctx"])
    i, ii = d3_i(m, syn), d3_ii(m, syn)
    rec = {
        "truth": {
            "theta": syn["truth_theta"],
            "neptune_1846": T,
            "note": "Uranus and Neptune from the cached Horizons states at the A1 epoch; "
            "group offsets zero; no noise",
        },
        "synthetic_places": {
            "years": syn["years"],
            "jd": syn["jd"],
            "oc_true_arcsec": syn["oc_true_arcsec"],
            "oc_mean_arcsec": syn["oc_mean_arcsec"],
            "his_oc_mean_arcsec": {p["year"]: p["oc_mean"] for p in m.places},
        },
        "c_fit": syn["c_fit"],
        "i": i,
        "ii": ii,
    }
    for key, lon in (
        ("i_d1c", i["d1c"]["longitudes"]["true_lon_1846_ecl_j2000"]),
        (
            "i_his_procedure",
            i["his_procedure"].get("longitudes", {}).get("true_lon_1846_ecl_j2000"),
        ),
        ("ii", ii["best"]["longitudes"]["ecl_j2000"]),
    ):
        rec[f"O_{key}"] = None if lon is None else _wrap(lon - T["ecl_j2000"])
    return rec


# ------------------------------------------------------------------ D3(iii): the nested job


def d3iii_dataset(with_box: bool = True):
    """Set D's layout with synthetic RA and Dec from the truth (``obs.predict`` of the step-3
    observation object with the values substituted), H1(G)'s error model, the Uranus box by the
    M5 rule on the synthetic data. Returns (obs, ctx, inc_raan, pq, box, truth θ)."""
    from fork.inference.fit_radec import uranus_box_radec
    from fork.model.observation import synthetic_values

    obs, ctx, inc_raan, pq, _box, _test = pa.dataset("h1_g")
    th, _ctx2, _ = truth()
    val, _off = synthetic_values(obs, obs.predict(th, ctx), np.random.default_rng(0), noise=False)
    obs_syn = obs.with_values(val)
    box = uranus_box_radec(obs_syn, ctx, inc_raan) if with_box else None
    return obs_syn, ctx, inc_raan, pq, box, th


def d3iii_space(cell: int, pq, box):
    from fork.inference.priors import PosteriorSpace

    held = PosteriorSpace(
        (None, math.log(HELD_A_AU), None, None, None),
        e_max=pa.E_MAX,
        pq=pq,
        uranus_box=box,
        inclined=True,
    )
    return held.with_e_slice(*pa.E_CELLS_HELD[cell])


@dataclass(frozen=True)
class D3Job:
    cell: int

    @property
    def id(self) -> str:
        return f"d3iii/c{self.cell}"

    def path(self, suffix: str) -> Path:
        return CACHE / f"d3iii__c{self.cell}{suffix}"

    def rng_seed(self) -> list[int]:
        return [*SEED_NAMESPACE, pa.ERRORS.index("g"), self.cell + 1]


def d3iii_jobs() -> list[D3Job]:
    return [D3Job(k) for k in range(len(pa.E_CELLS_HELD))]


def do_d3iii_job(job: D3Job, n_workers: int, maxiter: int | None = None) -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult, run_nested

    pa.check_frozen()
    out = job.path(".pkl")
    if out.exists():
        print(f"{job.id}: cached")
        return
    obs, ctx, _inc, pq, box, _th = d3iii_dataset()
    space = d3iii_space(job.cell, pq, box)
    resumed = job.path(".ckpt").exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(
        obs,
        ctx,
        space,
        np.random.default_rng(job.rng_seed()),
        nlive=NLIVE,
        n_workers=n_workers,
        maxiter=maxiter,
        progress_log=PROGRESS_LOG,
        checkpoint_file=None if maxiter else job.path(".ckpt"),
        label=f"p1b_adams_repro/{job.id}",
    )
    if maxiter is not None:
        print(
            f"{job.id}: smoke run, {res.niter} it, {res.ncall} calls, {res.wall_s:.1f} s, "
            f"ndim {space.ndim}, names {space.names}"
        )
        return
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append(
        f"p1b_adams_repro/{job.id}",
        "nested",
        n_likelihood_calls=res.ncall,
        n_integrations=res.ncall,
        sim_years=res.ncall * yrs,
        wall_s=time.time() - t0,
        cpu_s=ledger.cpu_now() - c0,
        n_workers=n_workers,
        estimated=resumed,
        resumed=resumed,
        nlive=NLIVE,
        niter=res.niter,
    )
    print(
        f"{job.id}: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, "
        f"{res.wall_s / 60:.1f} min"
    )


def run_d3iii(n_workers: int) -> None:
    """Both cells in parallel (step 3's scheduler: one retry, resuming from the checkpoint)."""
    pa.check_frozen()
    pending = [j for j in d3iii_jobs() if not j.path(".pkl").exists()]
    _log(f"run d3iii: {len(pending)} jobs pending, budget {n_workers} workers")
    running: dict[str, tuple[subprocess.Popen, float, int]] = {}
    failures: dict[str, int] = {}
    per = max(1, min(WORKERS_PER_CELL, n_workers // max(1, len(pending) or 1)))
    while pending or running:
        for jid, (proc, t0, _) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[jid]
            job = D3Job(int(jid.rsplit("c", 1)[1]))
            if rc == 0 and job.path(".pkl").exists():
                _log(f"done {jid} ({(time.time() - t0) / 60:.1f} min)")
            else:
                failures[jid] = failures.get(jid, 0) + 1
                _log(f"FAILED {jid} (exit {rc}, attempt {failures[jid]})")
                if failures[jid] < 2:
                    pending.insert(0, job)
        busy = sum(w for *_, w in running.values())
        while pending and (not running or busy + per <= n_workers):
            job = pending.pop(0)
            log = job.path(".log")
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a") as fh:
                proc = subprocess.Popen(
                    [
                        sys.executable,
                        "-m",
                        "fork.export.p1b_adams_repro",
                        "job",
                        str(job.cell),
                        "--n-workers",
                        str(per),
                    ],
                    stdout=fh,
                    stderr=subprocess.STDOUT,
                    cwd=ROOT,
                )
            running[job.id] = (proc, time.time(), per)
            busy += per
            _log(f"start {job.id} (nlive {NLIVE}, {per} workers, pid {proc.pid})")
        time.sleep(30)
    _log(f"run d3iii finished; failures: {failures or 'none'}")


def d3iii_summary() -> dict | None:
    """Median and 68 % of the perturber's heliocentric true longitude at 1846.762 over the
    equal-weight samples (step 3's subset size), the best-fit sample beside; O_iii against the
    truth's Neptune from the same model."""
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
    from fork.model.bodies import M_EARTH_MSUN
    from fork.model.forward import predict_states_ecl_j2000
    from fork.model.hidden_mass import I_PERTURBER

    jobs = d3iii_jobs()
    if not all(j.path(".pkl").exists() for j in jobs):
        return None
    res = [NestedResult.load(j.path(".pkl")) for j in jobs]
    obs, ctx, _inc, _pq, _box, th_truth = d3iii_dataset(with_box=False)
    T = perturber_lon_1846(th_truth, ctx)
    rng = np.random.default_rng([*SUMMARY_RNG, D3III_RNG_INDEX])
    part = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    eq = part.equal_weight(rng)
    space = res[0].space
    idx = rng.choice(len(eq), size=min(pa.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    lon_j2000 = np.empty(len(idx))
    lon_date = np.empty(len(idx))
    thetas = np.array([space.theta(v) for v in eq])
    for i, v in enumerate(eq[idx]):
        pos, _ = predict_states_ecl_j2000(space.theta(v), [JD_1846], ctx)
        p = pos[0, I_PERTURBER]
        lon_j2000[i] = math.degrees(math.atan2(p[1], p[0])) % 360.0
        lon_date[i] = lon_of_date_deg(p, JD_1846)
    best = max(res, key=lambda r: r.max_logl)
    th_best = space.theta(best.samples[int(np.argmax(best.results["logl"]))])
    best_lon = perturber_lon_1846(th_best, ctx)
    resid = obs.residuals_arcsec(*obs.predict(th_best, ctx))
    summ = m6c.lon_summary(lon_date, {"neptune_same_model": T["of_date"]})
    truth_vs = summ["vs"]["neptune_same_model"]
    return {
        "jobs": [
            {
                "id": j.id,
                "seed": j.rng_seed(),
                "nlive": r.nlive,
                "logz": r.logz,
                "logzerr": r.logzerr,
                "niter": r.niter,
                "ncall": r.ncall,
                "wall_min": r.wall_s / 60,
                "max_logl": r.max_logl,
            }
            for j, r in zip(jobs, res, strict=True)
        ],
        "evidence": {"logz": part.logz, "logzerr": part.logzerr, "shares": part.shares.tolist()},
        "n_samples": len(eq),
        "n_derived": len(idx),
        "true_lon_1846_of_date": summ,
        "true_lon_1846_ecl_j2000": m6c._quantiles(lon_j2000, 3),
        "O_iii_median": truth_vs["median_minus_target_deg"],
        "O_iii_68": [_wrap(summ["lo68"] - T["of_date"]), _wrap(summ["hi68"] - T["of_date"])],
        "best_fit": {
            "theta": th_best,
            "minus2_lnl": -2 * float(best.max_logl),
            "longitudes": best_lon,
            "O_iii_best": _wrap(best_lon["of_date"] - T["of_date"]),
            "mass_earth": float(10.0 ** th_best[6]),
            "e": float(th_best[8]),
            "inc_rel_uranus_deg": math.degrees(th_best[11]),
            "rms_arcsec": float(np.sqrt(np.mean(resid**2))),
        },
        "mass_earth": m6c._quantiles(10.0 ** thetas[:, 6], 3),
        "mass_vs_truth_earth": float(10.0 ** th_truth[6]),
        "mass_msun": m6c._quantiles(10.0 ** thetas[:, 6] * M_EARTH_MSUN, 8),
        "e": m6c._quantiles(thetas[:, 8], 4),
        "e_truth": float(th_truth[8]),
        "inc_rel_uranus_deg": m6c._quantiles(np.degrees(thetas[:, 11]), 3),
        "neptune_same_model": T,
        "box": "the M5 rule on the synthetic data (the Uranus-only RA/Dec fit)",
    }


# ------------------------------------------------------------------ R7: D4


def d4(m: A.Memoir, c_table: dict, d1_table: dict, T_date: float) -> dict:
    """Step 3's fixed-a posteriors re-derived at 1846.762 from the cached samples (the same
    subset and draws as step 3's summary), with the solutions placed by the step-3 rule."""
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
    from fork.model.forward import predict_states_ecl_j2000
    from fork.model.hidden_mass import I_PERTURBER

    out = {"epoch_jd_tdb": JD_1846, "neptune_same_model_of_date": T_date, "variants": {}}
    ctx = None
    for v in STEP3_HELD_VARIANTS:
        jobs = pa.all_jobs((v,))
        if not all(j.path(".pkl").exists() for j in jobs):
            out["variants"][v] = None
            continue
        if ctx is None:
            _obs, ctx, *_rest = pa.dataset("h1_g")
        res = [NestedResult.load(j.path(".pkl")) for j in jobs]
        rng = np.random.default_rng([*SUMMARY_RNG, pa.VARIANTS.index(v)])
        part = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
        eq = part.equal_weight(rng)
        space = res[0].space
        idx = rng.choice(len(eq), size=min(pa.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
        lon = np.empty(len(idx))
        for i, q in enumerate(eq[idx]):
            pos, _ = predict_states_ecl_j2000(space.theta(q), [JD_1846], ctx)
            lon[i] = lon_of_date_deg(pos[0, I_PERTURBER], JD_1846)
        hyp = 1 if v.startswith("h1") else 2
        H = f"H{hyp}"
        pr = d1_table[H]["printed"]
        br = A.BROOKES_1970[hyp]
        targets = {
            "adams_printed": pr["true_lon_1846"],
            "brookes_1970": A.true_longitude(
                br["mean_lon_1846"], br["eccentricity"], br["varpi_prime_1846"]
            ),
            "C2": c_table[f"{H}/C2"]["true_lon_1846"],
            "D1a": d1_table[H]["schemes"]["a"]["true_lon_1846"],
            "D1b": d1_table[H]["schemes"]["b"]["true_lon_1846"],
            "D1c": d1_table[H]["schemes"]["c"]["true_lon_1846"],
            "neptune_same_model": T_date,
        }
        summ = m6c.lon_summary(lon, targets)
        placements = {k: pa.status_of(summ["vs"][k]) for k in targets}
        out["variants"][v] = {
            "hypothesis": H,
            "n_derived": len(idx),
            "true_lon_1846_of_date": summ,
            "placements": placements,
            "gap_printed_minus_median_deg": -summ["vs"]["adams_printed"]["median_minus_target_deg"],
            "shortfall_median_minus_neptune_deg": summ["vs"]["neptune_same_model"][
                "median_minus_target_deg"
            ],
        }
    return out


# ------------------------------------------------------------------ R8: the readings of § 5


def q1_reading(c_table: dict, hyp: int) -> dict:
    H = f"H{hyp}"
    t = {s: c_table[f"{H}/{s}"]["vs_print"]["tier"] for s in ("C0", "C1", "C2")}
    tb = c_table[f"{H}/C2"]["vs_brookes"]["tier"]
    if all(v in (1, 2) for v in t.values()):
        row, text = (
            1,
            (
                "His printed solution follows from his printed data by his printed "
                "method; his arithmetic is reproducible. Proceed to Q2, Q3."
            ),
        )
    elif t["C0"] in (1, 2) and t["C2"] == 0:
        row = 2
        if tb in (1, 2):
            text = (
                "His data and method as printed do not give his answer; C2 is within tier "
                f"{tb} of Brookes's figures, so his re-derivation is confirmed and the "
                "located equation is the Q1 result: the (n) normal equation of the second "
                "group (its coefficients, rounded sums of t sin mt, through the "
                "ill-conditioned equation (2)). Q2 and Q3 run on C2's chain with Adams's "
                "printed solution carried beside it."
            )
        else:
            text = (
                "Not reproduced from the data and near neither Adams nor Brookes; Q2, Q3 "
                "run on his printed equations (C0's inputs), said in every table."
            )
    else:
        row, text = (
            3,
            (
                "C0 not reproduced: a transcription error in (1)–(4) and §§ 29/46 or a "
                "misreading of §§ 28/45; the step stops at Q1."
            ),
        )
    return {"row": row, "tiers_vs_print": t, "C2_tier_vs_brookes": tb, "text": text}


def q2_reading(d1_h: dict, d2_table: dict, step3_median_true_lon: float | None) -> dict:
    c = d1_h["schemes"]["c"]
    his_true = d1_h["printed"]["true_lon_1846"]
    a_c = abs(d1_h["a_minus_c_true_lon_deg"])
    b_c = abs(d1_h["b_minus_c_true_lon_deg"])
    a_his = abs(d1_h["schemes"]["a"]["true_lon_minus_printed_deg"])
    b_his = abs(d1_h["schemes"]["b"]["true_lon_minus_printed_deg"])
    move = c["true_lon_minus_printed_deg"]
    proc = {
        "c_tier_vs_print": c["vs_print"]["tier"],
        "a_minus_his_deg": a_his,
        "b_minus_his_deg": b_his,
        "a_minus_c_deg": a_c,
        "b_minus_c_deg": b_c,
        "c_minus_his_deg": move,
    }
    if c["vs_print"]["tier"] in (1, 2) and a_his <= NEAR_DEG and b_his <= NEAR_DEG:
        proc["row"], proc["text"] = (
            1,
            (
                "His shortcut did not move the answer: the procedure is "
                "not the cause of the step-3 gap."
            ),
        )
    elif c["vs_print"]["tier"] == 0 and step3_median_true_lon is not None:
        gap = _wrap(his_true - step3_median_true_lon)
        toward = (abs(move) >= NEAR_DEG) and (math.copysign(1, move) == -math.copysign(1, gap))
        if toward:
            proc["row"], proc["share"] = 2, move / -gap
            proc["text"] = (
                "The procedure accounts for a share of the gap, stated as "
                "(his − D1(c)) / (his − step-3 median)."
            )
        else:
            proc["row"], proc["text"] = (
                3,
                ("The procedure matters but does not explain the gap; reported, no share."),
            )
    else:
        proc["row"], proc["text"] = (
            1,
            (
                "D1(c) within tier 2 of his solution; (a) or (b) more "
                "than 3° from it: see the weighting row."
            ),
        )
    if a_c > NEAR_DEG or b_c > NEAR_DEG:
        proc["weighting_row"] = (
            "(a) or (b) more than 3° from (c): the weighting matters; (c)'s "
            "conditioning on his residuals is material and the three are "
            "reported as a bracket, no single procedure reading."
        )
    mn = d2_table["minimum"]
    d = d2_table["delta_chi2_at_0_50"]
    amin = mn["alpha_outside_bands"]
    if d >= 4 and amin <= 0.57:
        drow, dtext = (
            1,
            (
                "His data in his model disfavour his distance within the range the "
                f"series can reach; the minimum is at α = {amin:.2f}; Neptune's α is "
                "beyond the pole."
            ),
        )
    elif d <= 1:
        drow, dtext = (
            2,
            (
                "They do not; the preference step 3 found is not in his data as he "
                "reduced them, or not visible to his model."
            ),
        )
    else:
        drow, dtext = 3, "Weak or truncation-limited; stated, no statement."
    return {
        "procedure": proc,
        "distance": {
            "row": drow,
            "delta_chi2_at_0_50": d,
            "alpha_min": amin,
            "alpha_min_any": mn["alpha"],
            "flags_at_min_any": mn["flags"],
            "text": dtext,
        },
    }


def _classify(off: float, O3: float) -> str:
    if abs(off) <= NEAR_DEG:
        return "near"
    if off > NEAR_DEG:
        return "past"
    if abs(off - O3) <= NEAR_DEG:
        return "reproduces"
    if off < O3 - NEAR_DEG:
        return "beyond"
    return "between"


def q3_reading(O_i: float, O_ii: float, O_iii: float, O3: float) -> dict:
    """The § 5 Q3 table, rows in order; the shares as the table states them."""
    ci, cii, ciii = _classify(O_i, O3), _classify(O_ii, O3), _classify(O_iii, O3)
    out = {
        "O_i": O_i,
        "O_ii": O_ii,
        "O_iii": O_iii,
        "O_3": O3,
        "class": {"i": ci, "ii": cii, "iii": ciii},
        "rows": [],
    }

    def split_on_ii(prefix: str) -> str:
        if cii == "reproduces":
            return prefix + " Split on (ii): (ii) reproduces the shortfall → his model."
        if cii == "near":
            return prefix + " Split on (ii): (ii) near the truth → his layout, not his model."
        sl, sm = (O_iii - O_ii) / O_iii, (O_ii - O_i) / O_iii
        out["shares"] = {"layout": sl, "model": sm}
        return prefix + (
            f" Split on (ii): between → shares {sl:.2f} to the layout and {sm:.2f} to the model."
        )

    if ciii == "near":
        out["rows"].append(
            {
                "row": 1,
                "text": (
                    "The step-3 shortfall is in set D (this reconstruction of his data), not in "
                    "the assumed distance; A1–A5 are the candidates. (i) reported as a check: "
                    + (
                        "his method also finds the planet on clean data"
                        if ci == "near"
                        else f"his method does not find the planet on clean data ({ci}), "
                        "which does not change this row"
                    )
                    + "."
                ),
            }
        )
    elif ciii == "reproduces" and ci == "near":
        out["rows"].append(
            {
                "row": 2,
                "text": split_on_ii(
                    "The shortfall is a property of the assumed distance in the full model, not "
                    "of set D, and his method on his layout does not share it."
                ),
            }
        )
    elif ciii == "reproduces" and ci == "reproduces":
        out["rows"].append(
            {
                "row": 3,
                "text": (
                    "Both methods fall short on clean data at that distance: the shortfall is a "
                    "property of the assumed distance, and his nearness on real data owes to his "
                    "data or to compensating errors (Brookes's correction moves 325° 7′ to "
                    "323° 39′, away from the truth). A1–A5 and his reductions remain the "
                    "candidates."
                ),
            }
        )
    elif ciii == "reproduces" and ci == "between":
        share = (O_iii - O_i) / O_iii
        out["shares"] = {"layout_and_model": share}
        out["rows"].append(
            {
                "row": 4,
                "text": split_on_ii(
                    f"Partial: his layout and model together remove {share:.2f} of the "
                    "shortfall and his data the rest."
                ),
            }
        )
    elif ciii == "between":
        out["rows"].append(
            {
                "row": 6,
                "share_distance": O_iii / O3,
                "text": (
                    f"Partial: the assumed distance explains {O_iii / O3:.2f} of the shortfall "
                    "and set D the rest; (i) read against (iii) as in the rows above."
                ),
            }
        )
    elif ciii in ("beyond", "past"):
        out["rows"].append(
            {
                "row": 7,
                "text": (
                    "The noise-free full-model fit does not behave as step 3's did; reported "
                    "with the numbers, no attribution."
                ),
            }
        )
    if ci == "past":
        out["rows"].append(
            {
                "row": 5,
                "bias_deg": O_i,
                "text": (
                    f"His method/model carries a positive bias of {O_i:+.2f}° on clean data at "
                    "38.4 AU, compared with his real-data solutions' +2 to +3° beyond Neptune "
                    "(step 3 § 6)"
                    + (
                        ": the two agree within 3°, so his overshoot on real data is the method's"
                        if abs(O_i - 2.5) <= NEAR_DEG
                        else ""
                    )
                    + "."
                ),
            }
        )
    return out


# ------------------------------------------------------------------ the figure


def draw(payload: dict, path: Path) -> None:
    """Four panels: his § 9 data with the reproduced curve (C2) and D1(c)'s for each hypothesis;
    the D2 profile with the 3:1 pole inside the grid and the 2:1 pole marked at the edge; the D3
    offsets. Two hues (step 3's), ink for data; one series per role."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from fork.export.p1b_adams_plot import BLUE, GRID, INK, MUTED, ORANGE, _style

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6))
    fig.subplots_adjust(hspace=0.42, wspace=0.5)
    places = payload["his_data"]
    for j, hyp in enumerate((1, 2)):
        ax = axes[0, j]
        _style(ax)
        yrs = [p["year"] for p in places]
        ax.scatter(
            yrs,
            [p["oc_mean"] for p in places],
            s=14,
            color=INK,
            zorder=3,
            label="his § 9 residuals",
        )
        for key, color, ls, label in (
            ("C2", BLUE, "-", "C2 (his chain from the data)"),
            ("D1c", ORANGE, "--", "D1(c) least squares"),
        ):
            curve = payload["curves"][f"H{hyp}"][key]
            ax.plot(
                [c["year"] for c in curve],
                [c["total"] for c in curve],
                color=color,
                linestyle=ls,
                linewidth=1.6,
                label=label,
            )
        ax.set_title(
            f"Hypothesis {'I' if hyp == 1 else 'II'}: elements + perturbation, arcsec",
            fontsize=9.5,
            color=INK,
        )
        ax.set_xlabel("year", fontsize=8.5)
        if j == 0:
            ax.legend(fontsize=7.5, frameon=False)
    ax = axes[1, 0]
    _style(ax)
    prof = payload["d2"]["profile"]
    clean = [p for p in prof if not p["flags"]]
    flagged = [p for p in prof if p["flags"]]
    ax.plot([p["alpha"] for p in prof], [p["chi2"] for p in prof], color=GRID, linewidth=1.0)
    ax.plot(
        [p["alpha"] for p in clean],
        [p["chi2"] for p in clean],
        color=BLUE,
        marker="o",
        markersize=4,
        linestyle="none",
        label="χ², D1(c) weights",
    )
    ax.plot(
        [p["alpha"] for p in flagged],
        [p["chi2"] for p in flagged],
        color=BLUE,
        marker="o",
        markersize=4,
        markerfacecolor="white",
        linestyle="none",
        label="inside a band (pole or truncation)",
    )
    drop = payload["d2"]["profile_second_order_dropped"]
    ax.plot(
        [p["alpha"] for p in drop],
        [p["chi2"] for p in drop],
        color=ORANGE,
        linewidth=1.2,
        linestyle=":",
        label="second-order terms dropped",
    )
    ymax = max(p["chi2"] for p in clean) * 1.9
    ax.set_ylim(0, ymax)
    for a_, lab, dx in (
        (A.POLES["3:1"], "3:1 pole", -0.004),
        (0.5, "α = 0.5", 0.002),
        (0.515, "0.515", 0.002),
    ):
        ax.axvline(a_, color=MUTED, linewidth=0.8, linestyle="--")
        ax.text(
            a_ + dx,
            ymax * 0.62,
            lab,
            fontsize=7.5,
            color=MUTED,
            ha="right" if dx < 0 else "left",
            va="top",
            rotation=90,
        )
    ax.axvspan(0.58, 0.605, color=GRID, alpha=0.5, lw=0)
    ax.text(
        0.603,
        ymax * 0.03,
        "truncation-limited;\n2:1 pole at 0.630 →",
        fontsize=7,
        color=MUTED,
        ha="right",
        va="bottom",
    )
    ax.set_xlabel("α = a / a′", fontsize=8.5)
    ax.set_title(
        "D2: his equations rebuilt at each α, Hypothesis I weights", fontsize=9.5, color=INK
    )
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    ax = axes[1, 1]
    _style(ax)
    d3 = payload["d3"]
    rows = [
        ("(i) his layout, his method (D1(c))", d3["O_i_d1c"], None),
        ("(i) his own procedure", d3.get("O_i_his_procedure"), None),
        ("(ii) his layout, project model", d3["O_ii"], None),
    ]
    if d3.get("iii"):
        rows.append(
            (
                "(iii) set D layout, project sampler",
                d3["iii"]["O_iii_median"],
                d3["iii"]["O_iii_68"],
            )
        )
    ys = list(range(len(rows)))[::-1]
    for y, (_lab, off, band) in zip(ys, rows, strict=True):
        if off is None:
            continue
        if band is not None:
            ax.plot(band, [y, y], color=BLUE, linewidth=4, alpha=0.35, solid_capstyle="butt")
        ax.plot([off], [y], marker="o", color=BLUE, markersize=6)
        ax.text(off, y + 0.18, f"{off:+.1f}°", fontsize=7.5, color=INK, ha="center")
    ax.axvline(0, color=INK, linewidth=0.9)
    O3 = d3.get("O_3_step3")
    if O3 is not None:
        ax.axvline(O3, color=ORANGE, linewidth=1.2, linestyle="--")
        ax.text(
            O3 + 0.15,
            -0.55,
            f"step-3 shortfall {O3:+.1f}°",
            fontsize=7.5,
            color=ORANGE,
            ha="left",
            va="bottom",
        )
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlabel("longitude at 1846.762 − Neptune's, degrees (one truth)", fontsize=8.5)
    ax.set_title("D3: data against method", fontsize=9.5, color=INK)
    ax.set_ylim(-0.9, len(rows) - 0.3)
    ax.set_xlim(min(-12.5, (O3 or 0) - 2), 1.5)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ orchestration


def run_cheap() -> dict:
    """B, C, the location, D1, D2, D3 (i) and (ii): minutes; cached as JSON."""
    from fork import ledger

    t0, c0 = time.time(), ledger.cpu_now()
    m = A.load_memoir()
    out = {
        "memoir_sha256": memoir_hashes(),
        "stage_b": stage_b(m),
        "stage_c": stage_c(m),
        "locate": locate(m),
    }
    _log("B, C, locate done")
    out["d1"] = d1(m)
    _log("D1 done")
    out["d2"] = d2(m)
    _log("D2 done")
    out["d3"] = d3_cheap(m)
    _log("D3 (i), (ii) done")
    ledger.append(
        "p1b_adams_repro/cheap",
        "fit",
        n_likelihood_calls=None,
        n_integrations=None,
        sim_years=None,
        wall_s=time.time() - t0,
        cpu_s=ledger.cpu_now() - c0,
        n_workers=1,
        estimated=True,
        note="B, C, D1, D2 (2 × 21 + 4 least-squares profiles), D3 (i), (ii)",
    )
    _dump(CHEAP_JSON, out)
    return out


def curves(m: A.Memoir, c_table_raw: dict, d1_table: dict) -> dict:
    """For the figure: elements + perturbation at every normal place for C2's solution and
    D1(c)'s, both hypotheses."""
    out = {}
    for hyp in (1, 2):
        H = f"H{hyp}"
        run = A.run_hypothesis_1 if hyp == 1 else A.run_hypothesis_2
        ch = A.chain_with_print_fallback(m, hyp, run(m, "C2"))
        corr = A.element_corrections(ch, m, ch.solution)
        res = A.residuals(m, hyp, ch.solution, corr)
        s = A.setup_from_print(m, hyp)
        x = np.array([d1_table[H]["schemes"]["c"]["x"][k] for k in A.PARAMS])
        rows = A.place_rows(m, s, drop_1690=False)
        out[H] = {
            "C2": [{"year": r["year"], "total": r["total"]} for r in res],
            "D1c": [{"year": y, "total": A.model_value(eq, t, x, s)} for y, _k, t, eq in rows],
        }
    return out


def summarize(png: bool = True) -> dict:
    if not CHEAP_JSON.exists():
        raise SystemExit("run `cheap` first")
    cheap = json.loads(CHEAP_JSON.read_text(encoding="utf-8"))
    m = A.load_memoir()
    iii = d3iii_summary()
    if iii is not None:
        _dump(D3III_JSON, iii)
    T_date = cheap["d3"]["truth"]["neptune_1846"]["of_date"]
    d4_table = d4(m, cheap["stage_c"], cheap["d1"], T_date)
    h1g = d4_table["variants"].get("h1_g")
    step3_median = h1g["true_lon_1846_of_date"]["median"] if h1g else None
    O3 = h1g["shortfall_median_minus_neptune_deg"] if h1g else None
    d3 = dict(cheap["d3"])
    d3["iii"] = iii
    d3["O_3_step3"] = O3
    readings = {
        "Q1": {f"H{h}": q1_reading(cheap["stage_c"], h) for h in (1, 2)},
        "Q2": {
            f"H{h}": q2_reading(cheap["d1"][f"H{h}"], cheap["d2"], step3_median) for h in (1, 2)
        },
        "Q3": (
            q3_reading(d3["O_i_d1c"], d3["O_ii"], iii["O_iii_median"], O3)
            if iii is not None and O3 is not None
            else None
        ),
        "Q3_with_his_procedure": (
            q3_reading(d3["O_i_his_procedure"], d3["O_ii"], iii["O_iii_median"], O3)
            if iii is not None and O3 is not None and d3.get("O_i_his_procedure") is not None
            else None
        ),
        "Q3_best_fit_sample": (
            q3_reading(d3["O_i_d1c"], d3["O_ii"], iii["best_fit"]["O_iii_best"], O3)
            if iii is not None and O3 is not None
            else None
        ),
    }
    payload = {
        "description": "Phase 1b step 3b: Adams's 1846 solution reproduced from the memoir "
        "(docs/notes/P1B_adams_reproduction.md). Report order R1–R10 of § 6.",
        "epoch": {
            "year": A.REPORT_YEAR,
            "jd_tdb": JD_1846,
            "note": "heliocentric true longitude at 1846.762, of date unless a key says "
            "ecl_j2000; his orbits lie in the ecliptic, step 3's in Uranus's plane",
        },
        "R1_transcription": {
            "memoir_sha256": cheap["memoir_sha256"],
            "passes": "two independent passes, 1988 common values, 15 disagreements, all decided "
            "on the named render (docs/notes/adams_memoir_reconciliation.md)",
            "secondary": {
                "brookes_1970": "obtained (ADS), read in full, hash in MANIFEST.md",
                "smart_1947": "sought once, not obtained (print only)",
            },
        },
        "R2_stage_b": cheap["stage_b"],
        "R3_stage_c": cheap["stage_c"],
        "R3_locate": cheap["locate"],
        "R4_d1": cheap["d1"],
        "R5_d2": cheap["d2"],
        "R6_d3": d3,
        "R7_d4": d4_table,
        "R8_readings": readings,
        "R9_unseparated": (
            "His data fitted by the project's model would need Bouvard's 1821 tables with his "
            "§ 7 corrections at his 30 epochs (not transcribed). D2 cannot reach Neptune's "
            "distance: the 2:1 pole at α = 0.630 bounds the series, and the 3:1 pole at α = 0.481 "
            "found when the grid ran sits 0.02 inside Hypothesis I's α."
        ),
        "R10_prior_work": (
            "Brookes (1970) re-derived Hypothesis I from the memoir's equations (0.1509, 0.8267, "
            "311° 35′, 323° 39′) and attributed the difference from Adams to an unnamed numerical "
            "error; C2 confirms his figures and locates the difference in the (n) normal equation "
            "of the second group. His Table X and conclusion give the α dependence at sin 29°–32°; "
            "D2's points sit beside them."
        ),
        "his_data": [
            {"year": p["year"], "kind": p["kind"], "t": p["t_h1"], "oc_mean": p["oc_mean"]}
            for p in m.places
        ],
        "curves": curves(m, cheap["stage_c"], cheap["d1"]),
        "d2": cheap["d2"],
        "d3": d3,
    }
    _dump(OUT_JSON, payload)
    if png:
        draw(payload, OUT_PNG)
    print(json.dumps(readings, indent=1, ensure_ascii=False, default=_json_default))
    print(f"wrote {OUT_JSON}")
    return payload


def plan() -> None:
    for j in d3iii_jobs():
        state = (
            "done"
            if j.path(".pkl").exists()
            else "resumable"
            if j.path(".ckpt").exists()
            else "pending"
        )
        print(f"{j.id:10s} nlive {NLIVE}  workers {WORKERS_PER_CELL}  seed {j.rng_seed()}  {state}")
    print(f"cheap stages cached: {CHEAP_JSON.exists()}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    sub.add_parser("cheap")
    p_run = sub.add_parser("d3iii")
    p_run.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    p_job = sub.add_parser("job")
    p_job.add_argument("cell", type=int)
    p_job.add_argument("--n-workers", type=int, default=4)
    p_job.add_argument("--maxiter", type=int, default=None, help="smoke run: not saved")
    p_sum = sub.add_parser("summarize")
    p_sum.add_argument("--no-png", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "plan":
        plan()
    elif args.cmd == "cheap":
        run_cheap()
    elif args.cmd == "d3iii":
        run_d3iii(args.n_workers)
    elif args.cmd == "job":
        do_d3iii_job(D3Job(args.cell), args.n_workers, args.maxiter)
    else:
        summarize(png=not args.no_png)


if __name__ == "__main__":
    main()
