"""M10: profile likelihood in the semi-major axis under Student-t errors (docs/notes/M10.md,
pre-registered 2026-11-10, amended on approval the same day).

    uv run python -m fork.export.m10 profile --dataset all [--n-workers 12]
    uv run python -m fork.export.m10 report      # docs/notes/M10_profile.json and .png

A diagnostic: it changes no published result. Every setting is M10.md's; what M10.md leaves to
M9 is imported from ``fork.export.m9`` (grid, reference points, start counts, continuation,
candidate and reading-rule constants and functions). Nothing in ``m9`` is edited.

**What changes against M9:**

- **The likelihood:** r_i − b_g ~ t_ν(scale s_g), ν = 4, with S3's scales (the frozen error
  model's ``sigma_t4``) and the offsets' Gaussian prior τ_g. The profiled quantity is
  Q = Σ_i (ν + 1) ln(1 + (r_i − b_g)²/(ν s_g²)) + Σ_g (b_g/τ_g)², minimised over the perturber
  (a held), Uranus's four in-plane elements and the eight offsets. Δ ln L = −ΔQ/2. S3's marginal
  (``t_marginal_m2lnl``, offsets integrated) is evaluated at every optimum and reported beside Q.
- **The optimiser:** ``least_squares`` (trf, bounds, the e disk by ``_fit_bounded``) with the
  Student-t loss ρ(z) = (ν + 1) ln(1 + z/ν) on f_i = (r_i − b̂_g)/s_g; b̂ by IRLS at each
  evaluation (``t_offset_modes``); each group's prior term as one pseudo-residual f_g with
  ρ(f_g²) = (b̂_g/τ_g)². scipy's cost is then ½ Q. Ten variables: the six perturber parameters
  and Uranus's four (M9 projected Uranus out through the Gaussian marginal).
- **The search model:** M8's first-order Uranus response for the perturber, and Uranus's own
  elements linear about u_ref through the arcsec Jacobian of the Uranus-only model.
- **Convergence, pre-registered:** a polish is converged when scipy's status is positive and Q
  moved by less than 0.01 over its last five iterations; otherwise it restarts from its end,
  up to four rounds, every round recorded.
- **Six data sets** (``DATASETS``): all, post-1781, and all minus each pre-1781 observer. Seeds
  ``[2026, 10, d, k]``.
- **Pre-run changes:** the Uranus box from two starts (``uranus_box_radec_with_fits``), and the
  cache-only Horizons fetcher (``FORK_HORIZONS_CACHE_ONLY``), which the launch sets.

Records keep M9's key ``chi2`` for Q (the χ² scale), so M9's continuation, candidate and reading
helpers apply unchanged.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import least_squares

from fork.export import m6c, m9
from fork.export.m9 import (
    A_GRID_AU,
    CONT_TOL_CHI2,
    DISTINCT_E,
    DISTINCT_LON_DEG,
    DLNL_THRESHOLD,
    E_MAX,
    FLAT_RANGE_AU,
    INC_MAX_DEG,
    MAX_NFEV,
    MAX_PASSES,
    MAX_POLISH,
    MIN_EFF_SAMPLES,
    N_GRID,
    N_NEAREST,
    N_POSTERIOR,
    N_RANDOM,
    N_RERUN,
    PENALTY,
    POINTS_AU,
    POLISH_WINDOW_CHI2,
    REFERENCE_AU,
    STATUS_CODES,
    TIGHT_AU,
    U_FLAG_SIGMA,
    X_AU,
    _inc_raan_pq,
    _log,
    _path,
    _read,
    _start_of,
    _write,
    candidates,
    crossings,
    parabola_max,
    phi_of,
    posterior_marginal_a,
    posterior_starts,
    read_profile,
    replacement_posterior,
)
from fork.export.to_json import A1_EPOCH_JD_TDB
from fork.model.student_t import t_marginal_m2lnl, t_offset_modes

ROOT = m6c.ROOT
CACHE = ROOT / "data" / "cache" / "m10"
NOTES = ROOT / "docs" / "notes"
OUT_JSON = NOTES / "M10_profile.json"
OUT_PNG = NOTES / "M10_profile.png"
M9_JSON = NOTES / "M9_profile.json"
U_REF_FILE = m9.U_REF_FILE
CONJ_JSON = m9.CONJ_JSON

FROZEN = {  # M10.md, "Frozen inputs (SHA-256)"
    "data/historical/leverrier_1846.csv":
        "2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4",
    "data/historical/m6_error_model.json":
        "0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277",
    "docs/notes/M8_u_ref.json":
        "4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5",
    "docs/notes/M9_profile.json":
        "adbb9317b0e790b400996f7b64697e8fb44f3db8c3261e2235dca1050113883b",
}  # fmt: skip

NU = 4.0  # S3's ν, fixed
SEED_NAMESPACE = (2026, 10)  # M9's was (2026, 9)
DATASETS = ("all", "post1781", "no_flamsteed", "no_bradley", "no_mayer", "no_lemonnier")
OBSERVER = {"no_flamsteed": "Flamsteed", "no_bradley": "Bradley", "no_mayer": "Mayer",
            "no_lemonnier": "Le Monnier"}  # fmt: skip
LABEL = {"all": "A, all observations", "post1781": "B, the observations from 1781 on",
         "no_flamsteed": "A − Flamsteed", "no_bradley": "A − Bradley", "no_mayer": "A − Mayer",
         "no_lemonnier": "A − Le Monnier"}  # fmt: skip
EXPECTED_MEASUREMENTS = {"all": 529, "post1781": 492, "no_flamsteed": 519, "no_bradley": 528,
                         "no_mayer": 527, "no_lemonnier": 505}  # fmt: skip
S3_VARIANT = "s3_t4"
MAX_ROUNDS = 4  # polish restarts (M10.md, "Convergence, pre-registered")
MARGINAL_SPREAD_TOL = 1.0  # (marginal − Q) varying by more than this across a grid is reported
BORDERLINE_AU = 0.2  # amendment 1: |95 % width − 8 AU| ≤ this is "borderline"
PSEUDO_CAP = 700.0  # see t_residuals
SCREENED_ROWS = ("LV1846-T77-p131-l13-ra", "LV1846-T77-p135-l03-ra", "LV1846-T77-p135-l08-ra")


# ---------------------------------------------------------------------------- problem


def check_frozen() -> None:
    from fork.export import m8_runs

    m8_runs.check_frozen("headline")
    for rel, digest in FROZEN.items():
        if m6c.sha256(ROOT / rel) != digest:
            raise RuntimeError(f"{rel} differs from M10's frozen input")


def t_scales() -> np.ndarray:
    """S3's scales s_g in ``GROUPS`` order (the frozen error model's ``sigma_t4``)."""
    from fork.model.observation import GROUPS

    t = m6c.t_scales_arcsec()
    return np.array([t[g] for g in GROUPS], dtype=float)


def observations(dataset: str):
    """The data set's measurements with the headline's frozen (Gaussian) error model: the t
    scales ride beside it in the worker. Rows are dropped as in M9 (``m9.drop_measurements``);
    an observer's rows are those whose CSV ``observer`` column names him."""
    from fork.data.calibration import real_observations
    from fork.model import observation as om

    if dataset not in DATASETS:
        raise ValueError(f"dataset must be one of {DATASETS}")
    if dataset in ("all", "post1781"):
        return m9.observations(dataset)
    rows = om.read_rows()
    em = om.ErrorModel.load(ROOT / m6c.VARIANT_CONFIG["headline"]["error_model"])
    obs = real_observations(rows, em)
    if obs.n_measurements != len(rows):
        raise RuntimeError("measurements are not in row order")
    keep = np.array([r["observer"] != OBSERVER[dataset] for r in rows])
    return m9.drop_measurements(obs, keep)


def row_ids(dataset: str) -> list[str]:
    """The CSV ids of the data set's measurements, in measurement order."""
    from fork.model import observation as om

    rows = om.read_rows()
    if dataset == "all":
        return [r["id"] for r in rows]
    if dataset == "post1781":
        return [r["id"] for r in rows if r["era"] != "ancient"]
    return [r["id"] for r in rows if r["observer"] != OBSERVER[dataset]]


def u_reference() -> tuple[np.ndarray, np.ndarray]:
    u = json.loads(U_REF_FILE.read_text())["headline"]
    return np.array(u["u_ref"], dtype=float), np.array(u["u_std"], dtype=float)


def problem(dataset: str):
    """(obs, ctx, space, PerturbativeProblem, u_ref, u_std, box fits). The Uranus box is built
    from two starts (pre-run change (a)); u_ref and its posterior std are the headline's."""
    from fork.export import m8_runs
    from fork.inference.fit_radec import uranus_box_radec_with_fits
    from fork.inference.priors import free_a_space
    from fork.model.forward import ForwardContext
    from fork.model.perturbative import PerturbativeProblem

    check_frozen()
    obs = observations(dataset)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan, pq = _inc_raan_pq()
    u_ref, u_std = u_reference()
    box, fits = uranus_box_radec_with_fits(obs, ctx, inc_raan, extra_starts=[u_ref])
    space = free_a_space(pq, E_MAX, box, inclined=True)
    pr = PerturbativeProblem.build(obs, ctx, space, u_ref, step_days=m8_runs.STEP_DAYS,
                                   quadrature=m8_runs.QUADRATURE)  # fmt: skip
    return obs, ctx, space, pr, u_ref, u_std, fits


# ---------------------------------------------------------------------------- the objective


def t_loss(z: np.ndarray) -> np.ndarray:
    """scipy's ``loss`` for Student-t errors: ρ(z) = (ν + 1) ln(1 + z/ν) on z = f², with its
    first two derivatives. ½ Σ ρ(f_i²) is then ½ Q (the data part), exactly."""
    z = np.asarray(z, dtype=float)
    rho = (NU + 1.0) * np.log1p(z / NU)
    d1 = (NU + 1.0) / (NU + z)
    d2 = -(NU + 1.0) / (NU + z) ** 2
    return np.vstack([rho, d1, d2])


def t_residuals(resid_arcsec: np.ndarray, group_index: np.ndarray, scale: np.ndarray,
                tau: np.ndarray, nu: float = NU) -> np.ndarray:  # fmt: skip
    """f = (r − b̂_g)/s_g with b̂ the IRLS offset modes, followed by one pseudo-residual per
    group, f_g = √(ν (exp((b̂_g/τ_g)²/(ν + 1)) − 1)), so that ρ(f_g²) = (b̂_g/τ_g)² and
    Σ ρ(f²) over the whole vector is Q."""
    b = t_offset_modes(resid_arcsec, group_index, scale, tau, nu)
    f = (resid_arcsec - b[group_index]) / scale[group_index]
    # exp overflows past x ≈ 709, so the mapping is capped at PSEUDO_CAP: beyond |b/τ| ≈ 59
    # (120″ modern, 300″ ancient; no optimum is within a factor 20 of it) the prior term is
    # held at (ν + 1) × PSEUDO_CAP ≈ 3500 instead of becoming infinite, which scipy refuses
    # at an initial point. Any other non-finite entry becomes PENALTY.
    x = np.minimum((b / tau) ** 2 / (nu + 1.0), PSEUDO_CAP)
    pseudo = np.sqrt(nu * np.expm1(x))
    out = np.concatenate([f, pseudo])
    return np.where(np.isfinite(out), out, PENALTY)


def q_of(f: np.ndarray, nu: float = NU) -> float:
    """Q = Σ (ν + 1) ln(1 + f²/ν) over the vector ``t_residuals`` returns."""
    return float(np.sum((nu + 1.0) * np.log1p(np.asarray(f, dtype=float) ** 2 / nu)))


def irls_weights(resid_arcsec: np.ndarray, group_index: np.ndarray, scale: np.ndarray,
                 tau: np.ndarray, nu: float = NU) -> np.ndarray:  # fmt: skip
    """(ν + 1)/(ν + f²) per measurement at the IRLS offsets: 1 for a residual well inside the
    scale, small for one the t model discounts."""
    f = t_residuals(resid_arcsec, group_index, scale, tau, nu)[: resid_arcsec.size]
    return (nu + 1.0) / (nu + f**2)


# ---------------------------------------------------------------------------- starts


def random_starts(d: int, k: int, n: int = N_RANDOM, rerun: bool = False) -> list[dict]:
    """M9's draws from the headline prior with the M10 seed namespace [2026, 10, d, k] (+1 for
    the rerun). Uranus starts at the headline's u_ref."""
    from fork.inference.priors import PRIOR_LOG10_M_EARTH

    rng = np.random.default_rng([*SEED_NAMESPACE, d, k] + ([1] if rerun else []))
    u_ref = u_reference()[0].tolist()
    out = []
    for i in range(n):
        log10_m = rng.uniform(*PRIOR_LOG10_M_EARTH)
        e, varpi = E_MAX * rng.random(), 2.0 * math.pi * rng.random()
        lam0 = 2.0 * math.pi * rng.random()
        cos_i = 1.0 - (1.0 - math.cos(math.radians(INC_MAX_DEG))) * rng.random()
        node = 2.0 * math.pi * rng.random()
        out.append({"log10_m": float(log10_m), "k": math.sqrt(e) * math.cos(varpi),
                    "h": math.sqrt(e) * math.sin(varpi), "lam0": float(lam0),
                    "inc": math.acos(cos_i), "node": float(node), "u": u_ref,
                    "origin": f"{'rerun' if rerun else 'random'}:{i}"})  # fmt: skip
    return out


# ---------------------------------------------------------------------------- fits (workers)

_W: dict = {}


def _init_worker(dataset: str) -> None:
    from fork.inference.fit import uranus_x_to_theta
    from fork.inference.profile import GalleGeometry
    from fork.model.linearized import FD_REL_STEP, fd_jacobian

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    obs, ctx, space, pr, u_ref, u_std, box_fits = problem(dataset)

    def r_u(u: np.ndarray) -> np.ndarray:  # the Uranus-only model's residuals, arcsec
        theta = uranus_x_to_theta(np.concatenate([u, space.pq]))
        return obs.residuals_arcsec(*obs.predict(theta, ctx))

    j_u = fd_jacobian(r_u, u_ref, FD_REL_STEP * pr.sigma_raw, "central")
    geo = GalleGeometry.from_horizons()
    conj = json.loads(CONJ_JSON.read_text())["galle_night_jd_tdb"]
    if abs(geo.jd_tdb - conj) > 1e-9:
        raise RuntimeError("Galle's night differs between the geometry and m6c_conjunction.json")
    _W.update(dataset=dataset, obs=obs, ctx=ctx, space=space, pr=pr, u_ref=u_ref, u_std=u_std,
              j_u=j_u, geo=geo, scale=t_scales(), tau=np.asarray(obs.offset_prior_arcsec),
              g=np.asarray(obs.group_index), box_fits=box_fits,
              calls={"pert": 0, "exact": 0})  # fmt: skip


def _n_out() -> int:
    return _W["obs"].n_measurements + _W["scale"].size


def _pert_resid_arcsec(phi: np.ndarray, u: np.ndarray, count: bool = True) -> np.ndarray | None:
    """Observed − model in arcsec at (φ, u) from the first-order search model: Uranus's
    response to the perturber at u_ref, plus the Uranus-only Jacobian times (u − u_ref)."""
    from fork.model.apparent import apparent_ra_dec

    _W["calls"]["pert"] += int(count)
    pr, obs = _W["pr"], _W["obs"]
    try:
        dx = pr.delta_x(phi, second_order=False)
        k = pr.obs_node
        ra, dec = apparent_ra_dec(pr.x_ref[k] + dx, pr.v_ref[k], obs.geometry)
        r = obs.residuals_arcsec(ra, dec) + _W["j_u"] @ (np.asarray(u) - _W["u_ref"])
    except (ValueError, ZeroDivisionError, FloatingPointError, OverflowError):
        return None
    return r if np.all(np.isfinite(r)) else None


def _exact_resid_arcsec(q: np.ndarray, count: bool = True) -> np.ndarray | None:
    _W["calls"]["exact"] += int(count)
    obs, space = _W["obs"], _W["space"]
    if q[0] <= 0.0 or math.hypot(q[1], q[2]) >= 0.95:
        return None
    try:
        r = obs.residuals_arcsec(*obs.predict(space.theta(q), _W["ctx"]))
    except (ValueError, ZeroDivisionError, FloatingPointError, OverflowError):
        return None
    return r if np.all(np.isfinite(r)) else None


def _f_of(r: np.ndarray | None) -> np.ndarray:
    if r is None:
        return np.full(_n_out(), PENALTY)
    return t_residuals(r, _W["g"], _W["scale"], _W["tau"])


def fit_at(a_au: float, start: dict, exact: bool, trace: bool = False) -> dict:
    """One bounded Student-t fit at the held a from ``start`` (log10_m, k, h, lam0, inc, node,
    u). As ``m9.fit_at`` with the t loss, Uranus's four elements always variables, and Q in
    place of χ². With ``trace``, Q at every accepted iterate, the final gradient norm and
    scipy's optimality measure are recorded."""
    from fork.inference.fit import _fit_bounded, _prior_box, le_verrier_space

    sp = le_verrier_space(a_au, E_MAX)
    inc_max = math.radians(INC_MAX_DEG)
    inc0 = np.array([start["inc"], start["node"]], dtype=float)
    u0 = np.asarray(start["u"], dtype=float)
    c0 = dict(_W["calls"])

    def lm(layout, z, warm):
        z = np.asarray(z, dtype=float)
        lo_p, hi_p = _prior_box(layout)
        inc_s = warm.inc if warm is not None else inc0
        x0 = np.concatenate([z[:4], inc_s, z[4:]])
        lo = np.concatenate([np.full(4, -np.inf), [0.0, -np.inf], lo_p])
        hi = np.concatenate([np.full(4, np.inf), [inc_max, np.inf], hi_p])
        x0 = np.clip(x0, lo, hi)
        seen: dict[bytes, float] = {}

        def fun(x):
            phi = phi_of(layout.y_from_z(x[6:]), x[4:6])
            if exact:
                r = _exact_resid_arcsec(np.concatenate([x[:4], phi]))
            else:
                r = _pert_resid_arcsec(phi, x[:4])
            f = _f_of(r)
            if trace:
                seen[np.asarray(x, dtype=float).tobytes()] = q_of(f)
            return f

        hist: list[float] = []

        def cb(ir):
            x = np.asarray(getattr(ir, "x", ir), dtype=float)
            hist.append(seen.get(x.tobytes()) or q_of(fun(x)))

        extra = {"callback": cb} if trace else {}
        if trace:
            hist.append(q_of(fun(x0)))
        sol = least_squares(fun, x0, method="trf", loss=t_loss, bounds=(lo, hi), x_scale="jac",
                            xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=MAX_NFEV,
                            **extra)  # fmt: skip
        ns = SimpleNamespace(x=np.concatenate([sol.x[:4], sol.x[6:]]), cost=float(sol.cost),
                             success=bool(sol.success), status=int(sol.status),
                             nfev=int(sol.nfev), inc=sol.x[4:6].copy(),
                             message=str(sol.message))  # fmt: skip
        if trace:
            last5 = hist[-6] - hist[-1] if len(hist) >= 6 else hist[0] - hist[-1]
            ns.trace = {"n_iter": len(hist) - 1, "chi2_by_iter": hist,
                        "chi2_last5_change": float(last5),
                        "grad_norm": float(np.linalg.norm(sol.grad)),
                        "optimality": float(sol.optimality), "njev": int(sol.njev),
                        "message": str(sol.message)}  # fmt: skip
        return ns

    z0 = np.concatenate([u0, sp.z_from_y([start["log10_m"], a_au, start["k"], start["h"],
                                          start["lam0"]])])  # fmt: skip
    best, layout = _fit_bounded(lm, z0, sp)
    rec = _optimum(a_au, layout.y_from_z(best.x[4:]), best.inc, best, layout.on_e_bound, exact,
                   best.x[:4])  # fmt: skip
    rec["origin"] = start.get("origin")
    rec["calls"] = {k: _W["calls"][k] - c0[k] for k in c0}
    if trace:
        rec["trace"] = best.trace
    return rec


def _optimum(a_au, y, inc, sol, on_e_bound: bool, exact: bool, u) -> dict:
    """The fit record. ``chi2`` holds Q (M9's key, so M9's helpers apply); ``marginal_m2lnl``
    is S3's normalised marginal at the same point (offsets integrated)."""
    from fork.inference.priors import PRIOR_LOG10_M_EARTH
    from fork.inference.profile import galle_lon_of_date_rad

    phi = phi_of(y, inc)
    u = np.asarray(u, dtype=float)
    space, box = _W["space"], _W["space"].uranus_box
    theta = space.theta(np.concatenate([u, phi]))
    lon = math.degrees(galle_lon_of_date_rad(theta, _W["ctx"], _W["geo"]))
    if exact:  # bookkeeping evaluations, not counted as search calls
        r = _exact_resid_arcsec(np.concatenate([u, phi]), count=False)
    else:
        r = _pert_resid_arcsec(phi, u, count=False)
    marg = (t_marginal_m2lnl(r, _W["g"], _W["scale"], _W["tau"], NU) if r is not None
            else math.inf)  # fmt: skip
    e = float(y[2] ** 2 + y[3] ** 2)
    bounds = []
    if y[0] <= PRIOR_LOG10_M_EARTH[0] + 1e-6:
        bounds.append("mass_lo")
    if y[0] >= PRIOR_LOG10_M_EARTH[1] - 1e-6:
        bounds.append("mass_hi")
    if on_e_bound:
        bounds.append("e_max")
    if inc[0] >= math.radians(INC_MAX_DEG) - 1e-9:
        bounds.append("inc_max")
    return {"a_au": float(a_au), "chi2": 2.0 * float(sol.cost), "exact": exact,
            "marginal_m2lnl": float(marg),
            "log10_m": float(y[0]), "mass_earth": float(10.0 ** y[0]), "e": e,
            "varpi_deg": math.degrees(math.atan2(y[3], y[2]) % (2.0 * math.pi)),
            "k": float(y[2]), "h": float(y[3]), "lam0": float(y[4] % (2.0 * math.pi)),
            "inc": float(inc[0]), "node": float(inc[1] % (2.0 * math.pi)),
            "inc_deg": math.degrees(inc[0]), "u": u.tolist(), "galle_lon_deg": lon,
            "success": bool(sol.success), "status": int(sol.status), "nfev": int(sol.nfev),
            "on_e_bound": bool(on_e_bound), "active_bounds": bounds,
            "u_box_half_widths": float(np.max(np.abs((u - np.asarray(box.centre))
                                                     / np.asarray(box.half_width)))),
            "u_vs_u_ref_std": float(np.max(np.abs((u - _W["u_ref"]) / _W["u_std"]))),
            }  # fmt: skip


def _task(t: dict) -> dict:
    rec = fit_at(t["a_au"], t["start"], t["exact"], t.get("trace", False))
    rec["point"] = t["point"]
    return rec


def _run(tasks: list[dict], pool) -> list[dict]:
    if pool is None:
        return [_task(t) for t in tasks]
    return list(pool.map(_task, tasks, chunksize=1))


def _polish_rounds(t: dict) -> dict:
    """One candidate's exact polish with the pre-registered convergence rule: restart from the
    end point while Q moved by ≥ ``CONT_TOL_CHI2`` over the last five iterations, up to
    ``MAX_ROUNDS``. Returns the final record with every round summarised."""
    rounds = [fit_at(t["a_au"], t["start"], True, trace=True)]
    while (rounds[-1]["trace"]["chi2_last5_change"] >= CONT_TOL_CHI2
           or rounds[-1]["status"] <= 0) and len(rounds) < MAX_ROUNDS:  # fmt: skip
        s2 = _start_of(rounds[-1]) | {"origin": f"repolish:{t['start']['origin']}"}
        rounds.append(fit_at(t["a_au"], s2, True, trace=True))
    final = dict(rounds[-1])
    final["point"] = t["point"]
    final["origin"] = t["start"]["origin"]
    final["rounds"] = [{"chi2": x["chi2"], "status": x["status"], "nfev": x["nfev"],
                        **{k: v for k, v in x["trace"].items() if k != "chi2_by_iter"}}
                       for x in rounds]  # fmt: skip
    final["converged"] = (final["status"] > 0
                          and final["trace"]["chi2_last5_change"] < CONT_TOL_CHI2)  # fmt: skip
    final["calls"] = {k: sum(x["calls"][k] for x in rounds) for k in rounds[0]["calls"]}
    return final


def _run_polish(tasks: list[dict], pool) -> list[dict]:
    if pool is None:
        return [_polish_rounds(t) for t in tasks]
    return list(pool.map(_polish_rounds, tasks, chunksize=1))


# ---------------------------------------------------------------------------- the profile run


def profile(dataset: str, n_workers: int, cache: Path = CACHE, *, smoke: bool = False) -> dict:
    """The whole M10 profile for one data set, resumable per stage as M9's relaunch: the
    independent starts one point at a time, the continuation, the exact polishes with their
    convergence rounds, and the rerun of a failed point."""
    from multiprocessing import get_context

    from fork import ledger
    from fork.inference.profile import continuation_sweeps

    d = DATASETS.index(dataset)
    log_path = cache / "progress.log"
    points = list(range(len(POINTS_AU)))
    if smoke:
        points = [5, 13]
    t0, c0 = time.time(), ledger.cpu_now()
    _init_worker(dataset)
    calls = {"pert": 0, "exact": 0}

    def account(recs: list[dict]) -> None:
        for r in recs:
            for k in calls:
                calls[k] += r["calls"][k]

    n_meas = _W["obs"].n_measurements
    if not smoke and n_meas != EXPECTED_MEASUREMENTS[dataset]:
        raise RuntimeError(f"{dataset}: {n_meas} measurements, expected "
                           f"{EXPECTED_MEASUREMENTS[dataset]}")  # fmt: skip
    _log(f"{dataset}: {n_meas} measurements, {_W['obs'].obs_id.size} observations, "
         f"{n_workers} workers, smoke {smoke}, cache-only Horizons "
         f"{os.environ.get('FORK_HORIZONS_CACHE_ONLY', '')!r}", log_path)  # fmt: skip
    _log(f"{dataset}: Uranus box fits " + "; ".join(
        f"{f['start']}: χ² {f['chi2']:.1f}/{f['dof']}" + (" (centres the box)" if f["centres_box"]
                                                         else "") for f in _W["box_fits"]),
        log_path)  # fmt: skip
    _write(
        _path(dataset, "box.json", cache),
        {
            "fits": _W["box_fits"],
            "centre": list(_W["space"].uranus_box.centre),
            "half_width": list(_W["space"].uranus_box.half_width),
        },
    )
    pool = (get_context("spawn").Pool(n_workers, initializer=_init_worker, initargs=(dataset,))
            if n_workers > 1 else None)  # fmt: skip
    try:
        q, _, logl = replacement_posterior()
        pending = [k for k in points if not _path(dataset, f"starts_{k}.json", cache).exists()]
        if pending:
            _log(f"{dataset}: independent starts, {len(pending)} points pending", log_path)
        for k in pending:
            n_r, n_p = (2, 1) if smoke else (N_RANDOM, N_POSTERIOR)
            starts = random_starts(d, k, n_r) + posterior_starts(q, logl, POINTS_AU[k], n_p)
            tasks = [{"point": k, "a_au": POINTS_AU[k], "start": s, "exact": False}
                     for s in starts]  # fmt: skip
            recs = _run(tasks, pool)
            account(recs)
            _write(_path(dataset, f"starts_{k}.json", cache),
                   {"point": k, "a_au": POINTS_AU[k], "n_random": n_r, "n_posterior": n_p,
                    "seed": [*SEED_NAMESPACE, d, k], "optima": recs})  # fmt: skip
            ok = [r["chi2"] for r in recs if r["success"]]
            _log(f"{dataset}: point {k} (a = {POINTS_AU[k]} AU) starts done, {len(recs)} fits, "
                 f"best Q {min(ok):.3f}" if ok else
                 f"{dataset}: point {k} (a = {POINTS_AU[k]} AU) starts done, none succeeded",
                 log_path)  # fmt: skip
        optima = {k: _read(_path(dataset, f"starts_{k}.json", cache))["optima"] for k in points}

        cont_path = _path(dataset, "continuation.json", cache)
        cont = _read(cont_path)
        order = sorted(points, key=lambda k: POINTS_AU[k])
        if cont is None:
            best = {k: min((o for o in optima[k] if o["success"]), key=lambda o: o["chi2"])
                    for k in points}  # fmt: skip
            changed, passes, n_runs = set(points), 0, 0
            while changed and passes < MAX_PASSES:
                passes += 1
                tasks = []
                for i, k in enumerate(order):
                    for j in (
                        order[i - 1] if i > 0 else None,
                        order[i + 1] if i + 1 < len(order) else None,
                    ):
                        if j is not None and j in changed:
                            s = _start_of(best[j]) | {"origin": f"cont:{j}"}
                            tasks.append({"point": k, "a_au": POINTS_AU[k], "start": s,
                                          "exact": False})  # fmt: skip
                recs = _run(tasks, pool)
                account(recs)
                n_runs += len(recs)
                changed = set()
                for r in recs:
                    k = r["point"]
                    optima[k].append(r)
                    if r["success"] and r["chi2"] < best[k]["chi2"] - CONT_TOL_CHI2:
                        best[k] = r
                        changed.add(k)
                _log(f"{dataset}: neighbour pass {passes}: {len(recs)} fits, {len(changed)} "
                     f"points improved", log_path)  # fmt: skip
            converged = not changed
            best_list = [best[k] for k in order]

            def attempt(i: int, j: int):
                s = _start_of(best_list[j]) | {"origin": f"sweep:{order[j]}"}
                r = fit_at(POINTS_AU[order[i]], s, False)
                r["point"] = order[i]
                account([r])
                optima[order[i]].append(r)
                return SimpleNamespace(chi2=r["chi2"], rec=r) if r["success"] else None

            wrapped = [SimpleNamespace(chi2=b["chi2"], rec=b) for b in best_list]
            pairs, improved = continuation_sweeps(
                wrapped, attempt, tol=CONT_TOL_CHI2, log=lambda m: _log(f"{dataset}: {m}", log_path)
            )
            cont = {"passes": passes, "converged": converged, "n_neighbour_runs": n_runs,
                    "sweep_pairs": pairs, "sweep_improved": sorted(order[i] for i in improved),
                    "best": {str(k): w.rec for k, w in zip(order, wrapped, strict=True)},
                    "added": {str(k): [o for o in optima[k] if o["origin"].split(":")[0]
                                       in ("cont", "sweep")] for k in points}}  # fmt: skip
            _write(cont_path, cont)
            _log(f"{dataset}: continuation {'converged' if converged else 'NOT converged'} after "
                 f"{passes} passes, {pairs} sweep pairs, improved {cont['sweep_improved']}",
                 log_path)  # fmt: skip
        else:
            for k in points:
                optima[k] = optima[k] + cont["added"].get(str(k), [])

        for k in points:
            p_path = _path(dataset, f"polish_{k}.json", cache)
            if p_path.exists():
                continue
            attempts = []
            pool_optima = list(optima[k])
            for attempt_no in (0, 1):
                cands = candidates(pool_optima)
                tasks = [{"point": k, "a_au": POINTS_AU[k],
                          "start": _start_of(c) | {"origin": f"polish:{c['origin']}"}}
                         for c in cands]  # fmt: skip
                recs = _run_polish(tasks, pool)
                account(recs)
                for c, r in zip(cands, recs, strict=True):
                    r["chi2_perturbative"] = c["chi2"]
                    r["from"] = c["origin"]
                passed = bool(recs) and all(r["converged"] for r in recs)
                attempts.append({"rerun": attempt_no == 1, "candidates": cands,
                                 "polished": recs, "pass": passed})  # fmt: skip
                _log(f"{dataset}: point {k} polish{' (rerun)' if attempt_no else ''}: "
                     f"{len(recs)} candidates, rounds {[len(r['rounds']) for r in recs]}, "
                     f"statuses {[r['status'] for r in recs]}, exact Q "
                     f"{[round(r['chi2'], 3) for r in recs]}, marginal "
                     f"{[round(r['marginal_m2lnl'], 3) for r in recs]}", log_path)  # fmt: skip
                if passed:
                    break
                extra = random_starts(d, k, 2 if smoke else N_RERUN, rerun=True)
                more = _run([{"point": k, "a_au": POINTS_AU[k], "start": s, "exact": False}
                             for s in extra], pool)  # fmt: skip
                account(more)
                pool_optima += more
                attempts[-1]["rerun_optima"] = more
            _write(p_path, {"point": k, "a_au": POINTS_AU[k], "attempts": attempts,
                            "pass": attempts[-1]["pass"]})  # fmt: skip
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    wall, cpu = time.time() - t0, ledger.cpu_now() - c0
    summary = {"dataset": dataset, "n_measurements": n_meas, "calls": calls, "wall_s": wall,
               "cpu_s": cpu, "smoke": smoke}  # fmt: skip
    if not smoke and cache == CACHE:
        yrs = ledger.sim_years_per_integration(_W["obs"].jd_tdb, A1_EPOCH_JD_TDB)
        ledger.append(f"m10/profile/{dataset}", "profile",
                      n_likelihood_calls=calls["pert"] + calls["exact"],
                      n_integrations=calls["exact"], sim_years=calls["exact"] * yrs,
                      wall_s=wall, cpu_s=cpu, n_workers=n_workers, n_perturbative=calls["pert"],
                      n_exact=calls["exact"])  # fmt: skip
    _log(f"{dataset}: done, {calls['pert']} perturbative and {calls['exact']} exact calls, "
         f"{wall / 60:.1f} min wall, {cpu / 3600:.2f} CPU-h", log_path)  # fmt: skip
    return summary


# ---------------------------------------------------------------------------- after the run


def _best_polished(pol: dict) -> dict | None:
    ok = [r for r in pol["attempts"][-1]["polished"] if r["converged"]]
    return min(ok, key=lambda r: r["chi2"]) if ok else None


def posterior_tightness(width_au: float) -> str:
    """Amendment 1: "tight" below 8 − 0.2 AU, "borderline" within ±0.2 AU of 8, else "not
    tight"."""
    if abs(width_au - TIGHT_AU) <= BORDERLINE_AU + 1e-9:  # 0.2 inclusive, float-safe
        return "borderline"
    return "tight" if width_au < TIGHT_AU else "not tight"


def s3_posterior():
    """(q, weights, ln L) of S3's posterior: its four cells combined with their prior weights."""
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult

    cells = [NestedResult.load(m6c.Job(S3_VARIANT, k).path(".pkl"))
             for k in range(len(m6c.variant_cells(S3_VARIANT)))]  # fmt: skip
    q, w, _ = PartitionedResult(tuple(cells), err_scale=ERR_INFLATION).samples_and_weights()
    logl = np.concatenate([np.asarray(c.results["logl"]) for c in cells])
    return q, w, logl


def _point_diagnostics(obs, ctx, space, b: dict, a: float, ids: list[str]) -> dict:
    """At the best exact point: Q per group, the IRLS weights of the pre-1781 measurements and
    of the three rows M9's screen listed, and the ancient measurements' share of Q."""
    from fork.model.observation import GROUPS

    phi = phi_of([b["log10_m"], a, b["k"], b["h"], b["lam0"]], [b["inc"], b["node"]])
    q = np.concatenate([b["u"], phi])
    r = obs.residuals_arcsec(*obs.predict(space.theta(q), ctx))
    g, scale, tau = np.asarray(obs.group_index), t_scales(), np.asarray(obs.offset_prior_arcsec)
    f = t_residuals(r, g, scale, tau)
    n = r.size
    per = (NU + 1.0) * np.log1p(f**2 / NU)
    group_q = {}
    for i, name in enumerate(GROUPS):
        group_q[name] = float(np.sum(per[:n][g == i]) + per[n + i])
    w = irls_weights(r, g, scale, tau)
    anc = [i for i, name in enumerate(GROUPS) if name.startswith("ancient")]
    sel = np.isin(g, anc)
    return {"group_q": group_q, "q_total": float(np.sum(per)),
            "ancient_q": float(sum(group_q[GROUPS[i]] for i in anc)),
            "ancient_weights": [{"id": ids[i], "resid_arcsec": float(r[i]), "weight": float(w[i])}
                                for i in np.flatnonzero(sel)],
            "screened_weights": [{"id": ids[i], "resid_arcsec": float(r[i]),
                                  "weight": float(w[i])}
                                 for i in range(n) if ids[i] in SCREENED_ROWS]}  # fmt: skip


def summarize_dataset(dataset: str, cache: Path = CACHE) -> dict | None:
    """One data set's profile from the cache: per point the exact Q, S3's marginal there and
    the fit; the reading on Δ ln L = −ΔQ/2; the maximum, the intervals, the reference points;
    the marginal's own curve and its spread against Q."""
    from fork.inference.priors import free_a_space
    from fork.model.forward import ForwardContext

    cont = _read(_path(dataset, "continuation.json", cache))
    if cont is None:
        return None
    obs = observations(dataset)
    ids = row_ids(dataset)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    space = free_a_space(_inc_raan_pq()[1], E_MAX, None, inclined=True)
    pts = []
    for k, a in enumerate(POINTS_AU):
        pol = _read(_path(dataset, f"polish_{k}.json", cache))
        starts = _read(_path(dataset, f"starts_{k}.json", cache))
        rec = {"k": k, "a_au": a, "reference": k >= N_GRID,
               "name": None if k < N_GRID else list(REFERENCE_AU)[k - N_GRID]}  # fmt: skip
        if pol is None or starts is None:
            rec |= {"status": "not run"}
            pts.append(rec)
            continue
        last = pol["attempts"][-1]
        b = _best_polished(pol) if pol["pass"] else None
        rec |= {"status": "pass" if pol["pass"] else "missing", "rerun": len(pol["attempts"]) > 1,
                "n_starts": len(starts["optima"]),
                "perturbative_best_q": cont["best"][str(k)]["chi2"],
                "candidates": [{"origin": c["origin"], "q": c["chi2"], "e": c["e"],
                                "galle_lon_deg": c["galle_lon_deg"]} for c in last["candidates"]],
                "polished": [{"from": r["from"], "q": r["chi2"],
                              "q_perturbative": r["chi2_perturbative"],
                              "marginal_m2lnl": r["marginal_m2lnl"], "status": r["status"],
                              "n_rounds": len(r["rounds"]), "converged": r["converged"],
                              "last5_change": r["rounds"][-1]["chi2_last5_change"],
                              "grad_norm": r["rounds"][-1]["grad_norm"],
                              "optimality": r["rounds"][-1]["optimality"], "nfev": r["nfev"],
                              "mass_earth": r["mass_earth"], "e": r["e"],
                              "inc_deg": r["inc_deg"], "galle_lon_deg": r["galle_lon_deg"],
                              "active_bounds": r["active_bounds"],
                              "u_box_half_widths": r["u_box_half_widths"],
                              "u_vs_u_ref_std": r["u_vs_u_ref_std"]}
                             for r in last["polished"]]}  # fmt: skip
        if b is not None:
            rec |= {"q": b["chi2"], "ln_l": -0.5 * b["chi2"], "marginal_m2lnl": b["marginal_m2lnl"],
                    "mass_earth": b["mass_earth"], "e": b["e"], "inc_deg": b["inc_deg"],
                    "galle_lon_deg": b["galle_lon_deg"], "active_bounds": b["active_bounds"],
                    "u_box_half_widths": b["u_box_half_widths"],
                    "u_vs_u_ref_std": b["u_vs_u_ref_std"],
                    "u_flagged": bool(b["u_vs_u_ref_std"] > U_FLAG_SIGMA),
                    "q_perturbative_minus_exact": b["chi2_perturbative"] - b["chi2"],
                    "n_rounds": len(b["rounds"]),
                    "diagnostics": _point_diagnostics(obs, ctx, space, b, a, ids)}  # fmt: skip
        pts.append(rec)
    grid = [p for p in pts if not p["reference"]]
    have = [p for p in grid if "ln_l" in p]
    ln_max = max(p["ln_l"] for p in have) if have else None
    m_min = min(p["marginal_m2lnl"] for p in have) if have else None
    for p in pts:
        p["dlnl"] = None if ("ln_l" not in p or ln_max is None) else p["ln_l"] - ln_max
        p["dlnl_marginal"] = (None if ("marginal_m2lnl" not in p or m_min is None)
                              else -0.5 * (p["marginal_m2lnl"] - m_min))  # fmt: skip
    a_g = [p["a_au"] for p in grid]
    d_g = [p["dlnl"] for p in grid]
    reading = read_profile(a_g, d_g)
    reading_marginal = read_profile(a_g, [p["dlnl_marginal"] for p in grid])
    used = [(a, d) for a, d in zip(a_g, d_g, strict=True) if d is not None]
    au, du = [a for a, _ in used], [d for _, d in used]
    intervals = {
        "dlnl_minus_0p5": crossings(au, du, -0.5) if used else [None, None],
        "dlnl_minus_2": crossings(au, du, DLNL_THRESHOLD) if used else [None, None],
    }
    diff = [p["marginal_m2lnl"] - p["q"] for p in have]
    spread = (max(diff) - min(diff)) if diff else None
    a_max = max(have, key=lambda p: p["ln_l"]) if have else None
    refs = {}
    for p in pts:
        if p["reference"] and "ln_l" in p:
            refs[p["name"]] = {"a_au": p["a_au"], "dlnl": p["dlnl"],
                               "dlnl_marginal": p["dlnl_marginal"],
                               "above_grid_maximum": p["dlnl"] > 0, "mass_earth": p["mass_earth"],
                               "galle_lon_deg": p["galle_lon_deg"], "e": p["e"]}  # fmt: skip
    return {"dataset": dataset, "label": LABEL[dataset], "n_measurements": obs.n_measurements,
            "n_observations": int(obs.obs_id.size), "points": pts, "reading": reading,
            "reading_marginal": reading_marginal,
            "marginal_minus_q": {"values": diff, "spread": spread,
                                 "exceeds_tolerance": None if spread is None
                                 else spread > MARGINAL_SPREAD_TOL,
                                 "tolerance": MARGINAL_SPREAD_TOL},
            "a_max": None if a_max is None else {"a_au": a_max["a_au"],
                                                 "mass_earth": a_max["mass_earth"]},
            "parabola": parabola_max(au, du) if len(used) >= 3 else None,
            "intervals": intervals, "reference_points": refs,
            "continuation": {k: cont[k] for k in ("passes", "converged", "n_neighbour_runs",
                                                  "sweep_pairs", "sweep_improved")},
            "failed_points": [p["a_au"] for p in pts if p.get("status") == "missing"],
            "box": _read(_path(dataset, "box.json", cache)),
            "ancient_q_by_a": {p["a_au"]: p["diagnostics"]["ancient_q"] for p in have},
            "convergence": {"n_polishes": sum(len(p["polished"]) for p in pts if "polished" in p),
                            "n_repolished": sum(1 for p in pts if "polished" in p
                                                for r in p["polished"] if r["n_rounds"] > 1),
                            "max_last5_change": max((r["last5_change"] for p in pts
                                                     if "polished" in p for r in p["polished"]),
                                                    default=None),
                            "all_converged": all(r["converged"] for p in pts if "polished" in p
                                                 for r in p["polished"])}}  # fmt: skip


def _m9_row(ds: dict | None) -> dict | None:
    if ds is None or ds["a_max"] is None:
        return None
    pt = {p["a_au"]: p for p in ds["points"] if not p["reference"]}
    refs = ds["reference_points"]
    return {"reading": ds["reading"]["reading"], "a_max_au": ds["a_max"]["a_au"],
            "mass_at_a_max": ds["a_max"]["mass_earth"],
            "dlnl_minus_2": ds["intervals"]["dlnl_minus_2"],
            "dlnl_minus_0p5": ds["intervals"]["dlnl_minus_0p5"],
            "dlnl_neptune": refs["neptune"]["dlnl"], "dlnl_le_verrier": refs["le_verrier"]["dlnl"],
            "dlnl_26": pt[26.0]["dlnl"]}  # fmt: skip


def comparison_table(data: dict, m9_out: dict) -> dict:
    """The M10.md comparison table: A and B under both models (M9's values from the frozen
    M9_profile.json), and the four leave-one-out sets under Student-t."""
    return {"gaussian_m9": {ds: _m9_row(m9_out["datasets"].get(ds)) for ds in ("all", "post1781")},
            "student_t_m10": {ds: _m9_row(data.get(ds)) for ds in DATASETS}}  # fmt: skip


INTERPRETATION = {
    "i": "A constrained under both models: the distance result is robust to tails.",
    "ii": "A flat under Student-t (constrained under Gaussian): the distance result depends on "
    "Gaussian tails, that is, on a few old measurements; under heavier tails the data allow any "
    "a from 26 to 44 AU.",
    "iii": "A leave-one-out set is flat (A constrained under Student-t): that observer carries "
    "the distance.",
    "iv": "Every leave-one-out set is constrained (A constrained under Student-t): no single "
    "pre-1781 observer carries the distance on his own.",
}


def interpret(data: dict, m9_out: dict) -> dict:
    """M10.md's interpretation rows (i)–(iv), read after the three-way reading of every set;
    (iii)/(iv) only if A is constrained under Student-t (amendment 2 keeps every leave-one-out
    profile reported as description regardless)."""
    r = {ds: (data[ds]["reading"]["reading"] if data.get(ds) else None) for ds in DATASETS}
    a_gauss = m9_out["datasets"]["all"]["reading"]["reading"]
    rows: list[dict] = []
    loo = [ds for ds in DATASETS if ds in OBSERVER]
    if r["all"] == "constrained" and a_gauss == "constrained":
        rows.append({"row": "i", "text": INTERPRETATION["i"]})
    elif r["all"] == "flat" and a_gauss == "constrained":
        rows.append({"row": "ii", "text": INTERPRETATION["ii"]})
    elif r["all"] in (None, "partial"):
        rows.append({"row": None, "text": "A is partial or missing under Student-t: no row is "
                     "used; the curves and the failed clause are reported."})  # fmt: skip
    if r["all"] == "constrained":
        flat = [ds for ds in loo if r[ds] == "flat"]
        partial = [ds for ds in loo if r[ds] in (None, "partial")]
        if flat:
            rows.append({"row": "iii", "text": INTERPRETATION["iii"], "observers":
                         [OBSERVER[ds] for ds in flat]})  # fmt: skip
        elif all(r[ds] == "constrained" for ds in loo):
            rows.append({"row": "iv", "text": INTERPRETATION["iv"]})
        if partial:
            rows.append(
                {
                    "row": None,
                    "text": "Leave-one-out sets reported as partial, named in "
                    "neither (iii) nor (iv): " + ", ".join(OBSERVER[ds] for ds in partial),
                }
            )
    else:
        rows.append({"row": None, "text": "A is not constrained under Student-t: (iii)/(iv) are "
                     "not read; the leave-one-out sets are reported as curves only."})  # fmt: skip
    ab = (r["all"], r["post1781"])
    return {"readings": r, "a_gaussian": a_gauss, "rows": rows,
            "a_x_b_student_t": {"a": ab[0], "b": ab[1],
                                "text": m9.AB_TABLE.get(ab, "A profile is partial: the table is "
                                                        "not used.")}}  # fmt: skip


SITE_OUT_NAME = "m10_profile.json"
ERROR_MODEL = "Student-t errors (ν = 4) with the scales of the Student-t sensitivity run"


def _site_points(ds: dict) -> list[dict]:
    keys = ("a_au", "dlnl", "dlnl_marginal", "mass_earth", "e", "inc_deg", "galle_lon_deg",
            "active_bounds", "reference", "name", "status", "u_flagged")  # fmt: skip
    return [{k: p.get(k) for k in keys} for p in ds["points"]]


def site_payload(out: dict) -> dict:
    """The compact site copy of the report: what the page's § 5 text and Figs. 7 and 9 read."""
    ds = {}
    for name, d in out["datasets"].items():
        if d is None:
            continue
        ds[name] = {"label": d["label"], "observer": OBSERVER.get(name),
                    "measurements_removed": (out["datasets"]["all"]["n_measurements"]
                                             - d["n_measurements"]),
                    "n_measurements": d["n_measurements"],
                    "n_observations": d["n_observations"], "points": _site_points(d),
                    "reading": d["reading"], "reading_marginal": d["reading_marginal"]["reading"],
                    "marginal_minus_q_spread": d["marginal_minus_q"]["spread"],
                    "a_max": d["a_max"], "parabola": d["parabola"], "intervals": d["intervals"],
                    "reference_points": d["reference_points"],
                    "failed_points": d["failed_points"],
                    "convergence": d["convergence"]}  # fmt: skip
    marg = out["s3_posterior_marginal_a"]
    return {"description": "M10, the pre-registered profile likelihood in the semi-major axis "
            "under Student-t errors (docs/notes/M10.md): the site's copy of "
            "docs/notes/M10_profile.json. Q is on the χ² scale; Δ ln L = −ΔQ/2.",
            "pre_registration": "docs/notes/M10.md", "code_commit": out["code_commit"],
            "settings": out["settings"], "error_model": ERROR_MODEL, "datasets": ds,
            "leave_one_out": [ds_ for ds_ in DATASETS if ds_ in OBSERVER],
            "s3_posterior_marginal_a": {k: marg[k] for k in (
                "a_au", "width_95_au", "tightness", "tight_below_au", "borderline_au",
                "published_a_au")},
            "profile_vs_s3_posterior": out.get("profile_vs_s3_posterior"),
            "comparison_with_m9": out["comparison_with_m9"],
            "interpretation": out.get("interpretation"),
            "neptune_osculating_a_1781_au": out["neptune_osculating_a_1781_au"],
            "ledger": out["ledger"]}  # fmt: skip


def report(cache: Path = CACHE, out_json: Path = OUT_JSON, out_png: Path | None = OUT_PNG,
           site: bool = True) -> dict:  # fmt: skip
    """docs/notes/M10_profile.json (and the figure) from the cached runs of the six data sets,
    and the site's copy ``site/data/m10_profile.json`` once every data set is in."""
    from fork.export.to_json import write_json
    from fork.inference.priors import PRIOR_LOG10_M_EARTH
    from fork.model.observation import GROUPS

    check_frozen()
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()  # fmt: skip
    m9_out = json.loads(M9_JSON.read_text())
    data = {ds: summarize_dataset(ds, cache) for ds in DATASETS}
    q, w, _ = s3_posterior()
    marg = posterior_marginal_a(q, w)
    width = marg["a_au"]["hi95"] - marg["a_au"]["lo95"]
    marg["width_95_au"] = width
    marg["tightness"] = posterior_tightness(width)
    marg["tight"] = marg["tightness"] == "tight"
    marg["borderline_au"] = BORDERLINE_AU
    marg["published_a_au"] = json.loads((ROOT / "site" / "data" / "m6_fits.json").read_text())[
        "fits"][S3_VARIANT]["a_au"]  # fmt: skip
    out: dict = {"description": "M10: profile likelihood in the semi-major axis under Student-t "
                 "errors (ν = 4, S3's scales), pre-registered in docs/notes/M10.md (2026-11-10). "
                 "Diagnostic only; no published number changes. Q is on the χ² scale; "
                 "Δ ln L = −ΔQ/2.",
                 "code_commit": head, "frozen": FROZEN,
                 "settings": {"nu": NU, "t_scales_arcsec": dict(zip(GROUPS, t_scales().tolist(),
                                                                    strict=True)),
                              "seed_namespace": SEED_NAMESPACE, "datasets": DATASETS,
                              "grid_au": A_GRID_AU, "reference_au": REFERENCE_AU,
                              "n_random": N_RANDOM, "n_posterior": N_POSTERIOR,
                              "n_nearest": N_NEAREST, "n_rerun": N_RERUN, "e_max": E_MAX,
                              "log10_m_range": list(PRIOR_LOG10_M_EARTH),
                              "inc_max_deg": INC_MAX_DEG, "continuation_tol_q": CONT_TOL_CHI2,
                              "polish_window_q": POLISH_WINDOW_CHI2, "max_polish": MAX_POLISH,
                              "max_rounds": MAX_ROUNDS,
                              "distinct": {"galle_lon_deg": DISTINCT_LON_DEG, "e": DISTINCT_E},
                              "reading": {"dlnl": DLNL_THRESHOLD, "x_au": X_AU,
                                          "flat_range_au": FLAT_RANGE_AU},
                              "marginal_spread_tolerance": MARGINAL_SPREAD_TOL,
                              "borderline_au": BORDERLINE_AU, "status_codes": STATUS_CODES},
                 "datasets": data, "s3_posterior_marginal_a": marg,
                 "comparison_with_m9": comparison_table(data, m9_out),
                 "neptune_osculating_a_1781_au":
                     m9_out["neptune_osculating_a_1781_au"]}  # fmt: skip
    a_all = data["all"]
    if a_all is not None:
        vol = []
        for p in a_all["points"]:
            if p["reference"]:
                continue
            i = A_GRID_AU.index(p["a_au"])
            ok = (p["dlnl"] is not None and marg["ln_p_rel"][i] is not None
                  and marg["n_eff"][i] >= MIN_EFF_SAMPLES)  # fmt: skip
            vol.append({"a_au": p["a_au"], "volume_term": (marg["ln_p_rel"][i] - p["dlnl"])
                        if ok else None, "n_eff": marg["n_eff"][i]})  # fmt: skip
        out["volume_term"] = vol
        cmp = m9.compare_with_posterior(a_all["reading"], a_all["intervals"], marg)
        cmp["posterior_tightness"] = marg["tightness"]
        if marg["tightness"] == "borderline" and a_all["reading"]["reading"] == "flat":
            cmp["reading"] = (
                "flat profile; the posterior's 95 % width is borderline (within "
                f"{BORDERLINE_AU} AU of {TIGHT_AU} AU): no tight/loose claim"
            )
        out["profile_vs_s3_posterior"] = cmp
    if all(data.get(ds) for ds in ("all", "post1781")):
        out["interpretation"] = interpret(data, m9_out)
    ledger_rows = [json.loads(line) for line in (ROOT / "data" / "cache" / "compute_ledger.jsonl")
                   .read_text().splitlines() if '"m10/' in line]  # fmt: skip
    out["ledger"] = [{k: r.get(k) for k in ("job", "n_likelihood_calls", "n_integrations",
                                            "wall_s", "cpu_s", "n_workers")}
                     for r in ledger_rows]  # fmt: skip
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=1, default=float))
    if site and all(data.values()):
        write_json(SITE_OUT_NAME, site_payload(out))
    if out_png is not None:
        from fork.export.m10_plot import draw

        draw(out, m9_out, out_png)
    print(f"wrote {out_json}")
    return out


# ---------------------------------------------------------------------------- CLI


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profile")
    p.add_argument("--dataset", choices=DATASETS, required=True)
    p.add_argument("--n-workers", type=int, default=12)
    p.add_argument("--cache", type=Path, default=CACHE)
    p.add_argument("--smoke", action="store_true", help="2 points, 3 starts; scratch --cache only")
    r = sub.add_parser("report")
    r.add_argument("--cache", type=Path, default=CACHE)
    r.add_argument("--no-png", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "profile":
        cache = a.cache.resolve()
        if a.smoke and cache == CACHE.resolve():
            raise SystemExit("smoke runs go to a scratch --cache, never the run cache")
        profile(a.dataset, a.n_workers, cache, smoke=a.smoke)
    else:
        report(a.cache.resolve(), out_png=None if a.no_png else OUT_PNG)


if __name__ == "__main__":
    main()
