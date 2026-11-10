"""M9: profile likelihood in the semi-major axis (docs/notes/M9.md, pre-registered 2026-11-09).

    uv run python -m fork.export.m9 profile --dataset {all,post1781} [--n-workers 12]
    uv run python -m fork.export.m9 report      # docs/notes/M9_profile.json and .png from the cache
    uv run python -m fork.export.m9 screen      # the transcription screen (report only)

A diagnostic: it changes no published result. Every setting below is M9.md's; nothing is retuned
after a result is seen.

- **Grid:** ``A_GRID_AU`` (13 points) plus two reference points, Neptune's a and Le Verrier's,
  profiled the same way but outside the grid and the reading rule.
- **Datasets:** ``all`` (Table 77's 529 measurements) and ``post1781`` (the ``ancient`` group's
  37 measurements removed: the observations before Herschel's 1781 discovery). The headline's
  frozen error model in both; the ancient offsets drop out of the augmented χ² with their rows.
- **Profiled:** ln L(a) = −½ χ²_min(a), the headline's augmented χ² (``RaDecObservations.whitened``)
  minimised over the perturber's other six parameters (bounded to the headline's prior support:
  1–100 M⊕, e ≤ 0.3, i ≤ 10°), Uranus's four in-plane elements (free, docs/RULES.md rule 2) and
  the eight group offsets (closed form).
- **Search:** the M8 perturbative likelihood (``PerturbativeProblem``, first order, Δs = 60 d,
  cubic, u_ref from ``M8_u_ref.json``), in which Uranus and the offsets are already solved; the
  residual handed to ``scipy.optimize.least_squares`` (trust-region reflective) is the whitened
  vector projected off J_u's column space, whose squared norm is ``GaussianMarginal.chi2_min``.
  The e ≤ 0.3 disk is ``fork.inference.fit._fit_bounded``'s two-phase guard; log₁₀ m and i are
  box bounds of the optimizer. 32 random prior draws + 8 replacement-posterior starts per point,
  then neighbour continuation passes and sweeps (``fork.inference.profile``'s conventions).
- **Confirmation:** the best perturbative optimum, and up to two more distinct ones within 4 in
  χ², are polished with the exact likelihood (``fork.model.forward``, Wisdom–Holman) over all ten
  parameters at the held a. The lowest exact χ² is the profile value.
- **Checks:** every polish ends with scipy status > 0 and the continuation converged; a failed
  point is rerun once with 40 more random starts, then reported missing.
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

from fork.export import m6c
from fork.export.to_json import A1_EPOCH_JD_TDB

ROOT = m6c.ROOT
CACHE = ROOT / "data" / "cache" / "m9"
PROGRESS_LOG = CACHE / "progress.log"
NOTES = ROOT / "docs" / "notes"
OUT_JSON = NOTES / "M9_profile.json"
OUT_PNG = NOTES / "M9_profile.png"
SCREEN_JSON = NOTES / "M9_screen.json"
U_REF_FILE = NOTES / "M8_u_ref.json"
CONJ_JSON = ROOT / "site" / "data" / "m6c_conjunction.json"
REPLACEMENT_JSON = ROOT / "site" / "data" / "m6c_replacement.json"
SPOTCHECK_ROWS = ROOT / "docs" / "review" / "spotcheck_rows.csv"

FROZEN = {  # M9.md, "Frozen inputs (SHA-256)"
    "data/historical/leverrier_1846.csv":
        "2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4",
    "data/historical/m6_error_model.json":
        "0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277",
    "docs/notes/M8_u_ref.json":
        "4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5",
}  # fmt: skip

A_GRID_AU = (24.0, 26.0, 28.0, 29.0, 29.5, 30.0, 30.5, 31.0, 32.0, 34.0, 36.0, 40.0, 44.0)
REFERENCE_AU = {"neptune": 30.07, "le_verrier": 36.15}  # k = 13, 14; not in the reading rule
POINTS_AU = A_GRID_AU + tuple(REFERENCE_AU.values())
N_GRID = len(A_GRID_AU)
DATASETS = ("all", "post1781")  # d = 0, 1 in the start rng
N_ANCIENT_MEASUREMENTS = 37
E_MAX = 0.3  # the headline's prior
INC_MAX_DEG = 10.0
N_RANDOM = 32
N_POSTERIOR = 8
N_NEAREST = 200
N_RERUN = 40
CONT_TOL_CHI2 = 0.01
MAX_PASSES = 60
POLISH_WINDOW_CHI2 = 4.0
MAX_POLISH = 3
DISTINCT_LON_DEG = 2.0
DISTINCT_E = 0.03
U_FLAG_SIGMA = 3.0  # B: Uranus farther than this (headline posterior std) from u_ref is flagged
DLNL_THRESHOLD = -2.0
X_AU = 4.0
FLAT_RANGE_AU = (26.0, 44.0)
TIGHT_AU = 8.0  # 2X: the posterior's 95 % interval in a is "tight" below this width
MIN_EFF_SAMPLES = 50
SCREEN_SIGMA = 4.0
PENALTY = 1e3
MAX_NFEV = 3000

# q of the uncut inclined free-a space: u (4), log10 m, ln a, k, h, λ₀, sin i cos Ω, sin i sin Ω
Q_LN_A, Q_LOG10_M = 5, 4


# ---------------------------------------------------------------------------- problem


def check_frozen() -> None:
    from fork.export import m8_runs

    m8_runs.check_frozen("headline")  # the frozen data, the pinned M8 code and u_ref
    for rel, digest in FROZEN.items():
        if m6c.sha256(ROOT / rel) != digest:
            raise RuntimeError(f"{rel} differs from M9's frozen input")


def _log(msg: str, path: Path = PROGRESS_LOG) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with path.open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def drop_measurements(obs, keep: np.ndarray):
    """``obs`` without the measurements where ``keep`` is False; epochs with nothing left go
    too. Groups keep their σ and τ: an empty group's best offset is 0 and its prior term 0."""
    from fork.model.apparent import ApparentGeometry
    from fork.model.observation import RaDecObservations

    keep = np.asarray(keep, dtype=bool)
    ep = np.unique(obs.epoch_index[keep])
    new = np.full(obs.obs_id.size, -1)
    new[ep] = np.arange(ep.size)
    g = obs.geometry
    geom = ApparentGeometry(
        g.jd_tt[ep],
        g.earth_pos_ecl_j2000_au[ep],
        g.earth_vel_ecl_j2000_au_d[ep],
        g.ecl_j2000_to_true_of_date[ep],
    )
    return RaDecObservations(geom, obs.obs_id[ep], new[obs.epoch_index[keep]], obs.is_dec[keep],
                             obs.value_rad[keep], obs.group_index[keep], obs.sigma_arcsec,
                             obs.offset_prior_arcsec)  # fmt: skip


def observations(dataset: str):
    """Table 77 with the headline's frozen error model, exactly as ``m6c.dataset`` builds it
    (``all``), or without the ancient group's measurements (``post1781``)."""
    from fork.data.calibration import real_observations
    from fork.model import observation as om

    if dataset not in DATASETS:
        raise ValueError(f"dataset must be one of {DATASETS}")
    em = om.ErrorModel.load(ROOT / m6c.VARIANT_CONFIG["headline"]["error_model"])
    obs = real_observations(om.read_rows(), em)
    if dataset == "all":
        return obs
    ancient = [i for i, g in enumerate(om.GROUPS) if g.startswith("ancient")]
    return drop_measurements(obs, ~np.isin(obs.group_index, ancient))


def _inc_raan_pq() -> tuple[tuple[float, float], tuple[float, float]]:
    """Uranus's held (i, Ω) at the epoch and the (p, q) pair, as ``m6c.dataset`` sets them."""
    from fork.inference.fit import uranus_theta_to_x
    from fork.model.forward import uranus_inc_raan_at_epoch

    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    return inc_raan, tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])


def problem(dataset: str):
    """(obs, ctx, uncut space, PerturbativeProblem, u_ref, u_std) for a data set. The space is
    the headline's (e ≤ 0.3, inclined, a free) on the data set's own Uranus box; u_ref and its
    posterior std are the headline's (``M8_u_ref.json``) for both."""
    from fork.export import m8_runs
    from fork.inference.fit_radec import uranus_box_radec
    from fork.inference.priors import free_a_space
    from fork.model.forward import ForwardContext
    from fork.model.perturbative import PerturbativeProblem

    check_frozen()
    obs = observations(dataset)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan, pq = _inc_raan_pq()
    box = uranus_box_radec(obs, ctx, inc_raan)
    space = free_a_space(pq, E_MAX, box, inclined=True)
    u = json.loads(U_REF_FILE.read_text())["headline"]
    u_ref, u_std = np.array(u["u_ref"], dtype=float), np.array(u["u_std"], dtype=float)
    pr = PerturbativeProblem.build(obs, ctx, space, u_ref, step_days=m8_runs.STEP_DAYS,
                                   quadrature=m8_runs.QUADRATURE)  # fmt: skip
    return obs, ctx, space, pr, u_ref, u_std


def phi_of(y: np.ndarray, inc: np.ndarray) -> np.ndarray:
    """(log10 m, a, k, h, λ₀) and (i, Ω) → the 7 perturber q entries."""
    si = math.sin(inc[0])
    return np.array([y[0], math.log(y[1]), y[2], y[3], y[4], si * math.cos(inc[1]),
                     si * math.sin(inc[1])])  # fmt: skip


# ---------------------------------------------------------------------------- starts


def random_starts(d: int, k: int, n: int = N_RANDOM, rerun: bool = False) -> list[dict]:
    """``n`` draws from the headline prior (at the held a) with rng [2026, 9, d, k] (plus 1 for
    the rerun): log10 m uniform on [0, 2]; (√e cos ϖ, √e sin ϖ) uniform on the disk e ≤ 0.3;
    λ₀ uniform; cos i uniform on [cos 10°, 1] and Ω uniform."""
    from fork.inference.priors import PRIOR_LOG10_M_EARTH

    rng = np.random.default_rng([2026, 9, d, k] + ([1] if rerun else []))
    out = []
    for i in range(n):
        log10_m = rng.uniform(*PRIOR_LOG10_M_EARTH)
        e, varpi = E_MAX * rng.random(), 2.0 * math.pi * rng.random()
        lam0 = 2.0 * math.pi * rng.random()
        cos_i = 1.0 - (1.0 - math.cos(math.radians(INC_MAX_DEG))) * rng.random()
        node = 2.0 * math.pi * rng.random()
        out.append({"log10_m": float(log10_m), "k": math.sqrt(e) * math.cos(varpi),
                    "h": math.sqrt(e) * math.sin(varpi), "lam0": float(lam0),
                    "inc": math.acos(cos_i), "node": float(node),
                    "origin": f"{'rerun' if rerun else 'random'}:{i}"})  # fmt: skip
    return out


def replacement_posterior():
    """(q, weights, ln L) of the M6c replacement posterior, every cell's weighted samples."""
    from fork.export import m6c_repeat
    from fork.inference.nested import ERR_INFLATION, PartitionedResult

    cells = m6c_repeat.replacement_cells()
    q, w, _ = PartitionedResult(tuple(cells), err_scale=ERR_INFLATION).samples_and_weights()
    logl = np.concatenate([np.asarray(c.results["logl"]) for c in cells])
    return q, w, logl


def posterior_starts(q: np.ndarray, logl: np.ndarray, a_au: float,
                     n: int = N_POSTERIOR, n_nearest: int = N_NEAREST) -> list[dict]:  # fmt: skip
    """Of the ``n_nearest`` posterior samples nearest the grid point in ln a, the ``n`` with the
    highest likelihood, with a set to the grid value. Ties in ln a are broken by sample order."""
    near = np.argsort(np.abs(q[:, Q_LN_A] - math.log(a_au)), kind="stable")[:n_nearest]
    best = near[np.argsort(-logl[near], kind="stable")[:n]]
    out = []
    for j, i in enumerate(best):
        s = q[i]
        out.append({"log10_m": float(s[4]), "k": float(s[6]), "h": float(s[7]),
                    "lam0": float(s[8]), "inc": math.asin(min(math.hypot(s[9], s[10]), 1.0)),
                    "node": math.atan2(s[10], s[9]) % (2.0 * math.pi), "u": s[:4].tolist(),
                    "origin": f"posterior:{j}", "sample_ln_l": float(logl[i])})  # fmt: skip
    return out


# ---------------------------------------------------------------------------- fits (workers)

_W: dict = {}


def _init_worker(dataset: str) -> None:
    from fork.inference.profile import GalleGeometry

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    obs, ctx, space, pr, u_ref, u_std = problem(dataset)
    q_u, _ = np.linalg.qr(pr.jac_u)  # orthonormal basis of J_u's columns
    geo = GalleGeometry.from_horizons()
    conj = json.loads(CONJ_JSON.read_text())["galle_night_jd_tdb"]
    if abs(geo.jd_tdb - conj) > 1e-9:
        raise RuntimeError("Galle's night differs between the geometry and m6c_conjunction.json")
    _W.update(dataset=dataset, obs=obs, ctx=ctx, space=space, pr=pr, u_ref=u_ref, u_std=u_std,
              q_u=q_u, geo=geo, calls={"pert": 0, "exact": 0})  # fmt: skip


def _pert_residual(phi: np.ndarray) -> np.ndarray:
    """The whitened vector at u_ref projected off J_u's columns: ‖·‖² = χ² minimised over
    Uranus (``GaussianMarginal.chi2_min``) and, inside ``whitened``, over the offsets."""
    _W["calls"]["pert"] += 1
    pr = _W["pr"]
    try:
        z = pr.z_from_delta(pr.delta_x(phi, second_order=False))
    except (ValueError, ZeroDivisionError, FloatingPointError, OverflowError):
        return np.full(pr.jac_u.shape[0], PENALTY)
    if not np.all(np.isfinite(z)):
        return np.full(pr.jac_u.shape[0], PENALTY)
    q_u = _W["q_u"]
    return z - q_u @ (q_u.T @ z)


def _exact_residual(q: np.ndarray) -> np.ndarray:
    _W["calls"]["exact"] += 1
    obs, space = _W["obs"], _W["space"]
    n = obs.n_measurements + obs.sigma_arcsec.size
    if q[0] <= 0.0 or math.hypot(q[1], q[2]) >= 0.95:
        return np.full(n, PENALTY)
    try:
        z = obs.whitened(space.theta(q), _W["ctx"])
    except (ValueError, ZeroDivisionError, FloatingPointError, OverflowError):
        return np.full(n, PENALTY)
    return z if np.all(np.isfinite(z)) else np.full(n, PENALTY)


def fit_at(a_au: float, start: dict, exact: bool) -> dict:
    """One bounded least-squares fit at the held a from ``start`` (log10_m, k, h, lam0, inc,
    node; u for ``exact``). ``_fit_bounded`` runs the interior and, if needed, boundary layouts
    of the e disk; log10 m and i sit in the optimizer's box bounds; (i, Ω) ride along as
    variables ``_fit_bounded`` does not see, warm-started from the previous run."""
    from fork.inference.fit import _fit_bounded, _prior_box, le_verrier_space

    sp = le_verrier_space(a_au, E_MAX)
    inc_max = math.radians(INC_MAX_DEG)
    inc0 = np.array([start["inc"], start["node"]], dtype=float)
    nh = 4 if exact else 0
    head0 = np.asarray(start["u"], dtype=float) if exact else np.zeros(4)
    c0 = dict(_W["calls"])

    def lm(layout, z, warm):
        z = np.asarray(z, dtype=float)
        lo_p, hi_p = _prior_box(layout)
        inc_s = warm.inc if warm is not None else inc0
        x0 = np.concatenate([z[:4] if exact else [], inc_s, z[4:]])
        lo = np.concatenate([np.full(nh, -np.inf), [0.0, -np.inf], lo_p])
        hi = np.concatenate([np.full(nh, np.inf), [inc_max, np.inf], hi_p])
        x0 = np.clip(x0, lo, hi)

        def fun(x):
            phi = phi_of(layout.y_from_z(x[nh + 2 :]), x[nh : nh + 2])
            return _exact_residual(np.concatenate([x[:4], phi])) if exact else _pert_residual(phi)

        sol = least_squares(fun, x0, method="trf", bounds=(lo, hi), x_scale="jac", xtol=1e-10,
                            ftol=1e-10, gtol=1e-10, max_nfev=MAX_NFEV)  # fmt: skip
        return SimpleNamespace(x=np.concatenate([sol.x[:4] if exact else z[:4], sol.x[nh + 2:]]),
                               cost=float(sol.cost), success=bool(sol.success),
                               status=int(sol.status), nfev=int(sol.nfev),
                               inc=sol.x[nh:nh + 2].copy(), message=str(sol.message))  # fmt: skip

    z0 = np.concatenate([head0, sp.z_from_y([start["log10_m"], a_au, start["k"], start["h"],
                                             start["lam0"]])])  # fmt: skip
    best, layout = _fit_bounded(lm, z0, sp)
    rec = _optimum(a_au, layout.y_from_z(best.x[4:]), best.inc, best, layout.on_e_bound, exact,
                   best.x[:4] if exact else None)  # fmt: skip
    rec["origin"] = start.get("origin")
    rec["calls"] = {k: _W["calls"][k] - c0[k] for k in c0}
    return rec


def _optimum(a_au, y, inc, sol, on_e_bound: bool, exact: bool, u) -> dict:
    from fork.inference.priors import PRIOR_LOG10_M_EARTH
    from fork.inference.profile import galle_lon_of_date_rad

    phi = phi_of(y, inc)
    if u is None:  # the perturbative search: Uranus from the linear marginal at this φ
        u = _W["pr"].evaluate(phi).marginal.u_star
    u = np.asarray(u, dtype=float)
    space, box = _W["space"], _W["space"].uranus_box
    theta = space.theta(np.concatenate([u, phi]))
    lon = math.degrees(galle_lon_of_date_rad(theta, _W["ctx"], _W["geo"]))
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
    rec = fit_at(t["a_au"], t["start"], t["exact"])
    rec["point"] = t["point"]
    return rec


def _start_of(rec: dict) -> dict:
    """A fit record as a start for another fit (same or a neighbouring a)."""
    return {k: rec[k] for k in ("log10_m", "k", "h", "lam0", "inc", "node", "u")}


def _run(tasks: list[dict], pool) -> list[dict]:
    if pool is None:
        return [_task(t) for t in tasks]
    return list(pool.map(_task, tasks, chunksize=1))


def _distinct(a: dict, b: dict) -> bool:
    dlon = abs((a["galle_lon_deg"] - b["galle_lon_deg"] + 180.0) % 360.0 - 180.0)
    return dlon > DISTINCT_LON_DEG or abs(a["e"] - b["e"]) > DISTINCT_E


def candidates(optima: list[dict]) -> list[dict]:
    """The best successful perturbative optimum and the distinct others within
    ``POLISH_WINDOW_CHI2`` of it, best first, at most ``MAX_POLISH``."""
    ok = sorted((o for o in optima if o["success"]), key=lambda o: o["chi2"])
    out: list[dict] = []
    for o in ok:
        if out and o["chi2"] > out[0]["chi2"] + POLISH_WINDOW_CHI2:
            break
        if all(_distinct(o, p) for p in out):
            out.append(o)
        if len(out) == MAX_POLISH:
            break
    return out


# ---------------------------------------------------------------------------- the profile run


def _path(dataset: str, name: str, cache: Path) -> Path:
    return cache / dataset / name


def _write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload, indent=1, default=float))
    tmp.replace(path)


def _read(path: Path) -> dict | None:
    return json.loads(path.read_text()) if path.exists() else None


def profile(dataset: str, n_workers: int, cache: Path = CACHE, *, smoke: bool = False) -> dict:
    """The whole M9 profile for one data set, resumable: the independent starts per point, the
    continuation, the exact polishes and the rerun are each cached as they finish."""
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
    _log(f"{dataset}: {n_meas} measurements, {_W['obs'].obs_id.size} observations, "
         f"{n_workers} workers, smoke {smoke}", log_path)  # fmt: skip
    pool = (get_context("spawn").Pool(n_workers, initializer=_init_worker, initargs=(dataset,))
            if n_workers > 1 else None)  # fmt: skip
    try:
        # 1. independent starts, one point at a time, each written as it finishes (a run that
        #    dies loses at most one point's 40 fits)
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
                   {"point": k, "a_au": POINTS_AU[k], "n_random": N_RANDOM,
                    "n_posterior": N_POSTERIOR, "optima": recs})  # fmt: skip
            ok = [r["chi2"] for r in recs if r["success"]]
            _log(f"{dataset}: point {k} (a = {POINTS_AU[k]} AU) starts done, {len(recs)} fits, "
                 f"best χ² {min(ok):.3f}" if ok else
                 f"{dataset}: point {k} (a = {POINTS_AU[k]} AU) starts done, none succeeded",
                 log_path)  # fmt: skip
        optima = {k: _read(_path(dataset, f"starts_{k}.json", cache))["optima"] for k in points}

        # 2. continuation: parallel neighbour passes, then sequential sweeps (profile.py's)
        cont_path = _path(dataset, "continuation.json", cache)
        cont = _read(cont_path)
        order = sorted(points, key=lambda k: POINTS_AU[k])  # the reference points join in order
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

        # 3. exact polishes of the candidates, per point; 4. the rerun of a failed point
        for k in points:
            p_path = _path(dataset, f"polish_{k}.json", cache)
            if p_path.exists():
                continue
            attempts = []
            pool_optima = list(optima[k])
            for attempt_no in (0, 1):
                cands = candidates(pool_optima)
                tasks = [{"point": k, "a_au": POINTS_AU[k], "exact": True,
                          "start": _start_of(c) | {"origin": f"polish:{c['origin']}"}}
                         for c in cands]  # fmt: skip
                recs = _run(tasks, pool)
                account(recs)
                for c, r in zip(cands, recs, strict=True):
                    r["chi2_perturbative"] = c["chi2"]
                    r["from"] = c["origin"]
                passed = bool(recs) and all(r["status"] > 0 for r in recs)
                attempts.append({"rerun": attempt_no == 1, "candidates": cands,
                                 "polished": recs, "pass": passed})  # fmt: skip
                _log(f"{dataset}: point {k} polish{' (rerun)' if attempt_no else ''}: "
                     f"{len(recs)} candidates, statuses {[r['status'] for r in recs]}, exact χ² "
                     f"{[round(r['chi2'], 3) for r in recs]}", log_path)  # fmt: skip
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
        ledger.append(f"m9/profile/{dataset}", "profile",
                      n_likelihood_calls=calls["pert"] + calls["exact"],
                      n_integrations=calls["exact"], sim_years=calls["exact"] * yrs,
                      wall_s=wall, cpu_s=cpu, n_workers=n_workers, n_perturbative=calls["pert"],
                      n_exact=calls["exact"])  # fmt: skip
    _log(f"{dataset}: done, {calls['pert']} perturbative and {calls['exact']} exact calls, "
         f"{wall / 60:.1f} min wall, {cpu / 3600:.2f} CPU-h", log_path)  # fmt: skip
    return summary


# ---------------------------------------------------------------------------- reading


def read_profile(a_au: list[float], dlnl: list[float | None]) -> dict:
    """The pre-registered reading of one profile over the grid points (``None`` = missing, not
    used): constrained, flat or partial, with every clause reported."""
    pts = [(a, d) for a, d in zip(a_au, dlnl, strict=True) if d is not None]
    if not pts:
        return {"reading": "partial", "a_max_au": None, "clauses": {}, "failed": ["no points"]}
    top = max(d for _, d in pts)
    pts = [(a, d - top) for a, d in pts]  # Δ ln L relative to the grid maximum
    a_max = max(pts, key=lambda p: p[1])[0]
    t = DLNL_THRESHOLD
    clauses = {
        "a_max_not_at_grid_end": a_max not in (A_GRID_AU[0], A_GRID_AU[-1]),
        "falls_below_within_x_left": any(a_max - X_AU <= a < a_max and d < t for a, d in pts),
        "falls_below_within_x_right": any(a_max < a <= a_max + X_AU and d < t for a, d in pts),
        "no_second_branch": all(d < t for a, d in pts if abs(a - a_max) > X_AU),
        "flat_26_44": all(d > t for a, d in pts if FLAT_RANGE_AU[0] <= a <= FLAT_RANGE_AU[1]),
    }
    constrained_keys = ("a_max_not_at_grid_end", "falls_below_within_x_left",
                        "falls_below_within_x_right", "no_second_branch")  # fmt: skip
    if all(clauses[k] for k in constrained_keys):
        reading, failed = "constrained", []
    elif clauses["flat_26_44"]:
        reading, failed = "flat", []
    else:
        reading = "partial"
        failed = [k for k in constrained_keys if not clauses[k]] + ["flat_26_44"]
    where = {
        "second_branch_points": [a for a, d in pts if abs(a - a_max) > X_AU and d >= t],
        "below_threshold_in_26_44": [
            a for a, d in pts if FLAT_RANGE_AU[0] <= a <= FLAT_RANGE_AU[1] and d <= t
        ],
    }
    thresholds = {"dlnl": t, "x_au": X_AU, "flat_range_au": FLAT_RANGE_AU}
    missing = [a for a, d in zip(a_au, dlnl, strict=True) if d is None]
    return {"reading": reading, "a_max_au": a_max, "clauses": clauses, "failed": failed,
            "where": where, "thresholds": thresholds, "missing_points": missing}  # fmt: skip


def crossings(a_au: list[float], dlnl: list[float], level: float) -> list[float | None]:
    """[lo, hi] where the profile crosses ``level`` on each side of its maximum, by linear
    interpolation between grid points; None where it stays above the level to the grid's end."""
    a, d = np.asarray(a_au, float), np.asarray(dlnl, float)
    i = int(np.argmax(d))
    out: list[float | None] = [None, None]
    for side, rng_ in enumerate((range(i, 0, -1), range(i, len(a) - 1))):
        for j in rng_:
            k = j - 1 if side == 0 else j + 1
            if d[k] < level <= d[j]:
                f = (d[j] - level) / (d[j] - d[k])
                out[side] = float(a[j] + f * (a[k] - a[j]))
                break
    return out


def parabola_max(a_au: list[float], dlnl: list[float]) -> dict | None:
    """The vertex of the parabola through the grid maximum and its two neighbours (reported
    only); None when the maximum is at an end."""
    a, d = np.asarray(a_au, float), np.asarray(dlnl, float)
    i = int(np.argmax(d))
    if i == 0 or i == len(a) - 1:
        return None
    c = np.polyfit(a[i - 1 : i + 2], d[i - 1 : i + 2], 2)
    if c[0] >= 0:
        return None
    av = -c[1] / (2.0 * c[0])
    return {"a_au": float(av), "dlnl": float(np.polyval(c, av))}


def posterior_marginal_a(q: np.ndarray, w: np.ndarray) -> dict:
    """The replacement posterior's density in a on the grid: a histogram in ln a with bins one
    grid step wide (edges at the midpoints between grid points), ln p − max, effective samples
    per bin, and the posterior's 68 % and 95 % intervals in a."""
    from fork.inference.nested import weighted_quantiles

    ln_a = q[:, Q_LN_A]
    g = np.log(np.asarray(A_GRID_AU))
    mid = 0.5 * (g[1:] + g[:-1])
    edges = np.concatenate([[g[0] - (mid[0] - g[0])], mid, [g[-1] + (g[-1] - mid[-1])]])
    mass, _ = np.histogram(ln_a, bins=edges, weights=w)
    sq, _ = np.histogram(ln_a, bins=edges, weights=w * w)
    n_eff = np.where(sq > 0, mass**2 / np.where(sq > 0, sq, 1.0), 0.0)
    dens = mass / np.diff(edges)
    with np.errstate(divide="ignore"):
        ln_p = np.where(dens > 0, np.log(np.where(dens > 0, dens, 1.0)), -np.inf)
    ln_p = ln_p - ln_p.max()
    a = np.exp(ln_a)
    q16, q50, q84 = weighted_quantiles(a, w, (0.16, 0.5, 0.84))
    q025, q975 = weighted_quantiles(a, w, (0.025, 0.975))
    return {"grid_au": list(A_GRID_AU), "bin_edges_ln_a": edges.tolist(),
            "ln_p_rel": [None if not np.isfinite(v) else float(v) for v in ln_p],
            "n_eff": n_eff.tolist(), "mass": mass.tolist(),
            "a_au": {"median": q50, "lo68": q16, "hi68": q84, "lo95": q025, "hi95": q975},
            "tight": bool(q975 - q025 < TIGHT_AU), "tight_below_au": TIGHT_AU,
            "min_eff_samples": MIN_EFF_SAMPLES}  # fmt: skip


def compare_with_posterior(reading: dict, intervals: dict, marg: dict) -> dict:
    """M9.md's profile-against-posterior table for dataset A."""
    post = marg["a_au"]
    lo2, hi2 = intervals["dlnl_minus_2"]
    out = {
        "profile": reading["reading"],
        "posterior_tight": marg["tight"],
        "posterior_95": [post["lo95"], post["hi95"]],
        "profile_dlnl_minus_2": [lo2, hi2],
    }
    if reading["reading"] == "constrained" and lo2 is not None and hi2 is not None:
        within = lo2 <= post["lo95"] and post["hi95"] <= hi2
        ratio = (post["hi95"] - post["lo95"]) / (hi2 - lo2)
        out |= {"posterior_95_within_profile": within, "width_ratio_posterior_over_profile": ratio}
        if within and 0.5 <= ratio <= 2.0:
            out["reading"] = (
                "constrained, posterior within and widths within a factor 2: the "
                "data set the posterior's a; the headline interval stands as a "
                "measurement"
            )
        elif ratio > 2.0:
            out["reading"] = (
                "constrained, posterior more than twice as wide: the volume of the "
                "other parameters broadens the posterior beyond what the data fix "
                "at the best fit; see the volume term"
            )
        else:
            out["reading"] = "constrained; neither row of the table applies as written"
    elif reading["reading"] == "flat":
        out["reading"] = "flat profile" + (
            ", tight posterior: the tightness comes from prior "
            "volume or sampling, not the data; the volume term "
            "says which"
            if marg["tight"]
            else ", posterior not tight: no row of the table"
        )
    else:
        out["reading"] = "partial: no reading beyond the curves"
    return out


AB_TABLE = {
    ("constrained", "constrained"): "The modern series alone fix the distance. The ancient "
    "sightings are not what pins a, and Le Verrier's 36 AU is excluded by the data either way.",
    ("constrained", "flat"): "The distance comes from the pre-1781 sightings: the 155-year lever "
    "arm. Without them the data allow any a from 26 to 44 AU at Δ ln L > −2.",
    ("flat", "flat"): "The data do not fix the distance. The headline's a interval then comes "
    "from the prior and the posterior's volume, not from a likelihood peak.",
    ("flat", "constrained"): "Removing data sharpened the profile: the ancient rows pull against "
    "the modern ones. No physical reading is offered; a tension to examine.",
}


def neptune_osculating_a_1781_au() -> float:
    """Neptune's osculating a at the 1781 epoch, about the inner barycentre with the model's
    μ = G(M_central + m), from the cached Horizons state. Context only."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.model.bodies import M_CENTRAL_MSUN, NEPTUNE, G

    st = about("neptune", [A1_EPOCH_JD_TDB])[0]
    r = float(np.linalg.norm(st[:3]))
    v2 = float(st[3:] @ st[3:])
    return 1.0 / (2.0 / r - v2 / (G * (M_CENTRAL_MSUN + NEPTUNE.mass_msun)))


def _group_chi2(obs, theta, ctx) -> dict:
    from fork.model.observation import GROUPS

    z = obs.whitened(theta, ctx)
    m = obs.n_measurements
    out = {}
    for i, g in enumerate(GROUPS):
        sel = obs.group_index == i
        out[g] = float(np.sum(z[:m][sel] ** 2) + z[m + i] ** 2)
    return out


def summarize_dataset(dataset: str, cache: Path = CACHE) -> dict | None:
    """One data set's profile from the cache: per point the exact value and its fit, the
    reading, the maximum, the intervals and the reference points."""
    from fork.inference.priors import free_a_space
    from fork.model.forward import ForwardContext

    cont = _read(_path(dataset, "continuation.json", cache))
    if cont is None:
        return None
    obs = observations(dataset)
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
        ok = [r for r in last["polished"] if r["status"] > 0]
        rec |= {"status": "pass" if pol["pass"] else "missing", "rerun": len(pol["attempts"]) > 1,
                "n_starts": len(starts["optima"]),
                "perturbative_best_chi2": cont["best"][str(k)]["chi2"],
                "candidates": [{"origin": c["origin"], "chi2": c["chi2"], "e": c["e"],
                                "galle_lon_deg": c["galle_lon_deg"]} for c in last["candidates"]],
                "polished": [{kk: r[kk] for kk in ("from", "chi2", "chi2_perturbative", "status",
                                                   "nfev", "mass_earth", "e", "inc_deg",
                                                   "galle_lon_deg", "active_bounds",
                                                   "u_box_half_widths", "u_vs_u_ref_std")}
                             for r in last["polished"]]}  # fmt: skip
        if pol["pass"] and ok:
            b = min(ok, key=lambda r: r["chi2"])
            phi = phi_of([b["log10_m"], a, b["k"], b["h"], b["lam0"]], [b["inc"], b["node"]])
            theta = space.theta(np.concatenate([b["u"], phi]))
            rec |= {"chi2": b["chi2"], "ln_l": -0.5 * b["chi2"], "mass_earth": b["mass_earth"],
                    "e": b["e"], "inc_deg": b["inc_deg"], "galle_lon_deg": b["galle_lon_deg"],
                    "active_bounds": b["active_bounds"],
                    "u_box_half_widths": b["u_box_half_widths"],
                    "u_vs_u_ref_std": b["u_vs_u_ref_std"],
                    "u_flagged": bool(b["u_vs_u_ref_std"] > U_FLAG_SIGMA),
                    "chi2_perturbative_minus_exact": b["chi2_perturbative"] - b["chi2"],
                    "group_chi2": _group_chi2(obs, theta, ctx)}  # fmt: skip
        pts.append(rec)
    grid = [p for p in pts if not p["reference"]]
    have = [p for p in grid if "ln_l" in p]
    ln_max = max(p["ln_l"] for p in have) if have else None
    for p in pts:
        p["dlnl"] = None if ("ln_l" not in p or ln_max is None) else p["ln_l"] - ln_max
    a_g = [p["a_au"] for p in grid]
    d_g = [p["dlnl"] for p in grid]
    reading = read_profile(a_g, d_g)
    used = [(a, d) for a, d in zip(a_g, d_g, strict=True) if d is not None]
    au, du = [a for a, _ in used], [d for _, d in used]
    intervals = {
        "dlnl_minus_0p5": crossings(au, du, -0.5) if used else [None, None],
        "dlnl_minus_2": crossings(au, du, DLNL_THRESHOLD) if used else [None, None],
    }
    a_max = max(have, key=lambda p: p["ln_l"]) if have else None
    refs = {}
    for p in pts:
        if p["reference"] and "ln_l" in p:
            refs[p["name"]] = {"a_au": p["a_au"], "dlnl": p["dlnl"],
                               "above_grid_maximum": p["dlnl"] > 0, "mass_earth": p["mass_earth"],
                               "galle_lon_deg": p["galle_lon_deg"], "e": p["e"]}  # fmt: skip
    ancient = {p["a_au"]: p["group_chi2"]["ancient:ra"] + p["group_chi2"]["ancient:dec"]
               for p in have if "group_chi2" in p}  # fmt: skip
    return {"dataset": dataset, "n_measurements": obs.n_measurements,
            "n_observations": int(obs.obs_id.size), "points": pts, "reading": reading,
            "a_max": None if a_max is None else {"a_au": a_max["a_au"],
                                                 "mass_earth": a_max["mass_earth"]},
            "parabola": parabola_max(au, du) if len(used) >= 3 else None,
            "intervals": intervals, "reference_points": refs,
            "continuation": {k: cont[k] for k in ("passes", "converged", "n_neighbour_runs",
                                                  "sweep_pairs", "sweep_improved")},
            "failed_points": [p["a_au"] for p in pts if p.get("status") == "missing"],
            "ancient_group_chi2_by_a": ancient}  # fmt: skip


def report(cache: Path = CACHE, out_json: Path = OUT_JSON, out_png: Path | None = OUT_PNG) -> dict:
    """docs/notes/M9_profile.json (and the figure) from the cached runs of both data sets."""
    check_frozen()
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()  # fmt: skip
    data = {ds: summarize_dataset(ds, cache) for ds in DATASETS}
    q, w, _ = replacement_posterior()
    marg = posterior_marginal_a(q, w)
    pub = json.loads(REPLACEMENT_JSON.read_text())["a_au"]
    marg["published_a_au"] = pub
    out: dict = {"description": "M9: profile likelihood in the semi-major axis, pre-registered in "
                 "docs/notes/M9.md (2026-11-09). Diagnostic only; no published number changes.",
                 "code_commit": head, "frozen": FROZEN,
                 "settings": {"grid_au": A_GRID_AU, "reference_au": REFERENCE_AU,
                              "n_random": N_RANDOM, "n_posterior": N_POSTERIOR,
                              "n_nearest": N_NEAREST, "n_rerun": N_RERUN, "e_max": E_MAX,
                              "inc_max_deg": INC_MAX_DEG, "continuation_tol_chi2": CONT_TOL_CHI2,
                              "polish_window_chi2": POLISH_WINDOW_CHI2, "max_polish": MAX_POLISH,
                              "distinct": {"galle_lon_deg": DISTINCT_LON_DEG, "e": DISTINCT_E},
                              "reading": {"dlnl": DLNL_THRESHOLD, "x_au": X_AU,
                                          "flat_range_au": FLAT_RANGE_AU}},
                 "datasets": data, "posterior_marginal_a": marg,
                 "neptune_osculating_a_1781_au": neptune_osculating_a_1781_au()}  # fmt: skip
    a_all = data["all"]
    if a_all is not None:
        vol = []
        for p in a_all["points"]:
            if p["reference"]:
                continue
            i = A_GRID_AU.index(p["a_au"])
            ok = (
                p["dlnl"] is not None
                and marg["ln_p_rel"][i] is not None
                and marg["n_eff"][i] >= MIN_EFF_SAMPLES
            )
            vol.append({"a_au": p["a_au"], "volume_term": (marg["ln_p_rel"][i] - p["dlnl"])
                        if ok else None, "n_eff": marg["n_eff"][i]})  # fmt: skip
        out["volume_term"] = vol
        out["profile_vs_posterior"] = compare_with_posterior(
            a_all["reading"], a_all["intervals"], marg
        )
    if a_all is not None and data["post1781"] is not None:
        key = (a_all["reading"]["reading"], data["post1781"]["reading"]["reading"])
        out["interpretation"] = {
            "a": key[0],
            "b": key[1],
            "text": AB_TABLE.get(
                key,
                "A profile is partial: the table is "
                "not used; both curves are reported "
                "and no conclusion is drawn.",
            ),
        }
    ledger_rows = [json.loads(line) for line in (ROOT / "data" / "cache" / "compute_ledger.jsonl")
                   .read_text().splitlines() if '"m9/' in line]  # fmt: skip
    out["ledger"] = [
        {
            k: r.get(k)
            for k in ("job", "n_likelihood_calls", "n_integrations", "wall_s", "cpu_s", "n_workers")
        }
        for r in ledger_rows
    ]
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(out, indent=1, default=float))
    if out_png is not None:
        from fork.export.m9_plot import draw

        draw(out, out_png)
    print(f"wrote {out_json}")
    return out


# ---------------------------------------------------------------------------- the screen


def _page_line(row_id: str) -> tuple[str, int]:
    parts = row_id.split("-")  # LV1846-T77-p129-l01-ra
    return parts[2].removeprefix("p"), int(parts[3].removeprefix("l"))


def screen_table(group_index: np.ndarray, resid_arcsec: np.ndarray, rows: list[dict],
                 sigma_by_group: dict[str, float], hand_checked: set[tuple[str, int]],
                 threshold: float = SCREEN_SIGMA) -> list[dict]:  # fmt: skip
    """Every measurement more than ``threshold`` × σ_g from its group's median O − JPL, with
    its page and line, source, date, quantity, O − JPL, deviation, whether it was among the
    hand-checked rows, any logged misprint, and Le Verrier's own O − C beside its neighbours'
    (the previous and next rows of the same quantity in table order)."""
    from fork.model.observation import GROUPS

    g = np.asarray(group_index)
    r = np.asarray(resid_arcsec, dtype=float)
    med = {i: float(np.median(r[g == i])) for i in np.unique(g)}
    out = []
    for i, row in enumerate(rows):
        gi = int(g[i])
        s = sigma_by_group[GROUPS[gi]]
        dev = (r[i] - med[gi]) / s
        if abs(dev) <= threshold:
            continue
        page, line = _page_line(row["id"])
        same = [j for j, rr in enumerate(rows) if rr["quantity"] == row["quantity"]]
        pos = same.index(i)
        nb = [rows[same[j]]["lv_residual"] for j in (pos - 1, pos + 1) if 0 <= j < len(same)]
        out.append({"id": row["id"], "page": page, "line": line, "source": row["source_ref"],
                    "date": row["epoch"], "quantity": row["quantity"],
                    "o_minus_jpl_arcsec": float(r[i]), "group": GROUPS[gi],
                    "group_median_arcsec": med[gi], "sigma_g_arcsec": s, "deviation_sigma": dev,
                    "hand_checked": (page, line) in hand_checked,
                    "misprint_logged": row["misprints_applied"].strip() or None,
                    "le_verrier_o_minus_c_arcsec": row["lv_residual"],
                    "neighbours_o_minus_c_arcsec": nb})  # fmt: skip
    return out


def screen(out_json: Path = SCREEN_JSON) -> list[dict]:
    """The transcription screen on the real residuals (report only; the fit is unchanged)."""
    from fork.data.calibration import real_observations
    from fork.data.spotcheck import read_selection
    from fork.model import observation as om

    check_frozen()
    rows = om.read_rows()
    em = om.ErrorModel.load(ROOT / m6c.VARIANT_CONFIG["headline"]["error_model"])
    obs = real_observations(rows, em)
    resid = obs.residuals_arcsec(*om.true_uranus_ra_dec(obs))
    picks = {(p.page, int(p.line)) for p in read_selection(SPOTCHECK_ROWS)}
    table = screen_table(obs.group_index, resid, rows, em.sigma_arcsec, picks)
    payload = {"description": "M9 transcription screen: measurements more than 4 σ_g from their "
               "group's median O − JPL (docs/notes/M9.md). Report only; nothing is changed.",
               "threshold_sigma": SCREEN_SIGMA, "n_measurements": obs.n_measurements,
               "n_hand_checked_rows": len(picks), "flagged": table}  # fmt: skip
    out_json.write_text(json.dumps(payload, indent=1, default=float))
    cols = ("id", "source", "date", "quantity", "o_minus_jpl_arcsec", "deviation_sigma",
            "hand_checked", "misprint_logged", "le_verrier_o_minus_c_arcsec",
            "neighbours_o_minus_c_arcsec")  # fmt: skip
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))
    for t in table:
        print(
            "| "
            + " | ".join(f"{t[c]:.2f}" if isinstance(t[c], float) else str(t[c]) for c in cols)
            + " |"
        )
    print(f"wrote {out_json}")
    return table


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
    sub.add_parser("screen")
    a = ap.parse_args(argv)
    if a.cmd == "profile":
        cache = a.cache.resolve()
        if a.smoke and cache == CACHE.resolve():
            raise SystemExit("smoke runs go to a scratch --cache, never the run cache")
        profile(a.dataset, a.n_workers, cache, smoke=a.smoke)
    elif a.cmd == "report":
        report(a.cache.resolve(), out_png=None if a.no_png else OUT_PNG)
    else:
        screen()


if __name__ == "__main__":
    main()
