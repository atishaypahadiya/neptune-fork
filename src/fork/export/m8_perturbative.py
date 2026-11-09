"""M8 amendment, A8a (error-budget diagnostic) and A8e (per call) for the perturbative likelihood.

    uv run python -m fork.export.m8_perturbative a8a [--n-workers 8]   # resumable
    uv run python -m fork.export.m8_perturbative a8e
    uv run python -m fork.export.m8_perturbative summarize   # -> site/data/m8_perturbative.json

Spec: docs/PHASE1_SPEC.md, "M8 amendment" (5fc089f). Pointwise A8a is a diagnostic there, not
the gate; the gate is the reweighting ESS, measured once there is a fast posterior to reweight.
Nothing here samples.

**Points.** Stage 2's A8a points and their full-likelihood references (data/cache/m8/a8a: R2
core and ridge, S1, the whole prior and the ridge-branch region), plus 300 S2 points drawn here
from the S2 posterior with references computed the same way (full likelihood minimised over u,
exact-Hessian Laplace marginal).

**References u_ref**, one per data set: the median u of its M6c posterior (headline: the
replacement posterior, R2 as cell 2). Δs = 60 d, cubic quadrature.

**Terms at each point** (fast = the perturbative likelihood; each with the second-order switch
off and on; Δχ² and Δ(−2 ln L) for each, signed, fast-with-the-approximation minus without):
- total: fast against the full references;
- E_K: Keplerian perturber against the full integrator's perturber (tabulated at the nodes);
- E_JS: the forcing left out because Jupiter and Saturn respond to the perturber in the full
  model: the tabulated perturber without against with that forcing (not in the spec's list;
  found while building, docs/notes/M8.md);
- E_Δs: Δs = 60 d against 15 d;
- E_ref: u_ref against u_ref + 1 posterior std in each coordinate (all four kept);
- E₂ (switch-independent, full integrator only): the full model's z linearized in the mass,
  z(0) + (m/ε)[z(ε) − z(0)] with ε = 10⁻³ m, minimised over u, against the full χ²_min;
- the rest: total − E_K − E_JS, i.e. what is left once the perturber's and Jupiter's and
  Saturn's trajectories are exact: second order in mass (E₂ with the switch off), the u
  linearization about u_ref, and the quadrature.

**ESS forecast.** At the posterior points (R2, S1, S2: draws from the full posterior) with
Δ = ln L_(b) − ln L_fast (mode (b), the reweighting's full likelihood): the expected reweighting
ESS/N is 1 / (E_full[e^Δ] · E_full[e^−Δ]), an identity, not a Gaussian approximation. A forecast
of the gate, not the gate.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from fork.export import m8
from fork.export.to_json import write_json

CACHE = m8.CACHE / "perturbative"
OUT_NAME = "m8_perturbative.json"
S2_RNG = (8, 1846, 2)
N_S2 = 300
CHUNK = 25
FINE_STEP_DAYS = 15.0
E2_EPS = 1e-3
VARIANT_OF = {"s1": "s1_e05", "s2": "s2_sigma_data"}
VARIANTS = ("headline", "s1_e05", "s2_sigma_data")
INNER_OUTER_A = 30.5  # ridge branches (M6.md): inner a < 30.5 AU, outer a ≥ 30.5 AU


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with (CACHE / "progress.log").open("a") as fh:
        fh.write(line + "\n")


def references() -> dict:
    """{variant: (u_ref, u_std)}, cached (data-only: each variant's M6c posterior)."""
    from fork.model.perturbative import reference_from_m6c

    path = CACHE / "references.json"
    if path.exists():
        d = json.loads(path.read_text())
        return {k: (np.array(v["u_ref"]), np.array(v["u_std"])) for k, v in d.items()}
    out = {v: reference_from_m6c(v) for v in VARIANTS}
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({k: {"u_ref": a.tolist(), "u_std": b.tolist()}
                                for k, (a, b) in out.items()}, indent=1))  # fmt: skip
    return out


def s2_points() -> list[dict]:
    from fork.export import m6c
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult

    rng = np.random.default_rng(list(S2_RNG))
    res = PartitionedResult(tuple(NestedResult.load(m6c.Job("s2_sigma_data", k).path(".pkl"))
                                  for k in range(len(m6c.variant_cells("s2_sigma_data")))),
                            err_scale=ERR_INFLATION)  # fmt: skip
    eq = res.equal_weight(rng)
    _, _, space, _ = m8.problem("headline")
    out = [{"source": "s2", "q": eq[i].tolist()} for i in rng.choice(len(eq), N_S2, replace=False)]
    for p in out:
        p["class"] = m8.classify(space.theta(np.array(p["q"])))
    return out


def all_points() -> list[dict]:
    """Stage 2's records (with their references), then the S2 points (references to compute)."""
    return m8.load_a8a() + s2_points()


_W: dict = {}


def _init_worker(refs: dict) -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from fork.model.perturbative import PerturbativeProblem

    for v in VARIANTS:
        obs, ctx, space, lp = m8.problem(v)
        u_ref, u_std = refs[v]
        shifts = []
        for i in range(4):
            d = np.zeros(4)
            d[i] = u_std[i]
            shifts.append(PerturbativeProblem.build(obs, ctx, space, u_ref + d))
        _W[v] = {"obs": obs, "ctx": ctx, "space": space, "lp": lp,
                 "pr": PerturbativeProblem.build(obs, ctx, space, u_ref),
                 "fine": PerturbativeProblem.build(obs, ctx, space, u_ref,
                                                   step_days=FINE_STEP_DAYS),
                 "shifts": shifts}  # fmt: skip


def _pair(a, b) -> list[float]:
    """[Δχ², Δ(−2 ln L)] of marginal a minus marginal b."""
    return [a.chi2_min - b.chi2_min, a.m2lnl - b.m2lnl]


def _e2_chi2(w: dict, phi: np.ndarray, u_start: np.ndarray) -> tuple[float, int]:
    """χ²_min over u of the full model's z linearized in the perturber's mass (E₂)."""
    from fork.inference.fit import uranus_x_to_theta
    from fork.model.linearized import profile_min

    obs, ctx, space, lp = w["obs"], w["ctx"], w["space"], w["lp"]
    log10_eps = math.log10(E2_EPS)

    def z_lin(u):
        th = space.theta(np.concatenate([u, phi]))
        th_eps = th.copy()
        th_eps[6] += log10_eps
        z0 = obs.whitened(uranus_x_to_theta(np.concatenate([u, space.pq])), ctx)
        ze = obs.whitened(th_eps, ctx)
        return z0 + (ze - z0) / E2_EPS

    _, chi2, n = profile_min(z_lin, u_start, lp.sigma_raw)
    return chi2, 3 * n


def _evaluate(point: dict) -> dict:
    variant = point.get("variant") or VARIANT_OF.get(point["source"], "headline")
    w = _W[variant]
    space, lp, pr = w["space"], w["lp"], w["pr"]
    q = np.array(point["q"])
    phi = q[4:]
    t0 = time.process_time()
    rec = {k: point[k] for k in ("source", "class", "q", "cell") if k in point}
    th = space.theta(q)
    rec |= {"variant": variant, "a_au": float(th[7]), "e": float(th[8]),
            "mass_earth": float(10.0 ** th[6])}  # fmt: skip
    box = space.uranus_box
    if "chi2_full_min" in point:
        chi2_full, lap = point["chi2_full_min"], point["m2lnl_laplace"]
        offset = np.asarray(point["u_min_box_sigma"]) * np.asarray(box.sigma_scaled)
        u_min = np.asarray(box.centre) + offset
        m2lnl_b = point["b"]["m2lnl"]
    else:
        b = lp.marginal(phi, "b")
        u_min, chi2_full, _ = lp.profile(phi, b.u_star)
        lap = lp.laplace(phi, u_min)
        m2lnl_b = b.m2lnl
    rec |= {"chi2_full_min": chi2_full, "m2lnl_laplace": lap, "m2lnl_b": m2lnl_b}

    x_p, dx_j, dx_s = pr.integrated_states(phi)
    extra = pr.js_response_accel(dx_j, dx_s)
    for so in (False, True):
        tag = "so" if so else "fo"
        f = pr.evaluate(phi, second_order=so)
        fm = f.marginal
        tab = pr.evaluate(phi, second_order=so, perturber_pos=x_p).marginal
        tab_js = pr.evaluate(phi, second_order=so, perturber_pos=x_p, extra_accel=extra).marginal
        rec[tag] = {
            "chi2_min": fm.chi2_min, "m2lnl": fm.m2lnl,
            "total": [fm.chi2_min - chi2_full, fm.m2lnl - lap],
            "e_k": _pair(fm, tab), "e_js": _pair(tab, tab_js),
            "e_ds": _pair(fm, w["fine"].evaluate(phi, second_order=so).marginal),
            "e_ref": [_pair(fm, s.evaluate(phi, second_order=so).marginal) for s in w["shifts"]],
            "u_star_box_sigma": pr.box_sigma_offset(fm.u_star).tolist(),
            "outside_box": f.outside_box, "truncation_log10_bound": f.truncation_log10,
        }  # fmt: skip
    chi2_lin, n_int = _e2_chi2(w, phi, u_min)
    rec["e2_chi2"] = chi2_lin - chi2_full
    rec["e2_n_integrations"] = n_int
    rec["cpu_s"] = time.process_time() - t0
    return rec


def _chunk_path(k: int, sub: str = "a8a") -> Path:
    return CACHE / sub / f"chunk_{k:03d}.json"


def _run_chunk(chunk: list[dict]) -> list[dict]:
    return [_evaluate(p) for p in chunk]


def _run_points(points: list[dict], sub: str, n_workers: int) -> None:
    """Evaluate ``points`` in resumable chunks under CACHE/``sub``."""
    from fork import ledger

    refs = references()
    chunks = [points[i : i + CHUNK] for i in range(0, len(points), CHUNK)]
    todo = [k for k in range(len(chunks)) if not _chunk_path(k, sub).exists()]
    _log(f"{sub}: {len(points)} points, {len(chunks)} chunks, {len(todo)} to do, "
         f"{n_workers} workers")  # fmt: skip
    _chunk_path(0, sub).parent.mkdir(parents=True, exist_ok=True)
    t0, c0 = time.time(), ledger.cpu_now()
    with ProcessPoolExecutor(n_workers, initializer=_init_worker, initargs=(refs,)) as ex:
        futs = {k: ex.submit(_run_chunk, chunks[k]) for k in todo}
        for i, (k, f) in enumerate(futs.items()):
            recs = f.result()
            tmp = _chunk_path(k, sub).with_suffix(".tmp")
            tmp.write_text(json.dumps(recs))
            tmp.replace(_chunk_path(k, sub))  # atomic: a resumed run never reads half a chunk
            _log(f"{sub}: chunk {k} done ({i + 1}/{len(todo)})")
    if todo:
        ledger.append(f"m8/perturbative_{sub}", "diagnostic", n_likelihood_calls=None,
                      n_integrations=None, sim_years=None, wall_s=time.time() - t0,
                      cpu_s=ledger.cpu_now() - c0, n_workers=n_workers,
                      n_points=len(points))  # fmt: skip


def _load_points(n_points: int, sub: str) -> list[dict]:
    n = math.ceil(n_points / CHUNK)
    missing = [k for k in range(n) if not _chunk_path(k, sub).exists()]
    if missing:
        raise SystemExit(f"{sub} incomplete: chunks {missing} missing")
    return [r for k in range(n) for r in json.loads(_chunk_path(k, sub).read_text())]


def run_a8a(n_workers: int) -> None:
    _run_points(all_points(), "a8a", n_workers)


def load_a8a() -> list[dict]:
    return _load_points(len(all_points()), "a8a")


# ---------------------------------------------------------------------------- per-cell draws

CELLS_RNG = (8, 1846, 3)
N_PER_CELL_DRAW = 300
CELL_TAG = {"headline": "h", "s1_e05": "s1", "s2_sigma_data": "s2"}
CELLS_DRAWN = (0, 1, 3)  # cell 2 is covered by the R2, S1 and S2 points above


def cell_index(variant: str, a_au: float, e: float) -> int:
    """Which of the variant's M6c cells (a slice × e slice) holds (a, e)."""
    from fork.export import m6c

    for k, ((a_lo, a_hi), e_rng) in enumerate(m6c.variant_cells(variant)):
        if a_lo <= a_au < a_hi and (e_rng is None or e_rng[0] <= e < e_rng[1]):
            return k
    return -1


def cell_points() -> list[dict]:
    """300 equal-weight draws from each saved M6c run of cells 0, 1 and 3, for the headline,
    S1 and S2 (full-likelihood posteriors, one cell each)."""
    from fork.export import m6c
    from fork.inference.nested import NestedResult

    rng = np.random.default_rng(list(CELLS_RNG))
    _, _, space, _ = m8.problem("headline")
    out = []
    for v in VARIANTS:
        for k in CELLS_DRAWN:
            eq = NestedResult.load(m6c.Job(v, k).path(".pkl")).equal_weight(rng)
            idx = rng.choice(len(eq), N_PER_CELL_DRAW, replace=False)
            for i in idx:
                q = eq[i]
                out.append({"source": f"{CELL_TAG[v]}_c{k}", "variant": v, "cell": k,
                            "q": q.tolist(), "class": m8.classify(space.theta(q))})  # fmt: skip
    return out


def run_cells(n_workers: int) -> None:
    _run_points(cell_points(), "cells", n_workers)


def load_cells() -> list[dict]:
    return _load_points(len(VARIANTS) * len(CELLS_DRAWN) * N_PER_CELL_DRAW, "cells")


# ---------------------------------------------------------------------------- summaries


def _best_by_variant() -> dict:
    from fork.export import m6c
    from fork.inference.nested import NestedResult

    best = {"headline": m8.best_m2lnl()}
    best["s1_e05"] = best["headline"]  # same data and error model (Stage 2's domain)
    best["s2_sigma_data"] = float(min(
        -2.0 * np.max(NestedResult.load(m6c.Job("s2_sigma_data", k).path(".pkl")).results["logl"])
        for k in range(len(m6c.variant_cells("s2_sigma_data")))))  # fmt: skip
    return best


def _group(r: dict) -> str:
    if r["source"].startswith("prior"):
        return "prior_tail"
    if r["class"] == "ridge":
        return "ridge_inner" if r["a_au"] < INNER_OUTER_A else "ridge_outer"
    return r["class"]  # core or other_cells


def _ess_forecast(delta: np.ndarray, weights: np.ndarray | None = None) -> float | None:
    """1 / (E[e^Δ] E[e^−Δ]) over full-posterior draws (optionally weighted)."""
    if delta.size == 0:
        return None
    w = np.ones(delta.size) if weights is None else weights
    w = w / w.sum()
    s = delta - np.median(delta)  # the identity is shift-invariant; this keeps exp finite
    return float(1.0 / ((w @ np.exp(s)) * (w @ np.exp(-s))))


def _r2_stratum_mass() -> float:
    """The R2 posterior's mass with e ≥ 0.06 (to weight the stratified core / ridge draws)."""
    from fork.export import m6c_repeat
    from fork.inference.nested import NestedResult

    r2 = m6c_repeat.merged_cell2([NestedResult.load(j.path(".pkl"))
                                  for j in m6c_repeat.JOBS if j.group == "r2"])  # fmt: skip
    e = np.array([r2.space.theta(q)[8] for q in r2.samples])
    w = r2.weights
    return float(w[e >= m8.E_RIDGE].sum() / w.sum())


def summarize_a8a(recs: list[dict]) -> dict:
    best = _best_by_variant()
    for r in recs:
        r["in_domain"] = (r["chi2_full_min"] - best[r["variant"]]) / 2.0 <= m8.DOMAIN_DLNL
        r["group"] = _group(r)
    groups = ("core", "ridge_inner", "ridge_outer", "other_cells", "prior_tail")
    out = {"best_m2lnl": {k: round(v, 3) for k, v in best.items()},
           "n_points": len(recs),
           "points_by_source": {s: sum(r["source"] == s for r in recs)
                                for s in ("r2", "s1", "s2", "prior", "prior_ridge_region")},
           "tolerances_diagnostic": {"abs_d_chi2": m8.TOL_CHI2, "abs_d_m2lnl": m8.TOL_M2LNL,
                                     "domain_delta_lnl": m8.DOMAIN_DLNL},
           "groups": {}}  # fmt: skip

    def finite(r):
        return math.isfinite(r["m2lnl_laplace"])

    for g in groups:
        rs = [r for r in recs if r["group"] == g and finite(r)]
        dom = [r for r in rs if r["in_domain"]]
        entry = {"n": len(rs), "n_in_domain": len(dom),
                 "laplace_nan": sum(r["group"] == g and not finite(r) for r in recs)}  # fmt: skip
        for label, sub in (("in_domain", dom), ("all", rs)):
            if not sub:
                entry[label] = None
                continue
            e = {"e2_chi2": m8._dist(np.array([r["e2_chi2"] for r in sub]))}
            for tag in ("fo", "so"):
                t = {}
                for term in ("total", "e_k", "e_js", "e_ds"):
                    arr = np.array([r[tag][term] for r in sub])
                    t[term] = {"d_chi2": m8._dist(arr[:, 0]), "d_m2lnl": m8._dist(arr[:, 1])}
                rest = np.array([np.array(r[tag]["total"]) - r[tag]["e_k"] - r[tag]["e_js"]
                                 for r in sub])  # fmt: skip
                t["rest_after_k_js"] = {"d_chi2": m8._dist(rest[:, 0]),
                                        "d_m2lnl": m8._dist(rest[:, 1])}  # fmt: skip
                ref = np.array([r[tag]["e_ref"] for r in sub])  # (n, 4, 2)
                t["e_ref_max_over_shifts"] = {
                    "d_chi2": m8._dist(np.max(np.abs(ref[:, :, 0]), axis=1)),
                    "d_m2lnl": m8._dist(np.max(np.abs(ref[:, :, 1]), axis=1))}  # fmt: skip
                tot = np.array([r[tag]["total"] for r in sub])
                t["n_over_chi2_tol"] = int(np.sum(np.abs(tot[:, 0]) >= m8.TOL_CHI2))
                t["n_over_m2lnl_tol"] = int(np.sum(np.abs(tot[:, 1]) >= m8.TOL_M2LNL))
                t["n_outside_box"] = sum(r[tag]["outside_box"] for r in sub)
                t["max_truncation_log10_bound"] = round(max(r[tag]["truncation_log10_bound"]
                                                            for r in sub), 1)  # fmt: skip
                e[tag] = t
            entry[label] = e
        out["groups"][g] = entry

    # blind-spot screen (spec check 1): in-domain points where the fast ln L is below the full
    # by more than 3, i.e. Δ(−2 ln L) > 6
    out["blind_spot_candidates"] = {
        tag: [{"source": r["source"], "group": r["group"], "a_au": round(r["a_au"], 2),
               "e": round(r["e"], 3), "d_m2lnl": round(r[tag]["total"][1], 2)}
              for r in recs if r["in_domain"] and finite(r) and r[tag]["total"][1] > 6.0]
        for tag in ("fo", "so")}  # fmt: skip

    # ESS forecast at the full-posterior draws
    f_ridge = _r2_stratum_mass()
    fc: dict = {"r2_ridge_mass": round(f_ridge, 4)}
    for tag in ("fo", "so"):
        d = {}
        post = [r for r in recs if r["source"] in ("r2", "s1", "s2")]

        def delta(rs, tag=tag):
            return np.array([-0.5 * (r["m2lnl_b"] - r[tag]["m2lnl"]) for r in rs])

        r2 = [r for r in post if r["source"] == "r2"]
        core = [r for r in r2 if r["e"] < m8.E_RIDGE]
        ridge = [r for r in r2 if r["e"] >= m8.E_RIDGE]
        wts = np.array([(1 - f_ridge) / len(core) if r["e"] < m8.E_RIDGE
                        else f_ridge / len(ridge) for r in r2])  # fmt: skip
        d["r2_overall"] = _ess_forecast(delta(r2), wts)
        d["r2_core"] = _ess_forecast(delta(core))
        for name, sel in (("r2_ridge_inner", lambda r: r["a_au"] < INNER_OUTER_A),
                          ("r2_ridge_outer", lambda r: r["a_au"] >= INNER_OUTER_A)):  # fmt: skip
            sub = [r for r in ridge if sel(r)]
            d[name] = _ess_forecast(delta(sub))
            d[name + "_n"] = len(sub)
        for src in ("s1", "s2"):
            d[src + "_overall"] = _ess_forecast(delta([r for r in post if r["source"] == src]))
        d["delta_std"] = {k: round(float(np.std(delta([r for r in post if r["source"] == k]))), 4)
                          for k in ("r2", "s1", "s2")}  # fmt: skip
        fc[tag] = {k: (round(v, 4) if isinstance(v, float) else v) for k, v in d.items()}
    out["ess_forecast"] = fc
    out["ess_gates"] = {"overall": 0.5, "per_ridge_branch": 0.3}
    return out


# ---------------------------------------------------------------------------- A8e


def run_a8e() -> dict:
    """Per-call CPU time, one thread, on Stage 2's A8e points: full (M6c path), mode (b) (the
    reweighting's full likelihood), fast first order, fast second order."""
    from fork.export import m6c_repeat
    from fork.inference.nested import NestedResult
    from fork.model.perturbative import PerturbativeProblem

    refs = references()
    obs, ctx, space, lp = m8.problem("headline")
    t0 = time.process_time()
    pr = PerturbativeProblem.build(obs, ctx, space, refs["headline"][0])
    build_s = time.process_time() - t0
    r2 = m6c_repeat.merged_cell2([NestedResult.load(j.path(".pkl"))
                                  for j in m6c_repeat.JOBS if j.group == "r2"])  # fmt: skip
    rng = np.random.default_rng(list(m8.A8E_RNG))
    eq = r2.equal_weight(rng)
    qs = eq[rng.choice(len(eq), m8.A8E_N_PHI, replace=False)]
    calls = {"full": lambda q: -0.5 * obs.chi2(space.theta(q), ctx),
             "b": lambda q: lp.log_likelihood(q[4:], "b"),
             "fast_fo": lambda q: pr.log_likelihood(q[4:]),
             "fast_so": lambda q: pr.log_likelihood(q[4:], second_order=True)}  # fmt: skip
    for f in calls.values():
        f(qs[0])
    out = {"n_phi": m8.A8E_N_PHI, "repeats": 3, "ms_per_call": {},
           "precompute_cpu_s": round(build_s, 3), "n_nodes": int(pr.nodes_jd.size)}  # fmt: skip
    for name, f in calls.items():
        best = math.inf
        for _ in range(3):
            t0 = time.process_time()
            for q in qs:
                f(q)
            best = min(best, (time.process_time() - t0) / len(qs) * 1e3)
        out["ms_per_call"][name] = round(best, 4)
    ms = out["ms_per_call"]
    out["ratio_to_full"] = {k: round(v / ms["full"], 4) for k, v in ms.items()}
    out["per_call_target_ms"] = round(ms["full"] / 10.0, 3)
    n_cells, n_rw = 4, 6000
    out["reweighting_cpu_s_per_fit"] = round(n_cells * n_rw * ms["b"] / 1e3, 1)
    m6c_cpu_h = 28.5
    for tag in ("fast_fo", "fast_so"):
        sampling = m6c_cpu_h * ms[tag] / ms["full"]
        out[f"projected_fit_cpu_h_{tag}"] = {
            "sampling": round(sampling, 3),
            "reweighting": round(out["reweighting_cpu_s_per_fit"] / 3600.0, 3),
            "total": round(sampling + out["reweighting_cpu_s_per_fit"] / 3600.0, 3),
            "assumes": (
                "M6c's 28.5 CPU-h scaled by the per-call ratio (same number of "
                "calls); a projection, not the A8e measurement, which needs sampling"
            ),
        }
    out["a8e_target_cpu_h"] = 2.85
    _log(f"a8e: {json.dumps(ms)}")
    (CACHE / "a8e.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def summarize_cells(cell_recs: list[dict], a8a_recs: list[dict]) -> dict:
    """Per cell of the headline, S1 and S2: the error budget and the ESS forecast (the per-cell
    gate, ESS/N ≥ 0.3), plus the decision 2 switch rule. Cell 2 uses the A8a points: R2's
    stratified draws (weighted by the ridge mass) for the headline, S1's and S2's own draws
    that fall in their cell 2."""
    best = _best_by_variant()
    f_ridge = _r2_stratum_mass()
    per: dict = {}
    for v in VARIANTS:
        for k in range(4):
            if k == 2 and v == "headline":
                rs = [r for r in a8a_recs if r["source"] == "r2"]
                wts = np.array([(1 - f_ridge) if r["e"] < m8.E_RIDGE else f_ridge for r in rs])
                n_c = sum(r["e"] < m8.E_RIDGE for r in rs)
                wts = wts / np.where(np.array([r["e"] < m8.E_RIDGE for r in rs]), n_c,
                                     len(rs) - n_c)  # fmt: skip
                src = "R2 merged cell 2 (stratified core/ridge draws, weighted)"
            elif k == 2:
                tag = CELL_TAG[v]
                rs = [r for r in a8a_recs if r["source"] == tag
                      and cell_index(v, r["a_au"], r["e"]) == 2]  # fmt: skip
                wts = None
                src = f"{tag} posterior draws in cell 2"
            else:
                rs = [r for r in cell_recs if r.get("variant") == v and r.get("cell") == k]
                wts = None
                src = f"M6c {v} cell {k} run"
            e: dict = {"source": src, "n": len(rs)}
            dom = [r for r in rs if (r["chi2_full_min"] - best[v]) / 2.0 <= m8.DOMAIN_DLNL]
            e["n_in_domain"] = len(dom)
            for tag in ("fo", "so"):
                d = np.array([-0.5 * (r["m2lnl_b"] - r[tag]["m2lnl"]) for r in rs])
                t = {"ess_forecast": round(_ess_forecast(d, wts), 4),
                     "delta_std": round(float(np.std(d)), 4)}  # fmt: skip
                for term in ("total", "e_k", "e_js"):
                    t[term] = m8._dist(np.array([r[tag][term][0] for r in rs]))
                t["rest_after_k_js"] = m8._dist(
                    np.array(
                        [r[tag]["total"][0] - r[tag]["e_k"][0] - r[tag]["e_js"][0] for r in rs]
                    )
                )
                t["blind_spot_candidates_in_domain"] = sum(r[tag]["total"][1] > 6.0 for r in dom)
                t["n_outside_box"] = sum(r[tag]["outside_box"] for r in rs)
                t["max_truncation_log10_bound"] = round(
                    max(r[tag]["truncation_log10_bound"] for r in rs), 1)  # fmt: skip
                e[tag] = t
            per[f"{v}/c{k}"] = e
    rule = {"overall_min": 0.9, "branch_min": 0.5,
            "text": "switch off if the switch-off ESS forecast is >= 0.9 overall and >= 0.5 per "
                    "ridge branch, else on (decision 2, 2026-11-09)"}  # fmt: skip
    a8a = summarize_a8a(a8a_recs)["ess_forecast"]["fo"]
    overall = min(a8a["r2_overall"], a8a["s1_overall"], a8a["s2_overall"])
    branch = min(a8a["r2_ridge_inner"], a8a["r2_ridge_outer"])
    rule["applied"] = {"overall_forecast_min": overall, "branch_forecast_min": branch,
                       "switch": "off" if overall >= 0.9 and branch >= 0.5 else "on"}  # fmt: skip
    cell_fail_off = [c for c, e in per.items() if e["fo"]["ess_forecast"] < 0.3]
    rule["per_cell_gate"] = {
        "threshold": 0.3,
        "cells_failing_with_switch_off": cell_fail_off,
        "cells_failing_with_switch_on": [
            c for c, e in per.items() if e["so"]["ess_forecast"] < 0.3
        ],
    }
    return {"cells": per, "switch_rule": rule}


def summarize() -> None:
    from fork.model.perturbative import PerturbativeProblem

    refs = references()
    pre = {}
    for v in VARIANTS:
        obs, ctx, space, lp = m8.problem(v)
        pr = PerturbativeProblem.build(obs, ctx, space, refs[v][0])
        pre[v] = {"u_ref_box_sigma": lp.box_sigma_offset(refs[v][0]).round(3).tolist(),
                  "u_std_box_sigma": (refs[v][1] / np.asarray(space.uranus_box.sigma_scaled))
                  .round(3).tolist(),
                  "symplectic_rel_dev": float(f"{pr.symplectic_dev:.3g}"),
                  "n_nodes": int(pr.nodes_jd.size)}  # fmt: skip
    a8e_path = CACHE / "a8e.json"
    out = {"description": "M8 amendment: perturbative likelihood, A8a error budget (diagnostic) "
                          "and A8e per call; no sampling",
           "spec": "docs/PHASE1_SPEC.md, M8 amendment (5fc089f)",
           "precompute": pre,
           "a8a": summarize_a8a(load_a8a()),
           "per_cell": summarize_cells(load_cells(), load_a8a()),
           "a8e": json.loads(a8e_path.read_text()) if a8e_path.exists() else None}  # fmt: skip
    write_json(OUT_NAME, out)
    print(json.dumps({k: out[k] for k in ("precompute", "a8e")}, indent=1))


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="python -m fork.export.m8_perturbative")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("a8a", "cells"):
        a = sub.add_parser(name)
        a.add_argument("--n-workers", type=int, default=8)
    sub.add_parser("a8e")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "a8a":
        run_a8a(args.n_workers)
    elif args.cmd == "cells":
        run_cells(args.n_workers)
    elif args.cmd == "a8e":
        run_a8e()
    else:
        summarize()


if __name__ == "__main__":
    main()
