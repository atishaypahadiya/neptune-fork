"""M8 validation, Stage 2 (no sampling): A8a and the per-call part of A8e.

    uv run python -m fork.export.m8 a8a [--n-workers 8]   # resumable; progress in data/cache/m8
    uv run python -m fork.export.m8 a8e
    uv run python -m fork.export.m8 summarize             # -> site/data/m8_validation.json

Plan: docs/notes/M8.md. Spec: docs/PHASE1_SPEC.md, "M8". Reads the M6c inputs through
``fork.export.m6c`` without changing them (``check_frozen`` first) and the cached M6c / R2
nested results. Writes nothing under the M6c cache.

**A8a points** (Stage 2 instruction: ≥ 500 points from R2, S1 and prior draws):
- R2 posterior (merged cell 2): 300 core (e < 0.06) and 300 ridge (e ≥ 0.06), stratified;
- S1 posterior (its four cells combined): 300;
- headline prior: 150 over the whole prior, 150 in the ridge-branch region (a 28–33 AU,
  e 0.06–0.12).
At each: mode (a), mode (b), the full likelihood minimised over u (``profile_min``, started at
mode (b)'s u*), and the exact-Hessian Laplace marginal there. Spec tolerances: |Δχ²| < 0.1 for
the profile, |Δ(−2 ln L)| < 0.2 for the marginal, on points within Δ ln L = 30 of the best.
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

from fork.export.to_json import write_json

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "m8"
PROGRESS_LOG = CACHE / "progress.log"
OUT_NAME = "m8_validation.json"
A8A_RNG = (8, 1846, 1)
A8E_RNG = (8, 1846, 5)
N_R2_CORE, N_R2_RIDGE, N_S1, N_PRIOR, N_PRIOR_RIDGE = 300, 300, 300, 150, 150
CHUNK = 25
E_RIDGE = 0.06
RIDGE_REGION = ((28.0, 33.0), (0.06, 0.12))
TOL_CHI2, TOL_M2LNL, DOMAIN_DLNL = 0.1, 0.2, 30.0
A8E_N_PHI = 100


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with PROGRESS_LOG.open("a") as fh:
        fh.write(line + "\n")


def problem(variant: str = "headline"):
    """(obs, ctx, space, LinearizedProblem) on the M6c variant's data and box (uncut space)."""
    from fork.export import m6c
    from fork.inference.priors import free_a_space
    from fork.model.linearized import LinearizedProblem

    m6c.check_frozen(variant)
    obs, ctx, _, pq, box = m6c.dataset(variant=variant)
    space = free_a_space(pq, m6c.VARIANT_CONFIG[variant]["e_max"], box, inclined=True)
    return obs, ctx, space, LinearizedProblem.build(obs, ctx, space)


def classify(theta: np.ndarray) -> str:
    a, e = float(theta[7]), float(theta[8])
    if 28.0 <= a < 60.0 and e < 0.12:
        return "core" if e < E_RIDGE else "ridge"
    return "other_cells"


def best_m2lnl() -> float:
    """Lowest −2 ln L any M6c headline cell or R2 run reached (the Δ ln L = 30 domain)."""
    from fork.export import m6c_repeat
    from fork.inference.nested import NestedResult

    runs = list(m6c_repeat.replacement_cells())
    runs += [NestedResult.load(j.path(".pkl")) for j in m6c_repeat.JOBS if j.group == "r2"]
    return float(min(-2.0 * np.max(r.results["logl"]) for r in runs))


def a8a_points() -> list[dict]:
    """The A8a point list (deterministic): source, class and the 11-entry q of each."""
    from fork.export import m6c, m6c_repeat
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult

    rng = np.random.default_rng(list(A8A_RNG))
    _, _, space, _ = problem("headline")
    out: list[dict] = []

    r2 = m6c_repeat.merged_cell2([NestedResult.load(j.path(".pkl"))
                                  for j in m6c_repeat.JOBS if j.group == "r2"])  # fmt: skip
    eq = r2.equal_weight(rng)
    e = np.array([space.theta(q)[8] for q in eq])
    for sel, n in ((e < E_RIDGE, N_R2_CORE), (e >= E_RIDGE, N_R2_RIDGE)):
        idx = rng.choice(np.flatnonzero(sel), n, replace=False)
        out += [{"source": "r2", "q": eq[i].tolist()} for i in idx]

    s1 = PartitionedResult(tuple(NestedResult.load(m6c.Job("s1_e05", k).path(".pkl"))
                                 for k in range(len(m6c.variant_cells("s1_e05")))),
                           err_scale=ERR_INFLATION)  # fmt: skip
    eq = s1.equal_weight(rng)
    out += [{"source": "s1", "q": eq[i].tolist()}
            for i in rng.choice(len(eq), N_S1, replace=False)]  # fmt: skip

    (a_lo, a_hi), (e_lo, e_hi) = RIDGE_REGION
    ridge_space = space.with_a_slice(a_lo, a_hi).with_e_slice(e_lo, e_hi)
    for src, sp, n in (("prior", space, N_PRIOR),
                       ("prior_ridge_region", ridge_space, N_PRIOR_RIDGE)):  # fmt: skip
        out += [{"source": src, "q": sp.prior_transform(rng.random(sp.ndim)).tolist()}
                for _ in range(n)]  # fmt: skip
    for p in out:
        p["class"] = classify(space.theta(np.array(p["q"])))
    return out


_WORKER: dict = {}


def _init_worker() -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _WORKER["headline"] = problem("headline")
    _WORKER["s1_e05"] = problem("s1_e05")


def _evaluate(point: dict) -> dict:
    """Modes (a) and (b) against the full profile and the exact-Hessian Laplace marginal."""
    from fork.model.linearized import truncation_log10_bound

    variant = "s1_e05" if point["source"] == "s1" else "headline"
    _, _, space, prob = _WORKER[variant]
    q = np.array(point["q"])
    phi = q[4:]
    t0 = time.process_time()
    rec = dict(point)
    th = space.theta(q)
    rec |= {"a_au": float(th[7]), "e": float(th[8]), "mass_earth": float(10.0 ** th[6])}
    m = {mode: prob.marginal(phi, mode) for mode in ("a", "b")}
    u_min, chi2_full, n_prof = prob.profile(phi, m["b"].u_star)
    lap = prob.laplace(phi, u_min)
    box = space.uranus_box
    rec |= {"chi2_full_min": chi2_full, "m2lnl_laplace": lap, "n_profile_calls": n_prof,
            "u_min_box_sigma": prob.box_sigma_offset(u_min).tolist()}  # fmt: skip
    for mode, mm in m.items():
        rec[mode] = {"chi2_min": mm.chi2_min, "m2lnl": mm.m2lnl,
                     "d_chi2": mm.chi2_min - chi2_full, "d_m2lnl": mm.m2lnl - lap,
                     "u_star_box_sigma": prob.box_sigma_offset(mm.u_star).tolist(),
                     "truncation_log10_bound": truncation_log10_bound(mm.u_star, mm.cov_u,
                                                                      box.lo, box.hi),
                     "n_integrations": mm.n_integrations}  # fmt: skip
    rec["cpu_s"] = time.process_time() - t0
    return rec


def _chunk_path(k: int) -> Path:
    return CACHE / "a8a" / f"chunk_{k:03d}.json"


def run_a8a(n_workers: int) -> None:
    from fork import ledger

    points = a8a_points()
    chunks = [points[i : i + CHUNK] for i in range(0, len(points), CHUNK)]
    todo = [k for k in range(len(chunks)) if not _chunk_path(k).exists()]
    _log(f"a8a: {len(points)} points, {len(chunks)} chunks, {len(todo)} to do, "
         f"{n_workers} workers")  # fmt: skip
    _chunk_path(0).parent.mkdir(parents=True, exist_ok=True)
    t0, c0 = time.time(), ledger.cpu_now()
    n_int = 0
    with ProcessPoolExecutor(n_workers, initializer=_init_worker) as ex:
        futs = {k: ex.submit(_run_chunk, chunks[k]) for k in todo}
        for i, (k, f) in enumerate(futs.items()):
            recs = f.result()
            n_int += sum(r["a"]["n_integrations"] + r["b"]["n_integrations"]
                         + r["n_profile_calls"] + 33 for r in recs)  # fmt: skip
            tmp = _chunk_path(k).with_suffix(".tmp")
            tmp.write_text(json.dumps(recs))
            tmp.replace(_chunk_path(k))  # atomic: a resumed run never reads half a chunk
            _log(f"a8a: chunk {k} done ({i + 1}/{len(todo)})")
    if todo:
        ledger.append("m8/a8a", "diagnostic", n_likelihood_calls=None, n_integrations=n_int,
                      sim_years=None, wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0,
                      n_workers=n_workers, n_points=len(points))  # fmt: skip


def _run_chunk(chunk: list[dict]) -> list[dict]:
    return [_evaluate(p) for p in chunk]


def load_a8a() -> list[dict]:
    n = math.ceil(len(a8a_points()) / CHUNK)
    missing = [k for k in range(n) if not _chunk_path(k).exists()]
    if missing:
        raise SystemExit(f"a8a incomplete: chunks {missing} missing")
    return [r for k in range(n) for r in json.loads(_chunk_path(k).read_text())]


def _dist(x: np.ndarray) -> dict:
    a = np.abs(x)
    return {"n": int(x.size),
            "abs_q50": round(float(np.percentile(a, 50)), 5) if x.size else None,
            "abs_q90": round(float(np.percentile(a, 90)), 5) if x.size else None,
            "abs_q99": round(float(np.percentile(a, 99)), 5) if x.size else None,
            "abs_max": round(float(a.max()), 5) if x.size else None,
            "signed_median": round(float(np.median(x)), 5) if x.size else None}  # fmt: skip


def summarize_a8a(recs: list[dict]) -> dict:
    best = best_m2lnl()
    groups = {"core": lambda r: r["source"] in ("r2", "s1") and r["class"] == "core",
              "ridge": lambda r: r["source"] in ("r2", "s1") and r["class"] == "ridge",
              "other_cells": lambda r: r["source"] in ("r2", "s1")
              and r["class"] == "other_cells",
              "prior_tail": lambda r: r["source"].startswith("prior")}  # fmt: skip
    out = {"best_m2lnl": round(best, 3), "n_points": len(recs),
           "points_by_source": {s: sum(r["source"] == s for r in recs)
                                for s in ("r2", "s1", "prior", "prior_ridge_region")},
           "tolerances": {"abs_d_chi2": TOL_CHI2, "abs_d_m2lnl": TOL_M2LNL,
                          "domain_delta_lnl": DOMAIN_DLNL},
           "groups": {}, "modes": {}}  # fmt: skip
    for name, sel in groups.items():
        rs = [r for r in recs if sel(r)]
        dom = [r for r in rs if (r["chi2_full_min"] - best) / 2.0 <= DOMAIN_DLNL]
        g = {"n": len(rs), "n_in_domain": len(dom),
             "laplace_nan": sum(not math.isfinite(r["m2lnl_laplace"]) for r in rs)}  # fmt: skip
        for mode in ("a", "b"):
            for label, sub in (("in_domain", dom), ("all", rs)):
                ok = [r for r in sub if math.isfinite(r["m2lnl_laplace"])]
                dc = np.array([r[mode]["d_chi2"] for r in ok])
                dm = np.array([r[mode]["d_m2lnl"] for r in ok])
                us = np.array([np.abs(r[mode]["u_star_box_sigma"]) for r in ok]).reshape(-1, 4)
                g[f"{mode}_{label}"] = {
                    "d_chi2": _dist(dc), "d_m2lnl": _dist(dm),
                    "n_over_chi2_tol": int(np.sum(np.abs(dc) >= TOL_CHI2)),
                    "n_over_m2lnl_tol": int(np.sum(np.abs(dm) >= TOL_M2LNL)),
                    "max_abs_u_star_box_sigma": (us.max(axis=0).round(2).tolist()
                                                 if len(us) else None),
                    "max_truncation_log10_bound": (round(max(r[mode]["truncation_log10_bound"]
                                                             for r in ok), 1) if ok else None),
                }  # fmt: skip
        out["groups"][name] = g
    spec = [r for r in recs if r["source"] in ("r2", "s1")
            and (r["chi2_full_min"] - best) / 2.0 <= DOMAIN_DLNL
            and math.isfinite(r["m2lnl_laplace"])]  # fmt: skip
    for mode in ("a", "b"):
        dc = np.array([r[mode]["d_chi2"] for r in spec])
        dm = np.array([r[mode]["d_m2lnl"] for r in spec])
        out["modes"][mode] = {"n_spec_points": len(spec),
                              "max_abs_d_chi2": round(float(np.abs(dc).max()), 5),
                              "max_abs_d_m2lnl": round(float(np.abs(dm).max()), 5),
                              "pass": bool(np.abs(dc).max() < TOL_CHI2
                                           and np.abs(dm).max() < TOL_M2LNL)}  # fmt: skip
    out["mode_chosen"] = ("a" if out["modes"]["a"]["pass"]
                          else "b" if out["modes"]["b"]["pass"] else None)  # fmt: skip
    return out


def run_a8e() -> dict:
    """Per-call CPU time, one thread: the M6c full likelihood, mode (a), mode (b)."""
    from fork.export import m6c_repeat
    from fork.inference.nested import NestedResult

    obs, ctx, space, prob = problem("headline")
    r2 = m6c_repeat.merged_cell2([NestedResult.load(j.path(".pkl"))
                                  for j in m6c_repeat.JOBS if j.group == "r2"])  # fmt: skip
    rng = np.random.default_rng(list(A8E_RNG))
    eq = r2.equal_weight(rng)
    qs = eq[rng.choice(len(eq), A8E_N_PHI, replace=False)]
    calls = {"full": lambda q: -0.5 * obs.chi2(space.theta(q), ctx),
             "a": lambda q: prob.log_likelihood(q[4:], "a"),
             "b": lambda q: prob.log_likelihood(q[4:], "b")}  # fmt: skip
    for f in calls.values():  # warm the numba caches
        f(qs[0])
    out = {"n_phi": A8E_N_PHI, "repeats": 3, "ms_per_call": {}}
    for name, f in calls.items():
        best = math.inf
        for _ in range(3):
            t0 = time.process_time()
            for q in qs:
                f(q)
            best = min(best, (time.process_time() - t0) / len(qs) * 1e3)
        out["ms_per_call"][name] = round(best, 3)
    full = out["ms_per_call"]["full"]
    out["ratio_to_full"] = {k: round(v / full, 2) for k, v in out["ms_per_call"].items()}
    out["integrations_per_call"] = {"full": 1, "a": 1, "b": 6}
    out["fit_level"] = ("not measured in Stage 2: A8e's pass condition is per fit (<= 1/10 of "
                        "28.5 CPU-h) and needs sampling")  # fmt: skip
    _log(f"a8e: {json.dumps(out['ms_per_call'])}")
    (CACHE / "a8e.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def summarize() -> None:
    a8a = summarize_a8a(load_a8a())
    a8e_path = CACHE / "a8e.json"
    a8e = json.loads(a8e_path.read_text()) if a8e_path.exists() else None
    payload = {"description": "M8 Stage 2 validation: A8a (linearized vs full likelihood) and "
                              "the per-call part of A8e; no sampling. docs/notes/M8.md",
               "a8a": a8a, "a8e": a8e}  # fmt: skip
    write_json(OUT_NAME, payload)
    print(json.dumps({"modes": a8a["modes"], "mode_chosen": a8a["mode_chosen"]}, indent=1))


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="fork.export.m8")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("a8a")
    p.add_argument("--n-workers", type=int, default=8)
    sub.add_parser("a8e")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "a8a":
        run_a8a(args.n_workers)
    elif args.cmd == "a8e":
        run_a8e()
    else:
        summarize()


if __name__ == "__main__":
    main()
