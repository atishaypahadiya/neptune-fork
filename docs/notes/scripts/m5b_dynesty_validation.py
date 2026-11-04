"""Validate dynesty against the two emcee chains on seed 0 (free a, e ≤ 0.3).

    uv run python docs/notes/scripts/m5b_dynesty_validation.py [--n-workers 6]

Two nested runs (rng seeds 1 and 2; nlive 1000, multi/rslice, Δln Z < 0.1) on seed 0 with the
shared prior (``fork.inference.priors``; Uranus box from the Uranus-only fit). The weighted
samples are split at a = 25 AU into mode B (a ≈ 22 AU, e ≈ 0.28) and mode A (a ≈ 30 AU). Each
mode is compared with the emcee chain that sat in it (chain A: started from the 24-start point
fit; chain B: started on the 343° branch): medians and 16/84 % points, differences in units of
the emcee σ (half the 68 % width). Also: mode masses and ln Z with dynesty's error and the
spread over 100 jittered runs; run-to-run agreement. Output: docs/notes/M5b_dynesty_validation.json.
Nested results are cached in data/cache/nested/ (gitignored); re-runs reuse them.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from fork.data.synthetic import make_synthetic_observations
from fork.inference.fit import uranus_theta_to_x
from fork.inference.mcmc import circular_interval_deg
from fork.inference.nested import (
    NestedResult,
    a_au_of,
    mode_masses,
    run_nested,
)
from fork.inference.priors import UranusBox, free_a_space
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

EPOCH = 2371628.5
ROOT = Path(__file__).resolve().parents[3]
CHAINS = ROOT / "data" / "cache" / "chains"
CACHE = ROOT / "data" / "cache" / "nested"
OUT = ROOT / "docs" / "notes" / "M5b_dynesty_validation.json"
RNG_SEEDS = (1, 2)
SPLIT_A_AU = 25.0
EMCEE = {"A": "seed_0_free_a_e0.3", "B": "check_seed_0_free_a_e0.3_init343"}
N_LON = 3000  # samples per mode on which the 1846 longitude is computed
TOL_SIGMA = 0.5


def physical(space, q: np.ndarray) -> dict[str, np.ndarray]:
    th = np.array([space.theta(v) for v in q])
    return {"u_a_au": q[:, 0], "u_e_cos_varpi": q[:, 1], "u_e_sin_varpi": q[:, 2],
            "u_lam0_deg": np.degrees(q[:, 3]), "mass_earth": 10.0 ** th[:, 6],
            "a_au": th[:, 7], "e": th[:, 8], "lam0_deg": np.degrees(th[:, 10])}  # fmt: skip


def pct(x: np.ndarray, angle: bool) -> list[float]:
    if angle:
        lo, med, hi = circular_interval_deg(x)
        return [lo, med, hi]
    return [float(v) for v in np.percentile(x, [16, 50, 84])]


def lon_deg(space, q: np.ndarray, ctx, geo, rng) -> np.ndarray:
    idx = rng.choice(len(q), size=min(N_LON, len(q)), replace=False)
    return np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo)) for v in q[idx]])


def rel_deg(v: float, ref: float) -> float:
    return (v - ref + 180.0) % 360.0 - 180.0


def compare(dyn: list[float], emc: list[float], angle: bool) -> dict:
    sig = 0.5 * ((rel_deg(emc[2], emc[1]) - rel_deg(emc[0], emc[1])) if angle else emc[2] - emc[0])
    d = [rel_deg(a, b) if angle else a - b for a, b in zip(dyn, emc, strict=True)]
    return {"dynesty": [round(v, 6) for v in dyn], "emcee": [round(v, 6) for v in emc],
            "diff_sigma_16_50_84": [round(v / sig, 3) for v in d],
            "median_within_tol": bool(abs(d[1] / sig) <= TOL_SIGMA)}  # fmt: skip


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-workers", type=int, default=6)
    args = ap.parse_args()
    ctx = ForwardContext.from_horizons(EPOCH)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    geo = GalleGeometry.from_horizons()
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    space = free_a_space(pq, 0.3, box)
    rng = np.random.default_rng(7)

    # emcee reference chains, per mode (flat Uranus prior; check they sit inside the box)
    emc = {}
    for mode, name in EMCEE.items():
        q = np.load(CHAINS / f"{name}.npz")["samples"]
        inside = float(np.mean([box.contains(v[:4]) for v in q]))
        ph = physical(space, q)
        ph["lon_1846_deg"] = lon_deg(space, q, ctx, geo, rng)
        a = a_au_of(space, q)
        emc[mode] = {"ph": ph, "n": len(q), "inside_box": inside,
                     "frac_a_lt_split": float(np.mean(a < SPLIT_A_AU))}  # fmt: skip

    out = {"dataset": "seed 0", "space": "free a, e ≤ 0.3", "split_a_au": SPLIT_A_AU,
           "tol_sigma": TOL_SIGMA,
           "uranus_box": {"centre": box.centre, "half_width": box.half_width,
                          "sigma_scaled": box.sigma_scaled, "chi2_per_dof": box.chi2_per_dof,
                          "n_sigma": 50, "log_volume": box.log_volume},
           "emcee": {m: {"chain": EMCEE[m], "n_samples": v["n"], "inside_box": v["inside_box"],
                         "frac_a_lt_split": v["frac_a_lt_split"]} for m, v in emc.items()},
           "runs": {}}  # fmt: skip
    for seed in RNG_SEEDS:
        path = CACHE / f"seed_0_free_a_e0.3_rng{seed}.pkl"
        if path.exists():
            res = NestedResult.load(path)
        else:
            res = run_nested(obs, ctx, space, np.random.default_rng(seed),
                             n_workers=args.n_workers)  # fmt: skip
            res = NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s, seed)
            res.save(path)
        a = a_au_of(space, res.samples)
        masses = mode_masses(res, lambda s: a_au_of(space, s) < SPLIT_A_AU,
                             np.random.default_rng(100 + seed))  # fmt: skip
        eq = res.equal_weight(np.random.default_rng(200 + seed))
        a_eq = a_au_of(space, eq)
        run = {"wall_s": round(res.wall_s, 1), "niter": res.niter, "ncall": res.ncall,
               "logz": round(res.logz, 4), "logzerr": round(res.logzerr, 4),
               "logz_jitter_std": round(masses["logz_jitter_std"], 4),
               "chi2_min": round(-2.0 * res.max_logl, 4),
               "mode_mass_B_a_lt_25": round(masses["mass"], 4),
               "mode_mass_B_jitter_std": round(masses["mass_jitter_std"], 4),
               "n_equal_weight": len(eq), "modes": {}}  # fmt: skip
        for mode, sel in (("A", a_eq >= SPLIT_A_AU), ("B", a_eq < SPLIT_A_AU)):
            q = eq[sel]
            rec = {"n_equal_weight": int(sel.sum())}
            if sel.sum() >= 50:
                ph = physical(space, q)
                ph["lon_1846_deg"] = lon_deg(space, q, ctx, geo, rng)
                for k, x in ph.items():
                    angle = k.endswith("_deg")
                    rec[k] = compare(pct(x, angle), pct(emc[mode]["ph"][k], angle), angle)
                rec["all_medians_within_tol"] = all(v["median_within_tol"] for v in rec.values()
                                                    if isinstance(v, dict))  # fmt: skip
            run["modes"][mode] = rec
        # weighted (not resampled) mode-B mass as a cross-check of the equal-weight split
        run["mode_mass_B_weighted_check"] = round(float(res.weights[a < SPLIT_A_AU].sum()), 4)
        out["runs"][f"rng{seed}"] = run
        print(f"rng {seed}: ln Z {res.logz:.3f} ± {res.logzerr:.3f} (jitter "
              f"{masses['logz_jitter_std']:.3f}), mode B mass {masses['mass']:.4f} ± "
              f"{masses['mass_jitter_std']:.4f}, chi2_min {-2 * res.max_logl:.3f}, "
              f"{res.niter} it, {res.ncall} calls, {res.wall_s:.0f} s", flush=True)
        for mode in ("A", "B"):
            r = run["modes"][mode]
            print(f"  mode {mode}: n {r['n_equal_weight']}, medians within {TOL_SIGMA}σ: "
                  f"{r.get('all_medians_within_tol')}; " + ", ".join(
                      f"{k} {v['diff_sigma_16_50_84'][1]:+.2f}σ" for k, v in r.items()
                      if isinstance(v, dict)), flush=True)  # fmt: skip
    if len(out["runs"]) == 2:
        r1, r2 = (out["runs"][f"rng{s}"] for s in RNG_SEEDS)
        err = math.hypot(r1["logzerr"], r2["logzerr"])
        out["run_to_run"] = {
            "delta_logz": round(r2["logz"] - r1["logz"], 4), "combined_logzerr": round(err, 4),
            "logz_within_error": bool(abs(r2["logz"] - r1["logz"]) <= err),
            "delta_mode_mass_B": round(r2["mode_mass_B_a_lt_25"] - r1["mode_mass_B_a_lt_25"], 4),
        }  # fmt: skip
        print("run to run:", out["run_to_run"])
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
