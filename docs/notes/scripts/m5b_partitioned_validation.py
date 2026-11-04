"""Validate partitioned nested sampling on seed 0 (free a, e ≤ 0.3).

    uv run python docs/notes/scripts/m5b_partitioned_validation.py --run K RNG [--n-workers N]
    uv run python docs/notes/scripts/m5b_partitioned_validation.py --summarize

Slices in a: [20, 25], [25, 28], [28, 60] AU (``nested.A_SLICES_AU``), each sampled
separately by dynesty (nlive 1000, multi/rslice, Δln Z < 0.1), rng seeds 1 and 2. ``--run``
does one (slice, rng) run. It logs progress to data/cache/nested/progress.log, checkpoints to
data/cache/nested/*.ckpt (an interrupted run resumes), and caches the result as a .pkl.
``--summarize`` combines Z = Σ w_k Z_k per rng seed and checks the gate:

1. per-slice ln Z_k and mass shares agree across rng seeds within the combined error;
2. combined ln Z agrees with the unpartitioned runs (−73.046 ± 0.222, −73.313 ± 0.222);
3. slice [20, 25] recovers mode B and matches emcee chain B, and slice [28, 60] matches emcee
   chain A restricted to a ≥ 28 (medians within 0.5σ of the emcee σ);
4. the 1846-09-23 longitude posterior (median, 68 %), combined and per slice.

Output: docs/notes/M5b_partitioned_validation.json.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from fork.data.synthetic import make_synthetic_observations
from fork.inference.fit import uranus_theta_to_x
from fork.inference.mcmc import circular_interval_deg
from fork.inference.nested import A_SLICES_AU, NestedResult, PartitionedResult, run_nested
from fork.inference.priors import UranusBox, free_a_space
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

EPOCH = 2371628.5
ROOT = Path(__file__).resolve().parents[3]
CHAINS = ROOT / "data" / "cache" / "chains"
CACHE = ROOT / "data" / "cache" / "nested"
OUT = ROOT / "docs" / "notes" / "M5b_partitioned_validation.json"
RNG_SEEDS = (1, 2)
UNPARTITIONED = {"rng1": (-73.046, 0.222), "rng2": (-73.313, 0.222)}
EMCEE = {0: "check_seed_0_free_a_e0.3_init343", 2: "seed_0_free_a_e0.3"}  # slice → chain
N_LON = 3000
TOL_SIGMA = 0.5


def setup():
    ctx = ForwardContext.from_horizons(EPOCH)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    return ctx, obs, free_a_space(pq, 0.3, box)


def run_path(k: int, rng: int, suffix: str) -> Path:
    return CACHE / f"seed_0_free_a_e0.3_slice{k}_rng{rng}{suffix}"


def do_run(k: int, rng: int, n_workers: int) -> None:
    ctx, obs, full = setup()
    lo, hi = A_SLICES_AU[k]
    space = full.with_a_slice(lo, hi)
    out = run_path(k, rng, ".pkl")
    if out.exists():
        print(f"cached: {out}")
        return
    res = run_nested(obs, ctx, space, np.random.default_rng([rng, k]), n_workers=n_workers,
                     progress_log=CACHE / "progress.log",
                     checkpoint_file=run_path(k, rng, ".ckpt"),
                     label=f"slice{k} a[{lo:g},{hi:g}] rng{rng}")  # fmt: skip
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s, rng).save(out)
    print(f"slice {k} rng {rng}: ln Z_k {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, "
          f"{res.wall_s / 60:.1f} min")  # fmt: skip


def physical(space, q: np.ndarray) -> dict[str, np.ndarray]:
    th = np.array([space.theta(v) for v in q])
    return {"u_a_au": q[:, 0], "u_e_cos_varpi": q[:, 1], "u_e_sin_varpi": q[:, 2],
            "u_lam0_deg": np.degrees(q[:, 3]), "mass_earth": 10.0 ** th[:, 6],
            "a_au": th[:, 7], "e": th[:, 8], "lam0_deg": np.degrees(th[:, 10])}  # fmt: skip


def pct(x: np.ndarray, angle: bool) -> list[float]:
    if angle:
        return list(circular_interval_deg(x))
    return [float(v) for v in np.percentile(x, [16, 50, 84])]


def rel_deg(v: float, ref: float) -> float:
    return (v - ref + 180.0) % 360.0 - 180.0


def lon_deg(space, q, ctx, geo, rng) -> np.ndarray:
    idx = rng.choice(len(q), size=min(N_LON, len(q)), replace=False)
    return np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo)) for v in q[idx]])


def lon_summary(x: np.ndarray) -> dict:
    lo, med, hi = circular_interval_deg(x)
    return {"median": round(med, 3), "lo68": round(lo, 3), "hi68": round(hi, 3),
            "width68": round(rel_deg(hi, med) - rel_deg(lo, med), 3)}  # fmt: skip


def compare(dyn: dict, emc: dict) -> dict:
    out, ok = {}, True
    for name in dyn:
        angle = name.endswith("_deg")
        d, e = pct(dyn[name], angle), pct(emc[name], angle)
        sig = 0.5 * ((rel_deg(e[2], e[1]) - rel_deg(e[0], e[1])) if angle else e[2] - e[0])
        diff = [(rel_deg(a, b) if angle else a - b) / sig for a, b in zip(d, e, strict=True)]
        width = ((rel_deg(d[2], d[1]) - rel_deg(d[0], d[1])) if angle else d[2] - d[0]) / (2 * sig)
        within = abs(diff[1]) <= TOL_SIGMA
        ok &= within
        out[name] = {"dynesty_16_50_84": [round(v, 6) for v in d],
                     "emcee_16_50_84": [round(v, 6) for v in e],
                     "diff_sigma_16_50_84": [round(v, 3) for v in diff],
                     "width_ratio": round(width, 3), "median_within_tol": bool(within)}  # fmt: skip
    out["all_medians_within_tol"] = bool(ok)
    return out


def summarize() -> None:
    ctx, _, full = setup()
    geo = GalleGeometry.from_horizons()
    rng = np.random.default_rng(11)
    parts = {r: PartitionedResult(tuple(NestedResult.load(run_path(k, r, ".pkl"))
                                        for k in range(len(A_SLICES_AU))))
             for r in RNG_SEEDS}  # fmt: skip
    emc = {}
    for k, name in EMCEE.items():
        q = np.load(CHAINS / f"{name}.npz")["samples"]
        a = np.exp(q[:, 5])
        lo, hi = A_SLICES_AU[k]
        inside = (a >= lo) & (a <= hi)
        ph = physical(full, q[inside])
        ph["lon_1846_deg"] = lon_deg(full, q[inside], ctx, geo, rng)
        emc[k] = {"ph": ph, "chain": name, "frac_in_slice": round(float(inside.mean()), 4)}
    out = {"dataset": "seed 0", "space": "free a, e ≤ 0.3", "slices_au": A_SLICES_AU,
           "slice_prior_weights": [round(full.with_a_slice(*s).slice_weight, 6)
                                   for s in A_SLICES_AU],
           "unpartitioned_logz": UNPARTITIONED, "tol_sigma": TOL_SIGMA,
           "emcee_reference": {k: {"chain": v["chain"], "frac_in_slice": v["frac_in_slice"]}
                               for k, v in emc.items()},
           "rng": {}}  # fmt: skip
    for r, p in parts.items():
        rec = {"logz": round(p.logz, 4), "logzerr": round(p.logzerr, 4),
               "shares": [round(float(v), 5) for v in p.shares],
               "share_errs": [round(float(v), 5) for v in p.share_errs], "slices": []}  # fmt: skip
        for k, res in enumerate(p.slices):
            eq = res.equal_weight(np.random.default_rng(100 * r + k))
            lon = lon_deg(res.space, eq, ctx, geo, rng)
            th_best = res.space.theta(res.samples[int(np.argmax(res.results["logl"]))])
            s = {"a_au": A_SLICES_AU[k], "logz_k": round(res.logz, 4),
                 "logzerr_k": round(res.logzerr, 4), "niter": res.niter, "ncall": res.ncall,
                 "wall_min": round(res.wall_s / 60, 1), "chi2_min": round(-2 * res.max_logl, 4),
                 "best": {"a_au": round(float(th_best[7]), 3), "e": round(float(th_best[8]), 4),
                          "mass_earth": round(float(10 ** th_best[6]), 3)},
                 "lon_1846_deg": lon_summary(lon)}  # fmt: skip
            if k in emc:
                ph = physical(res.space, eq)
                ph["lon_1846_deg"] = lon
                s["vs_emcee"] = compare(ph, emc[k]["ph"])
            rec["slices"].append(s)
        rec["lon_1846_deg_combined"] = lon_summary(lon_deg(full, p.equal_weight(
            np.random.default_rng(300 + r)), ctx, geo, rng))  # fmt: skip
        rec["vs_unpartitioned"] = {
            f"rng{u_r}": {"delta": round(p.logz - uz, 4),
                          "combined_err": round(math.hypot(p.logzerr, uerr), 4),
                          "delta_in_sigma": round((p.logz - uz) / math.hypot(p.logzerr, uerr), 3)}
            for u_r, (uz, uerr) in ((1, UNPARTITIONED["rng1"]), (2, UNPARTITIONED["rng2"]))
        }  # fmt: skip
        out["rng"][f"rng{r}"] = rec
    r1, r2 = (out["rng"][f"rng{r}"] for r in RNG_SEEDS)
    gate1 = []
    for k in range(len(A_SLICES_AU)):
        a, b = r1["slices"][k], r2["slices"][k]
        dz, ez = b["logz_k"] - a["logz_k"], math.hypot(a["logzerr_k"], b["logzerr_k"])
        ds, es = r2["shares"][k] - r1["shares"][k], math.hypot(r1["share_errs"][k],
                                                                r2["share_errs"][k])  # fmt: skip
        gate1.append({"slice": k, "delta_logz_k": round(dz, 4), "err": round(ez, 4),
                      "logz_k_within": bool(abs(dz) <= ez), "delta_share": round(ds, 5),
                      "share_err": round(es, 5), "share_within": bool(abs(ds) <= es)})  # fmt: skip
    out["gate"] = {
        "1_per_slice_across_rng": gate1,
        "2_combined_vs_unpartitioned": {k: v["vs_unpartitioned"] for k, v in out["rng"].items()},
        "3_slice0_vs_chainB": [out["rng"][f"rng{r}"]["slices"][0]["vs_emcee"][
            "all_medians_within_tol"] for r in RNG_SEEDS],
        "3_slice2_vs_chainA": [out["rng"][f"rng{r}"]["slices"][2]["vs_emcee"][
            "all_medians_within_tol"] for r in RNG_SEEDS],
    }  # fmt: skip
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["gate"], indent=1))
    for r, rec in out["rng"].items():
        print(r, "ln Z", rec["logz"], "±", rec["logzerr"], "shares", rec["shares"], "±",
              rec["share_errs"], "lon", rec["lon_1846_deg_combined"])  # fmt: skip
        for s in rec["slices"]:
            print(
                "  ",
                s["a_au"],
                "ln Z_k",
                s["logz_k"],
                "±",
                s["logzerr_k"],
                "chi2_min",
                s["chi2_min"],
                "best",
                s["best"],
                "lon",
                s["lon_1846_deg"],
                f"{s['wall_min']} min",
            )
    print(f"wrote {OUT}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", nargs=2, type=int, metavar=("SLICE", "RNG"))
    ap.add_argument("--summarize", action="store_true")
    ap.add_argument("--n-workers", type=int, default=4)
    args = ap.parse_args()
    if args.run:
        do_run(*args.run, args.n_workers)
    if args.summarize:
        summarize()


if __name__ == "__main__":
    main()
