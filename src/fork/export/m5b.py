"""M5b: emcee posteriors on synthetic data → ``site/data/m5b_posteriors.json`` (summary only).

    uv run python -m fork.export.m5b            # all datasets (hours; local only, not CI)
    uv run python -m fork.export.m5b --seeds 0 1  # a subset (merged into the existing JSON)

Per dataset (noise seeds 0–9, plus noise-free truth), three posteriors:
  - ``free_a_e0.3``: free a, e ∈ [0, 0.3] (the M5 prior);
  - ``free_a_e0.1``: free a, e ∈ [0, 0.1] (prior-sensitivity run);
  - ``le_verrier_shape``: a, e, ϖ held at Le Verrier's published 1846 values; m, λ₀ and Uranus
    sampled.
Raw chains go to ``data/cache/chains/`` (gitignored) and are reused on re-runs, so an interrupted
run resumes. Reported per posterior: median and 16/84 % for m, a, e and the derived 1846-09-23
geocentric longitude (light-time corrected, ecliptic of date), whether truth lies inside each
68 % interval, τ, steps, acceptance, and the best χ² (chain maximum, or the bounded multistart
point fit if lower). A3b's calibration (truth inside the 68 % longitude interval in 50–85 % of
noise seeds) is computed over the noisy seeds.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np

from fork.export.to_json import A1_EPOCH_JD_TDB, SITE_DATA, write_json

CHAIN_DIR = Path(__file__).resolve().parents[3] / "data" / "cache" / "chains"
OUT_NAME = "m5b_posteriors.json"
NOISE_SEEDS = tuple(range(10))
RUNS = ("free_a_e0.3", "free_a_e0.1", "le_verrier_shape")
N_WALKERS = 48
N_STEPS = 20000
MAX_STEPS = 60000
DERIVED_MAX_SAMPLES = 6000  # 1846 longitude is computed on at most this many thinned samples
A3B_RANGE = (0.50, 0.85)


def _pct(x: np.ndarray) -> list[float]:
    return [float(v) for v in np.percentile(x, [16.0, 50.0, 84.0])]


def _summary_1d(x: np.ndarray, truth: float | None, nd: int) -> dict:
    lo, med, hi = _pct(x)
    out = {"median": round(med, nd), "lo68": round(lo, nd), "hi68": round(hi, nd),
           "minus": round(med - lo, nd), "plus": round(hi - med, nd)}  # fmt: skip
    if truth is not None:
        out["truth"] = round(truth, nd)
        out["truth_in_68"] = bool(lo <= truth <= hi)
    return out


def _heliocentric_lon_of_date_deg(theta: np.ndarray, ctx, jd: float) -> float:
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.model.forward import predict_states_ecl_j2000

    pos, _ = predict_states_ecl_j2000(theta, [jd], ctx)
    r = precession_matrix_ecl_j2000_to_ecl_of_date(jd) @ pos[0, 3]
    return math.degrees(math.atan2(r[1], r[0])) % 360.0


def _dataset(seed: int | None):
    from fork.data.synthetic import LongitudeObservations, make_synthetic_observations

    obs, truth = make_synthetic_observations(np.random.default_rng(0 if seed is None else seed))
    if seed is None:
        obs = LongitudeObservations(obs.jd_tdb, truth, obs.sigma_rad, obs.source)
    return obs


def _spaces(pq):
    from fork.inference import fit, mcmc
    from fork.model.hidden_mass import LE_VERRIER_1846 as LV

    varpi = LV.varpi_ecl_j2000_rad()
    s = math.sqrt(LV.e)
    return {
        "free_a_e0.3": (mcmc.free_a_space(pq, 0.3), dict(e_max=0.3)),
        "free_a_e0.1": (mcmc.free_a_space(pq, 0.1), dict(e_max=0.1)),
        "le_verrier_shape": (
            mcmc.held_shape_space(pq, LV.a_au, LV.e, varpi),
            dict(space=fit.PerturberSpace((None, LV.a_au, s * math.cos(varpi),
                                           s * math.sin(varpi), None)),
                 mode="le_verrier_shape"),
        ),
    }  # fmt: skip


def run_dataset(seed: int | None, ctx, inc_raan, n_workers: int, truth_p: dict) -> dict:
    from multiprocessing import get_context

    from fork.inference import mcmc
    from fork.inference.fit import fit_perturber_multistart, uranus_theta_to_x
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
    from fork.model.hidden_mass import LE_VERRIER_1846 as LV

    obs = _dataset(seed)
    key = "noise_free" if seed is None else f"seed_{seed}"
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    geo = GalleGeometry.from_horizons()
    out = {"seed": seed, "noise": seed is not None, "runs": {}}
    for run_name, (space, fit_kw) in _spaces(pq).items():
        t0 = time.time()
        point = fit_perturber_multistart(obs, ctx, inc_raan, n_workers=n_workers, **fit_kw).best
        path = CHAIN_DIR / f"{key}_{run_name}.npz"
        if path.exists():
            z = np.load(path)
            samples, lp = z["samples"], z["log_prob"]
            meta = json.loads(str(z["meta"]))
        else:
            q0 = space.q_from_theta(point.theta)
            # per-dataset, per-run sampler seed: reproducible and independent of run order
            rng = np.random.default_rng([2026, -1 if seed is None else seed, RUNS.index(run_name)])
            with get_context("spawn").Pool(n_workers) as pool:
                ch = mcmc.run_chain(obs, ctx, space, q0, rng, n_walkers=N_WALKERS,
                                    n_steps=N_STEPS, max_steps=MAX_STEPS, pool=pool)  # fmt: skip
            samples, lp = ch.samples, ch.log_prob
            meta = {"tau": [round(float(t), 2) for t in ch.tau], "n_steps": ch.n_steps,
                    "n_walkers": ch.n_walkers, "burn": ch.burn, "thin": ch.thin,
                    "acceptance": round(ch.acceptance, 4), "converged_50tau": ch.converged,
                    "chi2_chain_min": round(-2.0 * ch.max_log_prob, 4),
                    "names": list(space.names)}  # fmt: skip
            CHAIN_DIR.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, samples=samples, log_prob=lp, meta=json.dumps(meta))
        idx = np.linspace(0, len(samples) - 1, min(len(samples), DERIVED_MAX_SAMPLES)).astype(int)
        thetas = np.array([space.theta(q) for q in samples[idx]])
        lon = np.array([math.degrees(galle_lon_of_date_rad(th, ctx, geo)) for th in thetas])
        lo, med, hi = mcmc.circular_interval_deg(lon)

        def rel(v: float, med: float = med) -> float:
            return (v - med + 180.0) % 360.0 - 180.0

        nep = truth_p["galle_lon_ecl_of_date_deg"]
        chi2_best = min(meta["chi2_chain_min"], point.chi2)
        rec = {
            "n_samples": len(samples), "n_derived": len(idx), **meta,
            "chi2_best": round(chi2_best, 4), "dof": point.dof,
            "point_fit": {"chi2": round(point.chi2, 4), "mass_earth": round(point.mass_earth, 3),
                          "a_au": round(float(point.theta[7]), 3),
                          "e": round(float(point.theta[8]), 4)},
            "mass_earth": _summary_1d(10.0 ** thetas[:, 6], truth_p["mass_earth"], 3),
            "galle_lon_ecl_of_date_deg": {
                "median": round(med, 3), "lo68": round(lo, 3), "hi68": round(hi, 3),
                "minus": round(-rel(lo), 3), "plus": round(rel(hi), 3),
                "width68": round(rel(hi) - rel(lo), 3), "truth": nep,
                "median_minus_truth": round(-rel(nep), 3),
                "truth_in_68": bool(rel(lo) <= rel(nep) <= rel(hi)),
            },
        }  # fmt: skip
        if space.free[1]:
            rec["a_au"] = _summary_1d(thetas[:, 7], truth_p["a_au"], 3)
            rec["e"] = _summary_1d(thetas[:, 8], truth_p["e"], 4)
        if run_name == "le_verrier_shape":
            sub = thetas[:: max(1, len(thetas) // 2000)]
            helio = np.array([_heliocentric_lon_of_date_deg(th, ctx, LV.epoch_jd) for th in sub])
            hl, hm, hh = mcmc.circular_interval_deg(helio)
            rec["heliocentric_true_lon_1847_jan1_deg"] = {
                "median": round(hm, 3), "lo68": round(hl, 3), "hi68": round(hh, 3),
                "le_verrier_published": round(LV.true_lon_ecl_of_date_deg, 4)}  # fmt: skip
        out["runs"][run_name] = rec
        print(
            f"{key} {run_name}: lon {med:.2f} (-{rec['galle_lon_ecl_of_date_deg']['minus']:.2f}"
            f" +{rec['galle_lon_ecl_of_date_deg']['plus']:.2f}), m {rec['mass_earth']['median']}"
            f", chi2_best {chi2_best:.2f}, tau_max {max(meta['tau']):.0f}, "
            f"steps {meta['n_steps']}, {time.time() - t0:.0f} s",
            flush=True,
        )
    fa = out["runs"]["free_a_e0.3"]["chi2_best"]
    lvs = out["runs"]["le_verrier_shape"]["chi2_best"]
    out["le_verrier_shape_vs_free_a"] = {
        "delta_chi2": round(lvs - fa, 4),
        "delta_ln_likelihood": round(-0.5 * (lvs - fa), 4),
    }
    return out


def calibration(datasets: dict, run_name: str) -> dict:
    noisy = [d for d in datasets.values() if d["noise"] and run_name in d["runs"]]
    hits = [d["runs"][run_name]["galle_lon_ecl_of_date_deg"]["truth_in_68"] for d in noisy]
    widths = [d["runs"][run_name]["galle_lon_ecl_of_date_deg"]["width68"] for d in noisy]
    frac = float(np.mean(hits)) if hits else float("nan")
    return {"n_seeds": len(hits), "n_truth_in_68": int(sum(hits)), "fraction": round(frac, 3),
            "pass_50_85": bool(len(hits) >= 10 and A3B_RANGE[0] <= frac <= A3B_RANGE[1]),
            "width68_deg_median": round(float(np.median(widths)), 3) if widths else None,
            "width68_deg_range": [round(min(widths), 3), round(max(widths), 3)]
            if widths else None}  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.model import bodies
    from fork.model.forward import (
        ForwardContext,
        perturber_theta_from_state,
        uranus_inc_raan_at_epoch,
    )
    from fork.model.hidden_mass import true_neptune_geocentric_longitude

    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--seeds",
        nargs="*",
        default=None,
        help="noise seeds and/or 'none' (noise-free); default 0–9 and none",
    )
    ap.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    seeds = [*NOISE_SEEDS, None] if args.seeds is None else \
        [None if s.lower() == "none" else int(s) for s in args.seeds]  # fmt: skip
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    th_nep = perturber_theta_from_state(about("neptune", [A1_EPOCH_JD_TDB])[0],
                                        bodies.NEPTUNE.mass_msun)  # fmt: skip
    truth_p = {"mass_earth": round(10 ** th_nep[0], 4), "a_au": round(th_nep[1], 4),
               "e": round(th_nep[2], 5), "galle_lon_ecl_of_date_deg": round(math.degrees(
                   true_neptune_geocentric_longitude().lon_ecl_of_date_rad), 4)}  # fmt: skip
    path = SITE_DATA / OUT_NAME
    existing = json.loads(path.read_text())["datasets"] if path.exists() else {}
    for seed in seeds:
        key = "noise_free" if seed is None else f"seed_{seed}"
        existing[key] = run_dataset(seed, ctx, inc_raan, args.n_workers, truth_p)
        payload = {
            "description": "M5b emcee posteriors on synthetic data (Horizons Uranus incl. real "
            "Neptune + noise). Intervals are 16/50/84 % marginal percentiles; the 1846 longitude "
            "is the perturber's light-time-corrected geocentric ecliptic-of-date longitude on "
            "1846-09-23. Raw chains are not committed (data/cache/chains/).",
            "priors": {"log10_m_earth": [0.0, 2.0], "a_au_log_uniform": [20.0, 60.0],
                       "e_uniform_max": {"free_a_e0.3": 0.3, "free_a_e0.1": 0.1},
                       "angles": "uniform"},
            "le_verrier_shape": "a = 36.15 AU, e = 0.10761, ϖ = 284°45′ (equinox 1847.0) held",
            "sampler": {"package": "emcee", "moves": "0.8 DEMove + 0.2 DESnookerMove",
                        "n_walkers": N_WALKERS, "min_steps": N_STEPS, "max_steps": MAX_STEPS,
                        "target_steps_per_tau": 50, "burn_tau": 10},
            "truth": truth_p,
            "datasets": dict(sorted(existing.items())),
            "a3b_calibration": {r: calibration(existing, r) for r in RUNS},
        }  # fmt: skip
        write_json(OUT_NAME, payload)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
