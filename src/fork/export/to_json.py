"""Export site data to ``site/data/*.json``.

    uv run python -m fork.export.to_json

Currently exported:
- ``wh_convergence.json``: Uranus heliocentric-longitude error vs step size for both
  Wisdom–Holman compositions, 1781-03-13 → 1846-12-31, against an RK4 dt = 1 d reference
  (4th order; its own error is ~1e-10″). Feeds a future "integrator accuracy" panel.
- ``m5a_fits.json``: residuals after the M5a point fits on synthetic data (seed 0): no
  perturber, Le Verrier mode (a = 38.4 AU) and free a, perturber e ≤ 0.3 (the M5 prior);
  best-fit elements, local minima and the 1846-09-23 geocentric longitude vs true Neptune.
  Takes ~30 s (multistart LM).
- ``m5a_profile.json``: profile likelihood of the perturber's 1846-09-23 geocentric longitude
  (300°–350°, 1° steps, plus the true value), free a and Le Verrier mode, for noise seeds 0–3
  and noise-free data, with continuation sweeps and an e-bound log. Slow (about an hour: the ten
  curves run in parallel, one per process); regenerate alone with ``--only m5a_profile.json``.

M7 adds sightlines, residuals, sky map, degeneracy and energy exports here.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from fork.inference.fit import E_PRIOR_MAX
from fork.integrators import DEFAULT_DT_DAYS, rk4, wisdom_holman
from fork.integrators.wisdom_holman import SCHEME_DEFAULT, SCHEMES
from fork.model import bodies

SITE_DATA = Path(__file__).resolve().parents[3] / "site" / "data"
RAD_TO_ARCSEC = math.degrees(1.0) * 3600.0

A1_EPOCH_JD_TDB = 2371628.5  # 1781-03-13.0 TDB (Herschel's discovery; A1 initial conditions)
JD_1846_END_TDB = 2395296.5  # 1846-12-31.0 TDB
A1_BODIES = ("jupiter", "saturn", "uranus")
# Every dt divides 280 d, so the output grid below is on every integrator's step grid.
CONVERGENCE_DT_DAYS = (56.0, 40.0, 28.0, 20.0, 14.0, 12.0, 10.0, 7.0, 5.0, 3.5, 2.5)
CONVERGENCE_OUTPUT_STEP_DAYS = 280.0
REFERENCE_DT_DAYS = 1.0

SCHEME_LABELS = {
    "wh_drift_outer": "Wisdom–Holman, Kepler drift outermost (default; = REBOUND WHFast)",
    "wh_kick_outer": "Wisdom–Holman, kick outermost",
}


def _sig(x: float, digits: int = 4) -> float:
    return float(f"{x:.{digits}g}")


def a1_initial_conditions() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Heliocentric ecl_j2000 states of Jupiter, Saturn, Uranus at the A1 epoch (Horizons)."""
    from fork.data.horizons import state_vectors_ecl_j2000

    st = np.array([state_vectors_ecl_j2000(n, [A1_EPOCH_JD_TDB])[0] for n in A1_BODIES])
    masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in A1_BODIES])
    return st[:, :3].copy(), st[:, 3:].copy(), masses


def wh_convergence_data(
    x0: np.ndarray,
    v0: np.ndarray,
    masses: np.ndarray,
    dt_days: tuple[float, ...] = CONVERGENCE_DT_DAYS,
) -> dict:
    """Max |Δλ_Uranus| [″] and max |Δr_Uranus| [AU] over 1781–1846, per composition and dt."""
    t_out = np.arange(0.0, JD_1846_END_TDB - A1_EPOCH_JD_TDB, CONVERGENCE_OUTPUT_STEP_DAYS)
    ref = rk4.integrate(x0, v0, masses, t_out, REFERENCE_DT_DAYS)[:, 2]
    lon_ref_ecl_j2000_rad = np.unwrap(np.arctan2(ref[:, 1], ref[:, 0]))
    series = []
    for name, scheme in SCHEMES.items():
        dlon_arcsec, dr_au = [], []
        for dt in dt_days:
            pos = wisdom_holman.integrate(x0, v0, masses, t_out, dt, scheme)[:, 2]
            lon_ecl_j2000_rad = np.unwrap(np.arctan2(pos[:, 1], pos[:, 0]))
            dlon_arcsec.append(_sig(np.max(np.abs(lon_ecl_j2000_rad - lon_ref_ecl_j2000_rad))
                                    * RAD_TO_ARCSEC))  # fmt: skip
            dr_au.append(_sig(np.max(np.linalg.norm(pos - ref, axis=1))))
        series.append(
            {
                "scheme": name,
                "label": SCHEME_LABELS[name],
                "is_default": scheme == SCHEME_DEFAULT,
                "max_abs_dlon_uranus_arcsec": dlon_arcsec,
                "max_abs_dr_uranus_au": dr_au,
            }
        )
    return {
        "description": (
            "Uranus heliocentric ecliptic-J2000 longitude error vs Wisdom–Holman step size, "
            "Sun+Jupiter+Saturn+Uranus from JPL Horizons at 1781-03-13, integrated to "
            "1846-12-31. Error = max over outputs every 280 d of |WH − RK4(dt = 1 d)|."
        ),
        "span_jd_tdb": [A1_EPOCH_JD_TDB, JD_1846_END_TDB],
        "reference": {"integrator": "rk4", "dt_days": REFERENCE_DT_DAYS},
        "default_dt_days": DEFAULT_DT_DAYS,
        "observational_sigma_arcsec": 3.0,
        "dt_days": list(dt_days),
        "series": series,
    }


M5A_SEED = 0
M5A_YEARS = (1690, 1845)  # residual-curve grid, Jan 1.0 TDB


def _r(x, nd: int):
    return [round(float(v), nd) for v in np.atleast_1d(x)]


def _dlon_deg(a_rad: float, b_rad: float) -> float:
    return math.degrees((a_rad - b_rad + math.pi) % (2 * math.pi) - math.pi)


def m5a_fits_data(n_workers: int = 1) -> dict:
    """M5a site export. Deterministic given ``M5A_SEED`` (serial and pooled runs agree)."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about_bary
    from fork.data.synthetic import jd_tdb_jan1, make_synthetic_observations
    from fork.data.synthetic import true_uranus_lon_ecl_j2000_rad as true_lon
    from fork.inference.fit import fit_perturber_multistart, wrap_rad
    from fork.model.forward import (
        ForwardContext,
        perturber_theta_from_state,
        predict_uranus_lon,
        uranus_inc_raan_at_epoch,
    )
    from fork.model.hidden_mass import (
        BODE_A_AU,
        galle_jd_tdb,
        perturber_geocentric_longitude,
        true_neptune_geocentric_longitude,
    )

    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    obs, _truth = make_synthetic_observations(np.random.default_rng(M5A_SEED))
    post = obs.source == "synthetic_post_discovery"
    years = np.arange(M5A_YEARS[0], M5A_YEARS[1] + 1)
    grid = jd_tdb_jan1(years)
    truth_grid = true_lon(grid)
    nep = true_neptune_geocentric_longitude()

    def curve(theta):
        return _r(wrap_rad(truth_grid - predict_uranus_lon(theta, grid, ctx)) * RAD_TO_ARCSEC, 4)

    runs = {
        "le_verrier": fit_perturber_multistart(obs, ctx, inc_raan, a_fixed_au=BODE_A_AU,
                                               n_workers=n_workers),
        "free_a": fit_perturber_multistart(obs, ctx, inc_raan, n_workers=n_workers),
    }  # fmt: skip
    uranus_only = runs["free_a"].uranus_only
    curves = [{"id": "no_perturber", "label": "No perturber (A2)",
               "residual_arcsec": curve(uranus_only.theta)}]  # fmt: skip
    obs_resid = {"no_perturber": _r(uranus_only.residual_rad * RAD_TO_ARCSEC, 4)}
    fits = {}
    labels = {"le_verrier": "Le Verrier mode, a = 38.4 AU", "free_a": "Free a"}
    for mode, run in runs.items():
        best = run.best
        curves.append({"id": mode, "label": labels[mode], "residual_arcsec": curve(best.theta)})
        obs_resid[mode] = _r(best.residual_rad * RAD_TO_ARCSEC, 4)
        minima = []
        for m in run.minima:
            g = perturber_geocentric_longitude(m.best.theta, ctx)
            minima.append({
                "chi2": round(m.best.chi2, 4), "chi2_per_dof": round(m.best.chi2 / m.best.dof, 4),
                "mass_earth": round(m.best.mass_earth, 4), "a_au": round(float(m.best.theta[7]), 4),
                "e": round(float(m.best.theta[8]), 5),
                "lam0_ecl_j2000_deg": round(math.degrees(m.best.theta[10]), 4),
                "in_m5_prior": m.best.in_m5_prior, "n_starts": m.n_starts,
                "galle_dlon_vs_neptune_deg": round(_dlon_deg(g.lon_ecl_of_date_rad,
                                                             nep.lon_ecl_of_date_rad), 4),
            })  # fmt: skip
        g = perturber_geocentric_longitude(best.theta, ctx)
        fits[mode] = {
            "theta": _r(best.theta, 12),
            "mass_earth": round(best.mass_earth, 4), "a_au": round(float(best.theta[7]), 4),
            "e": round(float(best.theta[8]), 5), "chi2": round(best.chi2, 4), "dof": best.dof,
            "chi2_per_dof": round(best.chi2 / best.dof, 4),
            "rms_post_1781_arcsec": round(float(np.sqrt(np.mean(best.residual_rad[post] ** 2))
                                                * RAD_TO_ARCSEC), 4),
            "n_starts": len(run.fits), "status_counts": run.status_counts(),
            "n_distinct_minima": len(run.minima), "minima": minima,
            "galle": {"lon_ecl_of_date_deg": round(math.degrees(g.lon_ecl_of_date_rad), 4),
                      "lon_ecl_j2000_deg": round(math.degrees(g.lon_ecl_j2000_rad), 4),
                      "dlon_vs_neptune_deg": round(_dlon_deg(g.lon_ecl_of_date_rad,
                                                             nep.lon_ecl_of_date_rad), 4)},
        }  # fmt: skip
    th_nep = perturber_theta_from_state(about_bary("neptune", [A1_EPOCH_JD_TDB])[0],
                                        bodies.NEPTUNE.mass_msun)  # fmt: skip
    return {
        "description": "M5a point fits on synthetic data: Horizons Uranus (incl. real Neptune) + "
        "Gaussian noise; residual = noise-free truth − best fit, yearly Jan 1.0 TDB.",
        "seed": M5A_SEED, "dt_days": DEFAULT_DT_DAYS, "epoch_jd_tdb": A1_EPOCH_JD_TDB,
        "perturber_e_max": E_PRIOR_MAX,
        "galle_jd_tdb": galle_jd_tdb(),
        "years": years.tolist(),
        "curves": curves,
        "observations": {
            "year": _r(2000.0 + (obs.jd_tdb - 2451545.0) / 365.25, 4),
            "sigma_arcsec": _r(obs.sigma_rad * RAD_TO_ARCSEC, 4),
            "residual_arcsec": obs_resid,
        },
        "fits": fits,
        "true_neptune": {
            "mass_earth": round(10 ** th_nep[0], 4), "a_au_osculating_1781": round(th_nep[1], 4),
            "e_osculating_1781": round(th_nep[2], 5),
            "galle": {"lon_ecl_of_date_deg": round(math.degrees(nep.lon_ecl_of_date_rad), 4),
                      "lon_ecl_j2000_deg": round(math.degrees(nep.lon_ecl_j2000_rad), 4)},
        },
    }  # fmt: skip


def _profile_minimum(lon_deg: np.ndarray, chi2: np.ndarray) -> dict:
    """Lowest point of a sorted profile; if interior, refined by the parabola through it and its
    two neighbours. ``at_grid_edge`` flags a minimum that may lie outside the grid."""
    k = int(np.argmin(chi2))
    out = {"lon_ecl_of_date_deg": round(float(lon_deg[k]), 4), "chi2": round(float(chi2[k]), 4),
           "at_grid_edge": k in (0, lon_deg.size - 1)}  # fmt: skip
    if not out["at_grid_edge"]:
        c = np.polyfit(lon_deg[k - 1 : k + 2] - lon_deg[k], chi2[k - 1 : k + 2], 2)
        if c[0] > 0:
            x = -c[1] / (2 * c[0])
            if abs(x) <= max(lon_deg[k + 1] - lon_deg[k], lon_deg[k] - lon_deg[k - 1]):
                out["lon_ecl_of_date_deg"] = round(float(lon_deg[k] + x), 4)
                out["chi2"] = round(float(np.polyval(c, x)), 4)
    return out


M5A_PROFILE_SEEDS: tuple[int | None, ...] = (0, 1, 2, 3, None)  # None: noise-free truth
M5A_PROFILE_GRID_DEG = (300.0, 350.0, 1.0)
M5A_PROFILE_N_STARTS = 8


def _profile_interval_summary(lon_deg, chi2, best_lon, best_chi2, i_true, delta_chi2_intervals):
    """Intervals, overall minimum and Δχ²(truth) for one profile curve (grid + truth point +
    the unconstrained best, sorted). Δχ² is relative to the lowest χ² anywhere on the curve."""
    chi2_min = min(best_chi2, float(chi2.min()))
    curve_lon = np.concatenate([lon_deg, [best_lon]])
    curve_chi2 = np.concatenate([chi2, [best_chi2]])
    order = np.argsort(curve_lon)
    intervals = delta_chi2_intervals(curve_lon[order], curve_chi2[order] - chi2_min)
    return {
        "chi2_min": chi2_min,
        "intervals": [[round(a, 3), round(b, 3)] for a, b in intervals],
        "overall": _profile_minimum(curve_lon[order], curve_chi2[order]),
        "dchi2_true": round(float(chi2[i_true] - chi2_min), 4),
    }


def _profile_curve(seed: int | None, mode: str, grid_deg: tuple[float, float, float]) -> dict:
    """One profile curve (dataset × mode), serially in this process (module level: pickled)."""
    from fork.data.synthetic import LongitudeObservations, make_synthetic_observations
    from fork.inference.fit import fit_perturber_multistart, le_verrier_space
    from fork.inference.profile import (
        GalleGeometry,
        delta_chi2_intervals,
        galle_lon_of_date_rad,
        profile_galle_longitude,
    )
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch
    from fork.model.hidden_mass import BODE_A_AU, true_neptune_geocentric_longitude

    a_fixed = BODE_A_AU if mode == "le_verrier" else None
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    geo = GalleGeometry.from_horizons()
    nep_rad = true_neptune_geocentric_longitude().lon_ecl_of_date_rad
    nep_deg = math.degrees(nep_rad)
    lo, hi, step = grid_deg
    grid = np.arange(lo, hi + 0.5 * step, step)
    targets = np.sort(np.radians(np.concatenate([grid, [nep_deg]])))
    i_true = int(np.argmin(np.abs(targets - nep_rad)))
    obs, truth = make_synthetic_observations(np.random.default_rng(0 if seed is None else seed))
    if seed is None:
        obs = LongitudeObservations(obs.jd_tdb, truth, obs.sigma_rad, obs.source)
    best = fit_perturber_multistart(obs, ctx, inc_raan, a_fixed_au=a_fixed).best
    best_lon = math.degrees(galle_lon_of_date_rad(best.theta, ctx, geo))
    prof = profile_galle_longitude(obs, ctx, inc_raan, le_verrier_space(a_fixed), targets,
                                   mode=mode, n_starts=M5A_PROFILE_N_STARTS, geo=geo)  # fmt: skip
    lon_deg = np.degrees(targets)
    after = _profile_interval_summary(lon_deg, prof.chi2, best_lon, best.chi2, i_true,
                                      delta_chi2_intervals)  # fmt: skip
    before = _profile_interval_summary(lon_deg, prof.chi2_before_sweeps, best_lon, best.chi2,
                                       i_true, delta_chi2_intervals)  # fmt: skip
    chi2, chi2_min, overall = prof.chi2, after["chi2_min"], after["overall"]
    th = np.array([f.theta for f in prof.fits])
    on_bound = th[:, 8] >= E_PRIOR_MAX - 1e-9
    on_grid = np.ones(targets.size, dtype=bool)
    on_grid[i_true] = False
    drop = prof.chi2_before_sweeps - chi2
    return {
        "global_min": {
            "lon_ecl_of_date_deg": round(best_lon, 4), "chi2": round(best.chi2, 4),
            "dof": best.dof, "dlon_vs_neptune_deg": round(_dlon_deg(
                math.radians(best_lon), nep_rad), 4),
            "mass_earth": round(best.mass_earth, 3), "a_au": round(float(best.theta[7]), 3),
            "e": round(float(best.theta[8]), 4),
        },
        "overall_min": {**overall, "found_by": "multistart" if best.chi2 <= chi2_min
                        else "profile",
                        "dlon_vs_neptune_deg": round(overall["lon_ecl_of_date_deg"] - nep_deg, 4)},
        "profile_min_chi2_on_grid": round(float(chi2[on_grid].min()), 4),
        "profile_argmin_lon_on_grid_deg": round(float(lon_deg[on_grid][
            np.argmin(chi2[on_grid])]), 4),
        "dchi2_le_1_intervals_deg": after["intervals"],
        "dchi2_at_true_neptune": after["dchi2_true"],
        "before_sweeps": {"dchi2_le_1_intervals_deg": before["intervals"],
                          "dchi2_at_true_neptune": before["dchi2_true"],
                          "overall_min_lon_ecl_of_date_deg":
                              before["overall"]["lon_ecl_of_date_deg"]},
        "sweep_improved": [{"lon_ecl_of_date_deg": round(float(lon_deg[i]), 4),
                            "chi2_drop": round(float(drop[i]), 4)} for i in prof.sweep_improved],
        "n_sweep_pairs": prof.n_sweep_pairs,
        "on_e_bound_lon_deg": _r(lon_deg[on_grid & on_bound], 4),
        "profile": {
            "lon_ecl_of_date_deg": _r(lon_deg[on_grid], 4),
            "chi2": _r(chi2[on_grid], 4),
            "dchi2": _r(chi2[on_grid] - chi2_min, 4),
            "chi2_before_sweeps": _r(prof.chi2_before_sweeps[on_grid], 4),
            "mass_earth": _r(10.0 ** th[on_grid, 6], 3),
            "a_au": _r(th[on_grid, 7], 3),
            "e": _r(th[on_grid, 8], 4),
            "on_e_bound": [bool(v) for v in on_bound[on_grid]],
        },
        "at_true_neptune": {"chi2": round(float(chi2[i_true]), 4),
                            "mass_earth": round(float(10.0 ** th[i_true, 6]), 3),
                            "a_au": round(float(th[i_true, 7]), 3),
                            "e": round(float(th[i_true, 8]), 4),
                            "on_e_bound": bool(on_bound[i_true])},
        "max_constraint_miss_arcsec": _sig(max(f.constraint_miss_arcsec for f in prof.fits), 2),
        "n_lm_runs": prof.n_runs, "n_continuation_passes": prof.n_continuation_passes,
    }  # fmt: skip


def m5a_profile_data(
    n_workers: int = 1,
    seeds: tuple[int | None, ...] = M5A_PROFILE_SEEDS,
    grid_deg: tuple[float, float, float] = M5A_PROFILE_GRID_DEG,
) -> dict:
    """Profile χ² vs the perturber's 1846-09-23 geocentric longitude (ecliptic of date).

    For each dataset and mode: the unconstrained 24-start best fit, the profile on the grid
    (independent starts + neighbour passes + left→right / right→left continuation sweeps), and
    one profile point at true Neptune's longitude. Δχ² is relative to the smallest χ² found
    anywhere. The sweeps are sequential, so curves (dataset × mode) run in parallel, one per
    process. The pre-sweep intervals are kept for comparison, as is a log of grid points whose
    best fit sits on the e = 0.3 prior bound (kinks there are genuine).
    """
    from fork.inference.fit import run_pool
    from fork.inference.profile import GalleGeometry
    from fork.model.hidden_mass import BODE_A_AU, true_neptune_geocentric_longitude

    modes = ("le_verrier", "free_a")
    jobs = [(seed, mode, grid_deg) for seed in seeds for mode in modes]
    results = run_pool(_profile_curve, jobs, min(n_workers, len(jobs)))
    datasets = []
    for seed in seeds:
        entry = {"seed": seed, "noise": seed is not None, "modes": {}}
        for mode in modes:
            m = results[jobs.index((seed, mode, grid_deg))]
            entry["modes"][mode] = m
            print(f"seed {seed} {mode}: min {m['overall_min']['lon_ecl_of_date_deg']:.2f}° "
                  f"({m['overall_min']['found_by']}), Δχ²(true) {m['dchi2_at_true_neptune']:.2f}"
                  f", intervals {m['dchi2_le_1_intervals_deg']} (before sweeps "
                  f"{m['before_sweeps']['dchi2_le_1_intervals_deg']}), sweep-improved "
                  f"{len(m['sweep_improved'])} points, on e bound at "
                  f"{m['on_e_bound_lon_deg']}", flush=True)  # fmt: skip
        datasets.append(entry)
    geo = GalleGeometry.from_horizons()
    nep_deg = math.degrees(true_neptune_geocentric_longitude().lon_ecl_of_date_rad)
    return {
        "description": "Profile likelihood of the perturber's light-time-corrected geocentric "
        "ecliptic-of-date longitude on 1846-09-23 (Galle). At each longitude every other free "
        "parameter (Uranus a, e, ϖ, λ₀; perturber m, [a], e ≤ 0.3, ϖ; λ₀ solved from the "
        "constraint) is re-optimised: multistart LM, neighbour continuation, then sequential "
        "left→right / right→left sweeps; dchi2 = chi2 − min over all fits for that dataset.",
        "epoch_jd_tdb": A1_EPOCH_JD_TDB, "galle_jd_tdb": geo.jd_tdb, "dt_days": DEFAULT_DT_DAYS,
        "grid_deg": list(grid_deg), "n_starts_per_point": M5A_PROFILE_N_STARTS,
        "perturber_e_max": E_PRIOR_MAX,
        "le_verrier_a_au": BODE_A_AU,
        "true_neptune_lon_ecl_of_date_deg": round(nep_deg, 4),
        "datasets": datasets,
    }  # fmt: skip


def write_json(name: str, payload: dict, directory: Path = SITE_DATA) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    return path


def main(argv: list[str] | None = None) -> None:
    import argparse
    import os
    import time

    ap = argparse.ArgumentParser(description="Regenerate site/data/*.json")
    ap.add_argument("--only", action="append", help="write only this file (repeatable)")
    args = ap.parse_args(argv)
    n_workers = min(12, os.cpu_count() or 1)
    jobs = {
        "wh_convergence.json": lambda: wh_convergence_data(*a1_initial_conditions()),
        "m5a_fits.json": lambda: m5a_fits_data(n_workers=n_workers),
        "m5a_profile.json": lambda: m5a_profile_data(n_workers=n_workers),
    }
    for name in args.only or list(jobs):
        t0 = time.time()
        path = write_json(name, jobs[name]())
        print(f"wrote {path} ({path.stat().st_size} bytes, {time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
