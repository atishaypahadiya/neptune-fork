"""Export site data to ``site/data/*.json``.

    uv run python -m fork.export.to_json

Currently exported:
- ``wh_convergence.json``: Uranus heliocentric-longitude error vs step size for both
  Wisdom–Holman compositions, 1781-03-13 → 1846-12-31, against an RK4 dt = 1 d reference
  (4th order; its own error is ~1e-10″). Feeds a future "integrator accuracy" panel.
- ``m5a_fits.json``: residuals after the M5a point fits on synthetic data (seed 0): no
  perturber, Le Verrier mode (a = 38.4 AU) and free a; best-fit elements, local minima and the
  1846-09-23 geocentric longitude vs true Neptune. Takes ~30 s (multistart LM).

M7 adds sightlines, residuals, sky map, degeneracy and energy exports here.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

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


def write_json(name: str, payload: dict, directory: Path = SITE_DATA) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    return path


def main() -> None:
    import os

    x0, v0, masses = a1_initial_conditions()
    path = write_json("wh_convergence.json", wh_convergence_data(x0, v0, masses))
    print(f"wrote {path} ({path.stat().st_size} bytes)")
    path = write_json("m5a_fits.json", m5a_fits_data(n_workers=min(12, os.cpu_count() or 1)))
    print(f"wrote {path} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
