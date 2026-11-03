"""M5a diagnostics of the 1846 longitude miss → docs/notes/M5a_diagnostics.json.

    uv run python docs/notes/scripts/m5a_diagnostics.py

Datasets: noise-free Horizons truth, and synthetic noise seeds 0–4 (same schedule).
  1. Noise-free recovery: free-a and Le Verrier-mode 24-start multistart.
  3. Eccentricity sensitivity: free-a multistart with the perturber's e capped at 0.05 and 0.1
     (smooth cap, ``PerturberSpace.e_max``).
  4. Le Verrier's published 31 Aug 1846 solution (``LE_VERRIER_1846``: m, a, e, ϖ held; only
     its λ₀ and Uranus's a, e, ϖ, λ₀ free), with the alternative mass 1/9300 M☉ as a sensitivity
     case. Also compared: the fitted perturber's heliocentric mean and true longitude on 1847 Jan 1
     (ecliptic of date) against Le Verrier's own 318°47′ and 326°32′.
The profile likelihood (item 2) is ``python -m fork.export.to_json --only m5a_profile.json``.
"""

import json
import math
import os
import time
from pathlib import Path

import numpy as np

from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date, state_to_elements
from fork.data.synthetic import LongitudeObservations, make_synthetic_observations
from fork.inference.fit import (
    PerturberSpace,
    fit_perturber_multistart,
    held_orbit_space,
    wrap_rad,
)
from fork.model.bodies import M_CENTRAL_MSUN, M_EARTH_MSUN, G
from fork.model.forward import ForwardContext, predict_states_ecl_j2000, uranus_inc_raan_at_epoch
from fork.model.hidden_mass import (
    BODE_A_AU,
    LE_VERRIER_1846,
    LE_VERRIER_1846_MASS_ALT_MSUN,
    perturber_geocentric_longitude,
    true_neptune_geocentric_longitude,
)

EPOCH = 2371628.5
OUT = Path(__file__).resolve().parents[1] / "M5a_diagnostics.json"
ARC = math.degrees(1.0) * 3600.0
SEEDS = (None, 0, 1, 2, 3, 4)  # None: noise-free
N_WORKERS = min(12, os.cpu_count() or 1)


def lon_of_date_deg(lon_ecl_j2000_rad: float, jd: float) -> float:
    r = precession_matrix_ecl_j2000_to_ecl_of_date(jd)
    w = r @ np.array([math.cos(lon_ecl_j2000_rad), math.sin(lon_ecl_j2000_rad), 0.0])
    return math.degrees(math.atan2(w[1], w[0])) % 360.0


def heliocentric_1847(theta: np.ndarray, ctx: ForwardContext) -> dict:
    """Perturber's osculating mean and true heliocentric longitude at 1847 Jan 1.0, of date."""
    jd = LE_VERRIER_1846.epoch_jd
    pos, vel = predict_states_ecl_j2000(theta, [jd], ctx)
    r, v = pos[0, 3], vel[0, 3]
    m = 10.0 ** theta[6] * M_EARTH_MSUN
    _a, _e, _i, raan, argp, mean = state_to_elements(r, v, G * (M_CENTRAL_MSUN + m))
    true_lon = math.atan2(r[1], r[0])
    return {
        "mean_lon_ecl_of_date_deg": round(lon_of_date_deg(raan + argp + mean, jd), 3),
        "true_lon_ecl_of_date_deg": round(lon_of_date_deg(true_lon, jd), 3),
    }


def summary(run, obs, ctx, nep_rad: float, with_1847: bool = False) -> dict:
    b = run.best
    post = obs.source == "synthetic_post_discovery"
    g = perturber_geocentric_longitude(b.theta, ctx)
    out = {
        "mass_earth": round(b.mass_earth, 3), "a_au": round(float(b.theta[7]), 3),
        "e": round(float(b.theta[8]), 4), "varpi_ecl_j2000_deg": round(math.degrees(b.theta[9]), 2),
        "chi2": round(b.chi2, 3), "dof": b.dof, "chi2_per_dof": round(b.chi2 / b.dof, 4),
        "rms_post_1781_arcsec": round(float(np.sqrt(np.mean(b.residual_rad[post] ** 2)) * ARC), 3),
        "galle_lon_ecl_of_date_deg": round(math.degrees(g.lon_ecl_of_date_rad), 3),
        "dlon_vs_neptune_deg": round(math.degrees(wrap_rad(g.lon_ecl_of_date_rad - nep_rad)), 3),
        "n_distinct_minima": len(run.minima), "status_counts": run.status_counts(),
    }  # fmt: skip
    if with_1847:
        out["heliocentric_1847_jan1"] = heliocentric_1847(b.theta, ctx)
    return out


def main() -> None:
    t0 = time.time()
    ctx = ForwardContext.from_horizons(EPOCH)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    nep_rad = true_neptune_geocentric_longitude().lon_ecl_of_date_rad
    lv = LE_VERRIER_1846
    varpi = lv.varpi_ecl_j2000_rad()
    log10_m_alt = math.log10(LE_VERRIER_1846_MASS_ALT_MSUN / M_EARTH_MSUN)
    held = held_orbit_space(lv.log10_m_earth, lv.a_au, lv.e, varpi)
    lv_spaces = {
        "le_verrier_1846": held,
        "le_verrier_1846_m_1_9300": held_orbit_space(log10_m_alt, lv.a_au, lv.e, varpi),
        # context only (not item 4's spec): his a, e, ϖ held, mass free as well as λ₀
        "le_verrier_1846_m_free": PerturberSpace((None, *held.held[1:])),
    }
    results = {}
    for seed in SEEDS:
        obs, truth = make_synthetic_observations(np.random.default_rng(0 if seed is None else seed))
        if seed is None:
            obs = LongitudeObservations(obs.jd_tdb, truth, obs.sigma_rad, obs.source)
        key = "noise_free" if seed is None else f"seed_{seed}"
        res = {}
        runs = {
            "free_a": dict(),
            "le_verrier": dict(a_fixed_au=BODE_A_AU),
            "free_a_e_max_0.1": dict(e_max=0.1),
            "free_a_e_max_0.05": dict(e_max=0.05),
        }
        for name, kw in runs.items():
            run = fit_perturber_multistart(obs, ctx, inc_raan, n_workers=N_WORKERS, **kw)
            res[name] = summary(run, obs, ctx, nep_rad)
        for name, space in lv_spaces.items():
            run = fit_perturber_multistart(obs, ctx, inc_raan, space=space, mode=name,
                                           n_workers=N_WORKERS)  # fmt: skip
            res[name] = summary(run, obs, ctx, nep_rad, with_1847=True)
        results[key] = res
        print(key, {k: (v["dlon_vs_neptune_deg"], v["chi2"]) for k, v in res.items()}, flush=True)
    payload = {
        "description": "M5a diagnostics of the 1846-09-23 longitude miss; best of 24 λ₀ starts "
        "per fit. dlon = perturber − true Neptune, geocentric, light-time corrected, ecliptic "
        "of date.",
        "true_neptune_lon_ecl_of_date_deg": round(math.degrees(nep_rad), 4),
        "le_verrier_1846": {
            "a_au": lv.a_au,
            "e": lv.e,
            "mass_msun": lv.mass_msun,
            "mass_alt_msun": LE_VERRIER_1846_MASS_ALT_MSUN,
            "mass_earth": round(10**lv.log10_m_earth, 3),
            "varpi_ecl_of_date_deg": round(lv.varpi_ecl_of_date_deg, 4),
            "varpi_ecl_j2000_deg": round(math.degrees(lv.varpi_ecl_j2000_rad()), 4),
            "published_1847_jan1": {
                "mean_lon_ecl_of_date_deg": round(lv.mean_lon_ecl_of_date_deg, 4),
                "true_lon_ecl_of_date_deg": round(lv.true_lon_ecl_of_date_deg, 4),
            },
            "sources": list(lv.sources),
        },
        "results": results,
    }
    OUT.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {OUT} in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
