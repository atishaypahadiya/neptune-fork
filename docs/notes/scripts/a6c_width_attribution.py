"""Why are A6c's 68 % widths (median 5.1°) narrower than M5b's (9.1°)? A local (Fisher)
ablation at the true θ; no sampling, no data values (the Jacobian of the whitened residuals
does not depend on them).

    uv run python docs/notes/scripts/a6c_width_attribution.py

For each variant: Fisher F = JᵀJ of the whitened residual vector with respect to the free
parameters, and the delta-method σ of Neptune's 1846-09-23 geocentric ecliptic longitude of date,
σ_λ = √(gᵀ F⁻¹ g). Free parameters are Uranus (a, e, ϖ, λ₀) and the perturber
(log10 m, a, k = e cos ϖ, h = e sin ϖ, λ₀), plus (sin i cos Ω, sin i sin Ω) relative to Uranus's
plane when inclined. This is the curvature at the truth only. It ignores the priors and
their boundaries, which truncate the mass–distance direction, so these widths are larger than
the posterior's (M5b 20.5° against 9.1°; A6c 7.6° against 5.1°). Use the ratios between
variants, not the absolute numbers.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
from fork.data.synthetic import make_synthetic_observations
from fork.export.to_json import A1_EPOCH_JD_TDB as E
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model import bodies
from fork.model import observation as om
from fork.model.forward import (
    ForwardContext,
    inclined_perturber_theta_from_state,
    perturber_theta_from_state,
    predict_uranus_lon,
    uranus_theta_from_state,
)

# steps for central differences, per free parameter
STEP_U = np.array([1e-6, 1e-7, 1e-7, 1e-8])  # a, e, ϖ, λ₀
STEP_P = np.array([1e-4, 1e-4, 1e-5, 1e-5, 1e-6])  # log10 m, a, k, h, λ₀
STEP_I = np.array([1e-5, 1e-5])  # sin i cos Ω, sin i sin Ω


def to_theta(z: np.ndarray, base: np.ndarray, inclined: bool) -> np.ndarray:
    th = base.copy()
    th[0], th[1], th[4], th[5] = z[0], z[1], z[2], z[3]
    lm, a, k, h, lam = z[4:9]
    th[6], th[7], th[8], th[9], th[10] = lm, a, math.hypot(k, h), math.atan2(h, k), lam
    if inclined:
        x, y = z[9], z[10]
        th[11], th[12] = math.asin(math.hypot(x, y)), math.atan2(y, x)
    return th


def to_z(th: np.ndarray, inclined: bool) -> np.ndarray:
    z = [th[0], th[1], th[4], th[5], th[6], th[7], th[8] * math.cos(th[9]),
         th[8] * math.sin(th[9]), th[10]]  # fmt: skip
    if inclined:
        z += [math.sin(th[11]) * math.cos(th[12]), math.sin(th[11]) * math.sin(th[12])]
    return np.array(z)


def sigma_lon_deg(resid_fn, base: np.ndarray, inclined: bool, ctx, geo, hold=()) -> float:
    z0 = to_z(base, inclined)
    steps = np.concatenate([STEP_U, STEP_P, STEP_I if inclined else []])
    free = [k for k in range(z0.size) if k not in hold]
    cols, g = [], []
    for k in free:
        dz = np.zeros_like(z0)
        dz[k] = steps[k]
        tp, tm = to_theta(z0 + dz, base, inclined), to_theta(z0 - dz, base, inclined)
        cols.append((resid_fn(tp) - resid_fn(tm)) / (2 * steps[k]))
        d = galle_lon_of_date_rad(tp, ctx, geo) - galle_lon_of_date_rad(tm, ctx, geo)
        g.append(((d + math.pi) % (2 * math.pi) - math.pi) / (2 * steps[k]))
    jac, g = np.array(cols).T, np.array(g)
    cov = np.linalg.pinv(jac.T @ jac, rcond=1e-14)
    return math.degrees(math.sqrt(g @ cov @ g))


def subset(obs: om.RaDecObservations, keep: np.ndarray) -> om.RaDecObservations:
    return dataclasses.replace(obs, epoch_index=obs.epoch_index[keep], is_dec=obs.is_dec[keep],
                               value_rad=obs.value_rad[keep],
                               group_index=obs.group_index[keep])  # fmt: skip


def main() -> None:
    ctx = ForwardContext.from_horizons(E)
    geo = GalleGeometry.from_horizons()
    th_u = uranus_theta_from_state(about("uranus", [E])[0])
    nep = about("neptune", [E])[0]
    planar = np.concatenate([th_u, perturber_theta_from_state(nep, bodies.NEPTUNE.mass_msun)])
    inclined = np.concatenate([th_u, inclined_perturber_theta_from_state(
        nep, bodies.NEPTUNE.mass_msun, (th_u[2], th_u[3]))])  # fmt: skip

    lon_obs, _ = make_synthetic_observations(np.random.default_rng(0))

    def m5b_resid(th):
        return predict_uranus_lon(th, lon_obs.jd_tdb, ctx) / lon_obs.sigma_rad

    st = om.structure_from_rows(om.read_rows(), om.ErrorModel.load())
    full = st.with_values(om.synthetic_values(st, om.true_uranus_ra_dec(st), None,
                                              noise=False)[0])  # fmt: skip
    ancient = np.array([om.GROUPS[g].startswith("ancient") for g in full.group_index])

    def m5b_sigma(o):  # M5b's noise level on the A6c layout: 12″ ancient, 3″ modern
        return dataclasses.replace(o, sigma_arcsec=np.array(
            [12.0 if g.startswith("ancient") else 3.0 for g in om.GROUPS]))  # fmt: skip

    def tau(o, t):
        return dataclasses.replace(o, offset_prior_arcsec=np.full(len(om.GROUPS), t))

    def radec(o):
        return lambda th: o.whitened(th, ctx)

    ra = subset(full, ~full.is_dec)

    def sig(o, anc, mod=None):  # ancient σ, modern σ (None: calibrated)
        cal = o.sigma_arcsec
        return dataclasses.replace(o, sigma_arcsec=np.array(
            [anc if g.startswith("ancient") else (cal[k] if mod is None else mod)
             for k, g in enumerate(om.GROUPS)]))  # fmt: skip

    cal_anc = full.sigma_arcsec[[om.GROUPS.index("ancient:ra"), om.GROUPS.index("ancient:dec")]]
    hold = (9, 10)
    path = [
        ("M5b", m5b_resid, planar, False, ()),
        ("1 Table 77 epochs: 279 RA, M5b σ 12″/3″, no offsets",
         radec(tau(sig(ra, 12.0, 3.0), 1e-6)), inclined, True, hold),
        ("2 + modern σ_g calibrated (1.4–2.9″)", radec(tau(sig(ra, 12.0), 1e-6)), inclined, True,
         hold),
        ("3 + ancient σ_g calibrated (4.7–5.1″)", radec(tau(ra, 1e-6)), inclined, True, hold),
        ("4 + group offsets, τ 5″/2″", radec(ra), inclined, True, hold),
        ("5 + Dec (250)", radec(full), inclined, True, hold),
        ("6 + i, Ω free = A6c", radec(full), inclined, True, ()),
    ]  # fmt: skip
    print("path M5b → A6c (each step adds one change):")
    prev = None
    for name, fn, base, inc, hd in path:
        s = sigma_lon_deg(fn, base, inc, ctx, geo, hd)
        ratio = "" if prev is None else f"{prev / s:6.2f}×"
        print(f"  {name:58s} {s:7.3f} {ratio}")
        prev = s
    print("calibrated ancient σ:", cal_anc)
    variants = [
        ("M5b: 87 heliocentric longitudes, 12″/3″, 9 params", m5b_resid, planar, False),
        ("A6c: 529 RA+Dec, calibrated σ_g, τ 5″/2″, 11 params", radec(full), inclined, True),
        ("  same, i and Ω held (9 params)", radec(full), inclined, True, (9, 10)),
        ("  RA only (279)", radec(ra), inclined, True),
        ("  RA only, i and Ω held", radec(ra), inclined, True, (9, 10)),
        ("  no offsets (τ → 0)", radec(tau(full, 1e-6)), inclined, True),
        ("  no offsets, i and Ω held", radec(tau(full, 1e-6)), inclined, True, (9, 10)),
        ("  offsets unconstrained (τ → ∞)", radec(tau(full, 1e6)), inclined, True),
        ("  M5b σ (12″/3″) per measurement", radec(m5b_sigma(full)), inclined, True),
        ("  RA only, M5b σ, no offsets, i Ω held", radec(tau(m5b_sigma(ra), 1e-6)), inclined,
         True, (9, 10)),
        ("  modern only (no ancient)", radec(subset(full, ~ancient)), inclined, True),
    ]  # fmt: skip
    print(f"{'variant':58s} {'σ_λ (°)':>8s} {'2σ (°)':>7s}")
    for name, fn, base, inc, *hold in variants:
        s = sigma_lon_deg(fn, base, inc, ctx, geo, hold[0] if hold else ())
        print(f"{name:58s} {s:8.3f} {2 * s:7.3f}")
    w = 1.0 / full.sigma_arcsec[full.group_index] ** 2
    m5b_w = 1 / (lon_obs.sigma_rad * om.ARCSEC_PER_RAD) ** 2
    print("Σ 1/σ² [arcsec⁻²]: M5b", round(float(m5b_w.sum()), 2),
          "| A6c all", round(float(w.sum()), 2), "| A6c RA", round(float(w[~full.is_dec].sum()), 2),
          "| A6c modern RA", round(float(w[~full.is_dec & ~ancient].sum()), 2))  # fmt: skip


if __name__ == "__main__":
    main()
