"""Least-squares fits of forward models to longitude observations.

M4: ``fit_uranus_only`` — Uranus's 6 elements, no perturber (test A2).
M5 will add the multistart perturber fits here.

Uranus is fitted in non-singular variables

    x = (a, k, h, λ₀, p, q),  k + i·h = e·exp(iϖ),  p + i·q = tan(i/2)·exp(iΩ)

so circular and planar orbits are interior points. Per CLAUDE.md rule 2 the in-plane four
(a, k, h, λ₀) are always free and start from the data alone (mean motion from the post-1781
longitudes, e = 0); (p, q), i.e. i and Ω, are held at their epoch values.

Why hold i, Ω: from longitudes alone they enter only through the reduction to the ecliptic,
≈ −tan²(i/2)·sin 2u (~10″ for Uranus), which is quadratic in (p, q) and invariant under
Ω → Ω + π. Freed, they act as a spare second harmonic that absorbs part of a perturber's signal.
``free_inc_raan=True`` frees them anyway, as a diagnostic only; since the gradient in (p, q)
vanishes at i = 0, that mode starts from i = 1° at several Ω and keeps the best.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from fork.data.synthetic import LongitudeObservations
from fork.model.bodies import M_CENTRAL_MSUN, URANUS, G

TWO_PI = 2.0 * math.pi
_MU_URANUS = G * (M_CENTRAL_MSUN + URANUS.mass_msun)
DEFAULT_OMEGA_STARTS_DEG = (0.0, 45.0, 90.0, 135.0)  # Ω and Ω + 180° are equivalent
START_INC_DEG = 1.0


def wrap_rad(d: np.ndarray) -> np.ndarray:
    """Wrap-safe angle difference into [−π, π) (CLAUDE.md rule 6)."""
    return (np.asarray(d) + math.pi) % TWO_PI - math.pi


def uranus_x_to_theta(x: np.ndarray) -> np.ndarray:
    a, k, h, lam0, p, q = x
    e = math.hypot(k, h)
    varpi = math.atan2(h, k) % TWO_PI
    tan_half_i = math.hypot(p, q)
    inc = 2.0 * math.atan(tan_half_i)
    raan = math.atan2(q, p) % TWO_PI if tan_half_i > 0 else 0.0
    return np.array([a, e, inc, raan, varpi, lam0 % TWO_PI])


def uranus_theta_to_x(theta: np.ndarray) -> np.ndarray:
    a, e, inc, raan, varpi, lam0 = theta
    t = math.tan(0.5 * inc)
    return np.array([a, e * math.cos(varpi), e * math.sin(varpi), lam0, t * math.cos(raan),
                     t * math.sin(raan)])  # fmt: skip


@dataclass(frozen=True)
class FitResult:
    theta: np.ndarray
    residual_rad: np.ndarray  # wrap(obs − model), unweighted
    chi2: float
    dof: int
    success: bool
    message: str
    nfev: int
    start: int  # index of the winning start

    @property
    def rms_arcsec(self) -> float:
        return float(np.sqrt(np.mean(self.residual_rad**2)) * 206264.80624709636)


def initial_uranus_theta(obs: LongitudeObservations, epoch_jd_tdb: float) -> np.ndarray:
    """Data-only starting point: circular, planar orbit with the observed mean motion.

    Linear fit λ(t) = λ₀ + n (t − epoch) to the unwrapped post-1781 longitudes (σ = 3″ era),
    then a = (μ / n²)^(1/3) from Kepler's third law.
    """
    sel = obs.sigma_rad <= np.median(obs.sigma_rad)
    order = np.argsort(obs.jd_tdb[sel])
    t = obs.jd_tdb[sel][order] - epoch_jd_tdb
    lon = np.unwrap(obs.lon_ecl_j2000_rad[sel][order])
    n, lam0 = np.polyfit(t, lon, 1)
    a = (_MU_URANUS / n**2) ** (1.0 / 3.0)
    return np.array([a, 0.0, 0.0, 0.0, 0.0, lam0 % TWO_PI])


def fit_uranus_only(
    obs: LongitudeObservations,
    model,
    inc_raan_ecl_j2000_rad: tuple[float, float],
    *,
    free_inc_raan: bool = False,
    omega_starts_deg: tuple[float, ...] = DEFAULT_OMEGA_STARTS_DEG,
) -> FitResult:
    """Levenberg–Marquardt fit of Uranus's elements (no perturber).

    a, e, ϖ, λ₀ are free; (i, Ω) are held at ``inc_raan_ecl_j2000_rad`` (the epoch values, see
    ``fork.model.forward.uranus_inc_raan_at_epoch``). ``free_inc_raan=True`` fits all six from a
    data-only multistart — diagnostic only (CLAUDE.md rule 2).
    """
    epoch = model.context.epoch_jd_tdb
    theta0 = initial_uranus_theta(obs, epoch)
    free = np.ones(6, dtype=bool)
    if not free_inc_raan:
        free[4:] = False
    sigma = obs.sigma_rad

    def residuals(x_free: np.ndarray, x_full: np.ndarray) -> np.ndarray:
        x = x_full.copy()
        x[free] = x_free
        pred = model.predict(uranus_x_to_theta(x), obs.jd_tdb)
        return wrap_rad(obs.lon_ecl_j2000_rad - pred) / sigma

    starts = []
    if free_inc_raan:
        for om in omega_starts_deg:
            th = theta0.copy()
            th[2], th[3] = math.radians(START_INC_DEG), math.radians(om)
            starts.append(uranus_theta_to_x(th))
    else:
        th = theta0.copy()
        th[2], th[3] = inc_raan_ecl_j2000_rad
        starts.append(uranus_theta_to_x(th))

    best = None
    for k, x0 in enumerate(starts):
        sol = least_squares(residuals, x0[free], args=(x0,), method="lm", x_scale="jac",
                            xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=2000)  # fmt: skip
        x = x0.copy()
        x[free] = sol.x
        chi2 = float(np.sum(sol.fun**2))
        if best is None or chi2 < best[0]:
            best = (chi2, k, x, sol)
    chi2, k, x, sol = best
    theta = uranus_x_to_theta(x)
    resid = wrap_rad(obs.lon_ecl_j2000_rad - model.predict(theta, obs.jd_tdb))
    return FitResult(theta, resid, chi2, obs.jd_tdb.size - int(free.sum()), bool(sol.success),
                     str(sol.message), int(sol.nfev), k)  # fmt: skip
