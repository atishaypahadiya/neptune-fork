"""Analytic secular rates used to validate the Phase 2 forces, and the measurement of a
perihelion advance from integrated states.

All rates are in radians per day; ``ARCSEC_PER_CENTURY`` converts to arcseconds per Julian
century.
"""

from __future__ import annotations

import math

import numpy as np

from fork.data.frames import state_to_elements
from fork.integrators import wh_extra as W
from fork.model.bodies import G

ARCSEC_PER_CENTURY: float = math.degrees(1.0) * 3600.0 * 36525.0  # (rad/day) → ″/century


def mean_motion(a: float, mu: float) -> float:
    return math.sqrt(mu / a**3)


def gr_rate(
    a: float, e: float, mu: float, *, gr: float = 1.0, beta: float = 1.0, gamma: float = 1.0
) -> float:
    """1PN perihelion advance: gr (2 + 2γ − β)/3 · 6πμ/(c² a (1 − e²)) per orbit."""
    per_orbit = gr * (2.0 + 2.0 * gamma - beta) / 3.0 * 6.0 * math.pi * mu
    per_orbit /= W.C_AU_PER_DAY**2 * a * (1.0 - e * e)
    return per_orbit * mean_motion(a, mu) / (2.0 * math.pi)


def j2_rates(
    a: float, e: float, inc: float, mu: float, j2: float, r_sun_au: float = W.R_SUN_AU
) -> tuple[float, float]:
    """Secular (dω/dt, dΩ/dt) from J2, inclination ``inc`` to the equator of the central body:
    dω/dt = (3/4) n J₂ (R/p)² (5 cos² i − 1), dΩ/dt = −(3/2) n J₂ (R/p)² cos i (e.g. Murray &
    Dermott 1999, eqs. 6.249–6.250, to first order in J₂)."""
    p = a * (1.0 - e * e)
    k = mean_motion(a, mu) * j2 * (r_sun_au / p) ** 2
    return 0.75 * k * (5.0 * math.cos(inc) ** 2 - 1.0), -1.5 * k * math.cos(inc)


def laplace_coefficient(s: float, j: int, alpha: float, n: int = 20000) -> float:
    """b_s^{(j)}(α) = (1/π) ∫₀^{2π} cos jψ (1 − 2α cos ψ + α²)^{−s} dψ (periodic trapezoid)."""
    psi = 2.0 * math.pi * (np.arange(n) + 0.5) / n
    f = np.cos(j * psi) / (1.0 - 2.0 * alpha * np.cos(psi) + alpha * alpha) ** s
    return float(2.0 * f.mean())


def laplace_lagrange_rate(a: float, mu: float, mass: float, radius: float) -> float:
    """Lowest-order (in e) secular perihelion rate from a coplanar interior mass at ``radius``
    (a point planet on a circle, or a ring): (n/4) (m/M) α b_{3/2}^{(1)}(α), α = radius/a < 1
    (Murray & Dermott 1999, eq. 7.9 with ᾱ = 1 for an interior perturber)."""
    alpha = radius / a
    if not alpha < 1.0:
        raise ValueError("the mass must be interior")
    b = laplace_coefficient(1.5, 1, alpha)
    return 0.25 * mean_motion(a, mu) * (mass * G / mu) * alpha * b


def ring_gauss_rate(
    a: float, e: float, mu: float, mass: float, radius: float, n: int = 4000
) -> float:
    """Orbit-averaged dϖ/dt of a coplanar orbit of any eccentricity in the ring's exact field,
    from Gauss's equation dω/dt = √(1 − e²)/(n a e) · [−R cos f + T (1 + r/p) sin f] averaged
    over the mean anomaly (first order in the ring's mass). The field is radial in the ring's
    plane, so T = 0."""
    extra = W.make_extra(ring_mass=mass, ring_radius_au=radius)
    ecc = 2.0 * math.pi * (np.arange(n) + 0.5) / n
    r = a * (1.0 - e * np.cos(ecc))
    cos_f = (np.cos(ecc) - e) / (1.0 - e * np.cos(ecc))
    acc = np.empty(3)
    radial = np.empty(n)
    for i in range(n):
        W.position_accel(np.array([r[i], 0.0, 0.0]), mu, extra, acc)
        radial[i] = acc[0]
    weight = (1.0 - e * np.cos(ecc)) / n  # dM = (1 − e cos E) dE
    mean = float(np.sum(-radial * cos_f * weight))
    return math.sqrt(1.0 - e * e) / (mean_motion(a, mu) * a * e) * mean


def perihelion_longitude(x: np.ndarray, v: np.ndarray, mu: float) -> np.ndarray:
    """ϖ = Ω + ω [rad], unwrapped, from heliocentric states (T, 3)."""
    out = np.empty(x.shape[0])
    for i in range(x.shape[0]):
        el = state_to_elements(x[i], v[i], mu)
        out[i] = el[3] + el[4]
    return np.unwrap(out)


def elements_series(x: np.ndarray, v: np.ndarray, mu: float) -> np.ndarray:
    """(T, 6) of (a, e, i, Ω, ω, M) with Ω and ω unwrapped."""
    out = np.array([state_to_elements(x[i], v[i], mu) for i in range(x.shape[0])])
    out[:, 3] = np.unwrap(out[:, 3])
    out[:, 4] = np.unwrap(out[:, 4])
    return out


def slope(t: np.ndarray, y: np.ndarray) -> float:
    """Least-squares slope of y against t."""
    return float(np.polyfit(t - t.mean(), y, 1)[0])
