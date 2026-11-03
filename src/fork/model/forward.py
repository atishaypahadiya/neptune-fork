"""Forward model: θ → Uranus heliocentric ecliptic-J2000 longitudes.

Generic interface (Phase 2 adds other models beside this one)::

    model.param_names            # names of θ's entries, in order
    model.predict(theta, t_jd)   # -> lon_ecl_j2000_rad, shape like t_jd

Bodies integrated: Sun (with inner planets folded in, centred on their common barycenter),
Jupiter, Saturn, Uranus, and optionally one planar perturber. Jupiter and Saturn start from their
Horizons states at the model epoch and are **not** fitted. Uranus's six osculating heliocentric
elements at the epoch are all free (project rule 2): nothing about modern Uranus enters
except through θ.

Uranus θ (6), all angles in ecl_j2000:
    a_au, e, inc_rad, raan_rad (Ω), varpi_rad (ϖ = Ω + ω), lam0_rad (mean longitude λ = ϖ + M)
Perturber θ (5, optional, planar: i = 0 in ecl_j2000):
    log10_m_earth, a_au, e, varpi_rad, lam0_rad

Elements are converted with μ = G (M_central + m_body), the body's own heliocentric two-body μ.
The integrator is ``wisdom_holman`` with its project defaults (dt = 14 d, REBOUND composition).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numba import njit
from numpy.typing import ArrayLike

from fork.data.frames import elements_to_state, state_to_elements
from fork.integrators import DEFAULT_DT_DAYS
from fork.integrators.wisdom_holman import SCHEME_DEFAULT
from fork.integrators.wisdom_holman import integrate as wh_integrate
from fork.model.bodies import (
    JUPITER,
    M_CENTRAL_MSUN,
    M_EARTH_MSUN,
    SATURN,
    URANUS,
    G,
)

TWO_PI = 2.0 * math.pi
# plain floats: numba freezes module-level floats as compile-time constants
_M_JUPITER = JUPITER.mass_msun
_M_SATURN = SATURN.mass_msun
_M_URANUS = URANUS.mass_msun
_MU_URANUS = G * (M_CENTRAL_MSUN + URANUS.mass_msun)

URANUS_PARAM_NAMES: tuple[str, ...] = (
    "a_au",
    "e",
    "inc_ecl_j2000_rad",
    "raan_ecl_j2000_rad",
    "varpi_ecl_j2000_rad",
    "lam0_ecl_j2000_rad",
)
PERTURBER_PARAM_NAMES: tuple[str, ...] = (
    "p_log10_m_earth",
    "p_a_au",
    "p_e",
    "p_varpi_ecl_j2000_rad",
    "p_lam0_ecl_j2000_rad",
)


@njit(cache=True)
def _uranus_lon_kernel(
    uranus_el: np.ndarray,
    perturber_el: np.ndarray,
    x_js: np.ndarray,
    v_js: np.ndarray,
    t_rel_days: np.ndarray,
    dt: float,
    scheme: int,
) -> np.ndarray:
    """Heliocentric ecl_j2000 longitude of Uranus [rad, in [0, 2π)] at epoch + t_rel_days."""
    n_body = 3 if perturber_el.size == 0 else 4
    x0 = np.empty((n_body, 3))
    v0 = np.empty((n_body, 3))
    masses = np.empty(n_body + 1)
    masses[0] = M_CENTRAL_MSUN
    masses[1] = _M_JUPITER
    masses[2] = _M_SATURN
    masses[3] = _M_URANUS
    x0[0] = x_js[0]
    v0[0] = v_js[0]
    x0[1] = x_js[1]
    v0[1] = v_js[1]

    a, e, inc, raan, varpi, lam0 = (
        uranus_el[0], uranus_el[1], uranus_el[2], uranus_el[3], uranus_el[4], uranus_el[5]
    )  # fmt: skip
    r, v = elements_to_state(a, e, inc, raan, varpi - raan, lam0 - varpi, _MU_URANUS)
    x0[2] = r
    v0[2] = v

    if n_body == 4:
        m_p = 10.0 ** perturber_el[0] * M_EARTH_MSUN
        masses[4] = m_p
        a_p, e_p, varpi_p, lam0_p = (
            perturber_el[1],
            perturber_el[2],
            perturber_el[3],
            perturber_el[4],
        )
        rp, vp = elements_to_state(
            a_p, e_p, 0.0, 0.0, varpi_p, lam0_p - varpi_p, G * (M_CENTRAL_MSUN + m_p)
        )
        x0[3] = rp
        v0[3] = vp

    pos = wh_integrate(x0, v0, masses, t_rel_days, dt, scheme)
    out = np.empty(t_rel_days.size)
    for k in range(t_rel_days.size):
        lon = math.atan2(pos[k, 2, 1], pos[k, 2, 0])
        out[k] = lon + TWO_PI if lon < 0.0 else lon
    return out


@dataclass(frozen=True)
class ForwardContext:
    """Everything a prediction needs besides θ: the epoch and the (unfitted) Jupiter/Saturn ICs."""

    epoch_jd_tdb: float
    x_js_ecl_j2000: np.ndarray  # (2, 3) heliocentric positions of Jupiter, Saturn [AU]
    v_js_ecl_j2000: np.ndarray  # (2, 3) heliocentric velocities [AU/day]
    dt_days: float = DEFAULT_DT_DAYS
    scheme: int = SCHEME_DEFAULT

    @classmethod
    def from_states(cls, epoch_jd_tdb: float, jupiter_state: ArrayLike, saturn_state: ArrayLike,
                    **kw) -> ForwardContext:  # fmt: skip
        st = np.array([jupiter_state, saturn_state], dtype=np.float64)
        return cls(float(epoch_jd_tdb), st[:, :3].copy(), st[:, 3:].copy(), **kw)

    @classmethod
    def from_horizons(cls, epoch_jd_tdb: float, **kw) -> ForwardContext:
        """Jupiter/Saturn ICs relative to the Sun + inner-planet barycenter (see bodies.py)."""
        from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about_bary

        jup = about_bary("jupiter", [epoch_jd_tdb])[0]
        sat = about_bary("saturn", [epoch_jd_tdb])[0]
        return cls.from_states(epoch_jd_tdb, jup, sat, **kw)


@dataclass(frozen=True)
class UranusModel:
    """Uranus + Jupiter + Saturn, no perturber. θ = the 6 Uranus elements."""

    context: ForwardContext
    param_names: tuple[str, ...] = field(default=URANUS_PARAM_NAMES, init=False)

    def predict(self, theta: ArrayLike, t_jd: ArrayLike) -> np.ndarray:
        theta = np.asarray(theta, dtype=np.float64)
        if theta.shape != (6,):
            raise ValueError(f"UranusModel expects 6 parameters, got shape {theta.shape}")
        return predict_uranus_lon(theta, t_jd, self.context)


def predict_uranus_lon(theta: ArrayLike, t_jd: ArrayLike, context: ForwardContext) -> np.ndarray:
    """Uranus heliocentric longitude ``lon_ecl_j2000_rad`` at TDB JDs ``t_jd``.

    ``theta`` has 6 entries (Uranus only) or 11 (Uranus + planar perturber); see module docstring.
    """
    theta = np.asarray(theta, dtype=np.float64)
    if theta.shape == (6,):
        uranus, pert = theta, np.empty(0)
    elif theta.shape == (11,):
        uranus, pert = theta[:6], theta[6:]
    else:
        raise ValueError(f"theta must have 6 or 11 entries, got shape {theta.shape}")
    t = np.asarray(t_jd, dtype=np.float64)
    t_rel = np.ascontiguousarray(np.atleast_1d(t).ravel() - context.epoch_jd_tdb)
    lon = _uranus_lon_kernel(
        np.ascontiguousarray(uranus),
        np.ascontiguousarray(pert),
        context.x_js_ecl_j2000,
        context.v_js_ecl_j2000,
        t_rel,
        context.dt_days,
        context.scheme,
    )
    return lon.reshape(t.shape)


def uranus_theta_from_state(state_ecl_j2000: ArrayLike) -> np.ndarray:
    """6-element Uranus θ from a heliocentric state. For tests/diagnostics only — a fit must never
    be initialised or held at modern Uranus (project rule 2)."""
    s = np.asarray(state_ecl_j2000, dtype=np.float64)
    a, e, inc, raan, argp, mean = state_to_elements(s[:3], s[3:], _MU_URANUS)
    return np.array([a, e, inc, raan, (raan + argp) % TWO_PI, (raan + argp + mean) % TWO_PI])


def perturber_theta_from_state(state_ecl_j2000: ArrayLike, mass_msun: float) -> np.ndarray:
    """Planar 5-element perturber θ from a 3-D state: the osculating orbit is rotated into the
    ecliptic (a, e, ϖ = Ω + ω, λ = ϖ + M kept; i dropped). Tests/diagnostics only."""
    s = np.asarray(state_ecl_j2000, dtype=np.float64)
    a, e, _inc, raan, argp, mean = state_to_elements(s[:3], s[3:], G * (M_CENTRAL_MSUN + mass_msun))
    varpi = (raan + argp) % TWO_PI
    return np.array([math.log10(mass_msun / M_EARTH_MSUN), a, e, varpi, (varpi + mean) % TWO_PI])
