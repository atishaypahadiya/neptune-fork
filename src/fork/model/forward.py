"""Forward model: θ → Uranus heliocentric ecliptic-J2000 longitudes.

Generic interface (Phase 2 adds other models beside this one)::

    model.param_names            # names of θ's entries, in order
    model.predict(theta, t_jd)   # -> lon_ecl_j2000_rad, shape like t_jd

Bodies integrated: Sun (with inner planets folded in, centred on their common barycenter),
Jupiter, Saturn, Uranus, and optionally one planar perturber. Jupiter and Saturn start from their
Horizons states at the model epoch and are **not** fitted. θ carries all six osculating
heliocentric Uranus elements at the epoch; fits free the in-plane four (a, e, ϖ, λ₀) and hold
i, Ω at their epoch values (CLAUDE.md rule 2, ``uranus_inc_raan_at_epoch``).

Uranus θ (6), all angles in ecl_j2000:
    a_au, e, inc_rad, raan_rad (Ω), varpi_rad (ϖ = Ω + ω), lam0_rad (mean longitude λ = ϖ + M)
Perturber θ (5, optional, planar: i = 0 in ecl_j2000):
    log10_m_earth, a_au, e, varpi_rad, lam0_rad
Inclined perturber θ (7, M6): the same five, then inc_rel_rad, raan_rel_rad, its inclination
and node relative to Uranus's orbital plane. All of (ϖ, λ₀, Ω_rel) are then measured in the
Uranus-plane frame F = Rz(Ω_U)·Rx(i_U)·Rz(−Ω_U) (Ω_U, i_U from θ's Uranus part). F maps
the ecliptic onto Uranus's plane about their line of nodes, so it tends to ecl_j2000 as
i_U → 0, and longitudes in F differ from ecliptic ones by O(i_U²) (i_U = 0.77°). Why
inclined: declinations see the perturber's out-of-plane pull, ~2″ by the 1840s for Neptune,
which heliocentric longitudes hardly do (0.017″). See ``tests/test_observation.py``.

Elements are converted with μ = G (M_central + m_body), the body's own heliocentric two-body μ.
The integrator is ``wisdom_holman`` with its project defaults (dt = 12 d, REBOUND composition).
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
from fork.integrators.wisdom_holman import integrate_states as wh_integrate_states
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
INCLINED_PERTURBER_PARAM_NAMES: tuple[str, ...] = (
    "p_log10_m_earth",
    "p_a_au",
    "p_e",
    "p_varpi_uranus_frame_rad",
    "p_lam0_uranus_frame_rad",
    "p_inc_rel_uranus_plane_rad",
    "p_raan_rel_uranus_plane_rad",
)


@njit(cache=True)
def uranus_frame_matrix(inc_u: float, raan_u: float) -> np.ndarray:
    """F = Rz(Ω)·Rx(i)·Rz(−Ω): x_ecl_j2000 = F @ x_uranus_frame."""
    ci, si = math.cos(inc_u), math.sin(inc_u)
    co, so = math.cos(raan_u), math.sin(raan_u)
    rz = np.array([[co, -so, 0.0], [so, co, 0.0], [0.0, 0.0, 1.0]])
    rx = np.array([[1.0, 0.0, 0.0], [0.0, ci, -si], [0.0, si, ci]])
    rzi = np.array([[co, so, 0.0], [-so, co, 0.0], [0.0, 0.0, 1.0]])
    out = np.zeros((3, 3))
    tmp = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                tmp[i, j] += rx[i, k] * rzi[k, j]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                out[i, j] += rz[i, k] * tmp[k, j]
    return out


@njit(cache=True)
def _initial_conditions(
    uranus_el: np.ndarray, perturber_el: np.ndarray, x_js: np.ndarray, v_js: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(x0, v0, masses) for [J, S, U, (P)] from θ; heliocentric about the inner barycenter."""
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
        mu_p = G * (M_CENTRAL_MSUN + m_p)
        if perturber_el.size == 7:
            inc_r, node_r = perturber_el[5], perturber_el[6]
            rf, vf = elements_to_state(a_p, e_p, inc_r, node_r, varpi_p - node_r,
                                       lam0_p - varpi_p, mu_p)  # fmt: skip
            f = uranus_frame_matrix(inc, raan)
            for i in range(3):
                x0[3, i] = f[i, 0] * rf[0] + f[i, 1] * rf[1] + f[i, 2] * rf[2]
                v0[3, i] = f[i, 0] * vf[0] + f[i, 1] * vf[1] + f[i, 2] * vf[2]
        else:
            rp, vp = elements_to_state(a_p, e_p, 0.0, 0.0, varpi_p, lam0_p - varpi_p, mu_p)
            x0[3] = rp
            v0[3] = vp
    return x0, v0, masses


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
    x0, v0, masses = _initial_conditions(uranus_el, perturber_el, x_js, v_js)
    pos = wh_integrate(x0, v0, masses, t_rel_days, dt, scheme)
    out = np.empty(t_rel_days.size)
    for k in range(t_rel_days.size):
        lon = math.atan2(pos[k, 2, 1], pos[k, 2, 0])
        out[k] = lon + TWO_PI if lon < 0.0 else lon
    return out


@njit(cache=True)
def _states_kernel(
    uranus_el: np.ndarray,
    perturber_el: np.ndarray,
    x_js: np.ndarray,
    v_js: np.ndarray,
    t_rel_days: np.ndarray,
    dt: float,
    scheme: int,
) -> tuple[np.ndarray, np.ndarray]:
    x0, v0, masses = _initial_conditions(uranus_el, perturber_el, x_js, v_js)
    return wh_integrate_states(x0, v0, masses, t_rel_days, dt, scheme)


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


def _split_theta(theta: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    theta = np.asarray(theta, dtype=np.float64)
    if theta.shape == (6,):
        return np.ascontiguousarray(theta), np.empty(0)
    if theta.shape in ((11,), (13,)):
        return np.ascontiguousarray(theta[:6]), np.ascontiguousarray(theta[6:])
    raise ValueError(f"theta must have 6, 11 or 13 entries, got shape {theta.shape}")


def predict_uranus_lon(theta: ArrayLike, t_jd: ArrayLike, context: ForwardContext) -> np.ndarray:
    """Uranus heliocentric longitude ``lon_ecl_j2000_rad`` at TDB JDs ``t_jd``.

    ``theta`` has 6 entries (Uranus only) or 11 (Uranus + planar perturber); see module docstring.
    """
    uranus, pert = _split_theta(theta)
    t = np.asarray(t_jd, dtype=np.float64)
    t_rel = np.ascontiguousarray(np.atleast_1d(t).ravel() - context.epoch_jd_tdb)
    lon = _uranus_lon_kernel(uranus, pert, context.x_js_ecl_j2000, context.v_js_ecl_j2000, t_rel,
                             context.dt_days, context.scheme)  # fmt: skip
    return lon.reshape(t.shape)


def predict_states_ecl_j2000(
    theta: ArrayLike, t_jd: ArrayLike, context: ForwardContext
) -> tuple[np.ndarray, np.ndarray]:
    """(pos [AU], vel [AU/day]) of [Jupiter, Saturn, Uranus, (perturber)] at TDB JDs ``t_jd``,
    shape (n_t, n_body, 3), heliocentric about the inner barycenter, ecl_j2000."""
    uranus, pert = _split_theta(theta)
    t_rel = np.ascontiguousarray(np.atleast_1d(np.asarray(t_jd, dtype=np.float64)).ravel()
                                 - context.epoch_jd_tdb)  # fmt: skip
    return _states_kernel(uranus, pert, context.x_js_ecl_j2000, context.v_js_ecl_j2000, t_rel,
                          context.dt_days, context.scheme)  # fmt: skip


def uranus_theta_from_state(state_ecl_j2000: ArrayLike) -> np.ndarray:
    """6-element Uranus θ from a heliocentric state. For tests/diagnostics only — a fit must never
    be initialised or held at Horizons in-plane elements (CLAUDE.md rule 2)."""
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


def inclined_perturber_theta_from_state(
    state_ecl_j2000: ArrayLike, mass_msun: float, uranus_inc_raan: tuple[float, float]
) -> np.ndarray:
    """7-element inclined perturber θ from a 3-D state, in the Uranus-plane frame of
    ``uranus_inc_raan`` (see module docstring). Tests/diagnostics only."""
    s = np.asarray(state_ecl_j2000, dtype=np.float64)
    f = uranus_frame_matrix(float(uranus_inc_raan[0]), float(uranus_inc_raan[1]))
    r, v = f.T @ s[:3], f.T @ s[3:]
    a, e, inc, raan, argp, mean = state_to_elements(r, v, G * (M_CENTRAL_MSUN + mass_msun))
    varpi = (raan + argp) % TWO_PI
    return np.array([math.log10(mass_msun / M_EARTH_MSUN), a, e, varpi, (varpi + mean) % TWO_PI,
                     inc, raan % TWO_PI])  # fmt: skip


def uranus_inc_raan_at_epoch(epoch_jd_tdb: float) -> tuple[float, float]:
    """(inc_ecl_j2000_rad, raan_ecl_j2000_rad): Uranus's osculating i, Ω at the epoch, from
    Horizons (about the inner barycenter). The only Horizons Uranus information a fit may use
    (CLAUDE.md rule 2); the leak is bounded by ``tests/test_forward.py::test_inc_raan_leak``."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about_bary

    theta = uranus_theta_from_state(about_bary("uranus", [epoch_jd_tdb])[0])
    return float(theta[2]), float(theta[3])
