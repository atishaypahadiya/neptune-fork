"""Forward model for the inner solar system (Phase 2): θ → Mercury's heliocentric longitude.

Sun plus the eight planets, integrated with ``fork.integrators.wh_extra`` (the Phase 1
Wisdom–Holman map, fourth-order composition, with the optional forces). The seven planets
other than Mercury start from given states; Mercury's six elements are always in θ. Each
optional block of θ is present only when its flag is set, so one interface serves the
competing models:

    θ = (a, e, i, Ω, ω, M)                Mercury at the epoch, ecliptic J2000 [AU, rad]
        + (gr,)                           ``Flags.gr``: scale of the 1PN term (1 = GR)
        + (j2,)                           ``Flags.j2``: the Sun's J2, spin axis from the context
        + (ring_mass, ring_radius)        ``Flags.ring``: M_sun, AU; normal from the context
        + (mass, a, λ₀)                   ``Flags.planet``: a point planet on a circle in the
                                          plane with the context's normal, as an extra body
                                          (the Phase 1 unseen-planet treatment)

With every flag off this is the Newtonian N-body problem. Not used by any Phase 1 code.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from fork.data.frames import elements_to_state, state_to_elements
from fork.integrators import wh_extra as W
from fork.model.bodies import G, masses_msun

PLANETS = ("mercury", "venus", "earth_moon", "mars", "jupiter", "saturn", "uranus", "neptune")
DEFAULT_DT_DAYS = 1.0  # with order 4; docs/phase2/validation.md, V5
DEFAULT_ORDER = 4

OBLIQUITY_J2000_RAD = math.radians(84381.406 / 3600.0)  # IAU 2006


def solar_spin_axis_ecl_j2000() -> np.ndarray:
    """Unit vector of the Sun's north pole in ecliptic J2000: α₀ = 286.13°, δ₀ = 63.87°
    (IAU Working Group on Cartographic Coordinates, Archinal et al. 2018), rotated from the
    equator to the ecliptic of J2000."""
    ra, dec = math.radians(286.13), math.radians(63.87)
    x, y, z = math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)
    ce, se = math.cos(OBLIQUITY_J2000_RAD), math.sin(OBLIQUITY_J2000_RAD)
    return np.array([x, ce * y + se * z, -se * y + ce * z])


@dataclass(frozen=True)
class Flags:
    gr: bool = False
    j2: bool = False
    ring: bool = False
    planet: bool = False

    def n_theta(self) -> int:
        return 6 + self.gr + self.j2 + 2 * self.ring + 3 * self.planet


NO_FLAGS = Flags()


@dataclass(frozen=True)
class InnerSystemContext:
    """Epoch, the Sun-centred states (7, 6) of Venus … Neptune at it (ecliptic J2000), and
    the fixed directions."""

    epoch_jd_tdb: float
    other_states: np.ndarray
    spin_axis: np.ndarray = field(default_factory=solar_spin_axis_ecl_j2000)
    plane_normal: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.0, 1.0]))
    dt: float = DEFAULT_DT_DAYS
    order: int = DEFAULT_ORDER


def _plane_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = normal / np.linalg.norm(normal)
    ref = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = ref - n * (ref @ n)
    u /= np.linalg.norm(u)
    return u, np.cross(n, u)


def mercury_theta_from_state(state_ecl_j2000: ArrayLike) -> np.ndarray:
    """Mercury's six elements from a Sun-centred state (the inverse of what ``predict`` does
    with the first six entries of θ)."""
    s = np.asarray(state_ecl_j2000, dtype=np.float64)
    mu = G * (1.0 + masses_msun(("mercury",))[0])
    return np.array(state_to_elements(s[:3].copy(), s[3:].copy(), mu))


def states(
    theta: ArrayLike, t_jd: ArrayLike, context: InnerSystemContext, flags: Flags = NO_FLAGS
) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric (x, v) of Mercury at t_jd [JD TDB], each (T, 3), ecliptic J2000."""
    th = np.asarray(theta, dtype=np.float64)
    if th.size != flags.n_theta():
        raise ValueError(f"theta has {th.size} entries, the flags need {flags.n_theta()}")
    m = masses_msun(PLANETS)
    mu = G * (1.0 + m[0])
    r, v = elements_to_state(th[0], th[1], th[2], th[3], th[4], th[5], mu)
    x0 = np.vstack([r, context.other_states[:, :3]])
    v0 = np.vstack([v, context.other_states[:, 3:]])
    masses = np.concatenate([[1.0], m])
    k = 6
    kw: dict = {}
    if flags.gr:
        kw["gr"] = th[k]
        k += 1
    if flags.j2:
        kw["j2"] = th[k]
        kw["spin_axis"] = tuple(context.spin_axis)
        k += 1
    if flags.ring:
        kw["ring_mass"], kw["ring_radius_au"] = th[k], th[k + 1]
        kw["ring_normal"] = tuple(context.plane_normal)
        k += 2
    if flags.planet:
        mp, ap, lam = th[k], th[k + 1], th[k + 2]
        u, w = _plane_basis(context.plane_normal)
        speed = math.sqrt(G * (1.0 + mp) / ap)
        xp = ap * (math.cos(lam) * u + math.sin(lam) * w)
        vp = speed * (-math.sin(lam) * u + math.cos(lam) * w)
        x0, v0 = np.vstack([x0, xp]), np.vstack([v0, vp])
        masses = np.append(masses, mp)
    t = np.atleast_1d(np.asarray(t_jd, dtype=np.float64)) - context.epoch_jd_tdb
    x, vv = W.integrate_states(x0, v0, masses, t, context.dt, W.make_extra(**kw), context.order)
    return x[:, 0], vv[:, 0]


def predict(
    theta: ArrayLike, t_jd: ArrayLike, context: InnerSystemContext, flags: Flags = NO_FLAGS
) -> np.ndarray:
    """Mercury's heliocentric ecliptic longitude [rad, J2000] at t_jd: the project's generic
    ``predict(theta, t_jd) -> lon_rad`` with the flags and context bound by the caller."""
    x, _ = states(theta, t_jd, context, flags)
    return np.arctan2(x[:, 1], x[:, 0]) % (2.0 * math.pi)
