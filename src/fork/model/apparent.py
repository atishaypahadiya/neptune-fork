"""Geocentric apparent RA/Dec of date: the quantity Le Verrier's Table n° 77 prints.

The model frame is ``ecl_j2000`` (Horizons' ecliptic and mean equinox of J2000, obliquity
84381.448″), with positions about the Sun + inner-planet barycenter. For an observation at TT
(≈ TDB, < 2 ms) time t:

1. **Light time.** ρ = r_U(t − τ) − r_E(t) with r_U(t − τ) ≈ r_U(t) − τ·v_U(t) (the neglected
   ½·a·τ² is ~1e-8 AU for Uranus), τ = |ρ|/c, iterated.
2. **Annual aberration**, from the Earth's velocity: u' ∝ u + β − (u·β)·u, β = v_E/c. Neglected
   terms are O(β²) ≈ 1e-8 rad (0.002″). v_E is taken about the inner barycenter, which moves
   ≲ 15 m/s about the solar-system barycenter: ≲ 0.01″.
3. **Rotation to the true equator and equinox of date:** ecl_j2000 → ICRS (−ε₀ about x), then
   IAU 2006/2000A bias-precession-nutation (``erfa.pnm06a``).

Not modelled: gravitational light deflection by the Sun (milliarcseconds at Uranus's
elongations in these data), diurnal parallax (≤ 0.5″; the printed places are already
geocentric), and frame-tie differences between Horizons' J2000 and ICRS (~0.02″).
``tests/test_apparent.py`` checks this against astropy's ICRS → TETE (which includes the
deflection) to < 0.05″.

Geometry that does not depend on θ (Earth's state and the rotation matrices at each epoch)
is computed once in ``ApparentGeometry``; a likelihood evaluation only adds the vector algebra.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import erfa
import numpy as np
from numpy.typing import ArrayLike

from fork.data.frames import C_AU_PER_DAY

OBLIQUITY_J2000_ARCSEC = 84381.448  # Horizons' ecliptic J2000 frame (IAU 1976 value)
TWO_PI = 2.0 * math.pi


def ecl_j2000_to_icrs_matrix() -> np.ndarray:
    """R with x_icrs = R @ x_ecl_j2000 (rotation by −ε₀ about the x axis)."""
    eps = math.radians(OBLIQUITY_J2000_ARCSEC / 3600.0)
    c, s = math.cos(eps), math.sin(eps)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


@dataclass(frozen=True)
class ApparentGeometry:
    """Per-epoch quantities independent of θ."""

    jd_tt: np.ndarray  # (n,)
    earth_pos_ecl_j2000_au: np.ndarray  # (n, 3), about the inner barycenter
    earth_vel_ecl_j2000_au_d: np.ndarray  # (n, 3)
    ecl_j2000_to_true_of_date: np.ndarray  # (n, 3, 3)

    @classmethod
    def from_earth_states(cls, jd_tt: ArrayLike, earth_state: ArrayLike) -> ApparentGeometry:
        jd = np.ascontiguousarray(np.asarray(jd_tt, dtype=np.float64))
        st = np.asarray(earth_state, dtype=np.float64)
        npb = erfa.pnm06a(jd, np.zeros_like(jd))  # GCRS → true equator and equinox of date
        m = np.ascontiguousarray(npb @ ecl_j2000_to_icrs_matrix())
        return cls(jd, np.ascontiguousarray(st[:, :3]), np.ascontiguousarray(st[:, 3:]), m)

    @classmethod
    def from_horizons(cls, jd_tt: ArrayLike) -> ApparentGeometry:
        """Earth from (cached) Horizons, about the inner barycenter like the model's bodies.
        TT is used as TDB (they differ by < 2 ms)."""
        from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about

        jd = np.asarray(jd_tt, dtype=np.float64)
        return cls.from_earth_states(jd, about("earth", jd))


def apparent_ra_dec(
    pos_ecl_j2000_au: np.ndarray, vel_ecl_j2000_au_d: np.ndarray, geom: ApparentGeometry
) -> tuple[np.ndarray, np.ndarray]:
    """(ra_true_of_date_rad in [0, 2π), dec_true_of_date_rad) of a body whose state at each
    epoch of ``geom`` is (pos, vel), each (n, 3), same origin as the Earth's."""
    r = np.asarray(pos_ecl_j2000_au, dtype=np.float64)
    v = np.asarray(vel_ecl_j2000_au_d, dtype=np.float64)
    rho = r - geom.earth_pos_ecl_j2000_au
    for _ in range(2):
        tau = np.linalg.norm(rho, axis=1) / C_AU_PER_DAY
        rho = r - v * tau[:, None] - geom.earth_pos_ecl_j2000_au
    u = rho / np.linalg.norm(rho, axis=1)[:, None]
    beta = geom.earth_vel_ecl_j2000_au_d / C_AU_PER_DAY
    u = u + beta - np.sum(u * beta, axis=1)[:, None] * u
    w = np.einsum("nij,nj->ni", geom.ecl_j2000_to_true_of_date, u)
    w /= np.linalg.norm(w, axis=1)[:, None]
    ra = np.mod(np.arctan2(w[:, 1], w[:, 0]), TWO_PI)
    dec = np.arcsin(np.clip(w[:, 2], -1.0, 1.0))
    return ra, dec
