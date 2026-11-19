"""The contact time of a transit of Mercury with light time and aberration: the model a
real-data version of the Mercury study needs (docs/phase2/identifiability_full.md, § 8:
"light time and ΔT are the first required additions"). ΔT is in ``fork.data.delta_t``.

A contact is the instant the apparent separation of Mercury's centre from the Sun's, seen by
the observer, equals R☉ ∓ r☿ (interior, exterior), with the frozen plan's radii at 1 AU
(959.63″, 3.36″) scaled by the geometric distances.

- **Light time.** Mercury is seen where it was at t − τ, τ = |x☿(t − τ) − x_obs(t)| / c,
  solved by iteration. The Sun is the origin of the heliocentric frame, so its own light time
  moves nothing here.
- **Aberration.** Both directions are shifted by the observer's velocity, u' ∝ u + v/c (first
  order). The shift is about 20″ for each, but they nearly cancel: what is left in the
  separation is of order (v/c) × separation, under 0.1″.
- **Left out, with bounds.** The Sun's barycentric motion during the light time (about 13 m/s ×
  500 s, under 0.01″); the Sun's light deflection of Mercury, which lies in front of the Sun
  (under 10⁻⁴″); second-order aberration (under 0.01″ in the separation).
- **The observer** is any point whose heliocentric state is given: the geocentre, the
  Earth–Moon barycentre of the synthetic study, or (later) an observatory.

States are heliocentric, in AU and AU/day, at one reference epoch per contact (TDB); both
bodies are carried from it by their Keplerian acceleration about the Sun, as in the frozen
study. Nothing here reads or fits an observation.
"""

from __future__ import annotations

import math

import numpy as np

from fork.model.bodies import G
from fork.phase2 import contact_root as CROOT

C_AU_PER_DAY = 173.1446326846693  # IAU 2012: c × 86400 s / au
ARCSEC = math.pi / (180.0 * 3600.0)
R_SUN = 959.63 * ARCSEC  # at 1 AU, as the frozen plan
R_MERCURY = 3.36 * ARCSEC
MU_SUN = G  # Keplerian carry about the Sun (the bodies' own masses change it by < 2 × 10⁻⁷)
LIGHT_ITER = 4
ROOT_TOL = 1e-12  # rad
MAX_OFFSET_DAYS = 1.0


def carry(x0: np.ndarray, v0: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Position and velocity at offsets t [d] from the reference states, to second order."""
    a0 = -MU_SUN * x0 / np.linalg.norm(x0, axis=1)[:, None] ** 3
    tt = np.asarray(t, dtype=float)[:, None]
    return x0 + v0 * tt + 0.5 * a0 * tt**2, v0 + a0 * tt


def apparent(xm, vm, xo, vo, t, light_time: bool = True, aberration: bool = True):
    """Directions (unit vectors) and geometric distances of Mercury and of the Sun's centre
    from the observer at offsets t [d], and Mercury's light time τ [d]."""
    obs, v_obs = carry(xo, vo, t)
    tau = np.zeros(len(t))
    for _ in range(LIGHT_ITER if light_time else 1):
        merc, _ = carry(xm, vm, np.asarray(t) - tau)
        rho = merc - obs
        if light_time:
            tau = np.linalg.norm(rho, axis=1) / C_AU_PER_DAY
    sun = -obs
    d_m, d_s = np.linalg.norm(rho, axis=1), np.linalg.norm(sun, axis=1)
    u_m, u_s = rho / d_m[:, None], sun / d_s[:, None]
    if aberration:
        beta = v_obs / C_AU_PER_DAY
        u_m = (u_m + beta) / np.linalg.norm(u_m + beta, axis=1)[:, None]
        u_s = (u_s + beta) / np.linalg.norm(u_s + beta, axis=1)[:, None]
    return u_m, d_m, u_s, d_s, tau


def gap(xm, vm, xo, vo, exterior, t, light_time=True, aberration=True) -> np.ndarray:
    """Apparent separation of the centres minus the contact radius [rad] at offsets t [d]."""
    u_m, d_m, u_s, d_s, _ = apparent(xm, vm, xo, vo, t, light_time, aberration)
    sep = np.arctan2(np.linalg.norm(np.cross(u_m, u_s), axis=1), np.einsum("ij,ij->i", u_m, u_s))
    radius = R_SUN / d_s + np.where(exterior, 1.0, -1.0) * R_MERCURY / d_m
    return sep - radius


def contact_offsets(xm, vm, xo, vo, exterior, side, light_time=True,
                    aberration=True) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Offsets [s] of each contact (``side``: −1 ingress, contact I or II; +1 egress, III or IV)
    from its reference epoch and the root's residual [rad]; NaN where that contact does not
    exist within a day (``contact_root``)."""
    ext = np.asarray(exterior, dtype=bool)
    xm, vm, xo, vo = (np.atleast_2d(np.asarray(a, dtype=float)) for a in (xm, vm, xo, vo))

    def f(i, t):
        return gap(xm[i], vm[i], xo[i], vo[i], ext[i], t, light_time, aberration)

    t, res = CROOT.contact_roots(f, np.asarray(side), ROOT_TOL, MAX_OFFSET_DAYS)
    return t * 86400.0, res
