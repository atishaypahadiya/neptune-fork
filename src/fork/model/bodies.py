"""Physical constants and body definitions.

Units throughout: AU, days, solar masses. G = k² with k the Gaussian gravitational constant, so
G·M_sun = k² AU³/day² exactly.

Planet masses are *system* masses (planet + satellites), given as Sun/planet ratios (IAU 2009
system of astronomical constants). Those are the masses that belong with Horizons barycenter
states (Jupiter barycenter "5", etc.).

Inner planets are not integrated; their mass is folded into the central body
(`M_CENTRAL_MSUN`). Over 1690–1846 the inner planets act on Uranus essentially as extra solar
mass, so this is the standard approximation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Gaussian gravitational constant [AU^{3/2} day^{-1} M_sun^{-1/2}].
K_GAUSS: float = 0.01720209895
G: float = K_GAUSS**2  # [AU^3 day^-2 M_sun^-1]

# Sun / planet mass ratios (planet systems). Outer planets: required by PHASE1_SPEC M1.
MASS_RATIO_JUPITER: float = 1047.348644
MASS_RATIO_SATURN: float = 3497.9018
MASS_RATIO_URANUS: float = 22902.98
MASS_RATIO_NEPTUNE: float = 19412.26

# Inner planets (IAU 2009), only used to fold into the central mass.
MASS_RATIO_MERCURY: float = 6023600.0
MASS_RATIO_VENUS: float = 408523.719
MASS_RATIO_EARTH_MOON: float = 328900.56
MASS_RATIO_MARS: float = 3098703.59

M_SUN_MSUN: float = 1.0
M_INNER_MSUN: float = (
    1.0 / MASS_RATIO_MERCURY
    + 1.0 / MASS_RATIO_VENUS
    + 1.0 / MASS_RATIO_EARTH_MOON
    + 1.0 / MASS_RATIO_MARS
)
M_CENTRAL_MSUN: float = M_SUN_MSUN + M_INNER_MSUN

# 1 Earth mass (Earth alone, no Moon) in solar masses, for perturber priors quoted in M⊕.
# GM_sun / GM_earth = 332946.0487 (IAU 2009).
M_EARTH_MSUN: float = 1.0 / 332946.0487


@dataclass(frozen=True)
class Body:
    """A massive body: name, JPL Horizons target id, and mass in solar masses."""

    name: str
    horizons_id: str
    mass_msun: float


JUPITER = Body("jupiter", "5", 1.0 / MASS_RATIO_JUPITER)  # Jupiter system barycenter
SATURN = Body("saturn", "6", 1.0 / MASS_RATIO_SATURN)
URANUS = Body("uranus", "7", 1.0 / MASS_RATIO_URANUS)
NEPTUNE = Body("neptune", "8", 1.0 / MASS_RATIO_NEPTUNE)
# Earth geocenter (not the Earth–Moon barycenter): used only as an observing site, never as a
# perturber, so its mass is irrelevant to dynamics and is already folded into M_CENTRAL_MSUN.
EARTH = Body("earth", "399", 0.0)

BODIES: dict[str, Body] = {b.name: b for b in (JUPITER, SATURN, URANUS, NEPTUNE, EARTH)}


def masses_msun(names: list[str] | tuple[str, ...]) -> np.ndarray:
    """Masses [M_sun] of the named bodies, in order, as a float64 array."""
    return np.array([BODIES[n].mass_msun for n in names], dtype=np.float64)


def mu_heliocentric(mass_msun: float) -> float:
    """G·(M_central + m): the two-body μ for a planet's heliocentric Keplerian orbit."""
    return G * (M_CENTRAL_MSUN + mass_msun)
