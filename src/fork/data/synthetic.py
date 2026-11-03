"""Synthetic Uranus longitude observations (M4; CLAUDE.md rule 3: synthetic before historical).

Truth = JPL Horizons heliocentric (Sun-centred) ecliptic-J2000 longitude of Uranus, which
includes the real Neptune. Gaussian noise is added: σ = 12″ for pre-discovery sightings
(1690–1771) and σ = 3″ from 1781 onward.

Epoch schedule (a stand-in until M6 brings the real dates):
- **Pre-discovery**: 22 sightings with year counts approximating the commonly cited record
  (Flamsteed 6, Bradley 3, Mayer 1, Le Monnier 12; 1690–1771). Each is placed at Jan 1.0 TDB of
  its year, so several share an epoch; their noise is independent.
- **Post-discovery**: one normal place per year at Jan 1.0 TDB, 1781–1845 (65 epochs).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

ARCSEC_TO_RAD = math.radians(1.0 / 3600.0)
SIGMA_PRE_DISCOVERY_ARCSEC = 12.0
SIGMA_POST_DISCOVERY_ARCSEC = 3.0
DISCOVERY_YEAR = 1781

# (year, number of sightings): approximate, see module docstring.
PRE_DISCOVERY_SIGHTINGS: tuple[tuple[int, int], ...] = (
    (1690, 1),
    (1712, 1),
    (1715, 4),
    (1748, 1),
    (1750, 3),
    (1753, 1),
    (1756, 1),
    (1764, 1),
    (1768, 6),
    (1769, 2),
    (1771, 1),
)
POST_DISCOVERY_YEARS: tuple[int, ...] = tuple(range(1781, 1846))


@dataclass(frozen=True)
class LongitudeObservations:
    """Heliocentric Uranus longitudes. Never carries the noise-free truth."""

    jd_tdb: np.ndarray
    lon_ecl_j2000_rad: np.ndarray
    sigma_rad: np.ndarray
    source: np.ndarray  # str per point

    def __post_init__(self) -> None:
        n = self.jd_tdb.size
        for name in ("lon_ecl_j2000_rad", "sigma_rad", "source"):
            if getattr(self, name).shape != (n,):
                raise ValueError(f"{name} must have shape ({n},)")


def jd_tdb_jan1(years: ArrayLike) -> np.ndarray:
    """JD (TDB) of Jan 1.0 of each (Gregorian) year."""
    from fork.data.frames import jd_tdb_from_iso

    return np.asarray(jd_tdb_from_iso([f"{int(y)}-01-01T00:00:00" for y in np.atleast_1d(years)]))


def synthetic_epochs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(jd_tdb, sigma_arcsec, source) of the synthetic schedule, time-ordered."""
    years: list[int] = []
    source: list[str] = []
    for year, count in PRE_DISCOVERY_SIGHTINGS:
        years += [year] * count
        source += ["synthetic_pre_discovery"] * count
    years += list(POST_DISCOVERY_YEARS)
    source += ["synthetic_post_discovery"] * len(POST_DISCOVERY_YEARS)
    year_arr = np.array(years)
    sigma = np.where(
        year_arr < DISCOVERY_YEAR, SIGMA_PRE_DISCOVERY_ARCSEC, SIGMA_POST_DISCOVERY_ARCSEC
    )
    return jd_tdb_jan1(year_arr), sigma.astype(np.float64), np.array(source)


StateFn = Callable[[str, np.ndarray], np.ndarray]


def true_uranus_lon_ecl_j2000_rad(jd_tdb: ArrayLike, state_fn: StateFn | None = None) -> np.ndarray:
    """Noise-free Horizons heliocentric longitude of Uranus in [0, 2π)."""
    if state_fn is None:
        from fork.data.horizons import state_vectors_ecl_j2000 as state_fn
    jd = np.asarray(jd_tdb, dtype=np.float64)
    uniq, inverse = np.unique(jd, return_inverse=True)
    st = state_fn("uranus", uniq)
    lon = np.mod(np.arctan2(st[:, 1], st[:, 0]), 2.0 * math.pi)
    return lon[inverse]


def make_synthetic_observations(
    rng: np.random.Generator, state_fn: StateFn | None = None
) -> tuple[LongitudeObservations, np.ndarray]:
    """(observations, truth_lon_ecl_j2000_rad). Keep the truth away from any fit."""
    jd, sigma_arcsec, source = synthetic_epochs()
    truth = true_uranus_lon_ecl_j2000_rad(jd, state_fn)
    sigma_rad = sigma_arcsec * ARCSEC_TO_RAD
    obs = np.mod(truth + rng.normal(0.0, 1.0, jd.size) * sigma_rad, 2.0 * math.pi)
    return LongitudeObservations(jd, obs, sigma_rad, source), truth
