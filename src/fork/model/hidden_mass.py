"""Hidden-mass model: Uranus + one planar perturber (M5).

Forward-model θ (11) = 6 Uranus elements + perturber (log10 m [M⊕], a, e, ϖ, λ₀); see
``fork.model.forward``. Fits sample the perturber's eccentricity as

    (k, h) = (√e·cos ϖ, √e·sin ϖ)

(PHASE1_SPEC M5): a flat prior on the (k, h) disk is uniform in e and ϖ, and ϖ stays defined as
e → 0. Derived observable: the perturber's light-time-corrected geocentric ecliptic longitude on
Galle's discovery night, 1846-09-23.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from fork.data.frames import (
    geocentric_lon_lat_ecl_j2000_light_time,
    jd_tdb_from_iso,
    precession_matrix_ecl_j2000_to_ecl_of_date,
)
from fork.model.forward import (
    PERTURBER_PARAM_NAMES,
    URANUS_PARAM_NAMES,
    ForwardContext,
    predict_states_ecl_j2000,
    predict_uranus_lon,
)

TWO_PI = 2.0 * math.pi
I_PERTURBER = 3  # body index in predict_states_ecl_j2000 output: [J, S, U, P]
BODE_A_AU = 38.4  # Le Verrier's assumed distance (Bode's law)
PERTURBER_FIT_NAMES: tuple[str, ...] = (
    "p_log10_m_earth",
    "p_a_au",
    "p_sqrt_e_cos_varpi",
    "p_sqrt_e_sin_varpi",
    "p_lam0_ecl_j2000_rad",
)


def galle_jd_tdb() -> float:
    from fork.data.horizons import GALLE_EPOCH_ISO_TDB

    return float(jd_tdb_from_iso(GALLE_EPOCH_ISO_TDB))


def perturber_fit_to_theta(y: ArrayLike) -> np.ndarray:
    """(log10 m, a, √e cos ϖ, √e sin ϖ, λ₀) → forward-model (log10 m, a, e, ϖ, λ₀)."""
    log10_m, a, k, h, lam0 = np.asarray(y, dtype=np.float64)
    return np.array([log10_m, a, k * k + h * h, math.atan2(h, k) % TWO_PI, lam0 % TWO_PI])


def perturber_theta_to_fit(theta_p: ArrayLike) -> np.ndarray:
    log10_m, a, e, varpi, lam0 = np.asarray(theta_p, dtype=np.float64)
    s = math.sqrt(e)
    return np.array([log10_m, a, s * math.cos(varpi), s * math.sin(varpi), lam0])


@dataclass(frozen=True)
class HiddenMassModel:
    """Uranus + Jupiter + Saturn + one planar perturber. θ = 6 Uranus + 5 perturber elements."""

    context: ForwardContext
    param_names: tuple[str, ...] = field(
        default=URANUS_PARAM_NAMES + PERTURBER_PARAM_NAMES, init=False
    )

    def predict(self, theta: ArrayLike, t_jd: ArrayLike) -> np.ndarray:
        theta = np.asarray(theta, dtype=np.float64)
        if theta.shape != (11,):
            raise ValueError(f"HiddenMassModel expects 11 parameters, got shape {theta.shape}")
        return predict_uranus_lon(theta, t_jd, self.context)


@dataclass(frozen=True)
class GeocentricLongitude:
    jd_tdb: float
    lon_ecl_j2000_rad: float
    lon_ecl_of_date_rad: float
    light_time_days: float


def _geocentric(pos, vel, earth_pos, jd: float) -> GeocentricLongitude:
    lon, lat, tau = geocentric_lon_lat_ecl_j2000_light_time(pos, vel, earth_pos)
    u = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
    w = precession_matrix_ecl_j2000_to_ecl_of_date(jd) @ u
    return GeocentricLongitude(jd, lon, math.atan2(w[1], w[0]) % TWO_PI, tau)


def earth_pos_about_inner_barycenter(jd_tdb: float) -> np.ndarray:
    """Earth geocentre at ``jd_tdb`` from (cached) Horizons, in the model's origin."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000

    return states_about_inner_barycenter_ecl_j2000("earth", [jd_tdb])[0, :3]


def perturber_geocentric_longitude(
    theta: ArrayLike, context: ForwardContext, jd_tdb: float | None = None
) -> GeocentricLongitude:
    """Light-time-corrected geocentric ecliptic longitude of the model perturber (default: on
    1846-09-23, Galle's night). Earth from Horizons; both about the inner barycenter."""
    jd = galle_jd_tdb() if jd_tdb is None else float(jd_tdb)
    pos, vel = predict_states_ecl_j2000(theta, [jd], context)
    if pos.shape[1] <= I_PERTURBER:
        raise ValueError("theta has no perturber")
    return _geocentric(
        pos[0, I_PERTURBER], vel[0, I_PERTURBER], earth_pos_about_inner_barycenter(jd), jd
    )


def true_neptune_geocentric_longitude(jd_tdb: float | None = None) -> GeocentricLongitude:
    """The same observable for the real Neptune (Horizons), for comparison."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000

    jd = galle_jd_tdb() if jd_tdb is None else float(jd_tdb)
    nep = states_about_inner_barycenter_ecl_j2000("neptune", [jd])[0]
    return _geocentric(nep[:3], nep[3:], earth_pos_about_inner_barycenter(jd), jd)


# ------------------------------------------------------------------ Le Verrier's 1846 solution
JD_1847_JAN1_TDB = 2395662.5  # epoch of Le Verrier's elements (1847 Jan 1.0; TDB ≈ UT here)


@dataclass(frozen=True)
class PublishedSolution:
    """A historical perturber solution as published: a Keplerian orbit, angles in the ecliptic
    and mean equinox of ``epoch_jd`` (assumed; the sources do not state the equinox)."""

    label: str
    a_au: float
    e: float
    mass_msun: float
    varpi_ecl_of_date_deg: float
    mean_lon_ecl_of_date_deg: float  # at epoch_jd
    true_lon_ecl_of_date_deg: float  # heliocentric, at epoch_jd
    epoch_jd: float
    sources: tuple[str, ...]

    @property
    def log10_m_earth(self) -> float:
        from fork.model.bodies import M_EARTH_MSUN

        return math.log10(self.mass_msun / M_EARTH_MSUN)

    def varpi_ecl_j2000_rad(self) -> float:
        """ϖ referred to ecl_j2000: the direction (ϖ, β = 0) of date, rotated back."""
        w = math.radians(self.varpi_ecl_of_date_deg)
        r = precession_matrix_ecl_j2000_to_ecl_of_date(self.epoch_jd)
        v = r.T @ np.array([math.cos(w), math.sin(w), 0.0])
        return math.atan2(v[1], v[0]) % TWO_PI


# Le Verrier, Comptes rendus 23, 428 (31 Aug 1846), the solution Galle searched with.
# a, e, mass: Krajnović (2016), "The contrivance of Neptune", A&G 57, 5.28, arXiv:1610.06424,
# Table 2 (from Grosser 1962): a = 36.15 AU, e = 0.10761, mass 0.00011 M☉, P = 217.4 yr.
# ϖ, mean and true longitude at 1847 Jan 1: Encyclopædia Britannica (1911), "Neptune (planet)"
# (a = 36.154, e = 0.1076, ϖ = 284°45′, mean longitude 318°47′, true longitude 326°32′).
# Krajnović gives the mass to 2 figures only; the commonly quoted 1/9300 M☉ rounds to it and
# is used as a sensitivity case (LE_VERRIER_1846_MASS_ALT_MSUN).
LE_VERRIER_1846 = PublishedSolution(
    label="Le Verrier 1846 (31 Aug)",
    a_au=36.15,
    e=0.10761,
    mass_msun=0.00011,
    varpi_ecl_of_date_deg=284.0 + 45.0 / 60.0,
    mean_lon_ecl_of_date_deg=318.0 + 47.0 / 60.0,
    true_lon_ecl_of_date_deg=326.0 + 32.0 / 60.0,
    epoch_jd=JD_1847_JAN1_TDB,
    sources=(
        "Krajnović 2016, A&G 57, 5.28 (arXiv:1610.06424), Table 2, after Grosser 1962",
        "Encyclopædia Britannica 1911, 'Neptune (planet)', table of elements",
    ),
)
LE_VERRIER_1846_MASS_ALT_MSUN = 1.0 / 9300.0
