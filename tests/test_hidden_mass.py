"""M5a: hidden-mass model, perturber parameterization, 1846-09-23 geocentric longitude."""

import math

import numpy as np
import pytest

from fork.data.frames import C_AU_PER_DAY, geocentric_lon_lat_ecl_j2000_light_time
from fork.data.horizons import GALLE_EPOCH_ISO_TDB, states_about_inner_barycenter_ecl_j2000
from fork.data.synthetic import jd_tdb_jan1
from fork.inference.fit import wrap_rad
from fork.model import bodies
from fork.model.forward import (
    ForwardContext,
    perturber_theta_from_state,
    predict_states_ecl_j2000,
    predict_uranus_lon,
    uranus_theta_from_state,
)
from fork.model.hidden_mass import (
    HiddenMassModel,
    galle_jd_tdb,
    perturber_fit_to_theta,
    perturber_geocentric_longitude,
    perturber_theta_to_fit,
    true_neptune_geocentric_longitude,
)

ARC = math.degrees(1.0) * 3600.0
EPOCH = 2371628.5
# JPL Horizons astrometric (light-time corrected) geocentric RA/Dec of Neptune (id 8), ICRF, at
# 1846-09-23 22:59:54 UT = 23:00:00 TDB (fetched 2026-11-03, extra_precision).
HORIZONS_NEPTUNE_GALLE_RA_DEG = 330.378621006
HORIZONS_NEPTUNE_GALLE_DEC_DEG = -12.672344563
OBLIQUITY_J2000_RAD = math.radians(84381.448 / 3600.0)  # Horizons' ecliptic

pytestmark = pytest.mark.usefixtures("offline_horizons")


def _equatorial_to_ecl_lon(ra_deg: float, dec_deg: float) -> float:
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    v = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    c, s = math.cos(OBLIQUITY_J2000_RAD), math.sin(OBLIQUITY_J2000_RAD)
    w = np.array([v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]])
    return math.atan2(w[1], w[0]) % (2 * math.pi)


@pytest.fixture(scope="module")
def ctx(offline_horizons):
    return ForwardContext.from_horizons(EPOCH)


@pytest.fixture(scope="module")
def theta_planted(offline_horizons):
    u = uranus_theta_from_state(states_about_inner_barycenter_ecl_j2000("uranus", [EPOCH])[0])
    n = states_about_inner_barycenter_ecl_j2000("neptune", [EPOCH])[0]
    return np.concatenate([u, perturber_theta_from_state(n, bodies.NEPTUNE.mass_msun)])


def test_sqrt_e_parameterization_round_trip():
    for th in ([1.23, 30.2, 0.0039, 1.27, 3.27], [0.5, 38.4, 0.3, 5.9, 0.1]):
        back = perturber_fit_to_theta(perturber_theta_to_fit(th))
        np.testing.assert_allclose(back, th, rtol=0, atol=1e-13)
    y = perturber_theta_to_fit([1.0, 38.4, 0.25, math.radians(30.0), 0.0])
    assert math.hypot(y[2], y[3]) == pytest.approx(0.5)  # √e
    circ = perturber_fit_to_theta([1.0, 38.4, 0.0, 0.0, 7.0])
    assert circ[2] == 0.0 and 0.0 <= circ[4] < 2 * math.pi


def test_hidden_mass_model_interface(ctx, theta_planted):
    model = HiddenMassModel(ctx)
    assert len(model.param_names) == 11 and model.param_names[6] == "p_log10_m_earth"
    jd = jd_tdb_jan1([1781, 1822, 1845])
    np.testing.assert_array_equal(model.predict(theta_planted, jd),
                                  predict_uranus_lon(theta_planted, jd, ctx))  # fmt: skip
    with pytest.raises(ValueError):
        model.predict(theta_planted[:6], jd)


def test_state_prediction_matches_longitude_prediction(ctx, theta_planted):
    jd = jd_tdb_jan1(np.arange(1690, 1847, 13))
    for theta in (theta_planted[:6], theta_planted):
        pos, vel = predict_states_ecl_j2000(theta, jd, ctx)
        assert pos.shape == vel.shape == (jd.size, 3 if theta.size == 6 else 4, 3)
        lon = np.arctan2(pos[:, 2, 1], pos[:, 2, 0])
        assert np.max(np.abs(wrap_rad(lon - predict_uranus_lon(theta, jd, ctx)))) == 0.0


def test_light_time_geometry():
    """Static target: geometric direction; τ = distance / c. Moving target: shifted by v·τ."""
    r_t, r_e = np.array([30.0, 4.0, 0.3]), np.array([0.9, -0.4, 0.0])
    lon, lat, tau = geocentric_lon_lat_ecl_j2000_light_time(r_t, np.zeros(3), r_e)
    assert lat == pytest.approx(math.atan2(0.3, math.hypot(29.1, 4.4)), abs=1e-15)
    rho = r_t - r_e
    assert tau == pytest.approx(np.linalg.norm(rho) / C_AU_PER_DAY, rel=1e-15)
    assert lon == pytest.approx(math.atan2(rho[1], rho[0]) % (2 * math.pi), abs=1e-15)
    v = np.array([0.0, 0.003, 0.0])
    lon2, _, tau2 = geocentric_lon_lat_ecl_j2000_light_time(r_t, v, r_e)
    rho2 = r_t - v * tau2 - r_e
    assert lon2 == pytest.approx(math.atan2(rho2[1], rho2[0]), abs=1e-15) and lon2 < lon
    assert C_AU_PER_DAY == pytest.approx(173.14463267424034, rel=1e-15)


def test_true_neptune_galle_longitude_matches_horizons_astrometric():
    """Our light-time-corrected geocentric longitude of Neptune on Galle's night agrees with
    Horizons' own astrometric position (measured 0.001″; Horizons also applies light deflection,
    ~mas here)."""
    g = true_neptune_geocentric_longitude()
    assert g.jd_tdb == galle_jd_tdb() == pytest.approx(2395563.4583333335)
    want = _equatorial_to_ecl_lon(HORIZONS_NEPTUNE_GALLE_RA_DEG, HORIZONS_NEPTUNE_GALLE_DEC_DEG)
    assert abs(wrap_rad(g.lon_ecl_j2000_rad - want)) * ARC < 0.05
    assert g.light_time_days == pytest.approx(0.1685, abs=1e-3)
    # of date (mean equinox of 1846.7): 2.14° of precession less than J2000
    assert math.degrees(g.lon_ecl_of_date_rad) == pytest.approx(325.873, abs=1e-3)


def test_planted_neptune_galle_longitude_matches_truth(ctx, theta_planted):
    """The model's planar Neptune, integrated 1781 → 1846-09-23, lands on the real one.

    Measured 31″ (0.009°): flattening Neptune's 1.77° inclination into the ecliptic moves its
    ecliptic longitude by up to tan²(i/2) ≈ 49″. Negligible against the 1–2° A3 targets."""
    g = perturber_geocentric_longitude(theta_planted, ctx)
    nep = true_neptune_geocentric_longitude()
    assert abs(wrap_rad(g.lon_ecl_of_date_rad - nep.lon_ecl_of_date_rad)) * ARC < 60.0
    with pytest.raises(ValueError):
        perturber_geocentric_longitude(theta_planted[:6], ctx)


@pytest.mark.network
def test_live_horizons_astrometric_neptune_galle():
    from astroquery.jplhorizons import Horizons

    from fork.data.frames import jd_tdb_from_iso

    jd = float(jd_tdb_from_iso(GALLE_EPOCH_ISO_TDB))
    eph = Horizons(id="8", location="500", epochs=[jd - 6.0 / 86400.0]).ephemerides(
        quantities="1", extra_precision=True
    )
    assert float(eph["RA"][0]) == pytest.approx(HORIZONS_NEPTUNE_GALLE_RA_DEG, abs=1e-7)
    assert float(eph["DEC"][0]) == pytest.approx(HORIZONS_NEPTUNE_GALLE_DEC_DEG, abs=1e-7)
