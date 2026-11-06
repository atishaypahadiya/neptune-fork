"""Apparent geocentric RA/Dec of date (M6 observation model) against astropy (fast)."""

import math

import numpy as np
import pytest

from fork.data.horizons import state_vectors_ecl_j2000
from fork.model.apparent import ApparentGeometry, apparent_ra_dec, ecl_j2000_to_icrs_matrix

# Cached Horizons epochs (synthetic M4 schedule: Jan 1.0 TDB) spanning 1690–1845.
YEARS = (1690, 1715, 1750, 1769, 1781, 1800, 1820, 1830, 1845)


def _epochs():
    from fork.data.synthetic import jd_tdb_jan1

    return np.asarray(jd_tdb_jan1(YEARS))


def test_matches_astropy_icrs_to_tete():
    """Independent path: astropy places Uranus's barycentric ICRS position at the emission
    time (light time applied by us), then ICRS → TETE (true equator, true equinox of date)
    applies the Earth's position, aberration and solar light deflection from its own
    ephemeris and IAU 2006/2000A. Agreement to 0.05″ covers the deflection and frame-tie
    terms the model leaves out."""
    import astropy.units as u
    from astropy.coordinates import (
        ICRS,
        TETE,
        CartesianRepresentation,
        get_body_barycentric,
    )
    from astropy.time import Time

    from fork.data.frames import C_AU_PER_DAY

    jd = _epochs()
    uranus = state_vectors_ecl_j2000("uranus", jd)  # heliocentric
    earth = state_vectors_ecl_j2000("earth", jd)
    geom = ApparentGeometry.from_earth_states(jd, earth)
    ra, dec = apparent_ra_dec(uranus[:, :3], uranus[:, 3:], geom)

    t = Time(jd, format="jd", scale="tdb")
    rot = ecl_j2000_to_icrs_matrix()
    tau = np.linalg.norm(uranus[:, :3] - earth[:, :3], axis=1) / C_AU_PER_DAY
    for _ in range(3):
        rho = uranus[:, :3] - uranus[:, 3:] * tau[:, None] - earth[:, :3]
        tau = np.linalg.norm(rho, axis=1) / C_AU_PER_DAY
    t_emit = Time(jd - tau, format="jd", scale="tdb")
    sun_bary = get_body_barycentric("sun", t_emit).xyz.to_value(u.au).T
    r_emit_icrs = (rot @ (uranus[:, :3] - uranus[:, 3:] * tau[:, None]).T).T + sun_bary
    c = ICRS(CartesianRepresentation(r_emit_icrs.T * u.au))
    app = c.transform_to(TETE(obstime=t))
    d_ra = ((ra - app.ra.rad + math.pi) % (2 * math.pi) - math.pi) * np.cos(dec)
    d_dec = dec - app.dec.rad
    sep_arcsec = np.degrees(np.hypot(d_ra, d_dec)) * 3600
    assert sep_arcsec.max() < 0.05, sep_arcsec


def test_aberration_and_light_time_are_in_the_expected_range():
    """Sanity of the two corrections' sizes: aberration ≤ 20.5″, light time ~ 2.6–2.9 h."""
    jd = _epochs()
    uranus = state_vectors_ecl_j2000("uranus", jd)
    earth = state_vectors_ecl_j2000("earth", jd)
    geom = ApparentGeometry.from_earth_states(jd, earth)
    still = ApparentGeometry(geom.jd_tt, geom.earth_pos_ecl_j2000_au,
                             np.zeros_like(geom.earth_vel_ecl_j2000_au_d),
                             geom.ecl_j2000_to_true_of_date)  # fmt: skip
    ra, dec = apparent_ra_dec(uranus[:, :3], uranus[:, 3:], geom)
    ra0, dec0 = apparent_ra_dec(uranus[:, :3], uranus[:, 3:], still)
    ab = np.degrees(np.hypot((ra - ra0) * np.cos(dec), dec - dec0)) * 3600
    assert np.all(ab < 20.6) and ab.max() > 15.0
    ra1, dec1 = apparent_ra_dec(uranus[:, :3], np.zeros_like(uranus[:, 3:]), still)
    lt = np.degrees(np.hypot((ra0 - ra1) * np.cos(dec0), dec0 - dec1)) * 3600
    assert np.all(lt < 15.0) and lt.max() > 1.0  # Uranus's apparent motion over ~2.7 h


@pytest.mark.parametrize("year", [1690, 1845])
def test_true_of_date_differs_from_j2000_by_precession(year):
    """Rule 4 guard: an of-date RA differs from the J2000 RA by roughly 50″/yr × (2000 − year)
    · (cos ε + sin ε sin α tan δ); here only its order of magnitude, > 1.5°."""
    from fork.data.synthetic import jd_tdb_jan1

    jd = np.asarray(jd_tdb_jan1([year]))
    uranus = state_vectors_ecl_j2000("uranus", jd)
    earth = state_vectors_ecl_j2000("earth", jd)
    geom = ApparentGeometry.from_earth_states(jd, earth)
    ra, _ = apparent_ra_dec(uranus[:, :3], uranus[:, 3:], geom)
    j2000 = ApparentGeometry(jd, geom.earth_pos_ecl_j2000_au, geom.earth_vel_ecl_j2000_au_d,
                             ecl_j2000_to_icrs_matrix()[None])  # fmt: skip
    ra_j, _ = apparent_ra_dec(uranus[:, :3], uranus[:, 3:], j2000)
    assert abs((ra - ra_j + math.pi) % (2 * math.pi) - math.pi)[0] > math.radians(1.5)
