"""M1: precession and state ↔ elements conversion."""

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from fork.data.frames import (
    ecl_j2000_to_ecl_of_date,
    ecl_of_date_to_ecl_j2000,
    elements_to_state,
    jd_tdb_from_iso,
    precession_matrix_ecl_j2000_to_ecl_of_date,
    solve_kepler_elliptic,
    state_to_elements,
)
from fork.model.bodies import URANUS, mu_heliocentric

RAD_TO_ARCSEC = math.degrees(1.0) * 3600.0
TWO_PI = 2.0 * math.pi
MU = mu_heliocentric(URANUS.mass_msun)


def wrap(d):
    """Wrap-safe angle difference (project rule 6)."""
    return (np.asarray(d) + math.pi) % TWO_PI - math.pi


# ------------------------------------------------------------------ precession


def test_precession_fixed_direction_100yr_is_5029_arcsec():
    """A fixed direction's ecliptic-of-date longitude advances ≈ 5029″ per century near J2000."""
    t0 = jd_tdb_from_iso("1950-01-01T12:00:00")
    t1 = jd_tdb_from_iso("2050-01-01T12:00:00")
    lon_ecl_j2000_rad = np.linspace(0.0, TWO_PI, 12, endpoint=False)
    lat_ecl_j2000_rad = np.zeros_like(lon_ecl_j2000_rad)
    lon0_ecl_of_date_rad, _ = ecl_j2000_to_ecl_of_date(lon_ecl_j2000_rad, lat_ecl_j2000_rad, t0)
    lon1_ecl_of_date_rad, _ = ecl_j2000_to_ecl_of_date(lon_ecl_j2000_rad, lat_ecl_j2000_rad, t1)
    dlon_arcsec = wrap(lon1_ecl_of_date_rad - lon0_ecl_of_date_rad) * RAD_TO_ARCSEC
    np.testing.assert_allclose(dlon_arcsec, 5029.0, atol=1.0)


@pytest.mark.parametrize("year", [1690, 1781, 1822, 1846])
def test_precession_matches_iau2006_pA_in_historical_era(year):
    """J2000 equinox direction expressed in ecliptic of date ≈ IAU 2006 general precession p_A.

    Tighter than the spec's ±1″ and evaluated where our data live (1690–1846), where the rate is
    not 5029″/cy: the quadratic term matters at the arcsecond level.
    """
    t = jd_tdb_from_iso(f"{year}-01-01T12:00:00")
    tc = (t - 2451545.0) / 36525.0
    p_a_arcsec = 5028.796195 * tc + 1.1054348 * tc**2 + 0.00007964 * tc**3
    lon_ecl_of_date_rad, _ = ecl_j2000_to_ecl_of_date(0.0, 0.0, t)
    assert wrap(lon_ecl_of_date_rad) * RAD_TO_ARCSEC == pytest.approx(p_a_arcsec, abs=0.05)


def test_precession_round_trip_and_orthonormal():
    rng = np.random.default_rng(1)
    n = 200
    t = rng.uniform(jd_tdb_from_iso("1690-01-01"), jd_tdb_from_iso("1850-01-01"), n)
    lon_ecl_j2000_rad = rng.uniform(0, TWO_PI, n)
    lat_ecl_j2000_rad = rng.uniform(-0.3, 0.3, n)
    lon_d, lat_d = ecl_j2000_to_ecl_of_date(lon_ecl_j2000_rad, lat_ecl_j2000_rad, t)
    lon_back, lat_back = ecl_of_date_to_ecl_j2000(lon_d, lat_d, t)
    assert np.max(np.abs(wrap(lon_back - lon_ecl_j2000_rad))) < 1e-12
    assert np.max(np.abs(lat_back - lat_ecl_j2000_rad)) < 1e-12

    rot = precession_matrix_ecl_j2000_to_ecl_of_date(t)
    assert rot.shape == (n, 3, 3)
    err = np.einsum("nij,nkj->nik", rot, rot) - np.eye(3)
    assert np.max(np.abs(err)) < 1e-14
    np.testing.assert_allclose(np.linalg.det(rot), 1.0, atol=1e-14)


def test_precession_scalar_shape():
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(2451545.0)
    assert rot.shape == (3, 3)


# ------------------------------------------------------------------ Kepler's equation


@settings(max_examples=500, deadline=None)
@given(m=st.floats(-20.0, 20.0), e=st.floats(0.0, 0.99))
def test_solve_kepler_elliptic_residual(m, e):
    ecc = solve_kepler_elliptic(m, e)
    assert abs(wrap(ecc - e * math.sin(ecc) - m)) < 1e-14


# ------------------------------------------------------------------ elements ↔ state


def _assert_elements_close(got, want, tol=1e-12):
    a, e, inc, raan, argp, mean = got
    a0, e0, inc0, raan0, argp0, mean0 = want
    assert abs(a - a0) / a0 < tol
    assert abs(e - e0) < tol
    assert abs(inc - inc0) < tol
    for x, x0 in ((raan, raan0), (argp, argp0), (mean, mean0)):
        assert abs(wrap(x - x0)) < tol


# Individual Ω, ω, M are only *defined* to ~ε/sin i and ~ε/e: merely rounding (r, v) to doubles
# perturbs them by that much. 1e-12 per element is therefore asserted where it is attainable
# (e, sin i ≥ 0.01, margin ~50×); near-circular / near-planar orbits are covered below by the
# well-conditioned combinations instead. At e = 1e-3 the observed ω error is ~1.0e-12.
@settings(max_examples=1000, deadline=None)
@given(
    a=st.floats(0.3, 100.0),
    e=st.floats(0.01, 0.95),
    inc=st.floats(0.01, math.pi - 0.01),
    raan=st.floats(0.0, TWO_PI, exclude_max=True),
    argp=st.floats(0.0, TWO_PI, exclude_max=True),
    mean=st.floats(0.0, TWO_PI, exclude_max=True),
)
def test_elements_state_elements_round_trip(a, e, inc, raan, argp, mean):
    r, v = elements_to_state(a, e, inc, raan, argp, mean, MU)
    _assert_elements_close(state_to_elements(r, v, MU), (a, e, inc, raan, argp, mean))


@settings(max_examples=1000, deadline=None)
@given(
    a=st.floats(0.3, 100.0),
    e=st.floats(0.0, 0.01),
    inc=st.floats(0.0, 0.01),
    raan=st.floats(0.0, TWO_PI, exclude_max=True),
    argp=st.floats(0.0, TWO_PI, exclude_max=True),
    mean=st.floats(0.0, TWO_PI, exclude_max=True),
)
def test_near_circular_near_planar_round_trip(a, e, inc, raan, argp, mean):
    """Small e, i: (a, e, i), λ, eccentricity vector and node vector all round-trip to 1e-12."""
    r, v = elements_to_state(a, e, inc, raan, argp, mean, MU)
    a1, e1, inc1, raan1, argp1, mean1 = state_to_elements(r, v, MU)
    assert abs(a1 - a) / a < 1e-12
    assert abs(e1 - e) < 1e-12
    assert abs(inc1 - inc) < 1e-12
    assert abs(wrap((raan1 + argp1 + mean1) - (raan + argp + mean))) < 1e-12  # λ
    varpi, varpi1 = raan + argp, raan1 + argp1
    assert abs(e1 * math.cos(varpi1) - e * math.cos(varpi)) < 1e-12
    assert abs(e1 * math.sin(varpi1) - e * math.sin(varpi)) < 1e-12
    assert abs(math.sin(inc1) * math.cos(raan1) - math.sin(inc) * math.cos(raan)) < 1e-12
    assert abs(math.sin(inc1) * math.sin(raan1) - math.sin(inc) * math.sin(raan)) < 1e-12


def test_uranus_like_round_trip():
    want = (19.2, 0.047, math.radians(0.77), math.radians(74.0), math.radians(96.0), 1.0)
    r, v = elements_to_state(*want, MU)
    assert 18.0 < np.linalg.norm(r) < 20.5
    _assert_elements_close(state_to_elements(r, v, MU), want)


@pytest.mark.parametrize("e", [0.0, 1e-9, 0.047])
@pytest.mark.parametrize("inc", [0.0, 1e-12, 0.0134])
def test_degenerate_orbits_preserve_varpi_and_lambda(e, inc):
    """Planar and/or circular orbits: Ω, ω, M individually undefined; ϖ and λ must survive."""
    raan, argp, mean = 1.3, 2.1, 0.4
    a = 19.2
    r, v = elements_to_state(a, e, inc, raan, argp, mean, MU)
    a1, e1, inc1, raan1, argp1, mean1 = state_to_elements(r, v, MU)
    assert abs(a1 - a) / a < 1e-12
    assert abs(e1 - e) < 1e-12
    assert abs(inc1 - inc) < 1e-12
    # λ = Ω + ω + M is always well defined
    assert abs(wrap((raan1 + argp1 + mean1) - (raan + argp + mean))) < 1e-12
    if e > 1e-6:  # ϖ = Ω + ω is defined whenever e is not ~0
        assert abs(wrap((raan1 + argp1) - (raan + argp))) < 1e-9
    r1, v1 = elements_to_state(a1, e1, inc1, raan1, argp1, mean1, MU)
    np.testing.assert_allclose(r1, r, rtol=0, atol=1e-12 * a)
    np.testing.assert_allclose(v1, v, rtol=0, atol=1e-12 * np.linalg.norm(v))


def test_state_elements_state_round_trip_random():
    rng = np.random.default_rng(7)
    for _ in range(500):
        a = rng.uniform(0.5, 60.0)
        el = (a, rng.uniform(0, 0.95), rng.uniform(0, math.pi), *rng.uniform(0, TWO_PI, 3))
        r, v = elements_to_state(*el, MU)
        r1, v1 = elements_to_state(*state_to_elements(r, v, MU), MU)
        np.testing.assert_allclose(r1, r, rtol=0, atol=1e-12 * a)
        np.testing.assert_allclose(v1, v, rtol=0, atol=1e-12 * np.linalg.norm(v))


def test_unbound_raises():
    r = np.array([1.0, 0.0, 0.0])
    v = np.array([0.0, 1.0, 0.0])  # far above escape speed for μ ≈ 3e-4
    with pytest.raises(ValueError):
        state_to_elements(r, v, MU)
