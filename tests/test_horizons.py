"""M1: Horizons access and cache. Offline tests use a fake fetcher; network tests hit JPL."""

import math

import numpy as np
import pytest

from fork.data import horizons
from fork.data.frames import jd_tdb_from_iso


class FakeFetcher:
    """Deterministic stand-in for Horizons that records every epoch requested."""

    def __init__(self):
        self.calls: list[np.ndarray] = []

    def __call__(self, horizons_id: str, jd_tdb: np.ndarray) -> np.ndarray:
        self.calls.append(jd_tdb.copy())
        base = float(horizons_id)
        z = np.zeros_like(jd_tdb)
        return np.column_stack([z + base, jd_tdb, -jd_tdb, z, z, z])


def test_cache_never_refetches(tmp_path):
    fake = FakeFetcher()
    jd = np.array([2400000.5, 2400010.5, 2400020.5])
    st1 = horizons.state_vectors_ecl_j2000("uranus", jd, directory=tmp_path, fetcher=fake)
    assert len(fake.calls) == 1
    assert (tmp_path / "horizons_uranus.npz").exists()

    st2 = horizons.state_vectors_ecl_j2000("uranus", jd[::-1], directory=tmp_path, fetcher=fake)
    assert len(fake.calls) == 1, "fully cached request must not call the fetcher"
    np.testing.assert_array_equal(st2, st1[::-1])


def test_partial_cache_fetches_only_missing(tmp_path):
    fake = FakeFetcher()
    horizons.state_vectors_ecl_j2000("saturn", [10.0, 20.0], directory=tmp_path, fetcher=fake)
    st = horizons.state_vectors_ecl_j2000(
        "saturn", [30.0, 10.0, 30.0, 20.0], directory=tmp_path, fetcher=fake
    )
    assert len(fake.calls) == 2
    np.testing.assert_array_equal(fake.calls[1], [30.0])
    np.testing.assert_array_equal(st[:, 1], [30.0, 10.0, 30.0, 20.0])
    assert np.all(st[:, 0] == 6.0)  # Saturn barycenter id


def test_bodies_cached_separately(tmp_path):
    fake = FakeFetcher()
    horizons.state_vectors_ecl_j2000("earth", [1.0], directory=tmp_path, fetcher=fake)
    horizons.state_vectors_ecl_j2000("jupiter", [1.0], directory=tmp_path, fetcher=fake)
    assert len(fake.calls) == 2
    st = horizons.state_vectors_ecl_j2000("earth", 1.0, directory=tmp_path, fetcher=fake)
    assert st.shape == (1, 6) and st[0, 0] == 399.0


def test_unknown_body(tmp_path):
    with pytest.raises(KeyError):
        horizons.state_vectors_ecl_j2000("pluto", [1.0], directory=tmp_path, fetcher=FakeFetcher())


def test_default_fetch_plan():
    plan = horizons.default_fetch_plan()
    # outer planets + Earth, plus the inner planets needed for the inner-barycenter correction
    assert set(plan) == {"jupiter", "saturn", "uranus", "neptune", "earth",
                         "mercury", "venus", "earth_moon", "mars"}  # fmt: skip
    a1 = jd_tdb_from_iso(horizons.A1_EPOCH_ISO_TDB)
    assert all(a1 in v for v in plan.values())
    assert plan["uranus"].size == 1846 - 1690 + 1 + 1


# ------------------------------------------------------------------ live JPL Horizons


@pytest.mark.network
def test_live_uranus_state_and_cache(tmp_path):
    jd = np.array([jd_tdb_from_iso("1781-03-13T00:00:00"), jd_tdb_from_iso("1846-09-23T00:00:00")])
    st = horizons.state_vectors_ecl_j2000("uranus", jd, directory=tmp_path)
    r = np.linalg.norm(st[:, :3], axis=1)
    v = np.linalg.norm(st[:, 3:], axis=1)
    assert np.all((18.2 < r) & (r < 20.2))
    assert np.all((3.2e-3 < v) & (v < 4.0e-3))  # ≈ 6.8 km/s in AU/day
    assert abs(math.degrees(math.atan2(st[0, 2], math.hypot(*st[0, :2])))) < 0.8  # near ecliptic

    # Cached now: a second call must not touch the network.
    def boom(*_):
        raise AssertionError("refetched a cached epoch")

    st2 = horizons.state_vectors_ecl_j2000("uranus", jd, directory=tmp_path, fetcher=boom)
    np.testing.assert_array_equal(st2, st)


def _sep_arcsec(a: np.ndarray, b: np.ndarray) -> float:
    """Angle between two vectors, accurate at tiny angles (acos is not: ~3 mas floor)."""
    return math.degrees(math.atan2(np.linalg.norm(np.cross(a, b)), a @ b)) * 3600.0


@pytest.mark.network
@pytest.mark.filterwarnings("ignore::erfa.ErfaWarning")  # pre-1900 Sun ephemeris; cancels
def test_live_horizons_ecliptic_is_astropy_iau76_frame():
    """Horizons 'ecliptic J2000' must equal astropy HeliocentricEclipticIAU76 (our ecl_j2000).

    Fetch the same Uranus vector in ICRF and in ecliptic, rotate the ICRF one into each candidate
    astropy ecliptic frame, and compare directions. IAU76 must match to < 1 mas, while
    HeliocentricMeanEcliptic(J2000) must *not* (frame bias + IAU2006 obliquity, ~tens of mas),
    which shows the test can tell the two apart.
    """
    import astropy.units as u
    from astropy.coordinates import (
        ICRS,
        CartesianRepresentation,
        HeliocentricEclipticIAU76,
        HeliocentricMeanEcliptic,
    )
    from astropy.time import Time
    from astroquery.jplhorizons import Horizons

    jd = jd_tdb_from_iso("1822-01-01T00:00:00")
    t = Time(jd, format="jd", scale="tdb")

    def vec(refplane):
        tab = Horizons(id="7", location="500@10", epochs=[jd]).vectors(
            refplane=refplane, aberrations="geometric", cache=False
        )
        return np.array([float(tab[c][0]) for c in ("x", "y", "z")])

    def rot_frame_to_icrs(frame):
        """Pure rotation part of frame → ICRS (subtract the image of the origin)."""
        cols = frame.realize_frame(CartesianRepresentation(np.eye(3) * u.AU)).transform_to(ICRS())
        origin = frame.realize_frame(CartesianRepresentation([0, 0, 0] * u.AU)).transform_to(ICRS())
        o = origin.cartesian.xyz.to_value(u.AU)
        return cols.cartesian.xyz.to_value(u.AU) - o[:, None]

    r_ecl_horizons = vec("ecliptic")
    r_icrf = vec("earth")  # Horizons 'earth' refplane = ICRF equatorial

    r_iau76 = rot_frame_to_icrs(HeliocentricEclipticIAU76(obstime=t)).T @ r_icrf
    r_hme = rot_frame_to_icrs(HeliocentricMeanEcliptic(equinox=Time("J2000"), obstime=t)).T @ r_icrf

    sep_iau76 = _sep_arcsec(r_ecl_horizons, r_iau76)
    sep_hme = _sep_arcsec(r_ecl_horizons, r_hme)
    print(f"Horizons ecl vs IAU76: {sep_iau76:.2e}″; vs HME(J2000): {sep_hme:.2e}″")
    assert sep_iau76 < 1e-3
    assert sep_hme > 5e-3


@pytest.mark.network
def test_live_a1_fixture_matches_horizons(tmp_path):
    """tests/data/horizons_a1_1781-03-13.json is exactly what Horizons serves."""
    import json
    from pathlib import Path

    data = json.loads((Path(__file__).parent / "data" / "horizons_a1_1781-03-13.json").read_text())
    for name, row in data["state_ecl_j2000"].items():
        live = horizons.state_vectors_ecl_j2000(name, [data["jd_tdb"]], directory=tmp_path)[0]
        np.testing.assert_allclose(live, row, rtol=1e-15, atol=0)


def test_inner_barycenter_is_mass_weighted_inner_planets(tmp_path):
    from fork.model.bodies import INNER_BODIES, M_CENTRAL_MSUN

    fake = FakeFetcher()  # state rows = (id, jd, -jd, 0, 0, 0)
    c = horizons.inner_barycenter_state_ecl_j2000([5.0], directory=tmp_path, fetcher=fake)
    want_x = sum(b.mass_msun * float(b.horizons_id) for b in INNER_BODIES) / M_CENTRAL_MSUN
    assert c.shape == (1, 6)
    assert c[0, 0] == pytest.approx(want_x, rel=1e-14)
    rel = horizons.states_about_inner_barycenter_ecl_j2000(
        "uranus", [5.0], directory=tmp_path, fetcher=fake
    )
    assert rel[0, 0] == pytest.approx(7.0 - want_x, rel=1e-14)


def test_match_rows_tolerates_echo_just_below_request():
    """Horizons can echo a non-round JD ~5e-10 d below the request; each requested epoch must
    still map to its own row, not the next one (seen with a linspace grid, 2026-11-08)."""
    req = np.array([2379044.364285714, 2379135.675, 2379226.9857142857])
    got = np.array([req[2], req[0] - 4.7e-10, req[1] + 4.7e-10])  # reordered, rounded
    idx = horizons._match_rows(req, got)
    np.testing.assert_array_equal(idx, [1, 2, 0])


def test_match_rows_rejects_missing_or_extra_epochs():
    req = np.array([2400000.5, 2400010.5])
    with pytest.raises(RuntimeError, match="do not match"):
        horizons._match_rows(req, np.array([2400000.5, 2400011.5]))
    with pytest.raises(RuntimeError, match="do not match"):
        horizons._match_rows(req, np.array([2400000.5]))
