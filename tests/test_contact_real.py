"""ΔT and the light-time contact model against published values (no observation of the
historical transits is read or fitted).

Tolerances were set from the error budgets in the docstrings before the comparisons were run:
- light time against Horizons' astrometric vectors: the Sun's barycentric motion during the
  light time, under 0.01″, so 0.02″;
- contacts against Espenak's published geocentric times: the times are printed to 1 s, the ΔT
  used for them is not printed (a predicted value; IERS gives the one used here), and his solar
  and planetary radii are taken to be the plan's, so 2 s.
"""

import math
from pathlib import Path

import numpy as np
import pytest

from fork.data import delta_t as DT
from fork.phase2 import contact_real as CR

FIXTURE = Path(__file__).parent / "fixtures" / "phase2" / "contact_real_horizons.npz"
TOL_LIGHT_ARCSEC = 0.02
TOL_CONTACT_S = 2.0


# ---------------------------------------------------------------------------- ΔT


@pytest.mark.parametrize(("year", "value", "sigma"), DT.PUBLISHED)
def test_delta_t_reproduces_the_published_table(year, value, sigma):
    """Each polynomial value lies within the published standard error of the tabulated value
    (rounded to whole seconds, so within σ + 0.5 s)."""
    assert abs(float(DT.delta_t(year)) - value) <= sigma + 0.5
    assert float(DT.delta_t_sigma(year)) == pytest.approx(sigma)


def test_delta_t_segments_join():
    for (start, _), _ in zip(DT.SEGMENTS[1:], range(3), strict=False):
        below = float(DT.delta_t(start - 1e-9))
        above = float(DT.delta_t(start))
        assert abs(below - above) < 0.5, start


def test_delta_t_over_the_transit_span_is_the_plans():
    """The frozen plan's § 8: about 9 to 18 s over 1677 to 1800."""
    y = np.linspace(1677.0, 1800.0, 124)
    d = DT.delta_t(y)
    assert d.min() == pytest.approx(8.8, abs=0.2) and d.max() == pytest.approx(17.7, abs=0.2)
    assert float(DT.delta_t(1677.0)) == pytest.approx(17.69, abs=0.01)
    assert float(DT.delta_t_sigma(1677.0)) == pytest.approx(20.0 - 0.77 * 15.0, abs=1e-9)


def test_delta_t_refuses_outside_its_span():
    for y in (1599.9, 1900.1, math.nan):
        with pytest.raises(ValueError):
            DT.delta_t(y)


@pytest.mark.filterwarnings("ignore::erfa.ErfaWarning")  # pre-1960 UTC: not used here
def test_tdb_minus_tt_against_astropy():
    from astropy.time import Time

    jd = np.linspace(2333000.5, 2415020.5, 7)  # 1675 to 1900
    t = Time(jd, format="jd", scale="tt")
    full = (t.tdb.jd1 - t.tt.jd1 + t.tdb.jd2 - t.tt.jd2) * 86400.0
    np.testing.assert_allclose(DT.tdb_minus_tt(jd), full, rtol=0, atol=5e-5)


def test_ut_to_tdb_and_back():
    jd_ut = np.array([2333500.25, 2374000.75, 2408000.5])
    jd_tdb = DT.ut_to_tdb(jd_ut)
    dt = (jd_tdb - jd_ut) * 86400.0
    np.testing.assert_allclose(dt, DT.delta_t(DT.decimal_year(jd_ut)), rtol=0, atol=0.002)
    np.testing.assert_allclose(DT.tdb_to_ut(jd_tdb), jd_ut, rtol=0, atol=1e-10)


# ---------------------------------------------------------------------------- light time


@pytest.fixture(scope="module")
def horizons():
    with np.load(FIXTURE) as z:
        return {k: z[k] for k in z.files}


def _angle(a, b):
    return np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.einsum("ij,ij->i", a, b))


def test_light_time_matches_horizons_astrometric(horizons):
    m, e = horizons["mercury"], horizons["earth"]
    t = np.zeros(len(m))
    u, _, _, _, tau = CR.apparent(m[:, :3], m[:, 3:], e[:, :3], e[:, 3:], t, aberration=False)
    ref = horizons["mercury_astrometric_geocentric"]
    off = _angle(u, ref) / CR.ARCSEC
    assert off.max() < TOL_LIGHT_ARCSEC
    np.testing.assert_allclose(tau * 86400.0, horizons["light_time_days"] * 86400.0, atol=1e-3)
    geo, _, _, _, _ = CR.apparent(m[:, :3], m[:, 3:], e[:, :3], e[:, 3:], t, light_time=False,
                                  aberration=False)  # fmt: skip
    assert (_angle(geo, ref) / CR.ARCSEC).min() > 5.0  # light time is not small


# ---------------------------------------------------------------------------- contacts


def _offsets(horizons, **kw):
    m, e = horizons["mercury"], horizons["earth"]
    k = np.array([int(s.split()[1]) for s in horizons["labels"]])
    off, res = CR.contact_offsets(m[:, :3], m[:, 3:], e[:, :3], e[:, 3:], np.isin(k, (1, 4)), **kw)
    return off, res


def test_contacts_match_the_published_geocentric_times(horizons):
    """Espenak's geocentric contacts of 2016 and 2019 (UT1), at the IERS ΔT: the model's
    contacts fall within 2 s of the printed times."""
    off, res = _offsets(horizons)
    assert np.all(res < CR.ROOT_TOL)
    assert np.abs(off).max() < TOL_CONTACT_S, dict(
        zip(horizons["labels"], off.round(2), strict=True)
    )


def test_without_light_time_the_contacts_are_far_off(horizons):
    off, _ = _offsets(horizons, light_time=False)
    assert np.abs(off).max() > 30.0


# ---------------------------------------------------------------------------- the solver


def test_light_time_delays_a_uniformly_moving_planet(monkeypatch):
    """Mercury in uniform motion before a fixed Sun and observer: light time moves each contact
    later by τ (to second order in the speed)."""
    monkeypatch.setattr(CR, "MU_SUN", 0.0)
    xo, vo = np.array([[1.0, 0.0, 0.0]]), np.zeros((1, 3))
    for exterior in (False, True):
        xm = np.array([[0.4, -0.003, 0.0005]])
        vm = np.array([[0.0, 0.02, 0.0]])
        a, _ = CR.contact_offsets(xm, vm, xo, vo, [exterior], light_time=False, aberration=False)
        b, _ = CR.contact_offsets(xm, vm, xo, vo, [exterior], aberration=False)
        tau = 0.6 / CR.C_AU_PER_DAY * 86400.0
        assert b[0] - a[0] == pytest.approx(tau, abs=0.05)


def test_aberration_barely_moves_the_separation(horizons):
    on, _ = _offsets(horizons)
    off, _ = _offsets(horizons, aberration=False)
    assert np.abs(on - off).max() < 3.0  # (v/c) × separation, under 0.1″
