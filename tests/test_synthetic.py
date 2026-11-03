"""M4: synthetic observations."""

import math

import numpy as np
import pytest

from fork.data import synthetic
from fork.inference.fit import wrap_rad

ARC = math.degrees(1.0) * 3600.0
pytestmark = pytest.mark.usefixtures("offline_horizons")


def test_schedule():
    jd, sigma, source = synthetic.synthetic_epochs()
    pre = source == "synthetic_pre_discovery"
    assert pre.sum() == 22 and (~pre).sum() == 65
    assert np.all(sigma[pre] == 12.0) and np.all(sigma[~pre] == 3.0)
    assert np.all(np.diff(jd) >= 0)
    years = 2000.0 + (jd - 2451545.0) / 365.25
    assert years[pre].min() < 1691 and years[pre].max() < 1772
    assert 1780.9 < years[~pre].min() and years[~pre].max() < 1846


def test_seeded_and_noise_level():
    a, truth = synthetic.make_synthetic_observations(np.random.default_rng(42))
    b, _ = synthetic.make_synthetic_observations(np.random.default_rng(42))
    c, _ = synthetic.make_synthetic_observations(np.random.default_rng(43))
    np.testing.assert_array_equal(a.lon_ecl_j2000_rad, b.lon_ecl_j2000_rad)
    assert not np.array_equal(a.lon_ecl_j2000_rad, c.lon_ecl_j2000_rad)
    z = wrap_rad(a.lon_ecl_j2000_rad - truth) / a.sigma_rad
    assert 0.6 < np.std(z) < 1.4 and abs(np.mean(z)) < 0.5


def test_truth_is_sun_centred_horizons_longitude():
    from fork.data.horizons import state_vectors_ecl_j2000

    jd = synthetic.jd_tdb_jan1([1781, 1822])
    st = state_vectors_ecl_j2000("uranus", jd)
    want = np.mod(np.arctan2(st[:, 1], st[:, 0]), 2 * math.pi)
    np.testing.assert_array_equal(synthetic.true_uranus_lon_ecl_j2000_rad(jd), want)


def test_observations_do_not_carry_truth():
    obs, _ = synthetic.make_synthetic_observations(np.random.default_rng(0))
    assert set(vars(obs)) == {"jd_tdb", "lon_ecl_j2000_rad", "sigma_rad", "source"}
