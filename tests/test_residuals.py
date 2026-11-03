"""M4 / A2: fit Uranus alone (no perturber) to synthetic data; inspect the residual curve."""

import math

import numpy as np
import pytest

from fork.data.synthetic import (
    LongitudeObservations,
    jd_tdb_jan1,
    make_synthetic_observations,
    synthetic_epochs,
    true_uranus_lon_ecl_j2000_rad,
)
from fork.inference.fit import (
    fit_uranus_only,
    initial_uranus_theta,
    uranus_theta_to_x,
    uranus_x_to_theta,
    wrap_rad,
)
from fork.model.forward import ForwardContext, UranusModel

ARC = math.degrees(1.0) * 3600.0
EPOCH = 2371628.5
YEARS = np.arange(1690, 1846)

pytestmark = pytest.mark.usefixtures("offline_horizons")


@pytest.fixture(scope="module")
def model(offline_horizons):
    return UranusModel(ForwardContext.from_horizons(EPOCH))


def residual_curve(model, theta):
    """Noise-free residual curve: Horizons truth − fitted model, yearly 1690–1845 [″]."""
    grid = jd_tdb_jan1(YEARS)
    return wrap_rad(true_uranus_lon_ecl_j2000_rad(grid) - model.predict(theta, grid)) * ARC


def slope_sign_change_years(curve):
    s = np.sign(np.diff(curve))
    idx = np.nonzero(s[1:] != s[:-1])[0]
    return YEARS[1:-1][idx] + 0.5


def run_a2(model, seed):
    obs, _truth = make_synthetic_observations(np.random.default_rng(seed))
    fit = fit_uranus_only(obs, model)  # all 6 Uranus elements free, data-only start
    curve = residual_curve(model, fit.theta)
    return fit, curve, slope_sign_change_years(curve)


def test_a2_residual_curve(model):
    """A2: peak-to-peak ≈ 130″ (±25″) and a slope sign change within 1815–1830.

    Measured (seed 0): p2p 142.6″; slope sign changes at 1738.5, 1770.5, 1791.5, 1804.5,
    1818.5, 1833.5 → 1818.5 is the one in the window. See docs/notes/M4.md for why there are
    several.
    """
    fit, curve, flips = run_a2(model, seed=0)
    assert fit.success
    assert 105.0 <= np.ptp(curve) <= 155.0
    assert np.any((flips >= 1815) & (flips <= 1830))


@pytest.mark.slow
@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_a2_robust_across_noise_seeds(model, seed):
    """Measured p2p 137.8–141.6″, flip at 1818.5 for seeds 1–4."""
    fit, curve, flips = run_a2(model, seed)
    assert fit.success
    assert 105.0 <= np.ptp(curve) <= 155.0
    assert np.any((flips >= 1815) & (flips <= 1830))


def test_fit_recovers_model_generated_truth_exactly(model):
    """Self-consistency: data generated *by the model* (no Neptune) with no noise are fitted back
    to the generating elements — including i and Ω (up to Ω ≡ Ω + π) — from the data-only start."""
    theta_true = np.array([19.21, 0.046, math.radians(0.77), math.radians(74.0),
                           math.radians(170.0), math.radians(96.0)])  # fmt: skip
    jd, sigma, source = synthetic_epochs()
    obs = LongitudeObservations(jd, model.predict(theta_true, jd), sigma / ARC, source)
    fit = fit_uranus_only(obs, model)
    assert fit.rms_arcsec < 1e-4
    a, e, inc, raan, varpi, lam0 = fit.theta
    assert abs(a - theta_true[0]) < 1e-6 and abs(e - theta_true[1]) < 1e-7
    assert (
        abs(wrap_rad(varpi - theta_true[4])) < 1e-6 and abs(wrap_rad(lam0 - theta_true[5])) < 1e-8
    )
    assert abs(inc - theta_true[2]) < 1e-4
    assert (
        min(abs(wrap_rad(raan - theta_true[3])), abs(wrap_rad(raan - theta_true[3] - math.pi)))
        < 1e-2
    )


def test_initial_guess_is_data_only(model):
    """The starting point is built from the observations alone and is circular and planar."""
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    th0 = initial_uranus_theta(obs, EPOCH)
    assert th0[1] == 0.0 and th0[2] == 0.0 and th0[3] == 0.0
    assert 18.5 < th0[0] < 20.0


def test_x_theta_round_trip():
    th = np.array([19.2, 0.047, 0.0135, 1.29, 3.02, 1.67])
    np.testing.assert_allclose(uranus_x_to_theta(uranus_theta_to_x(th)), th, rtol=0, atol=1e-14)
    planar = uranus_x_to_theta(np.array([19.2, 0.0, 0.0, 1.0, 0.0, 0.0]))
    assert planar[1] == 0.0 and planar[2] == 0.0 and planar[3] == 0.0
