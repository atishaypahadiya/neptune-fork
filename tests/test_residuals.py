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
from fork.model.forward import ForwardContext, UranusModel, uranus_inc_raan_at_epoch

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
    # a, e, ϖ, λ₀ free from a data-only start; i, Ω held at epoch values (project rule 2)
    fit = fit_uranus_only(obs, model, uranus_inc_raan_at_epoch(EPOCH))
    curve = residual_curve(model, fit.theta)
    return fit, curve, slope_sign_change_years(curve)


def test_a2_residual_curve(model):
    """A2: peak-to-peak ≈ 130″ (±25″) and a slope sign change within 1815–1830.

    Measured (seed 0, dt = 14 d, i and Ω held): p2p 142.7″; slope sign changes at 1698.5,
    1734.5, 1783.5, 1806.5, 1828.5. See docs/notes/M4.md.
    """
    fit, curve, flips = run_a2(model, seed=0)
    assert fit.success
    assert 105.0 <= np.ptp(curve) <= 155.0
    assert np.any((flips >= 1815) & (flips <= 1830))


@pytest.mark.slow
@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_a2_robust_across_noise_seeds(model, seed):
    """Measured (dt = 14 d) p2p 142.7–145.4″, in-window flip at 1828.5 for seeds 1–4."""
    fit, curve, flips = run_a2(model, seed)
    assert fit.success
    assert 105.0 <= np.ptp(curve) <= 155.0
    assert np.any((flips >= 1815) & (flips <= 1830))


THETA_MODEL_TRUTH = np.array([19.21, 0.046, math.radians(0.77), math.radians(74.0),
                              math.radians(170.0), math.radians(96.0)])  # fmt: skip


def _model_generated_obs(model):
    jd, sigma, source = synthetic_epochs()
    return LongitudeObservations(jd, model.predict(THETA_MODEL_TRUTH, jd), sigma / ARC, source)


def test_fit_recovers_model_generated_truth_exactly(model):
    """Self-consistency: noise-free data generated *by the model* (no Neptune) are fitted back to
    the generating in-plane elements from the data-only start, with i, Ω held at their truth."""
    th = THETA_MODEL_TRUTH
    fit = fit_uranus_only(_model_generated_obs(model), model, (th[2], th[3]))
    assert fit.rms_arcsec < 1e-4
    a, e, inc, raan, varpi, lam0 = fit.theta
    assert (inc, raan) == (th[2], th[3])  # held exactly
    assert abs(a - th[0]) < 1e-6 and abs(e - th[1]) < 1e-7
    assert abs(wrap_rad(varpi - th[4])) < 1e-6 and abs(wrap_rad(lam0 - th[5])) < 1e-8


def test_free_inc_raan_diagnostic_recovers_model_truth(model):
    """The diagnostic all-six-free mode (data-only multistart) also recovers i and Ω, up to
    Ω ≡ Ω + π, when there is no unmodelled perturber to absorb."""
    th = THETA_MODEL_TRUTH
    fit = fit_uranus_only(_model_generated_obs(model), model, (0.0, 0.0), free_inc_raan=True)
    assert fit.rms_arcsec < 1e-4
    inc, raan = fit.theta[2], fit.theta[3]
    assert abs(inc - th[2]) < 1e-4
    assert min(abs(wrap_rad(raan - th[3])), abs(wrap_rad(raan - th[3] - math.pi))) < 1e-2


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
