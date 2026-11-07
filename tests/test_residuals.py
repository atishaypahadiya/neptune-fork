"""M4 / A2, A2b: fit Uranus to synthetic data without and with the true Neptune planted.

Every fit frees a, e, ϖ, λ₀ from a data-only start and holds i, Ω at epoch values
(docs/RULES.md rule 2).
"""

import math

import numpy as np
import pytest

from fork.data.horizons import states_about_inner_barycenter_ecl_j2000
from fork.data.synthetic import (
    DISCOVERY_YEAR,
    SIGMA_POST_DISCOVERY_ARCSEC,
    LongitudeObservations,
    jd_tdb_jan1,
    make_synthetic_observations,
    synthetic_epochs,
    true_uranus_lon_ecl_j2000_rad,
)
from fork.inference.fit import (
    fit_uranus_elements,
    fit_uranus_only,
    initial_uranus_theta,
    uranus_theta_to_x,
    uranus_x_to_theta,
    wrap_rad,
)
from fork.model import bodies
from fork.model.forward import (
    ForwardContext,
    UranusModel,
    perturber_theta_from_state,
    uranus_inc_raan_at_epoch,
)

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


def slope_sign_change_years(curve, years=YEARS):
    s = np.sign(np.diff(curve))
    idx = np.nonzero(s[1:] != s[:-1])[0]
    return years[1:-1][idx] + 0.5


def run_a2(model, seed):
    obs, _truth = make_synthetic_observations(np.random.default_rng(seed))
    fit = fit_uranus_only(obs, model, uranus_inc_raan_at_epoch(EPOCH))
    curve = residual_curve(model, fit.theta)
    post = YEARS >= DISCOVERY_YEAR
    return fit, curve, slope_sign_change_years(curve[post], YEARS[post])


def assert_a2(fit, curve, post_flips):
    assert fit.success
    assert 105.0 <= np.ptp(curve) <= 155.0  # peak-to-peak over the whole 1690–1845 curve
    assert np.any((post_flips >= 1815) & (post_flips <= 1830))  # post-1781 part only


def test_a2_residual_curve(model):
    """A2: peak-to-peak ≈ 130″ (±25″); the post-1781 residual slope changes sign in 1815–1830.

    Measured (seed 0, dt = 12 d): p2p 142.7″; post-1781 slope sign changes at 1783.5, 1806.5,
    **1828.5**. See docs/notes/M4.md.
    """
    assert_a2(*run_a2(model, seed=0))


@pytest.mark.slow
@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_a2_robust_across_noise_seeds(model, seed):
    """Measured p2p 142.7–145.4″; post-1781 in-window flip at 1828.5 for seeds 1–4."""
    assert_a2(*run_a2(model, seed))


@pytest.fixture(scope="module")
def neptune_theta(offline_horizons):
    """True Neptune at the epoch as a planar perturber θ (test only: A2b plants the truth)."""
    state = states_about_inner_barycenter_ecl_j2000("neptune", [EPOCH])[0]
    return perturber_theta_from_state(state, bodies.NEPTUNE.mass_msun)


def post_1781_rms_arcsec(obs, fit):
    post = obs.source == "synthetic_post_discovery"
    return float(np.sqrt(np.mean(fit.residual_rad[post] ** 2)) * ARC)


def run_a2b(model, neptune_theta, seed):
    obs, _ = make_synthetic_observations(np.random.default_rng(seed))
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    planted = fit_uranus_elements(obs, model.context, inc_raan, fixed_perturber_theta=neptune_theta)
    alone = fit_uranus_only(obs, model, inc_raan)
    return post_1781_rms_arcsec(obs, planted), post_1781_rms_arcsec(obs, alone), planted


def test_a2b_planted_neptune_fits_to_noise(model, neptune_theta):
    """A2b: with the true Neptune planted (not fitted) and Uranus refitted, the post-1781
    residual RMS is within 1.5 σ_noise (4.5″); without Neptune the same check fails.

    Measured (seed 0): 2.85″ = 0.95σ planted; 13.65″ = 4.55σ without Neptune.
    """
    rms_planted, rms_alone, fit = run_a2b(model, neptune_theta, seed=0)
    assert fit.success
    assert rms_planted <= 1.5 * SIGMA_POST_DISCOVERY_ARCSEC
    assert rms_alone > 1.5 * SIGMA_POST_DISCOVERY_ARCSEC


@pytest.mark.slow
@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_a2b_robust_across_noise_seeds(model, neptune_theta, seed):
    """Measured planted post-1781 RMS 0.93–1.02σ for seeds 1–4."""
    rms_planted, rms_alone, _ = run_a2b(model, neptune_theta, seed)
    assert rms_planted <= 1.5 * SIGMA_POST_DISCOVERY_ARCSEC
    assert rms_alone > 1.5 * SIGMA_POST_DISCOVERY_ARCSEC


def test_fixed_perturber_validation(model):
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    with pytest.raises(ValueError):
        fit_uranus_elements(obs, model.context, (0.0, 0.0), fixed_perturber_theta=np.zeros(4))


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
