"""M5b sampler: CI smoke tests only. Real chains run locally (``python -m fork.export.m5b``)."""

import math

import numpy as np
import pytest

from fork.data.synthetic import LongitudeObservations, synthetic_epochs
from fork.inference import mcmc
from fork.model.forward import ForwardContext, predict_uranus_lon

ARC = math.degrees(1.0) * 3600.0
EPOCH = 2371628.5
THETA_MODEL = np.array([19.21, 0.046, math.radians(0.77), math.radians(74.0), math.radians(170.0),
                        math.radians(96.0), math.log10(20.0), 32.0, 0.05, 1.0, 3.3])  # fmt: skip
PQ = (math.tan(THETA_MODEL[2] / 2) * math.cos(THETA_MODEL[3]),
      math.tan(THETA_MODEL[2] / 2) * math.sin(THETA_MODEL[3]))  # fmt: skip

pytestmark = pytest.mark.usefixtures("offline_horizons")


@pytest.fixture(scope="module")
def ctx(offline_horizons):
    return ForwardContext.from_horizons(EPOCH)


@pytest.fixture(scope="module")
def model_obs(ctx):
    jd, sigma, source = synthetic_epochs()
    return LongitudeObservations(jd, predict_uranus_lon(THETA_MODEL, jd, ctx), sigma / ARC, source)


def test_space_transforms_and_priors():
    sp = mcmc.free_a_space(PQ)
    assert sp.ndim == 9 and sp.names[5] == "p_ln_a_au"
    q = sp.q_from_theta(THETA_MODEL)
    assert sp.theta(q) == pytest.approx(THETA_MODEL, abs=1e-12)
    assert sp.log_prior(q) == 0.0
    for i, v in ((4, 2.1), (4, -0.1), (5, math.log(19.0)), (5, math.log(61.0)), (6, 0.56)):
        bad = q.copy()
        bad[i] = v
        assert sp.log_prior(bad) == -math.inf  # m, a outside the prior; e = 0.3136 > 0.3
    tight = mcmc.free_a_space(PQ, e_max=0.1)
    q2 = q.copy()
    q2[6] = 0.33  # e = 0.1089 + h²
    assert tight.log_prior(q2) == -math.inf and sp.log_prior(q2) == 0.0
    lv = mcmc.held_shape_space(PQ, 36.15, 0.10761, math.radians(286.9))
    assert lv.ndim == 6 and lv.names[4:] == ("p_log10_m_earth", "p_lam0_ecl_j2000_rad")
    th = lv.theta(np.concatenate([q[:4], [1.4, 2.0]]))
    assert th[7:10] == pytest.approx([36.15, 0.10761, math.radians(286.9)])


def test_circular_interval_wraps():
    x = np.array([355.0, 357.0, 359.0, 1.0, 3.0, 5.0])
    lo, med, hi = mcmc.circular_interval_deg(x)

    def signed(v):
        return (v + 180.0) % 360.0 - 180.0

    assert signed(med) == pytest.approx(0.0, abs=1e-9)
    assert -5.0 < signed(lo) < -1.0 and 1.0 < signed(hi) < 5.0  # not the 1°…359° wrong way


def test_tiny_chain_smoke(ctx, model_obs):
    """2 × (10 steps × 20 walkers), ~2 s: shapes, prior respected, finite log-prob, reproducible."""
    sp = mcmc.free_a_space(PQ)
    q0 = sp.q_from_theta(THETA_MODEL)
    runs = [mcmc.run_chain(model_obs, ctx, sp, q0, np.random.default_rng(7), n_walkers=20,
                           n_steps=10, max_steps=10, burn=2, thin=1) for _ in range(2)]  # fmt: skip
    r = runs[0]
    assert r.samples.shape == (8 * 20, 9) and r.n_steps == 10 and not r.converged
    assert np.all(np.isfinite(r.log_prob))
    assert all(sp.log_prior(q) == 0.0 for q in r.samples)
    assert (
        r.max_log_prob > -50.0
    )  # ball around the χ² = 0 truth (0.02 in m, a, k, h); measured −8.4
    assert np.array_equal(r.samples, runs[1].samples)
