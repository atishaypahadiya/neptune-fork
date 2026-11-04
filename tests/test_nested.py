"""dynesty wrapper (``fork.inference.nested``): CI smoke test only; real runs are local."""

import math

import numpy as np

from fork.data.synthetic import make_synthetic_observations
from fork.inference.fit import uranus_theta_to_x
from fork.inference.nested import NestedResult, a_au_of, mode_masses, run_nested
from fork.inference.priors import UranusBox, free_a_space
from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

EPOCH = 2371628.5


def test_nested_smoke_run_and_mode_bookkeeping(offline_horizons, tmp_path):
    """A few iterations on seed 0 with the real prior and likelihood: runs serially, every
    sample lies in the prior support, the result round-trips through the cache, and the mode
    split / jitter bookkeeping returns a mass in [0, 1]."""
    ctx = ForwardContext.from_horizons(EPOCH)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    space = free_a_space(pq, 0.3, box)
    res = run_nested(obs, ctx, space, np.random.default_rng(1), nlive=12, maxiter=15)
    assert res.niter >= 15 and math.isfinite(res.logz) and res.logzerr > 0
    assert all(space.log_prior(q) == 0.0 for q in res.samples)
    path = tmp_path / "r.pkl"
    res.save(path)
    again = NestedResult.load(path)
    assert np.array_equal(again.samples, res.samples) and again.logz == res.logz
    m = mode_masses(res, lambda s: a_au_of(space, s) < 25.0, np.random.default_rng(0), n_jitter=3)
    assert 0.0 <= m["mass"] <= 1.0 and m["logz_jitter_std"] >= 0.0
    assert abs(res.weights.sum() - 1.0) < 1e-9
