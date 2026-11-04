"""dynesty wrapper (``fork.inference.nested``): CI smoke test only; real runs are local."""

import math

import numpy as np
import pytest

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


class _Crash(RuntimeError):
    pass


class _CrashingLikelihood:
    """Module-level (so dynesty's checkpoint can pickle it): raises once ``calls`` exceeds the
    module flag ``_CRASH_AFTER``; harmless once the flag is None."""

    def __init__(self):
        self.calls = 0

    def __call__(self, q):
        import fork.inference.nested as nested
        from fork.inference.mcmc import chi2

        self.calls += 1
        if _CRASH_AFTER is not None and self.calls > _CRASH_AFTER:
            raise _Crash
        space, obs, ctx = nested._PROBLEM
        return -0.5 * chi2(q, space, obs, ctx)


_CRASH_AFTER: int | None = None


@pytest.mark.slow
def test_nested_checkpoint_resume_and_progress_log(offline_horizons, tmp_path, monkeypatch):
    """A run that dies mid-loop (here: the likelihood raises) leaves a checkpoint; calling again
    with the same file resumes from it (iterations continue rather than restart) and finishes.
    The progress log has start, periodic status, resume and done lines. (A run that *ends*,
    e.g. by maxiter, is finished and does not resume; only interruptions do.) ~15 s."""
    global _CRASH_AFTER
    import fork.inference.nested as nested

    ctx = ForwardContext.from_horizons(EPOCH)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    space = free_a_space(pq, 0.3, box)
    ckpt, log = tmp_path / "run.ckpt", tmp_path / "run.log"
    kw = dict(nlive=12, checkpoint_file=ckpt, checkpoint_every_s=0.2, progress_log=log,
              progress_every_s=0.0)  # fmt: skip
    monkeypatch.setattr(nested, "_log_likelihood", _CrashingLikelihood())
    _CRASH_AFTER = 300
    try:
        with pytest.raises(_Crash):
            run_nested(obs, ctx, space, np.random.default_rng(1), **kw)
    finally:
        _CRASH_AFTER = None
    assert ckpt.exists()
    import dynesty

    it_saved = dynesty.NestedSampler.restore(str(ckpt)).it
    assert it_saved > 1
    res = run_nested(obs, ctx, space, np.random.default_rng(99), maxiter=it_saved + 5, **kw)
    assert res.niter > it_saved
    text = log.read_text()
    assert "start:" in text and "resumed from" in text and "iter " in text and "done:" in text
