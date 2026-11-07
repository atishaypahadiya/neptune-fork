"""M5b: emcee posteriors for the hidden-mass model (PHASE1_SPEC M5).

Sampled vector q (docs/RULES.md rule 2: Uranus's in-plane four always free, i and Ω held):

    Uranus   (a, e·cos ϖ, e·sin ϖ, λ₀)
    perturber (log10 m [M⊕], ln a, √e·cos ϖ, √e·sin ϖ, λ₀)   — any of these may be held

Priors: ``fork.inference.priors`` (the one definition, shared with the nested sampler). λ₀ is
sampled on the real line here (periodic in the likelihood). Chains run before the Uranus box
existed used flat Uranus priors with an e < 0.95 guard; with a box set, emcee and dynesty
sample the identical prior. Likelihood: Gaussian, ln L = −χ²/2 with the observation σ.

Moves: a mix of differential-evolution moves (80 % DE, 20 % DE-snooker), which follow the
curved mass–distance–eccentricity ridge better than the default stretch move.
Walkers start in a small ball around a bounded multistart point fit. Chains run until n ≥ 50 τ
(spec) and the first 10 τ_max are discarded as burn-in.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from fork.data.synthetic import LongitudeObservations
from fork.inference.fit import wrap_rad
from fork.inference.priors import (  # the single prior definition, shared with ``nested``
    PosteriorSpace,
    UranusBox,
    free_a_space,
    held_shape_space,
)
from fork.model.forward import ForwardContext, predict_uranus_lon

__all__ = ["PosteriorSpace", "UranusBox", "free_a_space", "held_shape_space", "run_chain"]


def chi2(q: np.ndarray, space: PosteriorSpace, obs: LongitudeObservations,
         context: ForwardContext) -> float:  # fmt: skip
    if not isinstance(obs, LongitudeObservations):  # M6: RaDecObservations (own χ²)
        return obs.chi2(space.theta(q), context)
    pred = predict_uranus_lon(space.theta(q), obs.jd_tdb, context)
    return float(np.sum((wrap_rad(obs.lon_ecl_j2000_rad - pred) / obs.sigma_rad) ** 2))


def log_prob(q: np.ndarray, space: PosteriorSpace, obs: LongitudeObservations,
             context: ForwardContext) -> float:  # fmt: skip
    lp = space.log_prior(q)
    if not math.isfinite(lp):
        return -math.inf
    return lp - 0.5 * chi2(q, space, obs, context)


def initial_walkers(q0: np.ndarray, space: PosteriorSpace, n_walkers: int,
                    rng: np.random.Generator) -> np.ndarray:  # fmt: skip
    """Small ball around q0, inside the prior. Scales: Uranus 1e-4 (a, k, h), 1e-5 rad (λ₀);
    perturber 0.02 in log10 m, ln a, k, h and 0.01 rad in λ₀."""
    scale_p = np.array([0.02, 0.02, 0.02, 0.02, 0.01])[space.free]
    scale = np.concatenate([[1e-4, 1e-4, 1e-4, 1e-5], scale_p])
    out = np.empty((n_walkers, q0.size))
    k = 0
    while k < n_walkers:
        q = q0 + scale * rng.standard_normal(q0.size)
        if math.isfinite(space.log_prior(q)):
            out[k] = q
            k += 1
    return out


@dataclass(frozen=True)
class ChainResult:
    space: PosteriorSpace
    samples: np.ndarray  # (n, ndim) post burn-in, thinned, flattened over walkers
    log_prob: np.ndarray  # (n,)
    tau: np.ndarray  # integrated autocorrelation time per parameter [steps]
    n_steps: int
    n_walkers: int
    burn: int
    thin: int
    acceptance: float
    converged: bool  # n_steps ≥ TAU_RATIO · τ_max
    max_log_prob: float  # over the whole chain (incl. burn-in)
    q_max_log_prob: np.ndarray
    chain: np.ndarray | None = None  # (n_steps, n_walkers, ndim) raw, if kept


TAU_RATIO = 50.0  # PHASE1_SPEC: run ≥ 50 τ
BURN_TAU = 10.0  # discard the first 10 τ_max (walkers start in a small ball)


def run_chain(
    obs: LongitudeObservations,
    context: ForwardContext,
    space: PosteriorSpace,
    q0: np.ndarray,
    rng: np.random.Generator,
    *,
    n_walkers: int = 48,
    n_steps: int = 20000,
    max_steps: int = 60000,
    chunk: int = 5000,
    burn: int | None = None,
    thin: int | None = None,
    pool=None,
    keep_chain: bool = False,
    progress: bool = False,
) -> ChainResult:
    """Run emcee from a ball around ``q0`` for ``n_steps``, then extend in ``chunk`` steps until
    n ≥ 50 τ_max or ``max_steps``. Default burn = 10 τ_max, thin = τ_max / 2."""
    import emcee

    p0 = initial_walkers(q0, space, n_walkers, rng)
    moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    sampler = emcee.EnsembleSampler(n_walkers, space.ndim, log_prob, args=(space, obs, context),
                                    moves=moves, pool=pool)  # fmt: skip
    # emcee draws its moves from a legacy RandomState; seed it from rng for reproducibility
    sampler.random_state = np.random.RandomState(int(rng.integers(2**31))).get_state()
    sampler.run_mcmc(p0, n_steps, progress=progress)
    tau = sampler.get_autocorr_time(tol=0)
    while sampler.iteration < TAU_RATIO * np.max(tau) and sampler.iteration < max_steps:
        sampler.run_mcmc(None, min(chunk, max_steps - sampler.iteration), progress=progress)
        tau = sampler.get_autocorr_time(tol=0)
    t_max = float(np.max(tau))
    n = sampler.iteration
    # never discard more than half the chain (only matters for chains that hit max_steps)
    burn = min(math.ceil(BURN_TAU * t_max), n // 2) if burn is None else burn
    thin = max(1, min(int(t_max // 2), (n - burn) // 4)) if thin is None else thin
    lp_all = sampler.get_log_prob()
    i_max = np.unravel_index(np.argmax(lp_all), lp_all.shape)
    return ChainResult(
        space=space,
        samples=sampler.get_chain(discard=burn, thin=thin, flat=True),
        log_prob=sampler.get_log_prob(discard=burn, thin=thin, flat=True),
        tau=tau, n_steps=n, n_walkers=n_walkers, burn=burn, thin=thin,
        acceptance=float(np.mean(sampler.acceptance_fraction)),
        converged=bool(n >= TAU_RATIO * t_max),
        max_log_prob=float(lp_all[i_max]),
        q_max_log_prob=sampler.get_chain()[i_max].copy(),
        chain=sampler.get_chain() if keep_chain else None,
    )  # fmt: skip


def circular_interval_deg(lon_deg: np.ndarray, q: tuple[float, ...] = (16.0, 50.0, 84.0)
                          ) -> tuple[float, ...]:  # fmt: skip
    """Percentiles of angles [deg], unwrapped around their circular mean (results in [0, 360))."""
    z = np.exp(1j * np.radians(lon_deg))
    centre = math.degrees(math.atan2(z.imag.mean(), z.real.mean()))
    rel = (np.asarray(lon_deg) - centre + 180.0) % 360.0 - 180.0
    return tuple(float((centre + v) % 360.0) for v in np.percentile(rel, q))
