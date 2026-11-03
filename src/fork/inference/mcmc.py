"""M5b: emcee posteriors for the hidden-mass model (PHASE1_SPEC M5).

Sampled vector q (CLAUDE.md rule 2: Uranus's in-plane four always free, i and Ω held):

    Uranus   (a, e·cos ϖ, e·sin ϖ, λ₀)
    perturber (log10 m [M⊕], ln a, √e·cos ϖ, √e·sin ϖ, λ₀)   — any of these may be held

Priors: log10 m uniform on [0, 2] (1–100 M⊕); ln a uniform with a ∈ [20, 60] AU (i.e. log-uniform
in a, the spec's "sample in log a"); (k, h) = √e·(cos ϖ, sin ϖ) uniform on the disk e ≤ e_max,
which is uniform in e and ϖ; λ₀ uniform (sampled on the real line, periodic in the likelihood).
Uranus's elements have flat priors, with e < 0.95 as a physical guard. Likelihood: Gaussian,
ln L = −χ²/2 with the observation σ.

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
from fork.inference.fit import E_PRIOR_MAX, TWO_PI, uranus_x_to_theta, wrap_rad
from fork.model.forward import ForwardContext, predict_uranus_lon

PRIOR_LOG10_M_EARTH = (0.0, 2.0)
PRIOR_A_AU = (20.0, 60.0)
_URANUS_MAX_E = 0.95
PERTURBER_Q_NAMES = ("p_log10_m_earth", "p_ln_a_au", "p_sqrt_e_cos_varpi", "p_sqrt_e_sin_varpi",
                     "p_lam0_ecl_j2000_rad")  # fmt: skip
URANUS_Q_NAMES = ("u_a_au", "u_e_cos_varpi", "u_e_sin_varpi", "u_lam0_ecl_j2000_rad")


@dataclass(frozen=True)
class PosteriorSpace:
    """Which perturber entries of q are sampled; held ones carry their value (in q units)."""

    held: tuple[float | None, ...] = (None, None, None, None, None)
    e_max: float = E_PRIOR_MAX
    pq: tuple[float, float] = (0.0, 0.0)  # Uranus (p, q) = tan(i/2)(cos Ω, sin Ω), held

    @property
    def free(self) -> np.ndarray:
        return np.array([v is None for v in self.held])

    @property
    def ndim(self) -> int:
        return 4 + int(self.free.sum())

    @property
    def names(self) -> tuple[str, ...]:
        return URANUS_Q_NAMES + tuple(n for n, f in zip(PERTURBER_Q_NAMES, self.free,
                                                         strict=True) if f)  # fmt: skip

    def full_perturber_q(self, q: np.ndarray) -> np.ndarray:
        y = np.array([0.0 if v is None else v for v in self.held])
        y[self.free] = q[4:]
        return y

    def theta(self, q: np.ndarray) -> np.ndarray:
        """q → 11-entry forward-model θ."""
        theta_u = uranus_x_to_theta(np.concatenate([q[:4], self.pq]))
        log10_m, ln_a, k, h, lam0 = self.full_perturber_q(q)
        return np.concatenate([theta_u, [log10_m, math.exp(ln_a), k * k + h * h,
                                         math.atan2(h, k) % TWO_PI, lam0 % TWO_PI]])  # fmt: skip

    def q_from_theta(self, theta: np.ndarray) -> np.ndarray:
        from fork.inference.fit import uranus_theta_to_x

        xu = uranus_theta_to_x(theta[:6])[:4]
        s = math.sqrt(theta[8])
        y = np.array([theta[6], math.log(theta[7]), s * math.cos(theta[9]),
                      s * math.sin(theta[9]), theta[10]])  # fmt: skip
        return np.concatenate([xu, y[self.free]])

    def log_prior(self, q: np.ndarray) -> float:
        log10_m, ln_a, k, h, _ = self.full_perturber_q(q)
        if self.free[0] and not PRIOR_LOG10_M_EARTH[0] <= log10_m <= PRIOR_LOG10_M_EARTH[1]:
            return -math.inf
        if self.free[1] and not math.log(PRIOR_A_AU[0]) <= ln_a <= math.log(PRIOR_A_AU[1]):
            return -math.inf
        if (self.free[2] or self.free[3]) and k * k + h * h > self.e_max:
            return -math.inf
        if math.hypot(q[1], q[2]) >= _URANUS_MAX_E or q[0] <= 0.0:
            return -math.inf
        return 0.0


def free_a_space(pq: tuple[float, float], e_max: float = E_PRIOR_MAX) -> PosteriorSpace:
    return PosteriorSpace(e_max=e_max, pq=pq)


def held_shape_space(pq: tuple[float, float], a_au: float, e: float, varpi_rad: float
                     ) -> PosteriorSpace:  # fmt: skip
    """a, e, ϖ held (e.g. Le Verrier's published shape); m, λ₀ and Uranus sampled."""
    s = math.sqrt(e)
    return PosteriorSpace((None, math.log(a_au), s * math.cos(varpi_rad), s * math.sin(varpi_rad),
                           None), e_max=max(e, E_PRIOR_MAX), pq=pq)  # fmt: skip


def chi2(q: np.ndarray, space: PosteriorSpace, obs: LongitudeObservations,
         context: ForwardContext) -> float:  # fmt: skip
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
