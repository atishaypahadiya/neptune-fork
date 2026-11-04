"""M5b: dynesty nested sampling of the hidden-mass posterior (and its evidence).

Same q, prior and likelihood as the emcee sampler (``fork.inference.priors``, ``mcmc.chi2``);
the prior enters through ``PosteriorSpace.prior_transform``. Why nested sampling: on seed 0 the
free-a posterior has two modes (a ≈ 30 AU with e ≲ 0.05, and a ≈ 22 AU with e ≈ 0.28 against
the e ≤ 0.3 bound) that emcee does not move between, so each chain reports the mode it started
in. Nested sampling with multi-ellipsoid bounds explores the whole prior, weighs the modes, and
returns ln Z, which Phase 2's model comparison needs anyway.

Settings (PHASE1_SPEC M5b): bound='multi', sample='rslice' (dynesty's default 3 + ndim slices),
nlive=1000, stop at Δln Z < 0.1; the perturber's ϖ and λ₀ are periodic in the unit cube.
"""

from __future__ import annotations

import math
import pickle
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.data.synthetic import LongitudeObservations
from fork.inference.mcmc import chi2
from fork.inference.priors import PosteriorSpace
from fork.model.forward import ForwardContext

NLIVE = 1000
DLOGZ = 0.1
BOUND = "multi"
SAMPLE = "rslice"

# per-process problem, set once by ``_init_worker`` so dynesty's pool never pickles the data
_PROBLEM: tuple | None = None


def _init_worker(space: PosteriorSpace, obs: LongitudeObservations, ctx: ForwardContext) -> None:
    global _PROBLEM
    _PROBLEM = (space, obs, ctx)


def _log_likelihood(q: np.ndarray) -> float:
    space, obs, ctx = _PROBLEM  # type: ignore[misc]
    return -0.5 * chi2(q, space, obs, ctx)


def _prior_transform(u: np.ndarray) -> np.ndarray:
    return _PROBLEM[0].prior_transform(u)  # type: ignore[index]


@dataclass(frozen=True)
class NestedResult:
    space: PosteriorSpace
    results: object  # dynesty Results (samples in q units, logwt, logl, logz, logzerr, ...)
    nlive: int
    dlogz: float
    wall_s: float
    rng_seed: int | None = None

    @property
    def samples(self) -> np.ndarray:
        return np.asarray(self.results["samples"])

    @property
    def weights(self) -> np.ndarray:
        return np.asarray(self.results.importance_weights())

    @property
    def logz(self) -> float:
        return float(self.results["logz"][-1])

    @property
    def logzerr(self) -> float:
        return float(self.results["logzerr"][-1])

    @property
    def niter(self) -> int:
        return int(self.results["niter"])

    @property
    def ncall(self) -> int:
        return int(np.sum(self.results["ncall"]))

    @property
    def max_logl(self) -> float:
        return float(np.max(self.results["logl"]))

    def equal_weight(self, rng: np.random.Generator) -> np.ndarray:
        from dynesty.utils import resample_equal

        return resample_equal(self.samples, self.weights, rstate=rng)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".tmp")
        with tmp.open("wb") as f:
            pickle.dump(self, f)
        tmp.replace(path)

    @staticmethod
    def load(path: Path) -> NestedResult:
        # pickle: only for this project's own cache (data/cache/nested/, gitignored), written
        # by ``save`` above; never point this at files from elsewhere.
        with path.open("rb") as f:
            return pickle.load(f)


def run_nested(
    obs: LongitudeObservations,
    context: ForwardContext,
    space: PosteriorSpace,
    rng: np.random.Generator,
    *,
    nlive: int = NLIVE,
    dlogz: float = DLOGZ,
    n_workers: int = 1,
    maxiter: int | None = None,
    print_progress: bool = False,
) -> NestedResult:
    """Static nested sampling over ``space``'s prior (which must carry a Uranus box)."""
    from multiprocessing import get_context

    import dynesty

    space._check_transformable()
    _init_worker(space, obs, context)
    kw = dict(nlive=nlive, bound=BOUND, sample=SAMPLE, periodic=space.periodic_unit_dims,
              rstate=rng)  # fmt: skip
    t0 = time.time()
    if n_workers > 1:
        with get_context("spawn").Pool(n_workers, initializer=_init_worker,
                                       initargs=(space, obs, context)) as pool:  # fmt: skip
            sampler = dynesty.NestedSampler(_log_likelihood, _prior_transform, space.ndim,
                                            pool=pool, queue_size=n_workers, **kw)  # fmt: skip
            sampler.run_nested(dlogz=dlogz, maxiter=maxiter, print_progress=print_progress)
    else:
        sampler = dynesty.NestedSampler(_log_likelihood, _prior_transform, space.ndim, **kw)
        sampler.run_nested(dlogz=dlogz, maxiter=maxiter, print_progress=print_progress)
    return NestedResult(space, sampler.results, nlive, dlogz, time.time() - t0)


def weighted_quantiles(x: np.ndarray, w: np.ndarray, q=(0.16, 0.5, 0.84)) -> list[float]:
    from dynesty.utils import quantile

    return [float(v) for v in quantile(np.asarray(x), list(q), weights=np.asarray(w))]


def mode_masses(res: NestedResult, in_mode, rng: np.random.Generator, n_jitter: int = 100
                ) -> dict:  # fmt: skip
    """Posterior mass of the samples selected by ``in_mode(samples) -> bool array``, with the
    spread over ``n_jitter`` jittered runs (dynesty's resampling of the prior-volume shrinkage;
    the statistical error of one run). Also ln Z's jitter spread."""
    from dynesty.utils import jitter_run

    mask = np.asarray(in_mode(res.samples), dtype=bool)
    base = float(res.weights[mask].sum())
    fr, lz = [], []
    for _ in range(n_jitter):
        j = jitter_run(res.results, rstate=rng)
        w = np.exp(j["logwt"] - j["logz"][-1])
        fr.append(float(w[mask].sum() / w.sum()))
        lz.append(float(j["logz"][-1]))
    return {"mass": base, "mass_jitter_std": float(np.std(fr)),
            "logz_jitter_std": float(np.std(lz))}  # fmt: skip


def a_au_of(space: PosteriorSpace, samples: np.ndarray) -> np.ndarray:
    """Perturber a [AU] of q samples (free-a spaces: q[5] = ln a)."""
    if not space.free[1]:
        return np.full(len(samples), math.exp(space.held[1]))  # type: ignore[arg-type]
    return np.exp(samples[:, 4 + int(space.free[0])])
