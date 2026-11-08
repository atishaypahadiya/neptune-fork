"""M5b: dynesty nested sampling of the hidden-mass posterior (and its evidence).

Same q, prior and likelihood as the emcee sampler (``fork.inference.priors``, ``mcmc.chi2``);
the prior enters through ``PosteriorSpace.prior_transform``. Why nested sampling: on seed 0 the
free-a posterior has two modes (a ≈ 30 AU with e ≲ 0.05, and a ≈ 22 AU with e ≈ 0.28 against
the e ≤ 0.3 bound) that emcee does not move between, so each chain reports the mode it started
in. Nested sampling with multi-ellipsoid bounds explores the whole prior, weighs the modes, and
returns ln Z, which Phase 2's model comparison needs anyway.

Settings (PHASE1_SPEC M5b): bound='multi', sample='rslice' (dynesty's default 3 + ndim slices),
nlive=1000, stop at Δln Z < 0.1; the perturber's ϖ and λ₀ are periodic in the unit cube.

Long runs are monitorable and resumable: ``progress_log`` gets one status line per
``progress_every_s`` (iteration, calls, ln Z, Δln Z remaining, elapsed), and
``checkpoint_file`` is saved every ``checkpoint_every_s``; if it exists when ``run_nested`` is
called, the run resumes from it.
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


class _ProgressLog:
    """dynesty print_func: appends a status line to ``path`` at most every ``every_s`` s."""

    def __init__(self, path: Path, every_s: float, label: str) -> None:
        self.path, self.every_s, self.label = Path(path), every_s, label
        self.t0 = self.last = time.time()
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def line(self, text: str) -> None:
        stamp = time.strftime("%Y-%m-%d %H:%M:%S")
        with self.path.open("a") as f:
            f.write(f"{stamp} [{self.label}] {text}\n")

    def __call__(self, itresult, niter, ncall, **kw) -> None:
        now = time.time()
        if now - self.last < self.every_s:
            return
        self.last = now
        from dynesty.utils import get_print_fn_args

        a = get_print_fn_args(itresult, niter, ncall, **kw)
        self.line(f"iter {niter} | " + " | ".join(a.long_str) +
                  f" | elapsed {(now - self.t0) / 60:.1f} min")  # fmt: skip


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
    progress_log: Path | None = None,
    progress_every_s: float = 60.0,
    checkpoint_file: Path | None = None,
    checkpoint_every_s: float = 300.0,
    label: str = "nested",
    slices: int | None = None,
) -> NestedResult:
    """Static nested sampling over ``space``'s prior (which must carry a Uranus box).

    With ``checkpoint_file`` set, the sampler state is saved every ``checkpoint_every_s`` and a
    run whose checkpoint exists resumes from it (same settings assumed; ``rng`` is then ignored,
    the checkpoint carries its own generator state). ``slices``: rslice slices per proposal;
    None keeps dynesty's default (3 + ndim). ``wall_s`` counts this call only."""
    from contextlib import nullcontext
    from multiprocessing import get_context

    import dynesty

    space._check_transformable()
    _init_worker(space, obs, context)
    kw = dict(nlive=nlive, bound=BOUND, sample=SAMPLE, periodic=space.periodic_unit_dims,
              rstate=rng)  # fmt: skip
    if slices is not None:
        kw["slices"] = int(slices)
    log = None if progress_log is None else _ProgressLog(progress_log, progress_every_s, label)
    ckpt = None if checkpoint_file is None else Path(checkpoint_file)
    if ckpt is not None:
        ckpt.parent.mkdir(parents=True, exist_ok=True)
    run_kw = dict(dlogz=dlogz, maxiter=maxiter, print_progress=log is not None, print_func=log,
                  checkpoint_file=None if ckpt is None else str(ckpt),
                  checkpoint_every=checkpoint_every_s)  # fmt: skip
    t0 = time.time()
    pool_cm = (get_context("spawn").Pool(n_workers, initializer=_init_worker,
                                         initargs=(space, obs, context))
               if n_workers > 1 else nullcontext(None))  # fmt: skip
    with pool_cm as pool:
        pool_kw = {} if pool is None else dict(pool=pool, queue_size=n_workers)
        if ckpt is not None and ckpt.exists():
            # dynesty's own checkpoint (pickle), written by this function for this run
            sampler = dynesty.NestedSampler.restore(str(ckpt), pool=pool)
            if log is not None:
                log.line(f"resumed from {ckpt} at iteration {sampler.it}")
            sampler.run_nested(resume=True, **run_kw)
        else:
            if log is not None:
                sl = "" if slices is None else f", {slices} slices"
                log.line(f"start: ndim {space.ndim}, nlive {nlive}, {BOUND}/{SAMPLE}{sl}, "
                         f"dlogz {dlogz}, {n_workers} workers")  # fmt: skip
            sampler = dynesty.NestedSampler(_log_likelihood, _prior_transform, space.ndim,
                                            **pool_kw, **kw)  # fmt: skip
            sampler.run_nested(**run_kw)
    res = NestedResult(space, sampler.results, nlive, dlogz, time.time() - t0)
    if log is not None:
        log.line(f"done: {res.niter} it, {res.ncall} calls, ln Z {res.logz:.3f} ± "
                 f"{res.logzerr:.3f}, {res.wall_s / 60:.1f} min")  # fmt: skip
    return res


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


# ------------------------------------------------------------------ partitioned (cells)

# Round-1 partition (a only); kept for the record of that validation.
A_SLICES_AU: tuple[tuple[float, float], ...] = ((20.0, 25.0), (25.0, 28.0), (28.0, 60.0))

# M5b partition: cells of (a [AU], e) with e = None for all of [0, e_max]. [28, 60] AU is
# split by eccentricity (circular vs eccentric), not by a: it holds mode A (e ≈ 0.03) and
# an eccentric branch (a ≳ 33, e ≈ 0.16–0.3) that one sub-run weighs inconsistently, and
# every alternative branch found so far is eccentric. The e edge assumes e_max = 0.3.
PARTITION_CELLS: tuple[tuple[tuple[float, float], tuple[float, float] | None], ...] = (
    ((20.0, 25.0), None),
    ((25.0, 28.0), None),
    ((28.0, 60.0), (0.0, 0.12)),
    ((28.0, 60.0), (0.12, 0.3)),
)


# Empirical ln Z error inflation (PHASE1_SPEC M5b, "Evidence errors"). On seed 0's validation
# the single-run ln Z scatter between rng seeds was ≈ 1.8× dynesty's quoted logzerr (pooled
# over the four cells, 5 degrees of freedom), so every ln Z error that enters a comparison is
# σ_eff = ERR_INFLATION · σ_quoted. Applies to Phase 2 evidence too.
ERR_INFLATION = 1.8


def cells_for(e_max: float) -> tuple:
    """``PARTITION_CELLS`` restricted to e ≤ e_max: e ranges are clipped to [0, e_max] and
    cells left empty are dropped (e ≤ 0.1 keeps three cells, without [28, 60] × e ≥ 0.12)."""
    out = []
    for a_rng, e_rng in PARTITION_CELLS:
        if e_rng is None:
            out.append((a_rng, None))
            continue
        lo, hi = e_rng[0], min(e_rng[1], e_max)
        if lo < hi:
            out.append((a_rng, (lo, hi)))
    return tuple(out)


def cell_space(space: PosteriorSpace, cell) -> PosteriorSpace:
    """The restricted prior of one partition cell ((a_lo, a_hi), e range or None)."""
    (a_lo, a_hi), e_rng = cell
    out = space.with_a_slice(a_lo, a_hi)
    return out if e_rng is None else out.with_e_slice(*e_rng)


@dataclass(frozen=True)
class PartitionedResult:
    """Nested runs on disjoint cells (a slice × e slice) of one prior, combined:
    Z = Σ_k w_k Z_k with w_k the cell's prior mass (``PosteriorSpace.slice_weight``). Errors
    are propagated from each run's ln Z error, assumed independent: σ(ln Z)² = Σ_k s_k² σ_k²,
    with s_k = w_k Z_k / Z the cell's posterior mass share;
    var(s_k) = Σ_j s_k² (δ_kj − s_j)² σ_j². Each σ_k is the run's quoted logzerr times
    ``err_scale`` (``ERR_INFLATION`` for reported results; 1 reproduces dynesty's own)."""

    slices: tuple[NestedResult, ...]
    err_scale: float = 1.0

    def __post_init__(self) -> None:
        w = sum(r.space.slice_weight for r in self.slices)
        if abs(w - 1.0) > 1e-9:
            raise ValueError(f"cells must partition the prior (weights sum to {w})")

    @property
    def ln_terms(self) -> np.ndarray:
        return np.array([math.log(r.space.slice_weight) + r.logz for r in self.slices])

    @property
    def logz(self) -> float:
        t = self.ln_terms
        return float(t.max() + math.log(np.exp(t - t.max()).sum()))

    @property
    def shares(self) -> np.ndarray:
        return np.exp(self.ln_terms - self.logz)

    @property
    def logzerr(self) -> float:
        sig = self.err_scale * np.array([r.logzerr for r in self.slices])
        return float(math.sqrt(np.sum((self.shares * sig) ** 2)))

    @property
    def share_errs(self) -> np.ndarray:
        s = self.shares
        sig = self.err_scale * np.array([r.logzerr for r in self.slices])
        eye = np.eye(len(s))
        return np.sqrt([np.sum(s[k] ** 2 * (eye[k] - s) ** 2 * sig**2) for k in range(len(s))])

    def samples_and_weights(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """All cells' samples with combined posterior weights (summing to 1) and the cell
        index of each sample."""
        q = np.concatenate([r.samples for r in self.slices])
        w = np.concatenate([s * r.weights for s, r in zip(self.shares, self.slices, strict=True)])
        k = np.concatenate([np.full(len(r.samples), i) for i, r in enumerate(self.slices)])
        return q, w / w.sum(), k

    def equal_weight(self, rng: np.random.Generator) -> np.ndarray:
        from dynesty.utils import resample_equal

        q, w, _ = self.samples_and_weights()
        return resample_equal(q, w, rstate=rng)
