"""Fast-then-reweight (PHASE1_SPEC, M8 amendment, "Fast then reweight"; 5fc089f).

Sample with a fast likelihood, then evaluate the full one on N equal-weight posterior samples
per cell and importance-reweight:

- w_i ∝ exp(Δ_i), Δ_i = ln L_full(φ_i) − ln L_fast(φ_i);
- ln Z_k = ln Z_k,fast + ln[(1/N_k) Σ_i exp(Δ_i)], whose Monte Carlo error is
  √((N_k/ESS_k − 1)/N_k) (delta method on the mean weight);
- cells combined as in ``PartitionedResult``: Z = Σ_k c_k Z_k with c_k the cell's prior mass;
  only the sampler's part of each cell's ln Z error carries ERR_INFLATION;
- ESS = (Σ w)² / Σ w² over the combined weighted samples, and within any subset (the ridge
  branches), with the subset's own N reported beside it.

Precedents: relative binning (Zackay, Dai & Venumadhav 2018, arXiv:1806.08792); likelihood
reweighting in gravitational-wave inference (Payne, Talbot & Thrane 2019, PRD 100, 123017).
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

N_PER_CELL = 6000  # M6c's derived-sample count (m6c.DERIVED_MAX_SAMPLES)
ESS_GATE_OVERALL = 0.5
ESS_GATE_BRANCH = 0.3


def ess(log_w: np.ndarray) -> float:
    """Kish effective sample size of (unnormalised, log) weights; −inf weights count as 0."""
    lw = np.asarray(log_w, dtype=float)
    if lw.size == 0 or not np.any(np.isfinite(lw)):
        return 0.0
    return float(math.exp(2.0 * logsumexp(lw) - logsumexp(2.0 * lw)))


@dataclass(frozen=True)
class CellReweight:
    delta: np.ndarray  # (N,) ln L_full − ln L_fast at the cell's equal-weight samples
    ln_z_fast: float
    ln_z_fast_err: float  # the sampler's error (not yet inflated)
    ln_correction: float  # ln[(1/N) Σ exp(Δ)]
    ln_correction_err: float  # Monte Carlo
    ess: float

    @property
    def n(self) -> int:
        return int(self.delta.size)

    @property
    def ln_z(self) -> float:
        return self.ln_z_fast + self.ln_correction


def reweight_cell(delta: np.ndarray, ln_z_fast: float, ln_z_fast_err: float) -> CellReweight:
    d = np.asarray(delta, dtype=float)
    n = d.size
    corr = float(logsumexp(d) - math.log(n))
    e = ess(d)
    corr_err = math.sqrt(max(n / e - 1.0, 0.0) / n) if e > 0 else math.inf
    return CellReweight(d, float(ln_z_fast), float(ln_z_fast_err), corr, corr_err, e)


@dataclass(frozen=True)
class Reweighted:
    """The reweighted fit: cells combined, samples with their final (log) weights."""

    cells: tuple[CellReweight, ...]
    log_prior_mass: np.ndarray  # (n_cells,) ln c_k
    samples: np.ndarray  # (Σ N_k, d) the cells' equal-weight samples, stacked
    cell_of: np.ndarray  # (Σ N_k,) which cell each sample came from
    log_w: np.ndarray  # (Σ N_k,) final log weights (unnormalised)
    ln_z: float
    ln_z_err: float  # √(Σ_k s_k² [(ERR_INFLATION σ_fast,k)² + σ_MC,k²]), s_k = c_k Z_k / Z

    @property
    def ess(self) -> float:
        return ess(self.log_w)

    @property
    def ess_frac(self) -> float:
        return self.ess / self.log_w.size

    def ess_within(self, mask: np.ndarray) -> tuple[float, int]:
        """(ESS / N, N) of the weights restricted to the samples in ``mask``."""
        mask = np.asarray(mask, dtype=bool)
        n = int(mask.sum())
        return (ess(self.log_w[mask]) / n if n else math.nan), n

    def normalised_weights(self) -> np.ndarray:
        return np.exp(self.log_w - logsumexp(self.log_w))

    def equal_weight(self, rng: np.random.Generator, n: int | None = None) -> np.ndarray:
        """Resample to equal weights (systematic resampling)."""
        w = self.normalised_weights()
        n = w.size if n is None else n
        pos = (rng.random() + np.arange(n)) / n
        return self.samples[np.minimum(np.searchsorted(np.cumsum(w), pos), w.size - 1)]


def combine(cells: Sequence[CellReweight], samples: Sequence[np.ndarray],
            log_prior_mass: Sequence[float], err_inflation: float) -> Reweighted:  # fmt: skip
    """Combine reweighted cells: ln Z = logsumexp(ln c_k + ln Z_k); a sample of cell k weighs
    c_k Z_k,fast · exp(Δ_i) / N_k (its cell's share of the fast posterior, times its weight)."""
    lc = np.asarray(log_prior_mass, dtype=float)
    lz = np.array([c.ln_z for c in cells])
    ln_z = float(logsumexp(lc + lz))
    share = np.exp(lc + lz - ln_z)
    var = sum(s**2 * ((err_inflation * c.ln_z_fast_err) ** 2 + c.ln_correction_err**2)
              for s, c in zip(share, cells, strict=True))  # fmt: skip
    log_w = np.concatenate([lc[k] + c.ln_z_fast + c.delta - math.log(c.n)
                            for k, c in enumerate(cells)])  # fmt: skip
    cell_of = np.concatenate([np.full(c.n, k) for k, c in enumerate(cells)])
    return Reweighted(tuple(cells), lc, np.concatenate(list(samples)), cell_of, log_w, ln_z,
                      math.sqrt(var))  # fmt: skip


def reweight_fit(cell_samples: Sequence[np.ndarray], ln_z_fast: Sequence[float],
                 ln_z_fast_err: Sequence[float], log_prior_mass: Sequence[float],
                 log_l_fast: Callable[[np.ndarray], float],
                 log_l_full: Callable[[np.ndarray], float], err_inflation: float,
                 map_fn: Callable = map) -> Reweighted:  # fmt: skip
    """The whole step: both likelihoods at every cell's equal-weight samples (``map_fn`` may be
    a pool's map), then ``reweight_cell`` and ``combine``."""
    cells = []
    for k, s in enumerate(cell_samples):
        full = np.fromiter(map_fn(log_l_full, list(s)), dtype=float, count=len(s))
        fast = np.fromiter(map(log_l_fast, list(s)), dtype=float, count=len(s))
        cells.append(reweight_cell(full - fast, ln_z_fast[k], ln_z_fast_err[k]))
    return combine(cells, cell_samples, log_prior_mass, err_inflation)
