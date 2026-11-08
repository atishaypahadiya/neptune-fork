"""Student-t measurement errors with the M6 group offsets marginalised (M6c sensitivity S3).

Each measurement i in group g: r_i − b_g ~ t_ν(scale s_g), with b_g ~ N(0, τ_g²) as in the
Gaussian model (``fork.model.observation``). There is no closed form for the marginal over b_g,
but each group's integral is one-dimensional:

  L_g(θ) = ∫ Π_{i∈g} t_ν((r_i − b)/s_g)/s_g · N(b; 0, τ_g²) db.

It is evaluated by Gauss–Hermite quadrature (``N_NODES`` nodes) centred on the mode b̂_g, found
by iteratively reweighted least squares (the EM fixed point for the t), with node spacing from
the expected information n_g(ν+1)/((ν+3)s_g²) + 1/τ_g². ``tests/test_student_t.py`` checks the
quadrature against adaptive integration, and the ν → ∞ limit against the Gaussian marginal.

Unlike the Gaussian path, −2 ln L here is **fully normalised** (no constant dropped), so the
ln Z of a Student-t fit is compared with a Gaussian fit's ln Z only after adding back the
Gaussian constant (``gaussian_norm_m2lnl``).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from fork.model.forward import ForwardContext
from fork.model.observation import GROUPS, RaDecObservations

N_NODES = 48
IRLS_MAX_ITER = 100
IRLS_TOL_ARCSEC = 1e-10
_X, _W = np.polynomial.hermite.hermgauss(N_NODES)
_LOG_W_PLUS_X2 = np.log(_W) + _X**2


def _log_t_const(nu: float) -> float:
    return math.lgamma((nu + 1) / 2) - math.lgamma(nu / 2) - 0.5 * math.log(nu * math.pi)


def t_offset_modes(resid: np.ndarray, group_index: np.ndarray, scale: np.ndarray,
                   tau: np.ndarray, nu: float) -> np.ndarray:  # fmt: skip
    """b̂_g maximising Σ_i ln t_ν((r_i − b)/s_g) − b²/(2τ_g²), per group (IRLS from the
    Gaussian solution)."""
    ng = len(scale)
    s2 = scale[group_index] ** 2
    n = np.bincount(group_index, minlength=ng).astype(float)
    b = np.bincount(group_index, weights=resid / s2, minlength=ng) / (n / scale**2 + 1 / tau**2)
    for _ in range(IRLS_MAX_ITER):
        z2 = (resid - b[group_index]) ** 2 / s2
        w = (nu + 1.0) / (nu + z2) / s2
        b_new = (np.bincount(group_index, weights=w * resid, minlength=ng)
                 / (np.bincount(group_index, weights=w, minlength=ng) + 1 / tau**2))  # fmt: skip
        done = np.max(np.abs(b_new - b)) < IRLS_TOL_ARCSEC
        b = b_new
        if done:
            break
    return b


def t_marginal_m2lnl(resid: np.ndarray, group_index: np.ndarray, scale: np.ndarray,
                     tau: np.ndarray, nu: float) -> float:  # fmt: skip
    """−2 ln L, normalised, with every group's offset integrated out (see module docstring)."""
    ng = len(scale)
    b_hat = t_offset_modes(resid, group_index, scale, tau, nu)
    n = np.bincount(group_index, minlength=ng).astype(float)
    h = 1.0 / np.sqrt(n * (nu + 1.0) / ((nu + 3.0) * scale**2) + 1.0 / tau**2)
    c = _log_t_const(nu)
    total = 0.0
    for g in range(ng):
        r = resid[group_index == g]
        b = b_hat[g] + math.sqrt(2.0) * h[g] * _X  # (K,)
        z2 = (r[None, :] - b[:, None]) ** 2 / scale[g] ** 2
        log_data = r.size * (c - math.log(scale[g])) - 0.5 * (nu + 1.0) * np.sum(
            np.log1p(z2 / nu), axis=1)  # fmt: skip
        log_prior = -0.5 * b**2 / tau[g] ** 2 - 0.5 * math.log(2 * math.pi * tau[g] ** 2)
        terms = _LOG_W_PLUS_X2 + log_data + log_prior
        m = float(np.max(terms))
        total += m + math.log(float(np.sum(np.exp(terms - m))))
        total += math.log(math.sqrt(2.0) * h[g])
    return -2.0 * total


def gaussian_norm_m2lnl(obs: RaDecObservations) -> float:
    """The constant the Gaussian path drops (``fork.model.observation``): the normalised
    −2 ln L is ``obs.chi2`` + this."""
    sig = obs.sigma_arcsec[obs.group_index]
    n = obs.group_counts()
    ratio = obs.offset_prior_arcsec**2 / obs.sigma_arcsec**2
    return float(np.sum(np.log(2 * math.pi * sig**2)) + np.sum(np.log1p(n * ratio)))


@dataclass(frozen=True)
class StudentTObservations(RaDecObservations):
    """``RaDecObservations`` with Student-t errors: ``chi2`` is the normalised marginal −2 ln L
    above; ``sigma_arcsec`` is kept (unused here) and ``t_scale_arcsec`` is the t scale."""

    nu: float
    t_scale_arcsec: np.ndarray  # (len(GROUPS),)

    @classmethod
    def from_gaussian(cls, obs: RaDecObservations, nu: float,
                      t_scale_arcsec: dict[str, float]) -> StudentTObservations:  # fmt: skip
        scale = np.array([t_scale_arcsec[g] for g in GROUPS], dtype=float)
        return cls(obs.geometry, obs.obs_id, obs.epoch_index, obs.is_dec, obs.value_rad,
                   obs.group_index, obs.sigma_arcsec, obs.offset_prior_arcsec, float(nu),
                   scale)  # fmt: skip

    def chi2(self, theta: np.ndarray, context: ForwardContext) -> float:
        resid = self.residuals_arcsec(*self.predict(theta, context))
        v = t_marginal_m2lnl(resid, self.group_index, self.t_scale_arcsec,
                             self.offset_prior_arcsec, self.nu)  # fmt: skip
        return v if math.isfinite(v) else math.inf

    def best_offsets_arcsec(self, resid_arcsec: np.ndarray) -> np.ndarray:
        return t_offset_modes(resid_arcsec, self.group_index, self.t_scale_arcsec,
                              self.offset_prior_arcsec, self.nu)  # fmt: skip
