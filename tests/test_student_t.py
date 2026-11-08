"""Student-t marginal likelihood (M6c sensitivity S3): quadrature against adaptive integration,
the Gaussian limit, and the offset modes."""

import math

import numpy as np
import pytest
from scipy import integrate, stats

from fork.model import student_t as st


def _group_data(rng):
    r = np.concatenate([rng.standard_t(4, 60) * 1.5 + 0.8, [9.0, -7.5]])  # two outliers
    return r, np.zeros(r.size, dtype=int)


def test_quadrature_matches_adaptive_integration():
    rng = np.random.default_rng(3)
    for nu, s, tau in [(4.0, 1.5, 2.0), (4.0, 4.0, 5.0), (2.5, 0.8, 0.5)]:
        r, g = _group_data(rng)

        def f(b, r=r, nu=nu, s=s, tau=tau):
            return math.exp(np.sum(stats.t.logpdf(r - b, nu, scale=s))
                            + stats.norm.logpdf(b, 0, tau) + 300.0)  # fmt: skip

        b_hat = st.t_offset_modes(r, g, np.array([s]), np.array([tau]), nu)[0]
        val, _ = integrate.quad(f, b_hat - 30 * s, b_hat + 30 * s, points=[b_hat], limit=400,
                                  epsabs=0, epsrel=1e-12)  # fmt: skip
        ref = -2.0 * (math.log(val) - 300.0)
        got = st.t_marginal_m2lnl(r, g, np.array([s]), np.array([tau]), nu)
        assert got == pytest.approx(ref, abs=1e-7)


def test_large_nu_reproduces_the_gaussian_marginal():
    """ν → ∞: the t marginal equals the Gaussian path's −2 ln L (whitened² + the dropped
    constant). ν = 1e7, because a larger ν loses digits in the lgamma difference."""
    from fork.model.observation import GROUPS

    rng = np.random.default_rng(5)
    ng = len(GROUPS)
    g = np.repeat(np.arange(ng), 25)
    sig = rng.uniform(1.0, 5.0, ng)
    tau = np.where(np.arange(ng) < 2, 5.0, 2.0)
    r = rng.normal(0.0, sig[g]) + rng.normal(0.0, tau)[g]
    b = (np.bincount(g, weights=r / sig[g] ** 2) / (25 / sig**2 + 1 / tau**2))
    chi2 = np.sum(((r - b[g]) / sig[g]) ** 2) + np.sum((b / tau) ** 2)
    const = np.sum(np.log(2 * math.pi * sig[g] ** 2)) + np.sum(np.log1p(25 * tau**2 / sig**2))
    got = st.t_marginal_m2lnl(r, g, sig, tau, 1e7)
    assert got == pytest.approx(chi2 + const, abs=1e-4)


def test_offset_mode_is_robust_to_an_outlier():
    r = np.array([1.0] * 30 + [40.0])
    g = np.zeros(r.size, dtype=int)
    b_t = st.t_offset_modes(r, g, np.array([1.0]), np.array([1e6]), 4.0)[0]
    assert b_t == pytest.approx(1.0, abs=0.05)  # the Gaussian mean would be ≈ 2.26


def test_gaussian_norm_matches_the_observation_model():
    """gaussian_norm_m2lnl + obs.chi2 is the normalised Gaussian −2 ln L: check it against the
    ν → ∞ t marginal on the real layout."""
    from fork.data.calibration import real_observations
    from fork.model import observation as om

    obs = real_observations(om.read_rows(), om.ErrorModel.load())
    rng = np.random.default_rng(2)
    r = rng.normal(0.0, obs.sigma_arcsec[obs.group_index])
    chi2 = float(np.sum(obs.whitened_from_residuals(r) ** 2))
    got = st.t_marginal_m2lnl(r, obs.group_index, obs.sigma_arcsec, obs.offset_prior_arcsec, 1e7)
    assert got == pytest.approx(chi2 + st.gaussian_norm_m2lnl(obs), abs=1e-3)
