"""Fast-then-reweight on a toy pair of Gaussian likelihoods (PHASE1_SPEC, M8 amendment).

A uniform prior on [−V/2, V/2] (V ≫ the widths), L_fast = N(x; μ₁, s₁), L_full = c·N(x; μ₂, s₂).
Then Z_full / Z_fast = c, so the ln Z correction is ln c, and with w = L_full / L_fast under the
fast posterior N(μ₁, s₁):
    ESS / N → E[w]² / E[w²] = 1 / ∫ N₂² / N₁
            = s₂ √(2s₁² − s₂²) / s₁ · exp(−(μ₂ − μ₁)² / (2s₁² − s₂²)).
"""

import math

import numpy as np
import pytest
from scipy import integrate, stats

from fork.inference import reweight as rw

MU1, S1, MU2, S2, LN_C = 0.0, 1.0, 0.3, 0.8, 1.7


def _expected_ess_frac() -> float:
    d = 2 * S1**2 - S2**2
    return S2 * math.sqrt(d) / S1 * math.exp(-((MU2 - MU1) ** 2) / d)


def _delta(x):
    return LN_C + stats.norm.logpdf(x, MU2, S2) - stats.norm.logpdf(x, MU1, S1)


def test_closed_form_ess_matches_quadrature():
    val, _ = integrate.quad(lambda x: stats.norm.pdf(x, MU2, S2) ** 2 / stats.norm.pdf(x, MU1, S1),
                            -30, 30)  # fmt: skip
    assert 1.0 / val == pytest.approx(_expected_ess_frac(), rel=1e-10)


def test_toy_gaussian_ess_and_ln_z_correction():
    rng = np.random.default_rng(42)
    x = rng.normal(MU1, S1, 200_000)
    cell = rw.reweight_cell(_delta(x), ln_z_fast=-3.0, ln_z_fast_err=0.1)
    assert cell.ess / cell.n == pytest.approx(_expected_ess_frac(), rel=0.02)
    assert abs(cell.ln_correction - LN_C) < 4 * cell.ln_correction_err
    assert cell.ln_z == pytest.approx(-3.0 + cell.ln_correction, abs=1e-12)


def test_monte_carlo_error_matches_the_spread_over_repeats():
    rng = np.random.default_rng(7)
    corr, errs = [], []
    for _ in range(300):
        c = rw.reweight_cell(_delta(rng.normal(MU1, S1, 2000)), 0.0, 0.0)
        corr.append(c.ln_correction)
        errs.append(c.ln_correction_err)
    assert np.std(corr) == pytest.approx(np.mean(errs), rel=0.15)
    assert abs(np.mean(corr) - LN_C) < 3 * np.std(corr) / math.sqrt(len(corr)) + 1e-3


def test_two_cells_combine_like_a_partitioned_result():
    """Two cells (prior masses 0.3, 0.7); in cell 2 the full and fast likelihoods agree
    (Δ = 0); in cell 1 Δ is the toy's. The combined ln Z and the cells' weight shares follow."""
    rng = np.random.default_rng(3)
    x1 = rng.normal(MU1, S1, 100_000)
    x2 = rng.normal(5.0, 1.0, 50_000)
    c1 = rw.reweight_cell(_delta(x1), -2.0, 0.05)
    c2 = rw.reweight_cell(np.zeros(x2.size), -2.5, 0.05)
    assert c2.ln_correction == pytest.approx(0.0, abs=1e-12) and c2.ess == pytest.approx(x2.size)
    res = rw.combine([c1, c2], [x1[:, None], x2[:, None]], [math.log(0.3), math.log(0.7)],
                     err_inflation=1.8)  # fmt: skip
    want = math.log(0.3 * math.exp(-2.0 + c1.ln_correction) + 0.7 * math.exp(-2.5))
    assert res.ln_z == pytest.approx(want, abs=1e-12)
    w = res.normalised_weights()
    share1 = 0.3 * math.exp(-2.0 + c1.ln_correction) / math.exp(want)
    assert w[res.cell_of == 0].sum() == pytest.approx(share1, rel=1e-10)
    frac, n = res.ess_within(res.cell_of == 1)
    assert n == x2.size and frac == pytest.approx(1.0)
    s1, s2 = share1, 1 - share1
    expect_err = math.sqrt(s1**2 * ((1.8 * 0.05) ** 2 + c1.ln_correction_err**2)
                           + s2**2 * (1.8 * 0.05) ** 2)  # fmt: skip
    assert res.ln_z_err == pytest.approx(expect_err, rel=1e-12)
    eq = res.equal_weight(np.random.default_rng(0), 20_000)
    assert np.mean(eq[:, 0] > 2.5) == pytest.approx(s2, abs=0.02)


def test_reweight_fit_runs_both_likelihoods():
    rng = np.random.default_rng(5)
    s = rng.normal(MU1, S1, (4000, 1))
    res = rw.reweight_fit([s], [0.0], [0.0], [0.0],
                          lambda q: float(stats.norm.logpdf(q[0], MU1, S1)),
                          lambda q: float(LN_C + stats.norm.logpdf(q[0], MU2, S2)),
                          err_inflation=1.8)  # fmt: skip
    assert res.ess_frac == pytest.approx(_expected_ess_frac(), abs=0.05)
    assert res.ln_z == pytest.approx(LN_C, abs=4 * res.ln_z_err)
