"""The shared hidden-mass prior (``fork.inference.priors``): emcee's log_prior and dynesty's
prior_transform must describe the same distribution."""

import math

import numpy as np
import pytest

from fork.inference.fit import uranus_theta_to_x
from fork.inference.priors import (
    PRIOR_A_AU,
    PRIOR_LOG10_M_EARTH,
    PosteriorSpace,
    UranusBox,
    free_a_space,
    held_shape_space,
)

EPOCH = 2371628.5
BOX = UranusBox((19.14, -0.047, 0.0056, 1.6734), (6e-3, 4e-4, 4e-4, 1e-3))
PQ = (0.0066, 0.0007)


def _draws(space: PosteriorSpace, n: int, seed: int = 0) -> np.ndarray:
    u = np.random.default_rng(seed).random((n, space.ndim))
    return np.array([space.prior_transform(v) for v in u])


@pytest.mark.parametrize("e_max", [0.3, 0.1])
def test_prior_transform_is_the_log_prior_distribution(e_max):
    """Every transformed draw is inside the log_prior support (flat there), and the marginals are
    the spec's: log10 m and ln a uniform, e uniform on [0, e_max], ϖ and λ₀ uniform; Uranus fills
    its box uniformly."""
    space = free_a_space(PQ, e_max, BOX)
    q = _draws(space, 20000)
    assert all(space.log_prior(v) == 0.0 for v in q[:2000])
    log10_m, ln_a, k, h, lam0 = q[:, 4:].T
    e, varpi = k * k + h * h, np.arctan2(h, k) % (2 * math.pi)
    tol = 0.015  # ~2σ binomial at n = 20000 for these fractions
    assert np.mean(e < e_max / 2) == pytest.approx(0.5, abs=tol)  # uniform in e, not in √e
    assert np.mean(e < e_max / 4) == pytest.approx(0.25, abs=tol)
    assert e.max() <= e_max
    geo_mid = math.sqrt(PRIOR_A_AU[0] * PRIOR_A_AU[1])  # log-uniform: half the mass below √(lo·hi)
    assert np.mean(np.exp(ln_a) < geo_mid) == pytest.approx(0.5, abs=tol)
    assert np.exp(ln_a).min() >= PRIOR_A_AU[0] and np.exp(ln_a).max() <= PRIOR_A_AU[1]
    assert np.mean(log10_m < 1.0) == pytest.approx(0.5, abs=tol)
    assert log10_m.min() >= PRIOR_LOG10_M_EARTH[0] and log10_m.max() <= PRIOR_LOG10_M_EARTH[1]
    for ang in (varpi, lam0):
        assert np.mean(ang < math.pi / 2) == pytest.approx(0.25, abs=tol)
    rel = (q[:, :4] - np.array(BOX.centre)) / np.array(BOX.half_width)
    assert np.abs(rel).max() <= 1.0 and np.mean(rel[:, 0] < 0) == pytest.approx(0.5, abs=tol)


def test_log_prior_support_edges():
    space = free_a_space(PQ, 0.3, BOX)
    q = space.prior_transform(np.full(space.ndim, 0.5))
    assert space.log_prior(q) == 0.0
    for i, step in ((0, 1.01 * BOX.half_width[0]), (3, -1.01 * BOX.half_width[3])):
        bad = q.copy()
        bad[i] += step
        assert space.log_prior(bad) == -math.inf
    for i, v in ((4, 2.01), (4, -0.01), (5, math.log(19.9)), (5, math.log(60.1))):
        bad = q.copy()
        bad[i] = v
        assert space.log_prior(bad) == -math.inf
    bad = q.copy()
    bad[6], bad[7] = math.sqrt(0.301), 0.0
    assert space.log_prior(bad) == -math.inf
    assert BOX.log_volume == pytest.approx(float(np.sum(np.log(2 * np.array(BOX.half_width)))))


def test_unit_cube_layout_and_periodic_dims():
    free = free_a_space(PQ, 0.3, BOX)
    assert free.unit_names == ("u_a", "u_k", "u_h", "u_lam0", "log10_m", "ln_a", "e", "varpi",
                               "lam0")  # fmt: skip
    assert free.periodic_unit_dims == [7, 8]  # perturber ϖ, λ₀ (Uranus λ₀ is boxed)
    held = held_shape_space(PQ, 36.15, 0.10761, 5.0, BOX)
    assert held.unit_names == ("u_a", "u_k", "u_h", "u_lam0", "log10_m", "lam0")
    assert held.periodic_unit_dims == [5]
    q = _draws(held, 500)
    assert q.shape == (500, 6) and all(held.log_prior(v) == 0.0 for v in q)
    th = held.theta(q[0])
    assert th[7] == pytest.approx(36.15) and th[8] == pytest.approx(0.10761)


def test_prior_transform_requires_a_box():
    with pytest.raises(ValueError, match="Uranus box"):
        free_a_space(PQ).prior_transform(np.full(9, 0.5))


def test_uranus_box_contains_the_perturbed_uranus(offline_horizons):
    """Seed 0. The box is centred on the Uranus-only fit with 50 × (σ·√(χ²/dof)) half-widths.
    Uranus refitted with true Neptune planted sits well inside it, although it is 53 and 57 raw σ
    from the centre in a and e·cos ϖ, which an unscaled 50σ box would exclude. Measured:
    χ²/dof 21.6; the planted fit lies at ≤ 0.25 of the half-width in every coordinate (emcee
    samples of both seed-0 modes: ≤ 0.27)."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000
    from fork.data.synthetic import make_synthetic_observations
    from fork.inference.fit import fit_uranus_elements
    from fork.model import bodies
    from fork.model.forward import (
        ForwardContext,
        perturber_theta_from_state,
        uranus_inc_raan_at_epoch,
    )

    ctx = ForwardContext.from_horizons(EPOCH)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    assert box.chi2_per_dof > 10
    th_nep = perturber_theta_from_state(
        states_about_inner_barycenter_ecl_j2000("neptune", [EPOCH])[0], bodies.NEPTUNE.mass_msun
    )
    planted = fit_uranus_elements(obs, ctx, inc_raan, fixed_perturber_theta=th_nep)
    x = uranus_theta_to_x(planted.theta)[:4]
    rel = np.abs(x - np.array(box.centre)) / np.array(box.half_width)
    assert box.contains(x) and rel.max() < 0.35
    raw_sigma = np.array(box.sigma_scaled) / math.sqrt(box.chi2_per_dof)
    assert np.all(np.abs(x - np.array(box.centre))[:2] / raw_sigma[:2] > 50)


def test_a_slices_partition_the_prior_and_renormalise():
    """The three M5b slices carry prior masses ln(hi/lo)/ln 3 that sum to 1; inside a slice
    the transform stays on the slice and is log-uniform there; log_prior rejects a outside."""
    from fork.inference.nested import A_SLICES_AU

    full = free_a_space(PQ, 0.3, BOX)
    spaces = [full.with_a_slice(lo, hi) for lo, hi in A_SLICES_AU]
    w = [s.slice_weight for s in spaces]
    assert w == pytest.approx([math.log(hi / lo) / math.log(3.0) for lo, hi in A_SLICES_AU])
    assert sum(w) == pytest.approx(1.0) and full.slice_weight == 1.0
    sp = spaces[1]  # [25, 28]
    a = np.exp(_draws(sp, 5000)[:, 5])
    assert a.min() >= 25.0 and a.max() <= 28.0
    assert np.mean(a < math.sqrt(25.0 * 28.0)) == pytest.approx(0.5, abs=0.03)
    q = sp.prior_transform(np.full(sp.ndim, 0.5))
    q[5] = math.log(24.9)
    assert sp.log_prior(q) == -math.inf and full.log_prior(q) == 0.0
    with pytest.raises(ValueError):
        full.with_a_slice(15.0, 25.0)


def test_partitioned_combination_and_error_propagation():
    """Z = Σ w_k Z_k; shares s_k = w_k Z_k / Z; σ(ln Z)² = Σ s_k² σ_k²; delta-method share
    errors; combined weights sum to 1 and follow the shares."""
    from types import SimpleNamespace

    from fork.inference.nested import A_SLICES_AU, PartitionedResult

    full = free_a_space(PQ, 0.3, BOX)
    lnz, sig = [-75.0, -74.0, -73.2], [0.2, 0.25, 0.22]
    runs = []
    for i, ((lo, hi), z, e) in enumerate(zip(A_SLICES_AU, lnz, sig, strict=True)):
        runs.append(SimpleNamespace(space=full.with_a_slice(lo, hi), logz=z, logzerr=e,
                                    samples=np.full((4, 9), float(i)),
                                    weights=np.full(4, 0.25)))  # fmt: skip
    p = PartitionedResult(tuple(runs))
    w = np.array([r.space.slice_weight for r in runs])
    z = np.sum(w * np.exp(lnz))
    assert p.logz == pytest.approx(math.log(z))
    s = w * np.exp(lnz) / z
    assert p.shares == pytest.approx(s) and p.shares.sum() == pytest.approx(1.0)
    assert p.logzerr == pytest.approx(math.sqrt(np.sum((s * np.array(sig)) ** 2)))
    # share error, slice 0 by hand: s0² [(1−s0)² σ0² + s1² σ1² + s2² σ2²]
    e0 = s[0] * math.sqrt((1 - s[0]) ** 2 * sig[0] ** 2 + s[1] ** 2 * sig[1] ** 2
                          + s[2] ** 2 * sig[2] ** 2)  # fmt: skip
    assert p.share_errs[0] == pytest.approx(e0)
    _, wq, k = p.samples_and_weights()
    assert wq.sum() == pytest.approx(1.0)
    assert [wq[k == i].sum() for i in range(3)] == pytest.approx(s)
    with pytest.raises(ValueError, match="partition"):
        PartitionedResult(tuple(runs[:2]))


def test_e_cells_partition_the_prior_and_renormalise():
    """M5b cells: [28, 60] AU split at e = 0.12. Cell weights are w_a · (e_hi − e_lo)/e_max and
    sum to 1 over ``PARTITION_CELLS``; inside an e cell the transform stays in it and is
    uniform in e; log_prior rejects e outside; an e slice needs (k, h) free and e_hi ≤ e_max."""
    from fork.inference.nested import PARTITION_CELLS, cell_space

    full = free_a_space(PQ, 0.3, BOX)
    spaces = [cell_space(full, c) for c in PARTITION_CELLS]
    w_a = math.log(60.0 / 28.0) / math.log(3.0)
    assert [s.slice_weight for s in spaces[2:]] == pytest.approx([w_a * 0.4, w_a * 0.6])
    assert sum(s.slice_weight for s in spaces) == pytest.approx(1.0)
    for sp, (e_lo, e_hi) in ((spaces[2], (0.0, 0.12)), (spaces[3], (0.12, 0.3))):
        q = _draws(sp, 5000)
        e, a = q[:, 6] ** 2 + q[:, 7] ** 2, np.exp(q[:, 5])
        assert e.min() >= e_lo and e.max() <= e_hi and a.min() >= 28.0
        assert np.mean(e < 0.5 * (e_lo + e_hi)) == pytest.approx(0.5, abs=0.03)  # uniform in e
        assert all(sp.log_prior(v) == 0.0 for v in q[:500])
    q = spaces[3].prior_transform(np.full(spaces[3].ndim, 0.5))
    q[6], q[7] = math.sqrt(0.1), 0.0  # e = 0.1: in the circular cell, not the eccentric one
    assert spaces[3].log_prior(q) == -math.inf and spaces[2].log_prior(q) == 0.0
    with pytest.raises(ValueError):
        full.with_e_slice(0.12, 0.31)
    with pytest.raises(ValueError):
        held_shape_space(PQ, 36.15, 0.10761, 5.0, BOX).with_e_slice(0.0, 0.1)


def test_partitioned_result_accepts_the_four_cells():
    from types import SimpleNamespace

    from fork.inference.nested import PARTITION_CELLS, PartitionedResult, cell_space

    full = free_a_space(PQ, 0.3, BOX)
    runs = tuple(SimpleNamespace(space=cell_space(full, c), logz=-74.0, logzerr=0.2,
                                 samples=np.zeros((2, 9)), weights=np.full(2, 0.5))
                 for c in PARTITION_CELLS)  # fmt: skip
    p = PartitionedResult(runs)
    assert p.logz == pytest.approx(-74.0)  # equal Z_k: Z = Σ w_k Z_k = Z_k
    assert p.shares == pytest.approx([r.space.slice_weight for r in runs])
    with pytest.raises(ValueError, match="partition"):
        PartitionedResult(runs[:3])
