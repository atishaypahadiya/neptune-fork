"""Profile starts: the a = 22 AU backstop and posterior-seeded (full-θ) starts."""

import math

import numpy as np
import pytest

from fork.inference import profile as pr
from fork.inference.fit import _unpack_perturber_z, le_verrier_space, uranus_theta_to_x


def test_start_from_theta_round_trips_the_full_theta():
    theta = np.array([19.2, 0.047, 0.0135, 1.29, 2.97, 4.1,  # Uranus a, e, i, Ω, ϖ, λ₀
                      1.25, 31.2, 0.12, 2.2, 5.9])  # perturber log10 m, a, e, ϖ, λ₀  # fmt: skip
    space = le_verrier_space(None)
    pq = uranus_theta_to_x(theta[:6])[4:]
    lon = 5.7
    z, lam0 = pr.start_from_theta(space, theta, lon, lon)  # target = its own longitude
    assert lam0 == pytest.approx(theta[10], abs=1e-12)
    back = _unpack_perturber_z(np.concatenate([z, [lam0]]), pq, space.interior())
    np.testing.assert_allclose(back, theta, rtol=0, atol=1e-12)
    _, shifted = pr.start_from_theta(space, theta, lon, lon + 0.1)
    assert shifted == pytest.approx(theta[10] + 0.1, abs=1e-12)


def test_backstop_is_free_a_only():
    with pytest.raises(ValueError, match="free-a"):
        pr._start_vectors(le_verrier_space(38.4), np.zeros(4), None, None, 0.0, 8, backstop=True)


@pytest.mark.slow
def test_backstop_reaches_the_a22_branch_seed4_345():
    """Seed 4, free a, 1846 longitude held at 345°: the standard starts stop at χ² 93.48
    (a ≈ 26 AU on the e bound); the a = 22 AU backstop reaches the lower branch at 91.69
    (a ≈ 22 AU, e = 0.3), found by the pre-bounded-fit export (2026-11-04 diagnostics)."""
    from fork.data.synthetic import make_synthetic_observations
    from fork.export.to_json import A1_EPOCH_JD_TDB as epoch
    from fork.inference.fit import fit_uranus_elements
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    ctx = ForwardContext.from_horizons(epoch)
    obs, _ = make_synthetic_observations(np.random.default_rng(4))
    geo = pr.GalleGeometry.from_horizons()
    space = le_verrier_space(None)
    x = uranus_theta_to_x(fit_uranus_elements(obs, ctx, uranus_inc_raan_at_epoch(epoch)).theta)
    xu, pq = x[:4].copy(), x[4:].copy()
    t = math.radians(345.0)

    def best(starts):
        fits = [pr._run_constrained(z, lam0, pq, space, obs, ctx, geo, t) for z, lam0 in starts]
        return min((f for f in fits if f.success), key=lambda f: f.chi2)

    std = pr._start_vectors(space, xu, ctx, geo, t, 8)
    extra = pr._start_vectors(space, xu, ctx, geo, t, 8, backstop=True)[len(std) :]
    assert len(extra) == 2
    b = best(extra)
    assert b.chi2 == pytest.approx(91.69, abs=0.05)
    assert b.theta[7] == pytest.approx(22.0, abs=0.5) and b.theta[8] == pytest.approx(0.3, abs=1e-6)
    assert best(std).chi2 > b.chi2 + 1.0  # the gap the backstop closes


def test_seeded_starts_best_chi2_in_window_else_nearest():
    from fork.export.profile_seeds import seeded_starts

    lon = np.radians([320.0, 321.5, 322.5, 340.0])
    theta = np.arange(4.0)[:, None] * np.ones((4, 11))
    cell = {"theta": theta, "chi2": np.array([5.0, 3.0, 4.0, 1.0]), "lon_rad": lon}
    starts = seeded_starts([cell, cell], np.radians([321.0, 330.0, 359.0]))
    assert len(starts) == 3 and all(len(row) == 2 for row in starts)
    assert starts[0][0][0][0] == 1.0  # 320, 321.5, 322.5 in ±2°: lowest χ² is 321.5
    assert starts[1][0][0][0] == 2.0  # none within ±2° of 330: nearest is 322.5 (7.5° vs 10°)
    assert starts[2][0][0][0] == 3.0  # 359: nearest is 340
    assert starts[0][0][1] == pytest.approx(math.radians(321.5))
