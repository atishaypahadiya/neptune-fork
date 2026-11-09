"""M8 linearized likelihood (PHASE1_SPEC M8; plan in docs/notes/M8.md).

Toy problems check the marginal against brute-force numerical integration and the truncation
bound against brute-force box mass. The M6 data structure checks the finite-difference
Jacobian against independent difference schemes, and a synthetic seed runs A8a with small N.
"""

import math

import numpy as np
import pytest
from scipy import integrate, stats

from fork.model import linearized as lz

# ------------------------------------------------------------------ toy problems (no integrator)


def _toy(eps: float, rng_seed: int = 0):
    """z(u) = A u − y + eps · (u₀², u₁², u₀u₁, …) on 2-D u: linear for eps = 0."""
    rng = np.random.default_rng(rng_seed)
    a = rng.normal(0.0, 1.0, (40, 2)) * np.array([3.0, 0.5])
    y = a @ np.array([0.7, -1.3]) + rng.normal(0.0, 1.0, 40)
    c = rng.normal(0.0, 1.0, (40, 3))

    def z(u):
        u = np.asarray(u, dtype=float)
        return a @ u - y + eps * (c @ np.array([u[0] ** 2, u[1] ** 2, u[0] * u[1]]))

    return z


def _brute_m2lnl(z, lo, hi):
    """−2 ln (∫_box exp(−½‖z‖²) du / V_box), by adaptive quadrature, shifted for range."""
    grid = [(x, y) for x in np.linspace(lo[0], hi[0], 81) for y in np.linspace(lo[1], hi[1], 81)]
    shift = min(float(np.sum(z(np.array(p)) ** 2)) for p in grid)
    f = lambda y, x: math.exp(-0.5 * (float(np.sum(z(np.array([x, y])) ** 2)) - shift))  # noqa: E731
    val, err = integrate.dblquad(f, lo[0], hi[0], lo[1], hi[1], epsabs=1e-13, epsrel=1e-10)
    assert err < 1e-6 * val
    vol = float(np.prod(np.asarray(hi) - np.asarray(lo)))
    return shift - 2.0 * math.log(val / vol)


def test_linear_marginal_equals_brute_force_integral_on_a_linear_toy():
    z = _toy(0.0)
    u0 = np.array([5.0, 5.0])  # far from the optimum: exact for a linear z anyway
    jac = lz.fd_jacobian(z, u0, np.array([1e-3, 1e-3]), "central")
    lo, hi = np.array([-2.0, -4.0]), np.array([3.0, 2.0])
    m = lz.linear_marginal(z(u0), jac, u0, float(np.sum(np.log(hi - lo))))
    assert m.m2lnl == pytest.approx(_brute_m2lnl(z, lo, hi), abs=1e-6)
    assert m.chi2_min == pytest.approx(float(np.sum(z(m.u_star) ** 2)), abs=1e-9)


def test_refinement_and_laplace_approach_the_brute_force_integral_on_a_nonlinear_toy():
    z = _toy(0.02)
    lo, hi = np.array([-3.0, -6.0]), np.array([4.0, 4.0])
    logv = float(np.sum(np.log(hi - lo)))
    truth = _brute_m2lnl(z, lo, hi)
    u0, h = np.array([2.5, 2.0]), np.array([1e-4, 1e-4])
    jac0 = lz.fd_jacobian(z, u0, h, "central")
    a = lz.gauss_newton_marginal(z, u0, jac0, logv, refine=False, step=h)
    b = lz.gauss_newton_marginal(z, u0, jac0, logv, refine=True, step=h)
    u_min, chi2_min, _ = lz.profile_min(z, b.u_star, np.array([1.0, 1.0]))
    lap = lz.laplace_m2lnl(z, u_min, np.array([0.05, 0.05]), logv)
    assert a.n_integrations == 1 and b.n_integrations == 2 + 2
    assert abs(b.m2lnl - truth) < abs(a.m2lnl - truth)
    assert abs(b.chi2_min - chi2_min) < abs(a.chi2_min - chi2_min)
    assert abs(lap - truth) < 0.01  # Laplace error is second order in the curvature


def test_truncation_bound_bounds_the_brute_force_outside_mass():
    cov = np.array([[1.0, 0.6], [0.6, 2.0]])
    mu = np.array([0.0, 0.0])
    lo, hi = np.array([-1.5, -6.0]), np.array([6.0, 2.5])  # near edges on two faces
    pdf = stats.multivariate_normal(mu, cov).pdf
    inside, _ = integrate.dblquad(lambda y, x: pdf([x, y]), lo[0], hi[0], lo[1], hi[1],
                                  epsabs=1e-12)  # fmt: skip
    p_out, m2_err = lz.truncation_bound(mu, cov, lo, hi)
    assert 1.0 - inside <= p_out <= 1.5 * (1.0 - inside)
    assert m2_err == pytest.approx(-2.0 * math.log1p(-p_out), rel=1e-12)
    # far from every edge (28σ and 40σ) the probability underflows, but its log is kept:
    # two faces at 40/1 σ and two at 40/√2 σ
    expect = np.logaddexp(
        math.log(2.0) + stats.norm.logsf(40.0),
        math.log(2.0) + stats.norm.logsf(40.0 / math.sqrt(2.0)),
    ) / math.log(10)
    got = lz.truncation_log10_bound(mu, cov, mu - 40.0, mu + 40.0)
    assert got == pytest.approx(expect, abs=1e-9) and got < -170
    assert lz.truncation_bound(mu, cov, mu - 40.0, mu + 40.0)[1] < 1e-170


# ------------------------------------------------------------------ M6 data (integrator)


@pytest.fixture(scope="module")
def synthetic():
    """A synthetic seed on Table 77's layout: real Uranus (with the real Neptune) + offsets +
    noise from the headline error model; the box by the M6c rule on these values."""
    from fork.export.to_json import A1_EPOCH_JD_TDB
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.inference.priors import free_a_space
    from fork.model import observation as om
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    structure = om.structure_from_rows(om.read_rows(), om.ErrorModel.load())
    vals, _ = om.synthetic_values(structure, om.true_uranus_ra_dec(structure),
                                  np.random.default_rng(808))  # fmt: skip
    obs = structure.with_values(vals)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    box = uranus_box_radec(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    space = free_a_space(pq, 0.3, box, inclined=True)
    return obs, ctx, space, lz.LinearizedProblem.build(obs, ctx, space)


def _neptune_phi(space, ctx):
    """Neptune's own φ at the epoch (a test target only; the fit path never reads it)."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.model.bodies import NEPTUNE
    from fork.model.forward import inclined_perturber_theta_from_state

    box = space.uranus_box
    theta_u = space.theta(np.concatenate([box.centre, np.zeros(7)]))[:6]
    st = about("neptune", [ctx.epoch_jd_tdb])[0]
    theta_p = inclined_perturber_theta_from_state(st, NEPTUNE.mass_msun, (theta_u[2], theta_u[3]))
    return space.q_from_theta(np.concatenate([theta_u, theta_p]))[4:]


def test_jacobian_agrees_with_independent_difference_schemes(synthetic):
    _, ctx, space, prob = synthetic
    phi = _neptune_phi(space, ctx)
    z = prob.z_fn(phi)
    u = prob.u_ref + 3.0 * prob.sigma_raw
    h = lz.FD_REL_STEP * prob.sigma_raw
    j_c = lz.fd_jacobian(z, u, h, "central")
    j_2 = lz.fd_jacobian(z, u, 2.0 * h, "central")
    j_rich = (4.0 * j_c - j_2) / 3.0
    j_f = lz.fd_jacobian(z, u, h, "forward")
    norm = np.linalg.norm(j_rich, axis=0)
    assert np.max(np.linalg.norm(j_c - j_rich, axis=0) / norm) < 1e-6
    assert np.max(np.linalg.norm(j_f - j_rich, axis=0) / norm) < 1e-3
    # the linear prediction along random 1-raw-σ directions
    rng = np.random.default_rng(5)
    z0 = z(u)
    for _ in range(3):
        d = rng.normal(0.0, 1.0, 4) * prob.sigma_raw
        lin = j_c @ d
        assert np.linalg.norm(z(u + d) - z0 - lin) < 1e-3 * np.linalg.norm(lin)


def test_reference_jacobian_is_the_uranus_only_one_at_the_box_centre(synthetic):
    obs, ctx, space, prob = synthetic
    np.testing.assert_allclose(prob.u_ref, space.uranus_box.centre, rtol=0, atol=0)
    from fork.inference.fit import uranus_x_to_theta

    def z_u(u):
        return obs.whitened(uranus_x_to_theta(np.concatenate([u, space.pq])), ctx)

    j = lz.fd_jacobian(z_u, prob.u_ref, lz.FD_REL_STEP * prob.sigma_raw, "central")
    np.testing.assert_allclose(prob.jac_ref, j, rtol=1e-12, atol=0)


def test_student_t_observations_are_refused(synthetic):
    from fork.model.observation import GROUPS
    from fork.model.student_t import StudentTObservations

    obs, ctx, space, _ = synthetic
    t_obs = StudentTObservations.from_gaussian(obs, 4.0, dict.fromkeys(GROUPS, 1.0))
    with pytest.raises(TypeError, match="Student-t"):
        lz.LinearizedProblem.build(t_obs, ctx, space)


def test_a8a_small_n_on_a_synthetic_seed(synthetic):
    """A8a with small N: near the truth, mode (b) meets the spec's tolerances against the
    full likelihood minimised over u and its exact-Hessian Laplace marginal, and the box
    truncation is measured negligible."""
    _, ctx, space, prob = synthetic
    phi0 = _neptune_phi(space, ctx)
    rng = np.random.default_rng(11)
    jitter = np.array([0.02, 0.01, 0.01, 0.01, 0.002, 0.001, 0.001])
    for k in range(3):
        phi = phi0 + (0.0 if k == 0 else 1.0) * rng.normal(0.0, 1.0, 7) * jitter
        b = prob.marginal(phi, "b")
        u_min, chi2_full, _ = prob.profile(phi, b.u_star)
        lap = prob.laplace(phi, u_min)
        assert abs(b.chi2_min - chi2_full) < 0.1
        assert abs(b.m2lnl - lap) < 0.2
        assert prob.truncation(b)[1] < 1e-6


# ------------------------------------------------------------------ real posteriors (slow)


@pytest.mark.slow
def test_a8a_mode_b_on_real_r2_posterior_points():
    """A8a's slow part on the real data: R2 posterior draws, mode (b) (the mode A8a chose)
    within the spec's tolerances, and the box truncation measured negligible."""
    from fork.export import m6c_repeat, m8
    from fork.inference.nested import NestedResult

    paths = [j.path(".pkl") for j in m6c_repeat.JOBS if j.group == "r2"]
    if not all(p.exists() for p in paths):
        pytest.skip("the R2 nested results are not in data/cache")
    _, _, space, prob = m8.problem("headline")
    eq = m6c_repeat.merged_cell2([NestedResult.load(p) for p in paths]).equal_weight(
        np.random.default_rng(31)
    )
    rng = np.random.default_rng(32)
    for q in eq[rng.choice(len(eq), 12, replace=False)]:
        b = prob.marginal(q[4:], "b")
        u_min, chi2_full, _ = prob.profile(q[4:], b.u_star)
        assert abs(b.chi2_min - chi2_full) < 0.1
        assert abs(b.m2lnl - prob.laplace(q[4:], u_min)) < 0.2
        box = space.uranus_box
        assert lz.truncation_log10_bound(b.u_star, b.cov_u, box.lo, box.hi) < -100
