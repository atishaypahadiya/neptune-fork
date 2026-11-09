"""M8 amendment: the perturbative likelihood (PHASE1_SPEC, "M8 amendment", 5fc089f).

On the synthetic seed of tests/test_linearized.py (Table 77's layout, real Uranus with the real
Neptune, offsets and noise from the headline error model):
- Φ against direct finite-difference propagation, and its symplectic identity;
- the Kepler propagation against the closed-form elements;
- δx against the full integrator for a toy perturber, with the perturber's own trajectory
  from the integrator (isolating the response) and with the Kepler orbit;
- the second-order switch reducing the error at Neptune's mass;
- the quadrature's convergence in Δs;
- with the forcing off, the marginal is Stage 1's linear marginal at u_ref;
- A8a-style comparison against the full likelihood on a few points (a regression check: the
  spec's gate is the reweighting ESS, not these pointwise numbers).
"""

import math

import numpy as np
import pytest

from fork.model import linearized as lz
from fork.model import perturbative as pt


@pytest.fixture(scope="module")
def synthetic():
    from fork.export.to_json import A1_EPOCH_JD_TDB
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.inference.priors import free_a_space
    from fork.model import observation as om
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch
    from test_linearized import _neptune_phi

    structure = om.structure_from_rows(om.read_rows(), om.ErrorModel.load())
    vals, _ = om.synthetic_values(structure, om.true_uranus_ra_dec(structure),
                                  np.random.default_rng(808))  # fmt: skip
    obs = structure.with_values(vals)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    box = uranus_box_radec(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    space = free_a_space(pq, 0.3, box, inclined=True)
    lp = lz.LinearizedProblem.build(obs, ctx, space)
    phi0 = _neptune_phi(space, ctx)
    # no posterior for a synthetic seed: the spec's fallback, a full-likelihood point fit
    u_ref = lp.marginal(phi0, "b").u_star
    return obs, ctx, space, lp, phi0, pt.PerturbativeProblem.build(obs, ctx, space, u_ref)


def _true_delta_x(space, obs, ctx, u_ref, phi):
    from fork.model.forward import predict_states_ecl_j2000

    th = space.theta(np.concatenate([u_ref, phi]))
    pos, _ = predict_states_ecl_j2000(th, obs.geometry.jd_tt, ctx)
    pos0, _ = predict_states_ecl_j2000(th[:6], obs.geometry.jd_tt, ctx)
    return pos[:, 2] - pos0[:, 2]


def _rel(a, b):
    return float(np.abs(a - b).max() / np.abs(b).max())


def test_phi_matches_direct_propagation_and_is_symplectic(synthetic):
    from fork.inference.fit import uranus_x_to_theta
    from fork.integrators.wisdom_holman import integrate_states
    from fork.model.forward import _initial_conditions

    _, ctx, space, _, _, pr = synthetic
    assert pr.symplectic_dev < 1e-6
    x0, v0, masses = _initial_conditions(uranus_x_to_theta(np.concatenate([pr.u_ref, space.pq])),
                                         np.empty(0), ctx.x_js_ecl_j2000,
                                         ctx.v_js_ecl_j2000)  # fmt: skip
    t = np.ascontiguousarray(pr.nodes_jd - ctx.epoch_jd_tdb)
    p0, _ = integrate_states(x0, v0, masses, t, ctx.dt_days, ctx.scheme)
    rng = np.random.default_rng(3)
    n = np.linalg.norm(v0[2]) / np.linalg.norm(x0[2])
    d0 = rng.normal(size=6)
    d0 /= np.linalg.norm(d0)
    errs = []
    for amp in (1e-3, 5e-4):  # the same direction, two amplitudes
        dd = amp * d0 * np.r_[np.ones(3), np.full(3, n)]
        x1, v1 = x0.copy(), v0.copy()
        x1[2] += dd[:3]
        v1[2] += dd[3:]
        p1, _ = integrate_states(x1, v1, masses, t, ctx.dt_days, ctx.scheme)
        errs.append(_rel(np.einsum("kij,j->ki", pr.phi_pos, dd), p1[:, 2] - p0[:, 2]))
    # the linear prediction's error is the orbit's own nonlinearity: O(amplitude)
    assert errs[0] < 1e-3
    assert errs[1] / errs[0] == pytest.approx(0.5, abs=0.05)


def test_kepler_propagation_matches_closed_form_elements():
    from fork.data.frames import elements_to_state

    mu = 2.9e-4
    a, e, inc, raan, argp, m0 = 30.0, 0.15, 0.1, 1.0, 2.0, 0.5
    r0, v0 = elements_to_state(a, e, inc, raan, argp, m0, mu)
    n = math.sqrt(mu / a**3)
    dt = np.array([-40000.0, -1234.5, 0.0, 777.0, 60000.0])
    got = pt.kepler_positions(np.asarray(r0), np.asarray(v0), mu, dt)
    for k, d in enumerate(dt):
        r, _ = elements_to_state(a, e, inc, raan, argp, m0 + n * d, mu)
        np.testing.assert_allclose(got[k], r, rtol=0, atol=1e-10)


def test_delta_x_matches_full_integration_for_a_toy_perturber(synthetic):
    """A toy perturber at a = 30 AU, e = 0.01. At 0.01 M⊕ only first order survives: with the
    integrator's own perturber trajectory the left-out Jupiter/Saturn response (E_JS) is what
    remains, and adding it back leaves ~10⁻⁶. The Kepler orbit adds E_K. At 1 M⊕ the second
    order is visible and the switch removes it (measured: 5.9e-5 → 1.9e-6)."""
    obs, ctx, space, _, _, pr = synthetic
    for log10_m, first_bound in ((-2.0, 5e-6), (0.0, None)):
        toy = np.array([log10_m, math.log(30.0), 0.1, 0.0, 1.0, 0.0, 0.0])
        true = _true_delta_x(space, obs, ctx, pr.u_ref, toy)
        x_p, dx_j, dx_s = pr.integrated_states(toy)
        extra = pr.js_response_accel(dx_j, dx_s)
        tab = _rel(pr.delta_x(toy, perturber_pos=x_p), true)
        tab_js = _rel(pr.delta_x(toy, perturber_pos=x_p, extra_accel=extra), true)
        tab_js_so = _rel(pr.delta_x(toy, perturber_pos=x_p, extra_accel=extra,
                                    second_order=True), true)  # fmt: skip
        assert _rel(pr.delta_x(toy), true) < 3e-3  # E_K
        if first_bound is not None:
            assert tab_js < first_bound and tab > 10 * tab_js  # E_JS dominates the rest
        else:
            assert tab_js_so < tab_js / 10 and tab_js_so < 5e-6


def test_second_order_switch_reduces_the_error(synthetic):
    obs, ctx, space, _, phi0, pr = synthetic  # Neptune's own φ: 17 M⊕
    true = _true_delta_x(space, obs, ctx, pr.u_ref, phi0)
    x_p, dx_j, dx_s = pr.integrated_states(phi0)
    extra = pr.js_response_accel(dx_j, dx_s)
    first = _rel(pr.delta_x(phi0, perturber_pos=x_p, extra_accel=extra), true)
    second = _rel(pr.delta_x(phi0, perturber_pos=x_p, extra_accel=extra, second_order=True),
                  true)  # fmt: skip
    assert second < first / 3
    # and on the model as specified (Kepler, no J/S response)
    assert _rel(pr.delta_x(phi0, second_order=True), true) < _rel(pr.delta_x(phi0), true)


def test_quadrature_converges_in_the_grid_step(synthetic):
    obs, ctx, space, _, phi0, pr = synthetic
    fine = pt.PerturbativeProblem.build(obs, ctx, space, pr.u_ref, step_days=15.0)
    trap = pt.PerturbativeProblem.build(obs, ctx, space, pr.u_ref, quadrature="trapezoid")
    ref = fine.delta_x(phi0)
    assert _rel(pr.delta_x(phi0), ref) < 1e-7
    assert 1e-7 < _rel(trap.delta_x(phi0), ref) < 1e-4  # second order: visible, still small


def test_interval_weights_are_exact_for_cubics():
    rng = np.random.default_rng(1)
    s = np.sort(np.concatenate([np.arange(0.0, 600.0, 60.0), rng.uniform(0, 540, 7)]))
    idx, w = pt.interval_weights(s, "cubic")
    f = lambda x: 2.0 - x + 3e-3 * x**2 - 4e-6 * x**3  # noqa: E731
    big_f = lambda x: 2.0 * x - x**2 / 2 + 1e-3 * x**3 - 1e-6 * x**4  # noqa: E731
    got = np.sum(w * f(s[idx]), axis=1)
    np.testing.assert_allclose(got, big_f(s[1:]) - big_f(s[:-1]), rtol=1e-9, atol=1e-9)


def test_no_forcing_gives_the_stage_one_marginal_at_the_reference(synthetic):
    from fork.inference.fit import uranus_x_to_theta

    obs, ctx, space, _, _, pr = synthetic
    got = pr.marginal_from_z(pr.z_from_delta(np.zeros((pr.obs_node.size, 3))))
    z0 = obs.whitened(uranus_x_to_theta(np.concatenate([pr.u_ref, space.pq])), ctx)
    want = lz.linear_marginal(z0, pr.jac_u, pr.u_ref, space.uranus_box.log_volume)
    # positions from the tabulated reference vs a fresh integration to the epochs
    assert got.m2lnl == pytest.approx(want.m2lnl, abs=1e-6)
    np.testing.assert_allclose(got.u_star, want.u_star, rtol=0, atol=1e-9)


def test_student_t_observations_are_refused(synthetic):
    from fork.model.observation import GROUPS
    from fork.model.student_t import StudentTObservations

    obs, ctx, space, _, _, pr = synthetic
    t_obs = StudentTObservations.from_gaussian(obs, 4.0, dict.fromkeys(GROUPS, 1.0))
    with pytest.raises(TypeError, match="Student-t"):
        pt.PerturbativeProblem.build(t_obs, ctx, space, pr.u_ref)


def test_against_the_full_likelihood_near_the_truth(synthetic):
    """A regression check on the measured size (docs/notes/M8.md), not the spec's gate."""
    _, _, _, lp, phi0, pr = synthetic
    rng = np.random.default_rng(11)
    jitter = np.array([0.02, 0.01, 0.01, 0.01, 0.002, 0.001, 0.001])
    for k in range(3):
        phi = phi0 + (0.0 if k == 0 else 1.0) * rng.normal(0.0, 1.0, 7) * jitter
        b = lp.marginal(phi, "b")
        u_min, chi2_full, _ = lp.profile(phi, b.u_star)
        lap = lp.laplace(phi, u_min)
        for so in (False, True):
            r = pr.evaluate(phi, second_order=so)
            assert abs(r.marginal.chi2_min - chi2_full) < 1.0
            assert abs(r.marginal.m2lnl - lap) < 1.0
            assert not r.outside_box and r.truncation_log10 < -6
        assert pr.log_likelihood(phi) == pytest.approx(-0.5 * pr.evaluate(phi).marginal.m2lnl,
                                                       abs=1e-9)  # fmt: skip
