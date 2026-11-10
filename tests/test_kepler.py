"""M2: universal-variable Kepler propagator."""

import math

import numpy as np
import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

from fork.data.frames import elements_to_state
from fork.integrators.kepler import (
    kepler_drift_inplace,
    kepler_fg_diag,
    kepler_step,
    stumpff_c2_c3,
)
from fork.model.bodies import M_CENTRAL_MSUN, G

MU = G * M_CENTRAL_MSUN
TWO_PI = 2.0 * math.pi


def state_alpha(r, v):
    """1/a implied by the (rounded) state itself — the quantity the propagator actually uses."""
    return 2.0 / np.linalg.norm(r) - (v @ v) / MU


def state_period(r, v):
    al = state_alpha(r, v)
    return TWO_PI / math.sqrt(MU * al**3)


def random_elements(rng, a, e):
    return (a, e, rng.uniform(0, math.pi), *rng.uniform(0, TWO_PI, 3))


# ------------------------------------------------------------------ Stumpff functions


@pytest.mark.parametrize(
    "psi", [-30.0, -1.0 - 1e-12, -1.0 + 1e-12, -0.3, 0.3, 1.0 - 1e-12, 1.0 + 1e-12, 5.0, 39.4]
)
def test_stumpff_matches_closed_form(psi):
    """Series branch (|ψ| < 1) and closed-form branch agree with the definitions across ψ = ±1."""
    c2, c3 = stumpff_c2_c3(psi)
    s = math.sqrt(abs(psi))
    if psi > 0:
        want2, want3 = (1 - math.cos(s)) / psi, (s - math.sin(s)) / s**3
    else:
        want2, want3 = (math.cosh(s) - 1) / (-psi), (math.sinh(s) - s) / s**3
    assert c2 == pytest.approx(want2, rel=2e-15)
    assert c3 == pytest.approx(want3, rel=2e-15)


def test_stumpff_at_zero():
    assert stumpff_c2_c3(0.0) == (0.5, 1.0 / 6.0)


# ------------------------------------------------------------------ full period


@pytest.mark.parametrize("a", [1.0, 5.2, 9.5, 19.2])
@pytest.mark.parametrize("e", [0.0, 1e-9, 0.05, 0.1, 0.3])
def test_full_period_returns_to_start(a, e):
    """Spec M2: one full period returns to start to < 1e-12 AU.

    Exercised three ways, all with the period implied by the state itself:
    a single Δt = P call (period-reduction path), two genuine half-orbit solves, and ten random
    sub-steps. Covers every planet and the M5 perturber prior (e ≤ 0.3); e ≥ 0.6 is below.
    """
    rng = np.random.default_rng(int(a * 100 + e * 1000))
    for _ in range(50):
        r0, v0 = elements_to_state(*random_elements(rng, a, e), MU)
        p = state_period(r0, v0)

        r, v = kepler_step(r0, v0, MU, p)
        assert np.linalg.norm(r - r0) < 1e-12

        r, v = kepler_step(r0, v0, MU, 0.5 * p)
        r, v = kepler_step(r, v, MU, 0.5 * p)
        assert np.linalg.norm(r - r0) < 1e-12
        assert np.linalg.norm(v - v0) < 1e-12 * np.linalg.norm(v0) * 10

        frac = rng.uniform(0.05, 1.0, 10)
        frac /= frac.sum()
        r, v = r0.copy(), v0.copy()
        for f in frac:
            r, v = kepler_step(r, v, MU, f * p)
        assert np.linalg.norm(r - r0) < 1e-12


@pytest.mark.parametrize("a", [1.0, 19.2, 40.0])
@pytest.mark.parametrize("e", [0.6, 0.8, 0.95])
def test_full_period_high_eccentricity(a, e):
    """e ≥ 0.6 misses 1e-12 AU at large a; this is double precision, not the solver.

    At pericenter 1/a = 2/r − v²/μ cancels ~(1+e)/(1−e) ≈ 40×, so every rounded state's energy,
    and thus its period, is uncertain at ~40 ε. Times v_peri·P ≈ 39a that is ~5e-12 AU at
    a = 19. Measured worst (300 trials): 1.3e-12 AU at e = 0.6, a = 19.2; 1.2e-11 AU at e = 0.95,
    a = 19.2; 7.8e-11 AU at e = 0.95, a = 40. The bound
    asserted here is 1e-11 · a; see docs/notes/M2.md.
    """
    rng = np.random.default_rng(int(a * 10 + e * 100))
    for _ in range(100):
        r0, v0 = elements_to_state(*random_elements(rng, a, e), MU)
        p = state_period(r0, v0)
        r, v = kepler_step(r0, v0, MU, 0.5 * p)
        r, v = kepler_step(r, v, MU, 0.5 * p)
        assert np.linalg.norm(r - r0) < 1e-11 * a


# ------------------------------------------------------------------ vs analytic ellipse


@settings(max_examples=1000, deadline=None)
# The one case that failed on CI's x86 runner (2026-11-10, run 42816309814) with the former
# fixed bound of 1e-12·a: 4.87e-12 = 1.14e-12·a, at pericentre (mean anomaly 0.002) of an
# e = 0.875 orbit after 8 periods. Pinned so every run checks it.
@example(a=4.265625, e=0.875, inc=0.0, raan=0.0, argp=1.0, mean=0.001953125, n_periods=8.0)
@given(
    a=st.floats(0.3, 100.0),
    e=st.floats(0.0, 0.95),
    inc=st.floats(0.0, math.pi),
    raan=st.floats(0.0, TWO_PI),
    argp=st.floats(0.0, TWO_PI),
    mean=st.floats(0.0, TWO_PI),
    n_periods=st.floats(-20.0, 20.0),
)
def test_agrees_with_analytic_ellipse(a, e, inc, raan, argp, mean, n_periods):
    """Spec M2: universal-variable solution == Kepler-equation ellipse, 1000 cases, e ∈ [0, 0.95].

    The reference ellipse uses a = 1/α of the rounded state (otherwise the test measures the
    rounding of the initial state rather than the propagator).

    The bound scales with 1/(1 − e). Both sides carry a mean anomaly of up to 20 revolutions
    (about 126 rad) through libm's sin and cos, whose last bit differs between platforms; a
    last-bit difference δM in the mean anomaly becomes a difference δE = δM / (1 − e cos E) in
    the eccentric anomaly (the conditioning of Kepler's equation), which is 1/(1 − e) at
    pericentre, and then a position difference of about a·δE. So the platform-dependent part
    of the discrepancy is up to a/(1 − e) times the rounding of M, not a times it. Measured
    worst on Apple Silicon: 1.2e-13·a over 20 000 random cases (6e-14·a/(1 − e)); on CI's x86
    runner the pinned example reached 1.14e-12·a = 1.4e-13·a/(1 − e). The bound 1e-12·a/(1 − e)
    is 7× above that case and 8× above the former bound at e = 0.875.
    """
    r0, v0 = elements_to_state(a, e, inc, raan, argp, mean, MU)
    a_state = 1.0 / state_alpha(r0, v0)
    n = math.sqrt(MU / a_state**3)
    dt = n_periods * TWO_PI / n
    r_ref, v_ref = elements_to_state(a_state, e, inc, raan, argp, mean + n * dt, MU)
    r, v = kepler_step(r0, v0, MU, dt)
    tol = 1e-12 * a / (1.0 - e)  # the conditioning of Kepler's equation near pericentre
    assert np.linalg.norm(r - r_ref) < tol
    assert np.linalg.norm(v - v_ref) < tol * n * 50  # v ≤ ~6 n a at e = 0.95 pericenter


def test_very_long_dt():
    """10⁴ periods: reduction keeps the solve as cheap and accurate as a short one."""
    rng = np.random.default_rng(3)
    r0, v0 = elements_to_state(*random_elements(rng, 19.2, 0.047), MU)
    a_state = 1.0 / state_alpha(r0, v0)
    n = math.sqrt(MU / a_state**3)
    dt = 1.0e4 * TWO_PI / n + 1234.5
    r, _v = kepler_step(r0, v0, MU, dt)
    # Reference: same elements advanced by n·Δt. Phase error ~ ε·n·Δt ≈ 1e-11 rad → ×a.
    _a0, e0, inc0, raan0, argp0, mean0 = random_elements(np.random.default_rng(3), 19.2, 0.047)
    r_ref, _ = elements_to_state(a_state, e0, inc0, raan0, argp0, mean0 + n * dt, MU)
    assert np.linalg.norm(r - r_ref) < 1e-9


# ------------------------------------------------------------------ invariants & edge cases


def test_conserves_energy_and_angular_momentum_and_reverses():
    rng = np.random.default_rng(5)
    for _ in range(500):
        a = rng.uniform(0.5, 60.0)
        r0, v0 = elements_to_state(*random_elements(rng, a, rng.uniform(0, 0.95)), MU)
        dt = rng.uniform(-3, 3) * state_period(r0, v0)
        r, v = kepler_step(r0, v0, MU, dt)
        e_0 = 0.5 * v0 @ v0 - MU / np.linalg.norm(r0)
        e_1 = 0.5 * v @ v - MU / np.linalg.norm(r)
        assert abs(e_1 - e_0) < 1e-11 * abs(e_0)
        h0, h1 = np.cross(r0, v0), np.cross(r, v)
        assert np.linalg.norm(h1 - h0) < 1e-12 * np.linalg.norm(h0)
        rb, _vb = kepler_step(r, v, MU, -dt)
        assert np.linalg.norm(rb - r0) < 1e-11 * a


def test_dt_zero_is_identity():
    r0 = np.array([19.0, 1.0, 0.1])
    v0 = np.array([-1e-4, 3.9e-3, 1e-5])
    r, v = kepler_step(r0, v0, MU, 0.0)
    np.testing.assert_array_equal(r, r0)
    np.testing.assert_array_equal(v, v0)


def test_exactly_circular_matches_rotation():
    a = 19.2
    vc = math.sqrt(MU / a)
    r0 = np.array([a, 0.0, 0.0])
    v0 = np.array([0.0, vc, 0.0])
    n = vc / a
    for dt in (20.0, 1000.0, 12345.6, -777.0):
        r, _ = kepler_step(r0, v0, MU, dt)
        want = a * np.array([math.cos(n * dt), math.sin(n * dt), 0.0])
        assert np.linalg.norm(r - want) < 1e-12 * a


def test_high_e_pericenter_short_step_converges():
    """Pericenter of e = 0.95 with a tiny step: the a-based guess is off by ×20; must converge."""
    for a in (1.0, 30.0):
        r0, v0 = elements_to_state(a, 0.95, 0.1, 0.2, 0.3, 0.0, MU)
        p = state_period(r0, v0)
        for frac in (1e-9, 1e-6, 1e-3, 0.01, 0.49, 0.5, 0.51):
            r, v = kepler_step(r0, v0, MU, frac * p)
            rb, _ = kepler_step(r, v, MU, -frac * p)
            assert np.linalg.norm(rb - r0) < 1e-11 * a


def test_hyperbolic_invariants():
    r0 = np.array([1.0, 0.0, 0.0])
    v0 = np.array([0.0, 2.0 * math.sqrt(MU), 0.3 * math.sqrt(MU)])  # well above escape speed
    assert state_alpha(r0, v0) < 0
    e_0 = 0.5 * v0 @ v0 - MU
    for dt in (1.0, 50.0, 500.0, -300.0):
        r, v = kepler_step(r0, v0, MU, dt)
        e_1 = 0.5 * v @ v - MU / np.linalg.norm(r)
        assert abs(e_1 - e_0) < 1e-12 * abs(e_0)
        h0 = np.cross(r0, v0)
        assert np.linalg.norm(np.cross(r, v) - h0) < 1e-12 * np.linalg.norm(h0)
        rb, _ = kepler_step(r, v, MU, -dt)
        np.testing.assert_allclose(rb, r0, rtol=0, atol=1e-11 * max(1.0, np.linalg.norm(r)))


def test_drift_inplace_matches_kepler_step():
    rng = np.random.default_rng(9)
    x = np.empty((4, 3))
    v = np.empty((4, 3))
    for i, a in enumerate((5.2, 9.5, 19.2, 30.1)):
        x[i], v[i] = elements_to_state(*random_elements(rng, a, 0.05), MU)
    want = [kepler_step(x[i], v[i], MU, 20.0) for i in range(4)]
    kepler_drift_inplace(x, v, MU, 20.0)
    for i in range(4):
        np.testing.assert_allclose(x[i], want[i][0], rtol=1e-13, atol=1e-15)
        np.testing.assert_allclose(v[i], want[i][1], rtol=1e-13, atol=1e-15)


def test_short_steps_converge_by_newton_without_bisection():
    """Regression: a converged Newton step landing on the bracket edge used to trigger a
    spurious bisection jump (12.8k bisections per 60k calls, even for tiny steps).

    WH hot path: |Δt| ≤ 5% of P (dt = 20 d is 0.5% of Jupiter's P), e ≤ 0.6 → pure Newton, ≤ 6
    iterations. At e ≈ 0.9 a 5%-of-P step through pericenter can genuinely overshoot; there the
    safeguard may act, but rarely.
    """
    rng = np.random.default_rng(11)
    worst_iter = 0
    hi_e_calls = hi_e_bisect = 0
    for _ in range(4000):
        a = rng.uniform(0.3, 100.0)
        e = rng.uniform(0, 0.95)
        r0, v0 = elements_to_state(*random_elements(rng, a, e), MU)
        p = state_period(r0, v0)
        for dt in (20.0, rng.uniform(-0.05, 0.05) * p):
            if abs(dt) > 0.05 * p:
                continue
            out = kepler_fg_diag(*r0, *v0, MU, dt)
            n_iter, n_bisect = out[6], out[7]
            if e <= 0.6:
                assert n_bisect == 0
                worst_iter = max(worst_iter, n_iter)
            else:
                hi_e_calls += 1
                hi_e_bisect += n_bisect > 0
                assert n_iter <= 12
    assert worst_iter <= 6
    assert hi_e_bisect < 0.01 * hi_e_calls


def test_long_steps_rarely_need_bisection():
    rng = np.random.default_rng(12)
    calls_with_bisection = 0
    for _ in range(2000):
        a = rng.uniform(0.3, 100.0)
        r0, v0 = elements_to_state(*random_elements(rng, a, rng.uniform(0, 0.95)), MU)
        out = kepler_fg_diag(*r0, *v0, MU, rng.uniform(-20, 20) * state_period(r0, v0))
        calls_with_bisection += out[7] > 0
        assert out[6] < 60
    assert calls_with_bisection < 0.05 * 2000  # measured ≈ 2.4%
