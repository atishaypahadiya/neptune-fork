"""M3: integrators, including acceptance test A1 (vs REBOUND, energy, runtime)."""

import itertools
import json
import math
import time
from pathlib import Path

import numpy as np
import pytest

from fork.data.frames import elements_to_state
from fork.integrators import DEFAULT_DT_DAYS, leapfrog, rk4, wisdom_holman
from fork.integrators._nbody import energy_helio
from fork.integrators.kepler import kepler_step
from fork.integrators.wisdom_holman import (
    SCHEME_DEFAULT,
    SCHEME_DRIFT_OUTER,
    SCHEME_KICK_OUTER,
    SCHEMES,
)
from fork.model import bodies

RAD_TO_ARCSEC = math.degrees(1.0) * 3600.0
DT = DEFAULT_DT_DAYS  # 14 d
FIXTURE = Path(__file__).parent / "data" / "horizons_a1_1781-03-13.json"
A1_BODIES = ("jupiter", "saturn", "uranus")
I_URANUS = 2
JD_1846_END = 2395296.5  # 1846-12-31.0 TDB


def a1_initial_conditions():
    """Sun + Jupiter + Saturn + Uranus, heliocentric ecl_j2000, at JD 1781-03-13 (TDB)."""
    data = json.loads(FIXTURE.read_text())
    st = np.array([data["state_ecl_j2000"][n] for n in A1_BODIES])
    masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in A1_BODIES])
    return st[:, :3].copy(), st[:, 3:].copy(), masses, data["jd_tdb"]


def lon_ecl_j2000_rad(x_h: np.ndarray) -> np.ndarray:
    """Unwrapped heliocentric ecliptic-J2000 longitude of positions [..., 3]."""
    return np.unwrap(np.arctan2(x_h[..., 1], x_h[..., 0]))


def a1_output_times(jd0: float, step: float = DT) -> np.ndarray:
    """1781-03-13 → end of 1846 on the step grid (multiples of `step`, so no partial steps)."""
    return np.arange(0.0, JD_1846_END - jd0, step)


def test_defaults():
    assert DEFAULT_DT_DAYS == 14.0
    assert SCHEME_DEFAULT == SCHEME_DRIFT_OUTER == SCHEMES["wh_drift_outer"]
    assert SCHEMES["wh_kick_outer"] == SCHEME_KICK_OUTER


@pytest.fixture(scope="module")
def rebound_whfast_uranus():
    """REBOUND WHFast, democratic heliocentric, same ICs and dt: Uranus heliocentric positions."""
    rebound = pytest.importorskip("rebound")
    x0, v0, masses, jd0 = a1_initial_conditions()
    t_out = a1_output_times(jd0)
    sim = rebound.Simulation()
    sim.G = bodies.G
    sim.add(m=masses[0])
    for i in range(3):
        sim.add(m=masses[i + 1], x=x0[i, 0], y=x0[i, 1], z=x0[i, 2],
                vx=v0[i, 0], vy=v0[i, 1], vz=v0[i, 2])  # fmt: skip
    sim.integrator = "whfast"
    sim.integrator.coordinates = "democraticheliocentric"
    sim.dt = DT
    out = np.empty((t_out.size, 3))
    for k, t in enumerate(t_out):
        sim.integrate(t, exact_finish_time=1)
        assert sim.t == t
        out[k] = np.array(sim.particles[3].xyz) - np.array(sim.particles[0].xyz)
    return t_out, out


# ------------------------------------------------------------------ A1: agreement with REBOUND


def test_a1_default_wh_matches_rebound_whfast(rebound_whfast_uranus):
    """A1: Uranus heliocentric longitude, default WH vs REBOUND WHFast-DH, same ICs and dt,
    1781–1846: < 0.01″.

    The default composition is REBOUND's, so agreement is at roundoff (measured ~5e-9″). The
    second, much tighter assertion makes any real divergence (sign, mass, μ, ordering) fail loudly.
    """
    t_out, uranus_rb = rebound_whfast_uranus
    x0, v0, masses, _ = a1_initial_conditions()
    ours = wisdom_holman.integrate(x0, v0, masses, t_out, DT)[:, I_URANUS]  # default scheme
    dlon = (lon_ecl_j2000_rad(ours) - lon_ecl_j2000_rad(uranus_rb)) * RAD_TO_ARCSEC
    assert np.max(np.abs(dlon)) < 0.01
    assert np.max(np.abs(dlon)) < 1e-6


def test_kick_outer_differs_from_rebound_by_truncation_error_only(rebound_whfast_uranus):
    """wh_kick_outer is a different 2nd-order map: it differs from REBOUND by the sum of the two
    truncation errors (~0.024″ at dt = 14 d), not by anything larger."""
    t_out, uranus_rb = rebound_whfast_uranus
    x0, v0, masses, _ = a1_initial_conditions()
    ours = wisdom_holman.integrate(x0, v0, masses, t_out, DT, SCHEME_KICK_OUTER)[:, I_URANUS]
    dlon = (lon_ecl_j2000_rad(ours) - lon_ecl_j2000_rad(uranus_rb)) * RAD_TO_ARCSEC
    assert 0.005 < np.max(np.abs(dlon)) < 0.05


@pytest.mark.parametrize("scheme", [SCHEME_DRIFT_OUTER, SCHEME_KICK_OUTER], ids=list(SCHEMES))
def test_wh_convergence_order_is_two(scheme):
    """Permanent convergence-order test, both compositions.

    Uranus heliocentric position error over 1781–1846 vs an independent dt = 1 d reference
    (RK4: 4th order, its own error ~1e-10″) must drop ≈ 4× per halving of dt:
    ratio ∈ [3.5, 4.5] for dt = 28 → 14 → 7 → 3.5 d. Measured: 4.00, 4.00, 4.00 for both.
    """
    x0, v0, masses, jd0 = a1_initial_conditions()
    t_out = a1_output_times(jd0, step=28.0)  # on every grid below
    ref = rk4.integrate(x0, v0, masses, t_out, 1.0)[:, I_URANUS]
    errs = []
    for dt in (28.0, 14.0, 7.0, 3.5):
        pos = wisdom_holman.integrate(x0, v0, masses, t_out, dt, scheme)[:, I_URANUS]
        errs.append(np.max(np.linalg.norm(pos - ref, axis=1)))
    for coarse, fine in itertools.pairwise(errs):
        assert 3.5 <= coarse / fine <= 4.5


# ------------------------------------------------------------------ A1: energy over 1000 yr


def _energy_errors(integrate_states, x0, v0, masses, t, *args):
    e0 = energy_helio(x0, v0, masses)
    xs, vs = integrate_states(x0, v0, masses, t, DT, *args)
    return np.array([(energy_helio(xs[k], vs[k], masses) - e0) / abs(e0) for k in range(t.size)])


@pytest.fixture(scope="module")
def energy_1000yr():
    """Relative energy error at *every* step over 1000 yr (no sampling gaps)."""
    x0, v0, masses, _ = a1_initial_conditions()
    t = np.arange(0.0, 1000 * 365.25, DT)
    return t, {
        "wh_drift_outer": _energy_errors(wisdom_holman.integrate_states, x0, v0, masses, t),
        "wh_kick_outer": _energy_errors(
            wisdom_holman.integrate_states, x0, v0, masses, t, SCHEME_KICK_OUTER
        ),
        "rk4": _energy_errors(rk4.integrate_states, x0, v0, masses, t),
        "leapfrog": _energy_errors(leapfrog.integrate_states, x0, v0, masses, t),
    }


def test_a1_wh_energy_below_1e8(energy_1000yr):
    """A1: default WH relative energy error stays < 1e-8 over 1000 yr at dt = 14 d.

    Measured max 9.93e-9 (every step, 26,090 samples) — a thin 0.7% margin: the WH error scales
    as dt², so this test pins dt ≤ 14 d. REBOUND WHFast-DH behaves identically."""
    _, d_e = energy_1000yr
    assert np.max(np.abs(d_e["wh_drift_outer"])) < 1e-8


@pytest.mark.parametrize("key", ["wh_drift_outer", "wh_kick_outer", "leapfrog"])
def test_symplectic_energy_error_is_bounded(energy_1000yr, key):
    """No secular growth: the last century's max |ΔE/E| is no larger than the first's (±20%),
    and the fitted linear trend over 1000 yr is < 1% of the oscillation amplitude."""
    t, d_e = energy_1000yr
    de = d_e[key]
    n = t.size // 10
    first, last = np.max(np.abs(de[:n])), np.max(np.abs(de[-n:]))
    assert last < 1.2 * first
    slope_per_kyr = np.polyfit(t / 365250.0, de, 1)[0]
    assert abs(slope_per_kyr) < 0.01 * np.max(np.abs(de))


def test_kick_outer_energy_error_magnitude(energy_1000yr):
    _, d_e = energy_1000yr
    assert np.max(np.abs(d_e["wh_kick_outer"])) < 2.5e-8  # measured 2.0e-8 at dt = 14 d


def test_rk4_energy_drifts_secularly(energy_1000yr):
    """RK4 is not symplectic: its energy error grows ~linearly (A1 comparison panel)."""
    t, d_e = energy_1000yr
    de = d_e["rk4"]
    n = t.size // 10
    assert np.max(np.abs(de[-n:])) > 5 * np.max(np.abs(de[:n]))
    assert np.max(np.abs(de)) > 3 * np.max(np.abs(d_e["wh_drift_outer"]))


# ------------------------------------------------------------------ A1: runtime


def test_a1_wh_65yr_runtime_under_10ms():
    """A1: a 65-yr WH run (Sun + 3 planets, default dt and scheme) takes < 10 ms after JIT."""
    x0, v0, masses, _ = a1_initial_conditions()
    t_out = np.arange(0.0, 65 * 365.25, DT)
    wisdom_holman.integrate(x0, v0, masses, t_out, DT)  # warm-up / JIT
    best = math.inf
    for _ in range(20):
        t0 = time.perf_counter()
        wisdom_holman.integrate(x0, v0, masses, t_out, DT)
        best = min(best, time.perf_counter() - t0)
    assert best < 0.010  # measured ~1.6 ms


# ------------------------------------------------------------------ general integrator behaviour


def test_wh_with_massless_planets_is_exact_kepler():
    """With planet masses → 0 the kick and Sun-momentum terms vanish: WH must equal Kepler."""
    mu = bodies.G * bodies.M_CENTRAL_MSUN
    x0 = np.empty((2, 3))
    v0 = np.empty((2, 3))
    x0[0], v0[0] = elements_to_state(5.2, 0.048, 0.02, 1.7, 4.8, 0.3, mu)
    x0[1], v0[1] = elements_to_state(19.2, 0.047, 0.013, 1.29, 1.68, 2.5, mu)
    masses = np.array([bodies.M_CENTRAL_MSUN, 0.0, 0.0])
    t_out = np.array([0.0, 20.0, 1000.0, 7777.7, -3000.0])
    for scheme in (SCHEME_KICK_OUTER, SCHEME_DRIFT_OUTER):
        pos = wisdom_holman.integrate(x0, v0, masses, t_out, DT, scheme)
        for j, t in enumerate(t_out):
            for i in range(2):
                want, _ = kepler_step(x0[i], v0[i], mu, t)
                assert np.linalg.norm(pos[j, i] - want) < 1e-10


@pytest.mark.parametrize(
    "integrate_states",
    [
        lambda *a: wisdom_holman.integrate_states(*a, SCHEME_KICK_OUTER),
        lambda *a: wisdom_holman.integrate_states(*a, SCHEME_DRIFT_OUTER),
        leapfrog.integrate_states,
    ],
    ids=["wh_kick_outer", "wh_drift_outer", "leapfrog"],
)
def test_time_reversible(integrate_states):
    """Forward 65 yr, then integrate the end state back 65 yr: returns to the start."""
    x0, v0, masses, _ = a1_initial_conditions()
    span = 1187 * DT
    xs, vs = integrate_states(x0, v0, masses, np.array([span]), DT)
    xb, _ = integrate_states(xs[0], vs[0], masses, np.array([-span]), DT)
    assert np.max(np.linalg.norm(xb[0] - x0, axis=1)) < 1e-10


def test_output_times_do_not_change_trajectory():
    """Extra, unsorted, off-grid and negative outputs leave on-grid results bit-identical."""
    x0, v0, masses, _ = a1_initial_conditions()
    grid = np.array([0.0, 400.0, 4000.0])
    base = wisdom_holman.integrate(x0, v0, masses, grid, DT)
    t_mixed = np.array([4000.0, 133.3, -500.0, 0.0, 3999.99, 400.0, -20.0])
    mixed = wisdom_holman.integrate(x0, v0, masses, t_mixed, DT)
    np.testing.assert_array_equal(mixed[3], base[0])
    np.testing.assert_array_equal(mixed[5], base[1])
    np.testing.assert_array_equal(mixed[0], base[2])
    np.testing.assert_array_equal(mixed[3], x0)
    # an off-grid output is close to its on-grid neighbour, and backward output is meaningful
    assert np.linalg.norm(mixed[4] - mixed[0]) < 1e-3
    back = wisdom_holman.integrate(x0, v0, masses, np.array([-20.0]), DT)
    np.testing.assert_array_equal(mixed[6], back[0])


def test_integrators_agree_on_uranus():
    """WH (dt = 20 d), RK4 (dt = 5 d) and leapfrog (dt = 2 d) all describe the same orbit."""
    x0, v0, masses, _ = a1_initial_conditions()
    t_out = np.arange(0.0, 65 * 365.25, 400.0)
    ref = lon_ecl_j2000_rad(rk4.integrate(x0, v0, masses, t_out, 5.0)[:, I_URANUS])
    wh = lon_ecl_j2000_rad(wisdom_holman.integrate(x0, v0, masses, t_out, DT)[:, I_URANUS])
    lf = lon_ecl_j2000_rad(leapfrog.integrate(x0, v0, masses, t_out, 2.0)[:, I_URANUS])
    assert np.max(np.abs(wh - ref)) * RAD_TO_ARCSEC < 0.05
    assert np.max(np.abs(lf - ref)) * RAD_TO_ARCSEC < 1.0


@pytest.mark.parametrize(("module", "order"), [(rk4, 4), (leapfrog, 2)])
def test_convergence_order(module, order):
    """Self-convergence: errors vs a dt/8 run shrink by ≈ 2^order per halving of dt."""
    x0, v0, masses, _ = a1_initial_conditions()
    t_out = np.array([20 * 365.25])
    ref = module.integrate(x0, v0, masses, t_out, 1.25)[0, 0]  # Jupiter: the fast body
    errs = [np.linalg.norm(module.integrate(x0, v0, masses, t_out, h)[0, 0] - ref)
            for h in (20.0, 10.0)]  # fmt: skip
    assert errs[0] / errs[1] == pytest.approx(2**order, rel=0.15)


def test_input_validation():
    x0, v0, masses, _ = a1_initial_conditions()
    with pytest.raises(ValueError):
        wisdom_holman.integrate(x0, v0, masses[:-1], np.array([1.0]), DT)
    with pytest.raises(ValueError):
        wisdom_holman.integrate(x0, v0, masses, np.array([1.0]), -DT)
    with pytest.raises(ValueError):
        wisdom_holman.integrate(x0, v0, masses, np.array([1.0]), DT, 7)


def test_horizons_a1_fixture_matches_cache_format():
    data = json.loads(FIXTURE.read_text())
    assert data["refplane"] == "ecliptic" and data["center"] == "500@10"
    assert data["jd_tdb"] == 2371628.5
    for name in (*A1_BODIES, "neptune"):
        r = np.linalg.norm(data["state_ecl_j2000"][name][:3])
        assert 4.9 < r < 31.0
