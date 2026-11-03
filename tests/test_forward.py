"""M4: forward model θ → Uranus longitude."""

import math

import numpy as np
import pytest

from fork.data.frames import state_to_elements
from fork.data.horizons import state_vectors_ecl_j2000, states_about_inner_barycenter_ecl_j2000
from fork.data.synthetic import jd_tdb_jan1
from fork.inference.fit import wrap_rad
from fork.integrators import DEFAULT_DT_DAYS, wisdom_holman
from fork.model import bodies
from fork.model.forward import (
    PERTURBER_PARAM_NAMES,
    URANUS_PARAM_NAMES,
    ForwardContext,
    UranusModel,
    perturber_theta_from_state,
    predict_uranus_lon,
    uranus_inc_raan_at_epoch,
    uranus_theta_from_state,
)

ARC = math.degrees(1.0) * 3600.0
EPOCH = 2371628.5  # 1781-03-13.0 TDB
YEARS = np.arange(1690, 1847)

pytestmark = pytest.mark.usefixtures("offline_horizons")


@pytest.fixture(scope="module")
def setup(offline_horizons):
    ctx = ForwardContext.from_horizons(EPOCH)
    grid = jd_tdb_jan1(YEARS)
    truth = state_vectors_ecl_j2000("uranus", grid)
    lon_true = np.mod(np.arctan2(truth[:, 1], truth[:, 0]), 2 * math.pi)
    th_u = uranus_theta_from_state(states_about_inner_barycenter_ecl_j2000("uranus", [EPOCH])[0])
    th_n = perturber_theta_from_state(
        states_about_inner_barycenter_ecl_j2000("neptune", [EPOCH])[0], bodies.NEPTUNE.mass_msun
    )
    return ctx, grid, lon_true, th_u, th_n


def test_param_names_and_defaults(setup):
    ctx, *_ = setup
    model = UranusModel(ctx)
    assert model.param_names == URANUS_PARAM_NAMES and len(PERTURBER_PARAM_NAMES) == 5
    assert ctx.dt_days == DEFAULT_DT_DAYS and ctx.scheme == wisdom_holman.SCHEME_DEFAULT


def test_planted_real_neptune_reproduces_horizons_uranus(setup):
    """Forward model with Uranus + (planar) real Neptune from Horizons at 1781 must reproduce the
    Horizons Uranus longitude over 1690–1846. Measured max 0.55″, rms 0.12″ (the residual is
    mostly Neptune's dropped 1.8° inclination)."""
    ctx, grid, lon_true, th_u, th_n = setup
    r = wrap_rad(lon_true - predict_uranus_lon(np.concatenate([th_u, th_n]), grid, ctx)) * ARC
    assert np.max(np.abs(r)) < 1.0
    assert np.sqrt(np.mean(r**2)) < 0.25


def test_without_neptune_the_signal_is_large(setup):
    ctx, grid, lon_true, th_u, _ = setup
    r = wrap_rad(lon_true - UranusModel(ctx).predict(th_u, grid)) * ARC
    assert np.ptp(r) > 200.0  # measured 266″ before any fitting


def test_kernel_matches_direct_integration(setup):
    """Elements → state inside the kernel is exact: same as integrating the state directly."""
    ctx, grid, _, th_u, _ = setup
    names = ("jupiter", "saturn", "uranus")
    s = np.array([states_about_inner_barycenter_ecl_j2000(n, [EPOCH])[0] for n in names])
    masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in names])
    pos = wisdom_holman.integrate(
        s[:, :3].copy(), s[:, 3:].copy(), masses, grid - EPOCH, DEFAULT_DT_DAYS
    )
    lon = np.arctan2(pos[:, 2, 1], pos[:, 2, 0])
    assert np.max(np.abs(wrap_rad(lon - UranusModel(ctx).predict(th_u, grid)))) * ARC < 1e-6


def test_inner_barycenter_initial_conditions_matter(offline_horizons):
    """Regression for the M4 discovery: with the inner planets folded into the Sun, initial
    conditions must be relative to the Sun + inner-planet barycenter. Sun-centred ICs carry the
    Sun's ~7 cm/s reflex and shift Uranus by ~35″ over 1690–1846; barycentric ICs: < 0.2″."""
    grid = jd_tdb_jan1(YEARS)
    truth = state_vectors_ecl_j2000("uranus", grid)
    lon_true = np.arctan2(truth[:, 1], truth[:, 0])
    names = ("jupiter", "saturn", "uranus", "neptune")
    masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in names])
    err = {}
    for label, fn in (
        ("helio", state_vectors_ecl_j2000),
        ("bary", states_about_inner_barycenter_ecl_j2000),
    ):
        s = np.array([fn(n, [EPOCH])[0] for n in names])
        pos = wisdom_holman.integrate(
            s[:, :3].copy(), s[:, 3:].copy(), masses, grid - EPOCH, DEFAULT_DT_DAYS
        )
        err[label] = (
            np.max(np.abs(wrap_rad(lon_true - np.arctan2(pos[:, 2, 1], pos[:, 2, 0])))) * ARC
        )
    assert err["bary"] < 0.2  # measured 0.09″
    assert err["helio"] > 10.0  # measured 34.6″


def test_negligible_perturber_equals_no_perturber(setup):
    ctx, grid, _, th_u, th_n = setup
    tiny = th_n.copy()
    tiny[0] = -12.0  # 1e-12 Earth masses
    d = predict_uranus_lon(np.concatenate([th_u, tiny]), grid, ctx) - predict_uranus_lon(
        th_u, grid, ctx
    )
    assert np.max(np.abs(wrap_rad(d))) * ARC < 1e-6


def test_shapes_and_validation(setup):
    ctx, grid, _, th_u, _ = setup
    model = UranusModel(ctx)
    assert model.predict(th_u, grid[5]).shape == ()
    assert model.predict(th_u, grid[:6].reshape(2, 3)).shape == (2, 3)
    out = model.predict(th_u, grid)
    assert np.all((out >= 0) & (out < 2 * math.pi))
    with pytest.raises(ValueError):
        model.predict(np.zeros(5), grid)
    with pytest.raises(ValueError):
        predict_uranus_lon(np.zeros(7), grid, ctx)


def _lon_ecl(orbit_lon, inc, raan):
    """Ecliptic longitude of a point at orbital longitude Ω + u on an orbit with (i, Ω)."""
    u = orbit_lon - raan
    return raan + np.arctan2(np.cos(inc) * np.sin(u), np.cos(u))


def test_inc_raan_leak(offline_horizons):
    """project rule 2: holding Uranus's i, Ω at their 1781 values leaks < 1″ of Neptune.

    Integrate J+S+U with and without Neptune from the same 1781 ICs; Neptune's effect on Uranus's
    osculating (i, Ω) over 1781–1846 is converted to longitude as the worst case, over every
    orbital longitude, of the change in ecliptic longitude it would cause. Measured: |Δi| ≤ 2.4″,
    |ΔΩ| ≤ 64″, longitude-equivalent ≤ 0.017″ (at 1846)."""
    t_rel = np.append(np.arange(0.0, 2395296.5 - EPOCH, 30.0), 2395296.5 - EPOCH)
    mu = bodies.G * (bodies.M_CENTRAL_MSUN + bodies.URANUS.mass_msun)
    inc_raan = {}
    for label, names in (("without", ("jupiter", "saturn", "uranus")),
                         ("with", ("jupiter", "saturn", "uranus", "neptune"))):  # fmt: skip
        s = np.array([states_about_inner_barycenter_ecl_j2000(n, [EPOCH])[0] for n in names])
        masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in names])
        pos, vel = wisdom_holman.integrate_states(
            s[:, :3].copy(), s[:, 3:].copy(), masses, t_rel, DEFAULT_DT_DAYS
        )
        el = np.array([state_to_elements(pos[k, 2], vel[k, 2], mu) for k in range(t_rel.size)])
        inc_raan[label] = el[:, 2], el[:, 3]
    (i1, o1), (i0, o0) = inc_raan["with"], inc_raan["without"]
    orbit_lon = np.linspace(0.0, 2 * math.pi, 3601)[:, None]
    dlon = wrap_rad(_lon_ecl(orbit_lon, i1, o1) - _lon_ecl(orbit_lon, i0, o0))
    leak_arcsec = np.max(np.abs(dlon)) * ARC
    assert np.max(np.abs(i1 - i0)) * ARC < 10.0 and np.max(np.abs(wrap_rad(o1 - o0))) * ARC < 300.0
    assert leak_arcsec < 1.0
    assert leak_arcsec < 0.05  # measured 0.017″; a jump here means the leak test itself changed


def test_inc_raan_at_epoch_matches_horizons_state(offline_horizons):
    inc, raan = uranus_inc_raan_at_epoch(EPOCH)
    assert math.degrees(inc) == pytest.approx(0.777, abs=1e-3)
    assert math.degrees(raan) == pytest.approx(73.8, abs=0.1)
