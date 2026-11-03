"""M5a diagnostics: perturber spaces (e cap, held orbits), the 1846-longitude profile likelihood,
and Le Verrier's published solution."""

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from fork.data.synthetic import LongitudeObservations, synthetic_epochs
from fork.export.to_json import _profile_minimum
from fork.inference.fit import (
    PerturberSpace,
    fit_perturber_multistart,
    held_orbit_space,
    le_verrier_space,
    uranus_theta_to_x,
    wrap_rad,
)
from fork.inference.profile import (
    GalleGeometry,
    _run_constrained,
    _start_vectors,
    delta_chi2_intervals,
    galle_lon_of_date_rad,
)
from fork.model.forward import ForwardContext, predict_uranus_lon
from fork.model.hidden_mass import LE_VERRIER_1846, perturber_geocentric_longitude

ARC = math.degrees(1.0) * 3600.0
EPOCH = 2371628.5
ROOT = Path(__file__).resolve().parents[1]
PROFILE_JSON = ROOT / "site" / "data" / "m5a_profile.json"
N_WORKERS = min(12, os.cpu_count() or 1)
# Same model-generated truth as test_perturber_fit: Uranus-like + 20 M⊕ at 32 AU, e = 0.05.
THETA_MODEL = np.array([19.21, 0.046, math.radians(0.77), math.radians(74.0), math.radians(170.0),
                        math.radians(96.0), math.log10(20.0), 32.0, 0.05, 1.0, 3.3])  # fmt: skip

pytestmark = pytest.mark.usefixtures("offline_horizons")


@pytest.fixture(scope="module")
def ctx(offline_horizons):
    return ForwardContext.from_horizons(EPOCH)


@pytest.fixture(scope="module")
def model_obs(ctx):
    jd, sigma, source = synthetic_epochs()
    return LongitudeObservations(jd, predict_uranus_lon(THETA_MODEL, jd, ctx), sigma / ARC, source)


# ---------------------------------------------------------------------- perturber spaces


def test_perturber_space_free_held_and_e_cap():
    sp = le_verrier_space(38.4)
    assert sp.n_free == 4 and sp.held[1] == 38.4
    th = sp.theta_from_z(np.array([1.3, 0.3, -0.4, 2.0]))
    assert th[1] == 38.4 and th[2] == pytest.approx(0.25) and th[4] == pytest.approx(2.0)
    assert th[3] == pytest.approx(math.atan2(-0.4, 0.3) % (2 * math.pi))
    y = np.array([1.3, 30.0, 0.2, 0.1, 1.0])
    assert np.array_equal(PerturberSpace().y_from_z(PerturberSpace().z_from_y(y)), y)

    capped = le_verrier_space(None, e_max=0.05)
    small = capped.theta_from_z(np.array([1.0, 30.0, 0.05, 0.0, 1.0]))[2]
    assert small == pytest.approx(0.0025, rel=0.03)  # ≈ k² + h² well inside the cap
    for r in (0.3, 1.0, 10.0):
        assert capped.theta_from_z(np.array([1.0, 30.0, r, 0.0, 1.0]))[2] <= 0.05

    held = held_orbit_space(1.5, 36.15, 0.10761, math.radians(286.9))
    assert held.n_free == 1
    th = held.theta_from_z(np.array([4.0]))
    assert th[:4] == pytest.approx([1.5, 36.15, 0.10761, math.radians(286.9)])


def test_held_orbit_mode_recovers_lam0_and_uranus(ctx, model_obs):
    """Model-generated data, perturber m, a, e, ϖ held at truth: λ₀ from a start 20° off."""
    th = THETA_MODEL
    space = held_orbit_space(th[6], th[7], th[8], th[9])
    start = np.array([th[10] + math.radians(20.0)])
    r = fit_perturber_multistart(model_obs, ctx, (th[2], th[3]), space=space, mode="held",
                                 lam0_starts_rad=start)  # fmt: skip
    f = r.fits[0]
    assert r.mode == "held" and f.dof == model_obs.jd_tdb.size - 5
    assert f.status == "converged" and f.chi2 < 1e-6
    assert abs(wrap_rad(f.theta[10] - th[10])) < 1e-5 and abs(f.theta[0] - th[0]) < 1e-5


@pytest.mark.slow
def test_e_cap_recovers_model_perturber_inside_cap(ctx, model_obs):
    """~30 s: the start is 10° off and the capped (k, h) surface takes ~700 LM steps."""
    th = THETA_MODEL  # e = 0.05, inside a 0.1 cap
    start = np.array([th[10] - math.radians(10.0)])
    r = fit_perturber_multistart(model_obs, ctx, (th[2], th[3]), e_max=0.1, lam0_starts_rad=start)
    f = r.fits[0]
    # measured χ² 6e-14 (a tanh cap, tried first, stalled on its plateau at e = 0.09999)
    assert f.status == "converged" and f.chi2 < 1e-4
    assert abs(f.theta[8] - th[8]) < 1e-3 and abs(f.theta[7] - th[7]) < 0.05


def test_le_verrier_published_solution_constants():
    """Sources' internal consistency: mean longitude 318°47′ + equation of centre for e = 0.1076
    gives the published true longitude 326°32′ (to the elements' rounding)."""
    lv = LE_VERRIER_1846
    assert (lv.a_au, lv.e, lv.mass_msun) == (36.15, 0.10761, 0.00011)
    m = math.radians(lv.mean_lon_ecl_of_date_deg - lv.varpi_ecl_of_date_deg)
    ecc_anom = m
    for _ in range(20):  # Kepler's equation, Newton
        ecc_anom -= (ecc_anom - lv.e * math.sin(ecc_anom) - m) / (1 - lv.e * math.cos(ecc_anom))
    nu = 2 * math.atan2(math.sqrt(1 + lv.e) * math.sin(ecc_anom / 2),
                        math.sqrt(1 - lv.e) * math.cos(ecc_anom / 2))  # fmt: skip
    true_lon = (lv.varpi_ecl_of_date_deg + math.degrees(nu)) % 360
    assert true_lon == pytest.approx(lv.true_lon_ecl_of_date_deg, abs=0.1)
    # 1847 equinox → J2000: general precession ≈ 50.3″/yr × 153 yr ≈ 2.14°
    assert math.degrees(lv.varpi_ecl_j2000_rad()) - lv.varpi_ecl_of_date_deg == pytest.approx(
        2.137, abs=0.01
    )
    assert 10**lv.log10_m_earth == pytest.approx(36.6, abs=0.1)


# ------------------------------------------------------------------- profile likelihood


def test_delta_chi2_intervals_and_profile_minimum():
    x = np.arange(300.0, 351.0)
    parab = ((x - 320.3) / 2.0) ** 2  # Δχ² = 1 at 318.3 and 322.3
    (iv,) = delta_chi2_intervals(x, parab)
    assert iv == pytest.approx((318.3, 322.3), abs=0.1)  # linear interpolation on a 1° grid
    m = _profile_minimum(x, parab)
    assert m["lon_ecl_of_date_deg"] == pytest.approx(320.3, abs=1e-6) and not m["at_grid_edge"]
    two = np.minimum(((x - 305) / 2) ** 2, ((x - 340) / 2) ** 2 + 0.5)
    ivs = delta_chi2_intervals(x, two)
    assert len(ivs) == 2 and ivs[0][0] < 305 < ivs[0][1] and ivs[1][0] < 340 < ivs[1][1]
    edge = delta_chi2_intervals(x, (x - 350.0) ** 2 / 100)
    assert edge[0][1] == 350.0 and _profile_minimum(x, (x - 350.0) ** 2)["at_grid_edge"]


def test_constrained_fit_meets_constraint_and_recovers_truth(ctx, model_obs):
    """Fixed-a profile point at the model perturber's own 1846 longitude recovers it (χ² ≈ 0);
    one 4° away is met to ≪ 1″ (checked independently) at a higher χ². Only Δχ² = 1.84 even
    though these data are noise-free and a is held at truth: m, e and ϖ absorb most of a 4°
    shift (the degeneracy behind the M5a longitude scatter)."""
    th = THETA_MODEL
    geo = GalleGeometry.from_horizons()
    lon_true = galle_lon_of_date_rad(th, ctx, geo)
    space = le_verrier_space(th[7])
    x = uranus_theta_to_x(th[:6])
    xu, pq = x[:4] + np.array([1e-3, 2e-3, -1e-3, 1e-4]), x[4:]  # Uranus start slightly off
    chi2 = {}
    for dlon_deg in (0.0, 4.0):
        target = lon_true + math.radians(dlon_deg)
        z0, lam0 = _start_vectors(space, xu, ctx, geo, target, 1)[0]
        f = _run_constrained(z0, lam0, pq, space, model_obs, ctx, geo, target)
        assert f.success and f.constraint_miss_arcsec < 1e-6
        g = perturber_geocentric_longitude(f.theta, ctx)
        assert abs(wrap_rad(g.lon_ecl_of_date_rad - target)) * ARC < 1e-5
        chi2[dlon_deg] = f.chi2
    assert chi2[0.0] < 1e-3 and chi2[4.0] > 1.0  # measured 7e-15, 1.84


def test_committed_profile_json_schema():
    d = json.loads(PROFILE_JSON.read_text())
    lo, hi, step = d["grid_deg"]
    n_grid = round((hi - lo) / step) + 1
    assert (lo, hi, step) == (300.0, 350.0, 1.0) and d["n_starts_per_point"] >= 8
    assert d["true_neptune_lon_ecl_of_date_deg"] == pytest.approx(325.873, abs=1e-3)
    seeds = [ds["seed"] for ds in d["datasets"]]
    assert seeds[:4] == [0, 1, 2, 3] and None in seeds
    for ds in d["datasets"]:
        for mode in ("le_verrier", "free_a"):
            m = ds["modes"][mode]
            p = m["profile"]
            assert len(p["lon_ecl_of_date_deg"]) == n_grid == len(p["dchi2"]) == len(p["a_au"])
            assert min(p["dchi2"]) >= -1e-3 and m["dchi2_at_true_neptune"] >= -1e-3
            assert m["max_constraint_miss_arcsec"] < 1e-3
            assert m["overall_min"]["chi2"] <= m["global_min"]["chi2"] + 1e-6
            for a, b in m["dchi2_le_1_intervals_deg"]:
                assert lo <= a <= b <= hi or a == b
            if mode == "le_verrier":
                assert set(p["a_au"]) == {d["le_verrier_a_au"]}


@pytest.mark.slow
def test_committed_profile_point_reproduces(ctx):
    """Seed 0, Le Verrier mode, at true Neptune's longitude: a fresh 8-start constrained fit
    reproduces the committed profile value (all 8 starts agree on it; measured 82.475). ~10 s."""
    from fork.data.synthetic import make_synthetic_observations
    from fork.inference.fit import fit_uranus_elements, run_pool
    from fork.model.forward import uranus_inc_raan_at_epoch
    from fork.model.hidden_mass import BODE_A_AU, true_neptune_geocentric_longitude

    committed = json.loads(PROFILE_JSON.read_text())["datasets"][0]
    assert committed["seed"] == 0
    obs, _ = make_synthetic_observations(np.random.default_rng(0))
    geo = GalleGeometry.from_horizons()
    target = true_neptune_geocentric_longitude().lon_ecl_of_date_rad
    space = le_verrier_space(BODE_A_AU)
    x = uranus_theta_to_x(fit_uranus_elements(obs, ctx, uranus_inc_raan_at_epoch(EPOCH)).theta)
    jobs = [(z, lam0, x[4:], space, obs, ctx, geo, target)
            for z, lam0 in _start_vectors(space, x[:4], ctx, geo, target, 8)]  # fmt: skip
    best = min(f.chi2 for f in run_pool(_run_constrained, jobs, N_WORKERS) if f.success)
    assert best == pytest.approx(committed["modes"]["le_verrier"]["at_true_neptune"]["chi2"],
                                 abs=0.01)  # fmt: skip


# ------------------------------------ slow: diagnostics (docs/notes/M5a_diagnostics.json)


@pytest.fixture(scope="module")
def seed0_obs(ctx):
    from fork.data.synthetic import make_synthetic_observations

    return make_synthetic_observations(np.random.default_rng(0))[0]


def _galle_dlon_deg(theta, ctx):
    from fork.model.hidden_mass import true_neptune_geocentric_longitude

    g = perturber_geocentric_longitude(theta, ctx)
    nep = true_neptune_geocentric_longitude()
    return math.degrees(wrap_rad(g.lon_ecl_of_date_rad - nep.lon_ecl_of_date_rad))


@pytest.mark.slow
def test_le_verrier_published_solution_on_seed0(ctx, seed0_obs):
    """Le Verrier's 31 Aug 1846 m, a, e, ϖ held, only λ₀ (24 starts) and Uranus free. Measured
    χ² 330.44 / 82 and 1846-09-23 longitude +8.04° from Neptune: with his mass the data put the
    planet ~8° ahead of where he did (fitted 1847.0 mean longitude 327.1° vs his 318°47′). ~10 s."""
    from fork.model.forward import uranus_inc_raan_at_epoch

    lv = LE_VERRIER_1846
    space = held_orbit_space(lv.log10_m_earth, lv.a_au, lv.e, lv.varpi_ecl_j2000_rad())
    r = fit_perturber_multistart(seed0_obs, ctx, uranus_inc_raan_at_epoch(EPOCH), space=space,
                                 mode="le_verrier_1846", n_workers=N_WORKERS)  # fmt: skip
    b = r.best
    assert b.dof == 82 and b.theta[6] == pytest.approx(lv.log10_m_earth)
    assert b.chi2 == pytest.approx(330.44, abs=0.05)
    assert _galle_dlon_deg(b.theta, ctx) == pytest.approx(8.04, abs=0.01)


@pytest.mark.slow
@pytest.mark.parametrize(("e_max", "chi2", "dlon"), [(0.1, 72.78, 2.79), (0.05, 72.38, 5.20)])
def test_e_capped_free_a_on_seed0(ctx, seed0_obs, e_max, chi2, dlon):
    """Free a with e capped (24 starts). Measured values pinned; note 0.05 beats 0.1 in χ²,
    which a global optimiser could not do — the multistart misses minima at Δχ² ~ 0.5. ~15 s."""
    from fork.model.forward import uranus_inc_raan_at_epoch

    r = fit_perturber_multistart(seed0_obs, ctx, uranus_inc_raan_at_epoch(EPOCH), e_max=e_max,
                                 n_workers=N_WORKERS)  # fmt: skip
    b = r.best
    assert b.theta[8] <= e_max + 1e-12
    assert b.chi2 == pytest.approx(chi2, abs=0.01)
    assert _galle_dlon_deg(b.theta, ctx) == pytest.approx(dlon, abs=0.01)
