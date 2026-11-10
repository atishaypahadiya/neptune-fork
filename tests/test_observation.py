"""M6 observation model: groups, marginal likelihood, synthetic data, and the whole chain
(integrator → apparent place) against Horizons (fast; cached Horizons states)."""

import math

import numpy as np
import pytest

from fork.model import observation as om

SIG = {g: (10.0 if g.startswith("ancient") else 3.0) for g in om.GROUPS}
TAU = {g: (20.0 if g.startswith("ancient") else 5.0) for g in om.GROUPS}


@pytest.fixture(scope="module")
def structure():
    return om.structure_from_rows(om.read_rows(), om.ErrorModel(SIG, TAU, "test"))


def test_structure_follows_the_table_and_reads_no_values(structure):
    assert structure.geometry.jd_tt.size == 279
    assert structure.n_measurements == 529
    assert np.all(np.isnan(structure.value_rad))
    counts = dict(zip(om.GROUPS, structure.group_counts().tolist(), strict=True))
    assert counts == {"ancient:ra": 19, "ancient:dec": 18, "greenwich:ra": 87,
                      "greenwich:dec": 80, "paris_1801_1827:ra": 107, "paris_1801_1827:dec": 95,
                      "paris_1835_1845:ra": 66, "paris_1835_1845:dec": 57}  # fmt: skip
    per_epoch = np.bincount(structure.epoch_index)
    assert np.sum(per_epoch == 1) == 29 and np.sum(per_epoch == 2) == 250  # RA-only rows
    ra_only = structure.epoch_index[np.isin(structure.epoch_index, np.flatnonzero(per_epoch == 1))]
    assert not np.any(structure.is_dec[np.isin(structure.epoch_index, ra_only)])


def test_augmented_chi2_is_the_offset_marginal(structure):
    rng = np.random.default_rng(1)
    r = rng.normal(0.0, 1.0, structure.n_measurements) * 4.0 + 2.0
    chi2 = float(np.sum(structure.whitened_from_residuals(r) ** 2))
    g, s, t = structure.group_index, structure.sigma_arcsec, structure.offset_prior_arcsec
    expect = 0.0
    for k in range(len(om.GROUPS)):
        rk = r[g == k]
        expect += np.sum(rk**2) / s[k] ** 2 - (np.sum(rk) / s[k] ** 2) ** 2 / (
            rk.size / s[k] ** 2 + 1.0 / t[k] ** 2)  # fmt: skip
    assert chi2 == pytest.approx(expect, rel=1e-12)
    # brute force: the minimum over offsets of the augmented sum
    from scipy.optimize import minimize

    def aug(b):
        return float(np.sum((r - b[g]) ** 2 / s[g] ** 2) + np.sum(b**2 / t**2))

    assert minimize(aug, np.zeros(len(om.GROUPS)), tol=1e-12).fun == pytest.approx(chi2, rel=1e-6)


def test_synthetic_values_add_offsets_and_noise_on_the_sky(structure):
    truth = om.true_uranus_ra_dec(structure)
    bare, off0 = om.synthetic_values(structure, truth, None, noise=False)
    assert np.all(off0 == 0.0)
    obs = structure.with_values(bare)
    assert np.abs(obs.residuals_arcsec(*truth)).max() < 1e-9
    a, off_a = om.synthetic_values(structure, truth, np.random.default_rng(3))
    b, off_b = om.synthetic_values(structure, truth, np.random.default_rng(3))
    assert np.array_equal(a, b) and np.array_equal(off_a, off_b)  # same seed, same machine
    r = structure.with_values(a).residuals_arcsec(*truth) - off_a[structure.group_index]
    z = r / structure.sigma_arcsec[structure.group_index]
    assert abs(np.mean(z)) < 0.2 and 0.85 < np.std(z) < 1.15


def test_chain_matches_horizons_with_a_three_dimensional_neptune(structure):
    """Our integrator (J, S, U, N from their 1781 Horizons states, inner planets folded in)
    followed by the apparent-place model reproduces Horizons' Uranus to < 0.05″ in every
    group over 1690–1845."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.integrators import DEFAULT_DT_DAYS
    from fork.integrators.wisdom_holman import SCHEME_DEFAULT, integrate_states
    from fork.model.apparent import apparent_ra_dec
    from fork.model.bodies import JUPITER, M_CENTRAL_MSUN, NEPTUNE, SATURN, URANUS

    truth = om.true_uranus_ra_dec(structure)
    obs = structure.with_values(om.synthetic_values(structure, truth, None, noise=False)[0])
    s0 = np.array([about(b, [E])[0] for b in ("jupiter", "saturn", "uranus", "neptune")])
    masses = np.array([M_CENTRAL_MSUN, JUPITER.mass_msun, SATURN.mass_msun, URANUS.mass_msun,
                       NEPTUNE.mass_msun])  # fmt: skip
    x, v = integrate_states(np.ascontiguousarray(s0[:, :3]), np.ascontiguousarray(s0[:, 3:]),
                            masses, structure.geometry.jd_tt - E, DEFAULT_DT_DAYS,
                            SCHEME_DEFAULT)  # fmt: skip
    r = obs.residuals_arcsec(*apparent_ra_dec(x[:, 2], v[:, 2], structure.geometry))
    assert np.abs(r).max() < 0.05


def test_planar_perturber_misfits_declination_by_arcseconds(structure):
    """Known limitation, pinned: the forward model's perturber is planar (i = 0 in ecl_j2000).
    With Neptune flattened that way and everything else true, RA stays < 0.7″ but Dec misses
    by ~2″ in 1835–1845: Neptune's out-of-plane pull (mutual inclination 1.5°) tilts
    Uranus's orbit, which heliocentric longitudes hardly see (0.017″, M4) but Dec does."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.model import bodies
    from fork.model.forward import (
        ForwardContext,
        perturber_theta_from_state,
        uranus_theta_from_state,
    )

    truth = om.true_uranus_ra_dec(structure)
    obs = structure.with_values(om.synthetic_values(structure, truth, None, noise=False)[0])
    ctx = ForwardContext.from_horizons(E)
    theta = np.concatenate([uranus_theta_from_state(about("uranus", [E])[0]),
                            perturber_theta_from_state(about("neptune", [E])[0],
                                                       bodies.NEPTUNE.mass_msun)])  # fmt: skip
    r = obs.residuals_arcsec(*obs.predict(theta, ctx))
    late_dec = obs.group_index == om.GROUPS.index("paris_1835_1845:dec")
    assert 1.5 < math.sqrt(np.mean(r[late_dec] ** 2)) < 2.5
    assert np.abs(r[~obs.is_dec]).max() < 0.7


def test_uranus_only_fit_from_a_data_only_start(structure):
    """The RA/Dec Uranus-only fit starts from the data alone (docs/RULES.md rule 2) and lands
    near the real orbit; the box built from it contains the true in-plane elements."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import fit_uranus_radec, uranus_box_radec
    from fork.model.forward import (
        ForwardContext,
        uranus_inc_raan_at_epoch,
        uranus_theta_from_state,
    )

    truth = om.true_uranus_ra_dec(structure)
    vals, _ = om.synthetic_values(structure, truth, np.random.default_rng(0))
    obs = structure.with_values(vals)
    ctx = ForwardContext.from_horizons(E)
    inc_raan = uranus_inc_raan_at_epoch(E)
    theta, chi2, dof = fit_uranus_radec(obs, ctx, inc_raan)
    true_u = uranus_theta_from_state(about("uranus", [E])[0])
    assert abs(theta[0] - true_u[0]) < 0.05 and abs(theta[1] - true_u[1]) < 0.01
    assert 1.0 < chi2 / dof < 20.0  # Neptune's signal is in the residuals
    box = uranus_box_radec(obs, ctx, inc_raan)
    assert box.contains(uranus_theta_to_x(true_u)[:4])


def test_inclined_perturber_reproduces_horizons(structure):
    """The M6 forward model with an inclined perturber (i, Ω relative to Uranus's plane) and
    θ = the true 1781 states reproduces Horizons' Uranus to < 0.05″ in every group, where the
    planar perturber misses Dec by ~2″."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.model import bodies
    from fork.model.forward import (
        ForwardContext,
        inclined_perturber_theta_from_state,
        predict_states_ecl_j2000,
        uranus_theta_from_state,
    )

    truth = om.true_uranus_ra_dec(structure)
    obs = structure.with_values(om.synthetic_values(structure, truth, None, noise=False)[0])
    ctx = ForwardContext.from_horizons(E)
    th_u = uranus_theta_from_state(about("uranus", [E])[0])
    nep = about("neptune", [E])[0]
    th_p = inclined_perturber_theta_from_state(nep, bodies.NEPTUNE.mass_msun, (th_u[2], th_u[3]))
    assert math.degrees(th_p[5]) == pytest.approx(1.507, abs=0.01)  # Neptune–Uranus mutual i
    theta = np.concatenate([th_u, th_p])
    pos, vel = predict_states_ecl_j2000(theta, [E], ctx)
    np.testing.assert_allclose(pos[0, 3], nep[:3], rtol=0, atol=1e-12)
    np.testing.assert_allclose(vel[0, 3], nep[3:], rtol=0, atol=1e-15)
    assert np.abs(obs.residuals_arcsec(*obs.predict(theta, ctx))).max() < 0.05


def test_uranus_frame_reduces_to_the_ecliptic_when_uranus_is_planar():
    """With Uranus at i = 0 the Uranus-plane frame is ecl_j2000, so an inclined perturber with
    i_rel = 0 starts exactly where the planar one does."""
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.model.forward import ForwardContext, predict_states_ecl_j2000

    ctx = ForwardContext.from_horizons(E)
    th_u = np.array([19.14, 0.047, 0.0, 0.0, 3.02, 1.67])
    planar = np.concatenate([th_u, [1.2, 30.0, 0.01, 1.0, 3.0]])
    inclined = np.concatenate([planar, [0.0, 0.7]])
    a, _ = predict_states_ecl_j2000(planar, [E], ctx)
    b, _ = predict_states_ecl_j2000(inclined, [E], ctx)
    np.testing.assert_allclose(a[0, 3], b[0, 3], rtol=0, atol=1e-12)


def test_uranus_box_from_two_starts_takes_the_lower_chi2(structure, monkeypatch):
    """M10.md, pre-run change (a): the box fit is started from the data-only line start and
    from a second in-plane start; the lower χ² centres the box. With the line start made bad
    (a halved, from where LM stops in a far worse solution, as it does on the real rows), the
    true-orbit start wins and the box contains the truth; with no extra start the original rule
    is unchanged, bit for bit."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.inference import fit_radec
    from fork.inference.fit import uranus_theta_to_x
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch, uranus_theta_from_state

    truth = om.true_uranus_ra_dec(structure)
    vals, _ = om.synthetic_values(structure, truth, np.random.default_rng(0))
    obs = structure.with_values(vals)
    ctx = ForwardContext.from_horizons(E)
    inc_raan = uranus_inc_raan_at_epoch(E)
    true_x = uranus_theta_to_x(uranus_theta_from_state(about("uranus", [E])[0]))[:4]

    plain = fit_radec.uranus_box_radec(obs, ctx, inc_raan)
    same, fits = fit_radec.uranus_box_radec_with_fits(obs, ctx, inc_raan)
    assert same == plain and [f["start"] for f in fits] == ["line"] and fits[0]["centres_box"]

    good = fit_radec.initial_uranus_theta_radec

    def bad_start(obs_, epoch, ir):
        th = good(obs_, epoch, ir)
        th[0] *= 0.5
        return th

    monkeypatch.setattr(fit_radec, "initial_uranus_theta_radec", bad_start)
    box, fits = fit_radec.uranus_box_radec_with_fits(obs, ctx, inc_raan, extra_starts=[true_x])
    by = {f["start"]: f for f in fits}
    assert by["line"]["chi2"] > by["start:0"]["chi2"] * 10  # the bad start is stuck elsewhere
    assert by["start:0"]["centres_box"] and not by["line"]["centres_box"]
    assert box.contains(true_x)
    assert np.allclose(box.centre, by["start:0"]["u"])
