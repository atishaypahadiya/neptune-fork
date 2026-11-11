"""M10 (fork.export.m10): the pre-registered settings, the Student-t objective Q against the
likelihood module and against M9's Gaussian objective in the ν → ∞ limit, the loss's cost, the
six data sets' sizes, the seeding and the amendment-1 tightness rule. No fitting here."""

import math

import numpy as np
import pytest
from scipy.optimize import least_squares

from fork.export import m9, m10
from fork.model.student_t import _log_t_const, t_offset_modes


def test_settings_are_the_pre_registered_ones():
    assert m10.NU == 4.0 and m10.SEED_NAMESPACE == (2026, 10)
    assert m10.DATASETS == ("all", "post1781", "no_flamsteed", "no_bradley", "no_mayer",
                            "no_lemonnier")  # fmt: skip
    assert m10.EXPECTED_MEASUREMENTS == {"all": 529, "post1781": 492, "no_flamsteed": 519,
                                         "no_bradley": 528, "no_mayer": 527,
                                         "no_lemonnier": 505}  # fmt: skip
    assert (m10.MAX_ROUNDS, m10.MARGINAL_SPREAD_TOL, m10.BORDERLINE_AU) == (4, 1.0, 0.2)
    # what M10 takes from M9 unchanged
    assert m10.POINTS_AU == m9.POINTS_AU and m10.CONT_TOL_CHI2 == 0.01
    assert (m10.DLNL_THRESHOLD, m10.X_AU, m10.FLAT_RANGE_AU, m10.TIGHT_AU) == (-2, 4, (26, 44), 8)
    assert m10.FROZEN["docs/notes/M9_profile.json"].startswith("adbb9317b0e7")
    m10.check_frozen()
    s = m10.t_scales()
    assert s.shape == (8,) and abs(s[0] - 4.796) < 1e-3 and abs(s[-1] - 1.076) < 1e-3


def _toy(rng, n=60, ng=3):
    g = rng.integers(0, ng, n)
    r = rng.standard_t(4, n) * 2.0 + np.array([0.5, -1.0, 2.0])[g]
    r[:3] += 30.0  # outliers
    scale = np.array([1.5, 2.0, 2.5])
    tau = np.array([2.0, 2.0, 5.0])
    return r, g, scale, tau


def test_q_is_the_t_data_term_plus_the_prior_at_the_irls_offsets():
    r, g, scale, tau = _toy(np.random.default_rng(1))
    f = m10.t_residuals(r, g, scale, tau)
    assert f.shape == (r.size + scale.size,)
    b = t_offset_modes(r, g, scale, tau, m10.NU)
    data = np.sum((m10.NU + 1) * np.log1p((r - b[g]) ** 2 / (m10.NU * scale[g] ** 2)))
    prior = np.sum((b / tau) ** 2)
    assert m10.q_of(f) == pytest.approx(data + prior, rel=1e-12)
    # the normalised marginal differs from Q by a bounded amount (the Laplace-like term)
    from fork.model.student_t import t_marginal_m2lnl

    marg = t_marginal_m2lnl(r, g, scale, tau, m10.NU)
    const = -2 * r.size * _log_t_const(m10.NU) + 2 * np.sum(np.log(scale[g]))
    assert abs((marg - const) - m10.q_of(f)) < 10.0
    w = m10.irls_weights(r, g, scale, tau)
    assert w.shape == (r.size,) and np.all(w <= (m10.NU + 1) / m10.NU + 1e-12)
    assert np.all(w[:3] < 0.1)  # the planted outliers are discounted


def test_residuals_stay_finite_for_absurd_offsets():
    """A start far from the data (offsets of thousands of arcsec) must give finite residuals,
    or scipy refuses the initial point (the first launch died on this): the pseudo-residual's
    exponent is capped, and anything non-finite becomes PENALTY."""
    # 200 coherent residuals of 150″ in one group: the IRLS mode follows them (b ≈ 150″, 75 τ),
    # and the uncapped pseudo-residual would be √(ν expm1(1125)) = ∞
    g = np.zeros(200, dtype=int)
    r = np.full(200, 150.0)
    f = m10.t_residuals(r, g, np.array([2.0]), np.array([2.0]))
    assert np.all(np.isfinite(f)) and abs(f[:200].max()) < 1.0
    assert (m10.NU + 1) * np.log1p(f[200] ** 2 / m10.NU) == pytest.approx(
        (m10.NU + 1) * m10.PSEUDO_CAP, rel=1e-9
    )
    # the cap is far outside anything an optimum reaches: at |b/τ| ≈ 10 the mapping is exact
    g2 = np.array([0, 0, 1, 1])
    r2 = np.array([20.0, 20.0, 50.0, 50.0])
    f2 = m10.t_residuals(r2, g2, np.array([1.0, 1.0]), np.array([2.0, 5.0]))
    b = t_offset_modes(r2, g2, np.array([1.0, 1.0]), np.array([2.0, 5.0]), m10.NU)
    assert (m10.NU + 1) * np.log1p(f2[4:] ** 2 / m10.NU) == pytest.approx((b / [2.0, 5.0]) ** 2)


def test_loss_cost_is_half_q_including_the_pseudo_residuals():
    r, g, scale, tau = _toy(np.random.default_rng(2))
    f = m10.t_residuals(r, g, scale, tau)
    rho = m10.t_loss(f**2)
    assert rho.shape == (3, f.size)
    assert 0.5 * rho[0].sum() == pytest.approx(0.5 * m10.q_of(f), rel=1e-12)
    # derivatives of ρ by finite differences
    z = np.linspace(0.0, 50.0, 7)
    h = 1e-5
    d1 = (m10.t_loss(z + h)[0] - m10.t_loss(z - h)[0]) / (2 * h)
    assert np.allclose(d1, m10.t_loss(z)[1], rtol=1e-6, atol=1e-8)
    d2 = (m10.t_loss(z + h)[1] - m10.t_loss(z - h)[1]) / (2 * h)
    assert np.allclose(d2, m10.t_loss(z)[2], rtol=1e-5, atol=1e-8)
    # scipy's cost with the loss is ½ Q at a fixed residual
    sol = least_squares(lambda x: f + 0.0 * x.sum(), np.zeros(2), method="trf", loss=m10.t_loss,
                        max_nfev=2)  # fmt: skip
    assert sol.cost == pytest.approx(0.5 * m10.q_of(f), rel=1e-12)


def test_nu_to_infinity_with_sigma_scales_reproduces_m9s_augmented_chi2():
    """At the true 1781 state (Uranus and Neptune from Horizons), Q with ν → ∞ and s_g = σ_g
    equals M9's exact objective, the headline's augmented χ² (‖whitened‖²)."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.to_json import A1_EPOCH_JD_TDB as E
    from fork.model import bodies
    from fork.model.forward import (
        ForwardContext,
        perturber_theta_from_state,
        uranus_theta_from_state,
    )

    obs = m10.observations("all")
    ctx = ForwardContext.from_horizons(E)
    theta = np.concatenate([uranus_theta_from_state(about("uranus", [E])[0]),
                            perturber_theta_from_state(about("neptune", [E])[0],
                                                       bodies.NEPTUNE.mass_msun)])  # fmt: skip
    r = obs.residuals_arcsec(*obs.predict(theta, ctx))
    chi2 = float(np.sum(obs.whitened(theta, ctx) ** 2))
    nu = 1e9
    f = m10.t_residuals(r, obs.group_index, obs.sigma_arcsec, obs.offset_prior_arcsec, nu=nu)
    assert m10.q_of(f, nu=nu) == pytest.approx(chi2, rel=1e-6)
    assert chi2 > 500.0


def test_six_data_sets_have_the_planned_measurement_counts():
    for ds, n in m10.EXPECTED_MEASUREMENTS.items():
        obs = m10.observations(ds)
        assert obs.n_measurements == n, ds
        assert len(m10.row_ids(ds)) == n
    anc = m10.observations("all").group_counts()[:2].sum()
    assert anc == 37
    # the four observers' rows add up to the ancient group
    assert sum(529 - m10.EXPECTED_MEASUREMENTS[d] for d in m10.OBSERVER) == 37


def test_random_starts_are_deterministic_and_differ_from_m9s():
    a = m10.random_starts(0, 5)
    b = m10.random_starts(0, 5)
    assert a == b and len(a) == 32
    assert all(0.0 <= s["log10_m"] <= 2.0 and s["k"] ** 2 + s["h"] ** 2 <= 0.3 for s in a)
    assert all(0.0 <= s["inc"] <= math.radians(10.0) and len(s["u"]) == 4 for s in a)
    g = m9.random_starts(0, 5)
    assert a[0]["log10_m"] != g[0]["log10_m"] and a[0]["lam0"] != g[0]["lam0"]
    assert m10.random_starts(0, 5, rerun=True)[0]["origin"] == "rerun:0"
    assert m10.random_starts(1, 5)[0]["log10_m"] != a[0]["log10_m"]


def test_posterior_tightness_rule_of_amendment_1():
    assert m10.posterior_tightness(7.5) == "tight"
    assert m10.posterior_tightness(7.79) == "tight"
    assert m10.posterior_tightness(7.8) == "borderline"
    assert m10.posterior_tightness(7.99) == "borderline"
    assert m10.posterior_tightness(8.2) == "borderline"
    assert m10.posterior_tightness(8.21) == "not tight"


def test_interpretation_rows():
    def fake(reading):
        return None if reading is None else {"reading": {"reading": reading}}

    m9_out = {"datasets": {"all": fake("constrained")}}
    data = {"all": fake("constrained"), "post1781": fake("flat"),
            "no_flamsteed": fake("constrained"), "no_bradley": fake("constrained"),
            "no_mayer": fake("constrained"), "no_lemonnier": fake("constrained")}  # fmt: skip
    rows = [r["row"] for r in m10.interpret(data, m9_out)["rows"]]
    assert rows == ["i", "iv"]
    data["no_lemonnier"] = fake("flat")
    out = m10.interpret(data, m9_out)
    assert [r["row"] for r in out["rows"]] == ["i", "iii"]
    assert out["rows"][1]["observers"] == ["Le Monnier"]
    data["all"] = fake("flat")
    out = m10.interpret(data, m9_out)
    assert [r["row"] for r in out["rows"]] == ["ii", None]
    data["all"] = fake("partial")
    assert [r["row"] for r in m10.interpret(data, m9_out)["rows"]] == [None, None]


def test_comparison_table_reads_the_frozen_m9_results():
    import json

    m9_out = json.loads(m10.M9_JSON.read_text())
    t = m10.comparison_table({}, m9_out)
    g = t["gaussian_m9"]
    assert g["all"]["reading"] == "constrained" and g["post1781"]["reading"] == "flat"
    assert g["all"]["a_max_au"] == 30.0 and abs(g["all"]["dlnl_le_verrier"] + 3.15) < 0.01
    assert abs(g["all"]["dlnl_26"] + 17.96) < 0.01 and g["post1781"]["dlnl_minus_2"] == [None, None]
    assert all(v is None for v in t["student_t_m10"].values())


def test_site_payload_is_compact_and_matches_the_report():
    import json

    out = json.loads(m10.OUT_JSON.read_text())
    site = m10.site_payload(out)
    assert set(site["datasets"]) == set(m10.DATASETS)
    assert site["leave_one_out"] == ["no_flamsteed", "no_bradley", "no_mayer", "no_lemonnier"]
    for name, d in site["datasets"].items():
        full = out["datasets"][name]
        assert len(d["points"]) == len(full["points"])
        assert d["reading"]["reading"] == full["reading"]["reading"]
        assert "candidates" not in d["points"][0] and "diagnostics" not in d["points"][0]
        assert d["measurements_removed"] == 529 - d["n_measurements"]
        assert (d["observer"] is None) == (name in ("all", "post1781"))
    assert site["s3_posterior_marginal_a"]["tightness"] in ("tight", "borderline", "not tight")
    assert site["comparison_with_m9"]["gaussian_m9"]["all"]["reading"] == "constrained"
    assert len(json.dumps(site)) < 200_000


def test_committed_site_payload_is_current():
    """site/data/m10_profile.json must equal site_payload(docs/notes/M10_profile.json)."""
    import json

    from fork.export.to_json import SITE_DATA

    out = json.loads(m10.OUT_JSON.read_text())
    committed = json.loads((SITE_DATA / m10.SITE_OUT_NAME).read_text())
    assert committed == json.loads(json.dumps(m10.site_payload(out), default=float))
