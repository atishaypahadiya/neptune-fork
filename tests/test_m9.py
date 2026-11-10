"""M9 (fork.export.m9): the pre-registered settings, the reading rule on synthetic curves, the
start seeding, dataset B's size, the candidate selection and the screen on a planted row. No
fitting here; the projection identity behind the perturbative search is a slow test."""

import math

import numpy as np
import pytest

from fork.export import m9

GRID = list(m9.A_GRID_AU)


def _curve(fn):
    """fn on the grid, relative to its maximum (Δ ln L)."""
    v = [fn(a) for a in m9.A_GRID_AU]
    return [x - max(v) for x in v]


def test_settings_are_the_pre_registered_ones():
    assert m9.A_GRID_AU == (24, 26, 28, 29, 29.5, 30, 30.5, 31, 32, 34, 36, 40, 44)
    assert m9.REFERENCE_AU == {"neptune": 30.07, "le_verrier": 36.15}
    assert len(m9.POINTS_AU) == 15 and m9.N_GRID == 13
    assert (m9.N_RANDOM, m9.N_POSTERIOR, m9.N_NEAREST, m9.N_RERUN) == (32, 8, 200, 40)
    assert (m9.E_MAX, m9.INC_MAX_DEG) == (0.3, 10.0)
    assert (m9.DLNL_THRESHOLD, m9.X_AU, m9.FLAT_RANGE_AU, m9.TIGHT_AU) == (-2.0, 4.0, (26, 44), 8)
    assert (m9.POLISH_WINDOW_CHI2, m9.MAX_POLISH) == (4.0, 3)
    assert (m9.DISTINCT_LON_DEG, m9.DISTINCT_E, m9.CONT_TOL_CHI2) == (2.0, 0.03, 0.01)
    assert m9.SCREEN_SIGMA == 4.0 and m9.U_FLAG_SIGMA == 3.0 and m9.MIN_EFF_SAMPLES == 50
    assert m9.FROZEN["data/historical/leverrier_1846.csv"].startswith("2f546b9979f6")
    assert m9.FROZEN["docs/notes/M8_u_ref.json"].startswith("4bf4537c88ab")
    m9.check_frozen()


def test_reading_rule_constrained():
    d = _curve(lambda a: -0.5 * ((a - 30.0) / 1.2) ** 2)  # −2 at 30 ± 2.4 AU
    r = m9.read_profile(GRID, d)
    assert r["reading"] == "constrained" and r["a_max_au"] == 30.0 and r["failed"] == []
    assert all(
        r["clauses"][k]
        for k in (
            "a_max_not_at_grid_end",
            "falls_below_within_x_left",
            "falls_below_within_x_right",
            "no_second_branch",
        )
    )
    assert not r["clauses"]["flat_26_44"]


def test_reading_rule_flat_even_with_24_au_below():
    d = _curve(lambda a: -0.1 * ((a - 30.0) / 4.0) ** 2)  # −0.2 at 26, −1.2 at 44
    assert m9.read_profile(GRID, d)["reading"] == "flat"
    d[0] = -3.0  # 24 AU is outside the flat range
    r = m9.read_profile(GRID, d)
    assert r["reading"] == "flat" and r["where"]["below_threshold_in_26_44"] == []


def test_reading_rule_partial_each_clause():
    # a_max at a grid end (and not flat: 26 AU is 5.4 below the maximum at 44 AU)
    r = m9.read_profile(GRID, _curve(lambda a: 0.3 * a))
    assert r["reading"] == "partial" and "a_max_not_at_grid_end" in r["failed"]
    # a_max at the end but every point from 26 to 44 AU within 2: flat, by the rule's order
    assert m9.read_profile(GRID, _curve(lambda a: 0.1 * a))["reading"] == "flat"
    # falls below −2 only farther than 4 AU out
    r = m9.read_profile(GRID, _curve(lambda a: -0.5 * ((a - 30.0) / 2.9) ** 2))  # −2.14 at 24, 36
    assert r["reading"] == "partial"
    assert {"falls_below_within_x_left", "falls_below_within_x_right"} <= set(r["failed"])
    assert "no_second_branch" not in r["failed"]
    # a second branch rises above −2
    d = _curve(lambda a: -0.5 * ((a - 30.0) / 1.2) ** 2)
    d[GRID.index(40.0)] = -1.0
    r = m9.read_profile(GRID, d)
    assert r["reading"] == "partial" and r["failed"] == ["no_second_branch", "flat_26_44"]
    assert r["where"]["second_branch_points"] == [40.0]
    # a missing point is not used
    d[GRID.index(40.0)] = None
    r = m9.read_profile(GRID, d)
    assert r["reading"] == "constrained" and r["missing_points"] == [40.0]


def test_crossings_and_parabola_on_a_quadratic():
    d = _curve(lambda a: -0.5 * ((a - 30.0) / 1.2) ** 2)
    lo, hi = m9.crossings(GRID, d, -2.0)
    assert abs(lo - 27.6) < 0.3 and abs(hi - 32.4) < 0.3
    lo, hi = m9.crossings(GRID, d, -0.5)
    assert abs(lo - 28.8) < 0.3 and abs(hi - 31.2) < 0.3
    p = m9.parabola_max(GRID, d)
    assert abs(p["a_au"] - 30.0) < 1e-9 and abs(p["dlnl"]) < 1e-9
    assert m9.crossings([24.0, 26.0, 28.0], [-1.0, -0.5, 0.0], -2.0) == [None, None]
    assert m9.parabola_max(GRID, _curve(lambda a: 0.1 * a)) is None


def test_random_starts_are_deterministic_and_inside_the_prior():
    s1, s2 = m9.random_starts(0, 3), m9.random_starts(0, 3)
    assert s1 == s2 and len(s1) == 32
    assert s1 != m9.random_starts(0, 4) and s1 != m9.random_starts(1, 3)
    assert m9.random_starts(0, 3, rerun=True) != s1 and len(m9.random_starts(0, 3, 40, True)) == 40
    for s in s1 + m9.random_starts(1, 14, 40, rerun=True):
        assert 0.0 <= s["log10_m"] <= 2.0
        assert s["k"] ** 2 + s["h"] ** 2 <= 0.3 + 1e-12
        assert 0.0 <= s["lam0"] < 2 * math.pi and 0.0 <= s["node"] < 2 * math.pi
        assert 0.0 <= s["inc"] <= math.radians(10.0) + 1e-12


def test_posterior_starts_take_the_nearest_200_then_the_8_most_likely():
    rng = np.random.default_rng(1)
    n = 1000
    q = np.zeros((n, 11))
    q[:, m9.Q_LN_A] = np.log(rng.uniform(20.0, 60.0, n))
    q[:, 4] = rng.uniform(0.0, 2.0, n)
    q[:, 9:11] = 0.01 * rng.standard_normal((n, 2))
    logl = rng.standard_normal(n)
    picks = m9.posterior_starts(q, logl, 30.0)
    assert len(picks) == 8 and [p["origin"] for p in picks] == [f"posterior:{j}" for j in range(8)]
    near = set(np.argsort(np.abs(q[:, m9.Q_LN_A] - math.log(30.0)))[:200])
    chosen = [int(np.flatnonzero(logl == p["sample_ln_l"])[0]) for p in picks]
    assert set(chosen) <= near
    assert sorted(logl[list(near)])[-8:] == sorted(p["sample_ln_l"] for p in picks)
    assert picks[0]["sample_ln_l"] >= picks[-1]["sample_ln_l"]
    for p in picks:
        assert "u" in p and len(p["u"]) == 4 and 0 <= p["inc"] <= 0.1


def test_post1781_is_table77_without_the_ancient_group():
    from fork.model.observation import GROUPS

    a = m9.observations("all")
    b = m9.observations("post1781")
    assert (a.n_measurements, a.obs_id.size) == (529, 279)
    assert (b.n_measurements, b.obs_id.size) == (492, 260)
    counts = b.group_counts()
    assert all(counts[i] == 0 for i, g in enumerate(GROUPS) if g.startswith("ancient"))
    assert counts.sum() == 492 and a.n_measurements - b.n_measurements == 37
    # the epochs and values kept are the modern ones, in order, with the same geometry rows
    keep = ~np.isin(a.group_index, [0, 1])
    np.testing.assert_allclose(b.value_rad, a.value_rad[keep])
    np.testing.assert_allclose(b.geometry.jd_tt, np.unique(a.geometry.jd_tt[a.epoch_index[keep]]))
    assert np.all(b.geometry.jd_tt[b.epoch_index] == a.geometry.jd_tt[a.epoch_index][keep])
    # an empty group's best offset is 0, so its prior term vanishes from the augmented vector
    z = b.whitened_from_residuals(np.ones(492))
    assert z.size == 492 + len(GROUPS) and abs(z[492]) < 1e-15 and abs(z[493]) < 1e-15


def test_candidates_are_distinct_and_within_the_window():
    def o(chi2, lon, e, ok=True, origin="x"):
        return {"chi2": chi2, "galle_lon_deg": lon, "e": e, "success": ok, "origin": origin}

    opt = [
        o(10.0, 300.0, 0.1),
        o(9.0, 300.5, 0.11),
        o(11.0, 310.0, 0.1),
        o(12.5, 300.0, 0.2),
        o(13.5, 320.0, 0.1),
        o(8.0, 300.0, 0.1, ok=False),
        o(12.0, 300.2, 0.12),
    ]
    c = m9.candidates(opt)
    assert [x["chi2"] for x in c] == [9.0, 11.0, 12.5]  # best, then distinct within 4, max 3
    assert m9.candidates([o(1.0, 0.0, 0.0, ok=False)]) == []


def test_screen_flags_a_planted_row():
    from fork.model.observation import GROUPS

    gi = GROUPS.index("greenwich:ra")
    rows = [
        {
            "id": f"LV1846-T77-p129-l{i + 1:02d}-ra",
            "quantity": "ra",
            "source_ref": "src",
            "epoch": f"1790 Jan {i + 1}",
            "misprints_applied": "",
            "lv_residual": f"{i:+.1f}",
        }
        for i in range(12)
    ]
    rows[5]["misprints_applied"] = "mean_time: 1 -> 2"
    resid = np.array([0.3, -0.2, 0.1, 0.0, -0.1, 0.2, 0.0, 0.1, -0.3, 0.2, 0.05, -0.05])
    resid[5] += 5.0 * 1.9  # 5σ in a group with σ = 1.9″
    out = m9.screen_table(np.full(12, gi), resid, rows, {"greenwich:ra": 1.9}, {("129", 6)})
    assert len(out) == 1
    t = out[0]
    assert t["id"] == "LV1846-T77-p129-l06-ra" and (t["page"], t["line"]) == ("129", 6)
    assert abs(t["deviation_sigma"] - 5.0) < 0.2 and t["hand_checked"] is True
    assert t["misprint_logged"] == "mean_time: 1 -> 2"
    assert t["le_verrier_o_minus_c_arcsec"] == "+5.0"
    assert t["neighbours_o_minus_c_arcsec"] == ["+4.0", "+6.0"]
    assert (
        m9.screen_table(
            np.full(12, gi),
            resid - resid[5] * (np.arange(12) == 5),
            rows,
            {"greenwich:ra": 1.9},
            set(),
        )
        == []
    )


@pytest.mark.slow
def test_perturbative_residual_is_the_linear_marginal_profile():
    """‖z − Q Qᵀ z‖² equals ``GaussianMarginal.chi2_min`` (χ² minimised over Uranus)."""
    m9._init_worker("all")
    pr = m9._W["pr"]
    for phi in (
        m9.phi_of([math.log10(17.0), 30.0, 0.05, 0.02, 2.0], [0.02, 1.0]),
        m9.phi_of([1.5, 24.0, 0.3, 0.3, 4.0], [0.1, 3.0]),
    ):
        r = m9._pert_residual(phi)
        assert abs(float(r @ r) - pr.evaluate(phi).marginal.chi2_min) < 1e-6 * max(1.0, r @ r)
