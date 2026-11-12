"""Phase 1b step 3b (fork.export.p1b_adams_repro): stage D's pieces that need no run. The
nonlinear model agrees with the chain's arithmetic, the theory rebuilt at his α agrees with his
printed arguments and constants, the least squares cannot do worse than his point, the α grid
and its bands, the D3(iii) job's seeds and space, the synthetic data are noise-free, and the
Q3 table's rows on constructed offsets."""

import math

import numpy as np
import pytest

from fork.export import p1b_adams as pa
from fork.export import p1b_adams_repro as R
from fork.historical import adams_1846 as A


@pytest.fixture(scope="module")
def memoir():
    return A.load_memoir()


@pytest.mark.parametrize("hyp", [1, 2])
def test_the_nonlinear_model_is_the_chain_s_arithmetic(memoir, hyp):
    """§§ 51–52 totals at his printed solution with its § 50 corrections, from both code paths."""
    run = A.run_hypothesis_1 if hyp == 1 else A.run_hypothesis_2
    ch = A.chain_with_print_fallback(memoir, hyp, run(memoir, "C0"))
    corr = A.element_corrections(ch, memoir, ch.solution)
    x = A.x_of(ch.solution, corr)
    s = A.setup_from_print(memoir, hyp)
    totals = {r["year"]: r["total"] for r in A.residuals(memoir, hyp, ch.solution, corr)}
    for year, _kind, t, eq in A.place_rows(memoir, s):
        assert A.model_value(eq, t, x, s) == pytest.approx(totals[year], abs=1e-4)
    sol = ch.solution
    hk, P, Q = A.hkpq_from_elements(sol.theta, sol.m_prime, sol.e_prime_units, sol.beta_prime, s)
    assert P == pytest.approx(ch.solution.P, abs=1e-6)
    assert Q == pytest.approx(ch.solution.Q, abs=1e-6)
    for u in A.HKPQ:
        assert hk[u] == pytest.approx(ch.solution.hkpq[u], abs=1e-6)


@pytest.mark.parametrize("hyp,alpha", [(1, 0.5), (2, 0.515)])
def test_the_theory_rebuilt_at_his_alpha_matches_his_print(memoir, hyp, alpha):
    s = A.setup_from_print(memoir, hyp)
    sa = A.setup_at_alpha(memoir, alpha, hyp)
    for u, rate in s.th.rate.items():
        assert abs(sa.th.rate[u] - rate) * 60 <= 0.25, u  # his arguments are printed to 0.1′
    assert sa.n_prime_deg == pytest.approx(s.n_prime_deg, abs=0.5 / 60)
    for k in ("c_sin1", "c_sin2", "c_p1", "c_p2", "c_p3e"):
        assert abs(sa.th.c[k] - s.th.c[k]) <= 0.09  # B2's achieved level
    # the second-order constants scale with c_p3e: the printed pairs within 0.02″
    assert sa.th.c["c_sin3_e2"] == pytest.approx(memoir.const(hyp, 12 if hyp == 1 else 33,
                                                              "c_sin3_e2"), abs=0.02)  # fmt: skip
    assert sa.th.c["c_sin3_e2e"] == pytest.approx(memoir.const(hyp, 12 if hyp == 1 else 33,
                                                               "c_sin3_e2e"), abs=0.02)  # fmt: skip
    assert A.setup_at_alpha(memoir, alpha, hyp, drop_second_order=True).th.c["c_sin3_e2"] == 0.0


def test_the_arguments_follow_his_rule_for_any_n_prime():
    r = A.rates_at(13.01, 4.6)
    assert r["h1"] == pytest.approx(8.41) and r["h3"] == pytest.approx(25.23)
    assert r["p1"] == pytest.approx(4.6) and r["p2"] == pytest.approx(3.81)
    assert r["p3"] == pytest.approx(12.22) and r["dx2"] == pytest.approx(26.02)


def test_least_squares_does_not_do_worse_than_his_point(memoir):
    """Weighting (c) on Hypothesis I: the minimum χ² is at or below his solution's."""
    hyp = 1
    s = A.setup_from_print(memoir, hyp)
    rows = A.place_rows(memoir, s)
    w = A.weights_for(memoir, hyp, "c", rows)
    ch = A.chain_with_print_fallback(memoir, hyp, A.run_hypothesis_1(memoir, "C0"))
    his = A.chi2_of(rows, w, s, A.x_of(ch.solution, A.element_corrections(ch, memoir, ch.solution)))
    r = A.solve_least_squares(rows, w, s, A.printed_solution(memoir, hyp))
    assert r["chi2"] <= his + 1e-6
    assert r["dof"] == 29 - 8 and r["n_rows"] == 29
    assert r["solution"].m_prime > 0
    assert all(math.isfinite(v) for v in r["sigma"].values())


def test_weights_and_the_printed_rms(memoir):
    rows = A.place_rows(memoir, A.setup_from_print(memoir, 1))
    assert len(rows) == 29 and sum(k == "ancient" for _, k, _, _ in rows) == 8
    wa, wb, wc = (A.weights_for(memoir, 1, sch, rows) for sch in "abc")
    assert np.all(wa == 1.0)
    assert np.sum(wb > 0) == 21 and np.all(wb[[k == "ancient" for _, k, _, _ in rows]] == 0)
    rms = A.residual_rms(memoir, 1)
    assert rms["modern"] < rms["ancient"]
    first_modern = [k for _, k, _, _ in rows].index("modern")
    assert wc[first_modern] == pytest.approx(1 / rms["modern"] ** 2)
    assert wc[0] == pytest.approx(1 / rms["ancient"] ** 2)


def test_alpha_grid_and_its_bands():
    grid = A.alpha_grid()
    assert len(grid) == 21 and grid[0] == 0.40 and grid[-1] == 0.60 and 0.5 in grid
    assert A.POLES["2:1"] == pytest.approx(0.6300, abs=5e-5)
    assert A.POLES["3:1"] == pytest.approx(0.4807, abs=5e-5)
    assert A.alpha_flags(0.5) == [] and A.alpha_flags(0.515) == []
    assert A.alpha_flags(0.48) == ["pole_limited_3_1"]
    assert A.alpha_flags(0.47) == ["pole_limited_3_1"] and A.alpha_flags(0.49) == [
        "pole_limited_3_1"
    ]
    assert A.alpha_flags(0.46) == [] and A.alpha_flags(0.58) == ["truncation_limited"]
    # the pole is real: the i = 3 eccentric constants change sign across it
    e = 0.046679
    assert A.series_constants(0.47, 0.47**1.5, e)["c_p3e"] > 0
    assert A.series_constants(0.49, 0.49**1.5, e)["c_p3e"] < 0


def test_d3iii_jobs_seeds_and_space():
    jobs = R.d3iii_jobs()
    assert [j.rng_seed() for j in jobs] == [[2026, 11, 3, 7, 0, 1], [2026, 11, 3, 7, 0, 2]]
    assert {tuple(j.rng_seed()) for j in jobs}.isdisjoint(
        {tuple(j.rng_seed()) for j in pa.all_jobs()}
    )
    assert R.NLIVE == 1000 and R.HELD_A_AU == 38.4
    assert R.jd_of_year(1846.762) == pytest.approx(2395574.82, abs=0.01)


def test_epoch_conventions():
    assert R.JD_1810 == pytest.approx(2382267.30, abs=0.01)
    # precession from J2000 back to 1846.762 is about −2.14° in longitude
    assert R.lon_j2000_to_of_date(100.0, R.JD_1846) - 100.0 == pytest.approx(-2.14, abs=0.02)


@pytest.mark.slow
def test_synthetic_set_d_is_noise_free_and_the_space_is_step_3_s():
    obs, ctx, _inc, pq, _box, th = R.d3iii_dataset(with_box=False)
    resid = obs.residuals_arcsec(*obs.predict(th, ctx))
    assert np.max(np.abs(resid)) < 1e-6
    assert obs.n_measurements == pa.EXPECTED_FIT["total"][1]
    assert th.shape == (13,) and th[7] == pytest.approx(30.2, abs=0.3)
    from fork.inference.fit_radec import UranusBox

    box = UranusBox(np.zeros(4), np.ones(4))
    sp = R.d3iii_space(0, pq, box)
    assert sp.ndim == 10 and math.exp(sp.held[1]) == pytest.approx(38.4)


def test_q3_table_rows_on_constructed_offsets():
    O3 = -10.5
    assert R.q3_reading(1.0, -9.0, -1.0, O3)["rows"][0]["row"] == 1
    r = R.q3_reading(-1.0, -9.5, -10.0, O3)
    assert r["rows"][0]["row"] == 2 and "his model" in r["rows"][0]["text"]
    r = R.q3_reading(-1.0, -1.5, -10.0, O3)
    assert r["rows"][0]["row"] == 2 and "his layout" in r["rows"][0]["text"]
    r = R.q3_reading(-1.0, -5.0, -10.0, O3)
    assert r["rows"][0]["row"] == 2 and r["shares"] == pytest.approx({"layout": 0.5, "model": 0.4})
    assert R.q3_reading(-9.0, -9.0, -10.0, O3)["rows"][0]["row"] == 3
    r = R.q3_reading(-6.0, -4.0, -10.0, O3)
    assert r["rows"][0]["row"] == 4 and r["shares"]["layout"] == pytest.approx(0.6)
    r = R.q3_reading(-1.0, -1.0, -5.0, O3)
    assert r["rows"][0]["row"] == 6 and r["rows"][0]["share_distance"] == pytest.approx(5 / 10.5)
    assert R.q3_reading(-1.0, -1.0, -15.0, O3)["rows"][0]["row"] == 7
    rows = R.q3_reading(4.0, -1.0, -10.0, O3)["rows"]
    assert [x["row"] for x in rows] == [3, 5] or [x["row"] for x in rows] == [5]
    assert any(x["row"] == 5 and x["bias_deg"] == 4.0 for x in rows)


def test_q1_reading_rows():
    def table(t0, t1, t2, tb):
        mk = lambda t, b: {"vs_print": {"tier": t}, "vs_brookes": {"tier": b}}  # noqa: E731
        return {"H1/C0": mk(t0, 0), "H1/C1": mk(t1, 0), "H1/C2": mk(t2, tb)}

    assert R.q1_reading(table(1, 2, 2, 0), 1)["row"] == 1
    r = R.q1_reading(table(1, 2, 0, 2), 1)
    assert r["row"] == 2 and "Brookes" in r["text"] and "(n)" in r["text"]
    assert R.q1_reading(table(1, 0, 0, 0), 1)["row"] == 2
    assert R.q1_reading(table(0, 0, 0, 0), 1)["row"] == 3
