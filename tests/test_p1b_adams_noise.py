"""Phase 1b step 3c (docs/notes/P1B_adams_noise.md § 7): the seeds, the noise scales and the set
D window counts, the one-start solve against the multistart solve, the flags, P and the reading
rows on constructed inputs; the slow test regenerates step 3b's synthetic column and checks that
the noise-free draw reproduces D3(i)."""

import math

import numpy as np
import pytest

from fork.export import p1b_adams as pa
from fork.export import p1b_adams_noise as N
from fork.export import p1b_adams_repro as R
from fork.historical import adams_1846 as A


@pytest.fixture(scope="module")
def memoir():
    return A.load_memoir()


def test_seeds_are_the_plan_s_and_disjoint_from_earlier_steps():
    assert N.seed_for("a", 1) == [2026, 11, 3, 8, 0, 1]
    assert N.seed_for("b", 2) == [2026, 11, 3, 8, 1, 2]
    assert N.seed_for("c", 1) == [2026, 11, 3, 8, 2, 1]
    ours = {tuple(N.seed_for(mo, h)) for mo in N.MODELS for h in N.HYPS}
    assert ours.isdisjoint({tuple(j.rng_seed()) for j in pa.all_jobs()})
    assert ours.isdisjoint({tuple(j.rng_seed()) for j in R.d3iii_jobs()})
    assert N.N_DRAWS == 10_000 and N.N_MULTI == 1_000 and N.NU_T == 4


def test_noise_matrix_is_reproducible_and_scaled():
    sig = np.array([1.0, 2.0, 0.5])
    a1 = N.noise_matrix("a", 1, sig, 2000)
    a2 = N.noise_matrix("a", 1, sig, 2000)
    assert a1.shape == (2000, 3) and np.array_equal(a1, a2)
    assert np.std(a1[:, 1]) == pytest.approx(2.0, rel=0.1)
    # (c) is heavier-tailed than (b) on the same scales and a different stream
    b = N.noise_matrix("b", 1, sig, 20000)
    c = N.noise_matrix("c", 1, sig, 20000)
    assert not np.array_equal(b[:5], c[:5])
    assert np.mean(np.abs(c[:, 0]) > 3) > 2 * np.mean(np.abs(b[:, 0]) > 3)


def test_model_a_sigmas_are_his_rms_times_the_factor(memoir):
    assert A.SCALE_A == pytest.approx(math.sqrt(29 / 21), abs=1e-12)
    for hyp in (1, 2):
        rms = A.residual_rms(memoir, hyp)
        sig = A.noise_sigmas(memoir, "a", hyp=hyp)
        assert len(sig) == 29 and 1690 not in sig
        assert sig[1783] == pytest.approx(rms["modern"] * A.SCALE_A)
        assert sig[1712] == pytest.approx(rms["ancient"] * A.SCALE_A)
    assert A.residual_rms(memoir, 1)["modern"] == pytest.approx(1.24, abs=0.01)
    assert A.residual_rms(memoir, 1)["ancient"] == pytest.approx(6.21, abs=0.01)


def test_model_b_sigmas_follow_the_counts_with_the_floor(memoir):
    grp = {"ancient": 5.48, "greenwich": 2.22, "paris_1835_1845": 2.18}
    counts = {1783: 16, 1712: 4, 1840: 9}
    sig = A.noise_sigmas(memoir, "b", counts, 1, grp)
    assert sig[1783] == pytest.approx(2.22 / 4)
    assert sig[1712] == pytest.approx(5.48 / 2)
    assert sig[1840] == pytest.approx(2.18 / 3)
    assert sig[1780] == pytest.approx(2.22)  # empty window: floored at one Greenwich measurement
    assert sig[1813] == pytest.approx(2.22)
    assert sig[1769] == pytest.approx(5.48)
    assert A.noise_sigmas(memoir, "c", counts, 1, grp) == sig
    with pytest.raises(ValueError):
        A.noise_sigmas(memoir, "b")


def test_set_d_window_counts_place_every_observation_once(memoir):
    c = N.set_d_window_counts(memoir)
    assert c["n_observations"] == pa.EXPECTED_FIT["total"][0] == 145
    assert c["n_assigned"] == c["n_observations"]
    places = c["places"]
    assert len(places) == 29
    assert sorted(y for y, v in places.items() if v["n"] == 0) == list(N.PLAN_EMPTY_WINDOWS)
    assert all(v["groups"] in ([], [v["group_used"] + ":ra"]) for v in places.values())
    assert sum(v["n"] for y, v in places.items() if v["group_used"] == "ancient") == 18
    assert sum(v["n"] for y, v in places.items() if v["group_used"] == "greenwich") == 87
    assert sum(v["n"] for y, v in places.items() if v["group_used"] == "paris_1835_1845") == 40
    # the plan's month-rounded list differs at six Greenwich windows; totals agree
    diff = {y: v["n"] for y, v in places.items() if v["n"] != N.PLAN_COUNTS.get(y, 0)}
    assert set(diff) == {1783, 1786, 1789, 1792, 1795, 1798}
    assert sum(N.PLAN_COUNTS.values()) == 145


def test_warm_start_reproduces_the_multistart_minimum_on_his_data(memoir):
    s = A.setup_from_print(memoir, 1)
    rows = A.place_rows(memoir, s)
    w = A.weights_for(memoir, 1, "c", rows)
    full = A.solve_least_squares(rows, w, s, A.printed_solution(memoir, 1))
    x0 = np.array([full["x"][k] for k in A.PARAMS])
    warm = A.solve_warm(rows, w, s, x0)
    assert warm["converged"] and not warm["bound"]
    assert abs(warm["chi2"] - full["chi2"]) < 1e-6
    assert np.max(np.abs(warm["x"] - x0)) < 1e-3
    assert np.all(np.isfinite(warm["sigma_unscaled"]))
    # the formal σ with the χ²/dof factor is the unscaled one times √(χ²/dof)
    f = math.sqrt(full["chi2"] / full["dof"])
    assert warm["sigma_unscaled"][4] * f == pytest.approx(full["sigma"]["theta"], rel=1e-3)
    multi = A.solve_multistart(rows, w, s, x0)
    assert len(A.theta_grid_starts(x0)) == 12
    assert multi["converged"] and multi["chi2"] <= full["chi2"] + 1e-6
    assert N._wrap(multi["x"][4] - x0[4]) == pytest.approx(0.0, abs=0.01)


def test_flags_on_a_bound_start(memoir):
    s = A.setup_from_print(memoir, 1)
    rows = A.place_rows(memoir, s)
    w = A.weights_for(memoir, 1, "c", rows)
    r = A.solve_warm(rows, w, s, np.array([0, 0, 0, 0, 310.0, 0.8, 2.0, 310.0]), max_nfev=2)
    assert r["converged"] is False  # stopped on max_nfev, not on a tolerance
    r = A.solve_warm(rows, w, s, np.array([0, 0, 0, 0, 310.0, -1.0, 25.0, 310.0]))
    assert math.isfinite(r["chi2"])  # the start is clipped into the bounds, not rejected


def test_p_value_and_the_upper_bound():
    d = np.array([0.0, 5.0, 10.0, 12.0, 15.0])
    p = N.p_value(d, 11.0)
    assert p["n_exceeding"] == 2 and p["n"] == 5 and p["P"] == pytest.approx(0.4)
    assert "P_upper_bound" not in p
    p = N.p_value(d, 11.0, n_failed=5)
    assert p["P_upper_bound"] == pytest.approx(0.7)
    assert N.p_value(d, 20.0)["P"] == 0.0 and N.p_value(d, 20.0)["se"] == pytest.approx(0.2)


def test_reading_rows_on_constructed_p():
    assert N.reading(0.10, 0.20)["row"] == 1
    assert N.reading(0.05, 0.05)["row"] == 1
    assert N.reading(0.0, 0.002)["row"] == 2
    assert N.reading(0.0029, 0.0)["row"] == 2
    assert N.reading(0.10, 0.001)["row"] == 3 and "(a) admits" in N.reading(0.10, 0.001)["text"]
    assert N.reading(0.001, 0.10)["row"] == 3 and "(b) admits" in N.reading(0.001, 0.10)["text"]
    assert N.reading(0.01, 0.02)["row"] == 3
    assert N.reading(0.001, 0.01)["row"] == 3


def test_multistart_fraction_on_constructed_arrays():
    k = len(N._FIELDS) + 2
    warm = np.zeros((4, k))
    multi = np.zeros((3, k))
    i_lon = N._FIELDS.index("true_lon")
    i_chi = N._FIELDS.index("chi2")
    warm[:, -2] = 1.0
    multi[:, -2] = 1.0
    warm[:, i_lon] = [300.0, 300.0, 300.0, 300.0]
    multi[:, i_lon] = [300.5, 302.0, 300.0]
    warm[:, i_chi] = [1.0, 1.0, 1.0, 1.0]
    multi[:, i_chi] = [1.0, 0.5, 1.0]
    f = N.multistart_fraction(warm, multi)
    assert f["n_checked"] == 3 and f["n_differing"] == 1
    assert f["fraction_differing"] == pytest.approx(1 / 3)
    assert f["fraction_multistart_lower_chi2"] == pytest.approx(1 / 3)


@pytest.mark.slow
def test_the_noise_free_draw_reproduces_d3_i(memoir):
    """The regenerated synthetic column through the worker with zero noise lands on step 3b's
    D3(i) solution (O_i = −9.51° by D1(c), −9.92° by his chain) to 1e-6 of itself."""
    syn = R.synthetic_normal_places(memoir)
    col = {int(k): float(v) for k, v in syn["oc_mean_arcsec"].items()}
    T = R.perturber_lon_1846(syn["truth_theta"], syn["ctx"])
    s = A.setup_from_print(memoir, 1)
    rows = A.place_rows(memoir, s, oc=col)
    w = A.weights_for(memoir, 1, "c", rows)
    full = A.solve_least_squares(rows, w, s, A.printed_solution(memoir, 1))
    x0 = np.array([full["x"][k] for k in A.PARAMS])
    eps_c = syn["c_fit"]["eps_1810_ecl_j2000_deg"]
    O_full = N._wrap(
        N.longitudes_of(x0, eps_c, s.n_prime_deg)["true_lon_1846_ecl_j2000"] - T["ecl_j2000"]
    )
    assert O_full == pytest.approx(-9.51, abs=0.02)
    warm = A.solve_warm(rows, w, s, x0)
    O_warm = N._wrap(
        N.longitudes_of(warm["x"], eps_c, s.n_prime_deg)["true_lon_1846_ecl_j2000"] - T["ecl_j2000"]
    )
    assert abs(O_warm - O_full) < 1e-6
    import copy

    ms = copy.deepcopy(memoir)
    for p in ms.places:
        p["oc_mean"] = col[p["year"]]
    ch = N.chain_solve(ms, 1, rows, w, s)
    assert ch["converged"] and ch["n_roots"] == 6 and ch["n_positive"] == 2
    O_ch = N._wrap(
        N.longitudes_of(ch["x"], eps_c, s.n_prime_deg)["true_lon_1846_ecl_j2000"] - T["ecl_j2000"]
    )
    assert O_ch == pytest.approx(-9.92, abs=0.02)


def test_reading_with_flagged_intervals():
    assert N.reading_interval([0.10, 0.10], [0.20, 0.20])["row"] == 1
    assert N.reading_interval([0.001, 0.001], [0.0, 0.002])["row"] == 2
    r = N.reading_interval([0.001, 0.02], [0.0, 0.0])
    assert r["row"] == 3 and "(a)" in r["text"] and "0.003" in r["text"]
    r = N.reading_interval([0.10, 0.10], [0.04, 0.06])
    assert r["row"] == 3 and "(b)" in r["text"] and "0.05" in r["text"]
    assert N.reading_interval([0.01, 0.03], [0.01, 0.02])["row"] == 3


def test_branch_summary_on_constructed_arrays():
    k = len(N._FIELDS) + 2
    i_lon, i_th = N._FIELDS.index("true_lon"), N._FIELDS.index("theta")
    warm, multi = np.zeros((5, k)), np.zeros((5, k))
    warm[:, -2] = multi[:, -2] = 1.0
    warm[:, i_lon] = [100.0, 100.0, 100.0, 100.0, 100.0]
    multi[:, i_lon] = [100.0, 100.5, 115.0, 130.0, 60.0]  # Δ = 0, 0.5, 15, 30, −40 (T = 100)
    multi[:, i_th] = [10.0, 10.0, 10.0, 200.0, 10.0]
    b = N.branch_summary(warm, multi, 100.0, 0.0, 12.0, 195.0, 10.0)
    assert b["n_branch_changes"] == 3 and b["n_exceeding"] == 2
    assert b["n_exceeding_branch_changes"] == 2
    assert b["fraction_of_exceeding_that_are_branch_changes"] == pytest.approx(1.0)
    assert b["P_same_branch_only"] == 0.0
    assert b["exceeding_branch_within_10deg_of_real_theta"] == 1


def test_post_hoc_on_branch_counts_only_the_real_data_branch():
    # four draws on the branch (θ near 310°), two on the second minimum (θ near 95°)
    f = list(N._FIELDS)
    arr = np.zeros((6, len(f) + 2))
    arr[:, -2] = 1.0
    arr[:, f.index("theta")] = [310.0, 315.0, 320.0, 300.0, 95.0, 96.0]
    arr[:, f.index("true_lon")] = [320.0, 332.0, 310.0, 345.0, 130.0, 131.0]
    r = N.post_hoc_on_branch(arr, 329.0, -9.0, {"primary": 11.0}, 310.0)
    assert r["post_hoc"] is True
    assert r["n_on_branch"] == 4
    assert r["fraction_on_branch"] == pytest.approx(4 / 6)
    # Δ = (lon − 329) + 9: 0, 12, −10, 25 on the branch; the far draws (Δ ≈ +170) are left out
    assert r["P"]["primary"]["n_exceeding"] == 2
    assert r["P"]["primary"]["P"] == pytest.approx(0.5)
