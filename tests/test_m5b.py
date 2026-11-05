"""M5b job planning and the Le Verrier-shape single-mode check (fast; the runs themselves are
local and take a day, ``python -m fork.export.m5b run``)."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from fork.export import m5b


def test_job_plan_matches_the_m5b_decision():
    """10 seeds × 4 cells at e ≤ 0.3 (cell [20, 25] at nlive 2000), seeds 0–4 × 3 cells at
    e ≤ 0.1, Le Verrier shape unpartitioned on seeds 0–4 (seed 0 twice, for the check); the
    seed-0 e ≤ 0.3 cells 1–3 reuse the validation runs; ids round-trip; one rng per job."""
    jobs = m5b.all_jobs()
    by_run = {r: [j for j in jobs if j.run == r] for r in m5b.RUNS}
    assert len(by_run["free_a_e0.3"]) == 44 and len(by_run["free_a_e0.1"]) == 15
    assert len(by_run["le_verrier_shape"]) == 6
    assert len(m5b.all_jobs(with_lv=False)) == 61  # only seed 0's two check runs
    nf = [j for j in jobs if j.seed is None]
    assert [(j.run, j.cell) for j in nf] == [("free_a_e0.3", k) for k in range(4)]
    assert [j.seed for j in m5b.ordered(jobs)[-4:]] == [None] * 4  # after everything A3b needs
    assert {j.nlive for j in by_run["free_a_e0.3"] if j.cell == 0} == {2000}
    assert {j.nlive for j in jobs if not (j.run == "free_a_e0.3" and j.cell == 0)} == {1000}
    assert [j.cell for j in jobs if j.reuse_source is not None] == [1, 2, 3]
    assert all(m5b.Job.from_id(j.id) == j for j in jobs)
    draws = {j.id: j.rng().random() for j in jobs}
    assert len(set(draws.values())) == len(jobs)
    assert m5b.ordered(jobs)[0].run == "le_verrier_shape"  # the check decides 4 more jobs


def test_noise_free_check_offsets_in_sigma_on_the_side_facing_truth():
    lon = {"median": 326.5, "lo68": 325.0, "hi68": 330.5, "minus": 1.5, "plus": 4.0,
           "truth": 325.873, "median_minus_truth": 0.627, "truth_in_68": True}  # fmt: skip
    one = {"median": 17.0, "lo68": 16.0, "hi68": 19.0, "minus": 1.0, "plus": 2.0, "truth": 17.15,
           "truth_in_68": True}  # fmt: skip
    nf = {"runs": {"free_a_e0.3": {"galle_lon_ecl_of_date_deg": lon, "mass_earth": one,
                                    "a_au": one, "e": one}}}  # fmt: skip
    rec = m5b.noise_free_check(nf)
    assert rec["galle_lon_ecl_of_date_deg"]["offset_sigma"] == pytest.approx(0.627 / 1.5, abs=1e-3)
    assert rec["mass_earth"]["offset_sigma"] == pytest.approx(-0.15 / 2.0, abs=1e-3)
    assert m5b.noise_free_check({"runs": {}}) == {"status": "missing"}


def _fake_run(rng, logz, lam_centre, frac_second=0.0, n=20000):
    q = np.zeros((n, 6))
    q[:, 4] = rng.normal(1.39, 0.02, n)  # log10 m
    lam = rng.normal(lam_centre, 0.01, n)
    k = int(frac_second * n)
    lam[:k] = rng.normal(lam_centre + 1.0, 0.01, k)  # a second mode ~57° away
    q[:, 5] = lam % (2 * math.pi)
    return SimpleNamespace(samples=q, weights=np.full(n, 1.0 / n), logz=logz, logzerr=0.2)


@pytest.mark.parametrize(("logz2", "lam2", "frac", "single"), [
    (-80.1, 6.27, 0.0, True),    # agree; λ₀ near the 0/2π wrap is handled
    (-80.1, 6.27, 0.02, False),  # 2 % of the mass in a second λ₀ mode
    (-81.0, 6.27, 0.0, False),   # ln Z differs by 1.0 > 1.8 · 0.28
    (-80.1, 6.30, 0.0, False),   # λ₀ medians 3σ apart
])  # fmt: skip
def test_le_verrier_single_mode_check(logz2, lam2, frac, single):
    rng = np.random.default_rng(0)
    r1 = _fake_run(rng, -80.0, 6.27)
    r2 = _fake_run(rng, logz2, lam2, frac)
    rec = m5b.lv_single_mode(r1, r2)
    assert rec["single_mode"] is single


def test_summary_rng_accepts_noise_free_and_is_per_dataset():
    from fork.export.m5b import summary_rng

    a = summary_rng(None, "free_a_e0.3").random(3)
    np.testing.assert_array_equal(a, summary_rng(None, "free_a_e0.3").random(3))  # same seed
    assert not np.allclose(a, summary_rng(0, "free_a_e0.3").random(3))


def test_a3b_binomial_and_per_seed_offsets():
    b = m5b.binomial_calibration(7, 10)
    assert b["spec_range_counts"] == [5, 8]
    assert b["p_calibrated_in_spec_range"] == pytest.approx(0.8157, abs=1e-3)
    assert b["p_value_two_sided"] == pytest.approx(1.0, abs=1e-6)  # 7 is the mode
    assert m5b.binomial_calibration(3, 10)["p_value_two_sided"] < 0.05

    def ds(seed, med, hit):
        lon = {"median": med, "lo68": med - 2, "hi68": med + 4, "minus": 2.0, "plus": 4.0,
               "width68": 6.0, "median_minus_truth": med - 326.0, "truth_in_68": hit}  # fmt: skip
        return {"seed": seed, "noise": True, "runs": {"free_a_e0.3": {
            "galle_lon_ecl_of_date_deg": lon}}}  # fmt: skip

    cal = m5b.calibration({"a": ds(0, 327.0, True), "b": ds(1, 322.0, False)}, "free_a_e0.3")
    assert [r["offset_sigma"] for r in cal["per_seed"]] == pytest.approx([0.5, -1.0])
    assert cal["n_truth_in_68"] == 1 and cal["binomial"]["p_true"] == pytest.approx(0.68)


def test_profile_vs_posterior_flags_ends_beyond_two_degrees():
    lon = {"lo68": 318.0, "hi68": 334.0}
    ok = m5b.profile_vs_posterior([[317.4, 333.0]], lon)
    assert not ok["flag"] and ok["delta_lo_deg"] == pytest.approx(-0.6)
    bad = m5b.profile_vs_posterior([[317.4, 320.0], [322.0, 350.0]], lon, [300.0, 350.0])
    assert bad["flag"] and bad["n_profile_pieces"] == 2
    assert bad["profile_end_at_grid_edge"] == [False, True]
    assert bad["delta_hi_deg"] == pytest.approx(16.0)
