"""The low-distance solution's runner against its plan (docs/notes/S1_mode.md)."""

import json
import math

import numpy as np
import pytest

from fork import ledger
from fork.export import m8_runs
from fork.export import s1_mode as S

PLAN = (S.ROOT / "docs" / "notes" / "S1_mode.md").read_text()


def test_inputs_are_the_plans():
    S.check_inputs()
    for rel, digest in S.INPUTS.items():
        assert f"`{rel}` | `{digest}`" in PLAN


def test_region_masses_are_the_plans():
    w = S.prior_masses()
    S.check_masses(w)
    assert sum(w.values()) == pytest.approx(1.0, abs=1e-12)
    for name, (a, _e) in S.REGIONS.items():
        if name in ("c1", "c2", "c3"):
            continue
        assert f"a in [{a[0]:g}" in PLAN or f"a in [{a[0]:g}, {a[1]:g}" in PLAN


def test_settings_are_the_plans():
    for text in ("nsteps 28", "at least 4,000 live points", "8,000\nlive points and nsteps 56",
                 "`dlogz` 0.1, `frac_remain` 0.001", "s = 101 and 102",
                 "rng `[2026, 15, 1846, 1, r, s]`", "**ESS/N ≥ 0.5**", "6,000 equal-weight",
                 "chunks of 250"):  # fmt: skip
        assert text in PLAN, text
    assert (S.NLIVE, S.NSTEPS, S.X_NLIVE, S.X_NSTEPS) == (4000, 28, 8000, 56)
    assert (S.DLOGZ, S.FRAC_REMAIN, S.N_DRAWS, S.CHUNK) == (0.1, 0.001, 6000, 250)
    assert S.SEEDS == tuple(range(1, 9)) and S.X_SEEDS == (101, 102)
    assert (S.ESS_M, S.ESS_ALL, S.TRUNC_GATE) == (0.5, 0.3, 1e-6)
    assert S.rng_list("M", 3) == [2026, 15, 1846, 1, 1, 3]
    assert [S.REGION_R[n] for n in ("L1", "L2", "c1", "c2", "c3")] == [2, 3, 4, 5, 6]


def test_the_runs():
    r = S.runs()
    assert len(r) == 6 * 8 + 2 and len(set(r)) == len(r)
    assert [s for n, s in r if s > 8] == [101, 102] and all(n == "M" for n, s in r if s > 8)
    jobs = S.all_jobs(S.CACHE, [], None, [])
    assert len(jobs) == 50 + 16
    b = [j for j in jobs if j[0].startswith("b/")]
    assert {j[0] for j in b} == {f"b/c{k}/s{s}" for k in range(4) for s in (5, 6, 7, 8)}
    assert all("--campaign" in j[2] and "m8b" in j[2] for j in b)
    assert all(j[1].parent == m8_runs.CAMPAIGNS["m8b"].cache / "jobs" for j in b)


def test_g0_truth_by_an_independent_integral():
    from scipy.stats import norm

    x1 = np.linspace(0.0, 1.0, 400_001)
    m = 0.5 + S.G0_B * ((x1 - 0.5) ** 2 - S.G0_S1**2)
    inner = norm.cdf((1 - m) / S.G0_S2) - norm.cdf(-m / S.G0_S2)
    z = np.trapezoid(norm.pdf(x1, 0.5, S.G0_S1) * inner, x1)
    z *= (norm.cdf(0.5 / S.G0_SG) - norm.cdf(-0.5 / S.G0_SG)) ** 5
    assert S.g0_true_ln_z() == pytest.approx(math.log(z), abs=1e-6)
    # the density is normalised in R^7: at the centre of the ridge it is the Gaussian's peak
    x = np.full(7, 0.5)
    x[1] = 0.5 - S.G0_B * S.G0_S1**2
    peak = -math.log(S.G0_S1 * S.G0_S2 * S.G0_SG**5) - 3.5 * math.log(2 * math.pi)
    assert S.g0_logl(x) == pytest.approx(peak, abs=1e-12)


def test_g0_verdict():
    assert S.g0_verdict([0.1, -0.2, 0.25], 0.0)["pass"]
    assert not S.g0_verdict([0.1, -0.31], 0.0)["pass"]
    assert not S.g0_verdict([0.29, -0.29, 0.0], 0.3)["pass"]  # -0.59 from the truth
    assert not S.g0_verdict([0.28, -0.28], 0.0)["pass"]  # range 0.56 > 0.5
    assert S.g0_verdict([0.28, -0.28], 0.0)["range"] == pytest.approx(0.56)


def test_words_and_between():
    assert S.verdict_word(0.01, 0.001) == "not at all"
    assert S.verdict_word(0.2, 0.01) == "moved partly"
    assert S.verdict_word(0.7, 0.01) == "moved"
    assert S.verdict_word(0.045, 0.004) == "between not at all and moved partly"
    assert S.verdict_word(0.49, 0.01) == "between moved partly and moved"


def test_p_below_25_and_its_error():
    w = {"M": 0.057542, "L1": 0.107211, "L2": 0.038361, "c1": 0.103156, "c2": 0.166495,
         "c3": 0.527235}  # fmt: skip
    lnz = {"M": -362.0, "L1": -370.0, "L2": -369.0, "c1": -365.0, "c2": -360.0, "c3": -365.1}
    err = {n: 0.1 for n in w}
    p = S.p_below_25(lnz, err, w)
    z = {n: w[n] * math.exp(lnz[n] + 360) for n in w}
    assert p["p"] == pytest.approx((z["M"] + z["L1"] + z["L2"]) / sum(z.values()), rel=1e-12)
    rng = np.random.default_rng(0)
    draws = []
    for _ in range(40_000):
        lz = {n: lnz[n] + rng.normal(0, err[n]) for n in w}
        draws.append(S.p_below_25(lz, {n: 0.0 for n in w}, w)["p"])
    assert p["p_err"] == pytest.approx(np.std(draws), rel=0.05)
    from fork.export import m8_arbiter

    with_cells = m8_arbiter.implied_p_below_25(p["ln_z0"], [lnz["c1"], lnz["c2"], lnz["c3"]])
    assert with_cells == pytest.approx(p["p"], rel=1e-5)  # the plan's masses vs the space's


def test_merge_is_pooled_and_its_ln_z_the_mean_of_z():
    rng = np.random.default_rng(1)
    ds = [{"phi": rng.normal(size=(100, 7)), "lnl_fast": rng.normal(size=100),
           "lnl_b": rng.normal(size=100), "b_trunc_m2lnl": np.zeros(100),
           "galle_deg": rng.uniform(320, 330, 100), "ln_z_fast": -360.0 + i}
          for i in range(4)]  # fmt: skip
    m = S.pooled(ds, 25, np.random.default_rng(2))
    assert len(m["lnl_fast"]) == 100
    assert m["ln_z_fast"] == pytest.approx(math.log(np.mean(np.exp(np.arange(4.0)))) - 360.0)
    assert m["ln_z_fast_err"] == pytest.approx(np.std(np.arange(4.0), ddof=1) / 2)
    assert S.seed_error([359.9, 0.1, 0.0, 359.95], circular=True) < 0.1


def _b(holds, lo=320.6, hi=330.1):
    return {"holds": holds, "merged_95": [lo, hi], "sampler_error_95": {"lo95": 0.15, "hi95": 0.18},
            "conditions": {"1": True, "4p_std_over_sqrt_n_le_0p2": holds}}  # fmt: skip


def _a(resolved, p=None, w=None):
    r = {"resolved": resolved}
    if resolved:
        r.update(p=p, word=w)
    return {
        "reading": r,
        "regions": {"M": {"merged": {"galle_median": 315.2, "a_au": [23.4, 23.9, 24.5]}}},
    }


def test_the_table_of_section_6():
    t = S.range_reading(_b(True), _a(True, 0.01, "not at all"))
    assert t["final"] and t["label"] is None and t["range"] == [320.6, 330.1]
    t = S.range_reading(_b(True), _a(True, 0.12, "moved partly"))
    assert t["final"] and "e <= 0.3" in t["scope"] and t["beside"]["a_au"] == 23.9
    t = S.range_reading(_b(True), _a(False))
    assert not t["final"] and t["range"] == [320.6, 330.1] and t["label"] == "provisional"
    t = S.range_reading(_b(False), _a(True, 0.01, "not at all"))
    assert not t["final"] and t["range"] == list(S.M6C_95) and "4p" in t["why"]
    t = S.range_reading(_b(False), _a(False))
    assert t["part_a"] == "unresolved"
    t = S.range_reading(_b(True), _a(True, 0.04, "between not at all and moved partly"))
    assert t["final"] and "beside" in t  # "between" is not "under 5 %"


def test_part_b_reproduces_m8b_with_its_four_seeds():
    """Part B's summary code on seeds 1 to 4 and M8b's own merge gives site/data/m8b.json."""
    if not (S.b_cache() / "reweight" / "headline__c0__merged" / "draws.npz").exists():
        pytest.skip("M8b's cache is not here")
    b = S.part_b(
        S.CACHE,
        seeds=(1, 2, 3, 4),
        merged=True,
        load=lambda k, s: m8_runs.load_reweighted("headline", k, s, S.b_cache()),
    )
    m = json.loads((S.ROOT / "site" / "data" / "m8b.json").read_text())
    for key in ("lo95", "hi95"):
        assert b["std_95_bounds"][key] == pytest.approx(m["rule_95"]["std_95_bounds"][key],
                                                        abs=1e-9)  # fmt: skip
    for n in ("s1", "s2", "s3", "s4", "merged"):
        mine, pub = b["fits"][n], m["fits"][f"headline/{n}"]
        assert mine["ln_z"] == pytest.approx(pub["ln_z"], abs=1e-9)
        for key in ("lo95", "lo68", "median", "hi68", "hi95"):
            assert mine["galle_lon_ecl_of_date_deg"][key] == pytest.approx(
                pub["galle_lon_ecl_of_date_deg"][key], abs=1e-9
            )
    for r, q in zip(b["a8f"], m["tests"]["a8f"], strict=True):
        assert r["max_rel_dev"] == pytest.approx(q["max_rel_dev"], abs=1e-12)


def test_smoke_options_are_refused_in_the_real_cache():
    for argv in (["job-a", "M", "1", "--maxiter", "5"], ["job-a", "M", "1", "--synthetic", "1"],
                 ["pipeline", "--max-draws", "10"]):  # fmt: skip
        with pytest.raises(SystemExit):
            S.main(argv)
    with pytest.raises(SystemExit):
        S.do_job_a("L1", 101, S.CACHE)  # X's seeds are M's only


def test_pipeline_needs_g0(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "CACHE", tmp_path)
    with pytest.raises(SystemExit, match="G0"):
        S.pipeline(1, tmp_path)


@pytest.mark.slow
def test_a_capped_synthetic_run_writes_its_result_and_row(tmp_path):
    S.do_job_a("M", 1, tmp_path, maxiter=30, synthetic=1)
    r = np.load(tmp_path / "a" / "M" / "s1" / "result.npz")
    assert bool(r["capped"]) and int(r["nlive"]) == 4000 and r["points"].shape[1] == 7
    row = ledger.read()[-1]
    assert row["job"] == "s1_mode/a/M/s1" and row["outcome"] == "done"
    assert row["capped"] and row["synthetic"] and row["scratch"]
    assert row["n_likelihood_calls"] >= 4000
    S.do_job_a("M", 1, tmp_path, maxiter=30, synthetic=1)  # cached: no second row
    assert ledger.read()[-1]["recorded"] == row["recorded"] and len(ledger.read()) == 1
