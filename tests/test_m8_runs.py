"""The M8 A8b–A8f runner (fork.export.m8_runs): the job list, seeds and prior as pre-registered
in docs/notes/M8.md. No sampling here."""

import numpy as np
import pytest

from fork.export import m6c
from fork.export import m8_runs as mr


def test_jobs_and_seeds_are_the_pre_registered_ones():
    jobs = mr.all_jobs()
    assert len(jobs) == 48
    assert len({j.id for j in jobs}) == 48
    j = mr.RunJob.from_id("s2_sigma_data/c0/s3")
    assert j.rng_seed() == [2026, 8, 1846, 2, 1, 3]
    assert j.nlive == 2000 and mr.RunJob("headline", 3, 1).nlive == 1000
    assert mr.SEEDS == (1, 2, 3, 4) and mr.SECOND_ORDER is False
    assert mr.STEP_DAYS == pytest.approx(60.0) and mr.QUADRATURE == "cubic"
    # 3 data sets × 4 cells × (4 seeds + the merge)
    assert len(mr.posteriors()) == 60


def test_frozen_inputs_and_code_check_passes():
    for v in mr.VARIANTS:
        mr.check_frozen(v)
    assert m6c.sha256(mr.U_REF_FILE) == mr.U_REF_SHA256


def test_phi_prior_is_the_m6c_cells_restricted_to_phi():
    from fork.export import m8
    from fork.inference.nested import cell_space

    space = m8.problem("headline")[2]
    rng = np.random.default_rng(0)
    for k, cell in enumerate(m6c.variant_cells("headline")):
        cs = cell_space(space, cell)
        mr._S["space"] = cs
        for _ in range(20):
            u = rng.random(11)
            q = cs.prior_transform(u)
            np.testing.assert_allclose(mr._phi_transform(u[4:]), q[4:], rtol=0, atol=1e-15)
            th = cs.theta(q)
            (a_lo, a_hi), e_rng = cell
            assert a_lo <= th[7] <= a_hi
            if e_rng is not None:
                assert e_rng[0] - 1e-12 <= th[8] <= e_rng[1] + 1e-12
        assert [i - 4 for i in cs.periodic_unit_dims] == [3, 4, 6], k


def test_smoke_runs_cannot_write_to_the_run_cache():
    with pytest.raises(SystemExit, match="scratch"):
        mr.do_job(mr.RunJob("headline", 2, 1), maxiter=10)


def test_fit_ess_separates_allocation_cap_from_reweighting():
    """Gate 1 as amended: with Δ ≡ 0 the reweighting share is 1 whatever the cell shares,
    while ESS/N as written equals the allocation cap 1 / (K Σ p²)."""
    from fork.inference.reweight import combine, reweight_cell

    n, lz = 600, [-5.0, -5.0, -1.0, -4.0]
    cells = [reweight_cell(np.zeros(n), z, 0.1) for z in lz]
    res = combine(cells, [np.zeros((n, 1))] * 4, [0.0] * 4, 1.8)
    p = np.exp(np.array(lz) - np.logaddexp.reduce(lz))
    as_written, at0, share = mr.fit_ess(res)
    assert as_written == pytest.approx(1 / (4 * np.sum(p**2)), rel=1e-9)
    assert at0 == pytest.approx(as_written, rel=1e-9)
    assert share == pytest.approx(1.0, rel=1e-9)
    # a spread of Δ in one cell lowers only the share
    rng = np.random.default_rng(1)
    cells[2] = reweight_cell(rng.normal(0.0, 1.0, n), lz[2], 0.1)
    res = combine(cells, [np.zeros((n, 1))] * 4, [0.0] * 4, 1.8)
    _, _, share = mr.fit_ess(res)
    assert 0.2 < share < 0.6


def test_follow_up_campaigns_are_the_pre_registered_ones():
    """docs/notes/M8.md, "Follow-up pre-registration": the A8f rerun and M8b; main unchanged."""
    a8f = mr.all_jobs("a8f")
    assert len(a8f) == 12
    assert {(j.variant, j.cell) for j in a8f} == {("s1_e05", 1), ("s1_e05", 3),
                                                  ("s2_sigma_data", 2)}  # fmt: skip
    assert all(j.nlive == 4000 and j.slices == 20 for j in a8f)
    assert mr.RunJob("s1_e05", 1, 4, "a8f").rng_seed() == [2026, 8, 1846, 1, 2, 4, 2]
    m8b = mr.all_jobs("m8b")
    assert len(m8b) == 16 and {j.variant for j in m8b} == {"headline"}
    assert {j.cell: j.nlive for j in m8b} == {0: 8000, 1: 4000, 2: 4000, 3: 4000}
    assert all(j.slices == 28 for j in m8b)
    assert mr.RunJob("headline", 0, 2, "m8b").rng_seed() == [2026, 8, 1846, 0, 1, 2, 3]
    main = mr.all_jobs()
    assert len(main) == 48 and all(j.slices is None for j in main)
    assert len(mr.posteriors("a8f")) == 15 and len(mr.posteriors("m8b")) == 20
    caches = {c.cache for c in mr.CAMPAIGNS.values()}
    assert len(caches) == 3 and mr.CAMPAIGNS["main"].cache == mr.CACHE


def test_arbiter_is_the_pre_registered_one():
    from fork.export import m8_arbiter as arb

    assert (arb.VARIANT, arb.CELL, arb.NLIVE, arb.SLICES) == ("s1_e05", 0, 4000, 28)
    assert arb.rng_seed(1) == [2026, 6, 1846, 1, 1, 3, 1]
    assert arb.A8B_TOL == {"median": 0.35, "width_rel": 0.10, "ln_z": 0.5}
    # P(a < 25 AU) is cell 0's share: equal ln Z per prior mass gives the prior's share
    p = arb.implied_p_below_25(-360.0, [-360.0, -360.0, -360.0])
    assert 0.0 < p < 1.0
