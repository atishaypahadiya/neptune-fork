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
