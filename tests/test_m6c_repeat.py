"""M6c cell-2 repeatability study (pre-registered in docs/notes/M6.md): jobs, seeds, the cell-2
space, and the decision rule."""

import numpy as np
import pytest

from fork.export import m6c
from fork.export import m6c_repeat as rp


def test_jobs_settings_and_seeds():
    assert [j.id for j in rp.JOBS] == ["r1_1", "r1_2", "r1_3", "r1_4", "r2_1", "r2_2"]
    r1 = [j for j in rp.JOBS if j.group == "r1"]
    r2 = [j for j in rp.JOBS if j.group == "r2"]
    head = m6c.Job("headline", 2)
    assert all(j.nlive == head.nlive == 1000 and j.slices is None for j in r1)
    assert all(j.nlive == 4000 and j.slices == 2 * (3 + 11) == 28 for j in r2)
    assert rp.HEADLINE_SEED == head.rng_seed()
    seeds = [j.rng_seed() for j in rp.JOBS]
    assert seeds == [[*head.rng_seed(), 1, r] for r in range(1, 5)] + [
        [*head.rng_seed(), 2, r] for r in range(1, 3)
    ]
    assert not {tuple(s) for s in seeds} & {tuple(j.rng_seed()) for j in m6c.all_jobs()}
    assert sum(j.workers for j in r2) == 12 and sum(j.workers for j in r1) == 12
    assert all(rp.RepeatJob.from_id(j.id) == j for j in rp.JOBS)


def test_cell2_space_is_the_headline_cell_with_the_preregistered_dims():
    from fork.inference.priors import UranusBox

    box = UranusBox((19.2, 0.0, 0.0, 0.0), (0.01, 0.001, 0.001, 0.001))
    space = m6c.job_space(m6c.Job("headline", rp.HEADLINE_CELL), (0.0, 0.0), box)
    assert space.ndim == rp.CELL2_NDIM == 11
    assert rp.DEFAULT_SLICES == 3 + space.ndim


def test_spread_against_run_errors():
    s = rp.spread([1.0, 2.0, 3.0], [0.5, 0.5, 0.5])
    assert s["std"] == pytest.approx(1.0) and s["rms_run_error"] == pytest.approx(0.5)
    assert s["chi2"] == pytest.approx(8.0) and s["dof"] == 2
    assert s["p"] == pytest.approx(np.exp(-4.0))  # χ²₂ survival: exp(−χ²/2)
    assert s["implementation_error"] == pytest.approx(np.sqrt(0.75))


def test_decision_rule():
    big = [3.4, 5.7, 4.0, 4.5]  # R1 widths far apart against 0.05° errors
    calm = [4.00, 4.05, 3.97, 4.02]  # R1 widths within dynesty's errors
    d = rp.decide(big, [0.05] * 4, [4.6, 4.8])
    assert d["replace_headline_cell2"] and d["outcome"] == "replace"
    # (a) holds, (b) fails: R2's merged width is a cross-check, not a replacement
    d = rp.decide(calm, [0.05] * 4, [4.6, 4.8])
    assert d["r2_consistent"] and not d["r1_spread_exceeds_run_error"]
    assert not d["replace_headline_cell2"] and d["outcome"] == "cross_check"
    # (a) fails: R2 disagree by more than 10 % of their mean, whatever R1 does
    for r1 in (big, calm):
        d = rp.decide(r1, [0.05] * 4, [4.0, 4.5])
        assert not d["r2_consistent"] and d["outcome"] == "systematic_only"
