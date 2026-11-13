"""S1 cell 0 by sub-cells (docs/notes/S1_cell0.md): the cut, the seeds, the spaces and the
reading of § 5 on made-up inputs. No sampling here."""

import math
from pathlib import Path

import numpy as np
import pytest

from fork.export import m6c
from fork.export import s1_cell0 as S


@pytest.fixture(scope="module")
def space():
    return S.problem()[2]


def test_the_subcells_partition_cell_0(space):
    from fork.inference.nested import cell_space

    w = S.prior_masses(space)
    assert len(w) == 6
    assert w.sum() == pytest.approx(cell_space(space, m6c.S1_CELLS[0]).slice_weight, rel=1e-12)
    np.testing.assert_allclose(
        w, [0.042884, 0.032163, 0.032163, 0.038361, 0.028771, 0.028771], atol=1e-6
    )


def test_seeds_are_distinct_and_as_planned():
    seeds = [tuple(S.rng_a(k, s)) for k in range(6) for s in S.SEEDS]
    seeds += [tuple(S.rng_b(k, r)) for k in range(6) for r in S.B_RUNS]
    assert len(set(seeds)) == len(seeds) == 36
    assert S.rng_a(2, 3) == [2026, 14, 1846, 1, 3, 3]
    assert S.rng_b(4, 2) == [2026, 14, 1846, 1, 5, 102]
    assert (S.NLIVE, S.SLICES, S.N_DRAWS, S.CHUNK) == (4000, 28, 6000, 250)


def test_the_full_likelihood_space_is_m6c_s_cut_to_the_subcell(space):
    """Sub-cell spaces draw inside their a and e ranges, and the union of the cut is M6c's
    cell-0 job space: the same prior transform outside the two cut coordinates."""
    _, _, _, pq, box = m6c.dataset(variant="s1_e05")
    job = m6c.job_space(m6c.Job("s1_e05", 0), pq, box)
    rng = np.random.default_rng(5)
    for k, ((a0, a1), (e0, e1)) in enumerate(S.SUBCELLS):
        cs = S.space_b(k)
        assert cs.unit_names == job.unit_names
        for _ in range(20):
            u = rng.random(cs.ndim)
            th = cs.theta(cs.prior_transform(u))
            assert a0 <= th[7] <= a1 and e0 <= th[8] <= e1


def test_spread_and_shares():
    assert S.spread([359.9, 0.2, 0.1], circular=True) == pytest.approx(0.3)
    assert S.spread([-362.3, -362.8, -362.4, -362.5]) == pytest.approx(0.5)
    sh = S.evidence_shares(np.array([-360.0, -360.0]), np.array([0.1, 0.3]))
    np.testing.assert_allclose(sh, [0.25, 0.75])


def _stage_a(ln_z0: float, t1: bool = True, g2: bool = True) -> dict:
    return {"ln_z0": ln_z0, "ln_z0_err": 0.1, "T1": {"pass": t1},
            "gates": {"G2": {"pass": g2}, "G3": {"pass": True}}}  # fmt: skip


def test_the_check_values_and_the_bands():
    from fork.export import m8_arbiter

    assert S.p_below_25(-361.94, S.OTHER_FAST) == pytest.approx(
        m8_arbiter.implied_p_below_25(-361.94, list(S.OTHER_FAST)), rel=1e-12
    )
    assert S.p_below_25(-366.8432, S.OTHER_PUBLISHED) == pytest.approx(0.00104, abs=2e-5)
    assert S.p_below_25(-362.565, S.OTHER_FAST) == pytest.approx(0.109, abs=1e-3)
    assert [S.word(p) for p in (0.001, 0.049, 0.05, 0.3, 0.5, 0.9)] == [
        "not at all", "not at all", "moved partly", "moved partly", "moved", "moved"
    ]  # fmt: skip


@pytest.mark.parametrize(
    ("a", "b", "g1", "outcome", "verdict"),
    [
        (_stage_a(-362.0), {"pass": True}, True, "resolved", "moved partly"),
        (_stage_a(-368.0), {"pass": True}, True, "resolved", "not at all"),
        (_stage_a(-362.0, g2=False), {"pass": True}, True, "a gate", "stands"),
        (_stage_a(-362.0), {"pass": True}, False, "a gate", "stands"),
        (_stage_a(-362.0, t1=False), {"pass": True}, True, "T1", "stands"),
        (_stage_a(-362.0), {"pass": False}, True, "T2", "stands"),
        (_stage_a(-362.0), None, True, "T2", "stands"),
    ],
)
def test_the_reading_follows_the_plan(a, b, g1, outcome, verdict):
    r = S.reading(a, b, g1)
    assert outcome in r["outcome"] and verdict in r["verdict"]
    assert r["resolved"] == (outcome == "resolved")
    assert r["label_95"].startswith("provisional")  # under every outcome (§ 5)
    low = r["resolved"] and max(r["p_a_below_25_au"].values()) < 0.05
    assert ("no longer a reason" in r["label_95"]) == low


def test_a_resolved_reading_between_two_words():
    lo, hi = -364.0, -362.0
    for _ in range(60):  # a ln Z₀ whose two P(a < 25 AU) fall either side of 0.05
        mid = 0.5 * (lo + hi)
        p = (S.p_below_25(mid, S.OTHER_PUBLISHED), S.p_below_25(mid, S.OTHER_FAST))
        if min(p) < 0.05 <= max(p):
            break
        lo, hi = (mid, hi) if max(p) < 0.05 else (lo, mid)
    r = S.reading(_stage_a(mid), {"pass": True}, True)
    assert r["verdict"] == "between moved partly and not at all"


def test_within_marks_values_at_the_tolerance():
    assert S.within(0.49, "ln_z") == {"value": 0.49, "tolerance": 0.5, "pass": True,
                                      "at_the_tolerance": True}  # fmt: skip
    assert not S.within(-0.36, "median")["pass"]
    assert not S.within(0.2, "ln_z")["at_the_tolerance"]


def test_smoke_options_are_refused_in_the_run_cache():
    for kw in ({"maxiter": 10}, {"synthetic": 1}):
        with pytest.raises(SystemExit):
            S._guard(S.CACHE, kw.get("maxiter"), kw.get("synthetic"))
    S._guard(Path("/nonexistent/scratch"), 10, 1)


def test_the_plan_names_these_settings():
    plan = (m6c.ROOT / "docs" / "notes" / "S1_cell0.md").read_text(encoding="utf-8")
    frozen = plan.split("\n---\n")[0]
    for text in ("nlive 4000, 28 slices", "`[2026, 14, 1846, 1, k + 1, s]`", "100 + r",
                 "6,000 equal-weight draws", "ESS/N ≥ 0.3", "at least 0.01"):  # fmt: skip
        assert text in frozen, text
    assert math.isclose(S.MATERIAL_SHARE, 0.01) and S.ESS_GATE == 0.3
