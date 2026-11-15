"""The Mercury study with the full forward model, against its plan
(docs/phase2/identifiability_full.md)."""

import csv
import json
import math

import numpy as np
import pytest

from fork.model import inner_system as IS
from fork.phase2 import identifiability as ID
from fork.phase2 import identifiability_full as F
from fork.phase2 import identifiability_study as IDS


def test_inputs_and_plan_are_the_frozen_ones():
    F.check_all_inputs()
    plan = F.PLAN.read_text()
    for rel, digest in F.INPUTS.items():
        assert digest in plan, rel


def test_truths_file_is_the_plans_table():
    t = F.truths()
    plan = F.PLAN.read_text()
    assert t["b"]["extra"][0] == pytest.approx(1.9811152404862693e-07, rel=1e-15)
    assert "1.9811152404862693 × 10⁻⁷" in plan and "3.399891252673649 × 10⁻⁴" in plan
    assert [t[k]["model"] for k in "abcd"] == ["gr", "ring", "planet", "j2"]


def test_settings_are_the_plans():
    plan = F.PLAN.read_text()
    for text in ("800 full integrations", "5 degrees of freedom", "1.5 × the Laplace",
                 "within 0.1 in ln", "below 0.03", "less than 0.05", "below 0.5",
                 "within 20 iterations", "at most 0.1 % of calls", "Δχ² = 25",
                 "959.63″", "3.36″", "more than 1 h", "within 30 min", "10⁻¹² rad",
                 "15 points over ±6 posterior widths", "20 times finer", "0.1 decade",
                 "and ±4 decades from the best value", "e⁻²⁵",
                 "[2026, 21, 1882, truth, seed]"):  # fmt: skip
        assert text in plan, text
    assert (F.V2_DRAWS, F.V2_NU, F.V2_SCALE, F.V2_TOL, F.V2_ERR, F.V2_FRAC) == (
        800,
        5.0,
        1.5,
        0.1,
        0.03,
        0.95,
    )
    assert (F.V3_TOL, F.V4_TOL, F.MAX_IT, F.V6_FRAC, F.V6_DCHI) == (0.05, 0.5, 20, 0.001, 25.0)
    assert (F.FINE_N, F.FINE_HALF, F.SUB, F.COARSE_STEP, F.COARSE_HALF) == (15, 6.0, 20, 0.1, 4.0)
    assert F.TAIL_DECADES == (1.0, 2.0, 4.0) and F.TAIL_LN == 25.0
    assert (F.LEFT_OUT_S, F.V0_MAX_S, F.ROOT_TOL, F.DCHI_TOL) == (3600.0, 1800.0, 1e-12, 1e-3)


def test_tasks():
    s1, s2 = F.stage1_tasks(), F.stage2_tasks()
    assert len(s1) == 4 * (3 + 3) + 1 + 1 + 4 and len({t.id for t in s1}) == len(s1)
    assert len(s2) == 8 * 4 * 6 and all(t.seed > 0 and t.level == 0 for t in s2)
    for t in s1 + s2:
        assert F.Task.from_id(t.id) == t
    half = [t for t in s1 if t.half]
    assert {(t.truth, t.model) for t in half} == {(2, "R"), (3, "P"), (4, "J")}
    assert all(t.distance == 1 for t in half if t.model == "P")  # the planet's true distance
    lam = F.Task(1, 0, 0, "P", 0, 0).second_grid()
    assert np.allclose(lam, 2 * math.pi * np.arange(48) / 48)
    lam2 = np.concatenate([F.Task(3, 0, 0, "P", 1, b, True).second_grid() for b in range(4)])
    assert np.allclose(lam2, 2 * math.pi * np.arange(96) / 96)
    assert np.allclose(F.Task(2, 0, 0, "R", half=True).second_grid(), np.linspace(0.1, 0.3, 21))
    assert np.allclose(F.Task(2, 0, 0, "R").second_grid(), IDS.RING_GRID)
    assert F.Task(1, 0, 0, "J").levels() == IDS.LEVELS and F.Task(1, 3, 0, "J").levels() == (1.0,)


def test_the_integration_is_inner_systems():
    theta6, ctx = F.base_problem()
    jd = np.array([2385000.5, 2400000.5, 2410000.5])
    for model, extra in (("gr", (1.0,)), ("planet", (2e-7, 0.2, 0.3))):
        x, v, _, _ = F.states(model, theta6, extra, ctx, jd)
        xr, vr = IS.states(np.concatenate([theta6, extra]), jd, ctx, ID.flags_for(model))
        np.testing.assert_allclose(x, xr, rtol=0, atol=1e-15)
        np.testing.assert_allclose(v, vr, rtol=0, atol=1e-17)


def test_contact_solver_on_a_straight_crossing(monkeypatch):
    """Mercury crossing the Sun's disc in a straight line seen from a fixed Earth: the root is
    found to 10⁻¹² rad, at the contact the small-angle geometry predicts (to a few seconds,
    the size of the small-angle error)."""
    monkeypatch.setattr(F, "MU_MERCURY", 0.0)  # straight lines
    monkeypatch.setattr(F, "MU_EMB", 0.0)
    e, ve = np.array([[1.0, 0.0, 0.0]]), np.zeros((1, 3))
    d, b, speed = 0.6, 0.001, 1e-3  # distance from the Earth [AU], impact [rad], rad per day
    for exterior, sign in ((False, -1.0), (True, 1.0)):
        rc = F.R_SUN + sign * F.R_MERCURY / d
        s_contact = -math.sqrt(rc**2 - b**2)
        s0 = s_contact + 0.002 * speed  # the printed epoch 0.002 d after the contact
        xm = e + d * np.array([[-1.0, s0, b]])
        vm = d * np.array([[0.0, speed, 0.0]])
        off, res = F.contact_offsets(xm, vm, e, ve, np.array([exterior]))
        assert res[0] < 1e-12
        assert off[0] == pytest.approx(-0.002 * 86400.0, abs=5.0)


def test_lm_fit_cold_and_warm_reach_the_same_optimum():
    t = np.linspace(0.0, 4.0, 40)
    truth = np.array([2.0, -0.7, 0.3])
    y = truth[0] * np.exp(truth[1] * t) + truth[2]

    def resid(p):
        return (p[0] * np.exp(p[1] * t) + p[2] - y) / 0.01

    cold = F.lm_fit(resid, np.array([1.5, -0.5, 0.0]), np.full(3, 1e-7))
    assert cold.converged and cold.chi2 < 1e-6
    np.testing.assert_allclose(cold.theta6, truth, atol=1e-6)
    warm = F.lm_fit(resid, cold.theta6 + 0.01, np.full(3, 1e-7), jac0=cold.jac)
    assert warm.converged and warm.chi2 < 1e-6
    _, s, _ = np.linalg.svd(cold.jac)
    assert cold.ln_s == pytest.approx(float(np.sum(np.log(s))), rel=1e-12)


def test_no_contact_gives_zero_likelihood():
    def resid(p):
        raise F.NoContact

    c = F.lm_fit(resid, np.zeros(6), F.EL_STEPS)
    assert c.no_contact and c.ln_l(1.0) == -math.inf and not c.converged


def _profile(x0=-6.7, w=0.013, ln_s=3.0):
    xs_c = x0 + np.arange(-4.0, 4.0001, 0.1)
    xs_f = x0 + w * np.linspace(-14, 14, 33)
    pts = [{"x": float(x), "kind": "coarse", "chi2": ((x - x0) / w) ** 2, "ln_s": ln_s,
            "no_contact": False} for x in xs_c]  # fmt: skip
    pts += [{"x": float(x), "kind": "fine", "chi2": ((x - x0) / w) ** 2, "ln_s": ln_s,
             "no_contact": False} for x in xs_f]  # fmt: skip
    return pts, xs_f[1] - xs_f[0]


@pytest.mark.parametrize("level", (1.0, 0.5, 2.0))
def test_strength_quadrature_on_a_gaussian(level):
    w, ln_s = 0.013, 3.0
    pts, h = _profile(w=w, ln_s=ln_s)
    got = F.ln_z_strength(pts, level, (-9.0, -5.0), h)
    want = math.log(math.sqrt(2 * math.pi) * level * w) - ln_s + 6 * math.log(level) - math.log(4.0)
    assert got == pytest.approx(want, abs=1e-6)


def test_strength_quadrature_respects_the_box():
    pts, h = _profile(x0=-6.7)
    full = F.ln_z_strength(pts, 1.0, (-9.0, -5.0), h)
    half = F.ln_z_strength(pts, 1.0, (-6.7, -5.0), h)  # cuts the peak in two
    assert half == pytest.approx(full + math.log(0.5) + math.log(4.0 / 1.7), abs=1e-3)
    assert F.ln_z_strength(pts, 1.0, (-12.0, -11.0), h) < full - 1e3  # far from the peak


def _fake(results_scans):
    return {"scans": [{"p2": p2, "ln_z": {"1.0": {"prior": z}}} for p2, z in results_scans]}


def test_combining_radii_and_longitudes():
    radii = np.linspace(0.1, 0.3, 11)
    res = {F.Task(2, 0, 0, "R").id: _fake([(r, -10.0 + 5 * r) for r in radii])}
    got = F.ln_z_model(res, 2, 0, 0, "R", 1.0)
    z = np.exp(-10.0 + 5 * radii)
    assert got == pytest.approx(math.log(np.trapezoid(z, radii) / 0.2), rel=1e-12)
    for d in range(3):
        lam = F.Task(2, 0, 0, "P", d, 0).second_grid()
        res[F.Task(2, 0, 0, "P", d, 0).id] = _fake([(x, float(d) + 0.0 * x) for x in lam])
    got = F.ln_z_model(res, 2, 0, 0, "P", 1.0)
    assert got == pytest.approx(math.log(np.mean(np.exp([0.0, 1.0, 2.0]))), rel=1e-12)


def test_data_depend_only_on_epochs_kinds_weights_truth_and_seed(tmp_path):
    """V5: every other printed column can change without changing what the study reads."""
    rows = list(csv.DictReader(IDS.CONTACTS.open(encoding="utf-8")))
    for r in rows:
        r.update(k_sign="+", k_coeff="99", second_value="1", flag="x", printed_page="1",
                 contact="")  # fmt: skip
    p = tmp_path / "c.csv"
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    a, b = F.contact_rows(), F.contact_rows(p)
    for k in ("jd", "weight", "exterior"):
        np.testing.assert_array_equal(a[k], b[k])
    y = np.arange(5.0)
    s = np.full(5, 2.0)
    assert np.array_equal(F.synthetic_data(y, s, 1, 0), y)
    assert np.array_equal(F.synthetic_data(y, s, 2, 3), F.synthetic_data(y, s, 2, 3))
    assert not np.array_equal(F.synthetic_data(y, s, 2, 3), F.synthetic_data(y, s, 2, 4))


def test_exterior_contacts_are_tables_three_and_four():
    r = F.contact_rows()
    assert r["exterior"].sum() == 15 + 6 and len(r["jd"]) == 58
    assert set(r["table"][r["exterior"]]) == {"III", "IV"}


def test_stage2_refuses_without_stage1(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "CACHE", tmp_path)
    with pytest.raises(SystemExit, match="Stage 1"):
        F.stage2(1)
    (tmp_path / "stage1.json").write_text(json.dumps({"gates": {"pass": False}}))
    with pytest.raises(SystemExit, match="Stage 1"):
        F.stage2(1)


def _write_points(task, rows):
    task.dir.mkdir(parents=True)
    pts = F.Points(task.dir / "points.jsonl")
    for x, chi2, conv, nc in rows:
        pts.calls[pts.key(0.0, x)] = {"p2": 0.0, "x": x, "kind": "fine", "chi2": chi2,
                                      "converged": conv, "no_contact": nc}  # fmt: skip
    (task.dir / "points.jsonl").write_text(
        "\n".join(json.dumps(v) for v in pts.calls.values()) + "\n"
    )


def test_v6_frozen_and_amended(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "CACHE", tmp_path)
    monkeypatch.setattr(F, "V6_AMENDED", False)
    j = F.Task(4, 0, 0, "J")
    rows = [(float(i), 1.0 + i, True, False) for i in range(2000)]
    rows += [(-1.0, 500.0, False, False), (-2.0, math.inf, False, True)]
    _write_points(j, rows)
    ng = F.Task(1, 0, 0, "NG")
    _write_points(ng, [(0.0, 1107.5, True, False), (1.0, 0.0, True, False)])
    results = {j.id: {"scans": [{"p2": 0.0, "chi2_min": 1.0}]}, ng.id: {}}
    frozen = F.v6_check(results)
    assert frozen["no_contact_calls"] == 1 and frozen["calls"] == 2003
    assert frozen["not_converged"] == 1 and frozen["far_not_converged"] == 1
    assert frozen["fraction_all"] == pytest.approx(1 / 2003) and not frozen["near_best"]
    assert frozen["pass"] is True  # 1/2003 < 0.1 %
    monkeypatch.setattr(F, "V6_AMENDED", True)
    amended = F.v6_check(results)
    assert amended["not_converged_within_far"] == 0 and amended["pass"] is True
    _write_points(F.Task(2, 0, 0, "NG"), [(0.0, 900.0, False, False)])  # N fails at its own best
    results[F.Task(2, 0, 0, "NG").id] = {}
    worse = F.v6_check(results)
    assert len(worse["near_best"]) == 1 and worse["pass"] is False


def test_amendment_one_is_in_force():
    assert F.JAC_MODE == "whitened" and F.V6_AMENDED and F.V6_FAR == 200.0
    below = F.PLAN.read_text()[len(F._plan_frozen_part()) :]
    assert "## Amendment 1" in below and "Δχ² = 200" in below
