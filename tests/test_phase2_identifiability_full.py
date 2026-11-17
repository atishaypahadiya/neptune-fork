"""The Mercury study with the full forward model, against its plan
(docs/phase2/identifiability_full.md)."""

import csv
import json
import math

import numpy as np
import pytest

from fork.data.frames import elements_to_state
from fork.model import inner_system as IS
from fork.phase2 import contact_root as CRT
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
    # distance from the Earth [AU], impact [rad], rad per day (a real transit's rate: the disc,
    # 0.0093 rad, in about 7 h; amendment 2 looks for mid-transit within a day of the epoch)
    d, b, speed = 0.6, 0.001, 0.03
    for exterior, sign in ((False, -1.0), (True, 1.0)):
        rc = F.R_SUN + sign * F.R_MERCURY / d
        s_contact = -math.sqrt(rc**2 - b**2)
        s0 = s_contact + 0.002 * speed  # the printed epoch 0.002 d after the contact
        xm = e + d * np.array([[-1.0, s0, b]])
        vm = d * np.array([[0.0, speed, 0.0]])
        off, res = F.contact_offsets(xm, vm, e, ve, np.array([exterior]), np.array([-1]))
        assert res[0] < 1e-12
        assert off[0] == pytest.approx(-0.002 * 86400.0, abs=5.0)


def _crossing(b: float, t_mid: float, speed: float = 0.03, d: float = 0.6):
    """Straight-line crossing seen from a fixed Earth: impact b [rad], mid-transit at t_mid
    [d] from the epoch, at a real transit's rate. Returns the states for ``contact_offsets``."""
    e, ve = np.array([[1.0, 0.0, 0.0]]), np.zeros((1, 3))
    xm = e + d * np.array([[-1.0, -t_mid * speed, b]])
    vm = d * np.array([[0.0, speed, 0.0]])
    return xm, vm, e, ve


def test_contact_solver_mid_transit_at_the_epoch(monkeypatch):
    """Amendment 2's regression: mid-transit exactly at the printed epoch, where the gap's slope
    is zero (Newton's first step left the window there). Both contacts are found, ingress
    before and egress after, symmetric about the epoch."""
    monkeypatch.setattr(F, "MU_MERCURY", 0.0)
    monkeypatch.setattr(F, "MU_EMB", 0.0)
    st = _crossing(b=0.001, t_mid=0.0)
    for ext in (False, True):
        t_in, r_in = F.contact_offsets(*st, np.array([ext]), np.array([-1]))
        t_out, r_out = F.contact_offsets(*st, np.array([ext]), np.array([1]))
        assert r_in[0] < 1e-12 and r_out[0] < 1e-12
        assert t_in[0] < 0.0 < t_out[0]
        assert t_in[0] == pytest.approx(-t_out[0], abs=0.5)  # small-angle asymmetry only


def test_contact_solver_chooses_the_side_not_the_nearer_root(monkeypatch):
    """Mid-transit 2 h after the epoch: the ingress contact is the one asked for even where the
    egress contact lies nearer the epoch, and vice versa (Newton took the nearer one)."""
    monkeypatch.setattr(F, "MU_MERCURY", 0.0)
    monkeypatch.setattr(F, "MU_EMB", 0.0)
    for t_mid in (2.0 / 24.0, -2.0 / 24.0):
        st = _crossing(b=0.001, t_mid=t_mid)
        t_in, _ = F.contact_offsets(*st, np.array([False]), np.array([-1]))
        t_out, _ = F.contact_offsets(*st, np.array([False]), np.array([1]))
        assert t_in[0] < t_mid * 86400.0 < t_out[0]


def test_contact_solver_no_contact_only_without_a_root(monkeypatch):
    """A miss gives no contact of either kind; a grazing transit that crosses the limb but never
    lies wholly inside gives exterior contacts and no interior ones."""
    monkeypatch.setattr(F, "MU_MERCURY", 0.0)
    monkeypatch.setattr(F, "MU_EMB", 0.0)
    d = 0.6
    miss = _crossing(b=F.R_SUN + 2.0 * F.R_MERCURY / d, t_mid=0.0, d=d)
    graze = _crossing(b=F.R_SUN, t_mid=0.0, d=d)
    for side in (-1, 1):
        for ext in (False, True):
            off, _ = F.contact_offsets(*miss, np.array([ext]), np.array([side]))
            assert np.isnan(off[0])
        off_i, _ = F.contact_offsets(*graze, np.array([False]), np.array([side]))
        off_e, _ = F.contact_offsets(*graze, np.array([True]), np.array([side]))
        assert np.isnan(off_i[0]) and np.isfinite(off_e[0])


def test_contact_root_refuses_an_unstated_side():
    with pytest.raises(ValueError):
        CRT.contact_roots(lambda i, t: t**2 - 1e-4, np.array([0]), 1e-12)


def test_contact_root_matches_a_fine_grid_on_a_parabola():
    """g(t) = (t − m)² − w²: roots at m ± w, found to the double's precision."""
    m, w = np.array([0.3, -0.7, 0.0]), np.array([0.05, 0.01, 0.2])

    def gap(i, t):
        return (t - m[i]) ** 2 - w[i] ** 2

    t_in, r_in = CRT.contact_roots(gap, np.full(3, -1), 1e-12)
    t_out, _ = CRT.contact_roots(gap, np.full(3, 1), 1e-12)
    np.testing.assert_allclose(t_in, m - w, rtol=0, atol=1e-13)
    np.testing.assert_allclose(t_out, m + w, rtol=0, atol=1e-13)
    assert np.all(r_in < 1e-12)


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
            "no_contact": False, "converged": True} for x in xs_c]  # fmt: skip
    pts += [{"x": float(x), "kind": "fine", "chi2": ((x - x0) / w) ** 2, "ln_s": ln_s,
             "no_contact": False, "converged": True} for x in xs_f]  # fmt: skip
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
    for k in ("jd", "weight", "exterior", "contact", "side"):
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


# ---------------------------------------------------------------------------- domain's edge


def _position(theta):
    """A smooth residual that, like the model's, refuses elements outside a > 0, 0 ≤ e < 1."""
    if not F._in_domain(theta):
        raise F.OutsideDomain
    return elements_to_state(*theta, F.MU_MERCURY)[0]


def test_difference_is_the_old_forward_difference_inside_the_domain():
    theta = np.array([0.387, 0.2056, 0.12, 0.85, 0.51, 1.8])
    r = _position(theta)
    for d in (np.array([0.0, 2e-7, 0, 0, 0, 0]), np.array([1e-4, -3e-3, 2e-3, 1e-3, -1e-3, 5e-3])):
        n0 = F.N_EDGE[0]
        assert np.array_equal(F._difference(_position, theta, r, d), _position(theta + d) - r)
        assert F.N_EDGE[0] == n0
    steps = np.array([1e-10, 2e-7, 2e-7, 2e-7, 2e-7, 2e-7])
    old = np.array([(_position(theta + h * np.eye(6)[k]) - r) / h for k, h in enumerate(steps)]).T
    assert np.array_equal(F._jacobian(_position, theta, r, steps), old)


def test_near_circular_orbit_takes_the_backward_step():
    """e = 1e-4 and a step that lowers e by 2e-4: the forward trial has e < 0, which crashed
    the old code; the step is taken backward, into the domain, at its full length."""
    theta = np.array([0.387, 1e-4, 0.12, 0.85, 0.51, 1.8])
    d = np.array([0.0, -2e-4, 1e-5, 0.0, 0.0, 0.0])
    r = _position(theta)
    with pytest.raises(ValueError):
        elements_to_state(*(theta + d), F.MU_MERCURY)  # what the old forward difference hit
    n0 = F.N_EDGE[0]
    diff = F._difference(_position, theta, r, d)
    assert F.N_EDGE[0] == n0 + 1
    assert np.array_equal(diff, r - _position(theta - d))
    # a backward difference of the same length: agrees with the derivative to first order
    h = 1e-9
    deriv = (_position(theta + h * d / np.linalg.norm(d)) - r) / h * np.linalg.norm(d)
    np.testing.assert_allclose(diff, deriv, rtol=0.05, atol=1e-9)


def test_collapsed_direction_is_shortened_into_the_domain():
    """The Stage 1 crash (t3_s0_l0_P_d0_b0): a Jacobian whose smallest singular value had
    fallen to 1.3e-6 made one whitened unit a step of 4e5 in e and -1e3 AU in a, outside the
    domain on both sides. The step is halved along d until a trial lies inside."""
    theta = np.array([0.38709795748, 0.20375610234, 0.12415470373, 0.84568510339,
                      0.50881257509, 1.81576156722])  # fmt: skip
    d = np.array([-987.545, 373144.709, -54546.781, 407388.744, 282047.555, 426424.088])
    r = _position(theta)
    for sign in (1, -1):
        with pytest.raises(ValueError):
            elements_to_state(*(theta + sign * d), F.MU_MERCURY)
    n0 = F.N_EDGE[0]
    diff = F._difference(_position, theta, r, d)
    assert F.N_EDGE[0] == n0 + 1
    assert np.all(np.isfinite(diff))
    f = 1.0
    while not (F._in_domain(theta + f * d) or F._in_domain(theta - f * d)):
        f *= 0.5
    trial = theta + f * d if F._in_domain(theta + f * d) else theta - f * d
    sign = 1.0 if F._in_domain(theta + f * d) else -1.0
    np.testing.assert_allclose(diff, sign * (_position(trial) - r) / f, rtol=0, atol=0)


def test_whitened_jacobian_survives_a_collapsed_basis():
    theta = np.array([0.387, 0.2038, 0.124, 0.846, 0.509, 1.816])
    rng = np.random.default_rng(0)
    q, _ = np.linalg.qr(rng.normal(size=(6, 6)))
    basis = q @ np.diag([5e10, 4e9, 7e5, 2.5e5, 1.2e5, 1.3e-6]) @ q.T
    r = _position(theta)
    jac = F._whitened_jacobian(_position, theta, r, basis)
    assert jac.shape == (3, 6) and np.all(np.isfinite(jac))


def test_contact_identity_file_against_the_table():
    """Amendment 2: the identity file covers every row, repeats each printed number except the
    misprint I 11 (contact II), and makes every exterior row contact IV."""
    rows = F.contact_rows()
    lab = [f"{t} {r}" for t, r in zip(rows["table"], rows["row"], strict=True)]
    with IDS.CONTACTS.open(newline="", encoding="utf-8") as fh:
        printed = [r["contact"] for r in csv.DictReader(fh)]
    with F.CONTACT_IDENTITY.open(newline="", encoding="utf-8") as fh:
        assert [r["printed_contact"] for r in csv.DictReader(fh)] == printed
    changed = [n for n, p, c in zip(lab, printed, rows["contact"], strict=True) if p and p != c]
    assert changed == ["I 11"]
    assert set(rows["contact"][rows["exterior"]]) == {"IV"}
    assert set(rows["contact"][~rows["exterior"]]) == {"II", "III"}


@pytest.mark.slow
def test_contact_identity_agrees_with_truth_a():
    """Amendment 2's check, not its rule: on truth a's noise-free model, the contact nearer each
    printed epoch is the identity file's on every row but III 12, which is also the only row
    V0 leaves out (its printed day is one before table I's for the same transit, so both of
    truth a's contacts lie 18 to 24 h after it and the nearer one is ingress)."""
    rows = F.contact_rows()
    theta6, ctx = F.base_problem()
    a = F.truths()["a"]
    st = F.states(a["model"], theta6, a["extra"], ctx, rows["jd"])
    n = len(rows["jd"])
    off_in, _ = F.contact_offsets(*st, rows["exterior"], np.full(n, CRT.INGRESS))
    off_out, _ = F.contact_offsets(*st, rows["exterior"], np.full(n, CRT.EGRESS))
    near = F.nearer_side(off_in, off_out)
    lab = np.array([f"{t} {r}" for t, r in zip(rows["table"], rows["row"], strict=True)])
    assert lab[near != rows["side"]].tolist() == ["III 12"]
    off, _ = F.contact_offsets(*st, rows["exterior"], rows["side"])
    left = lab[~(np.isfinite(off) & (np.abs(off) <= F.LEFT_OUT_S))].tolist()
    assert left == ["III 12"]


def test_carry_is_the_exact_two_body_motion():
    """Amendment 2: the contact solver carries Mercury and the Earth from the epoch by the exact
    two-body motion about the Sun at every offset in its window (to a day), checked against an
    independent integration of the two-body equations."""
    from scipy.integrate import solve_ivp

    x0 = np.array([[0.30, -0.12, 0.03], [-0.6, 0.78, 0.0]])
    v0 = np.array([[0.008, 0.028, 0.002], [-0.014, -0.011, 0.0]])
    for mu in (F.MU_MERCURY, F.MU_EMB):
        for t in (2.0 / 24.0, 0.5, 1.0):

            def f(_, y, mu=mu):
                return np.concatenate([y[3:], -mu * y[:3] / np.linalg.norm(y[:3]) ** 3])

            got = F._carried(x0, v0, mu, np.full(2, t))
            for i in range(2):
                ref = solve_ivp(f, (0.0, t), np.concatenate([x0[i], v0[i]]), rtol=1e-13,
                                atol=1e-15).y[:3, -1]  # fmt: skip
                np.testing.assert_allclose(got[i], ref, rtol=0, atol=1e-12)


@pytest.mark.slow
def test_iii_12_is_left_out_of_the_likelihood(tmp_path, monkeypatch):
    """III 12 is the only row V0 leaves out, and it is not among the observation's rows, so no
    likelihood call sees it."""
    monkeypatch.setattr(F, "CACHE", tmp_path)
    prep = F.prepare()
    assert prep["V0"]["left_out"] == ["III 12"]
    rows = F.contact_rows()
    obs = F.observation(prep)
    k = [f"{t} {r}" for t, r in zip(rows["table"], rows["row"], strict=True)].index("III 12")
    assert len(obs.jd) == len(rows["jd"]) - 1 == prep["n_kept"]
    assert not np.any(obs.jd == rows["jd"][k])


class _StubData:
    """A likelihood in the strength alone: χ²(x) given, elements returned as started, and no
    contact from any start whose first element is in ``bad_a``."""

    def __init__(self, chi2_of_x, bad_a=()):
        self.chi2_of_x, self.bad_a, self.starts = chi2_of_x, tuple(bad_a), []

    def call(self, model, extra, start, jac0=None):
        x = math.log10(extra[0])
        self.starts.append(np.array(start))
        if start[0] in self.bad_a:
            return F.Call(math.inf, 0.0, np.array(start), None, 0, 0, False, no_contact=True)
        return F.Call(float(self.chi2_of_x(x)), 0.0, np.array(start), None, 1, 1, True)


def _stub_scanner(tmp_path, monkeypatch, data, model="P"):
    monkeypatch.setattr(F, "CACHE", tmp_path)
    task = F.Task(3, 0, 0, model, 0, 0) if model == "P" else F.Task(4, 0, 0, model)
    start = np.array([0.387, 0.2056, 0.12, 0.84, 1.35, 4.4])
    return F.Scanner(task, data, start), task


def test_scan_stops_at_the_box_when_the_data_push_out_of_it(tmp_path, monkeypatch):
    """Amendment 3 (A, a): with χ² still falling below the box (best value at x = −15), no
    point is computed outside the box, the best point is its lower edge, reported so, and the
    fine points there are no wider apart than the coarse scan's."""
    sc, _ = _stub_scanner(tmp_path, monkeypatch, _StubData(lambda x: 0.5 * (x + 15.0) ** 2))
    s = sc.scan_second(0.0, None, None)
    xs = sorted(p["x"] for p in sc.points.calls.values())
    lo, hi = F.box("P", max(F.R3_WIDTHS))
    assert lo - 1e-12 <= xs[0] and xs[-1] <= hi + 1e-12
    assert s["best_at_box_edge"] == "low" and s["best_point"]["x"] == pytest.approx(lo)
    fine = sorted(p["x"] for p in sc.points.calls.values() if p["kind"] == "fine")
    assert fine[0] == pytest.approx(lo) and fine[1] - fine[0] <= F.COARSE_STEP + 1e-12
    assert len(sc.points.calls) < 400


def test_scan_inside_the_box_is_unchanged(tmp_path, monkeypatch):
    """Amendment 3 (A): where the fine points stay inside the box they are the plan's,
    x0 + k·h, and no best value is reported at an edge."""
    sc, _ = _stub_scanner(tmp_path, monkeypatch, _StubData(lambda x: 400.0 * (x + 7.2) ** 2))
    s = sc.scan_second(0.0, None, None)
    assert s["best_at_box_edge"] is None
    fine = np.array(sorted(p["x"] for p in sc.points.calls.values() if p["kind"] == "fine"))
    k = (fine - fine[0]) / s["h"]
    np.testing.assert_allclose(k, np.round(k), atol=1e-6)
    assert s["x_best"] == pytest.approx(-7.2, abs=1e-3)


def test_a_start_with_no_contact_is_retried(tmp_path, monkeypatch):
    """Amendment 3 (B): a call whose start has no contact is retried from the nearest
    converged call of the scan with contacts, then the base elements; it has no contact only
    if every start has none, and its retries are recorded."""
    data = _StubData(lambda x: (x + 7.0) ** 2, bad_a=(9.0,))
    sc, _ = _stub_scanner(tmp_path, monkeypatch, data)
    good = sc.call(0.0, -7.0, "coarse")
    assert good["retries"] == 0 and not good["no_contact"]
    bad = np.array([9.0, 0.2, 0.1, 0.8, 1.3, 4.4])
    d = sc.call(0.0, -6.5, "tail", (bad, None))
    assert d["retries"] == 1 and not d["no_contact"] and d["chi2"] == pytest.approx(0.25)
    np.testing.assert_array_equal(data.starts[-1], good["theta6"])

    data2 = _StubData(lambda x: (x + 7.0) ** 2, bad_a=(9.0, 0.387))  # the base fails too
    sc2, _ = _stub_scanner(tmp_path / "b", monkeypatch, data2)
    d2 = sc2.call(0.0, -6.0, "coarse")
    assert d2["no_contact"] and d2["retries"] == 0  # the start was the base: nothing else
    d3 = sc2.call(0.0, -5.0, "tail", (bad, None))
    assert d3["no_contact"] and d3["retries"] == 1  # the base, after the given start


def test_no_contact_inside_the_edge_fallback_is_skipped():
    """The smoke-test crash of t3_s0_l0_P_d1_b2_half (amendment 3, addendum): the forward
    trial of a whitened column fell outside the elements' domain, which sends the column to
    the fallback, and the backward trial there had no contact, which the fallback did not
    catch. Now a trial with no contact inside the fallback is skipped like one outside the
    domain, the step is halved until a trial is allowed, and the column is counted."""
    theta = np.array([0.387, 0.2056, 0.12, 0.84, 1.35, 4.4])
    d = np.array([-0.5, 0.0, 0.0, 0.0, 0.0, 3.0])  # forward: a = −0.113, outside the domain

    def resid(t):
        if not F._in_domain(t):
            raise F.OutsideDomain
        if t[5] < theta[5] - 1.0:  # backward: λ moved by more than 1 rad, no contact
            raise F.NoContact
        return _position(t)

    r = resid(theta)
    with pytest.raises(F.OutsideDomain):
        resid(theta + d)
    with pytest.raises(F.NoContact):
        resid(theta - d)
    n0 = F.N_EDGE[0]
    diff = F._difference(resid, theta, r, d)  # crashed with NoContact before the addendum
    assert F.N_EDGE[0] == n0 + 1
    f = 0.5  # the first allowed trial: forward at half the step (a = 0.137)
    np.testing.assert_allclose(diff, (_position(theta + f * d) - r) / f, rtol=0, atol=0)

    def nowhere(t):  # forward of theta outside the domain, elsewhere no contact
        raise F.OutsideDomain if t[0] < theta[0] else F.NoContact

    with pytest.raises(ValueError, match="no step along d reaches an allowed trial"):
        F._difference(nowhere, theta, r, d)


def _oscillating_fine_points():
    """The shape of t3_s0_l0_P_d1_b2_half's fine points in the third smoke test: χ² rising from
    the box's lower edge, a cluster of points 4e-5 apart with χ² noise of ±0.2 (refinement
    chasing a false minimum), a gap, one unconverged call at χ² 3.3e7, one with no contact."""
    rng = np.random.default_rng(0)

    def chi2(x):
        return 485.0 + 120.0 * (x + 11.0) ** 2

    pts = []

    def add(x, c, conv=True, nc=False):
        pts.append({"x": float(x), "kind": "fine", "chi2": float(c), "ln_s": 0.0,
                    "converged": conv, "no_contact": nc})  # fmt: skip

    for x in np.arange(-11.0, -6.85, 0.054):
        add(x, chi2(x))
    for x in np.arange(-6.8418, -6.8294, 4e-5):
        add(x, chi2(x) + rng.uniform(-0.2, 0.2))
    for x in np.arange(-6.8256, -6.76, 0.0054):
        add(x, chi2(x))
    add(-6.6128, chi2(-6.6128))
    add(-6.064, 3.3e7, conv=False)
    add(-5.516, math.inf, conv=False, nc=True)
    return pts


def test_spline_does_not_oscillate_on_unconverged_and_crowded_points():
    """Amendment 3, second addendum: the spline through every finite fine point (the code
    before it) swings to χ² ≈ −6e4 near x = −6.66, as in the smoke test; from the knots alone
    (converged, with contacts, no two closer than KNOT_TOL) the best value is the lowest point,
    at the box's edge, and the spline and trapezoid integrals agree."""
    from scipy.interpolate import CubicSpline

    pts = _oscillating_fine_points()
    fin = sorted((p for p in pts if math.isfinite(p["chi2"])), key=lambda p: p["x"])
    x = np.array([p["x"] for p in fin])
    sub = np.linspace(x[0], x[-1], (len(x) - 1) * F.SUB + 1)
    assert CubicSpline(x, [p["chi2"] for p in fin])(sub).min() < -1e4  # the old failure

    kn = F.knots(pts)
    assert all(p["converged"] and not p["no_contact"] for p in kn)
    assert np.all(np.diff([p["x"] for p in kn]) >= F.KNOT_TOL)
    xb, _, smin = F.Scanner._best(pts, with_min=True)
    lowest = min(kn, key=lambda p: p["chi2"])
    assert xb == pytest.approx(-11.0) and lowest["chi2"] - smin <= 1.0
    for lv in (0.5, 1.0, 2.0):
        a = F.ln_z_strength(pts, lv, (-11.0, -3.0), 0.054)
        b = F.ln_z_strength(pts, lv, (-11.0, -3.0), 0.054, method="trapezoid")
        assert abs(a - b) < 0.05


def test_fine_spacing_has_a_floor(tmp_path, monkeypatch):
    """Amendment 3, second addendum: however narrow the width, fine points of one grid are at
    least H_MIN apart."""
    sc, _ = _stub_scanner(tmp_path, monkeypatch, _StubData(lambda x: 1e8 * (x + 7.2) ** 2))
    s = sc.scan_second(0.0, None, None)
    assert s["h"] >= F.H_MIN - 1e-15
    fine = sorted(p["x"] for p in sc.points.calls.values() if p["kind"] == "fine")
    assert len(fine) < 2000


def _flank_points():
    """A box edge on a steep flank, as in t3_s0_l0_P_d1_b2_half's 1-decade box: χ² a parabola
    with its minimum at x = −11, below the box (−7.5, −6.5); knots every 0.05 decade, the one at
    the box's lower edge stored 1e-11 outside it, as the coarse grid's floating point left it,
    and one uneven knot near the minimum, as refinement leaves them."""
    pts = [{"x": -10.97, "kind": "fine", "chi2": 485.0 + 120.0 * 0.03**2, "ln_s": 0.0,
            "converged": True, "no_contact": False}]  # fmt: skip
    for k in range(161):
        x = -11.0 + 0.05 * k
        if abs(x + 7.5) < 1e-9:
            x = -7.5 - 1e-11
        pts.append({"x": x, "kind": "fine", "chi2": 485.0 + 120.0 * (x + 11.0) ** 2,
                    "ln_s": 0.0, "converged": True, "no_contact": False})  # fmt: skip
    return pts


def test_box_edges_are_points_of_the_integral():
    """Amendment 3, third addendum: a knot within EDGE_EPS of a box edge is the knot at it, and
    the edge is a point of the integral. Before it, the edge was not a point of the grid and
    ln Z of a box cut on a steep flank depended on where the sub-grid fell: 0.18 too low here
    (0.53 and 0.76 in the fourth smoke test's 1-decade boxes); now it matches the exact Gaussian
    integral within 0.01, as does the box containing the peak."""
    from scipy.special import log_ndtr

    pts = _flank_points()
    lo, hi = -7.5, -6.5
    sig = 1.0 / math.sqrt(120.0)  # ln L = −242.5 − 60 (x + 11)²  at level 1

    def exact(a, b):
        la, lb = log_ndtr(-(a + 11.0) / sig), log_ndtr(-(b + 11.0) / sig)
        return (-242.5 + 0.5 * math.log(2 * math.pi) + math.log(sig)
                + la + math.log(-math.expm1(lb - la)) - math.log(b - a))  # fmt: skip

    got = F.ln_z_strength(pts, 1.0, (lo, hi), 0.05)
    assert got == pytest.approx(exact(lo, hi), abs=0.01)
    kn = F._box_knots(sorted(pts, key=lambda p: p["x"]), (lo, hi))
    assert any(p["x"] == lo for p in kn)  # exactly the edge, by construction
    wide = F.ln_z_strength(pts, 1.0, (-11.0, -3.0), 0.05)
    assert wide == pytest.approx(exact(-11.0, -3.0), abs=0.01)


class _StubP:
    """A planet likelihood in strength and longitude: χ² = 400 (x − x₀(λ))², and no contact
    from any start where ``nc(x, λ)``."""

    def __init__(self, x0_of, nc):
        self.x0_of, self.nc = x0_of, nc

    def call(self, model, extra, start, jac0=None):
        x, lam = math.log10(extra[0]), extra[2]
        if self.nc(x, lam):
            return F.Call(math.inf, 0.0, np.array(start), None, 0, 0, False, no_contact=True)
        chi2 = 400.0 * (x - self.x0_of(lam)) ** 2
        return F.Call(chi2, 0.0, np.array(start), None, 1, 1, True)


def test_a_guess_with_no_knot_falls_back_to_the_coarse_scan(tmp_path, monkeypatch):
    """Amendment 3, fourth addendum: the crash of t3_s0_l0_P_d0_b0 at λ = 7.5°. The fine points
    placed from the previous longitude's best (x = −6.26) all have no contact; before it, the
    best value of a grid with no knot raised IndexError. Now the longitude is scanned as the
    first one is, and the restart is counted."""
    data = _StubP(lambda lam: -6.26 if lam == 0.0 else -7.5,
                  lambda x, lam: lam > 0.0 and x > -6.8)  # fmt: skip
    sc, _ = _stub_scanner(tmp_path, monkeypatch, data)
    r0 = sc.scan_second(0.0, None, None)
    assert r0["x_best"] == pytest.approx(-6.26, abs=0.01)
    fine0 = [p for p in sc.points.calls.values() if p["kind"] == "fine"]
    r1 = sc.scan_second(0.1309, (r0["x_best"], r0["width"]), sc._warm_of(r0["best_point"]))
    assert r1["x_best"] == pytest.approx(-7.5, abs=0.01) and "flag" not in r1
    assert sc.coarse_restarts == [0.1309]
    at = [p for p in sc.points.calls.values() if F.Points.key(p["p2"], 0)[0] == 0.1309]
    assert any(p["kind"] == "coarse" for p in at)
    assert any(p["kind"] == "fine" and p["no_contact"] for p in at)  # the guess's points stay
    assert all(p["commit"] for p in fine0)  # every call stamped with the code's commit


def test_no_contact_anywhere_in_the_box_is_recorded_and_the_scan_goes_on(tmp_path, monkeypatch):
    """Amendment 3, fourth addendum: a longitude with no contact anywhere in the box is recorded
    as zero-likelihood points, flagged "no contact in box", and the next longitude starts from
    the coarse scan."""
    data = _StubP(lambda lam: -6.26, lambda x, lam: abs(lam - 0.2618) < 1e-6)
    sc, _ = _stub_scanner(tmp_path, monkeypatch, data)
    monkeypatch.setattr(F.Task, "second_grid", lambda self: np.array([0.0, 0.2618, 0.3927]))
    out = sc.run()
    assert [r.get("flag") for r in out] == [None, "no contact in box", None]
    assert out[1]["best_point"] is None and out[2]["x_best"] == pytest.approx(-6.26, abs=0.01)
    pts = [p for p in sc.points.calls.values() if F.Points.key(p["p2"], 0)[0] == 0.2618]
    assert pts and all(p["no_contact"] for p in pts)
    assert F.ln_z_strength(pts, 1.0, F.box("P"), out[1]["h"]) == -math.inf
    assert sc.coarse_restarts == [0.2618]


def test_a_task_given_up_is_reported(tmp_path, monkeypatch):
    """Amendment 3, fourth addendum: a task failed twice is listed by every status check, from
    its marker or, for a driver started before markers, from the driver's log."""
    monkeypatch.setattr(F, "CACHE", tmp_path)
    tasks = [F.Task(3, 0, 0, "P", 0, 0), F.Task(1, 0, 0, "J"), F.Task(2, 0, 0, "J")]
    for t in tasks:
        t.dir.mkdir(parents=True)
    (tmp_path / "driver.log").write_text(
        "2026-11-16 06:07:20 start t3_s0_l0_P_d0_b0\n"
        "2026-11-16 06:30:00 FAILED t3_s0_l0_P_d0_b0 (exit 1, attempt 1)\n"
        "2026-11-16 06:30:00 start t3_s0_l0_P_d0_b0\n"
        "2026-11-16 06:30:20 FAILED t3_s0_l0_P_d0_b0 (exit 1, attempt 2)\n"
        "2026-11-16 07:00:00 start t1_s0_l0_J\n"
        "2026-11-16 07:10:00 FAILED t1_s0_l0_J (exit 1, attempt 1)\n"
        "2026-11-16 07:10:00 start t1_s0_l0_J\n"
    )
    assert F.failed_tasks(tasks) == ["t3_s0_l0_P_d0_b0"]
    (tasks[2].dir / "log.txt").write_text("Traceback\nValueError: x\n")
    F._gave_up(tasks[2], 1, 2)
    assert F.failed_tasks(tasks) == ["t2_s0_l0_J", "t3_s0_l0_P_d0_b0"]


def test_no_knot_after_the_coarse_scan_is_integrated_linearly_and_counted_apart(
    tmp_path, monkeypatch
):
    """Amendment 3, fourth addendum (clarification): contacts only at the coarse scan's
    strengths, so the fine points from its best have no knot. The value is flagged, counted
    apart from the coarse restarts (its integral is ln L linear between the coarse points, not
    the spline), and finite."""
    on_grid = lambda x: abs(x * 10.0 - round(x * 10.0)) < 1e-6  # noqa: E731
    data = _StubP(lambda lam: -6.26 if lam == 0.0 else -7.53,
                  lambda x, lam: lam > 0.0 and not on_grid(x))  # fmt: skip
    sc, _ = _stub_scanner(tmp_path, monkeypatch, data)
    monkeypatch.setattr(F.Task, "second_grid", lambda self: np.array([0.0, 0.1309]))
    out = sc.run()
    assert out[0].get("flag") is None and out[1]["flag"] == F.NO_KNOT
    assert sc.coarse_restarts == [0.1309]
    pts = [p for p in sc.points.calls.values() if F.Points.key(p["p2"], 0)[0] == 0.1309]
    assert not F.knots(pts) and any(p["kind"] == "coarse" and not p["no_contact"] for p in pts)
    assert math.isfinite(F.ln_z_strength(pts, 1.0, F.box("P"), out[1]["h"]))


def test_the_sensitivity_drops_linear_values_only():
    """The gate reading's line with and without the values integrated linearly: those rows are
    dropped; a task whose only value is such is kept whole and listed; none, no change."""
    rows = [{"p2": 0.0}, {"p2": 0.1309, "flag": F.NO_KNOT}, {"p2": 0.2618, "flag": F.NO_CONTACT}]
    res = {
        "P": {"scans": rows, "no_knot_linear": [0.1309]},
        "J": {"scans": [{"p2": 0.0, "flag": F.NO_KNOT}], "no_knot_linear": [0.0]},
        "R": {"scans": [{"p2": 0.1}], "no_knot_linear": []},
    }
    out, kept = F.drop_no_knot(res)
    assert [s["p2"] for s in out["P"]["scans"]] == [0.0, 0.2618]
    assert kept == ["J"] and out["J"] is res["J"] and out["R"] is res["R"]
    assert F.no_knot_sensitivity({"R": res["R"]}, ())["n_values"] == 0
