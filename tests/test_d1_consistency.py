"""The study cache's deviation 1 results against the scratch copy's (the conditional go, (c))."""

import json

from fork.phase2 import d1_consistency as C

LEVELS = ("1.0", "0.5", "2.0")
BOXES = ("prior", "r3_1")


def _row(p2, z, tol=0.0025, **extra):
    d1 = {"tolerance": dict.fromkeys(BOXES, tol)}
    return {"p2": p2, "ln_z": {lv: dict.fromkeys(BOXES, z) for lv in LEVELS},
            "ln_z_before_d1": {}, "chi2_min_calls": 10.0, "x_best": -7.0,
            "x_best_argmin": -7.0, "d1": d1, **extra}  # fmt: skip


def _cache(root, tasks, refine=None):
    for tid, rows in tasks.items():
        d = root / "scans" / tid
        d.mkdir(parents=True)
        (d / "result.json").write_text(json.dumps({"scans": rows}))
        for (t, p2), xs in (refine or {}).items():
            if t != tid:
                continue
            (d / "refine").mkdir(exist_ok=True)
            lines = "".join(json.dumps({"x": x, "placed_by": "refine"}) + "\n" for x in xs)
            (d / "refine" / f"{p2:.10f}.jsonl").write_text(lines)


def _run(tmp_path, real_rows, scratch_rows, real_refine=None, scratch_refine=None):
    _cache(tmp_path / "real", real_rows, real_refine)
    _cache(tmp_path / "scratch", scratch_rows, scratch_refine)
    return C.compare(tmp_path / "real", tmp_path / "scratch")


def test_identical_caches_pass(tmp_path):
    rows = {"t1_s0_l0_P_d0_b0": [_row(0.5, -10.0)]}
    r = _run(tmp_path, rows, rows)
    assert r["pass"] and r["compared_rows"] == 1 and not r["differences"]


def test_the_read_back_rows_may_differ(tmp_path):
    r = _run(tmp_path, {"t2_s0_l0_J": [_row(0.0, -9.0)]}, {"t2_s0_l0_J": [_row(0.0, -10.0)]})
    assert r["pass"] and len(r["readback_rows"]) == 1


def test_a_difference_outside_the_allowed_groups_blocks(tmp_path):
    r = _run(tmp_path, {"t1_s0_l0_P_d0_b0": [_row(0.5, -10.001)]},
             {"t1_s0_l0_P_d0_b0": [_row(0.5, -10.0)]},
             real_refine={("t1_s0_l0_P_d0_b0", 0.5): [-7.1]})  # fmt: skip
    assert not r["pass"] and len(r["differences"]) == 1


def test_a_group_row_within_tau_and_traced_to_new_fits_is_allowed(tmp_path):
    tid = "t3_s0_l0_P_d1_b0"
    r = _run(tmp_path, {tid: [_row(0.5, -10.001)]}, {tid: [_row(0.5, -10.0)]},
             real_refine={(tid, 0.5): [-7.0, -7.05]},
             scratch_refine={(tid, 0.5): [-7.0]})  # fmt: skip
    assert r["pass"] and len(r["group_rows"]) == 1
    g = r["group_rows"][0]
    assert g["new_fits"] == [-7.05] and g["largest_over_tau"] < 1.0


def test_a_group_row_beyond_its_scratch_tolerance_blocks(tmp_path):
    tid = "t4_s0_l0_P_d2_b0"
    r = _run(tmp_path, {tid: [_row(0.5, -10.003)]}, {tid: [_row(0.5, -10.0)]},
             real_refine={(tid, 0.5): [-7.05]})  # fmt: skip
    assert not r["pass"] and r["differences"][0]["reason"] == "beyond the scratch tolerance"


def test_a_group_row_with_no_new_fit_blocks(tmp_path):
    tid = "t4_s0_l0_P_d2_b0"
    r = _run(tmp_path, {tid: [_row(0.5, -10.001)]}, {tid: [_row(0.5, -10.0)]},
             real_refine={(tid, 0.5): [-7.0]}, scratch_refine={(tid, 0.5): [-7.0]})  # fmt: skip
    assert not r["pass"] and r["differences"][0]["reason"] == "no fit added in this run"


def test_half_tasks_are_not_in_the_planet_groups(tmp_path):
    tid = "t3_s0_l0_P_d1_b0_half"
    r = _run(tmp_path, {tid: [_row(0.5, -10.001)]}, {tid: [_row(0.5, -10.0)]},
             real_refine={(tid, 0.5): [-7.05]})  # fmt: skip
    assert not r["pass"]


def test_guard_only_rows_take_the_largest_tolerance(tmp_path):
    tid = "t3_s0_l0_P_d2_b0"
    real = _row(0.5, -10.05)
    scratch = _row(0.5, -10.0)
    scratch["d1"]["tolerance"] = None
    r = _run(tmp_path, {tid: [real]}, {tid: [scratch]}, real_refine={(tid, 0.5): [-7.05]})
    assert r["pass"] and r["group_rows"][0]["tau"]["prior"] == C.TAU_MAX


def test_a_missing_row_blocks(tmp_path):
    r = _run(tmp_path, {"t1_s0_l0_R": []}, {"t1_s0_l0_R": [_row(0.14, -5.0)]})
    assert not r["pass"] and r["missing"]
