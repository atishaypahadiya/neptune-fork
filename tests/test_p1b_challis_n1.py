"""C2 under the named input N1 itself (docs/notes/P1B_challis_detection.md, the section of
12 November 2026): the committed result, its checks, and the step 4 JSON left as it was."""

import hashlib
import json
from pathlib import Path

import pytest

from fork.export import p1b_challis_n1 as N

ROOT = Path(__file__).resolve().parents[1]
R = json.loads((ROOT / "docs" / "notes" / "P1B_challis_detection_n1.json").read_text())


def test_n1_limits_are_the_sheets():
    lo, hi = N.N1["ra_limits_h"]
    assert lo == pytest.approx(21 + 56 / 60) and hi == pytest.approx(23 + 4 / 60)
    assert N.N1["dec_limits_deg"] == [-15.0, 15.0] and N.N1["equinox"] == 1800
    assert R["n1"] == N.N1


def test_the_run_is_on_the_step_4_draws_and_equals_the_first_fallback():
    c = R["checks"]
    assert c["n1_chart_test_equals_first_fallback_on_every_night"] is True
    assert c["acceptance_equals_step_4"] is True
    assert c["max_abs_A0_minus_step_4_A0"] == pytest.approx(0.0, abs=1e-12)
    assert c["max_abs_n1_minus_step_4_first_fallback"] == pytest.approx(0.0, abs=1e-12)
    assert R["n_accepted"] == 100_000


def test_reading_and_the_untouched_step_4_json():
    s = R["case"]["held"]["primary"]
    assert s["mean"] == pytest.approx(0.0, abs=1e-12) and s["q95"] == pytest.approx(0.0, abs=1e-12)
    assert R["case"]["reading_row"] == "would very probably not have"
    step4 = ROOT / "docs" / "notes" / "P1B_challis_detection.json"
    assert hashlib.sha256(step4.read_bytes()).hexdigest() == R["step_4_json_sha256"]


def test_the_planet_is_west_of_the_edge_on_the_two_recorded_nights():
    ra = R["neptune_ra_1800_h"]
    lo = N.N1["ra_limits_h"][0]
    assert ra["1846-08-04"] < lo and ra["1846-08-12"] < lo
    assert ra["1846-07-29"] > lo and ra["1846-07-30"] > lo
