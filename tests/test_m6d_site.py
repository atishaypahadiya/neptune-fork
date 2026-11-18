"""The main page's data on the shape model (site/data/m6d.json): every value is a copy of a
named field of docs/notes/M6d.json, and the committed file is what the exporter builds."""

import json

import pytest

from fork.export import m6d_site as S


def _site() -> dict:
    return json.loads((S.SITE_DATA / "m6d.json").read_text(encoding="utf-8"))


def test_the_committed_file_is_what_the_exporter_builds():
    assert _site() == json.loads(json.dumps(S.build(), ensure_ascii=False))


def test_every_value_is_its_source_field():
    d = _site()
    assert set(d["sources"]) == set(d["values"])
    cache: dict = {}
    for key, source in d["sources"].items():
        assert source[0] == S.NOTE, key
        assert d["values"][key] == S.resolve(source, cache), key


def test_the_page_states_what_the_readings_say():
    # § 4.5 reads the primary fit (38.4 AU) and states the secondary beside it; both fits were
    # read, neither put a draw within Le Verrier's miss, and the bands are as registered.
    v = _site()["values"]
    assert v["shape38.read"] is True and v["shape36.read"] is True
    assert v["shape38.held_a_au"] == 38.4 and v["shape36.held_a_au"] == 36.154
    assert v["shape38.ln_b_wording"] == "a strong preference for free distance"
    assert v["shape36.ln_b_wording"] == "a weak preference for free distance"
    assert v["shape38.p_shape"] == 0.0 and v["shape36.p_shape"] == 0.0
    assert v["shape38.n_draws"] == 6000


def test_the_window_is_the_one_the_plan_fixed():
    w = _site()["window"]
    plan = (S.ROOT / "docs" / "notes" / "M6d.md").read_text(encoding="utf-8")
    question = plan.split("## 1 Question", 1)[1].split("## 2", 1)[0]
    assert f"P(|λ − 325.873°| ≤ {w['deg']}°)" in question
    assert f"actual miss of {w['arcmin']}′" in question
    assert w["arcmin"] / 60 == pytest.approx(w["deg"], abs=0.005)


def test_the_held_values_are_described_in_the_plans_words():
    h = _site()["held"]
    plan = (S.ROOT / "docs" / "notes" / "M6d.md").read_text(encoding="utf-8")
    held = plan.split("- **Held value, two runs:**", 1)[1].split("- **Data,", 1)[0]
    bullet = " ".join(held.split())
    assert "a = 38.4 AU.** " in bullet and f"({h['shape38']})" in bullet
    assert f"a = 36.154 AU,** {h['shape36']}." in bullet
