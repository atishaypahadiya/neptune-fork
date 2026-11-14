"""The main page's data on the S1 cell-0 rerun (site/data/s1_cell0.json): every value is a copy
of a named field of docs/notes/S1_cell0.json, and the committed file is what the exporter
builds."""

import json

from fork.export import s1_site as S


def _site() -> dict:
    return json.loads((S.SITE_DATA / "s1_cell0.json").read_text(encoding="utf-8"))


def test_the_committed_file_is_what_the_exporter_builds():
    assert _site() == json.loads(json.dumps(S.build(), ensure_ascii=False))


def test_every_value_is_its_source_field():
    d = _site()
    assert set(d["sources"]) == set(d["values"])
    cache: dict = {}
    for key, source in d["sources"].items():
        assert source[0] == S.NOTE, key
        assert d["values"][key] == S.resolve(source, cache), key


def test_the_clause_the_page_prints_matches_the_reading():
    # The page prints its clause while the rerun is unresolved; the note says the 95 % label
    # stays provisional for both reasons in that case.
    v = _site()["values"]
    assert v["resolved"] is False
    assert v["label_95"].startswith("provisional")
    assert v["a_max_au"] == 25.0
