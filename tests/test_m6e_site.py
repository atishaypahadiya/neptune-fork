"""The main page's data for § 4.4 (site/data/m6e.json): every value is a copy of a named
field of docs/notes/M6e.json, and the committed file is what the exporter builds."""

import ast
import json
from pathlib import Path

from fork.export import m6e_site as S


def _site() -> dict:
    return json.loads((S.SITE_DATA / "m6e.json").read_text(encoding="utf-8"))


def test_the_committed_file_is_what_the_exporter_builds():
    assert _site() == json.loads(json.dumps(S.build(), ensure_ascii=False))


def test_every_value_is_its_source_field():
    d = _site()
    assert set(d["sources"]) == set(d["values"])
    cache: dict = {}
    for key, source in d["sources"].items():
        assert source[0] == "docs/notes/M6e.json", key
        assert d["values"][key] == S.resolve(source, cache), key


def test_the_histogram_is_of_the_posterior_the_text_quotes():
    v = _site()["values"]
    edges, counts = v["histogram.edges_deg"], v["histogram.counts"]
    assert len(edges) == len(counts) + 1 and sum(counts) == v["n_draws"]
    assert edges[0] < v["headline.1795-05-08.lon.median"] < edges[-1]


def test_the_exporter_does_not_read_the_held_out_table():
    tree = ast.parse(Path(S.__file__).read_text(encoding="utf-8"))
    imports = [n for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))]
    names = {a.name for n in imports for a in n.names} | {
        n.module for n in imports if isinstance(n, ast.ImportFrom)
    }
    assert names == {"__future__", "annotations", "json", "pathlib", "Path"}
