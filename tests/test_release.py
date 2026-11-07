"""Release package release/leverrier-table77/: Frictionless validation, counts, provenance."""

import csv
import json
import re
import shutil
from collections import Counter

import pytest
import yaml
from frictionless import Package, validate

from fork.data import transcription as tr
from fork.export import release as rel

PKG = rel.RELEASE / "datapackage.json"


def test_release_files_are_byte_copies_of_the_committed_project_files():
    """The release is a snapshot: a change to a source file must come with a new release."""
    for name, src in rel.SOURCES.items():
        assert (rel.RELEASE / name).read_bytes() == src.read_bytes(), name


def test_datapackage_json_is_the_builder_output():
    expected = json.dumps(rel.descriptor(), indent=2, ensure_ascii=False) + "\n"
    assert PKG.read_text(encoding="utf-8") == expected


def test_package_validates():
    """Schemas (types, patterns, enums), primary and foreign keys, and SHA-256/byte counts."""
    report = validate(str(PKG))
    errors = [e.message for e in report.errors]
    errors += [f"{t.name}: {e.message}" for t in report.tasks for e in t.errors]
    assert report.valid, errors[:10]
    assert {t.name for t in report.tasks} == {r["name"] for r in rel.RESOURCES}


def test_table77_has_279_rows_with_the_printed_page_counts():
    rows = Package(str(PKG)).get_resource("table77").read_rows()
    assert len(rows) == 279
    assert Counter(r["page"] for r in rows) == rel.PAGE_COUNTS
    assert sum(rel.PAGE_COUNTS.values()) == 279
    for name in ("pass1", "pass2"):  # both passes have the same row structure
        p = Package(str(PKG)).get_resource(name).read_rows()
        assert [(r["page"], r["line"]) for r in p] == [(r["page"], r["line"]) for r in rows]
    years = [r["year"] for r in rows]
    assert years.index(1781) == 19  # 19 ancient rows (1690–1771), 260 modern


@pytest.mark.parametrize(
    "file, edit, expected_error",
    [
        # one digit changed: the hash no longer matches
        ("leverrier_1846_table77.csv", ("55 49 19.7", "55 49 19.8"), "hash-count"),
        # a malformed sexagesimal value
        ("leverrier_1846_table77.csv", ("55 49 19.7", "55.49 19.7"), "constraint-error"),
        # a misprint entry pointing at a row that does not exist
        ("leverrier_1846_misprints.csv", ("\n135,3,dec_obs", "\n135,99,dec_obs"),
         "foreign-key"),
    ],
)  # fmt: skip
def test_validation_catches_tampering(tmp_path, file, edit, expected_error):
    """Negative controls: the validator is not vacuous."""
    d = tmp_path / "pkg"
    shutil.copytree(rel.RELEASE, d)
    p = d / file
    text = p.read_text(encoding="utf-8")
    assert text.count(edit[0]) == 1
    p.write_text(text.replace(*edit), encoding="utf-8")
    if expected_error != "hash-count":  # isolate the schema/key check from the hash check
        desc = json.loads((d / "datapackage.json").read_text(encoding="utf-8"))
        for r in desc["resources"]:
            r.pop("hash"), r.pop("bytes")
        (d / "datapackage.json").write_text(json.dumps(desc), encoding="utf-8")
    report = validate(str(d / "datapackage.json"))
    types = {e.type for t in report.tasks for e in t.errors}
    assert not report.valid and expected_error in types, types


def test_spotcheck_record_is_for_this_table():
    rec = json.loads((rel.RELEASE / "spotcheck_results.json").read_text(encoding="utf-8"))
    assert rec["table77_sha256"] == rel.sha256(rel.RELEASE / "leverrier_1846_table77.csv")
    assert len(rec["rows"]) == 28
    assert sum(len(r["cells"]) for r in rec["rows"]) == 336
    assert "not an independent human" in rec["reviewer"]


def test_data_dictionary_documents_every_column():
    """Each field of each tabular resource has a row `| \\`name\\` |` in its file's section."""
    text = (rel.RELEASE / "DATA_DICTIONARY.md").read_text(encoding="utf-8")
    sections = re.split(r"^## ", text, flags=re.M)
    for r in rel.RESOURCES:
        if "schema" not in r:
            continue
        sec = next((s for s in sections if f"`{r['path']}`" in s.splitlines()[0]), None)
        assert sec is not None, f"no section headed with `{r['path']}`"
        with (rel.RELEASE / r["path"]).open(newline="", encoding="utf-8") as fh:
            header = next(csv.reader(fh))
        assert header == [f["name"] for f in r["schema"]["fields"]], r["path"]
        missing = [c for c in header if f"| `{c}` |" not in sec]
        assert not missing, f"{r['path']}: undocumented {missing}"


def test_citation_metadata_agree():
    """datapackage.json, CITATION.cff and .zenodo.json carry one title, version and licence."""
    cff = yaml.safe_load((rel.RELEASE / "CITATION.cff").read_text(encoding="utf-8"))
    zen = json.loads((rel.RELEASE / ".zenodo.json").read_text(encoding="utf-8"))
    assert cff["title"] == zen["title"] == rel.TITLE
    assert cff["version"] == zen["version"] == rel.VERSION
    assert cff["license"] == rel.LICENSE["name"]
    assert zen["license"] == rel.LICENSE["name"].lower()
    author = cff["authors"][0]
    assert (author["family-names"], author["given-names"]) == ("Pahadiya", "Atishay")
    assert zen["creators"][0]["name"] == "Pahadiya, Atishay"
    assert zen["upload_type"] == "dataset"
    assert set(zen["keywords"]) == set(cff["keywords"]) == set(rel.KEYWORDS)


def test_manifest_hashes_match_the_project_manifest():
    """The release MANIFEST repeats the source-scan hashes of data/historical/sources."""
    hexes = re.compile(r"`([0-9a-f]{64})`")
    project = set(hexes.findall((tr.HIST / "sources" / "MANIFEST.md").read_text("utf-8")))
    release = set(hexes.findall((rel.RELEASE / "MANIFEST.md").read_text("utf-8")))
    assert len(project) == 3 and project <= release
