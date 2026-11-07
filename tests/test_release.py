"""Release package release/leverrier-table77/: Frictionless validation, counts, provenance."""

import csv
import json
import re
import shutil
from collections import Counter

import pytest
import yaml
from frictionless import Package, validate

from fork.data import spotcheck as sc
from fork.data import transcription as tr
from fork.export import release as rel

PKG = rel.RELEASE / "datapackage.json"
DOCS = {"README.md", "MANIFEST.md", "DATA_DICTIONARY.md", "LICENSE", "CITATION.cff", ".zenodo.json",
        "datapackage.json"}  # fmt: skip
PENDING = "human spot check is pending"


def test_release_files_are_byte_copies_of_the_committed_project_files():
    """The release is a snapshot: a change to a source file must come with a new release."""
    for name, src in rel.sources().items():
        assert (rel.RELEASE / name).read_bytes() == src.read_bytes(), name


def test_release_contains_exactly_the_built_files():
    """No stale or stray file (e.g. a human record copied in by hand, or left behind)."""
    assert {p.name for p in rel.RELEASE.iterdir()} == set(rel.sources()) | DOCS


def test_datapackage_json_is_the_builder_output():
    expected = json.dumps(rel.descriptor(), indent=2, ensure_ascii=False) + "\n"
    assert PKG.read_text(encoding="utf-8") == expected


def test_package_validates():
    """Schemas (types, patterns, enums), primary and foreign keys, and SHA-256/byte counts."""
    report = validate(str(PKG))
    errors = [e.message for e in report.errors]
    errors += [f"{t.name}: {e.message}" for t in report.tasks for e in t.errors]
    assert report.valid, errors[:10]
    assert {t.name for t in report.tasks} == {r["name"] for r in rel.resources()}


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


def _record(name):
    return json.loads((rel.RELEASE / name).read_text(encoding="utf-8"))


def test_spotcheck_records_are_complete_for_this_table():
    """The automated record: the reviewed sample, every cell, this table, no mismatch, and not
    marked human. The human record, once present, must pass the same check as human."""
    table = tr.read_pass(rel.RELEASE / "leverrier_1846_table77.csv")
    picks = sc.read_selection(rel.RELEASE / "spotcheck_rows.csv")
    sha = rel.sha256(rel.RELEASE / "leverrier_1846_table77.csv")
    rec = _record("spotcheck_results.json")
    assert sc.check_results(rec, table, picks, sha, human=False) == []
    assert sum(len(r["cells"]) for r in rec["rows"]) == 336
    assert "not an independent human" in rec["reviewer"]
    with pytest.raises(ValueError, match="reviewer = human"):
        sc.check_results(rec, table, picks, sha, human=True)
    assert rel.human_checked()  # the release requires the human record
    human = _record(sc.HUMAN_RESULTS_FILE.name)
    assert sc.check_results(human, table, picks, sha, human=True) == []
    assert human["reviewer_name"] == "Atishay Pahadiya" and human["all_match_rows"] == []


def _flat(path):
    """Lowercase text with all whitespace runs collapsed, so wrapped lines still match."""
    return " ".join(path.read_text(encoding="utf-8").split()).lower()


def test_verified_wording_and_no_pending_wording():
    """The human spot check is done: every description of the dataset states the verified
    result, and none still says it is pending. (This guard used to require "pending" while
    no human record existed; it now requires the record, see the test above.)"""
    files = ("README.md", "DATA_DICTIONARY.md", "CITATION.cff", ".zenodo.json",
             "datapackage.json", "../CHECKLIST_leverrier-table77.md")  # fmt: skip
    assert rel.VERIFIED == (
        "spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies)"
    )
    for f in files:
        text = _flat(rel.RELEASE / f)
        assert rel.VERIFIED.lower() in text, f"{f}: verified wording missing"
        if f != "../CHECKLIST_leverrier-table77.md":  # the checklist quotes the old wording
            assert PENDING not in text, f"{f}: still says pending"


def test_manifest_records_the_spotcheck_hashes():
    manifest = (rel.RELEASE / "MANIFEST.md").read_text(encoding="utf-8")
    for name in ("spotcheck_results.json", sc.HUMAN_RESULTS_FILE.name):
        assert f"`{rel.sha256(rel.RELEASE / name)}`" in manifest, name


def test_licence_is_cc0():
    text = (rel.RELEASE / "LICENSE").read_text(encoding="utf-8")
    assert "CC0 1.0 Universal" in text and "Statement of Purpose" in text
    assert "Attribution 4.0" not in text and "DRAFT" not in text
    readme = (rel.RELEASE / "README.md").read_text(encoding="utf-8")
    assert "CC0 1.0 Universal" in readme and "Please cite" in readme


def test_data_dictionary_documents_every_column():
    """Each field of each tabular resource has a row `| \\`name\\` |` in its file's section."""
    text = (rel.RELEASE / "DATA_DICTIONARY.md").read_text(encoding="utf-8")
    sections = re.split(r"^## ", text, flags=re.M)
    for r in rel.resources():
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
    assert zen["license"] == rel.LICENSE["name"].lower() == "cc0-1.0"
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


def test_transit_check_covers_every_row_as_the_readme_says():
    """README: every row's date and mean time fit the meridian-transit identity within 3.6 s,
    except the logged 1768 Dec 30 time (+9.5 s); a misread day or minute would miss by >= 60 s."""
    with (rel.RELEASE / "leverrier_1846.csv").open(newline="", encoding="utf-8") as fh:
        ra = [r for r in csv.DictReader(fh) if r["quantity"] == "ra"]
    assert len(ra) == 279
    resid = {r["obs_id"]: abs(float(r["transit_resid_s"])) for r in ra}
    assert resid.pop("LV1846-T77-p129-l12") == pytest.approx(9.5, abs=0.05)
    assert max(resid.values()) <= 3.6
