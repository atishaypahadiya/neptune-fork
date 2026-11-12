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
    """The first check's record: the reviewed sample, every cell, this table, no mismatch, and not
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
    assert author["orcid"] == "https://orcid.org/0009-0004-1264-1662"
    assert zen["creators"][0] == {"name": "Pahadiya, Atishay", "orcid": "0009-0004-1264-1662"}
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


def test_doi_is_recorded_consistently_and_never_fetched():
    """From v1.0.1 the package cites the concept DOI (all versions) with the version number: a
    version's own DOI exists only once Zenodo has the files, which are frozen by then. The same
    DOI is in CITATION.cff (doi + identifiers, with v1.0.0's version DOI beside it), the README's
    How to cite section and the Data Package id. It is NOT in .zenodo.json. Nothing here
    fetches it. No file of the release package carries a badge: the package is frozen with
    the published record (the project README, outside the package, has one)."""
    doi = "10.5281/zenodo.24044791"
    assert rel.DOI == doi and rel.DOI_V1_0_0 == "10.5281/zenodo.24044792"
    cff = yaml.safe_load((rel.RELEASE / "CITATION.cff").read_text(encoding="utf-8"))
    assert cff["doi"] == doi
    ids = [{k: i[k] for k in ("type", "value")} for i in cff["identifiers"]]
    assert {"type": "doi", "value": doi} in ids
    assert {"type": "doi", "value": rel.DOI_V1_0_0} in ids
    assert json.loads(PKG.read_text(encoding="utf-8"))["id"] == f"https://doi.org/{doi}"
    readme = (rel.RELEASE / "README.md").read_text(encoding="utf-8")
    how = readme.split("## How to cite", 1)[1].split("\n## ", 1)[0]
    assert f"https://doi.org/{doi}" in how and f"doi       = {{{doi}}}" in how
    for p in rel.RELEASE.iterdir():
        if p.suffix in ("", ".md", ".cff", ".json"):
            assert "badge" not in p.read_text(encoding="utf-8").lower(), p.name
    zen = (rel.RELEASE / ".zenodo.json").read_text(encoding="utf-8")
    assert "doi" not in json.loads(zen) and "2404479" not in zen
    assert "v1.0.1 rewords descriptive strings in spotcheck_results.json" in readme
    for p in rel.RELEASE.iterdir():
        if p.suffix in ("", ".md", ".cff", ".json"):
            assert "XXXXXXX" not in p.read_text(encoding="utf-8"), p.name


# MD5 of every file of version 1.0.1. PROVISIONAL until the Zenodo record of v1.0.1 is published:
# these are the MD5s of the files as built for the upload (2026-11-12); once the record exists
# they are compared with the MD5s it lists and this comment names the record. The package is
# frozen: a change to any of these files needs a new Zenodo version (and new values here),
# never an edit in place. Version 1.0.0 is Zenodo record 24044792 (published 2026-11-07).
PUBLISHED_MD5 = {
    ".zenodo.json": "fc450ae2d7119c6097fe03a5ef61d2aa",
    "CITATION.cff": "17c571e4551bb5fc375fa527a2cd707f",
    "DATA_DICTIONARY.md": "b7c65478a7b426711895bb0b414be7c4",
    "datapackage.json": "cd7b52fc1cf8af07529d164a8fb5a137",
    "leverrier_1846.csv": "39532d973859330c1818097d20054958",
    "leverrier_1846_misprints.csv": "712fc6c34df3cbb463bc908433fc2ea5",
    "leverrier_1846_pass1.csv": "45a261923cb0eaac392f4dc02cf568a2",
    "leverrier_1846_pass2.csv": "107b337807c9b55a3fbfeeef1cb5a03f",
    "leverrier_1846_resolution_log.csv": "b5104ab9f6e8b8d439f1f82191259ef1",
    "leverrier_1846_table77.csv": "a35258c41b7444e26e35daec711ada51",
    "LICENSE": "968a9a067b63eb6e76f8d076c00f31be",
    "MANIFEST.md": "aff5d209ba337b985ee010e9c347a655",
    "README.md": "be9604d649e875d528e4f3473e92f930",
    "spotcheck_human_results.json": "ddc60e3b7e3897244b329b10d7a925a8",
    "spotcheck_results.json": "840ebd9a98589eedf965fccc885afc6a",
    "spotcheck_rows.csv": "70e9f2c07a09796d7bda4b666955c9bc",
}


def test_release_files_match_the_published_zenodo_record():
    """Offline: every file in release/leverrier-table77/ is byte-identical to the published
    record, and no file has been added or removed."""
    import hashlib

    files = {p.name: p for p in rel.RELEASE.iterdir()}
    assert set(files) == set(PUBLISHED_MD5)
    changed = [n for n, p in files.items()
               if hashlib.md5(p.read_bytes()).hexdigest() != PUBLISHED_MD5[n]]  # fmt: skip
    assert not changed, f"published files changed (publish a new Zenodo version): {changed}"
