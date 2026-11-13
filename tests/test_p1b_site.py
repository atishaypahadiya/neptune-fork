"""The companion page's data (site/data/p1b_*.json) and the texts that quote Phase 1b numbers.

Every value the page shows is a copy of a named field of a committed JSON in docs/notes/
(docs/DESIGN.md, "Companion page"); these tests check each copy against its source, the track
and the recordings against their files, and the numbers of the README and the write-up note
against the same values.
"""

import csv
import json
import re

import pytest

from fork.export import p1b_site as S

SITE = S.SITE_DATA
FILES = ("p1b_adams.json", "p1b_challis.json", "p1b_page.json")


def _load(name: str) -> dict:
    return json.loads((SITE / name).read_text(encoding="utf-8"))


def test_committed_site_files_are_what_the_exporter_builds():
    built = S.build()
    assert set(built) == set(FILES)
    for name in FILES:
        assert _load(name) == json.loads(json.dumps(built[name], ensure_ascii=False)), name


@pytest.mark.parametrize("name", FILES[:2])
def test_every_value_is_its_source_field(name):
    d = _load(name)
    assert set(d["sources"]) == set(d["values"])
    cache: dict = {}
    for key, source in d["sources"].items():
        assert source[0].startswith("docs/notes/P1B_"), key
        assert d["values"][key] == S.resolve(source, cache), key


def test_track_and_recordings_match_their_files():
    d = _load("p1b_challis.json")
    pre = json.loads((S.ROOT / S.PRERUN).read_text())
    assert d["track"] == S._track()
    assert len(d["track"]) >= 60 and d["track"][0]["date"] == "1846-07-29"
    assert any(t["date"] == d["map_night"] for t in d["track"])
    assert pre  # the pre-run file the track is read from exists and parses
    with (S.ROOT / S.RECORDINGS).open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    counted = sum(int(r["nights"]) for r in rows if r["counted"] == "yes")
    assert counted == sum(r["nights"] for r in d["recordings"] if r["counted"] == "yes") == 8


def test_cited_constants_carry_a_citation():
    for key, c in _load("p1b_challis.json")["constants"].items():
        if key != "galle_night":
            assert c.get("citation"), key


def test_page_file_hashes_are_current():
    page = _load("p1b_page.json")
    for name, sha in page["built_from_sha256"].items():
        assert S._sha(name) == sha, name
    for note in page["notes"]:
        assert (S.ROOT / "docs" / "notes" / note).exists(), note


def test_the_copy_states_the_deviations():
    copy = (S.ROOT / "site" / "js" / "copy_p1b.js").read_text(encoding="utf-8")
    copy = re.sub(r"` \+\s+`", "", copy)  # sentences run over several template lines
    for phrase in ("post hoc", "by construction", "prior-dependent", "not halved",
                   "Outside the pre-registered reading", "depends on the method",
                   "an interpretation and not a reading",
                   "Post hoc, and defined after the cases were seen"):  # fmt: skip
        assert phrase in copy, phrase


def test_main_page_is_not_renumbered():
    html = (S.ROOT / "site" / "index.html").read_text(encoding="utf-8")
    assert re.findall(r"§ (\d)</p>", html) == [str(k) for k in range(1, 9)]
    assert html.count('data-v="companion_lead"') == 1


def _quoted() -> dict[str, str]:
    """Phase 1b numbers as the README and the write-up print them, from the site values."""
    a, c = _load("p1b_adams.json")["values"], _load("p1b_challis.json")["values"]
    nep = a["neptune_same_model"]
    return {
        "fit short, first hypothesis": f"{abs(a['H1_fit_true_lon_median'] - nep):.1f}°",
        "fit short, second hypothesis": f"{abs(a['H2_fit_true_lon_median'] - nep):.1f}°",
        "clean, his least squares": f"{abs(a['o_clean_least_squares']):.1f}°",
        "clean, project model on his places": f"{abs(a['o_clean_project_on_his_places']):.1f}°",
        "clean, project fit": f"{abs(a['o_clean_project_fit']):.1f}°",
        "post hoc, first hypothesis": f"{a['H1_post_hoc_P_a']:.3f}",
        "comparison carried through": f"{c['C1a_mean']:.2f} ({c['C1a_q05']:.2f} to "
        f"{c['C1a_q95']:.2f})",
        "S4": f"{c['C3_S4_mean']:.2f}",
        "S4 prompt": f"{c['C3_S4_prompt_mean']:.2f}",
        "Le Verrier June, prompt": f"{c['C4a_prompt_mean']:.2f}",
        "accepted": f"{c['n_accepted']:,}",
    }  # fmt: skip


@pytest.mark.parametrize("path", ["README.md", "docs/notes/P1B_writeup.md"])
def test_quoted_phase_1b_numbers(path):
    text = (S.ROOT / path).read_text(encoding="utf-8")
    for what, s in _quoted().items():
        assert s in text, f"{path}: {what} ({s})"


def test_section_3_5_matches_the_note():
    """The S8 sentence's premises, the ranges as ruled, and the manuscript line."""
    d = _load("p1b_adams.json")
    a = d["values"]
    assert a["bv_H1_d1c_row"] != a["bv_H1_chain_row"]
    assert not a["bv_s8a_H1_one_row"] and a["bv_s8b_H1_one_row"] and a["bv_s8c_H1_one_row"]
    got = [f"{a[f'bv_s8_bc_{k}']:+.1f}" for k in ("Ob_min", "Ob_max", "R_min", "R_max")]
    assert got == ["+8.6", "+12.0", "+0.5", "+3.1"]
    assert a["bv_H1_d1c_R"] + a["bv_H1_d1c_Ob"] == pytest.approx(a["bv_H1_d1c_shift"])
    ms = d["manuscript"]
    assert ms["source"] == "St John's College, Cambridge, MS W16 (arranged by R. A. Sampson, 1901)"
    assert (ms["n_values_agree"], ms["n_values_contradict"]) == (11, 0)
    assert "Maria Giovanna De Simone" in _load("p1b_page.json")["acknowledgements"]


def test_no_manuscript_image_on_the_site():
    names = [p.name.lower() for p in (S.ROOT / "site").rglob("*") if p.is_file()]
    assert not [n for n in names if "w16" in n or "stjohn" in n or n.startswith("168")]
