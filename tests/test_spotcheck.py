"""M6a spot-check sample: fixed seed, stratified, reproducible (fast; no scan needed)."""

import csv
import hashlib
import json

import numpy as np
import pytest

from fork.data import spotcheck as sc
from fork.data import transcription as tr


def test_sample_meets_the_strata_and_includes_every_logged_row():
    """Every row logged at review time (resolution log + the three original misprint
    entries) is in the sample."""
    table = tr.read_pass(tr.RECONCILED_FILE)
    picks = sc.read_selection(sc.SELECTION_FILE)
    forced = {(p, ln) for p, ln, _ in tr.read_resolutions(tr.RESOLUTION_FILE)}
    forced |= {("129", "10"), ("135", "3")}  # the misprint table at review time
    keys = [(p.page, p.line) for p in picks]
    assert len(keys) == sc.N_ROWS == 28 and len(set(keys)) == 28
    assert set(forced) <= set(keys)
    ancient = {(r["page"], r["line"]) for r in table[: [r["year"] for r in table].index("1781")]}
    assert len(set(keys) & ancient) >= sc.MIN_ANCIENT
    assert {p for p, _ in keys} == {r["page"] for r in table}


def test_sample_is_fixed_by_the_seed_and_matches_the_committed_selection():
    """The committed selection is the reviewed sample: drawing again with the seed and the
    forced rows it records gives it back exactly. (Misprint entries added after the review,
    e.g. the 1804 and 1807 times, are deliberately not in it.)"""
    table = tr.read_pass(tr.RECONCILED_FILE)
    committed_picks = sc.read_selection(sc.SELECTION_FILE)
    forced = sc.forced_in_selection(committed_picks)
    assert len(forced) == 16
    a = sc.select_rows(table, forced, np.random.default_rng(sc.SPOTCHECK_SEED))
    b = sc.select_rows(table, forced, np.random.default_rng(sc.SPOTCHECK_SEED))
    assert a == b == committed_picks
    with sc.SELECTION_FILE.open(newline="", encoding="utf-8") as fh:
        committed = [
            (r["page"], r["line"], r["reason"], int(r["seed"])) for r in csv.DictReader(fh)
        ]
    assert committed == [(p.page, p.line, p.reason, sc.SPOTCHECK_SEED) for p in a]


def test_all_match_is_disabled_for_logged_rows_only():
    table = tr.read_pass(tr.RECONCILED_FILE)
    picks = sc.read_selection(sc.SELECTION_FILE)
    forced = sc.forced_in_selection(picks)
    flagged = [p for p in picks if sc.is_flagged(p)]
    assert {(p.page, p.line) for p in flagged} == set(forced) and len(flagged) == 16
    page = sc.render_html(table, picks, {(p.page, p.line): "#" for p in picks}, "0" * 64)
    assert page.count('<button class="all" type="button" disabled>') == 16
    assert page.count('<button class="all" type="button">All match') == 12
    assert page.count('data-flagged="1"') == 16


def test_human_page_is_separate_from_the_first_review():
    """The human page stores ticks under its own key (so nothing from the automated review is
    pre-ticked), asks for a name, and exports a separate file marked reviewer = human."""
    table = tr.read_pass(tr.RECONCILED_FILE)
    picks = sc.read_selection(sc.SELECTION_FILE)
    links = {(p.page, p.line): "#" for p in picks}
    first = sc.render_html(table, picks, links, "0" * 64)
    human = sc.render_html(table, picks, links, "0" * 64, human=True)
    assert 'data-reviewer="human"' in human and "data-reviewer" not in first
    assert 'data-export="spotcheck_human_results.json"' in human
    assert 'id="reviewer-name"' in human and 'id="reviewer-name"' not in first
    assert '(REVIEWER ? REVIEWER + ":" : "")' in sc.SCRIPT  # the storage key
    assert "out.reviewer = REVIEWER" in sc.SCRIPT


def _first_record():
    return json.loads(sc.RESULTS_FILE.read_text(encoding="utf-8"))


def test_check_results_flags_mismatches_and_rejects_incomplete_records():
    table = tr.read_pass(tr.RECONCILED_FILE)
    picks = sc.read_selection(sc.SELECTION_FILE)
    sha = hashlib.sha256(tr.RECONCILED_FILE.read_bytes()).hexdigest()
    rec = _first_record()
    assert sc.check_results(rec, table, picks, sha, human=False) == []
    human = {**_first_record(), "reviewer": "human", "reviewer_name": "A. Reviewer"}
    human["rows"][0]["cells"]["day"]["matches"] = False
    assert sc.check_results(human, table, picks, sha, human=True) == ["p129 l4 day: '10'"]
    with pytest.raises(ValueError, match="reviewer = human"):
        sc.check_results({**human, "reviewer_name": " "}, table, picks, sha, human=True)
    with pytest.raises(ValueError, match="not the current"):
        sc.check_results(rec, table, picks, "f" * 64, human=False)
    with pytest.raises(ValueError, match="reviewed sample"):
        sc.check_results({**rec, "rows": rec["rows"][1:]}, table, picks, sha, human=False)
    wrong = _first_record()
    wrong["rows"][0]["cells"]["day"]["value"] = "11"
    with pytest.raises(ValueError, match="table has"):
        sc.check_results(wrong, table, picks, sha, human=False)
