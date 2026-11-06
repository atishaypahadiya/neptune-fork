"""M6a spot-check sample: fixed seed, stratified, reproducible (fast; no scan needed)."""

import csv

import numpy as np

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
