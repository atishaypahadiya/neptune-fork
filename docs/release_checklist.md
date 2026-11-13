# Release checklist: v1.0.0

**No tag before M6 is complete.** The version was `0.9.0` through M0–M5b and the M6 data and
observation model; `1.0.0` was set on 2026-11-10 with M6 complete (M6d, M6e and the variant
fits deferred, below) for the v1.0.0 tag. Versions follow
semver: the tag is `vX.Y.Z`, and `pyproject.toml` and `CITATION.cff` hold `X.Y.Z`.
`.zenodo.json` has no `version`: Zenodo takes it from the release tag (decided 2026-11-09).

## Decided (2026-11-07)

- No tag yet. Tag `v1.0.0` after M6 is complete.
- `0.9.0` until then, with semver everywhere and no "-phase1" suffix.
- Licenses: MIT for the code, and CC0 1.0 for `data/historical/`. The split is stated in
  `LICENSE`, `README.md`, `CITATION.cff` (`license: [MIT, CC0-1.0]`, CC0 also in the abstract)
  and `.zenodo.json` (`license: MIT`, one licence only; CC0 stated in `notes`).
  Zenodo's `license` field takes one value, so it says MIT, and the CC0 dedication is in
  `notes` (since 2026-11-09 the description covers only the method and the result).
- ORCID: 0009-0004-1264-1662 (https://orcid.org/0009-0004-1264-1662), filled in 2026-11-07: the URL in `CITATION.cff`, the bare iD
  in `.zenodo.json`.

## Before tagging: M6 complete

- [x] A6c passed and committed (`site/data/a6c_coverage.json`, with a note in
  `docs/notes/M6.md`). If it fails, stop and report, as the spec says. Don't tag.
- [x] M6c fits are done, with `site/data/m6_fits.json` and the figures the spec lists. (M6d and
  M6e moved to "Deferred to v1.1" below, 2026-11-09, and the M6c variant fits with them.)
  2026-11-10: S1 unresolved per outcome 3 (the A8d(S1) arbiter's pre-registered reading,
  `docs/notes/M8.md`); the published S1 verdict stands.
- [x] All A6 tests reported against their pre-registered thresholds (A6a–A6e).
  - A6a–A6d reported in `docs/notes/A6.md` (2026-11-09): A6a, A6c and A6d pass, and A6b's
    loader test passes. A6e is deferred with M6e.
  - A6b is report-only: Comparison against per-group precision not made: no precision figure was fixed before the data were seen. Medians reported.
- [x] README results and status updated to include the real-data results, keeping the
  synthetic and real numbers clearly apart.

## Deferred to v1.1 (decided 2026-11-09)

- **M6d, the Le Verrier shape-model comparison:** the Peirce probability is computed for every
  M6c fit (0.366 headline), but the ratio against the Le Verrier shape model and the Bayes
  factor between the two need a new pre-registered real-data nested run (hours of compute).
- **M6e, the held-out precoveries:** the Lalande 1795 and Galileo 1612–13 positions are not
  yet sourced or transcribed, and the predictive step is not built; A6e depends on it.
- **M6c variant fits** (`bouvard`, `adams`, the 1846 reduction track, the Le Verrier shape
  model; moved here 2026-11-09): not yet pre-registered; to be planned and pushed before
  launch; the `adams` variant opens the Adams–Challis study.
- **S1 cell 0 resolution:** new pre-registered plan with a mode-robust sampler. The A8d(S1)
  arbiter (2026-11-10) ended in its third outcome, "fail, agreeing with neither": the published
  verdict ("not at all") stands and the question is unresolved (`docs/notes/M8.md`).
- **Student-t profile in the semi-major axis:** M9's profile likelihood under S3's Student-t
  error model (ν = 4), to be pre-registered; the Gaussian profile is reported in
  `docs/notes/M9.md` (2026-11-10).
- **Cache-only fetcher for pre-registered runs:** future pre-registered runs use a cache-only
  fetcher that raises on miss. Today a run reads Horizons through the default fetcher, which
  goes live on a cache miss. The A8d(S1) arbiter (launched 2026-11-09) made no network calls,
  but only because every epoch it needs was cached.

## Repository state

- [x] Branch merged to `main`, working tree clean, and nothing in `data/cache/` is needed to
  read the release. 2026-11-10: all work is on `main`; every result the README and site read
  is under `site/data/` or `docs/notes/` (the M9 report, the arbiter summary, the fits).
- [x] `uv sync --locked` succeeds on a fresh clone. 2026-11-10: clone of b799fd1, `uv sync
  --locked` in 0.5 s (packages cached), `fork` 1.0.0 importable.
- [x] `uv run ruff check . && uv run ruff format --check .` passes. 2026-11-10: "All checks
  passed!", 100 files already formatted.
- [x] `uv run pytest -m "not slow and not network"` passes locally and on CI (the badge is
  green on `main`). 2026-11-09: a0ec860, CI run
  https://github.com/atishaypahadiya/neptune-fork/actions/runs/42740104130 (lint and test
  green, 322 passed) and locally 322 passed.
- [x] `uv run pytest -m slow` passes locally. Record the wall time. 2026-11-10: 22 passed in
  5 min 59 s wall (1344 s CPU), while the layout check ran alongside.
- [x] `uv run python examples/toy_fit.py` finishes in under 5 min on an idle machine. The time
  in `docs/reproduce.md` (4.6 min) was measured under A6c load, so re-time it and update.
  2026-11-10: 2 min 58 s on the idle machine; `docs/reproduce.md` updated.
- [x] Committed `site/data/*.json` matches a fresh export where tests check that
  (`test_export.py`, `test_ledger.py`). 2026-11-10: in the fast suite's 332 passed.
- [x] `uv run python -m fork.export.ledger` re-run, and the CPU-hour totals in `README.md` and
  `docs/reproduce.md` updated (they will include A6c and M6c). 2026-11-10: 1087 CPU-h (1004
  live, 83 estimated), through M9; both files say so.
- [x] The A6d static check holds: no M6a–M6d code path reads Neptune's real ephemeris or the
  M6e precovery file (`tests/test_a6d.py`, 2026-11-09; runtime half in `tests/test_m6c.py`).
- [x] Whitmore's *discoverneptune* is cited in the README. If any of its data are reused,
  check its license first. 2026-11-09: `README.md` line 165; none of its data are reused (it
  appears only as a citation, here and in `CITATION.cff`, `.zenodo.json` and
  `release/leverrier-table77/README.md`).

## Metadata

- [x] Version set to `1.0.0` in `pyproject.toml`, then `uv lock`, and the same value in
  `CITATION.cff` and `.zenodo.json`. 2026-11-10: `pyproject.toml` and `CITATION.cff` say
  `1.0.0`, `uv.lock` records `fork` 1.0.0; `.zenodo.json` carries no version (Zenodo takes it
  from the tag, decided 2026-11-09).
- [x] ORCID placeholder replaced in both files (2026-11-07):
  - `CITATION.cff` had `https://orcid.org/0000-0000-0000-0000`, now `https://orcid.org/0009-0004-1264-1662`;
  - `.zenodo.json` had `PLACEHOLDER-ORCID-replace-before-release`, now `0009-0004-1264-1662`.

  `grep -n "0000-0000-0000-0000\|PLACEHOLDER" CITATION.cff .zenodo.json` must print nothing.
- [x] `CITATION.cff`: add `date-released: YYYY-MM-DD` and validate with
  `uvx cffconvert --validate`. 2026-11-10: `date-released: "2026-11-10"`; "Citation metadata
  are valid according to schema version 1.2.0".
- [x] `.zenodo.json` is valid JSON and agrees with `CITATION.cff`. Zenodo uses `.zenodo.json`
  over `CITATION.cff` when both exist. 2026-11-10: parses; same title, creator, ORCID,
  affiliation, description text, dataset DOI and Whitmore reference; licence MIT with the CC0
  note, as decided.
- [x] Zenodo GitHub integration enabled for `atishaypahadiya/neptune-fork`. This has to happen
  **before** the GitHub release is published, or no DOI is minted. 2026-11-09: the repository
  has an active webhook to `zenodo.org/api/hooks/receivers/github` on `release` events
  (`gh api repos/atishaypahadiya/neptune-fork/hooks`); `docs/notes/M8.md`, Operational notes.

## Tag and release

- [x] `git tag -a v1.0.0 -m "v1.0.0: Phase 1 through M6"` 2026-11-10: on 1ac45a3, with CI (test
  and pages) green on it.
- [x] `git push origin v1.0.0` 2026-11-10: pushed over SSH.
- [x] GitHub release from the tag. 2026-11-10 22:36:57 UTC, not a draft:
  https://github.com/atishaypahadiya/neptune-fork/releases/tag/v1.0.0, with the notes below in
  the paper's voice (what was found, what is provisional or unresolved, what is deferred to
  v1.1, the dataset DOIs, the licence split, the ledger total). Release notes:
  - what is included;
  - the headline numbers, each marked synthetic or real:
    - A1: 2e-9″ against REBOUND;
    - A3a: the noise-free *point fit* at 0.07°;
    - the noise-free *posterior*: median −0.75σ from the truth;
    - A3b: 8/10;
    - A6c;
    - the M6c headline fit;
  - the ×1.8 evidence-error inflation;
  - the license split (MIT code, CC0 1.0 for `data/historical/`);
  - the compute ledger total.
- [x] Zenodo DOI minted. Add the DOI badge to `README.md` and `doi:` to `CITATION.cff`, in a
  follow-up commit. 2026-11-10: version DOI 10.5281/zenodo.24121060, concept DOI
  10.5281/zenodo.24121059 (https://zenodo.org/records/24121060), minted within 20 minutes of
  the release; badge and availability bullet in `README.md`, `doi:` and both identifiers in
  `CITATION.cff`, and the site's byline and Data availability read them from
  `site/data/release.json`.
- [x] The Pages site still deploys from `main`. 2026-11-10: the pages workflow deployed 1ac45a3
  and the live page carries the release byline ("code v1.0.0").

## v1.1.0 (2026-11-10)

Released on the user's explicit go-ahead ("Release v1.1.0"). What it adds to v1.0.0:

- **M10**, the pre-registered profile likelihood in the semi-major axis under Student-t errors
  (`docs/notes/M10.md`, results appended 2026-11-10): the v1.0.0 distance claim held only under
  the Gaussian error model. The abstract's distance sentences, § 5 (Fig. 7 with both error
  models, Fig. 9 with the leave-one-out profiles) and Limitations are updated; every number is
  read from `site/data/m10_profile.json`.
- The two pre-run changes M10 recorded: the Uranus-only box centred on the better of two
  starts (1c192e2) and the cache-only ephemeris fetcher, `FORK_HORIZONS_CACHE_ONLY=1`
  (41f7a12), which closes that item of the list above.
- Still deferred, now to v1.2: M6d, M6e, the M6c variant fits, the S1 cell-0 resolution
  (`site/data/m7_timeline.json` "deferred" says 1.2).
- Steps, all done 2026-11-10 (PST): version 1.1.0 in `pyproject.toml`, `uv.lock` and
  `CITATION.cff` (its `doi:` removed until Zenodo minted the new version DOI; the concept DOI
  stayed); `CHANGELOG.md` added; checks green (346 fast tests, ruff, prerender audit 0
  unmatched, layout check at 360–1440 px and print); CI (test and pages) green on ad5ada0; tag
  `v1.1.0` on ad5ada0; GitHub release published 2026-11-11 01:41:22 UTC (17:41 PST), not a
  draft, notes in the paper's voice; Zenodo minted version DOI 10.5281/zenodo.24122940 within
  10 minutes (https://zenodo.org/records/24122940; the record's publication date is the UTC
  day, 2026-11-11), and the concept DOI 10.5281/zenodo.24121059 resolves to that record;
  follow-up commit with the DOI in `CITATION.cff` (`doi:` and a v1.1.0 identifier beside
  v1.0.0's), `README.md` and `site/data/release.json`.

## v1.2.0 list (opened 2026-11-11; no release without the explicit word; closed below)

- Carried from v1.1: M6d, M6e, the M6c variant fits, the S1 cell-0 resolution.
- `docs/reproduce.md`: write sections for the M6c headline, S1–S3, R1/R2, M7, M8/M8b and the
  A8d(S1) arbiter from the notes and the ledger. No version of the file has had them; those
  runs appear only in its run-identifier table.
- At the site step: regenerate the site ledger (`site/data/ledger.json`) with the Phase 1b rows
  (`p1b_adams/*`, `p1b_adams_repro/*`, `p1b_adams_noise/draws`; 103 CPU-h, held in the local
  ledger until then).
- Release v1.2.0 so the concept DOI resolves to the files reworded in f7fecdc
  ("Describe transcription passes by method").
- With that release, the dataset tree: the same rewording in `src/fork/export/release.py`
  (the descriptor's description and the spot-check resource) was held back on 2026-11-11
  because `tests/test_release.py` pins `release/leverrier-table77/datapackage.json` to the
  builder's output, and `release/` is not edited outside a release. The dataset's own README,
  CITATION.cff, `.zenodo.json`, DATA_DICTIONARY, MANIFEST and CHECKLIST carry the same earlier
  wording; decide their text when the dataset version is cut.

## v1.2.0 (prepared 2026-11-12; not tagged until the explicit word)

Plan approved 2026-11-12 with four decisions: M6d, M6e, the S1 cell-0 resolution and the
remaining M6c variants go to v1.3; dataset v1.0.1 is cut with this release; the descriptive
strings of `spotcheck_results.json` are reworded with values unchanged; the Phase 1b tables are
in the code record only.

- [x] Dataset tree v1.0.1 (5988a01): documents, the builder and the spot-check record's strings
  reworded; `datapackage.json` regenerated; the package cites the concept DOI with the version
  number; `tests/test_release.py` pins the MD5s of the files as built, marked provisional.
- [x] `docs/reproduce.md`: the Phase 1 sections listed above.
- [x] Metadata: version 1.2.0 in `pyproject.toml`, `uv.lock` and `CITATION.cff` (its `doi:`
  removed until Zenodo mints the version DOI); `.zenodo.json`; `CHANGELOG.md`;
  `site/data/release.json`; the companion page's status; `site/data/m7_timeline.json`
  "deferred" says 1.3.
- [x] Checks: ruff, tests, prerender, number audit on both pages, layout in Chromium and WebKit.
- [x] Dataset v1.0.1 uploaded to Zenodo by hand as a new version of record 24044792:
  published 2026-11-12, https://zenodo.org/records/24158992, version DOI 10.5281/zenodo.24158992.
- [x] After the upload: the record's 16 MD5s equal the pins, and the test's comment names the
  record; the v1.0.1 version DOI is in the project README and `site/data/m7_dataset.json`.
- [x] Before the tag: the epoch of the held-distance fits' shortfall corrected to 1846.762
  (a262a35), and the prior-dependent figures marked in the changelog, README and release notes.
- [x] Tag `v1.2.0` on 5df9f42 (test and pages green), GitHub release published 2026-11-12
  11:14:57 UTC, not a draft; Zenodo minted version DOI 10.5281/zenodo.24159832
  (https://zenodo.org/records/24159832), and the concept DOI 10.5281/zenodo.24121059 resolves
  to that record.
- [x] Follow-up commit: that DOI in `CITATION.cff`, `README.md` and `site/data/release.json`.

## v1.2.1 (released on the author's word, 2026-11-13 UTC)

A patch. Contents in `CHANGELOG.md`: the chart case corrected, Krajnović's book, the Lille
link, public names and voice, the manuscript check, Bouvard's tables, the CI checks, the grid.

- [x] Metadata: version 1.2.1 in `pyproject.toml`, `uv.lock` and `CITATION.cff` (its `doi:`
  removed until Zenodo mints the version DOI); `site/data/release.json`; the companion
  page's status.
- [x] Checks: ruff, tests, prerender, number audit on both pages, layout in Chromium and
  WebKit, citation validation, one link pass.
- [x] Tag `v1.2.1` on 4729be9 (test and pages green), GitHub release published 2026-11-13
  03:35:34 UTC, not a draft; Zenodo minted version DOI 10.5281/zenodo.24183833
  (https://zenodo.org/records/24183833).
- [x] Follow-up commit: that DOI in `CITATION.cff`, `README.md` and `site/data/release.json`.
- [x] The v1.2.0 GitHub release notes corrected in place (title and five phrases).
- Not in this release: the dataset record is unchanged (its Lille handle is corrected by
  hand on Zenodo); the v1.2.0 release notes and Zenodo description are corrected in place.

## v1.3.0 (released on the author's word, 2026-11-13 UTC)

A minor release. Contents in `CHANGELOG.md`: the attribution of Adams's nearness to his
observations or to Bouvard's tables (pre-registered; "depends on the method", not
attributed), the companion page's § 3.5, the two findings on Bouvard's tables, the
manuscript check's spot check, the first Mercury data set and forces, the soupsieve
advisory. The list below was the v1.3 list; M6d, M6e, the M6c variants and S1 cell 0 were
not done in v1.3.0 and are carried to v1.4.

- [x] Metadata: version 1.3.0 in `pyproject.toml`, `uv.lock` and `CITATION.cff` (its `doi:`
  removed until Zenodo mints the version DOI); one sentence added to the abstract in
  `CITATION.cff` and `.zenodo.json`; `site/data/release.json`; the companion page's status;
  the deferred fits' target version (1.4) in `site/data/m7_timeline.json` and the README.
- [x] Checks: ruff, tests, prerender, number audit on both pages, layout in Chromium and
  WebKit, citation validation, one link pass.
- [ ] Tag `v1.3.0`, GitHub release, Zenodo version DOI.
- [ ] Follow-up commit: that DOI in `CITATION.cff`, `README.md` and `site/data/release.json`.

## v1.4 list (opened 2026-11-12 as the v1.3 list; no release without the explicit word)

- Carried from v1.3: M6d, M6e, the remaining M6c variants, the S1 cell-0 resolution.
- A human spot check of the Phase 1b tables (`data/historical/adams_challis/`), then a separate
  dataset record for them.

Status, 2026-11-13 (the preparation is merged into `main` and pushed, 736bcd3; nothing
launched, no plan frozen; this table replaces the one of 12 November, which described
local branches):

| item | status | next step |
|---|---|---|
| M6d | no plan existed; draft in `docs/notes/V13_drafts.md` | settle the open choices, move the plan into `M6.md`, push, then run (about 7 CPU-h) |
| M6c variants, modern track | no plan existed; draft in the same file | as above (about 55 CPU-h with the full likelihood) |
| M6c variants, 1846 track | blocked: the 1846 reduction track is not built | decide whether v1.4 includes it |
| M6e | Lalande's two places of 1795 sourced and transcribed (`data/historical/precoveries/`); Galileo's not sourced; draft in the same file | source Galileo 1612–13 or drop it, then the plan |
| S1 cell 0 | no plan existed; draft in the same file (sub-cell partition) | choose the sub-cells, push, then run |
| Spot check of the Phase 1b tables | kit ready: 42 rows drawn with seed 20261112 (`docs/review/p1b_spotcheck_rows.csv`), 36 from the nine tables read from printed pages; 0 rows in error among the 36 gives an exact two-sided 95 % upper bound of 9.7 % on the row error rate; results template `docs/review/p1b_spotcheck_results_template.json` | the review itself, by a person; the HTML kit is built outside the repository by `docs/notes/scripts/p1b_spotcheck.py kit` |
| Dataset record for the Phase 1b tables | not started | after the spot check |
| Named input N1 of Phase 1b step 4 | done: merged and released in v1.2.1, where the companion page's chart case was corrected | none |
| Manuscript check of Adams's solutions | merged (29fa136), with the author's spot check of 13 Nov 2026 | none |
| Phase 2 scoping, physics, validation and first data set | merged (16c78d2, 10ca793, 407f4dc); nothing fitted; no plan frozen | a pre-registered plan before any fit of Mercury's observations |
| The attribution to Bouvard's tables | done and released in v1.3.0 (`docs/notes/P1B_bouvard_attribution.md`) | none |
| Changelog | the Unreleased section covers `main` since v1.3.0 | date it at the release |
