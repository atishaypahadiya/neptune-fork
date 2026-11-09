# Release checklist: v1.0.0

**No tag before M6 is complete.** The current version is `0.9.0`, covering M0–M5b plus the M6
data and observation model. `v1.0.0` is tagged only when M6 is complete. Versions follow
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
- [ ] M6c fits are done, with `site/data/m6_fits.json` and the figures the spec lists. (M6d and
  M6e moved to "Deferred to v1.1" below, 2026-11-09.)
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

## Repository state

- [ ] Branch merged to `main`, working tree clean, and nothing in `data/cache/` is needed to
  read the release.
- [ ] `uv sync --locked` succeeds on a fresh clone.
- [ ] `uv run ruff check . && uv run ruff format --check .` passes.
- [ ] `uv run pytest -m "not slow and not network"` passes locally and on CI (the badge is
  green on `main`).
- [ ] `uv run pytest -m slow` passes locally. Record the wall time.
- [ ] `uv run python examples/toy_fit.py` finishes in under 5 min on an idle machine. The time
  in `docs/reproduce.md` (4.6 min) was measured under A6c load, so re-time it and update.
- [ ] Committed `site/data/*.json` matches a fresh export where tests check that
  (`test_export.py`, `test_ledger.py`).
- [ ] `uv run python -m fork.export.ledger` re-run, and the CPU-hour totals in `README.md` and
  `docs/reproduce.md` updated (they will include A6c and M6c).
- [x] The A6d static check holds: no M6a–M6d code path reads Neptune's real ephemeris or the
  M6e precovery file (`tests/test_a6d.py`, 2026-11-09; runtime half in `tests/test_m6c.py`).
- [ ] Whitmore's *discoverneptune* is cited in the README. If any of its data are reused,
  check its license first.

## Metadata

- [ ] Version set to `1.0.0` in `pyproject.toml`, then `uv lock`, and the same value in
  `CITATION.cff` and `.zenodo.json`.
- [x] ORCID placeholder replaced in both files (2026-11-07):
  - `CITATION.cff` had `https://orcid.org/0000-0000-0000-0000`, now `https://orcid.org/0009-0004-1264-1662`;
  - `.zenodo.json` had `PLACEHOLDER-ORCID-replace-before-release`, now `0009-0004-1264-1662`.

  `grep -n "0000-0000-0000-0000\|PLACEHOLDER" CITATION.cff .zenodo.json` must print nothing.
- [ ] `CITATION.cff`: add `date-released: YYYY-MM-DD` and validate with
  `uvx cffconvert --validate`.
- [ ] `.zenodo.json` is valid JSON and agrees with `CITATION.cff`. Zenodo uses `.zenodo.json`
  over `CITATION.cff` when both exist.
- [ ] Zenodo GitHub integration enabled for `atishaypahadiya/neptune-fork`. This has to happen
  **before** the GitHub release is published, or no DOI is minted.

## Tag and release

- [ ] `git tag -a v1.0.0 -m "v1.0.0: Phase 1 through M6"`
- [ ] `git push origin v1.0.0`
- [ ] GitHub release from the tag. Release notes:
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
- [ ] Zenodo DOI minted. Add the DOI badge to `README.md` and `doi:` to `CITATION.cff`, in a
  follow-up commit.
- [ ] The Pages site still deploys from `main`.
