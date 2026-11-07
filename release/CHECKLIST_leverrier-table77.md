# Release checklist: Le Verrier Table n° 77 dataset

This is a maintainer file, not part of the package. Do not upload it to Zenodo.

Do the steps in order. Run the commands from the repository root.

**Status (2026-11-07):** step 2 is complete. The dataset was spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies), by Atishay Pahadiya:
`docs/review/spotcheck_human_results.json`, SHA-256
`777515a11d0e2586767d3103835ef80b2beb0b73fa432650d9262bc1d1ff7a2e`. The package is out of
the "pending" state. The ORCID is filled in (step 1). Next: Zenodo (step 3).

## 1. Decisions (done 2026-11-07)

- [x] Keep both raw passes and the derived file `leverrier_1846.csv` in the package.
- [x] Licence: CC0 1.0 for the data and documentation, with a "please cite" request in the
  README.
- [x] Publish only after the human spot check.
- [x] ORCID: 0009-0004-1264-1662 (https://orcid.org/0009-0004-1264-1662), in `leverrier-table77/CITATION.cff` (as the URL) and
  `.zenodo.json` (as the bare iD, Zenodo's format). Done 2026-11-07.

## 2. Human spot check (done 2026-11-07)

Result: spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies), by the dataset author. `--check` printed `28 rows, 336 cells, 0 marked as
not matching`. The record's SHA-256 is in `leverrier-table77/MANIFEST.md`.

- [x] **Have the scan.** Put the Gallica scan at
  `data/historical/sources/cdt_1849_additions.pdf`. Check
  `shasum -a 256 data/historical/sources/cdt_1849_additions.pdf` against
  `leverrier-table77/MANIFEST.md`.
- [x] **Generate the page.** Run
  `uv run --with pymupdf python -I -m fork.data.spotcheck --human`. This writes
  `docs/review/spotcheck_human.html` and the row images in `docs/review/crops/`. The images
  are not committed.
- [x] **Review.** Open `docs/review/spotcheck_human.html` in a browser.
  - Enter your name.
  - Tick each cell that matches the print. Judge from the row image, or from the
    full-resolution Gallica link when an image is unclear.
  - Leave a cell unticked if it does not match, and say why in that row's notes.
  - Don't open `spotcheck_results.json` (the automated spot check) or the transcription notes until
    you have exported.
  - "All match" is disabled on the 16 logged rows, so tick those cell by cell.
- [x] **Export.** Click **Export results**. The browser downloads
  `spotcheck_human_results.json`. Move it, unedited, to
  `docs/review/spotcheck_human_results.json`.
- [x] **Record your spot-check results.** Run
  `uv run python -I -m fork.data.spotcheck --human --check docs/review/spotcheck_human_results.json`.
  It must print `28 rows, 336 cells, 0 marked as not matching` and exit 0. It checks that:
  - the file is for the current table (SHA-256);
  - every cell value is the table's;
  - the rows are the reviewed sample;
  - the file is marked `reviewer: human` with a name.
- [x] (Not needed: 0 discrepancies.) **If any cell is marked as not matching: stop.** Don't release.
  - A correction goes through the pass files and a logged entry, never a hand edit of the
    table.
  - It changes the table's SHA-256, so both spot-check records then need redoing for the new
    table. Plan that before fixing.
- [x] **Rebuild the package.** Run `uv run python -m fork.export.release`. It should print
  `human spot check: included`. The build refuses a human record that fails the check above.
- [x] **Take the dataset out of the "pending" state.** `tests/test_release.py` fails until
  you do. Replace the "A human spot check is pending" wording with "spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies)" in:
  - `leverrier-table77/README.md`: the "Read this first" paragraph, the first Limitations
    bullet, the "How it was made" spot-check step (add the human result), and the Contents
    table (add `spotcheck_human_results.json`);
  - `leverrier-table77/CITATION.cff`, the abstract;
  - `leverrier-table77/.zenodo.json`, the description and notes;
  - `leverrier-table77/MANIFEST.md`: add a row for `spotcheck_human_results.json`.

  `datapackage.json` updates itself. Keep stating that transcription and reconciliation were
  machine-assisted (an AI model), not done by an independent human.
- [x] **Test and commit.** Run `uv run pytest -q tests/test_release.py tests/test_spotcheck.py`.
  Then commit `docs/review/spotcheck_human_results.json`, `docs/review/spotcheck_human.html`
  and the release folder.

## 3. Zenodo

- [x] Create an account at https://zenodo.org. Sign in with GitHub or ORCID.
- [ ] Optional dry run on https://sandbox.zenodo.org. It needs a separate account and issues
  no real DOIs.
- [ ] Choose **New upload**. Add every visible file in `release/leverrier-table77/`, flat (no
  zip). `.zenodo.json` is not read by the web form: use it as the text to paste.
- [x] Under "Do you already have a DOI?" choose **No**, then **Get a DOI now**, which reserves
  the DOI. **Reserved: 10.5281/zenodo.24044792** (https://doi.org/10.5281/zenodo.24044792), 2026-11-07. It does not resolve until the record
  is published, so nothing (tests included) may fetch it.
  - [x] DOI in `README.md` ("How to cite" section, with a BibTeX entry; no badge yet).
  - [x] DOI in `CITATION.cff` (`doi` and an `identifiers` entry of type doi). It is validated
    with cffconvert.
  - [x] DOI as the Data Package `id` in `datapackage.json` (set by `fork.export.release`).
  - [x] `.zenodo.json` left without a DOI: Zenodo assigns the reserved DOI to the draft
    itself.
  - [ ] Upload, or replace in the draft if already uploaded, `README.md`, `CITATION.cff` and
    `datapackage.json`. The last one changed too, because it now carries the `id`.
  - [ ] At publication: set `date-released` in `CITATION.cff` to the publication date.
    Changing the README or CITATION.cff after that needs a new version.
- [ ] Fill in the form from `.zenodo.json`:
  - resource type Dataset, the title, and the creator with ORCID 0009-0004-1264-1662;
  - paste the description;
  - licence **CC0 1.0 Universal**, the keywords, version 1.0.0, language English;
  - related works: "is supplement to" the GitHub repository, and "is derived from" each of the
    three scan URLs;
  - the Le Verrier reference.
- [ ] **Publish.** Files are frozen after this, and any change becomes a new version.

## 4. Afterwards

- [ ] Push and merge `data-release`. Add the DOI badge for 10.5281/zenodo.24044792 to the dataset and project
  READMEs, only after publication, since before that it points at a 404. Then remove the
  no-badge guard in `tests/test_release.py`.
- [ ] Zenodo also mints a concept DOI (all versions) at publication. Consider adding it to
  `CITATION.cff` `identifiers` in the next version.
- [ ] A later change to any file in `data/historical/` or `docs/review/` that the package
  copies fails `tests/test_release.py`. Then rebuild, and publish a new Zenodo version under
  the same concept DOI.
