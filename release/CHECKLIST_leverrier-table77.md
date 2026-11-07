# Release checklist: Le Verrier Table n° 77 dataset

This is a maintainer file, not part of the package. Do not upload it to Zenodo.

Do the steps in order. **Do not publish before step 2 is complete.** Run the commands from the
repository root.

## 1. Decisions (done 2026-11-07)

- [x] Keep both raw passes and the derived file `leverrier_1846.csv` in the package.
- [x] Licence: CC0 1.0 for the data and documentation, with a "please cite" request in the
  README.
- [x] Publish only after the human spot check.
- [ ] ORCID: fill in the placeholder in `leverrier-table77/CITATION.cff` (and add `"orcid"` to
  the creator in `.zenodo.json`), or delete the placeholder.

## 2. Human spot check

- [ ] **Have the scan.** Put the Gallica scan at
  `data/historical/sources/cdt_1849_additions.pdf`. Check
  `shasum -a 256 data/historical/sources/cdt_1849_additions.pdf` against
  `leverrier-table77/MANIFEST.md`.
- [ ] **Generate the page.** Run
  `uv run --with pymupdf python -I -m fork.data.spotcheck --human`. This writes
  `docs/review/spotcheck_human.html` and the row images in `docs/review/crops/`. The images
  are not committed.
- [ ] **Review.** Open `docs/review/spotcheck_human.html` in a browser.
  - Enter your name.
  - Tick each cell that matches the print. Judge from the row image, or from the
    full-resolution Gallica link when an image is unclear.
  - Leave a cell unticked if it does not match, and say why in that row's notes.
  - Don't open `spotcheck_results.json` (the automated review) or the transcription notes until
    you have exported.
  - "All match" is disabled on the 16 logged rows, so tick those cell by cell.
- [ ] **Export.** Click **Export results**. The browser downloads
  `spotcheck_human_results.json`. Move it, unedited, to
  `docs/review/spotcheck_human_results.json`.
- [ ] **Record your spot-check results.** Run
  `uv run python -I -m fork.data.spotcheck --human --check docs/review/spotcheck_human_results.json`.
  It must print `28 rows, 336 cells, 0 marked as not matching` and exit 0. It checks that:
  - the file is for the current table (SHA-256);
  - every cell value is the table's;
  - the rows are the reviewed sample;
  - the file is marked `reviewer: human` with a name.
- [ ] **If any cell is marked as not matching: stop.** Don't release.
  - A correction goes through the pass files and a logged entry, never a hand edit of the
    table.
  - It changes the table's SHA-256, so both spot-check records then need redoing for the new
    table. Plan that before fixing.
- [ ] **Rebuild the package.** Run `uv run python -m fork.export.release`. It should print
  `human spot check: included`. The build refuses a human record that fails the check above.
- [ ] **Take the dataset out of the "pending" state.** `tests/test_release.py` fails until
  you do. Replace the "A human spot check is pending" wording in:
  - `leverrier-table77/README.md`: the "Read this first" paragraph, the first Limitations
    bullet, the "How it was made" spot-check step (add the human result), and the Contents
    table (add `spotcheck_human_results.json`);
  - `leverrier-table77/CITATION.cff`, the abstract;
  - `leverrier-table77/.zenodo.json`, the description and notes;
  - `leverrier-table77/MANIFEST.md`: add a row for `spotcheck_human_results.json`.

  `datapackage.json` updates itself. Keep stating that transcription and reconciliation were
  done by an AI model.
- [ ] **Test and commit.** Run `uv run pytest -q tests/test_release.py tests/test_spotcheck.py`.
  Then commit `docs/review/spotcheck_human_results.json`, `docs/review/spotcheck_human.html`
  and the release folder.

## 3. Zenodo

- [ ] Create an account at https://zenodo.org. Sign in with GitHub or ORCID.
- [ ] Optional dry run on https://sandbox.zenodo.org. It needs a separate account and issues
  no real DOIs.
- [ ] Choose **New upload**. Add every visible file in `release/leverrier-table77/`, flat (no
  zip). `.zenodo.json` is not read by the web form: use it as the text to paste.
- [ ] Under "Do you already have a DOI?" choose **No**, then **Get a DOI now**, which reserves
  the DOI. Then:
  - Put the DOI in `README.md`, replacing both `zenodo.XXXXXXX`.
  - Put it in `CITATION.cff`: uncomment `doi` and set `date-released`.
  - Commit, and replace those two files in the Zenodo draft. `datapackage.json` does not hash
    them, so no rebuild is needed.
- [ ] Fill in the form from `.zenodo.json`:
  - resource type Dataset, the title, and the creator (with ORCID if you have one);
  - paste the description;
  - licence **CC0 1.0 Universal**, the keywords, version 1.0.0, language English;
  - related works: "is supplement to" the GitHub repository, and "is derived from" each of the
    three scan URLs;
  - the Le Verrier reference.
- [ ] **Publish.** Files are frozen after this, and any change becomes a new version.

## 4. Afterwards

- [ ] Push and merge `data-release`. Add the DOI badge to the project README.
- [ ] A later change to any file in `data/historical/` or `docs/review/` that the package
  copies fails `tests/test_release.py`. Then rebuild, and publish a new Zenodo version under
  the same concept DOI.
