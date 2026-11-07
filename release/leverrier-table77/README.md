# Le Verrier (1846), Table n° 77: observations of Uranus, 1690–1845

In 1846 Urbain Le Verrier predicted the position of an unseen planet from the way Uranus
strayed from its predicted path. Galle found Neptune near that position on 23 September 1846.
The evidence was a table: 279 observed positions of Uranus from 1690 to 1845, each set beside
the position computed from the best theory then available, with the differences.

This dataset is a transcription of that table, Table n° 77 of his memoir. Every cell is kept
exactly as printed. Each step that produced it is recorded: two independent passes, a third
witness for every disagreement, an arithmetic check of every row, a log of misprints in the
source, and a spot check.

**Read this first:** the transcription and the reconciliation were done by an AI model
, not by an independent human. The dataset was then spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies). That check
was done by the dataset author. See [Limitations](#limitations).

## Contents

| file | content |
|---|---|
| `leverrier_1846_table77.csv` | **The dataset.** All 279 rows (8 printed pages), as printed |
| `leverrier_1846_misprints.csv` | 7 cells that are wrong in the print itself, with the argument for each |
| `leverrier_1846_resolution_log.csv` | The 14 cells where the two passes disagreed, and how each was resolved |
| `leverrier_1846_pass1.csv`, `leverrier_1846_pass2.csv` | The two independent passes, before reconciliation |
| `spotcheck_rows.csv`, `spotcheck_results.json` | The spot-check sample and the automated check: 28 rows, 336 cells |
| `spotcheck_human_results.json` | The human spot check of the same 28 rows: 336 cells, 0 discrepancies |
| `leverrier_1846.csv` | Derived: 529 RA and Dec measurements in degrees, with Julian Dates and observatories. Contains inferences |
| `DATA_DICTIONARY.md` | Every column: meaning, units, sign conventions, formats |
| `MANIFEST.md` | Source scans (URLs, SHA-256), PDF-to-printed page mapping, rows per page |
| `datapackage.json` | [Frictionless Data Package](https://specs.frictionlessdata.io/) descriptor, with typed schemas and the SHA-256 of every data file |

## What the table contains

Each row is one observation:
- the date and Paris mean time;
- the observed apparent right ascension and declination of Uranus;
- the seconds of the RA and Dec computed from Le Verrier's "théorie admise" (Uranus perturbed
  by Jupiter and Saturn, no unknown planet);
- the differences computed − observed (C − O) in RA, Dec, geocentric longitude and
  geocentric latitude.

There are two groups of observations:
- **19 "ancient" observations, 1690–1771**: sightings of Uranus made before Herschel
  discovered it in 1781, by Flamsteed, Bradley, Mayer and Le Monnier, re-reduced
  by Le Verrier. They print RA in arc.
- **260 meridian observations, 1781–1845**, from Greenwich and Paris. They print RA in time.

Things that are easy to get wrong (details in `DATA_DICTIONARY.md`):
- **The table prints C − O, not O − C.** The derived file's `lv_residual` is O − C.
- **RA units depend on the era.** The 1690–1771 rows give RA in degrees (° ′ ″) and their RA
  C − O in arcseconds. The 1781–1845 rows give RA in h m s and their RA C − O in seconds of
  time.
- **Computed values are printed as their seconds only.**
- **Blank year and month cells mean "same as above."**
- **`»` means no Dec was observed.**
- **The sign belongs to the whole sexagesimal string.** The 1843 Sep 20 declination is printed
  "− 0,33. 2,7" and stored as `-0 33 2.7`: −0° 33′ 2.7″, negative although its degree field
  is 0.

A minimal reader:

```python
import csv


def sexagesimal(s):  # "-0 33 2.7" -> -0.55075 (top unit: degrees, or hours for modern RA)
    sign = -1.0 if s.startswith("-") else 1.0
    a, b, c = (float(x) for x in s.lstrip("+-").split())
    return sign * (a + b / 60 + c / 3600)


rows, year, month = [], None, None
with open("leverrier_1846_table77.csv", encoding="utf-8") as fh:
    for r in csv.DictReader(fh):
        year, month = r["year"] or year, r["month"] or month  # blank = same as above
        rows.append({**r, "year": int(year), "month": month})
```

## Provenance

**Source.** U. J. Le Verrier, *Recherches sur les mouvements de la planète Herschel (dite
Uranus)*, in *Connaissance des Temps pour l'an 1849*, Additions, pp. 3–254 (Paris, Bachelier,
1846). It was also issued as an offprint with the same page numbers. Table n° 77, "Tableau de
la comparaison des observations d'Uranus avec la théorie admise", is on pp. 129–136. The text
describing its columns is on p. 127.

**Scans used.** None is redistributed here. URLs, SHA-256 hashes and page mappings are in
`MANIFEST.md`.
1. The 1846 offprint, Google Books scan on archive.org (`bub_gb_HMi0NyDVkhQC`): pass 1.
2. The 1846 offprint, Université de Lille (hdl:1908/3096): pass 2.
3. *Connaissance des Temps pour 1849*, BnF Gallica (`ark:/12148/bpt6k6505616w`): third
   witness and spot check.

Scans 1 and 2 are two physical copies of the same printing. Scan 3 is the volume the memoir
appeared in. All three are public domain.

## How it was made

1. **Two independent passes.** Two separate AI agents each transcribed the whole table.
   - Each worked from the page images of one scan: pass 1 from archive.org, pass 2 from
     Lille.
   - Neither used the scan's embedded OCR, and neither saw the other scan or the other pass.
   - Both recorded the cells exactly as printed, and wrote `?` for a digit they could not
     read.
   - Both found 279 rows with the same row structure (31, 37, 34, 37, 37, 34, 39 and 30 per
     page). The total matches the observation counts of Le Verrier's own §131 (p. 222).
2. **Comparison.** The passes were aligned row by row and compared cell by cell.
   - 14 of the 3348 data cells (279 × 12) disagreed.
   - In 10 of them, one pass had `?` and the other a digit.
   - 4 were real conflicts.
3. **Third witness.** Every disagreement was resolved from the Gallica page image of the
   *Connaissance des Temps* volume, and logged with the page
   (`leverrier_1846_resolution_log.csv`).
   - The third witness confirmed the other pass in all 10 `?` cases.
   - In one conflict neither pass had the right value (1792 Nov 12).
   - One conflict is a real difference between copies: the Lille copy prints 33′ where the
     other two print 22′ (1819 Jun 20).
   - One cell was undecidable from the glyph. The year on p. 135 l. 4 has a bold 3 and 5 that
     are near-identical in all three scans. It was resolved as 1835 from the table itself
     (chronological order and the RA progression between oppositions), and the dataset author
     confirmed it on the image.
4. **Arithmetic check of every row.** In each row the printed C − O must equal the computed
   seconds minus the observed seconds. The computed minute is not printed, so the nearest
   consistent minute is used.
   - All 279 RA rows and 250 Dec rows pass to within 0.05, except one: 1764 Jan 15, Dec. That
     row is printed the same way in all three witnesses and is listed as a misprint.
   - The check catches misread seconds and C − O digits. It cannot see an error in the degrees,
     hours or minutes.
5. **Misprint policy.** The transcription never corrects the source. A cell judged wrong in
   the print goes into `leverrier_1846_misprints.csv`, with:
   - the printed value and, where it can be determined, the correction;
   - an action: `apply` (the derived file uses the correction) or `log_only`;
   - an argument built from the table itself: neighbouring rows, daily motion, the row's
     arithmetic, ecliptic geometry, and the meridian-transit identity for times. No argument
     depends on a modern ephemeris.

   There are 7 entries, 4 `apply` and 3 `log_only`:
   - two wrong degree digits in a declination (1771 Dec 18, 1830 Nov 11; `apply`);
   - three mean times (1804 Mar 31 and 1807 Apr 25, `apply`; 1768 Dec 30, `log_only`, cause
     not found);
   - one Dec C − O whose sign contradicts its own row (1764 Jan 15, `log_only`);
   - one probable slip in a computed RA (1830 Nov 11, `log_only`).

   An AI model checked every entry as printed in the Gallica image. The dataset author approved
   each `apply` entry.
6. **Spot checks.** 28 rows (10 %) were checked cell by cell against the Gallica images,
   twice.
   - The sample is stratified with a fixed seed: every resolution-log row, the misprint rows
     known at review time, at least two ancient rows, every page, then random rows.
   - **First by an AI model:** 336 cells, 0 mismatches (`spotcheck_results.json`).
   - **Then by a human, the dataset author:** the dataset was spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies)
     (`spotcheck_human_results.json`, 2026-11-07).
   - The human review page did not show the automated results, and nothing on it was pre-ticked.
   - Neither review used the "all match" shortcut: every cell was ticked individually.

**Afterwards: comparison with a modern ephemeris.** Each observed RA and Dec was compared
with Uranus's apparent place from JPL Horizons at the same epoch. This is not part of the
transcription method: it was done for the analysis this dataset was made for, and it is how
the 1771 misprint was found.
- After the four corrections, all 529 measurements lie within 13.5″ of JPL. The median
  |O − JPL| is 0.9–1.9″ for the modern series and 3.1–3.5″ for the ancient ones.
- A misread degree or minute digit in an observed value, or a misread seconds digit in a
  modern RA (1 s of time is 15″), would show as a residual far larger than any seen. A
  misread tens-of-arcseconds digit in a modern Dec (10″) would also stand out against the
  1–3″ scatter of those series.
- The two largest modern residuals (+9.1″ and +9.2″ in RA) are matched by a similar jump in Le
  Verrier's own C − O. So they are in the observations, not in the transcription.

## Limitations

- **AI transcription; human check of a sample only.**
  - Both passes, the reconciliation and the spot check were done by an AI model, in separate
    sessions and from different scans, but by one model family.
  - Errors that model makes systematically are not independent between passes. Reading the
    same ambiguous glyph the same way in both passes is the main risk. Such a cell produces no
    disagreement, so it never reaches the third witness, and the automated spot check may share the
    same bias.
  - The dataset was spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies). That check guards against such errors only in the rows it covers.
  - It covers the same 28 rows as the automated check. **The other 251 rows have had no human
    check.**
  - The human reviewer is the dataset author, who also directed the work, approved the
    misprint corrections and confirmed one disputed glyph. It is not an independent
    third-party review.
- **Low-confidence agreements.** In 33 rows a transcriber noted a low-confidence reading
  (`low_conf:` in `notes`). Where both passes agreed, the cell was not sent to the third
  witness unless it fell in the spot-check sample.
- **Unequal checks across columns.**
  - The observed RA and Dec are covered by the double entry, the arithmetic check (seconds)
    and the ephemeris comparison (larger digits).
  - The computed seconds and the RA and Dec C − O are covered by the double entry and the
    arithmetic check.
  - **The longitude and latitude C − O columns have no internal check.** Only the double
    entry and the two spot checks (28 rows) stand behind them.
  - The ephemeris comparison cannot exclude errors below the observational scatter: tenths of
    a second of time, units of arcseconds, and tens of arcseconds in an ancient row.
- **The witnesses share one typesetting.** The two offprint scans are copies of the same
  printing, and the *Connaissance des Temps* volume appears to share its type: the misprints
  found occur in all three. The three witnesses guard against damaged copies and poor scans,
  not against errors in the type. Undetected misprints may remain, especially where a wrong
  computed value and a wrong C − O are consistent with each other.
- **262 or 260.** Le Verrier's text (p. 126) says 262 modern observations were chosen, but the
  table has 260. The §131 counts agree with the table, and nothing in the table marks two
  rejected rows. This is not resolved.
- **The derived file contains interpretation.** `leverrier_1846.csv` adds:
  - the time meridian (Paris, inferred from the transit identity except where stated);
  - observatories (from the text; the two 1828 rows inferred);
  - ΔT, from a polynomial.

  Use `leverrier_1846_table77.csv` for what was printed.
- **Le Verrier's own numbers carry his choices.** The ancient positions are his
  re-reductions, and the computed column is his 1846 theory. Some of his computed values are
  probably slips: see the `log_only` entry for 1830 Nov 11.

## Validation

`datapackage.json` validates with the Frictionless framework:

```sh
pip install frictionless
frictionless validate datapackage.json
```

This checks the column types and formats, primary keys, the foreign keys from the misprint
table, resolution log and spot-check sample to the table, and every file's SHA-256 and size.
The project repository runs the same validation in its tests, together with the row counts
(279; per page 31, 37, 34, 37, 37, 34, 39, 30).

## Citation

The dataset is in the public domain (see [Licence](#licence)), but please cite it and Le
Verrier's memoir:

> Pahadiya, A. (2026). *Le Verrier (1846), Table n° 77: 279 observations of Uranus,
> 1690–1845, transcribed and reconciled* (Version 1.0.0) [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.XXXXXXX

> Le Verrier, U. J. (1846). Recherches sur les mouvements de la planète Herschel (dite Uranus).
> *Connaissance des Temps pour l'an 1849*, Additions, 3–254. Paris: Bachelier.

Machine-readable citation metadata is in `CITATION.cff`.

## Licence

The 1846 source is in the public domain, and so is this dataset: the transcription, its logs
and its documentation are dedicated to the public domain under
[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) (`LICENSE`). You may
copy, modify and redistribute them without asking.

**Please cite anyway.** CC0 removes the legal requirement, not the scholarly one. If you use
this dataset, cite it by its DOI (https://doi.org/10.5281/zenodo.XXXXXXX; the entry above,
or `CITATION.cff`), and cite Le Verrier (1846) for the observations. Citing the version you
used also tells readers which transcription and which checks stand behind your numbers.

## Related work

- Bertrand (2025, *Gazette des Mathématiciens* 186) counts the same 279 data points.
- J. Whitmore, *discoverneptune* (2026), https://jbwhit.github.io/discoverneptune/, is a
  REBOUND-based reconstruction of the prediction.
- Lai, Lam & Young (1990) treat the Uranus perturbations as a forward and an inverse problem.

This dataset was made for [neptune-fork](https://github.com/atishaypahadiya/neptune-fork),
which rebuilds the prediction from these observations as a Bayesian inference.

## Acknowledgements

Page images: Internet Archive / Google Books; Université de Lille (Lillonum); Bibliothèque
nationale de France (Gallica).
