# Phase 2 data: Newcomb's observed contact times, transcribed

12 November 2026. Nothing is fitted to these data here or anywhere yet: fitting Mercury's
observations to competing models is the Phase 2 experiment and needs a frozen plan first.

## Choice of set

`data_inventory.md` names set D as the cheapest primary set for a first
blind test: the observed columns of Newcomb's comparison table, *Astronomical Papers of the
American Ephemeris* vol. 1, part VI (1882), pp. 451–452. It is the data behind the 43″, it
spans 1677 to 1881 in one time scale, already reduced to the Earth's centre, with weights,
on two pages. Le Verrier's 1859 chapter was not used here: its transit timings are spread
over 16 pages of narrative and its scan has no text layer to locate them by.

What this set is not: raw. Each time is Newcomb's combination of many observers' timings,
with his choices of phase, longitude and weight, made knowing the anomaly. Any plan built on
it has to say so.

## Method (as for Table 77)

1. **Scan and images.** The archive.org scan, identified by hash in
   `data/historical/mercury/MANIFEST.md`. The embedded page scans of pp. 451 and 452 were
   cut into 14 horizontal bands at full resolution
   (`docs/phase2/scripts/mercury_crops.py`), each showing only the columns to be
   transcribed. The first cut stopped short of the foot of tables I and II; three bands
   were added before either entry was complete.
2. **First entry** from the bands.
3. **Second entry**, made from the bands and the written column specification alone, with
   no sight of the first entry, of any text layer or of any other file.
4. **Comparison**, cell by cell: 58 rows × 12 fields = 696 cells.
5. **Reconciliation** on the named band for anything in doubt.
6. **Misprint log**: inconsistencies in the print are recorded, never corrected.
7. **Checks the table allows** (`tests/test_mercury_transcription.py`).

## Result

| | |
|---|---|
| rows | 58 (25 + 12 + 15 + 6) |
| cells compared | 696 |
| disagreements between the two entries | 0 |
| cells marked doubtful by an entry | 2 (table I row 1, seconds: 1 or 2, read 1; table IV row 4, a broken plus sign) |
| inconsistencies in the print, logged | 3 |

Both entries share one method of reading and one set of images, so their errors may be
correlated; agreement is weaker evidence than two unrelated readings would give. No second
printing or scan was compared (the second archive.org item appears to be the same file). A
person's spot check has not been made.

Logged, not corrected:

- 1782 November 12: both rows are printed as contact III. The earlier one (2h 42m 6s) is
  contact II by the table's pattern.
- 1848: the exterior contact is dated November 8 at 4h 29m 40s; table I dates contact III
  of the same transit November 9 at 4h 28m 8s. One of the two dates is a day off.
- 1848, exterior: "0k" without a sign, against "− 0k" in table I. No effect.

Checks that pass on the table as printed: one k term per transit across all four tables (23
transit years); within a transit the two interior contacts are 1 to 8 hours apart; each
exterior contact follows the last interior contact by 1 to 8 minutes (1782, the transit
near the limb, is the longest at 6.5 minutes); weights lie in (0, 6], the range the text
states.

## What a user of the file must know

- Times are Greenwich mean time on the astronomical day, which begins at noon: "1677 Nov. 6,
  21h 34m" is the morning of 7 November, civil.
- Each time carries a term in k, a correction to the adopted longitude or clock scale that
  Newcomb leaves symbolic (explained on p. 450 according to the data inventory; that
  page was not read for this note): the printed value is the time
  for k = 0, and the coefficient is common to a transit.
- Table II row 3 (1786 May 3, contact II) prints two times in a brace, 14h 59m 25s and
  15h 0m 8s. The first is in the time columns, the second in `second_value`.
- The tabular times, the symbolic corrections and the residuals on the pages that follow were not
  transcribed and are not in the repository. They should stay out of reach of any fitting
  code until a result is frozen.
- Mean solar time needs ΔT for 1677 to 1881 before comparison with a dynamical model.

## Reproduce

```
curl -L -o data/cache/mercury/astronomicalpap00newcgoog.pdf \
  https://archive.org/download/astronomicalpap00newcgoog/astronomicalpap00newcgoog.pdf
shasum -a 256 data/cache/mercury/astronomicalpap00newcgoog.pdf      # see the MANIFEST
uv run --with pymupdf --with pillow python -I docs/phase2/scripts/mercury_crops.py
uv run pytest -q tests/test_mercury_transcription.py
```

## Spot check by the author, 13 Nov 2026

Added below the note of 12 November, which stands as written. This section supersedes its
sentence "A person's spot check has not been made".

I checked 20 cells of the table against the page images
(`docs/review/newcomb_spotcheck_human_results.json`, exported 2026-11-13 06:34 UTC, seed
20261113, SHA-256 of `newcomb_1882_contacts.csv` `70d4d0d1…`). A cell is one field of one
row; the table has 560. All 20 match.

| how the cell was chosen | cells | match |
|---|---|---|
| marked doubtful by an entry | 2 | 2 |
| in the misprint log | 3 | 3 |
| the one braced second value | 1 | 1 |
| the one fractional second | 1 | 1 |
| random | 13 | 13 |

**Confirmed.** The two doubtful glyphs are as entered: the seconds of table I row 1 are 1,
and the sign of table IV row 4 is plus. The three logged inconsistencies are in the print
as logged: contact III for the first 1782 row, the day 8 for the exterior contact of 1848,
and its "0k" without a sign. The second time of the 1786 brace is 15h 0m 8s, and the
seconds of the 1868 contact III are 9.8.

**Scope.** 7 of the 20 cells were chosen, not drawn: they are the cells above. 13 were
random among the remaining 553. No discrepancy in 13 random cells bounds the cell error
rate only weakly: the two-sided 95 % upper limit is 1 − 0.025^(1/13) ≈ 25 %. The check is
therefore not an estimate of the error rate of the table. Its value is that every doubtful
or irregular reading was confirmed by a person. The images were whole tables cut from the
same scan, with the transcribed columns only.
