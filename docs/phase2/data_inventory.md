# Phase 2 data inventory: the observations of Mercury used by Le Verrier and Newcomb

Scoping note, 2026-11-12. Not a pre-registration. All URLs were requested on 2026-11-12 and
answered HTTP 200 unless stated. Nothing was downloaded into the repository.

## Unit of effort

One Phase 1 unit is the work on Le Verrier's Table 77: 279 rows by 12 data columns (3,348 cells)
on 8 printed pages, entered twice, reconciled against a third scan, with a 28-row human spot
check. Two conversions are given because Mercury's sources differ in layout:

- by cells: cells / 3,348;
- by pages: printed pages / 8.

Where the two disagree the range is shown. Free-text columns and mixed time scales are not
captured by either number and are noted in words.

## How statements were checked

SCAN: page images or OCR text read. UNVERIFIED: not confirmed. Row counts tagged "OCR tally" were
counted from OCR text of tables, which is poor; they need a recount from page images.

## Summary table

| id | set | printed in | scan | size | units |
|---|---|---|---|---|---|
| A | Le Verrier, transit timings 1631-1848 | Annales tome V, pp. 34-49 | URL 1 | a few hundred timings (UNVERIFIED) | about 2 |
| B | Le Verrier, 21 transit equations | Annales tome V, p. 82 | URL 1 | 21 rows | 0.07 |
| C | Le Verrier, Paris meridian places 1801-1842 | Annales tome V, pp. 51-61 | URL 1 | 397 rows, 5 columns | 0.6 to 1.4 |
| D | Newcomb, observed geocentric contacts 1677-1881 | APAE 1, part VI, pp. 451-452 | URL 2 | 55 to 60 rows (OCR tally) | 0.1 to 0.2 |
| E | Newcomb, individual contact observations | APAE 1, part VI, pp. 384-421 | URL 2 | of order 1,000 rows (UNVERIFIED) | 1.5 to 5 |
| F | Newcomb, station longitudes | APAE 1, part VI, pp. 373-377 | URL 2 | not counted | below 0.2 |
| G | Newcomb, transit equations of condition | APAE 1, part VI, from p. 453 | URL 2 | not counted | below 0.5 |
| H | Meridian observations used in 1895 | observatory annual volumes | not searched | Mercury 5,421 | 15 or more |
| I | Newcomb, secular variations table | 1895 volume, pp. 109-110 | URL 3 | 15 rows | 0.03 |
| J | Morrison and Ward, timings 1677-1973 | location not established | URL 4 (paper) | about 2400 | unknown |

"Annales tome V" is Annales de l'Observatoire imperial de Paris, tome V, chapter XV. "APAE 1" is
Astronomical Papers prepared for the use of the American Ephemeris and Nautical Almanac, vol. 1.

How each row was checked:

- A: pp. 34, 35 and 49 read as page images; the pages between were seen only as thumbnails. The
  section is 16 pages of narrative with small per-observer tables. The 1848 transit alone has 18
  timings (p. 49). Units by pages; the cells cannot be counted without reading every page.
- B: p. 82 read. About 11 coefficients per row. Not usable alone (see the notes).
- C: pp. 49-52 and 61 read. 240 observations for 1801-1828 and 157 for 1836-1842 (p. 49), in 195
  numbered groups (p. 61). Columns: mean time, apparent longitude, table error in longitude,
  apparent latitude, table error in latitude. 0.6 units by cells, 1.4 by pages.
- D: page range and layout from OCR text. Rows are November and May transits, interior and
  exterior contacts. Columns: contact, date, observed Greenwich mean time, weight, tabular time,
  5 or 6 symbolic coefficients. 0.1 units for the observed columns, 0.2 with the coefficients.
- E: page range from the memoir's contents (OCR). 38 pages. One row per observer and contact:
  station, observer, phase description (free text), observed local time, reduction to the Earth's
  centre, concluded Greenwich time. 1.5 units by cells at 1,000 rows, 5 by pages.
- F, G: page ranges from the memoir's contents (OCR). G is residual-type data.
- H: the 1895 volume prints counts, normal equations and results, not the observations (p. 3,
  OCR). Counts: Sun 40,176; Mercury 5,421; Venus 12,319; Mars 4,114; total 62,030. Observatories:
  Greenwich 1750-1892, Palermo 1791-1813, Paris 1801-1889, Konigsberg 1814-1845, Dorpat
  1823-1838, Cambridge 1828-1844, Berlin 1838-1842, Oxford (Radcliffe) 1840-1887, Pulkowa
  1842-1875, Washington 1846-1891, Leiden 1863-1871, Strassburg 1884-1887, Cape 1884-1890.
  Mercury alone is at least 6 units by cells at 4 columns per row, spread over more than a
  hundred volumes. Not costed further.
- I: pp. 109-110 from OCR. 15 rows by 3 values with mean errors. A summary, not observations.
- J: MNRAS 173, 183-206; pp. 184-185 read as page images. About 2400 internal-contact timings,
  transcribed from the original publications except 50 from Newcomb 1882. Whether the timings
  themselves are printed anywhere was not established.

URLs:

1. https://articles.adsabs.harvard.edu/pdf/1859AnPar...5....1L (PDF, 195 pages, about 10.5 MB,
   image-only: no text layer)
2. https://archive.org/details/astronomicalpap00newcgoog (526 images; OCR text at
   https://archive.org/download/astronomicalpap00newcgoog/astronomicalpap00newcgoog_djvu.txt).
   A second copy: https://archive.org/details/01Astronomicalpap00newcgoog
3. https://archive.org/details/elementsoffourin00newcuoft (234 images). Other copies:
   elementsoffourin00newcrich, cihm_16774.
4. https://articles.adsabs.harvard.edu/pdf/1975MNRAS.173..183M (PDF, 24 pages)

Licences of the scans were not assessed. Phase 1 committed stable URLs and page references, not
images, where reuse was unclear; the same rule would apply.

## Notes by set

### A, B: Le Verrier's transits

- Le Verrier begins his list with the transit of 7 November 1631 (Gassendi), which he describes
  and does not compute, then 1651 and 1661 (pp. 34-35, SCAN). His equations of condition start
  with 1697 (p. 82, SCAN). Whether he prints usable timings for 1677 was not checked.
- Times are given in Paris mean time for 1848 (p. 49) with the longitude and latitude of each
  station. Earlier transits may be in apparent time (not checked).
- The 21 equations (set B) are residuals against his own provisional tables, in terms of
  corrections to his elements. They cannot be used without those tables (Section V of the
  chapter, roughly pp. 110-190). This is the Phase 1 stop condition "only residuals against his
  own tables": set B fails it, set A passes it.

### C: Le Verrier's meridian positions

- These are positions, not only residuals: geocentric apparent longitude and latitude at a stated
  Paris mean time (civil, from midnight; p. 51, SCAN). The reduction details are not given (p. 50).
  The Sun's place came from Bessel's tables with corrections (p. 51).
- From 1836 the illuminated limb was observed and reduced to the centre with the apparent
  diameter; before that the centre was observed and no correction was applied (p. 50). This is a
  phase-defect systematic of the order of arcseconds.
- The printed table errors scatter by several arcseconds (pp. 51, 52, 61).
- Le Verrier states that the perihelion excess comes from the transits alone (p. 98). So set C is
  the direct analogue of Table 77 in format, but it is not where the signal is.

### D, E, F, G: Newcomb 1882

- Set D is already reduced: each row is a Greenwich mean time of geocentric contact, combined by
  Newcomb from set E with his own choices of phase, longitude and weight (sections 3 to 7 of
  Part I, pp. 369-379). No precise weight formula was used; five accordant internal contacts get
  weight 1 and the maximum is 6 (p. 451).
- All observations of first external contact were rejected; fourth contacts get low weight
  (p. 378).
- The observed times carry a symbolic clock term, minus k times a lunar-derived correction, which
  must be stripped or modelled (p. 450). It appears in the OCR as "- 33k", "- 26k" beside the
  times.
- The transit of 1782 needed a rigorous reduction to the Earth's centre because Mercury passed
  near the limb (p. 379).
- Transit years seen in the OCR of the comparison table: 1677, 1690, 1697, 1723, 1736, 1740,
  1743, 1753, 1756, 1769, 1782, 1786, 1789, 1799, 1802, 1809, 1822, 1832, 1845, 1848, 1861, 1868,
  1878, 1881 (OCR tally; some years appear only as exterior contacts).
- The archive.org OCR of the tables is unusable for transcription (digits and signs are lost).
  Transcription has to be done from page images, as in Phase 1.

### H: the meridian observations of 1895

- The 1895 volume does not print them. Newcomb reports that meridian observations of Mercury
  show systematic discordances that vary around the relative orbit (sections 30 and 33), and that
  he gave the transit equations more weight for that reason (p. 64, SCAN of OCR).
- The transit input to the 1895 solution is the 1882 result: V, W and N "from all the transits up
  to 1881" (pp. 61-62, SCAN of OCR). The transits of 1891 and 1894 were not found in that section.

### I: the secular-variation table

- Values are listed in `prior_art.md`, section 3. Newcomb warns that the theoretical mean errors
  are not independent because they come from the same masses (p. 109).
- Useful as a cheap multi-planet consistency check. It is a digest made by Newcomb with full
  knowledge of the anomaly. A test on it is not blind.

## The cheapest primary set for a first blind test

**Set D, the observed columns of Newcomb's comparison table, pp. 451-452.** About 55 to 60 rows
on two pages, about 0.1 unit.

Reasons:

1. It is the data that produced the 43 arcsec. Le Verrier and Newcomb both say the meridian
   observations do not carry the result.
2. It spans 204 years with one row per contact, already at the Earth's centre and in one time
   scale. The observation model is then one function: predicted time of geocentric contact.
3. It has weights, which give a first error model without reading 38 pages.
4. It is two pages. The whole Phase 1 transcription method (double entry, third source, spot
   check of every row) costs an afternoon.

What it costs in rigour, to be stated in any plan built on it:

1. It is not raw. Phase choices, longitudes, rejected observations and weights are Newcomb's.
   He knew the anomaly when he made them.
2. Blindness needs care. The same table prints the tabular times and the page after it the
   residuals. Only the columns contact, date, observed time and weight may be transcribed into the
   file the fit reads. The rest stays out of the pipeline's reach until the result is frozen.
3. The time scale is mean solar time. A clock correction (Delta T) for 1677-1881 is needed and
   its early values are uncertain at the level of the timing errors (see `design_sketch.md`).
4. The contact definition (which limb phase, which solar and planetary semidiameters) enters as
   nuisance parameters.

Second step, if the first is informative: set E for the transits that dominate the weight, to
replace Newcomb's combinations with an explicit per-observer error model (1.5 to 5 units), with
Morrison and Ward's method as the guide. Set A gives an independent partial check to 1848.
Set C is the cheapest positional set (about 1 unit) but should be expected to add little.
Set H is out of reach for a first phase.

## Not established

- Row counts for sets A, D and E (marked above).
- Whether Le Verrier used the transit of 1677.
- The image (leaf) numbers of pp. 451-452 in the archive.org scans.
- Whether a machine-readable collection of the timings exists (Morrison and Ward 1975; later
  solar-radius studies such as Sveshnikov, Astronomy Letters 28, 115-120 (2002), whose abstract
  says it analyses internal-contact timings for 1631-1973).
- The first and last printed pages of Newcomb's memoir as a whole (363 and about 485; the OCR is
  damaged at both).
