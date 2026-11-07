# Manifest

## Package files

| file | what it is |
|---|---|
| `leverrier_1846_table77.csv` | **The dataset.** Table n° 77, 279 rows, reconciled, exactly as printed |
| `leverrier_1846_misprints.csv` | Misprints in the source (7 cells), each with its argument and action |
| `leverrier_1846_resolution_log.csv` | The 14 cells where the two passes disagreed, and how each was resolved |
| `leverrier_1846_pass1.csv` | Pass 1 as transcribed (archive.org scan), before reconciliation |
| `leverrier_1846_pass2.csv` | Pass 2 as transcribed (Lille scan), before reconciliation |
| `spotcheck_rows.csv` | The 28-row spot-check sample and why each row was drawn |
| `spotcheck_results.json` | Automated spot-check result, cell by cell (336 cells) |
| `spotcheck_human_results.json` | The human spot-check result of the same 28 rows, cell by cell (336 cells, 0 discrepancies) |
| `leverrier_1846.csv` | Derived convenience file: one row per measured RA or Dec, in degrees, with epochs |
| `datapackage.json` | Frictionless Data Package descriptor: schemas, SHA-256 and size of every data file |
| `DATA_DICTIONARY.md` | Every column of every file: meaning, units, sign conventions, formats |
| `README.md` | What the dataset is, how it was made, its limitations |
| `MANIFEST.md` | This file: package contents, source scans, page mapping |
| `LICENSE`, `CITATION.cff`, `.zenodo.json` | Licence text and citation metadata |

The SHA-256 of each data file is in `datapackage.json` (`resources[].hash`).

## Verification records

Both spot checks reviewed the same 28 rows of `leverrier_1846_table77.csv` (SHA-256
`897f3cde1564998ac1aa1036c8665d91b24dda142d686e2fdf1ec5902e419127`) against the Gallica scan.

| record | reviewer | date | result | SHA-256 |
|---|---|---|---|---|
| `spotcheck_results.json` | Automated (AI model) | 2026-11-06 | 336 cells, 0 mismatches | `da90486b7905a8a2bf1796bca966f8a5bb9a5bad0c1c5a61e5d041e96a104ef9` |
| `spotcheck_human_results.json` | human: Atishay Pahadiya (dataset author) | 2026-11-07 | 336 cells, 0 discrepancies | `777515a11d0e2586767d3103835ef80b2beb0b73fa432650d9262bc1d1ff7a2e` |

## Source

U. J. Le Verrier, *Recherches sur les mouvements de la planète Herschel (dite Uranus)*. It was
printed in *Connaissance des Temps pour l'an 1849*, Additions, pp. 3–254 (Paris, Bachelier,
1846), and issued as an offprint (Bachelier, Paris, 1846) with the same page numbers. Table
n° 77, "Tableau de la comparaison des observations d'Uranus avec la théorie admise", is on
printed pp. 129–136. The columns are described in the text on p. 127. The source is in the
public domain.

## Source scans

The scans are not redistributed here. Each is identified by its SHA-256. To check one,
download it from the URL and run `shasum -a 256 <file>`. The file names are the local names
used during transcription.

| file | edition | source | SHA-256 |
|---|---|---|---|
| `leverrier_1846_archive_org.pdf` | 1846 offprint (Bachelier); Google Books scan, 264 PDF pages, page images ~1600 × 2500 px | https://archive.org/details/bub_gb_HMi0NyDVkhQC (`bub_gb_HMi0NyDVkhQC.pdf`) | `d17963fd996b3ae32e7ebe6e2b97b1964638ea718e280cf2bd8e78ff907ec042` |
| `leverrier_1846_lille.pdf` | 1846 offprint ("Extrait de la Connaissance des temps, pour 1849"); Université de Lille, 262 PDF pages, ~1890 × 2835 px | http://hdl.handle.net/1908/3096 → https://lillonum.univ-lille.fr/files/original/ffcb92a065dce94c1880b6ce3180498ec3449541.pdf | `298d2b57fdf5119b47a60c6449addb0d42964e4f7e4fe4e8484f0810d0ee51c6` |
| `cdt_1849_additions.pdf` | *Connaissance des Temps pour 1849*, Additions; BnF Gallica, views 468–725, 258 PDF pages, ~1000 × 1540 px, colour | https://gallica.bnf.fr/ark:/12148/bpt6k6505616w (views 468–725) | `bbb8bbcacd2d57abf25a5d75af623431e42925169b314bdd591cc28999c805d9` |

Roles in this dataset:
- Pass 1 was transcribed from the archive.org scan.
- Pass 2 was transcribed from the Lille scan.
- The Gallica scan of the *Connaissance des Temps* volume was the third witness. It was used
  for every disagreement and for the spot check.

The two offprint scans are two physical copies of the same 1846 offprint, not two
typesettings. One cell differs between copies: p. 133 l. 36, where the Lille copy prints 33′
and the other two print 22′ (see the resolution log).

A PDF downloaded from Gallica or Lille today may not reproduce these hashes, because a
library can re-render its files. Then the page mapping below still identifies the pages.

## Page mapping (printed page → PDF page, 1-based)

| file | PDF page | Table n° 77 (printed 129–136) | §78 equations (printed 137–141) | §131 33 conditions (printed 221–222) |
|---|---|---|---|---|
| `leverrier_1846_archive_org.pdf` | printed + 6 | 135–142 | 143–147 | 227–228 |
| `leverrier_1846_lille.pdf` | printed − 2 | 127–134 | 135–139 | 219–220 |
| `cdt_1849_additions.pdf` | printed + 2 | 131–138 | 139–143 | 223–224 |

The offsets were checked against the printed page headers. archive.org: PDF 135 = printed
129 and PDF 227 = printed 221. Lille: 127 = 129 and 219 = 221. Gallica: 131 = 129 and
223 = 221. The memoir ends on Gallica PDF p. 256, and Liouville's paper starts on p. 257.

The `witness` column of the resolution log and the reasons in the misprint table cite Gallica
**PDF** pages (printed + 2). Every `page` column in the CSV files is the **printed** page.

## Rows per printed page

| page | 129 | 130 | 131 | 132 | 133 | 134 | 135 | 136 | total |
|---|---|---|---|---|---|---|---|---|---|
| rows | 31 | 37 | 34 | 37 | 37 | 34 | 39 | 30 | 279 |

Rows 1–19 of p. 129 are the 19 ancient observations (1690–1771). The other 260 are modern
(1781–1845). Both passes counted these independently. The total agrees with the observation
counts of Le Verrier's 33 conditions (§131, p. 222).
