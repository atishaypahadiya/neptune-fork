# Historical source scans

The PDFs are not committed (`.gitignore`). Each one is identified by its SHA-256 below.
Re-download from the URL and check the hash with `shasum -a 256 <file>`.

All three are the same text: Le Verrier, *Recherches sur les mouvements de la planète Herschel
(dite Uranus)*. It was printed in *Connaissance des Temps pour 1849*, Additions, pp. 3–254
(Paris, 1846), and issued as an offprint (Bachelier, Paris, 1846) with the same page numbers.
All three are in the public domain.

| file | edition | source | SHA-256 |
|---|---|---|---|
| `leverrier_1846_archive_org.pdf` | 1846 offprint (Bachelier); Google Books scan, 264 PDF pages, page images ~1600 × 2500 px | https://archive.org/details/bub_gb_HMi0NyDVkhQC (`bub_gb_HMi0NyDVkhQC.pdf`) | `d17963fd996b3ae32e7ebe6e2b97b1964638ea718e280cf2bd8e78ff907ec042` |
| `leverrier_1846_lille.pdf` | 1846 offprint ("Extrait de la Connaissance des temps, pour 1849"); Université de Lille, 262 PDF pages, ~1890 × 2835 px, ABBYY FineReader 9 | http://hdl.handle.net/1908/3096 → https://lillonum.univ-lille.fr/files/original/ffcb92a065dce94c1880b6ce3180498ec3449541.pdf | `298d2b57fdf5119b47a60c6449addb0d42964e4f7e4fe4e8484f0810d0ee51c6` |
| `cdt_1849_additions.pdf` | *Connaissance des Temps pour 1849*, Additions; BnF Gallica, views 468–725, 258 PDF pages, ~1000 × 1540 px, colour | https://gallica.bnf.fr/ark:/12148/bpt6k6505616w (views 468–725) | `bbb8bbcacd2d57abf25a5d75af623431e42925169b314bdd591cc28999c805d9` |

## Page mapping (printed page → PDF page, 1-based)

| file | PDF page | Table n° 77 (printed 129–136) | §78 equations (printed 137–141) | §131 33 conditions (printed 221–222) | memoir ends |
|---|---|---|---|---|---|
| `leverrier_1846_archive_org.pdf` | printed + 6 | 135–142 | 143–147 | 227–228 | — |
| `leverrier_1846_lille.pdf` | printed − 2 | 127–134 | 135–139 | 219–220 | — |
| `cdt_1849_additions.pdf` | printed + 2 | 131–138 | 139–143 | 223–224 | 256 (Liouville starts at 257) |

These were checked against the printed page headers: archive.org 135 = 129 and 227 = 221;
Lille 127 = 129 and 219 = 221; Gallica 131 = 129 and 223 = 221. The rest follow from the
offsets.

## Use in M6a

- Pass 1 is from the archive.org offprint and pass 2 from the Lille offprint, transcribed
  independently.
- The Gallica CdT scan is the third witness. It resolves every pass 1 / pass 2 disagreement.
- The embedded OCR is not used for tables. Transcription is from page images only.

## Verification record

- **Double entry.** Pass 1 and pass 2 were transcribed by two separate AI agents. Each
  worked from one scan, from page images, under the brief in the session scratchpad. That
  gave 14 disagreeing cells, all resolved from the Gallica images
  (`leverrier_1846_resolution_log.csv`).
- **Spot check** (2026-11-06): 28 rows, 336 cells, 0 mismatches
  (`docs/review/spotcheck_results.json`).
  - The reviewer was an AI model in a chat session, against the Gallica images only.
  - Transcription, reconciliation and this first check were all done by the same model family.
- **Human spot check** (2026-11-07): the same 28 rows were spot-checked by a human against the
  Gallica scan (28 rows, 336 cells, 0 discrepancies)
  (`docs/review/spotcheck_human_results.json`).
  - The reviewer was Atishay Pahadiya (the dataset author). Every cell was ticked individually.
  - It was a verification against the image, not a blind re-transcription.
  - The other 251 rows have had no human check.
- **Misprints in the source** (`leverrier_1846_misprints.csv`) were each checked against the
  Gallica image by an AI model. Two misprint rows (129/10 and 135/3) are in the human sample, where
  they were checked as printed. No human has reviewed the misprint judgements themselves.
