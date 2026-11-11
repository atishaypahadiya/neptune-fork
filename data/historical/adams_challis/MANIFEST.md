# Adams–Challis study: source scans (Phase 1b, step 0)

Inventory date: 2026-11-10. The reading of each source is in `docs/notes/P1B_sources.md`.

Only one scan is committed here, `challis_1846_mnras_zenodo.pdf` (396 KB, CC0 on Zenodo). Every other
scan sits in `data/cache/adams_challis/` (git-ignored) and is identified below by its SHA-256.
Re-download from the URL and check with `shasum -a 256 <file>`. The `.txt` files beside them in the
cache are pypdf text extractions (archive.org `_djvu.txt` where noted) and are not sources.

Licence column: the 1846–47 texts are in the public domain. The scan providers' terms differ, which is
why most files are not redistributed from this repository:

- Zenodo records 1431903 and 1424621: CC0 1.0 (the record's licence field).
- ADS legacy scans ("© Royal Astronomical Society • Provided by the NASA Astrophysics Data System"):
  served by ADS with the publisher's permission; not redistributed here.
- archive.org Google scans: "possible-copyright-status: NOT_IN_COPYRIGHT"; Google's front matter asks
  for non-commercial use and attribution; not redistributed here (18 MB and 1.4 MB).
- Numdam: public-domain texts may be republished with the author's moral rights respected
  (numdam.org/reuse); bulk download and modification of the files are not allowed
  (numdam.org/conditions); the files carry the Gallica/Mathdoc cover page; not redistributed here.
- arXiv 1610.06424: arXiv non-exclusive licence (author's copyright; A&G version © RAS/OUP).
- arXiv 2108.06305: CC BY 4.0.
- Kollerstrom 2006 Part 1: © the author / SAGE (History of Science); the cached copy is from a
  third-party mirror; cache only.

## Primary printings

| file (cache) | source | printed pages | PDF pages | text layer | SHA-256 |
|---|---|---|---|---|---|
| `challis_1846_mnras_zenodo.pdf` (**committed**) | Zenodo 1431903, `article.pdf`, md5 861660a543ef180fc291dd2308965571; https://zenodo.org/records/1431903 | MNRAS 7 (9), 145–149 | 5 | yes | `ba2a100bd5e5f43687f4f1f83e15c1296dbe63950cc1d29e1c31bdac6f0c5121` |
| `ads_1846MNRAS...7..145C.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1846MNRAS...7..145C/ADS_PDF | MNRAS 7, 145–149 | 5 | yes | `a209fe102ef1aabf8ede28114557601f77f4a34572cfe0b43f3b591a5b6373ef` |
| `ads_1846MNRAS...7..149A.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1846MNRAS...7..149A/ADS_PDF | MNRAS 7, 149–152 | 4 | yes | `ba3cbb41ad73f2e7871ace9cb899432430d714b175e32fa45e7f6119e428f817` |
| `ads_1846MNRAS...7..121A.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1846MNRAS...7..121A/ADS_PDF | MNRAS 7, 121–144 (Airy's account; the whole item is the meeting report of 13 Nov 1846) | 24 | yes | `a53ac6482734791f28db90cdbbdbd082ab58a03055e7c3ae6610a71954720c1e` |
| `ads_1847MmRAS..16..385A.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1847MmRAS..16..385A/ADS_PDF | Memoirs RAS 16, 385–414 (Airy) | 30 | no | `e9b8f870b74eb405b3b6668bbaf3008606ea3551acc42ae88856cd9f80d30493` |
| `ads_1847MmRAS..16..415C.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1847MmRAS..16..415C/ADS_PDF | Memoirs RAS 16, 415–426 (Challis) | 12 | no | `83d02bf05ae65f1ce17d253271766891e3e827968aa06bb67a2e66b770fc28ff` |
| `ads_1847MmRAS..16..427A.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1847MmRAS..16..427A/ADS_PDF | Memoirs RAS 16, 427–460 (Adams) | 34 | no | `4ac10e1def2e7cdb69c100952923c4fd543f8270941d05618c326a422d31e0c9` |
| `memoirs_ras_16_1847_archive_org.pdf` | archive.org `memoirsroyalast01socigoog` (Google scan of Memoirs RAS vol. 16, 1847; source books.google.com id MqsRAAAAYAAJ); https://archive.org/download/memoirsroyalast01socigoog/memoirsroyalast01socigoog.pdf | whole volume; printed 385 = PDF 420, 415 = PDF 450, 427 = PDF 464 (checked on the page images; leaf 462 repeats printed 426) | 643 | no (OCR in `_djvu.txt`, SHA-256 `8d7351f1afc2fb48529716d9a7a8be30d8c2eeee988db23d7b4c0c18f3cba899`) | `5260eae02aa78692e3623e6671949996e7e2242f048118a554fdebb2cdd4ded2` |
| `airy_1847_an_reprint_zenodo_1424621.pdf` | Zenodo 1424621, `article.pdf`, md5 67f76508e28d261aaa2bc64790313df1; https://zenodo.org/records/1424621 (DOI 10.1002/asna.18470251002; the record's date "1970-01-01" is a placeholder) | Astronomische Nachrichten 25, no. 585 (1847): reprint of Airy's account | 14 | yes | `6c36020bb91944fb3d13571e38947e4805fcc6ab29d78edbde5261ab0c239591` |
| `adams_1846_nautical_almanac_appendix_archive_org.pdf` | archive.org `explanationobse00Adam` (Smithsonian Libraries; "Not in copyright"); https://archive.org/details/explanationobse00Adam | separate printing, London: W. Clowes & Sons, 1846, "From the Appendix to the Nautical Almanac for the year 1851", own pagination 1–31; §§ 1–61 | 44 | no (OCR in `_djvu.txt`, SHA-256 `7cfd0aee7274e0d856629c32a06696aa2e3fae7799c579a4dfa29a0d03228ca5`) | `451b7668369911d6241128616a4dc6a108320e84a0b34bf9f4fbbfd6c6e873ec` |
| `adams_1876_jmpa_part1_numdam.pdf` | Numdam `JMPA_1876_3_2__5_0`; https://www.numdam.org/item/JMPA_1876_3_2__5_0.pdf | J. Math. Pures Appl. 3e sér. 2 (1876), 5–32 (French translation, part 1, with an *avertissement*) | 29 | yes | `6d5b98090fe9e747eb75a32a2226b2699ba5c06d4fa4b0ab05e95f52fb94958b` |
| `adams_1876_jmpa_part2_numdam.pdf` | Numdam `JMPA_1876_3_2__69_0`; https://www.numdam.org/item/JMPA_1876_3_2__69_0.pdf | J. Math. Pures Appl. 3e sér. 2 (1876), 69–86 (part 2, *suite*) | 19 | yes | `812dae096e9ac5cc6930b6444aedca6cdc6f4c8952d8d7cdaadba6044c4bb77a` |
| `ads_1904MmRAS..54..143S.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1904MmRAS..54..143S/ADS_PDF | Memoirs RAS 54 (1904), 143 ff. (48 scanned pages including 20 plates): R. A. Sampson, description of Adams's manuscripts on the perturbations of Uranus (see the note for the check of this identification) | 48 | no | `2cf5fba128ed0ca74b61c643ed456f62fdb8306f0cf1ca871c3b974e0409e27c` |

## Secondary works (cross-checking only)

| file (cache) | source | SHA-256 |
|---|---|---|
| `krajnovic_2016_arxiv_1610.06424v2.pdf` | arXiv:1610.06424v2 (published as A&G 57 (5), 5.28–5.34, DOI 10.1093/astrogeo/atw183) | `0b07a6c4fc01d085bef2bbde7089523a3a5bf956761b8195e373f8696fe8009b` |
| `krajnovic_2021_arxiv_2108.06305.pdf` | arXiv:2108.06305 (chapter in *The Discovery of Neptune*, Springer 2021, DOI 10.1007/978-3-030-54218-4_6), CC BY 4.0 | `cd62724e87aa18593b67dc8c78d9cfcc68329aa542952f5c404f625629bde744` |
| `kollerstrom_2006_hiatus_part1_columbuslandfall.pdf` | History of Science 44 (1), 1–28, DOI 10.1177/007327530604400101; cached from a third-party mirror (columbuslandfall.com, expired certificate) | `ae54d89208021641674694e46eb86b7e2e8148d8ab93efc87eaba7bc3b3b82f0` |

Also in the cache, not sources: saved catalogue pages (`cul_*.html`, `cudl_rgo6.html`, `numdam_jmpa_*.html`,
`hathi_000552896.html`), `memoirs_scandata.xml` (archive.org leaf-to-page map used for the mapping above), and
low-resolution page renders used for the page checks (`memoirs_pdfpage_*.png`, `pamphlet_page_*.png`,
`sampson_p*.png`).

## Derived tables (committed, step 1)

- `adams_solutions.csv`: one row per dated solution (S1–S8), every value primary unless its cell says
  derived; sources and page references in the row.
- `adams_solutions_statements.csv`: the 44 printed statements the rows rest on, each with its page
  image name (renders in `data/cache/adams_challis/renders/`, made from the PDFs above with pymupdf).
- Transcription passes and reconciliation: `docs/notes/adams_solutions_pass1.csv`,
  `docs/notes/adams_solutions_pass2.csv`, `docs/notes/adams_solutions_reconciliation.md`.

## Derived tables (committed, step 2)

- `challis_sweeps.csv`: one row per observing session with the plan and the total; verdicts on
  Neptune's presence marked D where they rest on the computation.
- `challis_sweeps_statements.csv`: the 111 printed statements behind them, with page images.
- `neptune_on_challis_dates.csv`, `neptune_track_1846.csv`: computed from the cached Horizons states
  by `docs/notes/scripts/p1b_neptune_on_challis_dates.py` (no fetch).
- Passes and reconciliation: `docs/notes/challis_sweeps_pass1.csv`, `challis_sweeps_pass2.csv`,
  `challis_sweeps_reconciliation.md`.
