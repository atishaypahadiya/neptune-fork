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

## Primary printings (added in step 3, amendment 2)

Gallica refused every request on 2026-11-10 (HTTP 403 on the SRU search, the periodical's date
page and the IIIF manifests, four identifiers tried), so the volume was taken from archive.org
instead, once, online, before the launch.

| file (cache) | source | printed pages | PDF pages | text layer | SHA-256 |
|---|---|---|---|---|---|
| `cr22/comptesrendusheb22acad.pdf` | archive.org `comptesrendusheb22acad` (Comptes rendus hebdomadaires des séances de l'Académie des sciences, t. 22, janvier–juin 1846; contributor Natural History Museum Library, London; no licence field in the item metadata, an 1846 printing); https://archive.org/download/comptesrendusheb22acad/comptesrendusheb22acad.pdf | Le Verrier, "Recherches sur les mouvements d'Uranus", séance du 1er juin 1846, 907–918; the longitude statement at 917 (PDF 931): "en assignant 325 degrés de longitude héliocentrique à la planète, au 1er janvier 1847, on ne commet pas une erreur de 10 degrés", and the formula ν = 314°,5 + 12°,25 α + (1/m){20°,82 − 10°,79 α − 1°,14 α²} | whole volume, 1228; printed 907 = PDF 921 (checked on the page image) | yes (OCR) | `555936e73efaca9448b38ecd9c7f4aa62f405d232c44e582a637f92aa693dac9` |

The LV-June target of `docs/notes/P1B_adams_fit.md` therefore uses the printed value, marked P;
Airy's summary (MNRAS 7, 131–132) is no longer the fallback. Page renders in
`data/cache/adams_challis/renders/cr22_p921.png` … `cr22_p933.png`.

## Secondary works (added in step 3b, stage A, 2026-11-11)

Sought once each, as the step-3b plan requires, before any stage-B code.

| file (cache) | source | printed pages | PDF pages | text layer | SHA-256 |
|---|---|---|---|---|---|
| `ads_1970CeMec...3...67B.pdf` | https://ui.adsabs.harvard.edu/link_gateway/1970CeMec...3...67B/ADS_PDF (served from articles.adsabs.harvard.edu; ADS scan, publisher's copyright notice on every page, not redistributed) | C. J. Brookes, "On the prediction of Neptune", Celestial Mechanics 3 (1970), 67–80 (DOI 10.1007/BF01230433) | 14 | no (image scan; read from renders) | `054b6285a75b0b82803d405bb88e520c5d1a7a7aabcc25511f4f084af5d8da9e` |

Not obtained: W. M. Smart, "John Couch Adams and the discovery of Neptune", Occasional Notes of the
RAS No. 11 (vol. 2, 33–38 as Brookes cites it; a 56-page pamphlet, 1947). The ADS gateway has no
entry under the guessed bibcode (HTTP 404) and Nature's page for the 9 November 1946 abridgement
(Nature 158, 648) serves only a paywall landing page; the pamphlet exists in print only (St John's
College archives list a copy). Brookes 1970 acknowledges Smart for suggesting the problem.

What Brookes 1970 contains (read in full on 2026-11-11, the same day it was obtained): a restatement
of Adams's linearised model (his eqs 2.1–2.4, Tisserand's form of F_i); a re-derivation of the
Hypothesis I solution from the memoir's equations giving e′ 0.1509, m′ 0.8267, ϖ′ 311° 35′,
l′ 323° 39′ against Adams's 0.1610, 0.8282, 315° 57′, 325° 4′ [sic], the true longitude at
discovery 328° 4′ against 328° 42′, "due mainly to a small numerical error occurring in one of the
equations derived by Adams" (his 1969 thesis; the equation is not named in the paper); Hypothesis II
by Hypothesis I's procedure (0.1387, 0.7393, 303° 52′, 324° 47′, 331° 31′ against Adams's 0.1206,
0.7502, 299° 11′, 322° 58′, 329° 25′); the effect of Flamsteed 1690 (Table III); solutions with 61,
31 and 21 equations (Table V); a first-order-only variant (p. 75); solutions at α = sin 29°, 30°, 31°,
32° with the same theory (Table VI, X) and at sin 34°, 35°, 37° with F₄, G, H terms added
(Tables X–XIV; a′ = 30 excluded for the near-commensurability of n and n′); the conclusion that the
true longitude at discovery is insensitive to a′ between 30 and 40 AU for any epoch after 1830, and
that the mean longitude for 1846 rises then falls with decreasing a′ (323.1, 324.3, 325.5, 326.7,
321.4, 322.3, 322.0 for a′ 39.56 … 31.86), against Adams's § 54 expectation.

## Derived tables (committed, step 3b stage A)

Transcription of Adams's memoir §§ 8–55 (pamphlet, PDF 12–37) in two independent passes,
reconciled in `docs/notes/adams_memoir_reconciliation.md` (1988 common values, 15 disagreements,
all decided on the named renders; both passes kept under `docs/notes/adams_memoir_passes/`).

- `adams_memoir_normal_places.csv` (30 rows): §§ 8–9 residuals, §§ 51–52 theory and residuals.
- `adams_memoir_series.csv` (140 rows): §§ 12–13, 22, 28, 33–34, 45 constants and coefficients.
- `adams_memoir_equations.csv` (1407 rows): every printed equation of §§ 15–27, 29, 35–46.
- `adams_memoir_solution.csv` (104 rows): §§ 30–32, 47–50, 53, 55.
- `adams_memoir_statements.csv` (56 rows): the printed blocks behind them, with page numbers.

Renders used for the decisions: `data/cache/adams_challis/renders/pamphlet_hires/` (220 dpi
halves of PDF 18, 24–37 from the PDF's embedded scans, and 600 dpi crops `z018_*.png`,
`z024_formula.png`), regenerated by `uv run --with pymupdf python` from the cached PDF; not
committed.

## Named input N1, the Berlin sheet (fetched 2026-11-12; cached, not committed)

Under `data/cache/adams_challis/berlin/`. Used in `docs/notes/P1B_challis_detection.md`,
"Named input N1: search of 12 November 2026".

| file | source | SHA-256 |
|---|---|---|
| `krajnovic_2021_arxiv_2108.06305.pdf` | https://arxiv.org/pdf/2108.06305 (Fig. 6.10 on p. 28: the lower part of Hora XXII, Potsdam copy) | `cd62724e87aa18593b67dc8c78d9cfcc68329aa542952f5c404f625629bde744` |
| `jones_2002_highlights.pdf` | Jones 2002, DOI 10.1017/S1539299600013770 | `87b3f2c80bd968a667da22df56333bf0ef84c9cdb88741e0118bc28fcb7b52c2` |
| `glam_IAS_K178_58271_manifest.json` | https://glam.uwr.edu.pl/iiif/IAS_K178_58271/manifest (the 1859 catalogue, Wrocław University, shelfmark K178, 598 images) | `f06d6cde433c3564bf9149c25fb792d2cd20f550fde6a3b03b9a73dd5a89d34a` |
| `wroclaw/p007_w2600.jpg` | image 7 of that manifest at width 2600 (title page) | `434352692928ef5c812d2640f4794afbd20c610b9d84d48c9c843ca87130a0b5` |
| `wroclaw/p010_w2600.jpg` | image 10 (Vorbericht p. IV) | `34de3b35c1ae1dccb1fe662f20e8597f14626e00b713a3140214187fbceea823` |
| `wroclaw/p012_w2600.jpg` | image 12 (p. VI) | `7a2bc356b6e110d2e8ccdaa8b416aa03959afc966cc1fb6495a9381d51e85d84` |
| `wroclaw/p013_w2600.jpg` | image 13 (p. VII) | `d73e6f9c9ddf74e653dea86156a33ab52f653373007b07e99fac4db9527fed12` |
| `an_bub_gb_v1Q_AAAAcAAJ_djvu.txt` | archive.org `bub_gb_v1Q_AAAAcAAJ`, text layer; it is *Astronomische Nachrichten* 5 (1827), not the volume 4 wanted | `3b12add4d1967c329e5e40ecd54c1f7b4e2209aa51ccd6329399f07fe96f20b9` |

## Bouvard's 1821 tables: source search of 12 November 2026 (located, not obtained)

Bouvard, A., *Tables astronomiques publiées par le Bureau des longitudes de France, contenant
les tables de Jupiter, de Saturne et d'Uranus, construites d'après la théorie de la Mécanique
céleste* (Paris: Bachelier et Huzard, 1821; XXVIII + 110 pp. by the catalogue record of the
French defence libraries). These are the reference Adams subtracted (his § 7; see
`docs/notes/P1B_adams_reproduction.md`, the last named input).

No file was downloaded, so there is no hash and nothing in the cache. Every copy below
refused a scripted request once and was not asked again. Each needs a browser.

| copy | URL | answer to one scripted request (12 Nov 2026) |
|---|---|---|
| Gallica (BnF) | https://gallica.bnf.fr/ark:/12148/bpt6k9668766k | catalogue search (SRU) answered 403; the ark itself was not requested |
| Google Books, Harvard copy | https://books.google.com/books?id=SfZ4Eiw9yR0C | 429, a browser check page |
| Google Books, Austrian National Library copy | https://books.google.com/books?id=LPRVAAAAcAAJ | 429, a browser check page |
| Google Books, Oxford copy | https://books.google.com/books?id=E2xbAAAAQAAJ | 429, a browser check page |
| Google Books, a fourth record | https://books.google.com/books?id=eMFGAAAAcAAJ | 429, a browser check page |
| French defence libraries, record 77d0faf0-47b0-485f-8c1d-121335baa70f | https://bibliotheques-numeriques.defense.gouv.fr/document/77d0faf0-47b0-485f-8c1d-121335baa70f | no connection (timed out twice) |

Searched without a copy: archive.org (advanced search, eight queries: only the 1808 edition,
`bub_gb_Yuu0Iel8Jc8C`, Jupiter and Saturn in the decimal division, which has no Uranus
tables); HathiTrust (catalogue search answered 403 with a browser check); e-rara (browser
check); ADS (the search interface needs a token; a book of 1821 is not expected there); the
Munich digital library (its search needs a browser). The Google Books volumes API was over
its daily quota. That the four Google Books records and the Gallica ark are full scans of
the 1821 volume rests on a web search's result list; it has not been seen on a page.

Not done, because no page was seen: rendering the Uranus tables, recording their structure
(arguments, epochs, perturbation terms, masses) and estimating the transcription effort for
Adams's 30 normal-place epochs. Once a PDF is saved by hand, put it under
`data/cache/adams_challis/bouvard/` and add its SHA-256 here.

## St John's College, Cambridge, MS W16 (received 2026-11-12; cached, not committed)

Source: St John's College, Cambridge, MS W16 (Adams's papers on the perturbations of Uranus,
in the arrangement of Sampson 1901). 28 image files of seven sheets, supplied by the
College's Special Collections and received on 12 November 2026. Each sheet came in four
sizes. Under `data/cache/adams_challis/stjohns_W16/` (ignored by git).

**Rights condition.** The images may not be committed, put on the site or published without
prior notice to St John's College. Any reproduction carries the line "By permission of the
Master and Fellows of St John's College, Cambridge". Values transcribed from the images may
appear in the notes and tables.

| file | pixels | what it shows | date on the page | shelfmark stamp | SHA-256 |
|---|---|---|---|---|---|
| `16867.jpg` | 2695 × 4373 | the note of results: elements and mass of the supposed planet, signed J. C. Adams (the sheet behind S4) | undated (Sampson and the facsimile date it 18 September 1845) | E III 14 | `3ac508e20af6f0f745c66029ad87ba17011b703c470db1db11eeded1054dd6f6` |
| `16867_Large (1500 x 1500).jpg` | 924 × 1499 | the supplier's smaller copy of 16867.jpg |  |  | `ebcde29f054a6f55989d3308aace3e3c6b1c15c2bde07bf6769499576e4f2d9f` |
| `16867_Medium (1000 x 1000).jpg` | 616 × 1000 | the supplier's smaller copy of 16867.jpg |  |  | `b1d3e9574ee7e117497f88889c83a25c33355ee24ec321bd64406d7f33c8b2c6` |
| `16867_Small (500 x 500).jpg` | 308 × 500 | the supplier's smaller copy of 16867.jpg |  |  | `64c7898ef2a8b6b09ec5486493bfd8a42f43fe39514a205b464cdac44e1efa86` |
| `16868.jpg` | 5110 × 6636 | calculation for two trial values of θ (−45°, −43°) ending in θ very nearly −42° 51′ | August 20th 1846 | E VII 15 | `7b0e988c0b87ee57c5a85e2db14e64d3ca597c8082e1c165d095214fecdd8f6f` |
| `16868_Large (1500 x 1500).jpg` | 1155 × 1500 | the supplier's smaller copy of 16868.jpg |  |  | `6194d71a2d32b3b9e3a9f3064b3f7cf8201370f8bbb119b25b6f124a9e3bea4f` |
| `16868_Medium (1000 x 1000).jpg` | 770 × 1000 | the supplier's smaller copy of 16868.jpg |  |  | `063e7f1d5544f7e653a765972736eb6ad1a3f338ae34ae8749bfa04090b61c0a` |
| `16868_Small (500 x 500).jpg` | 385 × 500 | the supplier's smaller copy of 16868.jpg |  |  | `76b09804c42c9c626deceefba2356736218aea5406f43d9d3ff0abe69d847373` |
| `16869.jpg` | 5087 × 6561 | calculation of the additional perturbation of the radius vector for 1834, 1840 and 1846 | 1st Septr 1846 | F II 6 | `b58a61f90e275fadc83af90ca95fbcec59eda598d90d0d1ec5d361a85a24ecc4` |
| `16869_Large (1500 x 1500).jpg` | 1163 × 1500 | the supplier's smaller copy of 16869.jpg |  |  | `d827c3cba4aab3ba37397054f6529f44cd5c27ef30e31cea18a7c2650ba54524` |
| `16869_Medium (1000 x 1000).jpg` | 775 × 1000 | the supplier's smaller copy of 16869.jpg |  |  | `ae4716e078b9adb76ce4f1738e460e3952fcdecda1bd1baf72cd780cf032e6a1` |
| `16869_Small (500 x 500).jpg` | 387 × 499 | the supplier's smaller copy of 16869.jpg |  |  | `8ffec92f3b68c11dfcb30ebe2fa1f6e96e109088fc15c0875e333b328c4f652e` |
| `16870.jpg` | 4978 × 6456 | first approximation, θ = −51° 10′, ending in θ = −50° 34′ to be substituted | Septr 18th 1845 | E III 6 | `42ba3345cce776feb6a08872f6791b97654b97fa16ca064e3d980df3c2cee786` |
| `16870_Large (1500 x 1500).jpg` | 1156 × 1499 | the supplier's smaller copy of 16870.jpg |  |  | `614265c3e9b5cab8c790599a2b26a1d23c624181c3d42b435e036155293bd1b9` |
| `16870_Medium (1000 x 1000).jpg` | 771 × 1000 | the supplier's smaller copy of 16870.jpg |  |  | `74521c8ead975100a820ef134545a91e3f814c8bd29c864bc338d00e56e6f225` |
| `16870_Small (500 x 500).jpg` | 385 × 499 | the supplier's smaller copy of 16870.jpg |  |  | `88e874eac7a312965f127eca8b601836d0ff12a15a5e7efe47f6e644e97d99f1` |
| `16871.jpg` | 5141 × 6452 | calculation of errors for 1843 and 1846, with a calculation of the radius vector | undated | not visible on the image | `675a3339ef98e72b3a82c96b8deb5fe6676009dfe34df71f89fb97aceb3bb60c` |
| `16871_Large (1500 x 1500).jpg` | 1195 × 1500 | the supplier's smaller copy of 16871.jpg |  |  | `9cf688237561f2f7ac768e83ab941efe8530b5b91c65db517e6661cec21e17e0` |
| `16871_Medium (1000 x 1000).jpg` | 796 × 999 | the supplier's smaller copy of 16871.jpg |  |  | `9c4d0cc82acae54ec64a9bc636353dcaf84c0bc0c5eb26e029b2b14c296a7274` |
| `16871_Small (500 x 500).jpg` | 398 × 499 | the supplier's smaller copy of 16871.jpg |  |  | `c946a5322bd8afd8cabdc79365a65f38b805bed4f5259b2f8f655201a20c0e13` |
| `16872.jpg` | 3764 × 2579 | printed facsimile of the memorandum of 3 July 1841 and of the note of results | printed heading: 1841 July 3 and 1845 September 18 | reproduces the stamp E III 14 | `67556f81503a10c218026cfb837ca5db162ceb5c3f8395980f79f35c6dd10317` |
| `16872_Large (1500 x 1500).jpg` | 1499 × 1027 | the supplier's smaller copy of 16872.jpg |  |  | `0dae1d5b900ba03d9275dd4bc091a70a7da87bc9135a0c165c48f0960c9a43da` |
| `16872_Medium (1000 x 1000).jpg` | 1000 × 685 | the supplier's smaller copy of 16872.jpg |  |  | `4630b539316e9a093cba61084ff8e19d947dedaf664a68ea9dda533da222da9d` |
| `16872_Small (500 x 500).jpg` | 499 × 342 | the supplier's smaller copy of 16872.jpg |  |  | `4979e4dee9346ef92b622c66e4b676bde46f378fc546aece91ba15395a1573b8` |
| `16873.jpg` | 3837 × 5132 | letter from Robert Main, Royal Observatory, Greenwich, to Adams with a table of the errors of Uranus for six dates of 1845 | 1846 Septr 7 | A VII 3 | `93a88ddd0ac355ed2d2a064da5f7e049292fa45b800a5cd1fe0ac498ef1c72e1` |
| `16873_Large (1500 x 1500).jpg` | 1121 × 1499 | the supplier's smaller copy of 16873.jpg |  |  | `f1ef118e86a3238362d3bf631d4806a8aefdb71ffb542ba0777c9d1cd35e5318` |
| `16873_Medium (1000 x 1000).jpg` | 747 × 999 | the supplier's smaller copy of 16873.jpg |  |  | `9dabae92033881c0cee66befaab199e4e7dfff21a25cb04c539c0bb79bfe1601` |
| `16873_Small (500 x 500).jpg` | 373 × 499 | the supplier's smaller copy of 16873.jpg |  |  | `379c7aa6d26c76dfd28c20f4a6f15402cc0d90c7c7918c465be9c44e9bb50dab` |

Transcription (13 November 2026 amendment of `docs/notes/P1B_adams_solutions.md`): two
entries of 200 cells from 17 bands cut from the full-size images (under
`data/cache/adams_challis/stjohns_bands/`, not committed).

| file | content | SHA-256 |
|---|---|---|
| `adams_solutions_ms.csv` | 20 manuscript values set beside the committed ones | `d2dd0e07d2b4ad8a9d91c28e07d65d026b6c051f6889123a9a379205eefeadf5` |
| `adams_solutions_ms_reconciliation.csv` | the four differences between the entries and the recut | `dff231d49845255e406ac8d875007f2590b9d213e64020f4c78c7efff0c63bca` |
| `passes_ms/stjohns_w16_pass1.csv` | first entry | `3de30e11d667b66741477395b6b2e60efe09e88fb485ce6944a4c812abb4c3bd` |
| `passes_ms/stjohns_w16_pass2.csv` | second entry, from the bands and a list of labels alone | `4f8184c4515d1ec304156e95eeaba52a34fc36c73add2e845fd0a233701381ee` |
