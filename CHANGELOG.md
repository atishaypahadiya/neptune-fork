# Changelog

Versions follow semver. Dates are the release dates; the DOIs are in `CITATION.cff`.

## Unreleased

Everything on `main` since 1.2.1. No number already on either page changes; the companion
page gains § 3.5 and the text listed under Changed.

### Security

- soupsieve is 4.1, for advisory GHSA-mmhr-3c5m-f64f (a denial of service in selector
  matching). It comes in through astroquery. No patched soupsieve supports Python 3.11.0
  to 3.11.4, so the minimum Python is now 3.11.5.

### Changed

- The companion page: a new § 3.5 with Table 2 (the four cases), the pre-registered
  reading, the post hoc variants and an interpretation labelled as such; the two findings
  on Bouvard's tables and the manuscript check in § 2 (values only, no image); an
  acknowledgement; two references. The abstract, the end of § 3.4 and one limitation no
  longer say that Bouvard's tables were not used.

### Added

- The companion study's attribution step, pre-registered: Adams's four cases (his
  observations or noise-free ones, against Bouvard's reference or a clean one) by both
  methods and both hypotheses. Reading: "depends on the method"; the shift is not
  attributed. Post hoc variants of the 1780 place point the other way from the case as
  run (`docs/notes/P1B_bouvard_attribution.md`). The rows of Bouvard's tables of 1821
  that his places touch are transcribed in two entries, with a misprint log
  (`data/historical/adams_challis/bouvard_1821/`).
- Two findings on Bouvard's tables: the leading coefficient of the equation of the centre
  on p. xv differs from the tables in one digit, and Table I's mean longitude falls behind
  the introduction's mean motion, as a handwritten line in the Harvard copy also says
  (`docs/notes/P1B_bouvard_1821.md`).
- The manuscript check of Adams's solutions (St John's College, Cambridge, MS W16) has a
  spot check by the author: 20 cells read against the images, all matching. The two
  readings the entries left doubtful are resolved. 11 of the 20 cells were chosen and 9
  were random, so the check confirms every disputed reading but bounds the cell error rate
  only to about 34 % (`docs/notes/P1B_adams_solutions.md`).
- For a later study of Mercury's perihelion, the first data set: Newcomb's observed times
  of contact for the transits of 1677 to 1881 (1882, pp. 451–452), 58 rows in two entries
  that agree in all 696 cells, with a misprint log (`data/historical/mercury/`). The
  author's spot check of 20 cells, 7 chosen and 13 random, found no discrepancy; the random
  part bounds the cell error rate to about 25 %. Nothing is fitted to these data. A
  blinding record says what has already been seen (`docs/phase2/blinding.md`); the
  tabular times and the observed-minus-computed table are not in the repository.
- Scoping of that study, with nothing frozen: a survey of prior art that tags each
  statement by how it was checked, an inventory of the observations, a design sketch and a
  draft list of milestones (`docs/phase2/`). The survey is not a checked literature review.
- The forces that study needs, all off by default: the Sun's post-Newtonian term, the
  Sun's oblateness, a circular ring and a planet inside Mercury's orbit, with a forward
  model to Mercury's longitude. The validation reports five checks: the relativistic advance
  of 42.98″ per century is reproduced within 0.00001″, the oblateness and interior-mass rates agree
  with their formulas, and against the JPL ephemeris of 1800 to 2000 the longitude stays
  within 0.054″ with the term and drifts by 43.01″ per century in the perihelion without
  it (`docs/phase2/validation.md`). No existing source file is edited.
- A tool that generates Mercury's longitude under each competing model and reports how
  the predictions differ. Built and tested only; it reads no observation.
- Draft plans for the next release's runs, with their costs, and the code they need,
  exercised on synthetic data only (`docs/notes/V13_drafts.md`, `fork.export.v13`). No
  draft is frozen and none has been run on Le Verrier's table.
- Lalande's two places of Neptune of May 1795, in two entries with a misprint log
  (`data/historical/precoveries/`). Galileo's are not sourced.
- A spot-check kit for the companion study's tables: the script, 42 drawn rows and a
  results template (`docs/review/`).
- The README links the site and the companion study under its title.

### Changed

- The weekly link check also reads `docs/phase2/`.

## 1.2.1 (2026-11-13)

### Corrected

- The companion page's chart case. v1.2.0 printed two rows for the Berlin sheet Challis had,
  0.00 with its limits for the equinox of 1800 and 0.92 (0.77 to 1.00) with the same limits
  applied to right ascension of date, and read the case as depending on the sheet's limits.
  The second row came from a fallback applied in the wrong frame and is withdrawn. Under the
  sheet's own limits (21h 56m to 23h 04m, 15° south to 15° north, for 1800; read from the
  sheet and the Academy's rules of 1859, and checked against the sheet images) the planet
  lay off the sheet on every night Challis recorded it: 0.00, would very probably not have.
  Table 1, Figs. 3 and 4, the limitations and the write-up follow. Fig. 3 loses the second
  edge line, the sheet's edge read as right ascension of date, with the row.
- Krajnović 2021 is a chapter of *Neptune: From Grand Discovery to a World Revealed*
  (Springer), pp. 185–244, doi:10.1007/978-3-030-54218-4_6. v1.2.0 gave another title.
- The link to the Lille scan of Le Verrier's memoir: its handle no longer resolves, and the
  notes now give the file's address. The released dataset is unchanged.
- Public text names the studies descriptively and is written in the first person singular.

### Added

- A check of Adams's solutions against his manuscripts (St John's College, Cambridge,
  MS W16): 11 values agree with the printed sources, 9 are new, none contradicts
  (`docs/notes/P1B_adams_solutions.md`, amendment of 12 November 2026). No image is in the
  repository.
- Bouvard's tables of 1821: the copy, its structure and the printed formula for Uranus in
  two entries (`docs/notes/P1B_bouvard_1821.md`). No place is computed from it.
- A weekly link check, validation of the citation metadata, and tests of the public text.

### Changed

- Both pages are laid out on one grid of three widths.

## 1.2.0 (2026-11-12)

A companion study of the British search for the planet, read from the printed record
alone. It asks two questions.

How good was Adams? At the distances he assumed, the fit to his data choice is 10.5°
and 9.1° short of Neptune at epoch 1846.762 (5 October 1846). His two solutions are reproduced from his memoir. On
a noise-free Neptune at his first distance every method is short (his least squares 9.5°, the
model on his places 7.8°, the fit 10.4°), so the shortfall belongs to the
assumed distance, and his nearness to the truth (+1.7° and +2.5°) to his data or their
reduction. Whether noise alone accounts for that nearness is not settled: the pre-registered
statistic measured jumps to a second minimum of his least-squares problem, and the post hoc
statistic (0.017 under his first hypothesis) reads weak.

Could Challis have recognised the planet before 23 September 1846? A detection model
conditioned on his record gives, with the realised events held fixed: his actual procedure 0 by
construction (the cases measure what a change would have added); the comparison he broke off,
carried through, 0.71 (0.50 to 0.97); zones laid from Adams's first place or from Le Verrier's June place 0.41 and 0.41 (both prior-dependent), or 0.79 and 0.85 with prompt comparison. Each figure depends on an
input the record does not fix, three readings are prior-dependent, and none reaches "very
probably".

- `docs/notes/P1B_*.md`: sources, Adams's solutions, Challis's search log, the fit at his
  distances, the reproduction, noise against a systematic cause, the detection model, and the
  write-up with each claim's JSON path and deviation (`docs/notes/P1B_writeup.md`). Each
  pre-registration stands unedited above its results.
- `data/historical/adams_challis/`: the transcribed tables of the study (CC0). They are in the
  code archive only; they have had no human spot check and are not a dataset record.
- `fork.historical.adams_1846`, `fork.historical.challis_1846`; `fork.export.p1b_adams`,
  `p1b_adams_repro`, `p1b_adams_noise`, `p1b_challis`, `p1b_site`; their tests.
- Site: the companion page `adams-challis/`, with Table 1 and Figs. 1 to 4, every number a
  named field of a committed JSON; the main page gains one sentence at the end of its abstract.
  The prerender, the number audit and the layout check cover both pages.
- `docs/reproduce.md`: sections for the headline fit and the sensitivity runs, the
  repeatability runs, the site's data, the fast-likelihood validation, its reruns and the
  full-likelihood check of the eccentricity test, and the companion study.
- Wording: the notes, the source manifest and the documentation describe the transcription and
  its checks by their method. The caveats stand: both passes and the first spot check share one
  method, so their errors may be correlated, and 251 rows of Table 77 have had no human check.
- Dataset v1.0.1 (concept DOI 10.5281/zenodo.24044791): rewords descriptive strings in
  `spotcheck_results.json` and the documentation; no table value or check result changes.
- Compute ledger: 1193 CPU-hours (the main study: 1090; the companion study: 103).
- Deferred from v1.2 to v1.3: the Peirce shape-model comparison, the held-out
  precoveries, the remaining variants of the headline fit (the Bouvard and Adams data variants,
  the 1846 reduction track, the Le Verrier shape model), and the resolution of the eccentricity
  test's low-distance cell.

## 1.1.0 (2026-11-10)

The distance claim of v1.0.0 held only under the Gaussian error model. v1.0.0 said that the
pre-1781 observations place the planet at 28.0–34.0 AU and prefer Neptune's distance over
Le Verrier's 36 AU by Δ ln L = 3.15. A second pre-registered profile likelihood in the
semi-major axis repeats the Gaussian-error profile under the Student-t errors (ν = 4) that the
evidence favours,
on six data sets. Its reading is "partial" and draws no conclusion by the rule; what the curves
show is that the near side is the same under both error models (Δ ln L = −2 at 28.0 and
28.7 AU) while the far side depends on it (34.0 AU under Gaussian errors, 40.6 AU under
Student-t, with Le Verrier's 36.15 AU at −0.76 instead of −3.15), and that the near side rests on
Flamsteed's 10 measurements (without them the profile reads as flat; without Bradley's, Mayer's
or Le Monnier's it is almost unchanged). The headline longitude, mass and semi-major axis are
unchanged.

- `docs/notes/M10.md`: the pre-registration, its dated amendment, the code hashes, the relaunch
  note and the results in the pre-registered order; `docs/notes/M10_profile.{json,png}`.
- `fork.export.m10`, `fork.export.m10_plot`, `tests/test_m10.py`; the Uranus-only box is centred
  on the better of two starts; `FORK_HORIZONS_CACHE_ONLY=1` makes an ephemeris cache miss raise
  instead of fetching.
- Site and README: the abstract's distance sentences replaced; § 5 gains the Student-t
  paragraph, Fig. 7 shows both error models side by side, and Fig. 9 shows the four
  leave-one-out profiles; two Limitations bullets replaced (the far side depends on the error
  model; the near side rests on Flamsteed's measurements). Every number is read from
  `site/data/m10_profile.json`.
- Compute ledger: 1090 CPU-hours (the profile under Student-t errors: 3.1).
- Code DOI 10.5281/zenodo.24122940 (concept 10.5281/zenodo.24121059).
- Deferred from v1.1 to v1.2: the Peirce shape-model comparison, the held-out
  precoveries, the Bouvard and Adams data variants, the 1846 reduction track, the
  Le Verrier shape model, and the resolution of the eccentricity test's low-distance cell.

## 1.0.0 (2026-11-10)

First release: the main study through the headline fit on the historical data. The pre-registered fit to Le Verrier's 1846 Table 77 places
the unseen planet at 325.8° on the night Galle found it (68 %: 323.6°–327.8°; Neptune 325.87°),
with the 95 % range provisional and the eccentricity test unresolved; the Gaussian profile
likelihood in the semi-major axis; the site as a short paper with plates; the transcribed
Table 77 as a separate CC0 dataset (concept DOI 10.5281/zenodo.24044791). Code DOI
10.5281/zenodo.24121060 (concept 10.5281/zenodo.24121059).
