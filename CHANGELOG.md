# Changelog

Versions follow semver. Dates are the release dates; the DOIs are in `CITATION.cff`.

## 1.2.0 (2026-11-12)

Phase 1b, a companion study of the British search for the planet, read from the printed record
alone. It asks two questions.

How good was Adams? At the distances he assumed, the project's fit to his data choice is 10.5°
and 9.1° short of Neptune on Galle's night. His two solutions are reproduced from his memoir. On
a noise-free Neptune at his first distance every method is short (his least squares 9.5°, the
project's model on his places 7.8°, the project's fit 10.4°), so the shortfall belongs to the
assumed distance, and his nearness to the truth (+1.7° and +2.5°) to his data or their
reduction. Whether noise alone accounts for that nearness is not settled: the pre-registered
statistic measured jumps to a second minimum of his least-squares problem, and the post hoc
statistic (0.017 under his first hypothesis) reads weak.

Could Challis have recognised the planet before 23 September 1846? A detection model
conditioned on his record gives, with the realised events held fixed: his actual procedure 0 by
construction (the cases measure what a change would have added); the comparison he broke off,
carried through, 0.71 (0.50 to 0.97); zones from Adams's first place 0.41, or 0.79 with prompt
comparison; from Le Verrier's June place 0.85 with prompt comparison. Each figure depends on an
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
  full-likelihood check of the eccentricity test, and Phase 1b.
- Wording: the notes, the source manifest and the documentation describe the transcription and
  its checks by their method. The caveats stand: both passes and the first spot check share one
  method, so their errors may be correlated, and 251 rows of Table 77 have had no human check.
- Dataset v1.0.1 (concept DOI 10.5281/zenodo.24044791): rewords descriptive strings in
  `spotcheck_results.json` and the documentation; no table value or check result changes.
- Compute ledger: 1193 CPU-hours (Phase 1: 1090; Phase 1b: 103).
- Deferred from v1.2 to v1.3: the Peirce shape-model comparison (M6d), the held-out
  precoveries (M6e), the remaining M6c variants (the Bouvard and Adams data variants, the 1846
  reduction track, the Le Verrier shape model), and the resolution of the eccentricity test's
  low-distance cell (S1 cell 0).

## 1.1.0 (2026-11-10)

The distance claim of v1.0.0 held only under the Gaussian error model. v1.0.0 said that the
pre-1781 observations place the planet at 28.0–34.0 AU and prefer Neptune's distance over
Le Verrier's 36 AU by Δ ln L = 3.15. A second pre-registered profile likelihood in the
semi-major axis, M10, repeats M9 under the Student-t errors (ν = 4) that the evidence favours,
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
- Compute ledger: 1090 CPU-hours (M10: 3.1).
- Code DOI 10.5281/zenodo.24122940 (concept 10.5281/zenodo.24121059).
- Deferred from v1.1 to v1.2: the Peirce shape-model comparison (M6d), the held-out
  precoveries (M6e), the Bouvard and Adams data variants, the 1846 reduction track, the
  Le Verrier shape model, and the resolution of the eccentricity test's low-distance cell.

## 1.0.0 (2026-11-10)

First release: Phase 1 through M6. The pre-registered fit to Le Verrier's 1846 Table 77 places
the unseen planet at 325.8° on the night Galle found it (68 %: 323.6°–327.8°; Neptune 325.87°),
with the 95 % range provisional and the eccentricity test unresolved; the Gaussian profile
likelihood in the semi-major axis (M9); the site as a short paper with plates; the transcribed
Table 77 as a separate CC0 dataset (concept DOI 10.5281/zenodo.24044791). Code DOI
10.5281/zenodo.24121060 (concept 10.5281/zenodo.24121059).
