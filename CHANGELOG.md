# Changelog

Versions follow semver. Dates are the release dates; the DOIs are in `CITATION.cff`.

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
