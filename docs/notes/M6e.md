# M6e: do the fits predict Lalande's places of 1795? Pre-registration

Committed and pushed before any code is run at the 1795 epochs. Written 13 November 2026.
It replaces the draft in `docs/notes/V13_drafts.md` ("M6e: held-out precoveries"). Results
are appended below the line at the end; nothing above it is edited afterwards.

## 0 Scope and question

Lalande recorded a "star" on 8 and 10 May 1795 that Mauvais identified in 1847 as Neptune.
The fits of M6c use Uranus's observations of 1690 to 1845 only and never saw these two
places. **Question: does the posterior of the unseen planet, carried back to May 1795,
predict the places Lalande recorded?**

In scope: the two places of `data/historical/precoveries/lalande_1795.csv`. Out of scope:
Galileo's drawings of 1612 and 1613, which need their own pre-registration once Kowal and
Drake (1980) has been read.

This is A6e of the spec: reported, not pass or fail.

## 1 Inputs, fixed

**The observed places** (Mauvais 1847, *Comptes rendus* 24, 668, two entries, times as
corrected by the volume's erratum; SHA-256 of the table
`22484a07c4e8ee568ca3b84d4f35b16a9006aeb16d935a65f5ce0b095885849a`):

| | Paris mean time | apparent right ascension | apparent declination |
|---|---|---|---|
| 8 May 1795 | 11h 5m 11s | 212° 59′ 35″.0 | −11° 20′ 39″.1 |
| 10 May 1795 | 10h 57m 9s | 212° 56′ 36″.3 | −11° 19′ 38″.8 |

**Two posteriors, as published:**

- **The headline (Gaussian errors):** the reported posterior of M6c after the cell-2
  repeatability study, outcome `replace`: the headline run's cells 0, 1 and 3 with the
  merge of the two R2 runs as cell 2 (`fork.export.m6c_repeat.replacement_cells`; summary
  `site/data/m6c_replacement.json`, Galle-night median 325.815°).
- **Student-t errors, ν = 4 (S3, `s3_t4`):** its four cells (summary `fits.s3_t4` of
  `site/data/m6_fits.json`, Galle-night median 325.166°).

The cached nested results they are read from (`data/cache/nested/`, not in the repository):

| file | SHA-256 |
|---|---|
| `m6c/headline__c0.pkl` | `529407e435cadf99aa86889934e4fa5b525f2bc1b502973a8150aa720a665787` |
| `m6c/headline__c1.pkl` | `1cc6ff7b42bddcd4da80437aeb24d26ba511ec53dc975eabfe897c7a37812788` |
| `m6c/headline__c3.pkl` | `3fe3ee5434808d276d4db5aded4c4afc9ee6189bfb42af62cf223953d4f7e6b7` |
| `m6c_repeat/r2_1.pkl` | `72401426b028bbe77c4504d2f815c61705612d5a50896111eff5ac3508d34fe0` |
| `m6c_repeat/r2_2.pkl` | `33df8e2d5bafd97c3b68ee4c39fdb9b1e8470400a8f59d01454d62ee0e871b4f` |
| `m6c/s3_t4__c0.pkl` | `6b7f2c50ca356b32b615240f160d2eab91edd85baa482c10bd96bd5b6774d4d0` |
| `m6c/s3_t4__c1.pkl` | `8f6203df183399be1ac9f39f43cb6630b5f7e7534671da9417b9a48225100b4f` |
| `m6c/s3_t4__c2.pkl` | `d813abb56ff29af396248f099ccdcffea497e879d346684d929c41275e5c60f0` |
| `m6c/s3_t4__c3.pkl` | `ef7e0e6efa479235baf047de523c3f9a2910ba12816e1d3380f7a9a12b7c8ecb` |

The runner refuses to run if a hash differs. The published summaries:
`site/data/m6c_replacement.json` `d4d1b5c862169984b6152c38b57809da36404e8f84043c5fdf658c05b3f42e32`,
`site/data/m6_fits.json` `2d111e7fcb7a71100a72a8f963f6c8f0f26aea9da90ff22a246893618bb74fff`.
As a check that the right posterior is loaded, the Galle-night median of each sample set
drawn below must agree with the published median within 0.2° (a different random
subsample; the published 68 % widths are 4.2° and 6.3°).

**Samples.** Per model, the equal-weight samples of the partitioned posterior, then 6,000
of them without replacement; generator seeds `[20261113, 0]` (headline) and
`[20261113, 1]` (Student-t).

## 2 The predicted quantity

For each sample: the perturber's state integrated to the instant; its geocentric direction
with light time; the **apparent place referred to the true equator and equinox of date**
(precession, nutation and annual aberration applied), which is what Mauvais printed. From
it, right ascension, declination and the ecliptic longitude of date (true obliquity). The
observed places are converted to ecliptic longitude by the same rotation.

| step | handling | test, written before the run |
|---|---|---|
| time | Paris mean time, astronomical reckoning from mean noon: 11h 5m 11s on 8 May is 23h 5m 11s civil. UT = Paris mean time − 9m 20.9s (the Observatory's longitude, 2° 20′ 14″ E). TT = UT + 15 s; TDB taken equal to TT. Instants: 1795-05-08 22:56:05 and 1795-05-10 22:48:03 TDB | the observations are meridian transits: the apparent sidereal time at Paris at each instant must equal the observed right ascension within 1 minute of time. If it does not, the work stops and is reported |
| Earth | its place at the two instants about the inner barycenter, as in the fit, from JPL Horizons through the project's cache (one fetch, added to the CI fixture) | the existing ephemeris tests |
| light time, precession | the existing code (`geocentric_lon_lat_ecl_j2000_light_time`, the precession of `fork.data.frames`) | existing tests |
| nutation, aberration | new, with the IAU routines of `erfa` | the whole chain, applied to Uranus's ephemeris vector at the two instants, must reproduce JPL's geocentric apparent place of Uranus within 1″ in each coordinate. Uranus only; Neptune's ephemeris is not read |
| observer's parallax | not applied (under 0″.3 for Neptune) | |

An error of 5 s in ΔT moves Neptune by 0″.005. Mauvais reduced with the constants of his
day; a place of date depends on them by under 1″, which is not corrected.

No observational error is folded into the prediction. Mauvais states none; Lalande's zone
places are good to a few arcseconds, and the predictive region is expected to be degrees
wide. If the 68 % width in longitude comes out under 2′, that is said in the results and
the reading is marked as limited by the unknown error of the observation.

## 3 Reading, fixed in advance

**Primary: the central credible level of Lalande's place in ecliptic longitude**, per model
and per night. With F the fraction of the 6,000 draws whose longitude is below the observed
one (differences taken on the circle about the draws' circular median), the level is
L = |2F − 1|.

| L | class |
|---|---|
| ≤ 0.68 | inside 68 % |
| over 0.68 and ≤ 0.95 | inside 95 % |
| over 0.95 | outside |

A level within 0.01 of 0.68 or 0.95 is reported as "at the boundary" with the class on each
side named (the Monte Carlo error of L at 6,000 draws is under 0.006).

What each outcome means:

| outcome | meaning |
|---|---|
| inside 68 % | the posterior, which never saw 1795, predicts Lalande's longitude as a typical draw. Its strength is read from the width beside it (§ 4): a wide region makes this weak |
| inside 95 % | consistent, in the tail |
| outside | the orbit the fit infers does not carry back 51 years, although its longitude in 1846 is right. This is possible and is said here in advance: the fit constrains the planet's distance and eccentricity poorly (headline a = 30.2 AU, 95 % 28.7 to 32.8), and the place in 1795 depends on both. It is not a failure of the 1846 prediction |

- The two nights are two days apart (Neptune moved 3′). They are one test seen twice, not
  two. If their classes differ, both are given and the reading is "at a boundary between
  nights".
- The two models are read separately. If they differ, both are stated; no combined verdict.

## 4 Reported, not read

- **The tangent-plane level:** the fraction of draws whose Mahalanobis distance from the
  sample mean on the tangent plane (right ascension and declination) is no larger than the
  observed place's (`fork.export.v13.region_level`).
- The median and the 68 % and 95 % intervals in right ascension, declination and ecliptic
  longitude; the miss of the median in each.
- The predicted motion between the two nights against the observed (−2′ 58″.7, +1′ 0″.3).
- **The prior predictive:** 6,000 draws from the fit's prior (a cell chosen with its prior
  mass, then that cell's prior transform of a uniform draw; seeds `[20261113, 10]` and
  `[20261113, 11]`), carried to the same instants, with the same levels L and tangent-plane
  levels of Lalande's place.
- **The ratio of widths:** the posterior's 68 % and 95 % widths in ecliptic longitude
  divided by the prior predictive's, per model and night. The prior of the mean longitude
  is uniform, so the prior widths are expected near 245° and 342°.

## 5 What this cannot show

- Lalande's positional error. The comparison is at the scale of the predictive region.
- That the star was Neptune. The identification was made in 1847 from the modern orbit and
  is taken as given.
- An independent confirmation of the 1846 result: orbits that pass 1846 correctly are
  correlated at 1795.
- A blind analysis. Neptune's place in 1795 is public knowledge and the observed places
  were read while they were sourced. Only the fits are blind. Everything in §§ 1 to 4 is
  fixed before any prediction at these instants is computed.

## 6 Leak check

- The fits' inputs are Table 77 (Uranus only) and the error model, frozen at their
  pre-registration of 7 November 2026 (`fork.export.m6c.FROZEN`); the runner checks the
  hashes.
- The fits finished on 7 November 2026 (headline) and 8 November (S3 and the R2 runs;
  `docs/notes/M6.md`). The Lalande table was sourced on 12 November 2026 (`data/historical/precoveries/MANIFEST.md`).
- Existing tests, kept: no module on the fit path names a precovery file
  (`test_no_reference_to_a_precovery_file`) or reads Neptune's ephemeris (`tests/test_a6d.py`).
- New tests: `fork.export.m6e` is the only module under `src/` that names the precovery
  directory, and it is not reachable from the fit's entry modules.
- Disclosure: the predictive code (`fork.export.v13.predictive_places`) has been run on
  prior draws of synthetic data at an epoch of 1846 only. No posterior sample has been
  carried to 1795.

## 7 Order, cost, layout

Order: this plan committed and pushed; the code and its tests, run on prior draws or
synthetic data only; the Earth fetch and the JPL apparent place of Uranus; the time and
frame tests; the run; results below the line; a report before the results are pushed.

Cost: 24,000 integrations (two posteriors and two prior sets of 6,000, each to the two
instants) at about 5 ms, about 2 CPU-minutes, plus loading the cached results. Not
detached. One ledger row.

Layout: `src/fork/export/m6e.py`; `tests/test_m6e.py`; `docs/notes/M6e.{json,png}`. The
site is not touched.

---
