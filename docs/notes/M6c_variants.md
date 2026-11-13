# M6c: the two remaining data choices on the modern track (pre-registration)

Written 2026-11-13, before any code of this plan computed a quantity its reading uses.
Draft: `docs/notes/V13_drafts.md` ("M6c: the remaining data choices"). Everything above
the line at the end is frozen by the commit that adds this file; amendments and results go
below it, dated.

## 1 Question

The spec's M6c asks for the free-distance fit under three data choices. `leverrier` (all of
Table 77) is the published headline. This plan runs the other two on the modern reduction
track:

- `bouvard`: the modern observations only, 1781 to 1845 (492 measurements);
- `adams`: all but Flamsteed 1690 (527 measurements).

Question: does the choice of data move the predicted longitude on Galle's night, and by how
much, measured against the headline's interval?

**Not in this plan:** the 1846 reduction track (three fits and the shape model on it). It
has no specification: which constants and tables stand in for Le Verrier's reduction is not
written down anywhere in the spec or the notes. I do not invent one. Those four fits stay
blocked until a specification exists.

## 2 Inputs

| input | SHA-256 |
|---|---|
| `data/historical/leverrier_1846.csv` | 2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4 |
| `data/historical/m6_error_model.json` | 0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277 |
| `docs/notes/M8_u_ref.json` (the headline's reference orbit of Uranus) | 4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5 |
| `src/fork/model/perturbative.py` | 092f673b9352f9d2b0ae709f2a4d3f258dd09806245b767fd83b86e8ecccd68c |
| `src/fork/model/linearized.py` | 3e02e78bdd9ee74fb96a6666f192188c355113addc92f43e376729f0b616d15e |
| `src/fork/inference/reweight.py` | d0c5fbff6066ebaa531bc24668fa807a8f7bf0aa96cf77eead5c07aca9cd021b |
| `site/data/m8b.json` (the comparison fit) | 355104b24966ce216f53ff603931a060a4a5e738a996a683db5fac5c5e962027 |

- Row selections: `fork.export.v13.select_rows` (commit fbef030), 492 and 527 measurements.
- Comparison fit: the headline at the higher-resolution settings (`m8b.json`,
  `headline/merged`): median 325.735°, 68 % 323.511° to 327.672° (width 4.161°), 95 %
  319.593° to 330.213° (provisional), ln Z −361.523 ± 0.078.
- No Neptune ephemeris and no precovery file enters a fit. The comparison target 325.873°
  is used in reporting only.

## 3 Method

Each fit is made the way the published headline was confirmed (M8b), with nothing changed
but the rows:

- **Model and priors:** the headline's. Free distance, e ≤ 0.3, inclined perturber, the
  same priors on every parameter. Uranus's (i, Ω) held, its in-plane four free in the box
  given by the M5 rule on the chosen rows.
- **Error model:** σ_g and τ_g as frozen for the headline. In `bouvard` the ancient groups
  have no rows, so their offsets drop out.
- **Partition:** the headline's four cells: a 20 to 25 AU; 25 to 28 AU; 28 to 60 AU with
  e < 0.12; 28 to 60 AU with e 0.12 to 0.3.
- **Likelihood:** the fast (perturbative) likelihood, second order off, 60 d nodes, cubic
  quadrature, about the headline's reference orbit; then every draw reweighted exactly to
  the linearised marginal with refinement (mode (b)), as in M8.
- **Sampler:** dynesty static, multi-ellipsoid bounds, rslice with 28 slices, dlogz 0.1;
  nlive 8000 in cell 0 and 4000 in the others; four seeds per cell; one process per run.
  rng `[2026, 15, 1846, f, k + 1, s]`, f = 1 for `bouvard`, 2 for `adams`, k the cell,
  s = 1..4. That is 32 runs.
- **Reweighting:** 6,000 equal-weight draws per run and per four-seed merge (rng
  `[..., s, 7]`, s = 0 for the merge; chunks of 250 with rng `[..., s, 9, c]`). Each draw's
  Galle-night longitude through the full model, with Uranus drawn from mode (b)'s Gaussian.
- **The fit:** the four merged cells combined with their prior masses
  (`reweight.combine`), 6,000 draws per cell; and the four single-seed fits likewise.
  ln Z errors: the sampler's × 1.8, plus the reweighting's Monte Carlo error.

## 4 Gates and tests, in this order

- **G0, on synthetic data, before any run on the real rows.** One synthetic data set per
  choice on its own rows (`v13.synthetic_dataset`, seed 1; no observed value is read). The
  whole chain at a quarter of the live points, two seeds. Pass: every cell posterior has
  ESS/N ≥ 0.3 and a truncation bound below 10⁻⁶ in −2 ln L, and G4 passes on 50 draws. If
  G0 fails for a choice, that choice is not run on the real rows and is recorded as stopped.
- **G1:** every run stops on dlogz.
- **G2:** ESS/N ≥ 0.3 for each of the 40 cell posteriors (32 runs and 8 merges); ESS/N ≥
  0.5 of the between-cell-corrected weights for each fit (the fit-level gate of M8).
- **G3:** the truncation bound of mode (b) is below 10⁻⁶ in −2 ln L at every draw.
- **G4, the target itself.** Mode (b) was checked against the full likelihood on the
  headline's rows only. On 200 draws from each merged fit: the profile χ² of mode (b)
  within 0.1 of the full likelihood minimised over Uranus, and its marginal −2 ln L within
  0.2 of the Laplace value with the exact Hessian (A8a's two tolerances). Pass: both hold
  at 95 % of the draws at least, and no draw is beyond five times a tolerance.
- **T1, seeds.** In each cell: the four seeds' 68 % widths within 10 % of their mean (A8f)
  and their ln Z within a range of 0.5. Across the four single-seed fits: the range of the
  median ≤ 0.35°. Cells with less than 1 % of the fit's evidence are reported and not read
  in the width test.

## 5 What is reported and how it is read

Per fit, the spec's list: the longitude on Galle's night (median, 68 % and 95 %
intervals), the miss from 325.873° in degrees and in half-widths, mass and distance, ln Z
with its error, the cells' shares, and P(|λ − 325.873°| ≤ 0.87°).

Reading, fixed now. It is made only if G1 to G4 and T1 pass for the fit; otherwise the
fit's numbers are reported as not repeatable or not validated, and nothing is read.

| quantity | value | wording |
|---|---|---|
| shift of the median from the headline's 325.735° | below 0.35° | the choice does not move the prediction |
| | 0.35° or more, median inside the headline's 68 % interval | it moves the prediction within the headline's 68 % interval |
| | median outside the 68 % and inside the 95 % interval | it moves the prediction outside the headline's 68 % interval |
| | median outside the 95 % interval | it moves the prediction outside the headline's 95 % interval |
| 68 % width over the headline's 4.161° | 0.9 to 1.1 | the same width |
| | above 1.1, below 0.9 | wider, narrower, with the factor |
| 325.873° in the fit's 68 % and 95 % intervals | in or out | stated as such |

- 0.35° is the tolerance the headline's own repeated runs were held to; a shift below it
  cannot be told from sampler noise.
- A value within 0.02° of a threshold, or a ratio within 0.005, is reported as "at the
  threshold".
- No Bayes factor between data choices is formed: the data differ.
- The 95 % intervals of these fits are reported with the same "provisional" label as the
  headline's unless the std of the four single-seed fits' 95 % bounds is 0.4° or less.

## 6 What this cannot show

- It cannot say which choice is right. It measures how far the answer depends on the choice.
- `adams` here is Table 77 without 1690. It is not Adams's own data: that set was fitted
  separately (the Adams–Challis study, 278 measurements) and is not comparable row by row.
- `bouvard` here is the modern rows of Table 77 as Le Verrier reduced them, not Bouvard's
  1821 tables.
- Both fits use the modern reduction track and the headline's error model, which was
  calibrated on all the rows.
- Nothing here bears on the 1846 reduction track.

## 7 Choices in the draft, and how each is resolved

| open in the draft | resolution |
|---|---|
| full likelihood or fast likelihood with exact reweighting | the fast likelihood with exact reweighting, at the settings and with the gates of the run that confirmed the headline, four seeds per cell. A single full-likelihood run per cell was shown not to be repeatable in the headline's main cell, and four seeds of it would cost about 115 CPU-h a fit. G4 covers what the earlier validation did not: mode (b) on other rows |
| the reference orbit of the fast likelihood | the headline's, since no posterior exists for these rows. A poorer reference lowers ESS and cannot bias a reweighted result; G0 and G2 catch it |
| the reading rule | § 5, in units of the headline's interval |
| whether the 1846 track is included | no: blocked, no specification |
| the M6d probability per fit | reported, not read here |

## 8 Cost

From the ledger: the headline at these settings took 5.6 CPU-h of sampling and 0.3 CPU-h
of reweighting for 16 runs. Two fits: about 12 CPU-h, plus G0 (about 2) and G4 (under 1).
About 1.5 h on 12 workers. Runs start after the S1 cell 0 runs free the workers.

## 9 Layout

- Code: `src/fork/export/v14_fits.py`, tests `tests/test_v14_fits.py`, branch `v14-runs`.
- Cache: `data/cache/nested/v14_fits/<fit>/`. Results: `docs/notes/M6c_variants.json`.
- Ledger rows: `v14/<fit>/c<k>/s<s>`, `v14/<fit>/reweight`.

---
