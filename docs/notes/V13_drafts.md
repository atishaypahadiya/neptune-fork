# v1.3: draft plans for the carried items

These drafts are carried to the v1.4 list.

**DRAFT — not frozen, not pushed.** Written 2026-11-12 on a local branch. Nothing here is a
pre-registration until it has been reviewed, moved into the note of its milestone, committed
and pushed, and no run below may start before that. No data that these runs would produce
was looked at while writing this. Choices marked **open** are not settled.

Why drafts and no launch: none of the four items has a plan that was committed and pushed
before 12 November 2026. `docs/release_checklist.md` ("Deferred to v1.1") says so for each:
M6d "needs a new pre-registered real-data nested run", M6e's positions "are not yet sourced
or transcribed, and the predictive step is not built", the M6c variants are "not yet
pre-registered", and S1 cell 0 needs a "new pre-registered plan". Each also needs new code.

Cost figures are projections from the ledger of earlier runs (`docs/notes/M6.md`, `M8.md`,
`P1B_adams_fit.md`): a free-a fit on Table 77 with the full likelihood took 28.5 CPU-h and
3 h 8 min on 12 workers; a held-distance fit (two cells) 6.8 to 7.3 CPU-h and 40 to 46 min;
one fit with the fast likelihood of M8, reweighted, 0.8 CPU-h (seed 1 of the headline) to
5.9 CPU-h.

| item | plan before tonight | new code | CPU-h, full likelihood | wall, 12 workers | CPU-h, fast likelihood |
|---|---|---|---|---|---|
| M6d shape-model comparison | none | a held-distance job on Table 77 | about 7 (14 with the Student-t twin) | about 45 min (1.5 h) | about 1 |
| M6c variants, modern track (`bouvard`, `adams`) | none | two data selections in the runner | about 55 | about 6.5 h | 2 to 12 |
| M6c variants, 1846 track (3 fits and the shape model) | none | the 1846 reduction track itself | about 92 | about 10 h | 3 to 19 |
| M6e precoveries | none | sourcing, a table, the predictive step | under 2 | under 15 min | not applicable |
| S1 cell 0 | none | a sub-cell partition of cell 0 | 60 to 110 | 6 to 12 h | 10 to 25, plus the full-likelihood check below |

## M6d: the Le Verrier shape model against free distance (draft)

- **Question** (spec, M6d): the probability each model puts on landing within 0.87° of
  Neptune on Galle's night, the ratio of the two, and the Bayes factor between them.
- **Fit to add:** the shape model on the headline's data (`leverrier`, modern track, the
  frozen error model), semi-major axis held at Le Verrier's 36.154 AU, e ≤ 0.3, the
  headline's priors on every other parameter, cells e < 0.12 and e 0.12 to 0.3 (the
  held-distance layout of Phase 1b step 2), nlive 1000, the headline's sampler settings,
  ln Z errors × 1.8. **Open:** the held value (36.154 AU from his elements of 31 August
  1846, or the Bode-rule 38.4 AU of his first paper); the seed namespace.
- **Compared with:** the published headline (free a). No rerun of it.
- **Reported:** P(|λ − 325.87°| ≤ 0.87°) under each, their ratio with its sampling error
  (bootstrap over the posterior samples), ln B(shape : free) with its error; the same
  three under the Student-t error model if S3's twin is run (**open**).
- **Reading rule to fix before launch:** bands for the ratio and for ln B, in the wording
  of step 2 of Phase 1b (|ln B| < 1 "no preference", and so on). **Open.**
- **Known limit, stated in advance:** the free-a posterior depends on the distance prior
  and on cell 2's repeatability (M6c), so the ratio inherits both.

**Code ready** (12 November 2026, commit fbef030): `fork.export.v13`:
`shape_space` (ln a held, two e cells, 10 dimensions), `peirce_probability`, `peirce_ratio`
(bootstrap interval), `log_bayes_factor`. Tests: `tests/test_v13.py`
(`test_shape_space_holds_the_distance`, `test_peirce_probability_and_ratio`). Smoke run on
synthetic data only, 300 iterations per cell, nothing saved. Measured cost: 5.3 ms per
likelihood call (529 measurements). Still to write when the plan is frozen: the job list,
seeds, cache and ledger rows, which depend on the open choices.

## M6c: the remaining data choices (draft)

- **Fits** (spec, M6c): `bouvard` (1781 to 1845 only) and `adams` (all but Flamsteed 1690)
  on the modern track; all three on the 1846 track; the shape model on `leverrier`, both
  tracks (the modern one is M6d's fit above).
- **Blocked part:** the 1846 reduction track does not exist in the code. It has to be
  specified (which constants and tables stand in for Le Verrier's reduction), built and
  tested on synthetic data before any plan that uses it. **Open:** whether v1.3 includes it.
- **Runnable part after review:** the two modern-track variants, with the headline's
  sampler, partition, priors and error model (σ_g and the group offsets as frozen for the
  headline; the ancient group's offset drops out of `bouvard`).
- **Likelihood: open.** Full (about 28 CPU-h a fit), or the fast likelihood with exact
  reweighting as validated in M8 (A8b to A8f), which failed only in S1's cell 0.
- **Reported, per fit:** the spec's list (1846 longitude: median and 68 % interval; the
  miss in degrees and in half-widths; mass and a; ln Z), and the M6d probability.
- **Reading rule to fix before launch:** what counts as "the data choice moves the
  prediction", in units of the headline's interval. **Open.**
- **Overlap to state:** Phase 1b step 2 already fitted Adams's own data set D (278
  measurements). `adams` here is a different set (Table 77 without 1690).

**Code ready** (commit fbef030): `fork.export.v13.select_rows` for `bouvard` (492 measurements)
and `adams` (527); `model.observation.structure_from_rows` gains `allow_empty_groups`, off
by default, so a choice can drop the ancient series. Tests: `test_data_choices_select_the_rows_of_the_spec`,
`test_empty_groups_are_refused_unless_asked_for`,
`test_synthetic_layouts_give_a_finite_likelihood`. Smoke runs on synthetic data laid out on
each choice's rows. Measured cost per likelihood call: 2.4 ms for `bouvard`, 4.6 ms for
`adams`, against 5.3 ms for all rows; the `bouvard` fit should then cost about half the
projection above. **Not ready:** the 1846 reduction track (no specification yet), and so the
four fits on it.

## M6e: held-out precoveries (draft)

- **Sourcing first, no fit touched.** Lalande, 8 and 10 May 1795 (*Histoire céleste
  française*; Mauvais 1847, *Comptes rendus* 24, 666), and Galileo, 28 December 1612 and
  28 January 1613 (Kowal and Drake 1980, *Nature* 287, 311; only its abstract has been
  read here). Record each position with its stated uncertainty, equinox and source page,
  transcribed twice, in a file that no fitting module imports (A6d's static check extended
  to name it).
- **Predictive step:** for each published posterior (headline, S1 to S3, and the M6c
  variants when they exist), 6,000 equal-weight samples integrated exactly to each date;
  the perturber's geocentric place, with 68 % and 95 % regions on the sky.
- **Reported, not pass or fail** (A6e): whether each observed place lies in the 95 %
  region; the angular miss of the posterior median. Read with Standish's caution on
  Lalande and with Galileo's as sketches.
- **Stated in advance:** the fits use 1690 to 1845; 1612 is an extrapolation of 78 years
  before the first datum, and the regions are expected to be wide.
- **Open:** the uncertainty to assign to Galileo's drawings; whether the test is on
  longitude only.

**Code ready** (commit fbef030): `fork.export.v13.predictive_places` (the perturber's
geocentric right ascension and declination of date per posterior sample, one integration per
sample, about 5 ms), `region_summary` (68 % and 95 % intervals, plane mean and covariance)
and `region_level` (the smallest credible region that contains an observed place). Tests:
`test_region_summary_and_level`, `test_predictive_places_on_prior_draws` (prior draws on
synthetic data, at an epoch already in the cache). **Not ready:** Earth's place at the 1795 epochs (not fetched).

**Positions: Lalande sourced, Galileo not** (12 November 2026;
`data/historical/precoveries/MANIFEST.md`). Lalande, 8 and 10 May 1795: Mauvais's reduced
apparent places from *Comptes rendus* 24, 668, with the observation times as corrected by
the volume's erratum, checked against his letter in *Astronomische Nachrichten* 26, col. 97;
two entries, 117 cells, 0 disagreements; equinox of date; no stated uncertainty. Galileo,
1612–13: not sourced. Kowal and Drake 1980 is not readable here beyond its abstract, and by
Standish's summary the record is a drawing that gives at best one dimension (Neptune on the
line from Jupiter through a star), not a place. The existing A6d check
(`test_no_reference_to_a_precovery_file`) already fails on any fit-path string that names a
precovery file, so it covers the new directory without a change. **Open, added:** whether
Galileo enters the test at all, and if so as a one-dimensional constraint.

## S1 cell 0 (draft)

- **Where it stands:** published single run, ln Z₀ = −366.84, P(a < 25 AU) = 0.001; four
  fast-likelihood seeds, −362.26 to −362.84, 0.109; the full-likelihood arbiter, two seeds
  at −361.03 and −362.84, disagreeing with each other beyond tolerance (`M8.md`). Outcome
  3: unresolved, 0.001 to about 0.2.
- **Proposal:** treat cell 0 (a 20 to 25 AU, e ≤ 0.5) the way the whole space was treated:
  cut it into sub-cells so that each holds one mode, sample each separately, and add the
  evidences. Sub-cells fixed before any sampling from what is already published (the
  arbiter's two medians, 313.3° and 315.0°, and the e ≈ 0.3 ridge): **open**, for example
  e in [0, 0.2), [0.2, 0.35), [0.35, 0.5] crossed with a in [20, 22.5), [22.5, 25].
- **Runs:** each sub-cell with the fast likelihood, 4 seeds, nlive 4000, exact
  reweighting; then one full-likelihood run of the sub-cell that carries the most
  evidence, as the check.
- **Test to fix before launch:** seeds agree within A8b's tolerances in every sub-cell
  (ln Z within 0.5, median within 0.35°), and the full-likelihood check agrees with the
  fast one in its sub-cell. If so, the summed ln Z₀ and P(a < 25 AU) replace the range; if
  not, the question stays unresolved and is reported so.
- **Alternative, open:** a sampler built for separated modes without a partition (for
  example parallel tempering), which would need its own validation on synthetic data
  first and is the larger job.
- **Stated in advance:** whichever way it falls, S1's published verdict is corrected only
  by an erratum line that keeps the old number visible.

**Code ready** (commit fbef030): `fork.export.v13.S1_SUBCELLS` (the example cut: six
sub-cells), `subcell_space`, `combine_subcells` (evidence of the union from the sub-cells'
evidences and prior masses, with its error). Tests: `test_subcells_partition_cell_0` (the
prior masses add up to cell 0's), `test_combine_subcells`. Smoke runs on synthetic data in
two sub-cells. Measured cost: 5.3 ms per full-likelihood call; the fast likelihood was not
exercised here. The cut itself is still open.
