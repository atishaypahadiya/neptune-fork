# Phase 2 milestones

**DRAFT, scoping only, not a pre-registration.** 2026-11-12. Nothing here is frozen. No
threshold in this file has been committed anywhere. Numbers are estimates to be replaced.

Case: Mercury's perihelion. Question: given timings or positions only, does the engine prefer an
unseen mass, a change in the law, or neither, and how often would it be wrong?

## Cost basis

From `docs/reproduce.md` (ledger figures) and `docs/notes/M8.md`:

| Phase 1 figure | value |
|---|---|
| all of Phase 1 | 1090 CPU-h |
| headline fit, full likelihood, 4 cells | 28.5 CPU-h, 17.3 million calls |
| full likelihood per call | 5.52 ms |
| fast likelihood, 16 runs | 2.1 CPU-h, 27.7 million calls (0.13 CPU-h per run) |
| reweighting against the full likelihood | 4.8 CPU-h |
| synthetic coverage, 40 jobs (A6c) | 277 CPU-h |
| repeatability runs R1 and R2 | 103.6 CPU-h |
| full-likelihood arbiter, 2 runs | 98 CPU-h |
| M5b with validation | about 435 CPU-h |

Phase 2 unit costs assumed (`design_sketch.md`, section 6): full likelihood 0.13 s per call,
about 24 times Phase 1; fast likelihood 0.5 ms per call; 4 million calls per nested run. So one
full-likelihood nested run is about 145 CPU-h and one fast run about 0.6 CPU-h. Reweighting
20,000 posterior samples against the full likelihood is 0.7 CPU-h per run. Every figure below
inherits a factor of 2 either way from the per-call estimate.

## Milestones

### P2-M0: Rules, generic interfaces, pre-registration template
- Deliverable: Phase 2 additions to `docs/RULES.md` (Mercury's six elements always free; no JPL
  Mercury state held; no observed-minus-computed column read before the freeze; inputs allowed
  from modern knowledge listed one by one). A model-agnostic problem object replacing the
  Uranus-specific `PosteriorSpace` and `LongitudeObservations`. A pre-registration template with
  the decision rule left blank.
- Acceptance idea: every Phase 1 test passes unchanged through the generic interface; the Phase 1
  headline likelihood is reproduced to 1e-9 in ln L.
- CPU: below 1 CPU-h. Basis: Phase 1 test suites run in minutes.

### P2-M1: Inner-system integrator
- Deliverable: eight-planet Wisdom-Holman configuration with Mercury; step chosen by a
  convergence study; inner planets no longer folded into the Sun.
- Acceptance idea (A1 analogue): Mercury's heliocentric longitude agrees with REBOUND at the same
  step to 0.001 arcsec over 1677-1881; against a high-accuracy reference the Newtonian perihelion
  rate is converged to 0.01 arcsec per century; the rate itself is compared with the published
  values (Park et al.'s total of 575.31 less 42.98 is 532.3; Newcomb 1882, p. 473, reads 531.83
  in the OCR; Narlikar and Rana give 528.95), with the differences explained in the note;
  convergence order 2.
- CPU: 2 CPU-h. Basis: a few hundred 204-year integrations at 0.1 to 1 s each.

### P2-M2: Force-law family
- Deliverable: position-only kick terms: effective 1PN potential with factor xi, Hall's
  exponent, circular ring, solar J2. `model/modified_law.py`; ring and J2 in
  `model/hidden_mass.py`.
- Acceptance idea: two-body perihelion and node rates match the closed forms to 1e-4 relative
  (42.98 for xi = 1 at Mercury; 8.62, 3.84, 1.35 for Venus, Earth, Mars; (delta/2) n for Hall;
  the Laplace-coefficient rate for a ring at small eccentricity). The effective potential agrees
  with a full 1PN reference (REBOUNDx, tests only) to 0.01 arcsec over 204 years. Energy error
  stays bounded with each term switched on.
- CPU: 2 CPU-h. Basis: as P2-M1.

### P2-M3: Transit-contact observation model
- Deliverable: predicted time of geocentric contact with two-body light time, semidiameters,
  Earth-Moon offset, Delta T to 1881; the contact-time likelihood with weights, offsets and a
  Student-t variant; optional apparent ecliptic places for meridian data.
- Acceptance idea: with JPL states, the four contacts of each of four recent transits agree with
  JPL Horizons to 1 s. Internal: a 1 arcsec shift of Mercury along the track moves a November
  contact by about 10 s and a May contact by about 16 s.
- CPU: 1 CPU-h.

### P2-M4: Transcription of the first data set
- Deliverable: Newcomb 1882, pp. 451-452, observed columns only (contact, date, observed time,
  weight), by the Phase 1 method: two passes, reconciliation against a second scan, spot check
  of every row (there are only 55 to 60). Provenance file. The tabular times, coefficients and
  residuals are transcribed, if at all, into a separate file that the pipeline cannot read
  before the freeze.
- Acceptance idea (A6a analogue): the two passes agree after reconciliation; interior ingress
  and egress times of each transit differ by a plausible duration; dates match the known list of
  transits.
- Effort: 0.1 to 0.2 Phase 1 units. CPU: none.
- Optional second stage: section 8 of the memoir (1.5 to 5 units), only after P2-M6 shows it
  would change the answer.

### P2-M5: Synthetic recovery and the identifiability map
- Deliverable: synthetic contact times at the real epochs and weights from (i) JPL Mercury,
  which contains general relativity, (ii) a Newtonian integration with a planted ring, (iii) a
  Newtonian integration with nothing planted. Point fits and profiles for each model on each
  truth. A map of which ring planes and radii are distinguishable from the law models at the
  real noise level. The fast (perturbative) likelihood for all models, validated.
- Acceptance idea (A3 analogue): planted xi recovered within its 68 per cent interval in 5 to 8
  of 10 seeds; planted ring mass likewise; with nothing planted, no model beats the null. Fast
  against full likelihood: ln L differences below 0.1 over the posterior bulk, as in M8.
- CPU: about 320 CPU-h. Basis: 20 fast nested runs (12); reweighting them (14); two
  full-likelihood reference runs at 145 each (290). In Phase 1 the corresponding work (M5b with
  validation) took 435 CPU-h. If the reweighting diagnostics pass, the two reference runs can be
  dropped and the cost falls to about 30.

### P2-M6: Calibration of the decision rule
- Deliverable: the fork's error rates. For each of the three synthetic truths, 10 seeds, three
  models: evidences, Bayes factors, and the verdict under a candidate rule. Sensitivity of the
  verdict to the priors on the ring's plane and radius and on xi or delta. A written decision
  rule with an explicit "undecided" outcome.
- Acceptance idea: under the candidate rule the false-decision rate on synthetic data is below
  a stated level, and "undecided" is returned where `design_sketch.md` section 5 predicts it.
  If the rule returns "undecided" on every realistic synthetic case, report that and stop: the
  historical run is then a demonstration, not a test.
- CPU: about 125 CPU-h. Basis: 90 fast runs (54) plus reweighting (63), rounded. Phase 1's
  synthetic coverage test took 277 CPU-h for 40 full-likelihood jobs.

### P2-M7: Pre-registration and the blind historical run
- Deliverable: a committed plan (models, priors, error model, decision rule, sensitivity runs,
  stopping rules) before the pipeline reads the transcribed times. Then the run and its note.
- Acceptance idea (A6 analogue): thresholds written before the run; the perihelion excess
  implied by the xi posterior is reported against Newcomb's 43 only after the verdict is frozen.
- CPU: about 650 CPU-h in the cautious version. Basis: 12 fast runs, three models by two error
  models by two data variants (8), their reweighting (9); full-likelihood confirmation of the
  three headline models (435); a reserve of 200 for repeatability and an arbiter, which Phase 1
  needed (103.6 and 98). In the lean version, reweighting replaces the confirmation runs and the
  cost is about 50, with the reserve kept: 250.

### P2-M8: Multi-planet digest (optional)
- Deliverable: the same three models against Newcomb's 1895 table of secular discrepancies
  (pp. 109-110), as a linear-Gaussian problem with analytic evidence.
- Acceptance idea: none blind is possible; the scoping note already looked at this table. Report
  as descriptive.
- CPU: below 5 CPU-h.

### P2-M9: The Uranus case through the same fork (optional)
- Deliverable: a modified-law branch (Hall's exponent and an inverse-cube term) fitted to
  Table 77 with the Phase 1 observation model, and its evidence against the existing hidden-mass
  evidence. One engine, two cases, opposite expected answers.
- Acceptance idea: pre-registered like M6c; expected outcome stated in advance (hidden mass).
- CPU: about 25 CPU-h. Basis: one-parameter models at the Phase 1 full-likelihood cost, four
  runs at about 6 CPU-h (4 million calls at 5.5 ms).

### P2-M10: Export, notes and site
- Deliverable: `site/data` exports, ledger rows, milestone notes, a results page.
- Acceptance idea: number audit as in Phase 1.
- CPU: 5 CPU-h.

## Totals

| version | CPU-h |
|---|---|
| cautious: full-likelihood references and confirmations (P2-M5 320, P2-M7 650) | about 1135 |
| lean: reweighting in place of full-likelihood nested runs, reserve kept | about 445 |
| range allowing for the factor of 2 in the per-call cost | 250 to 2200 |

The cautious total is close to Phase 1's 1090 CPU-h. Three quarters of it is full-likelihood
nested sampling at about 145 CPU-h per run. Whether those runs are needed is decided by P2-M5.

## Order and stop conditions

1. P2-M0 to P2-M3 can proceed without data.
2. P2-M4 does not touch the pipeline.
3. P2-M5 and P2-M6 use synthetic data only (the Phase 1 rule "synthetic before historical").
4. Stop after P2-M6 if the decision rule cannot separate the branches on synthetic data at the
   real noise level. `design_sketch.md` section 5 expects that outcome for a ring in Mercury's
   plane. The stop is itself the result to report.
5. P2-M7 needs explicit approval of the plan, as every Phase 1 run did.
