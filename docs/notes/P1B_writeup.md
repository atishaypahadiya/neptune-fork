# Phase 1b write-up: Adams's solutions and Challis's search

Date: 2026-11-12. This note is the longer form of the companion page (`site/adams-challis/`,
specified in `docs/DESIGN.md`, "Companion page"). It states the two questions, the answers the
steps gave, and where each number comes from. It adds no result: every number below is a field
of a committed JSON in `docs/notes/`, named in the claims table, and the page reads the same
fields through `site/data/p1b_adams.json` and `site/data/p1b_challis.json`
(`src/fork/export/p1b_site.py`; `tests/test_p1b_site.py` checks each copy against its source).

The plans, runs, results and deviations of each step are in their own notes and are not
repeated: `P1B_sources.md` (step 0), `P1B_adams_solutions.md` (1), `P1B_challis_sweeps.md` (2),
`P1B_adams_fit.md` (3), `P1B_adams_reproduction.md` (3b), `P1B_adams_noise.md` (3c),
`P1B_challis_detection.md` (4).

## Q1. How good was Adams?

**Finding 1 (step 3). At his distances the project's fit is 9 to 10° short of Neptune.**
Adams assumed the planet's mean distance: 38.4 AU under his first hypothesis, 37.25 AU under
his second. The project's model, fitted to the 145 observations of Table 77 comparable to his
data with the distance held at those values, puts the true heliocentric longitude at
1846.762 (5 October 1846, the epoch at which every solution is compared) 10.5° short of Neptune under the first (68 % range −14.8° to −6.9°) and 9.1° short under
the second (−12.6° to −4.9°). Neptune is outside the 95 % range of the geocentric longitude in
both. With the distance free the same data put Neptune inside the 68 % range. His printed
solutions are +1.7° and +2.5° from Neptune.

**Finding 2 (step 3b). His solutions are reproduced from the memoir.** From his printed normal
equations the planet's longitude at epoch agrees with his to 0.011° and 0.003°. From his
tabulated data the same method gives a mean longitude 1.44° and 1.76° from his print and
within 0.02° and 0.01° of Brookes (1970). The difference sits in one normal equation, whose
coefficients he rounded to four decimals: one unit of the fourth decimal of one coefficient
moves the longitude by 2.9° (first hypothesis) and 1.3° (second).

**Finding 3 (step 3b). The shortfall belongs to the assumed distance.** On one noise-free
Neptune, sampled at his epochs and solved at 38.4 AU, every method is short: his least squares
by 9.5°, his own elimination by 9.9°, the project's model on his normal places by 7.8°, the
step 3 fit on synthetic data by 10.4° (the step 3 fit on the real data: 10.5°). His real data
through his least squares land +1.6° from Neptune. His nearness is therefore in his data or
their reduction, not in his method or his distance.

**Finding 4 (step 3c). Whether noise alone accounts for the nearness is not settled.** With
independent noise added to the noise-free places, 10,000 draws per noise model, the
pre-registered statistic (the share of draws that move toward Neptune by at least his 11.2°)
is, by hypothesis and noise model, H I: 0.183 (a), 0.074 (b); H II: 0.182 (a), 0.084 (b),
which the pre-registered rule reads "noise suffices". Model (a) is the scatter of his own
printed residuals; model (b) is the scatter of the project's data carried to his normal places. Deviation: the statistic treated the shift as a number on a line. The least-squares
problem has a second minimum (eccentricity 0.45, 165° from Neptune, 6.9 in χ² above the best),
and 92 % of the draws counted as exceeding his shift are jumps to it. Post hoc, among draws
whose solution stays within 30° of the real-data solution, the probability is H I: 0.017 (a), 0.011 (b);
H II: 0.046 (a), 0.037 (b): "weak" by the same thresholds, below 0.05 and
above 0.003. Neither luck nor a systematic cause is established. Bouvard's tables, on which
his reductions rested, are the route to an attribution; they were not transcribed.

## Q2. Could Challis's search have recognised the planet before 23 September?

**The model (step 4).** Each session is three events: the planet's place lies in the zone, it
is recorded, it is recognised (by a noted disc, by comparing two series, or against a chart).
What the record does not print is drawn from declared priors; a draw is kept only if his
actual procedure then reproduces the record. 100,000 of 1,003,853 draws are kept (10 %).

**The reading.** Primary, by ruling before the run: the realised events are held. A recording
that a case shares with the history keeps what happened (4 and 12 August: recorded, no disc
noted); only events the history does not contain are drawn.

| case | probability before 23 September (5 to 95 % range) | reading | note |
|---|---|---|---|
| his actual procedure | 0.00 | very probably not | by construction: the realised events are held fixed |
| the 12 August comparison with 30 July, carried through | 0.71 (0.50 to 0.97) | depends on p_follow | complete by 16 August in 99 % of the draws that flag the planet; the record does not inform p_follow (prior mean 0.75, posterior 0.75) |
| the 4 August sweep compared with 12 August | 0.28 (0.00 to 0.93) | depends on p_follow | |
| the Berlin sheet he had, equinox of 1800 | 0.00 | very probably not | the planet is off the sheet |
| the same sheet, right ascension of date | 0.92 (0.77 to 1.00) | depends on p_follow | prior-dependent |
| zones from Adams's first place (S4) | 0.41 (0.04 to 0.81) | depends on p_disc_routine | prior-dependent |
| from S4, with prompt comparison | 0.79 (0.22 to 0.99) | depends on the number of later sessions | |
| zones from Adams's second place (S5), with or without prompt comparison | 0.00 | very probably not | by timing |
| zones from Le Verrier's June place | 0.41 (0.04 to 0.81) | depends on p_disc_routine | prior-dependent |
| from Le Verrier's June place, with prompt comparison | 0.85 (0.63 to 0.99) | depends on p_follow | |
| the 29 September procedure, begun on 29 July | 0.35 (0.02 to 0.90) | depends on the breadth of the sweep | |

The chart case as a whole depends on the Berlin sheet's limits, taken from a secondary source.
No case reaches "very probably".

**Deviations that bear on these claims.**

1. Under the held reading his actual procedure is 0 by construction (the realised events are held
   fixed, so the cases measure what a change would have added), so the pre-registered
   paired-difference rule against it ("improves when the range excludes 0") is uninformative.
   Identified after reading the results. Its statements are not quoted on the page.
2. Three readings change row when every prior interval is halved and are marked
   prior-dependent. The priors on the zone's setting were not halved. That check used 20,000
   kept draws; the ranking of inputs used the first 10,000.
3. Declared choices: the Northumberland's aperture from Airy (1844); the start places at
   ecliptic latitude 0; zones alternating outward from the start place; Neptune's magnitude
   extrapolated from Mallama and Hilton (2018).

**Outside the pre-registered reading: the replay.** If every event is drawn afresh, his
procedure taken as a plan, the probabilities are (mean, 5 to 95 % range):

| case | replay |
|---|---|
| his actual procedure | 0.34 (0.03 to 0.70) |
| the 12 August comparison, carried through | 0.81 (0.56 to 0.98) |
| the 4 August sweep compared with 12 August | 0.52 (0.05 to 0.96) |
| the Berlin sheet, equinox of 1800 | 0.34 (0.03 to 0.70) |
| the Berlin sheet, right ascension of date | 0.94 (0.83 to 1.00) |

The start-place cases and the 29 September procedure share no recording with the history and
are the same under both readings. The replay answers a different question (how good the plan
was before it was carried out) and is reported, not read against the rule.

**Post hoc, description only: pairs of counterfactuals on the shared draws.**

| pair | mean difference | 5 to 95 % range | share of draws above 0, below 0, equal |
|---|---|---|---|
| 12 August comparison minus 4 August comparison | +0.43 | 0.00 to +0.96 | 0.59, 0.02, 0.39 |
| 12 August comparison minus the 29 September procedure | +0.37 | −0.27 to +0.85 | 0.82, 0.18, 0.00 |
| S4 with prompt comparison minus 12 August comparison | +0.08 | −0.50 to +0.37 | 0.86, 0.14, 0.00 |

The ties in the first pair are the draws where both comparisons cover the planet or neither
does (the dated correction in `P1B_challis_detection.md`).

## Claims, sources, deviations

`fit` = `P1B_adams_fit.json`, `repro` = `P1B_adams_reproduction.json`, `noise` =
`P1B_adams_noise.json`, `det` = `P1B_challis_detection.json`, all in `docs/notes/`.

| claim | JSON path | deviation that bears on it |
|---|---|---|
| printed solutions +1.7°, +2.5° from Neptune | `fit` `targets.S5/S6.places.1846-09-23.helio_true` minus `neptune_places.1846-09-23.helio_true` | none |
| fit at his distances 10.5°, 9.1° short | `repro` `R7_d4.variants.h1_g/h2_g.true_lon_1846_of_date` minus `R7_d4.neptune_same_model_of_date` | none |
| Neptune inside the 68 % range with the distance free | `fit` `reading.g.statement_1.F.status` | none |
| reproduced to 0.011°, 0.003° | `repro` `R3_stage_c.H1/C0, H2/C0.vs_print.rows.theta.delta` | none |
| 1.44°, 1.76° from his print; 0.02°, 0.01° from Brookes | `repro` `R3_stage_c.H1/C2, H2/C2.vs_print` and `.vs_brookes` `.rows.mean_lon_1846.delta` | none |
| 2.9°, 1.3° per unit of the fourth decimal | `repro` `R3_locate.H1/H2.n_coefficients.dy1.dtheta_per_step` | none |
| every method 7.8° to 10.4° short on a noise-free Neptune | `repro` `R6_d3.O_i_d1c`, `.O_i_his_procedure`, `.O_ii`, `.iii.O_iii_median` | none |
| his data through his least squares, +1.6° | `noise` `multistart_baselines.H1.real.O` | none |
| pre-registered 0.183, 0.074: "noise suffices" | `noise` `readings.H1/d1c/primary` | the statistic is linear in an angle and measures branch switching |
| post hoc 0.017, 0.011 (0.046, 0.037): "weak" | `noise` `post_hoc_readings.H1/d1c/primary` (`H2/…`) | post hoc, defined after the results |
| second minimum: e′ 0.45, 165° from Neptune, Δχ² 6.9 | `noise` `multistart_baselines.H1.second_minimum.real` | found during the run |
| Fig. 2's histograms | `P1B_adams_noise_hist.json` | binned from the cached draws after the results; the summary JSON is not edited |
| 100,000 kept of 1,003,853 | `det` `n_accepted`, `acceptance` | none |
| Table 1's probabilities and ranges | `det` `cases.<case>.held.primary` | his actual procedure is 0 by construction (the realised events are held fixed, so the cases measure what a change would have added); the paired rule is uninformative |
| "depends on" and its input | `det` `cases.<case>.reading_row`, `.top_input` | the ranking of inputs used 10,000 draws |
| prior-dependent rows | `det` `cases.<case>.prior_dependent`, `.halved_priors` | zone priors not halved; 20,000 draws |
| replay, his actual procedure 0.34 | `det` `cases.A0.replay.primary.mean` | outside the pre-registered reading |
| pairs of counterfactuals | `det` `post_hoc_pairs_description_only` | post hoc, description only |

## References for the pre-discovery recordings

The companion page's § 7 lists every reference. The three short communications behind the
count of earlier recordings, with their full titles:

- J. R. Hind, "Auszug aus einem Briefe des Herrn Hind", *Astronomische Nachrichten* 30,
  257–258 (1850), DOI 10.1002/asna.18500301706 (Lamont's recordings of 25 October 1845 and
  7 September 1846). The Crossref record appends the letter's dateline (Mr. Bishop's
  Observatory, Regents-Park, London) to the title; it is not part of it.
- V. Mauvais, "Sur une observation inédite de la nouvelle planète", *Comptes rendus* 24,
  666–669 (1847) (Lalande's recordings of 8 and 10 May 1795).
- E. Schönfeld, "Notiz über eine Beobachtung des Neptun in den Lamont'schen Zonen",
  *Astronomische Nachrichten* 114, 59–60 (1886), DOI 10.1002/asna.18861140406 (Lamont's
  recording of 11 September 1846).

Titles and page ranges: Hind and Schönfeld from the Crossref records of the journal; Mauvais
from the heading on p. 666 of the scanned volume.

## What is not claimed

- That Adams was lucky, or that his nearness has a systematic cause.
- That Challis ought to have done otherwise. The model covers his nights and his instrument,
  conditional on the printed record and the declared priors; the Cambridge observing books
  were not consulted.
- Any probability for the British search as a whole.
