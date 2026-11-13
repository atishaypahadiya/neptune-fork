# S1 cell 0: the low-distance solution of the eccentricity test, by sub-cells

**Pre-registration, 13 November 2026.** Written before any code of this plan exists and
before any sampling in a sub-cell. Everything below the line at the end is added after
this plan is pushed: amendments and results, dated, never edits to the text above it.
The draft is `docs/notes/V13_drafts.md`, "S1 cell 0"; § 7 lists how each of its open
choices is resolved.

## 0 Scope and question

S1 is the sensitivity run of M6c with the eccentricity allowed up to 0.5. Its reading rule
(`M6.md`) is on P(a < 25 AU), the posterior weight of the low-distance solution, which is
the evidence share of cell 0 (a from 20 to 25 AU): above 0.5 "moved", 0.05 to 0.5 "moved
partly", below 0.05 "not at all". The runs of cell 0 so far do not agree (`M8.md`, "A8d(S1)
arbiter: results"):

| run of cell 0 | ln Z₀ | P(a < 25 AU) |
|---|---|---|
| published, one full-likelihood run, nlive 2000 | −366.84 | 0.001 |
| fast likelihood, four seeds merged, reweighted | −362.56 | 0.109 |
| full likelihood, nlive 4000, 28 slices, seed 1 | −361.03 | |
| the same, seed 2 | −362.84 | |
| those two merged | −361.94 | 0.123 to 0.186 |

The published verdict "not at all" stands and the question is recorded as unresolved,
0.001 to about 0.2.

**Question: cut into sub-cells that are sampled separately, do repeated runs of cell 0
agree, and if they do, what are ln Z₀ and P(a < 25 AU)?**

Not in scope: the headline fit (e ≤ 0.3), its intervals, any other cell of S1, any other
sensitivity run, the site.

## 1 Inputs, fixed

| input | file | SHA-256 |
|---|---|---|
| observations | `data/historical/leverrier_1846.csv` | `2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4` |
| error model | `data/historical/m6_error_model.json` | `0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277` |
| reference Uranus for the fast likelihood | `docs/notes/M8_u_ref.json` | `4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5` |
| fast likelihood | `src/fork/model/perturbative.py` | `092f673b9352f9d2b0ae709f2a4d3f258dd09806245b767fd83b86e8ecccd68c` |
| full likelihood with Uranus marginalised | `src/fork/model/linearized.py` | `3e02e78bdd9ee74fb96a6666f192188c355113addc92f43e376729f0b616d15e` |
| reweighting | `src/fork/inference/reweight.py` | `d0c5fbff6066ebaa531bc24668fa807a8f7bf0aa96cf77eead5c07aca9cd021b` |
| the other cells, published | `site/data/m6_fits.json` | `2d111e7fcb7a71100a72a8f963f6c8f0f26aea9da90ff22a246893618bb74fff` |
| the other cells, fast likelihood | `site/data/m8_runs.json` | `4d50867272ab7ce6faec99ef3db1edb30e6b25ba61223f8adc0304c65f88b6cb` |

The runner refuses to start if any of the first six differs. The data set, the priors and
the space are S1's, unchanged: `m6c.dataset("s1_e05")`, the free-distance space with
e ≤ 0.5 and an inclined perturber.

**The other cells' ln Z** (cells 1, 2, 3), used only to turn ln Z₀ into P(a < 25 AU):
published −366.0354, −359.7935, −365.1071; fast likelihood −365.682, −360.290, −365.258.
Prior masses of the four cells: 0.203114, 0.103156, 0.166495, 0.527235.

## 2 The sub-cells

Cell 0 is a in [20, 25] AU, e in [0, 0.5]. It is cut into six sub-cells, the product of

- a: [20, 22.5) and [22.5, 25] AU;
- e: [0, 0.2), [0.2, 0.35) and [0.35, 0.5].

They are numbered 0 to 5 in the order of `fork.export.v13.S1_SUBCELLS` (a outer, e
inner). Their prior masses are 0.042884, 0.032163, 0.032163, 0.038361, 0.028771, 0.028771
and add up to cell 0's 0.203114. Each sub-cell's prior is cell 0's, renormalised on the
sub-cell (`nested.cell_space`).

The cut is the draft's example, taken unchanged. It was fixed from published summaries
only: the merged full-likelihood posterior has a = 24.08 AU (23.51 to 24.60) and
e = 0.324 (0.305 to 0.341). No sample of any earlier run was looked at to place it. It is
not tuned to separate the two longitudes the two full-likelihood seeds found (313.3° and
315.0°); if both lie in one sub-cell the cut may not help, and the test of § 4 will say so.

## 3 The runs

**Stage A: the fast likelihood, every sub-cell, four seeds.**

- Likelihood: the perturbative likelihood of M8, as pre-registered there (second order
  off, step 60 d, cubic quadrature), over the perturber's 7 parameters, Uranus
  marginalised.
- Sampler: dynesty static, bound `multi`, sample `rslice`, nlive 4000, 28 slices, dlogz
  0.1, no iteration limit. One process per run, checkpoint every 300 s.
- Seeds: rng `[2026, 14, 1846, 1, k + 1, s]` for sub-cell k and seed s = 1 to 4.
- 24 runs.

**Reweighting to the full likelihood.** For each run and for each sub-cell's four-seed
merge (30 posteriors): 6,000 equal-weight draws (rng `[2026, 14, 1846, 1, k + 1, s, 7]`,
s = 0 for the merge), each evaluated with the fast likelihood and with the full
likelihood with Uranus marginalised (mode (b) of M8), in chunks of 250 (chunk rng
`[..., s, 9, c]`). The reweighted evidence is ln Z = ln Z_fast + ln mean exp(Δ), with Δ
the difference of the two log-likelihoods, as `reweight.reweight_cell`. Galle-night
longitudes are computed through the full model for each draw as in M8.

**Stage B: the full-likelihood check.** In the sub-cell k* with the largest w_k Z_k among
the merged, reweighted sub-cells: two runs with the full likelihood of M6c (11
dimensions, Uranus sampled), nlive 4000, 28 slices, dlogz 0.1, no iteration limit, 6
workers each, rng `[2026, 14, 1846, 1, k* + 1, 100 + r]`, r = 1, 2. The draft asks for one
run; two are made because the earlier full-likelihood check of the whole cell failed on
the disagreement of its own two seeds, which one run cannot show.

**Quantities, per posterior:** ln Z and its error (the sampler's, times 1.8 where a
reported error is formed, as everywhere since M6c); the reweighting's effective sample
fraction ESS/N; the median and the 68 % width of the perturber's geocentric longitude on
Galle's night (`m6c.lon_summary` on resampled draws for Stage A, `m6c_repeat.run_stats`
for Stage B).

**Evidence of cell 0:** Z₀ = Σ w_k Z_k / Σ w_k over the six merged, reweighted sub-cells
(`v13.combine_subcells`), error propagated to first order from the six errors × 1.8.

**P(a < 25 AU):** cell 0's share with M6c's prior masses, once against the published
other cells and once against the fast-likelihood other cells
(`m8_arbiter.implied_p_below_25`). The two check values of M8, 0.00104 and 0.109, must be
reproduced by the same function before it is applied.

## 4 Gates and tests, fixed in advance

**Gates** (if one fails, there is no replacement and the failure is reported with its
numbers):

- G1: all 24 runs of Stage A and both runs of Stage B stop on the dlogz criterion.
- G2: ESS/N ≥ 0.3 for each of the 30 reweighted posteriors (M8's gate 2).
- G3: the truncation of the marginal over Uranus is below 10⁻⁶ in −2 ln L at every draw
  where the full likelihood is finite (M8's gate 4).

**T1, the seeds agree in every sub-cell.** With the four seeds' reweighted values:

- ln Z: largest minus smallest ≤ 0.5, in every sub-cell.
- Galle-night median: largest minus smallest ≤ 0.35°, in every sub-cell whose share of
  cell 0's evidence, w_k Z_k / Σ w_j Z_j from the merges, is at least 0.01. The result of
  the same test in the remaining sub-cells is reported beside it and not read: a sub-cell
  with under 1 % of the evidence cannot move P(a < 25 AU) by more than 1 % of itself,
  and a flat posterior there has no stable median.

**T2, the full likelihood agrees with the fast one in k*.**

- The two full-likelihood runs agree with each other: ln Z within 0.5, medians within
  0.35°, each 68 % width within 10 % of their mean.
- Their merge (`dynesty.utils.merge_runs`) against the fast four-seed merge of k*,
  reweighted: median within 0.35°, 68 % width within 10 %, ln Z within 0.5 (A8b's
  tolerances).

## 5 Reading, fixed in advance

| outcome | reading |
|---|---|
| G1 to G3, T1 and T2 all hold | **resolved.** ln Z₀ is the sum of § 3 and P(a < 25 AU) its two values. S1's rule is applied to each. If both give the same word, that word replaces "unresolved", and S1's published verdict is corrected by an erratum line that keeps the old number visible. If the two give different words, both are reported and the verdict is "between" the two words |
| a gate fails | **not resolved, the method failed.** The published verdict stands, the question stays unresolved, the failed gate is reported with its numbers |
| T1 fails | **not resolved, the sub-cells do not make the runs repeatable.** Reported per sub-cell. The published verdict stands |
| T1 holds and T2 fails | **not resolved, the fast likelihood and the full one disagree in k*.** Both are reported. The published verdict stands |

A value within 0.02 of a tolerance (ln Z), 0.02° (median) or 0.5 points (width) is
reported as "at the tolerance" and counted as it falls.

**The 95 % range of the headline.** Its "provisional" label was set by the fourth
condition of M8b's rule (the spread of each 95 % bound over four single-seed fits of the
headline, e ≤ 0.3; `M8.md`). This plan does not run the headline and does not test that
condition. What this run can change is fixed here:

- The label may be dropped only if both of these hold: (i) this plan resolves with
  P(a < 25 AU) < 0.05 against both sets of other cells; (ii) the fourth condition of
  M8b's rule is met. (ii) is not met and is not tested here, so **the label stays
  "provisional" under every outcome of this plan.**
- If (i) holds, the low-distance solution is no longer a reason for the label, and the
  only one left is (ii).
- If this plan resolves with P(a < 25 AU) ≥ 0.05, or does not resolve, the label keeps
  both reasons: (ii), and a solution near 314° under e ≤ 0.5 whose weight is the resolved
  P(a < 25 AU), or still 0.001 to about 0.2.

No interval of the headline is recomputed. The site is not touched by this plan.

## 6 What this cannot show

- Anything about the headline's posterior: S1 is a different prior (e ≤ 0.5).
- That the fast likelihood is right outside k*: it is reweighted everywhere, and checked
  against full-likelihood sampling in one sub-cell only.
- That a different cut would give the same answer. One cut is run.
- Whether the low-distance solution is physical. It is a statement about the evidence
  under this prior and error model.

## 7 The draft's open choices, each resolved

| open in the draft | resolved as | why |
|---|---|---|
| the cut | the draft's example, six sub-cells | it is the one the code and its tests already hold, and it was written before any run |
| slices | 28 | the setting of the earlier full-likelihood check and of the higher-resolution reruns |
| "seeds agree within A8b's tolerances" | range of four seeds ≤ 0.5 in ln Z in every sub-cell; ≤ 0.35° in the median where the sub-cell holds at least 1 % of the evidence | the range is the strictest reading for four values; the 1 % rule is stated with its reason in § 4 and the unrestricted result is reported |
| one full-likelihood run | two, with their mutual agreement tested | § 3 |
| the alternative sampler (parallel tempering) | not run | it needs its own validation first and is the larger job |
| which other cells | both sets, as the earlier check | no reason to prefer one |

## 8 Order, cost, layout

- Order: this plan pushed; the runner and its tests written; a smoke run on synthetic
  data in a scratch cache (300 iterations, two sub-cells, nothing kept); Stage A; the
  reweighting; Stage B in k*, chosen by the rule of § 3 by the runner, with no reading in
  between; the summary.
- Cost, from the ledger: a fast-likelihood run at nlive 4000 and 28 slices took 0.97 to
  1.14 CPU-h per cell of the headline, so Stage A is about 24 to 30 CPU-h, 2 to 2.5 h on
  12 processes. The reweighting is about 2.5 CPU-h. A full-likelihood run of the whole
  cell at these settings took 47.8 and 49.7 CPU-h; a sub-cell is at most that, so Stage B
  is at most about 100 CPU-h and 10 h on 12 workers. In all at most about 135 CPU-h.
- Runs are detached, resumable from checkpoints, logged, one ledger row per job.
- Layout: `src/fork/export/s1_cell0.py`, `tests/test_s1_cell0.py`, cache
  `data/cache/nested/s1_cell0/`, results `docs/notes/S1_cell0.json`.

---
