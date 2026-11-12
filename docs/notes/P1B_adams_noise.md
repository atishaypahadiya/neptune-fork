# Phase 1b, step 3c: is the shift between clean data and Adams's real-data solution within noise?

Pre-registered plan. Written 2026-11-11 after step 3b's results (commit 4a3301d); nothing below
§ 1–7 runs until the plan is committed, and the pre-registered sections are then byte-identical
(amendments and results only below the lines). Code and tests cite this note.

## Pre-registration (before any code or run)

### 0 Where step 3b left it

On one noise-free Neptune truth, his layout (30 normal places at his epochs) and his method
(the eight-parameter least squares with his § 52 weights, "D1(c)") at α = 0.5 placed the planet
at a heliocentric true longitude of 319.66° (ecl_j2000) for 1846.762, **9.51° short of
Neptune's 329.17° in the same model** (his own elimination chain: 9.92° short). His real data
through the same method (C2, which is Brookes's 1970 correction of Hypothesis I) place it at
328.05° of date against Neptune's 327.03°, **1.02° past**; his printed solution 328.77°, 1.74°
past. Step 3b read this as "his nearness on real data owes to his data or to compensating
errors" (Q3 row 3) and named the candidates A1–A5 and his reductions without ranking them.

### 1 Question

Is the shift between the clean-data solution and his real-data solution, **each method against
its own clean baseline** on the 1846.762 heliocentric true longitude,

- D1(c): Δ_obs = O(D1(c) on his real § 9 data, step 3b R4: 328.68° − 327.03° = +1.65°) − O_clean
  (−9.51°) = **+11.16°** (Hypothesis I);
- his chain: Δ_obs = O(C2, his chain on his real data = Brookes's correction: +1.02°) − O_clean
  (−9.92°) = **+10.94°**;
- secondary target, his printed solution: +1.74° − (−9.51°) = +11.25° against D1(c)'s baseline;
- Hypothesis II at α = 0.515: the same two differences from step 3b's D1(c) (332.00°, +4.97°) and
  C2 (331.48°, +4.45°) and the noise-free (i) solutions at 0.515, which step 3b did not compute
  and the pre-run record will (same code, his t_h2 and § 52 Hypothesis II weights),

within what independent observational noise produces in his layout through his method, or does
it require a systematic? "Systematic" means anything a smooth error in the data could be:
Bouvard's tables and his § 7 corrections, the construction of the normal places, a reduction
error common to a series. The step does not attribute further than that (Bouvard's tables are
not in the cache).

### 2 Design

Start from step 3b's noise-free his-layout normal places (D3(i): the truth's Uranus longitudes
at his 30 epochs minus the Bouvard-base "C", converted by r²/(ab); `cheap.json`'s
`d3.synthetic_places.oc_mean_arcsec`, recomputed by the same code and checked against it to
10⁻⁶″). To each of the 29 places used (1690 set aside as he did) add independent Gaussian noise,
ε_j ~ N(0, σ_j²), one draw per place per realisation; solve; repeat **10,000 realisations per
noise model**, seed `np.random.default_rng([2026, 11, 3, 8, model, hyp])` (model 0 = (a),
1 = (b); hyp 1, 2); the same 10,000 noise vectors (in units of σ_j) are reused across methods
within a model so that the two methods are compared on identical data.

**Two noise models, bracketing:**

- **(a) his own residual scatter.** σ_j = RMS of his printed § 52 residuals, modern and ancient
  separately, **scaled by √(N/(N − p)) = √(29/21) = 1.175** (N = 29 places, p = 8 parameters, the
  one factor applied to both groups): Hypothesis I raw 1.24″ (21 modern) and 6.21″ (8 ancient),
  scaled 1.46″ and 7.30″; Hypothesis II raw 1.22″ and 6.29″, scaled 1.44″ and 7.39″. The raw
  values are noted in every table. **Conditioned on his own fit**: residuals of his solution, so
  they understate the noise by more than the factor if his model is not exact; said in every
  table.
- **(b) set D's scatter propagated to 3-year normal places.** Per group, the RMS of step 3's
  H1(G) best-fit residuals about the fitted group offset, RA measurements only (his normal places
  are longitudes): ancient 5.48″ (18), Greenwich 2.22″ (87), Paris 1835–1845 2.18″ (40). Per
  place, σ_j = RMS_g / √N_j with N_j the number of set D observations whose epoch lies within
  ±1.518 yr (half his interval) of a modern place or within ±1 yr of an ancient one, and g the
  group those observations belong to. The counts from set D (each observation falls in exactly
  one window): ancient 1712: 1, 1715: 3, 1750: 2, 1753: 1, 1756: 1, 1764: 1, 1769: 8, 1771: 1;
  modern 1783: 19, 1786: 6, 1789: 15, 1792: 14, 1795: 11, 1798: 6, 1801: 2, 1807: 2, 1810: 1,
  1828: 6, 1831: 5 (Greenwich); 1834: 13, 1837: 17, 1840: 10 (Paris). **Seven modern windows
  hold no set D observation** (1780, 1804, 1813, 1816, 1819, 1822, 1825; Table 77 carries no
  Greenwich rows for 1811–1827 and set D excludes Paris 1801–1827), so N_j is **floored at 1**
  there and σ_j is the single-measurement scatter of Greenwich (2.22″). This makes (b) wider
  than the truth at those places (his normal places there rest on observations Table 77 does not
  carry), which biases (b) toward "noise suffices"; a "systematic required" reading that survives
  (b) is therefore the robust one. (b) treats his normal places as means of set D's measurements,
  which they are not exactly (Bouvard's and Greenwich's reductions); that is the point of
  bracketing with (a).

Neither model carries a group offset (step 3's τ_g): a common shift of a whole series is a
systematic of the kind the question sets against noise, not noise, and is the first candidate
if the reading says "systematic required".

**Sensitivity (c), reported and not read:** Student-t noise with ν = 4 and model (b)'s scales
(ε_j = σ_j · t₄), seed `[2026, 11, 3, 8, 2, hyp]`, the same 10,000 realisations through both
methods; it shows what heavier tails do to P and enters no reading row.

**Methods** (both at α = 0.5 for Hypothesis I, 0.515 for II, the printed series constants and
§ 28/45 forms, as D3(i)):

- primary, **D1(c)**: the eight-parameter weighted least squares (weights 1/RMS² of his § 52
  residuals, as in step 3b) started from the noise-free solution (one start: the method that gave
  −9.51°; a warm start is what an iterated solution does), m′ ≥ 0 and e′ ≤ 1 enforced;
- secondary, **his own chain** (C2 run on the noisy column: his normal equations, eliminations,
  § 25 multipliers, logarithmic equations, root choice by the smallest D1(c)-weighted χ² among the
  positive-mass roots, as D3(i) did).

Per realisation: θ, m′, e′, β′, the mean and the true longitude at 1846.762 (ecl_j2000, the
data's frame; the offset is frame-free), O = true longitude − T (T = 329.17°, the truth's
Neptune), and Δ = O − O_clean, where O_clean is the noise-free solution's offset by the same
method (−9.51° / −9.92°; Hypothesis II's from the pre-run record). **Escapes are outcomes**: a
draw that converges to another minimum is kept with its actual Δ (counted as `escaped`, θ more
than 30° from the clean θ, for description). Only non-convergence (`failed`: no convergence of
the least squares, or no positive-mass root of the chain) is excluded; failures are counted,
and if they exceed 1 % a second P that counts every failure as exceeding Δ_obs (an upper bound)
is stated beside the first. `bound` (m′ or e′ on a bound) is counted, not excluded.

**Multistart check.** For the first 1,000 realisations of each model × hypothesis, D1(c) is also
solved from a grid of 12 θ starts (every 30°, m′ 0.8, e′ 0.1, β′ the clean value) and the lowest
χ² kept; the fraction of realisations where this global minimum differs from the warm-start
solution by more than 1° in the 1846.762 true longitude is reported. **If it exceeds 1 %, all
draws are rerun with the multistart before the reading** (the warm start then stays as a
description of what an iterated solution would have done).

### 3 What is computed and reported

For each hypothesis × noise model × method: the distribution of Δ (quantiles 2.5, 16, 50, 84,
97.5 %, the mean and SD), the distribution of O, of the mean longitude, of m′ and e′ (as
description), the flagged fractions, and

**P = P(Δ ≥ Δ_obs), one-sided**, with its binomial standard error, for the primary target
(C2) and the secondary (his printed solution); the number of exceeding draws stated (at
P = 0.003 that is 30 of 10,000, so the threshold is resolved).

P is over the converged realisations (escapes included). Also: the SD of Δ against the
least-squares formal σ of the clean solution (step 3b: σ_θ 1.49° for (i)), as a check that the
noise propagation and the covariance agree under (a); and the multistart fraction of § 2.

### 4 Reading, fixed in advance

Read on Hypothesis I, primary method (D1(c)) against its own Δ_obs (+11.16°); his chain against
its own (+10.94°), Hypothesis II and the printed target reported beside, and any disagreement
between the two methods said per row. Sensitivity (c) enters no row.

| outcome | reading |
|---|---|
| P ≥ 0.05 under **both** (a) and (b) | **Noise suffices**: a shift as large as his is produced by independent observational noise in his layout at least one time in twenty; his nearness to Neptune on real data is within luck, and no systematic is needed to explain it (none is excluded, see below). |
| P < 0.003 under **both** | **A systematic is required**: independent noise of either size does not carry his method 10.5° from where clean data put it. Named candidates, not ranked and not attributed further without Bouvard's tables: Bouvard's 1821 tables and his § 7 corrections to them (what the "C" he subtracted contained), the construction of the normal places (which observations, which weights, which reductions), a reduction error common to a series. |
| otherwise | **Weak**, stated per model: the P under each, which model admits it and which does not, and no single statement. |

**What independent noise cannot do.** The draws are independent per place; a smooth systematic
(a slow error in Bouvard's elements, a series-wide offset) is correlated across places and is
exactly what the noise models do not contain, so "noise suffices" does not exclude a systematic:
it says one is not needed. Conversely "systematic required" does not say which; the candidates
above are where to look, and the first of them needs Bouvard's tables, which this study does not
have. The reading is about the real-data solution's distance from the clean one, not about the
truth: both solutions are measured from Neptune only because that is the fixed point step 3b
used.

### 5 Report order

R1 the σ_j tables for (a) (raw and scaled) and (b) with the counts and the frozen inputs'
hashes; R2 the clean solutions and Δ_obs per hypothesis, method and target; R3 the Δ
distributions (table and figure: histograms of Δ per model and method with each method's Δ_obs
marked, Hypothesis I above II, (c) beside); R4 P with standard errors and exceeding counts; R5
the failed, escaped and bound fractions, the multistart fraction and the SD–covariance check;
R6 the reading rows; R7 what remains unseparated. Results appended below the lines of this note;
`docs/notes/P1B_adams_noise.{json,png}`; a ledger row `p1b_adams_noise/draws`; the site
untouched; nothing in step 3b's note edited.

### 6 Cost and order of work

One warm-start least-squares solve takes 25 ms, one run of his chain 34 ms (measured on the
D3(i) column with noise). 10,000 realisations × 3 noise models × 2 hypotheses × 2 methods is
60,000 + 60,000 solves, plus 12 × 1,000 × 6 multistart solves, about 25 + 34 + 30 CPU-minutes,
run in a process pool (12 workers) in about 8 minutes wall; **minutes, not seconds; not
detached**, no checkpoint (a rerun is cheaper than one). A multistart rerun of all draws, if
§ 2 calls for it, adds about 5 CPU-hours (12 × 60,000 solves), still not detached. Order:
the pre-run record (σ_j tables, the Hypothesis II clean solution at 0.515, Δ_obs), the tests, the
draws, the summary and figure, the results below the lines; one commit for the plan, one for the
pre-run record and code, one for the results; no trailers.

### 7 Layout

`src/fork/historical/adams_1846.py` gains `noise_sigmas(m, model, counts)` and
`solve_warm(rows, w, s, x0)` (the one-start least squares with the flags); the runner
`src/fork/export/p1b_adams_noise.py` (`prerun`, `draws`, `summarize`), caches in
`data/cache/p1b_adams_noise/`; tests in `tests/test_p1b_adams_noise.py`: the seeds, the noise-free
draw reproduces D3(i)'s solution to 10⁻⁶, the σ tables and counts, the reading rows on
constructed P, the frozen inputs (`cheap.json`'s synthetic column is regenerated, not read, and
compared).

### 8 Amendments made on approval (2026-11-11, before the commit; nothing runs before it)

1. Escapes are outcomes: draws converging to another minimum are kept with their actual Δ; only
   non-convergence is excluded (counted; upper-bound P if > 1 %). Warm start stays primary, with
   a multistart check on 1,000 draws per model × hypothesis (12 θ starts, lowest χ² kept; the
   fraction differing from the warm start by > 1° in the 1846.762 longitude reported); if > 1 %,
   every draw is rerun with the multistart before the reading.
2. Δ_obs per method against its own clean baseline: D1(c) on his real § 9 data (+1.65°) minus
   D1(c) clean (−9.51°) = +11.16°; his chain on real data (C2, +1.02°) minus his chain clean
   (−9.92°) = +10.94°. Both reported; the reading uses each method's own.
3. Model (a): σ scaled by √(N/(N − p)) = 1.175, the one factor on both groups; raw RMS noted.
4. Sensitivity (c): Student-t ν = 4 with model (b)'s scales, same seed namespace; reported, not
   part of the reading rule.
5. The empty-window count corrected: seven modern windows (1780, 1804, 1813, 1816, 1819, 1822,
   1825), not eight.
