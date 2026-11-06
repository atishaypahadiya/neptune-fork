# Phase 1 spec: rebuild Neptune

**Goal:** from Uranus's heliocentric longitudes (1690–1845), with Neptune absent from the model,
infer a posterior over an unseen perturber and render where it sits in the sky on
**1846-09-23** (Galle's discovery night, Berlin), versus the true Neptune.

**Done when:** acceptance tests A1–A5 pass and the site is live.

---

## M0: Scaffold
- `uv init`, package `fork` under `src/`, ruff + pytest config, `.gitignore` (`data/cache/`).
- GitHub Actions: `test.yml` (fast tests on push), `pages.yml` (build + deploy `site/` on main).
- README: one-paragraph pitch, prior-art credits, build status badge.
**Exit:** `uv run pytest` runs (zero tests OK), CI green.

## M1: Frames and data access
- `data/horizons.py`: fetch heliocentric state vectors (ecliptic J2000, AU, AU/day) via
  `astroquery.jplhorizons` for Sun-relative Jupiter, Saturn, Uranus, Neptune, Earth at any JD list.
  Cache to `data/cache/*.npz`; never refetch if cached.
- `data/frames.py`: precession between mean ecliptic/equinox of date and J2000 (use astropy
  `GeocentricMeanEcliptic` / `HeliocentricMeanEcliptic` with `equinox=`), plus state-vector ↔
  orbital-elements conversion.
- Constants in `model/bodies.py`: mass ratios (Sun/planet) Jupiter 1047.348644, Saturn 3497.9018,
  Uranus 22902.98, Neptune 19412.26. Inner planets folded into the Sun's mass.
**Tests:** elements→state→elements round-trip to 1e-12; precession of a fixed direction over 100 yr
≈ 5029″ (±1″).

## M2: Kepler solver
- `integrators/kepler.py`: universal-variable two-body propagator (f and g functions, Stumpff
  functions), Newton iteration with bisection fallback. Must handle e≈0 and long Δt.
**Tests:** propagate one full period returns to start (< 1e-12 AU); agrees with an analytic ellipse
over 1000 random cases (hypothesis); e from 0 to 0.95.

## M3: Integrators
- `rk4.py` (baseline, not symplectic), `leapfrog.py` (KDK), `wisdom_holman.py`
  (democratic heliocentric coordinates). The default composition is REBOUND WHFast-DH's, with the
  Kepler drift outermost: Kepler drift (h/2) → Sun-momentum drift (h/2) → interaction kick (h) →
  Sun-momentum drift (h/2) → Kepler drift (h/2). The kick-outermost composition (half-kick →
  Sun-momentum drift → Kepler drift → Sun-momentum drift → half-kick) stays available as
  `wh_kick_outer`. All `numba.njit`, signature
  `integrate(x0, v0, masses, t_out, dt) -> positions[t_out, body, 3]`.
- Default dt = 12 days (`fork.integrators.DEFAULT_DT_DAYS`). *Amended:* the original 20 d, with
  the kick-outermost default, gave 0.049″ against REBOUND and 4e-8 energy error; 14 d then met
  the 1e-8 energy bound by only 0.7%, so the default was lowered to 12 d. See `docs/notes/M3.md`.
**A1 (acceptance):**
  - Sun+Jupiter+Saturn+Uranus from Horizons @ JD of 1781-03-13: Uranus heliocentric longitude agrees
    with REBOUND WHFast (same ICs, same dt) to **< 0.01″** over 1781–1846.
  - Relative energy error over 1000 yr stays bounded **< 1e-8** for WH; RK4 shows secular drift
    (save the comparison plot; it becomes a site panel).
  - WH run of 65 yr completes in **< 10 ms** after JIT warmup.
  - Convergence order (permanent test): Uranus position error vs a dt = 1 d reference drops by
    a ratio in **[3.5, 4.5]** per halving of dt (24 → 12 → 6 → 3 d), for both compositions.
    The data is exported to `site/data/wh_convergence.json` for a future panel.

## M4: Forward model and synthetic data
- `model/forward.py`: `predict_uranus_lon(theta, t_jd)` where
  θ = Uranus elements (a, e, ϖ, λ₀, plus i, Ω held at Horizons values in Phase 1 planar runs)
  + optional perturber (log10 m [M⊕], a, e, ϖ, λ₀), planar.
  Initial conditions are relative to the Sun + inner-planet barycenter, since the inner planets
  are folded into the central mass (`docs/notes/M4.md`).
  Holding i, Ω at their epoch values is checked by a leak test: Neptune's effect on Uranus's i, Ω
  over 1781–1846 (with vs without Neptune, same ICs) must be < 1″ in longitude.
  Jupiter and Saturn start from Horizons at the epoch and are not fitted.
- `data/synthetic.py`: true Uranus longitudes from Horizons (which include real Neptune) at the
  historical epochs; add Gaussian noise σ = 3″ (1781 onward) and σ = 12″ (pre-discovery sightings,
  1690–1771). Seeded.
**A2 (acceptance):** fit Uranus alone (no perturber) to synthetic data with `least_squares`.
Residual curve shows **peak-to-peak ≈ 130″ (±25″)** and, in its **post-1781** part, a **sign
change of the slope within 1815–1830** (the 1822 Uranus–Neptune conjunction). *Amended:* the
slope check was restricted to post-1781; over 1690–1845 the curve turns several times in the
sparse pre-discovery era.
**A2b (acceptance):** with the true Neptune planted (fixed, not fitted) and Uranus refitted, the
post-1781 residual RMS is **≤ 1.5 σ_noise** (4.5″); without Neptune the same check must fail.

## M5: Inverse problem on synthetic data
*Split:* **M5a** = the point fits below (no MCMC); **M5b** = MCMC, derived quantity and A3.
- `inference/fit.py`: multistart Levenberg–Marquardt over a 24-point grid in perturber λ₀ (0–360°),
  two modes: (a) **Le Verrier mode**, a fixed at 38.4 AU (Bode's law); (b) **free a**.
  Free: Uranus a, e, ϖ, λ₀ (i, Ω held, rule 2) and perturber (log10 m, a or fixed,
  √e·cos ϖ, √e·sin ϖ, λ₀). M5a reports per mode: best m, a, e, χ²/dof, post-1781 residual RMS,
  number of distinct local minima, and the 1846-09-23 geocentric longitude vs true Neptune;
  exports residuals after fit to `site/data/m5a_fits.json`. Results: `docs/notes/M5.md`.
- `inference/mcmc.py`: emcee, 32+ walkers, sample in
  (log10 m, log a, √e·cos ϖ, √e·sin ϖ, λ₀, Uranus elements).
  Priors: log10 m ∈ [0, 2] (1–100 M⊕), a ∈ [20, 60] AU, angles uniform. The perturber's
  eccentricity is sampled as (k, h) = (√e·cos ϖ, √e·sin ϖ), not as (e, ϖ): a flat prior on the
  (k, h) disk is exactly uniform in e and ϖ (dk dh = ½ de dϖ), and ϖ stays well defined as e → 0
  (see `docs/notes/M1.md`). Reject (log-prior = −∞) any proposal with e = k² + h² > 0.3.
  Report autocorrelation time; run ≥ 50τ. Save chains to `data/cache/chains/`.
  *Amended (M5a review):* point fits also enforce the prior bound e ≤ 0.3. M5b samples free a
  under **two priors**, e ∈ [0, 0.3] and e ∈ [0, 0.1], and reports prior sensitivity as a
  headline result. It also samples **Le Verrier's shape**: a, e, ϖ held at his published 1846
  values, with m, λ₀ and Uranus free. It reports Δ log-likelihood against free a, written up in
  `docs/notes/M5_le_verrier.md`. MCMC runs locally (`python -m fork.export.m5b`). CI runs only
  a tiny smoke test. A summary JSON (`site/data/m5b_posteriors.json`) is committed; raw chains
  are not.
  *Amended (M5b review, multimodality):* on seed 0 the free-a posterior has two modes. Mode A has
  a ≈ 30 AU and e ≲ 0.05. Mode B has a ≈ 22 AU and e ≈ 0.28, against the e ≤ 0.3 bound, and
  holds the best χ² (71.22 vs 71.91). An emcee chain started in either mode stays there for
  20–40k steps while passing the 50τ test, so emcee reports the mode it started in. **M5b
  posteriors therefore use dynesty** (`inference/nested.py`): static nested sampling,
  `bound='multi'`, `sample='rslice'`, `nlive=1000`, stop at Δln Z < 0.1, with the perturber's
  ϖ and λ₀ periodic. It reports ln Z ± error and the mode masses, found by splitting the
  weighted samples at a = 25 AU. Phase 2's evidence needs it anyway. emcee stays, for checks
  within one mode.
  Validation gate, before the full run: on seed 0 (free a, e ≤ 0.3), both modes are recovered,
  with within-mode medians and 68% points matching the two emcee chains to ≲ 0.5σ. A second rng
  seed gives ln Z within the stated error and mode masses within a few %.
  *Amended (dynesty validation, `docs/notes/M5b_dynesty_validation.md`):* unpartitioned dynesty
  failed the gate. ln Z was reproducible (−73.05 ± 0.22 vs −73.31 ± 0.22), but mode B was never
  found, and region C (a 24–28 AU) got 1% vs 39% of the mass. M5b therefore uses
  **partitioned nested sampling**:
  - the a prior is split into slices [20, 25], [25, 28] and [28, 60] AU, and each slice is
    sampled separately (same dynesty settings), with the prior renormalised on the slice;
  - the slice's prior mass is w_k = ln(a_hi/a_lo) / ln 3, and Z = Σ w_k Z_k;
  - σ(ln Z)² = Σ s_k² σ_k², where s_k = w_k Z_k / Z is the slice's posterior mass share
    (reported, with delta-method errors);
  - the posterior is the slices' weighted samples, each slice scaled by s_k.

  **Slice edges affect efficiency, not correctness.** Every slice is a valid restricted prior
  and the weights are exact, so any partition gives the same Z and posterior in the limit. The
  edges were chosen so that no sub-run has to discover a thin pocket next to a large mode. They
  are fixed in advance and the same for every dataset.
  Validation gate (seed 0, e ≤ 0.3, two rng seeds per slice):
  - per-slice ln Z_k and mass shares agree across rng seeds within error;
  - combined ln Z agrees with the unpartitioned runs;
  - slice [20, 25] recovers mode B and matches emcee chain B, and slice [28, 60] matches emcee
    chain A;
  - report the 1846 longitude median and 68% interval, combined and per slice.

  *Amended (partitioned validation round 1, `docs/notes/M5b_dynesty_validation.md`):* the
  three a slices failed the gate. [28, 60] AU holds mode A (e ≈ 0.03) and an eccentric branch
  (a ≳ 33, e 0.16–0.26) that the two rng seeds weighed 3% vs 18%, so its ln Z_k differed by
  1.2σ. The unpartitioned reference proved unreliable (mode B absent, [25, 28] at 1% vs 39%).
  **Partition axis for [28, 60]: circular vs eccentric**, because every alternative branch
  found so far is eccentric. The cells (`nested.PARTITION_CELLS`) are [20, 25] and [25, 28] AU
  (all e), and [28, 60] AU × e ∈ [0, 0.12) and × e ∈ [0.12, 0.3]. A cell's weight is
  w = w_a · (e_hi − e_lo)/0.3, with e uniform on the cell. Revised gate:
  - per-cell ln Z_k and mass shares agree across rng seeds within error;
  - combined ln Z and all cell shares agree across rng seeds within error (the unpartitioned
    comparison is dropped);
  - cell [20, 25] matches emcee chain B, both restricted to e > 0.25 (mode B, like-for-like),
    and cell [28, 60] × e < 0.12 matches emcee chain A, same tolerances;
  - report the 1846 longitude median and 68% interval, combined and per cell.

  *Amended (round 2, decision): the gate passes with measured errors.* **Evidence errors:**
  dynesty's quoted ln Z error understates the run-to-run scatter on these posteriors. On seed
  0 the pooled single-run SD of ln Z_k across rng seeds was 0.40 against a quoted 0.22
  (four cells, cell [20, 25] on three seeds: 5 degrees of freedom, so the factor's own 68%
  range is about 1.5–2.6). Every ln Z error that is reported or compared is therefore
  σ_eff = 1.8 · σ_quoted (`nested.ERR_INFLATION`), and share errors scale with it. **This
  applies to Phase 2 evidence too.** Seed-0 comparisons, quoted and corrected:

  | comparison | Δ/σ quoted | Δ/σ with ×1.8 |
  |---|---|---|
  | ln Z_k [20, 25], rng 1 vs 2 | 0.03 | 0.02 |
  | ln Z_k [25, 28], rng 1 vs 2 | 0.13 | 0.07 |
  | ln Z_k [28, 60] × e < 0.12, rng 1 vs 2 | 0.79 | 0.44 |
  | ln Z_k [28, 60] × e ≥ 0.12, rng 1 vs 2 | 1.39 | 0.77 |
  | share [20, 25], rng 1 vs 2 | 0.29 | 0.16 |
  | share [25, 28], rng 1 vs 2 | 0.17 | 0.09 |
  | share [28, 60] × e < 0.12, rng 1 vs 2 | 1.36 | 0.76 |
  | share [28, 60] × e ≥ 0.12, rng 1 vs 2 | 1.53 | 0.85 |
  | combined ln Z, rng 1 vs 2 | 0.44 | 0.25 |
  | ln Z_k [20, 25], rng 3 vs rng 1 / rng 2 | 3.08 / 3.12 | **1.71 / 1.73** |

  The rng 3 comparison stays above 1σ after correction: rng 3's ln Z_k is 1.0 higher with the
  same posterior shape, which points to the thin mode-B pocket being under-sampled at nlive
  1000. Hence nlive 2000 for that cell in the full run.

  **Full run** (`python -m fork.export.m5b run`, then `summarize`):
  - free a, e ≤ 0.3: noise seeds 0–9 × the four cells, one rng seed per cell; cell [20, 25] at
    nlive 2000, the others at 1000. The seed-0 validation runs of the other three cells are
    reused.
  - free a, e ≤ 0.1: seeds 0–4, only the cells that can hold e ≤ 0.1 (`nested.cells_for`):
    [20, 25] and [25, 28] AU, and [28, 60] AU × e ∈ [0, 0.1]. The [28, 60] × e ≥ 0.12 cell is
    outside this prior.
  - Le Verrier shape: seeds 0–4. The shape is fixed, so one unpartitioned run per seed, if seed
    0 passes a single-mode check: two rng seeds agree on ln Z within σ_eff and on the medians of
    log10 m and λ₀ within 0.5σ, and in each run < 0.5% of the mass lies beyond 8 half-widths of
    the median. If it fails, seeds 1–4 are not run until the shape is partitioned.
  - Noise-free truth, free a, e ≤ 0.3, the four cells: scheduled last, after everything A3b
    needs. The summary computes A3b first, then appends the noise-free check, a prior-bias
    measurement: truth inside the 68% intervals (1846 longitude, m, a, e), and the posterior
    median's offset from truth in σ (the 68% half-width on the side facing truth).

  Every long run writes a progress log and checkpoints, so it can be monitored and resumed.
- **Priors** (one definition, `inference/priors.py`, used by emcee and dynesty alike). The
  sampled vector is q = (Uranus a, e·cos ϖ, e·sin ϖ, λ₀; perturber log10 m, ln a, √e·cos ϖ,
  √e·sin ϖ, λ₀). The density is flat in q on the support below:
  - perturber log10 m uniform on [0, 2] (1–100 M⊕);
  - perturber a **log-uniform** on [20, 60] AU, i.e. ln a uniform. The emcee chains sample
    ln a with a constant density, so they used exactly this prior;
  - perturber e uniform on [0, e_max] and ϖ uniform on [0, 2π), equivalently (√e·cos ϖ,
    √e·sin ϖ) uniform on the disk of radius √e_max. The nested-sampling transform is
    e = e_max·u₁, ϖ = 2π·u₂;
  - perturber λ₀ uniform on [0, 2π);
  - Uranus (a, e·cos ϖ, e·sin ϖ, λ₀): uniform in a box centred on the Uranus-only fit, with
    half-widths of 50σ. Here σ is that fit's formal σ (finite-difference Jacobian) scaled by
    √(χ²/dof). The scaling is needed because the Uranus-only fit leaves the perturber's signal
    in its residuals (χ²/dof = 21.6 on seed 0). Its raw σ is ~4.6× too small, and the
    perturbed posteriors sit 40–63 raw σ from its centre, so an unscaled 50σ box would cut
    them. Scaled, every emcee sample of both seed-0 modes lies within 27% of the box half-width.
    Seed 0 half-widths: a 6.1e-3 AU, e·cos ϖ 4.3e-4, e·sin ϖ 4.1e-4, λ₀ 1.1e-3 rad (box ln
    volume −24.66). Rationale: the box only has to contain every posterior of interest. It is
    shared across all Phase 2 models fitted to the same data (hidden mass, modified law,
    mundane), so its log-volume largely cancels in Bayes factors. A model whose Uranus
    posterior came near the box edge would show up as samples piling at the edge, so that is
    checked.
- Derived quantity per sample: perturber geocentric ecliptic longitude on 1846-09-23
  (light-time corrected, Earth position from Horizons).
**A3 (acceptance, slow):**
  *Amended (M5a review):* the fixed targets ("free-a 1846 longitude within 1°", "Le Verrier mode
  within ~2°") are removed. The M5a profile likelihood showed that these data constrain the free-a
  1846 longitude only to about ±10° (Δχ² = 1). A fixed tolerance would therefore test the noise
  draw, not the code. The replacements:
  - **A3a (bias):** on noise-free data, the free-a fit recovers true Neptune's 1846-09-23
    geocentric longitude to within **0.5°**.
  - **A3b (calibration):** across **N ≥ 10 noise seeds**, the true 1846-09-23 longitude falls
    inside the 68% posterior interval (free a, e ∈ [0, 0.3]; **dynesty** posteriors, all modes)
    in **50–85%** of seeds. Report the interval width, and beside it the M5a profile's Δχ² ≤ 1
    interval for each seed.
  - Mass–distance posterior shows the expected ridge (plot saved).

## M6: Historical data
M5b is closed. Its methods (partitioned nested sampling, ERR_INFLATION ×1.8, free a, e ≤ 0.3 headline prior) carry over unchanged unless this section says otherwise.

### Goal

Run the validated pipeline on the observations Le Verrier actually used, and answer three questions with calibrated numbers:

1. Where does the real 1781–1845 Uranus data place the perturber, and does the 68% interval contain Neptune's true 1846 longitude (325.87°)?
2. How much did the outcome depend on the data choices of 1846: Bouvard's, Adams's and Le Verrier's?
3. Peirce's question (1847): was Galle's find a "happy accident"? Measure this as the probability that a given model points within Le Verrier's actual miss of the true position.

### Ground rules (in addition to the project rules)

- **Synthetic before historical, again.** Before fitting any real observation, repeat the A3b coverage test on synthetic data generated with the M6 observation model: real epochs, per-observer error structure, the same gaps. If coverage fails, stop and report. Don't tune.
- **Pre-register.** Write the A6 thresholds into this file and commit them before the first real-data fit. Never change them after seeing real results.
- **Keep the held-out data held out.** The Neptune precovery positions (M6e) are never used in fitting or model selection. Load them only in the A6e test.
- **Provenance for every number.** No observation enters the pipeline without a source volume, page and row reference.

### M6a: Transcription and provenance

**Stop condition first.** Before writing any loader, determine what Le Verrier's memoir actually tabulates and write it up in `docs/notes/M6.md`:
- observed positions (geocentric or heliocentric, RA/Dec or longitude/latitude), or
- residuals O − C against his own theory, or
- normal places (averages of several observations).

If it gives only residuals against his own tables, positions can't be recovered without those tables. Report this and stop; the source plan changes.

Sources, in priority order:
1. Le Verrier, "Recherches sur les mouvements de la planète Herschel (dite Uranus)", *Connaissance des Temps pour 1849* (Additions), 1846. Bertrand (2025, *Gazette des Mathématiciens* 186) counts 279 observational data points.
2. Le Verrier's notes in *Comptes Rendus* vols 21–23 (1845–46), on Gallica.
3. For the "what was true" variant: the original Greenwich and Paris meridian volumes, plus Airy's Greenwich planetary reductions.
4. The ancient sightings (1690–1771: Flamsteed, Bradley, Mayer, Lemonnier), as Le Verrier reduced them.

Process:
- Get the scans onto the local machine (the user can download them if your sandbox can't). Commit page images or stable URLs under `data/historical/sources/`.
- **Double-entry transcription.** Do two independent passes into `data/historical/leverrier_1846_pass{1,2}.csv`, then diff them. Every disagreement is resolved against the page image and logged. The user spot-checks 10% of rows chosen at random.
- Schema for `leverrier_1846.csv`, one row per observation or normal place:
  - `id`, `epoch` (as printed), `epoch_tt_jd`
  - `observer`, `observatory`, `instrument` (if known)
  - `quantity` (what was measured), `value`, `unit`, `frame_as_printed`
  - `lv_residual` (Le Verrier's O − C, if printed), `lv_weight` (if printed)
  - `era` (ancient / modern), `source_ref` (volume:page:row), `notes`

### M6b: Observation model

- **Frames.** Convert each row from its printed frame into the pipeline's frame.
  - Keep two reduction tracks:
    - **1846 track:** Le Verrier's constants for precession, aberration and nutation, as recorded in his memoir.
    - **Modern track:** IAU 2006 constants.
  - Store both reduced values. Report the difference per epoch as a figure.
- **Empirical error calibration.** Modern JPL ephemerides fix Uranus's true historical positions to far better than 1″. So for each observer and era, "observed minus JPL Uranus" *is* the real observation error.
  - Fit per-group offsets and scales. Groups are observer × instrument × era, falling back to coarser groups if counts are small.
  - This uses only Uranus truth, not Neptune, so it isn't leakage. State this explicitly in M6.md.
- **Likelihood.** The current likelihood plus per-group nuisance offsets, and error scales fixed from the calibration.
  - **Offset priors (pre-registered 2026-11-06, for A6c and the M6c headline):** b_g ~ N(0, τ²), τ = 2″ for modern groups and 5″ for the ancient group. The offsets are free nuisance parameters, marginalised, never fixed to O − JPL. A6c's synthetic truth draws the group offsets from these same priors.
  - **τ_ancient sensitivity (M6c, not A6c):** τ_ancient = 10″ and τ_ancient → ∞ (the latter should approach the `bouvard` variant). Report the posterior 1846 longitude and ln B against τ_ancient.
  - Run a sensitivity case with the scales free, using weak priors.
  - Robustness check: a Student-t likelihood (ν = 4) to catch outliers such as Flamsteed 1690.

### M6c: Fits under three data choices

Use the same sampler settings and partition as M5b, with free a and e ≤ 0.3.

| Variant | Data used | Historical basis |
|---|---|---|
| `bouvard` | modern observations only (1781–1845) | Bouvard 1821 set the ancient sightings aside |
| `adams` | all except Flamsteed 1690 | Adams treated 1690 as uncertain |
| `leverrier` | all | Le Verrier's choice |

Run each on both reduction tracks, which gives 6 fits. Also run the Le Verrier shape model (fixed a, as in M5) on the `leverrier` variant, both tracks.

Report for each fit:
- the posterior 1846 longitude: median and 68% interval
- the miss from 325.87°, in degrees and in units of the interval half-width (σ)
- mass and a posteriors
- ln Z, with the ×1.8 error inflation

### M6d: The Peirce number

For each fit, compute P(|λ − 325.87°| ≤ 0.87°), the probability the posterior puts on landing at least as close as Le Verrier's actual 52′ miss.
- Compare the Le Verrier shape model against free a.
- Express Peirce's "happy accident" as the ratio of these probabilities, together with the Bayes factor between the two models.

### M6e: Held-out test on Neptune precoveries

Use each M6c posterior to predict Neptune's geocentric position on:
- Lalande, 1795 May 8 and May 10
- Galileo, 1612 Dec 28 and 1613 Jan 28. The sources disagree on 27 vs 28 Jan; use Kowal & Drake 1980, *Nature* 287:311.

Record the predicted position and its 68% and 95% regions, and the observed positions with their stated uncertainties. Galileo's are large: his drawings were sketches, not measurements.

### Acceptance tests (pre-registered)

| Test | Pass condition |
|---|---|
| A6a transcription | Double-entry disagreements all resolved and logged; 10% spot-check shows zero uncorrected errors |
| A6b reduction | Modern-track residuals against JPL have median \|O−JPL\| within each group's published or expected precision. Report per group; no single pass/fail. Also the equinox/precession loader test below, which must pass |
| A6c synthetic coverage | On synthetic data with the M6 observation model, 10 seeds: truth within the 68% interval in 5–8 of 10 (same rule as A3b) |
| A6d leakage | Static check: no code path in M6a–M6d reads Neptune's real ephemeris or the M6e precovery file |
| A6e held-out | Reported, not pass/fail. Whether the 1795 positions fall inside the 95% predictive region. Interpret with Standish's caution that Lalande's observations are inconsistent with the others |

Headline result: the `leverrier` variant on the modern track, free a, e ≤ 0.3. Everything else is sensitivity.

**A6b loader test (fast, synthetic; carried over from the earlier M6 section).**
- Every observation records two things, with a citation:
  - its **equinox reference**: `of_date`, or `fixed` with the epoch (e.g. mean equinox of
    1800.0);
  - the **precession constant the source used**, e.g. Bessel's or Struve's value in ″/yr, plus
    any rate term.
- If a source pre-reduced its positions to a fixed epoch, first undo that reduction with *the
  source's own* constant, to recover the position of mean equinox of date. Then convert of date
  → J2000 with IAU 2006 (`frames.py`). Never apply IAU 2006 directly to a position another
  precession model produced.
- Observations missing either field fail to load.

The test:
1. Take known J2000 longitudes at 1690–1846 epochs, convert them to of date with IAU 2006, then
   pre-reduce them to a fixed epoch (1800.0) with a non-IAU constant, as a source would.
2. Use 50.10″/yr. Real 19th-century constants are within about 0.03″/yr of IAU 2006, too close
   to tell the paths apart reliably.
3. The loader must recover the J2000 longitudes to < 0.01″.
4. Two wrong paths must each be off by > 1″ at some epoch, so the test can't pass by accident:
   (a) treating the pre-reduced values as of date, and (b) undoing the reduction with IAU 2006
   instead of the source's constant.
5. A record missing equinox or constant metadata raises.

The earlier M6 numeric targets (A4: residual RMS ≲ 3″, 1846 longitude within ~2°) are dropped.

### Outputs

- `data/historical/leverrier_1846.csv`, with sources and a resolution log
- `site/data/m6_fits.json`: all variants, both tracks, the Peirce numbers, and the held-out predictions
- figures:
  - residuals against JPL by observer
  - posterior longitude for the three data choices
  - the Peirce ratio
  - the 1795 sky chart: predicted region, Lalande's positions, true Neptune
- `docs/notes/M6.md`: what the memoir tabulates (from the M6a stop condition), the transcription log, the error model, results, and limitations
- ledger entries as usual

### Out of scope for M6

- Re-deriving Le Verrier's own perturbation theory or his tables
- Adams's full dataset (a possible M6+ follow-up)
- Any change to the M5b sampler, partition or error inflation

## M7: Export and site
- `export/to_json.py` writes `site/data/`:
  - `sightlines.json`: [{jd, earth:[x,y], uranus:[x,y]}] 1690–1846
  - `residuals.json`: {jd[], obs_minus_no_perturber[], obs_minus_fitted[], sigma[]}
  - `skymap_1846.json`: 2D histogram (lon, lat) of posterior + markers {le_verrier, adams, galle_found, true}
  - `degeneracy.json`: thinned (log m, a) samples
  - `energy.json`: {t_yr[], rk4_dE[], wh_dE[]}
- `site/`: five panels, in order: sightline animation → residual curve (toggle with/without
  perturber) → integrator energy panel → degeneracy scatter → 1846 sky map. Each panel has one
  plain-language caption sentence.
**A5 (acceptance):** Pages deploy succeeds from CI; all five panels render with no console errors;
total `site/data` under 5 MB.

---

## Out of scope for Phase 1
Inclination of the perturber, fitting Jupiter/Saturn masses, raw geocentric RA/Dec forward model,
Mercury, model comparison, any LLM usage. These belong to Phase 2/3.

## Numbers to record in `docs/notes/` as you go
Integrator agreement vs REBOUND (″), WH runtime per 65-yr run, no-perturber peak-to-peak (″),
fitted RMS (″), MCMC τ and wall time, posterior median and 68% interval for m, a, and the
1846-09-23 longitude, versus Le Verrier's and Adams's published predictions.
