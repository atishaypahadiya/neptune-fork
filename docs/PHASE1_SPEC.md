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
  *Amended (M5a review):* point fits also enforce the prior bound e ≤ 0.3.
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
    inside the 68% posterior interval (free a, e ∈ [0, 0.3]) in **50–85%** of seeds. Report
    the interval width.
  - Mass–distance posterior shows the expected ridge (plot saved).

## M6: Historical data
- `data/observations.py`: load historical Uranus heliocentric longitudes (Whitmore's set if license
  permits, else transcribe from Le Verrier's memoirs / Bouvard's tables). Convert to J2000 via
  `frames.py`. Store with per-point σ and source column.
- Every observation records its **equinox reference** — `of_date` or `fixed` with the epoch
  (e.g. mean equinox of 1800.0) — and the **precession constant the source used** (e.g. Bessel's
  or Struve's value, in ″/yr, plus any rate term), with a citation. If a source pre-reduced its
  longitudes to a fixed epoch, first undo that reduction with *the source's own* constant to
  recover longitude of mean equinox of date, then convert of-date → J2000 with IAU 2006
  (`frames.py`). Never apply IAU 2006 directly to a longitude another precession model produced.
  Observations missing either field fail to load.
**Test (fast, synthetic):** take known J2000 longitudes at 1690–1846 epochs, convert to of-date
with IAU 2006, then pre-reduce to a fixed epoch (1800.0) with a non-IAU constant, as a source would.
Use 50.10″/yr: real 19th-century constants are within ~0.03″/yr of IAU 2006, too close to tell the
paths apart reliably. The loader must recover the J2000 longitudes to < 0.01″. Two wrong paths
must each be off by > 1″ at some epoch, so the test can't pass by accident: (a) treating the
pre-reduced values as of-date, and (b) undoing the reduction with IAU 2006 instead of the source's
constant. Also: a record missing equinox or constant metadata raises.
- Re-run M4 and M5 on historical data.
**A4 (acceptance, slow):** residual RMS with fitted perturber **≲ 3″**; posterior median
1846-09-23 longitude within **~2°** of true Neptune. Document outliers dropped and why.

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
