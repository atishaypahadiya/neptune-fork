# Method

This is the short version of the milestone notes (`docs/notes/M0.md` to `M6.md`). It covers
the main study up to M5b, plus the M6 data and observation model. It reports no real-data results:
the M6c headline fit on Table 77 is in `docs/notes/M6.md`. Every number here is from synthetic
data or from a validation against JPL Horizons or REBOUND. Each section names the note it condenses.

Units are AU, days and solar masses, with G = k² and k = 0.01720209895 (Gaussian). Angles are
radians inside the code and arcseconds only at input, output and in plots.

## The problem

Uranus drifted about 2′ off any Keplerian orbit that Jupiter and Saturn could explain. The
reason was Neptune. I fit Uranus's observed positions with Neptune left out of the model and
an unseen perturber put in, then ask where that perturber is on **1846-09-23**, the night Galle
found Neptune. The answer is a posterior distribution, not a single best fit. The data allow a
lighter planet closer in or a heavier one farther out (the mass–distance degeneracy), and the
posterior has to show that spread.

## Frames (M1)

- **`ecl_j2000`** is JPL Horizons' ecliptic, which is astropy's `HeliocentricEclipticIAU76`:
  ICRF rotated by ε₀ = 84381.448″. It matches Horizons to 1.5e-11″. Astropy's
  `HeliocentricMeanEcliptic(J2000)` is a different frame (0.04″ off), and a test checks that
  the two frames differ.
- **`ecl_of_date`** uses IAU 2006 precession through astropy. Precession is about 50″/yr, much
  larger than the ~130″ signal, so a frame bug could fake a planet. Every angle variable
  carries its frame in its name. A linear 50.29″/yr rate would already be 11″ wrong at 1690,
  so the full model is used everywhere.
- **Wrap-safe residuals.** Every angular residual is
  `((obs − model + π) mod 2π) − π`.
- **Central body.** The inner planets' masses are added to the Sun's. The central body then has
  to sit at the barycentre of the Sun and inner planets, not at the Sun's centre. Starting from
  Sun-centred states left a 7.4 cm/s velocity error. That error grew into a ±35″ drift in
  Uranus over 1690–1846. With the starting states moved to that barycentre, the drift is
  0.09″ (M4).

## Integrator (M2, M3)

`src/fork/integrators/` holds the integrators, written for this study and compiled with numba. REBOUND is used only
in tests.

- **Kepler step.** A universal-variable f–g propagator with Stumpff functions. It uses Newton
  iteration inside a guaranteed bracket and falls back to bisection. An elliptic step is first
  reduced modulo one period. Over 1000 random cases it agrees with an analytic ellipse to
  1.2e-13·a.
- **Wisdom–Holman** in democratic heliocentric coordinates. The Hamiltonian splits as

  H = H_Kepler + H_Sun + H_interaction,

  and one step of length h is Kepler(h/2) → Sun(h/2) → kick(h) → Sun(h/2) → Kepler(h/2). This
  is the same ordering as REBOUND's WHFast-DH. The default step is dt = 12 days.
- **Validation (A1)** uses Sun, Jupiter, Saturn and Uranus from Horizons at 1781-03-13.

  | check | target | result |
  |---|---|---|
  | Uranus longitude vs REBOUND WHFast-DH, 1781–1846 | < 0.01″ | **1.65e-9″** |
  | max \|ΔE/E\| over 1000 yr, checked at every step | < 1e-8 | 7.3e-9 |
  | one 65-yr run after compilation | < 10 ms | ~1.9 ms |
  | error ratio per halving of dt (24 → 3 d) | 3.5–4.5 | 4.00, 4.00, 4.00 |

  RK4 over the same 1000 years shows secular energy drift (`docs/notes/M3_energy.png`).

## Forward model (M4)

The parameter vector θ holds Uranus's osculating elements at 1781-03-13 plus the perturber's.
Jupiter and Saturn start from Horizons and are not fitted. The model integrates the system and
returns Uranus's heliocentric longitude (M5), or its full state for the apparent place (M6).

- **Uranus's a, e, ϖ, λ₀ are always free.** Its i and Ω are held at their 1781 values. From
  longitudes alone, i and Ω only change the reduction to the ecliptic, a ~9″ second harmonic.
  If they are left free, they absorb the perturber's signal (χ²/dof drops from 21.6 to 6.5
  with no perturber at all). Holding them costs very little: Neptune moves Uranus's i and Ω by
  amounts worth only **0.017″** in longitude over 1781–1846. A test pins that number.
- **Perturber:** log₁₀ m [M⊕], a, e, ϖ, λ₀. It is planar for the M5 longitude fits. For M6,
  where declinations enter, it gains an inclination i and node Ω relative to Uranus's plane
  (below).
- **A2/A2b:** a Uranus-only fit to synthetic data leaves a residual of 143″ peak to peak, with
  the slope turning in 1828. With the true Neptune planted and Uranus refitted, the post-1781
  RMS is 0.93–1.02σ, which is pure noise.

## Observation models

### Synthetic longitudes (M4, used through M5b)

The truth is Horizons' Uranus, which includes the real Neptune. Gaussian noise is added: 22
pre-discovery sightings between 1690 and 1771 at σ = 12″, and one normal place a year for
1781–1845 at σ = 3″. The likelihood is

  −2 ln L = χ² = Σᵢ (wrap(λᵢ,obs − λᵢ,model) / σᵢ)².

### Le Verrier's observations, apparent RA and Dec (M6a, M6b)

- **Data.** Le Verrier's 1846 memoir, Table n° 77: 279 observations (19 ancient, 260 modern),
  as apparent geocentric RA and Dec of date. 29 rows give RA only, so there are 529
  measurements. The table was transcribed twice, from two different scans. The 14 cells where
  the passes disagreed were resolved against a third witness (the Gallica scan of the Connaissance des Temps volume), and each resolution is logged.
  Source misprints are kept as printed in the transcription and corrected in the loader from a
  separate misprint table, each with a reason (`data/historical/`).
- **Times.** Every row is in Paris mean time. This was established by a transit check: for a
  meridian observation the local sidereal time must equal the RA. Greenwich rows miss by
  9 m 21 s under Greenwich time and fit to about 1 s under Paris time.
- **Model chain.** The integrated Uranus state goes through:
  1. light time, r − τv;
  2. annual aberration, from Earth's Horizons velocity;
  3. ecl_j2000 → ICRS → true equator and equinox of date (`erfa.pnm06a`, IAU 2006/2000A).

  Run with the true 1781 states and a 3-D Neptune, the chain matches Horizons to **< 0.05″**
  in every observer group.
- **Inclined perturber.** A planar Neptune misfits Dec by 1.9″ rms in 1835–1845. The perturber
  therefore gets (sin i cos Ω, sin i sin Ω) relative to Uranus's plane, with i ≤ 10° and a
  prior uniform in cos i.
- **Errors.** Residuals are on the sky: Δα·cos δ and Δδ. Measurements are grouped by observer
  series × quantity (ancient; Greenwich; Paris 1801–27; Paris 1835–45). Each group g has an
  error σ_g and an unknown offset b_g ~ N(0, τ_g²). The offset is integrated out analytically:

  −2 ln L = Σ_g [ Σ_{i∈g} rᵢ²/σ_g² − (Σ_{i∈g} rᵢ/σ_g²)² / (n_g/σ_g² + 1/τ_g²) ] + const.

- **Calibration.** σ_g comes from observed minus the real Uranus (Horizons), not from the fit.
  This uses no Neptune information, so nothing leaks. Values range from 1.2″ (Paris 1835–45
  Dec) to about 5″ (ancient). τ was fixed before any fit: 2″ for the modern groups and 5″ for
  the ancient group.

## Priors (M5)

The model is sampled in the vector

  q = (Uranus a, e cos ϖ, e sin ϖ, λ₀; perturber log₁₀ m, ln a, √e cos ϖ, √e sin ϖ, λ₀).

The prior density is flat in q on this support:

- log₁₀ m ∈ [0, 2], i.e. 1–100 M⊕;
- a log-uniform on [20, 60] AU;
- e uniform on [0, e_max], with e_max = 0.3 for the headline results and 0.1 as a sensitivity
  case;
- ϖ and λ₀ uniform. A flat prior on the (√e cos ϖ, √e sin ϖ) disk is exactly uniform in e and
  ϖ, and it stays well defined as e → 0;
- Uranus inside a box of ±50 scaled σ around the Uranus-only fit. The σ are scaled by
  √(χ²/dof), because that fit still contains the perturber's signal. The box's log-volume is
  the same for every model fitted to the same data, so it cancels in Bayes factors.

One definition (`inference/priors.py`) serves emcee, dynesty and the profile fits.

## Partitioned nested sampling (M5b)

Simpler samplers failed in two ways:

- **emcee** stayed in whichever mode it started in. Two seed-0 modes, a ≈ 30 AU and
  a ≈ 22 AU with e ≈ 0.28, never mixed in 20–40k steps, even though each chain passed the 50τ
  test.
- **Unpartitioned dynesty** got a reproducible ln Z but split the posterior mass between thin
  pockets at random. One region held 1 % of the mass in one run and 39 % in another.

The fix is to split the prior into fixed **cells** and sample each cell separately:

| cell | a [AU] | e | prior weight w_k |
|---|---|---|---|
| 0 | [20, 25] | all | 0.203 |
| 1 | [25, 28] | all | 0.103 |
| 2 | [28, 60] | [0, 0.12) | 0.277 |
| 3 | [28, 60] | [0.12, 0.3] | 0.416 |

w_k = ln(a_hi/a_lo)/ln 3 × (e_hi − e_lo)/0.3. Each cell is an ordinary restricted prior, so the
cell results combine exactly:

  Z = Σ_k w_k Z_k,  s_k = w_k Z_k / Z,  σ(ln Z)² = Σ_k s_k² σ_k².

The posterior is each cell's weighted samples scaled by its share s_k. The cell edges affect
efficiency, not correctness: any partition gives the same answer in the limit. The edges were
chosen so that no single run has to find a thin pocket next to a large mode. They were fixed
in advance and are the same for every dataset.

Settings: `bound='multi'`, `sample='rslice'`, stop at Δln Z < 0.1, perturber ϖ and λ₀
periodic, nlive 1000, and nlive 2000 in cell 0, whose thin pocket was under-sampled at 1000.
Every long run writes a progress log and checkpoints so it can resume.

## Evidence errors: the ×1.8 inflation

dynesty's quoted ln Z error understates the run-to-run scatter on these posteriors. On seed 0,
each cell was run with different sampler seeds (two seeds, three for cell 0). The pooled
single-run standard deviation of ln Z_k was **0.40**, against a quoted **0.22**: a factor of
1.8. With 5 degrees of freedom, that factor's own 68 % range is about 1.5–2.6. So every ln Z
error reported or compared is

  σ_eff = 1.8 × σ_quoted  (`nested.ERR_INFLATION`),

and cell-share errors scale with it. Later model comparison uses the same factor.

## Calibration tests

A fixed tolerance on the 1846 longitude would test the noise draw, not the code. The M5a
profile likelihood shows these data pin that longitude only to about ±10° (Δχ² = 1). So the
acceptance tests check bias and coverage instead.

- **A3a, bias.** On noise-free data, the free-a *point fit* (M5a, least squares) recovers true
  Neptune's 1846 longitude to **0.07°**. The target was within 0.5°. This is a check on the
  model and the fitter, not on the posterior. For the noise-free posterior, see below.
- **A3b, coverage.** With free a, e ≤ 0.3 and 10 noise seeds, the true longitude (325.87°) lies
  inside the posterior's central 68 % interval for **8 of 10 seeds**. A calibrated 68 %
  interval gives a two-sided binomial p = 0.52, and the pass range was 5–8. The intervals are
  7.2–10.7° wide (median 9.1°). The largest miss, seed 5 at −1.74σ, comes from prior volume
  inside a flat likelihood band, not from a missed mode (`docs/notes/M5.md`).
- **Prior sensitivity.** With e ≤ 0.1 the intervals are a third narrower (median 6.1°), and
  coverage drops to 5/10 (p = 0.31). ln Z prefers e ≤ 0.1 on all 10 seeds, by
  Δln Z = +0.2 to +1.5. That is the Occam factor of the smaller prior: the data don't need
  e > 0.1.
- **Noise-free posterior (M5b).** The median 1846 longitude is 322.05°, against the true
  325.87°: an offset of −3.8°, or −0.75σ, with the truth inside the 68 % interval. The offset
  lies along the mass–distance direction (cosine 0.99), where the likelihood is flat. So the
  posterior settles where more of the prior volume fits well. That is a lean, not a likelihood
  bias.
- **Le Verrier's published orbit** (a = 36.15 AU, e = 0.108, ϖ held) is disfavoured against free
  a by ln B = −11.9 to −21.6 on seeds 0–4.
- **A6c, coverage with the M6 observation model (passed).** This repeats A3b on synthetic
  data built from Table 77's real epochs. The truth is the real Uranus (Horizons, so the real
  3-D Neptune), converted to apparent RA/Dec. Group offsets are drawn from the pre-registered
  τ, and the calibrated σ_g noise is added. The fit uses the inclined perturber and the same
  four cells and settings. The pass condition is the same 5–8 of 10, and it was written into
  the spec before the run.
  - **Result:** the truth is inside the 68 % interval for 8 of 10 seeds (two-sided binomial
    p = 0.52). The rms offset is 0.87σ and the median 68 % width is 5.1°.
  - The intervals are narrower than M5b's 9.1° mostly because Table 77 has more measurements
    and smaller calibrated σ than M5b assumed (`docs/notes/M6.md`).
  - Only after this pass may real observations be fitted.

## Not modelled, and known limits

- Light deflection (milliarcseconds), diurnal parallax (the printed places are geocentric), and
  the Horizons/ICRS frame tie (~0.02″).
- Jupiter's and Saturn's masses and orbits are fixed at modern values.
- The Kepler solver's full-period closure misses 1e-12 AU for e ≳ 0.6 at large a. The cause is
  double-precision conditioning, not the algorithm, and the prior's e ≤ 0.3 is well inside the
  range where the target is met (M2).
- The table was transcribed twice from the scan and the two transcriptions reconciled, each
  disagreement resolved against the page image. 28 rows were then checked by hand against the
  scan, with no discrepancies; the remaining 251 rows have had no human check
  (`docs/notes/M6.md`).
