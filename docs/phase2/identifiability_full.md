# Mercury: the identifiability study with the full forward model (pre-registration)

**Pre-registration, 15 November 2026.** Written before any code of this plan exists and before
any data set of it is drawn. Everything above the line at the end is fixed by the commit that
adds this file; amendments and results go below it, dated, never edits to the text above.

It replaces the method of `identifiability_study.md`, which stopped at its gate V2: the six
elements of Mercury, marginalised linearly about the true ones, did not reproduce the full
integration at the posterior's own draws in any of its 12 fits. **The question, the truths, the
models, the priors on the models' parameters, σ₁'s calibration rule and the readings are that
plan's** (its §§ 1, 4, 5 and 7, and the rule of § 3). Three things change: the observation is
the contact time itself rather than a longitude standing in for it; all six elements are fitted
in every likelihood call, with the full integration; and each evidence is integrated on a grid.

**No observation of Mercury is fitted here or anywhere.** Newcomb's contact table gives the
epochs, the kind of contact and the weights only.

## 0 Why the inclination and the node are free, and why that changes the observation

A transit is seen only when Mercury is near a node of its orbit, and the times of first to
fourth contact depend on how far from the Sun's centre it crosses: contact times constrain the
inclination and the node strongly. The node's motion is also where the models differ most: a
solar J2, a ring in the ecliptic and a planet in the ecliptic move the node of Mercury's orbit,
while the 1PN term does not. Holding the inclination and the node, as an earlier draft of this
plan did, would have removed exactly that information. With them free, the earlier plan's
stand-in observation, Mercury's heliocentric longitude at each epoch, cannot carry it: the
longitude does not see the latitude at transit. I measured this (§ 1): with longitudes, six
free elements and a J2 model fit 1PN data almost exactly; with contact times they do not.

## 1 What I measured before writing this

Timing and conditioning only, on data that no reading uses: noise-free data made by a 1PN
term at strength 0.7 (not a truth of this study), fitted by the J2 model at J2 = 0.7 × truth
d's value (not a truth either), at the 58 printed epochs, all six elements free, Levenberg–
Marquardt with every residual and Jacobian column from a full integration. No data set of the
study was drawn, no χ² of the study was computed, and nothing below was set from these runs
except the choice of the observation.

**Disclosure.** This test configuration, 1PN at strength 0.7 fitted by J2, was run and its
results seen before this plan was frozen: χ²_min 24.8 with contact times against 0.001 with
longitudes. It lies outside the study's truths and data sets, and it was used for one decision
only, the choice of the observation (§ 0, § 3). No grid, gate, tolerance or reading below was
set from it, and its χ² values are not part of any result.

| | longitude (earlier plan's observation) | contact time (this plan) |
|---|---|---|
| condition number of the Jacobian at the optimum | 6.6 × 10⁷ | 4.5 × 10⁵ |
| χ²_min of J2 against 1PN-made data (each observation at its own σ₁, § 3) | 0.001 | 24.8 |
| iterations from the data's own elements | 4 | 4 |
| full integrations per call, from cold | 29 | 29 |
| ln Laplace − ln importance sampling at the optimum (six elements, 400 draws) | −0.67 ± 0.22 | −0.011 ± 0.037 |

- One full integration (Sun and eight planets, 1677 to 1881 from 1900.0, the validated step):
  0.54 CPU-s (Newtonian) to 0.80 CPU-s (ring). One cold call: 17 to 20 CPU-s.
- **Warm starts.** Started from the neighbouring grid point's optimum with fresh Jacobians at
  every iteration, a call still took 4 iterations and about 30 integrations: the conditional
  optimum moves at every strength step. Reusing the previous point's Jacobian, updated by
  Broyden's rank-one rule, six successive strength steps took 9, 27, 9, 9, 18 and 9
  integrations (mean 13.5, about 9 CPU-s), and reached the same χ² as the fresh-Jacobian fits
  to 10⁻³ at every step.
- The Laplace covariance has a condition number near 10¹¹: it is formed from the singular
  values of the whitened Jacobian, never by inverting JᵀWJ, and the importance-sampling
  proposal is drawn in whitened coordinates. (A library routine that drops eigenvalues below
  a relative 10⁻¹⁰ gave a spurious 22 in ln before this was done.)
- One printed epoch is 18.5 h from any contact of its kind: table III row 12, the exterior
  contact of 1848 dated November 8, a logged misprint (the table dates the transit's contact
  III November 9; `transcription.md`). Every other row is within 30 min of its contact
  (median 3 min).

## 2 Inputs

| input | SHA-256 |
|---|---|
| `data/historical/mercury/newcomb_1882_contacts.csv` (epochs, contact, weights only) | 70d4d0d1ea242c55c227d593b1cf9238593ce52c3cd39ae0ce9c1fdbd41712f9 |
| `data/historical/mercury/newcomb_1882_misprints.csv` | 79c4720ccae315cc4948c0a4d08c6928ec521ff95dc546fa939c9d6a9388fbce |
| `tests/fixtures/phase2/v4_reference.npz` (planetary states at 1900.0) | f8032fc2f622112303e9e3912fc32443b54e8e206251bb91babd7c3d29f24f77 |
| `src/fork/model/inner_system.py` | a132da7aab691230d988ee3c2d42b2ec16655af7e65159856ff480425e6f512b |
| `src/fork/phase2/secular.py` | 55d0eb506eaff1899731d119324474287bcb6091598e6dcd72e98213afab880b |
| `src/fork/phase2/identifiability.py` | 14676bb8d02a945dae5ed66f8c75c13084857b4e899357cc0b5faf9fa3c4d3f0 |

**The truths are the earlier plan's, already tuned** (its V1 passed), copied into a tracked file
`docs/phase2/identifiability_full_truths.json` before any data set is drawn:

| truth | cause | value | integrated rate (″ a century) |
|---|---|---|---|
| a | the 1PN term | strength 1 | 42.989 |
| b | ring in the ecliptic at 0.2 AU | 1.9811152404862693 × 10⁻⁷ M☉ | 42.980 |
| c | planet on a circle in the ecliptic at 0.2 AU, longitude 0 at epoch | 1.980728405006917 × 10⁻⁷ M☉ | 42.980 |
| d | solar J2 | 3.399891252673649 × 10⁻⁴ | 42.980 |

## 3 The observation: the geocentric contact time

For each row of the contact table, the synthetic datum is **the time of that contact as seen
from the Earth's centre**, computed from the full integration:

- the contact is the instant when the apparent separation of Mercury's centre from the Sun's
  centre, seen from the geocentre, equals R☉ − r☿ for an interior contact (II, III) and
  R☉ + r☿ for an exterior one (tables III and IV, whose rows have no contact number); R☉ =
  959.63″ and r☿ = 3.36″ at 1 AU, scaled by the distances;
- it is found from the integrated states of Mercury and the Earth at the printed epoch by
  Newton's method on the offset, the two bodies carried over the offset (under 30 min) by
  their heliocentric Keplerian acceleration; the root is required to 10⁻¹² rad;
- the residual is model contact time minus datum, in seconds, with error σ_i = σ₁ / √w_i (w_i
  Newcomb's weight). There is no ρ_i: the datum is already a time.

Not modelled, the same in data and model, so without effect on synthetic data: light time,
aberration, ΔT, the observers' places (the table's times are already reduced to the
geocentre), the limb's appearance. The semidiameters' values do not enter any reading for the
same reason.

**The row that is not a contact.** A row whose printed epoch lies more than 1 h from the
model's contact of its kind (on truth a's noise-free model) is left out, and listed. Fixed now,
this leaves out table III row 12 alone (§ 1), weight 0.2: 57 contacts remain. The misprint is
not corrected, as the transcription's rule says; the row is set aside because its printed time
is not the time of a contact, which is what this observation needs.

**σ₁** follows the earlier plan's rule unchanged: the timing error of a unit-weight contact at
which this observation, with the six elements free and the 1PN strength free, measures that
strength to 3.0 % on truth a's noise-free data (from the curvature of χ² by full
integrations). The measurement of § 1 gave about 17 s (the earlier observation's value was
13.4 s); the value is fixed in Stage 1, before any noisy data set is drawn. The noise of each
data set uses the earlier plan's rng, `[2026, 21, 1882, truth, seed]`.

## 4 The likelihood and the evidence

### 4.1 Six elements fitted in every likelihood call

- **Fitted:** a, e, i, Ω, ω and M of Mercury at 1900.0, by Levenberg–Marquardt with **every
  residual from a full integration**. A call is warm-started from the nearest grid point
  already done, with that point's Jacobian updated by Broyden's rank-one rule after each step;
  a fresh Jacobian (forward differences, the earlier plan's steps, each column a full
  integration) replaces it when a step fails twice, and at the first point of a scan. A call
  has converged when χ² drops by less than 10⁻³ and the Gauss–Newton step from a fresh Jacobian
  at the end predicts a further drop below 10⁻³ (at most 20 iterations); that last Jacobian
  is the one Laplace's method uses. Each model therefore has up to eight free
  parameters: the six elements, fitted in every call, and its own one or two (strength; radius
  or longitude), on the grid of § 4.2.
- **Marginalised by Laplace's method at that optimum:** ln L(p) = −½ χ²_min(p) − Σ ln s_k
  + a constant, s_k the singular values of the whitened Jacobian, with flat priors on the six
  elements over a box that is the same for every model, so the constant cancels in every Bayes
  factor. Nothing is linearised about the true elements; the one approximation left, that the
  likelihood is Gaussian in the six elements near their conditional optimum, is tested (V2′).

### 4.2 The evidence: quadrature on a fixed grid, no sampler

- **Strength** (log₁₀ J2 or log₁₀ mass): a fine scan of 15 points over ±6 posterior widths
  about the best value, χ² interpolated by a cubic spline and the integral of exp(−χ²/2) taken
  on a sub-grid 20 times finer. The best value is found by a coarse scan at 0.1 decade over the
  widest prior of R3 (8 decades about its centre) at the first radius or longitude of a scan,
  and from the neighbouring one elsewhere. Every strength scan also evaluates 6 points at ±1, ±2
  and ±4 decades from the best value, to show that the integrand outside the fine scan is below
  e⁻²⁵ of the peak; where it is not, the fine scan is widened until it is.
- **Ring radius:** the earlier plan's 11 radii, 0.10 to 0.30 AU, a strength scan at each.
- **Planet:** 48 longitudes at epoch at each of the three distances, a strength scan at each;
  the three distances averaged with equal weight, as before.
- **N and G** have no parameter beyond the elements: one call each.

Error levels: noise-free data sets at different levels differ only in σ₁, which scales every
σ_i alike, so one scan serves all levels. Noisy data sets are scanned at their own level.

## 5 Validation, before anything is read

| gate | pass condition |
|---|---|
| V0 contacts | on truth a's noise-free model, every kept row's contact is found to 10⁻¹² rad and lies within 30 min of its printed epoch; exactly the rows of § 3 are left out |
| V1 rates | as before, on the tracked truths file: each truth's rate within 0.05″ of 42.98″ a century |
| V2′ Laplace | for each of the 12 noise-free fits of the models with parameters (J, R, P; four truths), at 10 points of its scan (the best point and points at Δχ² = 1, 4, 9 and 25 on either side along the strength, at the best radius or longitude where there is one), the Laplace value agrees within 0.1 in ln with an importance-sampling marginal over the six elements (800 full integrations from a Student-t with 5 degrees of freedom, centred at the optimum with 1.5 × the Laplace covariance, drawn in whitened coordinates; its own error below 0.03), at 95 % of the points |
| V3′ grid | on each truth's own model and noise-free data: halving every grid spacing (strength, radius, longitude; for the planet, at its true distance) changes ln Z by less than 0.05 at each of the three levels |
| V4 closure | the truth's own model on its noise-free data: minimum χ² below 0.5 |
| V5 no observed value | a test shows the synthetic data depend only on the epochs, the contact kinds, the weights, the truth and the seed |
| V6 the fits converge | every likelihood call converges within 20 iterations, except at most 0.1 % of calls, none of them within Δχ² = 25 of a scan's best point |

V2′ is the earlier plan's gate on Laplace's method, kept with two changes forced by six
dimensions: 800 draws instead of 400 (at 400 the importance-sampling error was 0.037 in § 1,
above the gate's 0.03), and whitened coordinates (§ 1). The tolerance, the points and the
95 % stay.

If a gate fails, the study stops there; the failure is the result.

## 6 Stages and readings

The readings R1 to R4 are those of the earlier plan § 7, word for word, with its grades. Each is
read when the data it needs exist:

| stage | data sets | readings |
|---|---|---|
| 1 | the four noise-free sets, all five models, all three levels (one scan each); σ₁ calibrated; gates V0 to V6 | **R3** complete; **R2** on noise-free data; **R1 and R4 on noise-free data only**, reported as such and not as R1 or R4 |
| 2 | seeds 1 to 8 at the calibrated σ₁ (32 sets) | **R1** complete; R2 on all sets |

**Stage 3** (seeds 1 to 8 at σ₁/2 and 2σ₁, 64 sets, which R4 needs) **is not in this plan.**
Whether to run it is decided after Stage 2, in a plan of its own; until then R4 is reported on
noise-free data only.

**Continuation, fixed now so that no result decides it:** Stage 2 runs after Stage 1 if
Stage 1's gates pass, whatever its readings.

### Cost

From § 1. A likelihood call costs about 9 CPU-s warm-started and 18 CPU-s cold (the first point
of a scan); one integration 0.54 to 0.80 CPU-s.

| model | calls per fit | CPU per fit |
|---|---|---|
| N, G | 1 (cold) | 18 s |
| J | about 50 | about 8 min |
| R | about 550 (11 radii) | about 1.4 h |
| P | about 3,600 (3 distances × 48 longitudes) | about 9 h |
| one data set, all five | about 4,200 | about 10.5 h |

| stage | work | CPU (estimate) | on 12 workers |
|---|---|---|---|
| 1 | 4 noise-free scans, about 42; V2′ (120 points × 800 integrations), about 17; V3′ (the truths' own models at half spacing, about 5,000 calls), about 13; σ₁ and V0, V1, V4 under 1 | about 75 CPU-h | about 6.5 h |
| 2 | 32 noisy scans | about 340 CPU-h | about 28 h |

The earlier draft, with four elements, estimated about 50 and 220 CPU-h. The planet model is
86 % of every scan. The earlier plan's grid (48 longitudes, 3 distances) is kept, so that the
readings stay that plan's.

## 7 What this can show that the earlier plan could not

Nothing new about Mercury: the question is the same. The differences: the observation carries
the latitude at transit, where the models differ (the earlier plan's § 8 named this as its main
limit); no likelihood rests on a linearisation about the true elements; no evidence rests on a
sampler's noise.

## 8 What this cannot show

The earlier plan's § 8 still holds (nothing about the real Mercury; no time-scale or observer
errors, no outliers; bodies in the ecliptic; σ₁ from Newcomb's summary), except its limit of
the longitude-only observation, which this plan removes. In addition:

- **The contact model is geometric and geocentric.** It omits light time, aberration, ΔT and
  the observers' places. For synthetic data this is acceptable: data and model share the same
  omissions, so they cancel. For any version with real data, light time and ΔT are the first
  required additions: ΔT over 1677 to 1800 is about 9 to 18 s by the usual reconstructions,
  the same order as σ₁ (about 17 s), and it changes across the span, so it cannot be absorbed
  as a common offset of the time scale.
- One printed row (§ 3) is left out.
- Laplace's method in the six elements is tested at 10 points per fit, not everywhere.
- The planet's distance stays on three values.
- R4 is not complete: its noisy part is Stage 3.

## 9 Layout

- Code `src/fork/phase2/identifiability_full.py` (the earlier study's code and its pre-registered
  inputs unchanged), tests `tests/test_phase2_identifiability_full.py`, cache
  `data/cache/phase2/identifiability_full/`, results `docs/phase2/identifiability_full.json`.
- One process per scan, 12 at a time, each resumable (one file per finished grid point), every
  step writing a ledger row.

---

## Amendment 1, 15 November 2026, before any Stage 1 scan

Written after the runner's smoke test and before any scan in the study's cache. Items 1 to 7
and 10 record choices the plan left open. Items 8 and 9 change text above the line (the Jacobian
of § 4.1 and gate V6); both were approved before this amendment was committed. Nothing in
this amendment was set from a Bayes factor or a reading.

**What the smoke test computed.** In a scratch cache, on the study's own noise-free data
sets: truth d fitted by J at nine strengths (once with each Jacobian of item 8), truth a
fitted by N and G (with each Jacobian), and one complete J scan on truth d with the N and G
calls on truth a. None of it is reused; Stage 1 computes everything again in the study's
cache. That scan gave J's ln Z on truth d's noise-free data, which readings R1 and R3 use; I
saw it, and nothing here was set from it.

1. **The Earth is the Earth–Moon barycentre.** The integration carries the barycentre, not
   the geocentre; data and model share it, so it cancels for synthetic data. A version with
   real data needs the geocentre's offset from the barycentre (about 4,700 km, which moves
   Mercury against the Sun by several arcseconds, of order a minute in a contact time),
   together with light time and ΔT (§ 8).
2. **V2′ has 9 points per fit, not 10.** The gate names the best point and Δχ² = 1, 4, 9 and
   25 on either side: 9 points. That is 108 points in all, so the gate's 95 % allows at most
   5 outside the tolerance.
3. **A call with no contact has likelihood zero.** Where the model gives no contact within a
   day of a printed epoch, or the root is not found to 10⁻¹² rad, the call's likelihood is
   zero. Such calls are counted and listed separately; they are not failures to converge
   under V6.
4. **What is splined.** The fine scan splines ln L itself (−½ χ²/L² minus the Laplace term),
   not χ² alone, since the Laplace term changes along the scan too. Outside the fine points,
   ln L is interpolated linearly between the coarse and tail points on the same sub-grid
   step. The integral is clipped to the prior box.
5. **The planet's scans run one process per distance**, its 48 longitudes in sequence. V3′'s
   halved grid (96 longitudes at the true distance) runs as 4 processes of 24 longitudes, and
   the first longitude of each starts from its own coarse scan.
6. **Refining and widening the fine scan.** The 15 points are recentred and their spacing
   reset (at most six times) while the spline's best value moves by more than one width or
   its width changes by more than a factor 0.7 to 1.5. The points are then extended at the
   same spacing until both ends are at least Δχ² = 50 L² above the minimum for the largest
   level the scan serves (L = 2 for noise-free data), so that ln L there is at least 25
   below the peak. Where a tail point is not below e⁻²⁵ of the peak, the fine points are
   extended toward it (at most 1,500 more).
7. **The coarse scan** runs upward from the low end of its 8 decades. Its first point starts
   from the base elements; each later point starts from its neighbour.

8. **Fresh Jacobians along the singular directions (changes § 4.1).** The first Jacobian of
   a scan stays the plan's: forward differences with the earlier plan's element steps. Every
   later fresh Jacobian is taken by forward differences along the singular directions of the
   Jacobian in hand, each step one whitened unit (θ + v_k/s_k). This covers the one after
   two failed steps, the one for the convergence test and the one Laplace's method uses. It
   still costs six integrations.

   The reason is the forward model's own numerical noise. The scatter of χ² at an optimum is
   about 2 × 10⁻³ at χ² ≈ 1,100, which is about 3 × 10⁻⁵ σ in each residual. Divided by the
   element steps (10⁻¹⁰ to 2 × 10⁻⁷), that noise enters the Jacobian. Multiplied by large
   residuals, it then dominates the predicted drop of the convergence test. The predicted
   drop scales as the steps predict: 1.8 at a quarter of the steps, 0.038 at the plan's
   steps, 0.0012 at four times them.

   | | the plan's Jacobian | along the singular directions |
   |---|---|---|
   | converged, of 11 calls (J on truth d at 9 strengths; N, G on truth a) | 4 | 10 |
   | χ² of the calls that did not converge | 281 to 48,114 | 48,114 |
   | integrations, all 11 calls | 503 | 184 |
   | χ²_min where both converged | the same to 10⁻³ | |
   | Σ ln s_k | | 0.002 to 0.006 lower, at every point |

9. **V6 counts calls that the evidence uses (changes § 5).** V6's 0.1 % is read on the calls
   within Δχ² = 200 of the best point of their scan; the calls beyond are left out of it, and
   their number and how many of them did not converge are reported. Δχ² =
   200 is 50 L² at L = 2, so at every level their integrand is below e⁻²⁵ of the peak. The
   clause that no failed call may lie within Δχ² = 25 of a scan's best point is unchanged.

   The reason is that χ²'s numerical scatter grows with the residual, to a few × 10⁻² near χ²
   = 5 × 10⁴. There the test's 10⁻³ cannot be met with any Jacobian (the one failure in
   item 8's table). An error of that size in a call so far below the peak changes no ln Z.

10. **V4 reads the minimum of the fine scan's spline of χ²** (or the best grid point, if
    lower). The best grid point alone sits up to half a spacing from the minimum: in the
    smoke test, χ² 0.145 at the best point of truth d's own scan, against 1.6 × 10⁻⁵ at the
    truth.

The code also fixes a defect in the runner, which is not a change to the plan: a call that
starts at its optimum (the truth's own model on its noise-free data) now passes the
convergence test without taking a step.
