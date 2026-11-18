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

## Note, 15 November 2026: a runner failure in Stage 1, and its fix

Stage 1 started at 05:45 PST on 15 November and was stopped by hand at 13:28, with one of its
30 scans finished, 16 part done and 13 not started. No gate or reading had been computed.

**The failure.** Eleven scan processes stopped with `elements_to_state: need a > 0 and
0 <= e < 1`, in seven of the planet model's scans: t1_s0_l0_P_d0_b0, t1_s0_l0_P_d2_b0,
t3_s0_l0_P_d0_b0 and t3_s0_l0_P_d1_b2_half twice each (the runner retries a scan once, so
these four were left unfinished), and t1_s0_l0_P_d1_b0, t3_s0_l0_P_d1_b3_half and
t4_s0_l0_P_d0_b0 once each (they resumed). Every one stopped at the same place: the fresh
Jacobian taken along the singular directions of the current one (amendment 1, item 8), after
two rejected steps. The elements are Mercury's, the only body whose elements are fitted, and
its eccentricity was not near either bound: 0.204 at the two failures reproduced, and no
lower than 0.198 at any of the 45,640 points recorded before the stop. The bound was crossed
because the Jacobian being whitened, updated by Broyden's rule through the rejected steps,
had lost a direction: its smallest singular value was 1.3 × 10⁻⁶ in t3_s0_l0_P_d0_b0 and
6.9 × 10⁻⁷ in t3_s0_l0_P_d1_b2_half, against about 10⁵ for the next, so one whitened unit
along it was a step of about 4 × 10⁵ in e and 10³ AU in a, outside the domain on both sides.

**The fix.** A step of a finite-difference Jacobian is taken one-sided into the allowed
region: the elements' domain (a > 0, 0 ≤ e < 1) with contacts. Forward, as before, when the
forward trial is allowed; backward when it is not (as the runner already did for a trial
with no contact); and when neither trial is allowed, the step is halved along the same
direction until one is, the difference then scaled back to the full step. There is no change
of parameters. Each Jacobian column taken any way but the old one is counted in the call's
record (`edge_steps`), so that the readings can say how many such calls there are and how far
each lies from its scan's best point. Where the old code did not fail, the arithmetic is the
old one, bit for bit: at the 24 solutions of the finished scan (t3_s0_l0_P_d1_b3_half) the
Jacobians by the plan's steps and along the singular directions are identical with the old
and the new code. Tests: tests/test_phase2_identifiability_full.py (a near-circular orbit,
e = 10⁻⁴, and the collapsed direction of t3_s0_l0_P_d0_b0, both of which the old code failed
on). Rerun from their recorded points, the two reproduced failures now go on: the call at
the failing point of t3_s0_l0_P_d0_b0 took four columns at the edge and did not converge, at
χ² ≈ 2.96 × 10⁵, with two further unconverged calls beside it, all at the edge of the coarse
scan and far outside Δχ² = 200 of any best point.

The fix changes no model, prior, likelihood, data set, grid, tolerance or gate.

**The restart.** The four scans that failed twice restart from the beginning; their earlier
points are kept aside, not used. The other scans resume from their recorded points. A
resumed scan does not have the Jacobian it would have carried into its next point (it is not
recorded), so that point starts from a fresh one, as the first point of a scan does; each
point is still a converged fit, but a resumed scan is not bit-identical to an unbroken run.
The finished scan, and the three that resumed after one failure, are kept as they are.

## Amendment 2, 15 November 2026: the contact solver, and Stage 1 voided

Written after Stage 1's first launch was stopped and before it is launched again, and approved
with two additions (the carry's range, and the decision on 38bde9f from the measured result)
before it was committed. It changes
text above the line (§ 3's solver, its carry and its "under 30 min"; § 2's inputs gain one
file) and regenerates the truths' data. Nothing in it was set from a Bayes factor, a ln Z or a reading:
none had been computed.

**The defect.** § 3 finds each contact "by Newton's method on the offset" from the printed
epoch, with the bodies carried "over the offset (under 30 min)". The premise is false for the
study's own truths: truth a's contacts lie within 7.5 min of the printed epochs, but truth
b's lie up to 1.1 h from them and truth c's up to 2.0 h (the planet and the ring move the
transits). Over such offsets Newton's method from the epoch fails in two ways.

- Where a model's mid-transit falls near the printed epoch, the gap's slope there is about
  zero and the first step leaves the one-day window: the call is recorded as "no contact"
  though the transit is central. In t3_s0_l0_P_d0_b0 this happened at x = −6.4 (row I 1,
  1677): Mercury passes 12′ inside the limb, both interior contacts about 2.5 h from
  mid-transit.
- Which root Newton's method reaches depends on which side of the epoch mid-transit lies, so
  a row can take the other contact of the same transit: a jump of hours in one residual.

**What it did to Stage 1.** Truth c's data set took the other contact at three rows of 57
(I 5, 1736 November; I 11, 1782 November; II 1, 1740 May). Every planet
scan held spurious "no contact" points, and in seven of them χ² was still falling where such
a run began, within Δχ² = 200 of the row's best; in t3_s0_l0_P_d0_b0 the scan's best was
χ² = 100.0 at x = −6.5, against 32.9 at x ≈ −6.26 found by refitting in steps of 0.01. The
11 crashes of the note above (15 November) come from the same defect: re-run with the
corrected solver and data, under the old step code and the new, both scans pass the
crash points with no edge column, e = 0.20572 throughout, χ² ≈ 485 at x = −11 (it was
2.96 × 10⁵), and the smallest singular value of each call's final Jacobian 1.7 × 10⁴, 0.4
of the next (it was 10⁻⁶ in the collapsed Broyden Jacobians). The jumps in the residuals had
collapsed the rank-one updates.

**The fix.** A contact is a root of the gap (the separation minus the contact radius) as
before, found by bracketing instead of Newton's method (`src/fork/phase2/contact_root.py`):

1. mid-transit is the minimum of the gap within a day of the printed epoch, by
   golden-section search;
2. the row's own contact, ingress (I, II) or egress (III, IV), is the root between
   mid-transit and mid-transit ∓ 1 day, by bisection to a double's precision;
3. there is no contact only when the gap has no root there (no transit, or a grazing one
   whose limbs never reach this contact), or the root lies more than a day from the epoch,
   or |gap| at it is not below 10⁻¹² rad. Item 3 of amendment 1 (no contact, likelihood
   zero) applies to these cases only.

**The carry, and how far it reaches.** Contacts are sought within a day of the printed epoch
and no further: a contact more than a day from it is no contact (step 3), so no offset beyond
a day reaches the likelihood. Within that window, offsets do reach it well beyond 30 min: the
truths' own data lie up to 2.0 h from the printed epochs (truth a, the nearest to the real sky,
up to 7.5 min), and over the 467 states of validation (a), across the four truths' prior
ranges, contacts the solver returns lie up to 23.9 h from the epoch (b and c, at the largest
strengths) and 9.1 h (d); of the 24,739 returned on kept rows, 5,219 lie beyond 2 h and 2,280
beyond 6 h, almost all at strengths far from any truth's best fit. So every offset up to the
window's day reaches the likelihood. The plan's
carry, a second-order expansion in the epoch's Keplerian acceleration, is off by 0.009″ in the
direction from the Earth to Mercury at 2 h (about 0.1 s in a contact), 0.25″ at 6 h and 16″ at
a day, measured against exact two-body motion. A switch to exact motion past a bound such as
2 h would sit inside the range of truth c's own data, so the carry is now **the exact two-body
motion about the Sun at every offset** (Kepler's equation in universal variables, the repo's
`kepler_fg_scalar`). What it still leaves out is the pull of the other bodies over the offset:
against the full integration, at most 8.5 × 10⁻⁵″ in the separation at 2 h, 7.5 × 10⁻⁴″ at 6 h
and 0.013″ (about 0.13 s) at a day, the largest at the prior's edge (a planet of 10⁻⁵ M☉ at
0.15 AU); for truth a, 3 × 10⁻⁶″ at 2 h. The code enforces both halves: no contact is returned
beyond a day (`contact_root`), and a test checks the carry against an independent integration
of the two-body equations to 10⁻¹² AU at 2 h, 12 h and a day.

**III 12** is left out of V0 and of the likelihood: the observation is built from V0's kept
rows only, and a test asserts that III 12's epoch is not among the observation's 57.

**Which contact each row is.** Defined in data terms and frozen in
`data/historical/mercury/newcomb_1882_contact_identity.csv` (SHA-256
`6b03d75401d2bb2b15f31aed11f1f772cfac79871bbf54f6e8b898cd638c6193`, added to § 2's inputs), one row per row of the contact table, with the printed number, the contact
used and its basis:

- tables I and II: the printed number (II or III), except I 11, printed III, which is contact
  II: the misprints table already reads it so (the transit's two rows are both printed III;
  this one is the earlier);
- tables III and IV print no number. Nineteen of their 21 rows share a transit with rows of
  table I or II, and each follows that transit's printed contact III, by 1 to 3 min: contact
  IV. The other two are IV by Newcomb's Part I: 1690 (p. 385, Fontenay at Canton, "Sortie
  certaine et entière", 19h 28m 52s against the table's 19h 28m 54s) and 1756 (p. 390,
  Peking, "Egressus totus", "Sortie totale", "Fin de la sortie", 18h 54m 22s–49s against
  18h 54m 34s).

The solver reads this column and never chooses ingress or egress from a model. The contact
table itself is unchanged (it keeps printed values; its hash is unchanged). As a check, not
a rule, a test asserts that on truth a's noise-free model the contact nearer each printed
epoch is the file's on every row but III 12, and that III 12 is the only row V0 leaves out:
its printed day is one before table I's for the same transit, so both of truth a's contacts
lie 18 to 24 h after the printed epoch.

**Truth data regenerated.** With the file's contacts and the exact carry, truth c changes at
those three rows by 5,100 to 10,100 s each; every other datum moves by the carry alone, by at
most 1.6 × 10⁻⁵ s (truth a), 0.014 s (b), 0.096 s (c) and 0.0013 s (d). V0 is unchanged (III 12
alone left out; largest kept offset 447.2 s). σ₁ changes from 17.079679672 to 17.079680564 s. The real-data module (`contact_real`) uses the same root finder and states each
contact's side; its own carry is unchanged here and is for the real-data plan to set.

**Validation, before any run.**

(a) With the final code (solver, contact file and exact carry), on 467 states (truth a; truths b, c and d's model families over their prior strength
ranges in steps of 0.05, c at its three distances; and t3_s0_l0_P_d0_b0's elements across the
false edge, x = −6.6 to −6.2 in steps of 0.01), every row against an independent search: the
gap on a 1 s grid over ±1 day, every sign change refined by Brent's method, ingress the first
root and egress the second. 27,086 rows: 25,116 with two roots, where the solver's contact is
the grid's on the row's side to 1.5 × 10⁻⁹ s at worst (tolerance 10⁻⁵ s, a change in χ² below
10⁻⁴ for a 10 σ residual); 1,622 with none, where the solver gives none; 348 with one root
inside the window (transits moved about a day at extreme strengths), rechecked on ±2 days: 201
agree to 5.4 × 10⁻¹⁰ s and in 147 the row's contact lies beyond a day, where both give none. No
row has more than two roots. On the same states Newton's method took the other contact at
2,361 rows and gave a spurious "no contact" at 80.

(b) The 2016 and 2019 transits (Espenak's published geocentric contacts, with light time and
aberration, the real-data module): the new and the old solver agree to 0.1 ms, and the model
is within 1.7 s of the printed times (tolerance 2 s, set before).

(c) Every caller of the solver: every likelihood call of the study (Stage 1, Stage 2, V2′'s
importance sampling), `prepare` (V0, σ₁, the truths' data), and the real-data module's tests.
No real-data likelihood exists. With the final code all tests pass (the two modules' 47, slow
ones included, and the fast suite's 733) once this amendment, with the identity file's hash,
is in the plan; before it, the one failure was that check.

**Decision on the 15 November fix (38bde9f).** Measured first: the two scans that crashed
(t3_s0_l0_P_d0_b0 and t3_s0_l0_P_d1_b2_half), run from the start of their coarse scans past
the crash points with the corrected solver and data and the step code from before 38bde9f, run
clean: no crash, no column at the domain's edge, e = 0.20572 throughout, χ² ≈ 485 at
x = −11 (it was 2.96 × 10⁵), and the smallest singular value of each call's final Jacobian
1.7 × 10⁴, 0.4 of the next (10⁻⁶ in the collapsed Broyden Jacobians). The same run with
38bde9f's code gives the same rows, and so does the final code of this amendment (the exact
carry below moves χ² by 0.01). So the jumps in the residuals had collapsed the rank-one
updates, and 38bde9f's extension to calls with no contact was covering for the solver. It is
dropped: a forward trial with no contact takes the backward difference, as before 38bde9f, and
no contact on both sides of a step is an error (the scan stops and is reported), as before.
Step-halving is kept only as a fallback at the edge of the elements' domain (a > 0,
0 ≤ e < 1); where it does not fire the arithmetic is the plan's, bit for bit, and each column
it takes is counted (`edge_steps`) and reported at the gate reading. The note's account of why
the Jacobian collapsed is superseded by this amendment.

**Stage 1's first launch is void.** Started 05:45 PST and stopped at 15:22 PST on 15 November,
it is not used: its points, logs and the four scans set aside at the earlier restart are kept
outside the study's cache (`stage1_void_launch1/`), not deleted. An edge-handling change for
the false cliff was approved earlier on 15 November and then withdrawn before it ran; nothing
from it is used.

**Before the relaunch.** A smoke test in a scratch cache: two or three scans per truth,
including t3_s0_l0_P_d0_b0 (whose first row crosses the false edge) and the first row of the
other scan that crashed; nothing from it is reused. Stage
1 then starts again from nothing, with the regenerated data.

**Unchanged.** The models, priors, likelihood, evidence quadrature, grids, tolerances,
gates, readings and seeds.

## Amendment 3, 16 November 2026: the fine points stay in the box, and a start with no contact is retried

Written after amendment 2's smoke test and before Stage 1 is launched again; approved before
it was committed. It changes how items 3 and 6 of amendment 1 are carried out. Nothing in it
was set from a Bayes factor, a ln Z or a reading: none had been computed.

**What the smoke test showed.** Amendment 2's smoke test (a scratch cache; ten scans, each at
the first value of its second parameter) ran clean on what amendment 2 changed: no crash, no
edge column, the false edge gone (t3_s0_l0_P_d0_b0's best point at x = −6.27, χ² = 28.4). It
also showed two faults that predate amendment 2 and were in Stage 1's first launch too.

- *The fine points had no bound.* Item 6 extends them until both ends are Δχ² = 50 L² above
  the minimum. Where the data do not limit a weak perturber, χ² is flat at low strength and
  that never happens: the extension ran to its cap of 4,000 steps. The first launch's scans
  reached x = −686 (and +13 at the other end); in the smoke test, 10,093 of its 11,970 points
  lay outside the box. Every integral is clipped to the box (item 4), so none of those points
  entered one; they cost hours of computing each.
- *"No contact" was never the result of a fit.* Every no-contact call, 252 in the smoke
  test and 1,001 in the first launch, stopped at its first iteration: the call's starting
  elements gave no contact, so no fit was made, and item 3 set its likelihood to zero. A
  start is not a fit. In t1_s0_l0_P_d0_b0 all six tail points started from the best point's
  elements, fitted with the planet at x = −6.18, and stopped there; their zero likelihood
  passed the tail check, while the coarse point beside one of them (x = −10.2) is fitted, at
  χ² = 1,107.

**A. Nothing outside the box.** The box is the widest any integral uses, the coarse scan's 8
decades: x from −11 to −3 for the ring and the planet, −7.5 to 0.5 for J2. No strength outside
it is computed:

1. the fine points, their extension (item 6) and their widening toward a tail point stop at
   the box: when the next step would pass an edge, the edge itself is the last point;
2. tail points outside the box are not computed, and the tail check uses those inside;
3. the best value of a scan is the best inside the box.

Where the fine points stay inside the box, nothing changes: the same points, the same spline,
the same integrals. Where the old extension passed an edge, the integrals change only through
the spline's end: the fine points now end at the edge itself, and the spline is fitted to the
points inside. Measured on the smoke scans, recomputing every ln Z (each box, each level) from
the same points with those outside the box left out: the largest change is 1.6 × 10⁻⁴ (the
scan whose best value lies below the box, t3_s0_l0_P_d1_b2_half, at the widest box), at most
4 × 10⁻⁷ in the others, and none in the six whose points stayed inside. V3′'s tolerance is
0.05.

(a) *A best value at an edge.* Where the best value inside the box lies within one fine spacing
of an edge (the data pushing toward no perturber, at the lower edge), the fine points are
anchored at that edge, at a spacing no wider than the coarse scan's (0.1 decade; 0.05 for
V3′), so the integrand is dense where it is largest. A coarse point the fine points land on
counts as one of them. Each scan records whether its best point is at an edge
(`best_at_box_edge`), and the gate reading reports how many are, at which edge.

**B. A start with no contact is retried.** A call whose start gives no contact is called again
from the nearest call of the same scan (the same second parameter) that converged with
contacts, then from the base elements, skipping any start already tried. It has no contact,
and item 3 sets its likelihood to zero, only if every start gives none. The contact window
stays a day either side of the printed epoch (amendment 2); nothing here widens it. Each call
records how many starts it retried (`retries`); each result counts the calls retried and
those that then found contacts, and the gate reading reports both.

(b) *The smoke test's zero points under B.* Each of the 252 no-contact calls of amendment 2's
smoke test, called again under B with the first run's other calls as neighbours: 7 become
finite (4 of the 19 tail points, 3 of the 228 coarse points, none of the 5 fine points); the
other 245 have no contact from any start, at strengths that move a transit more than a day.

Where a call's first start gives contacts, nothing changes. A call retried and fitted is a
converged fit like any other; the scan's next point starts from it, so later points in the
same scan can differ from the old run's within the fit's tolerance.

**Tests.** A stub likelihood with its minimum outside the box (no point outside it, the best
point at the lower edge and reported so, fine points there no wider apart than 0.1); one with
its minimum inside (fine points on the plan's x₀ + k h, best value unchanged); a start with
no contact retried from the nearest converged call, and no contact only when the base
elements fail too.

**Before the relaunch.** The smoke test again, all ten scans to the end, including the two
stopped last time (t2_s0_l0_R and t3_s0_l0_P_d1_b2_half). Stage 1 then starts again from
nothing.

**Unchanged.** The models, priors, likelihood, contact window, evidence quadrature inside the
box, grids inside the box, tolerances, gates, readings and seeds.

### Addendum to amendment 3, 16 November 2026: a trial with no contact inside the edge fallback; smoke condition 2

Written after the second smoke test and before the third; approved before it was committed.

**What happened.** In the second smoke test (amendments 2 and 3, ten scans from scratch),
nine scans finished and t3_s0_l0_P_d1_b2_half stopped after 50 min, at x ≈ −7.15
(χ² ≈ 850). A column of a fresh whitened Jacobian had its forward trial outside the elements'
domain, which sends the column to the edge fallback kept by amendment 2; the fallback's
backward trial had no contact, and the fallback, as amendment 2 narrowed it, caught only
trials outside the domain, so the scan stopped. Before 38bde9f the column would have stopped
at the forward trial; 38bde9f's fallback skipped both kinds of trial.

**The fix.** Inside the edge fallback only, a trial with no contact is skipped like one
outside the domain, and the step is halved until a trial is allowed; each column taken there
is counted in `edge_steps`, and a column with no allowed trial is an error. Outside the
fallback nothing changes: a forward trial with no contact takes the backward difference, and
no contact on both sides of a full step is an error. Test: a residual whose forward trial is
outside the domain and whose backward trial has no contact (it raises NoContact on the code
before this addendum, as the scan did).

**Smoke condition 2, made precise before the rerun.** "Every no contact is a genuine miss
after retries" is checked as: every call with no contact was retried from every start B
allows, and none lies below (at a weaker strength than) a call of the same scan that has
contacts and whose fit converged. A call whose fit did not converge is not counted as a
contact neighbour. In the second smoke test the two calls with no contact below the last
call with contacts (x = −5.52 and −4.97, in the scan that stopped) had as neighbours only
calls that stopped at the iteration cap with χ² of 2 × 10⁵ to 5 × 10⁷; the nearest converged
call was at x ≤ −5.7.

**Unchanged.** Everything else in amendments 1 to 3.

### Second addendum to amendment 3, 16 November 2026: the spline's knots, a floor on the fine spacing; smoke conditions 6 and 7

Written after the third smoke test and before the fourth; approved before it was committed.

**What happened.** In the third smoke test all ten scans finished, but in
t3_s0_l0_P_d1_b2_half the cubic spline of item 4 and § 4.2, fitted through every fine point
with a finite χ² (one unconverged call at χ² = 3.3 × 10⁷, and points from several fine grids
as little as 1.6 × 10⁻⁶ decade apart), swung to χ² = −64,798 at x = −6.665. The scan took that
as its best value, with a width of 0.0002 decade, and its refinement chased it, adding points
about 10⁻⁵ decade apart (1,043 calls against about 130 in the other planet scans). The scan's
lowest converged point is at the box's lower edge, χ² = 485.0. Its integrals escaped only
because fine points with no contact switched them to linear interpolation; without those the
spline itself would have been integrated. The other nine scans' spline minima lay within
0.15 of their lowest points.

**The fix.**

1. *Knots.* The spline of the best value, the width, the refinement and the integral (item 4)
   is fitted to the knots only: fine points whose fit converged with contacts, of any two
   closer than 10⁻⁴ decade the one with the lower χ². An unconverged or no-contact fine point
   is never a knot; it is treated as the coarse and tail points are, entering only the linear
   interpolation outside the knots (a no-contact point there as zero likelihood, item 3). The
   best point of a scan is the knot with the lowest χ².
2. *Refinement guard.* The fine spacing of a grid is never below 10⁻³ decade (5 × 10⁻⁴ for
   V3′'s halved grids), and a scan value is refined at most six times (item 6's count, now
   enforced in code as a named limit).

**Smoke conditions 6 and 7, set before the fourth smoke test.**

6. In every scan, the spline's minimum is within Δχ² = 1.0 below the lowest knot's χ², and its
   position within one fine spacing (the scan's `h`) of that knot.
7. In every scan, at every level and for every box, ln Z by the spline integral and ln Z by
   the trapezoid rule on the same knots (the same treatment outside them) agree within 0.05.
   If they do not, the result is reported and no interpolant is changed without a further
   addendum.

Conditions 1 to 5 are as before (condition 2 as made precise in the first addendum; condition
4: edge steps zero or each reported).

**Tests.** A point set with the smoke test's shape: the spline through every finite point
reaches χ² below −10⁴, as it did; from the knots the best value is the lowest point, at the
box's edge, within Δχ² 1 of the spline's minimum, and the two integrals agree within 0.05 at
each level. A likelihood whose width would give a spacing below the floor keeps 10⁻³.

**Unchanged.** Everything else in amendments 1 to 3 and the first addendum.

### Third addendum to amendment 3, 16 November 2026: the box edges in the integral, and a ground-truth check

Written after the fourth smoke test and before the fifth; approved before it was committed.
Nothing in it was set from a Bayes factor or a reading: the only ln Z values computed are the
smoke tests' and the spot check's below, on synthetic data in a scratch cache.

**What the fourth smoke test showed.** Conditions 1 to 4 and 6 held; 5 and 7 failed, both in
the 1-decade box of t3_s0_l0_P_d1_b2_half, where the box leaves out the scan's best value and
its lower edge cuts the likelihood's flank about 230 in ln L below the peak. Two things were
found.

- *The edge was not a point of the integral.* Item 4's "clipped to the box" was carried out by
  dropping grid points outside the box, never by evaluating at the edge. Where ln L falls by
  hundreds per decade, the integral then depends on where the sub-grid happens to fall: in
  that scan the knot at the edge was stored at −7.50000000001 and dropped. Against a ground
  truth (below), the runner's ln Z was off by up to 0.76 in 5 of 54 scan × box × level cases,
  all in the 1-decade box, and by at most 0.011 in the prior box.
- *Linear ln L is not an accurate check at the peaks.* ln L linear between the knots and
  integrated exactly (piecewise-exponential) differed from the ground truth by up to 0.31, in
  23 of the 54 cases, at the peaks: the fine spacing is about 0.86 widths, 1.7 standard
  deviations at level 0.5, and a chord cuts under the parabola. The trapezoid rule on the knots
  is worse on the flanks (up to 2.8). Neither can be the condition.

The ground truth there: 96 to 117 real likelihood fits added to each of three smoke scans
(t3_s0_l0_P_d1_b2_half, t1_s0_l0_P_d0_b0, t3_s0_l0_P_d0_b0), placed as in item 3 below; ln Z
from the scan's points and the extra fits together. Its spline and piecewise-exponential values
agreed within 0.019 in every case. With the edge fix of item 1, the runner (the scan's own
points) agreed with it within 0.013 in all 54.

**1. The box edges are points of every integral.** A converged call with contacts within
10⁻⁹ decade of a box edge (the coarse grid has one on every edge) is that box's knot at exactly
the edge; fine knots within 10⁻⁴ decade of it give way to it. Both edges are points of the
quadrature grid, valued by the same interpolant as their neighbours (the spline inside the
knots, linear ln L outside them). The knot interval containing an edge is sub-divided 200 times
(20 elsewhere); outside the knots, the step at each edge is sub-divided 10 times more finely
(h/200). Nothing else in item 4 changes.

**2. Smoke conditions 5 and 7, restated.** The ground-truth spot check (item 3) is run on three
scans of the smoke test: t3_s0_l0_P_d1_b2_half, t1_s0_l0_P_d0_b0 and t3_s0_l0_P_d0_b0.

7. In each of the three, for every box (the prior box and the five R3 boxes) and every level,
   the runner's ln Z agrees with the ground truth within 0.05; each has at least 50 extra fits
   converged with contacts; and the ground truth is valid, its spline and piecewise-exponential
   values agreeing within 0.05. A case in which every rule gives zero likelihood in the box
   passes.
5. The box-bound integrals, measured against the same ground truth, within V3′'s 0.05: the same
   comparison as condition 7, read for the boxes' bounds; it holds when 7 does.

The piecewise-exponential and trapezoid values on the scan's own points are reported, not
conditions. Conditions 1 to 4 and 6 are unchanged.

**3. The ground-truth spot check.** In one scan (one value of the second parameter), real
likelihood fits, made as the scan's own calls are (started from the nearest converged call with
contacts, then from the base elements), at:

- 0.002, 0.005, 0.01 and 0.02 decade either side of every box edge (half those for V3′);
- 0.2, 0.4, 0.6 and 0.8 of the way along each of the six steepest knot intervals inside the
  widest box (steepest in ln L at the task's lowest level);
- 21 points across ± 2.5 widths of the best value, unless it sits on a box edge;

inside the widest box, none where the scan already has a point. If fewer than 50 converge with
contacts, four more in each of the next six steepest intervals, until 50 do or no interval is
left (then the check fails). The scan's own points are not changed. For every box and level:
the runner's ln Z from the scan's points; the ground truth, ln Z by the same rule from the
scan's points and the extra fits; and the piecewise-exponential value on those same points.
It passes where the ground truth is valid (its two values within 0.05) and the runner is within
0.05 of it.

**4. At Stage 1's gate reading**, the same spot check, with the same 0.05, on this sample,
fixed now:

- every scan whose best point inside the box sits on a box edge;
- every scan with a box that leaves out its best value and whose edge nearest that value lies
  within the knots' range with ln L there, by linear interpolation between the knots, more
  than 10 below the scan's peak ln L at some level;
- five scans drawn by `numpy.random.default_rng(20261116)` (`choice`, without replacement) from
  all of Stage 1's scans other than the no-perturber ones, listed by task identifier and then
  by scan index.

A scan in more than one part is checked once. Any failure is reported at the gate reading with
its scan, box, level and the three values.

**Tests.** A box edge on a steep flank with the edge's knot stored 10⁻¹¹ outside it: the code
before this addendum was 0.18 low there, and now ln Z matches the exact Gaussian integral
within 0.01, as does the box containing the peak.

**Before the relaunch.** The ten smoke scans again from nothing, then conditions 1 to 7 with the
spot check. Stage 1 starts again from nothing only if all hold.

**Unchanged.** Everything else in amendments 1 to 3 and the first two addenda.

### Fourth addendum to amendment 3, 16 November 2026: a deviation, after a crash in Stage 1

Written after Stage 1's relaunch (commit 7f180cf) and approved before it was committed. It is
a deviation from the plan, triggered by a crash: t3_s0_l0_P_d0_b0 failed twice at 06:30 on
16 November, at its second longitude (λ = 7.5°), and the driver gave up on it. It was decided
before this task had any result past 7.5°: its points then were its first longitude's
(computed in full) and the 15 that crashed. It changes the search procedure only; the
likelihood, the priors, the grids of the second parameter, the integrals, tolerances, gates
and readings are as before. The other Stage 1 tasks keep running with the code they started
with (below: what they compute is unchanged).

**What happened.** From the second value of the second parameter on, a scan does not repeat
the coarse scan: the fine points are placed about the previous value's best (item 6). In this
task that was x = −6.26, close to the strength above which a planet at 0° moves a transit more
than a day. At 7.5° all 15 fine points there had no contact from any start (amendment 3, B),
so the grid had no knot (second addendum), and the best value of a grid without a knot was not
defined: the code raised an IndexError. The smoke tests ran only each scan's first value, the
one that starts from the coarse scan, so this path was never run in them. The same task also
failed twice in the void first launch (15 November, 14:30); that was not followed up then.

**The change.**

1. *A guess with no knot.* When the fine points placed from the previous value's best contain
   no knot, that value is scanned as the first value is: the coarse scan over the box, then the
   fine points from its best. The points already computed stay and enter as any point does (a
   no-contact point as zero likelihood, item 3). Each task records the values at which this
   happened (`coarse_restarts`), and the gate reading reports them.
2. *No contact anywhere in the box.* If no coarse call at that value converged with contacts,
   the value is recorded with its points (all zero likelihood, so ln Z = −∞ in every box),
   flagged "no contact in box", without a best value; the next value starts from the coarse
   scan. If the coarse scan has such a call but the fine points from it still have no knot,
   the value is flagged "no knot after the coarse scan" and its ln Z is ln L interpolated
   linearly between its finite points. Each task lists its flagged values (`flagged`), and the
   gate reading reports them.
3. *What it does not change.* Where a fine grid has a knot, every call, point and number is as
   before: the change acts only where the code before it crashed. Checked before the rerun:
   the first three longitudes of a finished task (t3_s0_l0_P_d2_b0) rerun with the new code in
   a scratch cache must give every point and every result field bit-identical to Stage 1's,
   apart from the commit stamps of item 4. If they do not, nothing is relaunched.
4. *Commit stamps.* Every call records the git commit of the code that computed it ("+dirty"
   if the runner's source differed from it); each result row lists the commits among its
   longitude's points, and each task's result the commit that assembled it.
5. *Failures reported.* A task the driver gives up on (two failures) is marked in its directory
   and logged in capitals, and every status check (`identifiability_full plan`) lists it, also
   for a driver started before the marker existed (from its log).

**Before the rerun.** t3_s0_l0_P_d0_b0's first three longitudes (0°, 7.5°, 15°) in a scratch
cache with the new code, checked against smoke conditions 2 and 6. Then its directory is set
aside (not deleted) and the task is run again from nothing with the new code.

**Tests.** The crash's shape (fine points from the previous best with no contact): before
this addendum an IndexError, now the coarse scan, counted. A longitude with no contact anywhere
in the box: recorded, flagged, ln Z = −∞, the next longitude scanned. A task failed twice is
listed by the status check, from the driver's log and from its marker.

**Unchanged.** Everything else in amendments 1 to 3 and the first three addenda.

**Fourth addendum, clarification (16 November 2026, before any task ran with it).** The
addendum above says it changes the search procedure only. That holds for item 1 and for the
"no contact in box" case of item 2, which set nothing that was not already zero likelihood
(item 3). It does not hold for the "no knot after the coarse scan" case of item 2: there the
likelihood at that value is integrated by ln L linear between its finite points, not by the
spline on knots (item 4 and the second addendum), so that case changes how the likelihood is
integrated as well as the search. It is therefore counted apart: each task lists its coarse
restarts (`coarse_restarts`, the search), its values with no contact in the box
(`no_contact_in_box`) and its values integrated linearly (`no_knot_linear`) separately, and
Stage 1's result collects the three. At the gate reading, if any value was integrated
linearly, a sensitivity line gives every ln B against G of the readings with and without those
values (the values dropped from their tasks; a task whose every value is such is kept whole and
named); if none was, the line says so. Invariance (item 3) is checked with the code that
carries this clarification, and the rerun of t3_s0_l0_P_d0_b0 uses it.

### Fifth addendum to amendment 3, 16 November 2026: a second retry at fallback values, and condition 6 in a flat valley

Written after the fourth addendum's check run and approved before it was committed. That run
reran t3_s0_l0_P_d0_b0's first three longitudes in a scratch cache: invariance held
(t3_s0_l0_P_d2_b0's first three longitudes bit-identical), but smoke conditions 2 and 6 failed
at λ = 7.5°, the longitude that fell back to the coarse scan. The task was not relaunched.

**1. A second retry at fallback values (search only).** At 7.5°, 13 calls with no contact lay
among converged calls with contacts: they were the fine points placed from the previous
longitude's best, retried (rule B) before any call at 7.5° had converged, so from that best and
the base elements only; their neighbours fitted later, from this longitude's converged calls,
have contacts (χ² 2,221 to 3,867, against the best 485). Now, at a value that fell back to the
coarse scan, every call with no contact made before that value had a converged call with
contacts is called again under rule B once the fine points are done, from the nearest such call
then there. A call that then has contacts replaces the first; one that has none is recorded
again, with one retry more. Each call so made is marked (`retried_again`), and each task
records per value the calls retried again and the contacts newly found (`retried_again`), which
the gate reading reports. The fine points of that value are then all its fine points. Nothing
changes at a value that did not fall back, so invariance (fourth addendum, item 3) still holds.

**2. Smoke condition 6 in a flat valley (written after the λ = 7.5° result; prospective).** At
7.5° χ² is flat from x = −11 to −10.74 to within 0.01 (the fits' tolerance shows as jitter of
about 10⁻³): the spline's minimum (485.0116 at −10.74) is 0.0004 below the lowest knot (485.0120
at −11), so the Δχ² part holds, but it lies 0.26 decade from that knot against h = 0.054.
In a valley this flat the minimum's position is not determined. From now on: condition 6's
Δχ² ≤ 1.0 part always applies; its position part (within one fine spacing) is waived when every
knot between the spline's minimum and the lowest knot lies within Δχ² 0.1 of the lowest knot.
Every waived case is reported, at the smoke check and at the gate reading (Stage 1's result
lists condition 6's failures and waived cases for every scan value).

**Tests.** Calls with no contact before a fallback (a stub in which only a start from a
converged call at the same longitude has contacts) are retried again, find contacts, are marked
in the point file and counted, and no call with no contact is left below a converged call with
contacts. Condition 6: the 7.5° shape is waived and passes; a knot between rising by 0.5 fails;
a spline minimum 1.5 below the lowest knot fails; a minimum within h passes unwaived.

**Before the relaunch.** The three longitudes again from nothing, invariance on the same check,
and conditions 2 and 6; the task is relaunched as planned only if both hold.

**Unchanged.** Everything else in amendments 1 to 3 and the first four addenda.

**Fifth addendum, note: smoke condition 2 restated (16 November 2026, written after the λ = 7.5°
rerun).** The rerun of t3_s0_l0_P_d0_b0's first three longitudes with the fifth addendum (commit
c8ddc32) gave invariance bit-identical and condition 6 met (waived at 7.5° as a flat valley),
but one call still failed condition 2: the coarse point at x = −6.2 at 7.5°, made by the
fallback's own coarse scan and retried from the nearest converged call then there (x = −6.3),
with no contact; the second retry later found contacts on either side of it (χ² 3,435 at −6.2076
and 3,655 at −6.1904, against the best 485). Rule B retries once, from whatever converged
neighbour exists when the call is made, so whether such a call has contacts can depend on the
order of the calls. Rather than widen the retries again, condition 2 is restated, and applies
from now on: a call with no contact below a value's highest converged call with contacts blocks
only where its nearest converged neighbours with contacts (the lower χ² of the one on each side)
lie within Δχ² = 200 of that value's best. There any point's weight could matter. Beyond it a
point's weight relative to the best is below e⁻¹⁰⁰ at level 1 and e⁻²⁵ at the largest level, 2,
far below the 0.05 in ln Z the integrals are held to, so such calls are counted and reported
at the gate reading (Stage 1's result lists both kinds for every scan value) but do not block.
The call at −6.2 is of the second kind. Test: that case counted and not blocking; a no-contact
call between neighbours 115 above the best blocking.

**Amendment 3, sixth addendum: the spline across wide knot gaps (16 November 2026, written after
the Stage 1 ln Z values were seen and before any gate was computed).** Reading the 24 finished
Stage 1 tasks showed that where a scan's knots jump a steep stretch (at t2_s0_l0_P_d0_b0,
λ = 5.63, from χ² 491 at x = −7.68 to 52,751 at −6.26) the cubic spline through ln L rose far
above every knot inside the gap, and the integral took the bulge: ln Z at level 1 in the prior
box read +16,941 against −234 by the exponential-between-points rule. Since a model's ln Z sums
over its rows, one such row decides it. Smoke condition 6 had flagged these rows; none of the
smoke scans had such a gap.

The rule. A knot interval is integrated by the exponential-between-points rule (ln L linear
between its two knots, exp integrated exactly, on the same sub-grid), not the spline, when
either (a) flank: the interval's higher knot is more than 2 below the scan's best knot in ln L
and the spline rises more than 0.05 above that higher knot; or (b) anywhere: the spline exceeds
the scan's best knot by more than 0.5 in ln L. The scan's best knot is the largest ln L among its
knots (and a box-edge knot) at that level; the spline's largest value in an interval is taken at
its knots and at the roots of its derivative. Every other interval keeps the spline, so the peaks
stay as the ground-truth spot check validated them, and a row where the rule selects nothing is
integrated bit for bit as before (tested). The rule is applied per level and per box.

Also: (1) V4 takes its minimum from the corrected interpolant: the χ² spline on the knots, linear
in the intervals the same rule selects for ln L = −χ²/2 (with ln L = −χ²/2, (a) is a knot more
than 4 above the best in χ² and a dip of more than 0.1 below it, (b) a dip of more than 1.0 below
the best), never a spline value from such an interval. (2) A ground-truth check on the four rows
flagged (t2_s0_l0_P_d0_b0 at λ 5.63 and 6.15, t1_s0_l0_R at r 0.14, t1_s0_l0_P_d2_b0 at λ 1.70):
real fits inside each gap the rule selects (60 evenly spaced, from the gap's better knot; 60 more
at the midpoints, and once more, until at least 50 converge with contacts), then every box and
level: the corrected ln Z from the scan's points within 0.05 of the ground truth (the same rule
on the scan's points and the fits together, valid where it agrees with the exponential-between-
points rule within 0.05). Reported before the gates run. (3) A narrow R3 box with no knot inside
it is reported "not determined"; it is listed per row. (4) Every ln Z is recomputed from the
existing points (no scan is rerun); each row keeps the spline-only values beside the corrected
ones and the result as it was is kept beside it; every changed value is reported with its size.
Test: the λ = 5.63 row from Stage 1's points, the spline alone reproducing the stored +16,941.

**Amendment 3, seventh addendum: filling the gaps that carry weight (17 November 2026, written
after the sixth addendum's gap check and before any gate was computed).** The gap check (commit
671f469) found the corrected ln Z equal to its ground truth to 10⁻⁴ in the prior box where the
selected gap lies far below the best (t2_s0_l0_P_d0_b0 at λ 5.63, t1_s0_l0_R at r 0.14,
t1_s0_l0_P_d2_b0 at λ 1.70), but low by 0.35, 0.42 and 0.16 (levels 1, 0.5, 2) at
t2_s0_l0_P_d0_b0, λ 6.15, whose two gaps lie within Δχ² 6 of the best: there the straight line
under-integrates, and the check's fits brought the integral to the ground truth. Narrow boxes
inside the gaps were off by up to 12,685 either way.

The rule. A knot interval the sixth addendum's rule selects (any box, any level) whose higher
knot lies within 25 of the scan's best knot in ln L at that level carries weight, and is filled
with real fits: 60 evenly spaced in the interval, marching from its better knot, each started
from the nearest converged call with contacts (then the base); 60 more at the midpoints, and once
more, until at least 50 converge with contacts. The fits join the scan's points as fine points,
stamped "gap fill" with the commit that computed them, and ln Z is recomputed (the spline now
spans them). A selected interval whose higher knot is more than 25 below the best keeps the
straight line, which the gap check found exact in the prior box. Where the gap check has already
fitted a weighted gap, its fits are the fill: they were placed by the check's pre-set rule
(the one above, fixed in the sixth addendum before the check ran) and keep its commit. By this
rule the check's rows at r 0.14 and λ 6.15 are filled; λ 5.63 and λ 1.70 are not (their gaps'
higher knots lie more than 25 below the best at every level), and their check fits stay in the
check's record only. A fill call landing exactly on an existing point is not added (counted).

Validity: in every filled row, the corrected ln Z in the prior box within 0.05 of the
exponential-between-points rule on the filled points, at every level, and every gap with at
least 50 converged fits with contacts; any weighted interval the rule still selects after the
fill is reported. If any prior-box check fails after the fill, Mercury is paused and reported
for a design review of the fine-grid placement; no further addendum is drafted.

Narrow boxes. A narrow box overlapping an interval the rule selects, or a filled gap, has a value
at a level only where the spline alone and the exponential rule agree within 0.05 there;
otherwise it is "not determined" at that level (listed per row with the difference). Inside a
selected interval the corrected value is the exponential rule's, so the spline alone is the
comparison. The six tasks that finish after this addendum are recomputed and filled the same way
before the gates run.

## Amendment 4, 17 November 2026: the strength integral by adaptive refinement

Written after the seventh addendum's stop rule paused Stage 1 (two filled rows failed the
prior-box check) and before any gate was computed or any refinement fit made. It replaces the
integration rules of the sixth and seventh addenda (the straight-line fallback, the gap fill and
their checks); their records, fits and recomputed results are kept as they are, the values they
gave kept in each result as `ln_z_before_a4`. Nothing else in the plan changes.

### A. What each gate needs from ln Z (the tolerance)

A row is one strength scan at one radius or longitude; its integral is z_j. A model's ln Z is
ln Σ_j w_j e^{z_j} (the planet: equal weights over 3 distances × 48 longitudes; the ring: the
trapezoid over radii; J: one row). An error δ_j in each row moves the model by
ln Σ_j π_j e^{δ_j} ≈ Σ_j π_j δ_j, where π_j = w_j e^{z_j} / Σ w e^{z} is the row's share.

- V0, V1, V5: no evidence. V2′: Laplace against importance sampling at single points; no
  strength integral. V6: counts calls; the new fits are calls and count as such.
- V4: a χ² minimum, not an integral. It reads the smallest χ² among the row's calls; where the
  χ² spline's minimum lies more than 0.05 below the smallest call (a tenth of V4's 0.5), a call
  is made there first.
- V3′: |ln Z(full) − ln Z(half)| < 0.05 for each truth's own model at each level. Both sides
  carry an integration error of at most τ_m, so their difference carries at most 2τ_m. For
  integration error to move the measured difference by at most a fifth of the tolerance,
  2τ_m ≤ 0.05 / 5, so **τ_m = 0.005** in each model's ln Z. With a refined strength integral V3′
  tests the radius and longitude grids and the integration tolerance; the halved strength
  spacing of its tasks remains as their initial calls.
- Readings: grades change at |L| = 1, 3 and 5, and R1 needs a lead of 3. L is a difference of
  two models, so its error is at most 2τ_m = 0.01; a reading within 0.01 of a boundary is
  reported as at the boundary. Not binding.

Row tolerance. With N rows in a model's sum, row j is held to
τ_j = clip(τ_m / (2 N π_j), τ_m / 2, 0.1). Rows held to τ_m/2 contribute at most τ_m/2 in all
(Σπ ≤ 1); each relaxed row contributes at most 1.05 τ_m / (2N) (e^0.1 − 1 = 0.105), so the
model's error is at most 1.03 τ_m. A row's share is the largest over boxes and levels, from the
larger of its spline and exponential values (an upper bound), and is recomputed after each pass.

### B. Refinement

1. **Initial calls.** Every call already made at the row's value that converged with contacts
   is a knot, of any kind (fine, coarse, tail, gap fill), deduplicated within KNOT_TOL as now;
   the gap check's fits (sixth addendum, item 2) join their rows as calls. No scan is rerun.
2. **Error estimate.** For each knot interval inside a box, e_k = |∫ spline(all knots) −
   ∫ spline(every other knot)| over the interval, divided by the row's integral over the box,
   at each level; the largest over levels counts. (Every other knot counted from the first, the
   last always kept; the coarser spline is evaluated over the same interval.)
3. **Overshoot guard.** Any interval where the spline on all knots rises more than 0.05 above its
   higher knot while that knot is more than 2 below the row's best knot in ln L, or more than 0.5
   above the best knot, at any level, gets a real fit at its midpoint whatever its estimate, until
   the guard holds nowhere (the sixth addendum's limits, now a trigger for fits).
4. **Step.** Every interval with e_k > ε gets a real fit at its midpoint (warm-started from the
   nearer knot, from the other on no contact, then from the base elements). ε starts at τ_j / 4
   and halves when no interval exceeds it while Σ_k e_k ≥ τ_j. Repeat until Σ_k e_k < τ_j in
   every box and the guard holds nowhere.
5. **A midpoint that fails.** If it does not converge or finds no contact, the quarter points are
   tried once; if they fail too, the interval keeps its estimate and is listed.
6. **Cap: 150 fits per row.**
7. **The integral** is the spline on the refined knots, each knot interval in 64 steps on which
   the exponential of the chord is integrated exactly; outside the outermost knots the old
   linear rule on the other finite calls. No fallback interval remains.

### The share floor

A row whose share is below **10⁻⁶ in every box and level** gets the overshoot guard (B.3) only,
not the error-driven steps (B.4). After the guard pass every share is recomputed, and any row now
at or above 10⁻⁶ in some box and level gets full refinement; the shares are recomputed again
after each refinement pass (at most five times) and a row is refined again wherever its tolerance
tightened or it crossed the floor.

The bound. A model's sum has at most N = 144 rows, so the rows under the floor carry a share
below 1.44 × 10⁻⁴ in all. An error δ_j in such a row moves the model's ln Z by
ln(1 + Σ π_j (e^{δ_j} − 1)). An integral too large (δ_j > 0 from overshoot) is what the guard
removes, and whatever remains is already inside the computed share, so it moves the model by
less than 1.44 × 10⁻⁴; an integral too small moves it by less than 1.44 × 10⁻⁴ (e^{δ} − 1),
which is below **10⁻³** for every δ_j < 2.0 (a row underestimated by less than a factor 7.4). The
six named validation rows (D) all lie under the floor; their values after the guard pass alone
are compared with the ground truth and reported as the measured δ (reported, not a pass test).

### C. The plateau toward the lower edge

Where a knot lies on the box's lower edge (strength → 0, no perturber), the plateau is the
longest run of consecutive knots from that edge whose ln L at the smallest level lies within p of
the first knot's, if the run has at least 3 knots. It is integrated as L₀ × width (L₀ the mean of
L over the run's knots, width from the edge to the run's last knot), and only the intervals
beyond it are refined. The error of that product is at most (e^p − 1) times the plateau's share
of the row, so p = min(0.01, τ_j).

### D. Validation, before any gate

Ground truth by dense real fits, independent of the refinement: a fit every 0.01 decade where
ln L at the smallest level lies within 25 of the row's best knot, every 0.05 decade elsewhere in
the widest box, with the row's calls from before this amendment. It is valid where its spline and
exponential values agree within τ_j / 2; a row passes where the refined ln Z lies within τ_j of
it at every level in the prior box (where the ground truth must be valid), and in each R3 box
where both are determined. Rows: t2_s0_l0_P_d0_b0 at λ 6.02, t4_s0_l0_P_d1_b0 at λ 0.92, the
four gap rows (t2_s0_l0_P_d0_b0 at λ 5.63 and 6.15, t1_s0_l0_R at r 0.14, t1_s0_l0_P_d2_b0 at
λ 1.70), and 5 rows drawn by the seed 20261117 from the rows with a share of at least 10⁻³ in
some box and level after the refinement (sorted by task and value, the named rows left out).
The six named rows are refined in full whatever their share, so that the validation tests the
refinement itself.

### E. Stop rule, and who starts the gates

If a validation row fails, or a row with a share of at least 10⁻³ in some box and level ends its
refinement above its tolerance (at the cap, or with an interval whose fits failed), Stage 1
pauses and is reported; no further amendment is drafted. Rows with smaller shares that end above
tolerance are listed with their shares.

The gates run only on an explicit go: `stage1` computes no gate, schedules nothing and reads
nothing unless the file `go_stage1_gates` exists in the study's cache, written by hand after the
validation and the refinement are reported. No driver or queue writes it.

### F. Cost

Simulated on the 625 finished rows' knots (the workings are kept with the review records): with
the floor, about 2,100 error-driven fits in the 178 rows at or above 10⁻⁶ (none at the cap) and
about 900 guard fits in the rows below it, plus the six named rows in full; the 117 rows of the
three unfinished tasks add about a fifth. In all about 3,600 fits, at 9 to 21 CPU-s a fit
(21 s measured on a named row with other jobs running): 9 to 21 CPU-h, 1 to 2 h on 12 workers.
The validation: about 4,000 fits, 10 to 23 CPU-h. The three unfinished scans: about 11,000
calls, about 20 CPU-h, about 11 h (one process per task).

## Deviation 1, 17 November 2026: outcome of the design review required by the seventh addendum: adaptive integration

**Label.** The plan defines no kind of change made after results. Amendments 1 to 3 were all
written before results: amendment 1 "after the runner's smoke test and before any scan in the
study's cache", amendment 2 "after Stage 1's first launch was stopped and before it is launched
again", amendment 3 "after amendment 2's smoke test and before Stage 1 is launched again". From
the sixth addendum on, changes were written after Stage 1's per-row ln Z values were seen, and
the seventh addendum closes further addenda. This change comes after those values were seen, so
it is logged as a deviation, with the table in § 6, not as an amendment.

The section above headed "Amendment 4" (commit ff9a51b, 17 November, 11:42) is this rule's first
form. It was committed before this label was ruled and stands as committed; it is superseded by
this section and is not in force on its own. Where the two differ, this section holds.

### 1. Reported regardless of the outcome

The seventh addendum's prior-box check failed in two filled rows; its stop rule paused Stage 1
on 17 November at 01:01. The values (corrected ln Z minus the exponential-between-points rule,
over the prior box; tolerance 0.05):

| row | level 1 | level 0.5 | level 2 |
|---|---|---|---|
| t2_s0_l0_P_d0_b0, λ 6.0214 (345.0°) | −231.1537 − (−231.4837) = +0.3300 | −674.9857 − (−676.0386) = +1.0528 | −116.5174 − (−116.6108) = +0.0934 |
| t4_s0_l0_P_d1_b0, λ 0.9163 (52.5°) | −452.0893 − (−452.1874) = +0.0981 | −1553.4040 − (−1553.8124) = +0.4085 | −172.7674 − (−172.8101) = +0.0427 (passes) |

**At both longitudes the likelihood peaks below the prior box.** At every level ln L is largest at
x = −11.0, the scan box's lower edge, two decades below the prior box's −9 (smallest χ² 289.97 in
t2 and 709.13 in t4, both at the lowest strengths scanned; the calls stay within 20 of that
minimum from x = −11 up to −8.53 in t2 and −9.07 in t4). Inside the prior box the integrand is then
largest at its lower edge, and x = −9 is not a knot: the nearest knots are −9.068 and −8.527 in t2,
−9.065 and −8.943 in t4. In t4 the result's `x_best` field reads −7.178 with `best_at_box_edge`
"low"; the smallest χ² among its calls lies at −11.0. These facts are stated whatever the values
after this rule are.

### 2. The rule

Applied to every row of every Stage 1 task (every value of the second parameter of every model
with a strength scan), whether its checks passed or failed. The data, likelihood, models,
priors, boxes, levels and estimand are unchanged. The integral is unchanged too: the plan's
rule (`ln_z_strength`, the spline on the box's knots with amendment 3's addenda). The rule
only adds real fits to a row's points; ln Z is then recomputed by that rule.

1. **Pins at every integration boundary.** The boundaries are the edges of the prior box and of
   the five R3 boxes (`r3_1`, `r3_2`, `r3_4`, `r3_6`, `r3_8`): for the planet and the ring −11,
   −10, −9, −8, −7.5, −6.5, −6, −5, −4, −3 (ten; the prior box is `r3_4`); for J2 −7.5, −7,
   −6.5, −5.5, −4.5, −4, −3, −2.5, −2, −1.5, −0.5, 0.5 (twelve). A real fit is made exactly
   at each boundary unless a call within 10⁻⁹ decade of it converged with contacts, or a pin was
   already tried there. It is a fine point stamped "edge pin" with its commit; one that
   converges with contacts is a knot on the edge (the third addendum's edge rule then uses
   it); one with no contact from every start is recorded and not tried again.
2. **Refinement** as in amendment 4, B.2 to B.6 (the error estimate from every other knot, the
   overshoot guard, midpoint fits until Σ e < τ_j in every box and the guard holds nowhere,
   quarter points after a failed midpoint, a cap of 150 fits per row, the pins not counted),
   on the box's knots as the plan defines them (fine points that converged with contacts, and
   the box-edge calls). Refinement fits are fine points stamped "refine". The guard is the
   sixth addendum's selection itself, so once it holds nowhere that addendum's rule selects
   nothing and the integral is the plain spline.
3. **The convergence criterion** is amendment 4's tolerance: τ_m = 0.005 in each model's ln Z
   (from V3′), each row held to τ_j = clip(τ_m / (2 N π_j), τ_m / 2, 0.1).
4. **The share floor** as stated under amendment 4: below 10⁻⁶ in every box and level, the pins
   and the guard only; shares recomputed after that pass and after each refinement pass, and
   any row then at or above 10⁻⁶ gets full refinement. The six named validation rows are
   refined in full whatever their share.
5. **Kept from amendment 4:** the gap check's fits join their rows (B.1); V4 (A); the
   validation (D); the stop rule and the go file (E).
6. **Not used from amendment 4:** its own integrator (B.7), its knots of any kind (B.1: coarse
   and tail calls stay non-knots, as the plan has them), and the plateau (C).
7. **Unresolved boxes.** A box whose estimate ends at or above its tolerance, or where the guard
   still selects an interval, is listed per row and level with its values (`above_tolerance`,
   `guard_left`); an R3 box so listed is "not determined" under the seventh addendum's
   narrow-box rule. No value is changed by hand.

### 3. The smoke test (its conditions fixed now, before it runs)

From the committed code, into a fresh copy of the cache (the inputs only: `prepared.json`,
`scans/`, `gapcheck/`). Rows: t2_s0_l0_P_d0_b0 at λ 6.1523 (the filled row that passed the
seventh addendum's check), and the two that failed it (t2_s0_l0_P_d0_b0 at λ 6.0214,
t4_s0_l0_P_d1_b0 at λ 0.9163). All three are named validation rows, so each gets the pins and
the guard, then full refinement at its tolerance. The results are not rewritten.

- **(a) Bit-identity.** A stub row with a knot on every boundary and no interval above its
  estimate gets no pin and no fit, and its ln Z in all six boxes at every level is bit-identical
  to the value from the same points before the rule.
- **(b) The passing row.** After the rule it still passes the check in all six boxes at every
  level: the spline within 0.05 of the exponential-between-points rule on the same points.
- **(c) Corrections.** Any change in ln Z above 0.05 on a row that had not failed a check is
  flagged as a correction and listed by row, box and level. Every change is reported with its
  size. Flags are reported, not gated; the same holds in the full recompute.
- **The two failing rows:** the check in every box and level after the rule, the pins and fits
  used, and whether each ends below its tolerance within the cap.

The smoke test is reported and nothing further runs until the go.

### 4. Order after the smoke test, each step on the go

The full recompute (every finished row, by § 2; every change and correction reported); the
seventh addendum's check on every row in every box; resuming the three suspended scans
(t3_s0_l0_P_d0_b0, t4_s0_l0_P_d0_b0, t2_s0_l0_R_half), whose rows then get this rule; the
validation; and the gates, which `stage1` computes only when the go file exists.

### 5. Cost

On the 625 finished rows: 5,099 pins (counted on the cache); refinement about 3,000 fits (by
the simulation on the existing knots, before the pins); together about 8,100 fits, at 9 to 21
CPU-s each 20 to 47 CPU-h, 2 to 4 h on 12 workers. The smoke test: about 27 pins and up to 450 fits, under 3 h on 3 workers.

### 6. Deviation table (the format of Willroth & Atherton, 2024, *Advances in Methods and Practices in Psychological Science*)

| original plan | deviation | reason | timing | effect on the estimate |
|---|---|---|---|---|
| Item 6 and amendment 3 with its addenda: fine points placed from the previous value's best at a spacing from the width estimate; a box edge is a knot only where a call lies on it; the sixth addendum's straight-line rule and the seventh's gap fill where the spline overshoots. | Real fits pinned at every integration boundary of every row (ten for the planet and the ring, twelve for J2), then midpoint fits until a share-weighted error estimate lies below τ_j (τ_m = 0.005) and the overshoot guard selects nothing; rows below a share of 10⁻⁶ get the pins and the guard only. The integral itself is unchanged. | The seventh addendum's prior-box check failed in two rows, and its stop rule calls for a design review of the fine-grid placement. In both rows the likelihood peaks below the prior box and the box's lower edge is not a knot. | Written on 17 November 2026, after Stage 1's per-row ln Z values had been seen (from the sixth addendum, 16 November), after the check's failures (17 November, 01:01), after a simulation of the refinement on the finished rows' knots and their shares, after a rehearsal of the code on a scratch copy of the cache (§ 7), and after amendment 4's text was committed (ff9a51b). Before any integral by the committed rule on the study's cache, and before any gate, reading or Bayes factor. | Not known before the run. Every row's ln Z in every box may change; every change is reported with its size beside the value before, and changes above 0.05 on rows that had passed are flagged as corrections. Readings use only the values after. |
| The plan's change log: amendments 1 to 3, each written before results; then addenda to amendment 3, the seventh closing further addenda. | This change is recorded as Deviation 1, the outcome of the design review the seventh addendum requires; the section headed "Amendment 4" above is its superseded first form. | The plan defines no kind of change made after results. | As above; the label was ruled after ff9a51b was committed. | None on the estimate; a matter of labelling, recorded. |

### 7. What was run and seen before this was written

- The refinement simulation on the 625 finished rows' existing knots (no fits): fit counts by
  row and the rows' shares of their models' sums; 178 rows at or above 10⁻⁶ in some box and level.
  Model sums were formed only to compute shares; no model's ln Z was printed or compared with
  another's.
- A rehearsal on code not yet committed, on a scratch copy of the cache, 17 November 11:34 to
  11:38. A rehearsal, not validation; its outputs are kept, labelled as such, and were seen:
  the J row t1_s0_l0_J refined in full needed no fit and its ln Z moved by at most 10⁻⁴; on
  t4_s0_l0_P_d1_b0 at λ 0.9163, with the cap set to 8, 8 real fits were made (21 CPU-s each, one
  not converged, followed by its quarter points).
- No gate, reading or Bayes factor has been computed.

**Deviation 1, note (17 November 2026, written after its smoke test and before any recompute).**

*The smoke test* (commit 3d69dd2, a fresh copy of the cache, 11:55 to 12:05): (a) holds (no fit,
bit-identical); (b) holds (t2_s0_l0_P_d0_b0 at λ 6.1523: +0.044, +0.013, +0.016 in the prior box,
every box passing); (c) no correction (that row's largest change 0.0009). The two rows that
failed the seventh addendum's check still fail it, by less: t2 at λ 6.0214 +0.059 and +0.145
(levels 1 and 0.5), t4_s0_l0_P_d1_b0 at λ 0.9163 +0.070 (level 0.5) and in `r3_2`. Their spline
values moved by at most 0.018; the exponential-between-points values moved by 0.26 to 0.89
toward them. Both rows ended below their tolerance (0.1; shares near 10⁻¹⁵ and 10⁻³⁹). Pins
at −5, −4 and −3 found no contact in all three rows. In t4 the pin at −9 did not converge, and
its narrow boxes rose sharply (`r3_1` by +5,524 at level 0.5), the old values having been
interpolated linearly across far-apart points.

*Rulings, in force from this note:*

1. **A fix.** Reading a row's points back, a later record at the same point was dropped, so a
   pin made where the scan had a failed call was hidden and made again in the next pass (in the
   smoke test, at −3 in all three rows, every such pin failing too, so no value was affected).
   Now, at one point, the record that converged with contacts is kept; where none did, the
   later one.
2. **A failed pin is retried once at once, under rule B** (amendment 3, B): from the nearest
   call that converged with contacts on each side, nearer first, then the base elements, taking
   the first start that converges with contacts. It is stamped "edge pin retry". A boundary
   still without a converged call is listed as unresolved, per row.
3. **What "pass" means for a row:** the validation's dense-fit ground truth (D), not the
   seventh addendum's spline-against-exponential check, which is reported beside it.

*Facts stated for the record.* The prior box is log₁₀ mass −9 to −5 for the planet and the ring
and log₁₀ J2 −7 to −2 (the earlier plan's priors, adopted by this plan); −11 to −3 (−7.5 to 0.5
for J2) is the scan box, the widest R3 box. No prior-box integral includes strengths below its
lower edge: the integration grid and the exponential rule's nodes are cut at the box's edges.
The scan's best value (`x_best`) is the minimum of the χ² spline through the fine knots; in five
rows (t2_s0_l0_P_d0_b0 at λ 5.63 and 6.15, t1_s0_l0_R at r 0.14, t1_s0_l0_P_d2_b0 at λ 1.70,
t4_s0_l0_P_d1_b0 at λ 0.92) that spline dips far below every knot across a wide gap, and
`x_best` lies 3.3 to 4.1 decades from the lowest call (at −11.0). V4 never read it (it reads the
smallest knot and the sixth addendum's corrected minimum, −11.0 in all five), nor does any
part of deviation 1; it placed the next value's fine scan.

**Deviation 1, second note (17 November 2026, written after the first note and the retry of
t4's pins, before any recompute on the study's cache).**

*The deviation table's citation, checked.* Willroth, E. C., & Atherton, O. E. (2024). Best laid
plans: A guide to reporting preregistration deviations. *Advances in Methods and Practices in
Psychological Science*, 7(1), 1–14. https://doi.org/10.1177/25152459231213802. Its template is
Table 1 (p. 8), "Example Preregistration Deviations Table Template": one row per deviation, with
the columns No., Details (type, reason and timing, each from the categories of its pp. 9–10:
Tables 2 and 3 and the five timings), Original wording, Deviation description and Reader impact;
unregistered steps follow in a section of their own. The table in § 6 above has other columns
(no type, the plan paraphrased rather than quoted, "effect on the estimate" for reader impact),
so it does not follow that template. It stands as written; the deviations in that format:

| No. | Type | Reason | Timing | Original wording | Deviation description | Reader impact |
|---|---|---|---|---|---|---|
| 1 | Analytic approach | Preregistered plan not possible or inappropriate because of data | After results of the preregistered approach were known (Stage 1's per-row ln Z values and the seventh addendum's check; no gate, reading or Bayes factor) | "a fine scan of 15 points over ±6 posterior widths about the best value, χ² interpolated by a cubic spline and the integral of exp(−χ²/2) taken on a sub-grid 20 times finer" (§ 3, the strength) | Real fits pinned at every integration boundary of every row, then midpoint fits until a share-weighted error estimate lies below τ_j (τ_m = 0.005) and the sixth addendum's overshoot rule selects nothing; rows below a share of 10⁻⁶ get the pins and that guard only. The integral's formula is unchanged. | Every row's ln Z in every box may change; each change is reported with its size beside the value before, changes above 0.05 on rows that had passed are flagged as corrections, and readings use only the values after. Because the rule was written after per-row values were seen, a reader may weigh Stage 1's readings as less protected against data-dependent choices than those of a plan followed as written. |
| 2 | Inferential criteria | Researchers gained new knowledge | After results of the preregistered approach were known (the smoke test of deviation 1) | "in every filled row, the corrected ln Z in the prior box within 0.05 of the exponential-between-points rule on the filled points, at every level" (seventh addendum, validity) | A row passes against the validation's dense-fit ground truth (amendment 4, D, kept by deviation 1); the spline-against-exponential check is reported beside it. | In the smoke test the two rows that failed the old check moved toward passing it mostly because the exponential rule's value moved; judging against dense real fits tests the integral itself. Rows failing the old check are listed whatever the validation finds. |
| 3 | Variable operationalization or computation | Researchers gained new knowledge | After results of the preregistered approach were known | "The spline of the best value, the width, the refinement and the integral (item 4) is fitted to the knots only" (second addendum, item 1) | In the results, the best value reported (`x_best`) is the argmin of the corrected interpolant on the row's points (the sixth addendum's, as V4 already read); the scan's own value is kept as `x_best_placement`. Placement of the next value's scan is unchanged for the rest of Stage 1. | Five rows' reported best values move from −6.9 to −7.7 to −11.0, where their lowest calls lie; no integral, gate or reading used the value. |

The labelling of this change as a deviation (§ 6 above, second row) is not a deviation of the
plan's content and is not repeated here.

*Rulings, in force from this note.* (1) The three suspended Stage 1 scans (t3_s0_l0_P_d0_b0,
t4_s0_l0_P_d0_b0, t2_s0_l0_R_half) resume on their existing code (commit ca8ea04), with the
same placement rule as the other 27 tasks; their rows get deviation 1 later, with the rest.
(2) Validation, in this order, before anything touches the study's cache: deviation 1 over every
finished row in a scratch copy of the cache, then the 11-row validation exactly as registered
(amendment 4, D, as kept by deviation 1) on that copy; both are reported first. (3) Smoke
condition 6 and the spot-check placement keep reading the scan's own best value
(`x_best_placement`). The gates run only on the go file.

**Deviation 1, third note (17 November 2026, written after the second note, before any
recompute on the study's cache or its scratch copy).**

| No. | Type | Reason | Timing | Original wording | Deviation description | Reader impact |
|---|---|---|---|---|---|---|
| 4 | Variable operationalization or computation | Researchers gained new knowledge | After results of the preregistered approach were known (the per-row best values and the second note's count) | "The spline of the best value, the width, the refinement and the integral (item 4) is fitted to the knots only" (second addendum, item 1) | Applied to every row with a best value, whatever its result: where the argmin of the corrected interpolant lies at, or within one fine spacing h of, an edge of the row's scanned range (its lowest or highest call), the results report the best value as a bound, "≤ low edge" or "≥ high edge", flagged "edge" (`x_best` empty, `x_best_bound`, `x_best_flag`), not as a point; the argmin itself is kept as `x_best_argmin` and the scan's own value as `x_best_placement`. Elsewhere `x_best` is the argmin, as in the second note. | 265 of the 625 finished rows are reported as bounds, all at the low edge: 247 planet rows (≤ −11.0), 17 ring rows (≤ −11.0) and 1 J2 row (≤ −7.5); 360 keep a point value. The rows' scanned range, the scan box and the knots' range give the same 265 rows. No integral, gate or reading uses the value; a mass or J2 near the edge is stated as an upper limit, not an estimate. |

Counted on the points as they stand; the recompute adds fits and the count is given again with
its results. The five rows of the second note, before and after:

| row | value | x_best before | argmin of the corrected interpolant | reported | flag | h |
|---|---|---|---|---|---|---|
| t1_s0_l0_R | r 0.1400 | −6.907 | −11.000 | ≤ −11.0 | edge | 0.0374 |
| t1_s0_l0_P_d2_b0 | λ 1.7017 | −7.416 | −11.000 | ≤ −11.0 | edge | 0.0183 |
| t2_s0_l0_P_d0_b0 | λ 5.6287 | −6.865 | −11.000 | ≤ −11.0 | edge | 0.0073 |
| t2_s0_l0_P_d0_b0 | λ 6.1523 | −7.690 | −10.999 | ≤ −11.0 | edge | 0.9274 |
| t4_s0_l0_P_d1_b0 | λ 0.9163 | −7.178 | −11.000 | ≤ −11.0 | edge | 1.1215 |

## Deviation 5, 17 November 2026: the validation's ground truth for peaked rows

Written after the 11-row validation of deviation 1 failed on a scratch copy of the study's cache
(commit 4b03c4a), and before any new ground truth is computed.

**What failed.** The six named rows passed. The five rows drawn by the seed 20261117 from the
rows with a share of at least 10⁻³ failed in the prior box at every level, for one reason:
their ground truth was not valid by its own test (its spline and exponential-rule values within
τ/2). Its spline − exponential differences were +0.003 to +0.062, against limits of 0.00125
(0.0027 for t4_s0_l0_P_d1_b0 at λ 5.2360). In all five the refined ln Z lay within 0.0004 of
the dense fits' spline value. These rows' posterior widths are 0.013 to 0.020 decade, close to
the truth's 0.01-decade spacing, and the exponential rule, linear in ln L between points,
under-integrates a peak at that spacing.

| No. | Type | Reason | Timing | Original wording | Deviation description | Reader impact |
|---|---|---|---|---|---|---|
| 5 | Analytic approach | Preregistered plan not possible or inappropriate because of data | After results of the preregistered approach were known: the validation's failure on the scratch copy, its per-row values and the rows' widths; before any new ground truth is computed, and before any gate, reading or Bayes factor | "Ground truth by dense real fits, independent of the refinement: a fit every 0.01 decade where ln L at the smallest level lies within 25 of the row's best knot, every 0.05 decade elsewhere in the widest box … It is valid where its spline and exponential values agree within τ_j / 2; a row passes where the refined ln Z lies within τ_j of it at every level in the prior box" (amendment 4, D, kept by deviation 1) | For the five failing drawn rows only, a new ground truth across the prior box. The band is every x in the prior box where ln L ≥ (the row's maximum − 25) at some level, wherever it falls, a secondary region included: each run of calls that qualify, widened to the nearest call outside it on either side, and recomputed from every call until it no longer grows. Inside the band, real fits at spacing w/10 (w the row's posterior width estimate, the grid anchored at the box's lower edge); the w/5 grid is its every other point. Outside it, real fits on the registered 0.05-decade grid, with the row's calls before deviation 1. Two truths: ln Z on the w/10 grid and on the w/5 grid, each with the same points outside the band. The outside region's share of the w/10 truth's integral is measured per row and level (the larger of the spline and the exponential rule); it must be below 10⁻⁶, or the band is widened to ln L 35 and the truth recomputed. Valid where, at every level, the w/5 and w/10 spline values agree within τ_j/4 and the exponential rule lies nearer the spline on the w/10 grid than on the w/5 grid. A row passes where its refined ln Z lies within τ_j of the w/10 spline value at every level. The six named rows keep their result. If any of the five fails against a valid truth, Stage 1 stays paused and is reported, and nothing further is changed. | The validation's instrument is changed after it failed, so a pass under it is weaker evidence than a pass under the registered instrument would have been; the registered results are reported beside the new ones. The refined values are not changed by this deviation. |

Estimated on the rows' calls before the run: bands of 0.40 to 0.68 decade, about 2,000 fits in
all (about 400 a row), one row to a worker. The fits are not added to the rows.

**Code defect fixed before this deviation (ffbb963).** On read-back, of two records at one point
that both converged with contacts the earlier was kept, so a refinement fit landing on a
converged coarse call was dropped and the row's knots, estimate and guard differed from those
the refinement job had used. Now the later is kept, as the job held it. On the scratch copy this
touches 3 rows: t2_s0_l0_J, whose r3_1 estimate falls from 0.0098 to 0.0015 against its
tolerance of 0.0025 (its refinement rerun needed no fit; its ln Z is unchanged to 10⁻⁵);
t4_s0_l0_R at r 0.12 (share 2 × 10⁻³²) and t2_s0_l0_P_d1_b0 at λ 5.6287 (share 6 × 10⁻¹⁵), whose
ln Z move by up to 2.36 and 0.39 in some box. The study's cache has not been recomputed.
