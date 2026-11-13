# Phase 2 design sketch: the fork engine on Mercury's perihelion

Scoping note, 2026-11-12. Not a pre-registration, not a frozen plan. No project code was written
or run. Numbers marked "scratch" come from closed-form formulas evaluated numerically and from
one small three-body integration done outside the repository; they are orders of magnitude and
must be redone inside the project before any use.

Units as in docs/RULES.md rule 5: AU, days, solar masses. Rates are in arcseconds per Julian
century unless stated.

## 1. What Phase 1 already provides

| need | exists in `src/fork` | reuse |
|---|---|---|
| symplectic N-body map | `integrators/wisdom_holman.py`: `step`, `_kick`, `_jump`, `integrate_states` | yes, with an added kick term |
| Kepler drift | `integrators/kepler.py`: `kepler_drift_inplace` | yes |
| reference integrators | `integrators/rk4.py`, `leapfrog.py`; REBOUND in tests only (rule 1) | yes |
| elements and frames | `data/frames.py`: `elements_to_state`, `state_to_elements`, precession | yes |
| JPL states with cache | `data/horizons.py`: `state_vectors_ecl_j2000` | yes; new bodies |
| forward model interface | `model/forward.py`: `ForwardContext`, `predict_states_ecl_j2000`, `predict(theta, t_jd)` | interface yes, kernel new |
| light time, aberration, frame of date | `model/apparent.py`: `ApparentGeometry`, `apparent_ra_dec`; `frames.geocentric_lon_lat_ecl_j2000_light_time` | partly; Sun and contacts are new |
| grouped errors with marginalised offsets | `model/observation.py`: `ErrorModel`, `RaDecObservations.whitened`, `chi2` | pattern yes, class new |
| Student-t errors | `model/student_t.py` | yes |
| Delta T | `data/observations.py`: `delta_t_s` (1600-1860 only) | extend to 1881 |
| priors, unit-cube transform, cells | `inference/priors.py`: `PosteriorSpace`, `UranusBox` | pattern yes, spaces new |
| nested sampling with checkpoints | `inference/nested.py`: `run_nested`, `PartitionedResult` | yes |
| least squares, multistart, profiles | `inference/fit.py`, `inference/profile.py` | yes, after generalising names |
| fast likelihood and reweighting | `model/perturbative.py`: `PerturbativeProblem`; `model/linearized.py`; `inference/reweight.py` | yes; this is the main cost lever |
| double-entry transcription, spot check | `data/transcription.py`, `data/spotcheck.py` | yes |
| error calibration | `data/calibration.py` | pattern yes |
| compute ledger | `fork/ledger.py`: `append`, `read` | yes |

Phase 1's forward kernel is specific to Uranus: `_initial_conditions` builds Sun, Jupiter,
Saturn, Uranus and one perturber, and `model/bodies.py` folds the four inner planets into the
central mass (`M_CENTRAL_MSUN`). Neither survives for Mercury.

## 2. What must be built

### 2.1 Inner-system dynamics

- **Bodies.** Mercury, Venus, Earth-Moon barycentre, Mars, Jupiter, Saturn, and Uranus and
  Neptune for completeness (Le Verrier gives 0.1 for Uranus's share of Mercury's perihelion
  motion; p. 99). The central mass is the Sun alone. `bodies.py` needs the inner masses as
  bodies, and `horizons.py` needs heliocentric states for them.
  `states_about_inner_barycenter_ecl_j2000` no longer applies.
- **Step.** Phase 1 uses dt = 12 d for a 30,700 d orbit. Mercury's period is 88 d with e = 0.206.
  A first guess is dt = 2 d (44 steps per orbit). Rough scaling of the Wisdom-Holman error
  (periodic terms of relative size (k n dt)^2 / 24 on perturbations of a few arcseconds) suggests
  errors of 0.01 to 0.2 arcsec at that step, with the secular rates much better. This is
  UNVERIFIED and is the first thing to measure: perihelion rate against a high-accuracy reference,
  target 0.01 arcsec per century. If it fails: dt = 1 d, or symplectic correctors, or a
  higher-order composition. A different splitting is not needed for close encounters, since
  there are none, unless a discrete intra-Mercurial planet is integrated (section 3.1).
- **Initial conditions and the leak rule.** The analogue of rule 2: Mercury's six elements are
  always free. Its JPL state is never held, because the JPL ephemeris is integrated with general
  relativity and holding it would hand the answer to the modified-law branch. The other planets
  start from JPL states at the epoch and are integrated under the hypothesis being tested.
- **Earth for the geometry.** The observer's position can come from JPL at each epoch, as
  `ApparentGeometry` does now. It must be the geocentre, not the Earth-Moon barycentre: the
  offset, up to 4,670 km, subtends about 10 arcsec at Mercury's transit distance. The Moon is not
  part of any hypothesis, so a fixed JPL offset is legitimate.
- **Planetary masses.** Modern masses are known to parts in 10^8 and are independent of the
  anomaly. Le Verrier's own escape route was the mass of Venus (relation (A), p. 100). A
  "19th-century track" would give Venus's mass a wide prior. That belongs to the mundane-error
  branch.

### 2.2 Transit-contact observation model

New module, for example `model/transit.py` (name not decided):

1. Mercury and the Sun as seen from the geocentre at observation time t: light-time iteration
   for each body, then the angular separation of centres. Aberration is common to both to first
   order; the light-time difference is not. Mercury's light time is about 5 minutes, during which
   it moves about 30 arcsec relative to the Sun. `apparent.py` has the iteration for one body.
2. Contact condition: separation equals R_sun(t) + s R_merc(t), s = -1 for interior contacts and
   +1 for exterior ones. Solve for t by root finding from the tabulated date.
3. Predicted quantity: time of geocentric contact in the table's time scale. Residual in seconds,
   or converted to arcseconds along the track with the local rate.
4. Time scale: Greenwich mean time to TT through Delta T. `delta_t_s` stops at 1860. Mercury
   moves 0.098 arcsec per second relative to the Sun in November transits and 0.064 in May
   (scratch, from orbital speeds; consistent with central durations of about 5.5 h and 8 h).
   A Delta T error of 10 s is therefore 1 arcsec. Early Delta T values are uncertain at about
   that level (UNVERIFIED size). Options: adopt the lunar-derived Delta T as Morrison and Ward
   did; add a free quadratic term; or both. This choice affects the mean motion and, through
   the November-May difference, the perihelion.
5. Nuisance parameters: effective solar semidiameter offset (irradiation, black drop), possibly
   by era; Mercury's semidiameter (Newcomb used 3.34 arcsec at unit distance; 1895, p. 3); a
   contact-phase offset for exterior contacts.
6. Error model: Newcomb's weights as relative variances with one or two free or calibrated
   scales; Student-t as the robust alternative. Whether scales may be calibrated against
   "observed minus JPL", as Phase 1 did for Uranus, is a question for the pre-registration: it
   uses Mercury's true orbit for error sizes only, but the conservative route is scales from the
   data alone, as in Phase 1's S2 run.

Test targets: the contact times of modern transits (for example 2003, 2006, 2016, 2019) from
JPL Horizons, to about 1 s.

### 2.3 Meridian places (set C of the inventory)

Geocentric apparent ecliptic longitude and latitude at a Paris mean time. This is close to the
existing `apparent_ra_dec` path with an ecliptic-of-date output and a phase-defect nuisance term.
It is cheap to add and is expected to carry little weight (Le Verrier, p. 98).

## 3. The competing forward models

All share Mercury's six elements and the nuisance parameters. The interface stays
`predict(theta, t_jd)`; what is predicted is a contact time or an apparent place.

### 3.1 Hidden mass (`model/hidden_mass.py` grows these beside the Neptune case)

**(a) One intra-Mercurial planet.** Parameters: mass, a, e, i, node, perihelion, mean longitude.
Enters as one more body in `_initial_conditions`; the kick and drift are unchanged.
Difficulties: its period is 12 to 45 d, so the step falls to 0.5 d or less, and its mean
longitude is unconstrained by data that cannot see its short-period terms (section 5.5). In
practice it is indistinguishable from a ring unless a is above about 0.2 AU.

**(b) A ring.** Parameters: mass, radius, and the two angles of its plane; optionally a width.
Enters as a position-only acceleration on each planet, added in `_kick`. The potential of a
circular ring needs complete elliptic integrals, or N equal point masses on a fixed circle, or
the multipole series for radius well inside Mercury's perihelion distance (0.307 AU). An
exterior ring between Mercury and Venus (Newcomb 1895, section 59) is the same code with a
different radius prior.

**(c) Solar quadrupole J2.** One parameter if the Sun's spin axis is taken as known
(inclination 7.25 deg to the ecliptic, node about 75.8 deg; UNVERIFIED as sourced values).
Position-only acceleration in `_kick`. It is the small-radius limit of (b) with the plane fixed.

To lowest order in eccentricity and inclination, an interior ring or planet of mass m at radius
a' < a gives

    d(varpi)/dt = (n/4) (m/M_sun) alpha b(alpha),   alpha = a'/a,
    d(Omega)/dt = - d(varpi)/dt   (node measured on the ring's plane),

with b the Laplace coefficient b_{3/2}^{(1)}; for small alpha this is (3/4) n (m/M_sun) alpha^2.
For J2: d(varpi)/dt = (3/2) n J2 (R_sun/a)^2 / (1 - e^2)^2, with the node regressing at the same
rate about the solar equator. Both scale as n / a^2, that is a^(-7/2).

Mass needed for 43 arcsec per century at Mercury (scratch, lowest order):

| ring radius (AU) | mass (M_sun) | mass (Earth masses) | orbital period (d) |
|---|---|---|---|
| 0.10 | 1.4e-6 | 0.47 | 11.6 |
| 0.15 | 5.2e-7 | 0.17 | 21.2 |
| 0.20 | 2.2e-7 | 0.074 | 32.7 |
| 0.25 | 9.2e-8 | 0.031 | 45.7 |

The lowest-order formula overstates the mass as the radius nears Mercury's perihelion: a direct
integration with these masses gave 48, 49, 52 and 62 arcsec per century (scratch). J2 would need
to be 3.4e-4, which is an equatorial excess of the solar radius of 0.49 arcsec. Newcomb's "rather
more than half a second" (1895, p. 111) agrees. The measured J2 of 2.25e-7 gives 0.029.

### 3.2 Modified law (new `model/modified_law.py`)

**(a) First post-Newtonian term.** The test-particle 1PN Hamiltonian per unit mass is

    H = p^2/2 - mu/r - p^4/(8 c^2) - 3 mu p^2/(2 c^2 r) + mu^2/(2 c^2 r^2).

The correction mixes p and r, so it is not a kick (function of position) and not a drift. The
Wisdom-Holman split does not carry over unchanged. Three routes:

1. *Effective potential.* Writing p^2 = 2 (H_K + mu/r), with H_K the Kepler Hamiltonian, the
   correction becomes -(1/c^2) [H_K^2/2 + 4 H_K mu/r + 3 mu^2/r^2]. The first term is a function
   of H_K alone: it commutes with the Kepler drift and only rescales its step. The last is
   position-only: an exact kick with potential -3 mu^2/(c^2 r^2). The middle term, with
   H_K = -mu/(2a), acts as a rescaling of mu by 2 mu/(a c^2), about 5e-8 for Mercury, plus terms
   that do not move the perihelion at this order. Keeping only the kick
   dV = -3 mu^2/(c^2 r^2) reproduces the perihelion advance 6 pi mu / (c^2 a (1 - e^2)) per orbit
   exactly and shifts the mean motion by a constant that the free semi-major axis absorbs. This
   is the route of Nobili and Roxburgh (1986) and of Saha and Tremaine (1994) (both META; details
   not checked against the papers, UNVERIFIED). It stays symplectic and costs one extra term in `_kick`.
2. *Split with a modified drift.* Saha and Tremaine's full treatment keeps the H_K^2 term as a
   rescaled drift. More faithful, more code, no gain at 19th-century precision.
3. *Full Einstein-Infeld-Hoffmann forces with an implicit or non-symplectic step.* Needed only
   as a reference in tests (REBOUNDx, Tamayo et al. 2020, META; test dependency only, rule 1).

Route 1 differs from the full equations by periodic terms of order (v/c)^2 e, about 1e-8 rad or
0.005 arcsec (scratch), and by the planet-planet cross terms that Will (2018) puts at a few
parts in 10^6 of the main effect. Both are far below the data.

**(b) Parametrised perihelion term.** dV = -xi 3 mu^2/(c^2 r^2), with xi free. In the PPN
language xi = (2 + 2 gamma - beta)/3, and xi = 1 is general relativity. This also covers an
inverse-cube force of any origin, which Newcomb rejected on non-orbital grounds (1882, p. 476).
One parameter; the advance scales as n / (a (1 - e^2)), that is a^(-5/2).

**(c) Hall's exponent.** Force proportional to r^-(2 + delta). The difference from Newton is a
position-only potential, so it is again a kick. One parameter; the advance is (delta/2) n for
every planet, scaling as a^(-3/2). Whether the exponent applies to planet-planet forces as well
is a modelling choice; for the perihelia it makes no difference at this precision.

All three modified-law variants are position-only in practice. The honest statement on the
splitting is therefore: true 1PN dynamics is velocity-dependent and breaks the trivial split,
but the part of it that 19th-century data can see is carried by a potential. The general case
is deferred with a test that bounds what is dropped.

### 3.3 Mundane error (`model/mundane.py`, listed for completeness)

Venus's mass with a 19th-century prior; a Delta T polynomial; contact-phase and semidiameter
systematics by era. These are nuisance parameters in the other two branches and the whole model
in this one.

## 4. Likelihood and inference: what is touched

- `model/observation.py`: a sibling of `RaDecObservations` for contact times, with `whitened`
  and `chi2` of the same shape, so `inference/fit.py` and `inference/nested.py` need only a
  generalised problem object in place of `LongitudeObservations` and `ForwardContext`.
- `inference/priors.py`: `PosteriorSpace` is written around Uranus plus one perturber. It needs
  a model-agnostic version: a list of named blocks, each with a transform and a prior volume.
  The Mercury box plays the role of `UranusBox` and is common to all models, so its volume
  cancels in Bayes factors, as the Phase 1 docstring anticipates.
- `inference/nested.py`: unchanged in method. Partitioning by cells may be unnecessary; the
  ring's plane angles are the only place where multimodality is expected.
- `model/perturbative.py`: the first-order response of the target planet to a small extra
  acceleration, integrated along a reference orbit, is exactly the structure needed here. Every
  hypothesis is a perturbation of at most 1e-7 of the solar force. A linear response in the
  amplitude (xi, delta, ring mass, J2) is essentially exact; the ring's radius and plane enter
  non-linearly and need the quadrature per call or a precomputed grid.
- `inference/reweight.py`: importance reweighting of fast-likelihood samples against the full
  integration, as in M8.
- Evidence normalisation: Phase 1 dropped the error-model constant, so ln Z was comparable only
  between models on the same data and error model. That is still the case needed here, but the
  nuisance priors must be identical across branches and written down once.

## 5. Identifiability

**Stated plainly: an inner mass and a 1PN term both advance the perihelion, and the transit
data measure essentially that one number.** What follows quantifies what else differs.

### 5.1 What the transits determine

Newcomb (1882; restated 1895, pp. 61-62) reduces the transits to two functions of the
corrections to the elements, one from November and one from May transits, and a node function:

    V = 1.487 dl - 0.487 d(varpi) - 1.137 de - (Earth terms)
    W = 0.716 dl + 0.284 d(varpi) + 0.896 de - (Earth terms)
    dV/dt = -2.63 +/- 0.59,   dW/dt = +1.84 +/- 0.60,   dN/dt = +0.28 +/- 0.62   (arcsec/cy)

Two secular rates cannot fix three unknowns (mean motion, perihelion rate, eccentricity rate).
Newcomb takes the eccentricity rate from theory and solves for the other two. Propagating his
two errors through his own equations (scratch, errors taken as independent) gives

    perihelion excess = 43.0 +/- 1.3 arcsec per century.

So the transits detect the anomaly at about 33 sigma and measure its size to 3 per cent. The
node function is consistent with zero at 0.6.

### 5.2 Secular advance across planets

Formulas: 1PN, 3 n mu / (c^2 a (1 - e^2)); Hall, (delta/2) n; ring or J2, proportional to
n / a^2. All normalised to Mercury (scratch; Hall column as printed by Newcomb 1895, p. 119).

| planet | 1PN | Hall | ring or J2 | e | 1PN times e | Hall times e | ring times e | 1895 obs. minus theory, times e |
|---|---|---|---|---|---|---|---|---|
| Mercury | 42.98 | 42.34 | 42.98 | 0.2056 | 8.84 | 8.71 | 8.84 | +8.48 +/- 0.43 |
| Venus | 8.62 | 16.58 | 4.4 | 0.0068 | 0.058 | 0.11 | 0.030 | -0.05 +/- 0.25 |
| Earth | 3.84 | 10.20 | 1.4 | 0.0167 | 0.064 | 0.17 | 0.024 | +0.10 +/- 0.13 |
| Mars | 1.35 | 5.42 | 0.33 | 0.0934 | 0.126 | 0.51 | 0.031 | +0.75 +/- 0.35 |

The observable is e times the perihelion rate. For Venus and the Earth every model predicts
less than half of one mean error. For Mars, 1PN and the ring differ by 0.1, against an error of
0.35; Hall differs from 1PN by 0.38, about one mean error. **The scaling law across planets, the
cleanest theoretical discriminator, is invisible in 19th-century data.**

### 5.3 Nodes and inclinations

A ring or J2 makes each orbit's pole precess about the ring's axis at the perihelion rate. A
1PN term or Hall's exponent moves no node (the Lense-Thirring term is of order 0.001 arcsec per
century; UNVERIFIED size). Predicted rates when the hidden mass is tuned to 43 at Mercury
(scratch; quadrupole scaling for Venus; signs in the fixed ecliptic frame, not reduced to
Newcomb's conventions):

| hidden-mass plane | Mercury sin i dOmega/dt | Mercury di/dt | Venus sin i dOmega/dt | Venus di/dt |
|---|---|---|---|---|
| ecliptic | -5.2 | 0.0 | -0.26 | 0.0 |
| solar equator (J2) | -0.4 | -2.5 | +0.30 | 0.0 |
| Mercury's orbit plane | 0.0 | 0.0 | +0.21 | +0.26 |
| Venus's orbit plane | -3.0 | -1.2 | 0.0 | 0.0 |
| 1895 obs. minus theory | +0.61 +/- 0.52 | +0.38 +/- 0.80 | +0.60 +/- 0.17 | +0.38 +/- 0.33 |

Mutual inclinations used: Mercury to Venus 4.3 deg, Mercury to the solar equator 3.4 deg,
Mercury to the ecliptic 7.0 deg.

Reading:
- A ring in the ecliptic is excluded by Mercury's node at about 10 mean errors. Newcomb said
  so in 1882.
- A ring or planet in Mercury's own orbital plane leaves Mercury's node alone by construction.
  It then shows only in Venus's node and inclination at 0.2 to 0.3, one to two of Newcomb's mean
  errors.
- The solar-equator case is disfavoured by Mercury's inclination rate at about 3.5 mean errors.

### 5.4 Illustration on Newcomb's 1895 discrepancies

Sum of squared normalised differences over the eight quantities of sections 5.2 and 5.3
(scratch; errors taken as independent, which Newcomb warns they are not; each one-parameter
family fitted in amplitude; this is arithmetic on published numbers, not a run of the engine):

| model | free parameters | sum of squares |
|---|---|---|
| no anomaly | 0 | 410 |
| 1PN at its predicted value | 0 | 19.5 |
| law family, xi free | 1 | 18.9 |
| Hall's exponent | 1 | 16.6 |
| ring in the ecliptic | 1 | 122 |
| J2 or ring in the solar equator | 1 | 24.9 |
| ring in Mercury's orbit plane | 1 (plus the choice of plane) | 11.9 |

The best hidden-mass model fits Newcomb's numbers better than general relativity by about 7,
almost all of it from the Venus node (+0.60 +/- 0.17), with the Mars perihelion pulling toward
Hall. This is what Newcomb himself reported: the ring between Mercury and Venus fits best, and
he preferred Hall for other reasons. The Venus node excess is not real (UNVERIFIED as a sourced
statement). **A blind engine given the 1895 digest would probably not choose the modified law.**

### 5.5 Short-period signature of a discrete planet

One coplanar circular planet with the masses of section 3.1, Mercury on its real ellipse, 40
years, secular trends in e, perihelion and mean longitude removed (scratch, crude):

| planet a (AU) | periodic residual in Mercury's heliocentric longitude, rms (arcsec) |
|---|---|
| 0.10 | 0.05 |
| 0.15 | 0.16 |
| 0.20 | 1.0 |
| 0.25 | 2.6 |

Seen from the Earth at a transit these shrink by a factor of 0.5 to 0.8. One transit is good
to roughly 0.3 to 1 arcsec (Newcomb's V and W at epoch have errors of 0.31 and 0.25; Morrison
and Ward find a scatter of 11 s or more per timing, about 1 arcsec). About 25 transits sample
the perturbation at effectively random phases. So a planet inside 0.15 AU is undetectable by
its periodic terms; one at 0.2 to 0.25 AU would add a scatter comparable to the noise. It would
show as excess variance, not as a recoverable orbit. This is the opposite of the Neptune case,
where the periodic term was the whole signal.

### 5.6 Verdict

From 19th-century Mercury data alone the two branches are **not separable in the likelihood**:

1. The transits measure one combination well (the perihelion excess) and the node rate to 0.6.
2. The node rate removes hidden mass in or near the ecliptic. It cannot remove hidden mass in
   Mercury's orbital plane.
3. Every other difference (scaling with distance, Venus's node, Mars's perihelion, short-period
   terms) is at or below one mean error of the best 1895 determinations.

A Bayes factor from these data is therefore set by the priors: how much volume the hidden-mass
model wastes on planes and radii that the node excludes, against how the modified-law prior on
xi or delta is drawn. That is an Occam-factor result. It can be computed and reported honestly,
with its prior sensitivity, but it is not a discovery of the law from the data.

General relativity has one extra advantage that is not in the likelihood shape: it predicts the
size with no free parameter. Model (a) of section 3.2 with xi fixed at 1 has zero free
parameters against one or more for every rival. Using that requires knowing the theory, which
a blind engine of 1882 does not. The engine can still report the evidence for "xi = 1 exactly"
as a separate, clearly labelled line.

### 5.7 What additional data would separate them

| data | what it adds | era | cost |
|---|---|---|---|
| Venus: meridian places and the transits of 1761, 1769, 1874, 1882 | Venus's node to better than 0.17; a ring in Mercury's plane predicts 0.2 to 0.3 | 19th century | very high (set H scale); and the 1895 value was biased |
| Mars: perihelion | Hall against 1PN at about one mean error in 1895 | 19th century | very high |
| Mercury transits 1891-1973 | perihelion to 0.5 (Morrison and Ward); still one number | to 1973 | moderate if their collection is public |
| eclipse and transit searches for the matter itself | bounds on a discrete body or on scattered light | 1860-1908 | not orbital; enters as a prior |
| radar ranging to Mercury, Venus and Mars | the distance scaling; (2 + 2 gamma - beta)/3 to 0.7 per cent in 1972 | from the 1960s | available in ephemeris data sets; outside the historical brief |
| spacecraft ranging | the periodic 1PN signature; J2 separately (Park et al. 2017) | from 2011 | outside the brief |

Within the 19th century the only handle with real leverage is Venus's node, and the historical
determination of it was wrong in the direction of matter. The data that actually settle the
question are from the 20th century.

## 6. Cost model

From `docs/reproduce.md` and `docs/notes/M8.md`:

| quantity | Phase 1 | Phase 2 estimate |
|---|---|---|
| planets integrated | 3 or 4 | 8 |
| span and step | 156 yr at 12 d: 4,750 steps | 204 yr at 2 d: 37,250 steps |
| cost per step | 1.2 microseconds (5.52 ms per call) | about 3.5 microseconds (2 to 5; scales between bodies and pairs) |
| full likelihood, per call | 5.5 ms | 0.13 s (0.08 to 0.2); 0.26 s at dt = 1 d |
| fast likelihood, per call | 0.21 to 0.28 ms | 0.3 to 1 ms assumed |
| calls per nested run | 4.3 million per cell (17.3 million for 4 cells) | 4 million assumed (10 to 14 dimensions, one mode expected) |
| one nested run, full likelihood | 28.5 CPU-h (4 cells) | about 145 CPU-h (90 to 290) |
| one nested run, fast likelihood | 0.13 CPU-h per run | about 0.6 CPU-h |

A full-likelihood nested run costs about five Phase 1 headline fits. The fast likelihood is
therefore not an optimisation but a precondition, and its validation is a milestone of its own.

## 7. Open questions

1. Is the 1895 digest (section 5.4) in or out? Using it is cheap and it spans four planets, but
   it is not blind and it was looked at for this note.
2. May error scales be calibrated against JPL Mercury, as Phase 1 did against JPL Uranus?
3. Is the lunar-derived Delta T an allowed input?
4. Does "modified law" include the fixed-parameter relativistic prediction as its own model?
5. Is the Uranus data to be run through the same fork, with a modified-law branch added, so that
   the engine is shown on one case of each kind?

## 8. Two candidate questions (scoping, added 12 Nov 2026)

Neither is a plan. Nothing here has been run, and neither may be run on real data before a
plan is frozen. Both sources they lean on are already in `prior_art.md` (Newcomb 1895,
section 3; Seeliger 1906, section 4, read from the Academy's scan on 13 Nov 2026).

### 8.1 The Occam-factor question

Compare the relativistic term with no free parameter against a ring or dust model with
stated priors on mass and radius, and report how the Bayes factor depends on the prior
volume instead of reporting one number.

- Models: 1PN with its strength fixed (zero parameters); a circular ring (mass, radius);
  optionally Seeliger's form, two ellipsoids and a rotation of the reference system, as the
  flexible end.
- What is expected before any run: Mercury's advance fixes one combination of ring mass and
  radius to about 5 % (8.48 ± 0.43 in section 5.2's units). The evidence for the ring then
  falls roughly in proportion to the prior width of that combination; directions the data
  do not constrain cost nothing. A log-uniform mass prior over D decades costs a factor
  near D times 2.3 divided by 0.05 against a model that predicts the value outright.
- The deliverable would be a curve: log Bayes factor against prior width, for two or three
  prior shapes, with the widths at which the verdict changes sign or grade stated.
- What must be stated first: the prior families and their ranges, with a reason from the
  sources of the time (Le Verrier's and Newcomb's own mass limits, the visibility limits of
  section 5.5), and whether the 1PN strength is fixed or carries the PPN factor.
- Risk: the answer is known in outline (section 9 of `prior_art.md`). The new part would be
  the dependence on the prior, on period data, with the mundane-error model in the same
  comparison.

### 8.2 The cross-planet signature

1PN fixes the advance of every planet once Mercury's is given: it goes as n / (a (1 − e²)).
A ring's goes as n / a² outside the ring, and also moves the nodes. The question is whether
any period data can tell the two dependences apart.

- What in Newcomb 1895 would test it: the table of observed minus theoretical secular
  variations (pp. 109-110): e dϖ/dt for Mercury, Venus, the Earth and Mars with mean
  errors, and di/dt and sin i dΩ/dt for Mercury and Venus. Behind it, the 62,030 meridian
  observations of chapter I and the mass adjustments of chapters V and VI.
- What section 5.2 already shows: with those mean errors the perihelia of Venus and the
  Earth cannot separate the models (every prediction is under half a mean error), and for
  Mars 1PN and a ring differ by 0.1 against 0.35. The perihelia alone are underpowered.
- Where the leverage may be instead: the nodes. A ring in or near the ecliptic moves the
  nodes of Mercury and Venus; 1PN does not. Newcomb's sin i dΩ/dt for Venus is
  +0.60 ± 0.17. Whether that is signal or his mass of the Earth is the first thing to
  establish (section 5.3).
- What would be reported: for each model the predicted vector of the nine quantities, the
  covariance Newcomb gives or implies, and the expected log Bayes factor from synthetic
  data at his errors, before any real number is used.
- Constraint: the 1895 digest has been read (section 5.4), so a test on it is not blind.
  It can calibrate the method; it cannot be the result. See `blinding.md` on the data
  branch.
