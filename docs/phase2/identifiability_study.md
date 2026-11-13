# Mercury: a synthetic identifiability study (pre-registration)

Written 2026-11-13, before any code of this plan exists and before any synthetic data set
of it was drawn. Everything above the line at the end is frozen by the commit that adds
this file; amendments and results go below it, dated.

**No observation of Mercury is fitted here or anywhere.** Newcomb's contact table is used
for two things only: when the contacts happened and how he weighted them. The "data" of
this study are made by the forward model.

## 1 Question

Suppose Mercury had been observed when and as well as Newcomb's transit contacts were, and
suppose the excess advance of its perihelion had one of four causes. Could a comparison of
models by Bayesian evidence tell which? How much of any preference is the data's, and how
much is the width of the prior on an unseen mass?

## 2 Inputs

| input | SHA-256 |
|---|---|
| `data/historical/mercury/newcomb_1882_contacts.csv` (epochs and weights only) | 70d4d0d1ea242c55c227d593b1cf9238593ce52c3cd39ae0ce9c1fdbd41712f9 |
| `tests/fixtures/phase2/v4_reference.npz` (planetary states at 1900.0) | f8032fc2f622112303e9e3912fc32443b54e8e206251bb91babd7c3d29f24f77 |
| `src/fork/model/inner_system.py` | a132da7aab691230d988ee3c2d42b2ec16655af7e65159856ff480425e6f512b |
| `src/fork/phase2/secular.py` | 55d0eb506eaff1899731d119324474287bcb6091598e6dcd72e98213afab880b |
| `src/fork/phase2/identifiability.py` | 14676bb8d02a945dae5ed66f8c75c13084857b4e899357cc0b5faf9fa3c4d3f0 |

- The physics (1PN term, solar J2, ring, interior planet) is the code validated in
  `validation.md` (V1 to V5). Nothing is added to it.
- From the contact table, per row: the date and time (to place the epoch), the weight,
  and which contact it is. The k terms, the second value of the braced row and every flag
  are not read. The residual tables of pp. 453 to 454 are not in the repository and are
  not read.

## 3 The simplified observation, stated plainly

A contact time measures where Mercury is along its track at a known instant. The contact
model itself (limbs, semidiameters, light time, ΔT) is not built. I replace each contact
by **Mercury's heliocentric ecliptic longitude at the contact's epoch**, with a Gaussian
error σ_i = ρ_i σ₁ / √w_i:

- w_i is Newcomb's weight;
- ρ_i is the rate at which Mercury's heliocentric longitude gains on the Earth's at that
  epoch (from the Newtonian model; about 0.2″ per second in November and 0.08″ in May);
- σ₁, the timing error of a unit-weight contact, is not read from the source (the table
  gives relative weights only). It is fixed by calibration before any data are drawn: the
  value at which this observation, with Mercury's six elements free, measures the strength
  of the 1PN term to 3.0 %, the precision Newcomb's own errors give for the advance (43.0
  ± 1.3″ a century, `design_sketch.md` § 5.1). Sensitivity levels: σ₁/2 and 2σ₁.
- Epochs: the printed Greenwich mean time on the astronomical day (from noon), taken as
  TDB without ΔT. For synthetic data a common offset of the time scale changes nothing.

## 4 Truths

Four worlds, each integrated in full (no linearisation) from the same state at 1900.0:

| truth | cause | value |
|---|---|---|
| a | the 1PN term | strength 1 |
| b | a circular ring in the ecliptic at 0.2 AU | mass tuned to 42.98″ a century |
| c | a point planet on a circular orbit in the ecliptic at 0.2 AU, longitude 0 at epoch | mass tuned likewise |
| d | a large solar J2 | J2 tuned likewise (about 3 × 10⁻⁴) |

Tuning is on the integrated rate of the longitude of perihelion against the Newtonian
run, 1677 to 1881, by one linear rescaling of the analytic value, then checked (V1).

Per truth: one noise-free data set and eight with noise (seeds 1 to 8, rng `[2026, 21,
1882, truth, seed]`, truth = 1..4), at each of the three error levels.

## 5 Models fitted

Mercury's six elements are free in every model with flat priors, and are marginalised
analytically (the longitudes are linear in small element changes). They are the same for
every model, so they cancel in every Bayes factor.

| model | free parameters beyond the elements | prior |
|---|---|---|
| N, Newtonian | none | |
| G, 1PN | none (strength 1) | |
| J, solar J2 | J2 | log-uniform, 10⁻⁷ to 10⁻² |
| R, ring in the ecliptic | mass; radius | mass log-uniform, 10⁻⁹ to 10⁻⁵ M☉; radius uniform, 0.10 to 0.30 AU |
| P, planet in the ecliptic, circular | mass; longitude at epoch; distance | mass as the ring's; longitude uniform on the circle; distance one of 0.15, 0.20, 0.25 AU with equal weight |

- The planet's distance is not sampled continuously: over two centuries its phase at each
  contact turns many times for a small change of distance, and the likelihood is a comb.
  Three fixed distances, one of them the truth's, are a stated simplification in the
  planet model's favour when the truth is a planet.
- **Likelihood.** Each model's effect on the 58 longitudes is taken as linear in its
  strength (mass or J2), with the shape computed by full integrations: on 11 ring radii
  (cubic interpolation), and on 48 planet longitudes at each distance (trigonometric
  interpolation). Linearity and interpolation are tested (V2, V4).
- **Evidence.** dynesty, static, nlive 1000, dlogz 0.01, uniform sampling, over the
  parameters of the table (1 to 2 dimensions; the planet's three distances are three runs
  averaged). rng `[2026, 21, 1882, truth, seed, level, model]`. Models N and G have no
  parameter: their evidence is the marginal likelihood itself.

## 6 Validation, before anything is read

| gate | pass condition |
|---|---|
| V1 rates | each tuned truth's integrated perihelion rate is within 0.05″ of 42.98″ a century; truth a's own rate is within 0.05″ of 42.98″ |
| V2 linearity | at 30 posterior draws of each fit to a noise-free data set (models J, R, P; 12 fits), χ² recomputed by a full integration at the draw (parameters and conditional best elements) differs from the linear χ² by at most 0.5, for 95 % of the draws |
| V3 sampler | ln Z by nested sampling within 0.1 of direct quadrature on a grid, for every model with parameters on the four noise-free data sets |
| V4 closure | fitting the truth's own model to its noise-free data gives a minimum χ² below 0.5 |
| V5 no observed value | a test shows the synthetic data depend only on the epochs, the weights, the truth and the seed |

If a gate fails, the study stops there; the failure is the result.

## 7 Readings, fixed now

With L = ln B between two models. Grades: |L| < 1 none; 1 to 3 weak; 3 to 5 strong; 5 or
more decisive.

- **R1, is the truth recovered?** For each truth at the calibrated σ₁: "recovered" if the
  truth's model has the largest evidence and leads every other model by L ≥ 3, in the
  noise-free set and in at least 6 of the 8 noisy sets. "Preferred, not decided" if it is
  largest in at least 6 of 8 but leads by less than 3 somewhere. Otherwise "not
  recovered", with the model that wins instead. For truth a, model G against each of J,
  R, P is also given by grade.
- **R2, detection.** ln B(best of G, J, R, P : N) per truth: reported, expected decisive.
- **R3, the prior on the unseen mass.** ln B(R : G), ln B(P : G) and ln B(J : G) on each
  truth's noise-free data, with the mass (or J2) prior widened and narrowed about the
  centre of its range: 1, 2, 4, 6 and 8 decades (by quadrature). Reported as curves, with
  the widths at which the grade of R1 changes. Stated as: "of the preference for G over R
  on truth a, x comes from the data and y from the prior width", x being L at 1 decade.
- **R4, the error level.** R1 at σ₁/2 and 2σ₁. Reported as a table; a change of verdict
  is stated.

## 8 What this cannot show

- **Nothing about the real Mercury.** It says what data of this kind could decide, under
  the simplification of § 3.
- The longitude-only observation leaves out the latitude at transit, which is where a
  ring or planet in the ecliptic differs most from the 1PN term (it moves the node). A
  real contact model may separate the models better than this study finds.
- The fits linearise about the true elements; real data would not grant that.
- No error of the time scale, no systematic error by era or observer, no outliers.
- The ring and the planet lie in the ecliptic at fixed or coarsely gridded distances; an
  inclined or eccentric body is not tried. Hall's exponent and Seeliger's ellipsoids are
  not among the models.
- The calibration of σ₁ rests on a figure derived from Newcomb's summary, not on his
  stated error of a contact.

## 9 The cross-planet extension: scoping only

Not run. Synthetic perihelion data for Venus and the Earth would need those planets as
fitted bodies under each model, and that path has no validation: V1 to V5 test Mercury.
What is already on record (`design_sketch.md` § 5.2): with Newcomb's 1895 mean errors
every model predicts under half a mean error for Venus and the Earth. A scoping paragraph
is added with the results; no number is computed for it.

## 10 Cost and layout

- About 220 integrations of 204 years at under 1 s for the shapes, 500 more for V2, and
  540 one- or two-dimensional nested runs (4 truths, 9 data sets, 3 error levels, 5 runs)
  on a likelihood of microseconds: about 6 CPU-h, mostly sampler overhead. Run detached,
  12 processes, after the Neptune runs free the workers.
- Code `src/fork/phase2/identifiability_study.py`, tests
  `tests/test_phase2_identifiability_study.py`, branch `phase2-identifiability`. Cache
  `data/cache/phase2/identifiability_study/`. Results `docs/phase2/identifiability_study.json`.

---
