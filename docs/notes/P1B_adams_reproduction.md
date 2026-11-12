# Phase 1b, step 3b: Adams's own solution reproduced from his printed equations

Pre-registration. Drafted 2026-11-11, reviewed the same day and amended on review in six points
(§ 8) before freezing; the Q3 table of § 5 is shown to the user once more before the plan is
committed. On commit the plan block is frozen; the transcription record and the results go below
lines, as in `P1B_adams_fit.md`. Stage A (transcription) begins on the approval of 2026-11-11.

## Pre-registration (before any code or run)

### 0 Prior work

- **Sampson 1904**, Memoirs RAS 54, 143 ff. (in the cache, `ads_1904MmRAS..54..143S.pdf`):
  describes Adams's manuscripts and tabulates his determinations; it is the primary catalogue of
  the solutions (step 1) and is read again here for anything it says about his computation.
- **Smart 1947**, W. M. Smart, *John Couch Adams and the discovery of Neptune*, Occasional Notes
  of the RAS 2 (no. 11): the fullest account of his work from the manuscripts. Not yet obtained.
- **Brookes 1970**, C. J. Brookes, *On the prediction of Neptune*, Celestial Mechanics 3, 67–80:
  obtained from ADS and read in full on 2026-11-11 (hash and a summary in `MANIFEST.md`). It
  already does part of this step: from the memoir's equations he re-derives Hypothesis I as
  e′ 0.1509, m′ 0.8267, ϖ′ 311° 35′, mean longitude 323° 39′ at 1846.762 (Adams: 0.1610, 0.8282,
  315° 57′, 325° 7′) and the true longitude at discovery as 328° 4′ against Adams's 328° 42′,
  attributing the difference "mainly to a small numerical error occurring in one of the equations
  derived by Adams" (his 1969 thesis; the paper does not say which equation); Hypothesis II by the
  first hypothesis's procedure comes out 0.1387, 0.7393, 303° 52′, 324° 47′; he then varies the
  number of equations, the epoch and α (sin 29° … sin 37°, with F₄, G and H terms added from
  sin 34° on, and a′ = 30 AU excluded for the near-commensurability), finding the true longitude at
  discovery insensitive to a′ between 30 and 40 AU at any epoch after 1830.

Smart was sought once and not obtained (`MANIFEST.md`). **Consequences for this step, fixed before
freezing:** (i) Q1 has two targets, Adams's print and Brookes's re-derivation; the expected Q1
result is that C2 lands near Brookes and not near Adams, and the step's Q1 contribution is then to
*locate and quantify* the error Brookes reports but does not place (stage B's equation-by-equation
checks do exactly that); (ii) D2's point solutions at α = sin 29°, 30°, 31°, 32° are compared with
Brookes's Table X row (λ′ 323°.8, 327°.1, 330°.4, 333°.2), and the χ² profile and the pole are the
part he did not compute; (iii) the data–method split of Q3 is not in Brookes. **The novelty claim
of this step is limited to (i)'s location of the error, (ii)'s profile and (iii).** Sampson's
manuscript description is read again for the equation Brookes may have had in view.

### 1 Question and why it is the next step

Step 3 (`P1B_adams_fit.md`, results of 2026-11-11) ended in interpretation row 3: with data set D
(Table 77 cut to Adams's series) the planet is located when a is free, but at his distances the
posterior puts it short of Neptune, while his own solutions S4, S5, S6 lie beyond it. The
candidates for that gap are his solution method and the five places where this reconstruction does
not reproduce his data handling (A1 ancient weighting, A2 the 1831 normal place and the Paris
stand-in, A3 normal places and opposition selection, A4 his reductions, A5 Bessel 1823). Step 3
could not tell them apart because it never ran his method, and never saw his data.

**One comparison quantity** is used in Q2, Q3 and D4 and in every threshold: the planet's
heliocentric true longitude of date at 1846.762 (Adams's "October 6, 1846"; taken as
1846-10-06T00:00 TDB). On it, step 3's numbers at 1 October 1846 are: S5 328.75°, H1(G) median
316.50° [312.17, 320.12], Neptune 327.01°. The **step-3 gap** (S5 − H1(G) median) is 12.25° and the
**step-3 shortfall** (H1(G) median − Neptune) is −10.51°; the five days to 1846.762 move each by
less than 0.05°, so the thresholds below are set from 12.2°: **3° (a quarter) and 6° (a half)**. D4
recomputes both at 1846.762 from the cached samples; the thresholds are not revised by it. (Q1,
the reproduction of his print, is read on the quantities he printed, the mean longitude among them.)

Both his data and his method are in print. The memoir (*An explanation of the observed
irregularities in the motion of Uranus*, Nautical Almanac appendix, 1846; pamphlet in the cache,
printed pp. 1–31 = PDF 7–37; Memoirs RAS 16, 427 ff. is the same text) gives the 30 normal-place
residuals he fitted (§ 8 in true longitude, § 9 in mean longitude), his linearised model
(§§ 10–13), the normal equations and every elimination step (§§ 14–27), the final four equations
and their solution (§§ 28–31) for Hypothesis I, the same chain for Hypothesis II (§§ 32–48), and the
back-substituted values and residuals (§§ 49–53, 55). Step 3b runs that chain again and checks it
against his printed answers, then uses the same machinery to separate three things the step-3 gap
conflates: his *procedure* (how he solved his equations), his *model* (first-order perturbations of
the mean longitude at an assumed distance), and his *data* (normal places from Bouvard's residuals,
against set D).

Three questions, in this order:

- **Q1 Reproduction.** Do his printed inputs and his printed method give his printed Hypothesis I
  and II elements and longitudes, to the precision he states?
- **Q2 Procedure and distance.** Does a full least-squares solution of his own equations give the
  same answer as his elimination-and-summation procedure (§§ 21, 25), and do his data, in his
  model, prefer his assumed distance?
- **Q3 Data against method.** With one known truth (Neptune) generating three kinds of data, where
  does each method land at his distance: near the truth, where his solutions landed, or short of
  it, where step 3's fixed-distance posterior landed?

### 2 What Adams printed (inputs; page = pamphlet printed page, PDF = page + 6)

| item | section, page | content |
|---|---|---|
| I1 | § 8, p. 6 | Observation − theory in heliocentric true longitude, Bouvard's tables with his § 7 corrections: 9 ancient (1690, 1712, 1715, 1750, 1753, 1756, 1764, 1769, 1771) and 21 modern normal places (1780 to 1840 every third year), to 0.1″ and 0.01″. |
| I2 | § 9, p. 7 | The same converted to corrections of the mean longitude by the factor r²/(ab) (printed on the page; the conversion checks: 1840, −66.64 → −73.09 is the aphelion factor 1.0968). "These numbers form the basis of the subsequent investigations." |
| I3 | § 10, pp. 7–8 | Correction of the mean longitude from corrections of Uranus's elements: δε + t δn + δx₁ cos nt + δy₁ sin nt + δx₂ cos 2nt + δy₂ sin 2nt, with δx₂ = ¼e{cos(ε−ϖ) δx₁ + sin(ε−ϖ) δy₁}, δy₂ = −¼e{sin(ε−ϖ) δx₁ − cos(ε−ϖ) δy₁}; § 22 gives the numbers (e = 0.046679, ε − ϖ = 50° 15′.8: δx₂ = 0.007460 δx₁ + 0.008974 δy₁, δy₂ = −0.008974 δx₁ + 0.007460 δy₁). |
| I4 | § 11, p. 8 | The perturbation of the mean longitude in Pontécoulant's form, F_i, G_i, H_i as printed functions of z = i(n − n′), n, n′ and the Laplace coefficients A_i with α dA_i/dα, α² d²A_i/dα². |
| I5 | § 12, pp. 8–9 | For α = sin 30° = 0.5: log b₀..b₃, log α db/dα, log α² d²b/dα² (12 logarithms), and the eleven series coefficients (−36″.99 m′ sin(nt − n′t + ε − ε′) … −1″.08 m′e′ sin{3(…) − ϖ + ϖ′}), m′ in units of 1/5000 of the Sun's mass, e′ in units of 1/20. Then the h₁..₃, k₁..₃, p₁..₃, q₁..₃ form. |
| I6 | § 13, p. 9 | Epoch 1810.328 (mean opposition 1810); unit of time 3 synodic periods = 3.0362 yr; anomaly rate 13° 0′.5, n = 13° 0′.6, n′ = 4° 36′.0; the modern equation of condition with its twelve printed arguments; t = −10 … +10. |
| I7 | §§ 14–20, pp. 10–12 | The rule for the normal equations (§ 14, with closed-form sums); the first-group equations (ε), (x), (h₁), (h₂), (h₃), (p₂), (p₃) (§ 15), the second-group (n), (y), (k₁), (k₂), (k₃), (q₂), (q₃) (§ 18), and the two elimination stages (§§ 16–17, 19–20): about 300 coefficients to 4 decimals. |
| I8 | §§ 21–22, p. 13 | The summed equations of each group and, after eliminating δx₂, δy₂ by I3, equations (1) and (2). |
| I9 | § 23, p. 13 | The 9 ancient equations of condition with 16 printed coefficients each (the δn coefficient is the epoch in units: −39.31 for 1690, −32.30 for 1712, … −12.64 for 1771). |
| I10 | §§ 24–26, pp. 14–16 | The ancient equations after eliminating δε, δn, δx₁, δy₁; the multipliers 0.8, −0.6, 1.0, 1.0, 0.9, 0.6, 0.4, 0.3 and 1.0, 1.0, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1 (1690 set aside); equations (3) and (4). |
| I11 | §§ 27–29, pp. 16–17 | (2), (3), (4) with their left sides eliminated by (1); the expressions of h, k, p, q in θ = ε − ε′, β = ε − ϖ, β′ = ε − ϖ′, m′, e′ (§ 28); the three equations with logarithmic coefficients (§ 29). |
| I12 | §§ 30–31, p. 18 | The solution: θ = −51° 30′, p₃/m′ = 271″.57, q₃/m′ = −207″.24; ε = 217° 55′, ε′(1810.328) = 269° 25′; motion in 36 synodic periods 55° 12′ + precession 30′ → 325° 7′ at 1846.762 (6 Oct 1846); e′ = 3.2206 (×1/20 = 0.16103), β′ = 262° 28′, ϖ′ = 315° 27′ → 315° 57′ in 1846; m′ = 0.82816 → mass 0.0001656. |
| I13 | §§ 32–48, pp. 18–25 | Hypothesis II, α = sin 31° = 0.515: the 12 logarithms and eleven coefficients (§ 33), n′ = 4° 48′.5 and the arguments (§ 34), the normal equations (§§ 35–36; (p₂), (p₃), (q₂), (q₃) not formed, § 37), the eliminations (§§ 38–39, 42–43), the ancient equations (§ 40), the multipliers (§ 41), (3) and (4) (§ 43), the final equations (§§ 44–46), the solution θ = −46° 55′, 138″.92, −109″.83; ε′ = 264° 50′, 57° 42′ + 30′ → 323° 2′; e′ = 0.12062, β′ = 279° 14′, ϖ′ = 298° 41′ → 299° 11′; m′ = 0.75017 → 0.00015003 (§§ 47–48). |
| I14 | §§ 49–53, 55, pp. 26–29 | h, k, p, q per hypothesis (§ 49); δε, δn, δx₁, δy₁, δx₂, δy₂ and the element corrections (§ 50); the theoretical corrections split into elements and perturbation per normal place (§ 51); the residuals O − C per normal place and hypothesis (§ 52: modern within ±2.35″, ancient within ±12.8″); the 1843–1845 errors +7.11, +8.79, +12.40 and +5.77, +7.05, +10.18 (§ 53); Flamsteed 1690: 44″.5 and 50″.0 (§ 55). |

Everything is marked P (pamphlet image) once transcribed. The 1876 French translation (Journal de
mathématiques, 3e série, t. II, in the cache) is a second printing of the same numbers; it may decide
a disagreement between the two passes, the pamphlet image always wins over it.

### 3 Transcription (stage A)

Two passes as in steps 1 and 2: pass 1 from the renders `pamphlet_p012.png`–`pamphlet_p037.png`
(printed pp. 6–31); pass 2 by a subagent that sees only the renders and the column specification
below; every numeric token compared; disagreements decided on the named image in a reconciliation
file. About 600 tokens. Files under `data/historical/adams_challis/`, values kept as the printed
digit strings (no rounding, no conversion), angles as printed (`d m.m`):

- `adams_memoir_normal_places.csv`: one row per normal place (30): `year, kind, t_h1, t_h2`
  (modern: the integer t; ancient: the δn coefficient of §§ 23 and 40), `oc_true` (§ 8), `oc_mean`
  (§ 9), `elements_h1, perturbation_h1, total_h1, residual_h1` and the same `_h2` (§§ 51–52;
  blank for 1690, which has no row in §§ 51–52).
- `adams_memoir_series.csv`: long form `hypothesis, section, item, printed`: the 12 logarithms,
  the eleven coefficients (§§ 12, 33), the constants of §§ 13, 34 (epoch, unit, rates, arguments),
  the § 22 relations, and the numbers in the § 28 and § 45 expressions.
- `adams_memoir_equations.csv`: long form `hypothesis, section, equation, term, printed` for
  every printed equation of §§ 15–27, 29, 35–46; term names `lhs, de, dn, dx1, dy1, dx2, dy2,
  h1..h3, k1..k3, p1..p3, q1..q3`, and for the logarithmic equations of §§ 29, 46 `sin1, cos1,
  sin2, cos2, sin3, cos3, P3, Q3, P1c-Q1s, P1s+Q1c, P2c-Q2s, P2s+Q2c` with the sign and the
  bracketed logarithm kept as printed.
- `adams_memoir_solution.csv`: long form `hypothesis, section, item, printed` for §§ 30–31,
  47–50, 53, 55.
- `adams_memoir_statements.csv` and `adams_memoir_reconciliation.md`, as before.

### 4 What is computed

Code in a new subpackage `src/fork/historical/` (`adams_1846.py`, pure numpy, no numba; tests in
`tests/test_adams_1846.py`; runner `src/fork/export/p1b_adams_repro.py` writing
`docs/notes/P1B_adams_reproduction.{json,png}`). Layout approved on review.

**Stage B, his intermediate numbers from his inputs** (each a test with the tolerance stated):

- B1 Laplace coefficients: b^(i)_{1/2}(α) and their first two α-derivatives by quadrature, against the
  12 printed logarithms per hypothesis: |Δ log| ≤ 3 × 10⁻⁵ (his fifth decimal). If Pontécoulant's
  b differs from the modern b^(i)_{1/2} by a constant factor, the factor is found from b₀ and stated.
- B2 The eleven series coefficients from I4 with B1's values: |Δ| ≤ 0.02″ (two units of his last
  digit); the derived h, k, p, q forms of § 28 and § 45 reproduced symbolically from them.
- B3 The 21 modern equations of condition from I6 (coefficients are cosines and sines of the printed
  arguments) and the normal equations by his rule (§ 14): coefficients |Δ| ≤ 0.0005, right-hand
  sides |Δ| ≤ 0.05″ against §§ 15, 18, 35, 36.
- B4 The 9 ancient equations of condition from the epochs in I9: coefficients |Δ| ≤ 0.0005.
- B5 His eliminations (§§ 16–17, 19–20, 24, 38–39, 42) rerun from his printed normal equations:
  |Δ| ≤ 0.003 (round-off accumulates through two stages); the summed and weighted equations
  (§§ 21, 25, 41) and (1)–(4): |Δ| ≤ 0.003 on coefficients, 0.3″ on right-hand sides.

**Stage C, the solution** (both hypotheses, three entry points, so that a failure is located):

- C0 from his printed final equations ((1)–(4) and §§ 29/46): his iteration (§ 30) and a direct root
  find; this isolates the nonlinear solve.
- C1 from his printed normal equations through our rerun of his eliminations.
- C2 from his § 9 numbers alone through our equations of condition, normal equations, eliminations
  and solve: the full chain from data to elements.

Each entry point yields θ, m′, e′, β′ and from them ε′(1810.328), the mean longitude at 1846.762
(ε′ + 12 n′ + 30′, exactly his arithmetic), ϖ′ (1810 and 1846), the mass; then δε, δn, δx₁, δy₁
(§ 50), h, k, p, q (§ 49) and the 30 residuals (§ 52). The heliocentric true longitude at 1846.762
follows from the mean longitude, e′ and ϖ′ by Kepler's equation, the orbit in the ecliptic as his
model has it (step 3 laid his orbits in Uranus's plane; the difference is below 0.1° and is quoted).

**"Reproduced" means**, per hypothesis and entry point, with every quantity in the tier:

| tier | θ | mean longitude 1846.762 | e′ | ϖ′ | mass | residuals (§ 52) |
|---|---|---|---|---|---|---|
| 1, to his stated precision | ≤ 1′ | ≤ 2′ | ≤ 0.0005 | ≤ 10′ | ≤ 0.3 % | ≤ 0.05″ each |
| 2, within his rounding | ≤ 10′ | ≤ 15′ | ≤ 0.005 | ≤ 1° | ≤ 2 % | ≤ 0.5″ each |
| not reproduced | anything beyond tier 2 | | | | | |

The thresholds come from his print: angles to 1′, e′ to five figures, the mass to four, residuals
to 0.01″; tier 2 allows what two stages of 4-decimal elimination and an iteration stopped at 1′ can
move. The tier is reported per quantity and the worst quantity names the tier. C0 is expected to
reach tier 1 if the transcription is right; the result that matters is C2's tier, and where the
chain first leaves tier 1 if it does.

**Stage D, procedure, distance and the data–method split** (all on the heliocentric true longitude
at 1846.762 unless the item says otherwise):

- D1 *Procedure.* The same 30 (29 without 1690) equations and the same model solved as one
  weighted nonlinear least-squares problem in the eight parameters (δε, δn, δx₁, δy₁; θ, m′, e′, β′),
  at α = 0.5 and 0.515, in three weightings: (a) unit weights on all 29; (b) modern only; (c) two
  groups weighted by the § 52 residual RMS of the modern and the ancient places. (c) is the primary,
  **with the caveat that its weights come from the residuals of his own solution, so it is
  conditioned on his answer; (a) and (b) bracket it.** Reported per weighting: the four elements,
  the formal 1σ from the covariance with σ scaled to the residual RMS, the mean longitude at
  1846.762 against his printed value, and the true longitude.
- D2 *Distance, in his framework.* B1–B2 recomputed on a grid α = 0.40 … 0.60 by 0.01 (a from 48 to
  32 AU), n′ = n α^{3/2}; **the equations of condition rebuilt at each α from the fixed observation
  dates** (1810.328 + 3.0362 t in years for the modern places, the I9 epochs for the ancient), every
  rate in degrees per year, his unit of time dropped since it changes with n′; D1(c) solved at each
  α; the χ² profile, its minimum, and Δχ² at α = 0.5 and 0.515; the point solutions at α = sin 29°,
  30°, 31°, 32° beside Brookes's Table X. **Truncation and the small
  divisor:** his series keeps i = 1, 2, 3 of the F_i terms, the first-order e′ terms of four
  arguments and two second-order terms; F₂ carries the divisor (z² − n²) with z = 2(n − n′), which
  vanishes at n′ = n/2, α = 0.630. Neptune's α = 0.638 is on the far side of that pole, where no
  truncation of this series is meaningful, so the grid stops at 0.60, α ≥ 0.58 is marked
  truncation-limited in the reading, and **D2 cannot speak to Neptune's distance**; carrying the
  series further would not cure a small divisor, so no extended-series variant is run. Mean motions
  at each α from n′ = n α^{3/2}; the 1846.762 longitude from ε′ + n′ × 36.434 yr + 30′.
- D3 *Data against method, one truth.* The truth is the project's model with Uranus and the
  perturber started from the cached Horizons states of Uranus and Neptune at the integrator's epoch
  (no fitted quantity in the truth; the group offsets zero). Three noise-free data sets and three
  fits, all at a held at 38.4 AU (α = 0.5):
  - (i) *his data layout, his method*: Uranus's heliocentric true longitude at his 29 epochs from
    the truth; "C" = the same model without the perturber, Uranus's four coordinates fitted to the
    1781–1821 part only (Bouvard's base) with `fork.inference.fit.fit_uranus_elements`; O − C
    converted by r²/(ab) to mean-longitude corrections; solved by D1(c) and, if C2 reached tier 2,
    by his own procedure;
  - (ii) *his data layout, the project's model*: the same 29 noise-free heliocentric longitudes as
    `LongitudeObservations`, `fit_perturber_multistart` with `a_fixed_au=38.4` and the headline's λ₀
    starts (a best fit, no posterior);
  - (iii) *set D's layout, the project's model and sampler*: synthetic RA and Dec at set D's 278
    measurement epochs from the truth (`obs.predict` of the step-3 observation object with the
    synthetic values substituted), with H1(G)'s error model (the frozen σ_g, τ_g, offsets
    marginalised as in step 3), run through step 3's H1(G) job unchanged (two e cells, nlive 1000,
    seed namespace `[2026, 11, 3, 7, m, k+1]`, about 7 CPU-h, 45 min wall); reported as median and
    68 % as step 3 reports, with the best-fit sample beside.
  Reported for each: its offset from the truth O_x = (longitude at 1846.762) − (Neptune's from the
  same model), and the pairwise differences; (ii) against (iii) is the effect of the data layout in
  one model, (i) against (ii) the effect of the model on one layout; both enter the Q3 reading.
  Nothing stochastic enters the data (a noisy version would be a new pre-registration).
- D4 *Placement against step 3.* The step-3 fixed-a posteriors re-derived at 1846.762 from the cached
  samples (the same 6000-sample subset and `SUMMARY_RNG`; a heliocentric quantity needs no Earth
  state, so nothing is fetched); the reproduced solution (C2), D1(a–c) and his printed solution
  placed in them by the step-3 rule (inside ≤ 68 %, marginal ≤ 95 %, outside); the step-3 gap and
  shortfall restated at that epoch.

### 5 Reading, fixed in advance

Read on Hypothesis I first (S5 is the solution step 3 tested most), Hypothesis II beside it; a
reading that differs between the hypotheses is reported per hypothesis.

**Q1 (reproduction):**

| outcome | reading |
|---|---|
| C0, C1 and C2 all tier 1 or 2 | His printed solution follows from his printed data by his printed method; his arithmetic is reproducible. Proceed to Q2, Q3. |
| C0 tier 1–2, C2 not reproduced against Adams | His data and method as printed do not give his answer: the stage where B or C first fails is named and the discrepancy quantified; the transcription of that stage is rechecked once against the renders by both passes. If C2 is then within tier 2 of **Brookes's** figures (0.1509, 0.8267, 311° 35′, 323° 39′), his re-derivation is confirmed and the located equation is the Q1 result; Q2 and Q3 run on C2's (corrected) chain, with Adams's printed solution carried beside it. If C2 is near neither, the step reports "not reproduced from the data" and Q2, Q3 run on his *printed* equations (C0's inputs), with that said in every table. |
| C0 not reproduced | Either a transcription error in (1)–(4) and §§ 29/46 or a misreading of his nonlinear relations (§§ 28, 45); one recheck as above; if it stands, the step stops at Q1 and attributes nothing. |

**Q2 (procedure and distance):**

| outcome | reading |
|---|---|
| D1(c) within tier 2 of his solution, and (a), (b) within 3° of it on the true longitude | His shortcut (summing near-identical equations, § 21; ad hoc multipliers, § 25) did not move the answer: the procedure is not the cause of the step-3 gap. |
| D1(c) beyond tier 2, moving the 1846.762 true longitude by ≥ 3° toward the step-3 fixed-a median | The procedure accounts for a share of the gap, stated as (his − D1(c)) / (his − step-3 median). |
| D1(c) beyond tier 2, moving it away or sideways | The procedure matters but does not explain the gap; reported, no share. |
| (a) or (b) more than 3° from (c) | The weighting matters; (c)'s conditioning on his residuals is then material and the three are reported as a bracket, no single procedure reading. |
| D2: χ²(α = 0.5) − χ²_min ≥ 4, minimum at α ≤ 0.57 | His data in his model disfavour his distance within the range the series can reach; the minimum's α is reported, with the note that Neptune's α is beyond the pole. |
| D2: ≤ 1 | They do not; the preference step 3 found is not in his data as he reduced them, or not visible to his model. |
| between, or minimum at α ≥ 0.58 | weak or truncation-limited; stated, no statement. |

**Q3 (data against method, one truth).** Let T be Neptune's heliocentric true longitude at
1846.762 from the same model, O_x = (method x's longitude) − T, and O₃ the step-3 shortfall
(H1(G) median − Neptune, −10.5° on 1 Oct 1846, restated in D4). "Near the truth" means |O_x| ≤ 3°;
"reproduces the shortfall" means |O_x − O₃| ≤ 3°; "between" is everything else on the short side.
"More than 3° past the truth" is O_x > +3°. The reading depends first on where (iii) lands, then
on where (i) lands, and, where his method does better than the project's, on where (ii) lands to
split layout from model:

| (iii) set-D layout, project sampler | (i) his layout, his method | reading |
|---|---|---|
| near the truth | any | The step-3 shortfall is **in set D (this reconstruction of his data)**, not in the assumed distance: on clean data the project's own sampler at 38.4 AU finds the planet, so what pulled the real fit 10° short is in the measurements and their handling; **A1–A5 are the candidates**. His real-data solutions owe their nearness to his data. (i) is reported as a check that his method also finds the planet on clean data; if it does not, that is said, and it does not change this row. |
| reproduces the shortfall | near the truth | The shortfall is **a property of the assumed distance in the full model**, not of set D, and his method on his layout does not share it. Split on (ii): (ii) reproduces the shortfall → his **model** (the linearised perturbation theory at a wrong distance is less sensitive to that distance than the full model), A1–A5 not needed; (ii) near the truth → his **layout** (29 normal places at his epochs, the 1810.328 epoch), not his model; (ii) between → shares stated, (O_iii − O_ii)/O_iii to the layout and (O_ii − O_i)/O_iii to the model. |
| reproduces the shortfall | reproduces the shortfall | Both methods fall short on clean data at that distance: the shortfall is a property of the assumed distance, and his nearness on real data owes to his **data** (something in his normal places and reductions that clean data at that distance does not supply), or to compensating errors. Brookes (1970) bears on the second: his correction of Hypothesis I moves the mean longitude for 1846.762 from 325° 7′ to 323° 39′, **away from the truth**, so Adams's printed nearness contains at least one error that helped him; the residual gap after Brookes's correction (his 323° 39′ mean, true longitude at discovery 328° 4′ against Neptune's 327° 0′ in this model, both as C2 recomputes them) is stated beside the step-3 gap. A1–A5 and his reductions remain the candidates and the step cannot narrow them further. |
| reproduces the shortfall | between | Partial: his layout and model together remove (O_iii − O_i)/O_iii of the shortfall and his data the rest. Split on (ii) as in the row above: (ii) reproduces the shortfall → the removed share is his model's; (ii) near the truth → it is his layout's; (ii) between → (O_iii − O_ii)/O_iii to the layout, (O_ii − O_i)/O_iii to the model. Both candidates kept, the shares stated. |
| any of the above | more than 3° past the truth | His method/model carries a **positive bias on clean data** at 38.4 AU: O_i is stated as that bias and compared with his real-data solutions' +2 to +3° beyond Neptune (step 3, § 6); if the two agree within 3°, his overshoot on real data is the method's, not the data's. The (iii) column is still read by its own row for the shortfall. |
| between | any | Partial: the assumed distance explains O_iii/O₃ of the shortfall and set D the rest; (i) is read against (iii) as in the rows above and the shares are stated. |
| beyond the shortfall by more than 3°, or more than 3° past the truth | any | The noise-free full-model fit does not behave as step 3's did: the synthetic setup (no noise, Neptune's exact orbit, zero offsets) differs from set D in a way this plan did not foresee; reported with the numbers, no attribution. |

(ii) enters the rows above only where his method does better than the project's sampler; its
difference from (iii) is the effect of fitting 29 heliocentric normal places at his epochs instead
of 278 geocentric measurements in one model, and its difference from (i) the effect of the model on
one layout. Where (ii) is not called for by a row, both differences are still reported (R6).

**What this step cannot separate.** His data fitted by the project's model would need the "C" he
subtracted, Bouvard's 1821 tables with his § 7 corrections at his 30 epochs; that is not in the
cache and is not transcribed here. If Q3 reads "in set D", the step names this as the way to go
further and stops. Nor can D2 reach Neptune's distance: the pole of the F₂ divisor at α = 0.630 is
the 2:1 commensurability n′ = n/2 (Uranus making two revolutions to the perturber's one), and
Neptune's α = 0.638 lies beyond it, so a first-order theory of his form set up at an assumed α near
0.5 could not have represented Neptune's distance whatever the data said; this is description of
the method's reach, not a finding of this step, and the D2 reading is confined to α ≤ 0.60.

### 6 Report order (when it runs)

R1 transcription record (pass agreement, reconciliations, hashes; Smart and Brookes obtained or
not); R2 stage-B table; R3 stage-C table (hypothesis × entry point × quantity, tier); R4 D1 (a, b, c)
with 1σ; R5 D2 profile and figure; R6 D3 table (O_i, O_ii, O_iii and the pairwise differences);
R7 D4 placements and the restated gap and shortfall; R8 the readings of § 5 in order, one row each;
R9 what remains unseparated; R10 prior work read (§ 0) and what of this step it already contains.
Figure: his § 9 data with his reproduced curve and D1(c)'s (both hypotheses); the D2 profile with
the pole marked; the D3 offsets. Nothing is appended to `P1B_adams_fit.md`; the site is untouched;
`adams_solutions.csv` is not edited (any new printed value goes to the memoir files).

### 7 Cost and order of work

Transcription is the bulk (about 600 tokens, two passes; roughly the effort of step 1). Code: the
Laplace coefficients, the series, the equations, his eliminations, a least-squares and a profile
loop, the synthetic generation, all cheap; compute about 8 CPU-h, of which D3(iii)'s nested run is
7 (45 min wall on 6 workers) and the rest minutes. Order: A, then B as tests (a B failure stops the
step until resolved), then C, D1, D2, D3 (i, ii, then iii launched detached with progress log and
checkpoint), D4, summary. Ledger rows `p1b_adams_repro/<stage>` for D2 and D3.

### 8 Amendments made on review (2026-11-11, before freezing)

1. D3 gains (iii), the project's fixed-a fit on noise-free Neptune-truth data in set D's layout
   with H1(G)'s error model; each method's offset from the truth is reported, not only the
   differences; the Q3 table rewritten to read on where (iii) and (i) land.
2. D2 builds the equations of condition from the fixed observation dates at each α (his unit of
   time changes with n′); the truncation is stated, the F₂ small divisor at α = 0.630 identified,
   the grid cut at 0.60 and α ≥ 0.58 marked truncation-limited; no extended-series variant.
3. One comparison quantity in Q2, Q3 and D4: the heliocentric true longitude at 1846.762; the
   step-3 gap (12.25°) and shortfall (−10.51°) restated on it and the 3°/6° thresholds derived
   from it.
4. § 0 Prior work added (Sampson 1904, Smart 1947, Brookes 1970); no novelty claim until Brookes
   is read.
5. D1(c)'s weights noted as coming from his own solution's residuals; (a) and (b) bracket it.
6. The `src/fork/historical/` layout approved.
7. (After the review, before freezing.) Brookes 1970 obtained and read the same day; § 0 records
   what it contains and limits the novelty claim; Q1 gains Brookes's figures as a second target and
   the location of the error he reports as its contribution; D2 gains his Table X as a comparison.
   Smart 1947 sought once, not obtained.
8. (Q3 table review, 2026-11-11.) (ii) enters the attribution: where (iii) reproduces the shortfall
   and (i) is near the truth or between, the rows split on (ii) (shortfall → model; near the
   truth → layout; between → shares (O_iii − O_ii)/O_iii layout and (O_ii − O_i)/O_iii model).
   A row added for (i) more than 3° past the truth (positive bias of his method on clean data,
   compared with his real solutions' +2 to +3°). Row 1 reworded: the shortfall is in set D, this
   reconstruction of his data, A1–A5 the candidates. The compensating-errors row cites Brookes
   (1970): his correction moves Hypothesis I from 325° 7′ to 323° 39′, away from the truth, and
   the residual gap after it is stated beside the step-3 gap. The D2 pole named as the 2:1
   commensurability under "what this step cannot separate". D3(iii) confirmed as a rerun of the
   step-3 H1(G) nested job on the synthetic data, detached, best-fit sample beside the median.
