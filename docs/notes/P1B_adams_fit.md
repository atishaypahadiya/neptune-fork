# Phase 1b, step 3: the Adams-variant fit, pre-registration

Plan only. Nothing in this step runs until the user approves this note. It changes no published
result: the headline, S1–S3, M9 and M10 stand as they are. Written 2026-11-10, before any step-3
code.

## Pre-registration (before any code or run)

**The question.** Given Adams's own data choices and assumptions, what does the posterior say
about the planet's place, mass and orbit, and how do the solutions he communicated compare with
it? Three fits answer it: a free-a fit on his data, and two fits with a held at his two assumed
mean distances. The comparison targets are fixed below, before any run.

**Launch condition.** Nothing runs until all of these are true:
- this plan is committed and pushed;
- the user has approved it (amendments, if any, appended below the line, nothing above it edited);
- the pre-run changes below are committed, each with its hash recorded here;
- the step-3 code is committed and its hash recorded here.

**Frozen inputs (SHA-256):**
- `data/historical/leverrier_1846.csv` 2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4
- `data/historical/m6_error_model.json` 0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277
  (σ_g, τ_g and the `sigma_t4` column)
- `data/historical/adams_challis/adams_solutions.csv` 9952542fa6283d7e0f67db58714173a0a17d746f41aa3cae7c2cd199088347ea
  (the targets S4–S7)
- `data/historical/adams_challis/challis_sweeps.csv` f9061de02545168636832645cca164260063ccf50ec1e1498b06950fd450fafa
  (the dates Challis observed)

### 1. Adams's data, from the memoir

Source: J. C. Adams, *An Explanation of the Observed Irregularities in the Motion of Uranus*,
Nautical Almanac 1851 appendix, separate printing (1846), cached as
`adams_1846_nautical_almanac_appendix_archive_org.pdf` (page images `pamphlet_pNNN.png`, printed
page = PDF page − 6). Marks as in the Phase 1b notes: P primary-confirmed.

| # | Adams's choice | source | what we do |
|---|---|---|---|
| D1 | **Ancient observations: nine epochs**, 1690, 1712, 1715, 1750, 1753, 1756, 1764, 1769, 1771, one residual each (§ 8, p. 6; § 9, p. 7). "The errors of the Tables for the ancient Observations were taken from those given in the Equations of Condition of Bouvard's Tables" (§ 4, p. 4). | P, pamphlet pp. 4, 6–7 | Table 77's 19 pre-1781 observations cover every one of these epochs: Flamsteed 1690 (1), 1712 (1), 1715 (3); Le Monnier 1750 (2); Bradley 1753 (1); Mayer 1756 (1); Le Monnier 1764 (1), 1768–69 (8, Adams's "1769"), 1771 (1). The rows stand in for the epochs. |
| D2 | **1690 set aside in the final equations.** "It was thought unsafe to employ the first of the above equations, since it is derived from the single observation of Flamsteed, made in 1690, twenty-two years anterior to any other observation" (§ 25, p. 15). It stays in the residual tables (§ 8–9, § 52). | P, pamphlet pp. 15, 27 | Flamsteed 1690 dropped (1 observation, RA and Dec, 2 measurements). This is the spec's `adams` variant (`docs/PHASE1_SPEC.md`, M6c table). |
| D3 | **Ancient weights.** The 1769 and 1771 residuals "though not entitled to the same weight as the others, cannot, I think, be liable to much uncertainty" (§ 6, p. 5); the ancient equations are combined "by multiplying each equation by a quantity nearly proportional to the co-efficient of each of the unknown quantities p₃ and q₃" (§ 25, p. 15): a solution-driven weighting, not an error weighting. | P, pamphlet pp. 5, 15 | Not reproduced. The measurements keep the frozen error model's σ_g and τ_g (ancient group 4.66″/5.08″, τ 5″). Stated as approximation A1 below. |
| D4 | **Modern observations: Greenwich only to 1830,** "with the exception of an observation by Bessel, in 1823; and subsequently from the Cambridge and Greenwich Observations, and those given in various numbers of the *Astronomische Nachrichten*" (§ 4, p. 4). | P, pamphlet p. 4 | Table 77's modern rows are Greenwich 1781–1830 (87 observations, the same Greenwich series Adams used, in Le Verrier's reduction) and Paris 1801–1827 and 1835–1845 (Bouvard's and Arago's reductions). **Paris 1801–1827 is dropped** (107 observations: Adams did not use Paris). **Greenwich 1781–1830 is kept.** Bessel 1823 has no Table 77 row and is omitted. |
| D5 | **The series ends at 1840.** Twenty-one normal places 1780–1840 at three-year intervals (§ 6, p. 5; § 8–9). "In my last calculations, I might have used more recent observations, but ... the later observations might be used as a test of the theory" (§ 6, p. 5). The 1844–45 Greenwich normal errors reached him on 7 September 1846 (Airy, MNRAS 7, 140) and are used as a test in § 53. | P, pamphlet p. 5; Airy p. 140 | Rows after 1840-12-31 are **dropped from the fit** (Paris 1841–1845, 26 observations, 47 measurements) and **used as the test**: the posterior predictive residuals at those epochs are reported, as Adams's § 53 does with his 1843–45 errors. |
| D6 | **1831–1840: Cambridge, Greenwich and the AN.** Normal places 1831, 1834, 1837, 1840 (§ 8–9). | P, pamphlet pp. 4, 6–7 | Table 77 has no row in 1831–1834 and Paris rows only in 1835–1840 (40 observations). **Paris 1835–1840 stands in** for the 1834–1840 normal places; the 1831 normal place has no counterpart. Approximation A2. |
| D7 | **Three-year normal places near opposition** "so that the result should be nearly free from the effects of" the radius-vector error (§ 6, p. 5); 1780 interpolated from 1781–83 (§ 6). | P, pamphlet p. 5 | Not reproduced: every measurement enters individually, with its own date and full geometry, so the radius-vector effect is modelled, not avoided. Approximation A3. |
| D8 | **Reductions:** Bouvard's Tables with Adams's corrections for Jupiter's and Saturn's terms and Hansen's second-order terms (§ 7, pp. 5–6); residuals converted from true to mean longitude (§ 9, p. 7); heliocentric longitudes only, no latitudes, no radius vector. | P, pamphlet pp. 5–7 | Not reproduced. The fit is the project's modern track: apparent RA and Dec against the integrated Uranus (`docs/notes/M6.md`, M6b). Table 77's measured values are what Adams's observers measured; his and Le Verrier's reductions are both replaced by the forward model. Approximation A4. |
| D9 | **Uranus's elements free:** δε, δn, δe, eδϖ (§ 10, p. 7); a fixed by n. | P, pamphlet p. 7 | Uranus's four in-plane elements free (docs/RULES.md rule 2), i and Ω held, as every Phase 1 fit. |
| D10 | **Perturber in Uranus's plane, circular Uranus terms to first order in e′** (§ 4, § 11–12, pp. 4, 8). | P | The perturber's inclination keeps the headline's prior (i ≤ 10° relative to Uranus's plane), not held at zero: holding it would be a fourth change and Adams's choice was a computational one. Stated as a difference, not an approximation of his data. |

**The Adams data set (D):** Table 77 minus Flamsteed 1690, minus Paris 1801–1827, minus every
row after 1840. Counts, computed from the frozen CSV:

| group | observations | measurements (RA + Dec) |
|---|---|---|
| ancient (1712–1771) | 18 | 35 (18 + 17) |
| Greenwich 1781–1830 | 87 | 167 (87 + 80) |
| Paris 1835–1840 | 40 | 76 (40 + 36) |
| **total** | **145** | **278** |
| test rows (Paris 1841–1845, not fitted) | 26 | 47 |

The dropped Paris 1801–1827 group has no rows left: its offset is 0 and its prior term absent
(`m9.drop_measurements`). The other groups keep their σ_g, s_g and τ_g.

**Approximations, stated (not reproduced from the memoir):**
- **A1** his weighting of the ancient equations (D3): ours is the error model's.
- **A2** the 1831 normal place has no Table 77 counterpart, and 1834–1840 are Paris rows, not the
  Cambridge, Greenwich and AN observations he used (D6).
- **A3** normal places and the opposition selection (D7): we fit individual measurements.
- **A4** his reductions (D8): the modern track replaces both his and Le Verrier's.
- **A5** Bessel 1823 omitted (D4).
- The Greenwich 1781–1830 rows are the same observations he used but in Le Verrier's reduction,
  a different reduction of the same transits; this is what "Table 77 rows stand in" means for D4.

What the data set reproduces exactly: the epochs and observers of the ancient series, the
exclusion of 1690, the exclusion of the Paris series, the end of the fitted series at 1840 and
the use of the later rows as a test.

### 2. Assumptions: the three fits

| fit | perturber a | prior on a | cells |
|---|---|---|---|
| H1 | **held at 38.4 AU** (a/a′ = 0.5, a_Uranus = 19.2; memoir § 12, Airy No. 11, Sampson table V) | — | two, split in e at 0.12 (the headline's split for a ≥ 28 AU) |
| H2 | **held at 37.25 AU** (a/a′ = 0.515 = sin 31°, memoir § 32; 37.25 in Sampson table VI; 1.942 × 19.18 = 37.25) | — | two, as H1 |
| F | **free** | log-uniform on [20, 60] AU, the headline's | the headline's four ([20, 25]; [25, 28]; [28, 60] × e < 0.12; [28, 60] × e 0.12–0.3) |

Everything else is the headline's (`docs/notes/M6.md`, "M6c pre-registration"): log₁₀ m/M⊕
uniform on [0, 2]; e uniform on [0, 0.3]; ϖ and λ₀ uniform; i ≤ 10° relative to Uranus's plane,
cos i uniform, Ω uniform; Uranus's (a, e cos ϖ, e sin ϖ, λ₀) in the ±50 scaled-σ box centred on
the Uranus-only RA/Dec fit **of data set D** (the M5 rule applied to D, as every variant applies it
to its own data), i and Ω held at their 1781 values; eight group offsets b_g ~ N(0, τ_g²),
τ = 2″ modern, 5″ ancient, marginalised (analytically under Gaussian errors, by quadrature under
Student-t); observation model apparent RA/Dec of date with light time, aberration and
IAU 2006/2000A.

Held a is set through `PosteriorSpace.held` (ln a held, the other perturber entries free), the
mechanism `held_shape_space` already uses for Le Verrier's shape. Held a means the osculating a
at the 1781 epoch about the inner barycentre, the model's a, equated to Adams's mean distance;
the difference between that and his mean-distance definition is far below the 1.15 AU between
H1 and H2 and is not corrected.

**Error models: six runs.** Each fit under
- **G**, Gaussian, the headline's σ_g (the headline result of this step), and
- **T**, Student-t with ν = 4 and S3's scales `sigma_t4` (ancient 4.796″/3.584″, Greenwich
  1.784″/2.095″, Paris 1835–1845 1.366″/1.076″), the headline's τ_g, offsets integrated by
  48-node Gauss–Hermite quadrature (`StudentTObservations`, as S3 and M10).

The Gaussian run is the headline of the step; the Student-t run is reported beside it in every
table, as M10 is beside M9. The Uranus box is the Gaussian rule's for both (as S3 used the
headline's), so the prior is identical within each fit.

**Sampler:** dynesty as the headline: bound `multi`, sample `rslice`, dlogz 0.1, nlive 2000 in
cell 0 of F and 1000 elsewhere, ln Z errors ×1.8. Seeds, a new namespace: sampler rng
`[2026, 11, 3, f, m, k + 1]` with f = 0, 1, 2 for H1, H2, F, m = 0, 1 for G, T, k the cell index;
summary rng `[7, 1846, 11, 3]`. Derived quantities from 6000 equal-weight posterior samples per
run (the headline's count). Each job checks the frozen hashes, τ, dlogz and the i prior and
refuses to start if any differs. Checkpoint files per cell; a job resumes from its checkpoint.

### 3. Comparison targets, fixed now

From `adams_solutions.csv` (step 1, every value primary) and the ephemeris of this repository.
Each target is propagated as a Kepler orbit from its printed elements, in Uranus's plane
(Adams assumed the ecliptic; the difference is below 1′ in longitude), to each report date;
S4 has no mean distance of its own beyond 38.4 AU and its printed elements suffice.

| label | target | elements used | class |
|---|---|---|---|
| S4 | 18 Sep 1845, the paper to Challis | a 38.4, e 0.1428, ϖ 320° 30′, mean longitude 321° 40′ at end Sept 1845 (taken as 1845-09-30), m 0.000173 M☉ | communicated before the search |
| S5 | Hypothesis I, 21 Oct 1845 to Airy | a 38.4, e 0.16103, ϖ 315° 57′, mean longitude 325° 7′ on 1846-10-06, m 0.0001656 M☉ | communicated before the search |
| S6 | Hypothesis II, 2 Sep 1846 | a 37.25, e 0.120615, ϖ 299° 11′, mean longitude 323° 2′ on 1846-10-06, m 0.00015003 M☉ | final computed |
| S7 | 2 Sep 1846 extrapolation | mean longitude about 315° 20′ on 1846-10-01, a ≈ 33.7 AU (1/0.57), e → 0, no mass | **reported separately as an extrapolation, not a solution**; it enters no table of solutions and no summary sentence about Adams's prediction |
| LV-June | Le Verrier, Comptes rendus 1 June 1846 | heliocentric true longitude "about 325°" for the beginning of 1847, mean distance twice Uranus's, orbit in the ecliptic, no other elements (Airy's summary, MNRAS 7, 131–132, P for the summary; the Comptes rendus text is not held) | comparison |
| LV-Aug | Le Verrier, 31 Aug 1846 | `LE_VERRIER_1846`: a 36.15, e 0.10761, ϖ 284° 45′, mean longitude 318° 47′ at 1847.0, m 0.00011 M☉; true heliocentric longitude 326° 32′ at 1847.0 (Airy, MNRAS 7, 141) | comparison |
| Neptune | JPL Horizons, cached | the cached states; 325.873° geocentric on Galle's night (M6c constant), 327° 0′ heliocentric mean of date on 1 Oct 1846 (step 1) | truth |

LV-June is carried to the report dates as a circular orbit of 38.4 AU through 325° at
1847-01-01 (0.38° of motion over the quarter); its uncertainty, "an error of 10° in this place is
not probable", is quoted beside it. Kollerstrom's 23 September conversions (324° 35′, 325° 58′)
are secondary and are listed for the record only.

### 4. Report dates

- **1846-10-01**, Adams's date in the 2 September letter (and the search's end).
- **1846-09-23 23:00 TDB**, Galle's night, the headline's quantity, so that every table can be
  read against the headline (325.82°, 68 % 323.57–327.79°).
- **The nights Challis recorded or missed the planet**, from `challis_sweeps.csv`: 1846-07-29,
  1846-07-30, 1846-08-04, 1846-08-12, 1846-09-29. Each at 23:30 UT.
- **1845-09-30**, for S4 only: its own date, so that the paper to Challis is compared on the
  night it spoke of (geocentric 320° 30′, mean 321° 40′).

### 5. Reported, per run

For each of the six runs and each report date, from the 6000 equal-weight samples:
1. **Geocentric ecliptic longitude of date, light-time corrected** (the headline's Galle-night
   quantity, `GalleGeometry.lon_ecl_of_date_rad` with the date's Earth position): median, 68 %
   and 95 % circular intervals; for each target, median − target in degrees and in 68 %
   half-widths, and whether the 68 % and 95 % intervals contain it. The target's geocentric
   longitude is derived from its elements by the same code; S4's printed geocentric 320° 30′ is
   compared on 1845-09-30 as printed.
2. **Heliocentric true longitude of date** and **heliocentric mean longitude of date** (per
   sample: λ₀ + n t, in the model's frame, precessed to the date), the same statistics. Adams's
   printed numbers are mean longitudes, so the mean-longitude row is the like-for-like one; the
   true-longitude row is where the planet is.
3. **Geocentric RA and Dec of date** on the five Challis nights, median and 68 %, against the
   zone verdicts of `challis_sweeps.csv` (9′ zones: the posterior's declination spread against
   the 9′ breadth is reported as description; no new verdict on the sweeps is drawn here).
4. **Mass, e, ϖ, i:** median, 68 %, 95 %; for H1 and H2 also the mass against the target's
   printed mass (S5 0.0001656, S6 0.00015003 M☉) in posterior half-widths; **a** for F.
5. **ln Z** per cell and in total (errors ×1.8), cell shares, χ²_min (G) or Q_min (T), the
   best-fit group offsets.
6. **The test rows:** posterior-predictive residuals of the Paris 1841–1845 measurements (median
   and 68 % of O − C per epoch), beside Adams's own 1843–45 test errors (+7″.11, +8″.79,
   +12″.40 under H1; +5″.77, +7″.05, +10″.18 under H2; memoir § 53).
7. **Against the headline:** the Galle-night median shift in degrees and in headline 68 %
   half-widths, and the 68 % width ratio, for F(G) and F(T). No ln B against the headline (the
   data differ).
8. **H1 against H2 against F**, within each error model: ln B(H1 : F) and ln B(H2 : F) are not
   formed (held a against free a compares a point with a prior volume); instead Δ ln L_max and
   the posterior density of F at 38.4 and 37.25 AU (bins of ±0.5 AU, effective sample counts
   given) are reported.

### 6. Reading rule (fixed now)

Thresholds as M9/M10: a target is **inside** a posterior if it lies within the 68 % interval,
**marginal** between 68 % and 95 %, **outside** beyond 95 %; distances in 68 % half-widths are
given either way. The reading is made on the Gaussian run and reported for the Student-t run
beside it; where the two disagree on inside/marginal/outside for any target on 1 October 1846,
the report says so and gives both, and the Gaussian reading stands as the pre-registered one.

Three pre-registered statements, each read on the geocentric longitude of 1846-10-01 under F(G)
unless it names another run:
1. **Did Adams's data locate the planet?** Neptune inside F(G)'s 68 % interval: yes;
   marginal: partly; outside: no. The same for H1(G) and H2(G), reported as "under his assumed
   distance".
2. **Were his communicated solutions representative of his data?** For each of S4, S5 (under
   H1) and S6 (under H2): the target's true longitude inside, marginal or outside the fixed-a
   posterior's 68/95 % intervals. A solution outside its own fixed-a posterior means his point
   estimate was not where the data under his assumption put the planet; inside means it was.
   This is read on the true longitude (where the planet is) and the mean longitude (what he
   printed) separately; the report says if they differ.
3. **Did the data prefer his assumed distance?** F's posterior density in a at 38.4 and 37.25 AU
   relative to its maximum, and Δ ln L_max(H1) and Δ ln L_max(H2) against F's maximum. Read as
   M9's thresholds: Δ ln L ≤ −2 means the data disfavoured his distance; above −0.5 means they
   did not distinguish it.

S7 is excluded from all three. Le Verrier's June and August positions are reported against
every posterior in the same inside/marginal/outside terms, as comparison, and enter no statement.

Interpretation table (pre-registered; rows apply only where the reading rule's clauses hold):

| 1 | 2 (S5, S6) | 3 | reading |
|---|---|---|---|
| yes | inside | disfavoured | His data located the planet; his solutions sat where his data put it under his assumed distance; the distance was the error, and the data could have told him. |
| yes | inside | not distinguished | His data located the planet and his solutions represented them; the data of 1840 did not fix the distance, so his assumption was a free choice within what they allowed. |
| yes | outside | any | His data located the planet but his solutions did not represent its posterior under his own assumption: the gap is in his solution method or data handling (A1–A5 are the candidates), not in the data. |
| partly or no | any | any | His data, as reconstructed here, did not locate the planet; the comparison of his solutions is reported as description and no statement about his method is made. |

### 7. Pre-run changes (each committed before the launch, hash recorded here)

- **(a) Earth states at the report dates.** The cached Earth states hold no 1846 epoch but
  Galle's night. The seven report epochs (1845-09-30 23:30 UT; 1846-07-29, 07-30, 08-04, 08-12,
  09-29 at 23:30 UT; 1846-10-01 23:30 UT) are added to `fork.data.horizons`'s `--fetch` plan for
  Earth and Neptune and fetched **once, online, before the launch**; the two cache files' SHA-256
  after the fetch are recorded here. The run itself is `FORK_HORIZONS_CACHE_ONLY=1` and makes no
  network call (M10's fetcher). The CI fixture is regenerated with the same epochs by its build
  script so that the summary test runs offline.
- **(b) The data-set builder:** a `p1b_adams` variant in `fork.export.m6c.dataset`'s style,
  filtering Table 77's rows by observer, observatory and date as in § 1, with the test-row split;
  `tests/test_p1b_adams.py` checks the counts in § 1's table exactly, the held-a space (ln a held,
  nine free perturber-and-Uranus entries plus the two inclination entries), the seeding, and that
  the builder reads no Neptune state.
- **(c) The runner:** `src/fork/export/p1b_adams.py`, M6c's job structure (variant × cell, pkl +
  ckpt per job, `check_frozen`), six variants `h1_g, h1_t, h2_g, h2_t, f_g, f_t`, 2 + 2 + 4 cells
  each error model = 16 jobs; `summarize` writes `docs/notes/P1B_adams_fit.json` with every
  quantity of § 5 and `docs/notes/P1B_adams_fit.png` (longitude posteriors on 1 Oct 1846 with the
  targets marked; F's a marginal with 38.4 and 37.25 marked). No site payload unless the user
  puts the result on the site.

### 8. Launch and logistics

- Detached, as M10: its own process session, `nohup`, under `caffeinate -is`, parent not the
  session that starts it; `FORK_HORIZONS_CACHE_ONLY=1`; progress to
  `data/cache/p1b_adams/progress.log`, output to `data/cache/p1b_adams/run.out`; one ledger row
  per job, `p1b_adams/<variant>/<cell>`; resumable per job from its checkpoint.
- **Order:** H1(G), H2(G), F(G), then H1(T), H2(T), F(T); within a fit, cells in parallel on 12
  workers (F: 6 on cell 0, 2 on each other cell, as the headline).
- **Compute, projected** from the ledger: the headline (4 cells, 529 measurements) cost
  28.5 CPU-h and 3.1 h wall; S3 (Student-t) 30.4 CPU-h. D has 278 measurements, so a likelihood
  call costs about half. F(G) and F(T): about 15 CPU-h each. H1 and H2: one dimension fewer and
  two cells, about 6 CPU-h each. Total about 55 CPU-h, 7–9 h wall on 12 workers. Reported against
  the projection afterwards, as every milestone.
- **Smoke run** (`--maxiter`, not saved, not a result) of H1(G) cell 0 before the launch, to
  check the held-a space and the summary path end to end.

### 9. What this step does not do

- It does not reproduce Adams's reductions, weights or normal places (A1–A5). A later step could
  fit his printed normal-place residuals (§ 9) directly with his equations of condition; that is a
  different model (a linearised one about Bouvard's tables) and is not pre-registered here.
- It does not fit the 1845 state of his data (Hypothesis I was computed before he had the
  1831–1840 Cambridge and AN material in its final form); D is the memoir's data set, the one
  Hypothesis II and the restated Hypothesis I rest on.
- It does not model the sweeps; the Challis-night rows are positions for `challis_sweeps.csv`'s
  verdicts to be read against, not a coverage estimate.
- It draws no conclusion about priority or credit.

### 10. Outputs

- `docs/notes/P1B_adams_fit.json`, `docs/notes/P1B_adams_fit.png`, results appended to this note
  below the amendment line in the order of § 5 and § 6; ledger rows; `docs/reproduce.md` commands
  and the new ledger total.

## Amendment (2026-11-10, on approval; before any step-3 code)

The plan above is approved with four amendments. Nothing above this line is edited.

1. **Inclination:** as planned. The perturber keeps the headline's prior (i ≤ 10° relative to
   Uranus's plane); Adams's in-plane assumption stays listed as a difference (D10).
2. **LV-June:** before the launch, Le Verrier's paper in the *Comptes rendus* of 1 June 1846
   (vol. 22) is sought from Gallica, once, online; its URL and SHA-256 are recorded in
   `data/historical/adams_challis/MANIFEST.md`. If obtained, LV-June uses its printed value, marked
   P; if not, Airy's summary as planned, marked as the fallback. LV-June remains comparison only.
3. **Statement 3** gains a middle band: −2 < Δ ln L ≤ −0.5 reads "weakly disfavoured; reported, no
   statement". So: Δ ln L ≤ −2 disfavoured; −2 < Δ ln L ≤ −0.5 weakly disfavoured, no statement;
   above −0.5 not distinguished.
4. **Interpretation table**, two additions. (i) Statement 2 = marginal (a solution between the 68 %
   and 95 % intervals of its fixed-a posterior): reported as "at the edge of what his data allowed
   under his assumption"; no row of the table is applied and no statement about his method is made.
   (ii) S5 and S6 disagreeing (one inside, the other outside or marginal): the reading is reported
   per solution, and no single statement about "his solutions" is made.

The pre-run changes (a)–(c), the hash records, the smoke run and the detached launch follow this
commit, in § 7–8's order.

---

## Pre-run record (2026-11-10, after the amendment; nothing above the lines edited)

**Amendment 2, LV-June.** Gallica refused every request (HTTP 403 on the SRU search API, the
periodical's date page `ark:/12148/cb343481087/date1846` and four IIIF manifests), so the volume
was fetched once, online, from archive.org instead: `comptesrendusheb22acad` (Comptes rendus t.
22, janvier–juin 1846, Natural History Museum Library copy),
https://archive.org/download/comptesrendusheb22acad/comptesrendusheb22acad.pdf, SHA-256
`555936e73efaca9448b38ecd9c7f4aa62f405d232c44e582a637f92aa693dac9`, cached (not committed) as
`data/cache/adams_challis/cr22/comptesrendusheb22acad.pdf`; entry and page mapping in
`data/historical/adams_challis/MANIFEST.md`. Le Verrier's paper is at printed pages 907–918
(PDF 921–932). Page 917 (PDF 931) prints the value: "en assignant 325 degrés de longitude
héliocentrique à la planète, au 1er janvier 1847, on ne commet pas une erreur de 10 degrés", with
the formula ν = 314°,5 + 12°,25 α + (1/m){20°,82 − 10°,79 α − 1°,14 α²} (m the mass in
ten-thousandths of the Sun's, α undetermined). **LV-June therefore uses the printed value, mark
P**; Airy's summary agrees ("about 325°") and is no longer the fallback. The target is carried as
planned (§ 3): a circular orbit of 38.4 AU through 325° at 1847-01-01.

**Pre-run change (a), Earth and Neptune at the report epochs.** `fork.data.horizons` gained
`P1B_REPORT_EPOCHS_ISO_TDB` (the seven epochs of § 4, 23:30 UT written as TDB; the minute of
ΔT + TDB − UT is irrelevant at the arcminute level of the comparison) and `default_fetch_plan`
adds them, with Galle's epoch, for Earth, Neptune and the four inner bodies (the inner-barycenter
correction needs them). Fetched once, online, 2026-11-10; cache SHA-256 before → after:

| file | before | after |
|---|---|---|
| `data/cache/horizons_earth.npz` | `4f422ca0e3a0729418d7e625df6496cc662564fdd6b4af2ca3ef5db8b8611829` | `a73ea1c906131a01a184b16718ce5e4bc1a180b370a1f9d9411ce3eabd958332` |
| `data/cache/horizons_neptune.npz` | `d316ae71e79170b1362db1b7bd22afa63b1bcd5fa7d39f9c8c82b55691c4f0ad` | `559f361fc10bdde74f898a6e3c328b715b0790b48be447768592ef149146c6ee` |
| `data/cache/horizons_mercury.npz` | — | `3f102e3c29c94cbb02c7dec60f77a55c3a5dd644231158cbacf7bfc69eb68b77` |
| `data/cache/horizons_venus.npz` | — | `a62d5c340fc78b7b788f43df9a2d927f55bf7084d8a4587a2fa03c9c63de43e5` |
| `data/cache/horizons_earth_moon.npz` | — | `120a54df2c70069e34141a85905fd89b4715fa1dfbd3aa33339e8215f09cc1c9` |
| `data/cache/horizons_mars.npz` | — | `2f1ba5da00bc2d59f8a8ab9f5a8651b6e6f2f3a21921d2f02f45e343b0ff0bd9` |

Earth and Neptune now hold 166 epochs each (158 + Galle + the seven), the inner bodies 9 (A1 +
Galle + the seven). The CI fixture `tests/fixtures/horizons/` was regenerated by its build script
(hashes below). The run is `FORK_HORIZONS_CACHE_ONLY=1`.

**Pre-run changes (b) and (c), the code.** Committed before the launch; SHA-256 at the commit:

| file | role | SHA-256 |
|---|---|---|
| `src/fork/data/horizons.py` | (a) the report epochs in the fetch plan | `cfebb972725797f1b0ca54083a17a2f5d931f2a8d4b5b922618edc7fd85e875e` |
| `src/fork/export/p1b_adams.py` | (b) data set D and the test rows (`adams_rows`, `dataset`), (c) the runner: held-a spaces, 16 jobs, scheduler, targets, summary, reading rule | `c9c4d494ac441a4dfc77a8795e9c28da9fa603b74178bf5b381d8fb3c5949fcd` |
| `src/fork/export/p1b_adams_plot.py` | (c) the figure | `7362b44f765d0e41b7ecd25d0e2d91c213f37711544acbdfd49179c898ef6e3a` |
| `tests/test_p1b_adams.py` | (b) the checks of § 7 (b): counts 145/278 and 18/35, 87/167, 0/0, 40/76, test rows 26/47; the held-a space (ndim 10: Uranus 4 + log m, (k, h), λ₀ + the two inclination entries, ln a held at ln 38.4); seeds `[2026, 11, 3, f, m, k+1]`, disjoint from M6c's; no Neptune file opened by the data-set builder or the spaces; the targets against the printed places and Neptune; the reading rule's bands and rows | `2ec1512c47ee76bfc01c8ee5056c7c65311987f3aee5dc14472e0376d4619e77` |
| `tests/fixtures/horizons/horizons_earth.npz` | CI fixture, regenerated (296 epochs) | `707ac36efcc8af119bf5812b5287b500deafbd35bcec84be920c966c18c56ae8` |
| `tests/fixtures/horizons/horizons_neptune.npz` | CI fixture, regenerated (9 epochs) | `b13f703744bb6f3325f7fdfa7faffc6f5603c067c057a584ba5d3a1347f4704a` |

Implementation notes, none of them a change to the plan:
- Data set D is built as M9 built its data sets: Table 77 in full through the frozen error model,
  then the unused measurements dropped (`m9.drop_measurements`), because the observation layout
  requires every group populated at construction; the empty Paris 1801–1827 group keeps σ and τ
  and contributes nothing. The test rows are the same layout restricted to Paris 1841–1845.
- The held-a space is `PosteriorSpace(held=(None, ln a, None, None, None), inclined=True)` with
  an e slice per cell; its two cells' prior masses are 0.4 and 0.6 and sum to 1, as
  `PartitionedResult` requires.
- Targets are propagated from their printed elements in Uranus's plane (§ 3). Checks made by the
  tests: Neptune's geocentric longitude on Galle's night agrees with M6c's constant to 0.01°; its
  heliocentric true longitude on 1 October 1846 is 327.0° (step 1's figure; the osculating mean
  longitude is 328.2°, 1.2° ahead, e ≈ 0.01); S5's mean longitude carried from the memoir's
  325° 7′ (6 Oct) to 1 October is 325° 6′, 2′ from Airy's printed 325° 8′, the two printings
  themselves 2′ apart; S4's geocentric longitude from its elements on 30 September 1845 is
  320° 48′ against the sheet's "320° 30′ nearly, diminishing about 1′ daily", so both are reported
  (§ 5.1, `S4` and `S4_printed`); Neptune's mean-of-date RA/Dec on 12 August 1846 agrees with
  step 2's table (21h 57m 24s, −13° 02′). RA/Dec are referred to the mean equator and equinox of
  date (nutation and aberration, under 20″, are not applied; the 9′ zone is the comparison).
- The mean longitude per sample is Ω + ω + M of the osculating ellipse in the ecliptic of date
  (the perturber's state at the epoch, from the same integration as the other quantities).

**Smoke runs (not results, nothing saved under `data/cache/p1b_adams`).**
- `job h1_g/c0 --n-workers 4 --maxiter 400` with `FORK_HORIZONS_CACHE_ONLY=1`: ndim 10, the
  names `u_a_au, u_e_cos_varpi, u_e_sin_varpi, u_lam0, p_log10_m_earth, p_sqrt_e_cos_varpi,
  p_sqrt_e_sin_varpi, p_lam0, p_sin_i_cos_node, p_sin_i_sin_node` (a absent, as held); 401
  iterations, 1487 likelihood calls, 2.8 s wall on 4 workers.
- The summary path end to end on a scratch cache (H1(G), H2(G), F(G), nlive 200, maxiter 600,
  one worker, 15–21 s per cell; ledger disabled): `summarize` wrote the JSON with every § 5
  quantity for the three variants, the reading with statements 1–3 and an interpretation row,
  and the PNG, in 4.8 s. Those posteriors are unconverged by construction and were discarded.

**Checks before the launch.** `ruff check` and `ruff format --check` clean; the default test
job 353 passed from the rebuilt fixture (CI makes no live call). Launch time, commit and the
ledger rows are recorded with the results below.

---

## Results (2026-11-11; the run of commit d8d2990, launched 2026-11-10 21:01:37 PST, finished 04:54:43; nothing above the lines edited)

All 16 jobs completed, no failures, no retries (`data/cache/p1b_adams/driver.log`). Ledger rows
`p1b_adams/<variant>/c<k>` (16, kind `nested`) and `p1b_adams/summarize` (estimated; the summary
was run four times, the second to fourth only to redraw the figure after fixes to
`p1b_adams_plot.py`, so four rows). Summary: `docs/notes/P1B_adams_fit.json` (every quantity of
§ 5, the targets' places, the reading) and `docs/notes/P1B_adams_fit.png`. The site is untouched;
`site/data/ledger.json` is not regenerated here, so the rows are only in the local ledger and in
this record until the site step.

| | projected (§ 8) | actual |
|---|---|---|
| CPU-hours, 16 jobs | ≈ 55 | 70.6 (H1(G) 6.8, H2(G) 6.8, F(G) 20.7, H1(T) 7.2, H2(T) 7.3, F(T) 21.8) |
| wall, 12 workers | 7–9 h | 7 h 53 min (held-a cells 40–46 min; F cells 103–146 min) |
| likelihood calls | — | 51.9 M |

The free-a fits cost about 21 CPU-h each against the projected 15 (the cells ran longer than the
headline's, at ~1.5 % sampling efficiency), the held-a fits 7 against 6. Equal-weight samples: 79–81 k
per held-a run, 208–220 k per free-a run; 6000 of each were integrated for the derived quantities.

### § 5.1 Geocentric ecliptic longitude of date (light time)

**1 October 1846.** Neptune (cache) 325.71°. Targets from their elements: S4 322.89°, S5 327.59°,
S6 328.36°, LV-June 323.56°, LV-Aug 324.83°. Median, 68 % and 95 % intervals; for each target its
status (inside ≤ 68 %, marginal ≤ 95 %, outside) with median − target in 68 % half-widths:

| run | median [68 %] [95 %] | Neptune | S4 | S5 | S6 | LV-June | LV-Aug |
|---|---|---|---|---|---|---|---|
| H1(G) | 315.20 [310.80, 318.87] [306.42, 322.43] | outside, −2.9 | outside, −2.1 | outside, −3.4 | outside, −3.6 | outside, −2.3 | outside, −2.6 |
| H2(G) | 316.57 [313.00, 320.87] [310.31, 324.70] | outside, −2.1 | marginal, −1.5 | outside, −2.6 | outside, −2.7 | marginal, −1.6 | outside, −1.9 |
| F(G) | 324.55 [316.83, 342.15] [311.89, 358.34] | inside, −0.1 | inside, +0.2 | inside, −0.2 | inside, −0.2 | inside, +0.1 | inside, −0.0 |
| H1(T) | 317.93 [313.41, 322.39] [309.45, 327.09] | marginal, −1.7 | marginal, −1.1 | outside, −2.2 | outside, −2.3 | marginal, −1.3 | marginal, −1.5 |
| H2(T) | 319.87 [314.94, 324.45] [309.18, 328.37] | marginal, −1.3 | inside, −0.7 | marginal, −1.7 | marginal, −1.9 | inside, −0.8 | marginal, −1.1 |
| F(T) | 327.72 [321.44, 335.76] [316.40, 347.01] | inside, +0.3 | inside, +0.8 | inside, +0.0 | inside, −0.1 | inside, +0.7 | inside, +0.5 |

The fixed-a posteriors put the planet 9–10° short of Neptune; Adams's solutions put it 2–3° beyond
Neptune; the free-a posterior is centred on Neptune with a 68 % width of 25° (G) and 14° (T).

**Galle's night (1846-09-23 23:00 TDB).** Neptune 325.87° (M6c's constant). Medians 0.1° above the
1 October ones: H1(G) 315.31 [310.89, 318.99], H2(G) 316.68 [313.10, 321.00], F(G) 324.70 [316.94,
342.42], H1(T) 318.05 [313.51, 322.52], H2(T) 319.99 [315.04, 324.58], F(T) 327.89 [321.58, 335.96];
the Neptune statuses are those of 1 October.

**30 September 1845, S4's own date.** Neptune 323.47°; S4 from its elements 320.80°, the sheet's
printed "320° 30′ nearly" compared as printed. H1(G) 313.37 [308.98, 317.04]: S4 outside (−2.0),
S4 printed marginal (−1.9), Neptune outside (−2.8). H1(T) 316.11 [311.60, 320.57]: S4 marginal
(−1.1), S4 printed inside (−1.0), Neptune marginal (−1.6). F(G) 322.24 [314.84, 338.26]: all inside.

**The Challis nights** (geocentric longitude): the medians differ from the 1 October ones by the
planet's and the Earth's motion only; the statuses against Neptune (327.3° on 29 July, 326.9° on
12 August, 325.8° on 29 September) are those of 1 October for every run and target.

### § 5.2 Heliocentric true and mean longitude of date, 1 October 1846

Neptune: true 327.01°, osculating mean 328.23°. Targets, true then mean: S4 324.13 / 323.20; S5
328.75 / 325.10; S6 329.47 / 323.02; LV-June 324.62 / 324.62; LV-Aug 326.03 / 318.37; S7 mean
315.33 (mean longitude only).

| run | true: median [68 %] | statuses (Neptune, S4, S5, S6, LV-June, LV-Aug) | mean: median [68 %] | statuses (Neptune, S4, S5, S6, S7, LV-June, LV-Aug) |
|---|---|---|---|---|
| H1(G) | 316.50 [312.17, 320.12] | out, out, out, out, out, out | 312.22 [309.02, 314.95] | out, out, out, out, **marginal (−1.1)**, out, out |
| H2(G) | 317.87 [314.37, 322.10] | out, marg, out, out, marg, out | 314.10 [310.46, 317.15] | out, out, out, out, **inside (−0.4)**, out, marg |
| F(G) | 325.89 [318.23, 343.32] | all inside | 331.27 [321.51, 8.74] | in, in, in, in, marg (+1.6), in, marg |
| H1(T) | 319.19 [314.74, 323.56] | marg, marg, out, out, marg, marg | 312.18 [307.11, 315.57] | out, out, out, out, **inside (−0.9)**, out, marg |
| H2(T) | 321.10 [316.27, 325.59] | marg, in, marg, out, in, marg | 313.23 [308.13, 317.93] | out, out, out, out, **inside (−0.4)**, out, marg |
| F(T) | 329.04 [322.87, 336.95] | all inside | 343.07 [324.05, 354.94] | in, marg, in, marg, marg, in, marg |

The true-longitude statuses match § 5.1's target for target except H2(T)'s S6 (outside here,
marginal geocentrically). The mean-longitude intervals of the free-a runs are wide (F(G)'s 95 %
interval spans 307°–24°) because e and ϖ, hence the mean anomaly, are loosely constrained there;
they are reported, not read. Under the held distances the mean longitude is 3–5 half-widths below
every printed mean longitude of Adams's (S4 323° 12′, S5 325° 6′, S6 323° 2′). **S7's 315° 20′,
the 2 September 1846 extrapolation, falls inside the fixed-a mean-longitude intervals** of H2(G),
H1(T) and H2(T) and is marginal in H1(G)'s, 0.4–1.1 half-widths from their medians: description
only, S7 enters no statement.

### § 5.3 RA and Dec on the Challis nights (mean equator and equinox of date)

The posteriors' 68 % declination widths are 357′–582′ (6–10°) against the 9′ zones: the
posterior's spread is described, no verdict on the sweeps is drawn (plan § 5.3). On 12 August 1846,
Neptune 21h 57.4m, −13.04°:

| run | RA median [68 %] | Dec median [68 %] | Dec width 68 % |
|---|---|---|---|
| H1(G) | 21h 17.3m [20h 57.8m, 21h 33.6m] | −18.41 [−21.33, −14.89] | 386′ |
| H2(G) | 21h 21.4m [21h 05.2m, 21h 40.8m] | −16.36 [−20.18, −13.09] | 426′ |
| F(G) | 21h 52.6m [21h 19.5m, 22h 52.9m] | −14.71 [−18.82, −9.12] | 582′ |
| H1(T) | 21h 28.5m [21h 09.8m, 21h 48.4m] | −17.84 [−21.09, −14.18] | 414′ |
| H2(T) | 21h 37.1m [21h 15.0m, 21h 57.7m] | −17.48 [−20.39, −14.43] | 357′ |
| F(T) | 22h 09.4m [21h 45.5m, 22h 42.5m] | −15.44 [−18.97, −10.95] | 482′ |

On 29 and 30 July, 4 August and 29 September the medians move with the planet (Neptune 21h 58.8m,
−12.91° on 29 July; 21h 52.8m, −13.45° on 29 September) and the picture is the same: the held-a
medians lie 20–40 minutes of RA west and 3–5° south of the planet, Adams's 325° ephemeris
(21h 49m, −13° 10′ in step 2's table) between the two.

### § 5.4 Mass, e, ϖ, i; a for F

| run | mass M⊕: median [68 %] [95 %] | mass M☉ | vs printed mass (half-widths; status) | e median; one-sided 68 / 95 % upper | ϖ (Uranus frame) [68 %] | i rel. Uranus [68 %] (prior median 7.07°) | a (AU) |
|---|---|---|---|---|---|---|---|
| H1(G) | 40.3 [33.0, 50.0] [25.7, 62.3] | 1.21e-4 | S5 1.656e-4: −1.5, marginal | 0.091; 0.10 / 0.15 | 287 [263, 317] | 3.6 [1.5, 6.8] | 38.4 held |
| H2(G) | 37.0 [27.8, 43.9] [23.2, 49.4] | 1.11e-4 | S6 1.500e-4: −1.9, outside | 0.081; 0.09 / 0.15 | 288 [257, 319] | 3.5 [1.4, 6.9] | 37.25 held |
| F(G) | 26.0 [12.1, 36.0] [7.7, 48.2] | 7.8e-5 | — | 0.084; 0.13 / 0.29 | 5 [314, 55] | 3.7 [1.4, 6.6] | 32.35 [23.12, 34.74] [21.88, 40.37] |
| H1(T) | 35.1 [26.7, 44.5] [21.2, 54.1] | 1.05e-4 | S5: −2.1, outside | 0.107; 0.13 / 0.20 | 275 [249, 304] | 4.2 [1.6, 7.5] | 38.4 held |
| H2(T) | 30.1 [24.1, 41.5] [20.5, 57.8] | 9.1e-5 | S6: −1.7, marginal | 0.106; 0.13 / 0.18 | 267 [245, 317] | 4.2 [1.6, 7.4] | 37.25 held |
| F(T) | 24.9 [15.2, 33.6] [9.5, 40.8] | 7.5e-5 | — | 0.127; 0.20 / 0.26 | 32 [291, 41] | 4.6 [1.6, 7.4] | 30.12 [28.13, 33.52] [24.90, 38.09] |

Adams's masses are 55.1 M⊕ (S5) and 50.0 M⊕ (S6); Neptune's 17.1 M⊕ at 30.07 AU. F(G)'s marginal
in a is multimodal (figure): a narrow peak near 22 AU (cell [20, 25], share 0.17) and a broad one at
31–35 AU; F(T)'s peaks at 29.5 AU with the 20–25 AU cell at 0.03. The inclination posteriors sit
below the prior's median in every run, as the headline's did.

### § 5.5 Evidence, cells, best fits

Gaussian (constant dropped; comparable across H1, H2, F, which share data and error model):

| run | ln Z (± 1.8 σ) | cells: ln Z_k, share | −2 ln L_min (278 meas.) | best-fit group offsets, arcsec (ancient RA/Dec; Greenwich; Paris 1835–40) |
|---|---|---|---|---|
| H1(G) | −232.16 ± 0.34 | e < 0.12: −231.39, 0.864 ± 0.066; e ≥ 0.12: −233.64, 0.136 | 390.13 | +5.07 / −0.44; +1.62 / +0.65; +2.30 / +1.90 |
| H2(G) | −232.00 ± 0.33 | −231.26, 0.832 ± 0.077; −233.27, 0.168 | 389.42 | +3.89 / −0.64; +1.15 / +0.82; +1.63 / +2.53 |
| F(G) | −232.44 ± 0.26 | [20, 25]: −232.64, 0.167 ± 0.059; [25, 28]: −232.75, 0.076; [28, 60] e < 0.12: −231.58, 0.659 ± 0.099; e ≥ 0.12: −233.89, 0.098 | 382.25 | +1.76 / −1.09; +1.87 / +0.57; +2.41 / +1.55 |

Student-t (ν = 4, Q):

| run | ln Z | cells: ln Z_k, share | Q_min | offsets |
|---|---|---|---|---|
| H1(T) | −683.46 ± 0.28 | −683.01, 0.628 ± 0.129; −683.94, 0.371 | 1294.81 | +0.56 / −0.88; +2.43 / +0.64; +1.51 / +1.99 |
| H2(T) | −683.31 ± 0.28 | −682.94, 0.581 ± 0.133; −683.67, 0.419 | 1294.16 | +0.65 / −1.09; +2.16 / +0.75; +1.79 / +2.24 |
| F(T) | −682.77 ± 0.24 | [20, 25]: −684.84, 0.026 ± 0.009; [25, 28]: −682.60, 0.122; [28, 60] e < 0.12: −682.35, 0.422; e ≥ 0.12: −682.73, 0.430 | 1288.57 | −0.24 / −1.34; +2.36 / +0.69; +1.56 / +2.56 |

ln B(H1 : H2) = −0.16 ± 0.47 (G), −0.15 ± 0.40 (T): nothing between the two hypotheses. The held-a
evidences are also indistinguishable from F's (−232.16 and −232.00 against −232.44 with errors of
0.3; the plan does not form ln B(held : F) and none is formed here, § 5.8). The Paris 1801–1827 offsets
are 0 by construction (empty group).

### § 5.6 The test rows (Paris 1841–1845, 47 measurements)

Posterior-predictive O − C per measurement (median and 68 % over 6000 samples; per-measurement 68 %
widths ≈ 3″), summarised as the mean of the medians per year, RA then Dec, arcsec:

| run | 1841 (6 RA, 5 Dec) | 1842 (7, 7) | 1843 (3, 2) | 1844 (8, 5) | 1845 (2, 2) |
|---|---|---|---|---|---|
| H1(G) | +2.9 / +3.3 | +2.3 / +4.1 | +0.8 / +3.9 | +2.0 / +2.8 | −0.5 / +2.5 |
| H2(G) | +2.9 / +3.0 | +2.4 / +3.8 | +1.0 / +3.6 | +2.3 / +2.5 | −0.1 / +2.3 |
| F(G) | +2.6 / +3.2 | +2.1 / +4.0 | +0.5 / +3.9 | +1.7 / +2.8 | −1.0 / +2.5 |
| H1(T) | +2.4 / +3.2 | +2.0 / +4.1 | +0.8 / +4.1 | +2.2 / +3.1 | +0.2 / +3.0 |
| H2(T) | +2.4 / +3.3 | +2.1 / +4.2 | +1.0 / +4.3 | +2.5 / +3.3 | +0.6 / +3.3 |
| F(T) | +2.3 / +3.5 | +1.7 / +4.4 | +0.2 / +4.3 | +1.5 / +3.3 | −1.1 / +3.1 |

Adams's § 53 errors of the corrected tables at the 1843, 1844 and 1845 normal places: +7″.11, +8″.79,
+12″.40 (Hypothesis I); +5″.77, +7″.05, +10″.18 (Hypothesis II). Here the residuals stay at the level
of the group offset (≈ 2″, the fitted Paris offset) in every run and do not grow with time. The two
quantities are not like for like (his are errors in longitude of Bouvard's tables with his
perturbation applied, at normal places; these are RA/Dec residuals of the full fit with group offsets
not subtracted), so the comparison is description.

### § 5.7 Against the headline (M6c, all Table 77)

Galle-night geocentric longitude, headline 325.81 [323.91, 327.39], 68 % width 3.48°:

| run | median shift (°) | in headline 68 % half-widths | 68 % width ratio |
|---|---|---|---|
| F(G) | −1.11 | −0.64 | 7.33 |
| F(T) | +2.08 | +1.19 | 4.14 |

No ln B against the headline (the data differ). Dropping Flamsteed 1690, Paris 1801–1827 and the
rows after 1840 widens the Galle-night interval 4–7 times and moves the median by about one headline
half-width.

### § 5.8 H1 and H2 against F

Δ ln L_max (held − free) and F's posterior density in a (±0.5 AU bins; the fullest 1 AU bin is at
33.5 AU for F(G), 29.5 AU for F(T)):

| | Δ ln L_max, G | Δ ln L_max, T | F(G): fraction in bin; relative to max bin | F(T): fraction; relative |
|---|---|---|---|---|
| H1, 38.4 AU | −3.94 | −3.12 | 0.020 (4380 of 219 694); 0.13 | 0.011 (2373 of 207 744); 0.06 |
| H2, 37.25 AU | −3.58 | −2.80 | 0.037 (8099); 0.24 | 0.023 (4703); 0.12 |
| Neptune, 30.07 AU | — | — | 0.057 (12 602); 0.37 | 0.177 (36 700); 0.92 |

### § 6 The reading

Read on the Gaussian runs on 1 October 1846, the Student-t runs beside them.

**Statement 1, did Adams's data locate the planet?** **Yes.** Neptune is inside F(G)'s 68 % interval,
0.1 half-widths from the median. Under his assumed distance: **no** for both; Neptune is outside
H1(G) (−2.9 half-widths) and H2(G) (−2.1). Student-t: yes for F(T) (+0.3); marginal under H1(T)
(−1.7) and H2(T) (−1.3).

**Statement 2, were his communicated solutions representative of his data under his assumption?**
**No, for all three.** S4 and S5 are outside H1(G)'s intervals (−2.1 and −3.4 half-widths on the true
longitude, −4.0 and −4.7 on the mean) and S6 is outside H2(G)'s (−2.7 true, −2.9 mean); true and
mean longitude agree in every case. Student-t: S5 and S6 outside; S4 marginal on the true longitude
(−1.1) and outside on the mean (−3.3), the one target where the two quantities differ.

**Statement 3, did the data prefer his assumed distance?** **Disfavoured**, both hypotheses, both
error models: Δ ln L_max −3.94 (H1) and −3.58 (H2) under Gaussian errors, −3.12 and −2.80 under
Student-t, all below the −2 threshold; F(G)'s density at his distances is 13 % and 24 % of its
maximum. Amendment 3's weakly-disfavoured band is not entered. **Beside the evidence:** Δ ln L_max
disfavours his distances while the evidences of H1, H2 and F are indistinguishable (Gaussian
−232.16, −232.00, −232.44 with errors of 0.3; Student-t −683.46, −683.31, −682.77). The best fit at
a free a is better by 3.6–3.9 units of ln L, but the free-a prior spreads that gain over a volume of
a, so the marginal likelihoods tie; the data prefer his distances less than they prefer the best a,
and no more than that is said.

**Interpretation-table row.** Row 3 applies under both error models: his data located the planet, but
his solutions did not represent its posterior under his own assumption. The gap is therefore in the
solution method or the data handling and not in the data, **and the candidates explicitly include
this reconstruction's approximations A1–A5** (his ancient weighting, the 1831 normal place and the
Paris stand-in for Cambridge, normal places and opposition selection instead of individual
observations, his reductions, Bessel 1823) as much as his method. **The direction of the gap is
fixed:** Adams's solutions lie closer to Neptune (2–3° beyond it) than the reconstructed data allow at
his distances (9–10° short of it), so at a = 38.4 or 37.25 AU his point estimates are nearer the truth
than the posterior of data set D, not farther from it. **The gap is not yet attributable to his
method**: separating the data handling from the method is the business of a further step (3b, to be
pre-registered), which reproduces his solution from his printed normal places and equations of
condition. Neither amendment-4 case arises: S5 and S6 agree (both outside), and no marginal reading
enters the Gaussian row. S7's 315° 20′ is noted as description only: it falls inside the fixed-a
mean-longitude intervals (§ 5.2), where his computed solutions do not.

**Le Verrier, comparison only.** June (325° at 1847.0, printed, P) and August (`LE_VERRIER_1846`)
are both inside F(G) and F(T) (0.1 and 0.0 half-widths from F(G)'s median), outside H1(G), and for
H2(G) June marginal (−1.6) and August outside (−1.9).

**Gaussian against Student-t on 1 October 1846.** Ten disagreements, all in the held-a runs and all
one step softer under Student-t (outside → marginal for Neptune, S4, LV-June, LV-Aug under H1 and
for S5, S6, LV-Aug, Neptune under H2; marginal → inside for S4 and LV-June under H2); none for F,
none in the opposite sense. The Gaussian reading stands as pre-registered; under Student-t statement
1's "under his assumed distance" reads partly and statement 2 keeps S5 and S6 outside.

### Figure

`P1B_adams_fit.png`: the geocentric longitude of 1 October 1846 for H1, H2 and F, Gaussian solid and
Student-t dotted (M10's two hues), with Neptune, S4, S5, S6, LV-June and LV-Aug marked; below, F's
marginal in a with 38.4, 37.25 and Neptune's 30.07 AU marked; S7 is named in the caption and not drawn
as a solution. Three fixes to the plotting code after the first draw (longitudes unwrapped about their
circular mean so F's tail past 360° does not stretch the axis; target labels staggered; the legend
moved off the labels), hence the three further summarize runs.

### What this step established, in one paragraph

With Adams's data as reconstructed here, the planet's longitude is recovered when the distance is
left free, but at the distances he assumed the data place it 9–10° short of Neptune, while his own
solutions placed it 2–3° beyond; his three communicated solutions lie outside the fixed-distance
posteriors of these data, and the data disfavour his distances by Δ ln L_max of 3–4 while the
evidences tie. Whether the disagreement between his solutions and these posteriors comes from his
method or from the difference between his data handling and this reconstruction (A1–A5) is not
decided by this step.
