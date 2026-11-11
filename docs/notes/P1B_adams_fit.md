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
