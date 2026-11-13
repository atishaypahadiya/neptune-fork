# What of this repository's work was done before

Written 13 November 2026 from the nine reading notes in this directory. The claim numbers below
refer to the numbered claims in those notes, each of which is tagged there as read (with page
and form) or not read.

## What was and was not read

This table can only say what is in the pages I read. It cannot say that anything is absent
from the literature.

| source | read | not read |
|---|---|---|
| Lyttleton 1958 | all of it, pp. 551–559 (text layer; p. 557 on the image) | his 1960 *Vistas* paper with the numerical work |
| Sampson 1904 | pp. 143–154 as OCR text | the catalogue, pp. 154–170, beyond headings; the plates |
| Brookes 1970 | pp. 67, 78–80 as OCR; the rest through the repository's notes of 2026-11-11 | pp. 68–77 for this note; his 1969 thesis |
| Kollerstrom 2006 | Part 1 p. 1 and passages on pp. 7, 10, 11, 18–20, some notes; Part 2 pp. 349–352, 361, 365–367, parts of 357 and 364, some notes | about two thirds of each part, including most of the account of Challis's search |
| discoverneptune (Whitmore 2026) | five pages in full, two for prose and results, five by search | five pages; all figures; the source repository (404) |
| CoFI example | theory pages, `meta.yml`, the deterministic notebook's prose and printed results, part of the Bayesian notebook | the Bayesian results; figures |
| Lai, Lam and Young 1990 | the abstract | the article |
| Smart 1946–47 | one paragraph of the *Nature* article | both printings |
| Grosser 1962 | nothing | the book |
| Roseveare 1982 | nothing | the book |

Not among the nine items, and not read: the KTH thesis of 2018 that `README.md` names;
D. Rawlins, MNRAS 147, 177 (1970); E. W. Brown, MNRAS 92, 93 (1931); the histories by Standage
and by Sheehan and others.

## Table

Three grades only: **done before**, **partly done before**, **not found in what was read**.
The last means exactly that and no more.

| what this repository does | where | nearest prior work, and the claim it rests on | verdict |
|---|---|---|---|
| Its own symplectic integrator (Wisdom–Holman) | `src/fork/integrators/wisdom_holman.py` | discoverneptune integrates with REBOUND's WHFast, a symplectic integrator from a library (`discoverneptune_and_cofi.md` A3). The CoFI example has its own integrator, written with numba, but it is fourth-order Runge–Kutta (B3). | **partly done before.** A symplectic integrator was used; an own-written integrator was used; the two together I did not find. This is an implementation choice and not a result. |
| A full posterior, not a point fit | `docs/notes/M5.md`, `M6.md` | discoverneptune samples a posterior with emcee, with convergence diagnostics, calibration and predictive checks (A7, A9). CoFI samples one with emcee (B7). | **done before**, as a method, on modern ephemeris positions. On the historical observations themselves I did not find it: discoverneptune says its data are a JPL-derived proxy and that it holds no sourced transcription of the normal places (A2); CoFI fits Horizons vectors (B2). |
| The mass–distance degeneracy made explicit | `docs/notes/M9.md`, `M10.md` | discoverneptune traces the ridge by pinning the distance from 26 to 44 AU and refitting, gives the m/a² scaling, and shows the longitude barely moving (A6, A8). Brookes solves Adams's equations at seven distances from 31.86 to 39.56 AU and concludes the longitude is insensitive to a′ between 30 and 40 AU after 1830 (`brookes_1970.md` 6, 7). Sampson says the eccentricity and apse "compensate the excessive value attributed to the mean distance" (`sampson_1904.md` 6, 13). Lyttleton scans the distance with a circular orbit and finds a minimum of the residual near a′/a = 1.6 (`lyttleton_1958.md` 7). | **done before.** What I did not find is a profile likelihood over distance on the historical table under two error models; and Lyttleton's minimum is a published result that points to a determinable distance, which the notes here should address. |
| Bayesian evidence over a partitioned prior | `docs/notes/M5b_dynesty_validation.md`, `M6.md` | No evidence, Bayes factor or nested sampling in discoverneptune by a word search of its pages (A13), nor in the CoFI files read (B9). The nearest thing is a statement of multimodality: Lai, Lam and Young's inverse problem "admits a fairly good solution that places Neptune in a diametrically opposite position" (`lai_lam_young_1990.md` 5). | **not found in what was read.** The body of Lai et al. was not read, and that a second solution exists is in their abstract. |
| An error model with group offsets | `docs/notes/M6.md` | discoverneptune fits a stationary AR(1) discrepancy with free scale and correlation, because independent residuals are autocorrelated (A7). CoFI re-estimates a per-coordinate σ from the residuals (B4, B6). Neither has offsets by group of observations. | **partly done before.** A structured error model with free parameters exists; offsets by observing group I did not find. |
| The pre-registered check against Lalande's star of 1795, which is a check within the data span, not an extrapolation | `docs/notes/M6e.md` | discoverneptune runs a hindcast, fitting a prefix of 1781–1846 and scoring the held-out years, on its synthetic benchmark (A9). Kollerstrom records that Adams used Lalande's position in 1847 to fix the eccentricity of the orbit after the discovery (`kollerstrom_2006.md` 9), which is a use of the precovery and not a test of a prediction made without it. | **not found in what was read.** A held-out check as such was done before, on synthetic data; a check of a fit to the historical table against a recorded precovery I did not find, nor any pre-registration. |
| Reproduction of Adams's solutions from his memoir | `docs/notes/P1B_adams_reproduction.md` | Brookes re-derives Hypothesis I and II from the memoir's equations and tabulates the difference from Adams's print (`brookes_1970.md` 3, 4; from the repository's notes, not re-read on his pp. 68–77). Sampson tabulates the six solutions and when each was right (`sampson_1904.md` 4, 12). | **done before**, by Brookes; the repository already says so. |
| Three additions to that reproduction: placing the numerical difference in one normal equation; the χ² profile over α and its pole; the split between his data and his method | same note, § 0 and R10 | Brookes attributes the difference to "a small numerical error occurring in one of the equations" without naming it in the paper; his thesis is where it would be named (`brookes_1970.md` 3, and its last section). He gives point solutions over α, not a profile (6). Sampson states in words that the timing of the agreement was "an accidental circumstance" (`sampson_1904.md` 14). | **not found in what was read**, with two open doors: Brookes's 1969 thesis and Smart 1947, both unread, either of which may contain the first addition. |
| The study of Challis's search as a detection probability | `docs/notes/P1B_challis_detection.md` | Kollerstrom gives the documentary facts (start 29 July, Neptune seen 4 and 12 August, the "disc" entry) and no probability in the pages read (`kollerstrom_2006.md` 2, 7, 14). discoverneptune tells it as narrative: the comparison stopped at star 39, Neptune was entry 49 (A12). Lyttleton mentions the 30° by 10° belt Airy proposed (`lyttleton_1958.md` 5). | **not found in what was read.** Most of Kollerstrom's account of the search, and all of Grosser, were not read. |
| The attribution to Bouvard's tables | `docs/notes/P1B_bouvard_attribution.md` | Sampson records that Adams's data to 1821 were Bouvard's residuals and that Adams corrected Bouvard's tables in points due to Bessel and Hansen (`sampson_1904.md` 2, and p. 145). Brookes finds the solutions "particularly sensitive" to the discrepancies used (`brookes_1970.md` 8). Neither asks how much of Adams's shift the tables carry. | **not found in what was read.** Smart 1947, said to be the fullest account from the manuscripts, was not read. |

## What follows for the wording elsewhere

I changed nothing outside this directory. Three places would need a second look by whoever
edits them:

1. `README.md`, "What this project does differently": "full posteriors" and "the mass–distance
   degeneracy made explicit" are not differences from discoverneptune as it stood on
   2026-11-13. The differences that survive this reading are the data (the historical table,
   transcribed), the evidence over a partitioned prior, and the error model with group
   offsets.
2. `README.md` and `docs/notes/M6.md`, the benchmark "~2.5–2.9″": the figure I found is 2.5″
   (2.52″), and it is a misfit to JPL-derived annual longitudes, not to historical
   observations (`discoverneptune_and_cofi.md`, "Disagreement with the repository").
3. `docs/notes/P1B_sources.md` § 5: Kollerstrom's Part 2 is obtainable from ADS, and the
   appendices are in Part 2 (`kollerstrom_2006.md`). That note is append-only, so the
   correction belongs in a dated section at its end.

## Needs a browser / library

| item | what would have to be read | where |
|---|---|---|
| Grosser 1962 | the chapters on the computations and on Challis's search; any table of elements; any recomputation | a library |
| Smart 1947, Occasional Notes of the RAS 2 (no. 11) | any recomputation of Adams's solutions; p. 55 on the two hypotheses; the account of the 1845 papers | a library holding the Occasional Notes; print only as far as I found |
| Smart 1946, *Nature* 158, 648–652 | the article | subscription or library |
| Lai, Lam and Young 1990 | the data, method and residuals of the inverse problem; the opposite solution | subscription to the *American Journal of Physics*, or a library |
| Roseveare 1982 | any table setting the hypotheses for Mercury against one another | an archive.org loan or a library |
| Lyttleton 1960, *Vistas in Astronomy* 3, 25–46 | the seven solutions and their residuals | a library |
| Brookes 1969, London Ph.D. thesis | which of Adams's equations carries the numerical error | University of London theses |
| discoverneptune source repository | its licence and the provenance of the bundled data | GitHub, signed in; the address answered 404 without a login |

Reachable by script and simply not yet read: the rest of Kollerstrom's two parts and of
Sampson's catalogue (both cached); Rawlins 1970 and Brown 1931 (probably on ADS, not fetched).
