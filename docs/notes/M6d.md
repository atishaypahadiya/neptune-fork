# M6d: the Le Verrier shape model against free distance (pre-registration)

Written 2026-11-13, before any code of this plan computed a quantity its reading uses.
Draft: `docs/notes/V13_drafts.md` ("M6d"). Everything above the line at the end is frozen
by the commit that adds this file; amendments and results go below it, dated.

## 1 Question

Spec, M6d: for each model, the probability it puts on landing within 0.87° of Neptune on
Galle's night, P(|λ − 325.873°| ≤ 0.87°), which is Le Verrier's actual miss of 52′. Compare
the Le Verrier shape model (the distance held at his value) with the free-distance model:
the ratio of the two probabilities, and the Bayes factor between the models.

## 2 Inputs

| input | SHA-256 |
|---|---|
| `data/historical/leverrier_1846.csv` (all 529 measurements) | 2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4 |
| `data/historical/m6_error_model.json` | 0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277 |
| `docs/notes/M8_u_ref.json` | 4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5 |
| `src/fork/model/perturbative.py` | 092f673b9352f9d2b0ae709f2a4d3f258dd09806245b767fd83b86e8ecccd68c |
| `src/fork/model/linearized.py` | 3e02e78bdd9ee74fb96a6666f192188c355113addc92f43e376729f0b616d15e |
| `src/fork/inference/reweight.py` | d0c5fbff6066ebaa531bc24668fa807a8f7bf0aa96cf77eead5c07aca9cd021b |
| `site/data/m8b.json` (the free-distance fit) | 355104b24966ce216f53ff603931a060a4a5e738a996a683db5fac5c5e962027 |
| `site/data/m6c_replacement.json` (the reported posterior) | d4d1b5c862169984b6152c38b57809da36404e8f84043c5fdf658c05b3f42e32 |

- **Free distance, no rerun.** The headline at the higher-resolution settings
  (`m8b.json`, `headline/merged`): ln Z −361.523 ± 0.078, median 325.735°. Its probability
  within 0.87° is computed here from its cached draws (`data/cache/m8/runs_m8b/`), by the
  code that made `m8b.json`. Before use, the recomputed median must reproduce 325.735° to
  0.001°; if it does not, the item stops.
- The reported posterior's published probability, 0.3655 (`m6c_replacement.json`), is
  printed beside it, not used in the ratio: it comes from a different sampler path.
- No Neptune ephemeris enters a fit. 325.873° is used in reporting only.

## 3 Method

- **Shape model:** the headline's model and priors with the distance held: ln a fixed, e ≤
  0.3, everything else free (`v13.shape_space`). Two cells, e < 0.12 and e 0.12 to 0.3.
- **Held value, two runs:**
  - **primary, a = 38.4 AU.** The spec says "fixed a, as in M5", and M5's Le Verrier mode
    holds a at 38.4 AU (Bode's law, his first paper). It is also the value held in the
    Adams–Challis study.
  - **secondary, a = 36.154 AU,** from his elements of 31 August 1846. Run and reported
    the same way. The reading of § 5 is made on the primary; the secondary is stated
    beside it, and a different band for it is reported as a dependence on the held value.
- **Data, error model, track:** the headline's (all rows, modern track, frozen σ_g, τ_g).
- **Likelihood, sampler, reweighting:** as the fit it is compared with. The fast likelihood
  about the headline's reference orbit, reweighted exactly to mode (b); dynesty static,
  rslice with 28 slices, dlogz 0.1, nlive 4000 per cell, four seeds per cell; 6,000 draws
  per run and per merge. rng `[2026, 16, 1846, f, k + 1, s]`, f = 1 for 38.4 AU, 2 for
  36.154 AU, k the cell, s = 1..4; draws, chunks and fits seeded as in
  `M6c_variants.md` § 3. That is 16 runs. Both evidences are then integrals of the same
  likelihood over the same Uranus box, so their difference is a Bayes factor.
- **The ratio:** R = P_shape / P_free on 6,000 equal-weight draws of each fit, with a 68 %
  bootstrap interval (`v13.peirce_ratio`, 2,000 resamples, rng `[2026, 16, 1846, f, 0, 0,
  17]`). Also P_shape for each single-seed fit.
- **The Bayes factor:** ln B(shape : free) = ln Z_shape − ln Z_free, error the two in
  quadrature (`v13.log_bayes_factor`), each ln Z error being the sampler's × 1.8 plus the
  reweighting's.

## 4 Gates and tests, in this order

G0 to G4 and T1 of `M6c_variants.md` § 4, applied to each shape fit (G0 on a synthetic data
set on all rows, seed 1). In T1 the ln Z range of the four single-seed fits must also be
0.5 or less.

## 5 Reading, fixed now

Made only if the gates and T1 pass for the primary fit; otherwise the numbers are reported
as not validated or not repeatable and no band is read.

**Bayes factor.** With L = ln B(shape : free) and σ its error:

| L | wording |
|---|---|
| \|L\| < 1 | the data do not prefer either model |
| 1 ≤ \|L\| < 3 | a weak preference (for the shape model if L > 0, for free distance if L < 0) |
| 3 ≤ \|L\| < 5 | a strong preference |
| \|L\| ≥ 5 | a decisive preference |

A band is stated only if L − 2σ and L + 2σ both lie in it; otherwise "between" the two
bands they reach.

**Ratio.** If the 68 % interval of R contains 1: "holding the distance at his value did not
change the probability of a prediction as close as his". If it lies above 1: "made it more
probable, by a factor R"; below 1: "less probable, by a factor 1/R". The four single-seed
values of P_shape are listed; if their range exceeds 0.1 the ratio is reported as not
repeatable and not read.

**What "Peirce's accident" becomes.** One sentence, whichever way it falls: P_free is the
chance that the data alone, with no assumption about the distance, give a prediction as
good as Le Verrier's; P_shape is the same chance under his assumption; R compares them. No
statement about luck is made beyond these three numbers.

## 6 What this cannot show

- The free-distance posterior depends on the distance prior (log-uniform, 20 to 60 AU) and
  the Bayes factor depends on it directly: a wider prior on a lowers Z_free. The number
  answers "his assumption against this prior", nothing more general.
- The shape model here holds a only. Le Verrier's published elements also fix e and ϖ;
  that fully held model is M5's and is not rerun here.
- 0.87° is a threshold chosen after the event (his actual miss).
- Gaussian errors only. The Student-t twin is not run: its free-distance evidence comes
  from a single run whose repeatability was never tested, and the fast likelihood does not
  cover Student-t errors, so a ratio against it could not be read.
- The 1846 reduction track is not covered (no specification).

## 7 Choices in the draft, and how each is resolved

| open in the draft | resolution |
|---|---|
| held value, 36.154 or 38.4 AU | both; 38.4 AU primary because the spec says "as in M5" |
| seed namespace | tag 16, § 3 |
| Student-t twin | not run, § 6 |
| reading rule | § 5 |
| nlive 1000, one full-likelihood run per cell | four seeds at nlive 4000 with the fast likelihood and exact reweighting, to match the fit it is compared with and to test repeatability |
| free-distance comparison | the higher-resolution headline fit, not the first single-run fit |

## 8 Cost

16 runs at about 0.35 CPU-h, reweighting, G0 and G4: about 8 CPU-h, under 1 h on 12
workers. Runs start after the data-choice fits.

## 9 Layout

- Code: `src/fork/export/v14_fits.py` (shared with the data-choice fits), tests
  `tests/test_v14_fits.py`, branch `v14-runs`.
- Cache: `data/cache/nested/v14_fits/<fit>/`. Results: `docs/notes/M6d.json`.

---
