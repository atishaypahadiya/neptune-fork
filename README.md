# neptune-fork

[![test](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/test.yml/badge.svg)](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/test.yml)
[![pages](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/pages.yml/badge.svg)](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/pages.yml)
[![DOI, all versions](https://zenodo.org/badge/DOI/10.5281/zenodo.24044791.svg)](https://doi.org/10.5281/zenodo.24044791)

Refitting Le Verrier's prediction of Neptune from his observations of Uranus.

## Abstract

We refit the 279 observations of Uranus (529 measured coordinates, 1690–1845) that Le Verrier
published in 1846 (his Table 77), with a perturbed-orbit model that is never given Neptune: our
own Wisdom–Holman integrator, an apparent-place observation model, and partitioned nested
sampling run under a plan committed before the run. On the night Galle found Neptune,
23 September 1846, we find the unseen planet's geocentric longitude (ecliptic of date) at
**325.8°** (68 %: 323.6°–327.8°; 95 %: 320.2°–330.4°, provisional). Neptune was at 325.87°
(JPL), 0.03 half-widths from the median. The same posterior gives a mass of 17.3 M⊕ (68 %:
16.3–18.4; Neptune 17.15), a semi-major axis of 30.2 AU (29.6–30.9), osculating in 1781, and an
eccentricity below 0.028 (68 %, one-sided). Le Verrier's published elements, carried to that
night, give a geocentric longitude of 324.97°, inside our 68 % range. Phase 2 will extend the code to weigh three
explanations of an orbital anomaly against each other: a hidden mass, a modified law of
gravity, and error in the data.

The site (`site/`) presents the result as a short paper with plates; every number on it is read
from `site/data/`.

## Calibration on synthetic data

![Ten horizontal 68% intervals for the perturber's 1846 longitude, one per synthetic dataset, against a vertical line at true Neptune; two dashed intervals miss it](docs/notes/M5b_coverage.png)

*Coverage test (A3b) on synthetic data.* Each row is one noise realisation: the posterior
median and central 68 % interval for the inferred planet's longitude on 1846-09-23. True Neptune
(325.87°) lies inside 8 of the 10 intervals; for calibrated 68 % intervals that count gives a
two-sided binomial p = 0.52. The intervals are 7–11° wide.

Two noise-free checks:
- **Point fit (M5a, least squares):** recovers Neptune's 1846 longitude to 0.07°; the model and
  fitter are unbiased.
- **Posterior (M5b):** wider, with its median at 322.05° against the true 325.87°, a −0.75σ
  offset along the mass–distance degeneracy, where the likelihood is flat.

## Status

Phase 1. Milestones M0–M5b are complete:

- frames;
- the Kepler solver;
- the integrators, which match REBOUND to 2e-9″;
- the forward model;
- point fits and the profile likelihood;
- partitioned nested sampling with the A3 calibration tests.

The pre-registered real-data fit (M6c) is done and reported, and the A6 acceptance tests are
reported ([A6 report](docs/notes/A6.md)). M6d and M6e are deferred to v1.1. M8, a faster
likelihood that decides the provisional 95 % interval, is running. M7 is the public site.

**Data and model (M6).** Le Verrier's 279 observations are transcribed and double-checked
(`data/historical/`), and the apparent-place observation model and the error calibration are
built. The A6c coverage test on synthetic data passed: the truth is inside the 68 % interval for
8 of 10 seeds (two-sided binomial p = 0.52), with an rms offset of 0.87σ and a median 68 % width
of 5.1°. That is narrower than M5b's 9.1°, mostly because Table 77 has more measurements and
smaller calibrated errors than M5b assumed ([M6 notes](docs/notes/M6.md)).

**Result (M6c).** On Table 77 we find the planet at 325.8° on Galle's night (68 %:
323.6°–327.8°), with a sampler systematic of ±0.13° on the median. JPL's Neptune, 325.87°, is
0.03 half-widths from the median. Le Verrier's prediction, carried from his published elements
to that night (324.97°, geocentric, ecliptic of date), is also inside the interval. The posterior puts 37 % of its probability
within Le Verrier's 52′ miss of the true position (the Peirce probability). These figures are
derived after the fact from the pre-registered replacement posterior described below.

**Robustness.** JPL's Neptune is inside the 68 % interval in the headline fit and in all three
pre-registered sensitivity runs (e ≤ 0.5; error levels estimated from the data alone; Student-t
errors). The medians agree to within 0.4 headline half-widths. The evidence favours the
Student-t errors (ln B = +13.7), but the Gaussian error model was fixed by the
pre-registration and is kept; S3 shows the result is insensitive to it (the median moves 0.37
headline half-widths).

**Repeatability.** Repeated runs of the main sampling cell disagreed on the width by far more
than the sampler's own error estimate. Under a pre-registered repeatability study, two runs with
4× the live points agreed on the width, and their merge replaced that cell. This widened the
first published interval (323.9°–327.4°) by about 20 %, to the one above
([M6 notes](docs/notes/M6.md#m6c-cell-2-repeatability-results-2026-11-08)). The 68 % interval
is the headline. The 95 % interval (320.2°–330.4°) is provisional: it widened from 6.6° to 10.2°
because it depends on the small weights of three cells that are still single runs. M8, a faster
likelihood that makes repeated runs of every cell affordable, is specified in the
[Phase 1 spec](docs/PHASE1_SPEC.md#m8-fast-linearized-likelihood).

**Data availability.** Le Verrier's Table 77, transcribed, is on Zenodo (CC0).
- Concept DOI [10.5281/zenodo.24044791](https://doi.org/10.5281/zenodo.24044791) (the badge
  above) covers all versions; cite it for the dataset in general.
- Version DOI [10.5281/zenodo.24044792](https://doi.org/10.5281/zenodo.24044792) is v1.0.0;
  cite it for exact reproduction.

## Limitations

- A6b (the reduction check against JPL): Comparison against per-group precision not made: no precision figure was fixed before the data were seen. Medians reported. Details in the [A6 report](docs/notes/A6.md).
- A6a (the transcription check): 0 of 28 rows in error; 95 % upper bound on the per-row error rate 10.1 % (exact one-sided Clopper–Pearson, 0 of 28); the remaining 251 rows are unchecked.
- The observation model does not include diurnal parallax (up to 0.5″ for Uranus over 1781–1845).
- Error model: the evidence favours Student-t errors (S3, ln B = +13.7); the Gaussian error model
  was fixed by the pre-registration and is kept. S3 moves the median 0.37 headline half-widths.

## Quickstart

```sh
uv sync
uv run pytest -m "not slow and not network"   # the CI suite, offline
uv run python examples/toy_fit.py             # a few minutes, one core, offline
```

`examples/toy_fit.py` runs the real pipeline on one of the four prior cells at nlive 20. It is
a demo, not a result. Fitting anything else needs JPL Horizons states:
`uv run python -m fork.data.horizons --fetch` caches them in `data/cache/`. Commands and
runtimes for each milestone are in [docs/reproduce.md](docs/reproduce.md).

## Project structure

```
src/fork/
  integrators/  Kepler solver, Wisdom–Holman (democratic heliocentric), leapfrog, RK4
  data/         Horizons cache, frames, synthetic data, Table 77 transcription and loader,
                error calibration
  model/        bodies, forward model, hidden-mass model, apparent place, M6 observation model
  inference/    priors, least-squares fits, profile likelihood, emcee, dynesty (partitioned)
  export/       site/data/*.json writers, M5b and A6c run schedulers, compute ledger
examples/       toy_fit.py: one-cell demo
tests/          pytest; `slow` and `network` markers; offline Horizons fixtures in tests/data/
data/historical/  Le Verrier 1846 Table n° 77: both transcription passes, resolution log,
                  misprints, loader output, error model
docs/           method, reproduction commands, Phase 1 spec, milestone notes, figures
site/           static site (GitHub Pages) and its data
```

## Documentation

- [docs/method.md](docs/method.md): the method in one page. It covers frames, the integrator
  and its REBOUND validation, the observation models, the priors, partitioned nested sampling,
  the calibration tests and the ×1.8 evidence-error inflation.
- [docs/reproduce.md](docs/reproduce.md): commands and ledger runtimes for each milestone.
- [docs/PHASE1_SPEC.md](docs/PHASE1_SPEC.md): milestones and acceptance tests, with every
  amendment.
- Milestone notes, with the numbers achieved and the surprises:
  - [M0](docs/notes/M0.md): scaffold
  - [M1](docs/notes/M1.md): frames and data access
  - [M2](docs/notes/M2.md): Kepler solver
  - [M3](docs/notes/M3.md): integrators, A1
  - [M4](docs/notes/M4.md): forward model, A2
  - [M5](docs/notes/M5.md): point fits, profile, posteriors, A3
  - [M5b validation](docs/notes/M5b_dynesty_validation.md)
  - [M6](docs/notes/M6.md): historical data, observation model

## Compute

The ledger (`site/data/ledger.json`) records **934 CPU-hours** to date: 851 timed live, plus 83
estimated from logs. That covers 562 million likelihood evaluations, each one N-body
integration over 1690–1846. Nested sampling accounts for nearly all of it: M5b, and A6c with
277 CPU-hours. The M6c headline fit took 28.5, its sensitivity runs S1–S3 89, and the cell-2
repeatability study 104.

## Prior art

This is well-trodden ground. We credit and benchmark against:

- Jonathan Whitmore, *discoverneptune* (2026): <https://jbwhit.github.io/discoverneptune/>.
  A REBOUND-based reconstruction reaching ~2.5–2.9″ RMS, and our benchmark for the historical
  fit.
- Lai, H. M., Lam, C. C., & Young, K. (1990): forward and inverse treatment of Uranus's
  perturbations by Neptune.
- CoFI (Common Framework for Inference): the example "Finding Neptune using Uranus".
- KTH thesis (2018): Neptune from Uranus by Gauss–Newton, with an RK4 integrator.

What this project does differently:

- its own Wisdom–Holman integrator, with REBOUND used only in tests;
- full posteriors from partitioned nested sampling, with evidence;
- the mass–distance degeneracy made explicit;
- a forward-model interface built for comparing models in Phase 2.

## Citing and license

Citation metadata is in [CITATION.cff](CITATION.cff). The licenses are split ([LICENSE](LICENSE)):

- **Code:** MIT. This covers everything outside `data/historical/`.
- **Data:** CC0 1.0 (public-domain dedication) for everything in `data/historical/`: the
  Table n° 77 transcriptions, the resolution log, the misprint table, the loader output and the
  error model. Citing this project when you reuse them is appreciated but not required.

Le Verrier's 1846 memoir itself is in the public domain. Its scans are listed, with URLs and
hashes, in `data/historical/sources/MANIFEST.md`.
