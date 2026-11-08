# neptune-fork

[![test](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/test.yml/badge.svg)](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/test.yml)
[![pages](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/pages.yml/badge.svg)](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/pages.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24044791.svg)](https://doi.org/10.5281/zenodo.24044791)

This project rebuilds Le Verrier's 1846 prediction of Neptune from Uranus's observed motion, from
first principles: its own symplectic N-body integrator, and calibrated Bayesian posteriors for
where the unseen planet should be in place of a single best fit. Phase 2 will
generalise it into a "fork engine" that weighs hidden mass, a modified law of gravity and plain
data error as explanations for an orbital anomaly.

![Ten horizontal 68% intervals for the perturber's 1846 longitude, one per synthetic dataset, against a vertical line at true Neptune; two dashed intervals miss it](docs/notes/M5b_coverage.png)

*Coverage test (A3b) on synthetic data.* Each row is one noise realisation. It shows the
posterior median and central 68 % interval for where the inferred planet sits on 1846-09-23,
Galle's discovery night. True Neptune (325.87°) lies inside 8 of the 10 intervals. For
calibrated 68 % intervals that count gives a two-sided binomial p = 0.52, and the intervals
are 7–11° wide.

There are two separate noise-free checks:
- **Point fit (M5a, least squares):** recovers Neptune's 1846 longitude to 0.07°. This shows
  the model and fitter are unbiased.
- **Posterior (M5b):** wider, with its median at 322.05° against the true 325.87°, a −0.75σ
  offset. The shift runs along the mass–distance degeneracy, where the likelihood is flat.

Every number here comes from synthetic data. No real observations have been fitted yet.

## Status

Phase 1, **M6 in progress**. Milestones M0–M5b are complete:

- frames;
- the Kepler solver;
- the integrators, matching REBOUND to 2e-9″;
- the forward model;
- point fits and the profile likelihood;
- partitioned nested sampling with the A3 calibration tests.

For M6, Le Verrier's 279 observations are transcribed and double-checked (`data/historical/`),
and the apparent-place observation model and error calibration are built. The A6c coverage test
on synthetic data passed: the truth is inside the 68 % interval for 8 of 10 seeds (two-sided
binomial p = 0.52), with an rms offset of 0.87σ and a median 68 % width of 5.1°. That is
narrower than M5b's 9.1°, mostly because Table 77 has more measurements and smaller calibrated
errors than M5b assumed ([M6 notes](docs/notes/M6.md)). The real-data fit (M6c) is next.

Data: [Le Verrier Table 77, Zenodo, doi:10.5281/zenodo.24044792](https://doi.org/10.5281/zenodo.24044792)

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

The ledger (`site/data/ledger.json`) records **713 CPU-hours** to date: 631 timed live, plus 82
estimated from logs. That covers 428 million likelihood evaluations, each one N-body
integration over 1690–1846. Nested sampling accounts for nearly all of it: M5b, and A6c with
277 CPU-hours.

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
