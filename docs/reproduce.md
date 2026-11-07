# Reproducing each milestone

All commands run from the repository root after `uv sync`. Runtimes come from the compute
ledger (`site/data/ledger.json`, built from `data/cache/compute_ledger.jsonl`) where a job
recorded one. Rows marked *est.* were backfilled from saved results and logs, not timed live.
Ledger CPU-hours add up across jobs. Wall times are per job: many jobs ran side by side on one
14-core Apple Silicon machine, so summed wall time overstates the elapsed time. Jobs the ledger
does not cover took seconds to a few minutes and are marked "not in ledger".

**Ledger total so far: 435 CPU-hours** (353 timed live, 82 estimated), 257 million likelihood
calls. One likelihood call is one N-body integration over 1690–1846, so that is 40 billion
simulated planet-years. A6c is running now and is not in the total yet.

## Setup

```sh
uv sync
uv run python -m fork.data.horizons --fetch   # JPL Horizons → data/cache/ (network, once)
uv run pytest -m "not slow and not network"   # the CI suite: 240 tests, 72 s under A6c load
uv run python examples/toy_fit.py             # one-cell demo; see "Toy example" below
```

Most tests and the toy example run offline from the committed Horizons fixture
(`tests/data/horizons_m4.json`). The full exports need `data/cache/` from `--fetch`.

## M1–M2: frames and Kepler solver

```sh
uv run pytest tests/test_frames.py tests/test_kepler.py tests/test_horizons.py
```
Seconds, not in ledger. The `network` tests in `test_horizons.py` check the cache against live
Horizons. Add `-m "not network"` to skip them.

## M3: integrators (A1)

```sh
uv run pytest tests/test_integrators.py
uv run python -m fork.export.to_json --only wh_convergence.json   # site/data/wh_convergence.json
uv run python docs/notes/scripts/m3_energy_plot.py                # docs/notes/M3_energy.png
```
Seconds to a minute, not in ledger.

## M4: forward model and synthetic data (A2, A2b)

```sh
uv run pytest tests/test_forward.py tests/test_residuals.py tests/test_synthetic.py
uv run python docs/notes/scripts/m4_residuals_plot.py             # docs/notes/M4_residuals.png
```
Seconds, not in ledger.

## M5a: point fits and profile likelihood

| command | output | ledger wall (12 workers) |
|---|---|---|
| `uv run python -m fork.export.to_json --only m5a_fits.json` | `site/data/m5a_fits.json` | 80 s *est.* |
| `uv run python docs/notes/scripts/m5a_diagnostics.py` | `docs/notes/M5a_diagnostics.json` | 55 min *est.* |
| `uv run python -m fork.export.to_json --only m5a_profile.json` | `site/data/m5a_profile.json` | 3.4 h *est.* |
| `uv run python docs/notes/scripts/m5a_profile_plot.py` | `docs/notes/M5a_profile.png` | not in ledger |
| `uv run python docs/notes/scripts/m5a_residuals_plot.py` | `docs/notes/M5a_residuals.png` | not in ledger |

The profile export seeds starts from the M5b posteriors, so run it after M5b.

## M5b: partitioned nested sampling (A3)

Validation on seed 0 (`docs/notes/M5b_dynesty_validation.md`):

```sh
uv run python docs/notes/scripts/m5b_dynesty_validation.py --n-workers 6
uv run python docs/notes/scripts/m5b_partitioned_validation.py --run K RNG --n-workers N
uv run python docs/notes/scripts/m5b_partitioned_validation.py --summarize
```
13 runs, **79 CPU-h** *est.*, 49–165 min wall per run.

Full run, then the summary:

```sh
uv run python -m fork.export.m5b plan                          # jobs and their state
uv run python -m fork.export.m5b run --n-workers 12 --per-job 4
uv run python -m fork.export.m5b summarize                     # site/data/m5b_posteriors.json
```

| part | jobs | ledger CPU-h | wall per job (4 workers) |
|---|---|---|---|
| free a, e ≤ 0.3, seeds 0–9 × 4 cells | 37 (3 reused from validation) | 200 | 1.1–1.6 h; cell 0 at nlive 2000: 2.5–3.2 h |
| free a, e ≤ 0.1, seeds 0–9 × 3 cells | 30 | 117 | 1.0–1.9 h |
| Le Verrier shape, seeds 0–4 (+ one check rerun) | 6 | 16 | 0.7–0.8 h |
| noise-free truth, 4 cells | 4 | 21 | 1.1–2.9 h |
| `summarize` (4 passes) | 4 | 2.5 *est.* | 8–41 min on 12 workers |

Total for M5b including validation: **≈ 435 CPU-h**, about 36 h of wall time on 12 cores if
nothing else is running. The run checkpoints every 5 min and resumes from where it stopped.
Finished jobs are cached, so `run` can be restarted at any time. Progress is in
`data/cache/nested/m5b/progress.log`.

Then:

```sh
uv run pytest tests/test_recovery.py               # A3a (slow) and A3b on the committed summary
uv run python docs/notes/scripts/m5b_seed5_diagnostic.py   # docs/notes/M5b_seed5_lon.png, ~5 min
uv run python docs/notes/scripts/m5b_noise_free_degeneracy.py
uv run python -m fork.export.ledger                # site/data/ledger.json
```

## M6a: transcription and loader

```sh
uv run python -m fork.data.transcription    # diff pass 1 vs pass 2, per-row arithmetic check
uv run python -m fork.data.observations     # data/historical/leverrier_1846.csv
uv run pytest tests/test_transcription.py tests/test_historical.py tests/test_spotcheck.py
uv run --with pymupdf python -I -m fork.data.spotcheck   # docs/review/spotcheck.html
```
Seconds, not in ledger. The spot-check page needs the Gallica scan, which is not committed.
Its URL and SHA-256 are in `data/historical/sources/MANIFEST.md`.

## M6b: observation model and error calibration

```sh
uv run python -m fork.data.calibration      # data/historical/m6_error_model.json
uv run pytest tests/test_apparent.py tests/test_observation.py tests/test_calibration.py
```
Seconds, not in ledger. Calibration needs Horizons states for Uranus and Earth at the 279
epochs, which `--fetch` caches.

## A6c: synthetic coverage with the M6 model (running)

```sh
uv run python -m fork.export.a6c plan
uv run python -m fork.export.a6c run --n-workers 12 --per-job 4
uv run python -m fork.export.a6c summarize  # site/data/a6c_coverage.json
```
There are 40 jobs (10 seeds × 4 cells), with the same nlive as M5b. The ledger has no row yet.
The job count matches M5b's e ≤ 0.3 run (about 200 CPU-h), but each likelihood call costs more
here, because it computes the apparent place at 279 epochs. Treat 200 CPU-h as a lower bound
until the ledger records the real figure.

## Toy example (not a result)

```sh
uv run python examples/toy_fit.py
uv run pytest tests/test_toy_example.py     # capped at 15 iterations, a few seconds
```
This runs one cell (a 28–60 AU, e < 0.12) on synthetic seed 0 at nlive 20, single process,
offline. One run took 4.6 min on one core (674 iterations, 41k likelihood calls) while A6c
kept the machine at load ~19, so expect less on an idle machine. It shows the pipeline working.
Its interval and ln Z are not M5b results. That run's 68 % interval was about 5° wide, half the
width of the real M5b intervals, and its ln Z_k error was ±3. At nlive 20 the sampler
under-explores the cell.
