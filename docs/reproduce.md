# Reproducing each milestone

All commands run from the repository root after `uv sync`. Runtimes come from the compute
ledger (`site/data/ledger.json`, built from `data/cache/compute_ledger.jsonl`) where a job
recorded one. Rows marked *est.* were backfilled from saved results and logs, not timed live.
Ledger CPU-hours add up across jobs. Wall times are per job: many jobs ran side by side on one
14-core Apple Silicon machine, so summed wall time overstates the elapsed time. Jobs the ledger
does not cover took seconds to a few minutes and are marked "not in ledger".

**Ledger total so far: 1090 CPU-hours** (1007 timed live, 83 estimated), 1153 million likelihood
calls, 624 million of them N-body integrations over 1690–1846 (the rest are the fast linearized
likelihood of M8), or 96 billion simulated planet-years.

## Run identifiers

The site and the README name runs and checks descriptively. The milestone notes, the spec and
the commands below use these identifiers.

| name on the site and in the README | identifier | notes |
|---|---|---|
| the headline fit (its posterior: the headline posterior) | M6c headline; the replacement posterior, with cell 2 replaced by the R2 merge | `docs/notes/M6.md` |
| the sensitivity runs | M6c S1–S3 | `docs/notes/M6.md` |
| the eccentricity sensitivity test (eccentricity up to 0.5) | S1 (`s1_e05`) | `docs/notes/M6.md` |
| the error-level sensitivity test (error levels from the data alone) | S2 (`s2_sigma_data`) | `docs/notes/M6.md` |
| the Student-t sensitivity test | S3 (`s3_t4`) | `docs/notes/M6.md` |
| the repeatability runs | M6c cell-2 repeatability (R1, R2) | `docs/notes/M6.md` |
| the fast-likelihood validation | M8, A8b–A8f (48 runs) | `docs/notes/M8.md` |
| the fast-likelihood rerun | M8 A8f rerun | `docs/notes/M8.md` |
| the higher-resolution reruns | M8b | `docs/notes/M8.md` |
| the full-likelihood check of the eccentricity test | A8d(S1) arbiter (`m8_arbiter`) | `docs/notes/M8.md` |
| the profile likelihood in the semi-major axis under Gaussian errors (Fig. 7 left, Fig. 8) | M9 (`fork.export.m9`) | `docs/notes/M9.md` |
| the profile likelihood under Student-t errors, six data sets (Fig. 7 right, Fig. 9) | M10 (`fork.export.m10`) | `docs/notes/M10.md` |
| the 95 % stability condition | `std_95_bounds_le_0p4` (rule_95 in `m8b.json`) | `docs/notes/M8.md` |
| the convergence gates; the agreement with the headline fit; the seed-agreement condition | `gates_headline_and_single_seed_fits`; `a8b_a8c`; `a8f_headline_cells` | `docs/notes/M8.md` |
| the transcription check | A6a | `docs/notes/A6.md` |
| the reduction check against JPL | A6b | `docs/notes/A6.md` |
| the coverage test on synthetic data with the real error model | A6c | `docs/notes/A6.md` |
| the leak check (the fit never reads Neptune) | A6d | `docs/notes/A6.md` |
| the coverage test on synthetic data | A3b (M5b) | `docs/notes/M5.md` |

## Setup

```sh
uv sync
uv run python -m fork.data.horizons --fetch   # JPL Horizons → data/cache/ (network, once)
uv run pytest -m "not slow and not network"   # the CI suite: no network calls, about 1 min
uv run python examples/toy_fit.py             # one-cell demo; see "Toy example" below
```

The CI suite makes no network calls: every Horizons state it reads comes from the committed
fixtures (`tests/fixtures/horizons/`, see its README, and `tests/data/horizons_m4.json` for the
older `offline_horizons` modules), and a live fetch raises. The `network` tests, which exercise
the Horizons client, run in their own weekly workflow (`.github/workflows/network.yml`) or with
`uv run pytest -m network`. `slow` tests and the full exports read `data/cache/` from `--fetch`.
The toy example runs offline from `tests/data/horizons_m4.json`.

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

## A6c: synthetic coverage with the M6 model (passed: 8 of 10)

```sh
uv run python -m fork.export.a6c plan
uv run python -m fork.export.a6c run --n-workers 12 --per-job 4
uv run python -m fork.export.a6c summarize  # site/data/a6c_coverage.json
```
There are 40 jobs (10 seeds × 4 cells), with the same nlive as M5b. The recorded run took
**277 CPU-hours** and 171 million likelihood calls: 28.4 h wall on 12 workers (2026-11-06
10:04 to 2026-11-07 14:25), with each job on 4 workers taking 86–111 min (nlive 1000) or 202–221 min (nlive 2000). Each likelihood call
costs more than in M5b because it computes the apparent place at 279 epochs. The result is in
`site/data/a6c_coverage.json` and `docs/notes/M6.md`.

## M9: profile likelihood in the semi-major axis

```sh
uv run python -m fork.export.m9 profile --dataset all --n-workers 12       # 1.6 min wall, 0.09 CPU-h
uv run python -m fork.export.m9 profile --dataset post1781 --n-workers 12  # 2.3 min wall, 0.07 CPU-h
uv run python -m fork.export.m9 convergence --dataset all                  # iteration trace of every exact polish; seconds
uv run python -m fork.export.m9 report        # docs/notes/M9_profile.{json,png} and site/data/m9_profile.json
uv run python -m fork.export.m9 screen        # the transcription screen (docs/notes/M9_screen.json)
uv run python docs/notes/scripts/m9_mass_a.py # the after-the-fact corner plot, 23 s
```

Pre-registered in `docs/notes/M9.md` (code hash and relaunch recorded there); ledger rows
`m9/profile/*`. Resumable per grid point from `data/cache/m9/`; a run that dies loses at most one
point's 40 fits.

## M10: profile likelihood under Student-t errors

```sh
FORK_HORIZONS_CACHE_ONLY=1 uv run python -m fork.export.m10 profile --dataset all --n-workers 12   # 9.0 min wall, 0.38 CPU-h
# likewise post1781 (12.1 min, 0.68), no_flamsteed (10.0, 0.50), no_bradley (10.2, 0.49),
# no_mayer (11.1, 0.51), no_lemonnier (9.7, 0.49)
uv run python -m fork.export.m10 report       # docs/notes/M10_profile.{json,png} and site/data/m10_profile.json
uv run python docs/notes/scripts/m9_mass_a.py --variant s3_t4   # the Student-t posterior's mass–a plot, after the fact
```

Pre-registered in `docs/notes/M10.md` (amendment, code hashes and the relaunch recorded there);
ledger rows `m10/profile/*`, 3.1 CPU-h in all. `FORK_HORIZONS_CACHE_ONLY=1` makes a cache miss
raise instead of fetching. Resumable per grid point from `data/cache/m10/`.

## Toy example (not a result)

```sh
uv run python examples/toy_fit.py
uv run pytest tests/test_toy_example.py     # capped at 15 iterations, a few seconds
```
This runs one cell (a 28–60 AU, e < 0.12) on synthetic seed 0 at nlive 20, single process,
offline. One run took 3.0 min on one core (674 iterations, 41k likelihood calls) on an idle
machine (2026-11-10; 4.6 min under A6c's load on 2026-11-06). It shows the pipeline working.
Its interval and ln Z are not M5b results. That run's 68 % interval was about 5° wide, half the
width of the real M5b intervals, and its ln Z_k error was ±3. At nlive 20 the sampler
under-explores the cell.

## Site (Pages build)

```sh
uv sync --group site                                              # Playwright (uses Google Chrome)
uv run --group site python tools/prerender_site.py --out _site    # prerendered page + plate-i.png
uv run --group site python tools/check_prerender.py _site         # audit with JS off, = JS page
uv run --group site python tools/check_layout.py _site            # layout, print, interaction
```

The Pages workflow runs the same three steps and deploys `_site/`. Before a change to the site,
also run the layout check in WebKit, the engine of Safari on iPhones
(`uv run --group site playwright install webkit`, then `--engine webkit`); `--shots DIR` saves
screenshots of every plate and figure at each width. Without Google Chrome, run
`uv run --group site playwright install chromium` first. To preview, serve `_site/` (or
`site/`, which renders everything client-side) with any static server.

## Phase 1b

```
uv run python docs/notes/scripts/p1b_neptune_on_challis_dates.py   # data/historical/adams_challis/neptune_*.csv (cache only, seconds)
```

### Step 3: the Adams-variant fit (pre-registered in `docs/notes/P1B_adams_fit.md`)

```sh
uv run python -m fork.data.horizons --fetch          # once, online: Earth and Neptune at the report epochs (pre-run change (a))
uv run python tests/fixtures/build_horizons_fixture.py   # CI fixture with the same epochs
uv run pytest -q tests/test_p1b_adams.py              # data set D's counts, the held-a space, seeds, targets, reading rule
uv run python -m fork.export.p1b_adams plan           # the 16 jobs, their nlive, workers and seeds
FORK_HORIZONS_CACHE_ONLY=1 uv run python -m fork.export.p1b_adams job h1_g/c0 --n-workers 4 --maxiter 400   # smoke run, not saved
# the launch (detached: own process session, nohup, caffeinate -is); variants in § 8's order, cells in parallel on 12 workers
mkdir -p data/cache/p1b_adams
FORK_HORIZONS_CACHE_ONLY=1 nohup caffeinate -is uv run python -m fork.export.p1b_adams run --n-workers 12 \
    > data/cache/p1b_adams/run.out 2>&1 &
FORK_HORIZONS_CACHE_ONLY=1 uv run python -m fork.export.p1b_adams summarize   # docs/notes/P1B_adams_fit.{json,png}
```

Progress in `data/cache/p1b_adams/progress.log`, the driver's log in `driver.log`, one ledger row
per job (`p1b_adams/<variant>/c<k>`); resumable per job from its `.ckpt`. `FORK_HORIZONS_CACHE_ONLY=1`
makes a cache miss raise instead of fetching.

Run of 2026-11-10/11: 70.6 CPU-h over the 16 jobs (held-a variants 6.8–7.3 CPU-h each, free-a 20.7
and 21.8), 7 h 53 min wall on 12 workers, 51.9 M likelihood calls; results in
`docs/notes/P1B_adams_fit.md`. The ledger rows are in the local ledger only until the site step.

## Phase 1b step 3b, stage A (transcription)

No compute. Two transcription passes of the memoir (`docs/notes/adams_memoir_passes/`), compared by

```
python3 -I docs/notes/scripts/p1b_adams_memoir_compare.py   # run from docs/notes/adams_memoir_passes/
```

(the script expects `pass1/` and `pass2/` beside it); reconciliation in
`docs/notes/adams_memoir_reconciliation.md`, final tables in `data/historical/adams_challis/
adams_memoir_*.csv`. Hi-res renders for the decisions: see `MANIFEST.md`.

## Phase 1b step 3b, stages B–D (the reproduction)

```
uv run pytest -q tests/test_adams_1846.py          # stages B and C: his intermediates and his solution from three entry points (24 tests)
uv run pytest -q tests/test_p1b_adams_repro.py     # stage D's pieces (13 tests; one slow: the synthetic set D is noise-free)
uv run python -m fork.export.p1b_adams_repro plan  # the two D3(iii) jobs, seeds [2026, 11, 3, 7, 0, k+1], state
FORK_HORIZONS_CACHE_ONLY=1 uv run python -m fork.export.p1b_adams_repro cheap      # B, C, the location, D1, D2, D3 (i) and (ii): minutes, cached in data/cache/p1b_adams_repro/cheap.json
FORK_HORIZONS_CACHE_ONLY=1 uv run python -m fork.export.p1b_adams_repro job 0 --n-workers 4 --maxiter 300   # smoke run of a D3(iii) cell, not saved
# D3(iii), detached (own process session, nohup, caffeinate -is): step 3's H1(G) job on the synthetic set D, two cells on 6 workers each
mkdir -p data/cache/p1b_adams_repro
FORK_HORIZONS_CACHE_ONLY=1 nohup caffeinate -is uv run python -m fork.export.p1b_adams_repro d3iii --n-workers 12 \
    > data/cache/p1b_adams_repro/run.out 2>&1 &
FORK_HORIZONS_CACHE_ONLY=1 uv run python -m fork.export.p1b_adams_repro summarize   # D3(iii)'s summary, D4, the readings, docs/notes/P1B_adams_reproduction.{json,png}
```

Progress in `data/cache/p1b_adams_repro/progress.log`, the driver's log in `driver.log`, one ledger
row per nested job (`p1b_adams_repro/d3iii/c<k>`) and one for the cheap stages; the nested jobs
resume from their `.ckpt`. D4 needs step 3's cached results (`data/cache/p1b_adams/*.pkl`).
