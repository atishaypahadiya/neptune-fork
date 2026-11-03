# neptune-fork

Rebuild Le Verrier's 1846 prediction of Neptune from Uranus's observed motion, from first principles,
then (Phase 2) generalize it into a **fork engine**: given an orbital anomaly, decide between
*hidden mass*, *modified law*, or *mundane error* via Bayesian evidence.

Current phase: **Phase 1**. Full milestone spec: `docs/PHASE1_SPEC.md`. Work milestones in order.

## Prior art (credit it, don't duplicate it)
- Jonathan Whitmore, *discoverneptune* (2026): https://jbwhit.github.io/discoverneptune/ — REBOUND-based
  reconstruction, ~2.5–2.9″ RMS. Our benchmark for A4. Check its license before reusing any bundled data;
  cite it in README either way.
- Lai et al. (1990), forward + inverse treatment of Uranus perturbations.
- CoFI example "Finding Neptune using Uranus"; KTH thesis (2018, Gauss–Newton + RK4).
Our differences: own symplectic integrator, full posterior (not point fit), explicit mass–distance
degeneracy, and a code structure built for Phase 2 model comparison.

## Non-negotiable rules
1. **Own integrator in `src/fork/integrators/`.** REBOUND is a *test dependency only*, used to validate.
   Never import rebound from `src/`.
2. **Uranus's in-plane elements (a, e, ϖ, λ₀) are always free parameters** in any fit. Its i and Ω
   are held at their epoch (1781) values: longitude-only data can't constrain them (a second-order,
   ~10″ effect) and, freed, they absorb the perturber signal. The leak this introduces is pinned
   by a test (Neptune's effect on i, Ω over 1781–1846 is 0.017″ in longitude). Never hold Uranus's
   modern Horizons state (or any in-plane element) fixed when Neptune is removed — that leaks the
   answer.
3. **Synthetic before historical.** No code touches historical observations until test A3
   (recover a planted planet) passes.
4. **Frames are explicit.** Every angle carries its frame in the variable name or type
   (`lon_ecl_j2000_rad`, `lon_ecl_of_date_rad`). Precession (~50″/yr) dwarfs the ~130″ signal;
   a frame bug will fake a planet.
5. **Units:** AU, days, solar masses. G = k², k = 0.01720209895 (Gaussian). Angles internally in
   radians; arcseconds only at I/O and plotting boundaries.
6. **Wrap-safe residuals:** always `((obs - model + π) % 2π) - π`.
7. **No paid services.** Everything runs locally or on free GitHub Actions / Pages.
   API credits are reserved for Phase 3; do not add LLM calls in Phase 1.

## Stack
- Python ≥3.11, managed with `uv`. numpy, scipy, numba, emcee, astropy, astroquery, matplotlib.
- Tests: pytest (+ rebound, hypothesis). Lint/format: ruff. Types: annotate public functions.
- Site: static, `site/` — D3 + Observable Plot, reads `site/data/*.json`. Deployed to GitHub Pages.

## Commands
```
uv sync                       # install
uv run pytest -q              # all tests
uv run pytest -m "not slow"   # fast loop
uv run ruff check . && uv run ruff format .
uv run python -m fork.data.horizons --fetch   # cache JPL vectors to data/cache/ (gitignored)
uv run python -m fork.export.to_json          # regenerate site/data/
```

## Layout
```
src/fork/
  integrators/   leapfrog.py, kepler.py, wisdom_holman.py, rk4.py (baseline for comparison)
  data/          horizons.py, observations.py, frames.py, synthetic.py
  model/         bodies.py, forward.py (θ → Uranus longitudes), hidden_mass.py, likelihood.py
  inference/     fit.py (multistart least squares), mcmc.py (emcee)
  export/        to_json.py
tests/           test_kepler.py, test_integrators.py (A1), test_residuals.py (A2),
                 test_recovery.py (A3, slow), test_historical.py (A4, slow)
site/            index.html, js/, data/
docs/            PHASE1_SPEC.md, notes/
```
Phase 2 will add `model/mundane.py` and `model/modified_law.py` beside `hidden_mass.py`; keep the
forward-model interface generic: `predict(theta, t_jd) -> lon_rad` with a model-specific θ.

## Conventions
- Hot loops (`integrators/*`, `model/forward.py`) are `numba.njit`, pure-array, no Python objects.
- Time: Julian Date (TDB) floats internally; convert with astropy only at boundaries.
- Randomness: every stochastic function takes an explicit `rng: np.random.Generator`.
- Mark tests over ~5 s with `@pytest.mark.slow`.
- Each milestone ends with: tests green, a short note in `docs/notes/Mx.md` (what was built, numbers
  achieved, surprises), and one commit per milestone minimum.

## When stuck
Prefer adding a diagnostic plot to `docs/notes/` over guessing. If a numeric target in the spec is
missed, report the actual number and the suspected cause; don't loosen the test silently.
