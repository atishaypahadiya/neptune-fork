# neptune-fork

[![test](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/test.yml/badge.svg)](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/test.yml)
[![pages](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/pages.yml/badge.svg)](https://github.com/atishaypahadiya/neptune-fork/actions/workflows/pages.yml)

In 1846 Urbain Le Verrier looked at a century and a half of Uranus wandering ~2′ off its predicted
track and concluded that something unseen was pulling on it; Galle found Neptune within a degree of
his prediction on the first night of searching. This project rebuilds that inference from first
principles — our own symplectic N-body integrator, Uranus's orbit always fitted jointly with the
perturber, and a full Bayesian posterior rather than a point fit, so the famous mass–distance
degeneracy ("right direction, wrong distance") is visible instead of hidden. Phase 2 generalizes it
into a *fork engine*: given an orbital anomaly, weigh **hidden mass** vs **modified law** vs
**mundane error** by Bayesian evidence.

## Prior art

This is well-trodden ground; we credit and benchmark against:

- Jonathan Whitmore, *discoverneptune* (2026) — <https://jbwhit.github.io/discoverneptune/>.
  REBOUND-based reconstruction reaching ~2.5–2.9″ RMS; our benchmark for the historical fit.
- Lai, H. M., Lam, C. C., & Young, K. (1990), forward and inverse treatment of Uranus's
  perturbations by Neptune.
- CoFI (Common Framework for Inference) example, "Finding Neptune using Uranus".
- KTH thesis (2018), Neptune from Uranus via Gauss–Newton with an RK4 integrator.

What differs here: a home-grown Wisdom–Holman integrator (REBOUND used only to validate it), a full
posterior via MCMC, an explicit treatment of the mass–distance degeneracy, and a forward-model
interface designed for Phase 2 model comparison.

## Quick start

```sh
uv sync
uv run pytest -m "not slow"
uv run python -m fork.data.horizons --fetch   # cache JPL Horizons vectors to data/cache/
```

See `docs/PHASE1_SPEC.md` for milestones and acceptance tests.
