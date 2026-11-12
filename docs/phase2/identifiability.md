# Phase 2: synthetic identifiability tool

12 November 2026. Built and tested only. Running it as an experiment needs a plan, and it
reads no observation of Mercury.

`fork.phase2.identifiability` generates Mercury's longitude under each competing model
(`synthetic`) and reports how the models' predictions differ over a set of epochs
(`compare`). Each model gets the one parameter value whose analytic rate gives the same
perihelion advance, 42.98″ per century by default (`matched_parameters`):

| model | parameter matched | analytic rate used |
|---|---|---|
| `gr` | scale of the 1PN term | 6πμ/(c²a(1−e²)) per orbit |
| `j2` | the Sun's J2, spin axis from the context | dω/dt + dΩ/dt |
| `ring` | mass of a ring at 0.2 AU in the ecliptic | orbit-averaged Gauss rate |
| `planet` | mass of a point planet at 0.2 AU | Laplace–Lagrange, lowest order in e |

`compare` integrates each model and the Newtonian baseline from one initial state and
gives, relative to the baseline: the rates of ϖ, Ω and i (″ per century) and of e; the RMS
of the periodic part of the mean longitude; the RMS longitude difference; and the RMS
longitude difference for each pair of models.

## Example output

The CI fixture's 39 epochs, 1800 to 2000, initial state at 1900.0. This shows the shape of
the report. It is not a result: the epochs are not the observations', no noise is added
and no element is refitted, so every figure is an upper bound on what could separate the
models.

| model | parameter | ϖ rate | Ω rate | i rate | periodic part of λ |
|---|---|---|---|---|---|
| `gr` | 1.000 | 42.98″ | 0.00″ | 0.00″ | 0.008″ |
| `j2` | 3.40 × 10⁻⁴ | 43.03″ | −3.45″ | −2.50″ | 0.010″ |
| `ring` | 1.83 × 10⁻⁷ M☉ | 39.59″ | −42.18″ | −0.19″ | 0.006″ |
| `planet` | 2.22 × 10⁻⁷ M☉ | 48.22″ | −51.30″ | −0.23″ | 1.00″ |

Rates per century. The ring and the planet miss the target by 8 % and 12 %: both lie in the
ecliptic, so they also turn Mercury's inclined orbit about the ecliptic pole, and the
planet's rate formula ignores e = 0.2. The matching is by formula on purpose; a plan that
needs exact matching should match on the integrated rate.

## Limits

- No refit of Mercury's elements, no noise model, no time-scale error (ΔT).
- Heliocentric longitude only. Transit contact times and geocentric places are not built.
- One ring radius, one planet distance and phase per call; a plan must say which to scan.

## Reproduce

```
uv run python -m fork.phase2.identifiability
uv run pytest -q tests/test_phase2_identifiability.py
```
