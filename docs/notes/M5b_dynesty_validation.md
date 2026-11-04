# M5b: dynesty validation on seed 0 (free a, e ≤ 0.3): gate **not passed**

Setup: shared prior (`inference/priors.py`, Uranus box ±50 scaled σ, ln volume −24.66),
`bound='multi'`, `sample='rslice'` (12 slices), nlive 1000, Δln Z < 0.1, perturber ϖ and λ₀
periodic. Two rng seeds, 6 workers. Script: `scripts/m5b_dynesty_validation.py`; numbers:
`M5b_dynesty_validation.json`.
References: emcee chain A (started at the 24-start point fit, sat at a ≈ 30 AU) and chain B
(started on the 343° branch, a ≈ 22 AU, e ≈ 0.28).

| | rng 1 | rng 2 | gate |
|---|---|---|---|
| ln Z | −73.046 ± 0.222 | −73.313 ± 0.222 | Δ 0.27 < 0.31 combined: **pass** |
| iterations / calls / wall | 38 760 / 2.49 M / 56 min | 39 329 / 2.58 M / 49 min | |
| best χ² in run | 73.21 | 72.60 | (mode B peak 71.22, point fit 71.91) |
| mass in mode B (a < 25 AU) | 0.000 | 0.000 | **fail**: mode B not recovered |
| mass at a 24–28 AU (region C) | 0.013 | 0.388 | **fail**: not reproducible |
| mass at e > 0.2 | 0.044 | 0.351 | |
| a ≥ 28: median vs emcee A | ≤ 0.29σ | ≤ 0.32σ | pass |
| a ≥ 28: 68 % width / emcee A's | 1.1–2.3 | 1.0–1.5 | broader than emcee |

What the runs show:

- **Mode B is never found.** Each run drew ~2000 points at a < 24 AU, all with χ² ≥ 127. The
  pocket is thin, with e pressed against the 0.3 bound. A Laplace estimate from chain B's and
  chain A's sample covariances (flat prior in q) gives ln(M_B / M_A) ≈ −4.6, so B probably
  holds only ~1 % of the mass. Missing it costs little mass, but the gate requires it.
- **A third region, C**, sits around the point fit: a ≈ 26–28 AU, e ≈ 0.17–0.28. This is the
  branch the 340° constrained fits found (a 26.9, e 0.20, χ² 71.96). Run 1 gives it 1 % of the
  mass and run 2 gives it 39 %. Emcee chain A started there, left, and kept 2.5 % at a < 28 AU.
  Splitting the samples at a = 25 AU hides this, so the "mode A" comparison fails for run 2
  (e +1.1σ, longitude +1.0σ) only because region C is mixed in.
- ln Z still agrees between runs. Which pockets survive depends on the luck of where live
  points sit when the likelihood contour splits. With nlive 1000 and rslice, the sampler finds
  and weighs these small branches stochastically.

Conclusion: these settings are not adequate for M5b. Next step proposed in the review
message (partitioned nested sampling).
