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

## Partitioned validation, round 1 (three a slices): gate **not passed**

Slices [20, 25], [25, 28], [28, 60] AU, same settings, two rng seeds per slice (4 workers,
73–79 min per run). Script: `scripts/m5b_partitioned_validation.py` (at the commit that adds
this section); numbers: `M5b_partitioned_validation.json`.

| slice (a, AU) | ln Z_k rng 1 | ln Z_k rng 2 | Δ / σ | share rng 1 | share rng 2 |
|---|---|---|---|---|---|
| [20, 25] | −76.08 ± 0.23 | −76.09 ± 0.23 | 0.03 | 0.015 | 0.021 |
| [25, 28] | −73.86 ± 0.22 | −73.82 ± 0.22 | 0.13 | 0.072 | 0.105 |
| [28, 60] | −73.23 ± 0.22 | −73.60 ± 0.22 | **1.20** | 0.912 | 0.874 |
| combined | −73.50 ± 0.20 | −73.83 ± 0.20 | 1.17 | | |

- **Gate 1** fails on [28, 60] (1.2σ); shares agree.
- **Gate 2** (vs unpartitioned −73.05 / −73.31) fails: Δ −0.19 to −0.79, up to 2.7σ. The
  reference itself is unreliable. Split by slice, the unpartitioned runs give mode B zero
  mass in both seeds and [25, 28] 1.3 % vs 39 %, i.e. implied ln Z_k −75.1 vs −72.0, where
  the partitioned runs give −73.86 and −73.82.
- **Gate 3**: [28, 60] matches chain A on both seeds. [20, 25] matches chain B on rng 1; on
  rng 2 the a median is 0.55σ off (tolerance 0.5σ), and the a width is 1.68× chain B's: a
  tail toward 25 AU (98th percentile 23.8 vs 23.2) that the single-mode chain does not
  reach. The longitude median passes (0.41σ).
- **Gate 4**, 1846 longitude, median [68 %]: combined 329.1° [325.7, 332.4] (rng 1),
  328.6° [323.9, 333.0] (rng 2). Per slice in the JSON.

Cause: [28, 60] holds a second branch. Besides mode A (a ≈ 30, e ≈ 0.03) there is an
eccentric branch at a ≥ 33 (mean e 0.16–0.26), the e ≈ 0.3 branch of the profile fits. The
seeds weigh it 3 % vs 18 %; the a < 33 part has ln Z −73.26 vs −73.80. This is the same
mode-allocation scatter that broke the unpartitioned runs, now inside one slice, and it is
not in dynesty's quoted error.

Decision: split [28, 60] by eccentricity rather than by a, into e ∈ [0, 0.12) and
[0.12, 0.3], since every alternative branch found so far is eccentric. Gate 2 becomes
agreement of combined ln Z and all cell shares across seeds; gate 3 for [20, 25] compares
chain B with the slice's e > 0.25 samples.

## Partitioned validation, round 2 (four cells): gate **not passed** (shares, ln Z scatter)

Cells (`nested.PARTITION_CELLS`): [20, 25] and [25, 28] AU (round-1 runs reused), and
[28, 60] AU × e < 0.12 and × e ≥ 0.12 (new; 3 workers each, 149–157 min while sharing the
CPU). Cell [20, 25] has a third seed. Numbers: `M5b_partitioned_validation_cells.json`.

| cell | weight | ln Z_k rng 1 | ln Z_k rng 2 | Δ/σ | share rng 1 | share rng 2 | Δ/σ |
|---|---|---|---|---|---|---|---|
| a [20, 25] | 0.203 | −76.08 ± 0.23 | −76.09 ± 0.23 | 0.03 | 0.016 | 0.014 | 0.29 |
| a [25, 28] | 0.103 | −73.86 ± 0.22 | −73.82 ± 0.22 | 0.13 | 0.073 | 0.069 | 0.17 |
| [28, 60] × e < 0.12 | 0.277 | −72.62 ± 0.22 | −72.38 ± 0.22 | 0.79 | 0.678 | 0.781 | **1.36** |
| [28, 60] × e ≥ 0.12 | 0.416 | −74.10 ± 0.22 | −74.53 ± 0.22 | **1.39** | 0.233 | 0.136 | **1.53** |
| combined | | −73.52 ± 0.16 | −73.41 ± 0.18 | 0.44 | | | |

Cell [20, 25] rng 3: ln Z_k **−75.08 ± 0.23**, 1.0 above rng 1 and 2 (three-seed SD 0.58 vs
quoted 0.23).

- **Gate 1** fails on both e cells (eccentric ln Z_k 1.39σ; both shares).
- **Gate 2**: combined ln Z agrees (0.44σ); the shares do not (the two e cells).
- **Gate 3** passes. Cell [20, 25], e > 0.25: rng 3 matches chain B (all medians ≤ 0.29σ) and
  rng 1, not rng 2, so by the agreed rule rng 2 is an outlier draw. The circular cell matches
  chain A on both seeds (≤ 0.43σ; 97 % of chain A lies in the cell).
- **Gate 4**, 1846 longitude, median [68 %]: combined 329.2° [322.2, 333.1] (rng 1), 328.1°
  [322.5, 332.2] (rng 2). Circular cell 329.8° / 328.4°, eccentric cell 321.1° / 321.3°,
  [20, 25] 338.1° / 334.5°, [25, 28] 337.4° / 337.3°.

Cell [20, 25], three seeds: shapes agree well at the medians (seed-to-seed SD 0.09–0.39 of
chain B's σ) and worse in the upper tail of a, u_a, u_h and mass (SD 0.70–0.99σ at the 84th
percentile). The ln Z_k jump has the same posterior shape behind it. At fixed log L, rng 3's
ln X runs ~1.0 above rng 1 and 2 from log L ≈ −60 on (−29.3 vs −30.1/−29.6 at −60, −36.8 vs
−37.8 at −40), while all three agree to 0.2 down to log L −80. So the prior-volume estimates
part where the [20, 25] contour fragments, consistent with one or more pockets being lost
early (which biases ln Z low) in rng 1 and 2, or kept too long in rng 3.

Overall: on these cells dynesty's quoted ln Z error understates the run-to-run scatter. The
pooled single-run SD over the four cells is ≈ 0.40 against a quoted 0.22 (factor ≈ 1.8).
