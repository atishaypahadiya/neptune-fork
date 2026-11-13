# Two recent computational reconstructions: discoverneptune and the CoFI example

Reading note, 13 November 2026.

## A. discoverneptune (Whitmore 2026)

### Citation

Jonathan Whitmore, "The Rediscovery of Neptune: an executable reconstruction of Le Verrier's
1846 prediction", <https://jbwhit.github.io/discoverneptune/>. Pages dated on the site: "April
3, 2026" (tutorial; mass–distance explanation) and "April 29, 2026" (historical observations).
Built with Quarto 1.10.18; the server's `Last-Modified` for the index was Tue, 18 Aug 2026.
Footer on every page: "© Jonathan Whitmore".

**Verified** as a live site with that author and title. The repository it links to was not
reachable (below), so no commit or version can be cited.

### Licence

**Not stated in anything I could reach.** The site links its source as
`https://github.com/jbwhit/discoverneptune`. On 2026-11-13 that address answered 404, as did
the GitHub API for the repository, its licence endpoint and the raw `LICENSE` path; the
author's public repository list (84 entries) has no repository of that name. The seventeen site
pages I fetched contain no licence statement, only the copyright footer. So the licence "as
stated in its repository" cannot be recorded: the repository is private, renamed or removed.
Until that changes, the conservative reading is all rights reserved, and nothing from the site
(data files, figures, code) may be reused here. `docs/release_checklist.md` already records
that none of its data are reused.

### Access

Seventeen pages fetched with curl on 2026-11-13 into `data/cache/prior_work/discoverneptune/`;
URL of each is `https://jbwhit.github.io/discoverneptune/<name>.html` (the index is the bare
address). SHA-256 of those I cite:

| page | SHA-256 |
|---|---|
| `index.html` | `3b6ebe2c589b77259b4a2504a0c03b614e9855d3dec01c41b31487ac60bf4389` |
| `tutorial.html` | `8f4c25e384bc85bc069f50f12348e20ffa1d78a8e3988c55cfb5bf5ff3940e2c` |
| `neptune_discovery.html` | `eec028e5f39b407a7a3a2f5ff6eef1ecdbf03f4d05518a673ef2f04dd2493802` |
| `validation_historical.html` | `6abd95490339ab57a9962e238f953be8db3530cb94317cf71f23fc5b874ac542` |
| `leverrier-observations.html` | `4ade8a02c1e9a9d59e20595b23b3a71dbb3a1c2ceaf84aa846e2575ab3f6fd87` |
| `explanation-mass-distance-degeneracy.html` | `2b69afe9813338f1d0d4126d6391ca8e6635b019d77c40b596eabf19f9c29415` |
| `bayesian.html` | `224181e549ca5536b0622660ea7570259a53c7462a2274976ce99621b33f0f66` |
| `validation.html` | `1dfd4560813d60b13f69a1110511ce8b948618f00d5585cb1c692d3bb03bb205` |
| `sensitivity.html` | `a957bab259bc77a00b50094f16d74e2c2824ddd51fa058c503adfdc89573f3d9` |
| `mercury_vulcan.html` | `7a233c6140f45190a37b256179b664f7661fe054a1cf47deab845611db29b0f1` |
| `leverrier-legacy.html` | `6e72143a756d913f1aebe08c73fb05b80afc363af576bbf3e2fe20d2e07483f2` |

The other six (`leverrier-foundations`, `leverrier-uranus`, `leverrier-mercury`,
`rebound_neptune_tutorial`, `howto-observation-windows`, `reproducibility`) are hashed in the
cache's `MANIFEST.md`.

### What was read

The rendered pages with tags stripped, so prose, code cells and printed outputs but no
figures. Read in full: the index, "Notebook: Neptune Discovery", "Notebook: Historical
Comparison", "The Historical Observations", "Explanation: Mass-Distance Degeneracy". Read for
prose and printed results, skipping most code: "Notebook: Bayesian Inference", "Notebook:
Validation". Searched for particular terms and read the passages found: the tutorial,
sensitivity, Mercury and Vulcan, legacy, reproducibility. Not read: foundations, the Uranus
problem, the Mercury puzzle, the REBOUND tutorial, the observation-windows how-to. No source
file was read, since the repository is unreachable.

### Claims

1. The project "turns Uranus's 1781–1846 longitude anomalies into an inverse problem". [read:
   index]
2. Data: 66 annual longitudes of Uranus, 1781–1846, median cadence 365 days, in a bundled file;
   elsewhere called "the bundled annual JPL-derived proxy longitudes". They are not historical
   observations: "the later illustrative rows and the annual JPL-derived proxy are not
   historical observations", and "the
   repository does not yet contain a sourced transcription of those normal places". The only
   sourced historical rows are the "ancient" residuals of 1690–1771 from Le Verrier's Table I,
   which are displayed and not fitted. [read: neptune_discovery, validation,
   validation_historical, leverrier-observations]
3. Forward model: Jupiter, Saturn and Uranus start from JPL state vectors at 1781-01-01;
   Neptune has four free parameters (mass, semi-major axis, eccentricity, mean longitude),
   coplanar with pericentre fixed along the reference direction; integration with REBOUND's
   WHFast at dt = 0.05. [read: neptune_discovery; tutorial lines found by search]
4. Fit: `scipy.optimize.differential_evolution` over wide bounds (a from 25 to 45 AU), then a
   polish. [read: tutorial, by search; neptune_discovery]
5. Residual level: without Neptune, RMS 31.5″ and peak to peak 131.4″; with the fitted Neptune,
   "RMS after Neptune: 2.5 arcsec" (tutorial) and "JPL RMS (historical search domain): 2.52
   arcsec" (Bayesian notebook); a coarse grid in the discovery notebook reaches 3.16″. The index
   calls these "few-arcsecond, serially correlated residuals against modern JPL state-vector
   data" and says "a broader domain improves the fit, so this is a domain-specific result
   rather than an irreducible model floor". [read: validation, tutorial, bayesian,
   neptune_discovery, index]
6. Mass–distance degeneracy, shown two ways: a ridge traced by pinning the distance at values
   from 26 to 44 AU and refitting mass, eccentricity and longitude (on the synthetic benchmark:
   mass changes 11-fold, RMS under 0.1″ to about 1.3″, longitude within about ±16°); and three
   choices of bounds on the semi-major axis, with the longitude "insensitive to the assumed
   distance range". The leading-order scaling is given as m/a² with a caveat about Laplace
   coefficients. [read: explanation-mass-distance-degeneracy]
7. Posterior: emcee, 32 walkers, in two regimes. Synthetic, four parameters, 2″ injected noise,
   σ fixed: planted mass and eccentricity fall outside their 68 % intervals. JPL comparison,
   six parameters with a stationary AR(1) discrepancy (σ and ρ free): σ 2.37″ [1.66, 3.61]
   under a flat prior, 2.33″ under a prior centred near Bode's 38.8 AU; ρ sits at its upper
   bound; "modern Neptune lies outside the displayed orbital intervals. This is evidence of
   simplified-model discrepancy, not a physical recovery of Neptune's orbit." Convergence by
   the 50τ rule, ESS and split-R̂; chains checkpointed. [read: bayesian]
8. A flat prior against a Bode-centred prior on distance is compared: the Bode prior "pulls
   the semi-major axis estimate higher". [read: bayesian]
9. Validation: a hindcast (fit a prefix of 1781 to a cutoff, score the held-out years to 1846)
   run on "the synthetic observation window"; simulation-based calibration with four IID
   replicates and two AR(1) replicates, called "lightweight and demonstrative, not
   statistically rigorous"; posterior predictive checks, where the JPL fit's lag-1 tail
   probability is 0.00. [read: validation]
10. Sensitivity: injected noise of 0, 1, 2 and 5″ with ten trials each; truncated windows;
    bounds on distance of 25–35, 35–45 and 20–50 AU; all on the synthetic benchmark. [read:
    sensitivity, passages found by search]
11. Mercury: a six-body Newtonian integration over 500 years gives about 525″ per century for
    the perihelion, REBOUNDx adds about 43″, and a table gives the mass a Vulcan would need by
    first-order secular theory. No inference is run on Mercury observations. [read:
    mercury_vulcan, passages found by search]
12. Challis and Adams are treated in narrative: Adams's note reduced to 328° 41′ against
    326° 57′; Challis compared 4 and 12 August and "stopped at star number 39 … Neptune was
    entry 49". No probability of detection is computed. The page cites Smart as "Occasional
    Notes … 2 (1947): 33–55", a third page range for that pamphlet. [read: leverrier-legacy,
    passages found by search]
13. I found no nested sampling, no evidence or Bayes factor, no Lalande or Galileo precovery,
    and no reproduction of Adams's memoir, by searching all seventeen pages for those words.
    [read: search of the stripped text; absence in a search, not a reading of every page]

### Disagreement with the repository

`README.md` describes this work as "a REBOUND-based reconstruction reaching ~2.5–2.9″ RMS, and
the benchmark for the historical fit", and `docs/notes/M6.md` compares its own 2.85″–2.91″
with "about 2.5″–2.9″". In the pages as they stand on 2026-11-13 I find 2.5″ and 2.52″; I did
not find 2.9″. More important, that figure is the misfit of a four-parameter model to 66
annual JPL-derived longitudes, which the site itself says are not historical observations. The
repository's 2.85″–2.91″ is on Le Verrier's Table 77. The two numbers are close in size and
measure different things; calling one the benchmark of the other needs that sentence beside
it. The site may have changed since the README line was written; I cannot check, the
repository being unreachable.

## B. CoFI example "Finding Neptune using Uranus"

### Citation

InLab, CoFI development team, `cofi-examples`, directory
`examples/Finding_Neptune_Inversions/` with `theory/finding_neptune_deterministic.md` and
`theory/finding_neptune_bayesian.md`, <https://github.com/inlab-geo/cofi-examples>, at commit
`d81ceb78d768537457877b836315d587d474e147` (head of `main` on 2026-11-13, last push
2026-11-07). The deterministic notebook's first heading is "Finding Neptune using Uranus";
`meta.yml` titles the example "Neptune discovery".

**Verified:** the title in the README of this repository ("Finding Neptune using Uranus") is
the notebook's heading. No authors are named in the files I read beyond the team.

### Licence

BSD 2-Clause "Simplified" License, as GitHub reports for the repository and as the file
`LICENCE` at that commit reads: "Copyright 2021-2024 InLab, CoFI development team. All rights
reserved. Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met …". [read: `LICENCE`, first lines;
GitHub API licence field]

### Access

From `https://raw.githubusercontent.com/inlab-geo/cofi-examples/d81ceb78d768537457877b836315d587d474e147/`,
2026-11-13, into `data/cache/prior_work/cofi/`:

| file | SHA-256 |
|---|---|
| `LICENCE` | `2746a7c23b7e279cc745a6da55894fb6ebc6a68054c0994e0cf841c41b17f3a7` |
| `examples/Finding_Neptune_Inversions/meta.yml` | `c2294d7cf001ada4fc3f0966883fb94fa9d368a5bae5a32e9142cafad0628368` |
| `…/finding_neptune_via_deterministic_inv.ipynb` | `e6ff341f66876e8e50d6e9bb1e6a9427cea8eb89f447b72eaddf8905c0ded7ed` |
| `…/finding_neptune_via_bayesian_inv.ipynb` | `b25e87a7847410be494b6076084b09c8aae0254ebf8d522d81c1e95aabdd83d5` |
| `…/neptune_deterministic_methods.py` | `c2cbe03dae7c3d8af1a23afa9d03505c67ef79767d27c06eba44cf2817646185` |
| `…/neptune_bayesian_methods.py` | `a7487b898ee1d28b1280fe7995162b70a21a84bdceedaabd8e8cd98812d70d03` |
| `…/neptune_mcmc_helpers.py` | `ddca53c61c0cccc2fefd99a5c4f729a4ec8ad7197f3c64be50abbf663add3e7f` |
| `…/setup_inversion.py` | `cbdf4f0de2b35b6ec87412917fc8f575535113bf90b0084484ffb53c0da78dc8` |
| `theory/finding_neptune_deterministic.md` | `f84a7d385f95a5a61de8df69369b7ec2afb3c05444722211c7977664b8cebdfd` |
| `theory/finding_neptune_bayesian.md` | `ccbe6ab781136533f10ea4a61bf427ed23b1f79f17b2f6ee49afc654829728dd` |

The 9 MB results file `neptune_long_run_results.npz` was not downloaded.

### What was read

Both theory pages in full. `meta.yml` in full. The deterministic notebook: all prose cells to
section 4.2 and the printed outputs of 4.1 and 4.2. The Bayesian notebook: the prose cells and
the set-up cell of section 4.2; not its results summary or figures. The Python modules: only
the lines defining the data query and the step counts.

### Claims

1. Unknowns: seven, Neptune's mass and its position and velocity vectors at the start, "in the
   year 1775". [read: both theory pages; deterministic notebook, cell 4]
2. Data: the predicted and observed quantities are Uranus's Cartesian coordinates x, y, z, not
   longitudes. The "real data" are JPL Horizons vectors of Uranus relative to the solar-system
   barycentre, ecliptic J2000, from 1775-01-01 at a step of one year; T = 190 years, so N = 3T
   = 570 numbers. [read: theory pages; Bayesian notebook § 4; `neptune_deterministic_methods.py`
   line 754; deterministic notebook § 4]
3. Forward model: Newtonian N-body with the Sun and planets, integrated by a fourth-order
   Runge–Kutta step of one day compiled with numba. [read: theory pages; notebooks, cells 9–10]
4. Deterministic inversion: `scipy.optimize.least_squares` (`lm`, `trf`) through CoFI, with a
   regularisation weight chosen by Morozov's discrepancy principle (χ²/N nearest 1); on real
   data the per-coordinate σ is re-estimated from the residuals until it settles. [read:
   `meta.yml`; deterministic notebook § 4]
5. Results printed. Synthetic (noise 10⁻³ AU in x and y, 10⁻⁵ in z): mass recovered to 2.57 %,
   position components to 1–13 %. Real data: mass 5.129 × 10⁻⁵ against 5.151 × 10⁻⁵ solar
   masses (0.42 %), x and y to about 1 %, z off by 83 %, v_z by 42 %; the search starts from a
   "known-good seed model" close to the truth. [read: deterministic notebook, outputs of 4.1
   and 4.2]
6. Residual level on real data: effective σ of 7.8 × 10⁻² AU in x, 7.4 × 10⁻² AU in y and
   9.8 × 10⁻⁴ AU in z, described as "observational noise plus model mismatch for this forward
   model". No residual in arcseconds is given. [read: deterministic notebook 4.2 output;
   Bayesian notebook 4.2]
7. Bayesian inversion: emcee's stretch move through CoFI, 28 walkers on real data, a diagonal
   Gaussian likelihood, a uniform prior on synthetic data and on real data a Gaussian prior
   centred on the deterministic solution; the notebook runs a short validation chain by default
   and points to an 80 000-step chain for the real data. [read: `meta.yml`; Bayesian notebook
   §§ 2.1, 3.2, 4.2]
8. The notebook points to a thesis for the problem, at `www.diva-portal.org` (the KTH thesis
   named in this repository's README, by the address; I did not open it). [read: deterministic
   notebook, cell 4]
9. I found no evidence computation, no multimodal treatment and no historical observations in
   the files read. [read: the cells listed above; the Bayesian results cells were not read]

## What it means for this repository

Both use modern ephemeris positions as the data; neither fits the historical table. Within
that limit discoverneptune does several things this repository also does, and did them first
in time: a symplectic integrator (REBOUND's, not its own), a sampled posterior, the
mass–distance ridge shown explicitly, a prior-sensitivity comparison, a held-out check (on
synthetic data), calibration checks, and a Mercury counter-example. CoFI does a point fit and
a posterior on Cartesian positions with RK4. The README's four "does differently" lines should
be read against section A: "full posteriors" and "the mass–distance degeneracy made explicit"
are not differences from discoverneptune as it stands; the data, the evidence over a
partitioned prior and the error model are. Details are in `novelty.md`.

## Needs a browser / library

- The discoverneptune source repository and its licence: `https://github.com/jbwhit/discoverneptune`
  answered 404 to an unauthenticated request. Someone signed in to GitHub should check whether
  it is private or renamed, and read its `LICENSE` and the bundled data's provenance.
- The figures of both (they are images; I read the captions and printed numbers only).
