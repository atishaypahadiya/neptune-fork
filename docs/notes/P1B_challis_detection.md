# Phase 1b, step 4: Challis's search as a detection model

Pre-registered plan. Written 2026-11-11 after step 2 (`P1B_challis_sweeps.md`) and step 3c;
nothing below § 1–8 runs until the plan is committed, and the pre-registered sections are then
byte-identical (amendments and results only below the lines). Code and tests cite this note.

## Pre-registration (before any code or run)

### 0 Where step 2 left it

The printed record gives twelve rows of search log (`challis_sweeps.csv`). Challis recorded
Neptune three times: 4 August (sector sweep for reference stars), 12 August (No. 49 of a
fixed-telescope zone, "a star of the 8th magnitude") and 29 September ("it seems to have a
disk"). He did not recognise it before the news from Berlin reached him on 1 October. On
30 July the planet was outside his zone (primary); on 29 July it is not determinable. The
record prints no right-ascension limit and no declination setting for any night, and no dates
for the sessions between 12 August and 28 September. Neptune at about magnitude 8 was three
magnitudes brighter than his stated limit of the 11th, so depth is not what failed.

### 1 Question

With what probability would Challis have **recognised** Neptune, under his actual procedure and
under four named counterfactuals, given his actual nights, instrument and limit?

Two tiers of recognition, both reported for each case:

- **T1, flagged:** on some date he marks one object as a candidate planet, by a noted disc, by
  a star present in one series and absent from another, or by a star absent from a chart.
- **T2, confirmed:** T1 followed by a second observation showing motion, on the first usable
  night after T1.

**Primary outcome:** P(recognised before 23 September 1846, Galle's night), tier T1.
**Secondary:** by 1 October, the end of his search. T2 is reported beside both.

### 2 Frozen inputs

| input | file or source | mark |
|---|---|---|
| search log, twelve rows | `data/historical/adams_challis/challis_sweeps.csv`, sha256 f9061de02545… | P |
| supporting statements | `challis_sweeps_statements.csv` | P |
| Neptune's place per date, mean and apparent, equinox of date and of 1800 | `neptune_track_1846.csv` (7f688e20b1f7…), `neptune_on_challis_dates.csv` (815001ff41d6…) | D |
| Adams's 325° guidance column (the one Challis used) | as tabulated in `neptune_on_challis_dates.csv` | P |
| Adams's S4 and S5 solutions | `adams_solutions.csv` | P |
| Le Verrier's places of June and August 1846 | his own papers as cited in step 1; Kollerstrom's Appendix II only as a cross-check | P, to be transcribed in the pre-run |

Computed in the pre-run and frozen before any probability is computed:

- **Neptune's magnitude** on each date from its heliocentric and geocentric distances by the
  standard relation V = V(1,0) + 5 log₁₀(r Δ), the constant and its source recorded. Expected
  near 7.8; Challis's own estimate was the 8th.
- **Neptune's apparent diameter** on each date from the equatorial radius and Δ. Expected
  near 2.3″.
- **Neptune's motion** between each pair of his nights, in seconds of right ascension and
  minutes of declination, against his stated reading errors (≤ 4″ in declination with the comb,
  ≤ 3″ with the sector micrometer).
- **The geocentric places** implied on each date by S4, by S5 and by each of Le Verrier's two
  places, by the same conversion step 2 used for the 325° column.

### 3 Named inputs with fallbacks

Each named input is used if obtained before the run; otherwise the fallback is used and every
number that depends on it is flagged in the results.

| name | wanted | fallback, flagged [S] or [prior] |
|---|---|---|
| **N1 Berlin sheet limits** | Hora XXII's actual limits in right ascension and declination, its epoch, and its faintest plotted magnitude, from the sheet or the 1859 catalogue | Krajnović 2021: one hour of right ascension with a 4-minute overlap each side (21h 56m to 23h 04m), ±15°, equinox 1800. Faintest magnitude: the secondary value recorded in the pre-run with its page |
| **N2 Cambridge observing books** | per night: date, method, zone limits in right ascension and declination, star count, the cloud lines, the comparisons made (UA/Obsy, Northumberland series; UA/Obsy D.3) | the printed record of § 0, with the unknowns as the declared priors of § 4 |
| **N3 Airy's Suggestions of 13 July** | the declination intervals and their extents in right ascension (RGO 6/96a) | Challis's summary: 30° of longitude by 10° of latitude centred on 325°, zones in 1° declination intervals |
| **N4 Adams's July ephemeris, all five columns** | RGO 6/96a | the 325° column only, as reconstructed in step 2 |

N1 has a second fallback run, reported beside the first: the limits applied to right ascension
of date. Step 2 showed this is a frame mismatch if the sheets are in the equinox of 1800, but
it is what one secondary source assumes, and the counterfactual C2 turns on it.

### 4 The model

**Abduction: the draws are conditioned on the record before any counterfactual.** Each draw of
the priors below is run through his actual procedure A0 and kept only if it reproduces what
the printed record fixes:

1. the planet recorded on 4 August and on 12 August, and not on 30 July (geometry);
2. no recognition by 28 September 1846: no disc noted on either August recording, and nothing
   flagged by the comparison he ran, which stopped at No. 39;
3. a disc noted on 29 September, when he was looking for one;
4. the total number of stars recorded under the drawn sessions, durations, extents and rates
   falls in 3000 to 4000 ("above 3000 positions", printed).

Rejection sampling. Condition 1 is imposed through the zone priors, which are restricted to
settings consistent with it, so that the acceptance fraction measures conditions 2 to 4 and
is not swamped by the narrowness of a 9′ zone; the fraction of an unrestricted zone prior that
satisfies condition 1 is reported beside it. **His record enters through the abduction alone:**
no prior below is set from it, so the acceptance fraction is informative about whether the
independent priors are compatible with what happened, and it is reported overall and per
condition (2 and 3 together, 4 alone). Reported with it: for each behavioural parameter, and
for the session-count, duration, rate and extent priors, the posterior against the prior
(mean and 5 % to 95 % range, with a figure). **If the acceptance fraction is below 1 %, the model is reported as inconsistent
with history and the run stops before the counterfactuals.** C1 to C4 run on the accepted
draws only; a counterfactual changes the procedure and keeps the draw's parameters. Draws
continue until 100,000 are accepted.

**Per session, three events.** In zone (Z); recorded given in zone (R); recognised given
recorded (K). K has three routes, and a case states which are open:

- **K_disc:** the disc is noted at the eyepiece. One free probability per recording, with two
  values: p_disc_routine when he is recording positions only, p_disc_looking when he is
  examining the brighter stars for a disc, as on 29 September.
- **K_series:** two series covering the same ground on different nights are compared star by
  star. The planet is flagged if it is in one and the other series covered its place. Given
  the motion of § 2 against his reading errors, a compared planet cannot be matched to itself;
  the free probability is p_follow, that he treats an unmatched star as a candidate and not as
  an omission.
- **K_chart:** recorded stars are checked against a chart. The planet is flagged if the chart
  covers its place and plots stars to its magnitude; free probability p_follow again.

**Declared priors for what the record does not print** (used only where N2 is not obtained).
Each is sampled; seed `[2026, 11, 4, 1]` for the draws, the same accepted draws through every
case, and `[2026, 11, 4, 2, case]` for the events within a case.

| unknown | prior |
|---|---|
| declination of the 30 July and 12 August zone (9′ broad, the same on both nights) | uniform over the settings that contain the planet on 12 August and exclude it on 30 July by the declination edge; a second branch, weight 0.5, in which the planet entered by the right-ascension end and the declination is uniform over settings containing both places |
| right-ascension extent of each zone | hours of the session × 15° per hour for a fixed telescope (the sky's own rate); for a moving telescope the same length scaled by a factor uniform on 0.5 to 2 |
| the 29 July zone | same declination setting as 30 July (primary); extent as above |
| sessions between 12 August and 28 September | number uniform on 8 to 25; dates drawn without replacement from the nights not excluded by the Moon rule below; ground "earlier in right ascension" than the planet (primary), so Z is false on all of them in the actual procedure |
| usable nights | a night is excluded when the Moon is above the horizon and more than half lit within 30° of the zone during the session; weather: each remaining night usable with probability drawn uniform on 0.3 to 0.6 |
| p_disc_routine | uniform on 0.0 to 0.5. Basis, independent of Challis: the base rate of documented pre-discovery recordings of Neptune in which no disc was noted (Galileo 1612–13, Kowal & Drake 1980, Nature 287, 311; Lalande 1795; J. Herschel 1830; Lamont 1845–46). The count of such recordings, the instruments and the sources are tabulated in the pre-run; the interval is left wide because those observers used other apertures and powers |
| p_disc_looking | uniform on 0.1 to 1.0. Basis, independent of Challis: a disc of 2.3″ at power 166 on the 11¾-inch Northumberland refractor, against the resolution of that aperture and the seeing of an ordinary night; the arithmetic is recorded in the pre-run |
| stars recorded per hour of sweeping | uniform on 40 to 90 |
| hours per session | uniform on 1 to 3.5 |
| p_follow | uniform on 0.5 to 1.0 |
| R, recorded given in zone, at magnitude 8 against a limit of 11 | uniform on 0.9 to 1.0; clouds enter through the cloud-line rule of the case |
| desk hours to compare one zone with another, star by star | uniform on 1 to 6 hours per pair of zones (his own partial comparison of 39 stars was done "probably the next day") |
| desk hours he gives to comparison per day | uniform on 0.5 to 3 |
| delay before a comparison starts | 1 day after the later of the two sessions (C1's "promptly"); a second run at 7 days is reported beside it |

These priors are judgement, declared here so that they cannot be tuned after the results. The
record informs only some of them: conditions 2 and 3 update the two disc probabilities and
condition 4 the session count, duration and rate; none carries information on p_follow or on
the desk hours, whose posteriors equal their priors, and the results say so. Per accepted
draw, the probability of an outcome is the probability given that draw's inputs, computed in
closed form where the case allows and otherwise by 200 inner event draws with the same random
numbers across cases. Each P is given at the posterior mean with its 5 % to 95 % range over
the posterior of the judgement inputs.

### 5 The cases

- **A0, his actual procedure.** Pointing, nights and methods as in the log. K_disc open on
  every recording (routine on 4 and 12 August, looking on 29 September). K_series as he ran
  it: 12 August against 30 July, stopped at No. 39. K_chart closed. T2 after 29 September
  needs a usable night before 1 October; 30 September was lost to the Moon (primary).
- **C1, prompt comparison.** A0 with the comparisons carried through. The comparison has a
  **completion date**, modelled and not assumed: it starts after the declared delay, proceeds
  at the declared desk hours per day against the declared hours per pair of zones, and reaches
  the planet's entry at a point in the series proportional to its number (No. 49 on
  12 August). T1 is dated to the day the comparison reaches the planet; T2 to the first usable
  night after. The distribution of the T1 date is reported. Two parts, reported apart:
  - **C1a:** the 12 August series compared with 30 July through No. 49 and beyond. The planet
    is in the August series and the July zone covered that ground (primary: the August zone
    included the whole of the July one), so it is flagged with probability p_follow, unless
    its place lay in the part of the August zone beyond the July zone's extent; that fraction
    comes from the extent prior.
  - **C1b:** the 4 August sector stars compared with the 12 August zone. Flagged if the
    12 August zone covered the 4 August place, which depends on the declination prior (the
    two places are 4.5′ apart in declination and 49 s in right ascension).
- **C2, the Berlin Hora XXII chart.** A0 with every recorded star checked against the sheet
  Challis had, each session's check dated by the same desk-hours model as C1. K_chart open
  under N1. Run under N1's fallback,
  under its second fallback, and under N1 itself if obtained.
- **C3, starting at Adams's S4/S5 place.** The same nights, methods, zone breadth and rate,
  with the zones laid out from the geocentric place of S4 (and separately S5) for each date in
  place of the 325° column. The order of zones follows N3's fallback: declination intervals
  outward from the centre, each swept twice before the next. Z is computed, not conditioned.
  Recognition routes as in A0, and again with C1's prompt comparison.
- **C4, starting with Le Verrier's place.** Two parts:
  - **C4a:** as C3, centred on his June place, from 29 July.
  - **C4b:** the procedure of 29 September (a breadth in declination within his limits of
    longitude, brighter stars examined for a disc) adopted from 29 July with his June place.
    The limits of longitude are his printed ones, transcribed in the pre-run.

For C3 and C4 the zone layout needs a rule for where the first zone sits; the rule is: the
first zone is centred in declination on the tabulated place for that night, and its
right-ascension extent is centred on the tabulated right ascension.

### 6 Checks before reading

1. **Abduction (§ 4).** The acceptance fraction is at least 1 %, or the run stops. On the
   accepted draws A0 reproduces the record by construction; this is asserted in a test.
2. **A0 after 29 September.** On the accepted draws, A0's probability of T2 by 1 October is
   reported against the fact that he did not confirm before the news arrived.
3. **Geometry.** The computed planet-in-zone verdicts for A0 reproduce step 2's rows (inside
   on 4 and 12 August and 29 September, outside on 30 July).
4. **Tornado.** Per case, each judgement input and each named-input fallback is set in turn
   to its posterior 5 % and 95 % values with the others left to vary; the inputs are ranked by
   the change in the primary probability between the two settings, and the ranking is given
   as a table and a tornado figure.
5. **Prior width.** Each P is recomputed with every prior interval halved about its centre;
   a reading that changes row under this is marked "prior-dependent".

### 7 Reading, fixed in advance

Per case, on the primary outcome. **The 5 % to 95 % range is over the posterior of the
judgement inputs** (the accepted draws of § 4); it is not a sampling error and does not shrink
with more draws. Every reading names the case's top-ranked input from the tornado of § 6.4:

| condition | reading |
|---|---|
| lower end ≥ 0.8 | "would very probably have recognised the planet before Galle's night; most sensitive to [top input]" |
| upper end ≤ 0.2 | "would very probably not have; most sensitive to [top input]" |
| range inside 0.2 to 0.8 | "an open chance", with the mean stated and the top input named |
| range crossing 0.2 or 0.8 | "depends on [top input]", no single statement |

A case that rests on a fallback of N1 to N4 carries "conditional on [fallback]" in its reading.
C2 is read once per version of N1, and if the versions fall in different rows the reading is
"depends on the Berlin sheet's limits". T2 and the secondary outcome (by 1 October) are
reported beside every row and not read. No case is ranked against another on its marginal
range unless their ranges do not overlap.

**Paired differences.** For each counterfactual Cᵢ, the per-draw difference
P(Cᵢ) − P(A0) on the primary outcome, over the same accepted draws, with its mean and 5 % to
95 % range. "Cᵢ improves on A0" is stated only when that range excludes 0; "Cᵢ does worse
than A0" only when it lies wholly below 0; otherwise "no difference established". The
differences between counterfactuals are tabulated the same way and not read.

Not claimed in any outcome: what Challis ought to have done, or a probability for the British
search as a whole. The model covers his nights and instrument only.

### 8 Compute, outputs, tests

- Seconds of compute; one process; no detached run. A timed smoke run of 1,000 draws sets the
  estimate before the full run, and the estimate is logged.
- Code `src/fork/historical/challis_1846.py` and runner `src/fork/export/p1b_challis.py`
  (`prerun`, `draws`, `summarize`). Tests `tests/test_p1b_challis.py`: the geometry check of
  § 6.3 on the frozen files, the seeds, the reading rows on constructed inputs, the
  abduction (accepted draws reproduce the record; the stop below 1 %), the completion-date
  arithmetic of C1.
- Results below the lines of this note; `docs/notes/P1B_challis_detection.{json,png}`; one
  ledger row `p1b_challis/draws`, local only; the site untouched; nothing in the notes of
  steps 2 to 3c edited.
- Figures: per case, the T1 probability by date from 29 July to 1 October with its range, one
  panel per case, 23 September marked; the prior and posterior of each behavioural parameter;
  the tornado per case.
