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

---

## Dated lines (nothing above the lines edited)

- 2026-11-12: the Northumberland refractor's aperture is "a little more than 11½ inches", focal
  length 19 ft 5¾ in (Airy 1844, *Account of the Northumberland Equatoreal and Dome*), replacing
  the 11¾ inches of § 4's table, which was found in no source read. The prior on p_disc_looking
  is unchanged.
- 2026-11-12: the start places of C3 and C4 (S4, S5, Le Verrier's two places) lie at ecliptic
  latitude 0, as both predictors assumed orbits in the ecliptic. Neptune's latitude of about
  −0.5° is the source of most of the declination offsets in the table of P2.
- 2026-11-12: the zone order of C3 and C4 is alternating north and south outward from the
  centre, each zone swept twice (telescope fixed, then moving) before the next. This is the
  interpretation of N3's fallback.

## Pre-run (2026-11-12; frozen before any full draw)

Frozen files:

| file | sha256 |
|---|---|
| `docs/notes/P1B_challis_prerun.json` | `5d2dfad0842d164712014a526c1ee8b191cd2d59e06d21e963f9d06f89a1aba5` |
| `docs/notes/P1B_challis_smoke.json` | `c1b969c2f21ab2ddf9a239a5868dc3bd42f1c32dcc8180e652d1eb90e93c6bfb` |
| `data/historical/adams_challis/prediscovery_recordings.csv` | `9fb014c020fc5eeb0dc312866a60a4d748c6f2301cd7b8cd23cb1988c5ee9f47` |

`uv run python -m fork.export.p1b_challis prerun` rewrites the first; `smoke` the second, whose
timings change from run to run. The three frozen inputs of § 2 were checked against their
hashes before the pre-run.

### P1. Neptune per night

Constants: V = 5 log₁₀(r Δ) − 6.89 (Mallama & Hilton 2018, Astronomy and Computing 25, 10,
Eq. 16, the branch for dates before 1980.0; the fit rests on photometry of 1954 to 2014 and is
carried back to 1846, which is an extrapolation). Equatorial radius 24,764 km (IAU WGCCRE 2015:
Archinal et al. 2018, Celest. Mech. Dyn. Astron. 130, 22). 1 au = 149,597,870.7 km. Places are
geometric, mean equinox of date, at 23:30 UT, from two-body arcs through the cached Horizons
states of 1846-09-23 23:00 TDB; against step 2's track the right ascension differs by at most
6.4 s, which is the difference of session times and of the two-body arc. All 65 nights are in the JSON.

| night | RA | Dec | V | diameter | r (au) | Δ (au) |
|---|---|---|---|---|---|---|
| 1846-07-29 | 21h 58m 49.1s | −12° 54.4′ | 7.81 | 2.35″ | 30.013 | 29.063 |
| 1846-07-30 | 21h 58m 43.3s | −12° 54.9′ | 7.81 | 2.35″ | 30.013 | 29.057 |
| 1846-08-04 | 21h 58m 13.7s | −12° 57.7′ | 7.81 | 2.35″ | 30.013 | 29.032 |
| 1846-08-12 | 21h 57m 24.4s | −13° 02.2′ | 7.81 | 2.35″ | 30.013 | 29.007 |
| 1846-09-22 | 21h 53m 20.5s | −13° 24.0′ | 7.82 | 2.34″ | 30.011 | 29.172 |
| 1846-09-29 | 21h 52m 47.5s | −13° 26.9′ | 7.83 | 2.34″ | 30.011 | 29.245 |
| 1846-10-01 | 21h 52m 38.9s | −13° 27.6′ | 7.83 | 2.33″ | 30.011 | 29.268 |

Motion between his nights:

| from | to | ΔRA | ΔDec |
|---|---|---|---|
| 1846-07-29 | 1846-07-30 | -5.8 s | -0.5′ |
| 1846-07-30 | 1846-08-04 | -29.6 s | -2.7′ |
| 1846-08-04 | 1846-08-12 | -49.2 s | -4.5′ |
| 1846-08-12 | 1846-09-29 | -276.9 s | -24.7′ |

At power 166 a disc of 2.35″ subtends 6.5′ at the eye.

### P2. The start places

Transcribed from the papers. Le Verrier, 1 June 1846 (Comptes rendus 22, 917): heliocentric
longitude 325° at 1 January 1847, with "not an error of 10 degrees", read as 315° to 335°.
Le Verrier, 31 August 1846 (Comptes rendus 23, 428–438): semi-major axis 36,154; period
217,387 years; eccentricity 0,10761; longitude of perihelion 284° 45′; mean longitude at
1 January 1847 318° 47′; mass 1/9300; true heliocentric longitude 326° 32′ and distance 33,06
at 1 January 1847 (p. 432); search "entre 321 et 335 degrés de longitude héliocentrique"
starting from 326° 32′ (p. 436); mass between 1/14500 and 1/4700; predicted disc 3″,3
(p. 433). He prints no geocentric place. The elements and the limits were read on the page
images; Kollerstrom was not needed as a cross-check for any figure used.

Geocentric places, equinox of date, latitude 0:

| place | 29 July: RA, Dec | geocentric longitude | 22 September: RA, Dec | geocentric longitude |
|---|---|---|---|---|
| Adams S4 | 21h 46m 29.7s, −13° 25.7′ | 324.32° | 21h 41m 32.8s, −13° 50.9′ | 323.04° |
| Adams S5 | 22h 04m 49.3s, −11° 48.4′ | 329.07° | 21h 59m 47.4s, −12° 15.7′ | 327.76° |
| Le Verrier, 1 June | 21h 48m 29.3s, −13° 15.4′ | 324.83° | 21h 44m 07.2s, −13° 37.8′ | 323.70° |
| Le Verrier, 31 August | 21h 54m 05.1s, −12° 46.1′ | 326.28° | 21h 49m 06.8s, −13° 12.1′ | 324.99° |

Neptune minus the start place (ΔRA in minutes of time, ΔDec in arcminutes):

| night | 325° column | S4 | S5 | Le Verrier June | Le Verrier August |
|---|---|---|---|---|---|
| 1846-07-29 | +8.3 min, +9.6′ | +12.3 min, +31.2′ | -6.0 min, -66.0′ | +10.3 min, +21.0′ | +4.7 min, -8.3′ |
| 1846-08-12 | +8.1 min, +8.2′ | +12.2 min, +30.2′ | -6.1 min, -66.8′ | +10.1 min, +19.1′ | +4.6 min, -9.3′ |
| 1846-09-22 | +7.1 min, +3.4′ | +11.8 min, +26.9′ | -6.4 min, -68.3′ | +9.2 min, +13.8′ | +4.2 min, -11.9′ |
| 1846-09-29 | +7.1 min, +3.0′ | +11.7 min, +26.4′ | -6.5 min, -68.4′ | +9.1 min, +13.1′ | +4.2 min, -12.2′ |

Le Verrier's August place enters the pre-run as asked; the plan's C4 uses the June place.

### P3. Pre-discovery recordings (the basis of the prior on p_disc_routine)

`data/historical/adams_challis/prediscovery_recordings.csv`, reviewed before freezing.

| observer | nights | instrument | disc noted | identification | mark |
|---|---|---|---|---|---|
| Galileo | 28 Dec 1612; 27 and 28 Jan 1613 | Galilean refractor, aperture and power not verified | no | Kowal & Drake 1980, Nature 287, 311 | abstract read; wording from Standish and from Bernieri et al. 2022, secondary |
| M. Lefrançais de Lalande | 8 and 10 May 1795 | mural quadrant, 2½ pouces (about 68 mm), power not found | no | Mauvais 1847, Comptes rendus 24, 666–669 | primary, read |
| Lamont | 25 Oct 1845; 7 Sep 1846 | Reichenbach meridian circle, aperture and power not verified | no | Hind 1850, Astron. Nachr. 30, 257 | primary identification; instrument secondary |
| Lamont | 11 Sep 1846 | same | no | Schönfeld 1886, Astron. Nachr. 114, 59 | primary identification; instrument secondary |
| J. Herschel | 14 July 1830 | not verified | not established as a recording | Kollerstrom 2006, JAHH 9, 151: "nearly found it" | unverified; excluded from the count |

Count: 8 recordings on 8 nights by 3 observers, none with a disc noted. The instruments could
mostly not show a 2.35″ disc (a power of 20, for example, gives about 47″ at the eye), so the
base rate supports the width of the 0 to 0.5 prior and does not narrow it. Description only:
on 28 January 1613 Galileo noted that the star and its neighbour had seemed farther apart the
night before, an anomaly noted and not followed up; it bears on p_follow, whose frozen prior is
unchanged. Sources differ on Galileo's January nights (27 and 28, or 26/27 and 27/28).
Power 166 and the field of 9′ are Challis's own figures (MNRAS 7, 145).

### P4. The Moon rule per night

Lit fraction, separation from Neptune's place and altitude at Cambridge at 21:30, 23:30 and
01:30 UT are in the JSON. The rule of § 4 excludes: 08-06, 08-07, 08-08, 08-09, 09-02, 09-03, 09-04, 09-05, 09-29, 09-30, 10-01. Against the record: 4 August
is not excluded (lit 0.92, 48° away) and he observed; 29 September is excluded by the rule (lit
0.69, 28° away) though he observed for 3.5 hours; 30 September is excluded and he
records the Moon in the way. The nights of the record are therefore set as the record gives
them, and the rule acts only on the nights the record leaves open.

### P5. Timed smoke run

Seeds as § 5. 1000 accepted draws of 10061: acceptance 0.099 (conditions 2 and 3:
0.322; condition 4: 0.307), above the 1 % stop. Draws and abduction 0.3 s;
all 14 case variants in closed form 1.6 s; 1,000 inner event draws per draw and case
1.8 s. Estimate for 100,000 accepted draws: 3 min in closed form, 6 min with the
inner draws, on one core.

Implementation choice: every case has a closed form for the per-draw probability, so the
closed form is the reported value and 1,000 inner event draws per draw (the plan: 200) serve
as its check. The Monte Carlo share is the mean of p(1 − p)/n over draws against that plus the
variance of p between draws.

| case | primary, mean | 5–95 % | by 1 Oct | T2 primary | worst inner − exact (s.e.) | MC share, n = 200 | n = 1,000 |
|---|---|---|---|---|---|---|---|
| A0 | 0.325 | 0.026 to 0.690 | 0.786 | 0.325 | 3.0 | 0.020 | 0.004 |
| C1a | 0.831 | 0.616 to 0.987 | 0.947 | 0.831 | 3.3 | 0.046 | 0.010 |
| C1b | 0.513 | 0.036 to 0.955 | 0.847 | 0.513 | 3.2 | 0.009 | 0.002 |
| C1a_7d | 0.831 | 0.616 to 0.987 | 0.947 | 0.831 | 3.3 | 0.046 | 0.010 |
| C1b_7d | 0.513 | 0.036 to 0.955 | 0.847 | 0.513 | 2.9 | 0.009 | 0.002 |
| C2_1800 | 0.325 | 0.026 to 0.690 | 0.786 | 0.325 | 2.7 | 0.020 | 0.004 |
| C2_date | 0.944 | 0.829 to 1.000 | 0.983 | 0.944 | 3.1 | 0.076 | 0.016 |
| C3_S4 | 0.397 | 0.036 to 0.798 | 0.812 | 0.388 | 3.4 | 0.015 | 0.003 |
| C3_S4_prompt | 0.788 | 0.225 to 0.987 | 0.950 | 0.764 | 3.2 | 0.012 | 0.003 |
| C3_S5 | 0.000 | 0.000 to 0.000 | 0.683 | 0.000 | 0.0 | 0.000 | 0.000 |
| C3_S5_prompt | 0.000 | 0.000 to 0.000 | 0.683 | 0.000 | 0.0 | 0.000 | 0.000 |
| C4a | 0.401 | 0.032 to 0.805 | 0.810 | 0.401 | 3.3 | 0.015 | 0.003 |
| C4a_prompt | 0.846 | 0.625 to 0.988 | 0.952 | 0.844 | 2.8 | 0.042 | 0.009 |
| C4b | 0.353 | 0.018 to 0.900 | 0.781 | 0.353 | 3.0 | 0.009 | 0.002 |

These are smoke numbers on 1,000 draws, not results, and are not read against § 7.

Two things they show. C2 under the first fallback of N1 equals A0: the planet's right ascension
for the equinox of 1800 is west of the sheet's limit on 4 and 12 August, as step 2 found.
C3 from S5 gives 0 before 23 September by timing, not coverage: Neptune lies 66′ from S5's
place in declination, the 7th step on one side, the 14th or 15th zone, so the 27th session at
the earliest; by late September the offset is 68.4′, in the 8th step, and no draw has more than
28 sessions. No draw of the 1,000 reaches the band; had one done so, the extent in right
ascension (30 to 105 min each side against the 6 min needed) would cover the planet. From S4
the band is reached in 995 draws (959 before 23 September; the 8th to 13th session, median
8 September), and the extent in right ascension covers the planet in every one.

### P6. Implementation choices the plan does not fix (open for ruling before the full draws)

1. Replay reading of the abduction: an accepted draw keeps its parameters and the events are
   drawn afresh, so P(A0, primary | draw) = 1 − (1 − p_disc_routine)², not 0.
2. A0's zones: the plan places zones in right ascension only for C3 and C4. In A0, C1 and C2
   the record fixes the recordings of 4 and 12 August and 29 September, and the sessions
   between 12 August and 28 September are taken not to cover the planet (the record has no
   further recording).
3. C4b: "considerable breadth" is not quantified. The breadth in declination is uniform on
   30′ to 120′, the length in right ascension scales as 9′ over the breadth, one sweep per
   zone, p_disc_looking at each recording.
4. If fewer usable nights than drawn sessions fall between 13 August and 28 September, the
   sessions are cut to the usable nights.
5. Hours per session are drawn per session; stars per hour once per draw; the planet's place
   in the 12 August series is 49 over that session's stars, capped at 1.
6. Desk work: a comparison flagged on day t is confirmed on the first usable night from t; a
   disc noted at the eyepiece is confirmed on the next usable night.
