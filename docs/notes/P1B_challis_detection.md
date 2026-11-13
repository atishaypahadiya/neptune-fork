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

## Dated lines, 2026-11-12: rulings on P6

The rulings were decided, and the code implementing them was written, before the full draws
were launched (00:27 PST). These lines were written into the note about two minutes after the
launch, while the draws ran, and were committed alone and pushed before any result of the
full run was read (no chunk opened, no summary made). The smoke record's split of line 7 is in
the working tree and is committed with the results.

1. **Abduction reading.** Primary, and the one read against § 7: the counterfactual with the
   realised events held. A recording that a case shares with the actual history (the same
   night, the same sweep, the same mode) keeps its realised outcome: recorded on 4 and
   12 August with no disc noted, and recorded on 29 September with the disc noted. Events
   the history does not contain are drawn afresh. Under this reading A0 is 0 before
   23 September by construction. The paired differences of § 7 are taken under this reading.
   Secondary, reported and not read: the replay with every event redrawn, labelled "his
   procedure as a plan", A0 included. Both readings use the same accepted draws and seeds.
   What counts as shared, as implemented: in A0, C1 and C2 the three recordings of the record;
   in C3 and C4 the 29 September sweep only, since their zones and their 4 August sector
   sweep are laid from another place and are other recordings.
2. **The sessions of 13 August to 28 September do not cover the planet** in A0, C1 and C2.
   Basis: step 2's primary statement, row `1846-08-13-to-09-28` of `challis_sweeps.csv`, the
   ground "earlier in right ascension" than the planet (Challis's letter: "too early in right
   ascension"), and "did not again fall in with the planet".
3. **C4b, declared choice.** The breadth in declination of each sweep is uniform on 30′ to
   120′; the length in right ascension scales as 9′ over the breadth; one sweep per zone;
   p_disc_looking at each recording.
4. **Too few usable nights, declared choice.** If fewer usable nights than drawn sessions
   fall between 13 August and 28 September, the sessions are cut to the usable nights.
5. **Hours and rate, declared choice.** Hours are drawn per session; stars per hour once per
   draw; the planet's place in the 12 August series is 49 over that session's stars, capped
   at 1.
6. **Confirmation, declared choice.** A comparison that flags the planet on day t is
   confirmed on the first usable night from t; a disc noted at the eyepiece is confirmed on
   the next usable night.
7. **The smoke record is split.** `docs/notes/P1B_challis_smoke.json` now holds the content
   without timings, sha256
   `6c0c25969c707b6696a7852d7ef455184cc744f31297942dbab80f437a148d8a`; the timings are in
   `docs/notes/P1B_challis_smoke_timings.json`, not hashed. The numbers are those of the
   pre-run's code (commit d754436), under the replay reading; a rerun of `smoke` writes to
   the cache and leaves the record alone.

Fixed with these lines, before the full draws, where the pre-run's code fell short of § 5:
C1a now applies the plan's exception (the planet's 12 August place may lie in the part of the
August zone beyond the July zone's extent). In the branch where the planet entered by the
declination edge, the July zone lies within the August one at a setting uniform over the
possible ones, and the planet sits at 49 over the session's stars along the August zone. The
smoke record predates this and has C1a without the exception.

---

## Results (2026-11-12; pre-run d754436, rulings e6889d4; draws 00:27–00:47 PST on one core; nothing above the lines edited)

Order of events: the pre-run was frozen in d754436; the rulings on its open choices were
decided and implemented before the launch, written into the note about two minutes after it,
and committed alone and pushed as e6889d4 before any result was read. The results below were
read after that push.

Files: `docs/notes/P1B_challis_detection.json` (sha256 `54d15eafd18ca979729d9c3e7787816a5fc85a5c6440644ad89538909473b947`),
`P1B_challis_detection.png` (probability by date per case), `P1B_challis_detection_posterior.png`
(inputs before and after the abduction), `P1B_challis_detection_tornado.png`.
`uv run python -m fork.export.p1b_challis draws` then `summarize` reproduce them from the
frozen pre-run; the draws resume per chunk and log to `data/cache/p1b_challis/progress.log`.

### Summary

With the events of his own record held, Challis's actual procedure could not have recognised
the planet before 23 September: that is 0 by construction of the reading, and is stated as
such, not found. Four changes of procedure raise the probability with a range that excludes
0: carrying the 12 August comparison with 30 July through (0.71, range 0.50 to 0.97), a start
from Adams's S4 place (0.41, or 0.79 with prompt comparison), a start from Le Verrier's June
place (0.41, or 0.85 with prompt comparison), and the 29 September procedure adopted on
29 July (0.35). None of them reaches the row "very probably": every range crosses 0.2 or 0.8,
so each reads "depends on" its top input, which is p_follow for the comparisons and
p_disc_routine for the changes of start place. A start from S5 gives 0, by timing. The
Berlin sheet he had gives 0 if its limits are for the equinox of 1800 and 0.92 if they are
applied to right ascension of date: the chart case depends on the sheet's limits.

### R1. Abduction and checks (§ 6)

- **Acceptance** 0.0996 (100,000 of 1,003,853), above the 1 % stop. Conditions 2 and 3
  together 0.320; condition 4 alone 0.309. The independent priors are compatible with
  the record: about a third of prior draws reproduce the August silences and the September
  disc, and about a third the count of stars.
- **Condition 1 under an unrestricted zone prior:** a 9′ zone set uniformly in Airy's 10° of
  latitude contains the planet on 12 August in 0.015 of settings, and also excludes it on
  30 July by the declination edge in 0.012. This is arithmetic on the band, not a draw.
- **A0 after 29 September (§ 6.2):** held reading, T2 by 1 October is 0 in every draw. The
  disc was noted on 29 September and 30 September was lost to the Moon, which agrees with
  the fact that he did not confirm before the news arrived. Replay reading: 0.34, all of it
  from discs noted in August.
- **Geometry (§ 6.3):** inside on 4 and 12 August and 29 September, outside on 30 July, in
  every draw; asserted in a test on the frozen pre-run.
- **Closed form against 1,000 inner event draws** per draw, case and reading: worst mean
  difference 0.00007; worst single draw 4.8 standard errors, among 2.8 million comparisons.

Posterior against prior, mean (5 % to 95 %):

| input | prior | posterior |
|---|---|---|
| p disc routine | 0.25 (0.0249 to 0.475) | 0.197 (0.0146 to 0.448) |
| p disc looking | 0.55 (0.145 to 0.955) | 0.674 (0.247 to 0.975) |
| p follow | 0.75 (0.525 to 0.975) | 0.75 (0.525 to 0.975) |
| p recorded | 0.95 (0.905 to 0.995) | 0.95 (0.905 to 0.995) |
| stars per hour | 65 (42.5 to 87.5) | 70.6 (51.1 to 87.9) |
| moving extent factor | 1.25 (0.575 to 1.92) | 1.25 (0.575 to 1.93) |
| weather | 0.45 (0.315 to 0.585) | 0.453 (0.318 to 0.585) |
| desk hours per pair | 3.5 (1.25 to 5.75) | 3.5 (1.25 to 5.75) |
| desk hours per day | 1.75 (0.625 to 2.88) | 1.75 (0.625 to 2.88) |
| n unknown sessions | 14.8 (8 to 22) | 16.2 (11 to 22) |
| hours per session mean | 2.25 (1.97 to 2.53) | 2.27 (2.01 to 2.54) |
| hours total | 45.8 (30.6 to 63.6) | 49.1 (37.7 to 63.7) |
| stars total | 3.05e+03 (1.76e+03 to 4.75e+03) | 3.44e+03 (3.04e+03 to 3.93e+03) |
| extent jul30 deg | 42.2 (14.9 to 81.7) | 42.6 (15 to 82.1) |
| extent aug12 deg | 33.8 (16.9 to 50.6) | 34.1 (17 to 50.7) |

The record moves the two disc probabilities (routine down, looking up), the rate, and the
session count and hours through the count of stars. As the plan stated, p_follow, the desk
hours and the extents keep their priors: the record carries no information on them.

### R2. The reading (§ 7): realised events held, primary outcome before 23 September

Mean and 5 % to 95 % range over the posterior of the judgement inputs.

| case | P before 23 Sep | row of § 7 | top input | T2 before 23 Sep | T1 by 1 Oct | T2 by 1 Oct | range with priors halved |
|---|---|---|---|---|---|---|---|
| A0 actual procedure | 0.00 (0.00 to 0.00) | would very probably not have | none | 0.00 | 1.00 | 0.00 | 0.00 to 0.00 |
| C1a 12 Aug against 30 Jul | 0.71 (0.50 to 0.97) | depends on | p follow | 0.71 | 1.00 | 0.71 | 0.64 to 0.86 |
| C1b 4 Aug against 12 Aug | 0.28 (0.00 to 0.93) | depends on | p follow | 0.28 | 1.00 | 0.28 | 0.00 to 0.84 |
| C1a, 7-day delay | 0.71 (0.50 to 0.97) | depends on | p follow | 0.71 | 1.00 | 0.71 | 0.64 to 0.86 |
| C1b, 7-day delay | 0.28 (0.00 to 0.93) | depends on | p follow | 0.28 | 1.00 | 0.28 | 0.00 to 0.84 |
| C2 chart, equinox 1800 | 0.00 (0.00 to 0.00) | would very probably not have | none | 0.00 | 1.00 | 0.00 | 0.00 to 0.00 |
| C2 chart, RA of date | 0.92 (0.77 to 1.00) | depends on | p follow | 0.92 | 1.00 | 0.92 | 0.87 to 0.98, prior-dependent |
| C3 from S4 | 0.41 (0.04 to 0.81) | depends on | p disc routine | 0.40 | 1.00 | 0.42 | 0.33 to 0.71, prior-dependent |
| C3 from S4, prompt | 0.79 (0.22 to 0.99) | depends on | n unknown sessions | 0.76 | 1.00 | 0.84 | 0.77 to 0.95 |
| C3 from S5 | 0.00 (0.00 to 0.00) | would very probably not have | none | 0.00 | 1.00 | 0.00 | 0.00 to 0.00 |
| C3 from S5, prompt | 0.00 (0.00 to 0.00) | would very probably not have | none | 0.00 | 1.00 | 0.00 | 0.00 to 0.00 |
| C4a Le Verrier June | 0.41 (0.04 to 0.81) | depends on | p disc routine | 0.41 | 1.00 | 0.41 | 0.32 to 0.71, prior-dependent |
| C4a, prompt | 0.85 (0.63 to 0.99) | depends on | p follow | 0.85 | 1.00 | 0.85 | 0.78 to 0.95 |
| C4b looking from 29 Jul | 0.35 (0.02 to 0.90) | depends on | c4b breadth arcmin | 0.35 | 1.00 | 0.35 | 0.13 to 0.72 |

Readings, as the rule gives them:

- **A0 actual procedure:** would very probably not have; 0 in every draw, no input moves it; conditional on the printed record and the priors of § 4 (N2's fallback).
- **C1a 12 Aug against 30 Jul:** depends on p_follow; conditional on the printed record and the priors of § 4 (N2's fallback).
- **C1b 4 Aug against 12 Aug:** depends on p_follow; conditional on the printed record and the priors of § 4 (N2's fallback).
- **C1a, 7-day delay:** depends on p_follow; conditional on the printed record and the priors of § 4 (N2's fallback).
- **C1b, 7-day delay:** depends on p_follow; conditional on the printed record and the priors of § 4 (N2's fallback).
- **C2 chart, equinox 1800:** would very probably not have; 0 in every draw, no input moves it; conditional on N1's first fallback (the sheet's limits for the equinox of 1800).
- **C2 chart, RA of date:** depends on p_follow; conditional on N1's second fallback (the sheet's limits applied to right ascension of date).
- **C3 from S4:** depends on p_disc_routine; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C3 from S4, prompt:** depends on n_unknown_sessions; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C3 from S5:** would very probably not have; 0 in every draw, no input moves it; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C3 from S5, prompt:** would very probably not have; 0 in every draw, no input moves it; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C4a Le Verrier June:** depends on p_disc_routine; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C4a, prompt:** depends on p_follow; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C4b looking from 29 Jul:** depends on c4b_breadth_arcmin; conditional on N3's fallback (zones alternating outward from the start place) and N2's.
- **C2 over the versions of N1:** depends on the Berlin sheet's limits.

T1 by 1 October is 1 in every case under this reading, because the 29 September sweep and
its noted disc are held. It is reported and not read.

Prior-dependent (§ 6.5, 20,000 accepted draws with every interval of § 4 halved, acceptance
0.183): C2 chart, RA of date, C3 from S4, C4a Le Verrier June. With the priors halved the chart case under right ascension of date moves
to "very probably", and the two starts without prompt comparison to "an open chance".

### R3. Paired differences against A0 (held reading)

| case | P(C) − P(A0) | statement |
|---|---|---|
| C1a 12 Aug against 30 Jul | 0.71 (0.50 to 0.97) | improves on A0 |
| C1b 4 Aug against 12 Aug | 0.28 (0.00 to 0.93) | no difference established |
| C1a, 7-day delay | 0.71 (0.50 to 0.97) | improves on A0 |
| C1b, 7-day delay | 0.28 (0.00 to 0.93) | no difference established |
| C2 chart, equinox 1800 | 0.00 (0.00 to 0.00) | no difference established |
| C2 chart, RA of date | 0.92 (0.77 to 1.00) | improves on A0 |
| C3 from S4 | 0.41 (0.04 to 0.81) | improves on A0 |
| C3 from S4, prompt | 0.79 (0.22 to 0.99) | improves on A0 |
| C3 from S5 | 0.00 (0.00 to 0.00) | no difference established |
| C3 from S5, prompt | 0.00 (0.00 to 0.00) | no difference established |
| C4a Le Verrier June | 0.41 (0.04 to 0.81) | improves on A0 |
| C4a, prompt | 0.85 (0.63 to 0.99) | improves on A0 |
| C4b looking from 29 Jul | 0.35 (0.02 to 0.90) | improves on A0 |

Since A0 is 0 under this reading, each difference equals the case's own probability, and
"improves" means only that the case's 5 % point is above 0. C1b has a third of its draws at
0, where the 12 August zone does not cover the 4 August place. The differences between
counterfactuals are in the JSON (`pairs_between_counterfactuals_not_read`), tabulated and
not read.

### R4. C1's completion date

Mean over draws of the probability that the comparison has flagged the planet by the date:

| case | 13 Aug | 14 Aug | 16 Aug | 20 Aug | 22 Sep | 95 % of the flags by |
|---|---|---|---|---|---|---|
| C1a 12 Aug against 30 Jul | 0.50 | 0.66 | 0.71 | 0.71 | 0.71 | 1846-08-15 |
| C1b 4 Aug against 12 Aug | 0.16 | 0.23 | 0.27 | 0.28 | 0.28 | 1846-08-16 |
| C1a, 7-day delay | 0.00 | 0.00 | 0.00 | 0.66 | 0.71 | 1846-08-21 |
| C1b, 7-day delay | 0.00 | 0.00 | 0.00 | 0.23 | 0.28 | 1846-08-22 |

The desk work is short against the six weeks available: a comparison begun the next day is
done by 16 August in nearly every draw, and a 7-day delay moves the date by a week and the
probability before 23 September not at all.

### R5. Tornado (§ 6.4), top three inputs per case

Change in the mean probability before 23 September, input from its posterior 5 % to 95 %
value; (c) marks an input read by conditioning.

| case | first | second | third |
|---|---|---|---|
| C1a 12 Aug against 30 Jul | p follow +0.43 | moving extent factor +0.10 | branch ra +0.07 (c) |
| C1b 4 Aug against 12 Aug | p follow +0.17 | weather +0.02 (c) | hours per session mean -0.02 (c) |
| C1a, 7-day delay | p follow +0.43 | moving extent factor +0.10 | branch ra +0.07 (c) |
| C1b, 7-day delay | p follow +0.17 | weather +0.02 (c) | hours per session mean -0.02 (c) |
| C2 chart, RA of date | p follow +0.22 | hours per session mean -0.00 (c) | weather +0.00 (c) |
| C3 from S4 | p disc routine +0.75 | n unknown sessions +0.06 (c) | hours per session mean -0.05 (c) |
| C3 from S4, prompt | n unknown sessions +0.29 (c) | p disc routine +0.24 | p follow +0.23 |
| C4a Le Verrier June | p disc routine +0.76 | p recorded +0.03 | hours per session mean -0.03 (c) |
| C4a, prompt | p follow +0.26 | p disc routine +0.19 | p recorded +0.02 |
| C4b looking from 29 Jul | c4b breadth arcmin -0.46 | p disc routine +0.33 | p disc looking +0.18 |

A0, the chart under the equinox of 1800 and the two S5 cases are 0 in every draw and have no
ranking.

### R6. Secondary, reported and not read: his procedure as a plan (every event redrawn)

| case | P before 23 Sep | T1 by 1 Oct | T2 before 23 Sep | P(C) − P(A0) |
|---|---|---|---|---|
| A0 actual procedure | 0.34 (0.03 to 0.70) | 0.78 | 0.34 | 0.00 (0.00 to 0.00) |
| C1a 12 Aug against 30 Jul | 0.81 (0.56 to 0.98) | 0.94 | 0.81 | 0.47 (0.15 to 0.83) |
| C1b 4 Aug against 12 Aug | 0.52 (0.05 to 0.96) | 0.84 | 0.52 | 0.19 (0.00 to 0.73) |
| C1a, 7-day delay | 0.81 (0.56 to 0.98) | 0.94 | 0.81 | 0.47 (0.15 to 0.83) |
| C1b, 7-day delay | 0.52 (0.05 to 0.96) | 0.84 | 0.52 | 0.19 (0.00 to 0.73) |
| C2 chart, equinox 1800 | 0.34 (0.03 to 0.70) | 0.78 | 0.34 | 0.00 (0.00 to 0.00) |
| C2 chart, RA of date | 0.94 (0.83 to 1.00) | 0.98 | 0.94 | 0.61 (0.28 to 0.92) |
| C3 from S4 | 0.41 (0.04 to 0.81) | 0.81 | 0.40 | 0.07 (-0.02 to 0.14) |
| C3 from S4, prompt | 0.79 (0.22 to 0.99) | 0.95 | 0.76 | 0.45 (-0.02 to 0.82) |
| C3 from S5 | 0.00 (0.00 to 0.00) | 0.67 | 0.00 | -0.34 (-0.70 to -0.03) |
| C3 from S5, prompt | 0.00 (0.00 to 0.00) | 0.67 | 0.00 | -0.34 (-0.70 to -0.03) |
| C4a Le Verrier June | 0.41 (0.04 to 0.81) | 0.81 | 0.41 | 0.07 (-0.01 to 0.14) |
| C4a, prompt | 0.85 (0.63 to 0.99) | 0.95 | 0.85 | 0.51 (0.24 to 0.83) |
| C4b looking from 29 Jul | 0.35 (0.02 to 0.90) | 0.77 | 0.35 | 0.01 (-0.26 to 0.64) |

Under the replay A0 is 0.34: the chance that a disc is noted on one of two routine recordings,
1 − (1 − p_disc_routine)², on draws that were kept because no disc was noted. The starts
from S4, S5 and Le Verrier's place are the same under both readings before 23 September,
since they share no recording with the history before 29 September.

### R7. Declared choices and deviations

1. **Tornado.** Run on the first 10,000 accepted draws. Inputs that are a number on the draw
   are set to their posterior 5 % and 95 % values with everything else as drawn. Inputs that
   shape the sessions (weather, the number of sessions, mean hours per session) cannot be set
   without redrawing the sessions, so they are read by conditioning: the mean over the
   lowest tenth against the highest tenth of the draws. The zone order (north or south
   first) is set; the zone branch (declination edge or right-ascension end) is conditioned.
   Conditioned rows carry the other inputs' correlations with them.
2. **Prior width.** 20,000 accepted draws, not 100,000. The zone priors and the branch
   weight are not halved; every interval of § 4's table and the C4b breadth is.
3. **p_follow is drawn once per draw and applied independently to each comparison.** Where a
   case has two comparisons that can flag the planet (the chart case under right ascension
   of date; a prompt comparison on two zones), the chance is 1 − (1 − p_follow)².
4. **Event seeds** are `[2026, 11, 4, 2, case, chunk]`, so that a resumed run reproduces the
   inner draws chunk by chunk.
5. **Compute:** 20 min on one core against the plan's "seconds" and the pre-run's
   estimate of 6 min. The estimate covered one reading with inner draws; the run had two
   readings, the prior-width run and the tornado. One ledger row `p1b_challis/draws`, local
   only. The site is untouched.
6. **The paired statement is weak by construction.** With A0 at 0, "improves on A0" does not
   distinguish a case at 0.35 from one at 0.85; the rows of R2 do.

### R8. What this does and does not show

It shows which of the named changes could have mattered and what each rests on. The
comparison he began and broke off at No. 39 is the cheapest: no new observation, a few days
of desk work, and its probability is p_follow times the chance that the July zone reached
the planet's August place. The start place matters less than the method: from S4 or from
Le Verrier's June place the planet is recorded within weeks, but a routine recording is
recognised only through a noted disc, on which the record itself (two recordings, no disc)
and the eight earlier recordings by others both point low. It does not show what Challis
ought to have done, and it is not a probability for the British search as a whole. Every
number is conditional on the printed record standing in for the observing books (N2), and
the chart case on the limits of a sheet not yet examined (N1).

## After review, 2026-11-12 (results above not edited)

### Deviation: the paired-difference rule is uninformative under the held reading

Identified after reading the results. Under the reading with the realised events held, A0 is
0 before 23 September by construction, so the pre-registered rule of § 7 ("Cᵢ improves on A0"
when the 5 % to 95 % range of P(Cᵢ) − P(A0) excludes 0) tests only whether the case's own
5 % point is above 0. The statements of R3 stand as the rule gives them and carry no
information beyond R2.

**Post hoc, description only, not read against any rule:** paired differences between
counterfactuals on the shared accepted draws (held reading, before 23 September).

| pair | mean | 5 % to 95 % | draws above 0 | below 0 | equal |
|---|---|---|---|---|---|
| C1a minus C1b (the two comparisons) | +0.43 | 0.00 to +0.96 | 0.59 | 0.02 | 0.39 |
| C1a minus C4b (comparison against looking from 29 July) | +0.37 | −0.27 to +0.85 | 0.82 | 0.18 | 0.00 |
| C3 from S4 with prompt comparison minus C1a | +0.08 | −0.50 to +0.37 | 0.86 | 0.14 | 0.00 |

In words: carrying the 12 August comparison through beats the 4 August comparison in 59 % of
draws and ties in 39 % (the draws where neither zone setting matters or both flag alike); it
beats the 29 September procedure started in July in 82 % of draws; a start from S4 with
prompt comparison is ahead of it in 86 % of draws by a small margin and behind by a large
one in the rest, so the mean difference is +0.08 with a range that spans 0. The shares and
ranges are in the JSON under `post_hoc_pairs_description_only`; the JSON's sha256 is now
`eeb8f4e8bfbd7bb52d85d77e2479bf6106a8808cff67b45d6ad193f441e70079`, the earlier figures in it unchanged.

### Prior-width check (§ 6.5): the rows that changed

Three readings change row when every interval of § 4's table is halved about its centre, and
are marked **prior-dependent**:

| case | row, full priors | range | row, halved priors | range |
|---|---|---|---|---|
| C2 chart, RA of date | depends on | 0.77 to 1.00 | would very probably have | 0.87 to 0.98 |
| C3 from S4 | depends on | 0.04 to 0.81 | an open chance | 0.33 to 0.71 |
| C4a Le Verrier June | depends on | 0.04 to 0.81 | an open chance | 0.32 to 0.71 |

No other reading changes row. Two limits of the check as run: the zone priors (the
declination setting and the branch weight) were not halved, and the check used 20,000
accepted draws. The tornado of § 6.4 used the first 10,000 accepted draws.

### Correction, 2026-11-12

In the post hoc table above, the 39 % of draws in which C1a and C1b are equal are the draws
where both comparisons cover the planet or neither does (each case's probability is p_follow
times whether its comparison covers the planet). This replaces the gloss given there.

## Named input N1: search of 12 November 2026

N1 is the Berlin sheet's limits, epoch and faintest plotted magnitude (§ 3). The sheet in
question is Hora XXII, Argelander's, the one Challis had; Galle's sheet on 23 September was
Hora XXI, Bremiker's. Sources were fetched on 2026-11-12 and are cached under
`data/cache/adams_challis/berlin/` (ignored), with hashes in
`data/historical/adams_challis/MANIFEST.md`.

| what was found | where | class | URL or shelfmark | SHA-256 (first 16) | what it states |
|---|---|---|---|---|---|
| The sheet Hora XXII, lower part (declination −8° to −15°, the key, the title), photograph of the Potsdam copy | Krajnović 2021, Fig. 6.10 (preprint) | the primary sheet in reproduction | https://arxiv.org/pdf/2108.06305, p. 28; copy: Library of the Leibniz-Institut für Astrophysik Potsdam | `cd62724e87aa1859` (PDF) | Title: "Theil des Himmels, zwischen XXIIʰ und XXIIIʰ der geraden Aufsteigung, und 15° südlicher bis 15° nördlicher Abweichung, für 1800 … entworfen von Herrn Professor Fr. Argelander in Åbo". Grid lines every 4 minutes, lettered XXII, 4′ … 56′, XXIII, with one more column outside each: the outermost lines are 21h 56m and 23h 04m. Lowest line 15°. Key of magnitudes 1 to 9, the faintest class "9-10". Engraved "K. Kolbe sculps. Berlin" |
| The Academy's catalogue, Vorbericht p. VI, the rules for the sheets | *Cataloge zu den 24 Stunden der akademischen Sternkarten für 15° südlicher bis 15° nördlicher Abweichung, bearbeitet von verschiedenen Astronomen, herausgegeben von der Königlichen Akademie der Wissenschaften zu Berlin*, Berlin 1859 | primary | Wrocław University digital library, shelfmark K178, https://bibliotekacyfrowa.pl/dlibra/publication/113969 (images: https://glam.uwr.edu.pl/iiif/IAS_K178_58271/manifest, image 12) | `7a2bc356b6e110d2` | Rule 1: "Das Netz besteht aus Quadraten für die Grade der Aufsteigung und Abweichung … es geht von 4 Zeitminuten vor dem Anfange einer Stunde, bis zu 4 Zeitminuten nach dem Ende derselben, und enthält also 510 Quadrate"; the sheets run "von − 15° bis + 15° der Abweichung" |
| The same, p. VII | as above, image 13 | primary | as above | `d73e6f9c9ddf74e6` | Rule 2: the stars are drawn "für den Anfang des Jahres 1800". Rule 5: the sheet is completed with "alle Sterne … bis zu der aus dem Probeblatte hervorgehenden Grenze", those that "ein Fraunhoferscher Kometensucher von 34 Linien Öffnung und zehnmaliger Vergrößerung unter vortheilhaften Umständen noch sehen kann" |
| The same, p. IV, Bessel's letter to the Academy | as above, image 10 | primary | as above | `34de3b35c1ae1dcc` | Proposal 1: "Die Karten müssen in Blättern von 17° Breite und 30 Grad Höhe, für die Epoche von 1800, immer von Stunde zu Stunde der Rektascension und auf jeder Seite ein Grad mehr … entworfen" |
| The catalogue's title page | as above, image 7 | primary | as above | `434352692928ef5c` | The title as given above; "Berlin, gedruckt in der Druckerei der Königlichen Akademie der Wissenschaften, 1859. In Commission in F. Dümmler's Verlags-Buchhandlung". No editor is named on the title page; Jones 2002 says the introduction is signed by Encke and Dove (not read here) |
| Table of the 24 sheets | Jones 2002, *Highlights of Astronomy* 12, 367–370, DOI 10.1017/S1539299600013770 | secondary | publisher's PDF | `87b3f2c80bd968a6` | Hora XXII: Argelander in Åbo, 1832, 1436 catalogued stars and 65 more observed with a meridian circle. Nothing on limits or epoch |
| The design in a modern account | Krajnović 2021, pp. 10–11 of the preprint | secondary | as the first row | as the first row | 4-minute overlaps, 510 squares, epoch 1800, citing the 1859 catalogue pp. IV–V; the fallback of § 3 came from here. Fig. 6.10's caption places Neptune "between 21h 56m and 22h" on the sheet in August 1846, which holds for right ascension of date, not for the sheet's own equinox |

Not found or not reached: a scan of the whole sheet (the upper part, +15° to −8°, was not
seen; the title states the range); the catalogue's section for Hora XXII itself (it is in the
same Wrocław volume, 598 images, not paged through); Bessel 1826, *Astronomische Nachrichten*
4, 437 (the archive.org scan found was volume 5). Google Books' search interface refused with
"quota exceeded" and was left. e-rara, the Staatsbibliothek zu Berlin, the BBAW, Heidelberg,
Dresden, Munich, Gallica and HathiTrust were not searched once the sheet and the catalogue
had been read.

**N1 obtained 12 Nov 2026.** The sheet and the Academy's rules agree: Hora XXII runs from
21h 56m to 23h 04m of right ascension and from −15° to +15° of declination, for the equinox
of the beginning of 1800; its faintest class of stars in the key is 9-10. These are the
limits of N1's first fallback, now read from the sheet. The second fallback (the same limits
applied to right ascension of date) is the frame mismatch that § 3 said it would be if the
sheets were for 1800.

**C2 under N1 itself**, as § 5 provides, on the same accepted draws (seed (2026, 11, 4, 1),
100,000 accepted of 1,003,853; `uv run python -m fork.export.p1b_challis_n1 run`, 59 s;
`docs/notes/P1B_challis_detection_n1.json`; the step 4 JSON is unchanged, SHA-256 recorded in
the new file):

| case | held, before 23 Sep: mean (5 to 95 %) | § 7 row | by 1 Oct | halved priors | paired against A0 | replay (not read) |
|---|---|---|---|---|---|---|
| C2 chart, N1 itself | 0.00 (0.00 to 0.00) | would very probably not have | 1.00 | 0.00, same row | no difference established | 0.34 (0.03 to 0.70) |

Checks: N1's chart test equals the first fallback's on every night; the accepted draws are
those of step 4 (same acceptance counts; A0 equal draw by draw, largest difference 0); the
per-draw values equal the first fallback's (largest difference 0). So the number is the first
fallback's by construction; what changes is its footing. The event-sampling check of the
closed form was made in step 4 on those values and is not repeated.

**Reading (§ 7), N1 obtained 12 Nov 2026.** C2 under N1: would very probably not have; 0 in
every draw, no input moves it; no longer conditional on a fallback. The reason is the one
step 2 found: for the equinox of 1800 the planet was at 21h 55m 43s on 4 August and
21h 54m 54s on 12 August, west of the sheet's edge at 21h 56m, so a check of those nights'
stars against Hora XXII could not have met it. (On 29 and 30 July it was still on the sheet,
at 21h 56m 19s and 21h 56m 13s, within 19 s of the edge; C2 as planned checks the stars of
the two nights on which the record has the planet, 4 and 12 August.) By the letter of § 7,
"C2 is read once per version of N1", and the three versions do not fall in one row, so step
4's combined line "depends on the Berlin sheet's limits" is left as written; the version
that the sheet supports is the first. The companion page is not changed here.

Not settled by this: the sheet's faintest magnitude in modern terms (the key's last class is
9-10; Neptune was about 7.8, so the planet's absence from the sheet is not in question), and
whether the Cambridge copy was the same printing as the Potsdam one.

**N1, the author's verification, 12 Nov 2026.** The author checked Hora XXII against the
sheet images by eye on 12 November 2026. Read from the sheet:

| item | reading |
|---|---|
| title cartouche | "zwischen XXII^h und XXIII^h … 15° südlicher bis 15° nördlicher Abweichung, für 1800" |
| author | Argelander (Åbo) |
| plotted field | 21h 56m to 23h 04m: a 4-minute overlap column beyond each hour line |
| southern edge | −15° |
| key | down to class 9-10 |

The engraver's line is under the lower right corner. In the reproduction used for N1 its
letters are about 8 px high and the name cannot be confirmed letter by letter; the result
file's reading, "K. Kolbe sculps. Berlin", stands.

The nominal hour in the title (22h to 23h) is narrower than the plotted field. N1 uses the
plotted field, 21h 56m to 23h 04m, as above; the C2 figures and the § 7 reading are not
changed by this. With the nominal hour instead, the planet (21h 55m 43s and 21h 54m 54s for
1800 on 4 and 12 August) is further outside the sheet, so the reading would be the same.

## Correction, 13 November 2026: the withdrawn chart case, and one reference

- The case "C2 chart, RA of date" (0.92, 0.77 to 1.00; replay 0.94) was withdrawn in v1.2.1
  (§ "Named input N1" above): it applied the sheet's limits, which are for the equinox of 1800,
  to right ascensions of date. It still stands in the tables of § 7 and of the post hoc
  sections above, and those tables are left as written. Under N1 the chart case is 0.00,
  "would very probably not have" (`P1B_challis_detection_n1.json`, `case.held.primary`).
- The figure `P1B_challis_detection.png` still draws a panel for that case. It is not
  redrawn; this line is its correction, and the label in `src/fork/export/p1b_challis.py`
  carries the same note.
- Le Verrier's paper of 1 June 1846 is *Comptes rendus* 22, 907–918; the place quoted
  above is on p. 917. The citation "Comptes rendus 22, 917" above means that page.
