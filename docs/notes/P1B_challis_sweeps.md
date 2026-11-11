# Phase 1b, step 2: Challis's search log, reconstructed from the printed record

Date: 2026-11-10. No modelling. Deliverables in `data/historical/adams_challis/`:

- `challis_sweeps.csv`: one row per observing session (29 and 30 July, 4 and 12 August, the
  comparison of about 13 August, the two months after 12 August, 29 September, 30 September,
  1 October) plus the plan as Airy stated it, the plan as Challis modified it, and the total.
- `challis_sweeps_statements.csv`: the 111 printed statements the rows rest on, each with its page
  image, from Challis's account (Memoirs RAS 16, 415–426), the Society's abstract (MNRAS 7, 145–149)
  and Airy's account (MNRAS 7: Airy to Challis 9 July, the Suggestions of 13 July as Airy describes
  them, Challis to Airy 18 July, Challis to Main 7 August, Challis to Airy 2 September and 12 October).
- `neptune_on_challis_dates.csv` and `neptune_track_1846.csv`: the computation (below), made by
  `docs/notes/scripts/p1b_neptune_on_challis_dates.py`.

Transcription: two independent passes, compared on every numeric token, no disagreement; the
differences between printings (80′ against 70′ for the sector zone of 4 August, power 166 against 160,
"nearly 300" against 300, "above 3000" against 3000) are listed with the deciding images in
`challis_sweeps_reconciliation.md`. Every statement is primary (P); the session verdicts that rest on
the computation are marked D (derived).

## The log

| date 1846 | method | zone as printed | depth | stars | planet | Neptune in the zone |
|---|---|---|---|---|---|---|
| 29 July | trials, mostly telescope fixed | directed to Adams's 325° place; declination setting as on 30 July and 4 Aug | — | — | no | not determinable (RA extent not printed) |
| 30 July | telescope moving | zone 9′ broad | 11th mag | — | no | outside (P: "wandered into the zone between July 30 and August 12") |
| 4 Aug | sector, polar axis fixed | 80′ in NPD (70′ in the abstract and letter), brighter stars | bright | — | yes | inside (observed) |
| 12 Aug | telescope fixed | the 30 July zone, longer in RA | 11th mag | ≥ 49 | yes, No. 49 | inside (observed) |
| ~13 Aug | comparison | 30 July against 12 Aug to No. 39 | | 39 | missed by 10 | — |
| 13 Aug – 28 Sept | principally fixed, partly moving | "earlier in right ascension" | 10–11 mag | some thousands | no | outside (P) |
| 29 Sept | sweep in declination within Le Verrier's limits | breadth "considerable", not quantified | 11th mag; disks looked for | nearly 300 in 3½ h | yes, "seems to have a disk" | inside (observed) |
| 30 Sept | none (moon) | | | | | |
| 1 Oct | news of the discovery | | | | | |

The plan (Airy, 13 July): a zone centred at longitude 325° on the ecliptic, 15° each way and 5° of
latitude each way, divided into 1° declination intervals with the RA extent marked out for each,
power not under 120, zones 15′ broad, transits at one wire and NPD from horizontal wires 3′ apart,
three complete sweeps, 80 zones of an hour each, 300 hours. Challis's modifications: power 166 and
zones 9′ broad, declinations read on the micrometer comb (17″ per revolution, to a fourth), every
star visible in the illuminated field (11th magnitude), a first sweep with the telescope fixed and a
second with the telescope moving under clock drive, moonlit nights fixed only, clouds marked by a
line in the memorandum-book, and reference stars taken with the 80′ sector.

## The computation (D)

Neptune's geocentric place on each date is carried from the cached Horizons states of Neptune and
the Earth at 1846-09-23 23:00 TDB (the one 1846 epoch cached for both bodies) by two-body Kepler
arcs; the perturbations of Jupiter and Saturn over the 90 days spanned move the direction by well
under 1″. Light time, annual aberration and nutation are applied for the apparent place. No Horizons
fetch is made. Check against Challis's three reduced places (apparent, equinox of date):

| date | Challis | computed apparent | Challis − computed |
|---|---|---|---|
| 4 Aug, 13h 36m 25s GMT | 21h 58m 14.70s, NPD 102° 57′ 32.2″ | 21h 58m 15.0s, 102° 57′ 29″ | −0.3 s, +4″ |
| 12 Aug, 13h 3m 26s | 21h 57m 26.13s, 103° 2′ 0.2″ | 21h 57m 25.9s, 103° 1′ 58″ | +0.2 s, +2″ |
| 29 Sept, 10h 10m | 21h 52m 48.9s, 103° 26′ 30″ | 21h 52m 49.2s, 103° 26′ 38″ | −0.3 s, −8″ |

Challis's stated errors were 3″ and 4″ in NPD for the two August places; the 29 September place was
"not given as having any pretension to accuracy".

Against Adams's guidance ephemeris for 325° (the column Challis chose), Neptune stood 8.3 minutes of
time (2.0°) east of the tabulated right ascension and 10′ north of the tabulated declination at the
end of July, the offsets falling to 7.1 minutes and 3′ by 29 September. A 9′ zone centred on the
tabulated declination would have excluded the planet until mid September, but Challis's declination
settings are not printed, so this is not a verdict on 29 July.

Neptune's mean place of date moved from 22h 01.1m, −12° 49′ (1 July) to 21h 51.8m, −13° 31′
(15 October), westward about 6 s of RA a day; in the mean equinox of 1800, the epoch of the Berlin
charts, the same track is 21h 58.6m to 21h 49.3m, declination −12° 55′ to −13° 45′.

## The Berlin charts

The chart limits are secondary (one hour of RA per sheet with a 4-minute overlap on each side,
declination ±15°, epoch 1800, from Krajnović 2021; no sheet or the 1859 catalogue has been seen).
Under those limits Hora XXII covers 21h 56m to 23h 04m and Hora XXI 20h 56m to 22h 04m, both in the
equinox of 1800. Dates on which Neptune lay within each (D):

| assumption | Hora XXII (Challis had it) | Hora XXI (Bremiker's, in Berlin) |
|---|---|---|
| limits in equinox 1800 with the 4-minute overlap | until 2 August 1846 (crossed 21h 56m on 2–3 August) | the whole period |
| limits in equinox 1800 without overlap (22h–23h) | never during the search (Neptune was west of 22h before 1 July) | the whole period |
| limits applied to the RA of date (a frame mismatch) | until 26 August | the whole period |

Consequences, each conditional on the secondary limits: on 4 and 12 August, the two nights Challis
recorded the planet, it was already 11 s and 60 s of right ascension (2.7′ and 15′ of arc) west of
Hora XXII's edge; a comparison of those observations with the sheet he had would have met the
planet at or just beyond its margin, not inside it. Kollerstrom's "Neptune moved from Hour 22 onto
Hour 21 towards the end of August" matches the date-frame row and is therefore in the wrong frame if
the sheets are in the equinox of 1800, as the secondary sources state. The addendum to
`P1B_sources.md` corrects that note's earlier remark, which made the same comparison in the frame of
date. Hora XXI, which Challis did not find in the University Library, covered the planet on every
night of the search.

## What the printed record cannot determine

1. The zone boundaries of any session: no right-ascension limits are printed for any night, and no
   declination setting; only breadths (9′, 70′ or 80′, "a considerable breadth") and relations
   ("same position as July 30", "greater extent in RA", "earlier in right ascension").
2. The dates of the sessions between 12 August and 28 September, their number, and the division
   between fixed and moving nights.
3. Which edge of the 30 July zone the planet crossed before 12 August.
4. The per-night star counts other than 29 September (nearly 300) and 12 August (at least 49), and
   the magnitudes recorded.
5. The limits of the 29 September sweep ("within the limits of longitude marked out by
   M. Le Verrier"), which Challis does not quantify.
6. The comparison-star identities for the zones, hence the absolute positions of the zones.
7. The Suggestions themselves (the declination intervals and their RA extents) and Airy's note of
   21 July; both are in RGO 6/96a.
8. The Berlin sheets' actual limits and epoch, needed before any statement about the planet being on
   Hora XXII is more than conditional.

## What the observing books would add

Challis's memorandum-books (Cambridge University Library, UA/Obsy; the Northumberland observations
series) hold, per night: the date, the method, the declination-rod and sector readings, the
chronometer times and hour-circle readings of every star, the magnitudes, the reference stars, the
lines drawn at cloud interruptions, and the comparisons made. From them the zone of each session can
be placed on the sky to the minute of arc, the sessions after 12 August can be dated and counted, the
3000 positions can be reduced and compared with modern catalogues to give the actual completeness
of each sweep to the 11th magnitude, and the question of whether the planet fell inside the
29 July zone can be answered. Airy's papers (RGO 6/96a) would add the Suggestions, the 21 July note on
the Berlin maps, and the five columns of Adams's ephemeris.
