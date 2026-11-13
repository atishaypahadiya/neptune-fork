# Phase 1b, step 1: Adams's solutions, transcribed

Date: 2026-11-10. No fitting. Deliverable: `data/historical/adams_challis/adams_solutions.csv`
(one row per dated solution) with `adams_solutions_statements.csv` (every printed statement the rows
rest on, 44 statements, each with its page image). The two transcription passes and their
reconciliation are `adams_solutions_pass1.csv`, `adams_solutions_pass2.csv` and
`adams_solutions_reconciliation.md` in this directory. Sources and page renders are those of
`P1B_sources.md` and `data/historical/adams_challis/MANIFEST.md`, plus Sampson (1904), downloaded for
this step (ADS 1904MmRAS..54..143S, added to the MANIFEST).

Marks: **P** primary-confirmed, read on the page image of an 1846–47 printing, or of Sampson's 1904
catalogue of the St John's manuscripts, whose facsimile plates are the manuscripts themselves;
**S** secondary-only; **D** derived here by arithmetic from primary numbers (stated as such; nothing
derived is written into the solutions table except the two mean distances in AU flagged in their
cells).

## Method

1. Pass 1: the page renders were read directly (Sampson pp. 143–170 and the twenty plates; Airy's
   letters Nos. 11 and 20 in MNRAS 7 and again in Memoirs 16, pp. 395 and 406–407; the abstract
   MNRAS 7, 151–152; the memoir §§ 30–31, 47–48, 53–54 in the separate printing; Challis MNRAS 7, 148).
2. Pass 2: an independent reading of the same renders under a written brief, without access to pass 1
   or to the repository.
3. The passes were compared on every numeric token. There was no disagreement; the differences
   between printings of the same solution were each settled on the image named in the
   reconciliation note.
4. Sampson's six determinations were taken as the backbone (his Table, p. 151, epoch 1810.328); each
   was dated from his catalogue and tied to the published statements.

## The solutions

| id | Sampson | date | assumed a | e | ϖ | mean longitude | mass | mark |
|---|---|---|---|---|---|---|---|---|
| S1 | I | 1843-09-30 (sheet date) | 38.4 | 0 | — | 253° 28′ at 1810.328; 304°.0 at the 1843 opposition | 0.00023655 | P (plate and catalogue) |
| S2 | II | undated, 1843–45 | 38.4 | 0 | — | 263° 16′ at 1810.328 | 0.000410 | P (catalogue) |
| S3 | III | 1845-04-28 (note), comparison 1845-05-19 | 38.4 | 0.19717 | 248° 8′.3 | 271° 37′ at 1810.328; "Long. 1845 = 324° nearly" | 0.0001573 | P (catalogue) |
| S4 | IV | 1845-09-18 (θ); the September 1845 paper to Challis | 38.4 | 0.1428 | 320° 30′ | 267° 58′.2 at 1810.328; 321° 40′ at end Sept 1845; geocentric 320° 30′ "diminishing about 1′ daily" | 0.000173 | P (Sampson's transcription of the sheet) |
| S5 | V | 1845-10-21 (to Airy; "last days of October") | 38.4 (a/a′ = 0.5) | 0.16103 | 315° 57′ (315° 55′ in the Oct 1845 paper) | 323° 34′ at 1 Oct 1845; 325° 8′ at 1 Oct 1846; 325° 7′ at 6 Oct 1846 | 0.0001656 | P (four printings) |
| S6 | VI | 1846-08-20 (θ) to 1846-09-02 (letter) | 37.25 (a/a′ = 0.515) | 0.120615 | 299° 11′ | 323° 2′ at 1 and 6 Oct 1846 | 0.00015003 | P (four printings) |
| S7 | — | 1846-09-02 | 1/0.57 (≈ 33.7, D) | "very small" | — | about 315° 20′ at 1 Oct 1846 | — | P (letter No. 20; memoir § 53) |
| S8 | — | 1846-11-12 | 30.05 (from observations) | — | — | heliocentric 326° 39′ on 4 Aug 1846; node 309° 43′; i 1° 45′ | — | P (Challis) |

S7 is not a solution: no elements were computed. Adams took the ratio of the 1843 errors of S5 and S6
(6.84″ to 5.50″; 5 : 4 over 1843–45 in the memoir), inferred a/a′ = 0.57, and extrapolated the mean
longitude to "about 315° 20′, which I am inclined to think is not far from the truth". The memoir
(§ 54, 12 November 1846) says the direction of this extrapolation was wrong: with the mean motion
properly allowed for, the mean longitude at the present time "would probably soon arrive at a
minimum, and afterwards begin again to increase", and "I hastily inferred … that the effect of a
diminution in the mean distance would be to diminish the mean longitude". S8 is an orbit from Challis's
positions of the discovered planet, not a perturbation solution.

Sampson's own comparison (p. 152, P) converts I–VI to heliocentric longitudes for 1845 against Neptune's
323°.1: I 307°.0 (+16.1), II 316°.8 (+6.3), III 347°.2 (−24.1), IV 321°.0 (+2.1), V 325°.0 (−1.9),
VI 326°.0 (−2.9). His signs are (Neptune − predicted).

## Verification of the "final 2 Sep 1846 solution ≈ 315⅓°, ≈ 12° from Neptune" claim

Primary facts (letter No. 20, Adams to Airy, 2 September 1846, MNRAS 7, 137–139 and Memoirs 16,
405–407; memoir §§ 53–54):

- The letter communicates **two sets of elements**, Hypothesis I (S5) and Hypothesis II (S6), and then
  the extrapolation S7. The last *computed* solution before the discovery is S6: mean longitude 323° 2′
  on 1 October 1846, perihelion 299° 11′, e 0.12062.
- The figure 315° 20′ (= 315°.33, "315⅓°") is primary. It is a **mean** longitude for 1 October 1846,
  given for a/a′ = 0.57 with the remark that the eccentricity "would be very small"; it is an
  inference, not a solution, and Adams withdrew its direction in print ten weeks later.
- Neptune's heliocentric longitude on 1 October 1846, mean ecliptic and equinox of date, from the JPL
  Horizons vectors cached in this repository (interpolated between the cached epochs of 23 September
  and 1 October 1846): 327° 0′ (D). Kollerstrom gives 326° 57′ for 23 September (S); the two agree.
- Distance of 315° 20′ from 327° 0′: **11° 40′** (D). With e → 0, mean and true longitude coincide, so
  the comparison with the planet's true longitude is fair.

Verdict: the number and the date are correct and primary; the characterisation is not. 315° 20′ was
not Adams's final solution. It was an extrapolation beyond his final solution, offered with a hedge
and retracted in the memoir. Adams's final computed solution (S6) put the planet at a mean longitude
of 323° 2′, 4° 0′ short of Neptune's 327° 0′ in mean longitude, and at a true longitude of 329° 25′
(D, from the printed elements: equation of centre +6°.44 at anomaly 23° 48′), 2° 28′ beyond it.
Hypothesis I (S5) gives a true longitude of 328° 41′, 1° 44′ beyond. Any study that uses "Adams's
prediction" has to say which of S4, S5, S6 or S7 it means; the four differ by 14° in longitude and
by 5 AU in distance.

Derived true longitudes (D), for the cross-check below and nothing else: the mean longitude printed
for 6 October 1846 is carried to 23 September 1846 with the mean motion of the assumed distance
(238.0 yr for 38.4 AU, 227.3 yr for 37.25 AU), the equation of centre is taken to third order in e,
and Neptune is 326° 57′ on that date.

| | mean longitude 23 Sept 1846 | equation of centre | true longitude | Neptune − predicted |
|---|---|---|---|---|
| S5 | 325° 4′ | +3° 37′ | 328° 41′ | −1° 44′ |
| S6 | 322° 59′ | +6° 26′ | 329° 25′ | −2° 28′ |
| S7 | 315° 18′ | 0 | 315° 18′ | +11° 39′ |

## Cross-check against Kollerstrom (2006, Part 1, Appendices I–II and text), secondary

Agreements (his figure, then the primary): Adams 1 mean longitude 323° 34′ at 1 Oct 1845 (No. 11);
perihelion 315° 55′ (No. 11); eccentricity 0.161; mean distance 38.4; period 237.6 yr (No. 11's
1° 30′.9 per year gives 237.6); equation of centre 3° 38′ and true longitude 328° 41′ at discovery
(D above: 3° 37′, 328° 41′); error 1° 44′. Adams 2 mean longitude 323° 2′ at 1 Oct 1846 (No. 20);
perihelion 299° 11′; eccentricity 0.121 (0.12062 rounded); anomaly 23° 48′; equation of centre 6° 27′
(D: 6° 26′); true longitude 329° 27′ (D: 329° 25′); error 2° 30′. His statement that the 2 September
letter also gave a Hypothesis I mean longitude of 325° 8′ for 1 October 1846 is primary. His 315° 20′
passage with "which I am inclined to think is not far from the truth" is quoted correctly, and his
33.6 AU for a/a′ = 0.57 agrees with 33.7 (D).

Differences:

1. **Hypothesis II mean distance: 37.5 AU (Kollerstrom App. I) against 37.25 (Sampson table VI) and
   1.942 × Uranus (abstract) = 37.25–37.3.** His own period for Adams 2, 227.4 yr, is the period of
   37.25 AU, not of 37.5 (229.6 yr), so the table is internally inconsistent; his perihelion distance
   33.0 and discovery distance 33.4 descend from the 37.5.
2. **Anomaly of Hypothesis I: 11° in App. I against 9° 7′ in App. II** (the second is right for
   325° 2′ − 315° 55′).
3. **Dating of Hypothesis I "of uncertain date" (App. II).** Sampson (p. 168) dates the communication
   to Airy 21 October 1845 and Airy "one of the last days of October, 1845"; the solution was complete
   by 18 September 1845 apart from the 2′ correction of the perihelion (S4 → S5).
4. **"Adams MS E VII, p. 15: eccentricity, apse position and mean longitude were here cited for the
   hypothetical planet for the first time" (note 80).** Sampson's E.VII belongs to the second
   hypothesis of 1846, and his E.II (April–May 1845) and E.III (the September 1845 sheet of notepaper)
   already state eccentricity, perihelion and mean longitude. The sheet gives, in addition, a
   geocentric longitude and its daily rate, which is what an observer needs; whether Challis
   received that sheet is Sampson's inference ("very likely").
5. **Epoch of the September 1845 solution "3 May 1810" (note 77).** Sampson's epoch is 1810.328, the
   opposition of Uranus, which is 30 April 1810 (0.328 × 365.25 = 119.8 days after 1 January).
6. **"the solutions he communicated after this date [June 1846] ranged over twenty degrees of
   celestial longitude" (p. 2).** After June 1846 Adams communicated one solution (S6, 323° 2′), one
   extrapolation (S7, 315° 20′) and the July ephemeris for Challis, which tabulates positions for
   assumed heliocentric longitudes 315°–335° at 5° steps (Challis, Memoirs 16, 421). The 20° range is
   the ephemeris's grid of assumptions, not a spread of solutions.
7. **dioi.org "Adams Dated Computations": "Two such solutions are found in his manuscripts."**
   Sampson's catalogue lists six determinations (I–VI) with dates for five of them. The "six
   solutions" of the Phase 1b plan is therefore primary via Sampson; the two the page names
   (18 September 1845, −50° 34′; 20 August 1846, −42° 52′) are the θ values of S4 and S6, and the
   second is printed by Sampson as −42° 51′ "very nearly".
8. **Krajnović 2016** gives no elements and no count; nothing to compare beyond the search dates
   already noted in `P1B_sources.md`.

Kollerstrom's reading of the September 1845 solution (note 77, "the angular distance in longitude
between the Uranus and 'hypothetical planet' mean motions … for 3 May 1810") matches Sampson's θ
convention, θ = ε − ε′, the difference of the mean longitudes of Uranus and the planet at the epoch;
Adams's memoir uses the same quantity (§ 30: ε = 217° 55′, θ = −51° 30′, ε′ = 269° 25′).

## What this step does not settle

- The manuscripts themselves (St John's College) have not been seen; everything from them is through
  Sampson's catalogue and his twenty plates. The plates reproduce C.II (1843) in full and the sheet of
  notepaper is transcribed, not reproduced.
- The five columns of the July 1846 ephemeris (heliocentric longitudes 315°–335°) are in the Airy
  papers (RGO 6/96a) and are not transcribed here; only the 325° column is in print.
- Neptune's longitude is from the cached Horizons vectors at two bracketing epochs; the interpolation
  error over 8 days is below 1′.

## Amendment: manuscript check, 12 Nov 2026

Added after the step was closed; nothing above is changed, and `adams_solutions.csv` is not
edited. St John's College supplied images of seven sheets of MS W16 on 12 November 2026
(`data/historical/adams_challis/MANIFEST.md` lists them with the rights condition; the
images are not in the repository). The first bullet of "What this step does not settle" is
now settled for these seven sheets only.

**Method.** Bands were cut from the full-size images. A first entry was made from them. A
second entry was made from the bands and a list of labels alone, without the first entry,
the committed tables or Sampson. 200 cells were compared. The first cut of three bands lost
nine cells at an edge; the second entry read those from recut bands. The label list names
the two trial values of θ on the 1846 sheet, so those two cells are not independent
readings.

| | |
|---|---|
| cells compared | 200 |
| differences | 4 |
| of which the same reading written differently | 1 |
| settled by the sheet's own arithmetic | 1 (−23.27, not −23.77) |
| left doubtful | 2, both in the Greenwich letter's table, neither a result of Adams |

**What the manuscripts confirm** (`adams_solutions_ms.csv`, 11 values that agree):

- S4, the note of results (E III 14): mean distance 38.4, eccentricity 0.1428, longitude of
  perihelion 320° 30′, mean longitude about the end of September 1845 321° 40′, geocentric
  longitude 320° 30′ nearly and diminishing about 1′ daily, mass 0.000173. Sampson's
  transcription is exact in all of them. The sheet carries no date; the date is on E III 6.
- S4's date: E III 6 ends "Septr 18th 1845" with θ = −50° 34′, as Sampson p. 166.
- S6's first date: E VII 15 ends "August 20th 1846" with θ "very nearly = −42° 51′", as
  Sampson p. 167.
- The corrections of the radius vector for 1834, 1840 and 1846 on F II 6 (+0.004923,
  +0.006962, +0.008250) are the Hypothesis II figures of the letter to Airy of 2 September
  1846, digit for digit.

**What they contradict.** Nothing in the committed tables.

**What they add** (9 values):

- E III 6: the first approximation was θ = −51° 10′, where the equation left −2.8964.
- E VII 15: the equation's sum was −28.866 at θ = −45° and −2.077 at −43°. A straight line
  through the two crosses zero at −42° 50′.7, which is his −42° 51′. On this sheet
  θ − α = −95° 16′ at θ = −45°, so α = 50° 16′.
- F II 6 is dated "1st Septr 1846". The radius-vector corrections of Hypothesis II were
  therefore computed the day before the letter that communicates them.
- An undated sheet headed "Calcn of Errors for 1843 & 1846": sums −109″.76 and −148″.65,
  three totals −91″.65, −124″.23 and −159″.76, and corrections of heliocentric longitude
  −83″.59, −113″.55 and −146″.92. The sheet does not say which hypothesis it belongs to, and
  its first column has no year. Not attributed here.
- A VII 3 is a letter, not a calculation: Robert Main at Greenwich to Adams, "1846, Septr 7",
  forwarding "the Heliocentric Errors of Uranus for the year 1845", six dates, +7″.79 to
  +9″.55.

**Does 7 September 1846 change S6's dating?** No. The sheet of that date is the Greenwich
letter above. S6 stays 20 August to 2 September 1846, and F II 6 puts a dated Hypothesis II
calculation at 1 September inside that range. What the letter does bear on is S7 and the
memoir. The letter to Airy of 2 September rests its extrapolation on the errors of 1843
alone; the memoir (§ 53) uses 1843 to 1845. The 1845 errors reached Adams on or after
7 September, so the memoir's form of that argument was written after the letter of
2 September. This is an inference from two dates, not a statement on any sheet.

**Not done.** The trial with the 1843 sheet (S1) and S2, S3 and S5 have no manuscript image
here. The many logarithms on the sheets were not entered, only the results and the lists
that carry a total. No person has spot-checked the entries. Both entries share one method
of reading, so their agreement is weaker evidence than two unrelated readers would give.

### Spot check by the author, 13 Nov 2026

Added below the amendment of 12 November, which stands as written. This section supersedes
two of its statements: the table's row "left doubtful: 2" and the sentence "No person has
spot-checked the entries" under "Not done".

I checked 20 cells of the reconciled entries against the band images
(`docs/review/stjohns_spotcheck_human_results.json`, exported 2026-11-13 04:40 UTC, seed
20261112, first-entry file SHA-256 `3de30e11…`). All 20 match.

| how the cell was chosen | cells | match |
|---|---|---|
| the two entries differed | 3 | 3 |
| flagged doubtful in an entry or in the comparison | 4 | 4 |
| one per image not otherwise covered | 4 | 4 |
| random | 9 | 9 |

**Resolved, no longer doubtful.** Every reading the entries disputed or doubted is confirmed:

- the three differences in values: −23.27 (1846 column, term 4; the second entry's reading,
  which the sheet's sum had already settled), and the first entry's +114″.88 (Main's table,
  13 Jan 1845, error of geocentric longitude) and +2505 (18 Sep 1845, coefficient 4). The
  last two are the "left doubtful" of the table above; `adams_solutions_ms_reconciliation.csv`
  now gives them as resolved. The fourth difference was the same reading of the letter's
  head written two ways and needed no check.
- the cells flagged doubtful: in Main's table, coefficient 4 −9486 and the heliocentric
  error +7″.79 for 13 Jan 1845, and the geocentric error +124″.05 for 14 Dec 1845; on
  E VII 15, 2θ − α = −136° 16′ at θ = −43°.

**Scope.** 11 of the 20 cells were chosen, not drawn: they are the cells above. 9 were
random among the remaining cells. No discrepancy in 9 random cells bounds the cell error
rate only weakly: the two-sided 95 % upper limit is 1 − 0.025^(1/9) ≈ 34 %. The check is
therefore not an estimate of the error rate of the 200 cells. Its value is that every
disputed or doubtful reading was confirmed by a person.
