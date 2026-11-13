# Phase 1b, step 3d: is the shift carried by the reference Adams subtracted, or by the observations?

Pre-registered plan. Written 2026-11-13 after step 3c's results and the transcription of
Bouvard's printed formula (`P1B_bouvard_1821.md`, route A). Nothing below §§ 1–9 runs until
the plan is committed and pushed, and the pre-registered sections are then byte-identical
(amendments and results only below the lines). Code and tests cite this note.

## Pre-registration (before any code or run)

### 0 Where steps 3b and 3c left it

At his assumed distance (α = 0.5), noise-free data put Adams's planet about 10° short of
Neptune at 1846.762, and his real data put it about 1° past. The numbers, each method
against Neptune in its own frame (step 3c, R2):

| hypothesis | method | clean data, O_clean | his real data, O_real | shift Δ_obs |
|---|---|---|---|---|
| I, α = 0.5 | least squares, D1(c) | −9.512° | +1.647° | +11.158° |
| I | his chain (C2, which is Brookes's correction) | −9.923° | +1.021° | +10.944° |
| II, α = 0.515 | D1(c) | −7.621° | +4.964° | +12.585° |
| II | his chain | −8.127° | +4.447° | +12.574° |

Step 3c found the shift weakly beyond independent noise and could not say what carries it.
It named the reference he subtracted, Bouvard's 1821 tables with the corrections of his
§ 7, as the first place to look. The printed formula of those tables is now transcribed.

A note on the two figures usually quoted together. "−9.5°" is D1(c) on clean data and "+1°"
is his chain on real data; they are different methods. This plan keeps each method against
its own baseline, as step 3c did.

### 1 Question

Is the shift of about 11° carried by the reference C that Adams subtracted, by the
observations O, or by both?

### 2 Definitions

**His C** (memoir §§ 6–8). The memoir's numbers are "errors of the Tables in heliocentric
longitude" at the time of mean opposition, after "taking into account the several
corrections" of § 7. So C is the heliocentric true longitude of Uranus from Bouvard's 1821
tables plus five terms: two corrections to the tables (the equation in Saturn's mean
longitude minus twice that of Uranus, after Bessel; the principal Jupiter equation, for the
larger mass of Jupiter) and three terms of the second order in the disturbing force, after
Hansen and Delaunay. He says the Greenwich reductions already carried the Jupiter
correction. The coefficients and phases are to be transcribed in V2; the text layer reads
the last three as +32″.00, −8″.35 and −1″.49, not yet confirmed.

**The two references, at his 30 epochs** (the epochs of `adams_memoir_normal_places.csv`,
`t_h1` and `t_h2`; 1690 set aside in every solve, as he did):

- **C_Bouvard:** Bouvard's printed formula for the longitude (route A, groups I to XII and
  the reduction to the ecliptic), converted from the centesimal division, with Bouvard's
  precession from the fixed equinox of 1800, plus the five terms of § 7.
- **C_clean:** step 3b's D3(i) reference unchanged: the model without the perturber, its
  Uranus elements fitted to the clean truth over 1781–1821 (`fit_uranus_elements`).

**The two sets of observations:**

- **O_real** = C_Bouvard + his printed § 8 residuals (observation minus theory, true
  longitude).
- **O_clean:** step 3b's D3 truth unchanged: the model's Uranus with Neptune as perturber.

**What the design reduces to.** Write D(t) = C_Bouvard(t) − C_clean(t), both in one frame.
The solver takes only the column O − C, so the four cells are:

| cell | O | C | column solved |
|---|---|---|---|
| A | real | Bouvard | his § 8 residuals (his printed § 9 column, unchanged) |
| B | real | clean | § 8 residuals + D |
| C | clean | Bouvard | step 3b's clean O − C, minus D |
| D | clean | clean | step 3b's clean O − C (unchanged) |

Three consequences, stated before any run:

1. The only new ingredient is D(t). O_real is never compared with the sky as an absolute
   longitude, so an error in an epoch (Uranus moves 42″ a day) enters only through the
   slow change of D.
2. Cell A does not test C_Bouvard: Bouvard cancels in it by construction. The validity of
   C_Bouvard rests on V1 and V2 alone.
3. An error δ in my evaluation of C_Bouvard moves cells B and C in opposite directions. It
   is counted as "observations" with one sign and "reference" with the other, and it
   leaves their sum unchanged (§ 4).

**Frames and constants per cell.** D is formed in the mean ecliptic and equinox of date:
C_clean's longitude is carried from ecl_j2000 by the precession the earlier steps use
(`lon_j2000_to_of_date`); C_Bouvard is of date by his own precession. A constant or a
linear trend in D is the same as a change of Uranus's epoch and mean motion, which his
method solves for (δε, δn). V3 tests that this holds numerically. In each cell Uranus's ε
and the frame of the reported longitude are those of the cell's C, as in the two gate
cells (his 217° 55′ of date; the C_clean orbit's value in ecl_j2000); his series constants
(e, β) are kept in all four, as D3(i) did. The offset is the heliocentric true longitude
at 1846.762 minus Neptune's in the same frame (327.031° of date; 329.17° ecl_j2000).

**True to mean longitude.** Columns are formed in true longitude and converted by r²/(ab)
with the ellipse of the cell's C. Cell A uses his printed § 9 column as it stands.

### 3 Method

Every cell is solved by both methods of step 3c, with their code unchanged: the
eight-parameter least squares with his § 52 weights (D1(c), 12 θ starts, lowest χ² kept)
and his elimination chain (six roots, the positive-mass root of smaller χ² kept).
Hypothesis I at α = 0.5 is primary; Hypothesis II at α = 0.515 with `t_h2` and its weights
is reported beside it. **D1(c) on Hypothesis I is the method the reading uses**; his chain
is reported beside it and a disagreement between the two is said in the reading's row.

### 4 Attribution

With A, B, C, D the offsets of the four cells:

- effect of the reference, R = ½ [(A − B) + (C − D)];
- effect of the observations, Ob = ½ [(A − C) + (B − D)];
- interaction, I = ½ [(A − B) − (C − D)].

R + Ob = A − D exactly, which is the shift Δ_obs. The gates fix A and D, so the run
measures one new number, C − B = R − Ob, and the interaction.

### 5 Validation, before any cell runs

Each stage has a tolerance. If a stage fails, the work stops there, the failure is
reported with its numbers, and nothing is attributed.

**V1, the formula against the tables.** At five instants I compute Bouvard's longitude of
Uranus two ways: from the transcribed formula, and from his tables as a computer of the
time would, following his instructions (pp. xviii to xx). The table rows these instants
touch are transcribed in two entries, with a reconciliation file; nothing else of the
tables is. The instants are Paris mean midnights, the tables' own epoch, so that the tables
of hours and minutes are not needed: the midnight nearest each of his normal epochs of
1783, 1798, 1813 and 1828, and the midnight nearest 1846.762. By his unit of time
(1810.328 + 3.0362 k years) these fall on or next to 1 January 1783, 8 March 1798,
14 May 1813, 18 July 1828 and 6 October 1846; the exact dates are fixed from
`adams_memoir_normal_places.csv` before the rows are transcribed and recorded with them.

*The pass threshold is the rounding budget of the table path at each instant:* the number
of table entries summed, times half a unit of each entry's last printed digit. Every
longitude table seen so far is printed to 0″.1 centesimal (Tables I, IX, X, XI to XXI), so
half a unit is 0″.05 centesimal, 0″.0162. The entries a longitude takes:

| entry | tables | count |
|---|---|---|
| epoch of the mean longitude | I | 1 |
| motion for other centuries (before 1780 only) | II | 0 or 1 |
| inequality of the precession | Jupiter's Table III, p. 16 | 1 |
| motion for the month (none on 1 January) | III | 0 or 1 |
| motion for the day (none on the 1st) | IV | 0 or 1 |
| great inequality | IX | 1 |
| equation of the centre for 1800 | X | 1 |
| its secular variation | X | 1 |
| equations of the longitude | XI to XXI | 11 |
| reduction to the ecliptic | XXXII | 1 |

A value interpolated between two printed entries counts once: its rounding error is at
most that of one entry.

| instant | entries | budget, centesimal | budget |
|---|---|---|---|
| 1 January 1783 | 17 | 0″.85 | 0″.275 |
| March 1798, May 1813, July 1828, October 1846 (a month and a day each) | 19 | 0″.95 | 0″.308 |
| S4: the 1715 instant (a century, a month and a day) | 20 | 1″.00 | 0″.324 |

If an instant falls on the first of a month, its count drops by one and its budget by
0″.0162. If the precession table or Table XXXII is printed to another last digit than
0″.1, that entry's half unit is used in the sum; the rule does not change.

| largest difference, formula minus tables | consequence |
|---|---|
| within the instant's budget, at all five | pass: C_Bouvard is the formula |
| over the budget at any instant, and ≤ 2″ | stop and ask. The route may then need the table rows of all 30 epochs (route B), which is a new ruling |
| over 2″ at any | fail: the formula as transcribed is not the tables. Stop and report |

Stated beside the result and not part of the threshold: the tables' arguments II to XII
are whole numbers of 1,000 or 10,000 parts of the circle, so each equation also carries an
error of up to its amplitude × 2π/N × ½ from the rounded argument. This allowance is
computed from the transcribed amplitudes with V1 and reported. A difference between the
budget and the budget plus the allowance is still "stop and ask"; the report then says
that it lies within the argument rounding.

Also noted before transcription: the Harvard copy carries a handwritten line at the head
of Table I, "il faut ajouter à la longitude moyenne + 0″,10 t". It is not part of the
print. V1 uses the printed table; the note is transcribed and reported.

The conversion from the centesimal division is tested by itself (`tests/test_centesimal.py`
exists; the tests for table arguments in parts of the circle are added). Which
eccentricity the tables use (p. xv against p. ii) is read off this stage.

**V2, the corrections of § 7.** The five terms are transcribed from the pamphlet's page
image in two entries; the French printing of 1876 may decide a disagreement and the
pamphlet wins. Checks: the two entries agree; each argument is built from the mean
longitudes of Bouvard's p. ii; and, as a test of the true-to-mean conversion that cell C
needs, his § 8 column times r²/(ab) from his ellipse reproduces his § 9 column within
0″.02 (modern) and 0″.1 (ancient), his printing precision. The memoir's definition of C in
§ 2 above is checked against the page images and corrected below the lines if it is wrong.

**V3, the gates and the frame test.**

| check | must give | tolerance |
|---|---|---|
| cell A, D1(c) | +1.647° | 0.001° |
| cell A, his chain (C2) | +1.021° | 0.001° |
| cell D, D1(c), D3(i) | −9.512° | 0.001° |
| cell D, his chain | −9.923° | 0.001° |
| Hypothesis II, the four values of § 0 | as in § 0 | 0.001° |
| frame test: cell D with 10″ + 0″.2 per year added to its column | the offset of cell D | 0.05° |
| frame test: C_clean carried to the equinox of date against the JPL ephemeris's Uranus of date at the 30 epochs, after a constant and a trend are removed | the same curve as in ecl_j2000 | 0″.05 |

The first six are the same code on the same inputs, so they are regression gates: they
show that the pipeline of this step is the pipeline of steps 3b and 3c. They do not test
C_Bouvard (§ 2, consequence 2).

**The propagation check (part of V3; reported, not read).** Let d be the largest V1
difference in absolute value. Cells B and C are each solved again with D replaced by
D + d sin(φ) for four periodic errors: φ = Uranus's mean anomaly and φ = the argument of
the equation in Saturn's and Uranus's mean longitudes (argument III, φ′ − 2φ″), each at
phases 0 and 90°. The change of each cell's offset is reported for each, with the
largest, for both methods. It shows what an error of the size V1 could not exclude does
to a cell, so that the 3° and 6° classes of § 6 can be judged against it. It enters no
reading row. It needs cells B and C, so it is computed with the cells and reported with
V3; the report after V1 to V3 that precedes the cells states d and the four errors to be
applied.

### 6 Reading, fixed in advance

Read on Hypothesis I with D1(c). The shift is +11.16°. An effect is signed: positive is
toward and past Neptune, the direction of the shift.

| effect | class |
|---|---|
| ≥ +6° | carries the shift |
| over +3° and under +6° | partial |
| −3° to +3° | does not carry it |
| under −3° | works against it |

Because R + Ob = +11.16°, the combinations are limited:

| R, the reference | Ob, the observations | reading |
|---|---|---|
| carries | does not | **The reference carries the shift.** With the observations he had, a clean reference would have left his planet about as short as clean data do. His nearness on real data owes to what Bouvard's tables with his corrections contained |
| does not | carries | **The observations carry the shift.** Replacing his reference by a clean one changes little. The cause is in the observed places: their errors, their reductions or the construction of the normal places, which this step cannot tell apart (§ 8) |
| carries | partial | **Mostly the reference**, with a share from the observations, both stated in degrees |
| partial | carries | **Mostly the observations**, with a share from the reference |
| partial | partial | **Both**, in the shares found |
| carries | works against | **The reference carries more than the whole shift and the observations take some back.** Stated with both numbers; no simpler sentence is used |
| works against | carries | the mirror of the row above |

Overriding rows:

- **If the interaction |I| is 3° or more**, the two effects are not additive. The four
  cells are reported, the effect of each factor is given at each level of the other, and
  no single attribution is stated.
- **If D1(c) and his chain fall in different rows**, both rows are given and the reading
  is "depends on the method".
- **Hypothesis II** is reported with the same classes against its own shift (+12.59°). It
  enters no reading; a different class from Hypothesis I is said beside the reading.

Scale, for whoever reads the result: under step 3c's noise model (a), independent noise of
1″.5 on the modern places moved the 1846.762 longitude by several degrees (propagated
σ_θ 4.5°). An effect of 3° is therefore of the size that noise at his residuals' level
produces. The reading describes how the realised columns behave under his method. It is
not a significance statement.

### 7 Sensitivities: named, reported, not read

None of these enters a reading row.

- **S1, Saturn's mass.** The formula's coefficients are printed numbers, so the primary run
  does not use the mass at all, whichever figure is the misprint. If a primary source
  gives Bouvard's mass of Saturn as other than the printed 1/5512, cells B and C are rerun
  with the Saturn-argument terms of the formula scaled by the ratio of the two masses.
  This is what the error would be if the coefficients had been computed with the
  misprinted figure. If no primary source is found, S1 is reported as not run.
- **S2, the eccentricity.** Cells B and C with the equation of the centre rebuilt from
  e″ = 0.0466108 (p. ii) in place of the series of p. xv (about 0.04669).
- **S3, without § 7.** Cells B and C with C_Bouvard as the formula alone. It shows how much
  of R belongs to Adams's five terms and how much to Bouvard's tables.
- **S4, an ancient instant in V1.** The formula against the tables at the Paris midnight
  nearest his 1715 epoch, which uses the century table. Reported with V1 against its own
  budget; it changes no V1 outcome. The five instants of V1 do not test the tables before
  1780.

Expectation, recorded now: S2 should change little, because a change of Uranus's
eccentricity is of the form his δx₁, δy₁ absorb. I have no expectation for R against Ob.
Adams's own argument in § 7 is that terms of long period or of Uranus's period "would
little affect" the planet's position; R measures whether that held for what his C
actually contained.

### 8 What this cannot separate

- **Inside O_real:** the errors of the observations, the reductions made at Greenwich and
  by Adams, his grouping into normal places, and the interpolated place of 1780.
- **Inside C_Bouvard:** Bouvard's theory, Bouvard's fitted elements and Adams's five terms,
  except as S3 splits them.
- **The tables as used.** V1 compares the formula with the tables at five instants. The
  copy Adams or Greenwich used, later errata, and the interpolation habits of a computer
  are not tested. Any such difference is counted as "observations" with one sign and
  "reference" with the other (§ 2, consequence 3).
- **The clean side.** O_clean and C_clean are the model's. Where the model differs from the
  real solar system, the difference is counted as "observations".
- **Method.** The attribution is of his method at his assumed distances. It says nothing
  about a fit that frees the distance.
- **Cause against carrier.** "The reference carries the shift" says where the 11° enters
  his equations, not why Bouvard's tables were what they were.

### 9 Report order, cost, layout

R1 the transcription records (table rows of V1, the terms of § 7), with counts; R2 V1 with
the five differences, the budgets, the argument allowance and S4; R3 V2; R4 V3's gates and
frame tests. **The work stops there and is reported before any cell is run.** Then, only
if V1 to V3 pass and the author says to go on, R5 D(t) at the 29
places (table and figure, with the constant and trend removed shown beside); R6 the four
cells for both methods and both hypotheses; R7 R, Ob, I and the reading's row; R8 the
sensitivities; R9 what remains unseparated. Results are appended below the lines;
`docs/notes/P1B_bouvard_attribution.{json,png}`; the site and the earlier notes untouched.

Cost: each solve is tens of milliseconds; the whole run is under a minute, not detached.
The work is in the transcriptions: about 330 cells for V1's table rows and about 20 for
§ 7, each in two entries.

Layout: `src/fork/historical/bouvard_1821.py` (the formula, the table lookup for V1, the
five terms), the runner `src/fork/export/p1b_bouvard.py` (`validate`, `cells`,
`summarize`), tests in `tests/test_p1b_bouvard.py` (the entries agree; V1 to V3 as tests on
the committed numbers; the algebra R + Ob = A − D; the reading rows on constructed
effects). Order of work: this plan committed and pushed; the transcriptions; V1; V2; V3;
the report and the author's word; the cells with the propagation check; the results. One commit for the plan, one for the transcriptions and code, one
for the results; no trailers.

### 10 Rulings made on approval (2026-11-13, before the commit; nothing runs before it)

1. V1's pass threshold is the rounding budget of the table path at each instant, stated
   in § 5 before any table row is transcribed; over it and up to 2″ is "stop and ask"; over
   2″ fails.
2. The propagation check is added to V3: reported, not read.
3. Least squares on Hypothesis I is the method the reading uses, with his chain beside it.
4. The class "works against it" and the interaction override at 3° stand.
5. S3 and S4 are part of the plan. S1 runs only on a primary source for Saturn's mass.
6. The work stops after V1 to V3 and is reported before any cell runs.
7. Made with ruling 1, to keep the entry count small and certain: V1's instants are the
   Paris mean midnights nearest the approved instants, not the instants themselves. V1
   compares the formula with the tables at one and the same instant, so its purpose is
   unchanged.

---

## Results, 13 November 2026: V1 fails; the work stops before V2, V3 and every cell

Nothing above the line is changed. No cell has been run, and no gate of V3.

### Transcription

Two entries of every cell, the first from bands cut from the page images, the second made
independently from the same bands and a column specification alone.

| what | cells | differences between the entries |
|---|---|---|
| table rows of Uranus the six instants touch (Tables I to IV, IX, X, XI to XXI, XXXII) and four rows of the precession table | 764 | 0 |
| the great inequalities A, A′, A″ (pp. ii, iii) | 24 | 0 |
| Adams's § 7, five terms | 15 | 0 |

Instants (Paris civil date, the midnight that begins it): 1 January 1783, 7 March 1798,
13 May 1813, 18 July 1828, 5 October 1846; S4: 4 March 1715.

Deviations from the plan, all made before any comparison was read except the third:

1. The great inequalities were transcribed as well. The formula of p. xv is written in φ,
   φ′, φ″, which carry A, A′, A″ of pp. ii and iii; the 104 cells of the first
   transcription do not hold them, and without them the formula is not the theory.
2. The precession table is Jupiter's Table III on p. 12, not p. 16 as § 5 says.
3. Two readings of Table X are replaced, because the table's own printed first differences
   on both sides contradict them: row 324 reads 4.29.83,8 and is used as 4.28.83,8; row
   327 reads 4.39.52,6 and is used as 4.38.52,6. Both entries read the same ink. Only the
   1783 instant touches these rows. With the readings as they stand its difference is
   −0″.13; with the replacements −2″.93. The outcome below does not depend on it. One more
   reading the differences contradict (Table XIII, row 720, 2,6 against 2,5) is left as read.
4. Interpolation: four points in Tables IX, X, XII and the precession table, two elsewhere.

### V1

Formula minus tables, heliocentric longitude of Uranus.

| instant | budget | argument allowance | formula − tables | outcome |
|---|---|---|---|---|
| 1783 | 0″.275 | 0″.60 | −2″.93 | fail |
| 1798 | 0″.308 | 1″.04 | +0″.26 | pass |
| 1813 | 0″.308 | 0″.40 | +3″.40 | fail |
| 1828 | 0″.308 | 0″.80 | +2″.88 | fail |
| 1846 | 0″.308 | 0″.97 | −0″.38 | stop and ask |
| S4, 1715 | 0″.324 | 0″.98 | +1″.99 | (not in the threshold) |

**V1 fails: the largest difference is 3″.40, over 2″. The formula as transcribed is not the
tables.** By § 5 the work stops here. d for the propagation check would be 3″.40.

Where the difference sits (centesimal seconds, formula minus tables; × 0.324 for ″):

| instant | mean longitude | precession inequality | great inequality | equation of the centre | perturbations | reduction |
|---|---|---|---|---|---|---|
| 1783 | +0.23 | −0.32 | +0.03 | −8.91 | −0.04 | −0.03 |
| 1798 | +1.91 | −0.74 | +0.07 | −0.26 | −0.20 | +0.03 |
| 1813 | +3.46 | −1.33 | −0.04 | +8.66 | −0.29 | +0.02 |
| 1828 | +5.11 | −2.10 | −0.09 | +6.99 | −1.03 | +0.00 |
| 1846 | +7.06 | −3.25 | −0.15 | −4.59 | −0.26 | +0.01 |
| 1715 | +3.66 | −0.56 | −0.05 | +0.79 | +2.28 | +0.03 |

Three causes, read off the parts:

- **The equation of the centre.** Table X is the series of p. xv with a leading coefficient
  of 59417″.54, ten seconds under the printed 59427″.54. With 59417″.54 the equation of the
  centre agrees at every instant, and the eccentricity it implies is 0.046679, the value
  Adams's memoir uses (the printed coefficient gives 0.046687; p. ii prints 0.0466108).
  So the tables use the Uranus solution's eccentricity, and p. xv's first coefficient
  looks like a misprint of one digit. This is inferred from six instants, not shown.
- **The mean motion.** Table I's mean longitude falls behind the elements of p. ii by
  0″.105 centesimal a year from 1780 (0.076, 0.105, 0.104, 0.105, 0.106 at the five
  instants). The handwritten line in the Harvard copy at Table I says to add 0″,10 t to
  the mean longitude, t in years from 1780: the same amount.
- **The precession.** The tables apply the inequality of the precession from Jupiter's
  Table III; the formula has the uniform precession only.

Post hoc, and not a pass of V1: with 59417″.54 in the formula, the differences are +0″.05,
+0″.36, +0″.50, +0″.34, +1″.06 (1715: +1″.85). Taking out the mean-longitude and precession
differences as well leaves +0″.08, −0″.02, −0″.19, −0″.64, −0″.18 (1715: +0″.84): inside
the budget at four instants, and inside the budget plus the argument allowance at all
five. A constant and a trend fitted to the five differences with 59417″.54 leave −0″.02,
+0″.10, +0″.05, −0″.31, +0″.18. A constant, a trend, and a change of eccentricity are the
forms his δε, δn, δx₁, δy₁ take up, but that the differences reach no cell's offset is
not tested here: no cell has been run.

Also read off this stage:

- The worn sign before group III is −. With + the largest difference is 270″.
- The anomaly is the longitude of date less the perihelion of date.
- Table IX is A″ of p. iii (agreement within 0″.2 centesimal at all six instants).
- The sum of the constants added to the equations of the longitude is 717″.5 centesimal.
- 1715 is not explained to the same level: +2.3 centesimal seconds remain in the
  perturbations, with Table II's arguments for −100 years in the path.

Second-entry flags kept with the rows: three damaged digits in Tables I and IV, read alike
by both entries.

### Not done

V2's checks (the § 7 terms are transcribed; the arguments, the § 8 to § 9 conversion and
the definition of C are not checked), V3, the cells, the sensitivities. A ruling is needed
on what C_Bouvard is before anything further: the formula as printed has failed; the
alternatives are the tables themselves at the 30 epochs (route B), or the formula with the
three differences above built in and V1 run again on instants not used to find them.

Files: `data/historical/adams_challis/bouvard_1821/` (rows, passes),
`adams_memoir_s7_terms.csv`, `src/fork/historical/bouvard_1821.py`,
`src/fork/export/p1b_bouvard.py`, `tests/test_p1b_bouvard.py`,
`docs/notes/P1B_bouvard_attribution.json`.

## Amendment, 13 November 2026: what C_Bouvard is after the failure of V1, and V1′

Ruled by the author on 13 November 2026, after the results above were pushed (235bab4).
Written and pushed before any new row is transcribed. It supersedes the first row of the
consequence table of V1 in § 5 ("pass: C_Bouvard is the formula") and nothing else; no
earlier text is changed.

**V1 failed by its rule.** The formula as printed is not C_Bouvard.

**What C_Bouvard is to represent.** The tables. Adams's "errors of the Tables" were taken
against them, not against the formula of the introduction.

**The corrected formula** stands in for the tables only if V1′ below passes. It is the
formula of pp. ii, iii and xv as transcribed, with three changes, each fixed here:

1. The leading coefficient of the equation of the centre is 59417″.54, not the printed
   59427″.54.
2. The mean longitude of Uranus that enters the longitude and the anomaly is reduced by
   0″.105 centesimal × (Y − 1780), Y the date in years; for a date before 1780, which the
   tables reach by the row of Table I a century later and Table II, Y is the date plus
   100 years. The rate is the one V1 measured between Table I and the elements of p. ii.
   The handwritten line at the head of Table I in the Harvard copy ("il faut ajouter à la
   longitude moyenne + 0″,10 t, t étant le nombre d'années écoulées depuis 1780") names
   the same difference; it is cited as a witness, and the printed table is what is
   represented. The arguments II to XII are not changed: the tables carry them in their
   own columns.
3. The inequality of the precession is added to the longitude, the perihelion and the
   node as the tables apply it: interpolated in the rows of Jupiter's Table III (p. 12),
   by the same rule in the formula and in the table path. This part therefore cancels in
   V1′ by construction and is not tested by it.

Not known and not assumed: whether the copy of the tables that Greenwich or Adams used
carried the correction of the handwritten line.

**V1′.** The corrected formula against the tables at seven instants not used to find the
corrections, each the Paris mean midnight nearest one of his normal epochs:

| normal place | instant (Paris civil date, the midnight that begins it) | entries | budget |
|---|---|---|---|
| 1789 | 26 January 1789 | 18 | 0″.292 |
| 1804 | 3 April 1804 | 19 | 0″.308 |
| 1822 | 22 June 1822 | 19 | 0″.308 |
| 1837 | 27 August 1837 | 19 | 0″.308 |
| 1690 | 20 December 1690 | 20 | 0″.324 |
| 1750 | 5 November 1750 | 20 | 0″.324 |
| 1769 | 12 January 1769 | 19 | 0″.308 |

The table rows these touch are transcribed in two entries, as for V1, with the rows of the
precession table that bracket 1690. The table path, the interpolation and the count of
entries are those of V1 (`fork.historical.bouvard_1821.table_longitude`), unchanged.

Rule, as I read "the same budget-plus-allowance rule": an instant passes when the absolute
difference, corrected formula minus tables, is at most its rounding budget plus its
argument allowance, both computed as in V1; otherwise it fails. The difference against
the budget alone is reported beside.

| outcome | consequence |
|---|---|
| all seven pass | the corrected formula is C_Bouvard |
| any modern instant fails | route B: the table rows of all 30 epochs |
| the four modern pass, any ancient fails | route B for the ancient epochs only; the corrected formula for the modern |

A reading that the table's printed first differences on both sides contradict is replaced
only when both differences give the same value; the difference is reported both ways. If
the two ways give different outcomes at an instant, that instant is "stop and ask".

Reported with V1′ and not part of its rule: the corrected formula at the six instants of
V1 (found with them, so no test); the 1715 instant, whose +0″.84 after the three
corrections is unexplained and stays so unless V1′ explains it.

A correction of fact to the results above: the 1783 instant of V1, 1 January 1783, is the
epoch of Table I's row, 0.95 day after his normal epoch; the nearest midnight is that of
31 December 1782. V1 compares both paths at one instant, so nothing in it changes.

After V1′ the work stops and is reported before V2.

## Results, 13 November 2026: V1′; one instant is "stop and ask"; stopped before V2

Run after the amendment above was pushed (7306594). Nothing above is changed. No cell
has been run, and neither V2 nor V3.

### Transcription

609 cells of table rows for the seven instants, two entries (the first from bands cut from
the page images, the second made independently from the same bands and a column
specification alone). 0 differences. The row of Jupiter's Table III above 1700 is the row
of 1600; for the 1690 instant the precession inequality is interpolated in the rows of
1600, 1700, 1750 and 1800, for every other instant as in V1.

Readings the neighbouring figures contradict, read alike by both entries:

| table, row | reads | the printed differences on both sides give | used | instant it touches |
|---|---|---|---|---|
| X, 25 | 2.35.14,3 | 2.33.14,3 | 2.33.14,3 | 1804 |
| XII, 2600 | 3.41,6 | 3.31,6 | 3.31,6 | 1789 |
| XXXII, 80 | 17,7 | no differences printed; the neighbours (18,5; 17,8; 16,3) suggest 17,1 | 17,7 | 1750 |

### V1′

Corrected formula minus tables.

| instant | budget | budget + allowance | difference | with the readings as they stand | outcome |
|---|---|---|---|---|---|
| 26 January 1789 | 0″.292 | 1″.09 | −0″.50 | −2″.40 | **stop and ask** |
| 3 April 1804 | 0″.308 | 1″.33 | +0″.00 | +0″.72 | pass |
| 22 June 1822 | 0″.308 | 1″.18 | +0″.02 | +0″.02 | pass |
| 27 August 1837 | 0″.308 | 1″.32 | −0″.12 | −0″.12 | pass |
| 20 December 1690 | 0″.324 | 1″.60 | −1″.05 | −1″.05 | pass |
| 5 November 1750 | 0″.324 | 1″.20 | +0″.65 | +0″.65 | pass |
| 12 January 1769 | 0″.308 | 1″.17 | −0″.48 | −0″.48 | pass |

**Outcome by the rule of the amendment: stop and ask.** Six instants pass. At 1789 the
corrected formula passes with Table XII's row 2600 as its printed differences give it and
fails with the row as it reads, so the rule makes that instant "stop and ask" and no
consequence row applies until it is ruled. The formula as printed is 1″.6 to 4″.3 from
the tables at six of the seven instants (−2″.45, +1″.90, +4″.04, +1″.58, −4″.34, +4″.22,
−0″.25).

Against the budget alone: 1804, 1822 and 1837 are inside it; 1789, 1690, 1750 and 1769
are over it and inside the budget plus the allowance.

Where what is left sits (centesimal seconds, corrected formula minus tables):

| instant | mean longitude | great inequality | equation of the centre | perturbations | reduction |
|---|---|---|---|---|---|
| 1789 | −0.06 | −0.02 | +0.11 | −1.58 | −0.01 |
| 1804 | +0.00 | +0.06 | −0.03 | +0.01 | −0.03 |
| 1822 | +0.01 | −0.08 | −0.32 | +0.48 | −0.02 |
| 1837 | −0.03 | −0.13 | −0.01 | −0.18 | −0.03 |
| 1690 | −0.06 | −0.07 | +0.09 | −3.17 | −0.03 |
| 1750 | +0.02 | −0.05 | +0.04 | +1.71 | +0.30 |
| 1769 | +0.03 | +0.02 | +0.24 | −1.76 | −0.00 |

- The two corrections that V1′ can test hold at all seven fresh instants: the mean
  longitude agrees within 0.06 centesimal second (the rate of 0″.105 a year, and the
  rule of the row a century later for the ancient dates), and the equation of the centre
  within 0.32 (the coefficient 59417″.54, at seven new anomalies).
- What is left is in the equations of the perturbations, and it is larger at the ancient
  instants. There the tables' arguments, carried back by Table II, differ from the
  formula's by more than their rounding in places (argument III by 8 parts of 10,000 at
  1769, 2 at 1690 and 1750).
- The +0.30 in the reduction at 1750 is the reading 17,7 of Table XXXII, row 80. With
  17,1 the 1750 difference is +0″.56.

Not part of the rule: the corrected formula at the instants of V1 gives +0″.05, −0″.06,
−0″.23, −0″.57, −0″.04, and at 1715 +0″.73 (+0″.84 before the precession was interpolated
by the tables' rule). **1715 stays unexplained** in cause: V1′ shows the same kind of
difference at the other three ancient instants (in the perturbations, through the
arguments), within the budget plus the allowance, but does not say where it comes from.

A ruling is needed on the 1789 instant before any consequence is drawn, and on whether
the ancient instants are close enough.

## Amendment, 13 November 2026: rulings on V1′; S5; d for the propagation check

Ruled by the author on 13 November 2026, after the V1′ results above were pushed
(821b9c0). Written and pushed before V2 or V3 is run. No earlier text is changed.

1. **Table XII, row 2600, is 3.31,6.** Both printed difference columns agree on it. It
   and Table X, row 25 (2.33.14,3 by the differences) are logged as misprints of the
   tables in `data/historical/adams_challis/bouvard_1821/bouvard_1821_misprints.csv`.
   Table XXXII, row 80, stays as read (17,7), with 17,1 noted as the alternative.
2. **V1′ passes at all seven instants** under the rule of the amendment above (budget
   plus allowance). **The corrected formula is C_Bouvard for all 30 epochs**, the ancient
   ones included. This replaces "the printed formula (route A …)" in the definition of
   C_Bouvard in § 2; the five terms of § 7 are added as before.
3. **S5, a new named sensitivity, reported and not read.** Cells B and C are solved once
   more with C_Bouvard at the 1789 normal place changed to what a computer using row
   2600 as printed would have got: the equation of argument III interpolated by the
   tables' rule (four rows) with 3.41,6 in place of 3.31,6, at the argument of that
   epoch. At the V1′ instant of 1789 this is +1″.90 in the tables' longitude. The change
   of each cell's offset is reported for both methods. It bounds the question of the
   tables as used at that epoch. Any other normal place whose argument III falls on the
   same four rows is changed in the same way and listed.
4. **1715 stays unexplained** (+0″.73). The differences at the ancient instants (−1″.05,
   +0″.65, −0″.48, and +0″.73 at 1715) are small against the residuals his solution
   leaves at the ancient places, which reach ±12″.8.
5. **The propagation check of V3 uses d = 1″.05**, the largest V1′ difference under the
   adopted readings, in place of the largest V1 difference, at the two periods and two
   phases of § 5. S5 is separate from it.

Next: V2 and V3; then the work stops and is reported before any cell.

## Results, 13 November 2026: V2 and V3 pass; stopped before any cell

Run after the amendment above was pushed (ed463f4). No cell has been run: nothing has
been solved with D in its column, and no sensitivity.

### V2, the corrections of § 7

| check | result | tolerance | |
|---|---|---|---|
| two entries of the five terms | 15 cells, 0 differences | | pass |
| § 8 × r²/(ab) against § 9, modern places | largest 0″.009 | 0″.02 | pass |
| § 8 × r²/(ab) against § 9, ancient places | largest 0″.051 | 0″.1 | pass |

The factor uses his ellipse: e = 0.046679, ε − ϖ = 50° 15′.8 at 1810.328, n = 13°.01 a
unit. So the true-to-mean conversion that cell C needs is his.

The definition of C in § 2 was checked against the memoir's §§ 6 to 8 and stands: the
"errors of the Tables in heliocentric longitude" at the mean oppositions, after "the
corrections to be applied to Bouvard's Tables" (two terms) and "the terms to be added to
the longitude" (three); the Greenwich reductions already carried the Jupiter term.

The five terms at his 30 epochs, arguments from the mean longitudes of p. ii:

| term | coefficient | argument | largest after a constant and a trend are removed |
|---|---|---|---|
| 1 | +1″.918 | φ₁ − 2φ₂ − 13° 1′.5 | 2″.25 |
| 2 | +1″.085 | φ − φ₂ | 1″.07 |
| 3 | +32″.00 | 3φ₂ − 6φ₁ + 2φ + 22° 18′.8 | 0″.82 |
| 4 | −8″.35 | 2φ₂ − 6φ₁ + 2φ + 39° 10′.5 | 8″.56 |
| 5 | −1″.49 | 4φ₂ − 6φ₁ + 2φ + 34° 48′.4 | 1″.82 |

The sum runs from −29″.9 (1715) to −9″.9 (1780); with a constant and a trend removed it
is within ±10″.0. Term 3 is close to a constant and a trend over his 150 years, as he
says of the terms of long period; term 4 has Uranus's period and carries most of what
is left.

Two things the memoir does not settle, measured and not ruled:

- "Mean longitudes" with or without the great inequalities of Jupiter and Saturn. With
  them the sum differs by 1″.45 to 3″.66, and by up to 1″.09 after a constant and a trend
  are removed. The primary is without, as § 5 says (p. ii).
- Counted from the fixed equinox of 1800 (p. ii, the primary) or from the equinox of
  date: up to 0″.23 after a constant and a trend are removed.

The first is the same size as d of the propagation check.

### V3, the gates and the frame tests

| check | must give | gives | |
|---|---|---|---|
| cell A, D1(c), Hypothesis I | +1.647° | +1.6469° | pass |
| cell A, his chain | +1.021° | +1.0211° | pass |
| cell D, D1(c) | −9.512° | −9.5115° | pass |
| cell D, his chain | −9.923° | −9.9230° | pass |
| Hypothesis II, cell A, D1(c) and chain | +4.964°, +4.447° | +4.9641°, +4.4470° | pass |
| Hypothesis II, cell D, D1(c) and chain | −7.621°, −8.127° | −7.6206°, −8.1267° | pass |
| cell D with 10″ + 0″.2 a year added to its column, D1(c) | no change, 0.05° | +0.00002° | pass |
| the same, his chain | no change, 0.05° | 0.00000° | pass |
| ephemeris Uranus minus C_clean at the 30 epochs, of date against ecl_j2000, constant and trend removed (the curve runs −58″ to +45″) | 0″.05 | 0″.0009 | pass |

The pipeline of this step is the pipeline of steps 3b and 3c, and a constant or a trend
in D reaches neither method's offset.

**The propagation check**, to be computed with cells B and C: d = 1″.05; the four errors
are d sin M and d cos M (M Uranus's mean anomaly) and d sin III and d cos III (III the
argument φ′ − 2φ″), added to D.

**S5**, to be computed with the cells: the printed row 2600 changes the tables' longitude
by +1″.91 at the 1789 normal place and by −0″.02 at 1690, and at no other place.

Seen while preparing S5 and not ruled: Table X's row 25 as printed (2.35.14,3) is 200
centesimal seconds from the value its differences give. The 1804 normal place has an
anomaly of 26.9, where the four-row interpolation gives it a weight of +0″.72; a computer
interpolating between rows 26 and 27 alone would not have touched it.

Next, on approval: the cells.

## Amendment, 13 November 2026: the definition of C checked on the page images; S6, S7, S5b

Ruled by the author on 13 November 2026, after the V2 and V3 results above were pushed
(393dc5f). Written and pushed before any cell is run. No earlier text is changed.

1. **The definition of C, on the page images.** V2 above checked § 2's definition against
   the pamphlet's text layer; § 5 asks for the page images. Pages 5 and 6 of the pamphlet
   (§§ 6 to 8) were read as images and agree with the text layer and with § 2 in every
   point: "the error of the Tables in heliocentric longitude … corresponding to the time
   of mean opposition in the middle year of the group"; 21 normal errors between 1780
   and 1840; "the corrections to be applied to Bouvard's Tables" for the equation in
   Saturn's mean longitude minus twice that of Uranus (Bessel) and the principal
   equation of Jupiter (its mass); "In the Reduction of the Greenwich Observations, the
   latter correction was already taken into account"; three terms "depending on the
   square of the disturbing force" (Hansen, Delaunay), "to be added to the longitude";
   § 8: "Taking into account the several corrections above-mentioned, the residual
   differences between the theoretical and observed heliocentric longitudes". The five
   terms on the image are those transcribed. Nothing in § 2 is corrected.
2. **Three more named sensitivities, reported and not read.** Each solves cells B and C
   once more, by both methods, and reports the change of each offset.
   - **S6:** the mean longitudes in the arguments of § 7 include the great inequalities
     A, A′, A″ of pp. ii and iii.
   - **S7:** those mean longitudes are counted from the equinox of date (Bouvard's
     154″.63 a year added to each) instead of the fixed equinox of 1800.
   - **S5b:** C_Bouvard at the 1804 normal place as Table X's row 25 as printed
     (2.35.14,3) gives it under the four-row interpolation, at the anomaly of that
     epoch. At the V1′ instant of 1804 the difference, formula minus tables, was +0″.72
     with the row as printed, so the tables' longitude changes by −0″.72. Any other
     normal place whose anomaly falls on the same four rows is changed likewise and
     listed.
3. **The run.** Cells A to D, both methods, both hypotheses, with S2, S3, S5, S5b, S6, S7
   and the propagation check (d = 1″.05). S1 is not run: no primary source for Bouvard's
   mass of Saturn has been found. S4 was reported with V1.

## Results, 13 November 2026: the cells

Run after the amendment above was pushed (98e08c7). Runner:
`uv run python -m fork.export.p1b_bouvard` (about six minutes, not under a minute as § 9
expected: 4 cells and 12 variants, each by both methods and both hypotheses). Numbers:
`P1B_bouvard_attribution.json`; figure: `P1B_bouvard_attribution.png`.

**In one sentence:** by the rule fixed in advance the reading is "both" (reference +5.3°,
observations +5.9° of the +11.16° shift), but it depends on the method, on the
hypothesis, and on a few arcseconds at single places near 1780 to 1790, so no
attribution is supported.

### R5, D = C_Bouvard − C_clean at the 30 places

A constant of −14″.5 (at 1800) and a trend of −0″.81 a year, fitted to the 29 places
used. What is left runs from −10″.5 (1792) to +10″.5 (1750), a wave of about Uranus's
period (figure).

| year | D | less constant and trend | year | D | less constant and trend |
|---|---|---|---|---|---|
| 1690 | +75″.9 | +1″.9 | 1801 | −23″.2 | −7″.7 |
| 1712 | +57″.7 | +1″.0 | 1804 | −22″.5 | −4″.6 |
| 1715 | +58″.1 | +3″.8 | 1807 | −21″.3 | −0″.9 |
| 1750 | +35″.9 | +10″.5 | 1810 | −21″.5 | +1″.3 |
| 1753 | +32″.1 | +9″.2 | 1813 | −21″.6 | +3″.6 |
| 1756 | +28″.7 | +8″.0 | 1816 | −22″.4 | +5″.3 |
| 1764 | +18″.7 | +4″.0 | 1819 | −23″.0 | +7″.2 |
| 1769 | +10″.3 | −0″.3 | 1822 | −25″.8 | +6″.9 |
| 1771 | +5″.0 | −3″.3 | 1825 | −29″.1 | +6″.1 |
| 1780 | −7″.3 | −9″.1 | 1828 | −32″.5 | +5″.1 |
| 1783 | −7″.8 | −7″.2 | 1831 | −35″.5 | +4″.5 |
| 1786 | −12″.0 | −8″.9 | 1834 | −39″.1 | +3″.4 |
| 1789 | −16″.0 | −10″.4 | 1837 | −43″.0 | +2″.0 |
| 1792 | −18″.5 | −10″.5 | 1840 | −47″.4 | +0″.1 |
| 1795 | −20″.6 | −10″.1 | | | |
| 1798 | −22″.0 | −9″.0 | | | |

The 1780 place stands 2″ off its neighbours. Its epoch is 1779.96, so the rule of the
first amendment gives it the row of Table I a century later, whose mean longitude is
3″.4 from Table I's own branch. See the post hoc check below.

### R6, the four cells (offset from Neptune at 1846.762, degrees)

| cell | O | C | Hyp. I, D1(c) | Hyp. I, chain | Hyp. II, D1(c) | Hyp. II, chain |
|---|---|---|---|---|---|---|
| A | real | Bouvard | +1.647 | +1.021 | +4.964 | +4.447 |
| B | real | clean | −3.170 | −10.914 | −1.418 | −10.326 |
| C | clean | Bouvard | −3.801 | −1.187 | −0.714 | +3.722 |
| D | clean | clean | −9.512 | −9.923 | −7.621 | −8.127 |

Cells A and D are the gates. The chain found two positive-mass roots in every cell. In
cell B the two methods are 7.7° apart (χ² 22.5 for D1(c), 24.4 for the chain's root under
the same weights); in cell C 2.6° apart. D1(c)'s residuals: A 1″.23 modern and 6″.15
ancient; B 1″.08 and 5″.68; C 1″.00 and 2″.08; D 0″.43 and 1″.84.

### R7, attribution and reading

| | shift | R, the reference | Ob, the observations | I | classes |
|---|---|---|---|---|---|
| **Hyp. I, D1(c) (the reading)** | +11.16 | **+5.26** | **+5.89** | −0.45 | partial, partial |
| Hyp. I, chain | +10.94 | +10.34 | +0.61 | +1.60 | carries, does not carry |
| Hyp. II, D1(c) | +12.59 | +6.65 | +5.94 | −0.26 | carries, partial |
| Hyp. II, chain | +12.57 | +13.31 | −0.74 | +1.46 | carries, does not carry |

By § 6:

- The row for D1(c) on Hypothesis I is **"Both, in the shares found"**: +5.3° from the
  reference and +5.9° from the observations.
- The interaction is under 3° everywhere, so that override does not apply.
- **D1(c) and his chain fall in different rows.** The chain's row is "The reference
  carries the shift". By § 6 both rows are given and **the reading is "depends on the
  method"**.
- Hypothesis II with D1(c) falls in a different class for the reference (+6.65°,
  "carries"; the row "Mostly the reference").

Both partial effects of D1(c) lie within 1° of a class boundary (3° and 6°).

### R8, sensitivities: reported, not read

Change of the offset of cells B and C, degrees, Hypothesis I.

| | δ in C_Bouvard | B, D1(c) | C, D1(c) | B, chain | C, chain |
|---|---|---|---|---|---|
| S2, eccentricity of p. ii | ±28″ at Uranus's period | −0.09 | +0.04 | −0.16 | +0.03 |
| S3, without § 7 | +9″.9 to +29″.9 | −5.42 | +1.66 | −7.74 | −1.10 |
| S5, Table XII row 2600 as printed | +1″.91 at 1789 | −2.07 | +1.15 | +171.94 | +0.55 |
| S5b, Table X row 25 as printed | −0″.73 at 1804 | −0.68 | +0.36 | −1.05 | +0.18 |
| S6, § 7 with the great inequalities | +1″.4 to +3″.7 | −0.16 | +0.12 | −0.15 | +0.07 |
| S7, § 7 from the equinox of date | −0″.1 to +0″.7 | +0.10 | −0.08 | +0.09 | −0.08 |
| propagation, d sin M | ±1″.05 | −0.32 | +0.20 | −0.32 | +0.12 |
| propagation, d cos M | ±1″.05 | −0.11 | +0.05 | −0.14 | −0.02 |
| propagation, d sin III | ±1″.05 | −0.01 | +0.08 | +0.01 | +0.12 |
| propagation, d cos III | ±1″.05 | −0.07 | +0.12 | +0.00 | +0.18 |

S1 was not run: no primary source for Bouvard's mass of Saturn. S4 is with V1.

- **S2** changes little, as expected in § 7: a change of eccentricity is taken up by his
  δx₁, δy₁.
- **The propagation check:** a periodic error of the size V1′ could not exclude moves a
  cell by at most 0.32°. Smooth errors of that size do not threaten the 3° and 6° classes.
- **S3:** his five terms are worth several degrees. Without them R is +8.8° and Ob +2.4°
  (D1(c)).
- **S5:** one place changed by 1″.9 moves cell B by 2.1° and cell C by 1.2°; R becomes
  +6.9° ("carries") and Ob +4.3°. In the chain the same change makes cell B take its
  other root (offset +161°).
- **S5b, S6, S7** move a cell by under 0.7° (D1(c)).

### Post hoc, not in the plan

1. **The 1780 place.** Adams says its error was "found by interpolating between the errors
   of 1781, 1782, and 1783, and those given by the Ancient Observations of 1769 and
   1771". It was not computed from the tables at its own date, and the two table branches
   differ there by 3″.4 in mean longitude. With C_Bouvard at 1780 taken on Table I's own
   branch (+3″.40 at that one place), cell B moves by +4.81° and cell C by −2.76°
   (D1(c)); by +13.24° and −6.49° in the chain. The attribution becomes:

   | | R | Ob | I |
   |---|---|---|---|
   | Hyp. I, D1(c) | +1.48 | +9.68 | −1.47 |
   | Hyp. I, chain | +0.47 | +10.47 | −1.77 |
   | Hyp. II, D1(c) | +1.98 | +10.61 | −1.31 |
   | Hyp. II, chain | +0.60 | +11.97 | −1.84 |

   That is "the observations carry the shift" by both methods and both hypotheses, the
   opposite of what the chain gives under the rule as fixed. Which branch stands for an
   interpolated place the memoir does not say.
2. **A frame rate.** A trend of 0″.1 a year in true longitude (a difference of
   precession constants is of this size) moves cell B by +0.21° and cell C by −0.22°.
   Frame test 1 of V3 added the trend to the mean-longitude column, where it is taken up
   exactly; in true longitude a small part is not.

### R9, what this leaves

- The pre-registered reading stands as computed: "both", with "depends on the method"
  and a different class on Hypothesis II.
- It is not robust. Changes of 2″ to 3″ at one normal place (1789 in S5, 1780 post hoc)
  move a cell by 2° to 5° under D1(c) and switch roots or rows in the chain. The classes
  of § 6 are 3° wide. Step 3c found the same sensitivity to noise of this size.
- So this step does not attribute the shift. What it does show: smooth, periodic errors
  in the reference of the size the validation left open do not matter (under 0.4°); the
  five terms of § 7 matter (several degrees); and the result turns on single places
  between 1780 and 1790.
- Unseparated, as in § 8: the observations' errors, their reductions and the
  construction of the normal places. To this is added the treatment of the 1780 place.
