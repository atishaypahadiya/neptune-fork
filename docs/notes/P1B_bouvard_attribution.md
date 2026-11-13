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
