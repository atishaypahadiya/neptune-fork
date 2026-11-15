# Step 3e: the four cells of step 3d with Bouvard's tables themselves (route B). Pre-registration

Written 13 November 2026. Committed and pushed before any table longitude is computed at the
epochs below; the rows it uses were entered twice and committed first, and their hashes are
in § 1. No longitude has been formed from them.
Results are appended below the line at the end; nothing above it is edited afterwards.

## 0 Scope and question

Step 3d (`docs/notes/P1B_bouvard_attribution.md`) attributed nothing: its registered
reading was "depends on the method", and it turned on a few arcseconds at single places.
Its reference C_Bouvard was the corrected formula, which stands in for the tables within
the rounding budget plus the argument allowance at the instants tested (V1′), not at the
30 epochs themselves. Adams's "errors of the Tables" were taken against the tables.

**Question: with C_Bouvard computed from the printed tables at his epochs, as a computer
of 1821 to 1846 would, (R1) does it agree with the corrected formula within the tables'
rounding, (R2) what do the four cells of step 3d give, and (R3) is the reading the same?**

Out of scope: any change to the cells' method, to O_clean and C_clean, to the five terms of
§ 7, or to the classes of step 3d's § 6. The radius vector, the latitude and the
aberration are not used.

## 1 Inputs, fixed

| input | file | SHA-256 |
|---|---|---|
| the scan (Harvard copy, Google Books id `SfZ4Eiw9yR0C`; cached, not committed) | `data/cache/adams_challis/bouvard/bouvard_1821_SfZ4Eiw9yR0C.pdf` | `ebd3d9fbe1d44a9cc23eeb23225e15f44b226ca29b20119df0bd0c4ea4f6d1f6` |
| table rows already entered twice and reconciled (1,373 cells) | `data/historical/adams_challis/bouvard_1821/bouvard_1821_uranus_table_rows.csv` | `5f8919c5477e3bf959828697f3a0db876bf0ad376a44b2c7a512a6f55e4272c5` |
| readings contradicted by the tables' own figures, with the rulings | `…/bouvard_1821/bouvard_1821_misprints.csv` | `459eb0914691df29e9acaae265b2d85db94413ee91da0c196424fec0ae9c0a89` |
| the printed formula; the great inequalities | `…/bouvard_1821_uranus_formula.csv`; `…/bouvard_1821_great_inequalities.csv` | `42ceb111a5d0aece67cd43a7d57e5f0a46aaac889c15aa7b35a160ef2c7be1d3`; `01b494d839ac2fa6079da5d18d3a08fc59302ebfdad8c723918e3c180747f25b` |
| Adams's 30 normal places (epochs `t_h1`, `t_h2`; §§ 8, 9 columns) | `data/historical/adams_challis/adams_memoir_normal_places.csv` | `52aa0a9cb1459ab7181991e7441f22c303e87a8fbf4166564a8e6d1481f12963` |
| his five terms of § 7 | `data/historical/adams_challis/adams_memoir_s7_terms.csv` | `53339c28a8e9f9ce223fb5d92a2699f94720d36ca14298eb2ce5485cbf1888e6` |
| step 3d's numbers (the cells, S8) | `docs/notes/P1B_bouvard_attribution.json` | `3a8f7ea23449375110fe213d590575cfd0d08f196570b86043e2a08f851f51f3` |
| the new rows, as read and as used (1,855 cells; one difference between the entries, decided on the band) | `…/bouvard_1821/bouvard_1821_uranus_table_rows_route_b.csv` | `6b2efe05cf1f7844755bbc1e132ce441d204a2a1779c1ba8cfff1ee4d0cb07c1` |
| the reconciliation of the two entries | `…/bouvard_1821/bouvard_1821_route_b_reconciliation.csv` | `3fd4b42239d77cdcee976a2a572024febd011e444c30d82533c24956ccf5d098` |
| the two entries | `…/passes/bouvard_1821_table_rows_route_b_pass1.csv`; `…_pass2.csv` | `540984f284ba7d1e2c5578c68843edf54fffebfa026bc908415e689940b1be76`; `faea3e0e47f3b38adf7f6ea54b0ba6192ed76fcb4f96507eea5c8285544225e8` |

The new rows are in a file of their own, so that the rows file and the misprint log of step 3d
keep the hashes above. Readings of the new rows that the printed differences contradict are
logged in a new file, `…/bouvard_1821/bouvard_1821_misprints_route_b.csv`, under the rule of § 2,
before any longitude is formed; it is an output of this step, not an input.

The runner refuses to run if a hash differs.

**The epochs** are those of step 3d, unchanged: 1810.328 + 3.0362 t years with t the
`t_h1` (`t_h2` on Hypothesis II) of the normal places file, carried to Paris mean time by
`fork.historical.bouvard_1821.t_of_jd`. 1690 is set aside in every solve, as he did.

**The 1780 place is primary as the memoir has it: interpolated.** Adams says its error
was "found by interpolating between the errors of 1781, 1782, and 1783, and those given by
the Ancient Observations of 1769 and 1771". As in S8(c) of step 3d, D at the 1780 place is
the cubic through D at 1769, 1771, 1783 and 1786, evaluated at the epoch of the 1780 place,
with D now formed from the tables. **This choice is made after step 3d's results were
seen**: there S8 was a sensitivity defined post hoc, and variant (c) gave a different
reading from the registered variant (a). The readings below therefore cannot be presented
as blind to that outcome, and § 5 compares like with like for that reason.

## 2 C_tables: the table path

`C_tables(t)` = the heliocentric longitude of Uranus from the tables at t, of the mean
equinox of date, converted from the centesimal division, plus the five terms of § 7
exactly as in step 3d (arguments from the mean longitudes of p. ii, without the great
inequalities, from the fixed equinox of 1800).

The table path is `fork.historical.bouvard_1821.table_longitude` as validated in V1 and V1′
of step 3d, with its interpolation unchanged (four rows in Tables IX, X, XII and the
precession table, two elsewhere), the readings of the rows file's `value_used`, and one
extension, because the epochs are not midnights:

| step | tables | rule |
|---|---|---|
| epoch | I; II for a date before 1780 (the row of Table I a century later) | as in V1 |
| month, day | III, IV | as in V1 |
| time of day | V (decimal hours), VI (minutes; seconds as hundredths of its entries) | **new**: added to the mean longitude only; the tables give no motion of the perihelion, the node or the arguments for a fraction of a day |
| inequality of the precession | Jupiter's Table III, p. 12 | as in V1′ (interpolated by the same rule on both sides of R1, so R1 does not test it) |
| great inequality, corrections of the arguments | IX | as in V1 |
| equation of the centre, secular variation | X | as in V1 |
| equations of the longitude | XI to XXI | as in V1 |
| reduction to the ecliptic | XXXII | as in V1 |

If an epoch's argument, as the tables give it, falls on a row that was not entered, the
work stops at that epoch; the row is entered twice under a dated amendment; nothing is
filled in from the formula.

A reading in the new rows that the tables' printed differences contradict is handled by
the rule of step 3d's first amendment, quoted: "A reading that the table's printed first
differences on both sides contradict is replaced only when both differences give the same
value; the difference is reported both ways. If the two ways give different outcomes at an
instant, that instant is 'stop and ask'." Every such reading is logged in
`bouvard_1821_misprints_route_b.csv` before any longitude is formed.

## 3 R1: C_tables against C_formula at the epochs

C_formula is step 3d's C_Bouvard: the corrected formula (`corrected_formula_at`) plus § 7.
The § 7 terms are identical on both sides and cancel; R1 is the difference of the two
longitudes of Bouvard, formula minus tables, at each of the 29 epochs that have rows of
their own, and at the 1780 place computed at its own date on each branch (S8(a), S8(b)),
reported beside.

**The rounding budget, stated before any difference is evaluated.** As in V1 of step 3d:
"the number of table entries summed, times half a unit of each entry's last printed
digit". Every table in the path is printed to 0″.1 centesimal, so the half unit is 0″.05
centesimal = 0″.0162, and

    budget(epoch) = 0″.0162 × N,
    N = 17 + n_II + n_III + n_IV + n_V + n_VI,min + n_VI,sec

with 17 = the epoch (Table I), the inequality of the precession, the great inequality
(IX), the equation of the centre and its secular variation (X), the eleven equations (XI
to XXI) and the reduction (XXXII); and each n equal to 1 if that table is opened at the
epoch and 0 if not (II for a date before 1780; III unless the date is in January; IV
unless it is the first of the month; V unless the hour is 0; VI once for the minutes
unless 0 and once for the seconds unless 0). "A value interpolated between two printed
entries counts once." N runs from 17 to 23, the budget from 0″.275 to 0″.373 (step 3d's 17 for 1 January 1783 and 20 for the 1715 instant are the same count without the time of day). If a table
proves to be printed to another last digit, that entry's half unit is used in the sum; the
rule does not change.

**The argument allowance**, as in V1: for each of the arguments II to XII, the local slope
of its equation (the difference of the two bracketing entries over the step) × half a part
× the number of whole-number entries summed into the argument (I, II, III, IV, and the
correction of Table IX), summed over the eleven arguments. It is computed from the
entered rows by the code that computed it for V1 and V1′, unchanged.

**Rule**, that of V1′ (the amendment of step 3d), unchanged:

| at an epoch, the absolute difference is | outcome |
|---|---|
| ≤ budget | inside the rounding |
| over the budget and ≤ budget + allowance | inside the rounding of the arguments |
| over budget + allowance and ≤ 2″ | outside: the formula does not stand for the tables at that epoch |
| over 2″ | outside, and the entry of every row that epoch uses is checked against the page image again before anything else is done (a difference of that size is more likely a slip of entry or of row than a property of the tables) |

R1's reading: the count of epochs in each row, the largest difference, and the list of
epochs outside. R1 does not gate R2: the cells are run with C_tables whatever R1 gives,
because C_tables is what the plan is about. Stated in advance: 12 of the epochs lie within
a day of an instant already compared in V1 or V1′ and share most of their rows with it, and
the three corrections of the formula were found from five of those; R1 is a fresh test at
the other 17 only, and the table says which are which.

## 4 R2: the four cells with C_tables

Step 3d's § 2 design, with D(t) = C_tables(t) − C_clean(t) in place of
C_Bouvard(t) − C_clean(t) and the 1780 place as § 1 has it; everything else is step 3d's
code unchanged: cells A (his printed § 9 column) and D (step 3b's clean column) do not
contain Bouvard and must reproduce the gates of V3 within 0.001° (+1.647°, +1.021°,
−9.512°, −9.923°; Hypothesis II +4.964°, +4.447°, −7.621°, −8.127°), or the work stops.
Cells B and C are solved by both methods on both hypotheses. R, Ob and I as in § 4 of
step 3d.

**Reading rule: that of step 3d's § 6, quoted, not restated.**

> Read on Hypothesis I with D1(c). The shift is +11.16°. An effect is signed: positive is
> toward and past Neptune, the direction of the shift.
>
> | effect | class |
> |---|---|
> | ≥ +6° | carries the shift |
> | over +3° and under +6° | partial |
> | −3° to +3° | does not carry it |
> | under −3° | works against it |

with its table of rows (reference carries and observations do not: "The reference carries
the shift"; does not and carries: "The observations carry the shift"; carries and partial:
"Mostly the reference"; partial and carries: "Mostly the observations"; partial and
partial: "Both"; and the two rows with "works against") and its overriding rows:

> - **If the interaction |I| is 3° or more**, the two effects are not additive. The four
>   cells are reported, the effect of each factor is given at each level of the other, and
>   no single attribution is stated.
> - **If D1(c) and his chain fall in different rows**, both rows are given and the reading
>   is "depends on the method".
> - **Hypothesis II** is reported with the same classes against its own shift (+12.59°). It
>   enters no reading; a different class from Hypothesis I is said beside the reading.

As in step 3d, an effect within 1° of a class boundary (3° or 6°) is said to be so beside
the reading.

## 5 R3: the comparison with step 3d

Because the 1780 place is now primary as interpolated, like is compared with like:

| comparison | with tables (this step) | with the corrected formula (step 3d, as recorded) |
|---|---|---|
| **primary** | R2, the 1780 place interpolated | S8(c): Hyp. I, D1(c): R +2.18°, Ob +8.98°, I −1.26°; chain: R +2.39°, Ob +8.56° |
| beside it | the 1780 place on the century branch, S8(a) | the registered run: D1(c) R +5.26°, Ob +5.89°; chain R +10.34°, Ob +0.61°; reading "depends on the method" |
| beside it | the 1780 place on Table I's own branch, S8(b) | D1(c) R +1.48°, Ob +9.68°; chain R +0.47°, Ob +10.47° |

**"Same reading", fixed in advance.** For a given treatment of the 1780 place the reading
with the tables is *the same* as with the formula when all three hold on Hypothesis I:

1. D1(c)'s pair of classes (R, Ob) falls in the same row of § 6;
2. the override "depends on the method" applies in both or in neither (and, where it does
   not apply, his chain falls in that same row in both);
3. the interaction override (|I| ≥ 3°) applies in both or in neither.

Otherwise it is *a different reading*, and the report names which of the three differs.
Reported beside and not part of the rule: the change of each cell's offset (B and C, both
methods, both hypotheses), the change of R and Ob, and whether Hypothesis II changes
class.

What each outcome means:

| outcome on the primary comparison | meaning |
|---|---|
| same | the corrected formula stood for the tables well enough for step 3d's cells; step 3d's numbers stand as a description of what the tables give |
| different | differences between formula and tables of the size R1 reports move the attribution across a class. This is one more instance of the sensitivity step 3d found (2″ to 3″ at one place moved a cell by 2° to 5°), and it is not by itself evidence for either attribution |

Whatever R3 gives, step 3d's registered reading ("depends on the method"; the shift is not
attributed) is not replaced by this step: its plan fixed the century branch at 1780, and
the choice made here came after its results.

## 6 Reported, not read

- **S-a, S-b:** cells B and C with the 1780 place at its own date on the century branch
  and on Table I's own branch (the rows are in the supplement of the entries). If the
  supplement has not been entered, they are reported as not run.
- **S-printed:** cells B and C with every reading as printed (`value_as_printed`) in place
  of `value_used`, the seven logged cases and any new one. It answers S5 and S5b of step 3d
  for the tables as a computer would have met them.
- **S-two-row:** cells B and C with Tables IX and X interpolated as the book's own
  apparatus does it: two rows with the printed differences, Table IX's proportional part
  corrected for second differences by Table VII with the decimal of the year from Table
  VIII, Table X's by its printed second differences, Table XII between two rows. It needs
  the difference columns of the entries (the second supplement for the rows already in the
  file); not run if they are missing.
- **S-midnight:** cells B and C with each epoch moved to the nearest Paris midnight
  (Tables V and VI not used), as V1 and V1′ took their instants.
- R1's parts at each epoch (mean longitude, great inequality, equation of the centre,
  perturbations, reduction), as in step 3d's tables.

None of these enters a reading row.

## 7 What this cannot show

- Which copy of the tables Greenwich or Adams used, with what manuscript corrections (the
  Harvard copy carries a handwritten "+0″,10 t" at Table I, not applied here), or how
  their computers interpolated and rounded. C_tables is one careful reading of one copy.
- Whether the 1780 place was interpolated as the cubic of § 1 does it. His interpolation
  used errors of 1781 to 1783 that this study does not have.
- Anything about the observations. O_real is still C + his printed residuals; the
  errors of observation, the reductions and the construction of the normal places stay
  unseparated, as in § 8 of step 3d.
- An attribution of the shift. Step 3d's interpretation stands: changes of a few
  arcseconds at single early places move his solution by 2° to 5°. Agreement of tables
  and formula (R3 "same") would remove one doubt about step 3d's inputs, not that
  sensitivity.
- Significance. The classes describe how the realised columns behave under his method.

## 8 Cost, order, layout

- **Entries:** 1,531 cells for the 29 places (eight parts of 161 to 231), 76 for the 1780
  place at its own date, 248 printed differences beside rows already in the file; each
  twice, then compared and reconciled, with the misprint log. This is where the work is.
- **Compute:** step 3d records "about six minutes on one core" for "4 cells and 12
  variants, each by both methods and both hypotheses": 4 + 12 × 2 = 28 cell solves on each
  of two hypotheses, about 6 to 7 seconds a solve. R2 is the four cells on two hypotheses
  (8 solves, about one minute); the four sensitivities of § 6 add 16 solves (about two
  minutes). R1 is table lookups: under a second. **Under five minutes on one core, not
  detached, no sampling.**
- **Order:** the entries, their comparison and the reconciled rows file committed with
  hashes (done); this plan with those hashes committed and pushed; the table path's extension with tests (the six plus seven instants of V1 and V1′
  must reproduce their recorded table longitudes exactly as before); R1; the gates; R2; R3;
  the sensitivities. One commit for the plan, one for the entries, one for code and
  results.
- **Layout:** `fork.historical.bouvard_1821.table_longitude_at(t)` (the path of § 2 at any
  epoch), the runner `fork.export.p1b_bouvard_tables`, tests in
  `tests/test_p1b_bouvard_tables.py` (the entries agree; the V1 and V1′ instants unchanged;
  the budget formula on constructed dates; the "same reading" rule on constructed effects),
  results in `docs/notes/P1B_bouvard_tables.{md,json}`. The site and the earlier notes are
  untouched.

## 9 Open choices, each resolved conservatively

Each was open when this plan was drafted; each is resolved as the table says, and none was
resolved after any number of this step was seen.

| # | choice | resolution in this draft | why it is the conservative one | to reverse |
|---|---|---|---|---|
| 1 | the 1780 place | interpolated (S8(c)) is primary, declared as chosen after step 3d's results; (a) and (b) beside | it is what the memoir says, and R3 compares each treatment with its own counterpart so the post hoc choice cannot manufacture a change of reading | make (a) primary: swap the rows of § 5 |
| 2 | interpolation in the tables | the rule of V1 and V1′ unchanged (four rows in IX, X, XII); the book's two-row apparatus as S-two-row | the only rule already validated against the formula; no new freedom | make S-two-row primary; needs the second supplement entered. Tables VII and VIII (59 cells of the eight parts) serve S-two-row only |
| 3 | the instants | exactly the epochs of step 3d, with Tables V and VI for the time of day; the nearest midnight as S-midnight | D is formed at the same instants as C_clean and as step 3d | drop V and VI (91 fewer cells to enter) |
| 4 | arguments at a fraction of a day | no motion beyond Table IV | the tables print none | none needed |
| 5 | readings contradicted by the differences | step 3d's rule and rulings; S-printed beside | already ruled; both ways reported | — |
| 6 | R1's rule | V1′'s (budget, then budget + allowance), with 2″ as the level at which entries are rechecked | no new threshold | — |
| 7 | does R1 gate R2? | no | C_tables is the object; a gate would let the formula veto the tables | add a stop after R1 |
| 8 | "same reading" | same § 6 row for D1(c) on Hypothesis I, same state of the two overrides | uses only the registered classes; no new tolerance in degrees | add a numerical tolerance |
| 9 | the § 7 terms | as step 3d's primary | unchanged inputs | S6, S7 of step 3d could be rerun |
| 10 | the handwritten "+0″,10 t" | not applied; the printed table is represented | it is not print, and its date and hand are unknown | a sensitivity adding 0″.10 (Y − 1780) centesimal to the mean longitude |
| 11 | margin of rows entered | wider than 2 % for arguments II to XII (3 parts; 10 for III); stop and enter if a row is still missing | never fill a table value from the formula | — |
| 12 | Hypothesis II | reported, enters no reading | as step 3d | — |
| 13 | when the plan is frozen | after the entries exist (their hashes are inputs), before any longitude is formed from them | the entries are transcription; no number a reading uses exists before the freeze | freeze first with the template's hash only, and add the entries' hashes by amendment |
| 14 | Table VI for the seconds | entered and used, counted once in the budget | the book says the table serves for minutes and seconds; at most 0″.13 centesimal | round to the minute |
| 15 | where the new rows go | a file of their own; a misprint log of their own | step 3d's inputs keep their hashes, and its released results stay reproducible | append to the rows file of step 3d |
| 16 | a cell whose two bands print different signs (Table X, row 200, secular variation: −0″0 in the first, +0″0 in the second) | the first band, as both entries read it | the value is zero either way, so the sign changes no number; the rule for the band was fixed in the entry specification before any entry | none needed |

---

## Results, 13 November 2026

Run after the plan above was committed (f124d70). Nothing above the line is changed.
Runner: `uv run python -m fork.export.p1b_bouvard_tables` (167 s on one core, not
detached, no sampling); numbers in `P1B_bouvard_tables.json`. Commits: the misprint log
(df3d312), the table path and the runner with their tests (cde9266), these results.

**In one sentence:** the tables agree with the corrected formula within the rule of § 3 at
all 29 places (largest −1″.07, none outside), the gates pass, and with the tables the four
cells read "Mostly the observations, with a share from the reference" (R +3.16°, Ob +7.99°),
which by § 5 is *a different reading* from the formula's "The observations carry the shift"
(condition 1).

### The inputs

All twelve inputs of § 1 have their recorded SHA-256, the scan included; the runner checks
them before it does anything else.

### The new readings against their printed differences (before any longitude)

Every reading with a printed difference beside it was compared with that difference, by
magnitude (the tables print a sign on the first line of a run only), in tenths of a
centesimal second, before any table longitude was formed (df3d312).

| table | new readings checked on both sides | on one side (the row beyond was not entered) | not checkable | readings of step 3d's rows checked against the newly entered differences |
|---|---|---|---|---|
| IX | — | — | — | 19 |
| X | 33 | 32 | 0 | 48 |
| XI | 19 | 17 | 0 | 24 |
| XII | 41 | 6 | 0 | 41 |
| XIII | 21 | 11 | 0 | 24 |
| XIV | 26 | 19 | 0 | 26 |
| **total** | **140** | **85** | **0** | **182** |

By the rule § 2 quotes ("A reading that the table's printed first differences on both sides
contradict is replaced only when both differences give the same value"): **no reading of
the new rows is contradicted on both sides, and none is replaced.** `value_used` is
`value_as_printed` throughout the new rows, so no epoch is "stop and ask" under § 2.

`bouvard_1821_misprints_route_b.csv` (step 3d's columns) logs eleven cells:

| place | reads | by the evidence | what it is |
|---|---|---|---|
| IX, between 1790 and 1800, Différences | 25,4 | 25,6 | a printed difference both its readings contradict; each reading agrees with its other difference (21,5; 29,0) |
| X, between 84 and 85, first differences | 1.74,7 | 1.74,8 | the same kind; the printed second differences on 84 and 85 (15,0; 14,9) agree with 1.74,7 |
| X, between 142 and 143, first differences | 5.87,5 | 5.88,5 | the same kind |
| XIII, between 40 and 50, Différences | 9,8 | 9,6 | the same kind |
| XIII, row 720, Équation | 2,6 | 2,5 | a reading of step 3d's rows, contradicted on both sides by the differences now entered (2,3 and 1,4), both giving 2,5. Step 3d logged it and left it as read; that ruling stands (§ 9, choice 5) |
| X, second differences on 4, 83, 101, 141, 143, 189 | 1,2; 14,9; −14,5; 10,1; 9,9; 2,1 | 1,3; 14,8; 13,5; 10,3; 10,1; 2,0 | printed second differences the first differences on either side contradict (by 0,1 to 0,2, and by 1,0 on 101, where row 102 was not entered) |

None of the eleven is in the primary path: it uses readings, not differences, and XIII 720
enters it as step 3d used it (at 1771 and 1813). Four touch S-two-row only: IX 1790 to
1800 (at 1792, 1795, 1798), X 84 to 85 (1816), X 142 to 143 and the second difference on
143 (both 1828). They are used there as printed.

### R1: C_tables against C_formula at the epochs

Formula minus tables, heliocentric longitude of Uranus. N is the count of § 3; the budget
is 0″.0162 × N; the allowance is the argument allowance computed by V1's code.

| place | Paris date, decimal time | N | budget | budget + allowance | formula − tables | outcome | within a day of |
|---|---|---|---|---|---|---|---|
| 1690 | 1690-12-20, 1h 17m 21s | 23 | 0″.373 | 1″.65 | −1″.07 | inside the rounding of the arguments | 1690 (V1′) |
| 1712 | 1712-04-03, 0h 11m 28s | 22 | 0″.356 | 1″.51 | +0″.33 | inside the rounding | — |
| 1715 | 1715-03-03, 6h 24m 45s | 23 | 0″.373 | 1″.35 | +0″.77 | inside the rounding of the arguments | 1715 (V1) |
| 1750 | 1750-11-05, 0h 46m 03s | 22 | 0″.356 | 1″.23 | +0″.66 | inside the rounding of the arguments | 1750 (V1′) |
| 1753 | 1753-11-29, 1h 07m 80s | 23 | 0″.373 | 1″.74 | +0″.15 | inside the rounding | — |
| 1756 | 1756-08-23, 1h 82m 65s | 23 | 0″.373 | 1″.63 | +0″.12 | inside the rounding | — |
| 1764 | 1764-01-08, 9h 84m 73s | 22 | 0″.356 | 1″.03 | −0″.08 | inside the rounding | — |
| 1769 | 1769-01-11, 7h 88m 61s | 22 | 0″.356 | 1″.22 | −0″.43 | inside the rounding of the arguments | 1769 (V1′) |
| 1771 | 1771-12-12, 4h 01m 78s | 23 | 0″.373 | 0″.76 | −0″.49 | inside the rounding of the arguments | — |
| 1783 | 1782-12-31, 0h 60m 04s | 21 | 0″.340 | 1″.54 | +0″.04 | inside the rounding | 1783 (V1) |
| 1786 | 1786-01-13, 0h 32m 09s | 20 | 0″.324 | 0″.77 | −0″.20 | inside the rounding | — |
| 1789 | 1789-01-26, 0h 04m 14s | 20 | 0″.324 | 1″.13 | −0″.49 | inside the rounding of the arguments | 1789 (V1′) |
| 1792 | 1792-02-08, 9h 76m 19s | 22 | 0″.356 | 0″.96 | +0″.04 | inside the rounding | — |
| 1795 | 1795-02-21, 9h 48m 24s | 22 | 0″.356 | 1″.35 | −0″.03 | inside the rounding | — |
| 1798 | 1798-03-06, 9h 20m 29s | 22 | 0″.356 | 1″.40 | −0″.01 | inside the rounding | 1798 (V1) |
| 1801 | 1801-03-20, 8h 92m 34s | 22 | 0″.356 | 1″.20 | −0″.36 | inside the rounding of the arguments | — |
| 1804 | 1804-04-02, 8h 64m 39s | 22 | 0″.356 | 1″.38 | −0″.06 | inside the rounding | 1804 (V1′) |
| 1807 | 1807-04-16, 8h 36m 44s | 22 | 0″.356 | 0″.80 | −0″.02 | inside the rounding | — |
| 1810 | 1810-04-29, 8h 08m 49s | 22 | 0″.356 | 1″.32 | +0″.21 | inside the rounding | — |
| 1813 | 1813-05-12, 7h 80m 54s | 22 | 0″.356 | 0″.76 | −0″.23 | inside the rounding | 1813 (V1) |
| 1816 | 1816-05-25, 7h 52m 59s | 22 | 0″.356 | 1″.09 | −0″.48 | inside the rounding of the arguments | — |
| 1819 | 1819-06-08, 7h 24m 64s | 22 | 0″.356 | 0″.97 | −0″.94 | inside the rounding of the arguments | — |
| 1822 | 1822-06-21, 6h 96m 69s | 22 | 0″.356 | 1″.23 | −0″.21 | inside the rounding | 1822 (V1′) |
| 1825 | 1825-07-04, 6h 68m 74s | 22 | 0″.356 | 1″.39 | −0″.30 | inside the rounding | — |
| 1828 | 1828-07-17, 6h 40m 79s | 22 | 0″.356 | 1″.16 | −0″.45 | inside the rounding of the arguments | 1828 (V1) |
| 1831 | 1831-07-31, 6h 12m 84s | 22 | 0″.356 | 1″.50 | −0″.58 | inside the rounding of the arguments | — |
| 1834 | 1834-08-13, 5h 84m 89s | 22 | 0″.356 | 0″.74 | −0″.23 | inside the rounding | — |
| 1837 | 1837-08-26, 5h 56m 94s | 22 | 0″.356 | 1″.37 | −0″.08 | inside the rounding | 1837 (V1′) |
| 1840 | 1840-09-08, 5h 28m 99s | 22 | 0″.356 | 1″.02 | −0″.35 | inside the rounding | — |
| *beside:* 1780, S8(a), century branch | 1779-12-18, 0h 87m 99s | 22 | 0″.356 | 1″.27 | −0″.15 | inside the rounding | |
| *beside:* 1780, S8(b), Table I's own branch | 13 days 9h 12m 01s before 1780 | 21 | 0″.340 | 0″.89 | −0″.09 | inside the rounding | |

On Table I's own branch the formula is the corrected formula with the mean motion's slip
counted from 1780 without the century, as S8(b) of step 3d had it.

**R1's reading**, by § 3's rule ("≤ budget | inside the rounding"; "over the budget and ≤
budget + allowance | inside the rounding of the arguments"):

| epochs | inside the rounding | inside the rounding of the arguments | outside, ≤ 2″ | over 2″ | largest difference |
|---|---|---|---|---|---|
| **the 17 fresh** (1712, 1753, 1756, 1764, 1771, 1786, 1792, 1795, 1801, 1807, 1810, 1816, 1819, 1825, 1831, 1834, 1840) | 12 | 5 (1771, 1801, 1816, 1819, 1831) | 0 | 0 | −0″.94 (1819) |
| the 12 within a day of a V1 or V1′ instant | 6 | 6 (1690, 1715, 1750, 1769, 1789, 1828) | 0 | 0 | −1″.07 (1690) |
| all 29 | 18 | 11 | 0 | 0 | −1″.07 (1690) |

**No epoch is outside.** At the 17 fresh epochs, the only fresh test of § 3, the corrected
formula stands for the tables at every one; 12 of the 17 are inside the rounding alone. No
epoch is over 2″, so no entry had to be checked against the page images again.

R1's parts (centesimal seconds, formula minus tables; × 0.324 for ″):

| place | mean longitude | great inequality | equation of the centre | perturbations | reduction |
|---|---|---|---|---|---|
| 1690 | −0.16 | −0.07 | +0.10 | −3.15 | −0.03 |
| 1712 | +0.03 | −0.05 | −0.20 | +1.26 | −0.02 |
| 1715 | −0.00 | −0.05 | −0.01 | +2.41 | +0.03 |
| 1750 | +0.03 | −0.05 | +0.04 | +1.71 | +0.30 |
| 1753 | −0.04 | −0.05 | +0.08 | +0.49 | −0.03 |
| 1756 | +0.00 | −0.04 | −0.00 | +0.36 | +0.04 |
| 1764 | −0.05 | −0.01 | +0.14 | −0.28 | −0.05 |
| 1769 | +0.11 | +0.01 | +0.26 | −1.72 | −0.00 |
| 1771 | +0.08 | +0.03 | +0.21 | −1.82 | −0.02 |
| 1783 | −0.01 | +0.03 | +0.27 | −0.14 | −0.03 |
| 1786 | −0.12 | +0.00 | +0.20 | −0.68 | −0.02 |
| 1789 | −0.03 | −0.02 | +0.12 | −1.58 | −0.01 |
| 1792 | −0.04 | −0.00 | −0.09 | +0.20 | +0.05 |
| 1795 | −0.05 | +0.03 | +0.13 | −0.17 | −0.03 |
| 1798 | −0.06 | +0.07 | −0.14 | +0.07 | +0.03 |
| 1801 | +0.03 | +0.08 | −0.00 | −1.21 | +0.01 |
| 1804 | −0.08 | +0.06 | −0.08 | −0.04 | −0.03 |
| 1807 | −0.09 | +0.02 | −0.15 | +0.12 | +0.04 |
| 1810 | −0.10 | −0.02 | −0.31 | +1.13 | −0.04 |
| 1813 | −0.01 | −0.05 | −0.37 | −0.30 | +0.02 |
| 1816 | +0.08 | −0.06 | −0.44 | −1.11 | +0.04 |
| 1819 | +0.06 | −0.07 | −2.74 | −0.14 | −0.00 |
| 1822 | −0.04 | −0.08 | −0.31 | −0.18 | −0.02 |
| 1825 | −0.06 | −0.09 | −0.29 | −0.48 | −0.02 |
| 1828 | −0.06 | −0.09 | −0.63 | −0.61 | +0.00 |
| 1831 | −0.08 | −0.10 | −0.16 | −1.50 | +0.05 |
| 1834 | −0.09 | −0.12 | −0.05 | −0.69 | +0.22 |
| 1837 | +0.00 | −0.13 | +0.00 | −0.08 | −0.03 |
| 1840 | −0.11 | −0.14 | +0.02 | −0.90 | +0.04 |

Every part over 1 centesimal second is in the perturbations (the largest −3.15 at 1690),
except one: the equation of the centre at 1819 (−2.74), an anomaly of 99.2 whose four rows
(98 to 101) straddle the line of 100 that ends the first half of Table X.

Post hoc, beside R1 and not part of its rule: with step 3d's four replaced readings taken
as printed (Table X rows 324, 327, 25; XII 2600), formula minus tables is +2″.94 at 1783,
−2″.39 at 1789, +0″.67 at 1804 and −1″.08 at 1690. These are step 3d's rulings, not
readings of the new rows, so § 2's "stop and ask" does not apply to them; they enter
S-printed below.

### The gates

| check | must give | gives | |
|---|---|---|---|
| cell A, D1(c), Hypothesis I | +1.647° | +1.6469° | pass |
| cell A, his chain | +1.021° | +1.0211° | pass |
| cell D, D1(c) | −9.512° | −9.5115° | pass |
| cell D, his chain | −9.923° | −9.9230° | pass |
| Hypothesis II, cell A, D1(c) and chain | +4.964°, +4.447° | +4.9641°, +4.4470° | pass |
| Hypothesis II, cell D, D1(c) and chain | −7.621°, −8.127° | −7.6206°, −8.1267° | pass |

Also checked before the cells: step 3d's C_Bouvard put through this runner's columns gives
step 3d's recorded columns B and C exactly (largest difference 0″.000).

### R2: the four cells with C_tables

D with the tables less D with the formula (= tables minus formula, the § 7 terms
cancelling) runs from −0″.77 to +1″.07 at the 29 places (RMS 0″.43). At the 1780 place,
interpolated by the cubic of § 1 on both sides, the two differ by +0″.01.

Offset from Neptune at 1846.762, degrees:

| cell | O | C | Hyp. I, D1(c) | Hyp. I, chain | Hyp. II, D1(c) | Hyp. II, chain |
|---|---|---|---|---|---|---|
| A | real | Bouvard's tables | +1.647 | +1.021 | +4.964 | +4.447 |
| B | real | clean | −0.480 | −2.584 | +1.606 | −0.768 |
| C | clean | Bouvard's tables | −5.309 | −4.793 | −2.742 | −1.954 |
| D | clean | clean | −9.512 | −9.923 | −7.621 | −8.127 |

The chain found two positive-mass roots in every cell.

| | shift | R, the reference | Ob, the observations | I | classes | row of § 6 |
|---|---|---|---|---|---|---|
| **Hyp. I, D1(c) (the reading)** | +11.16 | **+3.16** | **+7.99** | −1.04 | partial, carries | Mostly the observations |
| Hyp. I, chain | +10.94 | +4.37 | +6.58 | −0.76 | partial, carries | Mostly the observations |
| Hyp. II, D1(c) | +12.58 | +4.12 | +8.47 | −0.76 | partial, carries | Mostly the observations |
| Hyp. II, chain | +12.57 | +5.69 | +6.88 | −0.48 | partial, carries | Mostly the observations |

**R2's reading**, by the rule § 4 quotes from step 3d's § 6: R is "over +3° and under +6° |
partial" and Ob is "≥ +6° | carries", so the row is **"partial and carries: 'Mostly the
observations'"**, with a share from the reference: +8.0° from the observations and +3.2°
from the reference.

- The interaction is −1.04°, under 3°: that override does not apply.
- D1(c) and his chain fall in the same row: "depends on the method" does not apply.
- Hypothesis II falls in the same classes.
- Beside the reading, as § 4 asks: D1(c)'s R (+3.16°) lies within 1° of the 3° boundary,
  0.16° above it. Also within 1° of a boundary, beside it: the chain's Ob on Hypothesis I
  (+6.58°) and both of the chain's effects on Hypothesis II (+5.69°, +6.88°), near 6°.

### R3: the comparison with step 3d

| comparison | tables (this step): D1(c); chain | formula (step 3d, as recorded): D1(c); chain | "same reading" (§ 5) |
|---|---|---|---|
| **primary: 1780 interpolated** | R +3.16, Ob +7.99 (partial, carries); R +4.37, Ob +6.58 (partial, carries) | S8(c): R +2.18, Ob +8.98 (does not carry, carries); R +2.39, Ob +8.56 (does not carry, carries) | **different**: condition 1 |
| beside: 1780 on the century branch | S-a: R +6.04, Ob +5.11 (carries, partial); R +11.41, Ob −0.47 (carries, does not carry) | the registered run: R +5.26, Ob +5.89 (partial, partial); R +10.34, Ob +0.61 (carries, does not carry) | different: condition 1 |
| beside: 1780 on Table I's own branch | S-b: R +2.34, Ob +8.82 (does not carry, carries); R +2.07, Ob +8.87 (does not carry, carries) | S8(b): R +1.48, Ob +9.68; R +0.47, Ob +10.47 (both does not carry, carries) | same |

**R3's reading on the primary comparison: a different reading.** Condition 1 differs:
D1(c)'s pair of classes falls in the row "Mostly the observations" with the tables and in
"The observations carry the shift" with the formula. Condition 2: "depends on the method"
applies in neither; his chain moves with D1(c) (from "does not carry, carries" to
"partial, carries"), so the clause of condition 2 on the chain's row differs as well, as it
must when condition 1 does and the override applies in neither. Condition 3: the
interaction override applies in neither.

The meaning § 5 fixed for this outcome, quoted: "differences between formula and tables of
the size R1 reports move the attribution across a class. This is one more instance of the
sensitivity step 3d found (2″ to 3″ at one place moved a cell by 2° to 5°), and it is not by
itself evidence for either attribution". The class crossed is the 3° boundary of R, crossed
by 0.16°.

On the comparisons beside it: on the century branch the reading also differs by condition
1 (with the tables "Mostly the reference", with the formula "Both"; "depends on the method"
in both); on Table I's own branch it is the same ("The observations carry the shift" in
both, by both methods).

As § 5 says, step 3d's registered reading ("depends on the method"; the shift is not
attributed) is not replaced by this step.

Reported beside, not part of the rule: change, tables minus formula, degrees.

| comparison | | B | C | R | Ob |
|---|---|---|---|---|---|
| primary | Hyp. I, D1(c) | −1.21 | +0.76 | +0.99 | −0.99 |
| | Hyp. I, chain | −2.53 | +1.44 | +1.98 | −1.98 |
| | Hyp. II, D1(c) | −1.58 | +0.97 | +1.28 | −1.28 |
| | Hyp. II, chain | −3.24 | +1.89 | +2.56 | −2.56 |
| century branch | Hyp. I, D1(c) | −1.00 | +0.56 | +0.78 | −0.78 |
| | Hyp. I, chain | −1.82 | +0.33 | +1.08 | −1.08 |
| | Hyp. II, D1(c) | −1.21 | +0.80 | +1.00 | −1.00 |
| | Hyp. II, chain | −1.99 | +0.75 | +1.37 | −1.37 |
| Table I's own branch | Hyp. I, D1(c) | −1.03 | +0.68 | +0.86 | −0.86 |
| | Hyp. I, chain | −1.92 | +1.27 | +1.60 | −1.60 |
| | Hyp. II, D1(c) | −1.41 | +0.85 | +1.13 | −1.13 |
| | Hyp. II, chain | −2.70 | +1.59 | +2.15 | −2.15 |

Hypothesis II changes class on the primary comparison (D1(c)'s R from +2.84°, "does not
carry", to +4.12°, "partial") and on Table I's own branch (+1.98° to +3.11°); not on the
century branch.

### § 6: reported, not read

Change of cells B and C from R2, degrees, Hypothesis I, and the attribution each variant
gives (D1(c); chain).

| | δ in D | B, D1(c) | C, D1(c) | B, chain | C, chain | R; Ob (D1(c)) | R; Ob (chain) | row of § 6 by D1(c) |
|---|---|---|---|---|---|---|---|---|
| S-a, 1780 at its own date, century branch | −2″.66 at 1780 | −3.69 | +2.07 | −10.15 | +3.94 | +6.04; +5.11 | +11.41; −0.47 | Mostly the reference; depends on the method |
| S-b, 1780 at its own date, Table I's branch | +0″.72 at 1780 | +1.08 | −0.57 | +2.99 | −1.61 | +2.34; +8.82 | +2.07; +8.87 | The observations carry the shift |
| S-printed, every reading as printed | −2″.91 at 1783, +1″.90 at 1789, −0″.73 at 1804, −3″.30 at 1780 (through the cubic) | −7.11 | +4.10 | +166.29 (the other root) | +3.58 | +8.77; +2.39 | −76.99; +87.93 | The reference carries the shift; depends on the method |
| S-two-row, the book's apparatus | −0″.07 to +0″.12 | +0.01 | −0.01 | +0.07 | −0.03 | +3.15; +8.01 | +4.32; +6.63 | Mostly the observations |
| S-midnight | not run | | | | | | | |

- **S-a and S-b** run: the supplement for the 1780 place at its own date was entered.
- **S-printed:** the readings as printed differ from those used only in step 3d's four
  replaced readings; the new rows have none. In the chain, cell B takes its other root
  (offset +163.70°), as in step 3d's S5.
- **S-two-row** run with the difference columns entered. It moves no cell by more than
  0.11° (Hypothesis II, chain, cell B).
- **S-midnight was not run.** At five epochs the nearest Paris midnight is the next civil
  day and that day's row of Table IV was not entered: 1764, 1792, 1819 and 1840 (Table IV,
  row 9), 1810 (Table IV, row 30). By § 2 the work stops at such an epoch and nothing is
  filled in from the formula; I did not enter the rows, so the variant has no value.
- R1's parts at each epoch are in the table of R1 above.

### Deviations from the plan, and choices where it was silent

Each was made before the number it touches was seen, except where marked.

1. **Rulings in a copy.** § 2 puts rulings into `value_used` "in a copy". Because the rows
   file's hash is an input, the copy is in memory: the loader applies the replacements the
   route B misprint log records. There are none.
2. **Comparison with the differences by magnitude.** The tables print a sign on the first
   line of a run only, and that line is often not among the rows entered; a reading is
   compared with the magnitude of each difference.
3. **What the log holds.** Besides readings, it logs printed first differences the
   readings at both their ends contradict, and Table X's second differences the first
   differences contradict, because S-two-row uses them. They are not readings of § 2's rule
   and nothing is replaced for them. XIII 720 is listed with step 3d's ruling, unchanged.
4. **Readings of step 3d's rows** checked against the differences now entered keep step
   3d's rulings (§ 9, choice 5). Only XIII 720 is contradicted on both sides.
5. **"Touches" in the log** was filled from the entry plan's rows by epoch (rows located
   with the corrected formula's arguments), so that no table longitude was formed before
   the log was committed.
6. **The time of day.** The decimal seconds are rounded to a whole second, as the entry
   plan did (no rounding reaches 100). The date used for Table IX, the secular variation of
   Table X and the precession is the Paris midnight of the civil date, as in V1: § 2 adds
   the time of day "to the mean longitude only". Table IX moves about 0.01 centesimal
   second in a day at most.
7. **Table I's own branch** for the 1780 place: Table I's row of 1780 less Table IV's row
   14 (13 days) and Tables V and VI for the 0.912 day left, all subtracted; Table III is not
   opened. Its formula counterpart is the corrected formula without the century in the
   mean motion's slip (+3″.40 at that place, S8(b)'s figure).
8. **"Stop and ask"** (§ 2) is applied to readings replaced under route B's rule; there are
   none. Step 3d's own replaced readings are reported both ways beside R1 (post hoc) and
   enter S-printed.
9. **Condition 2 of "same reading"** is reported in its two parts (the override; the
   chain's row where the override applies in neither). When condition 1 differs and the
   override applies in neither, the second part differs with it; the report names
   condition 1 as the cause. The outcome "same" or "different" does not depend on this.
10. **The boundary remark** of § 4 is checked against +6°, +3° and −3°, for D1(c) on
    Hypothesis I and, beside it, for the other three rows.
11. **S-two-row**, details the plan leaves open: Table IX's proportional part from its
    printed difference (the sign of the run carried down), corrected by Table VII
    interpolated between its rows (whole seconds) and columns (whole years), with the mean
    of the two second differences about the decade (Bessel's) and the year from Table VIII
    interpolated between its two dates; Table X from its printed first difference, with the
    sign the readings give, corrected by the mean of the two printed second differences of
    the interval, or the one printed where the other is not (at 1819, row 100; at 1756, row
    200); Table XII between two rows; the precession table by the rule of V1′, unchanged.
12. **The formula side of R3** is step 3d's JSON as recorded, not solved again; step 3d's
    columns were reproduced exactly through this runner (above).
13. **Hypothesis II's epochs**: `t_h2` equals `t_h1` at all 30 places in the normal places
    file, and step 3d's code forms D at `t_h1`; nothing changes.
14. **Cost**: 167 s on one core, as § 8 expected.

### What this leaves

- The questions of § 0: (R1) the corrected formula agrees with the tables within the
  tables' rounding at all 29 places, the 17 fresh ones included; (R2) with the tables the
  reading is "Mostly the observations", R +3.16° and Ob +7.99°; (R3) that is a different
  reading from the formula's, by condition 1, with R 0.16° over the 3° boundary.
- Step 3d's registered reading stands, as § 5 fixes. What § 7 says this cannot show is
  unchanged.
- S-midnight needs Table IV's rows 9 and 30 entered twice, under a dated amendment, if it
  is to be run.

## Spot check by a person, 14 November 2026

Nothing above is changed. Twenty of the 1,855 cells were checked by a person against the
page bands, in a kit drawn before any answer was given (seed [20261114, 1821]; the kit
carries page images and stays private). The answers are the kit's own export
(15 November 00:40 UTC), saved as
`docs/review/bouvard_route_b_spotcheck_human_results.json`. The reviewer field was left
empty by the export and was filled when the record was saved. The selection rule is in that
file.

| group | cells | matches |
|---|---|---|
| misprint log: every logged cell a computed path reads (all four are read by S-two-row), and one more drawn | 5 | 5 |
| the one cell on which the two entries disagreed (Table I, 1801, column VIII: 150) | 1 | 1 |
| marked doubtful (an alternative reading or a damaged figure in the entry note), drawn | 4 | 4 |
| read at the places from which the 1780 place is interpolated: one at each of 1769, 1771, 1783 and 1786, and one more | 5 | 5 |
| random, from the remaining cells | 5 | 5 |

**0 discrepancies in 20.** What this supports, and what it does not:

- The 15 chosen cells are where an error was likeliest or would have mattered most: every
  logged misprint that a computed path reads, the one disagreement and four doubtful
  readings all hold as entered.
- For the 1780 place, on which the reading turns, the check is a sample: 5 of the 124
  route-B cells the four neighbour places read. It is not a check of all of them.
- The 5 random cells bound the error rate of the whole transcription only weakly: with 0
  errors in 5, the one-sided 95 % upper bound is 45 %. The case for the rows rests mainly
  on the double entry and its reconciliation, as before.
