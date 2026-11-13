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
