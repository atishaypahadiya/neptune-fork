# Data dictionary

All CSV files are UTF-8, comma-separated, with one header row and LF line endings. A field is
quoted only when it contains a comma or a quote. Typed schemas for every file are in
`datapackage.json`.

## Conventions shared by the Table 77 files

**As printed.** Each cell of `leverrier_1846_table77.csv` holds the printed value as text.
Misprints are kept (they are listed in `leverrier_1846_misprints.csv`), and the printed number
of decimals is kept, because it shows the precision. The only normalization is typographic:

| printed | stored | rule |
|---|---|---|
| `6.51.40,38` | `6 51 40.38` | Unit separators (printed as periods) become one space; the decimal comma becomes a point |
| `+12.15,35,0` | `+12 15 35.0` | Separators are read by position (degrees, minutes, seconds), whatever mark is printed |
| `− 0,33. 2,7` | `-0 33 2.7` | 1843 Sep 20 Dec (p. 136, l. 18), printed with a comma after the degrees: −0° 33′ 2.7″ |
| `+`, `−` | `+`, `-` | Signs become ASCII |
| `» » »` | `»` | Printed in the Dec columns when there is no Dec observation |
| `+12,5` | `+12.5` | Decimal comma in single numbers |

Notes on these rules:
- **The sign belongs to the whole sexagesimal value.** Read it from the string, not from the
  first number: for `-0 33 2.7` the degree field is 0, and the value is negative.
- **Blank year and month.** The table prints the year and the month only when they change.
  A blank cell means "same as the nearest row above with a value". Every page (129–136)
  prints them on its first row, so you never have to look at the previous page.
- **`»` means no value, not "ditto".** It appears on 29 rows, all in the four Dec-related
  columns at once (`dec_obs`, `dec_calc_sec`, `dec_c_minus_o`, `lat_c_minus_o`). These are
  RA-only observations. Le Verrier still gives a longitude error for them, derived from the
  RA alone (p. 127–128).

**Two eras, two RA units.** The first 19 rows (p. 129, l. 1–19; 1690–1771, the "ancient"
observations) print RA in arc. The other 260 rows (1781–1845) print it in time.

| era | rows | `ra_obs` | `ra_calc_sec` | `ra_c_minus_o` |
|---|---|---|---|---|
| ancient | p. 129 l. 1–19 (1690–1771) | ° ′ ″, to 0.1″ (e.g. `55 49 19.7` = 55° 49′ 19.7″) | ″ of arc, 0.1 | ″ of arc, 0.1 |
| modern | the rest (1781–1845) | h m s, to 0.01 s (e.g. `6 12 26.63`) | s of time, 0.01 | s of time, 0.01 |

The Dec and longitude/latitude columns are in ° ′ ″ and ″ of arc in both eras. To convert a
modern RA difference from seconds of time to arcseconds of RA, multiply by 15. To get
arcseconds on the sky, also multiply by cos δ.

**Sign convention: the table prints C − O** (calculé moins observé, computed minus observed),
for RA, Dec, longitude and latitude. The derived file `leverrier_1846.csv` gives
`lv_residual` as **O − C**, the sign flipped (see below).

**"Computed" means Le Verrier's "théorie admise".** That is his ephemeris n° 75 (§72, p. 111):
his preliminary elliptic orbit of Uranus, perturbed by Jupiter and Saturn only. Neptune is not
in it. The C − O columns are therefore the anomaly that led to Neptune's prediction.

**Frame.** The observed and computed places are apparent geocentric RA and Dec: true equator
and equinox of date, with nutation and aberration (§72, p. 111; p. 127).

**Time.** `mean_time` is Paris mean time, in astronomical reckoning: day *D* begins at mean noon
of civil day *D*. So `Déc. 23, 9 41 25` is civil 23 December, 21:41:25. And `Sept. 25,
18 0 15` is civil 26 September, 06:00:15. The table does not name the meridian. The text
states Paris time for Mayer's 1756 row (p. 126). For every other row, Paris time is inferred
from the meridian-transit identity (see `transit_resid_s` below): Le Verrier converted the
Greenwich times to Paris mean time.

## `leverrier_1846_table77.csv`

The reconciled table. One row per printed row: 279 rows.

| column | type | meaning |
|---|---|---|
| `page` | integer | Printed page of the memoir, 129–136 (the same in the CdT volume and the offprint) |
| `line` | integer | Data row on that page, counted from 1. Headings and rules are not counted |
| `year` | integer, or blank | Année. Blank = same as the row above |
| `month` | text, or blank | Mois, as printed: `Janv.` `Fév.` `Mars` `Avril` `Mai` `Juin` `Juill.` `Août` `Sept.` `Oct.` `Nov.` `Déc.`. Two variants are kept as printed: `Juin.` (1 row) and `Juill` (1 row). Blank = same as above |
| `day` | integer | Jour: astronomical day (begins at mean noon of the civil date) |
| `mean_time` | text `h m s` | Temps moyen: Paris mean time, astronomical reckoning, to 1 s |
| `ra_obs` | text | Ascension droite observée: apparent RA. Ancient rows `° ′ ″` (0.1″), modern rows `h m s` (0.01 s) |
| `ra_calc_sec` | text (decimal) | Secondes de l'AR calculée: **seconds only** of the computed RA, unsigned. ″ of arc (ancient) or s of time (modern). The computed degrees/hours and minutes are not printed. They equal the observed ones, or differ by a minute or two where the C − O crosses a minute boundary (ancient C − O reach −66.9″) |
| `ra_c_minus_o` | text (signed decimal) | AR calculée moins AR observée: **C − O** in RA. ″ of arc (ancient) or s of time (modern) |
| `dec_obs` | text `±d m s` | Déclinaison observée: apparent Dec, to 0.1″. `»` = no Dec observed |
| `dec_calc_sec` | text (decimal) | Secondes de la déclinaison calculée: **seconds only** of the computed Dec, ″. A printed sign repeats the Dec's sign; it is not the sign of the seconds. Modern southern rows print `-` (155 rows). p. 136 l. 24–30 (1844–45, northern) print `+`. The other rows, including the three ancient southern rows, are unsigned. Use the magnitude: computed Dec = sign(`dec_obs`) × (its degrees, its minutes ± a minute, these seconds). `»` = none |
| `dec_c_minus_o` | text (signed decimal) | Décl. calculée moins décl. observée: **C − O** in Dec, ″. A difference of signed declinations: for observed −17° 6′ 13.3″ and computed seconds 1.2, computed is −17° 6′ 1.2″ and C − O = +12.1″. `»` = none |
| `lon_c_minus_o` | text (signed decimal) | Longitude calculée moins longitude observée: **C − O** in geocentric ecliptic longitude, ″. Derived by Le Verrier from the RA and Dec errors (formulas A and B, pp. 127–128), not measured. RA-only rows take it from the RA alone |
| `lat_c_minus_o` | text (signed decimal) | Latitude calculée moins latitude observée: **C − O** in geocentric latitude, ″, also derived. `»` = none |
| `notes` | text | The transcribers' notes, unedited, prefixed `p1:` or `p2:` by pass. `[comma]` stands for a printed comma. `low_conf:` marks a digit read with low confidence (33 rows). `check_fail_dec` marks a failed arithmetic check. These notes are working notes, not data, and they record what the transcriber thought at the time |

## `leverrier_1846_misprints.csv`

Cells that are wrong in the print itself, in all three editions consulted. The transcription
keeps the printed value. This file records the correction. 7 entries.

| column | type | meaning |
|---|---|---|
| `page` | integer | Printed page of the cell |
| `line` | integer | Row on that page |
| `field` | text | Table 77 column of the cell (a column name of `leverrier_1846_table77.csv`) |
| `printed` | text | The cell as printed. Equal to the value in `leverrier_1846_table77.csv` (the validation checks this) |
| `corrected` | text | The corrected value, in the same format. Blank for `log_only` |
| `action` | `apply` or `log_only` | `apply`: the derived file uses `corrected` (4 entries: two Dec degree digits, two mean times). `log_only`: recorded but not corrected, either because the right value cannot be determined or because the column is not used downstream (3 entries) |
| `reason` | text | The argument. It is built from the table itself: neighbouring rows, daily motion, the arithmetic of the row, ecliptic geometry, and the transit identity for times. Each entry also says which witnesses print the value |
| `uses_ephemeris` | `yes` or `no` | Whether the argument needs a modern ephemeris of Uranus. `no` for every entry. The 1771 entry was *found* by comparing with JPL, but its argument does not use the ephemeris |

## `leverrier_1846_resolution_log.csv`

Every cell where the two independent passes disagreed: 14 cells out of 279 × 12. Each was
resolved from the third witness.

| column | type | meaning |
|---|---|---|
| `page` | integer | Printed page |
| `line` | integer | Row on that page |
| `field` | text | Table 77 column |
| `pass1` | text | Pass 1's reading (archive.org scan). `?` = a digit the transcriber could not read |
| `pass2` | text | Pass 2's reading (Lille scan), same convention |
| `resolved` | text | The value adopted in `leverrier_1846_table77.csv` |
| `witness` | text | The file and **PDF** page of the third witness consulted (Gallica CdT 1849; PDF page = printed page + 2) |
| `reason` | text | Why the value was adopted. Usually "Gallica prints … clearly". One entry (p. 135 l. 4, year 1835 vs 1833) is resolved by argument, because the glyph is ambiguous in all three scans |

## `leverrier_1846_pass1.csv`, `leverrier_1846_pass2.csv`

The two independent transcriptions before reconciliation. Columns as in
`leverrier_1846_table77.csv`, except:
- A cell can contain `?` for an unreadable digit.
- `notes` is that pass's own note, without the `p1:`/`p2:` prefix.

`leverrier_1846_table77.csv` is exactly the merge of these two files with the resolution log.
In the project repository a test regenerates it and compares.

Columns of both files:

| column | type | meaning |
|---|---|---|
| `page` | integer | As in table77 |
| `line` | integer | As in table77 |
| `year` | text | As in table77; may contain `?` |
| `month` | text | As in table77 |
| `day` | text | As in table77; may contain `?` |
| `mean_time` | text | As in table77; may contain `?` |
| `ra_obs` | text | As in table77; may contain `?` |
| `ra_calc_sec` | text | As in table77; may contain `?` |
| `ra_c_minus_o` | text | As in table77; may contain `?` |
| `dec_obs` | text | As in table77; may contain `?` |
| `dec_calc_sec` | text | As in table77; may contain `?` |
| `dec_c_minus_o` | text | As in table77; may contain `?` |
| `lon_c_minus_o` | text | As in table77; may contain `?` |
| `lat_c_minus_o` | text | As in table77; may contain `?` |
| `notes` | text | This pass's notes |

## `spotcheck_rows.csv`

The spot-check sample: 28 rows. All 14 resolution-log rows and the 2 misprint rows that
existed at review time; then rows to cover at least 2 ancient rows and every page; then random
rows, with seed 20261106.

| column | type | meaning |
|---|---|---|
| `page` | integer | Printed page |
| `line` | integer | Row on that page |
| `reason` | text | Why the row was drawn: `resolution log: <field>`, `misprint table: <field>`, `stratum: …` or `random` |
| `seed` | integer | The sampling seed (20261106) |

## `spotcheck_results.json`

The spot-check record. It is JSON, not a table.

| key | meaning |
|---|---|
| `seed` | Sampling seed |
| `table77_sha256` | SHA-256 of the `leverrier_1846_table77.csv` that was checked. It equals the file in this package |
| `exported` | Time of the export (UTC) |
| `reviewer` | Who checked: "automated spot check, against Gallica scan only; not an independent human check" |
| `all_match_rows` | Rows confirmed with the review page's "all match" shortcut. Empty: every cell was ticked individually |
| `rows[]` | One object per sampled row: `page`, `line`, `reason`, `flagged` (drawn because of the log or the misprint table), `used_all_match`, `cells` (for each of the 12 data columns, the `value` checked and `matches`), and `notes` |

Result: 336 of 336 cells match the Gallica image.

## `spotcheck_human_results.json`

The human spot-check record. The dataset was spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies).
- It covers the same 28 rows as `spotcheck_results.json`.
- The reviewer was the dataset author, on 2026-11-07.
- Every cell was ticked individually; the "all match" shortcut was not used.
- The review page did not show the automated results, and its ticks were stored apart from that
  review's, so nothing was pre-ticked.

The keys are those of `spotcheck_results.json`, except `reviewer`:

| key | meaning |
|---|---|
| `reviewer` | `"human"` |
| `reviewer_name` | The name the reviewer entered on the review page |
| `seed`, `table77_sha256`, `exported`, `all_match_rows`, `rows[]` | As in `spotcheck_results.json`. `table77_sha256` equals the table in this package |

## `leverrier_1846.csv` (derived)

A convenience file for analysis, generated from `leverrier_1846_table77.csv` and the `apply`
entries of the misprint table by `fork.data.observations` in the project repository. It
contains **one row per measured quantity**: 279 RA rows and 250 Dec rows, 529 in all. The 29
RA-only observations have no Dec row.

Unlike the table, this file contains **inferences**: the epochs (time meridian, ΔT), the
observatories (the table has no such column), and the transit residuals. Where a value is not
printed, the column says how it was obtained. If you need only what Le Verrier printed, use
`leverrier_1846_table77.csv`.

| column | type | meaning |
|---|---|---|
| `id` | text | Unique: `<obs_id>-ra` or `<obs_id>-dec` |
| `obs_id` | text | The Table 77 row: `LV1846-T77-p<page>-l<line, 2 digits>`. The RA and Dec rows of one observation share it |
| `epoch` | text | Date and time as printed, with the year and month filled in: `1690 Déc. 23, 9 41 25 t.m.` (t.m. = temps moyen). It shows the **printed** time, before any misprint correction |
| `epoch_tt_jd` | number | `epoch_ut_jd` + ΔT, as a Julian Date in TT. ΔT comes from the Espenak & Meeus (2006) polynomials at mid-month, 6–14 s over 1690–1845 |
| `epoch_ut_jd` | number | Julian Date (UT): JD of civil noon on the printed date + mean time − Paris longitude (9 m 20.9 s). It uses the **corrected** time where a misprint is applied |
| `time_meridian` | text | `Paris` for every row (see Time, above) |
| `time_reckoning` | text | `astronomical (day begins at mean noon)` |
| `observer` | text | From §74 (pp. 124–126) for the ancient rows: Flamsteed (1690, 1712, 1715), Le Monnier (1750, 1764–1771), Bradley (1753), Mayer (1756). Blank for modern rows: the memoir does not name them |
| `observatory` | text | `Greenwich`, `Paris` or `Göttingen`. Ancient rows from §74. Modern rows from the date ranges on p. 127: Greenwich for 1781–1800, January 1807–1809 and 1828–1830; Paris for 1801–1827 and 1835–1845. The two 1828 rows are inferred (see next column) |
| `observatory_basis` | text | `stated`, or the inference. For the two 1828 rows (p. 134 l. 27–28): p. 127 lists 1828 under both observatories. The transit identity fits Greenwich with Paris mean time to about 2 s, and misses Paris by about 563 s |
| `instrument` | text | Always blank: not stated |
| `quantity` | `ra` or `dec` | The measured coordinate |
| `value` | number | The observed value in **degrees**, with the `apply` corrections. RA: ancient rows as printed in arc; modern rows 15 × hours. Dec: signed, with the sign taken from the string |
| `unit` | text | `deg` |
| `value_as_printed` | text | The `ra_obs` or `dec_obs` cell as printed, **before** correction |
| `unit_as_printed` | text | `deg arcmin arcsec` or `h m s` |
| `frame_as_printed` | text | `apparent geocentric; true equator and equinox of date (nutation, aberration)` |
| `equinox_ref` | text | `of_date` |
| `equinox_epoch` | text | Blank: not applicable to `of_date` |
| `equinox_citation` | text | Where the memoir states the frame (p. 111 §72; p. 127) |
| `precession_constant_arcsec_per_yr` | text | `not_stated`. §72 gives no precession constant, and none is needed for of-date places |
| `precession_rate_arcsec_per_yr2` | text | Blank |
| `precession_citation` | text | Citation for the two precession fields |
| `lv_residual` | number | **Le Verrier's O − C**, the negative of the printed C − O, in arcseconds of the coordinate itself. RA: ancient rows −(C − O) in ″; modern rows −15 × (C − O) with C − O in s. **Not** multiplied by cos δ. Dec: −(C − O) in ″. It is computed from the printed C − O, so for 1764 Jan 15 (Dec, a `log_only` misprint) it carries the inconsistent printed sign |
| `lv_residual_unit` | text | `arcsec of RA (O - C)` or `arcsec (O - C)` |
| `lv_c_minus_o_as_printed` | text | The printed `ra_c_minus_o` or `dec_c_minus_o` cell (C − O, in s or ″ as printed) |
| `lv_weight` | text | Always blank. Table 77 prints no weights; the only weights are the observation counts of §78 and §131 |
| `lv_dlon_c_minus_o_as_printed` | text | The printed `lon_c_minus_o`, **C − O**, ″. The same on the RA and Dec rows of an observation |
| `lv_dlat_c_minus_o_as_printed` | text | The printed `lat_c_minus_o`, **C − O**, ″; `»` = none |
| `transit_resid_s` | number | Local apparent sidereal time at the observatory (ERFA `gst06a`, IAU 2006/2000A) minus the observed RA, in seconds of time, wrapped to ±12 h, at the corrected epoch. A meridian observation timed correctly gives about 0. Modern rows: median \|residual\| 0.5 s, maximum 2.4 s. Ancient rows: median 1.4 s, maximum 3.6 s, except 1768 Dec 30 (+9.5 s, a `log_only` misprint). With Greenwich mean time instead, the Greenwich rows would miss by about 563 s. This column is the evidence for the time meridian |
| `era` | text | `ancient` (1690–1771) or `modern` (1781–1845) |
| `source_ref` | text | Citations. `CdT1849Add:p129:row1 (Table 77)` = *Connaissance des Temps pour 1849*, Additions, printed p. 129, row 1. The text citation for the attribution follows |
| `misprints_applied` | text | Corrections applied to this observation: `<field>: <printed> -> <corrected>`, separated by `; ` |
| `notes` | text | The source of the attribution, log-only misprints for the row, and for 1828 the transit residuals for all four site/time combinations |
