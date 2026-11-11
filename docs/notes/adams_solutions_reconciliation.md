# Adams's solutions: reconciliation of the two transcription passes

Date: 2026-11-10. Files: `adams_solutions_pass1.csv` (35 statements, read directly from the page
renders), `adams_solutions_pass2.csv` (44 statements, an independent reading of the same renders under
a written brief, with the pamphlet pages corrected twice during the pass because the first brief
named the wrong pages). Comparison: every numeric token in the element columns and in `other_values`
was extracted from both files and compared statement by statement (34 matched pairs; pass 2 has ten
further rows for residual tables and the Sampson narrative that pass 1 did not transcribe).

## Result

**No numeric disagreement between the passes.** Every differing cell was a difference of
presentation (raised dot against decimal point, `θ=` against `theta=`, ISO dates against the dates as
printed, numbers placed in `other_values` by one pass and in a named column by the other). Pass 2
reported no cell marked uncertain after the higher-resolution renders.

## Differences between printings of the same solution (not transcription errors)

Each was checked on the page image named; both passes read the same digits.

| solution | quantity | printing A | printing B | decided by |
|---|---|---|---|---|
| S5 | longitude of perihelion | 315° 55′ (Airy No. 11, MNRAS 7, 129; Memoirs 16, 395) | 315° 57′ (Airy No. 20; abstract; memoir § 31; Sampson table V) | `mnras7_airy_p009.png`, `memoirs16_p430.png`, `mnras7_airy_p018.png`, `pamphlet_p024.png`: both as printed; the October 1845 paper and the final value differ by 2′ (Adams's "slightly corrected", memoir § 4) |
| S5 | mass | 0.0001656 (No. 11; abstract; memoir) | 0.00016563 (No. 20) | `mnras7_airy_p018.png`, `memoirs16_p441.png`: the letter gives one more digit |
| S5 | eccentricity | 0.1610 (No. 11) | 0.16103 (elsewhere) | `mnras7_airy_p009.png`: four figures in the 1845 paper |
| S5 | mean longitude 1846 | 325° 8′ at 1 Oct 1846 (No. 20) | 325° 7′ at 6 Oct 1846 (memoir § 31; abstract) | `mnras7_airy_p018.png`, `pamphlet_p024.png`: as printed. Five days of mean motion is 1′, so 325° 7′ on 6 Oct implies 325° 6′ on 1 Oct, 2′ below the letter |
| S5 | mean longitude at epoch | 269° 25′ (memoir § 30) | 269° 25′.5 (Sampson table V) | `pamphlet_p024.png`, `sampson_p009.png`: Sampson keeps the half minute |
| S6 | eccentricity | 0.12062 (No. 20) | 0.120615 (memoir § 48; abstract "0·12 0615"; Sampson) | `mnras7_airy_p018.png`, `pamphlet_p031.png`, `mnras7_adams_p003.png`: the abstract has a space inside the number |
| S6 | mean longitude 1846 | 323° 2′ at 1 Oct 1846 (No. 20) | 323° 2′ at 6 Oct 1846 (memoir § 48) | same figure printed for both dates |
| S6 | θ | −46° 55′ (memoir § 47) | −46° 54′ 50″ (Sampson E.VIII) | `pamphlet_p031.png`, `sampson_p026.png` |
| S6 | perihelion 1810 | 298° 41′ (memoir § 47) | 298° 40′ 40″ (Sampson E.VIII) | same pages |
| S6 | assumed mean distance | a/a′ = 0.515 = sin 31° (memoir § 32); 1.942 × Uranus (abstract) | 37.25 (Sampson table VI) | `pamphlet_p024.png`, `mnras7_adams_p003.png`, `sampson_p009.png`: 1/0.515 = 1.9417; 1.942 × 19.18 = 37.25 |
| S4 | mass | 0.00017302 (Sampson table IV) | 0.000173 (sheet of notepaper, Sampson p. 166) | `sampson_p009.png`, `sampson_p024.png` |
| S4 | perihelion | 320° 30′.2 (table) | 320° 30′ (sheet) | same pages |
| S7 | distance ratio | a/a′ = 0.57 (No. 20) | sin 35° = 0.574 (memoir § 53; abstract) | `mnras7_airy_p019.png`, `pamphlet_p035.png` |

## Notes on the Sampson readings

- The table on Sampson p. 151 (`sampson_p009.png`) gives the six determinations at epoch 1810.328
  with mean distance, eccentricity, longitude of perihelion, mean longitude at epoch and mass; the
  eccentricities of I and II are assumed zero, and III–VI are the full-weight solutions.
- The 1843 result (S1) is read on the facsimile plate (`sampson_p038.png`, manuscript p. 10, lower
  right: "253·4 Long. Planet at oppn. 1810 for 33 yrs 50·6 304·0 Long. at oppn. of Uranus 1843") and in
  Sampson's catalogue (p. 160, `sampson_p018.png`). Both passes read the plate the same way.
- Sampson's dating of the six: I 1843 (the sheet dated September 30, 1843); II undated; III the note
  of April 28, 1845 and a comparison dated May 19, 1845; IV September 18, 1845 (θ) and the sheet of
  notepaper "about the end of Septr, 1845"; V communicated October 21, 1845 (Sampson p. 168) and
  left at Greenwich "on one of the last days of October" (Airy); VI the note of August 20, 1846 and
  the letter of September 2, 1846.
