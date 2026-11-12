# Adams's memoir, §§ 8–55: reconciliation of the two transcription passes

Date: 2026-11-11. Source: the pamphlet printing of Adams's *Explanation of the observed
irregularities in the motion of Uranus* (Nautical Almanac appendix, 1846; cache
`adams_1846_nautical_almanac_appendix_archive_org.pdf`, hash in `MANIFEST.md`), PDF pages 12–37,
printed pages 6–31 (printed = PDF − 6). Pass 1 was read from the archive.org renders
`pamphlet_p012.png`–`pamphlet_p037.png` except PDF 25–27, which archive.org serves at 364–367 ×
592–599 px against 546–587 × 888–958 px for the other pages; for those three pages pass 1 read
220 dpi renders in two halves made from the PDF's embedded scans (`renders/pamphlet_hires/p025a.png`
… `p027b.png`; the same was made for PDF 18, 24 and 28–37, used below). Pass 2 was an independent
reading of the 26 archive.org renders by a reader given only the renders and the column
specification, who upscaled slices of each page and marked every digit it could not decide with
"?". Files: `adams_memoir_passes/pass1/` and `pass2/` (six CSVs each, as delivered; pass 1's
normal places, series, equations, solution carry the printed values, its `statements.csv` the
blocks, its `notes.csv` the spec deviations; pass 2 likewise). Comparison by
`scripts/p1b_adams_memoir_compare.py`: every value keyed as the specification names it
(normal places by year and column; series and solution by hypothesis, section, item; equations by
hypothesis, section, equation, term), a leading "+" dropped on both sides, nothing else normalised.

| file | keys pass 1 | keys pass 2 | common | disagreements |
|---|---|---|---|---|
| normal places | 352 | 352 | 352 | 0 |
| series | 139 | 140 | 139 | 2 |
| equations | 1407 | 1407 | 1407 | 12 |
| solution | 91 | 104 | 90 | 1 |
| total | | | 1988 | 15 |

Of the 15: 12 are pass 2's "?" marks, all on PDF 25–27, where it could not separate 1 from 4 or
3 from 5 from 8 on the archive.org render; 3 are readings that differ. Each was decided on the
image named, at 220 dpi, and where that did not settle it, on a 600 dpi crop of the line.

## The twelve "?" tokens (PDF 25–27)

All settled on the 220 dpi renders, every one to pass 1's digit; the 600 dpi crops were not needed.

| key | pass 1 | pass 2 | decided | image |
|---|---|---|---|---|
| series 2/33 log_adb0 | 9.57333 | 9.?7333 | 9.57333 | `p025a.png` |
| series 2/33 log_adb3 | 9.50854 | 9.5?854 | 9.50854 | `p025a.png` |
| equations 2/35 (h₂) p₂ | 1.6727 | 1.67?7 | 1.6727 | `p026a.png` |
| equations 2/35 (h₃) δx₂ | 10.5558 | 10.55?8 | 10.5558 | `p026a.png` |
| equations 2/35 (h₃) p₃ | 3.5852 | 3.58?2 | 3.5852 | `p026a.png` |
| equations 2/36 (k₂) k₃ | 5.6294 | 5.629? | 5.6294 | `p026b.png` |
| equations 2/36 (k₃) k₂ | 5.6294 | 5.629? | 5.6294 | `p026b.png` (the two prints agree at 220 dpi) |
| equations 2/37 comb1 δx₁ | 6.0448 | 6.04?8 | 6.0448 | `p027a.png` |
| equations 2/37 comb1 h₂ | 13.6704 | 13.670? | 13.6704 | `p027a.png` |
| equations 2/37 comb1 h₃ | 19.2974 | 19.297? | 19.2974 | `p027a.png` |
| equations 2/37 comb1 p₃ | 2.4144 | 2.414? | 2.4144 | `p027a.png` |
| equations 2/37 comb2 k₂ | 7.3294 | 7.329? | 7.3294 | `p027a.png` |

## The three differing readings

| key | pass 1 | pass 2 | decided | how |
|---|---|---|---|---|
| equations 1/18 (k₁) q₃ (PDF 18) | 11.3316 | 11.3346 | **11.3346** (pass 2) | 600 dpi crop `z018_k1_q3b.png`: the fourth decimal glyph is a 4 with a thin stem, which reads as 1 at 220 dpi (`p018a.png`); the symmetric coefficient, (q₃) k₁ on the same page, prints 11.3346 (`z018_q3_k3b.png`), both passes agree on it, and a normal equation's matrix is symmetric by construction. |
| equations 2/35 (x) h₁ (PDF 26) | 7.3540 | 7.3510 | **7.3540** (pass 1) | `p026a.png`: 7.3540; the symmetric coefficient (h₁) δx₁ prints 7.3540 three lines below. |
| solution 1/30 p₃/m′ formula, first coefficient (PDF 24) | 48.55 | 18.55 | **48.55**, with the print noted | 600 dpi crop `z024_formula.png`: the glyph is a damaged 4 (its left stroke missing) and reads "18"; pass 2 recorded it as printed and said so. The q₃/m′ line directly below prints 48.55 with the same amplitude by the form of the expression, and the arithmetic settles it: with § 30's θ = −51° 30′, β = 50° 15′.8, e′ = 3.2206, β′ = 262° 28′ the pair gives p₃/m′ = 271.46, q₃/m′ = −207.38 against the printed 271.57, −207.24 for 48.55, and 258.9, −180.1 for 18.55. |

After these decisions the two passes agree on every one of the 1988 common values.

## Keys present in one pass only (no disagreement)

Pass 2 recorded 14 printed items the specification did not list: α = 0.5 under § 12; the epoch
1810.328 and the date 1846.762 under §§ 30 and 47; the q₃/m′ formula coefficients under §§ 30
(48.55, −93.01) and 47 (33.93, −63.41); α = sin 31° = 0.515 under § 32; α restated under §§ 49 and
50. All are kept in the final files (their values match pass 1's reading of the same pages where
pass 1 recorded the fact elsewhere). The § 53 statement α = sin 35° = 0.574 is one fact recorded
under hypothesis 0 by pass 1 and hypothesis 1 by pass 2; kept as hypothesis 0 (it belongs to
neither).

## Where the specification did not fit the page (both passes, merged)

1. §§ 23, 40: δε is printed with no numeral in the ancient equations; recorded as `de` = 1.
2. §§ 29, 46: the four bracketed terms of the logarithmic equations are the combinations
   p₃ cos θ − q₃ sin θ, p₃ sin θ + q₃ cos θ, p₃ cos 2θ − q₃ sin 2θ, p₃ sin 2θ + q₃ cos 2θ (all in
   p₃, q₃); the specification's names P1c-Q1s … P2s+Q2c are kept in order of appearance. The first
   term of L2 and L3 is printed without a sign (hyp. I: (1.65083), (1.79213); hyp. II: (1.65190),
   (1.91407)); recorded without a sign.
3. §§ 28, 45: h₃/m′ and k₃/m′ carry the terms 0.007460 p₃/m′ + 0.008974 q₃/m′ and
   −0.008974 p₃/m′ + 0.007460 q₃/m′; the printed layout interleaves h₁, h₂ with k₁, k₂ on shared
   lines; rows in specification order.
4. § 25, 41: the p₃ and q₃ multipliers are printed unsigned; the fourth q₃ multiplier (1753)
   prints "0 4" with the raised dot faint; recorded 0.4.
5. § 34 (hyp. II) prints no epoch, unit or anomaly rate, only the mean motions and arguments.
6. §§ 37, 38: the combined equations carry no printed label; recorded as comb1, comb2.
7. § 39: rel2 and rel4 have δy₁, δy₂ on the left (`lhs_dy1`, `lhs_dy2`); mult2's first
   multiplier (0.00165) is printed unsigned.
8. § 41: the combined equations still contain δε (2.80, 3.60) and δn (−27.82, −91.84), unlike
   § 25; recorded under `de`, `dn` (19 terms each); § 42 eliminates them.
9. § 43: the relation lines have four terms on the left (δx₁, δy₁, δx₂, δy₂); recorded as
   `lhs_dx1, lhs_dy1, lhs_dx2, lhs_dy2` plus `rhs_dx1, rhs_dy1`.
10. § 54 contains no numbers; the α = sin 35° statement is the last paragraph of § 53 (PDF 35).
11. § 55: the Flamsteed errors print +44″·5 and +50″·0; recorded +44.5, +50.0.
12. §§ 56–61 (radius-vector corrections, δζ tables) are outside the step and not transcribed
    (statements S54–S56).

Pass 2 also noted that § 18 (k₁)'s δn coefficient 88.4630 has a slightly less crisp 4, and that
§ 34's n′ = 4° 48′.5 could be read 48′.3 at archive.org's resolution; pass 1 read both at 220 dpi
(88.4630, 48′.5; `p018a.png`, `p025b.png`) and the two passes agree on them.

## What stage B inherits

Stage B's checks run on the final files. Two of the three decided readings were settled partly by
the symmetry of the normal-equation matrix; B1–B3 test that symmetry on every pair and will show
whether any other printed pair disagrees (none was found by either pass). The § 30 formula
coefficient is used as 48.55; the statement row S28 and the note above carry the print.
