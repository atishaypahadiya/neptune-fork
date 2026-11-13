# Brookes 1970, on the prediction of Neptune

Reading note, 13 November 2026, written from the notes already in the repository and checked
against four pages of the cached scan.

## Citation

C. J. Brookes, "On the prediction of Neptune", *Celestial Mechanics* 3 (1), 67–80 (1970).
DOI 10.1007/BF01230433. ADS bibcode 1970CeMec...3...67B. Received 5 February 1970; University
of Aston in Birmingham.

**Verified:** Crossref's record for the DOI (`data/cache/prior_work/crossref_brookes1970.json`,
SHA-256 `8ca9ae1802fdd8ed8802c4840f0f5b333838af3e2b8680b89b47a9a5b3f57092`) returns author
"Brookes, C. J.", that title, *Celestial Mechanics* volume 3, issue 1, pages 67–80, issued
1970-03. The scan's first page prints "Celestial Mechanics 3 (1970) 67-80" and "(Received 5
February, 1970)". Kollerstrom 2006 (Part 1, note 126) dates the issue "December 1970"; Crossref
says March. I log the difference and do not resolve it.

## Access

No new download. The scan is the one obtained on 2026-11-11:
`data/cache/adams_challis/ads_1970CeMec...3...67B.pdf`, from
`https://ui.adsabs.harvard.edu/link_gateway/1970CeMec...3...67B/ADS_PDF`, 14 pages, no text
layer, SHA-256 `054b6285a75b0b82803d405bb88e520c5d1a7a7aabcc25511f4f084af5d8da9e` (I recomputed
the hash on 2026-11-13; it matches `data/historical/adams_challis/MANIFEST.md`). Publisher's
copyright notice on every page; cache only.

## What was read

For this note I read pp. 67, 78, 79 and 80 (PDF pages 1, 12, 13, 14) as OCR text made from
200 dpi renders, with the tables' digits as the OCR gives them. Everything else below is taken
from the repository's notes, which record a full reading from renders on 2026-11-11:
`data/historical/adams_challis/MANIFEST.md` (the summary under the Brookes row),
`docs/notes/P1B_adams_reproduction.md` (§ 0, R5, R10) and `docs/notes/P1B_writeup.md`. Claims
from those notes are tagged with the note, not with a page I read.

## Claims

1. The abstract: an investigation of "the theoretical procedures involved in the theory of John
   Couch Adams, in order to determine the effect of variations in the time and mean distance on
   the accuracy of the predictions"; the work is part of a London Ph.D. thesis. [read: p. 67]
2. Table I gives Adams 37.25 AU, e 0.1206, perihelion 299° 11′, Sun/planet 6666, true longitude
   at discovery 329° 27′; Le Verrier 36.15, 0.1076, 284° 45′, 9300, 326° 0′; Neptune 30.07,
   0.0086, 44°, 19300, 326° 57′. [read: p. 67, OCR]
3. From the memoir's equations he re-derives Hypothesis I as e′ 0.1509, m′ 0.8267, ϖ′ 311° 35′,
   mean longitude 323° 39′ at 1846.762, true longitude at discovery 328° 4′ against Adams's
   328° 42′, and attributes the difference "mainly to a small numerical error occurring in one
   of the equations derived by Adams", without naming the equation. [not read for this note:
   `P1B_adams_reproduction.md` § 0 and the MANIFEST summary]
4. Hypothesis II by the first hypothesis's procedure comes out 0.1387, 0.7393, 303° 52′,
   324° 47′, true 331° 31′. [not read for this note: same two notes]
5. He varies the number of equations (61, 31, 21), the inclusion of Flamsteed 1690, a
   first-order-only variant, the epoch and α. [not read for this note: MANIFEST summary]
6. Table X, solutions for 1845/6 at epoch 1846.762, for α = sin 29°, 30°, 31°, 32°, 34°, 35°,
   37°: a′ 39.56, 38.38, 37.25, 36.19, 34.30, 33.40, 31.86; true longitude 323°.8, 327°.1,
   330°.4, 333°.2, 317°.6, 319°.1, 321°.8 against Neptune's 327°.0. Tables XI–XIV repeat this
   for 1840, 1835, 1830 and 1825. [read: pp. 78–79, OCR; the first four longitudes agree with
   those quoted in `P1B_adams_reproduction.md` § 0]
7. "The position of the planet Neptune could have been predicted with reasonable accuracy,
   using any value of a′ between 30 AU and 40 AU at any time after 1830." [read: p. 78]
8. The solutions "are particularly sensitive to any variation in the set of values representing
   the discrepancies in mean longitude", while the true longitude changes by "a little over
   1½°"; "there is little point in regarding the solution of Le Verrier as more accurate than
   that of Adams". [read: p. 79, OCR, the fraction checked on the page image]
9. "There is nevertheless an element of chance present": without the second-order terms no
   sensible solution is obtained at α = sin 34°, and with a fourth principal inequality none at
   sin 30°. [read: p. 79]
10. He thanks W. M. Smart for suggesting the problem and St John's College Library for copies
    of Adams's manuscripts; the references include his 1969 thesis, Lyttleton 1960 and Smart
    1947. [read: p. 80]
11. Kollerstrom reads Brookes p. 71 as saying "the first hypothesis was the more accurate of
    the two" and p. 78 as showing how the solutions vary with orbit radius, "317° at 34 AU".
    [not read on Brookes's p. 71: Kollerstrom 2006, Part 1, notes 126 and 128, OCR; the 317°.6
    at sin 34° is in claim 6]

## What it means for this repository

The repository's own account is already specific: the reproduction "agrees with" Brookes
(`README.md`), lands within 0.02° and 0.01° of him (`P1B_writeup.md`), and its novelty is
limited to locating the numerical difference he did not place, the χ² profile over α with its
pole, and the split between data and method (`P1B_adams_reproduction.md` § 0 and R10). What I
read supports that limit and adds one thing: claims 7 to 9 are a sensitivity study of Adams's
solution to distance, epoch and the terms kept, with the words "element of chance", so the
question of how much of his nearness was fortune was asked in print in 1970, by point solutions
and not by a noise study.

One disagreement inside the repository, already superseded: `docs/notes/P1B_sources.md` § 5
says "Paywalled; not obtained", written before the ADS scan was found on 2026-11-11; the later
notes are right.

## Needs a browser / library

- Brookes's 1969 London Ph.D. thesis, *The Discovery of Neptune: A Critical Examination of the
  Theory of Adams*, which is where the equation with the numerical error is said to be
  identified. Needs a library (University of London theses).
