# Bouvard's 1821 tables of Uranus: structure and what a transcription would cost

Post hoc scoping, 12 November 2026. No value is transcribed here and nothing is computed
from the tables. The figures quoted below were read once from the page images to describe
the book; none has had a second entry, and none may be used as data.

Source: Bouvard, A., *Tables astronomiques publiées par le Bureau des longitudes de France,
contenant les tables de Jupiter, de Saturne et d'Uranus* (Paris, 1821), the Harvard copy
(`data/historical/adams_challis/MANIFEST.md`: hash, completeness XXVIII + 110 pp.).

Why it matters: Adams's 30 normal places (memoir § 8) are observation minus Bouvard's
tables, with his § 7 corrections, in heliocentric true longitude. With Bouvard's longitude
at those 30 epochs, his data become places that the project's model can fit, which is what
the step 3c question (his shift of about 11°) needs before it can be attributed.

## The introduction's statement of method (pp. xiii to xvii)

- Two series of observations: 17 made before the discovery (one of Mayer, one of Bradley,
  three of Flamsteed, twelve of Lemonnier) and the oppositions and quadratures from 1781
  (Maskelyne and Pond to 1815; Paris 1800 to 1819).
- Perturbations by the formulas of the *Mécanique céleste*, with "les nouvelles masses de
  Jupiter et de Saturne".
- Provisional tables, then 77 equations of condition in four unknowns (corrections of the
  epoch, the mean motion, and the two parts of the equation of the centre):
  (1 + 2e cos φ) x + t (1 + 2e cos φ) y + z sin φ − u cos φ = dV″.
- One solution from both series left errors he would not accept (the table on p. xiv:
  +98″.7 in 1690, −124″.0 in 1750, +76″.2 in 1819, centesimal seconds as printed). He then
  fitted the modern series alone, and left it "aux tems à venir" to show whether the
  failure lies in the old observations "ou si elle dépend de quelque action étrangère et
  inaperçue, qui aurait agi sur la planète".
- The equations of condition are printed (pp. xvi to xvii) with the number of observations
  in each.

## Units and epochs

- The circle has 400 degrees, the degree 100 minutes, the minute 100 seconds; the day has
  10 hours of 100 minutes of 100 seconds (p. xviii). Every angle in the tables and in the
  printed formulas is in this division. Table LIV (p. 42) converts.
- Epochs are for the mean midnight between 31 December and 1 January, Paris.
- Mean elements are for 1800, counted from the fixed equinox of 1800 (p. ii). Precession is
  taken uniform at 154″.63 a year (centesimal); the correction is in Jupiter's Table III
  (p. 12), which the Uranus tables borrow (p. xx).

## Masses adopted (p. ii, Sun = 1)

| planet | as printed | remark |
|---|---|---|
| Jupiter | 1/1070,5 | |
| Saturn | 1/5512 | misprint candidate: Bouvard's mass of Saturn is usually quoted as 1/3512. Not settled here; a second source is needed before either is used |
| Uranus | 1/17918 | |

The masses come out of the Jupiter and Saturn solution; the Uranus theory takes those of
Jupiter and Saturn from it.

## The theory in closed form (pp. ii, xv, xvi)

The introduction prints the whole theory that the tables tabulate:

- mean longitudes of Jupiter, Saturn and Uranus (epoch value and annual motion), the
  perihelion and node of each, and the three eccentricities (p. ii);
- the heliocentric longitude V″ of Uranus in twelve groups, numbered I to XII as the
  tables' arguments are (p. xv): group I the equation of the centre, five terms with
  secular variations; II five terms in the Saturn minus Uranus argument; III the great
  inequality, argument (φ′ − 2φ″), with a secular term in coefficient and phase, and one
  term in (2φ′ − 4φ″); IV and V Jupiter terms; VI to XII one term each. About 40 printed
  numbers;
- the reduction to the ecliptic, one term;
- the radius vector: a constant with a secular term, four elliptic terms and seven
  perturbation terms (pp. xv to xvi);
- the polar distance: the elliptic term and three perturbation terms (p. xvi);
- the aberration of Uranus in longitude, three terms (p. xx).

## The tables of Uranus (pp. 83 to 110, and p. 41)

| table | pages | content | argument |
|---|---|---|---|
| I | 84 to 89 | epochs of the mean motions, 1780 to 1900, one row a year: mean longitude, perihelion, node, arguments II to XII | year |
| II | 90 | motions for other centuries, to add to Table I | century |
| III | 91 | motions for the months | month |
| IV | 92 | motions for the days | day |
| V to VIII | 93 | motions for hours, minutes, seconds; proportional parts; decimals of the year | |
| IX | 94 | the great inequality, every tenth year 1680 to 1900 (23 rows), with the corrections of arguments II to XII | year |
| X | 95 to 98 | equation of the centre for 1800 with its secular variation, every degree of 400, with first and second differences | I, mean anomaly |
| XI to XIV | 99 to 100 | equations of the longitude, "toujours additives" (a constant is added to each) | II to V; the circle in 10,000 or 1,000 parts |
| XV to XXI | 101 to 102 | the smaller equations of the longitude | VI to XII; the circle in 1,000 parts |
| XXII | 103 to 104 | radius vector for 1800 with secular variation | I |
| XXIII to XXVII | 105 to 106 | perturbations of the radius vector | II, III, IV, VI, VII−VI |
| XXVIII | 107 to 108 | distance from the north pole of the ecliptic for 1800 | XIII, longitude in the orbit minus the node |
| XXIX to XXXI | 109 | perturbations of the polar distance | III, VII, IX |
| XXXII | 110 | reduction to the ecliptic and log cosine of the latitude | polar distance |
| LII, LIII | 41 | aberration of Uranus in longitude | H − G and G |

The errata (p. xxviii) list one correction to the Uranus pages: p. 41, the signs opposite
175° in the aberration table.

Not established: whether the tables equal the printed formula to their last figure (the
constants added to make every equation additive, the secular terms, Table IX against group
III); whether the copy Adams or Greenwich used carried later corrections; the introduction
gives worked examples for Jupiter and Saturn only, none for Uranus.

## What a transcription would cost

Table 77 as the unit: 279 observations, about 3,350 cells at the 12 cells a row of its spot
check, entered twice. Cell counts below are estimates from the page layouts, not counts.

| route | what is entered | cells, one entry | in Table 77 units | what it gives |
|---|---|---|---|---|
| A, the printed formula | the mean elements, the twelve groups of V″, the reduction, the precession; the radius vector and polar distance only if wanted | about 60 for the longitude, about 100 with the rest | 0.02 to 0.03 | Bouvard's longitude at any epoch, by a short function with tests |
| B, the table rows the 30 epochs touch | Table I rows for the years (the nine old dates through Table II), months and days, Table IX, two rows of Table X and of each of XI to XXI per epoch, with differences | about 1,800 to 2,000 | 0.5 to 0.6 | the 30 places exactly as a computer of 1821 to 1846 got them |
| C, every Uranus longitude table | Tables I to XXI and LII, LIII in full | about 6,500 to 7,500 | about 2 | all of the above at any epoch, and a full check of A |
| D, all of pp. 83 to 110 | C plus radius vector, polar distance and reduction | about 10,000 | about 3 | the complete Uranus tables |

Proposed order, for a ruling: A first, then B at 5 of the 30 epochs (about 330 cells) as
the test of whether the formula reproduces the tables; B in full only if they differ by
more than the 0″.1 to which Adams printed his old places. C and D are not needed for the
question.

Also needed, and not in this book:

- the instant of each of the 30 places. The 21 modern ones are at mean oppositions, three
  synodic periods apart from 1810.328 (memoir § 13). The nine old ones need the dates of
  the observations (Table 77 has them);
- Adams's § 7 corrections to Bouvard, as printed in the memoir;
- conversion from the centesimal division, Paris midnight and the equinox of 1800 with
  Bouvard's precession to the project's frame;
- a ruling on the mass of Saturn above, if route A's check against the tables does not
  settle it.

This is preparation. A fit of Adams's data in the project's model is a new analysis and
needs its plan committed and pushed first.
