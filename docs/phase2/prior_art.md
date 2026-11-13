# Phase 2 prior art: Mercury's perihelion, hidden mass versus modified law

Scoping note, 2026-11-12. Not a pre-registration. All web access on 2026-11-12.

## How each source was checked

| tag | meaning |
|---|---|
| SCAN | page images or OCR text of the source were read for the statement made |
| ABS | the abstract was read (arXiv API, OpenAlex API or Semantic Scholar API) |
| META | only the bibliographic record was checked (Crossref API or archive.org metadata) |
| UNVERIFIED | stated without a source read for this note; nothing retrieved supports it |

A META tag means the citation exists as written. It says nothing about the content. Where a
content statement sits beside a META tag, the statement is UNVERIFIED unless marked otherwise.
The full source list with identifiers is in section 12.

## 1. Le Verrier 1859

**Source.** "Recherches astronomiques (suite), Chapitre XV, Theorie du mouvement de Mercure",
Annales de l'Observatoire imperial de Paris, tome V, pp. 1-195 (SCAN; the scan has 195 pages and
its last page is printed 195; the volume title page and the year 1859 were not seen, the year is
from the ADS identifier 1859AnPar...5....1L). The short announcement is the letter to Faye,
Comptes rendus 49, pp. 379-383 (SCAN of OCR text; the end page is approximate).

**What it did.** A full analytic theory of Mercury with tables, compared with two kinds of
observation: transits over the Sun and Paris meridian observations.

**Data.**
- Transits: narrative with per-observer contact times, pp. 34-49, from the transit of 1631
  (described, not computed) to that of 9 November 1848 (SCAN of pp. 34, 35, 49).
- The transits are condensed to 21 equations of condition (p. 82, SCAN): November transits of
  1697, 1723, 1736, 1743, 1769, 1782, 1789, 1802, 1848 and May transits of 1753, 1786, 1799, 1832,
  1845. Only internal contacts were used (p. 98, SCAN).
- Meridian: 397 Paris observations, 240 from 8 March 1801 to 22 October 1828 and 157 from
  20 April 1836 to 18 August 1842 (p. 49, SCAN), tabulated as geocentric apparent longitude and
  latitude with the error of the provisional tables, pp. 51-61, in 195 numbered groups.

**Result.** The computed secular motion of the perihelion is 526.7 arcsec per century (p. 99,
SCAN). The transits require an increase: relation (A) on p. 100 reads
288 nu' + 87 nu'' = 38.3 arcsec, with nu' and nu'' the fractional corrections to the masses of
Venus and the Earth (SCAN). Le Verrier states that the need for the increase "results exclusively
from the observations of the transits" (p. 98, SCAN, translated). The letter to Faye gives
"38 secondes" and proposes a planet, or a series of small bodies, between Mercury and the Sun
(SCAN of OCR).

**Settled.** The anomaly is real. Le Verrier's size is low by about 5 arcsec per century against
Newcomb.

**Gap for this project.** None in the physics. The raw contact times on pp. 34-49 are printed
and could be refitted; no machine-readable transcription was found (search limited to the
indexes in section 12).

## 2. Newcomb 1882

**Source.** "Discussion and results of observations on transits of Mercury, from 1677 to 1881",
Astronomical Papers prepared for the use of the American Ephemeris and Nautical Almanac, vol. 1,
part VI. The volume contents give the start as p. 363 (OCR reads "36-", last digit damaged); the
memoir's own contents run from p. 367 to the corrigenda on p. 484 (SCAN of OCR). The title above
is the usual citation form; the volume contents print "Discussion of Observed Transits of
Mercury, 1677-1881" (SCAN of OCR).

**What it did.** Collected every usable contact observation of the transits from 1677 to 1881,
reduced each to a Greenwich mean time of geocentric contact, and solved for corrections to
Le Verrier's elements. It also tested whether the Earth's rotation is uniform (a factor k on the
lunar-derived clock correction).

**Data.** Tabular summary of observations by transit, pp. 384-421. Comparison of observed and
tabular geocentric contact times, pp. 451-452. Equations of condition from p. 453 (SCAN of OCR).

**Result.** The transits fix two linear functions of the element corrections, V (November) and
W (May), and a node function N. With the secular variation of the eccentricity taken from
theory, "the observed centennial motion of the perihelion of Mercury is greater by 43 arcsec than
the theoretical motion" (p. 473, SCAN of OCR; the solved value is p = +42.95).

Section 5 of Part III (pp. 474-477) already lists the candidate causes and the tests between
them (SCAN of OCR):
- a planet or group of planets inside Mercury's orbit, rejected because it would be seen at
  eclipses and in transit, and because a group in the ecliptic would move the node of Mercury as
  much as the perihelion, which is not observed;
- matter like that of the zodiacal light, with the same node difficulty;
- an ellipticity of the Sun, rejected from measures of the solar disk;
- a ring, left open as "a field for investigation";
- an inverse third or fourth power term, rejected by its size at laboratory distances;
- Weber's velocity-dependent law, after Tisserand, which gives too little.

**Settled.** 43 arcsec per century from transits. The list of discriminators (node of Mercury,
node of Venus, visibility) is Newcomb's, not new.

## 3. Newcomb 1895

**Source.** "The Elements of the Four Inner Planets and the Fundamental Constants of Astronomy",
Supplement to the American Ephemeris and Nautical Almanac for 1897, Washington, 1895
(archive.org metadata and OCR text; SCAN of OCR for the statements below).

**What it did.** A joint solution for the elements and secular variations of Mercury, Venus, the
Earth and Mars from 62,030 meridian observations (Mercury 5,421; p. 3) and the transit results of
1882, then a comparison of observed and theoretical secular variations (table, pp. 109-110) and
an examination of hypotheses (chapter VI, sections 55-63).

**Result (pp. 109-110, arcsec per century, with Newcomb's mean errors, read from OCR).**

| quantity | observation minus theory |
|---|---|
| Mercury, e d(varpi)/dt | +8.48 +/- 0.43 |
| Mercury, de/dt | -0.88 +/- 0.50 |
| Mercury, di/dt | +0.38 +/- 0.80 |
| Mercury, sin i d(Omega)/dt | +0.61 +/- 0.52 |
| Venus, e d(varpi)/dt | -0.05 +/- 0.25 |
| Venus, di/dt | +0.38 +/- 0.33 |
| Venus, sin i d(Omega)/dt | +0.60 +/- 0.17 |
| Earth, e d(varpi)/dt | +0.10 +/- 0.13 |
| Mars, e d(varpi)/dt | +0.75 +/- 0.35 |

The digits were read from OCR text and need a check against the page images before any use.

Newcomb names four deviations that call for attention: the perihelion of Mercury, the node of
Venus, the perihelion of Mars and the eccentricity of Mercury (p. 110). He examines a
non-spherical Sun (needs an ellipticity of "rather more than half a second of arc", excluded by
heliometer measures), an intra-Mercurial ring, zodiacal-light matter, and a ring between Mercury
and Venus, which he calls the hypothesis that best fits and then sets aside as improbable. He
adopts Hall's hypothesis provisionally: exponent 2.000 000 1574, which gives perihelion motions
of 42.34 (Mercury), 16.58 (Venus), 10.20 (Earth) and 5.42 (Mars) arcsec per century (pp. 118-119).

**Settled.** The size of the Mercury anomaly. The Venus node excess did not survive; where and
when it was resolved was not verified here (UNVERIFIED).

**Gap.** Chapter VI is a hypothesis-by-hypothesis elimination by hand. Nobody is known to have
redone it as a formal model comparison with stated priors and error models (not found; see
section 11 for the limits of that search).

## 4. Hall 1894 and Seeliger

- Hall, "A suggestion in the theory of Mercury", Astronomical Journal 14, 49 (META). Newcomb's
  footnote in the 1895 OCR text reads "Vol. XIV, p. 7"; the OCR may have dropped a digit. The
  hypothesis: gravity varies as r to the power -(2 + delta). Every perihelion advances at
  (delta / 2) times the mean motion. What ended it is usually given as the Moon's perigee
  (UNVERIFIED); Newcomb states that the lunar excess would be about 140 arcsec per century
  (p. 119, SCAN of OCR).
- Seeliger, "Das Zodiakallicht und die empirischen Glieder in der Bewegung der innern Planeten",
  Sitzungsberichte der Bayerischen Akademie, 1906 (UNVERIFIED: volume, pages and exact title not
  confirmed by any index queried). His later paper "Ueber die Anomalien in der Bewegung der innern
  Planeten", Astronomische Nachrichten 201, 273-280 (1915) is META. The reply by Freundlich,
  Astronomische Nachrichten 201, 49-56 (1915) is META (Crossref spells the author "Freudlich").
  Content, UNVERIFIED: Seeliger fitted two ellipsoidal dust distributions that reproduced
  Newcomb's discrepancies; Freundlich argued the required mass exceeds what the zodiacal light
  allows.

**Gap.** Seeliger's model is the historical instance of a flexible hidden-mass model that fits
everything. It is the natural test of an Occam penalty. That use is not new in spirit (section 9).

## 5. Vulcan claims

- Lescarbault, "Passage d'une planete sur le disque du Soleil, observe a Orgeres", letter to
  Le Verrier, Comptes rendus 50, p. 40 onward (SCAN of OCR; start page from the OCR page marks,
  end page not checked), followed by Le Verrier's remarks.
- Newcomb, "On the supposed intra-Mercurial planets", Astronomical Journal 6, 162 (1860) (META).
  Newcomb 1882 cites it for the argument that a group massive enough would transit too often to
  escape notice (SCAN of OCR of the 1882 memoir).
- Eclipse of 29 July 1878: Watson, "Discovery of an intra-Mercurial planet", American Journal of
  Science (3) 16, 230-233, and "On the intra-Mercurial planets", 310-313 (META). Swift, "Discovery
  of Vulcan", Nature 18, 539 (META), and letters in Astronomische Nachrichten 93, 161-166 (META).
- Histories: Fontenrose, "In Search of Vulcan", Journal for the History of Astronomy 4, 145-158
  (1973) (META). Baum and Sheehan, "In Search of Planet Vulcan", 1997 (META via chapter records).

**Settled.** No intra-Mercurial planet exists at the mass required (UNVERIFIED as a sourced
statement; it follows from sections 6 and 7).

**Gap.** The non-detections are not orbital data. A dynamics-only engine cannot use them unless
they enter as a prior on a discrete body. That is a design choice, not a research gap.

## 6. Clemence 1947 and the mid-century position

- Clemence, "The Relativity Effect in Planetary Motions", Reviews of Modern Physics 19, 361-364
  (META). Content UNVERIFIED (no abstract available): it is the standard reference for the
  comparison of observed and relativistic perihelion advances of Mercury and the Earth.
- Clemence, "The motion of Mercury 1765-1937", Astronomical Journal 50, 126 (1943) (META); the
  full memoir in the Astronomical Papers series was not located (UNVERIFIED).
- Morrison and Ward, "An Analysis of the Transits of Mercury: 1677-1973", MNRAS 173, 183-206
  (1975) (ABS and SCAN of pp. 184-185). About 2400 timings of internal contact, transcribed from
  the original publications except 50 taken from Newcomb 1882. Clock corrections from lunar
  observations. Result: excess perihelion motion 41.9 +/- 0.5 arcsec per century, and a lunar
  tidal acceleration. The smallest standard deviation of the timings of any one transit is 11 s.
- Narlikar and Rana, MNRAS 213, 657-663 (1985) (ABS): a Newtonian N-body integration gives
  about 528.95 arcsec per century.
- Nobili and Will, "The real value of Mercury's perihelion advance", Nature 320, 39-41 (1986)
  (META).

**What would not be novel.** Refitting the transit timings with a modern ephemeris and getting
about 42 to 43 arcsec per century. Morrison and Ward did it with the full set of timings.

## 7. Solar oblateness

- Dicke and Goldenberg, Physical Review Letters 18, 313-316 (1967) (ABS): fractional difference of
  equatorial and polar radii (5.0 +/- 0.7) x 10^-5, implying a discrepancy of 8 per cent of the
  Einstein value. Full paper: Astrophysical Journal Supplement 27, 131 (1974) (META).
- Hill et al., Physical Review Letters 33, 1497-1500 (1974) (META). Hill and Stebbins,
  Astrophysical Journal 200, 471 (1975) (ABS): equatorial minus polar diameter
  18.4 +/- 12.5 milliarcsec; the abstract says this removes the consequence of the 1967 result for
  general relativity.
- Later shape and gravity measurements: Pijpers, MNRAS 297, L76-L80 (1998) (META); Fivian et al.,
  Science 322, 560-562 (2008) (META); Kuhn et al., Science 337, 1638-1640 (2012) (META); review by
  Rozelot and Damiani, European Physical Journal H 36, 407-436 (2011) (META).
- Dynamical values: J2 = (2.25 +/- 0.09) x 10^-7 (Park et al. 2017, ABS);
  J2 = (2.246 +/- 0.022) x 10^-7 (Genova et al. 2018, ABS).

**Settled.** The solar quadrupole contributes about 0.03 arcsec per century to Mercury's
perihelion (computed in `design_sketch.md` from J2 = 2.25 x 10^-7).

**Gap.** None. J2 is the modern form of the "hidden mass" branch and it is measured.

## 8. Modern ephemeris tests

- Shapiro et al., Physical Review Letters 28, 1594-1597 (1972) (ABS): radar, (2 + 2 gamma - beta)/3
  = 1.005 +/- 0.007.
- Pitjeva and Pitjev, MNRAS 432, 3431-3437 (2013) (ABS): EPM2011, about 678 thousand observations
  from 1913 to 2011; PPN parameters, solar oblateness, limits on dark matter inside Saturn's
  orbit. Pitjev and Pitjeva, Astronomy Letters 39, 141-149 (2013) (META).
- INPOP: Fienga et al., Celestial Mechanics and Dynamical Astronomy 111, 363-385 (2011) (META);
  Verma et al., Astronomy and Astrophysics 561, A115 (2014) (ABS); Fienga et al., Celestial
  Mechanics and Dynamical Astronomy 123, 325-349 (2015) (META); review by Fienga and Minazzoli,
  Living Reviews in Relativity 27, 1 (2024) (META).
- Park et al., Astronomical Journal 153, 121 (2017) (ABS): MESSENGER ranging; total perihelion
  precession 575.3100 +/- 0.0015 arcsec per century; the periodic 1PN signature in the argument
  of perihelion is measured, which separates beta from J2 when gamma is taken from Cassini.
- Genova et al., Nature Communications 9, 289 (2018) (ABS).
- Will, Physical Review Letters 120, 191101 (2018) (ABS): a further contribution from cross terms
  with the other planets at a few parts in 10^6 of 42.98 arcsec per century.
- Kolesnik and Masreliez, Astronomical Journal 128, 878-888 (2004) (ABS): about 240,000 optical
  observations of the Sun, Mercury and Venus, 1750-2000, against DE405.

**Settled.** General relativity holds at the 10^-5 level in the inner Solar System. Modern data
beat the transit data by about three orders of magnitude (0.0015 against roughly 1.3 arcsec per
century; the second figure is derived in `design_sketch.md`).

**What would not be novel.** Any constraint on PPN parameters, J2, dark matter or a fifth force.
A project using 19th-century timings cannot add to these.

## 9. Bayesian model selection: modified gravity against unseen mass

What the named papers did (all ABS):

| paper | method | evidence ratio, mass against law? |
|---|---|---|
| Brown and Mathur, AJ 166, 168 (2023) | secular theory with the MOND external field | no |
| Migaszewski, MNRAS 525, 805-826 (2023) | N-body integrations with the external field | no |
| Vokrouhlicky et al., ApJ 968, 47 (2024) | population simulations, Newton and AQUAL | no |

- Brown and Mathur show analytically that the external field can align the orbits of distant
  Kuiper belt objects, as an alternative to Planet Nine.
- Migaszewski integrates the known distant objects under the external field effect.
- Vokrouhlicky, Nesvorny and Tremaine simulate comet and detached-disk populations; the MOND runs
  fail to reproduce the observed binding energies. The comparison is between simulated and
  observed distributions.

None of the three computes Bayesian evidence between a hidden mass and a modified law.

Closest work found:

- Mariani et al., Physical Review D 108, 024047 (2023) (ABS) and Astronomy and Astrophysics 682,
  A175 (2024) (ABS): MCMC over a modified-gravity parameter with a full ephemeris rebuild (INPOP).
  Posterior for one family; no hidden-mass alternative; no evidence ratio.
- Hees et al., Physical Review D 89, 102002 (2014) (META); Park et al., Physical Review D 114,
  024066 (2026) (ABS); Desmond, Hees and Famaey, MNRAS 530, 1781-1795 (2024) (ABS): the MOND
  quadrupole from Cassini data and its tension with galaxy data. Bayesian inference in the 2024
  paper; no hidden-mass branch.
- Planet Nine from dynamics: Fienga et al., Astronomy and Astrophysics 587, L8 (2016) and 640, A6
  (2020); Holman and Payne, AJ 152, 80 and 152, 94 (2016); Brown and Batygin, AJ 162, 219 (2021)
  (all META). Hidden-mass branch only.
- Banik et al., MNRAS 527, 4573-4615 (2024) (ABS): wide binaries, Newton against MOND through an
  interpolating parameter, with undetected close companions as a nuisance population. This is the
  nearest structural match: a law parameter and a hidden-mass confound in one statistical model.
- Yee et al., Astrophysical Journal Letters 888, L5 (2020) (ABS): WASP-12b transit times; orbital
  decay against apsidal precession, with a Bayes factor of 70,000, and radial velocities used to
  rule out acceleration by a companion. An exoplanet instance of the three-way choice.
- Milgrom, MNRAS 399, 474-486 (2009) and Blanchet and Novak, MNRAS 412, 2530-2542 (2011) (META):
  MOND effects in the inner Solar System, including an anomalous perihelion precession (content
  UNVERIFIED). A modified law that predicts a perihelion term; no model comparison.
- Standish, AJ 105, 2000 (1993) (META): "Planet X: no dynamical evidence in the optical
  observations". The mundane-error branch winning (content UNVERIFIED).
- Berger and Jefferys, Journal of the Italian Statistical Society 1, 17-32 (1992) (META), and
  Jefferys and Berger, "Ockham's razor and Bayesian analysis", American Scientist, 1992
  (UNVERIFIED citation). Content UNVERIFIED but important: these are remembered as computing a
  Bayes factor between general relativity and a Newtonian law with an adjustable parameter, using
  one number, the observed perihelion advance of Mercury. If so, the summary-level Bayesian
  comparison for Mercury exists as a textbook illustration. It must be read before any novelty
  claim.

**What would not be novel.** A Bayes factor between Einstein's value and a one-parameter fudge,
computed from the single summary number 43 +/- something.

## 10. Roseveare 1982

Roseveare, "Mercury's Perihelion from Le Verrier to Einstein", Oxford, Clarendon Press, 1982,
ISBN 0198581742 (META from archive.org; the scan is lending-only and was not read). Reviews:
Isis 76, 246-247; American Journal of Physics 51, 478-479; Science 219, 1058 (META). Content
UNVERIFIED: the standard history of every hypothesis in sections 1 to 6, with the numbers.

**What would not be novel.** The narrative, and any qualitative statement of which hypothesis
failed on which observation.

## 11. Is there a gap at all?

Frankly: not in physics, and only a narrow one in method.

Not novel, and not to be claimed:
1. The size of the anomaly from transit data (Le Verrier, Newcomb, Morrison and Ward).
2. General relativity against alternatives (ephemeris groups, at precision out of reach here).
3. The list of discriminators: nodes of Mercury and Venus, other perihelia, visibility of the
   matter (Newcomb 1882 and 1895).
4. Bayesian Occam arguments on the single number (Berger and Jefferys, if confirmed).
5. The history (Roseveare; Baum and Sheehan).

What was not found, within the limits below:
1. A blind, pre-registered decision procedure that is given positions or timings only, and
   returns calibrated evidence for hidden mass, modified law or data error, with its error rates
   measured on synthetic anomalies first.
2. The same procedure applied to two historical anomalies with opposite answers: Uranus, where
   the mass was real, and Mercury, where the law was wrong. Phase 1 supplies the first.
3. A quantified statement of what the 19th-century data could decide. The arithmetic in
   `design_sketch.md` suggests the answer is "not much from Mercury alone", and that Newcomb's own
   numbers lean toward matter because of the Venus node. A formal version of that statement, with
   priors exposed, was not found.
4. An open, checked transcription of Newcomb's 1882 contact table with provenance (existence of
   one was not established either way; Morrison and Ward's collection may or may not be public).

Limits of the search. Bibliographic checks used Crossref, arXiv, OpenAlex, Semantic Scholar and
archive.org queries by title and keyword. ADS abstract pages, Gallica and HathiTrust refused
plain requests. No citation-graph search was possible. "Not found" means
not found by these means. Before any public novelty claim, a person should search ADS for papers
citing Morrison and Ward 1975 and Berger and Jefferys 1992.

So the gap is a methods and history-of-science contribution: a calibrated, pre-registered fork
tested on two cases. The likely honest result for Mercury is "the data of 1882 detect the anomaly
at high significance and cannot say what causes it; the verdict is set by the priors". That is a
publishable negative only if it is stated as such from the start.

## 12. Sources

| ref | identifier or location | tag |
|---|---|---|
| Le Verrier 1859, Annales tome V, pp. 1-195 | https://articles.adsabs.harvard.edu/pdf/1859AnPar...5....1L | SCAN |
| Le Verrier, letter to Faye, CR 49, 379-383 | https://archive.org/details/comptesrendusheb49acad | SCAN (OCR) |
| Lescarbault, CR 50, 40 onward | https://archive.org/details/comptesrendusheb50acad | SCAN (OCR) |
| Newcomb 1882, APAE 1, part VI | https://archive.org/details/astronomicalpap00newcgoog | SCAN (OCR) |
| Newcomb 1895 | https://archive.org/details/elementsoffourin00newcuoft | SCAN (OCR) |
| Newcomb 1860 | doi:10.1086/100803 | META |
| Hall 1894 | doi:10.1086/102055 | META |
| Seeliger 1915 | doi:10.1002/asna.19152011502 | META |
| Freundlich 1915 | doi:10.1002/asna.19152010302 | META |
| Seeliger 1906 | none found | UNVERIFIED |
| Watson 1878 | doi:10.2475/ajs.s3-16.93.230, doi:10.2475/ajs.s3-16.94.310 | META |
| Swift 1878 | doi:10.1038/018539a0, doi:10.1002/asna.18780931102 | META |
| Fontenrose 1973 | doi:10.1177/002182867300400301 | META |
| Baum and Sheehan 1997 | doi:10.1007/978-1-4899-6100-6_14 (one chapter) | META |
| Roseveare 1982 | https://archive.org/details/mercurysperiheli0000rose | META |
| Einstein 1915 | reprint doi:10.1002/3527608958.ch4; original pages not checked | META (reprint) |
| Clemence 1943 | doi:10.1086/105746 | META |
| Clemence 1947 | doi:10.1103/RevModPhys.19.361 | META |
| Morrison and Ward 1975 | doi:10.1093/mnras/173.1.183; https://articles.adsabs.harvard.edu/pdf/1975MNRAS.173..183M | ABS, SCAN |
| Narlikar and Rana 1985 | doi:10.1093/mnras/213.3.657 | ABS |
| Nobili and Will 1986 | doi:10.1038/320039a0 | META |
| Dicke and Goldenberg 1967 | doi:10.1103/PhysRevLett.18.313 | ABS |
| Dicke and Goldenberg 1974 | doi:10.1086/190292 | META |
| Hill et al. 1974 | doi:10.1103/PhysRevLett.33.1497 | META |
| Hill and Stebbins 1975 | doi:10.1086/153813 | ABS |
| Pijpers 1998 | doi:10.1046/j.1365-8711.1998.01801.x | META |
| Fivian et al. 2008 | doi:10.1126/science.1160863 | META |
| Kuhn et al. 2012 | doi:10.1126/science.1223231 | META |
| Rozelot and Damiani 2011 | doi:10.1140/epjh/e2011-20017-4 | META |
| Shapiro et al. 1972 | doi:10.1103/PhysRevLett.28.1594 | ABS |
| Pitjeva and Pitjev 2013 | doi:10.1093/mnras/stt695 | ABS |
| Pitjev and Pitjeva 2013 | doi:10.1134/S1063773713020060 | META |
| Fienga et al. 2011 | doi:10.1007/s10569-011-9377-8 | META |
| Verma et al. 2014 | doi:10.1051/0004-6361/201322124 | ABS |
| Fienga et al. 2015 | doi:10.1007/s10569-015-9639-y | META |
| Fienga and Minazzoli 2024 | doi:10.1007/s41114-023-00047-0 | META |
| Park et al. 2017 | doi:10.3847/1538-3881/aa5be2 | ABS |
| Genova et al. 2018 | doi:10.1038/s41467-017-02558-1 | ABS |
| Will 2018 | doi:10.1103/PhysRevLett.120.191101 | ABS |
| Kolesnik and Masreliez 2004 | doi:10.1086/422344 | ABS |
| Brown and Mathur 2023 | doi:10.3847/1538-3881/acef1e; arXiv:2304.00576 | ABS |
| Migaszewski 2023 | doi:10.1093/mnras/stad2250; arXiv:2303.13339 | ABS |
| Vokrouhlicky et al. 2024 | doi:10.3847/1538-4357/ad40a3; arXiv:2403.09555 | ABS |
| Mariani et al. 2023 | doi:10.1103/PhysRevD.108.024047 | ABS |
| Mariani et al. 2024 | doi:10.1051/0004-6361/202348082 | ABS |
| Hees et al. 2014 | doi:10.1103/PhysRevD.89.102002 | META |
| Park et al. 2026 | doi:10.1103/r7n8-kw38 | ABS |
| Desmond et al. 2024 | doi:10.1093/mnras/stae955 | ABS |
| Fienga et al. 2016, 2020 | doi:10.1051/0004-6361/201628227, doi:10.1051/0004-6361/202037919 | META |
| Holman and Payne 2016 | doi:10.3847/0004-6256/152/4/80, doi:10.3847/0004-6256/152/4/94 | META |
| Brown and Batygin 2021 | doi:10.3847/1538-3881/ac2056 | META |
| Banik et al. 2024 | doi:10.1093/mnras/stad3393 | ABS |
| Yee et al. 2020 | doi:10.3847/2041-8213/ab5c16 | ABS |
| Standish 1993 | doi:10.1086/116575 | META |
| Berger and Jefferys 1992 | doi:10.1007/BF02589047 | META |
| Jefferys and Berger 1992, American Scientist | none found | UNVERIFIED |
| Milgrom 2009; Blanchet and Novak 2011 | doi:10.1111/j.1365-2966.2009.15302.x; doi:10.1111/j.1365-2966.2010.18076.x | META |

Notes on names. The arXiv record of Brown and Mathur lists the first author as Katherine
Jones-Smith; the journal record lists Katherine Brown. Diacritics are dropped in this file.
