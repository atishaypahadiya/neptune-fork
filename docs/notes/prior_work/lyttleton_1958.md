# Lyttleton 1958, a short method for the discovery of Neptune

Reading note, 13 November 2026.

## Citation

R. A. Lyttleton, "A Short Method for the Discovery of Neptune", *Monthly Notices of the Royal
Astronomical Society* 118 (6), 551–559 (1958). DOI 10.1093/mnras/118.6.551. ADS bibcode
1958MNRAS.118..551L. Received 1958 September 25.

**Verified:** the citation as it was given to me (Lyttleton, MNRAS 118, 551, 1958, that title) is
correct. Crossref's record for the DOI
(<https://api.crossref.org/works/10.1093/mnras/118.6.551>) returns author "Lyttleton, R. A.",
title "A Short Method for the Discovery of Neptune", MNRAS volume 118, issue 6, pages 551–559,
issued 1958-12-01, publisher Oxford University Press. The ADS scan's running heads carry the
bibcode, "Vol. 118", "No. 6, 1958" and pages 551–559. This agrees with the entry `lyttleton1958`
in `site/data/m7_references.json`.

## Access

- <https://articles.adsabs.harvard.edu/pdf/1958MNRAS.118..551L>, 2026-11-13, 9 PDF pages,
  660757 bytes, SHA-256 `fddaa13c22e4cf642e4e69b64b5a19f45ce10a0a64c694817bbf8049c508b057`,
  saved as `data/cache/prior_work/ads_1958MNRAS.118..551L.pdf`. Each page carries "© Royal
  Astronomical Society • Provided by the NASA Astrophysics Data System"; cache only.
- Crossref record saved as `data/cache/prior_work/crossref_lyttleton1958.json`, SHA-256
  `717f4a437ac04d7e179122aa530dd19d0560e695466e6d4fad2bc8a77be7bd97`. It has no abstract.
- The ADS abstract page (`https://ui.adsabs.harvard.edu/abs/1958MNRAS.118..551L/abstract`)
  answered 405 to curl; not needed, since the scan has the summary.

## What was read

All nine pages (pp. 551–559) as the scan's own text layer. The layer drops raised decimal points
and degree signs (it prints "3290*4" for 329°·4), so I read p. 557, with Table III and Fig. 2, on
the page image as well; numbers quoted below from other pages are from the text layer and were
checked for sense against the surrounding sentence, not against the image. Figs. 1 and 3 were
not examined beyond their captions.

## Claims

1. The paper asks whether a much less laborious method could have reached the accuracy of Adams
   ("2½°") and Le Verrier ("1°"); the numerical application is "given in full elsewhere", in
   *Vistas in Astronomy* vol. III. [read: p. 551]
2. The work came out of discussions with J. E. Littlewood after the 1946 centenary; the
   conjunction step is "due principally to him". [read: p. 551]
3. Le Verrier's method as Lyttleton describes it: nine apparent unknowns; the mean distance
   assigned by Bode's law; 40 trial values of the epoch longitude at 9° steps, each with a
   least-squares solution of 18 consolidated equations of condition for 1690–1845 in seven
   variables. [read: pp. 551–552]
4. The time of heliocentric conjunction can be found from the longitude discrepancies alone: a
   second difference of the discrepancy about a trial instant, divided by (1 − cos nτ), is
   constant when the instant is conjunction. Applied to Adams's discrepancies for 1690–1840 it
   gives about 1822·3 against the actual 1821·74. [read: pp. 552–553]
5. With conjunction known and a circular orbit at Bode's distance (a′/a = 2), the longitude at
   the discovery date comes out about 13° behind the planet. [read: p. 554]
6. Regrouping the equations of condition as f(t+14) + f(t−14) − f(t), with 14 years being 60°
   of Uranus's motion, removes two unknowns and leaves three (δε, δn, m′) for each assumed
   a/a′; the same regrouping removes the singularity of the coefficient F₂ at the 2:1
   resonance, a/a′ ≈ 0·63. [read: pp. 555–556]
7. The equations were solved for seven values of a/a′ from 0·50 to 0·65, and the mean square
   residual plotted against a′/a (Fig. 2) has "a pronounced minimum at just about a′/a = 1·6",
   a mean distance of 30·71 a.u. against 30·07. [read: pp. 556–557, p. 557 on the image]
8. The longitude predicted for the discovery date is 329°·4 against 328°·4 (equinox 1950), and
   the mass 1/25 000 "in round figures" (Table III prints 25500, with Adams 6670, Le Verrier
   9350, Neptune 19300). [read: p. 557, image]
9. "Some element of good fortune entered into the predictions": Le Verrier's error changes sign
   about 1842, Adams's not until 1856; ±1° holds for about 5 years for Adams's orbit, 12 for
   Le Verrier's, 24 for the circular solution. [read: pp. 557–558]
10. Adams's conversion of corrections of true longitude into corrections of mean longitude by
    the factor r²/ab is called "an invalid step" for perturbed motion, also queried by E. W.
    Brown (M.N. 92, 93, 1931). [read: pp. 558–559]
11. Brown's 1931 method is said to put conjunction near 1840 and the discovery-date longitude
    some 27° ahead. [read: p. 559; Brown's paper itself not read]
12. The companion paper is R. A. Lyttleton, "The rediscovery of Neptune", *Vistas in Astronomy*
    3, 25–46 (1960), DOI 10.1016/0083-6656(60)90004-0. [not read: Crossref metadata only,
    `crossref_lyttleton_vistas_query.json`; Brookes 1970 cites it as "1960, Vistas in Astronomy
    3", read on Brookes p. 80]

## What it means for this repository

Claim 7 is a one-dimensional scan of fit quality over the assumed distance with a circular
orbit, on Adams's discrepancies, and it finds a minimum near 30·7 a.u. That is the nearest thing
I read to the profile over the semi-major axis in `docs/notes/M9.md` and `M10.md`, and it
points the other way from a flat mass–distance ridge: with his reduced equations the distance is
picked out. The two are not the same computation (his is a least-squares scan of three unknowns
on regrouped normal places of Adams; the profile here is over the full model on Le Verrier's
table with an error model), so I record it as a difference to explain, not a contradiction.
Claim 4 is the published form of the conjunction estimate that the site's reference list
credits him with. Claim 10 bears on the Adams reproduction (`P1B_adams_reproduction.md`), which
takes his r²/ab conversion as printed; that note does not discuss Lyttleton's objection.

## Needs a browser / library

- Lyttleton 1960 in *Vistas in Astronomy* 3 (Pergamon, now Elsevier): the numerical application
  (the seven solutions, the weights, the residual tables). Needs a library.
- E. W. Brown, MNRAS 92, 93 (1931), § 14: not sought; it is likely on ADS and would be the next
  thing to fetch.
