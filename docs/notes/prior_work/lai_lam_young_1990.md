# Lai, Lam and Young 1990, perturbation of Uranus by Neptune

Reading note, 13 November 2026. Only the abstract was read.

## Citation

H. M. Lai, C. C. Lam and K. Young, "Perturbation of Uranus by Neptune: A modern perspective",
*American Journal of Physics* 58 (10), 946–953 (October 1990). DOI 10.1119/1.16307.

**Corrected.** The title I was given to check, "Numerical simulation of the discovery of
Neptune", is wrong. Crossref's bibliographic search for the three authors, the journal and the
year (`data/cache/prior_work/crossref_lai1990_query.json`, SHA-256
`5c2de77b07df62700a14fdd845cb4d83b0488940c16e02bd56e05e035b50e24a`) returns as first hit DOI
10.1119/1.16307, authors "Lai, H. M.; Lam, C. C.; Young, K.", title "Perturbation of Uranus by
Neptune: A modern perspective", *American Journal of Physics* volume 58, issue 10, pages
946–953, issued 1990-10-01, publisher American Association of Physics Teachers. Authors,
journal, volume, first page and year were right as given. Two independent citations agree
with Crossref: Kollerstrom 2006, Part 1, a note on p. 24 ("H. M. Lai et al., 'Perturbation of
Uranus by Neptune: A modern perspective', American journal of physics, lviii (1990), 946-53"),
and the caption of his Fig. 6 in Part 2, p. 366 [both read as OCR]. The repository's README
gives authors and year without a title, so nothing there needs changing.

## Access

- Crossref record as above, 2026-11-13; it carries the abstract.
- The publisher's page (`https://doi.org/10.1119/1.16307`, resolving to `pubs.aip.org`)
  answered 403 to curl.
- ADS: `https://articles.adsabs.harvard.edu/pdf/1990AmJPh..58..946L` answered 403 (no scan
  served); the ADS abstract page answered 405.

## What was read

The abstract, as deposited with Crossref by the publisher. Nothing of the article's body,
equations, tables or figures, except that the lower panel of Kollerstrom's Fig. 6 is said in
its caption to be reprinted from p. 950 of this paper; the OCR of that figure is not readable
and I did not study the image.

## Claims

1. "The perturbation of Uranus by Neptune is analyzed in a simple model, in which the
   unperturbed orbits are regarded as circular and coplanar." [read: abstract]
2. "Both the forward problem of calculating the deviation in the position of Uranus and the
   inverse problem of using the observed deviations to infer the elements of Neptune are
   attempted." [read: abstract]
3. "The model accounts for the data quite well, offering a simple understanding of its
   features." [read: abstract]
4. "Because the perturbation is nearly resonant, the effect of Neptune is actually an order of
   magnitude larger than the historically reported remaining deviation." [read: abstract]
5. "The inverse problem also admits a fairly good solution that places Neptune in a
   diametrically opposite position." [read: abstract]
6. The paper is said to show that the perturbation "had reached a minimum over this historical
   period and not a maximum as had been supposed", with a principal term of 52 arcminutes over
   about 4200 years, and to give a reconstruction spanning a millennium (its p. 950). [not
   read: Kollerstrom 2006, Part 2, pp. 366–367, Appendix III and the caption of Fig. 6]
7. Which observed deviations are used (Adams's, Le Verrier's, or a modern set), how the inverse
   problem is solved, and what residual it leaves, I cannot say. [not read]

## What it means for this repository

Claim 5 is the one that matters most: a second, "fairly good" solution on the opposite side of
the Sun, found in a circular coplanar model. The repository's partitioned prior exists because
the posterior has separated modes (`docs/notes/M5b_dynesty_validation.md`, `M6.md`), and I
have not checked whether
any of its cells corresponds to Lai's opposite solution; a statement that the multimodality is
a finding of this work would have to be set against claim 5 first. Claim 4, with Sampson's
p. 147–148 (the long inequality of about 4220 years), says the same thing about what the
historical residuals are: a small remainder after the elements of Uranus absorb most of the
effect. The README's one-line description, "forward and inverse treatment of Uranus's
perturbations by Neptune", matches claims 1 and 2.

## Needs a browser / library

The whole article, pp. 946–953: a subscription to the *American Journal of Physics* or a
library. To be read for: the data used in the inverse problem; the method (least squares on
what unknowns); the residuals; the opposite solution, its longitude and how much worse it
fits; and whether a range of distances is explored.
