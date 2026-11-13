# Phase 2 blinding record

Written 12 November 2026, before any Phase 2 plan exists. It records what has already been
seen, so that a later plan can say honestly what its result is blind to.

## What has been seen

| what | where | when | by what route |
|---|---|---|---|
| Newcomb 1882, the table of observed minus computed contact times, pp. 453–454 (59 rows) | *Astronomical Papers … American Ephemeris* vol. 1, part VI (the source in `data/historical/mercury/MANIFEST.md`) | read 12 Nov 2026 | entered twice during the transcription work, then recognised as residuals and withdrawn |
| The published size of the anomaly: Le Verrier's 38″ a century (1859) and Newcomb's 43″ (1882) | the sources in `prior_art.md` | known before the project began | common knowledge of the subject |
| Newcomb 1895, observed minus theoretical secular variations for four planets (pp. 109–110) | `prior_art.md` section 3, `design_sketch.md` section 5 | read 12 Nov 2026, from OCR | scoping |
| The contact times themselves (58 rows) | `data/historical/mercury/newcomb_1882_contacts.csv` | transcribed 12 Nov 2026 | data, no model applied |

Consequences. The sign and size of Mercury's excess are known and cannot be made unknown.
The pattern of the 1882 residuals across transits has been looked at once, in table form,
without a fit, a plot or a summary statistic. No competing model has been fitted to any
Mercury observation.

## What the engine stays blind to

- Any fit of the real contact times, or of any real Mercury place, to any model (hidden
  mass, modified law, mundane error) until a plan is committed and pushed.
- Any evidence, Bayes factor or posterior from real Mercury data.
- Newcomb's and Le Verrier's residuals as inputs: the engine computes its own from the
  observed quantities. Their residual tables are not read by any code.
- The choice of priors, error model and data cuts: fixed on synthetic data and on the
  sources' stated accuracies, not on how the real residuals look.

## Held-back files

The two entries, the reconciliation and the misprint log of the pp. 453–454 table are kept
outside the repository. They stay out until a Phase 2 result is frozen. They may then be
added as a check on the engine's own residuals, labelled as seen on 12 November 2026.

## Rule for anything further

A source that prints residuals, corrections of elements or a derived excess is listed here,
with the date, before it is opened for transcription. If one is opened by mistake, as the
pp. 453–454 table was, the entry says so.

## 13 November 2026: Le Verrier's Paris places of Mercury, and Newcomb's p. 450

Written before either source is opened for transcription, under the rule above.

### Listed before opening

| what | where | prints residuals? | what is about to be opened, and why |
|---|---|---|---|
| Le Verrier, meridian places of Mercury observed at Paris, 1801–1842 (set C of `data_inventory.md`) | *Annales de l'Observatoire impérial de Paris*, tome V (1859), chapter XV, about pp. 51–61; scan https://articles.adsabs.harvard.edu/pdf/1859AnPar...5....1L | **Yes.** By `data_inventory.md` each row prints, beside the observed apparent longitude and latitude, the error of the tables in longitude and in latitude (observed minus computed) | The scan, to cut page bands for a double-entry transcription of the observed columns only: date and mean time, apparent longitude, apparent latitude, and whatever identifies the observation. It is the cheapest observation-side set among A, C and E of the inventory |
| Newcomb 1882, p. 450, the passage that defines the term in k attached to each contact time | the scan in `data/historical/mercury/MANIFEST.md` | Residual-adjacent: it is the page before the comparison table (pp. 451–452), and pp. 453–454 print the residuals seen on 12 November. Whether p. 450 itself prints a residual-type table is not known before opening | The page, to cut bands of the passage on k for an exact transcription of the text. `transcription.md` says of this page that it "was not read" |

### What is already seen of set C

`data_inventory.md` (12 November) records that pp. 49–52 and 61 of the Annales scan were
read as page images during scoping, and that "the printed table errors scatter by several
arcseconds (pp. 51, 52, 61)". So the size of the scatter of the table errors on three pages
has been seen once, in table form, without a fit, a plot or a summary statistic. That cannot
be undone and is stated here. No value from those columns was written down anywhere in the
repository.

### How the columns are kept out of sight

1. The page layout is found on renders of about 30–40 dpi, where digits cannot be read:
   table frames, column rules and row positions only.
2. The columns are identified from their headings: only the heading strip of a page is
   rendered at full resolution.
3. Pixel boundaries are taken from the rules. The two error columns are painted white or
   cut away, with a margin, before any band of the table body is looked at.
4. Only finished bands are viewed, a few of them, at reduced size first.
5. The error columns are never transcribed, and neither entry is given an image that
   shows them.

If a value of a table error is read by accident, this section will say what, where and how
many. What was in fact seen, by which route and in what form, is appended below when the
bands are made.
