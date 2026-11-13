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
