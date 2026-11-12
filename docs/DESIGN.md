# Site design (M7): a paper illustrated with plates

Written 2026-11-09. This replaces the earlier spec (42df69f), which described dark three.js
scenes. That direction is retired. The visual reference is `design-ref/neptune-fork site
directions.html`, which is local only, excluded from git, never committed, and only its
"Cyanotype" page is used. The reference sets the look. Its numbers and words are placeholders,
and none of them is used. Only its copy structure carries over: a plate label, a section mark,
a title that says what was done, prose that refers to lettered points, and a standalone caption.

## Concept

The site reads like a scientific paper illustrated with 19th-century astronomical plates. It
has two layers, and they never mix:

- **The archive layer (figures):** cyanotype ground, grain, dust and glow. Every figure is a
  plate: white line on blue.
- **The journal layer (text):** clean, sharp, flat paper, no texture. That covers the title
  block, prose, captions, Table 1, references and controls.

Texture is applied only to a figure's ground group. Labels, letters and numbers inside a figure
sit above the texture and stay sharp.

## The one rule

**No number or mark on the page that doesn't come from `site/data/` or a cited source.**
- Every number is read from a JSON file at load time and formatted in code. Nothing numeric is
  typed into HTML or CSS.
- A cited constant, such as Le Verrier's published elements, comes from a JSON field that
  carries its citation. The citation is shown in the caption or in the references.
- Every mark (an orbit, a point, a stipple dot, a tick) is drawn from data. Decorative marks are
  limited to the ground textures, which encode nothing and are named as texture in the code.
- A pending value says it is pending, taken from the JSON status (as now: `m8b.json`,
  `m8_arbiter.json`).
- **Check:** the screenshot script collects every number in the rendered text and fails on any
  that is not in an allowlist built from the loaded JSON and the citation list.

## Palette (cyanotype)

| token | hex | used for |
|---|---|---|
| `--deep` | #0f2a6b | darkest ground in figures; body text on paper |
| `--prussian` | #1c3f94 | figure ground base; links, Table 1 values |
| `--mid` | #2f62b8 | the light cloud layer in the ground only |
| `--wash` | #9db8dc | vignette at plate edges; the light wash clouds |
| `--paper` | #eef1f3 | journal ground |
| `--line` | #f4f7fb | every line, point, stipple and label in figures |

- **Contrast:** text on paper uses `--deep` (about 12:1) and `--prussian` for captions and
  values (about 8:1). Text on a plate is `--line` and is never placed over a `--wash` area.
- **Opacity in figures:** three steps only. 1 for data lines and points, 0.6 for secondary
  marks (the 95 % band edge, dotted orbits), 0.3 for context (the frame, gridlines in Fig. 3).
- **No other colours.** There is no status colour: pending and provisional are said in words
  and shown as dashed lines.

## Type

- **GFS Didot:** titles (the h1 on Plate I and the section h2s) and large numerals (the year
  on Plate III). Nothing else.
- **Source Serif 4** (optical sizes, roman and italic): everything else. That includes prose,
  captions, Table 1, figure labels, the labels along orbits and the controls.
- **Small caps:** only for "Plate", "Fig.", "Table" and "§". The reference also sets Table 1's
  source names in small caps; that is not carried over.
- **Italic:** the point letters (in figures and prose), mathematical variables, and titles of
  works in the references. Not captions, not emphasis.
- **Sizes, from the reference at 1280 px:**
  - title 112 px Didot;
  - h2 40/42 px Didot;
  - body 18/28 px;
  - captions 14/21 px;
  - Table 1 cells 20 px;
  - labels along orbits 13 px with 1.2 px letter spacing;
  - point letters 16–20 px italic;
  - plate label 18 px small caps with 2 px letter spacing.

  Phones scale these down with `clamp()`.
- **Handwriting:** only as faint texture, and only as an image of a real manuscript, with its
  source text cited beside it (archive, shelfmark, date, transcription). The reference's script
  is a modern typeface (Monsieur La Doulaise), not 19th-century handwriting, so it is not used.
  Until a scan with a citation and a usable licence is confirmed, there is no handwriting on
  the page.

## Figures

- **Medium:** SVG, drawn with D3 (pinned, from jsDelivr). three.js is removed.
- **Ground:** a stack of SVG filters with fixed seeds, so every load looks the same:
  - a dark low-frequency `feTurbulence` mottle;
  - a light cloud layer (`--mid`/`--wash`);
  - fine grain;
  - sparse dust specks;
  - a radial `--wash` vignette toward the plate edges;
  - a hairline plate mark inset from the edge (Plate I only).
- **Glow:** a small Gaussian blur under the bodies on Plate III and the highlighted point, and
  nowhere else.
- **Line:** white line `--line`.
  - Orbits are dotted (round caps, short dash).
  - Sightlines and the pull line are solid hairlines.
  - Model curves are 2 px solid.
  - Bodies are stippled: a cluster of dots, never a filled disc. The Sun's stipple is texture
    and encodes nothing. The perturber's stipple encodes the posterior (below).
- **Points lettered in italic:** *S* Sun, *U* Uranus, *N* Neptune (JPL), *C* the 1821
  conjunction, *L* Le Verrier's planet. Any other letter is defined in its caption.
  - The prose refers to points by letter. Each letter in the prose is a focusable element
    linked to its point; hover or focus highlights the point (glow, label at full weight).
  - The link works by keyboard and is announced (`aria-describedby` points to the caption).
- **Labels along orbits** (SVG `textPath`) carry only data text: observer names and years from
  the CSV, and fitted values from the JSON.
- **Captions:** every figure has a numbered caption ("Plate II.", "Fig. 3.") that stands alone.
  It states:
  - what is drawn;
  - which run (for example, the M6c replacement posterior);
  - the settings that matter (cells, nlive, frames);
  - the scale ("to scale", or "not to scale" with the factor);
  - any exaggeration factor, read from the JSON.
- **Posterior stipple (Plate II):** one dot per posterior sample at *C*. `m6c_conjunction.json`
  holds quantiles and histograms, not samples, so the restyle adds the samples' heliocentric
  longitudes at *C*. A script in `docs/notes/scripts/` reads the published M6c replacement
  posterior; `src/` is not changed. A seeded subsample of 600 is enough to read the density.

## Voice

A scientific paper.
- Titles say what was done. Every number carries its interval, or says why it has none.
- The verbs are "we find", "consistent with" and "we did not test". No promotional language,
  no rhetorical questions, no exclamation marks.

**Structure.** Sections are numbered, because the text cross-references them ("§ 4.2"):
1. Abstract (on Plate I)
2. Data: Table 1, observations by source; Fig. 1, the observations over time; the double
   entry and the spot check.
3. Method:
   - Fig. 2, the pre-registration timeline;
   - Plate III with Fig. 3, the perturbation;
   - the observation model, likelihood and priors in prose.
4. Results:
   - 4.1, the posterior on Galle's night (Fig. 4);
   - 4.2, the conjunction of 1821 (Plate II);
   - 4.3, the window in which Le Verrier's planet stays within the 68 % band (Fig. 5).
5. Robustness: Fig. 6, the repeatability correction and S1–S3, with the S1 verdict pending
   the arbiter.
6. Limitations: A6a, A6b, the provisional 95 % interval, S1 pending, the unchecked rows.
7. Data availability: the DOIs, CC0, the citations and the code licence.
8. References, with DOIs where they exist.

**Copy to rewrite.** The site's current copy and README.md are rewritten in this voice. The
numbers stay JSON-driven, and the claims stay the same as in A6.md, M6.md and M8.md.

## Plate I: title

- **Title block** (journal layer, over the plate):
  - a title that says what was done, drafted at the restyle from the project's own description;
  - the byline: author, year, the repository and the version DOI;
  - the abstract follows in § 1.
- **Figure:** *S*, Uranus's orbit and the outer orbit band:
  - the Uranus orbit (dotted) carries, along a `textPath`, the observers and years from the
    CSV;
  - the outer band shows the 68 % range of the perturber's longitude on Galle's night
    (`m6c_replacement.json`);
  - *N* is JPL's Neptune; *L* is Le Verrier's planet from his published elements (the cited
    sources in `m6c_conjunction.json`).
  - The caption states the frame, the epoch and that radii are not to scale, with the factor.

## Plate II: the conjunction of 1821 (§ 4.2)

From `m6c_conjunction.json`, heliocentric, ecliptic of date, at its stated epoch:
- *S*, *U* on its orbit, and the line *SU* produced to *C*;
- the posterior stipple at *C* (68 % and 95 % read from the JSON);
- *L* (Le Verrier's planet at the same epoch) and *N* (JPL's Neptune);
- the angle *CSL*, labelled with its value from the JSON;
- the dotted orbits.

The paragraph's numbers and the caption ("not to scale; radii ×…") come from the same file. The
file says it is an after-the-fact diagnostic, not pre-registered, and the caption says so too.

## Table 1: observations by source (§ 2)

- One row per source: Flamsteed, Bradley, Mayer, Le Monnier, Greenwich meridian
  observations, Paris 1801–1827 and Paris 1835–1845. The rows follow the CSV's `observer` and
  `observatory`, and the eras follow the error groups.
- **Columns:** source; observatory; years; observations (rows) and measurements (RA, Dec);
  median |O − JPL| in arcseconds, with N.
- **Built by** `docs/notes/scripts/m7_table1.py` → `site/data/m7_table1.json`, from
  `data/historical/leverrier_1846.csv`. The per-source |O − JPL| values come from
  `fork.data.calibration`, which the script imports; `src/` is not changed.
- **Footer row:** totals; "28 rows checked by hand against the scan"; the licence (CC0) and the version
  DOI from `m7_dataset.json`.
- **Caption:** "Table 1." in small caps, then what the residual is measured against.
- The values are typeset in Source Serif 4 tabular figures, with no rules except one hairline
  above the totals.

## Plate III: the pull on Uranus (§ 3), with Fig. 3

**Integration** (`docs/notes/scripts/m7_motion.py` → `site/data/m7_motion.json`):
- **Code:** REBOUND, IAS15. REBOUND stays a dev dependency, and the lint rule bans it only from
  `src/`, which this script is not in.
- **Bodies:** the Sun and all eight planets (Mercury, Venus, the Earth–Moon barycentre, Mars,
  Jupiter, Saturn, Uranus, Neptune).
- **Initial states:** JPL Horizons via `fork.data.horizons` (the cache), at 1781-03-13.0 TDB
  (Herschel's discovery; `A1_EPOCH_JD_TDB`). Masses come from `fork.model.bodies`, with its
  citations.
- **Two runs:** with Neptune, and with Neptune massless (it still moves, as a test particle).
- **Output:** every 0.1 yr, 1781–1850:
  - the heliocentric ecliptic x, y of *U*, *N* and the massless-run Uranus;
  - Neptune's heliocentric acceleration on Uranus, with its along-track sign;
  - the heliocentric longitude difference of Uranus (with minus without), in arcseconds;
  - the conjunction date, found from the run and checked against `m6c_conjunction.json`'s JPL
    date.
- **Validation:** the run with Neptune is compared with Horizons' Uranus in 1850, and the
  difference is recorded in the JSON and in M7.md.

**Drawing:**
- *S*, the dotted orbits of *U* and *N*, and the trail of *U*.
- The pull line *UN*, labelled with its along-track sign: "drawn ahead" or "held back".
- The massless-run Uranus as an open circle, its offset multiplied by a factor stored in the
  JSON and stated in the caption.
- *C* marked on the line from *S*.
- **Fig. 3:** the longitude difference in arcseconds against year, with a cursor at the current
  year and *C* marked.
- **Controls:** play/pause (a 44 px button) and a year scrubber (a range input, with arrow keys).
  The year is shown in Didot numerals.
- **Reduced motion:** starts paused at the conjunction year.

## Interaction

- **Letter highlight:** text ↔ figure, as above.
- **Plate III:** the controls above.
- **Hover:** a quiet probe on Figs. 1, 4, 5 and 6, showing values from the JSON. Nothing appears
  only on hover: every number a reader needs is in the text, the caption or the table.

## Layout

- **Desktop:** a plate is a full-bleed or half-bleed figure, with its text column on paper
  beside it (as in the reference's Plate II and Plate III). Sections scroll; there is no scroll
  snapping.
- **Phone (≤ 720 px):** the figure comes first, then the text. Labels along orbits are reduced
  to the outer one. Plate III's Fig. 3 goes below the plate, with the controls fixed under it
  while the plate is in view.

## Accessibility

- **Contrast:** as above.
- **Each figure:** an `aria-label` and a hidden long description built from the same JSON (as
  now).
- **Interaction:** focus rings on the letters and controls, and reduced motion respected.
- **Table 1:** a real `<table>` with a `<caption>` and header cells.

## Archival quotes

None is on the page until it is confirmed. The candidates, each with its exact source, are
listed in the restyle review. The reference's quotes are placeholders and are not used as
sources.

Confirmed on 2026-11-09 (`site/data/m7_quotes.json`). Each is printed in the journal layer as
it appears in Krajnović (2016), with the original reference he gives and "quoted in Krajnović
(2016)":
- § 3.2: Airy (1846), MNRAS 7, 121, on Le Verrier's letter to him of 28 June 1846.
- § 4.1: Galle to Le Verrier, 25 September 1846, with "Le [sic] planète" as printed, and
  Krajnović's translation.
- § 4.3: Airy to Le Verrier, private letter, 14 October 1846.

Arago's "au bout de sa plume" and Le Verrier's letter of 18 September 1846 are left off: their
wording is unverified.

## UI-rule audit (global rule), planned use

- **Small caps:** used only for "Plate", "Fig.", "Table" and "§". They are the period's
  labelling convention and mark cross-reference targets.
- **Section numbering:** used, because the prose cross-references sections.
- **Eyebrow:** the "§ 4.2" line above an h2 is the section number. It is the only eyebrow.
- **Italics:** used only for point letters, variables and titles of works, because the letters
  are the link between text and figure.
- **Gradient:** the radial vignette on plates is the cyanotype's uneven exposure, and it raises
  the bar for everything else. No other gradient, except a fade from a plate into paper if a
  confirmed quote band is added.
- **Several shades and opacities:** six palette colours and three line opacities, fixed above.
- **Grain, dust and cloud texture** (close to "dithering"): used on figure grounds only. They
  are the concept (the archive layer); text never carries them.
- **Bordered containers:** none. The hairline plate mark on Plate I is the plate's edge, not a
  box.
- **Periods:** "Plate II." and "Fig. 3." keep the period as the label convention. Titles have
  none.
- **Not used:** tags or pills, Instrument Serif, negative kerning, dividers (one hairline above
  the Table 1 totals only), monospace, ALL CAPS, code breadcrumbs, ASCII images, squircles,
  over-detailed diagrams (one idea per plate).

## Build

- Static, `site/` only. ES modules; an importmap pins D3 at an exact version from jsDelivr;
  fonts come from Google Fonts with serif fallbacks.
- **Prerendered for Pages (2026-11-09):** `tools/prerender_site.py` loads `site/` in headless
  Chrome and writes the rendered page (text, tables, captions, references, the static SVG of
  every plate and figure) as `_site/index.html`, with the head tags for link previews and
  `plate-i.png` (Plate I, 1200 × 630). Without JavaScript the page reads in full. With it,
  `main.js` redraws the same content (every render step replaces what it renders) and shows
  the controls that need JS (`data-js hidden`: Plate III's play and scrub, Copy). The Pages
  workflow runs the build, `tools/check_prerender.py` and `tools/check_layout.py`, and deploys
  `_site/`.
- **Text clear of marks (2026-11-10):** `site/js/clear.js`. After every draw, point letters take
  their drawn offset if it is clear, else the nearest clear spot around their point (Plate III's
  every frame). Lines that only guide the eye (gridlines, the launch line, Neptune's line in
  Fig. 6, Plate II's construction lines: class `yield`) break where they would cross a text.
  Labels along orbits are tried at a few positions and left out where none is clear; the
  caption carries the same words.
- **Kept from the current site:** `data.js` (loaders, the 95 % status, number formatting) and
  the per-panel error and pending states. Removed: the three.js scenes, `lanes.js` and the dark
  styles.
- **New data files:** `m7_table1.json`, `m7_motion.json` and `m7_plates.json` (the geometry of
  Plates I and II, with the posterior subsample; it takes the place of the
  `m7_conjunction_samples.json` planned here, so `m6c_conjunction.json` is not changed), plus
  `m7_references.json` (the reference list, so no cited number is typed into HTML). Each
  generated file is written by a script in `docs/notes/scripts/` and described in M7.md.

## Changes made while building (2026-11-09)

- **Plate II** is drawn with radii to scale and angles about *SC* magnified by a factor
  from `m7_plates.json` (the smallest of 1, 2, 5, 10, 20 that opens *CSL* to at least 10°:
  5). To scale, *C*, *N* and *L* fall within about 30 px of each other. The caption states
  the factor, and the angle label gives the true value.
- **Plate I** on wide screens (≥ 900 px) puts the title block at the bottom left and the
  figure on the right, so the figure can use the plate's full height. On phones the title sits
  below the figure, as specified.
- **Title:** "Refitting Le Verrier's prediction of Neptune", with a subtitle from the data
  ("A pre-registered Bayesian reconstruction from the 279 observations of Uranus he published
  in 1846").
- **Plates II and III** stay in view (sticky) while their text column scrolls past.

Review of the localhost page (2026-11-09):
- **Plate III's conjunction is *K*,** the conjunction of Uranus and Neptune (the run's, JPL's
  date to 0.2 h), not Plate II's *C* (the conjunction with the fitted planet, 3 October 1821).
  Fig. 3 uses *K* too. The 145.8 days in § 3.2 is computed from *K*'s date.
- **Split columns:** the in-text figure widening (`main > section > figure[data-fig]`) also
  caught the plates in Plates II and III, which made their tracks 707–1040 px wide. That pushed
  Plate III's side panel past the right edge at 768 px and cut Plate II's text column to
  104–288 px. The rule now skips `.split`. The two columns are `minmax(0, 1fr)` each, and
  below 900 px the plate and its text stack (at 900 px the narrower text column holds about
  47 characters).
- **Print** (`@media print`): one column; plates and figures full width at the aspect they were
  drawn at, at most 120 mm high (Plate I a full title page); captions under their figures;
  nothing sticky; no controls or copy buttons. Of the ground textures, print keeps the cloud
  and the vignette; mottle, grain and dust, rasterised at print resolution, take an A4 PDF from
  2.9 MB to over 20 MB.

## Acceptance

1. **Number audit:** no number on the rendered page outside the allowlist (the one rule).
2. **Layout** (`tools/check_layout.py`, run by the Pages workflow before deploy; phones as
   touch devices, and in WebKit as well as Chrome before a push): screenshots at 390, 430, 768,
   1024 and 1440 px for every plate and figure, saved to `docs/notes/m7_site/`. At each width
   (and at 360, the narrowest common phone, and at 899 and 900, either side of the stacking
   breakpoint):
   - no horizontal overflow;
   - no text box (SVG text, lettering along a path glyph by glyph, HTML text over a plate)
     intersects a stroke, a marker or another text, except its own leader line;
   - no label along a path is longer than its path (measured with its letter-spacing);
   - no emoji or emoji-presentation character in the rendered text, with or without JS (icons
     are inline SVG in the text colour);
   - nothing sticky or fixed below 900 px; blocks in the side panels of Plates II and III, and
     Plate III's controls, never overlap;
   - Table 1's headers and year ranges each on one line;
   - Plate II's and Plate III's side panels lie inside the viewport;
   - side by side, each text column holds at least 45 characters per line;
   - no vertical gap over 200 px between content blocks.
   In print (Save as PDF), from an 800 px window: nothing sticky, plates and figures at the full
   text width, every caption under its figure, controls hidden, and no text over a mark. The PDF
   is saved beside the screenshots.
3. **Interaction:** every point letter in the prose highlights its point. Plate III plays,
   pauses, scrubs and respects reduced motion. Tooltips open on tap as well as hover.
4. **Plate III:** the run with Neptune matches Horizons' Uranus in 1850 (the difference is
   reported), and the conjunction date agrees with `m6c_conjunction.json`'s JPL date (the
   difference is reported).
5. **Pending states:** while M8b and the arbiter have not reported, the provisional 95 % and
   the pending S1 verdict read as such.
6. **Guards:** no change to `src/`, `release/`, the pre-registered files or M8b; the
   design-ref file is never staged.
7. **Without JavaScript** (`tools/check_prerender.py`, run by the Pages workflow before
   deploy): the built page, loaded with JS off, passes the number audit, and its numbers, text
   and SVG labels match the JS-rendered page exactly (the JS-only controls left out of both).
   Title = the paper's title, description = the abstract's first sentence, canonical URL,
   Open Graph and Twitter card tags with `plate-i.png` at 1200 × 630.

## Figs. 7 and 8 (2026-11-10)

- **Fig. 7, the profile likelihood in a** (§ 5; beside Fig. 8 at the figure width in a two-column
  grid with a 32 px gap, one above the other below 900 px and in print): Δ ln L against the held
  a for all observations (2 px line) and for the observations from 1781 on (1 px line), dots at
  the grid values, the reading rule's threshold as a 0.6-opacity line among the 0.3 gridlines,
  *N* and *L* as open marks on the all-observations curve. The key is in the caption; neither
  curve is dashed, since dashed means provisional here. Tooltip: the nearest grid value's two
  readings.
- **Fig. 8, mass against a** (after the fact, and the caption says so): 600 seeded posterior
  samples as 0.6-opacity stipple, the 68 % and 95 % contours (2 px and 1 px), *N* and *L* as open
  marks. Tooltip: the marginal medians and 68 % ranges.
- Both read `site/data/m9_profile.json`; axis ticks are derived values. UI-rule audit: no new
  pattern. The figure pair is a grid, not bordered containers; no new colour, opacity step,
  gradient, divider, pill or eyebrow.

## Figs. 7 and 9 after M10 (2026-11-10)

- **Fig. 7 now has two panels** at the figure width (one above the other below 600 px of canvas):
  left the M9 profile on Gaussian errors, right the M10 profile on Student-t errors, each with
  the all-observations curve (2 px) over the from-1781-on curve (1 px), the dots, the threshold
  line, *N* and *L* as open marks (ids N and L on the left panel, Nt and Lt on the right, so the
  caption's letters light the left pair). One y axis, labelled on the left panel only; the
  panel titles are 14 px text above the plot area. Tooltip: the panel under the pointer and
  its nearest grid value.
- **Fig. 9, the leave-one-out profiles** (beside Fig. 8 in the two-column grid): four panels,
  one per pre-1781 observer removed, each the Student-t profile without that observer (2 px)
  over the all-observations curve (1 px). The y axis stops at −6 so that the −2 threshold and
  the shape near the maximum are legible in a 240 × 170 px panel; a point below −6 (the 26 AU
  dip) is drawn as an open circle on the bottom edge, and the caption says so.
- Both read `site/data/m10_profile.json` (Fig. 7's left panel still `m9_profile.json`); axis
  ticks and the −6 floor are derived values. The row-gap between panel rows is 68 px so that a
  panel title never meets the x labels of the row above (the layout check caught 50 px).
- UI-rule audit: no new pattern. Small multiples replace a four-curve overlay; no new colour,
  opacity step, dash, gradient, divider, pill or eyebrow; the panel titles are plain text, not
  eyebrows over headings.

## Companion page: the British search, Adams and Challis (2026-11-12)

Phase 1b gets its own page, `site/adams-challis/index.html`, served at `…/adams-challis/`. The
main page is not renumbered and gains one sentence: the last paragraph of its abstract ends
with a link to the companion study, its text read from `site/data/p1b_page.json`.

**Same system, second paper.** The concept, the one rule, the palette, the type, the figure
language, the voice and the interaction rules above apply unchanged. The page loads the same
`css/` and the shared modules of `js/` (`data.js`, `figs.js`, `plate.js`, `clear.js`); its own
text is in `js/copy_p1b.js` and its own figures in `js/figs_p1b.js`. It has no Plate I: the
title block sits on paper, with the title, the byline as on the main page, and one line
linking back to the main paper. Its sections, figures and tables are numbered from 1
(Figs. 1 to 4, Table 1); a cross-reference to the other page is written "Fig. N of the main
paper".

**Structure.**
1. Abstract: the two questions and their answers, each with its condition.
2. Sources: the printed record (Adams's memoir, Challis's reports, Le Verrier's two papers),
   the transcriptions and how each was checked, the frozen files and their hashes.
3. How good was Adams (four findings, one subsection each):
   - 3.1 the project's fit to his data at his distances falls short of Neptune (Fig. 1);
   - 3.2 his printed solutions reproduced from the memoir; his data through his method, and
     where the difference from his print sits;
   - 3.3 on a noise-free Neptune every method is short at his distance: the shortfall is the
     assumed distance (Fig. 1);
   - 3.4 noise against a systematic: neither is established (Fig. 2).
4. Challis's search as a detection model:
   - 4.1 the record, the model and the conditioning on the record (Fig. 3);
   - 4.2 the cases and their readings (Table 1, Fig. 4);
   - 4.3 what each reading rests on.
5. Limitations.
6. Data availability: the notes, the JSON files and the code paths; the DOIs as on the main
   page. Phase 1b is in the repository and not yet in a release, and the page says so.
7. References, the page's own list.

**Deviations stated where they bear on a claim** (each in the sentence that makes the claim,
not in a separate list):
- § 3.4: the pre-registered statistic treated an angle as linear and measured switching to
  the second minimum of his least-squares problem; the on-branch figures are post hoc and
  are labelled so wherever they appear.
- § 4.2: under the reading with the realised events held, his actual procedure is 0 before
  23 September by construction; the pre-registered paired rule is therefore uninformative,
  and the page says so in place of quoting its statements.
- § 4.2: three readings change row when the priors are halved and are marked
  prior-dependent in Table 1 and Fig. 4; the zone priors were not halved; the tornado used
  10,000 draws.
- § 4.2, one sentence: the replay reading (every event redrawn) is outside the
  pre-registered reading; its value for the actual procedure is given with its JSON path,
  and the note carries the rest.
- § 4.2, one sentence: paired differences between counterfactuals are post hoc and are in
  the note only; the sentence links to them.
- § 4.1: the Northumberland aperture per Airy 1844; the start places at ecliptic latitude 0;
  the zone order; Neptune's magnitude carried back from a modern fit.

**Table 1, the cases** (§ 4.2; journal layer, as Table 1 of the main paper): one row per case; columns: the
case, the probability before 23 September with its 5 to 95 % range, the row of the reading
rule, the input it is most sensitive to, what it is conditional on. Prior-dependent rows carry
a dagger defined in the caption. The two 7-day-delay variants are one sentence under the
table, not rows. Below 600 px the last column drops under each row's case name.

**Fig. 1, the shortfall follows the assumed distance** (§ 3.1, 3.3): one horizontal axis, the
heliocentric longitude minus Neptune's on Galle's night, Neptune a hairline at 0 lettered *N*.
One row per solution: his two printed solutions; his data through his method; the project's
fit at his two distances and with the distance free; each method on the noise-free Neptune at
his distance. A dot per value with its interval where the source has one (2 px line); rows
without an interval say so in the caption. Rows on the noise-free truth are open circles.
Tooltip: the row's value, interval and source step.

**Fig. 2, his shift against noise** (§ 3.4): the distribution over noise draws of the shift
toward Neptune, on an angular axis from −180° to 180° so that the second minimum is a
separate cluster and not a tail. Two panels (Hypothesis I, Hypothesis II). In each, noise
model (a) is the white 2 px line and model (b) the lower-opacity line, both as stepped
outlines of the same bins; the observed shift is a vertical hairline lettered and defined in
the caption. The caption states that the on-branch fraction is post hoc and gives the
pre-registered value beside it. The binned shifts are written from the cached draws of
step 3c to a new file, `docs/notes/P1B_adams_noise_hist.json`, with its hash and a dated line
in the step 3c note; `P1B_adams_noise.json` is not edited.

**Fig. 3, the search on the sky** (§ 4.1): right ascension against declination, equinox of
date, right ascension increasing to the left. Neptune's track from 29 July to 1 October as a
dotted line with the four nights of the record lettered; the 30 July and 12 August zone as a
band at its posterior median setting with its 5 to 95 % range as hairlines; the 4 August
sector sweep as a bar; the five start places (the 325° column, S4, S5, Le Verrier June and
August) as open marks lettered and defined in the caption; the edge of the Berlin sheet he had
as two vertical hairlines, one per equinox, labelled along the line. To scale, and the caption
says so. Tooltip: the nearest mark's place and source.

**Fig. 4, the cases against the reading rule** (§ 4.2): one row per case in Table 1's order;
the mean as a dot, the 5 to 95 % range as a 2 px line; the rule's two thresholds as
0.6-opacity lines among 0.3 gridlines; prior-dependent rows carry the same dagger as Table 1
and their halved-prior range as a 1 px line under the main one. Cases at 0 in every draw are
a dot at 0 with no line, and the caption gives the reason for each (by construction, by
timing, off the sheet). Tooltip: the row's numbers and what it is conditional on.

**Data.** `fork.export.to_json` writes `site/data/p1b_adams.json`, `p1b_challis.json` and
`p1b_page.json` as subsets of the committed `docs/notes/P1B_*.json`; every value keeps the
path of its source field, and a test asserts each equals its source. Cited constants (the
printed elements, the aperture, the magnitude constant) carry their citations. The site
ledger is regenerated with the Phase 1b rows.

**Build and checks, both pages.** `tools/prerender_site.py` prerenders `index.html` and
`adams-challis/index.html` into `_site/`. `tools/check_prerender.py` runs the number audit
and the JS-off comparison on each page, with one allowlist per page built from the files that
page loads. `tools/check_layout.py` takes the page as an argument and runs every width on
both, in Chromium and in WebKit; screenshots of the companion page go to
`docs/notes/p1b_site/`. The Pages workflow runs all three on both pages before the upload,
and one artifact deploys both.

**Acceptance, in addition to the list above:**
- the main page's rendered text differs from today's by the one sentence and nothing else;
- every number on the companion page is in its allowlist, and the allowlist's sources are
  committed files;
- the words "post hoc", "by construction" and "prior-dependent" appear in the sentences
  listed above, checked by a test on the copy;
- no change to `release/`, to any pre-registered block or to `src/` outside
  `fork/export/`; nothing is released.

**UI-rule audit:** no new pattern. The page uses what the main page uses, for the same
reasons: small caps on "Fig.", "Table" and "§"; section numbers, because the text
cross-references them; the one eyebrow per section; italics for point letters and titles of
works; the vignette on figure grounds. New on this page and justified: a dagger on three
rows, because the mark must sit in a table cell and in a figure row alike. Not used: pills,
bordered containers, dividers beyond Table 1's one hairline, monospace, ALL CAPS, any new
colour, opacity step or gradient. The approved "one hue plus grey" is carried inside the
palette: white line on the cyanotype ground, with the lower line opacity in the role of grey.

**Changes made while building (2026-11-12).**
- Order of the figures, by appearance: Fig. 1 the shortfall, Fig. 2 the noise draws (two
  panels, the first noise model as the 2 px white line and the second at the lower opacity),
  Fig. 3 the sky, Fig. 4 the cases.
- Fig. 1 has no row for the fit with the distance free: the committed JSON holds the true
  longitude on Galle's night for the held-distance fits only. The sentence of § 3.1 gives the
  free-distance result as the JSON states it (Neptune inside the 68 % range).
- Fig. 3 draws the zone of 30 July and 12 August at the middle of the settings the record
  allows, with hairlines at the two extremes; the posterior of the zone's setting is not in a
  committed JSON. The 4 August sector is not drawn: its centre is not in the record.
- Fig. 2 reads `docs/notes/P1B_adams_noise_hist.json`, binned from the cached draws (5° bins on
  the circle); `P1B_adams_noise.json` is not edited.
- Table 1's last column gives each case's condition in the page's words; the note's wording
  stays in `p1b_challis.json`. On phones the condition sits under the case, and Figs. 1 and 4
  use shorter row labels.
- Table 1 omits the two one-week-delay variants (identical to their cases to two decimals);
  the table's footnote says so.
- References give title, journal, volume and page as the notes record them. The titles and
  page ranges of three short communications were added from the publisher's records: Hind 1850
  and Schönfeld 1886 from the Crossref records of *Astronomische Nachrichten*
  (10.1002/asna.18500301706, 10.1002/asna.18861140406), Mauvais 1847 from the heading on
  p. 666 of the scanned volume (archive.org, `comptesrendusheb24acad`). ADS, the intended
  source, did not serve its records to a script.
- The companion page's link preview uses the main page's Plate I image.
- Screenshots (390 and 1440 px, Chromium and WebKit) are written to `docs/notes/p1b_site/`
  by `tools/check_layout.py --shots` and are not tracked: the folder is in `.gitignore`.
- The Pages workflow file is unchanged: its three commands now cover both pages.
