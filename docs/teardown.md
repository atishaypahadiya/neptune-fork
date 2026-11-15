# Teardown: neptune-fork site

ITTEN version: v0.1.0 (c1962e2)
Date: 2026-11-15

The cyanotype direction is approved (`docs/DESIGN.md`). This teardown does not reopen it. It
measures the reference the direction came from and the site as built, to find the closed
system the site is already trying to be, and which of its values are that system and which are
drift.

## References

| ID | Reference | Kind | Scale | Uncertainty | Notes |
|---|---|---|---|---|---|
| R1 | `design-ref/neptune-fork site directions.html`, "Cyanotype" page (local only, git-excluded, never committed) | live page | CSS px | exact | Four 1280 × 900 artboards in iframes: Plate I, Plate II, Table 1, Plate III. Each frame was rendered in Chromium, saved as a static snapshot (runtime scripts removed) and read at 1280 px. Its numbers and words are placeholders; only its geometry and type are measured. |
| R2 | The site as built: `site/index.html` (R2a) and `site/adams-challis/` (R2b), served locally | live site | CSS px | exact | Computed styles at 390, 1280 and 1440 px; inspector in distinct-value mode over 320 to 1920 px (`docs/itten/baseline.md`). |

Selectors for R1 are given as the element and its place in the artboard, because its markup
has no classes (inline styles on anonymous `div`s).

## Findings

### Type

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| Body | 18 px / 27.9 px (1.55), Source Serif 4 400 | measured | R1 Plate II, text column `p` @1280 |
| Section mark "§ 4.2" | 17 px, tracking 1.5 px | measured | R1 Plate II, first `p` @1280 |
| h2 | 40 px / 42 px, GFS Didot 400 | measured | R1 Plate II `h2` @1280 |
| Caption | 14 px / 21 px (1.5); "Plate II." label tracked 1 px | measured | R1 Plate II, last `p` @1280 |
| Table 1 cells | 20 px; source names 22 px; header 15 px tracked 1 px; title 16 px | measured | R1 Table 1 grid rows @1280 |
| Plate III year | 96 px / 86.4 px (0.9), GFS Didot | measured | R1 Plate III `p > span` @1280 |
| Plate III text column | 16/24, 14/20.3, Fig. 3 caption 13/18.85, control label 12 px | measured | R1 Plate III text column @1280 |
| Title (Plate I) | 112 px GFS Didot (SVG text) | measured | R1 Plate I `svg text` @1280 |
| Figure lettering | orbit labels 13 px tracked 1.2 px; point letters 16 to 20 px italic; plate label 18 px tracked 2 px; angle 15 px italic | measured | R1 Plates I to III, `svg text` @1280 |
| Steps on R1's text, ÷ body 18 | 0.72 (13), 0.78 (14), 0.83 (15), 0.94 (17), 1, 1.11 (20), 2.22 (40), 5.3 (96), 6.2 (112) | measured | R1 |
| Pattern | a tight cluster of text sizes from 13 to 20 px, then a jump of 2.2× to h2 and over 5× to display numerals | measured | R1 |
| Body (site) | `clamp(16.5px, 1.05vw + 10px, 18px)` / 1.56: 16.5/25.74 to 619 px, 18/28.08 from 762 px | measured | R2a `main > section > p` @390, @1440 |
| h2 (site) | `clamp(30px, 2.4vw + 12px, 42px)` / 1.05: 30 px to 750 px, 42 px from 1250 px | measured | R2a `h2` sweep |
| h3 (site) | body size, 600, line height 1.4 (25.2 px at 18) | measured | R2a `h3` @1440 |
| Plate III year (site) | `clamp(64px, 6vw, 104px)` / 1: 64 px to 1067 px, 86.4 px at 1440 | measured | R2a `#h-plate-iii` sweep |
| Plate I title (site) | `clamp(40px, 6.4vw, 100px)` below 900 px, `clamp(44px, 5.6vw, 92px)` from 900 px; 80.64 px at 1440 | measured | R2a `.titleblock h1` @390, @1440 |
| Subtitle, byline (site) | `clamp(15px, 1.25vw, 19px)`, `clamp(14px, 1.2vw, 17px)` | measured | R2a `.subtitle`, `.byline` |
| Companion title (site) | `clamp(32px, 3.6vw, 50px)` / 1.08 | measured | R2b `.paper-title h1` sweep |
| Table cells (site) | main `clamp(16px, 0.6vw + 12px, 20px)`; companion `clamp(15px, 0.4vw + 12px, 17px)` | measured | R2a, R2b `td` @390, @1440 |
| Fixed text sizes (site) | 14 px: captions, table headers, copy button, tooltip (74 elements at 390); 15 px: section marks, notes, references, citation (68) | measured | R2a @390 |
| Quote text (site) | `1.06em` of body: 17.49 to 19.08 px | measured | R2a `.q-text` sweep |
| Text sizes per width (site) | 6 on the main page, 4 to 6 on the companion | measured | `docs/itten/baseline.md` |
| Text that follows the browser text size (site) | none: 234 of 234 and 154 of 154 text elements keep their size at 2× | measured | `docs/itten/baseline.md` |
| Line heights (site) | one ratio, 1.56, inherited by every size that does not set its own, so 14 px gets 21.84 in table headers and 21 in captions, 15 px gets 23.4 | measured | R2a `th`, `figcaption`, `.sec` @1440 |

Of the six sizes the site uses at one width, four are intentional and two are drift:

- **Intentional:** caption 14 (R1 14), body 18 (R1 18), h2 42 (R1 40, enlarged on the site on
  purpose, `docs/DESIGN.md` "Type"), the Plate III year (R1 96).
- **15 px** is not in R1's text column; R1 uses 15 only for Table 1's header. On the site it is
  the size of section marks, notes and references, 68 elements, a deliberate second small step.
  It sits 1 px from the caption step (ratio 1.07), the kind of near-duplicate TYP-06 warns
  about. It is kept: merging it into 14 or into body moves 68 elements' size and is the kind of
  change this pass does not make.
- **Drift:** the quote at `1.06em` (19.08 px beside an 18 px body, a 6 % difference that reads
  as neither the same nor different), and the body's fractional phone size, 16.5 px.

### Spacing

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| Text column gap, R1 | 18 px between every block: section mark, h2, paragraphs (5 gaps) | measured | R1 Plate II text column (`display: flex; gap: 18px`) @1280 |
| Text column padding, R1 | 80 px top and left, 64 px right and bottom | measured | R1 Plate II text column @1280 |
| Table 1 rhythm, R1 | title to table 18 px; row padding 7 px (×7 rows); column gap 16 px; header padding 6 px; total row padding 22 px; 48 px below | measured | R1 Table 1 @1280 |
| Plate III text column, R1 | starts 64 px from the top; gaps 14, 10 and 4 px | measured | R1 Plate III @1280 |
| Distinct spacing values per width (site) | 22 to 24 (main), 16 to 17 (companion) | measured | `docs/itten/baseline.md` |
| Fixed spacing values (site) | 4, 6, 8, 9, 10, 12, 14, 15, 16, 22, 24, 28, 32, 36, 48, 64, 88, 96 (main); 6, 9, 10, 12, 21, 22, 24, 36, 48, 72, 81, 96 (companion) | measured | same |
| Paragraph gap (site) | `1em`: 18 px at 1440 (12 elements main, 16 companion), 16.5 at 390 | measured | R2a, R2b `section > p` |
| Section mark to h2 (site) | 6 px (`.sec` margin-bottom), 15 elements | measured | R2a, R2b `.sec` |
| h2 to text | 22 px | measured | R2a `h2` @1440 |
| h3 | 48 px above, 12 below | measured | R2a `h3` @1440 |
| Figure margins | 36 px above and below (`--fig-space`), 10 figures on each page | measured | R2a, R2b `main > section > .wide` |
| Caption under a figure | 10 px (figures), 12 px (tables) | measured | R2a `figcaption`, `.table > figcaption` |
| Space above a section mark | 88 px main, 72 px companion, 64 px on phones | measured | R2a, R2b `.sec` |
| Space around a half-bleed plate | 96 px (64 px on phones); main page ends with 96 px | measured | R2a `.split`, `main` |
| Plates II and III text column | 64 px top and bottom; sides `var(--gutter)`; 28 px top on phones, 32 px from 721 to 899 px | measured | R2a `.split > .side` |
| Values from `em` (site) | `.4em` (6.6 to 7.2), `.5em` (7.5 to 9), `.22em` (8.8 to 17.7), `.3em` (9.6 to 15), `.6em` (9), `.7em` (11.55 to 12.6), `1em` (15, 16.5 to 18), `1.1em` (18.15 to 19.8), `1.4em` (21), `1.6em` (24, 26.4 to 28.8) | measured | R2a, R2b sweep |
| Values from the window (site) | gutter `clamp(16px, 4vw, 56px)`; companion title `clamp(48px, 9vh, 112px)` (81 px in a 900 px window); Plate I title block `clamp(28px, 6vh, 64px)` | measured | R2a, R2b |
| Between-group ÷ within-group (site) | figure 36 ÷ paragraph 18 = 2.0; h3 48 ÷ 18 = 2.7; section 88 ÷ 18 = 4.9; caption 10 under its figure against 36 to the next block = 3.6 | measured | R2a @1440 |

**The scale the site is trying to be.** Eight of the site's values are exact multiples of 6 and
form a geometric series: 6, 12, 18, 24, 36, 48, 72, 96 (each step 1.5 or 1.33 times the one
below). They are the values the page's structure is built from: the section mark's 6, h3's
12 below, the paragraph's 18 (= the body size, and R1's text-column gap), 24, the figures'
36, h3's 48 above, the companion's 72 and the plates' 96. 6 is a third of the 18 px body.

The other values are drift from that series:

| Value | Where | Why it is drift |
|---|---|---|
| 2, 4, 5, 7, 8, 9, 10 | table cell and footer padding, citation row gap, list items, controls | single uses, or two values 1 to 2 px apart doing one job (table rows 5 px on the main page, 7 px on the companion) |
| 14, 15, 16 | copy button, Plate I caption, `.note`, controls gap, citation column gap | three values within 2 px for "a small gap" |
| 21, 22 | companion "back" line, h2 to text | 2 to 3 px off 24 |
| 28 | quote and controls bottom, phone text column | 4 px off 24 |
| 32 | main Table 1 margin, Figs. 8 and 9 gap, tablet text column | 4 px off the figures' 36, beside them on the same page; the companion's Table 1 already takes 36 |
| 64, 88 | phone section space, plate text column padding; main section space | each 8 px under a series step (72, 96) |
| 81 and the other window-height values | companion title, Plate I title block | follow the window height, so no two readers get the same gap |
| every `em` value | listed above | fractional, and they change with the type they sit beside |

### Edges and grid

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| R1 text column | 80 px from the left edge, 416 px wide, beside a 720 px plate (560 + 720 = 1280) | measured | R1 Plate II @1280 |
| R1 Plate III column | x = 880, 340 px wide, 60 px from the right edge | measured | R1 Plate III @1280 |
| Site grid | named lines full, wide, body; body = `min(40rem, 100% − 2 × gutter)`, wide = body + up to 12.5rem each side, gutter `clamp(16px, 4vw, 56px)` | measured | R2a `main` grid, `docs/DESIGN.md` "Widths" |
| Site lines at 1440 | body 400 to 1040, wide 200 to 1240 | measured | R2a @1440 (and DESIGN.md table) |
| Site lines at 390 | body = wide 16 to 374 | measured | R2a @390 |
| Plates II and III | two equal columns from 900 px; text measure 34rem inside the gutter; stacked below 900 px with the text on the body lines | measured | R2a `.split` |
| Near-miss edges | none on either page, Chromium and WebKit, 320 to 1920 px | measured | `docs/itten/baseline.md` |

### Density

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| Characters per line, R1 text column | 49 (mean of 8 lines, 44 to 52) | measured | R1 Plate II `p` @1280 |
| Characters per line, site prose | 75 at 1440 (10 lines, 66 to 81); 44 at 390 (38 to 50) | measured | R2a `main > section > p` |
| Characters per line, site plate text | 66 at 1440, 51 at 1024 | measured | R2a `.split .side p` |
| Plates per screen | one: Plates II and III are sticky and as tall as the window to 900 px | measured | R2a `.split > .plate` |

### Corners

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| Radii in R1 | none, except the Play button, a 44 px circle (radius 22) | measured | R1 Plate III `button` @1280 |
| Radii on the site | none, except the Play button, `border-radius: 50%` on 44 × 44 | measured | R2a `.play` |
| Nesting | no rounded box inside another | measured | R1, R2 |

### Colour

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| Plate ground | #1C3F94 (Plates I, II); #0F2A6B (Plate III) | measured (computed) | R1 artboard `div` backgrounds @1280 |
| Paper | #EEF1F3 | measured | R1 Plate II text column, Table 1 |
| Line on plates | #F4F7FB; Play button border at 70 % | measured | R1 `svg`, Plate III `button` |
| Text on paper, R1 | body #13235A (13.1 : 1 on paper); captions and headers #3B4A78 (7.6 : 1); table values #1C3F94 | measured | R1 Plate II, Table 1 |
| Text on paper, site | body `--deep` #0F2A6B (11.8 : 1); captions, labels, values, links `--prussian` #1C3F94 (8.5 : 1) | measured | R2a, palette in `docs/DESIGN.md` |
| Line on plates, site | `--line` #F4F7FB on `--prussian` 8.9 : 1, on `--deep` 12.5 : 1, on `--mid` 5.5 : 1 | measured | R2a |
| `--prussian` roles | figure ground, secondary text, links, Table 1 values, control borders, focus ring | measured | R2a `site.css` |
| Contrast findings | none on either page, light and dark scheme (the site is light only) | measured | `docs/itten/baseline.md` |

`--prussian` has several roles (COL-05). They are one family on purpose: the period's blue
print, where the plate's ink and the printed labels are the same blue. It is recorded, not
changed.

### Texture

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| Grain, mottle, cloud, dust, vignette | SVG filters on plate grounds only; no texture on paper or behind paper text | measured | R1, R2a `svg.ground` |

### Deliberately broken

| Finding | Value | Measured / estimated | Source |
|---|---|---|---|
| R1's controls in the default control font | the Play label renders in Times | measured | R1 Plate III, inspector `font-fallback` |
| R1's script lettering over the plate | Monsieur La Doulaise, 58 px | measured | R1 Plate I |
| Site's controls in the default control font | Play button and range input compute Arial 13.33 px (CRF-13) | measured | R2a, baseline i013 |

## Carry into the system

- The 6-based series 6, 12, 18, 24, 36, 48, 72, 96 (R2 spacing; 18 also R1's text-column gap)
  → `spacing.scale`.
- Figure ÷ paragraph 2.0 (R2a) → `spacing.groupRatioMin`: 2.
- Body 18/28 (R1 18/27.9, `docs/DESIGN.md` 18/28) and its phone size, 17/26 → `type.scale`
  `body` and `body-phone`, switched at 720 px, the site's phone breakpoint.
- h3, the body size at weight 600 (R2a) → `h3` 18/28 and `h3-phone` 17/26, switched with the body
  steps. Separate steps because ITTEN 0.2 matches a step on size, weight and family together.
- Caption 14/21 (R1, R2) → `caption`. The site's 15 px second small step → `note`, 15/23.
- h2, Plate III year, the two Plate I title ramps, companion title, subtitle, byline and the
  two Table 1 cell ramps (R2) → fluid steps, each the range the page renders now, with the
  limits in rem.
- One radius: the Play button's circle (R1, R2) → `radius.capsule`: true, no fixed radii.
- Palette roles (R1, `docs/DESIGN.md`) → `colour.roles`; pairs from the measured ratios.
- The glow's 160 ms fade on point highlight (R2a) → `motion`.

## Deliberately not carried

- R1's body and caption colours (#13235A, #3B4A78): `docs/DESIGN.md` fixed the palette at six
  colours; the site's `--deep` and `--prussian` pass at 11.8 and 8.5 : 1.
- R1's 80/64 text-column padding and its 7 and 22 px table padding: R1 is a 1280 px artboard,
  and the site's own repeated values cover the same jobs.
- R1's 112 px title and 96 px year: the site scales them with the window.
- R1's script lettering and default control font: the first is retired by `docs/DESIGN.md`, the
  second is a defect (CRF-13).
- The fluid gutter as a spacing step: it is the page's side padding and the outer grid track
  (`docs/DESIGN.md`, "Widths"), it reaches the window's edges at every width by construction,
  and stepping it would move the text edges at every width from 400 to 1400 px. It stays fluid
  and is logged in `docs/itten-overrides.md`.
- Figure heights (230 to 700 px) and SVG lettering inside figures (12.5, 13, 14, 19 px): drawing
  geometry, held to `tools/check_layout.py`'s text-over-marks check, not to the text scale.
- Print styles (pt): a separate medium, not inspected, unchanged.
