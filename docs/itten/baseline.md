# ITTEN baseline

ITTEN version: v0.1.0 (c1962e2)
Date: 2026-11-15
Branch: `itten-fixes`, from `itten-pass` at 3ce0015

The site as it stood before any ITTEN change, inspected without a system file (distinct-value
mode), so the counts are what the page uses, not what a scale allows. Reports and crops are in
`itten-output/` (ignored).

## Run

```sh
cd site && python3 -m http.server 8791 --bind 127.0.0.1
node itten/inspector/inspect.js http://127.0.0.1:8791/ --full --out itten-output/baseline-main
node itten/inspector/inspect.js http://127.0.0.1:8791/adams-challis/ --full --out itten-output/baseline-p1b
```

- Layout sweep: every 40 px from 320 to 1920 px (41 widths), 900 px tall.
- Contrast and fonts at 320, 768 and 1440 px, light and dark. Text doubled (browser default
  32 px). Device pixel ratio 2 at 320 and 1440 px.
- Engines: Chromium 159.0.8154.6 and WebKit 27.2 at 320, 768, 1440 and 1920 px. Firefox did not
  start on this machine ("Could not find profile folder"), so there is no Firefox opinion.

## Counts

| | Main page | Companion page |
|---|---|---|
| Errors | 0 | 0 |
| Warnings | 12 | 10 |
| Info | 2 | 0 |
| Distinct spacing values at one width | 22 to 24 | 16 to 17 |
| Distinct font sizes at one width | 6 | 4 to 6 |
| Distinct spacing values over the sweep | 68 | 49 |
| Distinct font sizes over the sweep | 48 | 40 |
| Fixed spacing values (px) | 4, 6, 8, 9, 10, 12, 14, 15, 16, 22, 24, 28, 32, 36, 48, 64, 88, 96 | 6, 9, 10, 12, 21, 22, 24, 36, 48, 72, 81, 96 |
| Fixed font sizes (px) | 14, 15 | 14, 15 |
| Radii | none | none |
| SVG text sizes (not judged) | 12.5, 13, 14, 19 | 12.5, 14, 19 |
| Text that ignores the browser text size | 234 of 234 | 154 of 154 |
| Elements / text elements measured | 315 / 234 | 187 / 154 |

The sweep counts are high because the body, h2, Plate III year and companion h1 sizes follow
the window (`clamp()` with `vw`), and every margin set in `em` follows them.

## Findings

### Main page (`/`)

| ID | Check | Severity | What |
|---|---|---|---|
| i013 | font-fallback | warning | Play button and range input computed in Arial 13.33 px, not Source Serif 4 |
| i003 | spacing-scale | warning | Paragraph `margin-bottom: 1em`: 16.5 to 18 px, 11 elements |
| i004 | spacing-scale | warning | `h3 .sc` `margin-right: 0.4em`: 6.6 to 7.2 px, 7 elements |
| i005 | spacing-scale | warning | `.quote` `padding-left: 1.6em`: 26.4 to 28.8 px, 3 elements |
| i006 | spacing-scale | warning | Plates II and III text column side padding `var(--gutter)`: 16 to 56 px, 26 values |
| i007 | spacing-scale | warning | `.limits` `padding-left: 1.1em`: 18.15 to 19.8 px |
| i008 | spacing-scale | warning | `.limits li` `margin-bottom: 0.7em`: 11.55 to 12.6 px, 10 elements |
| i009 | type-scale | warning | h2: 30 to 42 px, undeclared fluid range |
| i010 | type-scale | warning | Body: 16.5 to 18 px, undeclared fluid range, 82 elements |
| i011 | type-scale | warning | Quote text `1.06em`: 17.49 to 19.08 px |
| i012 | type-scale | warning | Plate III year: 64 to 104 px, undeclared fluid range |
| i014 | text-resize | warning | 234 of 234 text elements keep their size when the browser text size is doubled |
| i001 | hairline-subpixel | info | Play button border at y = 8943.2 px at 1× (could not be profiled) |
| i002 | hairline-subpixel | info | Copy buttons' borders at fractional positions; both engines snap them crisp |

### Companion page (`/adams-challis/`)

| ID | Check | Severity | What |
|---|---|---|---|
| i001 | spacing-scale | warning | h1 `margin-bottom: 0.3em`: 9.6 to 15 px |
| i002 | spacing-scale | warning | Paragraph `margin-bottom: 1em`: 16.5 to 18 px, 16 elements |
| i003 | spacing-scale | warning | `h3 .sc` `margin-right: 0.4em`: 6.6 to 7.2 px, 8 elements |
| i004 | spacing-scale | warning | `.limits` `padding-left: 1.1em` |
| i005 | spacing-scale | warning | `.limits li` `margin-bottom: 0.7em` |
| i006 | type-scale | warning | h1: 32 to 50 px, undeclared fluid range |
| i007 | type-scale | warning | Byline: 14 to 17 px, undeclared fluid range |
| i008 | type-scale | warning | h2: 30 to 42 px |
| i009 | type-scale | warning | Body: 16.5 to 18 px |
| i010 | text-resize | warning | 154 of 154 text elements keep their size when the browser text size is doubled |

No contrast, overflow, clipping, overlap, target-size, near-miss-edge, sibling-gaps or radius
finding on either page, in Chromium or WebKit.

## Not seen by the inspector

Everything inside a `<figure>` except its `<figcaption>` is skipped (21 figure parts on the main
page, 6 on the companion). On this site that hides more than drawings:

- **Table 1** (`figure.table`): its cells are `clamp(16px, 0.6vw + 12px, 20px)` on the main page
  (19.68 px at 1280) and `clamp(15px, 0.4vw + 12px, 17px)` on the companion; cell padding 5, 6,
  7 and 10 px; footer padding 2 and 8 px.
- **Plate I's title block** (`figure.plate`): h1 `clamp(40px, 6.4vw, 100px)`, then
  `clamp(44px, 5.6vw, 92px)` from 900 px (80.64 px at 1440); subtitle and byline `clamp()` in
  px and vw; h1 `margin-bottom: 0.22em` (8.8 to 17.74 px).

The teardown (`docs/teardown.md`) measures these with a separate computed-style probe.
