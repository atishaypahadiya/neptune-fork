# ITTEN 0.2 pass

ITTEN version: v0.2.0 (d998de0). System: `itten.system.json` version 1. Date: 2026-11-15.
Branch `itten-0.2`.

Every run is `--full` (320 to 1920 px every 40 px, 200 % text, light and dark, DPR 2, the
motion pass), with `--system itten.system.json`. Engines: Chromium 159.0.8154.6 for the sweep
and WebKit 27.2 for the layout checks at 320, 768, 1440 and 1920 px, with no finding that
differs between them. Firefox still does not start on this machine ("Could not find profile
folder"). Before is `2aee437` (the pin), served from a copy of that commit; after is `773928a`
plus the override markers committed with this note. Both served by a small Node static server.
`tools/check_prerender.py` and `tools/check_layout.py` (Chromium and WebKit) pass after.

## Counts

| Check | Main before | Main after | Companion before | Companion after |
|---|---|---|---|---|
| Errors | 0 | 0 | 0 | 0 |
| Warnings | 28 | 0 | 20 | 0 |
| `rag` | 23 | 0 (1 logged) | 16 | 0 |
| `widow` | 2 | 0 | 2 | 0 |
| `type-scale` | 2 | 0 | 2 | 0 |
| `spacing-scale` | 1 | 0 | 0 | 0 |
| `near-miss-edge` | 0 | 0 (2 logged) | 0 | 0 (1 logged) |
| Info | 59 | 58 | 52 | 48 |

Info: `text-wrap-pretty` 50 → 49 and 42 → 38; `text-wrap-balance` (6, 5), `optical-margin`
(1, 5) and `hairline-subpixel` (2, main) unchanged. At the default four widths the before
runs had 12 and 10 rag warnings (22); the full sweep finds 39.

"Logged" warnings are hidden by `data-itten-ignore` and have rows in
`docs/itten-overrides.md`. A run with those markers removed (`itten-output/v02-noignore-*`)
gives exactly those four: three `near-miss-edge`, an inspector artefact (`upstream.md`,
issue 4), and the rag at 920 px in Plate II's text column, above the phone breakpoint.

## What changed

1. `itten.system.json`: steps `h3` 18/28 and `h3-phone` 17/26 at weight 600 (ITTEN 0.2
   matches size, weight and family together).
2. `itten.system.json`: `spacing.fluid` "gutter", 16 to 56 px between 400 and 1400 px wide;
   its overrides row went.
3. `site/css/site.css`: on screens 720 px wide and narrower, `hyphens: auto` on text in
   `main`, with `hyphenate-limit-chars: 6 3 3` and (WebKit) `hyphenate-limit-lines: 2`.
   Headings, links, tables and Plate III's quote keep `hyphens: manual`. `lang="en"` was already
   on both pages.
4. `site/js/plate.js`, `keepWhole()`: at phone widths, URLs, DOIs, arXiv ids, paths, dotted field
   names and capitalised words go in a `span.whole` (`hyphens: manual`). Above 720 px the spans
   are removed, because WebKit does not kern across a span's edge and they moved line breaks at
   desktop widths.
5. Widows: a no-break space in "not pre-registered." (Plate II, twice) and in the translator
   credit "Krajnović (2016)", which hyphenation had left as "(2016))" alone at 360 px.

### Hyphenation, checked word by word

A script listed every word broken at a line end, both pages, Chromium and WebKit, at 320, 360,
390, 414, 600, 720, 721 and 768 px.

- No DOI, URL, arXiv id or path is hyphenated. Long DOIs and the longest path still break
  without a hyphen at 320 to 360 px, as before (`overflow-wrap`).
- Chromium does not hyphenate capitalised words. WebKit did (Nep-tune, Ver-rier's, Cam-bridge,
  Au-gust) until `keepWhole()`; after it, none.
- Without limits WebKit broke two-letter fragments (un-der, with-in); `6 3 3` stops that.
- Hyphenation narrowed Table 1's columns on the companion at 320 px, so protected names broke
  without a hyphen; tables are now exempt and unchanged.
- French and German titles in the references are hyphenated with English patterns (for example
  "observa-tion" in Mauvais 1847); the breaks seen are acceptable in either language.
- Without JavaScript (the prerendered page as a crawler or reader mode gets it) there are no
  spans, so WebKit on a phone would hyphenate names there.

## Visible changes

Text blocks whose line breaks differ, by width (Chromium / WebKit):

| Width | Main | Companion |
|---|---|---|
| 320 | 39 / 40 | 29 / 32 |
| 360 | 42 / 44 | 37 / 35 |
| 400 | 39 / 41 | 24 / 26 |
| 440 | 37 / 42 | 26 / 26 |
| 480 | 31 / 33 | 23 / 25 |
| 520 | 29 / 28 | 21 / 21 |
| 560 | 28 / 30 | 23 / 24 |
| 600 | 28 / 27 | 17 / 19 |
| 640 | 26 / 28 | 19 / 19 |
| 680 | 26 / 28 | 17 / 17 |
| 720 | 24 / 25 | 16 / 18 |
| 760 to 1920 | 1 / 1 at 920, 1160 and 1280 up; 0 elsewhere | 0 / 0 |

Above 720 px the one change is Plate II's "not pre-registered." moving down a line. Nothing
else changes: no colour, size, spacing or position, and the same blocks in the same order.

Crops: `itten-output/compare-0.2/` (ignored), 229 pairs, before left and after right, every
changed block at 320 and 720 px in both engines and the widow at 920, 1160 and 1440 px;
`linediff.json` there has the lines of each block before and after.
