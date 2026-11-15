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

## Second pass: spans with the text, hyphenation by column

Same ITTEN, engines and sweep. Before is `f2766e9` (the first pass); after is `dfb690b` plus the
override markers committed with this section.

1. `4ed9b53`: `keepWhole()` is gone. `appendText()` (`site/js/plate.js`) puts names and
   addresses in `span.whole` as `fill()`, the reference lists and the quotes write their text,
   so the prerender serializes them: 784 spans on the main page and 543 on the companion, the
   same with JavaScript on and off, and nothing is wrapped after load. WebKit does not kern
   across a span's edge, so each span takes what kerns with its neighbours: the punctuation
   around it, a trailing space, an elided prefix (l'Académie) and hyphenated initials (S.-F.).
   A canvas check in WebKit over the whole text finds 0 of 2452 span edges where the pair
   would kern. Line breaks above 720 px are unchanged in both engines; at 720 px and below
   0 to 4 blocks per width move by a word.
2. `dfb690b`: Plates II and III's text columns are inline-size containers and hyphenate under
   26em of their own width, 900 to about 1110 px. The main column has no box of its own (it is
   a grid track), and an inline-size container on `main`, `html` or a section changed WebKit's
   automatic column widths in the companion's Tables 1 and 2 at every width, so it keeps a
   media query, now on its measure: under 34em of the phone body, a window of 39.25rem
   (628 px). Both engines switch at the same widths (628 on, 632 off; 1110 on, 1120 off).

### Counts

| Check | Main first pass | Main now | Companion first pass | Companion now |
|---|---|---|---|---|
| Errors | 0 | 0 | 0 | 0 |
| Warnings | 0 | 0 | 0 | 0 |
| `rag` logged | 1 | 0 | 0 | 0 |
| `near-miss-edge` logged | 2 | 2 | 1 | 2 |
| `widow` logged | 0 | 0 | 0 | 2 |
| Info | 58 | 62 | 48 | 49 |

Info: `text-wrap-pretty` 49 → 53 and 38 → 39 (short last lines at 640 to 720 px, where text is
no longer hyphenated); the rest unchanged. No finding differs between Chromium and WebKit.

Without the markers (`itten-output/v03-noignore-*`) the runs give exactly the logged six:

- The rag at 920 px in Plate II is gone, and its marker and overrides row with it.
- `near-miss-edge`, all the inspector artefact of `upstream.md` issue 4, now that spans are in
  the text at every width: main, the citations at 400 px and the references at 520 px;
  companion, the references at 480 to 520 px. Measured over all inline content the edges are
  identical (16 and 16, 45 and 45, 43, 43 and 43 px). The quote-source pairs of the first pass
  no longer appear, so their marker in `main.js` went; the citation blocks gained one.
- `widow` on the companion: "arXiv:2108.06305." (360 px) and the path
  "docs/notes/P1B_challis_detection.json." (360 to 640 px). Both lines are as they were in
  the first pass, and both warnings were in the run before it. The first pass hid them by
  accident: ITTEN splits words at element edges, and its spans ended before the full stop, so
  the last line read as two words. Logged, with `upstream.md` issue 5 (an identifier as a
  widow).

### Hyphenation, checked word by word

The same script, both pages and engines, at 320, 360, 600, 628, 632, 720, 900, 920, 1000, 1080
and 1110 px. No name, DOI, URL, arXiv id or path is hyphenated; the long DOIs, the repository
URL and the longest path still break without a hyphen at 320 to 360 px (`overflow-wrap`), as
before. Nothing is hyphenated at 632 or 720 px. In Plates II and III, 900 to 1110 px, only
common words break (pos-terior, helio-centric, con-junction), at most two lines in a row. The
prerendered page now carries the spans, so WebKit keeps names whole without JavaScript too.

Not changed, seen in passing: in WebKit, Table 2 on the companion breaks words without a
hyphen in its narrow columns from 900 to 1110 px ("Bouvard|'s", "referenc|e"). It does the same
at `f2766e9`; the tables are kept out of hyphenation and their lines did not move.

### Visible changes

Text blocks whose line breaks differ from the first pass (Chromium / WebKit):

| Width | Main | Companion |
|---|---|---|
| 320 to 600 | 0 to 1 / 0 to 4 | 0 to 1 / 0 to 2 |
| 640 | 25 / 27 | 19 / 19 |
| 680 | 26 / 28 | 17 / 17 |
| 720 | 24 / 25 | 16 / 18 |
| 760 to 880 | 0 / 0 | 0 / 0 |
| 920 to 1080 | 3 to 4 / 3 to 4 | 0 / 0 |
| 1120 to 1920 | 0 / 0 | 0 / 0 |

At 640 to 720 px the text is back to its lines before ITTEN 0.2 (unhyphenated). From 920 to
1080 px the changes are Plate II's paragraphs, Fig. 3's caption and the captions of Plates II
and III, all in the side column. Nothing else changes: no colour, size, spacing or position.
`tools/check_prerender.py` and `tools/check_layout.py` (Chromium and WebKit) pass.

Crops: `itten-output/compare-0.2-r2/pairs/` (ignored), 198 pairs, before left and after right,
every changed block at 320, 640, 720, 920, 1000 and 1080 px in both engines;
`linediff-sweep.json` there has the lines of every changed block, 320 to 1920 px.

