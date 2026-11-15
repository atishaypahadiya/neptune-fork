# The system, with the old values mapped to it

ITTEN version: v0.1.0 (c1962e2). System: `itten.system.json` version 1, derived in
`docs/teardown.md`. This is the plan for step 4; nothing here is applied yet.

Δ is new minus old, in CSS px, at the width named or at every width.

## Spacing

Scale: **6, 12, 18, 24, 36, 48, 72, 96** (6 × 1, 2, 3, 4, 6, 8, 12, 16).

| Old | New | Δ | Where |
|---|---|---|---|
| 2 | 0 | −2 | Table 1, second footer row, padding-top |
| 4 | 6 | +2 | Citation block, row gap |
| 5 | 6 | +1 | Table 1 cells, padding top and bottom (main page) |
| 6 | 6 | 0 | Section mark below; quote text and translation below; phone cell padding-right |
| `0.4em` (6.6 to 7.2) | 6 | −0.6 to −1.2 | After "§ n" in an h3 |
| 7 | 6 | −1 | Companion tables, cell padding top and bottom |
| 8 | 6 | −2 | Table 1 footer padding-top; Plate III controls above; stacked plate text below |
| `0.5em` (7.5 to 9.5) | 6 | −1.5 to −3.5 | Plate I subtitle below (to the byline) |
| `0.6em` (9) | 12 | +3 | Between references |
| `0.22em` (8.8 to 17.7) | 12 below 900 px, 18 from 900 px | +3.2 at 390; +6.9 at 900; +0.3 at 1440 | Plate I title below |
| `0.3em` (9.6 to 15) | 12 | +2.4 to −3 | Companion title below |
| 10 | 12 | +2 | Figure caption above; table cell padding-right |
| `0.7em` (11.55 to 12.6) | 12 | −0.6 to +0.45 | Between Limitations items |
| 12 | 12 | 0 | h3 below; citation key; table caption; Plate III year below; control label gap |
| 14 | 12 | −2 | Copy button side padding; Plate I caption above |
| 15 (`1em` at 15 px) | 18 | +3 | Below the Plate III note |
| 16 | 18 | +2 | Plate III controls gap; citation column gap; Fig. 3 below |
| `1em` (16.5 to 18) | 18 | 0 from 762 px; +1.5 at 390 | Between paragraphs |
| `1.1em` (18.15 to 19.8) | 18 | −0.15 to −1.8 | Limitations list indent (bullets) |
| 21 (`1.4em` at 15 px) | 24 | +3 | Companion "A companion study to" line below |
| 22 | 24 | +2 | h2 below |
| 24 | 24 | 0 | Figs. 8 and 9 stacked gap; citation block; Fig. 3 above; references indent |
| `1.6em` (26.4 to 28.8) | 24 | −2.4 to −4.8 | Quote indent |
| 28 | 24 | −4 | Below a quote; below Plate III controls |
| 28 (phone), 32 (721 to 899 px) | 36 | +8, +4 | Stacked plate, space to its text (now the figures' 36) |
| 32 | 36 | +4 | Main Table 1 above and below (now the figures' 36, as on the companion); Figs. 8 and 9 side-by-side gap |
| 36 | 36 | 0 | Figures above and below |
| 48 | 48 | 0 | h3 above |
| 64 | 72 | +8 | Section mark above, on phones; plates above and below, on phones; Plates II and III text column, top and bottom padding |
| 72 | 72 | 0 | Section mark above, companion |
| 81 (`9vh` in a 900 px window; 48 to 112) | 72 below 900 px, 96 from 900 px | −9 to +15 at 900 px tall | Companion title above |
| 88 | 96 | +8 | Section mark above, main page |
| 96 | 96 | 0 | Plates above and below; end of page |

Positions on plates, held to the same scale:

| Old | New | Δ | Where |
|---|---|---|---|
| 22 | 24 | +2 | Plate I label from the top |
| 20 | 24 | +4 | Figure error message inset (only shown when data fail to load) |
| 6 × 9 | 6 × 12 | +3 wide | Tooltip padding (hover state) |
| 18, 24 | 18, 24 | 0 | Plate label from the bottom and right; plate error inset |

Not mapped (logged in `docs/itten-overrides.md`):

- **The gutter**, `clamp(16px, 4vw, 56px)`: the page's side padding and outer grid track,
  fluid by design (`docs/DESIGN.md`, "Widths"). It becomes `round(…, 1px)`, so it is whole
  pixels at every width (a change of at most 0.5 px); the range stays.
- **Plate I's title block placement**, `clamp(28px, 6vh, 64px)` from the bottom (40 to 80 px
  from 900 px): it is placed by the window height so the title sits on the plate's lower edge.

## Type

Every size in rem, so all text follows the browser's text size (TYP-03); at the default 16 px
the sizes below are what renders.

| Old | New step | Δ | Where |
|---|---|---|---|
| 14 px, line height 21 or 21.84 | caption 14/21 | line height −0.84 in table headers | Captions, table headers, copy button, tooltip, quote source, phone observatory line |
| 15 px / 23.4 | note 15/23 | line height −0.4 | Section marks, notes, translation, controls label, citation, references, companion "back" line, error |
| body `clamp(16.5px, 1.05vw + 10px, 18px)` / 1.56 | body-phone 17/26 to 720 px; body 18/28 from 721 px | +0.5 at 390 (16.5 → 17); −0.56 at 720 (17.56 → 17); +0.23 at 740 (17.77 → 18); line height 25.74 → 26, 28.08 → 28 | All prose |
| quote `1.06em` (17.49 to 19.08) | body | −0.5 at 390; −1.08 at 1440 | Archival quotes (3) |
| h3 `1em` / 1.4 | body (17/26, 18/28), weight 600 | line height +2.9 at 17 px, +2.8 at 18 px | h3 |
| Plate label 18 / 28.08 | body 18/28 | line height −0.08 | "Plate I." etc. |
| h2, Plate III year, Plate I title (both ramps), companion title, subtitle, byline, Table 1 cells (both pages) | fluid steps, the same ramps | 0 at the default text size | — |

Figure lettering inside SVG (12.5, 13, 14, 19 px) is drawing geometry and is not on the text
scale; it is logged, not changed.

## Would a reader notice

The largest changes are 8 px on gaps of 64 and 88 px (9 to 13 %), the quote indent moving by
about 5 px, quotes set at the body size instead of 1 px
larger, and the phone body going from 16.5 to 17 px. None changes a colour, a typeface, a
width, the order of anything, or a gap by more than 8 px, and no text size moves by more than
1.1 px at the default text size. My judgement is that a reader would not notice any of them,
so step 4 goes ahead; the before/after crops in the final report are the check on that.

The one change readers will see is the one asked for: with a larger default text size in the
browser, the text now grows.
