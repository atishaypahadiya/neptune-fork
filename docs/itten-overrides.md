# ITTEN overrides

ITTEN version: v0.2.0 (d998de0). System: `itten.system.json` version 1.

Every departure from the system or from a protocol rule, with its reason. An element marked
`data-itten-ignore` has a row here.

| Date | Where | Rule or check | System value | Used instead | Reason | Floors still met |
|---|---|---|---|---|---|---|
| 2026-11-15 | Companion title (`.p1b .paper-title h1`, `data-itten-ignore="type-scale"`) | `type-scale` line height | the `title-companion` step, 1.08 × size | the same, 1.08 × size | Not a departure: at about 833 and 1000 px wide this title and h2 are the same size (32 and 36 px), and the inspector judges a size against the first step that produces it, so it holds the title to h2's 1.05 leading. Reordering the steps moves the clash onto the eight h2s. Worth proposing to ITTEN: match fluid steps by element as well as size. | Yes |
| 2026-11-15 | Figure lettering inside SVG (`.fig .lbl` 13 px, `.lbl.row` 14 px, `.tick` 12.5 px, `.letter` 19 px; sizes set in `site/js/plates.js`) | CRF-06, CRF-07 (type only from the scale, no fractions), TYP-03 | caption 14 / note 15 / body steps, in rem | px sizes in the drawing's coordinates | Lettering is part of the drawing, placed against marks and held to `tools/check_layout.py`'s text-over-marks and label-length checks; it scales with the plate, not with the text setting, as in an engraved plate. The inspector does not judge inside graphics. | Yes (`--line` on `--prussian` 8.9 : 1) |
| 2026-11-15 | Plate I title block placement (`.titleblock` bottom) | CRF-01 (lengths from the system) | a spacing step | `clamp(28px, 6vh, 64px)`; from 900 px `clamp(40px, 8vh, 80px)` | It is a position on the plate, not a gap between blocks: the title sits on the plate's lower edge in proportion to the plate's height, which follows the window height. | Yes |
| 2026-11-15 | Tooltip placement (`site/js/plate.js`: 8 px from the plate's edge, 14 px right of and 10 px above the pointer) | CRF-01 | a spacing step | 8, 10, 14 px | Offsets from the pointer, set while hovering; not layout. | Yes |
| 2026-11-15 | Focus ring (`outline: 2px`, `outline-offset: 2px`) | CRF-01 | a spacing step | 2 px | A ring drawn outside the box, not a gap between blocks; 2 px keeps it next to the control it marks. | Yes (3 : 1 non-text, `--prussian` on `--paper` 8.5 : 1) |
| 2026-11-15 | Print (`@media print`) | the whole system | px steps | pt sizes and the em spacing print had before this pass | Print is a separate medium with its own 10.5 pt body; the inspector does not run in print, and `tools/check_layout.py` checks print layout separately. Unchanged by this pass. | n/a |
