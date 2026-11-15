# ITTEN after the fixes

ITTEN version: v0.1.0 (c1962e2). System: `itten.system.json` version 1. Date: 2026-11-15.

Same run as `docs/itten/baseline.md` (`--full`, 320 to 1920 px every 40 px, 200 % text, light
and dark, DPR 2), once with `--system itten.system.json` and once without it (distinct-value
mode, to compare with the baseline). Engines: Chromium 159.0.8154.6 and WebKit 27.2, with no
finding that differs between them. Firefox still does not start on this machine ("Could not
find profile folder"). Served by a small Node static server: Python's `http.server` (listen
backlog 5) reset some of the page's ~30 parallel data requests, so one earlier run measured
figures that had not drawn; these numbers come from runs where every figure drew.
`tools/check_prerender.py` and `tools/check_layout.py` (Chromium and
WebKit) pass.

## Counts

| | Main before | Main after | Companion before | Companion after |
|---|---|---|---|---|
| Errors (against the system) | — | 0 | — | 0 |
| Warnings (against the system) | — | 1 | — | 0 |
| Info | 2 | 2 | 0 | 0 |
| Warnings, distinct-value mode | 12 | 2 | 10 | 2 |
| Distinct spacing values at one width | 22 to 24 | 8 to 9 | 16 to 17 | 8 |
| Distinct spacing values over the sweep | 68 | 31 | 49 | 8 |
| Fixed spacing values (px) | 18 values, 4 to 96 | 6, 12, 18, 24, 36, 48, 72, 96 | 12 values, 6 to 96 | 6, 12, 18, 24, 36, 48, 72, 96 |
| Distinct font sizes at one width | 6 | 5 | 4 to 6 | 4 to 5 |
| Text that ignores the browser text size | 234 of 234 | 0 | 154 of 154 | 0 |
| Elements marked `data-itten-ignore` | 0 | 0 | 0 | 1 |

On the main page the 31 spacing values over the sweep are the eight steps plus the fluid side
padding (`--gutter`, whole pixels from 16 to 56). The font sizes over the sweep (38 main, 25
companion) are the declared fluid steps; at any one width there are at most five.

The distinct-value warnings that remain are fluid steps (h2 and the Plate III year on the main
page, h2 and the byline on the companion), which that mode cannot tell from undeclared ones;
with the system file they are declared and pass.

## What remains

| Page | Finding | Severity | Status |
|---|---|---|---|
| Main | Plates II and III text columns: side padding follows the window (16 to 56 px) | warning | Logged in `docs/itten-overrides.md` (the gutter) |
| Main | Play button border at y = 9120.48 px at 1×; could not be profiled | info | Fractional position from the fluid h2 and Table 1 heights above it; WebKit and Chromium both draw it, the line is a circle |
| Main | Copy buttons' borders at fractional positions | info | Both engines snap them crisp (86 of 86 samples, 8 of 8) |

The companion title's step clash is marked `data-itten-ignore="type-scale"` and logged.

## Before and after crops

In `itten-output/compare/` (ignored), each `<region>-<width>-pair.png` shows before on the left,
after on the right; the single crops are `-before.png` and `-after.png`.
