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
| Distinct spacing values at one width (table below) | 22 to 24 | 8 to 9 | 16 to 17 | 8 |
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

## Distinct spacing values per width

Margins, paddings and gaps the inspector measured at each width of the sweep, without the
system file. Before: the baseline (`5709398`). After: the build with the scale in rem.

| Width (px) | Main, before → after | Companion, before → after |
|---|---|---|
| 320 | 22 → 9 | 17 → 8 |
| 360 | 22 → 9 | 17 → 8 |
| 400 | 22 → 9 | 17 → 8 |
| 440 | 23 → 8 | 17 → 8 |
| 480 | 23 → 9 | 17 → 8 |
| 520 | 23 → 9 | 17 → 8 |
| 560 | 23 → 9 | 17 → 8 |
| 600 | 22 → 8 | 17 → 8 |
| 640 | 23 → 9 | 17 → 8 |
| 680 | 22 → 9 | 16 → 8 |
| 720 | 23 → 9 | 17 → 8 |
| 760 | 24 → 9 | 17 → 8 |
| 800 | 23 → 9 | 17 → 8 |
| 840 | 24 → 9 | 17 → 8 |
| 880 | 24 → 9 | 17 → 8 |
| 920 | 24 → 9 | 17 → 8 |
| 960 | 24 → 9 | 17 → 8 |
| 1000 | 24 → 9 | 17 → 8 |
| 1040 | 24 → 9 | 17 → 8 |
| 1080 | 24 → 9 | 17 → 8 |
| 1120 | 24 → 9 | 17 → 8 |
| 1160 | 24 → 9 | 17 → 8 |
| 1200 | 23 → 8 | 17 → 8 |
| 1240 | 24 → 9 | 17 → 8 |
| 1280 | 24 → 9 | 17 → 8 |
| 1320 | 24 → 9 | 17 → 8 |
| 1360 | 24 → 9 | 17 → 8 |
| 1400 | 24 → 9 | 17 → 8 |
| 1440 | 24 → 9 | 17 → 8 |
| 1480 | 24 → 9 | 17 → 8 |
| 1520 | 24 → 9 | 17 → 8 |
| 1560 | 24 → 9 | 17 → 8 |
| 1600 | 24 → 9 | 17 → 8 |
| 1640 | 24 → 9 | 17 → 8 |
| 1680 | 24 → 9 | 17 → 8 |
| 1720 | 24 → 9 | 17 → 8 |
| 1760 | 24 → 9 | 17 → 8 |
| 1800 | 24 → 9 | 17 → 8 |
| 1840 | 24 → 9 | 17 → 8 |
| 1880 | 24 → 9 | 17 → 8 |
| 1920 | 24 → 9 | 17 → 8 |

After, the main page has the eight steps plus one value of the gutter at each width; at 440,
600 and 1200 px the rounded gutter (18, 24, 48 px) is itself a step, so there are 8. The
companion page has no gutter padding and uses exactly the eight steps at every width.

## Spacing in rem

Every use of a scale step (margins, paddings, gaps, and the plate label and error insets) is
written in rem: 6 px = 0.375rem up to 96 px = 6rem, so gaps follow the browser's text size as
the type does. The gutter, the Plate I title block's placement, control sizes (44 px targets),
figure heights, borders and print stay as they were.

Nothing changes at the default text size. Full-page screenshots in Chromium, of both pages at
390 and 1440 px, before (`fc449e3`, scale in px) and after: 0 differing pixels in all four
(main 390 × 25756 and 1440 × 18866; companion 390 × 21308 and 1440 × 14144), with every
figure drawn and no failed request. A control pair (the px build loaded twice) also differs by
0 pixels. The inspector counts above and below are unchanged by the move, the 200 % pass has
no overflow, clipping or overlap, and `tools/check_layout.py` passes in Chromium and WebKit.

At 200 % text the gaps now double with the type: § 1 of the main page is 4893 px tall at
390 px (4737 with px gaps) and 1869 at 1440 (1689). The pairs are
`itten-output/compare/text-200pct-{390,1440}-pair.png` and
`text-200pct-companion-390-pair.png`: baseline, rem type with px gaps, rem type and gaps.

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
