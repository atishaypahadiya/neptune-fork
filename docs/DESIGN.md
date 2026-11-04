# neptune-fork site design (M7)

Mockups: https://claude.ai/artifact/HM1xYNzooiZxzS9M22oEVC (artboards Scene1–Scene4 are the approved direction; the others are rejected explorations, keep them out).

## The language in one line
Scenes, not interfaces: black space, a glowing horizon, perspective geometry, soft light against hairline lines, annotations like a pilot's manual, real equations as type. Greyscale only.

## Rules
- **One idea per scene, drawn as geometry.** The drift is the gap between two paths. The pull is a bent grid. Time is stacked planes. The answer is converging lines. Never a plain chart if geometry can carry it.
- **Every number is real.** Values come from `site/data/*.json` produced by the pipeline. Anything not yet computed shows a visible placeholder like `[M5B]`, never an invented value.
- **Every equation is real and belongs to its scene.** No decorative maths.
- **Greyscale only.** No accent color. Emphasis comes from brightness, weight and line style (solid = observed/true, dashed = predicted/assumed).
- **UI rule (project):** any pattern is allowed if it adds real value, visual or functional. Each scene must still pass a self-audit listing what it uses and why.

## Tokens
| Token | Value |
|---|---|
| Background | `#000000` (Scene 4 flips to a vertical black→`#EFEFEC` horizon gradient) |
| Ink on dark | `#FFFFFF` lines, `#EDEDED` text, `#B5B5B5` secondary, `#8F8F8F` tertiary |
| Ink on light | `#151515` lines/text, `#3A3A38` secondary |
| Labels | Barlow Condensed 500, caps, letter-spacing 0.12–0.22em, 13–15px |
| Headlines | Barlow Condensed 300, 60–64px, line-height 1.02, sentence case, no trailing period |
| Body | Barlow Condensed 400, 18px, line-height 1.45, max ~31rem |
| Equations | EB Garamond italic, 18–22px |
| Hairlines | 0.6–0.8px; main paths 2px solid; predicted 1.2px dash 6/6 |
| Grain | fractal noise (baseFrequency 0.9) over everything, opacity 0.10–0.14 |
| Glow | blurred radial ellipse on the horizon (stdDeviation ~18–24) |

## Components (shared by every scene)
- **Stage annotation:** 8px solid square on the point → 0.7px leader line → caps title (15px) + one or two secondary lines (13px). Pilot-manual style.
- **Readout (top right):** three stacked caps lines, tabular numbers; the most important one in full white.
- **Vertical edge title:** right edge, `writing-mode: vertical-rl`, letter-spacing 0.4em: `NEPTUNE-FORK — REDISCOVERING 1846`.
- **Chapter label (top left):** `0N · NAME` above the headline and a two-to-three sentence body.
- **Registration marks:** small `+` crosses scattered in empty sky, opacity ~0.55.
- **Legend + next (bottom right):** line samples (solid/dashed) and `NEXT · … →`.
- **Shaded sphere:** radial gradient white → mid grey → black with a dashed halo ring; light source matches the scene.

## Scenes
1. **Something unseen.** Perspective grid of the ecliptic plane to a glowing horizon. Observed vs predicted Uranus paths diverging (gap = residual, exaggerated and labelled as such). Stages 1690 / 1781 / 1822 / 1846. Unseen shaded sphere with a cone of lines onto the 1822 point. Equations: Kepler III, Newton, Le Verrier's disturbing function. Data: `residuals.json`.
2. **The pull.** Same grid, sinking into a 1/r gravity-well dish under the sphere. Observed path bends toward the well, predicted runs straight. Readout: 1822 separation ≈ 10.9 AU, pull ≈ 2.6 × 10⁻⁹ m/s². Equation: the planetary N-body equation of motion. Data: integrator states.
3. **Time as layers.** Isometric wireframe box on a dark base slab; 16 translucent planes 1690→1846; residual point on each plane joined into a rising curve; dashed zero wall. Equations: Δλ and χ². Data: `residuals.json`.
4. **The answer.** Inverted palette: black sky over a bright horizon. Sphere lit from below, cone of lines converging on Neptune's true 1846 position (diamond). Sky strip: free-distance posterior as a soft glow + crisp line, fixed-38.4 AU posterior dashed and narrow. Equation: Bayes. Data: `skymap_1846.json`, `m5a_profile.json`, M5b summary.

## Build
- Static site in `site/`, deployed by the existing Pages workflow.
- Scenes 1, 2 and 4 in **three.js** (real 3D: orbit the camera slightly on scroll/drag, animate the observed path drawing in). Scene 3 may stay SVG isometric if three.js adds nothing.
- Interaction, in priority order: (1) a time scrubber that moves Uranus along the path and updates the readout; (2) toggle free vs fixed distance in Scene 4; (3) hover a stack plane in Scene 3 to show that year's numbers.
- One scene per viewport; scroll or `NEXT` advances. No long scrolling page.
- Respect `prefers-reduced-motion`. Phone layout: scenes stack, labels shrink, the vertical edge title is hidden.

## Acceptance (adds to A5)
- All four scenes render from `site/data/` with no hard-coded science numbers.
- Self-audit per scene committed to `docs/notes/M7.md`.
- Lighthouse accessibility ≥ 90; every scene has a text alternative describing what it shows.
