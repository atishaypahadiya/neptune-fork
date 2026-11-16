# Issues for ITTEN

ITTEN version: v0.1.0 (c1962e2) for issues 1 to 3, v0.2.0 (d998de0) for issues 4 and 5. Found while
running ITTEN's loop on this site (`docs/itten/baseline.md`, `docs/itten/after.md`,
`docs/itten/after-0.2.md`). Written to be filed upstream; nothing in `itten/` was changed.

Each reproduction is a single page and, where needed, a system file. Serve the folder with
any static server that handles parallel requests, then run the command shown. The outputs
quoted are from Chromium 159.0.8154.6 on macOS, 2026-11-15.

---

## 1. Text is matched to a step by size alone

**What happens.** `type-scale` gives each font size the first step in `type.scale` that
produces it at that width (`stepAt` in `inspector/inspect.js`), then checks the line height
against that step. Weight is not compared, and nothing else is. Two steps that share a size,
or two fluid steps that cross at some width, cannot both be used: text on the second one is
held to the first one's line height.

**Reproduction.** `index.html`:

```html
<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<style>
  body { margin: 0; padding: 24px; font: 400 1.125rem/1.75rem Georgia, serif; }
  p { margin: 0 0 24px; }
  h3 { margin: 0 0 24px; font: 600 1.125rem/1.5625rem Georgia, serif; } /* the "h3" step: 18/25, bold */
</style>
<h3>A heading at the body size, bold, with tighter leading</h3>
<p>Body text, 18/28. The system declares both steps; they share a size and differ in line height and weight.</p>
```

`itten.system.json`:

```json
{ "name": "repro-1", "source": { "teardown": "n/a" },
  "spacing": { "scale": [24] },
  "type": { "families": [{ "name": "Georgia", "role": "text" }],
    "scale": [ { "name": "body", "size": 18, "lineHeight": 28, "weight": 400 },
               { "name": "h3", "size": 18, "lineHeight": 25, "weight": 600 } ] },
  "grid": { "breakpoints": [{ "minWidth": 0, "columns": 1, "gutter": 0, "margin": 24 }] },
  "radius": { "scale": [0], "nesting": "concentric" },
  "colour": { "roles": { "text": { "light": "#000" } } },
  "motion": { "durations": {}, "easings": {}, "reducedMotion": "instant" } }
```

```sh
node itten/inspector/inspect.js http://127.0.0.1:8795/1/ --system itten.system.json --no-engines --out out/1
```

**Observed.** `0 errors, 1 warnings`: *Line height 25 px at 18 px; the "body" step sets 28 px;
on 1 text element.* The h3 is exactly on its own declared step.

**Expected.** No finding: the h3 matches the `h3` step in size, line height and weight.

**Fluid form, as met on this site.** `h2` is `{min 30, max 42, minWidth 750, maxWidth 1250}`
with leading 1.05; the companion title is `{min 32, max 50, minWidth 888.89, maxWidth 1388.89}`
with leading 1.08. At about 833 px both are 32 px, and at 1000 px both are 36 px. With `h2`
listed first, the title gets *Line height 34.56 px at 32 px; the "h2" step sets 33.6 px*.
Listing the title first moves the same warning onto the eight h2s. We logged it and marked the
title `data-itten-ignore="type-scale"`, which also turns off its size check
(`docs/itten-overrides.md`).

**Suggestion.** When several steps fit the size at a width, accept the element if any of them
also fits its line height (and weight, where the step gives one); report only when none
does, naming the candidates.

---

## 2. Measuring a page before its figures have drawn, without saying so

**What happens.** The inspector loads the page, waits for `networkidle` (giving up silently
after 10 s), for `document.fonts.ready` and for 500 ms, then measures. Two things can leave it
measuring a page that has not finished drawing, and the report gives no sign of either:

- content drawn later than that (a timer, an animation frame chain, a slow request) is never
  measured;
- a request that fails (a 404, or a connection reset) is not recorded. A data-driven page
  catches the error and shows a message, so there is no page error either, and the inspector
  measures the error state as if it were the page.

**Reproduction.** `index.html` (no `data.json` beside it):

```html
<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<style>
  body { margin: 0; padding: 24px; font: 400 1.125rem/1.75rem Georgia, serif; }
  p { margin: 0 0 24px; }
  .late, .data { font-size: 13px; line-height: 1.5; margin: 0 0 7px; } /* off the scale on purpose */
</style>
<p>Two figures below. One draws from data.json, which is missing (404); the page catches the
error and shows a message, as a data-driven page does. The other draws 2 s after load.</p>
<div id="fig-data"></div>
<div id="fig-late"></div>
<script>
  fetch('data.json').then(r => { if (!r.ok) throw new Error('HTTP ' + r.status); return r.json(); })
    .then(d => { fig('fig-data', 'data', d.label); })
    .catch(e => { document.getElementById('fig-data').innerHTML = '<p>This figure could not be drawn: ' + e.message + '</p>'; });
  setTimeout(() => fig('fig-late', 'late', 'Drawn 2 s after load, at 13 px with a 7 px margin'), 2000);
  function fig(id, cls, text) { const p = document.createElement('p'); p.className = cls; p.textContent = text; document.getElementById(id).append(p); }
</script>
```

System: the file from issue 1 with only the `body` step.

```sh
node itten/inspector/inspect.js http://127.0.0.1:8795/2/ --system itten.system.json --no-engines --out out/2
```

**Observed.** `0 errors, 0 warnings, 0 info`. `report.json` has `"pageErrors": []` and no record
of the 404. The late paragraph (13 px, 7 px margin, both off the system) is not in
`distinct.fontSizes`, which lists only 18 px. In a browser, 2.5 s after load, the page shows
both "This figure could not be drawn: HTTP 404" and the late paragraph.

**Control.** The same page with the timer at 0 ms: `0 errors, 3 warnings` (*Font size 13 px is
not on the type scale*, *7 px is not on the spacing scale*, and `text-resize`). So the checks
work; they ran before the content existed.

**As met on this site.** Served by Python's `http.server` (listen backlog 5), the main page's
~30 parallel data requests were intermittently reset (`net::ERR_CONNECTION_RESET`). In one run
every figure showed "This figure could not be drawn: Failed to fetch", and the report measured
those error messages as part of the layout (`near-miss-edge` warnings on the error text) with no
mention of the failed requests. We caught it only because the findings named the error text.
That report was overwritten by the clean rerun, so this page's reproduction stands in for it.

**Suggestion.**
- Record failed requests and HTTP responses of 400 or more in `report.json` and at the top of
  `report.html`, as `pageErrors` already is, and exit non-zero with `--fail-on error`.
- Add `--wait-for <selector>` (and perhaps `--settle <ms>`) so a page can say when it has drawn.
- Say in the report when the `networkidle` wait timed out.

---

## 3. Fluid steps in distinct-value mode

**What happens.** Without `--system`, the inspector cannot be told which fluid values are
deliberate, and it treats type and spacing differently:

- a fluid font size is always a warning (*Declare the range as a fluid step in the system*),
  which cannot be done without a system file;
- a fluid spacing value is a warning only when some of its values are fractional: the same
  clamp rounded to whole pixels passes silently;
- every value a fluid step takes across the sweep is counted as a distinct size, so one
  heading adds 15 sizes to "Values in use", and a page with a few fluid steps looks far
  less consistent than it is.

**Reproduction.** `index.html`, no system file:

```html
<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<style>
  :root { --gutter: round(clamp(16px, 4vw, 56px), 1px); }
  body { margin: 0; font: 400 1.125rem/1.75rem Georgia, serif; }
  main { padding: 0 var(--gutter); }              /* fluid on purpose, whole px at every width */
  h2 { margin: 0 0 24px; font: 400 clamp(1.875rem, 0.75rem + 2.4vw, 2.625rem)/1.05 Georgia, serif; } /* fluid on purpose */
  p { margin: 0 0 24px; }
</style>
<main><h2>A heading that grows from 30 to 42 px between 750 and 1250 px</h2>
<p>Side padding that grows from 16 to 56 px, rounded to whole pixels.</p></main>
```

```sh
node itten/inspector/inspect.js http://127.0.0.1:8795/3/ --widths 320:1920:40 --no-engines --out out/3
```

**Observed.** `0 errors, 1 warnings`: *The font size of 1 text element scales with the
viewport, from 30 to 42 px across the sweep: 15 distinct values, 12 fractional. Declare the
range as a fluid step in the system…*. The padding, also fluid from 16 to 56 px, gets no
finding; it appears only in `distinct.fluid.spacing`. `distinct.fontSizes` lists 16 sizes for
a page with two styles of text.

**As met on this site.** The baseline and the final run in this mode report the h2, the
Plate III year and the companion byline as warnings, though all three are declared in
`itten.system.json` and pass with it; the gutter, declared fluid in the overrides, is not
reported. The main page's distinct font sizes over the sweep read 48 before and 38 after,
although at any one width there are 6 and 5. We report per-width counts for that reason.

**Suggestion.**
- In distinct-value mode, report each fluid range once, as info, with its range and the widths
  it spans (the data is already in `distinct.fluid`), and treat type and spacing alike.
- Count distinct values per width as the headline figure, with fluid ranges counted as one
  value each in the sweep totals.
- Optionally accept `--fluid-ok` or a minimal system file that declares only fluid ranges,
  for teardowns of references that have no system yet.

---

## 4. Text edges measured from a block's own text nodes only

Found with ITTEN v0.2.0 (d998de0). **Fixed in ITTEN 0.3** (`near-miss-edge` measures a text
block's edges over all its text, inline children included). With ITTEN v0.3.2 the reference
lists and the citation blocks clear without `data-itten-ignore`; the marker and the override
were removed on 2026-11-15.

**What happens.** `near-miss-edge` takes a text block's left and right edges from the client
rects of the block's own text nodes (`e.textNodes` in `inspector/collect.js`: direct children of
type text). Text inside inline children (`span`, `i`, `a`) is left out. When every line of a block
begins with an inline child, the measured edge falls somewhere inside the first line, and two
such blocks whose text starts on the same x can read as 0.5 to 3 px apart.

**Reproduction.** `4/index.html`:

```html
<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<style>
  body { margin: 0; padding: 24px; font: 400 1rem/1.5rem Georgia, serif; }
  ol { margin: 0; padding-left: 24px; }
  li { margin: 0 0 12px; }
</style>
<ol>
  <li><span>Brookes</span>, C. J. 1970.</li>
  <li><span>Bouvard</span>, A. 1821.</li>
</ol>
```

```sh
node itten/inspector/inspect.js http://127.0.0.1:8796/4/ --widths 400 --no-engines --out out/4
```

**Observed.** `near-miss-edge`, warning: *Left edges 2.56 px apart: li ", C. J. 1970." at
105.48, li ", A. 1821." at 108.05.* Both items' text starts at 48 px.

**As met on this site.** To keep WebKit from hyphenating names, `appendText()`
(`site/js/plate.js`) writes capitalised words in a span. Pairs of references and quote sources
then read as 1.4 to 2.6 px apart, in Chromium and WebKit, at widths where nothing moved
(measured with all inline content: 45 and 45 px, 88 and 88 px). Logged in
`docs/itten-overrides.md` with `data-itten-ignore="near-miss-edge"`.

**Suggestion.** Measure a text block's line boxes over all its text, inline descendants
included (a range over the block's contents, as `lines()` does per node), and keep the
own-text-node rule only for deciding which blocks are text.

---

## 5. An arXiv id or a path as a widow

Found with ITTEN v0.2.0 (d998de0). **Fixed in ITTEN 0.3** (`widow` does not count a last line
that is an arXiv id, another prefixed id or a path). With ITTEN v0.3.2 the companion's
reference list and data-availability paragraph clear without `data-itten-ignore`; the
markers and the override were removed on 2026-11-15.

**What happens.** `widow` does not count a last line that is a URL or a DOI
(`inspector/typography.js`: `^(https?://|www\.|doi:)` or a host followed by a path). Other
identifiers of the same kind, an arXiv id (`arXiv:2108.06305`) or a repository path
(`docs/notes/P1B_challis_detection.json`), are reported. None of them can be reworded or
hyphenated, and joining one to the word before it with a no-break space makes a longer
unbreakable run, which is worse at a narrow measure.

**Reproduction.** `5/index.html`:

```html
<!doctype html><html lang="en"><meta charset="utf-8"><link rel="icon" href="data:,"><meta name="viewport" content="width=device-width">
<style>
  body { margin: 0; padding: 1.5rem; overflow-wrap: break-word; font: 400 1rem/1.5rem Georgia, serif; }
  p { margin: 0 0 1.5rem; max-width: 20rem; text-wrap: pretty; }
</style>
<p>A chapter in an edited book on the discovery of the planet, with both of its identifiers listed at the end of the line: arXiv:2108.06305.</p>
<p>The figures are built from notes kept in the repository, each named for its step: docs/notes/P1B_adams_fit.md, docs/notes/P1B_challis_detection.json.</p>
<p>The same paragraph ending on an address the inspector already exempts from widows: https://doi.org/10.1007/978-3-030-54218-4_6.</p>
```

```sh
node itten/inspector/inspect.js http://127.0.0.1:8796/5/ --widths 400 --no-engines --out out/5
```

**Observed.** Two `widow` warnings: *The last line holds one word, "arXiv:2108.06305.", under 3
full lines, with text-wrap: pretty set.* and the same for *"docs/notes/P1B_challis_…"*. The
third paragraph, ending on a DOI URL, is not reported. (One `rag` warning on the same page is
not part of this issue.)

**As met on this site.** The companion's reference to Krajnović (2021) ends on its arXiv id,
alone on the last line at 360 px, and its first data-availability paragraph ends on the path
of a note, alone on the last line from 360 to 640 px. Logged in `docs/itten-overrides.md`
with `data-itten-ignore="widow"`.

**Suggestion.** Treat a last line that is one identifier as not a widow: besides URLs and
DOIs, a `scheme:` prefix with no spaces (`arXiv:`, `urn:`, `isbn:`) and a token holding a `/`
or `_` with no spaces.
