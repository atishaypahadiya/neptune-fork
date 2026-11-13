"""Layout, print and interaction acceptance checks for the site (docs/DESIGN.md, Acceptance).

Runs on every page of the site (the main paper and its companion, site/adams-challis/; the
checks of Plates II and III apply where the page has them). Loads the page with JavaScript at
360, 390, 430, 768, 899, 900, 1024 and 1440 px (phones as touch devices) and in print (A4 from an 800 px window), and fails on any of:

- **Text over marks:** in every plate and figure, a text box (SVG text, the lettering along
  orbits glyph by glyph, and the HTML text laid over a plate) that intersects a stroke, a
  marker or another text box. The one exception is a text's own leader line
  (``data-leader`` on the text, ``data-leader-of`` on the line).
- **Labels along a path** (``textPath``) longer than their path, so that their ends would be
  cut off. The length is measured with letter-spacing, which some engines leave out.
- Text leaving its figure; horizontal overflow of the page.
- **Grid lines** (docs/DESIGN.md, "Widths"): the six named column lines are read from the page
  at each width (and at 660 and 720 px, either side of where wide meets body). They must be
  the ones the page's own ``--measure``, ``--overhang`` and ``--gutter`` give, with wide
  centred on body and full on the window's edges. Every figure, plate and table has exactly
  one of the classes ``body``, ``wide``, ``full``, and its left and right edges are on the two
  lines of its class within 1 px; so are every block of the header and of the sections (body
  without a class) and every caption. None has an inline width or a transform. The exceptions
  are ``GRID_EXEMPT``, each with a rule of its own. In print the three widths are one.
- **Emoji:** any emoji or emoji-presentation character (pictographs, dingbats such as ❚,
  geometric shapes such as ▶, media symbols such as ⏸, variation selector 16) in the rendered
  text. Icons are inline SVG.
- **Sticky or fixed** elements below 900 px.
- Overlapping blocks in the text column of Plates II and III, and inside Plate III's controls.
- Table 1: a header or a year range on more than one line.
- Side by side, a text column under 45 characters a line; a vertical gap over 200 px.
- The number audit (tools/check_prerender.py's) at every width.
- Every point letter in the prose lights its point; Plate III's play, pause and scrub work
  and the button's icon is SVG.
- Print: nothing sticky, controls hidden, every caption under its figure.

    uv run --group site python tools/check_layout.py _site [--engine webkit] [--shots DIR]
"""

# The checks are JavaScript run in the page, kept whole here; its lines run past 100.
# ruff: noqa: E501

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_prerender import AUDIT
from prerender_site import PAGES, STRUCTURE, launch, open_rendered, serve

WIDTHS = [
    (360, 780),
    (390, 844),
    (430, 932),
    (768, 1024),
    (899, 900),
    (900, 900),
    (1024, 768),
    (1440, 900),
]
PHONE = 600  # below this width the page is loaded as a touch device
# Widths at which only the grid lines and the overflow are checked: either side of the width at
# which wide meets body (the measure plus twice the side padding, about 696 px).
GRID_WIDTHS = [(660, 900), (720, 900)]

# Figures that are not on two lines of the grid, each with the rule it is held to instead.
#   half:   one half of a full section: an outer line to the centre line (full when stacked).
#   column: fills the text column it stands in.
#   pair:   one half of the wide pair, less half the gap (wide when one is above the other).
GRID_EXEMPT = [
    {"sel": "#plate-ii > .plate", "rule": "half",
     "why": "Plate II is the right half of its full section, beside its text"},
    {"sel": "#plate-iii > .plate", "rule": "half",
     "why": "Plate III is the left half of its full section, beside its text"},
    {"sel": "#plate-iii .side > [data-fig]", "rule": "column",
     "why": "Fig. 3 stands in Plate III's text column and takes that column's width"},
    {"sel": ".pair > [data-fig]", "rule": "pair",
     "why": "Figs. 8 and 9 stand side by side from 900 px, each on half the wide pair"},
]  # fmt: skip

# The rendered text must hold no emoji or emoji-presentation character.
EMOJI = (
    r"[\p{Extended_Pictographic}\p{Emoji_Presentation}\u{FE0F}\u{20E3}"
    r"\u{2300}-\u{23FF}\u{25A0}-\u{25FF}\u{2600}-\u{27BF}\u{2B00}-\u{2BFF}]"
)

# Every text box in a figure against every stroke, marker and other text box. Runs in the page.
MARKS = r"""() => {
  const out = { collisions: [], paths: [], leaves: [] };
  const P = (x, y, m) => { const p = new DOMPoint(x, y).matrixTransform(m); return [p.x, p.y]; };
  const quad = (r, m) => [P(r.x, r.y, m), P(r.x + r.width, r.y, m),
                          P(r.x + r.width, r.y + r.height, m), P(r.x, r.y + r.height, m)];
  const aabb = (q) => {
    const xs = q.map((p) => p[0]), ys = q.map((p) => p[1]);
    return [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)];
  };
  const meets = (a, b, pad = 0) => a[0] < b[2] + pad && b[0] < a[2] + pad &&
                                   a[1] < b[3] + pad && b[1] < a[3] + pad;
  // Separating axis test for two convex quads, shrunk by `tol` px (touching is not crossing).
  const sat = (A, B, tol) => {
    for (const Q of [A, B]) for (let i = 0; i < 4; i++) {
      const [x1, y1] = Q[i], [x2, y2] = Q[(i + 1) % 4];
      const n = [y1 - y2, x2 - x1], len = Math.hypot(...n);
      if (!len) continue;
      const pa = A.map((p) => (p[0] * n[0] + p[1] * n[1]) / len);
      const pb = B.map((p) => (p[0] * n[0] + p[1] * n[1]) / len);
      if (Math.max(...pa) - tol <= Math.min(...pb) || Math.max(...pb) - tol <= Math.min(...pa)) return false;
    }
    return true;
  };
  // Distance from a point to a convex quad (0 inside).
  const dist = (Q, [x, y]) => {
    let inside = true, best = Infinity;
    for (let i = 0; i < 4; i++) {
      const [x1, y1] = Q[i], [x2, y2] = Q[(i + 1) % 4];
      const cross = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1);
      if (cross < 0) inside = false;
      const dx = x2 - x1, dy = y2 - y1, t = Math.max(0, Math.min(1, ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy || 1)));
      best = Math.min(best, Math.hypot(x - (x1 + t * dx), y - (y1 + t * dy)));
    }
    // The quads run clockwise on screen (y down); a point inside has every cross >= 0.
    return inside ? 0 : best;
  };
  const visible = (e) => e.checkVisibility ? e.checkVisibility({ opacityProperty: true, visibilityProperty: true })
    : getComputedStyle(e).display !== "none";
  const ctx = document.createElement("canvas").getContext("2d");
  const measure = (t) => {
    const cs = getComputedStyle(t);
    ctx.font = `${cs.fontStyle} ${cs.fontWeight} ${cs.fontSize} ${cs.fontFamily}`;
    const s = t.textContent, ls = parseFloat(cs.letterSpacing) || 0;
    return ctx.measureText(s).width + ls * s.length;
  };
  const name = (t) => `"${t.textContent.trim().slice(0, 40)}"`;

  for (const fig of document.querySelectorAll("[data-fig]")) {
    const svg = fig.querySelector(".fig");
    if (!svg || !visible(fig)) continue;
    const fb = (fig.querySelector(":scope > .canvas") ?? fig).getBoundingClientRect();
    const id = fig.dataset.fig;
    // Text items: each SVG text (glyph by glyph along a path) and each line of HTML text over
    // the plate. Each item: a list of screen quads.
    const texts = [];
    for (const t of svg.querySelectorAll("text")) {
      if (!visible(t) || !t.textContent.trim()) continue;
      const m = t.getScreenCTM();
      const tp = t.querySelector("textPath");
      let quads;
      if (tp) {
        quads = [];
        for (let i = 0; i < t.getNumberOfChars(); i++) {
          const r = t.getExtentOfChar(i);
          if (r.width > 0 && r.height > 0) quads.push(quad(r, m));
        }
        const path = svg.querySelector(tp.getAttribute("href"));
        const L = path.getTotalLength(), need = Math.max(measure(t), t.getComputedTextLength());
        if (need > L - 2) out.paths.push(`${id}: ${name(t)} needs ${need.toFixed(0)} px, its path is ${L.toFixed(0)} px`);
      } else {
        quads = [quad(t.getBBox(), m)];
      }
      texts.push({ el: t, quads, leader: t.dataset.leader ?? null, label: name(t) });
    }
    for (const e of fig.querySelectorAll(".plate-label, .titleblock h1, .titleblock p")) {
      if (!visible(e)) continue;
      const r = document.createRange();
      r.selectNodeContents(e);
      const quads = [...r.getClientRects()].filter((b) => b.width > 0)
        .map((b) => [[b.left, b.top], [b.right, b.top], [b.right, b.bottom], [b.left, b.bottom]]);
      texts.push({ el: e, quads, leader: null, label: name(e) });
    }
    for (const t of texts) {
      for (const q of t.quads) {
        const b = aabb(q);
        if (b[0] < fb.left - 1 || b[2] > fb.right + 1 || b[1] < fb.top - 1 || b[3] > fb.bottom + 1) {
          out.leaves.push(`${id}: ${t.label} leaves the figure`);
          break;
        }
      }
      t.box = t.quads.map(aabb).reduce((a, b) => [Math.min(a[0], b[0]), Math.min(a[1], b[1]), Math.max(a[2], b[2]), Math.max(a[3], b[3])]);
    }
    // Marks: strokes (sampled along their geometry) and markers (small filled shapes, by box).
    const marks = [];
    for (const e of svg.querySelectorAll("path, line, polyline, polygon, circle, ellipse, rect")) {
      if (e.closest("defs, text") || e.matches(".hit, .glow") || !visible(e)) continue;
      const cs = getComputedStyle(e);
      const stroked = cs.stroke !== "none" && parseFloat(cs.strokeWidth) > 0;
      const filled = cs.fill !== "none" && !/rgba\(.*,\s*0\)$/.test(cs.fill) && cs.fill !== "transparent";
      if (!stroked && !filled) continue;
      const m = e.getScreenCTM(), bb = e.getBBox();
      const scale = Math.hypot(m.a, m.b);
      const sw = stroked ? (parseFloat(cs.strokeWidth) * scale) / 2 : 0;
      const q = quad(bb, m), box = aabb(q).map((v, i) => v + (i < 2 ? -sw : sw));
      const small = bb.width * scale < 16 && bb.height * scale < 16;
      marks.push({ el: e, box, sw, q, small, filled, leaderOf: e.dataset.leaderOf ?? null, m });
    }
    const pts = (k) => {
      if (k.samples) return k.samples;
      const e = k.el, s = [];
      const len = e.getTotalLength ? e.getTotalLength() : 0;
      const n = Math.max(2, Math.ceil(len / 1));
      for (let i = 0; i <= n; i++) {
        const p = e.getPointAtLength((len * i) / n);
        s.push(P(p.x, p.y, k.m));
      }
      return (k.samples = s);
    };
    const tag = (k) => {
      const e = k.el;
      return `${e.tagName}.${e.getAttribute("class") ?? ""}${e.closest("[data-pt]") ? `[${e.closest("[data-pt]").dataset.pt}]` : ""}`;
    };
    for (const t of texts) {
      for (const k of marks) {
        if (!meets(t.box, k.box, 0.5)) continue;
        if (t.leader && t.leader === k.leaderOf) continue;
        let hit = false;
        if (k.small || (k.filled && k.el.tagName !== "path")) {
          hit = t.quads.some((q) => sat(q, k.q.map(([x, y]) => [x, y]), 0.5) &&
            meets(aabb(q), k.box, -0.5));
        } else {
          const s = pts(k);
          hit = t.quads.some((q) => { const b = aabb(q); return s.some((p) =>
            p[0] > b[0] - k.sw - 1 && p[0] < b[2] + k.sw + 1 && p[1] > b[1] - k.sw - 1 && p[1] < b[3] + k.sw + 1 &&
            dist(q, p) < Math.max(k.sw - 0.25, 0.25)); });
        }
        if (hit) out.collisions.push(`${id}: ${t.label} crosses ${tag(k)}`);
      }
    }
    for (let i = 0; i < texts.length; i++) for (let j = i + 1; j < texts.length; j++) {
      const a = texts[i], b = texts[j];
      if (!meets(a.box, b.box)) continue;
      if (a.quads.some((qa) => b.quads.some((qb) => sat(qa, qb, 0.5)))) {
        out.collisions.push(`${id}: ${a.label} overlaps ${b.label}`);
      }
    }
  }
  out.collisions = [...new Set(out.collisions)];
  return out;
}"""

PAGE = r"""({ emoji, narrow }) => {
  const vw = document.documentElement.clientWidth, out = {};
  out.overflowX = document.documentElement.scrollWidth > vw ? [`page ${document.documentElement.scrollWidth} px wide in a ${vw} px window`] : [];

  // Emoji in the rendered text, SVG labels, button text and the title.
  const re = new RegExp(emoji, "gu");
  const texts = [document.title, document.body.innerText,
    ...[...document.querySelectorAll("svg text, button")].map((e) => e.textContent)];
  out.emoji = [...new Set(texts.join("\n").match(re) ?? [])]
    .map((c) => `U+${c.codePointAt(0).toString(16).toUpperCase().padStart(4, "0")} ${c}`);

  // Nothing sticky or fixed below 900 px.
  out.sticky = [];
  if (narrow) for (const e of document.querySelectorAll("body *")) {
    const p = getComputedStyle(e).position;
    if ((p === "sticky" || p === "fixed") && e.getClientRects().length) out.sticky.push(`${e.tagName}.${e.className} is ${p}`);
  }

  // Blocks in the text column of a split plate, and Plate III's controls, never overlap.
  out.blocks = [];
  const over = (a, b) => Math.min(a.right, b.right) - Math.max(a.left, b.left) > 1 &&
                         Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top) > 1;
  const pairs = (els, where) => {
    const bs = els.filter((e) => e.getClientRects().length).map((e) => [e, e.getBoundingClientRect()]);
    for (let i = 0; i < bs.length; i++) for (let j = i + 1; j < bs.length; j++) {
      if (over(bs[i][1], bs[j][1])) out.blocks.push(`${where}: ${bs[i][0].tagName}.${bs[i][0].className} overlaps ${bs[j][0].tagName}.${bs[j][0].className}`);
    }
  };
  for (const s of document.querySelectorAll(".split")) {
    pairs([...s.children], s.id);
    pairs([...s.querySelector(".side").children], `${s.id} text column`);
  }
  const c = document.querySelector("#plate-iii .controls");
  if (c) pairs([c.querySelector(".play"), c.querySelector("label")], "Plate III controls");
  if (c) pairs([c.querySelector(".play"), c.querySelector("label > span") ?? c.querySelector("label"), c.querySelector(".scrub")].filter(Boolean), "Plate III controls");

  // Table 1: headers and year ranges on one line.
  out.table = [];
  const lines = (e) => {
    const r = document.createRange();
    r.selectNodeContents(e);
    return new Set([...r.getClientRects()].filter((b) => b.width > 0).map((b) => Math.round(b.top))).size;
  };
  for (const e of document.querySelectorAll("#table-1 th, #table-1 td.yrs")) {
    if (e.getClientRects().length && lines(e) > 1) out.table.push(`"${e.textContent}" wraps (${lines(e)} lines)`);
  }
  const tw = document.querySelector("#table-1 table").getBoundingClientRect().width;
  const fw = document.querySelector("#table-1").clientWidth;
  if (tw > fw + 1) out.table.push(`table ${Math.round(tw)} px wide in a ${fw} px column`);

  // Side panels inside the viewport; characters per line side by side; vertical gaps.
  out.columns = [];
  const ctx = document.createElement("canvas").getContext("2d");
  const cpl = (p) => {
    const cs = getComputedStyle(p);
    ctx.font = `${cs.fontStyle} ${cs.fontWeight} ${cs.fontSize} ${cs.fontFamily}`;
    const t = p.textContent.replace(/\s+/g, " ").trim();
    return Math.round(p.clientWidth / (ctx.measureText(t).width / t.length));
  };
  for (const s of document.querySelectorAll(".split")) {
    const side = s.querySelector(".side");
    for (const e of side.querySelectorAll("*")) {
      const b = e.getBoundingClientRect();
      if (b.width && (b.right > vw + 1 || b.left < -1)) out.columns.push(`${s.id}: ${e.tagName}.${e.className} outside the window`);
    }
    if (getComputedStyle(s).gridTemplateColumns.split(" ").length > 1) {
      const ps = [...side.querySelectorAll("p")].filter((p) => p.textContent.length > 120 && !p.classList.contains("caption"));
      const least = Math.min(...ps.map(cpl));
      if (least < 45) out.columns.push(`${s.id}: ${least} characters a line`);
    }
  }
  out.gaps = [];
  const body = [...document.querySelectorAll("main p, main li, main h2, main h3, main figure, main table, main blockquote, .plate, .canvas, .controls, header .caption")]
    .filter((e) => e.getBoundingClientRect().height > 0)
    .map((e) => { const b = e.getBoundingClientRect(); return [b.top + scrollY, b.bottom + scrollY, e]; })
    .sort((a, b) => a[0] - b[0]);
  let reach = body[0][1];
  for (const [t, bt, e] of body) {
    if (t - reach > 200) out.gaps.push(`${Math.round(t - reach)} px before ${e.tagName}.${e.className}`);
    reach = Math.max(reach, bt);
  }

  // Letter links: each prose letter finds its point and lights it on focus.
  out.links = [];
  let n = 0;
  for (const a of document.querySelectorAll(".pt")) {
    n++;
    const [fig, l] = a.dataset.link.split(":");
    const p = document.querySelector(`[data-fig="${fig}"] [data-pt="${l}"]`);
    if (!p) { out.links.push(`${a.dataset.link}: no point`); continue; }
    a.focus({ preventScroll: true });
    if (!p.classList.contains("on")) out.links.push(`${a.dataset.link}: not lit`);
    a.blur();
  }
  out.nLinks = n;
  return out;
}"""

# The named column lines, read from the page, and every block against them. Runs in the page.
GRID = r"""({ exempt, print }) => {
  const out = { fails: [], lines: {}, boxes: [] };
  const TOL = 1, near = (a, b) => Math.abs(a - b) <= TOL, f1 = (v) => (Math.round(v * 10) / 10).toFixed(1);
  const vw = document.documentElement.clientWidth;
  const NAMES = ["full-start", "wide-start", "body-start", "body-end", "wide-end", "full-end"];
  // The lines of a grid container: its content edge plus the resolved tracks, by line name.
  // In print there is no grid: one column of blocks, every line on a content edge.
  const read = (c) => {
    const cs = getComputedStyle(c), pos = {};
    let x = c.getBoundingClientRect().left + parseFloat(cs.borderLeftWidth) + parseFloat(cs.paddingLeft);
    if (print && cs.display === "block") {
      const r = c.getBoundingClientRect().right - parseFloat(cs.borderRightWidth) - parseFloat(cs.paddingRight);
      return Object.fromEntries(NAMES.map((n) => [n, n.endsWith("start") ? x : r]));
    }
    if (cs.display !== "grid") return pos;
    for (const t of cs.gridTemplateColumns.match(/\[[^\]]*\]|[\d.]+px/g) ?? []) {
      if (t[0] === "[") for (const n of t.slice(1, -1).split(/\s+/).filter(Boolean)) pos[n] = x;
      else x += parseFloat(t);
    }
    return pos;
  };
  const grids = [...document.querySelectorAll("body > header, main")];
  const L = read(grids.at(-1));
  out.lines = Object.fromEntries(NAMES.map((n) => [n, L[n]]));
  for (const c of grids) {
    const l = read(c), tag = c.tagName.toLowerCase();
    for (const n of NAMES) {
      if (!(n in l)) out.fails.push(`${tag}: no line ${n}`);
      else if (!near(l[n], L[n])) out.fails.push(`${tag}: line ${n} at ${f1(l[n])}, main's at ${f1(L[n])}`);
    }
  }
  if (NAMES.some((n) => !(n in L))) return out;

  // The lines are the ones the page's three lengths give.
  const px = (v) => {
    const d = document.createElement("div");
    d.style.cssText = `position:absolute;visibility:hidden;height:0;width:var(${v})`;
    document.body.append(d);
    const w = d.getBoundingClientRect().width;
    d.remove();
    return w;
  };
  const pad = px("--gutter"), over = px("--overhang"), measure = px("--measure");
  const body = print ? vw : Math.min(measure, vw - 2 * pad);
  const hang = print ? 0 : Math.max(0, Math.min(over, (vw - 2 * pad - body) / 2));
  const want = { "full-start": 0, "wide-start": (vw - body) / 2 - hang, "body-start": (vw - body) / 2,
                 "body-end": (vw + body) / 2, "wide-end": (vw + body) / 2 + hang, "full-end": vw };
  for (const n of NAMES) {
    if (!near(L[n], want[n])) out.fails.push(`line ${n} at ${f1(L[n])}, the page's lengths give ${f1(want[n])}`);
  }
  if (!near(L["body-start"] - L["wide-start"], L["wide-end"] - L["body-end"])) out.fails.push("wide is not centred on body");
  if (!print && L["wide-start"] - L["full-start"] < pad - TOL) out.fails.push("wide reaches into the side padding");

  const span = (c) => [L[`${c}-start`], L[`${c}-end`]];
  const name = (e) => e.id ? `#${e.id}` : e.dataset.fig ? `[data-fig=${e.dataset.fig}]`
    : `${e.parentElement.id ? `#${e.parentElement.id} > ` : ""}${e.tagName.toLowerCase()}${e.className ? `.${String(e.className).trim().split(/\s+/).join(".")}` : ""}`;
  const shown = (e) => e.getClientRects().length > 0;
  const classes = (e) => ["body", "wide", "full"].filter((c) => e.classList.contains(c));
  const inner = (e) => {  // content box; a display: contents parent has none, so its own parent's
    while (getComputedStyle(e).display === "contents") e = e.parentElement;
    const b = e.getBoundingClientRect(), cs = getComputedStyle(e);
    return [b.left + parseFloat(cs.paddingLeft), b.right - parseFloat(cs.paddingRight)];
  };
  // The rule of each exception: the edges it must have, or null where its class holds.
  const rules = {
    half(e) {
      const s = e.parentElement.getBoundingClientRect(), mid = (L["full-start"] + L["full-end"]) / 2;
      if (e.getBoundingClientRect().width > s.width - TOL) return span("full");
      return e === e.parentElement.firstElementChild ? [L["full-start"], mid] : [mid, L["full-end"]];
    },
    column: (e) => inner(e.parentElement),
    pair(e) {
      const p = e.parentElement, g = parseFloat(getComputedStyle(p).columnGap) || 0;
      if (e.getBoundingClientRect().width > p.getBoundingClientRect().width - TOL) return null;
      const mid = (L["wide-start"] + L["wide-end"]) / 2;
      return e === p.firstElementChild ? [L["wide-start"], mid - g / 2] : [mid + g / 2, L["wide-end"]];
    },
  };
  const FIGS = "[data-fig], .table, .pair, .split";
  const blocks = new Set([
    ...document.querySelectorAll("body > header > *, main > section, main > section:not(.split) > *"),
    ...document.querySelectorAll(FIGS), ...document.querySelectorAll(".body, .wide, .full"),
    ...document.querySelectorAll("[data-fig] > .canvas, [data-fig] > figcaption, .table > figcaption"),
  ]);
  const held = (e) => {  // the exception that holds this figure here, with the edges it gives
    const ex = exempt.find((x) => e.matches(x.sel)), edges = ex ? rules[ex.rule](e) : null;
    return edges ? [edges, `the exception "${ex.rule}"`] : null;
  };
  // The edges a block must have and what gives them, or [null, what is wrong with it].
  const expect = (e) => {
    const cs = classes(e), host = e.parentElement, hb = host.getBoundingClientRect();
    if (cs.length > 1) return [null, `has the classes ${cs.join(" and ")}`];
    if (held(e)) return held(e);
    if (e.matches("main > section:not(.split)")) return [span("full"), "full"];  // hands the lines on
    // A drawing follows its figure; so does the caption of a figure an exception holds.
    if (e.matches("[data-fig] > .canvas") || (e.matches("[data-fig] > figcaption") && held(host))) {
      return [[hb.left, hb.right], "its figure"];
    }
    if (e.matches("[data-fig] > figcaption, .table.wide > figcaption") && cs[0] !== "body") return [null, "is a caption under a figure on the grid without the class body"];
    if (e.matches(FIGS) && !cs.length) return [null, "has none of the classes body, wide, full"];
    return [span(cs[0] ?? "body"), cs[0] ?? "body"];
  };
  for (const e of blocks) {
    if (!shown(e)) continue;
    const b = e.getBoundingClientRect(), [edges, what] = expect(e), id = name(e);
    out.boxes.push({ name: id, on: what, left: b.left, right: b.right });
    if (!edges) { out.fails.push(`${id} ${what}`); continue; }
    if (!near(b.left, edges[0]) || !near(b.right, edges[1])) {
      out.fails.push(`${id} runs ${f1(b.left)} to ${f1(b.right)}; ${what} is ${f1(edges[0])} to ${f1(edges[1])}`);
    }
    for (const prop of ["width", "max-width", "min-width", "margin-left", "margin-right"]) {
      if (e.style.getPropertyValue(prop)) out.fails.push(`${id} sets its own ${prop}`);
    }
    if (getComputedStyle(e).transform !== "none") out.fails.push(`${id} is moved by a transform`);
  }
  return out;
}"""

PRINT = r"""() => {
  const out = { sticky: [], captions: [], controls: [] };
  for (const e of document.querySelectorAll("body *")) {
    const p = getComputedStyle(e).position;
    if ((p === "sticky" || p === "fixed") && e.getClientRects().length) out.sticky.push(`${e.tagName}.${e.className}`);
  }
  const pairs = [...document.querySelectorAll("[data-fig]")].map((e) => [e.dataset.fig, `cap-${e.dataset.fig}`]);
  for (const [f, c] of pairs) {
    const fe = document.querySelector(`[data-fig="${f}"]`);
    const fb = (fe.querySelector(":scope > .canvas") ?? fe).getBoundingClientRect();
    if (document.getElementById(c).getBoundingClientRect().top < fb.bottom - 1) out.captions.push(`${f}: caption not under it`);
  }
  for (const sel of [".controls", ".copy"]) {
    if ([...document.querySelectorAll(sel)].some((e) => getComputedStyle(e).display !== "none")) out.controls.push(`${sel} visible`);
  }
  return out;
}"""


def plate_iii_controls(page: Page) -> list[str]:
    """Play, pause and scrub change the frame; the button's icon is SVG, not text."""
    fails: list[str] = []
    page.eval_on_selector("#plate-iii", "e => e.scrollIntoView({block: 'start'})")
    page.wait_for_timeout(300)
    st = lambda: page.evaluate(  # noqa: E731
        "() => { const b = document.querySelector('.play'); return { v: document.querySelector('.scrub').value, "
        "label: b.getAttribute('aria-label'), text: b.textContent.trim(), "
        "icons: [...b.querySelectorAll('svg')].filter((s) => s.getClientRects().length).length }; }"
    )
    page.eval_on_selector(".play", "e => e.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    s0 = st()
    page.click(".play")
    page.wait_for_timeout(1200)
    s1 = st()
    page.click(".play")
    page.wait_for_timeout(200)
    s2 = st()
    page.focus(".scrub")
    page.keyboard.press("ArrowRight")
    page.keyboard.press("ArrowRight")
    s3 = st()
    for s in (s0, s1, s2):
        if s["text"]:
            fails.append(f"Plate III button has text {s['text']!r}")
        if s["icons"] != 1:
            fails.append(f"Plate III button shows {s['icons']} SVG icons")
    if not (s0["label"] == "Play" and s1["label"] == "Pause" and s2["label"] == "Play"):
        fails.append(f"Plate III button labels {s0['label']}, {s1['label']}, {s2['label']}")
    if s1["v"] == s0["v"]:
        fails.append("Plate III does not play")
    if int(s3["v"]) != int(s2["v"]) + 2:
        fails.append(f"Plate III arrows: {s2['v']} -> {s3['v']}")
    return fails


SHOTS = [("plate-i", "header"), ("table-1", "#table-1"), ("fig-1", "#fig-1"), ("fig-2", "#fig-2"),
         ("plate-iii", "#plate-iii"), ("plate-iii-side", "#plate-iii .side"), ("fig-4", "#fig-4"),
         ("plate-ii", "#plate-ii"), ("fig-5", "#fig-5"), ("fig-6", "#fig-6"), ("fig-7", "#fig-7"),
         ("fig-8", "#fig-8"), ("fig-9", "#fig-9"), ("fig-10", "#fig-10"), ("refs", "#s7")]  # fmt: skip
# The companion page: its title and abstract, each figure, Table 1 and the references.
SHOTS_COMPANION = [("title", "header"), ("fig-1", "#fig-1"), ("fig-2", "#fig-2"), ("fig-3", "#fig-3"),
                   ("table-1", "#table-1"), ("fig-4", "#fig-4"), ("refs", "#s7")]  # fmt: skip


def check(
    out: Path, engine: str, shots: Path | None, widths: list[tuple[int, int]], path: str = ""
) -> list[str]:
    """One page of the site: ``path`` under ``out`` ("" is the main paper)."""
    fails: list[str] = []
    main_page = path == ""
    with serve(out) as root, sync_playwright() as p:
        url = root + path
        browser = p.webkit.launch() if engine == "webkit" else launch(p)
        for w, h in widths:
            phone = w < PHONE
            ctx = browser.new_context(
                viewport={"width": w, "height": h}, reduced_motion="reduce",
                device_scale_factor=2 if phone else 1, is_mobile=phone, has_touch=phone,
            )  # fmt: skip
            problems: list[str] = []
            fetched: list[str] = []
            page = open_rendered(ctx, url, {"width": w, "height": h}, problems, fetched=fetched)
            page.wait_for_timeout(300)
            r = page.evaluate(PAGE, {"emoji": EMOJI, "narrow": w < 900})
            m = page.evaluate(MARKS)
            data = [json.loads((out / "data" / n).read_text()) for n in sorted(set(fetched))]
            derived = page.evaluate("globalThis.__derived.map((d) => d.x)")
            audit = page.evaluate(
                AUDIT, {"json": data, "derived": derived, "structure": STRUCTURE[path]}
            )
            g = page.evaluate(GRID, {"exempt": GRID_EXEMPT, "print": False})
            res = {
                "page": problems, "overflow": r["overflowX"], "grid lines": g["fails"],
                "emoji": r["emoji"],
                "sticky": r["sticky"], "blocks": r["blocks"], "table": r["table"],
                "columns": r["columns"], "gaps": r["gaps"], "links": r["links"],
                "text over marks": m["collisions"], "labels longer than their path": m["paths"],
                "text leaving its figure": m["leaves"], "number audit": audit["unmatched"],
            }  # fmt: skip
            if main_page and (w == 1440 or phone):
                res["Plate III controls"] = plate_iii_controls(page)
            bad = {k: v for k, v in res.items() if v}
            print(f"{engine} {w}px: {len(audit['numbers'])} numbers, {r['nLinks']} letter links, "
                  f"{'ok' if not bad else 'FAIL ' + ', '.join(bad)}")  # fmt: skip
            fails += [f"{w}px {k}: {x}" for k, v in bad.items() for x in v]
            if shots:
                shots.mkdir(parents=True, exist_ok=True)
                page.evaluate("window.scrollTo(0, 0)")
                for name, sel in SHOTS if main_page else SHOTS_COMPANION:
                    page.eval_on_selector(sel, "e => e.scrollIntoView({block: 'start'})")
                    page.evaluate("window.scrollBy(0, -16)")
                    page.wait_for_timeout(250)
                    page.screenshot(path=str(shots / f"{engine}_{w}_{name}.png" if not main_page
                                             else shots / f"{w}_{name}.png"))  # fmt: skip
                if not main_page and w in (390, 1440):
                    page.screenshot(path=str(shots / f"{engine}_{w}_page.png"), full_page=True)
            ctx.close()

        # Either side of the width at which wide meets body: the grid lines and the overflow.
        for w, h in GRID_WIDTHS if widths is WIDTHS else []:
            problems = []
            page = open_rendered(browser, url, {"width": w, "height": h}, problems)
            g = page.evaluate(GRID, {"exempt": GRID_EXEMPT, "print": False})
            over = page.evaluate(
                "document.documentElement.scrollWidth > document.documentElement.clientWidth"
            )
            res = {"page": problems, "overflow": ["page wider than the window"] if over else [],
                   "grid lines": g["fails"]}  # fmt: skip
            bad = {k: v for k, v in res.items() if v}
            print(
                f"{engine} {w}px: grid lines only, {'ok' if not bad else 'FAIL ' + ', '.join(bad)}"
            )
            fails += [f"{w}px {k}: {x}" for k, v in bad.items() for x in v]
            page.close()

        # Print: A4 from an 800 px window.
        ctx = browser.new_context(viewport={"width": 800, "height": 1000}, reduced_motion="reduce")
        problems = []
        page = open_rendered(ctx, url, {"width": 800, "height": 1000}, problems)
        if shots and engine != "webkit":  # Save as PDF from the 800 px window
            page.pdf(path=str(shots / "print_A4_from_800px.pdf"), format="A4",
                     print_background=True, prefer_css_page_size=True)  # fmt: skip
        # The checks: the print layout at the A4 text width.
        page.emulate_media(media="print")
        page.set_viewport_size({"width": 680, "height": 1000})
        page.wait_for_timeout(800)
        pr = page.evaluate(PRINT)
        m = page.evaluate(MARKS)
        g = page.evaluate(GRID, {"exempt": GRID_EXEMPT, "print": True})
        res = {"page": problems, "sticky": pr["sticky"], "captions": pr["captions"],
               "controls": pr["controls"], "grid lines": g["fails"],
               "text over marks": m["collisions"],
               "labels longer than their path": m["paths"], "text leaving its figure": m["leaves"]}  # fmt: skip
        bad = {k: v for k, v in res.items() if v}
        print(f"{engine} print: {'ok' if not bad else 'FAIL ' + ', '.join(bad)}")
        fails += [f"print {k}: {x}" for k, v in bad.items() for x in v]
        ctx.close()

        # The page without JavaScript: no emoji either.
        ctx = browser.new_context(
            java_script_enabled=False, viewport={"width": 1440, "height": 900}
        )
        page = ctx.new_page()
        page.goto(url)
        page.wait_for_load_state("load")
        r = page.evaluate(PAGE, {"emoji": EMOJI, "narrow": False})  # only its emoji count here
        print(
            f"{engine} without JS: {'ok' if not r['emoji'] else 'FAIL emoji ' + ' '.join(r['emoji'])}"
        )
        fails += [f"without JS emoji: {x}" for x in r["emoji"]]
        ctx.close()
        browser.close()
    return fails


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out", type=Path, help="the site to check: _site (built) or site")
    ap.add_argument("--engine", choices=["chromium", "webkit"], default="chromium")
    ap.add_argument("--shots", type=Path, help="save screenshots of every plate and figure here")
    ap.add_argument("--page", choices=PAGES, help="check one page only (default: every page)")
    ap.add_argument("--widths", help="comma-separated widths (default: all)")
    a = ap.parse_args()
    widths = WIDTHS
    if a.widths:
        want = {int(x) for x in a.widths.split(",")}
        widths = [wh for wh in WIDTHS if wh[0] in want] + [
            (x, 900) for x in want - {w for w, _ in WIDTHS}
        ]
    fails: list[str] = []
    for path in [a.page] if a.page is not None else PAGES:
        print(f"--- {path or 'main page'}")
        fails += [
            f"[{path or 'main'}] {f}"
            for f in check(a.out.resolve(), a.engine, a.shots, widths, path)
        ]
    if fails:
        print("\nFAILED:\n  " + "\n  ".join(fails))
        raise SystemExit(1)
    print("\nall checks pass")


if __name__ == "__main__":
    main()
