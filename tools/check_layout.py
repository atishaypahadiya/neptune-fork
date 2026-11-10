"""Layout, print and interaction acceptance checks for the site (docs/DESIGN.md, Acceptance).

Loads the page with JavaScript at 360, 390, 430, 768, 899, 900, 1024 and 1440 px (phones as touch
devices) and in print (A4 from an 800 px window), and fails on any of:

- **Text over marks:** in every plate and figure, a text box (SVG text, the lettering along
  orbits glyph by glyph, and the HTML text laid over a plate) that intersects a stroke, a
  marker or another text box. The one exception is a text's own leader line
  (``data-leader`` on the text, ``data-leader-of`` on the line).
- **Labels along a path** (``textPath``) longer than their path, so that their ends would be
  cut off. The length is measured with letter-spacing, which some engines leave out.
- Text leaving its figure; horizontal overflow of the page.
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
from prerender_site import launch, open_rendered, serve

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
  pairs([c.querySelector(".play"), c.querySelector("label")], "Plate III controls");
  pairs([c.querySelector(".play"), c.querySelector("label > span") ?? c.querySelector("label"), c.querySelector(".scrub")].filter(Boolean), "Plate III controls");

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

PRINT = r"""() => {
  const out = { sticky: [], captions: [], controls: [] };
  for (const e of document.querySelectorAll("body *")) {
    const p = getComputedStyle(e).position;
    if ((p === "sticky" || p === "fixed") && e.getClientRects().length) out.sticky.push(`${e.tagName}.${e.className}`);
  }
  const pairs = [["plate-i", "cap-plate-i"], ["plate-ii", "cap-plate-ii"], ["plate-iii", "cap-plate-iii"],
    ...["fig-1", "fig-2", "fig-3", "fig-4", "fig-5", "fig-6"].map((f) => [f, `cap-${f}`])];
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
         ("plate-ii", "#plate-ii"), ("fig-5", "#fig-5"), ("fig-6", "#fig-6"), ("refs", "#s7")]  # fmt: skip


def check(out: Path, engine: str, shots: Path | None, widths: list[tuple[int, int]]) -> list[str]:
    fails: list[str] = []
    with serve(out) as url, sync_playwright() as p:
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
            audit = page.evaluate(AUDIT, {"json": data, "derived": derived})
            res = {
                "page": problems, "overflow": r["overflowX"], "emoji": r["emoji"],
                "sticky": r["sticky"], "blocks": r["blocks"], "table": r["table"],
                "columns": r["columns"], "gaps": r["gaps"], "links": r["links"],
                "text over marks": m["collisions"], "labels longer than their path": m["paths"],
                "text leaving its figure": m["leaves"], "number audit": audit["unmatched"],
            }  # fmt: skip
            if w == 1440 or phone:
                res["Plate III controls"] = plate_iii_controls(page)
            bad = {k: v for k, v in res.items() if v}
            print(f"{engine} {w}px: {len(audit['numbers'])} numbers, {r['nLinks']} letter links, "
                  f"{'ok' if not bad else 'FAIL ' + ', '.join(bad)}")  # fmt: skip
            fails += [f"{w}px {k}: {x}" for k, v in bad.items() for x in v]
            if shots:
                shots.mkdir(parents=True, exist_ok=True)
                page.evaluate("window.scrollTo(0, 0)")
                for name, sel in SHOTS:
                    page.eval_on_selector(sel, "e => e.scrollIntoView({block: 'start'})")
                    page.evaluate("window.scrollBy(0, -16)")
                    page.wait_for_timeout(250)
                    page.screenshot(path=str(shots / f"{w}_{name}.png"))
            ctx.close()

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
        res = {"page": problems, "sticky": pr["sticky"], "captions": pr["captions"],
               "controls": pr["controls"], "text over marks": m["collisions"],
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
    ap.add_argument("--widths", help="comma-separated widths (default: all)")
    a = ap.parse_args()
    widths = WIDTHS
    if a.widths:
        want = {int(x) for x in a.widths.split(",")}
        widths = [wh for wh in WIDTHS if wh[0] in want] + [
            (x, 900) for x in want - {w for w, _ in WIDTHS}
        ]
    fails = check(a.out.resolve(), a.engine, a.shots, widths)
    if fails:
        print("\nFAILED:\n  " + "\n  ".join(fails))
        raise SystemExit(1)
    print("\nall checks pass")


if __name__ == "__main__":
    main()
