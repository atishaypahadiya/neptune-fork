// Keeps text clear of the marks in a figure (docs/DESIGN.md, Acceptance 2; checked by
// tools/check_layout.py): no text box crosses a stroke, a marker or another text.
//
// - Point letters find a free spot around their point: the drawn offset first, then nearby
//   positions in order of distance and turn.
// - Lines that only guide the eye (gridlines, the launch line, reference lines: class
//   "yield") break where they would cross a text.
// - A text's own leader line is not an obstacle to it (data-leader / data-leader-of).
//
// Geometry is in the data SVG's drawing space (w x h, as drawn).

const PAD = 1.5; // px kept clear around a text

/** Matrix from an element's user space to the drawing space of its data SVG. */
function toSvg(svg, el) {
  const a = svg.getScreenCTM(), b = el.getScreenCTM();
  return a && b ? a.inverse().multiply(b) : null;
}
const apply = (m, x, y) => { const p = new DOMPoint(x, y).matrixTransform(m); return [p.x, p.y]; };
const quadOf = (m, r) => [apply(m, r.x, r.y), apply(m, r.x + r.width, r.y),
  apply(m, r.x + r.width, r.y + r.height), apply(m, r.x, r.y + r.height)];
const boxOf = (q) => {
  const xs = q.map((p) => p[0]), ys = q.map((p) => p[1]);
  return [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)];
};
const meets = (a, b) => a[0] < b[2] && b[0] < a[2] && a[1] < b[3] && b[1] < a[3];
const grow = (b, d) => [b[0] - d, b[1] - d, b[2] + d, b[3] + d];

/** Does the segment (x1, y1)-(x2, y2) enter box b? Returns the entry and exit parameters. */
function clip(x1, y1, x2, y2, b) {
  let t0 = 0, t1 = 1;
  const dx = x2 - x1, dy = y2 - y1;
  for (const [p, q] of [[-dx, x1 - b[0]], [dx, b[2] - x1], [-dy, y1 - b[1]], [dy, b[3] - y1]]) {
    if (p === 0) { if (q < 0) return null; continue; }
    const r = q / p;
    if (p < 0) { if (r > t1) return null; if (r > t0) t0 = r; }
    else { if (r < t0) return null; if (r < t1) t1 = r; }
  }
  return [t0, t1];
}

/** Separating axis test for two convex quads. */
function sat(A, B) {
  for (const Q of [A, B]) for (let i = 0; i < 4; i++) {
    const [x1, y1] = Q[i], [x2, y2] = Q[(i + 1) % 4], n = [y1 - y2, x2 - x1];
    const pa = A.map((p) => p[0] * n[0] + p[1] * n[1]), pb = B.map((p) => p[0] * n[0] + p[1] * n[1]);
    if (Math.max(...pa) <= Math.min(...pb) || Math.max(...pb) <= Math.min(...pa)) return false;
  }
  return true;
}

/** Straight segments of a path in its own user space (M, L, H, V, Z; else sampled). */
function pathSegs(el) {
  const d = el.getAttribute("d") ?? "";
  const toks = d.match(/[MLHVZmlhvz]|[-+]?(?:\d*\.\d+|\d+\.?)(?:e[-+]?\d+)?/gi) ?? [];
  if (/[^MLHVZmlhvz\d\s,.eE+-]/.test(d) || toks.some((t) => /[mlhv]/.test(t))) {
    const len = el.getTotalLength(), n = Math.max(2, Math.ceil(len / 2)), pts = [];
    for (let i = 0; i <= n; i++) { const p = el.getPointAtLength((len * i) / n); pts.push([p.x, p.y]); }
    return pts.slice(1).map((p, i) => [...pts[i], ...p]);
  }
  const segs = [];
  let cmd = "M", x = 0, y = 0, sx = 0, sy = 0, i = 0;
  while (i < toks.length) {
    if (/[A-Za-z]/.test(toks[i])) { cmd = toks[i++].toUpperCase(); if (cmd === "Z") { segs.push([x, y, sx, sy]); x = sx; y = sy; } continue; }
    if (cmd === "M") { x = sx = +toks[i]; y = sy = +toks[i + 1]; i += 2; cmd = "L"; continue; }
    let nx = x, ny = y;
    if (cmd === "L") { nx = +toks[i]; ny = +toks[i + 1]; i += 2; }
    else if (cmd === "H") { nx = +toks[i]; i += 1; }
    else if (cmd === "V") { ny = +toks[i]; i += 1; }
    else { i += 1; continue; }
    segs.push([x, y, nx, ny]);
    x = nx; y = ny;
  }
  return segs;
}

const cache = new WeakMap();
const visible = (e) => (e.checkVisibility ? e.checkVisibility({ opacityProperty: true, visibilityProperty: true }) : true);

/**
 * The marks of a figure's data SVG: strokes as segments (with half their width), small marks
 * and filled shapes as boxes, texts as quads. `skip`: texts to leave out (those being placed).
 */
export function marks(svg, { skip = new Set() } = {}) {
  const out = { segs: [], boxes: [], texts: [] };
  for (const e of svg.querySelectorAll("path, line, polyline, polygon, circle, ellipse, rect")) {
    if (e.closest("defs, text") || e.matches(".hit, .glow") || !visible(e)) continue;
    const m = toSvg(svg, e);
    if (!m) continue;
    const key = [e.getAttribute("d"), e.getAttribute("cx"), e.getAttribute("cy"), e.getAttribute("r"),
      e.getAttribute("x"), e.getAttribute("y"), m.e, m.f, m.a].join("|");
    let c = cache.get(e);
    if (!c || c.key !== key) {
      const cs = getComputedStyle(e);
      const sw = cs.stroke !== "none" ? parseFloat(cs.strokeWidth) || 0 : 0;
      const filled = cs.fill !== "none" && cs.fill !== "transparent" && !/,\s*0\)$/.test(cs.fill);
      if (!sw && !filled) { cache.set(e, (c = { key, none: true })); continue; }
      const bb = e.getBBox(), scale = Math.hypot(m.a, m.b), hw = (sw * scale) / 2;
      const small = bb.width * scale < 16 && bb.height * scale < 16;
      if (small || e.tagName !== "path") {
        c = { key, box: grow(boxOf(quadOf(m, bb)), hw) };
      } else {
        c = { key, segs: pathSegs(e).map(([a, b, cc, d]) => [...apply(m, a, b), ...apply(m, cc, d), hw]) };
      }
      c.leaderOf = e.dataset.leaderOf ?? null;
      c.soft = e.classList.contains("yield");
      cache.set(e, c);
    }
    if (c.none) continue;
    if (c.box) out.boxes.push({ b: c.box, el: e, leaderOf: c.leaderOf, soft: c.soft });
    else for (const s of c.segs) out.segs.push({ s, el: e, leaderOf: c.leaderOf, soft: c.soft });
  }
  for (const t of svg.querySelectorAll("text")) {
    if (skip.has(t) || !visible(t) || !t.textContent.trim()) continue;
    for (const q of textQuads(svg, t)) out.texts.push({ q, b: boxOf(q), el: t });
  }
  // HTML text laid over the plate ("Plate II.", Plate I's title block), line by line.
  const inv = svg.getScreenCTM()?.inverse();
  const fig = svg.closest("[data-fig]");
  if (inv && fig) {
    for (const e of fig.querySelectorAll(".plate-label, .titleblock h1, .titleblock p")) {
      if (!e.getClientRects().length) continue;
      const r = document.createRange();
      r.selectNodeContents(e);
      for (const c of r.getClientRects()) {
        if (!c.width) continue;
        const q = [[c.left, c.top], [c.right, c.top], [c.right, c.bottom], [c.left, c.bottom]]
          .map(([x, y]) => apply(inv, x, y));
        out.texts.push({ q, b: boxOf(q), el: e });
      }
    }
  }
  return out;
}

/** The drawing's size (mount sets the viewBox to it). */
const drawing = (svg) => {
  const v = svg.viewBox?.baseVal;
  return v && v.width ? [v.x, v.y, v.x + v.width, v.y + v.height] : null;
};
const inside = (svg, quads) => {
  const d = drawing(svg);
  return !d || quads.every((q) => { const b = boxOf(q); return b[0] >= d[0] + 2 && b[1] >= d[1] + 2 && b[2] <= d[2] - 2 && b[3] <= d[3] - 2; });
};

/** A text's boxes in drawing space: one box, or one per glyph along a path. */
export function textQuads(svg, t) {
  const m = toSvg(svg, t);
  if (!m) return [];
  if (t.querySelector("textPath")) {
    const out = [];
    for (let i = 0; i < t.getNumberOfChars(); i++) {
      const r = t.getExtentOfChar(i);
      if (r.width > 0 && r.height > 0) out.push(quadOf(m, r));
    }
    return out;
  }
  return [quadOf(m, t.getBBox())];
}

/** Count what the quads of a text would cross. `soft`: count yielding lines too. */
export function crossings(quads, obs, { soft = true, leader = null, pad = PAD } = {}) {
  let n = 0;
  for (const q of quads) {
    const b = grow(boxOf(q), pad);
    for (const { s, soft: sf, leaderOf } of obs.segs) {
      if ((sf && !soft) || (leader && leaderOf === leader)) continue;
      const hb = grow(b, s[4]);
      if (Math.max(s[0], s[2]) < hb[0] || Math.min(s[0], s[2]) > hb[2] ||
          Math.max(s[1], s[3]) < hb[1] || Math.min(s[1], s[3]) > hb[3]) continue;
      if (clip(s[0], s[1], s[2], s[3], hb)) n++;
    }
    for (const { b: o, soft: sf, leaderOf } of obs.boxes) {
      if ((sf && !soft) || (leader && leaderOf === leader)) continue;
      if (meets(b, o)) n++;
    }
    for (const { q: o, b: ob } of obs.texts) if (meets(b, ob) && sat(q, o)) n++;
  }
  return n;
}

/** Candidate offsets around a point: the drawn one first, then by distance and turn. */
function candidates(dx, dy) {
  const r0 = Math.hypot(dx, dy), a0 = Math.atan2(dy, dx);
  const out = [[dx, dy, null]];
  const radii = [...new Set([Math.max(r0, 11), 12, 15, 19, 24, 30, 37, 45].map((r) => Math.round(r)))]
    .filter((r) => r >= Math.min(r0, 11)).sort((a, b) => a - b);
  const turns = Array.from({ length: 16 }, (_, k) => (k * Math.PI) / 8)
    .map((a) => [a, Math.abs(Math.atan2(Math.sin(a - a0), Math.cos(a - a0)))])
    .sort((x, y) => x[1] - y[1]).map(([a]) => a);
  for (const r of radii) for (const a of turns) {
    const c = Math.cos(a);
    out.push([r * c, r * Math.sin(a), c > 0.38 ? "start" : c < -0.38 ? "end" : "middle"]);
  }
  return out;
}

/**
 * Put each point letter where it crosses nothing: its drawn offset if that is clear, else the
 * nearest clear spot. Lines that yield are avoided where possible and break where not.
 * `letters`: the letters to place (default: all); the others stay where they are.
 */
export function placeLetters(svg, letters = null) {
  const list = letters ?? [...svg.querySelectorAll("[data-pt] > text.letter")];
  const obs = marks(svg, { skip: new Set(list) });
  for (const t of list) {
    if (!visible(t)) continue;
    if (!t.dataset.dx) {
      t.dataset.dx = t.getAttribute("x");
      t.dataset.dy = t.getAttribute("y");
      t.dataset.anchor = t.getAttribute("text-anchor") ?? "start";
    }
    const cands = candidates(+t.dataset.dx, +t.dataset.dy);
    const set = ([x, y, anchor]) => {
      t.setAttribute("x", x);
      t.setAttribute("y", y);
      t.setAttribute("text-anchor", anchor ?? t.dataset.anchor);
    };
    // Keep last frame's spot while it stays clear (Plate III moves its letters every frame).
    const order = t.__cand !== undefined ? [t.__cand, ...cands.keys()] : [...cands.keys()];
    let best = null;
    for (const soft of [true, false]) {
      for (const k of order) {
        set(cands[k]);
        const qs = textQuads(svg, t);
        const n = crossings(qs, obs, { soft }) + (inside(svg, qs) ? 0 : 1000); // stays in the figure
        if (n === 0) { best = k; break; }
        if (!soft && (best === null || n < best.n)) best = { k, n };
      }
      if (typeof best === "number") break;
    }
    const k = typeof best === "number" ? best : best.k;
    set(cands[k]);
    t.__cand = k;
    for (const q of textQuads(svg, t)) obs.texts.push({ q, b: boxOf(q), el: t });
  }
}

/** Break every yielding line where it crosses a text (with a little air). */
export function yieldLines(svg) {
  const texts = [...svg.querySelectorAll("text")].filter((t) => visible(t) && t.textContent.trim())
    .flatMap((t) => textQuads(svg, t).map((q) => grow(boxOf(q), 3)));
  for (const p of svg.querySelectorAll("path.yield")) {
    const m = toSvg(svg, p);
    if (!m || m.a !== 1 || m.d !== 1 || m.b || m.c) continue; // drawn unscaled
    const d0 = p.dataset.d0 ?? p.getAttribute("d");
    p.dataset.d0 = d0;
    const parts = [];
    for (const [x1, y1, x2, y2] of pathSegs(withD(p, d0))) {
      const X1 = x1 + m.e, Y1 = y1 + m.f, X2 = x2 + m.e, Y2 = y2 + m.f;
      const cuts = texts.map((b) => clip(X1, Y1, X2, Y2, b)).filter(Boolean).sort((a, b) => a[0] - b[0]);
      let t = 0;
      for (const [a, b] of cuts) {
        if (a > t) parts.push([t, a, x1, y1, x2, y2]);
        t = Math.max(t, b);
      }
      if (t < 1) parts.push([t, 1, x1, y1, x2, y2]);
    }
    const pt = (t, x1, y1, x2, y2) => `${x1 + t * (x2 - x1)},${y1 + t * (y2 - y1)}`;
    p.setAttribute("d", parts.filter(([a, b]) => b - a > 1e-6)
      .map(([a, b, ...s]) => `M${pt(a, ...s)}L${pt(b, ...s)}`).join("") || "M0,0");
    p.style.display = parts.length ? "" : "none";
  }
}

/** A detached path with the given d, for parsing. */
function withD(p, d) {
  const q = p.cloneNode(false);
  q.setAttribute("d", d);
  return q;
}

/** Everything in one pass, after a figure is drawn. */
export function settle(svg) {
  placeLetters(svg);
  yieldLines(svg);
}
