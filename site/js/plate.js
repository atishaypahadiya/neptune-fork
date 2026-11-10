// Shared drawing for every figure: the textured ground (archive layer), lettered points, the
// text <-> figure letter links, the hover probe, and redraw on resize.
//
// The ground is TEXTURE: SVG filters with fixed seeds (index.html #tex-*). It encodes nothing.
// It sits in its own SVG behind the data layer, so labels stay sharp and redraws of the data
// layer (Plate III's animation) never re-render the filters.

import * as d3 from "d3";

const SVGNS = "http://www.w3.org/2000/svg";

/** Insert the ground SVG behind a figure's data SVG. `mark`: the hairline plate mark (Plate I). */
export function ground(fig, { mark = false } = {}) {
  if (fig.querySelector(":scope > svg.ground")) return;
  const g = document.createElementNS(SVGNS, "svg");
  g.setAttribute("class", "ground");
  g.setAttribute("aria-hidden", "true");
  g.setAttribute("focusable", "false");
  const rect = (attrs) => {
    const r = document.createElementNS(SVGNS, "rect");
    for (const [k, v] of Object.entries({ x: 0, y: 0, width: "100%", height: "100%", ...attrs })) {
      r.setAttribute(k, v);
    }
    g.appendChild(r);
    return r;
  };
  rect({ fill: "var(--prussian)" });
  rect({ filter: "url(#tex-mottle)" });
  rect({ filter: "url(#tex-cloud)", opacity: 0.55 });
  rect({ filter: "url(#tex-grain)", opacity: 0.07, style: "mix-blend-mode: screen" });
  rect({ filter: "url(#tex-dust)", opacity: 0.5 });
  rect({ fill: "url(#tex-vignette)" });
  if (mark) {
    const r = rect({ x: 14, y: 14, fill: "none", stroke: "var(--line)", "stroke-opacity": 0.3,
                     "stroke-width": 1 });
    r.style.width = "calc(100% - 28px)";
    r.style.height = "calc(100% - 28px)";
  }
  fig.insertBefore(g, fig.firstChild);
}

/**
 * Mount a figure: ground + data layer, drawn by `draw(svg, w, h)` now and on every resize.
 * Returns a function that redraws on demand.
 */
export function mount(fig, draw, opts = {}) {
  const box = fig.querySelector(":scope > .canvas") ?? fig;
  ground(box, opts);
  const svg = d3.select(box.querySelector(":scope > svg.fig"));
  let last = "";
  const redraw = (force = false) => {
    const { width: w, height: h } = box.getBoundingClientRect();
    const key = `${Math.round(w)}x${Math.round(h)}`;
    if (!force && key === last) return;
    last = key;
    svg.attr("viewBox", `0 0 ${w} ${h}`).selectAll("*").remove();
    // Print keeps the drawn aspect (site.css @media print) and the viewBox scales the drawing.
    if (w > 0 && h > 0) box.style.setProperty("--ar", String(w / h));
    draw(svg, w, h);
    refreshLinks();
  };
  new ResizeObserver(() => redraw()).observe(box);
  redraw(true);
  return () => redraw(true);
}

/** A lettered point: glow (shown on highlight), mark, italic letter and a hit target. */
export function point(g, { id, x, y, letter, dx = 9, dy = -9, anchor = "start", mark = "dot",
                           r = 3.2, title = "" }) {
  const p = g.append("g").attr("data-pt", id).attr("transform", `translate(${x},${y})`);
  p.append("circle").attr("class", "glow").attr("r", 12);
  if (mark === "dot") p.append("circle").attr("class", "dot").attr("r", r);
  else if (mark === "open") p.append("circle").attr("class", "open").attr("r", r + 1.5);
  if (letter) {
    p.append("text").attr("class", "letter").attr("x", dx).attr("y", dy)
      .attr("text-anchor", anchor).attr("dominant-baseline", "middle").text(letter);
  }
  p.append("circle").attr("class", "hit").attr("r", 16);
  if (title) p.append("title").text(title);
  return p;
}

/** A stippled body: a cluster of dots. `pts` are pixel positions from data (or, for the Sun, a
 * fixed-seed texture that encodes nothing). */
export function stipple(g, pts, { r = 0.9, cls = "dot" } = {}) {
  g.selectAll(null).data(pts).join("circle").attr("class", cls)
    .attr("cx", (p) => p[0]).attr("cy", (p) => p[1]).attr("r", r);
}

/** TEXTURE: the Sun's stipple, a fixed-seed disc of dots. Encodes nothing. */
export function sunTexture(cx, cy, radius, n = 70) {
  const rnd = d3.randomLcg(0.1846);
  const out = [];
  for (let i = 0; i < n; i++) {
    const a = 2 * Math.PI * rnd(), rr = radius * Math.sqrt(rnd());
    out.push([cx + rr * Math.cos(a), cy + rr * Math.sin(a)]);
  }
  return out;
}

/** Dotted orbit path from points already in pixels. */
export function orbit(g, pts, { cls = "ln orbit", closed = false } = {}) {
  const line = closed ? d3.line().curve(d3.curveLinearClosed) : d3.line();
  return g.append("path").attr("class", cls).attr("d", line(pts));
}

/* ---------- prose: text with point letters ---------- */

/**
 * Fill an element from a template. `{X@fig}` becomes the point letter X, linked to point X of
 * figure `fig` (a focusable element; hover or focus highlights the point). `{i:word}` is an
 * italic variable. Everything else is text. No HTML from data is parsed.
 */
export function fill(el, template) {
  el.replaceChildren();
  const re = /\{([A-Za-z]{1,2})@([a-z0-9-]+)\}|\{i:([^}]+)\}/g;
  let at = 0, m;
  while ((m = re.exec(template))) {
    if (m.index > at) el.append(template.slice(at, m.index));
    if (m[3]) {
      const i = document.createElement("i");
      i.textContent = m[3];
      el.append(i);
    } else {
      const [, letter, fig] = m;
      const a = document.createElement("a");
      a.className = "pt";
      a.textContent = letter;
      a.href = `#${fig === "plate-i" ? "plate-i" : fig}`;
      a.dataset.link = `${fig}:${letter}`;
      a.setAttribute("aria-describedby", `cap-${fig}`);
      a.setAttribute("aria-label", `point ${letter} in ${figName(fig)}`);
      el.append(a);
    }
    at = re.lastIndex;
  }
  if (at < template.length) el.append(template.slice(at));
}

const figName = (fig) => fig.replace("plate-", "Plate ").replace("fig-", "Fig. ")
  .replace(/Plate (i+)$/, (_, r) => `Plate ${r.toUpperCase()}`);

/** Wire every prose letter to its point (called after each redraw). */
let wired = false;
export function refreshLinks() {
  if (wired) return;
  wired = true;
  const find = (link) => {
    const [fig, letter] = link.split(":");
    return document.querySelector(`[data-fig="${fig}"] [data-pt="${letter}"]`);
  };
  const on = (e, state) => {
    const a = e.target.closest?.(".pt");
    if (!a) return;
    find(a.dataset.link)?.classList.toggle("on", state);
  };
  for (const [ev, s] of [["mouseover", true], ["mouseout", false], ["focusin", true], ["focusout", false]]) {
    document.addEventListener(ev, (e) => on(e, s));
  }
  document.addEventListener("click", (e) => {
    const a = e.target.closest?.(".pt");
    if (!a) return;
    const fig = document.querySelector(`[data-fig="${a.dataset.link.split(":")[0]}"]`);
    if (!fig) return;
    const r = fig.getBoundingClientRect();
    if (r.top >= 0 && r.bottom <= innerHeight) e.preventDefault();
  });
}

/** Hover probe: a quiet text tip inside a figure. Nothing appears only here. */
export function probe(fig) {
  const box = fig.querySelector(":scope > .canvas") ?? fig;
  const tip = box.querySelector(":scope > .tip");
  return {
    show(x, y, text) {
      if (!tip) return;
      tip.textContent = text;
      tip.hidden = false;
      const w = box.clientWidth;
      tip.style.left = `${Math.min(Math.max(8, x + 14), w - tip.offsetWidth - 8)}px`;
      tip.style.top = `${Math.max(8, y - tip.offsetHeight - 10)}px`;
    },
    hide() { if (tip) tip.hidden = true; },
  };
}

/** Show a figure's error state (data missing or inconsistent). */
export function fail(fig, err) {
  const p = fig.querySelector(".error");
  if (p) {
    p.hidden = false;
    p.textContent = `This figure could not be drawn: ${err.message}`;
  }
  console.error(err);
}
