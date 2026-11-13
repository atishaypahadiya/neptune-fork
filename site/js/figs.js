// Figs. 1, 2, 4, 5 and 6: small plates (white line on the cyanotype ground). Every mark is drawn
// from site/data/; axis tick values are derived from the data's range and registered as such.

import * as d3 from "d3";
import { mount, point, probe, stipple } from "./plate.js";
import { fmt, drv, clock, dayMonth, duration, profileSegments } from "./data.js";
import { timelineFacts, SENS } from "./copy.js";

const L = d3.line();
const minus = (t) => `${t < 0 ? "−" : ""}${Math.abs(t)}`;
const years = (y) => (y[0] === y[1] ? `${y[0]}` : `${y[0]}–${y[1]}`);
// Tooltips follow a mouse and stay put after a tap: a touch "leave" comes with every lift.
const leave = (tip) => (e) => { if (e.pointerType === "mouse") tip.hide(); };

function xAxis(g, x, y0, ticks, fmtT = String) {
  for (const t of ticks) {
    drv(t, "axis tick");
    g.append("path").attr("class", "ln ctx").attr("d", L([[x(t), y0], [x(t), y0 + 5]]));
    g.append("text").attr("class", "tick").attr("x", x(t)).attr("y", y0 + 22)
      .attr("text-anchor", "middle").text(fmtT(t));
  }
}

/* ---------------- Fig. 1: the observations over time ---------------- */

export function fig1(fig, D) {
  const rows = D.table1.rows;
  const meta = new Map(D.dataset.observations.map((o) => [o.id, o]));
  const obs = D.table1.observations.map((o) => ({ ...o, ...meta.get(o.id), src: o.source }));
  // Each observation's place up or down within its row: the golden-ratio sequence in time
  // order, so neighbours in time sit apart. It spreads dense series; it encodes nothing.
  const kth = new Map();
  for (const o of obs) {
    const k = kth.get(o.src) ?? 0;
    kth.set(o.src, k + 1);
    o.spread = (((k * 0.6180339887) % 1) * 2) - 1;
  }
  // A row's name: the source, its observatory unless the name holds it, and its years when two
  // rows share a name (the two Paris series).
  const label = (r) => {
    const twin = rows.filter((x) => x.source === r.source).length > 1;
    return [r.source, ...(r.source.includes(r.observatory) ? [] : [r.observatory]),
      ...(twin ? [years(r.years)] : [])].join(", ");
  };
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: phone ? 14 : 190, r: 18, t: phone ? 34 : 22, b: 34 };
    const x = d3.scaleLinear().domain(d3.extent(obs, (o) => o.year)).nice().range([mg.l, w - mg.r]);
    const y = d3.scalePoint().domain(rows.map((r) => r.key)).range([mg.t, h - mg.b]).padding(0.5);
    const amp = Math.min(phone ? 10 : 8, y.step() * 0.2), r0 = phone ? 2 : 2.4;
    const cy = (o) => y(o.src) + amp * o.spread;
    // Dots shrink with a row's density: about the room each one has in its band.
    const rad = new Map(rows.map((r) => {
      const xs = obs.filter((o) => o.src === r.key).map((o) => x(o.year));
      const room = ((Math.max(...xs) - Math.min(...xs) + 2 * r0) * 2 * amp) / xs.length;
      return [r.key, Math.max(1.1, Math.min(r0, 0.5 * Math.sqrt(room)))];
    }));
    const rOf = (o) => (o.spot_checked ? Math.max(rad.get(o.src), 2) : rad.get(o.src));
    const g = svg.append("g");
    for (const r of rows) {
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[mg.l, y(r.key)], [w - mg.r, y(r.key)]]));
      // On phones the name sits above its row, clear of the spread dots.
      if (phone) g.append("text").attr("class", "lbl row").attr("x", mg.l).attr("y", y(r.key) - amp - 9).text(label(r));
      else g.append("text").attr("class", "lbl row").attr("x", mg.l - 12).attr("y", y(r.key) + 4).attr("text-anchor", "end").text(label(r));
    }
    xAxis(g, x, h - mg.b + 4, x.ticks(phone ? 4 : 8));
    const dots = g.append("g");
    dots.selectAll(null).data(obs).join("circle")
      .attr("class", (o) => (o.spot_checked ? "open" : "dot"))
      .attr("cx", (o) => x(o.year)).attr("cy", cy).attr("r", rOf);
    const q = d3.quadtree(obs, (o) => x(o.year), cy);
    svg.on("pointermove pointerdown", (e) => {
      const [px, py] = d3.pointer(e);
      const o = q.find(px, py, 24);
      if (!o) return tip.hide();
      const part = (k, n) => (o[k] === undefined ? "" : `${n} ${minus(o[k])}″`);
      tip.show(x(o.year), cy(o), `${o.observer}, ${o.observatory}, ${fmt(o.year, 3)}. O − JPL: ` +
        [part("ra", "RA"), part("dec", "Dec")].filter(Boolean).join(", ") +
        (o.spot_checked ? ". Checked by hand." : "."));
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 2: the pre-registration timeline ---------------- */

const SX = (s) => Math.sign(s) * Math.log10(1 + Math.abs(s) / 20);
const TICKS = [[-36000, "−10 h"], [-3600, "−1 h"], [-600, "−10 min"], [-60, "−1 min"], [0, "launch"],
  [60, "+1 min"], [600, "+10 min"], [3600, "+1 h"], [36000, "+10 h"]];

/** A Fig. 2 row label as one line, or two if it is longer than 28 characters. */
function rowLines(title) {
  if (title.length <= 28) return [title];
  const k = title.indexOf(" (");
  if (k > 0) return [title.slice(0, k), title.slice(k + 1)];
  const mid = title.length / 2;
  let best = -1;
  for (let i = title.indexOf(" "); i >= 0; i = title.indexOf(" ", i + 1)) {
    if (best < 0 || Math.abs(i - mid) < Math.abs(best - mid)) best = i;
  }
  return best < 0 ? [title] : [title.slice(0, best), title.slice(best + 1)];
}

export function fig2(fig, D) {
  const F = timelineFacts(D.timeline);
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: phone ? 14 : 300, r: 20, t: phone ? 34 : 20, b: 36 };
    const vals = F.flatMap((f) => [...f.pre.flatMap((p) => [p.dCommit, p.dPush]), f.dFinish, f.dNow])
      .filter((v) => v !== null && v !== undefined).map(SX);
    const x = d3.scaleLinear().domain([Math.min(...vals, SX(-36000)), Math.max(...vals, SX(36000))])
      .nice().range([mg.l, w - mg.r]);
    const y = d3.scalePoint().domain(F.map((f) => f.e.id)).range([mg.t, h - mg.b]).padding(0.5);
    const g = svg.append("g");
    // Ticks outward from launch, each kept only if its label clears the last one kept.
    const keep = new Set([0]);
    for (const side of [1, -1]) {
      let last = x(SX(0));
      for (const [s] of TICKS.filter(([s]) => Math.sign(s) === side).sort((a, b) => Math.abs(a[0]) - Math.abs(b[0]))) {
        if (Math.abs(x(SX(s)) - last) >= 62) { keep.add(s); last = x(SX(s)); }
      }
    }
    for (const [s, t] of TICKS) {
      if (!keep.has(s)) continue;
      drv(Math.abs(s) >= 3600 ? Math.abs(s) / 3600 : Math.abs(s) / 60, "Fig. 2 axis tick");
      g.append("path").attr("class", `ln yield ${s === 0 ? "" : "ctx"}`).attr("d", L([[x(SX(s)), mg.t - 8], [x(SX(s)), h - mg.b + 4]]));
      g.append("text").attr("class", "tick").attr("x", x(SX(s))).attr("y", h - mg.b + 20).attr("text-anchor", "middle").text(t);
    }
    for (const f of F) {
      const yy = y(f.e.id);
      const lane = g.append("g");
      // Long names take two lines: before a parenthesis, else at the space nearest the middle.
      const lines = rowLines(f.e.title);
      const t = phone
        ? lane.append("text").attr("class", "lbl row").attr("x", mg.l).attr("y", yy - 14 - 13 * (lines.length - 1))
        : lane.append("text").attr("class", "lbl row").attr("x", mg.l - 12).attr("y", yy + 4 - 6.5 * (lines.length - 1)).attr("text-anchor", "end");
      lines.forEach((s, i) => t.append("tspan").attr("x", t.attr("x")).attr("dy", i ? 13 : 0).text(s));
      if (!f.L) continue;
      const end = f.dFinish ?? f.dNow;
      lane.append("path").attr("class", `ln model${f.dFinish === null ? " dash" : ""}`)
        .attr("d", L([[x(0), yy], [x(SX(end)), yy]]));
      for (const p of f.pre) {
        lane.append("rect").attr("class", "open").attr("x", x(SX(p.dCommit)) - 4).attr("y", yy - 4).attr("width", 8).attr("height", 8);
        lane.append("rect").attr("class", "dot").attr("x", x(SX(p.dPush)) - 3.5).attr("y", yy - 3.5).attr("width", 7).attr("height", 7);
      }
      if (f.dSummary !== null) {
        const sx = x(SX(f.dSummary));
        lane.append("path").attr("class", "dot").attr("d", `M${sx},${yy - 6}L${sx + 6},${yy}L${sx},${yy + 6}L${sx - 6},${yy}Z`);
      }
      lane.append("rect").attr("class", "hit").attr("x", mg.l).attr("y", yy - 14).attr("width", w - mg.l - mg.r).attr("height", 28)
        .on("pointermove pointerdown", (e) => {
          const [px] = d3.pointer(e, svg.node());
          const p = f.pre[0];
          tip.show(px, yy, `${f.e.title}: plan committed ${dayMonth(p.committed)} ${clock(p.committed)}, ` +
            `pushed ${clock(p.pushed.time)} (${duration(drv(Math.round(Math.abs(p.dPush)), "s, push vs launch"))} ` +
            `${p.dPush < 0 ? "before" : "after"} launch); launched ${clock(f.e.launched.time)}` +
            `${f.e.finished ? `; finished ${clock(f.e.finished.time)}` : "; still running"}.`);
        })
        .on("pointerleave", leave(tip));
    }
    svg.on("pointerdown.away", (e) => { if (!e.target.closest(".hit")) tip.hide(); });
  });
}

/* ---------------- Fig. 4: the headline posterior ---------------- */

export function fig4(fig, D) {
  const hd = D.head, hist = hd.hist, i95 = hd.i95;
  const tip = probe(fig);
  const bins = hist.edges.slice(0, -1).map((e, k) => ({
    x0: hist.centre + e, x1: hist.centre + hist.edges[k + 1], n: hist.counts[k],
  }));
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: 20, r: 20, t: 40, b: 38 };
    const x = d3.scaleLinear().domain([bins[0].x0, bins.at(-1).x1]).range([mg.l, w - mg.r]);
    const y = d3.scaleLinear().domain([0, d3.max(bins, (b) => b.n)]).range([h - mg.b, mg.t + 10]);
    const g = svg.append("g");
    g.append("path").attr("class", "ln ctx").attr("d", L([[mg.l, h - mg.b], [w - mg.r, h - mg.b]]));
    xAxis(g, x, h - mg.b, x.ticks(phone ? 4 : 8), (t) => `${t}°`);
    const step = bins.flatMap((b) => [[x(b.x0), y(b.n)], [x(b.x1), y(b.n)]]);
    g.append("path").attr("class", "ln model").attr("d", L([[x(bins[0].x0), y(0)], ...step, [x(bins.at(-1).x1), y(0)]]));
    const vline = (v, cls, top = mg.t + 10) => g.append("path").attr("class", `ln ${cls}`).attr("d", L([[x(v), top], [x(v), h - mg.b]]));
    vline(hd.median, "");
    vline(hd.lo68, "sec");
    vline(hd.hi68, "sec");
    vline(i95.lo, "sec dash");
    vline(i95.hi, "sec dash");
    point(g, { id: "N", x: x(hd.jpl), y: mg.t - 6, letter: "N", dx: 0, dy: -14, anchor: "middle", r: 3, title: "Neptune (JPL)" });
    g.append("path").attr("class", "ln").attr("d", L([[x(hd.jpl), mg.t - 6], [x(hd.jpl), mg.t + 8]]));
    point(g, { id: "L", x: x(hd.leVerrier), y: mg.t - 6, letter: "L", mark: "open", dx: 0, dy: -14, anchor: "middle", r: 2.6, title: "Le Verrier's planet" });
    g.append("text").attr("class", "lbl").attr("x", x(hd.lo68) - 6).attr("y", mg.t + 26).attr("text-anchor", "end").text("68 %");
    g.append("text").attr("class", "lbl").attr("x", x(i95.lo) - 6).attr("y", mg.t + 26).attr("text-anchor", "end").text("95 %");
    svg.on("pointermove pointerdown", (e) => {
      const [px] = d3.pointer(e);
      const b = bins.find((bb) => x(bb.x0) <= px && px < x(bb.x1));
      if (!b) return tip.hide();
      tip.show(px, y(b.n), `${fmt(b.x0, 2)}–${fmt(b.x1, 2)}°: ${b.n} of ${D.conj.n_samples} samples`);
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 5: the window ---------------- */

export function fig5(fig, D) {
  const W = D.window, s = W.series, m = W.marks;
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: 50, r: phone ? 24 : 40, t: 30, b: 38 };
    const x = d3.scaleLinear().domain(d3.extent(s.year)).range([mg.l, w - mg.r]);
    const all = [...s.lo95_rel_deg, ...s.hi95_rel_deg, ...s.le_verrier_minus_median_deg];
    const y = d3.scaleLinear().domain(d3.extent(all)).nice(5).range([h - mg.b, mg.t]);
    const g = svg.append("g");
    for (const t of y.ticks(5)) {
      drv(t, "axis tick");
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[mg.l, y(t)], [w - mg.r, y(t)]]));
      g.append("text").attr("class", "tick").attr("x", mg.l - 8).attr("y", y(t) + 4).attr("text-anchor", "end").text(`${minus(t)}°`);
    }
    xAxis(g, x, h - mg.b, x.ticks(phone ? 4 : 8));
    const line = (key, cls) => g.append("path").attr("class", `ln ${cls}`)
      .attr("d", d3.line().x((_, k) => x(s.year[k])).y((v) => y(v))(s[key]));
    line("lo68_rel_deg", "sec");
    line("hi68_rel_deg", "sec");
    line("lo95_rel_deg", "sec dash");
    line("hi95_rel_deg", "sec dash");
    line("jpl_minus_median_deg", "");
    line("le_verrier_minus_median_deg", "model");
    const end = s.year.length - 1;
    point(g, { id: "L", x: x(s.year[end]), y: y(s.le_verrier_minus_median_deg[end]), letter: "L", mark: "none", dx: 6, dy: 0 });
    point(g, { id: "N", x: x(s.year[end]), y: y(s.jpl_minus_median_deg[end]), letter: "N", mark: "none", dx: 6, dy: 0 });
    for (const [k, letter] of [["conjunction", "C"], ["galle_night", "G"]]) {
      const xm = x(m[k].year);
      g.append("path").attr("class", "ln sec dash yield").attr("d", L([[xm, mg.t], [xm, h - mg.b]]));
      point(g, { id: letter, x: xm, y: mg.t - 2, letter, mark: "none", dx: 0, dy: -10, anchor: "middle" });
    }
    svg.on("pointermove pointerdown", (e) => {
      const [px] = d3.pointer(e);
      const k = d3.bisector((v) => v).center(s.year, x.invert(px));
      tip.show(x(s.year[k]), y(s.le_verrier_minus_median_deg[k]), `${fmt(s.year[k], 2)}: Le Verrier ` +
        `${minus(fmt(s.le_verrier_minus_median_deg[k], 2))}°, Neptune ${minus(fmt(s.jpl_minus_median_deg[k], 2))}°, ` +
        `68 % ${minus(fmt(s.lo68_rel_deg[k], 2))} to ${minus(fmt(s.hi68_rel_deg[k], 2))}°, relative to the median`);
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 6: robustness ---------------- */

export function fig6(fig, D) {
  const k = D.checks, names = SENS(k.fits);
  const g0 = (key) => k.fits[key].galle_lon_ecl_of_date_deg;
  const rows = [
    { key: "first", label: "Headline, first published", iv: k.first, dash95: false },
    { key: "reported", label: "Headline, reported", iv: k.reported, dash95: k.reported.i95.state !== "final" },
    ...Object.keys(names).map((key) => {
      const v = key === "s1_e05" ? k.verdicts.s1_e05 : null;
      const tag = v?.pending ? " (verdict pending)" : v?.unresolved ? " (unresolved)" : "";
      // A pending or unresolved verdict is drawn as provisional: dashed ranges and an open dot.
      return { key, label: names[key] + tag, iv: g0(key), dash95: false, provisional: !!(v?.pending || v?.unresolved) };
    }),
  ];
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: phone ? 14 : 340, r: 20, t: phone ? 60 : 56, b: 38 };
    const lo = d3.min(rows, (r) => r.iv.lo95), hi = d3.max(rows, (r) => r.iv.hi95);
    const x = d3.scaleLinear().domain([lo, hi]).nice().range([mg.l, w - mg.r]);
    const y = d3.scalePoint().domain(rows.map((r) => r.key)).range([mg.t, h - mg.b]).padding(0.6);
    const g = svg.append("g");
    xAxis(g, x, h - mg.b, x.ticks(phone ? 4 : 7), (t) => `${t}°`);
    g.append("path").attr("class", "ln yield").attr("d", L([[x(k.jpl), mg.t - 14], [x(k.jpl), h - mg.b]]));
    point(g, { id: "N", x: x(k.jpl), y: mg.t - 14, letter: "N", mark: "none", dx: 0, dy: -10, anchor: "middle" });
    for (const r of rows) {
      const yy = y(r.key);
      if (phone) g.append("text").attr("class", "lbl row").attr("x", mg.l).attr("y", yy - 13).text(r.label);
      else g.append("text").attr("class", "lbl row").attr("x", mg.l - 14).attr("y", yy + 4).attr("text-anchor", "end").text(r.label);
      g.append("path").attr("class", `ln sec${r.dash95 || r.provisional ? " dash" : ""}`).attr("d", L([[x(r.iv.lo95), yy], [x(r.iv.hi95), yy]]));
      g.append("path").attr("class", `ln model${r.key === "first" ? " sec" : ""}${r.provisional ? " dash" : ""}`)
        .attr("d", L([[x(r.iv.lo68), yy], [x(r.iv.hi68), yy]])).style("stroke-width", 3);
      g.append("circle").attr("class", r.provisional ? "open" : "dot").attr("cx", x(r.iv.median)).attr("cy", yy).attr("r", 4);
      g.append("rect").attr("class", "hit").attr("x", mg.l).attr("y", yy - 14).attr("width", w - mg.l - mg.r).attr("height", 28)
        .on("pointermove pointerdown", () => tip.show(x(r.iv.median), yy, `${r.label}: median ${fmt(r.iv.median, 2)}°, 68 % ` +
          `${fmt(r.iv.lo68, 2)}–${fmt(r.iv.hi68, 2)}°, 95 % ${fmt(r.iv.lo95, 2)}–${fmt(r.iv.hi95, 2)}°`))
        .on("pointerleave", leave(tip));
    }
    svg.on("pointerdown.away", (e) => { if (!e.target.closest(".hit")) tip.hide(); });
  });
}

/* ---------------- Fig. 7: the profile likelihood in a (M9) ---------------- */

/** One panel of a profile in a: the gridlines, the reading rule's threshold, the x axis, the
 * thin context curve (B) under the thick one (A), the dots, a title, and N and L as open marks
 * on A when `letters` gives their ids. A `floor` clips the y axis: a point below it is drawn
 * open at the bottom edge (Fig. 9). */
function profilePanel(g, { x, y, A, B, S, left, right, top, bottom, title, letters, yLabels,
                           floor = -Infinity, sparse }) {
  const grid = S.grid_au;
  const pts = (ds) => ds.points.filter((p) => !p.reference && p.dlnl != null);
  const yv = (p) => y(Math.max(p.dlnl, floor));
  for (const t of y.ticks(5)) {
    drv(t, "axis tick");
    g.append("path").attr("class", "ln ctx yield").attr("d", L([[left, y(t)], [right, y(t)]]));
    if (yLabels) {
      g.append("text").attr("class", "tick").attr("x", left - 8).attr("y", y(t) + 4).attr("text-anchor", "end").text(minus(t));
    }
  }
  g.append("path").attr("class", "ln sec yield").attr("d", L([[left, y(S.reading.dlnl)], [right, y(S.reading.dlnl)]]));
  xAxis(g, x, bottom, grid.filter((v) => Number.isInteger(v) && (!sparse || v % 4 === 0)));
  // Each run of grid points is one line; a dip between two families of solutions (the 26 AU
  // point of the all-observations profiles) is left unjoined, so the low-distance branch at
  // 24 AU and the main branch from 28 AU on are separate segments.
  const curve = (ds, cls) => {
    const P = pts(ds);
    for (const seg of profileSegments(P)) {
      if (seg.length > 1) g.append("path").attr("class", `ln ${cls}`).attr("d", L(seg.map((p) => [x(p.a_au), yv(p)])));
    }
    g.selectAll(null).data(P).join("circle").attr("class", (p) => (p.dlnl < floor ? "open" : "dot"))
      .attr("cx", (p) => x(p.a_au)).attr("cy", yv).attr("r", 2.6);
  };
  if (B) curve(B, "sec");
  curve(A, "model");
  if (title) g.append("text").attr("class", "lbl row").attr("x", left).attr("y", top - 12).text(title);
  if (letters) {
    for (const [name, letter, dy] of [["neptune", "N", -13], ["le_verrier", "L", 13]]) {
      const r = A.reference_points[name];
      point(g, { id: letters[letter], x: x(r.a_au), y: y(r.dlnl), letter, mark: "open", dx: 9, dy });
    }
  }
}

/** The panel under the pointer and the nearest grid value in it. */
function panelAt(panels, grid, px, py) {
  const P = panels.find((q) => px >= q.left - 20 && px <= q.right + 20 && py >= q.top - 20 && py <= q.bottom + 20)
    ?? panels[0];
  const a = grid[d3.bisector((v) => v).center(grid, P.x.invert(px))];
  return { P, a };
}

/* ---------------- Fig. 7: the profile likelihood in a under both error models ---------------- */

export function fig7(fig, D) {
  const M = D.m9, T = D.m10, S = M.settings, grid = S.grid_au;
  const spec = [
    { A: M.datasets.all, B: M.datasets.post1781, title: "Gaussian errors", letters: { N: "N", L: "L" } },
    { A: T.datasets.all, B: T.datasets.post1781, title: "Student-t errors", letters: { N: "Nt", L: "Lt" } },
  ];
  const pts = (ds) => ds.points.filter((p) => !p.reference && p.dlnl != null);
  const at = (ds, a) => ds.points.find((p) => !p.reference && p.a_au === a);
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const stacked = w < 600;
    const mg = { l: 54, r: stacked ? 16 : 28, t: 34, b: 38 }, gap = stacked ? 68 : 36;
    const low = d3.min(spec.flatMap((q) => [...pts(q.A), ...pts(q.B)]), (p) => p.dlnl);
    const y0 = Math.min(low, S.reading.dlnl) - 0.5;
    const g = svg.append("g");
    const panels = spec.map((q, i) => {
      const box = stacked
        ? { left: mg.l, right: w - mg.r, top: mg.t + i * ((h - mg.t - mg.b - gap) / 2 + gap),
            bottom: mg.t + i * ((h - mg.t - mg.b - gap) / 2 + gap) + (h - mg.t - mg.b - gap) / 2 }
        : { left: mg.l + i * ((w - mg.l - mg.r - gap) / 2 + gap),
            right: mg.l + i * ((w - mg.l - mg.r - gap) / 2 + gap) + (w - mg.l - mg.r - gap) / 2,
            top: mg.t, bottom: h - mg.b };
      const x = d3.scaleLinear().domain([grid[0] - 1, grid.at(-1) + 1]).range([box.left, box.right]);
      const y = d3.scaleLinear().domain([y0, 0.8]).nice(5).range([box.bottom, box.top]);
      profilePanel(g, { ...q, ...box, x, y, S, yLabels: stacked || i === 0, sparse: box.right - box.left < 400 });
      return { ...q, ...box, x, y };
    });
    svg.on("pointermove pointerdown", (e) => {
      const [px, py] = d3.pointer(e);
      const { P, a } = panelAt(panels, grid, px, py);
      const pa = at(P.A, a), pb = at(P.B, a);
      tip.show(P.x(a), P.y(pa.dlnl), `${P.title}, ${a} AU: all observations Δ ln L ${minus(fmt(pa.dlnl, 2))}, ` +
        `mass ${fmt(pa.mass_earth)} M⊕, longitude ${fmt(pa.galle_lon_deg, 1)}°; from 1781 on ` +
        `${minus(fmt(pb.dlnl, 2))}, ${fmt(pb.mass_earth)} M⊕, ${fmt(pb.galle_lon_deg, 1)}°`);
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 8: mass against a in the posterior (after the fact) ---------------- */

export function fig8(fig, D) {
  const Cn = D.m9.corner;
  const ae = Cn.a_edges_au, me = Cn.mass_edges_earth, nx = ae.length - 1, ny = me.length - 1;
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: 54, r: phone ? 22 : 36, t: 30, b: 38 };
    const x = d3.scaleLinear().domain([ae[0], ae.at(-1)]).range([mg.l, w - mg.r]);
    const y = d3.scaleLinear().domain([me[0], me.at(-1)]).range([h - mg.b, mg.t]);
    const g = svg.append("g");
    for (const t of y.ticks(5)) {
      drv(t, "axis tick");
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[mg.l, y(t)], [w - mg.r, y(t)]]));
      g.append("text").attr("class", "tick").attr("x", mg.l - 8).attr("y", y(t) + 4).attr("text-anchor", "end").text(String(t));
    }
    xAxis(g, x, h - mg.b, x.ticks(phone ? 4 : 6));
    // the posterior: a seeded subsample as stipple, and the density's 68 % and 95 % contours
    stipple(g, Cn.samples.map(([a, m]) => [x(a), y(m)]), { r: 1.1, cls: "dot sec" });
    const values = new Float64Array(nx * ny);
    for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) values[j * nx + i] = Cn.density[i][j];
    const da = ae[1] - ae[0], dm = me[1] - me[0];
    const path = d3.geoPath(d3.geoTransform({
      point(px, py) { this.stream.point(x(ae[0] + da * px), y(me[0] + dm * py)); },
    }));
    for (const c of d3.contours().size([nx, ny]).thresholds([Cn.levels.p95, Cn.levels.p68])(values)) {
      g.append("path").attr("class", `ln ${c.value === Cn.levels.p68 ? "model" : "sec"}`).attr("d", path(c));
    }
    const mk = Cn.marks;
    point(g, { id: "N", x: x(mk.neptune.a_au), y: y(mk.neptune.mass_earth), letter: "N", mark: "open", dx: 10, dy: -12 });
    point(g, { id: "L", x: x(mk.le_verrier.a_au), y: y(mk.le_verrier.mass_earth), letter: "L", mark: "open", dx: 10, dy: -12 });
    const qa = Cn.a_au, qm = Cn.mass_earth;
    svg.on("pointermove pointerdown", () => tip.show(x(qa.median), y(qm.median),
      `Posterior median a ${fmt(qa.median)} AU (68 %: ${fmt(qa.lo68)}–${fmt(qa.hi68)}), ` +
      `mass ${fmt(qm.median)} M⊕ (${fmt(qm.lo68)}–${fmt(qm.hi68)})`)).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 9: the Student-t profile without one pre-1781 observer ---------------- */

export function fig9(fig, D) {
  const T = D.m10, S = T.settings, grid = S.grid_au, A = T.datasets.all;
  const FLOOR = drv(-6, "Fig. 9's y-axis floor");
  const sets = T.leave_one_out.map((k) => T.datasets[k]);
  const at = (ds, a) => ds.points.find((p) => !p.reference && p.a_au === a);
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const mg = { l: 46, r: 12, t: 30, b: 34 }, gx = 26, gy = 68;
    const pw = (w - mg.l - mg.r - gx) / 2, ph = (h - mg.t - mg.b - gy) / 2;
    const g = svg.append("g");
    const panels = sets.map((ds, i) => {
      const col = i % 2, row = Math.floor(i / 2);
      const box = { left: mg.l + col * (pw + gx), right: mg.l + col * (pw + gx) + pw,
                    top: mg.t + row * (ph + gy), bottom: mg.t + row * (ph + gy) + ph };
      const x = d3.scaleLinear().domain([grid[0] - 1, grid.at(-1) + 1]).range([box.left, box.right]);
      const y = d3.scaleLinear().domain([FLOOR, 0.8]).range([box.bottom, box.top]);
      const title = `without ${ds.observer}`;
      profilePanel(g, { A: ds, B: A, S, ...box, x, y, title, yLabels: col === 0, floor: FLOOR, sparse: true });
      return { A: ds, B: A, title, ...box, x, y };
    });
    svg.on("pointermove pointerdown", (e) => {
      const [px, py] = d3.pointer(e);
      const { P, a } = panelAt(panels, grid, px, py);
      const pa = at(P.A, a), pb = at(P.B, a);
      tip.show(P.x(a), P.y(Math.max(pa.dlnl, FLOOR)), `${P.title}, ${a} AU: Δ ln L ${minus(fmt(pa.dlnl, 2))}, ` +
        `mass ${fmt(pa.mass_earth)} M⊕, longitude ${fmt(pa.galle_lon_deg, 1)}°; all observations ` +
        `${minus(fmt(pb.dlnl, 2))}, ${fmt(pb.mass_earth)} M⊕, ${fmt(pb.galle_lon_deg, 1)}°`);
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 10: the posterior on Lalande's nights ---------------- */

export function fig10(fig, D) {
  const M = D.m6e, v = M.values, night = M.histogram_night;
  const edges = v["histogram.edges_deg"], counts = v["histogram.counts"];
  const bins = counts.map((n, k) => ({ x0: edges[k], x1: edges[k + 1], n }));
  const q = (model, k) => v[`${model}.${night}.lon.${k}`];
  const lal = v[`lalande.${night}.ecl_lon_deg`]; // the night of the histogram only
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: 20, r: 20, t: 40, b: 78 };
    const x = d3.scaleLinear().domain([bins[0].x0, bins.at(-1).x1]).range([mg.l, w - mg.r]);
    const y = d3.scaleLinear().domain([0, d3.max(bins, (b) => b.n)]).range([h - mg.b, mg.t + 10]);
    const y0 = h - mg.b;
    const g = svg.append("g");
    g.append("path").attr("class", "ln ctx").attr("d", L([[mg.l, y0], [w - mg.r, y0]]));
    const [t0, t1] = x.domain();
    const ticks = d3.range(Math.ceil(t0), Math.floor(t1) + 1).filter((t) => !phone || t % 2 === 0);
    xAxis(g, x, y0, ticks, (t) => `${t}°`);
    const step = bins.flatMap((b) => [[x(b.x0), y(b.n)], [x(b.x1), y(b.n)]]);
    g.append("path").attr("class", "ln model").attr("d", L([[x(bins[0].x0), y(0)], ...step, [x(bins.at(-1).x1), y(0)]]));
    // The two 68 % ranges, each on its own line under the axis.
    [["headline", "ln model", 38, "Gaussian"], ["student_t", "ln", 62, "Student-t"]].forEach(([m, cls, dy, name]) => {
      const a = x(q(m, "lo68")), b = x(q(m, "hi68")), yy = y0 + dy;
      g.append("path").attr("class", cls).attr("d", L([[a, yy], [b, yy]]));
      for (const e of [a, b]) g.append("path").attr("class", "ln").attr("d", L([[e, yy - 4], [e, yy + 4]]));
      g.append("text").attr("class", "lbl").attr("x", b + 10).attr("y", yy).attr("dominant-baseline", "middle").text(name);
    });
    g.append("path").attr("class", "ln").attr("d", L([[x(lal), mg.t - 6], [x(lal), y0]]));
    point(g, { id: "La", x: x(lal), y: mg.t - 6, letter: "La", dx: 0, dy: -14, anchor: "middle", r: 2.6, title: "Lalande's star" });
    svg.on("pointermove pointerdown", (e) => {
      const [px] = d3.pointer(e);
      const b = bins.find((bb) => x(bb.x0) <= px && px < x(bb.x1));
      if (!b) return tip.hide();
      tip.show(px, y(b.n), `${fmt(b.x0, 2)}–${fmt(b.x1, 2)}°: ${b.n} of ${v.n_draws} draws`);
    }).on("pointerleave", leave(tip));
  });
}
