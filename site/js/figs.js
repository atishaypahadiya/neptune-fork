// Figs. 1, 2, 4, 5 and 6: small plates (white line on the cyanotype ground). Every mark is drawn
// from site/data/; axis tick values are derived from the data's range and registered as such.

import * as d3 from "d3";
import { mount, point, probe } from "./plate.js";
import { fmt, drv, clock, dayMonth, duration } from "./data.js";
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
