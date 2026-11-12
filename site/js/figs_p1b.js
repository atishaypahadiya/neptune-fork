// Figs. 1 to 4 of the companion page: small plates (white line on the cyanotype ground). Every
// mark is drawn from site/data/p1b_*.json; axis tick values are derived from the data's range
// and registered as such.

import * as d3 from "d3";
import { mount, point, probe } from "./plate.js";
import { fmt, drv } from "./data.js";
import { facts, caseRows, dayOf } from "./copy_p1b.js";

const L = d3.line();
const sgn = (t, n = 0) => `${t < 0 ? "−" : ""}${Math.abs(t).toFixed(n)}`;
const leave = (tip) => (e) => { if (e.pointerType === "mouse") tip.hide(); };
const tick = (g, x, y, text, anchor = "middle") =>
  g.append("text").attr("class", "tick").attr("x", x).attr("y", y).attr("text-anchor", anchor).text(text);

/* ---------------- Fig. 1: the shortfall follows the assumed distance ---------------- */

export function fig1(fig, D) {
  const { A, printedAt, fit } = facts(D);
  const f1 = fit("H1"), f2 = fit("H2"), c68 = A.o_clean_project_fit_68;
  const rows = [
    { id: "P", label: "Printed solution, first hypothesis", v: printedAt("H1") },
    { id: "Q", label: "Printed solution, second hypothesis", v: printedAt("H2") },
    { id: "R", label: "His data, his least squares", v: A.H1_o_real_least_squares },
    { id: "F", label: `Project fit, distance held at ${fmt(A.S5_a_au)} AU`, v: f1.med, lo: f1.lo, hi: f1.hi },
    { id: "G", label: `Project fit, distance held at ${fmt(A.S6_a_au, 2)} AU`, v: f2.med, lo: f2.lo, hi: f2.hi },
    { label: "Noise-free Neptune: his least squares", short: "Noise-free: his least squares", v: A.o_clean_least_squares, open: true },
    { label: "Noise-free Neptune: his elimination", short: "Noise-free: his elimination", v: A.o_clean_his_procedure, open: true },
    { label: "Noise-free Neptune: project model on his places", short: "Noise-free: project model, his places", v: A.o_clean_project_on_his_places, open: true },
    { label: "Noise-free Neptune: project fit", short: "Noise-free: project fit", v: A.o_clean_project_fit, lo: c68[0], hi: c68[1], open: true },
  ];
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: phone ? 16 : 370, r: 28, t: phone ? 58 : 40, b: 44 };
    const all = rows.flatMap((r) => [r.v, r.lo ?? r.v, r.hi ?? r.v]);
    const x = d3.scaleLinear().domain([Math.min(...all) - 1, Math.max(...all, 0) + 1]).nice()
      .range([mg.l, w - mg.r]);
    const y = d3.scalePoint().domain(rows.map((_, i) => i)).range([mg.t, h - mg.b]).padding(0.5);
    const g = svg.append("g");
    for (const t of x.ticks(phone ? 5 : 8)) {
      drv(t, "axis tick");
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[x(t), mg.t - 8], [x(t), h - mg.b + 6]]));
      tick(g, x(t), h - mg.b + 24, `${sgn(t)}°`);
    }
    g.append("path").attr("class", "ln sec yield").attr("d", L([[x(0), mg.t - 14], [x(0), h - mg.b + 6]]));
    point(g, { id: "N", x: x(0), y: mg.t - 14, letter: "N", mark: "none", dx: 8, dy: -6 });
    rows.forEach((r, i) => {
      const cy = y(i);
      if (phone) g.append("text").attr("class", "lbl row").attr("x", mg.l).attr("y", cy - 20).text(r.short ?? r.label);
      else g.append("text").attr("class", "lbl row").attr("x", mg.l - 16).attr("y", cy + 5).attr("text-anchor", "end").text(r.label);
      if (r.lo !== undefined) g.append("path").attr("class", "ln model").attr("d", L([[x(r.lo), cy], [x(r.hi), cy]]));
      // On phones the letter sits past the end of the row's range, on the row's own line.
      const past = phone ? (r.hi !== undefined ? x(r.hi) - x(r.v) : 0) + 10 : 10;
      if (r.id) point(g, { id: r.id, x: x(r.v), y: cy, letter: r.id, mark: "dot", dx: past, dy: phone ? 1 : -11 });
      else g.append("circle").attr("class", "open").attr("cx", x(r.v)).attr("cy", cy).attr("r", 4.5);
    });
    svg.on("pointermove pointerdown", (e) => {
      const [, py] = d3.pointer(e);
      const i = d3.minIndex(rows, (_, k) => Math.abs(y(k) - py));
      const r = rows[i];
      if (Math.abs(y(i) - py) > y.step() / 2) return tip.hide();
      tip.show(x(r.v), y(i), `${r.label}: ${sgn(r.v, 1)}° from Neptune` +
        (r.lo !== undefined ? ` (${sgn(r.lo, 1)}° to ${sgn(r.hi, 1)}°).` : "."));
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 2: his shift against noise ---------------- */

export function fig2(fig, D) {
  const { A } = facts(D);
  const edges = A.hist_edges_deg;
  const panels = [["H1", "First hypothesis"], ["H2", "Second hypothesis"]];
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const gap = phone ? 92 : 48;
    const mg = { l: 20, r: 20, t: 46, b: 40 };
    const pw = phone ? w - mg.l - mg.r : (w - mg.l - mg.r - gap) / 2;
    const ph = phone ? (h - mg.t - mg.b - gap) / 2 : h - mg.t - mg.b;
    const top = Math.max(...panels.flatMap(([H]) => ["a", "b"].flatMap((m) => A[`${H}_${m}_hist_counts`])));
    const g = svg.append("g");
    const hits = [];
    panels.forEach(([H, title], k) => {
      const x0 = mg.l + (phone ? 0 : k * (pw + gap)), y0 = mg.t + (phone ? k * (ph + gap) : 0);
      const x = d3.scaleLinear().domain([edges[0], edges.at(-1)]).range([x0, x0 + pw]);
      const y = d3.scaleLinear().domain([0, top]).range([y0 + ph, y0]);
      g.append("text").attr("class", "lbl row").attr("x", x0).attr("y", y0 - 18).text(title);
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[x0, y0 + ph], [x0 + pw, y0 + ph]]));
      for (const t of [-180, -90, 0, 90, 180]) {
        drv(t, "axis tick");
        g.append("path").attr("class", "ln ctx yield").attr("d", L([[x(t), y0 + ph], [x(t), y0 + ph + 5]]));
        tick(g, x(t), y0 + ph + 22, `${sgn(t)}°`, t === -180 ? "start" : t === 180 ? "end" : "middle");
      }
      const steps = (c) => c.flatMap((n, i) => [[x(edges[i]), y(n)], [x(edges[i + 1]), y(n)]]);
      g.append("path").attr("class", "ln sec").attr("d", L(steps(A[`${H}_b_hist_counts`])));
      g.append("path").attr("class", "ln model").attr("d", L(steps(A[`${H}_a_hist_counts`])));
      const dx = x(A[`${H}_delta_obs`]);
      g.append("path").attr("class", "ln").attr("d", L([[dx, y0 + ph], [dx, y0 + 6]]));
      if (k === 0) point(g, { id: "D", x: dx, y: y0 + 6, letter: "D", mark: "none", dx: 8, dy: 2 });
      hits.push({ H, title, x, x0, y0, pw, ph });
    });
    svg.on("pointermove pointerdown", (e) => {
      const [px, py] = d3.pointer(e);
      const p = hits.find((q) => px >= q.x0 && px <= q.x0 + q.pw && py >= q.y0 - 10 && py <= q.y0 + q.ph + 10);
      if (!p) return tip.hide();
      const i = Math.min(edges.length - 2, Math.max(0, d3.bisectRight(edges, p.x.invert(px)) - 1));
      tip.show(px, py, `${p.title}, ${sgn(edges[i])}° to ${sgn(edges[i + 1])}°: ` +
        `${A[`${p.H}_a_hist_counts`][i]} draws (first noise model), ` +
        `${A[`${p.H}_b_hist_counts`][i]} (second).`);
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 3: the search on the sky ---------------- */

const hm = (hours) => {
  let hh = Math.floor(hours), mm = Math.round((hours - hh) * 60);
  if (mm === 60) { hh += 1; mm = 0; }
  drv(hh, "hours of right ascension"); drv(mm, "minutes of right ascension");
  return `${hh}h ${String(mm).padStart(2, "0")}m`;
};

export function fig3(fig, D) {
  const { W, K } = facts(D);
  const track = D.challis.track, night = D.challis.map_night;
  const nep = track.find((t) => t.date === night);
  const places = [["C", "col325"], ["A", "S4"], ["B", "S5"], ["J", "LVJ"], ["L", "LVA"]]
    .map(([id, k]) => ({ id, ra: W[`place_${k}_ra_h`], dec: W[`place_${k}_dec_deg`] }));
  const half = drv(K.zone_breadth_arcmin.value / 60 / 2, "half the zone's breadth, degrees");
  const d12 = W.neptune_0812_dec_deg, d30 = W.neptune_0730_dec_deg;
  // The zone's north edge, by the record: the planet inside on 12 August, outside on 30 July
  // by the declination edge, or inside on both nights and entered by the right-ascension end.
  const north = { lo: d12, hi: d12 + 2 * half, mid: d12 + half };
  if (!(d30 > d12)) throw new Error("the planet's declinations on 30 July and 12 August");
  const edgeDate = drv(K.sheet_west_edge_ra.h + K.sheet_west_edge_ra.m / 60, "the sheet's edge, hours");
  const edge1800 = drv(edgeDate + (W.neptune_0812_ra_h - W.neptune_0812_ra1800_h),
    "the sheet's edge for the equinox of 1800, carried to the equinox of date");
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: 54, r: 18, t: 22, b: 42 };
    const ras = [...places.map((p) => p.ra), ...track.map((t) => t.ra_h), edgeDate, edge1800];
    const decs = [...places.map((p) => p.dec), ...track.map((t) => t.dec_deg), north.hi, north.lo - 2 * half];
    const dec0 = d3.mean(decs), k = 15 * Math.cos((dec0 * Math.PI) / 180); // degrees per hour
    const pad = 0.12;
    const spanX = (Math.max(...ras) - Math.min(...ras)) * k + 2 * pad;
    const spanY = Math.max(...decs) - Math.min(...decs) + 2 * pad;
    const s = Math.min((w - mg.l - mg.r) / spanX, (h - mg.t - mg.b) / spanY); // px per degree
    const cx = (mg.l + w - mg.r) / 2, cy = (mg.t + h - mg.b) / 2;
    const ra0 = (Math.max(...ras) + Math.min(...ras)) / 2, dc0 = (Math.max(...decs) + Math.min(...decs)) / 2;
    const X = (ra) => cx - (ra - ra0) * k * s, Y = (dec) => cy - (dec - dc0) * s;
    const xl = mg.l, xr = w - mg.r, yb = h - mg.b;
    const g = svg.append("g");
    g.append("path").attr("class", "ln ctx yield").attr("d", L([[xl, yb], [xr, yb]]));
    g.append("path").attr("class", "ln ctx yield").attr("d", L([[xl, mg.t], [xl, yb]]));
    const raAt = (px) => ra0 - (px - cx) / (k * s), decAt = (py) => dc0 - (py - cy) / s;
    const stepM = phone ? 10 : 4; // minutes of right ascension between ticks
    for (let m = Math.ceil((raAt(xr) * 60) / stepM) * stepM; m <= raAt(xl) * 60; m += stepM) {
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[X(m / 60), yb], [X(m / 60), yb + 5]]));
      tick(g, X(m / 60), yb + 22, hm(m / 60));
    }
    for (let d = Math.ceil(decAt(yb) * 2) / 2; d <= decAt(mg.t); d += 0.5) {
      const deg = Math.trunc(d), min = Math.round(Math.abs(d - deg) * 60);
      drv(Math.abs(deg), "degrees of declination"); drv(min, "minutes of declination");
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[xl - 5, Y(d)], [xl, Y(d)]]));
      tick(g, xl - 9, Y(d) + 4, `−${Math.abs(deg)}° ${String(min).padStart(2, "0")}′`, "end");
    }
    // the zone at the middle of the settings the record allows, and the range of settings
    for (const d of [north.mid, north.mid - 2 * half]) {
      g.append("path").attr("class", "ln sec yield").attr("d", L([[xl, Y(d)], [xr, Y(d)]]));
    }
    for (const d of [north.hi, north.lo - 2 * half]) {
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[xl, Y(d)], [xr, Y(d)]]));
    }
    g.append("text").attr("class", "lbl row").attr("x", xr - 8).attr("y", Y(north.hi) - 8)
      .attr("text-anchor", "end").text("zone");
    // the sheet's western edge, under each reading of its limits
    for (const [ra, label, side] of [[edge1800, String(K.sheet_west_edge_ra.equinox), -1],
                                     [edgeDate, "of date", 1]]) {
      g.append("path").attr("class", "ln").attr("d", L([[X(ra), mg.t], [X(ra), yb]]));
      g.append("text").attr("class", "lbl row").attr("x", X(ra) + 7 * side).attr("y", mg.t + 14)
        .attr("text-anchor", side < 0 ? "end" : "start").text(label);
    }
    g.append("path").attr("class", "ln orbit").attr("d", L(track.map((t) => [X(t.ra_h), Y(t.dec_deg)])));
    point(g, { id: "N", x: X(nep.ra_h), y: Y(nep.dec_deg), letter: "N", mark: "dot", dx: 10, dy: 12 });
    for (const p of places) {
      point(g, { id: p.id, x: X(p.ra), y: Y(p.dec), letter: p.id, mark: "open", dx: 10, dy: -10 });
    }
    const pts = [{ id: "N", ra: nep.ra_h, dec: nep.dec_deg }, ...places];
    svg.on("pointermove pointerdown", (e) => {
      const [px, py] = d3.pointer(e);
      const p = d3.least(pts, (q) => Math.hypot(X(q.ra) - px, Y(q.dec) - py));
      if (Math.hypot(X(p.ra) - px, Y(p.dec) - py) > 40) return tip.hide();
      tip.show(X(p.ra), Y(p.dec), `${p.id}: right ascension ${hm(p.ra)}, declination ` +
        `${sgn(p.dec, 2)}°, ${dayOf(night)}.`);
    }).on("pointerleave", leave(tip));
  });
}

/* ---------------- Fig. 4: the cases against the reading rule ---------------- */

export function fig4(fig, D) {
  const { K } = facts(D);
  const rows = caseRows(D);
  const tip = probe(fig);
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const mg = { l: phone ? 16 : 430, r: 26, t: phone ? 46 : 26, b: 42 };
    const x = d3.scaleLinear().domain([0, 1]).range([mg.l, w - mg.r]);
    const y = d3.scalePoint().domain(rows.map((r) => r.key)).range([mg.t, h - mg.b]).padding(0.5);
    const g = svg.append("g");
    for (const t of [0, 0.5, 1]) {
      drv(t, "axis tick");
      g.append("path").attr("class", "ln ctx yield").attr("d", L([[x(t), mg.t - 6], [x(t), h - mg.b + 6]]));
      tick(g, x(t), h - mg.b + 24, String(t));
    }
    for (const t of [K.rule_thresholds.low, K.rule_thresholds.high]) {
      g.append("path").attr("class", "ln sec yield").attr("d", L([[x(t), mg.t - 6], [x(t), h - mg.b + 6]]));
      tick(g, x(t), h - mg.b + 24, String(t));
    }
    for (const r of rows) {
      const cy = y(r.key), label = r.dagger ? `${r.name} †` : r.name;
      if (phone) g.append("text").attr("class", "lbl row").attr("x", mg.l).attr("y", cy - 14).text(r.dagger ? `${r.short} †` : r.short);
      else g.append("text").attr("class", "lbl row").attr("x", mg.l - 16).attr("y", cy + 5).attr("text-anchor", "end").text(label);
      if (r.hi > 0) g.append("path").attr("class", "ln model").attr("d", L([[x(r.lo), cy], [x(r.hi), cy]]));
      if (r.dagger) g.append("path").attr("class", "ln").attr("d", L([[x(r.hlo), cy + 7], [x(r.hhi), cy + 7]]));
      g.append("circle").attr("class", "dot").attr("cx", x(r.mean)).attr("cy", cy).attr("r", 4);
    }
    svg.on("pointermove pointerdown", (e) => {
      const [, py] = d3.pointer(e);
      const r = d3.least(rows, (q) => Math.abs(y(q.key) - py));
      if (Math.abs(y(r.key) - py) > y.step() / 2) return tip.hide();
      tip.show(x(r.mean), y(r.key), `${r.name}: ${r.p}. ${r.zeroWhy ? `${r.zeroWhy}. ` : ""}` +
        `Conditional on ${r.cond}.`);
    }).on("pointerleave", leave(tip));
  });
}
