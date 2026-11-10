// Plates I, II and III. All positions are heliocentric, in AU, from site/data/ (m7_plates.json,
// m7_motion.json); they are drawn to scale with y (ecliptic north view) pointing up and the
// equinox to the right. Stipples, orbits and points are data; the Sun's stipple is texture.

import * as d3 from "d3";
import { mount, point, stipple, sunTexture, orbit } from "./plate.js";
import { fmt, drv } from "./data.js";
import { marks, textQuads, crossings, placeLetters, yieldLines } from "./clear.js";

const TAU = 2 * Math.PI;
const lonOf = ([x, y]) => (Math.atan2(y, x) + TAU) % TAU;
const rOf = ([x, y]) => Math.hypot(x, y);

/** The last stretch of a time-ordered path that covers one full turn (or all of it). */
function lastTurn(path) {
  let total = 0;
  for (let i = path.length - 1; i > 0; i--) {
    let d = lonOf(path[i]) - lonOf(path[i - 1]);
    if (d > Math.PI) d -= TAU;
    if (d < -Math.PI) d += TAU;
    total += Math.abs(d);
    if (total >= TAU) return path.slice(i - 1);
  }
  return path;
}

/** Points along a circle of radius r (px) about (cx, cy), for a textPath, reading left to right. */
function labelArc(cx, cy, r, lonMid, spanRad) {
  const up = Math.sin(lonMid) >= 0; // upper half: run toward decreasing longitude
  const a0 = up ? lonMid + spanRad / 2 : lonMid - spanRad / 2;
  const a1 = up ? lonMid - spanRad / 2 : lonMid + spanRad / 2;
  const pts = d3.range(0, 121).map((k) => {
    const a = a0 + ((a1 - a0) * k) / 120;
    return [cx + r * Math.cos(a), cy - r * Math.sin(a)];
  });
  return { d: d3.line()(pts), up };
}

const ctx2d = document.createElement("canvas").getContext("2d");
/** A text's drawn width in the font and letter-spacing of `el`. Measured here because WebKit
 * leaves letter-spacing out of getComputedTextLength, which clipped labels on iPhones. */
function textWidth(el, text) {
  const cs = getComputedStyle(el);
  ctx2d.font = `${cs.fontStyle} ${cs.fontWeight} ${cs.fontSize} ${cs.fontFamily}`;
  return ctx2d.measureText(text).width + (parseFloat(cs.letterSpacing) || 0) * text.length;
}

let arcId = 0;
/**
 * Text along a circle about (cx, cy), at the first of `spots` where it fits its path and
 * crosses nothing. A spot is { r, lon } (the text's middle) or { r, from } (its start, running
 * toward increasing longitude). Shrinks once if it is too long for `maxSpan`; left out if no
 * spot is clear (the caption carries the same words).
 */
function arcText(g, cx, cy, text, spots, { maxSpan = 1.9 * Math.PI } = {}) {
  const svg = g.node().ownerSVGElement;
  const probe = g.append("text").attr("class", "lbl").text(text);
  const len0 = Math.max(textWidth(probe.node(), text), probe.node().getComputedTextLength());
  const size0 = parseFloat(getComputedStyle(probe.node()).fontSize);
  probe.remove();
  for (const spot of spots) {
    let size = null, span = (len0 * 1.08) / spot.r; // 8 % to spare at the ends
    if (span > maxSpan) { size = 11; span = (span * size) / size0; }
    if (span > maxSpan) continue;
    const lonMid = spot.from !== undefined ? spot.from + span / 2 : spot.lon;
    const id = `arc-${arcId++}`;
    const { d } = labelArc(cx, cy, up0(lonMid) ? spot.r : spot.r + 15, lonMid, span);
    const path = g.append("path").attr("id", id).attr("d", d).attr("fill", "none");
    const t = g.append("text").attr("class", "lbl");
    if (size) t.style("font-size", `${size}px`);
    t.append("textPath").attr("href", `#${id}`).attr("startOffset", "50%")
      .attr("text-anchor", "middle").text(text);
    const fits = textWidth(t.node(), text) <= path.node().getTotalLength() - 2;
    if (fits && crossings(textQuads(svg, t.node()), marks(svg, { skip: new Set([t.node()]) })) === 0) {
      return t;
    }
    t.remove();
    path.remove();
  }
  return null;
}
const up0 = (lon) => Math.sin(lon) >= 0;

/* ---------------- Plate I: Galle's night ---------------- */

export function plateI(fig, D) {
  const P = D.plates.plate_i_galle_night, b = P.bodies;
  const tb = fig.querySelector(".titleblock");
  mount(fig, (svg, w, h) => {
    const phone = w < 600, wide = w >= 900;
    const sub = P.posterior.subsample_xy;
    const rMax = Math.max(...sub.map(rOf), rOf(b.le_verrier.xy)) * 1.06;
    let s, cx, cy;
    if (wide) {
      // Title bottom-left (CSS); the drawing takes the right of the plate at full height.
      s = Math.min(h * 0.43, w * 0.3) / rMax;
      cx = w * 0.66; cy = h * 0.48;
    } else {
      const avail = h - tb.offsetHeight - (phone ? 40 : 64) - 40;
      s = Math.min(avail / 2, w * (phone ? 0.46 : 0.44)) / rMax;
      cx = w / 2; cy = 40 + avail / 2;
    }
    const X = ([x, y]) => [cx + s * x, cy - s * y];
    const g = svg.append("g");

    orbit(g, lastTurn(P.orbits.uranus_jpl).map(X));
    orbit(g, P.orbits.neptune_jpl.map(X));
    orbit(g, P.orbits.le_verrier.map(X), { closed: true });
    stipple(g.append("g"), sub.map(X), { r: phone ? 0.65 : 0.8 });
    stipple(g.append("g").attr("class", "ctx"), sunTexture(cx, cy, 6.5), { r: 0.8 });

    point(g, { id: "S", x: cx, y: cy, letter: "S", dx: -12, dy: -12, anchor: "end", mark: "none", title: "Sun" });
    const [ux, uy] = X(b.uranus_fit.xy);
    point(g, { id: "U", x: ux, y: uy, letter: "U", title: "Uranus" });
    const [nx, ny] = X(b.neptune_jpl.xy);
    point(g, { id: "N", x: nx, y: ny, letter: "N", dx: -12, dy: -2, anchor: "end", title: "Neptune (JPL)" });
    const [lx, ly] = X(b.le_verrier.xy);
    point(g, { id: "L", x: lx, y: ly, letter: "L", mark: "open", dx: 9, dy: 12,
               anchor: "start", title: "Le Verrier's planet" });

    // Labels along orbits, drawn after the points so they keep clear of them, and only where
    // they fit and cross nothing; the caption carries the same words. The planet's starts past
    // the stipple's 95 % edge. It carries no numbers: the plate's longitudes are heliocentric,
    // the headline's geocentric.
    const q = P.posterior.lon_deg, rad = (deg) => (deg * Math.PI) / 180;
    const rLab = s * P.posterior.r_au.hi95 + 8;
    arcText(g, cx, cy, "the planet",
      [rLab, rLab + 10].flatMap((r) => [4, 8, 12, 18, 24, 32].map((k) => ({ r, from: rad(q.hi95 + k) }))),
      { maxSpan: 1.2 * Math.PI });
    // Uranus's path: Table 1's sources and the span of their years, in the lower half so it
    // reads upright.
    const names = [...new Set(D.table1.rows.map((r) => r.source.replace(/ series$/, "")))];
    const ty = D.table1.total.years;
    const rU = s * rOf(b.uranus_fit.xy) + 9;
    arcText(g, cx, cy, `${names.join(" · ")}, ${ty[0]}–${ty[1]}`,
      [235, 245, 225, 255, 215, 265, 205].flatMap((a) => [{ r: rU, lon: rad(a) }, { r: rU + 8, lon: rad(a) }]),
      { maxSpan: 0.9 * Math.PI });
  }, { mark: true });
}

/* ---------------- Plate II: the conjunction ---------------- */

export function plateII(fig, D) {
  const P = D.plates.plate_ii_conjunction, b = P.bodies, s0 = D.conj.conjunction_helio;
  const K = P.angle_magnification;
  const med = (s0.all.median * Math.PI) / 180;
  const wrap = (a) => Math.atan2(Math.sin(a), Math.cos(a));
  // Magnify the angle about the line SC by K; keep the radius. Points more than 80° (drawn)
  // from SC are left out of the orbit paths.
  const mag = (p) => {
    const d = wrap(lonOf(p) - med) * K;
    return Math.abs(d) > (80 * Math.PI) / 180 ? null : [med + d, rOf(p)];
  };
  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const rLv = rOf(b.le_verrier.xy);
    const rMax = Math.max(rLv, P.posterior.r_au.hi95) * 1.04;
    const top = phone ? 40 : 64, bottom = phone ? 52 : 80;
    const s = (h - top - bottom) / rMax;
    const cx = w / 2 - s * rMax * 0.5 * Math.cos(med), cy = top;
    const XP = ([a, r]) => [cx + s * r * Math.cos(a), cy - s * r * Math.sin(a)];
    const X = (p) => XP(mag(p));
    const g = svg.append("g");
    const path = (pts, closed) => {
      const runs = [[]];
      for (const p of pts) {
        const m = mag(p);
        if (m) runs.at(-1).push(XP(m));
        else if (runs.at(-1).length) runs.push([]);
      }
      for (const r of runs) if (r.length > 1) orbit(g, r);
      return closed;
    };
    path(lastTurn(P.orbits.uranus_jpl));
    path(P.orbits.neptune_jpl);
    path(P.orbits.le_verrier);

    // The line SU produced to C, the line SL, and the angle CSL (labelled with its true value).
    const lvA = mag(b.le_verrier.xy)[0];
    const end = rMax;
    // Construction lines: they break where a letter needs the room (clear.js).
    g.append("path").attr("class", "ln yield").attr("d", d3.line()([[cx, cy], XP([med, end])]));
    g.append("path").attr("class", "ln sec yield").attr("d", d3.line()([[cx, cy], XP([lvA, end])]));
    const arc = d3.range(0, 25).map((k) => XP([med + ((lvA - med) * k) / 24, end]));
    g.append("path").attr("class", "ln sec yield").attr("d", d3.line()(arc));
    const [ax, ay] = XP([(med + lvA) / 2, end]);
    g.append("text").attr("class", "lbl").attr("x", ax).attr("y", ay + 20).attr("text-anchor", "middle")
      .text(`${fmt(Math.abs(s0.le_verrier.minus_median_deg), 2)}°`);

    stipple(g.append("g"), P.posterior.subsample_xy.map(X), { r: phone ? 0.65 : 0.8 });
    stipple(g.append("g").attr("class", "ctx"), sunTexture(cx, cy, 6.5), { r: 0.8 });

    point(g, { id: "S", x: cx, y: cy, letter: "S", dx: 12, dy: -4, mark: "none", title: "Sun" });
    const [ux, uy] = X(b.uranus_fit.xy);
    point(g, { id: "U", x: ux, y: uy, letter: "U", dx: -12, dy: 0, anchor: "end", title: "Uranus (fitted)" });
    const [mx, my] = XP([med, P.posterior.r_au.median]);
    const cOff = s * (P.posterior.r_au.hi68 - P.posterior.r_au.median);
    const sep = Math.max(cOff, 12);
    point(g, { id: "C", x: mx, y: my, letter: "C", dx: -(phone ? 26 : 40), dy: -sep, anchor: "end", mark: "none", title: "C: posterior median direction" });
    const [nx, ny] = X(b.neptune_jpl.xy);
    point(g, { id: "N", x: nx, y: ny, letter: "N", dx: -(phone ? 26 : 40), dy: my - ny + sep + 6, anchor: "end", r: 2.6, title: "Neptune (JPL)" });
    const [lx, ly] = X(b.le_verrier.xy);
    point(g, { id: "L", x: lx, y: ly, letter: "L", mark: "open", dx: 12, dy: -6, title: "Le Verrier's planet" });
  });
}

/* ---------------- Plate III: the pull, with Fig. 3 ---------------- */

export function plateIII(fig, fig3, side, D, onFrame) {
  const m = D.motion, F = m.frames, n = F.t_year.length, K = m.offset_exaggeration;
  const jdC = m.conjunction.grid_10_day.jd_tdb;
  const iC = d3.bisector((d) => d).center(F.jd_tdb, jdC);
  const state = { i: 0, update: null, update3: null };

  const lerp = (key) => {
    const j = d3.bisectLeft(F.jd_tdb, jdC);
    const t = (jdC - F.jd_tdb[j - 1]) / (F.jd_tdb[j] - F.jd_tdb[j - 1]);
    return F[key][j - 1] + t * (F[key][j] - F[key][j - 1]);
  };
  const nepC = [lerp("neptune_x"), lerp("neptune_y")];

  mount(fig, (svg, w, h) => {
    const phone = w < 600;
    const rMax = Math.max(...F.neptune_x.map((x, k) => Math.hypot(x, F.neptune_y[k]))) * 1.08;
    const s = (Math.min(w, h) / 2 - (phone ? 18 : 40)) / rMax;
    const cx = w / 2, cy = h / 2;
    const X = (x, y) => [cx + s * x, cy - s * y];
    const g = svg.append("g");
    const U = F.uranus_x.map((x, k) => X(x, F.uranus_y[k]));
    const N = F.neptune_x.map((x, k) => X(x, F.neptune_y[k]));
    orbit(g, U);
    orbit(g, N);
    // K on the line from S, at Neptune's position at the Uranus-Neptune conjunction of the run
    // (not Plate II's C, the conjunction with the fitted planet).
    const [qx, qy] = X(nepC[0] * 1.06, nepC[1] * 1.06);
    g.append("path").attr("class", "ln ctx yield").attr("d", d3.line()([[cx, cy], [qx, qy]]));
    point(g, { id: "K", x: qx, y: qy, letter: "K", dx: 0, dy: 16, anchor: "middle", mark: "none", title: "K: conjunction of Uranus and Neptune" });
    stipple(g.append("g").attr("class", "ctx"), sunTexture(cx, cy, 6.5), { r: 0.8 });
    point(g, { id: "S", x: cx, y: cy, letter: "S", dx: 12, dy: -10, mark: "none", title: "Sun" });

    const trail = g.append("path").attr("class", "ln model");
    const pull = g.append("path").attr("class", "ln").attr("data-leader-of", "pull");
    const pullLbl = g.append("text").attr("class", "lbl").attr("data-leader", "pull");
    const free = g.append("circle").attr("class", "open").attr("r", 4.5);
    const pu = point(g, { id: "U", x: 0, y: 0, letter: "U", r: 3.4, title: "Uranus" });
    const pn = point(g, { id: "N", x: 0, y: 0, letter: "N", r: 3.4, title: "Neptune" });
    pu.select(".glow").style("opacity", 0.6);
    pn.select(".glow").style("opacity", 0.6);

    const svgNode = svg.node(), letters = [pu, pn].map((p) => p.select("text.letter").node());
    state.update = (i) => {
      const [ux, uy] = U[i], [nx, ny] = N[i];
      trail.attr("d", d3.line()(U.slice(0, i + 1)));
      pu.attr("transform", `translate(${ux},${uy})`);
      pn.attr("transform", `translate(${nx},${ny})`);
      pull.attr("d", d3.line()([[ux, uy], [nx, ny]]));
      const [fx, fy] = X(F.uranus_x[i] + K * F.uranus_no_neptune_dx[i], F.uranus_y[i] + K * F.uranus_no_neptune_dy[i]);
      free.attr("cx", fx).attr("cy", fy);
      let ang = (Math.atan2(ny - uy, nx - ux) * 180) / Math.PI;
      if (ang > 90) ang -= 180;
      if (ang < -90) ang += 180;
      const mx = (ux + nx) / 2, my = (uy + ny) / 2;
      pullLbl.attr("transform", `translate(${mx},${my}) rotate(${ang}) translate(0,-7)`)
        .attr("text-anchor", "middle")
        .text(F.accel_along_track[i] >= 0 ? "drawn ahead" : "held back");
      // The letters move with their planets: each frame they take a clear spot (clear.js).
      pullLbl.style("display", "none");
      placeLetters(svgNode, letters);
      // The label only on a line long enough to carry it, and where it crosses nothing but
      // its own line.
      const len = pullLbl.node().getComputedTextLength();
      const room = Math.hypot(nx - ux, ny - uy) - len;
      const inside = mx - len / 2 > 8 && mx + len / 2 < w - 8 && my > 24 && my < h - 24;
      pullLbl.style("display", room > 64 && inside ? null : "none");
      if (room > 64 && inside) {
        const obs = marks(svgNode, { skip: new Set([pullLbl.node()]) });
        if (crossings(textQuads(svgNode, pullLbl.node()), obs, { leader: "pull" })) pullLbl.style("display", "none");
      }
    };
    state.update(state.i);
  });

  // Fig. 3: the longitude difference, with a cursor.
  mount(fig3, (svg, w, h) => {
    const mg = { l: 46, r: 14, t: 16, b: 26 };
    const x = d3.scaleLinear().domain(d3.extent(F.t_year)).range([mg.l, w - mg.r]);
    const y = d3.scaleLinear().domain(d3.extent(F.dlon_arcsec)).nice(4).range([h - mg.b, mg.t]);
    const g = svg.append("g");
    for (const t of y.ticks(4)) {
      drv(t, "Fig. 3 axis tick");
      g.append("path").attr("class", "ln ctx yield").attr("d", d3.line()([[mg.l, y(t)], [w - mg.r, y(t)]]));
      g.append("text").attr("class", "tick").attr("x", mg.l - 6).attr("y", y(t) + 4)
        .attr("text-anchor", "end").text(`${t < 0 ? "−" : ""}${Math.abs(t)}″`);
    }
    for (const t of x.ticks(w < 380 ? 3 : 5)) {
      drv(t, "Fig. 3 axis tick");
      g.append("text").attr("class", "tick").attr("x", x(t)).attr("y", h - 8)
        .attr("text-anchor", "middle").text(t);
    }
    g.append("path").attr("class", "ln model")
      .attr("d", d3.line().x((_, k) => x(F.t_year[k])).y((v) => y(v))(F.dlon_arcsec));
    const xc = x(F.t_year[iC]);
    g.append("path").attr("class", "ln sec dash yield").attr("d", d3.line()([[xc, mg.t], [xc, h - mg.b]]));
    point(g, { id: "K", x: xc, y: mg.t + 4, letter: "K", dx: 6, dy: 2, mark: "none", title: "K: conjunction of Uranus and Neptune" });
    const cur = g.append("path").attr("class", "ln yield");
    const dot = g.append("circle").attr("class", "dot").attr("r", 3.2);
    state.update3 = (i) => {
      const xi = x(F.t_year[i]);
      cur.attr("d", d3.line()([[xi, mg.t], [xi, h - mg.b]])).style("display", null);
      delete cur.node().dataset.d0;
      dot.attr("cx", xi).attr("cy", y(F.dlon_arcsec[i]));
      yieldLines(svg.node()); // the cursor breaks where it passes a label
    };
    state.update3(state.i);
  });

  // Controls: play/pause, scrubber (arrow keys work on the range input), reduced motion.
  const btn = side.querySelector(".play"), scrub = side.querySelector(".scrub");
  scrub.max = String(n - 1);
  const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
  let playing = false, raf = 0, lastT = 0, acc = 0, inView = false;
  const FPS = n / 24; // the whole run in about 24 s

  const show = (i) => {
    state.i = Math.max(0, Math.min(n - 1, i));
    scrub.value = String(state.i);
    state.update?.(state.i);
    state.update3?.(state.i);
    onFrame(state.i);
  };
  const tick = (t) => {
    if (!playing) return;
    acc += lastT ? ((t - lastT) / 1000) * FPS : 0;
    lastT = t;
    if (acc >= 1) {
      const steps = Math.floor(acc);
      acc -= steps;
      show(state.i + steps >= n ? 0 : state.i + steps);
    }
    raf = requestAnimationFrame(tick);
  };
  const setPlaying = (p) => {
    playing = p;
    btn.querySelector(".i-play").hidden = p;
    btn.querySelector(".i-pause").hidden = !p;
    btn.setAttribute("aria-label", p ? "Pause" : "Play");
    cancelAnimationFrame(raf);
    lastT = 0;
    if (p) raf = requestAnimationFrame(tick);
  };
  let userPaused = false;
  btn.addEventListener("click", () => { userPaused = playing; setPlaying(!playing); });
  scrub.addEventListener("input", () => { userPaused = true; setPlaying(false); show(Number(scrub.value)); });
  // Plays while Plate III's section is in view: on phones the controls and Fig. 3 sit a screen
  // below the plate, so watching the plate alone would stop the run as soon as it started.
  new IntersectionObserver(([e]) => {
    inView = e.isIntersecting;
    if (!inView) setPlaying(false);
    else if (!reduce && !userPaused) setPlaying(true);
  }, { threshold: 0.2 }).observe(fig.closest(".split") ?? fig);

  show(reduce ? iC : 0);
  setPlaying(false);
  return { show, iC, setPlaying };
}
