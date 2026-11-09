// Panel 2, what the data measured: the planet at its 1821 conjunction with Uranus against
// Galle's night, 25 years later (after-the-fact diagnostic, docs/notes/M6.md).
// Geometry (docs/DESIGN.md, "time as layers"): two translucent cards stand on a floor, the
// conjunction in front and Galle's night behind. Each carries the posterior's histograms
// (0.25° bins, each group normalised to unit area, as in docs/notes/M6c_conjunction.png),
// centred on that epoch's median. A ribbon on the floor joins the two 68 % ranges, and
// hairlines join each group's median, so the spreading is a shape. The camera is orthographic,
// so widths on the two cards compare directly. Every value comes from `d` (site/js/data.js).

import * as THREE from "three";
import { fmt, range, signed, whereIn } from "./data.js";
import { PHONE, REDUCED, drawAnnotations, fatLine, hairlines, watchVisibility } from "./common.js";

const SX = 0.68; // world units per degree
const SPAN = 8; // degrees either side of each card's median
const HY = 5.6; // world units per (unit-area density, per degree)
const CARD_H = 4.4;
const Z = { conj: 7, galle: -7 }; // the conjunction in front, Galle's night behind
const GROUPS = ["core", "inner", "outer"]; // "other cells" (a few %) is in the text, not drawn

/** Step outline of one group's histogram, unit area, relative to the card's median. */
function steps(hist, group, z) {
  const counts = hist.counts[group];
  const n = counts.reduce((s, c) => s + c, 0);
  const pts = [];
  hist.edges.slice(0, -1).forEach((lo, i) => {
    const hi = hist.edges[i + 1];
    if (Math.abs(lo) > SPAN || Math.abs(hi) > SPAN) return;
    const y = (counts[i] / n / (hi - lo)) * HY;
    pts.push(new THREE.Vector3(lo * SX, y, z), new THREE.Vector3(hi * SX, y, z));
  });
  const peak = pts.reduce((m, p) => (p.y > m.y ? p : m), pts[0]);
  pts.unshift(new THREE.Vector3(pts[0].x, 0, z));
  pts.push(new THREE.Vector3(pts.at(-1).x, 0, z));
  return { pts, peak, share: n };
}

export function mountConjunction(root, d) {
  const canvas = root.querySelector("canvas");
  const svg = root.querySelector("svg.leaders");
  let renderer;
  try {
    renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
  } catch {
    root.classList.add("no-webgl");
    return;
  }
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 1);
  const scene = new THREE.Scene();
  const camera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 200);
  const fat = [];
  const add = (o) => { scene.add(o); if (o.isLine2) fat.push(o); return o; };
  const rel = (card, lon) => (lon - card.all.median) * SX;
  const W = SPAN * SX + 0.4;

  // Floor: the base slab's outline and the time axis between the cards.
  add(hairlines([
    -W, 0, Z.conj + 2, W, 0, Z.conj + 2, W, 0, Z.conj + 2, W, 0, Z.galle - 2,
    W, 0, Z.galle - 2, -W, 0, Z.galle - 2, -W, 0, Z.galle - 2, -W, 0, Z.conj + 2,
  ], 0.22));

  // The ribbon: the 68 % range at the conjunction, carried to the 68 % range on Galle's night.
  const rib = [rel(d.conj, d.conj.all.lo68), rel(d.conj, d.conj.all.hi68),
               rel(d.galle, d.galle.all.lo68), rel(d.galle, d.galle.all.hi68)];
  const ribbonGeo = new THREE.BufferGeometry();
  const ribPos = new Float32Array(12);
  ribbonGeo.setAttribute("position", new THREE.BufferAttribute(ribPos, 3));
  ribbonGeo.setIndex([0, 1, 2, 1, 3, 2]);
  const ribbon = add(new THREE.Mesh(ribbonGeo, new THREE.MeshBasicMaterial({
    color: 0xffffff, transparent: true, opacity: 0.14, side: THREE.DoubleSide, depthWrite: false,
  })));
  const ribEdges = add(hairlines(new Array(12).fill(0), 0.75));
  const setRibbon = (p) => {
    const zg = Z.conj + (Z.galle - Z.conj) * p;
    const lerp = (a, b) => a + (b - a) * p;
    const gl = lerp(rib[0], rib[2]), gh = lerp(rib[1], rib[3]);
    ribPos.set([rib[0], 0, Z.conj, rib[1], 0, Z.conj, gl, 0, zg, gh, 0, zg]);
    ribbonGeo.attributes.position.needsUpdate = true;
    ribEdges.geometry.attributes.position.array.set(
      [rib[0], 0, Z.conj, gl, 0, zg, rib[1], 0, Z.conj, gh, 0, zg]);
    ribEdges.geometry.attributes.position.needsUpdate = true;
  };

  // Each group's median, joined across the 25 years.
  const medianLinks = [];
  for (const g of GROUPS) {
    medianLinks.push(rel(d.conj, d.conj.groups[g].median), 0, Z.conj,
                     rel(d.galle, d.galle.groups[g].median), 0, Z.galle);
  }
  const links = add(hairlines(medianLinks, 0.45));

  // The two cards, their histograms, Le Verrier (dashed) and Neptune (diamond, in the overlay).
  const cards = {};
  const peaks = {};
  for (const [name, card] of [["conj", d.conj], ["galle", d.galle]]) {
    const z = Z[name];
    const plane = add(new THREE.Mesh(new THREE.PlaneGeometry(2 * W, CARD_H),
      new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.045,
                                    side: THREE.DoubleSide, depthWrite: false })));
    plane.position.set(0, CARD_H / 2, z);
    const frame = add(hairlines([-W, 0, z, -W, CARD_H, z, -W, CARD_H, z, W, CARD_H, z,
                                 W, CARD_H, z, W, 0, z, W, 0, z, -W, 0, z], 0.35));
    const ticks = [];
    for (let t = -SPAN; t <= SPAN; t += 2) ticks.push(t * SX, 0, z, t * SX, -0.18, z);
    add(hairlines(ticks, 0.5));
    for (const g of GROUPS) {
      const st = steps(card.hist, g, z);
      add(fatLine(st.pts, g === "core" ? { width: 2 } : { width: 1.2, opacity: 0.6 }));
      peaks[`${name}-${g}`] = st.peak;
    }
    add(fatLine([new THREE.Vector3(rel(card, card.lv.value), 0, z),
                 new THREE.Vector3(rel(card, card.lv.value), CARD_H * 0.92, z)],
                { width: 1.2, dashed: true, opacity: 0.85 }));
    cards[name] = { plane, frame, card, z };
  }

  // Camera: orthographic, from the front left and above, so the cards step back to the right.
  const dir = { yaw: -0.28, yawGoal: -0.28, pitch: 0.5, pitchGoal: 0.5 };
  const target = new THREE.Vector3(0, 2, 0);
  let w = 1, h = 1, span = 15;
  const resize = () => {
    w = root.clientWidth; h = root.clientHeight;
    renderer.setSize(w, h, false);
    const phone = w / h < 1;
    span = phone ? 19 : 10; // half-height of the view, in world units
    const a = w / h;
    Object.assign(camera, { left: -span * a, right: span * a, top: span, bottom: -span });
    // Leave room for the copy: on a desktop it sits top left, on a phone across the top.
    if (phone) camera.setViewOffset(w, h, 0, -0.2 * h, w, h);
    else camera.setViewOffset(w, h, -0.17 * w, -0.03 * h, w, h);
    camera.updateProjectionMatrix();
    fat.forEach((l) => l.material.resolution.set(w, h));
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
  };
  new ResizeObserver(resize).observe(root);
  resize();
  const project = (v) => {
    const p = v.clone().project(camera);
    return { x: (p.x + 1) * 0.5 * w, y: (1 - p.y) * 0.5 * h };
  };

  // Drag orbits slightly; hovering a card shows that epoch's numbers.
  let drag = null;
  canvas.addEventListener("pointerdown", (e) => {
    drag = { x: e.clientX, y: e.clientY, yaw: dir.yawGoal, pitch: dir.pitchGoal };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointerup", () => { drag = null; });
  const tip = root.querySelector(".tip");
  const ray = new THREE.Raycaster();
  const ndc = new THREE.Vector2();
  canvas.addEventListener("pointermove", (e) => {
    const r = canvas.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    if (drag) {
      dir.yawGoal = Math.max(-0.55, Math.min(0, drag.yaw - (e.clientX - drag.x) * 0.0015));
      dir.pitchGoal = Math.max(0.3, Math.min(0.75, drag.pitch + (e.clientY - drag.y) * 0.0015));
      return;
    }
    ndc.set((px / w) * 2 - 1, -(py / h) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const hit = ray.intersectObjects([cards.conj.plane, cards.galle.plane])[0];
    for (const c of Object.values(cards)) c.frame.material.opacity = 0.35;
    if (!hit) { tip.hidden = true; return; }
    const c = Object.values(cards).find((x) => x.plane === hit.object);
    c.frame.material.opacity = 0.9;
    tip.hidden = false;
    tip.textContent = cardText(d, c.card);
    tip.style.transform = `translate(${Math.round(Math.min(px + 14, w - 380))}px, ${Math.round(py + 16)}px)`;
  });
  canvas.addEventListener("pointerleave", () => {
    tip.hidden = true;
    for (const c of Object.values(cards)) c.frame.material.opacity = 0.35;
  });

  // Annotations.
  const pt = (x, y, z) => new THREE.Vector3(x, y, z);
  const anchors = () => {
    const phone = PHONE.matches;
    return {
      conj: phone ? { point: pt(-W, 0, Z.conj), offset: [0, 40], marker: "none", leader: false }
        : { point: pt(-W, CARD_H * 0.55, Z.conj), offset: [-30, 0], marker: "none" },
      galle: phone ? { point: pt(-W, CARD_H, Z.galle), offset: [0, -16], marker: "none", leader: false }
        : { point: pt(W, CARD_H * 0.25, Z.galle), offset: [30, 0], marker: "none" },
      // On a phone the branch labels are hidden (CSS), so their markers go too.
      inner: { point: peaks["galle-inner"], offset: phone ? [22, -28] : [40, -40] },
      outer: { point: peaks["galle-outer"], offset: phone ? [-22, -28] : [-40, -40] },
      lv: { point: pt(rel(d.conj, d.conj.lv.value), CARD_H * 0.92, Z.conj), offset: [46, -30] },
      nep: { point: pt(rel(d.conj, d.conj.jpl.value), 0, Z.conj), offset: [-80, 62], marker: "diamond" },
      nepG: { point: pt(rel(d.galle, d.galle.jpl.value), 0, Z.galle), offset: [0, 0], marker: "diamond", leader: false },
      ribbon: { point: pt(rel(d.galle, d.galle.all.hi68) * 0.7 + rel(d.conj, d.conj.all.hi68) * 0.3, 0, Z.galle * 0.4), offset: [16, 0], marker: "none", leader: false },
    };
  };
  const tickEls = [];
  for (const [name, c] of Object.entries(cards)) {
    for (let t = -SPAN / 2; t <= SPAN / 2; t += SPAN / 2) {
      const el = document.createElement("span");
      el.className = "tick";
      el.textContent = t === 0 ? `${fmt(c.card.all.median)}°` : `${t > 0 ? "+" : "−"}${Math.abs(t)}°`;
      root.querySelector(".ticks").append(el);
      tickEls.push({ el, v: pt(t * SX, -0.18, c.z) });
    }
  }

  // Draw-in: the conjunction card is there first; the ribbon then runs back to 1846.
  let visible = false;
  let t0 = null; // the draw-in starts when the panel first scrolls into view
  watchVisibility(root, (v) => { visible = v; if (v && t0 === null) t0 = performance.now(); });
  const DRAW = REDUCED ? 0 : 1800;
  const frame = (now) => {
    requestAnimationFrame(frame);
    if (!visible) return;
    const p = !DRAW || t0 === null ? (DRAW ? 0 : 1) : Math.min(1, Math.max(0, (now - t0) / DRAW));
    const ease = 1 - (1 - p) ** 3;
    setRibbon(ease);
    links.material.opacity = 0.45 * ease;
    dir.yaw += (dir.yawGoal - dir.yaw) * 0.12;
    dir.pitch += (dir.pitchGoal - dir.pitch) * 0.12;
    camera.position.set(
      target.x + 60 * Math.sin(dir.yaw) * Math.cos(dir.pitch),
      target.y + 60 * Math.sin(dir.pitch),
      target.z + 60 * Math.cos(dir.yaw) * Math.cos(dir.pitch));
    camera.lookAt(target);
    renderer.render(scene, camera);
    drawAnnotations(svg, root, anchors(), project, w);
    for (const { el, v } of tickEls) {
      const q = project(v);
      el.style.transform = `translate(${Math.round(q.x)}px, ${Math.round(q.y + 6)}px) translateX(-50%)`;
    }
  };
  requestAnimationFrame(frame);
}

function cardText(d, c) {
  const date = c === d.conj ? d.conjDate.text : d.galleDate.text;
  const helio = c === d.galle ? ` Seen from the Sun the 68 % width is ${fmt(d.galleHelio.width68, 2)}°.` : "";
  return `${date}, ${c.frame}: median ${fmt(c.all.median, 2)}°, 68 % ${range(c.all.lo68, c.all.hi68, 2)}° ` +
    `(width ${fmt(c.all.width68, 2)}°), 95 % ${range(c.all.lo95, c.all.hi95, 2)}°. ` +
    `Neptune (JPL) ${signed(c.jpl.minus_median_deg, 2)}°; Le Verrier's orbit ` +
    `${signed(c.lv.minus_median_deg, 2)}°, ${whereIn(c.lv)} (our calculation).${helio}`;
}

// "inner branch (a ≈ 28.8 AU)" → title "Inner branch", detail "a ≈ 28.8 AU" (kept lowercase:
// a is the semi-major axis, so it must not go through the caps style of the titles).
const branchName = (s) => {
  const m = /^(.*?) \((.*)\)$/.exec(s) ?? [s, s, ""];
  return { title: m[1][0].toUpperCase() + m[1].slice(1), detail: m[2] };
};

/** The panel's text, from the data. */
export function fillConjunction(root, d) {
  const set = (sel, text) => { root.querySelector(sel).textContent = text; };
  const c = d.conj, g = d.galle;
  const yrs = Math.round(d.years);
  set("h2", `The data pin down ${d.conjDate.year}, not ${d.galleDate.year}`);
  set(".body",
    `Uranus overtook the unseen planet around ${d.conjDate.text}, when the pull between them ` +
    `was strongest. There the posterior fixes the planet's longitude to a 68 % width of ` +
    `${fmt(c.all.width68, 2)}°. Carried ${yrs} years forward to Galle's night, it spreads to ` +
    `${fmt(g.all.width68, 2)}°: a small spread in distance becomes a spread in speed. The two eccentric branches are ${fmt(d.split.conj)}° apart at the conjunction and ` +
    `${fmt(d.split.galle)}° apart by ${d.galleDate.year}.`);
  set(".note", `${d.note}.`);
  set("[data-anchor=conj] .t", `${d.conjDate.text} · conjunction`);
  set("[data-anchor=conj] .s", `Seen from the Sun · 68 % width ${fmt(c.all.width68, 2)}°`);
  set("[data-anchor=galle] .t", `${d.galleDate.text} · Galle's night`);
  set("[data-anchor=galle] .s", `Seen from Earth · 68 % width ${fmt(g.all.width68, 2)}°`);
  for (const b of ["inner", "outer"]) {
    const n = branchName(d.groupNames[b]);
    set(`[data-anchor=${b}] .t`, n.title);
    set(`[data-anchor=${b}] .s`, `${n.detail} · ${fmt(100 * c.groups[b].share)} % of the posterior`);
    set(`[data-anchor=${b}] .ts`, n.title.split(" ")[0]);
  }
  set("[data-anchor=lv] .t", `Le Verrier's orbit · ${signed(c.lv.minus_median_deg)}°`);
  set("[data-anchor=lv] .s", `${whereIn(c.lv)[0].toUpperCase()}${whereIn(c.lv).slice(1)} · our calculation`);
  set("[data-anchor=nep] .t", `Neptune (JPL) · ${signed(c.jpl.minus_median_deg, 2)}°`);
  set("[data-anchor=nep] .s", `Its own conjunction with Uranus: ${d.jplConjDate.text}`);
  set("[data-anchor=ribbon] .t", `68 % range, ${d.conjDate.year} → ${d.galleDate.year}`);
  set(".alt",
    `Figure: two translucent cards stand on a floor, ${d.conjDate.text} in front and ` +
    `${d.galleDate.text} behind. Each shows histograms of the planet's longitude (0.25° bins, ` +
    `each group scaled to unit area) for the core of the posterior and its two eccentric ` +
    `branches, centred on that epoch's median (${fmt(c.all.median, 2)}° heliocentric; ` +
    `${fmt(g.all.median, 2)}° geocentric). A ribbon on the floor joins the 68 % ranges, ` +
    `${fmt(c.all.width68, 2)}° wide in front and ${fmt(g.all.width68, 2)}° behind (` +
    `${fmt(d.galleHelio.width68, 2)}° seen from the Sun on Galle's night, so the change of ` +
    `viewpoint is not the reason). The branch ` +
    `medians are ${fmt(d.split.conj, 2)}° apart in ${d.conjDate.year} and ` +
    `${fmt(d.split.galle, 2)}° in ${d.galleDate.year}. Neptune (JPL) is ` +
    `${signed(c.jpl.minus_median_deg, 2)}° and ${signed(g.jpl.minus_median_deg, 2)}° from the ` +
    `medians. Le Verrier's orbit, by our calculation, is ${signed(c.lv.minus_median_deg, 2)}°, ` +
    `${whereIn(c.lv)}, in ${d.conjDate.year}, and ${signed(g.lv.minus_median_deg, 2)}°, ` +
    `${whereIn(g.lv)}, in ${d.galleDate.year}. The other cells (` +
    `${fmt(100 * c.groups.other.share)} % of the posterior) are not drawn. The conjunction date ` +
    `comes from the posterior alone; JPL's own Uranus–Neptune conjunction is ` +
    `${d.jplConjDate.text}. ${d.note}.`);
}
