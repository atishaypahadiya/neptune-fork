// Panel 1, the headline: Le Verrier's data → where the unseen planet should be on Galle's night.
// Geometry (docs/DESIGN.md): the ecliptic grid runs to a glowing horizon; the posterior for the
// planet's sky longitude stands on it as a lit curve; a cone of lines falls from the inferred
// planet onto the sky, and its footprint is the 68 % interval (the lines land on posterior
// quantiles 16 %…84 %); Neptune's true place is the diamond. Every value comes from `d`
// (site/js/data.js), read from site/data/ at load time.

import * as THREE from "three";
import { dm, fmt, range } from "./data.js";
import { PHONE, REDUCED, fatLine, glowTexture, hairlines, watchVisibility } from "./common.js";

const SX = 0.9; // world units per degree of longitude
const H = 3.4; // height of the posterior's peak
const SPAN = 9; // degrees shown either side of the histogram centre
const N_CONE = 17; // cone lines, at evenly spaced quantiles across the 68 % interval
const SMOOTH_BINS = 2; // display only: Gaussian kernel, σ = 2 bins (0.5°), on the drawn curve

/**
 * Posterior curve from the histogram: bin centres (deg), heights and the CDF. The CDF (cone,
 * hover) uses the raw counts; only the drawn height is smoothed, by SMOOTH_BINS.
 */
function posterior(hist) {
  const { centre, edges, counts } = hist;
  const total = counts.reduce((s, c) => s + c, 0);
  const kern = Array.from({ length: 6 * SMOOTH_BINS + 1 }, (_, j) =>
    Math.exp(-0.5 * ((j - 3 * SMOOTH_BINS) / SMOOTH_BINS) ** 2));
  const smooth = counts.map((_, i) => {
    let s = 0, w = 0;
    kern.forEach((k, j) => {
      const c = counts[i + j - 3 * SMOOTH_BINS];
      if (c !== undefined) { s += k * c; w += k; }
    });
    return s / w;
  });
  const peak = Math.max(...smooth);
  const bins = counts.map((c, i) => ({
    lon: centre + (edges[i] + edges[i + 1]) / 2, y: smooth[i] / peak,
    lo: centre + edges[i], hi: centre + edges[i + 1], p: c / total,
  }));
  let acc = 0;
  for (const b of bins) { b.cdfLo = acc; acc += b.p; b.cdfHi = acc; }
  const quantile = (q) => {
    const b = bins.find((x) => x.cdfHi >= q) ?? bins.at(-1);
    return b.lo + ((q - b.cdfLo) / Math.max(b.p, 1e-12)) * (b.hi - b.lo);
  };
  const cdfAt = (lon) => {
    const b = bins.find((x) => lon < x.hi) ?? bins.at(-1);
    return b.cdfLo + b.p * Math.min(1, Math.max(0, (lon - b.lo) / (b.hi - b.lo)));
  };
  const heightAt = (lon) => {
    const i = bins.findIndex((x) => x.lon >= lon);
    if (i <= 0) return bins[Math.max(i, 0)].y;
    const a = bins[i - 1], b = bins[i];
    return a.y + ((lon - a.lon) / (b.lon - a.lon)) * (b.y - a.y);
  };
  return { bins: bins.filter((b) => Math.abs(b.lon - centre) <= SPAN), quantile, cdfAt, heightAt };
}

export function mountHeadline(root, d) {
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
  scene.fog = new THREE.Fog(0x000000, 22, 80);
  const camera = new THREE.PerspectiveCamera(36, 1, 0.1, 300);
  const X = (lon) => (lon - d.hist.centre) * SX;
  const post = posterior(d.hist);
  const fatLines = [];

  // Ecliptic grid on y = 0, to a horizon at z = -FAR.
  const FAR = 80, NEAR = 18, W = 120, STEP = 2;
  const grid = [];
  for (let x = -W; x <= W; x += STEP) grid.push(x, 0, NEAR, x, 0, -FAR);
  for (let z = NEAR; z >= -FAR; z -= STEP) grid.push(-W, 0, z, W, 0, z);
  scene.add(hairlines(grid, 0.17));
  scene.add(hairlines([-400, 0, -FAR, 400, 0, -FAR], 0.75, { fog: false }));
  const glow = new THREE.Sprite(new THREE.SpriteMaterial({
    map: glowTexture(), fog: false, depthWrite: false, transparent: true,
    blending: THREE.AdditiveBlending,
  }));
  glow.scale.set(260, 46, 1);
  glow.position.set(0, 2, -FAR - 1);
  scene.add(glow);

  // The posterior, standing on the grid along z = 0: a soft glow, then a crisp 2px line.
  const curvePts = post.bins.map((b) => new THREE.Vector3(X(b.lon), b.y * H, 0));
  const fillPos = [], fillCol = [];
  for (const p of curvePts) {
    fillPos.push(p.x, 0, 0, p.x, p.y, 0);
    fillCol.push(1, 1, 1, 0.02, 1, 1, 1, 0.16 + 0.1 * (p.y / H));
  }
  const fillIdx = [];
  for (let i = 0; i < curvePts.length - 1; i++) {
    const a = 2 * i;
    fillIdx.push(a, a + 2, a + 1, a + 1, a + 2, a + 3);
  }
  const fillGeo = new THREE.BufferGeometry();
  fillGeo.setAttribute("position", new THREE.Float32BufferAttribute(fillPos, 3));
  fillGeo.setAttribute("color", new THREE.Float32BufferAttribute(fillCol, 4));
  fillGeo.setIndex(fillIdx);
  const fillMat = new THREE.MeshBasicMaterial({
    vertexColors: true, transparent: true, side: THREE.DoubleSide, depthWrite: false, opacity: 1,
  });
  scene.add(new THREE.Mesh(fillGeo, fillMat));
  const halo = fatLine(curvePts, { width: 9, opacity: 0.12 });
  const curve = fatLine(curvePts, { width: 2 });
  scene.add(halo, curve);
  fatLines.push(halo, curve);

  // Sky axis along the curve's base, ticks every 2°.
  const lonMin = post.bins[0].lon, lonMax = post.bins.at(-1).lon;
  const axis = [X(lonMin), 0, 0, X(lonMax), 0, 0];
  const ticks = [];
  for (let t = Math.ceil(lonMin / 2) * 2; t <= lonMax; t += 2) {
    axis.push(X(t), 0, 0, X(t), 0, 0.35);
    ticks.push(t);
  }
  scene.add(hairlines(axis, 0.8, { fog: false }));

  // Le Verrier's own elements, carried to that night: dashed (predicted).
  const lv = fatLine([new THREE.Vector3(X(d.leVerrier), 0, 0),
    new THREE.Vector3(X(d.leVerrier), H * 1.08, 0)], { width: 1.2, dashed: true, opacity: 0.85 });
  scene.add(lv);
  fatLines.push(lv);

  // Neptune's true place: solid hairline up to the curve, diamond drawn in the overlay.
  const nepY = post.heightAt(d.jpl) * H;
  scene.add(hairlines([X(d.jpl), 0, 0, X(d.jpl), nepY, 0], 0.9, { fog: false }));

  // The inferred planet: a shaded sphere lit from the horizon, a dashed halo (it is inferred).
  const R = 1.45;
  const planetPos = new THREE.Vector3(X(d.median) + 1.6, 10.6, -8);
  const sphere = new THREE.Mesh(
    new THREE.SphereGeometry(R, 64, 48),
    new THREE.MeshStandardMaterial({ color: 0xd8d8d8, roughness: 0.85, metalness: 0 }),
  );
  sphere.position.copy(planetPos);
  scene.add(sphere);
  scene.add(new THREE.AmbientLight(0xffffff, 0.05));
  const key = new THREE.DirectionalLight(0xffffff, 2.6);
  key.position.set(-0.25, -0.75, 0.6).multiplyScalar(20).add(planetPos);
  key.target = sphere;
  scene.add(key);
  const ringPts = [];
  for (let i = 0; i <= 128; i++) {
    const t = (i / 128) * Math.PI * 2;
    ringPts.push(new THREE.Vector3(Math.cos(t) * R * 1.32, Math.sin(t) * R * 1.32, 0));
  }
  const ring = new THREE.Line(new THREE.BufferGeometry().setFromPoints(ringPts),
    new THREE.LineDashedMaterial({ color: 0xffffff, dashSize: 0.06, gapSize: 0.1,
                                   transparent: true, opacity: 0.55, fog: false }));
  ring.computeLineDistances();
  ring.position.copy(planetPos);
  scene.add(ring);

  // The cone: from the planet to the posterior's 16 %…84 % quantiles on the sky axis.
  const coneEnds = Array.from({ length: N_CONE }, (_, i) =>
    new THREE.Vector3(X(post.quantile(0.16 + (0.68 * i) / (N_CONE - 1))), 0, 0));
  const coneGeo = new THREE.BufferGeometry();
  const conePos = new Float32Array(N_CONE * 6);
  coneGeo.setAttribute("position", new THREE.BufferAttribute(conePos, 3));
  const cone = new THREE.LineSegments(coneGeo, new THREE.LineBasicMaterial({
    color: 0xffffff, transparent: true, opacity: 0.4, fog: false }));
  scene.add(cone);
  const setCone = (p) => {
    coneEnds.forEach((e, i) => {
      const end = planetPos.clone().lerp(e, p);
      conePos.set([planetPos.x, planetPos.y, planetPos.z, end.x, end.y, end.z], i * 6);
    });
    coneGeo.attributes.position.needsUpdate = true;
  };

  // Hover: a hairline at the pointer's longitude, and its place in the posterior.
  const probe = hairlines([0, 0, 0, 0, 1, 0], 0.7, { fog: false });
  probe.visible = false;
  scene.add(probe);
  const tip = root.querySelector(".tip");

  // Overlay anchors (3D points) and HTML labels.
  const anchors = {
    neptune: new THREE.Vector3(X(d.jpl), nepY, 0),
    median: new THREE.Vector3(X(d.median), 0, 0),
    lv: new THREE.Vector3(X(d.leVerrier), H * 1.08, 0),
    planet: planetPos.clone().add(new THREE.Vector3(R * 1.32, 0, 0)),
  };
  const labels = Object.fromEntries(
    [...root.querySelectorAll("[data-anchor]")].map((el) => [el.dataset.anchor, el]));
  const tickEls = ticks.map((t) => {
    const el = document.createElement("span");
    el.className = "tick";
    el.textContent = `${t}°`;
    root.querySelector(".ticks").append(el);
    return { t, el };
  });

  // Camera: eye height e at distance r from the curve, looking just below level, so the
  // horizon sits mid-screen and the grid fills the lower half. Drag orbits it slightly.
  const TILT = 0.02;
  const view = { yaw: -0.05, e: 4.5, r: 22, yawGoal: -0.05, eGoal: 4.5 };
  let w = 1, h = 1;
  const resize = () => {
    w = root.clientWidth; h = root.clientHeight;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    const phone = w / h < 1;
    view.r = phone ? 44 : 22;
    // On a phone the copy takes the top third: shift the image down so the scene sits below it.
    if (phone) camera.setViewOffset(w, h, 0, -0.17 * h, w, h);
    else camera.clearViewOffset();
    camera.updateProjectionMatrix();
    fatLines.forEach((l) => l.material.resolution.set(w, h));
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
  };
  new ResizeObserver(resize).observe(root);
  resize();
  let drag = null;
  canvas.addEventListener("pointerdown", (e) => {
    drag = { x: e.clientX, y: e.clientY, yaw: view.yawGoal, e: view.eGoal };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointerup", () => { drag = null; });
  const project = (v) => {
    const p = v.clone().project(camera);
    return { x: (p.x + 1) * 0.5 * w, y: (1 - p.y) * 0.5 * h, z: p.z };
  };
  canvas.addEventListener("pointermove", (e) => {
    const r = canvas.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    if (drag) {
      const clamp = (v, a) => Math.max(-a, Math.min(a, v));
      view.yawGoal = clamp(drag.yaw - (e.clientX - drag.x) * 0.0012, 0.18);
      view.eGoal = Math.min(6, Math.max(3.6, drag.e + (e.clientY - drag.y) * 0.01));
      return;
    }
    const lvB = project(new THREE.Vector3(X(d.leVerrier), 0, 0));
    const lvT = project(anchors.lv);
    if (Math.abs(px - lvB.x) < 8 && py > lvT.y - 8 && py < lvB.y + 8) {
      probe.visible = false;
      showLvTip(px, py);
      return;
    }
    const a = project(new THREE.Vector3(X(lonMin), 0, 0));
    const b = project(new THREE.Vector3(X(lonMax), 0, 0));
    const top = project(new THREE.Vector3(0, H, 0)).y;
    const f = (px - a.x) / (b.x - a.x);
    if (f < 0 || f > 1 || py < top - 30 || py > a.y + 24) {
      probe.visible = false; tip.hidden = true; return;
    }
    const lon = lonMin + f * (lonMax - lonMin);
    const y = post.heightAt(lon) * H;
    probe.geometry.attributes.position.array.set([X(lon), 0, 0, X(lon), y, 0]);
    probe.geometry.attributes.position.needsUpdate = true;
    probe.visible = true;
    tip.hidden = false;
    tip.textContent = `${fmt(lon, 2)}° · ${Math.round(post.cdfAt(lon) * 100)} % of the posterior below`;
    const q = project(new THREE.Vector3(X(lon), y, 0));
    tip.style.transform = `translate(${Math.round(q.x + 10)}px, ${Math.round(q.y - 34)}px)`;
  });
  canvas.addEventListener("pointerleave", () => { probe.visible = false; tip.hidden = true; });
  // Le Verrier's printed prediction, on his line and on his label.
  function showLvTip(x, y) {
    tip.hidden = false;
    tip.textContent = lvTipText(d);
    tip.style.transform = `translate(${Math.round(Math.min(x + 12, w - 380))}px, ${Math.round(y + 14)}px)`;
  }
  const lvLabel = labels.lv;
  lvLabel.addEventListener("pointerenter", () => {
    const r = lvLabel.getBoundingClientRect(), o = root.getBoundingClientRect();
    showLvTip(r.left - o.left, r.bottom - o.top);
  });
  lvLabel.addEventListener("pointerleave", () => { tip.hidden = true; });

  // Leader lines and stage markers, re-projected each frame.
  const OFFSETS = {
    neptune: [70, -150], median: [-150, 46], lv: [-150, -40], planet: [48, -10],
  };
  const OFFSETS_PHONE = { ...OFFSETS, neptune: [36, -64], lv: [-34, -46] };
  const drawOverlay = () => {
    let s = "";
    for (const [k, v] of Object.entries(anchors)) {
      const p = project(v);
      const [dx, dy] = (PHONE.matches ? OFFSETS_PHONE : OFFSETS)[k];
      const lx = p.x + dx, ly = p.y + dy;
      const el = labels[k];
      // Hidden labels get no leader; on a phone the planet's label is unanchored (CSS).
      const shown = el && el.offsetParent !== null && !(k === "planet" && PHONE.matches);
      if (el && shown) {
        // Labels hang to the side their offset points; flip if that would leave the viewport.
        let left = dx < 0 ? lx - el.offsetWidth : lx;
        if (left + el.offsetWidth > w - 16) left = lx - el.offsetWidth;
        if (left < 16) left = lx;
        el.style.transform = `translate(${Math.round(left)}px, ${Math.round(ly - 9)}px)`;
      }
      if (k !== "planet" && shown) s += `<line x1="${p.x}" y1="${p.y}" x2="${lx}" y2="${ly}" />`;
      if (k === "neptune") {
        s += `<rect class="diamond" x="${p.x - 7}" y="${p.y - 7}" width="14" height="14" transform="rotate(45 ${p.x} ${p.y})" />`;
        s += `<rect class="dot" x="${p.x - 2.5}" y="${p.y - 2.5}" width="5" height="5" />`;
      } else if (k !== "planet") {
        s += `<rect class="dot" x="${p.x - 4}" y="${p.y - 4}" width="8" height="8" />`;
      }
    }
    svg.innerHTML = s;
    for (const { t, el } of tickEls) {
      const p = project(new THREE.Vector3(X(t), 0, 0.35));
      el.style.transform = `translate(${Math.round(p.x)}px, ${Math.round(p.y + 8)}px) translateX(-50%)`;
    }
  };

  // Draw-in: the curve, then the cone. Skipped with prefers-reduced-motion.
  const t0 = performance.now();
  const DRAW = REDUCED ? 0 : 1600;
  const segs = curvePts.length - 1;
  let visible = true;
  watchVisibility(root, (v) => { visible = v; });
  const frame = (now) => {
    if (!visible) { requestAnimationFrame(frame); return; }
    // rAF's timestamp can precede t0 on the first frame: clamp to [0, 1].
    const p = DRAW ? Math.min(1, Math.max(0, (now - t0) / DRAW)) : 1;
    const ease = 1 - (1 - p) ** 3;
    curve.geometry.instanceCount = halo.geometry.instanceCount = Math.ceil(segs * ease);
    fillMat.opacity = ease;
    setCone(DRAW ? Math.min(1, Math.max(0, (now - t0 - DRAW * 0.6) / DRAW)) : 1);
    view.yaw += (view.yawGoal - view.yaw) * 0.12;
    view.e += (view.eGoal - view.e) * 0.12;
    camera.position.set(view.r * Math.sin(view.yaw), view.e, view.r * Math.cos(view.yaw));
    camera.lookAt(0, view.e - view.r * TILT, 0);
    ring.quaternion.copy(camera.quaternion);
    renderer.render(scene, camera);
    drawOverlay();
    requestAnimationFrame(frame);
  };
  requestAnimationFrame(frame);
}

const lvTipText = (d) =>
  `Le Verrier's printed prediction: ${dm(d.leVerrierPrinted.lon)} heliocentric longitude for ` +
  `1 Jan 1847. Source: ${d.leVerrierPrinted.source}.`;

/** The 95 % range as one sentence, from its status (site/js/data.js, interval95). */
function sentence95(i) {
  const r = `${range(i.lo, i.hi)}°`;
  if (i.state === "final") return `The 95 % range is ${r}${i.source === "M6c" ? "" : ` (${i.source})`}.`;
  if (i.state === "pending") return `The 95 % range, ${r}, is provisional; a final check is still running.`;
  return `The 95 % range, ${r}, stays provisional: M8b did not meet the replacement rule.`;
}

/** The panel's text: headline, body, labels, and the text alternative. */
export function fillHeadline(root, d) {
  const set = (sel, text) => { root.querySelector(sel).textContent = text; };
  const m = d.mass, a = d.a, i95 = d.i95;
  set("h1", `Le Verrier's data point to ${fmt(d.median)}°`);
  set(".body",
    `We refitted the ${d.nObservations} observations of Uranus (${d.nMeasurements} measured ` +
    `coordinates) that Le Verrier published in 1846, from first principles. On the night ` +
    `Galle looked, ${d.galle.text}, the unseen planet comes out at ${fmt(d.median)}°, ` +
    `68 % range ${range(d.lo68, d.hi68)}°. Neptune was at ${fmt(d.jpl, 2)}°. ` +
    sentence95(i95));
  set("[data-anchor=neptune] .t", `Neptune, ${d.galle.text}`);
  set("[data-anchor=neptune] .s", `Where Galle found it · ${fmt(d.jpl, 2)}° (JPL)`);
  set("[data-anchor=median] .t", `This rebuild · ${fmt(d.median)}°`);
  set("[data-anchor=median] .s",
    `Posterior median · the cone lands on its 68 %, ${range(d.lo68, d.hi68)}°`);
  set("[data-anchor=lv] .t", `Le Verrier's own orbit · ${fmt(d.leVerrier)}°`);
  set("[data-anchor=lv] .s", "His published elements, carried to that night (our calculation)");
  set("[data-anchor=lv] .ts", "Le Verrier");
  set("[data-anchor=planet] .t", "The unseen planet");
  set("[data-anchor=planet] .s1",
    `${fmt(m.median)} M⊕ · 68 % ${range(m.lo68, m.hi68)} · Neptune ${fmt(d.neptune.mass, 2)}`);
  set("[data-anchor=planet] .s2",
    `${fmt(a.median)} AU · 68 % ${range(a.lo68, a.hi68)} · Neptune ${fmt(d.neptune.a, 2)}`);
  set(".axis-title", `Sky longitude on the night of ${d.galle.text} · ecliptic of date`);
  set(".alt",
    `Figure: the posterior for the unseen planet's geocentric sky longitude on ${d.galle.text}, ` +
    `drawn as a lit curve (smoothed over 0.5° for display) standing on a perspective grid. It peaks near ${fmt(d.median)}° ` +
    `(median); its 68 % range, ${range(d.lo68, d.hi68)}°, is the footprint of a cone of lines ` +
    `falling from a shaded sphere that stands for the planet. Neptune's true place, ` +
    `${fmt(d.jpl, 2)}°, is a diamond almost at the peak. Le Verrier's own published orbit ` +
    `gives ${fmt(d.leVerrier)}° by our calculation, shown as a dashed line. ${lvTipText(d)} ` +
    `${sentence95(i95)} Inferred ` +
    `mass ${fmt(m.median)} Earth masses (68 % ${range(m.lo68, m.hi68)}; Neptune ` +
    `${fmt(d.neptune.mass, 2)}); semi-major axis ${fmt(a.median)} AU, osculating in 1781 ` +
    `(68 % ${range(a.lo68, a.hi68)}; Neptune ${fmt(d.neptune.a, 2)}).`);
}
