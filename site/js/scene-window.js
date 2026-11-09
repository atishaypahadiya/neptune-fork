// Panel 3, the lucky window: Le Verrier's planet against the posterior, 1800–1870
// (after-the-fact diagnostic, docs/notes/M6.md). His published elements are held fixed and
// carried through time; the posterior band is measured from its median at each epoch.
// Geometry (docs/DESIGN.md, "the drift is the gap between two paths"): the band lies on the
// floor as a channel (time left to right, offset from the median into depth, ahead = back),
// Le Verrier's path runs across it dashed, Neptune's solid. Stage posts rise at 1836, the
// zero crossing, Galle's night and 1856; a gate cuts the channel at each. A bar in front marks
// the years his planet lies inside the 68 % band. Orthographic, so widths read true.
// Every value comes from `d` (site/js/data.js).

import * as THREE from "three";
import { fmt, signed, whereInBand } from "./data.js";

/** Short form for labels: inside the 68 %, outside the 68 %, or outside the 95 %. */
const bandWord = (v, b68, b95) => (b68[0] <= v && v <= b68[1] ? "inside the 68 %"
  : b95[0] <= v && v <= b95[1] ? "outside the 68 %" : "outside the 95 %");
import { PHONE, REDUCED, drawAnnotations, fatLine, hairlines, watchVisibility } from "./common.js";

const SXT = 0.3; // world units per year
const SZ = 0.3; // world units per degree of offset (ahead = −z, into the screen)
const POST = 2.6; // height of the stage posts

export function mountWindow(root, d) {
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
  // Everything lives in one group, so a phone can stretch the offset axis (z) to fit.
  const world = new THREE.Group();
  scene.add(world);
  const add = (o) => { world.add(o); if (o.isLine2) fat.push(o); return o; };

  const s = d.series;
  const y0 = s.year[0], y1 = s.year.at(-1), ymid = (y0 + y1) / 2;
  const X = (yr) => (yr - ymid) * SXT;
  const Zo = (deg) => -deg * SZ;
  const zMin = Zo(Math.max(...s.hi95_rel_deg)), zMax = Zo(Math.min(...s.lo95_rel_deg));
  const FRONT = zMax + 0.9, AXIS = zMax + 1.6;
  const V = (x, y, z) => new THREE.Vector3(x, y, z);

  // Floor: slab outline and time axis with a tick every decade.
  const W0 = X(y0) - 0.4, W1 = X(y1) + 0.4;
  add(hairlines([W0, 0, zMin - 0.6, W1, 0, zMin - 0.6, W1, 0, zMin - 0.6, W1, 0, AXIS,
                 W1, 0, AXIS, W0, 0, AXIS, W0, 0, AXIS, W0, 0, zMin - 0.6], 0.2));
  const decades = [];
  const ticks = [];
  for (let yr = Math.ceil(y0 / 10) * 10; yr <= y1; yr += 10) {
    ticks.push(X(yr), 0, AXIS, X(yr), 0, AXIS - 0.25);
    decades.push(yr);
  }
  add(hairlines(ticks, 0.6));

  // The bands as strips on the floor, and the median line.
  const strip = (lo, hi, opacity) => {
    const pos = [];
    const idx = [];
    s.year.forEach((yr, i) => {
      pos.push(X(yr), 0, Zo(s[lo][i]), X(yr), 0, Zo(s[hi][i]));
      if (i) idx.push(2 * i - 2, 2 * i - 1, 2 * i, 2 * i - 1, 2 * i + 1, 2 * i);
    });
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.Float32BufferAttribute(pos, 3));
    g.setIndex(idx);
    add(new THREE.Mesh(g, new THREE.MeshBasicMaterial({
      color: 0xffffff, transparent: true, opacity, side: THREE.DoubleSide, depthWrite: false })));
    const edge = (k) => s.year.flatMap((yr, i) => (i ? [X(s.year[i - 1]), 0, Zo(s[k][i - 1]), X(yr), 0, Zo(s[k][i])] : []));
    add(hairlines([...edge(lo), ...edge(hi)], opacity * 2.4));
  };
  strip("lo95_rel_deg", "hi95_rel_deg", 0.07);
  strip("lo68_rel_deg", "hi68_rel_deg", 0.13);
  add(hairlines([X(y0), 0, 0, X(y1), 0, 0], 0.35));

  // Neptune (JPL), solid; Le Verrier's planet, dashed (predicted) and drawn in over time.
  add(fatLine(s.year.map((yr, i) => V(X(yr), 0.01, Zo(s.jpl_minus_median_deg[i]))), { width: 1.2, opacity: 0.9 }));
  const nepAt = Math.round(0.33 * (s.year.length - 1)); // where Neptune's label attaches (display only)
  const lvPts = s.year.map((yr, i) => V(X(yr), 0.02, Zo(s.le_verrier_minus_median_deg[i])));
  const lv = add(fatLine(lvPts, { width: 2, dashed: true, dash: 0.16 }));

  // The window: years with his planet inside the 68 % band.
  const win = [];
  for (const [a, b] of d.in68) win.push(V(X(a), 0, FRONT), V(X(b), 0, FRONT));
  for (let i = 0; i < win.length; i += 2) {
    add(fatLine([win[i], win[i + 1]], { width: 2 }));
    add(hairlines([win[i].x, 0, FRONT - 0.2, win[i].x, 0, FRONT + 0.2,
                   win[i + 1].x, 0, FRONT - 0.2, win[i + 1].x, 0, FRONT + 0.2], 0.9));
  }

  // Stages: a post from his planet's place up, and a gate across the channel.
  const lvAt = (yr) => {
    const i = Math.min(s.year.length - 2, Math.max(0, Math.floor((yr - y0) / d.step)));
    const f = (yr - s.year[i]) / (s.year[i + 1] - s.year[i]);
    const at = (k) => s[k][i] + f * (s[k][i + 1] - s[k][i]);
    return { lv: at("le_verrier_minus_median_deg"), lo68: at("lo68_rel_deg"), hi68: at("hi68_rel_deg"),
             lo95: at("lo95_rel_deg"), hi95: at("hi95_rel_deg"), jpl: at("jpl_minus_median_deg") };
  };
  const stages = {
    s1836: { yr: d.marks.search_1836.year, z: Zo(d.marks.search_1836.le_verrier_minus_median_deg) },
    zero: { yr: d.zero[0], z: 0 },
    galle: { yr: d.marks.galle_night.year, z: Zo(d.marks.galle_night.le_verrier_minus_median_deg) },
    s1856: { yr: d.marks.search_1856.year, z: Zo(d.marks.search_1856.le_verrier_minus_median_deg) },
  };
  const posts = [];
  for (const st of Object.values(stages)) {
    const x = X(st.yr), b = lvAt(st.yr);
    posts.push(x, 0, st.z, x, POST, st.z);
    const g = add(new THREE.Mesh(new THREE.PlaneGeometry(Zo(b.lo95) - Zo(b.hi95), POST * 0.55),
      new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.05,
                                    side: THREE.DoubleSide, depthWrite: false })));
    g.rotation.y = Math.PI / 2;
    g.position.set(x, POST * 0.275, (Zo(b.lo95) + Zo(b.hi95)) / 2);
  }
  add(hairlines(posts, 0.8));

  // Hover: a probe across the channel at the pointer's year, with that year's numbers.
  const probe = add(hairlines([0, 0, zMin, 0, 0, zMax], 0.8));
  probe.visible = false;
  const tip = root.querySelector(".tip");

  // Camera: orthographic, from the front and above.
  const view = { yaw: 0, yawGoal: 0, pitch: 0.78, pitchGoal: 0.78 };
  const target = V(0, 0.6, 0);
  let w = 1, h = 1;
  const resize = () => {
    w = root.clientWidth; h = root.clientHeight;
    renderer.setSize(w, h, false);
    const phone = w / h < 1;
    const a = w / h;
    // Fit 1800–1870 across the width on a phone; there the offset axis is stretched to match.
    const half = phone ? (X(y1) + 0.9) / a : 10;
    lv.material.dashSize = lv.material.gapSize = phone ? 0.36 : 0.16;
    world.scale.z = phone ? 2.4 : 1;
    target.z = ((zMin + AXIS) / 2) * world.scale.z;
    Object.assign(camera, { left: -half * a, right: half * a, top: half, bottom: -half });
    // Desktop: the plot sits right of the equation (bottom left) and below the copy.
    if (phone) camera.setViewOffset(w, h, 0, -0.17 * h, w, h);
    else camera.setViewOffset(w, h, -0.085 * w, -0.17 * h, w, h);
    camera.updateProjectionMatrix();
    fat.forEach((l) => l.material.resolution.set(w, h));
    svg.setAttribute("viewBox", `0 0 ${w} ${h}`);
  };
  new ResizeObserver(resize).observe(root);
  resize();
  const project = (v) => {
    const p = v.clone().multiply(world.scale).project(camera);
    return { x: (p.x + 1) * 0.5 * w, y: (1 - p.y) * 0.5 * h };
  };

  let drag = null;
  canvas.addEventListener("pointerdown", (e) => {
    drag = { x: e.clientX, y: e.clientY, yaw: view.yawGoal, pitch: view.pitchGoal, moved: false };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointerup", () => { drag = null; });
  const ray = new THREE.Raycaster();
  const ndc = new THREE.Vector2();
  const floor = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0);
  const hit = new THREE.Vector3();
  canvas.addEventListener("pointermove", (e) => {
    const r = canvas.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    if (drag && e.pointerType === "mouse") {
      view.yawGoal = Math.max(-0.3, Math.min(0.3, drag.yaw - (e.clientX - drag.x) * 0.0015));
      view.pitchGoal = Math.max(0.5, Math.min(1.1, drag.pitch + (e.clientY - drag.y) * 0.0015));
      return;
    }
    ndc.set((px / w) * 2 - 1, -(py / h) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    if (!ray.ray.intersectPlane(floor, hit)) return;
    const yr = hit.x / SXT + ymid;
    const hz = hit.z / world.scale.z;
    if (yr < y0 || yr > y1 || hz < zMin - 0.6 || hz > AXIS) {
      probe.visible = false; tip.hidden = true; return;
    }
    const b = lvAt(yr);
    probe.geometry.attributes.position.array.set([X(yr), 0, Zo(b.hi95), X(yr), 0, Zo(b.lo95)]);
    probe.geometry.attributes.position.needsUpdate = true;
    probe.visible = true;
    tip.hidden = false;
    tip.textContent = `${fmt(yr, 2)}: posterior 68 % ${signed(b.lo68)}/${signed(b.hi68)}°, ` +
      `95 % ${signed(b.lo95)}/${signed(b.hi95)}° about its median. Le Verrier ${signed(b.lv)}° ` +
      `(${whereInBand(b.lv, [b.lo68, b.hi68], [b.lo95, b.hi95])}); Neptune ${signed(b.jpl, 2)}°.`;
    tip.style.transform = `translate(${Math.round(Math.min(px + 14, w - 380))}px, ${Math.round(py + 18)}px)`;
  });
  canvas.addEventListener("pointerleave", () => { probe.visible = false; tip.hidden = true; });

  // Annotations and axis labels.
  const anchors = () => {
    const phone = PHONE.matches;
    const top = (k) => V(X(stages[k].yr), POST, stages[k].z);
    return {
      s1836: { point: top("s1836"), offset: phone ? [-14, -22] : [-34, -36] },
      zero: { point: top("zero"), offset: phone ? [-6, -54] : [-26, -104] },
      galle: { point: top("galle"), offset: phone ? [12, -46] : [34, -70] },
      // On a phone the 1856 post is short and its top lands on Neptune's line, so the label
      // attaches to the post's foot, on Le Verrier's dashed path, and sits below it.
      s1856: phone ? { point: V(X(stages.s1856.yr), 0, stages.s1856.z), offset: [12, 44] }
        : { point: top("s1856"), offset: [34, -52] },
      lvp: { point: lvPts[20], offset: phone ? [8, -40] : [16, -46], marker: "none" },
      nep: { point: V(X(s.year[nepAt]), 0, Zo(s.jpl_minus_median_deg[nepAt])), offset: phone ? [10, 46] : [24, 64], marker: "none" },
      b68: { point: V(X(y1), 0, Zo(s.hi68_rel_deg.at(-1))), offset: [14, 0], marker: "none", leader: false },
      b95: { point: V(X(y1), 0, Zo(s.hi95_rel_deg.at(-1))), offset: [14, 0], marker: "none", leader: false },
      window: { point: V((win[0].x + win[1].x) / 2, 0, FRONT), offset: [0, 18], marker: "none", leader: false, centre: true },
      ahead: { point: V(X(y0), 0, zMin), offset: [-12, 0], marker: "none", leader: false },
      behind: { point: V(X(y0), 0, zMax), offset: [-12, 0], marker: "none", leader: false },
    };
  };
  const tickEls = decades.map((yr) => {
    const el = document.createElement("span");
    el.className = "tick";
    el.textContent = String(yr);
    root.querySelector(".ticks").append(el);
    return { el, v: V(X(yr), 0, AXIS) };
  });

  let visible = false;
  let t0 = null;
  watchVisibility(root, (v) => { visible = v; if (v && t0 === null) t0 = performance.now(); });
  const DRAW = REDUCED ? 0 : 2200;
  const segs = lvPts.length - 1;
  const frame = (now) => {
    requestAnimationFrame(frame);
    if (!visible) return;
    const p = !DRAW ? 1 : t0 === null ? 0 : Math.min(1, Math.max(0, (now - t0) / DRAW));
    lv.geometry.instanceCount = Math.max(1, Math.ceil(segs * p));
    view.yaw += (view.yawGoal - view.yaw) * 0.12;
    view.pitch += (view.pitchGoal - view.pitch) * 0.12;
    camera.position.set(
      target.x + 60 * Math.sin(view.yaw) * Math.cos(view.pitch),
      target.y + 60 * Math.sin(view.pitch),
      target.z + 60 * Math.cos(view.yaw) * Math.cos(view.pitch));
    camera.lookAt(target);
    renderer.render(scene, camera);
    const a = anchors();
    drawAnnotations(svg, root, a, project, w);
    const wl = root.querySelector('[data-anchor="window"]');
    if (wl.offsetParent) {
      const q = project(a.window.point);
      wl.style.transform = `translate(${Math.round(q.x - wl.offsetWidth / 2)}px, ${Math.round(q.y - wl.offsetHeight - 8)}px)`;
    }
    for (const [k, el] of [["ahead", root.querySelector('[data-anchor="ahead"]')],
                           ["behind", root.querySelector('[data-anchor="behind"]')]]) {
      const q = project(a[k].point);
      el.style.transform = `translate(${Math.round(q.x - el.offsetWidth - 12)}px, ${Math.round(q.y - 9)}px)`;
    }
    for (const { el, v } of tickEls) {
      const q = project(v);
      el.style.transform = `translate(${Math.round(q.x)}px, ${Math.round(q.y + 6)}px) translateX(-50%)`;
    }
  };
  requestAnimationFrame(frame);
}

/** "early 1842", "mid-1842" or "late 1842", from a decimal year. */
const partOfYear = (y) => {
  const f = y - Math.floor(y);
  return `${f < 1 / 3 ? "early " : f < 2 / 3 ? "mid-" : "late "}${Math.floor(y)}`;
};

/** The panel's text, from the data. */
export function fillWindow(root, d) {
  const set = (sel, text) => { root.querySelector(sel).textContent = text; };
  const m = d.marks;
  const [a, b] = d.in68[0] ?? [NaN, NaN];
  const years = Math.round(b - a);
  const miss = (k) => `${fmt(Math.abs(m[k].le_verrier_minus_jpl_deg))}°`;
  const lvWhere = (k) => bandWord(m[k].le_verrier_minus_median_deg, m[k].band_68_rel_deg, m[k].band_95_rel_deg);
  const zero = fmt(d.zero[0], 2);
  const galleIn = a <= m.galle_night.year && m.galle_night.year <= b;
  set("h2", d.in68.length === 1 && galleIn ? `A ${years}-year window, and ${m.galle_night.date.year} was in it`
    : `Le Verrier's planet and the posterior, ${Math.round(d.series.year[0])}–${Math.round(d.series.year.at(-1))}`);
  set(".body",
    `Hold Le Verrier's published orbit fixed and carry it through time. His planet crosses the ` +
    `posterior's median once, in ${partOfYear(d.zero[0])}, and stays inside its 68 % band from ${fmt(a, 2)} to ` +
    `${fmt(b, 2)}. On Galle's night it was ${miss("galle_night")} from Neptune. A search in ` +
    `${m.search_1836.date.year} would have been off by ${miss("search_1836")}; one in ` +
    `${m.search_1856.date.year}, by ${miss("search_1856")}.`);
  set(".note", `${d.note}. It does not model the search fields or star charts.`);
  const label = (k, title, short) => {
    set(`[data-anchor=${k}] .t`, title);
    set(`[data-anchor=${k}] .ts`, short);
  };
  label("s1836", `Search on ${m.search_1836.date.text}`, m.search_1836.date.text);
  set("[data-anchor=s1836] .s", `Le Verrier − Neptune ${signed(m.search_1836.le_verrier_minus_jpl_deg)}° · ${lvWhere("search_1836")}`);
  label("zero", `Crosses the median · ${zero}`, zero);
  set("[data-anchor=zero] .s", "His planet passes the posterior's centre");
  label("galle", `Galle's night · ${m.galle_night.date.text}`, "Galle");
  set("[data-anchor=galle] .s", `Le Verrier − Neptune ${signed(m.galle_night.le_verrier_minus_jpl_deg)}° · ${lvWhere("galle_night")}`);
  label("s1856", `Search on ${m.search_1856.date.text}`, m.search_1856.date.text);
  set("[data-anchor=s1856] .s", `Le Verrier − Neptune ${signed(m.search_1856.le_verrier_minus_jpl_deg)}° · ${lvWhere("search_1856")}`);
  label("lvp", "Le Verrier's planet", "Le Verrier");
  set("[data-anchor=lvp] .s", "His published orbit, held fixed (our calculation)");
  const nepMax = Math.max(...d.series.jpl_minus_median_deg.map(Math.abs));
  label("nep", "Neptune (JPL)", "Neptune");
  set("[data-anchor=nep] .s", `Within ${fmt(nepMax)}° of the median throughout`);
  set("[data-anchor=b68] .t", "68 %");
  set("[data-anchor=b95] .t", "95 %");
  label("window", `Inside the 68 % band · ${fmt(a, 2)}–${fmt(b, 2)}`, `68 % · ${fmt(a, 2)}–${fmt(b, 2)}`);
  set("[data-anchor=ahead] .t", "Ahead");
  set("[data-anchor=behind] .t", "Behind");
  const band = (k) => `${signed(m[k].band_68_rel_deg[0])}/${signed(m[k].band_68_rel_deg[1])}°`;
  set(".alt",
    `Figure: the posterior for the planet's heliocentric longitude, ${Math.round(d.series.year[0])}–` +
    `${Math.round(d.series.year.at(-1))}, measured from its median each quarter-year and laid on ` +
    `the floor as a channel (68 % band brighter, 95 % fainter; time left to right, ahead of the ` +
    `median towards the back). It is narrowest near the 1821 conjunction and widens after 1846. ` +
    `Le Verrier's planet, from his published elements held fixed (our calculation), is a dashed ` +
    `path that starts ${signed(d.series.le_verrier_minus_median_deg[0])}° ahead, crosses the median ` +
    `in ${zero} and ends ${signed(d.series.le_verrier_minus_median_deg.at(-1))}°. Neptune (JPL) ` +
    `is a solid line near the median. A bar marks ${fmt(a, 2)}–${fmt(b, 2)}, the years his planet ` +
    `lies inside the 68 % band` + (d.in95.length ? ` (inside the 95 % band from ` +
    `${fmt(d.in95[0][0], 2)})` : "") + `. Posts mark ${m.search_1836.date.year} (band ` +
    `${band("search_1836")}, Le Verrier ${signed(m.search_1836.le_verrier_minus_median_deg, 2)}°), ` +
    `Galle's night (band ${band("galle_night")}, ${signed(m.galle_night.le_verrier_minus_median_deg, 2)}°) ` +
    `and ${m.search_1856.date.year} (band ${band("search_1856")}, ` +
    `${signed(m.search_1856.le_verrier_minus_median_deg, 2)}°, ` +
    `${fmt(Math.abs(m.search_1856.le_verrier_in_68_half_widths), 2)} half-widths from the median, ` +
    `just past the 68 % edge). Frame: ${d.frame}. ${d.note}.`);
}
