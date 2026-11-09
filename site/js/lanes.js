// Lanes on a floor: the shared geometry of the honesty and dataset panels. Each lane is a
// hairline rail running left to right on an orthographic floor, back (top of the screen) to
// front. A lane carries bars (2px solid, or dashed for "not yet"), markers drawn in the overlay
// (commit = hollow square, push = solid square, summary = diamond, median = solid square) and,
// optionally, stems that rise from the floor. Walls are translucent planes standing across all
// lanes at one x (launch; Neptune's longitude). Every coordinate is given by the caller from
// the data; this module only draws.

import * as THREE from "three";
import { PHONE, REDUCED, drawAnnotations, fatLine, hairlines, watchVisibility } from "./common.js";

const GAP = 1.6; // world units between lanes
const WALL_H = 1.8;
const BAR_PTS = 48; // points per bar, so bars can draw in

export function mountLanes(root, spec) {
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
  // One group, so a phone can stretch the lane axis (z) and keep rows apart.
  const world = new THREE.Group();
  scene.add(world);
  const add = (o) => { world.add(o); if (o.isLine2) fat.push(o); return o; };
  const V = (x, y, z) => new THREE.Vector3(x, y, z);

  const n = spec.lanes.length;
  const Z = (i) => (i - (n - 1) / 2) * GAP;
  const { xMin, xMax } = spec;
  const zBack = Z(0) - GAP / 2, zFront = Z(n - 1) + GAP / 2;
  const AXIS = zFront + 0.35;

  // Floor outline, rails, axis ticks.
  add(hairlines([xMin, 0, zBack, xMax, 0, zBack, xMax, 0, zBack, xMax, 0, zFront,
                 xMax, 0, zFront, xMin, 0, zFront, xMin, 0, zFront, xMin, 0, zBack], 0.16));
  const rails = spec.lanes.map((_, i) => add(hairlines([xMin, 0, Z(i), xMax, 0, Z(i)], 0.3)));
  add(hairlines((spec.ticks ?? []).flatMap((t) => [t.x, 0, zFront, t.x, 0, zFront + 0.2]), 0.6));

  // Walls.
  for (const w of spec.walls ?? []) {
    const plane = add(new THREE.Mesh(new THREE.PlaneGeometry(zFront - zBack, WALL_H),
      new THREE.MeshBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.06,
                                    side: THREE.DoubleSide, depthWrite: false })));
    plane.rotation.y = Math.PI / 2;
    plane.position.set(w.x, WALL_H / 2, (zBack + zFront) / 2);
    add(fatLine([V(w.x, 0, zBack), V(w.x, 0, zFront)], { width: w.width ?? 1.2, opacity: 0.9 }));
    add(hairlines([w.x, WALL_H, zBack, w.x, WALL_H, zFront], 0.35));
  }

  // Bars: drawn in from x0 to x1 when the panel first shows.
  const bars = [];
  spec.lanes.forEach((lane, i) => {
    for (const b of lane.bars ?? []) {
      const pts = Array.from({ length: BAR_PTS }, (_, k) =>
        V(b.x0 + ((b.x1 - b.x0) * k) / (BAR_PTS - 1), 0.02, Z(i)));
      bars.push(add(fatLine(pts, { width: b.width ?? 2, opacity: b.opacity ?? 1,
                                   dashed: !!b.dashed, dash: b.dash ?? 0.12 })));
    }
  });

  // Stems: one hairline per item rising from the floor (the dataset's observations).
  const stems = [];
  spec.lanes.forEach((lane, i) => {
    for (const s of lane.stems ?? []) stems.push(s.x, 0, Z(i), s.x, s.h, Z(i));
  });
  if (stems.length) add(hairlines(stems, spec.stemOpacity ?? 0.55));

  // Hover: the lane under the pointer, with its text.
  const tip = root.querySelector(".tip");
  const ray = new THREE.Raycaster();
  const ndc = new THREE.Vector2();
  const floor = new THREE.Plane(new THREE.Vector3(0, 1, 0), 0);
  const hit = new THREE.Vector3();
  let hover = -1;

  // Camera: orthographic, high and from the front, so lanes read as rows.
  const view = { yaw: 0, yawGoal: 0, pitch: 0.92, pitchGoal: 0.92 };
  const target = V((xMin + xMax) / 2, 0, (zBack + AXIS) / 2);
  let w = 1, h = 1;
  const resize = () => {
    w = root.clientWidth; h = root.clientHeight;
    renderer.setSize(w, h, false);
    const phone = PHONE.matches;
    const a = w / h;
    const half = spec.half(phone, a, xMax - xMin);
    world.scale.z = phone ? spec.phoneZ ?? 1 : 1;
    target.z = ((zBack + AXIS) / 2) * world.scale.z;
    Object.assign(camera, { left: -half * a, right: half * a, top: half, bottom: -half });
    const [ox, oy] = spec.offset(phone);
    camera.setViewOffset(w, h, ox * w, oy * h, w, h);
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
    drag = { x: e.clientX, y: e.clientY, yaw: view.yawGoal, pitch: view.pitchGoal };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointerup", () => { drag = null; });
  canvas.addEventListener("pointermove", (e) => {
    const r = canvas.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    if (drag && e.pointerType === "mouse") {
      view.yawGoal = Math.max(-0.25, Math.min(0.25, drag.yaw - (e.clientX - drag.x) * 0.0015));
      view.pitchGoal = Math.max(0.6, Math.min(1.25, drag.pitch + (e.clientY - drag.y) * 0.0015));
      return;
    }
    ndc.set((px / w) * 2 - 1, -(py / h) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    hover = -1;
    if (ray.ray.intersectPlane(floor, hit) && hit.x >= xMin && hit.x <= xMax
        && hit.z / world.scale.z >= zBack && hit.z / world.scale.z <= zFront) {
      hover = Math.max(0, Math.min(n - 1, Math.round(hit.z / world.scale.z / GAP + (n - 1) / 2)));
    }
    rails.forEach((rl, i) => { rl.material.opacity = i === hover ? 0.85 : 0.3; });
    if (hover < 0 || !spec.lanes[hover].hover) { tip.hidden = true; return; }
    tip.hidden = false;
    tip.textContent = spec.lanes[hover].hover;
    tip.style.transform = `translate(${Math.round(Math.min(px + 14, w - 380))}px, ${Math.round(py + 18)}px)`;
  });
  canvas.addEventListener("pointerleave", () => {
    hover = -1; tip.hidden = true;
    rails.forEach((rl) => { rl.material.opacity = 0.3; });
  });

  // Overlay: markers per lane, lane labels, wall labels, ticks, extra anchors.
  const laneLabel = (lane) => root.querySelector(`[data-anchor="${lane.key}"]`);
  const tickEls = (spec.ticks ?? []).map((t) => {
    const el = document.createElement("span");
    el.className = t.minor ? "tick minor" : "tick";
    el.textContent = t.text;
    root.querySelector(".ticks").append(el);
    return { el, v: V(t.x, 0, zFront + 0.2) };
  });
  const marker = (p, kind) => {
    const { x, y } = p;
    if (kind === "push" || kind === "median") {
      return `<rect class="dot" x="${x - 4}" y="${y - 4}" width="8" height="8" />`;
    }
    if (kind === "commit") {
      return `<rect class="hollow" x="${x - 3.5}" y="${y - 3.5}" width="7" height="7" />`;
    }
    if (kind === "summary") {
      return `<rect class="diamond" x="${x - 6}" y="${y - 6}" width="12" height="12" transform="rotate(45 ${x} ${y})" />`;
    }
    if (kind === "spot") {
      return `<rect class="dot" x="${x - 3}" y="${y - 3}" width="6" height="6" />`;
    }
    if (kind === "tick") return `<line class="faint" x1="${x}" y1="${y - 6}" x2="${x}" y2="${y + 6}" />`;
    return "";
  };

  let visible = false;
  let t0 = null;
  watchVisibility(root, (v) => { visible = v; if (v && t0 === null) t0 = performance.now(); });
  const DRAW = REDUCED ? 0 : 1600;
  const frame = (now) => {
    requestAnimationFrame(frame);
    if (!visible) return;
    const p = !DRAW ? 1 : t0 === null ? 0 : Math.min(1, Math.max(0, (now - t0) / DRAW));
    for (const b of bars) b.geometry.instanceCount = Math.max(1, Math.ceil((BAR_PTS - 1) * p));
    view.yaw += (view.yawGoal - view.yaw) * 0.12;
    view.pitch += (view.pitchGoal - view.pitch) * 0.12;
    camera.position.set(
      target.x + 60 * Math.sin(view.yaw) * Math.cos(view.pitch),
      target.y + 60 * Math.sin(view.pitch),
      target.z + 60 * Math.cos(view.yaw) * Math.cos(view.pitch));
    camera.lookAt(target);
    renderer.render(scene, camera);

    const phone = PHONE.matches;
    const anchors = spec.anchors ? spec.anchors({ V, Z, phone, zBack, zFront, WALL_H }) : {};
    drawAnnotations(svg, root, anchors, project, w);
    let marks = "";
    if (p >= 1) {
      spec.lanes.forEach((lane, i) => {
        for (const m of lane.marks ?? []) marks += marker(project(V(m.x, m.y ?? 0, Z(i))), m.kind);
      });
    }
    svg.insertAdjacentHTML("beforeend", marks);
    spec.lanes.forEach((lane, i) => {
      const el = laneLabel(lane);
      if (!el || !el.offsetParent) return;
      const q = project(V(xMin, 0, Z(i)));
      if (phone) {
        el.style.transform = `translate(${Math.round(Math.max(16, q.x))}px, ${Math.round(q.y - el.offsetHeight - 4)}px)`;
      } else {
        el.style.transform = `translate(${Math.round(q.x - el.offsetWidth - 16)}px, ${Math.round(q.y - el.offsetHeight / 2)}px)`;
      }
    });
    for (const { el, v } of tickEls) {
      const q = project(v);
      el.style.transform = `translate(${Math.round(q.x)}px, ${Math.round(q.y + 6)}px) translateX(-50%)`;
    }
  };
  requestAnimationFrame(frame);
}
