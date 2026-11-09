// Shared three.js helpers for the panels (docs/DESIGN.md line styles).

import * as THREE from "three";
import { Line2 } from "three/addons/lines/Line2.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineGeometry } from "three/addons/lines/LineGeometry.js";

export const REDUCED = matchMedia("(prefers-reduced-motion: reduce)").matches;
export const PHONE = matchMedia("(max-width: 720px)");

/** A screen-space-width line: 2px solid for main paths, 1.2px dashed for predicted ones. */
export function fatLine(points, { width = 2, opacity = 1, dashed = false, dash = 0.14 } = {}) {
  const geo = new LineGeometry();
  geo.setPositions(points.flatMap((p) => [p.x, p.y, p.z]));
  const mat = new LineMaterial({
    color: 0xffffff, linewidth: width, transparent: opacity < 1, opacity,
    dashed, dashSize: dash, gapSize: dash, worldUnits: false,
  });
  const line = new Line2(geo, mat);
  line.computeLineDistances();
  return line;
}

/** 1px hairlines from a flat [x0, y0, z0, x1, y1, z1, …] segment list. */
export function hairlines(segments, opacity, { fog = true } = {}) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.Float32BufferAttribute(segments, 3));
  const mat = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity, fog });
  return new THREE.LineSegments(geo, mat);
}

/** Radial white-to-clear texture, for the horizon glow. */
export function glowTexture() {
  const c = document.createElement("canvas");
  c.width = c.height = 256;
  const g = c.getContext("2d");
  const grad = g.createRadialGradient(128, 128, 0, 128, 128, 128);
  grad.addColorStop(0, "rgba(255,255,255,0.62)");
  grad.addColorStop(0.35, "rgba(255,255,255,0.22)");
  grad.addColorStop(1, "rgba(255,255,255,0)");
  g.fillStyle = grad;
  g.fillRect(0, 0, 256, 256);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  return tex;
}

/** Calls `set(visible)` as the panel enters or leaves the viewport, so hidden panels skip rendering. */
export function watchVisibility(root, set) {
  new IntersectionObserver(([e]) => set(e.isIntersecting), { threshold: 0.01 }).observe(root);
}

/**
 * Stage annotations: an 8px square (or a diamond) on the point, a 0.7px leader, and the HTML
 * label at an offset. Labels flip sides rather than leave the viewport. `anchors` maps a key to
 * {point, offset:[dx, dy], marker:"dot"|"diamond"|"none", leader:bool}; a label is the element
 * with data-anchor=key inside root.
 */
export function drawAnnotations(svg, root, anchors, project, w) {
  let s = "";
  for (const [k, a] of Object.entries(anchors)) {
    const p = project(a.point);
    const [dx, dy] = a.offset;
    const lx = p.x + dx, ly = p.y + dy;
    const el = root.querySelector(`[data-anchor="${k}"]`);
    const shown = el && el.offsetParent !== null && !a.unanchored;
    if (shown) {
      let left = dx < 0 ? lx - el.offsetWidth : lx;
      if (left + el.offsetWidth > w - 16) left = lx - el.offsetWidth;
      if (left < 16) left = lx;
      el.style.transform = `translate(${Math.round(left)}px, ${Math.round(ly - 9)}px)`;
      if (a.leader !== false) s += `<line x1="${p.x}" y1="${p.y}" x2="${lx}" y2="${ly}" />`;
    }
    if (a.marker === "diamond") {
      s += `<rect class="diamond" x="${p.x - 7}" y="${p.y - 7}" width="14" height="14" transform="rotate(45 ${p.x} ${p.y})" />`;
      s += `<rect class="dot" x="${p.x - 2.5}" y="${p.y - 2.5}" width="5" height="5" />`;
    } else if (a.marker !== "none") {
      s += `<rect class="dot" x="${p.x - 4}" y="${p.y - 4}" width="8" height="8" />`;
    }
  }
  svg.innerHTML = s;
}
