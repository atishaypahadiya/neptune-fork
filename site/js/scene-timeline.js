// Panel 4a: when each plan went public. One lane per pre-registered run set, oldest at the
// back. Every lane is aligned on its own launch (the wall), so the question "was the plan
// public before the run?" is read as left or right of the wall. Time runs on a signed log
// scale (seconds to hours on one axis), labelled at every tick. Hollow square: plan committed;
// solid square: first pushed (public); bar: the run, launch to finish (dashed while running);
// diamond: first summary of the result. All times from m7_timeline.json.

import { clock, dayMonth, duration, secondsBetween } from "./data.js";
import { mountLanes } from "./lanes.js";

const SX = 2.6; // world units per decade of seconds
const S0 = 20; // seconds: the log scale is linear inside about ±S0
const X = (s) => Math.sign(s) * Math.log10(1 + Math.abs(s) / S0) * SX;
const TICKS = [[-36000, "−10 h"], [-3600, "−1 h"], [-600, "−10 min"], [-60, "−1 min"],
               [0, "launch"], [60, "+1 min"], [600, "+10 min"], [3600, "+1 h"], [36000, "+10 h"]];

/** Facts per entry, all from the JSON: offsets from launch in seconds. */
function facts(d) {
  return d.entries.map((e) => {
    const L = e.launched?.time ?? null;
    const rel = (iso) => (L && iso ? secondsBetween(L, iso) : null);
    return {
      e, L,
      pre: e.pre.map((p) => ({ ...p, dCommit: rel(p.committed), dPush: rel(p.pushed.time) })),
      dFinish: rel(e.finished?.time), dSummary: rel(e.first_summary?.time),
      dSmoke: rel(e.smoke?.time), dNow: rel(d.generated),
    };
  });
}

export function mountTimeline(root, d) {
  const F = facts(d);
  const lanes = F.map((f) => {
    const marks = [], bars = [];
    if (f.L) {
      for (const p of f.pre) {
        marks.push({ x: X(p.dCommit), kind: "commit" }, { x: X(p.dPush), kind: "push" });
      }
      if (f.dSmoke !== null) marks.push({ x: X(f.dSmoke), kind: "tick" });
      if (f.dSummary !== null) marks.push({ x: X(f.dSummary), kind: "summary" });
      if (f.dFinish !== null) bars.push({ x0: 0, x1: X(f.dFinish) });
      else bars.push({ x0: 0, x1: X(f.dNow), dashed: true, width: 1.2 });
    }
    return { key: `lane-${f.e.id}`, marks, bars, hover: laneText(f) };
  });
  const xs = F.flatMap((f) => [...f.pre.flatMap((p) => [p.dCommit, p.dPush]), f.dFinish, f.dNow])
    .filter((v) => v !== null).map(X);
  const xMin = Math.min(...xs, X(-36000)) - 0.3, xMax = Math.max(...xs, X(36000)) + 0.3;
  const head = F.findIndex((f) => f.e.id === "m6c_headline");
  const h = F[head];
  mountLanes(root, {
    lanes, xMin, xMax,
    walls: [{ x: 0 }],
    ticks: TICKS.map(([s, text]) => ({ x: X(s), text, minor: Math.abs(s) > 0 && Math.abs(s) < 3600 })),
    phoneZ: 1.85,
    half: (phone, a, span) => (phone ? (span / 2 + 0.6) / a : 11),
    offset: (phone) => (phone ? [0, -0.25] : [-0.21, -0.245]),
    anchors: ({ V, Z, phone, zBack, WALL_H }) => {
      const z = Z(head);
      const a = {
        wall: { point: V(0, WALL_H, zBack), offset: [0, -12], marker: "none", leader: false },
      };
      if (h?.L) {
        const pre = h.pre[0];
        a.hCommit = { point: V(X(pre.dCommit), 0, z), offset: phone ? [-8, -44] : [-60, -52], marker: "none" };
        a.hPush = { point: V(X(pre.dPush), 0, z), offset: phone ? [6, -64] : [20, -92], marker: "none" };
        if (h.dFinish !== null) {
          a.hFinish = { point: V(X(h.dFinish), 0, z), offset: phone ? [0, 0] : [-30, -50], marker: "none" };
        }
        if (h.dSummary !== null) {
          a.hSummary = { point: V(X(h.dSummary), 0, z), offset: phone ? [0, 0] : [24, -22], marker: "none" };
        }
      }
      return a;
    },
  });
}

function laneText(f) {
  const e = f.e;
  if (!f.L) {
    return `${e.title}: plan committed ${dayMonth(e.pre[0].committed)} ${clock(e.pre[0].committed)}, ` +
      `pushed ${clock(e.pre[0].pushed.time)}. Not launched yet.`;
  }
  const pre = f.pre.map((p) => `${p.role} (${p.sha}) committed ${clock(p.committed)}, pushed ` +
    `${clock(p.pushed.time)}, ${duration(p.dPush)} ${p.dPush < 0 ? "before" : "after"} launch`).join("; ");
  const fin = f.dFinish !== null ? `finished ${clock(e.finished.time)}` : "still running";
  const sum = f.dSummary !== null ? `; first summary ${clock(e.first_summary.time)}` : "";
  const smoke = f.dSmoke !== null ? ` A ${clock(e.smoke.time)} log line is a capped smoke job ` +
    `(31 iterations, no result saved), declared in the pre-registration.` : "";
  return `${e.title}, ${dayMonth(e.launched.time)}: ${pre}. Launched ${clock(e.launched.time)}, ` +
    `${fin}${sum}.${smoke} Sources: git, GitHub activity API, ${e.launched.source}.`;
}

/** The panel's text, from the data. Claims are computed, not asserted. */
export function fillTimeline(root, d) {
  const set = (sel, text) => { root.querySelector(sel).textContent = text; };
  const F = facts(d);
  const h = F.find((f) => f.e.id === "m6c_headline");
  const later = F.filter((f) => f !== h && f.L);
  const main = later.map((f) => f.pre[0]);
  const supplements = later.flatMap((f) => f.pre.slice(1).map((p) => ({ f, p })));
  const allMainBefore = main.every((p) => p.dPush < 0);
  // "Public before the result": every plan pushed before its run finished and before its
  // first summary, for every run set that has finished.
  const beforeResult = F.filter((f) => f.L && f.dFinish !== null).every((f) => f.pre.every((p) =>
    p.dPush < f.dFinish && (f.dSummary === null || p.dPush < f.dSummary)));

  set("h2", beforeResult ? "Every plan was public before its result" : "When each plan went public");
  const hp = h.pre[0];
  const lead = Math.min(...main.map((p) => -p.dPush)), most = Math.max(...main.map((p) => -p.dPush));
  let body =
    `The headline fit's plan was committed at ${clock(hp.committed)} on ${dayMonth(hp.committed)} ` +
    `and its run launched ${duration(hp.dCommit)} later, at ${clock(h.e.launched.time)}. The plan ` +
    `was first pushed at ${clock(hp.pushed.time)}, ${duration(hp.dPush)} after launch: public ` +
    `${duration(h.dFinish - hp.dPush)} before the run finished (${clock(h.e.finished.time)}) and ` +
    `${duration(h.dSummary - hp.dPush)} before the first summary of its result ` +
    `(${clock(h.e.first_summary.time)}).`;
  if (allMainBefore) {
    body += ` Every later plan was pushed before its runs launched, by ${duration(lead)} to ` +
      `${duration(most)}.`;
  }
  for (const { f, p } of supplements) {
    body += ` One supplement to the ${f.e.title.replace(/ \(.*\)$/, "")} plan (${p.role}) was committed ` +
      `${duration(p.dCommit)} before launch and pushed ${duration(p.dPush)} ` +
      `${p.dPush < 0 ? "before" : "after"} it` +
      (f.dFinish !== null && p.dPush < f.dFinish ? ", before the runs finished." : ".");
  }
  set(".body", body);
  set(".note", "Times are PST, from git, GitHub's activity API and the run logs.");

  for (const f of F) {
    const k = `[data-anchor="lane-${f.e.id}"]`;
    set(`${k} .t`, `${f.e.title} · ${dayMonth((f.e.launched ?? f.e.pre[0].pushed).time ?? f.e.pre[0].committed)}`);
    set(`${k} .ts`, f.e.title);
    const pushes = f.pre.map((p) => clock(p.pushed.time)).join(", ");
    set(`${k} .s`, f.L
      ? `Plan public ${pushes} · run ${clock(f.e.launched.time)}–${f.dFinish !== null ? clock(f.e.finished.time) : "running"}`
      : `Plan public ${pushes} · not launched`);
  }
  if (h.L) {
    set('[data-anchor="hCommit"] .t', `Committed ${clock(hp.committed)}`);
    set('[data-anchor="hPush"] .t', `First pushed ${clock(hp.pushed.time)}`);
    set('[data-anchor="hPush"] .s', `${duration(hp.dPush)} after launch`);
    set('[data-anchor="hFinish"] .t', `Finished ${clock(h.e.finished.time)}`);
    set('[data-anchor="hSummary"] .t', `First summary ${clock(h.e.first_summary.time)}`);
  }
  set('[data-anchor="wall"] .t', "Launch");

  set(".alt",
    `Figure: one row per pre-registered run set, oldest at the back, each aligned on its own ` +
    `launch (a wall marked "launch"); time from launch runs left to right on a log scale from ` +
    `ten hours before to ten hours after. A hollow square marks when the plan was committed, a ` +
    `solid square when it was first pushed (public), a bar the run from launch to finish ` +
    `(dashed while running) and a diamond the first summary of its result. ` +
    F.map((f) => laneText(f)).join(" "));
}
