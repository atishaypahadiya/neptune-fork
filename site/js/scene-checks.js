// Panel 4b: what we corrected, and what held. One lane per interval, all against Neptune's
// longitude on Galle's night (the wall, solid: true). Back to front: the first published
// headline (dimmed: superseded), the reported one after the repeatability study, and the three
// pre-registered sensitivity runs. 68 % is the 2px bar, 95 % the thin one (dashed where it is
// provisional), the median a square. Verdicts are the pre-registered reading rules applied to
// the published numbers (site/js/data.js loadChecks).

import { fmt, range } from "./data.js";
import { mountLanes } from "./lanes.js";

const SX = 1.15; // world units per degree

const NAMES = {
  s1_e05: ["S1 · eccentricity up to 0.5", "S1"],
  s2_sigma_data: ["S2 · error levels from the data alone", "S2"],
  s3_t4: ["S3 · heavy-tailed errors", "S3"],
};

function rows(d) {
  const f = d.fits;
  const g = (x) => x.galle_lon_ecl_of_date_deg;
  return [
    { key: "first", iv: d.first, dim: true },
    { key: "reported", iv: d.reported, provisional95: d.reported.i95.state !== "final" },
    ...Object.keys(NAMES).map((k) => ({ key: k, iv: g(f[k]) })),
  ];
}

export function mountChecks(root, d) {
  const R = rows(d);
  const all = R.flatMap((r) => [r.iv.lo95, r.iv.hi95, d.jpl]);
  const c = d.jpl;
  const X = (lon) => (lon - c) * SX;
  const xMin = X(Math.min(...all)) - 0.5, xMax = X(Math.max(...all)) + 0.5;
  const lanes = R.map((r) => ({
    key: `lane-${r.key}`,
    bars: [
      { x0: X(r.iv.lo95), x1: X(r.iv.hi95), width: 1.2, opacity: r.dim ? 0.35 : 0.7,
        dashed: !!r.provisional95, dash: 0.1 },
      { x0: X(r.iv.lo68), x1: X(r.iv.hi68), width: 2, opacity: r.dim ? 0.45 : 1 },
    ],
    marks: [{ x: X(r.iv.median), kind: "median" }],
    hover: hoverText(d, r),
  }));
  const ticks = [];
  for (let lon = Math.ceil(Math.min(...all) / 2) * 2; lon <= Math.max(...all); lon += 2) {
    ticks.push({ x: X(lon), text: `${lon}°` });
  }
  mountLanes(root, {
    lanes, xMin, xMax, ticks, phoneZ: 1.5,
    walls: [{ x: 0, width: 2 }],
    half: (phone, a, span) => (phone ? (span / 2 + 0.6) / a : 8.4),
    offset: (phone) => (phone ? [0, -0.27] : [-0.104, -0.23]),
    anchors: ({ V, zBack, WALL_H }) => ({
      jpl: { point: V(0, WALL_H, zBack), offset: [0, -12], marker: "none", leader: false },
    }),
  });
}

function hoverText(d, r) {
  const iv = r.iv;
  const base = `68 % ${range(iv.lo68, iv.hi68, 2)}° (${fmt(iv.width68, 2)}° wide), 95 % ` +
    `${range(iv.lo95, iv.hi95, 2)}°, median ${fmt(iv.median, 2)}°.`;
  if (r.key === "first") {
    return `First published headline: ${base} Main cell (cell 2) from one run at nlive 1000.`;
  }
  if (r.key === "reported") {
    const rp = d.repeat;
    return `Reported now: ${base} Cell 2 is the merge of two runs at 4 × the live points, which ` +
      `agree on the width to ${fmt(100 * rp.r2RelDiff)} % (rule: within ${fmt(100 * rp.tol, 0)} %). ` +
      (d.reported.i95.state === "final" ? `The 95 % range is final (${d.reported.i95.source}).`
        : `The 95 % range is provisional; ${d.reported.i95.state === "pending"
          ? "a final check is still running" : "M8b did not meet the replacement rule"}.`);
  }
  const v = d.verdicts[r.key];
  return `${NAMES[r.key][0]}: ${base} Pre-registered reading rule: ${v.rule}`;
}

function verdictLine(d, k) {
  const v = d.verdicts[k];
  const f = d.fits[k].vs_headline;
  const moved = `median moved ${fmt(Math.abs(f.galle_median_shift_headline_half_widths), 2)} half-widths`;
  if (k === "s1_e05") {
    const base = `Low-a peak: ${v.word} (P = ${fmt(v.p, 3)})`;
    return v.pending ? `${base} · pending the arbiter` : base;
  }
  if (k === "s2_sigma_data") return `${v.holds ? "Survives" : "Does not survive"} · ${moved}`;
  return `${v.holds ? "Headline unchanged" : "Changes the headline"} · ${moved} · ln B ` +
    `${f.ln_b > 0 ? "+" : "−"}${fmt(Math.abs(f.ln_b))} for it`;
}

/** The panel's text, from the data. */
export function fillChecks(root, d) {
  const set = (sel, text) => { root.querySelector(sel).textContent = text; };
  const rp = d.repeat;
  const inAll = Object.keys(NAMES).every((k) => d.fits[k].galle_lon_ecl_of_date_deg.vs.jpl.in_68)
    && d.reported.vs.jpl.in_68;
  const ratio = rp.r1Spread / rp.r1RunError;
  set("h2", "What we corrected, and what held");
  let body =
    `The first published 68 % interval, ${range(d.first.lo68, d.first.hi68, 2)}° ` +
    `(${fmt(d.first.width68, 2)}° wide), rested on one run of the main cell. Four repeats at the ` +
    `same settings gave widths from ${fmt(Math.min(...rp.r1Runs), 2)}° to ` +
    `${fmt(Math.max(...rp.r1Runs), 2)}°, a spread ${fmt(ratio, 0)} times the sampler's own error. ` +
    `Two runs with four times the live points agreed (${rp.r2Runs.map((x) => `${fmt(x, 2)}°`).join(" and ")}); ` +
    `under the pre-registered rule their merge replaced that cell, and the interval became ` +
    `${range(d.reported.lo68, d.reported.hi68, 2)}° (${fmt(d.reported.width68, 2)}°).`;
  if (inAll) {
    body += " Neptune is inside it, and inside the 68 % interval of all three sensitivity runs, " +
      "each read by a rule written before it ran.";
  }
  const s1 = d.verdicts.s1_e05;
  if (s1.pending) {
    body += ` S1's verdict is pending: the published run gives P(a < 25 AU) = ${fmt(s1.p, 3)}, ` +
      `M8's faster likelihood gives ${fmt(s1.m8p, 3)}, and a full-likelihood arbiter run decides.`;
  } else if (s1.arbiter) {
    body += " The arbiter for S1 has reported.";
  }
  set(".body", body);
  set(".note", "Verdicts are the pre-registered reading rules applied to the published numbers.");

  const label = (key, title, short, sub) => {
    set(`[data-anchor="lane-${key}"] .t`, title);
    set(`[data-anchor="lane-${key}"] .ts`, short);
    set(`[data-anchor="lane-${key}"] .s`, sub);
  };
  label("first", "First published", "First published",
    `68 % ${range(d.first.lo68, d.first.hi68, 2)}° · ${fmt(d.first.width68, 2)}° wide · superseded`);
  label("reported", "Reported now", "Reported now",
    `68 % ${range(d.reported.lo68, d.reported.hi68, 2)}° · ${fmt(d.reported.width68, 2)}° wide · ` +
    (d.reported.i95.state === "final" ? "95 % final" : "95 % provisional (dashed)"));
  for (const [k, [title, short]] of Object.entries(NAMES)) label(k, title, short, verdictLine(d, k));
  set('[data-anchor="jpl"] .t', `Neptune (JPL) · ${fmt(d.jpl, 2)}°`);

  set(".alt",
    `Figure: five rows against a wall at Neptune's longitude on Galle's night, ${fmt(d.jpl, 3)}°. ` +
    `Each row shows a 68 % interval (thick), a 95 % interval (thin) and the median. ` +
    rows(d).map((r) => `${r.key === "first" ? "First published (superseded)" : r.key === "reported"
      ? `Reported now (95 % ${d.reported.i95.state === "final" ? "final" : "provisional"})` : NAMES[r.key][0]}: 68 % ${range(r.iv.lo68, r.iv.hi68, 2)}°, ` +
      `95 % ${range(r.iv.lo95, r.iv.hi95, 2)}°` + (NAMES[r.key] ? `; ${verdictLine(d, r.key)}` : "")
    ).join(". ") + `. Reading rules: ${d.rulesSource}.`);
}
