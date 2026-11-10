// Loads site/data/, fills the text of the paper and draws its plates and figures. Each figure
// fails on its own (with a message in place) if its data are missing or inconsistent.
// Every render step replaces what it renders, so running over the prerendered page
// (tools/prerender_site.py) gives the same page plus the JS-only controls.

import {
  loadJSON, loadHeadline, loadConjunction, loadWindow, loadTimeline, loadChecks, fmt, drv,
} from "./data.js";
import * as C from "./copy.js";
import { fill, fail } from "./plate.js";
import { plateI, plateII, plateIII } from "./plates.js";
import { fig1, fig2, fig4, fig5, fig6 } from "./figs.js";

async function load() {
  const [head, conjL, window_, timeline, checks, rep, conj, m6, dataset, table1, plates, motion, refs,
    quotes, parallax] =
    await Promise.all([
      loadHeadline(), loadConjunction(), loadWindow(), loadTimeline(), loadChecks(),
      loadJSON("m6c_replacement.json"), loadJSON("m6c_conjunction.json"), loadJSON("m6_fits.json"),
      loadJSON("m7_dataset.json"), loadJSON("m7_table1.json"), loadJSON("m7_plates.json"),
      loadJSON("m7_motion.json"), loadJSON("m7_references.json"), loadJSON("m7_quotes.json"),
      loadJSON("m7_parallax.json"),
    ]);
  // The plates' posterior must be the published one (the script checks this too).
  const q = plates.plate_ii_conjunction.posterior.lon_deg, want = conj.conjunction_helio.all;
  if (Math.abs(q.median - want.median) > 1e-3) throw new Error("m7_plates.json does not match m6c_conjunction.json");
  return { head, conjL, window: window_, timeline, checks, rep, conj, m6, dataset, table1, plates, motion, refs,
    quotes, parallax };
}

function setSlots(values) {
  for (const [k, v] of Object.entries(values)) {
    for (const el of document.querySelectorAll(`[data-v="${k}"]`)) fill(el, String(v));
    for (const el of document.querySelectorAll(`[data-href="${k}"]`)) el.href = v;
  }
}

function table1(D) {
  const t = D.table1, tb = document.querySelector("#table-1 tbody"), tf = document.querySelector("#table-1 tfoot");
  tb.replaceChildren();
  tf.replaceChildren();
  const cell = (tr, text, cls = "") => {
    const td = document.createElement("td");
    td.textContent = text;
    if (cls) td.className = cls;
    tr.append(td);
  };
  const years = (y) => (y[0] === y[1] ? `${y[0]}` : `${y[0]}–${y[1]}`);
  for (const r of t.rows) {
    const tr = document.createElement("tr");
    cell(tr, r.source); cell(tr, r.observatory); cell(tr, years(r.years));
    cell(tr, r.observations, "n"); cell(tr, r.measurements.ra, "n"); cell(tr, r.measurements.dec, "n");
    cell(tr, `${fmt(r.median_abs_o_minus_jpl_arcsec, 2)}″`, "n v");
    tb.append(tr);
  }
  const tot = document.createElement("tr");
  cell(tot, "All sources"); cell(tot, ""); cell(tot, years(t.total.years));
  cell(tot, t.total.observations, "n"); cell(tot, t.total.measurements.ra, "n");
  cell(tot, t.total.measurements.dec, "n"); cell(tot, `${fmt(t.total.median_abs_o_minus_jpl_arcsec, 2)}″`, "n v");
  tf.append(tot);
  const z = D.dataset.zenodo, s = D.dataset.spot_checks.human;
  const note = document.createElement("tr");
  const td = document.createElement("td");
  td.colSpan = 7;
  td.textContent = `${s.rows} rows checked by hand against the scan (A6a). ` +
    `${z.license === "cc-zero" ? "CC0 1.0" : z.license}, doi:${z.version_doi}.`;
  note.append(td);
  tf.append(note);
}

/** Archival quotes (m7_quotes.json), each with its original reference and where it is quoted. */
function quotes(D) {
  for (const q of D.quotes.quotes) {
    const el = document.querySelector(`blockquote[data-quote="${q.id}"]`);
    if (!el) continue;
    el.replaceChildren();
    const add = (tag, cls, text, lang) => {
      const e = document.createElement(tag);
      if (cls) e.className = cls;
      if (lang) e.lang = lang;
      e.textContent = text;
      el.append(e);
    };
    add("p", "q-text", `“${q.text}”`, q.translation ? "fr" : "en");
    if (q.translation) add("p", "q-trans", `“${q.translation}” (translation: ${q.translation_by})`);
    add("footer", "q-src", `${q.who}, ${q.original}; quoted in ${q.quoted_in}.`);
  }
}

function references(D) {
  const ol = document.querySelector(".refs");
  ol.replaceChildren();
  for (const r of C.references(D)) {
    const li = document.createElement("li");
    li.id = `ref-${r.id}`;
    for (const p of r.parts) {
      if (typeof p === "string") li.append(p);
      else if (p.i) { const i = document.createElement("i"); i.textContent = p.i; li.append(i); }
      else { const a = document.createElement("a"); a.href = p.href; a.textContent = p.text; li.append(a); }
    }
    ol.append(li);
  }
}

function limits(D) {
  const ul = document.querySelector(".limits");
  ul.replaceChildren();
  for (const t of C.limits(D)) {
    const li = document.createElement("li");
    fill(li, t);
    ul.append(li);
  }
}

function copyButtons() {
  for (const b of document.querySelectorAll(".copy")) {
    b.addEventListener("click", async () => {
      const text = document.getElementById(b.dataset.copy).textContent;
      try { await navigator.clipboard.writeText(text); b.textContent = "Copied"; } catch { b.textContent = "Select and copy"; }
      setTimeout(() => { b.textContent = "Copy"; }, 1600);
    });
  }
}

function guard(sel, f) {
  const el = document.querySelector(sel);
  try { f(el); } catch (err) { fail(el, err); }
}

async function main() {
  let D;
  try {
    D = await load();
  } catch (err) {
    for (const el of document.querySelectorAll("[data-fig]")) fail(el, err);
    return;
  }
  D.head.galleJd ??= D.conj.galle_night_jd_tdb;
  setSlots({
    ...C.titleBlock(D),
    abstract: C.abstract(D),
    cap_plate_i: C.capPlateI(D),
    ...C.dataText(D),
    ...C.prereg(D),
    ...C.motionText(D),
    ...C.modelText(D),
    ...C.resultText(D),
    ...C.conjText(D),
    ...C.windowText(D),
    ...C.robustText(D),
    ...C.availability(D),
  });
  guard("#table-1", () => table1(D));
  limits(D);
  guard("main", () => quotes(D));
  references(D);
  copyButtons();

  guard('[data-fig="plate-i"]', (el) => plateI(el, D));
  guard('[data-fig="fig-1"]', (el) => fig1(el, D));
  guard('[data-fig="fig-2"]', (el) => fig2(el, D));
  const side = document.querySelector("#plate-iii .side");
  const yearEl = side.querySelector(".year"), pullEl = side.querySelector('[data-v="pull_text"]');
  guard('[data-fig="plate-iii"]', (el) => plateIII(el, document.querySelector('[data-fig="fig-3"]'), side, D, (i) => {
    yearEl.textContent = String(drv(Math.floor(D.motion.frames.t_year[i]), "Plate III year"));
    pullEl.textContent = C.pullAt(D.motion, i);
  }));
  guard('[data-fig="fig-4"]', (el) => fig4(el, D));
  guard('[data-fig="plate-ii"]', (el) => plateII(el, D));
  guard('[data-fig="fig-5"]', (el) => fig5(el, D));
  guard('[data-fig="fig-6"]', (el) => fig6(el, D));
  // Controls that only work with JS (Plate III's play and scrub, Copy) are hidden in the HTML.
  for (const el of document.querySelectorAll("[data-js]")) el.hidden = false;
  document.documentElement.dataset.ready = "1";
}

main();
