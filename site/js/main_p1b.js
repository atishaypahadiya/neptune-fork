// The companion page (site/adams-challis/): loads site/data/p1b_*.json, fills the text and
// draws Table 1 and Figs. 1 to 4. Each figure fails on its own, with a message in place.
// Every render step replaces what it renders, so running over the prerendered page gives the
// same page (tools/prerender_site.py).

import { loadJSON } from "./data.js";
import * as C from "./copy_p1b.js";
import { fill, fail } from "./plate.js";
import { fig1, fig2, fig3, fig4 } from "./figs_p1b.js";

async function load() {
  const [adams, challis, page, dataset, release] = await Promise.all([
    loadJSON("p1b_adams.json"), loadJSON("p1b_challis.json"), loadJSON("p1b_page.json"),
    loadJSON("m7_dataset.json"), loadJSON("release.json"),
  ]);
  return { adams, challis, page, dataset, release };
}

function setSlots(values) {
  for (const [k, v] of Object.entries(values)) {
    for (const el of document.querySelectorAll(`[data-v="${k}"]`)) fill(el, String(v));
    for (const el of document.querySelectorAll(`[data-href="${k}"]`)) el.href = v;
  }
}

function table1(D, note) {
  const tb = document.querySelector("#table-1 tbody"), tf = document.querySelector("#table-1 tfoot");
  tb.replaceChildren();
  tf.replaceChildren();
  const cell = (tr, text, cls = "") => {
    const td = document.createElement("td");
    td.textContent = text;
    if (cls) td.className = cls;
    tr.append(td);
    return td;
  };
  for (const r of C.caseRows(D)) {
    const tr = document.createElement("tr");
    const name = cell(tr, r.dagger ? `${r.name} †` : r.name);
    // The condition twice: its own column, and under the case on phones (p1b.css).
    const under = document.createElement("span");
    under.className = "cond-under";
    under.textContent = `Conditional on ${r.cond}`;
    name.append(under);
    cell(tr, r.p, "n");
    cell(tr, r.reading);
    cell(tr, r.cond, "cond");
    tb.append(tr);
  }
  const tr = document.createElement("tr");
  const td = cell(tr, note);
  td.colSpan = 4;
  tf.append(tr);
}

function list(sel, items) {
  const ul = document.querySelector(sel);
  ul.replaceChildren();
  for (const t of items) {
    const li = document.createElement("li");
    fill(li, t);
    ul.append(li);
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
      else { const i = document.createElement("i"); i.textContent = p.i; li.append(i); }
    }
    ol.append(li);
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
    for (const el of document.querySelectorAll("[data-fig], .paper-title")) fail(el, err);
    return;
  }
  const challis = C.challisText(D);
  setSlots({
    ...C.titleBlock(D), ...C.abstract(D), ...C.sourcesText(D), ...C.adamsText(D), ...challis,
    ...C.captions(D), ...C.availability(D),
  });
  guard("#table-1", () => table1(D, challis.table_note));
  list(".limits", C.limits(D));
  references(D);
  guard('[data-fig="fig-1"]', (el) => fig1(el, D));
  guard('[data-fig="fig-2"]', (el) => fig2(el, D));
  guard('[data-fig="fig-3"]', (el) => fig3(el, D));
  guard('[data-fig="fig-4"]', (el) => fig4(el, D));
  document.documentElement.dataset.ready = "1";
}

main();
