// Panel 5: the dataset. Its 279 observations stand as stems over 1690–1845, one lane per
// observatory (taller stems: RA and Dec both measured). The 28 rows checked against the scan,
// by machine and by hand, carry a square. Beside it: what the dataset is, its two DOIs, how to
// cite it and its licence. Everything from m7_dataset.json (built from the published package
// and Zenodo's records API).

import { fmt } from "./data.js";
import { mountLanes } from "./lanes.js";

const SX = 0.13; // world units per year

function lanesOf(d) {
  const groups = new Map();
  for (const o of d.observations) {
    if (!groups.has(o.observatory)) groups.set(o.observatory, []);
    groups.get(o.observatory).push(o);
  }
  // Back to front by first observation.
  return [...groups.entries()].sort((a, b) => a[1][0].year - b[1][0].year);
}

export function mountDataset(root, d) {
  const stage = root.querySelector(".stage");
  const [y0, y1] = d.counts.years;
  const mid = (y0 + y1) / 2;
  const X = (y) => (y - mid) * SX;
  const L = lanesOf(d);
  const H = (o) => (o.quantities.length > 1 ? 0.9 : 0.55);
  const lanes = L.map(([name, obs]) => ({
    key: `lane-${slug(name)}`,
    stems: obs.map((o) => ({ x: X(o.year), h: H(o) })),
    marks: obs.filter((o) => o.spot_checked).map((o) => ({ x: X(o.year), y: H(o), kind: "spot" })),
    hover: `${name}: ${obs.length} observations, ${fmt(obs[0].year, 0)}–${fmt(obs.at(-1).year, 0)}; ` +
      `${obs.filter((o) => o.spot_checked).length} of them in the spot-checked sample. Observers: ` +
      `${[...new Set(obs.map((o) => o.observer))].join(", ")}.`,
  }));
  const ticks = [];
  for (let y = Math.ceil(y0 / 20) * 20; y <= y1; y += 20) ticks.push({ x: X(y), text: String(y) });
  const spot = d.observations.filter((o) => o.spot_checked);
  const first = spot[Math.floor(spot.length / 2)];
  const li = L.findIndex(([name]) => name === first.observatory);
  mountLanes(stage, {
    lanes, ticks, xMin: X(y0) - 0.4, xMax: X(y1) + 0.4, stemOpacity: 0.6, phoneZ: 1.8,
    half: (phone, a, span) => (phone ? (span / 2 + 0.5) / a : 4.6),
    offset: (phone) => (phone ? [0, -0.06] : [-0.08, -0.04]),
    anchors: ({ V, Z, phone }) => ({
      spot: { point: V(X(first.year), H(first), Z(li)), offset: phone ? [-30, 64] : [40, -58],
              marker: "none" },
    }),
  });
}

const slug = (s) => s.normalize("NFD").replace(/[\u0300-\u036f]/g, "").replace(/[^\w]+/g, "-")
  .toLowerCase();

/** APA-style data citation, from the archive's own metadata. */
function citation(d, which) {
  const z = d.zenodo;
  const names = d.creators.map((c) => {
    const [family, given] = c.split(",").map((x) => x.trim());
    return `${family}, ${given.split(/\s+/).map((g) => `${g[0]}.`).join(" ")}`;
  }).join(", ");
  const year = z.publication_date.slice(0, 4);
  const doi = which === "version" ? z.version_doi : z.concept_doi;
  const ver = which === "version" ? ` (Version ${z.version})` : "";
  return `${names} (${year}). ${d.title}${ver} [Data set]. Zenodo. https://doi.org/${doi}`;
}

function sourceCitation(src) {
  const a = src.authors[0];
  const initials = a["given-names"].split(/\s+/).map((g) => `${g[0]}.`).join(" ");
  return `${a["family-names"]}, ${initials} (${src.year}). ${src.title}. ${src.journal}, ` +
    `${src.start}–${src.end}. ${src.location?.name ?? ""}${src.publisher ? `: ${src.publisher.name}` : ""}.`;
}

/** The panel's text and citation block, from the data. */
export function fillDataset(root, d) {
  const set = (sel, text) => { root.querySelector(sel).textContent = text; };
  const c = d.counts, s = d.spot_checks, z = d.zenodo;
  const span = `${Math.floor(c.years[0])} to ${Math.floor(c.years[1])}`;
  set("h2", "The observations, published on their own");
  set(".body",
    `Table 77 of Le Verrier's 1846 memoir: ${c.rows} observations of Uranus from ${span}, on ` +
    `${c.pages} printed pages, each beside the position he computed and the difference. Two ` +
    `independent transcription passes disagreed on ${c.pass_disagreements} cells, each resolved ` +
    `against a third scan; ${c.misprints} misprints in the print itself are logged with the ` +
    `argument for each. The transcription was machine-assisted, not done by an independent ` +
    `human. ${s.human.rows} rows (${s.human.cells} cells) were then checked against the scan, ` +
    `first automatically and then by hand: ${s.human.cells - s.human.matching} discrepancies.`);
  set(".note", `The other ${c.rows - s.human.rows} rows have had no human check.`);

  const cite = root.querySelector(".cite");
  const link = (sel, doi) => {
    const a = cite.querySelector(sel);
    a.href = `https://doi.org/${doi}`;
    a.textContent = doi;
  };
  link(".doi-concept", z.concept_doi);
  link(".doi-version", z.version_doi);
  cite.querySelector(".ver").textContent = `v${z.version}`;
  cite.querySelector("#cite-concept").textContent = citation(d, "concept");
  cite.querySelector("#cite-version").textContent = citation(d, "version");
  cite.querySelector(".cite-source").textContent = sourceCitation(d.citation.source);
  const rec = cite.querySelector(".record");
  rec.href = z.record_url;
  cite.querySelector(".licence").firstChild.textContent =
    `${z.license === "cc-zero" ? "CC0 1.0: public domain, no conditions" : z.license} · published ` +
    `${z.publication_date} · `;
  for (const b of cite.querySelectorAll("button[data-copy]")) {
    b.addEventListener("click", async () => {
      const text = cite.querySelector(`#${b.dataset.copy}`).textContent;
      try {
        await navigator.clipboard.writeText(text);
        b.textContent = "Copied";
      } catch {
        b.textContent = "Select and copy";
      }
      setTimeout(() => { b.textContent = "Copy"; }, 1600);
    });
  }

  for (const [name, obs] of lanesOf(d)) {
    const k = `[data-anchor="lane-${slug(name)}"]`;
    set(`${k} .t`, name);
    set(`${k} .ts`, `${name} · ${obs.length}`);
    set(`${k} .s`, `${obs.length} observation${obs.length > 1 ? "s" : ""}`);
  }
  set('[data-anchor="spot"] .t', `Checked by hand · ${s.human.rows} rows`);
  set('[data-anchor="spot"] .s', `${s.human.cells} cells, ${s.human.cells - s.human.matching} discrepancies`);
  set(".alt",
    `Figure: the ${c.rows} observations as stems over ${span}, one row per observatory (` +
    lanesOf(d).map(([n, o]) => `${n} ${o.length}`).join(", ") + `); taller stems have both right ` +
    `ascension and declination measured (${c.measurements} measurements in all). The ` +
    `${s.human.rows} spot-checked rows carry a square. Automated check: ${s.automated.cells} cells, ` +
    `${s.automated.cells - s.automated.matching} mismatches; human check: ${s.human.cells} cells, ` +
    `${s.human.cells - s.human.matching} discrepancies.`);
}
