// Loads the pipeline's JSON from site/data/ and turns it into what the panels show.
// No science number is written in this file: every value below is read from a JSON file at load
// time, so re-running the pipeline (M8b, the arbiter) updates the site with no code change.

const DATA = new URL("../data/", import.meta.url);

const cache = new Map();

/**
 * Fetch site/data/<name>, once per page (panels share files). With `optional`, a missing file
 * (404) returns null: a pending result. A network error (not an HTTP status) is retried twice.
 */
export function loadJSON(name, { optional = false } = {}) {
  if (!cache.has(name)) cache.set(name, fetchJSON(name, optional));
  return cache.get(name);
}

async function fetchJSON(name, optional) {
  for (let attempt = 0; ; attempt++) {
    let res;
    try {
      res = await fetch(new URL(name, DATA), { cache: "no-cache" });
    } catch (err) {
      if (attempt >= 2) throw err;
      await new Promise((r) => setTimeout(r, 300 * (attempt + 1)));
      continue;
    }
    if (!res.ok) {
      if (optional && res.status === 404) return null;
      throw new Error(`site/data/${name}: HTTP ${res.status}`);
    }
    return res.json();
  }
}

/** Calendar date (proleptic Gregorian) of a Julian Date, e.g. "23 Sep 1846". */
export function jdToDate(jd) {
  const z = Math.floor(jd + 0.5);
  const alpha = Math.floor((z - 1867216.25) / 36524.25);
  const a = z + 1 + alpha - Math.floor(alpha / 4);
  const b = a + 1524;
  const c = Math.floor((b - 122.1) / 365.25);
  const d = Math.floor(365.25 * c);
  const e = Math.floor((b - d) / 30.6001);
  const day = b - d - Math.floor(30.6001 * e);
  const month = e < 14 ? e - 1 : e - 13;
  const year = month > 2 ? c - 4716 : c - 4715;
  const names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
  return { day, month, year, text: `${day} ${names[month - 1]} ${year}` };
}

/**
 * Status of the 95 % interval. M6c's is provisional (`provisional` lists its bounds); M8 tried to
 * replace it and its rule did not hold (`m8_runs.json` rule_95.replace); M8b decides next
 * (`m8b.json`, absent until it reports).
 */
function interval95(rep, m8, m8b) {
  const g = rep.galle_lon_ecl_of_date_deg;
  const provisional = (rep.provisional ?? []).some((k) => k.endsWith("95"));
  const base = { lo: g.lo95, hi: g.hi95, source: "M6c" };
  if (!provisional) return { ...base, state: "final", note: "" };
  for (const [name, run, key] of [["M8b", m8b, "m8b_95"], ["M8", m8, "m8_95"]]) {
    const rule = run?.rule_95;
    if (rule?.replace) {
      return { lo: rule[key][0], hi: rule[key][1], source: name, state: "final",
               err: rule.sampler_error_std_over_2, note: `replaced by ${name}` };
    }
  }
  if (!m8b) {
    return { ...base, state: "pending",
             note: "provisional: the M8b reruns that decide it have not reported" };
  }
  const failed = Object.entries(m8b.rule_95?.conditions ?? {}).filter(([, ok]) => !ok)
    .map(([k]) => k);
  return { ...base, state: "provisional",
           note: `provisional: M8b did not meet the replacement rule (${failed.join(", ")})` };
}

/** Everything panel 1 (the headline) shows. */
export async function loadHeadline() {
  const [rep, conj, m5a, m6, m8, m8b] = await Promise.all([
    loadJSON("m6c_replacement.json"),
    loadJSON("m6c_conjunction.json"),
    loadJSON("m5a_fits.json"),
    loadJSON("m6_fits.json"),
    loadJSON("m8_runs.json"),
    loadJSON("m8b.json", { optional: true }),
  ]);
  const g = rep.galle_lon_ecl_of_date_deg;
  const h = conj.histograms.galle_night_geo;
  // The histogram must be of the same posterior as the reported interval.
  if (Math.abs(h.centre_deg - g.median) > 1e-6) {
    throw new Error("m6c_conjunction histogram is not centred on the reported median");
  }
  // The observation count is only in m6_fits.json's description ("279 observations, 529
  // measurements"); its measurement count must match the structured one.
  const nm = /\((\d+) observations, (\d+) measurements\)/.exec(m6.description);
  if (!nm || Number(nm[2]) !== rep.best_fit.n_measurements) {
    throw new Error("m6_fits.json: observation count not found or inconsistent");
  }
  const lvPub = m6.comparison_targets.le_verrier;
  const counts = h.edges_rel_deg.slice(0, -1).map((_, i) =>
    Object.values(h.counts).reduce((s, arr) => s + arr[i], 0));
  return {
    median: g.median, lo68: g.lo68, hi68: g.hi68,
    jpl: g.vs.jpl.target,
    leVerrier: g.vs.le_verrier_derived.target,
    mass: rep.mass_earth, a: rep.a_au,
    neptune: { mass: m5a.true_neptune.mass_earth, a: m5a.true_neptune.a_au_osculating_1781 },
    nMeasurements: rep.best_fit.n_measurements,
    nObservations: Number(nm[1]),
    // Le Verrier's printed prediction: heliocentric longitude at his elements' epoch, 1847.0.
    leVerrierPrinted: { lon: lvPub.helio_lon_ecl_of_date_1847_deg, source: lvPub.sources[0] },
    galle: jdToDate(conj.galle_night_jd_tdb),
    hist: { centre: h.centre_deg, edges: h.edges_rel_deg, counts },
    i95: interval95(rep, m8, m8b),
  };
}

export const fmt = (x, n = 1) => x.toFixed(n);
/** Degrees and arcminutes, e.g. 326°32′. */
export const dm = (x) => {
  let d = Math.floor(x), m = Math.round((x - d) * 60);
  if (m === 60) { d += 1; m = 0; }
  return `${d}°${String(m).padStart(2, "0")}′`;
};
export const range = (lo, hi, n = 1) => `${fmt(lo, n)}–${fmt(hi, n)}`;

/** Everything panel 2 (the 1821 conjunction) shows: m6c_conjunction.json. */
export async function loadConjunction() {
  const c = await loadJSON("m6c_conjunction.json");
  const card = (key) => {
    const s = c[key], h = c.histograms[key];
    if (Math.abs(h.centre_deg - s.all.median) > 1e-6) {
      throw new Error(`m6c_conjunction ${key}: histogram not centred on the median`);
    }
    return {
      key, frame: c.frames[key], all: s.all, groups: s.groups, jpl: s.jpl_neptune,
      lv: s.le_verrier, hist: { centre: h.centre_deg, edges: h.edges_rel_deg, counts: h.counts },
    };
  };
  const conj = card("conjunction_helio"), galle = card("galle_night_geo");
  const split = (k) => Math.abs(k.groups.outer.median - k.groups.inner.median);
  return {
    conj, galle,
    galleHelio: c.galle_night_helio.all,
    conjDate: jdToDate(c.epoch.jd_tdb),
    galleDate: jdToDate(c.galle_night_jd_tdb),
    jplConjDate: jdToDate(c.epoch.jpl_uranus_neptune_conjunction_jd_tdb),
    conjSpreadYears: c.epoch.sample_conjunction_years_68_width,
    years: (c.galle_night_jd_tdb - c.epoch.jd_tdb) / 365.25,
    groupNames: c.groups,
    split: { conj: split(conj), galle: split(galle) },
    nSamples: c.n_samples,
    note: c.description.split(":")[0],
  };
}

/** Where a value sits against an interval summary, in words. */
export function whereIn(v) {
  return v.in_68 ? "inside the 68 %" : v.in_95 ? "inside the 95 %, outside the 68 %"
    : "outside the 95 %";
}
export const signed = (x, n = 1) => `${x < 0 ? "−" : "+"}${Math.abs(x).toFixed(n)}`;
