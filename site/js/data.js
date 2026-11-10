// Loads the pipeline's JSON from site/data/ and turns it into what the panels show.
// No science number is written in this file: every value below is read from a JSON file at load
// time, so re-running the pipeline (M8b, the arbiter) updates the site with no code change.

const DATA = new URL("../data/", import.meta.url);

/**
 * Numbers computed from JSON values (differences, counts, calendar dates, unit changes) are
 * recorded here with a note, for the number audit (docs/DESIGN.md, the one rule): every number
 * in the rendered text must be a JSON value, a citation, or one of these.
 */
export const derived = [];
const seen = new Set();
export function drv(x, note) {
  const k = `${note}|${x}`;
  if (!seen.has(k)) { seen.add(k); derived.push({ x, note }); }
  return x;
}
globalThis.__derived = derived;

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
export function interval95(rep, m8, m8b) {
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
    galleJd: conj.galle_night_jd_tdb,
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

/** Everything panel 3 (the lucky window, 1800–1870) shows: m6c_offset_track.json. */
export async function loadWindow() {
  const t = await loadJSON("m6c_offset_track.json");
  const s = t.series;
  const n = s.year.length;
  for (const k of Object.keys(s)) {
    if (s[k].length !== n) throw new Error(`m6c_offset_track series ${k}: length ${s[k].length} ≠ ${n}`);
  }
  // Runs of samples where Le Verrier's planet lies inside a band (quarter-year resolution).
  const runs = (lo, hi) => {
    const out = [];
    let start = null;
    for (let i = 0; i <= n; i++) {
      const inside = i < n && s[lo][i] <= s.le_verrier_minus_median_deg[i]
        && s.le_verrier_minus_median_deg[i] <= s[hi][i];
      if (inside && start === null) start = i;
      if (!inside && start !== null) { out.push([s.year[start], s.year[i - 1]]); start = null; }
    }
    return out;
  };
  const mark = (m) => ({ ...m, date: jdToDate(m.jd_tdb) });
  return {
    series: s,
    step: (s.year.at(-1) - s.year[0]) / (n - 1),
    frame: t.frame,
    zero: t.le_verrier_zero_crossings_year,
    marks: Object.fromEntries(Object.entries(t.marks).map(([k, m]) => [k, mark(m)])),
    in68: runs("lo68_rel_deg", "hi68_rel_deg"),
    in95: runs("lo95_rel_deg", "hi95_rel_deg"),
    sources: t.le_verrier_sources,
    note: t.description.split(":")[0],
  };
}

/** Where an offset sits against [lo68, hi68] and [lo95, hi95], in words. */
export function whereInBand(v, b68, b95) {
  if (b68[0] <= v && v <= b68[1]) return "inside the 68 %";
  if (b95[0] <= v && v <= b95[1]) return "inside the 95 %, outside the 68 %";
  return "outside the 95 %";
}

// Times in m7_timeline.json are ISO strings already in PST ("2026-11-07T15:13:24-08:00");
// they are shown as written, never converted to the reader's time zone.
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
export const clock = (iso) => iso.slice(11, 19);
export const dayMonth = (iso) => `${Number(iso.slice(8, 10))} ${MONTHS[Number(iso.slice(5, 7)) - 1]}`;
export const secondsBetween = (a, b) => (Date.parse(b) - Date.parse(a)) / 1000;
/** A duration in words: "5 s", "16 min 48 s", "2 h 50 min". */
export function duration(s) {
  s = Math.round(Math.abs(s));
  if (s < 60) return `${s} s`;
  const m = Math.floor(s / 60);
  if (m < 60) return s % 60 ? `${m} min ${s % 60} s` : `${m} min`;
  const h = Math.floor(m / 60);
  return m % 60 ? `${h} h ${m % 60} min` : `${h} h`;
}

/** Panel 4a: the pre-registration timeline (m7_timeline.json). */
export async function loadTimeline() {
  const t = await loadJSON("m7_timeline.json");
  return { ...t, entries: t.entries.map((e) => ({ ...e, pre: e.preregistration })) };
}

/** Panel 4b: the repeatability correction and the sensitivity runs, with their reading rules. */
export async function loadChecks() {
  const [fits, rep, repeat, timeline, m8, m8b, arbiter] = await Promise.all([
    loadJSON("m6_fits.json"),
    loadJSON("m6c_replacement.json"),
    loadJSON("m6c_repeat.json"),
    loadJSON("m7_timeline.json"),
    loadJSON("m8_runs.json"),
    loadJSON("m8b.json", { optional: true }),
    loadJSON("m8_arbiter.json", { optional: true }),
  ]);
  const rules = timeline.sensitivity_reading_rules;
  const g = (f) => f.galle_lon_ecl_of_date_deg;
  const first = g(fits.fits.headline);
  // The reported interval's 95 % bounds and status follow panel 1 (pending / provisional / final).
  const i95 = interval95(rep, m8, m8b);
  const reported = { ...g(rep), lo95: i95.lo, hi95: i95.hi, i95 };
  const shift = (k) => fits.fits[k].vs_headline.galle_median_shift_headline_half_widths;
  const jplIn = (f) => g(f).vs.jpl.in_68;

  // Each verdict is the pre-registered rule applied to the published numbers.
  const s1p = fits.fits.s1_e05.p_a_below_25_au;
  const r1 = rules.s1_e05;
  const s1Word = s1p > r1.moved_above ? "moved" : s1p < r1.not_at_all_below ? "not at all" : "partly";
  const a8d = m8.tests?.a8d?.s1_e05;
  const keeps = (k) => jplIn(fits.fits[k]) === jplIn(fits.fits.headline)
    && Math.abs(shift(k)) < rules[k].max_shift_headline_half_widths;
  const verdicts = {
    s1_e05: {
      word: s1Word, p: s1p, rule: r1.text,
      // A8d failed for S1, which triggered the arbiter; until it reports the published verdict
      // stands and is marked pending.
      pending: !!(a8d && !a8d.pass) && !arbiter,
      m8p: a8d?.p_a_below_25_au ?? null,
      arbiter,
    },
    s2_sigma_data: { holds: keeps("s2_sigma_data"), rule: rules.s2_sigma_data.text },
    s3_t4: { holds: keeps("s3_t4"), rule: rules.s3_t4.text },
  };
  const r1Runs = ["r1_1", "r1_2", "r1_3", "r1_4"].map((k) => repeat.runs[k].width68);
  const r2Runs = ["r2_1", "r2_2"].map((k) => repeat.runs[k].width68);
  return {
    first, reported, jpl: fits.comparison_targets.jpl.galle_lon_ecl_of_date_deg,
    fits: fits.fits, verdicts, rulesSource: rules.source,
    repeat: {
      r1Runs, r2Runs,
      r1Spread: repeat.spread.r1.width68.std,
      r1RunError: repeat.spread.r1.width68.rms_run_error,
      r2RelDiff: repeat.decision.r2_width_rel_diff,
      tol: repeat.rule.r2_width_tol_rel,
      outcome: repeat.decision.outcome,
    },
  };
}

/** Panel 5: the dataset (m7_dataset.json). */
export const loadDataset = () => loadJSON("m7_dataset.json");
