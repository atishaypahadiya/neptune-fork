// The site's text, in the voice of a paper (docs/DESIGN.md, Voice). Every number is read from
// site/data/ at load time; a number computed from JSON values (a difference, a count, a date
// from a Julian Date) goes through drv(), which records it for the number audit.
// `{X@fig}` in a template is the point letter X linked to its point in that figure.

import {
  fmt, range, signed, jdToDate, clock, dayMonth, duration, secondsBetween, drv,
} from "./data.js";

const MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August",
  "September", "October", "November", "December"];

/** "23 September 1846" from a Julian Date (the parts are registered as derived). */
export function dateOfJd(jd) {
  const d = jdToDate(jd);
  drv(d.day, "calendar day of a JD"); drv(d.year, "calendar year of a JD");
  return `${d.day} ${MONTHS[d.month - 1]} ${d.year}`;
}
/** "28 September 1821, 07:22 TDB" from a JSON date string "1821-09-28 07:22". */
export function dateOfStr(s, { time = true } = {}) {
  const [y, m, d] = s.slice(0, 10).split("-").map(Number);
  const hm = s.slice(11, 16);
  return `${d} ${MONTHS[m - 1]} ${y}${time && hm && hm !== "00:00" ? `, ${hm} TDB` : time ? ", 0h TDB" : ""}`;
}
const deg = (x, n = 1) => `${fmt(x, n)}°`;
const rng = (lo, hi, n = 1) => `${range(lo, hi, n)}°`;
const yr = (x) => String(drv(Math.floor(x), "calendar year of a decimal year"));

/* ---------- Plate I and § 1 ---------- */

export function titleBlock(D) {
  const z = D.dataset.zenodo;
  return {
    author: D.dataset.creators.map((c) => c.split(", ").reverse().join(" ")).join(", "),
    pub_year: z.publication_date.slice(0, 4),
    doi_version: z.version_doi,
    doi_version_url: `https://doi.org/${z.version_doi}`,
    doi_concept: z.concept_doi,
    doi_concept_url: z.concept_url,
    version: z.version,
    subtitle: `A pre-registered Bayesian reconstruction from the ${D.head.nObservations} ` +
      `observations of Uranus he published in ${D.dataset.citation.source.year}`,
  };
}

/** The abstract in three paragraphs: the result in one sentence, the uncertainty (the
 * pre-registered Gaussian range and the Student-t one side by side), then mass, distance and
 * Le Verrier's planet, worded as what the headline fit gives. */
export function abstract(D) {
  const h = D.head, i95 = h.i95, w = D.window;
  const c = D.conj;
  const [a, b] = w.in68[0] ?? [NaN, NaN];
  const off = drv(Math.abs(h.median - h.jpl), "|median - JPL| on Galle's night");
  const t = D.checks.fits.s3_t4, g = t.galle_lon_ecl_of_date_deg;
  const inG = h.lo68 <= h.jpl && h.jpl <= h.hi68;
  const inT = g.vs.jpl.in_68 && g.lo68 <= g.vs.jpl.target && g.vs.jpl.target <= g.hi68;
  const inside = inG && inT ? "inside both" : inG ? "inside the first only" : inT ? "inside the second only" : "outside both";
  return {
    abstract: `Fitted to Le Verrier's own 1846 observations of Uranus, without any knowledge of ` +
      `Neptune, the model places the unseen planet at ${deg(h.median)} on the night Galle found ` +
      `it, ${deg(off, 2)} from Neptune's true position.`,
    abstract_p2: `The longitude is geocentric (ecliptic of date) on ` +
      `${dateOfJd(c.galle_night_jd_tdb)}, from partitioned nested sampling of the ` +
      `${h.nObservations} observations (${h.nMeasurements} measured coordinates, ` +
      `${D.table1.total.years[0]}–${D.table1.total.years[1]}) under a plan committed before the ` +
      `run; Neptune (JPL) was at ${deg(h.jpl, 2)}. The 68 % range is ${rng(h.lo68, h.hi68)} ` +
      `under the pre-registered Gaussian errors and ${rng(g.lo68, g.hi68)} under Student-t ` +
      `errors, which the evidence favours (ln {i:B} = ${signed(t.vs_headline.ln_b, 1)}); ` +
      `Neptune is ${inside}. The 95 % range, ${rng(i95.lo, i95.hi)}, is ` +
      `${i95.state === "final" ? "final" : "provisional"}.`,
    abstract_p3: `The headline fit gives a mass of ${fmt(h.mass.median)} Earth masses ` +
      `(68 %: ${range(h.mass.lo68, h.mass.hi68)}; Neptune ${fmt(h.neptune.mass, 2)}) and a ` +
      `semi-major axis of ${fmt(h.a.median)} AU (${range(h.a.lo68, h.a.hi68)}). Le Verrier's ` +
      `published orbit gives a geocentric longitude of ${deg(h.leVerrier, 2)} that night, inside ` +
      `the 68 % range. In heliocentric longitude (ecliptic of date) his planet lies outside the ` +
      `95 % range at the conjunction of ${yr(w.marks.conjunction.year)} (§ 4.2) and stays within ` +
      `the 68 % band only from ${fmt(a, 2)} to ${fmt(b, 2)} (§ 4.3).`,
  };
}

export function capPlateI(D) {
  const p = D.plates.plate_i_galle_night, q = p.posterior.lon_deg;
  // The same words as the lettering along Uranus's path (plates.js), which phones leave out.
  const names = [...new Set(D.table1.rows.map((r) => r.source.replace(/ series$/, "")))];
  const ty = D.table1.total.years;
  return `The planet on Galle's night, ${dateOfStr(p.date_tdb)}: heliocentric positions ` +
    `projected on the ecliptic of date, seen from the north with the equinox to the right; to ` +
    `scale. {S@plate-i} Sun. {U@plate-i} Uranus (median of the fitted positions) on its dotted ` +
    `path from JPL Horizons, lettered where there is room with the sources of Table 1 ` +
    `(${names.join(" · ")}, ${ty[0]}–${ty[1]}). Stipple: ` +
    `${p.posterior.subsample_xy.length} of the ${D.conj.n_samples} posterior samples of the ` +
    `planet (headline posterior). Their heliocentric longitude, 68 %: ` +
    `${rng(q.lo68, q.hi68)}, is not comparable with the geocentric ranges in the text. ` +
    `{N@plate-i} Neptune (JPL Horizons). {L@plate-i} Le Verrier's planet from his published ` +
    `elements (Krajnović 2016, after Grosser 1962).`;
}

/** A row title inside a sentence: first letter lower case ("Headline fit" → "headline fit"). */
const lc = (t) => t[0].toLowerCase() + t.slice(1);

/* ---------- § 2 Data ---------- */

export function dataText(D) {
  const c = D.dataset.counts, s = D.dataset.spot_checks.human;
  const unchecked = drv(c.rows - s.rows, "rows without a human check");
  const bad = drv(s.cells - s.matching, "cells in disagreement");
  return {
    data_p1: `Table 77 of Le Verrier's memoir (Le Verrier 1846) lists ${c.rows} observations of ` +
      `Uranus from ${yr(c.years[0])} to ${yr(c.years[1])} on ${c.pages} printed pages, each beside ` +
      `the position he computed. The table was transcribed twice from the scan and the two ` +
      `transcriptions reconciled, each disagreement resolved against the page image. ${s.rows} ` +
      `rows were then checked by hand against the scan, with ${bad === 0 ? "no" : bad} ` +
      `discrepancies; the remaining ${unchecked} rows have had no human check. The two ` +
      `transcriptions disagreed on ${c.pass_disagreements} cells, and ${c.misprints} misprints in ` +
      `the print are logged with the argument for each. The transcription is published as a dataset ` +
      `(${D.dataset.creators[0].split(",")[0]} ${D.dataset.zenodo.publication_date.slice(0, 4)}).`,
    data_p2: `Table 1 groups the observations by source: the four observers before Herschel, ` +
      `whom the error model pools into one group, and the three modern series that are its ` +
      `other groups. The last column is each source's median distance from JPL's Uranus taken ` +
      `through our observation model (§ 3.3). These residuals set the error level of each group ` +
      `in the likelihood; the error-level sensitivity test replaces them with levels estimated from the data ` +
      `alone (§ 5).`,
    cap_table_1: `Le Verrier's Table 77 by source. n: observations (rows). RA, Dec: measured ` +
      `coordinates. |O − JPL|: median absolute on-sky residual against JPL's Uranus at the ` +
      `apparent place of date, in arcseconds (right ascension as Δα cos δ). Built from ` +
      `the published transcription; the per-group medians reproduce the committed error model.`,
    cap_fig_1: `The ${c.rows} observations of Table 77 by source and year, one dot per ` +
      `observation, spread up and down within its row so that close dots stay apart (the ` +
      `spread encodes nothing); the ${s.rows} rows checked by hand are open circles. Tap or ` +
      `hover for the date and the residuals against JPL.`,
  };
}

/* ---------- § 3.1 Pre-registration (Fig. 2) ---------- */

export function timelineFacts(t) {
  return t.entries.map((e) => {
    const L = e.launched?.time ?? null;
    const rel = (iso) => (L && iso ? secondsBetween(L, iso) : null);
    return {
      e, L,
      pre: e.pre.map((p) => ({ ...p, dCommit: rel(p.committed), dPush: rel(p.pushed.time) })),
      dFinish: rel(e.finished?.time), dSummary: rel(e.first_summary?.time),
      dSmoke: rel(e.smoke?.time), dNow: rel(t.generated),
    };
  });
}

export function prereg(D) {
  const F = timelineFacts(D.timeline);
  const h = F.find((f) => f.e.id === "m6c_headline");
  const later = F.filter((f) => f !== h && f.L);
  const main = later.map((f) => f.pre[0]);
  const hp = h.pre[0];
  const dur = (s, note) => duration(drv(Math.round(Math.abs(s)), note));
  let p2 = `The headline fit's plan was committed at ${clock(hp.committed)} on ` +
    `${dayMonth(hp.committed)} and its run launched ${dur(hp.dCommit, "s, commit to launch")} later, ` +
    `at ${clock(h.e.launched.time)}. The plan was first pushed at ${clock(hp.pushed.time)}, ` +
    `${dur(hp.dPush, "s, launch to push")} after launch: public ` +
    `${dur(h.dFinish - hp.dPush, "s, push to finish")} before the run finished ` +
    `(${clock(h.e.finished.time)}) and ${dur(h.dSummary - hp.dPush, "s, push to summary")} before ` +
    `the first summary of its result (${clock(h.e.first_summary.time)}). The logs cannot rule out ` +
    `that one of its result files was opened by hand during the run.`;
  if (main.length && main.every((p) => p.dPush < 0)) {
    const lead = Math.min(...main.map((p) => -p.dPush)), most = Math.max(...main.map((p) => -p.dPush));
    p2 += ` Every later plan was pushed before its runs launched, by ` +
      `${dur(lead, "s, shortest push lead")} to ${dur(most, "s, longest push lead")}.`;
  }
  const running = F.filter((f) => f.L && f.dFinish === null).map((f) => f.e.title);
  const waiting = F.filter((f) => !f.L).map((f) => f.e.title);
  return {
    prereg_p1: `Every threshold, reading rule and run setting below was committed to the ` +
      `repository before the run it governs, and each commit's first push is a public record. ` +
      `Fig. 2 places each run set on its own launch time. Times are PST, from git, GitHub's ` +
      `activity API and the run logs.`,
    prereg_p2: p2,
    cap_fig_2: `When each plan went public. One row per pre-registered run set, each aligned on ` +
      `its own launch (the vertical line); time on a signed logarithmic scale. Open square: plan ` +
      `committed; filled square: first pushed; bar: the run, launch to finish` +
      `${running.length ? ` (dashed while running: the ${running.map(lc).join(", the ")})` : ""}; diamond: first ` +
      `summary of the result.${waiting.length ? ` Not yet launched: ${waiting.join(", ")}.` : ""} ` +
      `From the repository's timeline file, generated ${dayMonth(D.timeline.generated)} ` +
      `${clock(D.timeline.generated)}.`,
  };
}

/* ---------- § 3.2 Plate III and Fig. 3 ---------- */

export function motionText(D) {
  const m = D.motion, c = m.conjunction.grid_10_day, s = m.along_track_sign_change;
  const v = m.validation.uranus;
  const last = m.frames.t_year.at(-1);
  const kDays = drv(s.jd_tdb - c.jd_tdb, "days, K to the along-track sign change");
  return {
    pull_intro: `Plate III shows the perturbation Le Verrier inverted, computed directly. The Sun ` +
      `and the eight planets are integrated from JPL's states at Herschel's discovery, once with ` +
      `Neptune's mass and once without. {K@plate-iii} marks the heliocentric conjunction of ` +
      `Uranus and Neptune, ${dateOfStr(c.date_tdb)} (ecliptic of date; JPL's date to ` +
      `${fmt(Math.abs(c.minus_jpl_hours), 1)} h). It is Neptune's own conjunction, not ` +
      `§ 4.2's {C@plate-ii}, the conjunction with the fitted planet. Before {K@plate-iii} ` +
      `Neptune's pull has a component along Uranus's motion that draws it ahead; from ` +
      `${dateOfStr(s.date_tdb)}, ${fmt(kDays, 1)} days after {K@plate-iii}, it holds Uranus ` +
      `back. By ${fmt(last, 1)} Uranus is ${fmt(Math.abs(m.summary.dlon_end_arcsec), 1)}″ behind ` +
      `where it would be without Neptune.`,
    pull_open_circle: `Open circle: Uranus as it would move without Neptune, from the same ` +
      `${yr(m.frames.t_year[0])} state; its offset is drawn ${m.offset_exaggeration} times ` +
      `larger than life.`,
    cap_plate_iii: `Uranus ({U@plate-iii}) and Neptune ({N@plate-iii}) about the Sun ` +
      `({S@plate-iii}): heliocentric, ecliptic and mean equinox of J2000, seen from the north ` +
      `with the equinox to the right; to scale, except the open circle (offset ×` +
      `${m.offset_exaggeration}). ${m.integrator.code} ${m.integrator.version}, ` +
      `${m.integrator.scheme} (Rein & Liu 2012; Rein & Spiegel 2015): the Sun and eight planets ` +
      `from JPL Horizons states at ${dateOfStr(m.initial_states.date)}, run with and without ` +
      `Neptune's mass, drawn every ${m.step_years} yr to ${fmt(last, 2)}. The line ` +
      `{U@plate-iii}{N@plate-iii} carries Neptune's pull. {K@plate-iii}: the heliocentric ` +
      `conjunction of Uranus and Neptune in the run, ${dateOfStr(c.date_tdb)} (ecliptic of date), ` +
      `${fmt(Math.abs(c.minus_jpl_hours), 1)} h from JPL's; drawn at Neptune's place on that date. ` +
      `It is not Plate II's {C@plate-ii}. Horizons check: at ` +
      `${dateOfStr(v.last_date_tdb, { time: false })}, the last cached epoch (none reaches ` +
      `${D.motion.end.label}), the run's Uranus differs from Horizons' by ${v.final_arcsec}″` +
      (v.final_arcsec === v.max_arcsec ? `, the largest difference over all ${v.n_epochs} cached epochs checked.`
        : `; the largest over all ${v.n_epochs} cached epochs is ${v.max_arcsec}″.`),
    cap_fig_3: `Uranus's heliocentric longitude (ecliptic J2000) with Neptune minus without, in ` +
      `arcseconds, ${yr(m.frames.t_year[0])}–${fmt(last, 1)}. The cursor follows the plate's ` +
      `year; {K@fig-3} the conjunction of Uranus and Neptune ({K@plate-iii} in Plate III).`,
  };
}

/** The line beside Plate III for the frame on show. */
export function pullAt(m, i) {
  const f = m.frames, a = f.accel_along_track[i], dl = f.dlon_arcsec[i];
  const pull = a >= 0 ? "drawn ahead" : "held back";
  const where = Math.abs(dl) < 0.05 ? "on its Neptune-free path"
    : `${fmt(Math.abs(dl), 1)}″ ${dl > 0 ? "ahead of" : "behind"} its Neptune-free path`;
  return `Uranus is ${pull} by Neptune's pull, and ${where}.`;
}

/* ---------- § 3.3 Model ---------- */

export function modelText(D) {
  const fits = D.m6, p = fits.priors, em = fits.error_model, s = fits.sampler;
  const cells = D.rep.evidence.cells;
  const sig = Object.values(em.sigma_arcsec), off = Object.values(em.offset_prior_arcsec);
  const mlo = drv(10 ** p.log10_m_earth[0], "10^log10 m lower"), mhi = drv(10 ** p.log10_m_earth[1], "10^log10 m upper");
  return {
    model_p1: `Each measured right ascension and declination is predicted by integrating the Sun ` +
      `(with the inner planets folded in), Jupiter, Saturn, Uranus and one unseen planet with our ` +
      `own Wisdom–Holman integrator, then carrying Uranus through light time, annual aberration, ` +
      `precession and nutation (Capitaine et al. 2003) to the apparent place of date, the frame ` +
      `of Table 77. Jupiter and Saturn start from JPL states and are not fitted. Uranus's four ` +
      `in-plane elements are free; its inclination and node are held at their ` +
      `${yr(D.motion.frames.t_year[0])} values, which positions in longitude cannot constrain. The ` +
      `planet has seven parameters: mass, semi-major axis, eccentricity, two orbital angles, and ` +
      `its inclination and node relative to Uranus's plane.`,
    model_p2: `The likelihood is Gaussian in ${em.groups.length} error groups (the ancient ` +
      `observers, Greenwich, Paris 1801–1827 and Paris 1835–1845, each in RA and Dec), with σ ` +
      `from ${Math.min(...sig)}″ to ${Math.max(...sig)}″ calibrated against JPL and a free offset ` +
      `per group (prior σ ${Math.min(...off)}–${Math.max(...off)}″). Priors: mass log-uniform over ` +
      `${mlo}–${mhi} Earth masses; semi-major axis log-uniform over ` +
      `${p.a_au_log_uniform[0]}–${p.a_au_log_uniform[1]} AU; eccentricity uniform up to ` +
      `${p.e_uniform_max}; inclination to Uranus's plane up to ${p.inc_rel_uranus_deg_max}°, ` +
      `${p.inc}; the angles ${p.angles}. The posterior is sampled with dynesty (Speagle 2020) in ` +
      `${cells.length} cells (` +
      cells.map((c) => `${c.a_au[0]}–${c.a_au[1]} AU${c.e ? ` with {i:e} ${c.e[0]}–${c.e[1]}` : ""}`).join("; ") +
      `) with ` +
      `${s.nlive_default} live points (${Object.values(s.nlive_cells)[0]} in the first), and the ` +
      `cells are weighted by their evidence. The code path of the fit never reads Neptune; a ` +
      `static check of its call graph and a runtime test confirm it.`,
  };
}

/* ---------- § 4.1 ---------- */

export function resultText(D) {
  const h = D.head, i95 = h.i95, rep = D.rep;
  const g = rep.galle_lon_ecl_of_date_deg;
  const lvOff = g.vs.le_verrier_derived.median_minus_target_deg;
  const s95 = i95.state === "final" ? `The 95 % range is ${rng(i95.lo, i95.hi)}.`
    : `The 95 % range, ${rng(i95.lo, i95.hi)}, is ${i95.note}.`;
  return {
    res_p1: `Every longitude and difference in this section is geocentric, ecliptic of date. On ` +
      `Galle's night, ${dateOfJd(h.galleJd)}, we find the planet's longitude at ` +
      `${deg(h.median)}, 68 % range ${rng(h.lo68, h.hi68)}. ` +
      `Neptune was at ${deg(h.jpl, 2)} (JPL Horizons), ` +
      `${fmt(Math.abs(g.vs.jpl.median_minus_target_deg), 2)}° from the median, ` +
      `${fmt(Math.abs(g.vs.jpl.offset_sigma), 2)} of the 68 % half-width. ${s95} Le Verrier's ` +
      `published elements, carried to that night (our calculation; Kepler orbit, equinox 1847.0 ` +
      `assumed), give {L@fig-4} at ${deg(h.leVerrier, 2)}, ${fmt(Math.abs(lvOff), 2)}° from our ` +
      `median (geocentric), inside the 68 % range.`,
    res_p2: `The same posterior gives a mass of ${fmt(h.mass.median)} Earth masses (68 %: ` +
      `${range(h.mass.lo68, h.mass.hi68)}; Neptune ${fmt(h.neptune.mass, 2)}), a semi-major axis ` +
      `of ${fmt(h.a.median)} AU (${range(h.a.lo68, h.a.hi68)}), osculating in ` +
      `${yr(D.motion.frames.t_year[0])}, and an eccentricity below ${rep.e.upper68_one_sided} ` +
      `(68 %, one-sided). The inclination to Uranus's plane is dominated by its prior, and we ` +
      `do not report it as a result. The sampler's own systematic on the median is ` +
      `±${fmt(rep.sampler_systematic.galle_median_deg, 2)}°, not folded into the range.`,
    cap_fig_4: `Posterior of the planet's geocentric longitude on ${dateOfJd(h.galleJd)}, ` +
      `ecliptic of date: histogram of the ${D.conj.n_samples} derived samples of the ` +
      `headline posterior. Solid lines: median and 68 % range; dashed: the 95 % range ` +
      `(${i95.state === "final" ? "final" : "provisional"}). {N@fig-4} Neptune (JPL Horizons). ` +
      `{L@fig-4} Le Verrier's elements carried to that night.`,
  };
}

/* ---------- § 4.2 Plate II ---------- */

export function conjText(D) {
  const c = D.conj, p = D.plates.plate_ii_conjunction, s = c.conjunction_helio;
  const lv = s.le_verrier;
  const days = drv(Math.abs(c.epoch.jd_tdb - c.epoch.jpl_uranus_neptune_conjunction_jd_tdb), "days, C minus JPL conjunction");
  return {
    conj_year: c.epoch.date_tdb.slice(0, 4),
    conj_p1: `On ${dateOfStr(c.epoch.date_tdb)}, Uranus ({U@plate-ii}) and the planet were in ` +
      `heliocentric conjunction: the line {S@plate-ii}{U@plate-ii} produced passes through ` +
      `{C@plate-ii}, the median of the planet's heliocentric longitude, ${deg(s.all.median, 2)} ` +
      `(ecliptic of date; 68 %: ${rng(s.all.lo68, s.all.hi68, 2)}). Neptune ({N@plate-ii}, JPL) was ` +
      `at ${deg(s.jpl_neptune.value, 2)}. Le Verrier's orbit puts his planet at {L@plate-ii}, ` +
      `${deg(lv.value, 2)}: the angle {C@plate-ii}{S@plate-ii}{L@plate-ii} is ` +
      `${deg(Math.abs(lv.minus_median_deg), 2)}, ${fmt(Math.abs(lv.in_half_widths), 2)} ` +
      `half-widths of our 68 % range, ${lv.in_95 ? "inside" : "outside"} the 95 % range.`,
    conj_p2: `The date of {C@plate-ii} is the median, over the ${c.n_samples} posterior samples, of ` +
      `the date on which each sample's planet and Uranus share a heliocentric longitude (ecliptic ` +
      `of date); 68 % of those dates fall within ${c.epoch.sample_conjunction_years_68_width} yr. ` +
      `The same rule applied to JPL's Uranus and Neptune gives ` +
      `${dateOfStr(c.epoch.jpl_uranus_neptune_conjunction_date_tdb)}, ${fmt(days, 1)} days ` +
      `earlier, well inside that spread. This comparison was made after the fact; it was not ` +
      `pre-registered.`,
    cap_plate_ii: `Heliocentric positions on ${dateOfStr(p.date_tdb)}, projected on the ecliptic ` +
      `of date, seen from the north with the equinox to the right. Radii to scale; angles about ` +
      `the line {S@plate-ii}{C@plate-ii} drawn ${p.angle_magnification} times larger than life (the ` +
      `angle label gives the true value). {S@plate-ii} Sun. ` +
      `{U@plate-ii} Uranus (median of the fitted positions). Stipple: ` +
      `${p.posterior.subsample_xy.length} of the ${c.n_samples} posterior samples of the planet ` +
      `(headline posterior). {C@plate-ii} their median direction. {N@plate-ii} Neptune ` +
      `(JPL Horizons). {L@plate-ii} Le Verrier's planet from his published elements (Krajnović ` +
      `2016; Britannica 1911). Dotted: the paths of Uranus and Neptune from JPL Horizons and Le ` +
      `Verrier's ellipse. After-the-fact diagnostic, not pre-registered.`,
  };
}

/* ---------- § 4.3 Fig. 5 ---------- */

export function windowText(D) {
  const w = D.window, m = w.marks;
  const [a, b] = w.in68[0] ?? [NaN, NaN];
  const miss = (k) => deg(Math.abs(m[k].le_verrier_minus_jpl_deg));
  return {
    win_p1: `Hold Le Verrier's published orbit fixed and carry it through time. Every longitude ` +
      `and difference in this section is heliocentric, ecliptic of date. His planet crosses our ` +
      `posterior median once, at ${fmt(w.zero[0], 2)}, and stays inside the 68 % band from ` +
      `${fmt(a, 2)} to ${fmt(b, 2)}. On Galle's night it was ${miss("galle_night")} from ` +
      `Neptune in heliocentric longitude (the ${fmt(Math.abs(D.rep.galle_lon_ecl_of_date_deg.vs.le_verrier_derived.median_minus_target_deg), 2)}° ` +
      `of § 4.1 is geocentric, and from our median). A search on ` +
      `${dateOfJd(m.search_1836.jd_tdb)} would have been off by ${miss("search_1836")}; one on ` +
      `${dateOfJd(m.search_1856.jd_tdb)}, by ${miss("search_1856")}. We did not model the ` +
      `search fields or the star charts.`,
    cap_fig_5: `Heliocentric longitude (ecliptic of date) relative to the posterior median, ` +
      `${yr(w.series.year[0])}–${yr(w.series.year.at(-1))}: the 68 % band (solid edges) and 95 % ` +
      `band (dashed edges) of the headline posterior; {L@fig-5} Le Verrier's planet (heavy ` +
      `line); {N@fig-5} Neptune (JPL Horizons). {C@fig-5} the conjunction of § 4.2; ` +
      `{G@fig-5} Galle's night. After-the-fact diagnostic, not pre-registered.`,
  };
}

/* ---------- § 5 Fig. 6 ---------- */

export const SENS = (fits) => ({
  s1_e05: `Eccentricity up to ${fits.s1_e05.config.e_max}`,
  s2_sigma_data: "Error levels from the data alone",
  s3_t4: `Student-t errors (ν = ${fits.s3_t4.config.nu})`,
});

/** S3 against the headline: the median barely moves, the 68 % range widens (§ 5 and § 6). */
function s3Effect(D, it) {
  const f = D.checks.fits.s3_t4, holds = D.checks.verdicts.s3_t4.holds;
  const w0 = D.checks.reported.width68, w3 = f.galle_lon_ecl_of_date_deg.width68;
  const moves = fmt(Math.abs(f.vs_headline.galle_median_shift_headline_half_widths), 2);
  return `Under Student-t errors the median ${holds ? `is insensitive to ${it}` : `depends on ${it}`} (it moves ` +
    `${moves} headline half-widths), ${holds ? "but" : "and"} the 68 % range ` +
    `${w3 > w0 ? "widens" : "narrows"} from ${deg(w0, 2)} to ${deg(w3, 2)}`;
}

export function robustText(D) {
  const k = D.checks, rp = k.repeat, v = k.verdicts, f = k.fits;
  const ratio = drv(rp.r1Spread / rp.r1RunError, "R1 spread / run error");
  const shift = (key) => fmt(Math.abs(f[key].vs_headline.galle_median_shift_headline_half_widths), 2);
  const jplAll = ["headline", ...Object.keys(SENS(f))].every((key) => f[key].galle_lon_ecl_of_date_deg.vs.jpl.in_68);
  const s1 = v.s1_e05;
  const emax = f.s1_e05.config.e_max;
  let s1Text = `With eccentricity allowed up to ${emax}, P({i:a} < 25 AU) = ${fmt(s1.p, 3)}, which ` +
    `the rule reads as "${s1.word}".`;
  if (s1.pending) {
    s1Text = `With eccentricity allowed up to ${emax}, the published run gives P({i:a} < 25 AU) = ` +
      `${fmt(s1.p, 3)}, which the rule reads as "${s1.word}"; the faster approximate likelihood ` +
      `gives ${fmt(s1.m8p, 3)}, and the verdict is pending a full-likelihood check.`;
  }
  return {
    rob_p1: `The first published 68 % range, ${rng(k.first.lo68, k.first.hi68, 2)} ` +
      `(${fmt(k.first.width68, 2)}° wide), rested on one run of the main cell. Four repeats at ` +
      `the same settings gave widths from ${deg(Math.min(...rp.r1Runs), 2)} to ` +
      `${deg(Math.max(...rp.r1Runs), 2)}, a spread ${fmt(ratio, 0)} times the sampler's own error ` +
      `estimate. Two runs with four times the live points agreed on the width to ` +
      `${fmt(drv(100 * rp.r2RelDiff, "percent"), 1)} % (rule: within ${fmt(drv(100 * rp.tol, "percent"), 0)} %); ` +
      `under the pre-registered rule their merge replaced that cell, and the range became ` +
      `${rng(k.reported.lo68, k.reported.hi68, 2)} (${fmt(k.reported.width68, 2)}° wide).`,
    rob_p2: `Three sensitivity runs, with reading rules committed before they ran, test the ` +
      `model's assumptions. ${s1Text} With error levels estimated from the data alone the ` +
      `result ${v.s2_sigma_data.holds ? "survives" : "does not survive"}: the median moves ` +
      `${shift("s2_sigma_data")} headline half-widths. The evidence favours Student-t errors, ` +
      `ν = ${f.s3_t4.config.nu}, over the Gaussian ones: ln {i:B} = ` +
      `${signed(f.s3_t4.vs_headline.ln_b, 1)} ± ${fmt(f.s3_t4.vs_headline.ln_b_err, 1)}. The ` +
      `Gaussian error model was fixed by the pre-registration and stays the headline's. ` +
      `${s3Effect(D, "that choice")}. ` +
      `${jplAll ? "Neptune's longitude (JPL) is inside the 68 % range in the headline fit and in all three runs." : ""}`,
    cap_fig_6: `The planet's geocentric longitude on Galle's night (ecliptic of date): the ` +
      `headline as first published (superseded) and as reported, and the three sensitivity ` +
      `runs. Dot: median; thick line: 68 % range; thin line: 95 % range, dashed where ` +
      `provisional. {N@fig-6} Neptune (JPL Horizons). Reading rules: ${k.rulesSource}.`,
  };
}

/* ---------- § 6 Limitations ---------- */

/** "a", "a and b", "a, b and c" (no serial comma, as in the README). */
const andList = (xs) => (xs.length < 2 ? xs.join("") : `${xs.slice(0, -1).join(", ")} and ${xs.at(-1)}`);

/** The spec's M6c fits deferred to v1.1 (m7_timeline.json "deferred", from the release checklist). */
function deferredFits(d) {
  const variants = d.m6c_fits.filter((f) => f.kind === "variant").map((f) => f.name);
  const rest = d.m6c_fits.filter((f) => f.kind !== "variant").map((f) => `the ${f.name}`);
  const opens = Object.entries(d.opens).map(([key, what]) =>
    `the ${d.m6c_fits.find((f) => f.key === key).name} variant opens ${what}`);
  return `Only ${d.fitted} is fitted: ${andList([`the ${andList(variants)} variants`, ...rest])} ` +
    `are deferred to v${d.to_version} (${[d.status, ...opens].join("; ")}).`;
}

export function limits(D) {
  const s = D.dataset.spot_checks.human, c = D.dataset.counts, k = D.checks;
  const bad = drv(s.cells - s.matching, "cells in disagreement");
  // Exact one-sided Clopper-Pearson upper bound for 0 errors in n rows: 1 - 0.05^(1/n).
  const ub = drv(100 * (1 - 0.05 ** (1 / s.rows)), "Clopper-Pearson 95 % upper bound, percent");
  const out = [
    `Transcription check: ${bad === 0 ? 0 : "some"} of ${s.rows} rows in error; 95 % upper bound ` +
      `on the per-row error rate ${fmt(ub, 1)} % (exact one-sided Clopper–Pearson, 0 of ` +
      `${s.rows}); the remaining ${drv(c.rows - s.rows, "unchecked rows")} rows are unchecked. ` +
      `The table was transcribed twice from the scan and the two transcriptions reconciled, ` +
      `each disagreement resolved against the page image. ${s.rows} rows were then checked by ` +
      `hand against the scan, with ${bad === 0 ? "no" : bad} discrepancies; the remaining ` +
      `${drv(c.rows - s.rows, "unchecked rows")} rows have had no human check.`,
    `Reduction check against JPL: comparison against per-group precision not made: no ` +
      `precision figure was fixed before the data were seen. Medians reported (Table 1).`,
    `The 95 % range is ${k.reported.i95.state === "final" ? "final" : k.reported.i95.note}. The ` +
      `68 % range is the headline.`,
  ];
  if (k.verdicts.s1_e05.pending) {
    out.push(`The eccentricity test's verdict is pending: the published run and the faster ` +
      `approximate likelihood disagree, and a full-likelihood check decides.`);
  }
  out.push(
    deferredFits(D.timeline.deferred),
    `Error model: the evidence favours Student-t errors (ln B = ` +
      `${signed(D.checks.fits.s3_t4.vs_headline.ln_b, 1)}), but the Gaussian error model was ` +
      `fixed by the pre-registration and is kept. ${s3Effect(D, "it")} (§ 5).`,
    `Diurnal parallax is not modelled (up to ${fmt(D.parallax.max_arcsec, 1)}″).`,
    `The headline plan was first pushed after its run had launched; the logs cannot rule out ` +
      `that one of its result files was opened by hand during the run (§ 3.1).`,
    `Plate II and Fig. 5 are diagnostics made after the fact from the pre-registered posterior; ` +
      `they were not pre-registered.`,
  );
  return out;
}

/* ---------- § 7 Data availability ---------- */

function authorList(names) {
  return names.map((n) => {
    const [family, given] = n.split(", ");
    return `${family}, ${given.split(/\s+/).map((g) => `${g[0]}.`).join(" ")}`;
  }).join(", ");
}

export function availability(D) {
  const d = D.dataset, z = d.zenodo;
  const year = z.publication_date.slice(0, 4);
  const cite = (doi, ver) => `${authorList(d.creators)} (${year}). ${d.title}${ver ? ` (Version ${z.version})` : ""} ` +
    `[Data set]. Zenodo. https://doi.org/${doi}`;
  return {
    avail_p1: `The transcription of Table 77 (${d.counts.rows} rows, with both passes, the ` +
      `resolution log and the misprint table) is on Zenodo under ` +
      `${z.license === "cc-zero" ? "CC0 1.0" : z.license}. Cite the first DOI for the dataset in ` +
      `general and the second for exact reproduction, and cite Le Verrier (1846) for the ` +
      `observations themselves.`,
    cite_concept: cite(z.concept_doi, false),
    cite_version: cite(z.version_doi, true),
    avail_p2: `The code and documentation are at github.com/atishaypahadiya/neptune-fork under ` +
      `the MIT License; the files in the repository's data/historical/ folder, which hold the ` +
      `transcription, are under CC0 1.0. Every number on this page is read from the site's data ` +
      `files, which the repository's scripts regenerate.`,
  };
}

/* ---------- § 8 References ---------- */

export function references(D) {
  const src = D.dataset.citation.source, a = src.authors[0], z = D.dataset.zenodo;
  return D.refs.references.map((r) => {
    if (r.id === "lv1846") {
      return { id: r.id, parts: [`${a["family-names"]}, ${a["given-names"].split(/\s+/).map((g) => `${g[0]}.`).join(" ")} (${src.year}). `,
        { i: src.title }, `. ${src.journal}, ${src.start}–${src.end}. ${src.location.name}: ${src.publisher.name}.`] };
    }
    if (r.id === "dataset") {
      return { id: r.id, parts: [`${authorList(D.dataset.creators)} (${z.publication_date.slice(0, 4)}). `,
        { i: D.dataset.title }, ` (Version ${z.version}) [Data set]. Zenodo. `,
        { href: `https://doi.org/${z.version_doi}`, text: `doi:${z.version_doi}` }] };
    }
    // Each element ends with one full stop: none is added after text that already ends in one
    // ("et al."), and an element that is absent adds nothing.
    const parts = [];
    const stop = () => {
      const k = parts.length - 1;
      if (typeof parts[k] === "string") parts[k] = parts[k].replace(/\s+$/, "").replace(/([^.])$/, "$1.") + " ";
      else parts.push(". ");
    };
    parts.push(`${r.authors}${r.year ? ` (${r.year})` : ""}`); stop();
    if (r.book) { parts.push({ i: r.title }); stop(); parts.push(r.venue); stop(); }
    else {
      if (r.title) { parts.push(r.title); stop(); }
      if (r.venue) {
        parts.push({ i: r.venue });
        if (r.volume) parts.push(`, ${r.volume}${r.issue ? `(${r.issue})` : ""}`);
        if (r.pages) parts.push(`, ${r.pages}`);
        stop();
      }
    }
    if (r.doi) { parts.push({ href: `https://doi.org/${r.doi}`, text: `doi:${r.doi}` }); stop(); }
    if (r.url) { parts.push({ href: r.url, text: r.url_text }); stop(); }
    if (r.note) parts.push(r.note);
    return { id: r.id, parts };
  });
}
