// The companion page's text (site/adams-challis/), in the voice of a paper (docs/DESIGN.md,
// Voice and "Companion page"). Every number is read from site/data/p1b_*.json at load time; a
// number computed from JSON values goes through drv(), which records it for the number audit.
// `{X@fig}` in a template is the point letter X linked to its point in that figure.

import { fmt, drv } from "./data.js";

const MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August",
  "September", "October", "November", "December"];

/** "12 August" (or "12 August 1846") from "1846-08-12". */
export function dayOf(iso, { year = false } = {}) {
  const [y, m, d] = iso.slice(0, 10).split("-").map(Number);
  return `${d} ${MONTHS[m - 1]}${year ? ` ${y}` : ""}`;
}
/** 100,000: digit groups of a count (each group is registered as derived). */
export function count(n) {
  const s = String(n), groups = [];
  for (let i = s.length; i > 0; i -= 3) groups.unshift(s.slice(Math.max(0, i - 3), i));
  groups.forEach((g) => drv(Number(g), "digit group of a count"));
  return groups.join(",");
}
const signed = (x, n = 1) => `${x < 0 ? "−" : "+"}${Math.abs(x).toFixed(n)}`;
const abs = (x, n = 1, note = "magnitude") => fmt(drv(Math.abs(x), note), n);
const p2 = (x) => fmt(x, 2);
const rng = (lo, hi) => `${p2(lo)} to ${p2(hi)}`;
const INPUT = {
  p_follow: "the chance that he treats an unmatched star as a candidate",
  p_disc_routine: "the chance that a disc is noted in routine recording",
  n_unknown_sessions: "the number of his sessions after 12 August",
  c4b_breadth_arcmin: "the breadth of each sweep in declination",
};
export const CASE = {
  A0: "His actual procedure",
  C1a: "The comparison of 12 August with 30 July, carried through",
  C1b: "The 4 August sweep compared with 12 August",
  C2_1800: "The Berlin sheet he had, limits for the equinox of 1800",
  C2_date: "The same sheet, limits applied to right ascension of date",
  C3_S4: "Zones laid from Adams's first place (S4)",
  C3_S4_prompt: "From S4, with prompt comparison",
  C3_S5: "Zones laid from Adams's second place (S5)",
  C3_S5_prompt: "From S5, with prompt comparison",
  C4a: "Zones laid from Le Verrier's June place",
  C4a_prompt: "From Le Verrier's June place, with prompt comparison",
  C4b: "The 29 September procedure, begun on 29 July",
};
// What each case is conditional on, in the page's words; the note's wording is in
// p1b_challis.json (`*_conditional_on`).
const COND = {
  record: "the printed record and the declared priors",
  s1800: "the sheet's limits read for the equinox of 1800",
  sdate: "the sheet's limits read as right ascension of date",
  zones: "zones alternating outward from the start place; session dates from the priors",
};
const COND_OF = { A0: "record", C1a: "record", C1b: "record", C2_1800: "s1800", C2_date: "sdate" };
// The same cases in fewer words, for Fig. 4 on phones.
const CASE_SHORT = {
  A0: "His actual procedure",
  C1a: "12 August compared with 30 July",
  C1b: "4 August compared with 12 August",
  C2_1800: "Berlin sheet, equinox of 1800",
  C2_date: "Berlin sheet, right ascension of date",
  C3_S4: "From Adams's first place",
  C3_S4_prompt: "First place, prompt comparison",
  C3_S5: "From Adams's second place",
  C3_S5_prompt: "Second place, prompt comparison",
  C4a: "From Le Verrier's June place",
  C4a_prompt: "June place, prompt comparison",
  C4b: "29 September procedure from 29 July",
};
const ZERO_WHY = {
  A0: "by construction (the realised events are held fixed)",
  C2_1800: "the planet lies off the sheet",
  C3_S5: "by timing: the zones do not reach the planet",
  C3_S5_prompt: "by timing: the zones do not reach the planet",
};
export const galle = (K) => `${K.galle_night.day} ${K.galle_night.month}`;

/** O = a longitude minus Neptune's (registered as derived). */
const off = (x, ref, note) => drv(x - ref, note);

export function facts(D) {
  const A = D.adams.values, W = D.challis.values, K = D.challis.constants;
  const printed = (t) => off(A[`${t}_helio_true_galle`], A.neptune_helio_true_galle,
    `${t} minus Neptune, Galle's night`);
  const fit = (h) => ({
    med: off(A[`${h}_fit_true_lon_median`], A.neptune_same_model, `${h} fit minus Neptune`),
    lo: off(A[`${h}_fit_true_lon_lo68`], A.neptune_same_model, `${h} fit minus Neptune`),
    hi: off(A[`${h}_fit_true_lon_hi68`], A.neptune_same_model, `${h} fit minus Neptune`),
  });
  const clean = [A.o_clean_least_squares, A.o_clean_his_procedure,
    A.o_clean_project_on_his_places, A.o_clean_project_fit];
  return { A, W, K, printed, fit, cleanMin: Math.min(...clean.map(Math.abs)),
    cleanMax: Math.max(...clean.map(Math.abs)) };
}

export function titleBlock(D) {
  return {
    title: D.page.title,
    subtitle: "A companion study to ",
    main_title: D.page.main_title,
    author: D.dataset.creators.map((c) => c.split(", ").reverse().join(" ")).join(", "),
    status: `Phase 1b of the project: ${D.page.status}.`,
  };
}

export function abstract(D) {
  const { A, W, K, printed, fit, cleanMin, cleanMax } = facts(D);
  const R = K.range_percent;
  // The share of the draws that flag the planet in which the comparison is done by 16 August.
  const by = D.challis.sources.c1a_p_by_0816.at(-1);
  const share = drv(Math.round((100 * W.c1a_p_by_0816) / W.C1a_mean), "percent of flags by a date");
  return {
    abstract: `We put two questions to the British search for Neptune of ` +
      `${A.S4_date.slice(0, 4)} and ${K.galle_night.year}, using the printed record alone. How ` +
      `good were Adams's solutions? And could Challis's sweeps at Cambridge have recognised ` +
      `the planet before Galle found it on ${galle(K)} ${K.galle_night.year}?`,
    abstract_p2: `Adams's two final solutions place the planet ${abs(printed("S5"))}° and ` +
      `${abs(printed("S6"))}° from Neptune on that night. We reproduce both from his memoir to ` +
      `within ${fmt(Math.max(A.H1_c0_theta_delta, A.H2_c0_theta_delta), 2)}°. The project's ` +
      `model, fitted to comparable data at his assumed distance, puts the planet ` +
      `${abs(fit("H1").med)}° short of Neptune, and every method we tried on a noise-free ` +
      `Neptune at that distance is ${fmt(drv(cleanMin, "smallest clean shortfall"))}° to ` +
      `${fmt(drv(cleanMax, "largest clean shortfall"))}° short. The shortfall belongs to the ` +
      `assumed distance; his nearness belongs to his data or their reduction. Whether noise ` +
      `alone could account for it we could not settle: the pre-registered statistic measured ` +
      `something else, and the post hoc one reads weak.`,
    abstract_p3: `Under the reading we report, Challis's actual procedure has probability ` +
      `${p2(W.A0_mean)} of recognising the planet before ${galle(K)}. That is by ` +
      `construction and is not a finding: the realised events are held fixed, so the cases ` +
      `measure what a change would have added. Finishing the comparison of ` +
      `${dayOf(K.events.recorded_second)} with ${dayOf(K.events.zone_first_swept)} that he ` +
      `broke off (complete by ${dayOf(by)} in ${share} % of the draws that flag the planet) ` +
      `gives ${p2(W.C1a_mean)} (${R.low} to ${R.high} % range ` +
      `${rng(W.C1a_q05, W.C1a_q95)}) before ${galle(K)}. That figure rests on the chance ` +
      `that he treats an unmatched star as a candidate, which the record does not inform: ` +
      `its prior mean of ${p2(W.p_follow_prior_mean)} is ${p2(W.p_follow_posterior_mean)} ` +
      `after the conditioning. Zones laid from Adams's first place or from Le ` +
      `Verrier's give ${p2(W.C3_S4_mean)} and ${p2(W.C4a_mean)}, or ` +
      `${p2(W.C3_S4_prompt_mean)} and ${p2(W.C4a_prompt_mean)} with prompt comparison. Each ` +
      `depends on an input the record does not fix, and none reaches the pre-registered ` +
      `reading "very probably".`,
  };
}

export function sourcesText(D) {
  const { A } = facts(D);
  const passes = A.transcription_passes.split(", all decided")[0];
  const d = A.fit_data_set;
  return {
    src_p1: `Everything here rests on print: Adams's memoir, Challis's two accounts of his ` +
      `search, and Le Verrier's papers of June and August ${D.challis.constants.galle_night.year} ` +
      `(§ 7). The Cambridge observing books and Airy's papers were not consulted; § 5 says ` +
      `where that matters.`,
    src_p2: `The memoir's tables were transcribed in ${passes}, each decided on the page ` +
      `image. Challis's search is tabulated night by night from his own text, each row with ` +
      `its page. For the fit of § 3.1 we took, from the observations of the main paper, the ` +
      `groups comparable to what Adams used: ${d.total[0]} observations ` +
      `(${d.total[1]} measurements), of which ${d.ancient[0]} are the early chance ` +
      `observations, ${d.greenwich[0]} are from Greenwich and ` +
      `${d.paris_1835_1845[0]} from Paris.`,
    src_p3: `Each step was pre-registered before it ran. The plans, the results and every ` +
      `deviation are in the notes listed in § 6; where a deviation bears on a claim, the ` +
      `sentence that makes the claim says so.`,
  };
}

export function adamsText(D) {
  const { A, K, printed, fit } = facts(D);
  const f1 = fit("H1"), f2 = fit("H2"), F = K.fit_ranges_percent, T = K.noise_thresholds;
  const P = (x) => fmt(x, 3);
  return {
    a31_p1: `Adams did not solve for the planet's distance; he assumed it. His first ` +
      `hypothesis put the mean distance at ${fmt(A.S5_a_au)} AU and his second at ` +
      `${fmt(A.S6_a_au, 2)}; Neptune was then ${fmt(A.truth_r_au)} AU from the Sun. We fitted ` +
      `the main paper's model to the comparable data with the distance held at each value.`,
    a31_p2: `On Galle's night the posterior of the planet's true heliocentric longitude lies ` +
      `${abs(f1.med)}° short of Neptune under the first hypothesis (${F.inner} % range ` +
      `${signed(f1.lo)}° to ${signed(f1.hi)}°) and ${abs(f2.med)}° short under the second ` +
      `(${signed(f2.lo)}° to ${signed(f2.hi)}°); Neptune is ${A.H1_fit_status} the ` +
      `${F.outer} % range of the geocentric longitude in both. With the distance left free ` +
      `the same data put Neptune ${A.free_distance_status} the ${F.inner} % range. His ` +
      `printed solutions are much nearer than our fits at his own distances: ` +
      `${signed(printed("S5"))}° and ${signed(printed("S6"))}° (Fig. 1, {P@fig-1} and ` +
      `{Q@fig-1} against {F@fig-1} and {G@fig-1}).`,
    a32_p1: `To see whether that nearness is in his method, we rebuilt his calculation from ` +
      `the memoir: the series, the equations of condition, the normal equations and his ` +
      `elimination. From his printed normal equations our solution agrees with his in the ` +
      `planet's longitude at epoch to ${fmt(A.H1_c0_theta_delta, 3)}° and ` +
      `${fmt(A.H2_c0_theta_delta, 3)}° for the two hypotheses.`,
    a32_p2: `From his tabulated data, the same method gives a mean longitude ` +
      `${fmt(A.H1_c2_vs_print, 2)}° and ${fmt(A.H2_c2_vs_print, 2)}° away from his print, and ` +
      `within ${fmt(A.H1_c2_vs_brookes, 2)}° and ${fmt(A.H2_c2_vs_brookes, 2)}° of Brookes's ` +
      `re-derivation. The difference sits in one normal equation, whose coefficients he ` +
      `rounded to four decimals: one unit in the fourth decimal of one coefficient moves the ` +
      `longitude by ${abs(A.H1_dtheta_per_step)}° under the first hypothesis and ` +
      `${abs(A.H2_dtheta_per_step)}° under the second. His problem is that sensitive to its ` +
      `own arithmetic.`,
    a33_p1: `We then gave every method a Neptune without noise: Uranus's motion computed with ` +
      `the real Neptune, sampled at his epochs, and solved at his first distance. His least ` +
      `squares lands ${abs(A.o_clean_least_squares)}° short, his own elimination ` +
      `${abs(A.o_clean_his_procedure)}°, the project's model on his normal places ` +
      `${abs(A.o_clean_project_on_his_places)}°, and the fit of § 3.1 on synthetic data ` +
      `${abs(A.o_clean_project_fit)}° (Fig. 1, open marks). The fit of § 3.1 on the real ` +
      `data is ${abs(A.o_real_project_fit)}° short. So the shortfall is a property of the ` +
      `assumed distance, and his data through his method (${signed(A.H1_o_real_least_squares)}°, ` +
      `{R@fig-1}) are nearer Neptune than a noise-free Neptune allows.`,
    a34_p1: `Is that nearness luck? We added independent noise to the noise-free places, ` +
      `${count(A.noise_n_draws)} times per noise model, and asked how often the solution moves ` +
      `toward Neptune by as much as his did (${fmt(A.H1_delta_obs)}° under the first ` +
      `hypothesis, {D@fig-2}). The pre-registered statistic gives ${P(A.H1_prereg_P_a)} and ` +
      `${P(A.H1_prereg_P_b)} under the two noise models, which the pre-registered rule reads ` +
      `as "${A.H1_prereg_verdict}" (above ${T.suffices}).`,
    a34_p2: `That reading does not mean what it says, and we report it as a deviation. The ` +
      `statistic treated the shift as a number on a line, and it is an angle. His ` +
      `least-squares problem has a second minimum, at an eccentricity of ` +
      `${fmt(A.H1_second_minimum_e_prime, 2)} and ${fmt(A.H1_second_minimum_O, 0)}° from ` +
      `Neptune, only ${fmt(A.H1_second_minimum_chi2_minus_best)} in χ² above the best. Noise ` +
      `sends some draws there, and ${fmt(drv(100 * A.H1_a_exceeding_branch_changes, "percent"), 0)} % ` +
      `of the draws counted as exceeding his shift are such jumps (Fig. 2, the clusters at ` +
      `the edges). The pre-registered number measures how often noise switches minimum.`,
    a34_p3: `Post hoc, and labelled so wherever it appears: among the draws that stay on the ` +
      `minimum of the real data (${A.post_hoc_definition.split("among ")[1]}), a shift as ` +
      `large as his occurs with probability ${P(A.H1_post_hoc_P_a)} and ` +
      `${P(A.H1_post_hoc_P_b)} under the first hypothesis and ${P(A.H2_post_hoc_P_a)} and ` +
      `${P(A.H2_post_hoc_P_b)} under the second. By the same thresholds that is ` +
      `"${A.H1_post_hoc_verdict}": below ${T.suffices}, above ${T.systematic}. Neither luck ` +
      `nor a systematic cause is established. The tables of Bouvard that his reductions ` +
      `rested on are the route to an attribution, and we did not take it.`,
  };
}

export function challisText(D) {
  const { W, K } = facts(D);
  const R = K.range_percent, T = K.rule_thresholds, E = K.events;
  const recs = D.challis.recordings.filter((r) => r.counted === "yes");
  const nRec = drv(recs.reduce((s, r) => s + r.nights, 0), "pre-discovery recordings counted");
  const dep = D.challis.case_order.filter((c) => W[`${c}_prior_dependent`]);
  const names = (cs) => cs.map((c) => CASE[c].charAt(0).toLowerCase() + CASE[c].slice(1));
  const list = (xs) => (xs.length > 1 ? `${xs.slice(0, -1).join("; ")}; and ${xs.at(-1)}` : xs[0]);
  return {
    c41_p1: `Challis began on ${dayOf(E.search_start)}, sweeping zones ` +
      `${K.zone_breadth_arcmin.value}′ broad with the Northumberland refractor ` +
      `(${K.aperture.text} in aperture) at a power of ${K.power.value}, around the place ` +
      `Adams's work gave him ({C@fig-3}). He recorded the planet on ${dayOf(E.recorded_first)} ` +
      `and again on ${dayOf(E.recorded_second)}, as star No. ${K.planet_number_12_august.value} of ` +
      `that night's series, and did not know it: his comparison of that series with the ` +
      `earlier one stopped at No. ${K.comparison_stopped_at.value}. On ` +
      `${dayOf(E.disc_noted)}, looking for a disc after reading Le Verrier, he noted one.`,
    c41_p2: `We model each session as three events: the planet's place lies in the zone, it ` +
      `is recorded, and it is recognised, by a noted disc, by a comparison of two series, or ` +
      `against a chart. What the record does not print (the zone's setting, the dates and ` +
      `lengths of the later sessions, how readily a disc is noted) is drawn from declared ` +
      `priors, and a draw is kept only if his actual procedure then reproduces the record. ` +
      `Of ${count(W.n_drawn)} draws, ${count(W.n_accepted)} are kept ` +
      `(${fmt(drv(100 * W.acceptance, "percent"), 0)} %). The record lowers the chance of ` +
      `a disc noted in routine recording from a prior mean of ` +
      `${p2(W.p_disc_routine_prior_mean)} to ${p2(W.p_disc_routine_posterior_mean)}; the ` +
      `${nRec} earlier recordings of Neptune by other observers, none with a disc noted, are ` +
      `why that prior is wide (Kowal and Drake; Mauvais; Hind; Schönfeld, § 7).`,
    c41_p3: `Neptune's magnitude (${fmt(W.neptune_0729_V, 2)} on ${dayOf(E.search_start)}) is ` +
      `carried back from a modern fit (Mallama and Hilton) and its disc was ` +
      `${fmt(W.neptune_0729_diameter_arcsec, 2)}″ across. The places of Adams and Le Verrier ` +
      `are set in the ecliptic, as both assumed, which is where most of their offset from ` +
      `the planet in declination comes from (Fig. 3).`,
    c42_p1: `Table 1 gives, for each case, the probability that the planet is flagged ` +
      `before ${galle(K)}, as the mean and the ${R.low} to ${R.high} % range over the kept ` +
      `draws. The range is over what the record leaves open; it does not shrink with more ` +
      `draws. The pre-registered rule reads a case as "very probably" when the whole range ` +
      `is above ${T.high}, "very probably not" when it is below ${T.low}, and otherwise ` +
      `names the input it depends on (Fig. 4).`,
    c42_p2: `These are the reading with the realised events held: a recording that a case ` +
      `shares with the history keeps what happened (no disc noted in August), and only the ` +
      `events the history does not contain are drawn. Under it his actual procedure is ` +
      `${p2(W.A0_mean)} by construction: the realised events are held fixed, so the cases measure what a change would have added. The pre-registered rule for paired differences ` +
      `against his actual procedure is therefore uninformative, and we do not quote its ` +
      `statements.`,
    c42_p3: `Three readings change when every prior interval is halved, and are marked ` +
      `prior-dependent (†): ${list(names(dep))}. The priors on the zone's setting were not ` +
      `halved; that check used ${count(W.halved_n_accepted)} kept draws and the ranking of ` +
      `inputs the first ${count(W.tornado_n_draws)}.`,
    c42_p4: `Outside the pre-registered reading: if every event is drawn afresh, his ` +
      `procedure taken as a plan, the actual procedure's probability is ` +
      `${p2(W.replay_a0)} (cases.A0.replay.primary.mean in P1B_challis_detection.json); the ` +
      `note carries the rest. Paired differences between the cases are post hoc and are in ` +
      `the note only.`,
    c43_p1: `The cheapest change is the one he began. Carrying the ` +
      `${dayOf(E.recorded_second)} comparison through needs no new observation and a few days at ` +
      `the desk (done by ${dayOf(W.c1a_done_by)} in ${R.high} % of the draws that flag it), and its probability ` +
      `is the chance that he follows up an unmatched star, times the chance that the earlier ` +
      `zone reached the planet's place.`,
    c43_p2: `The place he started from matters less than the method. From Adams's first ` +
      `place or Le Verrier's the planet is recorded within weeks, but a routine recording ` +
      `is recognised only through a noted disc, and the evidence on that points low. From ` +
      `Adams's second place the zones do not reach the planet in time ({B@fig-3}). The ` +
      `chart he had helps only if its limits are read one way: the case ` +
      `${W.c2_combined}.`,
    table_note: `† Prior-dependent: the reading changes row when every prior interval is ` +
      `halved. A delay of a week before either comparison leaves its probability unchanged.`,
  };
}

/** Table 1's rows. */
export function caseRows(D) {
  const { W } = facts(D);
  return D.challis.case_order.map((c) => {
    const zero = W[`${c}_q95`] === 0;
    const top = W[`${c}_top_input`];
    const row = W[`${c}_row`];
    if (zero && !row.startsWith("would very probably not")) throw new Error(`row of ${c}: ${row}`);
    return {
      key: c, name: CASE[c], short: CASE_SHORT[c], dagger: W[`${c}_prior_dependent`],
      p: zero ? p2(W[`${c}_mean`]) : `${p2(W[`${c}_mean`])} (${rng(W[`${c}_q05`], W[`${c}_q95`])})`,
      mean: W[`${c}_mean`], lo: W[`${c}_q05`], hi: W[`${c}_q95`],
      hlo: W[`${c}_halved_q05`], hhi: W[`${c}_halved_q95`],
      reading: zero ? `very probably not: ${ZERO_WHY[c]}`
        : row === "depends on" ? `depends on ${INPUT[top] ?? top}` : row,
      cond: COND[COND_OF[c] ?? "zones"],
      zeroWhy: zero ? ZERO_WHY[c] : "",
    };
  });
}

export function captions(D) {
  const { A, W, K } = facts(D);
  const F = K.fit_ranges_percent, R = K.range_percent, E = K.events;
  return {
    cap_fig_1: `True heliocentric longitude on Galle's night minus Neptune's ({N@fig-1}), ` +
      `per solution. {P@fig-1}, {Q@fig-1}: Adams's printed solutions under his first and ` +
      `second hypotheses. {R@fig-1}: his data through his least squares. {F@fig-1}, ` +
      `{G@fig-1}: the project's fit with the distance held at ${fmt(A.S5_a_au)} and ` +
      `${fmt(A.S6_a_au, 2)} AU, with the ${F.inner} % range. Open marks: each method on a ` +
      `noise-free Neptune at the first distance; only the last has a range. Rows without a ` +
      `line are single solutions.`,
    cap_fig_2: `The shift of the least-squares solution toward Neptune over ` +
      `${count(A.noise_n_draws)} noise draws, on the circle, in bins of ` +
      `${A.hist_bin_deg}°. Heavy line: noise scaled to his own residuals; light line: the ` +
      `second noise model. {D@fig-2}: the shift of his real data. The clusters near the ` +
      `edges are draws on the second minimum. The on-branch probabilities quoted in § 3.4 ` +
      `are post hoc; the pre-registered values are ${fmt(A.H1_prereg_P_a, 3)} and ` +
      `${fmt(A.H2_prereg_P_a, 3)} for the heavy lines.`,
    cap_fig_3: `The search on the sky, to scale; right ascension increases to the left. ` +
      `Dotted: Neptune from ${dayOf(E.search_start)} to ${dayOf(E.track_end)}, with ` +
      `{N@fig-3} its place on ${dayOf(D.challis.map_night)}. The band is the ` +
      `${K.zone_breadth_arcmin.value}′ zone of ${dayOf(E.zone_first_swept)} and ` +
      `${dayOf(E.recorded_second)} at the middle of the settings the record allows, between the ` +
      `two hairlines. Open marks, on ${dayOf(D.challis.map_night)}: {C@fig-3} the place ` +
      `Adams's ephemeris gave Challis, {A@fig-3} and {B@fig-3} Adams's first and second ` +
      `places, {J@fig-3} and {L@fig-3} Le Verrier's of June and August. Vertical lines: the ` +
      `western edge of the Berlin sheet he had, with its limits read for the equinox of ` +
      `${K.sheet_west_edge_ra.equinox} and as right ascension of date.`,
    cap_fig_4: `Each case's probability of recognition before ${galle(K)}: the mean (dot) ` +
      `and the ${R.low} to ${R.high} % range (line), against the rule's thresholds at ` +
      `${K.rule_thresholds.low} and ${K.rule_thresholds.high}. † rows are prior-dependent; ` +
      `the thin line under each is its range with the priors halved. A dot at ` +
      `${p2(W.A0_mean)} with no line is a case that is ${p2(W.A0_mean)} in every draw: by ` +
      `construction for his actual procedure (the realised events are held fixed), off the sheet for the chart under the first ` +
      `limits, by timing for Adams's second place.`,
    cap_table_1: `Table 1. The cases: probability that the planet is flagged before ` +
      `${galle(K)} (mean, with the ${R.low} to ${R.high} % range), the row of the ` +
      `pre-registered reading rule, and what the case is conditional on.`,
  };
}

export function limits(D) {
  const { W } = facts(D);
  return [
    "The printed record stands in for the Cambridge observing books. Where it is silent the " +
      "model uses declared priors, and every probability in § 4 is conditional on them.",
    "Under the reading we report, his actual procedure cannot succeed before Galle's night by " +
      "construction: the realised events are held fixed, so the cases measure what a change would have added.",
    "The limits of the Berlin sheet he had are taken from a secondary source, and the chart " +
      `case ${W.c2_combined}.`,
    "The noise test of § 3.4 was read post hoc after its pre-registered statistic proved " +
      "unsuitable for an angle; its verdict is weak either way.",
    "The fit of § 3.1 uses the project's data and error model, not Adams's own reductions; " +
      "Bouvard's tables, on which those rested, were not transcribed.",
    "Neptune's magnitude in 1846 is extrapolated from photometry of the twentieth century.",
    "We do not say what Challis ought to have done, and the model is not a probability for " +
      "the British search as a whole: it covers his nights and his instrument only.",
  ];
}

export function availability(D) {
  const repo = D.release.repository;
  return {
    avail_p1: `The study is ${D.page.status}. Its plans, results and deviations are in the ` +
      `repository under docs/notes/: ${D.page.notes.join(", ")}. The numbers on this page are ` +
      `read from site/data/p1b_adams.json and p1b_challis.json, which copy named fields of ` +
      `the notes' JSON files.`,
    avail_p2: "The observations of Uranus and the code are those of the main paper, under " +
      "its licences and DOIs.",
    repo_url: repo,
  };
}

export function references(D) {
  return D.page.references.map((r) => ({
    id: r.id,
    parts: r.parts.map((s) => (s.startsWith("{") && s.endsWith("}") ? { i: s.slice(1, -1) } : s)),
  }));
}
