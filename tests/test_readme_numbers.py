"""Every number in the README, traced. Each result is computed here from the JSON it comes from
and must appear in README.md as printed; every other number in the README's prose (dates,
citations, levels, labels) is listed below with where it comes from. A number in the prose
that is neither fails the test."""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text(encoding="utf-8")


def _j(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _cp_upper(n: int, conf: float = 0.95) -> float:
    """Exact one-sided Clopper–Pearson upper bound for 0 successes in n trials."""
    return 1.0 - (1.0 - conf) ** (1.0 / n)


def _dm(deg: float) -> str:
    d = int(deg)
    return f"{d}° {round((deg - d) * 60)}′"


def _minus(s: str) -> str:
    return s.replace("-", "−")


def traced() -> list[tuple[str, str, str]]:
    """(what, the string printed in the README, its source)."""
    from fork.model.bodies import M_EARTH_MSUN
    from fork.model.hidden_mass import LE_VERRIER_1846_MASS_MEMOIR_MSUN as LV_MASS

    rep = _j("site/data/m6c_replacement.json")
    g, mass, a = rep["galle_lon_ecl_of_date_deg"], rep["mass_earth"], rep["a_au"]
    fits = _j("site/data/m6_fits.json")
    tgt, s3 = fits["comparison_targets"], fits["fits"]["s3_t4"]
    m9 = _j("site/data/m9_profile.json")
    m10 = _j("site/data/m10_profile.json")
    g9, d10 = m9["datasets"]["all"], m10["datasets"]
    m6d = _j("site/data/m6d.json")
    v6d = m6d["values"]
    ad = _j("site/data/p1b_adams.json")["values"]
    ch = _j("site/data/p1b_challis.json")["values"]
    fit = _j("docs/notes/P1B_adams_fit.json")
    rp = _j("docs/notes/P1B_adams_reproduction.json")["R3_stage_c"]
    cov = _j("site/data/a6c_coverage.json")
    s1 = _j("site/data/s1_mode.json")["values"]
    arb = _j("site/data/m8_arbiter.json")["merged"]
    n_spot = _j("site/data/m7_dataset.json")["spot_checks"]["human"]["rows"]
    par = _j("site/data/m7_parallax.json")
    tot = _j("site/data/ledger.json")["totals"]["all"]
    nep = ad["neptune_same_model"]
    n_obs = d10["all"]["n_observations"]
    n_draws = v6d["shape38.n_draws"]
    o_clean = [
        abs(ad[k])
        for k in ("o_clean_least_squares", "o_clean_project_on_his_places", "o_clean_project_fit")
    ]
    s7 = fit["targets"]["S7"]["mean_lon_1846_10_01_deg"]
    nep_oct1 = fit["neptune_places"]["1846-10-01"]["helio_true"]
    lv = tgt["le_verrier"]["galle_lon_ecl_of_date_deg_derived"]
    lv_a = m10["settings"]["reference_au"]["le_verrier"]
    nep_a = m10["settings"]["reference_au"]["neptune"]
    return [
        # the lead and the result table
        ("longitude, 68 %", f"{g['median']:.1f}° (68 % interval {g['lo68']:.1f}°–{g['hi68']:.1f}°)",
         "m6c_replacement.galle_lon_ecl_of_date_deg"),
        ("longitude, table", f"**{g['median']:.1f}°** (68 %: {g['lo68']:.1f}°–{g['hi68']:.1f}°)",
         "m6c_replacement.galle_lon_ecl_of_date_deg"),
        ("Neptune", f"Neptune was at {tgt['jpl']['galle_lon_ecl_of_date_deg']:.2f}°",
         "m6_fits.comparison_targets.jpl.galle_lon_ecl_of_date_deg"),
        ("Le Verrier, derived", f"{lv:.2f}°¹",
         "m6_fits.comparison_targets.le_verrier.galle_lon_ecl_of_date_deg_derived"),
        ("mass", f"**{mass['median']:.1f} Earths** ({mass['lo68']:.1f}–{mass['hi68']:.1f})",
         "m6c_replacement.mass_earth"),
        ("Neptune mass", f"| {cov['truth']['mass_earth']:.2f} |", "a6c_coverage.truth.mass_earth"),
        ("Le Verrier's mass", f"{LV_MASS / M_EARTH_MSUN:.1f} Earths (1/{round(1 / LV_MASS)} of the "
         "Sun)", "hidden_mass.LE_VERRIER_1846_MASS_MEMOIR_MSUN (Recherches § 137, p. 234)"),
        ("distance", f"**{a['median']:.1f} AU** ({a['lo68']:.1f}–{a['hi68']:.1f})",
         "m6c_replacement.a_au"),
        ("Neptune a, mean", f"| {nep_a:.2f} AU |", "m10_profile.settings.reference_au.neptune"),
        ("Le Verrier a", f"| {lv_a:.2f} AU |", "m10_profile.settings.reference_au.le_verrier"),
        ("observations", f"{n_obs} observations, {d10['all']['n_measurements']} measured "
         "coordinates", "m10_profile.datasets.all"),
        # the window (footnote 1 and the third finding)
        ("geocentric miss", f"{v6d['window.geocentric_miss_deg']:.2f}° short of Neptune",
         "M6d_window.geocentric_miss_deg, via m6d.values"),
        ("window", f"{m6d['window']['arcmin']}′ ({m6d['window']['deg']}°)", "m6d.window"),
        ("sensitivity window", f"a {v6d['window90.half_width_deg']:.2f}° window",
         "M6d_window.windows.sensitivity.half_width_deg"),
        ("sensitivity counts", f"its counts are {100 * v6d['window90.p_free']:.0f} % and none",
         "M6d_window.windows.sensitivity.p_free, shape38/36.p_shape"),
        # the findings
        ("Gaussian, 36 AU",
         _minus(f"(Δ ln L = {g9['reference_points']['le_verrier']['dlnl']:.2f})"),
         "m9_profile.datasets.all.reference_points.le_verrier.dlnl"),
        ("Student-t evidence", f"ln B = +{s3['vs_headline']['ln_b']:.1f}",
         "m6_fits.fits.s3_t4.vs_headline.ln_b"),
        ("Student-t, 36 AU",
         _minus(f"(Δ ln L = {d10['all']['reference_points']['le_verrier']['dlnl']:.2f})"),
         "m10_profile.datasets.all.reference_points.le_verrier.dlnl"),
        ("36 AU and 30", f"at {lv_a:.0f} AU, and the real Neptune is at {nep_a:.0f}",
         "m10_profile.settings.reference_au"),
        ("before 1781", f"the {d10['post1781']['measurements_removed']} measurements made before",
         "m10_profile.datasets.post1781.measurements_removed"),
        ("Flamsteed", f"Flamsteed's {d10['no_flamsteed']['measurements_removed']} alone",
         "m10_profile.datasets.no_flamsteed.measurements_removed"),
        ("held distances", f"({v6d['shape38.held_a_au']} AU, the Bode's-law distance",
         "m6d.values.shape38.held_a_au"),
        ("held distance 2", f"or {v6d['shape36.held_a_au']} AU from his elements",
         "m6d.values.shape36.held_a_au"),
        ("held fits short", f"{abs(v6d['shape38.median_minus_neptune_deg']):.1f}° and "
         f"{abs(v6d['shape36.median_minus_neptune_deg']):.1f}° short",
         "m6d.values.shape38/36.median_minus_neptune_deg"),
        ("draws", f"none of {n_draws:,} posterior draws land within {m6d['window']['deg']}°",
         "m6d.values.shape38.n_draws, p_shape; m6d.window"),
        ("free share", f"{100 * v6d['shape38.p_free']:.0f} % do", "m6d.values.shape38.p_free"),
        ("Adams, rounded", f"by {min(o_clean):.0f}° to {max(o_clean):.0f}°",
         "p1b_adams.values.o_clean_*, rounded"),
        ("Adams, noise-free", f"{min(o_clean):.1f}° to {max(o_clean):.1f}° short",
         "p1b_adams.values.o_clean_*"),
        ("Adams, final solutions", f"{ad['H1_printed_true_lon'] - nep:.1f}° and "
         f"{ad['H2_printed_true_lon'] - nep:.1f}° from Neptune",
         "p1b_adams.values.H1/H2_printed_true_lon - neptune_same_model"),
        # how it works
        ("coverage", f"{cov['a6c_coverage']['n_truth_in_68']} of "
         f"{cov['a6c_coverage']['n_seeds']} noise realisations", "a6c_coverage.a6c_coverage"),
        ("coverage p", f"p = {cov['a6c_coverage']['binomial']['p_value_two_sided']:.2f}",
         "a6c_coverage.a6c_coverage.binomial.p_value_two_sided"),
        ("scale", f"{tot['cpu_hours']:,.0f} CPU-hours in all", "ledger.totals.all.cpu_hours"),
        # the British search
        ("reproduction", "to {:.3f}° and {:.3f}°".format(
            *(abs(rp[h]["vs_print"]["rows"]["theta"]["delta"]) for h in ("H1/C0", "H2/C0"))),
         "P1B_adams_reproduction.R3_stage_c.H1/C0, H2/C0.vs_print.rows.theta.delta"),
        ("fits at his distances", f"{abs(ad['H1_fit_true_lon_median'] - nep):.1f}° and "
         f"{abs(ad['H2_fit_true_lon_median'] - nep):.1f}° short",
         "p1b_adams.values.H1/H2_fit_true_lon_median - neptune_same_model"),
        ("noise-free methods", "his least squares is {:.1f}° short, the model on his places "
         "{:.1f}° and the fit {:.1f}°".format(*o_clean), "p1b_adams.values.o_clean_*"),
        ("post hoc", f"gives {ad['H1_post_hoc_P_a']:.3f}", "p1b_adams.values.H1_post_hoc_P_a"),
        ("315° 20′", _dm(s7), "P1B_adams_fit.targets.S7.mean_lon_1846_10_01_deg"),
        ("12° reconciled", f"It lies {nep_oct1 - s7:.1f}° from Neptune",
         "P1B_adams_fit.neptune_places.1846-10-01.helio_true - S7"),
        ("Challis draws", f"({ch['n_accepted']:,} kept draws)", "p1b_challis.values.n_accepted"),
        ("Challis", f"{ch['C1a_mean']:.2f} ({ch['C1a_q05']:.2f} to {ch['C1a_q95']:.2f})",
         "p1b_challis.values.C1a_mean, q05, q95"),
        ("Challis zones", f"give {ch['C3_S4_mean']:.2f} and {ch['C4a_mean']:.2f}, or "
         f"{ch['C3_S4_prompt_mean']:.2f} and {ch['C4a_prompt_mean']:.2f}",
         "p1b_challis.values.C3_S4(_prompt)_mean, C4a(_prompt)_mean"),
        # limitations
        ("Student-t shift",
         f"moves {abs(s3['vs_headline']['galle_median_shift_headline_half_widths']):.2f} "
         "headline half-widths",
         "m6_fits.fits.s3_t4.vs_headline.galle_median_shift_headline_half_widths"),
        ("68 % widths",
         f"from {g['width68']:.2f}° to {s3['galle_lon_ecl_of_date_deg']['width68']:.2f}°",
         "m6c_replacement.width68; m6_fits.fits.s3_t4.galle_lon_ecl_of_date_deg.width68"),
        ("95 % interval", f"({s1['lo95']:.1f}°–{s1['hi95']:.1f}°)", "s1_mode.values.lo95, hi95"),
        ("low-distance solution", f"near {arb['median']:.0f}°", "m8_arbiter.merged.median"),
        ("its weight", "between 0.001 and about 0.2", "s1_mode.values.why"),
        ("spot check", f"0 errors in {n_spot} rows", "m7_dataset.spot_checks.human.rows"),
        ("spot-check bound", f"below {100 * _cp_upper(n_spot):.1f} %",
         "Clopper–Pearson on 0 of m7_dataset.spot_checks.human.rows"),
        ("unchecked rows", f"The other {n_obs - n_spot} rows", "279 - spot-checked rows"),
        ("parallax", f"up to {par['max_arcsec']:.1f}″", "m7_parallax.max_arcsec"),
    ]  # fmt: skip


def test_every_traced_number_is_printed():
    missing = [(w, s, src) for w, s, src in traced() if s not in README]
    assert not missing, missing


def test_readings_named_in_the_readme_match_the_json():
    m10 = _j("site/data/m10_profile.json")["datasets"]
    assert _j("site/data/m9_profile.json")["datasets"]["all"]["reading"]["reading"] == "constrained"
    assert m10["all"]["reading"]["reading"] == "partial"
    assert m10["post1781"]["reading"]["reading"] == "flat"
    assert m10["no_flamsteed"]["reading"]["reading"] == "flat"
    v = _j("site/data/m6d.json")["values"]
    assert v["shape38.p_shape"] == v["shape36.p_shape"] == 0.0
    assert v["window90.shape38.p_shape"] == v["window90.shape36.p_shape"] == 0.0
    assert v["window.post_hoc"] is True
    assert "0.001 to about 0.2" in _j("site/data/s1_mode.json")["values"]["why"]
    s7 = _j("docs/notes/P1B_adams_fit.json")["targets"]["S7"]
    assert s7["class"] == "extrapolation, not a solution"


def test_le_verrier_mass_on_the_site_is_the_readmes():
    from fork.model.bodies import M_EARTH_MSUN
    from fork.model.hidden_mass import LE_VERRIER_1846_MASS_MEMOIR_MSUN as LV_MASS

    mark = _j("site/data/m9_profile.json")["corner"]["marks"]["le_verrier"]
    assert mark["mass_earth"] == round(LV_MASS / M_EARTH_MSUN, 1)
    assert "§ 137, p. 234" in mark["mass_source"] and "§ 137, p. 234" in README


# Numbers in the README that are not results, each with where it comes from.
NOT_RESULTS = {
    # dates and spans of the record
    "1846": "the prediction and Galle's night", "1847": "1 January 1847, memoir p. 253",
    "1690": "Table 77's first year", "1845": "Table 77's last year",
    "1781": "Uranus recognised (Herschel)", "23": "23 September 1846, Galle's night",
    "31": "31 August 1846, Le Verrier's elements", "1": "1 January 1847",
    "2": "2 September 1846, Adams's place; 2 × 10⁻⁹″ (REBOUND, 1.65e-9″, docs/notes/M3.md)",
    "12": "12 August (Challis); the '12°' of Rawlins",
    "30": "30 July (Challis); Neptune at 30 AU (traced)",
    # citations and labels
    "77": "Table 77", "253": "memoir page", "137": "memoir section", "234": "memoir page",
    "1999": "Rawlins, DIO", "2026": "Whitmore", "1990": "Lai, Lam & Young", "2018": "KTH thesis",
    "4": "RK4, the integrator's name", "68": "68 % level", "95": "95 % level",
    "0": "0 by construction (Challis); 0 errors (traced)",
    "0.5": "eccentricity bound of the sensitivity run (docs/PHASE1_SPEC.md); parallax (traced)",
    "10": "× 10⁻⁹ (docs/notes/M3.md); 8 of 10 (traced)",
    # DOIs and versions, as in CITATION.cff
    "10.5281": "Zenodo's DOI prefix", "24044791": "dataset, all versions",
    "24158992": "dataset v1.0.1", "24121059": "code, all versions", "24357494": "code v1.5.1",
    "1.0.1": "dataset version", "1.5.1": "code version",
}  # fmt: skip


def _prose_numbers(text: str) -> set[str]:
    text = re.sub(r"```.*?```", " ", text, flags=re.DOTALL)
    text = re.sub(r"`[^`\n]*`", " ", text)
    text = re.sub(r"\]\([^)\s]*\)", "]", text)
    text = re.sub(r"https?://\S+", " ", text)
    return set(re.findall(r"(?<![\w.])\d+(?:,\d{3})*(?:\.\d+)*", text))


def test_every_number_in_the_prose_is_traced_or_listed():
    covered: set[str] = set(NOT_RESULTS)
    for _, s, _ in traced():
        covered |= _prose_numbers(s)
    stray = sorted(_prose_numbers(README) - covered)
    assert not stray, f"numbers in README.md with no source here: {stray}"
