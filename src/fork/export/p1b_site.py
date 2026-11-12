"""The companion page's data (site/adams-challis/): subsets of the committed Phase 1b JSONs.

Every entry of ``values`` is copied from one field of one committed file, named in ``sources``
as ``[file, key, key, ...]``; nothing is computed here, so the page's numbers trace to the
notes' JSONs (docs/DESIGN.md, the one rule). ``tests/test_p1b_site.py`` resolves every source.

    uv run python -m fork.export.p1b_site
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NOTES = "docs/notes/"
FIT = NOTES + "P1B_adams_fit.json"
REPRO = NOTES + "P1B_adams_reproduction.json"
NOISE = NOTES + "P1B_adams_noise.json"
HIST = NOTES + "P1B_adams_noise_hist.json"
DET = NOTES + "P1B_challis_detection.json"
PRERUN = NOTES + "P1B_challis_prerun.json"
RECORDINGS = "data/historical/adams_challis/prediscovery_recordings.csv"
SITE_DATA = ROOT / "site" / "data"
GALLE = "1846-09-23"
MAP_NIGHT = "1846-08-12"
CASES = ("A0", "C1a", "C1b", "C2_1800", "C2_date", "C3_S4", "C3_S4_prompt", "C3_S5",
         "C3_S5_prompt", "C4a", "C4a_prompt", "C4b")  # fmt: skip


def resolve(source: list, cache: dict | None = None):
    """The value a source names: ``[file, key, ...]`` followed through dicts and lists."""
    cache = {} if cache is None else cache
    name = source[0]
    if name not in cache:
        cache[name] = json.loads((ROOT / name).read_text())
    x = cache[name]
    for k in source[1:]:
        x = x[k]
    return x


def _collect(sources: dict) -> dict:
    cache: dict = {}
    return {k: resolve(s, cache) for k, s in sources.items()}


def fit_epoch() -> dict:
    """The epoch of the held-distance fits' longitudes (`fit_epoch_jd_tdb`, step 3b's comparison
    epoch) as a Julian year and a calendar date, so that the page states it and not Galle's
    night, 12 days earlier."""
    jd = resolve([REPRO, "R7_d4", "epoch_jd_tdb"], {})
    day = int(jd + 0.5) + 68569  # Fliegel and Van Flandern (1968), Gregorian date of a JD
    n = 4 * day // 146097
    day -= (146097 * n + 3) // 4
    i = 4000 * (day + 1) // 1461001
    day = day - 1461 * i // 4 + 31
    j = 80 * day // 2447
    d = day - 2447 * j // 80
    day = j // 11
    m, y = j + 2 - 12 * day, 100 * (n - 49) + i + day
    return {
        "julian_year": round(2000.0 + (jd - 2451545.0) / 365.25, 3),
        "date": f"{y}-{m:02d}-{d:02d}",
        "note": "the epoch at which step 3b compares every solution (TDB)",
    }


def adams_sources() -> dict:
    s: dict = {
        "fit_data_set": [FIT, "data_set", "fit"],
        "post_hoc_definition": [NOISE, "results", "H1/a", "post_hoc_on_branch", "definition"],
        "fit_n_observations": [FIT, "data_set", "fit", "total", 0],
        "fit_n_measurements": [FIT, "data_set", "fit", "total", 1],
        "neptune_helio_true_galle": [FIT, "neptune_places", GALLE, "helio_true"],
        "free_distance_status": [FIT, "reading", "g", "statement_1", "F", "status"],
        "truth_r_au": [REPRO, "R6_d3", "truth", "neptune_1846", "r_au"],
        "neptune_same_model": [REPRO, "R7_d4", "neptune_same_model_of_date"],
        "fit_epoch_jd_tdb": [REPRO, "R7_d4", "epoch_jd_tdb"],
        "o_clean_least_squares": [REPRO, "R6_d3", "O_i_d1c"],
        "o_clean_his_procedure": [REPRO, "R6_d3", "O_i_his_procedure"],
        "o_clean_project_on_his_places": [REPRO, "R6_d3", "O_ii"],
        "o_clean_project_fit": [REPRO, "R6_d3", "iii", "O_iii_median"],
        "o_clean_project_fit_68": [REPRO, "R6_d3", "iii", "O_iii_68"],
        "o_real_project_fit": [REPRO, "R6_d3", "O_3_step3"],
        "transcription_passes": [REPRO, "R1_transcription", "passes"],
        "n_coefficient_printed": [REPRO, "R3_locate", "H1", "n_coefficients", "dy1", "printed"],
        "n_coefficient_exact": [REPRO, "R3_locate", "H1", "n_coefficients", "dy1", "exact"],
        "noise_n_draws": [NOISE, "results", "H1/a", "n_draws"],
        "hist_edges_deg": [HIST, "edges_deg"],
        "hist_bin_deg": [HIST, "bin_deg"],
    }
    for t in ("S4", "S5", "S6"):
        s[f"{t}_helio_true_galle"] = [FIT, "targets", t, "places", GALLE, "helio_true"]
        s[f"{t}_a_au"] = [FIT, "targets", t, "elements", "a_au"]
        s[f"{t}_date"] = [FIT, "targets", t, "date"]
    for h in ("H1", "H2"):
        v = f"{h.lower()}_g"
        for q in ("median", "lo68", "hi68"):
            s[f"{h}_fit_true_lon_{q}"] = [REPRO, "R7_d4", "variants", v, "true_lon_1846_of_date", q]
        # His printed solution carried to the same epoch as the fits (Fig. 1 and § 3.1).
        s[f"{h}_printed_true_lon"] = [REPRO, "R7_d4", "variants", v, "true_lon_1846_of_date",
                                      "vs", "adams_printed", "target"]  # fmt: skip
        s[f"{h}_fit_status"] = [FIT, "reading", "g", "statement_1", h, "status"]
        s[f"{h}_c0_theta_delta"] = [REPRO, "R3_stage_c", f"{h}/C0", "vs_print", "rows", "theta",
                                    "delta"]  # fmt: skip
        s[f"{h}_c2_vs_print"] = [REPRO, "R3_stage_c", f"{h}/C2", "vs_print", "rows",
                                 "mean_lon_1846", "delta"]  # fmt: skip
        s[f"{h}_c2_vs_brookes"] = [REPRO, "R3_stage_c", f"{h}/C2", "vs_brookes", "rows",
                                   "mean_lon_1846", "delta"]  # fmt: skip
        s[f"{h}_dtheta_per_step"] = [REPRO, "R3_locate", h, "n_coefficients", "dy1",
                                     "dtheta_per_step"]  # fmt: skip
        b = [NOISE, "multistart_baselines", h]
        s[f"{h}_o_real_least_squares"] = [*b, "real", "O"]
        s[f"{h}_delta_obs"] = [*b, "delta_obs", "primary"]
        for q in ("theta", "e_prime", "chi2_minus_best", "O"):
            s[f"{h}_second_minimum_{q}"] = [*b, "second_minimum", "real", q]
        s[f"{h}_prereg_verdict"] = [NOISE, "readings", f"{h}/d1c/primary", "verdict"]
        s[f"{h}_post_hoc_verdict"] = [NOISE, "post_hoc_readings", f"{h}/d1c/primary", "verdict"]
        for m in ("a", "b"):
            s[f"{h}_prereg_P_{m}"] = [NOISE, "readings", f"{h}/d1c/primary", f"P_{m}"]
            s[f"{h}_post_hoc_P_{m}"] = [NOISE, "post_hoc_readings", f"{h}/d1c/primary", f"P_{m}"]
            r = [NOISE, "results", f"{h}/{m}"]
            s[f"{h}_{m}_fraction_on_branch"] = [*r, "post_hoc_on_branch", "fraction_on_branch"]
            s[f"{h}_{m}_exceeding_branch_changes"] = [
                *r, "branches", "fraction_of_exceeding_that_are_branch_changes"]  # fmt: skip
            s[f"{h}_{m}_hist_counts"] = [HIST, "cells", f"{h}/{m}", "counts"]
            s[f"{h}_{m}_hist_n"] = [HIST, "cells", f"{h}/{m}", "n"]
    return s


def challis_sources() -> dict:
    pre = json.loads((ROOT / PRERUN).read_text())
    idx = {n["date"]: i for i, n in enumerate(pre["nights"])}
    s: dict = {
        "n_accepted": [DET, "n_accepted"],
        "n_drawn": [DET, "acceptance", "n_drawn"],
        "acceptance": [DET, "acceptance", "all"],
        "acceptance_conditions_2_3": [DET, "acceptance", "conditions_2_3"],
        "acceptance_condition_4": [DET, "acceptance", "condition_4"],
        "halved_n_accepted": [DET, "acceptance_halved_priors", "n_accepted"],
        "tornado_n_draws": [DET, "tornado", "n_draws"],
        "c2_combined": [DET, "c2_combined_reading_row"],
        "replay_a0": [DET, "cases", "A0", "replay", "primary", "mean"],
        "c1a_done_by": [DET, "c1_completion", "C1a", "date_95_percent_of_flags"],
        "c1a_p_by_0816": [DET, "c1_completion", "C1a", "p_by", "1846-08-16"],
        "v10": [PRERUN, "constants", "V10_neptune", 0],
        "v10_source": [PRERUN, "constants", "V10_neptune", 1],
        "lv_august_true_lon": [PRERUN, "le_verrier_august", "true_lon_1847_deg"],
        "lv_august_limits": [PRERUN, "le_verrier_august", "limits_helio_lon_deg"],
        "lv_june_limits": [PRERUN, "le_verrier_june_limits_deg"],
    }
    for p in ("p_disc_routine", "p_disc_looking", "p_follow"):
        for w in ("prior", "posterior"):
            s[f"{p}_{w}_mean"] = [DET, "posterior_against_prior", p, w, "mean"]
    for c in CASES:
        b = [DET, "cases", c]
        for q in ("mean", "q05", "q95"):
            s[f"{c}_{q}"] = [*b, "held", "primary", q]
        s[f"{c}_row"] = [*b, "reading_row"]
        s[f"{c}_top_input"] = [*b, "top_input"]
        s[f"{c}_conditional_on"] = [*b, "conditional_on"]
        s[f"{c}_prior_dependent"] = [*b, "prior_dependent"]
        s[f"{c}_halved_q05"] = [*b, "halved_priors", "primary", "q05"]
        s[f"{c}_halved_q95"] = [*b, "halved_priors", "primary", "q95"]
        s[f"{c}_halved_row"] = [*b, "halved_priors", "reading_row"]
    for d in ("1846-07-29", "1846-07-30", "1846-08-04", "1846-08-12", "1846-09-22", "1846-09-29",
              "1846-10-01"):  # fmt: skip
        n = [PRERUN, "nights", idx[d]]
        key = d[5:].replace("-", "")
        for q in ("ra_h", "dec_deg", "ra1800_h", "V", "diameter_arcsec"):
            s[f"neptune_{key}_{q}"] = [*n, "neptune", q]
    for t in ("col325", "S4", "S5", "LVJ", "LVA"):
        n = [PRERUN, "nights", idx[MAP_NIGHT], "targets", t]
        s[f"place_{t}_ra_h"] = [*n, "ra_h"]
        s[f"place_{t}_dec_deg"] = [*n, "dec_deg"]
    return s


def _track() -> list[dict]:
    pre = json.loads((ROOT / PRERUN).read_text())
    return [{"date": n["date"], "ra_h": n["neptune"]["ra_h"], "dec_deg": n["neptune"]["dec_deg"]}
            for n in pre["nights"]]  # fmt: skip


def _recordings() -> list[dict]:
    with (ROOT / RECORDINGS).open() as f:
        return [{"observer": r["observer"], "nights": int(r["nights"]), "counted": r["counted"]}
                for r in csv.DictReader(f)]  # fmt: skip


# Cited constants: each carries its citation (docs/DESIGN.md, the one rule).
SWEEPS_CITE = "Challis 1847, Memoirs RAS 16, 415–426 (challis_sweeps.csv)"
CONSTANTS = {
    "zone_breadth_arcmin": {"value": 9, "citation": "Challis 1846, MNRAS 7, 145: the field at "
                            "power 166"},
    "power": {"value": 166, "citation": "Challis 1846, MNRAS 7, 145"},
    "aperture": {"text": "a little more than 11½ inches", "citation": "Airy 1844, Account of "
                 "the Northumberland Equatoreal and Dome"},
    "sheet_west_edge_ra": {"h": 21, "m": 56, "equinox": 1800,
                           "citation": "Krajnović 2021 (secondary): Hora "
                           "XXII, one hour of right ascension with a 4-minute overlap"},
    "rule_thresholds": {"low": 0.2, "high": 0.8, "citation": "the pre-registered reading rule, "
                        "docs/notes/P1B_challis_detection.md § 7"},
    "galle_night": {"day": 23, "month": "September", "year": 1846},
    "range_percent": {"low": 5, "high": 95, "citation": "the pre-registered reading rule, "
                      "docs/notes/P1B_challis_detection.md § 7"},
    "fit_ranges_percent": {"inner": 68, "outer": 95, "citation": "the pre-registered reading "
                           "rule, docs/notes/P1B_adams_fit.md"},
    "planet_number_12_august": {"value": 49, "citation": SWEEPS_CITE},
    "comparison_stopped_at": {"value": 39, "citation": SWEEPS_CITE},
    "events": {"search_start": "1846-07-29", "zone_first_swept": "1846-07-30",
               "recorded_first": "1846-08-04", "recorded_second": "1846-08-12",
               "disc_noted": "1846-09-29", "track_end": "1846-10-01", "citation": SWEEPS_CITE},
    "noise_thresholds": {"suffices": 0.05, "systematic": 0.003, "citation": "the "
                         "pre-registered reading rule, docs/notes/P1B_adams_noise.md"},
}  # fmt: skip


# The page's references: [text, italic title, text] parts; a string in braces is set in italics.
REFERENCES = [
    ["adams", "Adams, J. C. 1847. An explanation of the observed irregularities in the motion of "
     "Uranus, on the hypothesis of disturbances caused by a more distant planet. ",
     "{Memoirs of the Royal Astronomical Society}", " 16, 427–460."],
    ["airy", "Airy, G. B. 1844. ", "{Account of the Northumberland Equatoreal and Dome}",
     ". Cambridge: University Press."],
    ["brookes", "Brookes, C. J. 1970. On the prediction of Neptune. ", "{Celestial Mechanics}",
     " 3, 67–80."],
    ["challis-mn", "Challis, J. 1846. Account of observations at the Cambridge Observatory for "
     "detecting the planet exterior to Uranus. ",
     "{Monthly Notices of the Royal Astronomical Society}", " 7, 145–149."],
    ["challis-mem", "Challis, J. 1847. An account of observations undertaken in search of the "
     "planet discovered at Berlin on Sept. 23, 1846. ",
     "{Memoirs of the Royal Astronomical Society}", " 16, 415–426."],
    ["hind", "Hind, J. R. 1850. Auszug aus einem Briefe des Herrn Hind. ",
     "{Astronomische Nachrichten}", " 30, 257–258."],
    ["kowal", "Kowal, C. T. and Drake, S. 1980. Galileo's observations of Neptune. ", "{Nature}",
     " 287, 311."],
    ["krajnovic", "Krajnović, D. 2021. “That star is not on the map”: the German side of the "
     "discovery. In ", "{The Discovery of Neptune}", ". Springer. arXiv:2108.06305."],
    ["lv-june", "Le Verrier, U. J. 1846. Recherches sur les mouvements d'Uranus. ",
     "{Comptes rendus}", " 22, 917 (1 June)."],
    ["lv-aug", "Le Verrier, U. J. 1846. Sur la planète qui produit les anomalies observées dans "
     "le mouvement d'Uranus. ", "{Comptes rendus}", " 23, 428–438 (31 August)."],
    ["mallama", "Mallama, A. and Hilton, J. L. 2018. Computing apparent planetary magnitudes "
     "for The Astronomical Almanac. ", "{Astronomy and Computing}", " 25, 10–24."],
    ["mauvais", "Mauvais, V. 1847. Sur une observation inédite de la nouvelle planète. ",
     "{Comptes rendus}", " 24, 666–669."],
    ["schonfeld", "Schönfeld, E. 1886. Notiz über eine Beobachtung des Neptun in den "
     "Lamont'schen Zonen. ", "{Astronomische Nachrichten}", " 114, 59–60."],
]  # fmt: skip


def _sha(name: str) -> str:
    return hashlib.sha256((ROOT / name).read_bytes()).hexdigest()


def build() -> dict[str, dict]:
    a, c = adams_sources(), challis_sources()
    files = (FIT, REPRO, NOISE, HIST, DET, PRERUN, RECORDINGS)
    return {
        "p1b_adams.json": {
            "description": "Companion page, § 3: values copied from the committed step 3, 3b "
            "and 3c JSONs; `sources` names the field of each.",
            "sources": a, "values": _collect(a),
            "constants": {"fit_epoch": fit_epoch()},
        },
        "p1b_challis.json": {
            "description": "Companion page, § 4: values copied from the committed step 4 JSONs; "
            "`sources` names the field of each. `track` is Neptune's place per night from the "
            "pre-run; `recordings` the pre-discovery table's counts.",
            "sources": c, "values": _collect(c), "track": _track(),
            "recordings": _recordings(), "constants": CONSTANTS, "map_night": MAP_NIGHT,
            "case_order": list(CASES),
        },
        "p1b_page.json": {
            "description": "Companion page: where it lives and what it was built from.",
            "path": "adams-challis/",
            "title": "The British search for Neptune: Adams's solutions and Challis's sweeps",
            "status": "part of release v1.2.0",
            "main_title": "Refitting Le Verrier's prediction of Neptune",
            "main_link": "A companion study treats the British search of the same years: ",
            "references": [{"id": r[0], "parts": r[1:]} for r in REFERENCES],
            "notes": ["P1B_sources.md", "P1B_adams_solutions.md", "P1B_challis_sweeps.md",
                      "P1B_adams_fit.md", "P1B_adams_reproduction.md", "P1B_adams_noise.md",
                      "P1B_challis_detection.md"],
            "built_from_sha256": {f: _sha(f) for f in files},
        },
    }  # fmt: skip


def main() -> None:
    for name, payload in build().items():
        path = SITE_DATA / name
        path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"wrote {path} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
