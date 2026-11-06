"""M6a/A6a (transcription → schema) and the A6b loader test (equinox/precession handling)."""

import csv
import hashlib
import json
import math

import numpy as np
import pytest
from astropy.time import Time

from fork.data import frames
from fork.data import observations as obs
from fork.data import transcription as tr
from fork.data.spotcheck import SELECTION_FILE

RESULTS_FILE = tr.ROOT / "docs" / "review" / "spotcheck_results.json"


@pytest.fixture(scope="module")
def rows():
    return obs.load()


def _by_id(rows):
    return {r["id"]: r for r in rows}


# --- A6a ------------------------------------------------------------------------------------


def test_committed_file_is_the_build_of_the_transcription_and_misprint_table(rows):
    table = tr.read_pass(tr.RECONCILED_FILE)
    built = obs.build_rows(table, obs._misprints_by_cell(table))
    assert rows == built
    assert len({r["obs_id"] for r in rows}) == 279
    assert sum(r["quantity"] == "ra" for r in rows) == 279
    assert sum(r["era"] == "ancient" for r in rows if r["quantity"] == "ra") == 19


def test_spot_check_found_no_uncorrected_errors_in_the_current_table():
    """A6a: the 10% spot check covered the sampled rows of *this* table and found nothing."""
    res = json.loads(RESULTS_FILE.read_text(encoding="utf-8"))
    assert res["table77_sha256"] == hashlib.sha256(tr.RECONCILED_FILE.read_bytes()).hexdigest()
    with SELECTION_FILE.open(newline="", encoding="utf-8") as fh:
        sel = [(r["page"], r["line"]) for r in csv.DictReader(fh)]
    assert [(r["page"], r["line"]) for r in res["rows"]] == sel
    cells = [c for r in res["rows"] for c in r["cells"].values()]
    assert len(cells) == 28 * 12
    assert all(c["matches"] for c in cells)


def test_ancient_row_ra_in_arc_and_computed_seconds_in_arcsec(rows):
    """(a) 1690 Dec 23 (p. 129 row 1): RA printed 55° 49′ 19.7″ (arc, not time), computed
    seconds 12.8″, C − O −66.9″ (more than a minute of arc)."""
    ra = _by_id(rows)["LV1846-T77-p129-l01-ra"]
    assert ra["era"] == "ancient" and ra["unit_as_printed"] == "deg arcmin arcsec"
    assert float(ra["value"]) == pytest.approx(55 + 49 / 60 + 19.7 / 3600, abs=1e-9)
    assert float(ra["lv_residual"]) == pytest.approx(+66.9, abs=1e-9)  # O − C = −(C − O)
    # computed RA = observed + (C − O); its arcsecond part is the printed computed seconds
    calc_arcsec = float(ra["value"]) * 3600 - float(ra["lv_residual"])
    assert calc_arcsec % 60 == pytest.approx(12.8, abs=1e-6)
    assert ra["observer"] == "Flamsteed" and ra["observatory"] == "Greenwich"


def test_computed_seconds_reproduce_from_value_and_residual_in_both_eras(rows):
    """Units check over every row: ancient RA residuals are arcsec of arc, modern ones 15 × s
    of time. Rebuilding the printed computed seconds from value and lv_residual must work."""
    table = {
        f"LV1846-T77-p{r['page']}-l{int(r['line']):02d}": r
        for r in tr.read_pass(tr.RECONCILED_FILE)
    }
    checked = 0
    for r in rows:
        if r["quantity"] != "ra" or not r["lv_residual"]:
            continue
        printed = float(table[r["obs_id"]]["ra_calc_sec"])
        if r["era"] == "ancient":
            calc = (float(r["value"]) * 3600 - float(r["lv_residual"])) % 60  # arcsec
        else:
            calc = (float(r["value"]) * 240 - float(r["lv_residual"]) / 15) % 60  # s of time
        assert min(abs(calc - printed), 60 - abs(calc - printed)) < 1e-3, r["id"]
        checked += 1
    assert checked == 279


def test_negative_zero_degree_declination_keeps_its_sign(rows):
    """(b) 1843 Sept. 20 (p. 136 row 18): printed '− 0,33. 2,7' (comma as degree separator),
    transcribed '-0 33 2.7': −0° 33′ 2.7″, not +0° 33′ 2.7″ or −0.33°."""
    dec = _by_id(rows)["LV1846-T77-p136-l18-dec"]
    assert dec["value_as_printed"] == "-0 33 2.7"
    assert float(dec["value"]) == pytest.approx(-(33 / 60 + 2.7 / 3600), abs=1e-9)
    assert obs.parse_sexagesimal("+0 16 24.5") == pytest.approx(16 / 60 + 24.5 / 3600, abs=1e-12)
    for bad in ("-0.33 2.7", "-0,33 2.7", "-0 33", "1 6? 2.0", "» » »", "-0 61 0.0"):
        with pytest.raises(ValueError):
            obs.parse_sexagesimal(bad)


def test_ra_only_rows_have_no_dec_row(rows):
    ids = _by_id(rows)
    assert "LV1846-T77-p129-l08-ra" in ids  # Bradley, 1753 Dec 3
    assert "LV1846-T77-p129-l08-dec" not in ids
    assert ids["LV1846-T77-p129-l08-ra"]["observer"] == "Bradley"
    assert len(rows) - 279 == sum(tr.NO_VALUE not in r["dec_obs"]
                                  for r in tr.read_pass(tr.RECONCILED_FILE))  # fmt: skip


@pytest.mark.parametrize(
    ("year", "month", "observatory"),
    [(1756, 9, "Göttingen"), (1750, 10, "Paris"), (1781, 9, "Greenwich"), (1800, 3, "Greenwich"),
     (1801, 3, "Paris"), (1807, 1, "Greenwich"), (1807, 4, "Paris"), (1827, 7, "Paris"),
     (1828, 7, "ambiguous"), (1829, 7, "Greenwich"), (1830, 11, "Greenwich"),
     (1835, 7, "Paris"), (1845, 9, "Paris")],
)  # fmt: skip
def test_observatory_from_page_127(year, month, observatory):
    assert obs.attribute(year, month).observatory == observatory
    with pytest.raises(ValueError):
        obs.attribute(1832, 1)  # no source for 1831–1834 (p. 127)


def test_times_are_paris_mean_time_by_the_transit_condition(rows):
    """Meridian observations: apparent sidereal time at the site equals the observed RA.
    With Paris mean time (astronomical day) every row fits to < 4 s except the one logged
    1768 row; the 1828 rows fit Greenwich and miss Paris by the 9 m 21 s offset."""
    ra = [r for r in rows if r["quantity"] == "ra"]
    bad = {
        r["obs_id"] for r in ra if r["transit_resid_s"] and abs(float(r["transit_resid_s"])) >= 4
    }
    assert bad == {"LV1846-T77-p129-l12"}
    amb = [r for r in ra if r["observatory"] == "ambiguous"]
    assert len(amb) == 2 and all(
        "Greenwich +" in r["notes"] and "Paris +56" in r["notes"] for r in amb
    )
    # Greenwich rows timed in Greenwich mean time would miss by ~563 s: check one directly.
    g = next(r for r in ra if r["obs_id"] == "LV1846-T77-p129-l30")  # 1783 Oct 7
    ut_if_gmt = float(g["epoch_ut_jd"]) + obs.LON_EAST_H["Paris"] / 24
    off = obs.transit_residual_s(ut_if_gmt, 0.0, "Greenwich", float(g["value"]) / 15)
    assert abs(off) > 500


def test_epoch_astronomical_day_paris_mean_time(rows):
    """1783 Oct. 7, 17 54 6 Paris mean time (astronomical) = 1783-10-08 05:54:06 at Paris
    = 05:44:45.1 UT."""
    r = _by_id(rows)["LV1846-T77-p129-l30-ra"]
    expect = Time("1783-10-08 05:44:45.1", scale="ut1", precision=1).jd
    assert float(r["epoch_ut_jd"]) == pytest.approx(expect, abs=0.2 / 86400)
    dt = (float(r["epoch_tt_jd"]) - float(r["epoch_ut_jd"])) * 86400
    assert 5 < dt < 20


def test_delta_t_polynomials_are_continuous():
    assert obs.delta_t_s(1700.0) == pytest.approx(8.83, abs=1e-9)
    assert obs.delta_t_s(1800.0) == pytest.approx(13.72, abs=1e-9)
    for y in (1700.0, 1800.0):
        assert abs(obs.delta_t_s(y - 1e-6) - obs.delta_t_s(y)) < 0.2


def test_misprints_applied_only_where_tabled(rows):
    applied = {r["obs_id"]: r["misprints_applied"] for r in rows if r["misprints_applied"]}
    assert applied == {
        "LV1846-T77-p132-l05": "mean_time: 12 18 13 -> 12 17 13",
        "LV1846-T77-p132-l21": "mean_time: 11 32 34 -> 11 32 54",
        "LV1846-T77-p135-l03": "dec_obs: -18 12 6.9 -> -19 12 6.9",
    }
    dec = _by_id(rows)["LV1846-T77-p135-l03-dec"]
    assert dec["value_as_printed"] == "-18 12 6.9"
    assert float(dec["value"]) == pytest.approx(-(19 + 12 / 60 + 6.9 / 3600), abs=1e-9)


# --- A6b loader test --------------------------------------------------------------------------


def _record(t_jd, **kw):
    r = {"id": "syn", "obs_id": "syn", "epoch_tt_jd": f"{t_jd:.8f}", "quantity": "lon",
         "value": "0", "unit": "deg", "equinox_ref": "fixed", "equinox_epoch": "",
         "equinox_citation": "synthetic", "precession_constant_arcsec_per_yr": "50.10",
         "precession_citation": "synthetic", "source_ref": "synthetic"}  # fmt: skip
    r.update(kw)
    return r


def _sep_arcsec(a, b):
    return abs((a - b + math.pi) % (2 * math.pi) - math.pi) * obs.ARCSEC_PER_RAD


def test_a6b_fixed_epoch_reduction_is_undone_with_the_sources_constant():
    """J2000 longitudes at 1690–1846 → of date (IAU 2006) → pre-reduced to 1800.0 with
    50.10″/yr, as a source would. The loader recovers J2000 to < 0.01″; treating the values as
    of date, or undoing the reduction with IAU 2006, is off by > 1″ somewhere."""
    rng = np.random.default_rng(6)
    years = np.linspace(1690.0, 1846.0, 14)
    t = obs.J2000_JD + (years - 2000.0) * obs.JULIAN_YEAR_D
    epoch = obs.J2000_JD - 200.0 * obs.JULIAN_YEAR_D  # 1800.0
    lon_j = rng.uniform(0, 2 * math.pi, years.size)
    lat_j = rng.uniform(-0.01, 0.01, years.size)
    lon_d, lat_d = frames.ecl_j2000_to_ecl_of_date(lon_j, lat_j, t)
    lon_1800 = np.mod(lon_d - 50.10 * (t - epoch) / obs.JULIAN_YEAR_D / obs.ARCSEC_PER_RAD,
                      2 * math.pi)  # fmt: skip

    good, wrong_a, wrong_b = [], [], []
    for i in range(years.size):
        r = _record(t[i], equinox_epoch=f"{epoch:.1f}")
        lon, lat = obs.lon_lat_ecl_j2000(r, lon_1800[i], lat_d[i])
        good.append(max(_sep_arcsec(lon, lon_j[i]), abs(lat - lat_j[i]) * obs.ARCSEC_PER_RAD))
        la, _ = frames.ecl_of_date_to_ecl_j2000(lon_1800[i], lat_d[i], t[i])  # (a) as of date
        lb, _ = frames.ecl_of_date_to_ecl_j2000(lon_1800[i], lat_d[i], epoch)  # (b) IAU 2006
        wrong_a.append(_sep_arcsec(float(la), lon_j[i]))
        wrong_b.append(_sep_arcsec(float(lb), lon_j[i]))
    assert max(good) < 0.01
    assert max(wrong_a) > 1.0
    assert max(wrong_b) > 1.0


def test_a6b_records_without_equinox_or_constant_fail_to_load():
    t = obs.J2000_JD - 100 * obs.JULIAN_YEAR_D
    for kw in ({"equinox_ref": ""}, {"precession_constant_arcsec_per_yr": ""},
               {"equinox_citation": ""}, {"precession_citation": ""},
               {"equinox_epoch": ""},  # fixed without its epoch
               {"precession_constant_arcsec_per_yr": "not_stated", "equinox_epoch": "2378496.0"},
               {"equinox_ref": "mean_1800"}):  # fmt: skip
        with pytest.raises(ValueError):
            obs.check_record(_record(t, **{"equinox_epoch": "2378496.0", **kw}))
    obs.check_record(_record(t, equinox_epoch="2378496.0"))  # complete: loads
    obs.check_record(
        _record(t, equinox_ref="of_date", precession_constant_arcsec_per_yr="not_stated")
    )
