"""Le Verrier (1846) Table n° 77 → ``data/historical/leverrier_1846.csv`` (spec M6a schema),
and the observation loader with its equinox/precession rules (A6b loader test).

    uv run python -m fork.data.observations          # rebuild leverrier_1846.csv

Input: the reconciled transcription ``leverrier_1846_table77.csv`` (exactly as printed) and the
misprint table ``leverrier_1846_misprints.csv`` (its ``apply`` entries are substituted here, so
the transcription itself never changes).

One CSV row per *measured quantity*: an observation's RA and Dec are two rows sharing
``obs_id``, so an RA-only observation (e.g. Bradley 1753) simply has no Dec row and gets a 1-D
likelihood downstream. Values are degrees (``value``), with the printed string kept beside them.

What the memoir fixes (CdT 1849 Additions, page images; see docs/notes/M6.md):

* **Frame (p. 111 §72, p. 127).** Observed and computed places are apparent geocentric RA and
  Dec: true equator and equinox of date, with nutation and aberration. So ``equinox_ref`` is
  ``of_date``. §72 states no precession constant for the ephemeris; for an of-date position no
  fixed-epoch reduction has to be undone, so the constant is recorded as ``not_stated`` with
  that citation. Any record with ``equinox_ref = fixed`` must carry a numeric constant.
* **Time.** "Temps moyen" is Paris mean time, astronomical reckoning (the day begins at mean
  noon). p. 126 says so for Mayer's row; for every other row it is established by the transit
  check below (meridian observations: local apparent sidereal time = observed RA). Greenwich
  rows fit Paris time to ≤ 3.3 s and miss Greenwich time by ≈ 563 s (9 m 21 s sidereal).
* **Observers and observatories.** Ancient rows from §74 (pp. 124–126): Flamsteed (Greenwich)
  1690, 1712, 1715; Le Monnier (Paris) 1750, 1764–1771; Bradley (Greenwich) 1753; Mayer
  (Göttingen) 1756. Modern rows from p. 127: Greenwich 1781–1800, January 1807/1808/1809, and
  1828–1830; Paris 1801–1828 and 1835–1845. 1828 is in both ranges: ``ambiguous`` (pending
  Airy's *Reduction*), with the transit evidence in ``notes``.
* **Units as printed.** Ancient RA (before 1781) in arc, ° ′ ″, with the computed seconds and
  C − O in ″; modern RA in time, h m s, with computed seconds and C − O in s of time. Dec always
  ° ′ ″. Southern computed Dec seconds carry the Dec's sign. ``lv_residual`` is Le Verrier's
  O − C (the table prints C − O) in arcseconds of the coordinate itself (RA: 15 × s of time,
  not multiplied by cos δ).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import erfa
import numpy as np

from fork.data import frames
from fork.data import transcription as tr

OUT_FILE = tr.HIST / "leverrier_1846.csv"
SOURCE = "CdT1849Add"  # Connaissance des Temps pour 1849, Additions (Le Verrier 1846)
FIRST_MODERN_YEAR = 1781
JULIAN_YEAR_D = 365.25
J2000_JD = 2451545.0
ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi

MONTHS = {
    "Janv.": 1,
    "Fév.": 2,
    "Mars": 3,
    "Avril": 4,
    "Mai": 5,
    "Juin": 6,
    "Juin.": 6,
    "Juill.": 7,
    "Juill": 7,
    "Août": 8,
    "Sept.": 9,
    "Oct.": 10,
    "Nov.": 11,
    "Déc.": 12,
}  # fmt: skip  (as printed, including the two spellings without/with a stray dot)

# East longitudes in hours. Paris Observatory 2° 20′ 14″ E; Göttingen (Mayer's observatory)
# about 9° 56′ E. Only Paris enters the time conversion; site longitudes are used only by the
# transit check, where 1 s of time is the resolution.
LON_EAST_H = {"Greenwich": 0.0, "Paris": (9 * 60 + 20.9) / 3600, "Göttingen": (39 * 60 + 46) / 3600}
TIME_MERIDIAN = "Paris"

COLUMNS = (
    "id", "obs_id", "epoch", "epoch_tt_jd", "epoch_ut_jd", "time_meridian", "time_reckoning",
    "observer", "observatory", "instrument", "quantity", "value", "unit", "value_as_printed",
    "unit_as_printed", "frame_as_printed", "equinox_ref", "equinox_epoch", "equinox_citation",
    "precession_constant_arcsec_per_yr", "precession_rate_arcsec_per_yr2",
    "precession_citation", "lv_residual", "lv_residual_unit", "lv_c_minus_o_as_printed",
    "lv_weight", "lv_dlon_c_minus_o_as_printed", "lv_dlat_c_minus_o_as_printed",
    "transit_resid_s", "era", "source_ref", "misprints_applied", "notes",
)  # fmt: skip

FRAME_AS_PRINTED = "apparent geocentric; true equator and equinox of date (nutation, aberration)"
EQUINOX_CITATION = f"{SOURCE}:p111:§72 (ephemeris n° 75: apparent RA/Dec); p127:§75"
PRECESSION_CITATION = (
    f"{SOURCE}:p111:§72: no precession constant stated; of-date positions, so no fixed-epoch "
    "reduction to undo"
)


# --- parsing ------------------------------------------------------------------------------


def parse_sexagesimal(s: str) -> float:
    """'±a b c.f' → signed value in units of the first part (a + b/60 + c/3600).

    The sign is read from the string, so '-0 33 2.7' is −0.5507…, not +0.5507… (a −0°
    declination, printed '− 0,33. 2,7' in 1843). Exactly three parts; only the last may have
    a decimal point (a degree separator mis-normalized to '.', e.g. '-0.33 2.7', raises rather
    than loading as −0.33°). Unreadable ('?') or empty ('»') cells raise."""
    t = s.strip()
    if not t or "?" in t or tr.NO_VALUE in t:
        raise ValueError(f"no value in {s!r}")
    sign = -1.0 if t.startswith("-") else 1.0
    parts = t.lstrip("+-").split()
    if len(parts) != 3:
        raise ValueError(f"expected 3 sexagesimal parts in {s!r}")
    if any("." in p for p in parts[:2]) or not all(p.isdigit() for p in parts[:2]):
        raise ValueError(f"only the last part may be fractional in {s!r}")
    a, b, c = int(parts[0]), int(parts[1]), float(parts[2])
    if b >= 60 or c >= 60:
        raise ValueError(f"minutes or seconds out of range in {s!r}")
    return sign * (a + b / 60.0 + c / 3600.0)


# --- observatory, observer -----------------------------------------------------------------


@dataclass(frozen=True)
class Attribution:
    observer: str
    observatory: str  # "Greenwich", "Paris", "Göttingen", or "ambiguous"
    citation: str
    note: str = ""


ANCIENT = {
    1690: ("Flamsteed", "Greenwich", "p124:§74"),
    1712: ("Flamsteed", "Greenwich", "p124:§74"),
    1715: ("Flamsteed", "Greenwich", "p125:§74"),
    1750: ("Le Monnier", "Paris", "p125-126:§74"),
    1753: ("Bradley", "Greenwich", "p126:§74"),
    1756: ("Mayer", "Göttingen", "p126:§74"),
    1764: ("Le Monnier", "Paris", "p126:§74"),
    1768: ("Le Monnier", "Paris", "p126:§74"),
    1769: ("Le Monnier", "Paris", "p126:§74"),
    1771: ("Le Monnier", "Paris", "p126:§74"),
}


def attribute(year: int, month: int) -> Attribution:
    """Observer and observatory from the memoir's text (the table has no such column)."""
    if year < FIRST_MODERN_YEAR:
        if year not in ANCIENT:
            raise ValueError(f"no attribution for ancient year {year}")
        who, where, ref = ANCIENT[year]
        return Attribution(who, where, f"{SOURCE}:{ref}")
    ref = f"{SOURCE}:p127:§75"
    if year == 1828:
        return Attribution("", "ambiguous", ref, "p127 lists 1828 under both Greenwich and Paris")
    if 1781 <= year <= 1800 or (year in (1807, 1808, 1809) and month == 1) or year in (1829, 1830):
        return Attribution(
            "", "Greenwich", ref, "Greenwich, via Airy's Reduction (re-reduced by LV)"
        )
    if 1801 <= year <= 1827:
        return Attribution("", "Paris", ref, "Paris, CdT and two folio volumes")
    if 1835 <= year <= 1845:
        return Attribution("", "Paris", ref, "Paris, unpublished series given by Arago")
    raise ValueError(f"no observatory range on p127 covers {year}-{month:02d}")


# --- time ----------------------------------------------------------------------------------


def delta_t_s(year: float) -> float:
    """TT − UT in seconds, Espenak & Meeus (2006) polynomials (NASA Five Millennium Canon),
    1600–1860. Over 1690–1845 ΔT is 6–14 s: Uranus moves < 0.01″ in that time."""
    if 1600 <= year < 1700:
        t = year - 1600
        return 120 - 0.9808 * t - 0.01532 * t**2 + t**3 / 7129
    if 1700 <= year < 1800:
        t = year - 1700
        return 8.83 + 0.1603 * t - 0.0059285 * t**2 + 0.00013336 * t**3 - t**4 / 1174000
    if 1800 <= year < 1860:
        t = year - 1800
        c = (13.72, -0.332447, 0.0068612, 0.0041116, -0.00037436, 0.0000121272, -0.0000001699,
             0.000000000875)  # fmt: skip
        return sum(k * t**i for i, k in enumerate(c))
    raise ValueError(f"ΔT polynomial not set up for {year}")


def epoch_ut_jd(year: int, month: int, day: int, mean_time_h: float) -> float:
    """UT Julian Date of a Paris mean time on an astronomical day.

    The astronomical day ``day`` begins at mean noon of the civil date ``day``, which is where
    a Julian Date begins too: JD(civil date, 12:00) + mean time − Paris longitude."""
    j1, j2 = erfa.cal2jd(year, month, day)  # civil midnight
    noon = float(j1 + j2) + 0.5
    return noon + (mean_time_h - LON_EAST_H[TIME_MERIDIAN]) / 24.0


def transit_residual_s(ut_jd: float, delta_t: float, site: str, ra_h: float) -> float:
    """Local apparent sidereal time at ``site`` minus observed RA, in seconds of time, wrapped
    to ±12 h. Zero (to timing precision) for a meridian observation timed correctly."""
    tt = ut_jd + delta_t / 86400.0
    gast_h = erfa.gst06a(ut_jd, 0.0, tt, 0.0) * 12.0 / math.pi
    lst_h = (gast_h + LON_EAST_H[site]) % 24.0
    return ((lst_h - ra_h + 12.0) % 24.0 - 12.0) * 3600.0


# --- build ---------------------------------------------------------------------------------


def _misprints_by_cell(table: list[dict[str, str]]) -> dict[tuple[str, str, str], dict]:
    return {
        (m["page"], m["line"], m["field"]): m for m in tr.read_misprints(tr.MISPRINT_FILE, table)
    }


def build_rows(
    table: list[dict[str, str]], misprints: dict[tuple[str, str, str], dict[str, str]]
) -> list[dict[str, str]]:
    """The spec-schema rows for every Table 77 observation, in table order."""
    out: list[dict[str, str]] = []
    year = month_s = ""
    for r in table:
        year, month_s = r["year"] or year, r["month"] or month_s
        y, mo, d = int(year), MONTHS[month_s], int(r["day"])
        cell = dict(r)
        applied, logged = [], []
        for f in tr.DATA_FIELDS:
            m = misprints.get((r["page"], r["line"], f))
            if m is None:
                continue
            if m["action"] == "apply":
                cell[f] = m["corrected"]
                applied.append(f"{f}: {m['printed']} -> {m['corrected']}")
            else:
                logged.append(f"misprint (log only) {f}: {m['printed']}")
        ancient = y < FIRST_MODERN_YEAR
        era = "ancient" if ancient else "modern"
        att = attribute(y, mo)
        t_h = parse_sexagesimal(cell["mean_time"])
        ut = epoch_ut_jd(y, mo, d, t_h)
        dt = delta_t_s(y + (mo - 0.5) / 12.0)
        ra_top = parse_sexagesimal(cell["ra_obs"])
        ra_deg = ra_top if ancient else 15.0 * ra_top
        notes = [n for n in (att.note, *logged) if n]
        if att.observatory == "ambiguous":
            resid = ""
            both = {s: transit_residual_s(ut, dt, s, ra_deg / 15.0) for s in ("Greenwich", "Paris")}
            notes.append("transit check (Paris time): "
                         + "; ".join(f"{s} {v:+.1f} s" for s, v in both.items()))  # fmt: skip
        else:
            resid = f"{transit_residual_s(ut, dt, att.observatory, ra_deg / 15.0):+.1f}"
        obs_id = f"LV1846-T77-p{r['page']}-l{int(r['line']):02d}"
        base = {
            "obs_id": obs_id,
            "epoch": f"{year} {month_s} {r['day']}, {r['mean_time']} t.m.",
            "epoch_tt_jd": f"{ut + dt / 86400.0:.8f}",
            "epoch_ut_jd": f"{ut:.8f}",
            "time_meridian": TIME_MERIDIAN,
            "time_reckoning": "astronomical (day begins at mean noon)",
            "observer": att.observer,
            "observatory": att.observatory,
            "instrument": "",
            "frame_as_printed": FRAME_AS_PRINTED,
            "equinox_ref": "of_date",
            "equinox_epoch": "",
            "equinox_citation": EQUINOX_CITATION,
            "precession_constant_arcsec_per_yr": "not_stated",
            "precession_rate_arcsec_per_yr2": "",
            "precession_citation": PRECESSION_CITATION,
            "lv_weight": "",
            "lv_dlon_c_minus_o_as_printed": r["lon_c_minus_o"],
            "lv_dlat_c_minus_o_as_printed": r["lat_c_minus_o"],
            "transit_resid_s": resid,
            "era": era,
            "source_ref": f"{SOURCE}:p{r['page']}:row{r['line']} (Table 77); {att.citation}",
            "misprints_applied": "; ".join(applied),
            "notes": "; ".join(notes),
        }
        cmo = cell["ra_c_minus_o"]
        ra_scale = 1.0 if ancient else 15.0  # printed C − O unit → arcsec of RA
        out.append({
            **base,
            "id": f"{obs_id}-ra",
            "quantity": "ra",
            "value": f"{ra_deg:.9f}",
            "unit": "deg",
            "value_as_printed": r["ra_obs"],
            "unit_as_printed": "deg arcmin arcsec" if ancient else "h m s",
            "lv_residual": f"{-ra_scale * float(cmo):+.2f}" if cmo not in ("", tr.NO_VALUE) else "",
            "lv_residual_unit": "arcsec of RA (O - C)",
            "lv_c_minus_o_as_printed": cmo,
        })  # fmt: skip
        if tr.NO_VALUE not in cell["dec_obs"]:
            dcmo = cell["dec_c_minus_o"]
            out.append({
                **base,
                "id": f"{obs_id}-dec",
                "quantity": "dec",
                "value": f"{parse_sexagesimal(cell['dec_obs']):.9f}",
                "unit": "deg",
                "value_as_printed": r["dec_obs"],
                "unit_as_printed": "deg arcmin arcsec",
                "lv_residual": f"{-float(dcmo):+.2f}" if dcmo not in ("", tr.NO_VALUE) else "",
                "lv_residual_unit": "arcsec (O - C)",
                "lv_c_minus_o_as_printed": dcmo,
            })  # fmt: skip
    return out


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    tmp = path.with_suffix(".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


# --- loader --------------------------------------------------------------------------------

REQUIRED = ("id", "obs_id", "epoch_tt_jd", "quantity", "value", "unit", "equinox_ref",
            "equinox_citation", "precession_constant_arcsec_per_yr", "precession_citation",
            "source_ref")  # fmt: skip


def check_record(r: dict[str, str]) -> None:
    """Raise unless the record carries its equinox reference and the source's precession
    constant, each with a citation (A6b loader rule)."""
    missing = [k for k in REQUIRED if not (r.get(k) or "").strip()]
    if missing:
        raise ValueError(f"{r.get('id', '?')}: missing {missing}")
    ref = r["equinox_ref"]
    if ref == "fixed":
        if not (r.get("equinox_epoch") or "").strip():
            raise ValueError(f"{r['id']}: equinox_ref fixed without equinox_epoch")
        float(r["precession_constant_arcsec_per_yr"])  # must be the source's numeric constant
    elif ref == "of_date":
        c = r["precession_constant_arcsec_per_yr"]
        if c != "not_stated":
            float(c)
    else:
        raise ValueError(f"{r['id']}: equinox_ref {ref!r} is neither of_date nor fixed")


def load(path: Path = OUT_FILE) -> list[dict[str, str]]:
    """Rows of a historical observation file, each checked by :func:`check_record`."""
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        check_record(r)
    return rows


def _years(jd_a: np.ndarray | float, jd_b: np.ndarray | float) -> np.ndarray:
    return (np.asarray(jd_a, dtype=np.float64) - np.asarray(jd_b, dtype=np.float64)) / JULIAN_YEAR_D


def undo_fixed_epoch_reduction(
    lon_fixed_rad: np.ndarray | float,
    t_jd: np.ndarray | float,
    epoch_jd: float,
    constant_arcsec_per_yr: float,
    rate_arcsec_per_yr2: float = 0.0,
    rate_ref_jd: float | None = None,
) -> np.ndarray:
    """Mean longitude of date from one a source pre-reduced to the mean equinox of
    ``epoch_jd``, with *the source's own* annual general precession
    p(τ) = constant + rate·(τ − rate_ref) (years; rate_ref defaults to the epoch).

    Accumulated precession from the epoch to t is ∫ p dτ = constant·Δ + rate·[(t − ref)² −
    (epoch − ref)²]/2. Ecliptic latitude is unchanged by a reduction of this kind."""
    ref = epoch_jd if rate_ref_jd is None else rate_ref_jd
    dt = _years(t_jd, epoch_jd)
    acc = constant_arcsec_per_yr * dt + 0.5 * rate_arcsec_per_yr2 * (
        _years(t_jd, ref) ** 2 - _years(epoch_jd, ref) ** 2
    )
    return np.mod(np.asarray(lon_fixed_rad, dtype=np.float64) + acc / ARCSEC_PER_RAD, 2 * math.pi)


def lon_lat_ecl_j2000(r: dict, lon_rad: float, lat_rad: float) -> tuple[float, float]:
    """Ecliptic J2000 direction of an ecliptic (lon, lat) record, by its equinox metadata:
    ``fixed`` → undo the source's reduction with its own constant (to mean equinox of date),
    then of date → J2000 with IAU 2006 (``frames``). Never IAU 2006 straight from a position
    another precession model produced."""
    check_record(r)
    t = float(r["epoch_tt_jd"])
    if r["equinox_ref"] == "fixed":
        rate = r.get("precession_rate_arcsec_per_yr2") or "0"
        lon_rad = float(undo_fixed_epoch_reduction(
            lon_rad, t, float(r["equinox_epoch"]), float(r["precession_constant_arcsec_per_yr"]),
            float(rate)))  # fmt: skip
    lon, lat = frames.ecl_of_date_to_ecl_j2000(lon_rad, lat_rad, t)
    return float(lon), float(lat)


def main() -> None:
    table = tr.read_pass(tr.RECONCILED_FILE)
    rows = build_rows(table, _misprints_by_cell(table))
    write_csv(OUT_FILE, rows)
    n_obs = len({r["obs_id"] for r in rows})
    print(f"{n_obs} observations, {len(rows)} rows -> {OUT_FILE.relative_to(tr.ROOT)}")
    res = [abs(float(r["transit_resid_s"])) for r in rows if r["transit_resid_s"]]
    print(f"transit check: max |residual| {max(res):.1f} s over {len(res)} rows")


if __name__ == "__main__":
    main()
