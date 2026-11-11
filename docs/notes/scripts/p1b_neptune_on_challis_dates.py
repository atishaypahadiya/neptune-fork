"""Neptune's place on Challis's observing dates (Phase 1b step 2), from the cached states only.

Computation, not transcription. The Sun-centred states of Neptune and the Earth at the one cached
1846 epoch common to both (1846-09-23 23:00 TDB) are carried to each date as two-body Kepler arcs
(fork.integrators.kepler.kepler_drift_inplace, mu = k^2); Jupiter's and Saturn's pull on either body
over the 90 days spanned shifts Neptune's geocentric direction by well under 1". No Horizons
fetch is made (FORK_HORIZONS_CACHE_ONLY=1).

Outputs (data/historical/adams_challis/):
  neptune_on_challis_dates.csv  the sessions and Challis's three reduced places
  neptune_track_1846.csv        daily 1 July - 15 October 1846: mean place of date, equinox-1800
                                place, Adams's 325-degree ephemeris, and the Berlin chart tests

Frames: "mean of date" is the mean equator and equinox of the date (FK5 precession in astropy);
"apparent" adds light time, annual aberration and nutation (astropy TETE); "1800" is the mean
equinox of B1800, the epoch of the Berlin Akademische Sternkarten. Times: the sessions are taken at
23:30 UT (Cambridge night); Challis's places are at his Greenwich mean astronomical times (day
beginning at noon), converted to civil time.
"""

from __future__ import annotations

import csv
import os
from datetime import date
from pathlib import Path

os.environ["FORK_HORIZONS_CACHE_ONLY"] = "1"

import astropy.units as u
import numpy as np
from astropy.coordinates import FK5, GCRS, TETE, CartesianRepresentation, SkyCoord
from astropy.time import Time

from fork.data import horizons as h
from fork.integrators.kepler import kepler_drift_inplace

K = 0.01720209895
MU = K * K
C_AU_PER_DAY = 173.1446326846693
EPS = np.radians(23.4392911)  # obliquity J2000: ecliptic J2000 -> equator J2000
T0 = Time("1846-09-23T23:00:00", scale="tdb")
OUT = Path(__file__).resolve().parents[3] / "data" / "historical" / "adams_challis"

# Adams's ephemeris for 325 degrees of heliocentric longitude (Challis, Memoirs RAS 16, 421):
# RA minutes past 21h, declination in arcminutes.
ADAMS_325 = {
    "1846-07-20": (51.3, -(13 * 60 + 0)),
    "1846-08-09": (49.5, -(13 * 60 + 9)),
    "1846-08-29": (47.9, -(13 * 60 + 18)),
    "1846-09-18": (46.5, -(13 * 60 + 26)),
    "1846-10-08": (45.1, -(13 * 60 + 33)),
}


def _state(body: str) -> tuple[float, np.ndarray, np.ndarray]:
    jd, st = h._load_cache(h._cache_path(body, h.cache_dir()))
    i = int(np.argmin(np.abs(jd - T0.jd)))
    if abs(jd[i] - T0.jd) > 1e-3:
        raise RuntimeError(f"{body}: no cached state at {T0.iso}")
    return float(jd[i]), st[i, :3].copy(), st[i, 3:].copy()


JD0, XN, VN = _state("neptune")
_, XE, VE = _state("earth")


def _drift(x: np.ndarray, v: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    a, b = x.copy()[None, :], v.copy()[None, :]
    kepler_drift_inplace(a, b, MU, dt)
    return a[0], b[0]


def _ecl_to_eq(p: np.ndarray) -> np.ndarray:
    return np.array(
        [p[0], p[1] * np.cos(EPS) - p[2] * np.sin(EPS), p[1] * np.sin(EPS) + p[2] * np.cos(EPS)]
    )


def _fmt_ra(c) -> str:
    return c.ra.to_string(unit=u.hour, sep=("h ", "m ", "s"), precision=1, pad=True)


def _fmt_dec(c) -> str:
    return c.dec.to_string(
        unit=u.deg, sep=("° ", "′ ", "″"), precision=0, alwayssign=True, pad=True
    )


def place(t: Time) -> dict:
    """Neptune's geocentric place at ``t``: mean of date, apparent, equinox 1800, J2000."""
    dt = t.tdb.jd - JD0
    xn, _ = _drift(XN, VN, dt)
    xe, ve = _drift(XE, VE, dt)
    g = xn - xe
    dist = float(np.linalg.norm(g))
    # light time and annual aberration for the apparent place
    tau = dist / C_AU_PER_DAY
    xn_lt, _ = _drift(XN, VN, dt - tau)
    g_lt = xn_lt - xe
    u_ast = g_lt / np.linalg.norm(g_lt)
    u_app = u_ast + ve / C_AU_PER_DAY
    u_app /= np.linalg.norm(u_app)

    icrs = SkyCoord(CartesianRepresentation(_ecl_to_eq(g) * u.au), frame="icrs")
    mean_date = icrs.transform_to(FK5(equinox=t))
    eq1800 = icrs.transform_to(FK5(equinox=Time("B1800", scale="tt")))
    gcrs = SkyCoord(CartesianRepresentation(_ecl_to_eq(u_app)), frame=GCRS(obstime=t))
    app = gcrs.transform_to(TETE(obstime=t))
    return {
        "ra_mean_date": _fmt_ra(mean_date),
        "dec_mean_date": _fmt_dec(mean_date),
        "ra_apparent": _fmt_ra(app),
        "dec_apparent": _fmt_dec(app),
        "ra_1800": _fmt_ra(eq1800),
        "dec_1800": _fmt_dec(eq1800),
        "ra_J2000": _fmt_ra(icrs),
        "dec_J2000": _fmt_dec(icrs),
        "distance_au": round(dist, 3),
        "_ra_date_h": float(mean_date.ra.hour),
        "_dec_date_deg": float(mean_date.dec.deg),
        "_ra_1800_h": float(eq1800.ra.hour),
        "_ra_app_h": float(app.ra.hour),
        "_dec_app_deg": float(app.dec.deg),
    }


def adams_325(day: str) -> tuple[float, float]:
    xs = [date.fromisoformat(k).toordinal() for k in ADAMS_325]
    o = date.fromisoformat(day).toordinal()
    ra = float(np.interp(o, xs, [v[0] for v in ADAMS_325.values()]))
    de = float(np.interp(o, xs, [v[1] for v in ADAMS_325.values()]))
    return ra, de


def _hms_to_h(s: str) -> float:
    hh, mm, ss = s.split()
    return float(hh) + float(mm) / 60 + float(ss) / 3600


def _npd_to_deg(s: str) -> float:
    dd, mm, ss = s.split()
    return float(dd) + float(mm) / 60 + float(ss) / 3600


SESSIONS = [
    ("1846-07-29", "1846-07-29T23:30:00", "session (trial, telescope fixed)", None, None),
    ("1846-07-30", "1846-07-30T23:30:00", "session (telescope moving, zone 9')", None, None),
    (
        "1846-08-04",
        "1846-08-05T01:36:25",
        "session (sector, reference stars); Challis's place",
        "21 58 14.70",
        "102 57 32.2",
    ),
    (
        "1846-08-12",
        "1846-08-13T01:03:26",
        "session (telescope fixed, zone 9'); Challis's place",
        "21 57 26.13",
        "103 2 0.2",
    ),
    ("1846-09-23", "1846-09-23T23:30:00", "Berlin discovery (for reference)", None, None),
    (
        "1846-09-29",
        "1846-09-29T22:10:00",
        "session (breadth in declination, Le Verrier's limits); Challis's place",
        "21 52 48.9",
        "103 26 30",
    ),
    ("1846-09-30", "1846-09-30T23:30:00", "no observations (moon)", None, None),
    ("1846-10-01", "1846-10-01T23:30:00", "news of the discovery; search ended", None, None),
]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for day, iso, what, ra_ch, npd_ch in SESSIONS:
        t = Time(iso, scale="utc")
        p = place(t)
        ar, ad = adams_325(day)
        row = {
            "date": day,
            "time_ut": iso[11:16],
            "event": what,
            **{k: v for k, v in p.items() if not k.startswith("_")},
            "adams_325_ra": f"21h {ar:.1f}m",
            "adams_325_dec": f"-{int(-ad // 60)}° {(-ad) % 60:.0f}′",
            "neptune_minus_adams_ra_min_of_time": round(p["_ra_date_h"] * 60 - 21 * 60 - ar, 1),
            "neptune_minus_adams_dec_arcmin": round(p["_dec_date_deg"] * 60 - ad, 1),
            "challis_ra": ra_ch or "",
            "challis_npd": npd_ch or "",
            "challis_minus_apparent_ra_s": "",
            "challis_minus_apparent_npd_arcsec": "",
        }
        if ra_ch:
            dra = (_hms_to_h(ra_ch) - p["_ra_app_h"]) * 3600
            dnpd = (_npd_to_deg(npd_ch) - (90 - p["_dec_app_deg"])) * 3600
            row["challis_minus_apparent_ra_s"] = round(dra, 1)
            row["challis_minus_apparent_npd_arcsec"] = round(dnpd, 0)
        rows.append(row)
    with open(OUT / "neptune_on_challis_dates.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    track = []
    jd = Time("1846-07-01", scale="utc").jd
    end = Time("1846-10-15", scale="utc").jd
    while jd <= end + 0.5:
        t = Time(jd, format="jd", scale="utc")
        day = t.iso[:10]
        p = place(t)
        ar, ad = adams_325(day)
        r1800 = p["_ra_1800_h"]
        track.append(
            {
                "date": day,
                "ra_mean_date": p["ra_mean_date"],
                "dec_mean_date": p["dec_mean_date"],
                "ra_1800": p["ra_1800"],
                "dec_1800": p["dec_1800"],
                "adams_325_ra": f"21h {ar:.1f}m",
                "adams_325_dec": f"-{int(-ad // 60)}° {(-ad) % 60:.0f}′",
                "neptune_minus_adams_ra_min_of_time": round(p["_ra_date_h"] * 60 - 21 * 60 - ar, 1),
                "neptune_minus_adams_dec_arcmin": round(p["_dec_date_deg"] * 60 - ad, 1),
                "on_hora_xxii_epoch1800_4min_overlap": r1800 >= 21 + 56 / 60,
                "on_hora_xxii_epoch1800_no_overlap": r1800 >= 22.0,
                "on_hora_xxi_epoch1800_4min_overlap": r1800 <= 22 + 4 / 60,
                "on_hora_xxii_if_limits_applied_to_ra_of_date": p["_ra_date_h"] >= 21 + 56 / 60,
            }
        )
        jd += 1.0
    with open(OUT / "neptune_track_1846.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(track[0].keys()))
        w.writeheader()
        w.writerows(track)
    for r in rows:
        print(
            r["date"],
            r["ra_apparent"],
            r["dec_apparent"],
            "| 1800:",
            r["ra_1800"],
            r["dec_1800"],
            "| Challis - apparent:",
            r["challis_minus_apparent_ra_s"],
            r["challis_minus_apparent_npd_arcsec"],
        )


if __name__ == "__main__":
    main()
