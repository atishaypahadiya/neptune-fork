"""Phase 1b step 4: Challis's search as a detection model (docs/notes/P1B_challis_detection.md).

``prerun`` freezes the inputs the plan's § 2 names (Neptune's magnitude, diameter and motion per
night; the geocentric places of the guidance column, S4, S5 and Le Verrier's two places; the
Moon rule per night); ``smoke`` times a small number of accepted draws and checks the closed
form against inner event draws. The full draws, the summary and the figures follow the pre-run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import time
from datetime import date, timedelta
from pathlib import Path

import numpy as np

from fork.historical import challis_1846 as C

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "p1b_challis"
PRERUN_JSON = ROOT / "docs" / "notes" / "P1B_challis_prerun.json"
SMOKE_JSON = ROOT / "docs" / "notes" / "P1B_challis_smoke.json"
PROGRESS_LOG = CACHE / "progress.log"
HIST = ROOT / "data" / "historical" / "adams_challis"
FROZEN = {
    "challis_sweeps.csv": "f9061de02545",
    "neptune_track_1846.csv": "7f688e20b1f7",
    "neptune_on_challis_dates.csv": "815001ff41d6",
}
SEED_DRAWS = (2026, 11, 4, 1)
SEED_EVENTS = (2026, 11, 4, 2)
N_INNER = 1000  # inner event draws per accepted draw (check of the closed form; plan: 200)

# Constants and their sources (§ 2).
V10_NEPTUNE = -6.89  # Mallama & Hilton 2018, Astron. Comput. 25, 10, Eq. 16, t < 1980.0
R_EQ_NEPTUNE_KM = 24764.0  # IAU WGCCRE (Archinal et al. 2018, CeMDA 130, 22), 1-bar equatorial
V10_SOURCE = (
    "Mallama & Hilton 2018, Astronomy and Computing 25, 10, Eq. 16 (t < 1980.0); "
    "fitted to 1954–2014, carried to 1846"
)
R_EQ_SOURCE = "IAU WGCCRE 2015 (Archinal et al. 2018, Celest. Mech. Dyn. Astron. 130, 22)"
AU_KM = 149597870.7  # IAU 2012 Resolution B2
K_GAUSS = 0.01720209895
C_AU_PER_DAY = 173.1446326846693
T0_ISO = "1846-09-23T23:00:00"  # the one 1846 epoch cached for both Neptune and the Earth
SESSION_UT = "23:30:00"
CAMBRIDGE = (52.2139, 0.0956)  # Cambridge Observatory, degrees north and east

# Adams's ephemeris for 325° (Challis, Memoirs RAS 16, 421): RA minutes past 21h, Dec arcmin.
ADAMS_325 = {
    "1846-07-20": (51.3, -(13 * 60 + 0)),
    "1846-08-09": (49.5, -(13 * 60 + 9)),
    "1846-08-29": (47.9, -(13 * 60 + 18)),
    "1846-09-18": (46.5, -(13 * 60 + 26)),
    "1846-10-08": (45.1, -(13 * 60 + 33)),
}
# Le Verrier, Comptes rendus 23 (séance du 31 août 1846), 428–438; elements at 432, limits at
# 434–436 (archive.org comptesrendusheb23acad, PDF page = printed + 10).
LV_AUG = {
    "a_au": 36.154, "e": 0.10761, "varpi_deg": 284.0 + 45.0 / 60.0,
    "mean_lon_deg": 318.0 + 47.0 / 60.0, "mass_msun": 1.0 / 9300.0,
    "true_lon_1847_deg": 326.0 + 32.0 / 60.0, "r_au": 33.06,
    "limits_helio_lon_deg": (321.0, 335.0), "disc_arcsec": 3.3,
    "epoch_iso": "1847-01-01T00:00:00",
}  # fmt: skip
LV_JUNE_LIMITS_DEG = (315.0, 335.0)  # 325° "not an error of 10 degrees" (CR 22, 917)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    with PROGRESS_LOG.open("a") as f:
        f.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}\n")


def _dump(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")


def check_frozen() -> dict:
    out = {}
    for name, head in FROZEN.items():
        h = sha256(HIST / name)
        if not h.startswith(head):
            raise RuntimeError(f"{name}: sha256 {h[:12]} is not the plan's {head}")
        out[name] = h
    return out


def night_dates() -> list[str]:
    d0 = date.fromisoformat(C.NIGHT0)
    return [(d0 + timedelta(days=i)).isoformat() for i in range(C.N_NIGHTS)]


def _states_at_t0():
    from fork.data import horizons as h
    from fork.data.frames import jd_tdb_from_iso

    jd0 = float(jd_tdb_from_iso(T0_ISO))
    out = {}
    for body in ("neptune", "earth"):
        jd, st = h._load_cache(h._cache_path(body, h.cache_dir()))
        i = int(np.argmin(np.abs(jd - jd0)))
        if abs(jd[i] - jd0) > 1e-3:
            raise RuntimeError(f"{body}: no cached state at {T0_ISO}")
        out[body] = (st[i, :3].copy(), st[i, 3:].copy())
    return jd0, out


def _drift(x: np.ndarray, v: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    from fork.integrators.kepler import kepler_drift_inplace

    a, b = x.copy()[None, :], v.copy()[None, :]
    kepler_drift_inplace(a, b, K_GAUSS * K_GAUSS, dt)
    return a[0], b[0]


def adams_325(day: str) -> tuple[float, float]:
    """The guidance column interpolated to ``day``: RA hours, Dec degrees (equinox of date)."""
    xs = [date.fromisoformat(k).toordinal() for k in ADAMS_325]
    o = date.fromisoformat(day).toordinal()
    ra = float(np.interp(o, xs, [v[0] for v in ADAMS_325.values()]))
    de = float(np.interp(o, xs, [v[1] for v in ADAMS_325.values()]))
    return 21.0 + ra / 60.0, de / 60.0


def moon_table(days: list[str], ra_h: list[float], dec_deg: list[float]) -> list[dict]:
    """Per night: the Moon's lit fraction, its separation from Neptune's place and whether it
    is above the horizon at Cambridge at 21:30, 23:30 or 01:30 UT; the plan's rule excludes a
    night when it is up, more than half lit and within 30° of the zone."""
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body
    from astropy.time import Time

    loc = EarthLocation(lat=CAMBRIDGE[0] * u.deg, lon=CAMBRIDGE[1] * u.deg)
    out = []
    for day, ra, dec in zip(days, ra_h, dec_deg, strict=True):
        t = Time(f"{day}T{SESSION_UT}", scale="utc")
        moon, sun = get_body("moon", t), get_body("sun", t)
        lit = float((1.0 - math.cos(moon.separation(sun).rad)) / 2.0)
        zone = SkyCoord(ra=ra * 15.0 * u.deg, dec=dec * u.deg, frame="fk5", equinox=t)
        sep = float(moon.separation(zone).deg)
        alts = []
        for dh in (-2.0, 0.0, 2.0):
            tt = t + dh * u.hour
            alts.append(float(get_body("moon", tt, loc).transform_to(
                AltAz(obstime=tt, location=loc)).alt.deg))  # fmt: skip
        up = max(alts) > 0.0
        out.append({"date": day, "lit_fraction": round(lit, 3), "separation_deg": round(sep, 1),
                    "max_altitude_deg": round(max(alts), 1), "up": up,
                    "excluded": bool(up and lit > 0.5 and sep < 30.0)})  # fmt: skip
    return out


def prerun() -> dict:
    import warnings

    warnings.simplefilter("ignore")  # astropy's notices for dates before 1900
    from fork.data.frames import jd_tdb_from_iso, precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.export import p1b_adams as P
    from fork.model.bodies import M_CENTRAL_MSUN, NEPTUNE, G

    hashes = check_frozen()
    jd0, st = _states_at_t0()
    # The printed places lie in the ecliptic (latitude 0), as Adams's 325° column does: the
    # frame rotation of step 3 (Uranus's plane) is the identity here.
    inc_raan = (0.0, 0.0)
    tg = P.targets()
    src = "Comptes rendus 23, 432 (31 Aug 1846), archive.org comptesrendusheb23acad (P)"
    orbits = {
        "S4": tg["S4"]["orbit"],
        "S5": tg["S5"]["orbit"],
        "LVJ": tg["LV-June"]["orbit"],
        "LVA": P._published("Le Verrier, 31 August 1846 (Comptes rendus 23, 432)",
                            LV_AUG["a_au"], LV_AUG["e"], LV_AUG["mass_msun"], LV_AUG["varpi_deg"],
                            LV_AUG["mean_lon_deg"], LV_AUG["epoch_iso"], (src,)),
    }  # fmt: skip
    days = night_dates()
    nights = []
    for day in days:
        jd = float(jd_tdb_from_iso(f"{day}T{SESSION_UT}"))
        xn, vn = _drift(*st["neptune"], jd - jd0)
        xe, _ = _drift(*st["earth"], jd - jd0)
        geo = P.EpochGeometry(day, jd, xe, np.asarray(
            precession_matrix_ecl_j2000_to_ecl_of_date(jd)), P.mean_obliquity_deg(jd))  # fmt: skip
        mu_n = G * (M_CENTRAL_MSUN + NEPTUNE.mass_msun)
        nep = geo.derived(xn, vn, mu_n)
        r, delta = float(np.linalg.norm(xn)), float(np.linalg.norm(xn - xe))
        row = {
            "date": day,
            "neptune": {"ra_h": nep["ra_h"], "dec_deg": nep["dec_deg"],
                        "helio_true_deg": nep["helio_true"], "geo_lon_deg": nep["geo_lon"],
                        "r_au": r, "delta_au": delta,
                        "V": V10_NEPTUNE + 5.0 * math.log10(r * delta),
                        "diameter_arcsec": math.degrees(2.0 * R_EQ_NEPTUNE_KM / (delta * AU_KM))
                        * 3600.0},
            "targets": {"col325": dict(zip(("ra_h", "dec_deg"), adams_325(day), strict=True))},
        }  # fmt: skip
        for name, sol in orbits.items():
            mu = G * (M_CENTRAL_MSUN + sol.mass_msun)
            d = geo.derived(*P.published_state_ecl_j2000(sol, jd, inc_raan), mu)
            row["targets"][name] = {
                "ra_h": d["ra_h"], "dec_deg": d["dec_deg"],
                "helio_true_deg": d["helio_true"], "geo_lon_deg": d["geo_lon"],
            }  # fmt: skip
        nights.append(row)
    moon = moon_table(days, [n["neptune"]["ra_h"] for n in nights],
                      [n["neptune"]["dec_deg"] for n in nights])  # fmt: skip
    for n, m in zip(nights, moon, strict=True):
        n["moon"] = m
    # the equinox-1800 right ascension, for the Berlin sheet (N1)
    import astropy.units as u
    from astropy.coordinates import FK5, SkyCoord
    from astropy.time import Time

    for n in nights:
        t = Time(f"{n['date']}T{SESSION_UT}", scale="utc")
        c = SkyCoord(ra=n["neptune"]["ra_h"] * 15.0 * u.deg, dec=n["neptune"]["dec_deg"] * u.deg,
                     frame=FK5(equinox=t))  # fmt: skip
        c = c.transform_to(FK5(equinox=Time("B1800", scale="tt")))
        n["neptune"]["ra1800_h"] = float(c.ra.hour)
        n["neptune"]["dec1800_deg"] = float(c.dec.deg)
    idx = {d: i for i, d in enumerate(days)}
    known = ("1846-07-29", "1846-07-30", "1846-08-04", "1846-08-12", "1846-09-29")
    motion = []
    for a, b in itertools.pairwise(known):
        na, nb = nights[idx[a]]["neptune"], nights[idx[b]]["neptune"]
        motion.append({"from": a, "to": b,
                       "d_ra_s": (nb["ra_h"] - na["ra_h"]) * 3600.0,
                       "d_dec_arcmin": (nb["dec_deg"] - na["dec_deg"]) * 60.0})  # fmt: skip
    # check against step 2's table (mean of date at the session times it used)
    check = []
    with (HIST / "neptune_track_1846.csv").open() as f:
        track = {r["date"]: r for r in csv.DictReader(f)}
    for d in known:
        s = track[d]["ra_mean_date"].replace("h", "").replace("m", "").replace("s", "").split()
        ra = float(s[0]) + float(s[1]) / 60.0 + float(s[2]) / 3600.0
        check.append({"date": d, "ra_minus_step2_s": (nights[idx[d]]["neptune"]["ra_h"] - ra)
                      * 3600.0})  # fmt: skip
    out = {
        "plan": "docs/notes/P1B_challis_detection.md",
        "frozen_inputs_sha256": hashes,
        "constants": {
            "V10_neptune": [V10_NEPTUNE, V10_SOURCE],
            "equatorial_radius_km": [R_EQ_NEPTUNE_KM, R_EQ_SOURCE],
            "au_km": [AU_KM, "IAU 2012 Resolution B2"],
            "session_time_ut": SESSION_UT,
            "observatory_lat_lon_deg": CAMBRIDGE,
            "ephemeris": f"two-body arcs from the cached Horizons states at {T0_ISO} TDB",
        },
        "le_verrier_august": LV_AUG,
        "le_verrier_june_limits_deg": LV_JUNE_LIMITS_DEG,
        "nights": nights,
        "motion_between_his_nights": motion,
        "check_against_step2_track": check,
    }
    _dump(PRERUN_JSON, out)
    _log(f"prerun written, sha256 {sha256(PRERUN_JSON)[:16]}")
    return out


# Berlin sheet Hora XXII, N1's fallbacks: 21h 56m to 23h 04m (Krajnović 2021), read against the
# right ascension for the equinox of 1800 (first fallback) or of date (second fallback).
HORA_XXII_RA_H = (21.0 + 56.0 / 60.0, 23.0 + 4.0 / 60.0)


def geometry(pre: dict) -> dict:
    """The per-night arrays the model reads, from the frozen pre-run."""
    n = pre["nights"]
    ra = np.array([x["neptune"]["ra_h"] for x in n])
    ra1800 = np.array([x["neptune"]["ra1800_h"] for x in n])
    lo, hi = HORA_XXII_RA_H
    return {
        "ra_h": ra,
        "dec_deg": np.array([x["neptune"]["dec_deg"] for x in n]),
        "moon_excluded": np.array([x["moon"]["excluded"] for x in n]),
        "on_chart_1800": (ra1800 >= lo) & (ra1800 <= hi),
        "on_chart_date": (ra >= lo) & (ra <= hi),
        "targets": {k: {"ra_h": np.array([x["targets"][k]["ra_h"] for x in n]),
                        "dec_deg": np.array([x["targets"][k]["dec_deg"] for x in n])}
                    for k in n[0]["targets"]},
    }  # fmt: skip


def smoke(n_accept: int = 1000, n_inner: int = N_INNER) -> dict:
    """Timed smoke run: draws until ``n_accept`` are accepted, every case in closed form, and
    the same events by ``n_inner`` inner draws as a check; the Monte Carlo share of the
    per-draw spread that inner draws would carry."""
    pre = json.loads(PRERUN_JSON.read_text())
    geo = geometry(pre)
    rng = np.random.default_rng(list(SEED_DRAWS))
    ev = {c: np.random.default_rng([*SEED_EVENTS, i]) for i, c in enumerate(C.CASES)}
    t0, c0 = time.time(), time.process_time()
    n_drawn = n_c23 = n_c4 = 0
    acc = []
    while len(acc) < n_accept:
        d = C.draw(rng, geo)
        a = C.accept(d, rng)
        n_drawn += 1
        n_c23 += a["c23"]
        n_c4 += a["c4"]
        if a["all"]:
            acc.append(d)
    t_draw = time.time() - t0
    exact = {c: np.zeros(n_accept) for c in C.CASES}
    second = {c: np.zeros(n_accept) for c in C.CASES}
    t2 = {c: np.zeros(n_accept) for c in C.CASES}
    t1 = time.time()
    for i, d in enumerate(acc):
        for c in C.CASES:
            f1, f2 = C.cdfs(d, C.groups(d, geo, c))
            exact[c][i] = f1[C.PRIMARY_LAST]
            second[c][i] = f1[C.SECONDARY_LAST]
            t2[c][i] = f2[C.PRIMARY_LAST]
    t_exact = time.time() - t1
    t1 = time.time()
    inner = {c: np.zeros(n_accept) for c in C.CASES}
    for i, d in enumerate(acc):
        for c in C.CASES:
            s = C.sample(d, C.groups(d, geo, c), ev[c], n_inner)
            inner[c][i] = float(np.mean(s <= C.PRIMARY_LAST))
    t_inner = time.time() - t1
    cases = {}
    for c in C.CASES:
        p = exact[c]
        between = float(np.var(p))
        mc = {str(m): float(np.mean(p * (1 - p)) / m) for m in (200, 1000)}
        se = np.sqrt(p * (1 - p) / n_inner) + 1.0 / n_inner  # one count of slack
        cases[c] = {
            "primary_mean": float(p.mean()),
            "primary_q05_q95": [float(x) for x in np.quantile(p, [0.05, 0.95])],
            "secondary_mean": float(second[c].mean()),
            "t2_primary_mean": float(t2[c].mean()),
            "inner_minus_exact_mean": float(np.mean(inner[c] - p)),
            "inner_minus_exact_max_abs_in_se": float(np.max(np.abs(inner[c] - p) / se)),
            "variance_between_draws": between,
            "mc_share_of_spread": {m: (v / (v + between) if v + between > 0 else 0.0)
                                   for m, v in mc.items()},
        }  # fmt: skip
    out = {
        "n_accepted": n_accept, "n_drawn": n_drawn, "n_inner": n_inner,
        "acceptance": {"all": n_accept / n_drawn, "conditions_2_3": n_c23 / n_drawn,
                       "condition_4": n_c4 / n_drawn},
        "seconds": {"draws_and_abduction": t_draw, "closed_form_all_cases": t_exact,
                    "inner_draws_all_cases": t_inner, "cpu_total": time.process_time() - c0},
        "estimate_full_run_seconds": {
            "closed_form_100000_accepted": (t_draw + t_exact) * 100_000 / n_accept,
            "with_inner_draws": (t_draw + t_exact + t_inner) * 100_000 / n_accept},
        "cases": cases,
    }  # fmt: skip
    _dump(SMOKE_JSON, out)
    _log(f"smoke: {n_accept} accepted of {n_drawn}, closed form {t_exact:.1f} s, inner "
         f"{t_inner:.1f} s")  # fmt: skip
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prerun")
    sm = sub.add_parser("smoke")
    sm.add_argument("--n-accept", type=int, default=1000)
    sm.add_argument("--n-inner", type=int, default=N_INNER)
    a = ap.parse_args(argv)
    if a.cmd == "prerun":
        prerun()
    elif a.cmd == "smoke":
        smoke(a.n_accept, a.n_inner)


if __name__ == "__main__":
    main()
