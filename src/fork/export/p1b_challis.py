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
# The frozen smoke record in docs/notes (P1B_challis_smoke.json, timings in its sidecar) is the
# output of the pre-run's code; a rerun writes to the cache and leaves the record alone.
SMOKE_JSON = CACHE / "smoke.json"
SMOKE_TIMINGS_JSON = CACHE / "smoke_timings.json"
RESULTS_JSON = ROOT / "docs" / "notes" / "P1B_challis_detection.json"
FIGURES = {
    "curves": ROOT / "docs" / "notes" / "P1B_challis_detection.png",
    "posterior": ROOT / "docs" / "notes" / "P1B_challis_detection_posterior.png",
    "tornado": ROOT / "docs" / "notes" / "P1B_challis_detection_tornado.png",
}
N_ACCEPT = 100_000
N_ACCEPT_HALVED = 20_000  # the prior-width check of § 6.5
N_TORNADO = 10_000  # the first accepted draws, for the tornado of § 6.4
CHUNK = 5_000
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
        "cases": cases,
    }  # fmt: skip
    timings = {
        "seconds": {"draws_and_abduction": t_draw, "closed_form_all_cases": t_exact,
                    "inner_draws_all_cases": t_inner, "cpu_total": time.process_time() - c0},
        "estimate_full_run_seconds": {
            "closed_form_100000_accepted": (t_draw + t_exact) * 100_000 / n_accept,
            "with_inner_draws": (t_draw + t_exact + t_inner) * 100_000 / n_accept},
    }  # fmt: skip
    _dump(SMOKE_JSON, out)
    _dump(SMOKE_TIMINGS_JSON, timings)
    _log(f"smoke: {n_accept} accepted of {n_drawn}, closed form {t_exact:.1f} s, inner "
         f"{t_inner:.1f} s")  # fmt: skip
    return out


READINGS = ("held", "replay")
PARAMS = ("p_disc_routine", "p_disc_looking", "p_follow", "p_recorded", "stars_per_hour",
          "moving_extent_factor", "weather", "desk_hours_per_pair", "desk_hours_per_day",
          "n_unknown_drawn", "n_unknown_sessions", "hours_per_session_mean", "hours_total",
          "stars_total", "extent_jul30_deg", "extent_aug12_deg", "c4b_breadth_arcmin",
          "branch_ra", "north_first")  # fmt: skip
# the tornado's inputs: set on the draw (attribute, low and high from the posterior) or, for
# the inputs that shape the sessions, read by conditioning on the lowest and highest tenth
SETTABLE = {"p_disc_routine": "p_r", "p_disc_looking": "p_l", "p_follow": "p_f",
            "p_recorded": "p_rec", "stars_per_hour": "rate", "moving_extent_factor": "f_mov",
            "desk_hours_per_pair": "desk_pair", "desk_hours_per_day": "desk_day",
            "c4b_breadth_arcmin": "c4b_breadth"}  # fmt: skip
CONDITIONED = ("weather", "n_unknown_sessions", "hours_per_session_mean")
SWITCHES = {"north_first": "north_first", "branch_ra": "branch"}


def params_of(d: C.Draw) -> list[float]:
    h = [v for k, v in d.hours.items() if k != C.SEP29]
    return [d.p_r, d.p_l, d.p_f, d.p_rec, d.rate, d.f_mov, d.weather, d.desk_pair, d.desk_day,
            d.n_unknown, len(d.unknown), float(np.mean(h)), float(np.sum(h)) + C.SEP29_HOURS,
            C.total_stars(d), d.hours[C.JUL30] * d.f_mov * 15.0, d.hours[C.AUG12] * 15.0,
            d.c4b_breadth, float(d.branch == "ra"), float(d.north_first)]  # fmt: skip


def accepted_draws(geo: dict, n_accept: int, priors: dict | None = None):
    """The draws and the abduction, from the plan's seed: the accepted draws, the parameters
    of every draw and the counts per condition."""
    rng = np.random.default_rng(list(SEED_DRAWS))
    acc, prior_rows = [], []
    n_c23 = n_c4 = n = 0
    while len(acc) < n_accept:
        d = C.draw(rng, geo, priors)
        a = C.accept(d, rng)
        n += 1
        n_c23 += a["c23"]
        n_c4 += a["c4"]
        prior_rows.append(params_of(d))
        if a["all"]:
            acc.append(d)
    return acc, np.array(prior_rows), {"n_drawn": n, "n_accepted": n_accept,
                                       "all": n_accept / n, "conditions_2_3": n_c23 / n,
                                       "condition_4": n_c4 / n}  # fmt: skip


def _chunk(acc: list, geo: dict, k: int, n_inner: int) -> dict:
    """One chunk of accepted draws through every case and both readings."""
    n, nc = len(acc), len(C.CASES)
    t1 = np.zeros((2, nc, n, 2))  # reading, case, draw, (primary, secondary)
    t2 = np.zeros((2, nc, n, 2))
    curve = np.zeros((nc, n, C.N_NIGHTS), dtype=np.float32)  # held reading, T1 by night
    curve_replay_sum = np.zeros((nc, C.N_NIGHTS))
    inner_sum = np.zeros((2, nc))
    inner_maxz = np.zeros((2, nc))
    ev = [np.random.default_rng([*SEED_EVENTS, ci, k]) for ci in range(nc)]
    for i, d in enumerate(acc):
        for ci, case in enumerate(C.CASES):
            for ri, name in enumerate(READINGS):
                gs = C.groups(d, geo, case, held=name == "held")
                f1, f2 = C.cdfs(d, gs)
                t1[ri, ci, i] = f1[C.PRIMARY_LAST], f1[C.SECONDARY_LAST]
                t2[ri, ci, i] = f2[C.PRIMARY_LAST], f2[C.SECONDARY_LAST]
                if ri == 0:
                    curve[ci, i] = f1
                else:
                    curve_replay_sum[ci] += f1
                if n_inner:
                    p = f1[C.PRIMARY_LAST]
                    s = float(np.mean(C.sample(d, gs, ev[ci], n_inner) <= C.PRIMARY_LAST))
                    inner_sum[ri, ci] += s - p
                    z = abs(s - p) / (math.sqrt(p * (1 - p) / n_inner) + 1.0 / n_inner)
                    inner_maxz[ri, ci] = max(inner_maxz[ri, ci], z)
    return {"t1": t1, "t2": t2, "curve": curve, "curve_replay_sum": curve_replay_sum,
            "inner_sum": inner_sum, "inner_maxz": inner_maxz,
            "params": np.array([params_of(d) for d in acc])}  # fmt: skip


def _run(tag: str, geo: dict, n_accept: int, priors: dict | None, n_inner: int) -> list:
    out = CACHE / tag
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    acc, prior_rows, counts = accepted_draws(geo, n_accept, priors)
    np.savez_compressed(out / "prior.npz", params=prior_rows)
    _dump(out / "acceptance.json", counts)
    _log(f"{tag}: {n_accept} accepted of {counts['n_drawn']} in {time.time() - t0:.0f} s")
    for k in range(math.ceil(n_accept / CHUNK)):
        f = out / f"chunk_{k:03d}.npz"
        if f.exists():
            continue
        r = _chunk(acc[k * CHUNK : (k + 1) * CHUNK], geo, k, n_inner)
        tmp = out / f"chunk_{k:03d}.tmp.npz"
        np.savez_compressed(tmp, **r)
        tmp.rename(f)
        _log(f"{tag}: chunk {k + 1} of {math.ceil(n_accept / CHUNK)} at {time.time() - t0:.0f} s")
    return acc


def _tornado(acc: list, geo: dict) -> dict:
    """§ 6.4 on the held reading: each input at its posterior 5 % and 95 % values in turn, the
    others left to vary; the change in the mean primary probability per case."""
    import dataclasses

    par = np.array([params_of(d) for d in acc])
    col = {name: par[:, i] for i, name in enumerate(PARAMS)}

    def mean_primary(draws: list) -> np.ndarray:
        out = np.zeros(len(C.CASES))
        for d in draws:
            for ci, case in enumerate(C.CASES):
                out[ci] += C.cdfs(d, C.groups(d, geo, case, held=True))[0][C.PRIMARY_LAST]
        return out / len(draws)

    res = {}
    for name, attr in SETTABLE.items():
        lo, hi = np.quantile(col[name], [0.05, 0.95])
        a = mean_primary([dataclasses.replace(d, **{attr: float(lo)}) for d in acc])
        b = mean_primary([dataclasses.replace(d, **{attr: float(hi)}) for d in acc])
        res[name] = {"how": "set", "low": float(lo), "high": float(hi), "p_low": a.tolist(),
                     "p_high": b.tolist()}  # fmt: skip
        _log(f"tornado: {name}")
    for name, attr in SWITCHES.items():
        vals = (False, True) if attr == "north_first" else ("dec", "ra")
        a = mean_primary([dataclasses.replace(d, **{attr: vals[0]}) for d in acc
                          if attr == "north_first" or d.branch == vals[0]])  # fmt: skip
        b = mean_primary([dataclasses.replace(d, **{attr: vals[1]}) for d in acc
                          if attr == "north_first" or d.branch == vals[1]])  # fmt: skip
        res[name] = {"how": "set" if attr == "north_first" else "conditioned", "low": 0.0,
                     "high": 1.0, "p_low": a.tolist(), "p_high": b.tolist()}  # fmt: skip
    base = np.zeros((len(acc), len(C.CASES)))
    for i, d in enumerate(acc):
        for ci, case in enumerate(C.CASES):
            base[i, ci] = C.cdfs(d, C.groups(d, geo, case, held=True))[0][C.PRIMARY_LAST]
    for name in CONDITIONED:
        lo, hi = np.quantile(col[name], [0.1, 0.9])
        res[name] = {"how": "conditioned", "low": float(lo), "high": float(hi),
                     "p_low": base[col[name] <= lo].mean(axis=0).tolist(),
                     "p_high": base[col[name] >= hi].mean(axis=0).tolist()}  # fmt: skip
    return {"n_draws": len(acc), "cases": list(C.CASES), "inputs": res}


def draws(n_accept: int = N_ACCEPT, n_inner: int = N_INNER) -> None:
    """The full draws (resumable per chunk), the prior-width run and the tornado."""
    from fork import ledger

    geo = geometry(json.loads(PRERUN_JSON.read_text()))
    t0, c0 = time.time(), time.process_time()
    acc = _run("main", geo, n_accept, None, n_inner)
    _run("halved", geo, min(N_ACCEPT_HALVED, n_accept), C.halved(C.PRIORS), 0)
    tf = CACHE / "tornado.json"
    if not tf.exists():
        _dump(tf, _tornado(acc[: min(N_TORNADO, n_accept)], geo))
    wall, cpu = time.time() - t0, time.process_time() - c0
    _dump(CACHE / "timing.json", {"wall_s": wall, "cpu_s": cpu, "n_accept": n_accept})
    ledger.append("p1b_challis/draws", "monte_carlo", n_likelihood_calls=None,
                  n_integrations=None, sim_years=None, wall_s=wall, cpu_s=cpu, n_workers=1,
                  n_accepted=n_accept)  # fmt: skip
    _log(f"draws done: wall {wall:.0f} s, cpu {cpu:.0f} s")


CONDITIONAL_ON = {
    "A0": "the printed record and the priors of § 4 (N2's fallback)",
    "C1": "the printed record and the priors of § 4 (N2's fallback)",
    "C2_1800": "N1's first fallback (the sheet's limits for the equinox of 1800)",
    "C2_date": "N1's second fallback (the sheet's limits applied to right ascension of date)",
    "C3": "N3's fallback (zones alternating outward from the start place) and N2's",
    "C4": "N3's fallback (zones alternating outward from the start place) and N2's",
}
LABELS = {
    "A0": "A0 actual procedure", "C1a": "C1a 12 Aug against 30 Jul",
    "C1b": "C1b 4 Aug against 12 Aug", "C1a_7d": "C1a, 7-day delay",
    "C1b_7d": "C1b, 7-day delay", "C2_1800": "C2 chart, equinox 1800",
    "C2_date": "C2 chart, RA of date", "C3_S4": "C3 from S4",
    "C3_S4_prompt": "C3 from S4, prompt", "C3_S5": "C3 from S5",
    "C3_S5_prompt": "C3 from S5, prompt", "C4a": "C4a Le Verrier June",
    "C4a_prompt": "C4a, prompt", "C4b": "C4b looking from 29 Jul",
}  # fmt: skip


def _stats(x: np.ndarray) -> dict:
    q = np.quantile(x, [0.05, 0.95])
    return {"mean": float(x.mean()), "q05": float(q[0]), "q95": float(q[1])}


def _load(tag: str) -> dict:
    d = CACHE / tag
    files = sorted(d.glob("chunk_[0-9][0-9][0-9].npz"))
    parts = [np.load(f) for f in files]
    out = {k: np.concatenate([p[k] for p in parts], axis=2 if k in ("t1", "t2") else
                             (1 if k == "curve" else 0))
           for k in ("t1", "t2", "curve", "params")}  # fmt: skip
    out["curve_replay_sum"] = sum(p["curve_replay_sum"] for p in parts)
    out["inner_sum"] = sum(p["inner_sum"] for p in parts)
    out["inner_maxz"] = np.max([p["inner_maxz"] for p in parts], axis=0)
    out["prior"] = np.load(d / "prior.npz")["params"]
    out["acceptance"] = json.loads((d / "acceptance.json").read_text())
    return out


def _conditional(case: str) -> str:
    for k in (case, case[:2]):
        if k in CONDITIONAL_ON:
            return CONDITIONAL_ON[k]
    raise KeyError(case)


def summarize() -> dict:
    pre = json.loads(PRERUN_JSON.read_text())
    m, h = _load("main"), _load("halved")
    tor = json.loads((CACHE / "tornado.json").read_text())
    timing = json.loads((CACHE / "timing.json").read_text())
    n = m["params"].shape[0]
    days = night_dates()
    cases = {}
    ranking = {}
    for ci, case in enumerate(C.CASES):
        rows = sorted(((abs(v["p_high"][ci] - v["p_low"][ci]), name) for name, v in
                       tor["inputs"].items()), reverse=True)  # fmt: skip
        top = rows[0][1] if rows[0][0] > 1e-9 else None
        ranking[case] = [{"input": name, "how": tor["inputs"][name]["how"],
                          "low": tor["inputs"][name]["low"], "high": tor["inputs"][name]["high"],
                          "p_low": tor["inputs"][name]["p_low"][ci],
                          "p_high": tor["inputs"][name]["p_high"][ci], "change": ch}
                         for ch, name in rows]  # fmt: skip
        entry = {"label": LABELS[case], "conditional_on": _conditional(case), "top_input": top}
        for ri, name in enumerate(READINGS):
            entry[name] = {
                "primary": _stats(m["t1"][ri, ci, :, 0]),
                "by_1_october": _stats(m["t1"][ri, ci, :, 1]),
                "t2_primary": _stats(m["t2"][ri, ci, :, 0]),
                "t2_by_1_october": _stats(m["t2"][ri, ci, :, 1]),
                "paired_minus_A0": _stats(m["t1"][ri, ci, :, 0] - m["t1"][ri, 0, :, 0]),
                "inner_minus_exact_mean": float(m["inner_sum"][ri, ci] / n),
                "inner_minus_exact_max_abs_in_se": float(m["inner_maxz"][ri, ci]),
            }
        s = entry["held"]["primary"]
        row = C.reading(s["q05"], s["q95"])
        hs = _stats(h["t1"][0, ci, :, 0])
        entry["reading_row"] = row
        entry["halved_priors"] = {"primary": hs, "reading_row": C.reading(hs["q05"], hs["q95"])}
        entry["prior_dependent"] = entry["halved_priors"]["reading_row"] != row
        tail = f"most sensitive to {top}" if top else "0 in every draw, no input moves it"
        if row == "depends on":
            text = f"depends on {top}"
        elif row == "an open chance":
            text = f"an open chance (mean {s['mean']:.2f}); {tail}"
        else:
            text = f"{row}; {tail}"
        entry["reading"] = f"{text}; conditional on {entry['conditional_on']}"
        pd_ = entry["held"]["paired_minus_A0"]
        entry["paired_reading"] = None if case == "A0" else C.paired_reading(pd_["q05"], pd_["q95"])
        cases[case] = entry
    rows = {cases[c]["reading_row"] for c in ("C2_1800", "C2_date")}
    c2 = cases["C2_1800"]["reading_row"] if len(rows) == 1 else (
        "depends on the Berlin sheet's limits")  # fmt: skip
    pairs = {}
    for i, a in enumerate(C.CASES[1:], start=1):
        for j, b in enumerate(C.CASES[1:], start=1):
            if i < j:
                pairs[f"{a} minus {b}"] = _stats(m["t1"][0, i, :, 0] - m["t1"][0, j, :, 0])
    # post hoc (identified after reading; description only): paired differences between
    # counterfactuals on the shared draws, with the share of draws on each side
    post_hoc = {}
    for a_, b_ in (("C1a", "C1b"), ("C1a", "C4b"), ("C3_S4_prompt", "C1a")):
        diff = m["t1"][0, C.CASES.index(a_), :, 0] - m["t1"][0, C.CASES.index(b_), :, 0]
        post_hoc[f"{a_} minus {b_}"] = {
            **_stats(diff),
            "share_above_0": float((diff > 1e-12).mean()),
            "share_below_0": float((diff < -1e-12).mean()),
            "share_equal": float((np.abs(diff) <= 1e-12).mean()),
        }
    mean_curve = m["curve"].mean(axis=1)
    c1_dates = {}
    for case in ("C1a", "C1b", "C1a_7d", "C1b_7d"):
        cv = mean_curve[C.CASES.index(case)]
        final = cv[C.SEP29 - 1]
        c1_dates[case] = {
            "p_by": {d: float(cv[C.night_index(d)]) for d in (
                "1846-08-13", "1846-08-14", "1846-08-16", "1846-08-20", "1846-08-31",
                "1846-09-22")},
            "median_date_given_flagged": days[int(np.argmax(cv >= 0.5 * final))] if final > 0
            else None,
            "date_95_percent_of_flags": days[int(np.argmax(cv >= 0.95 * final))] if final > 0
            else None,
        }  # fmt: skip
    post = {name: {"prior": _stats(m["prior"][:, i]), "posterior": _stats(m["params"][:, i])}
            for i, name in enumerate(PARAMS)}  # fmt: skip
    nep = pre["nights"]
    d30, d12 = nep[C.JUL30]["neptune"]["dec_deg"], nep[C.AUG12]["neptune"]["dec_deg"]
    out = {
        "plan": "docs/notes/P1B_challis_detection.md",
        "prerun_sha256": sha256(PRERUN_JSON),
        "n_accepted": n,
        "seeds": {"draws": list(SEED_DRAWS), "events": [*SEED_EVENTS, "case", "chunk"]},
        "readings": {"held": "realised events held (read)",
                     "replay": "his procedure as a plan: every event redrawn (reported, not read)"},
        "acceptance": m["acceptance"],
        "acceptance_halved_priors": h["acceptance"],
        "condition_1_unrestricted": {
            "band_arcmin": 600.0,
            "fraction_containing_planet_12_aug": C.ZONE_ARCMIN / 600.0,
            "fraction_also_excluding_30_jul_by_declination": (d30 - d12) * 60.0 / 600.0,
            "note": "a 9′ zone set uniformly in Airy's 10° of latitude",
        },
        "checks": {
            "a0_t2_by_1_october_held": cases["A0"]["held"]["t2_by_1_october"],
            "a0_t2_by_1_october_replay": cases["A0"]["replay"]["t2_by_1_october"],
            "a0_primary_held_max": float(m["t1"][0, 0, :, 0].max()),
        },
        "cases": cases,
        "c2_combined_reading_row": c2,
        "tornado": {"n_draws": tor["n_draws"], "ranking": ranking},
        "pairs_between_counterfactuals_not_read": pairs,
        "post_hoc_pairs_description_only": post_hoc,
        "c1_completion": c1_dates,
        "posterior_against_prior": post,
        "compute": {"wall_s": round(timing["wall_s"], 1), "cpu_s": round(timing["cpu_s"], 1),
                    "n_workers": 1},
    }  # fmt: skip
    _dump(RESULTS_JSON, out)
    figures(out, m, days)
    return out


INK, MUTED, BLUE, BAND = "#1f2933", "#7b8794", "#2f6fb7", "#bcd3ec"


def _style(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=INK, labelsize=8, length=3)


def figures(out: dict, m: dict, days: list[str]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(C.N_NIGHTS)
    ticks = [C.night_index(d) for d in ("1846-08-01", "1846-09-01", "1846-10-01")]
    fig, axes = plt.subplots(3, 5, figsize=(13, 7.2), sharex=True, sharey=True)
    for ax, (ci, case) in zip(axes.flat, enumerate(C.CASES), strict=False):
        cv = m["curve"][ci]
        lo, hi = np.quantile(cv, [0.05, 0.95], axis=0)
        ax.fill_between(x, lo, hi, step="post", color=BAND, linewidth=0)
        ax.step(x, cv.mean(axis=0), where="post", color=BLUE, linewidth=2)
        ax.axvline(C.PRIMARY_LAST + 1, color=MUTED, linewidth=1)
        ax.set_title(LABELS[case], fontsize=9, color=INK, loc="left")
        ax.set_ylim(0, 1)
        ax.set_xticks(ticks, ["1 Aug", "1 Sep", "1 Oct"])
        _style(ax)
    axes.flat[0].text(C.PRIMARY_LAST, 0.5, "23 Sep", fontsize=8, color=MUTED, ha="right")
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0.0, 0.6, "Line: mean over the accepted draws\nBand: their 5 to 95 % "
                       "range\nRealised events held", fontsize=9, color=INK,
                       transform=axes.flat[-1].transAxes)  # fmt: skip
    fig.supylabel("Probability the planet is flagged by the date", fontsize=10, color=INK)
    fig.suptitle("Challis's search: probability of recognition by date, per case",
                 fontsize=12, color=INK, x=0.01, ha="left")  # fmt: skip
    fig.tight_layout()
    fig.savefig(FIGURES["curves"], dpi=150)
    plt.close(fig)

    names = [p for p in PARAMS if p not in ("branch_ra", "north_first", "n_unknown_drawn")]
    fig, axes = plt.subplots(4, 4, figsize=(13, 9))
    for ax, name in zip(axes.flat, names, strict=True):
        i = PARAMS.index(name)
        a, b = m["prior"][:, i], m["params"][:, i]
        lo, hi = min(a.min(), b.min()), max(a.max(), b.max())
        bins = np.arange(lo - 0.5, hi + 1.0) if name.startswith("n_") else np.linspace(lo, hi, 31)
        ax.hist(a, bins=bins, density=True, histtype="step", color=MUTED, linewidth=1.5)
        ax.hist(b, bins=bins, density=True, color=BLUE, alpha=0.75, linewidth=0)
        ax.set_title(name.replace("_", " "), fontsize=9, color=INK, loc="left")
        ax.set_yticks([])
        _style(ax)
        ax.spines["left"].set_visible(False)
    axes.flat[0].text(0.98, 0.9, "prior (outline)\nposterior (filled)", fontsize=8, color=INK,
                      ha="right", va="top", transform=axes.flat[0].transAxes)  # fmt: skip
    fig.suptitle("Inputs before and after conditioning on the record", fontsize=12, color=INK,
                 x=0.01, ha="left")  # fmt: skip
    fig.tight_layout()
    fig.savefig(FIGURES["posterior"], dpi=150)
    plt.close(fig)

    show = [c for c in C.CASES if out["cases"][c]["top_input"]]
    fig, axes = plt.subplots(2, math.ceil(len(show) / 2), figsize=(14, 6.5), sharex=True)
    for ax, case in zip(axes.flat, show, strict=False):
        rows = [r for r in out["tornado"]["ranking"][case] if r["change"] > 0.005][:6][::-1]
        y = np.arange(len(rows))
        ax.barh(y, [r["p_high"] - r["p_low"] for r in rows], height=0.6, color=BLUE)
        ax.axvline(0, color=MUTED, linewidth=1)
        ax.set_yticks(y, [r["input"].replace("_", " ") for r in rows], fontsize=8)
        ax.set_title(LABELS[case], fontsize=9, color=INK, loc="right")
        _style(ax)
    for ax in axes.flat[len(show) :]:
        ax.axis("off")
    fig.supxlabel("Change in the probability before 23 September, input from its 5 % to its "
                  "95 % value", fontsize=10, color=INK)  # fmt: skip
    fig.suptitle("Which inputs move each case", fontsize=12, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(FIGURES["tornado"], dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prerun")
    sm = sub.add_parser("smoke")
    sm.add_argument("--n-accept", type=int, default=1000)
    sm.add_argument("--n-inner", type=int, default=N_INNER)
    dr = sub.add_parser("draws")
    dr.add_argument("--n-accept", type=int, default=N_ACCEPT)
    dr.add_argument("--n-inner", type=int, default=N_INNER)
    sub.add_parser("summarize")
    a = ap.parse_args(argv)
    if a.cmd == "prerun":
        prerun()
    elif a.cmd == "smoke":
        smoke(a.n_accept, a.n_inner)
    elif a.cmd == "draws":
        draws(a.n_accept, a.n_inner)
    elif a.cmd == "summarize":
        summarize()


if __name__ == "__main__":
    main()
