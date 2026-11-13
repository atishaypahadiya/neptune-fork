"""M6e: do the fits predict Lalande's places of 8 and 10 May 1795?

The pre-registration is ``docs/notes/M6e.md``; nothing here may differ from it. This module is
the only one that reads the held-out precovery table, and no module on the fit path imports it
(``tests/test_m6e.py``, ``tests/test_a6d.py``).

    uv run python -m fork.export.m6e fetch   # JPL's apparent place of Uranus at the instants
    uv run python -m fork.export.m6e checks  # the time and frame tests of § 2
    uv run python -m fork.export.m6e run     # docs/notes/M6e.{json,png}
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
LALANDE = ROOT / "data" / "historical" / "precoveries" / "lalande_1795.csv"
LALANDE_SHA256 = "22484a07c4e8ee568ca3b84d4f35b16a9006aeb16d935a65f5ce0b095885849a"
OUT_JSON = ROOT / "docs" / "notes" / "M6e.json"
JPL_URANUS = ROOT / "tests" / "data" / "m6e_uranus_apparent_1795.json"
NESTED = ROOT / "data" / "cache" / "nested"

# § 2 of the plan: Paris mean time from mean noon; UT = PMT − 9m 20.9s; TT = UT + 15 s; TDB = TT.
INSTANTS_TDB = {"1795-05-08": "1795-05-08T22:56:05", "1795-05-10": "1795-05-10T22:48:03"}
PARIS_EAST_LONGITUDE_DEG = 2.0 + 20.0 / 60.0 + 14.0 / 3600.0
DELTA_T_S = 15.0
TOL_TRANSIT_MIN = 1.0
TOL_FRAME_ARCSEC = 1.0
TOL_GALLE_MEDIAN_DEG = 0.2

N_DRAWS = 6000
SEED = 20261113
MODELS = {"headline": 0, "student_t": 1}  # the seed's second entry; the prior's is 10 + this
PUBLISHED_GALLE_MEDIAN_DEG = {"headline": 325.815, "student_t": 325.166}
POSTERIOR_SHA256 = {
    "m6c/headline__c0.pkl": "529407e435cadf99aa86889934e4fa5b525f2bc1b502973a8150aa720a665787",
    "m6c/headline__c1.pkl": "1cc6ff7b42bddcd4da80437aeb24d26ba511ec53dc975eabfe897c7a37812788",
    "m6c/headline__c3.pkl": "3fe3ee5434808d276d4db5aded4c4afc9ee6189bfb42af62cf223953d4f7e6b7",
    "m6c_repeat/r2_1.pkl": "72401426b028bbe77c4504d2f815c61705612d5a50896111eff5ac3508d34fe0",
    "m6c_repeat/r2_2.pkl": "33df8e2d5bafd97c3b68ee4c39fdb9b1e8470400a8f59d01454d62ee0e871b4f",
    "m6c/s3_t4__c0.pkl": "6b7f2c50ca356b32b615240f160d2eab91edd85baa482c10bd96bd5b6774d4d0",
    "m6c/s3_t4__c1.pkl": "8f6203df183399be1ac9f39f43cb6630b5f7e7534671da9417b9a48225100b4f",
    "m6c/s3_t4__c2.pkl": "d813abb56ff29af396248f099ccdcffea497e879d346684d929c41275e5c60f0",
    "m6c/s3_t4__c3.pkl": "ef7e0e6efa479235baf047de523c3f9a2910ba12816e1d3380f7a9a12b7c8ecb",
}
LEVELS = (0.68, 0.95)
BOUNDARY = 0.01
ARCSEC = 3600.0
OBLIQUITY_J2000_RAD = math.radians(84381.448 / ARCSEC)  # IAU 1976, the frame of the states


# ------------------------------------------------------------------ the observed places


def lalande() -> dict[str, dict]:
    """The two places: apparent right ascension and declination in degrees."""
    if hashlib.sha256(LALANDE.read_bytes()).hexdigest() != LALANDE_SHA256:
        raise RuntimeError("the Lalande table is not the one of the pre-registration")
    out = {}
    with LALANDE.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            ra = (
                float(r["ra_app_deg"]) + float(r["ra_app_min"]) / 60 + float(r["ra_app_sec"]) / 3600
            )
            dec = (
                float(r["dec_app_deg"])
                + float(r["dec_app_min"]) / 60
                + float(r["dec_app_sec"]) / 3600
            )
            out[r["date"]] = {
                "paris_mean_time": r["paris_mean_time"],
                "ra_deg": ra,
                "dec_deg": -dec if r["dec_app_sign"] == "-" else dec,
            }
    return out


# ------------------------------------------------------------------ time and frame (§ 2)


def jd_of(date: str) -> float:
    from fork.data.frames import jd_tdb_from_iso

    return float(jd_tdb_from_iso(INSTANTS_TDB[date]))


def instant_from_paris_mean_time(date: str, hms: str) -> str:
    """The TDB instant of a Paris mean time reckoned from mean noon (ISO, whole seconds)."""
    from astropy.time import Time, TimeDelta

    h, m, s = (int(v) for v in hms.split(":"))
    civil = Time(f"{date}T12:00:00", scale="tt") + TimeDelta(h * 3600 + m * 60 + s, format="sec")
    tt = civil - TimeDelta(PARIS_EAST_LONGITUDE_DEG * 240.0 - DELTA_T_S, format="sec")
    return (tt + TimeDelta(0.5, format="sec")).isot[:19]


def paris_sidereal_time_deg(jd_tdb: float) -> float:
    """Apparent sidereal time at the Paris Observatory, degrees."""
    import erfa

    ut = jd_tdb - DELTA_T_S / 86400.0
    gst = erfa.gst06a(2400000.5, ut - 2400000.5, 2400000.5, jd_tdb - 2400000.5)
    return (math.degrees(gst) + PARIS_EAST_LONGITUDE_DEG) % 360.0


class Epoch:
    """What the apparent place needs at one instant: Earth's place and velocity (about the
    inner barycenter, ecl_j2000) and the rotation to the true equator and equinox of date."""

    def __init__(self, date: str):
        import erfa

        from fork.data.horizons import states_about_inner_barycenter_ecl_j2000

        self.date = date
        self.jd = jd_of(date)
        st = states_about_inner_barycenter_ecl_j2000("earth", [self.jd])[0]
        self.earth_pos, self.earth_vel = st[:3], st[3:]
        d1, d2 = 2400000.5, self.jd - 2400000.5
        self.npb = np.asarray(erfa.pnm06a(d1, d2))  # GCRS → true equator and equinox of date
        _, deps = erfa.nut06a(d1, d2)
        self.true_obliquity = float(erfa.obl06(d1, d2) + deps)

    def apparent(self, pos: np.ndarray, vel: np.ndarray) -> tuple[float, float, float]:
        """Apparent right ascension, declination and ecliptic longitude of date (degrees) of
        a body with state ``pos``, ``vel`` (ecl_j2000, the origin of ``earth_pos``)."""
        import erfa

        from fork.data.frames import C_AU_PER_DAY, geocentric_lon_lat_ecl_j2000_light_time

        lon, lat, _ = geocentric_lon_lat_ecl_j2000_light_time(pos, vel, self.earth_pos)
        u = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
        u, v = _to_equator(u), _to_equator(self.earth_vel) / C_AU_PER_DAY
        sun = float(np.linalg.norm(self.earth_pos))
        u = erfa.ab(u, v, sun, math.sqrt(1.0 - float(v @ v)))  # annual aberration
        w = self.npb @ u
        ra = math.degrees(math.atan2(w[1], w[0])) % 360.0
        dec = math.degrees(math.asin(w[2] / float(np.linalg.norm(w))))
        return ra, dec, ecliptic_longitude_deg(ra, dec, self.true_obliquity)


def _to_equator(v: np.ndarray) -> np.ndarray:
    c, s = math.cos(OBLIQUITY_J2000_RAD), math.sin(OBLIQUITY_J2000_RAD)
    return np.array([v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]])


def ecliptic_longitude_deg(ra_deg: float, dec_deg: float, obliquity_rad: float) -> float:
    """Ecliptic longitude of a place given on the equator of the same date."""
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    y = math.sin(a) * math.cos(d) * math.cos(obliquity_rad) + math.sin(d) * math.sin(obliquity_rad)
    return math.degrees(math.atan2(y, math.cos(a) * math.cos(d))) % 360.0


def fetch_jpl_uranus() -> dict:
    """JPL Horizons' geocentric apparent place of Uranus at the two instants (one request;
    Uranus only), written to ``tests/data`` for the frame test."""
    from astroquery.jplhorizons import Horizons

    out = {
        "source": "JPL Horizons, Uranus system barycenter (7, the body of the cached "
        "vectors), geocentric, apparent RA and DEC (airless), extra precision",
        "delta_t_s": DELTA_T_S,
        "places": {},
    }
    for date in INSTANTS_TDB:
        jd_ut = jd_of(date) - DELTA_T_S / 86400.0
        eph = Horizons(id="7", location="500@399", epochs=[jd_ut]).ephemerides(
            quantities="2", extra_precision=True
        )
        out["places"][date] = {
            "jd_ut": jd_ut,
            "ra_app_deg": float(eph["RA_app"][0]),
            "dec_app_deg": float(eph["DEC_app"][0]),
        }
    JPL_URANUS.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    return out


def checks() -> dict:
    """The two tests of § 2 that must pass before the run: the transit times and the frame."""
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000

    obs, out = lalande(), {"transit": {}, "frame": {}}
    jpl = json.loads(JPL_URANUS.read_text(encoding="utf-8"))["places"]
    for date in INSTANTS_TDB:
        iso = instant_from_paris_mean_time(date, obs[date]["paris_mean_time"])
        lst = paris_sidereal_time_deg(jd_of(date))
        ha_min = ((lst - obs[date]["ra_deg"] + 180.0) % 360.0 - 180.0) * 4.0
        out["transit"][date] = {
            "instant_from_the_table_tdb": iso,
            "instant_of_the_plan_tdb": INSTANTS_TDB[date],
            "hour_angle_min": ha_min,
            "pass": iso == INSTANTS_TDB[date] and abs(ha_min) <= TOL_TRANSIT_MIN,
        }
        ep = Epoch(date)
        st = states_about_inner_barycenter_ecl_j2000("uranus", [ep.jd])[0]
        ra, dec, _ = ep.apparent(st[:3], st[3:])
        dra = ((ra - jpl[date]["ra_app_deg"] + 180.0) % 360.0 - 180.0) * ARCSEC
        dra *= math.cos(math.radians(dec))
        ddec = (dec - jpl[date]["dec_app_deg"]) * ARCSEC
        out["frame"][date] = {
            "ra_minus_jpl_arcsec_on_sky": dra,
            "dec_minus_jpl_arcsec": ddec,
            "pass": abs(dra) <= TOL_FRAME_ARCSEC and abs(ddec) <= TOL_FRAME_ARCSEC,
        }
    out["pass"] = all(v["pass"] for part in ("transit", "frame") for v in out[part].values())
    return out


# ------------------------------------------------------------------ levels (§§ 3, 4)


def central_level(draws_deg: np.ndarray, observed_deg: float) -> float:
    """L = |2F − 1|, F the fraction of draws below the observed longitude, differences taken
    on the circle about the draws' circular median."""
    x = np.asarray(draws_deg, dtype=float)
    centre = math.degrees(math.atan2(np.sin(np.radians(x)).mean(), np.cos(np.radians(x)).mean()))
    d = (x - centre + 180.0) % 360.0 - 180.0
    centre += float(np.median(d))
    d = (x - centre + 180.0) % 360.0 - 180.0
    o = (observed_deg - centre + 180.0) % 360.0 - 180.0
    return abs(2.0 * float(np.mean(d < o)) - 1.0)


def level_class(level: float) -> str:
    name = (
        "inside 68 %" if level <= LEVELS[0] else "inside 95 %" if level <= LEVELS[1] else "outside"
    )
    for edge, pair in zip(
        LEVELS, ("inside 68 % or inside 95 %", "inside 95 % or outside"), strict=True
    ):
        if abs(level - edge) <= BOUNDARY:
            return f"at the boundary ({pair})"
    return name


def interval(x_deg: np.ndarray, centre_deg: float | None = None) -> dict:
    """Median and central 68 % and 95 % intervals of an angle, unwrapped about its centre."""
    x = np.asarray(x_deg, dtype=float)
    if centre_deg is None:
        centre_deg = math.degrees(
            math.atan2(np.sin(np.radians(x)).mean(), np.cos(np.radians(x)).mean())
        )
    d = (x - centre_deg + 180.0) % 360.0 - 180.0
    q = np.percentile(d, [2.5, 16.0, 50.0, 84.0, 97.5]) + centre_deg
    if centre_deg != 0.0:  # an angle on the circle: the median is given in [0°, 360°)
        q -= 360.0 * math.floor(q[2] / 360.0)
    return {
        "median": float(q[2]),
        "lo68": float(q[1]),
        "hi68": float(q[3]),
        "lo95": float(q[0]),
        "hi95": float(q[4]),
        "width68": float(q[3] - q[1]),
        "width95": float(q[4] - q[0]),
    }


def ecliptic_latitude_deg(ra_deg, dec_deg, obliquity_rad: float):
    """Ecliptic latitude of places given on the equator of the same date."""
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    return np.degrees(
        np.arcsin(
            np.sin(d) * math.cos(obliquity_rad) - np.cos(d) * math.sin(obliquity_rad) * np.sin(a)
        )
    )


def latitude_summary(p: dict, o: dict) -> dict:
    """Post hoc (docs/notes/M6e.md): the draws' ecliptic latitude of date and Lalande's."""
    eps = p["true_obliquity_rad"]
    q = np.percentile(
        ecliptic_latitude_deg(p["ra_deg"], p["dec_deg"], eps), [2.5, 16, 50, 84, 97.5]
    )
    return {
        "observed": float(ecliptic_latitude_deg(o["ra_deg"], o["dec_deg"], eps)),
        "median": float(q[2]),
        "lo68": float(q[1]),
        "hi68": float(q[3]),
        "lo95": float(q[0]),
        "hi95": float(q[4]),
    }


def summarize_draws(places: dict, obs: dict) -> dict:
    """Per night: the levels of Lalande's place and the intervals of the draws."""
    from fork.export import v13

    out = {}
    for date in INSTANTS_TDB:
        p, o = places[date], obs[date]
        o_lon = ecliptic_longitude_deg(o["ra_deg"], o["dec_deg"], p["true_obliquity_rad"])
        level = central_level(p["lon_deg"], o_lon)
        region = v13.region_summary(p["ra_deg"] / 15.0, p["dec_deg"])
        lon, ra, dec = interval(p["lon_deg"]), interval(p["ra_deg"]), interval(p["dec_deg"], 0.0)
        out[date] = {
            "observed": {"ra_deg": o["ra_deg"], "dec_deg": o["dec_deg"], "ecl_lon_deg": o_lon},
            "central_level_ecl_lon": level,
            "class": level_class(level),
            "tangent_plane_level": v13.region_level(
                region, p["ra_deg"] / 15.0, p["dec_deg"], o["ra_deg"], o["dec_deg"]
            ),
            "ecl_lon_deg": lon,
            "ra_deg": ra,
            "dec_deg": dec,
            "median_minus_observed_deg": {
                "ecl_lon": (lon["median"] - o_lon + 180.0) % 360.0 - 180.0,
                "ra": (ra["median"] - o["ra_deg"] + 180.0) % 360.0 - 180.0,
                "dec": dec["median"] - o["dec_deg"],
            },
        }
    a, b = (places[d] for d in INSTANTS_TDB)
    dra = ((b["ra_deg"] - a["ra_deg"] + 180.0) % 360.0 - 180.0) * ARCSEC
    ddec = (b["dec_deg"] - a["dec_deg"]) * ARCSEC
    oa, ob = (obs[d] for d in INSTANTS_TDB)
    out["motion_between_the_nights_arcsec"] = {
        "ra": {
            "median": float(np.median(dra)),
            "lo68": float(np.percentile(dra, 16)),
            "hi68": float(np.percentile(dra, 84)),
            "observed": (ob["ra_deg"] - oa["ra_deg"]) * ARCSEC,
        },
        "dec": {
            "median": float(np.median(ddec)),
            "lo68": float(np.percentile(ddec, 16)),
            "hi68": float(np.percentile(ddec, 84)),
            "observed": (ob["dec_deg"] - oa["dec_deg"]) * ARCSEC,
        },
    }
    return out


# ------------------------------------------------------------------ draws


def predictive(space, ctx, q: np.ndarray, epochs: dict[str, Epoch]) -> dict:
    """The perturber's apparent place at each instant for each row of ``q`` (one integration
    per row). Rows whose integration fails are dropped and counted."""
    from fork.model.forward import predict_states_ecl_j2000

    dates = list(epochs)
    jds = [epochs[d].jd for d in dates]
    vals = {d: [] for d in dates}
    failed = 0
    for v in q:
        try:
            pos, vel = predict_states_ecl_j2000(space.theta(v), jds, ctx)
            row = [epochs[d].apparent(pos[j, 3], vel[j, 3]) for j, d in enumerate(dates)]
        except (ValueError, FloatingPointError, ZeroDivisionError):
            failed += 1
            continue
        if not np.all(np.isfinite(row)):
            failed += 1
            continue
        for d, r in zip(dates, row, strict=True):
            vals[d].append(r)
    out = {}
    for d in dates:
        arr = np.array(vals[d])
        out[d] = {
            "ra_deg": arr[:, 0],
            "dec_deg": arr[:, 1],
            "lon_deg": arr[:, 2],
            "true_obliquity_rad": epochs[d].true_obliquity,
        }
    out["n_failed"] = failed
    return out


def cells_of(model: str) -> list:
    """The published posterior's per-cell nested results (§ 1 of the plan)."""
    from fork.export import m6c, m6c_repeat
    from fork.inference.nested import NestedResult

    if model == "headline":
        return m6c_repeat.replacement_cells()
    return [NestedResult.load(m6c.Job("s3_t4", k).path(".pkl")) for k in range(4)]


def posterior_draws(model: str, cells: list) -> np.ndarray:
    from fork.inference.nested import ERR_INFLATION, PartitionedResult

    rng = np.random.default_rng([SEED, MODELS[model]])
    eq = PartitionedResult(tuple(cells), err_scale=ERR_INFLATION).equal_weight(rng)
    return eq[rng.choice(len(eq), size=N_DRAWS, replace=False)]


def prior_draws(model: str, cells: list, n: int = N_DRAWS) -> np.ndarray:
    """Draws from the fit's prior: a cell with its prior mass, then its prior transform."""
    rng = np.random.default_rng([SEED, 10 + MODELS[model]])
    w = np.array([c.space.slice_weight for c in cells])
    which = rng.choice(len(cells), size=n, p=w / w.sum())
    return np.array(
        [cells[k].space.prior_transform(rng.random(cells[k].space.ndim)) for k in which]
    )


def check_hashes() -> None:
    from fork.export import m6c

    m6c.check_frozen("headline")
    m6c.check_frozen("s3_t4")
    for rel, digest in POSTERIOR_SHA256.items():
        if hashlib.sha256((NESTED / rel).read_bytes()).hexdigest() != digest:
            raise RuntimeError(f"{rel} is not the file of the pre-registration")


def run() -> dict:
    from fork import ledger
    from fork.export import m6c
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad

    t0, c0 = time.time(), ledger.cpu_now()
    check_hashes()
    pre = checks()
    if not pre["pass"]:
        raise SystemExit(f"a test of § 2 fails; stopped before the run: {json.dumps(pre)}")
    obs = lalande()
    epochs = {d: Epoch(d) for d in INSTANTS_TDB}
    geo = GalleGeometry.from_horizons()
    out: dict = {
        "plan": "docs/notes/M6e.md",
        "instants_tdb": INSTANTS_TDB,
        "checks": pre,
        "n_draws": N_DRAWS,
        "models": {},
    }
    for model in MODELS:
        cells = cells_of(model)
        space = cells[0].space
        _, ctx, *_ = m6c.dataset(variant="headline" if model == "headline" else "s3_t4")
        q = posterior_draws(model, cells)
        galle = np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo)) for v in q])
        med = interval(galle)["median"] % 360.0
        galle_width68 = interval(galle)["width68"]
        if abs(med - PUBLISHED_GALLE_MEDIAN_DEG[model]) > TOL_GALLE_MEDIAN_DEG:
            raise SystemExit(f"{model}: Galle-night median {med:.3f} is not the published one")
        post = predictive(space, ctx, q, epochs)
        prior = predictive(space, ctx, prior_draws(model, cells), epochs)
        entry = {
            "galle_night_median_deg": med,
            "published_galle_night_median_deg": PUBLISHED_GALLE_MEDIAN_DEG[model],
            "posterior": summarize_draws(post, obs) | {"n_failed": post["n_failed"]},
            "prior_predictive": summarize_draws(prior, obs) | {"n_failed": prior["n_failed"]},
        }
        entry["width_ratio_posterior_to_prior"] = {
            d: {
                k: entry["posterior"][d]["ecl_lon_deg"][k]
                / entry["prior_predictive"][d]["ecl_lon_deg"][k]
                for k in ("width68", "width95")
            }
            for d in INSTANTS_TDB
        }
        entry["histogram_8_may"] = _histogram(post["1795-05-08"]["lon_deg"])
        entry["post_hoc"] = {
            "ecl_lat_deg": {d: latitude_summary(post[d], obs[d]) for d in INSTANTS_TDB},
            "galle_night_ecl_lon_width68_deg": galle_width68,
        }
        classes = [entry["posterior"][d]["class"] for d in INSTANTS_TDB]
        entry["reading"] = (
            classes[0]
            if classes[0] == classes[1]
            else (f"at a boundary between nights: {classes[0]} (8 May), {classes[1]} (10 May)")
        )
        out["models"][model] = entry
        print(
            model,
            entry["reading"],
            {d: round(entry["posterior"][d]["central_level_ecl_lon"], 4) for d in INSTANTS_TDB},
            flush=True,
        )
    n = 4 * N_DRAWS
    out["cost"] = {"wall_s": time.time() - t0, "cpu_s": ledger.cpu_now() - c0, "n_integrations": n}
    ledger.append(
        "m6e/predictive",
        "summarize",
        n_likelihood_calls=0,
        n_integrations=n + 2 * N_DRAWS,
        sim_years=n * (1846.7 - 1795.35),
        wall_s=out["cost"]["wall_s"],
        cpu_s=out["cost"]["cpu_s"],
        n_workers=1,
        estimated=True,
        note="posterior and prior draws carried to 8 and 10 May 1795",
    )
    OUT_JSON.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    figure(out)
    return out


def _histogram(lon_deg: np.ndarray, bin_deg: float = 0.5) -> dict:
    centre = interval(lon_deg)["median"]
    d = (lon_deg - centre + 180.0) % 360.0 - 180.0
    lo, hi = math.floor(d.min() / bin_deg) * bin_deg, math.ceil(d.max() / bin_deg) * bin_deg
    counts, edges = np.histogram(d, bins=np.arange(lo, hi + bin_deg / 2, bin_deg))
    return {"edges_deg": (edges + centre).tolist(), "counts": counts.tolist()}


def figure(out: dict) -> None:
    """The posterior draws' ecliptic longitude on 8 May 1795 per model, with Lalande's place."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = {"headline": "Gaussian errors (headline)", "student_t": "Student-t errors"}
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)
    for ax, (model, e) in zip(axes, out["models"].items(), strict=True):
        h, p = e["histogram_8_may"], e["posterior"]["1795-05-08"]
        edges = np.array(h["edges_deg"])
        ax.stairs(h["counts"], edges, fill=True, color="#9db7d5")
        lon = p["ecl_lon_deg"]
        ax.axvspan(lon["lo68"], lon["hi68"], ymax=0.04, color="#1f4e79")
        ax.axvline(p["observed"]["ecl_lon_deg"], color="#b3261e", lw=1.4)
        ax.set_title(f"{names[model]}: {p['class']}", loc="left", fontsize=10)
        ax.set_ylabel("draws")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0].annotate(
        "Lalande, 8 May 1795",
        (out["models"]["headline"]["posterior"]["1795-05-08"]["observed"]["ecl_lon_deg"], 0.9),
        xycoords=("data", "axes fraction"),
        xytext=(5, 0),
        textcoords="offset points",
        color="#b3261e",
        fontsize=9,
    )
    axes[1].set_xlabel("predicted ecliptic longitude of date (degrees); bar: central 68 %")
    fig.tight_layout()
    fig.savefig(OUT_JSON.with_suffix(".png"), dpi=160)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    cmd = (argv or sys.argv[1:] or ["run"])[0]
    if cmd == "fetch":
        print(json.dumps(fetch_jpl_uranus(), indent=1))
    elif cmd == "checks":
        print(json.dumps(checks(), indent=1))
    elif cmd == "run":
        run()
    else:
        raise SystemExit(f"unknown command {cmd}")


if __name__ == "__main__":
    main()
