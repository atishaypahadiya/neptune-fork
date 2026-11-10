"""Uranus and Neptune, 1781–1850, integrated with and without Neptune's mass (site figure data).

    uv run python docs/notes/scripts/m7_motion.py

What it does. Two REBOUND IAS15 runs of the Sun and the eight planets from Herschel's discovery
epoch (1781-03-13.0 TDB, ``fork.export.to_json.A1_EPOCH_JD_TDB``) to J1850.0. Run 1 uses the real
masses; run 2 is identical except that Neptune's mass is zero (Neptune is still integrated, as a
test particle). Both start from the same heliocentric states, so run 2 is "the same planets on the
same day, without Neptune's pull". Outputs every 0.1 Julian year on a fixed grid from the start
epoch, plus a 1-day grid over 1810–1835 for the conjunction and every cached Horizons epoch of
Uranus and Neptune for validation (the integrator stops exactly at each requested time).

Inputs (read-only, no network):
- initial states: ``fork.data.horizons.state_vectors_ecl_j2000`` from the local cache in
  data/cache/; every call passes a fetcher that raises, so a cache miss fails instead of fetching;
- masses and G: ``fork.model.bodies`` (IAU 2009 Sun/planet mass ratios; G = k², Gaussian k). The
  Sun is 1 M_sun with the four inner planets as explicit bodies; ``M_CENTRAL_MSUN`` (which folds
  the inner planets into the Sun for the project's own model) is not used;
- JPL's conjunction epoch: site/data/m6c_conjunction.json.

Conjunction: as in docs/notes/scripts/m6c_conjunction.py (``lon_of_date_deg``,
``conjunction_jd``, ``jd_to_date``, copied below; ``lon_of_date_deg`` is vectorised over epochs,
same rotation matrices), on its 10-day 1810–1835 grid and on a 1-day grid.

Output: site/data/m7_motion.json, per-frame data as column arrays.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import rebound  # noqa: TID251 (docs/RULES.md rule 1 bans it from src/; this is a notes script)

from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
from fork.data.horizons import _cache_path, _load_cache, cache_dir, state_vectors_ecl_j2000
from fork.export.to_json import A1_EPOCH_JD_TDB
from fork.model import bodies

ROOT = Path(__file__).resolve().parents[3]
OUT_JSON = ROOT / "site" / "data" / "m7_motion.json"
M6C_JSON = ROOT / "site" / "data" / "m6c_conjunction.json"
SCRIPT = "docs/notes/scripts/m7_motion.py"

JD_J2000 = 2451545.0
JULIAN_YEAR_DAYS = 365.25
JD_END = JD_J2000 + (1850.0 - 2000.0) * JULIAN_YEAR_DAYS  # J1850.0 TDB
STEP_YEARS = 0.1
JD_1810, JD_1835 = 2382148.5, 2391279.5  # 1810 Jan 1.0 and 1835 Jan 1.0 TDB (as in m6c)
COARSE_STEP_DAYS, FINE_STEP_DAYS = 10.0, 1.0
ARCSEC_PER_DEG = 3600.0
ARCSEC_PER_RAD = math.degrees(1.0) * ARCSEC_PER_DEG

PLANETS = ("mercury", "venus", "earth_moon", "mars", "jupiter", "saturn", "uranus", "neptune")
I_URANUS, I_NEPTUNE = 1 + PLANETS.index("uranus"), 1 + PLANETS.index("neptune")  # 0 is the Sun
MASS_RATIO = {
    "mercury": bodies.MASS_RATIO_MERCURY,
    "venus": bodies.MASS_RATIO_VENUS,
    "earth_moon": bodies.MASS_RATIO_EARTH_MOON,
    "mars": bodies.MASS_RATIO_MARS,
    "jupiter": bodies.MASS_RATIO_JUPITER,
    "saturn": bodies.MASS_RATIO_SATURN,
    "uranus": bodies.MASS_RATIO_URANUS,
    "neptune": bodies.MASS_RATIO_NEPTUNE,
}
MASS_KIND = {
    "mercury": "planet",
    "venus": "planet",
    "earth_moon": "Earth + Moon (Earth-Moon barycentre)",
    "mars": "system (planet + satellites)",
    "jupiter": "system (planet + satellites)",
    "saturn": "system (planet + satellites)",
    "uranus": "system (planet + satellites)",
    "neptune": "system (planet + satellites)",
}
MASS_SOURCE = "IAU 2009 system of astronomical constants, Sun/planet mass ratio (fork.model.bodies)"


def no_fetch(horizons_id: str, jd_tdb: np.ndarray) -> np.ndarray:
    """Fetcher that refuses: a Horizons cache miss is an error, never a network call."""
    raise RuntimeError(f"Horizons cache miss for id={horizons_id} ({jd_tdb.size} epochs)")


def horizons(body: str, jd_tdb: np.ndarray) -> np.ndarray:
    """Sun-centred ecl_j2000 states (n, 6) [AU, AU/day] from the local cache only."""
    return state_vectors_ecl_j2000(body, jd_tdb, fetcher=no_fetch)


# ---- copied from docs/notes/scripts/m6c_conjunction.py (lon_of_date_deg vectorised) ----------


def wrap180(x: np.ndarray | float) -> np.ndarray:
    """Wrap degrees to [-180, 180)."""
    return (np.asarray(x) + 180.0) % 360.0 - 180.0


def lon_of_date_deg(pos_ecl_j2000: np.ndarray, jd: np.ndarray) -> np.ndarray:
    """Ecliptic-of-date longitude, degrees, of positions (n, 3) at epochs jd (n,)."""
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(np.asarray(jd, dtype=np.float64))
    v = np.einsum("nij,nj->ni", rot.reshape(-1, 3, 3), pos_ecl_j2000)
    return np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360.0


def conjunction_jd(jd: np.ndarray, lon_uranus: np.ndarray, lon_planet: np.ndarray) -> float:
    """First epoch where Uranus overtakes the planet in heliocentric longitude (NaN if none)."""
    d = wrap180(lon_uranus - lon_planet)
    k = np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0) & (np.abs(d[:-1]) < 90))
    if len(k) == 0:
        return float("nan")
    i = k[0]
    return float(jd[i] + (jd[i + 1] - jd[i]) * (-d[i]) / (d[i + 1] - d[i]))


def jd_to_date(jd: float) -> str:
    """'YYYY-MM-DD HH:MM' (TDB) for a TDB Julian date."""
    from astropy.time import Time

    return Time(jd, format="jd", scale="tdb").iso[:16]


# ---------------------------------------------------------------------------------------------


def jd_to_year(jd: np.ndarray | float) -> np.ndarray:
    """Decimal (Julian) year, 2000 + (jd - 2451545) / 365.25."""
    return 2000.0 + (np.asarray(jd) - JD_J2000) / JULIAN_YEAR_DAYS


def sig3(x: float) -> float:
    """Round to 3 significant figures."""
    return float(f"{x:.3g}")


def build_sim(neptune_mass: bool) -> rebound.Simulation:
    """Sun + 8 planets at the A1 epoch (time 0 = that epoch, days), centre-of-mass frame."""
    sim = rebound.Simulation()
    sim.G = bodies.G  # AU^3 / (M_sun day^2)
    sim.integrator = "ias15"
    sim.add(m=bodies.M_SUN_MSUN, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
    for name in PLANETS:
        st = horizons(name, np.array([A1_EPOCH_JD_TDB]))[0]
        m = 0.0 if (name == "neptune" and not neptune_mass) else 1.0 / MASS_RATIO[name]
        sim.add(m=m, x=st[0], y=st[1], z=st[2], vx=st[3], vy=st[4], vz=st[5])
    sim.move_to_com()
    return sim


def run(neptune_mass: bool, jd_out: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """Integrate to each epoch in jd_out (sorted); heliocentric positions and velocities
    (n, 9, 3) each, plus integrator diagnostics."""
    sim = build_sim(neptune_mass)
    e0 = sim.energy()
    n = len(jd_out)
    pos, vel = np.empty((n, 1 + len(PLANETS), 3)), np.empty((n, 1 + len(PLANETS), 3))
    buf = np.empty((1 + len(PLANETS), 6))
    for k, jd in enumerate(jd_out):
        sim.integrate(jd - A1_EPOCH_JD_TDB, exact_finish_time=1)
        for i, p in enumerate(sim.particles):
            buf[i] = (p.x, p.y, p.z, p.vx, p.vy, p.vz)
        pos[k], vel[k] = buf[:, :3] - buf[0, :3], buf[:, 3:] - buf[0, 3:]
    diag = {
        "relative_energy_error_end": sig3(abs((sim.energy() - e0) / e0)),
        "steps_done": int(sim.steps_done),
        "iterations_max_exceeded": int(sim.integrator.iterations_max_exceeded),
        "epsilon": float(sim.integrator.epsilon),
        "adaptive_mode": str(sim.integrator.adaptive_mode),
        "min_dt": float(sim.integrator.min_dt),
    }
    return pos, vel, diag


def ecl_lon_deg(xyz: np.ndarray) -> np.ndarray:
    """Ecliptic-J2000 longitude, degrees, of positions (..., 3)."""
    return np.degrees(np.arctan2(xyz[..., 1], xyz[..., 0]))


def angle_arcsec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Angle between position vectors (n, 3) and (n, 3), arcsec."""
    c = np.linalg.norm(np.cross(a, b), axis=1)
    return np.arctan2(c, np.einsum("ni,ni->n", a, b)) * ARCSEC_PER_RAD


def neptune_accel(r_u: np.ndarray, r_n: np.ndarray, gm_n: float) -> np.ndarray:
    """Heliocentric perturbing acceleration on Uranus by Neptune (direct + indirect), (n, 3)."""
    d = r_n - r_u
    direct = d / np.linalg.norm(d, axis=1, keepdims=True) ** 3
    indirect = r_n / np.linalg.norm(r_n, axis=1, keepdims=True) ** 3
    return gm_n * (direct - indirect)


def first_sign_change(jd: np.ndarray, f: np.ndarray, near_jd: float) -> float:
    """Epoch (linear interpolation) of the + to - crossing of f nearest near_jd (NaN if none)."""
    k = np.flatnonzero((f[:-1] > 0) & (f[1:] <= 0))
    if len(k) == 0:
        return float("nan")
    x = jd[k] + (jd[k + 1] - jd[k]) * f[k] / (f[k] - f[k + 1])
    return float(x[np.argmin(np.abs(x - near_jd))])


def exaggeration(max_offset_au: float) -> int | float:
    """Smallest K in {1, 2, 5} x 10^k with 0.8 <= K * max_offset <= 2 AU."""
    for k in range(-3, 9):
        for m in (1, 2, 5):
            K = m * 10.0**k
            if 0.8 <= K * max_offset_au <= 2.0:
                return int(K) if K >= 1 else K
    raise ValueError(f"no exaggeration factor for max offset {max_offset_au} AU")


def main() -> None:
    t0 = time.perf_counter()

    # -- output epochs --------------------------------------------------------------------------
    n_frames = math.floor((JD_END - A1_EPOCH_JD_TDB) / (STEP_YEARS * JULIAN_YEAR_DAYS)) + 1
    jd_frames = A1_EPOCH_JD_TDB + np.arange(n_frames) * STEP_YEARS * JULIAN_YEAR_DAYS
    jd_coarse = np.arange(JD_1810, JD_1835 + 1, COARSE_STEP_DAYS)  # as in m6c_conjunction.py
    jd_fine = np.arange(JD_1810, JD_1835 + 1, FINE_STEP_DAYS)
    jd_val = {}
    for body in ("uranus", "neptune"):
        cached, _ = _load_cache(_cache_path(body, cache_dir()))
        jd_val[body] = cached[(cached >= A1_EPOCH_JD_TDB) & (cached <= JD_END)]
    jd_all = np.unique(np.concatenate([jd_frames, jd_coarse, jd_fine, *jd_val.values()]))
    at = {
        name: np.searchsorted(jd_all, v)
        for name, v in (
            ("frames", jd_frames),
            ("coarse", jd_coarse),
            ("fine", jd_fine),
            ("uranus", jd_val["uranus"]),
            ("neptune", jd_val["neptune"]),
        )
    }

    # -- the two runs ---------------------------------------------------------------------------
    pos1, vel1, diag1 = run(True, jd_all)
    pos2, _, diag2 = run(False, jd_all)
    print(f"integrated {len(jd_all)} epochs twice in {time.perf_counter() - t0:.1f} s")

    # -- per-frame quantities -------------------------------------------------------------------
    f = at["frames"]
    r_u, v_u, r_n = pos1[f, I_URANUS], vel1[f, I_URANUS], pos1[f, I_NEPTUNE]
    off = pos2[f, I_URANUS] - r_u
    dlon = wrap180(ecl_lon_deg(r_u) - ecl_lon_deg(pos2[f, I_URANUS])) * ARCSEC_PER_DEG
    gm_n = bodies.G / bodies.MASS_RATIO_NEPTUNE
    acc = neptune_accel(r_u, r_n, gm_n)
    along = np.einsum("ni,ni->n", acc, v_u / np.linalg.norm(v_u, axis=1, keepdims=True))
    radial = np.einsum("ni,ni->n", acc, r_u / np.linalg.norm(r_u, axis=1, keepdims=True))
    mag = np.linalg.norm(acc, axis=1)
    off_xy = np.hypot(off[:, 0], off[:, 1])
    K = exaggeration(float(off_xy.max()))
    i_dlon = int(np.argmax(np.abs(dlon)))

    # -- conjunction (run 1, heliocentric, ecliptic of date) ------------------------------------
    jpl = json.loads(M6C_JSON.read_text())["epoch"]
    jd_jpl = float(jpl["jpl_uranus_neptune_conjunction_jd_tdb"])
    # JPL on the same 10-day grid, heliocentric (unrounded; m6c's value is about the inner
    # barycentre and rounded to 0.1 d)
    jd_c = jd_all[at["coarse"]]
    c_jpl_helio = conjunction_jd(
        jd_c,
        lon_of_date_deg(horizons("uranus", jd_c)[:, :3], jd_c),
        lon_of_date_deg(horizons("neptune", jd_c)[:, :3], jd_c),
    )
    conj = {}
    for grid_name, step in (("coarse", COARSE_STEP_DAYS), ("fine", FINE_STEP_DAYS)):
        g = at[grid_name]
        jd = jd_all[g]
        c = conjunction_jd(
            jd, lon_of_date_deg(pos1[g, I_URANUS], jd), lon_of_date_deg(pos1[g, I_NEPTUNE], jd)
        )
        conj[f"grid_{int(step)}_day"] = {
            "grid": f"1810-01-01.0 to 1835-01-01.0 TDB, every {int(step)} d, linear interpolation",
            "jd_tdb": round(c, 4),
            "date_tdb": jd_to_date(c),
            "minus_jpl_days": round(c - jd_jpl, 3),
            "minus_jpl_hours": round((c - jd_jpl) * 24.0, 1),
            "minus_jpl_heliocentric_unrounded_hours": round((c - c_jpl_helio) * 24.0, 2) + 0.0,
        }
    jd_fine_all = jd_all[at["fine"]]
    g = at["fine"]
    acc_fine = neptune_accel(pos1[g, I_URANUS], pos1[g, I_NEPTUNE], gm_n)
    vhat_fine = vel1[g, I_URANUS] / np.linalg.norm(vel1[g, I_URANUS], axis=1, keepdims=True)
    along_fine = np.einsum("ni,ni->n", acc_fine, vhat_fine)
    jd_flip = first_sign_change(jd_fine_all, along_fine, conj["grid_1_day"]["jd_tdb"])

    # -- validation against Horizons (cache only) -----------------------------------------------
    validation = {
        "method": "run 1 integrated to each cached Horizons epoch exactly (no interpolation); "
        "angle between heliocentric position vectors",
    }
    for body, idx in (("uranus", I_URANUS), ("neptune", I_NEPTUNE)):
        jd = jd_val[body]
        ang = angle_arcsec(pos1[at[body], idx], horizons(body, jd)[:, :3])
        validation[body] = {
            "n_epochs": len(jd),
            "first_jd_tdb": float(jd[0]),
            "last_jd_tdb": float(jd[-1]),
            "last_date_tdb": jd_to_date(float(jd[-1])),
            "max_arcsec": round(float(ang.max()), 3),
            "max_at_date_tdb": jd_to_date(float(jd[np.argmax(ang)])),
            "final_arcsec": round(float(ang[-1]), 3),
        }
    for body in ("uranus", "neptune"):
        v = validation[body]
        v["reaches_j1850"] = bool(v["last_jd_tdb"] >= JD_END - FINE_STEP_DAYS)
        if not v["reaches_j1850"]:
            v["note"] = f"no cached epoch reaches J1850.0; latest is {v['last_date_tdb']} TDB"

    payload = {
        "description": "Uranus and Neptune, 1781-1850: REBOUND IAS15, Sun and eight planets "
        "from JPL Horizons states at Herschel's discovery (1781-03-13.0 TDB), run with Neptune's "
        "mass (run 1) and with Neptune's mass set to zero (run 2, Neptune a test particle). "
        "Positions heliocentric; the no-Neptune Uranus is stored as an offset from run 1.",
        "script": SCRIPT,
        "integrator": {
            "code": "REBOUND",
            "version": rebound.__version__,
            "scheme": "IAS15",
            "epsilon": diag1["epsilon"],
            "adaptive_mode": diag1["adaptive_mode"],
            "min_dt": diag1["min_dt"],
            "exact_finish_time": True,
            "G": bodies.G,
            "units": "AU, day, M_sun; G = k^2 (Gaussian k = 0.01720209895)",
            "frame": "centre of mass (sim.move_to_com) after adding heliocentric states",
            "run_1": diag1,
            "run_2": diag2,
        },
        "bodies": [
            {
                "name": "sun",
                "mass_msun": bodies.M_SUN_MSUN,
                "source": "1 M_sun by definition; inner planets are explicit bodies",
            }
        ]
        + [
            {
                "name": n,
                "horizons_id": bodies.BODIES[n].horizons_id,
                "mass_msun": 1.0 / MASS_RATIO[n],
                "sun_to_body_mass_ratio": MASS_RATIO[n],
                "mass_kind": MASS_KIND[n],
                "source": MASS_SOURCE,
            }
            for n in PLANETS
        ],
        "initial_states": {
            "source": "JPL Horizons (cache)",
            "center": "Sun body centre (500@10), geometric",
            "jd_tdb": A1_EPOCH_JD_TDB,
            "date": jd_to_date(A1_EPOCH_JD_TDB) + " TDB",
            "same_in_both_runs": True,
        },
        "frame": "heliocentric, ecliptic and mean equinox of J2000; longitudes for the "
        "conjunction precessed to ecliptic of date",
        "step_years": STEP_YEARS,
        "end": {"jd_tdb": JD_END, "label": "J1850.0", "last_frame_jd_tdb": float(jd_frames[-1])},
        "offset_exaggeration": K,
        "offset_exaggeration_rule": "smallest K in {1, 2, 5} x 10^k with 0.8 <= K * max "
        "|(dx, dy)| <= 2 AU",
        "summary": {
            "max_offset_xy_au": round(float(off_xy.max()), 9),
            "max_abs_dlon_arcsec": round(float(abs(dlon[i_dlon])), 2),
            "max_abs_dlon_t_year": round(float(jd_to_year(jd_frames[i_dlon])), 3),
            "dlon_end_arcsec": round(float(dlon[-1]), 2),
        },
        "conjunction": {
            "definition": "heliocentric Uranus-Neptune conjunction from run 1: first epoch "
            "1810-1835 where Uranus overtakes Neptune in ecliptic-of-date longitude, as in "
            "docs/notes/scripts/m6c_conjunction.py",
            **conj,
            "jpl": {
                "jd_tdb": jd_jpl,
                "date_tdb": jpl["jpl_uranus_neptune_conjunction_date_tdb"],
                "source": "site/data/m6c_conjunction.json epoch."
                "jpl_uranus_neptune_conjunction_jd_tdb (Horizons, inner "
                "barycentre, 10-day grid, rounded to 0.1 d)",
                "heliocentric_10_day_grid_unrounded_jd_tdb": round(c_jpl_helio, 4),
            },
        },
        "along_track_sign_change": {
            "jd_tdb": round(jd_flip, 2),
            "date_tdb": jd_to_date(jd_flip),
            "t_year": round(float(jd_to_year(jd_flip)), 4),
            "minus_conjunction_days": round(jd_flip - conj["grid_1_day"]["jd_tdb"], 1),
            "method": "along_track on the 1-day 1810-1835 grid (run 1), + to - crossing nearest "
            "the 1-day-grid conjunction, linear interpolation",
        },
        "validation": validation,
        "fields": {
            "t_year": "2000 + (jd_tdb - 2451545) / 365.25",
            "jd_tdb": "Julian date, TDB",
            "uranus_x, uranus_y, neptune_x, neptune_y": "run 1, heliocentric ecl_j2000, AU",
            "uranus_no_neptune_dx, uranus_no_neptune_dy": "run 2 Uranus minus run 1 Uranus, "
            "AU (draw at offset_exaggeration times)",
            "dlon_arcsec": "Uranus heliocentric ecliptic-J2000 longitude, run 1 minus run 2",
            "accel_along_track": "Neptune's heliocentric perturbing acceleration on Uranus "
            "(direct + indirect), component along Uranus's heliocentric velocity, AU/day^2; "
            "+ = drawn ahead, - = held back",
            "accel_radial": "same, component along Uranus's heliocentric position (+ = outward)",
            "accel_magnitude": "same, magnitude",
        },
        "frames": {
            "t_year": [round(float(v), 4) for v in jd_to_year(jd_frames)],
            "jd_tdb": [round(float(v), 4) for v in jd_frames],
            "uranus_x": [round(float(v), 5) for v in r_u[:, 0]],
            "uranus_y": [round(float(v), 5) for v in r_u[:, 1]],
            "neptune_x": [round(float(v), 5) for v in r_n[:, 0]],
            "neptune_y": [round(float(v), 5) for v in r_n[:, 1]],
            "uranus_no_neptune_dx": [round(float(v), 9) for v in off[:, 0]],
            "uranus_no_neptune_dy": [round(float(v), 9) for v in off[:, 1]],
            "dlon_arcsec": [round(float(v), 2) for v in dlon],
            "accel_along_track": [sig3(float(v)) for v in along],
            "accel_radial": [sig3(float(v)) for v in radial],
            "accel_magnitude": [sig3(float(v)) for v in mag],
        },
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False) + "\n")

    show = {
        k: payload[k]
        for k in (
            "offset_exaggeration",
            "summary",
            "conjunction",
            "along_track_sign_change",
            "validation",
        )
    }
    print(json.dumps(show, indent=1, ensure_ascii=False))
    print(
        f"frames {n_frames}, wrote {OUT_JSON} ({OUT_JSON.stat().st_size} bytes) "
        f"in {time.perf_counter() - t0:.1f} s"
    )


if __name__ == "__main__":
    main()
