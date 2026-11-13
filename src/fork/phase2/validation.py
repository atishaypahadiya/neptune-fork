"""Validation V1–V5 of the Phase 2 forces (docs/phase2/validation.md). Engineering checks on
integrations against closed forms and against the JPL ephemeris; no observation is fitted.

    uv run python -m fork.phase2.validation run        # the long runs → docs/phase2/validation.json
    uv run python -m fork.phase2.validation plot       # figures under docs/phase2/

``run`` writes one JSON per check under ``data/cache/phase2/validation/`` as it finishes and
skips checks already there (``--fresh`` recomputes), logging to ``progress.log`` beside them.
The tests call the same functions on short spans.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from fork.data.frames import elements_to_state
from fork.integrators import wh_extra as W
from fork.integrators._nbody import energy_helio, helio_to_bary
from fork.model import inner_system as IS
from fork.model.bodies import MASS_RATIO_MERCURY, G, masses_msun
from fork.phase2 import reference as R
from fork.phase2 import secular as S

ROOT = Path(__file__).resolve().parents[3]
OUT_JSON = ROOT / "docs" / "phase2" / "validation.json"
CACHE = R.CACHE / "validation"
LOG = CACHE / "progress.log"

# Mercury's mean elements at J2000 (Standish's table in the Explanatory Supplement, 1992,
# Table 5.8.1): the "elements used" of V1–V3.
A_MERCURY, E_MERCURY = 0.38709893, 0.20563069
M_MERCURY = 1.0 / MASS_RATIO_MERCURY
J2_SUN = 2.2e-7
C = S.ARCSEC_PER_CENTURY
YEAR = 365.25

TOL_V1 = 0.1  # ″/century
TOL_V2 = 0.05
TOL_V3 = 0.05


def _two_body(e: float = E_MERCURY, inc: float = 0.0, raan: float = 0.5, argp: float = 1.0):
    masses = np.array([1.0, M_MERCURY])
    mu = G * (1.0 + M_MERCURY)
    r, v = elements_to_state(A_MERCURY, e, inc, raan, argp, 0.3, mu)
    return r[None, :].copy(), v[None, :].copy(), masses, mu


def _advance(extra: np.ndarray, years: float, *, dt: float = 1.0, order: int = 2, **orbit):
    """Slopes of ω and Ω [rad/day] of the two-body orbit under ``extra``, minus the same
    without it."""
    x0, v0, masses, mu = _two_body(**orbit)
    t = np.arange(0.0, years * YEAR, 3.3)
    out = []
    for ex in (extra, W.make_extra()):
        x, v = W.integrate_states(x0, v0, masses, t, dt, ex, order)
        el = S.elements_series(x[:, 0], v[:, 0], mu)
        out.append((S.slope(t, el[:, 4]), S.slope(t, el[:, 3])))
    return out[0][0] - out[1][0], out[0][1] - out[1][1]


def v1(years: float = 1000.0) -> dict:
    """The 1PN advance of Mercury's perihelion, Sun and Mercury only."""
    rows = []
    for gr, beta, gamma in ((1.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 1.2, 1.0), (1.0, 1.0, 0.8)):
        w, o = _advance(W.make_extra(gr=gr, beta=beta, gamma=gamma), years)
        ana = S.gr_rate(A_MERCURY, E_MERCURY, G, gr=gr, beta=beta, gamma=gamma) * C
        rows.append(
            {
                "gr": gr,
                "beta": beta,
                "gamma": gamma,
                "measured": (w + o) * C,
                "analytic": ana,
                "difference": (w + o) * C - ana,
            }
        )
    return {
        "years": years,
        "a_au": A_MERCURY,
        "e": E_MERCURY,
        "dt_days": 1.0,
        "order": 2,
        "unit": "arcsec per Julian century",
        "tolerance": TOL_V1,
        "rows": rows,
        "pass": all(abs(r["difference"]) <= TOL_V1 for r in rows),
    }


def v2(years: float = 1000.0) -> dict:
    """The J2 advance, at Mercury's inclination to the solar equator and at 30°."""
    rows = []
    for inc_deg in (3.38, 30.0):
        inc = math.radians(inc_deg)
        w, o = _advance(W.make_extra(j2=J2_SUN), years, inc=inc)
        aw, ao = S.j2_rates(A_MERCURY, E_MERCURY, inc, G, J2_SUN)
        rows.append(
            {
                "inclination_deg": inc_deg,
                "measured_argp": w * C,
                "analytic_argp": aw * C,
                "measured_node": o * C,
                "analytic_node": ao * C,
                "ratio_argp": w / aw,
                "ratio_node": o / ao,
            }
        )
    ok = all(abs(r[k] - 1.0) <= TOL_V2 for r in rows for k in ("ratio_argp", "ratio_node"))
    return {
        "years": years,
        "j2": J2_SUN,
        "r_sun_au": W.R_SUN_AU,
        "tolerance": TOL_V2,
        "unit": "arcsec per Julian century",
        "rows": rows,
        "pass": ok,
    }


def v3(years: float = 1000.0, point_planet: bool = True) -> dict:
    """The advance from an interior ring (and a point planet on a circle) in the orbit plane,
    at two masses and two radii, against the orbit-averaged rate in the ring's field; the
    lowest-order Laplace–Lagrange rate is checked on a nearly circular orbit (e = 0.02)."""
    rows = []
    for mass in (1e-7, 3e-7):
        for radius in (0.15, 0.25):
            ex = W.make_extra(ring_mass=mass, ring_radius_au=radius)
            w, o = _advance(ex, years)
            wl, ol = _advance(ex, years, e=0.02)
            row = {
                "mass_msun": mass,
                "radius_au": radius,
                "ring_measured": (w + o) * C,
                "ring_analytic": S.ring_gauss_rate(A_MERCURY, E_MERCURY, G, mass, radius) * C,
                "ring_measured_e002": (wl + ol) * C,
                "laplace_lagrange": S.laplace_lagrange_rate(A_MERCURY, G, mass, radius) * C,
            }
            row["ratio_ring"] = row["ring_measured"] / row["ring_analytic"]
            row["ratio_laplace_lagrange_e002"] = row["ring_measured_e002"] / row["laplace_lagrange"]
            if point_planet:
                row["planet_measured"] = _planet_advance(mass, radius, years) * C
                row["ratio_planet"] = row["planet_measured"] / row["ring_analytic"]
            rows.append(row)
    keys = [k for k in rows[0] if k.startswith("ratio")]
    return {
        "years": years,
        "tolerance": TOL_V3,
        "unit": "arcsec per Julian century",
        "rows": rows,
        "pass": all(abs(r[k] - 1.0) <= TOL_V3 for r in rows for k in keys),
    }


def _planet_advance(mass: float, radius: float, years: float) -> float:
    """Mercury's ϖ slope with a point planet of ``mass`` on a coplanar circle at ``radius``,
    minus the two-body slope (order 4, dt 0.25 d: the planet's period is 21 to 46 d)."""
    x0, v0, masses, mu = _two_body(inc=0.0, raan=0.0, argp=1.5)
    speed = math.sqrt(G * (1.0 + mass) / radius)
    xp, vp = np.array([[radius, 0.0, 0.0]]), np.array([[0.0, speed, 0.0]])
    t = np.arange(0.0, years * YEAR, 3.3)
    x, v = W.integrate_states(
        np.vstack([x0, xp]),
        np.vstack([v0, vp]),
        np.append(masses, mass),
        t,
        0.25,
        W.make_extra(),
        4,
    )
    xn, vn = W.integrate_states(x0, v0, masses, t, 0.25, W.make_extra(), 4)
    return S.slope(
        t,
        S.perihelion_longitude(x[:, 0], v[:, 0], mu)
        - S.perihelion_longitude(xn[:, 0], vn[:, 0], mu),
    )


def _lon(x: np.ndarray) -> np.ndarray:
    return np.arctan2(x[:, 1], x[:, 0])


def v4(
    ref: dict[str, np.ndarray] | None = None,
    dt: float = IS.DEFAULT_DT_DAYS,
    order: int = IS.DEFAULT_ORDER,
) -> dict:
    """Sun and eight planets from the Horizons states at 1900.0, Mercury's heliocentric
    longitude against Horizons 1800–2000, with the 1PN term off, on, and on with J2."""
    ref = R.load() if ref is None else ref
    masses = np.concatenate([[1.0], masses_msun(R.PLANETS)])
    x0, v0 = ref["epoch_states"][:, :3].copy(), ref["epoch_states"][:, 3:].copy()
    t = ref["jd"] - R.EPOCH_JD
    mu = G * (1.0 + masses[1])
    year = 1900.0 + t / YEAR
    truth = ref["mercury"]
    w_true = S.perihelion_longitude(truth[:, :3], truth[:, 3:], mu)
    cases = {
        "newtonian": W.make_extra(),
        "with_1pn": W.make_extra(gr=1.0),
        "with_1pn_and_j2": W.make_extra(
            gr=1.0, j2=J2_SUN, spin_axis=tuple(IS.solar_spin_axis_ecl_j2000())
        ),
    }
    out: dict = {
        "epoch_jd_tdb": R.EPOCH_JD,
        "n_epochs": int(t.size),
        "dt_days": dt,
        "order": order,
        "unit": "arcsec",
        "cases": {},
    }
    for name, ex in cases.items():
        x, v = W.integrate_states(x0, v0, masses, t, dt, ex, order)
        dlon = (
            np.degrees((_lon(x[:, 0]) - _lon(truth[:, :3]) + math.pi) % (2 * math.pi) - math.pi)
            * 3600.0
        )
        dw = S.perihelion_longitude(x[:, 0], v[:, 0], mu) - w_true
        decades = []
        for y0 in range(1800, 2000, 10):
            sel = (year >= y0) & (year < y0 + 10)
            if sel.any():
                decades.append(
                    {
                        "decade": y0,
                        "n": int(sel.sum()),
                        "mean": float(dlon[sel].mean()),
                        "rms": float(np.sqrt((dlon[sel] ** 2).mean())),
                        "max_abs": float(np.abs(dlon[sel]).max()),
                    }
                )
        out["cases"][name] = {
            "max_abs_dlon": float(np.abs(dlon).max()),
            "perihelion_residual_slope_arcsec_per_century": S.slope(t, dw) * C,
            "decades": decades,
            "series": {"year": year.round(4).tolist(), "dlon": dlon.round(5).tolist()},
        }
    return out


def _pn_energy_and_j(x: np.ndarray, v: np.ndarray, mu: float) -> tuple[np.ndarray, np.ndarray]:
    """Test-particle 1PN energy and angular momentum (β = γ = 1), per unit mass, from the
    Lagrangian L = v²/2 + μ/r + [v⁴/8 + 3μv²/(2r) − μ²/(2r²)]/c²."""
    c2 = W.C_AU_PER_DAY**2
    r = np.linalg.norm(x, axis=1)
    v2 = (v * v).sum(axis=1)
    e = 0.5 * v2 - mu / r + (0.375 * v2 * v2 + 1.5 * mu * v2 / r + 0.5 * mu * mu / (r * r)) / c2
    j = np.cross(x, v) * (1.0 + (0.5 * v2 + 3.0 * mu / r) / c2)[:, None]
    return e, j


def v5_conservation(
    years: float = 1000.0, dt: float = IS.DEFAULT_DT_DAYS, order: int = IS.DEFAULT_ORDER
) -> dict:
    """Relative drift of the conserved quantities with each term on.

    1PN: Sun and a massless Mercury, the test-particle 1PN energy and angular momentum.
    J2, ring, both: Sun and eight planets, total energy including the terms' potentials, and
    the total angular momentum about the common axis (the z axis here)."""
    t = np.linspace(0.0, years * YEAR, 2001)
    out: dict = {"years": years, "dt_days": dt, "order": order, "rows": []}
    x0, v0, _, _ = _two_body()
    masses = np.array([1.0, 0.0])
    x, v = W.integrate_states(x0, v0, masses, t, dt, W.make_extra(gr=1.0), order)
    e, j = _pn_energy_and_j(x[:, 0], v[:, 0], G)
    jn = np.linalg.norm(j, axis=1)
    out["rows"].append(
        {
            "term": "1pn",
            "system": "Sun and a massless Mercury",
            "max_rel_energy": float(np.abs(e / e[0] - 1.0).max()),
            "max_rel_angular_momentum": float(np.abs(jn / jn[0] - 1.0).max()),
        }
    )
    ref = R.load(R.FIXTURE)
    m8 = np.concatenate([[1.0], masses_msun(R.PLANETS)])
    x0, v0 = ref["epoch_states"][:, :3].copy(), ref["epoch_states"][:, 3:].copy()
    cases = {
        "none": W.make_extra(),
        "j2": W.make_extra(j2=J2_SUN),
        "ring": W.make_extra(ring_mass=1e-7, ring_radius_au=0.2),
        "j2_and_ring": W.make_extra(j2=J2_SUN, ring_mass=1e-7, ring_radius_au=0.2),
    }
    for name, ex in cases.items():
        x, v = W.integrate_states(x0, v0, m8, t, dt, ex, order)
        en = np.empty(t.size)
        lz = np.empty(t.size)
        for i in range(t.size):
            en[i] = energy_helio(x[i], v[i], m8)
            for k in range(x.shape[1]):
                en[i] += m8[k + 1] * W.position_potential(x[i, k], G, ex)
            xb, vb = helio_to_bary(x[i], v[i], m8)
            lz[i] = (m8[:, None] * np.cross(xb, vb)).sum(axis=0)[2]
        out["rows"].append(
            {
                "term": name,
                "system": "Sun and eight planets",
                "max_rel_energy": float(np.abs(en / en[0] - 1.0).max()),
                "max_rel_angular_momentum": float(np.abs(lz / lz[0] - 1.0).max()),
            }
        )
    return out


def v5_step(years: float = 100.0) -> dict:
    """Mercury's position against a reference (order 4, dt 0.125 d) after ``years`` forward
    and back from 1900, eight planets, each term on: the largest difference, in arcseconds
    as seen from the Sun, by step and order."""
    ref = R.load(R.FIXTURE)
    m8 = np.concatenate([[1.0], masses_msun(R.PLANETS)])
    x0, v0 = ref["epoch_states"][:, :3].copy(), ref["epoch_states"][:, 3:].copy()
    t = np.linspace(-years * YEAR, years * YEAR, 401)
    cases = {
        "newtonian": W.make_extra(),
        "1pn": W.make_extra(gr=1.0),
        "j2": W.make_extra(j2=J2_SUN),
        "ring": W.make_extra(ring_mass=1e-7, ring_radius_au=0.2),
    }
    grid = [(2, 1.0), (2, 0.25), (2, 0.0625), (4, 8.0), (4, 4.0), (4, 2.0), (4, 1.0), (4, 0.5)]
    out: dict = {
        "years_each_way": years,
        "reference": {"order": 4, "dt_days": 0.125},
        "unit": "arcsec",
        "rows": [],
    }
    for name, ex in cases.items():
        xr, _ = W.integrate_states(x0, v0, m8, t, 0.125, ex, 4)
        for order, dt in grid:
            t0 = time.perf_counter()
            x, _ = W.integrate_states(x0, v0, m8, t, dt, ex, order)
            sec = time.perf_counter() - t0
            d = np.linalg.norm(x[:, 0] - xr[:, 0], axis=1) / np.linalg.norm(xr[:, 0], axis=1)
            out["rows"].append(
                {
                    "term": name,
                    "order": order,
                    "dt_days": dt,
                    "max_arcsec": float(np.degrees(d.max()) * 3600.0),
                    "seconds": round(sec, 3),
                }
            )
    return out


CHECKS = {
    "v1": v1,
    "v2": v2,
    "v3": v3,
    "v4": v4,
    "v5_conservation": v5_conservation,
    "v5_step": v5_step,
}


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(f"{datetime.now().isoformat(timespec='seconds')} {msg}\n")


def run(fresh: bool = False) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    result = {}
    for name, fn in CHECKS.items():
        path = CACHE / f"{name}.json"
        if path.exists() and not fresh:
            _log(f"{name}: already done, kept")
        else:
            _log(f"{name}: start")
            t0 = time.perf_counter()
            res = fn()
            res["seconds"] = round(time.perf_counter() - t0, 1)
            path.write_text(json.dumps(res), encoding="utf-8")
            _log(f"{name}: done in {res['seconds']} s")
        result[name] = json.loads(path.read_text(encoding="utf-8"))
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    _log(f"wrote {OUT_JSON}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--fresh", action="store_true")
    sub.add_parser("plot")
    a = ap.parse_args()
    if a.cmd == "run":
        run(a.fresh)
    else:
        from fork.phase2.validation_plot import plot

        plot()


if __name__ == "__main__":
    main()
