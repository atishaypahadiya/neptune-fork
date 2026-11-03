"""Export site data to ``site/data/*.json``.

    uv run python -m fork.export.to_json

Currently exported:
- ``wh_convergence.json``: Uranus heliocentric-longitude error vs step size for both
  Wisdom–Holman compositions, 1781-03-13 → 1846-12-31, against an RK4 dt = 1 d reference
  (4th order; its own error is ~1e-10″). Feeds a future "integrator accuracy" panel.

M7 adds sightlines, residuals, sky map, degeneracy and energy exports here.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from fork.integrators import DEFAULT_DT_DAYS, rk4, wisdom_holman
from fork.integrators.wisdom_holman import SCHEME_DEFAULT, SCHEMES
from fork.model import bodies

SITE_DATA = Path(__file__).resolve().parents[3] / "site" / "data"
RAD_TO_ARCSEC = math.degrees(1.0) * 3600.0

A1_EPOCH_JD_TDB = 2371628.5  # 1781-03-13.0 TDB (Herschel's discovery; A1 initial conditions)
JD_1846_END_TDB = 2395296.5  # 1846-12-31.0 TDB
A1_BODIES = ("jupiter", "saturn", "uranus")
# Every dt divides 280 d, so the output grid below is on every integrator's step grid.
CONVERGENCE_DT_DAYS = (56.0, 40.0, 28.0, 20.0, 14.0, 12.0, 10.0, 7.0, 5.0, 3.5, 2.5)
CONVERGENCE_OUTPUT_STEP_DAYS = 280.0
REFERENCE_DT_DAYS = 1.0

SCHEME_LABELS = {
    "wh_drift_outer": "Wisdom–Holman, Kepler drift outermost (default; = REBOUND WHFast)",
    "wh_kick_outer": "Wisdom–Holman, kick outermost",
}


def _sig(x: float, digits: int = 4) -> float:
    return float(f"{x:.{digits}g}")


def a1_initial_conditions() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Heliocentric ecl_j2000 states of Jupiter, Saturn, Uranus at the A1 epoch (Horizons)."""
    from fork.data.horizons import state_vectors_ecl_j2000

    st = np.array([state_vectors_ecl_j2000(n, [A1_EPOCH_JD_TDB])[0] for n in A1_BODIES])
    masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in A1_BODIES])
    return st[:, :3].copy(), st[:, 3:].copy(), masses


def wh_convergence_data(
    x0: np.ndarray,
    v0: np.ndarray,
    masses: np.ndarray,
    dt_days: tuple[float, ...] = CONVERGENCE_DT_DAYS,
) -> dict:
    """Max |Δλ_Uranus| [″] and max |Δr_Uranus| [AU] over 1781–1846, per composition and dt."""
    t_out = np.arange(0.0, JD_1846_END_TDB - A1_EPOCH_JD_TDB, CONVERGENCE_OUTPUT_STEP_DAYS)
    ref = rk4.integrate(x0, v0, masses, t_out, REFERENCE_DT_DAYS)[:, 2]
    lon_ref_ecl_j2000_rad = np.unwrap(np.arctan2(ref[:, 1], ref[:, 0]))
    series = []
    for name, scheme in SCHEMES.items():
        dlon_arcsec, dr_au = [], []
        for dt in dt_days:
            pos = wisdom_holman.integrate(x0, v0, masses, t_out, dt, scheme)[:, 2]
            lon_ecl_j2000_rad = np.unwrap(np.arctan2(pos[:, 1], pos[:, 0]))
            dlon_arcsec.append(_sig(np.max(np.abs(lon_ecl_j2000_rad - lon_ref_ecl_j2000_rad))
                                    * RAD_TO_ARCSEC))  # fmt: skip
            dr_au.append(_sig(np.max(np.linalg.norm(pos - ref, axis=1))))
        series.append(
            {
                "scheme": name,
                "label": SCHEME_LABELS[name],
                "is_default": scheme == SCHEME_DEFAULT,
                "max_abs_dlon_uranus_arcsec": dlon_arcsec,
                "max_abs_dr_uranus_au": dr_au,
            }
        )
    return {
        "description": (
            "Uranus heliocentric ecliptic-J2000 longitude error vs Wisdom–Holman step size, "
            "Sun+Jupiter+Saturn+Uranus from JPL Horizons at 1781-03-13, integrated to "
            "1846-12-31. Error = max over outputs every 280 d of |WH − RK4(dt = 1 d)|."
        ),
        "span_jd_tdb": [A1_EPOCH_JD_TDB, JD_1846_END_TDB],
        "reference": {"integrator": "rk4", "dt_days": REFERENCE_DT_DAYS},
        "default_dt_days": DEFAULT_DT_DAYS,
        "observational_sigma_arcsec": 3.0,
        "dt_days": list(dt_days),
        "series": series,
    }


def write_json(name: str, payload: dict, directory: Path = SITE_DATA) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    return path


def main() -> None:
    x0, v0, masses = a1_initial_conditions()
    path = write_json("wh_convergence.json", wh_convergence_data(x0, v0, masses))
    print(f"wrote {path} ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
