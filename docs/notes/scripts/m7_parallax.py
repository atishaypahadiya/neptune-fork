"""Size of the diurnal parallax the observation model leaves out, for the site's Limitations
(site/data/m7_parallax.json).

Uranus's equatorial horizontal parallax, asin(R_earth / Δ), with Δ the geocentric distance from
JPL Horizons states (cache only; Uranus and the Earth at their shared cached epochs over
1781–1845). The parallax an observer sees is at most this (it scales with the cosine of the
altitude), so the maximum over the span bounds the omission.

Run from the repo root: ``uv run python docs/notes/scripts/m7_parallax.py``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import fork.data.horizons as hz

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "site" / "data" / "m7_parallax.json"
R_EARTH_KM = 6378.137  # IERS 2010 / WGS84 equatorial radius
AU_KM = 149597870.7  # IAU 2012
JD_1781, JD_1846 = 2371185.5, 2395297.5


def main() -> None:
    ju, su = hz._load_cache(hz._cache_path("uranus", hz.cache_dir()))
    je, se = hz._load_cache(hz._cache_path("earth", hz.cache_dir()))
    ku = {round(float(j), 4): k for k, j in enumerate(ju)}
    pairs = [
        (ku[round(float(j), 4)], k)
        for k, j in enumerate(je)
        if round(float(j), 4) in ku and JD_1781 <= j <= JD_1846
    ]
    iu = np.array([a for a, _ in pairs])
    ie = np.array([b for _, b in pairs])
    delta = np.linalg.norm(su[iu, :3] - se[ie, :3], axis=1)
    par = np.degrees(np.arcsin(R_EARTH_KM / AU_KM / delta)) * 3600.0
    k = int(par.argmax())
    payload = {
        "description": "Diurnal parallax of Uranus, not modelled by the observation model: "
        "equatorial horizontal parallax asin(R_earth / Delta) over 1781-1845, from JPL Horizons "
        "states (cache). Built by docs/notes/scripts/m7_parallax.py.",
        "n_epochs": len(pairs),
        "max_arcsec": round(float(par.max()), 3),
        "min_arcsec": round(float(par.min()), 3),
        "max_at_jd_tdb": round(float(ju[iu][k]), 4),
        "min_geocentric_distance_au": round(float(delta.min()), 3),
        "constants": {"r_earth_km": R_EARTH_KM, "au_km": AU_KM},
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(payload)


if __name__ == "__main__":
    main()
