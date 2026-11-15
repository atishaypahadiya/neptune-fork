"""Rebuild tests/fixtures/phase2/contact_real_horizons.npz: JPL Horizons states at the
published geocentric contacts of the 2016 and 2019 transits of Mercury (Espenak, eclipsewise.com
tm2016 and tm2019, Table 1, UT1), for tests/test_contact_real.py. Network, once.

For each of the eight contacts, at the published UT1 converted to TDB with the IERS values
bundled with astropy (no download):

- the Sun-centred geometric states of Mercury (199) and the Earth's centre (399), ecliptic
  and mean equinox of J2000, as ``fork.data.horizons`` fetches them;
- Mercury's astrometric (light-time corrected) position from the Earth's centre
  (``CENTER='500@399'``, ``VEC_CORR='LT'``), same frame, and its one-way light time.

Run from the repo root: ``uv run python tests/fixtures/build_contact_fixture.py``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "tests" / "fixtures" / "phase2" / "contact_real_horizons.npz"

# Espenak, "2016 Transit of Mercury" and "2019 Transit of Mercury", eclipsewise.com, Table 1
# (geocentric phases, UT1; ephemerides from JPL DE406)
PUBLISHED_UT1 = {
    "2016": ("2016-05-09", ("11:12:19", "11:15:31", "18:39:14", "18:42:26")),
    "2019": ("2019-11-11", ("12:35:27", "12:37:08", "18:02:33", "18:04:14")),
}


def main() -> None:
    from astropy.time import Time
    from astropy.utils import iers
    from astroquery.jplhorizons import Horizons

    from fork.data.horizons import fetch_from_horizons

    iers.conf.auto_download = False
    labels, ut1 = [], []
    for year, (date, times) in PUBLISHED_UT1.items():
        for k, hms in zip((1, 2, 3, 4), times, strict=True):
            labels.append(f"{year} {k}")
            ut1.append(f"{date} {hms}")
    t = Time(ut1, scale="ut1")
    jd_tdb = np.asarray(t.tdb.jd, dtype=np.float64)
    delta_t = np.asarray((t.tt.jd - t.jd) * 86400.0)
    mercury = fetch_from_horizons("199", jd_tdb)
    earth = fetch_from_horizons("399", jd_tdb)
    obj = Horizons(id="199", location="500@399", epochs=jd_tdb.tolist())
    tab = obj.vectors(refplane="ecliptic", aberrations="astrometric", cache=False)
    if not np.allclose(np.asarray(tab["datetime_jd"], float), jd_tdb, atol=1e-6, rtol=0):
        raise SystemExit("Horizons returned other epochs")
    astrometric = np.column_stack([np.asarray(tab[c], float) for c in ("x", "y", "z")])
    light_time = np.asarray(tab["lighttime"], float)  # days
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez(OUT, labels=np.array(labels), ut1=np.array(ut1), jd_tdb=jd_tdb,
             delta_t_iers_s=delta_t, mercury=mercury, earth=earth,
             mercury_astrometric_geocentric=astrometric, light_time_days=light_time)  # fmt: skip
    print(f"wrote {OUT.relative_to(ROOT)}: {len(labels)} contacts")


if __name__ == "__main__":
    main()
