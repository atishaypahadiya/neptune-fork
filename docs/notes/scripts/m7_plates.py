"""Geometry for the site's Plate I (Galle's night) and Plate II (the 1821 conjunction):
site/data/m7_plates.json.

Everything here is heliocentric (about the Sun + inner-planet barycentre, as in the model) and
projected onto the ecliptic of date of the plate's epoch. The x axis points to that date's
equinox, and longitudes are counted from it.

- **Posterior.** The 6000 derived samples of the M6c replacement posterior, drawn exactly as
  ``docs/notes/scripts/m6c_conjunction.py`` draws them (``m6c_repeat.REPORT_RNG``). The script
  checks that their conjunction longitudes reproduce ``site/data/m6c_conjunction.json`` (median,
  lo68, hi68 to 1e-3°) and their Galle-night longitudes reproduce the reported interval. A
  seeded subsample (``SUBSAMPLE`` samples, rng ``SUBSAMPLE_SEED``) is written as positions for
  the stipple.
- **Comparison bodies.** Uranus is the fitted one (the median of the samples' modelled positions;
  the cache has no Uranus state at these epochs). JPL's Neptune comes from the Horizons cache
  (read only; any cache miss raises instead of fetching). Le Verrier's planet comes from his
  1846 elements (``fork.model.hidden_mass.LE_VERRIER_1846``, cited there), propagated as in
  ``m6c.le_verrier_state_ecl_j2000``.
- **Orbits.** JPL's Uranus and Neptune paths over the epochs cached for them and for all inner
  planets, and Le Verrier's ellipse over one period, all projected onto the same plate frame.

Run from the repo root: ``uv run python docs/notes/scripts/m7_plates.py`` (a few minutes).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6c_conjunction as mc

import fork.data.horizons as hz
from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
from fork.export import m6c, m6c_repeat
from fork.export.to_json import SITE_DATA
from fork.inference.nested import ERR_INFLATION, PartitionedResult
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.bodies import INNER_BODIES
from fork.model.forward import predict_states_ecl_j2000
from fork.model.hidden_mass import I_PERTURBER, LE_VERRIER_1846
from fork.model.observation import I_URANUS

OUT = SITE_DATA / "m7_plates.json"
# Plate II magnifies angles about the line SC (radii stay to scale): the smallest factor in
# MAGNIFICATIONS that opens Le Verrier's angular offset from C to at least MIN_OPEN_DEG.
MAGNIFICATIONS = (1, 2, 5, 10, 20)
MIN_OPEN_DEG = 10.0
SUBSAMPLE = 600
SUBSAMPLE_SEED = 20261109


def offline(*_args: object, **_kwargs: object) -> np.ndarray:
    raise RuntimeError("m7_plates reads the Horizons cache only; a fetch was requested")


def about(body: str, jd: np.ndarray) -> np.ndarray:
    return hz.states_about_inner_barycenter_ecl_j2000(body, jd, fetcher=offline)


def cached_epochs(body: str) -> np.ndarray:
    """Epochs cached for ``body`` and for every inner planet (the inner-barycentre correction
    needs all of them), so a path is drawn in the model's frame without fetching."""

    def load(name: str) -> np.ndarray:
        return hz._load_cache(hz._cache_path(name, hz.cache_dir()))[0]

    jd = load(body)
    for inner in INNER_BODIES:
        other = load(inner.name)
        jd = jd[np.isin(np.round(jd, 6), np.round(other, 6))]
    return np.sort(jd)


def to_plate(pos_ecl_j2000: np.ndarray, jd_plate: float) -> np.ndarray:
    """(n, 2) x, y in AU on the ecliptic of date at ``jd_plate``."""
    m = precession_matrix_ecl_j2000_to_ecl_of_date(jd_plate)
    return (np.atleast_2d(pos_ecl_j2000) @ m.T)[:, :2]


def lon_deg(xy: np.ndarray) -> np.ndarray:
    return np.degrees(np.arctan2(xy[:, 1], xy[:, 0])) % 360.0


def linear_q(x: np.ndarray) -> dict:
    lo95, lo, med, hi, hi95 = np.percentile(x, (2.5, 16.0, 50.0, 84.0, 97.5))
    return {
        k: round(float(v), 3)
        for k, v in (("median", med), ("lo68", lo), ("hi68", hi), ("lo95", lo95), ("hi95", hi95))
    }


def r5(a: np.ndarray) -> list:
    return np.round(np.asarray(a, dtype=float), 5).tolist()


def main() -> None:
    conj = json.loads((SITE_DATA / "m6c_conjunction.json").read_text())
    jd_c = conj["epoch"]["jd_tdb"]
    jd_g = conj["galle_night_jd_tdb"]

    m6c.check_frozen("headline")
    rep = json.loads((SITE_DATA / m6c_repeat.OUT_NAME).read_text())
    _, ctx, *_ = m6c.dataset(variant="headline")
    geo = GalleGeometry.from_horizons()
    if abs(geo.jd_tdb - jd_g) > 1e-6:
        raise RuntimeError("Galle epoch differs from m6c_conjunction.json")

    res = m6c_repeat.replacement_cells()
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng(list(m6c_repeat.REPORT_RNG))
    eq = p.equal_weight(rng)
    idx = rng.choice(len(eq), size=min(m6c.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    th = np.array([res[0].space.theta(v) for v in eq[idx]])

    pos_c = np.empty((len(th), 3))
    pos_g = np.empty((len(th), 3))
    ura_c = np.empty((len(th), 3))
    ura_g = np.empty((len(th), 3))
    galle = np.empty(len(th))
    for i, t in enumerate(th):
        pos, _ = predict_states_ecl_j2000(t, [jd_c, jd_g], ctx)
        pos_c[i], pos_g[i] = pos[0, I_PERTURBER], pos[1, I_PERTURBER]
        ura_c[i], ura_g[i] = pos[0, I_URANUS], pos[1, I_URANUS]
        galle[i] = math.degrees(galle_lon_of_date_rad(t, ctx, geo)) % 360.0
        if i % 1000 == 0:
            print(f"samples {i}/{len(th)}", flush=True)

    xy_c, xy_g = to_plate(pos_c, jd_c), to_plate(pos_g, jd_g)
    lon_c = lon_deg(xy_c)
    qc, qg = mc.q(lon_c), mc.q(galle)
    want = conj["conjunction_helio"]["all"]
    if not all(abs(qc[k] - want[k]) < 1e-3 for k in ("median", "lo68", "hi68")):
        raise RuntimeError(f"conjunction longitudes {qc} do not reproduce {want}")
    if not all(abs(qg[k] - rep["reported"][k]) < 1e-3 for k in ("median", "lo68", "hi68")):
        raise RuntimeError(f"Galle-night longitudes {qg} do not reproduce {rep['reported']}")

    sub = np.random.default_rng(SUBSAMPLE_SEED).choice(len(th), size=SUBSAMPLE, replace=False)
    sub.sort()

    def bodies_at(jd: float, ura: np.ndarray) -> dict:
        # Uranus: the fitted one, the median of the posterior samples' modelled positions (the
        # Horizons cache has no Uranus state at these epochs, and nothing is fetched).
        u = np.median(to_plate(ura, jd), axis=0)
        n = to_plate(about("neptune", np.array([jd]))[:, :3], jd)[0]
        lv = to_plate(m6c.le_verrier_state_ecl_j2000(jd)[0], jd)[0]
        return {
            k: {
                "xy": r5(v),
                "lon_deg": round(float(lon_deg(v[None])[0]), 3),
                "r_au": round(float(np.hypot(*v)), 4),
            }
            for k, v in (("uranus_fit", u), ("neptune_jpl", n), ("le_verrier", lv))
        }

    def path(body: str, jd_plate: float) -> list:
        jd = cached_epochs(body)
        return r5(to_plate(about(body, jd)[:, :3], jd_plate))

    lv_period_days = 365.25 * LE_VERRIER_1846.a_au**1.5
    lv_jd = jd_c + np.linspace(0.0, lv_period_days, 361)
    lv_orbit = np.array([m6c.le_verrier_state_ecl_j2000(t)[0] for t in lv_jd])

    def plate(jd: float, xy: np.ndarray, ura: np.ndarray, date: str) -> dict:
        lon = lon_deg(xy)
        r = np.hypot(xy[:, 0], xy[:, 1])
        return {
            "jd_tdb": jd,
            "date_tdb": date,
            "bodies": bodies_at(jd, ura),
            "posterior": {
                "lon_deg": mc.q(lon),
                "r_au": linear_q(r),
                "subsample_xy": r5(xy[sub]),
            },
            "orbits": {
                "uranus_jpl": path("uranus", jd),
                "neptune_jpl": path("neptune", jd),
                "le_verrier": r5(to_plate(lv_orbit, jd)),
            },
        }

    lv_off = abs(conj["conjunction_helio"]["le_verrier"]["minus_median_deg"])
    k_mag = next(k for k in MAGNIFICATIONS if k * lv_off >= MIN_OPEN_DEG)

    payload = {
        "description": "Plate geometry for the site: Plate I (Galle's night, 1846 Sep 23) and "
        "Plate II (the 1821 Uranus-planet conjunction). Built by docs/notes/scripts/m7_plates.py.",
        "frame": "heliocentric (Sun + inner-planet barycentre), projected onto the ecliptic of "
        "date of each plate's epoch; x toward that date's equinox; AU",
        "posterior": "M6c replacement posterior, 6000 derived samples (as m6c_conjunction.py); "
        f"stipple: seeded subsample of {SUBSAMPLE} (rng {SUBSAMPLE_SEED})",
        "checks": "conjunction-longitude quantiles reproduce m6c_conjunction.json and Galle-night "
        "geocentric quantiles reproduce m6c_repeat.json's reported interval, to 1e-3 deg",
        "sources": {
            "jpl": "JPL Horizons (cache): Neptune's position and the Uranus and Neptune paths",
            "uranus_fit": "Uranus at the plate epoch: median of the posterior samples' modelled "
            "positions",
            "le_verrier": conj["le_verrier_sources"],
            "orbit_paths": "JPL paths over the cached epochs; Le Verrier's ellipse over one "
            "period from his 1846 elements",
        },
        "plate_ii_conjunction": {
            **plate(jd_c, xy_c, ura_c, conj["epoch"]["date_tdb"]),
            "angle_magnification": k_mag,
            "angle_magnification_rule": f"angles about SC drawn x K, radii to scale; K = smallest "
            f"of {list(MAGNIFICATIONS)} with K x |L - C| >= {MIN_OPEN_DEG} deg",
        },
        "plate_i_galle_night": plate(jd_g, xy_g, ura_g, mc.jd_to_date(jd_g)),
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False) + "\n")
    print("conjunction lon", qc)
    print("galle geo", qg)
    for k in ("plate_ii_conjunction", "plate_i_galle_night"):
        print(k, json.dumps(payload[k]["bodies"]), payload[k]["posterior"]["r_au"])
    print("bytes", OUT.stat().st_size)


if __name__ == "__main__":
    main()
