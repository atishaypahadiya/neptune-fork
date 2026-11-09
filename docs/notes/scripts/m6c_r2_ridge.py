"""After-the-fact diagnostic (not pre-registered): the e ≈ 0.06–0.12 ridge in the R2 runs.

    uv run python docs/notes/scripts/m6c_r2_ridge.py

The R2 runs (nlive 4000, 28 slices; M6.md "M6c cell-2 repeatability") reach a ridge at
e ≥ 0.06 that the headline run of cell 2 (a 28–60 AU, e < 0.12) hardly sampled. This script
characterises it and asks whether it is where Le Verrier-like solutions live.

1. From the merged R2 cell-2 posterior (``merge_runs`` of r2_1 and r2_2), equal-weight draws
   split into ridge (e ≥ 0.06) and core (e < 0.06), and the ridge into its two branches
   (a below and above 30.5 AU). For each part: posterior share, and the
   16/50/84 % of a, mass, e, the perturber's heliocentric longitude and its longitude of
   perihelion at 1847 Jan 1.0 (osculating heliocentric, ecliptic of date), and its Galle-night
   longitude.
2. From all of the R2 samples (dead points included): the lowest −2 ln L reached and the
   posterior mass by bin of a for e ≥ 0.06, and in a box around Le Verrier's elements. The
   lowest −2 ln L a sampler reached in a bin is an upper bound on the profile there, not the
   profile.
3. Le Verrier's elements from ``fork.model.hidden_mass.LE_VERRIER_1846`` (Krajnović 2016,
   A&G 57, 5.28, Table 2, after Grosser 1962; Encyclopædia Britannica 1911), with the often
   quoted mass 1/9300 M☉ (``LE_VERRIER_1846_MASS_ALT_MSUN``) beside Krajnović's 0.00011 M☉.

Reads data/cache/nested/m6c_repeat/r2_*.pkl. Output: docs/notes/M6c_r2_ridge.json.
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np

from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
from fork.export import m6c, m6c_repeat
from fork.inference.nested import NestedResult
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.bodies import M_CENTRAL_MSUN, M_EARTH_MSUN, G
from fork.model.forward import predict_states_ecl_j2000
from fork.model.hidden_mass import JD_1847_JAN1_TDB, LE_VERRIER_1846, LE_VERRIER_1846_MASS_ALT_MSUN

OUT = Path(__file__).resolve().parents[1] / "M6c_r2_ridge.json"
E_RIDGE = 0.06
N_DERIVED = 3000  # per part
A_EDGES = (28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 45, 60)
A_SPLIT = 30.5  # the ridge has two branches in a, with a gap at 30–31 AU
LV_BOX = {"a_au": (34.0, 38.0), "e": (0.08, 0.12)}
I_PERTURBER = 3


def q3(x, nd=3) -> list[float]:
    return [round(float(v), nd) for v in np.percentile(x, [16, 50, 84])]


def lon_q3(lon_deg: np.ndarray) -> list[float]:
    from fork.inference.mcmc import circular_interval_deg

    return [round(float(v), 2) for v in circular_interval_deg(lon_deg, (16.0, 50.0, 84.0))]


def state_1847(theta: np.ndarray, ctx) -> tuple[float, float]:
    """(heliocentric longitude, osculating longitude of perihelion), degrees, ecliptic of date
    at 1847 Jan 1.0, of the perturber."""
    pos, vel = predict_states_ecl_j2000(theta, [JD_1847_JAN1_TDB], ctx)
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(JD_1847_JAN1_TDB)
    r, v = rot @ pos[0, I_PERTURBER], rot @ vel[0, I_PERTURBER]
    mu = G * (M_CENTRAL_MSUN + 10.0 ** theta[6] * M_EARTH_MSUN)
    ecc = np.cross(v, np.cross(r, v)) / mu - r / np.linalg.norm(r)
    return (math.degrees(math.atan2(r[1], r[0])) % 360.0,
            math.degrees(math.atan2(ecc[1], ecc[0])) % 360.0)  # fmt: skip


def main() -> None:
    m6c.check_frozen("headline")
    _, ctx, _ = m6c_repeat.cell2_problem()
    geo = GalleGeometry.from_horizons()
    runs = [NestedResult.load(j.path(".pkl")) for j in m6c_repeat.JOBS if j.group == "r2"]
    merged = m6c_repeat.merged_cell2(runs)
    space = merged.space
    rng = np.random.default_rng([11, 1846, 6])

    eq = merged.equal_weight(rng)
    th = np.array([space.theta(q) for q in eq])
    ridge = th[:, 8] >= E_RIDGE
    parts = {"ridge_e_ge_0p06": ridge, "core_e_lt_0p06": ~ridge,
             "ridge_a_lt_30p5": ridge & (th[:, 7] < A_SPLIT),
             "ridge_a_ge_30p5": ridge & (th[:, 7] >= A_SPLIT)}  # fmt: skip
    out = {"description": __doc__.strip().splitlines()[0], "n_equal_weight": len(eq),
           "parts": {}}  # fmt: skip
    for name, sel in parts.items():
        idx = np.flatnonzero(sel)
        pick = idx[rng.choice(len(idx), size=min(N_DERIVED, len(idx)), replace=False)]
        galle = np.array([math.degrees(galle_lon_of_date_rad(th[i], ctx, geo)) for i in pick])
        hv = np.array([state_1847(th[i], ctx) for i in pick])
        out["parts"][name] = {
            "share": round(float(sel.mean()), 4), "n_derived": len(pick),
            "a_au_16_50_84": q3(th[sel, 7]), "e_16_50_84": q3(th[sel, 8], 4),
            "mass_earth_16_50_84": q3(10.0 ** th[sel, 6], 2),
            "helio_lon_1847_deg_16_50_84": lon_q3(hv[:, 0]),
            "varpi_1847_deg_16_50_84": lon_q3(hv[:, 1]),
            "galle_lon_deg_16_50_84": lon_q3(galle),
            "galle_lon_deg_2p5_97p5": [round(float(v), 2)
                                       for v in np.percentile(galle, [2.5, 97.5])],
        }  # fmt: skip
        print(name, json.dumps(out["parts"][name]), flush=True)

    res = merged.results
    m2 = -2.0 * np.asarray(res["logl"])
    w = np.asarray(res.importance_weights())
    ta = np.array([space.theta(q) for q in res["samples"]])
    ridge = ta[:, 8] >= E_RIDGE
    rows = []
    for lo, hi in itertools.pairwise(A_EDGES):
        s = ridge & (ta[:, 7] >= lo) & (ta[:, 7] < hi)
        rows.append({"a_lo": lo, "a_hi": hi, "posterior_mass": round(float(w[s].sum()), 5),
                     "min_m2lnl": round(float(m2[s].min()), 2) if s.any() else None,
                     "n_samples": int(s.sum())})  # fmt: skip
    box = ((ta[:, 7] >= LV_BOX["a_au"][0]) & (ta[:, 7] < LV_BOX["a_au"][1])
           & (ta[:, 8] >= LV_BOX["e"][0]) & (ta[:, 8] < LV_BOX["e"][1]))  # fmt: skip
    lv = LE_VERRIER_1846
    out |= {
        "cell2_min_m2lnl": round(float(m2.min()), 2),
        "ridge_by_a": rows,
        "le_verrier_box": LV_BOX | {"posterior_mass": round(float(w[box].sum()), 6),
                                    "min_m2lnl": round(float(m2[box].min()), 2) if box.any()
                                    else None, "n_samples": int(box.sum())},
        "le_verrier": {"a_au": lv.a_au, "e": lv.e,
                       "mass_earth": round(lv.mass_msun / M_EARTH_MSUN, 1),
                       "mass_earth_1_over_9300": round(LE_VERRIER_1846_MASS_ALT_MSUN
                                                       / M_EARTH_MSUN, 1),
                       "varpi_1847_deg": round(lv.varpi_ecl_of_date_deg, 3),
                       "helio_lon_1847_deg": round(lv.true_lon_ecl_of_date_deg, 3),
                       "galle_lon_deg_derived": round(m6c.le_verrier_galle_lon_deg(), 3),
                       "sources": list(lv.sources)},
    }  # fmt: skip
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("cell2_min_m2lnl", "le_verrier_box")}))
    for r in rows:
        print(r)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
