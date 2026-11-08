"""M6c S1–S3 diagnostics from the finished runs (no new sampling).

    uv run python docs/notes/scripts/m6c_cell2_repeatability.py

1. Per cell and variant: the Galle-night longitude (2.5/16/50/84/97.5 %) and a, e, mass
   (16/50/84 %) from that cell's own equal-weight samples.
2. Cell 2 (a 28–60 AU, e < 0.12) has the same prior and likelihood in the headline and in S1,
   so the two runs sample one distribution with different seeds. For each variant's cell 2:
   the run's own error on the median and 68 % width (200 ``resample_run`` and 200
   ``jitter_run`` draws), and, headline against S1, the posterior mass and the lowest
   −2 ln L reached in bins of a and of e.
3. The Uranus box edge check (as in M6.md, headline) on each variant's summary samples.

Reads data/cache/nested/m6c/*.pkl. Output: docs/notes/M6c_cell2_repeatability.json.
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np
from dynesty.utils import jitter_run, resample_run

from fork.export import m6c
from fork.inference.fit import uranus_theta_to_x
from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad

OUT = Path(__file__).resolve().parents[1] / "M6c_cell2_repeatability.json"
N_PER_CELL = 2000
N_BOOT = 200
LOGL_WINDOW = 30.0  # longitudes computed for samples within this of max ln L (weight beyond: ~1e-9)
A_EDGES = (28, 29, 29.5, 30, 30.5, 31, 31.5, 32, 33, 35, 40, 60)
E_EDGES = (0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12)


def wquantiles(x: np.ndarray, w: np.ndarray, q) -> np.ndarray:
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return np.interp(q, c / c[-1], x[o])


def r3(x) -> list[float]:
    return [round(float(v), 3) for v in x]


def load(variant: str) -> list[NestedResult]:
    n = len(m6c.variant_cells(variant))
    return [NestedResult.load(m6c.Job(variant, k).path(".pkl")) for k in range(n)]


def per_cell(variant: str, res: list[NestedResult], ctx, geo) -> list[dict]:
    space = res[0].space
    rows = []
    for k, r in enumerate(res):
        rng = np.random.default_rng([99, m6c.VARIANTS.index(variant), k])
        eq = r.equal_weight(rng)
        th = np.array([space.theta(q) for q in eq])
        pick = eq[rng.choice(len(eq), size=min(N_PER_CELL, len(eq)), replace=False)]
        g = np.array([math.degrees(galle_lon_of_date_rad(space.theta(q), ctx, geo)) for q in pick])
        p = np.percentile(g, [2.5, 16, 50, 84, 97.5])
        rows.append({
            "cell": k, "galle_2p5_16_50_84_97p5": r3(p), "galle_width68": round(p[3] - p[1], 3),
            "a_16_50_84": r3(np.percentile(th[:, 7], [16, 50, 84])),
            "e_16_50_84": [round(float(v), 4) for v in np.percentile(th[:, 8], [16, 50, 84])],
            "mass_16_50_84": r3(np.percentile(10.0 ** th[:, 6], [16, 50, 84])),
        })  # fmt: skip
    return rows


def cell2_internal_error(variant: str, r: NestedResult, ctx, geo) -> tuple[dict, dict]:
    res = r.results
    logl = np.asarray(res["logl"])
    assert np.all(np.diff(logl) >= 0)  # a static run's dead + live points, in ln L order
    keep = logl > logl.max() - LOGL_WINDOW
    g = np.full(len(logl), np.nan)
    for i in np.flatnonzero(keep):
        g[i] = math.degrees(galle_lon_of_date_rad(r.space.theta(res["samples"][i]), ctx, geo))
    w0 = np.asarray(res.importance_weights())
    q0 = wquantiles(g[keep], w0[keep], [0.16, 0.5, 0.84])
    rng = np.random.default_rng([5, m6c.VARIANTS.index(variant)])
    out = {"median": round(q0[1], 3), "width68": round(q0[2] - q0[0], 3),
           "weight_outside_window": float(w0[~keep].sum())}  # fmt: skip
    for name, f in (("resample_run", resample_run), ("jitter_run", jitter_run)):
        draws = []
        for _ in range(N_BOOT):
            rb = f(res, rstate=rng)
            idx = np.searchsorted(logl, rb["logl"])
            gb, wb = g[idx], np.asarray(rb.importance_weights())
            m = ~np.isnan(gb)
            q = wquantiles(gb[m], wb[m], [0.16, 0.5, 0.84])
            draws.append((q[1], q[2] - q[0]))
        d = np.array(draws)
        out[name] = {"std_median": round(float(d[:, 0].std()), 3),
                     "std_width68": round(float(d[:, 1].std()), 3)}  # fmt: skip
    th = np.array([r.space.theta(q) for q in res["samples"]])
    m2 = -2.0 * logl
    shape = {
        "min_m2lnl": round(float(m2.min()), 2),
        "posterior_mean_m2lnl_minus_min": round(float(np.sum(w0 * m2) - m2.min()), 2),
    }

    def bins(col: int, edges) -> list[dict]:
        rows = []
        for lo, hi in itertools.pairwise(edges):
            s = (th[:, col] >= lo) & (th[:, col] < hi)
            rows.append({"lo": lo, "hi": hi, "posterior_mass": round(float(w0[s].sum()), 4),
                         "min_m2lnl": round(float(m2[s].min()), 2) if s.any() else None,
                         "n_samples": int(s.sum())})  # fmt: skip
        return rows

    shape["a_bins"] = bins(7, A_EDGES)
    shape["e_bins"] = bins(8, E_EDGES)
    return out, shape


def box_edge(variant: str, res: list[NestedResult], box) -> dict:
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    eq = p.equal_weight(np.random.default_rng([*m6c.SUMMARY_RNG, m6c.VARIANTS.index(variant)]))
    space = res[0].space
    x = np.array([uranus_theta_to_x(space.theta(q)[:6])[:4] for q in eq])
    z = (x - np.asarray(box.centre)) / np.asarray(box.sigma_scaled)
    med, sd = np.median(z, axis=0), z.std(axis=0)
    edge = 50.0 - np.abs(med)
    return {"n_samples": len(eq), "chi2_per_dof": round(box.chi2_per_dof, 3),
            "order": ["a", "k", "h", "lambda0"], "median": r3(med), "std": r3(sd),
            "median_to_edge": r3(edge), "in_std": [round(float(v), 1) for v in edge / sd],
            "min": r3(z.min(axis=0)), "max": r3(z.max(axis=0)),
            "closest_sample_to_edge": round(float((50.0 - np.abs(z)).min()), 2)}  # fmt: skip


def main() -> None:
    geo = GalleGeometry.from_horizons()
    out = {"description": __doc__.strip().splitlines()[0], "variants": {}}
    for v in m6c.VARIANTS:
        m6c.check_frozen(v)
        res = load(v)
        _, ctx, _, _, box = m6c.dataset(variant=v)
        err, shape = cell2_internal_error(v, res[2], ctx, geo)
        out["variants"][v] = {"cells": per_cell(v, res, ctx, geo), "cell2_run_error": err,
                              "cell2_shape": shape, "box_edge": box_edge(v, res, box)}  # fmt: skip
        print(v, json.dumps(out["variants"][v]["cell2_run_error"]), flush=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
