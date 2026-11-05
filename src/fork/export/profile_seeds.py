"""Posterior-seeded starts for the M5a free-a profile (from the M5b partitioned posteriors).

At each grid longitude, one extra start per M5b cell: the best-χ² posterior draw of that cell
whose 1846-09-23 longitude is within ``SEED_WINDOW_DEG`` of the grid point, or, if none is,
the draw nearest in longitude. Draws are ``DRAWS_PER_CELL`` weighted posterior draws per cell
(duplicates dropped), with their 1846 longitudes computed directly. Seeding is deterministic
(explicit rng per dataset and cell).

Added after the bounded-fit profile re-export (2026-11-04) missed the a ≈ 22 AU, e = 0.3 branch
at high longitude on six free-a seeds; the posterior has already found every branch with
weight, so it is the natural start set (with ``FREE_A_BACKSTOP_AU`` as a backstop).
"""

from __future__ import annotations

import math

import numpy as np

SEED_WINDOW_DEG = 2.0
DRAWS_PER_CELL = 2000
SEED_RUN = "free_a_e0.3"  # the M5b posterior that matches the profile's model (free a, e ≤ 0.3)


def cell_draws(seed: int | None, ctx, geo, run: str = SEED_RUN,
               n_draws: int = DRAWS_PER_CELL) -> list[dict]:  # fmt: skip
    """Per M5b cell of ``run`` for this dataset: θ (n, 11), χ² (n,) and 1846 longitude (rad)
    of up to ``n_draws`` distinct weighted posterior draws. Raises if a cell has no result."""
    from fork.export.m5b import jobs_for
    from fork.inference.nested import NestedResult
    from fork.inference.profile import galle_lon_of_date_rad

    jobs = [j for j in jobs_for(seed) if j.run == run and j.rep == 1]
    out = []
    for j in jobs:
        if not j.path(".pkl").exists():
            raise FileNotFoundError(f"no M5b result for {j.id}")
        r = NestedResult.load(j.path(".pkl"))
        rng = np.random.default_rng([11, 10_000 if seed is None else seed, j.cell + 1])
        w = r.weights / r.weights.sum()
        idx = np.unique(rng.choice(w.size, size=n_draws, p=w))
        theta = np.array([r.space.theta(v) for v in r.samples[idx]])
        chi2 = -2.0 * np.asarray(r.results["logl"])[idx]
        lon = np.array([galle_lon_of_date_rad(th, ctx, geo) for th in theta])
        out.append({"cell": j.cell, "theta": theta, "chi2": chi2, "lon_rad": lon})
    return out


Start = tuple[np.ndarray, float]  # (θ, its 1846 longitude in rad)


def seeded_starts(cells: list[dict], targets_rad: np.ndarray,
                  window_deg: float = SEED_WINDOW_DEG) -> list[list[Start]]:  # fmt: skip
    """Per target: one (θ, 1846 longitude) per cell, the best-χ² draw within ``window_deg``,
    else the nearest draw in longitude."""
    win = math.radians(window_deg)
    starts = []
    for t in np.asarray(targets_rad, dtype=np.float64):
        row = []
        for c in cells:
            d = np.abs((c["lon_rad"] - t + math.pi) % (2.0 * math.pi) - math.pi)
            near = np.flatnonzero(d <= win)
            k = int(near[np.argmin(c["chi2"][near])]) if near.size else int(np.argmin(d))
            row.append((c["theta"][k], float(c["lon_rad"][k])))
        starts.append(row)
    return starts
