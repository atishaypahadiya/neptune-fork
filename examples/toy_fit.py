"""DEMO, NOT A RESULT: one partition cell of the M5b posterior, on synthetic data, small nlive.

    uv run python examples/toy_fit.py                  # a few minutes on one core, offline
    uv run python examples/toy_fit.py --seed 3 --nlive 40

What it runs: the real pipeline (own Wisdom–Holman integrator, Uranus a, e, ϖ, λ₀ free, the
shared prior, dynesty with bound='multi', sample='rslice') on synthetic Uranus longitudes
(noise seed ``--seed``), restricted to ONE of the four M5b prior cells, a ∈ [28, 60] AU ×
e ∈ [0, 0.12) (the cell true Neptune lies in), with nlive 20 instead of 1000–2000 and a loose
stopping tolerance. It prints the cell's ln Z and the posterior 1846-09-23 longitude of the
perturber beside true Neptune's.

What it is not: the M5b result. With one cell the posterior ignores the other three cells'
mass, and with nlive 20 the interval and ln Z carry large sampling error that the quoted error
understates. Use it to check an install and to see the moving parts. The real numbers are in
``site/data/m5b_posteriors.json`` and ``docs/notes/M5.md``.

Data: by default the Horizons states come from the committed test fixture
(``tests/data/horizons_m4.json``), so the demo needs no network. ``--live-cache`` uses
``data/cache/`` instead (fill it with ``python -m fork.data.horizons --fetch``).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "data" / "horizons_m4.json"
EPOCH_JD_TDB = 2371628.5  # 1781-03-13, the A1 / M4 epoch
CELL = 2  # nested.PARTITION_CELLS[2]: a [28, 60] AU × e [0, 0.12)
E_MAX = 0.3
NLIVE = 20
DLOGZ = 0.5
BANNER = "DEMO, NOT A RESULT"


def use_fixture_cache() -> Path:
    """Point ``fork.data.horizons`` at a temporary cache built from the committed fixture."""
    from fork.data import horizons

    directory = Path(tempfile.mkdtemp(prefix="fork_toy_horizons_"))
    data = json.loads(FIXTURE.read_text())
    for body, rec in data["bodies"].items():
        horizons._save_cache(directory / f"horizons_{body}.npz", np.array(rec["jd_tdb"]),
                             np.array(rec["state_ecl_j2000"]),
                             horizons.BODIES[body].horizons_id)  # fmt: skip
    os.environ["FORK_CACHE_DIR"] = str(directory)
    return directory


def run(seed: int = 0, *, nlive: int = NLIVE, dlogz: float = DLOGZ, rng_seed: int = 1,
        maxiter: int | None = None, n_draws: int = 400) -> dict:  # fmt: skip
    """One nested run on cell ``CELL``; returns a summary dict (degrees, Earth masses, AU)."""
    from fork.data.synthetic import make_synthetic_observations
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.nested import ERR_INFLATION, PARTITION_CELLS, cell_space, run_nested
    from fork.inference.priors import UranusBox, free_a_space
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch
    from fork.model.hidden_mass import true_neptune_geocentric_longitude

    t0 = time.time()
    ctx = ForwardContext.from_horizons(EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(EPOCH_JD_TDB)  # held (rule 2); a, e, ϖ, λ₀ are free
    obs, _ = make_synthetic_observations(np.random.default_rng(seed))
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    space = cell_space(free_a_space(pq, E_MAX, box), PARTITION_CELLS[CELL])
    res = run_nested(obs, ctx, space, np.random.default_rng(rng_seed), nlive=nlive,
                     dlogz=dlogz, maxiter=maxiter)  # fmt: skip

    rng = np.random.default_rng([rng_seed, 1])
    eq = res.equal_weight(rng)
    eq = eq[rng.choice(len(eq), size=min(n_draws, len(eq)), replace=False)]
    geo = GalleGeometry.from_horizons()
    thetas = np.array([space.theta(q) for q in eq])
    lon = np.array([math.degrees(galle_lon_of_date_rad(th, ctx, geo)) for th in thetas])
    # unwrap around the circular median so the 16/84 % points are not split across 0°/360°
    ref = math.degrees(math.atan2(np.sin(np.radians(lon)).mean(), np.cos(np.radians(lon)).mean()))
    lon = ref + ((lon - ref + 180.0) % 360.0) - 180.0
    lo, med, hi = np.percentile(lon, [16, 50, 84])
    truth = math.degrees(true_neptune_geocentric_longitude().lon_ecl_of_date_rad)
    (a_lo, a_hi), e_rng = PARTITION_CELLS[CELL]
    return {
        "seed": seed, "cell": {"a_au": [a_lo, a_hi], "e": list(e_rng)},
        "cell_prior_weight": space.slice_weight, "nlive": nlive, "dlogz": dlogz,
        "niter": res.niter, "ncall": res.ncall, "wall_s": time.time() - t0,
        "logz_k": res.logz, "logzerr_k": ERR_INFLATION * res.logzerr,
        "chi2_min": -2.0 * res.max_logl, "n_obs": int(obs.jd_tdb.size), "n_draws": len(eq),
        "lon_1846_deg": {"median": med % 360.0, "lo68": lo % 360.0, "hi68": hi % 360.0},
        "lon_1846_true_deg": truth,
        "mass_earth_median": float(np.median(10.0 ** thetas[:, 6])),
        "a_au_median": float(np.median(thetas[:, 7])),
        "e_median": float(np.median(thetas[:, 8])),
    }  # fmt: skip


def report(s: dict) -> str:
    lon = s["lon_1846_deg"]
    half = 0.5 * ((lon["hi68"] - lon["lo68"]) % 360.0)
    miss = ((lon["median"] - s["lon_1846_true_deg"] + 180.0) % 360.0) - 180.0
    return "\n".join([
        f"{BANNER}: one prior cell, nlive {s['nlive']} (M5b uses four cells, nlive 1000-2000).",
        f"synthetic seed {s['seed']}, {s['n_obs']} longitudes; cell a {s['cell']['a_au']} AU, "
        f"e {s['cell']['e']} (prior weight {s['cell_prior_weight']:.3f})",
        f"sampler: {s['niter']} iterations, {s['ncall']} likelihood calls, "
        f"{s['wall_s'] / 60:.1f} min",
        f"cell ln Z_k = {s['logz_k']:.2f} +/- {s['logzerr_k']:.2f} (dynesty error x1.8); "
        f"best chi2 {s['chi2_min']:.1f}",
        f"perturber, posterior medians: m {s['mass_earth_median']:.1f} M_earth, "
        f"a {s['a_au_median']:.1f} AU, e {s['e_median']:.3f}",
        f"1846-09-23 geocentric longitude (ecliptic of date), {s['n_draws']} draws: "
        f"median {lon['median']:.1f} deg, 68% [{lon['lo68']:.1f}, {lon['hi68']:.1f}]",
        f"true Neptune {s['lon_1846_true_deg']:.2f} deg; median - truth {miss:+.1f} deg "
        f"(half-width {half:.1f} deg)",
        f"{BANNER}. Real numbers: site/data/m5b_posteriors.json, docs/notes/M5.md.",
    ])  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=0, help="synthetic noise seed (default 0)")
    ap.add_argument("--nlive", type=int, default=NLIVE)
    ap.add_argument("--dlogz", type=float, default=DLOGZ)
    ap.add_argument("--rng-seed", type=int, default=1, help="sampler rng seed")
    ap.add_argument("--live-cache", action="store_true",
                    help="use data/cache/ instead of the committed Horizons fixture")  # fmt: skip
    args = ap.parse_args(argv)
    if not args.live_cache:
        use_fixture_cache()
    print(f"{BANNER}. Running one prior cell at nlive {args.nlive}...", file=sys.stderr)
    s = run(args.seed, nlive=args.nlive, dlogz=args.dlogz, rng_seed=args.rng_seed)
    print(report(s))


if __name__ == "__main__":
    main()
