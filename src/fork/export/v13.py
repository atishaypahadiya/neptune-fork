"""Code for the v1.3 draft plans (docs/notes/V13_drafts.md): M6d's shape-model space and
Peirce ratio, the M6c data choices, M6e's predictive step, the sub-cells of S1's cell 0.

Building blocks only. Nothing here is a runner for real data: no job list, no cache, no
ledger rows. The drafts are not frozen, and until they are these functions are exercised on
synthetic data alone (``synthetic_dataset``; tests/test_v13.py). No function reads a
Neptune state or a precovery position.

    uv run python -m fork.export.v13 smoke      # synthetic smoke runs and cost per call
"""

from __future__ import annotations

import argparse
import math
import time

import numpy as np

A1_EPOCH_JD_TDB = 2371628.5  # 1781-03-13.0 TDB, as fork.export.to_json
GALLE_LON_DEG = 325.873  # comparison target of M6d, as fork.export.m6c.JPL_GALLE_LON_DEG
PEIRCE_HALF_WIDTH_DEG = 0.87  # Le Verrier's 52′ miss (spec, M6d)
LE_VERRIER_A_AU = 36.154  # his elements of 31 August 1846; the draft leaves this value open
E_MAX = 0.3
SHAPE_E_CELLS = ((0.0, 0.12), (0.12, 0.3))  # the held-distance layout of Phase 1b step 2

DATA_CHOICES = ("leverrier", "bouvard", "adams")

# S1 cell 0 is a 20–25 AU, e ≤ 0.5. The draft's example cut, six sub-cells; open.
S1_E_MAX = 0.5
S1_CELL_0 = ((20.0, 25.0), (0.0, S1_E_MAX))
S1_SUBCELLS = tuple(((a0, a1), (e0, e1))
                    for a0, a1 in ((20.0, 22.5), (22.5, 25.0))
                    for e0, e1 in ((0.0, 0.2), (0.2, 0.35), (0.35, 0.5)))  # fmt: skip

DATA_RNG_TAG = 20261112  # synthetic data of this module only


# ------------------------------------------------------------------ M6c: data choices


def select_rows(rows: list[dict[str, str]], choice: str) -> list[dict[str, str]]:
    """Rows of ``leverrier_1846.csv`` under a data choice of the spec's M6c table:
    ``leverrier`` all; ``bouvard`` the modern observations only (1781–1845); ``adams`` all
    but Flamsteed 1690."""
    if choice == "leverrier":
        return list(rows)
    if choice == "bouvard":
        return [r for r in rows if r["era"] != "ancient"]
    if choice == "adams":
        return [r for r in rows
                if not (r["observer"] == "Flamsteed" and r["epoch"][:4] == "1690")]  # fmt: skip
    raise ValueError(f"unknown data choice {choice!r}; choose from {DATA_CHOICES}")


def synthetic_dataset(seed: int, choice: str = "leverrier", rows=None):
    """(observations, context, (p, q), Uranus box) for a synthetic data set on the layout of
    the chosen rows: epochs and groups from Table 77, values from the true Uranus with drawn
    errors and offsets (the A6c generator). No observed value is read."""
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.model import observation as om
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    rows = select_rows(om.read_rows() if rows is None else rows, choice)
    structure = om.structure_from_rows(rows, om.ErrorModel.load(), allow_empty_groups=True)
    truth = om.true_uranus_ra_dec(structure)
    values, _ = om.synthetic_values(structure, truth,
                                    np.random.default_rng([DATA_RNG_TAG, seed]))  # fmt: skip
    obs = structure.with_values(values)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    box = uranus_box_radec(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    return obs, ctx, pq, box


# ------------------------------------------------------------------ M6d: the shape model


def shape_space(pq, box, cell: int, a_au: float = LE_VERRIER_A_AU):
    """The Le Verrier shape model: ln a held at ``a_au``, everything else as the headline's
    space (inclined perturber, e ≤ 0.3), on one of the two e cells."""
    from fork.inference.priors import PosteriorSpace

    held = PosteriorSpace((None, math.log(a_au), None, None, None), e_max=E_MAX, pq=pq,
                          uranus_box=box, inclined=True)  # fmt: skip
    return held.with_e_slice(*SHAPE_E_CELLS[cell])


def peirce_probability(lon_deg: np.ndarray, target_deg: float = GALLE_LON_DEG,
                       half_width_deg: float = PEIRCE_HALF_WIDTH_DEG) -> float:  # fmt: skip
    """P(|λ − target| ≤ half width) over equal-weight posterior samples of the longitude."""
    d = (np.asarray(lon_deg) - target_deg + 180.0) % 360.0 - 180.0
    return float(np.mean(np.abs(d) <= half_width_deg))


def peirce_ratio(lon_shape_deg: np.ndarray, lon_free_deg: np.ndarray, rng: np.random.Generator,
                 n_boot: int = 2000, **kw) -> dict:  # fmt: skip
    """The M6d ratio P_shape / P_free with a bootstrap 68 % interval over the two sample
    sets (resampled independently). ``ratio`` is None when P_free is 0."""
    a, b = np.asarray(lon_shape_deg), np.asarray(lon_free_deg)
    pa, pb = peirce_probability(a, **kw), peirce_probability(b, **kw)
    boot = []
    for _ in range(n_boot):
        qa = peirce_probability(a[rng.integers(0, a.size, a.size)], **kw)
        qb = peirce_probability(b[rng.integers(0, b.size, b.size)], **kw)
        if qb > 0.0:
            boot.append(qa / qb)
    lo, hi = np.quantile(boot, [0.16, 0.84]).tolist() if boot else (None, None)
    return {"p_shape": pa, "p_free": pb, "ratio": pa / pb if pb > 0.0 else None,
            "ratio_lo68": lo, "ratio_hi68": hi, "n_boot_used": len(boot)}  # fmt: skip


def log_bayes_factor(lnz_a: float, err_a: float, lnz_b: float, err_b: float) -> tuple[float, float]:
    """ln B(a : b) and its error from two evidences with independent errors."""
    return lnz_a - lnz_b, math.hypot(err_a, err_b)


# ------------------------------------------------------------------ S1 cell 0: sub-cells


def subcell_space(pq, box, k: int):
    """The free-a space of S1 (e ≤ 0.5) restricted to sub-cell k of cell 0."""
    from fork.inference.nested import cell_space
    from fork.inference.priors import free_a_space

    return cell_space(free_a_space(pq, S1_E_MAX, box, inclined=True), S1_SUBCELLS[k])


def combine_subcells(lnz: np.ndarray, err: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """ln Z of the union of cells from per-cell ln Z_k (each under its own renormalised
    prior), their errors and the cells' prior masses w_k: Z = Σ w_k Z_k / Σ w_k. The error
    is propagated to first order, cells independent."""
    lnz, err, w = np.asarray(lnz), np.asarray(err), np.asarray(weights)
    terms = np.log(w / w.sum()) + lnz
    top = terms.max()
    total = top + math.log(np.exp(terms - top).sum())
    share = np.exp(terms - total)
    return float(total), float(np.sqrt(((share * err) ** 2).sum()))


# ------------------------------------------------------------------ M6e: predictive places


def predictive_places(space, ctx, eq: np.ndarray, epochs_iso_tdb: dict[str, str]) -> dict:
    """For equal-weight posterior samples ``eq`` (rows of q): the perturber's geocentric RA
    [h] and Dec [°] of date at each epoch, one integration per sample. Earth's place at the
    epochs comes from the cached ephemeris; the perturber's from the samples alone."""
    from fork.export.p1b_adams import EpochGeometry, derived_samples

    geos = {k: EpochGeometry.from_horizons(k, v) for k, v in epochs_iso_tdb.items()}
    d = derived_samples(space, ctx, eq, geos)
    return {k: {"ra_h": d[k]["ra_h"], "dec_deg": d[k]["dec_deg"]} for k in geos}


def region_summary(ra_h: np.ndarray, dec_deg: np.ndarray) -> dict:
    """Medians and 68 % / 95 % intervals per coordinate (RA unwrapped about its circular
    mean), and the sample mean and covariance on the tangent plane [°] for the 2-D region."""
    ra = np.asarray(ra_h) * 15.0
    centre = math.degrees(math.atan2(np.sin(np.radians(ra)).mean(),
                                     np.cos(np.radians(ra)).mean())) % 360.0  # fmt: skip
    dra = (ra - centre + 180.0) % 360.0 - 180.0
    dec = np.asarray(dec_deg)
    x = dra * math.cos(math.radians(float(np.median(dec))))
    q = lambda v: np.quantile(v, [0.025, 0.16, 0.5, 0.84, 0.975]).tolist()  # noqa: E731
    qr, qd = q(dra), q(dec)
    return {"ra_centre_deg": centre,
            "ra_deg": {"lo95": centre + qr[0], "lo68": centre + qr[1], "median": centre + qr[2],
                       "hi68": centre + qr[3], "hi95": centre + qr[4]},
            "dec_deg": dict(zip(("lo95", "lo68", "median", "hi68", "hi95"), qd, strict=True)),
            "plane_mean_deg": [float(x.mean()), float(dec.mean())],
            "plane_cov_deg2": np.cov(np.vstack([x, dec])).tolist(),
            "cos_dec": math.cos(math.radians(float(np.median(dec))))}  # fmt: skip


def region_level(summary: dict, ra_h: np.ndarray, dec_deg: np.ndarray, obs_ra_deg: float,
                 obs_dec_deg: float) -> float:  # fmt: skip
    """The smallest credible level whose region (by Mahalanobis distance from the sample
    mean) contains the observed place: the fraction of samples closer to the mean than it.
    Below 0.95, the place is inside the 95 % region."""
    cov = np.asarray(summary["plane_cov_deg2"])
    inv = np.linalg.inv(cov)
    mean = np.asarray(summary["plane_mean_deg"])
    dra = (np.asarray(ra_h) * 15.0 - summary["ra_centre_deg"] + 180.0) % 360.0 - 180.0
    pts = np.column_stack([dra * summary["cos_dec"], np.asarray(dec_deg)]) - mean
    d2 = np.einsum("ij,jk,ik->i", pts, inv, pts)
    o = np.array([((obs_ra_deg - summary["ra_centre_deg"] + 180.0) % 360.0 - 180.0)
                  * summary["cos_dec"], obs_dec_deg]) - mean  # fmt: skip
    return float(np.mean(d2 < o @ inv @ o))


# ------------------------------------------------------------------ smoke runs (synthetic)


def cost_per_call(obs, ctx, space, rng: np.random.Generator, n: int = 200) -> float:
    """Seconds per likelihood call at prior draws (after one warm-up call)."""
    from fork.inference import nested

    nested._init_worker(space, obs, ctx)
    u = rng.random((n + 1, space.ndim))
    nested._log_likelihood(space.prior_transform(u[0]))
    t0 = time.perf_counter()
    for row in u[1:]:
        nested._log_likelihood(space.prior_transform(row))
    return (time.perf_counter() - t0) / n


def smoke(maxiter: int = 300, nlive: int = 100) -> dict:
    """Short nested runs on synthetic data for each new space and data choice. Not results:
    a few hundred iterations, nothing saved."""
    from fork.inference.nested import run_nested

    out = {}
    cases = [("shape c0, leverrier layout", "leverrier", lambda pq, box: shape_space(pq, box, 0)),
             ("shape c1, leverrier layout", "leverrier", lambda pq, box: shape_space(pq, box, 1)),
             ("S1 sub-cell 0", "leverrier", lambda pq, box: subcell_space(pq, box, 0)),
             ("S1 sub-cell 4", "leverrier", lambda pq, box: subcell_space(pq, box, 4)),
             ("free a cell 2, bouvard layout", "bouvard", None),
             ("free a cell 2, adams layout", "adams", None)]  # fmt: skip
    for label, choice, make in cases:
        obs, ctx, pq, box = synthetic_dataset(1, choice)
        if make is None:
            from fork.inference.nested import cell_space, cells_for
            from fork.inference.priors import free_a_space

            space = cell_space(free_a_space(pq, E_MAX, box, inclined=True), cells_for(E_MAX)[2])
        else:
            space = make(pq, box)
        rng = np.random.default_rng([DATA_RNG_TAG, 99])
        per_call = cost_per_call(obs, ctx, space, rng)
        res = run_nested(obs, ctx, space, rng, nlive=nlive, maxiter=maxiter)
        out[label] = {"ndim": space.ndim, "n_measurements": int(obs.n_measurements)
                      if hasattr(obs, "n_measurements") else None,
                      "ms_per_call": round(1e3 * per_call, 3), "niter": int(res.niter),
                      "ncall": int(res.ncall), "wall_s": round(float(res.wall_s), 1)}  # fmt: skip
        print(label, out[label], flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=["smoke"])
    ap.parse_args()
    smoke()


if __name__ == "__main__":
    main()
