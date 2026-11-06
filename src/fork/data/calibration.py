"""M6b empirical error calibration: observed − JPL Uranus, per error group.

    uv run python -m fork.data.calibration      # table + data/historical/m6_error_model.json

Le Verrier's observed apparent places (``leverrier_1846.csv``, misprint corrections applied)
are compared with the real Uranus at the same epochs. The real Uranus is the Horizons state,
taken through ``fork.model.apparent``: light time, aberration, IAU 2006/2000A to the true
equator and equinox of date. This uses Uranus's true positions only: no Neptune and no
perturber model, so it is not leakage (spec M6b). It sets two things, and nothing else:

- **σ_g**, each group's error scale: the MAD-based σ (1.4826 · median |r − median r|), with a
  bootstrap 68 % interval. A Student-t (ν = 4) ML scale is reported beside it as a cross-check.
- **τ_g**, the offset prior width: the rms of the robust group means within the era (6 modern
  groups, 2 ancient), floored at 1″.

The measured group means are reported and not used. In every fit the offsets stay free
nuisance parameters with prior mean 0 (``fork.model.observation``). Fixing them to O − JPL
would correct 1846 data with modern truth.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

from fork.model import observation as om

N_BOOT = 4000
STUDENT_NU = 4.0
TAU_FLOOR_ARCSEC = 1.0
MAD_TO_SIGMA = 1.482602218505602  # 1 / Φ⁻¹(3/4)


def real_observations(rows: list[dict[str, str]] | None = None) -> om.RaDecObservations:
    """Table 77's measured values in the observation layout. σ and τ are placeholders (1″):
    calibration uses the residuals only."""
    rows = om.read_rows() if rows is None else rows
    unit = om.ErrorModel(dict.fromkeys(om.GROUPS, 1.0), dict.fromkeys(om.GROUPS, 1.0), "unit")
    st = om.structure_from_rows(rows, unit)
    order = {}
    jd = np.array([float(r["epoch_tt_jd"]) for r in rows])
    if np.any(np.diff(jd) < 0):
        raise ValueError("rows must be in time order")
    for r in rows:
        if r["unit"] != "deg":
            raise ValueError(f"{r['id']}: unit {r['unit']!r}")
        order[r["id"]] = math.radians(float(r["value"]))
    return st.with_values(np.array([order[r["id"]] for r in rows]))


def robust_sigma(x: np.ndarray) -> float:
    return float(MAD_TO_SIGMA * np.median(np.abs(x - np.median(x))))


def student_t_scale(x: np.ndarray, nu: float = STUDENT_NU) -> tuple[float, float]:
    """(location, scale) ML for a Student-t with fixed ν."""
    from scipy.stats import t

    loc, scale = t.fit(x, fix_df=nu)[1:]
    return float(loc), float(scale)


@dataclass(frozen=True)
class GroupCalibration:
    group: str
    n: int
    median: float
    mean: float
    sigma_mad: float
    sigma_mad_lo68: float
    sigma_mad_hi68: float
    sigma_t4: float
    median_abs: float
    worst_id: str
    worst_arcsec: float
    n_beyond_3sigma: int
    n_beyond_5sigma: int


def calibrate_group(name: str, r: np.ndarray, ids: np.ndarray, rng: np.random.Generator
                    ) -> GroupCalibration:  # fmt: skip
    med = float(np.median(r))
    s = robust_sigma(r)
    boot = np.array([robust_sigma(r[rng.integers(0, r.size, r.size)]) for _ in range(N_BOOT)])
    lo, hi = (float(v) for v in np.percentile(boot, [16.0, 84.0]))
    _, st = student_t_scale(r)
    k = int(np.argmax(np.abs(r - med)))
    z = np.abs(r - med) / s
    return GroupCalibration(name, int(r.size), med, float(np.mean(r)), s, lo, hi, st,
                            float(np.median(np.abs(r))), str(ids[k]), float(r[k]),
                            int(np.sum(z > 3)), int(np.sum(z > 5)))  # fmt: skip


def calibrate(rng: np.random.Generator) -> tuple[list[GroupCalibration], om.ErrorModel]:
    obs = real_observations()
    resid = obs.residuals_arcsec(*om.true_uranus_ra_dec(obs))
    ids = obs.obs_id[obs.epoch_index]
    table = [calibrate_group(g, resid[obs.group_index == k], ids[obs.group_index == k], rng)
             for k, g in enumerate(om.GROUPS)]  # fmt: skip
    tau = {}
    for era in ("ancient", "modern"):
        sel = [c for c in table if c.group.startswith("ancient") == (era == "ancient")]
        t = max(TAU_FLOOR_ARCSEC, math.sqrt(np.mean([c.median**2 for c in sel])))
        tau |= {c.group: round(t, 2) for c in sel}
    sigma = {c.group: round(c.sigma_mad, 2) for c in table}
    basis = ("O - JPL calibration (fork.data.calibration): sigma = MAD sigma per group; "
             "tau = rms of robust group means within the era, floor 1 arcsec")  # fmt: skip
    return table, om.ErrorModel(sigma, tau, basis)


def main() -> None:
    table, em = calibrate(np.random.default_rng(20261106))
    hdr = (f"{'group':22s} {'N':>4s} {'median':>7s} {'mean':>7s} {'σ_MAD':>6s} {'68% (boot)':>13s}"
           f" {'σ_t4':>6s} {'med|r|':>6s} {'>3σ':>4s} {'>5σ':>4s}  worst")  # fmt: skip
    print(hdr)
    for c in table:
        print(f"{c.group:22s} {c.n:4d} {c.median:7.2f} {c.mean:7.2f} {c.sigma_mad:6.2f} "
              f"[{c.sigma_mad_lo68:5.2f},{c.sigma_mad_hi68:5.2f}] {c.sigma_t4:6.2f} "
              f"{c.median_abs:6.2f} {c.n_beyond_3sigma:4d} {c.n_beyond_5sigma:4d}  "
              f"{c.worst_arcsec:+.1f}″ {c.worst_id}")  # fmt: skip
    print("tau:", em.offset_prior_arcsec)
    out = {"basis": em.basis, "sigma_arcsec": em.sigma_arcsec,
           "offset_prior_arcsec": em.offset_prior_arcsec,
           "note": "group means are reported only; offsets stay free nuisance parameters "
                   "with prior mean 0 in every fit",
           "table": [c.__dict__ for c in table]}  # fmt: skip
    tmp = om.ERROR_MODEL_FILE.with_suffix(".tmp")
    tmp.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(om.ERROR_MODEL_FILE)
    print(f"wrote {om.ERROR_MODEL_FILE.relative_to(om.ROOT)}")


if __name__ == "__main__":
    main()
