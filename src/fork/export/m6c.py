"""M6c: the real-data fit → ``site/data/m6_fits.json``.

    uv run python -m fork.export.m6c plan                 # jobs and their state
    uv run python -m fork.export.m6c run [--n-workers 12]  # every pending job, in the background
    uv run python -m fork.export.m6c job JOB_ID [--n-workers 4]
    uv run python -m fork.export.m6c summarize

**Pre-registered** in docs/notes/M6.md ("M6c pre-registration"), committed before any sampling.
Nothing below is retuned after seeing a result. ``FROZEN`` pins the inputs by SHA-256: a job
refuses to start if Table 77's CSV or the error model has changed.

**Headline** (spec: the `leverrier` variant on the modern track): all 279 observations of
Table 77 (529 measurements), the M6 observation model (apparent RA/Dec of date, IAU 2006/2000A),
calibrated σ_g, group offsets marginalised with τ = 2″ modern and 5″ ancient, the inclined
perturber, free a (log-uniform 20–60 AU), e ≤ 0.3, M5b's four (a, e) cells. Nested sampling
with A6c's settings: nlive 2000 in [20, 25] AU and 1000 elsewhere, dlogz 0.1, ``ERR_INFLATION``
on ln Z errors.

**No leakage (A6d).** The fit reads Table 77, the error model, and Horizons states for the
Sun–Jupiter–Saturn context and the Earth (apparent places). It never reads Neptune. The
comparison targets in the summary are constants (``JPL_*``), checked against Horizons by a test
outside this module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.export.to_json import A1_EPOCH_JD_TDB, SITE_DATA, write_json

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "cache" / "nested" / "m6c"
PROGRESS_LOG = CACHE / "progress.log"
DRIVER_LOG = CACHE / "driver.log"
OUT_NAME = "m6_fits.json"
# Sensitivity variants S1–S3 are pre-registered in docs/notes/M6.md ("M6c sensitivity runs
# S1–S3: pre-registration"). Each changes one thing against the headline; the order fixes the
# rng seeds, so it never changes.
VARIANTS = ("headline", "s1_e05", "s2_sigma_data", "s3_t4")
E_MAX = 0.3
S2_ERROR_MODEL = "data/historical/m6_error_model_s2_data_only.json"
VARIANT_CONFIG = {
    "headline": {"e_max": E_MAX, "error_model": "data/historical/m6_error_model.json",
                 "likelihood": "gaussian"},
    "s1_e05": {"e_max": 0.5, "error_model": "data/historical/m6_error_model.json",
               "likelihood": "gaussian"},
    "s2_sigma_data": {"e_max": E_MAX, "error_model": S2_ERROR_MODEL, "likelihood": "gaussian"},
    "s3_t4": {"e_max": E_MAX, "error_model": "data/historical/m6_error_model.json",
              "likelihood": "student_t", "nu": 4.0},
}  # fmt: skip
# S1: M5b's four cells with the last one's e range widened to the new bound, so the cells
# still partition the prior (``cells_for`` would clip it at 0.3).
S1_CELLS = (((20.0, 25.0), None), ((25.0, 28.0), None), ((28.0, 60.0), (0.0, 0.12)),
            ((28.0, 60.0), (0.12, 0.5)))  # fmt: skip
NLIVE_DEFAULT = 1000
NLIVE_CELL = {0: 2000}  # cell [20, 25] AU, as in M5b and A6c
WORKERS_DEFAULT = 2
WORKERS_CELL = {0: 6}  # the nlive-2000 cell gets more workers, so all four end together
DERIVED_MAX_SAMPLES = 6000
SUMMARY_RNG = (7, 1846, 6)

# Inputs, frozen at pre-registration (2026-11-07).
FROZEN = {
    "table77_csv": ("data/historical/leverrier_1846.csv",
                    "2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4"),
    "error_model": ("data/historical/m6_error_model.json",
                    "0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277"),
}  # fmt: skip
# Inputs of the sensitivity variants, frozen at their pre-registration (2026-11-07).
FROZEN_VARIANT = {
    "s2_sigma_data": {"s2_error_model": (
        S2_ERROR_MODEL, "79c8ddc69a61ba36ab44f0c790f166e25df8fee9a53f339e7dadc1aec490fa2d")},
}
TAU_ARCSEC = {"ancient": 5.0, "modern": 2.0}

# Comparison targets (reporting only; never enter a fit).
# JPL: Horizons Neptune, light-time-corrected geocentric ecliptic longitude of date on
# 1846-09-23 23:00 TDB (Galle's night), as in M5 (docs/notes/M5.md); and the heliocentric
# (inner-barycentre) ecliptic-of-date longitude at 1847 Jan 1.0, the epoch of Le Verrier's
# elements. tests/test_m6c.py recomputes both from Horizons.
JPL_GALLE_LON_DEG = 325.873
JPL_HELIO_1847_LON_DEG = 327.561
# Le Verrier, Comptes rendus 23, 428 (31 Aug 1846): heliocentric true longitude 326°32′ at
# 1847 Jan 1, from ``fork.model.hidden_mass.LE_VERRIER_1846`` (Krajnović 2016, A&G 57, 5.28,
# Table 2, after Grosser 1962; Encyclopædia Britannica 1911). His geocentric longitude on
# Galle's night is derived (``le_verrier_galle_lon_deg``), not printed in a source we hold.
INC_PRIOR_MAX_DEG = 10.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_frozen(variant: str = "headline") -> None:
    """Refuse to run on inputs other than the pre-registered ones (the headline's, plus the
    variant's own)."""
    from fork.inference.nested import DLOGZ
    from fork.inference.priors import INC_PRIOR_MAX_DEG as INC
    from fork.model.observation import ErrorModel

    for name, (rel, digest) in (FROZEN | FROZEN_VARIANT.get(variant, {})).items():
        if sha256(ROOT / rel) != digest:
            raise RuntimeError(f"{name} ({rel}) differs from the pre-registered file")
    em = ErrorModel.load(ROOT / VARIANT_CONFIG[variant]["error_model"])
    for g, t in em.offset_prior_arcsec.items():
        if abs(t - TAU_ARCSEC["ancient" if g.startswith("ancient") else "modern"]) > 1e-12:
            raise RuntimeError(f"tau for {g} is {t}, not the pre-registered value")
    if abs(DLOGZ - 0.1) > 1e-12 or abs(INC - INC_PRIOR_MAX_DEG) > 1e-12:
        raise RuntimeError("sampler or prior settings differ from the pre-registration")


@dataclass(frozen=True)
class Job:
    variant: str
    cell: int

    @property
    def id(self) -> str:
        return f"{self.variant}/c{self.cell}"

    @staticmethod
    def from_id(job_id: str) -> Job:
        variant, part = job_id.split("/")
        return Job(variant, int(part.removeprefix("c")))

    @property
    def nlive(self) -> int:
        return NLIVE_CELL.get(self.cell, NLIVE_DEFAULT)

    @property
    def workers(self) -> int:
        return WORKERS_CELL.get(self.cell, WORKERS_DEFAULT)

    def path(self, suffix: str) -> Path:
        return CACHE / (self.id.replace("/", "__") + suffix)

    def rng_seed(self) -> list[int]:
        return [2026, 6, 1846, VARIANTS.index(self.variant), self.cell + 1]

    def rng(self) -> np.random.Generator:
        return np.random.default_rng(self.rng_seed())


def variant_cells(variant: str) -> tuple:
    from fork.inference.nested import cells_for

    return S1_CELLS if variant == "s1_e05" else cells_for(VARIANT_CONFIG[variant]["e_max"])


def all_jobs(variants: tuple[str, ...] = VARIANTS) -> list[Job]:
    return [Job(v, k) for v in variants for k in range(len(variant_cells(v)))]


def t_scales_arcsec() -> dict[str, float]:
    """S3's Student-t scales: the ν = 4 ML scales of O − JPL per group, recorded in the frozen
    headline error model's table by the M6b calibration (``sigma_t4``)."""
    d = json.loads((ROOT / FROZEN["error_model"][0]).read_text(encoding="utf-8"))
    return {row["group"]: float(row["sigma_t4"]) for row in d["table"]}


def dataset(rows: list[dict[str, str]] | None = None, variant: str = "headline"):
    """(observations, context, Uranus (i, Ω), (p, q), Uranus box): Table 77's measured values
    with the variant's frozen error model. The box follows the M5 rule on the variant's
    Gaussian error model (S3's is the headline's). Reads no Neptune state."""
    from fork.data.calibration import real_observations
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.model import observation as om
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    cfg = VARIANT_CONFIG[variant]
    em = om.ErrorModel.load(ROOT / cfg["error_model"])
    obs = real_observations(om.read_rows() if rows is None else rows, em)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    box = uranus_box_radec(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    if cfg["likelihood"] == "student_t":
        from fork.model.student_t import StudentTObservations

        obs = StudentTObservations.from_gaussian(obs, cfg["nu"], t_scales_arcsec())
    return obs, ctx, inc_raan, pq, box


def job_space(job: Job, pq, box):
    from fork.inference.nested import cell_space
    from fork.inference.priors import free_a_space

    space = free_a_space(pq, VARIANT_CONFIG[job.variant]["e_max"], box, inclined=True)
    return cell_space(space, variant_cells(job.variant)[job.cell])


def do_job(job: Job, n_workers: int, maxiter: int | None = None) -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult, run_nested

    check_frozen(job.variant)
    out = job.path(".pkl")
    if out.exists():
        print(f"{job.id}: cached")
        return
    obs, ctx, _, pq, box = dataset(variant=job.variant)
    space = job_space(job, pq, box)
    resumed = job.path(".ckpt").exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(obs, ctx, space, job.rng(), nlive=job.nlive, n_workers=n_workers,
                     maxiter=maxiter, progress_log=PROGRESS_LOG,
                     checkpoint_file=None if maxiter else job.path(".ckpt"),
                     label=f"m6c/{job.id}")  # fmt: skip
    if maxiter is not None:  # smoke runs are not results
        print(f"{job.id}: smoke run, {res.niter} it, {res.ncall} calls, {res.wall_s:.1f} s")
        return
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append(f"m6c/{job.id}", "nested", n_likelihood_calls=res.ncall,
                  n_integrations=res.ncall, sim_years=res.ncall * yrs,
                  wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0, n_workers=n_workers,
                  estimated=resumed, resumed=resumed, nlive=job.nlive,
                  niter=res.niter)  # fmt: skip
    print(f"{job.id}: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, "
          f"{res.wall_s / 60:.1f} min")  # fmt: skip


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with DRIVER_LOG.open("a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def run_all(n_workers: int, variants: tuple[str, ...] = VARIANTS) -> None:
    """Scheduler: each job a subprocess with ``job.workers`` workers, within a budget of
    ``n_workers``; longest first; one retry per failed job (it resumes from its checkpoint)."""
    for v in variants:
        check_frozen(v)
    pending = sorted([j for j in all_jobs(variants) if not j.path(".pkl").exists()],
                     key=lambda j: (-j.nlive, j.variant, j.cell))  # fmt: skip
    _log(f"run {','.join(variants)}: {len(pending)} jobs pending, budget {n_workers} workers")
    running: dict[str, tuple[subprocess.Popen, float, int]] = {}
    failures: dict[str, int] = {}
    while pending or running:
        for jid, (proc, t0, _) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[jid]
            job = Job.from_id(jid)
            if rc == 0 and job.path(".pkl").exists():
                _log(f"done {jid} ({(time.time() - t0) / 60:.1f} min)")
            else:
                failures[jid] = failures.get(jid, 0) + 1
                _log(f"FAILED {jid} (exit {rc}, attempt {failures[jid]})")
                if failures[jid] < 2:
                    pending.insert(0, job)
        busy = sum(w for *_, w in running.values())
        while pending and (not running or busy + pending[0].workers <= n_workers):
            job = pending.pop(0)
            w = min(job.workers, n_workers)
            log = job.path(".log")
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a") as fh:
                proc = subprocess.Popen(
                    [sys.executable, "-m", "fork.export.m6c", "job", job.id, "--n-workers",
                     str(w)], stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)  # fmt: skip
            running[job.id] = (proc, time.time(), w)
            busy += w
            _log(f"start {job.id} (nlive {job.nlive}, {w} workers, pid {proc.pid})")
        time.sleep(30)
    _log(f"run finished; failures: {failures or 'none'}")


def plan() -> None:
    jobs = all_jobs()
    for j in jobs:
        state = ("done" if j.path(".pkl").exists() else
                 "resumable" if j.path(".ckpt").exists() else "pending")  # fmt: skip
        print(f"{j.id:16s} nlive {j.nlive:5d}  workers {j.workers:2d}  {state}")
    print(f"{len(jobs)} jobs, {sum(j.path('.pkl').exists() for j in jobs)} done")


# ------------------------------------------------------------------ comparison targets


def le_verrier_state_ecl_j2000(jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Le Verrier's published orbit (``LE_VERRIER_1846``) as a heliocentric Kepler ellipse in
    the ecliptic and mean equinox of 1847.0 (assumed, as in M5), rotated to ecl_j2000."""
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.model.bodies import M_CENTRAL_MSUN, G
    from fork.model.hidden_mass import LE_VERRIER_1846 as LV

    mu = G * (M_CENTRAL_MSUN + LV.mass_msun)
    a, e = LV.a_au, LV.e
    varpi = math.radians(LV.varpi_ecl_of_date_deg)
    mean_anom = math.radians(LV.mean_lon_ecl_of_date_deg) - varpi + math.sqrt(mu / a**3) * (
        jd_tdb - LV.epoch_jd)  # fmt: skip
    ecc = mean_anom
    for _ in range(50):
        ecc -= (ecc - e * math.sin(ecc) - mean_anom) / (1.0 - e * math.cos(ecc))
    half = ecc / 2.0
    nu = 2.0 * math.atan2(math.sqrt(1 + e) * math.sin(half), math.sqrt(1 - e) * math.cos(half))
    r, lon, p = a * (1 - e * math.cos(ecc)), nu + varpi, a * (1 - e * e)
    vr, vt = math.sqrt(mu / p) * e * math.sin(nu), math.sqrt(mu * p) / r
    pos = np.array([r * math.cos(lon), r * math.sin(lon), 0.0])
    vel = np.array([vr * math.cos(lon) - vt * math.sin(lon),
                    vr * math.sin(lon) + vt * math.cos(lon), 0.0])  # fmt: skip
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(LV.epoch_jd)
    return rot.T @ pos, rot.T @ vel


def le_verrier_galle_lon_deg() -> float:
    """Geocentric ecliptic-of-date longitude of Le Verrier's published planet on Galle's night,
    by the same light-time computation as the posterior's (derived, not printed)."""
    from fork.inference.profile import GalleGeometry

    geo = GalleGeometry.from_horizons()
    pos, vel = le_verrier_state_ecl_j2000(geo.jd_tdb)
    return math.degrees(geo.lon_ecl_of_date_rad(pos, vel)) % 360.0


def comparison_targets() -> dict:
    from fork.model.hidden_mass import LE_VERRIER_1846 as LV

    return {
        "jpl": {"galle_lon_ecl_of_date_deg": JPL_GALLE_LON_DEG,
                "helio_lon_ecl_of_date_1847_deg": JPL_HELIO_1847_LON_DEG,
                "source": "JPL Horizons (constants; tests/test_m6c.py recomputes them)"},
        "le_verrier": {"helio_lon_ecl_of_date_1847_deg": round(LV.true_lon_ecl_of_date_deg, 4),
                       "galle_lon_ecl_of_date_deg_derived": round(le_verrier_galle_lon_deg(), 4),
                       "sources": list(LV.sources),
                       "note": "Galle-night longitude derived from the published elements "
                               "(Kepler orbit, equinox 1847.0 assumed); not printed in a "
                               "source in the repo"},
    }  # fmt: skip


# ------------------------------------------------------------------ summary


def _quantiles(x: np.ndarray, nd: int) -> dict:
    q = np.percentile(x, [2.5, 16.0, 50.0, 84.0, 97.5])
    return {"median": round(float(q[2]), nd), "lo68": round(float(q[1]), nd),
            "hi68": round(float(q[3]), nd), "lo95": round(float(q[0]), nd),
            "hi95": round(float(q[4]), nd)}  # fmt: skip


def lon_summary(lon_deg: np.ndarray, targets: dict[str, float]) -> dict:
    """Circular median, 68 % and 95 % intervals; for each target, median − target in degrees
    and in 68 % half-widths (on the side facing the target), and whether each interval holds it."""
    from fork.inference.mcmc import circular_interval_deg

    lo95, lo, med, hi, hi95 = circular_interval_deg(lon_deg, (2.5, 16.0, 50.0, 84.0, 97.5))

    def rel(v: float) -> float:
        return (v - med + 180.0) % 360.0 - 180.0

    out = {"median": round(med, 3), "lo68": round(lo, 3), "hi68": round(hi, 3),
           "lo95": round(lo95, 3), "hi95": round(hi95, 3),
           "width68": round(rel(hi) - rel(lo), 3), "width95": round(rel(hi95) - rel(lo95), 3),
           "vs": {}}  # fmt: skip
    for name, t in targets.items():
        off = -rel(t)
        half = -rel(lo) if off > 0 else rel(hi)
        out["vs"][name] = {"target": round(t, 4), "median_minus_target_deg": round(off, 3),
                           "offset_sigma": round(off / half, 3) if half > 0 else None,
                           "in_68": bool(rel(lo) <= rel(t) <= rel(hi)),
                           "in_95": bool(rel(lo95) <= rel(t) <= rel(hi95))}  # fmt: skip
    return out


def inc_prior_quantiles_deg(q=(16.0, 50.0, 84.0)) -> list[float]:
    """i quantiles under the prior: cos i uniform on [cos 10°, 1]."""
    c = math.cos(math.radians(INC_PRIOR_MAX_DEG))
    return [round(math.degrees(math.acos(1.0 - p / 100.0 * (1.0 - c))), 3) for p in q]


def summarize_variant(variant: str) -> dict | None:
    from fork.export.m5b import _e_summary, _heliocentric_lon_of_date_deg
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
    from fork.model.hidden_mass import JD_1847_JAN1_TDB
    from fork.model.observation import GROUPS

    jobs = [Job(variant, k) for k in range(len(variant_cells(variant)))]
    if not all(j.path(".pkl").exists() for j in jobs):
        return None
    res = [NestedResult.load(j.path(".pkl")) for j in jobs]
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng([*SUMMARY_RNG, VARIANTS.index(variant)])
    eq = p.equal_weight(rng)
    obs, ctx, *_ = dataset(variant=variant)
    cfg = VARIANT_CONFIG[variant]
    geo = GalleGeometry.from_horizons()
    space = res[0].space
    idx = rng.choice(len(eq), size=min(DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    galle = np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo))
                      for v in eq[idx]])  # fmt: skip
    helio = np.array([_heliocentric_lon_of_date_deg(space.theta(v), ctx, JD_1847_JAN1_TDB)
                      for v in eq[idx]])  # fmt: skip
    thetas = np.array([space.theta(v) for v in eq])
    tg = comparison_targets()
    best = max(res, key=lambda r: r.max_logl)
    th_best = space.theta(best.samples[int(np.argmax(best.results["logl"]))])
    resid = obs.residuals_arcsec(*obs.predict(th_best, ctx))

    cells = variant_cells(variant)
    out = {
        "variant": variant, "n_samples": len(eq), "n_derived": len(galle),
        "evidence": {"logz": round(p.logz, 4), "logzerr": round(p.logzerr, 4),
                     "logzerr_inflation": ERR_INFLATION,
                     "cells": [{"a_au": list(c[0]), "e": None if c[1] is None else list(c[1]),
                                "nlive": r.nlive, "share": round(float(s), 5),
                                "share_err": round(float(se), 5), "logz_k": round(r.logz, 4),
                                "niter": r.niter, "ncall": r.ncall,
                                "wall_min": round(r.wall_s / 60, 1),
                                "chi2_min": round(-2 * r.max_logl, 3)}
                               for c, r, s, se in zip(cells, res, p.shares, p.share_errs,
                                                      strict=True)]},
        "galle_lon_ecl_of_date_deg": lon_summary(galle, {
            "jpl": tg["jpl"]["galle_lon_ecl_of_date_deg"],
            "le_verrier_derived": tg["le_verrier"]["galle_lon_ecl_of_date_deg_derived"]}),
        "helio_lon_ecl_of_date_1847_deg": lon_summary(helio, {
            "jpl": tg["jpl"]["helio_lon_ecl_of_date_1847_deg"],
            "le_verrier": tg["le_verrier"]["helio_lon_ecl_of_date_1847_deg"]}),
        "peirce_p_within_0p87_deg_of_jpl": round(float(np.mean(
            np.abs((galle - JPL_GALLE_LON_DEG + 180.0) % 360.0 - 180.0) <= 0.87)), 4),
        "mass_earth": _quantiles(10.0 ** thetas[:, 6], 3),
        "a_au": _quantiles(thetas[:, 7], 3),
        "e": _e_summary(thetas[:, 8], None, cfg["e_max"]),
        "inc_rel_uranus_deg": _quantiles(np.degrees(thetas[:, 11]), 3)
        | {"prior_16_50_84": inc_prior_quantiles_deg(),
           "note": "prior-dominated: compare with the prior quantiles"},
        "best_fit": {"chi2_augmented": round(-2 * best.max_logl, 3),
                     "n_measurements": obs.n_measurements,
                     "group_offsets_arcsec": dict(zip(GROUPS, (round(float(v), 3) for v in
                                                  obs.best_offsets_arcsec(resid)), strict=True))},
    }  # fmt: skip
    if variant == "headline":  # its record stays exactly as first published
        return out
    out["config"] = {k: v for k, v in cfg.items()} | {
        "cells": [{"a_au": list(a), "e": None if e is None else list(e)} for a, e in cells],
        "sigma_arcsec": dict(zip(GROUPS, (float(v) for v in obs.sigma_arcsec), strict=True)),
        "minus2_lnl_normalised": cfg["likelihood"] == "student_t"}
    if cfg["likelihood"] == "student_t":
        out["config"]["t_scale_arcsec"] = t_scales_arcsec()
    out["p_e_above_0p3"] = round(float(np.mean(thetas[:, 8] > 0.3)), 4)
    out["p_a_below_25_au"] = round(float(np.mean(thetas[:, 7] < 25.0)), 4)
    return out


def gaussian_norm_m2lnl_for(variant: str) -> float:
    """The constant the Gaussian likelihood drops, for the variant's error model."""
    from fork.data.calibration import real_observations
    from fork.model import observation as om
    from fork.model.student_t import gaussian_norm_m2lnl

    em = om.ErrorModel.load(ROOT / VARIANT_CONFIG[variant]["error_model"])
    return gaussian_norm_m2lnl(real_observations(om.read_rows(), em))


def vs_headline(variant: str, fit: dict, head: dict) -> dict:
    """Shift of the Galle-night median against the headline, and ln B against it (errors ×
    ERR_INFLATION, added in quadrature). ln Z of a Gaussian fit is normalised by adding back
    its dropped constant whenever the two likelihoods differ."""
    g, h = fit["galle_lon_ecl_of_date_deg"], head["galle_lon_ecl_of_date_deg"]
    shift = (g["median"] - h["median"] + 180.0) % 360.0 - 180.0
    zv, zh = fit["evidence"]["logz"], head["evidence"]["logz"]
    err = math.hypot(fit["evidence"]["logzerr"], head["evidence"]["logzerr"])
    lik = VARIANT_CONFIG[variant]["likelihood"]
    note = "same likelihood; the priors differ (e bound)"
    if variant == "s2_sigma_data":
        zv -= 0.5 * gaussian_norm_m2lnl_for(variant)
        zh -= 0.5 * gaussian_norm_m2lnl_for("headline")
        note = ("sigma_g estimated from these data (and the Uranus box rescaled with them): ln B "
                "favours S2 by construction; not a model comparison")
    elif lik == "student_t":
        zh -= 0.5 * gaussian_norm_m2lnl_for("headline")
        note = "normalised Student-t (nu 4) against normalised Gaussian; same data, prior, box"
    return {"galle_median_shift_deg": round(shift, 3),
            "galle_median_shift_headline_half_widths": round(
                shift / (0.5 * h["width68"]), 3),
            "width68_ratio": round(g["width68"] / h["width68"], 3),
            "ln_b": round(zv - zh, 3), "ln_b_err": round(err, 3), "note": note}  # fmt: skip


def summarize() -> None:
    from fork import ledger
    from fork.inference.nested import BOUND, DLOGZ, ERR_INFLATION, SAMPLE, cells_for
    from fork.model.observation import GROUPS, ErrorModel

    check_frozen()
    t0, c0 = time.time(), ledger.cpu_now()
    fits = {}
    for v in VARIANTS:
        if all(j.path(".pkl").exists() for j in all_jobs((v,))):
            check_frozen(v)
            fits[v] = summarize_variant(v)
    if "headline" in fits:
        for v, d in fits.items():
            if v != "headline":
                d["vs_headline"] = vs_headline(v, d, fits["headline"])
    em = ErrorModel.load()
    payload = {
        "description": "M6c: partitioned nested-sampling posteriors on Le Verrier's Table 77 "
        "(279 observations, 529 measurements), M6 observation model (apparent RA/Dec of date), "
        "calibrated sigma_g, group offsets marginalised (tau 2\" modern, 5\" ancient), inclined "
        "perturber, free a, e <= 0.3. Pre-registered in docs/notes/M6.md. ln Z errors x "
        "err_inflation.",
        "frozen_inputs": {k: {"path": p, "sha256": h} for k, (p, h) in FROZEN.items()},
        "error_model": {"groups": list(GROUPS), "sigma_arcsec": em.sigma_arcsec,
                        "offset_prior_arcsec": em.offset_prior_arcsec, "basis": em.basis},
        "sampler": {"bound": BOUND, "sample": SAMPLE, "dlogz": DLOGZ,
                    "nlive_default": NLIVE_DEFAULT, "nlive_cells": NLIVE_CELL,
                    "err_inflation": ERR_INFLATION,
                    "rng": {j.id: j.rng_seed() for j in all_jobs()},
                    "cells": [{"a_au": list(a), "e": None if e is None else list(e)}
                              for a, e in cells_for(E_MAX)]},
        "priors": {"log10_m_earth": [0.0, 2.0], "a_au_log_uniform": [20.0, 60.0],
                   "e_uniform_max": E_MAX, "inc_rel_uranus_deg_max": INC_PRIOR_MAX_DEG,
                   "inc": "uniform in cos i", "angles": "uniform",
                   "uranus": "box ±50 scaled σ around the Uranus-only RA/Dec fit"},
        "comparison_targets": comparison_targets(),
        "fits": fits,
        "missing_variants": [v for v in VARIANTS if v not in fits],
    }  # fmt: skip
    write_json(OUT_NAME, payload)
    print(json.dumps({v: d["galle_lon_ecl_of_date_deg"] for v, d in fits.items()}, indent=1))
    print(f"wrote {SITE_DATA / OUT_NAME}")
    n = 2 * sum(d["n_derived"] for d in fits.values())
    ledger.append("m6c/summarize", "summarize", n_likelihood_calls=0, n_integrations=n,
                  sim_years=n * (1847.0 - 1781.2), wall_s=time.time() - t0,
                  cpu_s=ledger.cpu_now() - c0, n_workers=1, estimated=True,
                  note="derived longitudes: Galle's night and 1847.0, per sample")  # fmt: skip


def make_s2_error_model() -> None:
    """S2's data-only σ_g (pre-registered procedure; run once, before any S2 sampling):
    the highest-likelihood sample of the headline run (all four cells), its Table 77 residuals,
    and per group the MAD σ (``robust_sigma``: 1.4826 · median |r − median r|, so a group
    offset drops out), rounded to 0.01″ as the headline's are. τ is unchanged. No JPL Uranus
    enters σ_g; the point itself was found under the headline's (JPL-calibrated) weights."""
    from fork.data.calibration import robust_sigma
    from fork.inference.nested import NestedResult
    from fork.model.observation import GROUPS, ErrorModel

    check_frozen()
    obs, ctx, *_ = dataset()
    res = [NestedResult.load(Job("headline", k).path(".pkl"))
           for k in range(len(variant_cells("headline")))]  # fmt: skip
    k = max(range(len(res)), key=lambda i: res[i].max_logl)
    i = int(np.argmax(res[k].results["logl"]))
    theta = res[k].space.theta(res[k].samples[i])
    resid = obs.residuals_arcsec(*obs.predict(theta, ctx))
    raw = {g: robust_sigma(resid[obs.group_index == j]) for j, g in enumerate(GROUPS)}
    head = ErrorModel.load()
    payload = {
        "basis": "sigma: MAD sigma per group of the Table 77 residuals at the headline M6c "
                 "run's highest-likelihood sample (no JPL Uranus); tau: pre-registered 2 arcsec "
                 "modern, 5 arcsec ancient",
        "sigma_arcsec": {g: round(v, 2) for g, v in raw.items()},
        "offset_prior_arcsec": head.offset_prior_arcsec,
        "sigma_unrounded_arcsec": raw,
        "source_point": {"run": "m6c headline", "cell": k,
                         "minus2_lnl": round(-2.0 * res[k].max_logl, 4),
                         "theta": [float(v) for v in theta]},
        "note": "M6c sensitivity S2 (docs/notes/M6.md). One iteration, no refit of sigma.",
    }  # fmt: skip
    out = ROOT / S2_ERROR_MODEL
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    tmp.replace(out)
    print(json.dumps(payload["sigma_arcsec"], indent=1))
    print(f"wrote {out}\nsha256 {sha256(out)}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    p_run = sub.add_parser("run")
    p_run.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    p_run.add_argument("--variant", action="append", choices=VARIANTS,
                       help="repeatable; variants run one after another, in the order given "
                            "(default: all, together)")  # fmt: skip
    p_job = sub.add_parser("job")
    p_job.add_argument("job_id")
    p_job.add_argument("--n-workers", type=int, default=4)
    p_job.add_argument("--maxiter", type=int, default=None, help="smoke run: not saved")
    sub.add_parser("summarize")
    sub.add_parser("make-s2-sigma", help="write S2's data-only error model (before sampling)")
    args = ap.parse_args(argv)
    if args.cmd == "plan":
        plan()
    elif args.cmd == "run":
        if args.variant:
            for v in args.variant:
                run_all(args.n_workers, (v,))
        else:
            run_all(args.n_workers)
    elif args.cmd == "job":
        do_job(Job.from_id(args.job_id), args.n_workers, args.maxiter)
    elif args.cmd == "make-s2-sigma":
        make_s2_error_model()
    else:
        summarize()


if __name__ == "__main__":
    main()
