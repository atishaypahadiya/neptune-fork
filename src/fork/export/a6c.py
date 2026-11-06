"""A6c: synthetic coverage with the M6 observation model → ``site/data/a6c_coverage.json``.

    uv run python -m fork.export.a6c plan                   # jobs and their state
    uv run python -m fork.export.a6c run [--n-workers 12 --per-job 4]   # all 40 jobs
    uv run python -m fork.export.a6c job JOB_ID [--n-workers 4]         # one nested run
    uv run python -m fork.export.a6c summarize              # the JSON, with the A6c count

The test (PHASE1_SPEC M6, pre-registered): on 10 synthetic data sets made with the M6
observation model, the truth's 1846 longitude lies inside the posterior 68 % interval in 5–8
of 10. This is the same rule as A3b.

**Data, per noise seed s** (``dataset``). Only Table 77's layout is used: the real 279
epochs, which rows are RA-only, the error groups, and so the same gaps. No observed value is
read. On that layout:
- truth = the real Uranus (Horizons, so including the real 3-D Neptune), as apparent RA/Dec
  of date;
- plus a group offset drawn from N(0, τ_g²), with τ pre-registered at 2″ modern and 5″
  ancient;
- plus Gaussian noise of the calibrated σ_g (``data/historical/m6_error_model.json``),
  rng ``[606, s]``.

**Fit.** M5b's sampler settings and partition: free a, e ≤ 0.3, the four (a, e) cells,
nlive 2000 in [20, 25] AU and 1000 elsewhere, ``ERR_INFLATION`` on ln Z errors. The
perturber is inclined (i ≤ 10° relative to Uranus's plane). Group offsets are marginalised,
and Uranus's box comes from a Uranus-only RA/Dec fit (``fit_radec``).

**Truth** for the coverage count: Neptune's light-time-corrected geocentric ecliptic-of-date
longitude on 1846-09-23, as in A3b. Jobs checkpoint every 5 min, log to
data/cache/nested/a6c/progress.log and resume; ``run`` can be restarted at any time.
"""

from __future__ import annotations

import argparse
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
CACHE = ROOT / "data" / "cache" / "nested" / "a6c"
PROGRESS_LOG = CACHE / "progress.log"
DRIVER_LOG = CACHE / "driver.log"
OUT_NAME = "a6c_coverage.json"
SEEDS = tuple(range(10))
E_MAX = 0.3
NLIVE_DEFAULT = 1000
NLIVE_CELL = {0: 2000}  # cell [20, 25] AU, as in M5b
DATA_RNG_TAG = 606
A6C_RANGE = (0.50, 0.85)  # 5–8 of 10, the A3b rule
DERIVED_MAX_SAMPLES = 6000


@dataclass(frozen=True)
class Job:
    seed: int
    cell: int

    @property
    def id(self) -> str:
        return f"seed_{self.seed}/c{self.cell}"

    @staticmethod
    def from_id(job_id: str) -> Job:
        key, part = job_id.split("/")
        return Job(int(key.removeprefix("seed_")), int(part.removeprefix("c")))

    @property
    def nlive(self) -> int:
        return NLIVE_CELL.get(self.cell, NLIVE_DEFAULT)

    def path(self, suffix: str) -> Path:
        return CACHE / (self.id.replace("/", "__") + suffix)

    def rng(self) -> np.random.Generator:
        return np.random.default_rng([2026, 6, self.seed, self.cell + 1])


def all_jobs() -> list[Job]:
    from fork.inference.nested import cells_for

    return [Job(s, k) for s in SEEDS for k in range(len(cells_for(E_MAX)))]


def dataset(seed: int, rows: list[dict[str, str]] | None = None):
    """(observations, context, Uranus (i, Ω), (p, q), Uranus box, drawn offsets [″])."""
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.model import observation as om
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    rows = om.read_rows() if rows is None else rows
    structure = om.structure_from_rows(rows, om.ErrorModel.load())
    truth = om.true_uranus_ra_dec(structure)
    values, offsets = om.synthetic_values(structure, truth,
                                          np.random.default_rng([DATA_RNG_TAG, seed]))  # fmt: skip
    obs = structure.with_values(values)
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    box = uranus_box_radec(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    return obs, ctx, inc_raan, pq, box, offsets


def job_space(job: Job, pq, box):
    from fork.inference.nested import cell_space, cells_for
    from fork.inference.priors import free_a_space

    space = free_a_space(pq, E_MAX, box, inclined=True)
    return cell_space(space, cells_for(E_MAX)[job.cell])


def do_job(job: Job, n_workers: int, maxiter: int | None = None) -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult, run_nested

    out = job.path(".pkl")
    if out.exists():
        print(f"{job.id}: cached")
        return
    obs, ctx, _, pq, box, _ = dataset(job.seed)
    space = job_space(job, pq, box)
    resumed = job.path(".ckpt").exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(obs, ctx, space, job.rng(), nlive=job.nlive, n_workers=n_workers,
                     maxiter=maxiter, progress_log=PROGRESS_LOG,
                     checkpoint_file=None if maxiter else job.path(".ckpt"),
                     label=f"a6c/{job.id}")  # fmt: skip
    if maxiter is not None:  # smoke runs are not results
        print(f"{job.id}: smoke run, {res.niter} it, {res.ncall} calls, {res.wall_s:.1f} s")
        return
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append(f"a6c/{job.id}", "nested", n_likelihood_calls=res.ncall,
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


def run_all(n_workers: int, per_job: int) -> None:
    """Scheduler: each job a subprocess with ``per_job`` workers, longest (nlive 2000) first,
    one retry per failed job (it resumes from its checkpoint)."""
    slots = max(1, n_workers // per_job)
    pending = sorted([j for j in all_jobs() if not j.path(".pkl").exists()],
                     key=lambda j: (-j.nlive, j.seed, j.cell))  # fmt: skip
    _log(f"run: {len(pending)} jobs pending, {slots} at a time × {per_job} workers")
    running: dict[str, tuple[subprocess.Popen, float]] = {}
    failures: dict[str, int] = {}
    while pending or running:
        for jid, (proc, t0) in list(running.items()):
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
        while pending and len(running) < slots:
            job = pending.pop(0)
            log = job.path(".log")
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a") as fh:
                proc = subprocess.Popen(
                    [sys.executable, "-m", "fork.export.a6c", "job", job.id, "--n-workers",
                     str(per_job)], stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)  # fmt: skip
            running[job.id] = (proc, time.time())
            _log(f"start {job.id} (nlive {job.nlive}, pid {proc.pid})")
        time.sleep(30)
    _log(f"run finished; failures: {failures or 'none'}")


def plan() -> None:
    jobs = all_jobs()
    for j in jobs:
        state = ("done" if j.path(".pkl").exists() else
                 "resumable" if j.path(".ckpt").exists() else "pending")  # fmt: skip
        print(f"{j.id:16s} nlive {j.nlive:5d}  {state}")
    print(f"{len(jobs)} jobs, {sum(j.path('.pkl').exists() for j in jobs)} done")


# ------------------------------------------------------------------ summary


def truth_params() -> dict:
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export.m5b import truth_params as m5b_truth
    from fork.model import bodies
    from fork.model.forward import inclined_perturber_theta_from_state, uranus_inc_raan_at_epoch

    th = inclined_perturber_theta_from_state(about("neptune", [A1_EPOCH_JD_TDB])[0],
                                             bodies.NEPTUNE.mass_msun,
                                             uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB))  # fmt: skip
    return m5b_truth() | {"inc_rel_uranus_deg": round(math.degrees(th[5]), 4)}


def summarize_seed(seed: int, truth: dict) -> dict | None:
    from fork.export.m5b import _e_summary, _lon_summary, _summary_1d
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult, cells_for
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad

    jobs = [Job(seed, k) for k in range(len(cells_for(E_MAX)))]
    if not all(j.path(".pkl").exists() for j in jobs):
        return None
    res = [NestedResult.load(j.path(".pkl")) for j in jobs]
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng([7, 606, seed])
    eq = p.equal_weight(rng)
    _, ctx, _, _, _, offsets = dataset(seed)
    geo = GalleGeometry.from_horizons()
    space = res[0].space
    idx = rng.choice(len(eq), size=min(DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    lon = np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo))
                    for v in eq[idx]])  # fmt: skip
    thetas = np.array([space.theta(v) for v in eq])
    cells = cells_for(E_MAX)
    return {
        "seed": seed, "n_samples": len(eq), "n_derived": len(lon),
        "drawn_offsets_arcsec": [round(float(v), 3) for v in offsets],
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
        "galle_lon_ecl_of_date_deg": _lon_summary(lon, truth["galle_lon_ecl_of_date_deg"]),
        "mass_earth": _summary_1d(10.0 ** thetas[:, 6], truth["mass_earth"], 3),
        "a_au": _summary_1d(thetas[:, 7], truth["a_au"], 3),
        "e": _e_summary(thetas[:, 8], truth["e"], E_MAX),
        "inc_rel_uranus_deg": _summary_1d(np.degrees(thetas[:, 11]), truth["inc_rel_uranus_deg"],
                                          3),
    }  # fmt: skip


def coverage(seeds: list[dict]) -> dict:
    from fork.export.m5b import binomial_calibration, offset_sigma

    lons = [d["galle_lon_ecl_of_date_deg"] for d in seeds]
    hits = [x["truth_in_68"] for x in lons]
    per_seed = []
    for d, x in zip(seeds, lons, strict=True):
        off, z = offset_sigma(x)
        per_seed.append({"seed": d["seed"], "truth_in_68": x["truth_in_68"],
                         "median_deg": x["median"], "width68_deg": x["width68"],
                         "median_minus_truth_deg": round(off, 3),
                         "offset_sigma": None if z is None else round(z, 3)})  # fmt: skip
    k, n = int(sum(hits)), len(hits)
    zs = [r["offset_sigma"] for r in per_seed if r["offset_sigma"] is not None]
    return {
        "n_seeds": n,
        "n_truth_in_68": k,
        "pass_5_to_8_of_10": bool(n == 10 and 5 <= k <= 8),
        "binomial": binomial_calibration(k, n) if n else None,
        "per_seed": per_seed,
        "offset_sigma_rms": round(float(np.sqrt(np.mean(np.square(zs)))), 3) if zs else None,
    }


def summarize() -> None:
    from fork import ledger
    from fork.inference.nested import BOUND, DLOGZ, ERR_INFLATION, SAMPLE, cells_for
    from fork.model.observation import GROUPS, ErrorModel

    t0, c0 = time.time(), ledger.cpu_now()
    truth = truth_params()
    seeds = [d for d in (summarize_seed(s, truth) for s in SEEDS) if d is not None]
    em = ErrorModel.load()
    payload = {
        "description": "A6c: partitioned nested-sampling posteriors on synthetic data with the "
        "M6 observation model (real Table 77 epochs and RA-only rows; truth = Horizons Uranus "
        "as apparent RA/Dec of date; group offsets ~ N(0, tau^2); Gaussian noise of the "
        "calibrated sigma_g). Inclined perturber, free a, e <= 0.3. Truth = Neptune's "
        "1846-09-23 geocentric ecliptic-of-date longitude. ln Z errors x err_inflation.",
        "pass_rule": "truth inside the 68 % interval in 5-8 of 10 seeds (pre-registered)",
        "error_model": {"groups": list(GROUPS), "sigma_arcsec": em.sigma_arcsec,
                        "offset_prior_arcsec": em.offset_prior_arcsec, "basis": em.basis},
        "sampler": {"bound": BOUND, "sample": SAMPLE, "dlogz": DLOGZ,
                    "nlive_default": NLIVE_DEFAULT, "nlive_cells": NLIVE_CELL,
                    "err_inflation": ERR_INFLATION,
                    "cells": [{"a_au": list(a), "e": None if e is None else list(e)}
                              for a, e in cells_for(E_MAX)]},
        "priors": {"log10_m_earth": [0.0, 2.0], "a_au_log_uniform": [20.0, 60.0],
                   "e_uniform_max": E_MAX, "inc_rel_uranus_deg_max": 10.0,
                   "inc": "uniform in cos i", "angles": "uniform",
                   "uranus": "box ±50 scaled σ around the Uranus-only RA/Dec fit"},
        "truth": truth, "seeds": seeds,
        "missing_seeds": [s for s in SEEDS if s not in {d["seed"] for d in seeds}],
        "a6c_coverage": coverage(seeds),
    }  # fmt: skip
    write_json(OUT_NAME, payload)
    print(json.dumps(payload["a6c_coverage"], indent=1))
    print(f"wrote {SITE_DATA / OUT_NAME}")
    n = sum(d["n_derived"] for d in seeds)
    ledger.append("a6c/summarize", "summarize", n_likelihood_calls=0, n_integrations=n,
                  sim_years=n * (1846.73 - 1781.2), wall_s=time.time() - t0,
                  cpu_s=ledger.cpu_now() - c0, n_workers=1, estimated=True,
                  note="derived 1846-longitude integrations (epoch to 1846-09-23)")  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    p_run = sub.add_parser("run")
    p_run.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    p_run.add_argument("--per-job", type=int, default=4)
    p_job = sub.add_parser("job")
    p_job.add_argument("job_id")
    p_job.add_argument("--n-workers", type=int, default=4)
    p_job.add_argument("--maxiter", type=int, default=None, help="smoke run: not saved")
    sub.add_parser("summarize")
    args = ap.parse_args(argv)
    if args.cmd == "plan":
        plan()
    elif args.cmd == "run":
        run_all(args.n_workers, args.per_job)
    elif args.cmd == "job":
        do_job(Job.from_id(args.job_id), args.n_workers, args.maxiter)
    else:
        summarize()


if __name__ == "__main__":
    main()
