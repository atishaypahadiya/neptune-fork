"""M5b: partitioned dynesty posteriors on synthetic data → ``site/data/m5b_posteriors.json``.

    uv run python -m fork.export.m5b plan                  # jobs and their state
    uv run python -m fork.export.m5b run [--n-workers 12 --per-job 4]   # all jobs (a day)
    uv run python -m fork.export.m5b job JOB_ID [--n-workers 4]         # one nested run
    uv run python -m fork.export.m5b summarize             # the JSON, with A3b

Posteriors (PHASE1_SPEC M5b, as amended after the partitioned validation):
  - ``free_a_e0.3``: free a, e ∈ [0, 0.3], noise seeds 0–9, partitioned into the four cells
    of ``nested.PARTITION_CELLS``, one rng seed per cell. Cell [20, 25] runs at nlive 2000
    (the thin mode-B pocket was under-sampled at 1000), the rest at 1000.
  - ``free_a_e0.1``: free a, e ∈ [0, 0.1], seeds 0–4: the three cells that can hold e ≤ 0.1
    (``nested.cells_for``), nlive 1000.
  - ``le_verrier_shape``: a, e, ϖ held at Le Verrier's published values; m, λ₀ and Uranus
    sampled. Seeds 0–4, one unpartitioned run each, provided seed 0 passes the single-mode
    check (``lv_single_mode``: two rng seeds agree, no mass far from the main peak); otherwise
    seeds 1–4 are not run and the summary says so.
  - noise-free truth, free a, e ≤ 0.3, the same four cells: scheduled after every job A3b needs
    and summarized after A3b, as a prior-bias check (``noise_free_check``: truth inside the
    68 % intervals, and the median's offset from truth in σ).
Every ln Z error that is reported or compared is ``nested.ERR_INFLATION`` × dynesty's quoted
error (the scatter measured on seed 0). The seed-0 e ≤ 0.3 cells [25, 28] and [28, 60] are
the validation runs (same prior and settings), reused.

``run`` is a scheduler: it runs each job as a subprocess with ``--per-job`` workers, longest
first, ``--n-workers`` in total. Jobs log to data/cache/nested/m5b/progress.log, checkpoint
every 5 min, and resume from the checkpoint; finished jobs are cached as .pkl and skipped, so
``run`` can be restarted at any time. Its own log is data/cache/nested/m5b/driver.log.
``summarize`` reports per posterior: ln Z and the cell shares (inflated errors), median and
16/84 % of m, a, e and the 1846-09-23 geocentric longitude (light-time corrected, ecliptic of
date), whether truth lies inside each 68 % interval, the best χ² (posterior maximum or the
bounded multistart point fit, if lower), and the M5a profile's Δχ² ≤ 1 interval beside the
posterior's. A3b: truth inside the 68 % longitude interval in 50–85 % of the noise seeds.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.export.to_json import A1_EPOCH_JD_TDB, SITE_DATA, write_json

ROOT = Path(__file__).resolve().parents[3]
NESTED_DIR = ROOT / "data" / "cache" / "nested"
CACHE = NESTED_DIR / "m5b"
PROGRESS_LOG = CACHE / "progress.log"
DRIVER_LOG = CACHE / "driver.log"
OUT_NAME = "m5b_posteriors.json"
NOISE_SEEDS = tuple(range(10))
# noise-free truth (seed None): free a, e ≤ 0.3, the same four cells; scheduled after every job
# A3b needs, summarized after A3b, as a prior-bias check
NOISE_FREE_RUNS = ("free_a_e0.3",)
RUNS = ("free_a_e0.3", "free_a_e0.1", "le_verrier_shape")
RUN_SEEDS: dict[str, tuple[int, ...]] = {
    "free_a_e0.3": NOISE_SEEDS,
    "free_a_e0.1": (0, 1, 2, 3, 4),
    "le_verrier_shape": (0, 1, 2, 3, 4),
}
E_MAX = {"free_a_e0.3": 0.3, "free_a_e0.1": 0.1}
NLIVE_DEFAULT = 1000
NLIVE_CELL = {("free_a_e0.3", 0): 2000}  # cell [20, 25] under e ≤ 0.3
LV = "le_verrier_shape"
LV_CHECK_SEED = 0  # Le Verrier shape: seed 0 runs twice (rep 1, 2) for the single-mode check
LV_CHECK_SIGMA = 8.0  # "far" = beyond median ± 8 × the 68 % half-width
LV_CHECK_MAX_FAR_MASS = 0.005
LV_CHECK_MEDIAN_TOL = 0.5  # rep medians agree within 0.5 of the 68 % half-width
# seed-0 e ≤ 0.3 validation runs (rng 1), reused for cells 1–3 (same prior, nlive, settings)
REUSE = {
    1: "seed_0_free_a_e0.3_slice1_rng1.pkl",
    2: "seed_0_free_a_e0.3_cell2_a28-60_e0-0.12_rng1.pkl",
    3: "seed_0_free_a_e0.3_cell3_a28-60_e0.12-0.3_rng1.pkl",
}
DERIVED_MAX_SAMPLES = 6000  # 1846 longitude is computed on at most this many posterior draws
CELL_LON_SAMPLES = 1000
A3B_RANGE = (0.50, 0.85)
# summarize's own compute, for its ledger row: derived-quantity integrations, point-fit nfev
_SUMMARY_COUNT = {"derived_integrations": 0, "point_fit_nfev": 0}
# rough cost of one job in "nlive-1000 cell" units, for longest-first ordering and the ETA
LV_COST = 0.7


# ------------------------------------------------------------------ jobs


@dataclass(frozen=True)
class Job:
    seed: int | None  # None: noise-free truth
    run: str
    cell: int | None  # index into nested.cells_for(E_MAX[run]); None: unpartitioned
    rep: int = 1  # Le Verrier shape on seed 0 also runs rep 2 (single-mode check)

    @property
    def id(self) -> str:
        part = "all" if self.cell is None else f"c{self.cell}"
        rep = "" if self.rep == 1 else f"_r{self.rep}"
        return f"{dataset_key(self.seed)}/{self.run}/{part}{rep}"

    @staticmethod
    def from_id(job_id: str) -> Job:
        key, run, part = job_id.split("/")
        part, _, rep = part.partition("_r")
        seed = None if key == "noise_free" else int(key.removeprefix("seed_"))
        return Job(seed, run,
                   None if part == "all" else int(part.removeprefix("c")),
                   int(rep) if rep else 1)  # fmt: skip

    @property
    def nlive(self) -> int:
        return NLIVE_CELL.get((self.run, self.cell), NLIVE_DEFAULT)

    @property
    def cost(self) -> float:
        return LV_COST if self.run == LV else self.nlive / NLIVE_DEFAULT

    def path(self, suffix: str) -> Path:
        return CACHE / (self.id.replace("/", "__") + suffix)

    @property
    def reuse_source(self) -> Path | None:
        if self.seed == 0 and self.run == "free_a_e0.3" and self.cell in REUSE:
            return NESTED_DIR / REUSE[self.cell]
        return None

    def rng(self) -> np.random.Generator:
        cell = 0 if self.cell is None else self.cell + 1  # SeedSequence entries must be ≥ 0
        seed = 10_000 if self.seed is None else self.seed
        return np.random.default_rng([2026, seed, RUNS.index(self.run), cell, self.rep])


def dataset_key(seed: int | None) -> str:
    return "noise_free" if seed is None else f"seed_{seed}"


def runs_for(seed: int | None) -> list[str]:
    if seed is None:
        return list(NOISE_FREE_RUNS)
    return [r for r in RUNS if seed in RUN_SEEDS[r]]


def jobs_for(seed: int | None, *, with_lv: bool = True) -> list[Job]:
    """This dataset's jobs. The Le Verrier-shape job for seeds other than ``LV_CHECK_SEED`` is
    only listed when ``with_lv`` (i.e. after seed 0 passed the single-mode check)."""
    from fork.inference.nested import cells_for

    out = []
    for run in runs_for(seed):
        if run == LV:
            if seed == LV_CHECK_SEED:
                out += [Job(seed, run, None, 1), Job(seed, run, None, 2)]
            elif with_lv:
                out.append(Job(seed, run, None))
        else:
            out += [Job(seed, run, k) for k in range(len(cells_for(E_MAX[run])))]
    return out


def all_jobs(*, with_lv: bool = True) -> list[Job]:
    return [j for s in (*NOISE_SEEDS, None) for j in jobs_for(s, with_lv=with_lv)]


# ------------------------------------------------------------------ one dataset's problem


def dataset(seed: int | None):
    """(observations, forward context, Uranus inc/raan, (p, q), Uranus box) for a noise seed,
    or for the noise-free truth (``seed=None``: seed 0's epochs and σ, true longitudes)."""
    from fork.data.synthetic import LongitudeObservations, make_synthetic_observations
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.priors import UranusBox
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    inc_raan = uranus_inc_raan_at_epoch(A1_EPOCH_JD_TDB)
    obs, truth = make_synthetic_observations(np.random.default_rng(0 if seed is None else seed))
    if seed is None:
        obs = LongitudeObservations(obs.jd_tdb, truth, obs.sigma_rad, obs.source)
    box = UranusBox.from_uranus_only_fit(obs, ctx, inc_raan)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    return obs, ctx, inc_raan, pq, box


def base_space(run: str, pq, box):
    from fork.inference.priors import free_a_space, held_shape_space
    from fork.model.hidden_mass import LE_VERRIER_1846 as LVE

    if run == LV:
        return held_shape_space(pq, LVE.a_au, LVE.e, LVE.varpi_ecl_j2000_rad(), box)
    return free_a_space(pq, E_MAX[run], box)


def job_space(job: Job, pq, box):
    from fork.inference.nested import cell_space, cells_for

    space = base_space(job.run, pq, box)
    return space if job.cell is None else cell_space(space, cells_for(E_MAX[job.run])[job.cell])


def do_job(job: Job, n_workers: int) -> None:
    from fork.inference.nested import NestedResult, run_nested

    out = job.path(".pkl")
    if out.exists():
        print(f"{job.id}: cached")
        return
    src = job.reuse_source
    if src is not None and src.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, out)
        print(f"{job.id}: reused {src.name}")
        return
    from fork import ledger

    obs, ctx, _, pq, box = dataset(job.seed)
    space = job_space(job, pq, box)
    resumed = job.path(".ckpt").exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(obs, ctx, space, job.rng(), nlive=job.nlive, n_workers=n_workers,
                     progress_log=PROGRESS_LOG, checkpoint_file=job.path(".ckpt"),
                     label=job.id)  # fmt: skip
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    # one likelihood call = one forward integration; a resumed job's calls include the work
    # before its checkpoint while wall/cpu cover this call only, so that row is estimated
    yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append(f"m5b/{job.id}", "nested", n_likelihood_calls=res.ncall,
                  n_integrations=res.ncall, sim_years=res.ncall * yrs,
                  wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0, n_workers=n_workers,
                  estimated=resumed, resumed=resumed, nlive=job.nlive,
                  niter=res.niter)  # fmt: skip
    print(f"{job.id}: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, "
          f"{res.wall_s / 60:.1f} min")  # fmt: skip


# ------------------------------------------------------------------ Le Verrier shape check


def _circ_rel(x: np.ndarray, ref: float) -> np.ndarray:
    return (x - ref + math.pi) % (2.0 * math.pi) - math.pi


def lv_single_mode(r1, r2) -> dict:
    """Is the Le Verrier-shape posterior on seed 0 single-mode, so that one unpartitioned run
    per seed suffices? Both rng seeds must agree on ln Z within the inflated error and on the
    medians of log10 m and λ₀ (≤ 0.5 of the 68 % half-width), and in each run < 0.5 % of the
    mass may lie beyond median ± 8 × half-width in either (a second mode would sit there)."""
    from fork.inference.nested import ERR_INFLATION, weighted_quantiles

    rec, ok = {"runs": []}, True
    meds = []
    for r in (r1, r2):
        q, w = r.samples, r.weights
        w = w / w.sum()
        m, lam = q[:, 4], q[:, 5]  # held-shape layout: Uranus 4, log10 m, λ₀
        c = math.atan2(np.sum(w * np.sin(lam)), np.sum(w * np.cos(lam)))
        lam_rel = _circ_rel(lam, c)
        one = {}
        for name, x in (("log10_m", m), ("lam0_rad", lam_rel)):
            lo, med, hi = weighted_quantiles(x, w)
            hw = 0.5 * (hi - lo)
            far = float(w[np.abs(x - med) > LV_CHECK_SIGMA * hw].sum())
            one[name] = {"median": med + (c if name == "lam0_rad" else 0.0), "half_width": hw,
                         "far_mass": round(far, 6)}  # fmt: skip
            ok &= far < LV_CHECK_MAX_FAR_MASS
        meds.append(one)
        rec["runs"].append(
            {
                "logz": round(r.logz, 4),
                "logzerr": round(r.logzerr, 4),
                **{k: {kk: round(vv, 6) for kk, vv in v.items()} for k, v in one.items()},
            }
        )
    dz = r2.logz - r1.logz
    ez = ERR_INFLATION * math.hypot(r1.logzerr, r2.logzerr)
    rec["delta_logz"], rec["err_eff"] = round(dz, 4), round(ez, 4)
    ok &= abs(dz) <= ez
    for name in ("log10_m", "lam0_rad"):
        a, b = meds[0][name], meds[1][name]
        d = a["median"] - b["median"]
        if name == "lam0_rad":
            d = float(_circ_rel(np.array(d), 0.0))
        sig = 0.5 * (a["half_width"] + b["half_width"])
        rec[f"{name}_median_diff_sigma"] = round(d / sig, 3)
        ok &= abs(d / sig) <= LV_CHECK_MEDIAN_TOL
    rec["single_mode"] = bool(ok)
    return rec


def lv_check_state() -> dict | None:
    """The single-mode check's record once both seed-0 runs exist (cached as lv_check.json)."""
    from fork.inference.nested import NestedResult

    path = CACHE / "lv_check.json"
    if path.exists():
        return json.loads(path.read_text())
    p1, p2 = (Job(LV_CHECK_SEED, LV, None, r).path(".pkl") for r in (1, 2))
    if not (p1.exists() and p2.exists()):
        return None
    rec = lv_single_mode(NestedResult.load(p1), NestedResult.load(p2))
    path.write_text(json.dumps(rec, indent=1) + "\n")
    return rec


# ------------------------------------------------------------------ scheduler


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with DRIVER_LOG.open("a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def ordered(jobs: list[Job]) -> list[Job]:
    """Le Verrier check runs first (they decide four more jobs), then longest first; the
    noise-free jobs come last, after everything A3b needs."""
    return sorted(
        jobs, key=lambda j: (j.seed is None, not (j.run == LV and j.seed == LV_CHECK_SEED), -j.cost)
    )


def run_all(n_workers: int, per_job: int) -> None:
    slots = max(1, n_workers // per_job)
    lv = lv_check_state()
    lv_known = lv is not None
    pending = ordered([j for j in all_jobs(with_lv=bool(lv and lv["single_mode"]))
                       if not j.path(".pkl").exists()])  # fmt: skip
    _log(f"run: {len(pending)} jobs pending, {slots} at a time × {per_job} workers")
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
                if failures[jid] < 2:  # one retry; resumes from its checkpoint
                    pending.insert(0, job)
        if not lv_known:
            lv = lv_check_state()
            if lv is not None:
                lv_known = True
                if lv["single_mode"]:
                    extra = [j for j in all_jobs() if j.run == LV and j.seed != LV_CHECK_SEED
                             and not j.path(".pkl").exists()]  # fmt: skip
                    pending = ordered(pending + extra)
                    _log(f"Le Verrier shape single-mode on seed {LV_CHECK_SEED}: "
                         f"queued {len(extra)} unpartitioned runs")  # fmt: skip
                else:
                    _log(
                        f"Le Verrier shape NOT single-mode on seed {LV_CHECK_SEED}: seeds "
                        f"1–4 not run (needs partitioning); see {CACHE / 'lv_check.json'}"
                    )
        while pending and len(running) < slots:
            job = pending.pop(0)
            log = job.path(".log")
            log.parent.mkdir(parents=True, exist_ok=True)
            with log.open("a") as fh:
                proc = subprocess.Popen(
                    [sys.executable, "-m", "fork.export.m5b", "job", job.id, "--n-workers",
                     str(per_job)], stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)  # fmt: skip
            running[job.id] = (proc, time.time(), proc.pid)
            _log(f"start {job.id} (nlive {job.nlive}, pid {proc.pid})")
        time.sleep(30)
    _log(f"run finished; failures: {failures or 'none'}")


def plan() -> None:
    lv = lv_check_state()
    jobs = all_jobs(with_lv=lv is None or lv["single_mode"])
    done = [j for j in jobs if j.path(".pkl").exists()]
    reused = [j for j in jobs if j.reuse_source is not None]
    cost = sum(j.cost for j in jobs if j not in done and j not in reused)
    for j in ordered(jobs):
        state = "done" if j in done else ("reuse" if j in reused else
                ("resumable" if j.path(".ckpt").exists() else "pending"))  # fmt: skip
        print(f"{j.id:36s} nlive {j.nlive:5d}  {state}")
    print(f"{len(jobs)} jobs, {len(done)} done, {len(reused)} reused; remaining cost "
          f"{cost:.1f} cell-units (nlive-1000 cell ≈ 75 min on 4 workers)")  # fmt: skip
    print("Le Verrier check:", "pending" if lv is None else lv["single_mode"])


# ------------------------------------------------------------------ summary


def _pct(x: np.ndarray) -> list[float]:
    return [float(v) for v in np.percentile(x, [16.0, 50.0, 84.0])]


def _summary_1d(x: np.ndarray, truth: float | None, nd: int) -> dict:
    lo, med, hi = _pct(x)
    out = {"median": round(med, nd), "lo68": round(lo, nd), "hi68": round(hi, nd),
           "minus": round(med - lo, nd), "plus": round(hi - med, nd)}  # fmt: skip
    if truth is not None:
        out["truth"] = round(truth, nd)
        out["truth_in_68"] = bool(lo <= truth <= hi)
    return out


def _lon_summary(lon: np.ndarray, truth: float | None = None) -> dict:
    from fork.inference.mcmc import circular_interval_deg

    lo, med, hi = circular_interval_deg(lon)

    def rel(v: float) -> float:
        return (v - med + 180.0) % 360.0 - 180.0

    out = {"median": round(med, 3), "lo68": round(lo, 3), "hi68": round(hi, 3),
           "minus": round(-rel(lo), 3), "plus": round(rel(hi), 3),
           "width68": round(rel(hi) - rel(lo), 3)}  # fmt: skip
    if truth is not None:
        out |= {"truth": truth, "median_minus_truth": round(-rel(truth), 3),
                "truth_in_68": bool(rel(lo) <= rel(truth) <= rel(hi))}  # fmt: skip
    return out


def _heliocentric_lon_of_date_deg(theta: np.ndarray, ctx, jd: float) -> float:
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.model.forward import predict_states_ecl_j2000

    pos, _ = predict_states_ecl_j2000(theta, [jd], ctx)
    r = precession_matrix_ecl_j2000_to_ecl_of_date(jd) @ pos[0, 3]
    return math.degrees(math.atan2(r[1], r[0])) % 360.0


def point_fit(seed: int | None, run: str, obs, ctx, inc_raan, n_workers: int) -> dict:
    """Bounded 24-start multistart point fit for this posterior (cached)."""
    from fork.inference import fit
    from fork.model.hidden_mass import LE_VERRIER_1846 as LVE

    path = CACHE / "points" / f"{dataset_key(seed)}__{run}.json"
    if path.exists():
        return json.loads(path.read_text())
    if run == LV:
        varpi, s = LVE.varpi_ecl_j2000_rad(), math.sqrt(LVE.e)
        kw = dict(space=fit.PerturberSpace((None, LVE.a_au, s * math.cos(varpi),
                                            s * math.sin(varpi), None)), mode=LV)  # fmt: skip
    else:
        kw = dict(e_max=E_MAX[run])
    ms = fit.fit_perturber_multistart(obs, ctx, inc_raan, n_workers=n_workers, **kw)
    _SUMMARY_COUNT["point_fit_nfev"] += sum(int(f.nfev) for f in ms.fits)
    b = ms.best
    rec = {"chi2": round(b.chi2, 4), "dof": b.dof, "mass_earth": round(b.mass_earth, 3),
           "a_au": round(float(b.theta[7]), 3), "e": round(float(b.theta[8]), 4)}  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rec) + "\n")
    return rec


def profile_intervals(seed: int | None) -> dict | None:
    """The committed M5a profile's Δχ² ≤ 1 intervals for this dataset (both modes), or None if
    the dataset was not profiled."""
    path = SITE_DATA / "m5a_profile.json"
    if not path.exists():
        return None
    for d in json.loads(path.read_text())["datasets"]:
        if d["seed"] == seed:
            return {mode: {"dchi2_le_1_intervals_deg": m["dchi2_le_1_intervals_deg"],
                           "overall_min_lon_ecl_of_date_deg":
                               m["overall_min"]["lon_ecl_of_date_deg"],
                           "dchi2_at_true_neptune": m["dchi2_at_true_neptune"]}
                    for mode, m in d["modes"].items()}  # fmt: skip
    return None


def summarize_run(seed: int | None, run: str, data, truth_p: dict, n_workers: int) -> dict | None:
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult, cells_for
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
    from fork.model.hidden_mass import LE_VERRIER_1846 as LVE

    obs, ctx, inc_raan, _, _ = data
    jobs = [j for j in jobs_for(seed) if j.run == run and j.rep == 1]
    if not jobs or not all(j.path(".pkl").exists() for j in jobs):
        return None
    res = [NestedResult.load(j.path(".pkl")) for j in jobs]
    geo = GalleGeometry.from_horizons()
    rng = np.random.default_rng([7, seed, RUNS.index(run)])
    nep = truth_p["galle_lon_ecl_of_date_deg"]

    def lon_of(space, q: np.ndarray, n: int) -> np.ndarray:
        idx = rng.choice(len(q), size=min(n, len(q)), replace=False)
        _SUMMARY_COUNT["derived_integrations"] += len(idx)
        return np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo))
                         for v in q[idx]])  # fmt: skip

    space = res[0].space
    if run == LV:
        logz, err_q = res[0].logz, res[0].logzerr
        evidence = {"logz": round(logz, 4), "logzerr_quoted": round(err_q, 4),
                    "logzerr": round(ERR_INFLATION * err_q, 4), "partitioned": False}  # fmt: skip
        eq = res[0].equal_weight(rng)
    else:
        p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
        cells = cells_for(E_MAX[run])
        cell_recs = []
        for k, (r, s, se) in enumerate(zip(res, p.shares, p.share_errs, strict=True)):
            lon_k = lon_of(r.space, r.equal_weight(rng), CELL_LON_SAMPLES)
            cell_recs.append({
                "a_au": list(cells[k][0]), "e": None if cells[k][1] is None else list(cells[k][1]),
                "prior_weight": round(r.space.slice_weight, 6), "nlive": r.nlive,
                "logz_k": round(r.logz, 4), "logzerr_k_quoted": round(r.logzerr, 4),
                "share": round(float(s), 5), "share_err": round(float(se), 5),
                "niter": r.niter, "ncall": r.ncall, "wall_min": round(r.wall_s / 60, 1),
                "chi2_min": round(-2 * r.max_logl, 4), "galle_lon_ecl_of_date_deg":
                    _lon_summary(lon_k)})  # fmt: skip
        evidence = {"logz": round(p.logz, 4), "logzerr": round(p.logzerr, 4),
                    "logzerr_quoted": round(PartitionedResult(tuple(res)).logzerr, 4),
                    "partitioned": True, "cells": cell_recs}  # fmt: skip
        eq = p.equal_weight(rng)
    thetas = np.array([space.theta(v) for v in eq])
    lon = lon_of(space, eq, DERIVED_MAX_SAMPLES)
    pt = point_fit(seed, run, obs, ctx, inc_raan, n_workers)
    chi2_post = min(-2 * r.max_logl for r in res)
    rec = {
        "evidence": evidence, "n_samples": len(eq), "n_derived": len(lon),
        "chi2_posterior_min": round(chi2_post, 4), "point_fit": pt,
        "chi2_best": round(min(chi2_post, pt["chi2"]), 4), "dof": pt["dof"],
        "mass_earth": _summary_1d(10.0 ** thetas[:, 6], truth_p["mass_earth"], 3),
        "galle_lon_ecl_of_date_deg": _lon_summary(lon, nep),
    }  # fmt: skip
    if run != LV:
        rec["a_au"] = _summary_1d(thetas[:, 7], truth_p["a_au"], 3)
        rec["e"] = _summary_1d(thetas[:, 8], truth_p["e"], 4)
    else:
        sub = thetas[rng.choice(len(thetas), size=min(2000, len(thetas)), replace=False)]
        helio = np.array([_heliocentric_lon_of_date_deg(th, ctx, LVE.epoch_jd) for th in sub])
        _SUMMARY_COUNT["derived_integrations"] += len(sub)
        rec["heliocentric_true_lon_1847_jan1_deg"] = {
            **{k: v for k, v in _lon_summary(helio).items() if k in ("median", "lo68", "hi68")},
            "le_verrier_published": round(LVE.true_lon_ecl_of_date_deg, 4)}  # fmt: skip
    return rec


def calibration(datasets: dict, run_name: str) -> dict:
    noisy = [d for d in datasets.values() if d["noise"] and run_name in d["runs"]]
    hits = [d["runs"][run_name]["galle_lon_ecl_of_date_deg"]["truth_in_68"] for d in noisy]
    widths = [d["runs"][run_name]["galle_lon_ecl_of_date_deg"]["width68"] for d in noisy]
    frac = float(np.mean(hits)) if hits else float("nan")
    return {"n_seeds": len(hits), "n_truth_in_68": int(sum(hits)), "fraction": round(frac, 3),
            "pass_50_85": bool(len(hits) >= 10 and A3B_RANGE[0] <= frac <= A3B_RANGE[1]),
            "width68_deg_median": round(float(np.median(widths)), 3) if widths else None,
            "width68_deg_range": [round(min(widths), 3), round(max(widths), 3)]
            if widths else None}  # fmt: skip


def truth_params() -> dict:
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.model import bodies
    from fork.model.forward import perturber_theta_from_state
    from fork.model.hidden_mass import true_neptune_geocentric_longitude

    th = perturber_theta_from_state(about("neptune", [A1_EPOCH_JD_TDB])[0],
                                    bodies.NEPTUNE.mass_msun)  # fmt: skip
    return {"mass_earth": round(10 ** th[0], 4), "a_au": round(th[1], 4), "e": round(th[2], 5),
            "galle_lon_ecl_of_date_deg": round(math.degrees(
                true_neptune_geocentric_longitude().lon_ecl_of_date_rad), 4)}  # fmt: skip


def summarize_dataset(seed: int | None, truth_p: dict, n_workers: int) -> dict:
    data = dataset(seed)
    out = {"seed": seed, "noise": seed is not None, "runs": {}, "missing": []}
    for run in runs_for(seed):
        rec = summarize_run(seed, run, data, truth_p, n_workers)
        if rec is None:
            out["missing"].append(run)
        else:
            out["runs"][run] = rec
    if LV in out["runs"] and "free_a_e0.3" in out["runs"]:
        a, b = out["runs"]["free_a_e0.3"], out["runs"][LV]
        out["le_verrier_shape_vs_free_a"] = {
            "delta_chi2": round(b["chi2_best"] - a["chi2_best"], 4),
            "delta_ln_likelihood": round(-0.5 * (b["chi2_best"] - a["chi2_best"]), 4),
            "ln_bayes_factor": round(b["evidence"]["logz"] - a["evidence"]["logz"], 4),
            "ln_bayes_factor_err": round(math.hypot(b["evidence"]["logzerr"],
                                                    a["evidence"]["logzerr"]), 4)}  # fmt: skip
    out["profile_m5a"] = profile_intervals(seed)
    print(
        f"{dataset_key(seed)}: "
        + ", ".join(
            f"{r} lon {v['galle_lon_ecl_of_date_deg']['median']:.2f} "
            f"[{v['galle_lon_ecl_of_date_deg']['lo68']:.2f}, "
            f"{v['galle_lon_ecl_of_date_deg']['hi68']:.2f}] "
            f"in68 {v['galle_lon_ecl_of_date_deg']['truth_in_68']}"
            for r, v in out["runs"].items()
        )
        + (f"; missing {out['missing']}" if out["missing"] else ""),
        flush=True,
    )
    return out


def noise_free_check(nf: dict) -> dict:
    """Prior-bias measurement on the noise-free posterior (free a, e ≤ 0.3): for the 1846
    longitude, m, a and e, whether truth lies in the 68 % interval, and the median's offset from
    truth in units of the 68 % half-width on the side facing truth (median − truth)/σ."""
    run = "free_a_e0.3"
    if run not in nf["runs"]:
        return {"status": "missing"}
    r = nf["runs"][run]
    out = {"status": "done", "run": run}
    for name in ("galle_lon_ecl_of_date_deg", "mass_earth", "a_au", "e"):
        x = r[name]
        off = x["median_minus_truth"] if "median_minus_truth" in x else x["median"] - x["truth"]
        sig = x["minus"] if off > 0 else x["plus"]  # half-width on the side facing truth
        out[name] = {"median": x["median"], "truth": x["truth"],
                     "truth_in_68": x["truth_in_68"], "median_minus_truth": round(off, 4),
                     "offset_sigma": round(off / sig, 3) if sig > 0 else None}  # fmt: skip
    return out


def summarize(n_workers: int) -> None:
    """A3b first (noise seeds; the JSON is written once they are summarized), then the
    noise-free posterior check is appended and the JSON rewritten."""
    from fork import ledger
    from fork.inference.nested import BOUND, DLOGZ, ERR_INFLATION, SAMPLE, cells_for

    t0, c0 = time.time(), ledger.cpu_now()
    truth_p = truth_params()
    lv = lv_check_state()
    datasets = {dataset_key(s): summarize_dataset(s, truth_p, n_workers) for s in NOISE_SEEDS}
    payload = {
        "description": "M5b partitioned nested-sampling (dynesty) posteriors on synthetic data "
        "(Horizons Uranus incl. real Neptune + noise; noise_free: the true longitudes). "
        "Intervals are 16/50/84 % marginal percentiles of the combined posterior; the 1846 "
        "longitude is the perturber's light-time-corrected geocentric ecliptic-of-date "
        "longitude on 1846-09-23. ln Z errors are dynesty's quoted errors × err_inflation "
        "(measured on seed 0). Raw runs are not committed (data/cache/nested/m5b/).",
        "priors": {"log10_m_earth": [0.0, 2.0], "a_au_log_uniform": [20.0, 60.0],
                   "e_uniform_max": E_MAX, "angles": "uniform",
                   "uranus": "box ±50 scaled σ around the Uranus-only fit"},
        "le_verrier_shape": "a = 36.15 AU, e = 0.10761, ϖ = 284°45′ (equinox 1847.0) held",
        "sampler": {"package": "dynesty", "bound": BOUND, "sample": SAMPLE, "dlogz": DLOGZ,
                    "nlive_default": NLIVE_DEFAULT,
                    "nlive_cells": {f"{r}/c{c}": n for (r, c), n in NLIVE_CELL.items()},
                    "err_inflation": ERR_INFLATION,
                    "cells": {r: [{"a_au": list(a), "e": None if e is None else list(e)}
                                  for a, e in cells_for(E_MAX[r])] for r in E_MAX},
                    "le_verrier_shape": "unpartitioned" if lv and lv["single_mode"] else
                    "seed 0 only (single-mode check not passed)" if lv else "check pending"},
        "le_verrier_single_mode_check": lv,
        "run_seeds": {**{r: list(v) for r, v in RUN_SEEDS.items()},
                      "noise_free": list(NOISE_FREE_RUNS)},
        "profile_m5a": "per dataset: M5a profile Δχ² ≤ 1 intervals (site/data/"
        "m5a_profile.json) where that dataset was profiled, else null",
        "truth": truth_p,
        "datasets": datasets,
        "a3b_calibration": {r: calibration(datasets, r) for r in ("free_a_e0.3", "free_a_e0.1")},
        "noise_free_check": {"status": "pending"},
    }  # fmt: skip
    write_json(OUT_NAME, payload)
    print(json.dumps(payload["a3b_calibration"], indent=1), flush=True)
    nf = summarize_dataset(None, truth_p, n_workers)
    datasets["noise_free"] = nf
    payload["noise_free_check"] = noise_free_check(nf)
    write_json(OUT_NAME, payload)
    print(json.dumps(payload["noise_free_check"], indent=1))
    print(f"wrote {SITE_DATA / OUT_NAME}")
    n_fit, n_der = _SUMMARY_COUNT["point_fit_nfev"], _SUMMARY_COUNT["derived_integrations"]
    yrs = ledger.sim_years_per_integration(dataset(0)[0].jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append("m5b/summarize", "summarize", n_likelihood_calls=n_fit,
                  n_integrations=n_fit + n_der, sim_years=(n_fit + n_der) * yrs,
                  wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0, n_workers=n_workers,
                  estimated=True, derived_integrations=n_der, point_fit_nfev=n_fit,
                  note="point-fit nfev excludes finite-difference Jacobian calls (lower "
                  "bound); cached point fits count 0")  # fmt: skip


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
    p_sum = sub.add_parser("summarize")
    p_sum.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    if args.cmd == "plan":
        plan()
    elif args.cmd == "run":
        run_all(args.n_workers, args.per_job)
    elif args.cmd == "job":
        do_job(Job.from_id(args.job_id), args.n_workers)
    else:
        summarize(args.n_workers)


if __name__ == "__main__":
    main()
