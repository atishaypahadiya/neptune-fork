"""M5b: partitioned dynesty posteriors on synthetic data → ``site/data/m5b_posteriors.json``.

    uv run python -m fork.export.m5b plan                  # jobs and their state
    uv run python -m fork.export.m5b run [--n-workers 12 --per-job 4]   # all jobs (a day)
    uv run python -m fork.export.m5b job JOB_ID [--n-workers 4]         # one nested run
    uv run python -m fork.export.m5b summarize             # the JSON, with A3b

Posteriors (PHASE1_SPEC M5b, as amended after the partitioned validation):
  - ``free_a_e0.3``: free a, e ∈ [0, 0.3], noise seeds 0–9, partitioned into the four cells
    of ``nested.PARTITION_CELLS``, one rng seed per cell. Cell [20, 25] runs at nlive 2000
    (the thin mode-B pocket was under-sampled at 1000), the rest at 1000.
  - ``free_a_e0.1``: free a, e ∈ [0, 0.1], seeds 0–9: the three cells that can hold e ≤ 0.1
    (``nested.cells_for``), nlive 1000. Seeds 5–9 were added after A3b so the prior-sensitivity
    result has its own 10-seed coverage count; it is a sensitivity result, not the headline.
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
    "free_a_e0.1": NOISE_SEEDS,
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


E_UPPER_PCT = (68.0, 95.0)  # one-sided upper limits on e
E_LIMIT_KEYS = ("upper68_one_sided", "upper95_one_sided", "truth_below_upper68",
                "truth_below_upper95", "prior_median")  # fmt: skip


def _e_summary(e: np.ndarray, truth: float | None, e_max: float) -> dict:
    """e is bounded below by 0, and the true e (0.0039) sits on that edge, so a central 16–84 %
    interval excludes it by construction. Alongside the central interval (kept for reference):
    one-sided 68 % and 95 % upper limits, and the prior median (e uniform on [0, e_max]) beside
    the posterior median, to show how weakly the data constrain e."""
    out = _summary_1d(e, truth, 4)
    up68, up95 = (float(v) for v in np.percentile(e, E_UPPER_PCT))
    out |= {"upper68_one_sided": round(up68, 4), "upper95_one_sided": round(up95, 4),
            "prior_median": round(e_max / 2.0, 4)}  # fmt: skip
    if truth is not None:
        out |= {"truth_below_upper68": bool(truth <= up68),
                "truth_below_upper95": bool(truth <= up95)}  # fmt: skip
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
    prof = json.loads(path.read_text())
    for d in prof["datasets"]:
        if d["seed"] == seed:
            return {mode: {"dchi2_le_1_intervals_deg": m["dchi2_le_1_intervals_deg"],
                           "overall_min_lon_ecl_of_date_deg":
                               m["overall_min"]["lon_ecl_of_date_deg"],
                           "dchi2_at_true_neptune": m["dchi2_at_true_neptune"],
                           "grid_deg": prof["grid_deg"][:2],
                           "dchi2_curve": {"lon_deg": m["profile"]["lon_ecl_of_date_deg"],
                                           "dchi2": m["profile"]["dchi2"],
                                           "chi2": m["profile"]["chi2"]}}
                    for mode, m in d["modes"].items()}  # fmt: skip
    return None


PROFILE_POSTERIOR_FLAG_DEG = 2.0
# a draw counts as below the profile if its χ² is lower than the profile's, linearly interpolated
# on the 1° grid, by more than this (the slack the interpolation needs)
DRAW_BELOW_PROFILE_TOL = 0.1


def profile_vs_posterior(intervals: list, lon: dict, grid_deg: list | None = None,
                         flag_deg: float = PROFILE_POSTERIOR_FLAG_DEG, *,
                         curve: dict | None = None, profile_min_deg: float | None = None,
                         sampled_bins: set[int] | None = None,
                         draws: dict | None = None) -> dict:  # fmt: skip
    """The profile's Δχ² ≤ 1 interval (hull of its pieces) next to the posterior 68 % interval
    of the 1846 longitude; flagged if either end differs by more than ``flag_deg``. Same model
    (free a, e ≤ 0.3); not the same quantity (profile likelihood vs marginal posterior). An end
    on the profile grid's edge is a bound, not a measured end (profile_end_at_grid_edge).

    Each flag is classified (``class``; needs ``curve``, ``profile_min_deg``, ``sampled_bins``;
    ``volume_beyond`` also needs ``draws``, the cells' posterior draws as {lon_deg, chi2}):
      - ``missed_mode``: the profile minimum is outside every posterior-sampled region, i.e. no
        cell's own posterior draws fall in its 1° bin or the two neighbouring bins;
      - ``volume_inside``: otherwise, if the profile is flat (Δχ² ≤ 1 at every grid point)
        across each flagged gap between the two ends, so the posterior leans within a flat band;
      - ``volume_beyond``: otherwise, if every gap that is not flat is one where the posterior's
        68 % end lies beyond the profile's Δχ² ≤ 1 end, and no posterior draw has χ² below the
        profile at its longitude (by more than ``DRAW_BELOW_PROFILE_TOL``): the posterior
        reaches where the best fit is worse but more prior volume fits, not a branch the
        profile missed;
      - ``other``: anything else.
    """

    def d(x: float, y: float) -> float:
        return (x - y + 180.0) % 360.0 - 180.0

    lo = intervals[0][0]
    hi = intervals[-1][1]
    dlo, dhi = d(lo, lon["lo68"]), d(hi, lon["hi68"])
    flag = bool(max(abs(dlo), abs(dhi)) > flag_deg)
    out = {"profile_dchi2_le_1_deg": [lo, hi], "n_profile_pieces": len(intervals),
           "posterior_68_deg": [lon["lo68"], lon["hi68"]],
           "delta_lo_deg": round(dlo, 3), "delta_hi_deg": round(dhi, 3),
           "flag_gt_deg": flag_deg, "flag": flag,
           "profile_end_at_grid_edge": None if grid_deg is None else
           [abs(lo - grid_deg[0]) < 1e-6, abs(hi - grid_deg[1]) < 1e-6],
           "class": None}  # fmt: skip
    if not flag or curve is None or profile_min_deg is None or sampled_bins is None:
        return out
    b = math.floor(profile_min_deg % 360.0)
    min_sampled = any(((b + k) % 360) in sampled_bins for k in (-1, 0, 1))
    g_lon, g_d = np.asarray(curve["lon_deg"]), np.asarray(curve["dchi2"])
    gaps, flat, beyond = [], [], []
    for prof_end, delta, sign in ((lo, dlo, 1.0), (hi, dhi, -1.0)):
        if abs(delta) <= flag_deg:
            continue
        a, c = sorted((prof_end, prof_end - delta))  # prof_end − delta = post_end, unwrapped
        sel = (g_lon >= a - 1e-9) & (g_lon <= c + 1e-9)
        gaps.append([round(a, 3), round(c, 3)])
        flat.append(bool(sel.any() and np.all(g_d[sel] <= 1.0)))
        beyond.append(bool(sign * delta > 0))  # posterior end outside the profile's interval
    out["profile_min_deg"] = profile_min_deg
    out["profile_min_sampled"] = min_sampled
    out["gaps_deg"] = gaps
    out["gap_flat"] = flat
    out["gap_posterior_beyond"] = beyond
    below = None
    if draws is not None and "chi2" in curve:
        dl = np.asarray(draws["lon_deg"], dtype=float) % 360.0
        inside = (dl >= g_lon[0]) & (dl <= g_lon[-1])
        margin = (np.asarray(draws["chi2"], dtype=float)[inside]
                  - np.interp(dl[inside], g_lon, np.asarray(curve["chi2"])))  # fmt: skip
        k = int(np.argmin(margin))
        below = bool(margin[k] < -DRAW_BELOW_PROFILE_TOL)
        out["draws_vs_profile"] = {
            "n_draws": int(inside.sum()), "min_draw_minus_profile_chi2": round(float(margin[k]), 3),
            "at_lon_deg": round(float(dl[inside][k]), 3), "tol": DRAW_BELOW_PROFILE_TOL,
            "any_draw_below_profile": below}  # fmt: skip
    if not min_sampled:
        cls = "missed_mode"
    elif all(flat):
        cls = "volume_inside"
    elif all(f or b for f, b in zip(flat, beyond, strict=True)) and below is False:
        cls = "volume_beyond"
    else:
        cls = "other"
    out["class"] = cls
    return out


def _sampled_bins(lon_deg: np.ndarray) -> list[int]:
    """Whole-degree bins (0–359) holding at least one posterior draw."""
    return sorted({int(v) for v in np.floor(np.asarray(lon_deg) % 360.0)})


def _chi2_of_rows(r, rows: np.ndarray) -> np.ndarray:
    """χ² = −2 ln L of equal-weight draws: they are copies of rows of ``r.samples``."""
    logl = dict(zip((v.tobytes() for v in r.samples), np.asarray(r.results["logl"]), strict=True))
    return np.array([-2.0 * logl[v.tobytes()] for v in rows])


def summary_rng(seed: int | None, run: str) -> np.random.Generator:
    """Subsampling rng for a run's summary (noise-free ↦ 10 000, as in ``Job.rng``)."""
    return np.random.default_rng([7, 10_000 if seed is None else seed, RUNS.index(run)])


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
    rng = summary_rng(seed, run)
    nep = truth_p["galle_lon_ecl_of_date_deg"]

    def lon_of(space, q: np.ndarray, n: int, rows_out: list | None = None) -> np.ndarray:
        idx = rng.choice(len(q), size=min(n, len(q)), replace=False)
        if rows_out is not None:
            rows_out.append(q[idx])
        _SUMMARY_COUNT["derived_integrations"] += len(idx)
        return np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo))
                         for v in q[idx]])  # fmt: skip

    space = res[0].space
    if run == LV:
        logz, err_q = res[0].logz, res[0].logzerr
        evidence = {"logz": round(logz, 4), "logzerr_quoted": round(err_q, 4),
                    "logzerr": round(ERR_INFLATION * err_q, 4),
                    "logzerr_inflation": ERR_INFLATION, "partitioned": False}  # fmt: skip
        eq = res[0].equal_weight(rng)
    else:
        p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
        cells = cells_for(E_MAX[run])
        cell_recs, draws_lon, draws_chi2 = [], [], []
        for k, (r, s, se) in enumerate(zip(res, p.shares, p.share_errs, strict=True)):
            rows: list = []
            lon_k = lon_of(r.space, r.equal_weight(rng), CELL_LON_SAMPLES, rows)
            draws_lon.append(lon_k)
            draws_chi2.append(_chi2_of_rows(r, rows[0]))
            cell_recs.append({
                "a_au": list(cells[k][0]), "e": None if cells[k][1] is None else list(cells[k][1]),
                "prior_weight": round(r.space.slice_weight, 6), "nlive": r.nlive,
                "logz_k": round(r.logz, 4), "logzerr_k_quoted": round(r.logzerr, 4),
                "logzerr_k": round(ERR_INFLATION * r.logzerr, 4),
                "share": round(float(s), 5), "share_err": round(float(se), 5),
                "niter": r.niter, "ncall": r.ncall, "wall_min": round(r.wall_s / 60, 1),
                "chi2_min": round(-2 * r.max_logl, 4), "galle_lon_ecl_of_date_deg":
                    _lon_summary(lon_k),
                "lon_bins_sampled_deg": _sampled_bins(lon_k)})  # fmt: skip
        evidence = {"logz": round(p.logz, 4), "logzerr": round(p.logzerr, 4),
                    "logzerr_quoted": round(PartitionedResult(tuple(res)).logzerr, 4),
                    "logzerr_inflation": ERR_INFLATION, "share_err_inflation": ERR_INFLATION,
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
        # not written to the JSON: summarize_dataset pops it for profile_vs_posterior
        rec["_cell_draws"] = {"lon_deg": np.concatenate(draws_lon),
                              "chi2": np.concatenate(draws_chi2)}  # fmt: skip
        rec["a_au"] = _summary_1d(thetas[:, 7], truth_p["a_au"], 3)
        rec["e"] = _e_summary(thetas[:, 8], truth_p["e"], E_MAX[run])
    else:
        sub = thetas[rng.choice(len(thetas), size=min(2000, len(thetas)), replace=False)]
        helio = np.array([_heliocentric_lon_of_date_deg(th, ctx, LVE.epoch_jd) for th in sub])
        _SUMMARY_COUNT["derived_integrations"] += len(sub)
        rec["heliocentric_true_lon_1847_jan1_deg"] = {
            **{k: v for k, v in _lon_summary(helio).items() if k in ("median", "lo68", "hi68")},
            "le_verrier_published": round(LVE.true_lon_ecl_of_date_deg, 4)}  # fmt: skip
    return rec


def offset_sigma(x: dict) -> tuple[float, float | None]:
    """(median − truth, that offset over the 68 % half-width on the side facing truth)."""
    off = x["median_minus_truth"] if "median_minus_truth" in x else x["median"] - x["truth"]
    sig = x["minus"] if off > 0 else x["plus"]
    return off, (off / sig if sig > 0 else None)


def binomial_calibration(k: int, n: int, p: float = 0.68) -> dict:
    """How surprising k of n hits is for a calibrated 68 % interval: exact two-sided binomial
    p-value (sum of outcomes no more likely than k, scipy ``binomtest``), and the probability a
    calibrated method lands in the spec's A3b range at all."""
    from scipy.stats import binom, binomtest

    lo, hi = math.ceil(A3B_RANGE[0] * n - 1e-9), math.floor(A3B_RANGE[1] * n + 1e-9)
    in_range = float(binom.cdf(hi, n, p) - binom.cdf(lo - 1, n, p))
    return {"p_true": p, "pmf_k": round(float(binom.pmf(k, n, p)), 4),
            "p_value_two_sided": round(float(binomtest(k, n, p).pvalue), 4),
            "spec_range_counts": [lo, hi],
            "p_calibrated_in_spec_range": round(in_range, 4)}  # fmt: skip


def calibration(datasets: dict, run_name: str) -> dict:
    noisy = [d for d in datasets.values() if d["noise"] and run_name in d["runs"]]
    lons = [d["runs"][run_name]["galle_lon_ecl_of_date_deg"] for d in noisy]
    hits = [x["truth_in_68"] for x in lons]
    widths = [x["width68"] for x in lons]
    frac = float(np.mean(hits)) if hits else float("nan")
    per_seed = []
    for d, x in zip(noisy, lons, strict=True):
        off, z = offset_sigma(x)
        per_seed.append({"seed": d["seed"], "truth_in_68": x["truth_in_68"],
                         "median_minus_truth_deg": round(off, 3),
                         "offset_sigma": None if z is None else round(z, 3)})  # fmt: skip
    zs = [r["offset_sigma"] for r in per_seed if r["offset_sigma"] is not None]
    return {"n_seeds": len(hits), "n_truth_in_68": int(sum(hits)), "fraction": round(frac, 3),
            "pass_50_85": bool(len(hits) >= 10 and A3B_RANGE[0] <= frac <= A3B_RANGE[1]),
            "binomial": binomial_calibration(int(sum(hits)), len(hits)) if hits else None,
            "per_seed": per_seed,
            "offset_sigma_rms": round(float(np.sqrt(np.mean(np.square(zs)))), 3) if zs else None,
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
    draws = {run: rec.pop("_cell_draws", None) for run, rec in out["runs"].items()}
    if LV in out["runs"] and "free_a_e0.3" in out["runs"]:
        a, b = out["runs"]["free_a_e0.3"], out["runs"][LV]
        out["le_verrier_shape_vs_free_a"] = {
            "delta_chi2": round(b["chi2_best"] - a["chi2_best"], 4),
            "delta_ln_likelihood": round(-0.5 * (b["chi2_best"] - a["chi2_best"]), 4),
            "ln_bayes_factor": round(b["evidence"]["logz"] - a["evidence"]["logz"], 4),
            "ln_bayes_factor_err": round(math.hypot(b["evidence"]["logzerr"],
                                                    a["evidence"]["logzerr"]), 4),
            "ln_bayes_factor_err_inflation": b["evidence"]["logzerr_inflation"]}  # fmt: skip
    out["profile_m5a"] = profile_intervals(seed)
    if out["profile_m5a"] and "free_a_e0.3" in out["runs"]:
        pf, run_rec = out["profile_m5a"]["free_a"], out["runs"]["free_a_e0.3"]
        bins = {b for c in run_rec["evidence"]["cells"] for b in c["lon_bins_sampled_deg"]}
        out["profile_vs_posterior_free_a"] = profile_vs_posterior(
            pf["dchi2_le_1_intervals_deg"], run_rec["galle_lon_ecl_of_date_deg"],
            pf["grid_deg"], curve=pf["dchi2_curve"],
            profile_min_deg=pf["overall_min_lon_ecl_of_date_deg"], sampled_bins=bins,
            draws=draws["free_a_e0.3"],
        )  # fmt: skip
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
    truth in units of the 68 % half-width on the side facing truth (median − truth)/σ. For e
    (truth on the e = 0 boundary) also the one-sided upper limits and the prior median."""
    run = "free_a_e0.3"
    if run not in nf["runs"]:
        return {"status": "missing"}
    r = nf["runs"][run]
    out = {"status": "done", "run": run}
    for name in ("galle_lon_ecl_of_date_deg", "mass_earth", "a_au", "e"):
        x = r[name]
        off, z = offset_sigma(x)  # half-width on the side facing truth
        out[name] = {"median": x["median"], "truth": x["truth"],
                     "truth_in_68": x["truth_in_68"], "median_minus_truth": round(off, 4),
                     "offset_sigma": None if z is None else round(z, 3)}  # fmt: skip
    e = r["e"]
    out["e"] |= {k: e[k] for k in E_LIMIT_KEYS}
    out["e"]["note"] = (
        "e ≥ 0 and the true e sits on that boundary: a central 68 % interval excludes a boundary "
        "truth by construction, so the one-sided upper limits are the check that applies. The "
        "posterior median against the prior median (uniform e) shows how weakly e is constrained."
    )
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
        "longitude on 1846-09-23. e also carries one-sided 68 / 95 % upper limits and its prior "
        "median (a central interval excludes a truth on the e = 0 boundary by construction). "
        "ln Z errors are dynesty's quoted errors × err_inflation "
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
