"""M8 A8b–A8f: the pre-registered runs (docs/notes/M8.md, "Pre-registration: A8b–A8f").

    uv run python -m fork.export.m8_runs laplace-check [--n-workers 12]  # pre-launch check
    uv run python -m fork.export.m8_runs plan
    uv run python -m fork.export.m8_runs run [--n-workers 12]   # the 48 sampling jobs
    uv run python -m fork.export.m8_runs job JOB_ID [--maxiter N --cache DIR]
    uv run python -m fork.export.m8_runs reweight [--n-workers 12]
    uv run python -m fork.export.m8_runs summarize   # gates, then A8b–A8f

Nothing here is retuned after a result is seen. Every setting below is the pre-registration's:

- **Fast likelihood:** ``PerturbativeProblem`` (a2b342a), ``second_order=False``, Δs = 60 d,
  cubic quadrature, u_ref from ``docs/notes/M8_u_ref.json`` (pinned by SHA-256 here).
- **Sampler:** M6c's: dynesty static, ``multi``/``rslice`` with dynesty's default 3 + ndim = 10
  slices, dlogz 0.1, nlive 2000 in cell 0 and 1000 elsewhere, over φ only (7 dimensions; the
  prior transform is the M6c cell's, restricted to φ, with the same periodic dimensions).
- **Jobs:** headline, S1, S2 × 4 cells × seeds 1–4; rng ``[2026, 8, 1846, v, cell + 1, s]``.
  Each job runs in one process (no pool: a 0.2 ms likelihood gains nothing from one) and
  checkpoints every 300 s.
- **Reweighting:** 6000 equal-weight draws per seed run and per 4-seed merge, each evaluated
  with the fast likelihood and with mode (b) (``LinearizedProblem``), in resumable chunks.

A saved run is a ``NestedResult`` whose ``space`` is the M6c cell's (11-entry q) and whose
samples are φ, the last 7 entries of q.

**Laplace check** (``laplace-check``, before launch; the user's, 2026-11-09): at 12 points, the
exact integral over u, ∫ exp(−½‖z(u, φ)‖²) du / V_box over the box, by importance sampling
from N(u*, 2²·Σ) with mode (b)'s u* and Σ, against mode (b)'s −2 ln L. At least 10 000 draws,
continued until the Monte Carlo error of −2 ln L is below 0.01. Gate: |difference| < 0.05 at
every point.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from fork.export import m8
from fork.export.to_json import write_json

ROOT = m8.ROOT
CACHE = m8.CACHE / "runs"
U_REF_FILE = ROOT / "docs" / "notes" / "M8_u_ref.json"
U_REF_SHA256 = "4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5"
VARIANTS = ("headline", "s1_e05", "s2_sigma_data")

# ---------------------------------------------------------------------------- Laplace check

LAPLACE_RNG = (8, 1846, 4)
LAPLACE_BATCH = 2000
LAPLACE_MIN_DRAWS = 10_000
LAPLACE_MAX_DRAWS = 600_000
LAPLACE_TARGET_ERR = 0.01  # MC error of −2 ln L
LAPLACE_GATE = 0.05
PROPOSAL_SCALE = 2.0  # proposal σ = 2 × mode (b)'s


def _log(msg: str, cache: Path = CACHE) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    with (cache / "progress.log").open("a") as fh:
        fh.write(line + "\n")


def laplace_points() -> list[dict]:
    """12 points (deterministic): 4 R2 core, 4 R2 ridge (2 inner, 2 outer branch), 2 headline
    cell 0, 2 S2; from the A8a and per-cell draws (full-posterior draws)."""
    from fork.export import m8_perturbative as mp

    rng = np.random.default_rng(list(LAPLACE_RNG))
    a8a, cells = mp.load_a8a(), mp.load_cells()

    def pick(recs, n, tag):
        idx = sorted(rng.choice(len(recs), n, replace=False).tolist())
        return [{"tag": tag, "source": recs[i]["source"], "variant": recs[i]["variant"],
                 "q": recs[i]["q"], "a_au": recs[i]["a_au"], "e": recs[i]["e"],
                 "m2lnl_laplace": recs[i]["m2lnl_laplace"]} for i in idx]  # fmt: skip

    r2 = [r for r in a8a if r["source"] == "r2"]
    inner = [r for r in r2 if r["class"] == "ridge" and r["a_au"] < mp.INNER_OUTER_A]
    outer = [r for r in r2 if r["class"] == "ridge" and r["a_au"] >= mp.INNER_OUTER_A]
    return (pick([r for r in r2 if r["class"] == "core"], 4, "r2_core")
            + pick(inner, 2, "ridge_inner") + pick(outer, 2, "ridge_outer")
            + pick([r for r in cells if r["source"] == "h_c0"], 2, "headline_cell0")
            + pick([r for r in a8a if r["source"] == "s2"], 2, "s2"))  # fmt: skip


def _laplace_one(args: tuple[int, dict]) -> dict:
    """Importance-sampling integral over u at one point, checkpointed per batch."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from scipy.stats import multivariate_normal

    k, p = args
    _, _, space, lp = m8.problem(p["variant"])
    phi = np.array(p["q"])[4:]
    t0 = time.process_time()
    b = lp.marginal(phi, "b")
    z = lp.z_fn(phi)
    box = space.uranus_box
    cov = PROPOSAL_SCALE**2 * b.cov_u
    prop = multivariate_normal(b.u_star, cov)
    rng = np.random.default_rng([*LAPLACE_RNG, k])
    ckpt = CACHE / "laplace" / f"point_{k:02d}.npz"
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    log_w = np.load(ckpt)["log_w"] if ckpt.exists() else np.empty(0)
    if log_w.size:  # advance the generator past the batches already done
        for _ in range(log_w.size // LAPLACE_BATCH):
            rng.multivariate_normal(b.u_star, cov, LAPLACE_BATCH)

    def estimate(lw):
        n = lw.size
        ln_i = float(logsumexp(lw) - math.log(n))
        e = float(math.exp(2.0 * logsumexp(lw) - logsumexp(2.0 * lw)))
        err = 2.0 * math.sqrt(max(n / e - 1.0, 0.0) / n)  # of −2 ln I (delta method)
        return -2.0 * ln_i, err, e

    while True:
        if log_w.size >= LAPLACE_MIN_DRAWS:
            m2, err, e = estimate(log_w)
            if err < LAPLACE_TARGET_ERR or log_w.size >= LAPLACE_MAX_DRAWS:
                break
        u = rng.multivariate_normal(b.u_star, cov, LAPLACE_BATCH)
        inside = np.all((u >= box.lo) & (u <= box.hi), axis=1)
        lw = np.full(LAPLACE_BATCH, -np.inf)
        for i in np.flatnonzero(inside):
            lw[i] = -0.5 * float(np.sum(z(u[i]) ** 2))
        lw -= prop.logpdf(u) + box.log_volume
        log_w = np.concatenate([log_w, lw])
        tmp = ckpt.with_name(ckpt.stem + ".tmp.npz")
        np.savez(tmp, log_w=log_w)
        tmp.replace(ckpt)
    m2, err, e = estimate(log_w)
    return {"k": k, "tag": p["tag"], "source": p["source"], "variant": p["variant"],
            "a_au": p["a_au"], "e": p["e"], "m2lnl_b": b.m2lnl,
            "m2lnl_laplace": p["m2lnl_laplace"], "m2lnl_is": m2, "mc_err": err,
            "diff_is_minus_b": m2 - b.m2lnl, "n_draws": int(log_w.size), "ess": e,
            "ess_frac": e / log_w.size, "n_outside_box": int(np.sum(~np.isfinite(log_w))),
            "converged": err < LAPLACE_TARGET_ERR,
            "cpu_s_this_session": time.process_time() - t0}  # fmt: skip


def laplace_check(n_workers: int) -> dict:
    from concurrent.futures import ProcessPoolExecutor

    from fork import ledger

    pts = laplace_points()
    _log(f"laplace-check: {len(pts)} points, proposal N(u*, {PROPOSAL_SCALE}² Σ_b), "
         f"≥{LAPLACE_MIN_DRAWS} draws until MC err < {LAPLACE_TARGET_ERR}, "
         f"{n_workers} workers")  # fmt: skip
    t0, c0 = time.time(), ledger.cpu_now()
    out = []
    with ProcessPoolExecutor(min(n_workers, len(pts))) as ex:
        for r in ex.map(_laplace_one, list(enumerate(pts))):
            out.append(r)
            _log(f"laplace-check {r['k']:2d} {r['tag']:15s} IS − b = {r['diff_is_minus_b']:+.4f} "
                 f"± {r['mc_err']:.4f}, N {r['n_draws']}, ESS/N {r['ess_frac']:.3f}")  # fmt: skip
    worst = max(abs(r["diff_is_minus_b"]) for r in out)
    res = {"points": out, "max_abs_diff": worst, "gate": LAPLACE_GATE,
           "passed": bool(worst < LAPLACE_GATE and all(r["converged"] for r in out)),
           "settings": {"rng": list(LAPLACE_RNG), "proposal_scale": PROPOSAL_SCALE,
                        "batch": LAPLACE_BATCH, "min_draws": LAPLACE_MIN_DRAWS,
                        "max_draws": LAPLACE_MAX_DRAWS,
                        "target_err": LAPLACE_TARGET_ERR}}  # fmt: skip
    (CACHE / "laplace_check.json").write_text(json.dumps(res, indent=1))
    n = sum(r["n_draws"] for r in out)
    ledger.append("m8/laplace_check", "diagnostic", n_likelihood_calls=n, n_integrations=n,
                  sim_years=None, wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0,
                  n_workers=n_workers, n_points=len(out))  # fmt: skip
    _log(f"laplace-check: max |IS − b| = {worst:.4f}, gate {LAPLACE_GATE}: "
         f"{'PASS' if res['passed'] else 'FAIL'}")  # fmt: skip
    return res


# ---------------------------------------------------------------------------- frozen inputs

# The pre-registration's code, pinned: the fast likelihood and the reweighting module as at
# a2b342a, mode (b) as at 5b7c80a. Checked before every job, as m6c.check_frozen is.
FROZEN_CODE = {
    "src/fork/model/perturbative.py":
        "092f673b9352f9d2b0ae709f2a4d3f258dd09806245b767fd83b86e8ecccd68c",
    "src/fork/model/linearized.py":
        "3e02e78bdd9ee74fb96a6666f192188c355113addc92f43e376729f0b616d15e",
    "src/fork/inference/reweight.py":
        "d0c5fbff6066ebaa531bc24668fa807a8f7bf0aa96cf77eead5c07aca9cd021b",
}  # fmt: skip
SECOND_ORDER = False  # the switch rule's outcome (docs/notes/M8.md)
STEP_DAYS = 60.0
QUADRATURE = "cubic"
SEEDS = (1, 2, 3, 4)
N_PHI = 7
CHECKPOINT_EVERY_S = 300.0
PROGRESS_EVERY_S = 60.0
REWEIGHT_CHUNK = 250
ESS_GATE_CELL = 0.3
TRUNC_GATE_M2LNL = 1e-6


def check_frozen(variant: str) -> None:
    from fork.export import m6c

    m6c.check_frozen(variant)
    for rel, digest in FROZEN_CODE.items():
        if m6c.sha256(ROOT / rel) != digest:
            raise RuntimeError(f"{rel} differs from the pre-registered code")
    if m6c.sha256(U_REF_FILE) != U_REF_SHA256:
        raise RuntimeError(f"{U_REF_FILE} differs from the pre-registered references")


def u_ref(variant: str) -> np.ndarray:
    return np.array(json.loads(U_REF_FILE.read_text())[variant]["u_ref"], dtype=float)


def fast_problem(variant: str):
    """(obs, ctx, uncut space, LinearizedProblem, PerturbativeProblem) for a data set."""
    from fork.model.perturbative import PerturbativeProblem

    obs, ctx, space, lp = m8.problem(variant)
    pr = PerturbativeProblem.build(obs, ctx, space, u_ref(variant), step_days=STEP_DAYS,
                                   quadrature=QUADRATURE)  # fmt: skip
    return obs, ctx, space, lp, pr


# ---------------------------------------------------------------------------- jobs


class RunJob:
    """One nested run: data set × cell × seed."""

    def __init__(self, variant: str, cell: int, seed: int) -> None:
        self.variant, self.cell, self.seed = variant, cell, seed

    @property
    def id(self) -> str:
        return f"{self.variant}/c{self.cell}/s{self.seed}"

    @staticmethod
    def from_id(job_id: str) -> RunJob:
        v, c, s = job_id.split("/")
        return RunJob(v, int(c.removeprefix("c")), int(s.removeprefix("s")))

    @property
    def nlive(self) -> int:
        from fork.export import m6c

        return m6c.NLIVE_CELL.get(self.cell, m6c.NLIVE_DEFAULT)

    def rng_seed(self) -> list[int]:
        from fork.export import m6c

        return [2026, 8, 1846, m6c.VARIANTS.index(self.variant), self.cell + 1, self.seed]

    def rng(self) -> np.random.Generator:
        return np.random.default_rng(self.rng_seed())

    def path(self, suffix: str, cache: Path = CACHE) -> Path:
        return cache / "jobs" / (self.id.replace("/", "__") + suffix)


def n_cells(variant: str) -> int:
    from fork.export import m6c

    return len(m6c.variant_cells(variant))


def all_jobs() -> list[RunJob]:
    return [RunJob(v, k, s) for v in VARIANTS for k in range(n_cells(v)) for s in SEEDS]


_S: dict = {}
_U_HALF = np.full(4, 0.5)


def _fast_logl(phi: np.ndarray) -> float:
    return _S["pr"].log_likelihood(phi, second_order=SECOND_ORDER)


def _phi_transform(u: np.ndarray) -> np.ndarray:
    """The M6c cell's prior transform on φ (the Uranus block is unused: u is marginalised)."""
    return _S["space"].prior_transform(np.concatenate([_U_HALF, u]))[4:]


def do_job(job: RunJob, maxiter: int | None = None, cache: Path = CACHE) -> None:
    import dynesty

    from fork import ledger
    from fork.export import m6c
    from fork.inference.nested import BOUND, DLOGZ, SAMPLE, NestedResult, _ProgressLog, cell_space

    if maxiter is not None and cache == CACHE:
        raise SystemExit("smoke runs (--maxiter) go to a scratch --cache, never the run cache")
    check_frozen(job.variant)
    out = job.path(".pkl", cache)
    if out.exists():
        print(f"{job.id}: cached")
        return
    *_, space, _, pr = fast_problem(job.variant)
    cs = cell_space(space, m6c.variant_cells(job.variant)[job.cell])
    _S.update(pr=pr, space=cs)
    periodic = [i - 4 for i in cs.periodic_unit_dims]
    if min(periodic) < 0 or len(cs.unit_names) - 4 != N_PHI:
        raise RuntimeError("unexpected prior layout")
    ckpt = job.path(".ckpt", cache)
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    log = _ProgressLog(cache / "progress.log", PROGRESS_EVERY_S, f"m8/{job.id}")
    run_kw = dict(dlogz=DLOGZ, maxiter=maxiter, print_progress=True, print_func=log,
                  checkpoint_file=str(ckpt), checkpoint_every=CHECKPOINT_EVERY_S)  # fmt: skip
    resumed = ckpt.exists()
    t0, c0 = time.time(), ledger.cpu_now()
    if resumed:
        sampler = dynesty.NestedSampler.restore(str(ckpt))
        log.line(f"resumed from {ckpt.name} at iteration {sampler.it}")
        sampler.run_nested(resume=True, **run_kw)
    else:
        log.line(f"start: ndim {N_PHI}, nlive {job.nlive}, {BOUND}/{SAMPLE} default slices, "
                 f"dlogz {DLOGZ}, rng {job.rng_seed()}, second_order {SECOND_ORDER}")  # fmt: skip
        sampler = dynesty.NestedSampler(_fast_logl, _phi_transform, N_PHI, nlive=job.nlive,
                                        bound=BOUND, sample=SAMPLE, periodic=periodic,
                                        rstate=job.rng())  # fmt: skip
        sampler.run_nested(**run_kw)
    res = NestedResult(cs, sampler.results, job.nlive, DLOGZ, time.time() - t0)
    res.save(out)
    ckpt.unlink(missing_ok=True)
    log.line(f"done: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, {res.ncall} calls")
    if cache == CACHE:
        ledger.append(f"m8/{job.id}", "nested", n_likelihood_calls=res.ncall, n_integrations=0,
                      sim_years=0.0, wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0,
                      n_workers=1, estimated=resumed, resumed=resumed, nlive=job.nlive,
                      niter=res.niter, likelihood="perturbative, second order off")  # fmt: skip
    print(f"{job.id}: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, "
          f"{res.ncall} calls, {(time.time() - t0) / 60:.1f} min")  # fmt: skip


def run_all(n_slots: int, cache: Path = CACHE, maxiter: int | None = None) -> bool:
    """Scheduler: each job a one-process subprocess, at most ``n_slots`` at once; cell 0
    (nlive 2000) first; one retry per failed job (it resumes from its checkpoint)."""
    import subprocess
    import sys

    for v in VARIANTS:
        check_frozen(v)
    pending = sorted([j for j in all_jobs() if not j.path(".pkl", cache).exists()],
                     key=lambda j: (-j.nlive, j.variant, j.cell, j.seed))  # fmt: skip
    _log(f"run: {len(pending)} jobs pending, {n_slots} slots"
         + (f", smoke maxiter {maxiter}" if maxiter else ""), cache)  # fmt: skip
    running: dict[str, tuple] = {}
    failures: dict[str, int] = {}
    while pending or running:
        for jid, (proc, t0) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[jid]
            job = RunJob.from_id(jid)
            if rc == 0 and job.path(".pkl", cache).exists():
                _log(f"done {jid} ({(time.time() - t0) / 60:.1f} min)", cache)
            else:
                failures[jid] = failures.get(jid, 0) + 1
                _log(f"FAILED {jid} (exit {rc}, attempt {failures[jid]})", cache)
                if failures[jid] < 2:
                    pending.insert(0, job)
        while pending and len(running) < n_slots:
            job = pending.pop(0)
            log = job.path(".log", cache)
            log.parent.mkdir(parents=True, exist_ok=True)
            cmd = [sys.executable, "-m", "fork.export.m8_runs", "job", job.id,
                   "--cache", str(cache)]  # fmt: skip
            if maxiter:
                cmd += ["--maxiter", str(maxiter)]
            with log.open("a") as fh:
                proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT,
                                        env=os.environ | {"OMP_NUM_THREADS": "1"})  # fmt: skip
            running[job.id] = (proc, time.time())
            _log(f"start {job.id} (nlive {job.nlive}, pid {proc.pid})", cache)
        time.sleep(5 if maxiter else 20)
    done = all(j.path(".pkl", cache).exists() for j in all_jobs())
    _log(f"run finished; failures: {failures or 'none'}; all done: {done}", cache)
    return done


def plan(cache: Path = CACHE) -> None:
    jobs = all_jobs()
    for j in jobs:
        state = ("done" if j.path(".pkl", cache).exists() else
                 "resumable" if j.path(".ckpt", cache).exists() else "pending")  # fmt: skip
        print(f"{j.id:24s} nlive {j.nlive:5d}  rng {j.rng_seed()}  {state}")
    print(f"{len(jobs)} jobs, {sum(j.path('.pkl', cache).exists() for j in jobs)} done")


# ---------------------------------------------------------------------------- reweighting

MERGED = 0  # the seed index of a cell's 4-seed merge


def posteriors() -> list[tuple[str, int, int]]:
    """(variant, cell, seed) of every cell posterior to reweight: seeds 1–4, then the merge."""
    return [(v, k, s) for v in VARIANTS for k in range(n_cells(v)) for s in (*SEEDS, MERGED)]


def pid(v: str, k: int, s: int) -> str:
    return f"{v}__c{k}__" + ("merged" if s == MERGED else f"s{s}")


def load_posterior(v: str, k: int, s: int, cache: Path = CACHE):
    from dynesty.utils import merge_runs

    from fork.inference.nested import NestedResult

    if s != MERGED:
        return NestedResult.load(RunJob(v, k, s).path(".pkl", cache))
    runs = [NestedResult.load(RunJob(v, k, t).path(".pkl", cache)) for t in SEEDS]
    merged = merge_runs([r.results for r in runs], print_progress=False)
    return NestedResult(runs[0].space, merged, sum(r.nlive for r in runs), runs[0].dlogz,
                        sum(r.wall_s for r in runs))  # fmt: skip


def _vi(v: str) -> int:
    from fork.export import m6c

    return m6c.VARIANTS.index(v)


def draws(v: str, k: int, s: int, cache: Path = CACHE, max_draws: int | None = None) -> dict:
    """The posterior's equal-weight draws (6000, or all if fewer), cached: φ, ln Z_fast and
    its sampler error. rng [2026, 8, 1846, v, cell + 1, seed (0: merge), 7]."""
    from fork.inference.reweight import N_PER_CELL

    path = cache / "reweight" / pid(v, k, s) / "draws.npz"
    if path.exists():
        d = np.load(path)
        return {key: d[key] for key in d.files}
    res = load_posterior(v, k, s, cache)
    rng = np.random.default_rng([2026, 8, 1846, _vi(v), k + 1, s, 7])
    eq = res.equal_weight(rng)
    n = min(N_PER_CELL if max_draws is None else max_draws, len(eq))
    phi = eq[rng.choice(len(eq), n, replace=False)]
    d = {"phi": phi, "ln_z_fast": np.array(res.logz), "ln_z_fast_err": np.array(res.logzerr),
         "ncall": np.array(res.ncall), "niter": np.array(res.niter)}  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("draws.tmp.npz")
    np.savez(tmp, **d)
    tmp.replace(path)
    return d


_RW: dict = {}


def _init_reweight() -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from fork.inference.profile import GalleGeometry

    for v in VARIANTS:
        _, ctx, space, lp, pr = fast_problem(v)
        _RW[v] = (ctx, space, lp, pr)
    _RW["geo"] = GalleGeometry.from_horizons()


def _reweight_chunk(task: tuple) -> tuple[str, int, float]:
    """Both likelihoods at draws [lo, hi) of one posterior, plus each draw's u from mode (b)'s
    N(u*, Σ) and its derived longitudes (Galle's night, 1847.0) through the full model."""
    from fork.export.m5b import _heliocentric_lon_of_date_deg
    from fork.inference.profile import galle_lon_of_date_rad
    from fork.model.hidden_mass import JD_1847_JAN1_TDB

    v, k, s, c, lo, hi, cache = task
    out = Path(cache) / "reweight" / pid(v, k, s) / f"chunk_{c:03d}.npz"
    if out.exists():
        return pid(v, k, s), c, 0.0
    t_start = time.process_time()
    ctx, space, lp, pr = _RW[v]
    geo = _RW["geo"]
    box = space.uranus_box
    phi = np.load(out.parent / "draws.npz")["phi"][lo:hi]
    n = len(phi)
    rng = np.random.default_rng([2026, 8, 1846, _vi(v), k + 1, s, 9, c])
    a = {"lnl_fast": np.full(n, -np.inf), "fast_outside": np.zeros(n, bool),
         "fast_trunc_log10": np.full(n, np.nan), "lnl_b": np.full(n, -np.inf),
         "b_u_star": np.full((n, 4), np.nan), "b_cov": np.full((n, 4, 4), np.nan),
         "b_outside": np.zeros(n, bool), "b_trunc_m2lnl": np.full(n, np.nan),
         "u_draw": np.full((n, 4), np.nan), "galle_deg": np.full(n, np.nan),
         "helio_1847_deg": np.full(n, np.nan)}  # fmt: skip
    cpu_b = 0.0
    for i, f in enumerate(phi):
        a["lnl_fast"][i] = pr.log_likelihood(f, second_order=SECOND_ORDER)
        try:
            r = pr.evaluate(f, second_order=SECOND_ORDER)
            a["fast_outside"][i], a["fast_trunc_log10"][i] = r.outside_box, r.truncation_log10
        except (ValueError, np.linalg.LinAlgError, ZeroDivisionError):
            pass
        t0 = time.process_time()
        try:
            m = lp.marginal(f, "b")
        except np.linalg.LinAlgError:
            cpu_b += time.process_time() - t0
            continue
        cpu_b += time.process_time() - t0
        if not math.isfinite(m.m2lnl):
            continue
        a["lnl_b"][i] = -0.5 * m.m2lnl
        a["b_u_star"][i], a["b_cov"][i] = m.u_star, m.cov_u
        a["b_outside"][i] = bool(np.any(m.u_star < box.lo) or np.any(m.u_star > box.hi))
        a["b_trunc_m2lnl"][i] = lp.truncation(m)[1]
        u = rng.multivariate_normal(m.u_star, m.cov_u)
        th = space.theta(np.concatenate([u, f]))
        a["u_draw"][i] = u
        a["galle_deg"][i] = math.degrees(galle_lon_of_date_rad(th, ctx, geo)) % 360.0
        a["helio_1847_deg"][i] = _heliocentric_lon_of_date_deg(th, ctx, JD_1847_JAN1_TDB)
    a["cpu_b_s"] = np.array(cpu_b)
    a["cpu_s"] = np.array(time.process_time() - t_start)
    tmp = out.with_name(out.stem + ".tmp.npz")
    np.savez(tmp, **a)
    tmp.replace(out)
    return pid(v, k, s), c, float(a["cpu_s"])


def reweight(n_workers: int, cache: Path = CACHE, max_draws: int | None = None) -> None:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    from fork import ledger

    for v in VARIANTS:
        check_frozen(v)
    missing = [j.id for j in all_jobs() if not j.path(".pkl", cache).exists()]
    if missing:
        raise SystemExit(f"runs missing: {missing}")
    t0, c0 = time.time(), ledger.cpu_now()
    tasks = []
    for v, k, s in posteriors():
        n = len(draws(v, k, s, cache, max_draws)["phi"])
        for c, lo in enumerate(range(0, n, REWEIGHT_CHUNK)):
            if not (cache / "reweight" / pid(v, k, s) / f"chunk_{c:03d}.npz").exists():
                tasks.append((v, k, s, c, lo, min(lo + REWEIGHT_CHUNK, n), str(cache)))
    _log(f"reweight: {len(posteriors())} posteriors, {len(tasks)} chunks to do, "
         f"{n_workers} workers", cache)  # fmt: skip
    cpu_chunks = 0.0
    with ProcessPoolExecutor(n_workers, initializer=_init_reweight) as ex:
        futs = [ex.submit(_reweight_chunk, t) for t in tasks]
        for i, f in enumerate(as_completed(futs)):
            p, c, cpu = f.result()
            cpu_chunks += cpu
            if (i + 1) % 24 == 0 or i + 1 == len(futs):
                _log(f"reweight: {i + 1}/{len(futs)} chunks (last {p} #{c})", cache)
    if cache == CACHE and tasks:
        n_draws = sum(t[5] - t[4] for t in tasks)
        ledger.append("m8/reweight", "reweight", n_likelihood_calls=n_draws, n_integrations=None,
                      sim_years=None, wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0,
                      n_workers=n_workers, cpu_chunks_s=round(cpu_chunks, 1),
                      note="mode (b) + fast + derived longitudes per draw")  # fmt: skip


def load_reweighted(v: str, k: int, s: int, cache: Path = CACHE) -> dict:
    d = draws(v, k, s, cache)
    n = len(d["phi"])
    parts = [np.load(cache / "reweight" / pid(v, k, s) / f"chunk_{c:03d}.npz")
             for c in range(math.ceil(n / REWEIGHT_CHUNK))]  # fmt: skip
    out = dict(d)
    for key in parts[0].files:
        if parts[0][key].ndim == 0:
            out[key] = float(sum(float(p[key]) for p in parts))
        else:
            out[key] = np.concatenate([p[key] for p in parts])
    return out


# ---------------------------------------------------------------------------- summary

OUT_NAME = "m8_runs.json"
INNER_OUTER_A = 30.5
E_RIDGE = 0.06
# Published comparison values, as written in the pre-registration (docs/notes/M8.md).
R2_CELL2 = {"median": 325.86, "width68": 4.06, "ln_z": -360.21}
REPLACEMENT = {"median": 325.82, "lo68": 323.57, "hi68": 327.79,
               "ln_z_k": (-364.87, -365.18, -359.77, -364.87)}  # fmt: skip
PUBLISHED = {"s1_e05": {"median": 326.21, "width68": 5.89, "ln_b": -0.53},
             "s2_sigma_data": {"median": 325.78, "width68": 4.41, "ln_b": -4.75}}  # fmt: skip
TOL = {"a8b_median": 0.35, "a8b_width_rel": 0.10, "a8b_ln_z": 0.5, "a8c_median": 0.35,
       "a8c_bound": 0.4, "a8c_ln_z_k": 1.5, "a8d_deg": 1.3, "a8d_ln_b": 2.2,
       "a8e_cpu_h": 2.85, "a8e_full_cpu_h": 28.5, "a8f_width_rel": 0.10,
       "rule_95_std": 0.4}  # fmt: skip
FIT_RNG_TAG = 11
CELL_RNG_TAG = 13


def _perturber(space, phi: np.ndarray) -> np.ndarray:
    """(mass M⊕, a AU, e) per draw: from φ alone (u does not enter them)."""
    th = np.array([space.theta(np.concatenate([space.uranus_box.centre, f])) for f in phi])
    return np.c_[10.0 ** th[:, 6], th[:, 7], th[:, 8]]


def _resampled(cells: list, log_c: list[float], rng: np.random.Generator, n: int):
    """Combine reweighted cells (``reweight.combine``) and draw n systematic-resampled indices
    into the stacked draws. Returns (Reweighted, indices)."""
    from fork.inference.nested import ERR_INFLATION
    from fork.inference.reweight import combine

    offs = np.cumsum([0] + [c.n for c in cells])
    idx = [np.arange(offs[i], offs[i + 1])[:, None] for i in range(len(cells))]
    res = combine(cells, idx, log_c, ERR_INFLATION)
    return res, res.equal_weight(rng, n)[:, 0].astype(int)


def _lon(galle: np.ndarray, targets: dict) -> dict:
    from fork.export import m6c

    return m6c.lon_summary(galle, targets)


def _wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / np.sum(w)
    return [float(x[o][min(np.searchsorted(c, p), len(x) - 1)]) for p in q]


def summarize(cache: Path = CACHE) -> dict:
    from fork import ledger
    from fork.export import m6c, m6c_repeat
    from fork.export import m8_perturbative as mp
    from fork.inference.nested import ERR_INFLATION
    from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
    from fork.inference.reweight import N_PER_CELL, ess, reweight_cell

    for v in VARIANTS:
        check_frozen(v)
    t0, c0 = time.time(), ledger.cpu_now()
    tg = m6c.comparison_targets()
    targets = {"jpl": tg["jpl"]["galle_lon_ecl_of_date_deg"],
               "le_verrier_derived":
                   tg["le_verrier"]["galle_lon_ecl_of_date_deg_derived"]}  # fmt: skip
    space_of = {v: m8.problem(v)[2] for v in VARIANTS}

    # per-cell posteriors
    post: dict = {}
    for v, k, s in posteriors():
        d = load_reweighted(v, k, s, cache)
        cr = reweight_cell(d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]),
                           float(d["ln_z_fast_err"]))  # fmt: skip
        d["pert"] = _perturber(space_of[v], d["phi"])
        d["cr"] = cr
        post[(v, k, s)] = d
    cells_out, gate2, gate3, gate4 = [], [], [], []
    for (v, k, s), d in post.items():
        cr = d["cr"]
        rng = np.random.default_rng([2026, 8, 1846, _vi(v), k + 1, s, CELL_RNG_TAG])
        _, ix = _resampled([cr], [0.0], rng, cr.n)
        g = _lon(d["galle_deg"][ix], targets)
        rec = {"variant": v, "cell": k, "seed": s, "n": cr.n, "ess_frac": cr.ess / cr.n,
               "ln_z_fast": cr.ln_z_fast, "ln_z_fast_err": cr.ln_z_fast_err,
               "ln_correction": cr.ln_correction, "ln_correction_err": cr.ln_correction_err,
               "ln_z": cr.ln_z, "galle_median": g["median"], "galle_width68": g["width68"],
               "ncall": int(d["ncall"]), "niter": int(d["niter"]),
               "fast_outside_box": int(d["fast_outside"].sum()),
               "b_outside_box": int(d["b_outside"].sum()),
               "b_trunc_m2lnl_max": float(np.nanmax(d["b_trunc_m2lnl"])),
               "n_b_failed": int(np.sum(~np.isfinite(d["lnl_b"])))}  # fmt: skip
        cells_out.append(rec)
        gate2.append({"posterior": pid(v, k, s), "ess_frac": rec["ess_frac"],
                      "pass": rec["ess_frac"] >= ESS_GATE_CELL})  # fmt: skip
        gate4.append({"posterior": pid(v, k, s), "max": rec["b_trunc_m2lnl_max"],
                      "pass": bool(np.all(d["b_trunc_m2lnl"][np.isfinite(d["lnl_b"])]
                                          < TRUNC_GATE_M2LNL))})  # fmt: skip
        if k == 2:
            a, e = d["pert"][:, 1], d["pert"][:, 2]
            for name, m in (("inner", (e >= E_RIDGE) & (a < INNER_OUTER_A)),
                            ("outer", (e >= E_RIDGE) & (a >= INNER_OUTER_A))):  # fmt: skip
                nb = int(m.sum())
                fr = ess(cr.delta[m]) / nb if nb else math.nan
                gate3.append({"posterior": pid(v, k, s), "branch": name, "n": nb,
                              "ess_frac": fr, "pass": bool(nb > 0 and fr >= 0.3)})  # fmt: skip

    # fits: the published one (merged cells) and the four single-seed fits, per data set
    fits, gate1 = {}, []
    for v in VARIANTS:
        for s in (MERGED, *SEEDS):
            cells = [post[(v, k, s)] for k in range(n_cells(v))]
            log_c = [math.log(_cell_space(v, k).slice_weight) for k in range(n_cells(v))]
            rng = np.random.default_rng([2026, 8, 1846, _vi(v), 0, s, FIT_RNG_TAG])
            res, ix = _resampled([c["cr"] for c in cells], log_c, rng, N_PER_CELL)
            galle = np.concatenate([c["galle_deg"] for c in cells])[ix]
            helio = np.concatenate([c["helio_1847_deg"] for c in cells])[ix]
            pert = np.concatenate([c["pert"] for c in cells])[ix]
            shares = np.exp(np.array(log_c) + np.array([c["cr"].ln_z for c in cells]) - res.ln_z)
            f = {"variant": v, "seed": s, "ln_z": res.ln_z, "ln_z_err": res.ln_z_err,
                 "ess_frac": res.ess_frac, "shares": shares.tolist(),
                 "galle_lon_ecl_of_date_deg": _lon(galle, targets),
                 "helio_lon_ecl_of_date_1847_deg": _lon(helio, {
                     "jpl": tg["jpl"]["helio_lon_ecl_of_date_1847_deg"],
                     "le_verrier": tg["le_verrier"]["helio_lon_ecl_of_date_1847_deg"]}),
                 "mass_earth": _q(pert[:, 0]), "a_au": _q(pert[:, 1]), "e": _q(pert[:, 2]),
                 "p_a_below_25_au": float(np.mean(pert[:, 1] < 25.0))}  # fmt: skip
            fits[f"{v}/{'merged' if s == MERGED else f's{s}'}"] = f
            gate1.append({"fit": f"{v}/{'merged' if s == MERGED else f's{s}'}",
                          "ess_frac": res.ess_frac, "pass": res.ess_frac >= 0.5})  # fmt: skip

    def gates_ok(v: str) -> bool:
        return all(g["pass"] for g in gate1 if g["fit"].startswith(v + "/")) and all(
            g["pass"] for g in gate2 + gate3 + gate4 if g["posterior"].startswith(v + "__")
        )

    # blind-spot check 3: core and branch masses in cell 2, fast / reweighted / R2
    d2 = post[("headline", 2, MERGED)]
    a, e = d2["pert"][:, 1], d2["pert"][:, 2]
    masks = {"core": e < E_RIDGE, "inner": (e >= E_RIDGE) & (a < INNER_OUTER_A),
             "outer": (e >= E_RIDGE) & (a >= INNER_OUTER_A)}  # fmt: skip
    w = np.exp(d2["cr"].delta - d2["cr"].delta.max())
    r2 = m6c_repeat.replacement_cells()[2]
    th_r2 = np.array([r2.space.theta(q) for q in r2.samples])
    w_r2 = r2.weights
    r2_masks = {"core": th_r2[:, 8] < E_RIDGE,
                "inner": (th_r2[:, 8] >= E_RIDGE) & (th_r2[:, 7] < INNER_OUTER_A),
                "outer": (th_r2[:, 8] >= E_RIDGE) & (th_r2[:, 7] >= INNER_OUTER_A)}  # fmt: skip
    check3 = {name: {"fast": float(m.mean()), "reweighted": float(w[m].sum() / w.sum()),
                     "r2": float(w_r2[r2_masks[name]].sum() / w_r2.sum())}
              for name, m in masks.items()}  # fmt: skip
    pert_json = json.loads((ROOT / "site" / "data" / mp.OUT_NAME).read_text())
    check1 = pert_json["a8a"]["blind_spot_candidates"]

    # tests
    head_ok = gates_ok("headline")
    cell = {(r["variant"], r["cell"], r["seed"]): r for r in cells_out}
    c2 = cell[("headline", 2, MERGED)]
    a8b = {"median": c2["galle_median"], "width68": c2["galle_width68"], "ln_z": c2["ln_z"],
           "d_median": c2["galle_median"] - R2_CELL2["median"],
           "width_rel": c2["galle_width68"] / R2_CELL2["width68"] - 1.0,
           "d_ln_z": c2["ln_z"] - R2_CELL2["ln_z"]}  # fmt: skip
    a8b["criteria"] = bool(abs(a8b["d_median"]) <= TOL["a8b_median"]
                           and abs(a8b["width_rel"]) <= TOL["a8b_width_rel"]
                           and abs(a8b["d_ln_z"]) <= TOL["a8b_ln_z"])  # fmt: skip
    a8b["pass"] = a8b["criteria"] and head_ok
    hg = fits["headline/merged"]["galle_lon_ecl_of_date_deg"]
    lnzk = [cell[("headline", k, MERGED)]["ln_z"] for k in range(4)]
    a8c = {"median": hg["median"], "lo68": hg["lo68"], "hi68": hg["hi68"], "ln_z_k": lnzk,
           "d_median": hg["median"] - REPLACEMENT["median"],
           "d_lo68": hg["lo68"] - REPLACEMENT["lo68"],
           "d_hi68": hg["hi68"] - REPLACEMENT["hi68"],
           "d_ln_z_k": [z - r for z, r in zip(lnzk, REPLACEMENT["ln_z_k"],
                                              strict=True)]}  # fmt: skip
    a8c["criteria"] = bool(abs(a8c["d_median"]) <= TOL["a8c_median"]
                           and abs(a8c["d_lo68"]) <= TOL["a8c_bound"]
                           and abs(a8c["d_hi68"]) <= TOL["a8c_bound"]
                           and all(abs(x) <= TOL["a8c_ln_z_k"]
                                   for x in a8c["d_ln_z_k"]))  # fmt: skip
    a8c["pass"] = a8c["criteria"] and head_ok
    a8d = {}
    hz = fits["headline/merged"]["ln_z"]
    for v in ("s1_e05", "s2_sigma_data"):
        fv = fits[f"{v}/merged"]
        g = fv["galle_lon_ecl_of_date_deg"]
        zv, zh = fv["ln_z"], hz
        if v == "s2_sigma_data":
            zv -= 0.5 * m6c.gaussian_norm_m2lnl_for(v)
            zh -= 0.5 * m6c.gaussian_norm_m2lnl_for("headline")
        ln_b = zv - zh
        rec = {"median": g["median"], "width68": g["width68"], "ln_b": ln_b,
               "ln_b_err": math.hypot(fv["ln_z_err"], fits["headline/merged"]["ln_z_err"]),
               "d_median": g["median"] - PUBLISHED[v]["median"],
               "d_width68": g["width68"] - PUBLISHED[v]["width68"],
               "d_ln_b": ln_b - PUBLISHED[v]["ln_b"],
               "jpl_in_68": g["vs"]["jpl"]["in_68"]}  # fmt: skip
        ok = (abs(rec["d_median"]) <= TOL["a8d_deg"] and abs(rec["d_width68"]) <= TOL["a8d_deg"]
              and abs(rec["d_ln_b"]) <= TOL["a8d_ln_b"] and rec["jpl_in_68"])  # fmt: skip
        if v == "s1_e05":
            rec["p_a_below_25_au"] = fv["p_a_below_25_au"]
            ok = ok and rec["p_a_below_25_au"] < 0.05
        else:
            shift = (g["median"] - hg["median"] + 180.0) % 360.0 - 180.0
            half = 0.5 * hg["width68"]
            rec["median_shift_headline_half_widths"] = shift / half if half > 0 else math.nan
            ok = ok and abs(rec["median_shift_headline_half_widths"]) < 0.5
        rec["criteria"] = bool(ok)
        rec["pass"] = bool(ok and gates_ok(v) and head_ok)
        a8d[v] = rec
    # A8e: headline seed 1, sampling CPU (ledger) plus 4 × 6000 mode (b) calls
    rows = [r for r in ledger.read() if r["job"] in {f"m8/headline/c{k}/s1" for k in range(4)}]
    b_cpu = [(post[key]["cpu_b_s"], len(post[key]["phi"])) for key in post]
    per_b = sum(c for c, _ in b_cpu) / sum(n for _, n in b_cpu)
    sampling_h = sum(r["cpu_s"] for r in rows) / 3600.0 if len(rows) == 4 else None
    rew_h = 4 * N_PER_CELL * per_b / 3600.0
    a8e = {"sampling_cpu_h": sampling_h, "reweight_cpu_h": rew_h,
           "mode_b_cpu_s_per_call": per_b,
           "ncall": [int(post[("headline", k, 1)]["ncall"]) for k in range(4)]}  # fmt: skip
    if sampling_h is not None:
        a8e["total_cpu_h"] = sampling_h + rew_h
        a8e["factor_vs_full"] = TOL["a8e_full_cpu_h"] / a8e["total_cpu_h"]
        a8e["pass"] = bool(a8e["total_cpu_h"] <= TOL["a8e_cpu_h"])
    a8f = []
    for v in VARIANTS:
        for k in range(n_cells(v)):
            wd = np.array([cell[(v, k, s)]["galle_width68"] for s in SEEDS])
            dev = np.abs(wd / wd.mean() - 1.0)
            a8f.append({"variant": v, "cell": k, "widths": wd.tolist(),
                        "max_rel_dev": float(dev.max()),
                        "pass": bool(dev.max() <= TOL["a8f_width_rel"]
                                     and gates_ok(v))})  # fmt: skip
    # the 95 % replacement rule
    single = [fits[f"headline/s{s}"]["galle_lon_ecl_of_date_deg"] for s in SEEDS]
    std95 = {b: float(np.std([g[b] for g in single], ddof=1)) for b in ("lo95", "hi95")}
    cond = {"gates_headline_and_single_seed_fits": head_ok,
            "a8b_a8c": a8b["pass"] and a8c["pass"],
            "a8f_headline_cells": all(r["pass"] for r in a8f if r["variant"] == "headline"),
            "std_95_bounds_le_0p4":
                all(x <= TOL["rule_95_std"] for x in std95.values())}  # fmt: skip
    rule95 = {"conditions": cond, "std_95_bounds": std95, "replace": all(cond.values()),
              "m8_95": [hg["lo95"], hg["hi95"]],
              "sampler_error_std_over_2": {b: x / 2.0 for b, x in std95.items()}}  # fmt: skip
    # hybrid: R2's cell 2 with M8's cells 0, 1, 3 (reported either way)
    ctx = m8.problem("headline")[1]
    geo = GalleGeometry.from_horizons()
    rng = np.random.default_rng([2026, 8, 1846, 0, 3, 0, 17])
    eq = r2.equal_weight(rng)
    q2 = eq[rng.choice(len(eq), min(N_PER_CELL, len(eq)), replace=False)]
    g2 = np.array([math.degrees(galle_lon_of_date_rad(r2.space.theta(q), ctx, geo)) % 360.0
                   for q in q2])  # fmt: skip
    hy_cells = [post[("headline", k, MERGED)]["cr"] for k in range(4)]
    hy_cells[2] = reweight_cell(np.zeros(len(q2)), r2.logz, r2.logzerr)
    hy_galle = [post[("headline", k, MERGED)]["galle_deg"] for k in range(4)]
    hy_galle[2] = g2
    log_c = [math.log(_cell_space("headline", k).slice_weight) for k in range(4)]
    res, ix = _resampled(hy_cells, log_c, np.random.default_rng([2026, 8, 1846, 0, 0, 0, 19]),
                         N_PER_CELL)  # fmt: skip
    hybrid = {"ln_z": res.ln_z, "ln_z_err": res.ln_z_err,
              "galle_lon_ecl_of_date_deg": _lon(np.concatenate(hy_galle)[ix], targets)}  # fmt: skip

    lap = cache / "laplace_check.json"
    payload = {
        "description": "M8 A8b-A8f (pre-registered in docs/notes/M8.md): perturbative fast "
        "likelihood (second order off), 48 nested runs (3 data sets x 4 cells x 4 seeds), "
        "reweighted with mode (b). Gates first; a test whose posteriors fail a gate is not "
        "passed. ln Z errors: sampler part x err_inflation, plus the reweighting MC error.",
        "settings": {"second_order": SECOND_ORDER, "step_days": STEP_DAYS,
                     "quadrature": QUADRATURE, "seeds": list(SEEDS), "n_draws": N_PER_CELL,
                     "err_inflation": ERR_INFLATION, "u_ref_sha256": U_REF_SHA256,
                     "frozen_code": FROZEN_CODE},
        "laplace_check": json.loads(lap.read_text()) if lap.exists() else None,
        "gates": {"1_fit_ess": gate1, "2_cell_ess": gate2, "3_ridge_branch_ess": gate3,
                  "4_truncation": gate4,
                  "5_blind_spot": {"check1_a8a_candidates": check1,
                                   "check2": "gate 3", "check3_cell2_masses": check3,
                                   "check4": "A8b"},
                  "all_pass": {v: gates_ok(v) for v in VARIANTS}},
        "tests": {"a8b": a8b, "a8c": a8c, "a8d": a8d, "a8e": a8e, "a8f": a8f},
        "rule_95": rule95, "hybrid": hybrid, "fits": fits, "cells": cells_out,
        "tolerances": TOL,
    }  # fmt: skip
    if cache == CACHE:
        write_json(OUT_NAME, payload)
        ledger.append("m8/summarize", "summarize", n_likelihood_calls=0,
                      n_integrations=len(q2), sim_years=None, wall_s=time.time() - t0,
                      cpu_s=ledger.cpu_now() - c0, n_workers=1)  # fmt: skip
    else:
        (cache / OUT_NAME).write_text(json.dumps(payload, indent=1, default=float))
    _log("summarize: gates " + json.dumps(payload["gates"]["all_pass"]), cache)
    return payload


def _cell_space(v: str, k: int):
    from fork.export import m6c
    from fork.inference.nested import cell_space

    return cell_space(m8.problem(v)[2], m6c.variant_cells(v)[k])


def _q(x: np.ndarray) -> dict:
    p = np.percentile(x, [2.5, 16.0, 50.0, 84.0, 97.5])
    return {"lo95": float(p[0]), "lo68": float(p[1]), "median": float(p[2]),
            "hi68": float(p[3]), "hi95": float(p[4])}  # fmt: skip


def pipeline(
    n_workers: int, cache: Path = CACHE, maxiter: int | None = None, max_draws: int | None = None
) -> None:
    """run → reweight → summarize; stops if any run is missing after the scheduler."""
    import subprocess

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()  # fmt: skip
    dirty = subprocess.run(["git", "status", "--porcelain", "src"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()  # fmt: skip
    _log(f"pipeline: commit {head}{' (src modified!)' if dirty else ''}", cache)
    if not run_all(n_workers, cache, maxiter):
        _log("pipeline: runs incomplete, stopping before reweighting", cache)
        return
    reweight(n_workers, cache, max_draws)
    summarize(cache)
    _log("pipeline: done", cache)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("laplace-check", "plan", "run", "job", "reweight", "summarize", "pipeline"):
        p = sub.add_parser(name)
        p.add_argument("--cache", type=Path, default=CACHE)
        if name in ("laplace-check", "run", "reweight", "pipeline"):
            p.add_argument("--n-workers", type=int, default=12)
        if name in ("job", "run", "pipeline"):
            p.add_argument("--maxiter", type=int, default=None, help="smoke runs only")
        if name in ("reweight", "pipeline"):
            p.add_argument("--max-draws", type=int, default=None, help="smoke runs only")
        if name == "job":
            p.add_argument("job_id")
    a = ap.parse_args(argv)
    cache = a.cache.resolve()
    if getattr(a, "max_draws", None) and cache == CACHE:
        raise SystemExit("--max-draws is for smoke runs in a scratch --cache")
    if a.cmd == "laplace-check":
        laplace_check(a.n_workers)
    elif a.cmd == "plan":
        plan(cache)
    elif a.cmd == "run":
        run_all(a.n_workers, cache, a.maxiter)
    elif a.cmd == "job":
        do_job(RunJob.from_id(a.job_id), a.maxiter, cache)
    elif a.cmd == "reweight":
        reweight(a.n_workers, cache, a.max_draws)
    elif a.cmd == "summarize":
        summarize(cache)
    else:
        pipeline(a.n_workers, cache, a.maxiter, a.max_draws)


if __name__ == "__main__":
    main()
