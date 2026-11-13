"""S1 cell 0 by sub-cells (pre-registered in docs/notes/S1_cell0.md; nothing here may differ).

    uv run python -m fork.export.s1_cell0 plan
    uv run python -m fork.export.s1_cell0 pipeline [--n-workers 12]   # A, reweight, B, summary
    uv run python -m fork.export.s1_cell0 job-a K S [--maxiter N --cache DIR --synthetic SEED]
    uv run python -m fork.export.s1_cell0 job-b K R [--n-workers 6 --maxiter N --cache DIR ...]
    uv run python -m fork.export.s1_cell0 reweight | summarize

- **Stage A:** the fast likelihood of M8 (``m8_runs``: second order off, 60 d, cubic) in each
  of the six sub-cells of ``v13.S1_SUBCELLS``, seeds 1 to 4, nlive 4000, 28 slices, one
  process per run. rng ``[2026, 14, 1846, 1, k + 1, s]``.
- **Reweighting:** 6,000 draws per run and per four-seed merge, the fast likelihood against
  mode (b), as ``m8_runs`` does it.
- **Stage B:** two full-likelihood runs (M6c's 11 dimensions) in the sub-cell with the
  largest w_k Z_k, nlive 4000, 28 slices, 6 workers each. rng ``[..., k + 1, 100 + r]``.
- **Summary:** gates G1 to G3, tests T1 and T2, the evidence of cell 0 and P(a < 25 AU), read
  by the plan's § 5.

``--synthetic SEED`` swaps the data for a synthetic set on Table 77's layout
(``v13.synthetic_dataset``) for smoke runs; it needs a scratch ``--cache``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from fork.export import m6c, m8_runs, v13

VARIANT = "s1_e05"
SUBCELLS = v13.S1_SUBCELLS
SEEDS = (1, 2, 3, 4)
MERGED = 0
B_RUNS = (1, 2)
NLIVE = 4000
SLICES = 28
B_WORKERS = 6
TAG = 14
N_PHI = m8_runs.N_PHI
N_DRAWS = 6000
CHUNK = 250
ROOT = m6c.ROOT
CACHE = ROOT / "data" / "cache" / "nested" / "s1_cell0"
OUT_JSON = ROOT / "docs" / "notes" / "S1_cell0.json"
ESS_GATE = 0.3
TRUNC_GATE = 1e-6
TOL = {"ln_z": 0.5, "median": 0.35, "width_rel": 0.10}
AT_TOL = {"ln_z": 0.02, "median": 0.02, "width_rel": 0.005}
MATERIAL_SHARE = 0.01
OTHER_PUBLISHED = (-366.0354, -359.7935, -365.1071)
OTHER_FAST = (-365.682, -360.290, -365.258)
CHECK_VALUES = {"published": (-366.8432, 0.00104, 2e-5), "fast": (-362.565, 0.109, 1e-3)}
BANDS = ((0.05, "not at all"), (0.5, "moved partly"), (math.inf, "moved"))


def rng_a(k: int, s: int) -> list[int]:
    return [2026, TAG, 1846, 1, k + 1, s]


def rng_b(k: int, r: int) -> list[int]:
    return [2026, TAG, 1846, 1, k + 1, 100 + r]


def path_a(k: int, s: int, suffix: str, cache: Path) -> Path:
    return cache / "a" / f"k{k}__s{s}{suffix}"


def path_b(k: int, r: int, suffix: str, cache: Path) -> Path:
    return cache / "b" / f"k{k}__r{r}{suffix}"


def pid(k: int, s: int) -> str:
    return f"k{k}__" + ("merged" if s == MERGED else f"s{s}")


def _log(msg: str, cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with (cache / "driver.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


# ---------------------------------------------------------------------------- the problem

_P: dict = {}


def problem(synthetic: int | None = None):
    """(obs, ctx, uncut space, LinearizedProblem, PerturbativeProblem): S1's, or the same
    space on a synthetic data set for smoke runs."""
    key = ("syn", synthetic)
    if key in _P:
        return _P[key]
    if synthetic is None:
        m8_runs.check_frozen(VARIANT)
        out = m8_runs.fast_problem(VARIANT)
    else:
        from fork.inference.priors import free_a_space
        from fork.model.linearized import LinearizedProblem
        from fork.model.perturbative import PerturbativeProblem

        obs, ctx, pq, box = v13.synthetic_dataset(synthetic)
        space = free_a_space(pq, v13.S1_E_MAX, box, inclined=True)
        lp = LinearizedProblem.build(obs, ctx, space)
        pr = PerturbativeProblem.build(
            obs,
            ctx,
            space,
            np.asarray(box.centre, dtype=float),
            step_days=m8_runs.STEP_DAYS,
            quadrature=m8_runs.QUADRATURE,
        )
        out = (obs, ctx, space, lp, pr)
    _P[key] = out
    return out


def subcell(space, k: int):
    from fork.inference.nested import cell_space

    return cell_space(space, SUBCELLS[k])


def prior_masses(space) -> np.ndarray:
    return np.array([subcell(space, k).slice_weight for k in range(len(SUBCELLS))])


def _real(cache: Path) -> bool:
    return cache.resolve() == CACHE.resolve()


def _guard(cache: Path, maxiter: int | None, synthetic: int | None) -> None:
    if _real(cache) and (maxiter is not None or synthetic is not None):
        raise SystemExit("smoke runs (--maxiter, --synthetic) go to a scratch --cache")


# ---------------------------------------------------------------------------- stage A


def do_job_a(
    k: int, s: int, cache: Path = CACHE, maxiter: int | None = None, synthetic: int | None = None
) -> None:
    import dynesty

    from fork import ledger
    from fork.inference.nested import BOUND, DLOGZ, SAMPLE, NestedResult, _ProgressLog

    _guard(cache, maxiter, synthetic)
    out = path_a(k, s, ".pkl", cache)
    if out.exists():
        print(f"a/k{k}/s{s}: cached")
        return
    *_, space, _, pr = problem(synthetic)
    cs = subcell(space, k)
    m8_runs._S.update(pr=pr, space=cs)
    periodic = [i - 4 for i in cs.periodic_unit_dims]
    if min(periodic) < 0 or len(cs.unit_names) - 4 != N_PHI:
        raise RuntimeError("unexpected prior layout")
    ckpt = path_a(k, s, ".ckpt", cache)
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    log = _ProgressLog(cache / "progress.log", m8_runs.PROGRESS_EVERY_S, f"s1_cell0/a/k{k}/s{s}")
    run_kw = dict(
        dlogz=DLOGZ,
        maxiter=maxiter,
        print_progress=True,
        print_func=log,
        checkpoint_file=str(ckpt),
        checkpoint_every=m8_runs.CHECKPOINT_EVERY_S,
    )
    resumed = ckpt.exists()
    t0, c0 = time.time(), ledger.cpu_now()
    if resumed:
        sampler = dynesty.NestedSampler.restore(str(ckpt))
        log.line(f"resumed from {ckpt.name} at iteration {sampler.it}")
        sampler.run_nested(resume=True, **run_kw)
    else:
        log.line(
            f"start: ndim {N_PHI}, nlive {NLIVE}, {BOUND}/{SAMPLE} {SLICES} slices, "
            f"dlogz {DLOGZ}, rng {rng_a(k, s)}, sub-cell {SUBCELLS[k]}"
        )
        sampler = dynesty.NestedSampler(
            m8_runs._fast_logl,
            m8_runs._phi_transform,
            N_PHI,
            nlive=NLIVE,
            bound=BOUND,
            sample=SAMPLE,
            periodic=periodic,
            slices=SLICES,
            rstate=np.random.default_rng(rng_a(k, s)),
        )
        sampler.run_nested(**run_kw)
    res = NestedResult(cs, sampler.results, NLIVE, DLOGZ, time.time() - t0)
    res.save(out)
    ckpt.unlink(missing_ok=True)
    log.line(f"done: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, {res.ncall} calls")
    if _real(cache):
        ledger.append(
            f"s1_cell0/a/k{k}/s{s}",
            "nested",
            n_likelihood_calls=res.ncall,
            n_integrations=0,
            sim_years=0.0,
            wall_s=time.time() - t0,
            cpu_s=ledger.cpu_now() - c0,
            n_workers=1,
            estimated=resumed,
            resumed=resumed,
            nlive=NLIVE,
            niter=res.niter,
            slices=SLICES,
            likelihood="perturbative, second order off",
        )


def _spawn(cmd: list[str], logfile: Path) -> subprocess.Popen:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    with logfile.open("a") as fh:
        return subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)


def _schedule(jobs: list[tuple], n_slots: int, cache: Path, extra: list[str]) -> bool:
    """Run ``jobs`` = (name, done-path, command, slots) with at most ``n_slots`` busy; one
    retry per failed job (it resumes from its checkpoint)."""
    pending = [j for j in jobs if not j[1].exists()]
    _log(f"schedule: {len(pending)} of {len(jobs)} jobs pending, {n_slots} slots", cache)
    running: dict[str, tuple] = {}
    failures: dict[str, int] = {}
    while pending or running:
        for name, (proc, t0, job) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[name]
            if rc == 0 and job[1].exists():
                _log(f"done {name} ({(time.time() - t0) / 60:.1f} min)", cache)
            else:
                failures[name] = failures.get(name, 0) + 1
                _log(f"FAILED {name} (exit {rc}, attempt {failures[name]})", cache)
                if failures[name] < 2:
                    pending.insert(0, job)
        busy = sum(j[3] for _, _, j in running.values())
        while pending and busy + pending[0][3] <= n_slots:
            job = pending.pop(0)
            proc = _spawn([*job[2], "--cache", str(cache), *extra], job[1].with_suffix(".log"))
            running[job[0]] = (proc, time.time(), job)
            busy += job[3]
            _log(f"start {job[0]} (pid {proc.pid}, {job[3]} slots)", cache)
        time.sleep(5 if extra else 30)
    ok = all(j[1].exists() for j in jobs)
    _log(f"schedule finished; failures: {failures or 'none'}; all done: {ok}", cache)
    return ok


def _cmd(*args: str) -> list[str]:
    return [sys.executable, "-m", "fork.export.s1_cell0", *args]


def run_stage_a(n_slots: int, cache: Path, extra: list[str], cells=None) -> bool:
    ks = range(len(SUBCELLS)) if cells is None else cells
    jobs = [
        (f"a/k{k}/s{s}", path_a(k, s, ".pkl", cache), _cmd("job-a", str(k), str(s)), 1)
        for s in SEEDS
        for k in ks
    ]
    return _schedule(jobs, n_slots, cache, extra)


# ---------------------------------------------------------------------------- reweighting


def load_posterior(k: int, s: int, cache: Path):
    from dynesty.utils import merge_runs

    from fork.inference.nested import NestedResult

    if s != MERGED:
        return NestedResult.load(path_a(k, s, ".pkl", cache))
    runs = [NestedResult.load(path_a(k, t, ".pkl", cache)) for t in SEEDS]
    merged = merge_runs([r.results for r in runs], print_progress=False)
    return NestedResult(
        runs[0].space,
        merged,
        sum(r.nlive for r in runs),
        runs[0].dlogz,
        sum(r.wall_s for r in runs),
    )


def draws(k: int, s: int, cache: Path, max_draws: int | None = None) -> dict:
    path = cache / "reweight" / pid(k, s) / "draws.npz"
    if path.exists():
        d = np.load(path)
        return {key: d[key] for key in d.files}
    res = load_posterior(k, s, cache)
    rng = np.random.default_rng([*rng_a(k, s), 7])
    eq = res.equal_weight(rng)
    n = min(N_DRAWS if max_draws is None else max_draws, len(eq))
    phi = eq[rng.choice(len(eq), n, replace=False)]
    d = {
        "phi": phi,
        "ln_z_fast": np.array(res.logz),
        "ln_z_fast_err": np.array(res.logzerr),
        "ncall": np.array(res.ncall),
        "niter": np.array(res.niter),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("draws.tmp.npz")
    np.savez(tmp, **d)
    tmp.replace(path)
    return d


_RW: dict = {}


def _init_reweight(synthetic: int | None) -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from fork.inference.profile import GalleGeometry

    _, ctx, space, lp, pr = problem(synthetic)
    _RW.update(ctx=ctx, space=space, lp=lp, pr=pr, geo=GalleGeometry.from_horizons())


def _reweight_chunk(task: tuple) -> tuple[str, int]:
    """Both likelihoods at draws [lo, hi) of one posterior, and each draw's Galle-night
    longitude through the full model with u drawn from mode (b)'s N(u*, Σ) (as m8_runs)."""
    from fork.inference.profile import galle_lon_of_date_rad

    k, s, c, lo, hi, cache = task
    out = Path(cache) / "reweight" / pid(k, s) / f"chunk_{c:03d}.npz"
    if out.exists():
        return pid(k, s), c
    ctx, space, lp, pr, geo = (_RW[x] for x in ("ctx", "space", "lp", "pr", "geo"))
    phi = np.load(out.parent / "draws.npz")["phi"][lo:hi]
    n = len(phi)
    rng = np.random.default_rng([*rng_a(k, s), 9, c])
    a = {
        "lnl_fast": np.full(n, -np.inf),
        "lnl_b": np.full(n, -np.inf),
        "b_trunc_m2lnl": np.full(n, np.nan),
        "galle_deg": np.full(n, np.nan),
    }
    for i, f in enumerate(phi):
        a["lnl_fast"][i] = pr.log_likelihood(f, second_order=m8_runs.SECOND_ORDER)
        try:
            m = lp.marginal(f, "b")
        except np.linalg.LinAlgError:
            continue
        if not math.isfinite(m.m2lnl):
            continue
        a["lnl_b"][i] = -0.5 * m.m2lnl
        a["b_trunc_m2lnl"][i] = lp.truncation(m)[1]
        u = rng.multivariate_normal(m.u_star, m.cov_u)
        th = space.theta(np.concatenate([u, f]))
        a["galle_deg"][i] = math.degrees(galle_lon_of_date_rad(th, ctx, geo)) % 360.0
    tmp = out.with_name(out.stem + ".tmp.npz")
    np.savez(tmp, **a)
    tmp.replace(out)
    return pid(k, s), c


def reweight(
    n_workers: int,
    cache: Path = CACHE,
    synthetic: int | None = None,
    max_draws: int | None = None,
    cells=None,
) -> None:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    from fork import ledger

    t0, c0 = time.time(), ledger.cpu_now()
    ks = range(len(SUBCELLS)) if cells is None else cells
    tasks = []
    for k in ks:
        for s in (*SEEDS, MERGED):
            n = len(draws(k, s, cache, max_draws)["phi"])
            for c, lo in enumerate(range(0, n, CHUNK)):
                if not (cache / "reweight" / pid(k, s) / f"chunk_{c:03d}.npz").exists():
                    tasks.append((k, s, c, lo, min(lo + CHUNK, n), str(cache)))
    _log(f"reweight: {len(tasks)} chunks to do, {n_workers} workers", cache)
    with ProcessPoolExecutor(n_workers, initializer=_init_reweight, initargs=(synthetic,)) as ex:
        futs = [ex.submit(_reweight_chunk, t) for t in tasks]
        for i, f in enumerate(as_completed(futs)):
            p, c = f.result()
            if (i + 1) % 48 == 0 or i + 1 == len(futs):
                _log(f"reweight: {i + 1}/{len(futs)} chunks (last {p} #{c})", cache)
    if _real(cache) and tasks:
        ledger.append(
            "s1_cell0/reweight",
            "reweight",
            n_likelihood_calls=sum(t[4] - t[3] for t in tasks),
            n_integrations=None,
            sim_years=None,
            wall_s=time.time() - t0,
            cpu_s=ledger.cpu_now() - c0,
            n_workers=n_workers,
            note="mode (b) + fast + Galle longitude per draw",
        )


def load_reweighted(k: int, s: int, cache: Path) -> dict:
    d = draws(k, s, cache)
    n = len(d["phi"])
    parts = [
        np.load(cache / "reweight" / pid(k, s) / f"chunk_{c:03d}.npz")
        for c in range(math.ceil(n / CHUNK))
    ]
    out = dict(d)
    for key in parts[0].files:
        out[key] = np.concatenate([p[key] for p in parts])
    return out


# ---------------------------------------------------------------------------- stage A summary


def posterior_stats(d: dict, rng: np.random.Generator) -> dict:
    """One reweighted posterior: ln Z, ESS/N, the truncation, the Galle-night summary."""
    from fork.inference.reweight import reweight_cell

    cr = reweight_cell(d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]), float(d["ln_z_fast_err"]))
    _, ix = m8_runs._resampled([cr], [0.0], rng, cr.n)
    g = m6c.lon_summary(d["galle_deg"][ix], {})
    finite = np.isfinite(d["lnl_b"])
    return {
        "n": cr.n,
        "ess_frac": cr.ess / cr.n,
        "ln_z_fast": cr.ln_z_fast,
        "ln_z_fast_err": cr.ln_z_fast_err,
        "ln_correction": cr.ln_correction,
        "ln_correction_err": cr.ln_correction_err,
        "ln_z": cr.ln_z,
        "galle_median": g["median"],
        "galle_width68": g["width68"],
        "ncall": int(d["ncall"]),
        "niter": int(d["niter"]),
        "b_trunc_m2lnl_max": float(np.nanmax(d["b_trunc_m2lnl"])) if finite.any() else None,
        "n_b_failed": int(np.sum(~finite)),
    }


def spread(values: list[float], circular: bool = False) -> float:
    v = np.asarray(values, dtype=float)
    if circular:
        v = (v - v[0] + 180.0) % 360.0 - 180.0
    return float(v.max() - v.min())


def evidence_shares(ln_z: np.ndarray, w: np.ndarray) -> np.ndarray:
    t = np.log(w) + np.asarray(ln_z)
    return np.exp(t - (t.max() + math.log(np.exp(t - t.max()).sum())))


def stage_a(cache: Path = CACHE, synthetic: int | None = None, cells=None) -> dict:
    """Per sub-cell: the four seeds and the merge, reweighted; shares; G2, G3 and T1."""
    from fork.inference.nested import ERR_INFLATION

    space = problem(synthetic)[2]
    ks = list(range(len(SUBCELLS)) if cells is None else cells)
    w = prior_masses(space)[ks]
    rows = []
    for k in ks:
        seeds = {}
        for s in (*SEEDS, MERGED):
            st = posterior_stats(
                load_reweighted(k, s, cache), np.random.default_rng([*rng_a(k, s), 13])
            )
            seeds["merged" if s == MERGED else f"s{s}"] = st
        four = [seeds[f"s{s}"] for s in SEEDS]
        rows.append(
            {
                "k": k,
                "a_au": list(SUBCELLS[k][0]),
                "e": list(SUBCELLS[k][1]),
                "prior_mass": float(subcell(space, k).slice_weight),
                "posteriors": seeds,
                "ln_z_range": spread([x["ln_z"] for x in four]),
                "median_range": spread([x["galle_median"] for x in four], True),
            }
        )
    shares = evidence_shares(np.array([r["posteriors"]["merged"]["ln_z"] for r in rows]), w)
    for r, sh in zip(rows, shares, strict=True):
        r["share"] = float(sh)
        r["material"] = bool(sh >= MATERIAL_SHARE)
        r["ln_z_pass"] = r["ln_z_range"] <= TOL["ln_z"]
        r["median_pass"] = r["median_range"] <= TOL["median"]
    posts = [(r["k"], name, p) for r in rows for name, p in r["posteriors"].items()]
    g2 = [
        {"posterior": f"k{k}/{n}", "ess_frac": p["ess_frac"]}
        for k, n, p in posts
        if p["ess_frac"] < ESS_GATE
    ]
    g3 = [
        {"posterior": f"k{k}/{n}", "max": p["b_trunc_m2lnl_max"]}
        for k, n, p in posts
        if p["b_trunc_m2lnl_max"] is None or p["b_trunc_m2lnl_max"] >= TRUNC_GATE
    ]
    t1 = {
        "ln_z_every_subcell": all(r["ln_z_pass"] for r in rows),
        "median_material_subcells": all(r["median_pass"] for r in rows if r["material"]),
        "median_every_subcell_reported_not_read": all(r["median_pass"] for r in rows),
    }
    t1["pass"] = t1["ln_z_every_subcell"] and t1["median_material_subcells"]
    lnz = np.array([r["posteriors"]["merged"]["ln_z"] for r in rows])
    err = np.array(
        [
            ERR_INFLATION
            * math.hypot(
                r["posteriors"]["merged"]["ln_z_fast_err"],
                r["posteriors"]["merged"]["ln_correction_err"],
            )
            for r in rows
        ]
    )
    ln_z0, ln_z0_err = v13.combine_subcells(lnz, err, w)
    return {
        "subcells": rows,
        "k_star": int(rows[int(np.argmax(shares))]["k"]),
        "gates": {
            "G2": {
                "min_ess_frac": min(p["ess_frac"] for *_, p in posts),
                "threshold": ESS_GATE,
                "failures": g2,
                "pass": not g2,
            },
            "G3": {"threshold": TRUNC_GATE, "failures": g3, "pass": not g3},
        },
        "T1": t1,
        "ln_z0": ln_z0,
        "ln_z0_err": ln_z0_err,
    }


# ---------------------------------------------------------------------------- stage B


def space_b(k: int, synthetic: int | None = None):
    """The 11-dimensional space of M6c's S1 job, cut to sub-cell k."""
    return subcell(problem(synthetic)[2], k)


def do_job_b(
    k: int,
    r: int,
    n_workers: int = B_WORKERS,
    cache: Path = CACHE,
    maxiter: int | None = None,
    synthetic: int | None = None,
) -> None:
    from fork import ledger
    from fork.export.to_json import A1_EPOCH_JD_TDB
    from fork.inference.nested import NestedResult, run_nested

    _guard(cache, maxiter, synthetic)
    out = path_b(k, r, ".pkl", cache)
    if out.exists():
        print(f"b/k{k}/r{r}: cached")
        return
    obs, ctx, *_ = problem(synthetic)
    ckpt = path_b(k, r, ".ckpt", cache)
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    resumed = ckpt.exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(
        obs,
        ctx,
        space_b(k, synthetic),
        np.random.default_rng(rng_b(k, r)),
        nlive=NLIVE,
        slices=SLICES,
        n_workers=n_workers,
        maxiter=maxiter,
        progress_log=cache / "progress.log",
        checkpoint_file=ckpt,
        label=f"s1_cell0/b/k{k}/r{r}",
    )
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    ckpt.unlink(missing_ok=True)
    if _real(cache):
        yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
        ledger.append(
            f"s1_cell0/b/k{k}/r{r}",
            "nested",
            n_likelihood_calls=res.ncall,
            n_integrations=res.ncall,
            sim_years=res.ncall * yrs,
            wall_s=time.time() - t0,
            cpu_s=ledger.cpu_now() - c0,
            n_workers=n_workers,
            estimated=resumed,
            resumed=resumed,
            nlive=NLIVE,
            niter=res.niter,
            slices=SLICES,
        )


def run_stage_b(k: int, n_workers: int, cache: Path, extra: list[str]) -> bool:
    w = max(1, min(B_WORKERS, n_workers // len(B_RUNS)))
    jobs = [
        (
            f"b/k{k}/r{r}",
            path_b(k, r, ".pkl", cache),
            _cmd("job-b", str(k), str(r), "--n-workers", str(w)),
            w,
        )
        for r in B_RUNS
    ]
    return _schedule(jobs, n_workers, cache, extra)


def within(d: float, key: str) -> dict:
    return {
        "value": d,
        "tolerance": TOL[key],
        "pass": abs(d) <= TOL[key],
        "at_the_tolerance": abs(abs(d) - TOL[key]) <= AT_TOL[key],
    }


def stage_b(k: int, a: dict, cache: Path = CACHE, synthetic: int | None = None) -> dict:
    """T2: the two full-likelihood runs against each other, and their merge against the fast
    four-seed merge of sub-cell k."""
    from fork.export import m6c_repeat
    from fork.inference.nested import ERR_INFLATION, NestedResult
    from fork.inference.profile import GalleGeometry

    ctx = problem(synthetic)[1]
    geo = GalleGeometry.from_horizons()
    rng = np.random.default_rng([2026, TAG, 1846, 1, k + 1, 0, 11])
    runs = [NestedResult.load(path_b(k, r, ".pkl", cache)) for r in B_RUNS]
    per = [m6c_repeat.run_stats(res, ctx, geo, rng) for res in runs]
    merged = m6c_repeat.run_stats(m6c_repeat.merged_cell2(runs), ctx, geo, rng)
    mean_w = 0.5 * (per[0]["width68"] + per[1]["width68"])
    mutual = {
        "ln_z": within(per[0]["logz"] - per[1]["logz"], "ln_z"),
        "median": within((per[0]["median"] - per[1]["median"] + 180.0) % 360.0 - 180.0, "median"),
        "width_rel": within(max(abs(p["width68"] / mean_w - 1.0) for p in per), "width_rel"),
    }
    fast = next(r for r in a["subcells"] if r["k"] == k)["posteriors"]["merged"]
    vs_fast = {
        "ln_z": within(merged["logz"] - fast["ln_z"], "ln_z"),
        "median": within(
            (merged["median"] - fast["galle_median"] + 180.0) % 360.0 - 180.0, "median"
        ),
        "width_rel": within(merged["width68"] / fast["galle_width68"] - 1.0, "width_rel"),
    }
    ok = all(x["pass"] for part in (mutual, vs_fast) for x in part.values())
    return {
        "k": k,
        "runs": {f"r{r}": p for r, p in zip(B_RUNS, per, strict=True)},
        "merged": merged,
        "merged_ln_z_err_inflated": ERR_INFLATION * merged["logzerr"],
        "fast_merge": fast,
        "mutual": mutual,
        "vs_fast": vs_fast,
        "pass": ok,
    }


# ---------------------------------------------------------------------------- reading


_LOG_MASS: list[float] = []


def p_below_25(ln_z0: float, other: tuple[float, ...]) -> float:
    """Cell 0's share with M6c's prior masses of S1's four cells: the function of
    ``m8_arbiter.implied_p_below_25`` with the masses computed once."""
    if not _LOG_MASS:
        _LOG_MASS.extend(math.log(m8_runs._cell_space(VARIANT, k).slice_weight) for k in range(4))
    t = np.array(_LOG_MASS) + np.array([ln_z0, *other])
    return float(math.exp(t[0] - (t.max() + math.log(np.exp(t - t.max()).sum()))))


def word(p: float) -> str:
    return next(name for edge, name in BANDS if p < edge)


def reading(a: dict, b: dict | None, g1: bool) -> dict:
    """The plan's § 5, in its order: gates, T1, T2, then the resolved values."""
    checks = {
        "published": p_below_25(CHECK_VALUES["published"][0], OTHER_PUBLISHED),
        "fast": p_below_25(CHECK_VALUES["fast"][0], OTHER_FAST),
    }
    if any(abs(checks[n] - CHECK_VALUES[n][1]) > CHECK_VALUES[n][2] for n in checks):
        raise RuntimeError(f"the check values of M8 are not reproduced: {checks}")
    p = {
        "against_published_other_cells": p_below_25(a["ln_z0"], OTHER_PUBLISHED),
        "against_fast_other_cells": p_below_25(a["ln_z0"], OTHER_FAST),
    }
    words = {n: word(x) for n, x in p.items()}
    gates_ok = g1 and a["gates"]["G2"]["pass"] and a["gates"]["G3"]["pass"]
    if not gates_ok:
        outcome = "not resolved, the method failed (a gate)"
    elif not a["T1"]["pass"]:
        outcome = "not resolved, the sub-cells do not make the runs repeatable (T1)"
    elif b is None or not b["pass"]:
        outcome = "not resolved, the fast likelihood and the full one disagree in k* (T2)"
    else:
        outcome = "resolved"
    resolved = outcome == "resolved"
    same = len(set(words.values())) == 1
    verdict = (
        "the published verdict stands; unresolved"
        if not resolved
        else next(iter(words.values()))
        if same
        else "between " + " and ".join(sorted(set(words.values())))
    )
    low = resolved and all(x < BANDS[0][0] for x in p.values())
    label = (
        "provisional: the low-distance solution is no longer a reason; the fourth "
        "condition of the higher-resolution reruns' rule is"
        if low
        else "provisional: for both reasons (the fourth condition of the higher-resolution "
        "reruns' rule, and the low-distance solution under e ≤ 0.5)"
    )
    return {
        "check_values_reproduced": checks,
        "ln_z0": a["ln_z0"],
        "ln_z0_err": a["ln_z0_err"],
        "p_a_below_25_au": p,
        "words": words,
        "G1": g1,
        "outcome": outcome,
        "resolved": resolved,
        "verdict": verdict,
        "label_95": label,
    }


def summarize(cache: Path = CACHE, synthetic: int | None = None, cells=None) -> dict:
    a = stage_a(cache, synthetic, cells)
    k = a["k_star"]
    have_b = all(path_b(k, r, ".pkl", cache).exists() for r in B_RUNS)
    try:
        b = stage_b(k, a, cache, synthetic) if have_b else None
    except IndexError:
        if _real(cache):
            raise
        b = None  # a smoke run cut short by --maxiter has no samples near its maximum
    ks = range(len(SUBCELLS)) if cells is None else cells
    g1 = all(path_a(kk, s, ".pkl", cache).exists() for kk in ks for s in SEEDS) and have_b
    out = {
        "plan": "docs/notes/S1_cell0.md",
        "settings": {
            "nlive": NLIVE,
            "slices": SLICES,
            "seeds": list(SEEDS),
            "n_draws": N_DRAWS,
            "tolerances": TOL,
            "material_share": MATERIAL_SHARE,
        },
        "stage_a": a,
        "stage_b": b,
        "reading": reading(a, b, g1),
    }
    target = OUT_JSON if _real(cache) else cache / "S1_cell0.json"
    target.write_text(
        json.dumps(out, indent=1, default=float, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    rd = out["reading"]
    _log(f"summarize: {rd['outcome']}; verdict: {rd['verdict']}", cache)
    return out


def pipeline(
    n_workers: int,
    cache: Path = CACHE,
    maxiter: int | None = None,
    synthetic: int | None = None,
    max_draws: int | None = None,
    cells=None,
) -> None:
    _guard(cache, maxiter, synthetic)
    extra = ([] if maxiter is None else ["--maxiter", str(maxiter)]) + (
        [] if synthetic is None else ["--synthetic", str(synthetic)]
    )
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    _log(f"pipeline: commit {head}, {n_workers} workers", cache)
    if not run_stage_a(n_workers, cache, extra, cells):
        _log("pipeline: stage A incomplete; stopped", cache)
        return
    reweight(n_workers, cache, synthetic, max_draws, cells)
    k = stage_a(cache, synthetic, cells)["k_star"]
    _log(f"pipeline: k* = {k} (largest w_k Z_k among the merged sub-cells)", cache)
    if not run_stage_b(k, n_workers, cache, extra):
        _log("pipeline: stage B incomplete; stopped", cache)
        return
    summarize(cache, synthetic, cells)
    _log("pipeline: done", cache)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "pipeline", "job-a", "job-b", "reweight", "summarize"):
        p = sub.add_parser(name)
        p.add_argument("--cache", type=Path, default=CACHE)
        p.add_argument("--maxiter", type=int, default=None, help="smoke runs only")
        p.add_argument("--synthetic", type=int, default=None, help="smoke runs only")
        p.add_argument("--max-draws", type=int, default=None, help="smoke runs only")
        p.add_argument("--cells", type=int, nargs="*", default=None, help="smoke runs only")
        p.add_argument("--n-workers", type=int, default=B_WORKERS if name == "job-b" else 12)
        if name in ("job-a", "job-b"):
            p.add_argument("k", type=int, choices=range(len(SUBCELLS)))
            p.add_argument("s", type=int)
    a = ap.parse_args(argv)
    cache = a.cache.resolve()
    if _real(cache) and (a.max_draws is not None or a.cells is not None):
        raise SystemExit("smoke options go to a scratch --cache")
    if a.cmd == "plan":
        for k, c in enumerate(SUBCELLS):
            done = sum(path_a(k, s, ".pkl", cache).exists() for s in SEEDS)
            print(f"k{k}: a {c[0]}, e {c[1]}: {done}/4 fast runs done")
    elif a.cmd == "pipeline":
        pipeline(a.n_workers, cache, a.maxiter, a.synthetic, a.max_draws, a.cells)
    elif a.cmd == "job-a":
        do_job_a(a.k, a.s, cache, a.maxiter, a.synthetic)
    elif a.cmd == "job-b":
        do_job_b(a.k, a.s, a.n_workers, cache, a.maxiter, a.synthetic)
    elif a.cmd == "reweight":
        reweight(a.n_workers, cache, a.synthetic, a.max_draws, a.cells)
    else:
        summarize(cache, a.synthetic, a.cells)


if __name__ == "__main__":
    main()
