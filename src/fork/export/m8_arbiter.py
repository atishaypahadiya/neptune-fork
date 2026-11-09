"""M8 A8d(S1) arbiter (spec, M8 "If A8d fails"; docs/notes/M8.md, "Follow-up
pre-registration"): a full-likelihood repeat of S1's cell 0 at R2's settings, 2 seeds,
compared with M8's S1 cell 0 under A8b's tolerances.

    uv run python -m fork.export.m8_arbiter plan
    uv run python -m fork.export.m8_arbiter run [--n-workers 12]   # both runs, 6 workers each
    uv run python -m fork.export.m8_arbiter job 1 [--n-workers 6] [--maxiter N --cache DIR]
    uv run python -m fork.export.m8_arbiter summarize

- **Problem:** exactly M6c's S1 cell-0 job: ``m6c.dataset("s1_e05")``, ``m6c.job_space``
  (11 dimensions, u sampled, full likelihood), ``run_nested``. M6c's code is imported, not
  changed.
- **Settings:** R2's, literally: nlive 4000 and 28 slices (published S1 cell 0: nlive 2000,
  dynesty's default 14 slices). Bound multi, rslice, dlogz 0.1. ERR_INFLATION as in M6c.
- **Seeds:** rng [2026, 6, 1846, 1, 1, 3, r] for r = 1, 2: M6c's S1 cell-0 seed, then 3 (the
  arbiter; R1 and R2 used 1 and 2), then the run.
- **Per run and for the 2-run merge** (``dynesty.utils.merge_runs``): the Galle-night median and
  68 % width with ``m6c_repeat.run_stats`` (as for R2, which set A8b's reference), ln Z.
- **Test:** the merge against M8's S1 cell 0 (4-seed merge, reweighted; site/data/m8_runs.json)
  under A8b's tolerances: median within 0.35°, 68 % width within 10 %, ln Z within 0.5. The
  published single run (ln Z −366.84) is reported beside it.
- **What it implies for the reading rule:** P(a < 25 AU) is cell 0's share (cell 0 is
  a = 20–25 AU, the other cells a ≥ 25 AU). Reported with the arbiter's ln Z₀ in place of
  M8's and of the published run's, the other cells unchanged.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from fork.export import m6c
from fork.export.to_json import A1_EPOCH_JD_TDB, write_json

VARIANT = "s1_e05"
CELL = 0
NLIVE = 4000
SLICES = 28
SEEDS = (1, 2)
GROUP_INDEX = 3  # after m6c_repeat's R1 (1) and R2 (2)
WORKERS = 6
CACHE = m6c.ROOT / "data" / "cache" / "nested" / "m8_arbiter"
OUT_NAME = "m8_arbiter.json"
A8B_TOL = {"median": 0.35, "width_rel": 0.10, "ln_z": 0.5}
PUBLISHED_LN_Z0 = -366.84  # site/data/m6_fits.json, s1_e05 cell 0 (logz_k)
SUMMARY_RNG = (2026, 8, 1846, 1, 1, GROUP_INDEX, 0)


def rng_seed(r: int) -> list[int]:
    return [*m6c.Job(VARIANT, CELL).rng_seed(), GROUP_INDEX, r]


def path(r: int, suffix: str, cache: Path = CACHE) -> Path:
    return cache / f"run{r}{suffix}"


def problem():
    """(observations, context, space) of S1's cell 0, exactly as M6c's job built it."""
    m6c.check_frozen(VARIANT)
    obs, ctx, _, pq, box = m6c.dataset(variant=VARIANT)
    return obs, ctx, m6c.job_space(m6c.Job(VARIANT, CELL), pq, box)


def _log(msg: str, cache: Path = CACHE) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with (cache / "driver.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def do_job(r: int, n_workers: int, maxiter: int | None = None, cache: Path = CACHE) -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult, run_nested

    if maxiter is not None and cache == CACHE:
        raise SystemExit("smoke runs (--maxiter) go to a scratch --cache, never the run cache")
    out = path(r, ".pkl", cache)
    if out.exists():
        print(f"run{r}: cached")
        return
    obs, ctx, space = problem()
    ckpt = path(r, ".ckpt", cache)
    resumed = ckpt.exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(obs, ctx, space, np.random.default_rng(rng_seed(r)), nlive=NLIVE,
                     slices=SLICES, n_workers=n_workers, maxiter=maxiter,
                     progress_log=cache / "progress.log", checkpoint_file=ckpt,
                     label=f"m8_arbiter/run{r}")  # fmt: skip
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    ckpt.unlink(missing_ok=True)
    if cache == CACHE:
        yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
        ledger.append(f"m8_arbiter/run{r}", "nested", n_likelihood_calls=res.ncall,
                      n_integrations=res.ncall, sim_years=res.ncall * yrs,
                      wall_s=time.time() - t0, cpu_s=ledger.cpu_now() - c0,
                      n_workers=n_workers, estimated=resumed, resumed=resumed, nlive=NLIVE,
                      niter=res.niter, slices=SLICES)  # fmt: skip
    print(f"run{r}: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, {res.ncall} calls, "
          f"{(time.time() - t0) / 60:.1f} min")  # fmt: skip


def run_all(n_workers: int, cache: Path = CACHE, maxiter: int | None = None) -> bool:
    """Both runs at once, ``WORKERS`` each within the budget; one retry from the checkpoint."""
    m6c.check_frozen(VARIANT)
    pending = [r for r in SEEDS if not path(r, ".pkl", cache).exists()]
    _log(f"run: {len(pending)} runs pending, budget {n_workers} workers", cache)
    running: dict[int, tuple] = {}
    failures: dict[int, int] = {}
    while pending or running:
        for r, (proc, t0, _) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[r]
            if rc == 0 and path(r, ".pkl", cache).exists():
                _log(f"done run{r} ({(time.time() - t0) / 60:.1f} min)", cache)
            else:
                failures[r] = failures.get(r, 0) + 1
                _log(f"FAILED run{r} (exit {rc}, attempt {failures[r]})", cache)
                if failures[r] < 2:
                    pending.insert(0, r)
        busy = sum(w for *_, w in running.values())
        while pending and (not running or busy + WORKERS <= n_workers):
            r = pending.pop(0)
            w = min(WORKERS, n_workers)
            cmd = [sys.executable, "-m", "fork.export.m8_arbiter", "job", str(r),
                   "--n-workers", str(w), "--cache", str(cache)]  # fmt: skip
            if maxiter:
                cmd += ["--maxiter", str(maxiter)]
            with path(r, ".log", cache).open("a") as fh:
                proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=m6c.ROOT)
            running[r] = (proc, time.time(), w)
            busy += w
            _log(f"start run{r} (nlive {NLIVE}, slices {SLICES}, {w} workers, rng "
                 f"{rng_seed(r)}, pid {proc.pid})", cache)  # fmt: skip
        time.sleep(5 if maxiter else 30)
    done = all(path(r, ".pkl", cache).exists() for r in SEEDS)
    _log(f"run finished; failures: {failures or 'none'}; all done: {done}", cache)
    return done


def implied_p_below_25(ln_z0: float, other_ln_z: list[float]) -> float:
    """Cell 0's share (= P(a < 25 AU)) with ``ln_z0`` for cell 0 and the other cells' ln Z_k,
    combined with M6c's prior masses (slice weights)."""
    from scipy.special import logsumexp

    from fork.export import m8_runs

    lc = np.array([math.log(m8_runs._cell_space(VARIANT, k).slice_weight) for k in range(4)])
    lz = np.array([ln_z0, *other_ln_z])
    return float(math.exp(lc[0] + lz[0] - logsumexp(lc + lz)))


def summarize(cache: Path = CACHE) -> dict:
    from fork import ledger
    from fork.export import m6c_repeat
    from fork.inference.nested import ERR_INFLATION, NestedResult
    from fork.inference.profile import GalleGeometry

    m6c.check_frozen(VARIANT)
    t0, c0 = time.time(), ledger.cpu_now()
    runs = [NestedResult.load(path(r, ".pkl", cache)) for r in SEEDS]
    _, ctx, _ = problem()
    geo = GalleGeometry.from_horizons()
    rng = np.random.default_rng(list(SUMMARY_RNG))
    per_run = {f"run{r}": m6c_repeat.run_stats(res, ctx, geo, rng)
               for r, res in zip(SEEDS, runs, strict=True)}  # fmt: skip
    merged = m6c_repeat.run_stats(m6c_repeat.merged_cell2(runs), ctx, geo, rng)
    m8 = json.loads((m6c.ROOT / "site" / "data" / "m8_runs.json").read_text())
    m8c = {(c["cell"], c["seed"]): c for c in m8["cells"] if c["variant"] == VARIANT}
    m8_c0 = m8c[(CELL, 0)]
    pub = json.loads((m6c.ROOT / "site" / "data" / "m6_fits.json").read_text())
    pub_lnz = [c["logz_k"] for c in pub["fits"][VARIANT]["evidence"]["cells"]]
    test = {"arbiter": {"median": merged["median"], "width68": merged["width68"],
                        "ln_z": merged["logz"],
                        "ln_z_err_inflated": ERR_INFLATION * merged["logzerr"]},
            "m8": {"median": m8_c0["galle_median"], "width68": m8_c0["galle_width68"],
                   "ln_z": m8_c0["ln_z"]},
            "d_median": merged["median"] - m8_c0["galle_median"],
            "width_rel": merged["width68"] / m8_c0["galle_width68"] - 1.0,
            "d_ln_z": merged["logz"] - m8_c0["ln_z"], "tolerances": A8B_TOL}  # fmt: skip
    test["pass"] = bool(abs(test["d_median"]) <= A8B_TOL["median"]
                        and abs(test["width_rel"]) <= A8B_TOL["width_rel"]
                        and abs(test["d_ln_z"]) <= A8B_TOL["ln_z"])  # fmt: skip
    m8_other = [m8c[(k, 0)]["ln_z"] for k in (1, 2, 3)]
    reading = {
        "published": {"ln_z0": pub_lnz[0], "p_a_below_25_au":
                      implied_p_below_25(pub_lnz[0], pub_lnz[1:])},
        "m8": {"ln_z0": m8_c0["ln_z"], "p_a_below_25_au":
               implied_p_below_25(m8_c0["ln_z"], m8_other)},
        "arbiter_with_m8_other_cells": implied_p_below_25(merged["logz"], m8_other),
        "arbiter_with_published_other_cells": implied_p_below_25(merged["logz"], pub_lnz[1:]),
        "rule": "P(a < 25 AU) < 0.05: not moved; 0.05-0.5: moved partly; > 0.5: moved"}  # fmt: skip
    payload = {"description": "M8 A8d(S1) arbiter: full-likelihood repeat of S1 cell 0 at R2's "
               "settings (nlive 4000, 28 slices), 2 seeds, merged, against M8's S1 cell 0 "
               "under A8b's tolerances; published single run beside it.",
               "settings": {"variant": VARIANT, "cell": CELL, "nlive": NLIVE, "slices": SLICES,
                            "seeds": [rng_seed(r) for r in SEEDS],
                            "err_inflation": ERR_INFLATION},
               "runs": per_run, "merged": merged, "published_ln_z0": PUBLISHED_LN_Z0,
               "test": test, "reading_rule": reading}  # fmt: skip
    if cache == CACHE:
        write_json(OUT_NAME, payload)
        ledger.append("m8_arbiter/summarize", "summarize", n_likelihood_calls=0,
                      n_integrations=0, sim_years=None, wall_s=time.time() - t0,
                      cpu_s=ledger.cpu_now() - c0, n_workers=1)  # fmt: skip
    else:
        (cache / OUT_NAME).write_text(json.dumps(payload, indent=1, default=float))
    _log(f"summarize: test pass {test['pass']}", cache)
    return payload


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "run", "job", "summarize", "pipeline"):
        p = sub.add_parser(name)
        p.add_argument("--cache", type=Path, default=CACHE)
        if name in ("run", "job", "pipeline"):
            p.add_argument("--n-workers", type=int, default=12 if name != "job" else WORKERS)
            p.add_argument("--maxiter", type=int, default=None, help="smoke runs only")
        if name == "job":
            p.add_argument("r", type=int, choices=SEEDS)
    a = ap.parse_args(argv)
    cache = a.cache.resolve()
    if a.cmd == "plan":
        for r in SEEDS:
            state = ("done" if path(r, ".pkl", cache).exists() else "resumable"
                     if path(r, ".ckpt", cache).exists() else "pending")  # fmt: skip
            print(f"run{r}: nlive {NLIVE}, slices {SLICES}, rng {rng_seed(r)}: {state}")
    elif a.cmd == "run":
        run_all(a.n_workers, cache, a.maxiter)
    elif a.cmd == "job":
        do_job(a.r, a.n_workers, a.maxiter, cache)
    elif a.cmd == "summarize":
        summarize(cache)
    else:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=m6c.ROOT, capture_output=True,
                              text=True).stdout.strip()  # fmt: skip
        _log(f"pipeline: commit {head}", cache)
        if run_all(a.n_workers, cache, a.maxiter):
            summarize(cache)
        _log("pipeline: done", cache)


if __name__ == "__main__":
    main()
