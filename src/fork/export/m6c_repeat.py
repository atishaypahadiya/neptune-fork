"""M6c cell-2 repeatability study → ``site/data/m6c_repeat.json``.

    uv run python -m fork.export.m6c_repeat plan
    uv run python -m fork.export.m6c_repeat run [--n-workers 12]
    uv run python -m fork.export.m6c_repeat job JOB_ID [--n-workers 6]
    uv run python -m fork.export.m6c_repeat summarize
    uv run python -m fork.export.m6c_repeat derived   # after the fact, outcome replace only

**Pre-registered** in docs/notes/M6.md ("M6c cell-2 repeatability: pre-registration"),
committed before any of its sampling. Nothing below is retuned after a result is seen.

Why: the headline's cell 2 (a 28–60 AU, e < 0.12) and S1's cell 2 sample the same prior and
likelihood, but their Galle-night 68 % widths are 3.37° and 5.75°, about 35 times the runs'
own resample errors (M6.md, "Cell 2 does not reproduce between runs").

Jobs, all on the headline's cell 2 (same data, frozen error model, box, prior, likelihood and
dlogz as ``fork.export.m6c``):
- **R1**, 4 runs at exactly the headline settings: nlive 1000, multi/rslice with dynesty's
  default 3 + ndim = 14 slices.
- **R2**, 2 runs with stronger exploration: nlive 4000 **and** 28 slices (2 × the default).

Seeds: the headline cell 2's ``[2026, 6, 1846, 0, 3]`` extended by (group, run):
``[2026, 6, 1846, 0, 3, 1, r]`` for R1 r = 1–4 and ``[2026, 6, 1846, 0, 3, 2, r]`` for R2
r = 1–2.

The analysis (``summarize``) and the decision rule (``decide``) are fixed here and in M6.md.
The replacement candidate is the merge of R2's runs only; R1 is diagnostic.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass

import numpy as np

from fork.export import m6c
from fork.export.to_json import A1_EPOCH_JD_TDB, SITE_DATA, write_json

CACHE = m6c.ROOT / "data" / "cache" / "nested" / "m6c_repeat"
PROGRESS_LOG = CACHE / "progress.log"
DRIVER_LOG = CACHE / "driver.log"
OUT_NAME = "m6c_repeat.json"

HEADLINE_CELL = 2
HEADLINE_SEED = [2026, 6, 1846, 0, HEADLINE_CELL + 1]  # m6c.Job("headline", 2).rng_seed()
CELL2_NDIM = 11  # inclined perturber: 4 Uranus + 7 perturber coordinates
DEFAULT_SLICES = 3 + CELL2_NDIM  # dynesty's rslice default
GROUP_INDEX = {"r1": 1, "r2": 2}
GROUP_RUNS = {"r1": 4, "r2": 2}
GROUP_NLIVE = {"r1": 1000, "r2": 4000}
GROUP_SLICES = {"r1": None, "r2": 2 * DEFAULT_SLICES}  # None: dynesty's default, as the headline
GROUP_WORKERS = {"r1": 3, "r2": 6}

# Analysis, fixed at pre-registration.
LOGL_WINDOW = 30.0  # Galle longitudes for samples within this of max ln L (weight beyond ~1e-9)
N_BOOT = 200  # resample_run draws per run: dynesty's own error on median and width
R2_WIDTH_TOL_REL = 0.10  # R2 consistent: |w1 − w2| ≤ 10 % of their mean
R1_CHI2_P = 0.01  # R1 spread exceeds dynesty's error: χ² of the widths, p < 0.01
REPORT_RNG = (*m6c.SUMMARY_RNG, 10)  # derived samples of the reported interval


@dataclass(frozen=True)
class RepeatJob:
    group: str
    r: int

    @property
    def id(self) -> str:
        return f"{self.group}_{self.r}"

    @staticmethod
    def from_id(job_id: str) -> RepeatJob:
        group, r = job_id.split("_")
        return RepeatJob(group, int(r))

    @property
    def nlive(self) -> int:
        return GROUP_NLIVE[self.group]

    @property
    def slices(self) -> int | None:
        return GROUP_SLICES[self.group]

    @property
    def workers(self) -> int:
        return GROUP_WORKERS[self.group]

    def rng_seed(self) -> list[int]:
        return [*HEADLINE_SEED, GROUP_INDEX[self.group], self.r]

    def path(self, suffix: str):
        return CACHE / (self.id + suffix)


JOBS = tuple(RepeatJob(g, r) for g in ("r1", "r2") for r in range(1, GROUP_RUNS[g] + 1))


def cell2_problem():
    """(observations, context, space) of the headline's cell 2, exactly as its job ran."""
    m6c.check_frozen("headline")
    obs, ctx, _, pq, box = m6c.dataset(variant="headline")
    space = m6c.job_space(m6c.Job("headline", HEADLINE_CELL), pq, box)
    if space.ndim != CELL2_NDIM:
        raise RuntimeError(f"cell-2 space has {space.ndim} dims, pre-registered {CELL2_NDIM}")
    return obs, ctx, space


def do_job(job: RepeatJob, n_workers: int, maxiter: int | None = None) -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult, run_nested

    out = job.path(".pkl")
    if out.exists():
        print(f"{job.id}: cached")
        return
    obs, ctx, space = cell2_problem()
    resumed = job.path(".ckpt").exists()
    t0, c0 = time.time(), ledger.cpu_now()
    res = run_nested(obs, ctx, space, np.random.default_rng(job.rng_seed()), nlive=job.nlive,
                     slices=job.slices, n_workers=n_workers, maxiter=maxiter,
                     progress_log=PROGRESS_LOG,
                     checkpoint_file=None if maxiter else job.path(".ckpt"),
                     label=f"m6c_repeat/{job.id}")  # fmt: skip
    if maxiter is not None:
        print(f"{job.id}: smoke run, {res.niter} it, {res.ncall} calls, {res.wall_s:.1f} s")
        return
    NestedResult(res.space, res.results, res.nlive, res.dlogz, res.wall_s).save(out)
    yrs = ledger.sim_years_per_integration(obs.jd_tdb, A1_EPOCH_JD_TDB)
    ledger.append(f"m6c_repeat/{job.id}", "nested", n_likelihood_calls=res.ncall,
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


def run_all(n_workers: int) -> None:
    """As ``m6c.run_all``: one subprocess per job within a worker budget, the costliest
    first (R2), one retry per failed job from its checkpoint."""
    m6c.check_frozen("headline")
    pending = sorted([j for j in JOBS if not j.path(".pkl").exists()],
                     key=lambda j: (-j.nlive, j.id))  # fmt: skip
    _log(f"run: {len(pending)} jobs pending, budget {n_workers} workers")
    running: dict[str, tuple[subprocess.Popen, float, int]] = {}
    failures: dict[str, int] = {}
    while pending or running:
        for jid, (proc, t0, _) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[jid]
            job = RepeatJob.from_id(jid)
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
                    [sys.executable, "-m", "fork.export.m6c_repeat", "job", job.id,
                     "--n-workers", str(w)], stdout=fh, stderr=subprocess.STDOUT,
                    cwd=m6c.ROOT)  # fmt: skip
            running[job.id] = (proc, time.time(), w)
            busy += w
            _log(f"start {job.id} (nlive {job.nlive}, slices {job.slices or 'default'}, "
                 f"{w} workers, pid {proc.pid})")  # fmt: skip
        time.sleep(30)
    _log(f"run finished; failures: {failures or 'none'}")


def plan() -> None:
    for j in JOBS:
        state = ("done" if j.path(".pkl").exists() else
                 "resumable" if j.path(".ckpt").exists() else "pending")  # fmt: skip
        print(f"{j.id:5s} nlive {j.nlive:5d}  slices {j.slices or DEFAULT_SLICES:3d}  "
              f"workers {j.workers:2d}  seed {j.rng_seed()}  {state}")  # fmt: skip
    print(f"{len(JOBS)} jobs, {sum(j.path('.pkl').exists() for j in JOBS)} done")


# ------------------------------------------------------------------ analysis (pre-registered)


def weighted_quantiles(x: np.ndarray, w: np.ndarray, q) -> np.ndarray:
    o = np.argsort(x)
    c = np.cumsum(w[o])
    return np.interp(q, c / c[-1], x[o])


def spread(values, errors) -> dict:
    """Run-to-run spread of one quantity against the runs' own error estimates (the
    comparison nestcheck makes): sample std (ddof 1), the rms of the per-run errors, their
    ratio, χ² = Σ (x_i − x̄)² / σ_i² with its p-value on n − 1 dof, and the implementation
    error √max(0, std² − mean σ²)."""
    from scipy import stats

    x, s = np.asarray(values, dtype=float), np.asarray(errors, dtype=float)
    n = len(x)
    std = float(x.std(ddof=1))
    rms = float(np.sqrt(np.mean(s**2)))
    chi2 = float(np.sum((x - x.mean()) ** 2 / s**2))
    return {"n": n, "mean": float(x.mean()), "std": std, "rms_run_error": rms,
            "ratio": std / rms, "chi2": chi2, "dof": n - 1,
            "p": float(stats.chi2.sf(chi2, n - 1)),
            "implementation_error": math.sqrt(max(0.0, std**2 - float(np.mean(s**2))))}  # fmt: skip


def decide(r1_widths, r1_width_errors, r2_widths) -> dict:
    """The pre-registered rule, on (a) R2's two widths agree to within ``R2_WIDTH_TOL_REL`` of
    their mean, and (b) R1's widths scatter more than dynesty's own errors allow
    (χ² p < ``R1_CHI2_P``).

    - (a) and (b), outcome ``replace``: the merge of R2's runs only replaces the headline's
      cell 2 in the reported interval. R1 is diagnostic and is not merged into it.
    - (a) without (b), outcome ``cross_check``: the headline is reported as it is, with R1's
      run-to-run std of the median and width as a systematic, and R2's merged width beside it
      as a cross-check, not a replacement.
    - otherwise, outcome ``systematic_only``: the headline as it is, with R1's spread as a
      systematic."""
    w2 = np.asarray(r2_widths, dtype=float)
    rel = float(abs(w2[0] - w2[1]) / w2.mean())
    s1 = spread(r1_widths, r1_width_errors)
    a = rel <= R2_WIDTH_TOL_REL
    b = s1["p"] < R1_CHI2_P
    outcome = "replace" if a and b else "cross_check" if a else "systematic_only"
    return {"r2_width_rel_diff": rel, "r2_consistent": a, "r1_width_chi2": s1["chi2"],
            "r1_width_p": s1["p"], "r1_spread_exceeds_run_error": b,
            "replace_headline_cell2": a and b, "outcome": outcome}  # fmt: skip


def run_stats(res, ctx, geo, rng: np.random.Generator) -> dict:
    """One cell-2 run: Galle-night median and 68 % width from all its weighted samples (within
    ``LOGL_WINDOW``), their std over ``N_BOOT`` resample_run draws, ln Z, a and e."""
    from dynesty.utils import resample_run

    from fork.inference.profile import galle_lon_of_date_rad

    r = res.results
    logl = np.asarray(r["logl"])
    if np.any(np.diff(logl) < 0):
        raise RuntimeError("ln L not sorted: not a static run")
    keep = logl > logl.max() - LOGL_WINDOW
    g = np.full(len(logl), np.nan)
    for i in np.flatnonzero(keep):
        g[i] = math.degrees(galle_lon_of_date_rad(res.space.theta(r["samples"][i]), ctx, geo))
    w = np.asarray(r.importance_weights())
    q = weighted_quantiles(g[keep], w[keep], [0.16, 0.5, 0.84])
    boot = []
    for _ in range(N_BOOT):
        rb = resample_run(r, rstate=rng)
        idx = np.searchsorted(logl, rb["logl"])
        gb, wb = g[idx], np.asarray(rb.importance_weights())
        m = ~np.isnan(gb)
        qb = weighted_quantiles(gb[m], wb[m], [0.16, 0.5, 0.84])
        boot.append((qb[1], qb[2] - qb[0]))
    boot = np.array(boot)
    th = np.array([res.space.theta(s) for s in r["samples"]])
    return {"median": float(q[1]), "width68": float(q[2] - q[0]),
            "median_err": float(boot[:, 0].std()), "width68_err": float(boot[:, 1].std()),
            "logz": res.logz, "logzerr": res.logzerr, "min_m2lnl": -2.0 * res.max_logl,
            "niter": res.niter, "ncall": res.ncall, "nlive": res.nlive,
            "weight_outside_window": float(w[~keep].sum()),
            "a_16_50_84": weighted_quantiles(th[:, 7], w, [0.16, 0.5, 0.84]).tolist(),
            "e_16_50_84": weighted_quantiles(th[:, 8], w, [0.16, 0.5, 0.84]).tolist()}  # fmt: skip


def merged_cell2(runs):
    """The runs merged into one (``dynesty.utils.merge_runs``), as a NestedResult."""
    from dynesty.utils import merge_runs

    from fork.inference.nested import NestedResult

    merged = merge_runs([r.results for r in runs])
    return NestedResult(runs[0].space, merged, sum(r.nlive for r in runs), runs[0].dlogz,
                        sum(r.wall_s for r in runs))  # fmt: skip


def interval_with_cell2(cell2, ctx, geo) -> dict:
    """The headline's partitioned posterior with its cell 2 replaced by ``cell2``: Galle-night
    longitude summarised as ``m6c.summarize_variant`` does (6000 derived samples)."""
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
    from fork.inference.profile import galle_lon_of_date_rad

    res = [NestedResult.load(m6c.Job("headline", k).path(".pkl")) for k in range(4)]
    res[HEADLINE_CELL] = cell2
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng(list(REPORT_RNG))
    eq = p.equal_weight(rng)
    idx = rng.choice(len(eq), size=min(m6c.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    space = res[0].space
    galle = np.array([math.degrees(galle_lon_of_date_rad(space.theta(v), ctx, geo))
                      for v in eq[idx]])  # fmt: skip
    tg = m6c.comparison_targets()
    targets = {
        "jpl": tg["jpl"]["galle_lon_ecl_of_date_deg"],
        "le_verrier_derived": tg["le_verrier"]["galle_lon_ecl_of_date_deg_derived"],
    }
    return {"logz": p.logz, "logzerr": p.logzerr, "shares": p.shares.tolist(),
            "galle_lon_ecl_of_date_deg": m6c.lon_summary(galle, targets)}  # fmt: skip


def summarize() -> None:
    from fork import ledger
    from fork.inference.nested import NestedResult
    from fork.inference.profile import GalleGeometry

    missing = [j.id for j in JOBS if not j.path(".pkl").exists()]
    if missing:
        raise SystemExit(f"not finished: {missing}")
    t0, c0 = time.time(), ledger.cpu_now()
    _, ctx, _ = cell2_problem()
    geo = GalleGeometry.from_horizons()
    runs = {j.id: NestedResult.load(j.path(".pkl")) for j in JOBS}
    # the two earlier runs of the same cell, for information (not in the decision)
    earlier = {f"{v}/c2": NestedResult.load(m6c.Job(v, HEADLINE_CELL).path(".pkl"))
               for v in ("headline", "s1_e05")}  # fmt: skip
    stats = {}
    for i, (name, r) in enumerate({**runs, **earlier}.items()):
        stats[name] = run_stats(r, ctx, geo, np.random.default_rng([5, 1846, i]))
        print(name, json.dumps({k: stats[name][k] for k in ("median", "width68", "logz")}))
    r1 = [j.id for j in JOBS if j.group == "r1"]
    r2 = [j.id for j in JOBS if j.group == "r2"]
    groups = {"r1": r1, "r2": r2, "r1_and_earlier": r1 + list(earlier)}
    spreads = {
        g: {q: spread([stats[n][q] for n in names], [stats[n][e] for n in names])
            for q, e in (("median", "median_err"), ("width68", "width68_err"),
                         ("logz", "logzerr"))}
        for g, names in groups.items()
    }  # fmt: skip
    decision = decide([stats[n]["width68"] for n in r1], [stats[n]["width68_err"] for n in r1],
                      [stats[n]["width68"] for n in r2])  # fmt: skip
    # R2 only is the replacement candidate; R1 only and R1 + R2 are reported for transparency
    merge_sets = {"r2": r2, "r1": r1, "r1_r2": r1 + r2}
    merged, merges = {}, {}
    for i, (name, ids) in enumerate(merge_sets.items()):
        merged[name] = merged_cell2([runs[n] for n in ids])
        merges[name] = {"runs": ids,
                        "cell2": run_stats(merged[name], ctx, geo,
                                           np.random.default_rng([5, 1846, 99, i])),
                        "interval": interval_with_cell2(merged[name], ctx, geo)}  # fmt: skip
    r2_interval = merges["r2"]["interval"]["galle_lon_ecl_of_date_deg"]
    head = json.loads((SITE_DATA / m6c.OUT_NAME).read_text())["fits"]["headline"]
    systematic = {"median_std_deg": spreads["r1"]["median"]["std"],
                  "width68_std_deg": spreads["r1"]["width68"]["std"]}  # fmt: skip
    if decision["outcome"] == "replace":
        reported = {"source": "headline cells 0, 1, 3 + merged R2 cell 2 (R1 not merged)",
                    **r2_interval}  # fmt: skip
    else:
        reported = {"source": "headline as published, R1 spread as a systematic",
                    **head["galle_lon_ecl_of_date_deg"], "systematic": systematic}  # fmt: skip
        if decision["outcome"] == "cross_check":
            reported["cross_check_not_replacement"] = {
                "source": "merged R2 cell 2",
                "cell2_median": merges["r2"]["cell2"]["median"],
                "cell2_width68": merges["r2"]["cell2"]["width68"],
                "interval_with_headline_cells_0_1_3": r2_interval,
            }
    payload = {
        "description": "M6c cell-2 repeatability (pre-registered in docs/notes/M6.md): R1, 4 "
        "runs at the headline settings; R2, 2 runs with nlive 4000 and 2x slices. Galle-night "
        "longitude of the cell-2 posterior per run, run-to-run spread against dynesty's own "
        "errors, merged posteriors (R2 only: the replacement candidate; R1 only and R1+R2: "
        "for transparency), and the reported interval by the decision rule.",
        "settings": {j.id: {"nlive": j.nlive, "slices": j.slices or DEFAULT_SLICES,
                            "seed": j.rng_seed(), "workers": j.workers} for j in JOBS},
        "rule": {"r2_width_tol_rel": R2_WIDTH_TOL_REL, "r1_chi2_p": R1_CHI2_P,
                 "n_boot": N_BOOT, "logl_window": LOGL_WINDOW},
        "runs": stats, "spread": spreads, "decision": decision,
        "merges": merges, "reported": reported,
    }  # fmt: skip
    write_json(OUT_NAME, payload)
    print(json.dumps({"decision": decision, "reported": reported}, indent=1))
    n = sum(int(np.sum(np.asarray(r.results["logl"]) > r.max_logl - LOGL_WINDOW))
            for r in [*runs.values(), *earlier.values(), *merged.values()])  # fmt: skip
    n += len(merged) * m6c.DERIVED_MAX_SAMPLES
    ledger.append("m6c_repeat/summarize", "summarize", n_likelihood_calls=0, n_integrations=n,
                  sim_years=n * (1846.73 - 1781.2), wall_s=time.time() - t0,
                  cpu_s=ledger.cpu_now() - c0, n_workers=1, estimated=True,
                  note="Galle-night longitude per weighted sample")  # fmt: skip


DERIVED_OUT_NAME = "m6c_replacement.json"


def replacement_cells() -> list:
    """The pre-registered replacement posterior's cells: the headline's cells 0, 1 and 3, and
    the merge of R2's runs as cell 2."""
    from fork.inference.nested import NestedResult

    res = [NestedResult.load(m6c.Job("headline", k).path(".pkl")) for k in range(4)]
    res[HEADLINE_CELL] = merged_cell2([NestedResult.load(j.path(".pkl"))
                                       for j in JOBS if j.group == "r2"])  # fmt: skip
    return res


def derived() -> None:
    """Mass, a, e, i, Peirce probability and the rest of the headline's summary record, from
    the replacement posterior (after the fact: the pre-registration fixed only its Galle-night
    interval). Same code as the headline record (``m6c.summarize_results``) with the rng of
    ``interval_with_cell2``, so the Galle-night numbers equal the reported interval."""
    from fork import ledger

    rep = json.loads((SITE_DATA / OUT_NAME).read_text())
    if rep["decision"]["outcome"] != "replace":
        raise SystemExit(f"outcome is {rep['decision']['outcome']}, not replace")
    m6c.check_frozen("headline")
    t0, c0 = time.time(), ledger.cpu_now()
    out = m6c.summarize_results("headline", replacement_cells(),
                                np.random.default_rng(list(REPORT_RNG)))  # fmt: skip
    reported = {k: v for k, v in rep["reported"].items() if k != "source"}
    if out["galle_lon_ecl_of_date_deg"] != reported:
        raise RuntimeError("Galle-night summary differs from the reported interval")
    r2 = rep["spread"]["r2"]["median"]
    payload = {
        "description": "M6c reported posterior after the pre-registered cell-2 repeatability "
        "study (outcome replace): the headline's cells 0, 1, 3 and the merged R2 cell 2. Every "
        "quantity other than the Galle-night 68 % interval is derived after the fact from this "
        "posterior with the headline's summary code; the 95 % bounds are provisional (they rest "
        "on the cell shares, and cells 0, 1, 3 are single runs).",
        "source": "headline cells 0, 1, 3 + merged R2 cell 2 (site/data/m6c_repeat.json)",
        "sampler_systematic": {
            "galle_median_deg": round(r2["std"] / math.sqrt(r2["n"]), 3),
            "basis": "R2 run-to-run std of the cell-2 Galle median (ddof 1) / sqrt(2 runs "
                     "merged); not folded into the interval",
            "r2_run_to_run_std_deg": round(r2["std"], 3),
        },
        "provisional": ["lo95", "hi95", "width95"],
        **out,
    }  # fmt: skip
    write_json(DERIVED_OUT_NAME, payload)
    keys = ("sampler_systematic", "galle_lon_ecl_of_date_deg", "peirce_p_within_0p87_deg_of_jpl",
            "mass_earth", "a_au", "e")  # fmt: skip
    print(json.dumps({k: payload[k] for k in keys}, indent=1))
    n = 2 * out["n_derived"]
    ledger.append("m6c_repeat/derived", "summarize", n_likelihood_calls=1, n_integrations=n,
                  sim_years=n * (1846.73 - 1781.2), wall_s=time.time() - t0,
                  cpu_s=ledger.cpu_now() - c0, n_workers=1, estimated=True,
                  note="replacement posterior: Galle, 1847 longitudes")  # fmt: skip


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="fork.export.m6c_repeat")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    p_run = sub.add_parser("run")
    p_run.add_argument("--n-workers", type=int, default=12)
    p_job = sub.add_parser("job")
    p_job.add_argument("job_id")
    p_job.add_argument("--n-workers", type=int, default=None)
    p_job.add_argument("--maxiter", type=int, default=None)
    sub.add_parser("summarize")
    sub.add_parser("derived")
    a = ap.parse_args(argv)
    if a.cmd == "plan":
        plan()
    elif a.cmd == "run":
        run_all(a.n_workers)
    elif a.cmd == "job":
        job = RepeatJob.from_id(a.job_id)
        do_job(job, a.n_workers or job.workers, a.maxiter)
    elif a.cmd == "derived":
        derived()
    else:
        summarize()


if __name__ == "__main__":
    main()
