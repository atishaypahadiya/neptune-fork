"""The low-distance solution by its own region, and the headline's 95 % bounds with eight seeds
(pre-registered in docs/notes/S1_mode.md; nothing here may differ).

    uv run python -m fork.export.s1_mode g0 [--n-workers 8]          # the sampler on a known answer
    uv run python -m fork.export.s1_mode pipeline [--n-workers 12]   # runs, reweighting, summary
    uv run python -m fork.export.s1_mode job-a REGION S [--maxiter N --cache DIR --synthetic SEED]
    uv run python -m fork.export.s1_mode reweight | summarize | plan

- **Part A:** S1's prior (e ≤ 0.5) in five regions: M (the low-distance solution), L (the rest
  of cell 0, two boxes L1 and L2) and cells 1 to 3. UltraNest's slice sampler on the fast
  likelihood of M8, eight seeds per region, and in M two runs at stronger settings (X). Each
  posterior reweighted to the full likelihood with Uranus marginalised (mode (b)).
- **Part B:** M8b's campaign with seeds 5 to 8, run by ``m8_runs job ... --campaign m8b`` into
  M8b's own cache; the eight-seed merge is drawn and reweighted here, into this plan's cache,
  so M8b's four-seed outputs are never rewritten.
- **Summary:** gates G0 to G4, R1, R2, X, rule B and the reading of § 6, into
  ``docs/notes/S1_mode.json``.

Implementation choices the plan does not fix, stated here and in the results:

- UltraNest draws from numpy's global generator, seeded with the plan's rng list
  (``np.random.seed([2026, 15, 1846, 1, r, s])``). A run resumed after a stop continues from
  its stored points with a fresh generator state, so it is not bit-identical to an
  uninterrupted one.
- A fast log-likelihood that is not finite is passed to UltraNest as −1e300.
- The eight seeds' merge of a region is the pooled posterior: 750 of each seed's 6,000
  reweighted draws (3,000 of each in X), so the merge needs no likelihood call of its own. Its
  ln Z is the log of the mean of the seeds' Z, as § 3 says; ESS/N, the median and the widths
  of the merge come from the pooled draws.
- The sampler error of P(a < 25 AU) is propagated to first order from each region's merged
  ln Z error (standard deviation over seeds / √n).
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
from pathlib import Path

import numpy as np

from fork import ledger
from fork.export import m6c, m8_runs, s1_cell0

ROOT = m6c.ROOT
CACHE = ROOT / "data" / "cache" / "nested" / "s1_mode"
OUT_JSON = ROOT / "docs" / "notes" / "S1_mode.json"
TAG = 15

# § 1: the inputs, fixed
INPUTS = {
    "data/historical/leverrier_1846.csv":
        "2f546b9979f613c62d9c2dc77cd92b8269a9dbea4d3f87eecf39d2145509ebc4",
    "data/historical/m6_error_model.json":
        "0856fa30d1675bf926d3921580de919b45fb72474e2f41ba58a5d153a6bec277",
    "docs/notes/M8_u_ref.json":
        "4bf4537c88ab2364c7f46247d87489de7beaa3031b677d8080435887caa4b3c5",
    "src/fork/model/perturbative.py":
        "092f673b9352f9d2b0ae709f2a4d3f258dd09806245b767fd83b86e8ecccd68c",
    "src/fork/model/linearized.py":
        "3e02e78bdd9ee74fb96a6666f192188c355113addc92f43e376729f0b616d15e",
    "src/fork/inference/reweight.py":
        "d0c5fbff6066ebaa531bc24668fa807a8f7bf0aa96cf77eead5c07aca9cd021b",
    "site/data/m6_fits.json":
        "2d111e7fcb7a71100a72a8f963f6c8f0f26aea9da90ff22a246893618bb74fff",
    "site/data/m8_runs.json":
        "4d50867272ab7ce6faec99ef3db1edb30e6b25ba61223f8adc0304c65f88b6cb",
    "site/data/m8b.json":
        "355104b24966ce216f53ff603931a060a4a5e738a996a683db5fac5c5e962027",
    "docs/notes/S1_cell0.json":
        "df29b5dc6c05882f319a341bea011b253858f3417d5cf2a9e48873acb7e327ec",
}  # fmt: skip

# § 2: the regions (a range AU, e range or None) and their prior masses as the plan states them
REGIONS = {
    "M": ((22.5, 25.0), (0.2, 0.5)),
    "L1": ((20.0, 22.5), (0.0, 0.5)),
    "L2": ((22.5, 25.0), (0.0, 0.2)),
    "c1": ((25.0, 28.0), None),
    "c2": ((28.0, 60.0), (0.0, 0.12)),
    "c3": ((28.0, 60.0), (0.12, 0.5)),
}
PLAN_MASS = {"M": 0.057542, "L": 0.145572, "c1": 0.103156, "c2": 0.166495, "c3": 0.527235}
CELL0 = ("M", "L1", "L2")
REGION_R = {"M": 1, "L1": 2, "L2": 3, "c1": 4, "c2": 5, "c3": 6}

# § 3: the method
SEEDS = tuple(range(1, 9))
X_SEEDS = (101, 102)
NLIVE, NSTEPS = 4000, 28
X_NLIVE, X_NSTEPS = 8000, 56
DLOGZ, FRAC_REMAIN = 0.1, 0.001
N_DRAWS = 6000
CHUNK = 250
BAD_LOGL = -1e300
PROGRESS_EVERY_S = 60.0

# § 4: gates and tests
G0_TOL, G0_RANGE = 0.3, 0.5
ESS_ALL, ESS_M = 0.3, 0.5
TRUNC_GATE = 1e-6
R1_TOL, R1_TOL_L, L_SHARE = 0.5, 1.0, 0.01
R2_MEDIAN_ERR, R2_WIDTH = 0.2, 0.10
X_TOL = {"ln_z": 0.5, "median": 0.35, "width_rel": 0.10}
AT_TOL = {"ln_z": 0.02, "median": 0.02, "width_rel": 0.005}

# § 5: Part B
B_SEEDS = tuple(range(1, 9))
B_NEW = (5, 6, 7, 8)
B_CELLS = range(4)
RULE_B_ERR = 0.2
M6C_95 = (320.18, 330.40)
B_MERGE_RNG = 8  # appended to M8b's rng for the eight-seed merge's draws and fit

# § 6: S1's reading
BANDS = ((0.05, "not at all"), (0.5, "moved partly"), (math.inf, "moved"))
PUBLISHED_WORD = "not at all"


def _sha(rel: str) -> str:
    return hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()


def check_inputs() -> None:
    bad = [rel for rel, h in INPUTS.items() if _sha(rel) != h]
    if bad:
        raise SystemExit(f"inputs differ from the plan (docs/notes/S1_mode.md § 1): {bad}")


def _log(msg: str, cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with (cache / "driver.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def _real(cache: Path) -> bool:
    return cache.resolve() == CACHE.resolve()


def _guard(cache: Path, maxiter: int | None, synthetic: int | None) -> None:
    if _real(cache) and (maxiter is not None or synthetic is not None):
        raise SystemExit("smoke runs (--maxiter, --synthetic) go to a scratch --cache")


# ---------------------------------------------------------------------------- the regions


def region_space(name: str, synthetic: int | None = None):
    from fork.inference.nested import cell_space

    return cell_space(s1_cell0.problem(synthetic)[2], REGIONS[name])


def prior_masses(synthetic: int | None = None) -> dict[str, float]:
    return {n: float(region_space(n, synthetic).slice_weight) for n in REGIONS}


def check_masses(w: dict[str, float]) -> None:
    got = {"M": w["M"], "L": w["L1"] + w["L2"], "c1": w["c1"], "c2": w["c2"], "c3": w["c3"]}
    bad = {k: v for k, v in got.items() if abs(v - PLAN_MASS[k]) > 5e-7}
    if bad or abs(sum(w.values()) - 1.0) > 1e-9:
        raise SystemExit(f"prior masses differ from the plan's § 2: {bad}")


def rng_list(name: str, s: int) -> list[int]:
    return [2026, TAG, 1846, 1, REGION_R[name], s]


def settings(s: int) -> tuple[int, int]:
    return (X_NLIVE, X_NSTEPS) if s in X_SEEDS else (NLIVE, NSTEPS)


def runs() -> list[tuple[str, int]]:
    return [(n, s) for n in REGIONS for s in SEEDS] + [("M", s) for s in X_SEEDS]


def run_dir(name: str, s: int, cache: Path) -> Path:
    return cache / "a" / name / f"s{s}"


def pid(name: str, s: int) -> str:
    return f"{name}__s{s}"


# ---------------------------------------------------------------------------- G0


G0_B, G0_S1, G0_S2, G0_SG, G0_DIM = 12.0, 0.06, 0.01, 0.03, 7


def g0_logl(x: np.ndarray) -> float:
    """The twisted Gaussian of § 4, normalised in R^7 (so ln Z on the cube is ln of the mass
    inside it)."""
    m = 0.5 + G0_B * ((x[0] - 0.5) ** 2 - G0_S1**2)
    r = np.concatenate([[(x[0] - 0.5) / G0_S1, (x[1] - m) / G0_S2], (x[2:] - 0.5) / G0_SG])
    return float(-0.5 * r @ r - math.log(G0_S1 * G0_S2 * G0_SG**5)
                 - 0.5 * G0_DIM * math.log(2 * math.pi))  # fmt: skip


def g0_true_ln_z() -> float:
    """ln of the density's mass in the unit cube: x₂ given x₁ in closed form, x₁ by quadrature."""
    from scipy.integrate import quad
    from scipy.stats import norm

    def f(x1: float) -> float:
        m = 0.5 + G0_B * ((x1 - 0.5) ** 2 - G0_S1**2)
        return norm.pdf(x1, 0.5, G0_S1) * (norm.cdf((1 - m) / G0_S2) - norm.cdf(-m / G0_S2))

    inner, _ = quad(f, 0.0, 1.0, points=[0.5], epsabs=1e-14, epsrel=1e-12, limit=200)
    side = norm.cdf(0.5 / G0_SG) - norm.cdf(-0.5 / G0_SG)
    return math.log(inner) + (G0_DIM - 2) * math.log(side)


def _g0_seed(s: int) -> dict:
    import tempfile

    import ultranest
    import ultranest.stepsampler as ss

    np.random.seed([2026, TAG, 1846, 0, 0, s])  # noqa: NPY002 (UltraNest draws from it)
    t0 = time.process_time()
    with tempfile.TemporaryDirectory() as d:
        sampler = ultranest.ReactiveNestedSampler(
            [f"x{i}" for i in range(G0_DIM)], g0_logl, lambda u: u, log_dir=d,
            resume="overwrite")  # fmt: skip
        sampler.stepsampler = ss.SliceSampler(
            nsteps=NSTEPS, generate_direction=ss.generate_mixture_random_direction
        )
        r = sampler.run(min_num_live_points=NLIVE, dlogz=DLOGZ, frac_remain=FRAC_REMAIN,
                        max_num_improvement_loops=0, show_status=False,
                        viz_callback=False)  # fmt: skip
    return {"seed": s, "ln_z": float(r["logz"]), "ln_z_err": float(r["logzerr"]),
            "ncall": int(r["ncall"]), "niter": int(r["niter"]),
            "cpu_s": time.process_time() - t0}  # fmt: skip


def g0_verdict(ln_z: list[float], truth: float) -> dict:
    dev = [z - truth for z in ln_z]
    rng_ = max(ln_z) - min(ln_z)
    ok = all(abs(d) <= G0_TOL for d in dev) and rng_ <= G0_RANGE
    return {"truth": truth, "deviations": dev, "range": rng_, "tol": G0_TOL,
            "range_tol": G0_RANGE, "pass": bool(ok)}  # fmt: skip


@ledger.recorded
def g0(n_workers: int = 8, cache: Path = CACHE) -> dict:
    from concurrent.futures import ProcessPoolExecutor

    rec = ledger.start("s1_mode/g0", "gate", n_workers=n_workers, synthetic=True,
                       scratch=not _real(cache))  # fmt: skip
    truth = g0_true_ln_z()
    _log(f"G0: true ln Z {truth:.4f}; eight seeds on {n_workers} workers", cache)
    with ProcessPoolExecutor(n_workers) as ex:
        rows = list(ex.map(_g0_seed, SEEDS))
    out = {"seeds": rows, **g0_verdict([r["ln_z"] for r in rows], truth),
           "settings": {"nlive": NLIVE, "nsteps": NSTEPS, "dlogz": DLOGZ,
                        "frac_remain": FRAC_REMAIN}}  # fmt: skip
    (cache / "g0.json").write_text(json.dumps(out, indent=1))
    _log(f"G0: range {out['range']:.3f}, largest |dev| "
         f"{max(abs(d) for d in out['deviations']):.3f}: {'pass' if out['pass'] else 'FAIL'}",
         cache)  # fmt: skip
    rec.done("gate passed" if out["pass"] else "gate failed",
             n_likelihood_calls=sum(r["ncall"] for r in rows), n_integrations=0,
             sim_years=0.0)  # fmt: skip
    return out


# ---------------------------------------------------------------------------- Part A runs

_S: dict = {"calls": 0}


def _transform(u: np.ndarray) -> np.ndarray:
    return _S["space"].prior_transform(np.concatenate([np.full(4, 0.5), u]))[4:]


def _logl(phi: np.ndarray) -> float:
    _S["calls"] += 1
    now = time.time()
    if now - _S["t_log"] >= PROGRESS_EVERY_S:
        _S["t_log"] = now
        _S["log"](f"{_S['calls']} calls")
    v = _S["pr"].log_likelihood(phi, second_order=m8_runs.SECOND_ORDER)
    return float(v) if math.isfinite(v) else BAD_LOGL


@ledger.recorded
def do_job_a(name: str, s: int, cache: Path = CACHE, maxiter: int | None = None,
             synthetic: int | None = None) -> None:  # fmt: skip
    import ultranest
    import ultranest.stepsampler as ss

    _guard(cache, maxiter, synthetic)
    if (name, s) not in runs():
        raise SystemExit(f"no run {name}/s{s} in the plan")
    d = run_dir(name, s, cache)
    out = d / "result.npz"
    if out.exists():
        print(f"a/{name}/s{s}: cached")
        return
    if synthetic is None:
        check_inputs()
    *_, pr = s1_cell0.problem(synthetic)
    cs = region_space(name, synthetic)
    periodic = [i - 4 for i in cs.periodic_unit_dims]
    if min(periodic) < 0 or len(cs.unit_names) - 4 != m8_runs.N_PHI:
        raise RuntimeError("unexpected prior layout")
    nlive, nsteps = settings(s)
    store = d / "un"
    resumed = (store / "results").exists() or (store / "chains").exists()
    d.mkdir(parents=True, exist_ok=True)

    def log(msg: str) -> None:
        with (cache / "progress.log").open("a") as fh:
            fh.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} s1_mode/a/{name}/s{s}: {msg}\n")

    _S.update(space=cs, pr=pr, calls=0, t_log=time.time(), log=log)
    rec = ledger.start(f"s1_mode/a/{name}/s{s}", "nested", n_workers=1, region=name, seed=s,
                       nlive=nlive, nsteps=nsteps, sampler="ultranest slice", resumed=resumed,
                       capped=maxiter is not None, synthetic=synthetic is not None,
                       scratch=not _real(cache), likelihood="perturbative, second order off",
                       n_integrations=0, sim_years=0.0)  # fmt: skip
    rec.set(n_likelihood_calls=lambda: _S["calls"])
    np.random.seed(rng_list(name, s))  # noqa: NPY002 (UltraNest draws from it)
    log(f"{'resume' if resumed else 'start'}: nlive {nlive}, nsteps {nsteps}, rng "
        f"{rng_list(name, s)}, region {REGIONS[name]}")  # fmt: skip
    wrapped = [i in periodic for i in range(m8_runs.N_PHI)]
    sampler = ultranest.ReactiveNestedSampler(
        list(cs.unit_names[4:]), _logl, _transform, wrapped_params=wrapped,
        log_dir=str(store), resume="resume")  # fmt: skip
    sampler.stepsampler = ss.SliceSampler(
        nsteps=nsteps, generate_direction=ss.generate_mixture_random_direction
    )
    r = sampler.run(min_num_live_points=nlive, dlogz=DLOGZ, frac_remain=FRAC_REMAIN,
                    max_num_improvement_loops=0, max_iters=maxiter, show_status=False,
                    viz_callback=False)  # fmt: skip
    ws = r["weighted_samples"]
    tmp = d / "result.tmp.npz"
    np.savez(tmp, points=ws["points"], weights=ws["weights"], logl=ws["logl"],
             ln_z=r["logz"], ln_z_err=r["logzerr"], ncall=r["ncall"], niter=r["niter"],
             capped=maxiter is not None, nlive=nlive, nsteps=nsteps)  # fmt: skip
    tmp.replace(out)
    log(f"done: ln Z {r['logz']:.3f} ± {r['logzerr']:.3f}, {r['niter']} it, {r['ncall']} calls")
    rec.done(niter=int(r["niter"]), ncall_total=int(r["ncall"]))


# ---------------------------------------------------------------------------- scheduling


def _spawn(cmd: list[str], logfile: Path) -> subprocess.Popen:
    logfile.parent.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1"}
    with logfile.open("a") as fh:
        return subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT, env=env)


def schedule(jobs: list[tuple[str, Path, list[str]]], n_slots: int, cache: Path) -> bool:
    """Run ``jobs`` = (name, done path, command), one slot each, at most ``n_slots`` at once;
    one retry per failed job (it resumes from its own store)."""
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
        while pending and len(running) < n_slots:
            job = pending.pop(0)
            proc = _spawn(job[2], cache / "logs" / (job[0].replace("/", "__") + ".log"))
            running[job[0]] = (proc, time.time(), job)
            _log(f"start {job[0]} (pid {proc.pid})", cache)
        time.sleep(30)
    ok = all(j[1].exists() for j in jobs)
    _log(f"schedule finished; failures: {failures or 'none'}; all done: {ok}", cache)
    return ok


def b_cache() -> Path:
    return m8_runs.CAMPAIGNS["m8b"].cache


def b_job_path(k: int, s: int, cache: Path | None = None) -> Path:
    return m8_runs.RunJob("headline", k, s, "m8b").path(".pkl", cache or b_cache())


def all_jobs(cache: Path, extra_a: list[str], b_cache_dir: Path | None, extra_b: list[str]):
    """Part A's 50 runs (X's two first, then M, then the rest) and Part B's 16, interleaved so
    the long cell-0 runs of Part B start early."""
    py = [sys.executable, "-m"]
    a = [(f"a/{n}/s{s}", run_dir(n, s, cache) / "result.npz",
          [*py, "fork.export.s1_mode", "job-a", n, str(s), "--cache", str(cache), *extra_a])
         for n, s in sorted(runs(), key=lambda r: (r[1] not in X_SEEDS, r[0] != "M"))]  # fmt: skip
    bc = b_cache_dir or b_cache()
    b = [(f"b/c{k}/s{s}", b_job_path(k, s, bc),
          [*py, "fork.export.m8_runs", "job", f"headline/c{k}/s{s}", "--campaign", "m8b",
           "--cache", str(bc), *extra_b])
         for k in B_CELLS for s in B_NEW]  # fmt: skip
    out = []
    while a or b:
        if b:
            out.append(b.pop(0))
        out.extend(a[:3])
        del a[:3]
    return out


# ---------------------------------------------------------------------------- reweighting


def draws_a(name: str, s: int, cache: Path) -> dict:
    """6,000 equal-weight draws of one run's posterior (systematic resampling of its weighted
    samples, then a choice without replacement, as for the dynesty runs), cached."""
    from dynesty.utils import resample_equal

    path = cache / "reweight" / pid(name, s) / "draws.npz"
    if path.exists():
        d = np.load(path)
        return {k: d[k] for k in d.files}
    r = np.load(run_dir(name, s, cache) / "result.npz")
    rng = np.random.default_rng([*rng_list(name, s), 7])
    w = r["weights"] / r["weights"].sum()
    eq = resample_equal(r["points"], w, rstate=rng)
    n = min(N_DRAWS, len(eq))
    d = {"phi": eq[rng.choice(len(eq), n, replace=False)], "ln_z_fast": r["ln_z"],
         "ln_z_fast_err": r["ln_z_err"], "ncall": r["ncall"], "niter": r["niter"],
         "capped": r["capped"]}  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("draws.tmp.npz")
    np.savez(tmp, **d)
    tmp.replace(path)
    return d


_RW: dict = {}


def _init_reweight(synthetic: int | None) -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from fork.inference.profile import GalleGeometry

    _, ctx, space, lp, pr = s1_cell0.problem(synthetic)
    _RW.update(ctx=ctx, space=space, lp=lp, pr=pr, geo=GalleGeometry.from_horizons())


def _reweight_chunk_a(task: tuple) -> tuple[str, int, float]:
    """Both likelihoods at draws [lo, hi) of one run, and each draw's Galle-night longitude
    through the full model with u drawn from mode (b)'s N(u*, Σ) (as ``s1_cell0``)."""
    from fork.inference.profile import galle_lon_of_date_rad

    name, s, c, lo, hi, cache = task
    out = Path(cache) / "reweight" / pid(name, s) / f"chunk_{c:03d}.npz"
    if out.exists():
        return pid(name, s), c, 0.0
    t0 = time.process_time()
    ctx, space, lp, pr, geo = (_RW[x] for x in ("ctx", "space", "lp", "pr", "geo"))
    phi = np.load(out.parent / "draws.npz")["phi"][lo:hi]
    n = len(phi)
    rng = np.random.default_rng([*rng_list(name, s), 9, c])
    a = {"lnl_fast": np.full(n, -np.inf), "lnl_b": np.full(n, -np.inf),
         "b_trunc_m2lnl": np.full(n, np.nan), "galle_deg": np.full(n, np.nan)}  # fmt: skip
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
    return pid(name, s), c, time.process_time() - t0


def b_merge_draws(k: int, cache: Path, max_draws: int | None = None) -> Path:
    """The eight-seed merge of Part B's cell k (dynesty's merge_runs of seeds 1 to 8), its
    6,000 draws written to this plan's cache, laid out as ``m8_runs`` lays out a merge."""
    from dynesty.utils import merge_runs

    from fork.inference.nested import NestedResult
    from fork.inference.reweight import N_PER_CELL

    bdir = cache / "b8"
    path = bdir / "reweight" / m8_runs.pid("headline", k, m8_runs.MERGED) / "draws.npz"
    if path.exists():
        return bdir
    rs = [NestedResult.load(b_job_path(k, s, _b_src(cache))) for s in B_SEEDS]
    merged = merge_runs([r.results for r in rs], print_progress=False)
    res = NestedResult(rs[0].space, merged, sum(r.nlive for r in rs), rs[0].dlogz,
                       sum(r.wall_s for r in rs))  # fmt: skip
    rng = np.random.default_rng([2026, 8, 1846, 0, k + 1, 0, 7, 3, B_MERGE_RNG])
    eq = res.equal_weight(rng)
    n = min(N_PER_CELL if max_draws is None else max_draws, len(eq))
    d = {"phi": eq[rng.choice(len(eq), n, replace=False)], "ln_z_fast": np.array(res.logz),
         "ln_z_fast_err": np.array(res.logzerr), "ncall": np.array(res.ncall),
         "niter": np.array(res.niter)}  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("draws.tmp.npz")
    np.savez(tmp, **d)
    tmp.replace(path)
    return bdir


def _b_src(cache: Path) -> Path:
    """Where Part B's runs are: M8b's cache for the real run, ``cache/b`` for a smoke run."""
    return b_cache() if _real(cache) else cache / "b"


@ledger.recorded
def reweight(n_workers: int, cache: Path = CACHE, synthetic: int | None = None,
             max_draws: int | None = None) -> None:  # fmt: skip
    from concurrent.futures import ProcessPoolExecutor, as_completed

    missing = [f"a/{n}/s{s}" for n, s in runs()
               if not (run_dir(n, s, cache) / "result.npz").exists()]  # fmt: skip
    bsrc = _b_src(cache)
    missing += [f"b/c{k}/s{s}" for k in B_CELLS for s in B_SEEDS
                if not b_job_path(k, s, bsrc).exists()]  # fmt: skip
    if missing:
        raise SystemExit(f"runs missing: {missing}")
    rec = ledger.start("s1_mode/reweight", "reweight", n_workers=n_workers,
                       scratch=not _real(cache), synthetic=synthetic is not None)  # fmt: skip
    ta, tb = [], []
    for n, s in runs():
        nd = len(draws_a(n, s, cache)["phi"])
        nd = nd if max_draws is None else min(nd, max_draws)
        for c, lo in enumerate(range(0, nd, CHUNK)):
            if not (cache / "reweight" / pid(n, s) / f"chunk_{c:03d}.npz").exists():
                ta.append((n, s, c, lo, min(lo + CHUNK, nd), str(cache)))
    for k in B_CELLS:
        for s in B_NEW:
            nd = len(m8_runs.draws("headline", k, s, bsrc, max_draws, "m8b")["phi"])
            for c, lo in enumerate(range(0, nd, m8_runs.REWEIGHT_CHUNK)):
                p = bsrc / "reweight" / m8_runs.pid("headline", k, s) / f"chunk_{c:03d}.npz"
                if not p.exists():
                    tb.append(("headline", k, s, c, lo, min(lo + m8_runs.REWEIGHT_CHUNK, nd),
                               str(bsrc)))  # fmt: skip
        bdir = b_merge_draws(k, cache, max_draws)
        mp = m8_runs.pid("headline", k, m8_runs.MERGED)
        nd = len(np.load(bdir / "reweight" / mp / "draws.npz")["phi"])
        for c, lo in enumerate(range(0, nd, m8_runs.REWEIGHT_CHUNK)):
            if not (bdir / "reweight" / mp / f"chunk_{c:03d}.npz").exists():
                tb.append(("headline", k, m8_runs.MERGED, c, lo,
                           min(lo + m8_runs.REWEIGHT_CHUNK, nd), str(bdir)))  # fmt: skip
    if not ta and not tb:
        rec.discard()
        return
    _log(f"reweight: {len(ta)} chunks of Part A, {len(tb)} of Part B, {n_workers} workers",
         cache)  # fmt: skip
    if synthetic is not None and tb:
        raise SystemExit("Part B has no synthetic data set; smoke runs reweight Part A only")
    done = 0
    with ProcessPoolExecutor(n_workers, initializer=_init_reweight,
                             initargs=(synthetic,)) as ex:  # fmt: skip
        futs = [ex.submit(_reweight_chunk_a, t) for t in ta]
        for f in as_completed(futs):
            f.result()
            done += 1
            if done % 48 == 0 or done == len(futs):
                _log(f"reweight A: {done}/{len(futs)} chunks", cache)
    if tb:
        with ProcessPoolExecutor(n_workers, initializer=m8_runs._init_reweight) as ex:
            futs = [ex.submit(m8_runs._reweight_chunk, t) for t in tb]
            for i, f in enumerate(as_completed(futs)):
                f.result()
                if (i + 1) % 48 == 0 or i + 1 == len(futs):
                    _log(f"reweight B: {i + 1}/{len(futs)} chunks", cache)
    rec.done(n_likelihood_calls=sum(t[4] - t[3] for t in ta) + sum(t[5] - t[4] for t in tb),
             n_integrations=None, sim_years=None,
             note="mode (b) + fast + Galle longitude per draw")  # fmt: skip


def load_a(name: str, s: int, cache: Path) -> dict:
    d = draws_a(name, s, cache)
    n = len(d["phi"])
    parts = [np.load(cache / "reweight" / pid(name, s) / f"chunk_{c:03d}.npz")
             for c in range(math.ceil(n / CHUNK))]  # fmt: skip
    if sum(len(p["lnl_fast"]) for p in parts) < n:  # a smoke run reweights fewer draws
        n = sum(len(p["lnl_fast"]) for p in parts)
    out = {k: d[k] for k in d}
    out["phi"] = d["phi"][:n]
    for key in parts[0].files:
        out[key] = np.concatenate([p[key] for p in parts])[:n]
    return out


# ---------------------------------------------------------------------------- statistics


def circ_dev(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (x - x[0] + 180.0) % 360.0 - 180.0


def spread(x, circular: bool = False) -> float:
    v = circ_dev(x) if circular else np.asarray(x, dtype=float)
    return float(v.max() - v.min())


def seed_error(x, circular: bool = False) -> float:
    """Sampler error of a merge from n seeds: standard deviation over seeds / √n."""
    v = circ_dev(x) if circular else np.asarray(x, dtype=float)
    return float(np.std(v, ddof=1) / math.sqrt(v.size))


def log_mean_exp(x) -> float:
    from scipy.special import logsumexp

    x = np.asarray(x, dtype=float)
    return float(logsumexp(x) - math.log(x.size))


def stats(d: dict, rng: np.random.Generator, space=None) -> dict:
    """One reweighted posterior: ln Z, ESS/N, truncation, the Galle-night summary, and with a
    space the perturber's a and e (16, 50, 84 %)."""
    from fork.inference.reweight import reweight_cell

    cr = reweight_cell(d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]),
                       float(d["ln_z_fast_err"]))  # fmt: skip
    _, ix = m8_runs._resampled([cr], [0.0], rng, cr.n)
    g = m6c.lon_summary(d["galle_deg"][ix], {})
    fin = np.isfinite(d["lnl_b"])
    out = {"n": cr.n, "ess_frac": cr.ess / cr.n, "ln_z_fast": cr.ln_z_fast,
           "ln_z_fast_err": cr.ln_z_fast_err, "ln_correction": cr.ln_correction,
           "ln_correction_err": cr.ln_correction_err, "ln_z": cr.ln_z,
           "galle": g, "galle_median": g["median"], "galle_width68": g["width68"],
           "b_trunc_m2lnl_max": float(np.nanmax(d["b_trunc_m2lnl"])) if fin.any() else None,
           "n_b_failed": int(np.sum(~fin))}  # fmt: skip
    if space is not None:
        p = m8_runs._perturber(space, d["phi"][ix])
        out["a_au"] = np.percentile(p[:, 1], [16, 50, 84]).tolist()
        out["e"] = np.percentile(p[:, 2], [16, 50, 84]).tolist()
    return out


def pooled(ds: list[dict], per: int, rng: np.random.Generator) -> dict:
    """The merge as the pooled posterior: ``per`` draws of each run, chosen without
    replacement; ln Z_fast as the log of the mean of the runs' Z_fast."""
    keys = ("phi", "lnl_fast", "lnl_b", "b_trunc_m2lnl", "galle_deg")
    pick = [rng.choice(len(d["lnl_fast"]), min(per, len(d["lnl_fast"])), replace=False)
            for d in ds]  # fmt: skip
    out = {k: np.concatenate([d[k][ix] for d, ix in zip(ds, pick, strict=True)]) for k in keys}
    zf = [float(d["ln_z_fast"]) for d in ds]
    out["ln_z_fast"] = log_mean_exp(zf)
    out["ln_z_fast_err"] = seed_error(zf)
    return out


def word(p: float) -> str:
    return next(w for lim, w in BANDS if p < lim)


def verdict_word(p: float, err: float) -> str:
    lo, hi = word(max(p - 2 * err, 0.0)), word(min(p + 2 * err, 1.0))
    return lo if lo == hi else f"between {lo} and {hi}"


def p_below_25(lnz: dict[str, float], err: dict[str, float], w: dict[str, float]) -> dict:
    """P(a < 25 AU) = cell 0's share over the six regions, with its first-order sampler error;
    and ln Z₀ of cell 0 (M and L's two boxes) with its error."""
    from fork.export import v13

    names = list(REGIONS)
    t = np.array([math.log(w[n]) + lnz[n] for n in names])
    share = np.exp(t - (t.max() + math.log(np.exp(t - t.max()).sum())))
    s = dict(zip(names, share, strict=True))
    p = sum(s[n] for n in CELL0)
    grad = {n: s[n] * (1 - p) if n in CELL0 else -p * s[n] for n in names}
    p_err = math.sqrt(sum((grad[n] * err[n]) ** 2 for n in names))
    ln_z0, ln_z0_err = v13.combine_subcells(np.array([lnz[n] for n in CELL0]),
                                            np.array([err[n] for n in CELL0]),
                                            np.array([w[n] for n in CELL0]))  # fmt: skip
    l_share = (s["L1"] + s["L2"]) / p
    return {"p": float(p), "p_err": p_err, "shares": {n: float(v) for n, v in s.items()},
            "ln_z0": float(ln_z0), "ln_z0_err": float(ln_z0_err),
            "l_share_of_cell0": float(l_share)}  # fmt: skip


def at_tol(value: float, tol: float, band: float) -> bool:
    return abs(abs(value) - tol) <= band


# ---------------------------------------------------------------------------- Part A summary


def part_a(cache: Path, synthetic: int | None = None) -> dict:
    from fork.export import m8_arbiter

    w = prior_masses(synthetic)
    space = s1_cell0.problem(synthetic)[2]
    reg: dict = {}
    for n in REGIONS:
        seeds = {s: load_a(n, s, cache) for s in SEEDS}
        st = {f"s{s}": stats(d, np.random.default_rng([*rng_list(n, s), 13]),
                             space if n == "M" else None)
              for s, d in seeds.items()}  # fmt: skip
        z = [st[f"s{s}"]["ln_z"] for s in SEEDS]
        pm = pooled(list(seeds.values()), N_DRAWS // len(SEEDS),
                    np.random.default_rng([*rng_list(n, 0), 21]))  # fmt: skip
        merged = stats(pm, np.random.default_rng([*rng_list(n, 0), 13]),
                       space if n == "M" else None)  # fmt: skip
        merged.update(ln_z_pooled=merged["ln_z"], ln_z=log_mean_exp(z), ln_z_err=seed_error(z))
        meds = [st[f"s{s}"]["galle_median"] for s in SEEDS]
        wds = np.array([st[f"s{s}"]["galle_width68"] for s in SEEDS])
        reg[n] = {"a_au": list(REGIONS[n][0]), "e": REGIONS[n][1], "prior_mass": w[n],
                  "seeds": st, "merged": merged, "ln_z_range": spread(z),
                  "median_error": seed_error(meds, circular=True),
                  "median_range": spread(meds, circular=True),
                  "width_rel_dev": (wds / wds.mean() - 1.0).tolist(),
                  "capped": [bool(seeds[s]["capped"]) for s in SEEDS],
                  "ncall": [int(seeds[s]["ncall"]) for s in SEEDS]}  # fmt: skip
    xs = {s: load_a("M", s, cache) for s in X_SEEDS}
    xst = {f"s{s}": stats(d, np.random.default_rng([*rng_list("M", s), 13]), space)
           for s, d in xs.items()}  # fmt: skip
    xz = [xst[f"s{s}"]["ln_z"] for s in X_SEEDS]
    xp = pooled(list(xs.values()), N_DRAWS // len(X_SEEDS),
                np.random.default_rng([*rng_list("M", 100), 21]))  # fmt: skip
    xm = stats(xp, np.random.default_rng([*rng_list("M", 100), 13]), space)
    xm.update(ln_z_pooled=xm["ln_z"], ln_z=log_mean_exp(xz), ln_z_err=seed_error(xz))

    lnz = {n: reg[n]["merged"]["ln_z"] for n in REGIONS}
    err = {n: reg[n]["merged"]["ln_z_err"] for n in REGIONS}
    p = p_below_25(lnz, err, w)

    posts = [(f"{n}/{k}", v) for n in REGIONS for k, v in
             [*reg[n]["seeds"].items(), ("merged", reg[n]["merged"])]]  # fmt: skip
    posts += [(f"X/{k}", v) for k, v in [*xst.items(), ("merged", xm)]]
    g1 = [f"{n}/s{s}" for n in REGIONS for s, c in zip(SEEDS, reg[n]["capped"], strict=True)
          if c] + [f"X/s{s}" for s in X_SEEDS if bool(xs[s]["capped"])]  # fmt: skip
    g2 = [{"posterior": k, "ess_frac": v["ess_frac"]} for k, v in posts
          if v["ess_frac"] < ESS_ALL]  # fmt: skip
    g3 = [{"posterior": k, "max": v["b_trunc_m2lnl_max"]} for k, v in posts
          if v["b_trunc_m2lnl_max"] is None or v["b_trunc_m2lnl_max"] >= TRUNC_GATE]  # fmt: skip
    g4_rows = [{"posterior": k, "ess_frac": v["ess_frac"]} for k, v in posts
               if k.startswith("M/") or k == "X/merged"]  # fmt: skip
    gates = {
        "G1": {"capped_runs": g1, "pass": not g1},
        "G2": {"threshold": ESS_ALL, "min_ess_frac": min(v["ess_frac"] for _, v in posts),
               "failures": g2, "pass": not g2},
        "G3": {"threshold": TRUNC_GATE, "failures": g3, "pass": not g3},
        "G4": {"threshold": ESS_M, "rows": g4_rows,
               "pass": all(r["ess_frac"] >= ESS_M for r in g4_rows)},
    }  # fmt: skip
    l_tol = R1_TOL_L if p["l_share_of_cell0"] < L_SHARE else R1_TOL
    r1 = {n: {"range": reg[n]["ln_z_range"], "tol": l_tol if n in ("L1", "L2") else R1_TOL}
          for n in REGIONS}  # fmt: skip
    for v in r1.values():
        v["pass"] = v["range"] <= v["tol"]
        v["at_tolerance"] = at_tol(v["range"], v["tol"], AT_TOL["ln_z"])
    m = reg["M"]
    r2 = {"median_error": m["median_error"], "median_tol": R2_MEDIAN_ERR,
          "median_pass": m["median_error"] <= R2_MEDIAN_ERR,
          "median_at_tolerance": at_tol(m["median_error"], R2_MEDIAN_ERR, AT_TOL["median"]),
          "width_rel_dev_max": float(np.max(np.abs(m["width_rel_dev"]))),
          "width_tol": R2_WIDTH}  # fmt: skip
    r2["width_pass"] = r2["width_rel_dev_max"] <= R2_WIDTH
    r2["width_at_tolerance"] = at_tol(r2["width_rel_dev_max"], R2_WIDTH, AT_TOL["width_rel"])
    r2["pass"] = r2["median_pass"] and r2["width_pass"]
    mm = m["merged"]
    x = {"d_ln_z": xm["ln_z"] - mm["ln_z"],
         "d_median": float(circ_dev([mm["galle_median"], xm["galle_median"]])[1]),
         "width_rel": xm["galle_width68"] / mm["galle_width68"] - 1.0}  # fmt: skip
    x["pass"] = bool(abs(x["d_ln_z"]) <= X_TOL["ln_z"] and abs(x["d_median"]) <= X_TOL["median"]
                     and abs(x["width_rel"]) <= X_TOL["width_rel"])  # fmt: skip
    x["at_tolerance"] = {k: at_tol(x[f"d_{k}" if k != "width_rel" else k], X_TOL[k], AT_TOL[k])
                         for k in X_TOL}  # fmt: skip

    checks = {}
    for key, other in (("published", s1_cell0.OTHER_PUBLISHED), ("fast", s1_cell0.OTHER_FAST)):
        z0, want, tol = s1_cell0.CHECK_VALUES[key]
        got = m8_arbiter.implied_p_below_25(z0, list(other))
        checks[key] = {"ln_z0": z0, "want": want, "got": got, "pass": abs(got - want) <= tol}
    beside = {key: m8_arbiter.implied_p_below_25(p["ln_z0"], list(other))
              for key, other in (("published", s1_cell0.OTHER_PUBLISHED),
                                 ("fast", s1_cell0.OTHER_FAST))}  # fmt: skip

    stop = not gates["G4"]["pass"]
    all_ok = (all(g["pass"] for g in gates.values()) and all(v["pass"] for v in r1.values())
              and r2["pass"] and x["pass"] and all(c["pass"] for c in checks.values()))  # fmt: skip
    g0 = _g0_result(cache)
    all_ok = all_ok and bool(g0 and g0["pass"])
    if stop:
        reading = {"status": "reading stopped (G4 failed)", "resolved": False}
    elif all_ok:
        reading = {"status": "resolved", "resolved": True, "p": p["p"], "p_err": p["p_err"],
                   "word": verdict_word(p["p"], p["p_err"]),
                   "erratum": verdict_word(p["p"], p["p_err"]) != PUBLISHED_WORD}  # fmt: skip
    else:
        reading = {"status": "not resolved", "resolved": False}
    return {
        "regions": reg,
        "X": {"seeds": xst, "merged": xm, "test": x},
        "G0": g0,
        "gates": gates,
        "R1": r1,
        "R2": r2,
        "p_below_25": p,
        "check_values": checks,
        "p_beside_not_read": beside,
        "reading": reading,
    }


def _g0_result(cache: Path) -> dict | None:
    f = cache / "g0.json"
    return json.loads(f.read_text()) if f.exists() else None


# ---------------------------------------------------------------------------- Part B summary


def load_b(k: int, s: int, cache: Path) -> dict:
    """A reweighted posterior of Part B: seeds 1 to 8 from M8b's cache, the eight-seed merge
    from this plan's."""
    if s == m8_runs.MERGED:
        return m8_runs.load_reweighted("headline", k, s, cache / "b8")
    return m8_runs.load_reweighted("headline", k, s, _b_src(cache))


def part_b(cache: Path, seeds: tuple[int, ...] = B_SEEDS, merged: bool = True,
           load=None) -> dict:  # fmt: skip
    """Rule B over ``seeds`` (the plan: 1 to 8, and their merge). With seeds 1 to 4 and
    M8b's own four-seed merge this reproduces ``site/data/m8b.json`` (a test)."""
    from fork.inference.reweight import N_PER_CELL, ess, reweight_cell

    load = load or (lambda k, s: load_b(k, s, cache))
    v, tag = "headline", m8_runs.CAMPAIGNS["m8b"].rng_tag
    tg = m6c.comparison_targets()
    targets = {
        "jpl": tg["jpl"]["galle_lon_ecl_of_date_deg"],
        "le_verrier_derived": tg["le_verrier"]["galle_lon_ecl_of_date_deg_derived"],
    }
    space = m8_runs.m8.problem(v)[2]
    names = [*seeds, m8_runs.MERGED] if merged else list(seeds)
    post, cells, g2, g3, g4 = {}, {}, [], [], []
    for k in B_CELLS:
        for s in names:
            d = load(k, s)
            cr = reweight_cell(d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]),
                               float(d["ln_z_fast_err"]))  # fmt: skip
            d["pert"], d["cr"] = m8_runs._perturber(space, d["phi"]), cr
            post[(k, s)] = d
            rng = np.random.default_rng([2026, 8, 1846, 0, k + 1, s, m8_runs.CELL_RNG_TAG, tag])
            _, ix = m8_runs._resampled([cr], [0.0], rng, cr.n)
            g = m8_runs._lon(d["galle_deg"][ix], targets)
            name = m8_runs.pid(v, k, s)
            cells[(k, s)] = {
                "ess_frac": cr.ess / cr.n,
                "ln_z": cr.ln_z,
                "galle_median": g["median"],
                "galle_width68": g["width68"],
            }
            g2.append({"posterior": name, "ess_frac": cr.ess / cr.n,
                       "pass": cr.ess / cr.n >= m8_runs.ESS_GATE_CELL})  # fmt: skip
            fin = np.isfinite(d["lnl_b"])
            g4.append({"posterior": name, "max": float(np.nanmax(d["b_trunc_m2lnl"])),
                       "pass": bool(np.all(d["b_trunc_m2lnl"][fin]
                                           < m8_runs.TRUNC_GATE_M2LNL))})  # fmt: skip
            if k == 2:
                a, e = d["pert"][:, 1], d["pert"][:, 2]
                for br, msk in (("inner", (e >= m8_runs.E_RIDGE) & (a < m8_runs.INNER_OUTER_A)),
                                ("outer", (e >= m8_runs.E_RIDGE)
                                 & (a >= m8_runs.INNER_OUTER_A))):  # fmt: skip
                    nb = int(msk.sum())
                    fr = ess(cr.delta[msk]) / nb if nb else math.nan
                    g3.append({"posterior": name, "branch": br, "n": nb, "ess_frac": fr,
                               "pass": bool(nb > 0 and fr >= 0.3)})  # fmt: skip
    fits, g1 = {}, []
    log_c = [math.log(m8_runs._cell_space(v, k).slice_weight) for k in B_CELLS]
    for s in names:
        cs = [post[(k, s)] for k in B_CELLS]
        rl = [2026, 8, 1846, 0, 0, s, m8_runs.FIT_RNG_TAG, tag]
        if s == m8_runs.MERGED and len(seeds) == len(B_SEEDS):
            rl.append(B_MERGE_RNG)
        res, ix = m8_runs._resampled([c["cr"] for c in cs], log_c, np.random.default_rng(rl),
                                     N_PER_CELL)  # fmt: skip
        galle = np.concatenate([c["galle_deg"] for c in cs])[ix]
        name = "merged" if s == m8_runs.MERGED else f"s{s}"
        fits[name] = {
            "ln_z": res.ln_z,
            "ln_z_err": res.ln_z_err,
            "galle_lon_ecl_of_date_deg": m8_runs._lon(galle, targets),
        }
        as_written, at0, share = m8_runs.fit_ess(res)
        g1.append({"fit": name, "ess_frac_as_written": as_written, "ess_frac_at_delta0": at0,
                   "reweighting_share": share, "threshold": m8_runs.ESS_GATE_FIT,
                   "pass": share >= m8_runs.ESS_GATE_FIT})  # fmt: skip

    def ok(rows: list[dict], key: str, keep) -> bool:
        return all(r["pass"] for r in rows if keep(r[key]))

    new = {f"s{s}" for s in B_NEW} | {"merged"}
    new_post = {m8_runs.pid(v, k, s) for k in B_CELLS for s in (*B_NEW, m8_runs.MERGED)}
    cond1 = ok(g1, "fit", lambda n: n in new) and all(
        ok(g, "posterior", lambda n: n in new_post) for g in (g2, g3, g4)
    )
    a8f = []
    for k in B_CELLS:
        wd = np.array([cells[(k, s)]["galle_width68"] for s in seeds])
        dev = np.abs(wd / wd.mean() - 1.0)
        pre = f"{v}__c{k}__"
        gates_ok = all(g["pass"] for g in g2 + g3 + g4 if g["posterior"].startswith(pre))
        a8f.append(
            {
                "cell": k,
                "widths": wd.tolist(),
                "max_rel_dev": float(dev.max()),
                "pass": bool(dev.max() <= m8_runs.TOL["a8f_width_rel"] and gates_ok),
            }
        )
    m8b = json.loads((ROOT / "site" / "data" / "m8b.json").read_text())
    single = [fits[f"s{s}"]["galle_lon_ecl_of_date_deg"] for s in seeds]
    std95 = {b: float(np.std([g[b] for g in single], ddof=1)) for b in ("lo95", "hi95")}
    err95 = {b: x / math.sqrt(len(seeds)) for b, x in std95.items()}
    cond = {
        "1_gates_new_seeds_and_merge": bool(cond1),
        "2_a8b_a8c": bool(m8b["tests"]["a8b"]["pass"] and m8b["tests"]["a8c"]["pass"]),
        "3_a8f_four_cells": all(r["pass"] for r in a8f),
        "4p_std_over_sqrt_n_le_0p2": all(x <= RULE_B_ERR for x in err95.values()),
    }
    out = {"seeds": list(seeds), "fits": fits,
           "gates": {"1_fit_ess": g1, "2_cell_ess": g2, "3_ridge_branch_ess": g3,
                     "4_truncation": g4},
           "a8f": a8f, "std_95_bounds": std95, "sampler_error_95": err95,
           "conditions": cond, "holds": all(cond.values())}  # fmt: skip
    if merged:
        h = fits["merged"]["galle_lon_ecl_of_date_deg"]
        out["merged_95"] = [h["lo95"], h["hi95"]]
    return out


# ---------------------------------------------------------------------------- § 6


def range_reading(b: dict, a: dict) -> dict:
    """The table of § 6: whether the 95 % range is final, and what it becomes."""
    r = a["reading"]
    if not b["holds"]:
        failed = [k for k, ok in b["conditions"].items() if not ok]
        return {
            "final": False,
            "range": list(M6C_95),
            "label": "provisional",
            "why": f"rule B fails ({', '.join(failed)})",
            "part_a": (f"P(a < 25 AU) = {r['p']:.3g}" if r["resolved"] else "unresolved"),
        }
    lo, hi = b["merged_95"]
    err = b["sampler_error_95"]
    rng = {"range": [lo, hi], "sampler_errors": [err["lo95"], err["hi95"]]}
    if not r["resolved"]:
        return {
            "final": False,
            **rng,
            "label": "provisional",
            "why": "the low-distance solution, weight 0.001 to about 0.2",
        }
    if r["p"] < 0.05 and r["word"] == "not at all":
        return {
            "final": True,
            **rng,
            "label": None,
            "s1": "S1's low-distance solution holds under 5 % of the posterior",
        }
    m = a["regions"]["M"]["merged"]
    return {
        "final": True,
        "scope": "the headline's prior (e <= 0.3)",
        **rng,
        "label": None,
        "beside": {"lon_deg": m["galle_median"], "a_au": m["a_au"][1], "p": r["p"]},
    }


@ledger.recorded
def summarize(cache: Path = CACHE, synthetic: int | None = None) -> dict:
    rec = ledger.start("s1_mode/summarize", "summarize", n_workers=1,
                       scratch=not _real(cache))  # fmt: skip
    if _real(cache):
        check_inputs()
    a = part_a(cache, synthetic)
    payload = {
        "plan": "docs/notes/S1_mode.md",
        "part_a": a,
        "settings": {
            "nlive": NLIVE,
            "nsteps": NSTEPS,
            "x_nlive": X_NLIVE,
            "x_nsteps": X_NSTEPS,
            "dlogz": DLOGZ,
            "frac_remain": FRAC_REMAIN,
            "seeds": list(SEEDS),
            "x_seeds": list(X_SEEDS),
            "n_draws": N_DRAWS,
            "ultranest": _version("ultranest"),
        },
    }
    if synthetic is None:
        b = part_b(cache)
        payload.update(part_b=b, range_95=range_reading(b, a))
    out = OUT_JSON if _real(cache) else cache / "S1_mode.json"
    out.write_text(json.dumps(payload, indent=1, default=_json) + "\n")
    _log(f"summary written: {out}", cache)
    rec.done(n_likelihood_calls=0, n_integrations=0, sim_years=None)
    return payload


def _json(x):
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    raise TypeError(type(x))


def _version(pkg: str) -> str:
    from importlib.metadata import version

    return version(pkg)


def plan() -> None:
    w = prior_masses()
    check_masses(w)
    for n, m in w.items():
        print(f"{n:3s} {REGIONS[n]}  prior mass {m:.6f}")
    jobs = all_jobs(CACHE, [], None, [])
    print(
        f"{len(jobs)} runs: {sum(j[0].startswith('a/') for j in jobs)} of Part A, "
        f"{sum(j[0].startswith('b/') for j in jobs)} of Part B"
    )
    for name, path, _ in jobs:
        print(f"  {name:14s} {'done' if path.exists() else 'pending'}")


def pipeline(n_workers: int, cache: Path = CACHE, maxiter: int | None = None,
             synthetic: int | None = None, max_draws: int | None = None) -> None:  # fmt: skip
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()  # fmt: skip
    dirty = subprocess.run(["git", "status", "--porcelain", "src"], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()  # fmt: skip
    _log(f"pipeline: commit {head}{' (src modified!)' if dirty else ''}", cache)
    if _real(cache):
        check_inputs()
        check_masses(prior_masses())
        g0r = _g0_result(cache)
        if not (g0r and g0r["pass"]):
            raise SystemExit("G0 has not passed: run `s1_mode g0` first")
    extra_a = ([] if maxiter is None else ["--maxiter", str(maxiter)]) + (
        [] if synthetic is None else ["--synthetic", str(synthetic)]
    )
    extra_b = [] if maxiter is None else ["--maxiter", str(maxiter)]
    bdir = None if _real(cache) else cache / "b"
    jobs = all_jobs(cache, extra_a, bdir, extra_b)
    if synthetic is not None:
        jobs = [j for j in jobs if j[0].startswith("a/")]
    if not schedule(jobs, n_workers, cache):
        _log("pipeline: runs incomplete, stopping before reweighting", cache)
        return
    reweight(n_workers, cache, synthetic, max_draws)
    summarize(cache, synthetic)
    _log("pipeline: done", cache)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("g0", "plan", "pipeline", "job-a", "reweight", "summarize"):
        p = sub.add_parser(name)
        p.add_argument("--cache", type=Path, default=CACHE)
        if name in ("g0", "pipeline", "reweight"):
            p.add_argument("--n-workers", type=int, default=8 if name == "g0" else 12)
        if name in ("pipeline", "job-a"):
            p.add_argument("--maxiter", type=int, default=None, help="smoke runs only")
        if name in ("pipeline", "job-a", "reweight", "summarize"):
            p.add_argument("--synthetic", type=int, default=None, help="smoke runs only")
        if name in ("pipeline", "reweight"):
            p.add_argument("--max-draws", type=int, default=None, help="smoke runs only")
        if name == "job-a":
            p.add_argument("region", choices=tuple(REGIONS))
            p.add_argument("seed", type=int, choices=(*SEEDS, *X_SEEDS))
    a = ap.parse_args(argv)
    cache = a.cache.resolve()
    if getattr(a, "max_draws", None) and _real(cache):
        raise SystemExit("--max-draws is for smoke runs in a scratch --cache")
    if getattr(a, "synthetic", None) is not None and _real(cache):
        raise SystemExit("--synthetic is for smoke runs in a scratch --cache")
    if a.cmd == "g0":
        g0(a.n_workers, cache)
    elif a.cmd == "plan":
        plan()
    elif a.cmd == "job-a":
        do_job_a(a.region, a.seed, cache, a.maxiter, a.synthetic)
    elif a.cmd == "reweight":
        reweight(a.n_workers, cache, a.synthetic, a.max_draws)
    elif a.cmd == "summarize":
        summarize(cache, a.synthetic)
    else:
        pipeline(a.n_workers, cache, a.maxiter, a.synthetic, a.max_draws)


if __name__ == "__main__":
    main()
