"""The data-choice fits and the shape model (pre-registered in docs/notes/M6c_variants.md
and docs/notes/M6d.md; nothing here may differ from those plans).

    uv run python -m fork.export.v14_fits plan
    uv run python -m fork.export.v14_fits g0 FIT [--n-workers 12]        # synthetic gate
    uv run python -m fork.export.v14_fits pipeline FIT [--n-workers 12]  # refuses without G0
    uv run python -m fork.export.v14_fits all [--n-workers 12]           # g0 + pipeline, each fit
    uv run python -m fork.export.v14_fits job FIT K S | reweight FIT | summarize FIT | report

Fits: ``bouvard`` and ``adams`` (free distance on a row selection of Table 77), ``shape38``
and ``shape36`` (all rows, the distance held at 38.4 and 36.154 AU). Each is made as the
headline was confirmed: the fast likelihood about the headline's reference orbit, four
seeds per cell at nlive 4000 (8000 in the low-distance cell), 28 slices, every draw
reweighted to mode (b).
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

from fork import ledger
from fork.export import m6c, m8_runs, v13

ROOT = m6c.ROOT
REAL_CACHE = ROOT / "data" / "cache" / "nested" / "v14_fits"
# Smoke runs only: a scratch cache and an iteration cap, both from the environment. The cap
# is refused in the real cache.
CACHE = Path(os.environ.get("FORK_V14_CACHE", REAL_CACHE))
_CAP = os.environ.get("FORK_V14_MAXITER")
SMOKE_MAXITER = None if _CAP is None else int(_CAP)
if SMOKE_MAXITER is not None and CACHE.resolve() == REAL_CACHE.resolve():
    raise SystemExit("FORK_V14_MAXITER needs a scratch FORK_V14_CACHE")
SEEDS = (1, 2, 3, 4)
MERGED = 0
SLICES = 28
N_DRAWS = 6000
CHUNK = 250
ESS_GATE_CELL = 0.3
ESS_GATE_FIT = 0.5
TRUNC_GATE = 1e-6
G4_DRAWS = 200
G4_TOL = {"chi2": 0.1, "m2lnl": 0.2}
G4_FRACTION = 0.95
G4_FAR = 5.0
G0 = {"seed": 1, "nlive_factor": 0.25, "seeds": (1, 2), "g4_draws": 50}
TOL = {"width_rel": 0.10, "ln_z": 0.5, "median": 0.35, "std_95": 0.4, "p_shape_range": 0.1}
MATERIAL_SHARE = 0.01
E_MAX = 0.3
JPL = v13.GALLE_LON_DEG
HEADLINE_FILE = ROOT / "site" / "data" / "m8b.json"
FROZEN = {
    "site/data/m8b.json": "355104b24966ce216f53ff603931a060a4a5e738a996a683db5fac5c5e962027",
    "site/data/m6c_replacement.json": (
        "d4d1b5c862169984b6152c38b57809da36404e8f84043c5fdf658c05b3f42e32"
    ),
}
SHIFT_AT = 0.02
RATIO_AT = 0.005
LN_B_BANDS = (
    (1.0, "the data do not prefer either model"),
    (3.0, "a weak preference"),
    (5.0, "a strong preference"),
    (math.inf, "a decisive preference"),
)


@dataclass(frozen=True)
class Fit:
    name: str
    tag: int  # rng namespace: 15 the data choices, 16 the shape model
    index: int
    choice: str
    held_a_au: float | None
    out: str  # results file under docs/notes/

    @property
    def n_cells(self) -> int:
        return 2 if self.held_a_au is not None else 4

    @property
    def n_phi(self) -> int:
        return 6 if self.held_a_au is not None else 7

    def nlive(self, k: int) -> int:
        return 4000 if self.held_a_au is not None else 4 * m6c.NLIVE_CELL.get(k, 1000)

    def rng(self, k: int, s: int, *more: int) -> list[int]:
        return [2026, self.tag, 1846, self.index, k + 1, s, *more]


FITS = {
    f.name: f
    for f in (
        Fit("bouvard", 15, 1, "bouvard", None, "M6c_variants.json"),
        Fit("adams", 15, 2, "adams", None, "M6c_variants.json"),
        Fit("shape38", 16, 1, "leverrier", 38.4, "M6d.json"),
        Fit("shape36", 16, 2, "leverrier", v13.LE_VERRIER_A_AU, "M6d.json"),
    )
}
ORDER = ("bouvard", "adams", "shape38", "shape36")


def _log(msg: str, cache: Path) -> None:
    cache.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with (cache / "driver.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def fit_cache(fit: Fit, synthetic: bool = False) -> Path:
    return CACHE / (fit.name + ("__g0" if synthetic else ""))


# ---------------------------------------------------------------------------- the problem


def real_observations(choice: str):
    """The chosen rows' measured values with the headline's frozen error model. A choice may
    leave a series without rows (``bouvard`` drops the ancient one)."""
    from fork.model import observation as om

    rows = v13.select_rows(om.read_rows(), choice)
    em = om.ErrorModel.load(ROOT / m6c.FROZEN["error_model"][0])
    st = om.structure_from_rows(rows, em, allow_empty_groups=True)
    jd = np.array([float(r["epoch_tt_jd"]) for r in rows])
    if np.any(np.diff(jd) < 0) or any(r["unit"] != "deg" for r in rows):
        raise ValueError("rows must be in time order and in degrees")
    return st.with_values(np.array([math.radians(float(r["value"])) for r in rows]))


_P: dict = {}


def problem(fit: Fit, synthetic: bool = False):
    """(obs, ctx, base space, LinearizedProblem, PerturbativeProblem). The base space is the
    uncut one (cell 0 of the shape model, whose φ → θ map does not depend on the cell)."""
    key = (fit.name, synthetic)
    if key in _P:
        return _P[key]
    from fork.inference.fit import uranus_theta_to_x
    from fork.inference.fit_radec import uranus_box_radec
    from fork.inference.priors import free_a_space
    from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch
    from fork.model.linearized import LinearizedProblem
    from fork.model.perturbative import PerturbativeProblem

    m8_runs.check_frozen("headline")
    for rel, digest in FROZEN.items():
        if m6c.sha256(ROOT / rel) != digest:
            raise RuntimeError(f"{rel} differs from the pre-registered file")
    if synthetic:
        obs, ctx, pq, box = v13.synthetic_dataset(G0["seed"], fit.choice)
    else:
        obs = real_observations(fit.choice)
        ctx = ForwardContext.from_horizons(v13.A1_EPOCH_JD_TDB)
        inc_raan = uranus_inc_raan_at_epoch(v13.A1_EPOCH_JD_TDB)
        box = uranus_box_radec(obs, ctx, inc_raan)
        pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, *inc_raan, 0.0, 0.0]))[4:])
    if fit.held_a_au is None:
        base = free_a_space(pq, E_MAX, box, inclined=True)
    else:
        base = v13.shape_space(pq, box, 0, fit.held_a_au)
    lp = LinearizedProblem.build(obs, ctx, base)
    pr = PerturbativeProblem.build(
        obs,
        ctx,
        base,
        m8_runs.u_ref("headline"),
        step_days=m8_runs.STEP_DAYS,
        quadrature=m8_runs.QUADRATURE,
    )
    _P[key] = (obs, ctx, base, lp, pr)
    return _P[key]


def cell_space(fit: Fit, k: int, synthetic: bool = False):
    from fork.inference import nested

    base = problem(fit, synthetic)[2]
    if fit.held_a_au is None:
        return nested.cell_space(base, nested.cells_for(E_MAX)[k])
    return v13.shape_space(base.pq, base.uranus_box, k, fit.held_a_au)


def cell_label(fit: Fit, k: int) -> dict:
    from fork.inference import nested

    if fit.held_a_au is None:
        a, e = nested.cells_for(E_MAX)[k]
        return {"a_au": list(a), "e": None if e is None else list(e)}
    return {"a_au": [fit.held_a_au, fit.held_a_au], "e": list(v13.SHAPE_E_CELLS[k])}


def log_prior_masses(fit: Fit, synthetic: bool = False) -> list[float]:
    return [math.log(cell_space(fit, k, synthetic).slice_weight) for k in range(fit.n_cells)]


# ---------------------------------------------------------------------------- runs


def job_path(k: int, s: int, suffix: str, cache: Path) -> Path:
    return cache / "jobs" / f"c{k}__s{s}{suffix}"


def pid(k: int, s: int) -> str:
    return f"c{k}__" + ("merged" if s == MERGED else f"s{s}")


_S: dict = {}


def _fast_logl(phi: np.ndarray) -> float:
    return _S["pr"].log_likelihood(phi, second_order=m8_runs.SECOND_ORDER)


def _phi_transform(u: np.ndarray) -> np.ndarray:
    return _S["space"].prior_transform(np.concatenate([np.full(4, 0.5), u]))[4:]


@ledger.recorded
def do_job(
    fit: Fit,
    k: int,
    s: int,
    synthetic: bool = False,
    maxiter: int | None = None,
    cache: Path | None = None,
) -> None:
    import dynesty

    from fork import ledger
    from fork.inference.nested import BOUND, DLOGZ, SAMPLE, NestedResult, _ProgressLog

    cache = fit_cache(fit, synthetic) if cache is None else cache
    maxiter = SMOKE_MAXITER if maxiter is None else maxiter
    if maxiter is not None and CACHE.resolve() == REAL_CACHE.resolve():
        raise SystemExit("capped runs go to a scratch cache")
    out = job_path(k, s, ".pkl", cache)
    if out.exists():
        print(f"{fit.name}/c{k}/s{s}: cached")
        return
    pr = problem(fit, synthetic)[4]
    cs = cell_space(fit, k, synthetic)
    _S.update(pr=pr, space=cs)
    n_phi = len(cs.unit_names) - 4
    periodic = [i - 4 for i in cs.periodic_unit_dims]
    if min(periodic) < 0 or n_phi != fit.n_phi:
        raise RuntimeError("unexpected prior layout")
    nlive = max(100, round(fit.nlive(k) * (G0["nlive_factor"] if synthetic else 1)))
    ckpt = job_path(k, s, ".ckpt", cache)
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    label = f"v14/{fit.name}{'/g0' if synthetic else ''}/c{k}/s{s}"
    log = _ProgressLog(cache / "progress.log", m8_runs.PROGRESS_EVERY_S, label)
    run_kw = dict(
        dlogz=DLOGZ,
        maxiter=maxiter,
        print_progress=True,
        print_func=log,
        checkpoint_file=str(ckpt),
        checkpoint_every=m8_runs.CHECKPOINT_EVERY_S,
    )
    resumed = ckpt.exists()
    t0 = time.time()
    real = CACHE.resolve() == REAL_CACHE.resolve()
    rec = ledger.start(f"v14/{fit.name}{'/g0' if synthetic else ''}/c{k}/s{s}", "nested",
                       n_workers=1, synthetic=synthetic, capped=maxiter is not None,
                       scratch=not real, n_integrations=0, sim_years=0.0)  # fmt: skip
    n0 = 0
    if resumed:
        sampler = dynesty.NestedSampler.restore(str(ckpt))
        n0 = int(sampler.ncall)
        log.line(f"resumed from {ckpt.name} at iteration {sampler.it}")
        rec.set(n_likelihood_calls=lambda: int(sampler.ncall) - n0)
        sampler.run_nested(resume=True, **run_kw)
    else:
        log.line(
            f"start: ndim {n_phi}, nlive {nlive}, {BOUND}/{SAMPLE} {SLICES} slices, "
            f"dlogz {DLOGZ}, rng {fit.rng(k, s)}"
        )
        sampler = dynesty.NestedSampler(
            _fast_logl,
            _phi_transform,
            n_phi,
            nlive=nlive,
            bound=BOUND,
            sample=SAMPLE,
            periodic=periodic,
            slices=SLICES,
            rstate=np.random.default_rng(fit.rng(k, s)),
        )
        rec.set(n_likelihood_calls=lambda: int(sampler.ncall))
        sampler.run_nested(**run_kw)
    res = NestedResult(cs, sampler.results, nlive, DLOGZ, time.time() - t0)
    res.save(out)
    ckpt.unlink(missing_ok=True)
    log.line(f"done: ln Z {res.logz:.3f} ± {res.logzerr:.3f}, {res.niter} it, {res.ncall} calls")
    # this segment's calls: a resumed run's earlier segments wrote their own rows
    rec.done(
        n_likelihood_calls=res.ncall - n0,
        ncall_total=res.ncall,
        estimated=False,
        resumed=resumed,
        nlive=nlive,
        niter=res.niter,
        slices=SLICES,
        likelihood="perturbative, second order off",
    )


def run_jobs(fit: Fit, n_slots: int, synthetic: bool = False) -> bool:
    """Every run of the fit as its own process, ``n_slots`` at a time, the long cell first;
    one retry per failed run (it resumes from its checkpoint)."""
    cache = fit_cache(fit, synthetic)
    seeds = G0["seeds"] if synthetic else SEEDS
    jobs = [(k, s) for k in range(fit.n_cells) for s in seeds]
    pending = [j for j in jobs if not job_path(*j, ".pkl", cache).exists()]
    _log(f"{fit.name}: {len(pending)} of {len(jobs)} runs pending, {n_slots} slots", cache)
    running: dict = {}
    failures: dict = {}
    while pending or running:
        for j, (proc, t0) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[j]
            if rc == 0 and job_path(*j, ".pkl", cache).exists():
                _log(f"done c{j[0]}/s{j[1]} ({(time.time() - t0) / 60:.1f} min)", cache)
            else:
                failures[j] = failures.get(j, 0) + 1
                _log(f"FAILED c{j[0]}/s{j[1]} (exit {rc}, attempt {failures[j]})", cache)
                if failures[j] < 2:
                    pending.insert(0, j)
        while pending and len(running) < n_slots:
            j = pending.pop(0)
            logf = job_path(*j, ".log", cache)
            logf.parent.mkdir(parents=True, exist_ok=True)
            cmd = [
                sys.executable,
                "-m",
                "fork.export.v14_fits",
                "job",
                fit.name,
                str(j[0]),
                str(j[1]),
                *(["--synthetic"] if synthetic else []),
            ]
            with logf.open("a") as fh:
                proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
            running[j] = (proc, time.time())
            _log(f"start c{j[0]}/s{j[1]} (pid {proc.pid})", cache)
        time.sleep(10)
    ok = all(job_path(*j, ".pkl", cache).exists() for j in jobs)
    _log(f"{fit.name}: runs finished; failures: {failures or 'none'}; all done: {ok}", cache)
    return ok


# ---------------------------------------------------------------------------- reweighting


def load_posterior(k: int, s: int, cache: Path, seeds=SEEDS):
    from dynesty.utils import merge_runs

    from fork.inference.nested import NestedResult

    if s != MERGED:
        return NestedResult.load(job_path(k, s, ".pkl", cache))
    runs = [NestedResult.load(job_path(k, t, ".pkl", cache)) for t in seeds]
    merged = merge_runs([r.results for r in runs], print_progress=False)
    return NestedResult(
        runs[0].space,
        merged,
        sum(r.nlive for r in runs),
        runs[0].dlogz,
        sum(r.wall_s for r in runs),
    )


def draws(fit: Fit, k: int, s: int, cache: Path, seeds=SEEDS) -> dict:
    path = cache / "reweight" / pid(k, s) / "draws.npz"
    if path.exists():
        d = np.load(path)
        return {key: d[key] for key in d.files}
    res = load_posterior(k, s, cache, seeds)
    rng = np.random.default_rng(fit.rng(k, s, 7))
    eq = res.equal_weight(rng)
    phi = eq[rng.choice(len(eq), min(N_DRAWS, len(eq)), replace=False)]
    d = {
        "phi": phi,
        "ln_z_fast": np.array(res.logz),
        "ln_z_fast_err": np.array(res.logzerr),
        "ncall": np.array(res.ncall),
        "niter": np.array(res.niter),
        "stopped_on_dlogz": np.array(_stopped(res)),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("draws.tmp.npz")
    np.savez(tmp, **d)
    tmp.replace(path)
    return d


def _stopped(res) -> bool:
    """G1: these runs have no cap on iterations or calls, so a saved run ended on dlogz. A
    capped smoke run is marked as not stopped."""
    return bool(res.niter > 0 and SMOKE_MAXITER is None)


_RW: dict = {}


def _init_worker(name: str, synthetic: bool) -> None:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from fork.inference.profile import GalleGeometry

    fit = FITS[name]
    _, ctx, space, lp, pr = problem(fit, synthetic)
    _RW.update(fit=fit, ctx=ctx, space=space, lp=lp, pr=pr, geo=GalleGeometry.from_horizons())


def _reweight_chunk(task: tuple) -> tuple[str, int]:
    from fork.inference.profile import galle_lon_of_date_rad

    k, s, c, lo, hi, cache = task
    out = Path(cache) / "reweight" / pid(k, s) / f"chunk_{c:03d}.npz"
    if out.exists():
        return pid(k, s), c
    fit, ctx, space, lp, pr, geo = (_RW[x] for x in ("fit", "ctx", "space", "lp", "pr", "geo"))
    phi = np.load(out.parent / "draws.npz")["phi"][lo:hi]
    n = len(phi)
    rng = np.random.default_rng(fit.rng(k, s, 9, c))
    a = {
        "lnl_fast": np.full(n, -np.inf),
        "lnl_b": np.full(n, -np.inf),
        "b_trunc_m2lnl": np.full(n, np.nan),
        "galle_deg": np.full(n, np.nan),
        "pert": np.full((n, 3), np.nan),
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
        th = space.theta(np.concatenate([rng.multivariate_normal(m.u_star, m.cov_u), f]))
        a["galle_deg"][i] = math.degrees(galle_lon_of_date_rad(th, ctx, geo)) % 360.0
        a["pert"][i] = (10.0 ** th[6], th[7], th[8])
    tmp = out.with_name(out.stem + ".tmp.npz")
    np.savez(tmp, **a)
    tmp.replace(out)
    return pid(k, s), c


def _pool(fit: Fit, synthetic: bool, n_workers: int):
    from concurrent.futures import ProcessPoolExecutor

    return ProcessPoolExecutor(n_workers, initializer=_init_worker, initargs=(fit.name, synthetic))


@ledger.recorded
def reweight(fit: Fit, n_workers: int, synthetic: bool = False) -> None:
    from concurrent.futures import as_completed

    from fork import ledger

    cache = fit_cache(fit, synthetic)
    seeds = G0["seeds"] if synthetic else SEEDS
    rec = ledger.start(f"v14/{fit.name}{'/g0' if synthetic else ''}/reweight", "reweight",
                       n_workers=n_workers, synthetic=synthetic,
                       scratch=CACHE.resolve() != REAL_CACHE.resolve())  # fmt: skip
    tasks = []
    for k in range(fit.n_cells):
        for s in (*seeds, MERGED):
            n = len(draws(fit, k, s, cache, seeds)["phi"])
            for c, lo in enumerate(range(0, n, CHUNK)):
                if not (cache / "reweight" / pid(k, s) / f"chunk_{c:03d}.npz").exists():
                    tasks.append((k, s, c, lo, min(lo + CHUNK, n), str(cache)))
    _log(f"{fit.name}: reweight, {len(tasks)} chunks to do, {n_workers} workers", cache)
    with _pool(fit, synthetic, n_workers) as ex:
        futs = [ex.submit(_reweight_chunk, t) for t in tasks]
        for i, f in enumerate(as_completed(futs)):
            f.result()
            if (i + 1) % 48 == 0 or i + 1 == len(futs):
                _log(f"{fit.name}: reweight {i + 1}/{len(futs)} chunks", cache)
    if not tasks:
        rec.discard()
        return
    rec.done(
        n_likelihood_calls=sum(t[4] - t[3] for t in tasks),
        n_integrations=None,
        sim_years=None,
        note="mode (b) + fast + Galle longitude per draw",
    )


def load_reweighted(fit: Fit, k: int, s: int, cache: Path, seeds=SEEDS) -> dict:
    d = dict(draws(fit, k, s, cache, seeds))
    n = len(d["phi"])
    parts = [
        np.load(cache / "reweight" / pid(k, s) / f"chunk_{c:03d}.npz")
        for c in range(math.ceil(n / CHUNK))
    ]
    for key in parts[0].files:
        d[key] = np.concatenate([p[key] for p in parts])
    return d


# ---------------------------------------------------------------------------- G4


def _g4_one(f: np.ndarray) -> list[float]:
    """[χ²(mode b) − χ²(full, minimised over Uranus), −2 ln L(mode b) − Laplace(exact Hessian)]."""
    lp = _RW["lp"]
    m = lp.marginal(f, "b")
    u_min, chi2_full, _ = lp.profile(f, m.u_star)
    return [float(m.chi2_min - chi2_full), float(m.m2lnl - lp.laplace(f, u_min))]


def g4_verdict(diffs: np.ndarray) -> dict:
    d = np.abs(np.asarray(diffs, dtype=float))
    tol = np.array([G4_TOL["chi2"], G4_TOL["m2lnl"]])
    inside = (d <= tol).mean(axis=0)
    far = int(np.any(d > G4_FAR * tol, axis=1).sum()) + int(np.any(~np.isfinite(d), axis=1).sum())
    return {
        "n": len(d),
        "tolerances": G4_TOL,
        "fraction_within": inside.tolist(),
        "max_abs": np.nanmax(d, axis=0).tolist(),
        "n_beyond_5x": far,
        "pass": bool(np.all(inside >= G4_FRACTION) and far == 0),
    }


def g4(fit: Fit, phi: np.ndarray, n_workers: int, synthetic: bool = False) -> dict:
    cache = fit_cache(fit, synthetic)
    path = cache / "g4.npz"
    if path.exists() and len(np.load(path)["diffs"]) == len(phi):
        return g4_verdict(np.load(path)["diffs"])
    _log(f"{fit.name}: G4 on {len(phi)} draws", cache)
    with _pool(fit, synthetic, n_workers) as ex:
        diffs = np.array(list(ex.map(_g4_one, phi, chunksize=2)))
    np.savez(path, diffs=diffs, phi=phi)
    return g4_verdict(diffs)


# ---------------------------------------------------------------------------- summary


def _q(x: np.ndarray) -> dict:
    lo95, lo, med, hi, hi95 = np.quantile(x, [0.025, 0.16, 0.5, 0.84, 0.975]).tolist()
    return {"median": med, "lo68": lo, "hi68": hi, "lo95": lo95, "hi95": hi95}


def _circ_range(values, circular: bool = True) -> float:
    v = np.asarray(values, dtype=float)
    if circular:
        v = (v - v[0] + 180.0) % 360.0 - 180.0
    return float(v.max() - v.min())


def assemble(fit: Fit, synthetic: bool = False) -> dict:
    """Cells, fits, gates G1 to G3 and T1. The longitudes of the merged fit and of each
    single-seed fit are kept under ``_galle`` (dropped before writing)."""
    from fork.inference.reweight import N_PER_CELL, reweight_cell

    cache = fit_cache(fit, synthetic)
    seeds = G0["seeds"] if synthetic else SEEDS
    log_c = log_prior_masses(fit, synthetic)
    post, cells = {}, []
    for k in range(fit.n_cells):
        for s in (*seeds, MERGED):
            d = load_reweighted(fit, k, s, cache, seeds)
            cr = reweight_cell(
                d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]), float(d["ln_z_fast_err"])
            )
            d["cr"] = cr
            post[(k, s)] = d
            _, ix = m8_runs._resampled([cr], [0.0], np.random.default_rng(fit.rng(k, s, 13)), cr.n)
            g = m6c.lon_summary(d["galle_deg"][ix], {})
            ok = np.isfinite(d["lnl_b"])
            cells.append(
                {
                    "cell": k,
                    **cell_label(fit, k),
                    "seed": s,
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
                    "b_trunc_m2lnl_max": float(np.nanmax(d["b_trunc_m2lnl"])) if ok.any() else None,
                    "n_b_failed": int(np.sum(~ok)),
                }
            )
    fits, galle_of, gate_fit = {}, {}, []
    for s in (MERGED, *seeds):
        cs = [post[(k, s)] for k in range(fit.n_cells)]
        res, ix = m8_runs._resampled(
            [c["cr"] for c in cs], log_c, np.random.default_rng(fit.rng(-1, s, 11)), N_PER_CELL
        )
        galle = np.concatenate([c["galle_deg"] for c in cs])[ix]
        pert = np.concatenate([c["pert"] for c in cs])[ix]
        shares = np.exp(np.array(log_c) + np.array([c["cr"].ln_z for c in cs]) - res.ln_z)
        name = "merged" if s == MERGED else f"s{s}"
        galle_of[name] = galle
        fits[name] = {
            "ln_z": res.ln_z,
            "ln_z_err": res.ln_z_err,
            "shares": shares.tolist(),
            "galle_lon_ecl_of_date_deg": m6c.lon_summary(galle, {"jpl": JPL}),
            "p_within_0p87_deg": v13.peirce_probability(galle),
            "mass_earth": _q(pert[:, 0]),
            "a_au": _q(pert[:, 1]),
            "e": _q(pert[:, 2]),
        }
        _, _, share = m8_runs.fit_ess(res)
        gate_fit.append({"fit": name, "reweighting_share": share, "pass": share >= ESS_GATE_FIT})
    shares = fits["merged"]["shares"]
    g2 = [c for c in cells if c["ess_frac"] < ESS_GATE_CELL]
    g3 = [
        c for c in cells if c["b_trunc_m2lnl_max"] is None or c["b_trunc_m2lnl_max"] >= TRUNC_GATE
    ]
    t1_cells = []
    for k in range(fit.n_cells):
        four = [c for c in cells if c["cell"] == k and c["seed"] != MERGED]
        w = np.array([c["galle_width68"] for c in four])
        rel = float(np.max(np.abs(w / w.mean() - 1.0)))
        lnz = _circ_range([c["ln_z"] for c in four], circular=False)
        material = shares[k] >= MATERIAL_SHARE
        t1_cells.append(
            {
                "cell": k,
                "share": shares[k],
                "material": bool(material),
                "widths": w.tolist(),
                "max_rel_dev": rel,
                "ln_z_range": lnz,
                "width_pass": bool(rel <= TOL["width_rel"]),
                "ln_z_pass": bool(lnz <= TOL["ln_z"]),
            }
        )
    single = [fits[f"s{s}"] for s in seeds]
    med_range = _circ_range([f["galle_lon_ecl_of_date_deg"]["median"] for f in single])
    fit_lnz_range = _circ_range([f["ln_z"] for f in single], circular=False)
    std95 = {
        b: float(np.std([f["galle_lon_ecl_of_date_deg"][b] for f in single], ddof=1))
        for b in ("lo95", "hi95")
    }
    t1 = {
        "cells": t1_cells,
        "median_range_single_seed_fits": med_range,
        "ln_z_range_single_seed_fits": fit_lnz_range,
        "pass": bool(
            all(c["ln_z_pass"] and (c["width_pass"] or not c["material"]) for c in t1_cells)
            and med_range <= TOL["median"]
            and (fit.held_a_au is None or fit_lnz_range <= TOL["ln_z"])
        ),
    }
    return {
        "fit": fit.name,
        "choice": fit.choice,
        "held_a_au": fit.held_a_au,
        "n_measurements": problem(fit, synthetic)[0].n_measurements,
        "settings": {
            "nlive": [fit.nlive(k) for k in range(fit.n_cells)],
            "slices": SLICES,
            "seeds": list(seeds),
            "n_draws": N_DRAWS,
            "rng": fit.rng(0, 1)[:4],
            "synthetic": synthetic,
        },
        "cells": cells,
        "fits": fits,
        "gates": {
            "G1": {"pass": bool(all(post[key]["stopped_on_dlogz"] for key in post))},
            "G2": {
                "cell_failures": [f"c{c['cell']}/s{c['seed']}" for c in g2],
                "min_cell_ess_frac": min(c["ess_frac"] for c in cells),
                "fits": gate_fit,
                "pass": not g2 and all(g["pass"] for g in gate_fit),
            },
            "G3": {"failures": [f"c{c['cell']}/s{c['seed']}" for c in g3], "pass": not g3},
        },
        "T1": t1,
        "std_95_bounds": std95,
        "label_95": "final" if max(std95.values()) <= TOL["std_95"] else "provisional",
        "_galle": galle_of,
        "_post": post,
    }


def g4_draws(a: dict, fit: Fit, n: int) -> np.ndarray:
    """n draws of φ from the merged fit (equal weight), fixed by the fit's rng."""
    from fork.inference.reweight import N_PER_CELL

    cs = [a["_post"][(k, MERGED)] for k in range(fit.n_cells)]
    _, ix = m8_runs._resampled(
        [c["cr"] for c in cs],
        log_prior_masses(fit, a["settings"]["synthetic"]),
        np.random.default_rng(fit.rng(-1, 0, 19)),
        N_PER_CELL,
    )
    phi = np.concatenate([c["phi"] for c in cs])[ix]
    return phi[np.random.default_rng(fit.rng(-1, 0, 23)).choice(len(phi), n, replace=False)]


def headline() -> dict:
    return json.loads(HEADLINE_FILE.read_text(encoding="utf-8"))["fits"]["headline/merged"]


def read_shift(fit_lon: dict, head_lon: dict) -> dict:
    """§ 5 of M6c_variants.md: the shift of the median and the width against the headline's."""
    d = (fit_lon["median"] - head_lon["median"] + 180.0) % 360.0 - 180.0

    def inside(lo: str, hi: str) -> bool:
        rel = lambda v: (v - head_lon["median"] + 180.0) % 360.0 - 180.0  # noqa: E731
        return rel(head_lon[lo]) <= d <= rel(head_lon[hi])

    if abs(d) < TOL["median"]:
        w = "the choice does not move the prediction"
    elif inside("lo68", "hi68"):
        w = "it moves the prediction within the headline's 68 % interval"
    elif inside("lo95", "hi95"):
        w = "it moves the prediction outside the headline's 68 % interval"
    else:
        w = "it moves the prediction outside the headline's 95 % interval"
    r = fit_lon["width68"] / head_lon["width68"]
    ww = "the same width" if 0.9 <= r <= 1.1 else ("wider" if r > 1.1 else "narrower")
    edges = [
        TOL["median"],
        *(
            abs((head_lon[b] - head_lon["median"] + 180.0) % 360.0 - 180.0)
            for b in ("lo68", "hi68", "lo95", "hi95")
        ),
    ]
    return {
        "shift_deg": d,
        "shift_wording": w,
        "shift_at_a_threshold": bool(any(abs(abs(d) - e) <= SHIFT_AT for e in edges)),
        "width_ratio": r,
        "width_wording": ww,
        "width_at_a_threshold": bool(min(abs(r - 0.9), abs(r - 1.1)) <= RATIO_AT),
        "jpl_in_68": fit_lon["vs"]["jpl"]["in_68"],
        "jpl_in_95": fit_lon["vs"]["jpl"]["in_95"],
    }


def band(x: float) -> str:
    return next(name for edge, name in LN_B_BANDS if abs(x) < edge)


def read_ln_b(ln_b: float, err: float) -> dict:
    lo, hi = band(ln_b - 2.0 * err), band(ln_b + 2.0 * err)
    side = "for the shape model" if ln_b > 0 else "for free distance"
    if lo == hi == band(ln_b):
        w = lo if lo == LN_B_BANDS[0][1] else f"{lo} {side}"
    else:
        w = f"between '{lo}' and '{hi}' ({side} where there is one)"
    return {"ln_b": ln_b, "err": err, "wording": w}


def read_ratio(r: dict, p_single: list[float]) -> dict:
    rng_p = max(p_single) - min(p_single)
    if rng_p > TOL["p_shape_range"]:
        w = "not repeatable; not read"
    elif r["ratio"] is None or r["ratio_lo68"] is None:
        w = "not defined (P_free is 0)"
    elif r["ratio_lo68"] <= 1.0 <= r["ratio_hi68"]:
        w = (
            "holding the distance at his value did not change the probability of a "
            "prediction as close as his"
        )
    elif r["ratio_lo68"] > 1.0:
        w = f"made it more probable, by a factor {r['ratio']:.2f}"
    else:
        w = (
            f"made it less probable, by a factor {1.0 / r['ratio']:.2f}"
            if r["ratio"] > 0
            else "made it less probable (P_shape is 0)"
        )
    return {**r, "p_shape_single_seed": p_single, "p_shape_range": rng_p, "wording": w}


def free_distance_galle() -> tuple[np.ndarray, dict]:
    """The published higher-resolution headline fit's longitudes, rebuilt from its cached
    draws exactly as ``m8_runs.summarize_followup`` built them; the median must reproduce."""
    from fork.inference.reweight import N_PER_CELL, reweight_cell

    camp = m8_runs.CAMPAIGNS["m8b"]
    crs, galle = [], []
    for k in range(4):
        d = m8_runs.load_reweighted("headline", k, MERGED, camp.cache)
        crs.append(
            reweight_cell(
                d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]), float(d["ln_z_fast_err"])
            )
        )
        galle.append(d["galle_deg"])
    log_c = [math.log(m8_runs._cell_space("headline", k).slice_weight) for k in range(4)]
    rng = np.random.default_rng([2026, 8, 1846, 0, 0, MERGED, m8_runs.FIT_RNG_TAG, camp.rng_tag])
    res, ix = m8_runs._resampled(crs, log_c, rng, N_PER_CELL)
    lon = np.concatenate(galle)[ix]
    pub = headline()
    med = m6c.lon_summary(lon, {})["median"]
    if (
        abs(med - pub["galle_lon_ecl_of_date_deg"]["median"]) > 0.001
        or abs(res.ln_z - pub["ln_z"]) > 1e-6
    ):
        raise RuntimeError(f"the published fit is not reproduced: median {med}, ln Z {res.ln_z}")
    return lon, pub


def _public(a: dict) -> dict:
    return {k: v for k, v in a.items() if not k.startswith("_")}


@ledger.recorded
def summarize(fit: Fit, n_workers: int, synthetic: bool = False) -> dict:
    cache = fit_cache(fit, synthetic)
    rec = ledger.start(f"v14/{fit.name}{'/g0' if synthetic else ''}/summarize", "summarize",
                       n_workers=n_workers, synthetic=synthetic,
                       scratch=CACHE.resolve() != REAL_CACHE.resolve(),
                       note="assembly and G4 (mode (b) against the full likelihood)")  # fmt: skip
    a = assemble(fit, synthetic)
    n4 = G0["g4_draws"] if synthetic else G4_DRAWS
    a["gates"]["G4"] = g4(fit, g4_draws(a, fit, n4), n_workers, synthetic)
    gates_ok = all(g["pass"] for g in a["gates"].values())
    if synthetic:
        cells_ok = not a["gates"]["G2"]["cell_failures"] and a["gates"]["G3"]["pass"]
        out = _public(a) | {
            "g0": {
                "pass": bool(cells_ok and a["gates"]["G4"]["pass"]),
                "rule": "every cell posterior ESS/N ≥ 0.3, truncation below 1e-6, G4 on 50 draws",
            }
        }
        (cache / "g0.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
        _log(f"{fit.name}: G0 {'passes' if out['g0']['pass'] else 'FAILS'}", cache)
        rec.done("gate passed" if out["g0"]["pass"] else "gate failed")
        return out
    out = _public(a)
    out["read"] = bool(gates_ok and a["T1"]["pass"])
    lon = out["fits"]["merged"]["galle_lon_ecl_of_date_deg"]
    if fit.held_a_au is None:
        out["reading"] = read_shift(lon, headline()["galle_lon_ecl_of_date_deg"])
    else:
        free, pub = free_distance_galle()
        ratio = v13.peirce_ratio(
            a["_galle"]["merged"], free, np.random.default_rng(fit.rng(-1, 0, 17))
        )
        ln_b, err = v13.log_bayes_factor(
            out["fits"]["merged"]["ln_z"],
            out["fits"]["merged"]["ln_z_err"],
            pub["ln_z"],
            pub["ln_z_err"],
        )
        repl = json.loads((ROOT / "site/data/m6c_replacement.json").read_text())
        out["reading"] = {
            "ratio": read_ratio(ratio, [out["fits"][f"s{s}"]["p_within_0p87_deg"] for s in SEEDS]),
            "bayes_factor": read_ln_b(ln_b, err),
            "free_distance": {
                "ln_z": pub["ln_z"],
                "ln_z_err": pub["ln_z_err"],
                "p_within_0p87_deg": ratio["p_free"],
                "reported_posterior_p_printed_beside": repl["peirce_p_within_0p87_deg_of_jpl"],
            },
        }
    if not out["read"]:
        out["reading"] = {
            "not_read": "a gate or T1 failed; numbers reported, no band read",
            "values": out["reading"],
        }
    (cache / "summary.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    _log(f"{fit.name}: summarized; gates {gates_ok}, T1 {a['T1']['pass']}", cache)
    return out


def _not_read_reasons(s: dict) -> tuple[str, list[str]]:
    """Why a fit is not read (§ 5): the gates that fail ("not validated") and the T1 tests that
    fail ("not repeatable"); a gate failure sets the label."""
    gates = [g for g in ("G1", "G2", "G3", "G4") if not s["gates"].get(g, {}).get("pass", False)]
    t1, why = s["T1"], [f"{g} fails" for g in gates]
    for c in t1["cells"]:
        if c["material"] and not c["width_pass"]:
            why.append(f"T1 fails: cell {c['cell']} 68 % widths {100 * c['max_rel_dev']:.1f} %"
                       f" from their mean vs {100 * TOL['width_rel']:g} %")  # fmt: skip
        if not c["ln_z_pass"]:
            why.append(f"T1 fails: cell {c['cell']} ln Z range {c['ln_z_range']:.2f}"
                       f" vs {TOL['ln_z']:g}")  # fmt: skip
    if t1["median_range_single_seed_fits"] > TOL["median"]:
        why.append(f"T1 fails: median range {t1['median_range_single_seed_fits']:.3f}°"
                   f" vs {TOL['median']:g}°")  # fmt: skip
    if s.get("held_a_au") is not None and t1["ln_z_range_single_seed_fits"] > TOL["ln_z"]:
        why.append(f"T1 fails: ln Z range of the single-seed fits"
                   f" {t1['ln_z_range_single_seed_fits']:.2f} vs {TOL['ln_z']:g}")  # fmt: skip
    return ("not validated" if gates else "not repeatable"), why


def unread_public(s: dict) -> dict:
    """A fit that is not read, as the results file carries it: its status and label (§ 5), the
    gates, T1 and the posterior summaries, and nothing that compares it with the headline or
    with Neptune's place, no band wording and no final/provisional label."""
    label, why = _not_read_reasons(s)
    out = {k: v for k, v in s.items() if k not in ("reading", "label_95")}
    out["fits"] = {
        name: {k: v for k, v in f.items() if k != "p_within_0p87_deg"}
        | {"galle_lon_ecl_of_date_deg":
           {k: v for k, v in f["galle_lon_ecl_of_date_deg"].items() if k != "vs"}}
        for name, f in s["fits"].items()
    }  # fmt: skip
    out["status"] = f"not read ({'; '.join(why) or 'a gate or T1 failed'})"
    out["label"] = label
    return out


def report() -> None:
    """Collect the finished fits' summaries into the two results files."""
    for out_name in sorted({f.out for f in FITS.values()}):
        fits = {}
        for name in ORDER:
            f = FITS[name]
            g0 = fit_cache(f, True) / "g0.json"
            s = fit_cache(f) / "summary.json"
            if f.out != out_name:
                continue
            fits[name] = {
                "g0": json.loads(g0.read_text())["g0"] if g0.exists() else None,
                "summary": json.loads(s.read_text()) if s.exists() else None,
            }
            sm = fits[name]["summary"]
            if sm is not None and not sm["read"]:
                fits[name]["summary"] = unread_public(sm)
        plan = "M6c_variants.md" if out_name.startswith("M6c") else "M6d.md"
        payload = {"plan": f"docs/notes/{plan}", "headline": headline(), "fits": fits}
        (ROOT / "docs" / "notes" / out_name).write_text(
            json.dumps(payload, indent=1, default=float, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"wrote docs/notes/{out_name}")


def g0(fit: Fit, n_workers: int) -> bool:
    """The synthetic gate: the whole chain at a quarter of the live points, two seeds."""
    path = fit_cache(fit, True) / "g0.json"
    if not path.exists():
        if not run_jobs(fit, n_workers, synthetic=True):
            return False
        reweight(fit, n_workers, synthetic=True)
        summarize(fit, n_workers, synthetic=True)
    return bool(json.loads(path.read_text())["g0"]["pass"])


def pipeline(fit: Fit, n_workers: int) -> None:
    cache = fit_cache(fit)
    path = fit_cache(fit, True) / "g0.json"
    if not path.exists() or not json.loads(path.read_text())["g0"]["pass"]:
        _log(f"{fit.name}: G0 has not passed; the real rows are not run", cache)
        return
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    _log(f"{fit.name}: pipeline, commit {head}, {n_workers} workers", cache)
    if not run_jobs(fit, n_workers):
        _log(f"{fit.name}: runs incomplete; stopped", cache)
        return
    reweight(fit, n_workers)
    summarize(fit, n_workers)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "g0", "pipeline", "all", "job", "reweight", "summarize", "report"):
        p = sub.add_parser(name)
        p.add_argument("--n-workers", type=int, default=12)
        if name in ("g0", "pipeline", "job", "reweight", "summarize"):
            p.add_argument("fit", choices=ORDER)
        if name == "job":
            p.add_argument("k", type=int)
            p.add_argument("s", type=int)
        if name in ("job", "reweight", "summarize"):
            p.add_argument("--synthetic", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "plan":
        for name in ORDER:
            f = FITS[name]
            done = sum(
                job_path(k, s, ".pkl", fit_cache(f)).exists()
                for k in range(f.n_cells)
                for s in SEEDS
            )
            print(
                f"{name}: G0 {'done' if (fit_cache(f, True) / 'g0.json').exists() else 'no'}"
                f", {done}/{f.n_cells * 4} runs, summary "
                f"{'yes' if (fit_cache(f) / 'summary.json').exists() else 'no'}"
            )
    elif a.cmd == "g0":
        g0(FITS[a.fit], a.n_workers)
    elif a.cmd == "pipeline":
        pipeline(FITS[a.fit], a.n_workers)
    elif a.cmd == "all":
        for name in ORDER:
            if g0(FITS[name], a.n_workers):
                pipeline(FITS[name], a.n_workers)
            else:
                _log(f"{name}: stopped at G0", fit_cache(FITS[name]))
        report()
    elif a.cmd == "job":
        do_job(FITS[a.fit], a.k, a.s, a.synthetic)
    elif a.cmd == "reweight":
        reweight(FITS[a.fit], a.n_workers, a.synthetic)
    elif a.cmd == "summarize":
        summarize(FITS[a.fit], a.n_workers, a.synthetic)
    else:
        report()


if __name__ == "__main__":
    main()
