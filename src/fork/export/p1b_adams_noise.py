"""Phase 1b step 3c: noise against systematic (docs/notes/P1B_adams_noise.md, pre-registered).

Is the shift between the clean-data solution of step 3b's D3(i) and Adams's real-data solution
within what independent observational noise produces in his layout through his method? Starting
from the noise-free his-layout normal places, add Gaussian noise per place under two bracketing
models ((a) his § 52 residual scatter scaled by √(N/(N − p)); (b) set D's step-3 scatter over √N
per 3-year window) and a Student-t sensitivity (c); solve each realisation by the warm-start
least squares (primary) and by his own elimination chain (secondary); read P(Δ ≥ Δ_obs) by the
plan's § 4 rows. ``prerun`` writes the frozen inputs, ``draws`` the realisations (a process
pool, minutes), ``summarize`` the note's JSON and figure. Caches under
data/cache/p1b_adams_noise/ (gitignored).
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from fork.export import p1b_adams as pa
from fork.export import p1b_adams_repro as R
from fork.historical import adams_1846 as A

ROOT = R.ROOT
CACHE = ROOT / "data" / "cache" / "p1b_adams_noise"
PRERUN_JSON = CACHE / "prerun.json"
PROGRESS_LOG = CACHE / "progress.log"
OUT_JSON = R.NOTES / "P1B_adams_noise.json"
OUT_PNG = R.NOTES / "P1B_adams_noise.png"

SEED_NAMESPACE = (2026, 11, 3, 8)  # [2026, 11, 3, 8, model index, hyp]
N_DRAWS = 10_000
N_MULTI = 1_000  # the multistart check's realisations per model × hypothesis
MODELS = ("a", "b", "c")
HYPS = (1, 2)
METHODS = ("d1c", "chain")
NU_T = 4  # sensitivity (c)
ESCAPE_DEG = 30.0  # θ farther than this from the clean θ: another minimum (kept, counted)
MULTI_DIFF_DEG = 1.0  # multistart check: a global minimum farther than this in the 1846.762 lon
MULTI_TRIGGER = 0.01  # ... in more than this fraction → every draw rerun with the multistart
FAIL_TRIGGER = 0.01  # failures above this fraction → the upper-bound P beside the first
P_NOISE_SUFFICES = 0.05
P_SYSTEMATIC = 0.003
WINDOW_MODERN_YR = A.UNIT_YR / 2  # ±1.518 yr about a modern place
WINDOW_ANCIENT_YR = 1.0
CHUNK = 250
CHUNK_DIR = CACHE / "chunks"
NFEV_DEFAULT = 2000  # per solve, (a) and (b)
NFEV_CAP_C = 500  # per solve in sensitivity (c); a capped start is skipped, capped draws counted
WATCHDOG_S = 60.0  # wall-clock guard per solve (12× for the multistart); a timeout is a failure
_FLAGS = ("warm_timeout", "warm_error", "warm_capped", "multi_timeout", "multi_error",
          "multi_n_capped", "chain_timeout", "chain_error")  # fmt: skip
QUANTILES = (2.5, 16, 50, 84, 97.5)
# § 2's listed counts (from a month-rounded parse when the plan was written); the rule is the
# window on the epoch, which the code applies; differences are reported in R1.
PLAN_COUNTS = {
    1712: 1,
    1715: 3,
    1750: 2,
    1753: 1,
    1756: 1,
    1764: 1,
    1769: 8,
    1771: 1,
    1783: 19,
    1786: 6,
    1789: 15,
    1792: 14,
    1795: 11,
    1798: 6,
    1801: 2,
    1807: 2,
    1810: 1,
    1828: 6,
    1831: 5,
    1834: 13,
    1837: 17,
    1840: 10,
}
PLAN_EMPTY_WINDOWS = (1780, 1804, 1813, 1816, 1819, 1822, 1825)

_wrap = R._wrap


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with PROGRESS_LOG.open("a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def seed_for(model: str, hyp: int) -> list[int]:
    return [*SEED_NAMESPACE, MODELS.index(model), hyp]


def noise_matrix(model: str, hyp: int, sigma: np.ndarray, n: int = N_DRAWS) -> np.ndarray:
    """n × 29 noise in arcsec: standard normal ((a), (b)) or Student-t ν = 4 ((c)) in units of
    σ_j, one generator per model × hypothesis."""
    rng = np.random.default_rng(seed_for(model, hyp))
    z = (
        rng.standard_t(NU_T, (n, len(sigma)))
        if model == "c"
        else rng.standard_normal((n, len(sigma)))
    )
    return z * sigma[None, :]


# ------------------------------------------------------------------ R1: the frozen inputs


def set_d_window_counts(m: A.Memoir) -> dict[int, dict]:
    """Set D observations (the step-3 fit set) per normal-place window: ±1.518 yr about a
    modern place, ±1 yr about an ancient one, the epoch from its TT Julian Date."""
    from fork.model import observation as om

    fit, _ = pa.adams_rows(om.read_rows())
    obs = {}
    for r in fit:
        obs.setdefault(
            r["obs_id"], (2000.0 + (float(r["epoch_tt_jd"]) - 2451545.0) / 365.25, om.group_of(r))
        )
    out = {}
    assigned = 0
    for p in m.places:
        if p["year"] == 1690:
            continue
        y = A.EPOCH_YEAR + A.UNIT_YR * float(p["t_h1"])
        half = WINDOW_MODERN_YR if p["kind"] == "modern" else WINDOW_ANCIENT_YR
        sel = [g for (yy, g) in obs.values() if abs(yy - y) <= half]
        assigned += len(sel)
        out[p["year"]] = {
            "epoch_year": y,
            "n": len(sel),
            "groups": sorted(set(sel)),
            "group_used": A.set_d_group_of_place(p),
        }
    return {"places": out, "n_observations": len(obs), "n_assigned": assigned}


def step3_group_rms() -> dict:
    """The RMS of step 3's H1(G) best-fit residuals about the fitted group offsets, RA only."""
    from fork.inference.nested import NestedResult
    from fork.model.observation import GROUPS

    obs, ctx, _inc, _pq, _box, _test = pa.dataset("h1_g")
    res = [NestedResult.load(j.path(".pkl")) for j in pa.all_jobs(("h1_g",))]
    best = max(res, key=lambda r: r.max_logl)
    th = best.space.theta(best.samples[int(np.argmax(best.results["logl"]))])
    resid = obs.residuals_arcsec(*obs.predict(th, ctx))
    r0 = resid - obs.best_offsets_arcsec(resid)[obs.group_index]
    out = {}
    for gi, g in enumerate(GROUPS):  # groups carry the quantity: "greenwich:ra"
        if not g.endswith(":ra"):
            continue
        sel = obs.group_index == gi
        if sel.sum():
            out[g[:-3]] = {"rms_arcsec": float(np.sqrt(np.mean(r0[sel] ** 2))), "n": int(sel.sum())}
    return out


def _rows_sigma_vector(m: A.Memoir, sig: dict[int, float]) -> np.ndarray:
    return np.array([sig[p["year"]] for p in m.places if p["year"] != 1690])


# ------------------------------------------------------------------ the two methods


def chain_solve(ms: A.Memoir, hyp: int, rows, w: np.ndarray, s: A.Setup) -> dict:
    """His own procedure on the § 9 column now in ``ms``: the C2 chain, every positive-mass
    root completed, the one with the smallest χ² under D1(c)'s weights kept (as D3(i))."""
    run = A.run_hypothesis_1 if hyp == 1 else A.run_hypothesis_2
    try:
        ch = run(ms, "C2")
    except (ValueError, ZeroDivisionError, KeyError, FloatingPointError, OverflowError):
        return {"converged": False, "n_roots": 0, "n_positive": 0}
    th_, _S, _e, beta, k, r_cos, r_sin, _R, eps, n_prime, prec = A._setup(ms, hyp, "print")
    eq1 = ch.steps[(22 if hyp == 1 else 39, "1")]
    best = None
    n_pos = 0
    for root in ch.roots:
        if root.get("m_prime", 0.0) <= 0:
            continue
        n_pos += 1
        try:
            cs = A.complete_solution(
                hyp,
                root["theta"],
                root["P"],
                root["Q"],
                th_.c,
                eq1,
                k,
                r_cos,
                r_sin,
                beta,
                eps,
                n_prime,
                prec,
            )
            corr = A.element_corrections(ch, ms, cs)
        except (ValueError, ZeroDivisionError, KeyError, FloatingPointError, OverflowError):
            continue
        x = A.x_of(cs, corr)
        chi2 = A.chi2_of(rows, w, s, x)
        if math.isfinite(chi2) and (best is None or chi2 < best[0]):
            best = (chi2, x)
    if best is None:
        return {"converged": False, "n_roots": len(ch.roots), "n_positive": n_pos}
    x = best[1].copy()
    x[4] %= 360.0
    x[7] %= 360.0
    return {
        "converged": True,
        "x": x,
        "chi2": best[0],
        "n_roots": len(ch.roots),
        "n_positive": n_pos,
    }


def longitudes_of(x: np.ndarray, eps_c: float, n_prime_deg: float) -> dict:
    """The solution's mean and true longitude at 1846.762 (ecl_j2000, the data's frame)."""
    return R._his_model_longitudes(float(x[4]), float(x[6]) / 20.0, float(x[7]), eps_c, n_prime_deg)


# ------------------------------------------------------------------ prerun


def _real_data_offsets(m: A.Memoir, hyp: int, T_of_date: float) -> dict:
    """O of his real-data solutions (of date, as step 3b reported them) by the three methods."""
    s = A.setup_from_print(m, hyp)
    rows = A.place_rows(m, s)
    w = A.weights_for(m, hyp, "c", rows)
    d1c = A.solve_least_squares(rows, w, s, A.printed_solution(m, hyp))
    run = A.run_hypothesis_1 if hyp == 1 else A.run_hypothesis_2
    ch = A.chain_with_print_fallback(m, hyp, run(m, "C2"))
    pr = A.printed_solution(m, hyp)
    out = {
        "d1c": d1c["true_lon_1846"],
        "C2": A.true_longitude(
            ch.solution.mean_lon_1846, ch.solution.eccentricity, ch.solution.varpi_prime_1846
        ),
        "printed": A.true_longitude(
            pr["mean_lon_1846"], pr["eccentricity"], pr["varpi_prime_1846"]
        ),
    }
    return {"true_lon_1846_of_date": out, "O": {k: _wrap(v - T_of_date) for k, v in out.items()}}


def prerun() -> dict:
    """The frozen inputs: the regenerated synthetic column, the counts, the group RMS, the σ
    tables, the clean solutions by both methods for both hypotheses, Δ_obs."""
    m = A.load_memoir()
    syn = R.synthetic_normal_places(m)
    T = R.perturber_lon_1846(syn["truth_theta"], syn["ctx"])
    col = {int(k): float(v) for k, v in syn["oc_mean_arcsec"].items()}
    check = None
    if R.CHEAP_JSON.exists():
        cached = json.loads(R.CHEAP_JSON.read_text())["d3"]["synthetic_places"]["oc_mean_arcsec"]
        check = max(abs(col[int(k)] - v) for k, v in cached.items())
    counts = set_d_window_counts(m)
    grp = step3_group_rms()
    group_rms = {g: v["rms_arcsec"] for g, v in grp.items()}
    n_by_year = {y: v["n"] for y, v in counts["places"].items()}
    eps_c = syn["c_fit"]["eps_1810_ecl_j2000_deg"]
    out = {
        "plan": "docs/notes/P1B_adams_noise.md",
        "memoir_sha256": R.memoir_hashes(),
        "frozen_inputs": {k: {"path": p, "sha256": h} for k, (p, h) in pa.FROZEN.items()},
        "truth": {"neptune_1846": T, "theta": syn["truth_theta"]},
        "synthetic_column_arcsec": col,
        "synthetic_column_max_abs_diff_vs_step3b_cache": check,
        "c_fit": {k: v for k, v in syn["c_fit"].items() if k != "theta"},
        "window_counts": counts,
        "counts_vs_plan": {
            y: {"rule": v["n"], "plan": PLAN_COUNTS.get(y, 0)}
            for y, v in counts["places"].items()
            if v["n"] != PLAN_COUNTS.get(y, 0)
        },
        "empty_windows": sorted(y for y, v in counts["places"].items() if v["n"] == 0),
        "step3_h1g_group_rms_ra": grp,
        "scale_a": A.SCALE_A,
        "sigma": {},
        "clean": {},
        "real_data": {},
        "delta_obs": {},
    }
    for hyp in HYPS:
        rms = A.residual_rms(m, hyp)
        out["sigma"][f"H{hyp}"] = {
            "a_raw_rms": rms,
            "a": A.noise_sigmas(m, "a", hyp=hyp),
            "b": A.noise_sigmas(m, "b", n_by_year, hyp, group_rms),
        }
        out["sigma"][f"H{hyp}"]["c"] = out["sigma"][f"H{hyp}"]["b"]
        s = A.setup_from_print(m, hyp)
        rows = A.place_rows(m, s, oc=col)
        w = A.weights_for(m, hyp, "c", rows)
        ls = A.solve_least_squares(rows, w, s, A.printed_solution(m, hyp))
        x_clean = np.array([ls["x"][k] for k in A.PARAMS])
        warm = A.solve_warm(rows, w, s, x_clean)
        lon = longitudes_of(x_clean, eps_c, s.n_prime_deg)
        ms = copy.deepcopy(m)
        for p in ms.places:
            p["oc_mean"] = col[p["year"]]
        chn = chain_solve(ms, hyp, rows, w, s)
        lon_ch = longitudes_of(chn["x"], eps_c, s.n_prime_deg) if chn["converged"] else None
        out["clean"][f"H{hyp}"] = {
            "alpha": 0.5 if hyp == 1 else 0.515,
            "n_prime_deg": s.n_prime_deg,
            "eps_c_1810_ecl_j2000": eps_c,
            "d1c": {
                "x": dict(zip(A.PARAMS, x_clean.tolist(), strict=True)),
                "chi2": ls["chi2"],
                "sigma_scaled": ls["sigma"],
                "sigma_unscaled": dict(zip(A.PARAMS, warm["sigma_unscaled"].tolist(), strict=True)),
                "warm_restart_chi2_diff": warm["chi2"] - ls["chi2"],
                "n_distinct_minima": ls["n_distinct_minima"],
                "longitudes": lon,
                "O_clean": _wrap(lon["true_lon_1846_ecl_j2000"] - T["ecl_j2000"]),
            },
            "chain": {
                "converged": chn["converged"],
                "n_roots": chn["n_roots"],
                "n_positive_mass_roots": chn["n_positive"],
                "x": dict(zip(A.PARAMS, chn["x"].tolist(), strict=True))
                if chn["converged"]
                else None,
                "chi2_c_weights": chn.get("chi2"),
                "longitudes": lon_ch,
                "O_clean": None
                if lon_ch is None
                else _wrap(lon_ch["true_lon_1846_ecl_j2000"] - T["ecl_j2000"]),
            },
        }
        real = _real_data_offsets(m, hyp, T["of_date"])
        out["real_data"][f"H{hyp}"] = real
        Oc = out["clean"][f"H{hyp}"]
        out["delta_obs"][f"H{hyp}"] = {
            "d1c": {
                "primary": real["O"]["d1c"] - Oc["d1c"]["O_clean"],
                "printed": real["O"]["printed"] - Oc["d1c"]["O_clean"],
            },
            "chain": None
            if Oc["chain"]["O_clean"] is None
            else {
                "primary": real["O"]["C2"] - Oc["chain"]["O_clean"],
                "printed": real["O"]["printed"] - Oc["chain"]["O_clean"],
            },
        }
    R._dump(PRERUN_JSON, out)
    _log(f"prerun written: {PRERUN_JSON}")
    return out


# ------------------------------------------------------------------ draws (process pool)

_STATE: dict = {}


def _init_worker(prerun_path: str) -> None:
    pre = json.loads(Path(prerun_path).read_text())
    m = A.load_memoir()
    col = {int(k): v for k, v in pre["synthetic_column_arcsec"].items()}
    for hyp in HYPS:
        s = A.setup_from_print(m, hyp)
        rows = A.place_rows(m, s, oc=col)
        w = A.weights_for(m, hyp, "c", rows)
        c = pre["clean"][f"H{hyp}"]
        ms = copy.deepcopy(m)
        for p in ms.places:
            p["oc_mean"] = col[p["year"]]
        _STATE[hyp] = {
            "m": m,
            "s": s,
            "rows": rows,
            "w": w,
            "ms": ms,
            "x_clean": np.array([c["d1c"]["x"][k] for k in A.PARAMS]),
            "eps_c": c["eps_c_1810_ecl_j2000"],
            "n_prime": c["n_prime_deg"],
        }


def _noisy_rows(rows, noise: np.ndarray):
    return [
        (y, k, t, dict(eq, lhs=eq["lhs"] + float(dn)))
        for (y, k, t, eq), dn in zip(rows, noise, strict=True)
    ]


_FIELDS = ("theta", "m_prime", "e_units", "beta_prime", "chi2", "mean_lon", "true_lon")


def _pack(res: dict | None, eps_c: float, n_prime: float) -> list[float]:
    if res is None or not res.get("converged"):
        return [float("nan")] * len(_FIELDS) + [0.0, 0.0]
    x = res["x"]
    lon = longitudes_of(x, eps_c, n_prime)
    return [
        x[4],
        x[5],
        x[6],
        x[7],
        res["chi2"],
        lon["mean_lon_1846_ecl_j2000"],
        lon["true_lon_1846_ecl_j2000"],
        1.0,
        float(res.get("bound", False)),
    ]


class _Watchdog(Exception):
    pass


def _alarm(_signum, _frame):
    raise _Watchdog()


def _guarded(fn, seconds: float):
    """Run ``fn`` under a wall-clock watchdog; returns (result | None, status) with status
    "ok", "timeout" or the exception's one-line description."""
    import signal

    signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return fn(), "ok"
    except _Watchdog:
        return None, "timeout"
    except Exception as exc:  # recorded per draw, never lost with the chunk
        return None, f"{type(exc).__name__}: {exc}"[:200]
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def chunk_path(hyp: int, model: str, i0: int) -> Path:
    return CHUNK_DIR / f"h{hyp}_{model}_{i0:05d}.npz"


def _task(args) -> tuple:
    """One chunk: every realisation through the warm start, the multistart (first ``n_multi``)
    and the chain, each stage under the watchdog; the chunk is written to disk by the worker."""
    hyp, model, i0, noise, n_multi = args
    out_path = chunk_path(hyp, model, i0)
    if out_path.exists():
        return hyp, model, i0, "cached", []
    st = _STATE[hyp]
    s, w = st["s"], st["w"]
    cap = NFEV_CAP_C if model == "c" else NFEV_DEFAULT
    nan_row = [float("nan")] * len(_FIELDS) + [0.0, 0.0]
    d1c, chain, multi, roots, flags, errors = [], [], [], [], [], []
    for i in range(noise.shape[0]):
        rows = _noisy_rows(st["rows"], noise[i])
        fl = [0.0] * len(_FLAGS)  # timeouts/errors per stage, capped solves
        r, status = _guarded(
            lambda rows=rows: A.solve_warm(rows, w, s, st["x_clean"], max_nfev=cap), WATCHDOG_S
        )
        if status != "ok":
            fl[_FLAGS.index("warm_timeout" if status == "timeout" else "warm_error")] = 1.0
            errors.append(f"draw {i0 + i} warm: {status}")
        elif not r["converged"] and r["nfev"] >= cap:
            fl[_FLAGS.index("warm_capped")] = 1.0
        d1c.append(_pack(r, st["eps_c"], st["n_prime"]) if r else nan_row)
        if i < n_multi:
            r, status = _guarded(
                lambda rows=rows: A.solve_multistart(rows, w, s, st["x_clean"], max_nfev=cap),
                12 * WATCHDOG_S,
            )
            if status != "ok":
                fl[_FLAGS.index("multi_timeout" if status == "timeout" else "multi_error")] = 1.0
                errors.append(f"draw {i0 + i} multistart: {status}")
            else:
                fl[_FLAGS.index("multi_n_capped")] = float(r.get("n_capped", 0))
            multi.append(_pack(r, st["eps_c"], st["n_prime"]) if r else nan_row)
        ms = st["ms"]
        j = 0
        for p in ms.places:
            if p["year"] == 1690:
                continue
            p["oc_mean"] = st["rows"][j][3]["lhs"] + float(noise[i][j])
            j += 1
        r, status = _guarded(lambda rows=rows, ms=ms: chain_solve(ms, hyp, rows, w, s), WATCHDOG_S)
        if status != "ok":
            fl[_FLAGS.index("chain_timeout" if status == "timeout" else "chain_error")] = 1.0
            errors.append(f"draw {i0 + i} chain: {status}")
            r = {"converged": False, "n_roots": 0, "n_positive": 0}
        chain.append(_pack(r, st["eps_c"], st["n_prime"]))
        roots.append([r["n_roots"], r["n_positive"]])
        flags.append(fl)
    CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tmp.npz")
    np.savez(tmp, d1c=np.array(d1c), chain=np.array(chain), roots=np.array(roots),
             multi=np.array(multi), flags=np.array(flags))  # fmt: skip
    os.replace(tmp, out_path)
    return hyp, model, i0, "done", errors


def draws_path(hyp: int, model: str) -> Path:
    return CACHE / f"draws_h{hyp}_{model}.npz"


def assemble(n_draws: int = N_DRAWS) -> list[tuple[int, str]]:
    """Join each complete cell's chunk files into its draws file; returns the cells written."""
    written = []
    for hyp in HYPS:
        for model in MODELS:
            paths = [chunk_path(hyp, model, i0) for i0 in range(0, n_draws, CHUNK)]
            if not all(p.exists() for p in paths):
                continue
            parts = [np.load(p) for p in paths]
            arrays = {k: np.concatenate([z[k] for z in parts if z[k].size])
                      for k in ("d1c", "chain", "roots", "multi", "flags")
                      if any(z[k].size for z in parts)}  # fmt: skip
            n_multi = arrays["multi"].shape[0] if "multi" in arrays else 0
            np.savez(draws_path(hyp, model), fields=np.array((*_FIELDS, "converged", "bound")),
                     flag_names=np.array(_FLAGS), multistart_all=np.array(n_multi == n_draws),
                     **arrays)  # fmt: skip
            written.append((hyp, model))
    return written


def run_draws(n_workers: int, multistart_all: bool = False, n_draws: int = N_DRAWS) -> None:
    """Every model × hypothesis: the noise matrix from its seed, both methods on every
    realisation, the multistart on the first N_MULTI (or on all with ``multistart_all``).
    Cells (a) and (b) of both hypotheses run first, (c) last; each chunk is a file on disk
    (a rerun skips the chunks already there); progress is logged per chunk."""
    from concurrent.futures import as_completed

    from fork import ledger

    pre = json.loads(PRERUN_JSON.read_text())
    m = A.load_memoir()
    t0 = time.time()
    tasks = []
    for model in MODELS:  # (a), (b) for both hypotheses, then (c)
        for hyp in HYPS:
            sig = _rows_sigma_vector(
                m, {int(k): v for k, v in pre["sigma"][f"H{hyp}"][model].items()}
            )
            noise = noise_matrix(model, hyp, sig, n_draws)
            n_multi_total = n_draws if multistart_all else min(N_MULTI, n_draws)
            for i0 in range(0, n_draws, CHUNK):
                tasks.append(
                    (hyp, model, i0, noise[i0 : i0 + CHUNK], max(0, min(CHUNK, n_multi_total - i0)))
                )
    _log(
        f"draws: {len(tasks)} chunks of {CHUNK} over {n_workers} workers "
        f"(multistart on {'all' if multistart_all else N_MULTI} per cell; order a, b, c)"
    )
    done = 0
    with ProcessPoolExecutor(
        n_workers, initializer=_init_worker, initargs=(str(PRERUN_JSON),)
    ) as ex:
        futures = {ex.submit(_task, t): t[:3] for t in tasks}
        for fut in as_completed(futures):
            hyp, model, i0 = futures[fut]
            done += 1
            try:
                _h, _m, _i, state, errors = fut.result()
            except Exception as exc:  # a chunk lost whole: logged, the others go on
                _log(f"CHUNK FAILED h{hyp}/{model}/{i0}: {type(exc).__name__}: {exc}")
                continue
            for e in errors[:5]:
                _log(f"h{hyp}/{model}: {e}")
            _log(f"chunk {done}/{len(tasks)} h{hyp}/{model}/{i0:05d} {state}, "
                 f"{(time.time() - t0) / 60:.1f} min")  # fmt: skip
    written = assemble(n_draws)
    ct = os.times()
    ledger.append(
        "p1b_adams_noise/draws",
        "fit",
        n_likelihood_calls=None,
        n_integrations=None,
        sim_years=None,
        wall_s=time.time() - t0,
        cpu_s=ct.children_user + ct.children_system + ct.user + ct.system,
        n_workers=n_workers,
        estimated=False,
        note=f"{n_draws} draws × {len(MODELS)} models × {len(HYPS)} hypotheses × 2 methods"
        f"{' + multistart on all' if multistart_all else f' + multistart on {N_MULTI}'}"
        "; per-chunk files",
    )
    _log(f"draws done in {(time.time() - t0) / 60:.1f} min; cells assembled: {written}")


# ------------------------------------------------------------------ summarize


def _col(arr: np.ndarray, name: str) -> np.ndarray:
    return arr[:, ((*_FIELDS, "converged", "bound")).index(name)]


def p_value(delta: np.ndarray, d_obs: float, n_failed: int = 0) -> dict:
    """One-sided P(Δ ≥ Δ_obs) over the converged draws with its binomial standard error; the
    upper bound counts every failed draw as exceeding."""
    n = int(delta.size)
    k = int(np.sum(delta >= d_obs))
    p = k / n if n else float("nan")
    out = {
        "n_exceeding": k,
        "n": n,
        "P": p,
        "se": math.sqrt(p * (1 - p) / n) if n and 0 < p < 1 else (1 / n if n else float("nan")),
    }
    if n_failed:
        out["P_upper_bound"] = (k + n_failed) / (n + n_failed)
    return out


def reading(P_a: float, P_b: float) -> dict:
    """§ 4's rows on the two bracketing models."""
    if P_a >= P_NOISE_SUFFICES and P_b >= P_NOISE_SUFFICES:
        return {
            "row": 1,
            "verdict": "noise suffices",
            "text": "A shift as large as his is produced by independent observational noise in "
            "his layout at least one time in twenty under both models; his nearness to Neptune "
            "on real data is within luck; no systematic is needed (none is excluded).",
        }
    if P_a < P_SYSTEMATIC and P_b < P_SYSTEMATIC:
        return {
            "row": 2,
            "verdict": "a systematic is required",
            "text": "Independent noise of either size does not carry his method from where clean "
            "data put it to where his real data put it. Candidates, not ranked: Bouvard's 1821 "
            "tables and his § 7 corrections, the construction of the normal places, a reduction "
            "error common to a series; not attributed further without Bouvard's tables.",
        }
    if P_a >= P_NOISE_SUFFICES > P_b:
        which = "(a) admits it, (b) does not"
    elif P_b >= P_NOISE_SUFFICES > P_a:
        which = "(b) admits it, (a) does not"
    else:
        which = "neither model reaches 0.05 and at least one is above 0.003"
    return {
        "row": 3,
        "verdict": "weak",
        "text": f"P = {P_a:.4f} under (a) and {P_b:.4f} under (b): {which}; no single statement.",
    }


def reading_interval(I_a: list[float], I_b: list[float]) -> dict:
    """§ 4's rows when a cell's P is an interval (flagged draws above 1 %): a cell whose
    interval crosses a reading threshold reads "weak"."""
    for name, (lo, hi) in (("a", I_a), ("b", I_b)):
        for thr in (P_SYSTEMATIC, P_NOISE_SUFFICES):
            if lo < thr <= hi:
                return {
                    "row": 3,
                    "verdict": "weak",
                    "text": f"P under ({name}) lies in [{lo:.4f}, {hi:.4f}] with the flagged draws "
                    f"counted either way, which crosses {thr}; no single statement.",
                }
    return reading(I_a[0], I_b[0])


def multistart_baselines(pre: dict) -> dict:
    """Δ_obs for D1(c) recomputed with the draws' own 12-start multistart on the clean column
    and on his real § 9 data (the reading's primary method, both ends by one solver)."""
    m = A.load_memoir()
    col = {int(k): v for k, v in pre["synthetic_column_arcsec"].items()}
    T = pre["truth"]["neptune_1846"]
    out = {}
    for hyp in HYPS:
        c = pre["clean"][f"H{hyp}"]
        s = A.setup_from_print(m, hyp)
        xc = np.array([c["d1c"]["x"][k] for k in A.PARAMS])
        rows = A.place_rows(m, s, oc=col)
        r = A.solve_multistart(rows, A.weights_for(m, hyp, "c", rows), s, xc)
        lon = longitudes_of(r["x"], c["eps_c_1810_ecl_j2000"], s.n_prime_deg)
        O_clean = _wrap(lon["true_lon_1846_ecl_j2000"] - T["ecl_j2000"])
        rows = A.place_rows(m, s)
        q = A.solve_multistart(rows, A.weights_for(m, hyp, "c", rows), s, xc)
        x = q["x"]
        mean = ((s.eps - x[4]) % 360.0 + 12.0 * s.n_prime_deg + s.prec) % 360.0
        true = A.true_longitude(mean, x[6] / 20.0, (s.eps - x[7] + s.prec) % 360.0)
        O_real = _wrap(true - T["of_date"])
        pr = pre["real_data"][f"H{hyp}"]["O"]["printed"]
        out[f"H{hyp}"] = {
            "clean": {"chi2": r["chi2"], "theta": float(r["x"][4]), "O": O_clean},
            "real": {"chi2": q["chi2"], "theta": float(x[4]), "O": O_real},
            "delta_obs": {"primary": O_real - O_clean, "printed": pr - O_clean},
            "minus_prerun_delta_obs": O_real
            - O_clean
            - pre["delta_obs"][f"H{hyp}"]["d1c"]["primary"],
        }
    return out


def _quant(v: np.ndarray) -> dict:
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {"n": 0}
    q = np.percentile(v, QUANTILES)
    return {
        "n": int(v.size),
        **{f"q{p}": float(x) for p, x in zip(QUANTILES, q, strict=True)},
        "mean": float(np.mean(v)),
        "sd": float(np.std(v, ddof=1)) if v.size > 1 else 0.0,
    }


def _method_summary(
    arr: np.ndarray,
    T_ecl: float,
    O_clean: float,
    theta_clean: float,
    d_obs: dict,
    n_total: int,
    flagged: np.ndarray | None = None,
) -> dict:
    """``flagged``: per draw, a timeout, an exception or a capped solve for this method; a
    non-converged draw is flagged too. Above FAIL_TRIGGER of the cell, P is the interval
    [P over the converged draws, P with every flagged draw counted as exceeding]."""
    conv = _col(arr, "converged") > 0.5
    flagged = ~conv if flagged is None else (flagged | ~conv)
    n_failed = int(n_total - conv.sum())
    true_lon = _col(arr, "true_lon")[conv]
    off = np.array([_wrap(v - T_ecl) for v in true_lon])
    delta = off - O_clean
    theta = _col(arr, "theta")[conv]
    escaped = np.abs(np.array([_wrap(v - theta_clean) for v in theta])) > ESCAPE_DEG
    bound = _col(arr, "bound")[conv] > 0.5
    out = {
        "n_total": n_total,
        "n_converged": int(conv.sum()),
        "fraction_failed": n_failed / n_total,
        "fraction_escaped": float(escaped.mean()) if conv.any() else float("nan"),
        "fraction_bound": float(bound.mean()) if conv.any() else float("nan"),
        "delta": _quant(delta),
        "O": _quant(off),
        "mean_lon_1846_ecl_j2000": _quant(_col(arr, "mean_lon")[conv]),
        "theta": _quant(theta),
        "theta_sd_deg": float(np.std(theta, ddof=1)) if conv.sum() > 1 else float("nan"),
        "m_prime": _quant(_col(arr, "m_prime")[conv]),
        "e_prime": _quant(_col(arr, "e_units")[conv] / 20.0),
        "P": {},
    }
    out["n_flagged"] = int(flagged.sum())
    out["fraction_flagged"] = float(flagged.mean())
    over = out["fraction_flagged"] > FAIL_TRIGGER
    for target, dv in d_obs.items():
        rec = {"delta_obs": dv, **p_value(delta, dv)}
        rec["P_interval"] = [rec["P"], rec["P"]]
        if over:
            k_clean = int(np.sum(delta[~flagged[conv]] >= dv))
            rec["P_interval"] = [rec["P"], (k_clean + int(flagged.sum())) / n_total]
        out["P"][target] = rec
    return out


def multistart_fraction(warm: np.ndarray, multi: np.ndarray) -> dict:
    n = multi.shape[0]
    cw = _col(warm[:n], "converged") > 0.5
    cm = _col(multi, "converged") > 0.5
    both = cw & cm
    diff = np.abs(
        np.array(
            [
                _wrap(a - b)
                for a, b in zip(
                    _col(multi, "true_lon")[both], _col(warm[:n], "true_lon")[both], strict=True
                )
            ]
        )
    )
    lower = _col(multi, "chi2")[both] < _col(warm[:n], "chi2")[both] - 1e-6
    return {
        "n_checked": int(both.sum()),
        "n_differing": int(np.sum(diff > MULTI_DIFF_DEG)),
        "fraction_differing": float(np.mean(diff > MULTI_DIFF_DEG)) if both.any() else float("nan"),
        "fraction_multistart_lower_chi2": float(np.mean(lower)) if both.any() else float("nan"),
        "max_diff_deg": float(diff.max()) if both.any() else float("nan"),
    }


def summarize(png: bool = True) -> dict:
    pre = json.loads(PRERUN_JSON.read_text())
    T_ecl = pre["truth"]["neptune_1846"]["ecl_j2000"]
    out = {
        "plan": pre["plan"],
        "prerun": {
            k: pre[k]
            for k in (
                "memoir_sha256",
                "frozen_inputs",
                "synthetic_column_max_abs_diff_vs_step3b_cache",
                "window_counts",
                "counts_vs_plan",
                "empty_windows",
                "step3_h1g_group_rms_ra",
                "scale_a",
                "sigma",
                "clean",
                "real_data",
                "delta_obs",
            )
        },
        "truth": pre["truth"]["neptune_1846"],
        "results": {},
        "multistart_check": {},
        "multistart_rerun_required": False,
        "readings": {},
        "cells_missing": [],
        "multistart_baselines": multistart_baselines(pre),
    }
    for hyp in HYPS:
        H = f"H{hyp}"
        c = pre["clean"][H]
        th_clean = {
            "d1c": c["d1c"]["x"]["theta"],
            "chain": c["chain"]["x"]["theta"] if c["chain"]["x"] else float("nan"),
        }
        O_clean = {"d1c": c["d1c"]["O_clean"], "chain": c["chain"]["O_clean"]}
        for model in MODELS:
            if not draws_path(hyp, model).exists():
                out["cells_missing"].append(f"{H}/{model}")
                continue
            z = np.load(draws_path(hyp, model))
            n_total = z["d1c"].shape[0]
            rec = {
                "seed": seed_for(model, hyp),
                "n_draws": n_total,
                "sigma_arcsec": pre["sigma"][H][model],
            }
            if model == "a":
                rec["sigma_raw_rms_arcsec"] = pre["sigma"][H]["a_raw_rms"]
                rec["note"] = "conditioned on his own fit; scaled by √(29/21)"
            if "flags" in z:
                fl = z["flags"]
                rec["flags"] = {
                    str(nm): int(np.sum(fl[:, k] > 0)) for k, nm in enumerate(z["flag_names"])
                }
            if "multi" in z:
                rec["multistart"] = multistart_fraction(z["d1c"], z["multi"][:N_MULTI])
                if z["multi"].shape[0] > N_MULTI:
                    rec["multistart_all_draws"] = multistart_fraction(z["d1c"], z["multi"])
                out["multistart_check"][f"{H}/{model}"] = rec["multistart"]
            use_multi = (
                bool(z["multistart_all"]) and "multi" in z and z["multi"].shape[0] == n_total
            )
            for method in METHODS:
                if O_clean[method] is None:
                    continue
                arr = z["multi"] if (method == "d1c" and use_multi) else z[method]
                d_obs = pre["delta_obs"][H][method]
                if method == "d1c" and use_multi:
                    d_obs = out["multistart_baselines"][H]["delta_obs"]
                fl_m = None
                if "flags" in z:
                    names = [str(v) for v in z["flag_names"]]
                    pref = (
                        "multi_"
                        if (method == "d1c" and use_multi)
                        else ("warm_" if method == "d1c" else "chain_")
                    )
                    cols = [k for k, nm in enumerate(names) if nm.startswith(pref)]
                    fl_m = np.any(z["flags"][:, cols] > 0, axis=1)
                oc = O_clean[method]
                if method == "d1c" and use_multi:
                    oc = out["multistart_baselines"][H]["clean"]["O"]
                rec[method] = _method_summary(
                    arr, T_ecl, oc, th_clean[method], d_obs, n_total, fl_m
                )
                if method == "d1c":
                    rec[method]["solver"] = "multistart (every draw)" if use_multi else "warm start"
                    if model == "a":
                        rec[method]["sd_check"] = {
                            "theta_sd_draws_deg": rec[method]["theta_sd_deg"],
                            "theta_sd_predicted_deg": pre["scale_a"]
                            * c["d1c"]["sigma_unscaled"]["theta"],
                            "step3b_sigma_theta_scaled_deg": c["d1c"]["sigma_scaled"]["theta"],
                        }
                if method == "chain":
                    roots = z["roots"]
                    rec[method]["roots"] = {
                        "mean_n_roots": float(roots[:, 0].mean()),
                        "mean_positive_mass": float(roots[:, 1].mean()),
                    }
            out["results"][f"{H}/{model}"] = rec
    out["multistart_rerun_required"] = any(
        v["fraction_differing"] > MULTI_TRIGGER for v in out["multistart_check"].values()
    )
    for hyp in HYPS:
        H = f"H{hyp}"
        for method in METHODS:
            if (
                f"{H}/a" not in out["results"]
                or f"{H}/b" not in out["results"]
                or method not in out["results"][f"{H}/a"]
            ):
                continue
            for target in ("primary", "printed"):
                Pa = out["results"][f"{H}/a"][method]["P"][target]["P"]
                Pb = out["results"][f"{H}/b"][method]["P"][target]["P"]
                Ia = out["results"][f"{H}/a"][method]["P"][target]["P_interval"]
                Ib = out["results"][f"{H}/b"][method]["P"][target]["P_interval"]
                rc = out["results"].get(f"{H}/c")
                Pc = rc[method]["P"][target]["P"] if rc else None
                out["readings"][f"{H}/{method}/{target}"] = {
                    "P_a": Pa,
                    "P_b": Pb,
                    "P_c_sensitivity": Pc,
                    "P_a_interval": Ia,
                    "P_b_interval": Ib,
                    **reading_interval(Ia, Ib),
                    "headline": hyp == 1 and method == "d1c" and target == "primary",
                }
    R._dump(OUT_JSON, out)
    if png:
        draw(out, pre, OUT_PNG)
    return out


def draw(out: dict, pre: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    T_ecl = pre["truth"]["neptune_1846"]["ecl_j2000"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), sharex=True)
    hues = {"a": "#1f6fd1", "b": "#d1521f", "c": "#777777"}
    labels = {
        "a": "(a) his § 52 scatter × 1.175",
        "b": "(b) set D scatter / √N",
        "c": "(c) Student-t ν = 4, (b)'s scales",
    }
    titles = {
        "d1c": "warm-start least squares (his method as D1(c))",
        "chain": "his own elimination chain",
    }
    for i, hyp in enumerate(HYPS):
        H = f"H{hyp}"
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            c = pre["clean"][H]
            O_clean = c[method]["O_clean"]
            if O_clean is None:
                ax.set_axis_off()
                continue
            for model in MODELS:
                if not draws_path(hyp, model).exists():
                    continue
                z = np.load(draws_path(hyp, model))
                arr = (
                    z["multi"]
                    if (
                        method == "d1c"
                        and bool(z["multistart_all"])
                        and "multi" in z
                        and z["multi"].shape[0] == z["d1c"].shape[0]
                    )
                    else z[method]
                )
                conv = _col(arr, "converged") > 0.5
                delta = np.array([_wrap(v - T_ecl) for v in _col(arr, "true_lon")[conv]]) - O_clean
                ax.hist(
                    delta,
                    bins=np.arange(-15, 20.01, 0.5),
                    histtype="step",
                    lw=1.6 if model != "c" else 1.0,
                    ls="-" if model != "c" else ":",
                    color=hues[model],
                    label=labels[model],
                    density=True,
                )
            d = pre["delta_obs"][H][method]
            ax.axvline(d["primary"], color="black", lw=1.2)
            ax.axvline(d["printed"], color="black", lw=0.8, ls="--")
            ps = out["results"]
            txt = "\n".join(
                f"P({model}) = {ps[f'{H}/{model}'][method]['P']['primary']['P']:.4f}"
                for model in MODELS
                if f"{H}/{model}" in ps
            )
            ax.text(0.02, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9)
            ax.text(
                d["primary"],
                ax.get_ylim()[1] * 0.98,
                f" Δ_obs {d['primary']:+.1f}°",
                fontsize=8,
                va="top",
            )
            ax.set_title(f"Hypothesis {'I' if hyp == 1 else 'II'}: {titles[method]}", fontsize=10)
            ax.grid(alpha=0.25)
            if i == 1:
                ax.set_xlabel("Δ = (offset from Neptune with noise) − (noise-free offset), degrees")
            if j == 0:
                ax.set_ylabel("density of realisations")
    axes[0, 0].legend(fontsize=8, loc="upper right", frameon=False)
    fig.suptitle(
        "Step 3c: does independent noise in his layout carry his method as far as his real "
        "data did? "
        "(solid line his real-data solution by the same method, dashed his printed one)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prerun")
    p_d = sub.add_parser("draws")
    p_d.add_argument("--n-workers", type=int, default=min(12, os.cpu_count() or 1))
    p_d.add_argument("--multistart-all", action="store_true")
    p_d.add_argument("--n-draws", type=int, default=N_DRAWS, help="smoke runs only")
    p_s = sub.add_parser("summarize")
    p_s.add_argument("--no-png", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "prerun":
        prerun()
    elif args.cmd == "draws":
        run_draws(args.n_workers, args.multistart_all, args.n_draws)
    else:
        summarize(png=not args.no_png)


if __name__ == "__main__":
    main()
