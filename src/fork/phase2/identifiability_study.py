"""Mercury: the synthetic identifiability study (pre-registered in
docs/phase2/identifiability_study.md; nothing here may differ from that plan).

    uv run python -m fork.phase2.identifiability_study prepare    # shapes, tuning, σ₁
    uv run python -m fork.phase2.identifiability_study run [--n-workers 12]
    uv run python -m fork.phase2.identifiability_study validate [--n-workers 12]
    uv run python -m fork.phase2.identifiability_study summarize
    uv run python -m fork.phase2.identifiability_study all [--n-workers 12]

No observation of Mercury is fitted. The contact table gives the epochs and the weights;
the data are the forward model's own longitudes under four truths, with drawn noise.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp

from fork import ledger
from fork.phase2 import identifiability as ID
from fork.phase2 import secular as S

ROOT = Path(__file__).resolve().parents[3]
CONTACTS = ROOT / "data" / "historical" / "mercury" / "newcomb_1882_contacts.csv"
REAL_CACHE = ROOT / "data" / "cache" / "phase2" / "identifiability_study"
CACHE = Path(os.environ.get("FORK_IDSTUDY_CACHE", REAL_CACHE))  # a scratch cache for smoke runs
_REAL = CACHE.resolve() == REAL_CACHE.resolve()
OUT_JSON = (ROOT / "docs" / "phase2" if _REAL else CACHE) / "identifiability_study.json"
COLUMNS_READ = ("year", "month", "day", "obs_h", "obs_m", "obs_s", "weight", "contact", "table")
ARCSEC = ID.ARCSEC
TARGET = ID.TARGET_ARCSEC_PER_CENTURY
TRUTHS = ("a", "b", "c", "d")
TRUTH_MODEL = {"a": "gr", "b": "ring", "c": "planet", "d": "j2"}
FIT_OF_TRUTH = {"a": "G", "b": "R", "c": "P", "d": "J"}
MODELS = ("N", "G", "J", "R", "P")
SEEDS = tuple(range(1, 9))
NOISE_FREE = 0
LEVELS = (1.0, 0.5, 2.0)  # × σ₁; index 0 is the calibrated level
PRECISION_1PN = 0.03
RING_RADIUS = PLANET_A = 0.2
RING_GRID = np.linspace(0.10, 0.30, 11)
PLANET_DISTANCES = (0.15, 0.20, 0.25)
N_PHASE = 48
LOG_MASS = (-9.0, -5.0)
LOG_J2 = (-7.0, -2.0)
R3_CENTRE = {"R": -7.0, "P": -7.0, "J": -3.5}  # J2: the amendment of 13 November 2026
R3_WIDTHS = (1.0, 2.0, 4.0, 6.0, 8.0)
NLIVE = 1000
DLOGZ = 0.01
N_EARTH = 2.0 * math.pi / 365.256363
V_TOL = {"rate": 0.05, "chi2_lin": 0.5, "frac": 0.95, "ln_z": 0.1, "closure": 0.5}
V2_DRAWS = 30
EL_STEPS = np.array([1e-10, 2e-7, 2e-7, 2e-7, 2e-7, 2e-7])
GRADES = ((1.0, "none"), (3.0, "weak"), (5.0, "strong"), (math.inf, "decisive"))


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with (CACHE / "progress.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def _wrap(d: np.ndarray) -> np.ndarray:
    return (np.asarray(d) + math.pi) % (2.0 * math.pi) - math.pi


def grade(x: float) -> str:
    return next(name for edge, name in GRADES if abs(x) < edge)


# ---------------------------------------------------------------------------- epochs


def julian_day_number(year: int, month: int, day: int) -> int:
    """JD at Greenwich noon of a Gregorian civil date (Fliegel and Van Flandern)."""
    a = (14 - month) // 12
    y, m = year + 4800 - a, month + 12 * a - 3
    return day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045


def contact_epochs(path: Path = CONTACTS) -> dict:
    """Epochs [JD, taken as TDB], weights and month (11 or 5) of the contacts. Only the
    columns of ``COLUMNS_READ`` are read: the date and time place the epoch, the weight
    sets the error; no other printed value is used."""
    months = {"Nov": 11, "May": 5}
    jd, w, mon = [], [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            row = {k: row[k] for k in COLUMNS_READ}
            m = months[row["month"]]
            hours = float(row["obs_h"]) + float(row["obs_m"]) / 60.0 + float(row["obs_s"]) / 3600.0
            # the astronomical day begins at noon of the civil date
            jd.append(julian_day_number(int(row["year"]), m, int(row["day"])) + hours / 24.0)
            w.append(float(row["weight"]))
            mon.append(m)
    return {"jd": np.array(jd), "weight": np.array(w), "month": np.array(mon)}


# ---------------------------------------------------------------------------- forward model


def base_problem():
    """(Mercury's six elements, context) at 1900.0 from the CI fixture."""
    theta, _, ctx = ID.fixture_problem()
    return np.asarray(theta, dtype=float)[:6], ctx


def lon(model: str, theta6: np.ndarray, jd: np.ndarray, ctx, extra=()) -> np.ndarray:
    return ID.synthetic(model, theta6, jd, ctx, tuple(extra))


def perihelion_rate(model: str, theta6, ctx, extra, span=(1677.0, 1881.0)) -> float:
    """Integrated rate of ϖ against the Newtonian run [″ per century] on an annual grid."""
    t = 2451545.0 + (np.arange(span[0], span[1] + 1.0) - 2000.0) * 365.25
    d = ID._elements(model, theta6, t, ctx, tuple(extra)) - ID._elements(
        "newtonian", theta6, t, ctx
    )
    return float(S.slope(t, np.unwrap(d[:, 3] + d[:, 4])) * ID.CENTURY_DAYS * ARCSEC)


def tuned_truths(theta6, ctx) -> dict:
    """Each truth's parameters with the integrated perihelion rate brought to the target by
    one linear rescaling of the analytic value (truth a keeps strength 1)."""
    m = ID.matched_parameters(theta6, ctx, ring_radius=RING_RADIUS, planet_a=PLANET_A)
    out = {}
    for t in TRUTHS:
        model = TRUTH_MODEL[t]
        extra = (1.0,) if t == "a" else m[model]
        if t != "a":
            r0 = perihelion_rate(model, theta6, ctx, extra)
            extra = (extra[0] * TARGET / r0, *extra[1:])
        out[t] = {
            "model": model,
            "extra": [float(v) for v in extra],
            "rate": perihelion_rate(model, theta6, ctx, extra),
        }
    return out


def element_jacobian(theta6, ctx, jd) -> np.ndarray:
    """∂λ/∂(six elements) at the epochs, central differences on the Newtonian model."""
    cols = []
    for i, h in enumerate(EL_STEPS):
        d = np.zeros(6)
        d[i] = h
        cols.append(
            _wrap(lon("newtonian", theta6 + d, jd, ctx) - lon("newtonian", theta6 - d, jd, ctx))
            / (2.0 * h)
        )
    return np.array(cols).T


def rates_vs_earth(theta6, ctx, jd, h: float = 0.01) -> np.ndarray:
    """ρ: the rate at which Mercury's heliocentric longitude gains on the Earth's [rad/day]."""
    d = _wrap(lon("newtonian", theta6, jd + h, ctx) - lon("newtonian", theta6, jd - h, ctx))
    return np.abs(d / (2.0 * h) - N_EARTH)


def shapes(theta6, ctx, jd) -> dict:
    """Each model's effect on the longitudes per unit of its strength, by full integrations
    at a strength that gives about the target advance (the linear regime the data live in)."""
    base = lon("newtonian", theta6, jd, ctx)

    def per_unit(model: str, extra: tuple) -> np.ndarray:
        return _wrap(lon(model, theta6, jd, ctx, extra) - base) / extra[0]

    m = ID.matched_parameters(theta6, ctx)
    out = {"G": per_unit("gr", (1.0,)), "J": per_unit("j2", m["j2"])}
    ring = []
    for r in RING_GRID:
        mr = ID.matched_parameters(theta6, ctx, ring_radius=float(r))["ring"]
        ring.append(per_unit("ring", mr))
    out["R"] = np.array(ring)  # (11, n)
    planet = []
    for a in PLANET_DISTANCES:
        mp = ID.matched_parameters(theta6, ctx, planet_a=a)["planet"][0]
        planet.append(
            [per_unit("planet", (mp, a, 2.0 * math.pi * k / N_PHASE)) for k in range(N_PHASE)]
        )
    out["P"] = np.array(planet)  # (3, 48, n)
    return out


def prepare() -> dict:
    """Everything fixed before data are drawn: epochs, Jacobian, shapes, tuned truths, σ₁."""
    path = CACHE / "prepared.npz"
    if path.exists():
        d = np.load(path, allow_pickle=False)
        out = {k: d[k] for k in d.files}
        out["truths"] = json.loads((CACHE / "truths.json").read_text())
        return out
    theta6, ctx = base_problem()
    ep = contact_epochs()
    jd = ep["jd"]
    _log(f"prepare: {len(jd)} epochs; Jacobian and shapes")
    jac = element_jacobian(theta6, ctx, jd)
    sh = shapes(theta6, ctx, jd)
    rho = rates_vs_earth(theta6, ctx, jd)  # rad / day
    unit = rho / 86400.0 / np.sqrt(ep["weight"])  # σ_i [rad] for σ₁ = 1 s
    a = np.c_[jac, sh["G"]]
    cov = np.linalg.inv(a.T @ (a / unit[:, None] ** 2))
    sigma1 = PRECISION_1PN / math.sqrt(cov[6, 6])
    _log(f"prepare: sigma_1 = {sigma1:.3f} s; tuning the truths")
    truths = tuned_truths(theta6, ctx)
    base = lon("newtonian", theta6, jd, ctx)
    y = np.array(
        [_wrap(lon(v["model"], theta6, jd, ctx, v["extra"]) - base) for v in truths.values()]
    )
    CACHE.mkdir(parents=True, exist_ok=True)
    (CACHE / "truths.json").write_text(json.dumps(truths, indent=1))
    out = {
        "jd": jd,
        "weight": ep["weight"],
        "month": ep["month"],
        "jac": jac,
        "G": sh["G"],
        "J": sh["J"],
        "R": sh["R"],
        "P": sh["P"],
        "rho": rho,
        "sigma_unit": unit,
        "sigma1_s": np.array(sigma1),
        "y_truth": y,
        "theta6": theta6,
    }
    np.savez(path, **out)
    out["truths"] = truths
    return out


# ---------------------------------------------------------------------------- data


def synthetic_data(y_truth: np.ndarray, sigma: np.ndarray, truth: int, seed: int) -> np.ndarray:
    """The truth's longitudes (as differences from the Newtonian reference) plus Gaussian
    noise of ``sigma``. Depends on nothing but its arguments: seed 0 is noise-free."""
    if seed == NOISE_FREE:
        return np.array(y_truth, dtype=float)
    z = np.random.default_rng([2026, 21, 1882, truth, seed]).standard_normal(len(y_truth))
    return y_truth + sigma * z


# ---------------------------------------------------------------------------- likelihood


class Fitter:
    """The marginal likelihood of each model for one data vector y at one error level: the
    six elements are marginalised analytically (flat priors, linear response), so
    −2 ln L = rᵀ P r with r = y − strength × shape and P the projector off the Jacobian."""

    def __init__(self, prep: dict, sigma: np.ndarray, y: np.ndarray) -> None:
        jac = prep["jac"]
        w = 1.0 / sigma**2
        jw = jac * w[:, None]
        self.jac, self.w = jac, w
        self.proj = np.diag(w) - jw @ np.linalg.solve(jac.T @ jw, jw.T)
        self.y = y
        self.py = self.proj @ y
        self.ypy = float(y @ self.py)
        self.shape_g, self.shape_j = prep["G"], prep["J"]
        self.ring = CubicSpline(RING_GRID, prep["R"], axis=0)
        self.planet_fft = np.fft.rfft(prep["P"], axis=1) / N_PHASE  # (3, 25, n)

    def planet_shape(self, ia: int, lam0: np.ndarray) -> np.ndarray:
        """Trigonometric interpolation in the planet's longitude at epoch; (..., n)."""
        c = self.planet_fft[ia]
        k = np.arange(c.shape[0])
        ph = np.exp(1j * np.multiply.outer(np.asarray(lam0, dtype=float), k))
        scale = np.full(c.shape[0], 2.0)
        scale[0] = 1.0
        scale[-1] = 1.0  # the Nyquist term of an even-length transform
        return np.real((ph * scale) @ c)

    def chi2(self, shape: np.ndarray, strength) -> np.ndarray:
        """rᵀ P r for r = y − strength × shape; shape (..., n), strength broadcast."""
        s = np.asarray(strength, dtype=float)
        ps = shape @ self.proj
        return self.ypy - 2.0 * s * (shape @ self.py) + s**2 * np.einsum("...i,...i", ps, shape)

    def lnl_fixed(self, model: str) -> float:
        if model == "N":
            return -0.5 * self.ypy
        return -0.5 * float(self.chi2(self.shape_g, 1.0))

    def lnl(self, model: str, x: np.ndarray, ia: int = 1) -> np.ndarray:
        """ln L at parameters x: J (log10 J2,), R (log10 m, radius), P (log10 m, λ₀)."""
        x = np.atleast_2d(x)
        s = 10.0 ** x[:, 0]
        if model == "J":
            return -0.5 * self.chi2(self.shape_j, s)
        shape = self.ring(x[:, 1]) if model == "R" else self.planet_shape(ia, x[:, 1])
        return -0.5 * self.chi2(shape, s)

    def best_elements(self, shape: np.ndarray, strength: float) -> np.ndarray:
        """The conditional best change of the six elements for residual y − strength × shape."""
        r = self.y - strength * shape
        jw = self.jac * self.w[:, None]
        return np.linalg.solve(self.jac.T @ jw, jw.T @ r)


def bounds(model: str, width: float | None = None) -> np.ndarray:
    """Prior box: the strength's log10 range (the plan's, or ``width`` decades about the R3
    centre) and the second parameter's range."""
    lo, hi = LOG_J2 if model == "J" else LOG_MASS
    if width is not None:
        lo, hi = R3_CENTRE[model] - width / 2.0, R3_CENTRE[model] + width / 2.0
    if model == "J":
        return np.array([[lo, hi]])
    second = (RING_GRID[0], RING_GRID[-1]) if model == "R" else (0.0, 2.0 * math.pi)
    return np.array([[lo, hi], second])


def ln_z_nested(f: Fitter, model: str, rng_seed: list[int], ia: int = 1):
    import dynesty

    b = bounds(model)
    kw = {"periodic": [1]} if model == "P" else {}
    sampler = dynesty.NestedSampler(
        lambda x: float(f.lnl(model, x, ia)[0]),
        lambda u: b[:, 0] + u * (b[:, 1] - b[:, 0]),
        len(b),
        nlive=NLIVE,
        bound="multi",
        sample="unif",
        rstate=np.random.default_rng(rng_seed),
        **kw,
    )
    sampler.run_nested(dlogz=DLOGZ, print_progress=False)
    r = sampler.results
    return float(r.logz[-1]), float(r.logzerr[-1]), r


def ln_z_quadrature(
    f: Fitter, model: str, width: float | None = None, ia: int = 1, n1: int = 2401, n2: int = 361
) -> float:
    """ln Z by the trapezoid rule on a grid over the prior box."""
    b = bounds(model, width)
    x1 = np.linspace(*b[0], n1)
    if model == "J":
        ll = f.lnl("J", x1[:, None])
        return float(logsumexp(ll, b=np.gradient(x1)) - math.log(b[0, 1] - b[0, 0]))
    x2 = np.linspace(*b[1], n2)
    rows = []
    for v in x2:
        ll = f.lnl(model, np.c_[x1, np.full(n1, v)], ia)
        rows.append(logsumexp(ll, b=_trap(x1)))
    return float(logsumexp(rows, b=_trap(x2)) - math.log(np.prod(b[:, 1] - b[:, 0])))


def _trap(x: np.ndarray) -> np.ndarray:
    w = np.full(len(x), x[1] - x[0])
    w[0] = w[-1] = 0.5 * (x[1] - x[0])
    return w


def planet_average(ln_z: list[float]) -> float:
    return float(logsumexp(ln_z) - math.log(len(ln_z)))


# ---------------------------------------------------------------------------- runs


def task_path(truth: int, seed: int, level: int) -> Path:
    return CACHE / "fits" / f"t{truth}_s{seed}_l{level}.json"


def run_task(task: tuple[int, int, int]) -> dict:
    """All five models on one data set: ln Z of each (nested sampling where a model has
    parameters). Cached per task, so a rerun resumes."""
    truth, seed, level = task
    path = task_path(truth, seed, level)
    if path.exists():
        return json.loads(path.read_text())
    prep = prepare()
    sigma = prep["sigma_unit"] * float(prep["sigma1_s"]) * LEVELS[level]
    y = synthetic_data(prep["y_truth"][truth - 1], sigma, truth, seed)
    f = Fitter(prep, sigma, y)
    out = {
        "truth": truth,
        "seed": seed,
        "level": level,
        "ln_z": {"N": f.lnl_fixed("N"), "G": f.lnl_fixed("G")},
        "ln_z_err": {},
    }
    for mi, model in ((2, "J"), (3, "R")):
        z, e, _ = ln_z_nested(f, model, [2026, 21, 1882, truth, seed, level, mi])
        out["ln_z"][model], out["ln_z_err"][model] = z, e
    zs = [
        ln_z_nested(f, "P", [2026, 21, 1882, truth, seed, level, 4, ia], ia)[:2]
        for ia in range(len(PLANET_DISTANCES))
    ]
    out["ln_z"]["P"] = planet_average([z for z, _ in zs])
    out["ln_z_planet_by_distance"] = [z for z, _ in zs]
    out["ln_z_err"]["P"] = max(e for _, e in zs)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(out))
    tmp.replace(path)
    return out


def all_tasks() -> list[tuple[int, int, int]]:
    return [
        (t, s, lv) for lv in range(len(LEVELS)) for t in range(1, 5) for s in (NOISE_FREE, *SEEDS)
    ]


def run(n_workers: int) -> None:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    prepare()
    todo = [t for t in all_tasks() if not task_path(*t).exists()]
    _log(f"run: {len(todo)} of {len(all_tasks())} data sets to fit, {n_workers} workers")
    with ProcessPoolExecutor(n_workers) as ex:
        futs = [ex.submit(run_task, t) for t in todo]
        for i, fu in enumerate(as_completed(futs)):
            fu.result()
            if (i + 1) % 12 == 0 or i + 1 == len(futs):
                _log(f"run: {i + 1}/{len(futs)} data sets")


# ---------------------------------------------------------------------------- validation


def _v2_one(task: tuple) -> list[float]:
    """χ² by a full integration at one posterior draw, and the linear χ²."""
    truth, model, ia, x = task
    prep = prepare()
    theta6, ctx = base_problem()
    sigma = prep["sigma_unit"] * float(prep["sigma1_s"])
    y = prep["y_truth"][truth - 1]
    f = Fitter(prep, sigma, y)
    s = 10.0 ** x[0]
    if model == "J":
        shape, name, extra = f.shape_j, "j2", (s,)
    elif model == "R":
        shape, name, extra = f.ring(x[1]), "ring", (s, x[1])
    else:
        shape, name = f.planet_shape(ia, x[1]), "planet"
        extra = (s, PLANET_DISTANCES[ia], x[1])
    d_el = f.best_elements(shape, s)
    full = _wrap(
        lon(name, theta6 + d_el, prep["jd"], ctx, extra) - lon("newtonian", theta6, prep["jd"], ctx)
    )
    return [float(np.sum(((y - full) / sigma) ** 2)), float(f.chi2(shape, s))]


def validate(n_workers: int) -> dict:
    """Gates V1 to V4 of the plan (V5 is a test in the test suite)."""
    from concurrent.futures import ProcessPoolExecutor

    path = CACHE / "validation.json"
    if path.exists():
        return json.loads(path.read_text())
    prep = prepare()
    truths = prep["truths"]
    v1 = {
        t: {"rate": v["rate"], "pass": abs(v["rate"] - TARGET) <= V_TOL["rate"]}
        for t, v in truths.items()
    }
    sigma = prep["sigma_unit"] * float(prep["sigma1_s"])
    v3, v4, tasks = [], {}, []
    for ti, t in enumerate(TRUTHS, start=1):
        f = Fitter(prep, sigma, prep["y_truth"][ti - 1])
        own = FIT_OF_TRUTH[t]
        if own == "G":
            v4[t] = -2.0 * f.lnl_fixed("G")
        for mi, model in ((2, "J"), (3, "R"), (4, "P")):
            for ia in range(len(PLANET_DISTANCES)) if model == "P" else (1,):
                seed = [2026, 21, 1882, ti, NOISE_FREE, 0, mi] + ([ia] if model == "P" else [])
                z, e, res = ln_z_nested(f, model, seed, ia)
                q = ln_z_quadrature(f, model, None, ia)
                v3.append(
                    {
                        "truth": t,
                        "model": model,
                        "distance": ia if model == "P" else None,
                        "nested": z,
                        "nested_err": e,
                        "quadrature": q,
                        "pass": abs(z - q) <= V_TOL["ln_z"],
                    }
                )
                if model == own and (model != "P" or ia == 1):
                    v4[t] = float(-2.0 * np.max(res.logl))
                eq = res.samples_equal(rstate=np.random.default_rng([*seed, 31]))
                n = V2_DRAWS if model != "P" else V2_DRAWS // len(PLANET_DISTANCES)
                for x in eq[:n]:
                    tasks.append((ti, model, ia, [float(v) for v in x]))
    _log(f"validate: V2 on {len(tasks)} draws by full integration, {n_workers} workers")
    with ProcessPoolExecutor(n_workers) as ex:
        pairs = np.array(list(ex.map(_v2_one, tasks, chunksize=4)))
    diff = np.abs(pairs[:, 0] - pairs[:, 1])
    v2 = []
    for ti, t in enumerate(TRUTHS, start=1):
        for model in ("J", "R", "P"):
            m = np.array([a[0] == ti and a[1] == model for a in tasks])
            frac = float(np.mean(diff[m] <= V_TOL["chi2_lin"]))
            v2.append(
                {
                    "truth": t,
                    "model": model,
                    "n": int(m.sum()),
                    "fraction_within": frac,
                    "max_abs_diff": float(diff[m].max()),
                    "pass": frac >= V_TOL["frac"],
                }
            )
    out = {
        "V1": {"truths": v1, "pass": all(v["pass"] for v in v1.values())},
        "V2": {"fits": v2, "pass": all(v["pass"] for v in v2)},
        "V3": {"fits": v3, "pass": all(v["pass"] for v in v3)},
        "V4": {"chi2_min": v4, "pass": all(v <= V_TOL["closure"] for v in v4.values())},
    }
    out["pass"] = all(out[k]["pass"] for k in ("V1", "V2", "V3", "V4"))
    path.write_text(json.dumps(out, indent=1))
    _log(f"validate: {'all gates pass' if out['pass'] else 'A GATE FAILS'}")
    return out


# ---------------------------------------------------------------------------- readings


def read_r1(fits: list[dict], truth: str) -> dict:
    """R1 on the nine data sets of one truth at one error level."""
    own = FIT_OF_TRUTH[truth]
    rows = []
    for f in fits:
        z = f["ln_z"]
        lead = min(z[own] - z[m] for m in MODELS if m != own)
        rows.append({"seed": f["seed"], "winner": max(z, key=z.get), "lead_of_truth": lead})
    noisy = [r for r in rows if r["seed"] != NOISE_FREE]
    free = next(r for r in rows if r["seed"] == NOISE_FREE)
    n_win = sum(r["winner"] == own for r in noisy)
    n_lead3 = sum(r["winner"] == own and r["lead_of_truth"] >= 3.0 for r in noisy)
    if n_lead3 >= 6 and free["winner"] == own and free["lead_of_truth"] >= 3.0:
        verdict = "recovered"
    elif n_win >= 6:
        verdict = "preferred, not decided"
    else:
        winners = [r["winner"] for r in noisy]
        verdict = "not recovered; most often " + max(set(winners), key=winners.count)
    return {
        "verdict": verdict,
        "wins_of_8": n_win,
        "wins_by_3_of_8": n_lead3,
        "noise_free": free,
        "noisy": noisy,
    }


def summarize() -> dict:
    prep = prepare()
    val = validate(1)
    out = {
        "plan": "docs/phase2/identifiability_study.md",
        "sigma1_s": float(prep["sigma1_s"]),
        "levels": list(LEVELS),
        "n_epochs": len(prep["jd"]),
        "truths": prep["truths"],
        "validation": val,
        "sigma_arcsec": {
            "november_median": float(
                np.median(
                    (prep["sigma_unit"] * float(prep["sigma1_s"]) * ARCSEC)[prep["month"] == 11]
                )
            ),
            "may_median": float(
                np.median(
                    (prep["sigma_unit"] * float(prep["sigma1_s"]) * ARCSEC)[prep["month"] == 5]
                )
            ),
        },
    }
    if not val["pass"]:
        out["readings"] = "not made: a validation gate failed"
        OUT_JSON.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
        return out
    fits = [json.loads(task_path(*t).read_text()) for t in all_tasks()]
    r1, r2, tables = {}, {}, {}
    for lv in range(len(LEVELS)):
        for ti, t in enumerate(TRUTHS, start=1):
            sel = [f for f in fits if f["truth"] == ti and f["level"] == lv]
            r1[f"{t}/level{lv}"] = read_r1(sel, t)
            noisy = [f for f in sel if f["seed"] != NOISE_FREE]
            free = next(f for f in sel if f["seed"] == NOISE_FREE)
            tables[f"{t}/level{lv}"] = {
                m: {
                    "noise_free": free["ln_z"][m] - free["ln_z"]["G"],
                    "median": float(np.median([f["ln_z"][m] - f["ln_z"]["G"] for f in noisy])),
                    "min": min(f["ln_z"][m] - f["ln_z"]["G"] for f in noisy),
                    "max": max(f["ln_z"][m] - f["ln_z"]["G"] for f in noisy),
                }
                for m in MODELS
            }
            if lv == 0:
                r2[t] = {
                    "ln_b_best_vs_newtonian": max(free["ln_z"][m] for m in MODELS[1:])
                    - free["ln_z"]["N"]
                }
    sigma = prep["sigma_unit"] * float(prep["sigma1_s"])
    r3 = {}
    for ti, t in enumerate(TRUTHS, start=1):
        f = Fitter(prep, sigma, prep["y_truth"][ti - 1])
        g = f.lnl_fixed("G")
        curves = {}
        for model in ("J", "R", "P"):
            vals = []
            for d in R3_WIDTHS:
                if model == "P":
                    z = planet_average(
                        [ln_z_quadrature(f, "P", d, ia) for ia in range(len(PLANET_DISTANCES))]
                    )
                else:
                    z = ln_z_quadrature(f, model, d)
                vals.append({"decades": d, "ln_b_vs_G": z - g, "grade": grade(z - g)})
            curves[model] = vals
        r3[t] = curves
    a_free = tables["a/level0"]
    out["readings"] = {
        "R1": {k: v for k, v in r1.items() if k.endswith("level0")},
        "R1_truth_a_G_against": {
            m: {"ln_b_G_over_m": -a_free[m]["noise_free"], "grade": grade(a_free[m]["noise_free"])}
            for m in ("J", "R", "P")
        },
        "R2": r2,
        "R3": r3,
        "R4": {k: v["verdict"] for k, v in r1.items()},
        "ln_b_vs_G": tables,
    }
    OUT_JSON.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    _log("summarize: written")
    return out


JOB = "phase2/identifiability_study"  # ledger rows (CLI steps; workers write none)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("prepare", "run", "validate", "summarize", "all"))
    ap.add_argument("--n-workers", type=int, default=12)
    a = ap.parse_args(argv)
    n = a.n_workers
    if a.cmd in ("prepare", "all"):
        with ledger.run(f"{JOB}/prepare", "prepare", n_workers=1):
            prepare()
    if a.cmd in ("validate", "all"):
        with ledger.run(f"{JOB}/validate", "gate", n_workers=n) as rec:
            ok = validate(n)["pass"]
            rec.done("gate passed" if ok else "gate failed")
    if a.cmd == "run" or (a.cmd == "all" and ok):
        with ledger.run(f"{JOB}/run", "nested", n_workers=n):
            run(n)
    if a.cmd in ("summarize", "all"):
        with ledger.run(f"{JOB}/summarize", "summarize", n_workers=1):
            summarize()


if __name__ == "__main__":
    main()
