"""Mercury: the identifiability study with the full forward model and contact times
(pre-registered in docs/phase2/identifiability_full.md; nothing here may differ).

    uv run python -m fork.phase2.identifiability_full prepare            # inputs, σ₁, V0, V1
    uv run python -m fork.phase2.identifiability_full stage1 [--n-workers 12]
    uv run python -m fork.phase2.identifiability_full stage2 [--n-workers 12]  # if Stage 1 passed
    uv run python -m fork.phase2.identifiability_full scan TASK          # one process per scan
    uv run python -m fork.phase2.identifiability_full plan | summarize

No observation of Mercury is read: the contact table gives the epochs, the kind of contact
and the weights; the data are the forward model's own output (plus noise from a fixed rng).

The pieces:

- **Forward model.** The inner-system integration of ``inner_system`` (Sun and eight
  planets, the validated step), returning Mercury and the Earth–Moon barycentre; each
  contact is the instant the apparent separation of Mercury's centre from the Sun's,
  seen from the barycentre, equals R☉ ∓ r☿ (interior, exterior), found by Newton's method
  on the offset from the printed epoch with both bodies carried by their Keplerian
  acceleration.
- **Likelihood call.** Levenberg–Marquardt on the six elements, every residual from a full
  integration, Broyden updates of the warm-start Jacobian, a fresh one when a step fails
  twice and at the end; Laplace's method at the optimum from the singular values of the
  whitened Jacobian.
- **Evidence.** Per second parameter (none for J; radius for R; longitude at epoch for P):
  a coarse strength scan, a fine scan of 15 points over ±6 widths widened to the levels it
  serves, six tail points; ln L splined on the fine scan and integrated on a sub-grid 20
  times finer, linearly between the other points.

ln Z is given up to a constant common to all models of one data set (the elements' flat
prior box and the data's normalisation), which cancels in every Bayes factor.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp

from fork import ledger
from fork.data.frames import elements_to_state
from fork.integrators import wh_extra as W
from fork.model import inner_system as IS
from fork.model.bodies import G, masses_msun
from fork.phase2 import identifiability as ID
from fork.phase2 import identifiability_study as IDS

ROOT = Path(__file__).resolve().parents[3]
PLAN = ROOT / "docs" / "phase2" / "identifiability_full.md"
PLAN_SHA256 = "0bf86af1ec4103046e56357e306a4950bc9a14a3c959b9843113337dedd03f7d"
TRUTHS_FILE = ROOT / "docs" / "phase2" / "identifiability_full_truths.json"
REAL_CACHE = ROOT / "data" / "cache" / "phase2" / "identifiability_full"
CACHE = Path(os.environ.get("FORK_IDFULL_CACHE", REAL_CACHE))  # scratch cache for smoke runs
OUT_JSON = ROOT / "docs" / "phase2" / "identifiability_full.json"
JOB = "phase2/identifiability_full"

# § 2: inputs
INPUTS = {
    "data/historical/mercury/newcomb_1882_contacts.csv":
        "70d4d0d1ea242c55c227d593b1cf9238593ce52c3cd39ae0ce9c1fdbd41712f9",
    "data/historical/mercury/newcomb_1882_misprints.csv":
        "79c4720ccae315cc4948c0a4d08c6928ec521ff95dc546fa939c9d6a9388fbce",
    "tests/fixtures/phase2/v4_reference.npz":
        "f8032fc2f622112303e9e3912fc32443b54e8e206251bb91babd7c3d29f24f77",
    "src/fork/model/inner_system.py":
        "a132da7aab691230d988ee3c2d42b2ec16655af7e65159856ff480425e6f512b",
    "src/fork/phase2/secular.py":
        "55d0eb506eaff1899731d119324474287bcb6091598e6dcd72e98213afab880b",
    "src/fork/phase2/identifiability.py":
        "14676bb8d02a945dae5ed66f8c75c13084857b4e899357cc0b5faf9fa3c4d3f0",
}  # fmt: skip
PLAN_TRUTHS = {
    "a": {"model": "gr", "extra": [1.0]},
    "b": {"model": "ring", "extra": [1.9811152404862693e-07, 0.2]},
    "c": {"model": "planet", "extra": [1.980728405006917e-07, 0.2, 0.0]},
    "d": {"model": "j2", "extra": [3.399891252673649e-04]},
}

# the earlier plan's choices, unchanged (§§ 1, 4, 5, 7 of identifiability_study.md)
TRUTHS, MODELS, SEEDS, NOISE_FREE = IDS.TRUTHS, IDS.MODELS, IDS.SEEDS, IDS.NOISE_FREE
LEVELS, PRECISION_1PN = IDS.LEVELS, IDS.PRECISION_1PN
FIT_OF_TRUTH, TRUTH_MODEL = IDS.FIT_OF_TRUTH, IDS.TRUTH_MODEL
RING_GRID, PLANET_DISTANCES, N_PHASE = IDS.RING_GRID, IDS.PLANET_DISTANCES, IDS.N_PHASE
LOG_MASS, LOG_J2, R3_CENTRE, R3_WIDTHS = IDS.LOG_MASS, IDS.LOG_J2, IDS.R3_CENTRE, IDS.R3_WIDTHS
EL_STEPS = IDS.EL_STEPS
STRENGTH_STEP_1PN = 0.01
MODEL_NAME = {"N": "newtonian", "G": "gr", "J": "j2", "R": "ring", "P": "planet"}

# § 3: the observation
ARCSEC = math.pi / (180.0 * 3600.0)
R_SUN = 959.63 * ARCSEC  # at 1 AU
R_MERCURY = 3.36 * ARCSEC
LEFT_OUT_S = 3600.0
V0_MAX_S = 1800.0
ROOT_TOL = 1e-12
MAX_OFFSET_DAYS = 1.0
NEWTON_ITER = 60

# § 4: the likelihood and the evidence
DCHI_TOL = 1e-3
MAX_IT = 20
COARSE_STEP, COARSE_HALF = 0.1, 4.0  # decades
FINE_N, FINE_HALF = 15, 6.0  # points, widths
SUB = 20
TAIL_DECADES = (1.0, 2.0, 4.0)
TAIL_LN = 25.0
HALF_BLOCKS = 4  # V3′'s halved planet grid (96 longitudes) runs as 4 processes of 24
MAX_WIDEN = 1500

# § 5: gates
V1_TOL = 0.05
V2_DRAWS, V2_NU, V2_SCALE = 800, 5.0, 1.5
V2_TOL, V2_ERR, V2_FRAC = 0.1, 0.03, 0.95
V2_DCHI = (1.0, 4.0, 9.0, 25.0)
V3_TOL = 0.05
V4_TOL = 0.5
V6_FRAC, V6_DCHI = 0.001, 25.0
V6_FAR = 200.0  # 50 L² at L = 2 (amendment 1, item 9)
V6_AMENDED = True  # amendment 1, item 9

# amendment 1, item 8: fresh Jacobians along the singular directions
JAC_MODE = os.environ.get("FORK_IDFULL_JAC", "whitened")
N_INT = [0]  # full integrations in this process


def _log(msg: str) -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    with (CACHE / "driver.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def _real() -> bool:
    return CACHE.resolve() == REAL_CACHE.resolve()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plan_frozen_part() -> str:
    text = PLAN.read_text()
    return text[: text.index("\n---\n") + 5]


def check_all_inputs() -> None:
    """The inputs of § 2, and the plan's text above its line as frozen (amendments below the
    line may grow)."""
    bad = [rel for rel, h in INPUTS.items() if _sha(ROOT / rel) != h]
    if hashlib.sha256(_plan_frozen_part().encode()).hexdigest() != PLAN_SHA256:
        bad.append("docs/phase2/identifiability_full.md, above its line")
    if bad:
        raise SystemExit(f"inputs differ from the plan: {bad}")


# ---------------------------------------------------------------------------- the table


def contact_rows(path: Path = IDS.CONTACTS) -> dict:
    """Epochs (JD, printed time taken as TDB, as the earlier plan), weights, and whether
    each contact is exterior (tables III and IV). Only the earlier plan's columns are read."""
    ep = IDS.contact_epochs(path)
    table, row = [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            table.append(r["table"])
            row.append(int(r["row"]))
    table = np.array(table)
    return {"jd": ep["jd"], "weight": ep["weight"], "exterior": np.isin(table, ("III", "IV")),
            "table": table, "row": np.array(row)}  # fmt: skip


# ---------------------------------------------------------------------------- forward model


_MASSES = masses_msun(IS.PLANETS)
MU_MERCURY = G * (1.0 + _MASSES[0])
MU_EMB = G * (1.0 + _MASSES[2])


def base_problem():
    return IDS.base_problem()


def states(model: str, theta6, extra, ctx, jd) -> tuple[np.ndarray, ...]:
    """Heliocentric (x, v) of Mercury and the Earth–Moon barycentre at ``jd``: the
    integration of ``inner_system.states`` with the same initial state and forces."""
    N_INT[0] += 1
    flags = ID.flags_for(model)
    th = np.concatenate([np.asarray(theta6, dtype=float), np.asarray(extra, dtype=float)])
    r, v = elements_to_state(*th[:6], MU_MERCURY)
    x0 = np.vstack([r, ctx.other_states[:, :3]])
    v0 = np.vstack([v, ctx.other_states[:, 3:]])
    masses = np.concatenate([[1.0], _MASSES])
    k, kw = 6, {}
    if flags.gr:
        kw["gr"], k = th[k], k + 1
    if flags.j2:
        kw["j2"], kw["spin_axis"], k = th[k], tuple(ctx.spin_axis), k + 1
    if flags.ring:
        kw["ring_mass"], kw["ring_radius_au"] = th[k], th[k + 1]
        kw["ring_normal"], k = tuple(ctx.plane_normal), k + 2
    if flags.planet:
        mp, ap, lam = th[k], th[k + 1], th[k + 2]
        u, w = IS._plane_basis(ctx.plane_normal)
        speed = math.sqrt(G * (1.0 + mp) / ap)
        x0 = np.vstack([x0, ap * (math.cos(lam) * u + math.sin(lam) * w)])
        v0 = np.vstack([v0, speed * (-math.sin(lam) * u + math.cos(lam) * w)])
        masses = np.append(masses, mp)
    t = np.asarray(jd, dtype=float) - ctx.epoch_jd_tdb
    x, vv = W.integrate_states(x0, v0, masses, t, ctx.dt, W.make_extra(**kw), ctx.order)
    return x[:, 0], vv[:, 0], x[:, 2], vv[:, 2]


def _carried(x0: np.ndarray, v0: np.ndarray, mu: float, t: np.ndarray) -> np.ndarray:
    a0 = -mu * x0 / np.linalg.norm(x0, axis=1)[:, None] ** 3
    return x0 + v0 * t[:, None] + 0.5 * a0 * t[:, None] ** 2


def _gap(xm, vm, xe, ve, exterior: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Separation of the centres minus the contact radius [rad], at offsets t [days]."""
    m, e = _carried(xm, vm, MU_MERCURY, t), _carried(xe, ve, MU_EMB, t)
    gm, gs = m - e, -e
    sep = np.arctan2(np.linalg.norm(np.cross(gm, gs), axis=1), np.einsum("ij,ij->i", gm, gs))
    rc = R_SUN / np.linalg.norm(gs, axis=1) + np.where(exterior, 1.0, -1.0) * R_MERCURY / (
        np.linalg.norm(gm, axis=1)
    )
    return sep - rc


def contact_offsets(xm, vm, xe, ve, exterior: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Offsets [s] of each contact from its printed epoch, and the root's residual [rad];
    NaN where Newton's method finds no contact within a day."""
    t = np.zeros(len(exterior))
    h = 1e-6
    ok = np.ones(len(t), bool)
    for _ in range(NEWTON_ITER):
        f = _gap(xm, vm, xe, ve, exterior, t)
        d = (_gap(xm, vm, xe, ve, exterior, t + h) - f) / h
        step = np.where(ok & (d != 0.0), f / np.where(d == 0.0, 1.0, d), 0.0)
        t = t - step
        ok &= np.isfinite(t) & (np.abs(t) <= MAX_OFFSET_DAYS)
        t = np.where(ok, t, 0.0)
        if np.all(np.abs(step[ok]) < 1e-9):
            break
    res = np.abs(_gap(xm, vm, xe, ve, exterior, t))
    ok &= res < ROOT_TOL
    return np.where(ok, t * 86400.0, np.nan), np.where(ok, res, np.nan)


class NoContact(Exception):
    """The model gives no contact within a day of a printed epoch: its likelihood is zero."""


@dataclass(frozen=True)
class Observation:
    jd: np.ndarray
    weight: np.ndarray
    exterior: np.ndarray
    ctx: object


def observe(model: str, theta6, extra, obs: Observation) -> np.ndarray:
    """The model's contact times [s from the printed epochs] at the kept rows."""
    off, _ = contact_offsets(*states(model, theta6, extra, obs.ctx, obs.jd), obs.exterior)
    if not np.all(np.isfinite(off)):
        raise NoContact
    return off


# ---------------------------------------------------------------------------- one call


@dataclass
class Call:
    chi2: float
    ln_s: float  # Σ ln s_k of the whitened Jacobian at the optimum
    theta6: np.ndarray
    jac: np.ndarray | None
    iterations: int
    n_int: int
    converged: bool
    no_contact: bool = False

    def ln_l(self, level: float = 1.0) -> float:
        """Laplace's marginal over the six elements at error level ``level`` (× σ₁), up to
        the constant common to every model: −χ²/(2 L²) − Σ ln s_k + 6 ln L."""
        if self.no_contact or not math.isfinite(self.chi2):
            return -math.inf
        return -0.5 * self.chi2 / level**2 - self.ln_s + 6.0 * math.log(level)

    def row(self) -> dict:
        return {"chi2": self.chi2, "ln_s": self.ln_s, "theta6": self.theta6.tolist(),
                "iterations": self.iterations, "n_int": self.n_int,
                "converged": self.converged, "no_contact": self.no_contact}  # fmt: skip


def _jacobian(resid, theta: np.ndarray, r: np.ndarray, steps: np.ndarray) -> np.ndarray:
    cols = []
    for k, h in enumerate(steps):
        d = np.zeros(len(theta))
        d[k] = h
        try:
            cols.append((resid(theta + d) - r) / h)
        except NoContact:
            cols.append((r - resid(theta - d)) / h)
    return np.array(cols).T


def _whitened_jacobian(resid, theta: np.ndarray, r: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Forward differences along the singular directions of ``basis``, each step one whitened
    unit (θ + v_k / s_k), returned in the elements' coordinates."""
    _, s, vt = np.linalg.svd(basis, full_matrices=False)
    cols = []
    for k in range(len(s)):
        d = vt[k] / s[k]
        try:
            cols.append(resid(theta + d) - r)
        except NoContact:
            cols.append(r - resid(theta - d))
    return np.array(cols).T @ (s[:, None] * vt)


def lm_fit(resid, start: np.ndarray, steps: np.ndarray, jac0: np.ndarray | None = None,
           tol: float = DCHI_TOL, max_it: int = MAX_IT,
           jac_mode: str = "plan") -> Call:  # fmt: skip
    """Levenberg–Marquardt on whitened residuals ``resid(θ)``. With ``jac0`` (the previous
    point's Jacobian) it is updated by Broyden's rank-one rule after each step; a fresh
    Jacobian replaces it when a step fails twice. Converged when the last step lowered χ² by
    less than ``tol`` (or none was taken) and the Gauss–Newton step from a fresh Jacobian
    predicts a further drop below ``tol``; that Jacobian gives Laplace's term. A fresh
    Jacobian is by forward differences with ``steps`` (``jac_mode="plan"``), or, given one to
    take directions from, along its singular directions (``"whitened"``)."""
    n0 = N_INT[0]
    th = np.array(start, dtype=float)
    try:
        r = resid(th)
    except NoContact:
        return Call(math.inf, math.nan, th, None, 0, N_INT[0] - n0, False, True)

    def fresh(theta, res, basis):
        if jac_mode == "whitened" and basis is not None:
            return _whitened_jacobian(resid, theta, res, basis)
        return _jacobian(resid, theta, res, steps)

    def predicted(jac, res) -> float:
        u = np.linalg.svd(jac, full_matrices=False)[0]
        return float(np.sum((u.T @ res) ** 2))

    chi = float(r @ r)
    is_fresh = jac0 is None
    jac = fresh(th, r, None) if is_fresh else np.array(jac0, dtype=float)
    lam = 1e-3 if is_fresh else 1e-6
    fails, it, converged, last_drop = 0, 0, False, 0.0
    while True:
        if is_fresh and last_drop < tol and predicted(jac, r) < tol:
            converged = True
            break
        if it >= max_it:
            break
        it += 1
        a = jac.T @ jac
        step = np.linalg.solve(a + lam * np.diag(np.diag(a)) + 1e-300 * np.eye(len(th)),
                               -jac.T @ r)  # fmt: skip
        try:
            r2 = resid(th + step)
            chi2 = float(r2 @ r2)
        except NoContact:
            r2, chi2 = None, math.inf
        if r2 is not None and step @ step > 0.0:
            jac = jac + np.outer(r2 - r - jac @ step, step) / (step @ step)
        is_fresh = False
        if chi2 < chi:
            last_drop = chi - chi2
            th, r, chi = th + step, r2, chi2
            lam, fails = max(lam / 10.0, 1e-9), 0
            if last_drop < tol:
                jac, is_fresh = fresh(th, r, jac), True
        else:
            lam *= 10.0
            fails += 1
            if fails >= 2:
                jac, is_fresh, fails, last_drop = fresh(th, r, jac), True, 0, 0.0
    final = jac if is_fresh else fresh(th, r, jac)
    s = np.linalg.svd(final, compute_uv=False)
    return Call(chi, float(np.sum(np.log(s))), th, final, it, N_INT[0] - n0, converged)


class Data:
    """One data set: contact times y [s] with errors σ_i = σ₁ / √w_i at level 1."""

    def __init__(self, obs: Observation, y: np.ndarray, sigma: np.ndarray) -> None:
        self.obs, self.y, self.sigma = obs, np.asarray(y, float), np.asarray(sigma, float)

    def resid_fn(self, model: str, extra):
        def resid(theta6: np.ndarray) -> np.ndarray:
            return (observe(MODEL_NAME[model], theta6, extra, self.obs) - self.y) / self.sigma

        return resid

    def call(self, model: str, extra, start, jac0=None) -> Call:
        return lm_fit(self.resid_fn(model, extra), start, EL_STEPS, jac0, jac_mode=JAC_MODE)


# ---------------------------------------------------------------------------- prepared


def truths() -> dict:
    t = json.loads(TRUTHS_FILE.read_text())
    if {k: {"model": v["model"], "extra": v["extra"]} for k, v in t.items()} != PLAN_TRUTHS:
        raise SystemExit("the truths file differs from the plan's table")
    return t


def prepare() -> dict:
    """Everything fixed before any data set is drawn: the kept rows (§ 3's rule, on truth a's
    noise-free model), σ₁ (§ 3), V0 and V1. Cached."""
    path = CACHE / "prepared.json"
    if path.exists():
        return json.loads(path.read_text())
    if _real():
        check_all_inputs()
    tr = truths()
    theta6, ctx = base_problem()
    rows = contact_rows()
    a = tr["a"]
    off, res = contact_offsets(*states(a["model"], theta6, a["extra"], ctx, rows["jd"]),
                               rows["exterior"])  # fmt: skip
    keep = np.isfinite(off) & (np.abs(off) <= LEFT_OUT_S)
    left = [f"{t} {r}" for t, r, k in zip(rows["table"], rows["row"], keep, strict=True) if not k]
    v0 = {"left_out": left, "max_kept_offset_s": float(np.max(np.abs(off[keep]))),
          "max_root_residual_rad": float(np.max(res[keep])),
          "offsets_s": [None if not math.isfinite(x) else float(x) for x in off]}  # fmt: skip
    v0["pass"] = bool(left == ["III 12"] and v0["max_kept_offset_s"] <= V0_MAX_S
                      and v0["max_root_residual_rad"] < ROOT_TOL)  # fmt: skip
    obs = Observation(rows["jd"][keep], rows["weight"][keep], rows["exterior"][keep], ctx)
    unit = 1.0 / np.sqrt(obs.weight)
    base = observe("gr", theta6, [1.0], obs)
    cols = []
    for k, h in enumerate(EL_STEPS):
        d = np.zeros(6)
        d[k] = h
        cols.append((observe("gr", theta6 + d, [1.0], obs) - base) / h)
    cols.append((observe("gr", theta6, [1.0 + STRENGTH_STEP_1PN], obs) - base) / STRENGTH_STEP_1PN)
    amat = np.array(cols).T / unit[:, None]
    _, s, vt = np.linalg.svd(amat, full_matrices=False)
    cov = (vt.T / s**2) @ vt
    sigma1 = PRECISION_1PN / math.sqrt(cov[6, 6])
    v1 = {}
    for t, v in tr.items():
        rate = IDS.perihelion_rate(v["model"], theta6, ctx, v["extra"])
        v1[t] = {"rate": rate, "pass": abs(rate - ID.TARGET_ARCSEC_PER_CENTURY) <= V1_TOL}
    out = {
        "keep": keep.tolist(),
        "sigma1_s": sigma1,
        "n_kept": int(keep.sum()),
        "V0": v0,
        "V1": {"truths": v1, "pass": all(x["pass"] for x in v1.values())},
    }
    y = {}
    for t, v in tr.items():
        y[t] = observe(v["model"], theta6, v["extra"], obs).tolist()
    out["y_truth"] = y
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    _log(
        f"prepare: {out['n_kept']} contacts kept, left out {left}; sigma_1 = {sigma1:.3f} s; "
        f"V0 {'pass' if v0['pass'] else 'FAIL'}, V1 {'pass' if out['V1']['pass'] else 'FAIL'}"
    )
    return out


def observation(prep: dict) -> Observation:
    _, ctx = base_problem()
    rows = contact_rows()
    keep = np.array(prep["keep"])
    return Observation(rows["jd"][keep], rows["weight"][keep], rows["exterior"][keep], ctx)


def synthetic_data(y_truth, sigma, truth: int, seed: int) -> np.ndarray:
    """The truth's contact times plus Gaussian noise of ``sigma``, rng [2026, 21, 1882, truth,
    seed]. Depends on nothing but its arguments: seed 0 is noise-free."""
    y = np.asarray(y_truth, dtype=float)
    if seed == NOISE_FREE:
        return y.copy()
    z = np.random.default_rng([2026, 21, 1882, truth, seed]).standard_normal(len(y))
    return y + np.asarray(sigma) * z


def dataset(prep: dict, truth: int, seed: int, level: float = 1.0) -> Data:
    obs = observation(prep)
    sigma = prep["sigma1_s"] / np.sqrt(obs.weight)
    y = synthetic_data(prep["y_truth"][TRUTHS[truth - 1]], level * sigma, truth, seed)
    return Data(obs, y, sigma)


# ---------------------------------------------------------------------------- scans


@dataclass(frozen=True)
class Task:
    """One scan: a data set (truth, seed, level index; seed 0 serves every level) and a model;
    for P one distance and a block of longitudes; ``half`` for V3′'s halved grids."""

    truth: int
    seed: int
    level: int
    model: str  # "NG", "J", "R" or "P"
    distance: int = -1
    block: int = -1
    half: bool = False

    @property
    def id(self) -> str:
        s = f"t{self.truth}_s{self.seed}_l{self.level}_{self.model}"
        if self.model == "P":
            s += f"_d{self.distance}_b{self.block}"
        return s + ("_half" if self.half else "")

    @staticmethod
    def from_id(text: str) -> Task:
        parts = text.split("_")
        t, s, lv, m = int(parts[0][1:]), int(parts[1][1:]), int(parts[2][1:]), parts[3]
        d = b = -1
        if m == "P":
            d, b = int(parts[4][1:]), int(parts[5][1:])
        return Task(t, s, lv, m, d, b, parts[-1] == "half")

    @property
    def dir(self) -> Path:
        return CACHE / "scans" / self.id

    @property
    def done_path(self) -> Path:
        return self.dir / "result.json"

    def levels(self) -> tuple[float, ...]:
        return LEVELS if self.seed == NOISE_FREE else (LEVELS[self.level],)

    def second_grid(self) -> np.ndarray:
        """Radii for R; this block's longitudes for P (halved spacing: twice as many)."""
        if self.model == "R":
            return np.linspace(RING_GRID[0], RING_GRID[-1], 21 if self.half else 11)
        if self.model == "P":
            n = 2 * N_PHASE if self.half else N_PHASE
            per = n // HALF_BLOCKS if self.half else n
            lam = 2.0 * math.pi * np.arange(n) / n
            return lam[self.block * per : (self.block + 1) * per]
        return np.array([0.0])

    def extra(self, x: float, p2: float) -> tuple[float, ...]:
        s = 10.0**x
        if self.model == "J":
            return (s,)
        if self.model == "R":
            return (s, float(p2))
        return (s, PLANET_DISTANCES[self.distance], float(p2))


def stage1_tasks() -> list[Task]:
    out = []
    for t in range(1, 5):
        out += [Task(t, 0, 0, "NG"), Task(t, 0, 0, "J"), Task(t, 0, 0, "R")]
        out += [Task(t, 0, 0, "P", d, 0) for d in range(3)]
    out += [Task(2, 0, 0, "R", half=True), Task(4, 0, 0, "J", half=True)]
    out += [Task(3, 0, 0, "P", 1, b, half=True) for b in range(HALF_BLOCKS)]
    return out


def stage2_tasks() -> list[Task]:
    out = []
    for s in SEEDS:
        for t in range(1, 5):
            out += [Task(t, s, 0, "NG"), Task(t, s, 0, "J"), Task(t, s, 0, "R")]
            out += [Task(t, s, 0, "P", d, 0) for d in range(3)]
    return out


class Points:
    """The calls of one scan, one JSON line per finished grid point (so a scan resumes)."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.calls: dict[tuple[float, float], dict] = {}
        if path.exists():
            for line in path.read_text().splitlines():
                if line.strip():
                    d = json.loads(line)
                    self.calls[self.key(d["p2"], d["x"])] = d

    @staticmethod
    def key(p2: float, x: float) -> tuple[float, float]:
        return (round(float(p2), 10), round(float(x), 10))

    def get(self, p2, x):
        return self.calls.get(self.key(p2, x))

    def add(self, p2, x, call: Call, kind: str) -> dict:
        d = {"p2": float(p2), "x": float(x), "kind": kind, **call.row()}
        self.calls[self.key(p2, x)] = d
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a") as fh:
            fh.write(json.dumps(d) + "\n")
        return d


class Scanner:
    """The strength scans of one task (§ 4.2)."""

    def __init__(self, task: Task, data: Data, start: np.ndarray) -> None:
        self.task, self.data, self.start = task, data, np.asarray(start, float)
        self.points = Points(task.dir / "points.jsonl")
        self.warm: tuple[np.ndarray, np.ndarray | None] = (self.start, None)
        self.jacs: dict[tuple[float, float], np.ndarray] = {}
        self.h_scale = 0.5 if task.half else 1.0
        self.max_level = max(task.levels())

    def call(self, p2: float, x: float, kind: str, warm=None) -> dict:
        got = self.points.get(p2, x)
        if got is not None:
            if not got["no_contact"]:
                self.warm = self._warm_of(got)
            return got
        th0, j0 = warm if warm is not None else self.warm
        c = self.data.call(self.task.model, self.task.extra(x, p2), th0, j0)
        if not c.no_contact:
            self.jacs[Points.key(p2, x)] = c.jac
            self.warm = (c.theta6, c.jac)
        else:
            self.warm = (th0, None)
        return self.points.add(p2, x, c, kind)

    def _warm_of(self, d: dict) -> tuple[np.ndarray, np.ndarray | None]:
        return (np.array(d["theta6"]), self.jacs.get(Points.key(d["p2"], d["x"])))

    def coarse(self, p2: float) -> list[dict]:
        c = R3_CENTRE["J" if self.task.model == "J" else self.task.model]
        step = COARSE_STEP * self.h_scale
        xs = c + np.arange(-COARSE_HALF, COARSE_HALF + 1e-9, step)
        self.warm = (self.start, None)
        return [self.call(p2, x, "coarse") for x in xs]

    @staticmethod
    def _parabola(pts: list[dict]) -> tuple[float, float]:
        """Best x and width from the lowest point and its neighbours: χ² ≈ χ²₀ + (x−x₀)²/w²."""
        ok = [p for p in sorted(pts, key=lambda p: p["x"]) if math.isfinite(p["chi2"])]
        i = int(np.argmin([p["chi2"] for p in ok]))
        i = min(max(i, 1), len(ok) - 2)
        x = np.array([ok[j]["x"] for j in (i - 1, i, i + 1)])
        y = np.array([ok[j]["chi2"] for j in (i - 1, i, i + 1)])
        a, b, _ = np.polyfit(x, y, 2)
        if a <= 0:
            return float(x[1]), float(abs(x[2] - x[0]) / 6.0)
        return float(-b / (2 * a)), float(1.0 / math.sqrt(a))

    def fine(self, p2: float, x0: float, w: float, warm) -> tuple[list[dict], float, float, float]:
        """Fine points about x0 at spacing 12 w / 14 (halved for V3′), recentred and refined
        until they cover ±6 widths of the best value, then widened until the integrand at
        both ends is below e⁻²⁵ of its peak at every level the scan serves."""
        for _ in range(6):
            h = 2.0 * FINE_HALF * w / (FINE_N - 1) * self.h_scale
            pts = self._grid(p2, x0, h, round(FINE_HALF * w / h), warm)
            xb, wb = self._best(pts)
            if wb < 0.7 * w or wb > 1.5 * w or abs(xb - x0) > w:
                x0, w = xb, wb
                warm = self._warm_of(min(pts, key=lambda p: abs(p["x"] - xb)))
                continue
            break
        need = 2.0 * TAIL_LN * self.max_level**2
        lo_n = hi_n = round(FINE_HALF * w / h)
        for _ in range(4000):
            fin = {round(p["x"], 10): p for p in pts}
            xs = sorted(fin)
            cmin = min(p["chi2"] for p in pts)
            grow = False
            if fin[xs[0]]["chi2"] - cmin < need:
                lo_n += 1
                grow = True
            if fin[xs[-1]]["chi2"] - cmin < need:
                hi_n += 1
                grow = True
            if not grow:
                break
            pts = self._grid(p2, x0, h, (lo_n, hi_n), None)
        return pts, x0, w, h, (lo_n, hi_n)

    def _grid(self, p2, x0, h, n, warm) -> list[dict]:
        lo_n, hi_n = (n, n) if isinstance(n, int) else n
        ks = list(range(0, hi_n + 1)) + list(range(-1, -lo_n - 1, -1))
        out, centre = [], None
        for k in ks:
            if k == -1 and centre is not None:
                self.warm = self._warm_of(centre)
            d = self.call(p2, x0 + k * h, "fine", warm if k == 0 else None)
            if k == 0:
                centre = d
            out.append(d)
        return out

    @staticmethod
    def _best(pts: list[dict], with_min: bool = False):
        ok = sorted((p for p in pts if math.isfinite(p["chi2"])), key=lambda p: p["x"])
        x = np.array([p["x"] for p in ok])
        cs = CubicSpline(x, [p["chi2"] for p in ok])
        sub = np.linspace(x[0], x[-1], (len(x) - 1) * SUB + 1)
        xb = float(sub[int(np.argmin(cs(sub)))])
        curv = float(cs(xb, 2))
        wb = math.sqrt(2.0 / curv) if curv > 0 else float(x[-1] - x[0]) / 12.0
        return (xb, wb, float(cs(xb))) if with_min else (xb, wb)

    def tails(self, p2: float, xb: float, best: dict) -> list[dict]:
        out = []
        for dx in TAIL_DECADES:
            for sgn in (1.0, -1.0):
                self.warm = self._warm_of(best)
                out.append(self.call(p2, xb + sgn * dx, "tail"))
        return out

    def widen_to(self, p2: float, x0: float, h: float, n: tuple[int, int], target: float,
                 best: dict) -> tuple[int, int]:  # fmt: skip
        """Extend the fine points at spacing h toward ``target`` (a tail point whose integrand
        is not negligible), until they reach it or MAX_WIDEN points are added."""
        lo_n, hi_n = n
        side = 1 if target > x0 else -1
        self.warm = self._warm_of(best)
        for _ in range(MAX_WIDEN):
            if (x0 + hi_n * h >= target) if side > 0 else (x0 - lo_n * h <= target):
                break
            if side > 0:
                hi_n += 1
                self.call(p2, x0 + hi_n * h, "fine")
            else:
                lo_n += 1
                self.call(p2, x0 - lo_n * h, "fine")
        return lo_n, hi_n

    def scan_second(self, p2: float, guess: tuple[float, float] | None, warm) -> dict:
        if guess is None:
            pts = self.coarse(p2)
            x0, w = self._parabola(pts)
            warm = self._warm_of(min((p for p in pts if math.isfinite(p["chi2"])),
                                     key=lambda p: p["chi2"]))  # fmt: skip
        else:
            x0, w = guess
        fine, x0, w, h, n = self.fine(p2, x0, w, warm)
        xb, wb = self._best(fine)
        best = min(fine, key=lambda p: p["chi2"])
        widened = []
        for t in self.tails(p2, xb, best):
            if any(_ln_l(t, lv) - _ln_l(best, lv) >= -TAIL_LN for lv in self.task.levels()):
                n = self.widen_to(p2, x0, h, n, t["x"], best)
                widened.append(t["x"])
        if widened:
            fine = [
                p
                for p in self.points.calls.values()
                if p["kind"] == "fine" and Points.key(p["p2"], 0)[0] == Points.key(p2, 0)[0]
            ]
            best = min(fine, key=lambda p: p["chi2"])
        xb, wb, spline_min = self._best(fine, with_min=True)
        return {"p2": p2, "x_best": xb, "width": wb, "h": h, "chi2_min": best["chi2"],
                "chi2_spline_min": spline_min,
                "best_point": best, "widened_to_tails": widened}  # fmt: skip

    def run(self) -> list[dict]:
        out, guess, warm = [], None, None
        for p2 in self.task.second_grid():
            r = self.scan_second(float(p2), guess, warm)
            out.append(r)
            guess = (r["x_best"], r["width"])
            warm = self._warm_of(r["best_point"])
        return out


# ---------------------------------------------------------------------------- quadrature


def ln_z_strength(pts: list[dict], level: float, box: tuple[float, float],
                  h_fine: float) -> float:  # fmt: skip
    """ln of (1/width of box) ∫ L dx over the box at one value of the second parameter: ln L
    splined on the fine points and integrated on a sub-grid 20 times finer; outside the fine
    points, ln L interpolated linearly between the coarse and tail points on the same step."""
    allp = sorted(pts, key=lambda p: p["x"])
    x = np.array([p["x"] for p in allp])
    ll = np.array([_ln_l(p, level) for p in allp])
    fin = np.isfinite(ll)
    fine = sorted((p for p in allp if p["kind"] == "fine"), key=lambda p: p["x"])
    xf = np.array([p["x"] for p in fine])
    llf = np.array([_ln_l(p, level) for p in fine])
    sub_f = np.linspace(xf[0], xf[-1], (len(xf) - 1) * SUB + 1)
    val_f = CubicSpline(xf, llf)(sub_f) if np.all(np.isfinite(llf)) else np.interp(sub_f, xf, llf)
    lo, hi = box
    step = h_fine / SUB
    left = np.arange(xf[0] - step, lo - step, -step)[::-1] if lo < xf[0] else np.empty(0)
    right = np.arange(xf[-1] + step, hi + step, step) if hi > xf[-1] else np.empty(0)

    def lin(xx: np.ndarray) -> np.ndarray:
        if xx.size == 0 or not fin.any():
            return np.full(xx.size, -np.inf)
        return np.interp(xx, x[fin], ll[fin], left=-np.inf, right=-np.inf)

    xs = np.concatenate([left, sub_f, right])
    vals = np.concatenate([lin(left), val_f, lin(right)])
    m = (xs >= lo) & (xs <= hi)
    xs, vals = xs[m], vals[m]
    if xs.size < 2:
        return -math.inf
    return float(logsumexp(vals, b=trapezoid_weights(xs)) - math.log(hi - lo))


def trapezoid_weights(x: np.ndarray) -> np.ndarray:
    """Trapezoid weights on an increasing, possibly uneven grid (half intervals at the ends)."""
    d = np.diff(x)
    w = np.zeros(x.size)
    w[:-1] += 0.5 * d
    w[1:] += 0.5 * d
    return w


def _ln_l(p: dict, level: float) -> float:
    if p["no_contact"] or not math.isfinite(p["chi2"]):
        return -math.inf
    return -0.5 * p["chi2"] / level**2 - p["ln_s"] + 6.0 * math.log(level)


def box(model: str, width: float | None = None) -> tuple[float, float]:
    if width is None:
        return LOG_J2 if model == "J" else LOG_MASS
    c = R3_CENTRE[model]
    return (c - width / 2.0, c + width / 2.0)


def tail_check(scan: dict, pts: list[dict], levels) -> bool:
    """Every tail point's integrand below e⁻²⁵ of the peak, at every level."""
    tails = [p for p in pts if p["kind"] == "tail"]
    peak = scan["best_point"]
    return all(_ln_l(p, lv) - _ln_l(peak, lv) < -TAIL_LN for p in tails for lv in levels)


# ---------------------------------------------------------------------------- a scan task


@ledger.recorded
def run_task(task: Task) -> dict:
    if task.done_path.exists():
        return json.loads(task.done_path.read_text())
    prep = prepare()
    data = dataset(prep, task.truth, task.seed, LEVELS[task.level])
    theta6, _ = base_problem()
    N_INT[0] = 0
    rec = ledger.start(f"{JOB}/{task.id}", "scan", n_workers=1, scratch=not _real(),
                       n_likelihood_calls=None, sim_years=None)  # fmt: skip
    rec.set(n_integrations=lambda: N_INT[0])
    t0 = time.time()
    if task.model == "NG":
        pts = Points(task.dir / "points.jsonl")
        res = {}
        for m in ("N", "G"):
            d = pts.get(MODELS.index(m), 0.0)
            if d is None:
                extra = () if m == "N" else (1.0,)
                d = pts.add(MODELS.index(m), 0.0, data.call(m, extra, theta6), "fixed")
            res[m] = d
        out = {
            "task": task.id,
            "calls": res,
            "ln_z": {m: {str(lv): _ln_l(res[m], lv) for lv in task.levels()} for m in res},
        }
    else:
        sc = Scanner(task, data, theta6)
        scans = sc.run()
        pts_by = {}
        for p in sc.points.calls.values():
            pts_by.setdefault(Points.key(p["p2"], 0)[0], []).append(p)
        rows = []
        model = task.model
        for s in scans:
            pts = pts_by[Points.key(s["p2"], 0)[0]]
            z = {str(lv): {"prior": ln_z_strength(pts, lv, box(model), s["h"]),
                           **{f"r3_{d:g}": ln_z_strength(pts, lv, box(model, d), s["h"])
                              for d in R3_WIDTHS}}
                 for lv in task.levels()}  # fmt: skip
            rows.append({**{k: v for k, v in s.items() if k != "best_point"},
                         "best_theta6": s["best_point"]["theta6"], "ln_z": z,
                         "tails_ok": tail_check(s, pts, task.levels()),
                         "n_points": len(pts)})  # fmt: skip
        out = {"task": task.id, "scans": rows}
    calls = list(Points(task.dir / "points.jsonl").calls.values())
    out["calls_summary"] = {
        "n": len(calls), "not_converged": sum(not c["converged"] and not c["no_contact"]
                                              for c in calls),
        "no_contact": sum(c["no_contact"] for c in calls),
        "max_iterations": max(c["iterations"] for c in calls),
        "integrations": sum(c["n_int"] for c in calls)}  # fmt: skip
    out["wall_s"] = time.time() - t0
    tmp = task.done_path.with_suffix(".tmp")
    tmp.write_text(json.dumps(out, indent=1, default=float))
    tmp.replace(task.done_path)
    rec.done(n_likelihood_calls=len(calls))
    return out


# ---------------------------------------------------------------------------- combining


def ln_z_task_second(res: dict, level: float, key: str = "prior") -> np.ndarray:
    return np.array([r["ln_z"][str(level)][key] for r in res["scans"]])


def ln_z_model(results: dict[str, dict], truth: int, seed: int, level_idx: int, model: str,
               level: float, key: str = "prior", half: bool = False,
               distances=None) -> float:  # fmt: skip
    """A model's ln Z on one data set, from its scan tasks: J alone; R by the trapezoid over
    radii; P by the mean over longitudes, then over distances with equal weight."""
    if model in ("N", "G"):
        return results[Task(truth, seed, level_idx, "NG").id]["ln_z"][model][str(level)]
    if model == "J":
        return float(ln_z_task_second(results[Task(truth, seed, level_idx, "J", half=half).id],
                                      level, key)[0])  # fmt: skip
    if model == "R":
        r = results[Task(truth, seed, level_idx, "R", half=half).id]
        z = ln_z_task_second(r, level, key)
        x = np.array([s["p2"] for s in r["scans"]])
        return float(logsumexp(z, b=IDS._trap(x)) - math.log(x[-1] - x[0]))
    per = []
    nb = HALF_BLOCKS if half else 1
    for d in range(3) if distances is None else distances:
        zs = np.concatenate([ln_z_task_second(results[Task(truth, seed, level_idx, "P", d, b,
                                                           half).id], level, key)
                             for b in range(nb)])  # fmt: skip
        per.append(float(logsumexp(zs) - math.log(zs.size)))
    return IDS.planet_average(per) if len(per) > 1 else per[0]


def fit_record(results, truth: int, seed: int, level_idx: int, level: float) -> dict:
    z = {m: ln_z_model(results, truth, seed, level_idx, m, level) for m in MODELS}
    return {"truth": truth, "seed": seed, "level": level_idx, "ln_z": z}


# ---------------------------------------------------------------------------- V2′


def v2_points(results) -> list[dict]:
    """For each truth's noise-free fits of J, R and P: the best point and the points at
    Δχ² = 1, 4, 9 and 25 on either side along the strength, at the best radius or
    longitude (the planet: at the best distance and longitude)."""
    out = []
    for t in range(1, 5):
        for model in ("J", "R", "P"):
            if model == "P":
                cands = [(Task(t, 0, 0, "P", d, 0), s) for d in range(3)
                         for s in results[Task(t, 0, 0, "P", d, 0).id]["scans"]]  # fmt: skip
            else:
                tk = Task(t, 0, 0, model)
                cands = [(tk, s) for s in results[tk.id]["scans"]]
            tk, s = min(cands, key=lambda c: c[1]["chi2_min"])
            pts = [p for p in Points(tk.dir / "points.jsonl").calls.values()
                   if p["kind"] == "fine" and round(p["p2"], 10) == round(s["p2"], 10)
                   and math.isfinite(p["chi2"])]  # fmt: skip
            pts.sort(key=lambda p: p["x"])
            x = np.array([p["x"] for p in pts])
            cs = CubicSpline(x, [p["chi2"] for p in pts])
            sub = np.linspace(x[0], x[-1], (len(x) - 1) * 200 + 1)
            c = cs(sub)
            ib = int(np.argmin(c))
            xs = [float(sub[ib])]
            for dchi in V2_DCHI:
                for side in (-1, 1):
                    rng_ = sub[:ib][::-1] if side < 0 else sub[ib:]
                    vv = c[:ib][::-1] if side < 0 else c[ib:]
                    j = np.nonzero(vv >= c[ib] + dchi)[0]
                    xs.append(float(rng_[j[0]]) if j.size else math.nan)
            for i, xv in enumerate(xs):
                out.append({"truth": t, "model": model, "task": tk.id, "p2": s["p2"],
                            "x": xv, "index": i})  # fmt: skip
    return out


def v2_one(point: dict) -> dict:
    """Laplace against importance sampling over the six elements at one point."""
    from scipy.stats import multivariate_t

    if not math.isfinite(point["x"]):
        return {**point, "pass": False, "note": "point outside the fine scan"}
    prep = prepare()
    tk = Task.from_id(point["task"])
    data = dataset(prep, tk.truth, 0, 1.0)
    pts = Points(tk.dir / "points.jsonl").calls.values()
    near = min((p for p in pts if round(p["p2"], 10) == round(point["p2"], 10)
                and math.isfinite(p["chi2"])), key=lambda p: abs(p["x"] - point["x"]))  # fmt: skip
    extra = tk.extra(point["x"], point["p2"])
    resid = data.resid_fn(tk.model, extra)
    c = lm_fit(resid, np.array(near["theta6"]), EL_STEPS)
    _, s, vt = np.linalg.svd(c.jac, full_matrices=False)
    tmat = (vt.T / s) * math.sqrt(V2_SCALE)
    ln_det_t = -float(np.sum(np.log(s))) + 3.0 * math.log(V2_SCALE)
    ln_lap = -0.5 * c.chi2 - c.ln_s + 3.0 * math.log(2.0 * math.pi)
    rng = np.random.default_rng([2026, 21, 1882, tk.truth, 0, MODELS.index(tk.model),
                                 point["index"], 2])  # fmt: skip
    prop = multivariate_t(loc=np.zeros(6), shape=np.eye(6), df=V2_NU)
    w = prop.rvs(size=V2_DRAWS, random_state=rng)
    lw = np.full(V2_DRAWS, -np.inf)
    for i in range(V2_DRAWS):
        try:
            r = resid(c.theta6 + tmat @ w[i])
            lw[i] = -0.5 * float(r @ r) - (prop.logpdf(w[i]) - ln_det_t)
        except NoContact:
            pass
    m = float(np.max(lw))
    wt = np.exp(lw - m)
    ln_is = m + math.log(float(wt.mean()))
    err = float(wt.std() / wt.mean() / math.sqrt(V2_DRAWS))
    diff = ln_lap - ln_is
    return {**point, "ln_laplace": ln_lap, "ln_is": ln_is, "diff": diff, "is_err": err,
            "converged": c.converged,
            "pass": bool(abs(diff) <= V2_TOL and err < V2_ERR)}  # fmt: skip


@ledger.recorded
def v2(n_workers: int, results: dict) -> dict:
    from concurrent.futures import ProcessPoolExecutor

    path = CACHE / "v2.json"
    if path.exists():
        return json.loads(path.read_text())
    rec = ledger.start(f"{JOB}/V2", "gate", n_workers=n_workers, scratch=not _real())
    pts = v2_points(results)
    _log(f"V2': {len(pts)} points, {V2_DRAWS} integrations each, {n_workers} workers")
    with ProcessPoolExecutor(n_workers) as ex:
        rows = list(ex.map(v2_one, pts))
    frac = float(np.mean([r["pass"] for r in rows]))
    out = {"points": rows, "fraction_pass": frac, "n_points": len(rows), "pass": frac >= V2_FRAC}
    path.write_text(json.dumps(out, indent=1, default=float))
    rec.done("gate passed" if out["pass"] else "gate failed",
             n_likelihood_calls=len(rows), n_integrations=len(rows) * V2_DRAWS,
             sim_years=None)  # fmt: skip
    return out


# ---------------------------------------------------------------------------- scheduling


def _spawn(task: Task) -> subprocess.Popen:
    task.dir.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1",
           "FORK_IDFULL_CACHE": str(CACHE)}  # fmt: skip
    cmd = [sys.executable, "-m", "fork.phase2.identifiability_full", "scan", task.id]
    with (task.dir / "log.txt").open("a") as fh:
        return subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT, env=env)


def schedule(tasks: list[Task], n_slots: int) -> bool:
    """One process per scan, the planet's first (the longest); one retry per failure."""
    order = sorted(tasks, key=lambda t: (t.model != "P", t.model != "R", not t.half))
    pending = [t for t in order if not t.done_path.exists()]
    _log(f"schedule: {len(pending)} of {len(tasks)} scans pending, {n_slots} slots")
    running: dict[str, tuple] = {}
    fails: dict[str, int] = {}
    while pending or running:
        for tid, (proc, t0, task) in list(running.items()):
            rc = proc.poll()
            if rc is None:
                continue
            del running[tid]
            if rc == 0 and task.done_path.exists():
                _log(f"done {tid} ({(time.time() - t0) / 60:.1f} min)")
            else:
                fails[tid] = fails.get(tid, 0) + 1
                _log(f"FAILED {tid} (exit {rc}, attempt {fails[tid]})")
                if fails[tid] < 2:
                    pending.insert(0, task)
        while pending and len(running) < n_slots:
            task = pending.pop(0)
            running[task.id] = (_spawn(task), time.time(), task)
            _log(f"start {task.id}")
        time.sleep(20)
    ok = all(t.done_path.exists() for t in tasks)
    _log(f"schedule finished; failures: {fails or 'none'}; all done: {ok}")
    return ok


def load_results(tasks: list[Task]) -> dict[str, dict]:
    return {t.id: json.loads(t.done_path.read_text()) for t in tasks}


# ---------------------------------------------------------------------------- stage 1


def gates_stage1(prep: dict, results: dict, n_workers: int) -> dict:
    lv1 = LEVELS[0]
    v3 = []
    for t, model, kw in ((2, "R", {}), (3, "P", {"distances": (1,)}), (4, "J", {})):
        for lv in LEVELS:
            a = ln_z_model(results, t, 0, 0, model, lv, **kw)
            b = ln_z_model(results, t, 0, 0, model, lv, half=True, **kw)
            v3.append({"truth": TRUTHS[t - 1], "model": model, "level": lv, "ln_z": a,
                       "ln_z_half": b, "pass": abs(a - b) < V3_TOL})  # fmt: skip
    v3.append({"truth": "a", "model": "G", "note": "no grid: one call", "pass": True})
    v4 = {}
    for t in range(1, 5):
        own = FIT_OF_TRUTH[TRUTHS[t - 1]]
        if own == "G":
            v4[TRUTHS[t - 1]] = results[Task(t, 0, 0, "NG").id]["calls"]["G"]["chi2"]
        else:
            ids = [k for k in results if k.startswith(f"t{t}_s0_l0_{own}") and "half" not in k]
            v4[TRUTHS[t - 1]] = min(min(s["chi2_spline_min"], s["chi2_min"])
                                    for k in ids for s in results[k]["scans"])  # fmt: skip
    v6 = v6_check(results)
    tails = [(k, s["p2"]) for k, r in results.items() for s in r.get("scans", [])
             if not s["tails_ok"]]  # fmt: skip
    out = {"V0": prep["V0"], "V1": prep["V1"],
           "V2'": v2(n_workers, results),
           "V3'": {"rows": v3, "pass": all(r["pass"] for r in v3)},
           "V4": {"chi2_min": v4, "pass": all(v < V4_TOL for v in v4.values())},
           "V6": v6, "tails_not_negligible": tails, "level_used": lv1}  # fmt: skip
    out["pass"] = all(out[k]["pass"] for k in ("V0", "V1", "V2'", "V3'", "V4", "V6"))
    return out


def v6_check(results: dict) -> dict:
    """V6 on every call (as frozen) and on the calls within Δχ² = V6_FAR of their scan's best
    point (amendment 1, item 9); ``V6_AMENDED`` says which decides the gate. A call with no
    contact is not a failure to converge (item 3); N and G are each their own best point."""
    n = n_in = bad = bad_in = nocontact = 0
    near, far = [], []
    for k in results:
        tk = Task.from_id(k)
        calls = list(Points(tk.dir / "points.jsonl").calls.values())
        best = {round(s["p2"], 10): s["chi2_min"] for s in results[k].get("scans", [])}
        for c in calls:
            if c["no_contact"]:
                nocontact += 1
                continue
            b = c["chi2"] if tk.model == "NG" else best.get(round(c["p2"], 10), math.inf)
            inside = c["chi2"] - b <= V6_FAR
            n += 1
            n_in += inside
            if c["converged"]:
                continue
            bad += 1
            bad_in += inside
            row = {"task": k, "p2": c["p2"], "x": c["x"], "chi2": c["chi2"], "dchi2": c["chi2"] - b}
            (near if c["chi2"] - b <= V6_DCHI else far if not inside else []).append(row)
    frac_all = bad / n if n else 0.0
    frac_in = bad_in / n_in if n_in else 0.0
    frac = frac_in if V6_AMENDED else frac_all
    return {"calls": n, "not_converged": bad, "fraction_all": frac_all,
            "calls_within_far": n_in, "not_converged_within_far": bad_in,
            "fraction_within_far": frac_in, "rule": "amended" if V6_AMENDED else "frozen",
            "no_contact_calls": nocontact, "near_best": near, "far_not_converged": len(far),
            "pass": frac <= V6_FRAC and not near}  # fmt: skip


def _noise_free_r1(f: dict, truth: str) -> dict:
    own = FIT_OF_TRUTH[truth]
    lead = min(f[own] - f[m] for m in MODELS if m != own)
    return {"noise_free_only": True, "winner": max(f, key=f.get), "lead_of_truth": lead,
            "grade": IDS.grade(abs(lead))}  # fmt: skip


def readings(results: dict, seeds: tuple[int, ...]) -> dict:
    """The earlier plan's § 7: R1 at the calibrated level (noise-free only without seeds),
    R2, R3 on noise-free data, R4 on noise-free data at the three levels (its noisy part is
    Stage 3, not in this plan)."""
    out: dict = {"R1": {}, "R2": {}, "R3": {}, "R4": {}, "ln_b_vs_G": {}}
    for ti, t in enumerate(TRUTHS, start=1):
        for li, lv in enumerate(LEVELS):
            free = fit_record(results, ti, 0, 0, lv)
            out["R4"][f"{t}/level{li}"] = _noise_free_r1(free["ln_z"], t)
            out["ln_b_vs_G"][f"{t}/level{li}"] = {m: free["ln_z"][m] - free["ln_z"]["G"]
                                                  for m in MODELS}  # fmt: skip
        fits = [fit_record(results, ti, 0, 0, LEVELS[0])]
        fits += [fit_record(results, ti, s, 0, LEVELS[0]) for s in seeds]
        out["R1"][t] = IDS.read_r1(fits, t) if seeds else _noise_free_r1(fits[0]["ln_z"], t)
        out["R2"][t] = [{"seed": f["seed"], "ln_b_best_vs_N":
                         max(f["ln_z"][m] for m in MODELS[1:]) - f["ln_z"]["N"]}
                        for f in fits]  # fmt: skip
        g = fits[0]["ln_z"]["G"]
        out["R3"][t] = {}
        for m in ("J", "R", "P"):
            rows = []
            for d in R3_WIDTHS:
                z = ln_z_model(results, ti, 0, 0, m, LEVELS[0], f"r3_{d:g}")
                rows.append({"decades": d, "ln_b_vs_G": z - g, "grade": IDS.grade(z - g)})
            out["R3"][t][m] = rows
    a_free = out["ln_b_vs_G"]["a/level0"]
    out["R1_truth_a_G_against"] = {m: {"ln_b_G_over_m": -a_free[m], "grade": IDS.grade(a_free[m])}
                                   for m in ("J", "R", "P")}  # fmt: skip
    return out


@ledger.recorded
def stage1(n_workers: int) -> dict:
    prep = prepare()
    tasks = stage1_tasks()
    if not schedule(tasks, n_workers):
        raise SystemExit("Stage 1: scans incomplete")
    results = load_results(tasks)
    rec = ledger.start(f"{JOB}/stage1/gates", "gate", n_workers=n_workers, scratch=not _real())
    gates = gates_stage1(prep, results, n_workers)
    out = {
        "plan": "docs/phase2/identifiability_full.md",
        "stage": 1,
        "sigma1_s": prep["sigma1_s"],
        "n_kept": prep["n_kept"],
        "gates": gates,
    }
    if gates["pass"]:
        out["readings"] = readings(results, ())
        out["readings"]["note"] = (
            "R1 and R4 on noise-free data only, reported as such and "
            "not as R1 or R4; R2 on noise-free data; R3 complete"
        )
    else:
        out["readings"] = "not made: a gate failed"
    path = CACHE / "stage1.json"
    path.write_text(json.dumps(out, indent=1, default=float))
    rec.done("gate passed" if gates["pass"] else "gate failed", n_likelihood_calls=None,
             n_integrations=None, sim_years=None)  # fmt: skip
    _log(f"Stage 1: gates {'pass' if gates['pass'] else 'FAIL'}")
    return out


@ledger.recorded
def stage2(n_workers: int) -> dict:
    s1 = CACHE / "stage1.json"
    if not s1.exists() or not json.loads(s1.read_text())["gates"]["pass"]:
        raise SystemExit("Stage 2 runs only if Stage 1's gates pass")
    tasks = stage2_tasks()
    if not schedule(tasks, n_workers):
        raise SystemExit("Stage 2: scans incomplete")
    results = load_results(stage1_tasks() + tasks)
    rec = ledger.start(f"{JOB}/stage2/summary", "summarize", n_workers=1, scratch=not _real())
    v6 = v6_check({t.id: results[t.id] for t in tasks})
    out = {"plan": "docs/phase2/identifiability_full.md", "stage": 2, "V6_stage2": v6,
           "readings": readings(results, SEEDS)}  # fmt: skip
    (CACHE / "stage2.json").write_text(json.dumps(out, indent=1, default=float))
    rec.done(n_likelihood_calls=None, n_integrations=None, sim_years=None)
    _log("Stage 2: summary written")
    return out


def summarize() -> None:
    out = {}
    for name in ("stage1", "stage2"):
        p = CACHE / f"{name}.json"
        if p.exists():
            out[name] = json.loads(p.read_text())
    target = OUT_JSON if _real() else CACHE / "identifiability_full.json"
    target.write_text(json.dumps(out, indent=1, default=float) + "\n")
    _log(f"summary: {target}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("prepare", "stage1", "stage2", "scan", "plan", "summarize"))
    ap.add_argument("task", nargs="?")
    ap.add_argument("--n-workers", type=int, default=12)
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        with ledger.run(f"{JOB}/prepare", "prepare", n_workers=1, scratch=not _real()) as rec:
            p = prepare()
            ok = p["V0"]["pass"] and p["V1"]["pass"]
            rec.done("gate passed" if ok else "gate failed")
    elif a.cmd == "stage1":
        stage1(a.n_workers)
    elif a.cmd == "stage2":
        stage2(a.n_workers)
    elif a.cmd == "scan":
        run_task(Task.from_id(a.task))
    elif a.cmd == "plan":
        for name, ts in (("Stage 1", stage1_tasks()), ("Stage 2", stage2_tasks())):
            done = sum(t.done_path.exists() for t in ts)
            print(f"{name}: {done}/{len(ts)} scans done")
    else:
        summarize()


if __name__ == "__main__":
    main()
