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
  seen from the barycentre, equals R☉ ∓ r☿ (interior, exterior): the row's own contact
  (ingress or egress), found as a bracketed root within a day of the predicted mid-transit
  (``contact_root``, amendment 2), with both bodies carried by their Keplerian
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
import itertools
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numba import njit
from scipy.interpolate import CubicSpline
from scipy.special import logsumexp

from fork import ledger
from fork.data.frames import elements_to_state
from fork.integrators import wh_extra as W
from fork.integrators.kepler import kepler_fg_scalar
from fork.model import inner_system as IS
from fork.model.bodies import G, masses_msun
from fork.phase2 import contact_root as CR
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
    # amendment 2: which contact each row is
    "data/historical/mercury/newcomb_1882_contact_identity.csv":
        "6b03d75401d2bb2b15f31aed11f1f772cfac79871bbf54f6e8b898cd638c6193",
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
MAX_OFFSET_DAYS = 1.0  # the window of a contact about its printed epoch (amendment 2)

# § 4: the likelihood and the evidence
DCHI_TOL = 1e-3
MAX_IT = 20
COARSE_STEP, COARSE_HALF = 0.1, 4.0  # decades
BOX_EPS = 1e-9  # decades: a point this close to a box edge is at it
# Amendment 3, second addendum: the spline's knots are converged fine points with contacts, no
# two closer than KNOT_TOL; fine points are never closer than H_MIN apart in one grid; a scan
# value is refined at most MAX_REFINE times.
KNOT_TOL = 1e-4  # decades
H_MIN = 1e-3  # decades (halved with the grid for V3′)
MAX_REFINE = 6
# Amendment 3, third addendum: the box edges are points of every integral; a call within
# EDGE_EPS of an edge is the call at it; the knot interval at an edge is sub-divided SUB_EDGE
# times (SUB elsewhere).
EDGE_EPS = 1e-9  # decades
SUB_EDGE = 200
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
N_EDGE = [0]  # Jacobian columns taken one-sided at the edge of the elements' domain (_difference)


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


CONTACT_IDENTITY = IDS.CONTACTS.with_name("newcomb_1882_contact_identity.csv")
SIDE_OF_CONTACT = {"I": CR.INGRESS, "II": CR.INGRESS, "III": CR.EGRESS, "IV": CR.EGRESS}


def contact_rows(path: Path = IDS.CONTACTS, identity: Path = CONTACT_IDENTITY) -> dict:
    """Epochs (JD, printed time taken as TDB, as the earlier plan), weights, and whether
    each contact is exterior (tables III and IV); and, from the contact-identity file
    (amendment 2), which contact each row is (I to IV) and so its side, ingress or egress."""
    ep = IDS.contact_epochs(path)
    table, row = [], []
    with path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            table.append(r["table"])
            row.append(int(r["row"]))
    table = np.array(table)
    with identity.open(newline="", encoding="utf-8") as fh:
        ident = {(r["table"], int(r["row"])): r for r in csv.DictReader(fh)}
    keys = list(zip(table.tolist(), row, strict=True))
    if sorted(ident) != sorted(keys):
        raise ValueError(f"{identity.name}: rows differ from {path.name}")
    which = np.array([ident[k]["contact"] for k in keys])
    exterior = np.isin(table, ("III", "IV"))
    if np.any(exterior != np.isin(which, ("I", "IV"))):
        raise ValueError(f"{identity.name}: a contact's kind disagrees with its table")
    return {"jd": ep["jd"], "weight": ep["weight"], "exterior": exterior, "table": table,
            "row": np.array(row), "contact": which,
            "side": np.array([SIDE_OF_CONTACT[c] for c in which])}  # fmt: skip


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


@njit(cache=True)
def _carried(x0: np.ndarray, v0: np.ndarray, mu: float, t: np.ndarray) -> np.ndarray:
    """Positions after offsets t [d] from the epoch states (x0, v0): the exact two-body motion
    about the Sun (amendment 2; the plan's carry by the epoch's Keplerian acceleration, a
    second-order expansion, is off by 0.009″ in the direction from the Earth at 2 h and 16″
    at a day)."""
    out = np.empty_like(x0)
    if mu == 0.0:  # straight lines (tests)
        for i in range(x0.shape[0]):
            out[i] = x0[i] + v0[i] * t[i]
        return out
    for i in range(x0.shape[0]):
        r = kepler_fg_scalar(x0[i, 0], x0[i, 1], x0[i, 2], v0[i, 0], v0[i, 1], v0[i, 2], mu, t[i])
        out[i, 0], out[i, 1], out[i, 2] = r[0], r[1], r[2]
    return out


def _gap(xm, vm, xe, ve, exterior: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Separation of the centres minus the contact radius [rad], at offsets t [days]."""
    m, e = _carried(xm, vm, MU_MERCURY, t), _carried(xe, ve, MU_EMB, t)
    gm, gs = m - e, -e
    sep = np.arctan2(np.linalg.norm(np.cross(gm, gs), axis=1), np.einsum("ij,ij->i", gm, gs))
    rc = R_SUN / np.linalg.norm(gs, axis=1) + np.where(exterior, 1.0, -1.0) * R_MERCURY / (
        np.linalg.norm(gm, axis=1)
    )
    return sep - rc


def contact_offsets(xm, vm, xe, ve, exterior: np.ndarray,
                    side: np.ndarray) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Offsets [s] of each row's contact (``side``: −1 ingress, +1 egress) from its printed
    epoch, and the root's residual [rad]; NaN where that contact does not exist within a day
    (``contact_root``: no transit, or a grazing one)."""
    exterior = np.asarray(exterior, dtype=bool)

    def gap(i: np.ndarray, t: np.ndarray) -> np.ndarray:
        return _gap(xm[i], vm[i], xe[i], ve[i], exterior[i], t)

    t, res = CR.contact_roots(gap, side, ROOT_TOL, MAX_OFFSET_DAYS)
    return t * 86400.0, res


def nearer_side(off_in: np.ndarray, off_out: np.ndarray) -> np.ndarray:
    """The side of whichever contact, ingress or egress, lies nearer the printed epoch, given a
    model's offsets [s] of both (NaN where none: then 0). Used only to check the contact
    identity against truth a (tests), never to choose a contact."""
    a, b = np.abs(off_in), np.abs(off_out)
    side = np.where(a <= b, CR.INGRESS, CR.EGRESS)
    side = np.where(np.isnan(a), np.where(np.isnan(b), 0, CR.EGRESS), side)
    return np.where(np.isnan(b) & ~np.isnan(a), CR.INGRESS, side)


class NoContact(Exception):
    """The model gives no contact within a day of a printed epoch: its likelihood is zero."""


class OutsideDomain(Exception):
    """Elements outside a > 0, 0 ≤ e < 1, where ``elements_to_state`` gives no orbit."""


@dataclass(frozen=True)
class Observation:
    jd: np.ndarray
    weight: np.ndarray
    exterior: np.ndarray
    ctx: object
    side: np.ndarray  # −1 ingress, +1 egress (amendment 2)


def observe(model: str, theta6, extra, obs: Observation) -> np.ndarray:
    """The model's contact times [s from the printed epochs] at the kept rows."""
    off, _ = contact_offsets(*states(model, theta6, extra, obs.ctx, obs.jd), obs.exterior,
                             obs.side)  # fmt: skip
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
    edge_steps: int = 0  # Jacobian columns of this call taken at the domain's edge (_difference)

    def ln_l(self, level: float = 1.0) -> float:
        """Laplace's marginal over the six elements at error level ``level`` (× σ₁), up to
        the constant common to every model: −χ²/(2 L²) − Σ ln s_k + 6 ln L."""
        if self.no_contact or not math.isfinite(self.chi2):
            return -math.inf
        return -0.5 * self.chi2 / level**2 - self.ln_s + 6.0 * math.log(level)

    def row(self) -> dict:
        return {"chi2": self.chi2, "ln_s": self.ln_s, "theta6": self.theta6.tolist(),
                "iterations": self.iterations, "n_int": self.n_int,
                "converged": self.converged, "no_contact": self.no_contact,
                "edge_steps": self.edge_steps}  # fmt: skip


def _in_domain(theta: np.ndarray) -> bool:
    """The elements' domain: a > 0 and 0 ≤ e < 1 (``elements_to_state``)."""
    return bool(theta[0] > 0.0 and 0.0 <= theta[1] < 1.0)


def _difference(resid, theta: np.ndarray, r: np.ndarray, d: np.ndarray) -> np.ndarray:
    """r(θ + d) − r(θ), the change of the residuals over the step d.

    By a forward difference, as the plan has it; a forward trial with no contact takes the
    backward difference r(θ) − r(θ − d), as before 38bde9f, and no contact on both sides is an
    error. Only at the edge of the elements' domain (a > 0, 0 ≤ e < 1, ``OutsideDomain``) is
    the step taken one-sided into it, shortened by halves along d until a trial is inside and
    the difference scaled back to the full step (38bde9f, kept as a fallback, amendment 2);
    within that fallback a trial with no contact is skipped like one outside the domain
    (amendment 3, addendum). Each such column is counted in ``N_EDGE``. Where it does not
    fire, the arithmetic is the plan's, bit for bit.
    """
    try:
        return resid(theta + d) - r
    except NoContact:
        return r - resid(theta - d)
    except OutsideDomain:
        pass
    f = 1.0
    for _ in range(80):
        for sign in (1.0, -1.0):
            try:
                rt = resid(theta + sign * f * d)
            except (OutsideDomain, NoContact):  # amendment 3, addendum: both are skipped here
                continue
            N_EDGE[0] += 1
            return sign * (rt - r) / f
        f *= 0.5
    raise ValueError("_difference: no step along d reaches an allowed trial")


def _jacobian(resid, theta: np.ndarray, r: np.ndarray, steps: np.ndarray) -> np.ndarray:
    cols = []
    for k, h in enumerate(steps):
        d = np.zeros(len(theta))
        d[k] = h
        cols.append(_difference(resid, theta, r, d) / h)
    return np.array(cols).T


def _whitened_jacobian(resid, theta: np.ndarray, r: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """Forward differences along the singular directions of ``basis``, each step one whitened
    unit (θ + v_k / s_k), returned in the elements' coordinates; one-sided into the elements'
    domain at its edge (``_difference``)."""
    _, s, vt = np.linalg.svd(basis, full_matrices=False)
    cols = [_difference(resid, theta, r, vt[k] / s[k]) for k in range(len(s))]
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
    n0, e0 = N_INT[0], N_EDGE[0]
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
    ln_s, n_edge = float(np.sum(np.log(s))), N_EDGE[0] - e0
    return Call(chi, ln_s, th, final, it, N_INT[0] - n0, converged, edge_steps=n_edge)


class Data:
    """One data set: contact times y [s] with errors σ_i = σ₁ / √w_i at level 1."""

    def __init__(self, obs: Observation, y: np.ndarray, sigma: np.ndarray) -> None:
        self.obs, self.y, self.sigma = obs, np.asarray(y, float), np.asarray(sigma, float)

    def resid_fn(self, model: str, extra):
        def resid(theta6: np.ndarray) -> np.ndarray:
            if not _in_domain(theta6):
                raise OutsideDomain
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
        out = json.loads(path.read_text())
        if "side" not in out:
            raise SystemExit(f"{path} predates amendment 2 (no contact sides): move it aside")
        return out
    if _real():
        check_all_inputs()
    tr = truths()
    theta6, ctx = base_problem()
    rows = contact_rows()
    a = tr["a"]
    side = rows["side"]
    off, res = contact_offsets(*states(a["model"], theta6, a["extra"], ctx, rows["jd"]),
                               rows["exterior"], side)  # fmt: skip
    keep = np.isfinite(off) & (np.abs(off) <= LEFT_OUT_S)
    left = [f"{t} {r}" for t, r, k in zip(rows["table"], rows["row"], keep, strict=True) if not k]
    v0 = {"left_out": left, "max_kept_offset_s": float(np.max(np.abs(off[keep]))),
          "max_root_residual_rad": float(np.max(res[keep])),
          "offsets_s": [None if not math.isfinite(x) else float(x) for x in off]}  # fmt: skip
    v0["pass"] = bool(left == ["III 12"] and v0["max_kept_offset_s"] <= V0_MAX_S
                      and v0["max_root_residual_rad"] < ROOT_TOL)  # fmt: skip
    obs = Observation(rows["jd"][keep], rows["weight"][keep], rows["exterior"][keep], ctx,
                      side[keep])  # fmt: skip
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
        "side": side.tolist(),
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
    if prep["side"] != rows["side"].tolist():
        raise SystemExit("prepared.json's contact sides differ from the contact-identity file")
    return Observation(rows["jd"][keep], rows["weight"][keep], rows["exterior"][keep], ctx,
                       rows["side"][keep])  # fmt: skip


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

    def relabel(self, d: dict, kind: str) -> dict:
        """The same call under another kind (the later line wins when the file is read)."""
        d = {**d, "kind": kind}
        self.calls[self.key(d["p2"], d["x"])] = d
        with self.path.open("a") as fh:
            fh.write(json.dumps(d) + "\n")
        return d

    def add(self, p2, x, call: Call, kind: str, retries: int = 0) -> dict:
        d = {"p2": float(p2), "x": float(x), "kind": kind, **call.row(), "retries": retries}
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
        # Amendment 3, A: no strength is computed outside the widest box any integral uses
        # (the coarse scan's 8 decades); the integrals are clipped to it.
        self.box = box(task.model, max(R3_WIDTHS)) if task.model != "NG" else None

    def call(self, p2: float, x: float, kind: str, warm=None) -> dict:
        got = self.points.get(p2, x)
        if got is not None:
            if kind == "fine" and got["kind"] == "coarse":
                # a coarse point the fine grid lands on (a grid anchored at the box's edge, on
                # the coarse step) is one of the fine points too (amendment 3)
                got = self.points.relabel(got, "fine")
            if not got["no_contact"]:
                self.warm = self._warm_of(got)
            return got
        th0, j0 = warm if warm is not None else self.warm
        extra = self.task.extra(x, p2)
        c = self.data.call(self.task.model, extra, th0, j0)
        retries = 0
        if c.no_contact:
            # Amendment 3, B: a start with no contact is not a fit. Retry from the nearest
            # converged call of this scan with contacts, then from the base elements; the
            # call has no contact only if every start has none.
            for th1, j1 in self._retry_starts(p2, x, th0):
                retries += 1
                c = self.data.call(self.task.model, extra, th1, j1)
                if not c.no_contact:
                    break
        if not c.no_contact:
            self.jacs[Points.key(p2, x)] = c.jac
            self.warm = (c.theta6, c.jac)
        else:
            self.warm = (th0, None)
        return self.points.add(p2, x, c, kind, retries)

    def _retry_starts(self, p2: float, x: float, tried: np.ndarray):
        """The starts amendment 3 (B) retries a call from, in order, skipping any already tried:
        the nearest call of this scan (same second parameter) that converged with contacts,
        then the base elements."""
        k2 = Points.key(p2, 0)[0]
        near = [d for d in self.points.calls.values()
                if Points.key(d["p2"], 0)[0] == k2
                and d["converged"] and not d["no_contact"]]  # fmt: skip
        starts = []
        if near:
            starts.append(self._warm_of(min(near, key=lambda d: (abs(d["x"] - x), d["x"]))))
        starts.append((self.start, None))
        seen = [np.asarray(tried)]
        for th, jac in starts:
            if any(np.array_equal(th, t) for t in seen):
                continue
            seen.append(th)
            yield th, jac

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
        for _ in range(MAX_REFINE):
            h = max(2.0 * FINE_HALF * w / (FINE_N - 1), H_MIN) * self.h_scale
            pts = self._grid(p2, x0, h, round(FINE_HALF * w / h), warm)
            xb, wb = self._best(pts)
            if wb < 0.7 * w or wb > 1.5 * w or abs(xb - x0) > w:
                x0, w = xb, wb
                warm = self._warm_of(min(pts, key=lambda p: abs(p["x"] - xb)))
                continue
            break
        lo, hi = self.box
        edge = self._edge_of(xb, h)
        if edge is not None:
            # Amendment 3 (a): a best value at the box's edge is a fine point, with the fine
            # points dense there (spacing no wider than the coarse scan's).
            x0 = edge
            h = max(min(h, COARSE_STEP * self.h_scale), H_MIN * self.h_scale)
            pts = self._grid(p2, x0, h, round(FINE_HALF * w / h), None)
        need = 2.0 * TAIL_LN * self.max_level**2
        lo_n = hi_n = round(FINE_HALF * w / h)
        for _ in range(4000):
            fin = {round(p["x"], 10): p for p in pts}
            xs = sorted(fin)
            cmin = min(p["chi2"] for p in knots(pts))
            grow = False
            if fin[xs[0]]["chi2"] - cmin < need and xs[0] > lo + BOX_EPS:
                lo_n += 1
                grow = True
            if fin[xs[-1]]["chi2"] - cmin < need and xs[-1] < hi - BOX_EPS:
                hi_n += 1
                grow = True
            if not grow:
                break
            pts = self._grid(p2, x0, h, (lo_n, hi_n), None)
        return pts, x0, w, h, (lo_n, hi_n)

    def _edge_of(self, x: float, h: float) -> float | None:
        """The box edge within one fine spacing of x, if any."""
        lo, hi = self.box
        if x - lo < h:
            return lo
        if hi - x < h:
            return hi
        return None

    def _side(self, x0: float, h: float, n: int, sign: int) -> list[float]:
        """x0 + sign·k·h for k = 1..n, stopping at the box: the last point is the box's edge
        when the step would pass it (amendment 3, A)."""
        lo, hi = self.box
        edge = hi if sign > 0 else lo
        out, prev = [], x0
        for k in range(1, n + 1):
            x = x0 + sign * k * h
            if (sign > 0 and x >= hi - BOX_EPS) or (sign < 0 and x <= lo + BOX_EPS):
                if abs(prev - edge) > BOX_EPS:
                    out.append(edge)
                break
            out.append(x)
            prev = x
        return out

    def _grid(self, p2, x0, h, n, warm) -> list[dict]:
        """Fine points x0 + k·h, k = −lo_n..hi_n, inside the box (amendment 3, A): the centre,
        then upward, then downward from the centre."""
        lo_n, hi_n = (n, n) if isinstance(n, int) else n
        lo, hi = self.box
        x0 = min(max(x0, lo), hi)
        centre = self.call(p2, x0, "fine", warm)
        out = [centre]
        for x in self._side(x0, h, hi_n, 1):
            out.append(self.call(p2, x, "fine"))
        down = self._side(x0, h, lo_n, -1)
        if down:
            self.warm = self._warm_of(centre)
        for x in down:
            out.append(self.call(p2, x, "fine"))
        return out

    @staticmethod
    def _best(pts: list[dict], with_min: bool = False):
        """Best value, width and (optionally) the spline's minimum, from the knots alone
        (``knots``: converged fine points with contacts, amendment 3's second addendum)."""
        ok = knots(pts)
        x = np.array([p["x"] for p in ok])
        if len(ok) < 2:
            xb = float(x[0])
            return (xb, COARSE_STEP, float(ok[0]["chi2"])) if with_min else (xb, COARSE_STEP)
        cs = CubicSpline(x, [p["chi2"] for p in ok])
        sub = np.linspace(x[0], x[-1], (len(x) - 1) * SUB + 1)
        xb = float(sub[int(np.argmin(cs(sub)))])
        curv = float(cs(xb, 2))
        wb = math.sqrt(2.0 / curv) if curv > 0 else float(x[-1] - x[0]) / 12.0
        return (xb, wb, float(cs(xb))) if with_min else (xb, wb)

    def tails(self, p2: float, xb: float, best: dict) -> list[dict]:
        """The tail points inside the box (amendment 3, A: outside it nothing enters an
        integral, so nothing is computed there)."""
        lo, hi = self.box
        out = []
        for dx in TAIL_DECADES:
            for sgn in (1.0, -1.0):
                x = xb + sgn * dx
                if lo - BOX_EPS <= x <= hi + BOX_EPS:
                    self.warm = self._warm_of(best)
                    out.append(self.call(p2, x, "tail"))
        return out

    def widen_to(self, p2: float, x0: float, h: float, n: tuple[int, int], target: float,
                 best: dict) -> tuple[int, int]:  # fmt: skip
        """Extend the fine points at spacing h toward ``target`` (a tail point whose integrand
        is not negligible), until they reach it or MAX_WIDEN points are added."""
        lo_n, hi_n = n
        lo, hi = self.box
        side = 1 if target > x0 else -1
        self.warm = self._warm_of(best)
        for _ in range(MAX_WIDEN):
            if (x0 + hi_n * h >= target) if side > 0 else (x0 - lo_n * h <= target):
                break
            if side > 0:
                hi_n += 1
                x = x0 + hi_n * h
                self.call(p2, min(x, hi), "fine")
                if x >= hi - BOX_EPS:
                    break
            else:
                lo_n += 1
                x = x0 - lo_n * h
                self.call(p2, max(x, lo), "fine")
                if x <= lo + BOX_EPS:
                    break
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
        best = min(knots(fine), key=lambda p: p["chi2"])
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
            best = min(knots(fine), key=lambda p: p["chi2"])
        xb, wb, spline_min = self._best(fine, with_min=True)
        lo, hi = self.box
        at_edge = ("low" if best["x"] <= lo + BOX_EPS else
                   "high" if best["x"] >= hi - BOX_EPS else None)  # fmt: skip
        return {"p2": p2, "x_best": xb, "width": wb, "h": h, "chi2_min": best["chi2"],
                "chi2_spline_min": spline_min,
                "best_point": best, "widened_to_tails": widened,
                "best_at_box_edge": at_edge}  # fmt: skip

    def run(self) -> list[dict]:
        out, guess, warm = [], None, None
        for p2 in self.task.second_grid():
            r = self.scan_second(float(p2), guess, warm)
            out.append(r)
            guess = (r["x_best"], r["width"])
            warm = self._warm_of(r["best_point"])
        return out


# ---------------------------------------------------------------------------- quadrature


def knots(pts: list[dict]) -> list[dict]:
    """The spline's knots (amendment 3, second addendum): fine points whose fit converged with
    contacts, in increasing x, and of any two closer than KNOT_TOL the one with the lower χ²
    (the first in x on a tie). Unconverged and no-contact points are never knots."""
    ok = sorted((p for p in pts if p["kind"] == "fine" and p["converged"]
                 and not p["no_contact"] and math.isfinite(p["chi2"])),
                key=lambda p: p["x"])  # fmt: skip
    out: list[dict] = []
    for p in ok:
        if out and p["x"] - out[-1]["x"] < KNOT_TOL:
            if p["chi2"] < out[-1]["chi2"]:
                out[-1] = p
            continue
        out.append(p)
    return out


def _box_knots(allp: list[dict], box: tuple[float, float]) -> list[dict]:
    """The knots for one box (amendment 3, third addendum): the scan's knots, plus each box
    edge that lies within their range as a knot of its own at exactly the edge, from the call
    computed there (the coarse scan has a point on every box edge) if it converged with
    contacts; fine knots within KNOT_TOL of such an edge give way to it."""
    kn = knots(allp)
    if not kn:
        return kn
    x0, x1 = kn[0]["x"], kn[-1]["x"]
    for e in box:
        if not x0 < e < x1:
            continue
        at = [p for p in allp if abs(p["x"] - e) <= EDGE_EPS and p["converged"]
              and not p["no_contact"] and math.isfinite(p["chi2"])]  # fmt: skip
        if at:
            kn = [p for p in kn if abs(p["x"] - e) >= KNOT_TOL] + [{**at[0], "x": float(e)}]
    return sorted(kn, key=lambda p: p["x"])


def ln_z_strength(pts: list[dict], level: float, box: tuple[float, float],
                  h_fine: float, method: str = "spline") -> float:  # fmt: skip
    """ln of (1/width of box) ∫ L dx over the box at one value of the second parameter.

    ``spline`` (the plan's rule, item 4, with amendment 3's addenda): ln L splined on the
    box's knots (``_box_knots``) and integrated on a sub-grid 20 times finer than the knots'
    spacing, SUB_EDGE times in the knot interval at each box edge; outside the knots, ln L
    interpolated linearly between the other points on steps of h/20 (h/SUB_EDGE in the step
    at each edge). The box edges themselves are always points of the grid, valued by the same
    interpolant, so nothing is dropped at an edge.

    ``pexp`` (smoke condition 7's independent rule): ln L linear between the same nodes (the
    box's knots, the other points inside the box, and the two edges) and exp integrated
    exactly on each interval. ``trapezoid``: the same nodes, L by the trapezoid rule (report
    only)."""
    lo, hi = box
    allp = sorted(pts, key=lambda p: p["x"])
    x = np.array([p["x"] for p in allp])
    ll = np.array([_ln_l(p, level) for p in allp])
    fin = np.isfinite(ll)
    kn = _box_knots(allp, box)
    if not kn:
        return -math.inf
    xk = np.array([p["x"] for p in kn])
    lk = np.array([_ln_l(p, level) for p in kn])
    cs = CubicSpline(xk, lk) if len(xk) > 1 else None

    def lin(xx: np.ndarray) -> np.ndarray:
        if xx.size == 0 or not fin.any():
            return np.full(xx.size, -np.inf)
        return np.interp(xx, x[fin], ll[fin], left=-np.inf, right=-np.inf)

    def value(xx: np.ndarray, inner: str) -> np.ndarray:
        inside = (xx >= xk[0]) & (xx <= xk[-1])
        if inner == "spline" and cs is not None:
            v_in = cs(np.clip(xx, xk[0], xk[-1]))
        else:
            v_in = np.interp(xx, xk, lk)
        return np.where(inside, v_in, lin(xx))

    if method == "spline":
        parts = []
        for a, b in itertools.pairwise(xk):
            n = SUB_EDGE if (a <= lo <= b or a <= hi <= b) else SUB
            parts.append(np.linspace(a, b, n + 1)[:-1])
        parts.append(xk[-1:])
        step = h_fine / SUB
        if lo < xk[0]:
            g = np.arange(xk[0] - step, lo - step, -step)[::-1]
            parts.append(g)
            parts.append(np.linspace(lo, lo + step, SUB_EDGE // SUB + 1))
        if hi > xk[-1]:
            g = np.arange(xk[-1] + step, hi + step, step)
            parts.append(g)
            parts.append(np.linspace(hi - step, hi, SUB_EDGE // SUB + 1))
        xs = np.unique(np.concatenate(parts))
        xs = np.unique(np.concatenate([xs[(xs > lo) & (xs < hi)], [lo, hi]]))
        vals = value(xs, "spline")
        return float(logsumexp(vals, b=trapezoid_weights(xs)) - math.log(hi - lo))
    if method not in ("pexp", "trapezoid"):
        raise ValueError(method)
    nodes = np.unique(np.concatenate([xk[(xk > lo) & (xk < hi)],
                                      x[fin & (x > lo) & (x < hi)
                                        & ((x < xk[0]) | (x > xk[-1]))],
                                      [lo, hi]]))  # fmt: skip
    vals = value(nodes, "linear")
    if method == "trapezoid":
        return float(logsumexp(vals, b=trapezoid_weights(nodes)) - math.log(hi - lo))
    terms = []
    for a, b, va, vb in zip(nodes[:-1], nodes[1:], vals[:-1], vals[1:], strict=True):
        if not (np.isfinite(va) and np.isfinite(vb)):
            continue  # an interval ending at a point of zero likelihood adds nothing
        d = abs(vb - va)
        terms.append(max(va, vb) + math.log(b - a)
                     + (math.log(-math.expm1(-d)) - math.log(d) if d > 1e-12 else 0.0))  # fmt: skip
    return float(logsumexp(terms) - math.log(hi - lo)) if terms else -math.inf


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
        "retried": sum(c.get("retries", 0) > 0 for c in calls),
        "retried_with_contacts": sum(c.get("retries", 0) > 0 and not c["no_contact"]
                                     for c in calls),
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
    at_edge = [(k, s["p2"], s["best_at_box_edge"]) for k, r in results.items()
               for s in r.get("scans", []) if s.get("best_at_box_edge")]  # fmt: skip
    retries = {k: (r["calls_summary"]["retried"], r["calls_summary"]["retried_with_contacts"])
               for k, r in results.items() if r["calls_summary"].get("retried")}  # fmt: skip
    out = {"V0": prep["V0"], "V1": prep["V1"],
           "V2'": v2(n_workers, results),
           "V3'": {"rows": v3, "pass": all(r["pass"] for r in v3)},
           "V4": {"chi2_min": v4, "pass": all(v < V4_TOL for v in v4.values())},
           "V6": v6, "tails_not_negligible": tails, "level_used": lv1,
           "best_at_box_edge": at_edge, "retried_calls": retries}  # fmt: skip
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


# --- Ground-truth spot check (amendment 3, third addendum) ---------------------------------
SPOT_TOL = 0.05  # ln Z: runner against ground truth, and ground truth's two rules
SPOT_MIN = 50  # extra fits converged with contacts, per scan
SPOT_EDGE_OFFSETS = (
    0.002,
    0.005,
    0.01,
    0.02,
)  # decades either side of each box edge (halved for V3′)
SPOT_STEEP, SPOT_PER = 6, (0.2, 0.4, 0.6, 0.8)  # steepest knot intervals; points inside each
SPOT_PEAK_N, SPOT_PEAK_HALF = 21, 2.5  # points across the peak, ± widths
SPOT_FLANK = 10.0  # Δ ln L below the peak at which a box edge on a flank enters the sample
SPOT_SEED, SPOT_RANDOM = 20261116, 5  # the gate reading's random part


def _boxes(model: str) -> list[tuple[str, tuple[float, float]]]:
    return [("prior", box(model))] + [(f"r3_{w:g}", box(model, w)) for w in R3_WIDTHS]


def _scan_points(task: Task, p2: float) -> list[dict]:
    k = Points.key(p2, 0.0)[0]
    pts = Points(task.dir / "points.jsonl").calls.values()
    return [p for p in pts if Points.key(p["p2"], 0.0)[0] == k]


def _steep_intervals(task: Task, kn: list[dict]) -> list[tuple[float, float]]:
    """Knot intervals inside the widest box, steepest first in ln L at the lowest level."""
    lo, hi = box(task.model, max(R3_WIDTHS))
    lv = min(task.levels())
    iv = [(abs(_ln_l(b, lv) - _ln_l(a, lv)), a["x"], b["x"])
          for a, b in itertools.pairwise(kn) if lo <= a["x"] and b["x"] <= hi]  # fmt: skip
    return [(a, b) for _, a, b in sorted(iv, reverse=True)]


def spot_xs(task: Task, scan: dict, pts: list[dict], n_steep: int = SPOT_STEEP) -> list[float]:
    """The extra strengths of a spot check: SPOT_EDGE_OFFSETS either side of every box edge,
    SPOT_PER of the way along each of the n_steep steepest knot intervals, and SPOT_PEAK_N
    across ± SPOT_PEAK_HALF widths of the best value unless it sits on a box edge; inside the
    widest box, none where the scan already has a point."""
    f = 0.5 if task.half else 1.0
    out = {round(e + sg * d * f, 6) for _, b in _boxes(task.model) for e in b
           for d in SPOT_EDGE_OFFSETS for sg in (-1, 1)}  # fmt: skip
    for a, b in _steep_intervals(task, knots(pts))[:n_steep]:
        out |= {round(a + q * (b - a), 6) for q in SPOT_PER}
    if not scan["best_at_box_edge"]:
        out |= {round(scan["x_best"] + k * scan["width"], 6)
                for k in np.linspace(-SPOT_PEAK_HALF, SPOT_PEAK_HALF, SPOT_PEAK_N)}  # fmt: skip
    lo, hi = box(task.model, max(R3_WIDTHS))
    have = {round(p["x"], 6) for p in pts}
    return sorted(x for x in out if lo <= x <= hi and x not in have)


def spot_check(task: Task, i: int) -> dict:
    """Real fits at the spot check's extra strengths in scan i of a finished task (the scan's
    own points are not touched), then, for every box and level, the runner's ln Z from the
    scan's points against the ground truth: ln Z by the same rule from the scan's points and
    the extra fits together, valid where that agrees with the piecewise-exponential rule on
    the same points within SPOT_TOL. Fewer than SPOT_MIN extra fits converged with contacts:
    SPOT_PER points in the next SPOT_STEEP steepest intervals, until there are or none is left.
    Written to the cache's spotcheck directory."""
    scan = json.loads(task.done_path.read_text())["scans"][i]
    pts = _scan_points(task, scan["p2"])
    data = dataset(prepare(), task.truth, task.seed, LEVELS[task.level])
    base, _ = base_problem()
    good = [p for p in pts if p["converged"] and not p["no_contact"]]
    extra: list[dict] = []

    def fit(x: float) -> None:
        near = sorted(good + [e for e in extra if e["converged"] and not e["no_contact"]],
                      key=lambda p: abs(p["x"] - x))[:1]  # fmt: skip
        starts = [np.array(near[0]["theta6"])] if near else []
        tried = 0
        for th in [*starts, base]:
            c = data.call(task.model, task.extra(x, scan["p2"]), th)
            tried += 1
            if not c.no_contact:
                break
        extra.append({"p2": scan["p2"], "x": x, "kind": "fine", **c.row(), "retries": tried - 1})

    n_steep = SPOT_STEEP
    for x in spot_xs(task, scan, pts, n_steep):
        fit(x)
    while sum(e["converged"] and not e["no_contact"] for e in extra) < SPOT_MIN:
        n_steep += SPOT_STEEP
        xs = [x for x in spot_xs(task, scan, pts, n_steep) if x not in {e["x"] for e in extra}]
        if not xs:
            break
        for x in xs:
            fit(x)
    dense = pts + extra
    rows = []
    for name, b in _boxes(task.model):
        for lv in task.levels():
            r = {"box": name, "level": lv,
                 "runner": ln_z_strength(pts, lv, b, scan["h"]),
                 "pexp": ln_z_strength(pts, lv, b, scan["h"], "pexp"),
                 "trapezoid": ln_z_strength(pts, lv, b, scan["h"], "trapezoid"),
                 "truth": ln_z_strength(dense, lv, b, scan["h"]),
                 "truth_pexp": ln_z_strength(dense, lv, b, scan["h"], "pexp")}  # fmt: skip
            both = [r["truth"], r["truth_pexp"], r["runner"]]
            if all(v == -math.inf for v in both):  # zero likelihood in the box by every rule
                r["valid"], r["pass"] = True, True
            else:
                r["valid"] = abs(r["truth"] - r["truth_pexp"]) < SPOT_TOL
                r["pass"] = r["valid"] and abs(r["runner"] - r["truth"]) < SPOT_TOL
            rows.append(r)
    n_ok = sum(e["converged"] and not e["no_contact"] for e in extra)
    out = {"task": task.id, "scan": i, "p2": scan["p2"], "extra_fits": len(extra),
           "extra_converged": n_ok,
           "extra_unconverged": sum(not e["converged"] and not e["no_contact"] for e in extra),
           "extra_no_contact": sum(e["no_contact"] for e in extra),
           "enough": n_ok >= SPOT_MIN, "rows": rows,
           "pass": n_ok >= SPOT_MIN and all(r["pass"] for r in rows), "extra": extra}  # fmt: skip
    d = CACHE / "spotcheck"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{task.id}__{i}.json").write_text(json.dumps(out, default=float))
    return out


def spot_sample(tasks: list[Task]) -> list[tuple[str, int, str]]:
    """The gate reading's spot-check sample: every scan whose best point inside the box sits on
    a box edge; every scan with a box that leaves out its best value and whose edge nearest it
    has ln L more than SPOT_FLANK below the scan's peak at some level; and SPOT_RANDOM scans
    drawn by default_rng(SPOT_SEED) from all scans (sorted by task, then scan index)."""
    picked: dict[tuple[str, int], str] = {}
    every = []
    for t in sorted(tasks, key=lambda t: t.id):
        if t.model == "NG" or not t.done_path.exists():
            continue
        for i, s in enumerate(json.loads(t.done_path.read_text())["scans"]):
            every.append((t.id, i))
            if s["best_at_box_edge"]:
                picked.setdefault((t.id, i), "best on a box edge")
                continue
            kn = knots(_scan_points(t, s["p2"]))
            if len(kn) < 2:
                continue
            xk = np.array([p["x"] for p in kn])
            for _, (lo, hi) in _boxes(t.model):
                if lo <= s["x_best"] <= hi:
                    continue
                e = lo if s["x_best"] < lo else hi
                for lv in t.levels():
                    lk = np.array([_ln_l(p, lv) for p in kn])
                    if xk[0] <= e <= xk[-1] and np.interp(e, xk, lk) < lk.max() - SPOT_FLANK:
                        picked.setdefault((t.id, i), "box edge on a flank")
    rng = np.random.default_rng(SPOT_SEED)
    for j in rng.choice(len(every), size=min(SPOT_RANDOM, len(every)), replace=False):
        picked.setdefault(every[int(j)], "random")
    return [(k[0], k[1], why) for k, why in sorted(picked.items())]


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "cmd",
        choices=(
            "prepare",
            "stage1",
            "stage2",
            "scan",
            "plan",
            "summarize",
            "spotcheck",
            "spotsample",
        ),
    )
    ap.add_argument("task", nargs="?")
    ap.add_argument("--scan", type=int, default=0)
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
    elif a.cmd == "spotcheck":
        r = spot_check(Task.from_id(a.task), a.scan)
        print(json.dumps({k: v for k, v in r.items() if k not in ("rows", "extra")}))
    elif a.cmd == "spotsample":
        for row in spot_sample(stage1_tasks()):
            print(*row)
    elif a.cmd == "plan":
        for name, ts in (("Stage 1", stage1_tasks()), ("Stage 2", stage2_tasks())):
            done = sum(t.done_path.exists() for t in ts)
            print(f"{name}: {done}/{len(ts)} scans done")
    else:
        summarize()


if __name__ == "__main__":
    main()
