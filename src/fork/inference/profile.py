"""Profile likelihood of the perturber's 1846-09-23 geocentric longitude (M5a diagnostics).

At each target longitude L (light-time-corrected geocentric, ecliptic of date, on Galle's night)
every other free parameter is re-optimised under the constraint lon(θ) = L, and the minimum χ²
is recorded.

The constraint is eliminated exactly. The perturber's λ₀ is not an LM variable: inside every
residual evaluation it is solved from the others by a warm-started secant iteration, to
|lon − L| < 1e-12 rad (~2e-7″). This is well defined because, at fixed Earth position and orbit
shape, the geocentric longitude of a body outside Earth's orbit is monotone in its orbital
phase. (A first attempt with a stiff penalty residual was so ill-conditioned that LM crawled and
hit its evaluation cap.)

Starts per grid point (``n_starts`` ≥ 8). Uranus always starts from the data-only no-perturber
fit and m from 10 M⊕. Then:
  - free-a mode: a ∈ {25, 30, 38.4, 50} AU × ϖ ∈ {0°, 180°}, e = 0.04;
  - fixed-a mode: ϖ on an 8-point grid, e = 0.1.
After the independent starts, continuation passes re-start each grid point from its neighbours'
best solutions until no point improves by more than 0.01 in χ², so one basin is followed
smoothly across the grid.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from fork.data.frames import (
    geocentric_lon_lat_ecl_j2000_light_time,
    precession_matrix_ecl_j2000_to_ecl_of_date,
)
from fork.data.synthetic import LongitudeObservations
from fork.inference.fit import (
    _PENALTY,
    START_LOG10_M_EARTH,
    TWO_PI,
    PerturberSpace,
    _outside_guards,
    _unpack_perturber_z,
    fit_uranus_elements,
    run_pool,
    uranus_theta_to_x,
    wrap_rad,
)
from fork.model.forward import ForwardContext, predict_states_ecl_j2000, predict_uranus_lon

FREE_A_STARTS_AU = (25.0, 30.0, 38.4, 50.0)
FREE_A_STARTS_VARPI_DEG = (0.0, 180.0)
FREE_A_START_E = 0.04
FIXED_A_STARTS_VARPI_DEG = tuple(45.0 * k for k in range(8))
FIXED_A_START_E = 0.1
CONTINUATION_TOL_CHI2 = 0.01
MAX_CONTINUATION_PASSES = 60  # stops earlier, once no point improves


@dataclass(frozen=True)
class GalleGeometry:
    """What the constraint needs besides θ: Galle's epoch, Earth there, and the precession
    rotation ecl_j2000 → ecl of date at that epoch."""

    jd_tdb: float
    earth_pos_ecl_j2000_au: np.ndarray
    r_ecl_j2000_to_of_date: np.ndarray

    @classmethod
    def from_horizons(cls) -> GalleGeometry:
        from fork.model.hidden_mass import earth_pos_about_inner_barycenter, galle_jd_tdb

        jd = galle_jd_tdb()
        return cls(jd, earth_pos_about_inner_barycenter(jd),
                   np.asarray(precession_matrix_ecl_j2000_to_ecl_of_date(jd)))  # fmt: skip

    def lon_ecl_of_date_rad(self, pos_ecl_j2000: np.ndarray, vel_ecl_j2000: np.ndarray) -> float:
        earth = self.earth_pos_ecl_j2000_au
        lon, lat, _ = geocentric_lon_lat_ecl_j2000_light_time(pos_ecl_j2000, vel_ecl_j2000, earth)
        u = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
        w = self.r_ecl_j2000_to_of_date @ u
        return math.atan2(w[1], w[0]) % TWO_PI


def galle_lon_of_date_rad(theta: np.ndarray, context: ForwardContext, geo: GalleGeometry) -> float:
    pos, vel = predict_states_ecl_j2000(theta, [geo.jd_tdb], context)
    return geo.lon_ecl_of_date_rad(pos[0, 3], vel[0, 3])


class _Lam0Solver:
    """Solve the perturber's λ₀ for the constraint, warm-started from the last solution."""

    def __init__(self, lam0_guess: float, tol_rad: float = 1e-12, max_iter: int = 12):
        self.lam0, self.slope = float(lam0_guess), 1.0
        self.tol, self.max_iter = tol_rad, max_iter

    def solve(self, z_rest, pq, space, context, geo, target_rad) -> tuple[float, float] | None:
        """(λ₀, |miss| [rad]) or None if the iteration fails / leaves the guards."""

        def miss(lam0: float) -> float | None:
            theta = _unpack_perturber_z(np.concatenate([z_rest, [lam0]]), pq, space)
            if _outside_guards(theta):
                return None
            return float(wrap_rad(galle_lon_of_date_rad(theta, context, geo) - target_rad))

        l0 = self.lam0
        f0 = miss(l0)
        if f0 is None:
            return None
        slope = self.slope
        for _ in range(self.max_iter):
            if abs(f0) < self.tol:
                self.lam0 = l0 % TWO_PI
                return self.lam0, abs(f0)
            l1 = l0 - f0 / slope
            f1 = miss(l1)
            if f1 is None:
                return None
            if abs(l1 - l0) > 1e-14 and f1 != f0:
                s_new = (f1 - f0) / (l1 - l0)
                if 0.05 < s_new < 20.0:
                    slope = s_new
            l0, f0 = l1, f1
            self.slope = slope
        return None


def _constrained_residuals(z_rest, pq, space, obs, context, geo, target_rad, solver) -> np.ndarray:
    n = obs.jd_tdb.size
    sol = solver.solve(z_rest, pq, space, context, geo, target_rad)
    if sol is None:
        return np.full(n, _PENALTY)
    theta = _unpack_perturber_z(np.concatenate([z_rest, [sol[0]]]), pq, space)
    lon_u = predict_uranus_lon(theta, obs.jd_tdb, context)
    return wrap_rad(obs.lon_ecl_j2000_rad - lon_u) / obs.sigma_rad


@dataclass(frozen=True)
class ConstrainedFit:
    target_lon_ecl_of_date_rad: float
    z_rest: np.ndarray  # LM vector: Uranus 4 + free perturber variables except λ₀
    theta: np.ndarray  # 11-entry forward-model θ
    chi2: float
    constraint_miss_arcsec: float  # |lon(θ) − target|
    success: bool
    nfev: int


def _run_constrained(z_rest0, lam0_guess, pq, space, obs, context, geo, target_rad
                     ) -> ConstrainedFit:  # fmt: skip
    """One constrained LM run; e ≤ space.e_max by the same two-phase scheme as the point fits
    (``fork.inference.fit._fit_bounded``). ``z_rest`` is returned in the interior layout."""
    solver = _Lam0Solver(lam0_guess)

    def lm(sp: PerturberSpace, z0):
        args = (pq, sp, obs, context, geo, target_rad, solver)
        return least_squares(_constrained_residuals, z0, args=args, method="lm", x_scale="jac",
                             xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=400)  # fmt: skip

    inner = space.interior()
    used, sol = inner, lm(inner, z_rest0)
    nfev = int(sol.nfev)
    if space.bounds_e:
        fin = solver.solve(sol.x, pq, inner, context, geo, target_rad)
        if fin is not None:
            y = inner.y_from_z(np.concatenate([sol.x[4:], [fin[0]]]))
            if y[2] ** 2 + y[3] ** 2 > space.e_max:
                used = space.boundary()
                sol = lm(used, np.concatenate([sol.x[:4], used.z_from_y(y)[:-1]]))
                nfev += int(sol.nfev)
    fin = solver.solve(sol.x, pq, used, context, geo, target_rad)
    if fin is None:
        return ConstrainedFit(float(target_rad), sol.x.copy(), np.full(11, np.nan), np.inf,
                              np.inf, False, nfev)  # fmt: skip
    theta = _unpack_perturber_z(np.concatenate([sol.x, [fin[0]]]), pq, used)
    lon_u = predict_uranus_lon(theta, obs.jd_tdb, context)
    chi2 = float(np.sum((wrap_rad(obs.lon_ecl_j2000_rad - lon_u) / obs.sigma_rad) ** 2))
    miss = fin[1] * 206264.80624709636
    ok = bool(sol.success) and not _outside_guards(theta)
    y_fin = used.y_from_z(np.concatenate([sol.x[4:], [fin[0]]]))
    z_rest = np.concatenate([sol.x[:4], inner.z_from_y(y_fin)[:-1]])
    return ConstrainedFit(float(target_rad), z_rest, theta, chi2, miss, ok, nfev)


def _start_vectors(space: PerturberSpace, xu: np.ndarray, context, geo, target_rad: float,
                   n_starts: int) -> list[tuple[np.ndarray, float]]:  # fmt: skip
    """(z_rest, λ₀ guess) per start. λ₀ guess: the target minus the epoch→1846 mean motion."""
    from fork.model.hidden_mass import BODE_A_AU

    free_a = bool(space.free[1])
    if free_a:
        combos = [(a, w, FREE_A_START_E) for a in FREE_A_STARTS_AU for w in FREE_A_STARTS_VARPI_DEG]
    else:
        combos = [(space.held[1], w, FIXED_A_START_E) for w in FIXED_A_STARTS_VARPI_DEG]
    if n_starts > len(combos):
        raise ValueError(f"at most {len(combos)} starts in this mode")
    out = []
    for a, w_deg, e in combos[:n_starts]:
        s, w = math.sqrt(e), math.radians(w_deg)
        a_au = BODE_A_AU if a is None else a
        n_deg_day = 0.9856076686 / a_au**1.5
        lam0 = (target_rad - math.radians(n_deg_day * (geo.jd_tdb - context.epoch_jd_tdb))) % TWO_PI
        y = np.array([START_LOG10_M_EARTH, a_au, s * math.cos(w), s * math.sin(w), lam0])
        out.append((np.concatenate([xu, space.z_from_y(y)[:-1]]), float(lam0)))
    return out


@dataclass(frozen=True)
class ProfileResult:
    mode: str
    targets_lon_ecl_of_date_rad: np.ndarray
    fits: tuple[ConstrainedFit, ...]  # best per target
    n_runs: int  # LM runs in total (independent starts + continuation)
    n_continuation_passes: int

    @property
    def chi2(self) -> np.ndarray:
        return np.array([f.chi2 for f in self.fits])


def profile_galle_longitude(
    obs: LongitudeObservations,
    context: ForwardContext,
    inc_raan_ecl_j2000_rad: tuple[float, float],
    space: PerturberSpace,
    targets_lon_ecl_of_date_rad: np.ndarray,
    *,
    mode: str = "custom",
    n_starts: int = 8,
    n_workers: int = 1,
    geo: GalleGeometry | None = None,
) -> ProfileResult:
    """Min χ² over all free parameters with the 1846-09-23 longitude held at each target.

    ``targets_lon_ecl_of_date_rad`` must be sorted ascending; continuation uses grid neighbours.
    """
    if space.held[4] is not None:
        raise ValueError("λ₀ must be free: it is the parameter the constraint mostly determines")
    geo = GalleGeometry.from_horizons() if geo is None else geo
    targets = np.asarray(targets_lon_ecl_of_date_rad, dtype=np.float64)
    uranus_only = fit_uranus_elements(obs, context, inc_raan_ecl_j2000_rad)
    x = uranus_theta_to_x(uranus_only.theta)
    xu, pq = x[:4].copy(), x[4:].copy()

    jobs, owner = [], []
    for i, t in enumerate(targets):
        for z0, lam0 in _start_vectors(space, xu, context, geo, float(t), n_starts):
            jobs.append((z0, lam0, pq, space, obs, context, geo, float(t)))
            owner.append(i)
    results = run_pool(_run_constrained, jobs, n_workers)
    n_runs = len(results)
    best: list[ConstrainedFit | None] = [None] * targets.size
    for i, f in zip(owner, results, strict=True):
        if f.success and (best[i] is None or f.chi2 < best[i].chi2):
            best[i] = f

    changed = set(range(targets.size))
    passes = 0
    while changed and passes < MAX_CONTINUATION_PASSES:
        passes += 1
        jobs, owner = [], []
        for i, t in enumerate(targets):
            for j in (i - 1, i + 1):
                if j in changed and 0 <= j < targets.size and best[j] is not None:
                    lam0 = float(best[j].theta[10] + t - targets[j])
                    jobs.append((best[j].z_rest, lam0, pq, space, obs, context, geo, float(t)))
                    owner.append(i)
        results = run_pool(_run_constrained, jobs, n_workers)
        n_runs += len(results)
        changed = set()
        for i, f in zip(owner, results, strict=True):
            if f.success and (best[i] is None or f.chi2 < best[i].chi2 - CONTINUATION_TOL_CHI2):
                best[i] = f
                changed.add(i)
    if any(b is None for b in best):
        raise RuntimeError("a profile point has no successful constrained fit")
    return ProfileResult(mode, targets, tuple(best), n_runs, passes)  # type: ignore[arg-type]


def delta_chi2_intervals(lon_deg: np.ndarray, dchi2: np.ndarray, level: float = 1.0
                         ) -> list[tuple[float, float]]:  # fmt: skip
    """Intervals where the (sorted, piecewise-linear) profile Δχ² ≤ ``level``, with edges
    linearly interpolated. An interval touching the grid edge is open there (edge value)."""
    lon_deg, dchi2 = np.asarray(lon_deg, float), np.asarray(dchi2, float)
    inside = dchi2 <= level
    out, k = [], 0
    while k < lon_deg.size:
        if not inside[k]:
            k += 1
            continue
        j = k
        while j + 1 < lon_deg.size and inside[j + 1]:
            j += 1

        def cross(i0: int, i1: int) -> float:
            f0, f1 = dchi2[i0] - level, dchi2[i1] - level
            return float(lon_deg[i0] + (lon_deg[i1] - lon_deg[i0]) * f0 / (f0 - f1))

        lo = lon_deg[k] if k == 0 else cross(k - 1, k)
        hi = lon_deg[j] if j == lon_deg.size - 1 else cross(j, j + 1)
        out.append((float(lo), float(hi)))
        k = j + 1
    return out
