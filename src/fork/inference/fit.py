"""Least-squares fits of forward models to longitude observations.

M4: ``fit_uranus_elements`` — Uranus's elements with no perturber (A2, via ``fit_uranus_only``)
or with a fixed, planted perturber (A2b).
M5a: ``fit_perturber_multistart`` — Uranus's in-plane elements + a planar perturber, from a
24-point grid of perturber λ₀ starts, in Le Verrier mode (a fixed) or with a free.

Uranus is fitted in non-singular variables

    x = (a, k, h, λ₀, p, q),  k + i·h = e·exp(iϖ),  p + i·q = tan(i/2)·exp(iΩ)

so circular and planar orbits are interior points. Per project rule 2 the in-plane four
(a, k, h, λ₀) are always free and start from the data alone (mean motion from the post-1781
longitudes, e = 0); (p, q), i.e. i and Ω, are held at their epoch values.

Why hold i, Ω: from longitudes alone they enter only through the reduction to the ecliptic,
≈ −tan²(i/2)·sin 2u (~10″ for Uranus), which is quadratic in (p, q) and invariant under
Ω → Ω + π. Freed, they act as a spare second harmonic that absorbs part of a perturber's signal.
``free_inc_raan=True`` frees them anyway, as a diagnostic only; since the gradient in (p, q)
vanishes at i = 0, that mode starts from i = 1° at several Ω and keeps the best.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares

from fork.data.synthetic import LongitudeObservations
from fork.model.bodies import M_CENTRAL_MSUN, URANUS, G
from fork.model.forward import ForwardContext, predict_uranus_lon

TWO_PI = 2.0 * math.pi
_MU_URANUS = G * (M_CENTRAL_MSUN + URANUS.mass_msun)
DEFAULT_OMEGA_STARTS_DEG = (0.0, 45.0, 90.0, 135.0)  # Ω and Ω + 180° are equivalent
START_INC_DEG = 1.0


def wrap_rad(d: np.ndarray) -> np.ndarray:
    """Wrap-safe angle difference into [−π, π) (project rule 6)."""
    return (np.asarray(d) + math.pi) % TWO_PI - math.pi


def uranus_x_to_theta(x: np.ndarray) -> np.ndarray:
    a, k, h, lam0, p, q = x
    e = math.hypot(k, h)
    varpi = math.atan2(h, k) % TWO_PI
    tan_half_i = math.hypot(p, q)
    inc = 2.0 * math.atan(tan_half_i)
    raan = math.atan2(q, p) % TWO_PI if tan_half_i > 0 else 0.0
    return np.array([a, e, inc, raan, varpi, lam0 % TWO_PI])


def uranus_theta_to_x(theta: np.ndarray) -> np.ndarray:
    a, e, inc, raan, varpi, lam0 = theta
    t = math.tan(0.5 * inc)
    return np.array([a, e * math.cos(varpi), e * math.sin(varpi), lam0, t * math.cos(raan),
                     t * math.sin(raan)])  # fmt: skip


@dataclass(frozen=True)
class FitResult:
    theta: np.ndarray
    residual_rad: np.ndarray  # wrap(obs − model), unweighted
    chi2: float
    dof: int
    success: bool
    message: str
    nfev: int
    start: int  # index of the winning start

    @property
    def rms_arcsec(self) -> float:
        return float(np.sqrt(np.mean(self.residual_rad**2)) * 206264.80624709636)


def initial_uranus_theta(obs: LongitudeObservations, epoch_jd_tdb: float) -> np.ndarray:
    """Data-only starting point: circular, planar orbit with the observed mean motion.

    Linear fit λ(t) = λ₀ + n (t − epoch) to the unwrapped post-1781 longitudes (σ = 3″ era),
    then a = (μ / n²)^(1/3) from Kepler's third law.
    """
    sel = obs.sigma_rad <= np.median(obs.sigma_rad)
    order = np.argsort(obs.jd_tdb[sel])
    t = obs.jd_tdb[sel][order] - epoch_jd_tdb
    lon = np.unwrap(obs.lon_ecl_j2000_rad[sel][order])
    n, lam0 = np.polyfit(t, lon, 1)
    a = (_MU_URANUS / n**2) ** (1.0 / 3.0)
    return np.array([a, 0.0, 0.0, 0.0, 0.0, lam0 % TWO_PI])


def fit_uranus_only(
    obs: LongitudeObservations,
    model,
    inc_raan_ecl_j2000_rad: tuple[float, float],
    *,
    free_inc_raan: bool = False,
    omega_starts_deg: tuple[float, ...] = DEFAULT_OMEGA_STARTS_DEG,
) -> FitResult:
    """Fit Uranus alone (no perturber) in ``model``'s context; see ``fit_uranus_elements``."""
    return fit_uranus_elements(
        obs, model.context, inc_raan_ecl_j2000_rad,
        free_inc_raan=free_inc_raan, omega_starts_deg=omega_starts_deg,
    )  # fmt: skip


def fit_uranus_elements(
    obs: LongitudeObservations,
    context: ForwardContext,
    inc_raan_ecl_j2000_rad: tuple[float, float],
    *,
    fixed_perturber_theta: np.ndarray | None = None,
    free_inc_raan: bool = False,
    omega_starts_deg: tuple[float, ...] = DEFAULT_OMEGA_STARTS_DEG,
) -> FitResult:
    """Levenberg–Marquardt fit of Uranus's elements, optionally with a fixed planted perturber.

    a, e, ϖ, λ₀ are free; (i, Ω) are held at ``inc_raan_ecl_j2000_rad`` (the epoch values, see
    ``fork.model.forward.uranus_inc_raan_at_epoch``). ``fixed_perturber_theta`` (5 entries) is
    integrated alongside but not fitted. ``free_inc_raan=True`` fits all six Uranus elements from
    a data-only multistart — diagnostic only (project rule 2).
    ``FitResult.theta`` holds the 6 Uranus elements only.
    """
    pert = np.empty(0) if fixed_perturber_theta is None else np.asarray(fixed_perturber_theta)
    if pert.shape not in ((0,), (5,)):
        raise ValueError("fixed_perturber_theta must have 5 entries")

    def predict(theta_u: np.ndarray, t_jd: np.ndarray) -> np.ndarray:
        return predict_uranus_lon(np.concatenate([theta_u, pert]), t_jd, context)

    theta0 = initial_uranus_theta(obs, context.epoch_jd_tdb)
    free = np.ones(6, dtype=bool)
    if not free_inc_raan:
        free[4:] = False
    sigma = obs.sigma_rad

    def residuals(x_free: np.ndarray, x_full: np.ndarray) -> np.ndarray:
        x = x_full.copy()
        x[free] = x_free
        pred = predict(uranus_x_to_theta(x), obs.jd_tdb)
        return wrap_rad(obs.lon_ecl_j2000_rad - pred) / sigma

    starts = []
    if free_inc_raan:
        for om in omega_starts_deg:
            th = theta0.copy()
            th[2], th[3] = math.radians(START_INC_DEG), math.radians(om)
            starts.append(uranus_theta_to_x(th))
    else:
        th = theta0.copy()
        th[2], th[3] = inc_raan_ecl_j2000_rad
        starts.append(uranus_theta_to_x(th))

    best = None
    for k, x0 in enumerate(starts):
        sol = least_squares(residuals, x0[free], args=(x0,), method="lm", x_scale="jac",
                            xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=2000)  # fmt: skip
        x = x0.copy()
        x[free] = sol.x
        chi2 = float(np.sum(sol.fun**2))
        if best is None or chi2 < best[0]:
            best = (chi2, k, x, sol)
    chi2, k, x, sol = best
    theta = uranus_x_to_theta(x)
    resid = wrap_rad(obs.lon_ecl_j2000_rad - predict(theta, obs.jd_tdb))
    return FitResult(theta, resid, chi2, obs.jd_tdb.size - int(free.sum()), bool(sol.success),
                     str(sol.message), int(sol.nfev), k)  # fmt: skip


# ---------------------------------------------------------------------------------- M5a

# Perturber starting point, identical for every λ₀ grid start and derived from nothing about
# Neptune: mass at the middle of the M5 prior (log10 m ∈ [0, 2] → 10 M⊕), a at Bode's 38.4 AU,
# and a small nonzero √e (the e-gradient vanishes at k = h = 0, as with i for Uranus).
START_LOG10_M_EARTH = 1.0
START_SQRT_E_K_H = (0.1, 0.0)  # e = 0.01, ϖ = 0
N_LAM0_STARTS = 24
# Outside these the forward model is unphysical (unbound / degenerate); residuals are set huge so
# LM rejects the step. Distinct from the M5 prior, which is reported, not enforced, here.
_MAX_E = 0.95
_MIN_A_AU = 1.0
_MAX_LOG10_M = 4.0
_PENALTY = 1e6
# Two converged starts are the same local minimum if all of these agree. Loose in a: in free-a
# mode LM stops at slightly different points along the flat mass–distance valley.
MINIMUM_TOL = {"chi2": 0.5, "lam0_deg": 5.0, "log10_m": 0.05, "a_au": 3.0}
# A run whose perturber improves χ² over the no-perturber fit by less than this has effectively
# switched the perturber off ("null"), not found a minimum.
NULL_DCHI2 = 1.0


@dataclass(frozen=True)
class PerturberFit:
    """One LM run. ``theta`` is the 11-entry forward-model θ (Uranus 6 + perturber 5)."""

    theta: np.ndarray
    residual_rad: np.ndarray
    chi2: float
    dof: int
    success: bool
    nfev: int
    start_lam0_rad: float
    status: str = "converged"  # converged | stalled_at_guard | null | failed

    @property
    def mass_earth(self) -> float:
        return float(10.0 ** self.theta[6])

    @property
    def in_m5_prior(self) -> bool:
        """log10 m ∈ [0, 2], a ∈ [20, 60] AU, e ≤ 0.3 (PHASE1_SPEC M5)."""
        lm, a, e = self.theta[6], self.theta[7], self.theta[8]
        return bool(0.0 <= lm <= 2.0 and 20.0 <= a <= 60.0 and e <= 0.3)


@dataclass(frozen=True)
class LocalMinimum:
    best: PerturberFit
    n_starts: int  # how many λ₀ starts converged here
    start_lam0_deg: tuple[float, ...]


@dataclass(frozen=True)
class MultistartResult:
    mode: str  # "le_verrier" or "free_a"
    fits: tuple[PerturberFit, ...]  # one per start, in grid order
    minima: tuple[LocalMinimum, ...]  # distinct converged minima, best first
    uranus_only: FitResult  # the no-perturber fit that seeds Uranus's elements

    @property
    def best(self) -> PerturberFit:
        return self.minima[0].best

    def status_counts(self) -> dict[str, int]:
        out = {k: 0 for k in ("converged", "stalled_at_guard", "null", "failed")}
        for f in self.fits:
            out[f.status] += 1
        return out


def _same_minimum(f: PerturberFit, g: PerturberFit, free_a: bool) -> bool:
    tol = MINIMUM_TOL
    dlam = abs(float(wrap_rad(f.theta[10] - g.theta[10])))
    return (abs(f.chi2 - g.chi2) <= tol["chi2"] and math.degrees(dlam) <= tol["lam0_deg"]
            and abs(f.theta[6] - g.theta[6]) <= tol["log10_m"]
            and (not free_a or abs(f.theta[7] - g.theta[7]) <= tol["a_au"]))  # fmt: skip


def classify(fit: PerturberFit, chi2_no_perturber: float) -> str:
    """``failed`` (LM did not converge), ``stalled_at_guard`` (stopped against the e / a / m
    guard that keeps the forward model physical), ``null`` (perturber explains nothing:
    Δχ² < ``NULL_DCHI2`` vs no perturber) or ``converged`` (a genuine local minimum)."""
    th = fit.theta
    if not fit.success:
        return "failed"
    if (th[8] >= _MAX_E - 0.01 or th[1] >= _MAX_E - 0.01 or th[7] <= _MIN_A_AU + 0.01
            or th[6] >= _MAX_LOG10_M - 0.01):  # fmt: skip
        return "stalled_at_guard"
    if fit.chi2 > chi2_no_perturber - NULL_DCHI2:
        return "null"
    return "converged"


def group_minima(fits: list[PerturberFit], free_a: bool) -> tuple[LocalMinimum, ...]:
    """Cluster ``converged`` fits into distinct local minima (``MINIMUM_TOL``), best χ² first."""
    groups: list[list[PerturberFit]] = []
    for f in sorted((f for f in fits if f.status == "converged"), key=lambda f: f.chi2):
        for g in groups:
            if _same_minimum(f, g[0], free_a):
                g.append(f)
                break
        else:
            groups.append([f])
    return tuple(LocalMinimum(g[0], len(g), tuple(math.degrees(f.start_lam0_rad) for f in g))
                 for g in groups)  # fmt: skip


@dataclass(frozen=True)
class PerturberSpace:
    """Map between the LM perturber vector and the forward-model perturber θ (5 entries).

    Fit variables y = (log10 m [M⊕], a [AU], k, h, λ₀) with (k, h) = √e·(cos ϖ, sin ϖ). Entries
    of ``held`` that are not None are fixed at that y value; the rest are free, in order.
    ``e_max`` caps the eccentricity smoothly: e = e_max·sin²(√(k² + h²)/√e_max). That is ≈ k² + h²
    for small e, so the cap changes little well inside it, and LM never meets a hard wall. The
    map is periodic, so there is no flat plateau to stall on (a tanh cap was tried first: LM ran
    into its saturated region and stopped at the cap even when the true e was inside it).
    """

    held: tuple[float | None, ...] = (None, None, None, None, None)
    e_max: float | None = None

    @property
    def free(self) -> np.ndarray:
        return np.array([v is None for v in self.held])

    @property
    def n_free(self) -> int:
        return int(self.free.sum())

    def y_from_z(self, z_p: np.ndarray) -> np.ndarray:
        y = np.array([0.0 if v is None else v for v in self.held])
        y[self.free] = z_p
        return y

    def theta_from_z(self, z_p: np.ndarray) -> np.ndarray:
        log10_m, a, k, h, lam0 = self.y_from_z(z_p)
        r2 = k * k + h * h
        e = r2 if self.e_max is None else self.e_max * math.sin(math.sqrt(r2 / self.e_max)) ** 2
        return np.array([log10_m, a, e, math.atan2(h, k) % TWO_PI, lam0 % TWO_PI])

    def z_from_y(self, y: np.ndarray) -> np.ndarray:
        return np.asarray(y, dtype=np.float64)[self.free]


def le_verrier_space(a_fixed_au: float | None, e_max: float | None = None) -> PerturberSpace:
    """Le Verrier mode (a held) or free-a mode (``a_fixed_au=None``)."""
    return PerturberSpace((None, a_fixed_au, None, None, None), e_max)


def held_orbit_space(log10_m_earth: float, a_au: float, e: float, varpi_rad: float
                     ) -> PerturberSpace:  # fmt: skip
    """Mass and orbit shape held (e.g. Le Verrier's published solution); only λ₀ free."""
    s = math.sqrt(e)
    return PerturberSpace((log10_m_earth, a_au, s * math.cos(varpi_rad), s * math.sin(varpi_rad),
                           None))  # fmt: skip


def _unpack_perturber_z(z: np.ndarray, pq: np.ndarray, space: PerturberSpace) -> np.ndarray:
    """LM vector z → 11-entry θ. z = Uranus (a, k, h, λ₀) + the free perturber variables."""
    theta_u = uranus_x_to_theta(np.concatenate([z[:4], pq]))
    return np.concatenate([theta_u, space.theta_from_z(z[4:])])


def _outside_guards(theta: np.ndarray) -> bool:
    return bool(theta[1] >= _MAX_E or theta[8] >= _MAX_E or theta[7] <= _MIN_A_AU
                or theta[6] >= _MAX_LOG10_M)  # fmt: skip


def _perturber_residuals(z, pq, space, obs: LongitudeObservations, context) -> np.ndarray:
    theta = _unpack_perturber_z(z, pq, space)
    if _outside_guards(theta):
        return np.full(obs.jd_tdb.size, _PENALTY)
    pred = predict_uranus_lon(theta, obs.jd_tdb, context)
    return wrap_rad(obs.lon_ecl_j2000_rad - pred) / obs.sigma_rad


def _run_perturber_start(z0, pq, space, obs, context, lam0: float) -> PerturberFit:
    """One LM run (module level so a process pool can pickle it)."""
    args = (pq, space, obs, context)
    sol = least_squares(_perturber_residuals, z0, args=args, method="lm", x_scale="jac",
                        xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=3000)  # fmt: skip
    theta = _unpack_perturber_z(sol.x, pq, space)
    chi2 = float(np.sum(sol.fun**2))
    n_obs = obs.jd_tdb.size
    ok = bool(sol.success) and chi2 < 0.5 * n_obs * _PENALTY**2
    resid = wrap_rad(obs.lon_ecl_j2000_rad - predict_uranus_lon(theta, obs.jd_tdb, context))
    return PerturberFit(theta, resid, chi2, n_obs - z0.size, ok, int(sol.nfev), float(lam0))


def run_pool(fn, jobs: list[tuple], n_workers: int) -> list:
    """``[fn(*job) for job in jobs]``, in a process pool when ``n_workers > 1`` (same results)."""
    if n_workers > 1 and len(jobs) > 1:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            return list(pool.map(fn, *zip(*jobs, strict=True)))
    return [fn(*job) for job in jobs]


def fit_perturber_multistart(
    obs: LongitudeObservations,
    context: ForwardContext,
    inc_raan_ecl_j2000_rad: tuple[float, float],
    *,
    a_fixed_au: float | None = None,
    e_max: float | None = None,
    space: PerturberSpace | None = None,
    mode: str | None = None,
    n_lam0_starts: int = N_LAM0_STARTS,
    lam0_starts_rad: np.ndarray | None = None,
    n_workers: int = 1,
) -> MultistartResult:
    """Multistart Levenberg–Marquardt: Uranus (a, e, ϖ, λ₀) + perturber, from λ₀ grid starts.

    Free variables: Uranus (a, e·cos ϖ, e·sin ϖ, λ₀) and perturber (log10 m, [a], √e·cos ϖ,
    √e·sin ϖ, λ₀); a is held at ``a_fixed_au`` in Le Verrier mode, and ``e_max`` optionally caps
    the perturber's e (``PerturberSpace``). ``space`` overrides both, e.g. ``held_orbit_space``
    for a published solution with only λ₀ free. Uranus's i, Ω are held at
    ``inc_raan_ecl_j2000_rad`` (project rule 2). Every start shares the data-only Uranus
    solution from ``fit_uranus_elements`` and the neutral perturber start above (for any variable
    not held), differing only in perturber λ₀ (``n_lam0_starts`` equally spaced over 0–360°, or
    ``lam0_starts_rad``). ``n_workers > 1`` runs starts in a process pool; results are identical
    to serial.
    """
    from fork.model.hidden_mass import BODE_A_AU

    if space is None:
        space = le_verrier_space(a_fixed_au, e_max)
        if mode is None:
            mode = "free_a" if a_fixed_au is None else "le_verrier"
    elif a_fixed_au is not None or e_max is not None:
        raise ValueError("pass either space or a_fixed_au / e_max")
    mode = mode or "custom"
    uranus_only = fit_uranus_elements(obs, context, inc_raan_ecl_j2000_rad)
    xu0 = uranus_theta_to_x(uranus_only.theta)  # (a, k, h, λ₀, p, q); p, q held
    pq = xu0[4:].copy()
    if lam0_starts_rad is None:
        lam0_starts_rad = np.arange(n_lam0_starts) * (TWO_PI / n_lam0_starts)
    k0, h0 = START_SQRT_E_K_H
    jobs = []
    for lam0 in np.asarray(lam0_starts_rad, dtype=np.float64):
        y0 = np.array([START_LOG10_M_EARTH, BODE_A_AU, k0, h0, lam0])
        jobs.append((np.concatenate([xu0[:4], space.z_from_y(y0)]), pq, space, obs, context,
                     float(lam0)))  # fmt: skip
    fits = run_pool(_run_perturber_start, jobs, n_workers)
    fits = [replace(f, status=classify(f, uranus_only.chi2)) for f in fits]
    free_a = bool(space.free[1])
    return MultistartResult(mode, tuple(fits), group_minima(fits, free_a), uranus_only)
