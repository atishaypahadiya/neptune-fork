"""Universal-variable two-body (Kepler) propagator.

Propagates a state (r, v) about a point mass μ by Δt using the universal anomaly χ with
Stumpff functions c₂, c₃ and Lagrange f, g coefficients (Danby 1988, §6.9; Vallado 2013, §2.3).
One code path covers elliptic, near-circular, and hyperbolic orbits.

Robustness:
- Elliptic Δt is reduced modulo the period first, so "long Δt" means the same work as short Δt
  and χ stays within a known bracket |χ| ≤ 2π/√α (since |ΔM| ≤ π ⇒ |ΔE| < π + 2e < 2π).
- Newton iteration is safeguarded by that bracket (f(χ) is monotone: f′ = r > 0). Any Newton
  step that leaves the bracket is replaced by bisection, so convergence is guaranteed.
- The new state is formed as r = r₀ + (f − 1) r₀ + g v₀ (and likewise for v) to avoid losing
  digits in the many short steps a Wisdom–Holman integrator takes.

All functions are ``numba.njit`` and allocation-free in the scalar core.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

TWO_PI = 2.0 * math.pi
_SERIES_PSI = 1.0  # |ψ| below this uses the Taylor series for c₂, c₃
_MAX_ITER = 200


@njit(cache=True)
def stumpff_c2_c3(psi: float) -> tuple[float, float]:
    """Stumpff functions c₂(ψ) = (1 − cos√ψ)/ψ and c₃(ψ) = (√ψ − sin√ψ)/√ψ³, all real ψ."""
    if abs(psi) < _SERIES_PSI:
        # c2 = Σ (−ψ)^k / (2k+2)!,  c3 = Σ (−ψ)^k / (2k+3)!  — 13 terms, truncation ≪ ε for |ψ|<1
        return _series(psi, 2), _series(psi, 3)
    if psi > 0.0:
        s = math.sqrt(psi)
        half = math.sin(0.5 * s)
        c2 = 2.0 * half * half / psi  # half-angle form: no cancellation in 1 − cos
        c3 = (s - math.sin(s)) / (psi * s)
        return c2, c3
    y = math.sqrt(-psi)
    half = math.sinh(0.5 * y)
    c2 = 2.0 * half * half / (-psi)
    c3 = (math.sinh(y) - y) / (-psi * y)
    return c2, c3


@njit(cache=True)
def _fact(n: int) -> float:
    out = 1.0
    for i in range(2, n + 1):
        out *= i
    return out


@njit(cache=True)
def _series(psi: float, offset: int) -> float:
    """Σ_{k=0}^{12} (−ψ)^k / (2k + offset)!, by Horner with term ratios (exact factorial ratios)."""
    # term_k / term_{k-1} = −ψ / ((2k+offset−1)(2k+offset))
    acc = 1.0
    for k in range(12, 0, -1):
        acc = 1.0 - psi * acc / ((2 * k + offset - 1) * (2 * k + offset))
    return acc / _fact(offset)


@njit(cache=True)
def kepler_fg_scalar(
    x: float, y: float, z: float, vx: float, vy: float, vz: float, mu: float, dt: float
) -> tuple[float, float, float, float, float, float]:
    """Propagate one state by ``dt`` about ``mu``. Returns (x, y, z, vx, vy, vz)."""
    out = kepler_fg_diag(x, y, z, vx, vy, vz, mu, dt)
    return out[0], out[1], out[2], out[3], out[4], out[5]


@njit(cache=True)
def kepler_fg_diag(
    x: float, y: float, z: float, vx: float, vy: float, vz: float, mu: float, dt: float
) -> tuple[float, float, float, float, float, float, int, int]:
    """:func:`kepler_fg_scalar` plus diagnostics: (..., n_iterations, n_bisection_steps)."""
    n_iter = 0
    n_bisect = 0
    if dt == 0.0:
        return x, y, z, vx, vy, vz, n_iter, n_bisect
    r0 = math.sqrt(x * x + y * y + z * z)
    v2 = vx * vx + vy * vy + vz * vz
    sqrt_mu = math.sqrt(mu)
    sigma0 = (x * vx + y * vy + z * vz) / sqrt_mu  # r₀·v₀/√μ
    alpha = 2.0 / r0 - v2 / mu  # 1/a (> 0 elliptic)
    one_minus_ar0 = 1.0 - alpha * r0

    if alpha > 0.0:
        period = TWO_PI / (sqrt_mu * alpha * math.sqrt(alpha))
        dt = dt - period * math.floor(dt / period + 0.5)  # now |dt| ≤ P/2
        if dt == 0.0:
            return x, y, z, vx, vy, vz, n_iter, n_bisect
        chi_max = TWO_PI / math.sqrt(alpha)
        lo, hi = -chi_max, chi_max
    else:
        lo, hi = -math.inf, math.inf

    target = sqrt_mu * dt

    # Short steps: second-order guess from f(χ) ≈ r₀χ + σ₀χ²/2 − √μ Δt.
    # Long elliptic steps (> 5% of P): χ ≈ √μ α Δt, exact for circular orbits.
    chi1 = target / r0
    chi = chi1 - 0.5 * sigma0 * chi1 * chi1 / r0
    if alpha > 0.0 and abs(dt) > 0.05 * period:
        chi = target * alpha
    if not (lo < chi < hi):
        chi = chi1 if lo < chi1 < hi else 0.0

    if alpha <= 0.0:
        # hyperbolic/parabolic: build a finite bracket by expansion (f is increasing)
        step = abs(chi) + abs(chi1) + 1e-300
        lo_b, hi_b = chi - step, chi + step
        for _ in range(200):
            c2, c3 = stumpff_c2_c3(alpha * lo_b * lo_b)
            f_lo = sigma0 * lo_b * lo_b * c2 + one_minus_ar0 * lo_b**3 * c3 + r0 * lo_b - target
            if f_lo <= 0.0:
                break
            lo_b -= 2.0 * step
            step *= 2.0
        step = abs(chi) + abs(chi1) + 1e-300
        for _ in range(200):
            c2, c3 = stumpff_c2_c3(alpha * hi_b * hi_b)
            f_hi = sigma0 * hi_b * hi_b * c2 + one_minus_ar0 * hi_b**3 * c3 + r0 * hi_b - target
            if f_hi >= 0.0:
                break
            hi_b += 2.0 * step
            step *= 2.0
        lo, hi = lo_b, hi_b

    c2 = 0.5
    c3 = 1.0 / 6.0
    psi = 0.0
    r = r0
    for _ in range(_MAX_ITER):
        n_iter += 1
        psi = alpha * chi * chi
        c2, c3 = stumpff_c2_c3(psi)
        chi2 = chi * chi
        fval = sigma0 * chi2 * c2 + one_minus_ar0 * chi2 * chi * c3 + r0 * chi - target
        r = chi2 * c2 + sigma0 * chi * (1.0 - psi * c3) + r0 * (1.0 - psi * c2)  # = f′(χ)
        if fval == 0.0:
            break
        delta = fval / r
        # Converged? Test *before* the bracket check: at convergence lo or hi equals χ, and a
        # last sub-ulp Newton step would otherwise be misread as leaving the bracket.
        done = abs(delta) <= 4.0e-16 * abs(chi)
        if not done:
            if fval < 0.0:
                lo = chi
            else:
                hi = chi
            new = chi - delta
            if lo < new < hi:
                chi = new
                continue
            new = 0.5 * (lo + hi)  # bisection fallback
            n_bisect += 1
            done = hi - lo <= 4.0e-16 * max(abs(lo), abs(hi))
            chi = new
            if not done:
                continue
        else:
            chi -= delta
        psi = alpha * chi * chi
        c2, c3 = stumpff_c2_c3(psi)
        chi2 = chi * chi
        r = chi2 * c2 + sigma0 * chi * (1.0 - psi * c3) + r0 * (1.0 - psi * c2)
        break

    chi2 = chi * chi
    fm1 = -chi2 * c2 / r0  # f − 1
    g = dt - chi2 * chi * c3 / sqrt_mu
    fdot = sqrt_mu / (r * r0) * chi * (psi * c3 - 1.0)
    gdm1 = -chi2 * c2 / r  # ġ − 1

    nx = x + fm1 * x + g * vx
    ny = y + fm1 * y + g * vy
    nz = z + fm1 * z + g * vz
    nvx = vx + fdot * x + gdm1 * vx
    nvy = vy + fdot * y + gdm1 * vy
    nvz = vz + fdot * z + gdm1 * vz
    return nx, ny, nz, nvx, nvy, nvz, n_iter, n_bisect


@njit(cache=True)
def kepler_step(
    r0: np.ndarray, v0: np.ndarray, mu: float, dt: float
) -> tuple[np.ndarray, np.ndarray]:
    """Two-body propagation of (r0 [AU], v0 [AU/day]) by dt [day] about mu [AU³/day²]."""
    out = kepler_fg_scalar(r0[0], r0[1], r0[2], v0[0], v0[1], v0[2], mu, dt)
    r = np.empty(3)
    v = np.empty(3)
    r[0], r[1], r[2] = out[0], out[1], out[2]
    v[0], v[1], v[2] = out[3], out[4], out[5]
    return r, v


@njit(cache=True)
def kepler_drift_inplace(x: np.ndarray, v: np.ndarray, mu: float, dt: float) -> None:
    """Kepler-drift every row of x, v (shape (N, 3)) about the same mu, in place."""
    for i in range(x.shape[0]):
        out = kepler_fg_scalar(x[i, 0], x[i, 1], x[i, 2], v[i, 0], v[i, 1], v[i, 2], mu, dt)
        x[i, 0], x[i, 1], x[i, 2] = out[0], out[1], out[2]
        v[i, 0], v[i, 1], v[i, 2] = out[3], out[4], out[5]
