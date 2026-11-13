"""Wisdom–Holman map with the Phase 2 extra forces: the Sun's 1PN term, the Sun's J2, a ring.

The map is ``wisdom_holman``'s drift-outer scheme with the extra accelerations added to the
kick. With every extra force off the step is that module's step, call for call, so the result
is bit-identical to ``wisdom_holman.integrate_states`` (tested). Nothing in Phase 1 imports
this module.

Extra forces, all centred on the central body and applied to every planet (the reaction on
the central body is neglected: it is of order planet mass / solar mass times terms that are
already 1e-7 of the Newtonian one):

1PN (``gr``)
    The Schwarzschild term in the PPN form of the IERS Conventions (2010), eq. 10.12 without
    the Lense–Thirring and de Sitter terms (the form of Nobili & Roxburgh 1986 and of
    REBOUNDx's ``gr`` for a single massive body)::

        a = gr · μ/(c² r³) · { [2(β+γ) μ/r − γ v²] r + 2(1+γ) (r·v) v },   μ = G m₀

    with r, v heliocentric. ``gr`` is a scale factor (0 off, 1 general relativity); β = γ = 1
    in general relativity. The perihelion advance per orbit is
    gr · (2 + 2γ − β)/3 · 6πμ / (c² a (1 − e²)). The term depends on velocity, so its kick is
    an implicit midpoint step in v at fixed r (time-symmetric; iterated to convergence).

J2 (``j2``)
    a = −(3/2) J₂ μ R² / r⁵ · { [1 − 5 (r·s)²/r²] r + 2 (r·s) s },  s the unit spin axis, R the
    reference radius. Potential per unit mass: Φ = μ J₂ R² (3 (r·s)²/r² − 1) / (2 r³).

Ring (``ring``)
    A uniform circular ring of mass m_r and radius R_r about the central body, with unit
    normal n. Potential Φ = −(2 G m_r/π) K(m) / A, A² = (R_r + ρ)² + z², m = 4 R_r ρ / A²
    (K the complete elliptic integral of the first kind, parameter m), ρ and z cylindrical
    about n. The acceleration is its exact gradient.

Parameters travel in one float array ``extra`` (``make_extra``), so the kernels stay
``numba.njit`` and pure-array.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from fork.integrators._nbody import check_inputs, plan_outputs
from fork.integrators.kepler import kepler_drift_inplace
from fork.integrators.wisdom_holman import _jump, _kick, dh_to_helio_into, helio_to_dh
from fork.model.bodies import G

# Speed of light [AU/day]: 299792.458 km/s × 86400 s / 149597870.7 km (IAU 2012 au).
C_AU_PER_DAY: float = 299792.458 * 86400.0 / 149597870.7
# Solar radius used with J2 [AU]: 696000 km, the reference radius of the J2 values quoted in
# the planetary ephemerides (e.g. Park et al. 2017, AJ 153, 121).
R_SUN_AU: float = 696000.0 / 149597870.7

# Layout of ``extra``.
I_GR, I_BETA, I_GAMMA = 0, 1, 2
I_J2, I_RSUN, I_SX, I_SY, I_SZ = 3, 4, 5, 6, 7
I_RING_M, I_RING_R, I_NX, I_NY, I_NZ = 8, 9, 10, 11, 12
N_EXTRA = 13
_GR_MAX_ITER = 30


def make_extra(
    *,
    gr: float = 0.0,
    beta: float = 1.0,
    gamma: float = 1.0,
    j2: float = 0.0,
    r_sun_au: float = R_SUN_AU,
    spin_axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
    ring_mass: float = 0.0,
    ring_radius_au: float = 0.0,
    ring_normal: tuple[float, float, float] = (0.0, 0.0, 1.0),
) -> np.ndarray:
    """The parameter array of the extra forces. Axes are normalised here. Masses in M_sun."""
    e = np.zeros(N_EXTRA)
    e[I_GR], e[I_BETA], e[I_GAMMA] = gr, beta, gamma
    s = np.asarray(spin_axis, dtype=np.float64)
    n = np.asarray(ring_normal, dtype=np.float64)
    e[I_J2], e[I_RSUN] = j2, r_sun_au
    e[I_SX : I_SZ + 1] = s / np.linalg.norm(s)
    e[I_RING_M], e[I_RING_R] = ring_mass, ring_radius_au
    e[I_NX : I_NZ + 1] = n / np.linalg.norm(n)
    if ring_mass != 0.0 and not ring_radius_au > 0.0:
        raise ValueError("a ring with mass needs a positive radius")
    return e


@njit(cache=True)
def any_extra(extra: np.ndarray) -> bool:
    return extra[I_GR] != 0.0 or extra[I_J2] != 0.0 or extra[I_RING_M] != 0.0


@njit(cache=True)
def ellip_ke(m: float) -> tuple[float, float]:
    """Complete elliptic integrals K(m), E(m), parameter m = k² in [0, 1), by the
    arithmetic–geometric mean (Abramowitz & Stegun 17.6)."""
    a = 1.0
    b = math.sqrt(1.0 - m)
    c = math.sqrt(m)
    s = 0.5 * c * c
    p = 0.5
    for _ in range(40):
        if abs(c) < 1e-17:
            break
        an = 0.5 * (a + b)
        c = 0.5 * (a - b)
        b = math.sqrt(a * b)
        a = an
        p *= 2.0
        s += p * c * c
    k = math.pi / (2.0 * a)
    return k, k * (1.0 - s)


@njit(cache=True)
def gr_accel(r: np.ndarray, v: np.ndarray, mu: float, extra: np.ndarray, out: np.ndarray) -> None:
    """The 1PN acceleration at heliocentric r, v (module docstring)."""
    beta, gamma = extra[I_BETA], extra[I_GAMMA]
    r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]
    rn = math.sqrt(r2)
    v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
    rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2]
    f = extra[I_GR] * mu / (C_AU_PER_DAY * C_AU_PER_DAY * r2 * rn)
    a = 2.0 * (beta + gamma) * mu / rn - gamma * v2
    b = 2.0 * (1.0 + gamma) * rv
    for k in range(3):
        out[k] = f * (a * r[k] + b * v[k])


@njit(cache=True)
def position_accel(r: np.ndarray, mu: float, extra: np.ndarray, out: np.ndarray) -> None:
    """J2 plus ring acceleration at heliocentric r, written into out."""
    out[0] = 0.0
    out[1] = 0.0
    out[2] = 0.0
    r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]
    if extra[I_J2] != 0.0:
        rs = r[0] * extra[I_SX] + r[1] * extra[I_SY] + r[2] * extra[I_SZ]
        rn = math.sqrt(r2)
        f = -1.5 * extra[I_J2] * mu * extra[I_RSUN] * extra[I_RSUN] / (r2 * r2 * rn)
        a = 1.0 - 5.0 * rs * rs / r2
        for k in range(3):
            out[k] += f * (a * r[k] + 2.0 * rs * extra[I_SX + k])
    if extra[I_RING_M] != 0.0:
        rr = extra[I_RING_R]
        z = r[0] * extra[I_NX] + r[1] * extra[I_NY] + r[2] * extra[I_NZ]
        px = r[0] - z * extra[I_NX]
        py = r[1] - z * extra[I_NY]
        pz = r[2] - z * extra[I_NZ]
        rho = math.sqrt(px * px + py * py + pz * pz)
        a2 = (rr + rho) * (rr + rho) + z * z
        an = math.sqrt(a2)
        m = 4.0 * rr * rho / a2
        kk, ee = ellip_ke(m)
        # dK/dm = [E − (1 − m) K] / (2 m (1 − m)); π/8 as m → 0
        dk = math.pi / 8.0 if m < 1e-10 else (ee - (1.0 - m) * kk) / (2.0 * m * (1.0 - m))
        c = 2.0 * G * extra[I_RING_M] / math.pi
        dm_drho = 4.0 * rr / a2 - 8.0 * rr * rho * (rr + rho) / (a2 * a2)
        dm_dz = -8.0 * rr * rho * z / (a2 * a2)
        g_rho = c * (dk * dm_drho / an - kk * (rr + rho) / (a2 * an))
        g_z = c * (dk * dm_dz / an - kk * z / (a2 * an))
        if rho > 0.0:
            out[0] += g_rho * px / rho
            out[1] += g_rho * py / rho
            out[2] += g_rho * pz / rho
        out[0] += g_z * extra[I_NX]
        out[1] += g_z * extra[I_NY]
        out[2] += g_z * extra[I_NZ]


@njit(cache=True)
def position_potential(r: np.ndarray, mu: float, extra: np.ndarray) -> float:
    """J2 plus ring potential per unit mass at heliocentric r."""
    phi = 0.0
    r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]
    rn = math.sqrt(r2)
    if extra[I_J2] != 0.0:
        rs = r[0] * extra[I_SX] + r[1] * extra[I_SY] + r[2] * extra[I_SZ]
        phi += mu * extra[I_J2] * extra[I_RSUN] ** 2 * (3.0 * rs * rs / r2 - 1.0) / (2.0 * r2 * rn)
    if extra[I_RING_M] != 0.0:
        rr = extra[I_RING_R]
        z = r[0] * extra[I_NX] + r[1] * extra[I_NY] + r[2] * extra[I_NZ]
        rho = math.sqrt(max(r2 - z * z, 0.0))
        a2 = (rr + rho) * (rr + rho) + z * z
        kk, _ = ellip_ke(4.0 * rr * rho / a2)
        phi -= 2.0 * G * extra[I_RING_M] * kk / (math.pi * math.sqrt(a2))
    return phi


@njit(cache=True)
def _kick_position(q: np.ndarray, vb: np.ndarray, mu: float, h: float, extra: np.ndarray) -> None:
    acc = np.empty(3)
    for i in range(q.shape[0]):
        position_accel(q[i], mu, extra, acc)
        for k in range(3):
            vb[i, k] += h * acc[k]


@njit(cache=True)
def _kick_gr(
    q: np.ndarray, vb: np.ndarray, masses: np.ndarray, h: float, extra: np.ndarray
) -> None:
    """Implicit midpoint in v for the velocity-dependent 1PN term, positions fixed. The Sun's
    barycentric velocity is taken at the start of the kick."""
    mu = G * masses[0]
    vs = np.zeros(3)
    for i in range(q.shape[0]):
        for k in range(3):
            vs[k] -= masses[i + 1] * vb[i, k] / masses[0]
    v0 = np.empty(3)
    vm = np.empty(3)
    acc = np.empty(3)
    for i in range(q.shape[0]):
        for k in range(3):
            v0[k] = vb[i, k] - vs[k]
            vm[k] = v0[k]
        for _ in range(_GR_MAX_ITER):
            gr_accel(q[i], vm, mu, extra, acc)
            change = 0.0
            for k in range(3):
                new = v0[k] + 0.5 * h * acc[k]
                change += abs(new - vm[k])
                vm[k] = new
            if change <= 1e-17 * (abs(vm[0]) + abs(vm[1]) + abs(vm[2])):
                break
        gr_accel(q[i], vm, mu, extra, acc)
        for k in range(3):
            vb[i, k] += h * acc[k]


@njit(cache=True)
def step(q: np.ndarray, vb: np.ndarray, masses: np.ndarray, h: float, extra: np.ndarray) -> None:
    """One drift-outer WH step on DH (Q, V) with the extra forces in the kick."""
    mu0 = G * masses[0]
    kepler_drift_inplace(q, vb, mu0, 0.5 * h)
    _jump(q, vb, masses, 0.5 * h)
    if not any_extra(extra):
        _kick(q, vb, masses, h)
    else:
        position = extra[I_J2] != 0.0 or extra[I_RING_M] != 0.0
        _kick(q, vb, masses, 0.5 * h)
        if position:
            _kick_position(q, vb, mu0, 0.5 * h, extra)
        if extra[I_GR] != 0.0:
            _kick_gr(q, vb, masses, h, extra)
        if position:
            _kick_position(q, vb, mu0, 0.5 * h, extra)
        _kick(q, vb, masses, 0.5 * h)
    _jump(q, vb, masses, 0.5 * h)
    kepler_drift_inplace(q, vb, mu0, 0.5 * h)


# Fourth-order composition of the second-order step (Yoshida 1990, Phys. Lett. A 150, 262):
# S4(h) = S2(w1 h) S2(w0 h) S2(w1 h), w1 = 1/(2 − 2^{1/3}), w0 = 1 − 2 w1 (negative).
_W1: float = 1.0 / (2.0 - 2.0 ** (1.0 / 3.0))
_W0: float = 1.0 - 2.0 * _W1


@njit(cache=True)
def step_order(
    q: np.ndarray, vb: np.ndarray, masses: np.ndarray, h: float, extra: np.ndarray, order: int
) -> None:
    """One step of the chosen order (2 or 4)."""
    if order == 2:
        step(q, vb, masses, h, extra)
    elif order == 4:
        step(q, vb, masses, _W1 * h, extra)
        step(q, vb, masses, _W0 * h, extra)
        step(q, vb, masses, _W1 * h, extra)
    else:
        raise ValueError("order must be 2 or 4")


@njit(cache=True)
def integrate_states(
    x0: np.ndarray,
    v0: np.ndarray,
    masses: np.ndarray,
    t_out: np.ndarray,
    dt: float,
    extra: np.ndarray,
    order: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric (positions, velocities) at t_out, each (len(t_out), N, 3). Conventions as
    ``fork.integrators._nbody``; ``extra`` from ``make_extra``; ``order`` 2 (the Phase 1 map)
    or 4 (its Yoshida composition, three second-order steps per step)."""
    check_inputs(x0, v0, masses)
    if extra.size != N_EXTRA:
        raise ValueError("extra must come from make_extra")
    idx_fwd, idx_bwd, n_steps, rem = plan_outputs(t_out, dt)
    out_x = np.empty((t_out.size, x0.shape[0], 3))
    out_v = np.empty_like(out_x)
    for sign, idx in ((1.0, idx_fwd), (-1.0, idx_bwd)):
        q, vb = helio_to_dh(x0, v0, masses)
        k = 0
        for j in idx:
            while k < n_steps[j]:
                step_order(q, vb, masses, sign * dt, extra, order)
                k += 1
            if rem[j] > 0.0:
                qc, vc = q.copy(), vb.copy()
                step_order(qc, vc, masses, sign * rem[j], extra, order)
                dh_to_helio_into(qc, vc, masses, out_x[j], out_v[j])
            else:
                dh_to_helio_into(q, vb, masses, out_x[j], out_v[j])
    return out_x, out_v
