"""Wisdom–Holman symplectic map in democratic heliocentric coordinates (Duncan, Levison & Lee 1998).

Canonical variables: Q_i = heliocentric positions, P_i = m_i·V_i with V_i the *barycentric*
velocities. The Hamiltonian splits into three exactly solvable parts:

    H_kepler = Σ_i [ m_i |V_i|²/2 − G m₀ m_i / |Q_i| ]   (μ = G·m₀ for every planet)
    H_sun    = |Σ_i m_i V_i|² / (2 m₀)                    (Sun-momentum "jump": Q_i += h·ΣmV/m₀)
    H_int    = −Σ_{i<j} G m_i m_j / |Q_i − Q_j|           (planet–planet "kick")

Two second-order compositions are provided:

``SCHEME_KICK_OUTER`` (default, PHASE1_SPEC M3)
    half-kick → Sun-momentum drift (h/2) → Kepler drift (h) → Sun-momentum drift (h/2) → half-kick

``SCHEME_DRIFT_OUTER`` (matches REBOUND WHFast, ``coordinates="democraticheliocentric"``)
    Kepler (h/2) → Sun-momentum (h/2) → kick (h) → Sun-momentum (h/2) → Kepler (h/2)

Both are symplectic, time-symmetric and second order. They are distinct maps whose errors differ
at O(ε·h²) (ε ~ planet/Sun mass), which at dt = 20 d is ~0.01″ in Uranus's longitude over a
century. ``SCHEME_DRIFT_OUTER`` exists so this code can be checked against REBOUND
algorithm-for-algorithm (A1); see docs/notes/M3.md.

Conventions for x0, v0, masses, t_out: see ``fork.integrators._nbody``.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from fork.integrators._nbody import check_inputs, plan_outputs
from fork.integrators.kepler import kepler_drift_inplace
from fork.model.bodies import G

SCHEME_KICK_OUTER = 0
SCHEME_DRIFT_OUTER = 1


@njit(cache=True)
def helio_to_dh(
    x_h: np.ndarray, v_h: np.ndarray, masses: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric (x, v) → democratic heliocentric (Q, V): Q = x, V = v − v_barycenter."""
    m_tot = masses.sum()
    vc = np.zeros(3)
    for i in range(x_h.shape[0]):
        vc += masses[i + 1] * v_h[i]
    vc /= m_tot
    q = x_h.copy()
    vb = v_h - vc
    return q, vb


@njit(cache=True)
def dh_to_helio_into(q: np.ndarray, vb: np.ndarray, masses: np.ndarray, out_x, out_v) -> None:
    """(Q, V) → heliocentric x = Q, v = V − V_sun with V_sun = −Σ m_i V_i / m₀."""
    p = np.zeros(3)
    for i in range(q.shape[0]):
        p += masses[i + 1] * vb[i]
    v_sun = -p / masses[0]
    for i in range(q.shape[0]):
        for k in range(3):
            out_x[i, k] = q[i, k]
            out_v[i, k] = vb[i, k] - v_sun[k]


@njit(cache=True)
def _kick(q: np.ndarray, vb: np.ndarray, masses: np.ndarray, h: float) -> None:
    """Planet–planet interaction: V_i += h Σ_{j≠i} G m_j (Q_j − Q_i)/|Q_j − Q_i|³."""
    n = q.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            dx = q[j, 0] - q[i, 0]
            dy = q[j, 1] - q[i, 1]
            dz = q[j, 2] - q[i, 2]
            r2 = dx * dx + dy * dy + dz * dz
            s = h * G / (r2 * np.sqrt(r2))
            si = s * masses[j + 1]
            sj = s * masses[i + 1]
            vb[i, 0] += si * dx
            vb[i, 1] += si * dy
            vb[i, 2] += si * dz
            vb[j, 0] -= sj * dx
            vb[j, 1] -= sj * dy
            vb[j, 2] -= sj * dz


@njit(cache=True)
def _jump(q: np.ndarray, vb: np.ndarray, masses: np.ndarray, h: float) -> None:
    """Sun-momentum drift: Q_i += h · Σ_j m_j V_j / m₀."""
    px = 0.0
    py = 0.0
    pz = 0.0
    for i in range(q.shape[0]):
        m = masses[i + 1]
        px += m * vb[i, 0]
        py += m * vb[i, 1]
        pz += m * vb[i, 2]
    s = h / masses[0]
    for i in range(q.shape[0]):
        q[i, 0] += s * px
        q[i, 1] += s * py
        q[i, 2] += s * pz


@njit(cache=True)
def step(q: np.ndarray, vb: np.ndarray, masses: np.ndarray, h: float, scheme: int) -> None:
    """One WH step of size h (may be negative) on DH (Q, V), in place."""
    mu0 = G * masses[0]
    if scheme == SCHEME_KICK_OUTER:
        _kick(q, vb, masses, 0.5 * h)
        _jump(q, vb, masses, 0.5 * h)
        kepler_drift_inplace(q, vb, mu0, h)
        _jump(q, vb, masses, 0.5 * h)
        _kick(q, vb, masses, 0.5 * h)
    elif scheme == SCHEME_DRIFT_OUTER:
        kepler_drift_inplace(q, vb, mu0, 0.5 * h)
        _jump(q, vb, masses, 0.5 * h)
        _kick(q, vb, masses, h)
        _jump(q, vb, masses, 0.5 * h)
        kepler_drift_inplace(q, vb, mu0, 0.5 * h)
    else:
        raise ValueError("unknown Wisdom–Holman scheme")


@njit(cache=True)
def integrate_states(
    x0: np.ndarray,
    v0: np.ndarray,
    masses: np.ndarray,
    t_out: np.ndarray,
    dt: float,
    scheme: int = SCHEME_KICK_OUTER,
) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric (positions, velocities) at t_out, each (len(t_out), N, 3)."""
    check_inputs(x0, v0, masses)
    idx_fwd, idx_bwd, n_steps, rem = plan_outputs(t_out, dt)
    out_x = np.empty((t_out.size, x0.shape[0], 3))
    out_v = np.empty_like(out_x)
    for sign, idx in ((1.0, idx_fwd), (-1.0, idx_bwd)):
        q, vb = helio_to_dh(x0, v0, masses)
        k = 0
        for j in idx:
            while k < n_steps[j]:
                step(q, vb, masses, sign * dt, scheme)
                k += 1
            if rem[j] > 0.0:
                qc, vc = q.copy(), vb.copy()
                step(qc, vc, masses, sign * rem[j], scheme)
                dh_to_helio_into(qc, vc, masses, out_x[j], out_v[j])
            else:
                dh_to_helio_into(q, vb, masses, out_x[j], out_v[j])
    return out_x, out_v


@njit(cache=True)
def integrate(
    x0: np.ndarray,
    v0: np.ndarray,
    masses: np.ndarray,
    t_out: np.ndarray,
    dt: float,
    scheme: int = SCHEME_KICK_OUTER,
) -> np.ndarray:
    """Heliocentric positions [t_out, body, 3]. Default dt for this project: 20 days."""
    return integrate_states(x0, v0, masses, t_out, dt, scheme)[0]
