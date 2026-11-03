"""Shared N-body helpers for the integrators (all ``numba.njit``).

Conventions shared by ``rk4``, ``leapfrog`` and ``wisdom_holman``:

- ``x0, v0``: shape (N, 3), **heliocentric** planet positions [AU] and velocities [AU/day]
  (Sun-relative, as Horizons provides them), in whatever inertial frame the caller uses
  (``ecl_j2000`` in this project).
- ``masses``: shape (N + 1,), [M_sun]; ``masses[0]`` is the central body, ``masses[1 + i]``
  belongs to planet i.
- ``t_out``: times [day] *relative to the epoch of x0*. Any order, any sign: negative times are
  reached by integrating backward from the epoch.
- Output: heliocentric positions (and velocities) with shape (len(t_out), N, 3).

Output times that are not on the step grid (t = k·dt from the epoch) are reached by a partial
step applied to a *copy* of the state, so requested outputs never change the trajectory.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from fork.model.bodies import G


@njit(cache=True)
def plan_outputs(
    t_out: np.ndarray, dt: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Visit order and step bookkeeping for the requested output times.

    Returns (idx_forward, idx_backward, n_steps, remainder):
    ``idx_forward`` are indices with t ≥ 0 sorted by increasing t; ``idx_backward`` those with
    t < 0 sorted by decreasing t (increasing |t|). For each output, ``n_steps`` full steps of
    size dt reach the last grid point before |t|, then a partial step of ``remainder`` days.
    """
    if not dt > 0.0:
        raise ValueError("dt must be positive")
    n = t_out.size
    order = np.argsort(t_out)
    n_fwd = 0
    for i in range(n):
        if t_out[i] >= 0.0:
            n_fwd += 1
    idx_forward = order[n - n_fwd :].copy()
    idx_backward = order[: n - n_fwd][::-1].copy()
    n_steps = np.empty(n, dtype=np.int64)
    remainder = np.empty(n)
    for i in range(n):
        t_abs = abs(t_out[i])
        k = np.int64(np.floor(t_abs / dt))
        # guard against t/dt rounding just below an integer
        if (k + 1) * dt <= t_abs:
            k += 1
        n_steps[i] = k
        remainder[i] = t_abs - k * dt
    return idx_forward, idx_backward, n_steps, remainder


@njit(cache=True)
def helio_to_bary(
    x_h: np.ndarray, v_h: np.ndarray, masses: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric planets (N, 3) → barycentric Sun + planets (N + 1, 3)."""
    n = x_h.shape[0]
    m_tot = masses.sum()
    rc = np.zeros(3)
    vc = np.zeros(3)
    for i in range(n):
        rc += masses[i + 1] * x_h[i]
        vc += masses[i + 1] * v_h[i]
    rc /= m_tot
    vc /= m_tot
    x_b = np.empty((n + 1, 3))
    v_b = np.empty((n + 1, 3))
    x_b[0] = -rc
    v_b[0] = -vc
    for i in range(n):
        x_b[i + 1] = x_h[i] - rc
        v_b[i + 1] = v_h[i] - vc
    return x_b, v_b


@njit(cache=True)
def bary_to_helio_into(
    x_b: np.ndarray, v_b: np.ndarray, out_x: np.ndarray, out_v: np.ndarray
) -> None:
    """Barycentric (N + 1, 3) → heliocentric planets, written into out_x, out_v (N, 3)."""
    for i in range(out_x.shape[0]):
        for k in range(3):
            out_x[i, k] = x_b[i + 1, k] - x_b[0, k]
            out_v[i, k] = v_b[i + 1, k] - v_b[0, k]


@njit(cache=True)
def accel_all_pairs(x: np.ndarray, masses: np.ndarray, out: np.ndarray) -> None:
    """Newtonian accelerations of every body from every other: out[i] = Σ_j G m_j r_ij/|r_ij|³."""
    n = x.shape[0]
    out[:] = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[j, 0] - x[i, 0]
            dy = x[j, 1] - x[i, 1]
            dz = x[j, 2] - x[i, 2]
            r2 = dx * dx + dy * dy + dz * dz
            inv_r3 = G / (r2 * np.sqrt(r2))
            fi = masses[j] * inv_r3
            fj = masses[i] * inv_r3
            out[i, 0] += fi * dx
            out[i, 1] += fi * dy
            out[i, 2] += fi * dz
            out[j, 0] -= fj * dx
            out[j, 1] -= fj * dy
            out[j, 2] -= fj * dz


@njit(cache=True)
def energy_helio(x_h: np.ndarray, v_h: np.ndarray, masses: np.ndarray) -> float:
    """Total (barycentric) energy of Sun + planets from heliocentric planet states."""
    x_b, v_b = helio_to_bary(x_h, v_h, masses)
    n = x_b.shape[0]
    kin = 0.0
    for i in range(n):
        kin += 0.5 * masses[i] * (v_b[i, 0] ** 2 + v_b[i, 1] ** 2 + v_b[i, 2] ** 2)
    pot = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = np.sqrt(
                (x_b[i, 0] - x_b[j, 0]) ** 2
                + (x_b[i, 1] - x_b[j, 1]) ** 2
                + (x_b[i, 2] - x_b[j, 2]) ** 2
            )
            pot -= G * masses[i] * masses[j] / d
    return kin + pot


@njit(cache=True)
def check_inputs(x0: np.ndarray, v0: np.ndarray, masses: np.ndarray) -> None:
    if x0.ndim != 2 or x0.shape[1] != 3 or v0.shape != x0.shape:
        raise ValueError("x0 and v0 must both have shape (N, 3)")
    if masses.size != x0.shape[0] + 1:
        raise ValueError("masses must have length N + 1 (central body first)")
