"""Classical 4th-order Runge–Kutta N-body integrator: the non-symplectic baseline.

Integrates Sun + planets in barycentric coordinates with direct all-pairs Newtonian gravity.
Its energy error drifts secularly; it exists to be compared against the symplectic integrators
(A1 energy panel). Conventions: see ``fork.integrators._nbody``.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from fork.integrators._nbody import (
    accel_all_pairs,
    bary_to_helio_into,
    check_inputs,
    helio_to_bary,
    plan_outputs,
)


@njit(cache=True)
def step(x: np.ndarray, v: np.ndarray, masses: np.ndarray, h: float) -> None:
    """One RK4 step of size h (may be negative) on barycentric x, v (N + 1, 3), in place."""
    a1 = np.empty_like(x)
    a2 = np.empty_like(x)
    a3 = np.empty_like(x)
    a4 = np.empty_like(x)
    accel_all_pairs(x, masses, a1)
    x2 = x + 0.5 * h * v
    v2 = v + 0.5 * h * a1
    accel_all_pairs(x2, masses, a2)
    x3 = x + 0.5 * h * v2
    v3 = v + 0.5 * h * a2
    accel_all_pairs(x3, masses, a3)
    x4 = x + h * v3
    v4 = v + h * a3
    accel_all_pairs(x4, masses, a4)
    x += (h / 6.0) * (v + 2.0 * v2 + 2.0 * v3 + v4)
    v += (h / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4)


@njit(cache=True)
def integrate_states(
    x0: np.ndarray, v0: np.ndarray, masses: np.ndarray, t_out: np.ndarray, dt: float
) -> tuple[np.ndarray, np.ndarray]:
    """Heliocentric (positions, velocities) at t_out, each (len(t_out), N, 3)."""
    check_inputs(x0, v0, masses)
    idx_fwd, idx_bwd, n_steps, rem = plan_outputs(t_out, dt)
    out_x = np.empty((t_out.size, x0.shape[0], 3))
    out_v = np.empty_like(out_x)
    for sign, idx in ((1.0, idx_fwd), (-1.0, idx_bwd)):
        x, v = helio_to_bary(x0, v0, masses)
        k = 0
        for j in idx:
            while k < n_steps[j]:
                step(x, v, masses, sign * dt)
                k += 1
            if rem[j] > 0.0:
                xc, vc = x.copy(), v.copy()
                step(xc, vc, masses, sign * rem[j])
                bary_to_helio_into(xc, vc, out_x[j], out_v[j])
            else:
                bary_to_helio_into(x, v, out_x[j], out_v[j])
    return out_x, out_v


@njit(cache=True)
def integrate(
    x0: np.ndarray, v0: np.ndarray, masses: np.ndarray, t_out: np.ndarray, dt: float
) -> np.ndarray:
    """Heliocentric positions [t_out, body, 3]."""
    return integrate_states(x0, v0, masses, t_out, dt)[0]
