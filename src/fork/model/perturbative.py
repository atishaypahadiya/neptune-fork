"""M8 amendment: first-order perturbative likelihood (PHASE1_SPEC, "M8 amendment", 5fc089f).

Uranus's motion = a reference (Sun with the inner planets folded in, Jupiter, Saturn, Uranus;
no perturber), integrated once, plus its linear response to the perturber's pull:

    δX(t) = Φ(t, t₀) ∫_{t₀}^{t} Φ(s, t₀)⁻¹ B g(s) ds,

with Φ Uranus's 6×6 state-transition matrix about the reference (central finite differences of
the reference integration in Uranus's initial Cartesian state), B the velocity rows, and g the
perturber's acceleration on Uranus in the integrator's heliocentric (inner-barycentre) frame,
direct + indirect: g = G m [(x_P − x_U)/|x_P − x_U|³ − x_P/|x_P|³]. The integral is a cumulative
quadrature over the nodes (``GRID_STEP_DAYS`` from t₀ ∪ every observation epoch ∪ t₀): on each
interval, the cubic through the interval's two nodes and their outer neighbours (shifted inward
at the ends), with weights that depend only on the grid and are precomputed. The mass factors
out of the first-order term.

**The perturber** starts from exactly the full model's initial state
(``fork.model.forward._initial_conditions``), moved to the Sun–Jupiter–Saturn barycentre, and
moves on a Kepler orbit about it with μ = G(M_central + M_J + M_S + m). Jupiter and Saturn act
on it only through that central mass and barycentre; Uranus not at all.

**Second order** (``second_order=True``): one Picard step. The forcing is re-evaluated along
the first-order path x_ref + δx: the perturber's pull there, plus the nonlinear remainder of the
Sun–Jupiter–Saturn field, F(x_ref + δx) − F(x_ref) − A δx (A = ∂F/∂x along the reference). The
remainder is as large as the change in the perturber's pull (both ~10⁻³ of it), so recomputing
only the perturber's pull would not be a second-order step.

**The likelihood**: Uranus at x_ref + δx with the reference velocity (the velocity change only
enters through light time, < 0.01″); exact apparent places; the group offsets marginalised as in
M6; Uranus's four coordinates marginalised analytically about u_ref with J_u = ∂z/∂u of the
Uranus-only model there (``fork.model.linearized.linear_marginal``). u_ref must be near the
posterior (spec: the M6c replacement posterior's median). Solutions outside the box are flagged.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import njit

from fork.inference.fit import uranus_x_to_theta
from fork.inference.priors import PosteriorSpace
from fork.integrators.wisdom_holman import integrate_states
from fork.model.apparent import apparent_ra_dec
from fork.model.bodies import JUPITER, M_CENTRAL_MSUN, SATURN, URANUS, G
from fork.model.forward import ForwardContext, _initial_conditions
from fork.model.linearized import (
    FD_REL_STEP,
    GaussianMarginal,
    fd_jacobian,
    linear_marginal,
    truncation_log10_bound,
)
from fork.model.observation import RaDecObservations

GRID_STEP_DAYS = 60.0
FD_STEP_POS_AU = 1e-5
FD_STEP_VEL_REL = 1e-5  # × |v₀| / |r₀| of Uranus, AU/day (the same relative size as positions)
QUADRATURES = ("cubic", "trapezoid")
_KEPLER_MAX_ITER = 50
_KEPLER_TOL = 1e-14
_M_J = JUPITER.mass_msun
_M_S = SATURN.mass_msun
_M_SJS = M_CENTRAL_MSUN + _M_J + _M_S
_MU_U = G * (M_CENTRAL_MSUN + URANUS.mass_msun)
_W = np.block([[np.zeros((3, 3)), np.eye(3)], [-np.eye(3), np.zeros((3, 3))]])
_B = np.vstack([np.zeros((3, 3)), np.eye(3)])


# ---------------------------------------------------------------------------- the grid


def node_grid(epochs_jd: np.ndarray, epoch_jd: float, step_days: float) -> np.ndarray:
    """Regular nodes at ``step_days`` from t₀ over the epochs' span, ∪ the epochs ∪ t₀."""
    lo, hi = float(np.min(epochs_jd)), float(np.max(epochs_jd))
    k_lo = math.floor((lo - epoch_jd) / step_days)
    k_hi = math.ceil((hi - epoch_jd) / step_days)
    regular = epoch_jd + step_days * np.arange(k_lo, k_hi + 1)
    return np.unique(np.concatenate([regular, epochs_jd, [epoch_jd]]))


def interval_weights(s: np.ndarray, method: str = "cubic") -> tuple[np.ndarray, np.ndarray]:
    """For every interval [s_j, s_{j+1}]: 4 node indices and weights with
    ∫_{s_j}^{s_{j+1}} f ds ≈ Σ_l w_{j,l} f(s_{idx_{j,l}}). ``cubic``: the interpolating cubic on
    nodes j−1 … j+2 (shifted inward at the ends; exact for cubics); ``trapezoid``: (j, j+1)."""
    if method not in QUADRATURES:
        raise ValueError(f"quadrature is one of {QUADRATURES}, not {method!r}")
    n = s.size
    idx = np.zeros((n - 1, 4), dtype=np.int64)
    w = np.zeros((n - 1, 4))
    for j in range(n - 1):
        if method == "trapezoid" or n < 4:
            h = s[j + 1] - s[j]
            idx[j] = (j, j + 1, j, j)
            w[j] = (0.5 * h, 0.5 * h, 0.0, 0.0)
            continue
        first = min(max(j - 1, 0), n - 4)
        st = np.arange(first, first + 4)
        x0, span = s[j], s[st[-1]] - s[st[0]]
        x = (s[st] - x0) / span  # scaled for conditioning
        b = (s[j + 1] - x0) / span
        moments = np.array([b ** (p + 1) / (p + 1) for p in range(4)])
        idx[j] = st
        w[j] = np.linalg.solve(np.vander(x, 4, increasing=True).T, moments) * span
    return idx, w


# ---------------------------------------------------------------------------- numba kernels


@njit(cache=True)
def kepler_positions(r0: np.ndarray, v0: np.ndarray, mu: float, dt: np.ndarray) -> np.ndarray:
    """Positions (n, 3) on the Kepler orbit through (r0, v0) after ``dt`` days (f and g in the
    eccentric-anomaly difference, well defined as e → 0). Elliptic orbits only."""
    r0n = math.sqrt(r0[0] ** 2 + r0[1] ** 2 + r0[2] ** 2)
    alpha = 2.0 / r0n - (v0[0] ** 2 + v0[1] ** 2 + v0[2] ** 2) / mu  # 1/a
    if alpha <= 0.0:
        raise ValueError("the perturber's barycentric orbit is not elliptic")
    a = 1.0 / alpha
    n = math.sqrt(mu * alpha**3)
    c1 = (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) * math.sqrt(alpha / mu)  # e sin E₀
    c2 = 1.0 - r0n * alpha  # e cos E₀
    out = np.empty((dt.size, 3))
    sq = math.sqrt(a**3 / mu)
    for k in range(dt.size):
        m = n * dt[k]
        de = m
        for _ in range(_KEPLER_MAX_ITER):
            sn, cs = math.sin(de), math.cos(de)
            f = de + c1 * (1.0 - cs) - c2 * sn - m
            de -= f / (1.0 + c1 * sn - c2 * cs)
            if abs(f) < _KEPLER_TOL:
                break
        ff = 1.0 - (a / r0n) * (1.0 - math.cos(de))
        gg = dt[k] + sq * (math.sin(de) - de)
        for i in range(3):
            out[k, i] = ff * r0[i] + gg * v0[i]
    return out


@njit(cache=True)
def unit_pull(x_u: np.ndarray, x_p: np.ndarray) -> np.ndarray:
    """The perturber's acceleration on Uranus per unit perturber mass (M☉), heliocentric frame:
    direct + indirect, G[(x_P − x_U)/|x_P − x_U|³ − x_P/|x_P|³]. Both (n, 3)."""
    out = np.empty_like(x_u)
    for k in range(x_u.shape[0]):
        d0, d1, d2 = x_p[k, 0] - x_u[k, 0], x_p[k, 1] - x_u[k, 1], x_p[k, 2] - x_u[k, 2]
        dn3 = (d0 * d0 + d1 * d1 + d2 * d2) ** 1.5
        pn3 = (x_p[k, 0] ** 2 + x_p[k, 1] ** 2 + x_p[k, 2] ** 2) ** 1.5
        out[k, 0] = G * (d0 / dn3 - x_p[k, 0] / pn3)
        out[k, 1] = G * (d1 / dn3 - x_p[k, 1] / pn3)
        out[k, 2] = G * (d2 / dn3 - x_p[k, 2] / pn3)
    return out


@njit(cache=True)
def _response_positions(m_k, phi_pos, accel, idx, w, i0, sel):
    """Position rows of δX(s_k) for k in ``sel``: Φ_pos(s_k) · ∫_{t₀}^{s_k} M(s) a(s) ds, with
    M = Φ⁻¹B (n, 6, 3), Φ_pos (n, 3, 6), a (n, 3); signed integral for s_k < t₀."""
    n = accel.shape[0]
    h = np.zeros((n, 6))
    for k in range(n):
        for i in range(6):
            h[k, i] = (m_k[k, i, 0] * accel[k, 0] + m_k[k, i, 1] * accel[k, 1]
                       + m_k[k, i, 2] * accel[k, 2])  # fmt: skip
    c = np.zeros((n, 6))
    for j in range(i0, n - 1):  # forward from t₀
        for i in range(6):
            s = 0.0
            for q in range(4):
                s += w[j, q] * h[idx[j, q], i]
            c[j + 1, i] = c[j, i] + s
    for j in range(i0 - 1, -1, -1):  # backward: ∫_{t₀}^{s_j} = −∫_{s_j}^{t₀}
        for i in range(6):
            s = 0.0
            for q in range(4):
                s += w[j, q] * h[idx[j, q], i]
            c[j, i] = c[j + 1, i] - s
    out = np.empty((sel.size, 3))
    for r in range(sel.size):
        k = sel[r]
        for i in range(3):
            s = 0.0
            for q in range(6):
                s += phi_pos[k, i, q] * c[k, q]
            out[r, i] = s
    return out


def _sjs_accel(x: np.ndarray, xj: np.ndarray, xs: np.ndarray) -> np.ndarray:
    """Heliocentric acceleration on Uranus from the Sun (inner planets folded in), Jupiter and
    Saturn, without the x-independent indirect terms (they cancel in differences)."""
    r = np.linalg.norm(x, axis=1)[:, None]
    out = -_MU_U * x / r**3
    for xp, mp in ((xj, _M_J), (xs, _M_S)):
        d = xp - x
        out = out + G * mp * d / np.linalg.norm(d, axis=1)[:, None] ** 3
    return out


def _sjs_jacobian(x: np.ndarray, xj: np.ndarray, xs: np.ndarray) -> np.ndarray:
    """A = ∂(``_sjs_accel``)/∂x along x, (n, 3, 3)."""
    eye = np.eye(3)[None]
    r = np.linalg.norm(x, axis=1)[:, None, None]
    out = -_MU_U * (eye / r**3 - 3.0 * np.einsum("ki,kj->kij", x, x) / r**5)
    for xp, mp in ((xj, _M_J), (xs, _M_S)):
        d = xp - x
        dn = np.linalg.norm(d, axis=1)[:, None, None]
        out = out - G * mp * (eye / dn**3 - 3.0 * np.einsum("ki,kj->kij", d, d) / dn**5)
    return out


def symplectic_rel_dev(phi: np.ndarray, n_mean: float) -> float:
    """max_k |Φ̃ᵀWΦ̃ − W| / max(1, ‖Φ̃‖²) with Φ̃ = SΦS⁻¹, S = diag(I, I/n): Φ in units where
    positions and velocities × 1/n are commensurate (Φ̃ is symplectic iff Φ is)."""
    s = np.concatenate([np.ones(3), np.full(3, 1.0 / n_mean)])
    pt = phi * s[None, :, None] / s[None, None, :]
    dev = np.abs(np.einsum("kji,jl,klm->kim", pt, _W, pt) - _W[None]).max(axis=(1, 2))
    scale = np.maximum(1.0, np.linalg.norm(pt, axis=(1, 2), ord=2) ** 2)
    return float(np.max(dev / scale))


# ---------------------------------------------------------------------------- the problem


@dataclass(frozen=True)
class PerturbativeResult:
    marginal: GaussianMarginal
    outside_box: bool
    truncation_log10: float  # the union bound on the box mass outside, log10
    delta_x_obs_au: np.ndarray  # (n_epoch, 3) Uranus's position change at the epochs


@dataclass(frozen=True, eq=False)
class PerturbativeProblem:
    """Everything precomputed once per data set, error model, reference and grid."""

    obs: RaDecObservations
    ctx: ForwardContext
    space: PosteriorSpace
    u_ref: np.ndarray  # (4,)
    nodes_jd: np.ndarray  # (n,)
    i0: int  # index of t₀ in the nodes
    obs_node: np.ndarray  # (n_epoch,) node index of each observation epoch
    x_ref: np.ndarray  # (n, 3) reference Uranus, heliocentric about the inner barycentre
    v_ref: np.ndarray
    x_j: np.ndarray  # (n, 3) reference Jupiter and Saturn
    x_s: np.ndarray
    phi_pos: np.ndarray  # (n, 3, 6) position rows of Φ(s_k, t₀)
    m_k: np.ndarray  # (n, 6, 3) Φ(s_k, t₀)⁻¹ B
    q_idx: np.ndarray  # (n−1, 4) quadrature stencils
    q_w: np.ndarray  # (n−1, 4) quadrature weights
    a_sjs: np.ndarray  # (n, 3, 3) ∂F_SJS/∂x along the reference
    jac_u: np.ndarray  # (m_aug, 4) ∂z/∂u, Uranus-only model at u_ref
    sigma_raw: np.ndarray  # (4,)
    symplectic_dev: float  # ``symplectic_rel_dev`` over the nodes
    quadrature: str
    step_days: float

    @classmethod
    def build(cls, obs: RaDecObservations, ctx: ForwardContext, space: PosteriorSpace,
              u_ref: np.ndarray, *, step_days: float = GRID_STEP_DAYS,
              quadrature: str = "cubic") -> PerturbativeProblem:  # fmt: skip
        from fork.model.student_t import StudentTObservations

        if isinstance(obs, StudentTObservations):
            raise TypeError("Student-t observations are out of scope for analytic M8")
        box = space.uranus_box
        if box is None:
            raise ValueError("the space needs its Uranus box")
        u_ref = np.asarray(u_ref, dtype=float)
        nodes = node_grid(obs.geometry.jd_tt, ctx.epoch_jd_tdb, step_days)
        i0 = int(np.flatnonzero(nodes == ctx.epoch_jd_tdb)[0])
        obs_node = np.searchsorted(nodes, obs.geometry.jd_tt)
        theta_u = uranus_x_to_theta(np.concatenate([u_ref, space.pq]))
        x0, v0, masses = _initial_conditions(theta_u, np.empty(0), ctx.x_js_ecl_j2000,
                                             ctx.v_js_ecl_j2000)  # fmt: skip
        t_rel = np.ascontiguousarray(nodes - ctx.epoch_jd_tdb)

        def run(x, v):
            return integrate_states(x, v, masses, t_rel, ctx.dt_days, ctx.scheme)

        pos, vel = run(x0, v0)
        n_mean = float(np.linalg.norm(v0[2]) / np.linalg.norm(x0[2]))
        steps = np.concatenate([np.full(3, FD_STEP_POS_AU), np.full(3, FD_STEP_VEL_REL * n_mean)])
        phi = np.empty((nodes.size, 6, 6))
        for c in range(6):
            xp, vp, xm, vm = x0.copy(), v0.copy(), x0.copy(), v0.copy()
            (xp if c < 3 else vp)[2, c % 3] += steps[c]
            (xm if c < 3 else vm)[2, c % 3] -= steps[c]
            (pp, vpp), (pm, vpm) = run(xp, vp), run(xm, vm)
            phi[:, :3, c] = (pp[:, 2] - pm[:, 2]) / (2.0 * steps[c])
            phi[:, 3:, c] = (vpp[:, 2] - vpm[:, 2]) / (2.0 * steps[c])
        phi_inv = -np.einsum("ij,kjl,lm->kim", _W, np.transpose(phi, (0, 2, 1)), _W)  # −WΦᵀW
        m_k = np.ascontiguousarray(phi_inv @ _B)
        q_idx, q_w = interval_weights(nodes, quadrature)
        sigma_raw = np.asarray(box.sigma_scaled) / math.sqrt(max(box.chi2_per_dof, 1.0))

        def z_u(u):
            return obs.whitened(uranus_x_to_theta(np.concatenate([u, space.pq])), ctx)

        jac_u = fd_jacobian(z_u, u_ref, FD_REL_STEP * sigma_raw, "central")
        x_ref, x_j, x_s = (np.ascontiguousarray(pos[:, i]) for i in (2, 0, 1))
        return cls(obs, ctx, space, u_ref, nodes, i0, obs_node, x_ref,
                   np.ascontiguousarray(vel[:, 2]), x_j, x_s,
                   np.ascontiguousarray(phi[:, :3, :]), m_k, q_idx, q_w,
                   _sjs_jacobian(x_ref, x_j, x_s), jac_u, sigma_raw,
                   symplectic_rel_dev(phi, n_mean), quadrature, step_days)  # fmt: skip

    # ------------------------------------------------------------------ pieces

    @property
    def log_volume(self) -> float:
        return self.space.uranus_box.log_volume

    def theta(self, phi: np.ndarray) -> np.ndarray:
        return self.space.theta(np.concatenate([self.u_ref, np.asarray(phi, dtype=float)]))

    def perturber_initial_state(self, phi: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        """(r0, v0) heliocentric at t₀ exactly as the full model builds them, and m in M☉."""
        th = self.theta(phi)
        x0, v0, masses = _initial_conditions(th[:6], th[6:], self.ctx.x_js_ecl_j2000,
                                             self.ctx.v_js_ecl_j2000)  # fmt: skip
        return x0[3], v0[3], float(masses[4])

    def perturber_positions(self, phi: np.ndarray) -> np.ndarray:
        """Heliocentric (inner-barycentre) positions of the Keplerian perturber at the nodes."""
        r0, v0, m = self.perturber_initial_state(phi)
        xj0, xs0 = self.ctx.x_js_ecl_j2000
        vj0, vs0 = self.ctx.v_js_ecl_j2000
        xb0 = (_M_J * xj0 + _M_S * xs0) / _M_SJS
        vb0 = (_M_J * vj0 + _M_S * vs0) / _M_SJS
        rel = kepler_positions(r0 - xb0, v0 - vb0, G * (_M_SJS + m),
                               self.nodes_jd - self.ctx.epoch_jd_tdb)  # fmt: skip
        return rel + (_M_J * self.x_j + _M_S * self.x_s) / _M_SJS

    def integrated_states(self, phi: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The full integrator with the perturber, from u_ref, at the nodes: the perturber's
        positions (Sun, J, S and Uranus all acting on it) and Jupiter's and Saturn's position
        changes from the reference. Diagnostics only (E_K and E_JS)."""
        th = self.theta(phi)
        x0, v0, masses = _initial_conditions(th[:6], th[6:], self.ctx.x_js_ecl_j2000,
                                             self.ctx.v_js_ecl_j2000)  # fmt: skip
        pos, _ = integrate_states(x0, v0, masses,
                                  np.ascontiguousarray(self.nodes_jd - self.ctx.epoch_jd_tdb),
                                  self.ctx.dt_days, self.ctx.scheme)  # fmt: skip
        return (np.ascontiguousarray(pos[:, 3]), pos[:, 0] - self.x_j, pos[:, 1] - self.x_s)

    def js_response_accel(self, dx_j: np.ndarray, dx_s: np.ndarray) -> np.ndarray:
        """The change in Jupiter's and Saturn's pull on the reference Uranus (direct and
        indirect) when they move by (dx_j, dx_s): the forcing the model leaves out (E_JS)."""

        def pull(xj, xs):
            out = np.zeros_like(self.x_ref)
            for xp, mp in ((xj, _M_J), (xs, _M_S)):
                d = xp - self.x_ref
                out += G * mp * (d / np.linalg.norm(d, axis=1)[:, None] ** 3
                                 - xp / np.linalg.norm(xp, axis=1)[:, None] ** 3)  # fmt: skip
            return out

        return pull(self.x_j + dx_j, self.x_s + dx_s) - pull(self.x_j, self.x_s)

    def _response(self, accel: np.ndarray, sel: np.ndarray) -> np.ndarray:
        return _response_positions(self.m_k, self.phi_pos, np.ascontiguousarray(accel),
                                   self.q_idx, self.q_w, self.i0, sel)  # fmt: skip

    def delta_x(self, phi: np.ndarray, *, second_order: bool = False,
                perturber_pos: np.ndarray | None = None,
                extra_accel: np.ndarray | None = None,
                at: np.ndarray | None = None) -> np.ndarray:  # fmt: skip
        """Uranus's position change (len(at), 3) at the nodes ``at`` (default: the observation
        epochs). Diagnostics: ``perturber_pos`` replaces the Keplerian perturber (E_K);
        ``extra_accel`` (n, 3) is added to the forcing in every pass (E_JS)."""
        at = self.obs_node if at is None else np.asarray(at, dtype=np.int64)
        x_p = self.perturber_positions(phi) if perturber_pos is None else perturber_pos
        _, _, m = self.perturber_initial_state(phi)
        extra = 0.0 if extra_accel is None else extra_accel
        pull = m * unit_pull(self.x_ref, x_p) + extra
        if not second_order:
            return self._response(pull, at)
        dx = self._response(pull, np.arange(self.nodes_jd.size))
        x1 = self.x_ref + dx
        rem = (_sjs_accel(x1, self.x_j, self.x_s) - _sjs_accel(self.x_ref, self.x_j, self.x_s)
               - np.einsum("kij,kj->ki", self.a_sjs, dx))  # fmt: skip
        return self._response(m * unit_pull(x1, x_p) + rem + extra, at)

    def z_from_delta(self, dx_obs: np.ndarray) -> np.ndarray:
        k = self.obs_node
        ra, dec = apparent_ra_dec(self.x_ref[k] + dx_obs, self.v_ref[k], self.obs.geometry)
        return self.obs.whitened_from_residuals(self.obs.residuals_arcsec(ra, dec))

    def marginal_from_z(self, z: np.ndarray) -> GaussianMarginal:
        return linear_marginal(z, self.jac_u, self.u_ref, self.log_volume)

    # ------------------------------------------------------------------ the likelihood

    def evaluate(self, phi: np.ndarray, *, second_order: bool = False,
                 perturber_pos: np.ndarray | None = None,
                 extra_accel: np.ndarray | None = None) -> PerturbativeResult:  # fmt: skip
        dx = self.delta_x(phi, second_order=second_order, perturber_pos=perturber_pos,
                          extra_accel=extra_accel)  # fmt: skip
        m = self.marginal_from_z(self.z_from_delta(dx))
        box = self.space.uranus_box
        outside = bool(np.any(m.u_star < box.lo) or np.any(m.u_star > box.hi))
        return PerturbativeResult(m, outside,
                                  truncation_log10_bound(m.u_star, m.cov_u, box.lo, box.hi),
                                  dx)  # fmt: skip

    def log_likelihood(self, phi: np.ndarray, *, second_order: bool = False) -> float:
        """ln L(φ) with u and the offsets marginalised: comparable with
        ``LinearizedProblem.log_likelihood`` and M6c's ln Z."""
        try:
            dx = self.delta_x(phi, second_order=second_order)
            v = self.marginal_from_z(self.z_from_delta(dx)).m2lnl
        except (ValueError, np.linalg.LinAlgError, ZeroDivisionError):
            return -math.inf
        return -0.5 * v if math.isfinite(v) else -math.inf

    def box_sigma_offset(self, u: np.ndarray) -> np.ndarray:
        box = self.space.uranus_box
        return (np.asarray(u) - np.asarray(box.centre)) / np.asarray(box.sigma_scaled)


def posterior_median_u(q: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted per-coordinate median of Uranus's four q entries (the spec's u_ref)."""
    w = np.asarray(weights, dtype=float) / np.sum(weights)
    out = np.empty(4)
    for i in range(4):
        o = np.argsort(q[:, i])
        out[i] = q[o, i][np.searchsorted(np.cumsum(w[o]), 0.5)]
    return out


def posterior_std_u(q: np.ndarray, weights: np.ndarray) -> np.ndarray:
    w = np.asarray(weights, dtype=float) / np.sum(weights)
    mu = w @ q[:, :4]
    return np.sqrt(w @ (q[:, :4] - mu) ** 2)


def reference_from_m6c(variant: str) -> tuple[np.ndarray, np.ndarray]:
    """(median u, std u) of the variant's M6c posterior: the replacement posterior (R2 as cell 2)
    for the headline, its own four cells for a sensitivity variant. Data-only quantities."""
    from fork.export import m6c, m6c_repeat
    from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult

    if variant == "headline":
        cells = tuple(m6c_repeat.replacement_cells())
    else:
        cells = tuple(NestedResult.load(m6c.Job(variant, k).path(".pkl"))
                      for k in range(len(m6c.variant_cells(variant))))  # fmt: skip
    q, w, _ = PartitionedResult(cells, err_scale=ERR_INFLATION).samples_and_weights()
    return posterior_median_u(q, w), posterior_std_u(q, w)
