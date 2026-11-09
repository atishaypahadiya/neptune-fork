"""M8: Uranus's four in-plane coordinates marginalised analytically (PHASE1_SPEC M8).

Plan and its reasons: docs/notes/M8.md. In short:

- u = q[:4], Uranus's (a, e·cos ϖ, e·sin ϖ, λ₀) in ecl_j2000, the coordinates of the M6c box
  (``fork.inference.priors.UranusBox``); i and Ω stay held (docs/RULES.md rule 2). φ = q[4:],
  the perturber's 7.
- The group offsets are already marginalised exactly, for any u, by
  ``RaDecObservations.whitened``: its output z(u, φ) has ‖z‖² = −2 ln L up to a constant of
  the error model. Linearising z in u, z ≈ z₀ + J (u − u₀), and integrating against the box
  prior 1/V_box over all of R⁴:

      −2 ln L(φ) = χ²_min(φ) + ln det(JᵀJ) − 4 ln 2π + 2 ln V_box

  With the same constant dropped as in M6c, ∫ L(φ) π(φ) dφ estimates M6c's Z.
- The box truncation the formula ignores is bounded by ``truncation_bound``.
- Mode (a): J₀ from the Uranus-only model at the box centre (the Uranus-only fit), the residual
  at the centre: 1 integration per call. Mode (b): one more Gauss–Newton step, re-integrating
  at mode (a)'s u* and re-evaluating J there by forward differences: 6 integrations per call.
- Jacobians are finite differences (docs/notes/M8.md says why not variational equations).

``profile_min`` and ``laplace_m2lnl`` are the A8a references: the full likelihood minimised
over u, and its Laplace marginal with the exact Hessian (Gauss–Newton term plus Σ z_k ∇²z_k).
Student-t observations (M6c S3) are refused: their marginal over u is not Gaussian.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp
from scipy.stats import norm

from fork.inference.fit import uranus_x_to_theta
from fork.inference.priors import PosteriorSpace
from fork.model.forward import ForwardContext
from fork.model.observation import RaDecObservations

N_U = 4
FD_REL_STEP = 0.01  # Jacobian step, in raw σ of each Uranus coordinate
HESS_REL_STEP = 0.5  # second-derivative step for the Laplace reference, in raw σ
_LN_2PI = math.log(2.0 * math.pi)

ZFn = Callable[[np.ndarray], np.ndarray]


class _Counted:
    """A z(u) that counts its calls (one integration each)."""

    def __init__(self, z_fn: ZFn) -> None:
        self.z_fn, self.n = z_fn, 0

    def __call__(self, u: np.ndarray) -> np.ndarray:
        self.n += 1
        return self.z_fn(u)


@dataclass(frozen=True)
class GaussianMarginal:
    """The linear-Gaussian marginal over u about one linearisation point."""

    m2lnl: float  # −2 ln ∫ exp(−½‖z‖²) du / V_box (the error model's constant dropped)
    chi2_min: float  # ‖z₀ + J Δu*‖², the linearised profile
    u_star: np.ndarray  # (n_u,)
    cov_u: np.ndarray  # (JᵀJ)⁻¹, the conditional covariance of u given φ
    logdet_jtj: float
    n_integrations: int = 0


def fd_jacobian(z_fn: ZFn, u: np.ndarray, step: np.ndarray, scheme: str = "central",
                z0: np.ndarray | None = None) -> np.ndarray:  # fmt: skip
    """∂z/∂u by finite differences: "central" (2 n_u calls) or "forward" (n_u calls, plus
    one for z₀ unless given)."""
    u = np.asarray(u, dtype=float)
    cols = []
    if scheme == "forward" and z0 is None:
        z0 = z_fn(u)
    for i, h in enumerate(np.asarray(step, dtype=float)):
        du = np.zeros_like(u)
        du[i] = h
        if scheme == "central":
            cols.append((z_fn(u + du) - z_fn(u - du)) / (2.0 * h))
        elif scheme == "forward":
            cols.append((z_fn(u + du) - z0) / h)
        else:
            raise ValueError(f"unknown scheme {scheme!r}")
    return np.array(cols).T


def linear_marginal(z0: np.ndarray, jac: np.ndarray, u0: np.ndarray,
                    log_volume: float) -> GaussianMarginal:  # fmt: skip
    """Exact marginal of exp(−½‖z₀ + J(u − u₀)‖²) / V over all of R^n (columns of J are
    rescaled to unit norm for the solve; the determinant is corrected for it)."""
    c = np.linalg.norm(jac, axis=0)
    js = jac / c
    du_s = np.linalg.lstsq(js, -z0, rcond=None)[0]
    du = du_s / c
    chi2 = float(np.sum((z0 + jac @ du) ** 2))
    jtj_s = js.T @ js
    sign, logdet_s = np.linalg.slogdet(jtj_s)
    if sign <= 0:
        raise np.linalg.LinAlgError("JᵀJ is not positive definite")
    logdet = float(logdet_s + 2.0 * np.sum(np.log(c)))
    cov = np.linalg.inv(jtj_s) / np.outer(c, c)
    n = jac.shape[1]
    return GaussianMarginal(chi2 + logdet - n * _LN_2PI + 2.0 * log_volume, chi2,
                            np.asarray(u0, dtype=float) + du, cov, logdet)  # fmt: skip


def gauss_newton_marginal(z_fn: ZFn, u0: np.ndarray, jac0: np.ndarray, log_volume: float, *,
                          refine: bool, step: np.ndarray) -> GaussianMarginal:  # fmt: skip
    """Mode (a) (``refine=False``): one evaluation at u₀ with the fixed J₀. Mode (b): then one
    more at u₁ = u*(a), J re-evaluated there by forward differences, and a second step."""
    z = _Counted(z_fn)
    m = linear_marginal(z(u0), jac0, u0, log_volume)
    if refine:
        u1 = m.u_star
        z1 = z(u1)
        m = linear_marginal(z1, fd_jacobian(z, u1, step, "forward", z0=z1), u1, log_volume)
    return replace(m, n_integrations=z.n)


def profile_min(z_fn: ZFn, u0: np.ndarray, scale: np.ndarray,
                rel_step: float = FD_REL_STEP) -> tuple[np.ndarray, float, int]:  # fmt: skip
    """(u_min, χ²_min, integrations): Levenberg–Marquardt on ‖z(u)‖² over unconstrained u,
    in coordinates u = u₀ + scale·s, central-difference Jacobian at ``rel_step`` in s."""
    z = _Counted(z_fn)
    u0, scale = np.asarray(u0, dtype=float), np.asarray(scale, dtype=float)

    def f(s):
        return z(u0 + scale * s)

    def jac(s):
        return fd_jacobian(f, s, np.full(s.size, rel_step), "central")

    sol = least_squares(f, np.zeros(u0.size), jac=jac, method="lm", xtol=1e-14, ftol=1e-14,
                        gtol=1e-14, max_nfev=200)  # fmt: skip
    u_min = u0 + scale * sol.x
    return u_min, float(np.sum(z(u_min) ** 2)), z.n


def laplace_m2lnl(z_fn: ZFn, u_min: np.ndarray, step: np.ndarray, log_volume: float) -> float:
    """Laplace −2 ln L about the minimum, with the exact Hessian of ½‖z‖²:
    JᵀJ + Σ_k z_k ∇²z_k, every derivative by central differences at ``step``."""
    u = np.asarray(u_min, dtype=float)
    h = np.asarray(step, dtype=float)
    n = u.size
    z0 = z_fn(u)
    e = np.eye(n) * h
    zp = [z_fn(u + e[i]) for i in range(n)]
    zm = [z_fn(u - e[i]) for i in range(n)]
    jac = np.array([(zp[i] - zm[i]) / (2.0 * h[i]) for i in range(n)]).T
    d2 = np.empty((n, n, z0.size))
    for i in range(n):
        d2[i, i] = (zp[i] - 2.0 * z0 + zm[i]) / h[i] ** 2
        for j in range(i + 1, n):
            v = (z_fn(u + e[i] + e[j]) - z_fn(u + e[i] - e[j]) - z_fn(u - e[i] + e[j])
                 + z_fn(u - e[i] - e[j])) / (4.0 * h[i] * h[j])  # fmt: skip
            d2[i, j] = d2[j, i] = v
    a = jac.T @ jac + np.einsum("ijk,k->ij", d2, z0)
    a_s = a * np.outer(h, h)  # in step units, for conditioning
    sign, logdet_s = np.linalg.slogdet(a_s)
    if sign <= 0:
        return math.nan
    logdet = float(logdet_s - 2.0 * np.sum(np.log(h)))
    return float(np.sum(z0**2)) + logdet - n * _LN_2PI + 2.0 * log_volume


def _log_face_masses(u_star, cov, lo, hi) -> np.ndarray:
    sd = np.sqrt(np.diag(cov))
    u = np.asarray(u_star, dtype=float)
    return np.concatenate([norm.logsf((np.asarray(hi) - u) / sd),
                           norm.logsf((u - np.asarray(lo)) / sd)])  # fmt: skip


def truncation_log10_bound(u_star, cov, lo, hi) -> float:
    """log10 of the union bound on the mass of N(u*, cov) outside the box [lo, hi]."""
    return float(logsumexp(_log_face_masses(u_star, cov, lo, hi)) / math.log(10.0))


def truncation_bound(u_star, cov, lo, hi) -> tuple[float, float]:
    """(P_out bound, bound on the error in −2 ln L): the union bound over the box's faces,
    Σ P(u_i beyond a face) with each coordinate's marginal σ, which is ≥ the true outside
    mass; ignoring the box overstates the integral by 1/(1 − P_out) at most."""
    p = float(min(1.0, math.exp(logsumexp(_log_face_masses(u_star, cov, lo, hi)))))
    return p, (math.inf if p >= 1.0 else -2.0 * math.log1p(-p))


@dataclass(frozen=True)
class LinearizedProblem:
    """The M6c likelihood with u marginalised, for one data set, context and prior space."""

    obs: RaDecObservations
    ctx: ForwardContext
    space: PosteriorSpace  # gives θ(q) and the Uranus box; cell slices do not matter
    sigma_raw: np.ndarray  # the Uranus-only fit's formal σ (box σ / √(χ²/dof))
    jac_ref: np.ndarray  # J₀ = ∂z/∂u, Uranus-only model, at the box centre

    @classmethod
    def build(cls, obs: RaDecObservations, ctx: ForwardContext,
              space: PosteriorSpace) -> LinearizedProblem:  # fmt: skip
        from fork.model.student_t import StudentTObservations

        if isinstance(obs, StudentTObservations):
            raise TypeError("Student-t observations have no analytic marginal over u "
                            "(docs/notes/M8.md, issue 2)")  # fmt: skip
        box = space.uranus_box
        if box is None:
            raise ValueError("the space needs its Uranus box")
        sigma_raw = np.asarray(box.sigma_scaled) / math.sqrt(max(box.chi2_per_dof, 1.0))
        u_ref = np.asarray(box.centre, dtype=float)

        def z_uranus_only(u):
            return obs.whitened(uranus_x_to_theta(np.concatenate([u, space.pq])), ctx)

        jac = fd_jacobian(z_uranus_only, u_ref, FD_REL_STEP * sigma_raw, "central")
        return cls(obs, ctx, space, sigma_raw, jac)

    @property
    def u_ref(self) -> np.ndarray:
        return np.asarray(self.space.uranus_box.centre, dtype=float)

    @property
    def log_volume(self) -> float:
        return self.space.uranus_box.log_volume

    def z_fn(self, phi: np.ndarray) -> ZFn:
        phi = np.asarray(phi, dtype=float)

        def z(u):
            return self.obs.whitened(self.space.theta(np.concatenate([u, phi])), self.ctx)

        return z

    def marginal(self, phi: np.ndarray, mode: str = "a") -> GaussianMarginal:
        if mode not in ("a", "b"):
            raise ValueError(f"mode is 'a' or 'b', not {mode!r}")
        return gauss_newton_marginal(self.z_fn(phi), self.u_ref, self.jac_ref, self.log_volume,
                                     refine=mode == "b",
                                     step=FD_REL_STEP * self.sigma_raw)  # fmt: skip

    def log_likelihood(self, phi: np.ndarray, mode: str = "a") -> float:
        try:
            v = self.marginal(phi, mode).m2lnl
        except np.linalg.LinAlgError:
            return -math.inf
        return -0.5 * v if math.isfinite(v) else -math.inf

    def profile(self, phi: np.ndarray, u_start: np.ndarray) -> tuple[np.ndarray, float, int]:
        return profile_min(self.z_fn(phi), u_start, self.sigma_raw)

    def laplace(self, phi: np.ndarray, u_min: np.ndarray) -> float:
        return laplace_m2lnl(self.z_fn(phi), u_min, HESS_REL_STEP * self.sigma_raw,
                             self.log_volume)  # fmt: skip

    def truncation(self, m: GaussianMarginal) -> tuple[float, float]:
        box = self.space.uranus_box
        return truncation_bound(m.u_star, m.cov_u, box.lo, box.hi)

    def box_sigma_offset(self, u: np.ndarray) -> np.ndarray:
        """(u − centre) in box-σ units, per coordinate (the spec's 40 box-σ record)."""
        box = self.space.uranus_box
        return (np.asarray(u) - np.asarray(box.centre)) / np.asarray(box.sigma_scaled)
