"""The hidden-mass prior, shared by every sampler (emcee in ``mcmc``, dynesty in ``nested``).

Sampled vector q (project rule 2: Uranus's in-plane four always free, i and Ω held):

    Uranus    (a [AU], e·cos ϖ, e·sin ϖ, λ₀ [rad])
    perturber (log10 m [M⊕], ln a [AU], √e·cos ϖ, √e·sin ϖ, λ₀ [rad])  — any may be held

Prior (PHASE1_SPEC M5, "Priors"), a product of independent pieces:

- perturber log10 m uniform on [0, 2] (1–100 M⊕);
- perturber a log-uniform on [20, 60] AU, i.e. ln a uniform (q holds ln a, so a flat density
  in q is exactly log-uniform in a);
- perturber e uniform on [0, e_max] and ϖ uniform on [0, 2π). In q this is (k, h) = √e·(cos ϖ,
  sin ϖ) uniform on the disk k² + h² ≤ e_max (area element k dk dϖ ∝ de dϖ), so the density
  is again flat in q;
- perturber λ₀ uniform on [0, 2π) (periodic);
- Uranus (a, k, h, λ₀) uniform in a box (``UranusBox``) centred on the Uranus-only fit, with
  half-widths ``URANUS_BOX_N_SIGMA`` = 50 times that fit's formal σ scaled by √(χ²/dof).
  The scaling matters: the Uranus-only fit leaves the perturber's signal in its residuals
  (χ²/dof ≈ 23 on seed 0), so its raw covariance is ~4.8× too tight, and the perturbed
  posteriors sit 40–63 raw σ from it (seed 0), which an unscaled 50σ box would truncate.
  Scaled, they sit ≲ 13 box-σ from the centre. The box is common to every Phase 2 model
  fitted to the same data, so its log-volume largely cancels in Bayes factors.

The density is flat in q on its support, so ``log_prior`` returns 0 or −∞ (unnormalised; emcee
needs nothing more), and ``prior_transform`` maps the unit cube onto exactly that support for
nested sampling, whose evidence is normalised by construction.

**Cells** (partitioned nested sampling): ``PosteriorSpace.a_au_range`` restricts a to
[a_lo, a_hi] ⊂ [20, 60] AU and ``e_range`` restricts e to [e_lo, e_hi] ⊂ [0, e_max], with
the prior renormalised on the cell. Both pieces are independent, so the cell's share of the
full prior is w = [ln(a_hi/a_lo) / ln(60/20)] · [(e_hi − e_lo) / e_max] (``slice_weight``),
and the full evidence is Z = Σ_k w_k Z_k. Cell edges change only how hard each sub-run is,
not the result.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from fork.data.synthetic import LongitudeObservations
from fork.inference.fit import E_PRIOR_MAX, TWO_PI, uranus_x_to_theta, wrap_rad
from fork.model.forward import ForwardContext, predict_uranus_lon

PRIOR_LOG10_M_EARTH = (0.0, 2.0)
PRIOR_A_AU = (20.0, 60.0)
URANUS_BOX_N_SIGMA = 50.0
_URANUS_MAX_E = 0.95  # physical guard, only used when no Uranus box is set
PERTURBER_Q_NAMES = ("p_log10_m_earth", "p_ln_a_au", "p_sqrt_e_cos_varpi", "p_sqrt_e_sin_varpi",
                     "p_lam0_ecl_j2000_rad")  # fmt: skip
URANUS_Q_NAMES = ("u_a_au", "u_e_cos_varpi", "u_e_sin_varpi", "u_lam0_ecl_j2000_rad")
# M6 inclined perturber: (sin i cos Ω, sin i sin Ω), i and Ω relative to Uranus's orbital plane
INCLINATION_Q_NAMES = ("p_sin_i_cos_node_rel_uranus", "p_sin_i_sin_node_rel_uranus")
INC_PRIOR_MAX_DEG = 10.0
# finite-difference steps for the Uranus-only Jacobian, in q units (a, k, h, λ₀)
_FD_STEP = (1e-6, 1e-7, 1e-7, 1e-8)


@dataclass(frozen=True)
class UranusBox:
    """Uniform box for Uranus's (a, e·cos ϖ, e·sin ϖ, λ₀)."""

    centre: tuple[float, float, float, float]
    half_width: tuple[float, float, float, float]
    sigma_scaled: tuple[float, float, float, float] = (math.nan,) * 4  # box-σ (for reporting)
    chi2_per_dof: float = math.nan  # of the Uranus-only fit the box was built from

    @property
    def lo(self) -> np.ndarray:
        return np.asarray(self.centre) - np.asarray(self.half_width)

    @property
    def hi(self) -> np.ndarray:
        return np.asarray(self.centre) + np.asarray(self.half_width)

    @property
    def log_volume(self) -> float:
        return float(np.sum(np.log(2.0 * np.asarray(self.half_width))))

    def contains(self, q_uranus: np.ndarray) -> bool:
        return bool(np.all((self.lo <= q_uranus) & (q_uranus <= self.hi)))

    @classmethod
    def from_uranus_only_fit(cls, obs: LongitudeObservations, context: ForwardContext,
                             inc_raan_ecl_j2000_rad: tuple[float, float],
                             n_sigma: float = URANUS_BOX_N_SIGMA) -> UranusBox:  # fmt: skip
        """Centre: the Uranus-only LM fit. σ: from its finite-difference Jacobian, scaled by
        √(χ²/dof) (see module docstring). Half-widths: ``n_sigma`` scaled σ."""
        from fork.inference.fit import fit_uranus_elements, uranus_theta_to_x

        fit = fit_uranus_elements(obs, context, inc_raan_ecl_j2000_rad)
        x = uranus_theta_to_x(fit.theta)

        def resid(x4: np.ndarray) -> np.ndarray:
            theta = uranus_x_to_theta(np.concatenate([x4, x[4:]]))
            pred = predict_uranus_lon(theta, obs.jd_tdb, context)
            return wrap_rad(obs.lon_ecl_j2000_rad - pred) / obs.sigma_rad

        cols = []
        for i, h in enumerate(_FD_STEP):
            dx = np.zeros(4)
            dx[i] = h
            cols.append((resid(x[:4] + dx) - resid(x[:4] - dx)) / (2.0 * h))
        jac = np.array(cols).T
        sigma = np.sqrt(np.diag(np.linalg.inv(jac.T @ jac)))
        chi2_nu = fit.chi2 / fit.dof
        s = sigma * math.sqrt(max(chi2_nu, 1.0))
        return cls(tuple(float(v) for v in x[:4]), tuple(float(v) for v in n_sigma * s),
                   tuple(float(v) for v in s), float(chi2_nu))  # fmt: skip


@dataclass(frozen=True)
class PosteriorSpace:
    """Which perturber entries of q are sampled (held ones carry their value, in q units), the
    eccentricity bound, Uranus's held (p, q) = tan(i/2)(cos Ω, sin Ω), and its prior box."""

    held: tuple[float | None, ...] = (None, None, None, None, None)
    e_max: float = E_PRIOR_MAX
    pq: tuple[float, float] = (0.0, 0.0)
    uranus_box: UranusBox | None = None  # required for prior_transform (nested sampling)
    a_au_range: tuple[float, float] = PRIOR_A_AU  # a slice of the a prior (renormalised on it)
    e_range: tuple[float, float] | None = None  # e slice of [0, e_max]; None = all of it
    # M6: perturber inclined relative to Uranus's plane; q gains (sin i cos Ω, sin i sin Ω).
    # Prior: cos i uniform on [cos inc_max, 1] and Ω uniform, i.e. density ∝ 1/cos i in q.
    inclined: bool = False
    inc_max_deg: float = INC_PRIOR_MAX_DEG

    @property
    def e_bounds(self) -> tuple[float, float]:
        return (0.0, self.e_max) if self.e_range is None else self.e_range

    @property
    def free(self) -> np.ndarray:
        return np.array([v is None for v in self.held])

    @property
    def ndim(self) -> int:
        return 4 + int(self.free.sum()) + 2 * int(self.inclined)

    @property
    def names(self) -> tuple[str, ...]:
        inc = INCLINATION_Q_NAMES if self.inclined else ()
        return URANUS_Q_NAMES + tuple(n for n, f in zip(PERTURBER_Q_NAMES, self.free,
                                                         strict=True) if f) + inc  # fmt: skip

    def with_box(self, box: UranusBox) -> PosteriorSpace:
        from dataclasses import replace

        return replace(self, uranus_box=box)

    def with_a_slice(self, a_lo_au: float, a_hi_au: float) -> PosteriorSpace:
        from dataclasses import replace

        if not (PRIOR_A_AU[0] <= a_lo_au < a_hi_au <= PRIOR_A_AU[1]) or not self.free[1]:
            raise ValueError("a slice must lie inside the a prior, with a free")
        return replace(self, a_au_range=(float(a_lo_au), float(a_hi_au)))

    def with_e_slice(self, e_lo: float, e_hi: float) -> PosteriorSpace:
        from dataclasses import replace

        if not (0.0 <= e_lo < e_hi <= self.e_max) or not (self.free[2] and self.free[3]):
            raise ValueError("an e slice must lie inside [0, e_max], with (k, h) free")
        return replace(self, e_range=(float(e_lo), float(e_hi)))

    @property
    def slice_weight(self) -> float:
        """Prior mass of this cell (a slice × e slice) within the full prior: log-uniform in a,
        uniform in e (1 for no slicing)."""
        lo, hi = self.a_au_range
        e_lo, e_hi = self.e_bounds
        return (math.log(hi / lo) / math.log(PRIOR_A_AU[1] / PRIOR_A_AU[0])
                * (e_hi - e_lo) / self.e_max)  # fmt: skip

    def full_perturber_q(self, q: np.ndarray) -> np.ndarray:
        y = np.array([0.0 if v is None else v for v in self.held])
        y[self.free] = q[4 : 4 + int(self.free.sum())]
        return y

    def theta(self, q: np.ndarray) -> np.ndarray:
        """q → 11-entry forward-model θ (13 when ``inclined``: + i, Ω relative to Uranus)."""
        theta_u = uranus_x_to_theta(np.concatenate([q[:4], self.pq]))
        log10_m, ln_a, k, h, lam0 = self.full_perturber_q(q)
        out = [log10_m, math.exp(ln_a), k * k + h * h, math.atan2(h, k) % TWO_PI, lam0 % TWO_PI]
        if self.inclined:
            sc, ss = q[-2], q[-1]
            out += [math.asin(min(math.hypot(sc, ss), 1.0)), math.atan2(ss, sc) % TWO_PI]
        return np.concatenate([theta_u, out])

    def q_from_theta(self, theta: np.ndarray) -> np.ndarray:
        from fork.inference.fit import uranus_theta_to_x

        xu = uranus_theta_to_x(theta[:6])[:4]
        s = math.sqrt(theta[8])
        y = np.array([theta[6], math.log(theta[7]), s * math.cos(theta[9]),
                      s * math.sin(theta[9]), theta[10]])  # fmt: skip
        inc = []
        if self.inclined:
            si = math.sin(theta[11])
            inc = [si * math.cos(theta[12]), si * math.sin(theta[12])]
        return np.concatenate([xu, y[self.free], inc])

    def log_prior(self, q: np.ndarray) -> float:
        """0 inside the support, −∞ outside (the density is flat in q; unnormalised). With
        ``inclined``, the (sin i cos Ω, sin i sin Ω) pair adds −ln cos i (uniform cos i)."""
        log10_m, ln_a, k, h, _ = self.full_perturber_q(q)
        if self.free[0] and not PRIOR_LOG10_M_EARTH[0] <= log10_m <= PRIOR_LOG10_M_EARTH[1]:
            return -math.inf
        lo_a, hi_a = self.a_au_range
        if self.free[1] and not math.log(lo_a) <= ln_a <= math.log(hi_a):
            return -math.inf
        if self.free[2] or self.free[3]:
            e_lo, e_hi = self.e_bounds
            if not e_lo <= k * k + h * h <= e_hi:
                return -math.inf
        if self.uranus_box is not None:
            if not self.uranus_box.contains(np.asarray(q[:4])):
                return -math.inf
        elif math.hypot(q[1], q[2]) >= _URANUS_MAX_E or q[0] <= 0.0:
            return -math.inf
        if self.inclined:
            s2 = q[-2] ** 2 + q[-1] ** 2
            if s2 > math.sin(math.radians(self.inc_max_deg)) ** 2:
                return -math.inf
            return -0.5 * math.log1p(-s2)
        return 0.0

    # ---------------------------------------------------------------- nested sampling

    def _check_transformable(self) -> None:
        if self.uranus_box is None:
            raise ValueError("prior_transform needs a Uranus box (PosteriorSpace.with_box)")
        if self.free[2] != self.free[3]:
            raise ValueError("(√e cos ϖ, √e sin ϖ) must be both free or both held")

    @property
    def unit_names(self) -> tuple[str, ...]:
        """Unit-cube coordinates, in order: Uranus 4, then the free perturber ones with the
        (k, h) pair replaced by (e, ϖ)."""
        names = ["u_a", "u_k", "u_h", "u_lam0"]
        for name, f in zip(("log10_m", "ln_a", "e", "varpi", "lam0"), self.free, strict=True):
            if f:
                names.append(name)
        if self.inclined:
            names += ["cos_i", "node"]
        return tuple(names)

    @property
    def periodic_unit_dims(self) -> list[int]:
        """Unit-cube indices that wrap: the perturber's ϖ and λ₀ (Uranus's λ₀ lives in a box
        of width ~1e-3 rad, so it is not periodic)."""
        return [i for i, n in enumerate(self.unit_names) if n in ("varpi", "lam0", "node")]

    def prior_transform(self, u: np.ndarray) -> np.ndarray:
        """Unit cube → q, pushing the uniform measure onto exactly the ``log_prior`` support:
        Uranus box; log10 m = 2u; ln a uniform; e = e_lo + (e_hi − e_lo)·u, ϖ = 2πu →
        (√e cos ϖ, √e sin ϖ); λ₀ = 2πu."""
        self._check_transformable()
        box = self.uranus_box
        u = np.asarray(u, dtype=float)
        q = np.empty(self.ndim)
        q[:4] = box.lo + 2.0 * np.asarray(box.half_width) * u[:4]
        names = self.unit_names
        j, out = 4, 4
        while j < len(names):
            n = names[j]
            if n == "log10_m":
                lo, hi = PRIOR_LOG10_M_EARTH
                q[out] = lo + (hi - lo) * u[j]
            elif n == "ln_a":
                lo, hi = math.log(self.a_au_range[0]), math.log(self.a_au_range[1])
                q[out] = lo + (hi - lo) * u[j]
            elif n == "e":  # consumes (e, ϖ) → (k, h)
                e_lo, e_hi = self.e_bounds
                s = math.sqrt(e_lo + (e_hi - e_lo) * u[j])
                varpi = TWO_PI * u[j + 1]
                q[out], q[out + 1] = s * math.cos(varpi), s * math.sin(varpi)
                j += 1
                out += 1
            elif n == "lam0":
                q[out] = TWO_PI * u[j]
            elif n == "cos_i":  # consumes (cos i, Ω) → (sin i cos Ω, sin i sin Ω)
                cos_i = 1.0 - (1.0 - math.cos(math.radians(self.inc_max_deg))) * u[j]
                s = math.sqrt(max(0.0, 1.0 - cos_i * cos_i))
                node = TWO_PI * u[j + 1]
                q[out], q[out + 1] = s * math.cos(node), s * math.sin(node)
                j += 1
                out += 1
            j += 1
            out += 1
        return q


def free_a_space(pq: tuple[float, float], e_max: float = E_PRIOR_MAX,
                 uranus_box: UranusBox | None = None, *,
                 inclined: bool = False) -> PosteriorSpace:  # fmt: skip
    return PosteriorSpace(e_max=e_max, pq=pq, uranus_box=uranus_box, inclined=inclined)


def held_shape_space(pq: tuple[float, float], a_au: float, e: float, varpi_rad: float,
                     uranus_box: UranusBox | None = None) -> PosteriorSpace:  # fmt: skip
    """a, e, ϖ held (e.g. Le Verrier's published shape); m, λ₀ and Uranus sampled."""
    s = math.sqrt(e)
    return PosteriorSpace((None, math.log(a_au), s * math.cos(varpi_rad), s * math.sin(varpi_rad),
                           None), e_max=max(e, E_PRIOR_MAX), pq=pq,
                          uranus_box=uranus_box)  # fmt: skip
