"""Uranus-only least-squares fit to RA/Dec measurements, and the Uranus prior box built from it
(the M6 counterpart of ``fit.fit_uranus_elements`` and ``priors.UranusBox.from_uranus_only_fit``).

docs/RULES.md rule 2, as in M4–M5: Uranus's in-plane (a, k, h, λ₀) are free and start from the
data alone; (p, q), i.e. i and Ω, are held at their epoch values. Declinations would constrain
i and Ω, but freeing them is a separate decision. The rule and its leak test are unchanged.

Data-only start: each measured direction is turned back into a geocentric ecliptic-J2000
longitude, using the inverse of the of-date rotation and ignoring aberration (~20″, harmless
for a start). A straight line is fitted to the post-1781 longitudes. Its slope is the mean
motion, since the annual parallax (±3°) averages out over 65 years, and Kepler's third law
gives a. The orbit starts circular, with λ₀ from the line's intercept at the epoch.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from scipy.optimize import least_squares

from fork.inference.fit import _MU_URANUS, TWO_PI, uranus_theta_to_x, uranus_x_to_theta
from fork.inference.priors import _FD_STEP, URANUS_BOX_N_SIGMA, UranusBox
from fork.model.forward import ForwardContext
from fork.model.observation import RaDecObservations

DISCOVERY_JD = 2371628.5  # 1781-03-13.0 TDB


def geocentric_lon_ecl_j2000_rad(obs: RaDecObservations) -> tuple[np.ndarray, np.ndarray]:
    """(jd, geocentric ecliptic-J2000 longitude) for every epoch that has both RA and Dec."""
    n = obs.geometry.jd_tt.size
    ra = np.full(n, np.nan)
    dec = np.full(n, np.nan)
    ra[obs.epoch_index[~obs.is_dec]] = obs.value_rad[~obs.is_dec]
    dec[obs.epoch_index[obs.is_dec]] = obs.value_rad[obs.is_dec]
    ok = np.isfinite(ra) & np.isfinite(dec)
    w = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=1)[ok]
    u = np.einsum("nji,nj->ni", obs.geometry.ecl_j2000_to_true_of_date[ok], w)  # Mᵀ w
    return obs.geometry.jd_tt[ok], np.mod(np.arctan2(u[:, 1], u[:, 0]), TWO_PI)


def initial_uranus_theta_radec(
    obs: RaDecObservations, epoch_jd_tdb: float, inc_raan: tuple[float, float]
) -> np.ndarray:
    jd, lon = geocentric_lon_ecl_j2000_rad(obs)
    sel = jd >= DISCOVERY_JD
    order = np.argsort(jd[sel])
    t = jd[sel][order] - epoch_jd_tdb
    n, lam0 = np.polyfit(t, np.unwrap(lon[sel][order]), 1)
    a = (_MU_URANUS / n**2) ** (1.0 / 3.0)
    return np.array([a, 0.0, inc_raan[0], inc_raan[1], 0.0, lam0 % TWO_PI])


def fit_uranus_radec(
    obs: RaDecObservations, context: ForwardContext, inc_raan: tuple[float, float],
    fixed_perturber_theta: np.ndarray | None = None,
    x0_in_plane: np.ndarray | None = None,
) -> tuple[np.ndarray, float, int]:  # fmt: skip
    """(Uranus θ (6), χ², dof) of the LM fit; χ² is the augmented (offset-marginal) one.
    dof = measurements − 4; the offsets are counted by their prior terms. The fit starts from
    the data-only line start, or from ``x0_in_plane`` (Uranus's (a, e cos ϖ, e sin ϖ, λ₀))
    when given (M10: a second start from the headline's u_ref)."""
    pert = np.empty(0) if fixed_perturber_theta is None else np.asarray(fixed_perturber_theta)
    x0 = uranus_theta_to_x(initial_uranus_theta_radec(obs, context.epoch_jd_tdb, inc_raan))
    if x0_in_plane is not None:
        x0 = np.concatenate([np.asarray(x0_in_plane, dtype=float), x0[4:]])

    def resid(x4: np.ndarray) -> np.ndarray:
        theta = uranus_x_to_theta(np.concatenate([x4, x0[4:]]))
        return obs.whitened(np.concatenate([theta, pert]), context)

    sol = least_squares(resid, x0[:4], method="lm", x_scale="jac", xtol=1e-12, ftol=1e-12,
                        gtol=1e-12, max_nfev=3000)  # fmt: skip
    theta = uranus_x_to_theta(np.concatenate([sol.x, x0[4:]]))
    return theta, float(np.sum(sol.fun**2)), obs.n_measurements - 4


def uranus_box_radec(
    obs: RaDecObservations, context: ForwardContext, inc_raan: tuple[float, float],
    n_sigma: float = URANUS_BOX_N_SIGMA,
) -> UranusBox:  # fmt: skip
    """The M5 box rule on RA/Dec data: centred on the Uranus-only fit, half-widths
    ``n_sigma`` × its formal σ scaled by √(χ²/dof)."""
    return uranus_box_radec_with_fits(obs, context, inc_raan, n_sigma)[0]


def uranus_box_radec_with_fits(
    obs: RaDecObservations, context: ForwardContext, inc_raan: tuple[float, float],
    n_sigma: float = URANUS_BOX_N_SIGMA, extra_starts: Sequence[np.ndarray] = (),
) -> tuple[UranusBox, list[dict]]:  # fmt: skip
    """``uranus_box_radec`` with the Uranus-only fit started from the data-only line start and
    from each of ``extra_starts`` (in-plane x, e.g. the headline's u_ref); the box is centred on
    the lowest-χ² solution. Returns the box and one record per start (name, χ², dof, the
    solution's in-plane x, whether it centres the box). With no extra starts this is the
    original rule, bit for bit (M10.md, pre-run change (a): M9 found the post-1781 line start
    converging to a poor local solution, χ² 2,900 for 488 dof)."""
    runs = []
    for name, start in [("line", None)] + [(f"start:{j}", s) for j, s in enumerate(extra_starts)]:
        try:
            runs.append((name, *fit_uranus_radec(obs, context, inc_raan, x0_in_plane=start)))
        except ValueError:  # the LM left the elements' domain (a ≤ 0 or e ≥ 1)
            if not extra_starts:
                raise
            runs.append((name, None, math.inf, 0))
    if all(r[1] is None for r in runs):
        raise ValueError("every start of the Uranus-only fit left the elements' domain")
    i_best = min(range(len(runs)), key=lambda i: runs[i][2])
    _, theta, chi2, dof = runs[i_best]
    x = uranus_theta_to_x(theta)
    fits = [{"start": name, "chi2": float(c2), "dof": int(df) if th is not None else None,
             "u": None if th is None else [float(v) for v in uranus_theta_to_x(th)[:4]],
             "centres_box": i == i_best}
            for i, (name, th, c2, df) in enumerate(runs)]  # fmt: skip

    def resid(x4: np.ndarray) -> np.ndarray:
        return obs.whitened(uranus_x_to_theta(np.concatenate([x4, x[4:]])), context)

    cols = []
    for i, h in enumerate(_FD_STEP):
        dx = np.zeros(4)
        dx[i] = h
        cols.append((resid(x[:4] + dx) - resid(x[:4] - dx)) / (2.0 * h))
    jac = np.array(cols).T
    sigma = np.sqrt(np.diag(np.linalg.inv(jac.T @ jac)))
    chi2_nu = chi2 / dof
    s = sigma * math.sqrt(max(chi2_nu, 1.0))
    box = UranusBox(tuple(float(v) for v in x[:4]), tuple(float(v) for v in n_sigma * s),
                    tuple(float(v) for v in s), float(chi2_nu))  # fmt: skip
    return box, fits
