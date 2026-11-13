"""Synthetic Mercury under each competing model, and how the models' predictions differ.

A tool, built and tested only. Running it as an experiment needs a plan, and nothing here
reads an observation of Mercury: the epochs are an argument and the "data" are the forward
model's own output.

Each model is given the one parameter value that produces the same secular perihelion
advance (``target``, 42.98″ per century by default) by its analytic rate:

    gr      scale of the 1PN term                       (secular.gr_rate)
    j2      the Sun's J2, spin axis from the context    (secular.j2_rates, dω/dt + dΩ/dt)
    ring    mass of a circular ring at ``ring_radius``  (secular.ring_gauss_rate)
    planet  mass of a point planet at ``planet_a``      (secular.laplace_lagrange_rate)

``compare`` integrates each against the Newtonian baseline from the same initial state and
reports, per model, the secular rates of ϖ, Ω, i and e relative to the baseline and the
size of the periodic part of the difference in mean longitude. No element is refitted: a
real test would absorb part of these differences into Mercury's elements, so the numbers
are upper bounds on what could tell the models apart.

    uv run python -m fork.phase2.identifiability        # the report on the CI fixture's epochs
"""

from __future__ import annotations

import json
import math

import numpy as np
from numpy.typing import ArrayLike

from fork.model import inner_system as IS
from fork.model.bodies import G, masses_msun
from fork.phase2 import secular as S

MODELS = ("gr", "j2", "ring", "planet")
TARGET_ARCSEC_PER_CENTURY = 42.98
RING_RADIUS_AU = 0.2
PLANET_A_AU = 0.2
ARCSEC = 180.0 * 3600.0 / math.pi
CENTURY_DAYS = 36525.0


def _mu() -> float:
    return G * (1.0 + masses_msun(("mercury",))[0])


def equator_inclination(theta_mercury: ArrayLike, spin_axis: ArrayLike) -> float:
    """Inclination [rad] of Mercury's orbit to the plane normal to ``spin_axis``."""
    _, _, inc, raan = np.asarray(theta_mercury, dtype=float)[:4]
    normal = np.array([math.sin(inc) * math.sin(raan), -math.sin(inc) * math.cos(raan),
                       math.cos(inc)])  # fmt: skip
    axis = np.asarray(spin_axis, dtype=float)
    return math.acos(float(normal @ axis) / float(np.linalg.norm(axis)))


def matched_parameters(theta_mercury: ArrayLike, context: IS.InnerSystemContext,
                       target: float = TARGET_ARCSEC_PER_CENTURY,
                       ring_radius: float = RING_RADIUS_AU, planet_a: float = PLANET_A_AU,
                       planet_lam0: float = 0.0) -> dict[str, tuple[float, ...]]:  # fmt: skip
    """The extra entries of θ per model that give ``target`` [″ per century] of perihelion
    advance by the analytic rate. Rates are linear in each parameter."""
    a, e = (float(v) for v in np.asarray(theta_mercury)[:2])
    mu, want = _mu(), target / S.ARCSEC_PER_CENTURY  # rad per day
    inc_eq = equator_inclination(theta_mercury, context.spin_axis)
    dw, dnode = S.j2_rates(a, e, inc_eq, mu, 1.0)
    planet = want / S.laplace_lagrange_rate(a, mu, 1.0, planet_a)
    return {
        "gr": (want / S.gr_rate(a, e, mu),),
        "j2": (want / (dw + dnode),),
        "ring": (want / S.ring_gauss_rate(a, e, mu, 1.0, ring_radius), ring_radius),
        "planet": (planet, planet_a, planet_lam0),
    }


def flags_for(model: str) -> IS.Flags:
    if model == "newtonian":
        return IS.NO_FLAGS
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; choose from {MODELS}")
    return IS.Flags(**{model: True})


def synthetic(model: str, theta_mercury: ArrayLike, t_jd: ArrayLike,
              context: IS.InnerSystemContext, extra: tuple[float, ...] = (),
              sigma_arcsec: float = 0.0,
              rng: np.random.Generator | None = None) -> np.ndarray:  # fmt: skip
    """Mercury's heliocentric longitude [rad] at t_jd under ``model``, with optional
    Gaussian noise of ``sigma_arcsec`` (then ``rng`` is required)."""
    theta = np.concatenate([np.asarray(theta_mercury, dtype=float)[:6], extra])
    lon = IS.predict(theta, t_jd, context, flags_for(model))
    if sigma_arcsec > 0.0:
        if rng is None:
            raise ValueError("noise needs an explicit rng")
        lon = (lon + rng.normal(0.0, sigma_arcsec / ARCSEC, lon.shape)) % (2.0 * math.pi)
    return lon


def _elements(model: str, theta_mercury, t_jd, context, extra=()) -> np.ndarray:
    theta = np.concatenate([np.asarray(theta_mercury, dtype=float)[:6], extra])
    x, v = IS.states(theta, t_jd, context, flags_for(model))
    return S.elements_series(x, v, _mu())


def _wrap(d: np.ndarray) -> np.ndarray:
    return (d + math.pi) % (2.0 * math.pi) - math.pi


def compare(theta_mercury: ArrayLike, t_jd: ArrayLike, context: IS.InnerSystemContext,
            **matched) -> dict:  # fmt: skip
    """Per model, relative to the Newtonian baseline over ``t_jd``: secular rates of ϖ, Ω
    and i [″ per century] and of e [per century]; the RMS of the difference in mean
    longitude λ = Ω + ω + M after a straight line is removed [″]; and the RMS of the raw
    difference in true longitude [″]. Then, per pair of models, the RMS of the difference
    of their longitudes. Epochs need not be evenly spaced; at least three are needed."""
    t = np.atleast_1d(np.asarray(t_jd, dtype=float))
    if t.size < 3:
        raise ValueError("at least three epochs")
    params = matched_parameters(theta_mercury, context, **matched)
    base = _elements("newtonian", theta_mercury, t, context)
    lon0 = synthetic("newtonian", theta_mercury, t, context)
    tc = t - t.mean()
    out: dict = {"n_epochs": int(t.size), "span_years": float((t.max() - t.min()) / 365.25),
                 "target_arcsec_per_century": matched.get("target", TARGET_ARCSEC_PER_CENTURY),
                 "models": {}, "pairs": {}}  # fmt: skip
    lons = {}
    for m in MODELS:
        el = _elements(m, theta_mercury, t, context, params[m])
        d = el - base
        varpi = d[:, 3] + d[:, 4]
        lam = _wrap(d[:, 3] + d[:, 4] + d[:, 5])
        lam = np.unwrap(lam)
        periodic = lam - np.polyval(np.polyfit(tc, lam, 1), tc)
        lons[m] = synthetic(m, theta_mercury, t, context, params[m])
        out["models"][m] = {
            "parameters": list(params[m]),
            "varpi_rate": S.slope(t, varpi) * CENTURY_DAYS * ARCSEC,
            "node_rate": S.slope(t, d[:, 3]) * CENTURY_DAYS * ARCSEC,
            "inclination_rate": S.slope(t, d[:, 2]) * CENTURY_DAYS * ARCSEC,
            "eccentricity_rate": S.slope(t, d[:, 1]) * CENTURY_DAYS,
            "mean_longitude_periodic_rms": float(np.sqrt((periodic**2).mean())) * ARCSEC,
            "longitude_rms_vs_newtonian": float(
                np.sqrt((_wrap(lons[m] - lon0) ** 2).mean())) * ARCSEC,
        }  # fmt: skip
    for i, a in enumerate(MODELS):
        for b in MODELS[i + 1 :]:
            out["pairs"][f"{a}:{b}"] = float(
                np.sqrt((_wrap(lons[a] - lons[b]) ** 2).mean())) * ARCSEC  # fmt: skip
    return out


def fixture_problem() -> tuple[np.ndarray, np.ndarray, IS.InnerSystemContext]:
    """(Mercury's θ, epochs, context) from the CI fixture of Part 5: 1900.0 states, epochs
    about every five years from 1800 to 2000."""
    from fork.phase2 import reference as R

    ref = R.load(R.FIXTURE)
    es = ref["epoch_states"]
    return IS.mercury_theta_from_state(es[0]), ref["jd"], IS.InnerSystemContext(R.EPOCH_JD, es[1:])


def main() -> None:
    theta, jd, ctx = fixture_problem()
    print(json.dumps(compare(theta, jd, ctx), indent=1))


if __name__ == "__main__":
    main()
