"""Reference frames, precession, and state-vector ↔ orbital-element conversion.

Frames used in this project (CLAUDE.md rule 4 — every angle names its frame):

``ecl_j2000``
    Heliocentric ecliptic and mean equinox of J2000 *as used by JPL Horizons*: ICRF rotated about
    x by the IAU 1976 obliquity ε₀ = 84381.448″, no frame bias. This is astropy's
    ``HeliocentricEclipticIAU76``. All Horizons vectors (``fork.data.horizons``) are in this
    frame, so it is the project's internal frame.

``ecl_of_date``
    Heliocentric mean ecliptic and mean equinox of date (IAU 2006 precession), astropy's
    ``HeliocentricMeanEcliptic(equinox=t)``. Historical longitudes (Bouvard, Le Verrier) are
    referred to this kind of frame.

``HeliocentricMeanEcliptic(equinox=J2000)`` differs from ``ecl_j2000`` by the ICRS frame bias and
the IAU76→IAU2006 obliquity change — tens of mas, see ``docs/notes/M1.md``. Going through astropy
for both ends of the precession keeps that difference accounted for rather than ignored.

Orbital elements are frame-agnostic: they are referred to whatever frame the state vector is in
(for us, ``ecl_j2000``). Angles are radians; a in AU; μ = G(M + m) in AU³/day².
"""

from __future__ import annotations

import math
import warnings

import astropy.units as u
import numpy as np
from astropy.coordinates import (
    CartesianRepresentation,
    HeliocentricEclipticIAU76,
    HeliocentricMeanEcliptic,
)
from astropy.time import Time
from erfa import ErfaWarning
from numba import njit
from numpy.typing import ArrayLike

TWO_PI = 2.0 * math.pi

# --------------------------------------------------------------------------------------------
# Time (boundary helper)
# --------------------------------------------------------------------------------------------


def jd_tdb_from_iso(iso: str | list[str]) -> np.ndarray | float:
    """Julian Date (TDB) of an ISO calendar date/time given *in TDB*.

    Pre-1960 UTC is not well defined, so historical epochs are specified directly in TDB; the
    difference from UT is ~10 s in the 18th–19th century, irrelevant at our 1″ scale for Uranus.
    """
    jd = Time(iso, scale="tdb").jd
    return float(jd) if np.ndim(jd) == 0 else np.asarray(jd, dtype=np.float64)


# --------------------------------------------------------------------------------------------
# Precession: ecliptic J2000 (Horizons) ↔ mean ecliptic of date
# --------------------------------------------------------------------------------------------


def precession_matrix_ecl_j2000_to_ecl_of_date(t_jd_tdb: ArrayLike) -> np.ndarray:
    """Rotation matrices R with ``x_ecl_of_date = R @ x_ecl_j2000``.

    Shape ``(*t.shape, 3, 3)``. Built by pushing the three basis vectors through astropy
    (``HeliocentricEclipticIAU76`` → ``HeliocentricMeanEcliptic(equinox=t)``) at a common obstime,
    so the Sun-offset legs of the affine transforms cancel and the result is a pure rotation.
    """
    t = np.asarray(t_jd_tdb, dtype=np.float64)
    flat = np.atleast_1d(t).ravel()
    time = Time(flat, format="jd", scale="tdb")

    eye = np.eye(3)
    # basis[k] = k-th unit vector, broadcast against the N epochs → data shape (3, N).
    xyz = np.broadcast_to(eye[:, :, None], (3, 3, flat.size))  # (component, basis k, N)
    rep = CartesianRepresentation(xyz * u.AU)
    src = HeliocentricEclipticIAU76(rep, obstime=time[None, :])
    with warnings.catch_warnings():
        # ERFA warns that pre-1900 dates are outside epv00's (Sun ephemeris) validity and that
        # UTC is "dubious" before 1960. Both feed only the Sun-offset legs, which cancel exactly
        # here because source and destination share obstime; precession itself (IAU 2006) is
        # valid for millennia. Tested against the IAU 2006 p_A polynomial back to 1690.
        warnings.simplefilter("ignore", ErfaWarning)
        dst = src.transform_to(
            HeliocentricMeanEcliptic(equinox=time[None, :], obstime=time[None, :])
        )
    out = dst.cartesian.xyz.to_value(u.AU)  # (component i, basis k, N): column k of R
    mats = np.moveaxis(out, -1, 0)  # (N, i, k)
    return mats.reshape((*t.shape, 3, 3))


def _unit_vector(lon_rad: np.ndarray, lat_rad: np.ndarray) -> np.ndarray:
    cl = np.cos(lat_rad)
    return np.stack([cl * np.cos(lon_rad), cl * np.sin(lon_rad), np.sin(lat_rad)], axis=-1)


def _lon_lat(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    lon = np.mod(np.arctan2(y, x), TWO_PI)
    lat = np.arctan2(z, np.hypot(x, y))
    return lon, lat


def ecl_j2000_to_ecl_of_date(
    lon_ecl_j2000_rad: ArrayLike, lat_ecl_j2000_rad: ArrayLike, t_jd_tdb: ArrayLike
) -> tuple[np.ndarray, np.ndarray]:
    """Precess heliocentric directions from ecliptic J2000 to mean ecliptic of date ``t``.

    Returns ``(lon_ecl_of_date_rad, lat_ecl_of_date_rad)``, lon in [0, 2π).
    """
    lon, lat, t = np.broadcast_arrays(
        np.asarray(lon_ecl_j2000_rad, dtype=np.float64),
        np.asarray(lat_ecl_j2000_rad, dtype=np.float64),
        np.asarray(t_jd_tdb, dtype=np.float64),
    )
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(t)
    vec_of_date = np.einsum("...ij,...j->...i", rot, _unit_vector(lon, lat))
    return _lon_lat(vec_of_date)


def ecl_of_date_to_ecl_j2000(
    lon_ecl_of_date_rad: ArrayLike, lat_ecl_of_date_rad: ArrayLike, t_jd_tdb: ArrayLike
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of :func:`ecl_j2000_to_ecl_of_date`.

    Returns ``(lon_ecl_j2000_rad, lat_ecl_j2000_rad)``, lon in [0, 2π).
    """
    lon, lat, t = np.broadcast_arrays(
        np.asarray(lon_ecl_of_date_rad, dtype=np.float64),
        np.asarray(lat_ecl_of_date_rad, dtype=np.float64),
        np.asarray(t_jd_tdb, dtype=np.float64),
    )
    rot = precession_matrix_ecl_j2000_to_ecl_of_date(t)
    vec_j2000 = np.einsum("...ji,...j->...i", rot, _unit_vector(lon, lat))  # Rᵀ v
    return _lon_lat(vec_j2000)


# --------------------------------------------------------------------------------------------
# Orbital elements (elliptic only), numba-compiled so forward models can call them.
# --------------------------------------------------------------------------------------------


# Speed of light [AU/day] (IAU 2012 au = 149 597 870 700 m exactly).
C_AU_PER_DAY: float = 299_792_458.0 * 86_400.0 / 149_597_870_700.0


def geocentric_lon_lat_ecl_j2000_light_time(
    target_pos_ecl_j2000_au: ArrayLike,
    target_vel_ecl_j2000_au_d: ArrayLike,
    earth_pos_ecl_j2000_au: ArrayLike,
    n_iter: int = 3,
) -> tuple[float, float, float]:
    """Astrometric (light-time corrected) geocentric (lon, lat) in ecl_j2000 [rad], and τ [day].

    Target and Earth positions share one origin at the observation time t. The target is moved
    back to t − τ along its velocity: r(t − τ) ≈ r(t) − v·τ, good to ~½·a·τ² ≈ 1e-8 AU for an
    outer planet (τ ≈ 0.17 d). No aberration, deflection or parallax — astrometric, geocentric.
    """
    r_t = np.asarray(target_pos_ecl_j2000_au, dtype=np.float64)
    v_t = np.asarray(target_vel_ecl_j2000_au_d, dtype=np.float64)
    r_e = np.asarray(earth_pos_ecl_j2000_au, dtype=np.float64)
    tau = 0.0
    for _ in range(n_iter):
        rho = r_t - v_t * tau - r_e
        tau = float(np.linalg.norm(rho)) / C_AU_PER_DAY
    rho = r_t - v_t * tau - r_e
    lon = math.atan2(rho[1], rho[0]) % (2.0 * math.pi)
    lat = math.atan2(rho[2], math.hypot(rho[0], rho[1]))
    return lon, lat, tau


@njit(cache=True)
def solve_kepler_elliptic(mean_anom_rad: float, e: float) -> float:
    """Eccentric anomaly E solving ``M = E - e sin E`` for 0 ≤ e < 1.

    Newton iteration from Danby's starter, with a bisection fallback on [-π, π] if Newton fails
    to converge (it should not for e < 1, but the fallback makes that a guarantee).
    """
    m = (mean_anom_rad + math.pi) % TWO_PI - math.pi  # in [-π, π)
    if e == 0.0:
        return m
    sign = 1.0 if math.sin(m) >= 0.0 else -1.0
    ecc = m + 0.85 * e * sign
    for _ in range(50):
        f = ecc - e * math.sin(ecc) - m
        fp = 1.0 - e * math.cos(ecc)
        d = f / fp
        ecc -= d
        if abs(d) < 1e-15:
            # one more step polishes the last ulp
            ecc -= (ecc - e * math.sin(ecc) - m) / (1.0 - e * math.cos(ecc))
            return ecc
    # Bisection fallback: f(E) = E - e sin E - M is monotone increasing on [-π, π].
    lo, hi = -math.pi, math.pi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid - e * math.sin(mid) - m > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


@njit(cache=True)
def elements_to_state(
    a: float,
    e: float,
    inc_rad: float,
    raan_rad: float,
    argp_rad: float,
    mean_anom_rad: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Keplerian elements → heliocentric (r [AU], v [AU/day]) in the elements' reference frame.

    ``raan`` = longitude of ascending node Ω, ``argp`` = argument of pericenter ω,
    ``mean_anom`` = mean anomaly M at the epoch. Requires a > 0, 0 ≤ e < 1.
    """
    if not (a > 0.0) or not (0.0 <= e < 1.0):
        raise ValueError("elements_to_state: need a > 0 and 0 <= e < 1")
    ecc = solve_kepler_elliptic(mean_anom_rad, e)
    cos_e, sin_e = math.cos(ecc), math.sin(ecc)
    b_over_a = math.sqrt((1.0 - e) * (1.0 + e))
    n = math.sqrt(mu / (a * a * a))
    denom = 1.0 - e * cos_e

    # perifocal frame: P toward pericenter, Q 90° ahead in the direction of motion
    xp = a * (cos_e - e)
    yp = a * b_over_a * sin_e
    vxp = -a * n * sin_e / denom
    vyp = a * n * b_over_a * cos_e / denom

    co, so = math.cos(raan_rad), math.sin(raan_rad)
    cw, sw = math.cos(argp_rad), math.sin(argp_rad)
    ci, si = math.cos(inc_rad), math.sin(inc_rad)
    # P and Q unit vectors in the reference frame (R3(-Ω) R1(-i) R3(-ω))
    px = co * cw - so * sw * ci
    py = so * cw + co * sw * ci
    pz = sw * si
    qx = -co * sw - so * cw * ci
    qy = -so * sw + co * cw * ci
    qz = cw * si

    r = np.empty(3)
    v = np.empty(3)
    r[0] = xp * px + yp * qx
    r[1] = xp * py + yp * qy
    r[2] = xp * pz + yp * qz
    v[0] = vxp * px + vyp * qx
    v[1] = vxp * py + vyp * qy
    v[2] = vxp * pz + vyp * qz
    return r, v


@njit(cache=True)
def state_to_elements(
    r: np.ndarray, v: np.ndarray, mu: float
) -> tuple[float, float, float, float, float, float]:
    """(r, v) → (a, e, inc_rad, raan_rad, argp_rad, mean_anom_rad), elliptic orbits only.

    Angles Ω, ω, M are returned in [0, 2π); inc in [0, π].
    Degenerate conventions: if the orbit is planar (sin i ≈ 0) Ω = 0 and ω is measured from the
    reference x-axis; if circular (e ≈ 0) ω = 0 and M is measured from the node. In both cases
    the well-defined combinations ϖ = Ω + ω and λ = ϖ + M are preserved.
    """
    rx, ry, rz = r[0], r[1], r[2]
    vx, vy, vz = v[0], v[1], v[2]
    rmag = math.sqrt(rx * rx + ry * ry + rz * rz)
    v2 = vx * vx + vy * vy + vz * vz

    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    hmag = math.sqrt(hx * hx + hy * hy + hz * hz)
    hxy = math.sqrt(hx * hx + hy * hy)

    inv_a = 2.0 / rmag - v2 / mu
    if not (inv_a > 0.0):
        raise ValueError("state_to_elements: orbit is not bound (elliptic)")
    a = 1.0 / inv_a

    # eccentricity vector: (v × h)/μ − r̂
    ex = (vy * hz - vz * hy) / mu - rx / rmag
    ey = (vz * hx - vx * hz) / mu - ry / rmag
    ez = (vx * hy - vy * hx) / mu - rz / rmag
    e = math.sqrt(ex * ex + ey * ey + ez * ez)

    inc = math.atan2(hxy, hz)

    # node direction n̂ = ẑ × ĥ normalized; planar → x-axis
    if hxy > 1e-15 * hmag:
        raan = math.atan2(hx, -hy)
    else:
        raan = 0.0
    nx, ny = math.cos(raan), math.sin(raan)
    # in-plane unit vector 90° ahead of the node: ĥ × n̂
    hhx, hhy, hhz = hx / hmag, hy / hmag, hz / hmag
    mx = -hhz * ny
    my = hhz * nx
    mz = hhx * ny - hhy * nx

    # argument of latitude of the body
    u_lat = math.atan2(rx * mx + ry * my + rz * mz, rx * nx + ry * ny)

    if e > 1e-14:
        argp = math.atan2(ex * mx + ey * my + ez * mz, ex * nx + ey * ny)
    else:
        argp = 0.0
    nu = u_lat - argp

    # true → eccentric → mean anomaly (half-angle form is well conditioned for all e < 1)
    ecc = 2.0 * math.atan2(
        math.sqrt(1.0 - e) * math.sin(0.5 * nu), math.sqrt(1.0 + e) * math.cos(0.5 * nu)
    )
    mean_anom = ecc - e * math.sin(ecc)

    return (
        a,
        e,
        inc,
        raan % TWO_PI,
        argp % TWO_PI,
        mean_anom % TWO_PI,
    )
