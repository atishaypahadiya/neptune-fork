"""ΔT = TT − UT for 1600 to 1900, its uncertainty, and the step from a printed UT to TDB.

**Series.** The polynomials of Espenak and Meeus, *Five Millennium Canon of Solar Eclipses:
−1999 to +3000* (NASA/TP-2006-214141, 2006), as printed in "Polynomial Expressions for Delta T"
(https://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html). For 1600 to 1900 they follow the
telescopic values of Morrison and Stephenson, "Historical values of the Earth's clock error ΔT
and the calculation of eclipses", *J. Hist. Astron.* 35, 327–336 (2004).

**Uncertainty.** The standard errors that page's companion table gives with the historical
values ("Delta T (ΔT)", https://eclipse.gsfc.nasa.gov/SEcat5/deltat.html): 20 s at 1600, 5 s at
1700, 2 s at 1750, 1 s at 1800, "< 1" s at 1850 and 1900 (taken here as 1 s), interpolated
linearly in between. The table gives no correlations between years; how to carry this error
into a likelihood is left to the plan that uses it.

Nothing here is fitted: the coefficients and standard errors are copied as printed.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

FIRST_YEAR, LAST_YEAR = 1600.0, 1900.0
J2000_JD = 2451545.0
SECONDS_PER_DAY = 86400.0

# (start year, ascending coefficients in t = y − start); each segment ends where the next begins
SEGMENTS = (
    (1600.0, (120.0, -0.9808, -0.01532, 1.0 / 7129.0)),
    (1700.0, (8.83, 0.1603, -0.0059285, 0.00013336, -1.0 / 1174000.0)),
    (1800.0, (13.72, -0.332447, 0.0068612, 0.0041116, -0.00037436, 0.0000121272,
              -0.0000001699, 0.000000000875)),
    (1860.0, (7.62, 0.5737, -0.251754, 0.01680668, -0.0004473624, 1.0 / 233174.0)),
)  # fmt: skip

# the published table: year, ΔT [s], standard error [s] ("< 1" entered as 1)
PUBLISHED = (
    (1600.0, 120.0, 20.0),
    (1700.0, 9.0, 5.0),
    (1750.0, 13.0, 2.0),
    (1800.0, 14.0, 1.0),
    (1850.0, 7.0, 1.0),
    (1900.0, -3.0, 1.0),
)


def decimal_year(jd: ArrayLike) -> np.ndarray:
    """Julian year from a Julian Date. (The Canon takes y = year + (month − 0.5)/12; the two
    differ by under half a month, which changes ΔT by under 0.05 s in 1600–1900.)"""
    return 2000.0 + (np.asarray(jd, dtype=float) - J2000_JD) / 365.25


def _check(y: np.ndarray) -> None:
    if np.any((y < FIRST_YEAR) | (y > LAST_YEAR)) or not np.all(np.isfinite(y)):
        raise ValueError(f"ΔT here covers {FIRST_YEAR:.0f} to {LAST_YEAR:.0f} only")


def delta_t(year: ArrayLike) -> np.ndarray:
    """ΔT = TT − UT [s] at a decimal year in 1600–1900."""
    y = np.asarray(year, dtype=float)
    _check(y)
    out = np.empty_like(y)
    starts = np.array([s for s, _ in SEGMENTS])
    seg = np.clip(np.searchsorted(starts, y, side="right") - 1, 0, len(SEGMENTS) - 1)
    for k, (start, coeffs) in enumerate(SEGMENTS):
        m = seg == k
        out[m] = np.polynomial.polynomial.polyval(y[m] - start, coeffs)
    return out


def delta_t_sigma(year: ArrayLike) -> np.ndarray:
    """The published standard error of ΔT [s], interpolated linearly between table years."""
    y = np.asarray(year, dtype=float)
    _check(y)
    yr, _, sig = np.array(PUBLISHED).T
    return np.interp(y, yr, sig)


def tdb_minus_tt(jd: ArrayLike) -> np.ndarray:
    """TDB − TT [s] by the usual two-term approximation, 1.657 ms sin g + 0.014 ms sin 2g with g
    the Sun's mean anomaly (a test checks it against astropy's full series to 50 μs)."""
    g = np.radians(357.53 + 0.98560028 * (np.asarray(jd, dtype=float) - J2000_JD))
    return 0.001657 * np.sin(g) + 0.000014 * np.sin(2.0 * g)


def ut_to_tdb(jd_ut: ArrayLike) -> np.ndarray:
    """JD (UT) → JD (TDB): add ΔT at that date, then TDB − TT."""
    jd = np.asarray(jd_ut, dtype=float)
    jd_tt = jd + delta_t(decimal_year(jd)) / SECONDS_PER_DAY
    return jd_tt + tdb_minus_tt(jd_tt) / SECONDS_PER_DAY


def tdb_to_ut(jd_tdb: ArrayLike, iterations: int = 3) -> np.ndarray:
    """The inverse of :func:`ut_to_tdb`, by fixed-point iteration (ΔT changes by under 1 s a
    year, so three iterations reach rounding)."""
    target = np.asarray(jd_tdb, dtype=float)
    jd = target.copy()
    for _ in range(iterations):
        jd = jd - (ut_to_tdb(jd) - target)
    return jd
