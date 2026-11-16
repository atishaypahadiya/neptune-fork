"""Contact times of a transit as bracketed roots (docs/phase2/identifiability_full.md,
amendment 2).

A contact is a root of a gap function g(t): the apparent separation of the centres minus the
contact radius, at an offset t [days] from a reference epoch. Over a day either side of an
inferior conjunction the separation has one minimum, at mid-transit, so a transit gives
exactly two roots, one before mid-transit (ingress, contact I or II) and one after (egress,
III or IV), and a miss gives none.

1. Mid-transit: the minimum of g over the window [−1, +1] day about the epoch, by
   golden-section search.
2. No contact if g at mid-transit is not below zero (no transit, or a grazing one in which
   the limbs never reach this contact).
3. The chosen root (``side``: −1 ingress, +1 egress) by bisection between mid-transit and
   mid-transit ∓ 1 day, where g must be positive; no contact otherwise.
4. No contact if the root lies more than a day from the epoch, or if |g| there is not below
   the root tolerance.

Newton's method from the epoch, which this replaces, fails where mid-transit falls near the
epoch (g′ ≈ 0 there, so the first step leaves the window) and can reach the other contact of
the same transit.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

INGRESS, EGRESS = -1, 1
WINDOW_DAYS = 1.0
GOLDEN_ITER = 80  # bracket 2 d × 0.618⁸⁰ ≈ 4 × 10⁻¹⁷ d
BISECT_ITER = 64  # bracket 1 d × 2⁻⁶⁴ ≈ 5 × 10⁻²⁰ d, below a double's spacing near the root
_INV_PHI = (np.sqrt(5.0) - 1.0) / 2.0

Gap = Callable[[np.ndarray, np.ndarray], np.ndarray]  # (row indices, offsets [d]) -> g [rad]


def mid_transit(gap: Gap, n: int, window: float = WINDOW_DAYS) -> tuple[np.ndarray, np.ndarray]:
    """Offset [d] of the minimum of g over [−window, +window] for each of n rows, and g there."""
    idx = np.arange(n)
    a, b = np.full(n, -window), np.full(n, window)
    c, d = b - _INV_PHI * (b - a), a + _INV_PHI * (b - a)
    fc, fd = gap(idx, c), gap(idx, d)
    for _ in range(GOLDEN_ITER):
        left = fc < fd  # the minimum lies in [a, d]
        b = np.where(left, d, b)
        a = np.where(left, a, c)
        c_new = np.where(left, b - _INV_PHI * (b - a), d)
        d_new = np.where(left, c, a + _INV_PHI * (b - a))
        ev = np.where(left, b - _INV_PHI * (b - a), a + _INV_PHI * (b - a))
        fe = gap(idx, ev)
        fc, fd = np.where(left, fe, fd), np.where(left, fc, fe)
        c, d = c_new, d_new
    t = 0.5 * (a + b)
    return t, gap(idx, t)


def contact_roots(gap: Gap, side: np.ndarray, root_tol: float,
                  window: float = WINDOW_DAYS) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """Offsets [d] of the chosen contact for each row (NaN where there is none) and |g| at it."""
    side = np.asarray(side)
    if not np.all(np.isin(side, (INGRESS, EGRESS))):
        raise ValueError("side must be −1 (ingress) or +1 (egress) for every row")
    n = len(side)
    idx = np.arange(n)
    tm, gm = mid_transit(gap, n, window)
    far = tm + side * window
    gf = gap(idx, far)
    ok = (gm < 0.0) & (gf > 0.0)
    lo, hi = tm.copy(), far.copy()  # g(lo) < 0 < g(hi) where ok
    for _ in range(BISECT_ITER):
        mid = 0.5 * (lo + hi)
        neg = gap(idx, mid) < 0.0
        lo, hi = np.where(neg, mid, lo), np.where(neg, hi, mid)
    g_lo, g_hi = gap(idx, lo), gap(idx, hi)
    t = np.where(np.abs(g_lo) <= np.abs(g_hi), lo, hi)
    res = np.minimum(np.abs(g_lo), np.abs(g_hi))
    ok &= (np.abs(t) <= window) & (res < root_tol)
    return np.where(ok, t, np.nan), np.where(ok, res, np.nan)
