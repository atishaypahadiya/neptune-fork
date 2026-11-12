"""Figures of the Phase 2 validation, from docs/phase2/validation.json.

    uv run python -m fork.phase2.validation plot

- ``docs/phase2/validation_v4.png``: Mercury's longitude against the JPL ephemeris, 1800–2000,
  without the 1PN term (top) and with it (bottom); two panels because the scales differ by
  three orders of magnitude.
- ``docs/phase2/validation_v5_step.png``: Mercury's position error against step size, second
  and fourth order, 1PN on.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fork.phase2.validation import OUT_JSON

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
BLUE, ORANGE = "#2a78d6", "#d9622b"  # checked for colour-vision separation on this surface


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(color=GRID, lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=INK_2, labelsize=9)


def plot() -> None:
    d = json.loads(OUT_JSON.read_text(encoding="utf-8"))
    v4 = d["v4"]["cases"]
    year = np.array(v4["newtonian"]["series"]["year"])
    fig, (a, b) = plt.subplots(2, 1, figsize=(9, 6), dpi=150, sharex=True, facecolor=SURFACE)
    for ax in (a, b):
        _style(ax)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.set_ylabel("model − JPL [″]", color=INK_2, fontsize=9)
    a.plot(year, v4["newtonian"]["series"]["dlon"], color=ORANGE, lw=0.6)
    slope = v4["newtonian"]["perihelion_residual_slope_arcsec_per_century"]
    a.set_title(f"Without the 1PN term: the perihelion falls behind by {abs(slope):.2f}″ per "
                "century", loc="left", color=INK, fontsize=10.5)  # fmt: skip
    b.plot(year, v4["with_1pn"]["series"]["dlon"], color=BLUE, lw=0.6, label="1PN")
    b.plot(year, v4["with_1pn_and_j2"]["series"]["dlon"], color=ORANGE, lw=0.6,
           label="1PN and J2")  # fmt: skip
    b.set_title(f"With it: within {v4['with_1pn']['max_abs_dlon']:.2f}″ over two centuries",
                loc="left", color=INK, fontsize=10.5)  # fmt: skip
    b.legend(frameon=False, fontsize=9, labelcolor=INK, loc="upper center", ncols=2)
    b.set_xlabel("year", color=INK_2, fontsize=9)
    fig.suptitle("Mercury's heliocentric longitude against the JPL ephemeris", x=0.01,
                 ha="left", color=INK, fontsize=12)  # fmt: skip
    fig.tight_layout()
    out = OUT_JSON.with_name("validation_v4.png")
    fig.savefig(out, facecolor=SURFACE)
    print(f"wrote {out}")

    rows = [r for r in d["v5_step"]["rows"] if r["term"] == "1pn"]
    fig, ax = plt.subplots(figsize=(7, 4.4), dpi=150, facecolor=SURFACE)
    _style(ax)
    for order, color in ((2, ORANGE), (4, BLUE)):
        r = sorted((x for x in rows if x["order"] == order), key=lambda x: x["dt_days"])
        dt, err = [x["dt_days"] for x in r], [x["max_arcsec"] for x in r]
        ax.loglog(dt, err, "o-", color=color, lw=2, ms=6)
        ax.annotate(f"order {order}", (dt[-1], err[-1]), xytext=(6, 0),
                    textcoords="offset points", va="center", color=INK, fontsize=9)  # fmt: skip
    ax.set_xlabel("step [days]", color=INK_2, fontsize=9)
    ax.set_ylabel("largest error in Mercury's place [″]", color=INK_2, fontsize=9)
    years = d["v5_step"]["years_each_way"]
    ax.set_title(f"Step size, 1PN on, {years:g} years each way from 1900", loc="left",
                 color=INK, fontsize=11)  # fmt: skip
    ax.margins(x=0.18)
    fig.tight_layout()
    out = OUT_JSON.with_name("validation_v5_step.png")
    fig.savefig(out, facecolor=SURFACE)
    print(f"wrote {out}")
