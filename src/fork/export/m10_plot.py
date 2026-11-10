"""The M10 figure (docs/notes/M10_profile.png) from docs/notes/M10_profile.json: Δ ln L, the
fitted mass and the Galle-night longitude against a for data sets A and B under Student-t errors
(solid) and, dotted in the same hue, M9's Gaussian curves, with S3's posterior marginal in a;
beside them, one small panel per leave-one-out set (Δ ln L, orange) over A's Student-t curve
(blue). Two hues, validated as a pair; the facets keep every panel to two series."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

BLUE, ORANGE = "#2a6fdb", "#d9480f"
INK, MUTED, GRID = "#1f1f1f", "#777777", "#dddddd"
JPL_GALLE_LON_DEG = 325.873  # fork.export.m6c.JPL_GALLE_LON_DEG, a comparison target only
LEFT = {"all": ("A, all rows", BLUE), "post1781": ("B, post-1781 rows", ORANGE)}
RIGHT = ("no_flamsteed", "no_bradley", "no_mayer", "no_lemonnier")


def _xy(ds: dict, key: str) -> tuple[list[float], list[float]]:
    pts = [p for p in ds["points"] if not p["reference"] and p.get(key) is not None]
    return [p["a_au"] for p in pts], [p[key] for p in pts]


def _style(ax) -> None:
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK, labelsize=8)


def _ticks(ax, grid: list[float]) -> None:
    ax.set_xticks(grid)
    ax.set_xticklabels([f"{v:g}" if float(v).is_integer() else "" for v in grid], fontsize=8)


def draw(out: dict, m9_out: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(4, 2, width_ratios=[3, 2])
    ax_l = fig.add_subplot(gs[0:2, 0])
    ax_m = fig.add_subplot(gs[2, 0], sharex=ax_l)
    ax_g = fig.add_subplot(gs[3, 0], sharex=ax_l)
    facets = [fig.add_subplot(gs[i, 1], sharex=ax_l) for i in range(4)]
    for ax in (ax_l, ax_m, ax_g, *facets):
        _style(ax)
    marg = out["s3_posterior_marginal_a"]
    grid = marg["grid_au"]
    pa = marg["a_au"]
    for lo, hi, alpha in ((pa["lo95"], pa["hi95"], 0.08), (pa["lo68"], pa["hi68"], 0.12)):
        for ax in (ax_l, ax_m, ax_g):
            ax.axvspan(lo, hi, color=MUTED, alpha=alpha, linewidth=0)
    lp = [np.nan if v is None else v for v in marg["ln_p_rel"]]
    ax_l.plot(grid, lp, color=MUTED, linestyle=":", linewidth=1.5, marker="s", markersize=3,
              label=f"S3 posterior, ln p(a) − max ({marg['tightness']})")  # fmt: skip
    for level, txt in ((-0.5, "−0.5"), (-2.0, "−2")):
        for ax in (ax_l, *facets):
            ax.axhline(level, color=MUTED, linewidth=0.8)
        ax_l.annotate(f"Δ ln L = {txt}", (44.5, level), fontsize=8, color=MUTED, va="bottom",
                      ha="right")  # fmt: skip
    lows = [v for v in lp if np.isfinite(v)]
    for name, (label, color) in LEFT.items():
        ds = out["datasets"].get(name)
        g9 = m9_out["datasets"].get(name)
        if g9 is not None:
            x, y = _xy(g9, "dlnl")
            ax_l.plot(x, y, color=color, linestyle=":", linewidth=1.2, marker="o", markersize=3,
                      markerfacecolor="none", label=f"{label}, Gaussian (M9)")  # fmt: skip
            lows += y
        if ds is None:
            continue
        for ax, key in ((ax_l, "dlnl"), (ax_m, "mass_earth"), (ax_g, "galle_lon_deg")):
            x, y = _xy(ds, key)
            if not x:
                continue
            ax.plot(x, y, color=color, linewidth=1.8, marker="o", markersize=4,
                    label=f"{label}, Student-t" if ax is ax_l else None)  # fmt: skip
            if ax is ax_l:
                lows += y
        for rname, r in ds["reference_points"].items():
            for ax, key in ((ax_l, "dlnl"), (ax_m, "mass_earth"), (ax_g, "galle_lon_deg")):
                ax.plot([r["a_au"]], [r[key]], marker="o", markersize=8, markerfacecolor="none",
                        markeredgecolor=color, linestyle="none")  # fmt: skip
            if name == "all":
                ax_l.annotate("Neptune's a" if rname == "neptune" else "Le Verrier's a",
                              (r["a_au"], r["dlnl"]), xytext=(0, 10), textcoords="offset points",
                              fontsize=8, color=INK, ha="center")  # fmt: skip
    ax_g.axhline(JPL_GALLE_LON_DEG, color=MUTED, linewidth=0.8)
    ax_g.annotate("JPL's Neptune", (44.5, JPL_GALLE_LON_DEG), fontsize=8, color=MUTED,
                  va="bottom", ha="right")  # fmt: skip
    ax_l.set_ylabel("Δ ln L(a) = −ΔQ/2")
    ax_m.set_ylabel("fitted mass (M⊕)")
    ax_g.set_ylabel("longitude on Galle's night (°)")
    ax_g.set_xlabel("semi-major axis a held (AU)")
    _ticks(ax_g, grid)
    ax_l.legend(loc="lower left", fontsize=8, frameon=False)
    lo = min(v for v in lows if np.isfinite(v)) if lows else -5.0
    ax_l.set_ylim(max(lo, -40.0) - 1.0, 1.5)

    a_t = out["datasets"].get("all")
    for ax, name in zip(facets, RIGHT, strict=True):
        ds = out["datasets"].get(name)
        if a_t is not None:
            x, y = _xy(a_t, "dlnl")
            ax.plot(x, y, color=BLUE, linewidth=1.2, marker="o", markersize=3, label="A")
        if ds is not None:
            x, y = _xy(ds, "dlnl")
            ax.plot(x, y, color=ORANGE, linewidth=1.8, marker="o", markersize=4,
                    label=f"{ds['label']} ({ds['reading']['reading']})")  # fmt: skip
        ax.legend(loc="lower left", fontsize=8, frameon=False)
        ax.set_ylim(max(lo, -40.0) - 1.0, 1.5)
        ax.set_ylabel("Δ ln L(a)")
    _ticks(facets[-1], grid)
    facets[-1].set_xlabel("semi-major axis a held (AU)")
    r = a_t["reading"]["reading"] if a_t else "not run"
    fig.suptitle("M10: profile likelihood in the semi-major axis under Student-t errors "
                 f"(pre-registered diagnostic; A reads {r})", fontsize=10)  # fmt: skip
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    import sys

    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/notes/M10_profile.json")
    draw(json.loads(src.read_text()), json.loads(Path("docs/notes/M9_profile.json").read_text()),
         src.with_suffix(".png"))  # fmt: skip
