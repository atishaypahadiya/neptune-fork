"""The M9 figure (docs/notes/M9_profile.png) from docs/notes/M9_profile.json: Δ ln L, the fitted
mass and the Galle-night longitude against a, both data sets, with the headline posterior's
marginal in a and the two reference points marked."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

SERIES = {"all": ("all rows", "#2a6fdb", "-"), "post1781": ("post-1781 rows", "#d9480f", "--")}
INK, MUTED, GRID = "#1f1f1f", "#777777", "#dddddd"
JPL_GALLE_LON_DEG = 325.873  # fork.export.m6c.JPL_GALLE_LON_DEG, a comparison target only


def _xy(ds: dict, key: str) -> tuple[list[float], list[float]]:
    pts = [p for p in ds["points"] if not p["reference"] and key in p and p[key] is not None]
    return [p["a_au"] for p in pts], [p[key] for p in pts]


def draw(out: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 9.0), sharex=True, constrained_layout=True)
    ax_l, ax_m, ax_g = axes
    marg = out["posterior_marginal_a"]
    for ax in axes:
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=INK, labelsize=9)
    # the posterior's 68 % and 95 % intervals in a, behind everything
    pa = marg["a_au"]
    for lo, hi, alpha in ((pa["lo95"], pa["hi95"], 0.08), (pa["lo68"], pa["hi68"], 0.12)):
        for ax in axes:
            ax.axvspan(lo, hi, color=MUTED, alpha=alpha, linewidth=0)
    # the marginal density in a (ln p − max) on the Δ ln L axis
    g = marg["grid_au"]
    lp = [np.nan if v is None else v for v in marg["ln_p_rel"]]
    ax_l.plot(g, lp, color=MUTED, linestyle=":", linewidth=1.5, marker="s", markersize=3,
              label="headline posterior, ln p(a) − max")  # fmt: skip
    for level, txt in ((-0.5, "−0.5"), (-2.0, "−2")):
        ax_l.axhline(level, color=MUTED, linewidth=0.8)
        ax_l.annotate(f"Δ ln L = {txt}", (44.5, level), fontsize=8, color=MUTED, va="bottom",
                      ha="right")  # fmt: skip
    for ds_name, (label, color, ls) in SERIES.items():
        ds = out["datasets"].get(ds_name)
        if ds is None:
            continue
        for ax, key in ((ax_l, "dlnl"), (ax_m, "mass_earth"), (ax_g, "galle_lon_deg")):
            x, y = _xy(ds, key)
            if not x:
                continue
            ax.plot(x, y, color=color, linestyle=ls, linewidth=1.6, marker="o", markersize=4,
                    label=label if ax is ax_l else None)  # fmt: skip
            ax.annotate(label, (x[-1], y[-1]), xytext=(4, 0), textcoords="offset points",
                        fontsize=8, color=INK, va="center")  # fmt: skip
        for name, r in ds["reference_points"].items():
            for ax, key in ((ax_l, "dlnl"), (ax_m, "mass_earth"), (ax_g, "galle_lon_deg")):
                ax.plot([r["a_au"]], [r[key]], marker="o", markersize=8, markerfacecolor="none",
                        markeredgecolor=color, linestyle="none")  # fmt: skip
            if ds_name == "all":
                ax_l.annotate("Neptune's a" if name == "neptune" else "Le Verrier's a",
                              (r["a_au"], r["dlnl"]), xytext=(0, 10), textcoords="offset points",
                              fontsize=8, color=INK, ha="center")  # fmt: skip
    ax_g.axhline(JPL_GALLE_LON_DEG, color=MUTED, linewidth=0.8)
    ax_g.annotate("JPL's Neptune", (44.5, JPL_GALLE_LON_DEG), fontsize=8, color=MUTED,
                  va="bottom", ha="right")  # fmt: skip
    ax_l.set_ylabel("Δ ln L(a)")
    ax_m.set_ylabel("fitted mass (M⊕)")
    ax_g.set_ylabel("longitude on Galle's night (°)")
    ax_g.set_xlabel("semi-major axis a held (AU)")
    ax_g.set_xticks(g)  # every grid point; the half-AU ones unlabelled so labels do not collide
    ax_g.set_xticklabels([f"{v:g}" if float(v).is_integer() else "" for v in g], fontsize=8)
    ax_l.legend(loc="lower left", fontsize=8, frameon=False)
    lo = min(
        v
        for v in lp
        + [
            p["dlnl"]
            for ds in out["datasets"].values()
            if ds
            for p in ds["points"]
            if p.get("dlnl") is not None
        ]
        if np.isfinite(v)
    )
    ax_l.set_ylim(max(lo, -40.0) - 1.0, 1.5)
    r = out["datasets"].get("all") or out["datasets"].get("post1781")
    reading = "" if r is None else f"; A reads {r['reading']['reading']}"
    fig.suptitle("M9: profile likelihood in the semi-major axis (pre-registered diagnostic"
                 f"{reading})", fontsize=10)  # fmt: skip
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    import sys

    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/notes/M9_profile.json")
    draw(json.loads(src.read_text()), src.with_suffix(".png"))
