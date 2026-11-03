"""M5a: residuals after the point fits → docs/notes/M5a_residuals.png.

    uv run python docs/notes/scripts/m5a_residuals_plot.py

Reads site/data/m5a_fits.json (no refitting). Top: noise-free residual curves (truth − best fit)
for no perturber (A2), Le Verrier mode (a = 38.4 AU) and free a. Bottom: the free-a fit's
observation residuals with ±σ. Right: 1846-09-23 geocentric longitude of each local minimum vs
true Neptune.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "site" / "data" / "m5a_fits.json"
OUT = Path(__file__).resolve().parents[1] / "M5a_residuals.png"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
COLORS = {"no_perturber": "#2a78d6", "le_verrier": "#eb6834", "free_a": "#1baf7a"}  # validated
SHORT = {"no_perturber": "no perturber", "le_verrier": "Le Verrier", "free_a": "free a"}
LABEL_DY = {"no_perturber": 0, "le_verrier": 7, "free_a": -7}  # end-label offsets [pt]


def main() -> None:
    d = json.loads(DATA.read_text())
    years = np.array(d["years"])
    fig = plt.figure(figsize=(11, 6.4), dpi=150, facecolor=SURFACE)
    gs = fig.add_gridspec(2, 2, width_ratios=[2.3, 1], height_ratios=[2, 1], wspace=0.55)
    ax, ax_obs, ax_g = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])
    for a in (ax, ax_obs, ax_g):
        a.set_facecolor(SURFACE)
        a.spines[["top", "right"]].set_visible(False)
        a.grid(axis="y", color=GRID, lw=0.8)
    for a in (ax, ax_obs):
        a.axhline(0, color=INK_MUTED, lw=0.8)
        a.axvline(1781, color=INK_MUTED, lw=0.8, ls=":")
        a.set_xlim(1688, 1858)

    for c in d["curves"]:
        r = np.array(c["residual_arcsec"])
        fit = d["fits"].get(c["id"])
        extra = "" if fit is None else (f": m {fit['mass_earth']:.1f} M⊕, a {fit['a_au']:.1f} AU, "
                                        f"e {fit['e']:.2f}, χ²/dof {fit['chi2_per_dof']:.2f}")  # fmt: skip
        ax.plot(years, r, color=COLORS[c["id"]], lw=2, label=c["label"] + extra)
        ax.annotate(SHORT[c["id"]], (years[-1], r[-1]), xytext=(5, LABEL_DY[c["id"]]),
                    textcoords="offset points",
                    va="center", fontsize=8.5, color=INK_2)  # fmt: skip
    ax.set_ylabel("truth − fit [″]")
    ax.set_title("M5a point fits on synthetic data (seed 0): residuals after fit", loc="left",
                 fontsize=11, color=INK)  # fmt: skip
    ax.legend(loc="lower left", frameon=False, fontsize=8, labelcolor=INK)

    o = d["observations"]
    yr, sig = np.array(o["year"]), np.array(o["sigma_arcsec"])
    res = np.array(o["residual_arcsec"]["free_a"])
    ax_obs.errorbar(yr, res, yerr=sig, fmt="o", ms=3.5, color=COLORS["free_a"], elinewidth=0.8,
                    capsize=0, alpha=0.8)  # fmt: skip
    ax_obs.set_ylabel("obs − free-a fit [″]")
    ax_obs.set_xlabel("year")

    true_lon = d["true_neptune"]["galle"]["lon_ecl_of_date_deg"]
    ax_g.axvline(0, color=INK, lw=1)
    rows = [(mode, k, m) for mode in ("le_verrier", "free_a")
            for k, m in enumerate(d["fits"][mode]["minima"])]  # fmt: skip
    labels = []
    for i, (mode, k, m) in enumerate(rows):
        ax_g.plot(m["galle_dlon_vs_neptune_deg"], i, "o", ms=8, mec=COLORS[mode], mew=2,
                  mfc=COLORS[mode] if m["in_m5_prior"] else SURFACE)  # fmt: skip
        tag = "best" if k == 0 else f"#{k + 1}"
        labels.append(f"{SHORT[mode]} {tag}  χ²/ν {m['chi2_per_dof']:.2f}  ×{m['n_starts']}")
        if k == 0:
            ax_g.annotate(f"{m['galle_dlon_vs_neptune_deg']:+.1f}°", (m["galle_dlon_vs_neptune_deg"], i),
                          xytext=(0, -11), textcoords="offset points", ha="center", fontsize=7.5,
                          color=INK_2)  # fmt: skip
    ax_g.set_yticks(range(len(rows)), labels, fontsize=7.5, color=INK_2)
    ax_g.tick_params(axis="y", length=0)
    ax_g.set_xlim(-180, 180)
    ax_g.set_xticks([-180, -90, 0, 90, 180])
    ax_g.set_ylim(len(rows) - 0.4, -1.2)
    ax_g.annotate(f"true Neptune {true_lon:.2f}° (of date)", (0, -0.9), ha="center", fontsize=7.5,
                  color=INK_2)  # fmt: skip
    ax_g.set_xlabel("1846-09-23 geocentric longitude\nminus true Neptune [°]")
    ax_g.set_title("local minima (filled: in M5 prior)", loc="left", fontsize=9,
                   color=INK)  # fmt: skip
    fig.text(0.01, 0.005, "Truth = JPL Horizons Uranus (incl. real Neptune) + noise σ = 12″ / 3″; "
             "Uranus i, Ω held at 1781 values; dt = 12 d", fontsize=7.5,
             color=INK_MUTED)  # fmt: skip
    fig.subplots_adjust(left=0.06, right=0.985, top=0.93, bottom=0.1, hspace=0.25)
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
