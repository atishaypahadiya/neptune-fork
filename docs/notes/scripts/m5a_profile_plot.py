"""Profile likelihood of the 1846-09-23 longitude → docs/notes/M5a_profile.png.

    uv run python docs/notes/scripts/m5a_profile_plot.py

Reads site/data/m5a_profile.json. One panel per mode (free a, Le Verrier a = 38.4 AU): Δχ² vs
the perturber's geocentric ecliptic-of-date longitude on Galle's night, for noise seeds 0–3 and
noise-free data. Dotted vertical: true Neptune. Horizontal: Δχ² = 1.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "site" / "data" / "m5a_profile.json"
OUT = Path(__file__).resolve().parents[1] / "M5a_profile.png"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
SEED_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")  # categorical slots 1–4, fixed order
Y_MAX = 12.0


def main() -> None:
    d = json.loads(DATA.read_text())
    nep = d["true_neptune_lon_ecl_of_date_deg"]
    titles = {"free_a": "Free a", "le_verrier": f"Le Verrier mode, a = {d['le_verrier_a_au']} AU"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), dpi=150, facecolor=SURFACE, sharey=True)
    for ax, mode in zip(axes, ("free_a", "le_verrier"), strict=True):
        ax.set_facecolor(SURFACE)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(1.0, color=INK_MUTED, lw=0.8)
        ax.axvline(nep, color=INK_MUTED, lw=1.0, ls=":")
        ax.annotate(f"true Neptune {nep:.2f}°", (nep, 0.55 * Y_MAX), xytext=(4, 0),
                    textcoords="offset points", fontsize=8, color=INK_2, zorder=6,
                    bbox={"fc": SURFACE, "ec": "none", "pad": 1.5})  # fmt: skip
        for k, ds in enumerate(d["datasets"]):
            m = ds["modes"][mode]
            p = m["profile"]
            noise_free = ds["seed"] is None
            color = INK if noise_free else SEED_COLORS[k]
            label = "noise-free" if noise_free else f"seed {ds['seed']}"
            ax.plot(p["lon_ecl_of_date_deg"], p["dchi2"], color=color, lw=2,
                    ls="--" if noise_free else "-", label=label, clip_on=True)  # fmt: skip
            mn = m["overall_min"]
            ax.plot(mn["lon_ecl_of_date_deg"], 0.0, "o", ms=8, mfc=color, mec=SURFACE, mew=2,
                    zorder=5, clip_on=False)  # fmt: skip
            if p["dchi2"][-1] < Y_MAX:  # direct label only where the line ends on the plot
                ax.annotate(label, (p["lon_ecl_of_date_deg"][-1], p["dchi2"][-1]), xytext=(4, 0),
                            textcoords="offset points", va="center", fontsize=8, color=INK_2,
                            annotation_clip=False)  # fmt: skip
        ax.set_xlim(300, 350)
        ax.set_ylim(0, Y_MAX)
        ax.set_title(titles[mode], loc="left", fontsize=10.5, color=INK)
        ax.set_xlabel("perturber geocentric longitude on 1846-09-23 [° ecliptic of date]")
    axes[0].set_ylabel("profile Δχ²  (other parameters re-optimised)")
    for ax, loc in zip(axes, ("upper left", "upper right"), strict=True):
        ax.legend(loc=loc, frameon=False, fontsize=8, labelcolor=INK)
    fig.text(0.01, 0.045, "Fixing a = 38.4 AU narrows the valley and excludes the truth: a wrong "
             "assumption buys false confidence.", fontsize=9, color=INK)  # fmt: skip
    fig.text(0.01, 0.01, "Synthetic data (σ = 12″ pre-1781, 3″ after), dt = 12 d; circles: overall "
             "minimum; Δχ² clipped at 12. Data: site/data/m5a_profile.json", fontsize=7.5,
             color=INK_MUTED)  # fmt: skip
    fig.tight_layout(rect=(0, 0.08, 0.97, 1))
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
