"""A3b coverage: one 68 % interval per noise seed → docs/notes/M5b_coverage.png.

    uv run python docs/notes/scripts/m5b_coverage_plot.py

Reads only the committed summary ``site/data/m5b_posteriors.json`` (free a, e ≤ 0.3): each
seed's posterior median and central 68 % interval of the perturber's 1846-09-23 geocentric
longitude, and whether true Neptune (325.87°) lies inside. A miss is drawn dashed, with a
hollow median marker and a direct label, so it reads without colour; the orange only
reinforces that.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "M5b_coverage.png"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
MISS = "#eb6834"  # categorical slot 2; passes the dataviz palette validator on SURFACE
RUN = "free_a_e0.3"


def main() -> None:
    post = json.loads((ROOT / "site" / "data" / "m5b_posteriors.json").read_text())
    truth = post["truth"]["galle_lon_ecl_of_date_deg"]
    cal = post["a3b_calibration"][RUN]
    inside = {r["seed"]: r["truth_in_68"] for r in cal["per_seed"]}
    rows = []
    for seed in post["run_seeds"][RUN]:
        lon = post["datasets"][f"seed_{seed}"]["runs"][RUN]["galle_lon_ecl_of_date_deg"]
        hit = lon["lo68"] <= truth <= lon["hi68"]
        assert hit == inside[seed], f"seed {seed}: interval and summary flag disagree"
        rows.append((seed, lon["lo68"], lon["median"], lon["hi68"], hit))

    fig, ax = plt.subplots(figsize=(8, 4.4), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for seed, lo, med, hi, hit in rows:
        y = -seed
        c = INK if hit else MISS
        ax.plot([lo, hi], [y, y], color=c, lw=2, solid_capstyle="round",
                ls="-" if hit else (0, (4, 2.5)))  # fmt: skip
        ax.plot(med, y, "o", ms=6.5, mec=c, mew=1.6, mfc=c if hit else SURFACE, zorder=3)
        if not hit:
            ax.text(hi + 0.6, y, "miss", color=INK_2, va="center", fontsize=9)
    ax.axvline(truth, color=INK_MUTED, lw=1, zorder=0)
    ax.text(truth + 0.3, 0.75, f"true Neptune {truth:.2f}°", color=INK_2, fontsize=9,
            va="bottom")  # fmt: skip

    ax.set_yticks([-s for s, *_ in rows], [f"seed {s}" for s, *_ in rows])
    ax.set_ylim(-len(rows) + 0.4, 1.4)
    ax.set_xlabel("perturber geocentric longitude on 1846-09-23 (ecliptic of date, °)",
                  color=INK_2)  # fmt: skip
    ax.set_title(f"Truth inside the 68 % interval in {cal['n_truth_in_68']} of "
                 f"{cal['n_seeds']} synthetic datasets", loc="left", color=INK,
                 fontsize=11)  # fmt: skip
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(INK_MUTED)
    ax.tick_params(colors=INK_2, length=0, labelsize=9)
    ax.tick_params(axis="x", length=3, color=INK_MUTED)
    fig.text(0.01, 0.01, "free a, e ≤ 0.3, partitioned nested sampling; dot = posterior "
             "median, line = central 68 % interval", color=INK_MUTED, fontsize=8)  # fmt: skip
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
