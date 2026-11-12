"""The step-3 figure (docs/notes/P1B_adams_fit.png): the posterior of the geocentric longitude
on 1 October 1846 for each fit (H1, H2, F), Gaussian solid and Student-t dotted in M10's two
hues, with the targets marked (Neptune, S4, S5, S6, Le Verrier June and August; S7 apart, as an
extrapolation); below, F's marginal in a with 38.4 and 37.25 AU and Neptune's distance marked.
Two hues, validated as a pair in M10; the facets keep every panel to two series."""

from __future__ import annotations

from pathlib import Path

import numpy as np

BLUE, ORANGE = "#2a6fdb", "#d9480f"
INK, MUTED, GRID = "#1f1f1f", "#777777", "#dddddd"
TITLES = {"h1": "H1, a held at 38.4 AU", "h2": "H2, a held at 37.25 AU", "f": "F, a free"}
MARKS = ("neptune", "S4", "S5", "S6", "LV-June", "LV-Aug")
LABELS = {"neptune": "Neptune", "LV-June": "LV June", "LV-Aug": "LV Aug"}


def _style(ax) -> None:
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK, labelsize=8)


def _density(ax, x: np.ndarray, lo: float, hi: float, color: str, ls: str, label: str) -> None:
    edges = np.linspace(lo, hi, 121)
    h, _ = np.histogram(x, bins=edges, density=True)
    mid = 0.5 * (edges[1:] + edges[:-1])
    ax.plot(mid, h, color=color, linewidth=1.6, linestyle=ls, label=label)


def draw(payload: dict, raw: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    epoch = payload["reading"]["epoch"]
    fig = plt.figure(figsize=(11, 7.2))
    gs = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.85), hspace=0.5, wspace=0.25)
    # longitudes unwrapped about the circular mean of all runs, so a posterior straddling 360°
    # is drawn in one piece
    lons = np.concatenate([raw[v]["lon_1oct_deg"] for v in raw])
    centre = float(np.degrees(np.angle(np.exp(1j * np.radians(lons)).mean()))) % 360.0

    def unwrap(x):
        return (np.asarray(x, dtype=float) - centre + 180.0) % 360.0 - 180.0 + centre

    marks = unwrap([payload["neptune_places"][epoch]["geo_lon"]] + [
        payload["targets"][n]["places"][epoch]["geo_lon"] for n in MARKS[1:]])  # fmt: skip
    u = unwrap(lons)
    lo = min(float(np.percentile(u, 0.2)), float(marks.min())) - 1.0
    hi = max(float(np.percentile(u, 99.8)), float(marks.max())) + 1.0
    for j, fit in enumerate(("h1", "h2", "f")):
        ax = fig.add_subplot(gs[0, j])
        _style(ax)
        any_series = False
        for err, color, ls, name in (("g", BLUE, "-", "Gaussian"), ("t", ORANGE, ":", "Student-t")):
            v = f"{fit}_{err}"
            if v in raw:
                _density(ax, unwrap(raw[v]["lon_1oct_deg"]), lo, hi, color, ls, name)
                any_series = True
        ax.set_title(TITLES[fit], fontsize=10, color=INK, loc="left")
        ax.set_xlabel(f"geocentric longitude of date, {epoch} (°)", fontsize=8, color=MUTED)
        ax.set_yticks([])
        top = ax.get_ylim()[1] if any_series else 1.0
        for k, name in enumerate(MARKS):
            x = float(marks[k])
            if not lo <= x <= hi:
                continue
            ax.axvline(x, color=INK if name == "neptune" else MUTED, linewidth=0.8,
                       ymax=0.9)  # fmt: skip
            ax.text(x, top * (0.98 - 0.065 * k), LABELS.get(name, name), fontsize=6.5,
                    color=INK if name == "neptune" else MUTED, ha="left" if k % 2 else "right",
                    va="top")  # fmt: skip
        if j == 0:
            ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax = fig.add_subplot(gs[1, :])
    _style(ax)
    for err, color, ls, name in (("g", BLUE, "-", "Gaussian"), ("t", ORANGE, ":", "Student-t")):
        v = f"f_{err}"
        if v in raw:
            _density(ax, raw[v]["a_au"], 20.0, 60.0, color, ls, name)
    for k, (x, label) in enumerate(
        ((38.4, "38.4 (H1)"), (37.25, "37.25 (H2)"), (30.07, "Neptune"))
    ):
        ax.axvline(x, color=INK if label == "Neptune" else MUTED, linewidth=0.8, ymax=0.9)
        ax.text(x, ax.get_ylim()[1] * (0.98 - 0.08 * (k % 2)), label, fontsize=7, ha="center",
                va="top", color=INK if label == "Neptune" else MUTED)  # fmt: skip
    ax.set_title("F, posterior in a", fontsize=10, color=INK, loc="left")
    ax.set_xlabel("a (AU)", fontsize=8, color=MUTED)
    ax.set_yticks([])
    ax.set_xlim(20, 60)
    s7 = payload["targets"]["S7"]["mean_lon_1846_10_01_deg"]
    fig.text(0.01, 0.005, f"S7 (2 Sept 1846 extrapolation, mean longitude {s7:.2f}° on {epoch}) "
             "is not drawn as a solution; see the JSON.", fontsize=7, color=MUTED)  # fmt: skip
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
