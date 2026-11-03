"""M5b posterior: mass–distance ridge and the 1846 longitude → docs/notes/M5b_posterior.png.

    uv run python docs/notes/scripts/m5b_ridge_plot.py [seed]      # default seed 0

Reads raw chains from data/cache/chains/ (written by ``python -m fork.export.m5b``; not
committed). Left: (a, m) samples under the two free-a priors (e ≤ 0.3, e ≤ 0.1) with true
Neptune marked. Right: the derived 1846-09-23 geocentric longitude under both priors and
Le Verrier's shape, with true Neptune marked.
"""

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import json

import matplotlib.pyplot as plt
import numpy as np

from fork.export.m5b import CHAIN_DIR, _spaces
from fork.inference.fit import uranus_theta_to_x
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.forward import ForwardContext, uranus_inc_raan_at_epoch

EPOCH = 2371628.5
ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "M5b_posterior.png"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"  # categorical slots 1–3 (all-pairs safe)
N_LON = 1500


def main() -> None:
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    key = f"seed_{seed}"
    post = json.loads((ROOT / "site" / "data" / "m5b_posteriors.json").read_text())
    truth = post["truth"]
    ctx = ForwardContext.from_horizons(EPOCH)
    inc, raan = uranus_inc_raan_at_epoch(EPOCH)
    pq = tuple(uranus_theta_to_x(np.array([19.2, 0.0, inc, raan, 0.0, 0.0]))[4:])
    spaces = _spaces(pq)
    geo = GalleGeometry.from_horizons()
    runs = (("free_a_e0.3", "free a, e ≤ 0.3", C1), ("free_a_e0.1", "free a, e ≤ 0.1", C2),
            ("le_verrier_shape", "Le Verrier's shape (a, e, ϖ held)", C3))  # fmt: skip

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), dpi=150, facecolor=SURFACE,
                                  gridspec_kw={"width_ratios": [1.1, 1]})  # fmt: skip
    for a in (ax, ax2):
        a.set_facecolor(SURFACE)
        a.spines[["top", "right"]].set_visible(False)
        a.grid(color=GRID, lw=0.8)
    for name, label, color in runs:
        z = np.load(CHAIN_DIR / f"{key}_{name}.npz")
        q = z["samples"]
        space = spaces[name][0]
        idx = np.linspace(0, len(q) - 1, min(len(q), N_LON)).astype(int)
        thetas = np.array([space.theta(v) for v in q[idx]])
        if name != "le_verrier_shape":
            ax.scatter(thetas[:, 7], 10 ** thetas[:, 6], s=4, color=color, alpha=0.35, lw=0,
                       label=label, rasterized=True)  # fmt: skip
        lon = np.array([math.degrees(galle_lon_of_date_rad(t, ctx, geo)) for t in thetas])
        ax2.hist(lon, bins=np.arange(290, 361, 1.5), histtype="step", lw=2, color=color,
                 density=True, label=label)  # fmt: skip
    ax.plot(truth["a_au"], truth["mass_earth"], "*", ms=14, mfc=INK, mec=SURFACE, mew=1.5,
            zorder=5, label="true Neptune")  # fmt: skip
    ax.set_xlabel("perturber a [AU]  (osculating, 1781)")
    ax.set_ylabel("perturber mass [M⊕]")
    ax.set_title(f"Mass–distance posterior (seed {seed})", loc="left", fontsize=10.5, color=INK)
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK, markerscale=3)
    nep = truth["galle_lon_ecl_of_date_deg"]
    ax2.axvline(nep, color=INK, lw=1.2, ls=":")
    ax2.annotate(f"true Neptune {nep:.2f}°", (nep, ax2.get_ylim()[1]), xytext=(4, -12),
                 textcoords="offset points", fontsize=8, color=INK_2)  # fmt: skip
    ax2.set_xlabel("perturber geocentric longitude, 1846-09-23 [° ecliptic of date]")
    ax2.set_ylabel("posterior density [1/°]")
    ax2.set_title("Derived 1846 longitude", loc="left", fontsize=10.5, color=INK)
    ax2.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK)
    fig.text(0.01, 0.01, "Synthetic data, emcee (48 walkers, ≥ 50 τ); priors log10 m ∈ [0, 2], "
             "a log-uniform on [20, 60] AU. Chains: data/cache/chains/ (not committed)",
             fontsize=7.5, color=INK_MUTED)  # fmt: skip
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
