"""A2 diagnostic: Uranus-only fit residuals → docs/notes/M4_residuals.png.

    uv run python docs/notes/scripts/m4_residuals_plot.py

Top: noise-free residual curve (Horizons truth − best no-perturber fit) for the A2 fit (all six
Uranus elements free) and, as a diagnostic, the fit with i and Ω held at Horizons values; noisy
observation residuals of the A2 fit as dots. Bottom: the slope of each curve (″/yr).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fork.data.synthetic import jd_tdb_jan1, make_synthetic_observations, true_uranus_lon_ecl_j2000_rad
from fork.inference.fit import fit_uranus_only, wrap_rad
from fork.model.forward import ForwardContext, UranusModel, uranus_inc_raan_at_epoch

ARC = 206264.80624709636
EPOCH = 2371628.5
OUT = Path(__file__).resolve().parents[1] / "M4_residuals.png"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
C_FREE, C_FIXED = "#2a78d6", "#eb6834"  # categorical slots 1, 2 (validated pair)


def jd_to_year(jd):
    return 2000.0 + (np.asarray(jd) - 2451545.0) / 365.25


def main() -> None:
    model = UranusModel(ForwardContext.from_horizons(EPOCH))
    obs, _truth = make_synthetic_observations(np.random.default_rng(0))
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    fits = {
        "free": fit_uranus_only(obs, model, inc_raan, free_inc_raan=True),
        "fixed": fit_uranus_only(obs, model, inc_raan),
    }
    years = np.arange(1690, 1846)
    grid = jd_tdb_jan1(years)
    truth_grid = true_uranus_lon_ecl_j2000_rad(grid)
    curves = {k: wrap_rad(truth_grid - model.predict(f.theta, grid)) * ARC for k, f in fits.items()}

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9, 6.4), dpi=150, sharex=True, facecolor=SURFACE,
                                  gridspec_kw={"height_ratios": [2.2, 1]})  # fmt: skip
    for a in (ax, ax2):
        a.set_facecolor(SURFACE)
        a.grid(axis="y", color=GRID, lw=0.8)
        a.spines[["top", "right"]].set_visible(False)
        a.axvspan(1815, 1830, color="#f1f0ec", zorder=0)
        a.axhline(0, color=INK_MUTED, lw=0.8)
    f = fits["free"]
    ax.errorbar(jd_to_year(obs.jd_tdb), f.residual_rad * ARC, yerr=obs.sigma_rad * ARC, fmt="o",
                ms=3.5, color=C_FREE, alpha=0.45, elinewidth=0.8, capsize=0,
                label="observations − fit (A2, σ = 12″ / 3″)")  # fmt: skip
    i_f, i_x = np.degrees(fits["free"].theta[2]), np.degrees(fits["fixed"].theta[2])
    o_f, o_x = np.degrees(fits["free"].theta[3]), np.degrees(fits["fixed"].theta[3])
    ax.plot(years, curves["free"], color=C_FREE, lw=2,
            label=f"truth − fit, all 6 elements free (A2): i = {i_f:.2f}°, Ω = {o_f:.0f}°")  # fmt: skip
    ax.plot(years, curves["fixed"], color=C_FIXED, lw=2,
            label=f"truth − fit, i, Ω held at Horizons (diagnostic): i = {i_x:.2f}°, Ω = {o_x:.0f}°")  # fmt: skip
    ax.set_ylabel("residual [″]")
    ax.set_title("Uranus without Neptune: residuals of the best Keplerian + Jupiter + Saturn fit",
                 loc="left", color=INK, fontsize=11)  # fmt: skip
    ax.legend(loc="lower left", frameon=False, fontsize=8.5, labelcolor=INK)
    ax.annotate("1815–1830", (1822.5, ax.get_ylim()[1]), xytext=(0, -12), textcoords="offset points",
                ha="center", fontsize=8.5, color=INK_2)  # fmt: skip
    for k, c in (("free", C_FREE), ("fixed", C_FIXED)):
        ax2.plot(years[1:] - 0.5, np.diff(curves[k]), color=c, lw=2)
    ax2.set_ylabel("slope [″/yr]")
    ax2.set_xlabel("year")
    ax2.set_xlim(1688, 1848)
    fig.text(0.01, 0.005, "Synthetic data (seed 0), truth = JPL Horizons (includes Neptune); "
             "fits never see Horizons Uranus", fontsize=7.5, color=INK_MUTED)  # fmt: skip
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
