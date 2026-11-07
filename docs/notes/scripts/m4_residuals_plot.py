"""A2 / A2b diagnostic → docs/notes/M4_residuals.png.

    uv run python docs/notes/scripts/m4_residuals_plot.py

Top: noise-free residual curves (Horizons truth − best fit), yearly 1690–1845:
  - A2: no perturber, a/e/ϖ/λ₀ free, i and Ω held at epoch values (docs/RULES.md rule 2), with the noisy
    observation residuals as dots;
  - diagnostic: no perturber, all six Uranus elements free;
  - A2b: true Neptune planted, Uranus refitted.
Bottom: slope of each curve (″/yr); post-1781 sign changes of the A2 slope are marked.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
from fork.data.synthetic import jd_tdb_jan1, make_synthetic_observations, true_uranus_lon_ecl_j2000_rad
from fork.inference.fit import fit_uranus_elements, wrap_rad
from fork.model import bodies
from fork.model.forward import (
    ForwardContext,
    perturber_theta_from_state,
    predict_uranus_lon,
    uranus_inc_raan_at_epoch,
)

ARC = 206264.80624709636
EPOCH = 2371628.5
OUT = Path(__file__).resolve().parents[1] / "M4_residuals.png"
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
# categorical slots 1–3 (validated; slot 3 is < 3:1 contrast, so every line is direct-labelled)
C_A2, C_FREE, C_A2B = "#2a78d6", "#eb6834", "#1baf7a"


def jd_to_year(jd):
    return 2000.0 + (np.asarray(jd) - 2451545.0) / 365.25


def main() -> None:
    ctx = ForwardContext.from_horizons(EPOCH)
    obs, _truth = make_synthetic_observations(np.random.default_rng(0))
    inc_raan = uranus_inc_raan_at_epoch(EPOCH)
    th_nep = perturber_theta_from_state(about("neptune", [EPOCH])[0], bodies.NEPTUNE.mass_msun)
    fits = {
        "a2": (fit_uranus_elements(obs, ctx, inc_raan), None),
        "free": (fit_uranus_elements(obs, ctx, inc_raan, free_inc_raan=True), None),
        "a2b": (fit_uranus_elements(obs, ctx, inc_raan, fixed_perturber_theta=th_nep), th_nep),
    }
    years = np.arange(1690, 1846)
    grid = jd_tdb_jan1(years)
    truth_grid = true_uranus_lon_ecl_j2000_rad(grid)
    curves = {}
    for k, (f, pert) in fits.items():
        theta = f.theta if pert is None else np.concatenate([f.theta, pert])
        curves[k] = wrap_rad(truth_grid - predict_uranus_lon(theta, grid, ctx)) * ARC

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(9.5, 6.6), dpi=150, sharex=True, facecolor=SURFACE,
                                  gridspec_kw={"height_ratios": [2.2, 1]})  # fmt: skip
    for a in (ax, ax2):
        a.set_facecolor(SURFACE)
        a.grid(axis="y", color=GRID, lw=0.8)
        a.spines[["top", "right"]].set_visible(False)
        a.axvspan(1815, 1830, color="#f1f0ec", zorder=0)
        a.axvline(1781, color=INK_MUTED, lw=0.8, ls=":")
        a.axhline(0, color=INK_MUTED, lw=0.8)

    a2 = fits["a2"][0]
    i_f, o_f = np.degrees(fits["free"][0].theta[2]), np.degrees(fits["free"][0].theta[3])
    ax.errorbar(jd_to_year(obs.jd_tdb), a2.residual_rad * ARC, yerr=obs.sigma_rad * ARC, fmt="o",
                ms=3.5, color=C_A2, alpha=0.4, elinewidth=0.8, capsize=0,
                label="observations − A2 fit (σ = 12″ / 3″)")  # fmt: skip
    series = (
        ("a2", C_A2, "-", "A2: no perturber, i, Ω held at epoch values", "A2"),
        ("free", C_FREE, "--",
         f"diagnostic: no perturber, all 6 free (i = {i_f:.2f}°, Ω = {o_f:.0f}°)", "all 6 free"),
        ("a2b", C_A2B, "-", "A2b: true Neptune planted, Uranus refitted", "A2b"),
    )  # fmt: skip
    for k, c, ls, label, short in series:
        ax.plot(years, curves[k], color=c, lw=2, ls=ls, label=label)
        ax.annotate(short, (years[-1], curves[k][-1]), xytext=(5, 0), textcoords="offset points",
                    va="center", fontsize=8.5, color=INK_2)  # fmt: skip
        ax2.plot(years[1:] - 0.5, np.diff(curves[k]), color=c, lw=2, ls=ls)

    post = years >= 1781
    s = np.sign(np.diff(curves["a2"][post]))
    flips = years[post][1:-1][np.nonzero(s[1:] != s[:-1])[0]] + 0.5
    for yr in flips:
        ax2.plot(yr, 0, "o", ms=7, mfc=SURFACE, mec=C_A2, mew=2, zorder=5)
        ax2.annotate(f"{yr:.1f}", (yr, 0), xytext=(0, 8), textcoords="offset points",
                     ha="center", fontsize=8, color=INK_2)  # fmt: skip

    ax.set_ylabel("truth − fit [″]")
    ax.set_title("Uranus residuals: without Neptune (A2) and with the true Neptune planted (A2b)",
                 loc="left", color=INK, fontsize=11)  # fmt: skip
    ax.legend(loc="lower left", frameon=False, fontsize=8.5, labelcolor=INK)
    ax.annotate("1815–1830", (1822.5, ax.get_ylim()[1]), xytext=(0, -12), textcoords="offset points",
                ha="center", fontsize=8.5, color=INK_2)  # fmt: skip
    ax2.set_ylabel("slope [″/yr]")
    ax2.set_xlabel("year  (dotted line: 1781; circles: post-1781 sign changes of the A2 slope)")
    ax2.set_xlim(1688, 1856)
    fig.text(0.01, 0.005, "Synthetic data (seed 0), dt = 12 d; truth = JPL Horizons (includes "
             "Neptune); fits never see Horizons in-plane Uranus elements", fontsize=7.5,
             color=INK_MUTED)  # fmt: skip
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")
    print("post-1781 A2 slope sign changes:", flips.tolist(), "p2p", round(float(np.ptp(curves["a2"])), 1))


if __name__ == "__main__":
    main()
