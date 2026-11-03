"""A1 energy-comparison figure → docs/notes/M3_energy.png.

    uv run python docs/notes/scripts/m3_energy_plot.py

Sun + Jupiter + Saturn + Uranus from Horizons at 1781-03-13, default dt (14 d), 1000 yr. Plots the
envelope (max over 25-yr windows) of |ΔE/E| on one log axis. This becomes a site panel in M7.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fork.integrators import DEFAULT_DT_DAYS, leapfrog, rk4, wisdom_holman
from fork.integrators._nbody import energy_helio
from fork.model import bodies

ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "tests" / "data" / "horizons_a1_1781-03-13.json"
OUT = ROOT / "docs" / "notes" / "M3_energy.png"
DT = DEFAULT_DT_DAYS

# Reference palette (dataviz skill), light mode; validated: CVD ΔE ≥ 9.1, normal ≥ 22.9.
# Aqua and yellow are < 3:1 on the surface → every line is also direct-labeled.
SURFACE, INK, INK_2, INK_MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984"
SERIES = [  # fixed categorical order: slot 1..4
    ("wh_drift_outer", "Wisdom–Holman (default)", "#2a78d6"),
    ("rk4", "RK4", "#eb6834"),
    ("wh_kick_outer", "Wisdom–Holman (kick outermost)", "#1baf7a"),
    ("leapfrog", "Leapfrog (KDK)", "#eda100"),
]


def main() -> None:
    data = json.loads(FIXTURE.read_text())
    names = ("jupiter", "saturn", "uranus")
    st = np.array([data["state_ecl_j2000"][n] for n in names])
    x0, v0 = st[:, :3].copy(), st[:, 3:].copy()
    masses = np.array([bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in names])
    t = np.arange(0.0, 1000 * 365.25, DT)
    e0 = energy_helio(x0, v0, masses)
    runs = {
        "wh_drift_outer": wisdom_holman.integrate_states(x0, v0, masses, t, DT),
        "wh_kick_outer": wisdom_holman.integrate_states(x0, v0, masses, t, DT, 0),
        "rk4": rk4.integrate_states(x0, v0, masses, t, DT),
        "leapfrog": leapfrog.integrate_states(x0, v0, masses, t, DT),
    }
    win = round(25 * 365.25 / DT)
    n_win = t.size // win
    t_yr = (np.arange(n_win) + 0.5) * win * DT / 365.25

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2})  # fmt: skip
    fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for key, label, color in SERIES:
        xs, vs = runs[key]
        de = np.abs([(energy_helio(xs[k], vs[k], masses) - e0) / abs(e0) for k in range(t.size)])
        env = de[: n_win * win].reshape(n_win, win).max(axis=1)
        ax.plot(t_yr, env, color=color, lw=2.0, label=label, solid_capstyle="round")
        ax.annotate(label, (t_yr[-1], env[-1]), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=9, color=INK)  # fmt: skip
    ax.axhline(1e-8, color=INK, lw=1.2, ls=(0, (4, 3)), zorder=5, label="A1 bound for WH (1e-8)")
    ax.set_yscale("log")
    ax.set_xlim(0, 1000)
    ax.set_xlabel("years after 1781-03-13")
    ax.set_ylabel("|ΔE / E|  (max per 25-yr window)")
    ax.set_title("Symplectic integrators keep energy error bounded; RK4 drifts",
                 loc="left", color=INK, fontsize=11.5)  # fmt: skip
    ax.grid(axis="y", which="major", color="#e6e5e1", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), ncol=2, frameon=False,
              fontsize=9, labelcolor=INK)  # fmt: skip
    fig.text(0.01, 0.005, f"Sun + Jupiter + Saturn + Uranus, Horizons ICs at 1781-03-13, dt = {DT:g} d",
             fontsize=8, color=INK_MUTED)  # fmt: skip
    fig.subplots_adjust(left=0.1, right=0.66, top=0.9, bottom=0.3)
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
