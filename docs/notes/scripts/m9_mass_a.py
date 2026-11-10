"""After-the-fact diagnostic (not pre-registered; listed in docs/notes/M9.md): the mass–a corner
plot of the M6c replacement posterior.

    uv run python docs/notes/scripts/m9_mass_a.py                    # -> docs/notes/M9_mass_a.png
    uv run python docs/notes/scripts/m9_mass_a.py --variant s3_t4    # -> docs/notes/M10_mass_a.png

Joint density of the perturber's mass and semi-major axis with the contours enclosing 68 % and
95 % of the posterior mass, each marginal beside it, and two marked points: Neptune (JPL;
17.15 M⊕, a = 30.07 AU) and Le Verrier's planet (36.6 M⊕, 36.15 AU; ``M6c_r2_ridge.json``'s
sources). Weighted samples from every cell via ``PartitionedResult``. With ``--variant s3_t4``
the same plot of S3's posterior (M10.md's after-the-fact diagnostic).
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT_PNG = ROOT / "docs" / "notes" / "M9_mass_a.png"
RIDGE_JSON = ROOT / "docs" / "notes" / "M6c_r2_ridge.json"
NEPTUNE = {"mass_earth": 17.15, "a_au": 30.07}
INK, MUTED, BLUE = "#1f1f1f", "#777777", "#2a6fdb"


def contour_levels(h: np.ndarray, fractions=(0.68, 0.95)) -> list[float]:
    flat = np.sort(h.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    return [float(flat[np.searchsorted(cum, f)]) for f in fractions]


def main(variant: str = "replacement") -> None:
    from fork.inference.nested import weighted_quantiles

    if variant == "replacement":
        from fork.export.m9 import replacement_posterior as posterior

        out_png, which = OUT_PNG, "the headline (replacement) posterior"
    else:
        from fork.export.m10 import s3_posterior as posterior

        out_png, which = OUT_PNG.with_name("M10_mass_a.png"), "S3's Student-t posterior"
    q, w, _ = posterior()
    a, m = np.exp(q[:, 5]), 10.0 ** q[:, 4]
    lv = json.loads(RIDGE_JSON.read_text())["le_verrier"]
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.0, 7.0),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [3, 1], "height_ratios": [1, 3]},
    )
    ax_j, ax_m, ax_a, ax_blank = axes[1, 0], axes[1, 1], axes[0, 0], axes[0, 1]
    ax_blank.axis("off")
    a_edges = np.linspace(20.0, 46.0, 105)
    m_edges = np.linspace(5.0, 60.0, 111)
    h, _, _ = np.histogram2d(a, m, bins=[a_edges, m_edges], weights=w)
    ac, mc = 0.5 * (a_edges[1:] + a_edges[:-1]), 0.5 * (m_edges[1:] + m_edges[:-1])
    levels = sorted(contour_levels(h))
    ax_j.contourf(ac, mc, h.T, levels=[levels[0], levels[1], h.max() * 1.01],
                  colors=[BLUE], alpha=0.35, antialiased=True)  # fmt: skip
    ax_j.contourf(ac, mc, h.T, levels=[levels[1], h.max() * 1.01], colors=[BLUE], alpha=0.5)
    ax_j.contour(ac, mc, h.T, levels=levels, colors=[BLUE], linewidths=1.0)
    for pt, label, mk in (
        (NEPTUNE, "Neptune (JPL)", "o"),
        ({"mass_earth": lv["mass_earth"], "a_au": lv["a_au"]}, "Le Verrier's planet", "s"),
    ):
        ax_j.plot(pt["a_au"], pt["mass_earth"], marker=mk, markersize=8, markerfacecolor="none",
                  markeredgecolor=INK, linestyle="none")  # fmt: skip
        ax_j.annotate(label, (pt["a_au"], pt["mass_earth"]), xytext=(8, 6),
                      textcoords="offset points", fontsize=9, color=INK)  # fmt: skip
    ax_a.hist(a, bins=a_edges, weights=w, color=BLUE, alpha=0.6)
    ax_m.hist(m, bins=m_edges, weights=w, color=BLUE, alpha=0.6, orientation="horizontal")
    for v in weighted_quantiles(a, w, (0.16, 0.5, 0.84)):
        ax_a.axvline(v, color=MUTED, linewidth=0.8)
    for v in weighted_quantiles(m, w, (0.16, 0.5, 0.84)):
        ax_m.axhline(v, color=MUTED, linewidth=0.8)
    ax_j.set_xlabel("semi-major axis a (AU, osculating 1781)")
    ax_j.set_ylabel("mass (M⊕)")
    ax_j.set_xlim(20, 46)
    ax_j.set_ylim(5, 60)
    ax_a.set_xlim(20, 46)
    ax_m.set_ylim(5, 60)
    ax_a.set_yticks([])
    ax_m.set_xticks([])
    for ax in (ax_j, ax_a, ax_m):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle(f"Mass against a in {which}: 68 % and 95 % contours; made after the fact",
                 fontsize=9)  # fmt: skip
    fig.savefig(out_png, dpi=160)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=("replacement", "s3_t4"), default="replacement")
    main(ap.parse_args().variant)
