"""Post hoc diagnostic of the low-distance plan (docs/notes/S1_mode.md, "Post hoc: cell 1, seed
7"): cell 1 (a 25 to 28 AU), seed 7 against the other seven seeds, from the stored samples only.
No run, rerun or new sampling: each seed's 6,000 reweighted draws are resampled exactly as the
summary resampled them (``s1_mode.stats``: the same rng list, ``[…, 13]``), so the medians agree
with docs/notes/S1_mode.json, which the script checks before it writes anything.

Writes docs/notes/S1_mode_cell1.json (the numbers the note quotes) and docs/notes/S1_mode_cell1.png.

    uv run python docs/notes/scripts/s1_mode_cell1.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from fork.export import m8_runs, s1_cell0, s1_mode
from fork.inference.reweight import reweight_cell

ROOT = Path(__file__).resolve().parents[3]
NOTES = ROOT / "docs" / "notes"
REGION, ODD = "c1", 7
# Seed 7's 84 % point; L1's merged median is 335.07°.
L1_PLACE_DEG = 333.9
GAP_DEG = 320.0  # inside the empty band between the two modes (about 311° to 330°) in every seed

INK, MUTED, GRID, ORANGE = "#1f1f1f", "#8a8a8a", "#dddddd", "#d9480f"


def seed_draws(s: int, space) -> dict:
    d = s1_mode.load_a(REGION, s, s1_mode.CACHE)
    cr = reweight_cell(d["lnl_b"] - d["lnl_fast"], float(d["ln_z_fast"]), float(d["ln_z_fast_err"]))
    rng = np.random.default_rng([*s1_mode.rng_list(REGION, s), 13])
    _, ix = m8_runs._resampled([cr], [0.0], rng, cr.n)
    lon = np.asarray(d["galle_deg"][ix], dtype=float)
    a = m8_runs._perturber(space, d["phi"][ix])[:, 1]
    return {
        "lon": lon,
        "a": a,
        "ln_z": cr.ln_z,
        "ln_z_err": math.hypot(cr.ln_z_fast_err, cr.ln_correction_err),
    }


def q(x: np.ndarray) -> list[float]:
    return [round(float(v), 3) for v in np.percentile(x, [16, 50, 84])]


def main() -> None:
    note = json.loads((NOTES / "S1_mode.json").read_text(encoding="utf-8"))
    published = note["part_a"]["regions"][REGION]["seeds"]
    space = s1_cell0.problem(None)[2]
    seeds = {s: seed_draws(s, space) for s in s1_mode.SEEDS}
    for s, d in seeds.items():  # the same draws as the summary's
        want = published[f"s{s}"]
        if (
            abs(float(np.median(d["lon"])) - want["galle_median"]) > 0.002
            or abs(d["ln_z"] - want["ln_z"]) > 1e-9
        ):
            raise SystemExit(f"seed {s}: draws do not reproduce docs/notes/S1_mode.json")

    rows = {}
    for s, d in seeds.items():
        lon = d["lon"]
        # The evidence split between the two modes at the empty band (GAP_DEG):
        # ln Z + ln f and ln Z + ln(1 − f), with f the share of the draws beyond it.
        f = round(float(np.mean(lon >= GAP_DEG)), 4)
        rows[f"s{s}"] = {
            "lon_deg": q(lon), "a_au": q(d["a"]),
            "frac_at_or_beyond_l1_place": round(float(np.mean(lon >= L1_PLACE_DEG)), 4),
            "frac_far_mode": f,
            "ln_z": round(d["ln_z"], 3), "ln_z_err": round(d["ln_z_err"], 3),
            "ln_z_far": round(d["ln_z"] + math.log(f), 3),
            "ln_z_rest": round(d["ln_z"] + math.log(1.0 - f), 3),
        }  # fmt: skip

    def span(key: str) -> float:
        v = [r[key] for r in rows.values()]
        return round(max(v) - min(v), 3)

    far = [d["a"][d["lon"] >= GAP_DEG] for d in seeds.values()]
    rest = [d["a"][d["lon"] < GAP_DEG] for d in seeds.values()]
    others = [r["frac_at_or_beyond_l1_place"] for k, r in rows.items() if k != f"s{ODD}"]
    out = {
        "description": (
            "Post hoc diagnostic of docs/notes/S1_mode.md: cell 1, seed 7 against the other seven "
            "seeds, from the stored reweighted draws (no new sampling). Quantiles are 16, 50, "
            "84 %. ln_z_err combines UltraNest's error and the reweighting's. The far mode is the "
            "draws beyond gap_deg; ln_z_far and ln_z_rest split each seed's ln Z between the modes."
        ),
        "post_hoc": True,
        "region": REGION, "odd_seed": ODD, "l1_place_deg": L1_PLACE_DEG, "gap_deg": GAP_DEG,
        "seeds": rows,
        "far_mode_a_au": q(np.concatenate(far)),
        "rest_a_au": q(np.concatenate(rest)),
        "range_ln_z": span("ln_z"), "range_ln_z_far": span("ln_z_far"),
        "range_ln_z_rest": span("ln_z_rest"),
        "others_frac_at_or_beyond_l1_place": [min(others), max(others)],
    }  # fmt: skip
    (NOTES / "S1_mode_cell1.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    plot(seeds)
    print(json.dumps(out, indent=1))


def plot(seeds: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_l, ax_a, ax_z) = plt.subplots(1, 3, figsize=(10.5, 3.3),
                                           gridspec_kw={"width_ratios": [1.5, 1, 0.8]})  # fmt: skip
    bins_l = np.arange(300.0, 345.01, 0.5)
    bins_a = np.linspace(25.0, 28.0, 61)
    for s, d in seeds.items():
        odd = s == ODD
        kw = {"color": ORANGE if odd else MUTED, "linewidth": 1.6 if odd else 0.8,
              "linestyle": "-" if odd else "--", "zorder": 3 if odd else 2}  # fmt: skip
        for ax, x, bins in ((ax_l, d["lon"], bins_l), (ax_a, d["a"], bins_a)):
            h, e = np.histogram(x, bins=bins, density=True)
            ax.step(e[:-1], h, where="post", **kw)
    ax_l.axvline(L1_PLACE_DEG, color=INK, linewidth=0.8, linestyle=":")
    top = ax_l.get_ylim()[1] * 0.92
    ax_l.text(L1_PLACE_DEG + 0.4, top, f"{L1_PLACE_DEG}°", color=INK, fontsize=8)
    ax_l.set_xlabel("longitude on Galle's night (°, geocentric)", fontsize=9, color=INK)
    ax_l.set_ylabel("density", fontsize=9, color=INK)
    ax_a.set_xlabel("semi-major axis (AU)", fontsize=9, color=INK)

    ys = list(seeds)
    for s in ys:
        d = seeds[s]
        odd = s == ODD
        ax_z.errorbar(d["ln_z"], s, xerr=d["ln_z_err"], fmt="o", markersize=4,
                      color=ORANGE if odd else INK, ecolor=ORANGE if odd else INK, elinewidth=1,
                      capsize=0)  # fmt: skip
        share = float(np.mean(d["lon"] >= L1_PLACE_DEG))
        ax_z.annotate(f"{100 * share:.1f} %", (d["ln_z"] + d["ln_z_err"], s), xytext=(4, 0),
                      textcoords="offset points", va="center", fontsize=7,
                      color=ORANGE if odd else INK)  # fmt: skip
    ax_z.set_yticks(ys, [f"seed {s}" for s in ys], fontsize=8)
    ax_z.invert_yaxis()
    ax_z.set_xlabel("ln Z (± sampler error); share at ≥ 333.9°", fontsize=9, color=INK)

    for ax in (ax_l, ax_a, ax_z):
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=INK, labelsize=8)
    ax_l.plot([], [], color=ORANGE, linewidth=1.6, label=f"seed {ODD}")
    ax_l.plot([], [], color=MUTED, linewidth=0.8, linestyle="--", label="seeds 1–6, 8")
    ax_l.legend(frameon=False, fontsize=8, loc="upper center")
    fig.tight_layout()
    fig.savefig(NOTES / "S1_mode_cell1.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
