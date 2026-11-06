"""Seed 5 (free a, e ≤ 0.3): why the posterior median sits −1.74σ from truth → M5b_seed5_lon.png.

    uv run python docs/notes/scripts/m5b_seed5_diagnostic.py

Reads the M5b cell results (data/cache/nested/m5b/seed_5__free_a_e0.3__c*.pkl). Draws N
combined-posterior samples, split over cells in proportion to their evidence shares, and
computes each draw's 1846-09-23 longitude. Top: the posterior histogram of that longitude,
stacked by cell group. Bottom: profile Δχ² (min over every other parameter, relative to the
global best fit): free a, e ≤ 0.3 from site/data/m5a_profile.json, and the circular branch
(e ≤ 0.12) computed here every 2° with the same machinery (standard starts + the cell's best
posterior draw near each longitude; ~5 min on 2 workers). Prints the numbers quoted in
docs/notes/M5.md.
"""

import json
import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fork.export.m5b import dataset, jobs_for, truth_params
from fork.export.profile_seeds import cell_draws, seeded_starts
from fork.export.to_json import A1_EPOCH_JD_TDB, SITE_DATA
from fork.inference.fit import le_verrier_space
from fork.inference.nested import ERR_INFLATION, NestedResult, PartitionedResult
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad, profile_galle_longitude
from fork.model.forward import ForwardContext

OUT = Path(__file__).resolve().parents[1] / "M5b_seed5_lon.png"
SEED, RUN, N_DRAWS = 5, "free_a_e0.3", 8000
SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
COLORS = ("#2a78d6", "#eb6834", "#1baf7a")  # categorical slots 1–3, fixed order
GROUPS = (("a 28–60 AU, e < 0.12", (2,)), ("a 28–60 AU, e ≥ 0.12", (3,)), ("a 20–28 AU", (0, 1)))
BINS = np.arange(300.0, 351.0, 1.0)
CIRC_E_MAX = 0.12  # the circular cell's e range
CIRC_GRID_DEG = np.arange(306.0, 337.0, 2.0)
CIRC_CACHE = (
    Path(__file__).resolve().parents[3] / "data" / "cache" / "profile" / "seed_5_circular.json"
)


def pct(x: np.ndarray) -> tuple[float, float, float]:
    lo, med, hi = np.percentile(x, [16, 50, 84])
    return float(lo), float(med), float(hi)


def offset_sigma(med: float, lo: float, hi: float, truth: float) -> float:
    off = med - truth
    return off / (med - lo if off > 0 else hi - med)


def main() -> None:
    ctx = ForwardContext.from_horizons(A1_EPOCH_JD_TDB)
    geo = GalleGeometry.from_horizons()
    truth = truth_params()["galle_lon_ecl_of_date_deg"]
    jobs = [j for j in jobs_for(SEED) if j.run == RUN and j.rep == 1]
    res = [NestedResult.load(j.path(".pkl")) for j in jobs]
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng([5, 2026])
    lon, chi2, cell = [], [], []
    for k, (r, share) in enumerate(zip(res, p.shares, strict=True)):
        n = round(N_DRAWS * float(share))
        if n == 0:
            continue
        w = r.weights / r.weights.sum()
        idx = rng.choice(w.size, size=n, p=w)
        q = r.samples[idx]
        lon += [math.degrees(galle_lon_of_date_rad(r.space.theta(v), ctx, geo)) for v in q]
        chi2 += list(-2.0 * np.asarray(r.results["logl"])[idx])
        cell += [k] * n
    lon, chi2, cell = np.array(lon) % 360.0, np.array(chi2), np.array(cell)

    lo, med, hi = pct(lon)
    print(f"seed {SEED}: N {lon.size}, lon median {med:.2f} [{lo:.2f}, {hi:.2f}], truth {truth}, "
          f"offset {offset_sigma(med, lo, hi, truth):+.2f}σ")  # fmt: skip
    for k, (s, se) in enumerate(zip(p.shares, p.share_errs, strict=True)):
        sel = cell == k
        if sel.sum() < 20:
            print(f"  cell {k}: share {s:.4f} ± {se:.4f} (×{ERR_INFLATION}), {sel.sum()} draws")
            continue
        clo, cmed, chi = pct(lon[sel])
        print(f"  cell {k}: share {s:.4f} ± {se:.4f} (×{ERR_INFLATION}), lon {cmed:.2f} "
              f"[{clo:.2f}, {chi:.2f}], offset {offset_sigma(cmed, clo, chi, truth):+.2f}σ, "
              f"χ² min {chi2[sel].min():.2f}")  # fmt: skip
    without = cell != 3
    wlo, wmed, whi = pct(lon[without])
    print(f"  without the e ≥ 0.12 cell: median {wmed:.2f} [{wlo:.2f}, {whi:.2f}], "
          f"offset {offset_sigma(wmed, wlo, whi, truth):+.2f}σ")  # fmt: skip
    h_all, _ = np.histogram(lon, BINS)
    dens = h_all / lon.size / np.diff(BINS)
    i_med, i_tru = np.searchsorted(BINS, med) - 1, np.searchsorted(BINS, truth) - 1
    print(f"  posterior density: {dens[i_med]:.4f}/° at the median's bin, {dens[i_tru]:.4f}/° at "
          f"truth's bin (ratio {dens[i_med] / dens[i_tru]:.1f})")  # fmt: skip

    prof = next(d for d in json.loads((SITE_DATA / "m5a_profile.json").read_text())["datasets"]
                if d["seed"] == SEED)["modes"]["free_a"]  # fmt: skip
    best_chi2 = prof["overall_min"]["chi2"]
    obs, _, inc_raan, _, _ = dataset(SEED)
    targets = np.radians(np.sort(np.concatenate([CIRC_GRID_DEG, [truth]])))
    if CIRC_CACHE.exists():
        cached = json.loads(CIRC_CACHE.read_text())
        circ_chi2, circ_theta = np.array(cached["chi2"]), np.array(cached["theta"])
    else:
        circ_cell = [c for c in cell_draws(SEED, ctx, geo) if c["cell"] == 2]
        space = le_verrier_space(None, e_max=CIRC_E_MAX)
        circ = profile_galle_longitude(obs, ctx, inc_raan, space, targets,
                                       mode="free_a_circular", n_workers=2, geo=geo, sweeps=False,
                                       seeded_starts=seeded_starts(circ_cell, targets))  # fmt: skip
        circ_chi2, circ_theta = circ.chi2, np.array([f.theta for f in circ.fits])
        CIRC_CACHE.parent.mkdir(parents=True, exist_ok=True)
        CIRC_CACHE.write_text(json.dumps({"lon_deg": np.degrees(targets).tolist(),
                                          "chi2": circ_chi2.tolist(),
                                          "theta": circ_theta.tolist()}))  # fmt: skip
    circ_lon = np.degrees(targets)
    circ_d = circ_chi2 - best_chi2
    k_tru = int(np.argmin(np.abs(circ_lon - truth)))
    k_min = int(np.argmin(circ_chi2))
    print(f"  global best χ² {best_chi2:.3f} at {prof['overall_min']['lon_ecl_of_date_deg']:.2f}°"
          f"; e ≤ 0.3 profile Δχ² at truth {prof['dchi2_at_true_neptune']:.2f} (fit there: "
          f"a {prof['at_true_neptune']['a_au']} AU, e {prof['at_true_neptune']['e']})")  # fmt: skip
    th_t = circ_theta[k_tru]
    print(f"  circular branch (e ≤ {CIRC_E_MAX}): Δχ² at truth {circ_d[k_tru]:.2f} (m "
          f"{10 ** th_t[6]:.1f}, a {th_t[7]:.2f}, e {th_t[8]:.3f}); its minimum Δχ² "
          f"{circ_d[k_min]:.2f} at {circ_lon[k_min]:.1f}°")  # fmt: skip

    ratios = {"height_ratios": [1.3, 1]}
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6.4), dpi=150, facecolor=SURFACE,
                                   sharex=True, gridspec_kw=ratios)  # fmt: skip
    centers = 0.5 * (BINS[1:] + BINS[:-1])
    bottom = np.zeros(centers.size)
    for (label, ks), color in zip(GROUPS, COLORS, strict=True):
        sel = np.isin(cell, ks)
        h, _ = np.histogram(lon[sel], BINS)
        dens = h / lon.size / np.diff(BINS)
        ax1.bar(centers, dens, width=0.92, bottom=bottom, color=color, label=label,
                linewidth=0)  # fmt: skip
        bottom += dens
    p = prof["profile"]
    ax2.plot(p["lon_ecl_of_date_deg"], np.array(p["chi2"]) - best_chi2, color=INK, lw=2,
             label="free a, e ≤ 0.3 (M5a profile)")  # fmt: skip
    ax2.plot(circ_lon, circ_d, "o-", color=COLORS[0], lw=2, ms=4, mec=SURFACE, mew=1,
             label=f"circular branch, e ≤ {CIRC_E_MAX}")  # fmt: skip
    for ax in (ax1, ax2):
        ax.set_facecolor(SURFACE)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.axvline(truth, color=INK_MUTED, lw=1.0, ls=":")
        ax.tick_params(colors=INK_2, labelsize=8)
    top = 1.25 * bottom.max()
    ax1.set_ylim(0, top)
    ax1.axvspan(lo, hi, color=INK_MUTED, alpha=0.12, lw=0)
    ax1.axvline(med, color=INK, lw=1.2)
    ax1.annotate(f"median {med:.1f}°, 68 % [{lo:.1f}, {hi:.1f}]", (lo, top), xytext=(-4, -12),
                 textcoords="offset points", ha="right", fontsize=8, color=INK_2)  # fmt: skip
    ax1.annotate(f"truth {truth:.1f}°", (truth, top), xytext=(4, -12),
                 textcoords="offset points", fontsize=8, color=INK_2)  # fmt: skip
    ax1.set_ylabel("posterior density (per °)", fontsize=9, color=INK_2)
    ax1.legend(frameon=False, fontsize=8, loc="upper right", labelcolor=INK_2)
    ax2.axhline(1.0, color=INK_MUTED, lw=0.8)
    ax2.set_ylim(-0.2, 4.0)
    ax2.set_ylabel("profile Δχ² (vs global best)", fontsize=9, color=INK_2)
    ax2.legend(frameon=False, fontsize=8, loc="upper center", labelcolor=INK_2)
    ax2.set_xlabel("perturber geocentric longitude on 1846-09-23 (ecliptic of date, °)",
                   fontsize=9, color=INK_2)  # fmt: skip
    ax2.set_xlim(BINS[0], BINS[-1])
    fig.suptitle(f"Seed {SEED}, free a, e ≤ 0.3: posterior mass vs best fit at each longitude",
                 fontsize=10, color=INK, x=0.07, ha="left")  # fmt: skip
    fig.tight_layout()
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
