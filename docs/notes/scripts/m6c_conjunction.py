"""After-the-fact diagnostic (not pre-registered): the hidden planet at the 1821–22 conjunction.

    uv run python docs/notes/scripts/m6c_conjunction.py compute   # integrations -> JSON (~10 min)
    uv run python docs/notes/scripts/m6c_conjunction.py plot      # JSON -> PNG (seconds)
    uv run python docs/notes/scripts/m6c_conjunction.py           # both

Uranus overtook the hidden planet around 1821–22, when the pull between them was strongest.
This script asks how well the M6c replacement posterior (the headline's cells 0, 1 and 3 with
the merged R2 cell 2; M6.md) fixes the planet's heliocentric longitude at that conjunction,
compared with its geocentric longitude on Galle's night (1846 Sep 23).

**The epoch, from the posterior alone.** For each derived sample, the conjunction is the
moment its modelled Uranus and its planet have equal heliocentric ecliptic longitude between
1810 and 1835 (a 10-day grid, interpolated linearly). The epoch used is the median of those
moments, rounded to the nearest 0h TDB. It is not taken from Neptune, so the fit's quantities
stay free of Neptune's ephemeris. JPL's own Uranus–Neptune conjunction is computed the same way
from Horizons and reported beside it, as a comparison target only.

At that fixed epoch, for every derived sample: the planet's heliocentric longitude, ecliptic of
date, about the inner barycentre (as for the 1847 longitude in ``m6c.summarize_results``).
Alongside: the Galle-night geocentric longitude, and the heliocentric one on that night.

Samples: the 6000 derived samples of the replacement posterior, drawn exactly as the reported
interval was (``m6c_repeat.REPORT_RNG``); their Galle-night quantiles are checked against it.
Groups, as in ``m6c_r2_ridge.py``: cell 2 (a 28–60 AU, e < 0.12) split into core (e < 0.06),
inner branch (e ≥ 0.06, a < 30.5) and outer branch (e ≥ 0.06, a ≥ 30.5); the rest (cells 0, 1
and 3) as "other cells".

Overlays: JPL's Neptune (Horizons) and Le Verrier's 1846 elements
(``fork.model.hidden_mass.LE_VERRIER_1846``; Krajnović 2016, A&G 57, 5.28, Table 2, after
Grosser 1962; Encyclopædia Britannica 1911), propagated as a Kepler ellipse with
``m6c.le_verrier_state_ecl_j2000``, as for its Galle-night longitude.

Output: site/data/m6c_conjunction.json (``compute``) and docs/notes/M6c_conjunction.png
(``plot``, drawn from that JSON alone, so layout changes never re-run the integrations).
"""

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
from fork.export import m6c, m6c_repeat
from fork.export.to_json import SITE_DATA
from fork.inference.mcmc import circular_interval_deg
from fork.inference.nested import ERR_INFLATION, PartitionedResult
from fork.inference.profile import GalleGeometry, galle_lon_of_date_rad
from fork.model.forward import predict_states_ecl_j2000
from fork.model.hidden_mass import I_PERTURBER, LE_VERRIER_1846
from fork.model.observation import I_URANUS

ROOT = Path(__file__).resolve().parents[3]
OUT_JSON = SITE_DATA / "m6c_conjunction.json"
OUT_PNG = ROOT / "docs" / "notes" / "M6c_conjunction.png"
JD_1810, JD_1835 = 2382148.5, 2391279.5  # 1810 Jan 1.0 and 1835 Jan 1.0 TDB
GRID_STEP_DAYS = 10.0
E_RIDGE, A_SPLIT = 0.06, 30.5
GROUPS = ("core", "inner", "outer", "other")
LABELS = {"core": "core (e < 0.06)", "inner": "inner branch (a ≈ 28.8 AU)",
          "outer": "outer branch (a ≈ 32.1 AU)", "other": "other cells"}  # fmt: skip
BIN_DEG = 0.25

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
# categorical slots 1–3 (validated all-pairs, light mode); "other cells" is not drawn as a shape
COLOR = {"core": "#2a78d6", "inner": "#eb6834", "outer": "#1baf7a"}


def wrap180(x):
    return (np.asarray(x) + 180.0) % 360.0 - 180.0


def lon_of_date_deg(pos_ecl_j2000: np.ndarray, jd: np.ndarray) -> np.ndarray:
    """Ecliptic-of-date longitude, degrees, of positions (n, 3) at epochs jd (n,)."""
    out = np.empty(len(jd))
    for i, (r, t) in enumerate(zip(pos_ecl_j2000, jd, strict=True)):
        v = precession_matrix_ecl_j2000_to_ecl_of_date(t) @ r
        out[i] = math.degrees(math.atan2(v[1], v[0])) % 360.0
    return out


def conjunction_jd(jd: np.ndarray, lon_uranus: np.ndarray, lon_planet: np.ndarray) -> float:
    """First epoch where Uranus overtakes the planet in heliocentric longitude (NaN if none)."""
    d = wrap180(lon_uranus - lon_planet)
    k = np.flatnonzero((d[:-1] < 0) & (d[1:] >= 0) & (np.abs(d[:-1]) < 90))
    if len(k) == 0:
        return float("nan")
    i = k[0]
    return float(jd[i] + (jd[i + 1] - jd[i]) * (-d[i]) / (d[i + 1] - d[i]))


def jd_to_date(jd: float) -> str:
    from astropy.time import Time

    return Time(jd, format="jd", scale="tdb").iso[:16]


def q(x: np.ndarray) -> dict:
    lo95, lo, med, hi, hi95 = circular_interval_deg(x, (2.5, 16.0, 50.0, 84.0, 97.5))
    r = lambda v: float(wrap180(v - med))  # noqa: E731
    return {"median": round(med, 3), "lo68": round(lo, 3), "hi68": round(hi, 3),
            "lo95": round(lo95, 3), "hi95": round(hi95, 3),
            "width68": round(r(hi) - r(lo), 3), "width95": round(r(hi95) - r(lo95), 3)}  # fmt: skip


def place(target: float, s: dict) -> dict:
    """target − median in degrees and in 68 % half-widths on the facing side; inside 68 / 95?"""
    off = float(wrap180(target - s["median"]))
    half = (s["hi68"] - s["median"]) % 360.0 if off > 0 else (s["median"] - s["lo68"]) % 360.0
    in68 = float(wrap180(s["lo68"] - s["median"])) <= off <= float(wrap180(s["hi68"] - s["median"]))
    in95 = float(wrap180(s["lo95"] - s["median"])) <= off <= float(wrap180(s["hi95"] - s["median"]))
    return {"value": round(target, 3), "minus_median_deg": round(off, 3),
            "in_half_widths": round(off / half, 2),
            "in_68": bool(in68), "in_95": bool(in95)}  # fmt: skip


def compute() -> None:
    m6c.check_frozen("headline")
    rep = json.loads((SITE_DATA / m6c_repeat.OUT_NAME).read_text())
    assert rep["decision"]["outcome"] == "replace"
    _, ctx, *_ = m6c.dataset(variant="headline")
    geo = GalleGeometry.from_horizons()

    # the replacement posterior's derived samples, drawn as m6c.summarize_results draws them
    res = m6c_repeat.replacement_cells()
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng(list(m6c_repeat.REPORT_RNG))
    eq = p.equal_weight(rng)
    idx = rng.choice(len(eq), size=min(m6c.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    th = np.array([res[0].space.theta(v) for v in eq[idx]])
    a, e = th[:, 7], th[:, 8]
    cell2 = (a >= 28.0) & (a < 60.0) & (e < 0.12)
    group = np.where(~cell2, "other", np.where(e < E_RIDGE, "core",
                     np.where(a < A_SPLIT, "inner", "outer")))  # fmt: skip

    # 1. per-sample conjunction epochs, 1810–1835
    grid = np.arange(JD_1810, JD_1835 + 1, GRID_STEP_DAYS)
    conj = np.empty(len(th))
    for i, t in enumerate(th):
        pos, _ = predict_states_ecl_j2000(t, grid, ctx)
        conj[i] = conjunction_jd(grid, lon_of_date_deg(pos[:, I_URANUS], grid),
                                 lon_of_date_deg(pos[:, I_PERTURBER], grid))  # fmt: skip
        if i % 1000 == 0:
            print(f"conjunction {i}/{len(th)}", flush=True)
    jd_c = float(np.round(np.nanmedian(conj) - 0.5) + 0.5)  # nearest 0h TDB

    jpl_u = lon_of_date_deg(about("uranus", grid)[:, :3], grid)
    jpl_n = lon_of_date_deg(about("neptune", grid)[:, :3], grid)
    jd_c_jpl = conjunction_jd(grid, jpl_u, jpl_n)

    # 2. longitudes at the fixed conjunction epoch and on Galle's night
    lon_c, galle, helio_g = (np.empty(len(th)) for _ in range(3))
    for i, t in enumerate(th):
        pos, _ = predict_states_ecl_j2000(t, [jd_c, geo.jd_tdb], ctx)
        lon_c[i], helio_g[i] = lon_of_date_deg(pos[:, I_PERTURBER], np.array([jd_c, geo.jd_tdb]))
        galle[i] = math.degrees(galle_lon_of_date_rad(t, ctx, geo)) % 360.0
        if i % 1000 == 0:
            print(f"longitudes {i}/{len(th)}", flush=True)
    reported = rep["reported"]
    sg = q(galle)
    assert all(abs(sg[k] - reported[k]) < 1e-3 for k in ("median", "lo68", "hi68")), sg

    lv_c = lon_of_date_deg(m6c.le_verrier_state_ecl_j2000(jd_c)[0][None], np.array([jd_c]))[0]
    lv_g = m6c.le_verrier_galle_lon_deg()
    jpl_c = lon_of_date_deg(about("neptune", [jd_c])[:, :3], np.array([jd_c]))[0]

    epochs = {
        "conjunction_helio": {"values": lon_c, "jpl": jpl_c, "le_verrier": lv_c},
        "galle_night_geo": {"values": galle, "jpl": m6c.JPL_GALLE_LON_DEG, "le_verrier": lv_g},
    }
    summary = {}
    for name, ep in epochs.items():
        allq = q(ep["values"])
        summary[name] = {
            "all": allq,
            "groups": {g: {"share": round(float(np.mean(group == g)), 4)}
                          | q(ep["values"][group == g]) for g in GROUPS},
            "jpl_neptune": place(ep["jpl"], allq),
            "le_verrier": place(ep["le_verrier"], allq),
        }  # fmt: skip
    summary["galle_night_helio"] = {"all": q(helio_g)}

    # histograms for the figure / site (per group, counts on BIN_DEG bins about each median)
    hist = {}
    for name in ("conjunction_helio", "galle_night_geo"):
        c = summary[name]["all"]["median"]
        rel = wrap180(epochs[name]["values"] - c)
        edges = np.arange(-12.0, 12.0 + BIN_DEG, BIN_DEG)
        hist[name] = {"centre_deg": c, "edges_rel_deg": edges.round(3).tolist(),
                      "counts": {g: np.histogram(rel[group == g], edges)[0].tolist()
                                 for g in GROUPS}}  # fmt: skip

    payload = {
        "description": __doc__.strip().splitlines()[0],
        "posterior": "M6c replacement: headline cells 0, 1, 3 + merged R2 cell 2 "
                     "(site/data/m6c_replacement.json)",
        "n_samples": len(th),
        "epoch": {
            "jd_tdb": jd_c, "date_tdb": jd_to_date(jd_c),
            "rule": "median over the derived samples of each one's Uranus-planet heliocentric "
                    "conjunction (1810-1835, 10-day grid), rounded to 0h TDB; no Neptune input",
            "sample_conjunction_jd_16_50_84": [round(float(v), 1) for v in
                                               np.nanpercentile(conj, [16, 50, 84])],
            "sample_conjunction_years_68_width": round(float(np.diff(
                np.nanpercentile(conj, [16, 84]))[0]) / 365.25, 2),
            "n_without_conjunction": int(np.isnan(conj).sum()),
            "jpl_uranus_neptune_conjunction_jd_tdb": round(jd_c_jpl, 1),
            "jpl_uranus_neptune_conjunction_date_tdb": jd_to_date(jd_c_jpl),
        },
        "galle_night_jd_tdb": geo.jd_tdb,
        "frames": {"conjunction_helio": "heliocentric (inner barycentre), ecliptic of date",
                   "galle_night_geo": "geocentric apparent, ecliptic of date (as reported)",
                   "galle_night_helio": "heliocentric (inner barycentre), ecliptic of date"},
        "groups": {g: LABELS[g] for g in GROUPS},
        "le_verrier_sources": list(LE_VERRIER_1846.sources),
        **summary,
        "histograms": hist,
    }  # fmt: skip
    OUT_JSON.write_text(json.dumps(payload, indent=1) + "\n")
    print(json.dumps({k: payload[k] for k in ("epoch",)}, indent=1))
    for name in ("conjunction_helio", "galle_night_geo", "galle_night_helio"):
        print(name, json.dumps(summary[name]))
    print(f"wrote {OUT_JSON}")


def plot() -> None:
    """Render the figure from the saved JSON only: no posterior, no integrations."""
    payload = json.loads(OUT_JSON.read_text())
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                         "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2})  # fmt: skip
    # y not shared: each panel is a density per degree on its own scale; the x spans are equal
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), facecolor=SURFACE)
    titles = {"conjunction_helio": f"At conjunction, {payload['epoch']['date_tdb'][:10]} "
                                   "(heliocentric)",
              "galle_night_geo": "Galle's night, 1846-09-23 (geocentric)"}  # fmt: skip
    for ax, name in zip(axes, ("conjunction_helio", "galle_night_geo"), strict=True):
        h, s = payload["histograms"][name], payload[name]
        c, edges = h["centre_deg"], np.array(h["edges_rel_deg"])
        ax.set_facecolor(SURFACE)
        ymax = 0.0
        for g in ("core", "inner", "outer"):
            n = np.array(h["counts"][g], dtype=float)
            dens = n / (n.sum() * BIN_DEG)  # each group normalised to unit area: shapes compare
            ymax = max(ymax, dens.max())
            nz = np.flatnonzero(n)  # draw only the occupied span, not a flat line at zero
            i0, i1 = max(nz[0] - 1, 0), min(nz[-1] + 2, len(n))
            ax.stairs(dens[i0:i1], c + edges[i0 : i1 + 1], color=COLOR[g], linewidth=2,
                      label=f"{LABELS[g]}, {s['groups'][g]['share']:.1%}")  # fmt: skip
        # the whole posterior's 68 % interval and median, as a bar under the shapes
        lo, hi = c + wrap180(s["all"]["lo68"] - c), c + wrap180(s["all"]["hi68"] - c)
        yb = -0.06 * ymax
        ax.plot([lo, hi], [yb, yb], color=INK, linewidth=3, solid_capstyle="butt")
        ax.plot([s["all"]["median"]], [yb], "o", color=INK, markersize=5)
        ax.annotate(f"all, 68 %: {s['all']['width68']:.2f}°", (lo, yb), xytext=(-5, 0),
                    textcoords="offset points", color=INK_2, fontsize=8, va="center",
                    ha="right")  # fmt: skip
        for key, ls, lab in (("jpl_neptune", "-", "JPL Neptune"),
                             ("le_verrier", "--", "Le Verrier")):  # fmt: skip
            x = c + float(wrap180(s[key]["value"] - c))
            ax.axvline(x, color=INK, linestyle=ls, linewidth=1)
            ax.annotate(lab, (x, ymax * 1.42), xytext=(3, 0), textcoords="offset points",
                        color=INK, fontsize=8, rotation=90, va="top")  # fmt: skip
        ax.set_xlim(c - 12, c + 12)  # same 24° span in both panels, so widths compare
        ax.set_ylim(yb * 2.2, ymax * 1.45)  # headroom above the shapes for the line labels
        ax.set_title(titles[name], color=INK, fontsize=10, loc="left")
        ax.set_xlabel("ecliptic longitude of date (°)")
        ax.grid(axis="x", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("density per degree (each group integrates to 1)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, loc="lower center", ncol=3,
               labelcolor=INK_2, bbox_to_anchor=(0.5, 0.0))  # fmt: skip
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(OUT_PNG, dpi=150, facecolor=SURFACE)
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else "all"
    if step not in ("compute", "plot", "all"):
        sys.exit(f"usage: {sys.argv[0]} [compute|plot|all]")
    if step in ("compute", "all"):
        compute()
    if step in ("plot", "all"):
        plot()
