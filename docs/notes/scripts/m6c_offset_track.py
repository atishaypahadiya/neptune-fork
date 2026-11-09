"""After-the-fact diagnostic (not pre-registered): Le Verrier's planet vs the posterior, 1800–1870.

    uv run python docs/notes/scripts/m6c_offset_track.py compute   # integrations -> JSON
    uv run python docs/notes/scripts/m6c_offset_track.py plot      # JSON -> PNG (seconds)
    uv run python docs/notes/scripts/m6c_offset_track.py           # both

Question: would Le Verrier's prediction have succeeded if the search had been made in 1836 or
1856 instead of 1846? His published elements are held fixed and carried to each year; this
does not ask what he would have predicted with only the data up to then.

At each epoch from 1800 to 1870 (281 epochs a quarter-year apart, plus the marked epochs) this
computes, from the M6c replacement posterior (the same 6000 derived samples as
``m6c_conjunction.py``, drawn with ``m6c_repeat.REPORT_RNG``):
- the planet's heliocentric longitude, ecliptic of date, about the inner barycentre: median and
  the 68 % and 95 % bounds, all relative to the median;
- Le Verrier's planet minus the median: his 1846 elements as a Kepler ellipse
  (``m6c.le_verrier_state_ecl_j2000``; ``fork.model.hidden_mass.LE_VERRIER_1846``: Krajnović
  2016, A&G 57, 5.28, Table 2, after Grosser 1962; Encyclopædia Britannica 1911);
- JPL's Neptune (Horizons) minus the median, as a comparison target only.

Marked epochs: the posterior's Uranus–planet conjunction (1821 Oct 3.0 TDB, read from
``site/data/m6c_conjunction.json``), Galle's night (1846 Sep 23), the same date in 1836 and
1856, and every year in which Le Verrier's offset from the median crosses zero.

The conjunction and Galle-night quantiles are checked against ``m6c_conjunction.json``.

Output: site/data/m6c_offset_track.json (``compute``) and docs/notes/M6c_offset_track.png
(``plot``, drawn from that JSON alone).
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SITE_DATA = ROOT / "site" / "data"
OUT_JSON = SITE_DATA / "m6c_offset_track.json"
OUT_PNG = ROOT / "docs" / "notes" / "M6c_offset_track.png"
CONJ_JSON = SITE_DATA / "m6c_conjunction.json"
# per-sample longitudes, so a failure after the integrations never repeats them (gitignored)
LON_CACHE = ROOT / "data" / "cache" / "m6c_offset_track_lon.npz"
JD_1800, JD_1870 = 2378496.5, 2404063.5  # 1800 Jan 1.0 and 1870 Jan 1.0 TDB
N_GRID = 281  # about a quarter-year apart
WHAT_IF = {"search_1836": "1836-09-23 00:00", "search_1856": "1856-09-23 00:00"}

SURFACE, INK, INK_2, INK_MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BAND = "#2a78d6"  # categorical slot 1: the posterior, one hue at two opacities (95 %, 68 %)


def wrap180(x):
    return (np.asarray(x) + 180.0) % 360.0 - 180.0


def compute() -> None:
    from astropy.time import Time

    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.export import m6c, m6c_repeat
    from fork.inference.mcmc import circular_interval_deg
    from fork.inference.nested import ERR_INFLATION, PartitionedResult
    from fork.model.forward import predict_states_ecl_j2000
    from fork.model.hidden_mass import I_PERTURBER, LE_VERRIER_1846

    m6c.check_frozen("headline")
    rep = json.loads((SITE_DATA / m6c_repeat.OUT_NAME).read_text())
    assert rep["decision"]["outcome"] == "replace"
    conj = json.loads(CONJ_JSON.read_text())
    _, ctx, *_ = m6c.dataset(variant="headline")

    # the replacement posterior's derived samples, drawn exactly as in m6c_conjunction.py
    res = m6c_repeat.replacement_cells()
    p = PartitionedResult(tuple(res), err_scale=ERR_INFLATION)
    rng = np.random.default_rng(list(m6c_repeat.REPORT_RNG))
    eq = p.equal_weight(rng)
    idx = rng.choice(len(eq), size=min(m6c.DERIVED_MAX_SAMPLES, len(eq)), replace=False)
    th = np.array([res[0].space.theta(v) for v in eq[idx]])

    marks = {"conjunction": conj["epoch"]["jd_tdb"], "galle_night": conj["galle_night_jd_tdb"]}
    marks |= {k: float(Time(v, scale="tdb").jd) for k, v in WHAT_IF.items()}
    grid = np.linspace(JD_1800, JD_1870, N_GRID)
    jd = np.unique(np.concatenate([grid, list(marks.values())]))
    rot = np.array([precession_matrix_ecl_j2000_to_ecl_of_date(t) for t in jd])

    def lon_deg(pos: np.ndarray) -> np.ndarray:  # (n_t, 3) ecl_j2000 -> (n_t,) of date
        v = np.einsum("tij,tj->ti", rot, pos)
        return np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360.0

    # network first: a Horizons failure should cost seconds, not the integrations
    jpl = lon_deg(about("neptune", jd)[:, :3])

    cached = np.load(LON_CACHE) if LON_CACHE.exists() else None
    if cached is not None and np.array_equal(cached["jd"], jd) and np.array_equal(
            cached["theta"], th):  # fmt: skip
        lon = cached["lon"]
        print(f"loaded {LON_CACHE}", flush=True)
    else:
        lon = np.empty((len(th), len(jd)))
        for i, t in enumerate(th):
            pos, _ = predict_states_ecl_j2000(t, jd, ctx)
            lon[i] = lon_deg(pos[:, I_PERTURBER])
            if i % 1000 == 0:
                print(f"samples {i}/{len(th)}", flush=True)
        LON_CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez(LON_CACHE, jd=jd, theta=th, lon=lon)

    pct = np.array([circular_interval_deg(lon[:, k], (2.5, 16.0, 50.0, 84.0, 97.5))
                    for k in range(len(jd))])  # fmt: skip
    med = pct[:, 2]
    rel = wrap180(pct - med[:, None])
    lv = lon_deg(np.array([m6c.le_verrier_state_ecl_j2000(t)[0] for t in jd]))
    lv_off, jpl_off = wrap180(lv - med), wrap180(jpl - med)
    year = Time(jd, format="jd", scale="tdb").decimalyear

    # consistency with m6c_conjunction.json (same samples, same frame)
    for key, name in (("conjunction", "conjunction_helio"), ("galle_night", "galle_night_helio")):
        k = int(np.flatnonzero(jd == marks[key])[0])
        ref = conj[name]["all"]
        got = {"median": med[k], "lo68": pct[k, 1], "hi68": pct[k, 3]}
        assert all(abs(wrap180(got[s] - ref[s])) < 2e-3 for s in got), (key, got, ref)

    # years in which Le Verrier's offset from the median changes sign (linear interpolation)
    s = np.flatnonzero(np.sign(lv_off[:-1]) != np.sign(lv_off[1:]))
    crossings = [round(float(year[i] + (year[i + 1] - year[i]) * lv_off[i]
                             / (lv_off[i] - lv_off[i + 1])), 2) for i in s]  # fmt: skip

    def at(jd_mark: float) -> dict:
        k = int(np.flatnonzero(jd == jd_mark)[0])
        half = rel[k, 3] if lv_off[k] > 0 else -rel[k, 1]
        return {"jd_tdb": jd_mark, "date_tdb": Time(jd_mark, format="jd", scale="tdb").iso[:16],
                "year": round(float(year[k]), 3), "median_deg": round(float(med[k]), 3),
                "band_68_rel_deg": [round(float(rel[k, 1]), 3), round(float(rel[k, 3]), 3)],
                "band_95_rel_deg": [round(float(rel[k, 0]), 3), round(float(rel[k, 4]), 3)],
                "le_verrier_minus_median_deg": round(float(lv_off[k]), 3),
                "le_verrier_in_68_half_widths": round(float(lv_off[k] / half), 2),
                "jpl_minus_median_deg": round(float(jpl_off[k]), 3),
                "le_verrier_minus_jpl_deg": round(float(wrap180(lv[k] - jpl[k])), 3)}  # fmt: skip

    payload = {
        "description": __doc__.strip().splitlines()[0],
        "posterior": "M6c replacement: headline cells 0, 1, 3 + merged R2 cell 2 "
                     "(site/data/m6c_replacement.json)",
        "n_samples": len(th),
        "frame": "heliocentric (inner barycentre for the posterior and JPL; Sun for Le Verrier's "
                 "Kepler ellipse), ecliptic of date",
        "le_verrier_sources": list(LE_VERRIER_1846.sources),
        "le_verrier_zero_crossings_year": crossings,
        "marks": {k: at(v) for k, v in marks.items()},
        "series": {
            "year": np.round(year, 4).tolist(),
            "median_deg": np.round(med, 3).tolist(),
            "lo95_rel_deg": np.round(rel[:, 0], 3).tolist(),
            "lo68_rel_deg": np.round(rel[:, 1], 3).tolist(),
            "hi68_rel_deg": np.round(rel[:, 3], 3).tolist(),
            "hi95_rel_deg": np.round(rel[:, 4], 3).tolist(),
            "le_verrier_minus_median_deg": np.round(lv_off, 3).tolist(),
            "jpl_minus_median_deg": np.round(jpl_off, 3).tolist(),
        },
    }  # fmt: skip
    OUT_JSON.write_text(json.dumps(payload, indent=1) + "\n")
    print(json.dumps({k: payload[k] for k in ("le_verrier_zero_crossings_year", "marks")},
                     indent=1))  # fmt: skip
    print(f"wrote {OUT_JSON}")


def plot() -> None:
    """Render the figure from the saved JSON only: no posterior, no integrations."""
    d = json.loads(OUT_JSON.read_text())
    s = {k: np.array(v) for k, v in d["series"].items()}
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                         "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK_2,
                         "xtick.color": INK_2, "ytick.color": INK_2})  # fmt: skip
    fig, ax = plt.subplots(figsize=(10, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    x = s["year"]
    ax.fill_between(x, s["lo95_rel_deg"], s["hi95_rel_deg"], color=BAND, alpha=0.16,
                    linewidth=0, label="posterior 95 %")  # fmt: skip
    ax.fill_between(x, s["lo68_rel_deg"], s["hi68_rel_deg"], color=BAND, alpha=0.36,
                    linewidth=0, label="posterior 68 %")  # fmt: skip
    ax.axhline(0.0, color=INK_MUTED, linewidth=0.8)
    ax.plot(x, s["jpl_minus_median_deg"], color=INK, linewidth=1.5, label="JPL Neptune")
    ax.plot(x, s["le_verrier_minus_median_deg"], color=INK, linewidth=2, linestyle="--",
            label="Le Verrier's planet")  # fmt: skip

    ylo, yhi = ax.get_ylim()
    m = d["marks"]
    vlines = [(m["conjunction"]["year"], "conjunction"), (m["galle_night"]["year"], "Galle")]
    vlines += [(y, "LV crosses 0") for y in d["le_verrier_zero_crossings_year"]]
    for xv, lab in vlines:
        ax.axvline(xv, color=INK_MUTED, linewidth=0.8, linestyle=":")
        ax.annotate(lab, (xv, yhi), xytext=(3, -3), textcoords="offset points", color=INK_2,
                    fontsize=8, rotation=90, va="top")  # fmt: skip
    # selective labels: Le Verrier minus Neptune where the question is asked
    # (the first says what the number is; the others sit clear of the JPL line)
    for key, prefix, off, ha in (("search_1836", "LV − Neptune ", (6, 6), "left"),
                                 ("galle_night", "", (-6, -14), "right"),
                                 ("search_1856", "", (6, 6), "left")):  # fmt: skip
        k = m[key]
        ax.plot([k["year"]], [k["le_verrier_minus_median_deg"]], "o", color=INK, markersize=5)
        ax.annotate(f"{prefix}{k['le_verrier_minus_jpl_deg']:+.1f}°",
                    (k["year"], k["le_verrier_minus_median_deg"]), xytext=off, ha=ha,
                    textcoords="offset points", color=INK, fontsize=8)  # fmt: skip
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel("year")
    ax.set_ylabel("heliocentric longitude minus posterior median (°)")
    ax.set_title("Le Verrier's planet against the posterior, 1800–1870", color=INK,
                 fontsize=10, loc="left")  # fmt: skip
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="lower left", labelcolor=INK_2, ncol=4)
    fig.tight_layout()
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
