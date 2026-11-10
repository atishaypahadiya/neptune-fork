"""Table 1 of the site (site/data/m7_table1.json): Table 77's observations of Uranus by source.

One row per source:
- the ancient observers named in Le Verrier's text: Flamsteed, Bradley, Mayer, Le Monnier;
- the modern series in the error model's groups: Greenwich 1781–1830, Paris 1801–1827 and
  Paris 1835–1845.

Per row: observatory, years, observations (Table 77 rows) and measurements (RA, Dec), and the
median absolute on-sky residual against JPL's Uranus, |O − JPL|. Per observation (for Fig. 1):
its source and signed residuals.

The residuals are the M6b calibration's own (``fork.data.calibration``: Table 77's values against
JPL's Uranus taken through the apparent-place model; RA as Δα·cos δ). Nothing in src/ is
changed. As a consistency check, the script recomputes the per-group median |O − JPL| from
these residuals and requires it to match the committed error model's table
(data/historical/m6_error_model.json) to 1e-6″.

Run from the repo root: `uv run python docs/notes/scripts/m7_table1.py` (cache only).
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from fork.data.calibration import real_observations
from fork.model import observation as om

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "site" / "data" / "m7_table1.json"
ERROR_MODEL = ROOT / "data" / "historical" / "m6_error_model.json"
SERIES_LABEL = {
    "greenwich": "Greenwich series",
    "paris_1801_1827": "Paris series",
    "paris_1835_1845": "Paris series",
}


def offline(*_args: object, **_kwargs: object) -> np.ndarray:
    raise RuntimeError("Table 1 reads the Horizons cache only; a fetch was requested")


def source_of(row: dict[str, str]) -> str:
    return row["observer"] if row["era"] == "ancient" else om.series_of(row)


def main() -> None:
    import fork.data.horizons as hz

    real = hz.state_vectors_ecl_j2000
    hz.state_vectors_ecl_j2000 = lambda body, jd, **kw: real(body, jd, **{**kw, "fetcher": offline})

    rows = om.read_rows()
    obs = real_observations(rows)
    resid = obs.residuals_arcsec(*om.true_uranus_ra_dec(obs))
    if len(resid) != len(rows):
        raise RuntimeError("residuals and rows are not aligned")

    # Consistency with the committed calibration, per error-model group.
    committed = {g["group"]: g["median_abs"] for g in json.loads(ERROR_MODEL.read_text())["table"]}
    by_group: dict[str, list[float]] = defaultdict(list)
    for r, x in zip(rows, resid, strict=True):
        by_group[om.group_of(r)].append(abs(x))
    for g, xs in by_group.items():
        if abs(float(np.median(xs)) - committed[g]) > 1e-6:
            raise RuntimeError(f"{g}: median |O-JPL| {np.median(xs)} != committed {committed[g]}")

    src: dict[str, dict] = {}
    for r, x in zip(rows, resid, strict=True):
        k = source_of(r)
        s = src.setdefault(
            k,
            {
                "years": [],
                "obs": set(),
                "ra": 0,
                "dec": 0,
                "abs": [],
                "observatory": r["observatory"],
            },
        )
        s["years"].append(int(r["epoch"][:4]))
        s["obs"].add(r["obs_id"])
        s[r["quantity"]] += 1
        s["abs"].append(abs(float(x)))
    out_rows = []
    for k, s in sorted(src.items(), key=lambda kv: min(kv[1]["years"])):
        out_rows.append(
            {
                "key": k,
                "source": SERIES_LABEL.get(k, k),
                "era": "ancient" if k not in SERIES_LABEL else "modern",
                "observatory": s["observatory"],
                "years": [min(s["years"]), max(s["years"])],
                "observations": len(s["obs"]),
                "measurements": {"ra": s["ra"], "dec": s["dec"]},
                "median_abs_o_minus_jpl_arcsec": round(float(np.median(s["abs"])), 2),
            }
        )
    per_obs: dict[str, dict] = {}
    for r, x in zip(rows, resid, strict=True):
        o = per_obs.setdefault(r["obs_id"], {"id": r["obs_id"], "source": source_of(r)})
        o[r["quantity"]] = round(float(x), 2)
    all_years = [y for s in src.values() for y in s["years"]]
    payload = {
        "description": "Table 1: Table 77's observations of Uranus by source. Built by "
        "docs/notes/scripts/m7_table1.py from data/historical/leverrier_1846.csv and the M6b "
        "calibration residuals (fork.data.calibration).",
        "residual": "median absolute on-sky residual against JPL's Uranus (apparent place of "
        "date; RA as Δα·cos δ), all of a source's measurements, arcsec",
        "grouping_note": "the error model pools the four ancient observers into one group "
        "(ancient RA, ancient Dec); the modern rows are its series",
        "rows": out_rows,
        "total": {
            "observations": len({r["obs_id"] for r in rows}),
            "measurements": {
                "ra": sum(s["ra"] for s in src.values()),
                "dec": sum(s["dec"] for s in src.values()),
            },
            "years": [min(all_years), max(all_years)],
            "median_abs_o_minus_jpl_arcsec": round(float(np.median(np.abs(resid))), 2),
        },
        "check": "per-group medians reproduce data/historical/m6_error_model.json to 1e-6 arcsec",
        "observations": list(per_obs.values()),
        "observations_note": "per observation (Table 77 row id): source key and the signed "
        "residual O - JPL of each measured coordinate, arcsec (RA as Δα·cos δ)",
    }
    OUT.write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")
    for r in out_rows:
        print(r)
    print(payload["total"])


if __name__ == "__main__":
    main()
