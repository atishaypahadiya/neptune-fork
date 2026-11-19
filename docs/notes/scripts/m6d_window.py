"""Post hoc sensitivity of the M6d window (docs/notes/M6d.md, note of 19 November 2026).

The registered window, 0.87°, is Le Verrier's own 52′: the difference in heliocentric longitude at
1847 January 1 between Galle's planet and the longitude he had published (Recherches, 1846,
p. 253). It is not his geocentric miss on Galle's night, which from his published elements is
0.90° (site/data/m6_fits.json, comparison_targets). This recomputes P_free, P_shape and the
one-sided 95 % bound at both half-widths, from the cached draws and by the code that made
docs/notes/M6d.json, after reproducing the registered values at 0.87°. No fit is rerun.

    uv run python docs/notes/scripts/m6d_window.py      # needs data/cache/m8 and nested/v14_fits
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from fork.export import v13
from fork.export import v14_fits as V

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "docs" / "notes" / "M6d_window.json"
REGISTERED = 0.87
SENSITIVITY = 0.90
HELD = ("shape38", "shape36")


def _cp_upper(n: int, conf: float = 0.95) -> float:
    """Exact one-sided Clopper–Pearson upper bound for 0 successes in n trials."""
    return 1.0 - (1.0 - conf) ** (1.0 / n)


def _row(w: float, free: np.ndarray, held: dict[str, dict[str, np.ndarray]], target: float) -> dict:
    row: dict = {"half_width_deg": w, "p_free": v13.peirce_probability(free, target, w)}
    for name, g in held.items():
        d = np.abs((g["merged"] - target + 180.0) % 360.0 - 180.0)
        row[name] = {
            "n_draws": int(g["merged"].size),
            "p_shape": v13.peirce_probability(g["merged"], target, w),
            "p_shape_single_seed": [
                v13.peirce_probability(x, target, w) for k, x in g.items() if k != "merged"
            ],
            "nearest_draw_deg": round(float(d.min()), 3),
        }
    n = row[HELD[0]]["n_draws"]
    row["cp_upper_95"] = _cp_upper(n)
    # As on the site: P_free over the bound, rounded down to the hundred.
    row["factor_at_least"] = math.floor(row["p_free"] / _cp_upper(n) / 100) * 100
    return row


def main() -> None:
    reg = json.loads((ROOT / "docs/notes/M6d.json").read_text(encoding="utf-8"))
    tgt = json.loads((ROOT / "site/data/m6_fits.json").read_text(encoding="utf-8"))
    tgt = tgt["comparison_targets"]
    free, _ = V.free_distance_galle()  # raises unless the published median and ln Z reproduce
    held = {name: V.assemble(V.FITS[name])["_galle"] for name in HELD}
    target = v13.GALLE_LON_DEG
    if v13.PEIRCE_HALF_WIDTH_DEG != REGISTERED:
        raise RuntimeError("the registered half-width is not 0.87°")
    rows = {"registered": _row(REGISTERED, free, held, target),
            "sensitivity": _row(SENSITIVITY, free, held, target)}  # fmt: skip
    for name in HELD:
        r = reg["fits"][name]["summary"]["reading"]
        if (
            abs(rows["registered"]["p_free"] - r["free_distance"]["p_within_0p87_deg"]) > 1e-12
            or rows["registered"][name]["p_shape"] != r["ratio"]["p_shape"]
        ):
            raise RuntimeError(f"0.87° does not reproduce docs/notes/M6d.json for {name}")
    geo = tgt["jpl"]["galle_lon_ecl_of_date_deg"]
    geo -= tgt["le_verrier"]["galle_lon_ecl_of_date_deg_derived"]
    out = {
        "description": (
            "Post hoc sensitivity of the M6d window: the registered 0.87° (Le Verrier's 52′, a "
            "heliocentric difference at 1847 January 1, Recherches p. 253) and 0.90° (his "
            "geocentric miss on Galle's night from his published elements). Built by "
            "docs/notes/scripts/m6d_window.py from the cached draws; no fit rerun."
        ),
        "post_hoc": True,
        "date": "2026-11-19",
        "target_deg": target,
        "window_52_arcmin": {
            "deg": REGISTERED,
            "measures": "heliocentric longitude at 1847 January 1, Galle's planet against "
            "the longitude Le Verrier published on 31 August 1846",
            "source": "Le Verrier, Recherches sur les mouvements de la planète Herschel (1846), "
            "p. 253: 327° 24′ against 326° 32′, difference 0° 52′",
        },
        "geocentric_miss_deg": round(geo, 4),
        "geocentric_miss_source": "site/data/m6_fits.json comparison_targets: jpl minus "
        "le_verrier.galle_lon_ecl_of_date_deg_derived",
        "n_free": int(free.size),
        "free_median_deg": round(float(np.median(free)), 3),
        "reproduces_registered": True,
        "windows": rows,
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
