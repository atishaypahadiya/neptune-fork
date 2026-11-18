"""The main page's data on Le Verrier's shape model against free distance (site/data/m6d.json,
a subset of docs/notes/M6d.json: the two shape fits' place on Galle's night, the Bayes factor
against free distance and the probability of landing as close as Le Verrier did, with their
registered wordings).

Every entry of ``values`` is copied from one field of the committed note JSON, named in
``sources`` as ``[file, key, key, ...]``; nothing is computed here (docs/DESIGN.md, the one
rule). ``tests/test_m6d_site.py`` resolves every source.

    uv run python -m fork.export.m6d_site
"""

from __future__ import annotations

import json
from pathlib import Path

from fork.export.m6e_site import resolve

ROOT = Path(__file__).resolve().parents[3]
NOTE = "docs/notes/M6d.json"
SITE_DATA = ROOT / "site" / "data"
FITS = ("shape38", "shape36")
WINDOW_DEG = 0.87
WINDOW_ARCMIN = 52
WINDOW_CITATION = "docs/notes/M6d.md, § 1: Le Verrier's actual miss of 52′"
# Where each held value comes from, in the plan's words (§ 3, "Held value, two runs").
HELD = {
    "shape38": "Bode's law, his first paper",
    "shape36": "from his elements of 31 August 1846",
    "citation": "docs/notes/M6d.md, § 3, Held value, two runs",
}


def sources() -> dict[str, list]:
    out: dict[str, list] = {}
    for f in FITS:
        s = [NOTE, "fits", f, "summary"]
        g = [*s, "fits", "merged", "galle_lon_ecl_of_date_deg"]
        r = [*s, "reading"]
        out |= {
            f"{f}.held_a_au": [*s, "held_a_au"],
            f"{f}.read": [*s, "read"],
            f"{f}.n_draws": [*s, "settings", "n_draws"],
            f"{f}.median": [*g, "median"],
            f"{f}.lo68": [*g, "lo68"],
            f"{f}.hi68": [*g, "hi68"],
            f"{f}.median_minus_neptune_deg": [*g, "vs", "jpl", "median_minus_target_deg"],
            f"{f}.ln_b": [*r, "bayes_factor", "ln_b"],
            f"{f}.ln_b_err": [*r, "bayes_factor", "err"],
            f"{f}.ln_b_wording": [*r, "bayes_factor", "wording"],
            f"{f}.p_shape": [*r, "ratio", "p_shape"],
            f"{f}.p_free": [*r, "ratio", "p_free"],
            f"{f}.ratio_wording": [*r, "ratio", "wording"],
        }
    return out


def build() -> dict:
    src, cache = sources(), {}
    return {
        "description": (
            "Le Verrier's shape model, the distance held at his values, against free distance "
            "(main page, § 4.5): copies of fields of docs/notes/M6d.json."
        ),
        "plan": "docs/notes/M6d.md",
        # The window is a constant of the plan, not a computed field: Le Verrier's actual miss.
        "window": {"deg": WINDOW_DEG, "arcmin": WINDOW_ARCMIN, "citation": WINDOW_CITATION},
        "held": HELD,
        "sources": src,
        "values": {k: resolve(v, cache) for k, v in src.items()},
    }


def main() -> None:
    path = SITE_DATA / "m6d.json"
    path.write_text(json.dumps(build(), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
