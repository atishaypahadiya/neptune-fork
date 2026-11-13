"""The main page's data for § 4.4 (site/data/m6e.json): a subset of docs/notes/M6e.json.

Every entry of ``values`` is copied from one field of the committed note JSON, named in
``sources`` as ``[file, key, key, ...]``; nothing is computed here, so the page's numbers
trace to the note (docs/DESIGN.md, the one rule). ``tests/test_m6e_site.py`` resolves every
source. This module reads the note's JSON only, never the held-out table.

    uv run python -m fork.export.m6e_site
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NOTE = "docs/notes/M6e.json"
SITE_DATA = ROOT / "site" / "data"
NIGHTS = ("1795-05-08", "1795-05-10")
MODELS = ("headline", "student_t")


def resolve(source: list, cache: dict | None = None):
    """The value a source names: ``[file, key, ...]`` followed through dicts and lists."""
    cache = {} if cache is None else cache
    name = source[0]
    if name not in cache:
        cache[name] = json.loads((ROOT / name).read_text(encoding="utf-8"))
    x = cache[name]
    for k in source[1:]:
        x = x[k]
    return x


def sources() -> dict[str, list]:
    s: dict[str, list] = {"n_draws": [NOTE, "n_draws"]}
    for night in NIGHTS:
        s[f"instant_tdb.{night}"] = [NOTE, "instants_tdb", night]
        obs = [NOTE, "models", "headline", "posterior", night, "observed"]
        s[f"lalande.{night}.ecl_lon_deg"] = [*obs, "ecl_lon_deg"]
        s[f"lalande.{night}.dec_deg"] = [*obs, "dec_deg"]
    for model in MODELS:
        m = [NOTE, "models", model]
        s[f"{model}.reading"] = [*m, "reading"]
        for night in NIGHTS:
            p, k = [*m, "posterior", night], f"{model}.{night}"
            s[f"{k}.level"] = [*p, "central_level_ecl_lon"]
            s[f"{k}.class"] = [*p, "class"]
            s[f"{k}.median_minus_lalande_deg"] = [*p, "median_minus_observed_deg", "ecl_lon"]
            for q in ("median", "lo68", "hi68", "lo95", "hi95", "width68"):
                s[f"{k}.lon.{q}"] = [*p, "ecl_lon_deg", q]
            for q in ("lo68", "hi68", "lo95", "hi95"):
                s[f"{k}.dec.{q}"] = [*p, "dec_deg", q]
            s[f"{k}.prior_width68"] = [*m, "prior_predictive", night, "ecl_lon_deg", "width68"]
            s[f"{k}.width_ratio68"] = [*m, "width_ratio_posterior_to_prior", night, "width68"]
    s["histogram.edges_deg"] = [NOTE, "models", "headline", "histogram_8_may", "edges_deg"]
    s["histogram.counts"] = [NOTE, "models", "headline", "histogram_8_may", "counts"]
    return s


def build() -> dict:
    src, cache = sources(), {}
    return {
        "description": (
            "The check against Lalande's places of May 1795 (main page, § 4.4): copies of "
            "fields of docs/notes/M6e.json, the pre-registered test's results. The histogram "
            "is of the headline posterior on the first night."
        ),
        "plan": "docs/notes/M6e.md",
        "nights": list(NIGHTS),
        "models": list(MODELS),
        "histogram_night": NIGHTS[0],
        "sources": src,
        "values": {k: resolve(v, cache) for k, v in src.items()},
    }


def main() -> None:
    path = SITE_DATA / "m6e.json"
    path.write_text(json.dumps(build(), indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
