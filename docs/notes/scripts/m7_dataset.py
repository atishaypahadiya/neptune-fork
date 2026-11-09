"""The dataset panel's data (site/data/m7_dataset.json), from primary sources only.

- The published package, release/leverrier-table77/ (read-only; frozen): the table, the derived
  measurements, the misprint and resolution logs, both spot-check records, its CITATION.cff.
- Zenodo's public records API, for what the archive itself says: version DOI, concept DOI,
  version, publication date, licence and record URL.

Run from the repo root: `uv run python docs/notes/scripts/m7_dataset.py` (network).
"""

from __future__ import annotations

import csv
import hashlib
import json
import urllib.request
from collections import Counter, OrderedDict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
PKG = ROOT / "release" / "leverrier-table77"
OUT = ROOT / "site" / "data" / "m7_dataset.json"
RECORD = "https://zenodo.org/api/records/24044792"


def rows(name: str) -> list[dict]:
    with open(PKG / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def jd_to_year(jd: float) -> float:
    return 2000.0 + (jd - 2451545.0) / 365.25


def main() -> None:
    table = rows("leverrier_1846_table77.csv")
    meas = rows("leverrier_1846.csv")
    misprints = rows("leverrier_1846_misprints.csv")
    resolved = rows("leverrier_1846_resolution_log.csv")
    auto = json.loads((PKG / "spotcheck_results.json").read_text())
    human = json.loads((PKG / "spotcheck_human_results.json").read_text())
    cff = yaml.safe_load((PKG / "CITATION.cff").read_text())

    def key(r: dict) -> tuple[str, str]:
        return (str(r["page"]), str(int(r["line"])))

    spot = {key(r) for r in human["rows"]}
    mis = Counter(key(r) for r in misprints)
    res = Counter(key(r) for r in resolved)

    # One record per observation, from the derived measurements (RA and/or Dec per row).
    obs: OrderedDict[str, dict] = OrderedDict()
    for m in meas:
        o = obs.setdefault(
            m["obs_id"],
            {
                "id": m["obs_id"],
                "year": round(jd_to_year(float(m["epoch_tt_jd"])), 3),
                "observatory": m["observatory"],
                "observer": m["observer"],
                "quantities": [],
            },
        )
        o["quantities"].append(m["quantity"])
    for o in obs.values():
        page, line = o["id"].split("-p")[1].split("-l")
        k = (page, str(int(line)))
        o.update(spot_checked=k in spot, misprints=mis[k], resolved=res[k])

    def cells(record: dict) -> tuple[int, int]:
        n = ok = 0
        for r in record["rows"]:
            for c in r["cells"].values():
                n += 1
                ok += bool(c["matches"])
        return n, ok

    n_auto, ok_auto = cells(auto)
    n_human, ok_human = cells(human)

    with urllib.request.urlopen(RECORD, timeout=30) as r:
        z = json.load(r)
    md = z["metadata"]

    out = {
        "description": "Le Verrier (1846) Table 77 dataset, for the site's dataset panel. Built by "
        "docs/notes/scripts/m7_dataset.py from release/leverrier-table77/ and Zenodo's "
        "records API.",
        "title": md["title"],
        "creators": [c["name"] for c in md["creators"]],
        "zenodo": {
            "version_doi": z["doi"],
            "concept_doi": z["conceptdoi"],
            "version": md.get("version"),
            "publication_date": md["publication_date"],
            "license": md["license"]["id"],
            "record_url": z["links"]["self_html"],
            "concept_url": f"https://doi.org/{z['conceptdoi']}",
        },
        "citation": {
            "message": " ".join(cff["message"].split()),
            "source": cff["references"][0],
        },
        "counts": {
            "rows": len(table),
            "pages": len({r["page"] for r in table}),
            "measurements": len(meas),
            "observations_with_measurements": len(obs),
            "years": [min(o["year"] for o in obs.values()), max(o["year"] for o in obs.values())],
            "pass_disagreements": len(resolved),
            "misprints": len(misprints),
        },
        "spot_checks": {
            "automated": {
                "rows": len(auto["rows"]),
                "cells": n_auto,
                "matching": ok_auto,
                "file": "spotcheck_results.json",
            },
            "human": {
                "rows": len(human["rows"]),
                "cells": n_human,
                "matching": ok_human,
                "reviewer": "the dataset author",
                "file": "spotcheck_human_results.json",
            },
            "same_rows": {key(r) for r in auto["rows"]} == spot,
        },
        "provenance": "Machine-assisted transcription and reconciliation (an AI model), not done "
        "by an independent human; 28 rows then spot-checked by a human against the Gallica scan.",
        "package_files": sorted(
            p.name for p in PKG.iterdir() if p.is_file() and not p.name.startswith(".")
        ),
        "table_sha256": hashlib.sha256(
            (PKG / "leverrier_1846_table77.csv").read_bytes()
        ).hexdigest(),
        "observations": list(obs.values()),
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print({k: out[k] for k in ("zenodo", "counts", "spot_checks")})
    print(Counter(o["observatory"] for o in obs.values()).most_common())


if __name__ == "__main__":
    main()
