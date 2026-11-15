"""The main page's data on the S1 cell-0 rerun (site/data/s1_cell0.json, a subset of
docs/notes/S1_cell0.json) and on the low-distance plan's results (site/data/s1_mode.json, a subset
of docs/notes/S1_mode.json: the 95 % range it sets and its status; never P(a < 25 AU), which
that plan's § 6 does not let it state while Part A is not resolved).

Every entry of ``values`` is copied from one field of the committed note JSON, named in
``sources`` as ``[file, key, key, ...]``; nothing is computed here (docs/DESIGN.md, the one
rule). The page uses it for one clause beside the 95 % range while the rerun is unresolved.
``tests/test_s1_site.py`` resolves every source.

    uv run python -m fork.export.s1_site
"""

from __future__ import annotations

import json
from pathlib import Path

from fork.export.m6e_site import resolve

ROOT = Path(__file__).resolve().parents[3]
NOTE = "docs/notes/S1_cell0.json"
MODE_NOTE = "docs/notes/S1_mode.json"
SITE_DATA = ROOT / "site" / "data"


def sources() -> dict[str, list]:
    r = [NOTE, "reading"]
    return {
        "resolved": [*r, "resolved"],
        "outcome": [*r, "outcome"],
        "verdict": [*r, "verdict"],
        "label_95": [*r, "label_95"],
        # The outer edge of the rerun's cell: the last sub-cell's upper distance.
        "a_max_au": [NOTE, "stage_a", "subcells", -1, "a_au", 1],
    }


def build() -> dict:
    src, cache = sources(), {}
    return {
        "description": (
            "The rerun of the eccentricity test's nearest cell in six sub-cells (main page, "
            "beside the 95 % range): copies of fields of docs/notes/S1_cell0.json."
        ),
        "plan": "docs/notes/S1_cell0.md",
        "sources": src,
        "values": {k: resolve(v, cache) for k, v in src.items()},
    }


def mode_sources() -> dict[str, list]:
    r = [MODE_NOTE, "range_95"]
    return {
        "lo95": [*r, "range", 0],
        "hi95": [*r, "range", 1],
        "err_lo95": [*r, "sampler_errors", 0],
        "err_hi95": [*r, "sampler_errors", 1],
        "final": [*r, "final"],
        "label": [*r, "label"],
        "why": [*r, "why"],
        "rule_b_holds": [MODE_NOTE, "part_b", "holds"],
        "part_a_resolved": [MODE_NOTE, "part_a", "reading", "resolved"],
        "part_a_status": [MODE_NOTE, "part_a", "reading", "status"],
    }


def build_mode() -> dict:
    src, cache = mode_sources(), {}
    return {
        "description": (
            "The low-distance plan's 95 % range and its status (main page, the 95 % range): "
            "copies of fields of docs/notes/S1_mode.json."
        ),
        "plan": "docs/notes/S1_mode.md",
        "sources": src,
        "values": {k: resolve(v, cache) for k, v in src.items()},
    }


def main() -> None:
    for name, data in (("s1_cell0.json", build()), ("s1_mode.json", build_mode())):
        path = SITE_DATA / name
        path.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
