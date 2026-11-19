"""Deviation 1's results on the study's cache against those on its scratch copy: the
consistency check of the conditional go for Stage 1's gates, (c).

Every row of the scratch copy is compared in every box and level of ln Z, and in the V4 minimum
and the reported best value. A row may differ only if it is

- one of the three rows the read-back fix (ffbb963) touches; or
- a row of truth c's or truth d's planet group (``t3_…_P_…`` or ``t4_…_P_…``, halves excluded),
  which gained a task after the scratch copy was made, and then only where the study's cache
  holds a deviation 1 fit at a strength the scratch copy's row lacks, and every ln Z lies within
  the row's tolerance on the scratch copy (TAU_MAX where that row had none).

Anything else, a row beyond its tolerance, or a row missing from the study's cache fails.
Standard library only, so it runs beside any commit of the study's code."""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

READBACK = {("t2_s0_l0_J", 0.0), ("t4_s0_l0_R", 0.12), ("t2_s0_l0_P_d1_b0", 5.6287)}
GROUPS = re.compile(r"^t[34]_s0_l0_P_d\d_b\d+$")
SAME = 1e-9
TAU_MAX = 0.1
FIELDS = ("chi2_min_calls", "x_best", "x_best_argmin")


def _same(a, b) -> bool:
    if a == b:
        return True
    if a is None or b is None:
        return False
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= SAME


def _fits(root: Path, tid: str, p2: float) -> set[float]:
    f = root / "scans" / tid / "refine" / f"{p2:.10f}.jsonl"
    if not f.exists():
        return set()
    return {round(json.loads(x)["x"], 9) for x in f.read_text().splitlines() if x.strip()}


def _row_diffs(s: dict, r: dict) -> list[dict]:
    out = []
    for lv, boxes in s["ln_z"].items():
        for b, v in boxes.items():
            w = r["ln_z"][lv][b]
            if not _same(v, w):
                d = w - v if math.isfinite(v) and math.isfinite(w) else math.inf
                out.append({"box": b, "level": lv, "scratch": v, "real": w, "diff": d})
    for f in FIELDS:
        if not _same(s.get(f), r.get(f)):
            out.append({"field": f, "scratch": s.get(f), "real": r.get(f)})
    if s.get("x_best_bound") != r.get("x_best_bound"):
        out.append({"field": "x_best_bound", "scratch": s.get("x_best_bound"),
                    "real": r.get("x_best_bound")})  # fmt: skip
    return out


def compare(real: Path, scratch: Path) -> dict:
    out: dict = {"compared_rows": 0, "readback_rows": [], "group_rows": [], "differences": [],
                 "missing": []}  # fmt: skip
    for d in sorted((Path(scratch) / "scans").iterdir()):
        sp = d / "result.json"
        if not sp.exists() or d.name.endswith("_NG"):
            continue
        rp = Path(real) / "scans" / d.name / "result.json"
        r_rows = {}
        if rp.exists():
            r_rows = {round(x["p2"], 9): x for x in json.loads(rp.read_text())["scans"]}
        for s in json.loads(sp.read_text())["scans"]:
            key = (d.name, round(s["p2"], 4))
            r = r_rows.get(round(s["p2"], 9))
            if r is None or "ln_z_before_d1" not in r:
                out["missing"].append(list(key))
                continue
            out["compared_rows"] += 1
            diffs = _row_diffs(s, r)
            if not diffs:
                continue
            item = {"task": d.name, "p2": s["p2"], "diffs": diffs}
            if key in READBACK:
                out["readback_rows"].append(item)
                continue
            if not GROUPS.match(d.name):
                out["differences"].append(item | {"reason": "outside the allowed rows"})
                continue
            new = sorted(_fits(Path(real), d.name, s["p2"]) - _fits(Path(scratch), d.name, s["p2"]))
            tol = (s.get("d1") or {}).get("tolerance") or {}
            tau = {b: tol.get(b) or TAU_MAX for b in s["ln_z"][next(iter(s["ln_z"]))]}
            over = max((abs(x["diff"]) / tau[x["box"]] for x in diffs if "box" in x), default=0.0)
            item |= {"new_fits": new, "tau": tau, "largest_over_tau": over}
            if not new:
                out["differences"].append(item | {"reason": "no fit added in this run"})
            elif over > 1.0:
                out["differences"].append(item | {"reason": "beyond the scratch tolerance"})
            else:
                out["group_rows"].append(item)
    out["pass"] = not out["differences"] and not out["missing"]
    return out


def main(argv: list[str]) -> int:
    """``d1_consistency.py <study cache> <scratch copy> <out.json>``"""
    real, scratch, dest = (Path(a) for a in argv)
    r = compare(real, scratch)
    dest.write_text(json.dumps(r, indent=1, default=float))
    print(f"compared {r['compared_rows']} rows; read-back rows differing "
          f"{len(r['readback_rows'])}; planet-group rows within tolerance "
          f"{len(r['group_rows'])}; other differences {len(r['differences'])}; missing "
          f"{len(r['missing'])}; pass {r['pass']}")  # fmt: skip
    for x in r["differences"][:20]:
        print("  ", x["task"], round(x["p2"], 4), x["reason"], x["diffs"][:2])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
