"""Compute ledger totals → ``site/data/ledger.json`` (and a one-time backfill).

    uv run python -m fork.export.ledger              # totals from data/cache/compute_ledger.jsonl
    uv run python -m fork.export.ledger --backfill   # add estimated rows for past runs, then total

Backfill (rows marked ``"estimated": true``; jobs already in the ledger are skipped):
  - dynesty results in data/cache/nested/*.pkl (the M5b validation runs): likelihood calls and
    wall time are exact (from the saved results); cpu_s is estimated as wall × workers, an
    upper bound; workers from the validation logs.
  - emcee chains in data/cache/chains/*.npz: likelihood calls = steps × walkers where the
    chain recorded its walkers; wall and workers were not recorded.
  - ``LOGGED_RUNS``: runs whose wall time is in this session's queue logs but that logged no
    evaluation counts (M5a diagnostics, point-fit regeneration, profile export, slow tests).
One likelihood call = one forward integration over the observation span (``fork.ledger``).
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

from fork import ledger
from fork.export.to_json import SITE_DATA, write_json

ROOT = ledger.ROOT
OUT_NAME = "ledger.json"
SPAN_YEARS = 155.0  # sim_years per integration on the synthetic epochs (1690.9–1845.9)
NESTED_DIR = ROOT / "data" / "cache" / "nested"
CHAIN_DIR = ROOT / "data" / "cache" / "chains"
# workers per validation run, from data/cache/nested/progress.log and M5b_dynesty_validation.md
VALIDATION_WORKERS = {"seed_0_free_a_e0.3_rng1": 6, "seed_0_free_a_e0.3_rng2": 6,
                      "seed_0_free_a_e0.3_slice0_rng3": 3}  # fmt: skip
VALIDATION_WORKERS_DEFAULT = {"slice": 4, "cell": 3}


def _t(s: str) -> datetime:
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S")


def _span_s(a: str, b: str) -> float:
    return (_t(b) - _t(a)).total_seconds()


# (job, kind, wall_s, n_workers, note) — wall times from this session's queue logs, 2026-11-04
LOGGED_RUNS: tuple[tuple[str, str, float, int, str], ...] = (
    ("m5a/diagnostics_rerun", "point_fits", _span_s("2026-11-04 06:07:00", "2026-11-04 07:01:34"),
     12, "54 multistart fits (6 datasets × 9), 24 starts; nfev not logged"),
    ("m5a/diagnostics_48_start/seed_2_free_a", "point_fits", 234.0, 12, "48-start rerun"),
    ("m5a/diagnostics_48_start/seed_3_le_verrier", "point_fits", 67.0, 12, "48-start rerun"),
    ("m5a/fits_regen", "point_fits", _span_s("2026-11-04 07:07:03", "2026-11-04 07:08:23"), 12,
     "site/data/m5a_fits.json, seed 0, two modes"),
    ("m5a/profile_reexport", "profile", _span_s("2026-11-04 07:14:49", "2026-11-04 10:41:22"),
     12, "site/data/m5a_profile.json, 16 curves × 51 grid points, 8 starts + continuation; "
     "LM evaluations not logged"),
    ("m5a/profile_reexport_post", "profile", 11629.0, 12,
     "site/data/m5a_profile.json, 11 free-a curves × 52 grid points (Le Verrier mode cached), "
     "8 starts + 2 a = 22 AU backstop + 4 M5b-posterior-seeded per point, + continuation; "
     "2026-11-05 17:36-20:50 (wall from the export log); LM evaluations not logged"),
    ("tests/slow_point_fit", "tests", 384.8, 12,
     "test_perturber_fit, test_export, test_profile, test_a3a (slow)"),
)  # fmt: skip


def backfill_rows() -> list[dict]:
    import numpy as np

    from fork.inference.nested import NestedResult

    rows = []
    for p in sorted(NESTED_DIR.glob("*.pkl")):
        r = NestedResult.load(p)
        stem = p.stem
        kind = "cell" if "_cell" in stem else "slice" if "_slice" in stem else None
        nw = VALIDATION_WORKERS.get(stem, VALIDATION_WORKERS_DEFAULT.get(kind or "", None))
        rows.append(
            dict(
                job=f"m5b_validation/{stem}",
                kind="nested",
                n_likelihood_calls=r.ncall,
                n_integrations=r.ncall,
                sim_years=r.ncall * SPAN_YEARS,
                wall_s=r.wall_s,
                cpu_s=None if nw is None else r.wall_s * nw,
                n_workers=nw,
                estimated=True,
                nlive=r.nlive,
                niter=r.niter,
                note="calls and wall exact; cpu_s = wall × workers (upper bound)",
            )
        )
    for p in sorted(CHAIN_DIR.glob("*.npz")):
        meta = json.loads(str(np.load(p)["meta"]))
        steps, walkers = meta.get("n_steps"), meta.get("n_walkers")
        calls = None if steps is None or walkers is None else int(steps) * int(walkers)
        rows.append(dict(job=f"emcee/{p.stem}", kind="emcee", n_likelihood_calls=calls,
                         n_integrations=calls,
                         sim_years=None if calls is None else calls * SPAN_YEARS,
                         wall_s=meta.get("wall_s"), cpu_s=None, n_workers=None, estimated=True,
                         note="calls = steps × walkers" if calls is not None else
                         "walkers not recorded; calls unknown"))  # fmt: skip
    for job, kind, wall, nw, note in LOGGED_RUNS:
        rows.append(dict(job=job, kind=kind, n_likelihood_calls=None, n_integrations=None,
                         sim_years=None, wall_s=wall, cpu_s=None, n_workers=nw,
                         estimated=True, note=note + "; wall from queue logs"))  # fmt: skip
    return rows


def backfill() -> int:
    have = {r["job"] for r in ledger.read()}
    n = 0
    for row in backfill_rows():
        if row["job"] in have:
            continue
        job, kind = row.pop("job"), row.pop("kind")
        ledger.append(job, kind, **row)
        n += 1
    return n


def _sum(rows: list[dict], key: str) -> float | int | None:
    vals = [r[key] for r in rows if r.get(key) is not None]
    if not vals:
        return None
    return round(sum(vals), 1) if any(isinstance(v, float) for v in vals) else sum(vals)


def totals(rows: list[dict]) -> dict:
    keys = ("n_likelihood_calls", "n_integrations", "sim_years", "wall_s", "cpu_s")

    def block(rs: list[dict]) -> dict:
        out = {"n_rows": len(rs), **{k: _sum(rs, k) for k in keys}}
        out["n_rows_missing"] = {k: sum(r.get(k) is None for r in rs) for k in keys}
        if out["cpu_s"] is not None:
            out["cpu_hours"] = round(out["cpu_s"] / 3600.0, 2)
        return out

    kinds = sorted({r["kind"] for r in rows})
    return {
        "all": block(rows),
        "measured": block([r for r in rows if not r.get("estimated")]),
        "estimated": block([r for r in rows if r.get("estimated")]),
        "by_kind": {k: block([r for r in rows if r["kind"] == k]) for k in kinds},
    }


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backfill", action="store_true")
    args = ap.parse_args(argv)
    if args.backfill:
        print(f"backfilled {backfill()} rows")
    rows = ledger.read()
    payload = {
        "description": "Compute ledger totals (data/cache/compute_ledger.jsonl). One likelihood "
        "call = one forward N-body integration over the observation span. Sums skip missing "
        "values; n_rows_missing counts them. 'estimated' rows were backfilled from saved "
        "results and logs.",
        "fields": ["job", "kind", "n_likelihood_calls", "n_integrations", "sim_years", "wall_s",
                   "cpu_s", "n_workers", "estimated"],
        "totals": totals(rows),
        "rows": rows,
    }  # fmt: skip
    write_json(OUT_NAME, payload)
    print(json.dumps(payload["totals"]["all"], indent=1))
    print(f"wrote {Path(SITE_DATA) / OUT_NAME}")


if __name__ == "__main__":
    main()
