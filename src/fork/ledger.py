"""Compute accounting: one JSON line per job in ``data/cache/compute_ledger.jsonl``.

Fields: job, kind, n_likelihood_calls, n_integrations, sim_years, wall_s, cpu_s, n_workers
(plus ``estimated`` and optional notes). One likelihood call is one forward-model evaluation,
i.e. one N-body integration from the 1781 epoch out to every observation epoch, so
``sim_years`` = n_integrations × the span the observations cover (1690–1846 ≈ 156 yr).
``cpu_s`` is this process plus its reaped children (pool workers). Totals:
``python -m fork.export.ledger`` → ``site/data/ledger.json``.
"""

from __future__ import annotations

import json
import os
import resource
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "data" / "cache" / "compute_ledger.jsonl"
DAYS_PER_YEAR = 365.25


def sim_years_per_integration(jd_tdb: np.ndarray, epoch_jd_tdb: float) -> float:
    """Years one forward evaluation integrates: from the epoch back to the first observation
    and forward to the last."""
    t = np.asarray(jd_tdb, dtype=np.float64)
    return float((max(t.max(), epoch_jd_tdb) - min(t.min(), epoch_jd_tdb)) / DAYS_PER_YEAR)


def cpu_now() -> float:
    """CPU seconds (user + system) of this process and its terminated, reaped children."""
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return s.ru_utime + s.ru_stime + c.ru_utime + c.ru_stime


def append(job: str, kind: str, *, n_likelihood_calls: int | None, n_integrations: int | None,
           sim_years: float | None, wall_s: float | None, cpu_s: float | None,
           n_workers: int | None, estimated: bool = False, path: Path | None = None,
           **extra) -> dict:  # fmt: skip
    """Append one row (a single O_APPEND write, so concurrent jobs do not interleave)."""
    row = {"job": job, "kind": kind, "n_likelihood_calls": n_likelihood_calls,
           "n_integrations": n_integrations,
           "sim_years": None if sim_years is None else round(sim_years, 1),
           "wall_s": None if wall_s is None else round(wall_s, 1),
           "cpu_s": None if cpu_s is None else round(cpu_s, 1), "n_workers": n_workers,
           "estimated": bool(estimated), "recorded": time.strftime("%Y-%m-%dT%H:%M:%S"),
           **extra}  # fmt: skip
    p = LEDGER if path is None else Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    line = (json.dumps(row, ensure_ascii=False) + "\n").encode()
    fd = os.open(p, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        os.write(fd, line)
    finally:
        os.close(fd)
    return row


def read(path: Path | None = None) -> list[dict]:
    p = LEDGER if path is None else Path(path)
    if not p.exists():
        return []
    return [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
