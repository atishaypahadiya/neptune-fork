"""Compute accounting: one JSON line per job in ``data/cache/compute_ledger.jsonl``.

Fields: job, kind, n_likelihood_calls, n_integrations, sim_years, wall_s, cpu_s, n_workers
(plus ``estimated`` and optional notes). One likelihood call is one forward-model evaluation,
i.e. one N-body integration from the 1781 epoch out to every observation epoch, so
``sim_years`` = n_integrations × the span the observations cover (1690–1846 ≈ 156 yr).
``cpu_s`` is this process plus its reaped children (pool workers). Totals:
``python -m fork.export.ledger`` → ``site/data/ledger.json``.

Runners open a record with :func:`start` (or the :func:`run` context manager) when a job begins
and close it with :meth:`Run.done`; a function that calls :func:`start` is decorated with
:func:`recorded`. A record that is not closed by ``done`` still writes its row, with the CPU spent
so far and an ``outcome``:

- ``done`` (or ``gate passed`` / ``gate failed`` for a gate's verdict step): closed normally;
- ``stopped``: SIGTERM or Ctrl-C. On SIGTERM the process's pool workers are ended and reaped
  first, so their CPU is counted; then the signal takes its default course (exit status −15);
- ``failed``: an uncaught exception (its type in ``error``);
- ``unfinished``: the interpreter exited with the record still open.

SIGKILL and power loss write nothing. Rows written before 2026-11-15 have no ``outcome``; they
were written only on success. Synthetic gates, capped smoke runs and runs into a scratch cache
write rows too (flagged ``synthetic``, ``capped`` or ``scratch`` where the runner knows it).

A record's CPU includes the children the process reaps, so a driver whose child processes write
their own rows (one process per run) opens no record of its own.
"""

from __future__ import annotations

import atexit
import contextlib
import functools
import json
import multiprocessing
import os
import resource
import signal
import sys
import threading
import time
from collections.abc import Iterator
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "data" / "cache" / "compute_ledger.jsonl"
DAYS_PER_YEAR = 365.25
OUTCOMES = ("done", "gate passed", "gate failed", "stopped", "failed", "unfinished")


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


# ---------------------------------------------------------------------------- records

_OPEN: list[Run] = []
_HOOKS = {"installed": False, "pid": None}


class Run:
    """One job's ledger record, open from :func:`start` until :meth:`done` (or until the
    process stops, fails or exits, which close it with that outcome)."""

    def __init__(self, job: str, kind: str, n_workers: int | None, path: Path | None,
                 fields: dict) -> None:  # fmt: skip
        self.job, self.kind, self.n_workers, self.path = job, kind, n_workers, path
        self.fields = dict(fields)
        self.t0, self.c0, self.pid = time.time(), cpu_now(), os.getpid()
        self.row: dict | None = None

    @property
    def closed(self) -> bool:
        return self.row is not None

    def set(self, **fields) -> None:
        """Add or replace fields. A callable value is called when the row is written (so a
        stopped run can report, for instance, the sampler's calls so far)."""
        self.fields.update(fields)

    def done(self, outcome: str = "done", **fields) -> dict:
        if outcome not in OUTCOMES:
            raise ValueError(f"unknown outcome {outcome!r}")
        return self._close(outcome, fields)

    def discard(self) -> None:
        """Close without a row: for a call that found nothing to do (all its work cached)."""
        self.row = {}
        if self in _OPEN:
            _OPEN.remove(self)

    def _close(self, outcome: str, fields: dict) -> dict:
        if self.closed:
            raise RuntimeError(f"ledger record {self.job} is already closed")
        vals = {}
        for k, v in {**self.fields, **fields}.items():
            if callable(v):
                try:
                    v = v()
                except Exception:  # a stopped run may lack what v reads
                    v = None
            vals[k] = v
        counts = {k: vals.pop(k, None) for k in ("n_likelihood_calls", "n_integrations",
                                                 "sim_years")}  # fmt: skip
        self.row = append(self.job, self.kind, **counts, wall_s=time.time() - self.t0,
                          cpu_s=cpu_now() - self.c0, n_workers=self.n_workers,
                          path=self.path, outcome=outcome, **vals)  # fmt: skip
        if self in _OPEN:
            _OPEN.remove(self)
        return self.row


def start(job: str, kind: str, *, n_workers: int | None, path: Path | None = None,
          **fields) -> Run:  # fmt: skip
    """Open a record for a job that starts now. Fields as for :func:`append`; any of them can
    also be given later with :meth:`Run.set` or :meth:`Run.done`."""
    _install_hooks()
    rec = Run(job, kind, n_workers, None if path is None else Path(path), fields)
    _OPEN.append(rec)
    return rec


@contextlib.contextmanager
def run(job: str, kind: str, *, n_workers: int | None, path: Path | None = None,
        **fields) -> Iterator[Run]:  # fmt: skip
    """``with ledger.run(...) as rec:`` — ``done`` on a normal exit (unless the body closed
    it), ``failed`` or ``stopped`` on an exception, which is re-raised."""
    rec = start(job, kind, n_workers=n_workers, path=path, **fields)
    try:
        yield rec
    except BaseException as e:
        if not rec.closed:
            rec._close(_outcome_of(e), {"error": type(e).__name__})
        raise
    if not rec.closed:
        rec.done()


def recorded(fn):
    """Decorator for a function that opens records with :func:`start`: when it returns, any of
    its records still open are closed as ``done`` (early returns, such as capped runs, keep
    their row); when it raises, as ``failed`` or ``stopped``. So a caller that catches the
    exception and goes on never leaves a record open."""

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        before = {id(r) for r in _OPEN}
        try:
            out = fn(*args, **kwargs)
        except BaseException as e:
            _close_new(before, _outcome_of(e), {"error": type(e).__name__})
            raise
        _close_new(before, "done", {})
        return out

    return wrapper


def _close_new(before: set[int], outcome: str, fields: dict) -> None:
    for rec in [r for r in _mine() if id(r) not in before]:
        rec._close(outcome, fields)


def _outcome_of(e: BaseException) -> str:
    return "stopped" if isinstance(e, KeyboardInterrupt) else "failed"


def _mine() -> list[Run]:
    return [r for r in _OPEN if r.pid == os.getpid() and not r.closed]


def _close_all(outcome: str, fields: dict) -> None:
    for rec in _mine():
        try:
            rec._close(outcome, fields)
        except Exception:  # never mask the original exit
            pass


def _end_children() -> None:
    """End and reap this process's multiprocessing children, so their CPU is counted."""
    kids = multiprocessing.active_children()
    for p in kids:
        p.terminate()
    for p in kids:
        p.join(timeout=10)


def _on_sigterm(signum: int, frame) -> None:
    if _mine():
        _end_children()
        _close_all("stopped", {"signal": "SIGTERM"})
    prev = _HOOKS.get("prev_sigterm", signal.SIG_DFL)
    if callable(prev):
        prev(signum, frame)
        return
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGTERM)


def _excepthook(etype, value, tb) -> None:
    _close_all(_outcome_of(value), {"error": etype.__name__})
    _HOOKS["prev_excepthook"](etype, value, tb)


def _at_exit() -> None:
    _close_all("unfinished", {})


def _install_hooks() -> None:
    """Once per process: SIGTERM handler (main thread only, and only over the default
    handler), exception hook and exit hook."""
    if _HOOKS["installed"] and _HOOKS["pid"] == os.getpid():
        return
    _HOOKS.update(installed=True, pid=os.getpid())
    if threading.current_thread() is threading.main_thread():
        prev = signal.getsignal(signal.SIGTERM)
        if prev in (signal.SIG_DFL, None) or (callable(prev) and prev is not _on_sigterm):
            _HOOKS["prev_sigterm"] = signal.SIG_DFL if prev is None else prev
            signal.signal(signal.SIGTERM, _on_sigterm)
    if sys.excepthook is not _excepthook:
        _HOOKS["prev_excepthook"] = sys.excepthook
        sys.excepthook = _excepthook
    atexit.register(_at_exit)
