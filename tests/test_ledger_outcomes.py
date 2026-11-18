"""Every run writes a ledger row with its CPU and outcome: done, gate passed or failed, stopped
(SIGTERM, Ctrl-C), failed (an exception) or unfinished (exit with the record open)."""

import ast
import json
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fork import ledger

SRC = Path(__file__).resolve().parents[1] / "src" / "fork"


def test_done_row_has_outcome_wall_and_cpu(tmp_path):
    p = tmp_path / "l.jsonl"
    rec = ledger.start("x/1", "nested", n_workers=1, path=p, nlive=10)
    sum(i * i for i in range(200_000))
    rec.done(n_likelihood_calls=5, niter=3)
    (row,) = ledger.read(p)
    assert row["outcome"] == "done" and row["job"] == "x/1" and row["n_workers"] == 1
    assert row["n_likelihood_calls"] == 5 and row["n_integrations"] is None
    assert row["nlive"] == 10 and row["niter"] == 3
    assert row["wall_s"] >= 0 and row["cpu_s"] >= 0
    with pytest.raises(RuntimeError):
        rec.done()
    bad = ledger.start("x/2", "nested", n_workers=1, path=p)
    with pytest.raises(ValueError):
        bad.done("finished")
    bad.discard()
    assert len(ledger.read(p)) == 1 and not ledger._mine()


def test_run_context_records_failed_and_stopped(tmp_path):
    p = tmp_path / "l.jsonl"
    with pytest.raises(ZeroDivisionError), ledger.run("x/f", "fit", n_workers=1, path=p):
        1 / 0  # noqa: B018
    with pytest.raises(KeyboardInterrupt), ledger.run("x/s", "fit", n_workers=1, path=p):
        raise KeyboardInterrupt
    with ledger.run("x/g", "gate", n_workers=2, path=p) as rec:
        rec.done("gate failed")
    with ledger.run("x/d", "fit", n_workers=1, path=p):
        pass
    rows = {r["job"]: r for r in ledger.read(p)}
    assert rows["x/f"]["outcome"] == "failed" and rows["x/f"]["error"] == "ZeroDivisionError"
    assert rows["x/s"]["outcome"] == "stopped"
    assert rows["x/g"]["outcome"] == "gate failed" and rows["x/d"]["outcome"] == "done"
    assert not ledger._mine()


def test_recorded_closes_on_early_return_and_on_a_caught_exception(tmp_path):
    p = tmp_path / "l.jsonl"

    @ledger.recorded
    def capped(n: int) -> int:
        rec = ledger.start("x/capped", "nested", n_workers=1, path=p, capped=True)
        if n < 10:
            return n  # an early return, as a smoke run's
        rec.done()
        return n

    @ledger.recorded
    def broken() -> None:
        ledger.start("x/broken", "nested", n_workers=1, path=p)
        raise RuntimeError("lost a chunk")

    assert capped(3) == 3
    try:  # a caller that goes on after a failure
        broken()
    except RuntimeError:
        pass
    rows = ledger.read(p)
    assert [(r["job"], r["outcome"]) for r in rows] == [
        ("x/capped", "done"),
        ("x/broken", "failed"),
    ]
    assert rows[0]["capped"] and rows[1]["error"] == "RuntimeError"
    assert not ledger._mine()


def test_callable_fields_are_read_when_the_row_is_written(tmp_path):
    p = tmp_path / "l.jsonl"
    box = {"n": 0}
    rec = ledger.start("x/c", "nested", n_workers=1, path=p)
    rec.set(n_likelihood_calls=lambda: box["n"], broken=lambda: box["missing"])
    box["n"] = 41
    row = rec.done()
    assert row["n_likelihood_calls"] == 41 and row["broken"] is None


def test_discard_writes_nothing(tmp_path):
    p = tmp_path / "l.jsonl"
    ledger.start("x/none", "reweight", n_workers=4, path=p).discard()
    assert ledger.read(p) == [] and not ledger._mine()


CHILD = textwrap.dedent(
    """
    import multiprocessing as mp, sys, time
    from pathlib import Path
    from fork import ledger

    def burn(_):
        t = time.process_time()
        while time.process_time() - t < 60:
            pass

    if __name__ == "__main__":
        mode, path, flag = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
        rec = ledger.start("child/" + mode, "nested", n_workers=2, path=path, synthetic=True)
        if mode == "fail":
            raise ValueError("boom")
        if mode == "exit":
            sys.exit(0)
        pool = mp.get_context("spawn").Pool(2)
        pool.map_async(burn, range(2))
        t = time.process_time()
        while time.process_time() - t < 0.3:
            pass
        flag.write_text("ready")
        time.sleep(60)
    """
)


def _child(tmp_path, mode: str) -> tuple[subprocess.Popen, Path, Path]:
    script = tmp_path / "child.py"
    script.write_text(CHILD)
    path, flag = tmp_path / "l.jsonl", tmp_path / "ready"
    proc = subprocess.Popen([sys.executable, str(script), mode, str(path), str(flag)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # fmt: skip
    return proc, path, flag


def test_sigterm_writes_a_stopped_row_with_the_pool_workers_cpu(tmp_path):
    proc, path, flag = _child(tmp_path, "stop")
    t_end = time.time() + 60
    while not flag.exists() and time.time() < t_end and proc.poll() is None:
        time.sleep(0.05)
    assert flag.exists(), proc.communicate()[1].decode()[-2000:]
    time.sleep(1.5)  # the two workers burn CPU
    proc.send_signal(signal.SIGTERM)
    assert proc.wait(timeout=30) == -signal.SIGTERM  # the signal still ends the process
    (row,) = ledger.read(path)
    assert row["outcome"] == "stopped" and row["signal"] == "SIGTERM" and row["synthetic"]
    # parent ≥ 0.3 s, each worker about 1.5 s: the workers were reaped and counted
    assert row["cpu_s"] > 1.5, row


@pytest.mark.parametrize(("mode", "outcome", "code"), [("fail", "failed", 1),
                                                      ("exit", "unfinished", 0)])  # fmt: skip
def test_uncaught_exception_and_exit_with_an_open_record(tmp_path, mode, outcome, code):
    proc, path, _ = _child(tmp_path, mode)
    assert proc.wait(timeout=60) == code
    (row,) = ledger.read(path)
    assert row["outcome"] == outcome
    if mode == "fail":
        assert row["error"] == "ValueError"


def test_totals_by_outcome_count_old_rows_as_done():
    from fork.export.ledger import totals

    rows = [
        {"job": "a", "kind": "nested", "cpu_s": 10.0},
        {"job": "b", "kind": "nested", "cpu_s": 5.0, "outcome": "stopped"},
        {"job": "c", "kind": "gate", "cpu_s": 1.0, "outcome": "gate failed"},
        {"job": "d", "kind": "nested", "cpu_s": 2.0, "outcome": "done"},
    ]
    t = totals(rows)["by_outcome"]
    assert list(t) == ["done", "gate failed", "stopped"]
    assert t["done"]["cpu_s"] == pytest.approx(12.0) and t["done"]["n_rows"] == 2
    assert t["stopped"]["cpu_s"] == pytest.approx(5.0)


# ---------------------------------------------------------------------------- runners


def test_v14_run_that_fails_records_its_segment(tmp_path, monkeypatch):
    """A synthetic-gate run whose sampler dies writes a failed row with the calls so far; a
    resumed run counts only the calls of its own segment."""
    import dynesty

    from fork.export import v14_fits as V
    from fork.inference import nested

    monkeypatch.setattr(V, "CACHE", tmp_path)
    fit = V.FITS["adams"]
    space = SimpleNamespace(unit_names=["u"] * (4 + fit.n_phi), periodic_unit_dims=[5])
    monkeypatch.setattr(V, "problem", lambda *a, **k: (None, None, None, None, object()))
    monkeypatch.setattr(V, "cell_space", lambda *a, **k: space)

    class Sampler:
        def __init__(self, *a, **k):
            self.ncall, self.it, self.results = 0, 0, {}

        def run_nested(self, resume=False, **kw):
            self.ncall += 700
            if not resume:
                raise RuntimeError("worker lost")

        @classmethod
        def restore(cls, path):
            s = cls()
            s.ncall = 300
            return s

    monkeypatch.setattr(dynesty, "NestedSampler", Sampler)
    fake = SimpleNamespace(logz=-1.0, logzerr=0.1, niter=9, ncall=1000, save=lambda p: None)
    monkeypatch.setattr(nested, "NestedResult", lambda *a, **k: fake)
    with pytest.raises(RuntimeError):
        V.do_job(fit, 0, 1, synthetic=True)
    ckpt = V.job_path(0, 1, ".ckpt", V.fit_cache(fit, True))
    ckpt.write_text("checkpoint")
    V.do_job(fit, 0, 1, synthetic=True)
    first, second = ledger.read()
    assert first["job"] == second["job"] == "v14/adams/g0/c0/s1"
    assert first["outcome"] == "failed" and first["n_likelihood_calls"] == 700
    assert first["synthetic"] and first["scratch"] and first["error"] == "RuntimeError"
    assert second["outcome"] == "done" and second["resumed"]
    assert second["n_likelihood_calls"] == 1000 - 300 and second["ncall_total"] == 1000


def test_mercury_study_writes_a_row_per_step(monkeypatch):
    from fork.phase2 import identifiability_study as S

    monkeypatch.setattr(S, "prepare", lambda: {})
    monkeypatch.setattr(S, "validate", lambda n: {"pass": False})
    monkeypatch.setattr(S, "run", lambda n: pytest.fail("run after a failed gate"))
    monkeypatch.setattr(S, "summarize", lambda: {})
    S.main(["all", "--n-workers", "3"])
    rows = [(r["job"].split("/")[-1], r["outcome"], r["n_workers"]) for r in ledger.read()]
    assert rows == [("prepare", "done", 1), ("validate", "gate failed", 3),
                    ("summarize", "done", 1)]  # fmt: skip


# ---------------------------------------------------------------------------- coverage

# Command-line modules that run no fit, sampler, integration or Monte Carlo, with the reason.
NOT_RUNNERS = {
    "data/horizons.py": "fetches and caches reference states",
    "data/transcription.py": "compares two entries of a table",
    "data/spotcheck.py": "builds a spot-check kit",
    "data/observations.py": "writes the observation table",
    "data/calibration.py": "residuals against reference positions, seconds",
    "phase2/reference.py": "fetches reference states",
    "export/to_json.py": "site data from saved results",
    "export/s1_site.py": "site data from saved results",
    "export/m6e_site.py": "site data from saved results",
    "export/m6d_site.py": "site data from saved results",
    "export/p1b_site.py": "site data from saved results",
    "export/ledger.py": "totals of the ledger itself",
}


# Code whose SHA-256 a frozen plan records: left byte-identical, so not converted.
FROZEN_CODE = {
    "export/p1b_adams.py": "docs/notes/P1B_adams_fit.md (writes its rows on success only)",
    "phase2/identifiability.py": "docs/phase2/identifiability_study.md § 2 (seconds, no row)",
}


def _modules() -> list[Path]:
    return sorted(p for p in SRC.rglob("*.py") if "__pycache__" not in p.parts)


def test_every_runner_records_and_only_the_ledger_appends():
    missing, direct, undecorated = [], [], []
    for p in _modules():
        rel = p.relative_to(SRC).as_posix()
        if rel in FROZEN_CODE:
            continue
        src = p.read_text()
        tree = ast.parse(src)
        if rel not in ("ledger.py", "export/ledger.py") and "ledger.append(" in src:
            direct.append(rel)
        has_main = any(isinstance(n, ast.FunctionDef) and n.name == "main" for n in tree.body)
        if (
            has_main
            and rel not in NOT_RUNNERS
            and not ("ledger.start(" in src or "ledger.run(" in src)
        ):
            missing.append(rel)
        for f in ast.walk(tree):
            if not isinstance(f, ast.FunctionDef):
                continue
            inner = {id(n) for g in ast.walk(f) if g is not f and isinstance(g, ast.FunctionDef)
                     for n in ast.walk(g)}  # fmt: skip
            opens = any(isinstance(n, ast.Call) and ast.unparse(n.func) == "ledger.start"
                        and id(n) not in inner for n in ast.walk(f))  # fmt: skip
            decorated = any(ast.unparse(d) == "ledger.recorded" for d in f.decorator_list)
            if opens and not decorated:
                undecorated.append(f"{rel}:{f.name}")
    assert not missing, f"runners with no ledger record: {missing}"
    assert not direct, f"append the ledger through ledger.start/run only: {direct}"
    assert not undecorated, f"functions that open a record need @ledger.recorded: {undecorated}"
    assert all((SRC / rel).exists() for rel in (*NOT_RUNNERS, *FROZEN_CODE))


def _rebinds(f: ast.FunctionDef, names: set[str], opens: list) -> list[int]:
    """Lines in ``f`` that bind one of ``names`` again (assignment, loop or with target)."""
    out = []
    for n in ast.walk(f):
        if isinstance(n, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            targets = n.targets if isinstance(n, ast.Assign) else [n.target]
        elif isinstance(n, (ast.For, ast.comprehension)):
            targets = [n.target]
        elif isinstance(n, ast.With):
            targets = [i.optional_vars for i in n.items if i.optional_vars is not None]
        else:
            continue
        if n in opens:
            continue
        for t in targets:
            if any(isinstance(x, ast.Name) and x.id in names for x in ast.walk(t)):
                out.append(n.lineno)
    return out


def test_a_record_is_not_overwritten_before_it_is_closed():
    """A name that holds a ledger record is bound once in its function: a later
    ``rec = {...}`` would make the closing ``rec.done(...)`` fail."""
    bad = []
    for p in _modules():
        rel = p.relative_to(SRC).as_posix()
        if rel in FROZEN_CODE:
            continue
        for f in ast.walk(ast.parse(p.read_text())):
            if not isinstance(f, ast.FunctionDef):
                continue
            opens = [n for n in ast.walk(f) if isinstance(n, ast.Assign)
                     and isinstance(n.value, ast.Call)
                     and ast.unparse(n.value.func) == "ledger.start"]  # fmt: skip
            names = {ast.unparse(t) for n in opens for t in n.targets}
            bad += [f"{rel}:{ln} {f.name}" for ln in _rebinds(f, names, opens)]
    assert not bad, f"ledger record overwritten: {bad}"


@pytest.mark.parametrize(("rel", "plan", "digest"), [
    ("export/p1b_adams.py", "docs/notes/P1B_adams_fit.md", "c9c4d494ac441a4d"),
    ("phase2/identifiability.py", "docs/phase2/identifiability_study.md", "14676bb8d02a945d"),
])  # fmt: skip
def test_frozen_code_is_unchanged(rel, plan, digest):
    import hashlib

    assert hashlib.sha256((SRC / rel).read_bytes()).hexdigest().startswith(digest)
    assert digest in (SRC.parents[1] / plan).read_text()


def test_rows_are_json_lines(tmp_path):
    p = tmp_path / "l.jsonl"
    ledger.start("x", "nested", n_workers=1, path=p, arr=int(np.int64(3))).done()
    assert json.loads(p.read_text().splitlines()[0])["arr"] == 3
    assert os.path.getsize(p) > 0
