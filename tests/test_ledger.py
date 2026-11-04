"""Compute ledger: rows appended per job, and the M5b job hook (run_nested faked)."""

from types import SimpleNamespace

import numpy as np
import pytest

from fork import ledger


def test_append_and_read_round_trip(tmp_path):
    p = tmp_path / "ledger.jsonl"
    ledger.append("a", "nested", n_likelihood_calls=10, n_integrations=10, sim_years=1560.0,
                  wall_s=2.0, cpu_s=7.5, n_workers=4, path=p)  # fmt: skip
    ledger.append("b", "profile", n_likelihood_calls=None, n_integrations=None, sim_years=None,
                  wall_s=60.0, cpu_s=None, n_workers=12, estimated=True, note="x", path=p)  # fmt: skip
    rows = ledger.read(p)
    assert [r["job"] for r in rows] == ["a", "b"]
    assert rows[0]["sim_years"] == pytest.approx(1560.0) and not rows[0]["estimated"]
    assert rows[1]["estimated"] and rows[1]["note"] == "x" and rows[1]["cpu_s"] is None
    assert ledger.cpu_now() > 0


def test_sim_years_spans_epoch_to_both_ends():
    epoch = 2371628.5
    jd = np.array([epoch - 365.25 * 91, epoch + 365.25 * 65])
    assert ledger.sim_years_per_integration(jd, epoch) == pytest.approx(156.0)


def test_m5b_job_appends_a_ledger_row(tmp_path, monkeypatch):
    from fork.export import m5b
    from fork.inference import nested

    monkeypatch.setattr(m5b, "CACHE", tmp_path)
    monkeypatch.setattr(ledger, "LEDGER", tmp_path / "ledger.jsonl")
    fake = SimpleNamespace(space=None, results={}, nlive=1000, dlogz=0.1, wall_s=1.0,
                           ncall=12345, niter=678, logz=-73.0, logzerr=0.2)  # fmt: skip
    monkeypatch.setattr(nested, "run_nested", lambda *a, **k: fake)
    monkeypatch.setattr(nested.NestedResult, "save", lambda self, path: path.write_text("x"))
    m5b.do_job(m5b.Job(3, "free_a_e0.1", 1), n_workers=4)
    (row,) = ledger.read()
    assert row["job"] == "m5b/seed_3/free_a_e0.1/c1" and row["kind"] == "nested"
    assert row["n_likelihood_calls"] == row["n_integrations"] == 12345
    assert 150 < row["sim_years"] / 12345 < 157  # 1690–1846 observation span, ≈ 155 yr
    assert row["n_workers"] == 4 and not row["estimated"] and row["niter"] == 678
