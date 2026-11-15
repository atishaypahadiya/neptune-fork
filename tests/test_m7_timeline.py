"""The pre-registration timeline's run status (docs/notes/scripts/m7_timeline.py)."""

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "docs" / "notes" / "scripts" / "m7_timeline.py"


@pytest.fixture(scope="module")
def T():
    spec = importlib.util.spec_from_file_location("m7_timeline", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["m7_timeline"] = mod
    spec.loader.exec_module(mod)
    return mod


LAUNCHED = {"time": "2026-11-09T16:27:37-08:00", "source": "x.log"}
FINISHED = {"time": "2026-11-10T02:40:42-08:00", "source": "x.log"}


def test_local_no_finish_line_is_running_or_stopped_by_process(T):
    seen = []

    def alive(pattern):
        seen.append(pattern)
        return pattern == "job-a"

    assert T.run_status(LAUNCHED, None, "job-a", local=True, alive=alive) == "running"
    assert T.run_status(LAUNCHED, None, "job-b", local=True, alive=alive) == "stopped"
    assert seen == ["job-a", "job-b"]
    assert T.run_status(LAUNCHED, FINISHED, "job-a", local=True, alive=alive) == "done"
    assert T.run_status(None, None, "job-a", local=True, alive=alive) == "not launched"


def test_without_process_view_no_finish_line_is_neither_running_nor_stopped(T):
    def alive(pattern):
        raise AssertionError("no process lookup where processes cannot be seen")

    assert T.run_status(LAUNCHED, None, "job", local=False, alive=alive) == "no finish line"
    assert T.run_status(LAUNCHED, FINISHED, "job", local=False, alive=alive) == "done"


@pytest.mark.parametrize("var", ["CI", "GITHUB_ACTIONS"])
def test_ci_cannot_see_processes(T, monkeypatch, var):
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.delenv("GITHUB_ACTIONS", raising=False)
    monkeypatch.setenv(var, "true")
    assert T.can_see_processes() is False


def test_arbiter_finish_is_in_its_driver_log(T, tmp_path, monkeypatch):
    """The arbiter's progress.log ends with each run's own "done:" line; "run finished" and
    "summarize" are written to driver.log (data/cache/nested/m8_arbiter, 10 November)."""
    d = tmp_path / "data" / "cache" / "nested" / "m8_arbiter"
    d.mkdir(parents=True)
    (d / "progress.log").write_text(
        "2026-11-09 16:27:37 [m8_arbiter/run2] start: ndim 11, nlive 4000\n"
        "2026-11-10 02:38:59 [m8_arbiter/run2] iter 197482 | dlogz:  0.108 >  0.100\n"
        "2026-11-10 02:39:54 [m8_arbiter/run2] done: 197792 it, ln Z -362.841 ± 0.126\n"
    )
    (d / "driver.log").write_text(
        "2026-11-09 16:27:36 run: 2 runs pending, budget 12 workers\n"
        "2026-11-10 02:40:12 done run2 (612.6 min)\n"
        "2026-11-10 02:40:42 run finished; failures: none; all done: True\n"
        "2026-11-10 03:10:32 summarize: test pass False\n"
        "2026-11-10 03:10:32 pipeline: done\n"
    )
    monkeypatch.setattr(T, "ROOT", tmp_path)
    launched = T.log_time(d / "progress.log", r"start")
    assert T.log_time(d / "progress.log", r"run finished") is None  # where it was looked for
    finished = T.log_time(d / "driver.log", r"^run finished")
    summary = T.log_time(d / "driver.log", r"^summarize")
    assert finished["time"] == "2026-11-10T02:40:42-08:00"
    assert summary["time"] == "2026-11-10T03:10:32-08:00"
    assert T.run_status(launched, finished, "x", local=False) == "done"
    # had it no finish line: its last log entry, from the log its launch was read from
    assert T.last_entry(launched)["time"] == "2026-11-10T02:39:54-08:00"
