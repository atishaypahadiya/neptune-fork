"""Shared fixtures.

Horizons in the default test job (``-m "not slow and not network"``): every Horizons state comes
from ``tests/fixtures/horizons/`` (a session copy, so nothing writes into the fixture), and any
live fetch raises. Session-wide, so module- and session-scoped fixtures are covered too.
- ``network`` tests get the live client and the normal cache back (they exercise the client).
- ``slow`` tests get the normal cache back: they also need other local caches (posteriors), and
  they are not part of the default job.
- Modules using ``offline_horizons`` keep their own M4 cache, as before.
Regenerate the fixture with ``uv run python tests/fixtures/build_horizons_fixture.py``.
"""

import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from fork.data import horizons

M4_FIXTURE = Path(__file__).parent / "data" / "horizons_m4.json"
HORIZONS_FIXTURE = Path(__file__).parent / "fixtures" / "horizons"
RECORD_ENV = "FORK_HORIZONS_FIXTURE_RECORD"  # "1": build the fixture instead of reading it
_LIVE_FETCH = horizons.fetch_from_horizons
_OUTER_CACHE_ENV = os.environ.get("FORK_CACHE_DIR")  # as outside pytest
_ID_TO_BODY = {b.horizons_id: name for name, b in horizons.BODIES.items()}


def _fixture_miss(horizons_id: str, jd_tdb: np.ndarray) -> np.ndarray:
    raise AssertionError(
        f"Horizons id={horizons_id}: {jd_tdb.size} epoch(s) not in tests/fixtures/horizons "
        f"(first {jd_tdb[:3]}). The default test job makes no live calls; regenerate the "
        "fixture with `uv run python tests/fixtures/build_horizons_fixture.py`."
    )


def _serve_from_local_cache(horizons_id: str, jd_tdb: np.ndarray) -> np.ndarray:
    """Record mode: the states from the local cache (data/cache), live only for a local miss."""
    body = _ID_TO_BODY[horizons_id]
    jd, st = horizons._load_cache(horizons._cache_path(body, horizons.DEFAULT_CACHE_DIR))
    if jd.size:
        i = np.clip(np.searchsorted(jd, jd_tdb), 0, jd.size - 1)
        if np.array_equal(jd[i], jd_tdb):
            return st[i]
    return _LIVE_FETCH(horizons_id, jd_tdb)


@pytest.fixture(scope="session", autouse=True)
def horizons_fixture_cache(tmp_path_factory):
    """FORK_CACHE_DIR -> a copy of tests/fixtures/horizons for the whole session; fetches raise.

    In record mode the copy starts empty, fetches are served from the local cache, and at the
    end the epochs the tests asked for become the new fixture.
    """
    directory = tmp_path_factory.mktemp("horizons_fixture")
    recording = os.environ.get(RECORD_ENV) == "1"
    if not recording:
        for f in sorted(HORIZONS_FIXTURE.glob("horizons_*.npz")):
            shutil.copyfile(f, directory / f.name)
    mp = pytest.MonkeyPatch()
    mp.setenv("FORK_CACHE_DIR", str(directory))
    mp.setattr(
        horizons, "fetch_from_horizons", _serve_from_local_cache if recording else _fixture_miss
    )
    yield directory
    mp.undo()
    if recording:
        HORIZONS_FIXTURE.mkdir(parents=True, exist_ok=True)
        for f in HORIZONS_FIXTURE.glob("horizons_*.npz"):
            f.unlink()
        for f in sorted(directory.glob("horizons_*.npz")):
            shutil.copyfile(f, HORIZONS_FIXTURE / f.name)


@pytest.fixture(autouse=True)
def _live_horizons_for_network_and_slow(request):
    """`network` and `slow` tests: the normal cache and the live client, as outside pytest."""
    node = request.node
    if "offline_horizons" in request.fixturenames or not (
        node.get_closest_marker("network") or node.get_closest_marker("slow")
    ):
        yield
        return
    mp = pytest.MonkeyPatch()
    if _OUTER_CACHE_ENV is None:
        mp.delenv("FORK_CACHE_DIR", raising=False)
    else:
        mp.setenv("FORK_CACHE_DIR", _OUTER_CACHE_ENV)
    mp.setattr(horizons, "fetch_from_horizons", _LIVE_FETCH)
    yield
    mp.undo()


@pytest.fixture(scope="module")
def offline_horizons(tmp_path_factory):
    """Horizons cache populated from tests/data/horizons_m4.json; any network fetch fails.

    Module-scoped so the no-network patch never leaks into the live `network` tests.

    Lets code that calls ``fork.data.horizons`` normally (ForwardContext.from_horizons,
    synthetic truth) run hermetically in CI.
    """
    directory = tmp_path_factory.mktemp("horizons_cache")
    data = json.loads(M4_FIXTURE.read_text())
    for body, rec in data["bodies"].items():
        horizons._save_cache(
            directory / f"horizons_{body}.npz",
            np.array(rec["jd_tdb"]),
            np.array(rec["state_ecl_j2000"]),
            horizons.BODIES[body].horizons_id,
        )

    def no_network(horizons_id, jd_tdb):
        raise AssertionError(f"offline test tried to fetch Horizons id={horizons_id} at {jd_tdb}")

    mp = pytest.MonkeyPatch()
    mp.setenv("FORK_CACHE_DIR", str(directory))
    mp.setattr(horizons, "fetch_from_horizons", no_network)
    yield directory
    mp.undo()
