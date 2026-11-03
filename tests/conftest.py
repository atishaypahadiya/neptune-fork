"""Shared fixtures."""

import json
from pathlib import Path

import numpy as np
import pytest

from fork.data import horizons

M4_FIXTURE = Path(__file__).parent / "data" / "horizons_m4.json"


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
