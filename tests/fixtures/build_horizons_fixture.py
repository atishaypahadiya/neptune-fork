"""Rebuild tests/fixtures/horizons/: the JPL Horizons states the default test job reads.

Runs the default job (``-m "not slow and not network"``) in record mode (``tests/conftest.py``,
``FORK_HORIZONS_FIXTURE_RECORD=1``): every Horizons state a test asks for is served from the
local cache (``data/cache/horizons_<body>.npz``, filled by ``python -m fork.data.horizons``),
or fetched live for an epoch the local cache lacks, and exactly those epochs are written to
``tests/fixtures/horizons/horizons_<body>.npz`` in the cache's own format. Then the job is run
again from the new fixture with live fetches disabled, to check it is complete.

Run from the repo root: ``uv run python tests/fixtures/build_horizons_fixture.py``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "tests" / "fixtures" / "horizons"
SELECT = ["-m", "not slow and not network", "-q", "-p", "no:cacheprovider"]


def main() -> None:
    env = {**os.environ, "FORK_HORIZONS_FIXTURE_RECORD": "1"}
    rec = subprocess.run([sys.executable, "-m", "pytest", *SELECT], cwd=ROOT, env=env)
    if rec.returncode:
        sys.exit(f"record run failed ({rec.returncode})")
    for f in sorted(FIXTURE.glob("horizons_*.npz")):
        with np.load(f) as z:
            jd = z["jd_tdb"]
            print(
                f"{f.name}: {jd.size} epochs, JD {jd.min():.4f}-{jd.max():.4f}, "
                f"id {z['horizons_id']}, center {z['center']}, refplane {z['refplane']}"
            )
    env.pop("FORK_HORIZONS_FIXTURE_RECORD")
    check = subprocess.run([sys.executable, "-m", "pytest", *SELECT], cwd=ROOT, env=env)
    sys.exit(check.returncode)


if __name__ == "__main__":
    main()
