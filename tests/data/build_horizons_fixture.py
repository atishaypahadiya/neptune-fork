"""Rebuild tests/data/horizons_m4.json from data/cache (run `--fetch` first).

    uv run python -m fork.data.horizons --fetch
    uv run python tests/data/build_horizons_fixture.py

Public-domain JPL Horizons states, committed so M4/M5 tests run without network.
"""

import json
from pathlib import Path

import numpy as np

from fork.data.frames import jd_tdb_from_iso
from fork.data.horizons import (
    A1_EPOCH_ISO_TDB,
    GALLE_EPOCH_ISO_TDB,
    HORIZONS_CENTER,
    HORIZONS_REFPLANE,
    state_vectors_ecl_j2000,
)
from fork.data.synthetic import jd_tdb_jan1
from fork.model.bodies import INNER_BODIES

OUT = Path(__file__).parent / "horizons_m4.json"


def plan() -> dict[str, np.ndarray]:
    epoch = np.array([jd_tdb_from_iso(A1_EPOCH_ISO_TDB)])
    galle = np.array([jd_tdb_from_iso(GALLE_EPOCH_ISO_TDB)])
    out = {"uranus": np.union1d(jd_tdb_jan1(np.arange(1690, 1847)), epoch)}
    out |= {b: epoch for b in ("jupiter", "saturn")}
    out |= {b: np.union1d(epoch, galle) for b in ("neptune", *(b.name for b in INNER_BODIES))}
    out["earth"] = galle
    return out


def main() -> None:
    data = {
        "description": "JPL Horizons Sun-centred geometric states, ecliptic J2000, AU and AU/day, "
        "at TDB JDs. Public-domain JPL data, committed so M4/M5 tests run without network. "
        "Rebuild with tests/data/build_horizons_fixture.py.",
        "center": HORIZONS_CENTER,
        "refplane": HORIZONS_REFPLANE,
        "bodies": {},
    }
    for body, jd in plan().items():
        st = state_vectors_ecl_j2000(body, jd)
        data["bodies"][body] = {"jd_tdb": jd.tolist(), "state_ecl_j2000": st.tolist()}
    OUT.write_text(json.dumps(data, indent=0))
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
