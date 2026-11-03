"""site/data exports."""

import json
from pathlib import Path

import numpy as np

from fork.export import to_json

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).parent / "data" / "horizons_a1_1781-03-13.json"


def _a1_ics():
    data = json.loads(FIXTURE.read_text())
    assert data["jd_tdb"] == to_json.A1_EPOCH_JD_TDB
    st = np.array([data["state_ecl_j2000"][n] for n in to_json.A1_BODIES])
    from fork.model import bodies

    masses = np.array(
        [bodies.M_CENTRAL_MSUN] + [bodies.BODIES[n].mass_msun for n in to_json.A1_BODIES]
    )
    return st[:, :3].copy(), st[:, 3:].copy(), masses


def test_wh_convergence_data_schema_and_order():
    d = to_json.wh_convergence_data(*_a1_ics(), dt_days=(24.0, 12.0))
    assert d["dt_days"] == [24.0, 12.0] and d["default_dt_days"] == 12.0
    names = [s["scheme"] for s in d["series"]]
    assert names == ["wh_drift_outer", "wh_kick_outer"]
    assert [s["is_default"] for s in d["series"]] == [True, False]
    for s in d["series"]:
        err = s["max_abs_dlon_uranus_arcsec"]
        assert 3.5 < err[0] / err[1] < 4.5
    assert d["series"][0]["max_abs_dlon_uranus_arcsec"][1] < 0.01  # default @ 12 d


def test_committed_wh_convergence_json_is_current():
    """site/data/wh_convergence.json must equal a fresh export (regenerate with
    `uv run python -m fork.export.to_json` if this fails)."""
    committed = json.loads((ROOT / "site" / "data" / "wh_convergence.json").read_text())
    fresh = json.loads(json.dumps(to_json.wh_convergence_data(*_a1_ics())))
    assert committed == fresh


def test_site_data_under_5mb():
    total = sum(p.stat().st_size for p in (ROOT / "site" / "data").glob("*") if p.is_file())
    assert total < 5 * 1024 * 1024
