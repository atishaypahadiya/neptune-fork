"""site/data exports."""

import json
from pathlib import Path

import numpy as np
import pytest

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


# ------------------------------------------------------------------------- M5a residuals


def _m5a_committed():
    return json.loads((ROOT / "site" / "data" / "m5a_fits.json").read_text())


def test_m5a_json_schema_and_summary():
    d = _m5a_committed()
    assert d["years"] == list(range(1690, 1846)) and d["seed"] == to_json.M5A_SEED
    assert [c["id"] for c in d["curves"]] == ["no_perturber", "le_verrier", "free_a"]
    assert all(len(c["residual_arcsec"]) == len(d["years"]) for c in d["curves"])
    n_obs = len(d["observations"]["year"])
    assert all(len(v) == n_obs for v in d["observations"]["residual_arcsec"].values())
    for mode in ("le_verrier", "free_a"):
        f = d["fits"][mode]
        assert len(f["theta"]) == 11 and f["n_starts"] == 24
        assert sum(f["status_counts"].values()) == 24 and f["n_distinct_minima"] == len(f["minima"])
        assert f["minima"][0]["chi2"] == f["chi2"]
        assert f["rms_post_1781_arcsec"] < 4.5  # A2b-style: within 1.5σ
    assert d["fits"]["le_verrier"]["a_au"] == 38.4


@pytest.mark.usefixtures("offline_horizons")
def test_m5a_json_curves_reproduce_from_committed_theta():
    """Each committed residual curve is truth − model(committed θ): cheap to recheck in CI."""
    from fork.data.synthetic import jd_tdb_jan1, true_uranus_lon_ecl_j2000_rad
    from fork.inference.fit import wrap_rad
    from fork.model.forward import ForwardContext, predict_uranus_lon

    d = _m5a_committed()
    ctx = ForwardContext.from_horizons(d["epoch_jd_tdb"])
    grid = jd_tdb_jan1(d["years"])
    truth = true_uranus_lon_ecl_j2000_rad(grid)
    for c in d["curves"][1:]:
        theta = np.array(d["fits"][c["id"]]["theta"])
        curve = wrap_rad(truth - predict_uranus_lon(theta, grid, ctx)) * to_json.RAD_TO_ARCSEC
        np.testing.assert_allclose(c["residual_arcsec"], curve, rtol=0, atol=2e-3)


@pytest.mark.slow
@pytest.mark.usefixtures("offline_horizons")
def test_m5a_json_is_current():
    """Fresh M5a export (~25 s) matches the committed file (regenerate with
    `uv run python -m fork.export.to_json`). Tolerances allow cross-platform LM roundoff."""
    import os

    fresh = json.loads(json.dumps(to_json.m5a_fits_data(n_workers=min(12, os.cpu_count() or 1))))
    committed = _m5a_committed()
    for mode in ("le_verrier", "free_a"):
        f, c = fresh["fits"][mode], committed["fits"][mode]
        assert f["n_distinct_minima"] == c["n_distinct_minima"]
        assert f["status_counts"] == c["status_counts"]
        np.testing.assert_allclose(f["theta"], c["theta"], rtol=1e-6, atol=1e-8)
        assert f["galle"]["dlon_vs_neptune_deg"] == pytest.approx(
            c["galle"]["dlon_vs_neptune_deg"], abs=1e-3
        )
    for cf, cc in zip(fresh["curves"], committed["curves"], strict=True):
        np.testing.assert_allclose(cf["residual_arcsec"], cc["residual_arcsec"], atol=1e-2)
