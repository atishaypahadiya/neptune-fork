"""M6e (docs/notes/M6e.md): the time and frame tests of § 2, the levels of § 3, the leak
check of § 6, and the recorded run. Posterior draws are not touched here."""

import ast
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fork.export import m6e

SRC = Path(m6e.__file__).resolve().parents[2]


def test_the_instants_follow_from_the_table():
    obs = m6e.lalande()
    for date, iso in m6e.INSTANTS_TDB.items():
        assert m6e.instant_from_paris_mean_time(date, obs[date]["paris_mean_time"]) == iso


def test_time_and_frame_checks_pass():
    """Meridian transits: sidereal time equals the observed right ascension within a minute.
    The chain on Uranus's cached vector reproduces JPL's apparent place within 1 arcsecond."""
    c = m6e.checks()
    for date in m6e.INSTANTS_TDB:
        assert abs(c["transit"][date]["hour_angle_min"]) <= m6e.TOL_TRANSIT_MIN
        f = c["frame"][date]
        assert abs(f["ra_minus_jpl_arcsec_on_sky"]) <= m6e.TOL_FRAME_ARCSEC
        assert abs(f["dec_minus_jpl_arcsec"]) <= m6e.TOL_FRAME_ARCSEC
    assert c["pass"]


def test_ecliptic_longitude_of_a_point_on_the_ecliptic():
    eps = math.radians(23.5)
    for lon in (10.0, 140.0, 215.0, 359.0):
        lam = math.radians(lon)
        ra = math.degrees(math.atan2(math.sin(lam) * math.cos(eps), math.cos(lam))) % 360.0
        dec = math.degrees(math.asin(math.sin(lam) * math.sin(eps)))
        assert m6e.ecliptic_longitude_deg(ra, dec, eps) == pytest.approx(lon, abs=1e-9)


def test_central_level_on_known_draws():
    x = np.linspace(-10.0, 10.0, 2001) + 359.0  # straddles 0°
    assert m6e.central_level(x, 359.0) == pytest.approx(0.0, abs=2e-3)
    assert m6e.central_level(x, 359.0 + 6.8) == pytest.approx(0.68, abs=2e-3)
    assert m6e.central_level(x, 359.0 - 9.5) == pytest.approx(0.95, abs=2e-3)
    assert m6e.central_level(x, 40.0) == pytest.approx(1.0)


@pytest.mark.parametrize(
    ("level", "name"),
    [
        (0.30, "inside 68 %"),
        (0.66, "inside 68 %"),
        (0.675, "at the boundary (inside 68 % or inside 95 %)"),
        (0.69, "at the boundary (inside 68 % or inside 95 %)"),
        (0.80, "inside 95 %"),
        (0.955, "at the boundary (inside 95 % or outside)"),
        (0.97, "outside"),
    ],
)
def test_level_classes(level, name):
    assert m6e.level_class(level) == name


def test_interval_unwraps_about_the_centre():
    x = (np.linspace(-5.0, 5.0, 1001) + 360.0) % 360.0
    q = m6e.interval(x)
    assert 0.0 <= q["median"] < 360.0 + 1e-9  # given on the circle, so 0° may read as 360°
    assert (q["median"] + 180.0) % 360.0 - 180.0 == pytest.approx(0.0, abs=1e-6)
    assert q["lo68"] < q["median"] < q["hi68"]
    assert q["width68"] == pytest.approx(6.8, abs=0.02)


def test_prior_draws_reach_1795():
    """The predictive step on draws from the fit's prior: every draw gives a place, and the
    prior's longitudes cover the circle."""
    from fork.export import m6c

    _, ctx, _, pq, box = m6c.dataset()
    cells = [
        SimpleNamespace(space=m6c.job_space(m6c.Job("headline", k), pq, box)) for k in range(4)
    ]
    w = [c.space.slice_weight for c in cells]
    assert sum(w) == pytest.approx(1.0)
    q = m6e.prior_draws("headline", cells, n=60)
    epochs = {d: m6e.Epoch(d) for d in m6e.INSTANTS_TDB}
    p = m6e.predictive(cells[0].space, ctx, q, epochs)
    lon = p["1795-05-08"]["lon_deg"]
    assert p["n_failed"] == 0 and len(lon) == 60
    assert m6e.interval(lon)["width95"] > 180.0
    s = m6e.summarize_draws(p, m6e.lalande())
    assert 0.0 <= s["1795-05-08"]["central_level_ecl_lon"] <= 1.0


def _names_the_table(path: Path) -> bool:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return any(isinstance(n, ast.Constant) and isinstance(n.value, str)
               and "precover" in n.value.lower() and len(n.value) < 200
               for n in ast.walk(tree))  # fmt: skip


def test_only_this_module_reads_the_precovery_table_and_the_fit_cannot_reach_it():
    from test_a6d import reachable_modules

    readers = sorted(
        p.relative_to(SRC).as_posix() for p in SRC.rglob("*.py") if _names_the_table(p)
    )
    assert readers == ["fork/export/m6e.py"]
    assert "fork.export.m6e" not in reachable_modules()


def test_the_run_loads_the_posteriors_of_the_plan():
    plan = (m6e.ROOT / "docs" / "notes" / "M6e.md").read_text(encoding="utf-8")
    for rel, digest in m6e.POSTERIOR_SHA256.items():
        assert f"`{rel}` | `{digest}`" in plan
    assert m6e.LALANDE_SHA256 in plan
    for iso in m6e.INSTANTS_TDB.values():
        assert iso.replace("T", " ") in plan
    assert json.loads(m6e.JPL_URANUS.read_text(encoding="utf-8"))["delta_t_s"] == m6e.DELTA_T_S


def test_the_recorded_run():
    r = json.loads(m6e.OUT_JSON.read_text(encoding="utf-8"))
    assert r["checks"]["pass"] and r["n_draws"] == m6e.N_DRAWS
    want = {"headline": (0.324, 0.324), "student_t": (0.2117, 0.2127)}
    for model, levels in want.items():
        e = r["models"][model]
        assert e["reading"] == "inside 68 %"
        assert abs(e["galle_night_median_deg"] - m6e.PUBLISHED_GALLE_MEDIAN_DEG[model]) <= 0.2
        for date, level in zip(m6e.INSTANTS_TDB, levels, strict=True):
            p = e["posterior"][date]
            assert p["central_level_ecl_lon"] == pytest.approx(level, abs=5e-4)
            assert m6e.level_class(p["central_level_ecl_lon"]) == p["class"]
            assert abs(p["median_minus_observed_deg"]["ecl_lon"]) < 0.35
            assert e["width_ratio_posterior_to_prior"][date]["width68"] < 0.01
        assert e["posterior"]["n_failed"] == 0 and e["prior_predictive"]["n_failed"] == 0
