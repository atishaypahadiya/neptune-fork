"""The synthetic identifiability tool (fork.phase2.identifiability): built and tested only."""

import numpy as np
import pytest

from fork.model import inner_system as IS
from fork.phase2 import identifiability as ID
from fork.phase2 import reference as R


@pytest.fixture(scope="module")
def problem():
    theta, _, ctx = ID.fixture_problem()
    jd = R.EPOCH_JD + np.linspace(-3650.0, 3650.0, 41)  # twenty years about the epoch
    return theta, jd, ctx


def test_matched_parameters_are_linear_and_gr_is_one(problem):
    theta, _, ctx = problem
    p = ID.matched_parameters(theta, ctx)
    assert p["gr"][0] == pytest.approx(1.0, abs=2e-3)  # 42.98″ is the 1PN rate itself
    assert 1e-4 < p["j2"][0] < 1e-3  # three orders above the Sun's 2.2e-7
    assert p["ring"][1] == ID.RING_RADIUS_AU and p["planet"][1] == ID.PLANET_A_AU
    q = ID.matched_parameters(theta, ctx, target=2 * ID.TARGET_ARCSEC_PER_CENTURY)
    for m in ID.MODELS:
        assert q[m][0] == pytest.approx(2 * p[m][0], rel=1e-12)


def test_synthetic_follows_the_forward_model_and_its_rng(problem):
    theta, jd, ctx = problem
    t = jd[18:23]
    clean = ID.synthetic("newtonian", theta, t, ctx)
    np.testing.assert_allclose(clean, IS.predict(theta, t, ctx), rtol=0, atol=1e-15)
    with pytest.raises(ValueError):
        ID.synthetic("newtonian", theta, t, ctx, sigma_arcsec=1.0)
    with pytest.raises(ValueError):
        ID.flags_for("other")
    a = ID.synthetic("gr", theta, t, ctx, (1.0,), 1.0, np.random.default_rng(3))
    b = ID.synthetic("gr", theta, t, ctx, (1.0,), 1.0, np.random.default_rng(3))
    assert np.array_equal(a, b)
    d = (a - ID.synthetic("gr", theta, t, ctx, (1.0,)) + np.pi) % (2 * np.pi) - np.pi
    assert 0 < np.abs(d).max() * ID.ARCSEC < 6.0


def test_compare_separates_the_models(problem):
    theta, jd, ctx = problem
    out = ID.compare(theta, jd, ctx)
    m = out["models"]
    target = ID.TARGET_ARCSEC_PER_CENTURY
    assert m["gr"]["varpi_rate"] == pytest.approx(target, rel=0.02)
    assert m["j2"]["varpi_rate"] == pytest.approx(target, rel=0.02)
    # a ring or planet in the ecliptic turns Mercury's inclined orbit about the ecliptic
    # pole as well, and the lowest-order planet rate ignores e = 0.2: within 20 %
    assert m["ring"]["varpi_rate"] == pytest.approx(target, rel=0.2)
    assert m["planet"]["varpi_rate"] == pytest.approx(target, rel=0.2)
    assert abs(m["gr"]["node_rate"]) < 0.01
    assert -6.0 < m["j2"]["node_rate"] < -2.0
    assert m["ring"]["node_rate"] < -30.0 and m["planet"]["node_rate"] < -30.0
    assert (
        m["planet"]["mean_longitude_periodic_rms"] > 10 * m["ring"]["mean_longitude_periodic_rms"]
    )
    assert set(out["pairs"]) == {"gr:j2", "gr:ring", "gr:planet", "j2:ring", "j2:planet",
                                 "ring:planet"}  # fmt: skip
    assert all(v > 0 for v in out["pairs"].values())
    with pytest.raises(ValueError):
        ID.compare(theta, jd[:2], ctx)
