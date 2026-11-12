"""Phase 2 forces (1PN, J2, ring, inner planet): bit-identity with the Phase 1 map when off,
the fields against independent calculations, and V1–V5 of docs/phase2/validation.md on short
spans. The long versions are run by ``python -m fork.phase2.validation run``."""

import math

import numpy as np
import pytest

from fork.integrators import wh_extra as W
from fork.integrators import wisdom_holman as WH
from fork.model import inner_system as IS
from fork.model.bodies import G, masses_msun
from fork.phase2 import reference as R
from fork.phase2 import secular as S
from fork.phase2 import validation as V

C = S.ARCSEC_PER_CENTURY


@pytest.fixture(scope="module")
def ref():
    return R.load(R.FIXTURE)


def test_off_is_the_phase_1_map_bit_for_bit(ref):
    """Documented invariant: with every extra force off, the step is the Phase 1 step."""
    m = np.concatenate([[1.0], masses_msun(R.PLANETS)])
    x0, v0 = ref["epoch_states"][:, :3].copy(), ref["epoch_states"][:, 3:].copy()
    t = np.array([-3000.0, -17.3, 0.0, 12.0, 2500.5])
    a = WH.integrate_states(x0, v0, m, t, 2.0)
    b = W.integrate_states(x0, v0, m, t, 2.0, W.make_extra())
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def test_elliptic_integrals():
    from scipy.special import ellipe, ellipk

    for m in (0.0, 1e-9, 0.1, 0.5, 0.9, 0.999):
        k, e = W.ellip_ke(m)
        assert k == pytest.approx(ellipk(m), rel=1e-13)
        assert e == pytest.approx(ellipe(m), rel=1e-13)


def test_ring_field_is_the_sum_over_the_ring_and_the_gradient_of_its_potential():
    normal = np.array([0.1, 0.2, 1.0])
    ex = W.make_extra(ring_mass=1e-7, ring_radius_au=0.2, ring_normal=tuple(normal))
    n = normal / np.linalg.norm(normal)
    u = np.cross(n, [1.0, 0.0, 0.0])
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    ph = 2 * np.pi * np.arange(20000) / 20000
    pts = 0.2 * (np.cos(ph)[:, None] * u + np.sin(ph)[:, None] * w)
    acc = np.empty(3)
    for p in (np.array([0.3, 0.1, 0.05]), np.array([0.05, -0.02, 0.3]), np.array([0.4, 0, 0])):
        W.position_accel(p, G, ex, acc)
        d = pts - p
        direct = (G * 1e-7 / 20000 * d / np.linalg.norm(d, axis=1)[:, None] ** 3).sum(axis=0)
        np.testing.assert_allclose(acc, direct, rtol=1e-9)
        fd = [-(W.position_potential(p + h, G, ex) - W.position_potential(p - h, G, ex)) / 2e-6
              for h in np.eye(3) * 1e-6]  # fmt: skip
        np.testing.assert_allclose(acc, fd, rtol=1e-6)


def test_j2_field_is_the_gradient_of_its_potential():
    ex = W.make_extra(j2=2.2e-7, spin_axis=(0.2, -0.1, 1.0))
    p = np.array([0.3, 0.1, 0.05])
    acc = np.empty(3)
    W.position_accel(p, G, ex, acc)
    fd = [-(W.position_potential(p + h, G, ex) - W.position_potential(p - h, G, ex)) / 2e-6
          for h in np.eye(3) * 1e-6]  # fmt: skip
    np.testing.assert_allclose(acc, fd, rtol=1e-6)


def test_the_analytic_gr_rate_is_42_98():
    assert S.gr_rate(V.A_MERCURY, V.E_MERCURY, G) * C == pytest.approx(42.98, abs=0.005)


def test_v1_gr_advance_short():
    r = V.v1(years=100.0)
    assert r["pass"], r["rows"]
    assert abs(r["rows"][0]["measured"] - 42.98) < 0.1


def test_v2_j2_short():
    r = V.v2(years=100.0)
    assert r["pass"], r["rows"]


def test_v3_ring_short():
    r = V.v3(years=100.0, point_planet=False)
    assert r["pass"], r["rows"]


@pytest.mark.slow
def test_v3_point_planet():
    r = V.v3(years=300.0, point_planet=True)
    assert r["pass"], r["rows"]


def test_v4_against_the_ephemeris_fixture(ref):
    r = V.v4(ref)
    off, on = r["cases"]["newtonian"], r["cases"]["with_1pn"]
    assert off["perihelion_residual_slope_arcsec_per_century"] == pytest.approx(-43.0, abs=0.3)
    assert abs(on["perihelion_residual_slope_arcsec_per_century"]) < 0.1
    assert on["max_abs_dlon"] < 0.2  # arcsec, over 1800–2000
    assert off["max_abs_dlon"] > 20.0


def test_v5_conservation_short():
    r = V.v5_conservation(years=50.0)
    for row in r["rows"]:
        assert row["max_rel_energy"] < 1e-9, row
        assert row["max_rel_angular_momentum"] < 1e-9, row


def test_forward_model_reproduces_the_ephemeris_and_its_blocks_switch(ref):
    ctx = IS.InnerSystemContext(R.EPOCH_JD, ref["epoch_states"][1:])
    th = IS.mercury_theta_from_state(ref["epoch_states"][0])
    jd = ref["jd"][::6]
    truth = np.arctan2(ref["mercury"][::6, 1], ref["mercury"][::6, 0])
    lon = IS.predict(np.append(th, 1.0), jd, ctx, IS.Flags(gr=True))
    d = (lon - truth + math.pi) % (2 * math.pi) - math.pi
    assert np.abs(np.degrees(d)).max() * 3600 < 0.2
    newton = IS.predict(th, jd, ctx)
    zero = IS.predict(np.concatenate([th, [0.0, 0.0, 0.0, 0.2]]), jd, ctx,
                      IS.Flags(gr=True, j2=True, ring=True))  # fmt: skip
    np.testing.assert_allclose(zero, newton, atol=1e-12)
    planet = IS.predict(np.concatenate([th, [1e-7, 0.2, 0.3]]), jd, ctx, IS.Flags(planet=True))
    assert np.abs(planet - newton).max() > 1e-7
    with pytest.raises(ValueError):
        IS.predict(th, jd, ctx, IS.Flags(gr=True))
