"""The synthetic identifiability study for Mercury (docs/phase2/identifiability_study.md):
what is read from the contact table, the data generator, the likelihood algebra and the
readings on made-up inputs. No integration over two centuries and no sampling here."""

import builtins
import csv
import math
from pathlib import Path

import numpy as np
import pytest

from fork.phase2 import identifiability_study as St


def test_only_the_epoch_and_weight_columns_are_read(monkeypatch):
    """V5, first half: a row is cut down to the allowed columns before anything is used, so
    the k terms, the second value of the braced row and the flags cannot enter."""
    seen: set[str] = set()

    class Spy(dict):
        def __getitem__(self, key):
            seen.add(key)
            return super().__getitem__(key)

    real = csv.DictReader

    def spying(fh, *a, **kw):
        return (Spy(row) for row in real(fh, *a, **kw))

    monkeypatch.setattr(St.csv, "DictReader", spying)
    ep = St.contact_epochs()
    assert seen == set(St.COLUMNS_READ)
    assert not seen & {"k_sign", "k_coeff", "second_value", "flag"}
    assert len(ep["jd"]) == 58 and set(ep["month"]) == {5, 11}
    assert np.all((ep["weight"] > 0) & (ep["weight"] <= 6))


def test_epochs_are_on_the_astronomical_day():
    """1677 Nov 6, 21h 34m 1s astronomical is the morning of 7 November, civil."""
    assert St.julian_day_number(2000, 1, 1) == 2451545
    ep = St.contact_epochs()
    jd0 = St.julian_day_number(1677, 11, 6) + (21 + 34 / 60 + 1 / 3600) / 24
    assert ep["jd"][0] == pytest.approx(jd0, abs=1e-9)
    assert ep["jd"].min() > St.julian_day_number(1677, 1, 1)
    assert ep["jd"].max() < St.julian_day_number(1882, 1, 1)


def test_synthetic_data_use_nothing_but_their_arguments(monkeypatch):
    """V5, second half: no file is opened while a data set is made; it is the truth's
    longitudes plus seeded noise of the given size, and nothing else."""

    def refuse(*a, **kw):
        raise AssertionError("a file was opened while making synthetic data")

    y = np.linspace(-1e-5, 1e-5, 58)
    sigma = np.full(58, 2e-6)
    monkeypatch.setattr(builtins, "open", refuse)
    monkeypatch.setattr(Path, "open", refuse)
    monkeypatch.setattr(Path, "read_text", refuse)
    free = St.synthetic_data(y, sigma, 2, St.NOISE_FREE)
    d1 = St.synthetic_data(y, sigma, 2, 3)
    d2 = St.synthetic_data(y, sigma, 2, 3)
    other = St.synthetic_data(y, sigma, 3, 3)
    wide = St.synthetic_data(y, 2.0 * sigma, 2, 3)
    np.testing.assert_allclose(free, y, rtol=0, atol=0)
    assert np.array_equal(d1, d2)  # same seed, same machine
    assert not np.allclose(d1, other)
    np.testing.assert_allclose(wide - y, 2.0 * (d1 - y), rtol=1e-12)
    assert abs(np.std((d1 - y) / sigma) - 1.0) < 0.35


def test_no_residual_table_is_named_in_the_module():
    src = Path(St.__file__).read_text(encoding="utf-8")
    for word in ("453", "454", "residual_table", "observed_minus", "o_minus_c"):
        assert word not in src


def _prep(rng: np.random.Generator, n: int = 58) -> dict:
    t = np.linspace(-1.0, 1.0, n)
    jac = np.c_[np.ones(n), t, np.sin(3 * t), np.cos(3 * t), t**2, np.sin(7 * t)]
    ring = np.array([(1.0 + r) * t**3 + 0.2 * r * np.cos(5 * t) for r in St.RING_GRID])
    lam = 2.0 * math.pi * np.arange(St.N_PHASE) / St.N_PHASE
    planet = np.array([[t**3 + 0.3 * np.sin(11 * t + ph) * (1 + ia) for ph in lam]
                       for ia in range(3)])  # fmt: skip
    return {"jac": jac, "G": t**3, "J": t**3 + 0.1 * t**4, "R": ring, "P": planet}


def test_marginal_likelihood_matches_a_direct_fit():
    """rᵀ P r is the χ² minimised over the six linear element changes."""
    rng = np.random.default_rng(4)
    prep = _prep(rng)
    sigma = rng.uniform(0.5, 2.0, 58)
    y = rng.normal(0.0, 1.0, 58) + 3.0 * prep["G"]
    f = St.Fitter(prep, sigma, y)
    for s in (0.0, 1.0, 2.5):
        r = y - s * prep["G"]
        a = prep["jac"] / sigma[:, None]
        coef, *_ = np.linalg.lstsq(a, r / sigma, rcond=None)
        direct = float(np.sum((r / sigma - a @ coef) ** 2))
        assert float(f.chi2(prep["G"], s)) == pytest.approx(direct, rel=1e-9, abs=1e-9)
        np.testing.assert_allclose(f.best_elements(prep["G"], s), coef, rtol=1e-8, atol=1e-10)
    assert f.lnl_fixed("N") == pytest.approx(-0.5 * float(f.chi2(prep["G"], 0.0)))
    assert f.lnl_fixed("G") == pytest.approx(-0.5 * float(f.chi2(prep["G"], 1.0)))


def test_planet_shape_interpolates_between_the_computed_longitudes():
    rng = np.random.default_rng(5)
    prep = _prep(rng)
    f = St.Fitter(prep, np.ones(58), np.zeros(58))
    lam = 2.0 * math.pi * np.arange(St.N_PHASE) / St.N_PHASE
    np.testing.assert_allclose(f.planet_shape(1, lam), prep["P"][1], atol=1e-12)
    t = np.linspace(-1.0, 1.0, 58)
    mid = 0.4321
    np.testing.assert_allclose(f.planet_shape(2, np.array([mid]))[0],
                               t**3 + 0.9 * np.sin(11 * t + mid), atol=1e-10)  # fmt: skip
    np.testing.assert_allclose(f.ring(0.2), prep["R"][5], atol=1e-12)


def test_quadrature_against_a_known_integral():
    """Model J on data with no signal: χ² is quadratic in the strength, so ln Z is a
    Gaussian integral over the strength with the log-uniform prior's Jacobian."""
    rng = np.random.default_rng(6)
    prep = _prep(rng)
    f = St.Fitter(prep, np.ones(58), rng.normal(0, 1, 58))
    x = np.linspace(*St.LOG_J2, 200001)
    exact = math.log(np.trapezoid(np.exp(f.lnl("J", x[:, None])), x) / (x[-1] - x[0]))
    assert St.ln_z_quadrature(f, "J") == pytest.approx(exact, abs=1e-3)
    z2 = St.ln_z_quadrature(f, "R", n1=801, n2=81)
    z3 = St.ln_z_quadrature(f, "R", n1=1601, n2=161)
    assert z2 == pytest.approx(z3, abs=0.01)


def test_prior_boxes_and_the_width_scan():
    np.testing.assert_allclose(St.bounds("J"), [[-7.0, -2.0]])
    np.testing.assert_allclose(St.bounds("R"), [[-9.0, -5.0], [0.10, 0.30]])
    np.testing.assert_allclose(St.bounds("P"), [[-9.0, -5.0], [0.0, 2.0 * math.pi]])
    np.testing.assert_allclose(St.bounds("J", 2.0)[0], [-4.5, -2.5])  # the amended centre
    np.testing.assert_allclose(St.bounds("R", 4.0)[0], [-9.0, -5.0])
    assert St.planet_average([0.0, 0.0, 0.0]) == pytest.approx(0.0)
    assert St.planet_average([math.log(3.0), -np.inf, -np.inf]) == pytest.approx(0.0)


def _fit(seed: int, **z: float) -> dict:
    return {"seed": seed, "ln_z": {"N": -500.0, "G": 0.0, "J": 0.0, "R": 0.0, "P": 0.0} | z}


def test_r1_verdicts():
    rec = [_fit(s, R=6.0) for s in range(9)]
    assert St.read_r1(rec, "b")["verdict"] == "recovered"
    weak = [_fit(s, R=1.5) for s in range(9)]
    assert St.read_r1(weak, "b")["verdict"] == "preferred, not decided"
    lost = [_fit(s, G=4.0) for s in range(9)]
    r = St.read_r1(lost, "b")
    assert r["verdict"] == "not recovered; most often G" and r["wins_of_8"] == 0
    six = [_fit(0, R=6.0)] + [_fit(s, R=6.0) for s in range(1, 7)] + [_fit(7, G=2.0), _fit(8, G=2)]
    assert St.read_r1(six, "b")["verdict"] == "recovered"
    five = [_fit(0, R=6.0)] + [_fit(s, R=6.0) for s in range(1, 6)] + [
        _fit(s, R=1.0) for s in (6, 7, 8)]  # fmt: skip
    assert St.read_r1(five, "b")["verdict"] == "preferred, not decided"
    assert [St.grade(x) for x in (0.5, -2.0, 3.0, -7.0)] == ["none", "weak", "strong", "decisive"]


def test_the_plan_names_these_settings():
    plan = (St.ROOT / "docs" / "phase2" / "identifiability_study.md").read_text(encoding="utf-8")
    frozen = plan.split("\n---\n")[0]
    for text in ("nlive 1000, dlogz 0.01", "`[2026, 21,\n1882, truth, seed]`", "10⁻⁹ to 10⁻⁵ M☉",
                 "0.15, 0.20, 0.25 AU", "11 ring radii", "48 planet longitudes", "3.0 %",
                 "1, 2, 4, 6 and 8 decades"):  # fmt: skip
        assert text in frozen, text
    assert "centred on 10⁻³·⁵" in plan.split("\n---\n")[1]
    assert (St.NLIVE, St.DLOGZ, St.N_PHASE, len(St.RING_GRID)) == (1000, 0.01, 48, 11)
    assert St.LEVELS == (1.0, 0.5, 2.0) and St.SEEDS == tuple(range(1, 9))
