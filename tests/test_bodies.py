import pytest

from fork.model import bodies


def test_gaussian_constant_and_g():
    assert bodies.K_GAUSS == 0.01720209895
    assert bodies.G == pytest.approx(0.01720209895**2, rel=0, abs=0)


@pytest.mark.parametrize(
    ("name", "ratio"),
    [("jupiter", 1047.348644), ("saturn", 3497.9018), ("uranus", 22902.98), ("neptune", 19412.26)],
)
def test_spec_mass_ratios(name, ratio):
    assert bodies.BODIES[name].mass_msun == 1.0 / ratio


def test_inner_planets_folded_into_central_mass():
    # Mercury 1.66e-7 + Venus 2.448e-6 + Earth–Moon 3.040e-6 + Mars 3.23e-7
    assert bodies.M_INNER_MSUN == pytest.approx(5.977e-6, rel=1e-3)
    assert bodies.M_CENTRAL_MSUN == 1.0 + bodies.M_INNER_MSUN
