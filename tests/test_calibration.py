"""M6b error calibration (O − JPL per group): reproducible, robust, and only σ and τ leave it."""

import json

import numpy as np
import pytest

from fork.data import calibration as cal
from fork.model import observation as om


def test_committed_error_model_is_the_calibration():
    table, em = cal.calibrate(np.random.default_rng(20261106))
    d = json.loads(om.ERROR_MODEL_FILE.read_text(encoding="utf-8"))
    assert d["sigma_arcsec"] == em.sigma_arcsec
    assert d["offset_prior_arcsec"] == em.offset_prior_arcsec
    assert [r["n"] for r in d["table"]] == [c.n for c in table] == [19, 18, 87, 80, 107, 95, 66, 57]
    loaded = om.ErrorModel.load()
    assert loaded.sigma_arcsec == em.sigma_arcsec


def test_error_model_carries_no_offsets():
    """Group means are reported, never used: the error model has σ and τ only, and the
    offsets stay free with prior mean 0 (the synthetic draws are centred on 0 too)."""
    em = om.ErrorModel.load()
    assert set(vars(em)) == {"sigma_arcsec", "offset_prior_arcsec", "basis"}


def test_robust_sigma_resists_an_outlier():
    rng = np.random.default_rng(4)
    x = rng.normal(0.0, 2.0, 2000)
    assert cal.robust_sigma(x) == pytest.approx(2.0, rel=0.05)
    y = x.copy()
    y[0] = 1e4
    assert cal.robust_sigma(y) == pytest.approx(cal.robust_sigma(x), rel=0.01)
