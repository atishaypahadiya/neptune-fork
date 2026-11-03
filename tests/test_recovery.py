"""A3: recovering a planted planet from synthetic data (PHASE1_SPEC M5, as amended after M5a).

A3a (bias) re-runs the noise-free point fit. A3b (calibration) checks the committed M5b posterior
summary, ``site/data/m5b_posteriors.json``: the chains behind it take hours and run locally
(``python -m fork.export.m5b``), never in CI.
"""

import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from fork.data.horizons import states_about_inner_barycenter_ecl_j2000
from fork.data.synthetic import LongitudeObservations, make_synthetic_observations
from fork.inference.fit import fit_perturber_multistart, wrap_rad
from fork.model import bodies
from fork.model.forward import ForwardContext, perturber_theta_from_state, uranus_inc_raan_at_epoch
from fork.model.hidden_mass import perturber_geocentric_longitude, true_neptune_geocentric_longitude

EPOCH = 2371628.5
N_WORKERS = min(12, os.cpu_count() or 1)
POSTERIORS = Path(__file__).resolve().parents[1] / "site" / "data" / "m5b_posteriors.json"

pytestmark = pytest.mark.usefixtures("offline_horizons")
# TEMPORARY: the M5b chains were still running when A3 was rewritten. This skip goes away in the
# commit that adds site/data/m5b_posteriors.json; after that a missing file must fail, not skip.
needs_posteriors = pytest.mark.skipif(
    not POSTERIORS.exists(), reason="m5b_posteriors.json not committed yet (chains running)"
)


@pytest.mark.slow
def test_a3a_noise_free_recovery_within_half_degree(offline_horizons):
    """A3a. Noise-free Horizons truth (real 3-D Neptune in the data, planar perturber in the
    model), free a, e ≤ 0.3, 24 λ₀ starts: the best fit's 1846-09-23 geocentric longitude is within
    0.5° of true Neptune's. Measured −0.07° (m 17.14 M⊕ vs 17.15, a 30.24 AU vs osculating 30.22,
    e 0.004). ~60 s: the slowest noise-free LM run sets the wall time."""
    ctx = ForwardContext.from_horizons(EPOCH)
    obs, truth = make_synthetic_observations(np.random.default_rng(0))
    clean = LongitudeObservations(obs.jd_tdb, truth, obs.sigma_rad, obs.source)
    best = fit_perturber_multistart(clean, ctx, uranus_inc_raan_at_epoch(EPOCH),
                                    n_workers=N_WORKERS).best  # fmt: skip
    g, nep = perturber_geocentric_longitude(best.theta, ctx), true_neptune_geocentric_longitude()
    dlon_deg = math.degrees(wrap_rad(g.lon_ecl_of_date_rad - nep.lon_ecl_of_date_rad))
    assert abs(dlon_deg) < 0.5
    th_nep = perturber_theta_from_state(
        states_about_inner_barycenter_ecl_j2000("neptune", [EPOCH])[0], bodies.NEPTUNE.mass_msun
    )
    assert abs(best.mass_earth - 10 ** th_nep[0]) < 0.5
    assert abs(best.theta[7] - th_nep[1]) < 0.3 and best.theta[8] < 0.02


@needs_posteriors
def test_a3b_calibration_of_committed_posteriors():
    """A3b. Across ≥ 10 noise seeds, true Neptune's 1846-09-23 longitude lies inside the 68 %
    posterior interval (free a, e ∈ [0, 0.3]) in 50–85 % of seeds. The summary is recomputed
    here from the per-seed records, not trusted from the file's own flag."""
    d = json.loads(POSTERIORS.read_text())
    noisy = [v for v in d["datasets"].values() if v["noise"]]
    assert len(noisy) >= 10
    for v in noisy:
        g = v["runs"]["free_a_e0.3"]["galle_lon_ecl_of_date_deg"]
        rel = [
            (x - g["median"] + 180.0) % 360.0 - 180.0 for x in (g["lo68"], g["truth"], g["hi68"])
        ]
        assert g["truth_in_68"] == (rel[0] <= rel[1] <= rel[2])
    hits = [v["runs"]["free_a_e0.3"]["galle_lon_ecl_of_date_deg"]["truth_in_68"] for v in noisy]
    frac = float(np.mean(hits))
    assert d["a3b_calibration"]["free_a_e0.3"]["fraction"] == pytest.approx(frac, abs=1e-3)
    assert 0.50 <= frac <= 0.85


@needs_posteriors
def test_committed_posteriors_converged():
    """Every committed chain ran ≥ 50 τ (spec) and the summary has all three runs per dataset."""
    d = json.loads(POSTERIORS.read_text())
    assert "noise_free" in d["datasets"]
    for v in d["datasets"].values():
        assert set(v["runs"]) == {"free_a_e0.3", "free_a_e0.1", "le_verrier_shape"}
        for r in v["runs"].values():
            assert r["converged_50tau"] and r["n_steps"] >= 50 * max(r["tau"])
            assert r["n_walkers"] >= 32
