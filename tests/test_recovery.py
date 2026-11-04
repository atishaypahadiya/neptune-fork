"""A3: recovering a planted planet from synthetic data (PHASE1_SPEC M5, as amended after M5a).

A3a (bias) re-runs the noise-free point fit. A3b (calibration) checks the committed M5b posterior
summary, ``site/data/m5b_posteriors.json``: the dynesty runs behind it take hours and run
locally (``python -m fork.export.m5b``), never in CI.
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
    posterior interval (free a, e ∈ [0, 0.3]; dynesty posteriors, all modes) in 50–85 % of
    seeds. The summary is recomputed here from the per-seed records, not trusted from the file's
    own flag."""
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
def test_committed_posteriors_complete():
    """Every noise seed has the posteriors it is assigned (``m5b.runs_for``: e ≤ 0.3 on seeds
    0–9, e ≤ 0.1 and Le Verrier shape on 0–4; Le Verrier shape beyond seed 0 only if the
    single-mode check passed). Each partitioned run has every cell of ``nested.cells_for`` at
    its nlive, ln Z errors carry the measured inflation, and the M5a profile interval is
    attached wherever the dataset was profiled."""
    from fork.export.m5b import E_MAX, NLIVE_CELL, NLIVE_DEFAULT, runs_for
    from fork.inference.nested import ERR_INFLATION, cells_for

    d = json.loads(POSTERIORS.read_text())
    assert d["sampler"]["err_inflation"] == ERR_INFLATION
    lv_ok = d["le_verrier_single_mode_check"]["single_mode"]
    assert sorted(v["seed"] for v in d["datasets"].values() if v["noise"]) == list(range(10))
    assert set(d["datasets"]["noise_free"]["runs"]) == {"free_a_e0.3"}
    assert d["noise_free_check"]["status"] == "done"
    for v in d["datasets"].values():
        want = {r for r in runs_for(v["seed"])
                if r != "le_verrier_shape" or lv_ok or v["seed"] == 0}  # fmt: skip
        assert set(v["runs"]) == want
        assert "profile_m5a" in v
        for name, r in v["runs"].items():
            ev = r["evidence"]
            assert ev["logzerr"] == pytest.approx(ERR_INFLATION * ev["logzerr_quoted"], abs=2e-4)
            if name == "le_verrier_shape":
                assert not ev["partitioned"]
                continue
            cells = cells_for(E_MAX[name])
            assert len(ev["cells"]) == len(cells)
            assert sum(c["share"] for c in ev["cells"]) == pytest.approx(1.0, abs=1e-4)
            for k, c in enumerate(ev["cells"]):
                assert c["nlive"] == NLIVE_CELL.get((name, k), NLIVE_DEFAULT)
