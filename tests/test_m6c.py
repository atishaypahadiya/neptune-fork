"""M6c runner: frozen inputs, the real data with the frozen error model, no Neptune in the fit
path (A6d), and the comparison targets."""

import math

import numpy as np
import pytest

from fork.export import m6c


def test_headline_jobs_and_worker_budget():
    jobs = m6c.all_jobs()
    assert [j.id for j in jobs] == ["headline/c0", "headline/c1", "headline/c2", "headline/c3"]
    assert [j.nlive for j in jobs] == [2000, 1000, 1000, 1000]
    assert sum(j.workers for j in jobs) == 12
    assert all(m6c.Job.from_id(j.id) == j for j in jobs)
    assert len({tuple(j.rng_seed()) for j in jobs}) == 4


def test_frozen_inputs_are_checked(monkeypatch):
    m6c.check_frozen()  # the committed files are the pre-registered ones
    bad = dict(m6c.FROZEN, error_model=(m6c.FROZEN["error_model"][0], "0" * 64))
    monkeypatch.setattr(m6c, "FROZEN", bad)
    with pytest.raises(RuntimeError, match="pre-registered"):
        m6c.check_frozen()


def test_dataset_is_table77_with_the_frozen_error_model_and_never_reads_neptune(monkeypatch):
    """A6d for the fit path: building the data, the context and the Uranus box reads no
    Neptune state (every Horizons state goes through ``state_vectors_ecl_j2000``)."""
    import fork.data.horizons as hz
    from fork.model.observation import GROUPS, ErrorModel

    real = hz.state_vectors_ecl_j2000

    def guarded(body, *args, **kwargs):
        if "neptune" in str(body).lower():
            raise AssertionError("the M6c fit path read Neptune's ephemeris")
        return real(body, *args, **kwargs)

    monkeypatch.setattr(hz, "state_vectors_ecl_j2000", guarded)
    obs, *_, box = m6c.dataset()
    em = ErrorModel.load()
    assert obs.n_measurements == 529 and np.all(np.isfinite(obs.value_rad))
    np.testing.assert_allclose(obs.sigma_arcsec, [em.sigma_arcsec[g] for g in GROUPS], rtol=0,
                               atol=1e-12)  # fmt: skip
    np.testing.assert_allclose(obs.offset_prior_arcsec,
                               [5.0 if g.startswith("ancient") else 2.0 for g in GROUPS],
                               rtol=0, atol=1e-12)  # fmt: skip
    assert box.chi2_per_dof > 1.0  # the Uranus-only fit leaves the perturber's signal


def test_jpl_galle_target_matches_horizons():
    from fork.model.hidden_mass import true_neptune_geocentric_longitude

    g = math.degrees(true_neptune_geocentric_longitude().lon_ecl_of_date_rad)
    assert g == pytest.approx(m6c.JPL_GALLE_LON_DEG, abs=5e-4)


@pytest.mark.network
def test_jpl_1847_heliocentric_target_matches_horizons():
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.data.horizons import states_about_inner_barycenter_ecl_j2000 as about
    from fork.model.hidden_mass import JD_1847_JAN1_TDB

    r = precession_matrix_ecl_j2000_to_ecl_of_date(JD_1847_JAN1_TDB) @ about(
        "neptune", [JD_1847_JAN1_TDB])[0, :3]  # fmt: skip
    lon = math.degrees(math.atan2(r[1], r[0])) % 360.0
    assert lon == pytest.approx(m6c.JPL_HELIO_1847_LON_DEG, abs=5e-4)


def test_le_verrier_published_orbit_reproduces_his_longitude_and_misses_by_54_arcmin():
    """The Kepler orbit from his published elements gives back his 1847.0 true longitude
    (326°32′) to < 0.02°; on Galle's night it lies 54′ short of JPL's Neptune, in line with the
    ~52′ miss the spec quotes (the remaining 2′ is within the elements' rounding)."""
    from fork.data.frames import precession_matrix_ecl_j2000_to_ecl_of_date
    from fork.model.hidden_mass import JD_1847_JAN1_TDB, LE_VERRIER_1846

    pos, _ = m6c.le_verrier_state_ecl_j2000(JD_1847_JAN1_TDB)
    r = precession_matrix_ecl_j2000_to_ecl_of_date(JD_1847_JAN1_TDB) @ pos
    lon = math.degrees(math.atan2(r[1], r[0])) % 360.0
    assert lon == pytest.approx(LE_VERRIER_1846.true_lon_ecl_of_date_deg, abs=0.02)
    g = m6c.le_verrier_galle_lon_deg()
    assert g == pytest.approx(324.973, abs=1e-3)
    assert (g - m6c.JPL_GALLE_LON_DEG) * 60.0 == pytest.approx(-54.0, abs=0.1)


def test_lon_summary_intervals_offsets_and_wrap():
    rng = np.random.default_rng(1)
    lon = (rng.normal(0.0, 2.0, 200_000) + 360.0) % 360.0  # straddles 0°/360°
    s = m6c.lon_summary(lon, {"near": 1.0, "far": 5.0})
    assert s["width68"] == pytest.approx(4.0, abs=0.05)
    assert s["width95"] == pytest.approx(7.84, abs=0.1)
    near, far = s["vs"]["near"], s["vs"]["far"]
    assert near["in_68"] and near["in_95"] and near["offset_sigma"] == pytest.approx(-0.5, abs=0.02)
    assert not far["in_68"] and not far["in_95"]
    assert far["median_minus_target_deg"] == pytest.approx(-5.0, abs=0.02)


def test_inclination_prior_quantiles():
    lo, med, hi = m6c.inc_prior_quantiles_deg()
    assert med == pytest.approx(7.07, abs=0.01) and lo < med < hi < 10.0
