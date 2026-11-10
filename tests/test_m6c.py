"""M6c runner: frozen inputs, the real data with the frozen error model, no Neptune in the fit
path (A6d), and the comparison targets."""

import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from fork.export import m6c

# Files opened while a test points the sink at a list (A6d rule 5, runtime). Audit hooks cannot
# be removed, so the hook is installed once and does nothing while the sink is None.
_AUDIT: dict[str, list | None] = {"sink": None}


def _audit(event: str, args: tuple) -> None:
    sink = _AUDIT["sink"]
    if sink is not None and event == "open" and isinstance(args[0], (str, bytes, os.PathLike)):
        sink.append(Path(os.fsdecode(args[0])).resolve())


sys.addaudithook(_audit)


def test_headline_jobs_and_worker_budget():
    jobs = m6c.all_jobs(("headline",))
    assert [j.id for j in jobs] == ["headline/c0", "headline/c1", "headline/c2", "headline/c3"]
    assert [j.nlive for j in jobs] == [2000, 1000, 1000, 1000]
    assert sum(j.workers for j in jobs) == 12
    assert all(m6c.Job.from_id(j.id) == j for j in jobs)
    assert [j.rng_seed() for j in jobs] == [[2026, 6, 1846, 0, k + 1] for k in range(4)]


def test_sensitivity_variants_change_one_thing_each():
    """S1–S3 (pre-registered): same cells, nlive, workers as the headline; distinct seeds; S1's
    cells partition its wider e prior; S2's error model is frozen; S3 keeps the headline's
    error model and Uranus box and swaps the likelihood."""
    from fork.inference.nested import cell_space
    from fork.inference.priors import free_a_space

    assert m6c.VARIANTS == ("headline", "s1_e05", "s2_sigma_data", "s3_t4")
    jobs = m6c.all_jobs()
    assert len(jobs) == 16 and len({tuple(j.rng_seed()) for j in jobs}) == 16
    for v in m6c.VARIANTS[1:]:
        vj = m6c.all_jobs((v,))
        assert [j.nlive for j in vj] == [2000, 1000, 1000, 1000]
        assert sum(j.workers for j in vj) == 12
        m6c.check_frozen(v)
    cfg = m6c.VARIANT_CONFIG
    head = cfg["headline"]

    def drop(d, *keys):
        return {k: v for k, v in d.items() if k not in keys}

    assert drop(cfg["s1_e05"], "e_max") == drop(head, "e_max")
    assert cfg["s1_e05"]["e_max"] == 0.5
    assert m6c.variant_cells("s1_e05")[3] == ((28.0, 60.0), (0.12, 0.5))
    w = sum(cell_space(free_a_space((0.0, 0.0), 0.5), c).slice_weight
            for c in m6c.variant_cells("s1_e05"))  # fmt: skip
    assert w == pytest.approx(1.0, abs=1e-12)
    assert drop(cfg["s2_sigma_data"], "error_model") == drop(head, "error_model")
    assert drop(cfg["s3_t4"], "likelihood", "nu") == drop(head, "likelihood")
    assert cfg["s3_t4"]["likelihood"] == "student_t" and cfg["s3_t4"]["nu"] == 4.0


def test_s2_error_model_is_data_only_with_the_headline_tau():
    import json

    from fork.model.observation import ErrorModel

    d = json.loads((m6c.ROOT / m6c.S2_ERROR_MODEL).read_text(encoding="utf-8"))
    assert "no JPL Uranus" in d["basis"] and d["source_point"]["run"] == "m6c headline"
    em, head = ErrorModel.load(m6c.ROOT / m6c.S2_ERROR_MODEL), ErrorModel.load()
    assert em.offset_prior_arcsec == head.offset_prior_arcsec
    for g, s in em.sigma_arcsec.items():
        assert s == pytest.approx(d["sigma_unrounded_arcsec"][g], abs=0.005)


def test_s3_dataset_is_student_t_on_the_headline_box():
    from fork.model.student_t import StudentTObservations

    *_, box_h = m6c.dataset()
    obs, *_, box = m6c.dataset(variant="s3_t4")
    assert isinstance(obs, StudentTObservations) and obs.nu == 4.0
    np.testing.assert_allclose(obs.t_scale_arcsec,
                               [m6c.t_scales_arcsec()[g] for g in m6c.t_scales_arcsec()],
                               rtol=0, atol=1e-12)  # fmt: skip
    np.testing.assert_allclose(box.centre, box_h.centre, rtol=0, atol=1e-12)
    np.testing.assert_allclose(box.half_width, box_h.half_width, rtol=0, atol=1e-12)


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


def test_fit_path_opens_nothing_under_tests_fixtures():
    """A6d rule 5, runtime half (the static half is tests/test_a6d.py): building the fit's data,
    context and Uranus box opens no file under tests/fixtures/. Its Horizons states come from
    FORK_CACHE_DIR, which the tests point at a temporary copy; outside pytest it is unset and
    the cache is data/cache/."""
    import fork.data.horizons as hz

    fixtures = (Path(__file__).parent / "fixtures").resolve()
    copy = Path(os.environ["FORK_CACHE_DIR"]).resolve()
    assert not copy.is_relative_to(fixtures)
    assert not hz.DEFAULT_CACHE_DIR.resolve().is_relative_to(fixtures.parent)
    opened: list[Path] = []
    _AUDIT["sink"] = opened
    try:
        m6c.dataset()
    finally:
        _AUDIT["sink"] = None
    assert not [p for p in opened if p.is_relative_to(fixtures)]
    horizons_reads = [p for p in opened if p.name.startswith("horizons_")]
    assert horizons_reads and all(p.parent == copy for p in horizons_reads)


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
