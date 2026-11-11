"""Phase 1b step 3 (fork.export.p1b_adams): the pre-registered data set D and its test rows, the
held-a spaces, the jobs and their seeds, the frozen-input check, the targets against the printed
places and Neptune's, the reading rule's bands and rows. No fitting here."""

import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from fork.export import m6c
from fork.export import p1b_adams as pa

_AUDIT: dict[str, list | None] = {"sink": None}


def _audit(event: str, args: tuple) -> None:
    sink = _AUDIT["sink"]
    if sink is not None and event == "open" and isinstance(args[0], (str, bytes, os.PathLike)):
        sink.append(Path(os.fsdecode(args[0])).resolve())


sys.addaudithook(_audit)


def test_data_set_d_has_the_pre_registered_counts():
    from fork.model import observation as om

    fit, test = pa.adams_rows(om.read_rows())
    assert pa.row_counts(fit) == pa.EXPECTED_FIT
    assert pa.row_counts(test)["total"] == pa.EXPECTED_TEST
    assert not any(r["observer"] == "Flamsteed" and r["epoch"].startswith("1690") for r in fit)
    assert all(om.series_of(r) != "paris_1801_1827" for r in fit + test)
    assert all(int(r["epoch"][:4]) > 1840 and r["observatory"] == "Paris" for r in test)
    assert max(int(r["epoch"][:4]) for r in fit) == 1840
    assert len(fit) + len(test) == 529 - 2 - sum(
        1 for r in om.read_rows() if om.series_of(r) == "paris_1801_1827"
    )


def test_jobs_seeds_and_worker_budget():
    jobs = pa.all_jobs()
    assert pa.VARIANTS == ("h1_g", "h2_g", "f_g", "h1_t", "h2_t", "f_t")
    assert len(jobs) == 16 and len({tuple(j.rng_seed()) for j in jobs}) == 16
    assert all(pa.Job.from_id(j.id) == j for j in jobs)
    for v in pa.VARIANTS:
        vj = pa.all_jobs((v,))
        assert sum(j.workers for j in vj) == 12
        fit, err = pa.split_variant(v)
        if fit == "f":
            assert [j.nlive for j in vj] == [2000, 1000, 1000, 1000]
            assert pa.variant_cells(v) == m6c.variant_cells("headline")
        else:
            assert [j.nlive for j in vj] == [1000, 1000]
            assert pa.variant_cells(v) == ((None, (0.0, 0.12)), (None, (0.12, 0.3)))
        assert [j.rng_seed() for j in vj] == [
            [2026, 11, 3, pa.FITS.index(fit), pa.ERRORS.index(err), k + 1] for k in range(len(vj))
        ]
    seeds = {tuple(j.rng_seed()) for j in jobs}
    assert not seeds & {tuple(j.rng_seed()) for j in m6c.all_jobs()}


def test_frozen_inputs_are_checked(monkeypatch):
    pa.check_frozen()
    bad = dict(pa.FROZEN)
    bad["table77_csv"] = (bad["table77_csv"][0], "0" * 64)
    monkeypatch.setattr(pa, "FROZEN", bad)
    with pytest.raises(RuntimeError, match="pre-registered"):
        pa.check_frozen()


def test_dataset_and_held_space_never_read_neptune():
    from fork.inference.nested import PartitionedResult
    from fork.model.student_t import StudentTObservations

    sink: list[Path] = []
    _AUDIT["sink"] = sink
    try:
        obs, _ctx, _inc_raan, pq, box, test = pa.dataset("h1_t")
        spaces = [pa.space_for("h1_t", k, pq, box) for k in range(2)]
        free0 = pa.space_for("f_g", 0, pq, box)
    finally:
        _AUDIT["sink"] = None
    assert not any("neptune" in p.name for p in sink), [p.name for p in sink if "neptune" in p.name]
    assert obs.n_measurements == 278 and test.n_measurements == 47
    assert isinstance(obs, StudentTObservations) and obs.nu == 4.0
    assert np.allclose(obs.t_scale_arcsec, [m6c.t_scales_arcsec()[g] for g in obs_groups()])
    s = spaces[0]
    assert s.ndim == 10 and s.held[1] == pytest.approx(math.log(38.4))
    assert len(s.names) == 10
    assert s.unit_names == (
        "u_a",
        "u_k",
        "u_h",
        "u_lam0",
        "log10_m",
        "e",
        "varpi",
        "lam0",
        "cos_i",
        "node",
    )
    u = np.full(10, 0.5)
    theta = s.theta(s.prior_transform(u))
    assert theta.shape == (13,) and theta[7] == pytest.approx(38.4)
    assert 0.0 <= theta[8] <= 0.12 and s.log_prior(s.prior_transform(u)) > -math.inf
    assert spaces[1].e_bounds == (0.12, 0.3)
    assert abs(sum(sp.slice_weight for sp in spaces) - 1.0) < 1e-12
    assert free0.ndim == 11 and free0.a_au_range == m6c.variant_cells("headline")[0][0]
    assert PartitionedResult.__name__  # the weights above are what it checks


def obs_groups():
    from fork.model.observation import GROUPS

    return GROUPS


def test_targets_reproduce_the_printed_places_and_neptune():
    from fork.model.forward import uranus_inc_raan_at_epoch

    geos = pa.report_geometries()
    inc_raan = uranus_inc_raan_at_epoch(pa.A1_EPOCH_JD_TDB)
    places = pa.target_places(geos, inc_raan)
    # Neptune: the headline's Galle-night constant and step 1's 1 October 1846 mean longitude
    assert places["neptune"]["1846-09-23"]["geo_lon"] == pytest.approx(
        m6c.JPL_GALLE_LON_DEG, abs=0.01
    )
    assert places["neptune"]["1846-10-01"]["helio_true"] == pytest.approx(327.0, abs=0.05)
    # osculating mean longitude: about 1.2° ahead of the true one on that date (e ≈ 0.01)
    assert 327.5 < places["neptune"]["1846-10-01"]["helio_mean"] < 329.0
    # S5: Airy's No. 20 prints 325° 8′ mean longitude for 1 October 1846, the memoir 325° 7′ for
    # 6 October (the two printings differ by 2′ against the planet's 1′ of motion between the
    # dates); carried from the memoir's value, 325° 6′ on 1 October: within 3′ of Airy's
    assert places["S5"]["1846-10-01"]["helio_mean"] == pytest.approx(325.0 + 8.0 / 60.0, abs=0.05)
    # S4: the sheet's geocentric "320° 30′ nearly ... diminishing about 1′ daily" for "about the
    # end of September 1845"; from its elements on 30 September, 320° 48′: the sheet's reduction
    # and date are loose to a few tenths of a degree, so both values are reported (plan § 5.1)
    assert places["S4"]["1845-09-30"]["geo_lon"] == pytest.approx(
        pa.S4_PRINTED_GEOCENTRIC_DEG, abs=0.35
    )
    # LV-June: 325° heliocentric at 1847.0 by construction, so 0.38° less on 1 October 1846
    assert places["LV-June"]["1846-10-01"]["helio_true"] == pytest.approx(325.0 - 0.38, abs=0.03)
    # LV-Aug: the published true longitude at 1847.0 (M6c's check) carried back a quarter
    lv = pa.targets()["LV-Aug"]["orbit"]
    assert lv.true_lon_ecl_of_date_deg == pytest.approx(326.0 + 32.0 / 60.0, abs=0.02)
    # RA/Dec of date on a Challis night: step 2's neptune_on_challis_dates.csv has Neptune at
    # 21h 57m 24s, −13° 02′ (mean equinox of date) on 12 August 1846
    n = places["neptune"]["1846-08-12"]
    assert 21.94 < n["ra_h"] < 21.97 and -13.1 < n["dec_deg"] < -13.0


def test_reading_rule_bands_and_rows():
    assert pa.band_of(-2.5) == "disfavoured"
    assert pa.band_of(-1.0) == "weakly disfavoured; reported, no statement"
    assert pa.band_of(-0.2) == "not distinguished"
    assert pa.status_of({"in_68": True, "in_95": True}) == "inside"
    assert pa.status_of({"in_68": False, "in_95": True}) == "marginal"
    assert pa.status_of({"in_68": False, "in_95": False}) == "outside"

    def s1(reading):
        return {"F": {"status": "x", "reading": reading}}

    def s2(s5, s6):
        return {"S5": {"true_lon": s5, "mean_lon": s5}, "S6": {"true_lon": s6, "mean_lon": s6}}

    def s3(b1, b2):
        return {"H1": {"band": b1}, "H2": {"band": b2}}

    dis = "disfavoured"
    assert pa.interpret(s1("partly"), s2("inside", "inside"), s3(dis, dis))["row"] == 4
    assert pa.interpret(s1("yes"), s2("inside", "inside"), s3(dis, dis))["row"] == 1
    nd = "not distinguished"
    assert pa.interpret(s1("yes"), s2("inside", "inside"), s3(nd, nd))["row"] == 2
    assert pa.interpret(s1("yes"), s2("outside", "outside"), s3(nd, nd))["row"] == 3
    weak = "weakly disfavoured; reported, no statement"
    assert pa.interpret(s1("yes"), s2("inside", "inside"), s3(weak, nd))["row"] == "amendment 3"
    assert "disagree" in pa.interpret(s1("yes"), s2("inside", "outside"), s3(nd, nd))["row"]
    assert "marginal" in pa.interpret(s1("yes"), s2("marginal", "marginal"), s3(nd, nd))["row"]
    assert "marginal" in pa.interpret(s1("yes"), s2("marginal", "outside"), s3(nd, nd))["row"]


def test_report_epochs_are_the_fetch_plan_s():
    from fork.data.frames import jd_tdb_from_iso
    from fork.data.horizons import P1B_REPORT_EPOCHS_ISO_TDB, default_fetch_plan

    plan = default_fetch_plan()
    for body in ("earth", "neptune"):
        for iso in pa.REPORT_EPOCHS.values():
            assert np.min(np.abs(plan[body] - jd_tdb_from_iso(iso))) < 1e-9
    galle = {"1846-09-23T23:00:00"}
    assert set(P1B_REPORT_EPOCHS_ISO_TDB) == set(pa.REPORT_EPOCHS.values()) - galle
