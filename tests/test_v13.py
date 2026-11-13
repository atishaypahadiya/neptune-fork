"""Building blocks for the v1.3 draft plans (fork.export.v13), on synthetic data only."""

import math

import numpy as np
import pytest

from fork.export import v13
from fork.model import observation as om


@pytest.fixture(scope="module")
def synth():
    return v13.synthetic_dataset(1, "leverrier")


def test_data_choices_select_the_rows_of_the_spec():
    rows = om.read_rows()
    n = {c: len(v13.select_rows(rows, c)) for c in v13.DATA_CHOICES}
    assert n["leverrier"] == len(rows) == 529
    assert n["bouvard"] == 492 and all(
        r["era"] != "ancient" for r in v13.select_rows(rows, "bouvard")
    )
    dropped = [r for r in rows if r not in v13.select_rows(rows, "adams")]
    assert dropped and all(
        r["observer"] == "Flamsteed" and r["epoch"][:4] == "1690" for r in dropped
    )
    with pytest.raises(ValueError):
        v13.select_rows(rows, "other")


def test_empty_groups_are_refused_unless_asked_for():
    modern = v13.select_rows(om.read_rows(), "bouvard")
    with pytest.raises(ValueError):
        om.structure_from_rows(modern, om.ErrorModel.load())
    st = om.structure_from_rows(modern, om.ErrorModel.load(), allow_empty_groups=True)
    assert st.n_measurements == 492 and st.group_counts()[:2].tolist() == [0, 0]


def test_synthetic_layouts_give_a_finite_likelihood(synth):
    for choice in ("bouvard", "adams"):
        obs, ctx, pq, box = v13.synthetic_dataset(1, choice)
        space = v13.subcell_space(pq, box, 0)
        q = space.prior_transform(np.full(space.ndim, 0.5))
        assert math.isfinite(obs.chi2(space.theta(q), ctx))


def test_shape_space_holds_the_distance(synth):
    _, _, pq, box = synth
    for cell, e_rng in enumerate(v13.SHAPE_E_CELLS):
        s = v13.shape_space(pq, box, cell)
        assert s.ndim == 10 and s.held[1] == pytest.approx(math.log(v13.LE_VERRIER_A_AU))
        assert s.e_bounds == e_rng and s.inclined
    assert v13.shape_space(pq, box, 0, a_au=38.4).held[1] == pytest.approx(math.log(38.4))


def test_subcells_partition_cell_0(synth):
    from fork.inference.nested import cell_space
    from fork.inference.priors import free_a_space

    _, _, pq, box = synth
    whole = cell_space(free_a_space(pq, v13.S1_E_MAX, box, inclined=True), v13.S1_CELL_0)
    parts = [v13.subcell_space(pq, box, k) for k in range(len(v13.S1_SUBCELLS))]
    assert sum(p.slice_weight for p in parts) == pytest.approx(whole.slice_weight, rel=1e-12)
    assert all(p.ndim == 11 for p in parts)


def test_combine_subcells():
    w = np.array([0.2, 0.3, 0.5])
    lnz, err = v13.combine_subcells(np.full(3, -360.0), np.full(3, 0.1), w)
    assert lnz == pytest.approx(-360.0) and err == pytest.approx(0.1 * math.sqrt((w**2).sum()))
    lnz, _ = v13.combine_subcells(np.array([-360.0, -400.0]), np.zeros(2), np.array([0.5, 0.5]))
    assert lnz == pytest.approx(-360.0 + math.log(0.5), abs=1e-9)


def test_peirce_probability_and_ratio():
    rng = np.random.default_rng(5)
    wide = rng.normal(v13.GALLE_LON_DEG, 4.0, 40_000)
    narrow = rng.normal(v13.GALLE_LON_DEG, 1.0, 40_000)
    exact = lambda s: math.erf(v13.PEIRCE_HALF_WIDTH_DEG / (s * math.sqrt(2)))  # noqa: E731
    assert v13.peirce_probability(wide) == pytest.approx(exact(4.0), abs=0.01)
    assert v13.peirce_probability(np.array([325.0, 326.7, 326.8, 685.9])) == pytest.approx(0.5)
    r = v13.peirce_ratio(narrow, wide, rng, n_boot=300)
    assert r["ratio"] == pytest.approx(exact(1.0) / exact(4.0), rel=0.05)
    assert r["ratio_lo68"] < r["ratio"] < r["ratio_hi68"]
    assert v13.peirce_ratio(narrow, wide + 90.0, rng, n_boot=20)["ratio"] is None
    assert v13.log_bayes_factor(-360.0, 0.3, -362.0, 0.4) == pytest.approx((2.0, 0.5))


def test_region_summary_and_level():
    rng = np.random.default_rng(6)
    ra_h = (359.5 + rng.normal(0, 1.0, 50_000)) % 360.0 / 15.0  # straddles 0h
    dec = rng.normal(-10.0, 0.5, 50_000)
    s = v13.region_summary(ra_h, dec)
    assert s["ra_deg"]["median"] == pytest.approx(359.5, abs=0.05)
    assert s["dec_deg"]["hi68"] - s["dec_deg"]["lo68"] == pytest.approx(1.0, abs=0.03)
    assert v13.region_level(s, ra_h, dec, 359.5, -10.0) < 0.01
    assert v13.region_level(s, ra_h, dec, 359.5, -8.0) > 0.99
    one_sigma = v13.region_level(s, ra_h, dec, 0.5 / math.cos(math.radians(10)), -10.0)
    assert one_sigma == pytest.approx(1 - math.exp(-0.5), abs=0.02)


def test_predictive_places_on_prior_draws(synth):
    _, ctx, pq, box = synth
    space = v13.shape_space(pq, box, 0)
    rng = np.random.default_rng(7)
    eq = np.array([space.prior_transform(u) for u in rng.random((3, space.ndim))])
    out = v13.predictive_places(space, ctx, eq, {"galle": "1846-09-23T23:00:00"})
    ra, dec = out["galle"]["ra_h"], out["galle"]["dec_deg"]
    assert ra.shape == dec.shape == (3,)
    assert np.all((ra >= 0) & (ra < 24)) and np.all(np.abs(dec) < 60)


@pytest.mark.slow
def test_smoke_runs_on_synthetic_data():
    out = v13.smoke(maxiter=40, nlive=50)
    assert len(out) == 6 and all(v["niter"] > 0 for v in out.values())


def test_lalande_precovery_table_matches_both_entries_and_the_printed_differences():
    import csv
    from pathlib import Path

    d = Path(__file__).resolve().parents[1] / "data" / "historical" / "precoveries"
    read = lambda p: list(csv.DictReader((d / p).open(encoding="utf-8")))  # noqa: E731
    p1, p2 = read("passes/lalande_1795_pass1.csv"), read("passes/lalande_1795_pass2.csv")
    assert len(p1) == len(p2) == 9
    for a, b in zip(p1, p2, strict=True):
        assert all(a[k] == b[k] for k in a if k != "flag"), a
    final = read("lalande_1795.csv")
    cr = [r for r in p1 if r["witness"] == "cr24_p668" and r["item"] == "position"]
    reads = [r for r in p1 if r["item"] == "reads"]
    for f, c, e in zip(final, cr, reads, strict=True):
        assert (f["ra_app_deg"], f["ra_app_min"], f["ra_app_sec"]) == (
            c["ra_deg"],
            c["ra_min"],
            c["ra_sec"],
        )
        assert (f["dec_app_deg"], f["dec_app_min"], f["dec_app_sec"]) == (
            c["dec_deg"],
            c["dec_min"],
            c["dec_sec"],
        )
        h, m, s = (int(x) for x in f["paris_mean_time"].split(":"))
        assert (h, m, s) == (int(e["time_h"]), int(e["time_m"]), int(e["time_s"]))
    sec = lambda r, q: int(r[f"{q}_min"]) * 60 + float(r[f"{q}_sec"])  # noqa: E731
    diff = next(r for r in p1 if r["item"] == "difference")
    assert sec(cr[0], "ra") - sec(cr[1], "ra") == pytest.approx(sec(diff, "ra"), abs=1e-9)
    assert sec(cr[0], "dec") - sec(cr[1], "dec") == pytest.approx(sec(diff, "dec"), abs=1e-9)
