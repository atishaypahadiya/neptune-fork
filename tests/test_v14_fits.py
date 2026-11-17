"""The data-choice fits and the shape model (docs/notes/M6c_variants.md, docs/notes/M6d.md):
rows, spaces, seeds and the readings on made-up inputs. No sampling here."""

import math

import numpy as np
import pytest

from fork.export import v14_fits as V


def test_fits_and_seeds_are_as_planned():
    assert V.ORDER == ("bouvard", "adams", "shape38", "shape36")
    assert V.FITS["bouvard"].rng(0, 1) == [2026, 15, 1846, 1, 1, 1]
    assert V.FITS["adams"].rng(3, 4, 7) == [2026, 15, 1846, 2, 4, 4, 7]
    assert V.FITS["shape38"].rng(1, 2) == [2026, 16, 1846, 1, 2, 2]
    assert V.FITS["shape36"].rng(-1, 0, 17) == [2026, 16, 1846, 2, 0, 0, 17]
    seeds = {tuple(f.rng(k, s)) for f in V.FITS.values() for k in range(f.n_cells) for s in V.SEEDS}
    assert len(seeds) == 4 * (4 + 4 + 2 + 2)
    assert [V.FITS["adams"].nlive(k) for k in range(4)] == [8000, 4000, 4000, 4000]
    assert [V.FITS["shape38"].nlive(k) for k in range(2)] == [4000, 4000]
    assert V.FITS["shape38"].held_a_au == pytest.approx(38.4)
    assert V.FITS["shape36"].held_a_au == pytest.approx(36.154)
    assert (V.SLICES, V.N_DRAWS, V.G4_DRAWS) == (28, 6000, 200)


@pytest.mark.parametrize(("name", "n"), [("bouvard", 492), ("adams", 527), ("shape38", 529)])
def test_rows_and_prior_masses(name, n):
    fit = V.FITS[name]
    obs, _, base, _, _ = V.problem(fit)
    assert obs.n_measurements == n
    w = np.exp(V.log_prior_masses(fit))
    assert w.sum() == pytest.approx(1.0, abs=1e-9)
    assert len(V.cell_space(fit, 0).unit_names) - 4 == fit.n_phi
    if fit.held_a_au is not None:
        th = base.theta(V.cell_space(fit, 1).prior_transform(np.full(4 + fit.n_phi, 0.3)))
        assert th[7] == pytest.approx(fit.held_a_au)
        assert 0.12 <= th[8] <= 0.3


def test_all_rows_reproduce_the_headline_observations():
    """The row reader with no rows dropped gives the headline's observation set."""
    from fork.export import m6c

    mine = V.real_observations("leverrier")
    ref = m6c.dataset()[0]
    assert mine.n_measurements == ref.n_measurements
    np.testing.assert_allclose(mine.value_rad, ref.value_rad, rtol=0, atol=1e-15)
    np.testing.assert_allclose(mine.sigma_arcsec, ref.sigma_arcsec)


HEAD = {"median": 325.735, "lo68": 323.511, "hi68": 327.672, "lo95": 319.593, "hi95": 330.213,
        "width68": 4.161}  # fmt: skip


def _lon(median: float, width: float) -> dict:
    return {"median": median, "width68": width, "vs": {"jpl": {"in_68": True, "in_95": True}}}


@pytest.mark.parametrize(
    ("median", "wording"),
    [
        (325.9, "does not move"),
        (326.5, "within the headline's 68 %"),
        (324.0, "within the headline's 68 %"),
        (328.5, "outside the headline's 68 %"),
        (320.5, "outside the headline's 68 %"),
        (331.0, "outside the headline's 95 %"),
        (318.0, "outside the headline's 95 %"),
    ],
)
def test_shift_reading(median, wording):
    assert wording in V.read_shift(_lon(median, 4.2), HEAD)["shift_wording"]


def test_width_reading_and_thresholds():
    assert V.read_shift(_lon(325.8, 4.161 * 1.2), HEAD)["width_wording"] == "wider"
    assert V.read_shift(_lon(325.8, 4.161 * 0.8), HEAD)["width_wording"] == "narrower"
    assert V.read_shift(_lon(325.8, 4.161), HEAD)["width_wording"] == "the same width"
    assert V.read_shift(_lon(325.735 + 0.36, 4.161 * 1.103), HEAD)["shift_at_a_threshold"]
    assert V.read_shift(_lon(325.735 + 0.36, 4.161 * 1.103), HEAD)["width_at_a_threshold"]
    assert not V.read_shift(_lon(326.6, 4.3), HEAD)["shift_at_a_threshold"]


@pytest.mark.parametrize(
    ("ln_b", "err", "wording"),
    [
        (0.2, 0.2, "the data do not prefer either model"),
        (2.0, 0.2, "a weak preference for the shape model"),
        (-2.0, 0.2, "a weak preference for free distance"),
        (-4.0, 0.3, "a strong preference for free distance"),
        (7.0, 0.5, "a decisive preference for the shape model"),
        (0.9, 0.2, "between"),
        (-3.1, 0.3, "between"),
    ],
)
def test_bayes_factor_reading(ln_b, err, wording):
    assert V.read_ln_b(ln_b, err)["wording"].startswith(wording)


def test_ratio_reading():
    r = {"p_shape": 0.4, "p_free": 0.37, "ratio": 1.08, "ratio_lo68": 0.98, "ratio_hi68": 1.2}
    assert "did not change" in V.read_ratio(r, [0.39, 0.4, 0.41, 0.4])["wording"]
    r = {"p_shape": 0.6, "p_free": 0.37, "ratio": 1.62, "ratio_lo68": 1.5, "ratio_hi68": 1.75}
    assert V.read_ratio(r, [0.6, 0.61, 0.59, 0.6])["wording"] == (
        "made it more probable, by a factor 1.62"
    )
    r = {"p_shape": 0.1, "p_free": 0.4, "ratio": 0.25, "ratio_lo68": 0.2, "ratio_hi68": 0.3}
    assert V.read_ratio(r, [0.1, 0.1, 0.1, 0.1])["wording"] == (
        "made it less probable, by a factor 4.00"
    )
    assert V.read_ratio(r, [0.1, 0.25, 0.1, 0.1])["wording"] == "not repeatable; not read"


def test_g4_verdict():
    ok = np.c_[np.full(200, 0.01), np.full(200, 0.05)]
    assert V.g4_verdict(ok)["pass"]
    some = ok.copy()
    some[:9, 0] = 0.3  # 4.5 % outside, none beyond five times the tolerance
    assert V.g4_verdict(some)["pass"]
    some[:11, 0] = 0.3
    assert not V.g4_verdict(some)["pass"]
    far = ok.copy()
    far[0, 1] = 1.5
    assert not V.g4_verdict(far)["pass"]


def test_the_published_free_distance_fit_is_reproduced():
    # Reads the published fit's reweighted draws from data/cache, which CI does not have.
    try:
        lon, pub = V.free_distance_galle()
    except FileNotFoundError as e:
        pytest.skip(f"the published fit's draws are not in data/cache ({e})")
    assert len(lon) == 6000
    assert pub["ln_z"] == pytest.approx(-361.5227, abs=1e-3)
    assert math.isfinite(float(np.mean(lon)))


def _walk_keys(o):
    if isinstance(o, dict):
        for k, v in o.items():
            yield k
            yield from _walk_keys(v)
    elif isinstance(o, list):
        for v in o:
            yield from _walk_keys(v)


def _unread_summary(gate_fails=False):
    lon = {"median": 316.204, "lo68": 312.606, "hi68": 320.775, "width68": 8.17,
           "vs": {"jpl": {"target": 325.873, "median_minus_target_deg": -9.669}}}  # fmt: skip
    fit = {"ln_z": -355.9, "galle_lon_ecl_of_date_deg": lon, "p_within_0p87_deg": 0.018,
           "mass_earth": {"median": 33.1}}  # fmt: skip
    cell = {"cell": 2, "material": True, "width_pass": True, "ln_z_pass": True,
            "max_rel_dev": 0.04, "ln_z_range": 0.1}  # fmt: skip
    return {
        "fit": "adams", "held_a_au": None, "read": False, "fits": {"merged": fit, "s1": fit},
        "gates": {g: {"pass": not (gate_fails and g == "G2")} for g in ("G1", "G2", "G3", "G4")},
        "T1": {"cells": [cell], "median_range_single_seed_fits": 0.353,
               "ln_z_range_single_seed_fits": 0.064, "pass": False},
        "std_95_bounds": {"lo95": 0.09, "hi95": 0.33}, "label_95": "final",
        "reading": {"not_read": "a gate or T1 failed", "values": {
            "shift_deg": -9.531, "shift_wording": "outside", "width_ratio": 1.963,
            "width_wording": "wider", "jpl_in_68": False, "jpl_in_95": False}},
    }  # fmt: skip


def test_an_unread_fit_carries_no_comparison_or_band():
    """§ 5: a fit not read is written with its status and label, its posterior summaries, and no
    headline comparison, band wording, true-place field or final/provisional label."""
    out = V.unread_public(_unread_summary())
    banned = {"reading", "label_95", "vs", "jpl", "target", "median_minus_target_deg",
              "p_within_0p87_deg", "shift_deg", "shift_wording", "width_ratio", "width_wording",
              "jpl_in_68", "jpl_in_95", "not_read", "values"}  # fmt: skip
    assert not banned & set(_walk_keys(out))
    assert out["status"] == "not read (T1 fails: median range 0.353° vs 0.35°)"
    assert out["label"] == "not repeatable"
    assert out["fits"]["merged"]["galle_lon_ecl_of_date_deg"]["median"] == 316.204
    assert out["fits"]["s1"]["mass_earth"] == {"median": 33.1}
    gate = V.unread_public(_unread_summary(gate_fails=True))
    assert gate["label"] == "not validated"
    assert gate["status"].startswith("not read (G2 fails; T1 fails")
