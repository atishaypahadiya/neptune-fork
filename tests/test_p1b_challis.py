"""Phase 1b step 4: the detection model's closed form, its sampler and the abduction."""

from __future__ import annotations

import numpy as np
import pytest

from fork.historical import challis_1846 as C


def _geo() -> dict:
    """A synthetic geometry: Neptune drifting west and south, one target at its place."""
    k = np.arange(C.N_NIGHTS)
    ra = 21.98 - 0.0016 * k
    dec = -12.9 - 0.009 * k
    return {
        "ra_h": ra,
        "dec_deg": dec,
        "moon_excluded": np.zeros(C.N_NIGHTS, dtype=bool),
        "on_chart_1800": np.zeros(C.N_NIGHTS, dtype=bool),
        "on_chart_date": np.ones(C.N_NIGHTS, dtype=bool),
        "targets": {
            name: {"ra_h": ra + off / 60.0, "dec_deg": dec + ddec / 60.0}
            for name, off, ddec in (("S4", -12.0, -30.0), ("S5", 6.0, 66.0), ("LVJ", 0.0, 0.0))
        },
    }


def _draws(n: int, seed: int = 5) -> list[C.Draw]:
    rng = np.random.default_rng(seed)
    return [C.draw(rng, _geo()) for _ in range(n)]


def test_night_indices() -> None:
    assert C.JUL29 == 0
    assert C.PRIMARY_LAST == 55
    assert C.SECONDARY_LAST == 63
    assert C.N_NIGHTS == 65


def test_draws_respect_the_priors_and_the_record() -> None:
    for d in _draws(300):
        for name, v in (("p_disc_routine", d.p_r), ("p_disc_looking", d.p_l),
                        ("stars_per_hour", d.rate), ("weather", d.weather)):  # fmt: skip
            lo, hi = C.PRIORS[name]
            assert lo <= v <= hi
        assert all(C.UNKNOWN_FIRST <= k <= C.UNKNOWN_LAST for k in d.unknown)
        assert len(d.unknown) == len(set(d.unknown)) <= d.n_unknown
        assert all(d.usable[list(C.KNOWN)])
        assert not d.usable[C.AUG4 + 1 : C.AUG12].any()
        assert not d.usable[C.SEP30]
        assert d.hours[C.SEP29] == C.SEP29_HOURS


def test_zone_priors_reproduce_condition_1() -> None:
    """Recorded on 12 August and not on 30 July, by declination or by right ascension."""
    geo = _geo()
    for d in _draws(300):
        assert C._zone_covers(d, geo, C.AUG12)
        assert not C._zone_covers(d, geo, C.JUL30)


def test_a0_closed_form() -> None:
    geo = _geo()
    for d in _draws(50):
        f1, f2 = C.cdfs(d, C.groups(d, geo, "A0"))
        assert f1[C.PRIMARY_LAST] == pytest.approx(1.0 - (1.0 - d.p_r) ** 2)
        assert f1[C.AUG4 - 1] == pytest.approx(0.0, abs=1e-15)
        assert f1[C.SEP29] == pytest.approx(1.0 - (1.0 - d.p_r) ** 2 * (1.0 - d.p_l))
        assert np.all(np.diff(f1) >= -1e-15) and np.all(np.diff(f2) >= -1e-15)
        assert np.all(f2 <= f1 + 1e-12)


@pytest.mark.parametrize("case", C.CASES)
def test_sampler_agrees_with_closed_form(case: str) -> None:
    geo = _geo()
    rng = np.random.default_rng(11)
    n = 20_000
    for d in _draws(8, seed=3):
        gs = C.groups(d, geo, case)
        f1, _ = C.cdfs(d, gs)
        t1 = C.sample(d, gs, rng, n)
        for night in (C.AUG12, C.PRIMARY_LAST, C.SECONDARY_LAST):
            p = f1[night]
            tol = 5.0 * np.sqrt(max(p * (1.0 - p), 1e-9) / n) + 1.0 / n
            assert abs(np.mean(t1 <= night) - p) < tol


def test_prompt_comparison_never_lowers_the_probability() -> None:
    geo = _geo()
    for d in _draws(100):
        a0, _ = C.cdfs(d, C.groups(d, geo, "A0"))
        for case in ("C1a", "C1b", "C2_date"):
            f, _ = C.cdfs(d, C.groups(d, geo, case))
            assert np.all(f >= a0 - 1e-12)
        late, _ = C.cdfs(d, C.groups(d, geo, "C1a_7d"))
        soon, _ = C.cdfs(d, C.groups(d, geo, "C1a"))
        assert np.all(soon >= late - 1e-12)
        off, _ = C.cdfs(d, C.groups(d, geo, "C2_1800"))  # planet off the sheet: A0
        np.testing.assert_allclose(off, a0, atol=1e-15)


def test_layout_from_the_true_place_records_the_planet_at_once() -> None:
    geo = _geo()
    for d in _draws(50):
        gs = C.groups(d, geo, "C4a")
        assert gs[0].rec[0][0] == C.JUL29
        assert gs[0].rec[0][1] == pytest.approx(d.p_rec)


def test_abduction_conditions() -> None:
    rng = np.random.default_rng(2)
    n_all = 0
    for d in _draws(2000, seed=9):
        a = C.accept(d, rng)
        assert a["all"] == (a["c23"] and a["c4"])
        assert a["c4"] == (C.STARS_TOTAL[0] <= C.total_stars(d) <= C.STARS_TOTAL[1])
        n_all += a["all"]
    assert 0.02 < n_all / 2000 < 0.3


def test_held_reading_keeps_the_realised_events() -> None:
    """A0 with the realised events held: nothing before 29 September, the disc on that night."""
    geo = _geo()
    for d in _draws(50):
        f1, f2 = C.cdfs(d, C.groups(d, geo, "A0", held=True))
        assert f1[C.PRIMARY_LAST] == 0.0
        assert f1[C.SEP29 - 1] == 0.0
        assert f1[C.SEP29] == 1.0
        assert f2[C.SECONDARY_LAST] == 0.0  # 30 September lost to the Moon
        for case in ("C1a", "C1b", "C2_date"):  # only the comparison is drawn
            g1, _ = C.cdfs(d, C.groups(d, geo, case, held=True))
            assert g1[C.PRIMARY_LAST] in (pytest.approx(0.0), pytest.approx(d.p_f),
                                          pytest.approx(1.0 - (1.0 - d.p_f) ** 2))  # fmt: skip
        for case in ("C3_S4", "C4a", "C4b"):  # before 29 September no event is shared
            a, _ = C.cdfs(d, C.groups(d, geo, case, held=True))
            b, _ = C.cdfs(d, C.groups(d, geo, case, held=False))
            np.testing.assert_allclose(a[: C.SEP29], b[: C.SEP29], atol=1e-15)


def test_c1_completion_date_arithmetic() -> None:
    geo = _geo()
    d = _draws(1)[0]
    d.desk_pair, d.desk_day = 6.0, 1.0
    assert C.desk_days(d, 0.5) == 3
    assert C.desk_days(d, 0.01) == 1
    d.rate, d.hours[C.AUG12] = 49.0, 2.0  # the planet half way through the series
    assert C._aug12_fraction(d) == pytest.approx(0.5)
    d.branch = "ra"
    g12 = C.groups(d, geo, "C1a", held=True)[1]
    assert g12.tc == C.AUG12 + 3  # starts the next day, three days of desk work
    assert C.groups(d, geo, "C1a_7d", held=True)[1].tc == C.AUG12 + 9
    f1, _ = C.cdfs(d, [g12])
    assert f1[g12.tc - 1] == 0.0
    assert f1[g12.tc] == pytest.approx(d.p_f)


def test_reading_rows() -> None:
    assert C.reading(0.85, 0.99).startswith("would very probably have")
    assert C.reading(0.0, 0.2) == "would very probably not have"
    assert C.reading(0.25, 0.75) == "an open chance"
    assert C.reading(0.1, 0.5) == "depends on"
    assert C.reading(0.5, 0.9) == "depends on"
    assert C.paired_reading(0.01, 0.9) == "improves on A0"
    assert C.paired_reading(0.0, 0.9) == "no difference established"
    assert C.paired_reading(-0.5, -0.1) == "does worse than A0"


def test_halved_priors() -> None:
    hp = C.halved(C.PRIORS)
    assert hp["p_disc_routine"] == pytest.approx((0.125, 0.375))
    assert hp["n_unknown_sessions"] == (13, 20)
    rng = np.random.default_rng(1)
    for _ in range(50):
        d = C.draw(rng, _geo(), hp)
        assert 0.125 <= d.p_r <= 0.375
        assert 13 <= d.n_unknown <= 20


def test_seeds_and_frozen_geometry() -> None:
    """The plan's seeds, and § 6.3 on the frozen pre-run: the 30 July / 12 August zone holds
    the planet on 12 August and not on 30 July; the sheet's first fallback misses it."""
    import json

    from fork.export import p1b_challis as X

    assert X.SEED_DRAWS == (2026, 11, 4, 1)
    assert X.SEED_EVENTS == (2026, 11, 4, 2)
    geo = X.geometry(json.loads(X.PRERUN_JSON.read_text()))
    rng = np.random.default_rng(list(X.SEED_DRAWS))
    for _ in range(200):
        d = C.draw(rng, geo)
        assert C._zone_covers(d, geo, C.AUG12)
        assert not C._zone_covers(d, geo, C.JUL30)
    assert not geo["on_chart_1800"][C.AUG4] and not geo["on_chart_1800"][C.AUG12]
    assert geo["on_chart_date"][C.AUG4] and geo["on_chart_date"][C.AUG12]
