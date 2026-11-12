"""Phase 1b step 3b, stages B and C (docs/notes/P1B_adams_reproduction.md § 4): Adams's printed
intermediate numbers from his printed inputs, and his solution from three entry points. Every
tolerance is the plan's; where a planned tolerance is missed the test asserts the achieved value
instead and names the miss here and in the note (the plan's "report the number, don't loosen
silently"): see ``KNOWN_MISSES`` below."""

import math

import numpy as np
import pytest

from fork.historical import adams_1846 as A


@pytest.fixture(scope="module")
def memoir():
    return A.load_memoir()


# ----------------------------------------------------------------------------- the print


def test_the_transcription_has_the_pre_registered_shape(memoir):
    assert len(memoir.places) == 30
    assert sum(p["kind"] == "ancient" for p in memoir.places) == 9
    assert len(A.modern_rows(memoir)) == 21
    assert [t for t, _ in A.modern_rows(memoir)] == list(range(-10, 11))
    assert len(memoir.equations) == 1407 - sum(len(v) for v in memoir.equations.values()) + len(
        memoir.equations
    )
    for hyp, secs in ((1, (15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29)),
                      (2, (35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46))):  # fmt: skip
        for sec in secs:
            assert memoir.eqs(hyp, sec), (hyp, sec)


def test_angles_and_logarithms_parse_as_printed():
    assert A.ang("13 0.5") == pytest.approx(13 + 0.5 / 60)
    assert A.ang("-51 30") == pytest.approx(-51.5)
    assert A.antilog("+(1.24782)") == pytest.approx(10**1.24782)
    assert A.antilog("-(9.48254)") == pytest.approx(-(10 ** (9.48254 - 10)))
    assert A.antilog("(1.65083)") > 0  # printed unsigned


def test_section_22_relations_are_e_over_4_times_cos_and_sin_beta(memoir):
    R = A.relations(memoir.const(1, 22, "e"), memoir.angle(1, 22, "eps_minus_varpi"))
    assert R["r11"] == pytest.approx(0.007460, abs=1e-6)
    assert R["r12"] == pytest.approx(0.008974, abs=1e-6)
    assert R["r21"] == pytest.approx(-0.008974, abs=1e-6)
    assert R["r22"] == pytest.approx(0.007460, abs=1e-6)


# ----------------------------------------------------------------------------- stage B

# Planned tolerances missed, with the achieved value (reported in the note's B tables).
KNOWN_MISSES = {
    # B1 (plan 3e-5 in log10): his second derivatives and every Hypothesis II value sit above it
    "b1_log10": 2.0e-4,
    # B2 (plan 0.02″): the i = 3 terms with the small divisor, and three Hypothesis II terms
    "b2_arcsec": 0.09,
    # B3 (plan 0.0005 on coefficients): his sums Σ t sin mt, which are the (n) row and, × 10, the
    # δn column of the second group, are low by 0.0005–0.0007 (q₁ by 0.004) in Hypothesis I and
    # the δy₁ one by 0.0005 in both; everything else holds the plan's 0.0005
    "b3_coef": 0.0075,
    # B4 (plan 0.0005): single entries of the 1764 equation
    "b4_coef": 0.0014,
    # B5 (plan 0.003 / 0.3″): § 24's 1764 h₁ (0.0082), § 44's e2 k₃ (0.025, his unrounded lhs of
    # (2)), § 41's δn (0.005, his rounding); § 37's comb1 prints 111.07 for 141.07 (a misprint
    # that § 38 does not carry)
    "b5_coef": 0.026,
    "b5_lhs_misprint": 30.05,
}


def test_b1_laplace_coefficients_against_his_twelve_logarithms(memoir):
    rows = A.b1_table(memoir)
    assert len(rows) == 24
    worst = max(abs(r["delta"]) for r in rows)
    assert worst <= KNOWN_MISSES["b1_log10"], worst
    # the plan's 3e-5 holds for b and α db/dα of Hypothesis I except b₃ (3.8e-5)
    h1_first = [r for r in rows if r["hyp"] == 1 and not r["item"].startswith("log_a2d2b")]
    assert max(abs(r["delta"]) for r in h1_first) <= 4e-5
    assert max(abs(r["delta"]) for r in h1_first if r["item"] != "log_b3") <= 3e-5


def test_b2_series_constants_from_section_11(memoir):
    rows = [r for r in A.b2_table(memoir) if not r["item"].startswith("§28")]
    assert len(rows) == 18
    worst = max(abs(r["delta"]) for r in rows)
    assert worst <= KNOWN_MISSES["b2_arcsec"], worst
    within = [r for r in rows if abs(r["delta"]) <= 0.02]
    assert len(within) >= 12
    # his printed Laplace values do not remove the larger misses: the cause is downstream of B1
    for r in rows:
        if abs(r["delta"]) > 0.03:
            assert abs(r["delta_his_logs"]) > 0.02, r


def test_b2_section_28_constants_follow_from_the_series(memoir):
    rows = [r for r in A.b2_table(memoir) if r["item"].startswith("§28")]
    assert len(rows) == 18
    for r in rows:
        # printed to 0.01 (the large constants), to 0.01 with one figure (0.18, 0.20), or to
        # 5–6 decimals (the P, Q multipliers; his 0.007460/0.008974 are e/4 cos β, e/4 sin β)
        tol = 0.02 if abs(r["printed"]) > 1 else 0.006 if abs(r["printed"]) >= 0.1 else 1e-4
        assert abs(r["delta"]) <= tol, r


@pytest.mark.parametrize("hyp", [1, 2])
def test_b3_normal_equations_from_the_21_equations_of_condition(memoir, hyp):
    rows = A.b3_table(memoir, hyp)
    coef = [r for r in rows if r["term"] != "lhs"]
    rhs = [r for r in rows if r["term"] == "lhs"]
    assert max(abs(r["delta"]) for r in rhs) <= 0.05
    worst = max(coef, key=lambda r: abs(r["delta"]))
    assert abs(worst["delta"]) <= KNOWN_MISSES["b3_coef"], worst
    # the plan's 0.0005 holds for every coefficient outside the (n) row and the δn column
    rest = [r for r in coef if not (r["equation"] == "n" or r["term"] == "dn")]
    assert max(abs(r["delta"]) for r in rest) <= 0.0005
    assert memoir.eq(hyp, 18 if hyp == 1 else 36, "n")["dy1"] == 9.3938  # exact 9.3943


@pytest.mark.parametrize("hyp", [1, 2])
def test_b4_ancient_equations_from_their_epochs(memoir, hyp):
    rows = [r for r in A.b4_table(memoir, hyp) if r["term"] not in ("lhs", "de", "dn")]
    worst = max(rows, key=lambda r: abs(r["delta"]))
    assert abs(worst["delta"]) <= KNOWN_MISSES["b4_coef"], worst
    rest = [r for r in rows if r["equation"] != "anc_1764"]
    assert max(abs(r["delta"]) for r in rest) <= 0.0005


@pytest.mark.parametrize("hyp", [1, 2])
def test_b5_each_printed_section_from_the_printed_ones_before_it(memoir, hyp):
    table = A.b5_table(memoir, hyp)
    for sec, c in table.items():
        if sec in (29, 46):
            assert c["sign_agreements"]
            assert c["log10_max"] <= 3e-4, (sec, c["log10_max"])
            continue
        assert c["coef_max"] <= KNOWN_MISSES["b5_coef"], (sec, c["coef_max"])
        if sec in (37, 38):
            assert c["lhs_max"] <= KNOWN_MISSES["b5_lhs_misprint"]
            assert c["rows"]["comb1"]["lhs"] > 29.9  # the misprint
        else:
            assert c["lhs_max"] <= 0.3, (sec, c["lhs_max"])
    clean = [c["coef_max"] for sec, c in table.items() if sec not in (24, 29, 41, 44, 46)]
    assert max(clean) <= 0.003


def test_b5_his_multipliers_are_recovered(memoir):
    c = A.compare_step(memoir, 1, 22)
    assert c["multipliers"]["1"] == pytest.approx((0.006768, 0.040287), abs=5e-6)
    assert c["multipliers"]["2"] == pytest.approx((-0.001869, 0.008187), abs=5e-6)
    c = A.compare_step(memoir, 2, 39)
    assert c["multipliers"]["1"] == pytest.approx((-1.7106, -0.03607), abs=5e-5)
    assert c["multipliers"]["2"] == pytest.approx((0.00165, -4.0487), abs=5e-4)


def test_b5_the_section_37_misprint_does_not_reach_section_38(memoir):
    """§ 37 prints 111.07 where −(h₁) + (h₂) + (h₃) gives 141.07; § 38's 268.02 follows from
    141.07."""
    c37 = A.compare_step(memoir, 2, 37)
    assert c37["rows"]["comb1"]["lhs"] == pytest.approx(30.0, abs=0.01)
    mine, _ = A.reproduce_step(memoir, 2, 38)
    from_141 = A.eliminate(
        {**memoir.eq(2, 37, "comb1"), "lhs": 141.07}, memoir.eq(2, 35, "eps"), "de"
    )
    assert from_141["lhs"] == pytest.approx(memoir.eq(2, 38, "comb1")["lhs"], abs=0.02)
    assert abs(mine["comb1"]["lhs"] - memoir.eq(2, 38, "comb1")["lhs"]) > 29


# ----------------------------------------------------------------------------- stage C


@pytest.fixture(scope="module")
def ctable(memoir):
    return A.c_table(memoir)


@pytest.mark.parametrize("hyp", [1, 2])
def test_c0_his_printed_final_equations_give_his_solution_at_tier_1(ctable, hyp):
    rec = ctable[f"H{hyp}/C0"]
    assert rec["vs_print"]["tier"] == 1, rec["vs_print"]
    assert rec["residual_tier"] in (1, 2)


def test_c1_his_printed_normal_equations_give_his_solution_within_his_rounding(ctable):
    """Hypothesis I at tier 2 (θ 8′, mean longitude 8′, ϖ′ 19′); Hypothesis II misses tier 2 on
    θ alone (13′ against 10′), the other quantities at tier 2."""
    rec = ctable["H1/C1"]
    assert rec["vs_print"]["tier"] == 2, rec["vs_print"]
    rec = ctable["H2/C1"]
    rows = rec["vs_print"]["rows"]
    assert rows["theta"]["delta"] < 15 / 60
    assert all(v["tier"] in (1, 2) for k, v in rows.items() if k != "theta"), rows


@pytest.mark.parametrize("hyp", [1, 2])
def test_c2_the_data_give_brookes_not_adams(ctable, hyp):
    rec = ctable[f"H{hyp}/C2"]
    assert rec["vs_print"]["tier"] == 0, rec["vs_print"]
    assert rec["vs_brookes"]["tier"] in (1, 2), rec["vs_brookes"]
    assert rec["vs_brookes"]["rows"]["mean_lon_1846"]["delta"] <= 2 / 60


def test_the_roots_include_one_physical_solution_per_hypothesis(ctable):
    for key, rec in ctable.items():
        positive = [r for r in rec["roots"] if r["m_prime"] > 0 and r["e_prime_units"] < 20]
        assert len(positive) >= 1, key


@pytest.mark.parametrize("hyp", [1, 2])
def test_the_gap_sits_in_the_second_group_normal_equation_n(memoir, hyp):
    """Swapping (n)'s printed coefficients for the exact sums moves θ by more than 1°; every
    other equation moves it by less than 0.15°; one unit of his fourth decimal on (n)'s δy₁
    coefficient moves θ by more than 1°."""
    s = A.sensitivity_table(memoir, hyp)
    n_shift = abs(s["equations"]["n"]["coefficients"]["theta"])
    assert n_shift > 1.0, n_shift
    others = max(abs(v["whole"]["theta"]) for k, v in s["equations"].items() if k != "n")
    assert others < 0.15, others
    assert abs(s["n_coefficients"]["dy1"]["dtheta_per_step"]) > 1.0
    assert abs(s["group2_exact"]["theta"] - s["base"]["theta"]) > 1.0
    assert abs(s["group1_exact"]["theta"] - s["base"]["theta"]) < 0.05


def test_true_longitude_reduces_to_the_mean_longitude_for_a_circle():
    assert A.true_longitude(123.4, 0.0, 300.0) == pytest.approx(123.4)
    # e = 0.16, M = 10°: the equation of the centre to third order
    e, M = 0.16, math.radians(10)
    series = (2 * e - e**3 / 4) * math.sin(M) + 1.25 * e**2 * math.sin(2 * M)
    series += 13 / 12 * e**3 * math.sin(3 * M) + 103 / 96 * e**4 * math.sin(4 * M)
    assert A.true_longitude(310.0, e, 300.0) - 310.0 == pytest.approx(
        math.degrees(series), abs=0.01
    )


def test_normal_equation_sums_use_his_closed_forms(memoir):
    """§ 14: Σ 2 t sin mt = ((t+1) sin mt − t sin m(t+1)) / (2 sin² ½m) at t = 10 equals the
    direct sum; his printed 9.3938 for Σ t sin(13° 0′.5 t)/10 is 9.3943."""
    m = math.radians(13 + 0.5 / 60)
    closed = (11 * math.sin(10 * m) - 10 * math.sin(11 * m)) / (2 * math.sin(m / 2) ** 2)
    t = np.arange(-10, 11)
    assert closed == pytest.approx(float(np.sum(t * np.sin(m * t))), abs=1e-9)
    assert closed / 10 == pytest.approx(9.3943, abs=5e-5)
    assert memoir.eq(1, 18, "n")["dy1"] == 9.3938
