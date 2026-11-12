"""Adams's 1846 theory of Uranus's perturbations, rebuilt from the memoir's printed numbers.

Phase 1b step 3b (docs/notes/P1B_adams_reproduction.md, pre-registered). The inputs are the
transcribed tables in ``data/historical/adams_challis/adams_memoir_*.csv``; nothing here is fitted
to the sky. Pure numpy, no numba: the problem is 30 equations in 16 unknowns.

Notation follows the memoir (pamphlet printed page = PDF page − 6):

- § 10: the correction of Uranus's mean longitude from the element corrections δε, δn, δx₁, δy₁
  with δx₂, δy₂ tied to δx₁, δy₁ by the factor e/4 (§ 22).
- § 11: the perturbation in Pontécoulant's form, F_i, G_i, H_i in the Laplace coefficients.
- § 12/§ 33: the eleven series constants; § 13/§ 34: the arguments and the unit of time.
- §§ 14–27 / 35–44: normal equations, eliminations, the ancient equations, the final four
  equations (1)–(4); § 28/§ 45: h, k, p, q in θ = ε − ε′, β = ε − ϖ, β′ = ε − ϖ′, m′, e′;
  § 29/§ 46: the three logarithmic equations; § 30/§ 47: the solution.

Equations are dicts ``term -> value`` with the right-hand side under ``"lhs"`` (the memoir writes
the observed quantity on the left); all values in arcseconds, coefficients dimensionless, time in
his unit (three synodic periods, 3.0362 yr) from 1810.328. m′ is in units of 1/5000 of the Sun,
e′ in units of 1/20, as he has them.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "historical" / "adams_challis"

GROUP1 = ("de", "dx1", "dx2", "h1", "h2", "h3", "p1", "p2", "p3")
GROUP2 = ("dn", "dy1", "dy2", "k1", "k2", "k3", "q1", "q2", "q3")
HKPQ = ("h1", "h2", "h3", "k1", "k2", "k3", "p1", "p2", "p3", "q1", "q2", "q3")
LOGTERMS = (
    "sin1",
    "cos1",
    "sin2",
    "cos2",
    "sin3",
    "cos3",
    "P3",
    "Q3",
    "P1c-Q1s",
    "P1s+Q1c",
    "P2c-Q2s",
    "P2s+Q2c",
)
ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
MASS_UNIT = 1.0 / 5000.0  # m′ = 1 is this fraction of the Sun
ECC_UNIT = 1.0 / 20.0  # e′ = 1 is this eccentricity
# the normal equations he formed (§ 15, § 18; § 35, § 36)
NORMAL_H1 = {
    1: ("eps", "x", "h1", "h2", "h3", "p2", "p3"),
    2: ("n", "y", "k1", "k2", "k3", "q2", "q3"),
}
NORMAL_H2 = {1: ("eps", "x", "h1", "h2", "h3"), 2: ("n", "y", "k1", "k2", "k3")}
PIVOT = {"eps": "de", "x": "dx1", "n": "dn", "y": "dy1"}
ANCIENT_YEARS = (1690, 1712, 1715, 1750, 1753, 1756, 1764, 1769, 1771)


# ----------------------------------------------------------------------------- the print


def _num(s: str) -> float:
    return float(
        s.replace("−", "-").replace("+", "", 1) if s.startswith("+") else s.replace("−", "-")
    )


def ang(s: str) -> float:
    """``"13 0.5"`` → 13.008333° ; ``"-51 30"`` → −51.5° ; ``"30"`` → 30° (callers scale)."""
    parts = s.replace("−", "-").split()
    d = float(parts[0])
    m = float(parts[1]) if len(parts) > 1 else 0.0
    sign = -1.0 if parts[0].startswith("-") else 1.0
    return sign * (abs(d) + m / 60.0)


def antilog(s: str) -> float:
    """``"+(1.24782)"`` → +10**1.24782 ; ``"(1.65083)"`` (unsigned) → positive."""
    s = s.strip().replace("−", "-")
    sign = -1.0 if s.startswith("-") else 1.0
    x = float(s.strip("+-()"))
    return sign * 10 ** (x - 10.0 if x >= 5.0 else x)  # his logarithms are augmented by 10


@dataclass(frozen=True)
class Memoir:
    """The transcribed tables, with printed strings kept beside parsed values."""

    places: list[dict]
    series: dict[tuple[int, int, str], str]
    equations: dict[tuple[int, int, str], dict[str, str]]
    solution: dict[tuple[int, int, str], str]

    def eq(self, hyp: int, sec: int, name: str) -> dict[str, float]:
        raw = self.equations[(hyp, sec, name)]
        return {
            k: (antilog(v) if v.lstrip("+-").startswith("(") else _num(v)) for k, v in raw.items()
        }

    def eqs(self, hyp: int, sec: int) -> dict[str, dict[str, float]]:
        names = [n for (h, s, n) in self.equations if h == hyp and s == sec]
        return {n: self.eq(hyp, sec, n) for n in names}

    def const(self, hyp: int, sec: int, item: str) -> float:
        return _num(self.series[(hyp, sec, item)])

    def angle(self, hyp: int, sec: int, item: str) -> float:
        return ang(self.series[(hyp, sec, item)])

    def sol(self, hyp: int, sec: int, item: str) -> str:
        return self.solution[(hyp, sec, item)]


def load_memoir(data: Path = DATA) -> Memoir:
    def rows(name):
        with open(data / f"adams_memoir_{name}.csv", newline="") as f:
            return list(csv.DictReader(f))

    places = []
    for r in rows("normal_places"):
        p = {"year": int(r["year"]), "kind": r["kind"]}
        for k, v in r.items():
            if k not in ("year", "kind"):
                p[k] = _num(v) if v != "" else None
        places.append(p)
    series = {
        (int(r["hypothesis"]), int(r["section"]), r["item"]): r["printed"] for r in rows("series")
    }
    equations: dict[tuple[int, int, str], dict[str, str]] = {}
    for r in rows("equations"):
        equations.setdefault((int(r["hypothesis"]), int(r["section"]), r["equation"]), {})[
            r["term"]
        ] = r["printed"]
    solution = {
        (int(r["hypothesis"]), int(r["section"]), r["item"]): r["printed"] for r in rows("solution")
    }
    return Memoir(places, series, equations, solution)


# ----------------------------------------------------------------------------- his theory


def laplace_b(i: int, alpha: float, n_quad: int = 8192) -> tuple[float, float, float]:
    """b^(i)_{1/2}(α) = (1/π)∫₀^{2π} cos iφ (1 − 2α cos φ + α²)^{-1/2} dφ and the printed
    companions α db/dα and α² d²b/dα², by the trapezoid rule on the periodic integrand."""
    ph = np.linspace(0.0, 2.0 * np.pi, n_quad, endpoint=False)
    c = np.cos(ph)
    d = 1.0 - 2.0 * alpha * c + alpha * alpha
    w = np.cos(i * ph)
    b0 = 2.0 * np.mean(w * d**-0.5)
    b1 = 2.0 * np.mean(w * (-0.5) * d**-1.5 * (2.0 * alpha - 2.0 * c))
    b2 = 2.0 * np.mean(w * (0.75 * d**-2.5 * (2.0 * alpha - 2.0 * c) ** 2 - d**-1.5))
    return float(b0), float(alpha * b1), float(alpha * alpha * b2)


@dataclass(frozen=True)
class Theory:
    """The constants of § 12–13 (or § 33–34) in working form."""

    alpha: float
    n_prime_over_n: float  # n′/n
    c: dict[str, float]  # the eleven series constants, arcsec per unit m′ (and per e′ unit)
    rate: dict[str, float]  # argument rates, degrees per unit of time, by unknown
    n_deg: float  # n, degrees per unit (δa/a needs it)

    @property
    def unknown_rates(self) -> dict[str, float]:
        return self.rate


def _A(i: int, alpha: float, b: tuple[float, float, float]) -> tuple[float, float, float]:
    """a A_i, a² dA_i/da, a³ d²A_i/da² with A_i = b_i/a′ (i ≥ 2) and (b₁ − α)/a′ (the indirect
    part)."""
    b0, adb, a2d2b = b
    if i == 1:
        return alpha * (b0 - alpha), alpha * (adb - alpha), alpha * a2d2b
    return alpha * b0, alpha * adb, alpha * a2d2b


def series_constants(
    alpha: float,
    n_prime_over_n: float,
    e_uranus: float,
    b: dict[int, tuple[float, float, float]] | None = None,
) -> dict[str, float]:
    """The nine first-order constants of § 12 from § 11 (F_i, G_i, H_i), in arcsec per unit m′;
    the e′ terms per unit e′ (1/20). ``b`` may carry his printed Laplace values; default is B1's."""
    n = 1.0
    bb = b or {i: laplace_b(i, alpha) for i in range(4)}
    out = {}
    k_f = ARCSEC_PER_RAD * MASS_UNIT  # m′/2 Σ_{±i} F_i sin = m′ F_i sin iψ
    k_g = ARCSEC_PER_RAD * MASS_UNIT * e_uranus
    k_h = ARCSEC_PER_RAD * MASS_UNIT * ECC_UNIT
    for i in (1, 2, 3):
        z = i * (n - n_prime_over_n)
        aA, adA, ad2A = _A(i, alpha, bb[i])
        f = (
            3 * i * n**4 / (z * z * (z * z - n * n)) + i * n * n / (z * z - n * n)
        ) * aA + 2 * n**3 / (z * (z * z - n * n)) * adA
        g = (
            -3 * i * (i - 1) * n**4 / ((z - n) ** 2 * z * (z - 2 * n))
            - i * (i + 1) * n * n / (z * (z - 2 * n))
            + i * n * n / (z * z - n * n)
            + 3 * i * n**3 / (z * (z - n) * (z - 2 * n))
        ) * aA
        g += (
            -1.5 * (i - 1) * n**4 / ((z - n) ** 2 * z * (z - 2 * n))
            - 0.5 * (i - 1) * n * n / (z * (z - 2 * n))
            - 0.5 * n * n / (z * z - n * n)
            - 2 * i * n**3 / (z * (z - n) * (z - 2 * n))
        ) * adA
        g -= n**3 / (z * (z - n) * (z - 2 * n)) * ad2A
        aA1, adA1, ad2A1 = (
            _A(i - 1, alpha, bb[i - 1])
            if i > 1
            else (alpha * bb[0][0], alpha * bb[0][1], alpha * bb[0][2])
        )
        h = (
            1.5 * (i - 1) * (2 * i - 1) * n**4 / ((z - n) ** 2 * z * (z - 2 * n))
            + 0.5 * (i - 1) * (2 * i - 1) * n * n / (z * (z - 2 * n))
        ) * aA1
        h += (
            1.5 * (i - 1) * n**4 / ((z - n) ** 2 * z * (z - 2 * n))
            + 0.5 * (i - 1) * n * n / (z * (z - 2 * n))
            + 2 * i * n**3 / (z * (z - n) * (z - 2 * n))
        ) * adA1
        h += n**3 / (z * (z - n) * (z - 2 * n)) * ad2A1
        out[f"c_sin{i}"] = k_f * f
        # i = 1: the argument is −(n′t + ε′ − ϖ), so the printed coefficient has the opposite sign
        sgn = -1.0 if i == 1 else 1.0
        out[f"c_p{i}"] = sgn * k_g * g
        out[f"c_p{i}e"] = sgn * k_h * h
    return out


def theory_from_print(m: Memoir, hyp: int) -> Theory:
    sec_c, sec_a = (12, 13) if hyp == 1 else (33, 34)
    c = {
        k: m.const(hyp, sec_c, k)
        for k in (
            "c_sin1",
            "c_sin2",
            "c_sin3",
            "c_p1",
            "c_p1e",
            "c_p2",
            "c_p2e",
            "c_p3",
            "c_p3e",
            "c_sin3_e2",
            "c_sin3_e2e",
        )
    }
    rate = {u: m.angle(hyp, sec_a, f"arg_{u}") for u in GROUP1[1:] + GROUP2[1:]}
    n = m.angle(hyp, sec_a, "n")
    n_prime = m.angle(hyp, sec_a, "n_prime")
    alpha = 0.5 if hyp == 1 else 0.515
    return Theory(alpha, n_prime / n, c, rate, n)


def hkpq_constants(c: dict[str, float]) -> dict[str, float]:
    """The constants of § 28 / § 45 from the eleven series constants (derivation in the note:
    the e′ terms are written through P = p₃/m′ and Q = q₃/m′; the second-order terms fold into
    h₃, k₃)."""
    A, B = c["c_p3"], c["c_p3e"]
    r = c["c_sin3_e2e"] / B  # the e e′ term through P, Q
    return {
        "k_dh1": c["c_sin1"],
        "k_dh2": c["c_sin2"],
        "k_dh3_a": c["c_sin3"] + c["c_sin3_e2"] - c["c_sin3_e2e"] * A / B,
        "k_dh3_r": r,  # h₃ += r (cos β P + sin β Q); k₃ += r (−sin β P + cos β Q)
        "k_dp1_a": -(c["c_p1"] - c["c_p1e"] * A / B),  # × sin(θ − β)
        "k_dp1_b": -c["c_p1e"] / B,  # printed −0.046247 {P cos 2θ − Q sin 2θ}
        "k_dp2_a": c["c_p2"] - c["c_p2e"] * A / B,  # × sin(2θ − β)
        "k_dp2_b": c["c_p2e"] / B,  # printed +0.13055 {P cos θ − Q sin θ}
    }


def hkpq_over_m(
    theta: float, P: float, Q: float, beta: float, k: dict[str, float], r_cos: float, r_sin: float
) -> dict[str, float]:
    """h, k, p, q divided by m′ at θ (deg), with P = p₃/m′, Q = q₃/m′; § 28 as printed, with the
    h₃/k₃ cross terms r_cos = 0.007460, r_sin = 0.008974 (his values) or derived."""
    t = math.radians(theta)
    b = math.radians(beta)
    s1, c1, s2, c2, s3, c3 = (
        math.sin(t),
        math.cos(t),
        math.sin(2 * t),
        math.cos(2 * t),
        math.sin(3 * t),
        math.cos(3 * t),
    )
    return {
        "h1": k["k_dh1"] * s1,
        "k1": k["k_dh1"] * c1,
        "h2": k["k_dh2"] * s2,
        "k2": k["k_dh2"] * c2,
        "h3": k["k_dh3_a"] * s3 + r_cos * P + r_sin * Q,
        "k3": k["k_dh3_a"] * c3 - r_sin * P + r_cos * Q,
        "p1": k["k_dp1_a"] * math.sin(t - b) + k["k_dp1_b"] * (P * c2 - Q * s2),
        "q1": -k["k_dp1_a"] * math.cos(t - b) - k["k_dp1_b"] * (P * s2 + Q * c2),
        "p2": k["k_dp2_a"] * math.sin(2 * t - b) + k["k_dp2_b"] * (P * c1 - Q * s1),
        "q2": k["k_dp2_a"] * math.cos(2 * t - b) + k["k_dp2_b"] * (P * s1 + Q * c1),
        "p3": P,
        "q3": Q,
    }


def print_constants(m: Memoir, hyp: int) -> tuple[dict[str, float], float, float]:
    """§ 28 / § 45 constants as printed, in the keys of :func:`hkpq_constants`."""
    sec = 28 if hyp == 1 else 45
    g = lambda item: m.const(hyp, sec, item)  # noqa: E731
    k = {
        "k_dh1": g("k_dh1"),
        "k_dh2": g("k_dh2"),
        "k_dh3_a": g("k_dh3_a"),
        "k_dp1_a": g("k_dp1_a"),
        "k_dp1_b": g("k_dp1_b"),
        "k_dp2_a": g("k_dp2_a"),
        "k_dp2_b": g("k_dp2_b"),
    }
    return k, g("k_dh3_b"), g("k_dh3_c")


# ----------------------------------------------------------------------------- equations


def modern_rows(m: Memoir) -> list[tuple[float, float]]:
    """(t, c) for the 21 modern normal places, c = the § 9 mean-longitude residual."""
    return [(float(p["t_h1"]), p["oc_mean"]) for p in m.places if p["kind"] == "modern"]


def ancient_rows(m: Memoir, hyp: int) -> list[tuple[int, float, float]]:
    return [
        (p["year"], float(p[f"t_h{hyp}"]), p["oc_mean"]) for p in m.places if p["kind"] == "ancient"
    ]


def condition_equation(t: float, c: float, th: Theory) -> dict[str, float]:
    """§ 13: c = δε + t δn + Σ (cos, sin) terms with his printed arguments."""
    eq = {"lhs": c, "de": 1.0, "dn": t}
    for u in GROUP1[1:]:
        eq[u] = math.cos(math.radians(th.rate[u] * t))
    for u in GROUP2[1:]:
        eq[u] = math.sin(math.radians(th.rate[u] * t))
    return eq


def normal_equations(
    rows: list[dict[str, float]], names: tuple[str, ...], group: tuple[str, ...]
) -> dict[str, dict[str, float]]:
    """§ 14: multiply each equation by the coefficient of the named unknown and add. The two
    groups separate exactly for t symmetric about 0; only the group's terms are kept."""
    out = {}
    for name in names:
        var = PIVOT.get(name, name)
        eq = {k: 0.0 for k in ("lhs", *group)}
        for r in rows:
            w = r[var] / (
                10.0 if name == "n" else 1.0
            )  # he scales (n) by t/10 (Σ t² = 770 → 77.0000)
            for k in eq:
                eq[k] += w * r[k]
        out[name] = eq
    return out


def eliminate(eq: dict[str, float], pivot: dict[str, float], var: str) -> dict[str, float]:
    """eq − (eq[var]/pivot[var]) × pivot, the variable dropped."""
    f = eq.get(var, 0.0) / pivot[var]
    out = {k: v - f * pivot.get(k, 0.0) for k, v in eq.items() if k != var}
    for k, v in pivot.items():
        if k not in out and k != var:
            out[k] = -f * v
    return out


def combine(eqs: dict[str, dict[str, float]], mult: dict[str, float]) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, w in mult.items():
        for k, v in eqs[name].items():
            out[k] = out.get(k, 0.0) + w * v
    return out


def relations(e: float, beta_deg: float) -> dict[str, float]:
    """§ 22: δx₂ = r11 δx₁ + r12 δy₁, δy₂ = r21 δx₁ + r22 δy₁ with e/4 and β = ε − ϖ."""
    q = e / 4.0
    b = math.radians(beta_deg)
    return {
        "r11": q * math.cos(b),
        "r12": q * math.sin(b),
        "r21": -q * math.sin(b),
        "r22": q * math.cos(b),
    }


def substitute_second(eq: dict[str, float], R: dict[str, float]) -> dict[str, float]:
    """Replace δx₂, δy₂ by their § 22 expressions in δx₁, δy₁."""
    out = {k: v for k, v in eq.items() if k not in ("dx2", "dy2")}
    x2, y2 = eq.get("dx2", 0.0), eq.get("dy2", 0.0)
    out["dx1"] = out.get("dx1", 0.0) + x2 * R["r11"] + y2 * R["r21"]
    out["dy1"] = out.get("dy1", 0.0) + x2 * R["r12"] + y2 * R["r22"]
    return out


def eliminate_first(
    eq: dict[str, float], x_eq: dict[str, float], y_eq: dict[str, float], R: dict[str, float]
) -> tuple[dict[str, float], tuple[float, float]]:
    """§ 22 / § 26 / § 39 / § 43: after the § 22 substitution, add multiples of (x) and (y) that
    cancel δx₁ and δy₁; returns the equation and the two multipliers."""
    a, bx, by = (substitute_second(q, R) for q in (eq, x_eq, y_eq))
    M = np.array([[bx["dx1"], by["dx1"]], [bx["dy1"], by["dy1"]]])
    mx, my = np.linalg.solve(M, -np.array([a["dx1"], a["dy1"]]))
    out = combine({"a": a, "x": bx, "y": by}, {"a": 1.0, "x": float(mx), "y": float(my)})
    out.pop("dx1")
    out.pop("dy1")
    return out, (float(mx), float(my))


def zero_lhs(eqs: dict[str, dict[str, float]], first: str = "1") -> dict[str, dict[str, float]]:
    """§ 27 / § 44: eliminate the left-hand members of the others by equation (1). His form
    scales the equation with the smaller |lhs| up to the other: (2) has lhs 3.695 against (1)'s
    89.641, so e2 = (1) − (lhs₁/lhs₂)(2); (3) and (4) have |lhs| ≈ 480, so
    e3 = (3) − (lhs₃/lhs₁)(1).
    The roots do not depend on this; the printed coefficients do."""
    one = eqs[first]
    out = {}
    for name, eq in eqs.items():
        if name == first:
            continue
        if abs(eq["lhs"]) >= abs(one["lhs"]):
            f = eq["lhs"] / one["lhs"]
            out["e" + name] = {k: (0.0 if k == "lhs" else eq[k] - f * one[k]) for k in eq}
        else:
            f = one["lhs"] / eq["lhs"]
            out["e" + name] = {k: (0.0 if k == "lhs" else one[k] - f * eq[k]) for k in eq}
    return out


def log_equations(
    zeros: dict[str, dict[str, float]],
    k: dict[str, float],
    r_cos: float,
    r_sin: float,
    beta_deg: float,
) -> dict[str, dict[str, float]]:
    """§ 29 / § 46: substitute § 28 into the three zero-sided equations; the coefficients of the
    twelve functions of θ, P, Q (returned as plain numbers; the print gives their logarithms)."""
    b = math.radians(beta_deg)
    cb, sb = math.cos(b), math.sin(b)
    out = {}
    for name, eq in zeros.items():
        a = {t: eq.get(t, 0.0) for t in HKPQ}
        # p1 = A1 sin(θ−β) + B1(P c2 − Q s2); q1 = −A1 cos(θ−β) − B1(P s2 + Q c2)
        # sin(θ−β) = sinθ cosβ − cosθ sinβ ; cos(θ−β) = cosθ cosβ + sinθ sinβ
        # p2 = A2 sin(2θ−β) + B2(P c1 − Q s1); q2 = A2 cos(2θ−β) + B2(P s1 + Q c1)
        A1, B1, A2, B2 = k["k_dp1_a"], k["k_dp1_b"], k["k_dp2_a"], k["k_dp2_b"]
        L = {
            "lhs": 0.0,
            "sin1": a["h1"] * k["k_dh1"] + a["p1"] * A1 * cb - a["q1"] * A1 * sb,
            "cos1": a["k1"] * k["k_dh1"] - a["p1"] * A1 * sb - a["q1"] * A1 * cb,
            "sin2": a["h2"] * k["k_dh2"] + a["p2"] * A2 * cb + a["q2"] * A2 * sb,
            "cos2": a["k2"] * k["k_dh2"] - a["p2"] * A2 * sb + a["q2"] * A2 * cb,
            "sin3": a["h3"] * k["k_dh3_a"],
            "cos3": a["k3"] * k["k_dh3_a"],
            "P3": a["p3"] + a["h3"] * r_cos - a["k3"] * r_sin,
            "Q3": a["q3"] + a["h3"] * r_sin + a["k3"] * r_cos,
            "P1c-Q1s": a["p2"] * B2,
            "P1s+Q1c": a["q2"] * B2,
            "P2c-Q2s": a["p1"] * B1,
            "P2s+Q2c": -a["q1"] * B1,
        }
        out["L" + str(int(name[1:]) - 1)] = L  # e2, e3, e4 → L1, L2, L3
    return out


def eval_log_equation(L: dict[str, float], theta: float, P: float, Q: float) -> float:
    t = math.radians(theta)
    s1, c1, s2, c2, s3, c3 = (
        math.sin(t),
        math.cos(t),
        math.sin(2 * t),
        math.cos(2 * t),
        math.sin(3 * t),
        math.cos(3 * t),
    )
    return (
        L["sin1"] * s1
        + L["cos1"] * c1
        + L["sin2"] * s2
        + L["cos2"] * c2
        + L["sin3"] * s3
        + L["cos3"] * c3
        + L["P3"] * P
        + L["Q3"] * Q
        + L["P1c-Q1s"] * (P * c1 - Q * s1)
        + L["P1s+Q1c"] * (P * s1 + Q * c1)
        + L["P2c-Q2s"] * (P * c2 - Q * s2)
        + L["P2s+Q2c"] * (P * s2 + Q * c2)
    )


def _pq_given_theta(
    L2: dict[str, float], L3: dict[str, float], theta: float
) -> tuple[float, float]:
    """L2 = L3 = 0 are linear in P, Q at fixed θ."""
    rows = []
    rhs = []
    for L in (L2, L3):
        f0 = eval_log_equation(L, theta, 0.0, 0.0)
        fP = eval_log_equation(L, theta, 1.0, 0.0) - f0
        fQ = eval_log_equation(L, theta, 0.0, 1.0) - f0
        rows.append([fP, fQ])
        rhs.append(-f0)
    P, Q = np.linalg.solve(np.array(rows), np.array(rhs))
    return float(P), float(Q)


def solve_log_equations(
    Ls: dict[str, dict[str, float]], theta_grid_deg: np.ndarray | None = None
) -> list[dict[str, float]]:
    """§ 30: every root of L1(θ, P(θ), Q(θ)) = 0 on (−180°, 180°], by sign changes on a grid
    and bisection; each with its P, Q."""
    from scipy.optimize import brentq

    L1, L2, L3 = Ls["L1"], Ls["L2"], Ls["L3"]

    def f(th):
        P, Q = _pq_given_theta(L2, L3, th)
        return eval_log_equation(L1, th, P, Q)

    grid = np.linspace(-180.0, 180.0, 3601) if theta_grid_deg is None else theta_grid_deg
    vals = np.array([f(x) for x in grid])
    roots = []
    for i in range(len(grid) - 1):
        if vals[i] == 0.0 or vals[i] * vals[i + 1] < 0.0:
            th = brentq(f, grid[i], grid[i + 1], xtol=1e-10) if vals[i] != 0.0 else grid[i]
            P, Q = _pq_given_theta(L2, L3, th)
            roots.append({"theta": float(th), "P": P, "Q": Q})
    return roots


@dataclass
class Solution:
    """θ, P, Q and what § 30 / § 47 derives from them."""

    hyp: int
    theta: float
    P: float
    Q: float
    e_prime_units: float  # e′ in units of 1/20
    beta_prime: float  # ε − ϖ′, degrees
    m_prime: float  # units of 1/5000
    eps: float
    eps_prime_1810: float
    mean_lon_1846: float
    varpi_prime_1810: float
    varpi_prime_1846: float
    hkpq: dict[str, float] = field(default_factory=dict)

    @property
    def eccentricity(self) -> float:
        return self.e_prime_units * ECC_UNIT

    @property
    def mass(self) -> float:
        return self.m_prime * MASS_UNIT

    def as_dict(self) -> dict:
        d = {
            k: getattr(self, k)
            for k in (
                "hyp",
                "theta",
                "P",
                "Q",
                "e_prime_units",
                "beta_prime",
                "m_prime",
                "eps",
                "eps_prime_1810",
                "mean_lon_1846",
                "varpi_prime_1810",
                "varpi_prime_1846",
            )
        }
        d["eccentricity"] = self.eccentricity
        d["mass"] = self.mass
        d["hkpq"] = dict(self.hkpq)
        return d


def complete_solution(
    hyp: int,
    theta: float,
    P: float,
    Q: float,
    c: dict[str, float],
    eq1: dict[str, float],
    k: dict[str, float],
    r_cos: float,
    r_sin: float,
    beta: float,
    eps: float,
    n_prime_deg: float,
    precession_deg: float,
) -> Solution:
    """§ 30: e′, β′ from P, Q; m′ from equation (1); the longitudes by his arithmetic
    (ε′ = ε − θ; mean longitude at 1846.762 = ε′ + 12 n′ + precession; ϖ′ = ε − β′)."""
    A, B = c["c_p3"], c["c_p3e"]
    t3 = math.radians(3 * theta - beta)
    es = (P - A * math.sin(t3)) / B  # e′ sin(3θ − β′), B signed (−93.01)
    ec = (Q - A * math.cos(t3)) / B
    e_prime = math.hypot(es, ec)
    beta_prime = (3 * theta - math.degrees(math.atan2(es, ec))) % 360.0
    over_m = hkpq_over_m(theta, P, Q, beta, k, r_cos, r_sin)
    m_prime = eq1["lhs"] / sum(eq1.get(t, 0.0) * over_m[t] for t in HKPQ)
    eps_prime = (eps - theta) % 360.0
    return Solution(
        hyp,
        theta,
        P,
        Q,
        e_prime,
        beta_prime,
        m_prime,
        eps,
        eps_prime,
        (eps_prime + 12.0 * n_prime_deg + precession_deg) % 360.0,
        (eps - beta_prime) % 360.0,
        (eps - beta_prime + precession_deg) % 360.0,
        {t: m_prime * v for t, v in over_m.items()},
    )


# ----------------------------------------------------------------------------- the chains


@dataclass
class Chain:
    """Every intermediate of one hypothesis, keyed as the print is, for the B-stage checks and
    the C entry points. ``source`` says where the chain started (C0, C1, C2)."""

    hyp: int
    source: str
    steps: dict[tuple[int, str], dict[str, float]] = field(default_factory=dict)
    multipliers: dict[str, tuple[float, float]] = field(default_factory=dict)
    roots: list[dict[str, float]] = field(default_factory=list)
    solution: Solution | None = None

    def put(self, sec: int, name: str, eq: dict[str, float]) -> dict[str, float]:
        self.steps[(sec, name)] = eq
        return eq


def _setup(m: Memoir, hyp: int, k_source: str):
    th = theory_from_print(m, hyp)
    sec = {
        1: dict(
            c=12,
            a=13,
            nrm=(15, 18),
            e1=16,
            e2=19,
            x2=17,
            y2=20,
            anc=23,
            anc2=24,
            mul=25,
            fin=26,
            zero=27,
            hk=28,
            log=29,
            sol=30,
        ),
        2: dict(
            c=33,
            a=34,
            nrm=(35, 36),
            comb=37,
            e1=38,
            rel=39,
            anc=40,
            mul=41,
            anc2=42,
            fin=43,
            zero=44,
            hk=45,
            log=46,
            sol=47,
        ),
    }[hyp]
    e = m.const(1, 22, "e")
    beta = m.angle(1, 22, "eps_minus_varpi")
    if k_source == "print":
        k, r_cos, r_sin = print_constants(m, hyp)
        R = {
            "r11": m.const(1, 22, "rel_dx2_dx1"),
            "r12": m.const(1, 22, "rel_dx2_dy1"),
            "r21": m.const(1, 22, "rel_dy2_dx1"),
            "r22": m.const(1, 22, "rel_dy2_dy1"),
        }
    else:
        k = hkpq_constants(th.c)
        r = k.pop("k_dh3_r")
        r_cos, r_sin = r * math.cos(math.radians(beta)), r * math.sin(math.radians(beta))
        R = relations(e, beta)
    eps = ang(m.sol(hyp, sec["sol"], "eps"))
    n_prime = m.angle(hyp, sec["a"], "n_prime")
    prec = ang(m.sol(hyp, sec["sol"], "precession")) / 60.0  # printed "30" is 30′
    return th, sec, e, beta, k, r_cos, r_sin, R, eps, n_prime, prec


def run_hypothesis_1(m: Memoir, start: str = "C2", k_source: str = "print") -> Chain:
    """Hypothesis I from the data (C2), from his printed normal equations (C1) or from his
    printed final equations (C0); ``k_source`` chooses printed or derived § 28 constants."""
    th, S, _e, beta, k, r_cos, r_sin, R, eps, n_prime, prec = _setup(m, 1, k_source)
    ch = Chain(1, start)
    if start == "C0":
        Ls = m.eqs(1, 29)
        eq1 = m.eq(1, 22, "1")
    else:
        if start == "C2":
            rows = [condition_equation(t, c, th) for t, c in modern_rows(m)]
            N1 = normal_equations(rows, NORMAL_H1[1], GROUP1)
            N2 = normal_equations(rows, NORMAL_H1[2], GROUP2)
            anc = {f"anc_{y}": condition_equation(t, c, th) for y, t, c in ancient_rows(m, 1)}
        else:  # C1
            N1, N2 = m.eqs(1, 15), m.eqs(1, 18)
            anc = m.eqs(1, 23)
        for name, eq in list(N1.items()) + list(N2.items()):
            ch.put(S["nrm"][0] if name in N1 else S["nrm"][1], name, eq)
        for name, eq in anc.items():
            ch.put(S["anc"], name, eq)
        # §§ 16–17, 19–20
        E1 = {n: ch.put(16, n, eliminate(N1[n], N1["eps"], "de")) for n in NORMAL_H1[1][1:]}
        E2 = {n: ch.put(19, n, eliminate(N2[n], N2["n"], "dn")) for n in NORMAL_H1[2][1:]}
        X2 = {n: ch.put(17, n, eliminate(E1[n], E1["x"], "dx1")) for n in NORMAL_H1[1][2:]}
        Y2 = {n: ch.put(20, n, eliminate(E2[n], E2["y"], "dy1")) for n in NORMAL_H1[2][2:]}
        # § 21: "after giving the unknown quantities the same sign in them all": each equation
        # enters with the sign that makes its δx₂ (δy₂) coefficient negative
        sum1 = ch.put(21, "sum1", combine(X2, {n: -math.copysign(1.0, X2[n]["dx2"]) for n in X2}))
        sum2 = ch.put(21, "sum2", combine(Y2, {n: -math.copysign(1.0, Y2[n]["dy2"]) for n in Y2}))
        # § 22
        eq1, ch.multipliers["1"] = eliminate_first(sum1, E1["x"], E2["y"], R)
        eq2, ch.multipliers["2"] = eliminate_first(sum2, E1["x"], E2["y"], R)
        ch.put(22, "1", eq1)
        ch.put(22, "2", eq2)
        # § 24: δε by (ε), δx₁ by § 16's (x), δn by (n), δy₁ by § 19's (y)
        anc2 = {}
        for name, eq in anc.items():
            q = eliminate(eq, N1["eps"], "de")
            q = eliminate(q, E1["x"], "dx1")
            q = eliminate(q, N2["n"], "dn")
            q = eliminate(q, E2["y"], "dy1")
            anc2[name] = ch.put(24, name, q)
        # § 25
        mp = m.eq(1, 25, "mult_p3")
        mq = m.eq(1, 25, "mult_q3")
        comb_p = ch.put(25, "comb_p3", combine(anc2, mp))
        comb_q = ch.put(25, "comb_q3", combine(anc2, mq))
        # § 26
        eq3, ch.multipliers["3"] = eliminate_first(comb_p, E1["x"], E2["y"], R)
        eq4, ch.multipliers["4"] = eliminate_first(comb_q, E1["x"], E2["y"], R)
        ch.put(26, "3", eq3)
        ch.put(26, "4", eq4)
        zeros = zero_lhs({"1": eq1, "2": eq2, "3": eq3, "4": eq4})
        for n, q in zeros.items():
            ch.put(27, n, q)
        Ls = log_equations(zeros, k, r_cos, r_sin, beta)
    for n, q in Ls.items():
        ch.put(29, n, q)
    ch.roots = solve_log_equations(Ls)
    ch.solution = _pick_root(ch, m, 1, th, eq1, k, r_cos, r_sin, beta, eps, n_prime, prec)
    return ch


def run_hypothesis_2(m: Memoir, start: str = "C2", k_source: str = "print") -> Chain:
    """Hypothesis II by his shortened procedure (§§ 35–47)."""
    th, _S, _e, beta, k, r_cos, r_sin, R, eps, n_prime, prec = _setup(m, 2, k_source)
    ch = Chain(2, start)
    if start == "C0":
        Ls = m.eqs(2, 46)
        eq1 = m.eq(2, 39, "1")
    else:
        if start == "C2":
            rows = [condition_equation(t, c, th) for t, c in modern_rows(m)]
            N1 = normal_equations(rows, NORMAL_H2[1], GROUP1)
            N2 = normal_equations(rows, NORMAL_H2[2], GROUP2)
            anc = {f"anc_{y}": condition_equation(t, c, th) for y, t, c in ancient_rows(m, 2)}
        else:
            N1, N2 = m.eqs(2, 35), m.eqs(2, 36)
            anc = m.eqs(2, 40)
        for name, eq in list(N1.items()) + list(N2.items()):
            ch.put(35 if name in N1 else 36, name, eq)
        for name, eq in anc.items():
            ch.put(40, name, eq)
        # § 37: −(h₁) + (h₂) + (h₃); −(k₁) + (k₂) + (k₃)
        comb1 = ch.put(37, "comb1", combine(N1, {"h1": -1.0, "h2": 1.0, "h3": 1.0}))
        comb2 = ch.put(37, "comb2", combine(N2, {"k1": -1.0, "k2": 1.0, "k3": 1.0}))
        # § 38
        x = ch.put(38, "x", eliminate(N1["x"], N1["eps"], "de"))
        y = ch.put(38, "y", eliminate(N2["y"], N2["n"], "dn"))
        c1 = ch.put(38, "comb1", eliminate(comb1, N1["eps"], "de"))
        c2 = ch.put(38, "comb2", eliminate(comb2, N2["n"], "dn"))
        # § 39
        eq1, ch.multipliers["1"] = eliminate_first(c1, x, y, R)
        eq2, ch.multipliers["2"] = eliminate_first(c2, x, y, R)
        ch.put(39, "1", eq1)
        ch.put(39, "2", eq2)
        # § 41: the combinations first, then § 42 δε, δn out, then § 43 δx₁, δy₁ out
        mp = m.eq(2, 41, "mult_p3")
        mq = m.eq(2, 41, "mult_q3")
        comb_p = ch.put(41, "comb_p3", combine(anc, mp))
        comb_q = ch.put(41, "comb_q3", combine(anc, mq))
        cp = ch.put(42, "comb_p3", eliminate(eliminate(comb_p, N1["eps"], "de"), N2["n"], "dn"))
        cq = ch.put(42, "comb_q3", eliminate(eliminate(comb_q, N1["eps"], "de"), N2["n"], "dn"))
        eq3, ch.multipliers["3"] = eliminate_first(cp, x, y, R)
        eq4, ch.multipliers["4"] = eliminate_first(cq, x, y, R)
        ch.put(43, "3", eq3)
        ch.put(43, "4", eq4)
        zeros = zero_lhs({"1": eq1, "2": eq2, "3": eq3, "4": eq4})
        for n, q in zeros.items():
            ch.put(44, n, q)
        Ls = log_equations(zeros, k, r_cos, r_sin, beta)
    for n, q in Ls.items():
        ch.put(46, n, q)
    ch.roots = solve_log_equations(Ls)
    ch.solution = _pick_root(ch, m, 2, th, eq1, k, r_cos, r_sin, beta, eps, n_prime, prec)
    return ch


def _pick_root(
    ch: Chain, m: Memoir, hyp: int, th: Theory, eq1, k, r_cos, r_sin, beta, eps, n_prime, prec
) -> Solution | None:
    """Complete every root; keep the one with positive mass nearest his printed θ (all roots stay
    in ``ch.roots`` with their masses)."""
    if not ch.roots:
        return None
    printed_theta = ang(m.sol(hyp, 30 if hyp == 1 else 47, "theta"))
    best = None
    for r in ch.roots:
        s = complete_solution(
            hyp, r["theta"], r["P"], r["Q"], th.c, eq1, k, r_cos, r_sin, beta, eps, n_prime, prec
        )
        r["m_prime"] = s.m_prime
        r["e_prime_units"] = s.e_prime_units
        d = abs((r["theta"] - printed_theta + 180.0) % 360.0 - 180.0)
        if s.m_prime > 0 and (best is None or d < best[0]):
            best = (d, s)
    return best[1] if best else None


# ----------------------------------------------------------------------------- back-substitution


def element_corrections(
    ch: Chain, m: Memoir, sol: Solution, R: dict[str, float] | None = None
) -> dict[str, float]:
    """§ 50: δx₁, δy₁ from (x), (y) with the perturbation terms known; δx₂, δy₂ by § 22;
    δε, δn from (ε), (n); δa/a, δe, e δϖ from them."""
    hyp = sol.hyp
    e = m.const(1, 22, "e")
    beta = m.angle(1, 22, "eps_minus_varpi")
    R = R or relations(e, beta)
    sx = 16 if hyp == 1 else 38
    sy = 19 if hyp == 1 else 38
    x = substitute_second(ch.steps[(sx, "x")], R)
    y = substitute_second(ch.steps[(sy, "y")], R)
    known = sol.hkpq

    def residual(eq):
        return eq["lhs"] - sum(eq.get(t, 0.0) * known[t] for t in HKPQ)

    M = np.array([[x["dx1"], x["dy1"]], [y["dx1"], y["dy1"]]])
    dx1, dy1 = np.linalg.solve(M, np.array([residual(x), residual(y)]))
    dx2 = R["r11"] * dx1 + R["r12"] * dy1
    dy2 = R["r21"] * dx1 + R["r22"] * dy1
    se, sn = (15, 18) if hyp == 1 else (35, 36)
    E, N = ch.steps[(se, "eps")], ch.steps[(sn, "n")]
    vals = {"dx1": dx1, "dy1": dy1, "dx2": dx2, "dy2": dy2, **known}
    de = (
        E["lhs"] - sum(E.get(t, 0.0) * vals.get(t, 0.0) for t in E if t not in ("lhs", "de"))
    ) / E["de"]
    dn = (
        N["lhs"] - sum(N.get(t, 0.0) * vals.get(t, 0.0) for t in N if t not in ("lhs", "dn"))
    ) / N["dn"]
    th = theory_from_print(m, hyp)
    n_arcsec = th.n_deg * 3600.0
    b = math.radians(beta)
    return {
        "de": float(de),
        "dn": float(dn),
        "dx1": float(dx1),
        "dy1": float(dy1),
        "dx2": float(dx2),
        "dy2": float(dy2),
        "da_over_a": float(-2.0 / 3.0 * dn / n_arcsec),
        "de_elem": float((dx1 * math.sin(b) + dy1 * math.cos(b)) / 2.0),
        "e_dvarpi": float((-dx1 * math.cos(b) + dy1 * math.sin(b)) / 2.0),
    }


def residuals(m: Memoir, hyp: int, sol: Solution, corr: dict[str, float]) -> list[dict[str, float]]:
    """§§ 51–52: for each normal place the elements part, the perturbation part, their sum and
    O − theory, with his arguments."""
    th = theory_from_print(m, hyp)
    out = []
    for p in m.places:
        t = float(p[f"t_h{hyp}"])
        eq = condition_equation(t, p["oc_mean"], th)
        el = (
            corr["de"] + t * corr["dn"] + sum(eq[u] * corr[u] for u in ("dx1", "dx2", "dy1", "dy2"))
        )
        pert = sum(eq[u] * sol.hkpq[u] for u in HKPQ)
        out.append(
            {
                "year": p["year"],
                "t": t,
                "elements": el,
                "perturbation": pert,
                "total": el + pert,
                "residual": p["oc_mean"] - el - pert,
            }
        )
    return out


def true_longitude(mean_lon_deg: float, e: float, varpi_deg: float) -> float:
    """Heliocentric true longitude from the mean longitude, e and ϖ (Kepler's equation)."""
    M = math.radians((mean_lon_deg - varpi_deg) % 360.0)
    E = M
    for _ in range(50):
        E = E - (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
    nu = 2.0 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    return (math.degrees(nu) + varpi_deg) % 360.0


def printed_solution(m: Memoir, hyp: int) -> dict[str, float]:
    """§ 30–31 / § 47–48 as numbers."""
    s = 30 if hyp == 1 else 47
    el = 31 if hyp == 1 else 48
    return {
        "theta": ang(m.sol(hyp, s, "theta")),
        "P": _num(m.sol(hyp, s, "p3_over_m")),
        "Q": _num(m.sol(hyp, s, "q3_over_m")),
        "eps_prime_1810": ang(m.sol(hyp, s, "eps_prime_1810")),
        "mean_lon_1846": ang(m.sol(hyp, s, "mean_lon_1846")),
        "e_prime_units": _num(m.sol(hyp, s, "e_prime_units")),
        "beta_prime": ang(m.sol(hyp, s, "beta_prime")),
        "varpi_prime_1810": ang(m.sol(hyp, s, "varpi_prime_1810")),
        "varpi_prime_1846": ang(m.sol(hyp, s, "varpi_prime_1846")),
        "m_prime": _num(m.sol(hyp, s, "m_prime")),
        "eccentricity": _num(m.sol(hyp, el, "eccentricity")),
        "mass": _num(m.sol(hyp, el, "mass")),
    }


# ------------------------------------------------------------------------ stage B, step by step


def _printed_relations(m: Memoir) -> dict[str, float]:
    return {
        "r11": m.const(1, 22, "rel_dx2_dx1"),
        "r12": m.const(1, 22, "rel_dx2_dy1"),
        "r21": m.const(1, 22, "rel_dy2_dx1"),
        "r22": m.const(1, 22, "rel_dy2_dy1"),
    }


def reproduce_step(m: Memoir, hyp: int, sec: int) -> tuple[dict[str, dict[str, float]], dict]:
    """Section ``sec`` of hypothesis ``hyp`` computed from the *printed* sections it rests on
    (B5); returns the equations and any multipliers found. Unknown sections raise KeyError."""
    R = _printed_relations(m)
    beta = m.angle(1, 22, "eps_minus_varpi")
    P = lambda s: m.eqs(hyp, s)  # noqa: E731
    mult: dict = {}
    if hyp == 1:
        if sec == 16:
            N = P(15)
            return {n: eliminate(N[n], N["eps"], "de") for n in NORMAL_H1[1][1:]}, mult
        if sec == 17:
            E = P(16)
            return {n: eliminate(E[n], E["x"], "dx1") for n in NORMAL_H1[1][2:]}, mult
        if sec == 19:
            N = P(18)
            return {n: eliminate(N[n], N["n"], "dn") for n in NORMAL_H1[2][1:]}, mult
        if sec == 20:
            E = P(19)
            return {n: eliminate(E[n], E["y"], "dy1") for n in NORMAL_H1[2][2:]}, mult
        if sec == 21:
            X2, Y2 = P(17), P(20)
            return {
                "sum1": combine(X2, {n: -math.copysign(1.0, X2[n]["dx2"]) for n in X2}),
                "sum2": combine(Y2, {n: -math.copysign(1.0, Y2[n]["dy2"]) for n in Y2}),
            }, mult
        if sec == 22:
            S, x, y = P(21), P(16)["x"], P(19)["y"]
            e1, mult["1"] = eliminate_first(S["sum1"], x, y, R)
            e2, mult["2"] = eliminate_first(S["sum2"], x, y, R)
            return {"1": e1, "2": e2}, mult
        if sec == 24:
            N1, N2, E1, E2 = P(15), P(18), P(16), P(19)
            out = {}
            for n, eq in P(23).items():
                q = eliminate(eq, N1["eps"], "de")
                q = eliminate(q, E1["x"], "dx1")
                q = eliminate(q, N2["n"], "dn")
                out[n] = eliminate(q, E2["y"], "dy1")
            return out, mult
        if sec == 25:
            A2 = P(24)
            return {
                "comb_p3": combine(A2, m.eq(1, 25, "mult_p3")),
                "comb_q3": combine(A2, m.eq(1, 25, "mult_q3")),
            }, mult
        if sec == 26:
            C, x, y = P(25), P(16)["x"], P(19)["y"]
            e3, mult["3"] = eliminate_first(C["comb_p3"], x, y, R)
            e4, mult["4"] = eliminate_first(C["comb_q3"], x, y, R)
            return {"3": e3, "4": e4}, mult
        if sec == 27:
            return zero_lhs({**P(22), **P(26)}), mult
        if sec == 29:
            k, rc, rs = print_constants(m, 1)
            return log_equations(P(27), k, rc, rs, beta), mult
    else:
        if sec == 37:
            N1, N2 = P(35), P(36)
            return {
                "comb1": combine(N1, {"h1": -1.0, "h2": 1.0, "h3": 1.0}),
                "comb2": combine(N2, {"k1": -1.0, "k2": 1.0, "k3": 1.0}),
            }, mult
        if sec == 38:
            N1, N2, C = P(35), P(36), P(37)
            return {
                "x": eliminate(N1["x"], N1["eps"], "de"),
                "y": eliminate(N2["y"], N2["n"], "dn"),
                "comb1": eliminate(C["comb1"], N1["eps"], "de"),
                "comb2": eliminate(C["comb2"], N2["n"], "dn"),
            }, mult
        if sec == 39:
            E = P(38)
            e1, mult["1"] = eliminate_first(E["comb1"], E["x"], E["y"], R)
            e2, mult["2"] = eliminate_first(E["comb2"], E["x"], E["y"], R)
            return {"1": e1, "2": e2}, mult
        if sec == 41:
            A0 = P(40)
            return {
                "comb_p3": combine(A0, m.eq(2, 41, "mult_p3")),
                "comb_q3": combine(A0, m.eq(2, 41, "mult_q3")),
            }, mult
        if sec == 42:
            N1, N2, C = P(35), P(36), P(41)
            return {
                n: eliminate(eliminate(C[n], N1["eps"], "de"), N2["n"], "dn")
                for n in ("comb_p3", "comb_q3")
            }, mult
        if sec == 43:
            C, E = P(42), P(38)
            e3, mult["3"] = eliminate_first(C["comb_p3"], E["x"], E["y"], R)
            e4, mult["4"] = eliminate_first(C["comb_q3"], E["x"], E["y"], R)
            return {"3": e3, "4": e4}, mult
        if sec == 44:
            return zero_lhs(
                {**{n: P(39)[n] for n in ("1", "2")}, **{n: P(43)[n] for n in ("3", "4")}}
            ), mult
        if sec == 46:
            k, rc, rs = print_constants(m, 2)
            return log_equations(P(44), k, rc, rs, beta), mult
    raise KeyError((hyp, sec))


STEPS = {1: (16, 17, 19, 20, 21, 22, 24, 25, 26, 27, 29), 2: (37, 38, 39, 41, 42, 43, 44, 46)}
INPUT_NAMES = {"mult_p3", "mult_q3", "rel1", "rel2", "rel3", "rel4", "mult1", "mult2"}


def compare_step(m: Memoir, hyp: int, sec: int) -> dict:
    """Max |Δ| of coefficients and of right-hand sides between :func:`reproduce_step` and the
    print (up to the overall sign of an equation, which he chooses freely), per equation."""
    mine, mult = reproduce_step(m, hyp, sec)
    printed = m.eqs(hyp, sec)
    rows = {}
    for name, eq in printed.items():
        if name in INPUT_NAMES:
            continue
        me = mine[name]
        sgn = -1.0 if me.get("lhs", 0.0) * eq.get("lhs", 0.0) < 0 else 1.0
        terms = {t: sgn * me.get(t, 0.0) - v for t, v in eq.items()}
        d_lhs = abs(terms.pop("lhs", 0.0))
        worst = max(terms, key=lambda t: abs(terms[t])) if terms else None
        rows[name] = {
            "coef_max": abs(terms[worst]) if worst else 0.0,
            "coef_at": worst,
            "lhs": d_lhs,
            "sign_flipped": sgn < 0,
        }
    return {
        "rows": rows,
        "multipliers": mult,
        "coef_max": max(r["coef_max"] for r in rows.values()),
        "lhs_max": max(r["lhs"] for r in rows.values()),
    }


# ----------------------------------------------------------------------------- stage B tables


def _log10_printed(m: Memoir, hyp: int, sec: int, item: str) -> float:
    x = m.const(hyp, sec, item)
    return x - 10.0 if x >= 5.0 else x


def b1_table(m: Memoir) -> list[dict]:
    """B1: the twelve printed logarithms per hypothesis against quadrature."""
    out = []
    for hyp, alpha, sec in ((1, 0.5, 12), (2, 0.515, 33)):
        for i in range(4):
            b = laplace_b(i, alpha)
            for item, v in zip(("log_b", "log_adb", "log_a2d2b"), b, strict=True):
                printed = _log10_printed(m, hyp, sec, f"{item}{i}")
                out.append(
                    {
                        "hyp": hyp,
                        "item": f"{item}{i}",
                        "printed": printed,
                        "computed": math.log10(v),
                        "delta": math.log10(v) - printed,
                    }
                )
    return out


def printed_laplace(m: Memoir, hyp: int) -> dict[int, tuple[float, float, float]]:
    sec = 12 if hyp == 1 else 33
    return {
        i: tuple(
            10 ** _log10_printed(m, hyp, sec, f"{it}{i}")
            for it in ("log_b", "log_adb", "log_a2d2b")
        )
        for i in range(4)
    }


def b2_table(m: Memoir) -> list[dict]:
    """B2: the nine first-order series constants from § 11 with B1's Laplace values and with his
    printed ones, against the print."""
    out = []
    e = m.const(1, 22, "e")
    for hyp in (1, 2):
        th = theory_from_print(m, hyp)
        mine = series_constants(th.alpha, th.n_prime_over_n, e)
        his = series_constants(th.alpha, th.n_prime_over_n, e, printed_laplace(m, hyp))
        for k in mine:
            out.append(
                {
                    "hyp": hyp,
                    "item": k,
                    "printed": th.c[k],
                    "computed": mine[k],
                    "with_his_logs": his[k],
                    "delta": mine[k] - th.c[k],
                    "delta_his_logs": his[k] - th.c[k],
                }
            )
        kd = hkpq_constants(th.c)
        kp, rc, rs = print_constants(m, hyp)
        r = kd.pop("k_dh3_r")
        b = math.radians(m.angle(1, 22, "eps_minus_varpi"))
        for k in kp:
            out.append(
                {
                    "hyp": hyp,
                    "item": f"§28 {k}",
                    "printed": kp[k],
                    "computed": kd[k],
                    "delta": kd[k] - kp[k],
                }
            )
        out.append(
            {
                "hyp": hyp,
                "item": "§28 k_dh3_b",
                "printed": rc,
                "computed": r * math.cos(b),
                "delta": r * math.cos(b) - rc,
            }
        )
        out.append(
            {
                "hyp": hyp,
                "item": "§28 k_dh3_c",
                "printed": rs,
                "computed": r * math.sin(b),
                "delta": r * math.sin(b) - rs,
            }
        )
    return out


def exact_normal_equations(m: Memoir, hyp: int):
    th = theory_from_print(m, hyp)
    rows = [condition_equation(t, c, th) for t, c in modern_rows(m)]
    NH = NORMAL_H1 if hyp == 1 else NORMAL_H2
    return normal_equations(rows, NH[1], GROUP1), normal_equations(rows, NH[2], GROUP2)


def _diff_rows(
    mine: dict[str, dict[str, float]], printed: dict[str, dict[str, float]]
) -> list[dict]:
    out = []
    for name, eq in printed.items():
        if name in INPUT_NAMES:
            continue
        me = mine[name]
        for t, v in eq.items():
            out.append(
                {
                    "equation": name,
                    "term": t,
                    "printed": v,
                    "computed": me.get(t, 0.0),
                    "delta": me.get(t, 0.0) - v,
                }
            )
    return out


def b3_table(m: Memoir, hyp: int) -> list[dict]:
    """B3: every coefficient and right-hand side of his modern normal equations against the sums
    over the 21 equations of condition formed from § 9 and § 13."""
    N1, N2 = exact_normal_equations(m, hyp)
    s1, s2 = (15, 18) if hyp == 1 else (35, 36)
    return _diff_rows(N1, m.eqs(hyp, s1)) + _diff_rows(N2, m.eqs(hyp, s2))


def b4_table(m: Memoir, hyp: int) -> list[dict]:
    """B4: the nine ancient equations of condition from the printed epochs (the δn coefficient)."""
    th = theory_from_print(m, hyp)
    anc = {f"anc_{y}": condition_equation(t, c, th) for y, t, c in ancient_rows(m, hyp)}
    return _diff_rows(anc, m.eqs(hyp, 23 if hyp == 1 else 40))


def b5_table(m: Memoir, hyp: int) -> dict[int, dict]:
    """B5: each printed section from the printed sections it rests on (:func:`compare_step`); the
    logarithmic equations also in log units."""
    out = {}
    for sec in STEPS[hyp]:
        c = compare_step(m, hyp, sec)
        if sec in (29, 46):
            mine, _ = reproduce_step(m, hyp, sec)
            pr = m.eqs(hyp, sec)
            c["log10_max"] = max(
                abs(math.log10(abs(mine[n][t])) - math.log10(abs(pr[n][t])))
                for n in pr
                for t in pr[n]
                if t != "lhs"
            )
            c["sign_agreements"] = all(
                math.copysign(1, mine[n][t]) == math.copysign(1, pr[n][t])
                for n in pr
                for t in pr[n]
                if t != "lhs"
            )
        out[sec] = c
    return out


# ----------------------------------------------------------------------------- stage C tables

TIERS = {  # quantity: (tier 1, tier 2); angles in degrees, e in eccentricity, mass relative
    "theta": (1 / 60, 10 / 60),
    "mean_lon_1846": (2 / 60, 15 / 60),
    "eccentricity": (0.0005, 0.005),
    "varpi_prime_1846": (10 / 60, 1.0),
    "mass_rel": (0.003, 0.02),
}
BROOKES_1970 = {  # his re-derivation from the memoir's equations (Celestial Mechanics 3, 71)
    1: {
        "eccentricity": 0.1509,
        "m_prime": 0.8267,
        "varpi_prime_1846": 311 + 35 / 60,
        "mean_lon_1846": 323 + 39 / 60,
        "true_lon_discovery": 328 + 4 / 60,
    },
    2: {
        "eccentricity": 0.1387,
        "m_prime": 0.7393,
        "varpi_prime_1846": 303 + 52 / 60,
        "mean_lon_1846": 324 + 47 / 60,
        "true_lon_discovery": 331 + 31 / 60,
    },
}


def _wrap(d: float) -> float:
    return (d + 180.0) % 360.0 - 180.0


def tier_of(sol: Solution, target: dict[str, float]) -> dict:
    """Per-quantity differences from a target and the tier each lands in; the worst names the
    tier (1, 2 or 0 for not reproduced). Quantities the target lacks are skipped."""
    rows = {}
    s = sol.as_dict()
    for q, (t1, t2) in TIERS.items():
        if q == "mass_rel":
            if "m_prime" not in target:
                continue
            d = abs(s["m_prime"] / target["m_prime"] - 1.0)
        elif q not in target:
            continue
        elif q in ("theta", "mean_lon_1846", "varpi_prime_1846"):
            d = abs(_wrap(s[q] - target[q]))
        else:
            d = abs(s[q] - target[q])
        rows[q] = {"delta": d, "tier": 1 if d <= t1 else 2 if d <= t2 else 0}
    tiers = [r["tier"] for r in rows.values()]
    worst = 0 if (not tiers or 0 in tiers) else max(tiers)
    return {"rows": rows, "tier": worst}


def printed_target(m: Memoir, hyp: int) -> dict[str, float]:
    p = printed_solution(m, hyp)
    return {
        "theta": p["theta"],
        "mean_lon_1846": p["mean_lon_1846"],
        "eccentricity": p["eccentricity"],
        "varpi_prime_1846": p["varpi_prime_1846"],
        "m_prime": p["m_prime"],
    }


def printed_residuals(m: Memoir, hyp: int) -> list[dict]:
    return [
        {"year": p["year"], "residual": p[f"residual_{'h1' if hyp == 1 else 'h2'}"]}
        for p in m.places
        if p["year"] != 1690
    ]


def chain_with_print_fallback(m: Memoir, hyp: int, ch: Chain) -> Chain:
    """For C0 the chain has no normal equations; § 50 needs (ε), (x), (n), (y): take the print."""
    need = {
        1: ((15, "eps"), (16, "x"), (18, "n"), (19, "y")),
        2: ((35, "eps"), (38, "x"), (36, "n"), (38, "y")),
    }
    for sec, name in need[hyp]:
        if (sec, name) not in ch.steps:
            ch.steps[(sec, name)] = m.eq(hyp, sec, name)
    return ch


def c_table(m: Memoir) -> dict:
    """C0, C1, C2 for both hypotheses: the solution, its tier against the print and against
    Brookes, the § 50 corrections and the § 52 residuals against the print."""
    out = {}
    for hyp, run in ((1, run_hypothesis_1), (2, run_hypothesis_2)):
        pt = printed_target(m, hyp)
        for start in ("C0", "C1", "C2"):
            ch = chain_with_print_fallback(m, hyp, run(m, start))
            s = ch.solution
            corr = element_corrections(ch, m, s)
            res = residuals(m, hyp, s, corr)
            pr = {r["year"]: r["residual"] for r in printed_residuals(m, hyp)}
            dres = [abs(r["residual"] - pr[r["year"]]) for r in res if r["year"] in pr]
            rec = {
                "solution": s.as_dict(),
                "true_lon_1846": true_longitude(
                    s.mean_lon_1846, s.eccentricity, s.varpi_prime_1846
                ),
                "corrections": corr,
                "residuals": res,
                "residual_max_delta": max(dres),
                "residual_tier": 1 if max(dres) <= 0.05 else 2 if max(dres) <= 0.5 else 0,
                "vs_print": tier_of(s, pt),
                "vs_brookes": tier_of(s, BROOKES_1970[hyp]),
                "roots": ch.roots,
                "multipliers": ch.multipliers,
            }
            rec["tier"] = (
                min(rec["vs_print"]["tier"], rec["residual_tier"]) if rec["vs_print"]["tier"] else 0
            )
            out[f"H{hyp}/{start}"] = rec
    return out


# ----------------------------------------------------------------------------- locating the gap


def with_equations(
    m: Memoir, hyp: int, override: dict[tuple[int, str], dict[str, float]]
) -> Memoir:
    """A copy of the memoir with some printed equations replaced (for the sensitivity table)."""
    import copy

    mm = copy.deepcopy(m)
    for (sec, name), eq in override.items():
        mm.equations[(hyp, sec, name)] = {k: repr(v) for k, v in eq.items()}
    return mm


def sensitivity_table(m: Memoir, hyp: int, step: float = 0.0005) -> dict:
    """Which printed normal equation carries the difference between C1 and C2: each equation
    swapped for the exact one (whole, right-hand side only, coefficients only), and each
    coefficient of (n) moved by one unit of his fourth decimal."""
    run = run_hypothesis_1 if hyp == 1 else run_hypothesis_2
    s1, s2 = (15, 18) if hyp == 1 else (35, 36)
    N1, N2 = exact_normal_equations(m, hyp)
    P1, P2 = m.eqs(hyp, s1), m.eqs(hyp, s2)

    def sol(override):
        s = run(with_equations(m, hyp, override), "C1").solution
        return {
            "theta": s.theta,
            "mean_lon_1846": s.mean_lon_1846,
            "e": s.eccentricity,
            "m_prime": s.m_prime,
            "varpi_prime_1846": s.varpi_prime_1846,
        }

    base = sol({})
    out = {"base": base, "equations": {}, "n_coefficients": {}}
    out["group1_exact"] = sol({(s1, n): N1[n] for n in N1})
    out["group2_exact"] = sol({(s2, n): N2[n] for n in N2})
    for name in list(N1) + list(N2):
        sec, mine, printed = (s1, N1[name], P1[name]) if name in N1 else (s2, N2[name], P2[name])
        rhs_only = dict(printed)
        rhs_only["lhs"] = mine["lhs"]
        coef_only = dict(mine)
        coef_only["lhs"] = printed["lhs"]
        rec = {}
        for lab, eq in (("whole", mine), ("rhs", rhs_only), ("coefficients", coef_only)):
            r = sol({(sec, name): eq})
            rec[lab] = {
                k: _wrap(r[k] - base[k]) if k != "e" and k != "m_prime" else r[k] - base[k]
                for k in r
            }
        out["equations"][name] = rec
    for t in P2["n"]:
        if t in ("lhs", "dn"):
            continue
        eq = dict(P2["n"])
        eq[t] = eq[t] + step
        r = sol({(s2, "n"): eq})
        out["n_coefficients"][t] = {
            "printed": P2["n"][t],
            "exact": N2["n"][t],
            "dtheta_per_step": _wrap(r["theta"] - base["theta"]),
            "dL_per_step": _wrap(r["mean_lon_1846"] - base["mean_lon_1846"]),
        }
    return out


# ------------------------------------------------------- stage D1/D2: his model, solved directly

UNIT_YR = 3.0362  # his unit of time, years (the interval between modern normal places)
EPOCH_YEAR = 1810.328  # t = 0
REPORT_YEAR = 1846.762  # ε′ + 12 n′: the comparison epoch of the whole step
POLES = {"3:1": 3.0 ** (-2.0 / 3.0), "2:1": 2.0 ** (-2.0 / 3.0)}  # α where n′ = n/3, n/2
PARAMS = ("de", "dn", "dx1", "dy1", "theta", "m_prime", "e_prime_units", "beta_prime")


@dataclass(frozen=True)
class Setup:
    """Everything his model needs at one α: the theory (constants and argument rates), the § 28
    forms, the § 22 relations, Uranus's β, ε and the precession to 1846.762."""

    hyp: int
    th: Theory
    k: dict[str, float]
    r_cos: float
    r_sin: float
    R: dict[str, float]
    beta: float
    eps: float
    n_prime_deg: float  # per unit of time
    prec: float

    @property
    def alpha(self) -> float:
        return self.th.alpha


def setup_from_print(m: Memoir, hyp: int, k_source: str = "print") -> Setup:
    th, _S, _e, beta, k, r_cos, r_sin, R, eps, n_prime, prec = _setup(m, hyp, k_source)
    return Setup(hyp, th, k, r_cos, r_sin, R, beta, eps, n_prime, prec)


def rates_at(n_deg: float, n_prime_deg: float) -> dict[str, float]:
    """His arguments (§ 13) for any n′: δx, δy at n and 2n; h_i, k_i at i(n − n′); p_i, q_i at
    n′, n − 2n′, 2n − 3n′ (the e′ terms of the i-th order carry (i − 1)n − i n′)."""
    n, n1 = n_deg, n_prime_deg
    return {
        "dx1": n, "dy1": n, "dx2": 2 * n, "dy2": 2 * n,
        "h1": n - n1, "k1": n - n1, "h2": 2 * (n - n1), "k2": 2 * (n - n1),
        "h3": 3 * (n - n1), "k3": 3 * (n - n1),
        "p1": n1, "q1": n1, "p2": n - 2 * n1, "q2": n - 2 * n1,
        "p3": 2 * n - 3 * n1, "q3": 2 * n - 3 * n1,
    }  # fmt: skip


def second_order_at(m: Memoir, c: dict[str, float], drop: bool = False) -> dict[str, float]:
    """The two second-order constants of § 12 (e² and e e′ terms of the third argument: +0.57,
    −1.08 at α = 0.5; +0.40, −0.74 at 0.515), which § 11 does not give. They are carried as the
    printed Hypothesis I values scaled by c_p3e(α)/c_p3e(0.5), the small-divisor factor the
    printed pairs follow within 2 % (0.682 against 0.70 and 0.685); ``drop`` sets them to 0 for
    the sensitivity profile. Not derived: a choice, stated in the results."""
    if drop:
        return {"c_sin3_e2": 0.0, "c_sin3_e2e": 0.0}
    f = c["c_p3e"] / m.const(1, 12, "c_p3e")
    return {"c_sin3_e2": f * m.const(1, 12, "c_sin3_e2"),
            "c_sin3_e2e": f * m.const(1, 12, "c_sin3_e2e")}  # fmt: skip


def setup_at_alpha(m: Memoir, alpha: float, hyp: int = 1, drop_second_order: bool = False) -> Setup:
    """His model rebuilt at α: B1–B2 recomputed, n′ = n α^{3/2}, the § 28 forms derived, the
    § 22 relations from his e and β, the arguments from :func:`rates_at`. Rates stay in degrees
    per his unit of time with t in that unit, which is the same as degrees per year with t in
    years (D2's wording); the dates are the fixed ones either way."""
    e = m.const(1, 22, "e")
    beta = m.angle(1, 22, "eps_minus_varpi")
    n = m.angle(1, 13, "n")
    ratio = alpha**1.5
    n1 = n * ratio
    c = series_constants(alpha, ratio, e)
    c = c | second_order_at(m, c, drop=drop_second_order)
    k = hkpq_constants(c)
    r = k.pop("k_dh3_r")
    b = math.radians(beta)
    th = Theory(alpha, ratio, c, rates_at(n, n1), n)
    sec = 30 if hyp == 1 else 47
    eps = ang(m.sol(hyp, sec, "eps"))
    prec = ang(m.sol(hyp, sec, "precession")) / 60.0
    return Setup(hyp, th, k, r * math.cos(b), r * math.sin(b), relations(e, beta), beta, eps, n1,
                 prec)  # fmt: skip


def hkpq_from_elements(
    theta: float, m_prime: float, e_units: float, beta_prime: float, s: Setup
) -> tuple[dict[str, float], float, float]:
    """h, k, p, q from (θ, m′, e′, β′) through P = p₃/m′ and Q = q₃/m′ (the inverse of § 30)."""
    A, B = s.th.c["c_p3"], s.th.c["c_p3e"]
    t3 = math.radians(3 * theta - s.beta)
    t3p = math.radians(3 * theta - beta_prime)
    P = A * math.sin(t3) + B * e_units * math.sin(t3p)
    Q = A * math.cos(t3) + B * e_units * math.cos(t3p)
    over_m = hkpq_over_m(theta, P, Q, s.beta, s.k, s.r_cos, s.r_sin)
    return {u: m_prime * v for u, v in over_m.items()}, P, Q


def place_rows(m: Memoir, s: Setup, oc: dict[int, float] | None = None, drop_1690: bool = True):
    """(year, kind, t, equation) for his normal places with the arguments of ``s``; ``oc``
    replaces the § 9 residuals (synthetic data), keyed by year."""
    out = []
    for p in m.places:
        if drop_1690 and p["year"] == 1690:
            continue
        t = float(p[f"t_h{s.hyp}"])
        c = p["oc_mean"] if oc is None else oc[p["year"]]
        out.append((p["year"], p["kind"], t, condition_equation(t, c, s.th)))
    return out


def model_value(eq: dict[str, float], t: float, x: np.ndarray, s: Setup) -> float:
    de, dn, dx1, dy1, theta, mp, eu, bp = x
    dx2 = s.R["r11"] * dx1 + s.R["r12"] * dy1
    dy2 = s.R["r21"] * dx1 + s.R["r22"] * dy1
    hk, _, _ = hkpq_from_elements(theta, mp, eu, bp, s)
    v = de + t * dn + eq["dx1"] * dx1 + eq["dx2"] * dx2 + eq["dy1"] * dy1 + eq["dy2"] * dy2
    return v + sum(eq[u] * hk[u] for u in HKPQ)


def weights_for(m: Memoir, hyp: int, scheme: str, rows) -> np.ndarray:
    """(a) unit weights; (b) modern only (ancient weight 0); (c) 1/RMS² of his printed § 52
    residuals for the modern and the ancient places (conditioned on his own solution)."""
    kinds = [r[1] for r in rows]
    if scheme == "a":
        return np.ones(len(rows))
    if scheme == "b":
        return np.array([1.0 if k == "modern" else 0.0 for k in kinds])
    col = f"residual_h{hyp}"
    rms = {}
    for kind in ("modern", "ancient"):
        v = [p[col] for p in m.places if p["kind"] == kind and p[col] is not None]
        rms[kind] = math.sqrt(sum(r * r for r in v) / len(v))
    return np.array([1.0 / rms[k] ** 2 for k in kinds])


def residual_rms(m: Memoir, hyp: int) -> dict[str, float]:
    col = f"residual_h{hyp}"
    out = {}
    for kind in ("modern", "ancient"):
        v = [p[col] for p in m.places if p["kind"] == kind and p[col] is not None]
        out[kind] = math.sqrt(sum(r * r for r in v) / len(v))
    return out


def _starts(printed: dict[str, float] | None) -> list[np.ndarray]:
    starts = []
    if printed is not None:
        starts.append(np.array([0.0, 0.0, 0.0, 0.0, printed["theta"], printed["m_prime"],
                                printed["e_prime_units"], printed["beta_prime"]]))  # fmt: skip
    for th in range(0, 360, 30):
        for bp in range(0, 360, 90):
            starts.append(np.array([0.0, 0.0, 0.0, 0.0, float(th), 0.8, 2.0, float(bp)]))
    return starts


def solve_least_squares(
    rows, w: np.ndarray, s: Setup, printed: dict[str, float] | None = None
) -> dict:
    """His 29 equations and his model as one weighted nonlinear least-squares problem in the
    eight parameters, from a θ × β′ grid of starts (and his printed solution when given);
    m′ ≥ 0 and e′ ≤ 1 enforced. Returns the solution, the formal 1σ with σ² scaled to χ²/dof,
    χ², the RMS per group and every start's end point."""
    from scipy.optimize import least_squares

    keep = w > 0
    R = [r for r, k in zip(rows, keep, strict=True) if k]
    sw = np.sqrt(w[keep])
    lhs = np.array([r[3]["lhs"] for r in R])

    def fun(x):
        return sw * (lhs - np.array([model_value(r[3], r[2], x, s) for r in R]))

    lo = np.array([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 0.0, 0.0, -np.inf])
    hi = np.array([np.inf, np.inf, np.inf, np.inf, np.inf, 50.0, 20.0, np.inf])
    best, ends = None, []
    for x0 in _starts(printed):
        try:
            sol = least_squares(fun, x0, bounds=(lo, hi), method="trf", x_scale="jac",
                                xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=4000)  # fmt: skip
        except ValueError:
            continue
        chi2 = float(np.sum(sol.fun**2))
        ends.append({"start_theta": float(x0[4]), "start_beta_prime": float(x0[7]),
                     "theta": float(sol.x[4] % 360.0), "chi2": chi2,
                     "m_prime": float(sol.x[5])})  # fmt: skip
        if best is None or chi2 < best[0] - 1e-9:
            best = (chi2, sol)
    chi2, sol = best
    x = sol.x.copy()
    x[4] %= 360.0
    x[7] %= 360.0
    dof = len(R) - len(x)
    J = sol.jac
    try:
        cov = np.linalg.inv(J.T @ J) * chi2 / dof
        sig = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    except np.linalg.LinAlgError:
        sig = np.full(len(x), float("nan"))
    hk, P, Q = hkpq_from_elements(x[4], x[5], x[6], x[7], s)
    eps_prime = (s.eps - x[4]) % 360.0
    solution = Solution(s.hyp, float(x[4]), P, Q, float(x[6]), float(x[7]), float(x[5]), s.eps,
                        eps_prime, (eps_prime + 12.0 * s.n_prime_deg + s.prec) % 360.0,
                        (s.eps - x[7]) % 360.0, (s.eps - x[7] + s.prec) % 360.0, hk)  # fmt: skip
    resid = lhs - np.array([model_value(r[3], r[2], x, s) for r in R])
    rms = {}
    for kind in ("modern", "ancient"):
        v = [float(d) for d, r in zip(resid, R, strict=True) if r[1] == kind]
        rms[kind] = math.sqrt(sum(q * q for q in v) / len(v)) if v else None
    return {
        "x": dict(zip(PARAMS, (float(v) for v in x), strict=True)),
        "sigma": dict(zip(PARAMS, (float(v) for v in sig), strict=True)),
        "solution": solution,
        "chi2": chi2,
        "dof": dof,
        "n_rows": len(R),
        "rms_arcsec": rms,
        "residuals": [
            {"year": r[0], "kind": r[1], "t": r[2], "residual": float(d)}
            for r, d in zip(R, resid, strict=True)
        ],
        "true_lon_1846": true_longitude(
            solution.mean_lon_1846, solution.eccentricity, solution.varpi_prime_1846
        ),
        "starts": ends,
        "n_distinct_minima": len({round(e["chi2"], 3) for e in ends}),
    }


def chi2_of(rows, w: np.ndarray, s: Setup, x: np.ndarray) -> float:
    keep = w > 0
    return float(sum(wi * (r[3]["lhs"] - model_value(r[3], r[2], x, s)) ** 2
                     for r, wi, k in zip(rows, w, keep, strict=True) if k))  # fmt: skip


def x_of(sol: Solution, corr: dict[str, float]) -> np.ndarray:
    """The eight parameters of a chain solution with its § 50 corrections."""
    return np.array([corr["de"], corr["dn"], corr["dx1"], corr["dy1"], sol.theta, sol.m_prime,
                     sol.e_prime_units, sol.beta_prime])  # fmt: skip


def alpha_grid(lo: float = 0.40, hi: float = 0.60, step: float = 0.01) -> list[float]:
    n = round((hi - lo) / step)
    return [round(lo + i * step, 4) for i in range(n + 1)]


def alpha_flags(alpha: float) -> list[str]:
    """D2's bands: truncation-limited at α ≥ 0.58 (the plan); pole-limited within 0.015 of the
    3:1 commensurability at α = 0.4807, where the i = 3 eccentric terms' divisor n − 3n′
    vanishes (found when the grid ran, reported, not pre-registered)."""
    flags = []
    if alpha >= 0.58 - 1e-9:
        flags.append("truncation_limited")
    if abs(alpha - POLES["3:1"]) <= 0.015:
        flags.append("pole_limited_3_1")
    return flags


# ----------------------------------------------------------------------- step 3c: noise draws
# docs/notes/P1B_adams_noise.md § 2: per-place noise scales and the one-start solve.

N_PLACES_USED = 29  # his 30 places less 1690, which he set aside
N_PARAMS = len(PARAMS)
SCALE_A = math.sqrt(N_PLACES_USED / (N_PLACES_USED - N_PARAMS))  # √(N/(N − p)) = 1.175
GREENWICH_LAST_YEAR = 1831  # modern places to here rest on Greenwich in set D; after, Paris


def set_d_group_of_place(place: dict) -> str:
    """The set D group whose observations a normal place's window holds (model (b)); an empty
    modern window takes Greenwich's single-measurement scatter (the plan's floor)."""
    if place["kind"] == "ancient":
        return "ancient"
    return "greenwich" if place["year"] <= GREENWICH_LAST_YEAR else "paris_1835_1845"


def noise_sigmas(
    m: Memoir,
    model: str,
    counts: dict[int, int] | None = None,
    hyp: int = 1,
    group_rms: dict[str, float] | None = None,
    scale_a: float = SCALE_A,
) -> dict[int, float]:
    """σ_j in arcsec per place (keyed by year; 1690 excluded). (a) his § 52 residual RMS per
    kind × √(N/(N − p)); (b) and (c) the set D group RMS over √max(N_j, 1) with N_j the set D
    observations in the place's window (step 3c § 2)."""
    out = {}
    if model == "a":
        rms = residual_rms(m, hyp)
        for p in m.places:
            if p["year"] != 1690:
                out[p["year"]] = rms[p["kind"]] * scale_a
        return out
    if model not in ("b", "c"):
        raise ValueError(model)
    if counts is None or group_rms is None:
        raise ValueError("models (b) and (c) need the window counts and the group RMS")
    for p in m.places:
        if p["year"] == 1690:
            continue
        n = max(int(counts.get(p["year"], 0)), 1)
        out[p["year"]] = group_rms[set_d_group_of_place(p)] / math.sqrt(n)
    return out


_LS_LO = np.array([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 0.0, 0.0, -np.inf])
_LS_HI = np.array([np.inf, np.inf, np.inf, np.inf, np.inf, 50.0, 20.0, np.inf])


def _ls_problem(rows, w: np.ndarray, s: Setup):
    keep = w > 0
    R = [r for r, k in zip(rows, keep, strict=True) if k]
    sw = np.sqrt(w[keep])
    lhs = np.array([r[3]["lhs"] for r in R])

    def fun(x):
        return sw * (lhs - np.array([model_value(r[3], r[2], x, s) for r in R]))

    return R, fun


def solve_warm(rows, w: np.ndarray, s: Setup, x0: np.ndarray, max_nfev: int = 2000) -> dict:
    """The eight-parameter weighted least squares from one start (step 3c's primary method:
    what an iterated solution does). Returns x (θ, β′ wrapped), χ², ``converged`` (the solver
    stopped on a tolerance, not on ``max_nfev``), ``bound`` (m′ or e′ on a bound), the unscaled
    formal σ (√diag (JᵀWJ)⁻¹, the noise-propagation prediction before any χ²/dof factor) and
    the number of function evaluations."""
    from scipy.optimize import least_squares

    _R, fun = _ls_problem(rows, w, s)
    x0 = np.asarray(x0, dtype=float).copy()
    x0[5] = min(max(x0[5], _LS_LO[5]), _LS_HI[5])
    x0[6] = min(max(x0[6], _LS_LO[6]), _LS_HI[6])
    try:
        sol = least_squares(fun, x0, bounds=(_LS_LO, _LS_HI), method="trf", x_scale="jac",
                            xtol=1e-10, ftol=1e-10, gtol=1e-10, max_nfev=max_nfev)  # fmt: skip
    except ValueError:
        return {"x": np.full(N_PARAMS, np.nan), "chi2": float("nan"), "converged": False,
                "bound": False, "sigma_unscaled": np.full(N_PARAMS, np.nan), "nfev": 0}  # fmt: skip
    x = sol.x.copy()
    x[4] %= 360.0
    x[7] %= 360.0
    bound = bool(x[5] <= _LS_LO[5] + 1e-9 or x[6] >= _LS_HI[6] - 1e-9 or x[5] >= _LS_HI[5] - 1e-9)
    try:
        sig = np.sqrt(np.clip(np.diag(np.linalg.inv(sol.jac.T @ sol.jac)), 0.0, None))
    except np.linalg.LinAlgError:
        sig = np.full(N_PARAMS, np.nan)
    return {
        "x": x,
        "chi2": float(np.sum(sol.fun**2)),
        "converged": bool(sol.status > 0 and np.all(np.isfinite(x))),
        "bound": bound,
        "sigma_unscaled": sig,
        "nfev": int(sol.nfev),
    }


def theta_grid_starts(x_clean: np.ndarray, step_deg: float = 30.0) -> list[np.ndarray]:
    """The multistart check's grid: θ every ``step_deg``, m′ 0.8, e′ 0.1 (2.0 units), β′ the
    clean value, the linear corrections zero (step 3c § 2)."""
    starts = []
    for th in np.arange(0.0, 360.0, step_deg):
        starts.append(np.array([0.0, 0.0, 0.0, 0.0, float(th), 0.8, 2.0, float(x_clean[7])]))
    return starts


def solve_multistart(
    rows, w: np.ndarray, s: Setup, x_clean: np.ndarray, max_nfev: int = 2000
) -> dict:
    """The lowest-χ² end point of ``theta_grid_starts`` (12 solves); ``n_capped`` counts the
    starts that stopped on ``max_nfev`` (skipped)."""
    best = None
    n_capped = 0
    for x0 in theta_grid_starts(x_clean):
        r = solve_warm(rows, w, s, x0, max_nfev=max_nfev)
        if not r["converged"] or not math.isfinite(r["chi2"]):
            n_capped += int(r["nfev"] >= max_nfev)
            continue
        if best is None or r["chi2"] < best["chi2"] - 1e-9:
            best = r
    if best is None:
        best = {"x": np.full(N_PARAMS, np.nan), "chi2": float("nan"), "converged": False,
                "bound": False, "sigma_unscaled": np.full(N_PARAMS, np.nan), "nfev": 0}  # fmt: skip
    return {**best, "n_capped": n_capped}
