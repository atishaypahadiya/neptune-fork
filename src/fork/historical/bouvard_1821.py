"""Bouvard's 1821 theory of Uranus: the printed formula and a table path for a few instants.

Two routes to the heliocentric longitude of Uranus, both in the book's centesimal division
(the circle in 400 degrees, a degree in 100 minutes, a minute in 100 seconds), reckoned from
the mean equinox of date with Bouvard's uniform precession of 154″.63 a year:

- ``formula``: the closed form of the introduction (pp. ii, iii, xv). The mean longitudes φ,
  φ′, φ″ of p. xv carry the great inequalities A, A′, A″ of pp. ii–iii.
- ``table_longitude``: the tables as a computer of the time would use them (pp. xviii–xx),
  from the rows transcribed for the validation of docs/notes/P1B_bouvard_attribution.md. Only
  the instants of that note can be computed this way.

Time is t in Julian years from the Paris mean midnight that begins 1 January 1800.
"""

from __future__ import annotations

import csv
import datetime
import math
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from fork.data.centesimal import CENTESIMAL_SEC_IN_ARCSEC

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "historical" / "adams_challis"
BOUVARD = DATA / "bouvard_1821"

CIRCLE_CS = 4_000_000.0  # centesimal seconds in the circle
ARGS = ("II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII")
ARG_PARTS = {a: (10_000 if a == "III" else 1_000) for a in ARGS}
# Table IX, head of each block: the sign of the correction of each argument
ARG_CORRECTION_SIGN = {
    "II": -1, "III": -1, "IV": 1, "V": 1, "VI": -1, "VII": -1,
    "VIII": 1, "IX": 1, "X": -1, "XI": -1, "XII": -1,
}  # fmt: skip
# argument -> (table of the transcription, column, step of the printed argument)
EQUATION_TABLE = {
    "II": ("XI", "eq", 10), "III": ("XII", "eq", 100), "IV": ("XIII", "eq", 10),
    "V": ("XIV", "eq", 10), "VI": ("XV-XVII", "VI", 10), "VII": ("XV-XVII", "VII", 10),
    "VIII": ("XV-XVII", "VIII", 10), "IX": ("XVIII-XXI", "IX", 10),
    "X": ("XVIII-XXI", "X", 10), "XI": ("XVIII-XXI", "XI", 10), "XII": ("XVIII-XXI", "XII", 10),
}  # fmt: skip
HALF_UNIT_CS = 0.05  # every longitude table used is printed to 0″.1 centesimal
EPOCH_1800 = datetime.date(1800, 1, 1).toordinal()


@dataclass(frozen=True)
class Instant:
    """A Paris mean midnight (the beginning of the civil date) and the table rows it uses."""

    name: str
    date: tuple[int, int, int]
    table_year: int  # the row of Table I
    century: str | None  # the row of Table II, for a year before 1780
    month: str | None  # the row of Table III; None on 1 January
    day: str | None  # the row of Table IV; None on the first of a month

    @property
    def t(self) -> float:
        return (datetime.date(*self.date).toordinal() - EPOCH_1800) / 365.25

    @property
    def year(self) -> float:
        y = self.date[0]
        return y + (datetime.date(*self.date).toordinal() - datetime.date(y, 1, 1).toordinal()) / (
            366 if y % 4 == 0 and (y % 100 != 0 or y % 400 == 0) else 365
        )


# The Paris midnight nearest each of his normal epochs of 1783, 1798, 1813, 1828, nearest
# 1846.762, and (S4) nearest his 1715 epoch; fixed from adams_memoir_normal_places.csv.
INSTANTS = (
    Instant("1783", (1783, 1, 1), 1783, None, None, None),
    Instant("1798", (1798, 3, 7), 1798, None, "Mars common", "7"),
    Instant("1813", (1813, 5, 13), 1813, None, "Mai common", "13"),
    Instant("1828", (1828, 7, 18), 1828, None, "Juillet leap", "18"),
    Instant("1846", (1846, 10, 5), 1846, None, "Octobre common", "5"),
    Instant("1715", (1715, 3, 4), 1815, "-100", "Mars common", "4"),
)
# V1′: the midnights nearest seven of his normal epochs not used to find the corrections
INSTANTS_V1P = (
    Instant("1789", (1789, 1, 26), 1789, None, None, "26"),
    Instant("1804", (1804, 4, 3), 1804, None, "Avril leap", "3"),
    Instant("1822", (1822, 6, 22), 1822, None, "Juin common", "22"),
    Instant("1837", (1837, 8, 27), 1837, None, "Août common", "27"),
    Instant("1690", (1690, 12, 20), 1790, "-100", "Décembre common", "20"),
    Instant("1750", (1750, 11, 5), 1850, "-100", "Novembre common", "5"),
    Instant("1769", (1769, 1, 12), 1869, "-100", None, "12"),
)


def centesimal_seconds(s: str) -> float:
    """A printed centesimal angle in centesimal seconds: '111.58.62,6' (d.m.s,t), '81.02.48'
    (d.m.s), '1.21,3' (m.s,t), '26,3' (s,t)."""
    p = s.replace(",", ".").split(".")
    if len(p) == 4:
        return int(p[0]) * 10_000 + int(p[1]) * 100 + int(p[2]) + int(p[3]) / 10
    if len(p) == 3 and len(p[2]) == 2:
        return int(p[0]) * 10_000 + int(p[1]) * 100 + int(p[2])
    if len(p) == 3:
        return int(p[0]) * 100 + int(p[1]) + int(p[2]) / 10
    if len(p) == 2:
        return int(p[0]) + int(p[1]) / 10
    raise ValueError(s)


def _num(s: str) -> float:
    return float(s.replace(",", ".").replace("−", "-"))


@cache
def formula_cells() -> dict[int, str]:
    with (BOUVARD / "bouvard_1821_uranus_formula.csv").open(encoding="utf-8") as f:
        return {int(r["n"]): r["value_as_printed"] for r in csv.DictReader(f)}


@cache
def great_inequality_cells() -> dict[tuple[str, str], float | str]:
    out: dict[tuple[str, str], float | str] = {}
    with (BOUVARD / "bouvard_1821_great_inequalities.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = (r["formula"].split(",")[0], r["cell"])
            v = r["value_as_printed"]
            out[key] = v if r["cell"] == "lead_sign" else _num(v)
    return out


@cache
def table_rows(used: bool = True) -> dict[tuple[str, str, str], str]:
    with (BOUVARD / "bouvard_1821_uranus_table_rows.csv").open(encoding="utf-8") as f:
        col = "value_used" if used else "value_as_printed"
        return {(r["table"], r["row"], r["column"]): r[col] for r in csv.DictReader(f)}


def _sin(grades: float) -> float:
    return math.sin(grades * math.pi / 200.0)


def formula(
    t: float,
    centre_coefficient: float | None = None,
    sign_group_iii: int | None = None,
    great_inequalities: bool = True,
    mean_slip: float = 0.0,
    precession: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> dict[str, float]:
    """The printed formula at t. Angles in centesimal seconds unless named ``_grades``.

    ``centre_coefficient`` replaces the leading coefficient of the equation of the centre
    (printed 59427″.54); ``sign_group_iii`` the sign before group III (printed on worn type).
    ``mean_slip`` (centesimal seconds) is taken from the mean longitude that enters the
    longitude and the anomaly, not from the arguments; ``precession`` is the inequality of the
    precession added to the longitude, the perihelion and the node.
    """
    c = formula_cells()
    g = great_inequality_cells()
    s = 1e-4  # centesimal seconds -> degrees of 400
    mj = _num(c[4]) + _num(c[5]) * t
    ms = _num(c[6]) + _num(c[7]) * t
    mu = _num(c[8]) + _num(c[9]) * t
    sg = {"+": 1.0, "-": -1.0}
    a1 = (
        5 * ms
        - 2 * mj
        + g["A", "arg_const"]
        + g["A", "arg_t"] * s * t
        + g["A", "arg_t2"] * s * t * t
    )
    big_a = sg[g["A", "lead_sign"]] * (
        g["A", "coef"] + g["A", "coef_t"] * t + g["A", "coef_t2"] * t * t
    ) * _sin(a1) + g["A", "second_coef"] * _sin(2 * a1)
    a2 = (
        5 * ms - 2 * mj + g["A′", "arg_const"] + g["A′", "arg_t"] * s * t
        + g["A′", "arg_t2"] * s * t * t
    )  # fmt: skip
    big_a1 = (
        sg[g["A′", "lead_sign"]]
        * (g["A′", "coef"] + g["A′", "coef_t"] * t + g["A′", "coef_t2"] * t * t)
        * _sin(a2)
        + (g["A′", "second_coef"] + g["A′", "second_coef_t"] * t) * _sin(2 * a2)
        + g["A′", "third_coef"] * _sin(3 * mu - ms + g["A′", "third_arg_const"])
    )
    big_a2 = (
        sg[g["A″", "lead_sign"]]
        * (g["A″", "coef"] + g["A″", "coef_t"] * t)
        * _sin(3 * mu - ms + g["A″", "arg_const"] + g["A″", "arg_t"] * s * t)
    )
    if not great_inequalities:
        big_a = big_a1 = big_a2 = 0.0
    f, f1, f2 = mj + big_a * s, ms + big_a1 * s, mu + big_a2 * s
    prec = _num(c[27]) * s * t
    peri = _num(c[20]) + _num(c[21]) * s * t
    node = _num(c[22]) + _num(c[23]) * s * t
    # the perihelion of p. ii moves with the equinox (Table I prints the same motion), so the
    # anomaly is the longitude of date less the perihelion
    own = f2 + prec + (precession[0] - mean_slip) * s  # the mean longitude of date
    anomaly = own - peri - precession[1] * s
    lead = _num(c[28]) if centre_coefficient is None else centre_coefficient
    centre = (lead - _num(c[29]) * t) * _sin(anomaly)
    for k, n in enumerate((30, 32, 34, 36), 2):
        centre += (_num(c[n]) - _num(c[n + 1]) * t) * _sin(k * anomaly)
    x = f1 - f2
    ii = _num(c[38]) * _sin(x + _num(c[39])) + sum(
        _num(c[n]) * _sin(k * x) for k, n in enumerate((40, 41, 42, 43), 2)
    )
    s3 = sg[c[44]] if sign_group_iii is None else float(sign_group_iii)
    iii = s3 * (_num(c[45]) - _num(c[46]) * t) * _sin(
        f1 - 2 * f2 + _num(c[47]) + _num(c[48]) * s * t
    ) + _num(c[49]) * _sin(2 * f1 - 4 * f2 + _num(c[50]))
    iv = _num(c[51]) * _sin(f - f2) + _num(c[52]) * _sin(2 * (f - f2))
    arg = {
        "II": f1 - f2, "III": f1 - 2 * f2, "IV": f - f2, "V": f - 2 * f2,
        "VI": 2 * f1 - 3 * f2, "VII": f1, "VIII": 2 * f - f2, "IX": f,
        "X": 2 * f1 - 5 * f2, "XI": 2 * f1 - f2, "XII": 3 * f1 - 4 * f2,
    }  # fmt: skip
    small = sum(
        _num(c[n]) * _sin(arg[a] + _num(c[n + 1]))
        for a, n in zip(ARGS[3:], (53, 55, 57, 59, 61, 63, 65, 67), strict=True)
    )
    in_orbit = own / s + centre + ii + iii + iv + small
    reduction = _num(c[69]) * _sin(2 * (in_orbit * s - node - precession[2] * s))
    return {
        "longitude": (in_orbit + reduction) % CIRCLE_CS,
        "in_orbit": in_orbit % CIRCLE_CS,
        "mean_longitude": ((mu + prec) / s - mean_slip) % CIRCLE_CS,
        "great_inequality": big_a2,
        "anomaly_grades": anomaly % 400.0,
        "centre": centre,
        "perturbations": ii + iii + iv + small,
        "reduction": reduction,
        **{f"arg_{a}": (arg[a] % 400.0) / 400.0 * ARG_PARTS[a] for a in ARGS},
    }


def _lagrange(xs: list[float], ys: list[float], x: float) -> float:
    total = 0.0
    for i, (xi, yi) in enumerate(zip(xs, ys, strict=True)):
        p = yi
        for j, xj in enumerate(xs):
            if j != i:
                p *= (x - xj) / (xi - xj)
        total += p
    return total


CENTRE_FROM_TABLES = 59417.54  # the leading coefficient Table X answers to
MEAN_SLIP_PER_YEAR = 0.105  # centesimal seconds: Table I against the mean motion of p. ii
MEAN_SLIP_FROM = 1780.0


def precession_inequality(year: float, used: bool = True) -> tuple[float, float, float]:
    """The inequality of the precession for the longitude, the perihelion and the node
    (Jupiter's Table III, p. 12), four rows about the date."""
    d = table_rows(used)
    rows = [1600, 1700, 1750, 1800] if year < 1700 else [1700, 1750, 1800, 1900]
    return tuple(
        _lagrange(rows, [centesimal_seconds(d["P", str(r), col]) for r in rows], year)
        for col in ("lon", "peri", "node")
    )


def corrected_formula(inst: Instant) -> dict[str, float]:
    """The formula with the three corrections of the amendment of 13 November 2026 to
    docs/notes/P1B_bouvard_attribution.md: Table X's leading coefficient, Table I's mean
    motion (for a date before 1780 the row of Table I is a century later), and the
    inequality of the precession as the tables apply it."""
    row_year = inst.year + (100.0 if inst.century is not None else 0.0)
    return formula(
        inst.t,
        centre_coefficient=CENTRE_FROM_TABLES,
        mean_slip=MEAN_SLIP_PER_YEAR * (row_year - MEAN_SLIP_FROM),
        precession=precession_inequality(inst.year),
    )


JD_1800_PARIS_MIDNIGHT = 2378496.5 - (9 * 60 + 21) / 86400.0  # Paris is 9 m 21 s east


def t_of_jd(jd: float) -> float:
    """Julian years from the Paris mean midnight that begins 1 January 1800."""
    return (jd - JD_1800_PARIS_MIDNIGHT) / 365.25


def corrected_formula_at(t: float) -> dict[str, float]:
    """``corrected_formula`` at any time t (Julian years from 1800)."""
    year = 1800.0 + t
    row_year = year + (100.0 if year < MEAN_SLIP_FROM else 0.0)
    return formula(
        t,
        centre_coefficient=CENTRE_FROM_TABLES,
        mean_slip=MEAN_SLIP_PER_YEAR * (row_year - MEAN_SLIP_FROM),
        precession=precession_inequality(year),
    )


def mean_longitudes_deg(t: float, great_inequalities: bool = False) -> tuple[float, float, float]:
    """The mean longitudes of Jupiter, Saturn and Uranus of p. ii in degrees of 360, from the
    fixed equinox of 1800; with the great inequalities of pp. ii and iii if asked."""
    c = formula_cells()
    out = [_num(c[4]) + _num(c[5]) * t, _num(c[6]) + _num(c[7]) * t, _num(c[8]) + _num(c[9]) * t]
    if great_inequalities:
        with_a = formula(t)
        without = formula(t, great_inequalities=False)
        # arguments VII, IX carry the longitudes of Saturn and Jupiter; Uranus by A″ itself
        out[0] += (with_a["arg_IX"] - without["arg_IX"]) * 400.0 / ARG_PARTS["IX"]
        out[1] += (with_a["arg_VII"] - without["arg_VII"]) * 400.0 / ARG_PARTS["VII"]
        out[2] += with_a["great_inequality"] * 1e-4
    return tuple(x * 0.9 for x in out)


def s7_correction(t: float, great_inequalities: bool = False, terms=None) -> float:
    """Adams's § 7: the sum of his five terms in arcseconds, the arguments built from the
    mean longitudes of p. ii."""
    lon = mean_longitudes_deg(t, great_inequalities)
    total = 0.0
    for k, (coef, mult, const) in enumerate(s7_terms()):
        if terms is not None and k not in terms:
            continue
        arg = sum(m * x for m, x in zip(mult, lon, strict=True)) + const
        total += coef * math.sin(math.radians(arg))
    return total


def c_bouvard_deg(t: float) -> float:
    """C_Bouvard: the corrected formula's heliocentric longitude of Uranus (mean equinox of
    date, degrees of 360) plus the five terms of § 7."""
    return (corrected_formula_at(t)["longitude"] * 0.9e-4 + s7_correction(t) / 3600.0) % 360.0


def misprint_row_2600_shift(t: float) -> float:
    """S5: the change of the tables' longitude (arcseconds) when Table XII's row 2600 is used
    as printed (3.41,6) instead of 3.31,6, by the tables' four-row interpolation at the
    argument III of t."""
    arg = corrected_formula_at(t)["arg_III"]
    lo = math.floor(arg / 100) * 100
    grid = [lo - 100, lo, lo + 100, lo + 200]
    if 2600 not in grid:
        return 0.0
    weight = _lagrange(grid, [1.0 if g == 2600 else 0.0 for g in grid], arg)
    return 10.0 * weight * CENTESIMAL_SEC_IN_ARCSEC


def _centre_row(d: dict, row: int) -> tuple[float, float]:
    """Table X: the equation (the constant of p. xix subtracted, as printed) and the secular
    variation at a whole degree. Rows print the degrees in full (39x, a negative equation) or
    their last digit only."""
    p = d["X", str(row), "eq"].replace(",", ".").split(".")
    deg = int(p[0])
    if deg >= 390:
        deg -= 400
    elif row > 200:
        deg -= 10
    return (
        deg * 10_000 + int(p[1]) * 100 + int(p[2]) + int(p[3]) / 10,
        centesimal_seconds(d["X", str(row), "sv"]),
    )


def centre_constant() -> float:
    """The sum of the constants added to the equations of the longitude, read off Table X at
    anomaly 400 (= 0), where the equation of the centre itself is zero."""
    return -_centre_row(table_rows(), 400)[0]


def table_longitude(inst: Instant, used: bool = True) -> dict:
    """The longitude of Uranus at ``inst`` from the transcribed table rows."""
    d = table_rows(used)
    out: dict = {"entries": 0}

    def summed(col: str, conv) -> float:
        v = conv(d["I", str(inst.table_year), col])
        for table, row in (("II", inst.century), ("III", inst.month), ("IV", inst.day)):
            if row is not None:
                v += conv(d[table, row, col])
        return v

    lon = summed("lon", centesimal_seconds)
    peri = summed("peri", centesimal_seconds)
    node = summed("node", centesimal_seconds)
    n_epoch = 1 + sum(r is not None for r in (inst.century, inst.month, inst.day))
    y = inst.year
    dec = int(y // 10) * 10
    yrs = [dec - 10, dec, dec + 10, dec + 20]
    great = _lagrange(yrs, [centesimal_seconds(d["IX", str(r), "eq"]) for r in yrs], y)
    prec, prec_peri, prec_node = precession_inequality(y, used)
    corrected = lon + great + prec
    anomaly = ((corrected - peri - prec_peri) % CIRCLE_CS) / 10_000
    k = math.floor(anomaly)
    rows = [k - 1, k, k + 1, k + 2]
    vals = [_centre_row(d, r if r > 0 else r + 400) for r in rows]
    centre = _lagrange(rows, [v[0] for v in vals], anomaly)
    centre_linear = _lagrange(rows[1:3], [v[0] for v in vals[1:3]], anomaly)
    sv = _lagrange(rows[1:3], [v[1] for v in vals[1:3]], anomaly)
    # p. xix: the secular variation is additive after 1800; its printed sign is − in the
    # first half of the circle of the anomaly and + in the second
    secular = (-1.0 if anomaly < 200 else 1.0) * sv * (y - 1800.0) / 100.0
    args, eqs, slope = {}, {}, {}
    for a in ARGS:
        raw = summed(a, int)
        c0, c1 = int(d["IX", str(dec), a]), int(d["IX", str(dec + 10), a])
        args[a] = (raw + ARG_CORRECTION_SIGN[a] * (c0 + (c1 - c0) * (y - dec) / 10)) % ARG_PARTS[a]
        table, col, step = EQUATION_TABLE[a]
        lo = math.floor(args[a] / step) * step
        grid = [lo - step, lo, lo + step, lo + 2 * step] if a == "III" else [lo, lo + step]
        ys = [centesimal_seconds(d[table, str(int(r)), col]) for r in grid]
        eqs[a] = _lagrange(grid, ys, args[a])
        pair = ys[1:3] if a == "III" else ys
        slope[a] = abs(pair[1] - pair[0]) / step
    in_orbit = corrected + centre + secular + sum(eqs.values())
    u = ((in_orbit - node - prec_node) % CIRCLE_CS) / 10_000
    folded = u % 200.0
    r = folded if folded < 100 else 200.0 - folded
    kk = min(math.floor(r), 98)
    red = (-1.0 if folded < 100 else 1.0) * _lagrange(
        [kk, kk + 1],
        [centesimal_seconds(d["XXXII", str(kk + j), "red"]) for j in (0, 1)],
        r,
    )
    # entries summed into the longitude: epoch (with century, month, day), inequality of the
    # precession, great inequality, equation of the centre, its secular variation, eleven
    # equations, reduction
    entries = n_epoch + 1 + 1 + 1 + 1 + len(ARGS) + 1
    # the arguments are whole numbers of parts: each entry summed into one can be half a part
    # off, and the equation follows at its local slope
    n_arg = n_epoch + 1
    allowance = sum(slope[a] * 0.5 * n_arg for a in ARGS)
    out.update(
        longitude=(in_orbit + red) % CIRCLE_CS,
        in_orbit=in_orbit % CIRCLE_CS,
        mean_longitude=lon % CIRCLE_CS,
        precession_inequality=prec,
        great_inequality=great,
        anomaly_grades=anomaly,
        centre=centre + centre_constant(),
        centre_linear=centre_linear + centre_constant(),
        secular=secular,
        perturbations=sum(eqs.values()) - centre_constant(),
        reduction=red,
        args=args,
        entries=entries,
        budget_cs=entries * HALF_UNIT_CS,
        budget_arcsec=entries * HALF_UNIT_CS * CENTESIMAL_SEC_IN_ARCSEC,
        argument_allowance_cs=allowance,
    )
    return out


@cache
def s7_terms() -> tuple[tuple[float, tuple[int, int, int], float], ...]:
    """Adams's § 7: (coefficient in arcseconds, multipliers of the mean longitudes of
    (Jupiter, Saturn, Uranus), constant in degrees) for his five terms."""
    cells: dict[str, dict[str, str]] = {}
    with (DATA / "adams_memoir_s7_terms.csv").open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            cells.setdefault(r["term"], {})[r["cell"]] = r["value_as_printed"]
    out = []
    for k in sorted(cells):
        mult = [0, 0, 0]
        expr = cells[k]["arg"].replace("-", "+-")
        for piece in filter(None, expr.split("+")):
            body = piece.lstrip("-")
            digits = "".join(ch for ch in body[: body.index("p")])
            which = {"phi": 0, "phi1": 1, "phi2": 2}[body[body.index("p") :]]
            mult[which] = (-1 if piece.startswith("-") else 1) * int(digits or 1)
        const = 0.0
        if cells[k]["const"]:
            sign = -1.0 if cells[k]["const"].startswith("-") else 1.0
            deg, minutes = cells[k]["const"].lstrip("+-").split()
            const = sign * (int(deg) + float(minutes) / 60.0)
        out.append((float(cells[k]["coef"]), (mult[0], mult[1], mult[2]), const))
    return tuple(out)
