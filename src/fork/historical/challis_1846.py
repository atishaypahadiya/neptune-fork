"""Challis's search of 1846 as a detection model (docs/notes/P1B_challis_detection.md, § 4–5).

A *draw* is one setting of the judgement inputs and of what the printed record does not give
(zone settings, session dates and lengths). ``accept`` is the abduction: the draw is kept if
his actual procedure, run with fresh event draws, reproduces the record. A *case* turns a draw
into independent *groups* of events; ``cdfs`` gives the probability of recognition by each
night in closed form and ``sample`` draws the same events, as a check of the closed form.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import date

import numpy as np

NIGHT0 = "1846-07-29"
N_NIGHTS = 65  # 29 July … 1 October 1846


def night_index(iso: str) -> int:
    return (date.fromisoformat(iso) - date.fromisoformat(NIGHT0)).days


JUL29, JUL30, AUG4, AUG12 = (night_index(f"1846-{d}") for d in ("07-29", "07-30", "08-04", "08-12"))
SEP29, SEP30 = night_index("1846-09-29"), night_index("1846-09-30")
PRIMARY_LAST = night_index("1846-09-22")  # recognised before 23 September, Galle's night
SECONDARY_LAST = SEP30  # by 1 October: the nights and days up to 30 September
UNKNOWN_FIRST, UNKNOWN_LAST = night_index("1846-08-13"), night_index("1846-09-28")
KNOWN = (JUL29, JUL30, AUG4, AUG12, SEP29)
NEVER = 10_000

ZONE_ARCMIN = 9.0  # breadth of a zone at power 166 (Memoirs 16, 417–421)
SECTOR_HALF_ARCMIN = 35.0  # the 4 August sector sweep: "through 70′" (MNRAS 7, 147)
SEP29_HOURS, SEP29_STARS = 3.5, 300.0  # printed (Memoirs 16, 423)
PLANET_NUMBER_AUG12 = 49.0  # the planet was No. 49 of the 12 August series
STARS_TOTAL = (3000.0, 4000.0)  # abduction condition 4
MAX_UNKNOWN = 25

PRIORS = {  # § 4's table: uniform on these intervals
    "p_disc_routine": (0.0, 0.5),
    "p_disc_looking": (0.1, 1.0),
    "p_follow": (0.5, 1.0),
    "p_recorded": (0.9, 1.0),
    "stars_per_hour": (40.0, 90.0),
    "hours_per_session": (1.0, 3.5),
    "moving_extent_factor": (0.5, 2.0),
    "weather": (0.3, 0.6),
    "desk_hours_per_pair": (1.0, 6.0),
    "desk_hours_per_day": (0.5, 3.0),
    "n_unknown_sessions": (8, 25),
    # not fixed by the plan (pre-run record): the breadth in declination of C4b's sweeps
    "c4b_breadth_arcmin": (30.0, 120.0),
}


@dataclass
class Draw:
    p_r: float
    p_l: float
    p_f: float
    p_rec: float
    rate: float
    f_mov: float
    weather: float
    desk_pair: float
    desk_day: float
    n_unknown: int
    branch: str  # "dec": the planet entered the 12 August zone by its declination edge; "ra"
    north_edge_deg: float  # declination of the north edge of the 30 July / 12 August zone
    east_end_h: float  # branch "ra": right ascension of the 30 July zone's east end
    hours: dict[int, float]  # hours of each session, by night
    unknown: list[int]  # nights of the sessions between 12 August and 28 September
    usable: np.ndarray  # nights on which he could observe
    north_first: bool  # C3, C4: the first zone off centre lies north
    u_zone: np.ndarray  # where in a series the planet falls, per comparison
    c4b_breadth: float
    extra: dict = field(default_factory=dict)


def _u(rng, name: str) -> float:
    lo, hi = PRIORS[name]
    return float(rng.uniform(lo, hi))


def draw(rng: np.random.Generator, geo: dict) -> Draw:
    """One draw of § 4's priors. ``geo``: per-night arrays ``ra_h``, ``dec_deg`` of Neptune,
    ``moon_excluded`` and the targets' ``ra_h``/``dec_deg``."""
    dec30, dec12 = geo["dec_deg"][JUL30], geo["dec_deg"][AUG12]
    half = ZONE_ARCMIN / 60.0
    if rng.random() < 0.5:
        branch = "dec"  # north edge between the planet's two declinations
        north = float(rng.uniform(dec12, dec30))
        east = float("nan")
    else:
        branch = "ra"  # the zone holds both declinations; its east end lay between the places
        north = float(rng.uniform(dec30, dec12 + half))
        east = float(rng.uniform(geo["ra_h"][AUG12], geo["ra_h"][JUL30]))
    weather = _u(rng, "weather")
    usable = (~geo["moon_excluded"]) & (rng.random(N_NIGHTS) < weather)
    usable[list(KNOWN)] = True
    usable[AUG4 + 1 : AUG12] = False  # moonlight and weather, 4 to 12 August (primary)
    usable[SEP30] = False  # "the moon being in the way" (primary)
    lo, hi = PRIORS["n_unknown_sessions"]
    n = int(rng.integers(lo, hi + 1))
    cand = np.flatnonzero(usable[UNKNOWN_FIRST : UNKNOWN_LAST + 1]) + UNKNOWN_FIRST
    unknown = sorted(rng.choice(cand, size=min(n, cand.size), replace=False).tolist())
    hours = {k: _u(rng, "hours_per_session") for k in (*KNOWN[:4], *unknown)}
    hours[SEP29] = SEP29_HOURS
    return Draw(
        p_r=_u(rng, "p_disc_routine"), p_l=_u(rng, "p_disc_looking"), p_f=_u(rng, "p_follow"),
        p_rec=_u(rng, "p_recorded"), rate=_u(rng, "stars_per_hour"),
        f_mov=_u(rng, "moving_extent_factor"), weather=weather,
        desk_pair=_u(rng, "desk_hours_per_pair"), desk_day=_u(rng, "desk_hours_per_day"),
        n_unknown=n, branch=branch, north_edge_deg=north, east_end_h=east, hours=hours,
        unknown=unknown, usable=usable, north_first=bool(rng.random() < 0.5),
        u_zone=rng.random(MAX_UNKNOWN + 4), c4b_breadth=_u(rng, "c4b_breadth_arcmin"),
    )  # fmt: skip


def total_stars(d: Draw) -> float:
    return SEP29_STARS + d.rate * sum(h for k, h in d.hours.items() if k != SEP29)


def accept(d: Draw, rng: np.random.Generator) -> dict:
    """The abduction (§ 4): conditions 2 and 3 by fresh event draws of his actual procedure,
    condition 4 on the total of stars. Condition 1 holds by the zone priors."""
    no_disc_august = bool(rng.random() >= d.p_r) and bool(rng.random() >= d.p_r)
    disc_sep29 = bool(rng.random() < d.p_l)
    c4 = STARS_TOTAL[0] <= total_stars(d) <= STARS_TOTAL[1]
    return {"c23": no_disc_august and disc_sep29, "c4": c4,
            "all": no_disc_august and disc_sep29 and c4}  # fmt: skip


# ------------------------------------------------------------------ groups of events


@dataclass
class Group:
    """Independent block of events: each recording happens with probability ``q``; a recorded
    planet has its disc noted that night with probability ``pd``; if anything was recorded, a
    comparison flags it with probability ``pc`` on day ``tc``."""

    rec: list[tuple[int, float, float]]  # (night, q, pd)
    pc: float = 0.0
    tc: int = NEVER


def next_usable(d: Draw, day: int, strictly_after: bool) -> int:
    start = day + 1 if strictly_after else day
    for k in range(max(start, 0), N_NIGHTS):
        if d.usable[k]:
            return k
    return NEVER


def desk_days(d: Draw, fraction: float) -> int:
    """Days of desk work to reach ``fraction`` of one pair of zones."""
    return max(1, math.ceil(fraction * d.desk_pair / d.desk_day))


def _aug12_fraction(d: Draw) -> float:
    return min(1.0, PLANET_NUMBER_AUG12 / (d.rate * d.hours[AUG12]))


def _history(d: Draw) -> list[Group]:
    """What the record fixes: recorded on 4 and 12 August (routine) and 29 September (looking)."""
    return [Group([(AUG4, 1.0, d.p_r)]), Group([(AUG12, 1.0, d.p_r)]),
            Group([(SEP29, 1.0, d.p_l)])]  # fmt: skip


def _zone_covers(d: Draw, geo: dict, night_of_place: int) -> bool:
    """Did the 30 July / 12 August zone cover Neptune's place of ``night_of_place``?"""
    dec, ra = geo["dec_deg"][night_of_place], geo["ra_h"][night_of_place]
    in_dec = d.north_edge_deg - ZONE_ARCMIN / 60.0 <= dec <= d.north_edge_deg
    in_ra = True if d.branch == "dec" else ra <= d.east_end_h
    return bool(in_dec and in_ra)


def _sweeps(d: Draw) -> list[int]:
    return [JUL29, JUL30, AUG12, *d.unknown]


def _layout_groups(
    d: Draw, geo: dict, target: str, prompt: bool, delay: int,
    breadth_arcmin: float = ZONE_ARCMIN, looking: bool = False,
) -> list[Group]:  # fmt: skip
    """C3, C4: zones laid out from ``target``'s place, declination intervals outward from the
    centre, each swept twice (fixed, then moving) before the next; the 4 August sector sweep
    and the 29 September sweep as in the record."""
    tra, tdec = geo["targets"][target]["ra_h"], geo["targets"][target]["dec_deg"]
    b = breadth_arcmin / 60.0
    sweeps = _sweeps(d)
    per_zone = 1 if looking else 2
    out = []
    for z in range(math.ceil(len(sweeps) / per_zone)):
        nights = sweeps[z * per_zone : (z + 1) * per_zone]
        step = (z + 1) // 2
        sign = 0 if z == 0 else (1 if (z % 2 == 1) == d.north_first else -1)
        centre = tdec[nights[0]] + sign * step * b
        rec = []
        for j, n in enumerate(nights):
            length = d.hours[n] * (d.f_mov if j == 1 else 1.0) * ZONE_ARCMIN / breadth_arcmin
            inside = (abs(geo["dec_deg"][n] - centre) <= b / 2.0
                      and abs(geo["ra_h"][n] - tra[n]) <= length / 2.0)  # fmt: skip
            rec.append((n, d.p_rec if inside else 0.0, d.p_l if looking else d.p_r))
        g = Group(rec)
        if prompt and len(nights) == 2:
            g.pc, g.tc = d.p_f, nights[1] + delay + desk_days(d, float(d.u_zone[z])) - 1
        out.append(g)
    inside = (abs(geo["dec_deg"][AUG4] - tdec[AUG4]) <= SECTOR_HALF_ARCMIN / 60.0
              and abs(geo["ra_h"][AUG4] - tra[AUG4]) <= d.hours[AUG4] / 2.0)  # fmt: skip
    out.append(Group([(AUG4, d.p_rec if inside else 0.0, d.p_r)]))
    out.append(Group([(SEP29, 1.0, d.p_l)]))
    return out


CASES = ("A0", "C1a", "C1b", "C1a_7d", "C1b_7d", "C2_1800", "C2_date", "C3_S4", "C3_S4_prompt",
         "C3_S5", "C3_S5_prompt", "C4a", "C4a_prompt", "C4b")  # fmt: skip


def groups(d: Draw, geo: dict, case: str) -> list[Group]:
    if case == "A0":
        return _history(d)
    if case.startswith("C1"):
        delay = 7 if case.endswith("_7d") else 1
        g4, g12, g29 = _history(d)
        if case.startswith("C1a"):  # 12 August against 30 July, carried through No. 49
            if _zone_covers(d, geo, AUG12):
                g12.pc = d.p_f
                g12.tc = AUG12 + delay + desk_days(d, _aug12_fraction(d)) - 1
        else:  # the 4 August sector stars against the 12 August zone
            if _zone_covers(d, geo, AUG4):
                g4.pc = d.p_f
                g4.tc = AUG12 + delay + desk_days(d, float(d.u_zone[0])) - 1
        return [g4, g12, g29]
    if case.startswith("C2"):
        key = "on_chart_1800" if case == "C2_1800" else "on_chart_date"
        g4, g12, g29 = _history(d)
        for g, n, u in ((g4, AUG4, d.u_zone[0]), (g12, AUG12, _aug12_fraction(d))):
            if geo[key][n]:
                g.pc, g.tc = d.p_f, n + 1 + desk_days(d, float(u)) - 1
        return [g4, g12, g29]
    if case.startswith("C3"):
        return _layout_groups(d, geo, case.split("_")[1], case.endswith("prompt"), 1)
    if case.startswith("C4a"):
        return _layout_groups(d, geo, "LVJ", case.endswith("prompt"), 1)
    if case == "C4b":
        return _layout_groups(d, geo, "LVJ", False, 1, d.c4b_breadth, looking=True)
    raise ValueError(case)


def _no_flag_by(g: Group, when: np.ndarray, dates_rec: list[int], date_c: int) -> np.ndarray:
    """P(no event of the group dated at or before each entry of ``when``)."""
    k = len(g.rec)
    out = np.zeros(when.shape)
    for mask in range(1 << k):
        p = 1.0
        term = np.ones(when.shape)
        for i, (_n, q, pd) in enumerate(g.rec):
            if mask >> i & 1:
                p *= q
                term = term * (1.0 - pd * (dates_rec[i] <= when))
            else:
                p *= 1.0 - q
        if p == 0.0:
            continue
        if mask:
            term = term * (1.0 - g.pc * (date_c <= when))
        out += p * term
    return out


def cdfs(d: Draw, gs: list[Group]) -> tuple[np.ndarray, np.ndarray]:
    """P(T1 by each night) and P(T2 by each night), closed form. T2 is the first usable night
    after a disc noted at the eyepiece, or the first usable night from the day a comparison
    flags the planet."""
    when = np.arange(N_NIGHTS)
    no1, no2 = np.ones(N_NIGHTS), np.ones(N_NIGHTS)
    for g in gs:
        no1 *= _no_flag_by(g, when, [n for n, _, _ in g.rec], g.tc)
        no2 *= _no_flag_by(g, when, [next_usable(d, n, True) for n, _, _ in g.rec],
                           next_usable(d, g.tc, False) if g.tc < NEVER else NEVER)  # fmt: skip
    return 1.0 - no1, 1.0 - no2


def sample(d: Draw, gs: list[Group], rng: np.random.Generator, n: int) -> np.ndarray:
    """``n`` event draws of the groups: the day of T1 (``NEVER`` if none)."""
    t1 = np.full(n, NEVER)
    for g in gs:
        any_rec = np.zeros(n, dtype=bool)
        for night, q, pd in g.rec:
            rec = rng.random(n) < q
            any_rec |= rec
            disc = rec & (rng.random(n) < pd)
            t1 = np.where(disc, np.minimum(t1, night), t1)
        if g.pc > 0.0:
            flag = any_rec & (rng.random(n) < g.pc)
            t1 = np.where(flag, np.minimum(t1, g.tc), t1)
    return t1
