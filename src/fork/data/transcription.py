"""Double-entry transcription of Le Verrier (1846), Table n° 77: diff and self-checks.

The two passes (``data/historical/leverrier_1846_pass{1,2}.csv``) are written independently
from two scans, exactly as printed (see ``docs/notes/M6.md``). This module only compares them
and checks each row's internal arithmetic. It does not interpret or convert anything: that is
the loader's job, which comes after the passes are reconciled.

    uv run python -m fork.data.transcription            # diff pass 1 vs pass 2
"""

from __future__ import annotations

import csv
import difflib
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HIST = ROOT / "data" / "historical"
PASS_FILES = (HIST / "leverrier_1846_pass1.csv", HIST / "leverrier_1846_pass2.csv")
COLUMNS = ("page", "line", "year", "month", "day", "mean_time", "ra_obs", "ra_calc_sec",
           "ra_c_minus_o", "dec_obs", "dec_calc_sec", "dec_c_minus_o", "lon_c_minus_o",
           "lat_c_minus_o", "notes")  # fmt: skip
DATA_FIELDS = COLUMNS[2:-1]  # compared between passes (notes are free text)
NO_VALUE = "»"


def read_pass(path: Path) -> list[dict[str, str]]:
    """Rows of one pass, as strings, header checked."""
    with path.open(newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        if tuple(rd.fieldnames or ()) != COLUMNS:
            raise ValueError(f"{path.name}: header {rd.fieldnames} != {COLUMNS}")
        return [{k: (v or "").strip() for k, v in r.items()} for r in rd]


@dataclass(frozen=True)
class Disagreement:
    page: str
    line1: str | None  # None: the row is missing from pass 1
    line2: str | None
    field: str  # a column name, or "row" for an unmatched row
    pass1: str
    pass2: str


def _key(r: dict[str, str]) -> str:
    """What identifies a printed row across passes: its day and mean time."""
    return f"{r['day']}|{r['mean_time']}"


def diff_passes(a: list[dict[str, str]], b: list[dict[str, str]]) -> list[Disagreement]:
    """Field-by-field differences, page by page. Rows are aligned by (day, mean time) with a
    sequence match, so a row one pass skipped shows up as one unmatched row, not as a cascade
    of shifted differences."""
    out: list[Disagreement] = []
    pages = sorted({r["page"] for r in a} | {r["page"] for r in b}, key=int)
    for pg in pages:
        ra = [r for r in a if r["page"] == pg]
        rb = [r for r in b if r["page"] == pg]
        sm = difflib.SequenceMatcher(a=[_key(r) for r in ra], b=[_key(r) for r in rb],
                                     autojunk=False)  # fmt: skip
        for tag, i0, i1, j0, j1 in sm.get_opcodes():
            if tag == "equal" or (tag == "replace" and i1 - i0 == j1 - j0):
                pairs = zip(ra[i0:i1], rb[j0:j1], strict=True)
                for x, y in pairs:
                    for f in DATA_FIELDS:
                        if x[f] != y[f]:
                            out.append(Disagreement(pg, x["line"], y["line"], f, x[f], y[f]))
                continue
            for x in ra[i0:i1]:
                out.append(Disagreement(pg, x["line"], None, "row", _key(x), ""))
            for y in rb[j0:j1]:
                out.append(Disagreement(pg, None, y["line"], "row", "", _key(y)))
    return out


def _sexa(s: str) -> tuple[float, int]:
    """'±d m s.f' → (signed value in the last unit's terms of the top unit, n parts)."""
    s = s.strip()
    sign = -1.0 if s.startswith("-") else 1.0
    parts = [float(p) for p in s.lstrip("+-").split()]
    v = 0.0
    for p in parts:
        v = v * 60.0 + p
    return sign * v, len(parts)


def arithmetic_check(r: dict[str, str]) -> dict[str, float | None]:
    """Printed C − O minus (computed seconds − observed seconds), for RA and Dec, in the
    printed unit (s or ″). The computed value is only printed as its seconds, so its minute is
    unknown: the observed minute −2…+2 that best fits is used (ancient C − O exceed 60″, e.g.
    1690: −66.9″). So this catches misread seconds and C − O digits, not misread minutes. Near
    zero when the row is self-consistent; None where a cell is empty or »."""
    out: dict[str, float | None] = {}
    for name in ("ra", "dec"):
        obs_f, calc_f, cmo_f = f"{name}_obs", f"{name}_calc_sec", f"{name}_c_minus_o"
        vals = (r[obs_f], r[calc_f], r[cmo_f])
        if any(v in ("", NO_VALUE) or "?" in v for v in vals):
            out[name] = None
            continue
        obs, _ = _sexa(vals[0])
        calc_s, cmo = float(vals[1]), float(vals[2])
        sign = -1.0 if obs < 0 else 1.0
        mag = abs(obs)
        base = mag - (mag % 60.0)  # observed value with its seconds removed
        res = (cmo - sign * (base + k * 60.0 + calc_s - mag) for k in range(-2, 3))
        out[name] = round(min(res, key=abs), 6)
    return out


def main() -> None:
    a, b = (read_pass(p) for p in PASS_FILES)
    d = diff_passes(a, b)
    print(f"pass 1: {len(a)} rows, pass 2: {len(b)} rows, {len(d)} disagreements")
    for x in d:
        print(f"p{x.page} l{x.line1}/{x.line2} {x.field}: {x.pass1!r} vs {x.pass2!r}")


if __name__ == "__main__":
    main()
