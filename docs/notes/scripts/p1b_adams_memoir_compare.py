"""Compare the two transcription passes of Adams's memoir token by token.

Run from ``docs/notes/adams_memoir_passes/`` (expects ``pass1/`` and ``pass2/`` beside it).
Keys follow the column specification: normal places by year and column, series and solution by
hypothesis, section and item, equations by hypothesis, section, equation and term. A leading "+"
is dropped on both sides; nothing else is normalised.
"""

import csv
import pathlib
import re

S = pathlib.Path.cwd()


def norm(v: str) -> str:
    v = v.strip().replace("−", "-").replace("–", "-").replace("·", ".").replace(" ", " ")
    v = re.sub(r"\s+", " ", v)
    return v[1:] if v.startswith("+") else v


def load(side: str, name: str, keyf) -> dict | None:
    p = S / side / name
    if not p.exists():
        return None
    out = {}
    with open(p, newline="") as f:
        for r in csv.DictReader(f):
            for k, v in keyf(r):
                out[k] = norm(v)
    return out


def _np(r):
    return [((r["year"], c), r[c]) for c in r if c not in ("year", "kind") and r[c] != ""]


SPECS = {
    "normal_places.csv": _np,
    "series.csv": lambda r: [((r["hypothesis"], r["section"], r["item"]), r["printed"])],
    "equations.csv": lambda r: [
        ((r["hypothesis"], r["section"], r["equation"], r["term"]), r["printed"])
    ],
    "solution.csv": lambda r: [((r["hypothesis"], r["section"], r["item"]), r["printed"])],
}


def main() -> None:
    tot = dis = 0
    for name, keyf in SPECS.items():
        a, b = load("pass1", name, keyf), load("pass2", name, keyf)
        if a is None or b is None:
            print(name, "missing on a side", a is None, b is None)
            continue
        ka, kb = set(a), set(b)
        common = ka & kb
        diffs = [(k, a[k], b[k]) for k in sorted(common) if a[k] != b[k]]
        tot += len(common)
        dis += len(diffs)
        print(f"== {name}: pass1 {len(ka)} keys, pass2 {len(kb)}, common {len(common)}, "
              f"disagree {len(diffs)}")  # fmt: skip
        for k in sorted(ka - kb):
            print("   only pass1:", k, a[k])
        for k in sorted(kb - ka):
            print("   only pass2:", k, b[k])
        for k, x, y in diffs:
            print("   DIFF", k, "pass1", x, "pass2", y)
    print("common tokens", tot, "disagreements", dis)


if __name__ == "__main__":
    main()
