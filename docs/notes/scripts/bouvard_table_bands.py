"""Cut bands of rows from the page images of Bouvard's 1821 tables of Uranus.

A band is the heading strip of a table (or of one half of it, where the page prints the
table in two halves) above the rows asked for, with one neighbouring row above and below,
at the resolution of the scan (600 dpi). The bands are what the table rows are entered
from; they are private (the scan is not committed) and go to the cache.

    uv run --with pymupdf --with pillow python -I docs/notes/scripts/bouvard_table_bands.py \
        --rows ROWS.csv [--pdf PDF] [--out DIR] [--check]

ROWS.csv has the columns ``table`` and ``row`` in the names of
``data/historical/adams_challis/bouvard_1821/bouvard_1821_uranus_table_rows.csv`` (other
columns are ignored). The scan is the Harvard copy recorded in
``data/historical/adams_challis/MANIFEST.md``; the script refuses any other file.

The rows are found on the page, not assumed: in each half of a table the lines of print in
the column of the argument are counted, and the count must equal the number of rows the
table has there. ``--check`` does only this and prints the counts.

Output: ``tNN_pPPP_bXX.png`` (NN the table's number, the first of a group printed together;
PPP the page; XX the band on that page) and ``bands_index.csv``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
PDF = ROOT / "data" / "cache" / "adams_challis" / "bouvard" / "bouvard_1821_SfZ4Eiw9yR0C.pdf"
PDF_SHA256 = "ebd3d9fbe1d44a9cc23eeb23225e15f44b226ca29b20119df0bd0c4ea4f6d1f6"
OUT = ROOT / "data" / "cache" / "adams_challis" / "bouvard" / "bands_route_b"
DPI = 600  # the resolution of the scan's own page images
INK = 128  # a pixel darker than this is print
LINE_PITCH = 80  # px: the tables' lines of figures are set this far apart at 600 dpi
MERGE_GAP = 2  # rows of clustered argument values closer than this many rows share a band

MONTHS = ("Janvier", "Février", "Mars", "Avril", "Mai", "Juin", "Juillet", "Août",
          "Septembre", "Octobre", "Novembre", "Décembre")  # fmt: skip


@dataclass(frozen=True)
class Panel:
    """One block of rows on a page. Positions are fractions of the page's width and height."""

    table: str  # the name in the rows file
    number: int  # NN of the band's name
    page: int  # printed page
    leaf: int  # 0-based index in the PDF
    x: tuple[float, float]  # the columns of the block
    head: tuple[tuple[float, float], ...]  # the heading strips, stacked
    body: tuple[float, float]  # where the rows are
    label_x: tuple[float, float]  # the column of the argument, where lines are counted
    labels: tuple[str, ...]  # the rows, top to bottom; empty for a block cut whole


def _range(a: int, b: int, step: int = 1) -> tuple[str, ...]:
    return tuple(str(v) for v in range(a, b + 1, step))


def _months(kind: str) -> tuple[str, ...]:
    return tuple(f"{m} {kind}" for m in MONTHS)


def _x_halves(page: int, leaf: int, first: int, left, right, head, body) -> tuple[Panel, ...]:
    """Table X, one page: two halves of 51 rows; the last row of one is the first of the next."""
    return tuple(
        Panel("X", 10, page, leaf, x, (head,), body, lab, _range(first + 50 * j, first + 50 * j + 50))
        for j, (x, lab) in enumerate((left, right))
    )  # fmt: skip


PANELS = (
    # Table I, two facing pages a row: the left one to argument IV, the right one V to XII
    Panel("I", 1, 84, 125, (0.150, 1.000), ((0.270, 0.349),), (0.349, 0.880), (0.172, 0.252),
          _range(1780, 1809)),
    Panel("I", 1, 85, 126, (0.000, 0.862), ((0.296, 0.373),), (0.373, 0.906), (0.022, 0.100),
          _range(1780, 1809)),
    Panel("I", 1, 86, 127, (0.150, 1.000), ((0.189, 0.267),), (0.267, 0.912), (0.176, 0.255),
          _range(1810, 1854)),
    Panel("I", 1, 87, 128, (0.000, 0.875), ((0.185, 0.264),), (0.264, 0.914), (0.044, 0.125),
          _range(1810, 1854)),
    Panel("I", 1, 88, 129, (0.145, 1.000), ((0.183, 0.261),), (0.261, 0.911), (0.172, 0.253),
          _range(1855, 1900)),
    Panel("I", 1, 89, 130, (0.000, 0.862), ((0.185, 0.264),), (0.264, 0.919), (0.034, 0.118),
          _range(1855, 1900)),
    # Table III: the common years above, the leap years below under their own title
    Panel("III", 3, 91, 132, (0.045, 0.860), ((0.177, 0.305),), (0.305, 0.582), (0.064, 0.131),
          _months("common")),
    Panel("III", 3, 91, 132, (0.045, 0.860), ((0.234, 0.305), (0.584, 0.634)), (0.634, 0.908),
          (0.064, 0.131), _months("leap")),
    Panel("IV", 4, 92, 133, (0.140, 0.955), ((0.255, 0.332),), (0.332, 0.914), (0.166, 0.204),
          _range(1, 31)),
    # the small tables of p. 93 and Table IX are cut whole
    Panel("V", 5, 93, 134, (0.025, 0.440), ((0.040, 0.352),), (0.0, 0.0), (0.0, 0.0), ()),
    Panel("VI", 6, 93, 134, (0.025, 0.440), ((0.352, 0.905),), (0.0, 0.0), (0.0, 0.0), ()),
    Panel("VII", 7, 93, 134, (0.440, 0.850), ((0.040, 0.508),), (0.0, 0.0), (0.0, 0.0), ()),
    Panel("VIII", 8, 93, 134, (0.440, 0.850), ((0.508, 0.905),), (0.0, 0.0), (0.0, 0.0), ()),
    Panel("IX", 9, 94, 135, (0.145, 0.965), ((0.040, 0.800),), (0.0, 0.0), (0.0, 0.0), ()),
    *_x_halves(95, 136, 0, ((0.048, 0.460), (0.065, 0.099)), ((0.448, 0.870), (0.459, 0.492)),
               (0.130, 0.185), (0.185, 0.893)),
    *_x_halves(96, 137, 100, ((0.150, 0.550), (0.167, 0.197)), ((0.540, 0.945), (0.550, 0.582)),
               (0.140, 0.192), (0.192, 0.893)),
    *_x_halves(97, 138, 200, ((0.050, 0.460), (0.066, 0.101)), ((0.448, 0.856), (0.459, 0.488)),
               (0.170, 0.222), (0.222, 0.931)),
    *_x_halves(98, 139, 300, ((0.158, 0.563), (0.176, 0.209)), ((0.552, 0.955), (0.561, 0.595)),
               (0.178, 0.230), (0.230, 0.933)),
    Panel("XI", 11, 99, 140, (0.040, 0.244), ((0.152, 0.207),), (0.207, 0.924), (0.060, 0.097),
          _range(0, 500, 10)),
    Panel("XI", 11, 99, 140, (0.236, 0.440), ((0.152, 0.207),), (0.207, 0.924), (0.246, 0.287),
          _range(500, 1000, 10)),
    Panel("XII", 12, 99, 140, (0.432, 0.636), ((0.152, 0.207),), (0.207, 0.924), (0.442, 0.486),
          _range(0, 5000, 100)),
    Panel("XII", 12, 99, 140, (0.628, 0.842), ((0.152, 0.207),), (0.207, 0.924), (0.638, 0.680),
          _range(5000, 10000, 100)),
    Panel("XIII", 13, 100, 141, (0.155, 0.365), ((0.158, 0.212),), (0.212, 0.925), (0.177, 0.214),
          _range(0, 500, 10)),
    Panel("XIII", 13, 100, 141, (0.353, 0.560), ((0.158, 0.212),), (0.212, 0.925), (0.368, 0.405),
          _range(500, 1000, 10)),
    Panel("XIV", 14, 100, 141, (0.548, 0.757), ((0.158, 0.212),), (0.212, 0.925), (0.563, 0.604),
          _range(0, 500, 10)),
    Panel("XIV", 14, 100, 141, (0.744, 0.958), ((0.158, 0.212),), (0.212, 0.925), (0.759, 0.797),
          _range(500, 1000, 10)),
    Panel("XV-XVII", 15, 101, 142, (0.035, 0.444), ((0.150, 0.203),), (0.203, 0.909),
          (0.053, 0.085), _range(0, 500, 10)),
    Panel("XV-XVII", 15, 101, 142, (0.434, 0.848), ((0.150, 0.203),), (0.203, 0.909),
          (0.445, 0.477), _range(500, 1000, 10)),
    Panel("XVIII-XXI", 18, 102, 143, (0.160, 0.563), ((0.152, 0.206),), (0.206, 0.911),
          (0.177, 0.202), _range(0, 500, 10)),
    Panel("XVIII-XXI", 18, 102, 143, (0.553, 0.958), ((0.152, 0.206),), (0.206, 0.911),
          (0.563, 0.595), _range(500, 1000, 10)),
    Panel("XXXII", 32, 110, 151, (0.138, 0.542), ((0.153, 0.207),), (0.207, 0.913),
          (0.158, 0.216), _range(0, 50)),
    Panel("XXXII", 32, 110, 151, (0.530, 0.945), ((0.153, 0.207),), (0.207, 0.913),
          (0.543, 0.600), _range(50, 100)),
)  # fmt: skip


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


class Scan:
    """The pages of the scan as arrays of grey levels, rendered once each."""

    def __init__(self, pdf: Path) -> None:
        import pymupdf

        if sha256(pdf) != PDF_SHA256:
            raise SystemExit(f"{pdf} is not the scan recorded in the MANIFEST (SHA-256 differs)")
        self.doc = pymupdf.open(pdf)
        self.matrix = pymupdf.Matrix(DPI / 72, DPI / 72)
        self.gray = pymupdf.csGRAY
        self.pages: dict[int, np.ndarray] = {}

    def page(self, leaf: int) -> np.ndarray:
        if leaf not in self.pages:
            pm = self.doc[leaf].get_pixmap(matrix=self.matrix, colorspace=self.gray)
            self.pages[leaf] = np.frombuffer(pm.samples, dtype=np.uint8).reshape(
                pm.height, pm.width
            )
        return self.pages[leaf]


def text_lines(img: np.ndarray, p: Panel) -> list[tuple[int, int]]:
    """The lines of print in the argument's column of a block: (top, bottom) in pixels.

    The tables are set on a regular pattern: lines a fixed distance apart, in groups (of
    five, of three in Table III) with a rule or a space between groups, the last group
    taking any odd line. Vertical and horizontal rules are taken out of the column first
    (columns and pixel rows inked over most of their length); the ink of each pixel row is
    smoothed; and the pattern (first line, spacing, extra space between groups) is the one
    under which the lines carry the most ink. Each line is then moved to the peak of ink
    nearest its place in the pattern."""
    h, w = img.shape
    y0, y1 = int(p.body[0] * h), int(p.body[1] * h)
    x0, x1 = int(p.label_x[0] * w), int(p.label_x[1] * w)
    ink = (img[y0:y1, x0:x1] < INK).copy()
    ink = ink[:, ink.mean(axis=0) < 0.35]
    ink[ink.mean(axis=1) > 0.6, :] = False
    smooth = np.convolve(ink.sum(axis=1).astype(float), np.ones(31) / 31, mode="same")
    n = len(p.labels)
    size = 3 if p.table == "III" else 5
    i = np.arange(n)
    group = np.minimum(i // size, max(n // size - 1, 0))
    best, best_score = None, -1.0
    for pitch in np.arange(LINE_PITCH - 3.0, LINE_PITCH + 3.01, 0.25):
        for extra in np.arange(0.0, 200.0, 2.0):
            rel = i * pitch + group * extra
            room = len(smooth) - 1 - rel[-1]
            if room < 0:
                break
            starts = np.arange(0.0, room + 1.0, 2.0)
            pos = (starts[:, None] + rel[None, :]).astype(int)
            score = smooth[pos].sum(axis=1)
            k = int(np.argmax(score))
            if score[k] > best_score:
                best_score, best = float(score[k]), pos[k]
    if best is None:
        return []
    half = int(0.45 * LINE_PITCH)
    reach = LINE_PITCH // 8
    out = []
    for c in best:
        lo, hi = max(c - reach, 0), min(c + reach + 1, len(smooth))
        c = lo + int(np.argmax(smooth[lo:hi]))
        out.append((c - half + y0, c + half + y0))
    return out


def row_boxes(img: np.ndarray, p: Panel) -> list[tuple[int, int]]:
    """The rows of a block, checked: as many lines as rows, in order, none too near another."""
    lines = text_lines(img, p)
    where = f"Table {p.table}, p. {p.page} ({p.labels[0]} to {p.labels[-1]})"
    if len(lines) != len(p.labels):
        raise SystemExit(f"{where}: {len(lines)} lines found, {len(p.labels)} rows expected")
    gaps = np.diff([a for a, _ in lines])
    if gaps.min() < 0.6 * LINE_PITCH:
        raise SystemExit(f"{where}: two lines {gaps.min()} px apart")
    return lines


def clusters(idx: list[int]) -> list[tuple[int, int]]:
    """Runs of wanted row positions, joined when their neighbours would overlap."""
    out: list[list[int]] = []
    for i in sorted(idx):
        if out and i - out[-1][1] <= MERGE_GAP + 1:
            out[-1][1] = i
        else:
            out.append([i, i])
    return [(a, b) for a, b in out]


def cut(scan: Scan, wanted: dict[str, set[str]], out: Path, check: bool) -> list[dict]:
    from PIL import Image

    index: list[dict] = []
    counter: dict[tuple[int, int], int] = {}
    for p in PANELS:
        img = scan.page(p.leaf)
        h, w = img.shape
        x0, x1 = int(p.x[0] * w), int(p.x[1] * w)
        want = wanted.get(p.table, set())
        if not p.labels:
            if check or not want:
                continue
            n = counter[p.number, p.page] = counter.get((p.number, p.page), 0) + 1
            name = f"t{p.number:02d}_p{p.page:03d}_b{n:02d}.png"
            y0, y1 = p.head[0]
            Image.fromarray(img[int(y0 * h) : int(y1 * h), x0:x1]).save(out / name)
            index.append(dict(band=name, table=p.table, page=p.page, leaf=p.leaf,
                              first_row="whole table", last_row="", rows=" | ".join(sorted(want))))  # fmt: skip
            continue
        boxes = row_boxes(img, p)
        if check:
            pitch = (boxes[-1][0] - boxes[0][0]) / (len(boxes) - 1)
            print(f"Table {p.table:9s} p. {p.page:3d}  {p.labels[0]:>15s} to {p.labels[-1]:<15s} "
                  f"{len(boxes):3d} rows, {pitch:5.1f} px apart")  # fmt: skip
            continue
        idx = [i for i, lab in enumerate(p.labels) if lab in want]
        gap = np.full((24, x1 - x0), 255, dtype=np.uint8)
        heads: list[np.ndarray] = []
        for y0, y1 in p.head:
            heads += [img[int(y0 * h) : int(y1 * h), x0:x1], gap]
        for a, b in clusters(idx):
            lo, hi = max(a - 1, 0), min(b + 1, len(boxes) - 1)
            pad = int(0.6 * (boxes[1][0] - boxes[0][0]))
            top = max(boxes[lo][0] - pad, int(p.body[0] * h))
            bottom = min(boxes[hi][1] + pad, h)
            body = img[top:bottom, x0:x1]
            n = counter[p.number, p.page] = counter.get((p.number, p.page), 0) + 1
            name = f"t{p.number:02d}_p{p.page:03d}_b{n:02d}.png"
            Image.fromarray(np.vstack([*heads, body])).save(out / name)
            index.append(dict(band=name, table=p.table, page=p.page, leaf=p.leaf,
                              first_row=p.labels[lo], last_row=p.labels[hi],
                              rows=" | ".join(p.labels[i] for i in idx if a <= i <= b)))  # fmt: skip
    return index


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rows", type=Path, help="CSV with the columns table and row")
    ap.add_argument("--pdf", type=Path, default=PDF)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--check", action="store_true", help="count the rows on each page and stop")
    a = ap.parse_args()
    scan = Scan(a.pdf)
    if a.check:
        cut(scan, {}, a.out, True)
        return
    if a.rows is None:
        ap.error("--rows is needed unless --check is given")
    wanted: dict[str, set[str]] = {}
    with a.rows.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            wanted.setdefault(r["table"], set()).add(r["row"])
    known = {p.table: set() for p in PANELS}
    for p in PANELS:
        known[p.table].update(p.labels)
    for table, rows in wanted.items():
        if table not in known:
            raise SystemExit(f"Table {table} has no block defined here")
        missing = rows - known[table] if known[table] else set()
        if missing:
            raise SystemExit(f"Table {table}: rows not on any block: {sorted(missing)}")
    a.out.mkdir(parents=True, exist_ok=True)
    index = cut(scan, wanted, a.out, False)
    with (a.out / "bands_index.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(index[0]))
        w.writeheader()
        w.writerows(index)
    print(f"{len(index)} bands in {a.out}")


if __name__ == "__main__":
    main()
