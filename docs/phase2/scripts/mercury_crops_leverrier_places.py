"""Page bands for the transcription of Le Verrier's Paris meridian places of Mercury, 1801–1842.

    uv run --with pymupdf --with pillow python -I \
        docs/phase2/scripts/mercury_crops_leverrier_places.py

Source: Annales de l'Observatoire impérial de Paris, tome V (1859), chapter XV, pp. 51–61, the
table "Longitudes et latitudes géocentriques de Mercure, déduites des observations méridiennes".
Reads the cached scan (data/cache/mercury/1859AnPar_5_1L.pdf, SHA-256 in
data/historical/mercury/MANIFEST.md); PDF page n is printed page n. The embedded page scans are
4480 px wide, black and white.

The table prints, beside each observed longitude and latitude, the error of the tables (observed
minus computed). Those two columns are painted white on every band before it is written
(docs/phase2/blinding.md). The cut is made without reading them:

- rows are found from the ink profile of the mean-time column only;
- the page is cut into slices at the blank lines between the printed blocks of five rows;
- longitude error: in each slice a white rectangle runs from the first ink-free run of pixel
  columns to the right of the longitudes to the last one to the left of the latitudes (their
  minus signs included), each at least MIN_GAP wide;
- latitude error: from the first such run to the right of the latitudes, every connected piece
  of ink less than BRACE_HEIGHT tall is removed as far as the column of group numbers. A glyph
  of the table is at most about 65 px tall; the braces that tie the rows of a group are taller
  and stay. The errors stand too close to the braces for a straight cut.

Bands: pNNN_b00.png is the heading of the page (no table row); pNNN_b01 … hold 5 to 14 rows
of their own (whole printed blocks) and repeat the last row of the band before and the first row of the band after, so that
a brace cut by a band edge can be followed. Output: data/cache/mercury/renders_leverrier_places/
(ignored).
"""

import hashlib
import itertools
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[3]
PDF = ROOT / "data" / "cache" / "mercury" / "1859AnPar_5_1L.pdf"
PDF_SHA256 = "9668f9a6e5d984ab35efd1a3cec76d2fab40b1e968e268d3d564729e020db3f4"
OUT = ROOT / "data" / "cache" / "mercury" / "renders_leverrier_places"

MIN_GAP = 16  # px, narrowest ink-free run of pixel columns accepted beside a column
BRACE_HEIGHT = 75  # px; ink shorter than this, right of the latitudes, is not a brace

# Per printed page, in pixels of the embedded scan:
#   head  (top, bottom) of the heading band; it ends above the unit marks of the first row
#   body  (top, bottom) of the table rows, unit marks included
#   x     (left, right) of every band
#   rows  (left, right) of the mean-time column, used to find the rows
#   gaps  (lo, hi) of the three searches for an ink-free run: right of the longitudes,
#         left of the latitudes, right of the latitudes; then the x up to which short ink is
#         removed (25 px left of the group numbers, right of every latitude error)
PAGES: dict[int, dict[str, tuple]] = {
    51: {"head": (2700, 3355), "body": (3357, 5165), "x": (280, 3570), "rows": (1151, 1465),
         "gaps": ((1948, 2118), (2292, 2469), (2864, 3051), 3416)},
    52: {"head": (354, 902), "body": (904, 5210), "x": (860, 4140), "rows": (1701, 2022),
         "gaps": ((2514, 2649), (2864, 3018), (3422, 3599), 3947)},
    53: {"head": (275, 823), "body": (825, 5140), "x": (215, 3550), "rows": (1062, 1403),
         "gaps": ((1892, 2032), (2232, 2391), (2812, 2984), 3346)},
    54: {"head": (389, 937), "body": (939, 5235), "x": (900, 4250), "rows": (1752, 2099),
         "gaps": ((2591, 2711), (2936, 3066), (3495, 3644), 4014)},
    55: {"head": (307, 855), "body": (857, 5225), "x": (150, 3465), "rows": (999, 1311),
         "gaps": ((1792, 1967), (2143, 2317), (2712, 2896), 3264)},
    56: {"head": (359, 907), "body": (909, 5270), "x": (950, 4270), "rows": (1798, 2120),
         "gaps": ((2606, 2750), (2970, 3108), (3524, 3679), 4070)},
    57: {"head": (347, 895), "body": (897, 5190), "x": (290, 3600), "rows": (1135, 1443),
         "gaps": ((1998, 2100), (2328, 2450), (2853, 3029), 3367)},
    58: {"head": (393, 941), "body": (943, 5237), "x": (880, 4200), "rows": (1689, 2019),
         "gaps": ((2518, 2646), (2896, 3048), (3454, 3635), 3972)},
    59: {"head": (265, 813), "body": (815, 5187), "x": (250, 3595), "rows": (1064, 1393),
         "gaps": ((1906, 2039), (2289, 2424), (2842, 3014), 3352)},
    60: {"head": (336, 884), "body": (886, 5237), "x": (870, 4185), "rows": (1664, 2016),
         "gaps": ((2511, 2634), (2886, 3013), (3440, 3603), 3935)},
    61: {"head": (270, 818), "body": (820, 3140), "x": (220, 3540), "rows": (1026, 1350),
         "gaps": ((1847, 1998), (2231, 2394), (2790, 2968), 3302)},
}  # fmt: skip


def runs(mask: np.ndarray, join: int) -> list[tuple[int, int]]:
    """Inclusive (start, end) of the runs of True in a 1-d mask; gaps up to `join` are bridged."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(idx) > join)
    starts = np.concatenate(([idx[0]], idx[breaks + 1]))
    ends = np.concatenate((idx[breaks], [idx[-1]]))
    return [(int(s), int(e)) for s, e in zip(starts, ends, strict=True)]


def page_ink(doc: pymupdf.Document, page: int) -> np.ndarray:
    """The embedded scan of a page as a boolean array, True where there is ink."""
    xref = doc[page - 1].get_images()[0][0]
    pix = pymupdf.Pixmap(doc, xref)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    return a[:, :, 0] > 128  # the scan is stored as a mask: high value is ink


def find_rows(ink: np.ndarray, spec: dict[str, tuple]) -> list[tuple[int, int]]:
    """(top, bottom) of each table row, from the mean-time column only."""
    (y0, y1), (x0, x1) = spec["body"], spec["rows"]
    profile = ink[y0:y1, x0:x1].sum(axis=1)
    lines = [(s + y0, e + y0) for s, e in runs(profile > 2, join=3) if e - s >= 20]
    # The unit marks (h m s) stand as a short line just above the first row on most pages.
    if len(lines) > 1 and lines[0][1] - lines[0][0] < 45 and lines[1][0] - lines[0][1] < 20:
        lines = lines[1:]
    return lines


def blocks_of(rows: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """Rows grouped as printed: a blank line separates blocks (of five rows, mostly)."""
    out = [[rows[0]]]
    for prev, row in itertools.pairwise(rows):
        if row[1] - prev[1] > 130:
            out.append([])
        out[-1].append(row)
    return out


def free_run(profile: np.ndarray, lo: int, last: bool, where: str) -> tuple[int, int]:
    """The first (or last) ink-free run at least MIN_GAP wide in a column profile from x=lo."""
    free = [(s + lo, e + lo) for s, e in runs(profile == 0, join=1) if e - s + 1 >= MIN_GAP]
    if not free:
        raise SystemExit(f"no ink-free run of {MIN_GAP} px {where}")
    return free[-1] if last else free[0]


def blank_error_columns(
    ink: np.ndarray, spec: dict[str, tuple], blocks: list[list[tuple[int, int]]]
) -> list[tuple[int, int, int, int, int]]:
    """Remove the two error columns, slice by slice. Returns (top, bottom, three cuts)."""
    y0, y1 = spec["body"]
    (a0, a1), (b0, b1), (c0, c1), x_end = spec["gaps"]
    bounds = [y0] + [(a[-1][1] + b[0][0]) // 2 for a, b in itertools.pairwise(blocks)]
    bounds.append(y1)
    cuts = []
    for top, bottom in itertools.pairwise(bounds):
        where = f"at y {top}..{bottom}"
        left = free_run(ink[top:bottom, a0:a1].sum(axis=0), a0, False, where)[0] + MIN_GAP // 2
        right = free_run(ink[top:bottom, b0:b1].sum(axis=0), b0, True, where)[1] - MIN_GAP // 2
        ink[top:bottom, left:right] = False
        start = free_run(ink[top:bottom, c0:c1].sum(axis=0), c0, False, where)[0] + MIN_GAP // 2
        cuts.append((top, bottom, left, right, start))
    # Latitude error: short pieces of ink between the cut and x_end go; braces stay.
    x_from = min(c[4] for c in cuts)
    labels, _ = ndimage.label(ink[y0:y1, x_from:], structure=np.ones((3, 3)))
    for i, (ys, xs) in enumerate(ndimage.find_objects(labels), start=1):
        mid = y0 + (ys.start + ys.stop) // 2
        start = next(c[4] for c in cuts if c[0] <= mid < c[1])
        x = x_from + xs.start
        if start <= x < x_end and ys.stop - ys.start < BRACE_HEIGHT:
            ink[y0:y1, x_from:][ys, xs][labels[ys, xs] == i] = False
    return cuts


def bands_of(
    blocks: list[list[tuple[int, int]]], body: tuple[int, int]
) -> list[tuple[int, int, int]]:
    """(top, bottom, rows of its own) of each band: whole blocks, 10 rows or more, plus the
    neighbouring row on each side."""
    groups: list[list[tuple[int, int]]] = [[]]
    for block in blocks:
        if len(groups[-1]) >= 10:
            groups.append([])
        groups[-1].extend(block)
    if len(groups) > 1 and len(groups[-1]) < 5 and len(groups[-2]) + len(groups[-1]) <= 15:
        groups[-2].extend(groups.pop())
    out = []
    for i, g in enumerate(groups):
        top = body[0] if i == 0 else groups[i - 1][-1][0] - 15
        bottom = body[1] if i == len(groups) - 1 else groups[i + 1][0][1] + 15
        out.append((top, bottom, len(g)))
    return out


def save(ink: np.ndarray, box: tuple[int, int, int, int], name: str) -> None:
    left, top, right, bottom = box
    Image.fromarray(np.where(ink[top:bottom, left:right], 0, 255).astype(np.uint8)).save(OUT / name)


def main() -> None:
    if hashlib.sha256(PDF.read_bytes()).hexdigest() != PDF_SHA256:
        raise SystemExit(f"{PDF.name}: not the file in the MANIFEST")
    OUT.mkdir(parents=True, exist_ok=True)
    doc = pymupdf.open(PDF)
    total = 0
    for page, spec in PAGES.items():
        ink = page_ink(doc, page).copy()
        rows = find_rows(ink, spec)
        # The heading band is not blanked, so it must end above the first row and its unit marks.
        if not spec["head"][1] <= spec["body"][0] < rows[0][0] - 40:
            raise SystemExit(f"p{page}: the heading band would reach the first row")
        blocks = blocks_of(rows)
        cuts = blank_error_columns(ink, spec, blocks)
        left, right = spec["x"]
        save(ink, (left, spec["head"][0], right, spec["head"][1]), f"p{page:03d}_b00.png")
        bands = bands_of(blocks, spec["body"])
        total += len(rows)
        print(f"p{page:03d}: {len(rows)} rows in blocks of {[len(b) for b in blocks]}")
        for cut in cuts:
            print("  slice y {}..{}: white x {}..{}, short ink removed from x {}".format(*cut))
        print(f"  p{page:03d}_b00.png heading, y {spec['head'][0]}..{spec['head'][1]}")
        for i, (top, bottom, n) in enumerate(bands, start=1):
            save(ink, (left, top, right, bottom), f"p{page:03d}_b{i:02d}.png")
            print(f"  p{page:03d}_b{i:02d}.png {n} rows of its own, y {top}..{bottom}")
    print(f"{total} rows on {len(PAGES)} pages")


if __name__ == "__main__":
    main()
