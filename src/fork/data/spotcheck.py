"""Spot check of the reconciled Table n° 77 against the Gallica CdT 1849 page images.

Picks a fixed-seed, stratified sample of rows and writes ``docs/review/spotcheck.html``. For
each row, the page shows the Gallica row image above the transcribed values, with a checkbox
per cell, a notes field and an export button. Row images go to ``docs/review/crops/``
(gitignored: they are cut from the Gallica scan, see MANIFEST.md); the page, the selection
and this generator are committed.

    uv run --with pymupdf python -I -m fork.data.spotcheck [--overlay DIR]
    uv run --with pymupdf python -I -m fork.data.spotcheck --human     # page for a human review
    uv run python -I -m fork.data.spotcheck [--human] --check RESULTS.json

``--human`` writes ``docs/review/spotcheck_human.html`` for the same reviewed sample. Its ticks
are stored separately from the first review's (so nothing is pre-ticked), export asks for the
reviewer's name, and it downloads ``spotcheck_human_results.json`` marked ``reviewer: human``.
``--check`` validates an exported results file against the table and the sample and lists the
cells marked as not matching (exit status 1 if there are any).

Sample (``select_rows``): every row touched by the resolution log or the misprint table, then
at least ``MIN_ANCIENT`` rows before 1781 and at least one row per printed page, then random
rows up to ``N_ROWS``; all random draws come from ``SPOTCHECK_SEED``.

Row images: rows are found in a wide band over the mean-time and RA columns, from ink runs
on the page's regular row grid (``row_centres``). The table edges and header line per page are
calibrated by hand in ``PAGES``. The row count per page must equal the table's and every gap
must be a whole number of row pitches, or this raises. Each row is then cut as a level strip
along the line through its left and right band positions (the right one from aligning the two
bands' ink profiles), so a slightly tilted scan line comes out level.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import io
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fork.data import transcription as tr

SPOTCHECK_SEED = 20261106
N_ROWS = 28
MIN_ANCIENT = 2
FIRST_MODERN_YEAR = "1781"

GALLICA_PDF = tr.HIST / "sources" / "cdt_1849_additions.pdf"
GALLICA_ARK = "ark:/12148/bpt6k6505616w"
PDF_PAGE_OFFSET = 2  # Gallica PDF page = printed page + 2 (MANIFEST.md)
GALLICA_VIEW_OFFSET = 467  # Gallica view f468 is PDF page 1 (views 468-725 = PDF 1-258)

REVIEW_DIR = tr.ROOT / "docs" / "review"
HTML_FILE = REVIEW_DIR / "spotcheck.html"
HUMAN_HTML_FILE = REVIEW_DIR / "spotcheck_human.html"
RESULTS_FILE = REVIEW_DIR / "spotcheck_results.json"  # the first spot check (2026-11-06)
HUMAN_RESULTS_FILE = REVIEW_DIR / "spotcheck_human_results.json"
SELECTION_FILE = REVIEW_DIR / "spotcheck_rows.csv"
CROP_DIR = REVIEW_DIR / "crops"

# Per printed page, native scan pixels: the table's left and right edges, and the y just below
# the header rule. Calibrated by eye on the embedded page images.
PAGES: dict[str, dict[str, int]] = {
    "129": {"x0": 45, "x1": 848, "y0": 346},
    "130": {"x0": 85, "x1": 880, "y0": 358},
    "131": {"x0": 48, "x1": 850, "y0": 347},
    "132": {"x0": 88, "x1": 882, "y0": 351},
    "133": {"x0": 50, "x1": 842, "y0": 354},
    "134": {"x0": 88, "x1": 882, "y0": 366},
    "135": {"x0": 62, "x1": 862, "y0": 345},
    "136": {"x0": 85, "x1": 882, "y0": 347},
}
CROP_MARGIN = (12, 25)  # px added left and right of the table edges in the row strip
LEFT_BAND = (0.17, 0.44)  # mean time + observed RA, as fractions of the table width
ALIGN_BANDS = ((0.55, 0.80), (0.84, 0.92), (0.93, 1.0))  # Dec; lon C − O; lat C − O
INK_BELOW_PAPER = 65  # grey levels below the local paper level that count as ink
TEXT_HEIGHT = 17.0  # px, a typical row's ink run (old-style figures, with ascenders/descenders)
HALF_HEIGHT = 18  # px above and below the row's centre line in the strip
SCALE = 2  # strip upscale for display

FIELD_LABELS = {
    "year": "year", "month": "month", "day": "day", "mean_time": "mean time",
    "ra_obs": "RA obs", "ra_calc_sec": "RA calc s", "ra_c_minus_o": "RA C−O",
    "dec_obs": "Dec obs", "dec_calc_sec": "Dec calc s", "dec_c_minus_o": "Dec C−O",
    "lon_c_minus_o": "lon C−O", "lat_c_minus_o": "lat C−O",
}  # fmt: skip


@dataclass(frozen=True)
class Pick:
    page: str
    line: str
    reason: str


def _forced(
    resolutions: dict[tuple[str, str, str], dict[str, str]], misprints: list[dict[str, str]]
) -> dict[tuple[str, str], list[str]]:
    out: dict[tuple[str, str], list[str]] = {}
    for page, line, field in resolutions:
        out.setdefault((page, line), []).append(f"resolution log: {field}")
    for m in misprints:
        out.setdefault((m["page"], m["line"]), []).append(f"misprint table: {m['field']}")
    return out


def select_rows(
    table: list[dict[str, str]],
    forced: dict[tuple[str, str], list[str]],
    rng: np.random.Generator,
    n: int = N_ROWS,
) -> list[Pick]:
    """The stratified sample, in table order. Draws happen in a fixed order (ancient top-up,
    then each uncovered page in page order, then the rest), so one seed gives one sample."""
    keys = [(r["page"], r["line"]) for r in table]
    years = [r["year"] for r in table]
    n_ancient = years.index(FIRST_MODERN_YEAR)
    ancient = set(keys[:n_ancient])
    missing = [k for k in forced if k not in set(keys)]
    if missing:
        raise ValueError(f"forced rows not in the table: {missing}")
    chosen: dict[tuple[str, str], str] = {k: "; ".join(v) for k, v in forced.items()}

    def draw(pool: list[tuple[str, str]], reason: str) -> None:
        free = [k for k in pool if k not in chosen]
        chosen[free[int(rng.integers(len(free)))]] = reason

    while sum(k in ancient for k in chosen) < MIN_ANCIENT:
        draw(keys[:n_ancient], "stratum: ancient")
    for pg in sorted({p for p, _ in keys}, key=int):
        if not any(p == pg for p, _ in chosen):
            draw([k for k in keys if k[0] == pg], f"stratum: page {pg}")
    if len(chosen) > n:
        raise ValueError(f"{len(chosen)} rows are forced, more than the sample size {n}")
    while len(chosen) < n:
        draw(keys, "random")
    return [Pick(p, ln, chosen[(p, ln)]) for p, ln in keys if (p, ln) in chosen]


def write_selection(path: Path, picks: list[Pick]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(("page", "line", "reason", "seed"))
        for p in picks:
            w.writerow((p.page, p.line, p.reason, SPOTCHECK_SEED))


def read_selection(path: Path) -> list[Pick]:
    """The committed sample: the record of what was reviewed (2026-11-06). Later additions to
    the resolution log or misprint table do not change it; ``--reselect`` draws a new one."""
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if {int(r["seed"]) for r in rows} != {SPOTCHECK_SEED}:
        raise ValueError(f"{path.name}: seed is not {SPOTCHECK_SEED}")
    return [Pick(r["page"], r["line"], r["reason"]) for r in rows]


def forced_in_selection(picks: list[Pick]) -> dict[tuple[str, str], list[str]]:
    """The forced rows as they were when the sample was drawn, from its recorded reasons."""
    return {(p.page, p.line): p.reason.split("; ") for p in picks if is_flagged(p)}


# --- row images ----------------------------------------------------------------------------


def _page_image(doc: object, printed_page: str) -> np.ndarray:  # doc: a pymupdf Document
    from PIL import Image

    page = doc[int(printed_page) + PDF_PAGE_OFFSET - 1]
    xref = page.get_images()[0][0]
    im = Image.open(io.BytesIO(doc.extract_image(xref)["image"])).convert("L")
    return np.asarray(im)


def _band(page: str, which: tuple[float, float]) -> tuple[int, int]:
    x0, x1 = PAGES[page]["x0"], PAGES[page]["x1"]
    return int(x0 + which[0] * (x1 - x0)), int(x0 + which[1] * (x1 - x0))


def _ink(img: np.ndarray, page: str, which: tuple[float, float]) -> np.ndarray:
    """Ink mask of a column band below the header: darker than the local paper level by
    ``INK_BELOW_PAPER`` (so show-through from the reverse page drops out), with the vertical
    column rules removed (a vertical opening keeps only strokes taller than any glyph)."""
    from scipy.ndimage import binary_dilation, binary_opening, percentile_filter

    xa, xb = _band(page, which)
    band = img[PAGES[page]["y0"] :, xa:xb].astype(float)
    paper = percentile_filter(np.percentile(band, 90, axis=1), 50, size=81)
    ink = band < (paper - INK_BELOW_PAPER)[:, None]
    rules = binary_opening(ink, structure=np.ones((31, 1), bool))
    return ink & ~binary_dilation(rules, structure=np.ones((1, 5), bool))


def _runs(ink: np.ndarray) -> list[tuple[int, int, int]]:
    """(start, end, ink mass) of each run of pixel rows holding ink, at least 7 px tall."""
    has = ink.sum(axis=1) >= 3
    out, start = [], None
    for y, v in enumerate(np.append(has, False)):
        if v and start is None:
            start = y
        elif not v and start is not None:
            if y - start >= 7:
                out.append((start, y, int(ink[start:y].sum())))
            start = None
    return out


def row_centres(img: np.ndarray, page: str, n_rows: int) -> list[float]:
    """y of each data row's centre in the left band.

    Ink runs are classified by height and position. Under 12 px: a units line (h m s) or a
    rule fragment, dropped. A low-ink first run: the units line under the header, dropped.
    Taller runs are rows; two or more rows that touch merge into one run, which is split by
    its height. Rows are centred from the bottom of the run, so a units line merged into the
    first row does not shift it. Runs after the last row (the page footer) must be clearly
    separated from it. Then every gap between consecutive rows must be a whole number of row
    pitches (blank lines between year groups), or this raises."""
    ink = _ink(img, page, LEFT_BAND)
    runs = [r for r in _runs(ink) if r[1] - r[0] >= 12]
    if not runs:
        raise ValueError(f"p.{page}: no text lines found")
    median_mass = float(np.median([m for *_, m in runs]))
    if runs[0][2] < 0.7 * median_mass and runs[0][1] - runs[0][0] < 26:
        runs = runs[1:]
    gaps = np.diff([(a + b) / 2 for a, b, _ in runs])
    pitch = float(np.median(gaps[(gaps > 18) & (gaps < 26)]))
    y0 = PAGES[page]["y0"]
    centres: list[float] = []
    for a, b, _ in runs:
        k = 1 + int((b - a - TEXT_HEIGHT + 6) // pitch)
        centres += [y0 + b - TEXT_HEIGHT / 2 - j * pitch for j in reversed(range(k))]
    if len(centres) < n_rows:
        raise ValueError(f"p.{page}: found {len(centres)} rows, table has {n_rows}")
    rows, extra = centres[:n_rows], centres[n_rows:]
    if extra and extra[0] - rows[-1] < 1.5 * pitch:
        raise ValueError(f"p.{page}: {len(centres)} rows found, table has {n_rows}, and the "
                         "extra line is not separated like a footer")  # fmt: skip
    steps = np.diff(rows) / pitch
    bad = np.flatnonzero((steps < 0.75) | (np.abs(steps - np.rint(steps)) > 0.3))
    if bad.size:
        raise ValueError(f"p.{page}: row spacing off the {pitch:.1f} px grid after line "
                         f"{bad[0] + 1}: {steps[bad[0]]:.2f} pitches")  # fmt: skip
    return rows


Line = tuple[tuple[float, float], ...]  # (x, y) knots of one row's centre line


MAX_STEP = 9  # px; under half the row pitch, so alignment cannot lock onto the next row


def _shift(pa: np.ndarray, pb: np.ndarray, c: int, base: int) -> int:
    """Offset d (within ``base`` ± ``MAX_STEP``) that best aligns profile pb with pa around
    row centre c."""
    lo = max(c - 30, MAX_STEP + abs(base))
    hi = min(c + 31, len(pa) - MAX_STEP - abs(base))
    win = np.arange(lo, hi)
    ds = range(base - MAX_STEP, base + MAX_STEP + 1)
    return max(ds, key=lambda d: float(np.dot(pa[win], pb[win + d])))


def _smooth_fit(y: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Robust quadratic in y through the per-row offsets: tilt and curvature vary smoothly
    over a page, single rows do not."""
    keep = np.ones(len(y), bool)
    for _ in range(3):
        c = np.polyfit(y[keep], d[keep], 2)
        keep = np.abs(np.polyval(c, y) - d) < 3
    return np.polyval(c, y)


def row_lines(img: np.ndarray, page: str, n_rows: int) -> list[Line]:
    """Each row's centre line as (x, y) knots: the left band centre, then each alignment band
    centre with y offset by the scan's tilt and curvature there. Offsets are found by aligning
    each band's ink profile with the previous band's around every row (chained left to right,
    each step within ±``MAX_STEP`` of the previous offset), then replaced by a robust
    quadratic fit in y over the page."""
    left = row_centres(img, page, n_rows)
    y0 = PAGES[page]["y0"]
    ys = np.asarray(left)
    prev = _ink(img, page, LEFT_BAND).sum(axis=1).astype(float)
    prev_d = np.zeros(n_rows, int)
    knots_x = [sum(_band(page, LEFT_BAND)) / 2]
    knots_dy = [np.zeros(n_rows)]
    for band in ALIGN_BANDS:
        prof = _ink(img, page, band).sum(axis=1).astype(float)
        steps = np.array([_shift(prev, prof, round(y) - y0 + int(b), 0)
                          for y, b in zip(left, prev_d, strict=True)])  # fmt: skip
        total = _smooth_fit(ys, prev_d + steps)
        knots_x.append(sum(_band(page, band)) / 2)
        knots_dy.append(total)
        prev, prev_d = prof, np.rint(total).astype(int)
    return [
        tuple((x, y + float(dy[i])) for x, dy in zip(knots_x, knots_dy, strict=True))
        for i, y in enumerate(left)
    ]


def _y_along(line: Line, xs: np.ndarray) -> np.ndarray:
    """y on the row's centre line at each x: piecewise linear through the knots, extended
    linearly beyond the end knots."""
    kx = np.array([k[0] for k in line])
    ky = np.array([k[1] for k in line])
    y = np.interp(xs, kx, ky)
    lo, hi = xs < kx[0], xs > kx[-1]
    y[lo] = ky[0] + (ky[1] - ky[0]) / (kx[1] - kx[0]) * (xs[lo] - kx[0])
    y[hi] = ky[-1] + (ky[-1] - ky[-2]) / (kx[-1] - kx[-2]) * (xs[hi] - kx[-1])
    return y


def row_strip(img: np.ndarray, page: str, line: Line) -> np.ndarray:
    """The row as a level strip: each column x is cut around the row's centre line there, so
    tilt and curvature are removed."""
    x0 = max(PAGES[page]["x0"] - CROP_MARGIN[0], 0)
    x1 = min(PAGES[page]["x1"] + CROP_MARGIN[1], img.shape[1])
    xs = np.arange(x0, x1)
    rows = np.rint(_y_along(line, xs.astype(float))).astype(int)[None, :]
    rows = np.clip(rows + np.arange(-HALF_HEIGHT, HALF_HEIGHT + 1)[:, None], 0, img.shape[0] - 1)
    return img[rows, xs[None, :]]


def _iiif_link(page: str, line: Line, shape: tuple[int, ...]) -> str:
    h, w = shape[:2]
    x0, x1 = PAGES[page]["x0"], PAGES[page]["x1"]
    ys = [y for _, y in line]
    y0, y1 = min(ys) - 3 * HALF_HEIGHT, max(ys) + 3 * HALF_HEIGHT
    region = (100 * x0 / w, 100 * y0 / h, 100 * (x1 - x0) / w, 100 * (y1 - y0) / h)
    view = int(page) + PDF_PAGE_OFFSET + GALLICA_VIEW_OFFSET
    pct = ",".join(f"{v:.2f}" for v in region)
    return f"https://gallica.bnf.fr/iiif/{GALLICA_ARK}/f{view}/pct:{pct}/full/0/native.jpg"


# --- page ----------------------------------------------------------------------------------

STYLE = """
body { font: 15px/1.45 system-ui, sans-serif; color: #1a1a1a; max-width: 1240px;
       margin: 2rem auto; padding: 0 1rem; }
h1 { font-size: 1.4rem; font-weight: 600; margin: 0 0 .5rem; }
section { margin: 2.5rem 0; }
h2 { font-size: 1rem; font-weight: 600; margin: 0 0 .4rem; }
p.meta { margin: 0 0 .5rem; color: #555; }
img { display: block; width: 100%; height: auto; }
table { border-collapse: collapse; width: 100%; margin-top: .5rem;
        font-variant-numeric: tabular-nums; }
th { font-weight: 500; color: #555; text-align: left; padding: .2rem .4rem; }
td { padding: .2rem .4rem; white-space: nowrap; }
.blank { color: #999; }
.meta { color: #555; }
p { margin: .4rem 0; }
textarea { width: 100%; min-height: 2.6rem; font: inherit; margin-top: .4rem; }
button { font: inherit; }
"""

SCRIPT = """
const REVIEWER = document.body.dataset.reviewer || "";
const KEY = "spotcheck:" + (REVIEWER ? REVIEWER + ":" : "") + document.body.dataset.sha;
const saved = JSON.parse(localStorage.getItem(KEY) || "{}");
document.querySelectorAll("input[type=checkbox], textarea").forEach(el => {
  if (el.id in saved) {
    if (el.type === "checkbox") el.checked = saved[el.id];
    else el.value = saved[el.id];
  }
  el.addEventListener("input", save);
});
document.querySelectorAll("section").forEach(sec => {
  if (saved[sec.id + "-all"]) sec.dataset.usedAllMatch = "1";
});
function save() {
  const s = {};
  document.querySelectorAll("input[type=checkbox]").forEach(el => s[el.id] = el.checked);
  document.querySelectorAll("textarea").forEach(el => s[el.id] = el.value);
  document.querySelectorAll("section").forEach(sec => {
    if (sec.dataset.usedAllMatch) s[sec.id + "-all"] = true;
  });
  localStorage.setItem(KEY, JSON.stringify(s));
}
document.querySelectorAll("button.all").forEach(b => b.addEventListener("click", () => {
  const sec = b.closest("section");
  if (sec.dataset.flagged) return;  // logged rows need per-cell ticks
  sec.querySelectorAll("input[type=checkbox]").forEach(c => c.checked = true);
  sec.dataset.usedAllMatch = "1";
  save();
}));
document.getElementById("export").addEventListener("click", () => {
  const nameEl = document.getElementById("reviewer-name");
  if (REVIEWER && !nameEl.value.trim()) { alert("Enter your name before exporting."); return; }
  const rows = [...document.querySelectorAll("section")].map(sec => {
    const cells = {};
    sec.querySelectorAll("input[type=checkbox]").forEach(c => {
      cells[c.dataset.field] = { value: c.dataset.value, matches: c.checked };
    });
    return { page: sec.dataset.page, line: sec.dataset.line, reason: sec.dataset.reason,
             flagged: Boolean(sec.dataset.flagged),
             used_all_match: Boolean(sec.dataset.usedAllMatch),
             cells, notes: sec.querySelector("textarea").value };
  });
  const out = { seed: Number(document.body.dataset.seed), table77_sha256: document.body.dataset.sha,
                exported: new Date().toISOString(),
                all_match_rows: rows.filter(r => r.used_all_match).map(r => r.page + "/" + r.line),
                rows };
  if (REVIEWER) { out.reviewer = REVIEWER; out.reviewer_name = nameEl.value.trim(); }
  const a = document.createElement("a");
  const blob = new Blob([JSON.stringify(out, null, 2)], { type: "application/json" });
  a.href = URL.createObjectURL(blob);
  a.download = document.body.dataset.export || "spotcheck_results.json";
  a.click();
});
"""


def is_flagged(p: Pick) -> bool:
    """In the sample because the resolution log or the misprint table touches it: these rows
    need per-cell ticks, so "All match" is disabled for them."""
    sources = ("resolution log", "misprint table")
    return any(part.startswith(sources) for part in p.reason.split("; "))


def _context(table: list[dict[str, str]]) -> dict[tuple[str, str], str]:
    """Year and month carried down from the last printed ones (for orientation only)."""
    out, year, month = {}, "", ""
    for r in table:
        year, month = r["year"] or year, r["month"] or month
        out[(r["page"], r["line"])] = f"{year} {month} {r['day']}"
    return out


def render_html(
    table: list[dict[str, str]], picks: list[Pick], links: dict, sha: str, human: bool = False
) -> str:
    """The review page. ``human=True``: the page for a human reviewer (separate storage,
    reviewer name required, exports ``spotcheck_human_results.json`` marked as human)."""
    rows = {(r["page"], r["line"]): r for r in table}
    ctx = _context(table)
    e = html.escape
    body = f'data-seed="{SPOTCHECK_SEED}" data-sha="{sha}"'
    if human:
        body += f' data-reviewer="human" data-export="{HUMAN_RESULTS_FILE.name}"'
    parts = [
        "<!DOCTYPE html>",
        '<html lang="en"><head><meta charset="utf-8">',
        f"<title>Table n° 77 {'human ' if human else ''}spot check</title>",
        f"<style>{STYLE}</style></head>",
        f"<body {body}>",
        f"<h1>Le Verrier 1846, Table n° 77: {'human ' if human else ''}spot check</h1>",
        f"<p>{len(picks)} rows, seed {SPOTCHECK_SEED}, from "
        f"<code>data/historical/leverrier_1846_table77.csv</code> (sha256 {sha[:12]}…). "
        "Each image is the row from the Gallica CdT 1849 scan. Tick a cell when the value "
        "below matches the print; an empty cell is correct when nothing is printed there. "
        "Headings give the date with year and month carried down from earlier rows; the year "
        "and month cells are blank where nothing is printed on the row. Ticks and notes are "
        "kept in this browser until you export.</p>",
    ]
    if human:
        parts += [
            "<p><strong>Human review.</strong> Judge each cell from the image (or the "
            "full-resolution Gallica link), not from memory of the earlier review. Don't open "
            f"<code>{RESULTS_FILE.name}</code> or the transcription notes until you have "
            "exported: this check is meant to be independent of the first one. Export saves "
            f"<code>{HUMAN_RESULTS_FILE.name}</code>.</p>",
            '<p><label>Your name <input id="reviewer-name" type="text" '
            'autocomplete="name"></label></p>',
        ]
    parts += ['<p><button id="export">Export results</button></p>']
    for p in picks:
        r = rows[(p.page, p.line)]
        sid = f"p{p.page}l{p.line}"
        cells_h = "".join(f"<th>{FIELD_LABELS[f]}</th>" for f in tr.DATA_FIELDS)
        blank = '<td class="blank">blank</td>'
        vals = "".join(f"<td>{e(r[f])}</td>" if r[f] else blank for f in tr.DATA_FIELDS)
        boxes = "".join(
            f'<td><input type="checkbox" id="{sid}-{f}" data-field="{f}" data-value="{e(r[f])}" '
            f'aria-label="{FIELD_LABELS[f]} matches"></td>'
            for f in tr.DATA_FIELDS
        )
        flagged = is_flagged(p)
        button = (
            '<button class="all" type="button" disabled>All match</button> '
            '<span class="meta">Tick each cell: this row is in the resolution log or the '
            "misprint table.</span>"
            if flagged
            else '<button class="all" type="button">All match</button>'
        )
        attrs = f'id="{sid}" data-page="{p.page}" data-line="{p.line}" data-reason="{e(p.reason)}"'
        if flagged:
            attrs += ' data-flagged="1"'
        parts += [
            f"<section {attrs}>",
            f"<h2>p. {p.page}, line {p.line}: {e(ctx[(p.page, p.line)])}</h2>",
            f'<p class="meta">In the sample by: {e(p.reason)}. Gallica PDF p. '
            f"{int(p.page) + PDF_PAGE_OFFSET}; "
            f'<a href="{e(links[(p.page, p.line)])}">full-resolution crop on Gallica</a></p>',
            f'<img src="crops/{sid}.png" alt="Gallica scan of p. {p.page}, line {p.line}">',
            f"<table><tr>{cells_h}</tr><tr>{vals}</tr><tr>{boxes}</tr></table>",
            f"<p>{button}</p>",
            f'<textarea id="{sid}-notes" placeholder="Notes"></textarea>',
            "</section>",
        ]
    parts += [f"<script>{SCRIPT}</script>", "</body></html>", ""]
    return "\n".join(parts)


def check_results(
    rec: dict, table: list[dict[str, str]], picks: list[Pick], sha: str, human: bool
) -> list[str]:
    """Validate an exported results record against the table (by SHA-256 and by every cell
    value) and the sample. Raises ValueError if it is not a complete record of this sample,
    or (``human=True``) not marked ``reviewer: human`` with a name. Returns the cells marked
    as not matching the print, as ``"p<page> l<line> <field>: <value>"``."""
    if rec.get("table77_sha256") != sha:
        raise ValueError(f"results are for table77 sha256 {rec.get('table77_sha256')!r}, "
                         f"not the current {sha}")  # fmt: skip
    if rec.get("seed") != SPOTCHECK_SEED:
        raise ValueError(f"seed {rec.get('seed')!r} is not {SPOTCHECK_SEED}")
    if human and (rec.get("reviewer") != "human" or not str(rec.get("reviewer_name", "")).strip()):
        raise ValueError("not marked reviewer = human with a reviewer_name")
    got = [(r.get("page"), r.get("line")) for r in rec.get("rows", [])]
    if got != [(p.page, p.line) for p in picks]:
        raise ValueError("the rows are not the reviewed sample, in order")
    cells = {(r["page"], r["line"]): r for r in table}
    bad = []
    for r, p in zip(rec["rows"], picks, strict=True):
        if set(r["cells"]) != set(tr.DATA_FIELDS):
            raise ValueError(f"p{p.page} l{p.line}: cells {sorted(r['cells'])}")
        if r["used_all_match"] and is_flagged(p):
            raise ValueError(f"p{p.page} l{p.line}: 'all match' used on a logged row")
        for f, c in r["cells"].items():
            if c["value"] != cells[(p.page, p.line)][f]:
                raise ValueError(f"p{p.page} l{p.line} {f}: checked {c['value']!r}, table has "
                                 f"{cells[(p.page, p.line)][f]!r}")  # fmt: skip
            if not isinstance(c["matches"], bool):
                raise ValueError(f"p{p.page} l{p.line} {f}: matches is not true/false")
            if not c["matches"]:
                bad.append(f"p{p.page} l{p.line} {f}: {c['value']!r}")
    return bad


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--overlay", type=Path, help="also write per-page row overlays here")
    ap.add_argument(
        "--reselect",
        action="store_true",
        help="draw a new sample from the current logs (replaces the reviewed one)",
    )
    ap.add_argument("--human", action="store_true", help="the page or check for a human review")
    ap.add_argument("--check", type=Path, help="validate an exported results file and exit")
    args = ap.parse_args(argv)

    table = tr.read_pass(tr.RECONCILED_FILE)
    sha = hashlib.sha256(tr.RECONCILED_FILE.read_bytes()).hexdigest()
    if args.check:
        rec = json.loads(args.check.read_text(encoding="utf-8"))
        picks = read_selection(SELECTION_FILE)
        bad = check_results(rec, table, picks, sha, human=args.human)
        n_cells = sum(len(r["cells"]) for r in rec["rows"])
        print(f"{args.check.name}: reviewer {rec.get('reviewer')!r} "
              f"({rec.get('reviewer_name', '')}); {len(picks)} rows, {n_cells} cells, "
              f"{len(bad)} marked as not matching")  # fmt: skip
        for b in bad:
            print("  " + b)
        raise SystemExit(1 if bad else 0)
    if args.human and args.reselect:
        ap.error("--human reviews the committed sample; --reselect would replace it")

    import pymupdf
    from PIL import Image

    if SELECTION_FILE.exists() and not args.reselect:
        picks = read_selection(SELECTION_FILE)  # the reviewed sample stays fixed
    else:
        forced = _forced(tr.read_resolutions(tr.RESOLUTION_FILE),
                         tr.read_misprints(tr.MISPRINT_FILE, table))  # fmt: skip
        picks = select_rows(table, forced, np.random.default_rng(SPOTCHECK_SEED))
        write_selection(SELECTION_FILE, picks)

    doc = pymupdf.open(GALLICA_PDF)
    CROP_DIR.mkdir(parents=True, exist_ok=True)
    links: dict[tuple[str, str], str] = {}
    pages = sorted({r["page"] for r in table}, key=int)
    for pg in pages:  # every page, so a calibration error anywhere is caught
        img = _page_image(doc, pg)
        n = sum(r["page"] == pg for r in table)
        lines = row_lines(img, pg, n)
        if args.overlay:
            _overlay(img, pg, lines, args.overlay)
        for p in picks:
            if p.page != pg:
                continue
            line = lines[int(p.line) - 1]
            strip = Image.fromarray(row_strip(img, pg, line))
            strip = strip.resize((strip.width * SCALE, strip.height * SCALE), Image.LANCZOS)
            strip.save(CROP_DIR / f"p{p.page}l{p.line}.png")
            links[(p.page, p.line)] = _iiif_link(pg, line, img.shape)
    out = HUMAN_HTML_FILE if args.human else HTML_FILE
    out.write_text(render_html(table, picks, links, sha, human=args.human), encoding="utf-8")
    print(f"{len(picks)} rows -> {out.relative_to(tr.ROOT)}; crops in "
          f"{CROP_DIR.relative_to(tr.ROOT)}/")  # fmt: skip
    print(json.dumps({p.page + "/" + p.line: p.reason for p in picks}, indent=0))


def _overlay(img: np.ndarray, page: str, lines: list[Line], out: Path) -> None:
    from PIL import Image, ImageDraw

    out.mkdir(parents=True, exist_ok=True)
    im = Image.fromarray(img).convert("RGB")
    dr = ImageDraw.Draw(im)
    xs = np.arange(PAGES[page]["x0"], PAGES[page]["x1"], 4, dtype=float)
    for i, line in enumerate(lines, 1):
        ys = _y_along(line, xs)
        dr.line(list(zip(xs.tolist(), ys.tolist(), strict=True)), fill=(255, 0, 0))
        dr.text((xs[-1] + 4, ys[-1] - 6), str(i), fill=(255, 0, 0))
    im.save(out / f"overlay_p{page}.png")


if __name__ == "__main__":
    main()
