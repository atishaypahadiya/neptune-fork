"""Human spot check of the Phase 1b tables: the sample, the review kit and the results template.

    uv run python docs/notes/scripts/p1b_spotcheck.py select     # docs/review/p1b_spotcheck_rows.csv
    uv run --with pymupdf python docs/notes/scripts/p1b_spotcheck.py kit --out DIR   # the kit
    uv run python docs/notes/scripts/p1b_spotcheck.py template          # the results template
    uv run python docs/notes/scripts/p1b_spotcheck.py check FILE        # read an exported result

The sample. The twelve tables of data/historical/adams_challis/ are strata. Nine were
transcribed from printed pages; from them N_PAGE = 36 rows are drawn, at least MIN_PER_TABLE
from each and the rest in proportion to the number of rows (largest remainders). With 36 rows
and 0 rows in error, the exact (Clopper-Pearson) two-sided 95 % upper bound on the row error
rate is 1 - 0.025**(1/36) = 9.7 % (one-sided: 8.0 %), treating the rows as a simple random
sample; the allocation is not proportional, so the bound is a statement about rows drawn this
way, not a weighted estimate for the tables together. Three tables have no printed page
behind a row (two are computed from the ephemeris, one summarises the literature); EXTRA_EACH
rows of each are drawn as well and shown with their provenance, outside the bound.

The kit shows each sampled row's cells beside the page image or images the row was read from.
The images are whole pages (or the pages of the memoir section), not cell crops: the reviewer
finds the value on the page. Page images come from the cached scans (data/cache/, ignored) and
are written outside the repository.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "historical" / "adams_challis"
CACHE = ROOT / "data" / "cache" / "adams_challis"
REVIEW = ROOT / "docs" / "review"
ROWS_FILE = REVIEW / "p1b_spotcheck_rows.csv"
TEMPLATE_FILE = REVIEW / "p1b_spotcheck_results_template.json"

SEED = 20261112
N_PAGE = 36
MIN_PER_TABLE = 2
EXTRA_EACH = 2
PAMPHLET = "adams_1846_nautical_almanac_appendix_archive_org.pdf"

# Tables read from printed pages, and how a row finds its page.
PAGE_TABLES = {
    "adams_memoir_equations.csv": "memoir section",
    "adams_memoir_normal_places.csv": "memoir sections 8, 9, 51, 52",
    "adams_memoir_series.csv": "memoir section",
    "adams_memoir_solution.csv": "memoir section",
    "adams_memoir_statements.csv": "pdf_page",
    "adams_solutions.csv": "the statements of the same solution",
    "adams_solutions_statements.csv": "image_file",
    "challis_sweeps.csv": "the statements of the same date",
    "challis_sweeps_statements.csv": "image_file",
}
# Tables with no printed page behind a row.
OTHER_TABLES = {
    "neptune_on_challis_dates.csv": "computed from the cached JPL ephemeris "
    "(docs/notes/scripts/p1b_neptune_on_challis_dates.py); Challis's own places in the last "
    "columns are from his printed papers",
    "neptune_track_1846.csv": "computed from the cached JPL ephemeris "
    "(docs/notes/scripts/p1b_neptune_on_challis_dates.py)",
    "prediscovery_recordings.csv": "a summary of the literature named in the row's `sources`; "
    "check against those papers",
}
# image_file prefix -> cached PDF (the number in the name is the PDF page)
PDF_OF = {
    "sampson_p": "ads_1904MmRAS..54..143S.pdf",
    "mnras7_airy_p": "ads_1846MNRAS...7..121A.pdf",
    "mnras7_adams_p": "ads_1846MNRAS...7..149A.pdf",
    "mnras7_challis_zenodo_p": "../../historical/adams_challis/challis_1846_mnras_zenodo.pdf",
    "mnras7_challis_p": "ads_1846MNRAS...7..145C.pdf",
    "memoirs16_p": "memoirs_ras_16_1847_archive_org.pdf",
    "pamphlet_p": PAMPHLET,
}


def read(name: str) -> list[dict[str, str]]:
    with open(DATA / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def allocation() -> dict[str, int]:
    """Rows per page-sourced table: MIN_PER_TABLE each, the rest by largest remainders."""
    n = {t: len(read(t)) for t in PAGE_TABLES}
    rest = N_PAGE - MIN_PER_TABLE * len(n)
    total = sum(n.values())
    quota = {t: rest * k / total for t, k in n.items()}
    out = {t: MIN_PER_TABLE + int(q) for t, q in quota.items()}
    left = N_PAGE - sum(out.values())
    for t in sorted(quota, key=lambda t: (-(quota[t] - int(quota[t])), t))[:left]:
        out[t] += 1
    return {t: min(k, n[t]) for t, k in out.items()}


def select() -> list[dict[str, str]]:
    """The sample: (table, row number from 1, counted without the header), by the fixed seed."""
    rng = np.random.default_rng(SEED)
    picks = []
    alloc = allocation()
    for t in sorted(PAGE_TABLES):
        idx = sorted(rng.choice(len(read(t)), size=alloc[t], replace=False).tolist())
        picks += [{"table": t, "row": str(i + 1), "group": "page"} for i in idx]
    for t in sorted(OTHER_TABLES):
        idx = sorted(rng.choice(len(read(t)), size=EXTRA_EACH, replace=False).tolist())
        picks += [{"table": t, "row": str(i + 1), "group": "no page"} for i in idx]
    for k, p in enumerate(picks, start=1):
        p["id"] = f"R{k:02d}"
    return picks


def write_selection() -> None:
    picks = select()
    with open(ROWS_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "table", "row", "group"], lineterminator="\n")
        w.writeheader()
        w.writerows(picks)
    print(f"wrote {ROWS_FILE}: {len(picks)} rows; allocation {allocation()}")


def upper_bound(n: int, level: float = 0.95, two_sided: bool = True) -> float:
    """Clopper-Pearson upper bound on the error rate after 0 errors in n rows."""
    a = (1 - level) / (2 if two_sided else 1)
    return 1 - a ** (1 / n)


# ------------------------------------------------------------------ pages of a row


def _section_pages() -> dict[str, list[int]]:
    """Memoir section -> PDF pages of the separate printing, from the statements table; a
    section the statements do not name takes the pages between its neighbours."""
    known: dict[int, set[int]] = {}
    for r in read("adams_memoir_statements.csv"):
        known.setdefault(int(r["section"]), set()).add(int(r["pdf_page"]))
    secs = sorted(known)
    out = {}
    for s in range(secs[0], secs[-1] + 1):
        if s in known:
            pages = set(known[s])
            later = [t for t in secs if t > s]
            if later:  # pages between this section's last and the next one's first
                pages |= set(range(max(pages) + 1, min(known[later[0]])))
        else:
            before = max(t for t in secs if t < s)
            after = min(t for t in secs if t > s)
            pages = set(range(max(known[before]), min(known[after]) + 1))
        out[str(s)] = sorted(pages)
    return out


def _images_of(names: str) -> list[tuple[str, int]]:
    out = []
    for name in (x.strip() for x in names.split(";") if x.strip()):
        stem = name.rsplit(".", 1)[0]
        for prefix, pdf in PDF_OF.items():
            if stem.startswith(prefix):
                out.append((pdf, int(stem[len(prefix) :])))
                break
    return out


def pages_of(table: str, row: dict[str, str]) -> list[tuple[str, int]]:
    """(cached PDF, PDF page) for the pages a row was read from."""
    sec = _section_pages()
    if table in ("adams_memoir_equations.csv", "adams_memoir_series.csv",
                 "adams_memoir_solution.csv"):  # fmt: skip
        return [(PAMPHLET, p) for p in sec.get(row["section"], [])]
    if table == "adams_memoir_normal_places.csv":
        return [(PAMPHLET, p) for s in ("8", "9", "51", "52") for p in sec[s]]
    if table == "adams_memoir_statements.csv":
        return [(PAMPHLET, int(row["pdf_page"]))]
    if table in ("adams_solutions_statements.csv", "challis_sweeps_statements.csv"):
        return _images_of(row["image_file"])
    if table == "adams_solutions.csv":
        return sorted({p for r in read("adams_solutions_statements.csv")
                       if r["solution_id"] == row["solution_id"]
                       for p in _images_of(r["image_file"])})  # fmt: skip
    if table == "challis_sweeps.csv":
        st = read("challis_sweeps_statements.csv")
        key = row["date_1846"] or row["row_id"]
        hit = [r for r in st if (r["date_1846"] == key if row["date_1846"]
                                 else r["statement_id"].startswith(row["row_id"]))]  # fmt: skip
        return sorted({p for r in hit for p in _images_of(r["image_file"])})
    return []


# ------------------------------------------------------------------ the kit

STYLE = """
body{font:16px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1100px;padding:16px;
  background:#fff;color:#111}
h1{font-size:1.4rem} h2{font-size:1.1rem;margin-top:2.2rem}
table{border-collapse:collapse;width:100%;margin:.5rem 0}
td,th{border-bottom:1px solid #ccc;padding:4px 8px;text-align:left;vertical-align:top}
td.v{font-family:ui-monospace,monospace;white-space:pre-wrap;word-break:break-word}
img{max-width:100%;border:1px solid #ccc;margin:.4rem 0}
textarea{width:100%} label{margin-right:1rem} .src{color:#444}
@media (prefers-color-scheme:dark){body{background:#111;color:#eee} .src{color:#bbb}
  td,th,img{border-color:#555}}
"""

SCRIPT = """
const KEY = "p1b_spotcheck_" + META.seed;
const state = (() => { try { return JSON.parse(localStorage.getItem(KEY)) || {}; }
  catch (e) { return {}; } })();
const save = () => { try { localStorage.setItem(KEY, JSON.stringify(state)); } catch (e) {} };
document.querySelectorAll("input[type=radio]").forEach((el) => {
  if (state[el.name] === el.value) el.checked = true;
  el.addEventListener("change", () => { state[el.name] = el.value; save(); count(); });
});
document.querySelectorAll("textarea").forEach((el) => {
  if (state[el.name]) el.value = state[el.name];
  el.addEventListener("input", () => { state[el.name] = el.value; save(); });
});
function count() {
  const names = new Set([...document.querySelectorAll("input[type=radio]")].map((e) => e.name));
  const done = [...names].filter((n) => state[n]).length;
  document.getElementById("progress").textContent = done + " of " + names.size + " cells judged";
}
count();
document.getElementById("export").addEventListener("click", () => {
  const rows = META.rows.map((r) => ({
    id: r.id, table: r.table, row: r.row, group: r.group,
    cells: Object.fromEntries(Object.entries(r.cells).map(([k, v]) =>
      [k, { value: v, matches: state[r.id + "|" + k] === "matches" ? true
            : state[r.id + "|" + k] === "differs" ? false : null,
            judgement: state[r.id + "|" + k] || null }])),
    notes: state[r.id + "|notes"] || "",
  }));
  const out = { seed: META.seed, tables_sha256: META.tables_sha256,
    exported: new Date().toISOString(), rows, reviewer: "human",
    reviewer_name: document.getElementById("who").value };
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([JSON.stringify(out, null, 2)],
    { type: "application/json" }));
  a.download = "p1b_spotcheck_human_results.json";
  a.click();
});
"""


def _cells(row: dict[str, str]) -> dict[str, str]:
    return {k: v for k, v in row.items() if v != ""}


def kit(out: Path, width: int = 1700) -> None:
    import pymupdf

    picks = list(csv.DictReader(open(ROWS_FILE, encoding="utf-8")))
    assert picks == select(), "p1b_spotcheck_rows.csv is not what the seed gives"
    (out / "pages").mkdir(parents=True, exist_ok=True)
    docs: dict[str, object] = {}
    made: dict[tuple[str, int], str] = {}

    def image(pdf: str, page: int) -> str | None:
        if (pdf, page) in made:
            return made[pdf, page]
        path = (CACHE / pdf).resolve()
        if not path.exists():
            return None
        doc = docs.setdefault(pdf, pymupdf.open(path))
        if not 1 <= page <= len(doc):
            return None
        name = f"pages/{Path(pdf).stem}_p{page:03d}.jpg"
        pg = doc[page - 1]
        zoom = width / pg.rect.width
        pg.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom)).save(out / name, jpg_quality=80)
        made[pdf, page] = name
        return name

    meta_rows, body = [], []
    for p in picks:
        row = read(p["table"])[int(p["row"]) - 1]
        cells = _cells(row)
        meta_rows.append({**p, "cells": cells})
        if p["group"] == "page":
            pages = pages_of(p["table"], row)
            imgs = [(pdf, pg, image(pdf, pg)) for pdf, pg in pages]
            src = "; ".join(f"{Path(pdf).name} p. {pg}" for pdf, pg, _ in imgs) or "no page found"
            pics = "".join(f'<img loading="lazy" src="{html.escape(n)}" alt="{html.escape(pdf)} '
                           f'page {pg}">' for pdf, pg, n in imgs if n)  # fmt: skip
            missing = [f"{pdf} p. {pg}" for pdf, pg, n in imgs if not n]
            if missing:
                pics += f"<p>Not in the cache: {html.escape('; '.join(missing))}</p>"
        else:
            src, pics = OTHER_TABLES[p["table"]], ""
        trs = "".join(
            f'<tr><th>{html.escape(k)}</th><td class="v">{html.escape(v)}</td><td>'
            + "".join(f'<label><input type="radio" name="{p["id"]}|{html.escape(k)}" '
                      f'value="{j}"> {j}</label>' for j in ("matches", "differs", "unreadable"))
            + "</td></tr>" for k, v in cells.items())  # fmt: skip
        body.append(
            f'<h2 id="{p["id"]}">{p["id"]}: {html.escape(p["table"])}, row {p["row"]}</h2>'
            f'<p class="src">Read from: {html.escape(src)}</p>'
            f"<table><tr><th>cell</th><th>transcribed</th><th>judgement</th></tr>{trs}</table>"
            f'<textarea name="{p["id"]}|notes" rows="2" placeholder="note"></textarea>{pics}'
        )
    meta = {
        "seed": SEED,
        "rows": meta_rows,
        "tables_sha256": {t: sha256(DATA / t) for t in sorted({p["table"] for p in picks})},
    }
    n = sum(p["group"] == "page" for p in picks)
    page = (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>Phase 1b spot check</title><style>{STYLE}</style></head><body>"
        "<h1>Phase 1b tables: human spot check</h1>"
        f"<p>{len(picks)} rows drawn with seed {SEED}: {n} from the tables read from printed "
        f"pages and {len(picks) - n} from the tables with no printed page. For each cell, find "
        "the value on the page image below the row and choose one judgement. A row is in error "
        "if any cell differs. With 0 rows in error among the "
        f"{n}, the exact two-sided 95 % upper bound on the row error rate is "
        f"{100 * upper_bound(n):.1f} %. Judgements are kept in this browser until exported.</p>"
        '<p><label>Reviewer <input id="who"></label> <button id="export">Export results'
        '</button> <span id="progress"></span></p>'
        + "".join(body)
        + f"<script>const META = {json.dumps(meta, ensure_ascii=False)};{SCRIPT}</script>"
        "</body></html>")  # fmt: skip
    (out / "index.html").write_text(page, encoding="utf-8")
    print(f"wrote {out / 'index.html'}: {len(picks)} rows, {len(made)} page images")


def template() -> None:
    """A results file with every judgement empty, in the shape of the Table 77 human record
    (seed, hash of what was checked, export time, rows with cells of value and matches, notes,
    reviewer, reviewer_name); `judgement` adds 'unreadable' beside matches true or false."""
    picks = select()
    rows = []
    for p in picks:
        row = read(p["table"])[int(p["row"]) - 1]
        rows.append({**p, "cells": {k: {"value": v, "matches": None, "judgement": None}
                                    for k, v in _cells(row).items()}, "notes": ""})  # fmt: skip
    out = {"seed": SEED,
           "tables_sha256": {t: sha256(DATA / t) for t in sorted({p["table"] for p in picks})},
           "exported": None, "rows": rows, "reviewer": "human", "reviewer_name": ""}  # fmt: skip
    TEMPLATE_FILE.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {TEMPLATE_FILE}")


def check(path: Path) -> None:
    d = json.loads(path.read_text(encoding="utf-8"))
    assert d["seed"] == SEED and [(r["table"], r["row"]) for r in d["rows"]] == [
        (p["table"], p["row"]) for p in select()], "not this sample"  # fmt: skip
    page = [r for r in d["rows"] if r["group"] == "page"]
    cells = [c for r in page for c in r["cells"].values()]
    bad = [r["id"] for r in page if any(c["matches"] is False for c in r["cells"].values())]
    open_ = sum(c["matches"] is None for c in cells)
    print(f"{len(page)} page-sourced rows, {len(cells)} cells: {len(bad)} rows in error "
          f"{bad}, {open_} cells not judged or unreadable")  # fmt: skip
    if not bad and not open_:
        two, one = upper_bound(len(page)), upper_bound(len(page), two_sided=False)
        print(f"upper bound on the row error rate: {100 * two:.1f} % (two-sided 95 %), "
              f"{100 * one:.1f} % one-sided")  # fmt: skip
    stale = [t for t, h in d["tables_sha256"].items() if sha256(DATA / t) != h]
    print("tables changed since the kit was built:", stale or "none")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("select")
    k = sub.add_parser("kit")
    k.add_argument("--out", type=Path, required=True)
    sub.add_parser("template")
    c = sub.add_parser("check")
    c.add_argument("file", type=Path)
    a = ap.parse_args()
    if a.cmd == "select":
        write_selection()
    elif a.cmd == "kit":
        kit(a.out)
    elif a.cmd == "template":
        template()
    else:
        check(a.file)


if __name__ == "__main__":
    main()
