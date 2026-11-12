"""Le Verrier (1846) Table n° 77 as a standalone dataset: ``release/leverrier-table77/``.

    uv run python -m fork.export.release      # copy the data files, write datapackage.json

The data files are byte-for-byte copies of the committed project files (``sources()``). The
human spot-check record (``docs/review/spotcheck_human_results.json``) joins the package once
it exists; the build refuses it unless it is a complete human record of the reviewed sample
with no cell marked as not matching. This module copies the files and writes the Frictionless
Data Package descriptor (``datapackage.json``, Data Package v1) with each file's SHA-256, size
and table schema. README.md, MANIFEST.md, DATA_DICTIONARY.md and the licence/citation metadata
are written by hand.
``tests/test_release.py`` checks that the committed release equals this module's output.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from fork.data import spotcheck as sc
from fork.data import transcription as tr

ROOT = tr.ROOT
RELEASE = ROOT / "release" / "leverrier-table77"
REVIEW = sc.REVIEW_DIR

NAME = "leverrier-1846-table77"
TITLE = (
    "Le Verrier (1846), Table n° 77: 279 observations of Uranus, 1690–1845, "
    "transcribed and reconciled"
)
VERSION = "1.0.1"
# The concept DOI (all versions; it resolves to the latest). From v1.0.1 the package cites it with
# the version number: a version's own DOI exists only once Zenodo has the files, and the files
# are frozen by then. Data Package v1 reserves `id` for a global identifier such as a DOI.
DOI = "10.5281/zenodo.24044791"
DOI_V1_0_0 = "10.5281/zenodo.24044792"
LICENSE = {
    "name": "CC0-1.0",
    "path": "https://creativecommons.org/publicdomain/zero/1.0/",
    "title": "CC0 1.0 Universal (public domain dedication)",
}
# How the package states the human spot check (the 0 is enforced by human_checked()).
VERIFIED = "spot-checked by a human against the Gallica scan (28 rows, 336 cells, 0 discrepancies)"
PAGE_COUNTS = {129: 31, 130: 37, 131: 34, 132: 37, 133: 37, 134: 34, 135: 39, 136: 30}

# release file name -> committed source (the human spot check is added by sources())
SOURCES = {
    "leverrier_1846_table77.csv": tr.RECONCILED_FILE,
    "leverrier_1846_misprints.csv": tr.MISPRINT_FILE,
    "leverrier_1846_resolution_log.csv": tr.RESOLUTION_FILE,
    "leverrier_1846_pass1.csv": tr.PASS_FILES[0],
    "leverrier_1846_pass2.csv": tr.PASS_FILES[1],
    "spotcheck_rows.csv": REVIEW / "spotcheck_rows.csv",
    "spotcheck_results.json": REVIEW / "spotcheck_results.json",
    "leverrier_1846.csv": tr.HIST / "leverrier_1846.csv",
}


# --- schemas -------------------------------------------------------------------------------
# One-line descriptions only; DATA_DICTIONARY.md is the full reference.


def _f(name: str, type_: str, desc: str, **constraints: object) -> dict:
    d: dict = {"name": name, "type": type_, "description": desc}
    if constraints:
        d["constraints"] = constraints
    return d


_PAGE = _f("page", "integer", "Printed page of the memoir (129–136).", required=True,
           minimum=129, maximum=136)  # fmt: skip
_LINE = _f("line", "integer", "Data row on that page, from 1.", required=True, minimum=1)
_MONTHS = [
    "Janv.",
    "Fév.",
    "Mars",
    "Avril",
    "Mai",
    "Juin",
    "Juin.",
    "Juill.",
    "Juill",
    "Août",
    "Sept.",
    "Oct.",
    "Nov.",
    "Déc.",
]  # as printed, both variant spellings
_TO_TABLE = [{"fields": ["page", "line"],
              "reference": {"resource": "table77", "fields": ["page", "line"]}}]  # fmt: skip

TABLE77_FIELDS = [
    _PAGE,
    _LINE,
    _f("year", "integer", "Year as printed; blank = same as the row above.",
       minimum=1690, maximum=1845),
    _f("month", "string", "Month as printed (French abbreviation); blank = same as above.",
       enum=_MONTHS),
    _f("day", "integer", "Day of the month (astronomical day, begins at mean noon).",
       required=True, minimum=1, maximum=31),
    _f("mean_time", "string", "Temps moyen, h m s (Paris mean time, astronomical reckoning).",
       required=True, pattern=r"\d{1,2} \d{1,2} \d{1,2}"),
    _f("ra_obs", "string", "Observed apparent RA: arc (° ′ ″) before 1781, time (h m s) after.",
       required=True, pattern=r"\d{1,3} \d{1,2} \d{1,2}\.\d{1,2}"),
    _f("ra_calc_sec", "string", "Seconds of Le Verrier's computed RA (″ before 1781, s after).",
       required=True, pattern=r"\d{1,2}\.\d{1,2}"),
    _f("ra_c_minus_o", "string", "Printed RA computed − observed (″ before 1781, s after).",
       required=True, pattern=r"[+-]\d{1,3}\.\d{1,2}"),
    _f("dec_obs", "string", "Observed apparent Dec, ± ° ′ ″; » = no Dec observed.",
       required=True, pattern=r"[+-]\d{1,2} \d{1,2} \d{1,2}\.\d|»"),
    _f("dec_calc_sec", "string", "Seconds of the computed Dec, ″ (sign as printed); » = none.",
       required=True, pattern=r"[+-]?\d{1,2}\.\d|»"),
    _f("dec_c_minus_o", "string", "Printed Dec computed − observed, ″; » = none.",
       required=True, pattern=r"[+-]\d{1,3}\.\d|»"),
    _f("lon_c_minus_o", "string", "Printed geocentric longitude computed − observed, ″.",
       required=True, pattern=r"[+-]\d{1,3}\.\d"),
    _f("lat_c_minus_o", "string", "Printed geocentric latitude computed − observed, ″; » = none.",
       required=True, pattern=r"[+-]\d{1,3}\.\d|»"),
    _f("notes", "string", "Transcribers' notes, prefixed p1: / p2: by pass."),
]  # fmt: skip

# The raw passes keep '?' for unreadable digits, so only page and line are typed.
PASS_FIELDS = [_PAGE, _LINE] + [
    _f(f["name"], "string", f["description"]) for f in TABLE77_FIELDS[2:]
]

MISPRINT_FIELDS = [
    _PAGE,
    _LINE,
    _f("field", "string", "Table 77 column of the misprinted cell.", required=True,
       enum=list(tr.DATA_FIELDS)),
    _f("printed", "string", "The cell as printed (equals table77).", required=True),
    _f("corrected", "string", "Corrected value; blank for log_only."),
    _f("action", "string", "apply: corrected in the derived file; log_only: recorded only.",
       required=True, enum=["apply", "log_only"]),
    _f("reason", "string", "The argument for the misprint, from the table itself.",
       required=True),
    _f("uses_ephemeris", "string", "Whether the argument relies on a modern ephemeris.",
       required=True, enum=["yes", "no"]),
]  # fmt: skip

RESOLUTION_FIELDS = [
    _PAGE,
    _LINE,
    _f("field", "string", "Table 77 column where the passes disagree.", required=True,
       enum=list(tr.DATA_FIELDS)),
    _f("pass1", "string", "Pass 1 reading (archive.org offprint); ? = unreadable digit.",
       required=True),
    _f("pass2", "string", "Pass 2 reading (Lille offprint); ? = unreadable digit.",
       required=True),
    _f("resolved", "string", "Value adopted in table77.", required=True),
    _f("witness", "string", "Third witness consulted: file and PDF page.", required=True),
    _f("reason", "string", "Why the resolved value was adopted.", required=True),
]  # fmt: skip

SPOTCHECK_FIELDS = [
    _PAGE,
    _LINE,
    _f("reason", "string", "Why the row is in the sample (log, misprint, ancient, page, random).",
       required=True),
    _f("seed", "integer", "Sampling seed.", required=True),
]  # fmt: skip

_STR = "string"
OBSERVATION_FIELDS = [
    _f("id", _STR, "One measured quantity: <obs_id>-ra or <obs_id>-dec.", required=True,
       unique=True),
    _f("obs_id", _STR, "One Table 77 row: LV1846-T77-p<page>-l<line>.", required=True),
    _f("epoch", _STR, "Date and mean time as printed, with the year and month filled in."),
    _f("epoch_tt_jd", "number", "Epoch, Julian Date (TT).", required=True),
    _f("epoch_ut_jd", "number", "Epoch, Julian Date (UT).", required=True),
    _f("time_meridian", _STR, "Meridian of the printed mean time.", enum=["Paris"]),
    _f("time_reckoning", _STR, "Day convention of the printed time."),
    _f("observer", _STR, "Observer, from the memoir's text (ancient rows only)."),
    _f("observatory", _STR, "Observatory, from the memoir's text (1828 inferred).",
       enum=["Greenwich", "Paris", "Göttingen"]),
    _f("observatory_basis", _STR, "'stated', or how the observatory was inferred."),
    _f("instrument", _STR, "Not stated in the memoir; always blank."),
    _f("quantity", _STR, "Measured coordinate.", required=True, enum=["ra", "dec"]),
    _f("value", "number", "Observed value in degrees, after applied misprint corrections.",
       required=True),
    _f("unit", _STR, "Unit of value.", enum=["deg"]),
    _f("value_as_printed", _STR, "The observed cell as printed (before any correction)."),
    _f("unit_as_printed", _STR, "Units of the printed cell."),
    _f("frame_as_printed", _STR, "Reference frame of the printed places."),
    _f("equinox_ref", _STR, "Equinox convention.", enum=["of_date"]),
    _f("equinox_epoch", _STR, "Epoch of a fixed equinox; blank for of_date."),
    _f("equinox_citation", _STR, "Where the memoir states the frame."),
    _f("precession_constant_arcsec_per_yr", _STR, "Precession constant; not_stated here."),
    _f("precession_rate_arcsec_per_yr2", _STR, "Precession rate; blank here."),
    _f("precession_citation", _STR, "Citation for the precession fields."),
    _f("lv_residual", "number", "Le Verrier's O − C (sign flipped from the print), ″."),
    _f("lv_residual_unit", _STR, "Unit and sign convention of lv_residual."),
    _f("lv_c_minus_o_as_printed", _STR, "The printed C − O cell (s or ″)."),
    _f("lv_weight", _STR, "Not printed per observation; always blank."),
    _f("lv_dlon_c_minus_o_as_printed", _STR, "Printed longitude C − O, ″."),
    _f("lv_dlat_c_minus_o_as_printed", _STR, "Printed latitude C − O, ″; » = none."),
    _f("transit_resid_s", "number", "Local apparent sidereal time − observed RA, s."),
    _f("era", _STR, "ancient (1690–1771) or modern (1781–1845).", enum=["ancient", "modern"]),
    _f("source_ref", _STR, "Page, row and text citations."),
    _f("misprints_applied", _STR, "Corrections applied to this row: field: printed -> corrected."),
    _f("notes", _STR, "Attribution notes, log-only misprints, 1828 transit evidence."),
]  # fmt: skip


def _table(name: str, path: str, title: str, desc: str, fields: list[dict], pk: list[str],
           fks: list[dict] | None = None) -> dict:  # fmt: skip
    schema: dict = {"fields": fields, "primaryKey": pk, "missingValues": [""]}
    if fks:
        schema["foreignKeys"] = fks
    return {"name": name, "path": path, "title": title, "description": desc,
            "profile": "tabular-data-resource", "format": "csv", "mediatype": "text/csv",
            "encoding": "utf-8", "schema": schema}  # fmt: skip


RESOURCES = [
    _table("table77", "leverrier_1846_table77.csv", "Table n° 77, reconciled",
           "All 279 rows of Table n° 77 (pp. 129–136) exactly as printed, reconciled from two "
           "independent passes and a third witness. Misprints are kept as printed.",
           TABLE77_FIELDS, ["page", "line"]),
    _table("misprints", "leverrier_1846_misprints.csv", "Misprints in the source",
           "Cells that are wrong in the print itself (all three editions), each with an "
           "argument and an action.", MISPRINT_FIELDS, ["page", "line", "field"], _TO_TABLE),
    _table("resolution_log", "leverrier_1846_resolution_log.csv", "Resolution log",
           "Every cell where pass 1 and pass 2 disagreed, resolved from the third witness.",
           RESOLUTION_FIELDS, ["page", "line", "field"], _TO_TABLE),
    _table("pass1", "leverrier_1846_pass1.csv", "Transcription pass 1 (archive.org offprint)",
           "Independent transcription from the archive.org scan, before reconciliation.",
           PASS_FIELDS, ["page", "line"]),
    _table("pass2", "leverrier_1846_pass2.csv", "Transcription pass 2 (Lille offprint)",
           "Independent transcription from the Lille scan, before reconciliation.",
           PASS_FIELDS, ["page", "line"]),
    _table("spotcheck_rows", "spotcheck_rows.csv", "Spot-check sample",
           "The 28 rows checked cell by cell against the Gallica scan.", SPOTCHECK_FIELDS,
           ["page", "line"], _TO_TABLE),
    {"name": "spotcheck_results", "path": "spotcheck_results.json",
     "title": "Spot-check results (first check)",
     "description": "Per-cell result of the first spot check (336 cells), with the "
                    "reviewer and the SHA-256 of the table checked.",
     "format": "json", "mediatype": "application/json", "encoding": "utf-8"},
    _table("observations", "leverrier_1846.csv", "Derived: one row per measured quantity",
           "Convenience file derived from table77 and the apply misprints: 279 RA and 250 Dec "
           "measurements in degrees, with epochs (JD), observatories and Le Verrier's O − C. "
           "Contains inferences; see DATA_DICTIONARY.md.",
           OBSERVATION_FIELDS, ["id"]),
]  # fmt: skip


HUMAN_RESOURCE = {
    "name": "spotcheck_human_results", "path": sc.HUMAN_RESULTS_FILE.name,
    "title": "Spot-check results (human)",
    "description": "Per-cell result of the human spot check of the same 28 rows (336 cells), "
                   "marked reviewer = human, with the SHA-256 of the table checked.",
    "format": "json", "mediatype": "application/json", "encoding": "utf-8",
}  # fmt: skip


def human_checked() -> bool:
    """Whether the human spot-check record exists. If it does, it must be a complete human
    record of the reviewed sample for the current table with no mismatch, or this raises."""
    path = sc.HUMAN_RESULTS_FILE
    if not path.exists():
        return False
    table = tr.read_pass(tr.RECONCILED_FILE)
    rec = json.loads(path.read_text(encoding="utf-8"))
    bad = sc.check_results(rec, table, sc.read_selection(sc.SELECTION_FILE),
                           sha256(tr.RECONCILED_FILE), human=True)  # fmt: skip
    if bad:
        raise ValueError(f"the human spot check marks {len(bad)} cells as not matching: {bad}; "
                         "resolve them before releasing")  # fmt: skip
    return True


def sources() -> dict[str, Path]:
    out = dict(SOURCES)
    if human_checked():
        out[sc.HUMAN_RESULTS_FILE.name] = sc.HUMAN_RESULTS_FILE
    return out


def resources() -> list[dict]:
    return [*RESOURCES, HUMAN_RESOURCE] if human_checked() else list(RESOURCES)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def descriptor() -> dict:
    """The datapackage.json descriptor for the files currently in ``RELEASE``."""
    review = (
        f"with a first spot check against the scan; then {VERIFIED}"
        if human_checked()
        else "with a first spot check against the scan; a human spot check is pending"
    )
    resources_ = []
    for r in resources():
        p = RELEASE / r["path"]
        resources_.append({**r, "bytes": p.stat().st_size, "hash": f"sha256:{sha256(p)}"})
    return {
        "profile": "data-package",
        "id": f"https://doi.org/{DOI}",
        "name": NAME,
        "title": TITLE,
        "description": (
            "Transcription of Table n° 77 of U. J. Le Verrier, 'Recherches sur les mouvements "
            "de la planète Herschel (dite Uranus)', Connaissance des Temps pour 1849, "
            "Additions, pp. 129–136 (Paris, 1846): the observed apparent right ascensions and "
            "declinations of Uranus that Le Verrier compared with theory before predicting "
            "Neptune, with his computed values and residuals. Transcribed twice, each pass "
            f"without access to the other, and reconciled, {review}. "
            "See README.md."
        ),
        "version": VERSION,
        "homepage": "https://github.com/atishaypahadiya/neptune-fork",
        "licenses": [LICENSE],
        "contributors": [{"title": "Atishay Pahadiya", "role": "author"}],
        "keywords": KEYWORDS,
        "sources": [
            {
                "title": "Le Verrier 1846 offprint (Bachelier), archive.org scan",
                "path": "https://archive.org/details/bub_gb_HMi0NyDVkhQC",
            },
            {
                "title": "Le Verrier 1846 offprint, Université de Lille scan",
                "path": "http://hdl.handle.net/1908/3096",
            },
            {
                "title": "Connaissance des Temps pour 1849, Additions, BnF Gallica scan",
                "path": "https://gallica.bnf.fr/ark:/12148/bpt6k6505616w",
            },
        ],
        "resources": resources_,
    }


KEYWORDS = [
    "Le Verrier",
    "Uranus",
    "Neptune",
    "history of astronomy",
    "astrometry",
    "historical observations",
    "meridian observations",
    "Connaissance des Temps",
    "transcription",
    "celestial mechanics",
]


def build() -> None:
    RELEASE.mkdir(parents=True, exist_ok=True)
    files = sources()
    for name, src in files.items():
        shutil.copyfile(src, RELEASE / name)
    text = json.dumps(descriptor(), indent=2, ensure_ascii=False) + "\n"
    (RELEASE / "datapackage.json").write_text(text, encoding="utf-8")
    print(f"wrote {len(files)} files and datapackage.json to {RELEASE.relative_to(ROOT)}"
          f" (human spot check: {'included' if human_checked() else 'pending'})")  # fmt: skip


if __name__ == "__main__":
    build()
