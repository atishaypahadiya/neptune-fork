"""site/data/release.json: the code's release version, date and archive DOI, read from
pyproject.toml and CITATION.cff so that no number is typed into the page.

    uv run python docs/notes/scripts/release_json.py

The software DOI (``doi:`` in CITATION.cff) and its concept DOI (an ``identifiers`` entry of
type doi whose description names it) are absent until Zenodo mints them for the release; the
site then says the archive is pending.
"""

import json
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "site" / "data" / "release.json"
REPO = "https://github.com/atishaypahadiya/neptune-fork"


def main() -> None:
    version = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    date = re.search(r'^date-released: "?([0-9-]+)"?', cff, re.M)
    doi = re.search(r'^doi: "?([^"\n]+)"?', cff, re.M)
    concept = re.search(r'value: "?(10\.5281/zenodo\.\d+)"?\n\s+description: "?[^"\n]*[Cc]oncept', cff)
    out = {
        "description": "The code's release, from pyproject.toml and CITATION.cff: version, tag, "
        "release date, and the Zenodo DOIs of the archived code once minted (null before).",
        "version": version,
        "tag": f"v{version}",
        "date_released": date.group(1) if date else None,
        "software_doi": doi.group(1) if doi else None,
        "software_concept_doi": concept.group(1) if concept else None,
        "repository": REPO,
        "release_url": f"{REPO}/releases/tag/v{version}",
        "license": "MIT (code); CC0 1.0 (data/historical/)",
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {OUT}: {out['tag']}, {out['date_released']}, DOI {out['software_doi']}")


if __name__ == "__main__":
    main()
