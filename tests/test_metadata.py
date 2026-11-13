"""Deposit metadata: .zenodo.json parses and carries what Zenodo needs to mint a version."""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = ("title", "creators", "upload_type", "license", "description")
UPLOAD_TYPES = {"software", "dataset", "publication", "poster", "presentation", "image",
                "video", "lesson", "physicalobject", "other"}  # fmt: skip


def test_zenodo_json_has_the_required_keys():
    meta = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    missing = [k for k in REQUIRED if not meta.get(k)]
    assert not missing, f".zenodo.json lacks {missing}"
    assert meta["upload_type"] in UPLOAD_TYPES
    assert isinstance(meta["title"], str) and isinstance(meta["description"], str)
    assert isinstance(meta["license"], str | dict)
    creators = meta["creators"]
    assert isinstance(creators, list) and creators
    assert all(isinstance(c, dict) and c.get("name") for c in creators)
    for rel in meta.get("related_identifiers", []):
        assert rel.get("identifier") and rel.get("relation")
