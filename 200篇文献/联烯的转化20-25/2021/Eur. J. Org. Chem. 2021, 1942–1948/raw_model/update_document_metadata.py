#!/usr/bin/env python3
"""Apply source-visible article metadata before v3 template initialization."""
import json
from pathlib import Path


path = Path("runs/current/dataset.json")
data = json.loads(path.read_text())
document = next(d for d in data["documents"] if d["document_id"] == "doc_a10334c5f49b8d92")
document["title"] = (
    "DMAP Mediated Efficient Construction of Functionalized Chromenes "
    "through One-Pot Reaction of para-Quinone Methides with Allenoates"
)
document["supplementary_material"] = {
    "status": "mentioned_not_supplied",
    "note": (
        "Page 1 explicitly states that Supporting Information is available online, "
        "but no SI file was supplied or explicitly associated in this workspace."
    ),
}
path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
