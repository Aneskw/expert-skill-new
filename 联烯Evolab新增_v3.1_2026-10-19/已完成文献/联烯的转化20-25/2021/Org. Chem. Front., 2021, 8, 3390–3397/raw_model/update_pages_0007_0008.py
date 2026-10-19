#!/usr/bin/env python3
"""Idempotently add the page-7--8 visual finding to the current v3.1 bundle.

Pages 7--8 are bibliography only, so this update creates no molecule or
reaction records.  It preserves the previously imported records verbatim and
adds one evidence-linked source-context entry before standard merge.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SUBMITTED = ROOT / "runs/current/submitted.json"
INIT_TEMPLATE = ROOT / "runs/current/raw_model/pages_0007_0008_init.json"
CONTEXT_ID = "doc_d31888e9faa3be67:bibliography_pages7_8"


bundle = json.loads(SUBMITTED.read_text())
init_template = json.loads(INIT_TEMPLATE.read_text())
before = {
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
}

# The prior page-1--6 submission snapshot omitted the rendered page-7 and
# page-8 image evidence objects.  Standard merge requires the prepared source
# block to match exactly, so use the freshly initialized authoritative block.
# This changes no molecule or reaction content.
bundle["documents"] = init_template["documents"]

entry = {
    "context_id": CONTEXT_ID,
    "kind": "literature_context",
    "description": (
        "Original PDF pages 7-8 are bibliography continuation only: page 7 "
        "contains the end of reference 11 and references 12-19(a), and page 8 "
        "completes reference 19. Cited prior-work experiments are not imported "
        "as current-paper experiments. No standalone characterization dataset "
        "appears on these pages."
    ),
    "evidence_ids": [
        "doc_d31888e9faa3be67:page7:text",
        "doc_d31888e9faa3be67:page7:image",
        "doc_d31888e9faa3be67:page8:text",
        "doc_d31888e9faa3be67:page8:image",
    ],
}

matches = [
    item for item in bundle["source_context"]
    if item.get("context_id") == CONTEXT_ID
]
if matches:
    if len(matches) != 1:
        raise RuntimeError(f"duplicate source_context id: {CONTEXT_ID}")
    matches[0].clear()
    matches[0].update(entry)
else:
    bundle["source_context"].append(entry)

after = {
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
}
if after["molecules"] != before["molecules"] or after["reactions"] != before["reactions"]:
    raise RuntimeError("bibliography-only update unexpectedly changed records")

SUBMITTED.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"before": before, "after": after}, ensure_ascii=False))
