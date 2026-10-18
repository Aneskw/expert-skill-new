#!/usr/bin/env python3
"""Record the visual review of PDF pages 7-8 for the current document.

Both pages are bibliography-only.  This script therefore adds source-context
entries and coverage/task metadata without inventing reaction or molecule
records.  It is idempotent and preserves the existing page 1-6 extraction.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_a3d4333b7270ebfb"


def dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


bundle_path = RUN / "submitted.json"
bundle = json.loads(bundle_path.read_text())
template = json.loads((RUN / "raw_model/pages_0007_0008_template.json").read_text())
# `standard merge` requires prepared document/evidence metadata to match the
# mechanical dataset exactly.  Reuse the freshly initialized template rather
# than reconstructing or partially copying that protected metadata.
bundle["documents"] = template["documents"]
assert any(document.get("document_id") == DOC for document in bundle["documents"])

contexts = {
    item.get("context_id"): item for item in bundle.get("source_context", [])
}
contexts[f"{DOC}:page7_references"] = {
    "context_id": f"{DOC}:page7_references",
    "document_id": DOC,
    "kind": "bibliography",
    "description": (
        "PDF page 7 contains only bibliographic references: the continuation "
        "of reference 10(g), references 11-24, and the beginning of reference "
        "25. These cited publications are not experiments performed in the "
        "current paper and were not converted into reaction records."
    ),
    "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
}
contexts[f"{DOC}:page8_references"] = {
    "context_id": f"{DOC}:page8_references",
    "document_id": DOC,
    "kind": "bibliography",
    "description": (
        "PDF page 8 contains only bibliographic references: the continuation "
        "of reference 25 and references 26-27. These cited publications are "
        "not experiments performed in the current paper and were not converted "
        "into reaction records."
    ),
    "evidence_ids": [f"{DOC}:page8:image", f"{DOC}:page8:text"],
}
bundle["source_context"] = list(contexts.values())
dump(bundle_path, bundle)

raw = {
    "document_id": DOC,
    "pages_inspected": [7, 8],
    "visual_inspection_completed": True,
    "chemical_accuracy_verified": False,
    "page_findings": [
        {
            "page": 7,
            "finding": (
                "Bibliography only: continuation of reference 10(g), references "
                "11-24, and the beginning of reference 25."
            ),
            "reaction_arrows": 0,
            "tables": 0,
            "current_paper_experiments": 0,
            "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
        },
        {
            "page": 8,
            "finding": (
                "Bibliography only: continuation of reference 25 and references "
                "26-27."
            ),
            "reaction_arrows": 0,
            "tables": 0,
            "current_paper_experiments": 0,
            "evidence_ids": [f"{DOC}:page8:image", f"{DOC}:page8:text"],
        },
    ],
    "new_molecules": 0,
    "new_reactions": 0,
    "new_source_context_entries": 2,
    "limitations": [
        "The bibliography was classified as source context rather than current-paper experimental evidence.",
        "No supporting-information file was supplied or explicitly associated.",
        "No chemistry record was source-reviewed; existing machine records remain needs_review.",
    ],
}
dump(RUN / "raw_model/pages_0007_0008_extraction.json", raw)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "two-column bibliography, full page",
                    "finding": (
                        "References only: continuation of reference 10(g), references "
                        "11-24, and the beginning of reference 25. No current-paper "
                        "reaction arrow, experimental table, experimental footnote, "
                        "or characterization peak table is present."
                    ),
                    "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                }
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "single-column bibliography, upper-left page area",
                    "finding": (
                        "References only: continuation of reference 25 and references "
                        "26-27. No current-paper reaction arrow, experimental table, "
                        "experimental footnote, or characterization peak table is present."
                    ),
                    "evidence_ids": [f"{DOC}:page8:image", f"{DOC}:page8:text"],
                }
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
dump(RUN / "coverage/pages_0007_0008.json", coverage)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") in {
        f"{DOC}:page7:text",
        f"{DOC}:page8:text",
    }:
        task["state"] = "attempted"
        task["reason"] = (
            "Original page image inspected and mapped in "
            "coverage/pages_0007_0008.json; bibliography only, so no reaction "
            "record was created."
        )
dump(tasks_path, tasks)

print(
    json.dumps(
        {
            "pages": [7, 8],
            "new_molecules": 0,
            "new_reactions": 0,
            "source_context_entries": len(bundle["source_context"]),
        },
        ensure_ascii=False,
    )
)
