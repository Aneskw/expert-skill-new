#!/usr/bin/env python3
"""Add the visually reviewed, non-reaction inventory for original PDF page 7."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/current"
DOC = "doc_93797a8609d47752"
TEXT_EVIDENCE = f"{DOC}:page7:text"
IMAGE_EVIDENCE = f"{DOC}:page7:image"
CONTEXT_ID = f"{DOC}:page7_reference_inventory"


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def main() -> None:
    submitted_path = OUT / "submitted.json"
    bundle = json.loads(submitted_path.read_text())
    init_template = json.loads(
        (OUT / "raw_model/pages_0007_0007_init_template.json").read_text()
    )
    # Refresh the immutable source/evidence envelope from this pass's standard-init
    # template.  The older bundle predates the prepared page-7 image evidence.
    bundle["documents"] = init_template["documents"]

    # Page 7 contains only the end of reference 37 and references 38--40.
    # Replace by stable context_id on rerun without touching reaction records.
    page_context = {
        "context_id": CONTEXT_ID,
        "document_id": DOC,
        "evidence_ids": [TEXT_EVIDENCE, IMAGE_EVIDENCE],
        "kind": "nonreaction_page",
        "location": {
            "page": 7,
            "region": "reference-list continuation (reference 37 ending and references 38–40)",
        },
        "note": (
            "Original PDF page 7 was visually inspected. It contains bibliographic "
            "references only and no current-paper experiment, reaction arrow, table, "
            "figure, footnote, structure drawing, or independent characterization result; "
            "the cited publications were not promoted to reaction records."
        ),
    }
    contexts = [
        item for item in bundle.get("source_context", [])
        if item.get("context_id") != CONTEXT_ID
    ]
    contexts.append(page_context)
    bundle["source_context"] = contexts
    write_json(submitted_path, bundle)

    raw_inventory = {
        "document_id": DOC,
        "page_range": [7, 7],
        "semantic_accuracy_verified": False,
        "note": (
            "Raw page inventory retained separately from normalized records. "
            "The page image was reviewed at original detail; no reaction candidate was created."
        ),
        "pages": [
            {
                "page": 7,
                "regions": [
                    "upper-left reference-list continuation",
                    "remaining blank page area and publication footer",
                ],
                "findings": [
                    "The upper-left text completes reference 37 and lists references 38–40.",
                    "No current-paper experiment, reaction arrow, structure, table, figure, footnote, or standalone characterization result is present.",
                    "The cited publications are bibliography only and were not promoted to literature-example or current-paper reaction records.",
                ],
                "record_ids": [],
                "evidence_ids": [TEXT_EVIDENCE, IMAGE_EVIDENCE],
            }
        ],
    }
    write_json(OUT / "raw_model/pages_0007_0007_extraction.json", raw_inventory)

    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 7,
                "state": "attempted",
                "regions": [
                    {
                        "region": "upper-left reference-list continuation",
                        "finding": (
                            "Reference 37 ends and references 38–40 are listed; these are cited works, "
                            "not current-paper experiments, and no reaction records were created."
                        ),
                        "evidence_ids": [TEXT_EVIDENCE, IMAGE_EVIDENCE],
                    },
                    {
                        "region": "remaining blank page area and publication footer",
                        "finding": (
                            "No reaction arrow, scheme, table, figure, footnote, experimental prose, "
                            "or independent characterization observation is present."
                        ),
                        "evidence_ids": [IMAGE_EVIDENCE],
                    },
                ],
                "record_ids": [],
                "unresolved": [],
            }
        ],
    }
    write_json(OUT / "coverage/pages_0007_0007.json", coverage)

    print(json.dumps({
        "page_range": [7, 7],
        "new_molecules": 0,
        "new_reactions": 0,
        "page_record_ids": [],
        "source_context_id": CONTEXT_ID,
        "total_molecules": len(bundle.get("molecules", [])),
        "total_reactions": len(bundle.get("reactions", [])),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
