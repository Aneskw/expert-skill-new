#!/usr/bin/env python3
"""Build the bounded v3.1 page-7 extraction delta.

Page 7 was visually inspected at original resolution.  It contains the end of
the acknowledgments and the bibliography, but no experiment, reaction arrow,
table, figure, footnote, or standalone characterization dataset.  The delta
therefore adds source context only and deliberately creates no reaction or
molecule records.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_975f54c401a61631"
TEXT_EVIDENCE = f"{DOC}:page7:text"
IMAGE_EVIDENCE = f"{DOC}:page7:image"
EVIDENCE = [TEXT_EVIDENCE, IMAGE_EVIDENCE]


template_path = RUN / "raw_model/submitted_pages_0007_0007.template.json"
bundle = json.loads(template_path.read_text())
assert bundle["molecules"] == []
assert bundle["reactions"] == []
assert bundle["source_context"] == []

bundle["source_context"] = [
    {
        "context_id": f"{DOC}:context:page7_acknowledgments_references",
        "document_id": DOC,
        "kind": "acknowledgments_and_bibliography",
        "evidence_ids": EVIDENCE,
        "note": (
            "Original PDF page 7 contains the closing acknowledgments and references 1-27. "
            "It contains no figure, table, reaction arrow, footnote, experimental procedure, "
            "or standalone characterization peak list. References 12-17 provide bibliography "
            "for prior-work reactions discussed and drawn earlier in the article; the citation "
            "list itself is not an independent synthesis experiment."
        ),
    }
]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left continuation above REFERENCES",
                    "finding": (
                        "Closing acknowledgments and funding statement only; no chemical "
                        "experiment, reaction condition, result, or molecular structure."
                    ),
                    "evidence_ids": EVIDENCE,
                },
                {
                    "region": "left column, REFERENCES 1-18",
                    "finding": (
                        "Bibliographic citations only. References 12-17 concern prior allene/"
                        "sulfone chemistry cited in the article, but this page supplies no "
                        "reaction arrow or experiment-specific conditions and is not counted "
                        "as a new literature-example experiment."
                    ),
                    "evidence_ids": EVIDENCE,
                },
                {
                    "region": "right column, REFERENCES 19-27",
                    "finding": (
                        "Bibliographic citations only; no figure, table, reaction arrow, "
                        "footnote, synthesis procedure, control experiment, or characterization "
                        "peak table is present."
                    ),
                    "evidence_ids": EVIDENCE,
                },
            ],
            "record_ids": [],
            "unresolved": [
                "Bibliography entries alone do not report enough source detail to create new "
                "reaction records; no reaction identity, quantities, outcomes, or structures "
                "were inferred from citation titles."
            ],
        }
    ],
}

assert len(coverage["pages"]) == 1
assert coverage["pages"][0]["page"] == 7
assert coverage["pages"][0]["regions"]

serialized_bundle = json.dumps(bundle, ensure_ascii=False, indent=2) + "\n"
(RUN / "submitted.json").write_text(serialized_bundle)
(RUN / "raw_model/pages_0007_0007_extraction.json").write_text(serialized_bundle)
(RUN / "coverage/pages_0007_0007.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n"
)

print(
    json.dumps(
        {
            "page_range": "7-7",
            "molecules": 0,
            "reactions": 0,
            "source_context": 1,
            "coverage_pages": 1,
        },
        ensure_ascii=False,
    )
)
