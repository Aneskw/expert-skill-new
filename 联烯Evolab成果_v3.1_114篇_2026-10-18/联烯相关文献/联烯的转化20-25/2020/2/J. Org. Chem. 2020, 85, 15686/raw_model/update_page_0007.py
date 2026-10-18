#!/usr/bin/env python3
"""Add the bounded page-7 visual evidence and non-reaction context.

Page 7 is article back matter (author information, notes, acknowledgments,
and references).  It contains no reaction scheme, experiment table, reaction
arrow, procedure, or characterization peak table, so this update deliberately
adds no molecule or reaction record.
"""
from __future__ import annotations

import json
from pathlib import Path


RUN = Path(__file__).resolve().parents[1]
DOC = "doc_496d854c8fcb9457"


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path: Path, value) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


bundle = read(RUN / "submitted.json")
template = read(RUN / "raw_model/pages_0007_0007_template.json")

document = bundle["documents"][0]
template_document = template["documents"][0]
known_evidence = {item["evidence_id"] for item in document["evidence"]}
for item in template_document["evidence"]:
    if item.get("page") == 7 and item["evidence_id"] not in known_evidence:
        document["evidence"].append(item)
        known_evidence.add(item["evidence_id"])

context_id = f"{DOC}:page7_back_matter"
if not any(item.get("context_id") == context_id for item in bundle["source_context"]):
    bundle["source_context"].append(
        {
            "context_id": context_id,
            "document_id": DOC,
            "evidence_ids": [
                f"{DOC}:page7:text",
                f"{DOC}:page7:image",
            ],
            "kind": "bibliographic_back_matter",
            "description": (
                "Original PDF page 7 contains the continuation of author "
                "information, contact information and notes, acknowledgments, "
                "and references 1-18. Visual review found no reaction scheme, "
                "experiment table, reaction arrow, synthetic procedure, or "
                "standalone characterization peak table on this page."
            ),
        }
    )

raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "bounded_page_range": [7, 7],
    "visual_pages_opened_original_detail": [7],
    "source_mapping": {
        "presentation_conversion": False,
        "supporting_information_association": "unestablished",
    },
    "page_inventory": [
        {
            "page": 7,
            "regions": [
                {
                    "region": "upper left",
                    "finding": (
                        "Continuation of author affiliations, Yan Qiao entry, "
                        "contact link, Notes, and no-competing-interest statement."
                    ),
                },
                {
                    "region": "middle left",
                    "finding": (
                        "Acknowledgments and funding text; no experimental "
                        "conditions or reaction outcome."
                    ),
                },
                {
                    "region": "lower left and full right column",
                    "finding": (
                        "References 1-18 only. Some titles discuss chemical "
                        "reactions, but these bibliography entries are not "
                        "experiment reports or reaction arrows in this paper."
                    ),
                },
            ],
            "reaction_arrows": 0,
            "experiment_table_rows": 0,
            "new_molecule_records": 0,
            "new_reaction_records": 0,
        }
    ],
    "decisions": {
        "record_creation": (
            "No molecule or reaction record added because page 7 contains "
            "bibliographic back matter only."
        ),
        "references": (
            "Bibliography titles were retained as source context and were not "
            "converted into literature-example experiments."
        ),
        "si": (
            "No SI file is present or explicitly associated; page 7 does not "
            "alter the unestablished SI status."
        ),
    },
    "formal_counts_after_update": {
        "molecules": len(bundle["molecules"]),
        "reactions": len(bundle["reactions"]),
        "source_context": len(bundle["source_context"]),
    },
}

write(RUN / "submitted.json", bundle)
write(RUN / "raw_model/pages_0007_0007_raw_extraction.json", raw)

print(
    json.dumps(
        {
            "page": 7,
            "new_molecules": 0,
            "new_reactions": 0,
            "molecules_total": len(bundle["molecules"]),
            "reactions_total": len(bundle["reactions"]),
            "source_context_total": len(bundle["source_context"]),
        },
        ensure_ascii=False,
    )
)
