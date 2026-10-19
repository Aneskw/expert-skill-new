#!/usr/bin/env python3
"""Build the v3.1 incremental submission for source PDF pages 19--21.

These pages were visually inspected and contain bibliography/back matter only.
The batch therefore contributes source context and page coverage, but no new
molecule or reaction records.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills/chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))

from standard_records import write  # noqa: E402


DOC = "doc_c32350973ac9b683"
PAGES = (19, 20, 21)


def evidence(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


init_path = ROOT / "runs/current/raw_model/pages_0019_0021_init.json"
bundle = json.loads(init_path.read_text())

if [document["document_id"] for document in bundle["documents"]] != [DOC]:
    raise RuntimeError("standard init document does not match the requested source")
if bundle["molecules"] or bundle["reactions"]:
    raise RuntimeError("standard init unexpectedly contains records")

source_context = [
    {
        "document_id": DOC,
        "evidence_ids": evidence(19),
        "kind": "references_only",
        "text": (
            "Original PDF page 19 (journal page 5044) was visually inspected. "
            "It contains the continuation of reference 2(f) and references 3--14(j). "
            "These are bibliography entries and cited prior literature, not current-paper "
            "experiments; no figure, table, reaction arrow, footnote, or standalone "
            "characterization block appears on the page."
        ),
    },
    {
        "document_id": DOC,
        "evidence_ids": evidence(20),
        "kind": "references_only",
        "text": (
            "Original PDF page 20 (journal page 5045) was visually inspected. "
            "It continues reference 14(j) and lists references 15--24(e). These are "
            "bibliographic descriptions of prior work, not current-paper experiments; "
            "no figure, table, reaction arrow, footnote, or standalone characterization "
            "block appears on the page."
        ),
    },
    {
        "document_id": DOC,
        "evidence_ids": evidence(21),
        "kind": "references_only",
        "text": (
            "Original PDF page 21 (journal page 5046) was visually inspected. "
            "It continues reference 24(e) and lists references 25--50. Parenthetical "
            "author comments in references 26 and 32 concern literature preparation or "
            "commercial cost context and do not constitute complete current-paper "
            "experiments. No figure, table, reaction arrow, footnote, or standalone "
            "characterization block appears on the page."
        ),
    },
]

bundle["source_context"] = source_context
write(ROOT / "runs/current/submitted.json", bundle)
write(ROOT / "runs/current/raw_model/pages_0019_0021_source_context.json", source_context)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 19,
            "state": "attempted",
            "regions": [
                {
                    "region": "full two-column page: continuation of reference 2(f) and references 3--14(j)",
                    "finding": (
                        "References/back matter only. Cited reviews and synthesis papers are "
                        "prior literature and were not converted into experiments of this paper. "
                        "No figure, table, reaction arrow, footnote, or standalone mass-spectrum "
                        "block is present."
                    ),
                    "evidence_ids": evidence(19),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Cited primary papers were not supplied as source documents, so their experiments were not expanded into records."
            ],
        },
        {
            "page": 20,
            "state": "attempted",
            "regions": [
                {
                    "region": "full two-column page: continuation of reference 14(j) and references 15--24(e)",
                    "finding": (
                        "References/back matter only. Named reactions and cited experimental work "
                        "describe prior literature rather than new experiments reported on this page. "
                        "No figure, table, reaction arrow, footnote, or standalone characterization "
                        "block is present."
                    ),
                    "evidence_ids": evidence(20),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Bibliographic summaries do not provide complete source procedures or row-level conditions and were not treated as reaction records."
            ],
        },
        {
            "page": 21,
            "state": "attempted",
            "regions": [
                {
                    "region": "left column: continuation of reference 24(e) and references 25--38(d)",
                    "finding": (
                        "References/back matter only. Reference 26 includes an author comment about "
                        "commercial cost and cites preparation of allenamide 15a, but it is not a "
                        "complete performed experiment on this page."
                    ),
                    "evidence_ids": evidence(21),
                },
                {
                    "region": "right column: continuation of reference 38(d) and references 39--50",
                    "finding": (
                        "References/back matter only. The reference 32 price comment and remaining "
                        "citations are contextual bibliography, not new current-paper experiments. "
                        "No figure, table, reaction arrow, footnote, or standalone characterization "
                        "block is present."
                    ),
                    "evidence_ids": evidence(21),
                },
            ],
            "record_ids": [],
            "unresolved": [
                "The cited preparation and pricing sources are not attached source documents and were not inferred as Supporting Information or expanded into experiments."
            ],
        },
    ],
}
write(ROOT / "runs/current/coverage/pages_0019_0021.json", coverage)

inventory = {
    "document_id": DOC,
    "pages": list(PAGES),
    "visual_inspection": "completed from original page images at original detail",
    "new_reaction_count": 0,
    "new_molecule_count": 0,
    "source_context_count": len(source_context),
    "record_ids": [],
    "decisions": [
        "All three pages are bibliography/back matter rather than synthetic or mechanistic experiment pages.",
        "Cited prior work was not promoted to experiments of the current paper.",
        "No Supporting Information association was inferred.",
        "No stereochemical molecule review was required because this increment creates no molecule records.",
        "Reagent-dose and outcome-target checks are not applicable to this zero-reaction increment.",
    ],
}
write(ROOT / "runs/current/raw_model/pages_0019_0021_inventory.json", inventory)

print(json.dumps(inventory, ensure_ascii=False))
