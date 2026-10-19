#!/usr/bin/env python3
"""Add visually reviewed nonexperimental context for PDF pages 13--15."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
DOC = "doc_dfa009e2c8318d59"


def evidence(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


bundle = json.loads((RUN / "submitted.json").read_text())
# The merge contract requires the prepared dataset's document/evidence registry
# byte-for-byte, while the molecule/reaction content remains the accumulated
# v3.1 submission from the preceding page batches.
bundle["documents"] = json.loads((RUN / "dataset.json").read_text())["documents"]

contexts = [
    {
        "context_id": "ctx_page13_references",
        "type": "references_and_literature_context",
        "evidence_ids": evidence(13),
        "summary": (
            "Reference-list continuation only (reference 1(f) through reference 6(i)). "
            "The page contains bibliographic citations to prior ferrocene, pyrrole, and "
            "materials literature, with no reaction scheme, table, experimental procedure, "
            "product outcome, or characterization block from the present paper."
        ),
    },
    {
        "context_id": "ctx_page14_references",
        "type": "references_and_literature_context",
        "evidence_ids": evidence(14),
        "summary": (
            "Reference-list continuation only (reference 6(j) through reference 13). "
            "Titles describing reactions and syntheses are citations to prior work, not "
            "experimental records of the present paper; no reaction arrows, tables, or "
            "present-paper outcome data occur on this page."
        ),
    },
    {
        "context_id": "ctx_page15_references",
        "type": "references_and_literature_context",
        "evidence_ids": evidence(15),
        "summary": (
            "Final reference-list page only (references 14 through 24), including cited "
            "reaction literature and crystallographic software references. No new experiment, "
            "reaction arrow, table row, footnote condition, yield, or spectrum from the present "
            "paper is reported."
        ),
    },
]

existing = {item.get("context_id") for item in bundle.get("source_context", [])}
for item in contexts:
    if item["context_id"] not in existing:
        bundle["source_context"].append(item)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page, two-column reference list",
                    "finding": (
                        "Reference-list continuation from item 1(f) through item 6(i). "
                        "Visual inspection found no present-paper reaction arrow, table, "
                        "experimental procedure, outcome, footnote condition, or standalone "
                        "characterization dataset."
                    ),
                    "evidence_ids": evidence(13),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "Cited prior-work experiments were not converted into records because this page supplies bibliographic references rather than source procedures or outcomes for the present paper."
            ],
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page, two-column reference list",
                    "finding": (
                        "Reference-list continuation from item 6(j) through item 13. "
                        "Reaction and synthesis terms occur only inside cited titles; visual "
                        "inspection found no present-paper scheme, table, experimental arrow, "
                        "yield row, or characterization block."
                    ),
                    "evidence_ids": evidence(14),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "The chemical details of cited works are outside this document-page extraction and were not inferred from citation titles."
            ],
        },
        {
            "page": 15,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper half, two-column final reference list; lower half blank",
                    "finding": (
                        "References 14 through 24 and blank page remainder. The entries include "
                        "prior reaction papers and SHELX/SADABS software citations, but no new "
                        "present-paper experiment, reaction arrow, table, footnote condition, "
                        "yield, or spectrum."
                    ),
                    "evidence_ids": evidence(15),
                }
            ],
            "record_ids": [],
            "unresolved": [
                "No reaction record is linked because the page is exclusively bibliography and blank space."
            ],
        },
    ],
}

inventory = {
    "document_id": DOC,
    "pages_reviewed": [13, 14, 15],
    "visual_review": "Original-resolution JPEG render inspected for each page.",
    "new_molecule_count": 0,
    "new_reaction_count": 0,
    "source_context_count": 3,
    "pages": {
        "13": {
            "content": "reference list, items 1(f)-6(i)",
            "reaction_arrows": 0,
            "tables": 0,
            "experimental_records": 0,
        },
        "14": {
            "content": "reference list, items 6(j)-13",
            "reaction_arrows": 0,
            "tables": 0,
            "experimental_records": 0,
        },
        "15": {
            "content": "reference list, items 14-24; lower page blank",
            "reaction_arrows": 0,
            "tables": 0,
            "experimental_records": 0,
        },
    },
}

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model" / "submitted_pages_0013_0015.json").write_text(
    json.dumps(bundle, ensure_ascii=False, indent=2) + "\n"
)
(RUN / "raw_model" / "page_inventory_0013_0015.json").write_text(
    json.dumps(inventory, ensure_ascii=False, indent=2) + "\n"
)
(RUN / "coverage" / "pages_0013_0015.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n"
)

print(
    json.dumps(
        {
            "molecules": len(bundle["molecules"]),
            "reactions": len(bundle["reactions"]),
            "source_context": len(bundle["source_context"]),
            "new_molecules": 0,
            "new_reactions": 0,
            "new_source_context": 3,
        },
        ensure_ascii=False,
    )
)
