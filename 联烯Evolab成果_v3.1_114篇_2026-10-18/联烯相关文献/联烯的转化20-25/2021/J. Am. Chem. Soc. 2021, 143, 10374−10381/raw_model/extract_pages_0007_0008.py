#!/usr/bin/env python3
"""Add the visually reviewed, non-experimental context for PDF pages 7-8.

The two pages contain conclusion/metadata/reference matter only.  Existing
pages 1-6 records are preserved solely to permit an atomic standard merge;
they are not re-extracted or claimed as reviewed by this page-range pass.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_9508ed0ca7d13d60"


def dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads((RUN / "dataset.json").read_text())
template = json.loads((RUN / "raw_model/init_pages_0007_0008.json").read_text())
standard = dataset["standard"]

# Keep the already merged registry/reactions intact.  This scoped pass adds no
# molecule or reaction records because neither page contains a present-paper
# experimental arrow, table row, procedure, or result.
template["documents"] = dataset["documents"]
template["molecules"] = standard["molecules"]
template["reactions"] = standard["reactions"]
template["source_context"] = list(standard.get("source_context", []))

new_context = [
    {
        "context_id": f"{DOC}:page7_back_matter",
        "document_id": DOC,
        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
        "kind": "article_back_matter",
        "note": (
            "Visually reviewed page 7: conclusion tail, Supporting Information "
            "availability statement, author/contact/contribution/competing-interest "
            "information, acknowledgments, and references 1-8. No present-paper "
            "reaction arrow, experimental table, procedure, yield, selectivity, or "
            "standalone characterization result appears on this page. The SI is "
            "mentioned but was not supplied or associated in this workspace."
        ),
    },
    {
        "context_id": f"{DOC}:page8_references",
        "document_id": DOC,
        "evidence_ids": [f"{DOC}:page8:image", f"{DOC}:page8:text"],
        "kind": "reference_list",
        "note": (
            "Visually reviewed page 8: references 9-29 only. Bibliographic titles "
            "mention peptide synthesis, coupling reagents, ynamides, active esters, "
            "and carfilzomib, but these cited works are not present-paper experiments "
            "and were not promoted to molecule or reaction records."
        ),
    },
]
existing_ids = {row.get("context_id") for row in template["source_context"]}
template["source_context"].extend(
    row for row in new_context if row["context_id"] not in existing_ids
)

dump(RUN / "submitted.json", template)
dump(
    RUN / "raw_model/pages_0007_0008_visual_review.json",
    {
        "document_id": DOC,
        "batch_id": "allene-stereo-v31-200-20261017",
        "pages_visually_inspected": [7, 8],
        "new_molecule_records": 0,
        "new_reaction_records": 0,
        "source_context_added": [row["context_id"] for row in new_context],
        "page_findings": {
            "7": [
                "conclusion tail",
                "Supporting Information availability statement",
                "author information, acknowledgments, references 1-8",
                "no present-paper experiment or standalone characterization record",
            ],
            "8": [
                "references 9-29 only",
                "no present-paper experiment or standalone characterization record",
            ],
        },
        "limitations": [
            "The explicitly mentioned Supporting Information was not supplied or associated.",
            "Cited literature was inventoried as bibliography, not re-extracted as current-paper experiments.",
            "Existing pages 1-6 records were preserved but were not re-reviewed in this scoped pass.",
        ],
    },
)
dump(
    RUN / "coverage/pages_0007_0008.json",
    {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 7,
                "state": "attempted",
                "regions": [
                    {
                        "region": "left column, conclusion continuation",
                        "finding": (
                            "Conclusion-level summary of the allenone-mediated amide/peptide "
                            "formation method; no new reaction arrow, experimental condition, "
                            "or product-specific outcome."
                        ),
                        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                    },
                    {
                        "region": "left column, Associated Content / Supporting Information",
                        "finding": (
                            "The article states that SI contains experimental procedures and "
                            "characterization data, but no SI file was supplied or explicitly "
                            "associated; no SI-derived facts were extracted."
                        ),
                        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                    },
                    {
                        "region": "left column lower half and right column upper section",
                        "finding": (
                            "Author information, contributions, competing-interest statement, "
                            "and acknowledgments; administrative metadata, not experiments."
                        ),
                        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                    },
                    {
                        "region": "right column, References 1-8",
                        "finding": (
                            "Bibliographic citations to prior peptide/amide synthesis work; cited "
                            "work was not promoted to present-paper experiment records."
                        ),
                        "evidence_ids": [f"{DOC}:page7:image", f"{DOC}:page7:text"],
                    },
                ],
                "record_ids": [],
                "unresolved": [
                    "Supporting Information is explicitly mentioned but was not supplied or associated.",
                    "No reaction or molecule record is warranted from the administrative and bibliographic material on this page.",
                ],
            },
            {
                "page": 8,
                "state": "attempted",
                "regions": [
                    {
                        "region": "both columns, References 9-29",
                        "finding": (
                            "Reference list only. Titles include automated peptide synthesis, "
                            "coupling reagents, ynamides, active esters, and carfilzomib, but no "
                            "current-paper experimental arrows, table rows, procedures, or results "
                            "are printed on the page."
                        ),
                        "evidence_ids": [f"{DOC}:page8:image", f"{DOC}:page8:text"],
                    }
                ],
                "record_ids": [],
                "unresolved": [
                    "Cited papers were not supplied as source documents and are not treated as current-paper experiments.",
                    "No reaction or molecule record is warranted from this reference-only page.",
                ],
            },
        ],
    },
)

print(
    json.dumps(
        {
            "preserved_molecules": len(template["molecules"]),
            "preserved_reactions": len(template["reactions"]),
            "source_context_total": len(template["source_context"]),
            "new_reactions": 0,
        },
        ensure_ascii=False,
    )
)
