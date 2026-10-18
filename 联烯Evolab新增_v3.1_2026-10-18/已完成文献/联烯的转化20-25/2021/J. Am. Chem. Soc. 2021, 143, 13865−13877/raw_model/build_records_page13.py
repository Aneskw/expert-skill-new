#!/usr/bin/env python3
"""Add source context and exact coverage for original PDF page 13.

Page 13 contains only references 78--86 and the article footer.  It does not
contain a current-paper synthesis experiment, table, scheme, spectrum, or
reaction arrow.  Existing records are preserved solely so standard merge can
update the same persistent dataset atomically; none is claimed as a page-13
result.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_dfedfa74f4eea49f"
EVIDENCE = [f"{DOC}:page13:image", f"{DOC}:page13:text"]


def write(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = copy.deepcopy(json.loads((RUN / "standard_records.json").read_text()))

# Use the current prepared document/evidence registry, which includes both the
# page-13 text and rendered page image.  Do not reconstruct evidence by hand.
bundle["documents"] = copy.deepcopy(dataset["documents"])

context_id = f"{DOC}:page13_references_78_86"
context = {
    "context_id": context_id,
    "document_id": DOC,
    "kind": "reference_list",
    "description": (
        "Original PDF page 13 is the continuation and end of the reference "
        "list (references 78--86). It contains no current-paper reaction "
        "arrow, experimental table, synthetic procedure, characterization "
        "table, or standalone mechanistic experiment. The cited publications "
        "are literature context and are not recast as experiments of this paper."
    ),
    "reference_numbers": list(range(78, 87)),
    "evidence_ids": EVIDENCE,
}
existing = {item.get("context_id") for item in bundle.get("source_context", [])}
if context_id not in existing:
    bundle.setdefault("source_context", []).append(context)

write(RUN / "submitted.json", bundle)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left reference-list column, references 78--86",
                    "finding": (
                        "Reference-list continuation only. References 78--86 cite prior "
                        "synthetic, computational, and catalyst-pocket studies; they are "
                        "not current-paper experiments and no reaction record is created."
                    ),
                    "evidence_ids": EVIDENCE,
                },
                {
                    "region": "remaining page area and footer",
                    "finding": (
                        "Blank page area plus journal/article footer and DOI; no figure, "
                        "table, reaction arrow, footnote experiment, spectrum, or other "
                        "extractable chemical experiment is present."
                    ),
                    "evidence_ids": [f"{DOC}:page13:image"],
                },
            ],
            "record_ids": [],
            "unresolved": [
                "No page-13 experimental record exists to associate; cited papers were not supplied as source documents and were not extracted.",
                "Supporting Information remains unassociated and was not inferred from the article's references or DOI.",
            ],
        }
    ],
}
write(RUN / "coverage/pages_0013_0013.json", coverage)

visual = {
    "document_id": DOC,
    "actual_pages": [13],
    "page_images_inspected": [f"assets/{DOC}/page13.jpg"],
    "regions": [
        {
            "region": "references 78--86",
            "classification": "reference_list",
            "new_reaction_records": 0,
            "existing_reaction_records_updated": 0,
        },
        {
            "region": "blank lower page and footer",
            "classification": "nonexperimental_bibliographic_matter",
            "new_reaction_records": 0,
            "existing_reaction_records_updated": 0,
        },
    ],
    "inventoried_features": {
        "reaction_arrows": 0,
        "tables": 0,
        "figures": 0,
        "experimental_footnotes": 0,
        "spectra_or_peak_tables": 0,
        "current_paper_experiments": 0,
        "cited_references": 9,
    },
    "source_context_id": context_id,
    "evidence_ids": EVIDENCE,
}
write(RUN / "raw_model/visual_transcription_page13.json", visual)

counts = {
    "document_id": DOC,
    "actual_pages": [13],
    "new_molecules_page13": 0,
    "new_reactions_page13": 0,
    "updated_existing_reactions_page13": 0,
    "new_source_context_items_page13": 1 if context_id not in existing else 0,
    "molecules_total_in_submission": len(bundle.get("molecules", [])),
    "reactions_total_in_submission": len(bundle.get("reactions", [])),
    "needs_review_total": sum(
        r.get("review", {}).get("status") == "needs_review"
        for r in bundle.get("reactions", [])
    ),
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
write(RUN / "raw_model/premerge_counts_page13.json", counts)

marker = "## Extension: original PDF page 13"
workflow_path = ROOT / "runs/workflow.md"
workflow = workflow_path.read_text() if workflow_path.exists() else "# Extraction workflow\n"
if marker not in workflow:
    workflow += f"""

{marker}

- Batch: `allene-stereo-v31-200-20261017`; document `{DOC}`; claimed visual scope is exactly original PDF page 13. The supplied full-page JPEG was opened at original detail, and the extracted page text was used only as corroborating evidence.
- Source mapping identifies the main article and leaves SI association unestablished. No SI, cited-paper PDF, or crystallographic file was inferred or associated.
- Page inventory: the upper-left column contains references 78–86; the rest is blank page area and the journal/DOI footer. There are zero figures, tables, reaction arrows, experimental footnotes, spectra/peak tables, and current-paper experiments on this page.
- The cited publications include synthetic and computational topics, but they remain reference-list context. They are not converted into current-paper `literature_example` or `specific_experiment` records because page 13 provides citations rather than source reaction schemes or experimental observations.
- Therefore this increment adds zero molecules and zero reactions, updates zero existing reactions, and adds one `source_context` item plus `coverage/pages_0013_0013.json`. Existing 67 reactions are carried through unchanged solely for atomic v3.1 merge into the persistent dataset; they are not claimed as page-13 results.
- No stereochemical molecule review is applicable to page 13 because no molecule drawing or new molecular identity occurs. Existing molecule stereochemistry fields are not altered.
- Validation and audit are mechanical checks only. Existing machine records remain `needs_review`; semantic and chemical accuracy flags remain false. Existing source gaps and warnings are not filled from the reference list.
"""
    workflow_path.write_text(workflow)

print(json.dumps(counts, ensure_ascii=False, indent=2))
