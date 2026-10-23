#!/usr/bin/env python3
"""Merge source-derived evidence from PDF page 7 into the existing v3.1 bundle."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_65d1cdf1cb09f398"
PAGE7_IMAGE = f"{DOC}:page7:image"
PAGE7_TEXT = f"{DOC}:page7:text"
PAGE7_EVIDENCE = [PAGE7_IMAGE, PAGE7_TEXT]
RXN_3H = f"{DOC}:table2_3h"
RXN_5C = f"{DOC}:table3_5c"


def read(path: Path):
    return json.loads(path.read_text())


def write(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def extend_unique(values: list[str], additions: list[str]) -> None:
    for value in additions:
        if value not in values:
            values.append(value)


bundle_path = RUN / "submitted.json"
bundle = read(bundle_path)
# The current standard-init template is the prepared dataset's authoritative
# documents/evidence block. The prior page-1--6 submission predates the page-7
# rendered-image evidence, so synchronize this block before standard merge.
init_template = read(RUN / "raw_model/init_template_pages_0007_0007.json")
bundle["documents"] = init_template["documents"]
before = {
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
}

# Page 7 maps CCDC deposits to the already registered products. It does not
# restate the absolute configurations; page 4 remains the evidence for those.
for molecule in bundle["molecules"]:
    if molecule["molecule_id"] not in {"m_3h", "m_5c"}:
        continue
    extend_unique(molecule["occurrences"][0]["evidence_ids"], PAGE7_EVIDENCE)
    extend_unique(molecule["occurrences"][0]["compound_label"]["evidence_ids"], PAGE7_EVIDENCE)
    extend_unique(molecule["stereochemistry"]["assignment"]["evidence_ids"], PAGE7_EVIDENCE)
    extend_unique(molecule["identity_note"]["evidence_ids"], PAGE7_EVIDENCE)

for reaction in bundle["reactions"]:
    if reaction["reaction_id"] not in {RXN_3H, RXN_5C}:
        continue
    extend_unique(reaction["evidence_ids"], PAGE7_EVIDENCE)
    for note in reaction["notes"]:
        normalized = note.get("normalized_value")
        if isinstance(normalized, dict) and normalized.get("CCDC") in {"1937754", "1968335"}:
            extend_unique(note["evidence_ids"], PAGE7_EVIDENCE)

context_id = f"{DOC}:page7_references_and_ccdc"
context = {
    "context_id": context_id,
    "type": "crystallographic_and_bibliographic_context",
    "description": (
        "Page 7 contains only the continuation of references [15], reference [16] mapping "
        "CCDC 1937754 to 3h and CCDC 1968335 to 5c, and reference [17]. It contains no "
        "reaction arrow, structure drawing, table row, procedure, or standalone spectrum. "
        "The CCDC availability statement does not establish an associated SI file in this workspace."
    ),
    "linked_reaction_ids": [RXN_3H, RXN_5C],
    "linked_molecule_ids": ["m_3h", "m_5c"],
    "evidence_ids": PAGE7_EVIDENCE,
}
bundle["source_context"] = [
    row for row in bundle["source_context"] if row.get("context_id") != context_id
] + [context]

write(bundle_path, bundle)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left column: continuation of reference [15]",
                    "finding": (
                        "Cited prior literature only; no current-paper experiment, reaction arrow, "
                        "structure, table row, or procedure was present."
                    ),
                    "evidence_ids": PAGE7_EVIDENCE,
                },
                {
                    "region": "reference [16] across both columns",
                    "finding": (
                        "Crystallographic-data notice mapping CCDC 1937754 to product 3h and "
                        "CCDC 1968335 to product 5c. This evidence was merged into the existing "
                        "Table 2 3h and Table 3 5c records; it is not a new experiment and does not "
                        "establish an SI association."
                    ),
                    "evidence_ids": PAGE7_EVIDENCE,
                },
                {
                    "region": "upper-right column: reference [17] and remainder of page",
                    "finding": (
                        "Cited medicinal-chemistry literature followed by journal footer and blank "
                        "page area; no current-paper experiment or independent characterization table."
                    ),
                    "evidence_ids": PAGE7_EVIDENCE,
                },
            ],
            "record_ids": [RXN_3H, RXN_5C],
            "unresolved": [
                "The referenced crystallographic files were not supplied and were not inspected.",
                "No Supporting Information file is explicitly associated in the workspace.",
            ],
        }
    ],
}
write(RUN / "coverage/pages_0007_0007.json", coverage)

raw_update = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "page_range": [7, 7],
    "source_image_inspected": "assets/doc_65d1cdf1cb09f398/page7.jpg",
    "evidence_ids": PAGE7_EVIDENCE,
    "new_reaction_records": [],
    "updated_reaction_records": [RXN_3H, RXN_5C],
    "updated_molecule_records": ["m_3h", "m_5c"],
    "new_source_context": context,
    "decision": (
        "No new synthesis, control, mechanism experiment, table row, reaction arrow, or spectrum "
        "occurs on page 7. References were not converted into current-paper experiments."
    ),
    "si_association": "unestablished",
    "review_status": "needs_review",
}
write(RUN / "raw_model/page_0007_evidence_update.json", raw_update)

workflow_path = ROOT / "runs/workflow.md"
workflow = workflow_path.read_text()
workflow = workflow.replace(
    "- Page 7 was not part of the bounded visual pass and is not claimed as covered.\n",
    "- Page 7 crystallographic files and SI were not supplied; their contents remain unreviewed.\n",
)
marker = "## Batch pass — PDF page 7 of 7"
section = f"""

## Batch pass — PDF page 7 of 7

- Scope: original PDF page 7 only. `page7.jpg` was opened at original detail and checked across both columns and the remaining blank page area. The source mapping confirms `presentation_conversion: false`.
- Inventory: continuation of reference [15], crystallographic-data reference [16], reference [17], journal footer, and blank space. No structure drawing, reaction arrow, table, experimental procedure, or standalone spectrum occurs on this page.
- Reference handling: cited papers [15] and [17] were kept as bibliography and were not converted into current-paper experiments.
- Existing-record merge: reference [16] maps CCDC 1937754 to 3h and CCDC 1968335 to 5c. Page-7 image/text evidence was added to `{RXN_3H}` and `{RXN_5C}` and to molecules `m_3h` and `m_5c`; no duplicate experiment was created.
- SI boundary: the statement that crystallographic data are obtainable from CCDC does not associate a local SI file. SI remains unestablished, and the CCDC files were not inspected.
- Output count for this pass: 0 new reactions, 2 existing reactions updated, 2 existing molecules updated, and 1 new `source_context` item. The full bundle remains at 51 reactions and all machine records remain `needs_review`.
- Coverage: `runs/current/coverage/pages_0007_0007.json`; raw decision trace: `runs/current/raw_model/page_0007_evidence_update.json`.
"""
if marker in workflow:
    workflow = workflow[: workflow.index(marker)].rstrip() + section
else:
    workflow = workflow.rstrip() + section
workflow_path.write_text(workflow.rstrip() + "\n")

after = {
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
}
print(json.dumps({"before": before, "after": after, **raw_update}, ensure_ascii=False))
