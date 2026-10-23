#!/usr/bin/env python3
"""Record the bounded visual review of original PDF pages 7--8.

The two pages contain administrative matter and bibliography only.  This
addendum deliberately creates no molecule or reaction records and preserves
the existing page-1--6 extraction unchanged.
"""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_e47326ab6eaa3ee5"


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


page7 = f"{DOC}:page7:image"
page8 = f"{DOC}:page8:image"

raw_review = {
    "document_id": DOC,
    "requested_original_pdf_pages": [7, 8],
    "pages_visually_inspected_at_original_detail": [7, 8],
    "new_molecule_records": 0,
    "new_reaction_records": 0,
    "pages": [
        {
            "page": 7,
            "evidence_id": page7,
            "visual_findings": [
                "Author affiliations and the complete-contact-information link occupy the upper left column.",
                "The Notes section states that the authors declare no competing financial interest.",
                "The Acknowledgments section lists funding support and contains no synthetic operation or measured chemical result.",
                "References 1--9 occupy the rest of the page; these are bibliographic citations to background, reviews, and prior examples, not experiments performed in the present paper.",
            ],
            "reaction_arrows": 0,
            "tables": 0,
            "figures": 0,
            "synthetic_experiments": 0,
            "characterization_only_observations": 0,
        },
        {
            "page": 8,
            "evidence_id": page8,
            "visual_findings": [
                "The page contains only the continuation of reference 9 and references 10--14 in the upper-left portion.",
                "The remainder of the page is blank apart from the printed page number.",
                "No reaction arrow, table, figure, experimental procedure, reaction result, or standalone characterization table is present.",
            ],
            "reaction_arrows": 0,
            "tables": 0,
            "figures": 0,
            "synthetic_experiments": 0,
            "characterization_only_observations": 0,
        },
    ],
    "decisions": [
        "No new v3.1 molecule or reaction record was created from pages 7--8.",
        "Bibliographic entries were not promoted to present-paper literature_example records because these pages provide citations only, not depicted reaction families or experimental conditions.",
        "The existing Scheme 1 prior-work records remain supported by their depicted source region on page 2; bibliography pages 7--8 were not added as reaction evidence.",
        "No SI association was inferred from references or publisher links.",
        "No stereochemical claim, selectivity outcome, or formation-bond classification arises from these two pages.",
    ],
}
dump(RUN / "raw_model/pages_0007_0008_visual_review.json", raw_review)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left author affiliations and contact-information block",
                    "finding": "Administrative author and contact information; no reaction, experimental condition, or chemical result.",
                    "evidence_ids": [page7],
                },
                {
                    "region": "left-column Notes and Acknowledgments",
                    "finding": "Competing-interest statement and funding acknowledgments; no synthetic experiment or characterization observation.",
                    "evidence_ids": [page7],
                },
                {
                    "region": "References 1--9 across the lower left and right columns",
                    "finding": "Bibliographic citations to natural products, dearomatization, and domino Heck literature. They are cited prior work and are not recorded as experiments of this paper.",
                    "evidence_ids": [page7],
                },
            ],
            "record_ids": [],
            "unresolved": [
                "The cited external papers were not supplied or inspected; their individual experiments and stereochemistry cannot be extracted from bibliography entries alone."
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left continuation of reference 9 and references 10--14",
                    "finding": "Bibliography continuation only. No present-paper reaction arrow, table row, procedure, outcome, mechanistic experiment, or standalone spectrum/peak table is shown.",
                    "evidence_ids": [page8],
                },
                {
                    "region": "remaining blank page area and footer",
                    "finding": "No additional scientific content beyond the reference list and printed page number.",
                    "evidence_ids": [page8],
                },
            ],
            "record_ids": [],
            "unresolved": [
                "References 10--14 identify prior publications but do not provide enough source evidence here to reconstruct their experiments or stereochemical records."
            ],
        },
    ],
}
dump(RUN / "coverage/pages_0007_0008.json", coverage)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("evidence_id") in {f"{DOC}:page7:text", f"{DOC}:page8:text"}:
        task["state"] = "attempted"
        task["reason"] = (
            "Original rendered page inspected at original detail; administrative/bibliographic "
            "findings and the no-experiment decision are recorded in coverage/pages_0007_0008.json."
        )
dump(tasks_path, tasks)

submitted_path = RUN / "submitted.json"
submitted = json.loads(submitted_path.read_text())
# `inspect` materialized the rendered page-7/page-8 image evidence in the
# prepared dataset.  The merge contract requires the submitted document
# registry to be the exact prepared registry, so carry it over without
# modifying any source metadata or hashes.
dataset = json.loads((RUN / "dataset.json").read_text())
submitted["documents"] = dataset["documents"]
new_context = [
    {
        "context_id": f"{DOC}:context:page7_administrative_bibliography",
        "evidence_ids": [page7],
        "kind": "administrative_and_bibliography",
        "value": "Author/contact information, no-competing-interest statement, acknowledgments, and references 1--9; no present-paper experiment or standalone characterization record.",
    },
    {
        "context_id": f"{DOC}:context:page8_bibliography",
        "evidence_ids": [page8],
        "kind": "bibliography",
        "value": "Continuation of reference 9 and references 10--14 only; no present-paper experiment, reaction result, or standalone characterization record.",
    },
]
by_id = {item.get("context_id"): item for item in submitted.get("source_context", [])}
for item in new_context:
    by_id[item["context_id"]] = item
submitted["source_context"] = list(by_id.values())
dump(submitted_path, submitted)

print(json.dumps({
    "pages": [7, 8],
    "new_molecules": 0,
    "new_reactions": 0,
    "source_context_added": 2,
}, ensure_ascii=False))
