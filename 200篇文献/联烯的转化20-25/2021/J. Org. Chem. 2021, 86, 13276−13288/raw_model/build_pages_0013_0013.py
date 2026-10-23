#!/usr/bin/env python3
"""Record the visual review of PDF page 13 for doc_5bda04d0a1e86544.

Page 13 is a bibliography continuation.  It contains no current-paper
reaction experiment, so this incremental extraction adds source context and
coverage only; it deliberately creates no molecule or reaction records.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_5bda04d0a1e86544"
BATCH = "allene-stereo-v31-200-20261017"
RUN = ROOT / "runs/current"
DATASET_PATH = RUN / "dataset.json"
BUNDLE_PATH = RUN / "submitted.json"
PAGE = 13
EVIDENCE = [f"{DOC}:page{PAGE}:image", f"{DOC}:page{PAGE}:text"]
CONTEXT_ID = f"{DOC}:page13_bibliography"


dataset = json.loads(DATASET_PATH.read_text())
bundle = json.loads(BUNDLE_PATH.read_text())

if dataset["documents"][0]["document_id"] != DOC:
    raise RuntimeError("dataset document_id does not match the bounded task")

available_evidence = {
    evidence["evidence_id"]
    for document in dataset["documents"]
    for evidence in document.get("evidence", [])
}
missing_evidence = [evidence_id for evidence_id in EVIDENCE if evidence_id not in available_evidence]
if missing_evidence:
    raise RuntimeError(f"page-13 evidence missing from dataset: {missing_evidence}")

# Refresh the document block from the mechanically prepared dataset so the
# formal bundle includes the supplied page-13 image evidence as well as text.
bundle["documents"] = copy.deepcopy(dataset["documents"])

context = {
    "context_id": CONTEXT_ID,
    "document_id": DOC,
    "kind": "bibliography",
    "evidence_ids": EVIDENCE,
    "note": (
        "PDF page 13 was visually inspected. It is a two-column continuation "
        "of the References (entries 4(c)-17) plus the journal footer. It "
        "contains no current-paper reaction scheme, reaction arrow, table, "
        "experimental procedure, conditions, outcome, or standalone "
        "characterization dataset. Reaction descriptions appearing in cited "
        "paper titles remain bibliography and were not imported as experiments."
    ),
}

existing = {item.get("context_id"): index for index, item in enumerate(bundle["source_context"])}
if CONTEXT_ID in existing:
    bundle["source_context"][existing[CONTEXT_ID]] = context
else:
    bundle["source_context"].append(context)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "full-page two-column References continuation (entries 4(c)-17)",
                    "finding": (
                        "Bibliographic citations and review/article titles only; no current-paper "
                        "experiment, reaction arrow, conditions, or result is present."
                    ),
                    "evidence_ids": EVIDENCE,
                },
                {
                    "region": "right-column references 11-17",
                    "finding": (
                        "Several titles concern allene/allenamide reactions, but they are cited prior "
                        "literature and are not extracted as experiments from this paper."
                    ),
                    "evidence_ids": EVIDENCE,
                },
                {
                    "region": "page footer and article citation line",
                    "finding": "Journal pagination/DOI footer only; no chemical record.",
                    "evidence_ids": EVIDENCE,
                },
            ],
            "record_ids": [],
            "unresolved": [
                "The cited papers' underlying experimental details are outside the supplied source page and were not imported."
            ],
        }
    ],
}

inventory = {
    "document_id": DOC,
    "batch_id": BATCH,
    "visually_inspected_pages": [13],
    "image_asset": f"assets/{DOC}/page13.jpg",
    "image_reviewed": True,
    "page_classification": "bibliography_only",
    "new_molecule_records": 0,
    "new_reaction_records": 0,
    "new_source_context_ids": [CONTEXT_ID],
    "notes": [
        "No current-paper reaction arrow, table row, experimental procedure, or outcome occurs on page 13.",
        "Reaction-bearing titles in references were not treated as source experiments.",
        "No stereochemistry review or standard stereo command was needed because no molecule was registered from this page.",
        "No SI document was supplied or associated.",
    ],
}

st.write(BUNDLE_PATH, bundle)
st.write(RUN / "coverage/pages_0013_0013.json", coverage)
st.write(RUN / "raw_model/pages_0013_0013_extraction.json", inventory)

print(
    json.dumps(
        {
            "page": 13,
            "new_molecules": 0,
            "new_reactions": 0,
            "source_context_id": CONTEXT_ID,
        },
        ensure_ascii=False,
    )
)
