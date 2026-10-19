#!/usr/bin/env python3
"""Merge visual review of source PDF pages 7-8 into the existing v3.1 bundle.

Pages 7-8 add source evidence and coverage, but no new reaction identity:
page 7 continues characterization/results for the existing Fig. 2 scope
experiments 2l-2q, and page 8 is references only.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
RAW = RUN / "raw_model"
DOC = "doc_0414a2e90dea4adf"
E7I = f"{DOC}:page7:image"
E7T = f"{DOC}:page7:text"
E8I = f"{DOC}:page8:image"
E8T = f"{DOC}:page8:text"

sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def append_unique(values: list[str], additions: list[str]) -> None:
    for value in additions:
        if value not in values:
            values.append(value)


def add_image_beside_text(node: object, text_evidence: str, image_evidence: str) -> None:
    """Recursively add the inspected page image wherever its page text is cited."""
    if isinstance(node, dict):
        evidence = node.get("evidence_ids")
        if isinstance(evidence, list) and text_evidence in evidence:
            append_unique(evidence, [image_evidence])
        for value in node.values():
            add_image_beside_text(value, text_evidence, image_evidence)
    elif isinstance(node, list):
        for value in node:
            add_image_beside_text(value, text_evidence, image_evidence)


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = json.loads((RUN / "submitted.json").read_text())

# Restore the complete prepared document/evidence registry. The earlier segment
# bundle stopped its image evidence list at page 6 even though prepare rendered
# all eight pages.
bundle["documents"] = dataset["documents"]

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}

# Page 7 explicitly names and characterizes products 2m/3m through 2q/3q.
# Its image is added beside the already retained page-text evidence. Candidate
# structures remain machine records (needs_review); this does not promote them.
for letter in "mnopq":
    for number in (2, 3):
        molecule_id = f"{DOC}:mol:{number}{letter}"
        add_image_beside_text(molecules[molecule_id], E7T, E7I)
    reaction_id = f"{DOC}:scope_2{letter}"
    add_image_beside_text(reactions[reaction_id], E7T, E7I)
    append_unique(reactions[reaction_id]["evidence_ids"], [E7I])
    # The previous segment script had a stale structure-gap note even though
    # candidate SMILES were present. Remove only that contradictory note.
    reactions[reaction_id]["review"]["unresolved"] = [
        item
        for item in reactions[reaction_id]["review"]["unresolved"]
        if not item.startswith("Complete candidate structures for the long-chain/benzyloxy")
    ]

# The top-left of page 7 is the continuation of 2l/3l characterization from
# page 6. It maps to the same scope experiment and is not a new synthesis.
scope_2l_id = f"{DOC}:scope_2l"
append_unique(reactions[scope_2l_id]["evidence_ids"], [E7I])

context_ids = {
    item.get("context_id") for item in bundle["source_context"] if isinstance(item, dict)
}
new_context = [
    {
        "context_id": f"{DOC}:context:characterization_page_7",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": [E7I, E7T],
        "record_ids": [f"{DOC}:scope_2{x}" for x in "lmnopq"],
        "text": (
            "Page 7 continues the 2l/3l characterization and gives product paragraphs "
            "for the existing Fig. 2 experiments 2m/3m through 2q/3q. NMR and MS peak "
            "lists are characterization context, not additional syntheses."
        ),
    },
    {
        "context_id": f"{DOC}:context:back_matter_page_7",
        "document_id": DOC,
        "kind": "back_matter",
        "evidence_ids": [E7I, E7T],
        "text": "Conflicts of interest, acknowledgements, and the beginning of the reference list; no new experiment.",
    },
    {
        "context_id": f"{DOC}:context:references_page_8",
        "document_id": DOC,
        "kind": "references",
        "evidence_ids": [E8I, E8T],
        "text": "Reference list only; no reaction arrow, experimental table, procedure, or new result is shown.",
    },
]
bundle["source_context"].extend(
    item for item in new_context if item["context_id"] not in context_ids
)

visual_inventory = {
    "document_id": DOC,
    "pages": [
        {
            "page": 7,
            "regions": [
                {
                    "region": "upper-left continuation above the 2m heading",
                    "finding": "Continuation of 2l/3l NMR and EI-MS characterization; same Fig. 2 scope experiment, not a new synthesis.",
                    "evidence_ids": [E7I],
                },
                {
                    "region": "left column product paragraphs 2m/3m, 2n/3n, 2o/3o, and 2p/3p",
                    "finding": "Existing substrate-scope experiments with isolated yields, regioisomer ratios, and characterization; no new reaction conditions are introduced.",
                    "evidence_ids": [E7I],
                },
                {
                    "region": "lower-left to upper-right product paragraph 2q/3q",
                    "finding": "Existing substrate-scope experiment with 75% yield and 83:17 2q:3q ratio, followed by characterization.",
                    "evidence_ids": [E7I],
                },
                {
                    "region": "right column back matter and references 1-4",
                    "finding": "Conflicts of interest, acknowledgements, and bibliography; no experimental record.",
                    "evidence_ids": [E7I],
                },
            ],
        },
        {
            "page": 8,
            "regions": [
                {
                    "region": "full two-column page",
                    "finding": "Continuation and completion of references 4-17 only; no figure, table, reaction arrow, footnote, procedure, or experimental result.",
                    "evidence_ids": [E8I],
                }
            ],
        },
    ],
}

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": visual_inventory["pages"][0]["regions"],
            "record_ids": [f"{DOC}:scope_2{x}" for x in "lmnopq"],
            "unresolved": [
                "The source does not report the measurement method for the 2m:3m through 2q:3q regioisomer ratios.",
                "The page-7 value 72:18 for 2n:3n conflicts with Fig. 2 (72:28); the source values are retained without repair.",
                "Page-7 yields/ratios for 2o and 2p conflict with the corresponding Fig. 2 labels; both source alternatives are retained.",
                "Absolute allene scale and numeric solvent volumes remain unreported in the supplied main text; no SI was associated.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": visual_inventory["pages"][1]["regions"],
            "record_ids": [],
            "unresolved": [
                "Bibliographic citations alone do not provide enough source detail to create literature-example experiments."
            ],
        },
    ],
}

raw_review = {
    "document_id": DOC,
    "review_scope": {"pages": [7, 8], "evidence_ids": [E7I, E8I]},
    "page_7": {
        "existing_record_results": {
            "2m/3m": {"isolated_yield": "34%", "ratio": "93:7"},
            "2n/3n": {"isolated_yield": "53%", "ratio": "72:18"},
            "2o/3o": {"isolated_yield": "48%", "ratio": "74:26"},
            "2p/3p": {"isolated_yield": "50%", "ratio": "84:16"},
            "2q/3q": {"isolated_yield": "75%", "ratio": "83:17"},
        },
        "new_reaction_records": [],
        "notes": [
            "The top-left text is a continuation of 2l/3l characterization.",
            "NMR/MS peak lists were not converted into synthesis records.",
            "Conflicting Fig. 2 versus experimental-paragraph values remain alternatives.",
        ],
    },
    "page_8": {"new_reaction_records": [], "finding": "references only"},
}

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and task["evidence_id"] in {E7T, E8T}:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; decisions mapped in coverage/pages_0007_0008.json."

write_json(RUN / "submitted.json", bundle)
write_json(RAW / "visual_inventory_pages_0007_0008.json", visual_inventory)
write_json(RAW / "page7_8_visual_review.json", raw_review)
write_json(
    RAW / "source_mentions_pages_0007_0008.json",
    st.source_mentions(dataset, ["2l", "2m", "2n", "2o", "2p", "2q"], context=260),
)
write_json(RUN / "coverage/pages_0007_0008.json", coverage)
write_json(RUN / "tasks.json", tasks)

print(
    json.dumps(
        {
            "pages": [7, 8],
            "new_reactions": 0,
            "existing_reactions_enriched": 6,
            "coverage_records": len(coverage["pages"][0]["record_ids"]),
        },
        ensure_ascii=False,
    )
)
