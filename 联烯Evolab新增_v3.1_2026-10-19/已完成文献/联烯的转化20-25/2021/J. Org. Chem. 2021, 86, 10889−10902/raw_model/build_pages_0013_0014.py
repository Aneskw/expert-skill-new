#!/usr/bin/env python3
"""Extend the v3.1 bundle with the visually reviewed PDF pages 13--14.

These pages are bibliography pages. Page 13 also contains several current-paper
footnotes. Two of those are duplicate evidence for existing records; vague
SI-dependent observations are retained as source context rather than expanded
into invented experiments. Page 14 contains cited literature only.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as sr  # noqa: E402


RUN = ROOT / "runs/current"
DOC = "doc_036d1c8f4ae80aff"
bundle = json.loads((RUN / "standard_records.json").read_text())
init_bundle = json.loads((RUN / "raw_model/pages_0013_0014_init.json").read_text())
bundle["documents"] = copy.deepcopy(init_bundle["documents"])


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def fact(raw, normalized=None, evidence=(), *, method="original_extraction",
         status="pending_review", verification="unreviewed", derivation=None):
    if normalized is None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   status=status, verification=verification, derivation=derivation)


def append_evidence(target: dict, evidence: list[str]) -> None:
    target["evidence_ids"] = list(dict.fromkeys(target.get("evidence_ids", []) + evidence))


def append_fact_evidence(target: dict, evidence: list[str]) -> None:
    if isinstance(target, dict) and "evidence_ids" in target:
        append_evidence(target, evidence)


def append_note_once(reaction: dict, note: dict) -> None:
    normalized = note["normalized_value"]
    if not any(existing.get("normalized_value") == normalized for existing in reaction["notes"]):
        reaction["notes"].append(note)


reactions = {r["reaction_id"]: r for r in bundle["reactions"]}
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}


# Footnote 23 repeats a result already represented by the page-2/page-6
# cyclization record. Only the explicitly repeated byproduct facts receive the
# new evidence; page 13 does not support copying the full reaction conditions.
rid_2c = f"{DOC}:scheme2_1c_2c"
r_2c = reactions[rid_2c]
append_evidence(r_2c, ev(13))
byproduct_outcome = next(o for o in r_2c["outcomes"] if o["outcome_id"] == "byproduct_yield")
for key in ("product_role", "value", "scope"):
    append_fact_evidence(byproduct_outcome[key], ev(13))
byproduct_mid = byproduct_outcome["target_molecule_ids"][0]
byproduct_mol = molecules[byproduct_mid]
if not any(set(ev(13)) <= set(o["evidence_ids"]) for o in byproduct_mol["occurrences"]):
    byproduct_mol["occurrences"].append({
        "document_id": DOC,
        "compound_label": fact("2c′", "2c′", ev(13)),
        "evidence_ids": ev(13),
    })
append_fact_evidence(byproduct_mol["stereochemistry"]["mixture"], ev(13))
append_note_once(r_2c, fact(
    "The six-membered 3-propenyl isoquinoline 2c′ was also isolated in a low 20% yield as a mixture of Z and E isomers.",
    "Page-13 footnote 23 repeats the existing 2c′ byproduct result (20%; Z/E mixture); it is additional evidence, not a new experiment.",
    ev(13),
))


# Footnote 29 is an author-proposed interpretation based on the X-ray structure
# of 2e. It augments the existing proposed-mechanism record and is not promoted
# to a mechanistic experiment.
rid_mech = f"{DOC}:proposed_inner_outer_mechanism"
r_mech = reactions[rid_mech]
append_evidence(r_mech, ev(13))
append_note_once(r_mech, fact(
    "From the X-ray structure of 2e, a tentative interpretation would be that ‘outer’ sphere pathway is perhaps more favored due to the tendency to establish attractive π−π interactions between both aryl rings on nitrogens of 3d−k.",
    "Tentative author interpretation favoring an outer-sphere pathway; no dedicated mechanistic experiment is reported on page 13.",
    ev(13),
))


# Bibliographic pages and sparse SI-dependent footnotes are source context.
# Footnotes 26 and 27 do indicate low-yield/negative observations, but they do
# not identify exact rows, structures, quantities, or the number of experiments;
# the details are explicitly deferred to SI that is not supplied here.
new_contexts = [
    {
        "context_id": f"{DOC}:context_page13_references",
        "document_id": DOC,
        "pages": [13],
        "kind": "bibliographic_context",
        "evidence_ids": ev(13),
        "note": "References 18-38 and literature discussion only; cited-paper reactions are not experiments in this paper.",
    },
    {
        "context_id": f"{DOC}:context_page13_current_footnotes",
        "document_id": DOC,
        "pages": [13],
        "kind": "current_work_footnotes_with_unavailable_SI",
        "evidence_ids": ev(13),
        "note": (
            "Footnotes 24-27 mention other ligand conditions in SI, CCDC 1983304 for 2e, "
            "low yields of 2e with rac-BNP acid, and failure of 1,1-disubstituted or "
            "trisubstituted allenes under standard conditions. No SI file, exact substrate "
            "structures, quantitative results, or countable experimental rows are supplied, "
            "so no new specific experiment is invented."
        ),
    },
    {
        "context_id": f"{DOC}:context_page14_references",
        "document_id": DOC,
        "pages": [14],
        "kind": "bibliographic_context",
        "evidence_ids": ev(14),
        "note": "Continuation of reference 38 and references 39-45 only; no current-paper experiment, scheme, table, figure, or characterization block is present.",
    },
]
existing_context_ids = {c.get("context_id") for c in bundle["source_context"]}
bundle["source_context"].extend(c for c in new_contexts if c["context_id"] not in existing_context_ids)


for reaction in bundle["reactions"]:
    reaction["review"]["status"] = "needs_review"
    reaction["review"]["method"] = "model"
    reaction["review"]["reviewer"] = None


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "both columns, references 18-22, 28, and 30-38",
                    "finding": "Cited literature only; reaction titles describe other papers and were not attributed to the current paper.",
                    "evidence_ids": ev(13),
                },
                {
                    "region": "right column, footnote 23",
                    "finding": "Repeats the current-paper 2c′ byproduct result (20% yield, Z/E mixture); merged as evidence into the existing 1c cyclization record.",
                    "evidence_ids": ev(13),
                },
                {
                    "region": "right column, footnotes 24-27",
                    "finding": "SI referrals, CCDC notice, a nonquantitative low-yield statement for 2e with rac-BNP acid, and unresolved negative scope classes; insufficient to enumerate new experiments without the unavailable SI.",
                    "evidence_ids": ev(13),
                },
                {
                    "region": "right column, footnote 29",
                    "finding": "Tentative author interpretation favoring an outer-sphere pathway from the X-ray structure of 2e; merged into the existing proposed-mechanism record, not treated as a mechanistic experiment.",
                    "evidence_ids": ev(13),
                },
            ],
            "record_ids": [rid_2c, rid_mech],
            "unresolved": [
                "Supporting Information was not supplied or associated, so footnotes 24, 26, and 27 cannot be resolved into exact conditions or distinct experiment counts.",
                "Footnote 26 reports only ‘low chemical yields’ of 2e with rac-BNP acid; no value, row identity, or full conditions are stated.",
                "Footnote 27 does not identify the 1,1-disubstituted or trisubstituted allene structures, number of tested substrates, or individual outcomes.",
                "CCDC 1983304 availability is bibliographic characterization context; the CIF was not provided or used.",
            ],
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper-left bibliography block and continuation from page 13",
                    "finding": "Continuation of reference 38 and references 39-45 only; no current-paper experiment, reaction arrow, table, figure, footnote result, or characterization data.",
                    "evidence_ids": ev(14),
                }
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}


raw_inventory = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "title": "Enantioenriched α-Vinyl 1,4-Benzodiazepines and 1,4-Benzoxazepines via Enantioselective Rhodium-Catalyzed Hydrofunctionalizations of Alkynes and Allenes",
    "pages": [13, 14],
    "visual_inspection": "Both supplied full-page JPEGs were opened and inspected at original detail.",
    "new_molecule_ids": [],
    "new_reaction_ids": [],
    "patched_existing_record_ids": [rid_2c, rid_mech],
    "new_source_context_ids": [c["context_id"] for c in new_contexts],
    "decisions": [
        "Page-13 footnote 23 is duplicate evidence for the existing 2c′ byproduct, not a new experiment.",
        "Page-13 footnote 29 is a proposed mechanistic interpretation, not a mechanistic-study experiment.",
        "Footnotes 26-27 are too underspecified and SI-dependent to support fabricated structures or experiment counts.",
        "All reaction titles in the page-13/14 reference list belong to cited literature.",
    ],
}


(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "raw_model").mkdir(parents=True, exist_ok=True)
sr.write(RUN / "coverage/pages_0013_0014.json", coverage)
sr.write(RUN / "raw_model/pages_0013_0014_extraction.json", raw_inventory)
sr.write(RUN / "submitted.json", bundle)

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": 0,
    "patched_existing": 2,
    "new_source_contexts": 3,
}, ensure_ascii=False, indent=2))
