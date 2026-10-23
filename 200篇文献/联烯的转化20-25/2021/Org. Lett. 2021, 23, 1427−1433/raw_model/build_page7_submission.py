#!/usr/bin/env python3
"""Page-7 evidence delta for doc_26447b9b5c1dd197.

The page contains bibliography plus endnotes.  It creates no new reaction,
but endnote 17 supplies additional evidence for two existing Scheme 4
records.  The delta is initialized from the standard v3.1 template and then
merged into the accumulated dataset.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


RUN = ROOT / "runs/current"
DOC = "doc_26447b9b5c1dd197"
P3 = f"{DOC}:page3:image"
T3 = f"{DOC}:page3:text"
P7 = f"{DOC}:page7:image"
T7 = f"{DOC}:page7:text"


def add_evidence(container: dict, *evidence_ids: str) -> None:
    for evidence_id in evidence_ids:
        if evidence_id not in container["evidence_ids"]:
            container["evidence_ids"].append(evidence_id)


def inferred_status(value: str, evidence_ids: tuple[str, ...], reason: str) -> dict:
    return st.fact(
        raw=None,
        normalized=value,
        evidence=evidence_ids,
        method="inference",
        derivation=reason,
    )


template = st.read(RUN / "raw_model/page7_init.json")
prior = st.read(RUN / "standard_records.json")
reaction_index = {r["reaction_id"]: r for r in prior["reactions"]}
molecule_index = {m["molecule_id"]: m for m in prior["molecules"]}


# Endnote 17 clarifies that the unstable N,O-aminal derived from 1b was not
# isolated and was used after a basic workup.  Correct the prior generic-scale
# inheritance rather than treating the footnote as a new experiment.
rid_4ba = f"{DOC}:scheme4_4ba"
r4ba = copy.deepcopy(reaction_index[rid_4ba])
add_evidence(r4ba, P7, T7)

one_b_id = f"{DOC}:mol:1b"
one_b = st.new_molecule(one_b_id, DOC, (T3, P7, T7))
one_b["identity_status"] = "source_identified"
one_b["occurrences"][0]["compound_label"] = st.fact("1b", "1b", (T3, P7, T7))
one_b["names"]["original"] = st.fact("1b", "allenamide 1b", (T3, T7))
one_b["names"]["english"] = st.fact(evidence=(T3, T7))
one_b["names"]["abbreviations"] = st.fact(["1b"], ["1b"], (T3, T7))
one_b["names"]["aliases"] = st.fact(evidence=(T3, T7))
one_b["cas_numbers"] = st.fact(evidence=(T3, T7))
one_b["representation_type"] = "unresolved"
one_b["smiles"] = st.fact(evidence=(T3, P7, T7))
one_b["generic_template"] = None
reason_1b = (
    "Page 7 identifies 1b only by compound label and as the precursor of an "
    "unstable N,O-aminal; no page-7 structure is available for safe encoding."
)
one_b["stereochemistry"]["chirality_status"] = inferred_status("unresolved", (T3, P7, T7), reason_1b)
one_b["stereochemistry"]["geometry_status"] = inferred_status("unresolved", (T3, P7, T7), reason_1b)
one_b["stereochemistry"]["representation_status"] = inferred_status("unresolved", (T3, P7, T7), reason_1b)
one_b["identity_note"] = st.fact(
    "obtained from 1b",
    "Label-level source identity only; structure and stereochemistry remain unresolved.",
    (T3, P7, T7),
    method="inference",
    derivation="The source identifies the precursor label but page 7 does not provide a structure or CAS number.",
)

starting = st.new_participant("starting_allenamide", one_b_id, "reactant", (T3, P7, T7))
r4ba["participants"] = [p for p in r4ba["participants"] if p["participant_id"] != "starting_allenamide"]
r4ba["participants"].insert(0, starting)

substrate = next(p for p in r4ba["participants"] if p["participant_id"] == "substrate")
substrate["charges"]["amount"] = st.fact(evidence=(P3, P7, T7))
r4ba["reaction_scale"] = st.fact(evidence=(P3, P7, T7))

preformation = st.new_step("intermediate_formation", (T3, P7, T7), "reaction")
preformation["participant_ids"] = ["starting_allenamide", "substrate"]
preformation["operation"] = st.fact(
    "N-tosyl N,O-aminals (obtained from 1b and 1l)",
    "formation of the N-tosyl N,O-aminal intermediate from 1b; detailed preparation unreported in the supplied main paper",
    (T3, P7, T7),
)
basic_workup = st.new_step("basic_workup", (P7, T7), "workup")
basic_workup["participant_ids"] = ["substrate"]
basic_workup["operation"] = st.fact("after basic workup", "basic workup", (P7, T7))
r4ba["steps"] = [
    s for s in r4ba["steps"]
    if s["step_id"] not in {"intermediate_formation", "basic_workup"}
]
r4ba["steps"] = [preformation, basic_workup] + r4ba["steps"]

footnote_17 = (
    "The halobenzoyloxylation products N-tosyl N,O-aminals (obtained from 1b and 1l) "
    "were directly used for alkynylation reaction without isolation (after basic workup) "
    "because of the instability (of aminal from 1b) and formation of an isomeric mixture (from 1l)."
)
r4ba["procedure"]["modifications"] = st.fact(
    footnote_17,
    {
        "intermediate_origin": "N-tosyl N,O-aminal obtained from 1b",
        "isolation": "not isolated",
        "pretreatment": "basic workup",
        "reason": "instability of the aminal from 1b",
    },
    (P7, T7),
)
r4ba["procedure"]["pretreatment"] = st.fact("after basic workup", "basic workup", (P7, T7))

yield_4ba = next(o for o in r4ba["outcomes"] if o["outcome_id"] == "yield")
yield_4ba["yield_type"] = st.fact(
    "yield mentioned is the overall yield",
    "overall yield",
    (P3,),
)
yield_4ba["scope"] = st.fact(
    {
        "scheme_label": "4ba, 73% (footnote b)",
        "endnote_17": footnote_17,
    },
    {
        "start_molecule_id": one_b_id,
        "intermediate_molecule_id": substrate["molecule_id"],
        "end_molecule_id": yield_4ba["target_molecule_ids"][0],
    },
    (P3, P7, T7),
)
instability_note = st.fact(
    "because of the instability (of aminal from 1b)",
    "The N,O-aminal from 1b was unstable and therefore not isolated before alkynylation.",
    (P7, T7),
)
r4ba["notes"] = [n for n in r4ba["notes"] if n["raw_value"] != instability_note["raw_value"]]
r4ba["notes"].append(instability_note)
unresolved_4ba = "Endnote 17 establishes an unisolated intermediate and basic workup, but the supplied main paper does not report the starting 1b scale or the basic-workup reagent, duration, temperature, or operations."
r4ba["review"]["unresolved"] = [u for u in r4ba["review"]["unresolved"] if u != unresolved_4ba]
r4ba["review"]["unresolved"].append(unresolved_4ba)


# The same endnote clarifies the existing 1l -> 6 sequence.  It does not add a
# second reaction and does not resolve the intermediate mixture composition.
rid_1l = f"{DOC}:scheme4_1l_elimination_6"
r1l = copy.deepcopy(reaction_index[rid_1l])
add_evidence(r1l, P7, T7)
r1l["procedure"]["modifications"] = st.fact(
    footnote_17,
    {
        "intermediate_origin": "N-tosyl N,O-aminal obtained from 1l",
        "isolation": "not isolated",
        "pretreatment": "basic workup",
        "intermediate_state": "isomeric mixture",
    },
    (P7, T7),
)
r1l["procedure"]["pretreatment"] = st.fact("after basic workup", "basic workup", (P7, T7))
workup = next(s for s in r1l["steps"] if s["step_id"] == "basic_workup")
workup["operation"] = st.fact("after basic workup", "basic workup; reagent and detailed operation unreported", (P3, P7, T7))
mixture_note = st.fact(
    "formation of an isomeric mixture (from 1l)",
    "The N,O-aminal intermediate derived from 1l was an isomeric mixture.",
    (P7, T7),
)
r1l["notes"] = [n for n in r1l["notes"] if n["raw_value"] != mixture_note["raw_value"]]
r1l["notes"].append(mixture_note)
unresolved_1l = [
    "Endnote 17 does not report the identities, configurations, ratio, or measurement method for the 1l-derived N,O-aminal isomeric mixture.",
    "The source does not explicitly label the printed 42% as an overall yield from 1l versus a post-workup alkynylation-stage yield.",
]
r1l["review"]["unresolved"] = [u for u in r1l["review"]["unresolved"] if u not in unresolved_1l]
r1l["review"]["unresolved"].extend(unresolved_1l)


template["reactions"] = [r4ba, r1l]
needed_molecule_ids = {
    p["molecule_id"]
    for reaction in template["reactions"]
    for p in reaction["participants"]
}
template["molecules"] = [
    one_b if molecule_id == one_b_id else copy.deepcopy(molecule_index[molecule_id])
    for molecule_id in sorted(needed_molecule_ids)
]
template["source_context"] = [
    {
        "context_id": f"{DOC}:context:page7_bibliography",
        "document_id": DOC,
        "kind": "bibliography_only",
        "note": "Upper page region continues references 9 and lists references 10-16; reaction language belongs to cited prior literature and is not a current-paper experiment.",
        "evidence_ids": [P7, T7],
    },
    {
        "context_id": f"{DOC}:context:page7_endnote15",
        "document_id": DOC,
        "kind": "source_availability",
        "note": "Endnote 15 refers preparation details to Supporting Information. No SI was supplied or explicitly associated, so no procedure was inferred.",
        "evidence_ids": [P7, T7],
    },
    {
        "context_id": f"{DOC}:context:page7_endnote17",
        "document_id": DOC,
        "kind": "current_experiment_clarification",
        "note": "Endnote 17 clarifies the unisolated, basic-worked-up N,O-aminal intermediates derived from 1b and 1l; it is linked to the existing 4ba and product-6 records rather than counted as a new experiment.",
        "record_ids": [rid_4ba, rid_1l],
        "evidence_ids": [P7, T7],
    },
    {
        "context_id": f"{DOC}:context:page7_endnotes18_20",
        "document_id": DOC,
        "kind": "computation_and_mechanistic_discussion",
        "note": "Endnotes 18-20 state computational levels, define conformational terminology, and discuss alkynylzinc aggregation limitations; they do not report new observed reactions.",
        "evidence_ids": [P7, T7],
    },
]


st.write(RUN / "submitted.json", template)
st.write(
    RUN / "raw_model/page7_extraction_inventory.json",
    {
        "document_id": DOC,
        "scope_pages": [7],
        "new_reaction_count": 0,
        "updated_reaction_ids": [rid_4ba, rid_1l],
        "new_molecule_ids": [one_b_id],
        "source_context_ids": [c["context_id"] for c in template["source_context"]],
        "chemical_accuracy_verified": False,
    },
)
print(json.dumps({
    "bundle_molecules": len(template["molecules"]),
    "new_reactions": 0,
    "updated_reactions": [rid_4ba, rid_1l],
    "source_context": len(template["source_context"]),
}, ensure_ascii=False, indent=2))
