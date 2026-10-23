#!/usr/bin/env python3
"""Author the bounded PDF-page-7 v3.1 delta after visual inspection."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_448eb03f81b6e3bb"
PAGE4 = [f"{DOC}:page4:image", f"{DOC}:page4:text"]
PAGE7 = [f"{DOC}:page7:image", f"{DOC}:page7:text"]


def write(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def fact(raw=None, normalized=None, evidence=(), *, method="original_extraction",
         derivation=None, status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def qty(number, unit):
    return {"number": number, "unit": unit}


def extend_unique(values, additions):
    for value in additions:
        if value not in values:
            values.append(value)


old = json.loads((RUN / "standard_records.json").read_text())
old_molecules = {m["molecule_id"]: m for m in old["molecules"]}
old_reactions = {r["reaction_id"]: r for r in old["reactions"]}
bundle = json.loads((RUN / "submitted.json").read_text())


def copy_molecule(key: str):
    molecule_id = f"{DOC}:mol:{key}"
    molecule = copy.deepcopy(old_molecules[molecule_id])
    bundle["molecules"].append(molecule)
    return molecule


needed = [
    "1b", "4a", "Na2CO3", "PPh3", "toluene", "5aprac",
    "1a", "4g", "K2CO3", "Cs2CO3", "P5", "5g",
]
molecules = {key: copy_molecule(key) for key in needed}


def add_page7_occurrence(molecule, label):
    occurrence = next(o for o in molecule["occurrences"] if o["document_id"] == DOC)
    extend_unique(occurrence["evidence_ids"], PAGE7)
    extend_unique(occurrence["compound_label"]["evidence_ids"], PAGE7)
    # Page 7 repeats only the source label/name; structure and stereochemical
    # drawing evidence remain tied to their original figure pages.
    if occurrence["compound_label"]["raw_value"] is None:
        occurrence["compound_label"] = fact(label, label, PAGE7)


for key in ("4g", "K2CO3", "Cs2CO3", "5g"):
    add_page7_occurrence(molecules[key], key)

# Footnote 13 supplies the exact X-ray assignment/deposition evidence for 5g.
m5g = molecules["5g"]
m5g["stereochemistry"]["assignment"] = fact(
    {
        "assignment": "The absolute configuration of 5g was assigned based on X-ray crystallographic analysis.",
        "deposition": "CCDC 1970725 contains the supplementary crystallographic data for this paper."
    },
    {
        "absolute_configuration_status": "assigned by X-ray crystallographic analysis",
        "ccdc": "1970725",
        "descriptors": "not stated in footnote 13"
    },
    [f"{DOC}:page3:text", *PAGE7],
)
extend_unique(m5g["stereochemistry"]["chirality_status"]["evidence_ids"], PAGE7)

# Footnote 11 names isoindigo but supplies no drawing or stereochemical detail.
iso = st.new_molecule(f"{DOC}:mol:isoindigo", DOC, PAGE7)
iso["identity_status"] = "source_identified"
iso["occurrences"][0]["compound_label"] = fact("isoindigo", "isoindigo", PAGE7)
iso["names"]["original"] = fact("isoindigo", "isoindigo", PAGE7)
iso["names"]["english"] = fact("isoindigo", "isoindigo", PAGE7)
iso["identity_note"] = fact(
    "Footnote 11 names isoindigo but does not display its structure or specify its stereochemical/geometric form.",
    "Source-named additional product; connectivity and stereochemical/geometric identity unresolved from page 7.",
    PAGE7,
)
for field_name in ("chirality_status", "geometry_status", "representation_status"):
    iso["stereochemistry"][field_name] = fact(
        "Footnote 11 provides only the name isoindigo.",
        "unresolved",
        PAGE7,
        method="inference",
        derivation="No structure drawing, configuration, E/Z assignment, or sample composition is provided on page 7; no stereochemical absence or racemic state was inferred.",
    )
bundle["molecules"].append(iso)

# Augment the already-existing Figure 4 1b experiment with its footnote-11
# additional-product observation. This is the same experiment, not a new row.
fig4 = copy.deepcopy(old_reactions[f"{DOC}:fig4_1b"])
extend_unique(fig4["evidence_ids"], PAGE7)
iso_participant = st.new_participant(
    "isoindigo_product", iso["molecule_id"], "product", PAGE7
)
iso_participant["product_role"] = fact(
    "the formation of isoindigo was observed",
    "additional observed product",
    PAGE7,
)
fig4["participants"].append(iso_participant)
if fig4["steps"]:
    extend_unique(fig4["steps"][0]["participant_ids"], ["isoindigo_product"])
iso_outcome = st.new_outcome(
    "isoindigo_observation", "other", [iso["molecule_id"]], PAGE7
)
iso_outcome["product_role"] = fact(
    "the formation of isoindigo was observed",
    "additional observed product",
    PAGE7,
)
iso_outcome["value"] = fact(
    "the formation of isoindigo was observed", "observed", PAGE7
)
iso_outcome["scope"] = fact(
    "the formation of isoindigo was observed",
    {"metric": "qualitative product observation", "amount": "not reported"},
    PAGE7,
)
fig4["outcomes"].append(iso_outcome)
fig4["notes"].append(fact(
    "presumably resulted from bisoxindole substrate",
    "The authors tentatively attribute isoindigo formation to the bisoxindole substrate.",
    PAGE7,
))
fig4["classification"]["limitations"].append(
    "The main annulation is partially classified; the separately observed isoindigo-formation pathway is not classified because page 7 gives no product drawing or atom mapping."
)
fig4["review"]["unresolved"] = [
    item for item in fig4["review"]["unresolved"]
    if item != "Full structures, atom mapping, and SI-only operational details require expert/source review."
] + [
    "Full main-product structures, atom mapping, and SI-only operational details require expert/source review.",
    "Isoindigo amount/yield, measurement method, exact structure/geometric form, and formation pathway are unreported; the substrate-origin statement is explicitly tentative."
]
bundle["reactions"].append(fig4)

# Include the existing general scheme so the two new experiment records can
# retain a valid general_scheme_id in the independently validated delta bundle.
generic_id = f"{DOC}:generic_asymmetric_42"
bundle["reactions"].append(copy.deepcopy(old_reactions[generic_id]))


def inherited(raw, normalized):
    return fact(
        raw,
        normalized,
        [*PAGE4, *PAGE7],
        method="inherited_general_condition",
        derivation="Inherited from the Figure 5 general conditions because footnote 12 identifies a base-variation experiment for the same annulation of substrate 4g; the changed base identity is taken directly from footnote 12. No unreported reaction time or base dose was inherited.",
    )


def add_participant(reaction, participant_id, molecule_key, role, amount=None,
                    volume=None, product_role=None):
    participant = st.new_participant(
        participant_id, molecules[molecule_key]["molecule_id"], role,
        [*PAGE4, *PAGE7]
    )
    if amount is not None:
        participant["charges"]["amount"] = inherited(amount[0], amount[1])
    if volume is not None:
        participant["charges"]["volume"] = inherited(volume[0], volume[1])
    if product_role is not None:
        participant["product_role"] = fact(product_role, product_role, PAGE7)
    reaction["participants"].append(participant)


def base_variant(base_key: str):
    rid = f"{DOC}:footnote12_4g_{base_key.lower()}"
    reaction = st.new_reaction(rid, DOC, [*PAGE4, *PAGE7])
    reaction["record_type"] = "specific_experiment"
    reaction["experiment_categories"] = fact(["control"], ["control"], PAGE7)
    reaction["location"] = fact(
        f"Footnote 12, 4g annulation with {base_key}",
        f"Footnote 12, 4g annulation with {base_key}",
        PAGE7,
    )
    reaction["general_scheme_ids"] = [generic_id]
    reaction["reaction_name"] = fact(
        f"{base_key} as the base for the annulation of substrate 4g",
        "Base-variation [4+2] annulation of substrate 4g",
        PAGE7,
    )
    reaction["transformation_type"] = fact(
        None,
        "enantioselective [4+2] annulation",
        PAGE7,
        method="inference",
        derivation="Normalized from footnote 12's explicit '[4 + 2] product 5g' statement and its Figure 5 context.",
    )

    add_participant(reaction, "allene", "1a", "reactant",
                    amount=("Figure 5: 1 (0.1 mmol)", qty(0.1, "mmol")))
    add_participant(reaction, "bisox", "4g", "reactant",
                    amount=("Figure 5: 4 (0.12 mmol)", qty(0.12, "mmol")))
    add_participant(reaction, "base", base_key, "reagent")
    add_participant(reaction, "cat", "P5", "catalyst",
                    amount=("Figure 5: P5 (0.02 mmol)", qty(0.02, "mmol")))
    add_participant(reaction, "solvent", "toluene", "solvent",
                    volume=("Figure 5: toluene (2 mL)", qty(2, "mL")))
    add_participant(reaction, "product", "5g", "product",
                    product_role="[4 + 2] product 5g")

    reaction["reaction_scale"] = inherited(
        "Figure 5: 1 (0.1 mmol)",
        {"reference_molecule_id": molecules["1a"]["molecule_id"],
         "amount": qty(0.1, "mmol")},
    )
    reaction["solvent_conditions"]["composition"] = inherited(
        "Figure 5: in toluene", ["toluene"]
    )
    reaction["solvent_conditions"]["total_volume"] = inherited(
        "Figure 5: toluene (2 mL)", qty(2, "mL")
    )

    step = st.new_step("annulation", [*PAGE4, *PAGE7])
    step["participant_ids"] = [
        "allene", "bisox", "base", "cat", "solvent", "product"
    ]
    step["operation"] = fact(
        f"employed {base_key} as the base for the annulation of substrate 4g",
        f"Conduct the 4g annulation with {base_key} as base",
        PAGE7,
    )
    step["temperature"] = inherited(
        "Figure 5: at room temperature", "r.t."
    )
    # Footnote 12 and the Figure 5 caption do not give a general time for 5g.
    reaction["steps"].append(step)

    rr = st.new_outcome("regioisomer_ratio", "selectivity",
                        [molecules["5g"]["molecule_id"]], PAGE7)
    rr["product_role"] = fact(
        "the [4 + 2] product 5g", "main reported product", PAGE7
    )
    rr["value"] = fact(
        "rr ratio of 4 to 1",
        {"ratio": "4:1", "qualifier": "exact"},
        PAGE7,
    )
    rr["scope"] = fact(
        "rr ratio of 4 to 1",
        {"metric": "regioisomer ratio",
         "isomer_order": "not restated in footnote 12"},
        PAGE7,
    )
    reaction["outcomes"].append(rr)

    # The same labelled substrate/product transformation as Figure 5 entry 5g
    # is defensibly C-C/C-C at element level, but complete mapping is absent.
    reaction["classification"] = copy.deepcopy(
        old_reactions[f"{DOC}:fig5_5g"]["classification"]
    )
    reaction["classification"]["limitations"] = [
        "Partial element-level classification inherited from the displayed Figure 5 4g-to-5g transformation; footnote 12 supplies no atom mapping, coproduct balance, or complete structures."
    ]
    reaction["notes"] = [fact(
        "Utilization of different bases did not have effects on the regioselectivity of the reaction.",
        "Author observation: changing the base to K2CO3 or Cs2CO3 did not change the reported 4:1 rr.",
        PAGE7,
    )]
    reaction["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [
            "Footnote 12 does not report the base amount, reaction time, isolated yield, ee, or rr measurement method.",
            "The 4:1 regioisomer order is not restated in footnote 12, and complete structures/atom mapping remain unresolved.",
            "Figure 5 conditions other than the explicitly changed base are inherited as contextual general conditions and remain needs_review."
        ],
    }
    reaction["completeness"] = "partial"
    bundle["reactions"].append(reaction)
    return rid


new_ids = [base_variant("K2CO3"), base_variant("Cs2CO3")]

bundle["source_context"].append({
    "context_id": f"{DOC}:page7_reference_inventory",
    "document_id": DOC,
    "kind": "reference_section_inventory",
    "note": "Page 7 otherwise contains continuations of references 9-10 and references 14-18. These cited works were not converted into experiments performed in this article. Footnotes 11-13 were handled separately as current-paper evidence.",
    "evidence_ids": PAGE7,
})

write(RUN / "submitted.json", bundle)

raw = {
    "document_id": DOC,
    "pages": [7],
    "semantic_accuracy_verified": False,
    "visual_review": {
        "asset": f"runs/current/assets/{DOC}/page7.jpg",
        "regions": [
            {
                "region": "left column, references 9(j)-10(e)",
                "finding": "Cited prior literature only; no experiment performed in the current article is displayed."
            },
            {
                "region": "left column, footnote 11",
                "finding": "Additional isoindigo formation observation for the existing Figure 4 experiment using allenoate 1b."
            },
            {
                "region": "left column, footnote 12",
                "finding": "Two distinct base-variation experiments for 4g, one with K2CO3 and one with Cs2CO3; each reports 5g at rr 4:1."
            },
            {
                "region": "left column, footnote 13",
                "finding": "5g absolute-configuration assignment by X-ray analysis and deposition CCDC 1970725; stereochemical evidence, not a new reaction."
            },
            {
                "region": "lower left and upper right, references 14-18",
                "finding": "Cited literature only; no current-paper experiment."
            }
        ]
    },
    "new_reaction_ids": new_ids,
    "augmented_reaction_ids": [f"{DOC}:fig4_1b"],
    "stereochemical_evidence_molecule_ids": [molecules["5g"]["molecule_id"]],
    "limitations": [
        "No Supporting Information was supplied or associated.",
        "Footnote 12 omits base charges, reaction times, yield, ee, and rr measurement method.",
        "Page 7 gives no structure drawing for isoindigo and no explicit R/S descriptors for 5g."
    ]
}
write(RUN / "raw_model/pages_0007_0007_extraction.json", raw)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [{
        "page": 7,
        "state": "attempted",
        "regions": [
            {
                "region": "left column, references 9(j)-10(e)",
                "finding": "Reference-list continuation describing cited prior work; no scheme, table, reaction arrow, or current-paper experiment is present.",
                "evidence_ids": PAGE7,
            },
            {
                "region": "left column, footnote 11",
                "finding": "Current-paper qualitative observation: when allenic ester 1b was used in the Figure 4 feasibility experiment, isoindigo formation was observed and tentatively attributed by the authors to bisoxindole substrate. Merged into the existing experiment rather than duplicated.",
                "evidence_ids": PAGE7,
            },
            {
                "region": "left column, footnote 12",
                "finding": "Two distinct current-paper base-variation controls were identified: annulation of 4g with K2CO3 and with Cs2CO3; both report [4+2] product 5g with rr 4:1.",
                "evidence_ids": PAGE7,
            },
            {
                "region": "left column, footnote 13",
                "finding": "Stereochemical attribution for molecule 5g: absolute configuration assigned by X-ray crystallographic analysis; CCDC 1970725. This is molecule-assignment evidence, not a separate synthesis experiment.",
                "evidence_ids": PAGE7,
            },
            {
                "region": "lower left and upper right, references 14-18",
                "finding": "Bibliographic references only; cited experiments were not converted into experiments of this paper.",
                "evidence_ids": PAGE7,
            },
        ],
        "record_ids": [
            f"{DOC}:fig4_1b",
            f"{DOC}:fig5_5g",
            *new_ids,
        ],
        "unresolved": [
            "Footnote 11 does not give an isoindigo amount/yield, measurement method, structure drawing, or exact geometric form.",
            "Footnote 12 does not report base amount, reaction time, yield, ee, rr measurement method, or explicitly restate the regioisomer ordering.",
            "Footnote 13 does not state explicit absolute-configuration descriptors for 5g on this page; the structure remains not encoded as isomeric SMILES."
        ],
    }]
}
write(RUN / "coverage/pages_0007_0007.json", coverage)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task.get("document_id") == DOC and task.get("evidence_id") == f"{DOC}:page7:text":
        task["state"] = "attempted"
write(RUN / "tasks.json", tasks)

print(json.dumps({
    "delta_molecules": len(bundle["molecules"]),
    "delta_reactions": len(bundle["reactions"]),
    "new_reactions": new_ids,
    "augmented_reaction": f"{DOC}:fig4_1b",
}, ensure_ascii=False))
