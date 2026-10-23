#!/usr/bin/env python3
"""Build the v3.1 delta bundle for inspected PDF pages 7--8.

Page 7 is references only.  Page 8 is a graphical abstract that repeats the
article's generic reaction scheme, so this pass enriches the existing generic
scheme instead of inventing a new concrete experiment.
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
DOC = "doc_439fe3cfd49f9dfc"
P3 = [f"{DOC}:page3:text", f"{DOC}:page3:image"]
P7 = [f"{DOC}:page7:text", f"{DOC}:page7:image"]
P8 = [f"{DOC}:page8:text", f"{DOC}:page8:image"]
SCHEME_EVIDENCE = P3 + P8


def fact(raw=None, normalized=None, evidence=(), *, method="original_extraction",
         derivation=None, status=None):
    return st.fact(
        raw,
        normalized,
        evidence,
        method=method,
        derivation=derivation,
        status=status,
    )


def add_evidence(value, evidence):
    """Append source evidence to all nested facts/occurrences without duplication."""
    if isinstance(value, dict):
        if "normalized_value" in value and "acquisition_method" in value:
            value["evidence_ids"] = list(dict.fromkeys(value["evidence_ids"] + evidence))
        elif "document_id" in value and "compound_label" in value and "evidence_ids" in value:
            value["evidence_ids"] = list(dict.fromkeys(value["evidence_ids"] + evidence))
        for item in value.values():
            add_evidence(item, evidence)
    elif isinstance(value, list):
        for item in value:
            add_evidence(item, evidence)


def variable(label, attachment, definition, constraints, evidence):
    return {
        "label": label,
        "attachment_points": fact(attachment, attachment, evidence),
        "definition": fact(definition, definition, evidence),
        "allowed_values": fact(evidence=evidence),
        "constraints": fact(constraints, constraints, evidence),
    }


def template(fixed_structure, variables, reason, evidence):
    return {
        "fixed_structure": fact(fixed_structure, fixed_structure, evidence),
        "variables": variables,
        "instance_mapping": [],
        "reason": fact(reason, reason, evidence),
    }


def set_stereo(molecule, chirality, geometry, representation, evidence, rationale):
    molecule["stereochemistry"]["chirality_status"] = fact(
        normalized=chirality,
        evidence=evidence,
        method="inference",
        derivation=rationale,
    )
    molecule["stereochemistry"]["geometry_status"] = fact(
        normalized=geometry,
        evidence=evidence,
        method="inference",
        derivation=rationale,
    )
    molecule["stereochemistry"]["representation_status"] = fact(
        normalized=representation,
        evidence=evidence,
        method="inference",
        derivation=rationale,
    )


def new_named_molecule(molecule_id, label, name, evidence):
    molecule = st.new_molecule(molecule_id, DOC, evidence)
    molecule["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    molecule["names"]["original"] = fact(name, name, evidence)
    molecule["names"]["english"] = fact(name, name, evidence)
    molecule["cas_numbers"] = fact(evidence=evidence)
    molecule["external_properties"] = []
    return molecule


template_bundle = json.loads((RUN / "raw_model/pages_0007_0008_template.json").read_text())
existing = json.loads((RUN / "standard_records.json").read_text())
molecules_by_id = {m["molecule_id"]: m for m in existing["molecules"]}
reactions_by_id = {r["reaction_id"]: r for r in existing["reactions"]}

bundle = st.empty(template_bundle["documents"])

# Re-express the three source-drawn R-group structures as templates.  No
# concrete SMILES is claimed, so the standard stereo command has no concrete
# graph to inspect in this page range.
substrate = copy.deepcopy(molecules_by_id["mol_generic_2_aminoacrylate"])
add_evidence(substrate, P8)
substrate["representation_type"] = "generic_template"
substrate["generic_template"] = template(
    "source-drawn 2-aminoacrylate framework bearing N-R1 and ester CO2R3",
    [
        variable("R1", "nitrogen substituent", "generic nitrogen protecting/substituent group", "identity not specified in the graphical abstract", SCHEME_EVIDENCE),
        variable("R3", "ester oxygen substituent", "generic ester substituent", "identity not specified in the graphical abstract", SCHEME_EVIDENCE),
    ],
    "The source uses R1 and R3 labels; expanding them to a single concrete molecule would invent identity.",
    SCHEME_EVIDENCE,
)
substrate["identity_note"] = fact(
    "R1/R3 generic 2-aminoacrylate in the source scheme.",
    "R1/R3 generic 2-aminoacrylate in the source scheme.",
    SCHEME_EVIDENCE,
)
set_stereo(
    substrate,
    "unresolved",
    "unresolved",
    "unresolved",
    SCHEME_EVIDENCE,
    "R1/R3 are not expanded and the graphical abstract does not assign stereochemistry; unspecified information is not treated as achiral or racemic.",
)

allenoate = copy.deepcopy(molecules_by_id["mol_generic_allenoate_2"])
add_evidence(allenoate, P8)
allenoate["representation_type"] = "generic_template"
allenoate["generic_template"] = template(
    "source-drawn terminal allenoate framework bearing ester CO2R2",
    [
        variable("R2", "ester oxygen substituent", "generic ester substituent", "identity not specified in the graphical abstract", SCHEME_EVIDENCE),
    ],
    "The source uses an R2 label; no concrete ester substituent is selected.",
    SCHEME_EVIDENCE,
)
allenoate["stereochemistry"]["axial"] = fact(
    "terminal allene framework shown; no axial descriptor assigned",
    "terminal allene framework shown; no axial descriptor assigned",
    SCHEME_EVIDENCE,
)
allenoate["identity_note"] = fact(
    "R2 generic allenoate in the source scheme; full identity remains unresolved.",
    "R2 generic allenoate in the source scheme; full identity remains unresolved.",
    SCHEME_EVIDENCE,
)
set_stereo(
    allenoate,
    "unresolved",
    "not_applicable",
    "unresolved",
    SCHEME_EVIDENCE,
    "The drawn terminal allene has no assigned axial configuration, but R2 remains generic and no concrete stereochemical graph is asserted.",
)

product = copy.deepcopy(molecules_by_id["mol_generic_3_pyrroline_3"])
add_evidence(product, P8)
product["representation_type"] = "generic_template"
product["generic_template"] = template(
    "source-drawn 3-pyrroline framework bearing R1, CO2R2 and CO2R3",
    [
        variable("R1", "ring nitrogen substituent", "generic nitrogen substituent", "identity not specified in the graphical abstract", SCHEME_EVIDENCE),
        variable("R2", "ester oxygen substituent", "generic ester substituent derived from the allenoate", "identity not specified in the graphical abstract", SCHEME_EVIDENCE),
        variable("R3", "ester oxygen substituent at the quaternary carbon", "generic ester substituent derived from the 2-aminoacrylate", "identity not specified in the graphical abstract", SCHEME_EVIDENCE),
    ],
    "The source product is a generic R-group drawing and cannot be represented as one concrete molecule.",
    SCHEME_EVIDENCE,
)
product["stereochemistry"]["atom_centers"] = fact(
    "An amino quaternary stereogenic center",
    "amino quaternary stereogenic center; absolute configuration unassigned",
    P8,
)
product["stereochemistry"]["assignment"] = fact(
    "An amino quaternary stereogenic center",
    "stereogenic center indicated; no absolute configuration or mixture identity reported",
    P8,
)
product["identity_note"] = fact(
    "Generic source-drawn 3-pyrroline; the highlighted amino quaternary stereogenic center has no absolute configuration assignment.",
    "Generic source-drawn 3-pyrroline; the highlighted amino quaternary stereogenic center has no absolute configuration assignment.",
    SCHEME_EVIDENCE,
)
set_stereo(
    product,
    "unspecified",
    "not_applicable",
    "not_encoded",
    SCHEME_EVIDENCE,
    "The source explicitly marks a stereogenic center but gives no absolute configuration; the generic template does not encode a stereoisomer.",
)

phosphine = new_named_molecule(
    "mol_generic_tertiary_phosphine_pr3",
    "cat. PR3",
    "generic tertiary phosphine catalyst PR3",
    SCHEME_EVIDENCE,
)
phosphine["representation_type"] = "generic_template"
phosphine["generic_template"] = template(
    "tertiary phosphine PR3",
    [
        variable("R", "three phosphorus substituent positions", "generic substituent on phosphorus", "catalyst substituent identity is not specified in the generic scheme", SCHEME_EVIDENCE),
    ],
    "The scheme states cat. PR3 rather than a specific phosphine; PPh3 is not substituted into the generic record.",
    SCHEME_EVIDENCE,
)
phosphine["identity_note"] = fact(
    "Generic catalyst notation cat. PR3; catalyst identity and loading are unreported in this scheme.",
    "Generic catalyst notation cat. PR3; catalyst identity and loading are unreported in this scheme.",
    SCHEME_EVIDENCE,
)
set_stereo(
    phosphine,
    "unresolved",
    "not_applicable",
    "unresolved",
    SCHEME_EVIDENCE,
    "The R substituents are not identified, so possible phosphorus stereochemistry and a concrete representation cannot be determined.",
)

solvent = new_named_molecule(
    "mol_generic_anhydrous_solvent",
    "anhydrous solvent",
    "unspecified anhydrous solvent",
    SCHEME_EVIDENCE,
)
solvent["identity_note"] = fact(
    "The generic scheme states only anhydrous solvent; composition and volume are unreported.",
    "The generic scheme states only anhydrous solvent; composition and volume are unreported.",
    SCHEME_EVIDENCE,
)
set_stereo(
    solvent,
    "unresolved",
    "unresolved",
    "unresolved",
    SCHEME_EVIDENCE,
    "Solvent identity is unreported, so molecular and stereochemical properties cannot be assigned.",
)

bundle["molecules"] = [substrate, allenoate, product, phosphine, solvent]

# Enrich the existing generic scheme.  Page 8 is duplicate evidence, not a new
# experiment, and its aggregate yield/example count is not made an outcome.
reaction_id = f"{DOC}:scheme1_this_work_generic"
reaction = copy.deepcopy(reactions_by_id[reaction_id])
add_evidence(reaction, P8)
reaction["location"] = fact(
    "Scheme 1 (This work) and page-8 graphical abstract",
    "Scheme 1 (This work) and page-8 graphical abstract",
    SCHEME_EVIDENCE,
)
reaction["participants"] = []
for participant_id, molecule_id, role in [
    ("substrate", substrate["molecule_id"], "reactant"),
    ("allenoate", allenoate["molecule_id"], "reactant"),
    ("catalyst", phosphine["molecule_id"], "catalyst"),
    ("solvent", solvent["molecule_id"], "solvent"),
    ("product", product["molecule_id"], "product"),
]:
    participant = st.new_participant(participant_id, molecule_id, role, SCHEME_EVIDENCE)
    if role == "product":
        participant["product_role"] = fact("generic main product", "generic main product", SCHEME_EVIDENCE)
    reaction["participants"].append(participant)

reaction["solvent_conditions"]["composition"] = fact(
    "anhydrous solvent",
    [solvent["molecule_id"]],
    SCHEME_EVIDENCE,
)
step = st.new_step("generic_reaction", SCHEME_EVIDENCE)
step["participant_ids"] = ["substrate", "allenoate", "catalyst", "solvent", "product"]
step["operation"] = fact(
    "phosphine-catalyzed formal [3 + 2] annulation",
    "phosphine-catalyzed formal [3 + 2] annulation",
    SCHEME_EVIDENCE,
)
step["temperature"] = fact("r.t.", "r.t.", SCHEME_EVIDENCE)
step["time"] = fact(evidence=SCHEME_EVIDENCE)
step["atmosphere"] = fact("Ar", "Ar", SCHEME_EVIDENCE)
reaction["steps"] = [step]
reaction["procedure"]["original_text"] = fact(
    "cat. PR3; anhydrous solvent; Ar, r.t.",
    "cat. PR3; anhydrous solvent; Ar, r.t.",
    SCHEME_EVIDENCE,
)
reaction["outcomes"] = []
reaction["notes"] = [
    fact("In situ generated ketoimine", "In situ generated ketoimine", SCHEME_EVIDENCE),
    fact("An amino quaternary stereogenic center", "An amino quaternary stereogenic center", P8),
    fact(
        "up to 94% yield; 30 examples",
        "article-level graphical-abstract aggregate; not a per-experiment outcome",
        P8,
    ),
]
reaction["evidence_ids"] = SCHEME_EVIDENCE
reaction["review"] = {
    "status": "needs_review",
    "method": "model",
    "reviewer": None,
    "checked_fields": [],
    "unresolved": [
        "Generic R1/R2/R3 groups are not expanded to concrete molecules.",
        "The generic PR3 catalyst identity/loading, solvent identity/volume and reaction time are unreported.",
        "The graphical-abstract statements 'up to 94% yield' and '30 examples' are aggregate claims and are not assigned to one experiment or product outcome.",
        "The highlighted product stereogenic center has no absolute configuration or mixture assignment.",
    ],
}
reaction["completeness"] = "partial"
bundle["reactions"] = [reaction]

bundle["source_context"] = [
    {
        "context_id": f"{DOC}:page7_references",
        "document_id": DOC,
        "evidence_ids": P7,
        "kind": "references",
        "summary": "Page 7 contains references 6-18 only; cited literature is not treated as experiments performed in this article.",
    },
    {
        "context_id": f"{DOC}:page8_graphical_abstract",
        "document_id": DOC,
        "evidence_ids": P8,
        "kind": "graphical_abstract",
        "summary": "Graphical abstract repeats the generic this-work annulation and gives aggregate claims 'up to 94% yield' and '30 examples'; these are not a new concrete experiment or row-specific outcome.",
        "linked_reaction_id": reaction_id,
    },
]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "two-column references section (references 6-18)",
                    "finding": "Reference list only; cited papers are prior literature and no experiment performed in this article is shown on the page.",
                    "evidence_ids": P7,
                }
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "left graphical reaction panel",
                    "finding": "Generic phosphine-catalyzed [3+2] annulation under cat. PR3, anhydrous solvent, Ar and room temperature; duplicate evidence for the existing this-work generic scheme, not a new concrete experiment.",
                    "evidence_ids": P8,
                },
                {
                    "region": "pink feature box and yield/example callouts",
                    "finding": "Article-level summary claims an amino quaternary stereogenic center, up to 94% yield and 30 examples; aggregate claims are not assigned as a product-specific outcome.",
                    "evidence_ids": P8,
                },
                {
                    "region": "right summary paragraph and author line",
                    "finding": "Graphical-abstract narrative and bibliographic content; no additional table row, reaction arrow or experimental procedure.",
                    "evidence_ids": P8,
                },
            ],
            "record_ids": [reaction_id],
            "unresolved": [
                "R1/R2/R3 identities, PR3 identity/loading, anhydrous-solvent identity/volume and reaction time are not specified in the graphical abstract.",
                "No absolute configuration or mixture identity is assigned to the highlighted stereogenic center.",
                "The 94% maximum is not tied to a single row on this page and its measurement basis is not reported here.",
            ],
        },
    ],
}

raw_extraction = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "inspected_original_pdf_pages": [7, 8],
    "image_evidence": [P7[1], P8[1]],
    "observations": {
        "page7": "References 6-18 only; no current-paper reaction, table, arrow, footnote condition or characterization block.",
        "page8": {
            "content_type": "graphical abstract",
            "reaction_text": "generic allenoate + generic 2-aminoacrylate; cat. PR3; anhydrous solvent; Ar, r.t.; generic 3-pyrroline",
            "stereochemistry_text": "An amino quaternary stereogenic center",
            "aggregate_text": ["up to 94% yield", "30 examples", "Gram-scale synthesis", "Mechanistic studies"],
        },
    },
    "decisions": [
        "Merge page-8 evidence into the existing generic this-work scheme rather than create a duplicate concrete experiment.",
        "Represent source R groups and PR3 as generic templates; leave concrete SMILES empty.",
        "Use generic PR3 rather than substituting PPh3 into the source's cat. PR3 notation.",
        "Correct the generic this-work atmosphere to Ar; N2 belongs to the adjacent previous-work schemes on page 3.",
        "Do not turn 'up to 94% yield; 30 examples' into an outcome because it is an article-level aggregate without a unique target experiment.",
        "Keep the marked product center unspecified; do not infer an absolute configuration or racemate.",
        "No SI association is inferred.",
    ],
    "new_concrete_experiments": 0,
    "new_reaction_records": 0,
    "updated_record_ids": [reaction_id],
    "new_molecule_ids": [phosphine["molecule_id"], solvent["molecule_id"]],
    "standard_stereo": "No new concrete SMILES occurs on pages 7-8; generic templates cannot be passed to concrete graph stereo inspection.",
    "semantic_accuracy_verified": False,
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "coverage/pages_0007_0008.json", coverage)
st.write(RUN / "raw_model/pages_0007_0008_extraction.json", raw_extraction)

print(json.dumps({
    "molecules_in_delta": len(bundle["molecules"]),
    "reactions_in_delta": len(bundle["reactions"]),
    "source_context_in_delta": len(bundle["source_context"]),
    "coverage_pages": [p["page"] for p in coverage["pages"]],
}, ensure_ascii=False))
