#!/usr/bin/env python3
"""Build the page-7-only v3.1 submission extension for DOI 10.1039/d1sc02681j.

The page contains references and three substantive footnotes.  Only footnote 19
is sufficiently distinct to support a new experiment record.  Footnotes 20 and
21 are retained as source context without inventing an SI association or a
separately identifiable reaction row.
"""

from copy import deepcopy
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


DOC = "doc_68d6021490685b69"
PAGE7_TEXT = f"{DOC}:page7:text"
PAGE7_IMAGE = f"{DOC}:page7:image"
E7 = [PAGE7_TEXT, PAGE7_IMAGE]
PAGE2_TEXT = f"{DOC}:page2:text"


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None, verification=None, alternatives=()):
    return st.fact(
        raw,
        normalized,
        evidence,
        method=method,
        derivation=derivation,
        status=status,
        verification=verification,
        alternatives=alternatives,
    )


def qty(number, unit, raw_text):
    return {"number": number, "unit": unit, "raw_text": raw_text}


def inherited(raw, normalized):
    return F(
        raw,
        normalized,
        [PAGE2_TEXT],
        method="inherited_general_condition",
        derivation=(
            "Footnote 19 explicitly identifies the experiment as Table 1, "
            "entry 1 carried out at 0 °C; the remaining Table 1 entry-1 "
            "conditions are inherited without treating unavailable ESI Table "
            "S1 entry 15 as supplied evidence."
        ),
    )


template = json.loads((RUN / "raw_model/page7_init_template.json").read_text())
old = json.loads((RUN / "standard_records.json").read_text())
old_molecules = {m["molecule_id"]: m for m in old["molecules"]}

# Reuse exact registered identities for the model substrates, catalyst and
# solvent.  They are copied unchanged into the incoming bundle only so the
# standalone submission validates before merge.
for molecule_id in ("m_1a", "m_2a", "a_NiOTf2", "a_CH2Cl2"):
    template["molecules"].append(deepcopy(old_molecules[molecule_id]))


# Page 7 does not print a compound number or structure for the detected [2+2]
# product.  A separate identity prevents conflating a product from "the racemic
# reaction" with the enantioenriched, X-ray-assigned m_3aa registry entry.
product_id = "m_note19_detected_2plus2_product"
product = st.new_molecule(product_id, DOC, E7)
product["identity_status"] = "source_identified"
product["occurrences"] = [{
    "document_id": DOC,
    "compound_label": F(
        "the [2+2] cycloaddition product",
        "[2+2] cycloaddition product detected in footnote 19",
        E7,
    ),
    "evidence_ids": E7,
}]
product["names"]["original"] = F(
    "the [2+2] cycloaddition product",
    "the [2+2] cycloaddition product",
    E7,
)
product["names"]["english"] = F(
    "the [2+2] cycloaddition product",
    "[2+2] cycloaddition product detected in the 0 °C racemic control",
    E7,
)
product["names"]["abbreviations"] = F()
product["names"]["aliases"] = F()
product["cas_numbers"] = F()
product["representation_type"] = "unresolved"
product["smiles"] = F()
product["generic_template"] = None
product["stereochemistry"]["atom_centers"] = F()
product["stereochemistry"]["double_bonds"] = F()
product["stereochemistry"]["axial"] = F()
product["stereochemistry"]["mixture"] = F()
product["stereochemistry"]["assignment"] = F(
    "the racemic reaction",
    {
        "reaction_context": "described by the authors as racemic",
        "product_stereochemical_assignment": "not explicitly reported",
    },
    E7,
    method="inference",
    derivation=(
        "The adjective 'racemic' modifies the reaction in footnote 19, not an "
        "explicit product ee, enantiomer ratio or absolute configuration; no "
        "product stereochemical composition is asserted."
    ),
)
product["stereochemistry"]["chirality_status"] = F(
    "product configuration not printed",
    "unresolved",
    E7,
    method="inference",
    derivation=(
        "The page reports a racemic reaction but neither labels the product "
        "racemic nor gives an enantiomeric composition."
    ),
)
product["stereochemistry"]["geometry_status"] = F(
    "product structure not printed",
    "unresolved",
    E7,
    method="inference",
    derivation="No product drawing or E/Z/cis/trans assignment appears in footnote 19.",
)
product["stereochemistry"]["representation_status"] = F(
    "product structure not printed",
    "not_encoded",
    E7,
    method="inference",
    derivation=(
        "No concrete structure is printed on page 7, so no SMILES or generic "
        "structure is introduced from context alone."
    ),
)
product["stereochemistry"]["isomeric_smiles"] = F()
product["stereochemistry"]["identity_key"] = F()
product["identity_note"] = F(
    "Footnote 19 does not print a compound number or structure.",
    (
        "Kept separate from m_3aa because the page-7 note does not explicitly "
        "identify the detected product as the single X-ray-assigned enantiomer."
    ),
    E7,
)
template["molecules"].append(product)


def participant(participant_id, molecule_id, role, *, charges=None,
                reference=None, product_role=None):
    p = st.new_participant(participant_id, molecule_id, role, E7)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, E7)
    for key, fact in (charges or {}).items():
        p["charges"][key] = fact
    if reference is not None:
        p["equivalence_reference"] = inherited(reference, reference)
    return p


reaction_id = f"{DOC}:footnote19_racemic_0C_2plus2_detection"
reaction = st.new_reaction(reaction_id, DOC, E7)
reaction["record_type"] = "specific_experiment"
reaction["experiment_categories"] = F(
    "racemic reaction; mechanistic inference",
    ["control", "mechanistic_study"],
    E7,
)
reaction["location"] = F(
    "footnote 19",
    {"page": 7, "location": "right column, footnote 19"},
    E7,
)
reaction["general_scheme_ids"] = []
reaction["reaction_name"] = F(
    "the racemic reaction (Table 1, entry 1) was carried out at 0 °C",
    "0 °C racemic [2+2] cycloaddition control",
    E7,
)
reaction["transformation_type"] = F(
    "[2+2] cycloaddition product was detected",
    "[2+2] cycloaddition",
    E7,
    method="inference",
    derivation="Normalized directly from the product class named in footnote 19.",
)

reaction["participants"] = [
    participant(
        "p_1a",
        "m_1a",
        "reactant",
        charges={"amount": inherited(
            "Table 1 entry 1: 1a (0.05 mmol)",
            qty(0.05, "mmol", "0.05 mmol"),
        )},
    ),
    participant(
        "p_2a",
        "m_2a",
        "reactant",
        charges={"equivalents": inherited(
            "Table 1 entry 1: 2a (1.0 equiv.)",
            qty(1.0, "equiv", "1.0 equiv."),
        )},
        reference="m_1a",
    ),
    participant(
        "p_ni",
        "a_NiOTf2",
        "catalyst",
        charges={"mol_percent": inherited(
            "Table 1 entry 1: Ni(OTf)2 (10 mol%)",
            qty(10, "mol%", "10 mol%"),
        )},
        reference="m_1a",
    ),
    participant("p_solvent", "a_CH2Cl2", "solvent"),
    participant(
        "p_detected_product",
        product_id,
        "product",
        product_role="detected [2+2] cycloaddition product",
    ),
]
reaction["reaction_scale"] = inherited(
    "Table 1 entry 1: 1a (0.05 mmol)",
    {
        "reference_molecule_id": "m_1a",
        "amount": qty(0.05, "mmol", "0.05 mmol"),
    },
)
reaction["solvent_conditions"]["composition"] = inherited(
    "Table 1 entry 1: CH2Cl2",
    "CH2Cl2",
)
reaction["solvent_conditions"]["concentration"] = inherited(
    "Table 1 entry 1: CH2Cl2 (0.1 M)",
    qty(0.1, "M", "0.1 M"),
)
reaction["solvent_conditions"]["concentration_basis"] = inherited(
    "Table 1 entry 1: concentration based on 1a",
    "m_1a",
)

step = st.new_step("racemic_low_temperature_control", E7)
step["participant_ids"] = [
    "p_1a", "p_2a", "p_ni", "p_solvent", "p_detected_product"
]
step["operation"] = F(
    "the racemic reaction (Table 1, entry 1) was carried out at 0 °C",
    "carry out Table 1 entry 1 under racemic conditions at 0 °C",
    E7,
)
step["temperature"] = F(
    "0 °C",
    qty(0, "°C", "0 °C"),
    E7,
)
step["time"] = inherited(
    "Table 1 entry 1: 16 h",
    qty(16, "h", "16 h"),
)
step["atmosphere"] = inherited(
    "Table 1 entry 1: under N2",
    "N2",
)
reaction["steps"] = [step]

detected = st.new_outcome("detected_2plus2_product", "other", [product_id], E7)
detected["product_role"] = F(
    "the [2+2] cycloaddition product",
    "detected product",
    E7,
)
detected["value"] = F(
    "was detected",
    "detected",
    E7,
)
detected["yield_type"] = F()
detected["measurement_basis"] = F()
detected["scope"] = F(
    "the [2+2] cycloaddition product",
    (
        "qualitative detection only; compound number, amount, yield, detection "
        "method and stereochemical composition are unreported"
    ),
    E7,
)
reaction["outcomes"] = [detected]

changes = []
for change_id in ("b1", "b2"):
    changes.append(F(
        "[2+2] cycloaddition: new C-C connection (partial)",
        {
            "change_id": change_id,
            "atoms": [
                {"atom_map": None, "element": "C"},
                {"atom_map": None, "element": "C"},
            ],
            "before_order": 0,
            "after_order": 1,
            "change_type": "formed",
        },
        E7,
        method="structure_derivation",
        derivation=(
            "The authors explicitly call the detected material a [2+2] "
            "cycloaddition product.  Two C-C connections are the partial net "
            "classification; no atom correspondence or concrete product "
            "structure is asserted from page 7."
        ),
    ))
reaction["classification"]["bond_changes"] = changes
reaction["classification"]["labels"] = F(
    "[2+2] cycloaddition: two C-C connections (partial)",
    [{
        "category_number": 2,
        "code": "C-C",
        "supporting_change_ids": ["b1", "b2"],
    }],
    E7,
    method="structure_derivation",
    derivation=(
        "Partial label derived from the source reaction class; all atom maps, "
        "bond-order changes and stereochemical correspondence remain unresolved."
    ),
)
reaction["classification"]["coverage"] = F(
    "partial",
    "partial",
    E7,
    method="inference",
    derivation=(
        "Only the two characteristic new C-C connections are classified; the "
        "page does not print the product structure or a mapped reaction."
    ),
)
reaction["classification"]["limitations"] = [
    "No product compound number, concrete structure or atom mapping is printed in footnote 19.",
    "The author describes the reaction as racemic, but no product ee/er or absolute configuration is reported.",
    "ESI Table S1 entry 15 is cited but no SI document is associated in this workspace.",
]
reaction["classification"]["model_votes"] = []
reaction["notes"] = [
    F(
        "Therefore, we deemed that the formal [4+2] adduct was generated from the transformation of the [2+2] cycloaddition product.",
        "Author's mechanistic interpretation; not treated as direct proof of the pathway.",
        E7,
    ),
    F(
        "See ESI, Table S1 entry 15 for details.",
        "Cited ESI details are unavailable and are not inferred or associated.",
        E7,
    ),
]
reaction["evidence_ids"] = E7
reaction["review"] = {
    "status": "needs_review",
    "method": "model",
    "reviewer": None,
    "checked_fields": [],
    "unresolved": [
        "The detected product has no printed compound number, structure, yield/amount, detection method, ee/er or configuration in footnote 19.",
        "The inherited Table 1 entry-1 conditions are explicit, but ESI Table S1 entry 15 is unavailable for any additional modifications or outcomes.",
    ],
}
reaction["completeness"] = "partial"
template["reactions"].append(reaction)


template["source_context"] = [
    {
        "context_id": f"{DOC}:page7_references",
        "document_id": DOC,
        "kind": "bibliographic_references",
        "evidence_ids": E7,
        "text": (
            "References 12–18 and 22–23 are cited prior literature/reviews; "
            "they do not constitute current-paper experiments on page 7."
        ),
    },
    {
        "context_id": f"{DOC}:page7_footnote19_interpretation",
        "document_id": DOC,
        "kind": "mechanistic_interpretation",
        "evidence_ids": E7,
        "related_reaction_ids": [reaction_id],
        "text": (
            "The authors use detection of the low-temperature [2+2] product to "
            "infer that the formal [4+2] adduct arises by transformation of the "
            "[2+2] product; this is an interpretation, not a separate experiment."
        ),
    },
    {
        "context_id": f"{DOC}:page7_footnote20_crystallography",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": E7,
        "text": (
            "Footnote 20 lists CCDC 2042864 (3aa), 2052253 (rac-4ka) and "
            "2079315 (rac-5aa); crystallographic identifiers are characterization "
            "context rather than new synthesis experiments."
        ),
    },
    {
        "context_id": f"{DOC}:page7_footnote21_onepot",
        "document_id": DOC,
        "kind": "qualitative_experimental_note",
        "evidence_ids": E7,
        "text": (
            "Footnote 21 states that performing 'the reaction' as a one-pot "
            "procedure resulted in product family 4 as the major product and "
            "points to ESI Table S5.  The page does not identify a substrate "
            "instance, conditions, ratio or yield, so it is not promoted to a "
            "distinct reaction record without the unassociated ESI."
        ),
    },
]


(RUN / "submitted.json").write_text(
    json.dumps(template, ensure_ascii=False, indent=2) + "\n"
)

raw = {
    "document_id": DOC,
    "page": 7,
    "page_visually_inspected": True,
    "evidence_ids": E7,
    "regions": [
        {
            "region": "references 12–18",
            "raw_finding": "bibliographic citations only; no table, scheme or current-paper reaction arrow",
        },
        {
            "region": "footnote 19",
            "raw_finding": (
                "Table 1 entry 1 was run at 0 °C under racemic conditions; an "
                "unnumbered [2+2] cycloaddition product was detected"
            ),
        },
        {
            "region": "footnote 20",
            "raw_finding": "CCDC identifiers for 3aa, rac-4ka and rac-5aa",
        },
        {
            "region": "footnote 21",
            "raw_finding": (
                "one-pot procedure gave product family 4 as major product; no "
                "instance, yield, ratio or detailed condition on this page"
            ),
        },
        {
            "region": "references 22–23 and footer",
            "raw_finding": "bibliographic citations and publication footer only",
        },
    ],
    "new_reaction_ids": [reaction_id],
    "not_promoted_to_reactions": [
        "footnote 20 crystallographic identifiers",
        "footnote 21 under-specified qualitative one-pot note",
        "bibliographic references",
    ],
    "semantic_accuracy_verified": False,
}
(RUN / "raw_model/page7_visual_extraction.json").write_text(
    json.dumps(raw, ensure_ascii=False, indent=2) + "\n"
)

# Page 7 has no structure drawing of its own.  These diagnostics document that
# the concrete structures reused from the explicit Table 1 reference retain
# their previously registered stereochemical graph; they do not certify the
# unnumbered detected product.
reused_stereo = {
    "m_1a": st.stereo_inspect(old_molecules["m_1a"]["smiles"]["normalized_value"]),
    "m_2a": st.stereo_inspect(old_molecules["m_2a"]["smiles"]["normalized_value"]),
    "page7_detected_product": {
        "standard_stereo_run": "not_applicable",
        "reason": "No concrete product structure or SMILES is printed in footnote 19.",
    },
    "chemical_accuracy_verified": False,
}
(RUN / "raw_model/page7_reused_stereo_diagnostics.json").write_text(
    json.dumps(reused_stereo, ensure_ascii=False, indent=2) + "\n"
)

print(json.dumps({
    "incoming_molecules": len(template["molecules"]),
    "new_molecules": 1,
    "new_reactions": 1,
    "reaction_id": reaction_id,
}, ensure_ascii=False))
