#!/usr/bin/env python3
"""Merge source facts visually extracted from PDF pages 7-8 into the v3.1 bundle."""
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st


RUN = ROOT / "runs/current"
DOC = "doc_a52c76191d5b058b"
PAGE7 = [f"{DOC}:page7:image", f"{DOC}:page7:text"]
PAGE8 = [f"{DOC}:page8:image", f"{DOC}:page8:text"]
PAGE5_7 = [f"{DOC}:page5:image", f"{DOC}:page5:text", *PAGE7]


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None):
    return st.fact(raw, normalized, evidence, method=method,
                   derivation=derivation, status=status)


def generic_template_molecule(mid, label, name, fixed_structure, variables,
                              stereo):
    m = st.new_molecule(mid, DOC, PAGE7)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, PAGE7)
    m["names"]["original"] = F(label, label, PAGE7)
    m["names"]["english"] = F(name, name, PAGE7)
    m["representation_type"] = "generic_template"
    m["smiles"] = F(
        "Scheme 3 generic drawing",
        None,
        PAGE7,
        "original_extraction",
        "A unique concrete SMILES is not assigned because Scheme 3 contains R and/or L_n variables.",
    )
    template_vars = []
    for variable in variables:
        template_vars.append({
            "label": variable["label"],
            "attachment_points": F(
                variable["attachment_raw"], variable["attachment"], PAGE7),
            "definition": F(variable["definition_raw"], variable["definition"], PAGE7),
            "allowed_values": F(
                variable.get("allowed_raw"), variable.get("allowed"), PAGE7,
                "original_extraction", status="unreported"
                if variable.get("allowed") is None else None),
            "constraints": F(
                variable.get("constraints_raw"), variable.get("constraints"), PAGE7,
                "original_extraction", status="unreported"
                if variable.get("constraints") is None else None),
        })
    m["generic_template"] = {
        "fixed_structure": F("Scheme 3 drawing", fixed_structure, PAGE7),
        "variables": template_vars,
        "instance_mapping": [],
        "reason": F(
            "R and/or L_n shown in Scheme 3",
            "The source is a proposed generic catalytic cycle, so variable groups and ligand sets are retained instead of expanded.",
            PAGE7,
        ),
    }
    common = (
        "Stereochemistry was reviewed from the Scheme 3 drawing. Generic R/L_n variables and "
        "metal coordination prevent a unique fully encoded molecular graph."
    )
    m["stereochemistry"]["atom_centers"] = F(
        stereo["atom_raw"], stereo["atom"], PAGE7, "inference", common)
    m["stereochemistry"]["double_bonds"] = F(
        stereo["double_raw"], stereo["double"], PAGE7, "inference", common)
    m["stereochemistry"]["axial"] = F(
        stereo["axial_raw"], stereo["axial"], PAGE7, "inference", common)
    m["stereochemistry"]["mixture"] = F(
        "No mixture composition is assigned in Scheme 3",
        "unresolved for the generic mechanistic entity",
        PAGE7,
        "inference",
        common,
    )
    m["stereochemistry"]["assignment"] = F(
        stereo["assignment_raw"], stereo["assignment"], PAGE7,
        "inference", common)
    m["stereochemistry"]["chirality_status"] = F(
        None, stereo["chirality_status"], PAGE7, "inference", common)
    m["stereochemistry"]["geometry_status"] = F(
        None, stereo["geometry_status"], PAGE7, "inference", common)
    m["stereochemistry"]["representation_status"] = F(
        None, "partial", PAGE7, "inference", common)
    m["stereochemistry"]["isomeric_smiles"] = F(
        None, None, PAGE7, "original_extraction",
        "No unique isomeric SMILES can represent the generic mechanistic entity.",
        "unreported",
    )
    m["stereochemistry"]["identity_key"] = F(
        None, None, PAGE7, "original_extraction",
        "No stereo graph key was generated without a unique molecular graph.",
        "unreported",
    )
    m["identity_note"] = F(
        f"Scheme 3 label {label}",
        "Proposed mechanistic entity, not an isolated compound or independently observed experiment on pages 7-8.",
        PAGE7,
    )
    return m


def add_occurrence(molecule, label, evidence):
    if not any(set(o["evidence_ids"]) == set(evidence)
               for o in molecule["occurrences"]):
        molecule["occurrences"].append({
            "document_id": DOC,
            "compound_label": F(label, label, evidence),
            "evidence_ids": list(evidence),
        })


def participant(pid, mid, role, product_role=None):
    p = st.new_participant(pid, mid, role, PAGE7)
    if product_role:
        p["product_role"] = F(product_role, product_role, PAGE7)
    return p


def proposed_step(step_id, participant_ids, raw, normalized):
    step = st.new_step(step_id, PAGE5_7)
    step["participant_ids"] = participant_ids
    step["operation"] = F(raw, normalized, PAGE5_7)
    return step


bundle = st.read(RUN / "submitted.json")

# Evidence declarations must match the dataset-prepared document exactly. The real
# title remains in source_context because standard init preserves document.title=null.
template = st.read(RUN / "raw_model/pages_0007_0008_template.json")
bundle["documents"] = template["documents"]

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}

# Page 7 supplies a new source-explicit absolute assignment for 3al.
m3al = molecules["m_3al"]
add_occurrence(m3al, "(2S,3S)-3al", PAGE7)
m3al["stereochemistry"]["atom_centers"] = F(
    "(2S,3S)-3al",
    [
        {"source_position": 2, "assignment": "S"},
        {"source_position": 3, "assignment": "S"},
    ],
    PAGE7,
)
m3al["stereochemistry"]["double_bonds"] = F(
    "terminal vinyl group in the 3al drawing",
    [{"source_locator": "terminal vinyl group", "assignment": "E/Z not applicable"}],
    [f"{DOC}:page5:image", *PAGE7],
    "structure_derivation",
    "The displayed alkene is terminal; no E/Z pair is defined for that double bond.",
)
m3al["stereochemistry"]["axial"] = F(
    None, "not_applicable", [f"{DOC}:page5:image", *PAGE7],
    "inference",
    "The product drawing is not an allene or other explicitly axially chiral framework.",
)
m3al["stereochemistry"]["assignment"] = F(
    "X-ray data for compounds (2S,3S)-3al (CCDC 1999447)",
    {
        "absolute_configuration": "(2S,3S)",
        "assignment_context": "X-ray data listed in Associated Content",
        "ccdc": "1999447",
    },
    PAGE7,
)
m3al["stereochemistry"]["chirality_status"] = F(
    "(2S,3S)-3al", "specified", PAGE7,
)
m3al["stereochemistry"]["geometry_status"] = F(
    None, "not_applicable", [f"{DOC}:page5:image", *PAGE7],
    "structure_derivation",
    "The displayed product contains a terminal vinyl group for which E/Z is not applicable.",
)
m3al["stereochemistry"]["representation_status"] = F(
    "(2S,3S)-3al; full X-ray/SI structure not supplied in this workspace",
    "not_encoded",
    PAGE7,
    "inference",
    "The absolute descriptors are source-explicit, but a full candidate structure was not encoded without the associated SI/X-ray file.",
)
m3al["stereochemistry"]["isomeric_smiles"] = F(
    None, None, PAGE7, "original_extraction",
    "The source-explicit (2S,3S) assignment is retained without guessing a full isomeric SMILES.",
    "unreported",
)
m3al["identity_note"] = F(
    "X-ray data for compounds (2S,3S)-3al (CCDC 1999447)",
    "Page 7 explicitly identifies 3al as (2S,3S) and lists CCDC 1999447; the X-ray/SI file itself is not associated in the workspace.",
    PAGE7,
)

# 1a is drawn again in Scheme 3; append occurrence evidence without creating a duplicate identity.
add_occurrence(molecules["m_1a"], "1a", PAGE7)

var_r = [{
    "label": "R",
    "attachment_raw": "RO shown in Scheme 3",
    "attachment": "substituent attached through oxygen",
    "definition_raw": "R",
    "definition": "unspecified alkyl/organic substituent in the proposed generic cycle",
    "allowed_raw": None,
    "allowed": None,
    "constraints_raw": None,
    "constraints": None,
}]
var_r_ln = var_r + [{
    "label": "L_n",
    "attachment_raw": "L_nPd shown in Scheme 3",
    "attachment": "ligand set bound to palladium",
    "definition_raw": "L_n",
    "definition": "unspecified palladium ligand set in the proposed cycle",
    "allowed_raw": None,
    "allowed": None,
    "constraints_raw": None,
    "constraints": None,
}]

generic_specs = [
    (
        "m_scheme3_allene_2", "2", "generic alkoxyallene 2",
        "RO-substituted terminal allene drawn as the allene input to intermediate A",
        var_r,
        {
            "atom_raw": "No tetrahedral stereocenter specified in the drawn fixed allene fragment",
            "atom": [],
            "double_raw": "terminal cumulated C=C=C unit",
            "double": [{"type": "cumulated allene", "E_Z": "not applicable"}],
            "axial_raw": "terminal allene with a CH2 terminus",
            "axial": "the allene axis is not stereogenic as drawn because one terminus is CH2; R remains unspecified",
            "assignment_raw": "No stereochemical descriptor printed",
            "assignment": "R-group stereochemistry unresolved",
            "chirality_status": "unresolved",
            "geometry_status": "not_applicable",
        },
    ),
    (
        "m_scheme3_intermediate_A", "A", "proposed palladium hydride intermediate A",
        "amino-acid carboxylate-bound L_nPd(II)-H intermediate A",
        var_r_ln[1:],
        {
            "atom_raw": "amino-acid-derived stereocenter drawn without an R/S label",
            "atom": [{"source_locator": "carbon bearing NHPh and Bn", "assignment": "unspecified"}],
            "double_raw": "carbonyl shown; no alkene E/Z assignment",
            "double": [{"type": "carbonyl", "E_Z": "not applicable"}],
            "axial_raw": "No axial descriptor shown",
            "axial": "unresolved/not indicated for the generic metal complex",
            "assignment_raw": "No R/S, R_a/S_a, P/M, or cis/trans descriptor printed",
            "assignment": "unresolved",
            "chirality_status": "unspecified",
            "geometry_status": "not_applicable",
        },
    ),
    (
        "m_scheme3_intermediate_B", "B", "proposed pi-allyl palladium intermediate B",
        "amino-acid carboxylate-bound L_nPd(II) pi-allyl species B with an RO substituent",
        var_r_ln,
        {
            "atom_raw": "amino-acid-derived center and allyl complex shown without stereochemical descriptors",
            "atom": [{"source_locator": "carbon bearing NHPh and Bn", "assignment": "unspecified"}],
            "double_raw": "pi-allyl depiction; no E/Z descriptor printed",
            "double": "coordination drawing not converted to a unique alkene geometry",
            "axial_raw": "No axial descriptor shown",
            "axial": "unresolved/not indicated",
            "assignment_raw": "No R/S, R_a/S_a, P/M, or cis/trans descriptor printed",
            "assignment": "unresolved",
            "chirality_status": "unspecified",
            "geometry_status": "unresolved",
        },
    ),
    (
        "m_scheme3_intermediate_C", "C", "proposed palladium(III) intermediate C",
        "amino-acid carboxylate-bound L_nPd(III) allyl species C with an RO substituent",
        var_r_ln,
        {
            "atom_raw": "amino-acid-derived center and allyl complex shown without stereochemical descriptors",
            "atom": [{"source_locator": "carbon bearing NHPh and Bn", "assignment": "unspecified"}],
            "double_raw": "allyl depiction; no E/Z descriptor printed",
            "double": "coordination drawing not converted to a unique alkene geometry",
            "axial_raw": "No axial descriptor shown",
            "axial": "unresolved/not indicated",
            "assignment_raw": "No R/S, R_a/S_a, P/M, or cis/trans descriptor printed",
            "assignment": "unresolved",
            "chirality_status": "unspecified",
            "geometry_status": "unresolved",
        },
    ),
    (
        "m_scheme3_intermediate_D", "D", "proposed palladium(III) intermediate D",
        "C-C-coupled L_nPd(III) allyl species D with Bn, NHPh, and RO substituents",
        var_r_ln,
        {
            "atom_raw": "carbon bearing Bn and NHPh is drawn without an R/S label",
            "atom": [{"source_locator": "carbon bearing Bn and NHPh", "assignment": "unspecified"}],
            "double_raw": "allyl depiction; no E/Z descriptor printed",
            "double": "coordination drawing not converted to a unique alkene geometry",
            "axial_raw": "No axial descriptor shown",
            "axial": "unresolved/not indicated",
            "assignment_raw": "No R/S, R_a/S_a, P/M, or cis/trans descriptor printed",
            "assignment": "unresolved",
            "chirality_status": "unspecified",
            "geometry_status": "unresolved",
        },
    ),
    (
        "m_scheme3_ester_5aa", "5aa", "generic proposed ester intermediate 5aa",
        "amino-acid carboxylate ester of an RO-substituted allylic fragment",
        var_r,
        {
            "atom_raw": "amino-acid-derived center shown without an R/S label",
            "atom": [{"source_locator": "carbon bearing NHPh and Bn", "assignment": "unspecified"}],
            "double_raw": "terminal vinyl group",
            "double": [{"type": "terminal alkene", "E_Z": "not applicable"}],
            "axial_raw": "No axial framework or descriptor shown",
            "axial": "not applicable to the fixed drawn fragment",
            "assignment_raw": "No R/S descriptor printed",
            "assignment": "unresolved",
            "chirality_status": "unspecified",
            "geometry_status": "not_applicable",
        },
    ),
    (
        "m_scheme3_product_3", "3", "generic branched hydroaminoalkylation product 3",
        "branched terminal-alkene product bearing Bn, NHPh, and OR substituents",
        var_r,
        {
            "atom_raw": "solid wedge shown at the carbon bearing NHPh; no absolute descriptor printed",
            "atom": [{"source_locator": "carbon bearing NHPh", "assignment": "stereochemistry depicted; absolute descriptor unspecified"}],
            "double_raw": "terminal vinyl group",
            "double": [{"type": "terminal alkene", "E_Z": "not applicable"}],
            "axial_raw": "No axial framework or descriptor shown",
            "axial": "not applicable to the fixed drawn fragment",
            "assignment_raw": "wedge bond in generic product 3; no R/S descriptor",
            "assignment": "partially specified drawing; absolute configuration unresolved",
            "chirality_status": "partially_specified",
            "geometry_status": "not_applicable",
        },
    ),
]

for spec in generic_specs:
    mid = spec[0]
    molecule = generic_template_molecule(*spec)
    if mid in molecules:
        index = next(i for i, existing in enumerate(bundle["molecules"])
                     if existing["molecule_id"] == mid)
        bundle["molecules"][index] = molecule
    else:
        bundle["molecules"].append(molecule)
    molecules[mid] = molecule

reaction_id = f"{DOC}:proposed_mechanism_text"
reaction = next(r for r in bundle["reactions"] if r["reaction_id"] == reaction_id)
reaction["experiment_categories"] = F(
    "Scheme 3. Proposal Mechanism", [], PAGE7,
)
reaction["location"] = F(
    ["page 5 prose describing Scheme 3", "page 7 Scheme 3. Proposal Mechanism"],
    {
        "pages": [5, 7],
        "regions": ["mechanistic prose", "Scheme 3 proposed catalytic cycle"],
    },
    PAGE5_7,
)
reaction["reaction_name"] = F(
    "Scheme 3. Proposal Mechanism",
    "proposed mechanism for decarboxylative hydroaminoalkylation",
    PAGE7,
)
reaction["transformation_type"] = F(
    None,
    "author-proposed dual palladium/photoredox catalytic cycle",
    PAGE5_7,
    "inference",
    "The source caption explicitly labels the drawing as a proposal mechanism, and page 5 describes it as a potential mechanism.",
)
reaction["participants"] = [
    participant("substrate_1a", "m_1a", "reactant"),
    participant("allene_2", "m_scheme3_allene_2", "reactant"),
    participant("intermediate_A", "m_scheme3_intermediate_A", "unknown"),
    participant("intermediate_B", "m_scheme3_intermediate_B", "unknown"),
    participant("intermediate_C", "m_scheme3_intermediate_C", "unknown"),
    participant("intermediate_D", "m_scheme3_intermediate_D", "unknown"),
    participant("ester_5aa", "m_scheme3_ester_5aa", "product", "proposed competing ester intermediate"),
    participant("product_3", "m_scheme3_product_3", "product", "proposed generic hydroaminoalkylation product"),
]
reaction["steps"] = [
    proposed_step(
        "proposal_1_pd_hydride",
        ["substrate_1a", "intermediate_A"],
        "1a + L_nPd(0) -> A",
        "proposed oxidative addition/Pd(II)-hydride formation from amino acid 1a",
    ),
    proposed_step(
        "proposal_2_allene_insertion",
        ["intermediate_A", "allene_2", "intermediate_B"],
        "A + 2 -> B",
        "proposed allene migratory insertion to pi-allyl palladium intermediate B",
    ),
    proposed_step(
        "proposal_3_ester_branch",
        ["intermediate_B", "ester_5aa"],
        "B -> 5aa + Pd(0)",
        "proposed competing reductive elimination to ester intermediate 5aa",
    ),
    proposed_step(
        "proposal_4_photoredox_set",
        ["intermediate_B", "intermediate_C"],
        "Scheme 3 shows Ir(III)* / Ir(II) electron-transfer coupling to the B/C region",
        "proposed reductive-quenching/electron-transfer branch; exact arrow-to-intermediate assignment remains schematic",
    ),
    proposed_step(
        "proposal_5_decarboxylation_rebound",
        ["intermediate_B", "intermediate_C", "intermediate_D"],
        "B/C -> D; -e-, -CO2",
        "proposed single-electron decarboxylation and radical rebound giving palladium(III) intermediate D",
    ),
    proposed_step(
        "proposal_6_cc_reductive_elimination",
        ["intermediate_D", "product_3"],
        "D -> 3 + L_nPd(I)",
        "proposed C-C-forming reductive elimination to generic product 3 and Pd(I)",
    ),
    proposed_step(
        "proposal_7_catalyst_regeneration",
        ["intermediate_D"],
        "L_nPd(I) / L_nPd(0) and Ir(II) / Ir(III) regeneration arrows",
        "proposed reduction of Pd(I) by Ir(II), regenerating Pd(0) and Ir(III)",
    ),
]
reaction["outcomes"] = []
reaction["classification"] = {
    "standard_id": st.STANDARD,
    "mapped_reaction_smiles": F(),
    "bond_changes": [],
    "labels": F(
        None, [], PAGE5_7, "inference",
        "Scheme 3 is an author-proposed catalytic cycle, not an independently executed net experiment; no formation taxonomy labels are indexed for it.",
    ),
    "coverage": F(
        None, "not_applicable", PAGE5_7, "inference",
        "Bond-formation classification applies to observed net experiments, not proposed elementary intermediates.",
    ),
    "limitations": [
        "Author-proposed mechanism; intermediates A-D are not reported as isolated products on pages 7-8.",
        "R and L_n are generic variables, so no complete atom mapping or unique concrete SMILES is asserted.",
        "The B/C electron-transfer arrow topology is schematic; page-5 prose and page-7 drawing are retained together without inventing a more precise sequence.",
    ],
    "model_votes": [],
}
reaction["notes"] = [
    F(
        "Scheme 3. Proposal Mechanism",
        "The entire cycle is one proposed_mechanism record; individual arrows are ordered proposal steps, not separate experimental reactions.",
        PAGE7,
    ),
    F(
        "light-bulb symbol and Ir(III)* / Ir(II) / Ir(III) cycle",
        "A photoredox redox cycle is drawn, but wavelength, power, quantities, and a standalone experimental procedure are not reported in Scheme 3.",
        PAGE7,
    ),
]
reaction["evidence_ids"] = PAGE5_7
reaction["review"] = {
    "status": "needs_review",
    "method": "assistant_visual",
    "reviewer": None,
    "checked_fields": [
        "record_type", "location", "participants", "steps",
        "classification", "stereochemical labels in Scheme 3",
    ],
    "unresolved": [
        "Supporting Information and the CCDC structure file are explicitly mentioned but are not provided or associated in this workspace.",
        "Intermediates A-D and generic entities 2/3/5aa contain R and/or L_n variables and therefore have no unique concrete SMILES.",
        "Scheme 3 is proposed rather than direct experimental proof; no mechanistic arrow is promoted to a specific_experiment.",
    ],
}
reaction["completeness"] = "partial"

# Add page-7 evidence to the existing SI-source-limit context without associating an SI file.
si_context = next(c for c in bundle["source_context"]
                  if c["context_id"] == f"{DOC}:si_status")
si_context["note"] = (
    "The article explicitly states on page 7 that Supporting Information is available and lists experimental procedures, "
    "characterization, and X-ray data for (2S,3S)-3al (CCDC 1999447), but no SI or CCDC file was provided or explicitly "
    "associated in this workspace."
)
si_context["evidence_ids"] = list(dict.fromkeys(si_context["evidence_ids"] + PAGE7))

new_contexts = [
    {
        "context_id": f"{DOC}:page7_xray_3al",
        "document_id": DOC,
        "kind": "stereochemical_assignment_context",
        "note": "Associated Content states: X-ray data for (2S,3S)-3al, CCDC 1999447. The explicit assignment is attached to m_3al; the external structure file was not supplied.",
        "evidence_ids": PAGE7,
    },
    {
        "context_id": f"{DOC}:pages7_8_references",
        "document_id": DOC,
        "kind": "cited_literature",
        "note": "The references on pages 7-8 cite reviews and prior experiments; they are not experiments performed in the current paper and were not converted into current-paper reaction records.",
        "evidence_ids": [*PAGE7, *PAGE8],
    },
]
existing_context_ids = {c["context_id"] for c in bundle["source_context"]}
for context in new_contexts:
    if context["context_id"] not in existing_context_ids:
        bundle["source_context"].append(context)

st.write(RUN / "submitted.json", bundle)

raw = {
    "document_id": DOC,
    "pages": [7, 8],
    "source": "assistant visual review of rendered original PDF pages with extracted-text cross-check",
    "page7": {
        "scheme3": {
            "caption": "Scheme 3. Proposal Mechanism",
            "classification": "proposed_mechanism",
            "merged_record_id": reaction_id,
            "drawn_entities": ["1a", "2", "A", "B", "C", "D", "5aa", "3"],
            "decision": "One proposed-mechanism record with ordered arrow inventory; not split into experimental reactions.",
        },
        "associated_content": {
            "raw_stereo_statement": "(2S,3S)-3al (CCDC 1999447)",
            "molecule_id": "m_3al",
            "si_association": "unestablished",
        },
        "other_regions": ["author information", "author contributions", "notes", "acknowledgments", "references 1-4"],
    },
    "page8": {
        "regions": ["references 5-21"],
        "decision": "Cited literature only; no current-paper experiment, reaction arrow, table, or standalone characterization record.",
    },
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/pages_0007_0008_extraction.json", raw)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {
                    "region": "upper half: Scheme 3, Proposal Mechanism",
                    "finding": "Author-proposed Pd/photoredox catalytic cycle with 1a, generic allene 2, intermediates A-D, competing ester 5aa, product 3, Pd(0/I/II/III) and Ir(II/III) redox arrows. Merged into the existing proposed_mechanism record, not treated as experimental proof.",
                    "evidence_ids": PAGE7,
                },
                {
                    "region": "lower left: Associated Content / Supporting Information",
                    "finding": "The article states that SI contains experimental procedures, characterization, and X-ray data for (2S,3S)-3al (CCDC 1999447). The absolute assignment was attached to m_3al; no SI/CCDC file was associated.",
                    "evidence_ids": PAGE7,
                },
                {
                    "region": "lower left: Author Information, Author Contributions, Notes, Acknowledgments",
                    "finding": "Bibliographic and administrative content; no synthetic experiment or standalone characterization table.",
                    "evidence_ids": PAGE7,
                },
                {
                    "region": "lower right: References 1-4",
                    "finding": "Cited reviews/prior literature, not experiments performed in the current article.",
                    "evidence_ids": PAGE7,
                },
            ],
            "record_ids": [reaction_id],
            "unresolved": [
                "The SI and CCDC 1999447 structure file are mentioned but not supplied or explicitly associated in the workspace.",
                "Generic R and L_n groups prevent unique SMILES and complete stereochemical encoding for Scheme 3 entities.",
                "Intermediates A-D are proposed and not claimed as isolated or directly observed on this page.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page: References 5-21",
                    "finding": "Reference list continuation only. It includes cited reviews, prior synthetic methods, the SI cross-reference, an er-derivatization explanation, and mechanistic prior art; none is a new experiment performed in the current paper on this page.",
                    "evidence_ids": PAGE8,
                },
            ],
            "record_ids": [],
            "unresolved": [
                "Cited literature was not re-extracted as current-paper experiments; no reaction arrows, tables, or independent spectral peak lists appear on the page.",
            ],
        },
    ],
}
st.write(RUN / "coverage/pages_0007_0008.json", coverage)

print({
    "updated_reaction_records": [reaction_id],
    "new_generic_mechanism_molecules": [spec[0] for spec in generic_specs],
    "updated_stereo_molecule": "m_3al",
    "coverage_pages": [7, 8],
})
