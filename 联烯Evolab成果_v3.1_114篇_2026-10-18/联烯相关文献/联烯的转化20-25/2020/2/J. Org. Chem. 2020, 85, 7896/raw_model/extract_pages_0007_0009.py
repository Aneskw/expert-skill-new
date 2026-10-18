#!/usr/bin/env python3
"""Build the incremental v3.1 submission for source PDF pages 7-9.

This file is a raw machine extraction artifact.  It updates nine already
enumerated scope experiments with their page-7 characterization evidence and
does not create experiments from author information or references on pages
8-9.  All records remain machine candidates requiring expert review.
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


RUN = ROOT / "runs/current"
DOC = "doc_b298156f3445f448"
E6 = f"{DOC}:page6:text"
I6 = f"{DOC}:page6:image"
E7 = f"{DOC}:page7:text"
I7 = f"{DOC}:page7:image"
E8 = f"{DOC}:page8:text"
I8 = f"{DOC}:page8:image"
E9 = f"{DOC}:page9:text"
I9 = f"{DOC}:page9:image"


def fact(raw, normalized, evidence, method="original_extraction", derivation=None, status=None):
    return st.fact(
        raw,
        normalized,
        evidence,
        method=method,
        derivation=derivation,
        status=status,
    )


def append_unique(values, additions):
    for value in additions:
        if value not in values:
            values.append(value)


existing = json.loads((RUN / "standard_records.json").read_text())
template = json.loads((RUN / "raw_model/pages_0007_0009_init.json").read_text())
old_molecules = {m["molecule_id"]: m for m in existing["molecules"]}
old_reactions = {r["reaction_id"]: r for r in existing["reactions"]}


products = {
    "3s": {
        "name": "(R,E)-3-Allyl-1-(3-(benzyloxy)allyl)-1-methylnaphthalen-2(1H)-one",
        "smiles": "C=CCC1=Cc2ccccc2[C@@](C)(C/C=C/OCc2ccccc2)C1=O",
        "configuration": "R",
        "mass": 61.4,
        "appearance": "Colorless oil",
        "purification": "petroleum ether/EtOAc = 40:1-10:1",
        "er": "93:8",
        "hplc": "Chiralcel OD, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 8.11 min, minor 7.23 min",
        "reaction": "scheme2_3s",
        "name_evidence": [I6, E6],
        "mass_evidence": [I6, E6],
        "er_evidence": [I7, E7],
        "yield_repeated_on_page7": False,
    },
    "3t": {
        "name": "(R,E)-1-(3-(Benzyloxy)allyl)-1-methyl-3-phenylnaphthalen-2(1H)-one",
        "smiles": "C[C@]1(C/C=C/OCc2ccccc2)C(=O)C(c2ccccc2)=Cc2ccccc21",
        "configuration": "R",
        "mass": 57.2,
        "appearance": "Colorless oil",
        "purification": "petroleum ether/EtOAc = 40:1-10:1",
        "er": "92.5:7.5",
        "hplc": "Chiralcel AD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 13.21 min, minor 12.24 min",
        "reaction": "scheme2_3t",
    },
    "3u": {
        "name": "(R,E)-1-(3-(Benzyloxy)allyl)-1,3-dimethyl-6-phenylnaphthalen-2(1H)-one",
        "smiles": "CC1=Cc2ccc(-c3ccccc3)cc2[C@@](C)(C/C=C/OCc2ccccc2)C1=O",
        "configuration": "R",
        "mass": 56.1,
        "appearance": "Colorless oil",
        "purification": "petroleum ether/EtOAc = 40:1-10:1",
        "er": "91.5:8.5",
        "hplc": "Chiralcel IA, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 11.40 min, minor 12.27 min",
        "reaction": "scheme2_3u",
    },
    "5a": {
        "name": "(S,E)-2-Allyl-2-(3-(benzyloxy)allyl)naphthalen-1(2H)-one",
        "smiles": "C=CC[C@]1(C/C=C/OCc2ccccc2)C=Cc2ccccc2C1=O",
        "configuration": "S",
        "mass": 43.7,
        "appearance": "Yellow oil",
        "purification": "petroleum ether/EtOAc = 40:1-8:1",
        "er": "91.5:8.5",
        "hplc": "Chiralcel OD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 8.43 min, minor 7.50 min",
        "reaction": "scheme3_5a",
    },
    "5b": {
        "name": "(S,E)-2-Allyl-2-(3-((4-ethylbenzyl)oxy)allyl)naphthalen-1(2H)-one",
        "smiles": "C=CC[C@]1(C/C=C/OCc2ccc(CC)cc2)C=Cc2ccccc2C1=O",
        "configuration": "S",
        "mass": 44.4,
        "appearance": "Yellow oil",
        "purification": "petroleum ether/EtOAc = 40:1-8:1",
        "er": "92:8",
        "hplc": "Chiralcel IA, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 11.16 min, minor 10.30 min",
        "reaction": "scheme3_5b",
    },
    "5c": {
        "name": "(S,E)-2-Allyl-2-(3-((3,5-dimethoxybenzyl)oxy)allyl)naphthalen-1(2H)-one",
        "smiles": "C=CC[C@]1(C/C=C/OCc2cc(OC)cc(OC)c2)C=Cc2ccccc2C1=O",
        "configuration": "S",
        "mass": 55.4,
        "appearance": "Yellow oil",
        "purification": "petroleum ether/EtOAc = 40:1-8:1",
        "er": "92:8",
        "hplc": "Chiralcel AD, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 11.41 min, minor 8.96 min",
        "reaction": "scheme3_5c",
    },
    "5d": {
        "name": "(S,E)-2-Allyl-2-(3-((4-fluorobenzyl)oxy)allyl)naphthalen-1(2H)-one",
        "smiles": "C=CC[C@]1(C/C=C/OCc2ccc(F)cc2)C=Cc2ccccc2C1=O",
        "configuration": "S",
        "mass": 45.2,
        "appearance": "Yellow oil",
        "purification": "petroleum ether/EtOAc = 40:1-8:1",
        "er": "91.5:8.5",
        "hplc": "Chiralcel AD, 99:1 hexanes/i-PrOH, 1 mL/min; tR major 21.04 min, minor 17.95 min",
        "reaction": "scheme3_5d",
    },
    "5e": {
        "name": "(S,E)-2-Allyl-2-(3-(benzyloxy)allyl)-4-methoxynaphthalen-1(2H)-one",
        "smiles": "C=CC[C@]1(C/C=C/OCc2ccccc2)C=C(OC)c2ccccc2C1=O",
        "configuration": "S",
        "mass": 56.9,
        "appearance": "Yellow oil",
        "purification": "petroleum ether/EtOAc = 40:1-8:1",
        "er": "94:6",
        "hplc": "Chiralcel OD, 90:10 hexanes/i-PrOH, 1 mL/min; tR major 9.44 min, minor 5.12 min",
        "reaction": "scheme3_5e",
    },
    "5f": {
        "name": "(S,E)-2-Allyl-2-(3-(benzyloxy)allyl)-4-ethoxynaphthalen-1(2H)-one",
        "smiles": "C=CC[C@]1(C/C=C/OCc2ccccc2)C=C(OCC)c2ccccc2C1=O",
        "configuration": "S",
        "mass": 59.1,
        "appearance": "Yellow oil",
        "purification": "petroleum ether/EtOAc = 20:1-4:1",
        "er": "93:7",
        "hplc": "Chiralcel OD, 95:5 hexanes/i-PrOH, 1 mL/min; tR major 6.60 min, minor 5.83 min",
        "reaction": "scheme3_5f",
    },
}


stereo_checks = {}
incoming_molecules = []
incoming_reactions = []

for label, spec in products.items():
    molecule_id = f"{DOC}:mol:{label}"
    reaction_id = f"{DOC}:{spec['reaction']}"
    molecule = copy.deepcopy(old_molecules[molecule_id])
    reaction = copy.deepcopy(old_reactions[reaction_id])
    name_evidence = spec.get("name_evidence", [I7, E7])
    mass_evidence = spec.get("mass_evidence", [I7, E7])
    er_evidence = spec.get("er_evidence", [I7, E7])
    combined_evidence = list(dict.fromkeys(name_evidence + er_evidence))

    diag = st.stereo_inspect(spec["smiles"])
    stereo_checks[molecule_id] = diag
    centers = diag["atom_centers"]
    named_center = next((c for c in centers if c.get("CIP") == spec["configuration"]), None)
    if named_center is None:
        raise ValueError(f"Candidate stereochemistry for {label} does not encode {spec['configuration']}")
    e_bond = next((b for b in diag["double_bonds"] if b.get("stereo") == "STEREOE"), None)
    if e_bond is None:
        raise ValueError(f"Candidate stereochemistry for {label} does not encode E geometry")

    occurrence = next(o for o in molecule["occurrences"] if o["document_id"] == DOC)
    append_unique(occurrence["evidence_ids"], combined_evidence)
    occurrence["compound_label"] = fact(label, label, name_evidence)
    molecule["names"]["original"] = fact(spec["name"], spec["name"], name_evidence)
    molecule["names"]["english"] = fact(spec["name"], spec["name"], name_evidence)
    structure_basis = (
        "Candidate connectivity was retained from the scheme drawing and source systematic name; "
        "the named R/S center and E alkene were encoded and checked with standard stereo. "
        "Toolkit parsing is not source or chemical verification."
    )
    molecule["representation_type"] = "smiles"
    molecule["smiles"] = fact(
        spec["name"],
        diag["canonical_isomeric_smiles"],
        name_evidence,
        "structure_derivation",
        structure_basis,
    )
    molecule["stereochemistry"]["atom_centers"] = fact(
        f"{spec['configuration']} in the source compound name",
        [
            {
                "atom_index": named_center["atom_index"],
                "configuration": spec["configuration"],
                "index_basis": "zero-based atom index in canonical isomeric SMILES",
            }
        ],
        name_evidence,
        "structure_derivation",
        structure_basis,
    )
    molecule["stereochemistry"]["double_bonds"] = fact(
        "E in the source compound name",
        [
            {
                "bond_index": e_bond["bond_index"],
                "configuration": "E",
                "index_basis": "RDKit bond index in canonical isomeric SMILES",
            }
        ],
        name_evidence,
        "structure_derivation",
        structure_basis,
    )
    molecule["stereochemistry"]["axial"] = fact(
        "No allene or separately named stereogenic axis in the product name/drawing",
        [],
        name_evidence,
        "structure_derivation",
        "The product candidate contains no C=C=C unit; absence of an encoded axis is limited to this product graph.",
    )
    molecule["stereochemistry"]["mixture"] = fact(
        f"Isolated sample reported at {spec['er']} e.r.; individual retention-peak configurations are not assigned",
        {
            "metric": "e.r.",
            "ratio": spec["er"],
            "order": "major:minor",
            "retention_peak_configurations": "unassigned",
        },
        er_evidence,
    )
    molecule["stereochemistry"]["assignment"] = fact(
        f"{spec['configuration']},E in the compound name; no ring cis/trans descriptor reported",
        {
            "atom_center": spec["configuration"],
            "alkene_geometry": "E",
            "ring_cis_trans": "not applicable to the depicted single-center product",
            "assignment_method": "not reported in this page-7 characterization paragraph",
        },
        name_evidence,
        "structure_derivation",
        "The R/S and E descriptors are source-derived; ring cis/trans applicability is assessed from the candidate graph.",
    )
    molecule["stereochemistry"]["chirality_status"] = fact(
        spec["configuration"], "specified", name_evidence
    )
    molecule["stereochemistry"]["geometry_status"] = fact("E", "specified", name_evidence)
    molecule["stereochemistry"]["representation_status"] = fact(
        f"Candidate isomeric SMILES encodes named {spec['configuration']},E stereochemistry",
        "encoded",
        name_evidence,
        "structure_derivation",
        structure_basis,
    )
    molecule["stereochemistry"]["isomeric_smiles"] = fact(
        spec["name"],
        diag["canonical_isomeric_smiles"],
        name_evidence,
        "structure_derivation",
        structure_basis,
    )
    molecule["stereochemistry"]["identity_key"] = fact(
        None,
        diag["stereo_graph_key"],
        name_evidence,
        "structure_derivation",
        "Computed from the candidate isomeric molecular graph; not a unique sample identifier.",
    )
    molecule["identity_note"] = fact(
        f"Source label {label} and source name {spec['name']}; structure remains a machine candidate.",
        f"Source-identified {label}; candidate structure and stereochemistry require expert review.",
        name_evidence,
    )

    append_unique(reaction["evidence_ids"], combined_evidence)
    # Page 7 corroborates the already enumerated yield and gives the product-specific
    # chiral HPLC method.  It is not a second synthesis.
    yield_outcome = next(o for o in reaction["outcomes"] if o["outcome_id"] == "yield")
    if spec.get("yield_repeated_on_page7", True):
        append_unique(yield_outcome["value"]["evidence_ids"], [I7, E7])
    er_outcome = next(o for o in reaction["outcomes"] if o["outcome_id"] == "er")
    er_outcome["value"] = fact(
        f"{spec['er']} e.r.",
        {"metric": "e.r.", "ratio": spec["er"], "order": "major:minor"},
        er_evidence,
    )
    er_outcome["measurement_basis"] = fact(spec["hplc"], spec["hplc"], er_evidence)
    er_outcome["scope"] = fact(
        "e.r.; tR reported as major and minor without assigning either retention peak to an absolute configuration",
        {
            "metric": "e.r.",
            "order": "major:minor",
            "target_molecule_id": molecule_id,
            "retention_peak_configurations": "unassigned",
        },
        er_evidence,
    )
    mass_outcome = st.new_outcome("mass", "product_amount", [molecule_id], mass_evidence)
    mass_outcome["value"] = fact(
        f"{spec['mass']} mg",
        {"number": spec["mass"], "unit": "mg"},
        mass_evidence,
    )
    mass_outcome["measurement_basis"] = fact(
        f"{spec['appearance']}; {spec['purification']}; {spec['mass']} mg, yield reported in the characterization paragraph",
        "reported product mass in the preparation/characterization paragraph",
        mass_evidence,
    )
    mass_outcome["scope"] = fact(
        "product amount for the same scope experiment; page 7 is a characterization continuation, not a second reaction",
        "product amount for the same scope experiment",
        mass_evidence,
    )
    reaction["outcomes"] = [o for o in reaction["outcomes"] if o["outcome_id"] != "mass"]
    reaction["outcomes"].append(mass_outcome)
    # Make reruns idempotent: replace this pass's explanatory note rather than
    # accumulating a second copy when the script reads an already merged dataset.
    reaction["notes"] = [
        note
        for note in reaction["notes"]
        if not str(note.get("raw_value") or "").startswith("Page 7")
    ]
    reaction["notes"].append(
        fact(
            (
                "Page 7 continues the 3s characterization paragraph with chiral HPLC and HRMS."
                if label == "3s"
                else "Page 7 contains the preparation/characterization paragraph (mass, repeated yield, optical rotation, NMR, chiral HPLC, HRMS)."
            ),
            (
                "The product-specific chiral HPLC was merged into the existing 3s scope record; its mass/name remain evidenced on page 6, and HRMS remains source context."
                if label == "3s"
                else "Mass and product-specific chiral HPLC were merged into this existing scope record; spectral/HRMS data remain source context and do not create another synthesis."
            ),
            combined_evidence,
            "inference",
            "The paragraph has no new reaction arrow or changed conditions and repeats the Scheme 2/3 product label and yield.",
        )
    )
    reaction["review"]["status"] = "needs_review"
    reaction["review"]["method"] = "model"
    reaction["review"]["reviewer"] = None
    append_unique(
        reaction["review"]["unresolved"],
        [
            "Page-7 NMR, optical rotation, and HRMS were retained as characterization context and were not independently checked peak-by-peak.",
            "The chiral HPLC retention peaks are labeled major/minor, but the page does not explicitly assign each peak to an absolute configuration.",
            "The isomeric SMILES is a structure-derived machine candidate; standard stereo parsing does not establish chemical correctness.",
        ],
    )
    reaction["completeness"] = "partial"

    incoming_molecules.append(molecule)
    incoming_reactions.append(reaction)


# A standard bundle is independently valid, so include unchanged generic-scheme
# dependencies and every molecule referenced by the nine updated reactions.
# Merge replaces records by stable IDs, therefore these copies do not create new
# experiments or duplicate the existing records.
dependency_reactions = [
    copy.deepcopy(old_reactions[f"{DOC}:scheme2_generic"]),
    copy.deepcopy(old_reactions[f"{DOC}:scheme3_generic"]),
]
bundle_reactions = dependency_reactions + incoming_reactions
updated_molecule_map = {m["molecule_id"]: m for m in incoming_molecules}
required_molecule_ids = {
    participant["molecule_id"]
    for reaction in bundle_reactions
    for participant in reaction["participants"]
}
template["molecules"] = [
    copy.deepcopy(updated_molecule_map.get(mid, old_molecules[mid]))
    for mid in sorted(required_molecule_ids)
]
template["reactions"] = bundle_reactions
template["source_context"] = [
    {
        "context_id": f"{DOC}:context:characterization_page7",
        "kind": "characterization",
        "evidence_ids": [I7, E7],
        "value": (
            "The opening lines continue 3s chiral-HPLC/HRMS characterization, followed by preparation/characterization paragraphs for 3t, 3u, and 5a-5f. Product masses and "
            "product-specific chiral-HPLC methods are merged into the existing Scheme 2/3 scope records. "
            "Optical rotations, NMR peaks, HRMS values, and the opening continuation of 3s are context, not separate syntheses."
        ),
    },
    {
        "context_id": f"{DOC}:context:si_notice_page7",
        "kind": "supporting_information_notice",
        "evidence_ids": [I7, E7],
        "value": (
            "The article advertises a Supporting Information PDF containing spectra, but no SI file was explicitly provided or associated; "
            "no SI evidence is used in this extraction."
        ),
    },
    {
        "context_id": f"{DOC}:context:author_information_page8",
        "kind": "bibliographic_metadata",
        "evidence_ids": [I8, E8],
        "value": "Author information, conflict-of-interest statement, acknowledgments, and the beginning of the reference list; no current-paper experiment is depicted.",
    },
    {
        "context_id": f"{DOC}:context:references_pages8_9",
        "kind": "references",
        "evidence_ids": [I8, E8, I9, E9],
        "value": "Bibliographic references only. Cited paper titles are not treated as experiments performed in the current article.",
    },
]


candidate_path = RUN / "raw_model/pages_0007_0009_submission.json"
st.write(candidate_path, template)
st.write(RUN / "raw_model/pages_0007_0009_stereo_checks.json", stereo_checks)
st.write(RUN / "submitted.json", template)
print(
    json.dumps(
        {
            "pages": [7, 8, 9],
            "updated_reactions": len(incoming_reactions),
            "updated_molecules": len(incoming_molecules),
            "new_reactions": 0,
            "unchanged_dependency_reactions": len(dependency_reactions),
            "source_context_entries": len(template["source_context"]),
        },
        ensure_ascii=False,
    )
)
