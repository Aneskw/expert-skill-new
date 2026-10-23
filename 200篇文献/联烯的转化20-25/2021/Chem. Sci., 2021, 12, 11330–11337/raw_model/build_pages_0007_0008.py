#!/usr/bin/env python3
"""Build the page-7--8 v3.1 incremental submission from direct visual transcription."""
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
DOC = "doc_b620dbe6e0542262"
PREFIX = "d1sc02167b"
raw = json.loads((RUN / "raw_model/page7_visual_transcription.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
bundle = json.loads((RUN / "raw_model/pages_0007_0008_template.json").read_text())


def evid(page: int):
    return [f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"]


def direct(raw_value, normalized=None, page=7):
    if normalized is None:
        normalized = raw_value
    return st.fact(raw_value, normalized, evid(page))


def derived(normalized, derivation, raw_value=None, page=7):
    return st.fact(raw_value, normalized, evid(page), method="structure_derivation", derivation=derivation)


def inferred(normalized, derivation, page=7):
    return st.fact(None, normalized, evid(page), method="inference", derivation=derivation)


def missing(page=7):
    return st.fact(evidence=evid(page))


def mid(label: str):
    return f"{PREFIX}_p7_mol_{label.lower().replace('-', '_').replace('·', '_').replace('(', '').replace(')', '').replace('/', '_').replace(' ', '_')}"


def rid(label: str):
    return f"{DOC}:{label}"


CANONICAL = {
    "N-methylmaleimide": ("CN1C(=O)C=CC1=O", "190263ccdd57d896472650be8aded4f146a3dbe55a59ca59cdcd0c68e4251ce4"),
    "diethyl acetylenedicarboxylate": ("CCOC(=O)C#CC(=O)OCC", "551fa135c408030686071e2028a73c2f9108f14ce2d4e4c817689f822b40472e"),
    "NaOEt": ("CC[O-].[Na+]", "31e739e8ca19eb08fe7f478d05c9356d240684b3ae53d70432fad169749a0e90"),
    "DMSO": ("CS(C)=O", "ef4a329fe2b9444e2ae17d5a7d6b0fd38e35d5eae829d6ca83527841c4d3898f"),
    "PhB(OH)2": ("OB(O)c1ccccc1", "8688c09a3a54367bd40d7d2c37c9f8fe5fca8d30504eb1c11832e4077cfc95c2"),
    "K3PO4": ("O=P([O-])([O-])[O-].[K+].[K+].[K+]", "a11342c3e792c6bea3f5b6732e26cfb995c153795495fa8ce3550ec8efd7c9bd"),
    "DME": ("COCCOC", "ffa88e79e5e17aaae45e234d5a2f2016cca9140f2e2d1be6cb0e7159fb0f6bba"),
    "H2O": ("O", "c4694f2e93d5c4e7d51f9c5deb75e6cc8be5e1114178c6a45b6fc2c566a0aa8c"),
    "methyl acrylate": ("C=CC(=O)OC", "7eef9856e5d090f08a9653a9e2a7e84e712302799ef5b58c552caf9661c56dd7"),
    "DMF": ("CN(C)C=O", "aeb2e86e450d0fab9900120c06f6f4a77a774c09aacea28aa2f6ecfa8650bc7f"),
}

NAMES = {
    "N-methylmaleimide": "N-methylmaleimide",
    "diethyl acetylenedicarboxylate": "diethyl acetylenedicarboxylate",
    "NaOEt": "sodium ethoxide",
    "DMSO": "dimethyl sulfoxide",
    "PhB(OH)2": "phenylboronic acid",
    "Pd(PPh3)4": "tetrakis(triphenylphosphine)palladium(0)",
    "K3PO4": "tripotassium phosphate",
    "DME": "1,2-dimethoxyethane",
    "H2O": "water",
    "methyl acrylate": "methyl acrylate",
    "[Cp*RhCl2]2": "[Cp*RhCl2]2",
    "Cu(OAc)2·H2O": "copper(II) acetate monohydrate",
    "DMF": "N,N-dimethylformamide",
    "Sc(OTf)3": "scandium triflate",
    "E-3an": "source compound E-3an",
    "14": "reported product 14",
    "15": "reported product 15",
    "16": "reported product 16 (Z/E mixture)",
    "17": "reported product 17 (Z/E mixture)",
    "18": "reported product 18 (Z/E mixture)",
    "cis-19": "reported product cis-19",
    "E-20": "reported product E-20",
}


def base_molecule(label: str):
    m = st.new_molecule(mid(label), DOC, evid(7))
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label)
    name = NAMES[label]
    m["names"]["original"] = direct(label, label)
    m["names"]["english"] = direct(name, name)
    m["identity_note"] = direct(
        f"Scheme 2 identifies this species as {label}; structure and stereochemistry are encoded only to the extent supported by the page image.",
        f"Scheme 2 identity {label}; machine record needs review",
    )
    return m


def simple_molecule(label: str):
    m = base_molecule(label)
    smiles, key = CANONICAL[label]
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, "Candidate connectivity derived from the explicit Scheme 2 reagent/solvent name and checked with standard stereo; parsing does not independently verify source identity.", f"Scheme 2 label {label}")
    m["stereochemistry"]["atom_centers"] = derived([], "The standard stereo diagnostic found no potential tetrahedral atom centers in this candidate graph.")
    m["stereochemistry"]["double_bonds"] = derived([], "No configurable alkene E/Z element was found in this candidate graph; carbonyl or constrained ring double bonds are not treated as E/Z assignments.")
    m["stereochemistry"]["axial"] = derived([], "No allene or other axial stereogenic unit occurs in the candidate graph.")
    m["stereochemistry"]["chirality_status"] = derived("achiral", "The candidate graph has no potential tetrahedral or axial stereogenic element.")
    m["stereochemistry"]["geometry_status"] = derived("not_applicable", "The candidate graph has no configurable alkene geometry.")
    m["stereochemistry"]["representation_status"] = derived("encoded", "The candidate graph contains no source-specified stereochemical element omitted from the isomeric SMILES.")
    m["stereochemistry"]["isomeric_smiles"] = derived(smiles, "Canonical isomeric SMILES returned by the standard stereo diagnostic.")
    m["stereochemistry"]["identity_key"] = derived(key, "Stereo graph key returned by the standard stereo diagnostic; it is not a sample identity proof.")
    return m


def unresolved_molecule(label: str, *, chirality="unresolved", geometry="unresolved", representation="unresolved", assignment=None, mixture=None):
    m = base_molecule(label)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["stereochemistry"]["chirality_status"] = inferred(chirality, "Center/relative stereochemistry was inspected in the Scheme 2 drawing but cannot be fully encoded or assigned from this page alone.")
    m["stereochemistry"]["geometry_status"] = inferred(geometry, "Alkene or ring-relative geometry is retained from the source label/result when explicit; remaining geometry is unresolved.")
    m["stereochemistry"]["representation_status"] = inferred(representation, "No fully verified stereochemical SMILES was asserted from the compact Scheme 2 drawing.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = direct(assignment, assignment)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = direct(mixture, mixture)
    return m


molecules = []
for label in CANONICAL:
    molecules.append(simple_molecule(label))
for label in ("Pd(PPh3)4", "[Cp*RhCl2]2", "Cu(OAc)2·H2O", "Sc(OTf)3"):
    molecules.append(unresolved_molecule(label, chirality="unresolved", geometry="not_applicable", representation="unresolved"))
molecules += [
    unresolved_molecule("E-3an", geometry="specified", representation="not_encoded", assignment={"alkene_configuration": "E", "source_label": "E-3an", "locants": "not restated in Scheme 2"}),
    unresolved_molecule("14", chirality="mixture", geometry="unresolved", representation="not_encoded", assignment={"major_diastereomer": "structure drawn", "descriptor": "not assigned"}, mixture={"type": "diastereomers", "composition": ">20:1", "order": ["major", "minor"]}),
    unresolved_molecule("15"),
    unresolved_molecule("16", geometry="mixture", representation="not_encoded", mixture={"type": "geometric_isomers", "metric": "Z/E", "composition": "96/4"}),
    unresolved_molecule("17", geometry="mixture", representation="not_encoded", mixture={"type": "geometric_isomers", "metric": "Z/E", "composition": "93/7"}),
    unresolved_molecule("18", geometry="mixture", representation="not_encoded", mixture={"type": "geometric_isomers", "metric": "Z/E", "composition": "95/5"}),
    unresolved_molecule("cis-19", chirality="partially_specified", geometry="unresolved", representation="not_encoded", assignment={"ring_relative_configuration": "cis", "source_label": "cis-19", "absolute_descriptors": "not reported"}),
    unresolved_molecule("E-20", geometry="specified", representation="not_encoded", assignment={"alkene_configuration": "E", "source_label": "E-20", "locants": "not restated in Scheme 2"}),
]

# Reuse only source molecule identities needed to prevent duplicate registration;
# no prior experiment values are used as page-7 results.
old_molecules = {m["molecule_id"]: m for m in dataset.get("standard", {}).get("molecules", [])}
for existing_id in (f"{PREFIX}_mol_3bf", f"{PREFIX}_mol_3am"):
    molecules.append(copy.deepcopy(old_molecules[existing_id]))

bundle["molecules"] = molecules


def participant(pid, molecule_id, role, *, stage="reaction", ref=None):
    p = st.new_participant(pid, molecule_id, role, evid(7), stage=stage)
    if role == "product":
        p["product_role"] = direct("reported product", "main_product")
    if ref:
        p["equivalence_reference"] = direct(f"relative to {ref}", ref)
    return p


def add_charge(p, key, raw_value, number, unit):
    p["charges"][key] = direct(raw_value, {"number": number, "unit": unit})


def outcome(r, oid, kind, target, raw_value, normalized, *, yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, [target], evid(7))
    o["value"] = direct(raw_value, normalized)
    if kind in {"yield", "product_amount"}:
        o["product_role"] = direct("reported product", "main_product")
    if yield_type is not None:
        o["yield_type"] = direct(yield_type, yield_type)
    if basis is not None:
        o["measurement_basis"] = direct(basis, basis)
    if scope is not None:
        o["scope"] = direct(scope, scope)
    r["outcomes"].append(o)


def classify(r, changes, labels, limitation):
    r["classification"]["bond_changes"] = [
        derived(change, "Net bond changes inferred from the Scheme 2 reactant/product drawings without a complete atom mapping; atom identities remain partial.")
        for change in changes
    ]
    r["classification"]["labels"] = derived(labels, "Formation labels follow the new connections and bond-order increases visible in Scheme 2; unchanged and bond-order-decreased bonds are excluded from labels.")
    r["classification"]["coverage"] = inferred("partial", limitation)
    r["classification"]["limitations"] = [limitation, "No complete atom-mapped reaction SMILES was asserted."]


def atom_change(cid, a, b, before, after, ctype):
    return {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}], "before_order": before, "after_order": after, "change_type": ctype}


def base_reaction(suffix, source_name, transformation, reactants, product, *, solvent=None, temp=None, time_h=None, atmosphere=None):
    r = st.new_reaction(rid(suffix), DOC, evid(7))
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = direct("derivatization", ["derivatization"])
    r["location"] = direct(f"PDF page 7, Scheme 2, arrow to {product}", {"page": 7, "region": f"Scheme 2 arrow to {product}"})
    r["reaction_name"] = direct(source_name, source_name)
    r["transformation_type"] = inferred(transformation, "Normalized from the Scheme 2 arrow and product drawing.")
    ps = []
    for i, molecule_id in enumerate(reactants, 1):
        ps.append(participant(f"p_reactant_{i}", molecule_id, "reactant"))
    product_id = mid(product)
    ps.append(participant("p_product", product_id, "product"))
    if solvent:
        solvent_ids = solvent if isinstance(solvent, list) else [solvent]
        for i, sid in enumerate(solvent_ids, 1):
            ps.append(participant(f"p_solvent_{i}", sid, "solvent"))
        r["solvent_conditions"]["composition"] = direct(
            "/".join(NAMES_BY_ID[s] for s in solvent_ids), solvent_ids
        )
    r["participants"] = ps
    step = st.new_step("reaction", evid(7))
    step["participant_ids"] = [p["participant_id"] for p in ps if p["role"] != "product"]
    step["operation"] = direct(f"Scheme 2 transformation to {product}", f"conduct reaction under displayed Scheme 2 conditions to give {product}")
    if isinstance(temp, (int, float)):
        step["temperature"] = direct(f"{temp} °C", {"number": temp, "unit": "°C"})
    elif temp == "r.t.":
        step["temperature"] = direct("r.t.", "r.t.")
    if time_h is not None:
        step["time"] = direct(f"{time_h} h", {"number": time_h, "unit": "h"})
    if atmosphere:
        step["atmosphere"] = direct(atmosphere, "argon" if atmosphere == "Ar" else atmosphere)
    r["steps"] = [step]
    r["procedure"]["original_text"] = direct(f"Scheme 2 arrow to {product}: displayed reagents/solvent, temperature, time, and atmosphere as separately structured.", f"Scheme 2 arrow to {product}")
    r["procedure"]["general_procedure_reference"] = missing()
    r["review"]["unresolved"] = [
        "No explicitly associated SI was supplied; reaction scale, addition order, workup, purification, and full characterization remain unreported here.",
        "The percentage yield is displayed, but its isolated/NMR measurement basis is not stated in Scheme 2.",
        "Compact product drawings do not justify a fully verified stereochemical SMILES; machine structure and bond-change assignments require expert review.",
    ]
    r["completeness"] = "partial"
    return r


NAMES_BY_ID = {mid(label): label for label in NAMES}
NAMES_BY_ID.update({f"{PREFIX}_mol_3bf": "3bf", f"{PREFIX}_mol_3am": "3am"})

reactions = []

# 3bf + N-methylmaleimide -> 14
r = base_reaction("scheme2_14", "Diels–Alder reaction of 3bf with N-methylmaleimide", "Diels–Alder cycloaddition", [f"{PREFIX}_mol_3bf", mid("N-methylmaleimide")], "14", solvent=mid("DMF"), temp=60, time_h=24, atmosphere="Ar")
add_charge(r["participants"][1], "equivalents", "1 equiv", 1, "equiv")
r["participants"][1]["equivalence_reference"] = direct("relative to 3bf", f"{PREFIX}_mol_3bf")
outcome(r, "o_yield", "yield", mid("14"), "14, 80%", {"number": 80, "unit": "%"}, yield_type="reported yield")
outcome(r, "o_dr", "selectivity", mid("14"), "dr > 20:1", {"relation": ">", "major": 20, "minor": 1, "raw_text": ">20:1"}, scope={"metric": "dr", "order": ["major", "minor"], "major_configuration": "not assigned by descriptor"})
r["notes"].append(direct("Starting material shown as 3bf, Z/E = 95/5.", {"input": "3bf", "metric": "Z/E", "value": "95/5"}))
classify(r, [atom_change("b1", "C", "C", 0, 1, "formed"), atom_change("b2", "C", "C", 0, 1, "formed"), atom_change("b3", "C", "C", 1, 2, "order_increased"), atom_change("b4", "C", "C", 2, 1, "order_decreased"), atom_change("b5", "C", "C", 2, 1, "order_decreased")], [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}, {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]}], "Diels–Alder net changes are inferred from the drawings; atom-level correspondence and implicit-H changes are not fully mapped.")
reactions.append(r)

# 3bf + diethyl acetylenedicarboxylate -> 15
r = base_reaction("scheme2_15", "Diels–Alder reaction of 3bf with diethyl acetylenedicarboxylate", "Diels–Alder cycloaddition", [f"{PREFIX}_mol_3bf", mid("diethyl acetylenedicarboxylate")], "15", solvent=mid("DMF"), temp=60, time_h=24, atmosphere="Ar")
add_charge(r["participants"][1], "equivalents", "3 equiv", 3, "equiv")
r["participants"][1]["equivalence_reference"] = direct("relative to 3bf", f"{PREFIX}_mol_3bf")
outcome(r, "o_yield", "yield", mid("15"), "15, 77%", {"number": 77, "unit": "%"}, yield_type="reported yield")
r["notes"].append(direct("Starting material shown as 3bf, Z/E = 95/5.", {"input": "3bf", "metric": "Z/E", "value": "95/5"}))
classify(r, [atom_change("b1", "C", "C", 0, 1, "formed"), atom_change("b2", "C", "C", 0, 1, "formed"), atom_change("b3", "C", "C", 1, 2, "order_increased"), atom_change("b4", "C", "C", 2, 1, "order_decreased"), atom_change("b5", "C", "C", 2, 1, "order_decreased"), atom_change("b6", "C", "C", 3, 2, "order_decreased")], [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}, {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]}], "Cycloaddition net changes are inferred from the drawings; the alkyne-to-alkene decrease is retained but is not itself a new double-bond label.")
reactions.append(r)

# Suzuki coupling -> 16
r = base_reaction("scheme2_16", "Suzuki coupling of 3bf with PhB(OH)2", "Suzuki aryl–aryl cross-coupling", [f"{PREFIX}_mol_3bf", mid("PhB(OH)2")], "16", solvent=[mid("DME"), mid("H2O")], temp=65, time_h=18, atmosphere="Ar")
add_charge(r["participants"][1], "equivalents", "1.2 equiv", 1.2, "equiv")
r["participants"][1]["equivalence_reference"] = direct("relative to 3bf", f"{PREFIX}_mol_3bf")
p = participant("p_palladium", mid("Pd(PPh3)4"), "catalyst", ref=f"{PREFIX}_mol_3bf"); add_charge(p, "mol_percent", "2 mol%", 2, "mol%"); r["participants"].append(p)
p = participant("p_base", mid("K3PO4"), "reagent", ref=f"{PREFIX}_mol_3bf"); add_charge(p, "equivalents", "2 equiv", 2, "equiv"); r["participants"].append(p)
r["steps"][0]["participant_ids"] += ["p_palladium", "p_base"]
r["solvent_conditions"]["mixture_ratio"] = direct("DME/H2O = 4:1", {"components": [mid("DME"), mid("H2O")], "ratio": [4, 1]})
r["solvent_conditions"]["total_volume"] = direct("2 mL", {"number": 2, "unit": "mL"})
outcome(r, "o_yield", "yield", mid("16"), "16, 89%", {"number": 89, "unit": "%"}, yield_type="reported yield")
outcome(r, "o_selectivity", "selectivity", mid("16"), "Z/E = 96/4", "96/4", scope={"metric": "Z/E", "order": ["Z", "E"], "bond": "not specified in Scheme 2"})
r["notes"].append(direct("Starting material shown as 3bf, Z/E = 95/5.", {"input": "3bf", "metric": "Z/E", "value": "95/5"}))
r["review"]["unresolved"].append("Scheme 2 does not state which product alkene the Z/E ratio refers to; the ratio is retained without assigning a bond.")
classify(r, [atom_change("b1", "C", "C", 0, 1, "formed"), atom_change("b2", "C", "Br", 1, 0, "broken"), atom_change("b3", "C", "B", 1, 0, "broken")], [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}], "The aryl C–C coupling and leaving-group bonds are inferred from the Scheme 2 drawings without a complete atom map.")
reactions.append(r)

# Rh-catalyzed C-H olefination -> 17
r = base_reaction("scheme2_17", "Subsequent C–H functionalization of 3bf with methyl acrylate", "directed indole C–H olefination", [f"{PREFIX}_mol_3bf", mid("methyl acrylate")], "17", solvent=mid("DMF"), temp=80, time_h=12)
add_charge(r["participants"][1], "equivalents", "1.5 equiv", 1.5, "equiv")
r["participants"][1]["equivalence_reference"] = direct("relative to 3bf", f"{PREFIX}_mol_3bf")
p = participant("p_rhodium", mid("[Cp*RhCl2]2"), "catalyst", ref=f"{PREFIX}_mol_3bf"); add_charge(p, "mol_percent", "5 mol%", 5, "mol%"); r["participants"].append(p)
p = participant("p_copper", mid("Cu(OAc)2·H2O"), "reagent", ref=f"{PREFIX}_mol_3bf"); add_charge(p, "equivalents", "1 equiv", 1, "equiv"); r["participants"].append(p)
r["steps"][0]["participant_ids"] += ["p_rhodium", "p_copper"]
r["solvent_conditions"]["total_volume"] = direct("DMF (1 mL)", {"number": 1, "unit": "mL"})
next(p for p in r["participants"] if p["role"] == "solvent")["charges"]["volume"] = direct("DMF (1 mL)", {"number": 1, "unit": "mL"})
outcome(r, "o_yield", "yield", mid("17"), "17, 65%", {"number": 65, "unit": "%"}, yield_type="reported yield")
outcome(r, "o_selectivity", "selectivity", mid("17"), "Z/E = 93/7", "93/7", scope={"metric": "Z/E", "order": ["Z", "E"], "bond": "not specified in Scheme 2"})
r["notes"].append(direct("Starting material shown as 3bf, Z/E = 95/5.", {"input": "3bf", "metric": "Z/E", "value": "95/5"}))
r["review"]["unresolved"].append("Scheme 2 does not state which of product 17's alkenes the Z/E ratio refers to; no bond-specific assignment was invented.")
classify(r, [atom_change("b1", "C", "C", 0, 1, "formed")], [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}], "The newly connected aryl/vinyl C–C bond is visible in the drawings; dehydrogenation details and implicit-H losses are not fully atom mapped.")
reactions.append(r)

# Directing-group removal -> 18
r = base_reaction("scheme2_18", "Removal of the 2-pyrimidinyl directing group from 3bf", "N-deprotection of indole", [f"{PREFIX}_mol_3bf"], "18", solvent=mid("DMSO"), temp="r.t.", time_h=36, atmosphere="Ar")
p = participant("p_naoet", mid("NaOEt"), "reagent", ref=f"{PREFIX}_mol_3bf"); add_charge(p, "equivalents", "4 equiv", 4, "equiv"); r["participants"].append(p); r["steps"][0]["participant_ids"].append("p_naoet")
outcome(r, "o_yield", "yield", mid("18"), "18, 66%", {"number": 66, "unit": "%"}, yield_type="reported yield")
outcome(r, "o_selectivity", "selectivity", mid("18"), "Z/E = 95/5", "95/5", scope={"metric": "Z/E", "order": ["Z", "E"], "bond": "not specified in Scheme 2"})
r["notes"].append(direct("Starting material shown as 3bf, Z/E = 95/5.", {"input": "3bf", "metric": "Z/E", "value": "95/5"}))
classify(r, [atom_change("b1", "N", "H", 0, 1, "formed"), atom_change("b2", "C", "N", 1, 0, "broken")], [{"category_number": 29, "code": "N-H", "supporting_change_ids": ["b1"]}], "The product drawing shows indole N–H after cleavage of the N–pyrimidinyl bond; the hydrogen source and full atom mapping remain unresolved.")
reactions.append(r)

# 3am cycloaddition -> cis-19
r = base_reaction("scheme2_cis19", "Diels–Alder reaction of 3am with diethyl acetylenedicarboxylate", "Lewis-acid-catalyzed Diels–Alder cycloaddition", [f"{PREFIX}_mol_3am", mid("diethyl acetylenedicarboxylate")], "cis-19", temp=100, time_h=24, atmosphere="Ar")
add_charge(r["participants"][1], "equivalents", "3 equiv", 3, "equiv")
r["participants"][1]["equivalence_reference"] = direct("relative to 3am", f"{PREFIX}_mol_3am")
p = participant("p_scandium", mid("Sc(OTf)3"), "catalyst", ref=f"{PREFIX}_mol_3am"); add_charge(p, "mol_percent", "10 mol%", 10, "mol%"); r["participants"].append(p); r["steps"][0]["participant_ids"].append("p_scandium")
outcome(r, "o_yield", "yield", mid("cis-19"), "cis-19, 55%", {"number": 55, "unit": "%"}, yield_type="reported yield")
outcome(r, "o_relative_configuration", "selectivity", mid("cis-19"), "cis-19", {"configuration": "cis", "ratio": None}, scope={"metric": "relative configuration", "order": None, "assignment": "cis"})
r["notes"].append(direct("Starting material 3am is shown with Z/E of C1=C2 = 99/1 and E/Z of C3=C4 = 100/0.", {"input": "3am", "C1=C2": {"metric": "Z/E", "value": "99/1"}, "C3=C4": {"metric": "E/Z", "value": "100/0"}}))
r["review"]["unresolved"].append("No solvent is displayed for the cis-19 arrow; it remains unreported rather than inferred from another Scheme 2 reaction.")
classify(r, [atom_change("b1", "C", "C", 0, 1, "formed"), atom_change("b2", "C", "C", 0, 1, "formed"), atom_change("b3", "C", "C", 1, 2, "order_increased"), atom_change("b4", "C", "C", 2, 1, "order_decreased"), atom_change("b5", "C", "C", 2, 1, "order_decreased"), atom_change("b6", "C", "C", 3, 2, "order_decreased")], [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}, {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b3"]}], "Cycloaddition changes are inferred from the drawings; cis relative configuration is retained separately and atom correspondence is partial.")
reactions.append(r)

# E-3an directing-group removal -> E-20
r = base_reaction("scheme2_E20", "Removal of the 2-pyrimidinyl directing group from E-3an", "N-deprotection of indole", [mid("E-3an")], "E-20", solvent=mid("DMSO"), temp=40, time_h=30, atmosphere="Ar")
p = participant("p_naoet", mid("NaOEt"), "reagent", ref=mid("E-3an")); add_charge(p, "equivalents", "4 equiv", 4, "equiv"); r["participants"].append(p); r["steps"][0]["participant_ids"].append("p_naoet")
outcome(r, "o_yield", "yield", mid("E-20"), "E-20, 58%", {"number": 58, "unit": "%"}, yield_type="reported yield")
outcome(r, "o_configuration_retention", "selectivity", mid("E-20"), "E-20 could be produced without isomerization", {"configuration": "E", "isomerization": "not observed", "ratio": None}, scope={"metric": "geometric configuration retention", "order": None, "assignment": "E"})
r["review"]["unresolved"].append("The phrase 'without isomerization' does not provide a numerical E/Z ratio or analytical method; neither was inferred.")
classify(r, [atom_change("b1", "N", "H", 0, 1, "formed"), atom_change("b2", "C", "N", 1, 0, "broken")], [{"category_number": 29, "code": "N-H", "supporting_change_ids": ["b1"]}], "The product drawing shows indole N–H after directing-group cleavage; the hydrogen source and complete atom mapping remain unresolved.")
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": rid("page7_nonexperimental_sections"),
        "document_id": DOC,
        "kind": "article_context",
        "evidence_ids": evid(7),
        "note": "Conclusions, data availability, author contributions, conflicts, acknowledgements, and references are non-synthetic context. The page states that ESI exists, but no SI file was explicitly associated in the workspace."
    },
    {
        "context_id": rid("page8_references"),
        "document_id": DOC,
        "kind": "bibliography",
        "evidence_ids": evid(8),
        "note": "Page 8 contains bibliography entries only. Cited literature was not converted into experiments of this paper."
    }
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and task["evidence_id"] in {f"{DOC}:page7:text", f"{DOC}:page8:text"}:
        task["state"] = "attempted"
        task["reason"] = "Page image visually inspected; Scheme 2 arrows or bibliography-only status represented in coverage with unresolved items explicit."
(RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "source_context": len(bundle["source_context"])}, ensure_ascii=False))
