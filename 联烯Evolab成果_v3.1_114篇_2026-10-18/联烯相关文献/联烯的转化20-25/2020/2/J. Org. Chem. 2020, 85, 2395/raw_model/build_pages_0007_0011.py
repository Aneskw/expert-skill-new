#!/usr/bin/env python3
"""Extend the v3.1 extraction with source facts visible on PDF pages 7--11.

The existing standard bundle is used only to de-duplicate and augment experiments
whose table rows appeared on earlier pages.  Every added fact below is sourced to
pages 7--11; the two 1 mmol procedures are new, independent experiments.
"""

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_fdf0107b6851842e"
DATASET = json.loads((RUN / "dataset.json").read_text())
OLD = json.loads((RUN / "standard_records.json").read_text())
E = {p: f"{DOC}:page{p}:text" for p in range(7, 12)}
I = {p: f"{DOC}:page{p}:image" for p in range(7, 12)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def inferred(normalized, evidence, why):
    return f(None, normalized, evidence, method="inference", derivation=why)


def derived(normalized, evidence, why):
    return f(None, normalized, evidence, method="structure_derivation", derivation=why)


def uniq_extend(seq, values):
    for value in values:
        if value not in seq:
            seq.append(value)


def add_record_evidence(record, evidence):
    uniq_extend(record["evidence_ids"], evidence)


def participant(pid, mid, role, evidence, stage="reaction", amount=None,
                mass=None, volume=None, equiv=None, molpct=None, ref=None,
                product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if amount is not None:
        p["charges"]["amount"] = f(str(amount), amount, evidence)
    if mass is not None:
        p["charges"]["mass"] = f(str(mass), mass, evidence)
    if volume is not None:
        p["charges"]["volume"] = f(str(volume), volume, evidence)
    if equiv is not None:
        p["charges"]["equivalents"] = f(str(equiv), equiv, evidence)
    if molpct is not None:
        p["charges"]["mol_percent"] = f(str(molpct), molpct, evidence)
    if ref is not None:
        p["equivalence_reference"] = f(ref, ref, evidence)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def add_step(record, sid, evidence, pids=(), operation=None, temperature=None,
             time=None, atmosphere=None, stage="reaction"):
    step = st.new_step(sid, evidence, stage=stage)
    step["participant_ids"] = list(pids)
    if operation is not None:
        step["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        step["temperature"] = f(str(temperature), temperature, evidence)
    if time is not None:
        step["time"] = f(str(time), time, evidence)
    if atmosphere is not None:
        step["atmosphere"] = f(atmosphere, atmosphere, evidence)
    record["steps"].append(step)
    return step


def add_outcome(record, oid, kind, targets, evidence, value, raw=None,
                basis=None, role="main product", scope=None):
    # Re-running the extension against its already-merged bundle must replace,
    # not duplicate, the same source-derived outcome.
    record["outcomes"] = [o for o in record["outcomes"] if o["outcome_id"] != oid]
    outcome = st.new_outcome(oid, kind, targets, evidence)
    outcome["value"] = f(raw if raw is not None else str(value), value, evidence)
    outcome["product_role"] = f(role, role, evidence)
    if kind == "yield":
        outcome["yield_type"] = f("isolated yield", "isolated", evidence)
    if basis is not None:
        outcome["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        outcome["scope"] = f(scope, scope, evidence)
    record["outcomes"].append(outcome)
    return outcome


def add_participant(record, item):
    record["participants"] = [
        p for p in record["participants"] if p["participant_id"] != item["participant_id"]
    ]
    record["participants"].append(item)


def base_reaction(rid, evidence, category, location, name, transformation):
    record = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    record["record_type"] = "specific_experiment"
    record["experiment_categories"] = f([category], [category], evidence)
    record["location"] = f(location, location, evidence)
    record["reaction_name"] = f(name, name, evidence)
    record["transformation_type"] = inferred(
        transformation, evidence,
        "Transformation description inferred from the source procedure and named products.")
    record["review"]["unresolved"] = []
    return record


CHANGE_CC = {
    "change_id": "b1",
    "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
    "before_order": 0,
    "after_order": 1,
    "change_type": "formed",
}
CHANGE_OH = {
    "change_id": "b1",
    "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}],
    "before_order": 0,
    "after_order": 1,
    "change_type": "hydrogen_gain",
    "hydrogen_count": 1,
}
CHANGE_DBCC = {
    "change_id": "b2",
    "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
    "before_order": 1,
    "after_order": 2,
    "change_type": "order_increased",
}


def set_classification(record, evidence, changes, labels, limitation):
    record["classification"]["mapped_reaction_smiles"] = f()
    record["classification"]["bond_changes"] = [
        derived(change, evidence,
                "Net change inferred from named/drawn reactant and product; no complete atom mapping asserted.")
        for change in changes
    ]
    record["classification"]["labels"] = derived(
        labels, evidence, "Taxonomy labels follow only the listed net formed/promoted bonds.")
    record["classification"]["coverage"] = inferred(
        "partial", evidence,
        "No complete source-verified atom mapping is available; classification is deliberately partial.")
    record["classification"]["limitations"] = [limitation]
    record["classification"]["model_votes"] = []


molecules = {m["molecule_id"]: copy.deepcopy(m) for m in OLD["molecules"]}
reactions = {r["reaction_id"]: copy.deepcopy(r) for r in OLD["reactions"]}
incoming_reaction_ids = []


def reaction(short_id):
    rid = f"{DOC}:{short_id}"
    if rid not in reactions:
        raise KeyError(rid)
    if rid not in incoming_reaction_ids:
        incoming_reaction_ids.append(rid)
    return reactions[rid]


def extend_molecule(mid, name, evidence, geometry=None, representation=None,
                    double_bonds=None):
    mol = molecules[mid]
    mol["names"]["original"] = f(name, name, evidence)
    mol["names"]["english"] = f(name, name, evidence)
    occurrence = mol["occurrences"][0]
    uniq_extend(occurrence["evidence_ids"], evidence)
    uniq_extend(occurrence["compound_label"]["evidence_ids"], evidence)
    if geometry is not None:
        mol["stereochemistry"]["geometry_status"] = inferred(
            geometry, evidence,
            "Geometry status follows the explicit stereochemical prefix in the source name.")
    if representation is not None:
        mol["stereochemistry"]["representation_status"] = inferred(
            representation, evidence,
            "The source assignment is retained textually; the current machine representation does not fully encode it.")
    if double_bonds is not None:
        mol["stereochemistry"]["double_bonds"] = f(
            double_bonds, double_bonds, evidence)


names_2 = {
    "da": "Ethyl 5-(5-(4-bromophenyl)-2-methylfuran-3-yl)-2-diazo-3-oxo-5-phenylpentanoate",
    "ea": "Ethyl 5-(5-(4-chlorophenyl)-2-methylfuran-3-yl)-2-diazo-3-oxo-5-phenylpentanoate",
    "fa": "Ethyl 2-diazo-5-(5-(4-fluorophenyl)-2-methylfuran-3-yl)-3-oxo-5-phenylpentanoate",
    "ga": "Ethyl 2-diazo-5-(2-methyl-5-(thiophen-3-yl)furan-3-yl)-3-oxo-5-phenylpentanoate",
    "ha": "Ethyl 5-(5-butyl-2-methylfuran-3-yl)-2-diazo-3-oxo-5-phenylpentanoate",
    "ia": "Ethyl 5-(5-(tert-butyl)-2-methylfuran-3-yl)-2-diazo-3-oxo-5-phenylpentanoate",
    "ja": "Ethyl 2-diazo-5-(2-methylfuran-3-yl)-3-oxo-5-phenylpentanoate",
    "ab": "2-Diazo-5-(2-methyl-5-phenylfuran-3-yl)-1,5-diphenylpentane-1,3-dione",
    "bb": "2-Diazo-5-(2-methyl-5-(p-tolyl)furan-3-yl)-1,5-diphenylpentane-1,3-dione",
    "cb": "2-Diazo-5-(5-(4-fluorophenyl)-2-methylfuran-3-yl)-1,5-diphenylpentane-1,3-dione",
    "db": "2-Diazo-5-(2-methyl-5-(thiophen-3-yl)furan-3-yl)-1,5-diphenylpentane-1,3-dione",
    "eb": "5-(5-Butyl-2-methylfuran-3-yl)-2-diazo-1,5-diphenylpentane-1,3-dione",
    "fb": "5-(5-(tert-butyl)-2-methylfuran-3-yl)-2-diazo-1,5-diphenylpentane-1,3-dione",
    "gb": "2-Diazo-5-(2-methylfuran-3-yl)-1,5-diphenylpentane-1,3-dione",
}
page_2 = {code: 7 for code in ["da", "ea", "fa", "ga", "ha", "ia", "ja", "ab", "bb", "cb", "db", "eb"]}
page_2.update({"fb": 8, "gb": 8})
for code, name in names_2.items():
    page = page_2[code]
    extend_molecule(f"{DOC}:mol:2{code}", name, [E[page], I[page]])


names_3 = {
    "aa": "(Z)-Ethyl 5-methyl-3-oxo-6-(2-oxo-2-phenylethylidene)-1,2,3,6-tetrahydro-[1,1'-biphenyl]-4-carboxylate",
    "fa": "(Z)-Ethyl 6-(2-(4-fluorophenyl)-2-oxoethylidene)-5-methyl-3-oxo-1,2,3,6-tetrahydro-[1,1'-biphenyl]-4-carboxylate",
    "ha": "(Z)-Ethyl 5-methyl-3-oxo-6-(2-oxohexylidene)-1,2,3,6-tetrahydro-[1,1'-biphenyl]-4-carboxylate",
    "ia": "(Z)-Ethyl 6-(3,3-dimethyl-2-oxobutylidene)-5-methyl-3-oxo-1,2,3,6-tetrahydro-[1,1'-biphenyl]-4-carboxylate",
    "ja": "(Z)-Ethyl 5-methyl-3-oxo-6-(2-oxoethylidene)-1,2,3,6-tetrahydro-[1,1'-biphenyl]-4-carboxylate",
    "ab": "(Z)-4-Benzoyl-5-methyl-6-(2-oxo-2-phenylethylidene)-1,6-dihydro-[1,1'-biphenyl]-3(2H)-one",
    "bb": "(Z)-4-Benzoyl-5-methyl-6-(2-oxo-2-(p-tolyl)ethylidene)-1,6-dihydro-[1,1'-biphenyl]-3(2H)-one",
    "cb": "(Z)-4-Benzoyl-4'-fluoro-5-methyl-6-(2-oxo-2-phenylethylidene)-1,6-dihydro-[1,1'-biphenyl]-3(2H)-one",
    "eb": "(Z)-4-Benzoyl-5-methyl-6-(2-oxohexylidene)-1,6-dihydro-[1,1'-biphenyl]-3(2H)-one",
    "gb": "(Z)-2-(4-benzoyl-3-methyl-5-oxo-5,6-dihydro-[1,1'-biphenyl]-2(1H)-ylidene)acetaldehyde",
}
page_3 = {code: 8 for code in ["aa", "fa", "ha", "ia", "ja", "ab", "bb", "cb"]}
page_3.update({"eb": 9, "gb": 9})
for code, name in names_3.items():
    page = page_3[code]
    extend_molecule(
        f"{DOC}:mol:3{code}", name, [E[page], I[page]], geometry="specified",
        representation="partial",
        double_bonds=[{"configuration": "Z", "basis": "explicit (Z) prefix in source name"}],
    )


names_4 = {
    "aa": "Ethyl 5-hydroxy-3-methyl-2-(2-oxo-2-phenylethyl)-[1,1'-biphenyl]-4-carboxylate",
    "ba": "Ethyl 5-hydroxy-3-methyl-2-(2-oxo-2-(p-tolyl)ethyl)-[1,1'-biphenyl]-4-carboxylate",
    "ca": "Ethyl 5-hydroxy-2-(2-(4-methoxyphenyl)-2-oxoethyl)-3-methyl-[1,1'-biphenyl]-4-carboxylate",
    "da": "Ethyl 2-(2-(4-bromophenyl)-2-oxoethyl)-5-hydroxy-3-methyl-[1,1'-biphenyl]-4-carboxylate",
    "ea": "Ethyl 2-(2-(4-chlorophenyl)-2-oxoethyl)-5-hydroxy-3-methyl-[1,1'-biphenyl]-4-carboxylate",
    "fa": "Ethyl 2-(2-(4-fluorophenyl)-2-oxoethyl)-5-hydroxy-3-methyl-[1,1'-biphenyl]-4-carboxylate",
    "ga": "Ethyl 5-hydroxy-3-methyl-2-(2-oxo-2-(thiophen-3-yl)ethyl)-[1,1'-biphenyl]-4-carboxylate",
    "ha": "Ethyl 5-hydroxy-3-methyl-2-(2-oxohexyl)-[1,1'-biphenyl]-4-carboxylate",
    "ia": "Ethyl 2-(3,3-dimethyl-2-oxobutyl)-5-hydroxy-3-methyl-[1,1'-biphenyl]-4-carboxylate",
    "ab": "2-(4-Benzoyl-5-hydroxy-3-methyl-[1,1'-biphenyl]-2-yl)-1-phenylethanone",
    "cb": "2-(4-Benzoyl-5-hydroxy-3-methyl-[1,1'-biphenyl]-2-yl)-1-(4-fluorophenyl)ethanone",
    "eb": "1-(4-Benzoyl-5-hydroxy-3-methyl-[1,1'-biphenyl]-2-yl)hexan-2-one",
}
page_4 = {code: 9 for code in ["aa", "ba", "ca", "da", "ea", "fa", "ga"]}
page_4.update({code: 10 for code in ["ha", "ia", "ab", "cb", "eb"]})
for code, name in names_4.items():
    page = page_4[code]
    extend_molecule(f"{DOC}:mol:4{code}", name, [E[page], I[page]])


# Sodium sulfate is explicitly introduced only by the page 9 workup.
na2so4_id = f"{DOC}:mol:na2so4"
na2so4 = st.new_molecule(na2so4_id, DOC, [E[9], I[9]])
na2so4["identity_status"] = "source_identified"
na2so4["occurrences"][0]["compound_label"] = f("Na2SO4", "Na2SO4", [E[9], I[9]])
na2so4["names"]["original"] = f("Na2SO4", "Na2SO4", [E[9], I[9]])
na2so4["names"]["english"] = f("sodium sulfate", "sodium sulfate", [E[9], I[9]])
na2so4["names"]["abbreviations"] = f("Na2SO4", ["Na2SO4"], [E[9], I[9]])
na2so4["names"]["aliases"] = f()
na2so4["cas_numbers"] = f()
na2so4["representation_type"] = "smiles"
na2so4_smiles = "[Na+].[Na+].[O-]S(=O)(=O)[O-]"
na2so4_diag = st.stereo_inspect(na2so4_smiles)
na2so4["smiles"] = derived(
    na2so4_smiles, [E[9], I[9]],
    "Simple named drying agent represented as a machine candidate; source review remains pending.")
na2so4["stereochemistry"]["atom_centers"] = f([], [], [E[9], I[9]])
na2so4["stereochemistry"]["double_bonds"] = f([], [], [E[9], I[9]])
na2so4["stereochemistry"]["axial"] = f([], [], [E[9], I[9]])
na2so4["stereochemistry"]["mixture"] = f(None, None, [E[9], I[9]])
na2so4["stereochemistry"]["assignment"] = f("achiral", "achiral", [E[9], I[9]])
na2so4["stereochemistry"]["chirality_status"] = inferred(
    "achiral", [E[9], I[9]], "No stereogenic feature is present in the candidate graph.")
na2so4["stereochemistry"]["geometry_status"] = inferred(
    "not_applicable", [E[9], I[9]], "No source-relevant E/Z or ring geometry applies.")
na2so4["stereochemistry"]["representation_status"] = inferred(
    "encoded", [E[9], I[9]], "Candidate connectivity is fully represented for this simple material.")
na2so4["stereochemistry"]["isomeric_smiles"] = derived(
    na2so4_diag["canonical_isomeric_smiles"], [E[9], I[9]],
    "RDKit canonical isomeric SMILES; source_verified=false.")
na2so4["stereochemistry"]["identity_key"] = derived(
    na2so4_diag["stereo_graph_key"], [E[9], I[9]],
    "RDKit stereo graph key from the candidate; not a reviewed sample identity.")
na2so4["identity_note"] = f(
    "Explicitly named drying agent in section 4.5; no CAS number inferred.",
    "Explicitly named drying agent in section 4.5; no CAS number inferred.",
    [E[9], I[9]],
)
molecules[na2so4_id] = na2so4


# Characterization-paragraph isolated masses for Table 2 products.
product2_masses = {
    "da": (7, 78.0), "ea": (7, 69.9), "fa": (7, 65.6), "ga": (7, 65.4),
    "ha": (7, 72.7), "ia": (7, 61.2), "ja": (7, 42.4), "ab": (7, 69.5),
    "bb": (7, 71.8), "cb": (7, 67.9), "db": (7, 69.6), "eb": (8, 68.0),
    "fb": (8, 70.5), "gb": (8, 56.6),
}
for code, (page, mass) in product2_masses.items():
    ev = [E[page], I[page]]
    record = reaction(f"table2_scope_2{code}")
    add_record_evidence(record, ev)
    add_outcome(
        record, f"product_amount_2{code}_characterization", "product_amount",
        [f"{DOC}:mol:2{code}"], ev, {"number": mass, "unit": "mg"},
        raw=f"{mass:g} mg",
        basis="isolated material reported in the compound characterization paragraph",
    )


# Section 4.4 gives the detailed common procedure for every Table 3 row.
table3_codes = ["aa", "ba", "ca", "da", "ea", "fa", "ga", "ha", "ia", "ja",
                "ab", "bb", "cb", "db", "eb", "fb", "gb"]
for code in table3_codes:
    ev = [E[8], I[8]]
    record = reaction(f"table3_scope_2{code}")
    add_record_evidence(record, ev)
    catalyst = next(p for p in record["participants"] if p["participant_id"] == "catalyst")
    catalyst["charges"]["mass"] = f("8.8 mg", {"number": 8.8, "unit": "mg"}, ev)
    catalyst["charges"]["amount"] = f("0.02 mmol", {"number": 0.02, "unit": "mmol"}, ev)
    add_participant(record, participant("dcm_substrate_solution", f"{DOC}:mol:dcm", "solvent", ev,
                                        volume={"number": 1, "unit": "mL"}))
    add_participant(record, participant("dcm_catalyst_solution", f"{DOC}:mol:dcm", "solvent", ev,
                                        volume={"number": 2, "unit": "mL"}))
    add_participant(record, participant("silica", f"{DOC}:mol:silica", "purification_agent", ev,
                                        stage="purification"))
    add_participant(record, participant("petroleum_ether", f"{DOC}:mol:pe", "purification_agent", ev,
                                        stage="purification"))
    add_participant(record, participant("ethyl_acetate", f"{DOC}:mol:ea", "purification_agent", ev,
                                        stage="purification"))
    product_ids = [p["participant_id"] for p in record["participants"] if p["role"] == "product"]
    record["steps"] = []
    add_step(record, "prepare_substrate_solution", ev, ["substrate", "dcm_substrate_solution"],
             "prepare diazo-furan solution in DCM")
    add_step(record, "dropwise_addition_and_reaction", ev,
             ["substrate", "catalyst", "dcm_substrate_solution", "dcm_catalyst_solution"],
             "syringe substrate solution dropwise into vigorously stirred Rh2(OAc)4 solution in a water bath",
             {"number": 25, "unit": "°C"}, {"number": 30, "unit": "min"}, "nitrogen")
    add_step(record, "concentrate", ev, [], "concentrate under reduced pressure", stage="workup")
    add_step(record, "column", ev, product_ids + ["silica", "petroleum_ether", "ethyl_acetate"],
             "purify crude product by silica-gel column chromatography with a suitable PE/EA mixture",
             stage="purification")
    record["procedure"]["original_text"] = f(
        "Section 4.4: 2 (0.2 mmol) in DCM (1 mL) was added dropwise to Rh2(OAc)4 (8.8 mg, 0.02 mmol) in anhydrous DCM (2 mL), 25 °C, 30 min; concentrate and purify by silica/PE-EA.",
        "Detailed section 4.4 general procedure", ev)
    record["procedure"]["general_procedure_reference"] = f(
        "Section 4.4 general procedure", "Section 4.4", ev)
    record["procedure"]["addition_order"] = f(
        "substrate solution in DCM into Rh2(OAc)4 solution in anhydrous DCM",
        ["substrate solution", "Rh2(OAc)4 solution"], ev)


# Table 3 product masses; direct aromatized products are also characterized later.
table3_product_masses = {
    "aa": (8, "3", 57.7), "fa": (8, "3", 64.4), "ha": (8, "3", 59.5),
    "ia": (8, "3", 42.5), "ja": (8, "3", 48.9), "ab": (8, "3", 34.1),
    "bb": (8, "3", 54.7), "cb": (9, "3", 67.9), "eb": (9, "3", 58.0),
    "gb": (9, "3", 34.8), "ba": (9, "4", 52.8), "ca": (9, "4", 66.3),
    "da": (9, "4", 61.6), "ea": (9, "4", 59.7), "ga": (9, "4", 60.9),
}
for code, (page, series, mass) in table3_product_masses.items():
    ev = [E[page], I[page]]
    record = reaction(f"table3_scope_2{code}")
    add_record_evidence(record, ev)
    add_outcome(
        record, f"product_amount_{series}{code}_characterization", "product_amount",
        [f"{DOC}:mol:{series}{code}"], ev, {"number": mass, "unit": "mg"},
        raw=f"{mass:g} mg",
        basis="isolated material reported in the compound characterization paragraph",
    )


# Page 9 prints 45% for 3gb while Table 3 prints 40%; retain both as a conflict.
gb_record = reaction("table3_scope_2gb")
gb_yield = next(o for o in gb_record["outcomes"] if o["outcome_id"] == "yield")
gb_yield["value"] = f(
    ["40% (Table 3)", "45% (section 4.4.10)"], None,
    [f"{DOC}:page3:text", f"{DOC}:page3:image", E[9], I[9]], status="conflict",
    alternatives=[
        f("40%", {"number": 40, "unit": "%"}, [f"{DOC}:page3:text", f"{DOC}:page3:image"]),
        f("45%", {"number": 45, "unit": "%"}, [E[9], I[9]]),
    ],
)
gb_note = "3gb printed yield conflict retained without selecting a value"
gb_record["notes"] = [
    note for note in gb_record["notes"] if note.get("normalized_value") != gb_note
]
gb_record["notes"].append(f(
    "Table 3 reports 40% for 3gb, whereas section 4.4.10 reports 34.8 mg, 45%; both printed yields are retained.",
    gb_note, [f"{DOC}:page3:text", f"{DOC}:page3:image", E[9], I[9]],
))
gb_unresolved = "3gb yield conflicts between Table 3 (40%) and section 4.4.10 (45%); no value selected."
if gb_unresolved not in gb_record["review"]["unresolved"]:
    gb_record["review"]["unresolved"].append(gb_unresolved)


# Section 4.5 gives common workup details for all ten Scheme 2 arrows.
scheme2_codes = ["aa", "fa", "ha", "ia", "ja", "ab", "bb", "cb", "eb", "gb"]
for code in scheme2_codes:
    ev = [E[9], I[9]]
    record = reaction(f"scheme2_isomerization_3{code}")
    add_record_evidence(record, ev)
    record["steps"] = [step for step in record["steps"] if step["step_id"] == "reaction"]
    bf3 = next(p for p in record["participants"] if p["participant_id"] == "bf3")
    bf3["charges"]["amount"] = f("0.2 mmol", {"number": 0.2, "unit": "mmol"}, ev)
    add_participant(record, participant("nahco3_quench", f"{DOC}:mol:nahco3", "workup_agent", ev,
                                        stage="workup"))
    add_participant(record, participant("dcm_extract", f"{DOC}:mol:dcm", "solvent", ev,
                                        stage="workup"))
    add_participant(record, participant("na2so4_dry", na2so4_id, "workup_agent", ev,
                                        stage="workup"))
    add_participant(record, participant("silica", f"{DOC}:mol:silica", "purification_agent", ev,
                                        stage="purification"))
    add_participant(record, participant("petroleum_ether", f"{DOC}:mol:pe", "purification_agent", ev,
                                        stage="purification"))
    add_participant(record, participant("ethyl_acetate", f"{DOC}:mol:ea", "purification_agent", ev,
                                        stage="purification"))
    add_step(record, "quench", ev, ["nahco3_quench"],
             "quench with saturated aqueous NaHCO3", stage="workup")
    add_step(record, "extract", ev, ["dcm_extract"],
             "extract the water layer with DCM and combine organic layers", stage="workup")
    add_step(record, "dry_and_concentrate", ev, ["na2so4_dry"],
             "dry combined organic layers over Na2SO4 and concentrate in vacuum", stage="workup")
    add_step(record, "column", ev, ["product", "silica", "petroleum_ether", "ethyl_acetate"],
             "purify residue by silica-gel column chromatography with a suitable PE/EA mixture",
             stage="purification")
    record["procedure"]["original_text"] = f(
        "Section 4.5: 3 (0.2 mmol) in DCM (2 mL) was added to BF3.OEt2 (0.2 mmol, 1.0 equiv) at 25 °C for several hours, quenched with saturated aqueous NaHCO3, extracted with DCM, dried over Na2SO4, concentrated, and chromatographed.",
        "Detailed section 4.5 general procedure", ev)
    record["procedure"]["general_procedure_reference"] = f(
        "Section 4.5 general procedure", "Section 4.5", ev)
    record["procedure"]["addition_order"] = f(
        "solution of compound 3 in DCM into BF3.OEt2",
        ["compound 3 solution", "BF3.OEt2"], ev)


scheme2_product_masses = {
    "aa": (9, 57.7), "fa": (9, 64.4), "ha": (10, 59.5), "ia": (10, 42.5),
    "ab": (10, 34.1), "cb": (10, 67.9), "eb": (10, 34.8),
}
for code, (page, mass) in scheme2_product_masses.items():
    ev = [E[page], I[page]]
    record = reaction(f"scheme2_isomerization_3{code}")
    add_record_evidence(record, ev)
    add_outcome(
        record, f"product_amount_4{code}_characterization", "product_amount",
        [f"{DOC}:mol:4{code}"], ev, {"number": mass, "unit": "mg"},
        raw=f"{mass:g} mg",
        basis="isolated material reported in the compound characterization paragraph",
    )


# Two independent 1 mmol experiments printed in section 4.5.1.
ev = [E[9], I[9]]
scale_rh = base_reaction(
    "section4_5_1_scaleup_2aa_to_3aa", ev, "scale_up",
    {"page": 9, "region": "Section 4.5.1 first 1 mmol procedure"},
    "1 mmol scale synthesis of 3aa", "Rh2(OAc)4-mediated cascade cyclization",
)
scale_rh["general_scheme_ids"] = [f"{DOC}:general_rh_cascade"]
scale_rh["participants"] = [
    participant("substrate", f"{DOC}:mol:2aa", "reactant", ev,
                amount={"number": 1, "unit": "mmol"}, mass={"number": 402.4, "unit": "mg"}),
    participant("product", f"{DOC}:mol:3aa", "product", ev, product_role="main product"),
    participant("catalyst", f"{DOC}:mol:rh2_oac4", "catalyst", ev,
                amount={"number": 0.1, "unit": "mmol"}, mass={"number": 44.2, "unit": "mg"},
                molpct={"number": 10, "unit": "mol%"}, ref=f"{DOC}:mol:2aa"),
    participant("dcm_substrate_solution", f"{DOC}:mol:dcm", "solvent", ev,
                volume={"number": 4, "unit": "mL"}),
    participant("dcm_catalyst_solution", f"{DOC}:mol:dcm", "solvent", ev,
                volume={"number": 2, "unit": "mL"}),
    participant("silica", f"{DOC}:mol:silica", "purification_agent", ev, stage="purification"),
    participant("petroleum_ether", f"{DOC}:mol:pe", "purification_agent", ev, stage="purification"),
    participant("ethyl_acetate", f"{DOC}:mol:ea", "purification_agent", ev, stage="purification"),
]
scale_rh["reaction_scale"] = f(
    "2aa (402.4 mg, 1 mmol)",
    {"reference_molecule_id": f"{DOC}:mol:2aa", "amount": {"number": 1, "unit": "mmol"}}, ev)
scale_rh["solvent_conditions"]["composition"] = f("anhydrous DCM", ["anhydrous DCM"], ev)
scale_rh["solvent_conditions"]["total_volume"] = f(
    "4 mL + 2 mL", {"number": 6, "unit": "mL"}, ev)
add_step(scale_rh, "prepare_substrate_solution", ev, ["substrate", "dcm_substrate_solution"],
         "prepare 2aa solution in DCM")
add_step(scale_rh, "dropwise_addition_and_reaction", ev,
         ["substrate", "catalyst", "dcm_substrate_solution", "dcm_catalyst_solution"],
         "syringe 2aa solution dropwise into vigorously stirred Rh2(OAc)4 solution in a water bath",
         {"number": 25, "unit": "°C"}, {"number": 30, "unit": "min"}, "nitrogen")
add_step(scale_rh, "concentrate", ev, [], "concentrate under reduced pressure", stage="workup")
add_step(scale_rh, "column", ev, ["product", "silica", "petroleum_ether", "ethyl_acetate"],
         "silica-gel column chromatography, ethyl acetate/petroleum ether = 1:6",
         stage="purification")
scale_rh["procedure"]["original_text"] = f(
    "2aa (402.4 mg, 1 mmol) in DCM (4 mL) was added dropwise to Rh2(OAc)4 (44.2 mg, 0.1 mmol) in anhydrous DCM (2 mL) at 25 °C for 30 min; concentrated and chromatographed (EA/PE 1:6).",
    "Section 4.5.1 first scale-up procedure", ev)
scale_rh["procedure"]["addition_order"] = f(
    "2aa solution into Rh2(OAc)4 solution", ["2aa solution", "Rh2(OAc)4 solution"], ev)
add_outcome(scale_rh, "yield_3aa", "yield", [f"{DOC}:mol:3aa"], ev,
            {"number": 75, "unit": "%"}, raw="75%", basis="isolated after column chromatography")
add_outcome(scale_rh, "product_amount_3aa_mmol", "product_amount", [f"{DOC}:mol:3aa"], ev,
            {"number": 0.75, "unit": "mmol"}, raw="0.75 mmol", basis="isolated material")
add_outcome(scale_rh, "product_amount_3aa_mass", "product_amount", [f"{DOC}:mol:3aa"], ev,
            {"number": 280.8, "unit": "mg"}, raw="280.8 mg", basis="isolated material")
set_classification(
    scale_rh, ev, [CHANGE_CC],
    [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]}],
    "Net ring-forming C-C change only; N2 loss and all bond-order changes are not completely mapped.")
reactions[scale_rh["reaction_id"]] = scale_rh
incoming_reaction_ids.append(scale_rh["reaction_id"])


scale_bf3 = base_reaction(
    "section4_5_1_scaleup_3aa_to_4aa", ev, "scale_up",
    {"page": 9, "region": "Section 4.5.1 second 1 mmol procedure"},
    "1 mmol scale transformation of 3aa into 4aa", "BF3-promoted isomerization/aromatization",
)
scale_bf3["participants"] = [
    participant("substrate", f"{DOC}:mol:3aa", "reactant", ev,
                amount={"number": 1, "unit": "mmol"}, mass={"number": 374.2, "unit": "mg"}),
    participant("product", f"{DOC}:mol:4aa", "product", ev, product_role="main product"),
    participant("bf3", f"{DOC}:mol:bf3_oet2", "reagent", ev,
                amount={"number": 1, "unit": "mmol"}, equiv={"number": 1, "unit": "equiv"},
                ref=f"{DOC}:mol:3aa"),
    participant("dcm_reaction", f"{DOC}:mol:dcm", "solvent", ev,
                volume={"number": 5, "unit": "mL"}),
    participant("nahco3_quench", f"{DOC}:mol:nahco3", "workup_agent", ev, stage="workup"),
    participant("dcm_extract", f"{DOC}:mol:dcm", "solvent", ev, stage="workup"),
    participant("na2so4_dry", na2so4_id, "workup_agent", ev, stage="workup"),
    participant("silica", f"{DOC}:mol:silica", "purification_agent", ev, stage="purification"),
    participant("petroleum_ether", f"{DOC}:mol:pe", "purification_agent", ev, stage="purification"),
    participant("ethyl_acetate", f"{DOC}:mol:ea", "purification_agent", ev, stage="purification"),
]
scale_bf3["reaction_scale"] = f(
    "3aa (374.2 mg, 1 mmol)",
    {"reference_molecule_id": f"{DOC}:mol:3aa", "amount": {"number": 1, "unit": "mmol"}}, ev)
scale_bf3["solvent_conditions"]["composition"] = f("DCM", ["DCM"], ev)
scale_bf3["solvent_conditions"]["total_volume"] = f(
    "5 mL", {"number": 5, "unit": "mL"}, ev)
add_step(scale_bf3, "addition_and_reaction", ev, ["substrate", "bf3", "dcm_reaction"],
         "syringe 3aa solution into BF3.OEt2 in a water bath",
         {"number": 25, "unit": "°C"}, "several hours")
add_step(scale_bf3, "quench", ev, ["nahco3_quench"],
         "quench with saturated aqueous NaHCO3", stage="workup")
add_step(scale_bf3, "extract", ev, ["dcm_extract"],
         "extract the water layer with DCM and combine organic layers", stage="workup")
add_step(scale_bf3, "dry_and_concentrate", ev, ["na2so4_dry"],
         "dry over Na2SO4 and concentrate in vacuum", stage="workup")
add_step(scale_bf3, "column", ev, ["product", "silica", "petroleum_ether", "ethyl_acetate"],
         "silica-gel column chromatography with a suitable PE/EA mixture", stage="purification")
scale_bf3["procedure"]["original_text"] = f(
    "3aa (374.2 mg, 1 mmol) in DCM (5 mL) was added to BF3.OEt2 (1 mmol, 1.0 equiv) at 25 °C for several hours; NaHCO3 quench, DCM extraction, Na2SO4 drying, concentration and column chromatography afforded 4aa.",
    "Section 4.5.1 second scale-up procedure", ev)
scale_bf3["procedure"]["addition_order"] = f(
    "3aa solution into BF3.OEt2", ["3aa solution", "BF3.OEt2"], ev)
add_outcome(scale_bf3, "yield_4aa", "yield", [f"{DOC}:mol:4aa"], ev,
            {"number": 60, "unit": "%"}, raw="60%", basis="isolated after column chromatography")
add_outcome(scale_bf3, "product_amount_4aa_mmol", "product_amount", [f"{DOC}:mol:4aa"], ev,
            {"number": 0.60, "unit": "mmol"}, raw="0.60 mmol", basis="isolated material")
add_outcome(scale_bf3, "product_amount_4aa_mass", "product_amount", [f"{DOC}:mol:4aa"], ev,
            {"number": 224.5, "unit": "mg"}, raw="224.5 mg", basis="isolated material")
set_classification(
    scale_bf3, ev, [CHANGE_OH, CHANGE_DBCC],
    [
        {"category_number": 24, "code": "O-H", "supporting_change_ids": ["b1"]},
        {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b2"]},
    ],
    "Hydrogen source/acceptor and complete aromatic bond-order changes are unresolved.")
reactions[scale_bf3["reaction_id"]] = scale_bf3
incoming_reaction_ids.append(scale_bf3["reaction_id"])


# Context blocks are kept separate from synthesis records.
source_context = [
    {
        "context_id": f"{DOC}:context:pages7_10_characterization",
        "document_id": DOC,
        "evidence_ids": [E[7], I[7], E[8], I[8], E[9], I[9], E[10], I[10]],
        "kind": "characterization",
        "note": "NMR/IR/HRMS, melting-point and Rf blocks support compound identity; spectral peaks and molecular-ion values are not synthesis experiments.",
    },
    {
        "context_id": f"{DOC}:context:page10_associated_content",
        "document_id": DOC,
        "evidence_ids": [E[10], I[10]],
        "kind": "associated_content_notice",
        "note": "The paper states that SI/CIF are available online, but no SI file was supplied in this workspace and none is associated.",
    },
    {
        "context_id": f"{DOC}:context:pages10_11_bibliography",
        "document_id": DOC,
        "evidence_ids": [E[10], I[10], E[11], I[11]],
        "kind": "bibliography",
        "note": "Author information, acknowledgments and cited references are nonexperimental context; cited literature is not represented as this paper's experiments.",
    },
]


# Page-specific coverage uses only records actually evidenced on that page.
page_records = {
    7: [f"{DOC}:table2_scope_2{c}" for c in
        ["da", "ea", "fa", "ga", "ha", "ia", "ja", "ab", "bb", "cb", "db", "eb"]],
    8: ([f"{DOC}:table2_scope_2{c}" for c in ["eb", "fb", "gb"]] +
        [f"{DOC}:table3_scope_2{c}" for c in table3_codes]),
    9: ([f"{DOC}:table3_scope_2{c}" for c in ["cb", "eb", "gb", "ba", "ca", "da", "ea", "ga"]] +
        [f"{DOC}:scheme2_isomerization_3{c}" for c in scheme2_codes] +
        [scale_rh["reaction_id"], scale_bf3["reaction_id"]]),
    10: ([f"{DOC}:table3_scope_2ga"] +
         [f"{DOC}:scheme2_isomerization_3{c}" for c in ["ha", "ia", "ab", "cb", "eb"]]),
    11: [],
}
regions = {
    7: [
        {"region": "Sections 4.3.4-4.3.15", "finding": "Compound-isolation masses/yields and characterization for 2da-2eb; these augment the corresponding Table 2 experiments rather than creating duplicate reactions.", "evidence_ids": [E[7], I[7]]},
        {"region": "NMR/IR/HRMS/mp/Rf blocks", "finding": "Characterization-only measurements retained in source_context; mass spectra and molecular-ion values are not synthesis experiments.", "evidence_ids": [E[7], I[7]]},
    ],
    8: [
        {"region": "Sections 4.3.15-4.3.17", "finding": "Completion of 2eb and isolation/characterization of 2fb and 2gb; linked to existing Table 2 rows.", "evidence_ids": [E[8], I[8]]},
        {"region": "Section 4.4 general procedure", "finding": "Detailed 0.2 mmol Rh2(OAc)4 procedure, including split DCM volumes, 25 °C/30 min, concentration and chromatography; applied to the 17 existing Table 3 rows.", "evidence_ids": [E[8], I[8]]},
        {"region": "Sections 4.4.1-4.4.8", "finding": "Isolation masses and explicit (Z) names for 3aa-3cb; characterization is not duplicated as separate experiments.", "evidence_ids": [E[8], I[8]]},
    ],
    9: [
        {"region": "Sections 4.4.8-4.4.10", "finding": "Isolation/characterization of 3cb, 3eb and 3gb; 3gb yield conflict (40% in Table 3 versus 45% here) retained unresolved.", "evidence_ids": [E[9], I[9]]},
        {"region": "Section 4.5 general procedure", "finding": "Detailed BF3.OEt2 addition, 25 °C reaction, NaHCO3 quench, DCM extraction, Na2SO4 drying and chromatography; applied to the ten existing Scheme 2 arrows.", "evidence_ids": [E[9], I[9]]},
        {"region": "Section 4.5.1 first procedure", "finding": "Independent 1 mmol scale 2aa-to-3aa experiment: 75% isolated, 0.75 mmol/280.8 mg.", "evidence_ids": [E[9], I[9]]},
        {"region": "Section 4.5.1 second procedure", "finding": "Independent 1 mmol scale 3aa-to-4aa experiment: several-hour duration, 60% isolated, 0.60 mmol/224.5 mg.", "evidence_ids": [E[9], I[9]]},
        {"region": "Sections 4.5.2-4.5.8", "finding": "Isolation masses and characterization for 4aa-4ga; linked to existing Table 3 or Scheme 2 experiments as appropriate.", "evidence_ids": [E[9], I[9]]},
    ],
    10: [
        {"region": "Sections 4.5.8-4.5.13", "finding": "Completion of 4ga and isolation/characterization of 4ha, 4ia, 4ab, 4cb and 4eb; linked to existing reactions.", "evidence_ids": [E[10], I[10]]},
        {"region": "Associated Content and Author Information", "finding": "SI/CIF availability notice and bibliographic context only; no SI file is supplied or associated in this run.", "evidence_ids": [E[10], I[10]]},
        {"region": "References 1-7 start", "finding": "Cited literature only, not experiments performed in this paper.", "evidence_ids": [E[10], I[10]]},
    ],
    11: [
        {"region": "References 7-22", "finding": "Bibliography only; no reaction arrow, table row, footnote, or experimental procedure from this paper is present.", "evidence_ids": [E[11], I[11]]},
    ],
}
unresolved = {
    7: ["Spectral assignments were not independently interpreted; characterization remains source context.", "No SI was supplied or associated."],
    8: ["The general procedure says a suitable PE/EA mixture but does not give a ratio.", "Structures and stereochemical candidates remain machine records requiring review."],
    9: ["3gb yield conflicts between Table 3 (40%) and section 4.4.10 (45%).", "The second 1 mmol procedure reports only 'several hours'; no exact duration is inferred.", "Several characterization mass/yield pairs are retained verbatim and are not silently recalculated.", "No SI was supplied or associated."],
    10: ["Associated-content text refers to files not present in the workspace; no SI relationship is inferred.", "Spectral blocks are not treated as independent synthesis experiments."],
    11: ["Reference-only page; cited methods were not imported as experiments of this paper."],
}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": page, "state": "attempted", "regions": regions[page],
         "record_ids": page_records[page], "unresolved": unresolved[page]}
        for page in range(7, 12)
    ],
}


stereo_checks = {}
page_molecule_ids = (
    [f"{DOC}:mol:2{c}" for c in names_2] +
    [f"{DOC}:mol:3{c}" for c in names_3] +
    [f"{DOC}:mol:4{c}" for c in names_4] + [na2so4_id]
)
for mid in page_molecule_ids:
    mol = molecules[mid]
    smiles = mol["smiles"]["normalized_value"]
    if smiles:
        stereo_checks[mid] = st.stereo_inspect(smiles)
    else:
        stereo_checks[mid] = {
            "canonical_isomeric_smiles": None,
            "source_assignment": mol["stereochemistry"]["double_bonds"]["normalized_value"],
            "source_verified": False,
            "limitations": [
                "Concrete SMILES withheld because the ring/exocyclic geometry transcription is not sufficiently secure.",
                "Explicit (Z) source assignment is preserved textually with partial representation status.",
            ],
        }


dependency_reaction_ids = [f"{DOC}:general_furanization", f"{DOC}:general_rh_cascade"]
bundle = st.empty(DATASET["documents"])
bundle["molecules"] = list(molecules.values())
bundle["reactions"] = [reactions[rid] for rid in dependency_reaction_ids + incoming_reaction_ids]
bundle["source_context"] = source_context

raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visual_pages_opened_original_detail": [7, 8, 9, 10, 11],
    "new_reaction_ids": [scale_rh["reaction_id"], scale_bf3["reaction_id"]],
    "submission_dependency_reaction_ids": dependency_reaction_ids,
    "updated_existing_reaction_ids": [
        rid for rid in incoming_reaction_ids
        if rid not in {scale_rh["reaction_id"], scale_bf3["reaction_id"]}
    ],
    "new_molecule_ids": [na2so4_id],
    "page_regions": regions,
    "decisions": [
        "Product characterization paragraphs augment the already-recorded Table 2/Table 3/Scheme 2 experiments and are not duplicate reactions.",
        "The two 1 mmol procedures in section 4.5.1 are separate scale-up experiments.",
        "NMR/IR/HRMS/mp/Rf blocks and bibliography are source_context, not synthesis experiments.",
        "The 3gb 40%/45% source conflict is preserved with alternatives and no selected normalized value.",
        "All explicit (Z) labels are retained; missing stereochemical assignments are not converted to achiral or racemic claims.",
        "No SI file was supplied or associated.",
    ],
}

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "coverage/pages_0007_0011.json", coverage)
st.write(RUN / "raw_model/pages_0007_0011_raw_extraction.json", raw)
st.write(RUN / "raw_model/pages_0007_0011_standard_stereo_checks.json", stereo_checks)
st.write(
    RUN / "raw_model/pages_0007_0011_source_mentions.json",
    st.source_mentions(DATASET, ["4.5.1", "402.4 mg", "374.2 mg", "3gb", "Na2SO4"], context=220),
)
print(json.dumps({
    "incoming_reactions": len(bundle["reactions"]),
    "updated_existing": len(incoming_reaction_ids) - 2,
    "unchanged_dependency_records": len(dependency_reaction_ids),
    "new_reactions": 2,
    "molecules_in_submission": len(bundle["molecules"]),
    "new_molecules": 1,
    "pages": [7, 8, 9, 10, 11],
}, ensure_ascii=False))
