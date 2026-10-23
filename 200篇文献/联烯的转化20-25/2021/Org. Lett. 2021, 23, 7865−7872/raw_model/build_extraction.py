#!/usr/bin/env python3
"""Build source-derived v3.1 records for doc_e47326ab6eaa3ee5.

This is an extraction artifact, not a modification of the shared tooling.
All structures that were not safely encoded from the supplied pages remain unresolved.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_e47326ab6eaa3ee5"
RUN = ROOT / "runs/current"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def F(raw=None, norm=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, norm, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def qty(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def raw_value(value):
    if isinstance(value, dict):
        if value.get("raw_text"):
            return value["raw_text"]
        if "number" in value and "unit" in value:
            return f"{value['number']} {value['unit']}"
        if value.get("ratio"):
            return value["ratio"]
    return str(value)


bundle = json.loads((RUN / "submitted.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
title = "Asymmetric Domino Heck/Dearomatization Reaction of β-Naphthols to Construct Indole–Terpenoid Frameworks"
dataset["documents"][0]["title"] = title
bundle["documents"] = dataset["documents"]

molecules = {}


def add_molecule(label, evidence, name=None, stereo="unresolved", smiles=None,
                 identity_note=None, assignment=None):
    mid = f"{DOC}:mol:{label}"
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name, name, evidence) if name else F()
    m["names"]["abbreviations"] = F(label, [label], evidence) if label != (name or label) else F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    deriv = "Status assigned from the inspected source drawing; absence of a label was not treated as proof of achirality."
    if stereo == "achiral":
        chirality, geometry, rep = "achiral", "not_applicable", "not_encoded"
        deriv = "The inspected scaffold has no tetrahedral stereocenter; for allenamides 2 the terminal CH2 end prevents an allene stereogenic axis."
    elif stereo == "product":
        chirality, geometry, rep = "partially_specified", "not_applicable", "not_encoded"
        m["stereochemistry"]["atom_centers"] = F(
            "product drawing uses a hashed bond at the newly formed quaternary center",
            [{"center": "new all-carbon quaternary center", "configuration": "drawn but absolute label unreported"}], evidence)
    elif stereo == "3aa":
        chirality, geometry, rep = "specified", "not_applicable", "not_encoded"
        m["stereochemistry"]["atom_centers"] = F(
            "the absolute configuration of 3aa was determined to be the S-configuration by X-ray crystallographic analysis",
            [{"center": "new all-carbon quaternary center", "configuration": "S"}], evidence)
        m["stereochemistry"]["assignment"] = F(
            "S-configuration by X-ray crystallographic analysis", {"configuration": "S", "method": "X-ray crystallography"}, evidence)
    elif stereo == "negative_3aa":
        chirality, geometry, rep = "partially_specified", "not_applicable", "not_encoded"
        m["stereochemistry"]["mixture"] = F(
            "negative ee sign for 3aa in Table 1", {"major_enantiomer": "opposite sign to positive-ee 3aa", "absolute_configuration": None}, evidence)
        m["stereochemistry"]["assignment"] = F(
            "negative ee sign; absolute configuration not stated for these rows",
            {"relative_assignment": "opposite chromatographic sign", "absolute_configuration": None}, evidence)
    elif stereo == "chiral_ligand":
        chirality, geometry, rep = "partially_specified", "not_applicable", "not_encoded"
        m["stereochemistry"]["axial"] = F(
            "chiral phosphoramidite ligand drawing", [{"axis": "ligand stereogenic element", "configuration": "not textually assigned"}], evidence)
    else:
        chirality, geometry, rep = "unresolved", "unresolved", "unresolved"
    for key, val in (("chirality_status", chirality), ("geometry_status", geometry), ("representation_status", rep)):
        m["stereochemistry"][key] = F(None, val, evidence, method="inference", derivation=deriv)
    if assignment and stereo not in {"3aa"}:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = F(None, smiles, evidence, method="structure_derivation",
                         derivation="Candidate connectivity transcribed from the inspected drawing and checked only for graph parsing; not source-certified.")
        m["stereochemistry"]["isomeric_smiles"] = F(None, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation", derivation="RDKit canonicalization of the candidate source-derived SMILES.")
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evidence,
            method="structure_derivation", derivation="RDKit stereo graph key for the candidate SMILES; not a sample identity.")
        if label == "2a":
            m["stereochemistry"]["representation_status"] = F(None, "partial", evidence, method="inference",
                derivation="Connectivity is encoded, while toolkit diagnostics are not used to certify allene stereochemistry.")
        else:
            m["stereochemistry"]["representation_status"] = F(None, "encoded", evidence, method="inference",
                derivation="The achiral candidate connectivity is encoded; this is not source structure certification.")
    if identity_note:
        m["identity_note"] = F(None, identity_note, evidence, method="inference", derivation=identity_note)
    molecules[label] = m
    return mid


# Reactants and products. Only three low-risk connectivity candidates are encoded.
for label in [f"1{x}" for x in "abcdefghijk"]:
    smi = "Cc1c(O)c(C)c2ccccc2c1" if label == "1a" else "Cc1c(O)cc2ccccc2c1" if label == "1f" else None
    add_molecule(label, [ev(3 if label == "1a" else 4)], name=f"β-naphthol substrate {label}", stereo="achiral", smiles=smi,
                 identity_note="Dense scope drawing retained by paper label; complete SMILES omitted unless separately encoded.")
for label in [f"2{x}" for x in "abcdefghijklmnopqrstu"]:
    smi = "C=C=CN(c1ccccc1I)S(=O)(=O)c1ccc(C)cc1" if label == "2a" else None
    add_molecule(label, [ev(3 if label == "2a" else 4)], name=f"N-(2-iodophenyl)allenamide {label}", stereo="achiral", smiles=smi,
                 identity_note="The allenamide has a terminal CH2 allene end; substituent details are retained by source label when no SMILES is supplied.")

product_labels = ["3aa", "3ab", "3ac", "3ad", "3ae", "3af", "3ag", "3ah", "3ai", "3aj", "3ak"]
product_labels += ["3ba", "3ca", "3da", "3ea", "3fa", "3ga", "3ha", "3ia", "3ja", "3ka"]
product_labels += ["3al", "3am", "3an", "3ao", "3ap", "3aq", "3ar", "3as", "3at", "3au"]
for label in product_labels:
    page = 3 if label == "3aa" else 4
    add_molecule(label, [ev(page)], name=f"dearomatized product {label}", stereo="3aa" if label == "3aa" else "product",
                 identity_note="Product structure is retained by source label and drawing; complete stereochemical SMILES was not safely encoded from the main-text image.")
add_molecule("3aa_negative_ee", [ev(3)], name="3aa, opposite-sign enantiomer-enriched sample", stereo="negative_3aa",
             identity_note="Kept distinct from positive-ee/S-assigned 3aa; the negative ee sign was not converted into an unsupported R assignment.")
for label in [str(x) for x in range(5, 13)]:
    add_molecule(label, [ev(5)], name=f"derivatized product {label}", stereo="product",
                 identity_note="Stereochemical drawing and experiment-specific ee/dr are retained; complete stereochemical SMILES was not safely encoded.")

# Catalysts, ligands, reagents, solvents, and generic Scheme 1 identities.
agent_specs = {
    "Pd_allyl_dimer": "[Pd(C3H5)Cl]2", "Cs2CO3": "cesium carbonate", "PhCH3": "toluene",
    "EtOAc": "ethyl acetate", "petroleum_ether": "petroleum ether", "NaBH4": "sodium borohydride",
    "MeOH": "methanol", "MeMgBr": "methylmagnesium bromide", "THF": "tetrahydrofuran",
    "TMSCH2Li": "((trimethylsilyl)methyl)lithium", "Pd_C": "palladium on carbon", "H2": "hydrogen",
    "EtOH": "ethanol", "Mg": "magnesium", "NH4Cl": "ammonium chloride", "benzene": "benzene",
    "PCC": "pyridinium chlorochromate", "CH2Cl2": "dichloromethane", "Co_acac2": "cobalt(II) acetylacetonate",
    "TBHP": "tert-butyl hydroperoxide", "acetone": "acetone", "KOH": "potassium hydroxide",
    "generic_allyl_reagent": "generic allylic reagent", "generic_naphthol": "generic β-naphthol",
    "generic_prior_product": "generic allylic dearomatization product", "generic_Nu": "generic nucleophile",
    "generic_prior_allenamide": "generic N-(2-iodophenyl)allenamide", "generic_prior_indole": "generic 3-substituted indole",
}
for label, name in agent_specs.items():
    page = 2 if label.startswith("generic_") else 3 if label in {"Pd_allyl_dimer", "Cs2CO3", "PhCH3", "EtOAc", "petroleum_ether"} else 5
    add_molecule(label, [ev(page)], name=name, stereo="unresolved")
for i in range(1, 13):
    add_molecule(f"L{i}", [ev(3)], name=f"ligand L{i}", stereo="chiral_ligand" if i >= 4 else "unresolved",
                 identity_note="Ligand structure retained by source drawing; no complete stereochemical SMILES asserted.")


def mid(label):
    return molecules[label]["molecule_id"]


def participant(pid, label, role, evidence, stage="reaction", amount=None, mass=None, volume=None,
                equiv=None, molpct=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid(label), role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    for key, val in (("amount", amount), ("mass", mass), ("volume", volume),
                     ("equivalents", equiv), ("mol_percent", molpct)):
        if val is not None:
            p["charges"][key] = F(raw_value(val), val, evidence)
    if reference:
        p["equivalence_reference"] = F(reference, mid(reference), evidence)
    return p


def manual_classification(r, evidence, changes, coverage="partial", limitation=None):
    facts = []
    labels = {}
    for i, (a, b, before, after, change_type) in enumerate(changes, 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": change_type}
        facts.append(F(None, change, evidence, method="structure_derivation",
                       derivation="Partial net bond change inferred from the inspected source structures; atom-complete mapping is unavailable."))
        if change_type in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(a, b, after)
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = F(None, list(labels.values()), evidence, method="structure_derivation",
        derivation="Taxonomy labels are generated from the explicitly listed partial net bond changes.")
    r["classification"]["coverage"] = F(None, coverage, evidence, method="inference",
        derivation="No atom-mapped complete reaction was asserted; only source-evident transformation classes are retained.")
    r["classification"]["limitations"] = [limitation or "Atom-complete mapping, coproducts, and full hydrogen balance remain unresolved."]


def base_reaction(suffix, evidence, record_type="specific_experiment", categories=None, location=None,
                  name="Asymmetric domino Heck/dearomatization reaction"):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence) if categories is not None else F()
    r["location"] = F(location, location, evidence) if location is not None else F()
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(None, "domino Heck cyclization/dearomatization", evidence,
        method="inference", derivation="Normalized from the article's named transformation.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["No explicitly associated SI: detailed operation, workup, purification, and spectra remain unavailable.",
                                  "Machine extraction and visual transcription require expert chemical review."]}
    return r


def add_common_domino(r, evidence, substrate1, substrate2, product, pd_pct=1, ligand="L5", ligand_pct=1,
                      temperature=0, time=12, scale1=0.1, scale2=0.12, volume=1.0, concentration=0.1,
                      atmosphere="N2"):
    r["participants"] = [
        participant("p_1", substrate1, "reactant", evidence, amount=qty(scale1, "mmol")),
        participant("p_2", substrate2, "reactant", evidence, amount=qty(scale2, "mmol")),
        participant("p_pd", "Pd_allyl_dimer", "catalyst", evidence, molpct=qty(pd_pct, "mol%"), reference=substrate1),
        participant("p_lig", ligand, "ligand", evidence, molpct=qty(ligand_pct, "mol%"), reference=substrate1),
        participant("p_base", "Cs2CO3", "reagent", evidence, equiv=qty(1.5, "equiv"), reference=substrate1),
        participant("p_solvent", "PhCH3", "solvent", evidence, volume=qty(volume, "mL")),
        participant("p_product", product, "product", evidence, product_role="main_product"),
    ]
    r["reaction_scale"] = F(f"{substrate1} ({scale1} mmol)", {"reference_molecule_id": mid(substrate1), "amount": qty(scale1, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = F("PhCH3", [mid("PhCH3")], evidence)
    r["solvent_conditions"]["total_volume"] = F(f"{volume} mL", qty(volume, "mL"), evidence)
    r["solvent_conditions"]["concentration"] = F(f"c = {concentration} M", qty(concentration, "M"), evidence)
    r["solvent_conditions"]["concentration_basis"] = F(substrate1, mid(substrate1), evidence)
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    step["operation"] = F("reaction under the tabulated conditions", "combine and react", evidence)
    step["temperature"] = F(f"{temperature} °C", qty(temperature, "°C"), evidence)
    step["time"] = F(f"{time} h", qty(time, "h"), evidence)
    step["atmosphere"] = F(atmosphere, atmosphere, evidence)
    r["steps"] = [step]
    r["procedure"]["general_procedure_reference"] = F("table/scheme shared reaction conditions", "shared conditions shown in the main-text table or scheme", evidence)
    manual_classification(r, evidence, [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")],
                          limitation="Two framework C–C formations are source-evident for the domino sequence; atom identities, bond-order changes, and hydrogen balance are unresolved.")


def outcome(oid, kind, target, value, evidence, basis=None, scope=None, yield_type=None, product_role="main_product"):
    o = st.new_outcome(oid, kind, [mid(target)] if target else [], evidence)
    o["product_role"] = F(product_role, product_role, evidence) if target else F()
    o["value"] = F(raw_value(value), value, evidence)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    return o


reactions = []

# Scheme 1: two prior-work families and one current-work generic scheme.
r = base_reaction("scheme1_prior_allylic_dearomatization", [ev(2)], "literature_example", None,
                  {"page": 2, "region": "Scheme 1a"}, "Allylic dearomatization reaction of naphthol (refs 7)")
r["participants"] = [participant("p1", "generic_allyl_reagent", "reactant", [ev(2)]),
                     participant("p2", "generic_naphthol", "reactant", [ev(2)]),
                     participant("p3", "generic_prior_product", "product", [ev(2)], product_role="generic_product")]
r["review"]["unresolved"].append("Generic cited-work scheme; no single experimental procedure or yield is reported here.")
manual_classification(r, [ev(2)], [("C", "C", 0, 1, "formed")], limitation="Generic cited-work drawing supports only a new C–C connection.")
reactions.append(r)

r = base_reaction("scheme1_prior_domino_heck", [ev(2)], "literature_example", None,
                  {"page": 2, "region": "Scheme 1b"}, "Domino Heck reactions of N-(2-iodophenyl)allenamide (refs 10–14)")
r["participants"] = [participant("p1", "generic_prior_allenamide", "reactant", [ev(2)]),
                     participant("p2", "generic_Nu", "reactant", [ev(2)]),
                     participant("p3", "generic_prior_indole", "product", [ev(2)], product_role="generic_product")]
r["review"]["unresolved"].append("Nu is generic and can represent different bond-forming classes; classification is unresolved.")
r["classification"]["coverage"] = F(None, "unresolved", [ev(2)], method="inference", derivation="Generic Nu prevents a unique bond classification.")
reactions.append(r)

r = base_reaction("scheme1_current_generic", [ev(1), ev(2)], "generic_scheme", None,
                  [{"page": 1, "region": "graphical abstract"}, {"page": 2, "region": "Scheme 1c"}])
r["participants"] = [participant("p1", "generic_prior_allenamide", "reactant", [ev(1), ev(2)]),
                     participant("p2", "generic_naphthol", "reactant", [ev(1), ev(2)]),
                     participant("p3", "generic_prior_product", "product", [ev(1), ev(2)], product_role="generic_product")]
manual_classification(r, [ev(1), ev(2)], [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")])
reactions.append(r)

# Table 1 optimization rows.
opt_rows = [
    (1, "L1", 5, 10, 0, "trace", None), (2, "L2", 5, 10, 0, "trace/complex", None),
    (3, "L3", 5, 10, 0, "trace", None), (4, "L4", 5, 10, 0, 80, 30),
    (5, "L5", 5, 10, 0, 81, 88), (6, "L5", 5, 5, 0, 88, 92),
    (7, "L6", 5, 5, 0, 51, 87), (8, "L7", 5, 5, 0, 82, 90),
    (9, "L8", 5, 5, 0, 25, 88), (10, "L9", 5, 5, 0, 47, 85),
    (11, "L10", 5, 5, 0, 76, {"number": 29, "unit": "%", "sign": "negative", "raw_text": "−29% ee"}),
    (12, "L11", 5, 5, 0, 87, {"number": 80, "unit": "%", "sign": "negative", "raw_text": "−80% ee"}),
    (13, "L12", 5, 5, 0, 85, 86), (14, "L5", 2.5, 2.5, 0, 90, 92),
    (15, "L5", 1, 1, 0, 89, 93), (16, "L5", 1, 1, 25, 84, 83),
    (17, "L5", 1, 1, -20, 37, 93), (18, "L5", 1, 1, 0, 71, {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99% ee"}),
]
for entry, lig, x, y, temp, yld, ee in opt_rows:
    evidence = [ev(3)] + ([ev(4)] if entry == 15 else [])
    loc = [{"page": 3, "region": f"Table 1 entry {entry}"}]
    if entry == 15:
        loc.append({"page": 4, "region": "Table 2 product 3aa duplicate report"})
    product = "3aa_negative_ee" if entry in (11, 12) else "3aa"
    r = base_reaction(f"table1_entry_{entry:02d}", evidence, categories=["optimization"], location=loc)
    add_common_domino(r, evidence, "1a", "2a", product, pd_pct=x, ligand=lig, ligand_pct=y, temperature=temp)
    if isinstance(yld, str):
        r["outcomes"] = [outcome("trace_product", "trace", product, yld, evidence, scope="Table 1 yield column")]
        r["review"]["unresolved"].append("Trace/complex is qualitative; detection limit and isolation basis for this row are unreported.")
    else:
        r["outcomes"] = [outcome("yield", "yield", product, qty(yld, "%"), evidence,
                                  basis="isolated", scope="reaction product", yield_type="isolated yield")]
    if ee is not None:
        ee_value = qty(ee, "%") if isinstance(ee, (int, float)) else ee
        r["outcomes"].append(outcome("ee", "selectivity", product, ee_value, evidence,
                                      basis="chiral HPLC", scope={"metric": "ee", "major_configuration": "S" if entry in (15, 18) else None}))
    if entry == 18:
        r["participants"] += [participant("p_etoac", "EtOAc", "purification_agent", evidence, stage="purification"),
                              participant("p_pet", "petroleum_ether", "purification_agent", evidence, stage="purification")]
        s = st.new_step("crystallization_enrichment", evidence, stage="purification")
        s["participant_ids"] = ["p_product", "p_etoac", "p_pet"]
        s["operation"] = F("enrichment via crystallization from EtOAc/petroleum ether", "crystallization enrichment", evidence)
        r["steps"].append(s)
        r["review"]["unresolved"].append("Crystallization solvent ratio, temperature, time, recovered racemate amount, and mother-liquor handling are unreported in the main text.")
    reactions.append(r)

# Table 2; 3aa is the same experiment as Table 1 entry 15 and is not duplicated.
table2 = [
    ("3ab", "2b", "H", 88, 92), ("3ac", "2c", "4-tBu", 92, 90), ("3ad", "2d", "4-MeO", 74, 90),
    ("3ae", "2e", "3,5-Me2", 72, 92), ("3af", "2f", "2,4,6-Me3", 67, 69), ("3ag", "2g", "4-NO2", 38, 83),
    ("3ah", "2h", "naphthalenesulfonyl", 90, 91), ("3ai", "2i", "Ac", 70, 49),
    ("3aj", "2j", "Boc", 63, 68), ("3ak", "2k", "Cbz", 50, 70),
]
for product, sub2, rgroup, yld, ee_pct in table2:
    evidence = [ev(4)] + ([ev(3)] if product in {"3af", "3ag", "3ah"} else [])
    r = base_reaction(f"table2_{product}", evidence, categories=["substrate_scope"],
                      location={"page": 4, "region": f"Table 2 {product}"})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    add_common_domino(r, evidence, "1a", sub2, product)
    r["notes"].append(F(f"N-protecting group shown as {rgroup}", {"N_protecting_group": rgroup}, evidence))
    r["outcomes"] = [outcome("yield", "yield", product, qty(yld, "%"), evidence, basis="isolated", yield_type="isolated yield"),
                     outcome("ee", "selectivity", product, qty(ee_pct, "%"), evidence, basis="chiral HPLC", scope={"metric": "ee", "absolute_configuration": None})]
    reactions.append(r)

# Table 3 naphthol scope and allenamide aryl scope.
table3_naph = [
    ("3ba", "1b", "Et", 91, 93), ("3ca", "1c", "allyl", 91, 92), ("3da", "1d", "Ph", 81, 90),
    ("3ea", "1e", "Bn", 62, 85), ("3fa", "1f", "H", 79, 82), ("3ga", "1g", "Cl", 82, 24),
    ("3ha", "1h", "Br", 81, 3), ("3ia", "1i", "6-Ph isomer", 72, 80), ("3ja", "1j", "6-Ph isomer", 64, 88),
    ("3ka", "1k", "2,6-naphthalenediol derivative", 86, 91),
]
table3_allene = [
    ("3al", "2l", "F", 86, 93), ("3am", "2m", "Cl", 90, 94), ("3an", "2n", "Br", 81, 91),
    ("3ao", "2o", "Me", 76, 87), ("3ap", "2p", "F", 84, 92), ("3aq", "2q", "Cl", 88, 93),
    ("3ar", "2r", "Br", 80, 93), ("3as", "2s", "CF3", 82, 85), ("3at", "2t", "MeO", 80, 86),
    ("3au", "2u", "CO2Me", 78, 85),
]
for product, sub1, substituent, yld, ee_pct in table3_naph:
    evidence = [ev(4)] + ([ev(5)] if product in {"3ga", "3ha", "3ia", "3ja", "3ka"} else [])
    r = base_reaction(f"table3_{product}", evidence, categories=["substrate_scope"], location={"page": 4, "region": f"Table 3 {product}"})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    add_common_domino(r, evidence, sub1, "2a", product)
    r["notes"].append(F(f"β-naphthol substituent shown as {substituent}", {"substituent": substituent}, evidence))
    r["outcomes"] = [outcome("yield", "yield", product, qty(yld, "%"), evidence, basis="isolated", yield_type="isolated yield"),
                     outcome("ee", "selectivity", product, qty(ee_pct, "%"), evidence, basis="chiral HPLC", scope={"metric": "ee", "absolute_configuration": None})]
    reactions.append(r)
for product, sub2, substituent, yld, ee_pct in table3_allene:
    evidence = [ev(4), ev(5)]
    r = base_reaction(f"table3_{product}", evidence, categories=["substrate_scope"], location={"page": 4, "region": f"Table 3 {product}"})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    add_common_domino(r, evidence, "1a", sub2, product)
    r["notes"].append(F(f"allenamide aryl substituent shown as {substituent}", {"substituent": substituent}, evidence))
    r["outcomes"] = [outcome("yield", "yield", product, qty(yld, "%"), evidence, basis="isolated", yield_type="isolated yield"),
                     outcome("ee", "selectivity", product, qty(ee_pct, "%"), evidence, basis="chiral HPLC", scope={"metric": "ee", "absolute_configuration": None})]
    reactions.append(r)

# Scheme 2 scale-up experiments; enrichment is a downstream purification stage in the same experiment.
scale_rows = [
    ("3aa", "1a", qty(1.03, "g"), 90, 93, 72, {"number": 99, "unit": "%", "qualifier": "gt", "raw_text": ">99% ee"}),
    ("3fa", "1f", qty(0.91, "g"), 83, 81, 56, qty(96, "%")),
]
for product, sub1, mass, yld, ee_before, enriched_yld, ee_after in scale_rows:
    evidence = [ev(5), ev(6)]
    r = base_reaction(f"scheme2_scale_{product}", evidence, categories=["scale_up"], location={"page": 5, "region": f"Scheme 2 {product}"})
    add_common_domino(r, evidence, sub1, "2a", product, scale1=2.5, scale2=3.0, volume=None if False else 25.0, concentration=0.1)
    # Scheme 2 reports c=0.1 M, so 25 mL is calculated rather than source-reported.
    r["participants"][5]["charges"]["volume"] = F()
    r["solvent_conditions"]["total_volume"] = F(None, qty(25.0, "mL"), evidence, method="inference",
        derivation="Calculated from 2.5 mmol substrate and c = 0.1 M; the scheme does not print a solvent volume.")
    r["participants"].append(participant("p_enrichment_solvent", "EtOAc", "purification_agent", evidence, stage="purification"))
    s = st.new_step("enrichment", evidence, stage="purification")
    s["participant_ids"] = ["p_product", "p_enrichment_solvent"]
    s["operation"] = F("after enrichment", "enrichment; method details unreported in Scheme 2", evidence)
    r["steps"].append(s)
    r["outcomes"] = [
        outcome("product_mass", "product_amount", product, mass, evidence, scope="before enrichment"),
        outcome("yield_before", "yield", product, qty(yld, "%"), evidence, scope="before enrichment", yield_type="reported yield"),
        outcome("ee_before", "selectivity", product, qty(ee_before, "%"), evidence, scope={"metric": "ee", "stage": "before enrichment"}),
        outcome("yield_after", "yield", product, qty(enriched_yld, "%"), evidence, scope="after enrichment", yield_type="reported yield"),
        outcome("ee_after", "selectivity", product, ee_after, evidence, scope={"metric": "ee", "stage": "after enrichment"}),
    ]
    for o in r["outcomes"]:
        if o["kind"] == "selectivity": o["measurement_basis"] = F("ee", "method not stated in Scheme 2; paper elsewhere uses chiral HPLC", evidence, status="pending_review")
        elif o["kind"] == "yield": o["measurement_basis"] = F()
    r["review"]["unresolved"].append("Scheme 2 does not state the yield measurement basis or the enrichment procedure; the calculated 25 mL solvent volume is not a reported charge.")
    reactions.append(r)


def derivative(suffix, start, products, reagents, solvents, temp=None, pressure=None, evidence=(ev(5), ev(6)),
               class_changes=(), name="Transformation of a dearomatized product"):
    r = base_reaction(suffix, list(evidence), categories=["derivatization"], location={"page": 5, "region": suffix}, name=name)
    r["transformation_type"] = F(None, name, evidence, method="inference", derivation="Normalized from Scheme 3/4 structures and reagents.")
    r["participants"] = [participant("p_start", start, "reactant", evidence)]
    for i, lab in enumerate(reagents, 1):
        role = "catalyst" if lab in {"Pd_C", "Co_acac2"} else "reagent"
        r["participants"].append(participant(f"p_reagent_{i}", lab, role, evidence))
    for i, lab in enumerate(solvents, 1):
        r["participants"].append(participant(f"p_solvent_{i}", lab, "solvent", evidence))
    for i, lab in enumerate(products, 1):
        r["participants"].append(participant(f"p_product_{i}", lab, "product", evidence, product_role="main_product" if i == 1 else "coproduct"))
    r["solvent_conditions"]["composition"] = F("/".join(solvents), [mid(s) for s in solvents], evidence) if solvents else F()
    if solvents == ["MeOH", "benzene"]:
        r["solvent_conditions"]["mixture_ratio"] = F("MeOH/C6H6 (1:1)", {"components": [mid("MeOH"), mid("benzene")], "ratio": "1:1"}, evidence)
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    step["operation"] = F("Scheme 3/4 arrow", name, evidence)
    if temp is not None: step["temperature"] = F(f"{temp} °C", qty(temp, "°C"), evidence)
    if pressure is not None: step["pressure"] = F(pressure, pressure, evidence)
    r["steps"] = [step]
    r["procedure"]["original_text"] = F("main-text scheme and prose; detailed procedure deferred to SI", "main-text-only transformation", evidence)
    manual_classification(r, evidence, list(class_changes))
    r["review"]["unresolved"].append("Reaction scale, reagent doses, duration, workup, purification, and yield measurement basis are not reported in the supplied main-text pages.")
    return r


# Scheme 3.
r = derivative("scheme3_3fa_to_5", "3fa", ["5"], ["NaBH4"], ["MeOH"],
               class_changes=[("C", "O", 2, 1, "order_decreased"), ("C", "H", 0, 1, "hydrogen_gain"), ("O", "H", 0, 1, "hydrogen_gain")],
               name="carbonyl reduction")
r["outcomes"] = [outcome("yield", "yield", "5", qty(92, "%"), [ev(5)], scope="product 5"),
                 outcome("ee", "selectivity", "5", qty(96, "%"), [ev(5)], scope={"metric": "ee"}),
                 outcome("dr", "selectivity", "5", {"ratio": ">20:1", "qualifier": "gt"}, [ev(5)], scope={"metric": "dr", "order": "major:minor"})]
reactions.append(r)
r = derivative("scheme3_3fa_to_6", "3fa", ["6"], ["MeMgBr"], ["THF"], temp=0,
               class_changes=[("C", "O", 2, 1, "order_decreased"), ("C", "C", 0, 1, "formed"), ("O", "H", 0, 1, "hydrogen_gain")],
               name="methylmagnesium bromide addition")
r["outcomes"] = [outcome("yield", "yield", "6", qty(90, "%"), [ev(5)], scope="product 6"), outcome("ee", "selectivity", "6", qty(93, "%"), [ev(5)], scope={"metric": "ee"}), outcome("dr", "selectivity", "6", {"ratio": ">20:1", "qualifier": "gt"}, [ev(5)], scope={"metric": "dr", "order": "major:minor"})]
reactions.append(r)
r = derivative("scheme3_3fa_to_7", "3fa", ["7"], ["TMSCH2Li"], ["THF"], temp=0,
               class_changes=[("C", "O", 2, 1, "order_decreased"), ("C", "C", 0, 1, "formed"), ("O", "H", 0, 1, "hydrogen_gain")],
               name="((trimethylsilyl)methyl)lithium addition")
r["outcomes"] = [outcome("yield", "yield", "7", qty(86, "%"), [ev(5)], scope="product 7"), outcome("ee", "selectivity", "7", qty(95, "%"), [ev(5)], scope={"metric": "ee"}), outcome("dr", "selectivity", "7", {"ratio": ">20:1", "qualifier": "gt"}, [ev(5)], scope={"metric": "dr", "order": "major:minor"})]
reactions.append(r)
r = derivative("scheme3_hydrogenation_3fa", "3fa", ["8", "9"], ["Pd_C", "H2"], ["EtOH"], pressure="1 atm",
               class_changes=[("C", "C", 2, 1, "order_decreased"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "O", 2, 1, "order_decreased"), ("O", "H", 0, 1, "hydrogen_gain")],
               name="catalytic hydrogenation")
r["outcomes"] = [outcome("yield_8", "yield", "8", qty(52, "%"), [ev(5)], scope="cyclohexanone 8"), outcome("ee_8", "selectivity", "8", qty(95, "%"), [ev(5)], scope={"metric": "ee"}), outcome("yield_9", "yield", "9", qty(39, "%"), [ev(5)], scope="cyclohexanol 9"), outcome("ee_9", "selectivity", "9", qty(95, "%"), [ev(5)], scope={"metric": "ee"}), outcome("dr_9", "selectivity", "9", {"ratio": ">20:1", "qualifier": "gt"}, [ev(5)], scope={"metric": "dr", "order": "major:minor"})]
reactions.append(r)

# Scheme 4.
r = derivative("scheme4_8_to_9", "8", ["9"], ["NaBH4"], ["MeOH"],
               class_changes=[("C", "O", 2, 1, "order_decreased"), ("C", "H", 0, 1, "hydrogen_gain"), ("O", "H", 0, 1, "hydrogen_gain")], name="ketone reduction")
r["outcomes"] = [outcome("yield", "yield", "9", qty(99, "%"), [ev(5)], scope="product 9"), outcome("ee", "selectivity", "9", qty(95, "%"), [ev(5)], scope={"metric": "ee"}), outcome("dr", "selectivity", "9", {"ratio": ">20:1", "qualifier": "gt"}, [ev(5)], scope={"metric": "dr", "order": "major:minor"})]
r["notes"].append(F("quantitative yield without erosion of enantiomeric purity", "author prose description corresponding to the 99% scheme value", [ev(6)]))
reactions.append(r)
r = derivative("scheme4_8_to_10", "8", ["10"], ["Mg", "NH4Cl"], ["MeOH", "benzene"],
               class_changes=[("N", "S", 1, 0, "broken"), ("N", "H", 0, 1, "hydrogen_gain")], name="detosylation")
r["outcomes"] = [outcome("yield", "yield", "10", qty(82, "%"), [ev(5)], scope="product 10"), outcome("ee", "selectivity", "10", qty(96, "%"), [ev(5)], scope={"metric": "ee"})]
reactions.append(r)
r = derivative("scheme4_9_to_8", "9", ["8"], ["PCC"], ["CH2Cl2"],
               class_changes=[("C", "O", 1, 2, "order_increased"), ("C", "H", 1, 0, "hydrogen_loss"), ("O", "H", 1, 0, "hydrogen_loss")], name="alcohol oxidation")
r["outcomes"] = [outcome("yield", "yield", "8", qty(97, "%"), [ev(5)], scope="product 8"), outcome("ee", "selectivity", "8", qty(95, "%"), [ev(5)], scope={"metric": "ee"})]
reactions.append(r)
r = derivative("scheme4_9_to_11", "9", ["11"], ["Co_acac2", "TBHP"], ["acetone"],
               class_changes=[("C", "O", 0, 2, "formed"), ("C", "H", 1, 0, "hydrogen_loss")], name="benzylic oxidation")
r["outcomes"] = [outcome("yield", "yield", "11", qty(71, "%"), [ev(5)], scope="product 11"), outcome("ee", "selectivity", "11", qty(96, "%"), [ev(5)], scope={"metric": "ee"})]
reactions.append(r)
r = derivative("scheme4_11_to_12", "11", ["12"], ["KOH"], ["EtOH"],
               class_changes=[("C", "C", 1, 2, "order_increased"), ("C", "H", 1, 0, "hydrogen_loss"), ("C", "O", 1, 0, "broken")], name="base-promoted dehydration")
r["outcomes"] = [outcome("yield", "yield", "12", qty(92, "%"), [ev(5)], scope="product 12"), outcome("ee", "selectivity", "12", qty(96, "%"), [ev(5)], scope={"metric": "ee"})]
reactions.append(r)

# Figure 2 is an actual mechanistic/linearity study, not a proposed mechanism.
r = base_reaction("figure2_ee_linearity", [ev(6)], categories=["mechanistic_study"],
                  location={"page": 6, "region": "Figure 2"}, name="Linearity relationship between ee of ligand L5 and product 3aa")
r["participants"] = [participant("p1", "1a", "reactant", [ev(6)]), participant("p2", "2a", "reactant", [ev(6)]),
                     participant("pl", "L5", "ligand", [ev(6)]), participant("pp", "3aa", "product", [ev(6)], product_role="measured_product")]
s = st.new_step("linearity_series", [ev(6)])
s["participant_ids"] = ["p1", "p2", "pl", "pp"]
s["operation"] = F("determination of linear correlation between the ee values of product 3aa and L5", "series with varied ligand ee", [ev(6)])
r["steps"] = [s]
o = st.new_outcome("regression", "other", [mid("3aa")], [ev(6)])
o["value"] = F("y = 0.9423x - 2.4568; R² = 0.9957", {"slope": 0.9423, "intercept": -2.4568, "R2": 0.9957,
                  "x": "ee of L5 (%)", "y": "ee of product 3aa (%)"}, [ev(6)])
o["measurement_basis"] = F("Figure 2 plotted ee values", "method for individual ee points not restated in Figure 2", [ev(6)])
o["scope"] = F("clear linear correlation suggesting a monomeric chiral catalyst might be involved", "author mechanistic interpretation", [ev(6)])
r["outcomes"] = [o]
manual_classification(r, [ev(6)], [("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")])
r["review"]["unresolved"].append("Exact numeric coordinates and reaction conditions for each plotted point are not tabulated on the supplied page; they were not reverse-engineered from pixels.")
reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:ccdc", "evidence_ids": [ev(4), ev(6)], "kind": "crystallography",
     "value": "CCDC 2096151–2096152; 3aa absolute S assignment is stated in the main text. This is not a synthetic experiment."},
    {"context_id": f"{DOC}:context:si_notice", "evidence_ids": [ev(6)], "kind": "availability_notice",
     "value": "The paper advertises Supporting Information, but no SI file was explicitly supplied or associated in this session."},
    {"context_id": f"{DOC}:context:summary", "evidence_ids": [ev(1), ev(6)], "kind": "author_summary",
     "value": "Reported scope reaches up to 92% yield and 94% ee; summary text is not counted as a separate experiment."},
]

# Operational task state records visual review attempts without claiming semantic completeness.
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    page = int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0])
    if 1 <= page <= 6:
        task["state"] = "attempted"
        task["reason"] = "Original rendered page inspected; findings and unresolved items are recorded in coverage/pages_0001_0006.json."

raw = {
    "document_id": DOC,
    "pages_visually_inspected": [1, 2, 3, 4, 5, 6],
    "table1_rows": opt_rows,
    "table2_rows_excluding_duplicate_3aa": table2,
    "table3_naphthol_rows": table3_naph,
    "table3_allenamide_rows": table3_allene,
    "scale_rows": scale_rows,
    "derivatization_arrows": ["3fa→5", "3fa→6", "3fa→7", "3fa→8+9", "8→9", "8→10", "9→8", "9→11", "11→12"],
    "decisions": [
        "Table 2 product 3aa duplicates Table 1 entry 15 and is evidence-linked, not duplicated.",
        "Each Table 1 row is an independent optimization experiment.",
        "Scheme 2 before/after enrichment values are outcomes/stages of two scale-up reactions.",
        "Figure 2 is retained as a mechanistic study series; plot coordinates were not digitized.",
        "Dense source structures without safely encoded stereochemical SMILES remain unresolved rather than guessed.",
    ],
}

raw["literal_mentions"] = st.source_mentions(dataset, ["1a", "2a", "3aa", "3fa", "3ha", "3au"], context=180)

def page_record_ids(page):
    eid = ev(page)
    return [r["reaction_id"] for r in reactions if eid in r["evidence_ids"]]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, authors, citation and abstract", "finding": "Article identity and current-work summary; no separate experiment created from the abstract.", "evidence_ids": [ev(1)]},
                {"region": "graphical abstract", "finding": "Current-work generic asymmetric domino Heck/dearomatization scheme; duplicate presentation of Scheme 1c, linked to one generic record.", "evidence_ids": [ev(1)]},
                {"region": "Figure 1 and background prose", "finding": "Natural-product examples and background discussion, not experiments from this paper.", "evidence_ids": [ev(1)]},
            ],
            "record_ids": page_record_ids(1),
            "unresolved": ["Graphical-abstract generic structures are not atom-mapped and do not establish a complete reaction mapping."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 1a", "finding": "Generic cited allylic dearomatization family (ref. 7), classified as prior literature rather than a current-paper experiment.", "evidence_ids": [ev(2)]},
                {"region": "Scheme 1b", "finding": "Generic cited domino Heck/Nu trapping family (refs. 10–14); Nu is variable and classification remains unresolved.", "evidence_ids": [ev(2)]},
                {"region": "Scheme 1c", "finding": "Current-work generic domino Heck/dearomatization concept; same generic record as the graphical abstract.", "evidence_ids": [ev(2)]},
                {"region": "lower-column optimization discussion", "finding": "Model substrates and entries 1–5 discussion; exact rows are recorded from Table 1 on page 3.", "evidence_ids": [ev(2)]},
            ],
            "record_ids": page_record_ids(2),
            "unresolved": ["The prose mentions additional Pd2(dba)3 ligand tests in Table S1, but SI was not provided; no unseen optimization rows were invented."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Table 1 reaction drawing and ligand structures", "finding": "Model reaction, desired 3aa, unobserved byproduct 4, and ligands L1–L12 inspected.", "evidence_ids": [ev(3)]},
                {"region": "Table 1 entries 1–18", "finding": "All 18 optimization rows recorded independently, including trace/complex, negative ee signs, −20 °C, and crystallization-enriched entry 18.", "evidence_ids": [ev(3)]},
                {"region": "Table 1 footnotes", "finding": "0.1/0.12 mmol charges, catalyst/ligand loadings, Cs2CO3, 1 mL PhCH3, N2, 12 h, isolated-yield basis, chiral-HPLC ee basis, and entry-18 enrichment captured.", "evidence_ids": [ev(3)]},
                {"region": "lower prose", "finding": "Optimized entry 15, crystallization behavior, S assignment of 3aa, and selected Table 2 results cross-linked.", "evidence_ids": [ev(3)]},
            ],
            "record_ids": page_record_ids(3),
            "unresolved": ["Entry 18 lacks crystallization ratio, temperature, duration, isolated racemate amount, and mother-liquor procedure.", "Negative ee signs identify opposite chromatographic sense but do not explicitly assign an absolute configuration."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "Table 2 shared arrow and footnote", "finding": "Shared conditions and isolated-yield/chiral-HPLC footnotes captured for the N-protecting-group scope.", "evidence_ids": [ev(4)]},
                {"region": "Table 2 products 3aa–3ak", "finding": "Eleven displayed products inspected; 3aa is the same 89%/93% experiment as Table 1 entry 15 and was evidence-merged, while 3ab–3ak are ten separate records.", "evidence_ids": [ev(4)]},
                {"region": "Table 3 shared arrow and footnote", "finding": "Shared conditions and measurement bases captured for β-naphthol and allenamide scope.", "evidence_ids": [ev(4)]},
                {"region": "Table 3 products 3ba–3ka and 3al–3au", "finding": "Twenty scope entries recorded independently with printed yields and ee values.", "evidence_ids": [ev(4)]},
                {"region": "lower prose", "finding": "Selected β-naphthol scope results cross-checked; continuation appears on page 5.", "evidence_ids": [ev(4)]},
            ],
            "record_ids": page_record_ids(4),
            "unresolved": ["Dense substrate/product drawings were not converted to complete stereochemical SMILES where manual encoding was unsafe.", "Absolute configurations other than textually assigned 3aa are not stated in the inspected main-text page."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2", "finding": "Two 2.5 mmol scale-up reactions (3aa and 3fa) captured with substrate charges, product masses, pre-enrichment yield/ee, and post-enrichment yield/ee.", "evidence_ids": [ev(5)]},
                {"region": "Scheme 3 arrows a–d", "finding": "Four distinct experiments captured: 3fa→5, 3fa→6, 3fa→7, and one hydrogenation giving 8 plus 9.", "evidence_ids": [ev(5)]},
                {"region": "Scheme 4a", "finding": "Two independent arrows captured: 8→9 reduction and 8→10 detosylation.", "evidence_ids": [ev(5)]},
                {"region": "Scheme 4b", "finding": "Three independent arrows captured: 9→8 oxidation, 9→11 benzylic oxidation, and 11→12 dehydration.", "evidence_ids": [ev(5)]},
                {"region": "lower scope prose", "finding": "Table 3 continuation cross-checks 3ha–3ka and 3al–3au ranges; no additional scope arrows occur here.", "evidence_ids": [ev(5)]},
            ],
            "record_ids": page_record_ids(5),
            "unresolved": ["Scheme 2 does not state enrichment details or yield measurement bases.", "Schemes 3–4 omit reaction scale, reagent amounts, durations, workup, purification, and yield measurement bases; these remain explicit missing fields."],
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [
                {"region": "left-column scale-up and derivatization prose", "finding": "Scheme 2–4 experiments cross-checked, including MeOH/benzene 1:1 for detosylation and qualitative descriptions such as quantitative yield.", "evidence_ids": [ev(6)]},
                {"region": "Figure 2 and adjacent prose", "finding": "Ligand-ee/product-ee linearity study recorded as a mechanistic-study series with printed regression and R²; individual point coordinates were not digitized.", "evidence_ids": [ev(6)]},
                {"region": "summary", "finding": "Author summary only; no duplicate reaction records created.", "evidence_ids": [ev(6)]},
                {"region": "Associated Content and Accession Codes", "finding": "SI availability notice and CCDC identifiers retained as source context, not as associated SI or synthetic experiments.", "evidence_ids": [ev(6)]},
            ],
            "record_ids": page_record_ids(6),
            "unresolved": ["Exact Figure 2 data-point coordinates and per-point reaction conditions are not tabulated.", "The advertised SI was not supplied, so procedures and characterization remain unavailable."],
        },
    ],
}

(RUN / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/extracted_rows.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"])}, ensure_ascii=False))
