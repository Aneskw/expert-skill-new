#!/usr/bin/env python3
"""Build the source-derived v3.1 submission for doc_6095e7fe387f2dc7.

This is an extraction artifact, not a chemistry validator.  Every record remains
machine-generated/needs_review after standard merge.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st
import reaction_tools as rt

DOC = "doc_6095e7fe387f2dc7"
E1 = f"{DOC}:page1:image"
E2 = f"{DOC}:page2:image"
E3 = f"{DOC}:page3:image"
E4 = f"{DOC}:page4:image"
E5 = f"{DOC}:page5:image"
T1 = f"{DOC}:page1:text"
T2 = f"{DOC}:page2:text"
T3 = f"{DOC}:page3:text"
T4 = f"{DOC}:page4:text"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", *, derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def missing(evidence=(), *, not_applicable=False):
    if not_applicable:
        return f(evidence=evidence, method="not_applicable", status="not_applicable")
    return f(evidence=evidence)


def mol(mid, label, evidence, *, name=None, smiles=None, chirality="unresolved",
        geometry="unresolved", representation="unresolved", stereo_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["names"]["english"] = f(name or label, name or label, evidence)
    m["names"]["abbreviations"] = f([], [], evidence)
    m["names"]["aliases"] = f([], [], evidence)
    m["cas_numbers"] = missing(evidence)
    m["external_properties"] = []
    if smiles:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "smiles"
        m["smiles"] = f({"source_drawing": label}, smiles, evidence,
                         method="structure_derivation",
                         derivation="Candidate connectivity transcribed from the cited source drawing; RDKit parsing does not verify source identity.")
        m["generic_template"] = None
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="RDKit stereo graph key calculated from the candidate SMILES; not a sample identity or source verification.")
        m["stereochemistry"]["isomeric_smiles"] = f(
            {"candidate": smiles}, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the source-derived candidate; unspecified stereo remains unspecified.")
    else:
        m["identity_status"] = "source_identified"
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = f(
        label, chirality, evidence, method="inference",
        derivation=stereo_note or "Stereo status conservatively assigned from the cited drawing and candidate graph; absence of a label is not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = f(
        label, geometry, evidence, method="inference",
        derivation=stereo_note or "Geometry status conservatively assigned from the cited drawing and candidate graph.")
    m["stereochemistry"]["representation_status"] = f(
        label, representation, evidence, method="inference",
        derivation=stereo_note or "Representation status reflects whether all source-visible stereochemical information is encoded.")
    m["stereochemistry"]["atom_centers"] = missing(evidence)
    m["stereochemistry"]["double_bonds"] = missing(evidence)
    m["stereochemistry"]["axial"] = missing(evidence)
    m["stereochemistry"]["mixture"] = missing(evidence)
    m["stereochemistry"]["assignment"] = missing(evidence)
    m["identity_note"] = f(
        {"source_label": label, "structure_status": "candidate" if smiles else "unresolved"},
        "Source-linked identity retained; any SMILES is a machine candidate and is not structure-verified.", evidence)
    return m


def template_mol(mid, label, evidence, fixed, variables, *, name=None):
    m = mol(mid, label, evidence, name=name or label)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f({"source_drawing": fixed}, {"description": fixed}, evidence),
        "variables": [
            {
                "label": var,
                "attachment_points": f(f"Source attachment point for {var}", {"label": var}, evidence),
                "definition": f(desc, desc, evidence),
                "allowed_values": missing(evidence),
                "constraints": missing(evidence),
            }
            for var, desc in variables
        ],
        "instance_mapping": [],
        "reason": f("Generic R/PG drawing", "The source uses variable groups, so no wildcard concrete SMILES is asserted.", evidence),
    }
    m["stereochemistry"]["chirality_status"] = f(label, "unresolved", evidence, method="inference",
        derivation="Generic substituents can alter chirality; the source overview does not define a single stereochemical identity.")
    m["stereochemistry"]["geometry_status"] = f(label, "unresolved", evidence, method="inference",
        derivation="Generic substituents can alter alkene geometry; the overview does not define one instance.")
    m["stereochemistry"]["representation_status"] = f(label, "unresolved", evidence, method="inference",
        derivation="Generic template is intentionally not a complete stereochemical representation.")
    return m


def simple_agent(mid, label, evidence, smiles=None, *, geometry="not_applicable"):
    if smiles:
        return mol(mid, label, evidence, name=label, smiles=smiles, chirality="achiral",
                   geometry=geometry, representation="not_applicable",
                   stereo_note="Candidate graph has no assigned stereocenter; geometry is not applicable unless explicitly marked otherwise. This is a graph derivation, not source identity verification.")
    return mol(mid, label, evidence, name=label, chirality="unresolved",
               geometry="unresolved", representation="unresolved")


def ar_hydrazone_smiles(ar):
    return f"Cc1ccc(S(=O)(=O)NN=C{ar})cc1"


def product_smiles(ar, diphenyl=False):
    scaffold = "C=CC1CN(S(=O)(=O)c2ccc(C)cc2)C=C1"
    if diphenyl:
        return scaffold + "C(=C(c1ccccc1)c1ccccc1)"  # intentionally replaced below if invalid
    return scaffold + "/C=C/" + ar


def make_product(mid, label, ar, evidence, *, diphenyl=False):
    if diphenyl:
        # The trisubstituted alkene is represented without an E/Z claim because
        # the hydrazone-derived carbon bears two identical phenyl groups.
        smi = "C=CC1CN(S(=O)(=O)c2ccc(C)cc2)C=C1C=C(c1ccccc1)c1ccccc1"
        m = mol(mid, label, evidence, name=label, smiles=smi, chirality="unspecified",
                geometry="not_applicable", representation="partial",
                stereo_note="The ring carbon bearing vinyl is a possible stereocenter but no absolute/relative assignment is printed; the diphenyl-substituted alkene has no E/Z descriptor.")
        return m
    smi = product_smiles(ar)
    m = mol(mid, label, evidence, name=label, smiles=smi, chirality="unspecified",
            geometry="specified", representation="partial",
            stereo_note="The scope drawing specifies the aryl/vinyl-pyrroline alkene as the drawn E geometry; a possible ring stereocenter is not assigned, so representation remains partial.")
    m["stereochemistry"]["double_bonds"] = f(
        {"source": "drawn trans alkene"}, [{"assignment": "E", "site": "hydrazone-derived exocyclic alkene"}], evidence,
        method="structure_derivation", derivation="E assignment is transcribed from the line-angle drawing; only 3a is additionally supported by the printed 1H NMR coupling statement.")
    return m


def participant(pid, mid, role, evidence, *, amount=None, volume=None, equivalents=None,
                mol_percent=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    if amount is not None:
        p["charges"]["amount"] = f(f"{amount:g} mmol", {"number": amount, "unit": "mmol"}, evidence)
    if volume is not None:
        p["charges"]["volume"] = f(f"{volume:g} mL", {"number": volume, "unit": "mL"}, evidence)
    if equivalents is not None:
        p["charges"]["equivalents"] = f(f"{equivalents:g} equiv", {"number": equivalents, "unit": "equiv"}, evidence)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(f"{mol_percent:g} mol%", {"number": mol_percent, "unit": "mol%"}, evidence)
    if reference is not None and (equivalents is not None or mol_percent is not None):
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def classified(reaction, evidence, limitation):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 2, "change_type": "formed"},
    ]
    reaction["classification"]["bond_changes"] = [
        f(None, c, evidence, method="structure_derivation",
          derivation="Partial bond change derived from the depicted reactant/product connectivity and the authors' statement that two C-C bonds are formed.")
        for c in changes
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b2"]},
    ]
    reaction["classification"]["labels"] = f(
        None, labels, evidence, method="structure_derivation",
        derivation="One new sigma C-C ring-closing connection and one newly connected C=C carbene-coupling bond are visible in the depicted overall transformation.")
    reaction["classification"]["coverage"] = f(
        None, "partial", evidence, method="inference",
        derivation="No complete atom mapping is asserted because eliminated hydrazone atoms, acetate, metal species, and hydrogen transfers are not fully represented.")
    reaction["classification"]["limitations"] = [limitation, "No complete atom mapping; implicit-H and coproduct accounting remain unresolved."]
    reaction["classification"]["model_votes"] = []


def base_reaction(rid, evidence, record_type, categories, location, name, transformation=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    if transformation:
        r["transformation_type"] = f(None, transformation, evidence, method="inference",
            derivation="Transformation description normalized from the cited source title, scheme, and product connectivity.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def optimization_reaction(row):
    n, catalyst, base, solvent, temp, hours, yld = row
    evidence = [E2]
    if n == 1:
        evidence = [E1, E2]
    if n == 18:
        evidence = [E2, E3]
    r = base_reaction(f"table1_entry_{n:02d}", evidence, "specific_experiment",
        ["optimization"] + (["substrate_scope"] if n == 18 else []),
        {"page": [2, 3] if n == 18 else ([1, 2] if n == 1 else [2]), "table": "Table 1", "entry": n,
         "duplicate_scope_location": "Table 2, 3a" if n == 18 else None},
        "Pd-catalyzed cascade metallo-ene cyclization/metallo-carbene coupling",
        "cascade metallo-ene cyclization/metallo-carbene coupling")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_overview"]
    cat_mid = catalyst_ids[catalyst]
    base_mid = "mol_tbuoli" if base == "tBuOLi" else "mol_ch3oli"
    solv_mid = solvent_ids[solvent]
    catalyst_amount = 0.01 if n == 21 else 0.02
    catalyst_mol = 5 if n == 21 else 10
    r["participants"] = [
        participant("p_substrate_1", "mol_1", "reactant", evidence, amount=0.2),
        participant("p_hydrazone_2a", "mol_2a", "reactant", evidence, amount=0.4),
        participant("p_catalyst", cat_mid, "catalyst", evidence, amount=catalyst_amount,
                    mol_percent=catalyst_mol, reference="mol_1"),
        participant("p_base", base_mid, "reagent", evidence, amount=0.5,
                    equivalents=2.5, reference="mol_1"),
        participant("p_solvent", solv_mid, "solvent", evidence, volume=2.0),
        participant("p_product_3a", "mol_3a", "product", evidence, product_role="target product"),
    ]
    if n == 1:
        r["participants"].append(participant("p_side_product_3aprime", "mol_3aprime", "product", evidence,
                                              product_role="speculated side product"))
    r["reaction_scale"] = f("1 (0.2 mmol)", {"reference_molecule_id": "mol_1", "amount": {"number": 0.2, "unit": "mmol"}}, evidence)
    r["solvent_conditions"]["composition"] = f(solvent, [{"molecule_id": solv_mid}], evidence)
    r["solvent_conditions"]["total_volume"] = f("2.0 mL", {"number": 2.0, "unit": "mL"}, evidence)
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = ["p_substrate_1", "p_hydrazone_2a", "p_catalyst", "p_base", "p_solvent"]
    step["operation"] = f(f"Table 1 entry {n}", "combine and heat under the tabulated conditions", evidence)
    step["temperature"] = f(f"{temp} °C", {"number": temp, "unit": "°C"}, evidence)
    step["time"] = f(f"{hours:g} h", {"number": hours, "unit": "h"}, evidence)
    if n == 1:
        step["atmosphere"] = f("under nitrogen atmosphere", "nitrogen", evidence)
    if n == 19:
        r["special_conditions"]["reactor"] = f("[c] The threaded tube reaction.", "threaded tube", evidence)
    r["steps"] = [step]
    r["procedure"]["original_text"] = f(
        {"table_row": n, "catalyst": catalyst, "base": base, "solvent": solvent,
         "temperature_C": temp, "time_h": hours,
         "charges": "1 (0.2 mmol), 2a (0.4 mmol), [Pd] (0.02 mmol unless footnote e), base (0.5 mmol), solvent (2.0 mL)"},
        "Table 1 row and footnote conditions", evidence)
    out = st.new_outcome("yield_3a", "yield", ["mol_3a"], evidence)
    out["product_role"] = f("target product 3a", "target product", evidence)
    out["value"] = f(f"{yld}%", {"number": yld, "unit": "%"}, evidence)
    out["yield_type"] = f("Yield (%) [b]", "isolated yield", evidence)
    out["measurement_basis"] = f("[b] Isolated.", "isolated", evidence)
    out["scope"] = f(f"Table 1 entry {n}", {"table": "Table 1", "entry": n}, evidence)
    r["outcomes"] = [out]
    if n == 1:
        side = st.new_outcome("gcms_side_product_3aprime", "other", ["mol_3aprime"], evidence)
        side["product_role"] = f("side-product 3a'", "speculated side product", evidence)
        side["value"] = f("we speculate that the reaction generated the side-product 3a' by GC-MS data",
                          "speculated present", evidence)
        side["measurement_basis"] = f("GC-MS data", "GC-MS", evidence)
        side["scope"] = f("Entry 1 discussion", {"entry": 1}, evidence)
        r["outcomes"].append(side)
        sel = st.new_outcome("e_configuration_3a", "selectivity", ["mol_3a"], evidence)
        sel["product_role"] = f("target product 3a", "target product", evidence)
        sel["value"] = f("excellent E-configuration", "E; ratio not reported", evidence)
        sel["measurement_basis"] = f("1H NMR; coupling constant of the proton on C=C bond is 16 Hz",
                                     {"method": "1H NMR", "J_Hz": 16}, evidence)
        sel["scope"] = f("configuration of the newly formed C=C bond", {"metric": "E/Z", "ratio": None}, evidence)
        r["outcomes"].append(sel)
    if n == 2:
        r["notes"].append(f("there were still reductive products present", "reductive products present; identities and amounts unreported", evidence))
    if n == 5:
        r["notes"].append(f("MeTHF as solvent could inhibit the formation of byproducts", "byproduct formation inhibited; identities and amounts unreported", evidence))
    if n == 20:
        r["notes"].append(f("[d] When other bases were used, the yield was reduced.", "footnote attached to the CH3OLi entry", evidence))
    if n == 21:
        r["notes"].append(f("[e] [Pd] (0.01 mmol).", "catalyst charge reduced to 0.01 mmol (5 mol% relative to 1)", evidence))
    classified(r, evidence, "The exact atom identities are unresolved; labels express the two source-depicted net C-C formations only.")
    r["review"]["unresolved"] = [
        "No associated SI was supplied; addition order, workup, purification details, and product mass remain unreported in the supplied main text.",
        "Candidate structures and partial bond classification require expert review; eliminated atoms and hydrogens are not completely mapped.",
    ]
    if n == 19:
        r["review"]["unresolved"].append("Threaded-tube pressure is not reported.")
    return r


def scope_reaction(label, yld, hyd_mid, prod_mid, conflict=False):
    evidence = [E3]
    r = base_reaction(f"table2_{label}", evidence, "specific_experiment", ["substrate_scope"],
        {"page": 3, "table": "Table 2", "product": label},
        "Pd-catalyzed cascade metallo-ene cyclization/metallo-carbene coupling",
        "cascade metallo-ene cyclization/metallo-carbene coupling")
    r["general_scheme_ids"] = [f"{DOC}:table2_generic_scope"]
    r["participants"] = [
        participant("p_substrate_1", "mol_1", "reactant", evidence, amount=0.2),
        participant(f"p_hydrazone_{label[1:]}", hyd_mid, "reactant", evidence, amount=0.4),
        participant("p_catalyst", "mol_pdcl2_dppf", "catalyst", evidence, amount=0.02, mol_percent=10, reference="mol_1"),
        participant("p_base", "mol_tbuoli", "reagent", evidence, amount=0.5, equivalents=2.5, reference="mol_1"),
        participant("p_solvent", "mol_2methf", "solvent", evidence, volume=2.0),
        participant(f"p_product_{label}", prod_mid, "product", evidence, product_role="scope product"),
    ]
    r["reaction_scale"] = f("1 (0.2 mmol)", {"reference_molecule_id": "mol_1", "amount": {"number": 0.2, "unit": "mmol"}}, evidence)
    r["solvent_conditions"]["composition"] = f("2-MeTHF", [{"molecule_id": "mol_2methf"}], evidence)
    r["solvent_conditions"]["total_volume"] = f("2-MeTHF (2.0 mL)", {"number": 2.0, "unit": "mL"}, evidence)
    step = st.new_step("reaction", evidence)
    step["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    step["operation"] = f("Table 2 condition", "combine and heat under the Table 2 conditions", evidence)
    step["temperature"] = f("90 °C", {"number": 90, "unit": "°C"}, evidence)
    step["time"] = f("3 h", {"number": 3, "unit": "h"}, evidence)
    r["steps"] = [step]
    r["procedure"]["original_text"] = f(
        "[a] The condition: 1 (0.2 mmol), 2 (0.4 mmol), PdCl2(dppf) (0.02 mmol), tBuOLi (0.5 mmol), 2-MeTHF (2.0 mL); scheme arrow: 90 °C, 3 h.",
        "Table 2 general condition", evidence)
    out = st.new_outcome(f"yield_{label}", "yield", [prod_mid], evidence)
    out["product_role"] = f(f"product {label}", "scope product", evidence)
    if conflict:
        a = f("3j, 80% (Table 2 drawing)", {"number": 80, "unit": "%"}, evidence)
        b = f("The desired product 3j ... was isolated in 75% (prose)", {"number": 75, "unit": "%"}, [T3])
        out["value"] = f(["3j, 80%", "isolated in 75%"], None, [E3, T3], status="conflict", alternatives=[a, b])
        out["yield_type"] = f("Table value versus prose isolated value", None, [E3, T3], status="conflict",
                              alternatives=[f("Table 2 value", "yield; basis not printed", [E3]),
                                            f("isolated in 75%", "isolated yield", [T3])])
        out["measurement_basis"] = f("prose: isolated in 75%; Table 2 basis not printed", "isolated stated only for the 75% prose alternative", [E3, T3])
        r["notes"].append(f("Source conflict: Table 2 prints 3j, 80%, while the adjacent prose says product 3j was isolated in 75%.",
                            "Both values retained without choosing.", [E3, T3]))
        r["review"]["unresolved"].append("Resolve the 3j 80% table versus 75% isolated prose conflict from authoritative source/SI; no value is selected here.")
    else:
        out["value"] = f(f"{label}, {yld}%", {"number": yld, "unit": "%"}, evidence)
        out["yield_type"] = missing(evidence)
        out["measurement_basis"] = missing(evidence)
    out["scope"] = f(f"Table 2 product {label}", {"table": "Table 2", "product": label}, evidence)
    r["outcomes"] = [out]
    classified(r, evidence, "Scope-member atom mapping is not asserted; substituent identity and the two depicted net C-C formations are retained.")
    r["review"]["unresolved"] += [
        "Table 2 does not print a yield measurement basis for this row except that the prose calls the conflicting 75% value for 3j isolated.",
        "Atmosphere, pressure, addition order, workup, purification, and product mass are unreported in the supplied main text.",
        "Candidate structures include source-derived connectivity and drawn alkene geometry but require expert structure review; the possible ring stereocenter is unassigned.",
    ]
    return r


def overview_reaction(rid, record_type, evidence, name, participants, location):
    r = base_reaction(rid, evidence, record_type, [], location, name)
    r["participants"] = participants
    r["classification"]["coverage"] = f(None, "unresolved", evidence, method="inference",
        derivation="The overview lacks sufficiently specific atom mapping for a bond-difference assertion.")
    r["review"]["unresolved"] = ["Overview/reference arrow only; experimental quantities, procedure, outcomes, and exact structures are not supplied here."]
    return r


dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
bundle = st.empty(dataset["documents"])

# Core substrate and reagents.
molecules = []
molecules.append(mol("mol_1", "1", [E2, E3], name="allenamide 1",
    smiles="C=C=CN(CC=CCOC(C)=O)S(=O)(=O)c1ccc(C)cc1",
    chirality="achiral", geometry="unspecified", representation="partial",
    stereo_note="The allene has a terminal CH2 and is not assigned an axial configuration; the tether alkene geometry is not explicitly labelled, so it remains unspecified."))
molecules.append(mol("mol_2a", "2a", [E1, E2, E3], name="benzaldehyde N-tosylhydrazone",
    smiles=ar_hydrazone_smiles("c2ccccc2"), chirality="achiral", geometry="unspecified", representation="partial",
    stereo_note="No stereogenic center is present in the candidate graph; hydrazone C=N geometry is not explicitly named and remains unspecified."))
molecules.append(make_product("mol_3a", "3a", "c1ccccc1", [E1, E2, E3]))
molecules[-1]["stereochemistry"]["assignment"] = f(
    "excellent E-configuration; 1H NMR coupling constant of the proton on C=C bond is 16 Hz",
    {"alkene": "E", "method": "1H NMR", "J_Hz": 16}, [T1, E1], method="original_extraction")
molecules.append(mol("mol_3aprime", "3a'", [E1, E2], name="side-product 3a'",
    chirality="unresolved", geometry="unresolved", representation="unresolved"))

hydrazone_defs = {
    "b": ("4-methylphenyl", "c2ccc(C)cc2"),
    "c": ("4-tert-butylphenyl", "c2ccc(C(C)(C)C)cc2"),
    "d": ("4-ethylphenyl", "c2ccc(CC)cc2"),
    "e": ("4-isopropylphenyl", "c2ccc(C(C)C)cc2"),
    "f": ("4-fluorophenyl", "c2ccc(F)cc2"),
    "g": ("4-chlorophenyl", "c2ccc(Cl)cc2"),
    "h": ("4-(trifluoromethyl)phenyl", "c2ccc(C(F)(F)F)cc2"),
    "i": ("furan-2-yl", "c2ccco2"),
    "k": ("3-methylphenyl", "c2cccc(C)c2"),
    "l": ("2-methylphenyl", "c2ccccc2C"),
    "m": ("3-methoxyphenyl", "c2cccc(OC)c2"),
    "n": ("2-methoxyphenyl", "c2ccccc2OC"),
    "o": ("3-chlorophenyl", "c2cccc(Cl)c2"),
    "p": ("2-chlorophenyl", "c2ccccc2Cl"),
    "q": ("[1,1'-biphenyl]-4-yl", "c2ccc(-c3ccccc3)cc2"),
    "r": ("3,5-dimethylphenyl", "c2cc(C)cc(C)c2"),
    "s": ("3,5-dimethoxyphenyl", "c2cc(OC)cc(OC)c2"),
    "t": ("3,5-difluorophenyl", "c2cc(F)cc(F)c2"),
}
for suffix, (name, ar) in hydrazone_defs.items():
    molecules.append(mol(f"mol_2{suffix}", f"2{suffix}", [E3], name=f"{name} aldehyde N-tosylhydrazone",
        smiles=ar_hydrazone_smiles(ar), chirality="achiral", geometry="unspecified", representation="partial",
        stereo_note="No stereogenic center is present in the candidate graph; hydrazone C=N geometry is not explicitly named and remains unspecified."))
    molecules.append(make_product(f"mol_3{suffix}", f"3{suffix}", ar, [E3]))

molecules.append(mol("mol_2j", "2j", [E3], name="benzophenone N-tosylhydrazone",
    smiles="Cc1ccc(S(=O)(=O)NN=C(c2ccccc2)c2ccccc2)cc1", chirality="achiral",
    geometry="not_applicable", representation="not_applicable",
    stereo_note="The hydrazone carbon bears two identical phenyl substituents; no C=N E/Z descriptor applies to this candidate graph."))
molecules.append(make_product("mol_3j", "3j", "", [E3], diphenyl=True))

# Reagents, catalysts, ligands/solvents as source-named entities.
for mid, label, smi in [
    ("mol_tbuoli", "tBuOLi", "[Li+].[O-]C(C)(C)C"),
    ("mol_ch3oli", "CH3OLi", "[Li+].[O-]C"),
    ("mol_toluene", "toluene", "Cc1ccccc1"),
    ("mol_iproh", "iPrOH", "CC(O)C"),
]:
    molecules.append(simple_agent(mid, label, [E1, E2, E3], smi))
molecules.append(mol("mol_2methf", "2-MeTHF", [E1, E2, E3], name="2-methyltetrahydrofuran",
    smiles="CC1CCCO1", chirality="unspecified", geometry="not_applicable", representation="partial",
    stereo_note="2-Methyltetrahydrofuran has a potential stereocenter at C2, but the source gives no enantiomeric identity or mixture composition; no racemic assumption is made."))
molecules.append(simple_agent("mol_dioxane", "dioxane", [E1, E2]))

catalyst_ids = {
    "PdCl2(PPh3)2": "mol_pdcl2_pph3_2", "PdCl2(dppf)": "mol_pdcl2_dppf",
    "Pd2(dba)3/L1": "mol_pd2dba3_l1", "Pd(OAc)2/L2": "mol_pdoac2_l2",
    "Pd(OAc)2/L3": "mol_pdoac2_l3", "Pd(OAc)2/L4": "mol_pdoac2_l4",
    "Pd(OAc)2/L5": "mol_pdoac2_l5", "Pd(OAc)2/L6": "mol_pdoac2_l6",
    "Pd(OAc)2/L7": "mol_pdoac2_l7", "Pd(OAc)2/dppm": "mol_pdoac2_dppm",
    "Pd(OAc)2/dppe": "mol_pdoac2_dppe", "Pd(OAc)2/dppp": "mol_pdoac2_dppp",
    "Pd(OAc)2/dppf": "mol_pdoac2_dppf", "Pd(OAc)2/dppb": "mol_pdoac2_dppb",
}
for label, mid in catalyst_ids.items():
    molecules.append(simple_agent(mid, label, [E2]))

solvent_ids = {"dioxane": "mol_dioxane", "toluene": "mol_toluene", "iPrOH": "mol_iproh", "MeTHF": "mol_2methf"}

# Generic overview/template identities.
molecules += [
    template_mol("mol_generic_A", "A", [E1], "allenamide A bearing a tethered allyl acetate and N-protecting group PG", [("PG", "source-defined nitrogen protecting group")]),
    template_mol("mol_generic_I", "I", [E1], "2,3-dihydropyrrole product bearing generic R from Suzuki coupling", [("R", "generic Suzuki-derived substituent"), ("PG", "nitrogen protecting group")]),
    template_mol("mol_generic_II", "II", [E1], "2,3-dihydropyrrole product bearing generic alkynyl R", [("R", "generic alkynyl substituent"), ("PG", "nitrogen protecting group")]),
    template_mol("mol_generic_III", "III", [E1], "2,3-dihydropyrrole triene product bearing R and R'", [("R", "generic carbene-derived substituent"), ("R'", "generic carbene-derived substituent"), ("PG", "nitrogen protecting group")]),
    template_mol("mol_generic_suzuki_partner", "Suzuki coupling partner", [E1], "generic coupling partner supplying R", [("R", "generic organic substituent")]),
    template_mol("mol_generic_alkyne_partner", "Sonogashira coupling partner", [E1], "generic terminal alkyne supplying alkynyl R", [("R", "generic organic substituent")]),
    template_mol("mol_generic_carbene_partner", "carbene precursor R/R'", [E1], "generic diazo/carbene precursor", [("R", "generic substituent"), ("R'", "generic substituent")]),
    template_mol("mol_generic_hydrazone_2", "2", [E3], "N-tosyl hydrazone bearing R1 and H at the hydrazone carbon", [("R1", "Ar or alkyl, as printed in Table 2")]),
    template_mol("mol_generic_product_3", "3", [E3], "2,3-dihydropyrrole triene scope product bearing R1", [("R1", "Ar or alkyl, as printed in Table 2")]),
    simple_agent("mol_generic_pd0", "Pd(0)/PdLn", [E1, E3]),
]
bundle["molecules"] = molecules

optimization_rows = [
    (1, "PdCl2(PPh3)2", "tBuOLi", "dioxane", 90, 4, 56),
    (2, "PdCl2(dppf)", "tBuOLi", "dioxane", 90, 4, 73),
    (3, "PdCl2(dppf)", "tBuOLi", "toluene", 90, 4, 71),
    (4, "PdCl2(dppf)", "tBuOLi", "iPrOH", 90, 4, 26),
    (5, "PdCl2(dppf)", "tBuOLi", "MeTHF", 90, 4, 80),
    (6, "Pd2(dba)3/L1", "tBuOLi", "MeTHF", 90, 4, 46),
    (7, "Pd(OAc)2/L2", "tBuOLi", "MeTHF", 90, 3.5, 31),
    (8, "Pd(OAc)2/L3", "tBuOLi", "MeTHF", 90, 3, 36),
    (9, "Pd(OAc)2/L4", "tBuOLi", "MeTHF", 90, 2.5, 22),
    (10, "Pd(OAc)2/L5", "tBuOLi", "MeTHF", 90, 2.5, 46),
    (11, "Pd(OAc)2/L6", "tBuOLi", "MeTHF", 90, 3, 61),
    (12, "Pd(OAc)2/L7", "tBuOLi", "MeTHF", 90, 3, 51),
    (13, "Pd(OAc)2/dppm", "tBuOLi", "MeTHF", 90, 3, 16),
    (14, "Pd(OAc)2/dppe", "tBuOLi", "MeTHF", 90, 3, 44),
    (15, "Pd(OAc)2/dppp", "tBuOLi", "MeTHF", 90, 3.5, 69),
    (16, "Pd(OAc)2/dppf", "tBuOLi", "MeTHF", 90, 3, 73),
    (17, "Pd(OAc)2/dppb", "tBuOLi", "MeTHF", 90, 3.5, 67),
    (18, "PdCl2(dppf)", "tBuOLi", "MeTHF", 90, 3, 82),
    (19, "PdCl2(dppf)", "tBuOLi", "MeTHF", 110, 3, 70),
    (20, "PdCl2(dppf)", "CH3OLi", "MeTHF", 90, 3, 10),
    (21, "PdCl2(dppf)", "tBuOLi", "MeTHF", 90, 3, 51),
]

reactions = [optimization_reaction(row) for row in optimization_rows]
scope_yields = {"b": 81, "c": 79, "d": 77, "e": 79, "f": 84, "g": 83, "h": 79, "i": 77,
                "j": 80, "k": 73, "l": 76, "m": 74, "n": 71, "o": 75, "p": 73, "q": 91,
                "r": 82, "s": 75, "t": 72}
for suffix, yld in scope_yields.items():
    reactions.append(scope_reaction(f"3{suffix}", yld, f"mol_2{suffix}", f"mol_3{suffix}", conflict=(suffix == "j")))

# Scheme 1: two referenced prior-work arrows and one present-work generic overview.
reactions.append(overview_reaction("scheme1_prior_suzuki", "literature_example", [E1],
    "Prior Suzuki coupling capture of palladium intermediate B/B'",
    [participant("p_A", "mol_generic_A", "reactant", [E1]),
     participant("p_partner", "mol_generic_suzuki_partner", "reactant", [E1]),
     participant("p_I", "mol_generic_I", "product", [E1], product_role="literature overview product")],
    {"page": 1, "scheme": "Scheme 1", "branch": "Suzuki coupling (I)"}))
reactions.append(overview_reaction("scheme1_prior_sonogashira", "literature_example", [E1],
    "Prior Sonogashira coupling capture of palladium intermediate B/B'",
    [participant("p_A", "mol_generic_A", "reactant", [E1]),
     participant("p_partner", "mol_generic_alkyne_partner", "reactant", [E1]),
     participant("p_II", "mol_generic_II", "product", [E1], product_role="literature overview product")],
    {"page": 1, "scheme": "Scheme 1", "branch": "Sonogashira coupling (II)"}))
r = overview_reaction("scheme1_current_overview", "generic_scheme", [E1],
    "Pd-catalyzed metallo-ene/metallo-carbene coupling overview",
    [participant("p_A", "mol_generic_A", "reactant", [E1]),
     participant("p_carbene", "mol_generic_carbene_partner", "reactant", [E1]),
     participant("p_III", "mol_generic_III", "product", [E1], product_role="generic target scaffold")],
    {"page": 1, "scheme": "Scheme 1", "branch": "This work (III)"})
classified(r, [E1], "Generic R/R'/PG groups and exact atom correspondence remain unresolved.")
reactions.append(r)

r = overview_reaction("table2_generic_scope", "generic_scheme", [E3],
    "Generic substrate-scope scheme for N-tosyl hydrazones",
    [participant("p_1", "mol_1", "reactant", [E3], amount=0.2),
     participant("p_2", "mol_generic_hydrazone_2", "reactant", [E3], amount=0.4),
     participant("p_cat", "mol_pdcl2_dppf", "catalyst", [E3], amount=0.02, mol_percent=10, reference="mol_1"),
     participant("p_base", "mol_tbuoli", "reagent", [E3], amount=0.5, equivalents=2.5, reference="mol_1"),
     participant("p_solv", "mol_2methf", "solvent", [E3], volume=2.0),
     participant("p_3", "mol_generic_product_3", "product", [E3], product_role="generic scope product")],
    {"page": 3, "table": "Table 2", "region": "generic header"})
step = st.new_step("reaction", [E3]); step["participant_ids"] = ["p_1", "p_2", "p_cat", "p_base", "p_solv"]
step["operation"] = f("Table 2 generic arrow", "combine and heat", [E3])
step["temperature"] = f("90 °C", {"number": 90, "unit": "°C"}, [E3])
step["time"] = f("3 h", {"number": 3, "unit": "h"}, [E3])
r["steps"] = [step]
r["reaction_scale"] = f("1 (0.2 mmol)", {"reference_molecule_id": "mol_1", "amount": {"number": 0.2, "unit": "mmol"}}, [E3])
r["solvent_conditions"]["composition"] = f("2-MeTHF", [{"molecule_id": "mol_2methf"}], [E3])
r["solvent_conditions"]["total_volume"] = f("2-MeTHF (2.0 mL)", {"number": 2.0, "unit": "mL"}, [E3])
classified(r, [E3], "Generic R1 prevents complete atom mapping; only the depicted net bond classes are retained.")
reactions.append(r)

# Scheme 2 is a proposed pathway, not a performed mechanistic experiment.
r = base_reaction("scheme2_proposed_mechanism", [E3, T3, E4, T4], "proposed_mechanism", [],
    {"pages": [3, 4], "scheme": "Scheme 2", "continuation": "page 4 first sentence"},
    "Proposed mechanism for the Pd-carbene reaction system",
    "proposed oxidative addition / metallo-ene cyclization / palladium-carbene insertion / beta-H elimination pathway")
r["participants"] = [
    participant("p_1", "mol_1", "reactant", [E3]),
    participant("p_2a", "mol_2a", "reactant", [E3]),
    participant("p_base", "mol_tbuoli", "reagent", [E3]),
    participant("p_pd", "mol_generic_pd0", "catalyst", [E3]),
    participant("p_3a", "mol_3a", "product", [E3, E4], product_role="proposed pathway product"),
]
r["procedure"]["original_text"] = f(
    "Oxidative addition gives 1-III/1-IV; hydrazone gives diazo 2a-I; palladium carbene 2a-II, migratory insertion 2a-III, then beta-H elimination gives 3a.",
    "author-proposed pathway, not an experimental procedure", [T3, T4])
classified(r, [E3, E4], "Elementary organopalladium intermediates are schematic and not fully atom-mapped; classification refers only to the proposed overall product-forming transformation.")
r["review"]["unresolved"] = ["Proposed mechanism only; no kinetic, trapping, isotope-labeling, or other independent mechanistic experiment is reported on the supplied pages."]
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:figure1_bioactive_examples", "evidence_ids": [E2], "kind": "nonreaction_figure",
     "note": "Figure 1 shows representative bioactive molecules and is not a synthesis experiment in this paper."},
    {"context_id": f"{DOC}:figure2_ligand_angle_plot", "evidence_ids": [E3], "kind": "derived_plot",
     "note": "Figure 2 plots P-Pd-P bond angle against yields already represented by Table 1 entries; it is not a separate experiment."},
    {"context_id": f"{DOC}:entry1_gcms_observation", "evidence_ids": [E1, E2], "kind": "characterization_observation",
     "note": "The authors speculate that side-product 3a' was generated based on GC-MS; no amount or spectrum is printed."},
    {"context_id": f"{DOC}:si_statement", "evidence_ids": [E1, T1], "kind": "document_metadata",
     "note": "The paper states SI is available online, but no SI file was supplied or explicitly associated in this session."},
    {"context_id": f"{DOC}:pages4_5_nonexperimental", "evidence_ids": [E4, E5], "kind": "coverage_note",
     "note": "Page 4 contains the mechanism continuation, conclusion, acknowledgments, and references; page 5 contains references and manuscript dates, with no additional current-paper synthesis arrows or tables."},
]

out = ROOT / "runs/current/submitted.json"
out.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

# Save every graph diagnostic separately; this does not elevate identity status.
diagnostics = []
for m in molecules:
    smi = m["smiles"]["normalized_value"]
    if smi:
        diagnostics.append({"molecule_id": m["molecule_id"], "candidate_smiles": smi, "diagnostic": st.stereo_inspect(smi)})
(ROOT / "runs/current/raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

mentions = st.source_mentions(dataset, ["1", "2a", "3a", "3a'", "3j", "tBuOLi", "PdCl2(dppf)"], context=220)
(ROOT / "runs/current/raw_model/source_mentions.json").write_text(json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")

page2_records = [f"{DOC}:table1_entry_{n:02d}" for n in range(1, 22)]
page3_scope = [f"{DOC}:table2_3{x}" for x in "bcdefghijklmnopqrst"]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "title, abstract, and lower-right optimization prose", "finding": "Article title and initial Table 1 experiment: 1 + 2a, PdCl2(PPh3)2, 2.5 equiv tBuOLi, 90 °C in anhydrous dioxane under nitrogen; 3a 56%, E assignment by 1H NMR J=16 Hz, and GC-MS-based speculation of side-product 3a'.", "evidence_ids": [E1, T1]},
                {"region": "Scheme 1", "finding": "Two cited prior-work branches (Suzuki product I and Sonogashira product II) and one present-work generic metallo-carbene target III; overview arrows, not additional current-paper concrete experiments.", "evidence_ids": [E1]},
                {"region": "supporting-information statement", "finding": "SI is stated to be available online, but no SI file was supplied or explicitly associated; association remains unestablished.", "evidence_ids": [E1, T1]},
            ],
            "record_ids": [f"{DOC}:table1_entry_01", f"{DOC}:scheme1_prior_suzuki", f"{DOC}:scheme1_prior_sonogashira", f"{DOC}:scheme1_current_overview"],
            "unresolved": ["Exact side-product 3a' connectivity and amount are not established by the supplied page.", "No detailed workup/purification procedure is printed on this page."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Figure 1", "finding": "Representative bioactive structures only; not synthesis experiments in this paper.", "evidence_ids": [E2]},
                {"region": "Table 1 header, entries 1-21, ligand drawings, and footnotes", "finding": "Twenty-one independent condition-screening rows transcribed with catalyst system, base, solvent, temperature, time, isolated yield, common charges, threaded-tube footnote, alternative-base note, and reduced catalyst charge in entry 21.", "evidence_ids": [E2]},
                {"region": "right-column discussion", "finding": "Confirms solvent/ligand/temperature/time/base trends and identifies entry 18 as the standard conditions; byproduct comments are retained without inventing identities or amounts.", "evidence_ids": [E2, T2]},
            ],
            "record_ids": page2_records,
            "unresolved": ["Ligand charge relative to Pd is not printed for slash catalyst/ligand systems.", "Pressure for the threaded-tube entry is not reported.", "Workup, purification, and isolated product mass are not available without associated SI."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Figure 2", "finding": "Plot of P-Pd-P angle versus yields already represented by Table 1; not separate experiments.", "evidence_ids": [E3]},
                {"region": "Table 2 generic header and products 3a-3t", "finding": "Generic scope conditions and twenty pictured products reviewed. 3a/82% duplicates Table 1 entry 18 and is merged; 3b-3t are nineteen additional scope experiments. Table prints 3j 80% while prose states isolated 75%, retained as an unresolved conflict.", "evidence_ids": [E3, T3]},
                {"region": "Scheme 2 and mechanism prose", "finding": "Author-proposed Pd pathway (oxidative addition, metallo-ene cyclization, diazo formation, Pd-carbene insertion, beta-H elimination); classified as proposed mechanism, not an experimental mechanistic study.", "evidence_ids": [E3, T3]},
            ],
            "record_ids": [f"{DOC}:table1_entry_18", f"{DOC}:table2_generic_scope", *page3_scope, f"{DOC}:scheme2_proposed_mechanism"],
            "unresolved": ["Yield basis is not printed for Table 2 rows except that prose calls the conflicting 75% value for 3j isolated.", "Possible product ring stereocenters are not assigned; scope drawings provide alkene geometry candidates only.", "No independent mechanistic control, isotope, kinetic, or trapping experiment is printed."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "top-left continuation", "finding": "Completes the proposed Scheme 2 sentence: intermediate 2a-III gives 3a by beta-H elimination; no new experiment.", "evidence_ids": [E4, T4]},
                {"region": "conclusion, acknowledgments, conflict statement, keywords, and references 1-19", "finding": "No additional current-paper synthesis table, arrow, or concrete experiment; references are contextual literature, not extracted as this paper's experiments.", "evidence_ids": [E4]},
            ],
            "record_ids": [f"{DOC}:scheme2_proposed_mechanism"],
            "unresolved": ["Conclusion supplies no additional experimental detail."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "references 20-21 and manuscript dates", "finding": "Reference list and received/revised/accepted dates only; no current-paper reaction arrow, table, or synthesis experiment.", "evidence_ids": [E5]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
(ROOT / "runs/current/coverage/pages_0001_0005.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(ROOT / "runs/current/raw_model/page1_5_visual_inventory.json").write_text(json.dumps({
    "document_id": DOC,
    "inspected_original_page_images": [E1, E2, E3, E4, E5],
    "page_inventory": coverage["pages"],
    "limitations": "Visual inspection is machine-assisted and remains needs_review; attempted page coverage is not semantic accuracy verification."
}, ensure_ascii=False, indent=2) + "\n")

tasks_path = ROOT / "runs/current/tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task.get("document_id") == DOC:
        task["state"] = "attempted"
        task["reason"] = "Original page image and extracted text inspected; see coverage/pages_0001_0005.json. Attempted does not certify semantic accuracy."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

# Audit the merged compatibility records when present.  This remains advisory
# and deliberately does not rewrite missing facts.
if dataset.get("experiments"):
    audit = rt.dataset_audit(dataset)
    (ROOT / "runs/current/audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({"molecules": len(molecules), "reactions": len(reactions)}, ensure_ascii=False))
