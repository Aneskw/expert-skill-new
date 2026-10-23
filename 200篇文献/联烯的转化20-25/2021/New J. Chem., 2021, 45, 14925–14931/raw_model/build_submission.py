#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1–6.

Page 7 is used only as contextual evidence for product names and explicit
figure/experimental-paragraph conflicts; it is not included in the coverage
claim produced by this script.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_0414a2e90dea4adf"
DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model"


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


E1I, E2I, E3I, E4I, E5I, E6I = [ev(i, "image") for i in range(1, 7)]
E1T, E2T, E3T, E4T, E5T, E6T, E7T = [ev(i) for i in range(1, 8)]


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def q(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="structure_derivation",
             derivation=derivation)


def not_applicable(evidence=()):
    return f(evidence=evidence, method="not_applicable", status="not_applicable")


data = json.loads(DATASET.read_text())
bundle = st.empty(data["documents"])
molecules = {}
stereo_diagnostics = {}


def set_status(m, key, value, evidence, derivation, method="inference", raw=None):
    m["stereochemistry"][key] = f(raw, value, evidence, method=method,
                                    derivation=derivation)


def make_template(label, fixed, variable_definition, evidence,
                  chirality="unresolved", geometry="unresolved",
                  representation="not_encoded"):
    mid = f"{DOC}:mol:{label}"
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "generic_template"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    m["names"]["english"] = f(label, label, evidence)
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, evidence),
        "variables": [{
            "label": "R",
            "attachment_points": f("R attached at the marked carbon", ["marked carbon"], evidence),
            "definition": f(variable_definition, variable_definition, evidence),
            "allowed_values": f("alkyl or aryl as drawn/stated", ["alkyl", "aryl"], evidence),
            "constraints": f("No unique R identity in the generic scheme", "source-defined variable", evidence),
        }],
        "instance_mapping": [],
        "reason": f("Generic R-group drawing", "The source uses an unexpanded R variable.", evidence),
    }
    set_status(m, "chirality_status", chirality, evidence,
               "Status preserves what can be concluded from the generic source drawing.")
    set_status(m, "geometry_status", geometry, evidence,
               "Status preserves the source drawing without inventing an isomer.")
    set_status(m, "representation_status", representation, evidence,
               "A generic R-group cannot be completely encoded as one concrete SMILES.")
    m["identity_note"] = f("Generic source structure; R is not expanded.",
                           "Generic source structure; R is not expanded.", evidence)
    molecules[mid] = m
    return mid


def make_molecule(label, name, evidence, smiles=None, stereo="nonstereo",
                  identity_note=None):
    mid = f"{DOC}:mol:{label}"
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f(smiles, smiles, evidence)
        m["stereochemistry"]["isomeric_smiles"] = inferred(
            diag["canonical_isomeric_smiles"], evidence,
            "RDKit canonical isomeric SMILES generated from the source-derived candidate graph.")
        m["stereochemistry"]["identity_key"] = inferred(
            diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key generated from the candidate SMILES; not a sample identity.")
        if stereo == "terminal_allene":
            m["stereochemistry"]["axial"] = inferred(
                {"allene": "terminal CH2 end; no R_a/S_a assignment"}, evidence,
                "The source drawing is a monosubstituted terminal allene; one terminal carbon has two hydrogens.")
            set_status(m, "chirality_status", "achiral", evidence,
                       "The terminal CH2 has two identical substituents, so this allene axis is not stereogenic.",
                       method="structure_derivation")
            set_status(m, "geometry_status", "not_applicable", evidence,
                       "E/Z notation is not assigned to the terminal cumulene drawing.", method="structure_derivation")
            set_status(m, "representation_status", "partial", evidence,
                       "The graph is encoded, but RDKit reports limited/unspecified cumulene stereo diagnostics.",
                       method="structure_derivation")
        elif stereo == "z_alkene":
            m["stereochemistry"]["double_bonds"] = f(
                "(Z) in the source product name", [{"configuration": "Z"}], evidence)
            m["stereochemistry"]["assignment"] = f(
                "(Z) product name", {"alkene_geometry": "Z", "assignment_source": "product name"}, evidence)
            set_status(m, "chirality_status", "achiral", evidence,
                       "No tetrahedral or axial stereogenic element is present in the encoded candidate graph.",
                       method="structure_derivation")
            set_status(m, "geometry_status", "specified", evidence,
                       "The product name explicitly assigns Z geometry.", method="original_extraction",
                       raw="(Z)")
            set_status(m, "representation_status", "encoded", evidence,
                       "The candidate isomeric SMILES encodes the source-assigned Z alkene.",
                       method="structure_derivation")
        elif stereo == "terminal_alkene":
            set_status(m, "chirality_status", "achiral", evidence,
                       "No tetrahedral or axial stereogenic element is present in the candidate graph.",
                       method="structure_derivation")
            set_status(m, "geometry_status", "not_applicable", evidence,
                       "The alkene is terminal CH2=C and has no E/Z descriptor.", method="structure_derivation")
            set_status(m, "representation_status", "encoded", evidence,
                       "The complete candidate connectivity is encoded; source correctness remains unreviewed.",
                       method="structure_derivation")
        else:
            set_status(m, "chirality_status", "not_applicable", evidence,
                       "No source stereochemical role is assigned to this non-product reagent/solvent entry.",
                       method="structure_derivation")
            set_status(m, "geometry_status", "not_applicable", evidence,
                       "No source E/Z or ring cis/trans assignment applies to this entry.",
                       method="structure_derivation")
            set_status(m, "representation_status", "encoded", evidence,
                       "The non-stereochemical candidate graph is fully represented.",
                       method="structure_derivation")
    else:
        m["representation_type"] = "unresolved"
        set_status(m, "chirality_status", "unresolved", evidence,
                   "The complete structure/stereochemical inventory was not safely recoverable.")
        set_status(m, "geometry_status", "unresolved", evidence,
                   "The complete structure or geometric assignment was not safely recoverable.")
        set_status(m, "representation_status", "not_encoded", evidence,
                   "No complete candidate SMILES is asserted for this source identity.")
    note = identity_note or ("Source identity with a machine candidate SMILES; parsing does not establish correctness."
                             if smiles else "Source identity retained without a guessed complete structure.")
    m["identity_note"] = f(note, note, evidence)
    molecules[mid] = m
    return mid


# Generic current-work structures.
g1 = make_template("generic_allene_1", "R-CH=C=CH2", "R is the allene substituent", [E3I, E4I],
                   chirality="achiral", geometry="not_applicable", representation="partial")
g2 = make_template("generic_product_2", "CH2=C(Bpin)-CH2-R", "R is alkyl or aryl", [E3I, E4I],
                   chirality="achiral", geometry="not_applicable", representation="not_encoded")
g3 = make_template("generic_product_3", "(Z)-R-CH=C(Bpin)-CH3", "R is alkyl or aryl", [E4I],
                   chirality="achiral", geometry="specified", representation="not_encoded")
molecules[g3]["stereochemistry"]["double_bonds"] = f(
    "generic minor product is drawn/named as Z", [{"configuration": "Z"}], [E4I])
molecules[g3]["stereochemistry"]["assignment"] = f(
    "3 minor; Z products named in experimental text", {"alkene_geometry": "Z"}, [E4I, E5T, E6T])


BPIN = "B1OC(C)(C)C(C)(C)O1"
scope = {
    "a": ("phenylallene", "C=C=Cc1ccccc1",
          f"C=C({BPIN})Cc1ccccc1", f"C/C(=C\\c1ccccc1){BPIN}"),
    "b": ("p-tolylallene", "C=C=Cc1ccc(C)cc1",
          f"C=C({BPIN})Cc1ccc(C)cc1", f"C/C(=C\\c1ccc(C)cc1){BPIN}"),
    "c": ("p-tert-butylphenylallene", "C=C=Cc1ccc(C(C)(C)C)cc1",
          f"C=C({BPIN})Cc1ccc(C(C)(C)C)cc1", f"C/C(=C\\c1ccc(C(C)(C)C)cc1){BPIN}"),
    "d": ("p-chlorophenylallene", "C=C=Cc1ccc(Cl)cc1",
          f"C=C({BPIN})Cc1ccc(Cl)cc1", f"C/C(=C\\c1ccc(Cl)cc1){BPIN}"),
    "e": ("m-chlorophenylallene", "C=C=Cc1cccc(Cl)c1",
          f"C=C({BPIN})Cc1cccc(Cl)c1", None),
    "f": ("o-chlorophenylallene", "C=C=Cc1ccccc1Cl",
          f"C=C({BPIN})Cc1ccccc1Cl", f"C/C(=C\\c1ccccc1Cl){BPIN}"),
    "g": ("p-trifluoromethylphenylallene", "C=C=Cc1ccc(C(F)(F)F)cc1",
          f"C=C({BPIN})Cc1ccc(C(F)(F)F)cc1", f"C/C(=C\\c1ccc(C(F)(F)F)cc1){BPIN}"),
    "h": ("p-fluorophenylallene", "C=C=Cc1ccc(F)cc1",
          f"C=C({BPIN})Cc1ccc(F)cc1", f"C/C(=C\\c1ccc(F)cc1){BPIN}"),
    "i": ("p-methoxyphenylallene", "C=C=Cc1ccc(OC)cc1",
          f"C=C({BPIN})Cc1ccc(OC)cc1", f"C/C(=C\\c1ccc(OC)cc1){BPIN}"),
    "j": ("p-benzyloxyphenylallene", "C=C=Cc1ccc(OCc2ccccc2)cc1",
          f"C=C({BPIN})Cc1ccc(OCc2ccccc2)cc1", None),
    "k": ("p-allyloxyphenylallene", "C=C=Cc1ccc(OCC=C)cc1",
          f"C=C({BPIN})Cc1ccc(OCC=C)cc1", f"C/C(=C\\c1ccc(OCC=C)cc1){BPIN}"),
    "l": ("cyclohexylallene", "C=C=CC1CCCCC1",
          f"C=C({BPIN})CC1CCCCC1", f"C/C(=C\\C1CCCCC1){BPIN}"),
    "m": ("allene 1m (source structure)", "C=C=CCCOCc1ccccc1",
          f"C=C({BPIN})CCCOCc1ccccc1", f"C/C(=C\\CCOCc1ccccc1){BPIN}"),
    "n": ("allene 1n (source structure)", "C=C=CCOCc1ccccc1",
          f"C=C({BPIN})CCOCc1ccccc1", f"C/C(=C\\COCc1ccccc1){BPIN}"),
    "o": ("allene 1o (source structure)", "C=C=" + "C" * 8,
          f"C=C({BPIN})" + "C" * 8, f"C/C(=C\\{'C' * 7}){BPIN}"),
    "p": ("allene 1p (source structure)", "C=C=" + "C" * 9,
          f"C=C({BPIN})" + "C" * 9, f"C/C(=C\\{'C' * 8}){BPIN}"),
    "q": ("allene 1q (source structure)", "C=C=" + "C" * 13,
          f"C=C({BPIN})" + "C" * 13, f"C/C(=C\\{'C' * 12}){BPIN}"),
}

product_names = {
    "a": ("4,4,5,5-tetramethyl-2-(3-phenylprop-1-en-2-yl)-1,3,2-dioxaborolane",
          "(Z)-4,4,5,5-tetramethyl-2-(1-phenylprop-1-en-2-yl)-1,3,2-dioxaborolane"),
    "b": ("4,4,5,5-tetramethyl-2-(3-(p-tolyl)prop-1-en-2-yl)-1,3,2-dioxaborolane",
          "(Z)-4,4,5,5-tetramethyl-2-(1-(p-tolyl)prop-1-en-2-yl)-1,3,2-dioxaborolane"),
    "c": ("2-(3-(4-(tert-butyl)phenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-(4-(tert-butyl)phenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "d": ("2-(3-(4-chlorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-(4-chlorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "e": ("2-(3-(3-chlorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "3e; source name says (Z)-2-(1-(4-chlorophenyl)prop-1-en-2-yl)-..."),
    "f": ("2-(3-(2-chlorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-(2-chlorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "g": ("4,4,5,5-tetramethyl-2-(3-(4-(trifluoromethyl)phenyl)prop-1-en-2-yl)-1,3,2-dioxaborolane",
          "(Z)-4,4,5,5-tetramethyl-2-(1-(4-(trifluoromethyl)phenyl)prop-1-en-2-yl)-1,3,2-dioxaborolane"),
    "h": ("2-(3-(4-fluorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-(4-fluorophenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "i": ("2-(3-(4-methoxyphenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-(4-methoxyphenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "j": ("2-(3-(4-(benzyloxy)phenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "3j; source name says (Z)-2-(1-(3-(benzyloxy)phenyl)prop-1-en-2-yl)-..."),
    "k": ("2-(3-(4-(allyloxy)phenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-(4-(allyloxy)phenyl)prop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "l": ("2-(3-cyclohexylprop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(1-cyclohexylprop-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "m": ("2-(5-(benzyloxy)pent-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(5-(benzyloxy)pent-2-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "n": ("2-(4-(benzyloxy)but-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(4-(benzyloxy)but-2-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "o": ("2-(dec-1-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane",
          "(Z)-2-(dec-2-en-2-yl)-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
    "p": ("4,4,5,5-tetramethyl-2-(undec-1-en-2-yl)-1,3,2-dioxaborolane",
          "(Z)-4,4,5,5-tetramethyl-2-(undec-2-en-2-yl)-1,3,2-dioxaborolane"),
    "q": ("4,4,5,5-tetramethyl-2-(pentadec-1-en-2-yl)-1,3,2-dioxaborolane",
          "(Z)-4,4,5,5-tetramethyl-2-(pentadec-2-en-2-yl)-1,3,2-dioxaborolane"),
}

scope_mids = {}
for letter, (sname, ssmi, p2smi, p3smi) in scope.items():
    sev = [E3I, E4I] if letter == "a" else [E4I]
    if letter in "ab": pev = [E4I, E5T]
    elif letter in "cdefghijkl": pev = [E4I, E6T]
    else: pev = [E4I, E7T]
    m1 = make_molecule(f"1{letter}", sname, sev, ssmi,
                       "terminal_allene" if ssmi else "unresolved")
    note3 = None
    if letter == "e":
        note3 = "Identity conflict retained: Fig. 2 implies the m-chloro series, while the experimental product name for 3e says 4-chlorophenyl."
    if letter == "j":
        note3 = "Identity conflict retained: Fig. 2 shows para-benzyloxy substitution, while the experimental product name for 3j says 3-(benzyloxy)phenyl."
    m2 = make_molecule(f"2{letter}", product_names[letter][0], pev, p2smi,
                       "terminal_alkene" if p2smi else "unresolved")
    m3 = make_molecule(f"3{letter}", product_names[letter][1], pev, p3smi,
                       "z_alkene" if p3smi else "unresolved", note3)
    if p3smi is None:
        molecules[m3]["stereochemistry"]["double_bonds"] = f(
            "(Z) in the product name", [{"configuration": "Z"}], pev)
        molecules[m3]["stereochemistry"]["assignment"] = f(
            "(Z) in the product name", {"alkene_geometry": "Z"}, pev)
        set_status(molecules[m3], "geometry_status", "specified", pev,
                   "The source product name explicitly assigns Z geometry.",
                   method="original_extraction", raw="(Z)")
        set_status(molecules[m3], "chirality_status", "unresolved", pev,
                   "Complete graph unavailable; no additional stereochemical claim is made.")
        set_status(molecules[m3], "representation_status", "not_encoded", pev,
                   "The source Z assignment is retained textually without a guessed full graph.")
    scope_mids[letter] = (m1, m2, m3)


# Common reagents, catalysts, solvents, and workup/purification agents.
common_ev = [E3I, E4I, E5T]
reagents = {}
def reg(label, name, smiles=None, evidence=common_ev, note=None):
    reagents[label] = make_molecule(label, name, evidence, smiles, "nonstereo", note)
    return reagents[label]

reg("B2pin2", "bis(pinacolato)diboron", "CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C")
reg("BnNH2", "benzylamine", "NCc1ccccc1")
reg("PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1")
reg("water", "water", "O")
reg("CuSO4", "copper(II) sulfate", "[Cu+2].[O-]S(=O)(=O)[O-]")
reg("toluene", "toluene", "Cc1ccccc1", [E3I, E4I])
reg("chloroform", "chloroform", "ClC(Cl)Cl", [E5T])
reg("ethyl_acetate", "ethyl acetate", "CCOC(C)=O", [E5T])
reg("silica_gel", "silica gel", None, [E5T])
reg("4_picoline", "4-picoline", "Cc1ccncc1", [E3I])
reg("proton_sponge", "proton sponge", None, [E3I])
reg("DMAP", "4-dimethylaminopyridine (DMAP)", "CN(C)c1ccncc1", [E3I])
reg("NaOAc", "sodium acetate", "CC(=O)[O-].[Na+]", [E3I])
reg("PCy3", "tricyclohexylphosphine", "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1", [E3I])
reg("DPEPhos", "DPEPhos", None, [E3I])
reg("SPhos", "SPhos", None, [E3I])
reg("Cu_acac2", "copper(II) acetylacetonate, Cu(acac)2", None, [E3I])
reg("Cu_BF4_2", "Cu(BF4)2 (as reported)", None, [E3I])
reg("Cu_pyrazinecarboxylate", "Cu(II) 2-pyrazinecarboxylate, C10H6CuN4O4", None, [E3I])
reg("Cu_OH2", "copper(II) hydroxide", None, [E3I])
reg("TPGS_750_M", "TPGS-750-M surfactant", None, [E3I, E4I],
    "Named surfactant retained without guessing a unique molecular composition.")
reg("CuCl", "copper(I) chloride", "[Cu]Cl", [E3I])
reg("NaOtBu", "sodium tert-butoxide", "CC(C)(C)[O-].[Na+]", [E3I])
reg("methanol", "methanol", "CO", [E3I])
reg("THF", "tetrahydrofuran", "C1CCOC1", [E3I])
reg("diethyl_ether", "diethyl ether", "CCOCC", [E3I])
reg("HBpin", "pinacolborane", "[BH]1OC(C)(C)C(C)(C)O1", [E3I])
reg("dioxane", "dioxane (isomer not specified in Fig. 1c label)", None, [E3I])
reg("tBuOH", "tert-butanol", "CC(C)(C)O", [E3I])
reg("P_ligand", "P-ligand (Ma 2012)", None, [E3I])
reg("NHC_ligand", "NHC ligand (Hoveyda 2013)", None, [E3I])
reg("CF3Ar_Xan", "CF3Ar-Xan ligand", None, [E3I])
reg("chiral_NHC_Cu", "chiral NHC-Cu complex", None, [E3I])

# Named mixtures/templates that must not be forced into one concrete molecule.
hexanes = make_template("hexanes", "mixture of hexane isomers", "mixture composition", [E5T],
                        chirality="unresolved", geometry="not_applicable", representation="not_encoded")
molecules[hexanes]["generic_template"]["variables"] = []


def charge(p, key, raw, normalized, evidence, reference=None, method="original_extraction"):
    p["charges"][key] = f(raw, normalized, evidence, method=method)
    if reference is not None:
        p["equivalence_reference"] = f(
            f"relative to {reference}", reference, evidence, method=method)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def set_classification(r, evidence, observed=True):
    c = r["classification"]
    c["mapped_reaction_smiles"] = f()
    if observed:
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "B"}],
             "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
             "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
             "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ]
        labels = [
            {"category_number": 11, "code": "C-B", "supporting_change_ids": ["b1"]},
            {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]},
        ]
        c["bond_changes"] = [inferred(x, evidence,
            "Net protoboration change derived from the source substrate/product drawings; no complete atom mapping asserted.") for x in changes]
        c["labels"] = inferred(labels, evidence,
            "Taxonomy labels derived from the partial source-visible net changes.")
        c["coverage"] = inferred("partial", evidence,
            "C-B formation, net C-H gain, and C=C order decrease are retained; full atom mapping and all partners are unresolved.")
        c["limitations"] = [
            "No reviewed atom-mapped reaction SMILES is asserted.",
            "Hydrogen source is consistent with the aqueous protoboration description but elementary atom provenance is unresolved.",
            "Spectator/catalyst changes and complete leaving-group accounting are not mapped.",
        ]
    else:
        c["bond_changes"] = []
        c["labels"] = f(None, [], evidence, method="inference",
                          derivation="The source reports 0% product; no observed bond formation is asserted.")
        c["coverage"] = f(None, "not_applicable", evidence, method="inference",
                            derivation="No observed transformation in this control entry.")
        c["limitations"] = ["Intended products are retained only to target the reported 0% yield."]
    c["model_votes"] = []


def reaction_base(rid, evidence, record_type, categories, location, name, transformation):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(transformation, transformation, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


reactions = []

# Fig. 1 literature examples. Each source arrow is a separate literature-example record.
lit_hov_sub = make_molecule("lit_hoveyda_2013_substrate", "dienyl-substituted terminal allene shown in Fig. 1b", [E3I])
lit_hov_p1 = make_molecule("lit_hoveyda_2013_product_internal", "internal-protoboration product shown in Fig. 1b", [E3I])
lit_hov_p2 = make_molecule("lit_hoveyda_2013_product_terminal", "terminal-protoboration product shown in Fig. 1b", [E3I])
lit_ma_p1 = make_molecule("lit_ma_2012_product_internal", "Ma 2012 internal-protoboration product", [E3I])
lit_ma_p2 = make_molecule("lit_ma_2012_product_terminal", "Ma 2012 terminal-protoboration product", [E3I])
lit_semba_p = make_molecule("lit_semba_2013_product", "Semba 2013 vinylboronate product", [E3I])
lit_chiral_sub = make_template("lit_hoveyda_2014_disubstituted_allene", "R-C(R1)=C=CH2", "R/R1 substituents", [E3I],
                               chirality="unspecified", geometry="not_applicable", representation="not_encoded")
lit_chiral_p = make_template("lit_hoveyda_2014_chiral_product", "enantioenriched branched vinyl-Bpin product", "R/R1 substituents", [E3I],
                             chirality="specified", geometry="not_applicable", representation="not_encoded")
molecules[lit_chiral_p]["stereochemistry"]["assignment"] = f(
    "wedge/dash configuration shown; no R/S label", {"assignment": "source wedge/dash only"}, [E3I])
molecules[lit_chiral_sub]["generic_template"]["variables"].append({
    "label": "R1",
    "attachment_points": f("R1 attached at the substituted allene terminus", ["marked allene carbon"], [E3I]),
    "definition": f("R1 is the second source-defined allene substituent", "second allene substituent", [E3I]),
    "allowed_values": f("not enumerated in Fig. 1d", None, [E3I], status="unreadable"),
    "constraints": f("R and R1 are distinct generic substituent positions", "source-defined variable", [E3I]),
})

literature_specs = [
    ("fig1a_ma_2012", scope_mids["a"][0], [lit_ma_p1, lit_ma_p2],
     [("CuCl", "catalyst", 5), ("P_ligand", "ligand", 6), ("NaOtBu", "reagent", 0.2),
      ("methanol", "reagent", 2.0), ("B2pin2", "reagent", 1.2), ("diethyl_ether", "solvent", None)],
     "Et2O, rt, 1 h", "Ma, 2012"),
    ("fig1b_hoveyda_2013", lit_hov_sub, [lit_hov_p1, lit_hov_p2],
     [("CuCl", "catalyst", 0.05), ("NHC_ligand", "ligand", 0.05), ("NaOtBu", "reagent", 0.4),
      ("methanol", "reagent", 2.0), ("B2pin2", "reagent", 1.1), ("THF", "solvent", None)],
     "THF, rt, 1 h", "Hoveyda, 2013"),
    ("fig1c_semba_2013", scope_mids["l"][0], [lit_semba_p],
     [("CuCl", "catalyst", 2), ("CF3Ar_Xan", "ligand", 2), ("NaOtBu", "reagent", 12),
      ("HBpin", "reagent", 1.2), ("dioxane", "solvent", None)],
     "dioxane, 28 °C, 2 h", "Semba, 2013"),
    ("fig1d_hoveyda_2014", lit_chiral_sub, [lit_chiral_p],
     [("chiral_NHC_Cu", "catalyst", None), ("B2pin2", "reagent", None), ("tBuOH", "reagent", None)],
     "chiral NHC-Cu complex (1–5 mol%); B2pin2, tBuOH", "Hoveyda, 2014"),
]

for rid, sub, prods, agents, condition, citation in literature_specs:
    r = reaction_base(rid, [E3I], "literature_example", ["other"],
                      f"Fig. 1 ({rid[4]})", "allene protoboration", "protoboration")
    ps = [participant("substrate", sub, "reactant", [E3I])]
    for i, pmid in enumerate(prods, 1):
        ps.append(participant(f"product_{i}", pmid, "product", [E3I], product_role="alternative regioisomer"))
    for label, role, amount in agents:
        p = participant(label, reagents[label], role, [E3I])
        if amount is not None:
            unit = "mol%" if role in ("catalyst", "ligand") and amount >= 1 else "equiv"
            if rid == "fig1b_hoveyda_2013" and role in ("catalyst", "ligand"):
                unit = "equiv"
            charge(p, "mol_percent" if unit == "mol%" else "equivalents",
                   f"{amount} {unit}", q(amount, unit), [E3I], sub)
        ps.append(p)
    r["participants"] = ps
    r["solvent_conditions"]["composition"] = f(condition, condition.split(",")[0], [E3I])
    step = st.new_step("literature_arrow", [E3I])
    step["participant_ids"] = [p["participant_id"] for p in ps]
    step["operation"] = f(condition, condition, [E3I])
    if "28 °C" in condition:
        step["temperature"] = f("28 °C", q(28, "°C"), [E3I])
        step["time"] = f("2 h", q(2, "h"), [E3I])
    elif "rt, 1 h" in condition:
        step["temperature"] = f("rt", "rt", [E3I])
        step["time"] = f("1 h", q(1, "h"), [E3I])
    r["steps"] = [step]
    r["procedure"]["original_text"] = f(condition, condition, [E3I])
    r["notes"] = [f(citation, {"citation_label": citation, "context": "prior work"}, [E3I])]
    r["review"]["unresolved"] = ["Literature example only; yields and complete operational details are not given in Fig. 1.",
                                      "Some ligand and product structures are retained unresolved rather than guessed."]
    set_classification(r, [E3I], observed=True)
    reactions.append(r)


# One current-work generic scheme, with duplicate source evidence merged.
scheme_id = f"{DOC}:current_generic_scheme"
r = reaction_base("current_generic_scheme", [E3I, E4I], "generic_scheme", ["other"],
                  "Fig. 1e and Fig. 2 generic arrow", "protoboration of allenes",
                  "Cu(II)-catalyzed protoboration")
r["participants"] = [
    participant("generic_allene", g1, "reactant", [E3I, E4I]),
    participant("generic_major", g2, "product", [E3I, E4I], product_role="major regioisomer"),
    participant("generic_minor", g3, "product", [E4I], product_role="minor regioisomer"),
]
for label, role in [("CuSO4", "catalyst"), ("BnNH2", "additive"), ("PPh3", "ligand"),
                    ("B2pin2", "reagent"), ("water", "solvent")]:
    p = participant(label, reagents[label], role, [E3I, E4I])
    amounts = {"CuSO4": ("mol_percent", 1, "mol%"), "BnNH2": ("mol_percent", 5, "mol%"),
               "PPh3": ("mol_percent", 10, "mol%"), "B2pin2": ("equivalents", 1.1, "equiv")}
    if label in amounts:
        key, number, unit = amounts[label]
        charge(p, key, f"{number} {unit}", q(number, unit), [E3I, E4I], g1)
    r["participants"].append(p)
r["solvent_conditions"]["composition"] = f("H2O", ["water"], [E3I, E4I])
step = st.new_step("generic_reaction", [E3I, E4I])
step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
step["operation"] = f("H2O, rt, 3 h; open to air", "stir under the shown conditions", [E3I])
step["temperature"] = f("rt", "rt", [E3I, E4I])
step["time"] = f("3 h", q(3, "h"), [E3I, E4I])
step["atmosphere"] = f("Open to air", "open air", [E3I])
r["steps"] = [step]
r["procedure"]["original_text"] = f("CuSO4, NH2Bn, PPh3, B2pin2; H2O, rt, 3 h", "generic scheme only", [E3I, E4I])
r["review"]["unresolved"] = ["Generic scheme has no absolute scale or solvent volume.",
                                  "R is source-defined and not expanded into one concrete molecule."]
set_classification(r, [E3I, E4I], observed=True)
reactions.append(r)


# Table 1: 21 independent optimization rows.
rows = [
    (1, "4_picoline", "CuSO4", None, 10, (34, 66), "a"),
    (2, "proton_sponge", "CuSO4", None, 13, (40, 60), "a"),
    (3, "DMAP", "CuSO4", None, 9, (30, 70), "a"),
    (4, "NaOAc", "CuSO4", None, 16, (38, 62), "a"),
    (5, "BnNH2", "CuSO4", None, 24, (35, 65), "a"),
    (6, "BnNH2", "CuSO4", "PCy3", 65, (23, 77), "a"),
    (7, "BnNH2", "CuSO4", "PCy3", 49, (21, 79), "ac"),
    (8, "BnNH2", None, "PCy3", 0, None, "a"),
    (9, None, "CuSO4", "PCy3", 68, (33, 67), "ad"),
    (10, "BnNH2", "CuSO4", "DPEPhos", 52, (46, 54), ""),
    (11, "BnNH2", "CuSO4", "SPhos", 63, (81, 19), ""),
    (12, "BnNH2", "CuSO4", "PPh3", 78, (87, 13), ""),
    (13, "BnNH2", "CuSO4", "PCy3", 63, (31, 69), ""),
    (14, "BnNH2", "Cu_acac2", "PPh3", 60, (85, 15), ""),
    (15, "BnNH2", "Cu_BF4_2", "PPh3", 53, (77, 23), ""),
    (16, "BnNH2", "Cu_pyrazinecarboxylate", "PPh3", 32, (81, 19), ""),
    (17, "BnNH2", "Cu_OH2", "PPh3", 19, (76, 24), ""),
    (18, "BnNH2", "CuSO4", "PPh3", 50, (74, 26), "e"),
    (19, "BnNH2", "CuSO4", "PCy3", 57, (27, 73), "e"),
    (20, "BnNH2", "CuSO4", "PCy3", 52, (30, 70), "f"),
    (21, "BnNH2", "CuSO4", "PPh3", 36, (64, 36), "f"),
]

for entry, base, cu, ligand, yld, ratio, flags in rows:
    rid = f"table1_entry_{entry:02d}"
    r = reaction_base(rid, [E3I, E3T], "specific_experiment", ["optimization"],
                      f"Table 1, entry {entry}", "allene protoboration", "protoboration")
    r["general_scheme_ids"] = [scheme_id]
    ps = [participant("allene_1a", scope_mids["a"][0], "reactant", [E3I])]
    # Products are intended table targets even for entry 8; entry 8 explicitly reports 0%.
    ps += [participant("product_2a", scope_mids["a"][1], "product", [E3I], product_role="regioisomer 2a"),
           participant("product_3a", scope_mids["a"][2], "product", [E3I], product_role="regioisomer 3a")]
    b2 = participant("B2pin2", reagents["B2pin2"], "reagent", [E3I])
    b2eq = 1.5 if "a" in flags else 1.1
    charge(b2, "equivalents", f"{b2eq} equiv.", q(b2eq, "equiv"), [E3I], scope_mids["a"][0])
    ps.append(b2)
    water = participant("water", reagents["water"], "solvent", [E3I])
    ps.append(water)
    if base:
        p = participant("base", reagents[base], "additive", [E3I])
        charge(p, "mol_percent", "5 mol%", q(5, "mol%"), [E3I], scope_mids["a"][0])
        ps.append(p)
    if cu:
        p = participant("copper_source", reagents[cu], "catalyst", [E3I])
        charge(p, "mol_percent", "1 mol%", q(1, "mol%"), [E3I], scope_mids["a"][0])
        ps.append(p)
    if ligand:
        p = participant("ligand", reagents[ligand], "ligand", [E3I])
        charge(p, "mol_percent", "10 mol%", q(10, "mol%"), [E3I], scope_mids["a"][0])
        ps.append(p)
    if "c" in flags:
        p = participant("surfactant", reagents["TPGS_750_M"], "additive", [E3I])
        charge(p, "amount", "1% surfactant TPGS-750-M", q(1, "%"), [E3I])
        ps.append(p)
    if flags in ("e", "f"):
        ps.append(participant("toluene", reagents["toluene"], "solvent", [E3I]))
    r["participants"] = ps
    if flags == "e":
        r["solvent_conditions"]["composition"] = f("4:1 water:toluene", ["water", "toluene"], [E3I])
        r["solvent_conditions"]["mixture_ratio"] = f("4:1 water:toluene", {"water": 4, "toluene": 1}, [E3I])
    elif flags == "f":
        r["solvent_conditions"]["composition"] = f("4:1 toluene:water", ["toluene", "water"], [E3I])
        r["solvent_conditions"]["mixture_ratio"] = f("4:1 toluene:water", {"toluene": 4, "water": 1}, [E3I])
    else:
        r["solvent_conditions"]["composition"] = f("H2O", ["water"], [E3I])
        r["solvent_conditions"]["mixture_ratio"] = not_applicable([E3I])
    step = st.new_step("reaction", [E3I])
    step["participant_ids"] = [p["participant_id"] for p in ps]
    step["operation"] = f(f"Table 1 entry {entry}; stir under listed conditions", "stir", [E3I])
    step["temperature"] = f("rt", "rt", [E3I])
    hours = 6 if "d" in flags else 3
    step["time"] = f(f"{hours} h", q(hours, "h"), [E3I])
    step["atmosphere"] = f("Open to air (Fig. 1e current-work conditions)", "open air", [E3I],
                             method="inherited_general_condition")
    r["steps"] = [step]
    r["procedure"]["original_text"] = f(
        {"entry": entry, "base": base or "None", "Cu source": cu or "None", "ligand": ligand or "None"},
        "Table 1 row plus header/footnotes", [E3I])
    oy = st.new_outcome("combined_yield", "yield", [scope_mids["a"][1], scope_mids["a"][2]], [E3I])
    oy["product_role"] = f("combined 2a + 3a", "combined regioisomer yield", [E3I])
    oy["value"] = f(f"{yld}%", q(yld, "%"), [E3I])
    oy["yield_type"] = f("% yield", "NMR yield", [E3I])
    oy["measurement_basis"] = f(
        "Yields were determined by 1H NMR analysis of the crude reaction mixture after extraction.",
        "1H NMR of crude after extraction", [E3I])
    oy["scope"] = f("2a + 3a", {"products": [scope_mids["a"][1], scope_mids["a"][2]]}, [E3I])
    r["outcomes"] = [oy]
    if ratio:
        osel = st.new_outcome("regioisomer_ratio", "selectivity",
                              [scope_mids["a"][1], scope_mids["a"][2]], [E3I])
        osel["product_role"] = f("2a:3a", "regioisomer ratio", [E3I])
        osel["value"] = f(f"{ratio[0]}:{ratio[1]}",
                           {"ratio": list(ratio), "order": [scope_mids["a"][1], scope_mids["a"][2]]}, [E3I])
        osel["scope"] = f("(2a:3a)", {"metric": "2a:3a ratio", "order": ["2a", "3a"]}, [E3I])
        r["outcomes"].append(osel)
    if entry == 8:
        ono = st.new_outcome("no_product_observed", "no_reaction", [], [E3I])
        ono["value"] = f("0", "0% combined yield", [E3I])
        ono["scope"] = f("control without CuSO4", "copper-free control", [E3I])
        r["outcomes"].append(ono)
        r["review"]["unresolved"] = ["Detection limit for the 0% copper-free control is not reported."]
        set_classification(r, [E3I], observed=False)
    else:
        r["review"]["unresolved"] = ["Absolute substrate scale and solvent volume are unreported in Table 1.",
                                          "Measurement method for the 2a:3a ratio is not explicitly stated."]
        set_classification(r, [E3I], observed=True)
    reactions.append(r)


scope_values = {
    "a": (78, (87, 13)), "b": (67, (89, 11)), "c": (47, (90, 10)),
    "d": (87, (90, 10)), "e": (62, (89, 11)), "f": (87, (90, 10)),
    "g": (19, (93, 7)), "h": (19, (83, 17)), "i": (46, (67, 33)),
    "j": (37, (66, 34)), "k": (49, (75, 25)), "l": (84, (78, 22)),
    "m": (34, (93, 7)), "n": (53, None), "o": (None, None),
    "p": (None, None), "q": (75, (83, 17)),
}


def conflict_value(raws, alternatives):
    return f(raws, None, sorted({eid for _, _, eids in alternatives for eid in eids}),
             status="conflict",
             alternatives=[f(raw, norm, eids) for raw, norm, eids in alternatives])


for letter in scope:
    m1, m2, m3 = scope_mids[letter]
    specific_ev = [E4I]
    if letter in "ab": specific_ev += [E5T]
    elif letter in "cdefghijkl": specific_ev += [E6T]
    else: specific_ev += [E7T]
    evidence = list(dict.fromkeys(specific_ev + [E5T]))
    r = reaction_base(f"scope_2{letter}", evidence, "specific_experiment", ["substrate_scope"],
                      f"Fig. 2, products 2{letter}/3{letter}", "protoboration of allene",
                      "Cu(II)-catalyzed protoboration")
    r["general_scheme_ids"] = [scheme_id]
    ps = [participant(f"allene_1{letter}", m1, "reactant", evidence),
          participant(f"product_2{letter}", m2, "product", evidence, product_role="major regioisomer"),
          participant(f"product_3{letter}", m3, "product", evidence, product_role="minor Z regioisomer")]
    charge(ps[0], "equivalents", "1 equiv.", q(1, "equiv"), [E4I], m1)
    for label, role, key, number, unit in [
        ("B2pin2", "reagent", "equivalents", 1.1, "equiv"),
        ("BnNH2", "additive", "mol_percent", 5, "mol%"),
        ("PPh3", "ligand", "mol_percent", 10, "mol%"),
        ("CuSO4", "catalyst", "mol_percent", 1, "mol%"),
    ]:
        p = participant(label, reagents[label], role, evidence)
        charge(p, key, f"{number} {unit}", q(number, unit), [E4I, E5T], m1,
               method="inherited_general_condition")
        ps.append(p)
    ps.append(participant("water", reagents["water"], "solvent", evidence))
    if letter == "j":
        p = participant("TPGS_750_M", reagents["TPGS_750_M"], "additive", [E4I])
        charge(p, "amount", "1% TPGS-750-M", q(1, "%"), [E4I])
        ps.append(p)
    if letter in "opq":
        ps.append(participant("toluene", reagents["toluene"], "solvent", [E4I]))
    # General workup and purification materials are explicit in the page-5 procedure.
    ps += [participant("chloroform", reagents["chloroform"], "workup_agent", [E5T], stage="workup"),
           participant("silica_gel", reagents["silica_gel"], "purification_agent", [E5T], stage="purification"),
           participant("ethyl_acetate", reagents["ethyl_acetate"], "purification_agent", [E5T], stage="purification"),
           participant("hexanes", hexanes, "purification_agent", [E5T], stage="purification")]
    r["participants"] = ps
    if letter in "opq":
        r["solvent_conditions"]["composition"] = f("3:1 toluene:water", ["toluene", "water"], [E4I])
        r["solvent_conditions"]["mixture_ratio"] = f("3:1 toluene:water", {"toluene": 3, "water": 1}, [E4I])
    else:
        r["solvent_conditions"]["composition"] = f("H2O", ["water"], [E4I, E5T],
                                                      method="inherited_general_condition")
        r["solvent_conditions"]["mixture_ratio"] = f(
            "Milli-Q water and CuSO4 solution were added in equal volume amounts",
            {"Milli-Q water": 1, "aqueous CuSO4 solution": 1}, [E5T],
            method="inherited_general_condition")
    reaction_step = st.new_step("reaction", [E4I, E5T])
    reaction_step["participant_ids"] = [p["participant_id"] for p in ps if p["stage"] == "reaction"]
    reaction_step["operation"] = f(
        "dark reaction mixture was allowed to stir until completion (typically 3 h) as monitored by GC or TLC",
        "stir; monitor by GC or TLC", [E5T], method="inherited_general_condition")
    reaction_step["temperature"] = f("rt", "rt", [E4I])
    reaction_step["time"] = f("3 h", q(3, "h"), [E4I])
    reaction_step["atmosphere"] = f("All borylation reactions were carried out under atmospheric conditions.",
                                      "open air / atmospheric conditions", [E5T],
                                      method="inherited_general_condition")
    reaction_step["pressure"] = f("atmospheric conditions", "atmospheric pressure", [E5T],
                                    method="inherited_general_condition")
    workup = st.new_step("workup", [E5T], stage="workup")
    workup["participant_ids"] = ["chloroform"]
    workup["operation"] = f(
        "chloroform was added; water layer extracted 3× with chloroform; solvent removed in vacuo",
        ["add chloroform", "extract aqueous layer three times", "remove solvent in vacuo"],
        [E5T], method="inherited_general_condition")
    purification = st.new_step("purification", [E5T], stage="purification")
    purification["participant_ids"] = ["silica_gel", "ethyl_acetate", "hexanes"]
    purification["operation"] = f("flash chromatography (0–4% EtOAc in hexanes)",
                                    "flash chromatography; 0–4% EtOAc in hexanes", [E5T],
                                    method="inherited_general_condition")
    r["steps"] = [reaction_step, workup, purification]
    r["procedure"]["original_text"] = f(
        "General procedure for synthesis of vinylboronic acid derivatives",
        "Use page-5 general procedure with the Fig. 2 substrate-specific result.", [E5T],
        method="inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f(
        "General procedure for synthesis of vinylboronic acid derivatives",
        "page 5 general procedure", [E5T], method="inherited_general_condition")
    r["procedure"]["addition_order"] = f(
        "allene, B2pin2, benzylamine, PPh3; then Milli-Q water and CuSO4 solution",
        ["allene", "B2pin2", "benzylamine", "PPh3", "Milli-Q water", "CuSO4 solution"],
        [E5T], method="inherited_general_condition")

    oy = st.new_outcome("combined_isolated_yield", "yield", [m2, m3], specific_ev)
    oy["product_role"] = f(f"2{letter} + 3{letter}", "combined regioisomer mixture", specific_ev)
    yld, ratio = scope_values[letter]
    if letter == "o":
        oy["value"] = conflict_value(
            ["Fig. 2: 50%", "Experimental paragraph: 48%"],
            [("50%", q(50, "%"), [E4I]), ("48%", q(48, "%"), [E7T])])
    elif letter == "p":
        oy["value"] = conflict_value(
            ["Fig. 2: 48%", "Experimental paragraph: 50%"],
            [("48%", q(48, "%"), [E4I]), ("50%", q(50, "%"), [E7T])])
    else:
        oy["value"] = f(f"{yld}%", q(yld, "%"), specific_ev)
    oy["yield_type"] = f("isolated", "isolated yield", [E5T] if letter in "ab" else [E6T] if letter in "cdefghijkl" else [E7T])
    oy["measurement_basis"] = f("isolated", "isolated product mixture", [E5T] if letter in "ab" else [E6T] if letter in "cdefghijkl" else [E7T])
    oy["scope"] = f(f"combined 2{letter}/3{letter} product mixture",
                     {"products": [m2, m3], "start": m1}, specific_ev)
    r["outcomes"] = [oy]

    if letter == "n":
        ratio_fact = conflict_value(
            ["Fig. 2: 72:28", "Experimental paragraph: 72:18"],
            [("72:28", {"ratio": [72, 28], "order": [m2, m3]}, [E4I]),
             ("72:18", {"ratio": [72, 18], "order": [m2, m3]}, [E7T])])
    elif letter == "o":
        ratio_fact = conflict_value(
            ["Fig. 2: 84:16", "Experimental paragraph: 74:26"],
            [("84:16", {"ratio": [84, 16], "order": [m2, m3]}, [E4I]),
             ("74:26", {"ratio": [74, 26], "order": [m2, m3]}, [E7T])])
    elif letter == "p":
        ratio_fact = conflict_value(
            ["Fig. 2: 74:26", "Experimental paragraph: 84:16"],
            [("74:26", {"ratio": [74, 26], "order": [m2, m3]}, [E4I]),
             ("84:16", {"ratio": [84, 16], "order": [m2, m3]}, [E7T])])
    else:
        ratio_fact = f(f"{ratio[0]}:{ratio[1]}",
                       {"ratio": list(ratio), "order": [m2, m3]}, specific_ev)
    osel = st.new_outcome("regioisomer_ratio", "selectivity", [m2, m3], specific_ev)
    osel["product_role"] = f(f"2{letter}:3{letter}", "major:minor regioisomer ratio", specific_ev)
    osel["value"] = ratio_fact
    osel["scope"] = f(f"2{letter}:3{letter}",
                       {"metric": "regioisomer ratio", "order": [f"2{letter}", f"3{letter}"]}, specific_ev)
    # No ratio measurement method is explicitly stated.
    r["outcomes"].append(osel)
    if letter == "h":
        onmr = st.new_outcome("combined_nmr_yield", "yield", [m2, m3], [E6T])
        onmr["product_role"] = f("2h + 3h", "combined regioisomer mixture", [E6T])
        onmr["value"] = f("50% NMR yield", q(50, "%"), [E6T])
        onmr["yield_type"] = f("NMR yield", "NMR yield", [E6T])
        onmr["measurement_basis"] = f("NMR yield", "NMR; further method not stated in product paragraph", [E6T])
        onmr["scope"] = f("same 2h/3h experiment", {"products": [m2, m3]}, [E6T])
        r["outcomes"].append(onmr)
    unresolved = ["Absolute allene scale and solvent volume are not reported in the supplied main-text procedure.",
                  "Measurement method for the regioisomer ratio is unreported.",
                  "No SI file was supplied or associated; omitted details remain missing."]
    if letter in "mno pq".replace(" ", ""):
        unresolved.append("Complete candidate structures for the long-chain/benzyloxy substrate-product set were not asserted from this segment.")
    if letter == "e":
        unresolved.append("3e naming conflict: the major series/Fig. 2 is meta-chloro, but the experimental minor-product name says 4-chlorophenyl.")
    if letter == "j":
        unresolved.append("3j naming conflict: Fig. 2 shows para-benzyloxy substitution, but the experimental minor-product name says 3-(benzyloxy)phenyl.")
    if letter in "nop":
        unresolved.append("Fig. 2 and the page-7 experimental paragraph disagree on the reported ratio and/or yield; both values are retained as conflict alternatives.")
    r["review"]["unresolved"] = unresolved
    r["notes"] = [f("* indicates minor isomer in the experimental section.",
                      "Asterisked characterization peaks belong to the minor isomer.", [E5T])]
    if letter in "nop":
        r["notes"].append(f(
            "Source conflict retained between Fig. 2 and the later experimental paragraph; no value was selected.",
            "Source conflict retained without correction.", [E4I, E7T]))
    set_classification(r, [E4I], observed=True)
    reactions.append(r)


# Fig. 3 is a proposed mechanism, not an observed preparative experiment.
r = reaction_base("fig3_proposed_mechanism", [E4T, E5I, E5T], "proposed_mechanism",
                  ["mechanistic_study"], "Fig. 3", "proposed mechanism for Cu(II)-catalyzed protoboration",
                  "proposed catalytic cycle")
r["general_scheme_ids"] = [scheme_id]
r["participants"] = [
    participant("generic_allene", g1, "reactant", [E5I]),
    participant("B2pin2", reagents["B2pin2"], "reagent", [E5I]),
    participant("water", reagents["water"], "reagent", [E4T, E5I]),
    participant("generic_major", g2, "product", [E5I], product_role="proposed major product 12"),
    participant("generic_minor", g3, "product", [E5I], product_role="proposed alternative product 14"),
]
step = st.new_step("proposed_cycle", [E4T, E5I])
step["participant_ids"] = [p["participant_id"] for p in r["participants"]]
step["operation"] = f(
    "base-assisted boronate activation, transmetallation, allene insertion, and protonation via proposed intermediates 4–14",
    "author-proposed catalytic sequence", [E4T, E5I])
r["steps"] = [step]
r["procedure"]["original_text"] = f("Proposed mechanism for the Cu(II)-catalyzed protoboration.",
                                      "proposed mechanism, not an experimental procedure", [E5I])
r["notes"] = [f(
    "The authors state that the mechanism proceeds as described previously by Hoveyda and coworkers.",
    "Mechanistic interpretation is an author proposal, not an observed preparative record.", [E4T])]
r["review"]["unresolved"] = ["Intermediates 4–14 are proposed and were not encoded as independently observed molecules.",
                                  "No kinetic/isotopic experiment is reported in the inspected pages."]
set_classification(r, [E4T, E5I], observed=True)
reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:metadata", "document_id": DOC,
     "kind": "bibliographic_metadata", "evidence_ids": [E2I, E2T],
     "title": "Copper(II)-catalyzed protoboration of allenes in aqueous media and open air",
     "authors": ["Russell F. Snead", "Jan Nekvinda", "Webster L. Santos"],
     "doi": "10.1039/d0nj02010a"},
    {"context_id": f"{DOC}:context:esi_statement", "document_id": DOC,
     "kind": "supplement_availability_statement", "evidence_ids": [E2I, E2T],
     "text": "The article states that ESI is available, but no SI file was supplied; no association was made."},
    {"context_id": f"{DOC}:context:characterization_pages_5_6", "document_id": DOC,
     "kind": "characterization", "evidence_ids": [E5I, E5T, E6I, E6T],
     "text": "NMR and MS characterization for product mixtures 2a–3l is retained as context, not separate synthesis records."},
    {"context_id": f"{DOC}:context:cover", "document_id": DOC,
     "kind": "cover_art", "evidence_ids": [E1I],
     "text": "Journal cover artwork; no concrete experimental record."},
]

st.write(OUT, bundle)
st.write(RAW / "stereo_diagnostics.json", stereo_diagnostics)
mentions = st.source_mentions(data, ["1a", "2a", "3a", "2h", "2n", "2o", "2p"], context=220)
st.write(RAW / "source_mentions.json", mentions)

visual_inventory = {
    "document_id": DOC,
    "pages": [
        {"page": 1, "regions": [{"region": "full page", "finding": "journal cover artwork; no experiment", "evidence_ids": [E1I]}]},
        {"page": 2, "regions": [{"region": "title/abstract/introduction", "finding": "bibliographic metadata, abstract, prior-work prose, ESI availability statement; no concrete arrow/table", "evidence_ids": [E2I]}]},
        {"page": 3, "regions": [
            {"region": "Fig. 1a–d", "finding": "four cited prior-work protoboration arrows", "evidence_ids": [E3I]},
            {"region": "Fig. 1e", "finding": "current-work generic scheme", "evidence_ids": [E3I]},
            {"region": "Table 1 rows 1–21 and footnotes a–f", "finding": "21 optimization experiments", "evidence_ids": [E3I]},
            {"region": "results/discussion prose", "finding": "interpretation and duplicate row references", "evidence_ids": [E3I, E3T]},
        ]},
        {"page": 4, "regions": [
            {"region": "Fig. 2 generic arrow/caption", "finding": "current-work scheme conditions and footnotes", "evidence_ids": [E4I]},
            {"region": "Fig. 2 products 2a–2q", "finding": "17 substrate-scope experiments with yields/ratios", "evidence_ids": [E4I]},
            {"region": "mechanism discussion", "finding": "author-proposed mechanism leading into Fig. 3", "evidence_ids": [E4I, E4T]},
        ]},
        {"page": 5, "regions": [
            {"region": "Fig. 3", "finding": "proposed mechanism; not observed experiments", "evidence_ids": [E5I]},
            {"region": "general experimental details/procedure", "finding": "general borylation operations, workup, purification", "evidence_ids": [E5I, E5T]},
            {"region": "2a/3a, 2b/3b, start of 2c/3c characterization", "finding": "same Fig. 2 experiments; isolated yields/ratios plus spectra, no duplicate reactions", "evidence_ids": [E5I, E5T]},
        ]},
        {"page": 6, "regions": [
            {"region": "2c/3c through 2l/3l product paragraphs", "finding": "same Fig. 2 experiments; yields/ratios and characterization", "evidence_ids": [E6I, E6T]},
            {"region": "NMR/MS peak lists", "finding": "characterization context, not synthesis records", "evidence_ids": [E6I, E6T]},
        ]},
    ],
}
st.write(RAW / "visual_inventory_pages_0001_0006.json", visual_inventory)

byid = {r["reaction_id"]: r for r in reactions}
prior_ids = [f"{DOC}:fig1a_ma_2012", f"{DOC}:fig1b_hoveyda_2013",
             f"{DOC}:fig1c_semba_2013", f"{DOC}:fig1d_hoveyda_2014"]
opt_ids = [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 22)]
scope_ids = [f"{DOC}:scope_2{x}" for x in scope]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": visual_inventory["pages"][0]["regions"],
         "record_ids": [], "unresolved": ["Cover artwork is illustrative and was not treated as an experimental scheme."]},
        {"page": 2, "state": "attempted", "regions": visual_inventory["pages"][1]["regions"],
         "record_ids": [], "unresolved": ["ESI is mentioned but was not supplied or associated."]},
        {"page": 3, "state": "attempted", "regions": visual_inventory["pages"][2]["regions"],
         "record_ids": prior_ids + [scheme_id] + opt_ids,
         "unresolved": ["Complete structures for several cited-work ligands/products are not safely recoverable from the compact figure.",
                        "Optimization-table absolute scale and solvent volumes are unreported."]},
        {"page": 4, "state": "attempted", "regions": visual_inventory["pages"][3]["regions"],
         "record_ids": [scheme_id] + scope_ids + [f"{DOC}:fig3_proposed_mechanism"],
         "unresolved": ["Fig. 2 does not state the ratio measurement method.",
                        "Fig. 2 values for 2o/2p and the later experimental paragraphs conflict; no value was selected."]},
        {"page": 5, "state": "attempted", "regions": visual_inventory["pages"][4]["regions"],
         "record_ids": [f"{DOC}:fig3_proposed_mechanism"] + scope_ids,
         "unresolved": ["General procedure gives equal water/CuSO4-solution volumes but no numeric volumes or substrate scale.",
                        "Proposed intermediates are not observed experiments."]},
        {"page": 6, "state": "attempted", "regions": visual_inventory["pages"][5]["regions"],
         "record_ids": [f"{DOC}:scope_2{x}" for x in "cdefghijkl"],
         "unresolved": ["3e and 3j product names conflict with their Fig. 2 substitution patterns.",
                        "Spectral peak lists were retained as context and not transformed into reaction records."]},
    ],
}
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5, 6]
assert all(p["regions"] for p in coverage["pages"])
assert all(rid in byid for p in coverage["pages"] for rid in p["record_ids"])
st.write(ROOT / "runs/current/coverage/pages_0001_0006.json", coverage)

# Mark only the six required page tasks attempted; later pages remain pending.
tasks_path = ROOT / "runs/current/tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task["document_id"] == DOC and task["evidence_id"] in {E1T, E2T, E3T, E4T, E5T, E6T}:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; decisions mapped in coverage/pages_0001_0006.json."
st.write(tasks_path, tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "scope_records": len(scope_ids), "optimization_records": len(opt_ids)}, ensure_ascii=False))
