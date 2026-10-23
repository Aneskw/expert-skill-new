#!/usr/bin/env python3
"""Build source-derived v3.1 candidates for PDF pages 1-10.

This file is an extraction trace, not a chemistry-validation program.  It uses
the repository's unchanged standard-record helpers and keeps all records in
needs_review/partial state.
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
DOC = "doc_2f5be1806b6542bb"
TITLE = (
    "Regio- and Stereoselective Synthesis of Highly Functionalized "
    "Tetrasubstituted Olefins by Iodine-Mediated Iodofunctionalization "
    "of Ferrocene-Containing Allenylphosphonates"
)


def ev(*pages: int) -> list[str]:
    result = []
    for page in pages:
        result.extend([f"{DOC}:page{page}:image", f"{DOC}:page{page}:text"])
    return result


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    if normalized is None and raw is not None:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def missing(evidence=(), status="unreported"):
    return st.fact(evidence=evidence, status=status)


def not_applicable(evidence=()):
    return st.fact(evidence=evidence, method="not_applicable",
                   status="not_applicable")


def q(number, unit, raw_text=None, qualifier=None):
    result = {"number": number, "unit": unit}
    if raw_text is not None:
        result["raw_text"] = raw_text
    if qualifier is not None:
        result["qualifier"] = qualifier
    return result


def add_occurrence(molecule, label, evidence):
    molecule["occurrences"] = [{
        "document_id": DOC,
        "compound_label": F(label, label, evidence),
        "evidence_ids": list(evidence),
    }]


stereo_checks = {}
molecules = {}


def add_concrete(mid, label, name, smiles, evidence, kind="ordinary",
                 source_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    add_occurrence(m, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "smiles"
    derivation = (
        "Candidate SMILES transcribed from the source drawing/name; RDKit "
        "parsing and canonicalization do not verify source identity."
    )
    m["smiles"] = F({"source": "drawn or named structure", "label": label},
                    smiles, evidence, "structure_derivation",
                    derivation=derivation)
    m["generic_template"] = None
    diag = st.stereo_inspect(smiles)
    stereo_checks[mid] = diag
    iso = diag["canonical_isomeric_smiles"]
    m["stereochemistry"]["isomeric_smiles"] = F(
        None, iso, evidence, "structure_derivation",
        derivation="RDKit canonical isomeric SMILES for the candidate graph."
    )
    m["stereochemistry"]["identity_key"] = F(
        None, diag["stereo_graph_key"], evidence, "structure_derivation",
        derivation="RDKit stereo graph key; not a sample identity claim."
    )
    m["stereochemistry"]["mixture"] = missing(evidence)
    unspecified_atoms = [p for p in diag["potential_stereo"]
                         if str(p["type"]).startswith("Atom_") and p["specified"] == "Unspecified"]
    unspecified_bonds = [p for p in diag["potential_stereo"]
                         if str(p["type"]).startswith("Bond_") and p["specified"] == "Unspecified"]
    if kind == "allene":
        m["stereochemistry"]["atom_centers"] = not_applicable(evidence)
        m["stereochemistry"]["double_bonds"] = F(
            "cumulated C=C=C system; no E/Z descriptor reported",
            [{"type": "allene", "geometry": "not assigned as E/Z"}], evidence)
        m["stereochemistry"]["axial"] = F(
            "allene axis present; no R_a/S_a or P/M assignment",
            [{"type": "allene", "configuration": "unspecified"}], evidence)
        m["stereochemistry"]["assignment"] = F(
            "no axial configuration stated", "unspecified", evidence)
        m["stereochemistry"]["chirality_status"] = F(
            None, "unspecified", evidence, "inference",
            derivation="Both allene termini bear non-equivalent substituents, but the source assigns no axial configuration.")
        m["stereochemistry"]["geometry_status"] = F(
            None, "not_applicable", evidence, "inference",
            derivation="The cumulated allene is not assigned an alkene E/Z descriptor.")
        m["stereochemistry"]["representation_status"] = F(
            None, "partial", evidence, "inference",
            derivation="Connectivity is encoded; allene axial stereochemistry is unspecified and not encoded.")
    elif kind == "E_product":
        if not any(d["stereo"] == "STEREOE" for d in diag["double_bonds"]):
            raise ValueError(f"E product candidate did not retain E geometry: {mid}")
        m["stereochemistry"]["atom_centers"] = F(
            "new tetrahedral carbon is drawn without an absolute descriptor",
            [{"site": "ferrocenyl-bearing carbon", "configuration": "unspecified"}], evidence)
        m["stereochemistry"]["double_bonds"] = F(
            f"{label} is explicitly labelled E",
            [{"site": "new tetrasubstituted alkene", "configuration": "E"}], evidence)
        m["stereochemistry"]["axial"] = not_applicable(evidence)
        m["stereochemistry"]["assignment"] = F(
            "E alkene; no R/S assignment for the tetrahedral carbon",
            {"alkene": "E", "tetrahedral_center": "unspecified"}, evidence)
        m["stereochemistry"]["chirality_status"] = F(
            None, "unspecified", evidence, "inference",
            derivation="The product graph has a potential tetrahedral stereocenter, but the source reports only E alkene geometry.")
        m["stereochemistry"]["geometry_status"] = F(
            "E", "specified", evidence)
        m["stereochemistry"]["representation_status"] = F(
            None, "partial", evidence, "inference",
            derivation="E alkene geometry is encoded; the potential tetrahedral center remains unspecified.")
    else:
        m["stereochemistry"]["atom_centers"] = F(
            None, diag["atom_centers"] or [], evidence, "structure_derivation",
            derivation="RDKit potential/assigned tetrahedral-center diagnostic for the candidate graph.")
        m["stereochemistry"]["double_bonds"] = F(
            None, diag["double_bonds"] or [], evidence, "structure_derivation",
            derivation="RDKit double-bond stereo diagnostic for the candidate graph.")
        m["stereochemistry"]["axial"] = not_applicable(evidence)
        m["stereochemistry"]["assignment"] = missing(evidence)
        chirality = "unspecified" if unspecified_atoms else "achiral"
        geometry = "unspecified" if unspecified_bonds else "not_applicable"
        representation = "partial" if (unspecified_atoms or unspecified_bonds) else "encoded"
        m["stereochemistry"]["chirality_status"] = F(
            None, chirality, evidence, "inference",
            derivation="Status follows candidate-graph stereo diagnostics; it is not inferred from absence of source wedges alone.")
        m["stereochemistry"]["geometry_status"] = F(
            None, geometry, evidence, "inference",
            derivation="Status follows candidate-graph double-bond diagnostics.")
        m["stereochemistry"]["representation_status"] = F(
            None, representation, evidence, "inference",
            derivation="Candidate graph encoding status after stereo diagnostics.")
    note = source_note or "Machine-transcribed candidate; source identity and stereochemical completeness require expert review."
    if "[Fe+2]" in smiles:
        note += " Ferrocene is represented as disconnected charged Cp/Fe components; this is a graph surrogate."
    m["identity_note"] = F(None, note, evidence, "inference", derivation=note)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def add_unresolved(mid, label, name, evidence, E_product=False, note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    add_occurrence(m, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(evidence)
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = missing(evidence)
    m["stereochemistry"]["axial"] = not_applicable(evidence) if E_product else missing(evidence)
    m["stereochemistry"]["mixture"] = missing(evidence)
    if E_product:
        m["stereochemistry"]["double_bonds"] = F(
            f"{label} is explicitly labelled E",
            [{"site": "new tetrasubstituted alkene", "configuration": "E"}], evidence)
        m["stereochemistry"]["assignment"] = F("E", "E", evidence)
        m["stereochemistry"]["chirality_status"] = F(
            None, "unresolved", evidence, "inference",
            derivation="No complete molecular graph was encoded, so other stereogenic elements remain unresolved.")
        m["stereochemistry"]["geometry_status"] = F("E", "specified", evidence)
        m["stereochemistry"]["representation_status"] = F(
            None, "not_encoded", evidence, "inference",
            derivation="The E assignment is retained textually, but the complete structure was not safely encoded.")
    else:
        m["stereochemistry"]["double_bonds"] = missing(evidence)
        m["stereochemistry"]["assignment"] = missing(evidence)
        for key in ("chirality_status", "geometry_status", "representation_status"):
            m["stereochemistry"][key] = F(
                None, "unresolved", evidence, "inference",
                derivation="Complete identity/structure is unresolved from the bounded source region.")
    m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
    m["stereochemistry"]["identity_key"] = missing(evidence)
    message = note or "Complete structure intentionally left unresolved rather than guessed."
    m["identity_note"] = F(None, message, evidence, "inference", derivation=message)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def add_generic(mid, label, name, fixed, variables, evidence):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    add_occurrence(m, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(evidence)
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": F(f"source attachment point for {v}", f"source attachment point for {v}", evidence),
            "definition": F(definition, definition, evidence),
            "allowed_values": missing(evidence),
            "constraints": missing(evidence),
        } for v, definition in variables],
        "instance_mapping": [],
        "reason": F("generic R-group drawing", "Source uses variable substituents; no single concrete molecule is asserted.", evidence),
    }
    m["stereochemistry"]["atom_centers"] = missing(evidence)
    m["stereochemistry"]["double_bonds"] = missing(evidence)
    m["stereochemistry"]["axial"] = missing(evidence)
    m["stereochemistry"]["mixture"] = missing(evidence)
    m["stereochemistry"]["assignment"] = missing(evidence)
    for key in ("chirality_status", "geometry_status", "representation_status"):
        m["stereochemistry"][key] = F(
            None, "unresolved", evidence, "inference",
            derivation="Generic variable structure cannot carry one concrete stereochemical assignment.")
    m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
    m["stereochemistry"]["identity_key"] = missing(evidence)
    m["identity_note"] = F(None, "Generic template from the source; variables are not expanded beyond stated scope.",
                           evidence, "inference", derivation="Generic source drawing.")
    m["external_properties"] = []
    molecules[mid] = m
    return mid


FC = "[Fe+2].[cH-]1cccc1.[c-]1({sub})cccc1"


def allene_smiles(r1):
    return FC.format(sub=f"C=C=C(P(=O)(OC)OC){r1}")


def product_smiles(oxy, r1):
    # Opposite slash directions encode the E relationship of I and P (CIP-high groups).
    return FC.format(sub=f"C({oxy})C(/I)=C(\\P(=O)(OC)OC){r1}")


def ensure_small_molecules():
    small = [
        ("m_i2", "I2", "iodine", "II", ev(1,2,3,4,5)),
        ("m_tbuoli", "t-BuOLi", "lithium tert-butoxide", "[Li+].CC(C)(C)[O-]", ev(1,2,3,4,5)),
        ("m_naoh", "NaOH", "sodium hydroxide", "[Na+].[OH-]", ev(2)),
        ("m_koh", "KOH", "potassium hydroxide", "[K+].[OH-]", ev(2)),
        ("m_cs2co3", "Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O", ev(2)),
        ("m_lioh", "LiOH", "lithium hydroxide", "[Li+].[OH-]", ev(2)),
        ("m_dbu", "DBU", "1,8-diazabicyclo[5.4.0]undec-7-ene", "C1CCC2=NCCCN2CC1", ev(2)),
        ("m_thf", "THF", "tetrahydrofuran", "C1CCOC1", ev(2)),
        ("m_mecn", "CH3CN", "acetonitrile", "CC#N", ev(2)),
        ("m_meoh", "CH3OH", "methanol", "CO", ev(2)),
        ("m_toluene", "toluene", "toluene", "Cc1ccccc1", ev(2)),
        ("m_dioxane", "dioxane", "1,4-dioxane", "O1CCOCC1", ev(2)),
        ("m_dcm", "CH2Cl2", "dichloromethane", "ClCCl", ev(2,3,4,5)),
        ("m_etoac", "EtOAc", "ethyl acetate", "CCOC(C)=O", ev(5)),
    ]
    for row in small:
        add_concrete(*row)
    add_unresolved("m_petroleum_ether", "petroleum ether", "petroleum ether", ev(5),
                   note="Named distillation-range mixture; composition is not specified, so no unique SMILES is assigned.")


ensure_small_molecules()

# Generic source drawings.
add_generic("m_generic_allene", "1", "ferrocene-containing allenylphosphonate 1",
            "Fc-CH=C=C(R)-P(=O)(OMe)2", [("R", "allene substituent")], ev(1,2,3,4))
add_generic("m_generic_oxygen_nucleophile", "R'OH", "carboxylic acid, alcohol, or phenol",
            "R'-OH", [("R'", "carboxylic-acid, alcohol, or phenol residue")], ev(1,2))
add_generic("m_generic_E_product", "E-product", "E-tetrasubstituted alkenylphosphonate",
            "Fc-CH(OR')-C(I)=C(R)-P(=O)(OMe)2, E",
            [("R", "allene substituent"), ("R'", "oxygen-nucleophile residue")], ev(1,2))
add_generic("m_generic_carboxylic_acid", "2", "carboxylic acid 2", "R2-C(=O)OH",
            [("R2", "carboxylic-acid substituent")], ev(3,5))
add_generic("m_generic_E3", "E-3a-t", "E ester product 3",
            "Fc-CH(O-C(=O)R2)-C(I)=C(R1)-P(=O)(OMe)2, E",
            [("R1", "allene substituent"), ("R2", "acid substituent")], ev(3))
add_generic("m_generic_alcohol", "4", "alcohol 4", "R2-OH",
            [("R2", "alcohol substituent")], ev(4))
add_generic("m_generic_E5", "E-5a-j", "E ether product 5",
            "Fc-CH(OR2)-C(I)=C(R1)-P(=O)(OMe)2, E",
            [("R1", "allene substituent"), ("R2", "alcohol substituent")], ev(4))
add_generic("m_generic_phenol", "6", "phenol 6", "Ar-OH",
            [("R", "phenyl-ring substituent")], ev(5))
add_generic("m_generic_E7", "E-7a-i", "E aryl ether product 7",
            "Fc-CH(OAr)-C(I)=C(Ph)-P(=O)(OMe)2, E",
            [("R", "phenoxy-ring substituent")], ev(5))

# Concrete allene substrates; the source does not assign allene axial configuration.
r1_fragments = {
    "ph": "c2ccccc2",
    "35me": "c2cc(C)cc(C)c2",
    "4me": "c2ccc(C)cc2",
    "4tbu": "c2ccc(C(C)(C)C)cc2",
    "4ome": "c2ccc(OC)cc2",
    "4oh": "c2ccc(O)cc2",
    "4biphenyl": "c2ccc(-c3ccccc3)cc2",
    "4br": "c2ccc(Br)cc2",
    "4cl": "c2ccc(Cl)cc2",
    "ethyl": "CC",
}
allene_defs = {
    "1a": ("ph", ev(2,3,4,5)),
    "1l": ("35me", ev(3)),
    "1m": ("4me", ev(3)),
    "1n": ("4tbu", ev(3)),
    "1o": ("4ome", ev(3)),
    "1p": ("4oh", ev(3)),
    "1q": ("4biphenyl", ev(3)),
    "1r": ("4br", ev(3)),
    "1s": ("4cl", ev(3)),
    "1t": ("ethyl", ev(3)),
}
for label, (rkey, evidence) in allene_defs.items():
    add_concrete(f"m_{label}", label, f"allenylphosphonate {label}",
                 allene_smiles(r1_fragments[rkey]), evidence, kind="allene")

# Carboxylic acids 2a-k.
acid_defs = {
    "2a": ("benzoic acid", "O=C(O)c1ccccc1"),
    "2b": ("4-bromobenzoic acid", "O=C(O)c1ccc(Br)cc1"),
    "2c": ("3-methoxybenzoic acid", "O=C(O)c1cccc(OC)c1"),
    "2d": ("4-chlorobenzoic acid", "O=C(O)c1ccc(Cl)cc1"),
    "2e": ("3-iodobenzoic acid", "O=C(O)c1cccc(I)c1"),
    "2f": ("3-chlorobenzoic acid", "O=C(O)c1cccc(Cl)c1"),
    "2g": ("3-fluorobenzoic acid", "O=C(O)c1cccc(F)c1"),
    "2h": ("4-cyanobenzoic acid", "N#Cc1ccc(C(=O)O)cc1"),
    "2i": ("indole-3-acetic acid", "O=C(O)Cc1c[nH]c2ccccc12"),
    "2j": ("hexanoic acid", "CCCCCC(=O)O"),
}
for label, (name, smiles) in acid_defs.items():
    add_concrete(f"m_{label}", label, name, smiles, ev(3, 5) if label == "2a" else ev(3))
add_unresolved("m_2k", "2k", "long-chain aliphatic carboxylic acid 2k", ev(3),
               note="The bounded drawing/printed product name does not establish the long-chain carbon count confidently; SMILES is withheld.")

# Alcohols 4a-g and phenols 6a-i.
alcohol_defs = {
    "4a": ("ethanol", "CCO"),
    "4b": ("2-bromoethanol", "OCCBr"),
    "4c": ("(4-iodophenyl)methanol", "OCc1ccc(I)cc1"),
    "4d": ("but-3-yn-1-ol", "OCCC#C"),
    "4e": ("octan-2-ol", "CC(O)CCCCCC"),
    "4f": ("cyclohexanol", "OC1CCCCC1"),
    "4g": ("tert-butanol", "CC(C)(C)O"),
}
for label, (name, smiles) in alcohol_defs.items():
    add_concrete(f"m_{label}", label, name, smiles, ev(4))
phenol_defs = {
    "6a": ("phenol", "Oc1ccccc1"),
    "6b": ("4-methylphenol", "Cc1ccc(O)cc1"),
    "6c": ("4-tert-butylphenol", "CC(C)(C)c1ccc(O)cc1"),
    "6d": ("3,5-dinitrophenol", "Oc1cc([N+](=O)[O-])cc([N+](=O)[O-])c1"),
    "6e": ("4-bromophenol", "Oc1ccc(Br)cc1"),
    "6f": ("3-fluorophenol", "Oc1cccc(F)c1"),
    "6g": ("4-chloro-3-methylphenol", "Oc1cc(C)c(Cl)cc1"),
    "6h": ("2-chlorophenol", "Oc1ccccc1Cl"),
    "6i": ("4-phenylphenol", "Oc1ccc(-c2ccccc2)cc1"),
}
for label, (name, smiles) in phenol_defs.items():
    add_concrete(f"m_{label}", label, name, smiles, ev(5))

# Positive E-products.  Product 3k remains unresolved because its acid chain is ambiguous.
acid_oxy = {
    "3a": "OC(=O)c2ccccc2",
    "3b": "OC(=O)c2ccc(Br)cc2",
    "3c": "OC(=O)c2cccc(OC)c2",
    "3d": "OC(=O)c2ccc(Cl)cc2",
    "3e": "OC(=O)c2cccc(I)c2",
    "3f": "OC(=O)c2cccc(Cl)c2",
    "3g": "OC(=O)c2cccc(F)c2",
    "3j": "OC(=O)CCCCC",
    "3l": "OC(=O)c2ccccc2",
    "3m": "OC(=O)c2ccccc2",
    "3n": "OC(=O)c2ccccc2",
    "3o": "OC(=O)c2ccccc2",
    "3p": "OC(=O)c2ccccc2",
    "3q": "OC(=O)c2ccccc2",
    "3r": "OC(=O)c2ccccc2",
    "3s": "OC(=O)c2ccccc2",
    "3t": "OC(=O)c2ccccc2",
}
product_r1 = {
    **{k: "ph" for k in ["3a","3b","3c","3d","3e","3f","3g","3j"]},
    "3l":"35me", "3m":"4me", "3n":"4tbu", "3o":"4ome", "3p":"4oh",
    "3q":"4biphenyl", "3r":"4br", "3s":"4cl", "3t":"ethyl",
}
product_extra_pages = {
    "3a": ev(3,5), "3b": ev(3,5), "3c": ev(3,5,6),
    "3d": ev(3,6), "3e": ev(3,6), "3f": ev(3,6), "3g": ev(3,6),
    "3j": ev(3,6), "3l": ev(3,6), "3m": ev(3,6), "3n": ev(3,6,7),
    "3o": ev(3,7), "3p": ev(3,7), "3q": ev(3,7),
    "3r": ev(3,7), "3s": ev(3,7), "3t": ev(3,7),
}
for label, oxy in acid_oxy.items():
    add_concrete(f"m_{label}", f"E-{label}", f"(E)-product {label}",
                 product_smiles(oxy, r1_fragments[product_r1[label]]),
                 product_extra_pages.get(label, ev(3)), kind="E_product")
add_unresolved("m_3k", "E-3k", "(E)-product 3k", ev(3,6), E_product=True,
               note="E geometry and product label are explicit, but the long-chain substituent length is ambiguous in the bounded source image/OCR.")

alcohol_oxy = {
    "5a": "OCC", "5b": "OCCBr", "5c": "OCc2ccc(I)cc2", "5d": "OCCC#C",
    "5e": "OC(C)CCCCCC", "5f": "OC2CCCCC2", "5g": "OC(C)(C)C",
    "5h": "OCC", "5i": "OCC", "5j": "OCC",
}
alcohol_r1 = {**{k:"ph" for k in ["5a","5b","5c","5d","5e","5f","5g"]},
              "5h":"35me", "5i":"4biphenyl", "5j":"4br"}
alcohol_product_evidence = {
    "5a": ev(4,7), "5b": ev(4,7), "5c": ev(4,7),
    "5d": ev(4,7,8), "5e": ev(4,8), "5f": ev(4,8),
    "5g": ev(4,8), "5h": ev(4,8), "5i": ev(4,8), "5j": ev(4,8),
}
for label, oxy in alcohol_oxy.items():
    add_concrete(f"m_{label}", f"E-{label}", f"(E)-product {label}",
                 product_smiles(oxy, r1_fragments[alcohol_r1[label]]),
                 alcohol_product_evidence[label], kind="E_product")

phenol_oxy = {
    "7a":"Oc2ccccc2", "7b":"Oc2ccc(C)cc2", "7c":"Oc2ccc(C(C)(C)C)cc2",
    "7e":"Oc2ccc(Br)cc2", "7f":"Oc2cccc(F)c2", "7g":"Oc2cc(C)c(Cl)cc2",
    "7i":"Oc2ccc(-c3ccccc3)cc2",
}
phenol_product_evidence = {
    "7a": ev(5,8), "7b": ev(5,8), "7c": ev(5,8),
    "7e": ev(5,8,9), "7f": ev(5,9), "7g": ev(5,9), "7i": ev(5,9),
}
for label, oxy in phenol_oxy.items():
    add_concrete(f"m_{label}", f"E-{label}", f"(E)-product {label}",
                 product_smiles(oxy, r1_fragments["ph"]),
                 phenol_product_evidence[label], kind="E_product")


def charge_fact(participant, key, raw, normalized, evidence, reference=None,
                method="original_extraction", derivation=None):
    participant["charges"][key] = F(raw, normalized, evidence, method,
                                     derivation=derivation)
    if key in ("equivalents", "mol_percent"):
        participant["equivalence_reference"] = F(
            reference, reference, evidence,
            method if method != "structure_derivation" else "inference",
            derivation=derivation)


def participant(pid, mid, role, evidence, stage="reaction", amounts=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    for key, (raw, normalized, method, derivation) in (amounts or {}).items():
        charge_fact(p, key, raw, normalized, evidence, reference, method, derivation)
    return p


def set_reaction_basics(r, categories, location, reaction_name, transformation,
                        evidence, scheme_ids=()):
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F(categories, categories, evidence)
    r["location"] = F(location, location, evidence)
    r["general_scheme_ids"] = list(scheme_ids)
    r["reaction_name"] = F(reaction_name, reaction_name, evidence)
    r["transformation_type"] = F(
        None, transformation, evidence, "inference",
        derivation="Transformation type inferred from the source reactant/product drawings and labels.")
    r["completeness"] = "partial"


def set_solvent(r, name, volume, evidence):
    r["solvent_conditions"]["composition"] = F(name, [name], evidence)
    r["solvent_conditions"]["mixture_ratio"] = missing(evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = F(
            f"{volume} mL {name}", q(volume, "mL"), evidence)
    else:
        r["solvent_conditions"]["total_volume"] = missing(evidence)
    r["solvent_conditions"]["concentration"] = missing(evidence)
    r["solvent_conditions"]["concentration_basis"] = missing(evidence)


def add_step(r, sid, pids, operation, temperature, time, atmosphere, evidence,
             stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence)
    s["temperature"] = (F(temperature, temperature, evidence)
                        if temperature is not None else missing(evidence))
    s["time"] = (F(time, time, evidence)
                 if time is not None else missing(evidence))
    s["atmosphere"] = (F(atmosphere, atmosphere, evidence)
                       if atmosphere is not None else missing(evidence))
    s["pressure"] = missing(evidence)
    r["steps"].append(s)


def add_yield(r, oid, product_mid, percent, basis, evidence, product_mass=None,
              product_mass_evidence=None, qualitative=None):
    if qualitative == "trace":
        o = st.new_outcome(oid, "trace", [product_mid], evidence)
        o["value"] = F("trace", "trace", evidence)
        o["product_role"] = F("desired product", "desired product", evidence)
        o["yield_type"] = missing(evidence)
        o["measurement_basis"] = missing(evidence)
        o["scope"] = F("trace desired product", {"product": product_mid}, evidence)
        r["outcomes"].append(o)
        return
    o = st.new_outcome(oid, "yield", [product_mid], evidence)
    o["value"] = F(f"{percent}%", q(percent, "%"), evidence)
    o["product_role"] = F("desired product", "desired product", evidence)
    o["yield_type"] = F("isolated yield", "isolated", evidence) if basis else missing(evidence)
    o["measurement_basis"] = F(basis, basis, evidence) if basis else missing(evidence)
    o["scope"] = F("yield of the depicted product", {"product": product_mid}, evidence)
    r["outcomes"].append(o)
    if product_mass is not None:
        mass_evidence = product_mass_evidence or evidence
        p = st.new_outcome(oid + "_mass", "product_amount", [product_mid], mass_evidence)
        p["value"] = F(f"{product_mass} mg", q(product_mass, "mg"), mass_evidence)
        p["product_role"] = F("desired product", "desired product", mass_evidence)
        p["yield_type"] = not_applicable(mass_evidence)
        p["measurement_basis"] = F("isolated product mass", "isolated", mass_evidence)
        p["scope"] = F("mass reported in the product preparation paragraph", {"product": product_mid}, mass_evidence)
        r["outcomes"].append(p)


def add_no_reaction(r, oid, raw, evidence, attempted_label):
    o = st.new_outcome(oid, "no_reaction", [], evidence)
    o["value"] = F(raw, q(0, "%") if raw == "0%" else raw, evidence)
    o["product_role"] = not_applicable(evidence)
    o["yield_type"] = not_applicable(evidence)
    o["measurement_basis"] = missing(evidence)
    o["scope"] = F(f"attempted formation of {attempted_label}",
                   {"attempted_product_label": attempted_label}, evidence)
    r["outcomes"].append(o)


def add_selectivity(r, oid, product_mid, raw, evidence):
    o = st.new_outcome(oid, "selectivity", [product_mid], evidence)
    o["value"] = F(raw, raw, evidence)
    o["product_role"] = F("desired E product", "desired product", evidence)
    o["yield_type"] = not_applicable(evidence)
    o["measurement_basis"] = missing(evidence)
    o["scope"] = F(
        "alkene geometry; E product is shown/reported, numerical E/Z ratio unreported",
        {"metric": "alkene geometry", "major": "E", "ratio": None,
         "target": product_mid}, evidence)
    r["outcomes"].append(o)


def set_classification(r, evidence, state="positive"):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(evidence)
    c["model_votes"] = []
    if state == "positive":
        changes = [
            {"change_id":"c1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"O"}],
             "before_order":0,"after_order":1,"change_type":"formed"},
            {"change_id":"c2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"I"}],
             "before_order":0,"after_order":1,"change_type":"formed"},
            {"change_id":"c3","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],
             "before_order":2,"after_order":1,"change_type":"order_decreased"},
            {"change_id":"c4","atoms":[{"atom_map":None,"element":"O"},{"atom_map":None,"element":"H"}],
             "before_order":1,"after_order":0,"change_type":"broken"},
        ]
        c["bond_changes"] = [F(None, x, evidence, "structure_derivation",
                               derivation="Partial net bond-change reading of the source reaction drawing; no complete atom map asserted.")
                             for x in changes]
        labels = [
            {"category_number":3,"code":"C-O","supporting_change_ids":["c1"]},
            {"category_number":10,"code":"C-I","supporting_change_ids":["c2"]},
        ]
        c["labels"] = F(None, labels, evidence, "structure_derivation",
                        derivation="Only newly formed C-O and C-I bonds are classified; the decreased allene C=C bond is not relabelled C-C.")
        c["coverage"] = F(None, "partial", evidence, "inference",
                          derivation="Source-derived partial net changes; atom mapping, base coproducts, and iodine balance are unresolved.")
        c["limitations"] = [
            "No complete atom-mapped reaction SMILES is asserted.",
            "Ferrocene connectivity is represented by a graph surrogate in candidate SMILES.",
            "Base/iodide coproducts and proton transfer partners are not fully mapped.",
        ]
    elif state == "no_reaction":
        c["bond_changes"] = []
        c["labels"] = F(None, [], evidence, "inference",
                        derivation="No desired product was observed; no intended bond formation is assigned.")
        c["coverage"] = F(None, "not_applicable", evidence, "inference",
                          derivation="Negative experiment; expected transformation is not classified as observed chemistry.")
        c["limitations"] = ["No product/detection-limit structure was reported for the failed experiment."]
    else:
        c["bond_changes"] = []
        c["labels"] = missing(evidence)
        c["coverage"] = F(None, "unresolved", evidence, "inference",
                          derivation="The source presents a proposed multistep mechanism rather than one experimentally mapped net reaction.")
        c["limitations"] = ["Proposed intermediates and elementary arrows are not treated as experimentally verified atom mappings."]


def finish_review(r, unresolved):
    r["review"] = {"status":"needs_review", "method":"model", "reviewer":None,
                   "checked_fields":[], "unresolved":list(unresolved)}
    r["completeness"] = "partial"


reactions = []


def add_generic_reaction(rid, page, participants_spec, name, location, procedure,
                         classification="positive"):
    evidence = ev(page)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = "generic_scheme"
    r["experiment_categories"] = F(["other"], ["other"], evidence)
    r["location"] = F(location, location, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(None, "iodofunctionalization of an allenylphosphonate", evidence,
                                  "inference", derivation="Inferred from the generic source arrow.")
    for pid, mid, role in participants_spec:
        r["participants"].append(participant(pid, mid, role, evidence,
                                             product_role="generic product" if role == "product" else None))
    r["procedure"]["original_text"] = F(procedure, procedure, evidence)
    r["procedure"]["general_procedure_reference"] = missing(evidence)
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = missing(evidence)
    r["procedure"]["pretreatment"] = missing(evidence)
    set_classification(r, evidence, classification)
    finish_review(r, ["Generic R-group scheme; not a single concrete experiment.",
                      "Variable identities, time, solvent volume, and workup are not all specified on this arrow."])
    reactions.append(r)
    return r["reaction_id"]


scheme1_id = add_generic_reaction(
    "scheme1_general", 2,
    [("gs1_allene","m_generic_allene","reactant"),
     ("gs1_nuc","m_generic_oxygen_nucleophile","reactant"),
     ("gs1_i2","m_i2","reagent"), ("gs1_base","m_tbuoli","reagent"),
     ("gs1_product","m_generic_E_product","product")],
    "Synthesis of Ferrocene-Containing Tetrasubstituted Alkenylphosphonates",
    "Scheme 1 (also summarized in the page-1 abstract graphic)",
    "I2 (1 equiv), t-BuOLi (1 equiv); R'OH = carboxylic acids, phenols and alcohols."
)
reactions[-1]["evidence_ids"] = ev(1,2)

scheme2_id = add_generic_reaction(
    "scheme2_general", 3,
    [("gs2_allene","m_generic_allene","reactant"),
     ("gs2_acid","m_generic_carboxylic_acid","reactant"),
     ("gs2_i2","m_i2","reagent"), ("gs2_base","m_tbuoli","reagent"),
     ("gs2_product","m_generic_E3","product")],
    "Iodine-Mediated Iodofunctionalization of Various Allenylphosphonates and Carboxylic Acids",
    "Scheme 2 generic arrow", "I2 (1.0 equiv), t-BuOLi (1.0 equiv), CH2Cl2, r.t."
)
scheme3_id = add_generic_reaction(
    "scheme3_general", 4,
    [("gs3_allene","m_generic_allene","reactant"),
     ("gs3_alcohol","m_generic_alcohol","reactant"),
     ("gs3_i2","m_i2","reagent"), ("gs3_base","m_tbuoli","reagent"),
     ("gs3_product","m_generic_E5","product")],
    "Highly Regio- and Stereoselective Iodoalkoxylation of Allenylphosphonates with Various Alcohols",
    "Scheme 3 generic arrow", "I2 (1.0 equiv), t-BuOLi (1.0 equiv), CH2Cl2, r.t."
)
scheme4_id = add_generic_reaction(
    "scheme4_general", 5,
    [("gs4_allene","m_1a","reactant"),
     ("gs4_phenol","m_generic_phenol","reactant"),
     ("gs4_i2","m_i2","reagent"), ("gs4_base","m_tbuoli","reagent"),
     ("gs4_product","m_generic_E7","product")],
    "Synthesis of Tetrasubstituted Alkenylphosphonates by the Reaction of Allenylphosphonate 1a with Various Phenols",
    "Scheme 4 generic arrow", "I2 (1.0 equiv), t-BuOLi (1.0 equiv), CH2Cl2, r.t."
)


# Table 1: each row is a distinct optimization experiment.
optimization_rows = [
    (1,"m_thf","THF","m_naoh","NaOH","r.t.",45),
    (2,"m_mecn","CH3CN","m_naoh","NaOH","r.t.",9),
    (3,"m_meoh","CH3OH","m_naoh","NaOH","r.t.",19),
    (4,"m_toluene","toluene","m_naoh","NaOH","r.t.",35),
    (5,"m_dioxane","dioxane","m_naoh","NaOH","r.t.",42),
    (6,"m_dcm","CH2Cl2","m_naoh","NaOH","r.t.",47),
    (7,"m_dcm","CH2Cl2","m_koh","KOH","r.t.",17),
    (8,"m_dcm","CH2Cl2","m_cs2co3","Cs2CO3","r.t.",39),
    (9,"m_dcm","CH2Cl2","m_tbuoli","t-BuOLi","r.t.",70),
    (10,"m_dcm","CH2Cl2","m_lioh","LiOH","r.t.",60),
    (11,"m_dcm","CH2Cl2","m_dbu","DBU","r.t.","trace"),
    (12,"m_dcm","CH2Cl2","m_tbuoli","t-BuOLi",q(40,"°C"),54),
    (13,"m_dcm","CH2Cl2","m_tbuoli","t-BuOLi","r.t.",78),
]
for entry, solvent_mid, solvent_name, base_mid, base_name, temp, outcome in optimization_rows:
    evidence = ev(2)
    r = st.new_reaction(f"{DOC}:table1_entry_{entry:02d}", DOC, evidence)
    set_reaction_basics(r, ["optimization"], f"Table 1 entry {entry}",
                        "Optimization of Reaction Conditions",
                        "iodoesterification of allenylphosphonate 1a with benzoic acid",
                        evidence, [scheme1_id])
    if entry < 13:
        amounts = {
            "allene": ("0.20 mmol", q(0.20,"mmol"), "original_extraction", None),
            "acid": ("0.10 mmol", q(0.10,"mmol"), "original_extraction", None),
            "i2": ("0.15 mmol", q(0.15,"mmol"), "original_extraction", None),
            "base": ("0.10 mmol", q(0.10,"mmol"), "original_extraction", None),
        }
    else:
        amounts = {
            "allene": ("0.10 mmol", q(0.10,"mmol"), "original_extraction", None),
            "acid": ("1a/2a = 1.0:2.0", q(2.0,"equiv"), "original_extraction", None),
            "i2": ("1a/iodine = 1.0:1.0", q(1.0,"equiv"), "original_extraction", None),
            "base": ("1a/base = 1.0:1.0", q(1.0,"equiv"), "original_extraction", None),
        }
    r["participants"] = [
        participant("allene", "m_1a", "reactant", evidence,
                    amounts={"amount": amounts["allene"]}),
        participant("acid", "m_2a", "reactant", evidence,
                    amounts={"amount" if entry < 13 else "equivalents": amounts["acid"]},
                    reference="m_1a" if entry == 13 else None),
        participant("iodine", "m_i2", "reagent", evidence,
                    amounts={"amount" if entry < 13 else "equivalents": amounts["i2"]},
                    reference="m_1a" if entry == 13 else None),
        participant("base", base_mid, "reagent", evidence,
                    amounts={"amount" if entry < 13 else "equivalents": amounts["base"]},
                    reference="m_1a" if entry == 13 else None),
        participant("solvent", solvent_mid, "solvent", evidence,
                    amounts={"volume": ("2.0 mL solvent", q(2.0,"mL"), "original_extraction", None)}),
        participant("product", "m_3a", "product", evidence, product_role="desired product"),
    ]
    r["reaction_scale"] = F(
        "1a 0.20 mmol" if entry < 13 else "1a 0.10 mmol",
        {"reference_molecule_id":"m_1a","amount":q(0.20 if entry < 13 else 0.10,"mmol")}, evidence)
    set_solvent(r, solvent_name, 2.0, evidence)
    time = q(1,"h") if entry == 1 else None
    atmosphere = "air" if entry == 1 else None
    add_step(r, "reaction", ["allene","acid","iodine","base","solvent"],
             "model reaction under the tabulated solvent/base/temperature",
             temp, time, atmosphere, evidence)
    add_step(r, "purification", ["product"],
             "purified by column chromatography on silica gel",
             None, None, None, evidence, stage="purification")
    r["procedure"]["original_text"] = F(
        "Table 1 footnotes a-c and row values", {"table":"Table 1","entry":entry}, evidence)
    for key in ("general_procedure_reference","modifications","addition_order","pretreatment"):
        r["procedure"][key] = missing(evidence)
    if outcome == "trace":
        add_yield(r, "trace_3a", "m_3a", None, None, evidence, qualitative="trace")
    else:
        add_yield(r, "yield_3a", "m_3a", outcome,
                  "isolated product purified by column chromatography on silica gel", evidence)
    table_selectivity = (
        "E-3a depicted and described as the only stereoisomer; numerical E/Z ratio unreported"
        if entry == 1 else
        "E-3a depicted; numerical E/Z ratio unreported for this optimization row"
    )
    add_selectivity(r, "selectivity_3a", "m_3a", table_selectivity, evidence)
    set_classification(r, evidence)
    unresolved = ["Numerical E/Z ratio and stereochemical measurement method are unreported."]
    if entry != 1:
        unresolved += ["Reaction time and atmosphere are not stated for this individual optimization row; entry-1 prose is not inherited."]
    if outcome == "trace":
        unresolved += ["Trace detection/measurement method is not reported."]
    finish_review(r, unresolved)
    reactions.append(r)


# Scheme 2 concrete rows.
scope3 = [
    ("3a","1a","2a",1,78,100), ("3b","1a","2b",1,77,110),
    ("3c","1a","2c",1,62,83), ("3d","1a","2d",1,64,87),
    ("3e","1a","2e",1,60,92), ("3f","1a","2f",1,42,57),
    ("3g","1a","2g",1.5,54,75), ("3h","1a","2h",4,0,None),
    ("3i","1a","2i",3,0,None), ("3j","1a","2j",2,46,60),
    ("3k","1a","2k",2,35,54), ("3l","1l","2a",1,78,98),
    ("3m","1m","2a",1,83,105), ("3n","1n","2a",1,94,115),
    ("3o","1o","2a",1,90,113), ("3p","1p","2a",1,82,104),
    ("3q","1q","2a",1,85,103), ("3r","1r","2a",1,76,92),
    ("3s","1s","2a",1,69,86), ("3t","1t","2a",2,53,72),
]
scope3_mass_page = {
    "3a":5, "3b":5, "3c":5, "3d":6, "3e":6, "3f":6, "3g":6,
    "3j":6, "3k":6, "3l":6, "3m":6, "3n":6,
    "3o":7, "3p":7, "3q":7, "3r":7, "3s":7, "3t":7,
}
for label, allene_label, acid_label, hours, percent, mass in scope3:
    extra = product_extra_pages.get(label, ev(3))
    evidence = list(dict.fromkeys(ev(3,5) + extra))
    conditions_evidence = ev(3,5)
    # Page 5 supplies the general carboxylic-acid procedure to every Scheme-2 row;
    # pages 6-7 are added only where a preparation/characterization paragraph is visible.
    r = st.new_reaction(f"{DOC}:scope_{label}", DOC, evidence)
    set_reaction_basics(r, ["substrate_scope"], f"Scheme 2, {label}",
                        "Iodine-Mediated Iodofunctionalization of Various Allenylphosphonates and Carboxylic Acids",
                        "iodoesterification of an allenylphosphonate",
                        evidence, [scheme2_id])
    r["participants"] = [
        participant("allene", f"m_{allene_label}", "reactant", conditions_evidence,
                    amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),
                             "equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)}, reference=f"m_{allene_label}"),
        participant("acid", f"m_{acid_label}", "reactant", conditions_evidence,
                    amounts={"amount":("0.40 mmol",q(0.40,"mmol"),"original_extraction",None),
                             "equivalents":("2.0 equiv",q(2.0,"equiv"),"original_extraction",None)}, reference=f"m_{allene_label}"),
        participant("iodine", "m_i2", "reagent", conditions_evidence,
                    amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),
                             "equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)}, reference=f"m_{allene_label}"),
        participant("base", "m_tbuoli", "reagent", conditions_evidence,
                    amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),
                             "equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)}, reference=f"m_{allene_label}"),
        participant("solvent", "m_dcm", "solvent", conditions_evidence,
                    amounts={"volume":("3.0 mL CH2Cl2",q(3.0,"mL"),"original_extraction",None)}),
    ]
    positive = percent > 0
    if positive:
        r["participants"].append(participant("product", f"m_{label}", "product", extra,
                                             product_role="desired product"))
        r["participants"].extend([
            participant("petroleum_ether", "m_petroleum_ether", "purification_agent", ev(5), stage="purification"),
            participant("ethyl_acetate", "m_etoac", "purification_agent", ev(5), stage="purification"),
        ])
    r["reaction_scale"] = F("0.20 mmol of 1", {"reference_molecule_id":f"m_{allene_label}","amount":q(0.20,"mmol")}, conditions_evidence)
    set_solvent(r, "CH2Cl2", 3.0, conditions_evidence)
    add_step(r, "addition_and_reaction", ["allene","acid","base","solvent","iodine"],
             "I2 was added to the solution of allenylphosphonate, carboxylic acid, and t-BuOLi; mixture stirred",
             "r.t.", q(hours,"h"), "N2", conditions_evidence)
    if positive:
        add_step(r, "solvent_removal", ["product"], "solvent removed in vacuo", None, None, None, ev(5), stage="workup")
        add_step(r, "column", ["product","petroleum_ether","ethyl_acetate"],
                 "silica gel column; petroleum ether/EtOAc = 5:1", None, None, None, ev(5), stage="purification")
    r["procedure"]["original_text"] = F(
        "Scheme 2 footnote plus the page-5 General Procedure for carboxylic acids",
        {"scheme":"Scheme 2","general_procedure":"page 5"}, evidence)
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure for the Electrophilic Addition of Ferrocene-Containing Allenylphosphonates with Carboxylic Acids",
        "page-5 carboxylic-acid general procedure", ev(5))
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = F(
        "I2 added to a solution of 1, 2, and t-BuOLi in CH2Cl2",
        ["allenylphosphonate","carboxylic acid","t-BuOLi","CH2Cl2","then I2"], ev(5))
    r["procedure"]["pretreatment"] = missing(evidence)
    if positive:
        add_yield(r, f"yield_{label}", f"m_{label}", percent,
                  "isolated product after chromatography", extra, product_mass=mass,
                  product_mass_evidence=ev(scope3_mass_page[label]) if mass is not None else None)
        add_selectivity(r, f"selectivity_{label}", f"m_{label}",
                        f"E-{label} depicted; numerical E/Z ratio unreported", ev(3))
        set_classification(r, evidence)
        unresolved = ["Numerical E/Z ratio and the method used to establish stereoselectivity are unreported in the bounded pages."]
        if label == "3k":
            unresolved.append("Long-chain acid/product carbon count is not confidently legible; structures remain unresolved.")
    else:
        add_no_reaction(r, f"no_reaction_{label}", "0%", ev(3), f"E-{label}")
        set_classification(r, evidence, "no_reaction")
        unresolved = ["No desired product detected/0%; detection method and limit are unreported, so no product molecule is fabricated."]
    finish_review(r, unresolved)
    reactions.append(r)


# Independent scale-up experiment on page 5.
evidence = ev(5)
r = st.new_reaction(f"{DOC}:large_scale_3a", DOC, evidence)
set_reaction_basics(r, ["scale_up"], "Large-Scale Preparation of 3a",
                    "Large-Scale Preparation of 3a",
                    "iodoesterification of allenylphosphonate 1a with benzoic acid",
                    evidence, [scheme2_id])
r["participants"] = [
    participant("allene","m_1a","reactant",evidence,amounts={"amount":("2.5 mmol",q(2.5,"mmol"),"original_extraction",None),"mass":("1.02 g",q(1.02,"g"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference="m_1a"),
    participant("acid","m_2a","reactant",evidence,amounts={"amount":("5.0 mmol",q(5.0,"mmol"),"original_extraction",None),"mass":("0.61 g",q(0.61,"g"),"original_extraction",None),"equivalents":("2.0 equiv",q(2.0,"equiv"),"original_extraction",None)},reference="m_1a"),
    participant("base","m_tbuoli","reagent",evidence,amounts={"amount":("2.5 mmol",q(2.5,"mmol"),"original_extraction",None),"mass":("0.20 g",q(0.20,"g"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference="m_1a"),
    participant("iodine","m_i2","reagent",evidence,amounts={"amount":("2.5 mmol",q(2.5,"mmol"),"original_extraction",None),"mass":("0.63 g",q(0.63,"g"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference="m_1a"),
    participant("solvent","m_dcm","solvent",evidence,amounts={"volume":("25.0 mL",q(25.0,"mL"),"original_extraction",None)}),
    participant("petroleum_ether","m_petroleum_ether","purification_agent",evidence,stage="purification"),
    participant("ethyl_acetate","m_etoac","purification_agent",evidence,stage="purification"),
    participant("product","m_3a","product",evidence,product_role="desired product"),
]
r["reaction_scale"] = F("1a: 2.5 mmol (1.02 g)", {"reference_molecule_id":"m_1a","amount":q(2.5,"mmol"),"mass":q(1.02,"g")}, evidence)
set_solvent(r,"CH2Cl2",25.0,evidence)
add_step(r,"addition_and_reaction",["allene","acid","base","solvent","iodine"],"I2 added to the solution; mixture stirred","r.t.",q(2,"h"),"N2",evidence)
add_step(r,"solvent_removal",["product"],"solvent removed in vacuo",None,None,None,evidence,stage="workup")
add_step(r,"column",["product","petroleum_ether","ethyl_acetate"],"silica gel column; petroleum ether/EtOAc = 5:1",None,None,None,evidence,stage="purification")
r["procedure"]["original_text"] = F("Large-Scale Preparation of 3a paragraph", "full page-5 scale-up paragraph", evidence)
r["procedure"]["general_procedure_reference"] = missing(evidence)
r["procedure"]["modifications"] = F("2.5 mmol scale", {"scale_up":True}, evidence)
r["procedure"]["addition_order"] = F("I2 added after 1a, 2a, and t-BuOLi were dissolved in CH2Cl2", ["1a","2a","t-BuOLi","CH2Cl2","I2"], evidence)
r["procedure"]["pretreatment"] = missing(evidence)
add_yield(r,"yield_3a","m_3a",66,"isolated after silica gel column chromatography",evidence)
pm = st.new_outcome("product_mass_3a","product_amount",["m_3a"],evidence)
pm["value"] = F("1.08 g",q(1.08,"g"),evidence); pm["product_role"] = F("desired product","desired product",evidence)
pm["yield_type"] = not_applicable(evidence); pm["measurement_basis"] = F("isolated product mass","isolated",evidence)
pm["scope"] = F("scale-up isolated mass",{"product":"m_3a"},evidence); r["outcomes"].append(pm)
add_selectivity(r,"selectivity_3a","m_3a","E-3a; numerical E/Z ratio unreported",evidence)
set_classification(r,evidence)
finish_review(r,["Numerical E/Z ratio and stereochemical measurement method are unreported."])
reactions.append(r)


# Scheme 3 alcohol scope.
scope5 = [
    ("5a","1a","4a",1,62,70), ("5b","1a","4b",1,57,74),
    ("5c","1a","4c",1,57,86), ("5d","1a","4d",1,63,75),
    ("5e","1a","4e",2,59,77), ("5f","1a","4f",1,50,62),
    ("5g","1a","4g",2.5,14,17), ("5h","1l","4a",1,66,74),
    ("5i","1q","4a",1,58,63), ("5j","1r","4a",1,54,58),
]
scope5_mass_page = {**{label:7 for label in ("5a","5b","5c","5d")},
                    **{label:8 for label in ("5e","5f","5g","5h","5i","5j")}}
for label, allene_label, alcohol_label, hours, percent, mass in scope5:
    extra = alcohol_product_evidence[label]
    evidence = list(dict.fromkeys(ev(4,7) + extra))
    conditions_evidence = ev(4,7)
    r = st.new_reaction(f"{DOC}:scope_{label}",DOC,evidence)
    set_reaction_basics(r,["substrate_scope"],f"Scheme 3, {label}",
                        "Highly Regio- and Stereoselective Iodoalkoxylation of Allenylphosphonates with Various Alcohols",
                        "iodoalkoxylation of an allenylphosphonate",evidence,[scheme3_id])
    r["participants"] = [
        participant("allene",f"m_{allene_label}","reactant",conditions_evidence,amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference=f"m_{allene_label}"),
        participant("alcohol",f"m_{alcohol_label}","reactant",conditions_evidence,amounts={"amount":("0.40 mmol",q(0.40,"mmol"),"original_extraction",None),"equivalents":("2.0 equiv",q(2.0,"equiv"),"original_extraction",None)},reference=f"m_{allene_label}"),
        participant("iodine","m_i2","reagent",conditions_evidence,amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference=f"m_{allene_label}"),
        participant("base","m_tbuoli","reagent",conditions_evidence,amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference=f"m_{allene_label}"),
        participant("solvent","m_dcm","solvent",conditions_evidence,amounts={"volume":("3.0 mL",q(3.0,"mL"),"original_extraction",None)}),
        participant("petroleum_ether","m_petroleum_ether","purification_agent",ev(7),stage="purification"),
        participant("ethyl_acetate","m_etoac","purification_agent",ev(7),stage="purification"),
        participant("product",f"m_{label}","product",extra,product_role="desired product"),
    ]
    r["reaction_scale"] = F("0.20 mmol of 1",{"reference_molecule_id":f"m_{allene_label}","amount":q(0.20,"mmol")},conditions_evidence)
    set_solvent(r,"CH2Cl2",3.0,conditions_evidence)
    add_step(r,"addition_and_reaction",["allene","alcohol","base","solvent","iodine"],
             "I2 was added to the solution of allenylphosphonate, alcohol, and t-BuOLi; mixture stirred",
             "r.t.",q(hours,"h"),"N2",conditions_evidence)
    add_step(r,"solvent_removal",["product"],"solvent removed in vacuo",None,None,None,ev(7),stage="workup")
    add_step(r,"column",["product","petroleum_ether","ethyl_acetate"],
             "silica gel column; petroleum ether/EtOAc = 5:1",None,None,None,ev(7),stage="purification")
    r["procedure"]["original_text"] = F(
        "Scheme 3 row/footnote plus the page-7 General Procedure for alcohols",
        {"scheme":"Scheme 3","product":label,"general_procedure":"page 7"},evidence)
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure for the Electrophilic Addition of Ferrocene-Containing Allenylphosphonates with Alcohols",
        "page-7 alcohol general procedure",ev(7))
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = F(
        "I2 added to a solution of 1, 4, and t-BuOLi in CH2Cl2",
        ["allenylphosphonate","alcohol","t-BuOLi","CH2Cl2","then I2"],ev(7))
    r["procedure"]["pretreatment"] = missing(evidence)
    add_yield(r,f"yield_{label}",f"m_{label}",percent,
              "isolated product after chromatography",extra,product_mass=mass,
              product_mass_evidence=ev(scope5_mass_page[label]))
    wording = "only E-isomer" if label == "5a" else "single stereoisomer" if label in {"5h","5i","5j"} else "E product depicted"
    add_selectivity(r,f"selectivity_{label}",f"m_{label}",f"{wording}; numerical ratio unreported",evidence)
    set_classification(r,evidence)
    finish_review(r,["Numerical E/Z ratio and stereochemical measurement method are unreported."])
    reactions.append(r)


# Scheme 4 phenol scope.
scope7 = [
    ("7a","6a",2,61,75), ("7b","6b",2,59,76), ("7c","6c",2,23,31),
    ("7d","6d",3,0,None), ("7e","6e",2,50,69), ("7f","6f",2,52,66),
    ("7g","6g",2,48,64), ("7h","6h",2,0,None), ("7i","6i",2,40,55),
]
scope7_mass_page = {**{label:8 for label in ("7a","7b","7c","7e")},
                    **{label:9 for label in ("7f","7g","7i")}}
for label, phenol_label, hours, percent, mass in scope7:
    extra = phenol_product_evidence.get(label, ev(5,8))
    evidence=list(dict.fromkeys(ev(5,8) + extra))
    conditions_evidence = ev(5,8)
    r=st.new_reaction(f"{DOC}:scope_{label}",DOC,evidence)
    set_reaction_basics(r,["substrate_scope"],f"Scheme 4, {label}",
                        "Synthesis of Tetrasubstituted Alkenylphosphonates by the Reaction of Allenylphosphonate 1a with Various Phenols",
                        "iodoaryloxylation of allenylphosphonate 1a",evidence,[scheme4_id])
    r["participants"]=[
        participant("allene","m_1a","reactant",conditions_evidence,amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference="m_1a"),
        participant("phenol",f"m_{phenol_label}","reactant",conditions_evidence,amounts={"amount":("0.40 mmol",q(0.40,"mmol"),"original_extraction",None),"equivalents":("2.0 equiv",q(2.0,"equiv"),"original_extraction",None)},reference="m_1a"),
        participant("iodine","m_i2","reagent",conditions_evidence,amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference="m_1a"),
        participant("base","m_tbuoli","reagent",conditions_evidence,amounts={"amount":("0.20 mmol",q(0.20,"mmol"),"original_extraction",None),"equivalents":("1.0 equiv",q(1.0,"equiv"),"original_extraction",None)},reference="m_1a"),
        participant("solvent","m_dcm","solvent",conditions_evidence,amounts={"volume":("3.0 mL",q(3.0,"mL"),"original_extraction",None)}),
    ]
    positive=percent>0
    if positive:
        r["participants"].extend([
            participant("petroleum_ether","m_petroleum_ether","purification_agent",ev(8),stage="purification"),
            participant("ethyl_acetate","m_etoac","purification_agent",ev(8),stage="purification"),
            participant("product",f"m_{label}","product",extra,product_role="desired product"),
        ])
    r["reaction_scale"]=F("0.20 mmol of 1a",{"reference_molecule_id":"m_1a","amount":q(0.20,"mmol")},conditions_evidence)
    set_solvent(r,"CH2Cl2",3.0,conditions_evidence)
    add_step(r,"addition_and_reaction",["allene","phenol","base","solvent","iodine"],
             "I2 was added to the solution of allenylphosphonate, phenol, and t-BuOLi; mixture stirred",
             "r.t.",q(hours,"h"),"N2",conditions_evidence)
    if positive:
        add_step(r,"solvent_removal",["product"],"solvent removed in vacuo",None,None,None,ev(8),stage="workup")
        add_step(r,"column",["product","petroleum_ether","ethyl_acetate"],
                 "silica gel column; petroleum ether/EtOAc = 5:1",None,None,None,ev(8),stage="purification")
    r["procedure"]["original_text"]=F(
        "Scheme 4 row/footnote plus the page-8 General Procedure for phenols",
        {"scheme":"Scheme 4","product":label,"general_procedure":"page 8"},evidence)
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure for the Electrophilic Addition of Ferrocene-Containing Allenylphosphonates with Phenols",
        "page-8 phenol general procedure",ev(8))
    r["procedure"]["modifications"] = missing(evidence)
    r["procedure"]["addition_order"] = F(
        "I2 added to a solution of 1, 6, and t-BuOLi in CH2Cl2",
        ["allenylphosphonate","phenol","t-BuOLi","CH2Cl2","then I2"],ev(8))
    r["procedure"]["pretreatment"] = missing(evidence)
    if positive:
        add_yield(r,f"yield_{label}",f"m_{label}",percent,
                  "isolated product after chromatography",extra,product_mass=mass,
                  product_mass_evidence=ev(scope7_mass_page[label]))
        add_selectivity(r,f"selectivity_{label}",f"m_{label}","complete stereoselectivity to the E product; numerical ratio unreported",evidence)
        set_classification(r,evidence)
        unresolved=["Numerical E/Z ratio and stereochemical measurement method are unreported."]
    else:
        add_no_reaction(r,f"no_reaction_{label}","0%",evidence,f"E-{label}")
        set_classification(r,evidence,"no_reaction")
        unresolved=["No desired product detected/0%; detection method and limit are unreported, so no product molecule is fabricated."]
    finish_review(r,unresolved)
    reactions.append(r)


# Scheme 5 is an author-proposed mechanism, not a mechanistic experiment.
evidence=ev(4,5)
r=st.new_reaction(f"{DOC}:scheme5_proposed_mechanism",DOC,evidence)
r["record_type"]="proposed_mechanism"
r["experiment_categories"]=F(["other"],["other"],evidence)
r["location"]=F("Scheme 5 and accompanying Results and Discussion paragraph","Scheme 5",evidence)
r["general_scheme_ids"]=[scheme1_id]
r["reaction_name"]=F("Rationalization of Regio- and Stereoselectivity","Rationalization of Regio- and Stereoselectivity",evidence)
r["transformation_type"]=F(None,"proposed iodonium/neighboring-group-participation mechanism",evidence,"inference",derivation="Classified as proposed mechanism from the source heading and prose wording 'may be interpreted'.")
r["participants"]=[participant("allene","m_1a","reactant",evidence),participant("iodine","m_i2","reagent",evidence),participant("oxygen_nucleophile","m_generic_oxygen_nucleophile","reactant",evidence),participant("product","m_generic_E_product","product",evidence,product_role="proposed E product")]
for sid,op in [("iodonium","electrophilic addition of I2 gives proposed iodonium intermediate A"),("neighboring_group","phosphate oxygen participation gives proposed five-membered intermediate B"),("nucleophilic_attack","oxygen nucleophile attack and C-O cleavage gives the E-isomer")]:
    add_step(r,sid,["allene","iodine","oxygen_nucleophile","product"],op,None,None,None,evidence)
r["procedure"]["original_text"]=F("Initially A; next B; finally oxygen nucleophile attack affords E-isomer","author-proposed three-stage mechanism",evidence)
for key in ("general_procedure_reference","modifications","addition_order","pretreatment"):
    r["procedure"][key]=missing(evidence)
set_classification(r,evidence,"mechanism")
finish_review(r,["This is an author-proposed mechanism, not a directly observed mechanistic experiment.","Intermediates A and B are not assigned experimental yields or complete atom mappings."])
reactions.append(r)


source_context = [
    {"context_id":f"{DOC}:page1_intro_context","document_id":DOC,"evidence_ids":ev(1),
     "kind":"introduction_and_abstract","note":"Title, abstract, TOC graphic, and Figure 1 bioactive examples are context, not separate experiments."},
    {"context_id":f"{DOC}:pages5_6_characterization","document_id":DOC,"evidence_ids":ev(5,6),
     "kind":"characterization_context","note":"General instrumentation plus IR/NMR/HRMS/melting-point/Rf data are not separate synthesis experiments; visible preparation masses/yields were merged into Scheme-2 records."},
    {"context_id":f"{DOC}:pages7_9_characterization","document_id":DOC,"evidence_ids":ev(7,8,9),
     "kind":"characterization_context","note":"IR/NMR/HRMS/melting-point/Rf peak lists are characterization, not separate synthesis experiments; preparation masses/yields and the alcohol/phenol general procedures were merged into their existing scope records."},
    {"context_id":f"{DOC}:pages9_10_back_matter","document_id":DOC,"evidence_ids":ev(9,10),
     "kind":"publication_back_matter","note":"Associated-content notice, author information, acknowledgments, and references are publication context or cited prior work, not current-paper experiments."},
    {"context_id":f"{DOC}:si_unavailable","document_id":DOC,"evidence_ids":ev(1),
     "kind":"source_limitation","note":"The article advertises Supporting Information, but source_mapping.json does not explicitly associate an SI file; none was used."},
]


dataset = json.loads((RUN / "dataset.json").read_text())
dataset["documents"][0]["title"] = TITLE
(RUN / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
bundle = json.loads((RUN / "submitted.json").read_text())
bundle["documents"] = copy.deepcopy(dataset["documents"])
bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = source_context
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_checks.json").write_text(json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/extraction_inventory.json").write_text(json.dumps({
    "document_id":DOC,
    "title":TITLE,
    "bounded_pages":[1,2,3,4,5,6,7,8,9,10],
    "current_batch_pages":[7,8,9,10],
    "record_counts":{
        "generic_schemes":sum(r["record_type"]=="generic_scheme" for r in reactions),
        "optimization_entries":13,
        "carboxylic_acid_scope":20,
        "scale_up":1,
        "alcohol_scope":10,
        "phenol_scope":9,
        "proposed_mechanisms":1,
        "total":len(reactions),
    },
    "molecule_count":len(molecules),
    "notes":[
        "All records are machine candidates and remain needs_review/partial.",
        "No SI file was associated.",
        "Characterization-only passages were retained as source_context, not reactions.",
        "Pages 7-9 extend existing scope records with reported product masses and explicit alcohol/phenol general procedures; page 10 contains references only.",
    ],
}, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules":len(molecules),"reactions":len(reactions),"stereo_checks":len(stereo_checks)},ensure_ascii=False))
