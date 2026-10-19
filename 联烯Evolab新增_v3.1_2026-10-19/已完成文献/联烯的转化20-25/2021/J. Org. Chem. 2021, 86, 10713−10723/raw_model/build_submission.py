#!/usr/bin/env python3
"""Build the page-1--11 v3.1 extraction bundle for doc_fdd3bd33ea5bf346.

This is a reproducible machine submission.  It deliberately remains
needs_review and does not claim chemical or semantic verification.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DATASET = json.loads((RUN / "dataset.json").read_text())
DOC = "doc_fdd3bd33ea5bf346"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def fe(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
       status=None, verification=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def missing(evidence=()):
    return fe(evidence=evidence)


def na(evidence=(), reason="Not applicable to this field."):
    return fe(evidence=evidence, method="not_applicable", status="not_applicable",
              verification="missing", derivation=reason)


def derived(normalized, evidence, reason, raw=None):
    return fe(raw, normalized, evidence, method="structure_derivation", derivation=reason)


def inferred(normalized, evidence, reason, raw=None):
    return fe(raw, normalized, evidence, method="inference", derivation=reason)


bundle = st.empty(DATASET["documents"])
molecules = {}


def molecule(mid, label, name, smiles=None, evidence=None, chirality="achiral",
             geometry="not_applicable", representation=None, axial=None,
             double_bonds=None, identity_note=None):
    evidence = list(evidence or ev(1))
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fe(label, label, evidence)
    m["names"]["original"] = fe(name, name, evidence)
    m["names"]["english"] = fe(name, name, evidence)
    m["names"]["abbreviations"] = missing(evidence)
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["external_properties"] = []
    m["identity_note"] = fe(identity_note, identity_note, evidence) if identity_note else missing(evidence)
    stereo = m["stereochemistry"]
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence,
            "Candidate SMILES transcribed from the named or drawn source structure; parsing is not source verification.",
            {"source_structure": name, "compound_label": label})
        diag = st.stereo_inspect(smiles)
        stereo["identity_key"] = derived(diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key calculated from the candidate SMILES.")
        stereo["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], evidence,
            "Canonical isomeric SMILES calculated from the candidate SMILES.")
        representation = representation or ("not_applicable" if chirality == "achiral" and geometry == "not_applicable" else "partial")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(evidence)
        stereo["identity_key"] = missing(evidence)
        stereo["isomeric_smiles"] = missing(evidence)
        representation = representation or "unresolved"
    m["generic_template"] = None
    stereo["chirality_status"] = inferred(chirality, evidence,
        "Status assigned conservatively from the displayed structure and absence/presence of explicit stereochemical annotation; unspecified is not treated as achiral.",
        {"displayed_label": label})
    stereo["geometry_status"] = inferred(geometry, evidence,
        "Status assigned from the displayed structure; unannotated configurable bonds remain unspecified.",
        {"displayed_label": label})
    stereo["representation_status"] = inferred(representation, evidence,
        "Compared the source stereochemical content with the candidate SMILES representation.")
    stereo["atom_centers"] = inferred([], evidence,
        "No atom-center assignment is asserted; an empty list does not convert unspecified stereochemistry to achiral.")
    stereo["double_bonds"] = inferred(double_bonds or [], evidence,
        "Double-bond assignments retained only when displayed or named explicitly.")
    stereo["axial"] = inferred(axial or [], evidence,
        "Allene axial configuration is retained as unassigned where the allene termini can support axial chirality.")
    stereo["mixture"] = missing(evidence)
    stereo["assignment"] = missing(evidence)
    molecules[mid] = m
    return mid


def unresolved_molecule(mid, label, name, evidence):
    return molecule(mid, label, name, None, evidence, chirality="unresolved",
                    geometry="unresolved", representation="unresolved")


# Common reagents, solvents, and workup materials. Structures are intentionally
# unresolved because the extraction target is source identity, not a database fill.
common = {
    "m_lta": ("LTA", "lead(IV) acetate"),
    "m_pida": ("PIDA", "phenyliodine diacetate"),
    "m_tfe": ("TFE", "2,2,2-trifluoroethanol"),
    "m_hfip": ("HFIP", "hexafluoroisopropanol"),
    "m_dcm": ("DCM", "dichloromethane"),
    "m_dce": ("DCE", "1,2-dichloroethane"),
    "m_benzene": ("benzene", "benzene"),
    "m_acoh": ("AcOH", "acetic acid"),
    "m_pb_obz4": ("Pb(OBz)4", "lead(IV) tetrabenzoate"),
    "m_pb_pfbz4": ("Pb(OCO-4-F-C6H4)4", "lead(IV) tetra(4-fluorobenzoate)"),
    "m_pb_pmbz4": ("Pb(OCO-4-MeO-C6H4)4", "lead(IV) tetra(4-methoxybenzoate)"),
    "m_pb_ptfbz4": ("Pb(OCO-4-CF3-C6H4)4", "lead(IV) tetra(4-trifluoromethylbenzoate)"),
    "m_pb_cinn4": ("Pb(OCOCH=CHPh)4", "lead(IV) tetracinnamate"),
    "m_pyridine": ("pyridine", "pyridine"),
    "m_methyl_chloroformate": ("methyl chloroformate", "methyl chloroformate"),
    "m_water": ("water", "water"),
    "m_brine": ("brine", "brine"),
    "m_na2so4": ("Na2SO4", "sodium sulfate"),
    "m_etoac": ("EtOAc", "ethyl acetate"),
    "m_hexanes": ("hexanes", "hexanes"),
    "m_pdcl2pph3": ("PdCl2(PPh3)2", "bis(triphenylphosphine)palladium(II) dichloride"),
    "m_cui": ("CuI", "copper(I) iodide"),
    "m_thf": ("THF", "tetrahydrofuran"),
    "m_et3n": ("NEt3", "triethylamine"),
    "m_nh4cl": ("NH4Cl", "ammonium chloride"),
}
for mid, (label, name) in common.items():
    unresolved_molecule(mid, label, name, ev(2) if mid.startswith("m_pb") or mid in {"m_lta", "m_pida", "m_tfe", "m_hfip", "m_dcm", "m_dce", "m_benzene", "m_acoh"} else ev(4) + ev(5))

# Generic identities used only by general-scheme records.
for mid, label, name, evidence in [
    ("m_generic_o_allenyl_aniline", "generic o-allenyl aniline", "carbamate-protected o-allenyl aniline", ev(1)),
    ("m_generic_indole_product", "generic indole product", "2-(acyloxyalkyl)indole carbamate", ev(1)),
    ("m_generic_iodoaniline", "2-iodoaniline derivative", "protected or substituted 2-iodoaniline derivative", ev(4)),
    ("m_generic_carbamate", "protected 2-iodoaniline", "carbamate-protected 2-iodoaniline derivative", ev(4)),
    ("m_generic_alkyne", "alkyne", "terminal propargylic alcohol reagent", ev(5)),
    ("m_generic_alkynyl_product", "o-alkynyl aniline", "carbamate-protected o-alkynyl aniline", ev(5)),
]:
    unresolved_molecule(mid, label, name, evidence)

# Core allene substrates and indole products shown on pages 2--4.
subs = {
    "m_1": ("1", "methyl [4-bromo-2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccc(Br)cc1C=C=C", 2),
    "m_s27": ("S27", "methyl [2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccccc1C=C=C", 2),
    "m_s28": ("S28", "methyl [4-fluoro-2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccc(F)cc1C=C=C", 2),
    "m_s29": ("S29", "methyl [4-chloro-2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccc(Cl)cc1C=C=C", 2),
    "m_s30": ("S30", "methyl [4-methyl-2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccc(C)cc1C=C=C", 2),
    "m_s31": ("S31", "methyl [4-(tert-butyl)-2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccc(C(C)(C)C)cc1C=C=C", 2),
    "m_s32": ("S32", "methyl [4-methoxy-2-(propa-1,2-dien-1-yl)phenyl]carbamate", "COC(=O)Nc1ccc(OC)cc1C=C=C", 2),
    "m_s33": ("S33", "methyl [2-(propa-1,2-dien-1-yl)-4-(trifluoromethyl)phenyl]carbamate", "COC(=O)Nc1ccc(C(F)(F)F)cc1C=C=C", 2),
    "m_s34": ("S34", "tert-butyl [2-(propa-1,2-dien-1-yl)phenyl]carbamate", "CC(C)(C)OC(=O)Nc1ccccc1C=C=C", 2),
    "m_s35": ("S35", "benzyl [2-(propa-1,2-dien-1-yl)phenyl]carbamate", "O=C(OCc1ccccc1)Nc1ccccc1C=C=C", 2),
    "m_s36": ("S36", "2,2,2-trichloroethyl [2-(propa-1,2-dien-1-yl)phenyl]carbamate", "O=C(OCC(Cl)(Cl)Cl)Nc1ccccc1C=C=C", 2),
}
for mid, (label, name, smi, page) in subs.items():
    molecule(mid, label, name, smi, ev(page), chirality="achiral")

prod_smiles = {
    2: "COC(=O)n1c(COC(C)=O)cc2cc(Br)ccc21",
    3: "COC(=O)n1c(COC(C)=O)cc2cc(Cl)ccc21",
    4: "COC(=O)n1c(COC(C)=O)cc2cc(F)ccc21",
    5: "COC(=O)n1c(COC(C)=O)cc2cc(C)ccc21",
    6: "COC(=O)n1c(COC(C)=O)cc2cc(C(C)(C)C)ccc21",
    7: "COC(=O)n1c(COC(C)=O)cc2cc(OC)ccc21",
    8: "COC(=O)n1c(COC(C)=O)cc2cc(C(F)(F)F)ccc21",
    9: "COC(=O)n1c(COC(C)=O)cc2ccccc21",
    10: "CC(C)(C)OC(=O)n1c(COC(C)=O)cc2ccccc21",
    11: "O=C(OCc1ccccc1)n1c(COC(C)=O)cc2ccccc21",
    12: "O=C(OCC(Cl)(Cl)Cl)n1c(COC(C)=O)cc2ccccc21",
    13: "COC(=O)n1c(COC(=O)c2ccccc2)cc2ccccc21",
    14: "COC(=O)n1c(COC(=O)c2ccc(F)cc2)cc2ccccc21",
    15: "COC(=O)n1c(COC(=O)c2ccc(OC)cc2)cc2ccccc21",
    16: "COC(=O)n1c(COC(=O)c2ccc(C(F)(F)F)cc2)cc2ccccc21",
    17: "COC(=O)n1c(COC(=O)/C=C/c2ccccc2)cc2ccccc21",
}
for n, smi in prod_smiles.items():
    name = f"source-drawn indole product {n}"
    geom = "specified" if n == 17 else "not_applicable"
    rep = "encoded" if n == 17 else "not_applicable"
    db = [{"group": "cinnamoyl C=C", "configuration": "E (drawn trans)"}] if n == 17 else []
    molecule(f"m_{n}", str(n), name, smi, ev(2 if n <= 12 else 3),
             chirality="achiral", geometry=geom, representation=rep, double_bonds=db)

# Terminally substituted allene scope.
molecule("m_18", "18", "methyl [4-bromo-2-(buta-1,2-dien-1-yl)phenyl]carbamate",
         "COC(=O)Nc1ccc(Br)cc1C=C=CC", ev(3), chirality="unspecified",
         representation="not_encoded", axial=[{"axis": "allene", "configuration": "unassigned"}])
molecule("m_19", "19", "methyl [2-(buta-1,2-dien-1-yl)phenyl]carbamate",
         "COC(=O)Nc1ccccc1C=C=CC", ev(3), chirality="unspecified",
         representation="not_encoded", axial=[{"axis": "allene", "configuration": "unassigned"}])
molecule("m_20", "20", "methyl [4-bromo-2-(3-methylbuta-1,2-dien-1-yl)phenyl]carbamate",
         "COC(=O)Nc1ccc(Br)cc1C=C=C(C)C", ev(3), chirality="achiral")
molecule("m_21", "21", "methyl [2-(3-methylbuta-1,2-dien-1-yl)phenyl]carbamate",
         "COC(=O)Nc1ccccc1C=C=C(C)C", ev(3), chirality="achiral")
for n, smi, ch in [
    (22, "COC(=O)n1c(C(C)OC(C)=O)cc2cc(Br)ccc21", "unspecified"),
    (23, "COC(=O)n1c(C(C)OC(C)=O)cc2ccccc21", "unspecified"),
    (24, "COC(=O)n1c(C(C)(C)OC(C)=O)cc2cc(Br)ccc21", "achiral"),
    (25, "COC(=O)n1c(C(C)(C)OC(C)=O)cc2ccccc21", "achiral"),
]:
    molecule(f"m_{n}", str(n), f"source-drawn indole product {n}", smi, ev(3),
             chirality=ch, representation="partial" if ch == "unspecified" else "not_applicable")

# Control substrates and products.
molecule("m_26", "26", "methyl o-tolylcarbamate", "COC(=O)Nc1ccccc1C", ev(3))
unresolved_molecule("m_27", "27", "6-[(methoxycarbonyl)amino]-1-methylcyclohexa-2,4-dien-1-yl acetate", ev(3))
molecule("m_28", "28", "phenylallene control", "c1ccccc1C=C=C", ev(3))
molecule("m_29", "29", "methyl (2-ethynylphenyl)carbamate", "COC(=O)Nc1ccccc1C#C", ev(3))
molecule("m_30", "30", "methyl 3-acetoxy-1H-indole-1-carboxylate", "COC(=O)n1cc(OC(C)=O)c2ccccc21", ev(3))
molecule("m_31", "31", "methyl (2-allylphenyl)carbamate", "COC(=O)Nc1ccccc1CC=C", ev(3))
molecule("m_32", "32", "methyl (2-vinylphenyl)carbamate", "COC(=O)Nc1ccccc1C=C", ev(3), geometry="not_applicable", representation="not_applicable")
molecule("m_33", "33", "2-{2-[(methoxycarbonyl)amino]phenyl}ethane-1,1-diyl diacetate", "COC(=O)Nc1ccccc1CC(OC(C)=O)OC(C)=O", ev(3))
molecule("m_34", "34", "methyl 1H-indole-1-carboxylate", "COC(=O)n1ccc2ccccc21", ev(3))
molecule("m_35", "35", "methyl [2-(propa-1,2-dien-1-yl-1-d)phenyl]carbamate", "COC(=O)Nc1ccccc1[C]([2H])=C=C", ev(3)+ev(4))
molecule("m_36", "36", "methyl 2-(acetoxymethyl)-1H-indole-1-carboxylate-3-d", "COC(=O)n1c(COC(C)=O)c([2H])c2ccccc21", ev(3)+ev(4))

# General-procedure A starting materials and products S1--S11.
prep_a = {
    "m_a1": ("2-iodoaniline", "2-iodoaniline", "Nc1ccccc1I"),
    "m_a2": ("4-fluoro-2-iodoaniline", "4-fluoro-2-iodoaniline", "Nc1ccc(F)cc1I"),
    "m_a3": ("4-chloro-2-iodoaniline", "4-chloro-2-iodoaniline", "Nc1ccc(Cl)cc1I"),
    "m_a4": ("4-bromo-2-iodoaniline", "4-bromo-2-iodoaniline", "Nc1ccc(Br)cc1I"),
    "m_a5": ("2-iodo-4-methylaniline", "2-iodo-4-methylaniline", "Nc1ccc(C)cc1I"),
    "m_a7": ("2-iodo-4-methoxyaniline", "2-iodo-4-methoxyaniline", "Nc1ccc(OC)cc1I"),
    "m_a8": ("2-iodo-4-(trifluoromethyl)aniline", "2-iodo-4-(trifluoromethyl)aniline", "Nc1ccc(C(F)(F)F)cc1I"),
    "m_2iodobenzoic_acid": ("2-iodobenzoic acid", "2-iodobenzoic acid", "O=C(O)c1ccccc1I"),
}
for mid, (label, name, smi) in prep_a.items():
    molecule(mid, label, name, smi, ev(4)+ev(5))

s_products = {
    1: "COC(=O)Nc1ccccc1I",
    2: "COC(=O)Nc1ccc(F)cc1I",
    3: "COC(=O)Nc1ccc(Cl)cc1I",
    4: "COC(=O)Nc1ccc(Br)cc1I",
    5: "COC(=O)Nc1ccc(C)cc1I",
    6: "COC(=O)Nc1ccc(C(C)(C)C)cc1I",
    7: "COC(=O)Nc1ccc(OC)cc1I",
    8: "COC(=O)Nc1ccc(C(F)(F)F)cc1I",
    9: "CC(C)(C)OC(=O)Nc1ccccc1I",
    10: "O=C(OCc1ccccc1)Nc1ccccc1I",
    11: "O=C(OCC(Cl)(Cl)Cl)Nc1ccccc1I",
}
s_names = {
    1:"methyl (2-iodophenyl)carbamate",2:"methyl (4-fluoro-2-iodophenyl)carbamate",
    3:"methyl (4-chloro-2-iodophenyl)carbamate",4:"methyl (4-bromo-2-iodophenyl)carbamate",
    5:"methyl (2-iodo-4-methylphenyl)carbamate",6:"methyl [4-(tert-butyl)-2-iodophenyl]carbamate",
    7:"methyl (2-iodo-4-methoxyphenyl)carbamate",8:"methyl [2-iodo-4-(trifluoromethyl)phenyl]carbamate",
    9:"tert-butyl (2-iodophenyl)carbamate",10:"benzyl (2-iodophenyl)carbamate",
    11:"2,2,2-trichloroethyl (2-iodophenyl)carbamate",
}
for n, smi in s_products.items():
    evidence = ev(4)+ev(5) if n <= 2 else ev(5)
    molecule(f"m_s{n}", f"S{n}", s_names[n], smi, evidence)

# General-procedure B products S12--S26 and alkyne reagents.
molecule("m_propargyl_alcohol", "propargyl alcohol", "propargyl alcohol", "C#CCO", ev(5)+ev(6))
molecule("m_but3yn2ol", "but-3-yn-2-ol", "but-3-yn-2-ol", "C#CC(O)C", ev(6), chirality="unspecified", representation="partial")
molecule("m_2methylbut3yn2ol", "2-methylbut-3-yn-2-ol", "2-methylbut-3-yn-2-ol", "C#CC(O)(C)C", ev(6))
s12_26 = {
12:("methyl [2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccccc1C#CCO","achiral"),
13:("methyl [4-fluoro-2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccc(F)cc1C#CCO","achiral"),
14:("methyl [4-chloro-2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccc(Cl)cc1C#CCO","achiral"),
15:("methyl [4-bromo-2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccc(Br)cc1C#CCO","achiral"),
16:("methyl [2-(3-hydroxyprop-1-yn-1-yl)-4-methylphenyl]carbamate","COC(=O)Nc1ccc(C)cc1C#CCO","achiral"),
17:("methyl [4-(tert-butyl)-2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccc(C(C)(C)C)cc1C#CCO","achiral"),
18:("methyl [2-(3-hydroxyprop-1-yn-1-yl)-4-methoxyphenyl]carbamate","COC(=O)Nc1ccc(OC)cc1C#CCO","achiral"),
19:("methyl [2-(3-hydroxyprop-1-yn-1-yl)-4-(trifluoromethyl)phenyl]carbamate","COC(=O)Nc1ccc(C(F)(F)F)cc1C#CCO","achiral"),
20:("tert-butyl [2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","CC(C)(C)OC(=O)Nc1ccccc1C#CCO","achiral"),
21:("benzyl [2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","O=C(OCc1ccccc1)Nc1ccccc1C#CCO","achiral"),
22:("2,2,2-trichloroethyl [2-(3-hydroxyprop-1-yn-1-yl)phenyl]carbamate","O=C(OCC(Cl)(Cl)Cl)Nc1ccccc1C#CCO","achiral"),
23:("methyl [2-(3-hydroxybut-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccccc1C#CC(O)C","unspecified"),
24:("methyl [2-(3-hydroxy-3-methylbut-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccccc1C#CC(O)(C)C","achiral"),
25:("methyl [4-bromo-2-(3-hydroxybut-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccc(Br)cc1C#CC(O)C","unspecified"),
26:("methyl [4-bromo-2-(3-hydroxy-3-methylbut-1-yn-1-yl)phenyl]carbamate","COC(=O)Nc1ccc(Br)cc1C#CC(O)(C)C","achiral"),
}
for n,(name,smi,ch) in s12_26.items():
    evidence = ev(5)+ev(6) if n == 14 else ev(5) if n <= 13 else ev(6)
    molecule(f"m_s{n}", f"S{n}", name, smi, evidence, chirality=ch,
             representation="partial" if ch == "unspecified" else "not_applicable")


def part(pid, mid, role, evidence, stage="reaction", eq=None, amount=None,
         volume=None, mass=None, mol_percent=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role:
        p["product_role"] = fe(product_role, product_role, evidence)
    else:
        p["product_role"] = na(evidence)
    vals = {"equivalents": eq, "amount": amount, "volume": volume,
            "mass": mass, "mol_percent": mol_percent}
    for key, val in vals.items():
        if val is not None:
            p["charges"][key] = fe(str(val), val, evidence)
    if eq is not None or mol_percent is not None:
        p["equivalence_reference"] = inferred(None, evidence, "Set by caller after all participants are assembled.")
    else:
        p["equivalence_reference"] = na(evidence)
    return p


def set_refs(parts, reference_mid):
    for p in parts:
        if p["charges"]["equivalents"]["normalized_value"] is not None or p["charges"]["mol_percent"]["normalized_value"] is not None:
            p["equivalence_reference"] = fe(reference_mid, reference_mid, p["charges"]["equivalents"]["evidence_ids"] or p["charges"]["mol_percent"]["evidence_ids"])


def step(sid, evidence, participants, operation, temperature=None, time=None,
         atmosphere=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage)
    s["participant_ids"] = list(participants)
    s["operation"] = fe(operation, operation, evidence)
    s["temperature"] = fe(str(temperature), temperature, evidence) if temperature is not None else missing(evidence)
    s["time"] = fe(str(time), time, evidence) if time is not None else missing(evidence)
    s["atmosphere"] = fe(atmosphere, atmosphere, evidence) if atmosphere is not None else missing(evidence)
    s["pressure"] = missing(evidence)
    return s


def classification(evidence, bonds=(), coverage="partial", limitation=None):
    if coverage in ("unresolved", "not_applicable"):
        labels = [] if coverage == "not_applicable" else None
        return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(evidence),
                "bond_changes": [],
                "labels": inferred(labels, evidence, "No observed transformation for this record." if labels == [] else "Insufficient atom-level source mapping."),
                "coverage": inferred(coverage, evidence, "Classification scope set conservatively."),
                "limitations": [limitation or "No complete atom mapping was asserted."], "model_votes": []}
    changes=[]; labs={}
    for i,(a,b) in enumerate(bonds,1):
        cid=f"b{i}"
        change={"change_id":cid,"atoms":[{"atom_map":None,"element":a},{"atom_map":None,"element":b}],
                "before_order":0,"after_order":1,"change_type":"formed"}
        changes.append(derived(change,evidence,"Partial net bond change inferred from the source-drawn transformation."))
        lab=st.bond_label(a,b,1)
        labs.setdefault(lab["code"],dict(lab,supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    return {"standard_id":st.STANDARD,"mapped_reaction_smiles":missing(evidence),
            "bond_changes":changes,
            "labels":derived(list(labs.values()),evidence,"Labels computed from the explicitly listed partial formed-bond changes."),
            "coverage":inferred(coverage,evidence,"Only source-supported net formed bonds are included; complete atom mapping is unavailable."),
            "limitations":[limitation or "Partial connectivity classification; hydrogen transfers, aromatic bond-order changes, coproducts, and complete atom mapping remain unresolved."],
            "model_votes":[]}


reactions=[]


def reaction(suffix, evidence, record_type="specific_experiment", categories=None,
             location=None, name=None, transformation=None, scheme_ids=(), bonds=(),
             class_coverage="partial", unresolved=None):
    rid=f"{DOC}:{suffix}"
    r=st.new_reaction(rid,DOC,evidence)
    r["record_type"]=record_type
    categories=categories or ["other"]
    r["experiment_categories"]=fe(categories,categories,evidence)
    r["location"]=fe(location or suffix, {"pages":sorted({int(x.split(':page')[1].split(':')[0]) for x in evidence if ':page' in x}),"region":location or suffix}, evidence)
    r["general_scheme_ids"]=list(scheme_ids)
    r["reaction_name"]=fe(name,name,evidence) if name else missing(evidence)
    r["transformation_type"]=inferred(transformation,evidence,"Transformation wording normalized from the source title/drawing.",transformation) if transformation else missing(evidence)
    r["classification"]=classification(evidence,bonds,class_coverage)
    r["review"]={"status":"needs_review","method":"model","reviewer":None,"checked_fields":[],
                 "unresolved":list(unresolved or ["Machine extraction requires expert source review."])}
    r["completeness"]="partial"
    reactions.append(r)
    return r


def yield_outcome(r, oid, target, pct, evidence, basis=None, scope=None):
    o=st.new_outcome(oid,"yield",[target],evidence)
    o["product_role"]=fe("product","product",evidence)
    o["value"]=fe(f"{pct}%",{"number":pct,"unit":"%"},evidence)
    o["yield_type"]=fe("reported yield","reported yield",evidence)
    o["measurement_basis"]=fe(basis,basis,evidence) if basis else missing(evidence)
    o["scope"]=fe(scope,scope,evidence) if scope else missing(evidence)
    r["outcomes"].append(o)


def no_reaction_outcome(r, oid, target, evidence, raw="no reaction"):
    o=st.new_outcome(oid,"no_reaction",[target],evidence)
    o["product_role"]=na(evidence)
    o["value"]=fe(raw,True,evidence)
    o["yield_type"]=na(evidence)
    o["measurement_basis"]=missing(evidence)
    o["scope"]=fe("reported control or 0% outcome","reported control or 0% outcome",evidence)
    r["outcomes"].append(o)


# Page 1: four prior-work examples and this-work generic transformation.
lit_specs=[
 ("literature_mukai_sn2prime","Mukai (2011) intramolecular SN2-prime cyclization","base-mediated indole formation"),
 ("literature_mukai_swamy_michael","Mukai (2004)/Swamy (2014) allene Michael cyclizations","intramolecular conjugate addition"),
 ("literature_kanai","Kanai (2014) copper-catalyzed asymmetric cyclization/carbonyl addition","asymmetric allene aminocupration and carbonyl addition"),
 ("literature_arisawa","Arisawa (2019) palladium-mediated cycloisomerization","palladium-mediated allene cyclization"),
]
for suffix,nm,tr in lit_specs:
    reaction(suffix,ev(1),"literature_example",["other"],"page 1 Scheme 1 prior-work row",nm,tr,
             class_coverage="unresolved",unresolved=["Prior literature example only; yields and full experimental details are not supplied in this article."])

generic=reaction("generic_oxidative_indole",ev(1)+ev(2),"generic_scheme",["substrate_scope"],
    "page 1 this-work row and page 2 Scheme 2 header","Intramolecular indole-forming oxidative amino-hydroxylation",
    "o-allenyl aniline to 2-(acyloxyalkyl)indole",bonds=(("C","N"),("C","O")),
    unresolved=["Generic R1/R2 identities vary; no single experiment is represented by this record."])
generic["participants"]=[part("generic_substrate","m_generic_o_allenyl_aniline","reactant",generic["evidence_ids"]),
                         part("generic_product","m_generic_indole_product","product",generic["evidence_ids"],product_role="generic target")]
generic["procedure"]["original_text"]=fe("Pb(OCOR2)4, CH2Cl2, rt", "Pb(OCOR2)4 in CH2Cl2 at room temperature", generic["evidence_ids"])

# Page 2 Table 1 optimization rows.
opt_rows=[
 (1,"m_pida",1.2,"m_tfe","rt",None,0,0.3),(2,"m_pida",1.2,"m_hfip","rt",None,0,0.3),
 (3,"m_pida",1.2,"m_dcm","rt",None,4,0.3),(4,"m_lta",0.5,"m_dcm","rt",None,21,0.3),
 (5,"m_lta",1.0,"m_dcm","rt",None,40,0.3),(6,"m_lta",1.2,"m_dcm","rt",None,51,0.3),
 (7,"m_lta",1.5,"m_dcm","rt",None,81,0.3),(8,"m_lta",1.2,"m_dcm",{"number":40,"unit":"°C"},{"number":4,"unit":"h"},63,0.3),
 (9,"m_lta",1.2,"m_dce",{"number":60,"unit":"°C"},{"number":4,"unit":"h"},52,0.3),
 (10,"m_lta",1.5,"m_dcm","rt",None,44,0.1),(11,"m_lta",1.5,"m_dcm","rt",None,50,0.2),
 (12,"m_lta",1.5,"m_dcm","rt",None,51,0.4),(13,"m_lta",1.5,"m_benzene","rt",None,33,0.3),
 (14,"m_lta",1.5,"m_dcm","rt",None,55,0.3),
]
for idx,oxidant,eqv,solv,temp,duration,yld,conc in opt_rows:
    r=reaction(f"table1_entry_{idx:02d}",ev(2),categories=["optimization"],location=f"page 2 Table 1 entry {idx}",
      name="Optimization of intramolecular indole-forming oxidative amino-hydroxylation",
      transformation="oxidative amino-hydroxylation and indole formation",scheme_ids=[generic["reaction_id"]],
      bonds=() if yld==0 else (("C","N"),("C","O")),class_coverage="not_applicable" if yld==0 else "partial",
      unresolved=["Yield measurement basis is not printed in Table 1.","Overnight duration is qualitative and not converted to hours."] if duration is None else ["Yield measurement basis is not printed in Table 1."])
    parts=[part("substrate_1","m_1","reactant",ev(2),eq={"number":1,"unit":"equiv"}),
           part("oxidant",oxidant,"reagent",ev(2),eq={"number":eqv,"unit":"equiv"}),
           part("solvent",solv,"solvent",ev(2))]
    if idx==14: parts.append(part("cosolvent_acoh","m_acoh","solvent",ev(2)))
    if yld>0: parts.append(part("product_2","m_2","product",ev(2),product_role="reported product"))
    set_refs(parts,"m_1");r["participants"]=parts
    r["reaction_scale"]=missing(ev(2))
    r["solvent_conditions"]["composition"]=fe("DCM/AcOH" if idx==14 else molecules[solv]["names"]["original"]["normalized_value"],
                                                  ["DCM","AcOH"] if idx==14 else [molecules[solv]["names"]["original"]["normalized_value"]],ev(2))
    r["solvent_conditions"]["mixture_ratio"]=missing(ev(2))
    r["solvent_conditions"]["concentration"]=fe(f"{conc} M",{"number":conc,"unit":"M"},ev(2))
    r["solvent_conditions"]["concentration_basis"]=fe("allene substrate 1","m_1",ev(2))
    dur=duration if duration is not None else "overnight"
    r["steps"]=[step("reaction",ev(2),[p["participant_id"] for p in parts],"stir reaction mixture",temp,dur)]
    if yld==0:no_reaction_outcome(r,"no_product","m_1",ev(2),"0% yield")
    else:yield_outcome(r,"yield_2","m_2",yld,ev(2),scope=f"Table 1 entry {idx}")

# Page 2 Scheme 2 scope.
scope2=[(2,"m_1","m_2",81),(3,"m_s29","m_3",71),(4,"m_s28","m_4",60),
        (5,"m_s30","m_5",61),(6,"m_s31","m_6",72),(7,"m_s32","m_7",64),
        (8,"m_s33","m_8",0),(9,"m_s27","m_9",75),(10,"m_s34","m_10",73),
        (11,"m_s35","m_11",65),(12,"m_s36","m_12",56)]
for num,sub,prod,yld in scope2:
    r=reaction(f"scheme2_product_{num}",ev(2),categories=["substrate_scope"],location=f"page 2 Scheme 2 product {num}",
      name="Scope and limitation for indole formation with unsubstituted allenes",
      transformation="oxidative amino-hydroxylation and indole formation",scheme_ids=[generic["reaction_id"]],
      bonds=() if yld==0 else (("C","N"),("C","O")),class_coverage="not_applicable" if yld==0 else "partial",
      unresolved=["LTA loading, reaction time, concentration, and yield measurement basis are not printed in Scheme 2."])
    parts=[part("substrate",sub,"reactant",ev(2)),part("lta","m_lta","reagent",ev(2)),part("dcm","m_dcm","solvent",ev(2))]
    if yld>0:parts.append(part("product",prod,"product",ev(2),product_role="reported product"))
    r["participants"]=parts
    r["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(2))
    r["steps"]=[step("reaction",ev(2),[p["participant_id"] for p in parts],"treat allene with LTA in DCM","rt")]
    if yld==0:no_reaction_outcome(r,"no_product",sub,ev(2),"0%")
    else:yield_outcome(r,"yield",prod,yld,ev(2),scope=f"Scheme 2 product {num}")

# Page 3 Scheme 3 lead(IV) carboxylates.
scheme3=[(13,"m_pb_obz4","m_13",57),(14,"m_pb_pfbz4","m_14",44),(15,"m_pb_pmbz4","m_15",31),
         (16,"m_pb_ptfbz4","m_16",54),(17,"m_pb_cinn4","m_17",67)]
for num,ox,prod,yld in scheme3:
    r=reaction(f"scheme3_product_{num}",ev(3),categories=["substrate_scope"],location=f"page 3 Scheme 3 product {num}",
      name="Oxidation with lead(IV) benzoate/cinnamate",transformation="oxidative amino-hydroxylation and indole formation",
      scheme_ids=[generic["reaction_id"]],bonds=(("C","N"),("C","O")),
      unresolved=["Oxidant loading, reaction time, concentration, and yield measurement basis are not printed in Scheme 3."])
    parts=[part("substrate","m_s27","reactant",ev(3)),part("oxidant",ox,"reagent",ev(3)),
           part("dcm","m_dcm","solvent",ev(3)),part("product",prod,"product",ev(3),product_role="reported product")]
    r["participants"]=parts;r["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(3))
    r["steps"]=[step("reaction",ev(3),[p["participant_id"] for p in parts],"treat allene with lead(IV) carboxylate","rt")]
    yield_outcome(r,"yield",prod,yld,ev(3),scope=f"Scheme 3 product {num}")

# Page 3 Scheme 4 terminal substitution.
for sub,prod,yld in [("m_18","m_22",48),("m_19","m_23",49),("m_20","m_24",14),("m_21","m_25",14)]:
    num=prod.split("_")[-1]
    r=reaction(f"scheme4_product_{num}",ev(3),categories=["substrate_scope"],location=f"page 3 Scheme 4 product {num}",
      name="Scope of indole formation with terminally substituted allenes",
      transformation="oxidative amino-hydroxylation and indole formation",scheme_ids=[generic["reaction_id"]],
      bonds=(("C","N"),("C","O")),unresolved=["LTA loading, reaction time, concentration, and yield measurement basis are not printed in Scheme 4.",
      "Products 22 and 23 contain an unassigned stereogenic carbon; no ee/dr/configuration is reported."] if num in {"22","23"} else ["LTA loading, reaction time, concentration, and yield measurement basis are not printed in Scheme 4."])
    parts=[part("substrate",sub,"reactant",ev(3)),part("lta","m_lta","reagent",ev(3)),part("dcm","m_dcm","solvent",ev(3)),part("product",prod,"product",ev(3),product_role="reported product")]
    r["participants"]=parts;r["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(3))
    r["steps"]=[step("reaction",ev(3),[p["participant_id"] for p in parts],"treat allene with LTA","rt")]
    yield_outcome(r,"yield",prod,yld,ev(3),scope=f"Scheme 4 product {num}")

# Page 3 Scheme 5 controls.
controls=[("control_26_to_27","m_26",[("m_27",44)],False,(("C","O"),),{"number":5,"unit":"equiv"}),
          ("control_28_no_reaction","m_28",[],True,(),None),
          ("control_29_to_30","m_29",[("m_30",17)],False,(("C","N"),("C","O")),None),
          ("control_31_no_reaction","m_31",[],True,(),None),
          ("control_32_products_33_34","m_32",[("m_33",36),("m_34",9)],False,(("C","O"),("C","O"),("C","N")),None)]
for suffix,sub,products,nr,bonds,ltaeq in controls:
    r=reaction(suffix,ev(3),categories=["control","mechanistic_study"],location=f"page 3 Scheme 5 {suffix}",
      name="Lead(IV) acetate control experiment",transformation="control oxidation",
      bonds=bonds,class_coverage="not_applicable" if nr else "partial",
      unresolved=["Reaction time, concentration, and outcome measurement basis are not printed in Scheme 5."])
    parts=[part("substrate",sub,"reactant",ev(3)),part("lta","m_lta","reagent",ev(3),eq=ltaeq),part("dcm","m_dcm","solvent",ev(3))]
    for j,(prod,yld) in enumerate(products,1):parts.append(part(f"product_{j}",prod,"product",ev(3),product_role="reported product" if j==1 else "minor product"))
    if ltaeq:set_refs(parts,sub)
    r["participants"]=parts;r["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(3))
    r["steps"]=[step("reaction",ev(3),[p["participant_id"] for p in parts],"treat control substrate with LTA","rt")]
    if nr:no_reaction_outcome(r,"no_reaction",sub,ev(3))
    else:
        for j,(prod,yld) in enumerate(products,1):yield_outcome(r,f"yield_{j}",prod,yld,ev(3),scope=f"Scheme 5 {suffix}")

# Pages 3--4: isotope/KIE experiment, plus the proposed mechanism.
kie=reaction("isotope_kie_35_to_36",ev(3)+ev(4),categories=["mechanistic_study"],location="page 4 Scheme 6 isotope experiment and page 3 discussion",
    name="Deuterium kinetic isotope effect experiment",transformation="oxidative amino-hydroxylation of deuterated allene",
    scheme_ids=[generic["reaction_id"]],bonds=(("C","N"),("C","O")),
    unresolved=["The displayed comparison gives kH/kD = 0.87 but does not print individual kinetic traces, exact charges, concentration, duration, or numerical product yield on pages 1--6."])
kie["participants"]=[part("substrate_35","m_35","reactant",kie["evidence_ids"]),part("lta","m_lta","reagent",kie["evidence_ids"]),part("dcm","m_dcm","solvent",kie["evidence_ids"]),part("product_36","m_36","product",kie["evidence_ids"],product_role="deuterated product")]
kie["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],kie["evidence_ids"])
kie["steps"]=[step("reaction",kie["evidence_ids"],[p["participant_id"] for p in kie["participants"]],"perform deuterated reaction for kinetic comparison","rt")]
o=st.new_outcome("kie","other",["m_36"],kie["evidence_ids"]);o["product_role"]=fe("kinetic comparison product","kinetic comparison product",kie["evidence_ids"]);o["value"]=fe("kH/kD = 0.87",{"metric":"kH/kD","number":0.87},kie["evidence_ids"]);o["yield_type"]=na(kie["evidence_ids"]);o["measurement_basis"]=fe("proton NMR kinetic comparison","proton NMR kinetic comparison",kie["evidence_ids"]);o["scope"]=fe("formation of 36 compared with protio product 9","36 versus 9",kie["evidence_ids"]);kie["outcomes"].append(o)

mechanism=reaction("scheme6_mechanistic_proposal",ev(4),"proposed_mechanism",["mechanistic_study"],"page 4 Scheme 6",
    "Mechanistic proposal","Wessely dearomatization, cyclization, carboxylate trapping, and rearomatization",
    scheme_ids=[generic["reaction_id"]],bonds=(("C","N"),("C","O")),
    unresolved=["A--F are proposed/observed intermediates in a mechanistic scheme, not separately isolated synthetic experiments; only C and E are described as observed by proton NMR."])

# General Procedure A.
gpA=reaction("general_procedure_A",ev(4),"generic_scheme",["other"],"page 4 General Procedure A",
    "General Procedure A: methyl carbamate protection","methyl carbamate protection",bonds=(("C","N"),),
    unresolved=["Generic derivative identity varies by experiment; the record is a procedure template, not one isolated product."])
gpA_parts=[part("iodoaniline","m_generic_iodoaniline","reactant",ev(4),eq={"number":1,"unit":"equiv"},amount={"number":10,"unit":"mmol"}),
 part("pyridine","m_pyridine","reagent",ev(4),eq={"number":3.3,"unit":"equiv"},amount={"number":33,"unit":"mmol"},volume={"number":2.66,"unit":"mL"}),
 part("methyl_chloroformate","m_methyl_chloroformate","reagent",ev(4),eq={"number":3,"unit":"equiv"},volume={"number":2.32,"unit":"mL"}),
 part("dcm","m_dcm","solvent",ev(4),volume={"number":20,"unit":"mL"}),part("product","m_generic_carbamate","product",ev(4),product_role="generic protected product"),
 part("water","m_water","workup_agent",ev(4),"workup",volume={"number":40,"unit":"mL"}),part("brine","m_brine","workup_agent",ev(4),"workup"),part("na2so4","m_na2so4","workup_agent",ev(4),"workup"),
 part("etoac","m_etoac","purification_agent",ev(4),"purification"),part("hexanes","m_hexanes","purification_agent",ev(4),"purification")]
set_refs(gpA_parts,"m_generic_iodoaniline");gpA["participants"]=gpA_parts
gpA["reaction_scale"]=fe("10.00 mmol 2-iodoaniline derivative",{"reference_molecule_id":"m_generic_iodoaniline","amount":{"number":10,"unit":"mmol"}},ev(4))
gpA["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(4));gpA["solvent_conditions"]["total_volume"]=fe("20 mL",{"number":20,"unit":"mL"},ev(4));gpA["solvent_conditions"]["concentration"]=fe("1.0 M",{"number":1.0,"unit":"M"},ev(4));gpA["solvent_conditions"]["concentration_basis"]=fe("2-iodoaniline derivative","m_generic_iodoaniline",ev(4))
gpA["steps"]=[step("addition_0c",ev(4),["iodoaniline","pyridine","methyl_chloroformate","dcm"],"add methyl chloroformate over 15 min",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"nitrogen"),step("warm_and_hold",ev(4),["iodoaniline","pyridine","methyl_chloroformate","dcm"],"allow to warm to room temperature and stir","room temperature",{"number":4,"unit":"h"},"nitrogen"),step("water_quench_extract",ev(4),["water","dcm","brine","na2so4"],"quench with water; extract three times with DCM; wash with brine; dry over Na2SO4","room temperature",None,None,"workup"),step("chromatography",ev(4),["etoac","hexanes","product"],"silica gel flash chromatography, 15:85 EtOAc/hexanes",None,None,None,"purification")]
gpA["procedure"]["original_text"]=fe("General Procedure A paragraph","Full operation retained in ordered steps.",ev(4))

# S1--S8, with explicit Procedure A inheritance; S6 is a separate literature procedure.
a_rows={1:("m_a1",10.00,{"number":2.64,"unit":"g"},95),2:("m_a2",10.00,{"number":2.77,"unit":"g"},94),3:("m_a3",10.00,{"number":2.99,"unit":"g"},96),4:("m_a4",10.00,{"number":3.41,"unit":"g"},96),5:("m_a5",10.00,{"number":2.75,"unit":"g"},95),7:("m_a7",10.00,{"number":2.91,"unit":"g"},95),8:("m_a8",10.00,{"number":3.26,"unit":"g"},95)}
for n,(start,scale,mass,yld) in a_rows.items():
    evidence=ev(4)+ev(5) if n==2 else ev(5) if n>=3 else ev(4)
    r=reaction(f"prep_S{n}",evidence,categories=["other"],location=f"pages 4--5 preparation of S{n}" if n==2 else f"page {4 if n==1 else 5} preparation of S{n}",
      name="Methyl carbamate protection",transformation="methyl carbamate formation",scheme_ids=[gpA["reaction_id"]],bonds=(("C","N"),),
      unresolved=["Machine structure and stereochemical status require expert review."] + (["S2 characterization values printed across pages 4--5 appear internally inconsistent with the named molecular formula; they are retained as source context, not used to alter the reaction identity."] if n==2 else []))
    parts=[part("starting_material",start,"reactant",evidence,eq={"number":1,"unit":"equiv"},amount={"number":scale,"unit":"mmol"}),part("pyridine","m_pyridine","reagent",ev(4),eq={"number":3.3,"unit":"equiv"},amount={"number":33,"unit":"mmol"},volume={"number":2.66,"unit":"mL"}),part("methyl_chloroformate","m_methyl_chloroformate","reagent",ev(4),eq={"number":3,"unit":"equiv"},volume={"number":2.32,"unit":"mL"}),part("dcm","m_dcm","solvent",ev(4),volume={"number":20,"unit":"mL"}),part("product",f"m_s{n}","product",evidence,mass=mass,product_role="isolated product"),part("water","m_water","workup_agent",ev(4),"workup",volume={"number":40,"unit":"mL"}),part("brine","m_brine","workup_agent",ev(4),"workup"),part("na2so4","m_na2so4","workup_agent",ev(4),"workup"),part("etoac","m_etoac","purification_agent",ev(4),"purification"),part("hexanes","m_hexanes","purification_agent",ev(4),"purification")]
    set_refs(parts,start);r["participants"]=parts
    r["reaction_scale"]=fe(f"{scale} mmol starting material",{"reference_molecule_id":start,"amount":{"number":scale,"unit":"mmol"}},evidence)
    r["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(4));r["solvent_conditions"]["total_volume"]=fe("20 mL",{"number":20,"unit":"mL"},ev(4));r["solvent_conditions"]["concentration"]=fe("1.0 M",{"number":1.0,"unit":"M"},ev(4));r["solvent_conditions"]["concentration_basis"]=fe("starting 2-iodoaniline derivative",start,ev(4))
    hold_temp={"number":45,"unit":"°C"} if n==8 else "room temperature";hold_time={"number":24,"unit":"h"} if n==8 else {"number":4,"unit":"h"}
    r["steps"]=[step("addition_0c",ev(4),["starting_material","pyridine","methyl_chloroformate","dcm"],"add methyl chloroformate over 15 min",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"nitrogen"),step("hold",evidence,["starting_material","pyridine","methyl_chloroformate","dcm"],"warm and stir" if n!=8 else "heat in an oil bath before quench",hold_temp,hold_time,"nitrogen"),step("workup",ev(4),["water","dcm","brine","na2so4"],"water quench, three DCM extractions, brine wash, Na2SO4 drying",None,None,None,"workup"),step("purification",ev(4),["etoac","hexanes","product"],"silica gel flash chromatography, 15:85 EtOAc/hexanes",None,None,None,"purification")]
    r["procedure"]["general_procedure_reference"]=fe("general procedure A","general procedure A",evidence);r["procedure"]["modifications"]=fe("45 °C for 24 h before quench","45 °C for 24 h before quench",evidence) if n==8 else na(evidence)
    yield_outcome(r,"yield",f"m_s{n}",yld,evidence,basis="product obtained after the stated purification",scope=f"preparation of S{n}")

literature_preps=[
 (6,"m_s1",7.22,{"number":819,"unit":"mg"},34,"known procedure"),
 (9,"m_a1",18.202,{"number":5.76,"unit":"g"},99,"literature procedure"),
 (10,"m_2iodobenzoic_acid",16.13,{"number":2.14,"unit":"g"},62,"literature procedure"),
 (11,"m_2iodobenzoic_acid",16.13,{"number":4.32,"unit":"g"},68,"literature procedure"),
]
for n,start,scale,mass,yld,proc in literature_preps:
    r=reaction(f"prep_S{n}",ev(5),categories=["other"],location=f"page 5 preparation of S{n}",name=f"Preparation of S{n} following a {proc}",transformation="carbamate preparation",bonds=(("C","N"),),unresolved=["The cited procedure, reagents, conditions, workup, purification, and yield measurement basis are not reproduced on pages 1--6."])
    r["participants"]=[part("starting_material",start,"reactant",ev(5),amount={"number":scale,"unit":"mmol"}),part("product",f"m_s{n}","product",ev(5),mass=mass,product_role="reported product")]
    r["reaction_scale"]=fe(f"{scale} mmol starting material",{"reference_molecule_id":start,"amount":{"number":scale,"unit":"mmol"}},ev(5));r["steps"]=[step("reported_preparation",ev(5),["starting_material","product"],f"follow {proc}; details not reproduced")]
    r["procedure"]["general_procedure_reference"]=fe(proc,proc,ev(5));yield_outcome(r,"yield",f"m_s{n}",yld,ev(5),basis="product obtained following the cited procedure",scope=f"preparation of S{n}")

# General Procedure B.
gpB=reaction("general_procedure_B",ev(5),"generic_scheme",["other"],"page 5 General Procedure B","General Procedure B: Sonogashira coupling of protected 2-iodoaniline derivatives","Sonogashira coupling",bonds=(("C","C"),),unresolved=["Generic protected iodoaniline and alkyne identities vary; this is a procedure template."])
gpB_parts=[part("iodoaniline","m_generic_carbamate","reactant",ev(5),eq={"number":1,"unit":"equiv"}),part("alkyne","m_generic_alkyne","reactant",ev(5),eq={"number":1.2,"unit":"equiv"}),part("pd","m_pdcl2pph3","catalyst",ev(5),eq={"number":0.02,"unit":"equiv"}),part("cui","m_cui","catalyst",ev(5),eq={"number":0.03,"unit":"equiv"}),part("thf","m_thf","solvent",ev(5)),part("et3n","m_et3n","solvent",ev(5)),part("product","m_generic_alkynyl_product","product",ev(5),product_role="generic coupled product"),part("nh4cl","m_nh4cl","workup_agent",ev(5),"workup"),part("etoac","m_etoac","workup_agent",ev(5),"workup"),part("brine","m_brine","workup_agent",ev(5),"workup"),part("na2so4","m_na2so4","workup_agent",ev(5),"workup"),part("hexanes","m_hexanes","purification_agent",ev(5),"purification")]
set_refs(gpB_parts,"m_generic_carbamate");gpB["participants"]=gpB_parts
gpB["solvent_conditions"]["composition"]=fe("THF/NEt3",["THF","NEt3"],ev(5));gpB["solvent_conditions"]["mixture_ratio"]=fe("4:1",{"THF":4,"NEt3":1},ev(5));gpB["solvent_conditions"]["concentration"]=fe("0.5 M",{"number":0.5,"unit":"M"},ev(5));gpB["solvent_conditions"]["concentration_basis"]=fe("protected 2-iodoaniline derivative","m_generic_carbamate",ev(5))
gpB["steps"]=[step("nitrogen_and_addition",ev(5),["iodoaniline","alkyne","pd","cui","thf","et3n"],"bubble nitrogen through solution, then add alkyne dropwise","room temperature",None,"nitrogen"),step("reaction_hold",ev(5),["iodoaniline","alkyne","pd","cui","thf","et3n"],"stir reaction mixture","room temperature",{"minimum":3,"maximum":6,"unit":"h"},"nitrogen"),step("workup",ev(5),["nh4cl","etoac","brine","na2so4"],"quench with saturated aqueous NH4Cl; extract three times with EtOAc; brine wash; Na2SO4 drying",None,None,None,"workup"),step("purification",ev(5),["etoac","hexanes","product"],"silica gel chromatography, 20:80 to 40:60 EtOAc/hexanes",None,None,None,"purification")]
gpB["procedure"]["original_text"]=fe("General Procedure B paragraph","Full operation retained in ordered steps.",ev(5))

b_rows={12:("m_s1","m_propargyl_alcohol",10.83,{"number":1.78,"unit":"g"},80),13:("m_s2","m_propargyl_alcohol",6.79,{"number":1.30,"unit":"g"},86),14:("m_s3","m_propargyl_alcohol",6.43,{"number":1.19,"unit":"g"},77),15:("m_s4","m_propargyl_alcohol",4.12,{"number":1.04,"unit":"g"},89),16:("m_s5","m_propargyl_alcohol",4.12,{"number":613,"unit":"mg"},67),17:("m_s6","m_propargyl_alcohol",2.45,{"number":552,"unit":"mg"},88),18:("m_s7","m_propargyl_alcohol",4.40,{"number":855,"unit":"g"},83),19:("m_s8","m_propargyl_alcohol",4.12,{"number":1.04,"unit":"g"},95),20:("m_s9","m_propargyl_alcohol",3.23,{"number":765,"unit":"mg"},97),21:("m_s10","m_propargyl_alcohol",4.25,{"number":951,"unit":"mg"},80),22:("m_s11","m_propargyl_alcohol",2.55,{"number":707,"unit":"mg"},86),23:("m_s1","m_but3yn2ol",5.56,{"number":1.21,"unit":"g"},99),24:("m_s1","m_2methylbut3yn2ol",10.83,{"number":2.49,"unit":"g"},98),25:("m_s4","m_but3yn2ol",4.21,{"number":1.23,"unit":"g"},98),26:("m_s4","m_2methylbut3yn2ol",1.40,{"number":275,"unit":"mg"},62)}
for n,(start,alkyne,scale,mass,yld) in b_rows.items():
    evidence=ev(5)+ev(6) if n==14 else ev(5) if n<=13 else ev(6)
    unresolved=["Individual reagent masses/volumes are not printed; catalyst/alkyne equivalents and 3--6 h duration are inherited explicitly from General Procedure B."]
    if n in (23,25):unresolved.append("The propargylic alcohol product contains an unassigned stereogenic center; no configuration or enantiomeric composition is reported.")
    if n==18:unresolved.append("The source prints product mass as 855 g at 4.40 mmol and 83%; the apparent internal inconsistency is preserved without silently changing it to mg.")
    r=reaction(f"prep_S{n}",evidence,categories=["other"],location=f"pages 5--6 preparation of S{n}" if n==14 else f"page {5 if n<=13 else 6} preparation of S{n}",name="Sonogashira coupling using General Procedure B",transformation="Sonogashira C-C coupling",scheme_ids=[gpB["reaction_id"]],bonds=(("C","C"),),unresolved=unresolved)
    parts=[part("iodoaniline",start,"reactant",evidence,eq={"number":1,"unit":"equiv"},amount={"number":scale,"unit":"mmol"}),part("alkyne",alkyne,"reactant",evidence,eq={"number":1.2,"unit":"equiv"}),part("pd","m_pdcl2pph3","catalyst",ev(5),eq={"number":0.02,"unit":"equiv"}),part("cui","m_cui","catalyst",ev(5),eq={"number":0.03,"unit":"equiv"}),part("thf","m_thf","solvent",ev(5)),part("et3n","m_et3n","solvent",ev(5)),part("product",f"m_s{n}","product",evidence,mass=mass,product_role="isolated product"),part("nh4cl","m_nh4cl","workup_agent",ev(5),"workup"),part("etoac","m_etoac","workup_agent",ev(5),"workup"),part("brine","m_brine","workup_agent",ev(5),"workup"),part("na2so4","m_na2so4","workup_agent",ev(5),"workup"),part("hexanes","m_hexanes","purification_agent",ev(5),"purification")]
    set_refs(parts,start);r["participants"]=parts;r["reaction_scale"]=fe(f"{scale} mmol protected 2-iodoaniline",{"reference_molecule_id":start,"amount":{"number":scale,"unit":"mmol"}},evidence)
    r["solvent_conditions"]["composition"]=fe("THF/NEt3",["THF","NEt3"],ev(5));r["solvent_conditions"]["mixture_ratio"]=fe("4:1",{"THF":4,"NEt3":1},ev(5));r["solvent_conditions"]["concentration"]=fe("0.5 M",{"number":0.5,"unit":"M"},ev(5));r["solvent_conditions"]["concentration_basis"]=fe("protected 2-iodoaniline derivative",start,ev(5))
    r["steps"]=[step("nitrogen_and_addition",ev(5),["iodoaniline","alkyne","pd","cui","thf","et3n"],"bubble nitrogen, then add alkyne dropwise","room temperature",None,"nitrogen"),step("reaction_hold",ev(5),["iodoaniline","alkyne","pd","cui","thf","et3n"],"stir reaction mixture","room temperature",{"minimum":3,"maximum":6,"unit":"h"},"nitrogen"),step("workup",ev(5),["nh4cl","etoac","brine","na2so4"],"saturated NH4Cl quench; three EtOAc extractions; brine wash; Na2SO4 drying",None,None,None,"workup"),step("purification",ev(5),["etoac","hexanes","product"],"silica gel chromatography, 20:80 to 40:60 EtOAc/hexanes",None,None,None,"purification")]
    r["procedure"]["general_procedure_reference"]=fe("general procedure B","general procedure B",evidence)
    yield_outcome(r,"yield",f"m_s{n}",yld,evidence,basis="product obtained after the stated chromatography",scope=f"preparation of S{n}")

# Pages 7--11: actual image review and experimental continuation.  Records that
# were drawn in Schemes 2--5 are enriched below rather than duplicated.
def extend_evidence(obj, evidence):
    """Add evidence IDs without changing any source value or review status."""
    if isinstance(obj, dict):
        if isinstance(obj.get("evidence_ids"), list):
            obj["evidence_ids"] = list(dict.fromkeys(obj["evidence_ids"] + list(evidence)))
        for value in obj.values():
            extend_evidence(value, evidence)
    elif isinstance(obj, list):
        for value in obj:
            extend_evidence(value, evidence)


def generic_template_molecule(mid, label, name, fixed_structure, variables, evidence):
    unresolved_molecule(mid, label, name, evidence)
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": fe(fixed_structure, fixed_structure, evidence),
        "variables": [{
            "label": variable,
            "attachment_points": fe(f"source-drawn attachment for {variable}", f"source-drawn attachment for {variable}", evidence),
            "definition": fe(f"variable substituent {variable}", f"variable substituent {variable}", evidence),
            "allowed_values": missing(evidence),
            "constraints": fe("varies among the explicitly instantiated compounds", "source-defined substituent", evidence),
        } for variable in variables],
        "instance_mapping": [],
        "reason": fe("General procedure covers a family of substituted compounds.", "generic source template", evidence),
    }
    return mid


def allene_classification(evidence):
    changes = [
        derived({"change_id":"b1","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"H"}],
                 "before_order":0,"after_order":1,"change_type":"formed"}, evidence,
                "Net hydrogen (or deuterium for compound 35) gain is implicit in conversion of the propargylic alcohol to the allene."),
        derived({"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"C"}],
                 "before_order":1,"after_order":2,"change_type":"order_increased"}, evidence,
                "The source structures show promotion of the propargylic C-C single bond to an allene double bond."),
    ]
    c_h = st.bond_label("C", "H", 1)
    db_cc = st.bond_label("C", "C", 2)
    return {"standard_id":st.STANDARD,"mapped_reaction_smiles":missing(evidence),
            "bond_changes":changes,
            "labels":derived([dict(c_h,supporting_change_ids=["b1"]),dict(db_cc,supporting_change_ids=["b2"])], evidence,
                             "Labels follow the explicit net hydrogen gain and C-C bond-order increase."),
            "coverage":inferred("partial",evidence,"The broken C-O bond and C-C triple-bond order decrease are retained as limitations, not new-bond labels."),
            "limitations":["Partial net classification: leaving groups, broken C-O bond, triple-to-double order decrease, and complete atom mapping remain unresolved."],
            "model_votes":[]}


# Reagents and generic structures introduced by General Procedures C--E and the
# page-10 control-substrate preparations.
for mid, label, name, evidence in [
    ("m_pph3","PPh3","triphenylphosphine",ev(7)),
    ("m_ipnbsh","IPNBSH","N-isopropylidene-N'-2-nitrobenzenesulfonyl hydrazine",ev(7)),
    ("m_diad","DIAD","diisopropyl azodicarboxylate",ev(7)),
    ("m_mgso4","MgSO4","magnesium sulfate",ev(7)),
    ("m_et2o","Et2O","diethyl ether",ev(7)),
    ("m_nitromethane","nitromethane","nitromethane",ev(7)),
    ("m_nbsh","NBSH","2-nitrobenzenesulfonylhydrazide",ev(7)),
    ("m_agotf","AgOTf","silver triflate",ev(7)),
    ("m_tfoh","TfOH","triflic acid",ev(7)),
    ("m_nahco3","NaHCO3","sodium bicarbonate",ev(7)),
    ("m_alumina_basic","basic alumina","basic alumina, activity 5",ev(8)),
    ("m_tfe_d1","TFE-d1","deuterated 2,2,2-trifluoroethanol",ev(8)),
    ("m_d2o","D2O","deuterium oxide",ev(8)),
    ("m_o_toluidine","o-toluidine","o-toluidine",ev(10)),
    ("m_pdp4","Pd(PPh3)4","tetrakis(triphenylphosphine)palladium(0)",ev(10)),
    ("m_dmf","DMF","N,N-dimethylformamide",ev(10)),
    ("m_tributyl_vinyl_tin","tributyl(vinyl)tin","tributyl(vinyl)tin",ev(10)),
    ("m_nh4oh_5pct","5% NH4OH","5% ammonium hydroxide",ev(10)),
    ("m_tms_acetylene","TMS/acetylene","trimethylsilylacetylene",ev(10)),
    ("m_k2co3","K2CO3","potassium carbonate",ev(10)),
    ("m_meoh","MeOH","methanol",ev(10)),
]:
    unresolved_molecule(mid,label,name,evidence)

generic_template_molecule("m_generic_propargylic_alcohol","generic propargylic alcohol derivative",
    "carbamate-protected o-alkynyl propargylic alcohol derivative","Ar-C#C-C(OH)(R1)(R2)",["Ar","R1","R2"],ev(7))
generic_template_molecule("m_generic_allene_c","generic o-allenyl aniline derivative",
    "carbamate-protected o-allenyl aniline derivative","Ar-C(R1)=C=C(R2)R3",["Ar","R1","R2","R3"],ev(7))
generic_template_molecule("m_generic_lead_tetracarboxylate","Pb(OCOR)4",
    "lead(IV) tetracarboxylate","Pb(OCOR)4",["R"],ev(8))
molecule("m_silyl_alkyne_intermediate","TMS-protected 29",
    "methyl [2-(2-trimethylsilylethynyl)phenyl]carbamate",
    "COC(=O)Nc1ccccc1C#C[Si](C)(C)C",ev(10))


# Every named molecule in these pages receives page-local structural and stereo
# evidence; an RDKit parse remains diagnostic rather than source verification.
page7_molecules = ["m_s12","m_s13","m_s14","m_s15","m_s16","m_s17","m_s18","m_s19","m_s20","m_s21",
                   "m_s27","m_s28","m_s29","m_1","m_s30","m_s31","m_s32","m_s33","m_s34","m_s35"]
for mid in page7_molecules:
    extend_evidence(molecules[mid],ev(7))
for mid in ["m_s22","m_s23","m_s24","m_s25","m_s26","m_s36","m_18","m_19","m_20","m_21","m_35"]:
    extend_evidence(molecules[mid],ev(8))
for mid in ["m_9","m_4","m_3","m_2","m_5"]:
    extend_evidence(molecules[mid],ev(8))
for mid in ["m_5","m_6","m_7","m_10","m_11","m_12","m_23","m_22","m_24","m_36","m_13","m_14","m_16"]:
    extend_evidence(molecules[mid],ev(9))
for mid in ["m_16","m_15","m_17","m_26","m_32","m_29","m_31","m_27","m_33","m_30"]:
    extend_evidence(molecules[mid],ev(10))


# General Procedure C.
gpC = reaction("general_procedure_C",ev(7),"generic_scheme",["other"],"page 7 General Procedure C",
    "General Procedure C: synthesis of carbamate-protected o-allenyl aniline derivatives",
    "reductive deoxyallenylation of propargylic alcohols",unresolved=[
        "Generic substituents vary among S27--S36 and compounds 1, 18, 19, and 35; this is a procedure template.",
        "Dropwise DIAD addition duration and warm-up ramp time are not reported."])
gpC["classification"] = allene_classification(ev(7))
gpC_parts = [
    part("propargylic_alcohol","m_generic_propargylic_alcohol","reactant",ev(7),eq={"number":1,"unit":"equiv"}),
    part("pph3","m_pph3","reagent",ev(7),eq={"number":1.5,"unit":"equiv"}),
    part("ipnbsh","m_ipnbsh","reagent",ev(7),eq={"number":1.5,"unit":"equiv"}),
    part("diad","m_diad","reagent",ev(7),eq={"number":1.5,"unit":"equiv"}),
    part("thf","m_thf","solvent",ev(7)),part("tfe","m_tfe","solvent",ev(7)),part("water","m_water","solvent",ev(7)),
    part("product","m_generic_allene_c","product",ev(7),product_role="generic allene product"),
    part("et2o","m_et2o","workup_agent",ev(7),"workup"),part("water_wash","m_water","workup_agent",ev(7),"workup"),
    part("brine","m_brine","workup_agent",ev(7),"workup"),part("mgso4","m_mgso4","workup_agent",ev(7),"workup"),
    part("etoac","m_etoac","purification_agent",ev(7),"purification"),part("hexanes","m_hexanes","purification_agent",ev(7),"purification")]
set_refs(gpC_parts,"m_generic_propargylic_alcohol");gpC["participants"]=gpC_parts
gpC["solvent_conditions"]["composition"]=fe("THF; then 1:1 TFE/H2O",["THF","TFE","water"],ev(7))
gpC["solvent_conditions"]["mixture_ratio"]=fe("TFE/H2O 1:1",{"TFE":1,"H2O":1},ev(7))
gpC["solvent_conditions"]["concentration"]=fe("0.2 M",{"number":0.2,"unit":"M"},ev(7))
gpC["solvent_conditions"]["concentration_basis"]=fe("propargyl alcohol derivative","m_generic_propargylic_alcohol",ev(7))
gpC["steps"]=[
    step("charge_and_dissolve",ev(7),["propargylic_alcohol","pph3","ipnbsh","thf"],"charge propargylic alcohol, PPh3, and IPNBSH; dissolve in THF",None,None,"nitrogen"),
    step("cool_and_add_diad",ev(7),["propargylic_alcohol","pph3","ipnbsh","diad","thf"],"cool to 0 °C and add DIAD dropwise",{"number":0,"unit":"°C"},None,"nitrogen"),
    step("hold_at_zero",ev(7),["propargylic_alcohol","pph3","ipnbsh","diad","thf"],"hold after DIAD addition",{"number":0,"unit":"°C"},{"number":5,"unit":"min"},"nitrogen"),
    step("warm_and_hold",ev(7),["propargylic_alcohol","pph3","ipnbsh","diad","thf"],"allow to warm to room temperature and stir","room temperature",{"number":2,"unit":"h"},"nitrogen"),
    step("diazene_formation",ev(7),["tfe","water"],"add 1:1 TFE/H2O and stir","room temperature",{"number":2,"unit":"h"},"nitrogen"),
    step("workup",ev(7),["et2o","water_wash","brine","mgso4"],"partition between Et2O and water; wash organic layer four times with water, then brine; dry over MgSO4",None,None,None,"workup"),
    step("purification",ev(7),["etoac","hexanes","product"],"silica gel chromatography, 15:85 EtOAc/hexanes",None,None,None,"purification")]
gpC["procedure"]["original_text"]=fe("General Procedure C paragraph","Ordered operations retained in the reaction steps.",ev(7))


# General Procedure D.
gpD = reaction("general_procedure_D",ev(7),"generic_scheme",["other"],"page 7 General Procedure D",
    "General Procedure D: synthesis of dimethyl-substituted carbamate-protected o-allenyl aniline derivatives",
    "acid/silver-promoted deoxyallenylation",unresolved=["Generic substrate identity varies; this is a procedure template."])
gpD["classification"] = allene_classification(ev(7))
gpD_parts=[part("propargylic_alcohol","m_generic_propargylic_alcohol","reactant",ev(7),eq={"number":1,"unit":"equiv"}),
    part("nbsh","m_nbsh","reagent",ev(7),eq={"number":1.2,"unit":"equiv"}),part("agotf","m_agotf","reagent",ev(7),eq={"number":0.2,"unit":"equiv"}),
    part("tfoh","m_tfoh","reagent",ev(7),eq={"number":0.1,"unit":"equiv"}),part("nitromethane","m_nitromethane","solvent",ev(7)),
    part("product","m_generic_allene_c","product",ev(7),product_role="generic dimethyl-substituted allene product"),
    part("nahco3","m_nahco3","workup_agent",ev(7),"workup"),part("et2o","m_et2o","workup_agent",ev(7),"workup"),
    part("brine","m_brine","workup_agent",ev(7),"workup"),part("mgso4","m_mgso4","workup_agent",ev(7),"workup"),
    part("etoac","m_etoac","purification_agent",ev(7),"purification"),part("hexanes","m_hexanes","purification_agent",ev(7),"purification")]
set_refs(gpD_parts,"m_generic_propargylic_alcohol");gpD["participants"]=gpD_parts
gpD["solvent_conditions"]["composition"]=fe("nitromethane",["nitromethane"],ev(7));gpD["solvent_conditions"]["concentration"]=fe("0.167 M",{"number":0.167,"unit":"M"},ev(7));gpD["solvent_conditions"]["concentration_basis"]=fe("propargylic alcohol","m_generic_propargylic_alcohol",ev(7))
gpD["steps"]=[step("charge_substrate",ev(7),["propargylic_alcohol","nitromethane"],"dissolve propargylic alcohol in nitromethane"),
    step("add_nbsh",ev(7),["nbsh"],"add NBSH and flush with argon; spin until dissolved",None,None,"argon"),
    step("add_agotf",ev(7),["agotf"],"add AgOTf and flush with argon",None,None,"argon"),
    step("add_tfoh",ev(7),["tfoh"],"add TfOH and replace septum with a Teflon cap",None,None,"argon"),
    step("heated_hold",ev(7),["propargylic_alcohol","nbsh","agotf","tfoh","nitromethane"],"heat in an oil bath and spin",{"number":35,"unit":"°C"},{"number":1,"unit":"h"},"argon"),
    step("workup",ev(7),["nahco3","et2o","brine","mgso4"],"quench with NaHCO3; extract three times with Et2O; wash with brine; dry over MgSO4",None,None,None,"workup"),
    step("purification",ev(7),["etoac","hexanes","product"],"silica gel chromatography, 10:90 EtOAc/hexanes",None,None,None,"purification")]
gpD["procedure"]["original_text"]=fe("General Procedure D paragraph","Ordered operations retained in the reaction steps.",ev(7))


# General Procedure E.
gpE = reaction("general_procedure_E",ev(8),"generic_scheme",["other"],"page 8 General Procedure E",
    "General Procedure E: indole synthesis","oxidative amino-hydroxylation and indole formation",
    bonds=(("C","N"),("C","O")),unresolved=["The lead tetracarboxylate and allene identities vary; this is a procedure template."])
gpE_parts=[part("lead_tetracarboxylate","m_generic_lead_tetracarboxylate","reagent",ev(8),eq={"number":1.5,"unit":"equiv"}),
    part("allene","m_generic_o_allenyl_aniline","reactant",ev(8),eq={"number":1,"unit":"equiv"}),part("dcm","m_dcm","solvent",ev(8)),
    part("product","m_generic_indole_product","product",ev(8),product_role="generic indole product"),
    part("alumina","m_alumina_basic","purification_agent",ev(8),"purification"),part("dcm_eluent","m_dcm","purification_agent",ev(8),"purification"),
    part("etoac","m_etoac","purification_agent",ev(8),"purification"),part("hexanes","m_hexanes","purification_agent",ev(8),"purification")]
set_refs(gpE_parts,"m_generic_o_allenyl_aniline");gpE["participants"]=gpE_parts
gpE["solvent_conditions"]["composition"]=fe("dichloromethane",["dichloromethane"],ev(8));gpE["solvent_conditions"]["concentration"]=fe("0.30 M",{"number":0.30,"unit":"M"},ev(8));gpE["solvent_conditions"]["concentration_basis"]=fe("2-allenylaniline derivative","m_generic_o_allenyl_aniline",ev(8))
gpE["steps"]=[step("prepare_oxidant",ev(8),["lead_tetracarboxylate"],"charge lead tetracarboxylate into a flame-dried 4 mL vial; purge and maintain under nitrogen",None,None,"nitrogen"),
    step("prepare_allene_solution",ev(8),["allene","dcm"],"dissolve the allene in DCM in a separate flame-dried 4 mL vial and load into a syringe",None,None,"nitrogen"),
    step("addition",ev(8),["lead_tetracarboxylate","allene","dcm"],"add allene solution to the stirred lead tetracarboxylate over 30 s","room temperature",{"number":30,"unit":"s"},"nitrogen"),
    step("reaction_hold",ev(8),["lead_tetracarboxylate","allene","dcm"],"continue stirring overnight","room temperature",{"number":15,"unit":"h","qualifier":"approx","raw_text":"approximately 15 h"},"nitrogen"),
    step("alumina_filtration",ev(8),["alumina","dcm_eluent"],"filter through a 1.5 cm diameter basic-alumina plug packed to 5 cm (activity 5), eluting with DCM",None,None,None,"purification"),
    step("silica_purification",ev(8),["etoac","hexanes","product"],"concentrate and immediately purify by silica gel chromatography, 20:80 EtOAc/hexanes",None,None,None,"purification")]
gpE["procedure"]["original_text"]=fe("General Procedure E paragraph","Ordered operations retained in the reaction steps.",ev(8))


def add_product_mass_basis(outcome, mass, evidence, procedure):
    outcome["measurement_basis"] = fe(
        f"product obtained as {mass['number']} {mass['unit']} after {procedure}",
        "product mass after the stated purification",evidence)


def build_allene_prep(suffix, start, product, scale, mass, yld, evidence, proc="C", deuterated=False):
    gp = gpC if proc == "C" else gpD
    r = reaction(suffix,evidence,categories=["other"],location=f"page {min(int(x.split(':page')[1].split(':')[0]) for x in evidence if ':page' in x)} preparation of {molecules[product]['occurrences'][0]['compound_label']['normalized_value']}",
        name=f"Allene synthesis using General Procedure {proc}",transformation="reductive deoxyallenylation of a propargylic alcohol",
        scheme_ids=[gp["reaction_id"]],unresolved=["The source reports product mass and percent yield; no independent yield-assay method is stated."])
    r["classification"] = allene_classification(evidence)
    parts=[part("substrate",start,"reactant",evidence,eq={"number":1,"unit":"equiv"},amount={"number":scale,"unit":"mmol"})]
    if proc == "C":
        parts += [part("pph3","m_pph3","reagent",ev(7),eq={"number":1.5,"unit":"equiv"}),part("ipnbsh","m_ipnbsh","reagent",ev(7),eq={"number":1.5,"unit":"equiv"}),part("diad","m_diad","reagent",ev(7),eq={"number":1.5,"unit":"equiv"}),part("thf","m_thf","solvent",ev(7))]
        if deuterated:
            parts += [part("tfe_d1","m_tfe_d1","solvent",evidence),part("d2o","m_d2o","solvent",evidence)]
        else:
            parts += [part("tfe","m_tfe","solvent",ev(7)),part("water","m_water","solvent",ev(7))]
        parts += [part("et2o","m_et2o","workup_agent",ev(7),"workup"),part("water_wash","m_water","workup_agent",ev(7),"workup"),part("brine","m_brine","workup_agent",ev(7),"workup"),part("mgso4","m_mgso4","workup_agent",ev(7),"workup"),part("etoac","m_etoac","purification_agent",ev(7),"purification"),part("hexanes","m_hexanes","purification_agent",ev(7),"purification")]
    else:
        parts += [part("nbsh","m_nbsh","reagent",ev(7),eq={"number":1.2,"unit":"equiv"}),part("agotf","m_agotf","reagent",ev(7),eq={"number":0.2,"unit":"equiv"}),part("tfoh","m_tfoh","reagent",ev(7),eq={"number":0.1,"unit":"equiv"}),part("nitromethane","m_nitromethane","solvent",ev(7)),part("nahco3","m_nahco3","workup_agent",ev(7),"workup"),part("et2o","m_et2o","workup_agent",ev(7),"workup"),part("brine","m_brine","workup_agent",ev(7),"workup"),part("mgso4","m_mgso4","workup_agent",ev(7),"workup"),part("etoac","m_etoac","purification_agent",ev(7),"purification"),part("hexanes","m_hexanes","purification_agent",ev(7),"purification")]
    parts.append(part("product",product,"product",evidence,mass=mass,product_role="reported product"));set_refs(parts,start);r["participants"]=parts
    r["reaction_scale"]=fe(f"{scale} mmol substrate",{"reference_molecule_id":start,"amount":{"number":scale,"unit":"mmol"}},evidence)
    if proc == "C":
        composition=["THF","TFE-d1","D2O"] if deuterated else ["THF","TFE","water"]
        r["solvent_conditions"]["composition"]=fe("THF; then 1:1 TFE-d1/D2O" if deuterated else "THF; then 1:1 TFE/H2O",composition,evidence+ev(7))
        r["solvent_conditions"]["mixture_ratio"]=fe("1:1",{"TFE-d1" if deuterated else "TFE":1,"D2O" if deuterated else "H2O":1},evidence+ev(7))
        r["solvent_conditions"]["concentration"]=fe("0.2 M",{"number":0.2,"unit":"M"},ev(7));r["solvent_conditions"]["concentration_basis"]=fe("propargylic alcohol derivative",start,ev(7))
        tfe_ids=["tfe_d1","d2o"] if deuterated else ["tfe","water"]
        r["steps"]=[step("cool_and_add_diad",ev(7),["substrate","pph3","ipnbsh","diad","thf"],"cool THF solution to 0 °C and add DIAD dropwise",{"number":0,"unit":"°C"},None,"nitrogen"),step("hold_at_zero",ev(7),["substrate","pph3","ipnbsh","diad","thf"],"hold after DIAD addition",{"number":0,"unit":"°C"},{"number":5,"unit":"min"},"nitrogen"),step("warm_and_hold",ev(7),["substrate","pph3","ipnbsh","diad","thf"],"allow to warm to room temperature and stir","room temperature",{"number":2,"unit":"h"},"nitrogen"),step("diazene_formation",evidence+ev(7),tfe_ids,"add 1:1 deuterated solvent mixture and stir" if deuterated else "add 1:1 TFE/H2O and stir","room temperature",{"number":2,"unit":"h"},"nitrogen"),step("workup",ev(7),["et2o","water_wash","brine","mgso4"],"Et2O/water partition; four water washes; brine wash; MgSO4 drying",None,None,None,"workup"),step("purification",ev(7),["etoac","hexanes","product"],"silica gel chromatography, 15:85 EtOAc/hexanes",None,None,None,"purification")]
    else:
        r["solvent_conditions"]["composition"]=fe("nitromethane",["nitromethane"],ev(7));r["solvent_conditions"]["concentration"]=fe("0.167 M",{"number":0.167,"unit":"M"},ev(7));r["solvent_conditions"]["concentration_basis"]=fe("propargylic alcohol",start,ev(7))
        r["steps"]=[step("add_nbsh",ev(7),["substrate","nbsh","nitromethane"],"add NBSH and flush with argon; spin until dissolved",None,None,"argon"),step("add_agotf",ev(7),["agotf"],"add AgOTf and flush with argon",None,None,"argon"),step("add_tfoh",ev(7),["tfoh"],"add TfOH and replace septum with a Teflon cap",None,None,"argon"),step("heated_hold",ev(7),["substrate","nbsh","agotf","tfoh","nitromethane"],"heat in oil bath and spin",{"number":35,"unit":"°C"},{"number":1,"unit":"h"},"argon"),step("workup",ev(7),["nahco3","et2o","brine","mgso4"],"NaHCO3 quench; three Et2O extractions; brine wash; MgSO4 drying",None,None,None,"workup"),step("purification",ev(7),["etoac","hexanes","product"],"silica gel chromatography, 10:90 EtOAc/hexanes",None,None,None,"purification")]
    r["procedure"]["general_procedure_reference"]=fe(f"general procedure {proc}",f"general procedure {proc}",evidence)
    if deuterated:r["procedure"]["modifications"]=fe("TFE and H2O were replaced by TFE-d1 and D2O","TFE-d1/D2O substitution",evidence)
    yield_outcome(r,"yield",product,yld,evidence,scope=f"preparation of {molecules[product]['occurrences'][0]['compound_label']['normalized_value']}")
    add_product_mass_basis(r["outcomes"][-1],mass,evidence,f"General Procedure {proc}")
    return r


allene_rows=[
    ("prep_S27","m_s12","m_s27",2.03,{"number":305,"unit":"mg"},79,ev(7),"C",False),
    ("prep_S28","m_s13","m_s28",1.35,{"number":181,"unit":"mg"},65,ev(7),"C",False),
    ("prep_S29","m_s14","m_s29",1.67,{"number":188,"unit":"mg"},50,ev(7),"C",False),
    ("prep_1","m_s15","m_1",0.42,{"number":24,"unit":"mg"},43,ev(7),"C",False),
    ("prep_S30","m_s16","m_s30",1.83,{"number":263,"unit":"mg"},71,ev(7),"C",False),
    ("prep_S31","m_s17","m_s31",1.38,{"number":257,"unit":"mg"},76,ev(7),"C",False),
    ("prep_S32","m_s18","m_s32",1.28,{"number":149,"unit":"mg"},73,ev(7),"C",False),
    ("prep_S33","m_s19","m_s33",0.85,{"number":190,"unit":"mg"},87,ev(7),"C",False),
    ("prep_S34","m_s20","m_s34",0.89,{"number":101,"unit":"mg"},49,ev(7),"C",False),
    ("prep_S35","m_s21","m_s35",1.71,{"number":306,"unit":"mg"},68,ev(7),"C",False),
    ("prep_S36","m_s22","m_s36",0.75,{"number":138,"unit":"mg"},61,ev(7)+ev(8),"C",False),
    ("prep_19","m_s23","m_19",0.51,{"number":42,"unit":"mg"},41,ev(8),"C",False),
    ("prep_21","m_s24","m_21",0.64,{"number":58,"unit":"mg"},41,ev(8),"D",False),
    ("prep_18","m_s25","m_18",1.21,{"number":134,"unit":"mg"},47,ev(8),"C",False),
    ("prep_20","m_s26","m_20",0.32,{"number":38,"unit":"mg"},40,ev(8),"D",False),
    ("prep_35","m_s12","m_35",1.35,{"number":200,"unit":"mg"},78,ev(8),"C",True),
]
for row in allene_rows:
    rr=build_allene_prep(*row)
    if row[2] in {"m_18","m_19"}:
        rr["review"]["unresolved"].append("The allene is axially stereogenic, but no configuration or enantiomeric composition is reported; axial stereochemistry remains unassigned and not encoded.")


def find_reaction(suffix):
    rid=f"{DOC}:{suffix}"
    return next(r for r in reactions if r["reaction_id"]==rid)


def enrich_gpe(suffix, sub, prod, scale_raw, scale_value, mass, yld, detail_evidence, oxidant="m_lta", conflict=None):
    r=find_reaction(suffix)
    r["evidence_ids"]=list(dict.fromkeys(r["evidence_ids"]+ev(8)+detail_evidence))
    pages=sorted({int(x.split(':page')[1].split(':')[0]) for x in r["evidence_ids"] if ':page' in x})
    r["location"]=fe(f"scheme entry and detailed experimental paragraph on pages {pages}",{"pages":pages,"region":"scheme entry merged with General Procedure E experimental paragraph"},r["evidence_ids"])
    if gpE["reaction_id"] not in r["general_scheme_ids"]:r["general_scheme_ids"].append(gpE["reaction_id"])
    parts=[part("substrate",sub,"reactant",detail_evidence,eq={"number":1,"unit":"equiv"},amount=scale_value),part("oxidant",oxidant,"reagent",ev(8),eq={"number":1.5,"unit":"equiv"}),part("dcm","m_dcm","solvent",ev(8)),part("product",prod,"product",detail_evidence,mass=mass,product_role="reported product"),part("alumina","m_alumina_basic","purification_agent",ev(8),"purification"),part("dcm_eluent","m_dcm","purification_agent",ev(8),"purification"),part("etoac","m_etoac","purification_agent",ev(8),"purification"),part("hexanes","m_hexanes","purification_agent",ev(8),"purification")]
    set_refs(parts,sub);r["participants"]=parts
    r["reaction_scale"]=fe(scale_raw,{"reference_molecule_id":sub,"amount":scale_value},detail_evidence)
    r["solvent_conditions"]["composition"]=fe("dichloromethane",["dichloromethane"],ev(8));r["solvent_conditions"]["concentration"]=fe("0.30 M",{"number":0.30,"unit":"M"},ev(8));r["solvent_conditions"]["concentration_basis"]=fe("2-allenylaniline derivative",sub,ev(8))
    r["steps"]=[step("prepare_oxidant",ev(8),["oxidant"],"charge oxidant into a flame-dried 4 mL vial; purge and maintain under nitrogen",None,None,"nitrogen"),step("prepare_allene_solution",ev(8),["substrate","dcm"],"dissolve substrate in DCM in a separate flame-dried 4 mL vial and load into a syringe",None,None,"nitrogen"),step("addition",ev(8),["substrate","oxidant","dcm"],"add substrate solution to stirred oxidant over 30 s","room temperature",{"number":30,"unit":"s"},"nitrogen"),step("reaction_hold",ev(8),["substrate","oxidant","dcm"],"continue stirring overnight","room temperature",{"number":15,"unit":"h","qualifier":"approx","raw_text":"approximately 15 h"},"nitrogen"),step("alumina_filtration",ev(8),["alumina","dcm_eluent"],"filter through a 1.5 cm diameter, 5 cm packed basic-alumina plug (activity 5), eluting with DCM",None,None,None,"purification"),step("silica_purification",ev(8),["etoac","hexanes","product"],"concentrate and immediately purify by silica gel chromatography, 20:80 EtOAc/hexanes",None,None,None,"purification")]
    r["procedure"]["general_procedure_reference"]=fe("general procedure E","general procedure E",detail_evidence);r["procedure"]["original_text"]=fe(f"Compound {molecules[prod]['occurrences'][0]['compound_label']['normalized_value']} was prepared from the stated substrate and lead tetracarboxylate using General Procedure E.","specific General Procedure E application",detail_evidence)
    target_outcome=next((o for o in r["outcomes"] if prod in o["target_molecule_ids"] and o["kind"]=="yield"),None)
    if target_outcome is None:
        yield_outcome(r,"yield",prod,yld,detail_evidence,scope=f"preparation of {molecules[prod]['occurrences'][0]['compound_label']['normalized_value']}")
        target_outcome=r["outcomes"][-1]
    else:
        extend_evidence(target_outcome,detail_evidence)
    target_outcome["value"]=fe(f"{yld}%",{"number":yld,"unit":"%"},detail_evidence)
    add_product_mass_basis(target_outcome,mass,detail_evidence,"General Procedure E")
    r["review"]["unresolved"]=["Machine-transcribed structure, stereochemical status, and atom-level bond classification require expert review."]
    if conflict:
        r["notes"].append(fe(conflict,conflict,detail_evidence));r["review"]["unresolved"].append(conflict)
    return r


gpe_rows=[
    ("scheme2_product_9","m_s27","m_9","0.215 mmol scale",{"number":0.215,"unit":"mmol"},{"number":40,"unit":"mg"},75,ev(8),"m_lta",None),
    ("scheme2_product_4","m_s28","m_4","0.096 mmol scale",{"number":0.096,"unit":"mmol"},{"number":15.3,"unit":"mg"},60,ev(8),"m_lta",None),
    ("scheme2_product_3","m_s29","m_3","0.088 mmol scale",{"number":0.088,"unit":"mmol"},{"number":17.4,"unit":"mg"},71,ev(8),"m_lta",None),
    ("scheme2_product_2","m_1","m_2","0.062 mmol scale",{"number":0.062,"unit":"mmol"},{"number":16.2,"unit":"mg"},81,ev(8),"m_lta",None),
    ("scheme2_product_5","m_s30","m_5","0.167 mmol scale",{"number":0.167,"unit":"mmol"},{"number":18.6,"unit":"mg"},61,ev(8)+ev(9),"m_lta",None),
    ("scheme2_product_6","m_s31","m_6","0.100 mmol scale",{"number":0.100,"unit":"mmol"},{"number":21.3,"unit":"mg"},72,ev(9),"m_lta","The product-6 paragraph names the starting material as methyl [4-methoxy-2-(propa-1,2-dien-1-yl)phenyl]carbamate while assigning label S31; the S31 preparation paragraph identifies the tert-butyl substrate. The record follows the explicit label S31 and product structure but preserves this source conflict."),
    ("scheme2_product_7","m_s32","m_7","0.187 mmol scale",{"number":0.187,"unit":"mmol"},{"number":33.2,"unit":"mg"},64,ev(9),"m_lta",None),
    ("scheme2_product_10","m_s34","m_10","0.113 mmol scale",{"number":0.113,"unit":"mmol"},{"number":23.8,"unit":"mg"},73,ev(9),"m_lta",None),
    ("scheme2_product_11","m_s35","m_11","0.100 mmol scale",{"number":0.100,"unit":"mmol"},{"number":21.1,"unit":"mg"},65,ev(9),"m_lta",None),
    ("scheme2_product_12","m_s36","m_12","0.077 mmol scale",{"number":0.077,"unit":"mmol"},{"number":14.9,"unit":"mg"},56,ev(9),"m_lta",None),
    ("scheme4_product_23","m_19","m_23","0.197 mmol scale",{"number":0.197,"unit":"mmol"},{"number":25,"unit":"mg"},49,ev(9),"m_lta",None),
    ("scheme4_product_22","m_18","m_22","0.071 mmol scale",{"number":0.071,"unit":"mmol"},{"number":11.6,"unit":"mg"},48,ev(9),"m_lta",None),
    ("scheme4_product_24","m_20","m_24","0.059 mmol scale",{"number":0.059,"unit":"mmol"},{"number":3,"unit":"mg"},14,ev(9),"m_lta",None),
    ("isotope_kie_35_to_36","m_35","m_36","0.110 mmol scale",{"number":0.110,"unit":"mmol"},{"number":18.7,"unit":"mg"},68,ev(9),"m_lta",None),
    ("scheme3_product_13","m_s27","m_13","0.108 mmol scale",{"number":0.108,"unit":"mmol"},{"number":20.1,"unit":"mg"},57,ev(9),"m_pb_obz4",None),
    ("scheme3_product_14","m_s27","m_14","0.159 mmol scale",{"number":0.159,"unit":"mmol"},{"number":22.8,"unit":"mg"},44,ev(9),"m_pb_pfbz4",None),
    ("scheme3_product_16","m_s27","m_16","0.159 mmol scale",{"number":0.159,"unit":"mmol"},{"number":32.4,"unit":"mg"},54,ev(9)+ev(10),"m_pb_ptfbz4",None),
    ("scheme3_product_15","m_s27","m_15","0.160 mol scale",{"number":0.160,"unit":"mol"},{"number":16.6,"unit":"mg"},31,ev(10),"m_pb_pmbz4","The source prints a 0.160 mol scale together with 16.6 mg product and 31% yield. The apparent unit inconsistency is retained exactly and not corrected to mmol."),
    ("scheme3_product_17","m_s27","m_17","0.169 mmol scale",{"number":0.169,"unit":"mmol"},{"number":38.1,"unit":"mg"},67,ev(10),"m_pb_cinn4",None),
]
for row in gpe_rows:enrich_gpe(*row)


# Page-10 starting-material preparations.
prep26=reaction("prep_26",ev(10),categories=["other"],location="page 10 preparation of 26",name="Methyl carbamate protection of o-toluidine",transformation="methyl carbamate formation",scheme_ids=[gpA["reaction_id"]],bonds=(("C","N"),),unresolved=["Specific absolute reagent volumes were not restated; only General Procedure A equivalents and concentration are associated."])
prep26_parts=[part("o_toluidine","m_o_toluidine","reactant",ev(10),eq={"number":1,"unit":"equiv"},amount={"number":9.33,"unit":"mmol"}),part("pyridine","m_pyridine","reagent",ev(4),eq={"number":3.3,"unit":"equiv"}),part("methyl_chloroformate","m_methyl_chloroformate","reagent",ev(4),eq={"number":3,"unit":"equiv"}),part("dcm","m_dcm","solvent",ev(4)),part("water","m_water","workup_agent",ev(4),"workup"),part("brine","m_brine","workup_agent",ev(4),"workup"),part("na2so4","m_na2so4","workup_agent",ev(4),"workup"),part("etoac","m_etoac","purification_agent",ev(4),"purification"),part("hexanes","m_hexanes","purification_agent",ev(4),"purification"),part("product","m_26","product",ev(10),mass={"number":1.49,"unit":"g"},product_role="reported product")]
set_refs(prep26_parts,"m_o_toluidine");prep26["participants"]=prep26_parts;prep26["reaction_scale"]=fe("9.33 mmol o-toluidine",{"reference_molecule_id":"m_o_toluidine","amount":{"number":9.33,"unit":"mmol"}},ev(10));prep26["solvent_conditions"]["composition"]=fe("DCM",["dichloromethane"],ev(4));prep26["solvent_conditions"]["concentration"]=fe("1.0 M",{"number":1.0,"unit":"M"},ev(4));prep26["solvent_conditions"]["concentration_basis"]=fe("o-toluidine","m_o_toluidine",ev(10));prep26["steps"]=[step("addition_0c",ev(4),["o_toluidine","pyridine","methyl_chloroformate","dcm"],"add methyl chloroformate over 15 min",{"number":0,"unit":"°C"},{"number":15,"unit":"min"},"nitrogen"),step("warm_and_hold",ev(4),["o_toluidine","pyridine","methyl_chloroformate","dcm"],"allow to warm to room temperature and stir","room temperature",{"number":4,"unit":"h"},"nitrogen"),step("workup",ev(4),["water","brine","na2so4"],"water quench, DCM extraction, brine wash, Na2SO4 drying",None,None,None,"workup"),step("purification",ev(4),["etoac","hexanes","product"],"silica gel chromatography, 15:85 EtOAc/hexanes",None,None,None,"purification")];prep26["procedure"]["general_procedure_reference"]=fe("general procedure A","general procedure A",ev(10));yield_outcome(prep26,"yield","m_26",97,ev(10),scope="preparation of 26");add_product_mass_basis(prep26["outcomes"][-1],{"number":1.49,"unit":"g"},ev(10),"General Procedure A")

prep32=reaction("prep_32",ev(10),categories=["other"],location="page 10 preparation of 32",name="Stille vinylation to prepare methyl (2-vinylphenyl)carbamate",transformation="palladium-catalyzed aryl-vinyl coupling",bonds=(("C","C"),),unresolved=["The printed CuI amount (0.058 mmol) and loading (0.06 equiv) are arithmetically inconsistent relative to 0.722 mmol substrate; both source values are retained."])
prep32_parts=[part("s1","m_s1","reactant",ev(10),eq={"number":1,"unit":"equiv"},amount={"number":0.722,"unit":"mmol"}),part("cui","m_cui","catalyst",ev(10),eq={"number":0.06,"unit":"equiv"},amount={"number":0.058,"unit":"mmol"}),part("pdp4","m_pdp4","catalyst",ev(10),eq={"number":0.04,"unit":"equiv"},amount={"number":0.029,"unit":"mmol"}),part("dmf","m_dmf","solvent",ev(10)),part("vinyl_tin","m_tributyl_vinyl_tin","reactant",ev(10),eq={"number":1.20,"unit":"equiv"},amount={"number":0.87,"unit":"mmol"}),part("nh4oh","m_nh4oh_5pct","workup_agent",ev(10),"workup"),part("etoac","m_etoac","workup_agent",ev(10),"workup"),part("water","m_water","workup_agent",ev(10),"workup"),part("brine","m_brine","workup_agent",ev(10),"workup"),part("na2so4","m_na2so4","workup_agent",ev(10),"workup"),part("hexanes","m_hexanes","purification_agent",ev(10),"purification"),part("product","m_32","product",ev(10),mass={"number":124,"unit":"mg"},product_role="reported product")]
set_refs(prep32_parts,"m_s1");prep32["participants"]=prep32_parts;prep32["reaction_scale"]=fe("0.722 mmol S1",{"reference_molecule_id":"m_s1","amount":{"number":0.722,"unit":"mmol"}},ev(10));prep32["solvent_conditions"]["composition"]=fe("DMF",["DMF"],ev(10));prep32["solvent_conditions"]["concentration"]=fe("0.1 M",{"number":0.1,"unit":"M"},ev(10));prep32["solvent_conditions"]["concentration_basis"]=fe("methyl (2-iodophenyl)carbamate","m_s1",ev(10));prep32["steps"]=[step("charge_and_purge",ev(10),["s1","cui","pdp4","dmf"],"charge S1, CuI, and Pd(PPh3)4; purge with nitrogen and add DMF",None,None,"nitrogen"),step("add_vinyl_tin",ev(10),["vinyl_tin"],"add tributyl(vinyl)tin dropwise",None,None,"nitrogen"),step("reaction_hold",ev(10),["s1","cui","pdp4","dmf","vinyl_tin"],"stir reaction mixture",None,{"number":16,"unit":"h"},"nitrogen"),step("workup",ev(10),["nh4oh","etoac","water","brine","na2so4"],"quench with 5% NH4OH; extract with EtOAc; wash combined organic layers five times with water and then brine; dry over Na2SO4",None,None,None,"workup"),step("purification",ev(10),["etoac","hexanes","product"],"silica gel chromatography, 15:85 EtOAc/hexanes",None,None,None,"purification")];prep32_note="Source conflict: 0.058 mmol CuI is printed together with 0.06 equiv on a 0.722 mmol scale; neither value was recalculated or replaced.";prep32["notes"].append(fe(prep32_note,prep32_note,ev(10)));yield_outcome(prep32,"yield","m_32",97,ev(10),scope="preparation of 32");add_product_mass_basis(prep32["outcomes"][-1],{"number":124,"unit":"mg"},ev(10),"silica gel chromatography")

prep29=reaction("prep_29",ev(10),categories=["other"],location="page 10 two-step preparation of 29",name="Two-step Sonogashira coupling and desilylation to prepare 29",transformation="aryl ethynylation followed by trimethylsilyl deprotection",scheme_ids=[gpB["reaction_id"]],bonds=(("C","C"),("C","H")),unresolved=["Amounts of trimethylsilylacetylene, K2CO3, and MeOH and the desilylation time/temperature/workup are not printed in the paragraph."])
prep29_parts=[part("s1","m_s1","reactant",ev(10),eq={"number":1,"unit":"equiv"},amount={"number":5.42,"unit":"mmol"}),part("tms_acetylene","m_tms_acetylene","reactant",ev(5),eq={"number":1.2,"unit":"equiv"}),part("pd","m_pdcl2pph3","catalyst",ev(5),eq={"number":0.02,"unit":"equiv"}),part("cui","m_cui","catalyst",ev(5),eq={"number":0.03,"unit":"equiv"}),part("thf","m_thf","solvent",ev(5)),part("et3n","m_et3n","solvent",ev(5)),part("intermediate","m_silyl_alkyne_intermediate","product",ev(10),mass={"number":1.36,"unit":"g"},product_role="stage-1 isolated intermediate"),part("k2co3","m_k2co3","reagent",ev(10)),part("meoh","m_meoh","solvent",ev(10)),part("product","m_29","product",ev(10),mass={"number":856,"unit":"mg"},product_role="final product")]
set_refs(prep29_parts,"m_s1");prep29["participants"]=prep29_parts;prep29["reaction_scale"]=fe("5.42 mmol S1",{"reference_molecule_id":"m_s1","amount":{"number":5.42,"unit":"mmol"}},ev(10));prep29["solvent_conditions"]["composition"]=fe("THF/NEt3, then MeOH",["THF","NEt3","MeOH"],ev(5)+ev(10));prep29["solvent_conditions"]["mixture_ratio"]=fe("THF/NEt3 4:1",{"THF":4,"NEt3":1},ev(5));prep29["solvent_conditions"]["concentration"]=fe("0.5 M for General Procedure B stage",{"number":0.5,"unit":"M"},ev(5));prep29["solvent_conditions"]["concentration_basis"]=fe("S1","m_s1",ev(10));prep29["steps"]=[step("sonogashira_addition",ev(5)+ev(10),["s1","tms_acetylene","pd","cui","thf","et3n"],"perform General Procedure B with trimethylsilylacetylene; add alkyne dropwise under nitrogen","room temperature",None,"nitrogen"),step("sonogashira_hold",ev(5)+ev(10),["s1","tms_acetylene","pd","cui","thf","et3n"],"stir Sonogashira reaction","room temperature",{"minimum":3,"maximum":6,"unit":"h"},"nitrogen"),step("isolate_intermediate",ev(10),["intermediate"],"isolate brown trimethylsilyl-protected intermediate",None,None,None,"purification"),step("desilylation",ev(10),["intermediate","k2co3","meoh","product"],"react intermediate with K2CO3 in MeOH to afford 29",None,None,None)];prep29["procedure"]["general_procedure_reference"]=fe("general procedure B for stage 1","general procedure B",ev(10));yield_outcome(prep29,"stage1_yield","m_silyl_alkyne_intermediate",97,ev(10),scope="stage 1: S1 to TMS-protected alkyne intermediate");add_product_mass_basis(prep29["outcomes"][-1],{"number":1.36,"unit":"g"},ev(10),"General Procedure B");yield_outcome(prep29,"stage2_yield","m_29",90,ev(10),scope="stage 2: isolated TMS-protected intermediate to 29");add_product_mass_basis(prep29["outcomes"][-1],{"number":856,"unit":"mg"},ev(10),"K2CO3/MeOH desilylation");yield_outcome(prep29,"overall_yield","m_29",87,ev(10),scope="overall two-step yield from S1 to 29");prep29["outcomes"][-1]["measurement_basis"]=fe("87% in two steps","reported overall yield from S1 to 29",ev(10))

prep31=reaction("prep_31",ev(10),categories=["other"],location="page 10 preparation of 31",name="Preparation of methyl (2-allylphenyl)carbamate following a known procedure",transformation="aryl allylation",bonds=(("C","C"),),unresolved=["The cited procedure, reagents, conditions, workup, purification, and explicit yield-assay method are not reproduced."])
prep31["participants"]=[part("s1","m_s1","reactant",ev(10),amount={"number":1.81,"unit":"mmol"}),part("product","m_31","product",ev(10),mass={"number":174,"unit":"mg"},product_role="reported product")];prep31["reaction_scale"]=fe("1.81 mmol S1",{"reference_molecule_id":"m_s1","amount":{"number":1.81,"unit":"mmol"}},ev(10));prep31["steps"]=[step("reported_preparation",ev(10),["s1","product"],"follow a known procedure; details are not reproduced")];prep31["procedure"]["general_procedure_reference"]=fe("known procedure (reference 26)","reference 26",ev(10));yield_outcome(prep31,"yield","m_31",50,ev(10),scope="preparation of 31");add_product_mass_basis(prep31["outcomes"][-1],{"number":174,"unit":"mg"},ev(10),"cited procedure")


# The detailed page-10 control paragraphs are the same experiments as Scheme 5.
def enrich_control(suffix, sub, product, mass, yld, evidence):
    r=find_reaction(suffix);r["evidence_ids"]=list(dict.fromkeys(r["evidence_ids"]+ev(8)+evidence));r["location"]=fe("page 3 Scheme 5 and page 10 detailed control paragraph",{"pages":[3,10],"region":"Scheme 5 control merged with detailed experimental paragraph"},r["evidence_ids"]);r["general_scheme_ids"].append(gpE["reaction_id"])
    products=[(product,mass,"reported product")]
    if suffix=="control_32_products_33_34":products.append(("m_34",None,"minor product"))
    parts=[part("substrate",sub,"reactant",evidence,eq={"number":1,"unit":"equiv"},amount={"number":0.60,"unit":"mmol"}),part("lta","m_lta","reagent",evidence,eq={"number":5,"unit":"equiv"}),part("dcm","m_dcm","solvent",ev(8))]
    for i,(mid,pmass,role) in enumerate(products,1):parts.append(part(f"product_{i}",mid,"product",evidence,mass=pmass,product_role=role))
    parts += [part("alumina","m_alumina_basic","purification_agent",ev(8),"purification"),part("dcm_eluent","m_dcm","purification_agent",ev(8),"purification"),part("etoac","m_etoac","purification_agent",ev(8),"purification"),part("hexanes","m_hexanes","purification_agent",ev(8),"purification")]
    set_refs(parts,sub);r["participants"]=parts;r["reaction_scale"]=fe("0.60 mmol substrate",{"reference_molecule_id":sub,"amount":{"number":0.60,"unit":"mmol"}},evidence);r["solvent_conditions"]["composition"]=fe("dichloromethane",["dichloromethane"],ev(8));r["solvent_conditions"]["concentration"]=fe("0.30 M",{"number":0.30,"unit":"M"},ev(8));r["solvent_conditions"]["concentration_basis"]=fe("control substrate",sub,evidence);r["steps"]=[step("addition",ev(8),["substrate","lta","dcm"],"add substrate solution to 5 equiv LTA over 30 s","room temperature",{"number":30,"unit":"s"},"nitrogen"),step("reaction_hold",ev(8),["substrate","lta","dcm"],"continue stirring overnight","room temperature",{"number":15,"unit":"h","qualifier":"approx","raw_text":"approximately 15 h"},"nitrogen"),step("alumina_filtration",ev(8),["alumina","dcm_eluent"],"filter through the stated basic-alumina plug with DCM",None,None,None,"purification"),step("silica_purification",ev(8),["etoac","hexanes"]+[f"product_{i}" for i in range(1,len(products)+1)],"concentrate and immediately purify by silica gel chromatography, 20:80 EtOAc/hexanes",None,None,None,"purification")];r["procedure"]["general_procedure_reference"]=fe("modified general procedure E using 5 equivalents of LTA","modified general procedure E",evidence);target=next(o for o in r["outcomes"] if product in o["target_molecule_ids"]);extend_evidence(target,evidence);target["value"]=fe(f"{yld}%",{"number":yld,"unit":"%"},evidence);add_product_mass_basis(target,mass,evidence,"modified General Procedure E");r["review"]["unresolved"]=["The source does not state whether the reported percent is gravimetric isolated yield, although a product mass after purification is given."]

enrich_control("control_26_to_27","m_26","m_27",{"number":60,"unit":"mg"},44,ev(10))
enrich_control("control_29_to_30","m_29","m_30",{"number":24,"unit":"mg"},17,ev(10))
enrich_control("control_32_products_33_34","m_32","m_33",{"number":62,"unit":"mg"},36,ev(10))


bundle["molecules"]=list(molecules.values())
bundle["reactions"]=reactions
bundle["source_context"]=[
 {"context_id":f"{DOC}:general_experimental","document_id":DOC,"evidence_ids":ev(4),"kind":"general_experimental_context","summary":"General atmosphere, solvent drying, chromatography, spectroscopy, HRMS, and X-ray instrumentation; not a separate reaction."},
 {"context_id":f"{DOC}:kinetics_observations","document_id":DOC,"evidence_ids":ev(3)+ev(4),"kind":"mechanistic_observation","summary":"Overall rate law reported as rate = k[allene][LTA]; intermediates C and E observed by proton NMR; lead NMR gave no useful signals. Exact kinetic traces are deferred to unprovided SI."},
 {"context_id":f"{DOC}:xray_product4","document_id":DOC,"evidence_ids":ev(2),"kind":"characterization","summary":"Scheme 2 includes an X-ray crystal depiction for product 4; it is characterization evidence, not an additional synthesis."},
 {"context_id":f"{DOC}:characterization_pages4_6","document_id":DOC,"evidence_ids":ev(4)+ev(5)+ev(6),"kind":"characterization","summary":"NMR/IR/HRMS blocks for S1--S26 are retained as characterization context and do not create extra reaction records."},
 {"context_id":f"{DOC}:source_anomalies","document_id":DOC,"evidence_ids":ev(4)+ev(5)+ev(6),"kind":"source_conflict","summary":"S2 characterization as printed appears inconsistent with the named formula, and S18 mass is printed as 855 g at 4.40 mmol/83%; neither was silently repaired."},
 {"context_id":f"{DOC}:characterization_pages7_10","document_id":DOC,"evidence_ids":ev(7)+ev(8)+ev(9)+ev(10),"kind":"characterization","summary":"NMR, IR, HRMS, physical-state, freezer-solidification, molecular-formula, and mass-spectral peak data accompanying S27--S36 and compounds 1--36 are characterization context; they were not converted into extra synthesis experiments."},
 {"context_id":f"{DOC}:source_conflicts_pages9_10","document_id":DOC,"evidence_ids":ev(9)+ev(10),"kind":"source_conflict","summary":"The product-6 paragraph names a methoxy substrate while citing S31, whose own paragraph is tert-butyl; the product-15 paragraph prints a 0.160 mol scale with 16.6 mg product/31%; both conflicts are preserved without correction."},
 {"context_id":f"{DOC}:supporting_information_notice","document_id":DOC,"evidence_ids":ev(10),"kind":"document_notice","summary":"The article points to an external Supporting Information PDF, but no SI file was explicitly supplied or associated in source_mapping.json; no SI content was imported."},
 {"context_id":f"{DOC}:bibliographic_back_matter","document_id":DOC,"evidence_ids":ev(11),"kind":"bibliographic_context","summary":"Author information, acknowledgments, and references only; no new reaction arrow, table row, or original experiment was identified on page 11."},
]

# Compatibility generation also requires evidence arrays to be sets.  Some
# cross-page expressions intentionally concatenate the procedure and specific
# paragraph evidence, so normalize duplicates once without dropping lineage.
extend_evidence(bundle,[])

st.write(RUN/"submitted.json",bundle)

# Raw visual inventory retained separately from the formal machine bundle.
inventory={"document_id":DOC,"page_range":[1,11],"semantic_accuracy_verified":False,
 "record_count":len(reactions),"molecule_count":len(molecules),
 "decisions":["Merged repeated evidence for S2 (pages 4--5), S14 (pages 5--6), and the KIE experiment (pages 3--4).",
              "Did not treat NMR/IR/HRMS peak lists, molecular formulas, or X-ray depiction as additional synthesis experiments.",
              "Did not associate or download SI; missing procedure and kinetic details remain unresolved.",
              "Kept all machine records needs_review and partial.",
              "Viewed the original page images for pages 7--11 and merged detailed experimental paragraphs into their existing Scheme 2--5 records.",
              "Preserved the product-6 substrate-name/S31 conflict and the product-15 0.160 mol scale without silent correction."],
 "records_by_page":{}}
for page in range(1,12):
    marker=f":page{page}:"
    inventory["records_by_page"][str(page)]=[r["reaction_id"] for r in reactions if any(marker in e for e in r["evidence_ids"])]
st.write(RUN/"raw_model/page_inventory.json",inventory)

coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[]}
regions={
1:[("Scheme 1, right column","Four prior-literature reaction rows and one current-work generic oxidative amino-hydroxylation row were distinguished.")],
2:[("Table 1, left column","All 14 optimization rows and footnotes were transcribed as separate experiments."),("Scheme 2, right column","Products 2--12 were inventoried as 11 separate scope outcomes, including product 8 at 0%."),("X-ray inset for 4","Characterization only; no extra synthesis record.")],
3:[("Scheme 3, upper left","Five lead(IV) benzoate/cinnamate experiments, products 13--17."),("Scheme 4, lower left","Four terminal-allene substitution experiments, products 22--25."),("Scheme 5, upper right","Five control experiments; the 32 experiment has two product outcomes."),("kinetics discussion, lower right","Rate law, failed lead-NMR observation, and deuterium/KIE comparison inventoried; only the KIE transformation is a specific reaction record.")],
4:[("Scheme 6, left column","Proposed mechanism separated from the real deuterium/KIE experiment."),("General Experimental, upper right","Non-reaction laboratory and characterization context."),("General Procedure A, right column","Procedure template plus preparations S1 and S2; S2 continues on page 5.")],
5:[("left and upper-right compound paragraphs","Preparations S3--S11 plus continuation of S2; spectra retained as characterization context."),("General Procedure B, right column","Sonogashira procedure template and preparations S12--S14; S14 continues on page 6.")],
6:[("two-column compound paragraphs","Continuation of S14 and preparations S15--S26; all spectral peak lists treated as characterization context."),("S18 source line","Printed mass 855 g retained and flagged as internally inconsistent rather than repaired.")],
}
unresolved_by_page={
1:["Prior-work rows omit yields and full procedures; no external papers were imported."],
2:["Table/scope yield measurement bases are not printed; Scheme 2 omits time, concentration, and LTA loading."],
3:["Schemes 3--5 omit several loadings/times/measurement bases; products 22/23 have unassigned stereogenic centers."],
4:["Kinetic details are deferred to unprovided SI; mechanism intermediates are not all isolated; S2 characterization is internally suspect."],
5:["Literature procedures for S6 and S9--S11 are not reproduced; S14 paragraph continues to page 6."],
6:["S18 printed mass is inconsistent with scale/yield; S23/S25 stereocenters are unassigned; no SI was supplied."],
}
for page in range(1,7):
    coverage["pages"].append({"page":page,"state":"attempted",
        "regions":[{"region":a,"finding":b,"evidence_ids":ev(page)} for a,b in regions[page]],
        "record_ids":inventory["records_by_page"][str(page)],"unresolved":unresolved_by_page[page]})
st.write(RUN/"coverage/pages_0001_0006.json",coverage)

# Required coverage artifact for the current bounded review.  Page numbers are
# original PDF pages; each region was checked against the supplied page image.
coverage_7_11={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[]}
regions_7_11={
7:[("upper-left General Procedure C","General Procedure C transcribed with 0 °C DIAD addition, 5 min hold, warm to room temperature for 2 h, TFE/H2O addition and a further 2 h, workup, and chromatography."),("upper-left General Procedure D","General Procedure D transcribed with sequential NBSH, AgOTf, and TfOH additions, 35 °C/1 h hold, workup, and chromatography."),("two-column compound paragraphs","Preparations S27--S35 and compound 1 inventoried; S36 begins at the bottom and continues on page 8. Spectral peak lists were classified as characterization context.")],
8:[("upper-left continuation and allene paragraphs","S36 continuation and preparations 19, 21, 18, 20, and deuterated allene 35 inventoried; axial configuration of 18/19 remains unassigned."),("General Procedure E","General Procedure E transcribed with 1.5 equiv oxidant, 0.30 M DCM, 30 s addition, approximately 15 h hold, alumina filtration, and silica purification."),("two-column product paragraphs","Detailed experiments for products 9, 4, 3, 2, and 5 inventoried; product 5 continues to page 9.")],
9:[("two-column product paragraphs","Detailed General Procedure E experiments for products 5, 6, 7, 10, 11, 12, 23, 22, 24, 36, 13, 14, and 16 inventoried; product 16 continues on page 10."),("product 6 paragraph","The printed methoxy-substrate name conflicts with label S31 and the tert-butyl product; both readings are retained for review."),("spectral blocks","NMR/IR/HRMS and physical descriptions retained only as characterization context, not separate syntheses.")],
10:[("upper-left continuation and products 15/17","Product 16 continuation plus detailed product 15 and 17 experiments inventoried; the printed 0.160 mol scale for 15 is preserved as a source conflict."),("control-substrate preparations","Preparations 26, 32, two-step 29, and 31 inventoried; stage and overall yields for 29 remain separately scoped."),("control products","Detailed modified-General-Procedure-E experiments 26→27, 32→33 (with scheme-only minor 34), and 29→30 merged with their Scheme 5 records."),("Associated Content notice","External SI notice recorded as document context only; no SI file was supplied or associated.")],
11:[("full page back matter","Author information, acknowledgments, and references only; no new reaction arrow, table row, general procedure, or original experiment was found.")],
}
unresolved_7_11={
7:["DIAD dropwise-addition duration and warm-up ramp time are not reported.","S36 paragraph continues onto page 8."],
8:["No configuration or ee is reported for axially chiral allenes 18 and 19.","Product 5 characterization continues onto page 9."],
9:["Product 6 substrate wording conflicts with compound label S31; no silent choice was made beyond linking the explicit label to the displayed product.","Product 16 paragraph continues onto page 10."],
10:["Product 15 scale is printed as 0.160 mol and is internally inconsistent with product mass/yield; retained unchanged.","For preparation 32, 0.058 mmol CuI conflicts arithmetically with the simultaneously printed 0.06 equiv on a 0.722 mmol scale.","External Supporting Information was not supplied or associated; omitted details remain missing."],
11:["Reference details were not imported as current-paper reaction experiments."],
}
for page in range(7,12):
    coverage_7_11["pages"].append({"page":page,"state":"attempted",
        "regions":[{"region":region,"finding":finding,"evidence_ids":ev(page)} for region,finding in regions_7_11[page]],
        "record_ids":inventory["records_by_page"][str(page)],"unresolved":unresolved_7_11[page]})
st.write(RUN/"coverage/pages_0007_0011.json",coverage_7_11)

# Deterministic stereo diagnostics for every concrete candidate SMILES.
checks=[]
for m in bundle["molecules"]:
    smi=m["smiles"]["normalized_value"]
    if smi:
        checks.append({"molecule_id":m["molecule_id"],"input_smiles":smi,"diagnostic":st.stereo_inspect(smi)})
st.write(RUN/"raw_model/stereo_checks.json",checks)
print(json.dumps({"molecules":len(molecules),"reactions":len(reactions),"stereo_checks":len(checks)},ensure_ascii=False))
