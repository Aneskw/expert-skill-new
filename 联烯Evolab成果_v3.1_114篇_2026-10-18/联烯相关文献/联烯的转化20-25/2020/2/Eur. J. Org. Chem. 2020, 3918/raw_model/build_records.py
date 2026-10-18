#!/usr/bin/env python3
"""Build source-derived v3.1 records for doc_518d3edb75482abc, PDF pages 1-6.

This is a transparent extraction artifact.  It calls only the bundled record
helpers and does not alter the schema, validators, quality rules, or tools.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st

DOC = "doc_518d3edb75482abc"
P = {i: f"{DOC}:page{i}:text" for i in range(1, 9)}
IMG = {i: f"{DOC}:page{i}:image" for i in range(1, 7)}


def with_page_images(evidence):
    """Pair page text with the supplied page image when that image exists."""
    expanded = []
    for evidence_id in evidence:
        if evidence_id not in expanded:
            expanded.append(evidence_id)
        for page, text_id in P.items():
            if evidence_id == text_id and page in IMG and IMG[page] not in expanded:
                expanded.append(IMG[page])
                break
    return expanded


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    evidence = with_page_images(evidence)
    if normalized is None and raw is not None and status != "conflict":
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def missing():
    return st.fact()


def qraw(value):
    """Render a normalized quantity back to source-like text, never Python repr."""
    if not isinstance(value, dict):
        return str(value)
    if value.get("raw_text"):
        return value["raw_text"]
    unit = value.get("unit", "")
    spacer = "" if unit in {"%", "mol%"} else " "
    if "number" in value:
        return f"{value['number']:g}{spacer}{unit}".strip()
    if "minimum" in value and "maximum" in value:
        return f"{value['minimum']:g}–{value['maximum']:g}{spacer}{unit}".strip()
    return json.dumps(value, ensure_ascii=False)


def inferred(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="inference", derivation=derivation)


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = json.loads((RUN / "raw_model/submitted_init.json").read_text())
assert bundle["documents"] == dataset["documents"]


def molecule(mid, label, name, ev, smiles=None, chirality="unspecified",
             geometry="not_applicable", representation="partial", mixture=None,
             stereo_note=None, identity_note=None):
    ev = with_page_images(ev)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f([], [], ev)
    m["names"]["aliases"] = f([], [], ev)
    m["cas_numbers"] = missing()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source": "structure drawing or printed compound name", "candidate_smiles": smiles},
            smiles,
            ev,
            method="structure_derivation",
            derivation="Candidate connectivity transcribed from the cited source drawing/name; parsing is not source verification.",
        )
        m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = f(
        stereo_note or "No absolute stereochemical assignment is printed for this identity.",
        chirality, ev)
    m["stereochemistry"]["geometry_status"] = f(
        "Alkene/allene geometry status retained from the cited drawing and ratio labels.",
        geometry, ev)
    m["stereochemistry"]["representation_status"] = f(
        "Candidate connectivity does not encode unassigned centers/axes or mixture components.",
        representation, ev)
    if chirality in {"unspecified", "partially_specified", "mixture"}:
        m["stereochemistry"]["atom_centers"] = f(
            "No R/S label or wedge/dash assignment is printed for the product alpha carbon.",
            [{"site": "alpha carbon attached to N and P", "configuration": "unspecified"}], ev)
    if geometry == "mixture":
        m["stereochemistry"]["double_bonds"] = f(
            "The source reports an E/Z or stereoisomer mixture without isolating both configured identities.",
            [{"site": "product alkene", "configuration": "mixture"}], ev)
    else:
        m["stereochemistry"]["double_bonds"] = missing()
    if mixture:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, ev)
    else:
        m["stereochemistry"]["mixture"] = missing()
    m["stereochemistry"]["axial"] = missing()
    m["stereochemistry"]["assignment"] = f(stereo_note, stereo_note, ev) if stereo_note else missing()
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = missing()
    m["identity_note"] = f(
        identity_note or "Machine candidate; CAS not reported and stereochemical silence was not converted to achiral or racemic.",
        identity_note or "Machine candidate; CAS unreported; stereochemistry retained as unresolved/unspecified.", ev)
    m["external_properties"] = []
    return m


def simple_agent(mid, name, smiles, ev):
    m = molecule(mid, name, name, ev, smiles=smiles, chirality="achiral",
                 geometry="not_applicable", representation="not_applicable",
                 stereo_note="No stereogenic element is present in the candidate molecular graph.")
    m["stereochemistry"]["atom_centers"] = missing()
    m["stereochemistry"]["assignment"] = inferred(
        "achiral candidate graph", ev,
        "Derived only from the candidate graph; does not certify source identity.")
    return m


def generic_molecule(mid, label, fixed, variables, ev, name):
    ev = with_page_images(ev)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(name, name, ev)
    m["names"]["english"] = f(name, name, ev)
    m["names"]["abbreviations"] = f([], [], ev)
    m["names"]["aliases"] = f([], [], ev)
    m["cas_numbers"] = missing()
    m["representation_type"] = "generic_template"
    m["smiles"] = missing()
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, ev),
        "variables": [
            {
                "label": label_,
                "attachment_points": f(defn.get("attachment"), defn.get("attachment"), ev),
                "definition": f(defn.get("definition"), defn.get("definition"), ev),
                "allowed_values": f(defn.get("allowed"), defn.get("allowed"), ev),
                "constraints": f(defn.get("constraints"), defn.get("constraints"), ev),
            }
            for label_, defn in variables.items()
        ],
        "instance_mapping": [],
        "reason": f("Scheme uses R/EWG/OR′ variables.", "Source generic structure cannot be a concrete SMILES.", ev),
    }
    for key in ("chirality_status", "geometry_status", "representation_status"):
        m["stereochemistry"][key] = f(
            "Generic drawing does not specify all substituents or configurations.", "unresolved", ev)
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing()
    m["identity_note"] = f("Generic Scheme 1 identity; variables retained.",
                            "Generic Scheme 1 identity; variables retained.", ev)
    m["external_properties"] = []
    return m


products = {
    "m5": ("5", "Diethyl 1-(4-methoxybenzylamino)-3-methylbut-2-enylphosphonate",
           "COc1ccc(CNC(P(=O)(OCC)OCC)C=C(C)C)cc1", [P[2], P[3], P[5]], "unspecified", "not_applicable", None),
    "m12": ("12", "Diethyl [1-(benzylamino)-3-methylbut-2-en-1-yl]phosphonate",
            "c1ccccc1CNC(P(=O)(OCC)OCC)C=C(C)C", [P[2], P[3], P[5]], "unspecified", "not_applicable", None),
    "m13": ("13", "Diethyl 1-[N-(4-methoxybenzyl)methylsulfonamido]-3-methylbut-2-enylphosphonate",
            "COc1ccc(CN(S(=O)(=O)C)C(P(=O)(OCC)OCC)C=C(C)C)cc1", [P[2], P[3], P[5]], "unspecified", "not_applicable", None),
    "m14": ("14", "tert-Butyl 1-(diethoxyphosphoryl)-3-methylbut-2-enyl(4-methoxybenzyl)carbamate",
            "COc1ccc(CN(C(=O)OC(C)(C)C)C(P(=O)(OCC)OCC)C=C(C)C)cc1", [P[3], P[5]], "unspecified", "not_applicable", None),
    "m15": ("15", "Diethyl [1-(diethoxyphosphoryl)-3-methylbut-2-en-1-yl]-(4-methoxybenzyl)phosphoramidate",
            "COc1ccc(CN(P(=O)(OCC)OCC)C(P(=O)(OCC)OCC)C=C(C)C)cc1", [P[3], P[5]], "unspecified", "not_applicable", None),
    "m16": ("16", "Dimethyl 1-(4-methoxybenzylamino)-3-methylbut-2-enylphosphonate",
            "COc1ccc(CNC(P(=O)(OC)OC)C=C(C)C)cc1", [P[3], P[5]], "unspecified", "not_applicable", None),
    "m17": ("17", "1-(5,5-Dimethyl-2-oxido-1,3,2-dioxaphosphinan-2-yl)-N-(4-methoxybenzyl)-3-methylbut-2-en-1-amine",
            "COc1ccc(CNC(P2(=O)OCC(C)(C)CO2)C=C(C)C)cc1", [P[3], P[5]], "unspecified", "not_applicable", None),
    "m18": ("18", "Ethyl phenyl {1-[(4-methoxybenzyl)amino]-3-methylbut-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)(Oc2ccccc2))C=C(C)C)cc1", [P[3], P[5]], "mixture", "not_applicable", "1:1 mixture of two inseparable diastereoisomers"),
    "m19": ("19", "(E,Z)-Diethyl {1-[(4-methoxybenzyl)amino]but-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=CC)cc1", [P[3], P[5], P[6]], "unspecified", "mixture", "E/Z = 1:1"),
    "m20": ("20", "Diethyl 1-(4-methoxybenzylamino)prop-2-en-1-ylphosphonate (drawn; decomposed)",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C)cc1", [P[3]], "unspecified", "unspecified", None),
    "m21": ("21", "Boc/PMB-protected dimethyl 1-aminoprop-2-en-1-ylphosphonate (drawn; decomposed)",
            "COc1ccc(CN(C(=O)OC(C)(C)C)C(P(=O)(OC)OC)C=C)cc1", [P[3]], "unspecified", "unspecified", None),
    "m22": ("22", "Diethyl {2-cyclopentylidene-1-[(4-methoxybenzyl)amino]ethyl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C2CCCC2)cc1", [P[3], P[6]], "unspecified", "not_applicable", None),
    "m23": ("23", "Diethyl {2-cyclohexylidene-1-[(4-methoxybenzyl)amino]ethyl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C2CCCCC2)cc1", [P[3], P[6]], "unspecified", "not_applicable", None),
    "m24": ("24", "Diethyl {1-[(4-methoxybenzyl)amino]-3-methylpent-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C(C)CC)cc1", [P[3], P[6]], "unspecified", "mixture", "E/Z = 1:1"),
    "m25": ("25", "Diethyl {1-[(4-methoxybenzyl)amino]-3-phenylbut-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C(C)c2ccccc2)cc1", [P[3], P[4], P[6]], "unspecified", "mixture", "E/Z = 9:1, E major"),
    "m25prime": ("25′", "Diethyl (3-phenylbut-2-en-1-yl)phosphonate",
                 "CCOP(=O)(CC=C(C)c1ccccc1)OCC", [P[4], P[6]], "achiral", "mixture", "E/Z = 9:1, E major"),
    "m26": ("26", "Diethyl {3-cyclopropyl-1-[(4-methoxybenzyl)amino]but-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C(C)C2CC2)cc1", [P[3], P[6]], "unspecified", "mixture", "E/Z = 1:1"),
    "m27": ("27", "Diethyl {4-(3,4-dimethoxyphenyl)-1-[(4-methoxybenzyl)amino]-3-methylbut-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C(C)Cc2ccc(OC)c(OC)c2)cc1", [P[3], P[6]], "unspecified", "mixture", "E/Z = 1:1"),
    "m28": ("28", "Diethyl 1-(2-methoxybenzylamino)-3-methylbut-2-enylphosphonate",
            "COc1ccccc1CNC(P(=O)(OCC)OCC)C=C(C)C", [P[4], P[6]], "unspecified", "not_applicable", None),
    "m29": ("29", "Diethyl [1-(allylamino)-3-methylbut-2-en-1-yl]phosphonate",
            "C=CCNC(P(=O)(OCC)OCC)C=C(C)C", [P[4]], "unspecified", "not_applicable", None),
    "m30": ("30", "Diethyl [1-(isobutylamino)-3-methylbut-2-en-1-yl]phosphonate",
            "CC(C)CNC(P(=O)(OCC)OCC)C=C(C)C", [P[4]], "unspecified", "not_applicable", None),
    "m31": ("31", "Diethyl {1-[(3,4-dimethoxybenzyl)amino]-3-methylbut-2-en-1-yl}phosphonate",
            "COc1ccc(CNC(P(=O)(OCC)OCC)C=C(C)C)cc1OC", [P[4]], "unspecified", "not_applicable", None),
    "m32": ("32", "Diethyl 1-[(4-trifluoromethylbenzyl)amino]-3-methylbut-2-enylphosphonate",
            "FC(F)(F)c1ccc(CNC(P(=O)(OCC)OCC)C=C(C)C)cc1", [P[4]], "unspecified", "not_applicable", None),
    "m33": ("33", "Diethyl 3-methyl-1-(phenylamino)but-2-enylphosphonate",
            "c1ccc(NC(P(=O)(OCC)OCC)C=C(C)C)cc1", [P[4]], "unspecified", "not_applicable", None),
    "m34": ("34", "Diethyl [1-({2-[(tert-butyldimethylsilyl)oxy]ethyl}amino)-3-methylbut-2-en-1-yl]phosphonate",
            "CC(C)(C)[Si](C)(C)OCCNC(P(=O)(OCC)OCC)C=C(C)C", [P[4]], "unspecified", "not_applicable", None),
}

for mid, (label, name, smiles, ev, chiral, geom, mix) in products.items():
    m = molecule(mid, label, name, ev, smiles, chiral, geom,
                 "partial" if chiral not in {"achiral"} or geom == "mixture" else "not_applicable",
                 mixture=mix)
    if mid == "m18":
        m["stereochemistry"]["atom_centers"] = f(
            "Compound 18 is reported as a 1:1 mixture of two inseparable diastereoisomers; no configuration is assigned.",
            [
                {"site": "alpha carbon attached to N and P", "configuration": "unspecified"},
                {"site": "phosphorus bearing C, =O, OEt and OPh substituents", "configuration": "unspecified"},
            ], ev)
    bundle["molecules"].append(m)

# Explicitly drawn/concrete allene precursors.  Allene axial configuration is left
# unspecified where substituent patterns permit it; no R_a/S_a is invented.
precursor_smiles = {
    "m4": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)C)cc1",
    "m4a": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)c2ccccc2)cc1",
    "pre12": "c1ccccc1CN(S(=O)(=O)c2ccc(OC)cc2)C(P(=O)(OCC)OCC)=C=C(C)C",
    "pre13": "COc1ccc(CN(S(=O)(=O)C)C(P(=O)(OCC)OCC)=C=C(C)C)cc1",
    "pre14": "COc1ccc(CN(C(=O)OC(C)(C)C)C(P(=O)(OCC)OCC)=C=C(C)C)cc1",
    "pre15": "COc1ccc(CN(P(=O)(OCC)OCC)C(P(=O)(OCC)OCC)=C=C(C)C)cc1",
    "pre16": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OC)OC)=C=C(C)C)cc1",
    "pre17": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P3(=O)OCC(C)(C)CO3)=C=C(C)C)cc1",
    "pre18": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)(Oc3ccccc3))=C=C(C)C)cc1",
    "pre19": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=CC)cc1",
    "pre20": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C)cc1",
    "pre21": "COc1ccc(CN(C(=O)OC(C)(C)C)C(P(=O)(OC)OC)=C=C)cc1",
    "pre22": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C3CCCC3)cc1",
    "pre23": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C3CCCCC3)cc1",
    "pre24": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)CC)cc1",
    "pre26": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)C3CC3)cc1",
    "pre27": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)Cc3ccc(OC)c(OC)c3)cc1",
    "pre28": "COc1ccccc1CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)C",
    "pre29": "C=CCN(S(=O)(=O)c1ccc(C)cc1)C(P(=O)(OCC)OCC)=C=C(C)C",
    "pre30": "CC(C)CN(S(=O)(=O)c1ccc(C)cc1)C(P(=O)(OCC)OCC)=C=C(C)C",
    "pre31": "COc1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)C)cc1OC",
    "pre32": "FC(F)(F)c1ccc(CN(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)C)cc1",
    "pre33": "c1ccc(N(S(=O)(=O)c2ccc(C)cc2)C(P(=O)(OCC)OCC)=C=C(C)C)cc1",
    "pre34": "CC(C)(C)[Si](C)(C)OCCN(S(=O)(=O)c1ccc(C)cc1)C(P(=O)(OCC)OCC)=C=C(C)C",
}

precursor_evidence = {
    "m4": [P[2], P[3], P[5]], "m4a": [P[3], P[4]],
    "pre12": [P[2], P[3]], "pre13": [P[2], P[3]], "pre14": [P[3]], "pre15": [P[3]],
    "pre16": [P[3]], "pre17": [P[3]], "pre18": [P[3]], "pre19": [P[3]],
    "pre20": [P[3]], "pre21": [P[3]], "pre22": [P[3]], "pre23": [P[3]],
    "pre24": [P[3]], "pre26": [P[3]], "pre27": [P[3]], "pre28": [P[4]],
    "pre29": [P[4]], "pre30": [P[4]], "pre31": [P[4]], "pre32": [P[4]],
    "pre33": [P[4]], "pre34": [P[4]],
}
axially_chiral_precursors = {"m4a", "pre19", "pre24", "pre26", "pre27"}
phosphorus_chiral_precursors = {"pre18"}

for mid, smiles in precursor_smiles.items():
    label = "4" if mid == "m4" else "4a" if mid == "m4a" else f"precursor to {mid[3:]}"
    name = "N-tosyl alpha-amino allenylphosphonate 4" if mid == "m4" else (
        "N-tosyl alpha-amino allenylphosphonate 4a" if mid == "m4a" else
        f"Corresponding alpha-amino allenylphosphonate precursor to product {mid[3:]}")
    has_unassigned_stereogenic_element = mid in axially_chiral_precursors | phosphorus_chiral_precursors
    m = molecule(
        mid, label, name, precursor_evidence[mid], smiles,
        chirality="unspecified" if has_unassigned_stereogenic_element else "achiral",
        geometry="not_applicable",
        representation="partial" if has_unassigned_stereogenic_element else "not_applicable",
        stereo_note=(
            "Allene axial configuration is not assigned in the cited source drawing."
            if mid in axially_chiral_precursors else
            "Phosphorus configuration is not assigned in the cited source drawing."
            if mid in phosphorus_chiral_precursors else None
        ),
        identity_note="Connectivity candidate transcribed from the scope drawing; allene axial configuration remains unspecified where applicable.",
    )
    # The allene terminal carbon is sp2-like and must not be represented as the
    # tetrahedral alpha-carbon center used for the allyl products.
    if mid in phosphorus_chiral_precursors:
        m["stereochemistry"]["atom_centers"] = f(
            "The unsymmetrical phosphonate substituents make phosphorus stereogenic; the source assigns no configuration.",
            [{"site": "phosphorus bearing allene carbon, =O, OEt and OPh substituents", "configuration": "unspecified"}],
            precursor_evidence[mid], method="structure_derivation",
            derivation="Potential phosphorus stereogenicity follows from the four different substituent paths in the candidate graph; source assignment is absent.")
    else:
        m["stereochemistry"]["atom_centers"] = missing()
    if mid in axially_chiral_precursors:
        m["stereochemistry"]["axial"] = f(
            "No R_a/S_a or P/M label is printed for the unsymmetrically substituted allene.",
            [{"site": "allene C=C=C axis", "configuration": "unspecified"}],
            precursor_evidence[mid], method="structure_derivation",
            derivation="Both terminal allene carbons have two different substituent paths in the candidate graph, so axial stereochemistry is possible but unassigned.")
    else:
        m["stereochemistry"]["axial"] = f(
            "The candidate allene has two identical substituent paths at one terminus; no stereogenic allene axis results.",
            [{"site": "allene C=C=C axis", "configuration": "not_applicable"}],
            precursor_evidence[mid], method="structure_derivation",
            derivation="Allene-axis applicability assessed from the evidence-linked candidate connectivity; this does not independently certify the structure.")
    if not has_unassigned_stereogenic_element:
        m["stereochemistry"]["assignment"] = inferred(
            "achiral candidate graph", precursor_evidence[mid],
            "No tetrahedral, phosphorus, or allene-axis stereogenic element remains in the candidate graph; source identity still needs review.")
    bundle["molecules"].append(m)

# Generic Scheme 1 identities.
variables = {
    "EWG": {"attachment": "on nitrogen", "definition": "electron-withdrawing N-protecting group", "allowed": "not enumerated in Scheme 1", "constraints": "source-defined"},
    "R": {"attachment": "on nitrogen", "definition": "N-substituent", "allowed": "not enumerated in Scheme 1", "constraints": "source-defined"},
    "R1": {"attachment": "allene terminus", "definition": "allene substituent", "allowed": "not enumerated in Scheme 1", "constraints": "source-defined"},
    "R2": {"attachment": "allene terminus", "definition": "allene substituent", "allowed": "not enumerated in Scheme 1", "constraints": "source-defined"},
    "Rprime": {"attachment": "phosphonate oxygens", "definition": "phosphonate ester substituent", "allowed": "not enumerated in Scheme 1", "constraints": "same or independently variable; source does not specify"},
}
bundle["molecules"] += [
    generic_molecule("generic_1", "1", "EWG/R substituted alpha-amino allenylphosphonate core with P(=O)(OR′)2", variables, [P[1]], "generic alpha-amino allenylphosphonate 1"),
    generic_molecule("generic_2", "2", "corresponding alpha-amino allylphosphonate core", variables, [P[1]], "generic alpha-amino allylphosphonate 2"),
    generic_molecule("generic_3", "3", "corresponding alpha-amino vinylphosphonate core", variables, [P[1]], "generic alpha-amino vinylphosphonate 3"),
]

# Reagents, solvents, work-up, and purification materials explicitly used in pages 1-6.
agents = {
    "a_et4nbf4": ("Tetraethylammonium tetrafluoroborate", "[B-](F)(F)(F)F.CC[N+](CC)(CC)CC", [P[2], P[3], P[4], P[5]]),
    "a_mecn": ("Acetonitrile", "CC#N", [P[2], P[3], P[4], P[5]]),
    "a_dmf": ("N,N-Dimethylformamide", "CN(C)C=O", [P[2], P[3], P[4], P[5]]),
    "a_malonate": ("Diethyl malonate", "CCOC(=O)CC(=O)OCC", [P[2], P[3], P[4], P[5]]),
    "a_succinimide": ("Succinimide", "O=C1CCC(=O)N1", [P[2]]),
    "a_phenol": ("Phenol", "Oc1ccccc1", [P[2]]),
    "a_benzoic": ("Benzoic acid", "O=C(O)c1ccccc1", [P[2]]),
    "a_nbu4nclo4": ("Tetrabutylammonium perchlorate", "CCCC[N+](CCCC)(CCCC)CCCC.[O-]Cl(=O)(=O)=O", [P[2]]),
    "a_liclo4": ("Lithium perchlorate", "[Li+].[O-]Cl(=O)(=O)=O", [P[2]]),
    "a_argon": ("Argon", "[Ar]", [P[4], P[5]]),
    "a_etoac": ("Ethyl acetate", "CCOC(C)=O", [P[4], P[5], P[6]]),
    "a_silica": ("Silica gel", "O=[Si]=O", [P[4], P[5], P[6]]),
    "a_celite": ("Celite", "O=[Si]=O", [P[5]]),
    "a_pe": ("Petroleum ether (composition not specified)", None, [P[4], P[5], P[6]]),
    "a_dcm": ("Dichloromethane", "ClCCl", [P[5], P[6]]),
    "a_meoh": ("Methanol", "CO", [P[5], P[6]]),
    "a_nh3": ("Aqueous ammonia used in chromatography eluent", "N", [P[6]]),
    "a_h2": ("Hydrogen", "[H][H]", [P[1]]),
    "a_pdc": ("Poisoned palladium catalyst (Pd/C, quinoline)", None, [P[1]]),
    "a_ethanol": ("Ethanol", "CCO", [P[1]]),
}
for mid, (name, smiles, ev) in agents.items():
    if smiles:
        bundle["molecules"].append(simple_agent(mid, name, smiles, ev))
    else:
        bundle["molecules"].append(molecule(mid, name, name, ev, identity_note="Named source material with unresolved composition/structure; no concrete SMILES asserted."))


def reaction(rid, ev, record_type="specific_experiment", categories=None, location=None,
             name="Electrochemical reduction of an alpha-amino allenylphosphonate",
             transformation="chemoselective partial reduction of an allene to an alkene",
             unresolved=None):
    ev = with_page_images(ev)
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories or [], categories or [], ev)
    r["location"] = f(location or rid, location or rid, ev)
    r["reaction_name"] = f(name, name, ev)
    r["transformation_type"] = inferred(
        transformation, ev,
        "Transformation description inferred from the depicted reactant/product relationship and article wording.",
        raw="source structures and article wording")
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": list(unresolved or ["Candidate structures and bond changes require expert review."])
    }
    r["completeness"] = "partial"
    return r


def add_participant(r, pid, mid, role, ev, stage="reaction", amount=None, mass=None,
                    volume=None, equivalents=None, mol_percent=None, reference=None,
                    product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    for key, val in (("amount", amount), ("mass", mass), ("volume", volume),
                     ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if val is not None:
            p["charges"][key] = f(qraw(val), val, ev)
    if reference is not None:
        p["equivalence_reference"] = f(reference, reference, ev)
    r["participants"].append(p)
    return pid


def add_step(r, sid, stage, pids, ev, operation, temperature=None, time=None, atmosphere=None):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = f(qraw(temperature), temperature, ev)
    if time is not None:
        s["time"] = f(qraw(time), time, ev)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, ev, value, role=None, yield_type=None,
            basis=None, scope=None, value_fact=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["product_role"] = f(role, role, ev) if role else missing()
    o["value"] = value_fact if value_fact is not None else f(qraw(value), value, ev)
    o["yield_type"] = f(yield_type, yield_type, ev) if yield_type else missing()
    o["measurement_basis"] = f(basis, basis, ev) if basis else missing()
    o["scope"] = f(scope, scope, ev) if scope else missing()
    r["outcomes"].append(o)
    return o


def classification(r, mode):
    ev = r["evidence_ids"]
    r["classification"]["mapped_reaction_smiles"] = missing()
    r["classification"]["model_votes"] = []
    if mode == "unresolved":
        r["classification"]["bond_changes"] = []
        r["classification"]["labels"] = missing()
        r["classification"]["coverage"] = inferred(
            "unresolved", ev, "Observed products or atom correspondence are insufficient for a safe bond-change assignment.")
        r["classification"]["limitations"] = ["No atom mapping asserted.", "Unknown decomposition or incompletely reported outcome."]
        return
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    if mode == "detosyl":
        changes += [
            {"change_id": "b4", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "b5", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "S"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ]
    r["classification"]["bond_changes"] = [
        f(None, ch, ev, method="structure_derivation",
          derivation="Net graph change inferred from the source allene/product drawings; atom mapping and proton provenance are not established.")
        for ch in changes
    ]
    labels = [{"category_number": 1, "code": "C-H", "supporting_change_ids": ["b1", "b2"]}]
    if mode == "detosyl":
        labels.append({"category_number": 29, "code": "N-H", "supporting_change_ids": ["b4"]})
    r["classification"]["labels"] = f(
        None, labels, ev, method="structure_derivation",
        derivation="Taxonomy labels include only inferred formed C-H/N-H bonds; allene bond-order decrease and N-S cleavage are retained but not labeled as formed bonds.")
    r["classification"]["coverage"] = inferred(
        "partial", ev, "Hydrogen gains are net changes; complete atom mapping and proton sources are unresolved.")
    r["classification"]["limitations"] = [
        "No complete atom mapping asserted.",
        "Hydrogen gains are net graph changes; proton donor identity per hydrogen is not established.",
        "The reduced C=C bond is an order decrease, not a newly formed C-C sigma bond.",
    ]


def set_scheme_conditions(r, substrate_mid, product_mids, ev, potential=-2.50,
                          fmol=4, solvent="mixed", malonate=True, malonate_equiv=4.5,
                          detailed=False, product_roles=None):
    product_roles = product_roles or {}
    sub_ev = [P[5]] if detailed else ev
    add_participant(r, "p_substrate", substrate_mid, "reactant", sub_ev,
                    amount={"number": 0.215, "unit": "mmol"} if detailed else None)
    for index, mid in enumerate(product_mids, 1):
        add_participant(r, f"p_product_{index}", mid, "product", ev,
                        product_role=product_roles.get(mid, "product"))
    add_participant(r, "p_electrolyte", "a_et4nbf4", "additive", ev)
    if solvent == "mixed":
        add_participant(r, "p_mecn", "a_mecn", "solvent", ev)
        add_participant(r, "p_dmf", "a_dmf", "solvent", ev)
        r["solvent_conditions"]["composition"] = f("anhydrous CH3CN/DMF", ["acetonitrile", "DMF"], ev)
        r["solvent_conditions"]["mixture_ratio"] = f("1:1", "1:1", ev)
    elif solvent == "mecn":
        add_participant(r, "p_mecn", "a_mecn", "solvent", ev)
        r["solvent_conditions"]["composition"] = f("CH3CN", ["acetonitrile"], ev)
        r["solvent_conditions"]["mixture_ratio"] = f("pure CH3CN", "not_applicable", ev)
    elif solvent == "dmf":
        add_participant(r, "p_dmf", "a_dmf", "solvent", ev)
        r["solvent_conditions"]["composition"] = f("DMF", ["DMF"], ev)
        r["solvent_conditions"]["mixture_ratio"] = f("pure DMF", "not_applicable", ev)
    r["solvent_conditions"]["concentration"] = f("Et4NBF4 0.1 M", {"number": 0.1, "unit": "M"}, ev)
    r["solvent_conditions"]["concentration_basis"] = f("Et4NBF4", "a_et4nbf4", ev)
    if detailed:
        r["solvent_conditions"]["total_volume"] = f("25 mL", {"number": 25, "unit": "mL"}, [P[5]])
        r["reaction_scale"] = f("N-protected alpha-AAP (0.215 mmol)",
                                  {"reference_molecule_id": substrate_mid, "amount": {"number": 0.215, "unit": "mmol"}}, [P[5]])
    if malonate:
        add_participant(r, "p_malonate", "a_malonate", "reagent", ev,
                        amount={"number": 0.966, "unit": "mmol"} if detailed and malonate_equiv == 4.5 else None,
                        equivalents={"number": malonate_equiv, "unit": "equiv"}, reference=substrate_mid)
    if detailed:
        add_participant(r, "p_argon", "a_argon", "reagent", [P[5]])
        add_participant(r, "p_etoac", "a_etoac", "workup_agent", [P[5]], stage="workup")
        add_participant(r, "p_celite", "a_celite", "workup_agent", [P[5]], stage="workup")
        add_participant(r, "p_silica", "a_silica", "purification_agent", [P[5]], stage="purification")
        add_step(r, "s1_setup", "reaction", ["p_substrate"] + (["p_malonate"] if malonate else []) + ["p_electrolyte"], [P[5]],
                 "Add substrate and proton source to electrolyte solution in a divided cell.")
        add_step(r, "s2_degas", "reaction", ["p_argon"], [P[5]], "Degas under argon before electrolysis.",
                 time={"number": 15, "unit": "min"}, atmosphere="argon")
        add_step(r, "s3_electrolysis", "reaction", ["p_substrate", "p_electrolyte"] + (["p_malonate"] if malonate else []), [P[5]],
                 "Apply constant potential; maintain argon until current is less than 5% of its initial value.",
                 temperature="room temperature", time={"number": 5, "unit": "h"}, atmosphere="argon")
        add_step(r, "s4_remove_solvent", "workup", [], [P[5]], "Remove solvent under vacuum.")
        add_step(r, "s5_wash_filter", "workup", ["p_etoac", "p_celite"], [P[5]],
                 "Wash residue with ethyl acetate and filter the suspension through Celite.")
        add_step(r, "s6_evaporate", "workup", [], [P[5]], "Evaporate solvent under vacuum.")
        add_step(r, "s7_chromatography", "purification", ["p_silica"], [P[5]], "Purify the residue by silica-gel chromatography.")
        r["procedure"]["original_text"] = f(
            "General Procedure for the Electrochemical Reduction of the alpha-AAP, page 5.",
            "General Procedure A as printed on page 5.", [P[5]])
        r["procedure"]["general_procedure_reference"] = f("general procedure A", "general procedure A", [P[5]])
        r["procedure"]["addition_order"] = f(
            "N-protected alpha-AAP and diethyl malonate were added to the electrolyte solution; degassed before electrolysis.",
            ["substrate", "diethyl malonate", "electrolyte solution", "argon degassing", "electrolysis"], [P[5]])
        r["procedure"]["pretreatment"] = f("degassed under argon for 15 min", "argon degassing for 15 min", [P[5]])
        r["special_conditions"]["electrodes"] = f(
            "GC working electrode ca. 120 x 10 x 2 mm3; SCE reference; Pt-gauze counter electrode; divided cell",
            {"working": "glassy carbon stick ca. 120 x 10 x 2 mm", "reference": "SCE", "counter": "Pt gauze"}, [P[4], P[5]])
        r["special_conditions"]["reactor"] = f("divided cell", "divided electrochemical cell", [P[5]])
    else:
        add_step(r, "s1_electrolysis", "reaction", ["p_substrate", "p_electrolyte"] + (["p_malonate"] if malonate else []), ev,
                 "Controlled-potential electrolysis in a divided cell.", temperature="room temperature")
        r["special_conditions"]["reactor"] = f("divided cell", "divided electrochemical cell", ev)
    r["special_conditions"]["voltage"] = f(f"{potential:.2f} V vs. SCE", {"number": potential, "unit": "V", "reference": "SCE"}, ev)
    r["special_conditions"]["current"] = f(f"{fmol} F/mol", {"number": fmol, "unit": "F/mol"}, ev)


def isolated_basis(ev):
    return "isolated; yields refer to chromatographically and spectroscopically pure compounds"


def add_yield(r, product, percent, ev, oid="yield", conflict=None):
    if conflict:
        alternatives = [f(raw, norm, evidence) for raw, norm, evidence in conflict]
        vf = f([x[0] for x in conflict], None, ev, status="conflict", alternatives=alternatives)
        outcome(r, oid, "yield", [product], ev, None, role="desired product",
                yield_type="isolated", basis=isolated_basis(ev), value_fact=vf)
    else:
        outcome(r, oid, "yield", [product], ev, {"number": percent, "unit": "%"},
                role="desired product", yield_type="isolated", basis=isolated_basis(ev))


# Scheme 1: cited prior work and current generic strategy are not concrete current-paper experiments.
r = reaction("scheme1_previous_hydrogenation", [P[1]], record_type="literature_example",
             location="Scheme 1, Previous work", name="Hydrogenation of alpha-amino allenylphosphonates to vinylphosphonates",
             transformation="literature-reported allene partial reduction to vinylphosphonate")
add_participant(r, "p1", "generic_1", "reactant", [P[1]])
add_participant(r, "p2", "generic_3", "product", [P[1]], product_role="Z-major vinylphosphonate")
add_participant(r, "p_h2", "a_h2", "reagent", [P[1]])
add_participant(r, "p_cat", "a_pdc", "catalyst", [P[1]])
add_participant(r, "p_etoh", "a_ethanol", "solvent", [P[1]])
add_step(r, "s1", "reaction", ["p1", "p_h2", "p_cat", "p_etoh"], [P[1]], "H2, Pd/C, quinoline in EtOH at room temperature.", temperature="room temperature", atmosphere="hydrogen")
classification(r, "protected")
bundle["reactions"].append(r)

r = reaction("scheme1_current_generic", [P[1]], record_type="generic_scheme",
             location="Scheme 1, This work", name="Electrochemical reduction of alpha-amino allenylphosphonates to allylphosphonates")
add_participant(r, "p1", "generic_1", "reactant", [P[1]])
add_participant(r, "p2", "generic_2", "product", [P[1]], product_role="allylphosphonate")
add_step(r, "s1", "reaction", ["p1", "p2"], [P[1]], "Controlled-potential reduction in a divided cell.")
r["special_conditions"]["reactor"] = f("divided cell", "divided electrochemical cell", [P[1]])
classification(r, "protected")
bundle["reactions"].append(r)


# Page-2 prose-only supporting-electrolyte and solvent trials.
for suffix, agent_mid, label in (("nbu4nclo4", "a_nbu4nclo4", "nBu4NClO4"), ("liclo4", "a_liclo4", "LiClO4")):
    r = reaction(f"optimization_supporting_electrolyte_{suffix}", [P[2]], categories=["optimization"],
                 location="Page 2, supporting-electrolyte prose", unresolved=[
                     "Exact potential, scale, conversion, products, time, and result metric are not reported for this electrolyte trial.",
                     "The source only calls the results less conclusive."])
    add_participant(r, "p_substrate", "m4", "reactant", [P[2]])
    add_participant(r, "p_electrolyte", agent_mid, "additive", [P[2]])
    add_step(r, "s1", "reaction", ["p_substrate", "p_electrolyte"], [P[2]], f"Investigated {label} as supporting electrolyte.")
    outcome(r, "o1", "other", [], [P[2]], "less conclusive results", scope="qualitative optimization observation")
    classification(r, "unresolved")
    bundle["reactions"].append(r)

r = reaction("optimization_dmf_no_external_proton", [P[2]], categories=["optimization"],
             location="Page 2, solvent discussion", unresolved=[
                 "Exact conversion, yield, reaction time, and product-specific analytical basis are unreported."])
set_scheme_conditions(r, "m4", ["m5"], [P[2]], solvent="dmf", malonate=False, detailed=False)
outcome(r, "o_conversion", "conversion", ["m4"], [P[2]], "incomplete", scope="qualitative starting-material conversion")
outcome(r, "o_time", "other", ["m5"], [P[2]], "much longer than in acetonitrile", scope="qualitative reaction-time comparison")
classification(r, "detosyl")
bundle["reactions"].append(r)


def table_record(entry, proton_mid, conversion, yield_value, malonate=False, mixed=False):
    rid = f"table1_entry{entry}" + ("_mixed_solvent" if mixed else "")
    ev = [P[2]] + ([P[3], P[5]] if mixed else ([P[3]] if entry == 1 else []))
    r = reaction(rid, ev, categories=["optimization"], location=f"Table 1 entry {entry}" + (" footnote a" if mixed else ""),
                 unresolved=["Yield measurement basis is unreported." ] if not mixed else [
                     "Table reports 40-120 min and 4 F/mol, while the product-5 General Procedure A reports usually 5 h and 4-5 F/mol; both are retained.",
                     "Candidate structures and partial bond changes require expert review."])
    solvent = "mixed" if mixed else "mecn"
    set_scheme_conditions(r, "m4", ["m5"], ev, solvent=solvent, malonate=malonate,
                          detailed=mixed, malonate_equiv=4.5)
    if proton_mid and not malonate:
        add_participant(r, "p_proton_source", proton_mid, "reagent", [P[2]])
        r["steps"][0]["participant_ids"].append("p_proton_source")
    if not mixed:
        r["steps"][0]["time"] = f("40 to 120 min", {"minimum": 40, "maximum": 120, "unit": "min", "raw_text": "40 to 120 min"}, [P[2]])
    else:
        a = f("40 to 120 min", {"minimum": 40, "maximum": 120, "unit": "min", "raw_text": "40 to 120 min"}, [P[2]])
        b = f("usually 5 h", {"number": 5, "unit": "h"}, [P[5]])
        r["steps"][2]["time"] = f(["40 to 120 min", "usually 5 h"], None, [P[2], P[5]], status="conflict", alternatives=[a, b])
        q1 = f("4 F/mol", {"number": 4, "unit": "F/mol"}, [P[2], P[3]])
        q2 = f("usually 4-5 F/mol", {"minimum": 4, "maximum": 5, "unit": "F/mol"}, [P[5]])
        r["special_conditions"]["current"] = f(["4 F/mol", "usually 4-5 F/mol"], None, [P[2], P[3], P[5]], status="conflict", alternatives=[q1, q2])
    outcome(r, "o_conversion", "conversion", ["m4"], [P[2]], {"number": conversion, "unit": "%"}, scope="Table 1 conversion")
    if yield_value == "ND":
        outcome(r, "o_yield_not_determined", "other", ["m5"], [P[2]], "ND (not determined)", scope="Table 1 yield cell")
    elif yield_value is not None:
        basis = "isolated yield" if mixed else None
        outcome(r, "o_yield", "yield", ["m5"], [P[2]] + ([P[5]] if mixed else []),
                {"number": yield_value, "unit": "%"}, role="desired product",
                yield_type="isolated" if mixed else None, basis=basis,
                scope="Table 1 product-5 yield")
    if mixed:
        outcome(r, "o_product_amount", "product_amount", ["m5"], [P[5]], {"number": 54, "unit": "mg"},
                role="isolated product", basis="after silica-gel chromatography")
    classification(r, "detosyl")
    return r


bundle["reactions"] += [
    table_record(1, None, 100, 66, malonate=False),
    table_record(2, "a_succinimide", 75, "ND"),
    table_record(3, "a_phenol", 100, 66),
    table_record(4, "a_benzoic", 9, "ND"),
    table_record(5, "a_malonate", 100, 58, malonate=True),
    table_record(5, "a_malonate", 100, 68, malonate=True, mixed=True),
]

# Scheme 2 is an author-proposed mechanism, not an independently executed mechanistic experiment.
r = reaction("scheme2_proposed_mechanism", [P[2]], record_type="proposed_mechanism",
             categories=["mechanistic_study"], location="Scheme 2", name="Proposed mechanism for reduction of 4 to 5",
             transformation="author-proposed electron/proton sequence with beta-fragmentation",
             unresolved=["No independent isotope, trapping, kinetic, or control experiment is shown on pages 1-6.",
                         "Intermediates 6-11 and Ts radical are proposed, not isolated products."])
add_participant(r, "p4", "m4", "reactant", [P[2]])
add_participant(r, "p5", "m5", "product", [P[2]], product_role="overall proposed product")
add_participant(r, "p_malonate", "a_malonate", "reagent", [P[2]])
add_step(r, "s1", "reaction", ["p4", "p_malonate", "p5"], [P[2]], "Proposed four-electron/protonation sequence with N-Ts beta-fragmentation.")
classification(r, "detosyl")
bundle["reactions"].append(r)

# Scheme 3 generic protecting-group arrow.
r = reaction("scheme3_protecting_group_generic", [P[3]], record_type="generic_scheme",
             location="Scheme 3 generic arrow", name="Effect of the nitrogen protecting group",
             unresolved=["The generic arrow spans cleavage and retention outcomes depending on EWG; a single complete classification is unsafe."])
add_participant(r, "p_substrate", "generic_1", "reactant", [P[3]])
add_participant(r, "p_product", "generic_2", "product", [P[3]], product_role="products 12-15")
add_participant(r, "p_electrolyte", "a_et4nbf4", "additive", [P[3]])
add_participant(r, "p_malonate", "a_malonate", "reagent", [P[3]])
add_participant(r, "p_mecn", "a_mecn", "solvent", [P[3]])
add_participant(r, "p_dmf", "a_dmf", "solvent", [P[3]])
add_step(r, "s1", "reaction", ["p_substrate", "p_product", "p_electrolyte", "p_malonate"], [P[3]],
         "Controlled-potential electrolysis: 2-4 F/mol, room temperature, -2.50 V; footnotes modify potential/solvent.", temperature="room temperature")
classification(r, "unresolved")
bundle["reactions"].append(r)


product_mass = {"m12": 44, "m13": 6, "m14": 77, "m15": 118, "m16": 58, "m17": 55,
                "m18": 24, "m19": 20, "m22": 65, "m23": 48, "m24": 43,
                "m25": 12, "m25prime": 48, "m26": 49, "m27": 50}


def scope_record(rid, precursor, products_, ev, yield_rows=(), conversion_rows=(), selectivity_rows=(),
                 potential=-2.50, fmol=4, solvent="mixed", malonate=True, malonate_equiv=4.5,
                 mode="detosyl", detailed=False, unresolved=None, mass_rows=()):
    r = reaction(rid, ev, categories=["substrate_scope"], location=rid.replace("_", " "),
                 unresolved=unresolved)
    set_scheme_conditions(r, precursor, products_, ev, potential=potential, fmol=fmol,
                          solvent=solvent, malonate=malonate, malonate_equiv=malonate_equiv,
                          detailed=detailed,
                          product_roles={m: ("major product" if i == 0 else "byproduct") for i, m in enumerate(products_)})
    for item in yield_rows:
        mid, percent, yev = item[:3]
        if len(item) == 4 and item[3] == "conflict13":
            add_yield(r, mid, None, yev, conflict=[
                ("56% in Scheme 3/prose", {"number": 56, "unit": "%"}, [P[2], P[3]]),
                ("55% in product paragraph", {"number": 55, "unit": "%"}, [P[5]]),
            ])
        else:
            add_yield(r, mid, percent, yev)
    for mid, value_, cev in conversion_rows:
        outcome(r, f"conversion_{mid}", "conversion", [precursor], cev,
                {"number": value_, "unit": "%"}, scope=f"conversion associated with product {mid[1:]}")
    for mid, ratio, sev, metric in selectivity_rows:
        outcome(r, f"selectivity_{mid}", "selectivity", [mid], sev, ratio,
                role="stereoisomer mixture", basis="crude 31P NMR" if P[3] in sev else None,
                scope={"metric": metric, "order": metric.replace(" ratio", "")})
    for mid, mg, mev in mass_rows:
        outcome(r, f"amount_{mid}", "product_amount", [mid], mev, {"number": mg, "unit": "mg"},
                role="isolated product", basis="after flash chromatography")
    classification(r, mode)
    return r


# Scheme 3 concrete executions.
bundle["reactions"].append(scope_record("scheme3_product12", "pre12", ["m12"], [P[2], P[3], P[5]],
    yield_rows=[("m12", 55, [P[2], P[3], P[5]])], detailed=True, mode="detosyl",
    mass_rows=[("m12", 44, [P[5]])]))
bundle["reactions"].append(scope_record("scheme3_product13", "pre13", ["m13"], [P[2], P[3], P[5]],
    yield_rows=[("m13", None, [P[2], P[3], P[5]], "conflict13")], detailed=True, fmol=2,
    malonate_equiv=2, mode="protected", mass_rows=[("m13", 6, [P[5]])],
    unresolved=["Scheme/prose report 56% while the product paragraph reports 55%; no value chosen.",
                "The printed 6 mg product amount appears internally difficult to reconcile with 0.215 mmol and 55-56%; it is preserved without correction.",
                "Candidate structures and partial bond changes require expert review."]))

for potential in (-2.70, -2.90):
    r = scope_record(f"scheme3_ms_potential_{str(abs(potential)).replace('.', '_')}V", "pre13", ["m13", "m5"], [P[2], P[3]],
                     potential=potential, fmol=2, malonate_equiv=2, mode="detosyl",
                     unresolved=["No individual yields, scale, time, or ratio measurement method are reported.",
                                 "The 5:13 ratio is preserved in the author's stated order."])
    outcome(r, "selectivity_5_vs_13", "selectivity", ["m5", "m13"], [P[2]], "1:2",
            role="minor 5 / major 13", scope={"metric": "product ratio", "order": ["5", "13"]})
    bundle["reactions"].append(r)

bundle["reactions"].append(scope_record("scheme3_product14", "pre14", ["m14"], [P[3], P[5]],
    yield_rows=[("m14", 77, [P[3], P[5]])], potential=-2.70, detailed=True, mode="protected",
    mass_rows=[("m14", 77, [P[5]])]))
bundle["reactions"].append(scope_record("scheme3_product14_pure_mecn", "pre14", ["m14"], [P[2], P[3]],
    yield_rows=[("m14", 55, [P[2], P[3]])], potential=-2.70, solvent="mecn", malonate=False,
    mode="protected", unresolved=["Scale, time, work-up, purification, and product amount are unreported for this pure-MeCN execution."]))
bundle["reactions"].append(scope_record("scheme3_product15", "pre15", ["m15"], [P[3], P[5]],
    yield_rows=[("m15", 83, [P[3], P[5]])], potential=-2.70, detailed=True, mode="protected",
    mass_rows=[("m15", 118, [P[5]])]))


# Scheme 4 concrete scope; product 5 is the Table-1 mixed-solvent record above.
scope_specs = [
    ("scope16", "pre16", "m16", 72, None, None, [P[3], P[5]], True),
    ("scope16_pure_mecn", "pre16", "m16", 35, None, None, [P[2], P[3]], False),
    ("scope17", "pre17", "m17", 77, None, None, [P[3], P[5]], True),
    ("scope17_pure_mecn", "pre17", "m17", 37, None, None, [P[2], P[3]], False),
    ("scope18", "pre18", "m18", 21, None, ("1:1", "dr ratio"), [P[3], P[5]], True),
    ("scope19", "pre19", "m19", 20, None, ("1:1", "E/Z ratio"), [P[3], P[5], P[6]], True),
    ("scope22", "pre22", "m22", 77, None, None, [P[3], P[6]], True),
    ("scope23", "pre23", "m23", 56, 88, None, [P[3], P[6]], True),
    ("scope24", "pre24", "m24", 46, 81, ("1:1", "E/Z ratio"), [P[3], P[6]], True),
    ("scope26", "pre26", "m26", 56, None, ("1:1", "E/Z ratio"), [P[3], P[6]], True),
    ("scope27", "pre27", "m27", 70, None, ("1:1", "E/Z ratio"), [P[3], P[6]], True),
]
for rid, pre, prod, y, conv, sel, ev, detailed in scope_specs:
    pure = rid.endswith("pure_mecn")
    r = scope_record(rid, pre, [prod], ev,
                     yield_rows=[(prod, y, ev)],
                     conversion_rows=[(prod, conv, [P[3]])] if conv is not None else [],
                     selectivity_rows=[(prod, sel[0], [P[3]], sel[1])] if sel else [],
                     solvent="mecn" if pure else "mixed", malonate=not pure,
                     detailed=detailed, mode="detosyl",
                     mass_rows=[(prod, product_mass[prod], [P[5]] if prod in {"m16", "m17", "m18"} else [P[6]])] if detailed else [],
                     unresolved=["Scale, reaction time, work-up, purification, and product amount are unreported for this pure-MeCN execution."] if pure else None)
    bundle["reactions"].append(r)

for prod, pre in (("m20", "pre20"), ("m21", "pre21")):
    r = scope_record(f"scope{prod[1:]}_decomposition", pre, [prod], [P[3]], yield_rows=[], mode="unresolved",
                     unresolved=["The expected product decomposed under the reaction conditions; decomposition products and analytical method are unreported.",
                                 "No yield or stable isolated product is asserted."])
    outcome(r, f"decomposition_{prod}", "decomposition", [prod], [P[3]], "decomposed under the reaction conditions",
            role="unstable expected product", scope="qualitative decomposition observation")
    bundle["reactions"].append(r)

# Product 25/25' is one experiment with two product-specific isolated yields.
r = scope_record("scope25_competing_products", "m4a", ["m25", "m25prime"], [P[3], P[4], P[6]],
                 yield_rows=[("m25", 12, [P[3], P[4], P[6]]), ("m25prime", 72, [P[4], P[6]])],
                 selectivity_rows=[("m25", "9:1, E major", [P[3], P[4]], "E/Z ratio"),
                                   ("m25prime", "9:1, E major", [P[4]], "E/Z ratio")],
                 potential=-2.40, detailed=True, mode="detosyl",
                 mass_rows=[("m25", 12, [P[6]]), ("m25prime", 48, [P[6]])],
                 unresolved=["The mechanistic origin and complete bond changes for byproduct 25′ require expert review.",
                             "The product-ratio 25:25′ is not explicitly reported; only separate isolated yields are given."])
bundle["reactions"].append(r)

r = scope_record("scope25_potential_minus2_30V", "m4a", ["m25", "m25prime"], [P[4]],
                 potential=-2.30, detailed=False, mode="detosyl",
                 unresolved=["The source says the product ratio was not modified, but gives no yields or explicit 25:25′ numerical ratio for this potential.",
                             "Exact scale, time, work-up, and measurement method are unreported."])
outcome(r, "o_ratio", "other", ["m25", "m25prime"], [P[4]], "Changing to -2.30 V did not modify the ratio of products formed.",
        scope={"metric": "25:25′ product ratio", "order": ["25", "25′"], "numeric_value": None})
bundle["reactions"].append(r)


# Scheme 6 N-substituent scope.
scheme6_specs = [
    ("scope28", "pre28", "m28", 61, 81, -2.50, "mixed", True, [P[4], P[6]]),
    ("scope28_pure_mecn", "pre28", "m28", 52, 100, -2.50, "mecn", False, [P[2], P[4]]),
    ("scope29", "pre29", "m29", 47, 81, -2.50, "mixed", True, [P[4]]),
    ("scope30", "pre30", "m30", 61, None, -2.50, "mixed", True, [P[4]]),
    ("scope31", "pre31", "m31", 47, 91, -2.50, "mixed", True, [P[4]]),
    ("scope32", "pre32", "m32", 41, 72, -2.50, "mixed", True, [P[4]]),
    ("scope33", "pre33", "m33", 78, None, -2.42, "mixed", True, [P[4]]),
    ("scope33_pure_mecn", "pre33", "m33", 45, None, -2.42, "mecn", False, [P[2], P[4]]),
    ("scope34", "pre34", "m34", 53, 92, -2.50, "mixed", True, [P[4]]),
]
for rid, pre, prod, y, conv, pot, solvent, malonate, ev in scheme6_specs:
    r = scope_record(rid, pre, [prod], ev, yield_rows=[(prod, y, ev)],
                     conversion_rows=[(prod, conv, [P[4]])] if conv is not None else [],
                     potential=pot, solvent=solvent, malonate=malonate, detailed=False, mode="detosyl",
                     unresolved=["The covered pages give scheme conditions/outcomes but not product-specific scale, time, work-up, purification, or product amount.",
                                 "No E/Z or alpha-center configuration is assigned."])
    bundle["reactions"].append(r)


# Source context: analytical blocks and apparatus are not extra synthesis records.
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:context_title_abstract",
        "document_id": DOC, "evidence_ids": [P[1]], "kind": "bibliographic_and_abstract",
        "text": "Electrochemical Reduction of alpha-Amino Allenylphosphonates to alpha-Amino Allylphosphonates; authors Pauline Adler, Pedro de Oliveira, Nicolas Rabasso."
    },
    {
        "context_id": f"{DOC}:context_mechanistic_discussion",
        "document_id": DOC, "evidence_ids": [P[2]], "kind": "author_proposed_mechanism",
        "text": "Radical-anion/protonation/beta-fragmentation sequence proposed in Scheme 2; not treated as an independently verified mechanistic experiment."
    },
    {
        "context_id": f"{DOC}:context_electrochemical_apparatus",
        "document_id": DOC, "evidence_ids": [P[4], P[5]], "kind": "general_instrumentation",
        "text": "EG&G 273A potentiostat; SCE reference; Pt-gauze counter; glassy-carbon working electrode; argon deaeration; divided cell."
    },
    {
        "context_id": f"{DOC}:context_characterization_pages5_6",
        "document_id": DOC, "evidence_ids": [P[5], P[6]], "kind": "compound_characterization",
        "text": "NMR, HRMS, Rf, physical descriptions, and chromatography details for products 5, 12-19, and 22-28. Spectral peak lists are context, not separate synthetic reactions."
    },
]

# Raw model provenance artifacts.
stereo_checks = {}
for m in bundle["molecules"]:
    smi = st.value(m["smiles"])
    if smi:
        try:
            stereo_checks[m["molecule_id"]] = st.stereo_inspect(smi)
        except Exception as exc:
            stereo_checks[m["molecule_id"]] = {"error": str(exc)}
dump(RUN / "raw_model/stereo_checks.json", stereo_checks)

mentions = st.source_mentions(dataset, ["4", "5", "13", "14", "16", "17", "18", "19", "25", "25′", "28", "33"], context=180)
dump(RUN / "raw_model/source_mentions.json", mentions)

page_inventory = {
    "document_id": DOC,
    "visually_inspected_pages": [1, 2, 3, 4, 5, 6],
    "image_detail": "original",
    "pages": [
        {"page": 1, "regions": ["title/abstract", "Scheme 1 previous work arrow", "Scheme 1 this-work arrow", "introduction/background"]},
        {"page": 2, "regions": ["optimization prose", "Table 1 rows 1-5 and footnote a", "Scheme 2 mechanism", "protecting-group prose"]},
        {"page": 3, "regions": ["Scheme 3 arrow/products/footnotes", "Scheme 4 products 5,16-27", "Scheme 4 caption/footnote", "scope prose"]},
        {"page": 4, "regions": ["Scheme 5 arrow and products 25/25′", "Scheme 6 products 28-34/footnotes", "conclusion", "general electrochemical apparatus"]},
        {"page": 5, "regions": ["General Procedure A", "product paragraphs 5,12-19", "NMR/HRMS/Rf characterization"]},
        {"page": 6, "regions": ["continuation product 19", "product paragraphs 22-28", "product 25′", "NMR/HRMS/Rf characterization"]},
    ],
}
dump(RUN / "raw_model/page_visual_inventory.json", page_inventory)

# Coverage is derived from the visual inventory and explicit record evidence, not from count alone.
page_records = {i: [r["reaction_id"] for r in bundle["reactions"] if P[i] in r["evidence_ids"]] for i in range(1, 7)}
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title and abstract", "finding": "Article identity and overall current-work claim; no standalone experiment.", "evidence_ids": [P[1]]},
            {"region": "Scheme 1, Previous work arrow", "finding": "Cited hydrogenation literature example 1 to 3; kept as literature_example.", "evidence_ids": [P[1]]},
            {"region": "Scheme 1, This work arrow", "finding": "Generic electrochemical route 1 to 2; kept as generic_scheme, not a concrete execution.", "evidence_ids": [P[1]]},
            {"region": "introduction prose", "finding": "Background and citations; no additional current-paper concrete experiment.", "evidence_ids": [P[1]]},
        ], "record_ids": page_records[1], "unresolved": ["Generic R/EWG/R′ identities and all stereochemical assignments remain unresolved."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "supporting-electrolyte and solvent prose", "finding": "nBu4NClO4, LiClO4, and pure-DMF trials retained as partial optimization experiments.", "evidence_ids": [P[2]]},
            {"region": "Table 1 rows 1-5 and footnote", "finding": "Five rows plus the footnoted mixed-solvent execution represented independently; ND preserved as not determined.", "evidence_ids": [P[2]]},
            {"region": "Scheme 2", "finding": "Author-proposed mechanism 4 to 5, not independent proof.", "evidence_ids": [P[2]]},
            {"region": "protecting-group prose", "finding": "Products 12/13 and -2.70/-2.90 V 5:13 trials linked to Scheme 3 records.", "evidence_ids": [P[2]]},
        ], "record_ids": page_records[2], "unresolved": ["Table yields other than footnoted 68% lack an explicit measurement basis.", "Less-conclusive electrolyte trials lack numeric outcomes."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Scheme 3 arrow/products/footnotes", "finding": "Generic protecting-group scheme plus concrete executions 12-15, including alternate pure-MeCN 14.", "evidence_ids": [P[3]]},
            {"region": "Scheme 4, products 5 and 16-18", "finding": "Main and footnoted pure-MeCN executions separated; 18 dr retained.", "evidence_ids": [P[3]]},
            {"region": "Scheme 4, products 19-27", "finding": "Scope products, decompositions, conversions, and E/Z outcomes each mapped to records.", "evidence_ids": [P[3]]},
            {"region": "Scheme 4 caption", "finding": "Shared electrolyte/solvent/proton-source/charge/potential and crude 31P-NMR ratio basis retained.", "evidence_ids": [P[3]]},
        ], "record_ids": page_records[3], "unresolved": ["Most source drawings do not assign alpha-center or allene axial stereochemistry.", "Products 20 and 21 decompose; decomposition identities are not given."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Scheme 5", "finding": "One 4a execution with product-specific yields for 25 and 25′; separate -2.30 V execution retained from prose.", "evidence_ids": [P[4]]},
            {"region": "Scheme 6", "finding": "Products 28-34 and footnoted pure-MeCN variants mapped independently with conversions.", "evidence_ids": [P[4]]},
            {"region": "conclusion", "finding": "Summary only; no new reaction record.", "evidence_ids": [P[4]]},
            {"region": "Experimental Section apparatus", "finding": "General instrumentation/context retained; cyclic voltammograms are not counted as synthesis records.", "evidence_ids": [P[4]]},
        ], "record_ids": page_records[4], "unresolved": ["The -2.30 V 25:25′ product ratio is described only as unchanged, without a number.", "Covered pages omit product-specific procedures for 29-34."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "General Procedure A", "finding": "Scale, doses, solvent volume/ratio, degassing, electrolysis endpoint, work-up, and chromatography mapped in operation order.", "evidence_ids": [P[5]]},
            {"region": "product paragraphs 5 and 12-19", "finding": "Product masses, yields, purification details, and mixture statements associated with existing reactions; not duplicated.", "evidence_ids": [P[5]]},
            {"region": "NMR/HRMS/Rf blocks", "finding": "Characterization retained as source_context rather than reactions.", "evidence_ids": [P[5]]},
        ], "record_ids": page_records[5], "unresolved": ["Product 13 prints 6 mg and 55%, while Scheme 3/prose prints 56%; values are preserved without correction.", "General Procedure A 5 h/4-5 F mol-1 differs from Table 1 product-5 40-120 min/4 F mol-1."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "product 19 continuation", "finding": "E/Z-mixture characterization continues from page 5; no new synthesis record.", "evidence_ids": [P[6]]},
            {"region": "product paragraphs 22-27 and 25′", "finding": "Masses/yields/purification and mixture characterization linked to Scheme 4/5 records.", "evidence_ids": [P[6]]},
            {"region": "product 28 heading/start", "finding": "Preparation paragraph begins and continues beyond bounded page 6; page-4 scheme supplies covered outcome.", "evidence_ids": [P[6]]},
            {"region": "NMR/HRMS/Rf blocks", "finding": "Characterization preserved in source_context, not counted as reactions.", "evidence_ids": [P[6]]},
        ], "record_ids": page_records[6], "unresolved": ["Product 28 procedure/characterization continues onto page 7 and is outside this visual-coverage claim.", "No SI file was supplied or associated."]},
    ],
}
for page_entry in coverage["pages"]:
    for region in page_entry["regions"]:
        region["evidence_ids"] = with_page_images(region["evidence_ids"])
dump(RUN / "coverage/pages_0001_0006.json", coverage)

dump(RUN / "submitted.json", bundle)

workflow = f"""# Extraction workflow log

## Batch, source, and bounded coverage

- Batch: `allene-stereo-v31-200-20261017`.
- Document: `{DOC}`, *Electrochemical Reduction of α-Amino Allenylphosphonates to α-Amino Allylphosphonates*, DOI `10.1002/ejoc.202000534`.
- Current visual coverage is exactly original PDF pages 1-6. All six supplied page images were opened at original detail. Other extracted pages were used only for literal mention searching; no visual or semantic coverage beyond page 6 is claimed.
- `runs/current/source_mapping.json` identifies only the main-text PDF. The article advertises SI, but no SI file was explicitly supplied or associated.
- The prior `runs/workflow.md` concerned a different article and was replaced completely. No old molecule, reaction, decision, count, or review state was reused.

## Workflow and source decisions

1. Read `AGENTS.md`, the complete extraction skill, `extraction_standard.md`, `stereochemistry.md`, `structured_fields.md`, `quality_rules.md`, `literature.schema.json`, and `bond_taxonomy.json`. Used only `/EvoLab/student-workspaces/xllin/OpenEvo-chem-integration/.venv/bin/python` for bundled tools. No schema, tool, test, quality rule, or control program was changed.
2. Confirmed the prepared dataset had zero candidates, zero experiments, and no v3 bundle. Added the exact title, ran `standard init`, and generated the source-derived submission from the empty template.
3. Page 1 Scheme 1 has two arrows: the hydrogen/Pd route is prior work and is `literature_example`; the electrochemical route is a `generic_scheme`. Neither is counted as a current concrete execution.
4. Page 2 Table 1 has five printed rows. The parenthesized 68% footnote changes the solvent to CH3CN/DMF and is a distinct execution, so the table contributes six records. `ND` is retained as not determined, not zero/no reaction. Only the 68% footnote explicitly says isolated; other table-yield bases remain missing.
5. The prose-only nBu4NClO4, LiClO4, and pure-DMF trials are separate partial optimization records because they are different executed conditions. Their missing numeric outcomes are not invented.
6. Scheme 2 is an author-proposed radical-anion/protonation/beta-fragmentation mechanism. It is a `proposed_mechanism`, not independent mechanistic proof; no isotope, trapping, kinetic, or control experiment appears on pages 1-6.
7. Scheme 3 has a generic protecting-group arrow and concrete executions for 12-15. The -2.70 and -2.90 V mesyl trials are independent records with the stated 5:13 = 1:2 order. Product 14's pure-MeCN/no-malonate 55% execution is separate from its 77% mixed-solvent execution.
8. Scheme 4 was expanded by actual product/condition. Scheme footnote pure-MeCN executions for 5, 16, and 17 were not merged with mixed-solvent executions. Product 5's mixed-solvent 68% and pure-MeCN 66% outcomes are the same experiments already represented by Table 1 and therefore only gain cross-page evidence, not duplicate reaction IDs.
9. Products 20 and 21 are decomposition outcomes. Their expected/transient identities are retained, but no stable isolated product, decomposition product, or yield is invented.
10. Scheme 5 is one -2.40 V execution producing both 25 (12%) and 25′ (72%); both yields and both 9:1 E/Z outcomes are product-targeted. The -2.30 V prose trial is a separate execution, but its exact 25:25′ ratio and yields are unreported.
11. Scheme 6 has seven main scope executions plus the two footnoted pure-MeCN variants for 28 and 33. Conversion is distinct from yield. Covered pages do not provide product-specific scale/time/work-up for 29-34, so those fields remain missing.
12. General Procedure A was decomposed into setup, 15 min argon degassing, electrolysis, solvent removal, EtOAc wash/Celite filtration, evaporation, and silica chromatography. Explicit substrate/malonate doses and 25 mL electrolyte solution were structured. Unknown mobile-phase volumes were not inferred.
13. Table 1 reports 40-120 min and 4 F/mol for product 5, while General Procedure A reports usually 5 h and 4-5 F/mol. The mixed-solvent product-5 record retains these as conflicts with both page sources; neither value is selected.
14. Scheme 3/prose reports product 13 at 56%, while the product paragraph reports 55% and 6 mg. The yield is a conflict fact and the mass is preserved verbatim; the source was not arithmetically corrected.
15. NMR, HRMS, Rf, physical descriptions, and spectral peak lists on pages 5-6 are `source_context` linked to preparations, not extra synthesis records. Cyclic-voltammetry apparatus is also context, not a reaction.

## Molecules, stereochemistry, and classification

- Molecules are independently registered and reactions refer only by `molecule_id`. Generic R/EWG/R′ structures are templates, not wildcard SMILES. CAS values are absent and were not guessed.
- Candidate concrete SMILES were transcribed for drawn/named products and precursors. Every candidate was passed through the bundled `standard stereo` implementation; diagnostics are in `runs/current/raw_model/stereo_checks.json`. Parsing is not claimed as source verification.
- Structure-derived facts cite both the supplied page text and page image for pages 1-6. The images, not OCR alone, support scheme/table structures, arrows, footnotes, and row values.
- No product alpha-center R/S assignment, allene R_a/S_a assignment, or racemic label was invented. The allene terminal carbon is not misrepresented as a tetrahedral center: unsymmetrically substituted allene axes remain unassigned/partial, while candidates with duplicate substituents at one allene terminus are marked axis-not-applicable by structure derivation. Compound 18 separately records its unassigned alpha-carbon and phosphorus stereogenic sites plus the reported 1:1 diastereomer mixture. E/Z mixtures remain mixture identities with experiment-specific selectivity outcomes.
- Bond changes are explicit but partial/inferred: two net C-H gains and the allene C=C bond-order decrease; detosylative paths additionally retain N-H formation and N-S cleavage. The C=C order decrease is not mislabeled as a new C-C sigma bond. No proton provenance, complete atom mapping, or three-model vote is claimed.
- All machine records remain `needs_review`; semantic and chemical accuracy are false.

## Raw artifacts and coverage

- `runs/current/raw_model/submitted_init.json`: exact `standard init` template.
- `runs/current/raw_model/build_records.py`: this transparent source-build script.
- `runs/current/raw_model/page_visual_inventory.json`: page-region inventory for the six viewed images.
- `runs/current/raw_model/stereo_checks.json`: graph/stereo diagnostics for candidate SMILES.
- `runs/current/raw_model/source_mentions.json`: literal cross-region hints only; it does not merge experiments.
- `runs/current/coverage/pages_0001_0006.json`: exact pages 1-6, every page attempted and every page with nonempty image-locatable regions.

## Unresolved limitations to preserve through validation/audit

- Final audit has no `reaction_identity`, `outcome_target`, or `structure_gap` warnings. Its remaining warnings are `measurement_basis` gaps for conversions and the non-footnoted Table 1 yields; the source gives no analytical method/basis for those values, so the facts remain missing and are not auto-filled.
- No associated SI; product-specific details that would normally be in SI remain missing.
- Table yield bases except the explicit isolated 68% footnote; conversion analytical bases; exact outcomes for alternative electrolytes; pure-DMF numerical conversion/time.
- Product-5 timing/charge conflict; product-13 56% versus 55% conflict and printed 6 mg amount.
- Exact 25:25′ product ratio at -2.30 V; decomposition identities for 20/21; product-specific scale/time/work-up for 29-34.
- Full source-verified stereochemical structures, allene axial configurations, separated E/Z identities, complete atom mapping, and proton sources.
- Schema validity and graph parsing do not establish chemistry correctness, semantic completeness, or reproducibility.
"""
(ROOT / "runs/workflow.md").write_text(workflow)

tasks_path = RUN / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    page = int(task["evidence_id"].split(":page", 1)[1].split(":", 1)[0])
    if 1 <= page <= 6:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected at original detail; findings and unresolved items are recorded in coverage/pages_0001_0006.json."
dump(tasks_path, tasks)

# On a post-merge rerun, persist the advisory audit alongside raw/provenance
# artifacts.  Warnings are not auto-repaired and never certify correctness.
if dataset.get("experiments"):
    import reaction_tools
    dump(RUN / "raw_model/audit.json", reaction_tools.dataset_audit(dataset))

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "page_record_counts": {str(k): len(v) for k, v in page_records.items()},
}, ensure_ascii=False, indent=2))
