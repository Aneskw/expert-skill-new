#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for the six supplied main-text pages.

This script deliberately leaves compact scope structures unresolved.  It only uses
the prepared evidence already present in runs/current/dataset.json.
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
DATASET = RUN / "dataset.json"
SUBMITTED = RUN / "submitted.json"
DOC = "doc_d044c13bc5ccf0c9"
TITLE = "Direct synthesis of pentasubstituted pyrroles and hexasubstituted pyrrolines from propargyl sulfonylamides and allenamides"


def ev(page: int, image: bool = True):
    ids = [f"{DOC}:page{page}:text"]
    if image:
        ids.append(f"{DOC}:page{page}:image")
    return ids


def F(raw=None, norm=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return st.fact(raw, norm, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


dataset = st.read(DATASET)
dataset["documents"][0]["title"] = TITLE
st.write(DATASET, dataset)
bundle = st.read(SUBMITTED)
bundle["documents"] = dataset["documents"]
molecules = {}


def molecule(mid, label, page, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, note=None,
             axial=None, assignment=None):
    evidence = ev(page)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name or label, name or label, evidence)
    m["names"]["english"] = F(name, name, evidence) if name else F(evidence=evidence)
    m["names"]["abbreviations"] = F(evidence=evidence)
    m["names"]["aliases"] = F(evidence=evidence)
    m["cas_numbers"] = F(evidence=evidence)
    m["generic_template"] = None
    if smiles:
        diag = st.stereo_inspect(smiles)
        canonical = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = F(None, canonical, evidence, method="structure_derivation",
                         derivation="Candidate connectivity transcribed from the supplied main-text drawing/name and canonicalized with RDKit; parsing does not certify source identity.")
        m["stereochemistry"]["isomeric_smiles"] = F(None, canonical, evidence, method="structure_derivation",
                         derivation="RDKit canonical isomeric SMILES for the candidate graph; unassigned potential stereo remains unassigned.")
        m["stereochemistry"]["identity_key"] = F(None, diag["stereo_graph_key"], evidence, method="structure_derivation",
                         derivation="Hash of the candidate canonical isomeric graph from standard stereo diagnostics.")
        rep = representation or ("partial" if diag["potential_stereo"] else "not_applicable")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(evidence=evidence)
        m["stereochemistry"]["isomeric_smiles"] = F(evidence=evidence)
        m["stereochemistry"]["identity_key"] = F(evidence=evidence)
        rep = representation or "not_encoded"
    for key in ("atom_centers", "double_bonds", "mixture"):
        m["stereochemistry"][key] = F(evidence=evidence)
    if axial is None:
        m["stereochemistry"]["axial"] = F(evidence=evidence)
    else:
        m["stereochemistry"]["axial"] = F(None, axial, evidence, method="structure_derivation",
            derivation="The drawn cumulene/allene has different substituents at both termini; no R_a/S_a or P/M assignment is printed.")
    m["stereochemistry"]["assignment"] = (
        F(assignment, assignment, evidence) if assignment else F(evidence=evidence)
    )
    status_derivation = "Status assigned after inspection of the supplied drawing; absence of a descriptor is not treated as achiral or racemic."
    m["stereochemistry"]["chirality_status"] = F(None, chirality, evidence, method="inference", derivation=status_derivation)
    m["stereochemistry"]["geometry_status"] = F(None, geometry, evidence, method="inference", derivation=status_derivation)
    m["stereochemistry"]["representation_status"] = F(None, rep, evidence, method="inference",
        derivation="Reports whether the candidate representation retains the source-visible stereochemical state; compact or unassigned drawings remain partial/not encoded.")
    m["identity_note"] = F(None, note, evidence, method="inference", derivation=note) if note else F(evidence=evidence)
    m["external_properties"] = []
    molecules[mid] = m
    return m


def add_occurrence(mid, label, page):
    evidence = ev(page)
    m = molecules[mid]
    if not any(o["compound_label"]["normalized_value"] == label and set(o["evidence_ids"]) == set(evidence) for o in m["occurrences"]):
        m["occurrences"].append({"document_id": DOC, "compound_label": F(label, label, evidence), "evidence_ids": evidence})


def generic_molecule(mid, label, page, fixed, variables, note):
    evidence = ev(page)
    m = molecule(mid, label, page, name=label, chirality="unresolved", geometry="unresolved",
                 representation="not_encoded", note=note)
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": F(f"Figure attachment point for {v}", {"label": v, "attachment": "as drawn"}, evidence),
            "definition": F(f"{v} as defined by the source scope drawing", f"source-defined {v}", evidence),
            "allowed_values": F(evidence=evidence),
            "constraints": F(evidence=evidence),
        } for v in variables],
        "instance_mapping": [],
        "reason": F(note, note, evidence),
    }
    return m


# Core structures and generic templates.
molecule("m_1a", "1a", 2, "N-(1,3-diphenylprop-2-yn-1-yl)-N-ethylbenzenesulfonamide",
         "CCN(S(=O)(=O)c1ccccc1)C(c1ccccc1)C#Cc1ccccc1",
         chirality="unspecified", geometry="not_applicable", representation="partial",
         note="Candidate name-derived graph; the propargylic stereocenter is unassigned in the source drawing.")
molecule("m_2a", "2a", 2, "2-amino-1-ethyl-3-cyano-4,5-diphenylpyrrole",
         "CCn1c(N)c(C#N)c(-c2ccccc2)c1-c1ccccc1",
         chirality="achiral", geometry="not_applicable", representation="not_applicable",
         note="Candidate graph derived from the clearly drawn Table 1 product; no stereogenic element is present in this candidate graph.")
molecule("m_3a", "3a", 3, "trisubstituted allenamide 3a",
         "CCN(S(=O)(=O)c1ccccc1)C(=C=C(c1ccccc1))c1ccccc1",
         chirality="unspecified", geometry="unspecified", representation="partial",
         axial={"axis": "allene C=C=C", "configuration": "unspecified"},
         note="Candidate connectivity from Fig. 3; axial configuration is not printed and is not encoded by the toolkit candidate.")
generic_molecule("m_generic_propargyl_1", "generic propargyl sulfonamide 1", 2,
                 {"motif": "R1-C#C-C(R2)(R3)-N(SO2R)-", "source": "Fig. 2 generic substrate 1"},
                 ["R", "R1", "R2", "R3"], "R-group source drawing; not a concrete SMILES.")
generic_molecule("m_generic_pyrrole_2", "generic pentasubstituted pyrrole 2", 2,
                 {"motif": "2-amino-3-cyanopyrrole bearing R1/R2/R3/N-substituent", "source": "Fig. 2 generic product 2"},
                 ["R", "R1", "R2", "R3"], "R-group source drawing; not a concrete SMILES.")
generic_molecule("m_generic_allene_4", "generic tetrasubstituted allenamide 4", 3,
                 {"motif": "R1R2C=C=C(R3)-N(R4)SO2Ph", "source": "Fig. 4 generic substrate 4"},
                 ["R1", "R2", "R3", "R4"], "Tetrasubstituted allene template; axial configuration is not assigned.")
generic_molecule("m_generic_pyrroline_5", "generic hexasubstituted pyrroline 5", 3,
                 {"motif": "amino/cyano pyrroline bearing R1-R4", "source": "Fig. 4 generic product 5"},
                 ["R1", "R2", "R3", "R4"], "R-group source drawing; product stereochemistry is not assigned.")


# Common reagents/solvents.  Simple candidates are structural aids, not source-certified identities.
simple = [
    ("m_tmscn", "TMSCN", "trimethylsilyl cyanide", "C[Si](C)(C)C#N"),
    ("m_cs2co3", "Cs2CO3", "cesium carbonate", "[Cs+].[Cs+].[O-]C([O-])=O"),
    ("m_k2co3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    ("m_k3po4", "K3PO4", "potassium phosphate", "[K+].[K+].[K+].[O-]P([O-])([O-])=O"),
    ("m_kf", "KF", "potassium fluoride", "[K+].[F-]"),
    ("m_koh", "KOH", "potassium hydroxide", "[K+].[OH-]"),
    ("m_kotbu", "KOtBu", "potassium tert-butoxide", "CC(C)(C)[O-].[K+]"),
    ("m_et3n", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("m_dmf", "DMF", "N,N-dimethylformamide", "CN(C)C=O"),
    ("m_acetonitrile", "CH3CN", "acetonitrile", "CC#N"),
    ("m_dme", "DME", "1,2-dimethoxyethane", "COCCOC"),
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("m_dce", "DCE", "1,2-dichloroethane", "ClCCCl"),
    ("m_dioxane", "dioxane", "1,4-dioxane", "O1CCOCC1"),
    ("m_dmso_d6", "DMSO-d6", "dimethyl sulfoxide-d6", "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]"),
    ("m_temppo", "TEMPO", "2,2,6,6-tetramethylpiperidin-1-oxyl", "CC1(C)CCCC(C)(C)N1[O]"),
    ("m_bzcl", "BzCl", "benzoyl chloride", "O=C(Cl)c1ccccc1"),
    ("m_pyridine", "pyridine", "pyridine", "n1ccccc1"),
    ("m_cubr", "CuBr", "copper(I) bromide", "[Cu]Br"),
    ("m_4f_benzaldehyde", "4-fluorobenzaldehyde", "4-fluorobenzaldehyde", "O=Cc1ccc(F)cc1"),
    ("m_phso2cl", "PhSO2Cl", "benzenesulfonyl chloride", "O=S(=O)(Cl)c1ccccc1"),
    ("m_dcm", "DCM", "dichloromethane", "ClCCl"),
    ("m_tfa", "TFA", "trifluoroacetic acid", "O=C(O)C(F)(F)F"),
    ("m_naoh", "NaOH", "sodium hydroxide", "[Na+].[OH-]"),
    ("m_methanol", "CH3OH", "methanol", "CO"),
    ("m_phenylacetylene", "phenylacetylene", "phenylacetylene", "C#Cc1ccccc1"),
]
for mid, label, name, smi in simple:
    molecule(mid, label, 2 if mid in {"m_tmscn","m_cs2co3","m_k2co3","m_k3po4","m_kf","m_koh","m_kotbu","m_et3n","m_dmf","m_acetonitrile","m_dme","m_toluene","m_dce","m_dioxane"} else 3,
             name, smi, chirality="achiral", geometry="not_applicable", representation="not_applicable")

for mid, label, name in [
    ("m_ni", "Ni(PPh3)2Cl2", "bis(triphenylphosphine)nickel(II) dichloride"),
    ("m_pd", "Pd(OAc)2", "palladium(II) acetate"),
    ("m_cu", "Cu(OAc)2", "copper(II) acetate"),
    ("m_co", "Co(OAc)2", "cobalt(II) acetate"),
]:
    molecule(mid, label, 2, name, chirality="unresolved", geometry="not_applicable", representation="not_encoded",
             note="Metal complex formula is reported; a complete coordination-state SMILES was not asserted.")


reactions = []


def reaction(suffix, page, record_type, categories, location, reaction_name, transformation, evidence=None):
    evidence = evidence or ev(page)
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(str(categories), categories, evidence)
    r["location"] = F(location, {"page": page, "location": location}, evidence)
    r["reaction_name"] = F(reaction_name, reaction_name, evidence)
    r["transformation_type"] = F(None, transformation, evidence, method="inference",
                                  derivation="Transformation class inferred from the source-drawn input and product labels.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": ["No explicitly associated SI was supplied; workup, purification and omitted procedural details remain unreported."]}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def participant(r, pid, mid, role, evidence, eq=None, amount=None, volume=None, mol_percent=None,
                reference=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    if eq is not None:
        p["charges"]["equivalents"] = F(f"{eq} equiv.", {"number": eq, "unit": "equiv"}, evidence)
    if amount is not None:
        p["charges"]["amount"] = F(f"{amount} mmol", {"number": amount, "unit": "mmol"}, evidence)
    if volume is not None:
        p["charges"]["volume"] = F(f"{volume} mL", {"number": volume, "unit": "mL"}, evidence)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = F(f"{mol_percent} mol%", {"number": mol_percent, "unit": "mol%"}, evidence)
    if (eq is not None or mol_percent is not None) and reference:
        p["equivalence_reference"] = F(reference, reference, evidence)
    r["participants"].append(p)
    return p


def scale(r, mid, amount, evidence):
    r["reaction_scale"] = F(f"{mid}: {amount} mmol", {"reference_molecule_id": mid, "amount": {"number": amount, "unit": "mmol"}}, evidence)


def solvent(r, mid, name, evidence, volume=None, method="original_extraction"):
    r["solvent_conditions"]["composition"] = F(name, [name], evidence, method=method)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = F(f"{volume} mL", {"number": volume, "unit": "mL"}, evidence, method=method)
    participant(r, f"p_{mid}", mid, "solvent", evidence, volume=volume)


def step(r, sid, evidence, pids, temperature=None, time=None, operation=None, method="original_extraction"):
    s = st.new_step(sid, evidence)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = F(operation, operation, evidence, method=method)
    if temperature is not None:
        norm = temperature if isinstance(temperature, dict) else temperature
        s["temperature"] = F(str(temperature), norm, evidence, method=method)
    if time is not None:
        norm = time if isinstance(time, dict) else time
        s["time"] = F(str(time), norm, evidence, method=method)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, norm, basis=None, yield_type=None, scope_text=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, norm, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = F(basis, basis, evidence)
    if scope_text:
        o["scope"] = F(scope_text, scope_text, evidence)
    r["outcomes"].append(o)
    return o


def classify(r, pairs, evidence, limitation):
    changes = []
    grouped = {}
    for i, (a, b) in enumerate(pairs, 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(change)
        lab = st.bond_label(a, b, 1)
        grouped.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    r["classification"]["bond_changes"] = [F(None, c, evidence, method="structure_derivation",
        derivation="Partial net bond inventory inferred from the source-drawn scaffold transformation; no complete atom map is asserted.") for c in changes]
    r["classification"]["labels"] = F(None, list(grouped.values()), evidence, method="structure_derivation",
        derivation="Taxonomy labels computed from the partial inferred bond changes.")
    r["classification"]["coverage"] = F(None, "partial", evidence, method="inference", derivation=limitation)
    r["classification"]["limitations"] = [limitation, "No complete atom mapping was asserted; omitted bond cleavages and proton sources remain unresolved."]
    r["classification"]["model_votes"] = []


def not_applicable_classification(r, evidence, reason):
    r["classification"]["bond_changes"] = []
    r["classification"]["labels"] = F(None, [], evidence, method="inference", derivation=reason)
    r["classification"]["coverage"] = F(None, "not_applicable", evidence, method="inference", derivation=reason)
    r["classification"]["limitations"] = [reason]
    r["classification"]["model_votes"] = []


def cascade_classification(r, evidence, kind="pyrrole"):
    pairs = [("C", "C"), ("C", "N"), ("N", "H")]
    classify(r, pairs, evidence, f"Partial {kind} cascade inventory; the compact main-text scheme does not support a unique complete atom correspondence.")


# Generic transformations supported by the page-1 overview and scope headers.
r = reaction("scheme_propargyl_to_pyrrole", 2, "generic_scheme", [], "Fig. 1c and Fig. 2 header",
             "base-promoted cascade addition–cyclization", "propargyl sulfonamide to pentasubstituted pyrrole", ev(1)+ev(2))
participant(r, "p_generic_1", "m_generic_propargyl_1", "reactant", r["evidence_ids"])
participant(r, "p_tmscn", "m_tmscn", "reagent", r["evidence_ids"], eq=3, reference="m_generic_propargyl_1")
participant(r, "p_generic_2", "m_generic_pyrrole_2", "product", r["evidence_ids"], product_role="main product class")
cascade_classification(r, r["evidence_ids"])

r = reaction("scheme_allene_to_pyrroline", 3, "generic_scheme", [], "Fig. 1c and Fig. 4 header",
             "base-promoted cascade addition–cyclization", "tetrasubstituted allenamide to hexasubstituted pyrroline", ev(1)+ev(3))
participant(r, "p_generic_4", "m_generic_allene_4", "reactant", r["evidence_ids"])
participant(r, "p_tmscn", "m_tmscn", "reagent", r["evidence_ids"], eq=3, reference="m_generic_allene_4")
participant(r, "p_generic_5", "m_generic_pyrroline_5", "product", r["evidence_ids"], product_role="main product class")
cascade_classification(r, r["evidence_ids"], "pyrroline")


# Table 1: sixteen independent optimization entries.
table_rows = [
    (1,"m_ni","m_cs2co3","m_dmf",67),(2,"m_pd","m_cs2co3","m_dmf",65),
    (3,"m_cu","m_cs2co3","m_dmf",65),(4,"m_co","m_cs2co3","m_dmf",63),
    (5,None,"m_cs2co3","m_dmf",65),(6,None,"m_kf","m_dmf","Trace"),
    (7,None,"m_k3po4","m_dmf","Trace"),(8,None,"m_k2co3","m_dmf",48),
    (9,None,"m_koh","m_dmf",52),(10,None,"m_kotbu","m_dmf",46),
    (11,None,"m_et3n","m_dmf","Trace"),(12,None,"m_cs2co3","m_acetonitrile",18),
    (13,None,"m_cs2co3","m_dme",23),(14,None,"m_cs2co3","m_toluene","Trace"),
    (15,None,"m_cs2co3","m_dce","Trace"),(16,None,"m_cs2co3","m_dioxane","Trace"),
]
names = {mid:molecules[mid]["names"]["original"]["normalized_value"] for mid in molecules}
for idx, cat, base, solv, y in table_rows:
    evidence = ev(2)
    r = reaction(f"table1_entry_{idx:02d}", 2, "specific_experiment", ["optimization"], f"Table 1 entry {idx}",
                 "cyclization reaction", "propargyl sulfonamide to pentasubstituted pyrrole", evidence)
    r["general_scheme_ids"] = [f"{DOC}:scheme_propargyl_to_pyrrole"]
    participant(r,"p_1a","m_1a","reactant",evidence,eq=1,amount=0.1,reference="m_1a")
    participant(r,"p_tmscn","m_tmscn","reagent",evidence,eq=3,amount=0.3,reference="m_1a")
    if cat:
        participant(r,"p_cat",cat,"catalyst",evidence,mol_percent=10,reference="m_1a")
    participant(r,"p_base",base,"reagent",evidence,eq=3,amount=0.3,reference="m_1a")
    solvent(r,solv,names[solv],evidence,volume=1)
    participant(r,"p_2a","m_2a","product",evidence,product_role="desired product")
    scale(r,"m_1a",0.1,evidence)
    step(r,"reaction",evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},time={"number":10,"unit":"h"},operation="heat reaction mixture")
    if isinstance(y,int):
        outcome(r,"yield_2a","yield",["m_2a"],evidence,f"{y}%",{"number":y,"unit":"%"},basis="isolated yield",yield_type="isolated yield")
    else:
        outcome(r,"trace_2a","trace",["m_2a"],evidence,"Trace","trace",basis="Table footnote: isolated yield")
    r["procedure"]["original_text"] = F("1a (0.1 mmol, 1 equiv.), TMSCN (0.3 mmol, 3 equiv.), cat. (0 or 10 mol%), base (0.3 mmol, 3 equiv.) and solvent (1 mL), at 80 °C for 10 h; isolated yield.",
                                            "Table 1 footnote conditions",evidence)
    cascade_classification(r,evidence)


# Figure 2 substrate scope (2b–2v); the compact product structures are retained by label.
scope2 = [("b",66),("c",78),("d",64),("e",64),("f",67),("g",80),("h",62),("i",53),
          ("j",81),("k",70),("l",79),("m",76),("n",72),("o",78),("p",61),("q",54),
          ("r",47),("s",67),("t",77),("u",74),("v",44)]
for letter, y in scope2:
    prod = f"m_2{letter}"
    sub = f"m_1{letter}" if letter in "tuv" else f"m_precursor_2{letter}"
    molecule(prod, f"2{letter}", 2, f"compound 2{letter}", chirality="unspecified", geometry="unspecified",
             representation="not_encoded", note="Compact Fig. 2 product drawing retained by source label; no stereochemical assignment is printed and no complete SMILES was asserted.")
    sublabel = f"1{letter}" if letter in "tuv" else f"Fig. 2 precursor corresponding to 2{letter}"
    molecule(sub, sublabel, 2, sublabel, chirality="unspecified", geometry="unspecified", representation="not_encoded",
             note="Scope precursor identity is linked to the printed product position; a complete main-text structure was not safely transcribed.")
    evidence = ev(2)
    r = reaction(f"fig2_scope_2{letter}",2,"specific_experiment",["substrate_scope"],f"Fig. 2 product 2{letter}",
                 "substrate scope of propargylamides","propargyl sulfonamide to substituted pyrrole",evidence)
    r["general_scheme_ids"]=[f"{DOC}:scheme_propargyl_to_pyrrole"]
    participant(r,"p_substrate",sub,"reactant",evidence,eq=1,amount=0.2,reference=sub)
    participant(r,"p_tmscn","m_tmscn","reagent",evidence,eq=3,amount=0.6,reference=sub)
    participant(r,"p_base","m_cs2co3","reagent",evidence,eq=3,amount=0.6,reference=sub)
    solvent(r,"m_dmf","N,N-dimethylformamide",evidence,volume=2)
    participant(r,"p_product",prod,"product",evidence,product_role="scope product")
    scale(r,sub,0.2,evidence)
    step(r,"reaction",evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},time={"number":10,"unit":"h"},operation="heat reaction mixture")
    outcome(r,f"yield_2{letter}","yield",[prod],evidence,f"{y}%",{"number":y,"unit":"%"},basis="isolated yield",yield_type="isolated yield")
    cascade_classification(r,evidence)
    r["review"]["unresolved"].append("Full substrate/product SMILES and substitution mapping require the explicitly associated SI or higher-resolution structure transcription.")

# Narrative failed scope class: exact aliphatic R2 substrate not specified.
molecule("m_r2_aliphatic_failed", "propargyl sulfonamide with aliphatic R2", 3,
         "propargyl sulfonamide with aliphatic R2", chirality="unresolved", geometry="unresolved",
         representation="not_encoded", note="Narrative reports the class only; no exact substrate label or structure is supplied.")
r = reaction("r2_aliphatic_failed",3,"specific_experiment",["substrate_scope"],"page 3 scope discussion",
             "failed substrate-scope case","no desired pyrrole formation",ev(2)+ev(3))
r["general_scheme_ids"]=[f"{DOC}:scheme_propargyl_to_pyrrole"]
participant(r,"p_substrate","m_r2_aliphatic_failed","reactant",r["evidence_ids"],eq=1,reference="m_r2_aliphatic_failed")
participant(r,"p_tmscn","m_tmscn","reagent",r["evidence_ids"],eq=3,reference="m_r2_aliphatic_failed")
participant(r,"p_base","m_cs2co3","reagent",r["evidence_ids"],eq=3,reference="m_r2_aliphatic_failed")
solvent(r,"m_dmf","N,N-dimethylformamide",r["evidence_ids"],volume=2,method="inherited_general_condition")
step(r,"reaction",r["evidence_ids"],[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},time={"number":10,"unit":"h"},operation="apply Fig. 2 scope conditions",method="inherited_general_condition")
# The narrative failure belongs to the Fig. 2 scope discussion; mark copied
# loadings as inherited rather than pretending they are reprinted on page 3.
for p in r["participants"]:
    for f in p["charges"].values():
        if f["normalized_value"] is not None:
            f["acquisition_method"] = "inherited_general_condition"
outcome(r,"no_desired_product","no_reaction",[],r["evidence_ids"],"reaction failed to provide the desired product",True,scope_text="R2 is an aliphatic group")
not_applicable_classification(r,r["evidence_ids"],"No observed product transformation is specified; no intended bond label is assigned.")
r["review"]["unresolved"].append("Exact aliphatic R2 substrate, amounts and analytical detection basis are not reported on the supplied pages.")


# Figure 3: allenamide 3a to 2a under explicitly described standard conditions.
evidence = ev(3)+ev(2)
r = reaction("fig3_3a_to_2a",3,"specific_experiment",["mechanistic_study"],"Fig. 3",
             "synthesis of substituted pyrroles from allenes","allenamide to pentasubstituted pyrrole",evidence)
r["general_scheme_ids"]=[f"{DOC}:scheme_propargyl_to_pyrrole"]
participant(r,"p_3a","m_3a","reactant",evidence,eq=1,amount=0.2,reference="m_3a")
participant(r,"p_tmscn","m_tmscn","reagent",evidence,eq=3,amount=0.6,reference="m_3a")
participant(r,"p_base","m_cs2co3","reagent",evidence,eq=3,amount=0.6,reference="m_3a")
solvent(r,"m_dmf","N,N-dimethylformamide",evidence,volume=2,method="inherited_general_condition")
participant(r,"p_2a","m_2a","product",evidence,product_role="isolated pyrrole")
scale(r,"m_3a",0.2,evidence)
step(r,"reaction",evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},time={"number":10,"unit":"h"},operation="subject allenamide to standard reaction conditions",method="inherited_general_condition")
outcome(r,"yield_2a","yield",["m_2a"],evidence,"82%",{"number":82,"unit":"%"},basis="isolated yield",yield_type="isolated yield")
cascade_classification(r,evidence)


# Figure 4: twelve tetrasubstituted allenamide scope experiments.
scope5=[("a",73),("b",76),("c",74),("d",77),("e",75),("f",69),("g",67),("h",73),("i",83),("j",62),("k",64),("l",82)]
for letter,y in scope5:
    sub=f"m_4{letter}";prod=f"m_5{letter}"
    molecule(sub,f"4{letter}",3,f"tetrasubstituted allenamide 4{letter}",chirality="unspecified",geometry="unspecified",representation="not_encoded",
             axial={"axis":"allene C=C=C","configuration":"unspecified"},note="Compact Fig. 4 scope substrate; axial configuration is not assigned and no complete SMILES was asserted.")
    molecule(prod,f"5{letter}",3,f"hexasubstituted pyrroline 5{letter}",chirality="unspecified",geometry="unspecified",representation="not_encoded",
             note="Compact Fig. 4 product drawing; no absolute/relative stereochemical descriptor is printed.",
             assignment="X-ray structure reported for 5a (CCDC 2024921); no R/S descriptor is printed." if letter=="a" else None)
    evidence=ev(3)
    r=reaction(f"fig4_scope_5{letter}",3,"specific_experiment",["substrate_scope"],f"Fig. 4 product 5{letter}",
               "substrate scope of tetrasubstituted allenamides","tetrasubstituted allenamide to hexasubstituted pyrroline",evidence)
    r["general_scheme_ids"]=[f"{DOC}:scheme_allene_to_pyrroline"]
    participant(r,"p_substrate",sub,"reactant",evidence,eq=1,amount=0.1,reference=sub)
    participant(r,"p_tmscn","m_tmscn","reagent",evidence,eq=3,amount=0.3,reference=sub)
    participant(r,"p_base","m_k2co3","reagent",evidence,eq=3,amount=0.3,reference=sub)
    solvent(r,"m_dmf","N,N-dimethylformamide",evidence,volume=1)
    participant(r,"p_product",prod,"product",evidence,product_role="scope product")
    scale(r,sub,0.1,evidence)
    step(r,"reaction",evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},time={"number":10,"unit":"h"},operation="heat reaction mixture")
    outcome(r,f"yield_5{letter}","yield",[prod],evidence,f"{y}%",{"number":y,"unit":"%"},basis="isolated yield",yield_type="isolated yield")
    cascade_classification(r,evidence,"pyrroline")
    r["review"]["unresolved"].append("Full allene/product connectivity and stereochemical assignment remain unresolved from the compact main-text graphic.")


# Figure 5 molecules and application arrows.
for mid,label,name in [
    ("m_6","6","bis-benzoylated pyrrole 6"),("m_7","7","monobenzoylated pyrrole 7"),
    ("m_8","8","lactam precursor 8"),("m_9","9","five-membered fused lactam 9"),
    ("m_10","10","lactam precursor 10"),("m_11","11","six-membered fused lactam 11"),
    ("m_amino_chain_start","Fig. 5c amino ester starting material","protected amino ester starting material"),
    ("m_propargylamine_intermediate","unlabeled 48% product","propargylamine intermediate"),
    ("m_12","12","propargyl sulfonamide 12"),("m_13","13","pyrrole analogue 13"),
    ("m_14","14","benzoylated pyrrole analogue 14"),("m_15","15","atorvastatin analogue 15"),
]:
    molecule(mid,label,3,name,chirality="unspecified",geometry="unspecified",representation="not_encoded",
             note="Structure is visible in compact Fig. 5, but a complete stereochemically faithful SMILES was not asserted from the main-text graphic.",
             assignment="X-ray structure reported (CCDC 2024922); no printed absolute stereochemical descriptor." if mid=="m_7" else None)

app_evidence=ev(3)+[f"{DOC}:page4:text"]
def application(suffix, loc, sub, prod, y, class_pairs, name_text="synthetic application"):
    r=reaction(suffix,3,"specific_experiment",["derivatization"],loc,name_text,"synthetic derivatization",app_evidence)
    participant(r,"p_substrate",sub,"reactant",app_evidence)
    participant(r,"p_product",prod,"product",app_evidence,product_role="reported product")
    outcome(r,"yield","yield",[prod],app_evidence,f"{y}%",{"number":y,"unit":"%"})
    classify(r,class_pairs,app_evidence,"Partial bond inventory from the source-drawn application arrow; no complete atom mapping was asserted.")
    r["review"]["unresolved"].append("Yield measurement basis and omitted workup/purification are not stated on the supplied pages.")
    return r

r=application("fig5a_product_6","Fig. 5a upper arrow","m_2a","m_6",76,[("C","N"),("C","N")],"protection of the free amino group")
participant(r,"p_bzcl","m_bzcl","reagent",app_evidence,eq=4,reference="m_2a")
participant(r,"p_et3n","m_et3n","reagent",app_evidence,eq=6,reference="m_2a")
solvent(r,"m_toluene","toluene",app_evidence)
step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature="reflux",operation="benzoylation")
r=application("fig5a_product_7","Fig. 5a lower arrow","m_2a","m_7",97,[("C","N")],"protection of the free amino group")
participant(r,"p_bzcl","m_bzcl","reagent",app_evidence,eq=2,reference="m_2a")
solvent(r,"m_pyridine","pyridine",app_evidence)
step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature="reflux",operation="benzoylation")

for sub,prod,y,suffix in [("m_8","m_9",44,"fig5b_lactam_9"),("m_10","m_11",47,"fig5b_lactam_11")]:
    r=application(suffix,f"Fig. 5b product {prod.split('_')[-1]}",sub,prod,y,[("C","C"),("C","N"),("N","H")],"one-pot synthesis of fused lactam")
    participant(r,"p_tmscn","m_tmscn","reagent",app_evidence,eq=3,reference=sub)
    participant(r,"p_base","m_cs2co3","reagent",app_evidence,eq=3,reference=sub)
    solvent(r,"m_dmf","N,N-dimethylformamide",app_evidence)
    step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},operation="one-pot cascade and lactam formation")

# Atorvastatin-analogue sequence: five distinct arrows.
r=application("fig5c_step1_propargylamine","Fig. 5c top arrow","m_amino_chain_start","m_propargylamine_intermediate",48,[("C","C"),("C","N")],"CuBr-mediated three-component coupling")
participant(r,"p_phenylacetylene","m_phenylacetylene","reactant",app_evidence)
participant(r,"p_4f_benzaldehyde","m_4f_benzaldehyde","reactant",app_evidence,eq=1.2,reference="m_amino_chain_start")
participant(r,"p_cubr","m_cubr","catalyst",app_evidence,mol_percent=20,reference="m_amino_chain_start")
solvent(r,"m_toluene","toluene",app_evidence)
step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature="reflux",operation="three-component coupling")

r=application("fig5c_step2_to_12","Fig. 5c vertical sulfonylation","m_propargylamine_intermediate","m_12",80,[("N","S")],"N-sulfonylation")
participant(r,"p_phso2cl","m_phso2cl","reagent",app_evidence,eq=1.5,reference="m_propargylamine_intermediate")
participant(r,"p_et3n","m_et3n","reagent",app_evidence,eq=3,reference="m_propargylamine_intermediate")
solvent(r,"m_dcm","dichloromethane",app_evidence)
step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature="r.t.",time={"number":12,"unit":"h"},operation="sulfonylation")

r=application("fig5c_step3_to_13","Fig. 5c 12 to 13","m_12","m_13",62,[("C","C"),("C","N"),("N","H")],"cascade pyrrole synthesis")
participant(r,"p_tmscn","m_tmscn","reagent",app_evidence,eq=3,reference="m_12")
participant(r,"p_base","m_cs2co3","reagent",app_evidence,eq=3,reference="m_12")
solvent(r,"m_dmf","N,N-dimethylformamide",app_evidence)
step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},operation="cascade addition–cyclization")

r=application("fig5c_step4_to_14","Fig. 5c 13 to 14","m_13","m_14",72,[("C","N")],"benzoylation")
participant(r,"p_bzcl","m_bzcl","reagent",app_evidence,eq=1.5,reference="m_13")
solvent(r,"m_pyridine","pyridine",app_evidence)
step(r,"reaction",app_evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":100,"unit":"°C"},operation="benzoylation")

r=application("fig5c_step5_to_15","Fig. 5c 14 to 15","m_14","m_15",70,[("O","H"),("O","H")],"two-stage deprotection and hydrolysis")
participant(r,"p_tfa","m_tfa","reagent",app_evidence,volume=0.5)
participant(r,"p_dcm","m_dcm","solvent",app_evidence,volume=0.5)
participant(r,"p_naoh","m_naoh","reagent",app_evidence,eq=3,reference="m_14")
participant(r,"p_methanol","m_methanol","solvent",app_evidence,volume=0.5)
r["solvent_conditions"]["composition"] = F("DCM then CH3OH",["dichloromethane","methanol"],app_evidence)
step(r,"stage_1_tfa",app_evidence,["p_substrate","p_tfa","p_dcm"],temperature="r.t.",time={"number":3,"unit":"h"},operation="TFA deprotection")
step(r,"stage_2_naoh",app_evidence,["p_substrate","p_naoh","p_methanol"],temperature="r.t.",time={"number":12,"unit":"h"},operation="NaOH hydrolysis")
r["review"]["unresolved"].append("The partial O-H classification does not claim complete atom mapping across the two-stage deprotection/hydrolysis.")


# Figure 6a radical-trapping experiment.
mech_evidence=ev(4)
r=reaction("fig6a_tempo_control",4,"specific_experiment",["control","mechanistic_study"],"Fig. 6a",
           "radical trapping reaction","propargyl sulfonamide to pyrrole in the presence of TEMPO",mech_evidence)
r["general_scheme_ids"]=[f"{DOC}:scheme_propargyl_to_pyrrole"]
participant(r,"p_1a","m_1a","reactant",mech_evidence,eq=1,reference="m_1a")
participant(r,"p_tmscn","m_tmscn","reagent",mech_evidence,eq=3,reference="m_1a")
participant(r,"p_base","m_cs2co3","reagent",mech_evidence,eq=3,reference="m_1a")
participant(r,"p_tempo","m_temppo","additive",mech_evidence,eq=3,reference="m_1a")
solvent(r,"m_dmf","N,N-dimethylformamide",mech_evidence)
participant(r,"p_2a","m_2a","product",mech_evidence,product_role="desired product")
step(r,"reaction",mech_evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},operation="radical-trapping reaction")
outcome(r,"yield_2a","yield",["m_2a"],mech_evidence,"62%",{"number":62,"unit":"%"})
cascade_classification(r,mech_evidence)
r["notes"].append(F("reaction was not inhibited", "TEMPO did not inhibit product formation", mech_evidence))
r["review"]["unresolved"].append("Reaction time and yield measurement basis are not printed for Fig. 6a.")

# Figure 6b in-situ NMR kinetic experiment.
r=reaction("fig6b_in_situ_nmr",4,"specific_experiment",["mechanistic_study"],"Fig. 6b",
           "1H NMR of crude reaction solution","in situ monitoring of 1a isomerization and pyrrole formation",mech_evidence)
participant(r,"p_1a","m_1a","reactant",mech_evidence,eq=1,reference="m_1a")
participant(r,"p_tmscn","m_tmscn","reagent",mech_evidence,eq=3,reference="m_1a")
participant(r,"p_base","m_cs2co3","reagent",mech_evidence,eq=3,reference="m_1a")
solvent(r,"m_dmso_d6","DMSO-d6",mech_evidence)
participant(r,"p_3a","m_3a","product",mech_evidence,product_role="observed intermediate")
participant(r,"p_2a","m_2a","product",mech_evidence,product_role="gradually formed product")
step(r,"monitoring",mech_evidence,[p["participant_id"] for p in r["participants"]],temperature={"number":80,"unit":"°C"},time={"number":60,"unit":"min"},operation="record in situ 1H NMR spectra at 0, 10, 20, 40 and 60 min")
outcome(r,"conversion_1a_10min","conversion",["m_1a"],mech_evidence,"1a was completely converted to 3a in 10 min",{"number":100,"unit":"%"},basis="in situ 1H NMR of crude reaction solution",scope_text="conversion at 10 min")
outcome(r,"time_course_products","other",["m_3a","m_2a"],mech_evidence,"2a was gradually produced; no other reaction intermediates were observed or identified","qualitative time-course observation",basis="in situ 1H NMR of crude reaction solution")
cascade_classification(r,mech_evidence)
r["review"]["unresolved"].append("The source does not report calibrated concentrations or a numerical product time-course on the supplied page.")

# Proposed mechanism is not a performed experiment.
r=reaction("fig6e_proposed_mechanism",4,"proposed_mechanism",[],"Fig. 6e and page 5 discussion",
           "proposed mechanism","author-proposed cyanide addition, desulfonylation, cyclization and protonation",ev(4)+ev(5))
participant(r,"p_generic_1","m_generic_propargyl_1","reactant",r["evidence_ids"])
participant(r,"p_generic_4","m_generic_allene_4","reactant",r["evidence_ids"])
participant(r,"p_generic_2","m_generic_pyrrole_2","product",r["evidence_ids"],product_role="proposed product class")
participant(r,"p_generic_5","m_generic_pyrroline_5","product",r["evidence_ids"],product_role="proposed product class")
r["notes"].append(F("Fig. 6e proposed mechanism", "Propargylamide first isomerizes to allenamide A; two cyanide additions, sulfonyl-group release, cyclization and protonation are proposed.", r["evidence_ids"]))
not_applicable_classification(r,r["evidence_ids"],"Proposed elementary steps are author interpretation, not a separate observed net experiment.")
r["review"]["unresolved"].append("Intermediates A/B/C/D and primed analogues are proposed rather than isolated or directly identified.")


# Non-synthetic context is retained separately.
bundle["source_context"] = [
    {"context_id":f"{DOC}:page1_bioactive_examples","document_id":DOC,"kind":"background_figure","evidence_ids":ev(1),"note":"Fig. 1a shows bioactive pyrrole-containing molecules; these are not experiments from this work."},
    {"context_id":f"{DOC}:page1_prior_work","document_id":DOC,"kind":"literature_context","evidence_ids":ev(1),"note":"Fig. 1b depicts cited propargylamide transformations (refs 10–12), not experiments from this work."},
    {"context_id":f"{DOC}:xray_5a_7","document_id":DOC,"kind":"characterization","evidence_ids":ev(3),"note":"X-ray depictions of 5a (CCDC 2024921) and 7 (CCDC 2024922) are characterization context, not separate synthesis experiments."},
    {"context_id":f"{DOC}:fig6c_nbo","document_id":DOC,"kind":"computational_observation","evidence_ids":ev(4),"note":"B3LYP-D3 calculations report internal-allene-carbon NBO charges 0.11 for 3b and 0.18 for 4a."},
    {"context_id":f"{DOC}:fig6d_dft","document_id":DOC,"kind":"computational_observation","evidence_ids":ev(4)+ev(5),"note":"DFT energy profile compares internal- versus terminal-carbon cyanide addition; it is not a synthetic experiment."},
    {"context_id":f"{DOC}:page5_data_availability","document_id":DOC,"kind":"source_limitation","evidence_ids":[f"{DOC}:page5:text"],"note":"The article says ESI contains experimental and analytical details, but no SI file is explicitly associated in this workspace."},
    {"context_id":f"{DOC}:page6_references","document_id":DOC,"kind":"bibliography","evidence_ids":ev(6),"note":"Page 6 contains references only and no new reaction arrow or experiment from this work."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
st.write(SUBMITTED, bundle)

# Raw visual inventory and stereo diagnostics remain separate from the formal bundle.
inventory = {
    "document_id": DOC,
    "pages": {
        "1": ["Fig. 1a background structures", "Fig. 1b cited prior work", "Fig. 1c two generic this-work transformations"],
        "2": ["Table 1 entries 1-16", "Fig. 2 scope products 2b-2v and common conditions"],
        "3": ["Fig. 3 one allenamide experiment", "Fig. 4 products 5a-5l", "Fig. 5 nine reaction arrows", "X-ray insets are characterization"],
        "4": ["Fig. 6a TEMPO control", "Fig. 6b in-situ NMR experiment", "Fig. 6c-d computational observations", "Fig. 6e proposed mechanism", "Fig. 5 prose is duplicate evidence"],
        "5": ["Proposed-mechanism continuation, conclusion, data availability, references; no new synthetic arrow"],
        "6": ["References only; no experiment from this work"],
    },
    "reaction_ids": [r["reaction_id"] for r in reactions],
    "limitations": ["Compact scope structures were not converted to guessed complete SMILES.", "No SI association was inferred."],
}
st.write(RUN/"raw_model/page1_6_visual_inventory.json", inventory)
diags = {mid: st.stereo_inspect(st.value(m["smiles"])) for mid,m in molecules.items() if st.value(m["smiles"])}
st.write(RUN/"raw_model/standard_stereo_diagnostics.json", diags)

page_records = {i: [] for i in range(1,7)}
for r in reactions:
    pages = {int(e.split(":page")[1].split(":")[0]) for e in r["evidence_ids"] if ":page" in e}
    for p in pages:
        page_records[p].append(r["reaction_id"])
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page":1,"state":"attempted","regions":[
            {"region":"Fig. 1a","finding":"Background bioactive pyrrole structures; no experiment from this work.","evidence_ids":ev(1)},
            {"region":"Fig. 1b","finding":"Cited prior-work transformations; retained as source context.","evidence_ids":ev(1)},
            {"region":"Fig. 1c","finding":"Two generic this-work transformations, represented by generic-scheme records.","evidence_ids":ev(1)}],
         "record_ids":page_records[1],"unresolved":["R-group structures are generic templates; no concrete instances are inferred from this overview."]},
        {"page":2,"state":"attempted","regions":[
            {"region":"Table 1","finding":"Sixteen optimization rows, each represented as a separate experiment.","evidence_ids":ev(2)},
            {"region":"Fig. 2","finding":"Twenty-one scope products 2b-2v with common conditions and isolated yields.","evidence_ids":ev(2)},
            {"region":"scope prose","finding":"Introduces scope and continues onto page 3; no duplicate records created.","evidence_ids":[f"{DOC}:page2:text"]}],
         "record_ids":page_records[2],"unresolved":["Most compact scope precursor/product structures lack safely transcribed complete SMILES."]},
        {"page":3,"state":"attempted","regions":[
            {"region":"scope prose continuation","finding":"Reports failed aliphatic-R2 class and narrative context for 2i-2v.","evidence_ids":[f"{DOC}:page3:text"]},
            {"region":"Fig. 3","finding":"One allenamide 3a to pyrrole 2a experiment.","evidence_ids":ev(3)},
            {"region":"Fig. 4","finding":"Twelve tetrasubstituted-allene scope experiments 5a-5l.","evidence_ids":ev(3)},
            {"region":"Fig. 5","finding":"Nine distinct synthetic-application arrows; X-ray insets are characterization, not synthesis.","evidence_ids":ev(3)}],
         "record_ids":page_records[3],"unresolved":["Fig. 5 omits several reaction times/scales and yield measurement bases; these remain unreported."]},
        {"page":4,"state":"attempted","regions":[
            {"region":"Fig. 5 discussion","finding":"Duplicate prose evidence for the application records; not counted again.","evidence_ids":[f"{DOC}:page4:text"]},
            {"region":"Fig. 6a","finding":"TEMPO radical-trapping control, product 2a reported at 62%.","evidence_ids":ev(4)},
            {"region":"Fig. 6b","finding":"In-situ 1H NMR kinetic experiment with spectra at 0-60 min.","evidence_ids":ev(4)},
            {"region":"Fig. 6c-d","finding":"NBO/DFT computational observations, retained in source_context.","evidence_ids":ev(4)},
            {"region":"Fig. 6e","finding":"Author-proposed mechanism, represented separately from performed experiments.","evidence_ids":ev(4)}],
         "record_ids":page_records[4],"unresolved":["TEMPO-control time and yield basis are not reported; kinetic curves are qualitative in the supplied graphic."]},
        {"page":5,"state":"attempted","regions":[
            {"region":"mechanism discussion and conclusions","finding":"Textual continuation of Fig. 6e and conclusion; no new synthetic experiment.","evidence_ids":[f"{DOC}:page5:text"]},
            {"region":"data availability","finding":"States that ESI exists, but no SI is explicitly associated in the workspace.","evidence_ids":[f"{DOC}:page5:text"]},
            {"region":"references","finding":"Bibliography begins; cited literature is not imported as this-work experiments.","evidence_ids":ev(5)}],
         "record_ids":page_records[5],"unresolved":["Experimental, NMR, HRMS and IR details said to be in ESI remain unavailable here."]},
        {"page":6,"state":"attempted","regions":[
            {"region":"references 7-21","finding":"Bibliography only; no reaction arrow, table row or this-work experiment.","evidence_ids":ev(6)}],
         "record_ids":page_records[6],"unresolved":[]},
    ]
}
(RUN/"coverage").mkdir(exist_ok=True)
st.write(RUN/"coverage/pages_0001_0006.json",coverage)
print(json.dumps({"molecules":len(bundle["molecules"]),"reactions":len(reactions),"source_context":len(bundle["source_context"])},ensure_ascii=False))
