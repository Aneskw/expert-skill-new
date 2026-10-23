#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for adsc.202100768 pages 1--6.

This script only assembles facts transcribed from the supplied source images/text.
It does not promote machine records to source-reviewed status.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_c690e758ccda5a7c"
MID = lambda s: f"{DOC}:mol_{s}"
RID = lambda s: f"{DOC}:{s}"
P = {i: (f"{DOC}:page{i}:image", f"{DOC}:page{i}:text") for i in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


bundle = json.loads((Path(__file__).with_name("standard_init.json")).read_text())
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


def molecule(key, label, name, evidence, smiles=None, chirality="unspecified",
             geometry="unspecified", rep_status=None, stereo_note=None):
    mid = MID(key)
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = f()
    m["names"]["aliases"] = f()
    m["cas_numbers"] = f()
    m["identity_note"] = f(
        None, f"Source identity retained as {label}; CAS was not reported on supplied pages 1–6.", evidence,
        method="inference", derivation="No CAS number was found in the inspected source regions; no external lookup was performed.")
    if smiles:
        m["representation_type"] = "smiles"
        deriv = "Candidate SMILES transcribed from the source drawing; parsing is not independent identity verification."
        m["smiles"] = f(f"source identity {label}", smiles, evidence,
                         method="structure_derivation", derivation=deriv)
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["isomeric_smiles"] = f(
            f"source identity {label}", diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation", derivation=deriv)
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="Computed from the candidate isomeric graph; not a sample identity.")
        m["stereochemistry"]["representation_status"] = f(
            rep_status or "partial", rep_status or "partial", evidence,
            method="inference", derivation="Graph was checked by the standard stereo diagnostic; axial and source-drawing correctness still require review.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        m["stereochemistry"]["isomeric_smiles"] = f()
        m["stereochemistry"]["identity_key"] = f()
        m["stereochemistry"]["representation_status"] = f(
            "not encoded", "not_encoded", evidence, method="inference",
            derivation="The compact source drawing was not converted to a sufficiently confident complete graph.")
    m["stereochemistry"]["chirality_status"] = f(
        chirality, chirality, evidence, method="inference",
        derivation="Status preserves source-specified or source-unspecified stereochemistry; absence of a descriptor was not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = f(
        geometry, geometry, evidence, method="inference",
        derivation="Status preserves source-specified or source-unspecified alkene/ring geometry.")
    m["stereochemistry"]["atom_centers"] = f()
    m["stereochemistry"]["double_bonds"] = f()
    m["stereochemistry"]["axial"] = f(
        stereo_note, stereo_note, evidence, method="inference",
        derivation="Assessment is limited to the visible source drawing and candidate graph.") if stereo_note else f()
    m["stereochemistry"]["mixture"] = f()
    m["stereochemistry"]["assignment"] = f()
    molecules[mid] = m
    return mid


def generic_molecule(key, label, name, evidence, fixed, variables):
    mid = MID(key)
    if mid in molecules:
        return mid
    molecule(key, label, name, evidence)
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, evidence),
        "variables": [{
            "label": var,
            "attachment_points": f("shown attachment point", "shown attachment point", evidence),
            "definition": f(desc, desc, evidence),
            "allowed_values": f(),
            "constraints": f(),
        } for var, desc in variables],
        "instance_mapping": [],
        "reason": f("R-group or mixture drawing", "Source gives an R-group family or non-unique mixture, not one concrete molecule.", evidence),
    }
    molecules[mid] = m
    return mid


def simple(key, label, smiles, evidence=P[3]):
    return molecule(key, label, label, evidence, smiles=smiles,
                    chirality="achiral", geometry="not_applicable", rep_status="partial")


# Recurrent concrete and operational species.
m1a = molecule("1a", "1a", "1,6-allenyne 1a", P[2] + P[3] + P[4] + P[5],
    smiles="C#Cc1ccccc1N(S(=O)(=O)c1ccc(C)cc1)C=C=C",
    chirality="achiral", geometry="not_applicable", rep_status="partial",
    stereo_note="Terminal-CH2 allene shown; this allene unit is not axially stereogenic, but toolkit cumulene diagnostics are limited.")
m2a = molecule("2a", "2a / TsN3", "tosyl azide", P[2] + P[3] + P[4] + P[5],
    smiles="Cc1ccc(S(=O)(=O)N=[N+]=[N-])cc1", chirality="achiral", geometry="not_applicable")
m2b = molecule("2b", "2b / MsN3", "methanesulfonyl azide", P[3],
    smiles="CS(=O)(=O)N=[N+]=[N-]", chirality="achiral", geometry="not_applicable")
m2c = molecule("2c", "2c / NsN3", "p-nitrobenzenesulfonyl azide", P[3],
    smiles="O=[N+]([O-])c1ccc(S(=O)(=O)N=[N+]=[N-])cc1", chirality="achiral", geometry="not_applicable")
m3a = molecule("3a", "3a", "allenamide-tethered N-sulfonyl-1,2,3-triazole 3a", P[2] + P[4])
m4a = molecule("4a", "4a", "N-allyl 1,2,3-triazoloquinazoline 4a", P[2] + P[3] + P[4] + P[5] + P[6],
    chirality="unspecified", geometry="not_applicable")

common = {
    "CuTc": generic_molecule("CuTc", "CuTc", "copper(I) thiophene-2-carboxylate", P[2] + P[3] + P[4] + P[5], "CuTc catalyst as printed", []),
    "toluene": simple("toluene", "toluene", "Cc1ccccc1"),
    "benzene": simple("benzene", "benzene", "c1ccccc1"),
    "MeCN": simple("MeCN", "MeCN", "CC#N"),
    "THF": simple("THF", "THF", "C1CCOC1"),
    "water": simple("water", "H2O", "O"),
    "SiO2": generic_molecule("SiO2", "SiO2", "silica gel", P[2] + P[3] + P[4] + P[5], "silica material; composition/particle form not uniquely specified", []),
    "MS4A": generic_molecule("MS4A", "4 Å MS", "4 Å molecular sieves", P[3], "4 Å molecular sieves; non-unique material", []),
    "xylenes": generic_molecule("xylenes", "xylenes", "xylenes", P[3], "xylene isomer mixture not specified", []),
    "H2SO4": simple("H2SO4", "H2SO4", "OS(=O)(=O)O"),
    "MsOH": simple("MsOH", "MsOH", "CS(=O)(=O)O"),
    "TsOH": simple("TsOH", "TsOH", "Cc1ccc(S(=O)(=O)O)cc1", P[3] + P[4] + P[5]),
    "AlCl3": simple("AlCl3", "AlCl3", "Cl[Al](Cl)Cl"),
    "NaHCO3": simple("NaHCO3", "saturated NaHCO3(aq)", "[Na+].O=C([O-])O", P[5]),
    "EtOAc": simple("EtOAc", "EtOAc", "CCOC(C)=O", P[5]),
    "hexanes": generic_molecule("hexanes", "hexanes", "hexanes", P[5], "commercial hexanes mixture; composition not specified", []),
    "DCM": simple("DCM", "DCM", "ClCCl", P[5]),
    "PdC": generic_molecule("PdC", "Pd/C", "palladium on carbon", P[5], "heterogeneous Pd/C catalyst; support loading not specified", []),
    "H2": simple("H2", "H2", "[H][H]", P[5]),
    "IPA": simple("IPA", "IPA", "CC(O)C", P[5]),
}

# Scope compounds; compact structures remain unresolved unless a safe graph is available.
substrates, products = {"a": m1a}, {"a": m4a}
for label in list("bcdefghijklmnopqrstuvwxyz") + ["a_prime"]:
    printed = "1a′" if label == "a_prime" else f"1{label}"
    prod = "4a′" if label == "a_prime" else f"4{label}"
    note = "Substituted allene (R2 = Me) has an unassigned allene axis." if label == "a_prime" else "Terminal-CH2 allene shown; no axial descriptor is reported."
    substrates[label] = molecule(f"1{label}", printed, f"1,6-allenyne {printed}", P[4], stereo_note=note)
    products[label] = molecule(f"4{label}", prod, f"N-allyl triazolo fused product {prod}", P[4],
        geometry="unspecified" if label == "a_prime" else "not_applicable")

# Other-scope and reduction compounds.
other = {}
for x, desc in [("b", "aliphatic allenamide-based"), ("c", "phenoxyallene-based"),
                ("d", "allenoate-based"), ("e", "unbiased-allene-based")]:
    other[f"1{x}"] = molecule(f"1{x}_prime", f"1{x}′", f"{desc} 1,6-allenyne 1{x}′", P[3] + P[4], stereo_note="Allene configuration is not assigned in the source.")
    other[f"3{x}"] = molecule(f"3{x}_prime", f"3{x}′", f"N-sulfonyl-1,2,3-triazole 3{x}′", P[3] + P[4])
    other[f"4{x}"] = molecule(f"4{x}_prime", f"4{x}′", f"desired cyclized product 4{x}′", P[4])

m3z = molecule("3z", "3z", "pyridine-derived N-sulfonyl-1,2,3-triazole 3z", P[3] + P[4])
for x in ("a", "i", "v"):
    molecule(f"5{x}", f"5{x}", f"N-alkyl triazolo fused product 5{x}", P[5])

# Generic families for cited work and the current general scheme.
g_triazole = generic_molecule("lit_triazole", "N-sulfonyl-1,2,3-triazole", "generic N-sulfonyl-1,2,3-triazole", P[1] + P[2], "1,2,3-triazole bearing R and N-sulfonyl substituents", [("R", "substituent R"), ("S", "sulfonyl substituent")])
g_olefin = generic_molecule("lit_olefin", "olefin", "generic olefin", P[2], "olefin with R1/R2 substituents", [("R1", "olefin substituent"), ("R2", "olefin substituent")])
g_iminium = generic_molecule("lit_iminium", "iminium", "generic iminium", P[2], "iminium with two R3 substituents", [("R3", "iminium substituent")])
g_vinyl_ether = generic_molecule("lit_vinyl_ether", "vinyl ether", "generic vinyl ether", P[2], "vinyl ether bearing R4/R5", [("R4", "alkene substituent"), ("R5", "alkoxy substituent")])
g_alkyne = generic_molecule("lit_alkyne", "alkyne", "generic alkyne", P[2], "alkyne bearing R6", [("R6", "alkyne substituent")])
for i, name in enumerate(["N2-alkyl triazole", "N-dialkylaminomethyl triazole", "N2-alkoxyalkyl triazole", "N1-vinyl triazole"], 1):
    generic_molecule(f"lit_product{i}", f"literature product {i}", name, P[2], f"generic {name} product drawing", [("R", "source R-group variables")])
g_1 = generic_molecule("family1", "1", "1,6-allenyne family 1", P[2] + P[4], "aniline-linked 1,6-allenyne with variables R/R1/R2", [("R", "N-protecting group"), ("R1", "aryl substituent"), ("R2", "allene-terminal substituent")])
g_2 = generic_molecule("family2", "2", "organic azide family 2", P[2] + P[3], "organic sulfonyl azide RN3", [("R", "organic sulfonyl group")])
g_4 = generic_molecule("family4", "4", "N-allyl 1,2,3-triazole-fused product family 4", P[2] + P[4], "triazolo fused product with R/R1/R2 variables", [("R", "N substituent"), ("R1", "aryl substituent"), ("R2", "allyl substituent")])


def reaction(key, evidence, record_type="specific_experiment", categories=(), location="", name=None,
             transformation=None, unresolved=None):
    r = st.new_reaction(RID(key), DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(list(categories), list(categories), evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else f()
    if transformation:
        r["transformation_type"] = f(None, transformation, evidence, method="inference",
            derivation="Transformation class inferred from the source arrow and accompanying description.")
    r["review"]["unresolved"] = list(unresolved or ["Complete structures/atom mapping and unavailable SI details require expert review."])
    return r


def participant(r, mid, role, evidence, stage="reaction", label=None, charges=None, eqref=None, product_role=None):
    pid = f"p{len(r['participants'])+1}_{label or mid.rsplit('_',1)[-1]}"
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, val in (charges or {}).items():
        raw, norm = val
        p["charges"][key] = f(raw, norm, evidence)
    if eqref:
        p["equivalence_reference"] = f(f"relative to {eqref.rsplit('_',1)[-1]}", eqref, evidence)
    r["participants"].append(p)
    return pid


def step(r, stage, evidence, operation, pids=(), temperature=None, time=None,
         atmosphere=None, reactor=None, time_status=None):
    s = st.new_step(f"s{len(r['steps'])+1}_{stage}", evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence, status=time_status)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    r["steps"].append(s)
    if reactor is not None:
        r["special_conditions"]["reactor"] = f(reactor, reactor, evidence)
    return s


def outcome(r, kind, targets, evidence, raw, norm, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(f"o{len(r['outcomes'])+1}_{kind}", kind, targets, evidence)
    o["value"] = f(raw, norm, evidence)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    r["outcomes"].append(o)
    return o


def classification(r, mode, count=1):
    ev = r["evidence_ids"]
    c = r["classification"]
    c["mapped_reaction_smiles"] = f()
    c["bond_changes"] = []
    c["model_votes"] = []
    if mode == "none":
        c["labels"] = f("no observed transformation", [], ev, method="inference", derivation="No desired/observable transformation was reported.")
        c["coverage"] = f(None, "not_applicable", ev, method="inference", derivation="No observed bond formation to classify.")
        c["limitations"] = ["Expected transformation was not assigned to a negative experiment."]
        return
    if mode == "unresolved":
        c["labels"] = f()
        c["coverage"] = f(None, "unresolved", ev, method="inference", derivation="Insufficient structural detail for a safe classification.")
        c["limitations"] = ["No bond-formation labels asserted."]
        return
    changes, labels = [], {}
    if mode == "CN":
        for i in range(count):
            cid = f"b{i+1}"
            changes.append({"change_id": cid, "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"})
            labels.setdefault("C-N", {"category_number": 5, "code": "C-N", "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    elif mode == "hydrogenation":
        changes.append({"change_id": "b1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"})
        for i in (2, 3):
            cid = f"b{i}"
            changes.append({"change_id": cid, "atoms": [{"element": "C"}, {"element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1})
            labels.setdefault("C-H", {"category_number": 1, "code": "C-H", "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    c["bond_changes"] = [f(None, ch, ev, method="inference", derivation="Class-level net change inferred from the source reaction drawing; exact atoms are not mapped.") for ch in changes]
    c["labels"] = f(None, list(labels.values()), ev, method="inference", derivation="Taxonomy labels derived from the class-level bond changes.")
    c["coverage"] = f(None, "partial", ev, method="inference", derivation="Exact atom identities, complete hydrogen accounting, coproducts, and atom mapping are unresolved.")
    c["limitations"] = ["Class-level inference only; no complete atom-mapped reaction is asserted."]


def set_scale_and_solvent(r, reference, amount, solvent_mid, volume, evidence, raw_amount=None):
    r["reaction_scale"] = f(raw_amount or f"{amount} mmol {reference.rsplit('_',1)[-1]}",
        {"reference_molecule_id": reference, "amount": q(amount, "mmol")}, evidence)
    r["solvent_conditions"]["composition"] = f(molecules[solvent_mid]["names"]["original"]["normalized_value"], [solvent_mid], evidence)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume} mL", q(volume, "mL"), evidence)


def add_general_workup(r, evidence=P[5]):
    existing = {p["molecule_id"]: p["participant_id"] for p in r["participants"]}
    ids = {}
    for key, role, stage in [("NaHCO3", "workup_agent", "workup"), ("EtOAc", "workup_agent", "workup"),
                             ("SiO2", "purification_agent", "purification"), ("hexanes", "purification_agent", "purification"),
                             ("DCM", "purification_agent", "purification")]:
        ids[key] = existing.get(common[key]) or participant(r, common[key], role, evidence, stage=stage, label=key)
    step(r, "workup", evidence, "Quench with saturated aqueous NaHCO3 and extract with EtOAc.", [ids["NaHCO3"], ids["EtOAc"]])
    step(r, "workup", evidence, "Concentrate the combined organic phase under reduced pressure.", [ids["EtOAc"]])
    step(r, "purification", evidence, "Flash chromatography on silica gel with hexanes/EtOAc/DCM eluent.", [ids["SiO2"], ids["hexanes"], ids["EtOAc"], ids["DCM"]])


reactions = []

# Four cited prior-work arrows in Scheme 1a.
lit_specs = [
    ("literature_eq1", g_olefin, MID("lit_product1"), "TsOH", "N2-selective addition to olefins"),
    ("literature_eq2", g_iminium, MID("lit_product2"), None, "N2-selective addition to iminiums"),
    ("literature_eq3", g_vinyl_ether, MID("lit_product3"), None, "gold-catalyzed N2-selective addition to vinyl ethers"),
    ("literature_eq4", g_alkyne, MID("lit_product4"), None, "gold-catalyzed N1-selective addition to alkynes"),
]
for key, co, prod, reagent_key, desc in lit_specs:
    r = reaction(key, P[1] + P[2], "literature_example", (), "Scheme 1a", desc, desc)
    participant(r, g_triazole, "reactant", P[2], label="triazole")
    participant(r, co, "reactant", P[2], label="co_reactant")
    participant(r, prod, "product", P[2], label="product", product_role="reported literature product family")
    if reagent_key:
        participant(r, common[reagent_key], "reagent", P[2], label=reagent_key)
    classification(r, "CN")
    reactions.append(r)

# Current-paper general transformation.
r = reaction("generic_tandem_annulation", P[2] + P[4], "generic_scheme", (), "Scheme 1b-c and Table 2 header",
             "one-pot tandem annulation of 1,6-allenynes with organic azides", "CuAAC/intramolecular annulation")
participant(r, g_1, "reactant", P[2] + P[4], label="family1")
participant(r, g_2, "reactant", P[2] + P[4], label="family2")
participant(r, g_4, "product", P[2] + P[4], label="family4", product_role="generic desired product family")
classification(r, "CN", 3)
reactions.append(r)
generic_id = r["reaction_id"]

# Table 1 optimization rows.
opt = {
  1:(m2a,"SiO2","toluene",68,"default"), 2:(m2a,None,"toluene",None,"default"),
  3:(m2a,"SiO2","toluene",None,"ms"), 4:(m2b,"SiO2","toluene",42,"default"),
  5:(m2c,"SiO2","toluene",None,"default"), 6:(m2a,"SiO2","benzene",56,"default"),
  7:(m2a,"SiO2","xylenes",61,"default"), 8:(m2a,"SiO2","MeCN",None,"default"),
  9:(m2a,"SiO2","THF",None,"default"), 10:(m2a,"H2SO4","toluene",None,"acid5"),
  11:(m2a,"MsOH","toluene",55,"acid5"), 12:(m2a,"TsOH","toluene",69,"acid5"),
  13:(m2a,"AlCl3","toluene",39,"acid5"), 14:(m2a,"TsOH","toluene",47,"low6h"),
  15:(m2a,"TsOH","toluene",71,"long12h"), 16:(m2a,"TsOH","toluene",0,"hot"),
}
for n, (azide, cat_key, solv_key, yld, variant) in opt.items():
    ev = P[3] + (P[2] if n in (1,2,3,4,5,8,9,10,11,12,13,14,15,16) else ())
    r = reaction(f"table1_entry{n}", ev, categories=("optimization",), location=f"Table 1, entry {n}",
                 name="optimization of one-pot tandem annulation", transformation="CuAAC/intramolecular annulation",
                 unresolved=["Exact atom mapping and full byproduct accounting are unavailable."])
    r["general_scheme_ids"] = [generic_id]
    p1 = participant(r, m1a, "reactant", ev, label="1a", charges={"amount":("0.1 mmol",q(0.1,"mmol"))})
    p2 = participant(r, azide, "reactant", ev, label="2", charges={"amount":("0.1 mmol",q(0.1,"mmol"))})
    pc = participant(r, common["CuTc"], "catalyst", ev, label="CuTc", charges={"mol_percent":("10 mol%",q(10,"mol%"))}, eqref=m1a)
    ps = participant(r, common[solv_key], "solvent", ev, label=solv_key, charges={"volume":("1.0 mL",q(1.0,"mL"))})
    stage1 = [p1,p2,pc,ps]
    step(r,"reaction",ev,"Combine 1a, organic azide 2, CuTc, and solvent; stir.",stage1,("rt","rt"),("12 h",q(12,"h")))
    stage2 = []
    if cat_key:
        charges = {}
        if cat_key == "SiO2": charges["mass"]=("200 mg",q(200,"mg"))
        elif variant == "low6h": charges["mol_percent"]=("20 mol%",q(20,"mol%"))
        else: charges["mol_percent"]=("100 mol%",q(100,"mol%"))
        stage2.append(participant(r,common[cat_key],"catalyst",ev,label=cat_key,charges=charges,eqref=m1a if "mol_percent" in charges else None))
    if variant == "ms":
        stage2.append(participant(r,common["MS4A"],"additive",ev,label="MS4A",charges={"mass":("100 mg",q(100,"mg"))}))
    else:
        stage2.append(participant(r,common["water"],"additive",ev,label="water",charges={"amount":("0.1 mmol",q(0.1,"mmol"))},eqref=m1a))
    time2 = ("24 h",q(24,"h")); temp2=("rt","rt")
    if variant == "acid5": time2=("5 min",q(5,"min"))
    elif variant == "low6h": time2=("6 h",q(6,"h"))
    elif variant == "long12h": time2=("12 h",q(12,"h"))
    elif variant == "hot": temp2=("100 °C",q(100,"°C"))
    step(r,"reaction",ev,"Add the second-stage catalyst/additive and continue the reaction.",stage2,temp2,time2)
    set_scale_and_solvent(r,m1a,0.1,common[solv_key],1.0,ev)
    if yld is not None:
        pp = participant(r,m4a,"product",ev,label="4a",product_role="desired product; 0% means not observed" if yld==0 else "main product")
        outcome(r,"yield",[m4a],ev,f"{yld}%",q(yld,"%"),basis="isolated",yield_type="isolated yield")
        if yld == 0:
            outcome(r,"no_reaction",[],ev,"no desired product","no desired product")
            r["notes"].append(f("the allenyne substrate is unstable at high temperature", "the allenyne substrate is unstable at high temperature", P[2]))
            classification(r,"none")
        else:
            classification(r,"CN",3)
    else:
        outcome(r,"no_reaction",[],ev,"table dash; prose reports no product/not suitable/did not work","no desired product")
        classification(r,"none")
    reactions.append(r)


def add_standard_scope_reaction(key, sub_mid, prod_mid, yld, label, extra_evidence=(),
                                stage2_catalyst="TsOH", stage2_time=None, temp="rt", amount=0.2,
                                isolated=True, categories=("substrate_scope",), include_workup=True):
    ev = P[4] + extra_evidence + (P[5] if include_workup and stage2_catalyst=="TsOH" else ())
    r = reaction(key, ev, categories=categories, location=label,
                 name="one-pot tandem annulation", transformation="CuAAC/intramolecular annulation")
    r["general_scheme_ids"]=[generic_id]
    p1=participant(r,sub_mid,"reactant",ev,label="substrate",charges={"amount":(f"{amount} mmol",q(amount,"mmol")),"equivalents":("1.0 equiv",q(1,"equiv"))},eqref=sub_mid)
    p2=participant(r,m2a,"reactant",ev,label="2a",charges={"amount":(f"{amount} mmol",q(amount,"mmol")),"equivalents":("1.0 equiv",q(1,"equiv"))},eqref=sub_mid)
    pc=participant(r,common["CuTc"],"catalyst",ev,label="CuTc",charges={"mol_percent":("10 mol%",q(10,"mol%"))},eqref=sub_mid)
    ps=participant(r,common["toluene"],"solvent",ev,label="toluene",charges={"volume":("2.0 mL",q(2,"mL"))})
    step(r,"reaction",ev,"Add tosyl azide to a solution of the 1,6-allenyne and CuTc in toluene; stir.",[p1,p2,pc,ps],("rt","rt"),("12 h",q(12,"h")))
    stage2=[]
    if stage2_catalyst=="TsOH":
        pw=participant(r,common["water"],"additive",ev,label="water",charges={"amount":(f"{amount} mmol",q(amount,"mmol")),"equivalents":("1.0 equiv",q(1,"equiv"))},eqref=sub_mid)
        pa=participant(r,common["TsOH"],"catalyst",ev,label="TsOH",charges={"mol_percent":("100 mol%",q(100,"mol%"))},eqref=sub_mid)
        stage2=[pw,pa]
        operation="Successively add H2O and TsOH and continue until TLC completion."
    else:
        pa=participant(r,common["SiO2"],"catalyst",ev,label="SiO2",charges={"mass":("200 mg",q(200,"mg"))})
        stage2=[pa]; operation="Add SiO2 and continue the comparison reaction."
    t = stage2_time
    if t is None:
        time_fact=("5 min–48 h",None)
        s2=step(r,"reaction",ev,operation,stage2,(temp,q(float(temp),"°C") if isinstance(temp,(int,float)) else temp),time_fact,time_status="unreadable")
        s2["time"]["derivation"]="The caption gives a collective 5 min–48 h range but does not assign an exact duration to this substrate."
    else:
        step(r,"reaction",ev,operation,stage2,(f"{temp} °C",q(temp,"°C")) if isinstance(temp,(int,float)) else (temp,temp),stage2_time)
    set_scale_and_solvent(r,sub_mid,amount,common["toluene"],2.0,ev)
    participant(r,prod_mid,"product",ev,label="product",product_role="main product")
    outcome(r,"yield",[prod_mid],ev,f"{yld}%",q(yld,"%"),basis="isolated" if isolated else None,yield_type="isolated yield" if isolated else None)
    if include_workup and stage2_catalyst=="TsOH":
        add_general_workup(r,P[5])
        r["procedure"]["general_procedure_reference"] = f("General Procedure for the Preparation of Compounds 4", "General Procedure for the Preparation of Compounds 4", P[5])
        r["procedure"]["original_text"] = f("Page 5 general procedure", "Tosyl azide was added to allenyne/CuTc in toluene; rt 12 h; H2O then TsOH; TLC; NaHCO3 quench; EtOAc extraction; concentration; flash chromatography.", P[5])
        r["procedure"]["addition_order"] = f("2a into 1/CuTc/toluene; then H2O; then TsOH", ["2a into 1/CuTc/toluene","H2O","TsOH"], P[5])
    classification(r,"CN",3)
    reactions.append(r)
    return r


# Table 2 principal yields, excluding z which is explicitly a two-arrow sequence.
scope_yields = {
 "a":69,"b":59,"c":83,"d":57,"e":89,"f":22,"g":41,"h":61,"i":60,"j":81,
 "k":79,"l":37,"m":71,"n":60,"o":53,"p":68,"q":76,"r":46,"s":47,"t":42,
 "u":57,"v":67,"w":53,"x":74,"y":61,"a_prime":75,
}
for lab,y in scope_yields.items():
    printed="4a′" if lab=="a_prime" else f"4{lab}"
    add_standard_scope_reaction(f"table2_{printed.replace('′','prime')}",substrates[lab],products[lab],y,
                                f"Table 2, {printed}")

# Four explicitly parenthesized SiO2 comparison experiments.
for lab,y in {"j":24,"k":26,"q":46,"x":16}.items():
    add_standard_scope_reaction(f"table2_4{lab}_sio2",substrates[lab],products[lab],y,
        f"Table 2, 4{lab}, parenthesized SiO2 comparison", stage2_catalyst="SiO2",
        stage2_time=("24 h",q(24,"h")), include_workup=False)

# 5 mmol scale experiment is distinct from the 0.2 mmol Table 2 4a row.
r = reaction("table2_4a_scaleup", P[3] + P[4], categories=("scale_up",), location="Table 2, 4a parenthesized 65% and footnote c",
             name="gram-scale tandem annulation", transformation="CuAAC/intramolecular annulation",
             unresolved=["Only the 5 mmol scale and standard-condition statement are explicit; scaled solvent volume and exact stage-2 time are not reported."])
r["general_scheme_ids"]=[generic_id]
p1=participant(r,m1a,"reactant",P[3]+P[4],label="1a",charges={"amount":("5 mmol scale",q(5,"mmol")),"equivalents":("1.0 equiv",q(1,"equiv"))},eqref=m1a)
p2=participant(r,m2a,"reactant",P[3]+P[4],label="2a",charges={"equivalents":("1.0 equiv under standard conditions",q(1,"equiv"))},eqref=m1a)
pc=participant(r,common["CuTc"],"catalyst",P[3]+P[4],label="CuTc",charges={"mol_percent":("10 mol% under standard conditions",q(10,"mol%"))},eqref=m1a)
pw=participant(r,common["water"],"additive",P[3]+P[4],label="water",charges={"equivalents":("1.0 equiv under standard conditions",q(1,"equiv"))},eqref=m1a)
pa=participant(r,common["TsOH"],"catalyst",P[3]+P[4],label="TsOH",charges={"mol_percent":("100 mol% under standard conditions",q(100,"mol%"))},eqref=m1a)
step(r,"reaction",P[3]+P[4],"Perform the reported 5 mmol scale synthesis under standard reaction conditions.",[p1,p2,pc,pw,pa],("rt","rt"))
r["reaction_scale"]=f("5 mmol scale",{"reference_molecule_id":m1a,"amount":q(5,"mmol")},P[4])
participant(r,m4a,"product",P[3]+P[4],label="4a",product_role="main product")
outcome(r,"yield",[m4a],P[3]+P[4],"65%",q(65,"%"),basis="isolated",yield_type="isolated yield")
classification(r,"CN",3);reactions.append(r)

# Pyridine row: two actual arrows, precursor isolation then elevated-temperature cyclization.
r = reaction("table2_3z_precursor", P[3]+P[4], categories=("substrate_scope",), location="Table 2, dashed 3z/4z box, first arrow",
             name="CuAAC preparation of 3z", transformation="CuAAC")
p1=participant(r,substrates["z"],"reactant",P[4],label="1z",charges={"amount":("0.2 mmol",q(0.2,"mmol"))})
p2=participant(r,m2a,"reactant",P[4],label="2a",charges={"amount":("0.2 mmol",q(0.2,"mmol"))})
pc=participant(r,common["CuTc"],"catalyst",P[4],label="CuTc",charges={"mol_percent":("10 mol%",q(10,"mol%"))},eqref=substrates["z"])
ps=participant(r,common["toluene"],"solvent",P[4],label="toluene",charges={"volume":("2.0 mL",q(2,"mL"))})
step(r,"reaction",P[4],"CuTc-catalyzed reaction of pyridine-derived allenyne 1z with 2a in toluene.",[p1,p2,pc,ps],("rt","rt"),("12 h",q(12,"h")))
set_scale_and_solvent(r,substrates["z"],0.2,common["toluene"],2.0,P[4])
participant(r,m3z,"product",P[3]+P[4],label="3z",product_role="isolated intermediate")
outcome(r,"yield",[m3z],P[3]+P[4],"91% (described as almost quantitative)",q(91,"%"),basis="isolated",yield_type="isolated yield")
classification(r,"CN",2);reactions.append(r)

r = reaction("table2_4z_cyclization", P[3]+P[4], categories=("substrate_scope",), location="Table 2, dashed 3z→4z arrow",
             name="elevated-temperature cyclization of 3z", transformation="intramolecular annulation",
             unresolved=["Reaction time and detailed procedure for the elevated-temperature step are not reported on pages 1–6."])
p3=participant(r,m3z,"reactant",P[4],label="3z")
pa=participant(r,common["TsOH"],"catalyst",P[4],label="TsOH",charges={"mol_percent":("100 mol%",q(100,"mol%"))},eqref=m3z)
ps=participant(r,common["toluene"],"solvent",P[4],label="toluene")
step(r,"reaction",P[4],"Treat isolated 3z with TsOH in toluene at elevated temperature.",[p3,pa,ps],("100 °C",q(100,"°C")))
participant(r,products["z"],"product",P[4],label="4z",product_role="main product")
outcome(r,"yield",[products["z"]],P[3]+P[4],"41%",q(41,"%"),basis="isolated",yield_type="isolated yield")
classification(r,"CN");reactions.append(r)

# Scheme 2: isolate the CuAAC intermediates and test each cyclization separately.
prec_yields={"b":90,"c":93,"d":83,"e":72}
for x,y in prec_yields.items():
    r=reaction(f"scheme2_3{x}_precursor",P[3]+P[4],categories=("substrate_scope",),location=f"Scheme 2{x}, first arrow",
               name=f"CuAAC preparation of 3{x}′",transformation="CuAAC")
    p1=participant(r,other[f"1{x}"],"reactant",P[4],label=f"1{x}prime")
    p2=participant(r,m2a,"reactant",P[4],label="2a")
    pc=participant(r,common["CuTc"],"catalyst",P[4],label="CuTc",charges={"mol_percent":("10 mol%",q(10,"mol%"))},eqref=other[f"1{x}"])
    ps=participant(r,common["toluene"],"solvent",P[4],label="toluene")
    step(r,"reaction",P[4],"CuTc-catalyzed CuAAC reaction in toluene.",[p1,p2,pc,ps],("rt","rt"))
    participant(r,other[f"3{x}"],"product",P[4],label=f"3{x}prime",product_role="isolated intermediate")
    basis="isolated" if x in ("b","c") else None
    outcome(r,"yield",[other[f"3{x}"]],P[3]+P[4],f"{y}%",q(y,"%"),basis=basis,yield_type="isolated yield" if basis else None)
    classification(r,"CN",2);reactions.append(r)

for x in ("b","c","d","e"):
    r=reaction(f"scheme2_4{x}_cyclization",P[3]+P[4]+(P[6] if x=="b" else ()),categories=("substrate_scope",),location=f"Scheme 2{x}, second arrow",
               name=f"TsOH treatment of 3{x}′",transformation="intramolecular annulation",
               unresolved=["Reaction time, scale, detailed workup, and product measurement method are not reported on pages 1–6."])
    p3=participant(r,other[f"3{x}"],"reactant",P[4],label=f"3{x}prime")
    pa=participant(r,common["TsOH"],"catalyst",P[4],label="TsOH",charges={"mol_percent":("100 mol%",q(100,"mol%"))},eqref=other[f"3{x}"])
    ps=participant(r,common["toluene"],"solvent",P[4],label="toluene")
    step(r,"reaction",P[4],"Treat isolated triazole with TsOH in toluene.",[p3,pa,ps],("rt","rt"))
    if x=="c":
        participant(r,other["4c"],"product",P[4],label="4cprime",product_role="main product")
        outcome(r,"yield",[other["4c"]],P[3]+P[4],"25%",q(25,"%"))
        classification(r,"CN")
    else:
        outcome(r,"no_reaction",[],P[3]+P[4],"desired product not detected","desired product not detected")
        if x=="b":
            outcome(r,"recovery",[other["3b"]],P[6],"recovered almost quantitatively","almost quantitative")
            outcome(r,"other",[],P[3],"complex mixture","complex mixture")
            r["notes"].append(f("no desired product or desulfonation product was detected", "no desired product or desulfonation product was detected", P[6]))
        classification(r,"none")
    reactions.append(r)

# Two named but detail-free mechanistic studies are retained without inventing SI facts.
for key,name in [("mechanistic_deuterium_study","deuterium-labeling experiment"),("mechanistic_insitu_nmr","in-situ NMR investigation")]:
    r=reaction(key,P[4],record_type="unresolved",categories=("mechanistic_study",),location="Page 4 prose; details said to be in Supporting Information",
               name=name,transformation=None,
               unresolved=["The paper mentions the study but the Supporting Information was not supplied; substrates, conditions, observations, and isotope positions are unreported in the provided pages."])
    r["notes"].append(f(name,name,P[4]))
    classification(r,"unresolved");reactions.append(r)

# Proposed mechanism is one non-executed record.
r=reaction("scheme3_proposed_mechanism",P[3]+P[4],record_type="proposed_mechanism",location="Scheme 3",
           name="possible reaction mechanisms (paths a and b)",transformation="author-proposed tandem CuAAC/cyclization mechanism",
           unresolved=["Intermediates I–IV and paths a/b are proposals, not isolated experiments; exact atom mapping is not asserted."])
participant(r,m1a,"reactant",P[4],label="1a");participant(r,m2a,"reactant",P[4],label="2a")
participant(r,m3a,"unknown",P[4],label="3a");participant(r,m4a,"product",P[4],label="4a",product_role="mechanistic endpoint")
classification(r,"unresolved");reactions.append(r)

# Scheme 4 contains three distinct hydrogenations.
for x,y in [("a",95),("i",97),("v",97)]:
    start=products[x]; end=MID(f"5{x}")
    r=reaction(f"scheme4_{x}_hydrogenation",P[5],categories=("derivatization",),location=f"Scheme 4, 4{x}→5{x}",
               name="Pd/C-catalyzed chemoselective hydrogenation",transformation="alkene hydrogenation",
               unresolved=["Reaction scale, time, workup, purification, and yield measurement basis are not reported on page 5."])
    ps=participant(r,start,"reactant",P[5],label=f"4{x}")
    pc=participant(r,common["PdC"],"catalyst",P[5],label="PdC",charges={"mol_percent":("50 mol%",q(50,"mol%"))},eqref=start)
    ph=participant(r,common["H2"],"reagent",P[5],label="H2")
    pi=participant(r,common["IPA"],"solvent",P[5],label="IPA")
    step(r,"reaction",P[5],"Hydrogenate under an H2 balloon with Pd/C in IPA.",[ps,pc,ph,pi],("rt","rt"),atmosphere="H2 (balloon)")
    r["solvent_conditions"]["composition"]=f("IPA",[common["IPA"]],P[5])
    participant(r,end,"product",P[5],label=f"5{x}",product_role="main product")
    outcome(r,"yield",[end],P[5],f"{y}%",q(y,"%"))
    classification(r,"hydrogenation");reactions.append(r)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
  {"context_id":RID("context_page1_biological_examples"),"document_id":DOC,"kind":"non_synthetic_figure",
   "description":"Figure 1 shows biological triazoloazine examples; it is not a current-paper synthesis experiment.","evidence_ids":list(P[1])},
  {"context_id":RID("context_4a_xray"),"document_id":DOC,"kind":"structure_assignment",
   "description":"Product 4a was assigned by X-ray diffraction; CCDC 2076608 is cited, but the CIF was not supplied.","evidence_ids":list(P[2]+P[6])},
  {"context_id":RID("context_characterization_inventory"),"document_id":DOC,"kind":"negative_inventory",
   "description":"Pages 1–6 contain no standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block to convert into a synthesis record.","evidence_ids":[P[i][0] for i in range(1,7)]},
  {"context_id":RID("context_references"),"document_id":DOC,"kind":"references",
   "description":"Pages 5–6 contain acknowledgements/references; cited papers and external CCDC data were not imported as current experiments.","evidence_ids":list(P[5]+P[6])},
]

out = ROOT / "runs/current/submitted.json"
st.write(out,bundle)

# Raw diagnostics and a page/arrow inventory remain separate from the formal bundle.
diagnostics={}
for mid,m in molecules.items():
    smi=st.value(m["smiles"])
    if smi:
        diagnostics[mid]=st.stereo_inspect(smi)
st.write(ROOT/"runs/current/raw_model/stereo_checks.json",diagnostics)

inventory={
 "document_id":DOC,"pages":[
  {"page":1,"items":["title/abstract/SI availability notice","Figure 1 biological context","introduction and four cited prior transformations described in prose"]},
  {"page":2,"items":["Scheme 1a four literature arrows","Scheme 1b current proposal","Scheme 1c discovery/entry-1 duplicate","optimization discussion"]},
  {"page":3,"items":["Table 1 entries 1–16 and footnotes","prose for Table 2, gram scale, Scheme 2, proposed mechanism"]},
  {"page":4,"items":["Table 2 27 named outcomes plus four SiO2 comparisons and one 5 mmol result","Scheme 2 eight arrows","Scheme 3 proposed mechanism","two SI-only mechanistic-study mentions"]},
  {"page":5,"items":["Scheme 4 three hydrogenation arrows","General Procedure for compounds 4","acknowledgements/references"]},
  {"page":6,"items":["references","CCDC 2076608 notice","footnotes 18–19 including 3b-prime recovery"]}
 ],"reaction_count":len(reactions),"molecule_count":len(molecules),
 "limitations":["SI was advertised but not supplied or associated.","Compact complex structures remain unresolved rather than guessed."]}
st.write(ROOT/"runs/current/raw_model/page_extraction.json",inventory)

def page_record_ids(page):
    image_id = P[page][0]
    return [r["reaction_id"] for r in reactions if image_id in r["evidence_ids"]]

coverage={"document_id":DOC,"semantic_accuracy_verified":False,"pages":[
 {"page":1,"state":"attempted","regions":[
   {"region":"title, abstract, keywords, and SI notice","finding":"Current-paper summary and an online-SI availability notice; no SI file was supplied or associated.","evidence_ids":[P[1][0]]},
   {"region":"Figure 1 and introduction, lower half","finding":"Biological scaffold context plus prose describing four cited prior transformations; Figure 1 is not a synthesis experiment.","evidence_ids":[P[1][0],P[1][1]]}],
  "record_ids":page_record_ids(1),"unresolved":["The advertised Supporting Information was not supplied."]},
 {"page":2,"state":"attempted","regions":[
   {"region":"Scheme 1a, upper left","finding":"Four cited prior-work reaction arrows, retained as literature examples.","evidence_ids":[P[2][0]]},
   {"region":"Scheme 1b-c, left center","finding":"Current generic proposal and discovery experiment; Scheme 1c is the same 68% SiO2 experiment as Table 1 entry 1 and was not duplicated.","evidence_ids":[P[2][0],P[2][1]]},
   {"region":"right-column optimization discussion","finding":"Narrative outcomes for Table 1 entries, including no-product controls and high-temperature instability.","evidence_ids":[P[2][0],P[2][1]]}],
  "record_ids":page_record_ids(2),"unresolved":["Scheme 1a cited examples do not provide full procedures or yields on this page."]},
 {"page":3,"state":"attempted","regions":[
   {"region":"Table 1 and footnotes, upper half","finding":"All 16 optimization rows and footnotes were separately recorded; dash cells remain qualitative no-product results.","evidence_ids":[P[3][0],P[3][1]]},
   {"region":"lower-left scope prose","finding":"Pyridine 3z/4z sequence, 4a-prime, gram-scale 4a, and Scheme 2 precursor/cyclization outcomes are described.","evidence_ids":[P[3][0],P[3][1]]},
   {"region":"lower-right mechanism prose","finding":"Two proposed paths are discussed; these are not executed intermediate syntheses.","evidence_ids":[P[3][0],P[3][1]]}],
  "record_ids":page_record_ids(3),"unresolved":["Individual Table 2 second-stage times are not assigned by the collective 5 min–48 h range."]},
 {"page":4,"state":"attempted","regions":[
   {"region":"Table 2, left column","finding":"All printed 4a–4z and 4a-prime outcomes, four parenthesized SiO2 comparisons, the 5 mmol result, and the separate 3z arrow were inventoried.","evidence_ids":[P[4][0],P[4][1]]},
   {"region":"Scheme 2, upper right","finding":"Four first-stage triazole isolations and four distinct TsOH-treatment arrows, including three failures, were recorded separately.","evidence_ids":[P[4][0],P[4][1]]},
   {"region":"Scheme 3 and lower-right prose","finding":"One proposed mechanism plus mentions of a deuterium-labeling study and in-situ NMR study; SI-only details remain unresolved.","evidence_ids":[P[4][0],P[4][1]]}],
  "record_ids":page_record_ids(4),"unresolved":["Compact complex structures and exact atom mappings remain unverified.","Mechanistic study details are unavailable without SI."]},
 {"page":5,"state":"attempted","regions":[
   {"region":"Scheme 4, upper left","finding":"Three distinct Pd/C-H2 hydrogenations (4a/4i/4v to 5a/5i/5v) with 95/97/97% yields.","evidence_ids":[P[5][0],P[5][1]]},
   {"region":"Experimental Section, left center","finding":"General procedure for compounds 4, including all charges, addition order, TLC endpoint, quench, extraction, concentration, and chromatography.","evidence_ids":[P[5][0],P[5][1]]},
   {"region":"acknowledgements and references","finding":"Non-experimental material; no synthesis record created.","evidence_ids":[P[5][0]]}],
  "record_ids":page_record_ids(5),"unresolved":["Scheme 4 does not report reaction time, scale, workup, purification, or yield measurement basis."]},
 {"page":6,"state":"attempted","regions":[
   {"region":"references [6]–[17]","finding":"Bibliography only; cited papers were not imported as current-paper experiments.","evidence_ids":[P[6][0],P[6][1]]},
   {"region":"reference [16] CCDC notice","finding":"CCDC 2076608 supports 4a structure assignment; no CIF was supplied.","evidence_ids":[P[6][0],P[6][1]]},
   {"region":"footnotes [18]–[19], lower right","finding":"Source limitations and near-quantitative recovery of 3b-prime with no desired/desulfonated product; linked to the existing Scheme 2b cyclization record.","evidence_ids":[P[6][0],P[6][1]]}],
  "record_ids":page_record_ids(6),"unresolved":["External CCDC data and Supporting Information were not supplied."]}
]}
st.write(ROOT/"runs/current/coverage/pages_0001_0006.json",coverage)

print(json.dumps({"reactions":len(reactions),"molecules":len(molecules)},ensure_ascii=False))
