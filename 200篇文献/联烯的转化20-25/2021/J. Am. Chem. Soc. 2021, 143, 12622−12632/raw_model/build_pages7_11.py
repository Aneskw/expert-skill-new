#!/usr/bin/env python3
"""Visual extraction for source PDF pages 7--11 (JACS 2021, 143, 12622--12632).

Only Scheme 8 contains new present-paper reaction arrows in this page range.
Dense source structures remain unresolved unless their identity is already in
the v3.1 registry; the script never treats RDKit parsing as source review.
"""
from __future__ import annotations

import copy
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_f2a3e97f1c7189a7"
DATASET = ROOT / "runs/current/dataset.json"
SUBMITTED = ROOT / "runs/current/submitted.json"


def ev(page: int) -> str:
    return f"{DOC}:page{page}:image"


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


dataset = st.read(DATASET)
bundle = st.empty(dataset["documents"])
existing = {m["molecule_id"]: m for m in dataset["standard"]["molecules"]}
molecules = {}


def copy_existing(key: str, evidence):
    mid = f"{DOC}:mol:{key}"
    if mid not in existing:
        raise KeyError(f"existing molecule not found: {mid}")
    if mid not in molecules:
        molecules[mid] = copy.deepcopy(existing[mid])
    m = molecules[mid]
    for occ in m["occurrences"]:
        if occ["document_id"] == DOC:
            occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + list(evidence)))
            occ["compound_label"]["evidence_ids"] = list(dict.fromkeys(
                occ["compound_label"]["evidence_ids"] + list(evidence)))
    return mid


def molecule(label, evidence, key=None, assignment=None, axial=None, centers=None,
             geometry="unresolved", geometry_raw=None, note=None):
    mid = f"{DOC}:mol:{key or safe(label)}"
    evidence = list(dict.fromkeys(evidence))
    if mid in molecules:
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    m["names"]["english"] = F(label, label, evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = F(evidence=evidence)
    m["generic_template"] = None
    m["cas_numbers"] = F(evidence=evidence)
    for field in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                  "isomeric_smiles", "identity_key"):
        m["stereochemistry"][field] = F(evidence=evidence)
    if assignment:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, evidence)
        m["stereochemistry"]["chirality_status"] = F(assignment, "specified", evidence)
        if axial:
            m["stereochemistry"]["axial"] = F(
                axial, {"axis": "allene", "configuration": axial}, evidence)
        if centers:
            m["stereochemistry"]["atom_centers"] = F(centers, centers, evidence)
        m["stereochemistry"]["representation_status"] = F(
            "source stereochemistry retained textually; complete isomeric SMILES not safely encoded",
            "not_encoded", evidence)
    else:
        m["stereochemistry"]["chirality_status"] = F(
            "configuration is not specified for this Scheme 8 structure", "unresolved", evidence)
        m["stereochemistry"]["representation_status"] = F(
            "complete structure and stereochemistry not safely encoded", "unresolved", evidence)
    graw = geometry_raw or (
        "no applicable alkene E/Z or ring cis/trans assignment in the displayed product"
        if geometry == "not_applicable" else
        "alkene geometry or ring relationship not independently assigned from the dense drawing")
    m["stereochemistry"]["geometry_status"] = F(graw, geometry, evidence)
    m["identity_note"] = F(
        note or "Source label and drawing retained; complete connectivity was not safely converted to SMILES from the main-text figure.",
        note or "Source label and drawing retained; complete connectivity was not safely converted to SMILES from the main-text figure.",
        evidence)
    m["external_properties"] = []
    molecules[mid] = m
    return mid


def reagent(label, smiles, evidence, key=None):
    mid = molecule(label, evidence, key=key, geometry="not_applicable",
                   note="Named Scheme 8 reagent/solvent; candidate connectivity is retained as a derived representation, not source identity verification.")
    m = molecules[mid]
    diag = st.stereo_inspect(smiles)
    derivation = ("Candidate connectivity generated from the unambiguous source reagent name and checked by RDKit; "
                  "parsing and graph diagnostics do not constitute source verification.")
    m["representation_type"] = "smiles"
    m["smiles"] = F(label, smiles, evidence, method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["chirality_status"] = F(
        "no potential stereocenter in candidate graph", "achiral", evidence,
        method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["geometry_status"] = F(
        "no applicable E/Z or ring cis/trans assignment", "not_applicable", evidence,
        method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["representation_status"] = F(
        "candidate graph stereochemistry fully represented", "encoded", evidence,
        method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["isomeric_smiles"] = F(
        label, diag["canonical_isomeric_smiles"], evidence,
        method="structure_derivation", derivation=derivation)
    m["stereochemistry"]["identity_key"] = F(
        label, diag["stereo_graph_key"], evidence,
        method="structure_derivation", derivation=derivation)
    return mid, diag


M = {}
for key in ("R_S_3aa", "rac_1e", "2o", "R_L3", "S_L3", "SSp_L6"):
    M[key] = copy_existing(key, [ev(7)])

# Scheme 8 complex structures: retain labels/configurations without guessed SMILES.
M["allenyl_acid"] = molecule("allenyl acid derived from (R,S)-3aa", [ev(7)], "allenyl_acid_RS_3aa",
                              assignment="R,S", axial="R",
                              centers=[{"source_order": 2, "configuration": "S"}])
M["SS_8"] = molecule("(S,S)-8", [ev(7)], "S_S_8", assignment="S,S",
                      centers=[{"source_order": 1, "configuration": "S"},
                               {"source_order": 2, "configuration": "S"}],
                      geometry="mixture",
                      geometry_raw=">20:1 Z/E; major alkene geometry is Z")
molecules[M["SS_8"]]["stereochemistry"]["double_bonds"] = F(
    ">20:1 Z/E", {"major": "Z", "minor": "E", "ratio": ">20:1", "isomer_order": "Z:E"}, [ev(7)])
M["3ba_generic"] = molecule("3ba (configuration not specified at Scheme 8 arrow c)", [ev(7)], "scheme8_3ba")
M["3ga_generic"] = molecule("3ga (configuration not specified at Scheme 8 arrow c)", [ev(7)], "scheme8_3ga")
M["9"] = molecule("9 (R = Me)", [ev(7)], "9", note="Scheme 8 intermediate 9 (R = Me); configuration is not printed beside the central structure.")
M["10"] = molecule("10 [R = -(CH2)4-]", [ev(7)], "10", note="Scheme 8 intermediate 10 [R = -(CH2)4-]; configuration is not printed beside the central structure.")

for num in (11, 12, 13, 14):
    pages = [ev(7), ev(8)] if num in (13, 14) else [ev(7)]
    for cfg, key in (("R,S", f"R_S_{num}"), ("S,S", f"S_S_{num}")):
        M[key] = molecule(f"({cfg})-{num}", pages, key, assignment=cfg,
                          centers=[{"source_order": 1, "configuration": cfg.split(',')[0]},
                                   {"source_order": 2, "configuration": "S"}],
                          geometry="not_applicable")
for cfg, key in (("R,S", "R_S_3eo"), ("S,S", "S_S_3eo")):
    M[key] = molecule(f"({cfg})-3eo", [ev(7), ev(8)], key, assignment=cfg,
                      axial=cfg.split(",")[0],
                      centers=[{"source_order": 2, "configuration": "S"}],
                      geometry="unresolved",
                      geometry_raw="allene axial configuration is specified, but other drawn double-bond geometry was not independently assigned")

stereo_diagnostics = {}
for key, label, smiles in (
    ("KOH", "KOH", "[K+].[OH-]"), ("EtOH", "EtOH", "CCO"),
    ("I2", "I2", "II"), ("hexane", "hexane", "CCCCCC"),
    ("TBAF", "TBAF", "CCCC[N+](CCCC)(CCCC)CCCC.[F-]"),
    ("THF_s8", "THF", "C1CCOC1"),
    ("AgOTf", "AgOTf", "[Ag+].O=S(=O)([O-])C(F)(F)F"),
    ("CH2Cl2", "CH2Cl2", "ClCCl"), ("PdCl2", "PdCl2", "Cl[Pd]Cl"),
    ("allyl_bromide", "allylic bromide", "C=CCBr"),
    ("DMA", "DMA", "CC(=O)N(C)C"),
    ("citric_acid", "citric acid", "O=C(O)CC(O)(CC(=O)O)C(=O)O"),
    ("TsCl", "TsCl", "Cc1ccc(S(=O)(=O)Cl)cc1"),
    ("pyridine", "pyridine", "n1ccccc1"),
):
    M[key], stereo_diagnostics[key] = reagent(label, smiles, [ev(7)], key=key)
M["Ph3PAuCl"] = molecule("Ph3PAuCl", [ev(7)], "Ph3PAuCl", geometry="unresolved",
                          note="Gold(I) catalyst is source-identified by formula/name; no complete coordination-complex SMILES is asserted.")
M["citric_10pct"] = molecule("10% citric acid solution (basis unreported)", [ev(7)], "citric_acid_10pct",
                              geometry="not_applicable",
                              note="The figure reports a 10% citric acid solution but does not state the percentage basis; the solution is retained as a source-labelled reagent.")


def reaction(rid, pages, name, category="derivatization", transform=None):
    evidence = [ev(p) for p in pages]
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = F([category], [category], evidence)
    r["location"] = F({"pdf_pages": pages, "scheme": 8, "arrow": rid},
                      {"pdf_pages": pages, "scheme": 8, "arrow": rid}, evidence)
    r["reaction_name"] = F(name, name, evidence)
    r["transformation_type"] = F(
        transform or name, transform or name, evidence,
        method="structure_derivation",
        derivation="Transformation type derived from the displayed Scheme 8 arrow and product class.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": [
                       "Supporting Information was mentioned but not supplied or explicitly associated; detailed procedure, quantities, times, workup and purification are unavailable.",
                       "Complete complex-substrate/product connectivity and stereochemical SMILES remain unresolved from the dense main-text drawing."]}
    r["completeness"] = "partial"
    return r


def addp(r, mid, role, evidence, stage="reaction", product_role=None):
    pid = f"p{len(r['participants'])+1}_{safe(mid.split(':')[-1])}"
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    r["participants"].append(p)
    return pid


def step(r, sid, evidence, participant_ids, operation, temperature=None, stage="reaction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = participant_ids
    s["operation"] = F(operation, operation, evidence)
    if temperature:
        s["temperature"] = F(temperature, temperature, evidence)
    s["time"] = F(evidence=evidence, status="unreported")
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, target, raw, normalized, evidence, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target], evidence)
    o["value"] = F(raw, normalized, evidence)
    o["product_role"] = F("reported product", "reported product", evidence)
    if scope is not None:
        o["scope"] = F(scope, scope, evidence)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    # Scheme 8 prints values but not the yield/selectivity analytical basis.
    o["measurement_basis"] = F(evidence=evidence, status="unreported")
    r["outcomes"].append(o)


def set_bonds(r, evidence, formed, reason):
    changes, labels = [], {}
    for i, (a, b, kind) in enumerate(formed, 1):
        cid = f"b{i}"
        change = {"change_id": cid,
                  "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": kind}
        changes.append(F(None, change, evidence, method="structure_derivation", derivation=reason))
        lab = st.bond_label(a, b, 1)
        labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(cid)
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F(None, list(labels.values()), evidence,
                                          method="structure_derivation", derivation=reason)
    r["classification"]["coverage"] = F(None, "partial", evidence,
                                            method="inference",
                                            derivation="Only visually evident newly formed bonds are classified; full atom mapping, broken bonds, bond-order changes and hydrogen balance are unresolved.")
    r["classification"]["limitations"] = [
        "Complete atom-mapped reaction SMILES are unavailable from the dense figure.",
        "Only visually evident newly formed bonds are asserted; implicit hydrogen and leaving-group balance remain unresolved."]


def add_selectivity_set(r, target, evidence, dr=None, ee=None, ze=None):
    if dr:
        qual = "gt" if dr.startswith(">") else None
        major, minor = dr.lstrip(">").split(":")
        outcome(r, "dr", "selectivity", target, f"{dr} dr",
                {"metric": "dr", "major": float(major), "minor": float(minor), "qualifier": qual,
                 "isomer_order": "displayed product : other diastereomer"}, evidence,
                scope={"metric": "dr", "target_entity_ids": [target],
                       "isomer_order": "displayed product : other diastereomer"})
    if ee:
        qual = "gt" if ee.startswith(">") else None
        value = float(ee.lstrip(">"))
        outcome(r, "ee", "selectivity", target, f"{ee}% ee",
                q(value, "%", qual, f"{ee}% ee"), evidence,
                scope={"metric": "ee", "target_entity_ids": [target],
                       "configuration": molecules.get(target, existing.get(target))["stereochemistry"]["assignment"]["normalized_value"]})
    if ze:
        qual = "gt" if ze.startswith(">") else None
        major, minor = ze.lstrip(">").split(":")
        outcome(r, "Z_E", "selectivity", target, f"{ze} Z/E",
                {"metric": "Z/E", "Z": float(major), "E": float(minor), "qualifier": qual,
                 "isomer_order": "Z:E"}, evidence,
                scope={"metric": "Z/E", "target_entity_ids": [target], "isomer_order": "Z:E"})


E7 = [ev(7)]

# (R,S)-3aa --a,b--> (S,S)-8. Preserve the reported two-step yield basis.
r = reaction("scheme8_RS_3aa_to_SS_8", [7], "hydrolysis and electrophilic iodocyclization",
             transform="two-step hydrolysis/iodocyclization with axial-to-central chirality transfer")
p_start = addp(r, M["R_S_3aa"], "reactant", E7)
p_koh = addp(r, M["KOH"], "reagent", E7)
p_etoh = addp(r, M["EtOH"], "solvent", E7)
p_int = addp(r, M["allenyl_acid"], "product", E7, product_role="step-a intermediate")
p_i2 = addp(r, M["I2"], "reagent", E7)
p_hex = addp(r, M["hexane"], "solvent", E7)
p_prod = addp(r, M["SS_8"], "product", E7, product_role="reported final product")
step(r, "a_hydrolysis", E7, [p_start, p_koh, p_etoh, p_int], "hydrolyzed with KOH in EtOH", "reflux")
step(r, "b_iodocyclization", E7, [p_int, p_i2, p_hex, p_prod], "electrophilic iodocyclization with I2 in hexane", "rt")
outcome(r, "yield", "yield", M["SS_8"], "47% yield for two steps", q(47, "%", raw_text="47% yield for two steps"), E7,
        yield_type="two-step overall yield from (R,S)-3aa to (S,S)-8")
add_selectivity_set(r, M["SS_8"], E7, dr=">20:1", ee=">99", ze=">20:1")
r["notes"].append(F("The structure of (S,S)-8 was confirmed by X-ray crystallography.",
                    "source reports X-ray confirmation of (S,S)-8", E7))
set_bonds(r, E7, [("C", "O", "formed"), ("C", "I", "formed")],
          "The displayed iodocyclization forms a lactone C-O bond and a C-I bond; the hydrolysis/iodocyclization sequence is not fully atom mapped.")
bundle["reactions"].append(r)

# Scheme 8 central c arrow has an explicit '3ba or 3ga' split and two products.
for source_key, product_key, rid in (
    ("3ba_generic", "9", "scheme8_3ba_to_9"),
    ("3ga_generic", "10", "scheme8_3ga_to_10"),
):
    r = reaction(rid, [7], "TBAF deprotection", transform="silyl ether deprotection")
    ps = addp(r, M[source_key], "reactant", E7)
    pt = addp(r, M["TBAF"], "reagent", E7)
    pthf = addp(r, M["THF_s8"], "solvent", E7)
    pp = addp(r, M[product_key], "product", E7, product_role="displayed deprotected intermediate")
    step(r, "c_deprotection", E7, [ps, pt, pthf, pp], "TBAF deprotection in THF", "rt")
    set_bonds(r, E7, [("O", "H", "hydrogen_gain")],
              "TBAF deprotection gives the displayed alcohol; the O-H bond is the newly formed bond and O-Si cleavage is not a formation label.")
    bundle["reactions"].append(r)

# Gold(I)-catalyzed rearrangements (four separately displayed arrows).
for start_key, prod_key, yld, dr, ee in (
    ("9", "R_S_11", 61, ">20:1", ">99"), ("9", "S_S_11", 43, "12:1", "99"),
    ("10", "R_S_12", 83, ">20:1", ">99"), ("10", "S_S_12", 40, "18:1", ">99"),
):
    rid = f"scheme8_{start_key}_to_{prod_key}"
    r = reaction(rid, [7], "gold(I)-catalyzed allene rearrangement",
                 transform="gold(I)-catalyzed cycloisomerization to dihydrofuran")
    ps = addp(r, M[start_key], "reactant", E7)
    pau = addp(r, M["Ph3PAuCl"], "catalyst", E7)
    pag = addp(r, M["AgOTf"], "reagent", E7)
    psol = addp(r, M["CH2Cl2"], "solvent", E7)
    pp = addp(r, M[prod_key], "product", E7, product_role="reported product")
    step(r, "d_gold_cyclization", E7, [ps, pau, pag, psol, pp], "Ph3PAuCl/AgOTf cyclization in CH2Cl2", "rt")
    outcome(r, "yield", "yield", M[prod_key], f"{yld}%", q(yld, "%"), E7)
    add_selectivity_set(r, M[prod_key], E7, dr=dr, ee=ee)
    set_bonds(r, E7, [("C", "O", "formed")],
              "The displayed allene-to-dihydrofuran rearrangement creates the ring C-O bond; other bond-order changes are unresolved.")
    bundle["reactions"].append(r)

# Pd-catalyzed cyclizative coupling of 9 (two separately displayed arrows).
for prod_key, yld, dr, ee in (
    ("R_S_13", 93, "13:1", "99"), ("S_S_13", 82, "5:1", "99"),
):
    rid = f"scheme8_9_to_{prod_key}"
    r = reaction(rid, [7, 8], "Pd-catalyzed cyclizative coupling",
                 transform="Pd-catalyzed allenol cyclization/allylation to dihydrofuran")
    evidence = [ev(7), ev(8)]
    ps = addp(r, M["9"], "reactant", evidence)
    ppd = addp(r, M["PdCl2"], "catalyst", evidence)
    pab = addp(r, M["allyl_bromide"], "reactant", evidence)
    psol = addp(r, M["DMA"], "solvent", evidence)
    pp = addp(r, M[prod_key], "product", evidence, product_role="reported product")
    step(r, "e_palladium_cyclization", evidence, [ps, ppd, pab, psol, pp],
         "PdCl2-catalyzed cyclizative coupling with allylic bromide in DMA", "rt")
    outcome(r, "yield", "yield", M[prod_key], f"{yld}%", q(yld, "%"), evidence)
    add_selectivity_set(r, M[prod_key], evidence, dr=dr, ee=ee)
    set_bonds(r, evidence, [("C", "O", "formed"), ("C", "C", "formed")],
              "The displayed cyclizative coupling forms the dihydrofuran C-O bond and attaches the allyl fragment through a new C-C bond.")
    bundle["reactions"].append(r)

# Two ligand-controlled rac-1e + 2o scope arrows. No unprinted common conditions are inherited.
for lig_key, prod_key, yld, dr in (
    ("R_L3", "R_S_3eo", 66, "17:1"), ("S_L3", "S_S_3eo", 65, ">20:1"),
):
    rid = f"scheme8_rac1e_2o_to_{prod_key}"
    r = reaction(rid, [7], "stereodivergent allenylation", category="substrate_scope",
                 transform="ligand-controlled allenylic C-C bond formation")
    p1 = addp(r, M["rac_1e"], "reactant", E7)
    p2 = addp(r, M["2o"], "reactant", E7)
    pl1 = addp(r, M[lig_key], "ligand", E7)
    pl2 = addp(r, M["SSp_L6"], "ligand", E7)
    pp = addp(r, M[prod_key], "product", E7, product_role="reported product")
    step(r, "allenylation", E7, [p1, p2, pl1, pl2, pp],
         "rac-1e and 2o reacted using the two displayed chiral ligands; other conditions are not printed on this arrow")
    outcome(r, "yield", "yield", M[prod_key], f"{yld}% yield", q(yld, "%"), E7)
    add_selectivity_set(r, M[prod_key], E7, dr=dr, ee="98")
    r["review"]["unresolved"].append(
        "Scheme 8 prints only the two chiral ligands on this arrow; catalyst precursors, base, solvent, loadings, temperature and time were not inherited from earlier tables.")
    set_bonds(r, E7, [("C", "C", "formed")],
              "The displayed allenylation joins the carbon nucleophile derived from 2o to the allenylic fragment from rac-1e.")
    bundle["reactions"].append(r)

# Four explicitly ordered operations f,g,c,e to each 14 stereoisomer.
for start_key, prod_key, yld, ee in (
    ("R_S_3eo", "R_S_14", 80, "99"), ("S_S_3eo", "S_S_14", 71, "98"),
):
    rid = f"scheme8_{start_key}_to_{prod_key}"
    r = reaction(rid, [7, 8], "four-operation conversion to dihydrofuran amino-acid derivative",
                 transform="deprotection, tosylation, desilylation and Pd-catalyzed cyclization")
    evidence = [ev(7), ev(8)]
    ps = addp(r, M[start_key], "reactant", evidence)
    pcit = addp(r, M["citric_10pct"], "reagent", evidence)
    pthf = addp(r, M["THF_s8"], "solvent", evidence)
    pts = addp(r, M["TsCl"], "reagent", evidence)
    ppyr = addp(r, M["pyridine"], "reagent", evidence)
    ptbaf = addp(r, M["TBAF"], "reagent", evidence)
    ppd = addp(r, M["PdCl2"], "catalyst", evidence)
    pab = addp(r, M["allyl_bromide"], "reactant", evidence)
    pdma = addp(r, M["DMA"], "solvent", evidence)
    pp = addp(r, M[prod_key], "product", evidence, product_role="reported final product")
    step(r, "f_citric_acid", evidence, [ps, pcit, pthf],
         "treated with 10% citric acid (percentage basis unreported) in THF", "rt")
    step(r, "g_tosylation", evidence, [pts, ppyr], "treated with TsCl in pyridine", "rt")
    step(r, "c_desilylation", evidence, [ptbaf, pthf], "treated with TBAF in THF", "rt")
    step(r, "e_palladium_cyclization", evidence, [ppd, pab, pdma, pp],
         "PdCl2-catalyzed cyclizative coupling with allylic bromide in DMA", "rt")
    outcome(r, "yield_step_e", "yield", M[prod_key], f"{yld}% yield for step e", q(yld, "%", raw_text=f"{yld}% yield for step e"), evidence,
            yield_type="yield for final operation e only; not an overall four-operation yield")
    add_selectivity_set(r, M[prod_key], evidence, dr=">20:1", ee=ee)
    r["notes"].append(F(
        f"{yld}% yield is explicitly labelled 'yield for step e'; it is not treated as the overall yield from {start_key}.",
        {"reported_yield_percent": yld, "yield_scope": "step e only", "overall_yield": None}, evidence))
    set_bonds(r, evidence, [("N", "S", "formed"), ("N", "H", "hydrogen_gain"),
                            ("C", "O", "formed"), ("C", "C", "formed")],
              "Across the displayed f,g,c,e sequence, tosylation, deprotection and cyclizative coupling form N-S, N-H, ring C-O and allyl C-C bonds; intermediates are not drawn and full atom provenance is unresolved.")
    bundle["reactions"].append(r)

bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id": f"{DOC}:page7_scheme7_dft", "document_id": DOC, "evidence_ids": [ev(7)],
     "kind": "theoretical_calculation", "note": "Scheme 7 presents four calculated nucleophilic-substitution pathways, not four performed experiments."},
    {"context_id": f"{DOC}:page7_figure2_dft", "document_id": DOC, "evidence_ids": [ev(7)],
     "kind": "theoretical_calculation", "note": "Figure 2 reports calculated relative transition-state energies 0.0, 10.27, 9.84 and 13.64 kcal/mol; it is not a synthetic experiment."},
    {"context_id": f"{DOC}:page8_conclusions_si", "document_id": DOC, "evidence_ids": [ev(8)],
     "kind": "article_context", "note": "Page 8 completes the Scheme 8 discussion, gives conclusions, and mentions SI. No SI file is supplied or explicitly associated in source_mapping.json."},
    {"context_id": f"{DOC}:page9_references", "document_id": DOC, "evidence_ids": [ev(9)],
     "kind": "references", "note": "Bibliographic references only; cited literature was not recorded as a present-paper experiment."},
    {"context_id": f"{DOC}:page10_references", "document_id": DOC, "evidence_ids": [ev(10)],
     "kind": "references", "note": "Bibliographic references only; cited literature was not recorded as a present-paper experiment."},
    {"context_id": f"{DOC}:page11_references_note", "document_id": DOC, "evidence_ids": [ev(11)],
     "kind": "references_and_correction", "note": "References plus Note Added After ASAP Publication: an incorrect Scheme 1 structure in the August 5 version was corrected and reposted August 9. This is source context, not a new experiment."},
]

st.write(SUBMITTED, bundle)
st.write(ROOT / "runs/current/raw_model/pages7_11_stereo_diagnostics.json", stereo_diagnostics)
st.write(ROOT / "runs/current/raw_model/pages7_11_source_mentions.json",
         st.source_mentions(dataset, ["3aa", "3ba", "3ga", "3eo", "14"], context=260))

inventory = {
    "document_id": DOC,
    "scope": "original PDF pages 7-11",
    "pages_visually_opened": [7, 8, 9, 10, 11],
    "page_inventory": {
        "7": ["Scheme 7 DFT pathways", "Figure 2 DFT energies", "Scheme 8 reaction arrows and footnote conditions"],
        "8": ["continuation of Scheme 8 discussion", "conclusions", "SI mention", "author information", "references begin"],
        "9": ["references only"], "10": ["references only"],
        "11": ["references", "Note Added After ASAP Publication correction"],
    },
    "scheme8_arrow_count": 13,
    "notes": [
        "The two c branches (3ba->9 and 3ga->10) were retained separately because the figure explicitly says '3ba or 3ga'.",
        "The two f,g,c,e sequences retain the printed operation order and final-step-only yield scope.",
        "All records remain needs_review; complex structures were not assigned guessed SMILES.",
    ],
}
st.write(ROOT / "runs/current/raw_model/pages7_11_visual_transcription.json", inventory)

by_page = {p: [] for p in range(7, 12)}
for rxn in bundle["reactions"]:
    for p in by_page:
        if ev(p) in rxn["evidence_ids"]:
            by_page[p].append(rxn["reaction_id"])

regions = {
    7: [
        {"region": "upper-left Scheme 7", "finding": "Four DFT pathways; theoretical calculation, not four performed reactions.", "evidence_ids": [ev(7)]},
        {"region": "upper-right Figure 2", "finding": "Calculated transition states and relative energies; retained as source_context, not synthesis records.", "evidence_ids": [ev(7)]},
        {"region": "lower-right Scheme 8 and footnote", "finding": "Thirteen displayed transformation branches were inventoried: one two-step 3aa->8 sequence, two c deprotections, four d cyclizations, two e cyclizations, two rac-1e/2o scope arrows, and two f,g,c,e sequences to 14.", "evidence_ids": [ev(7)]},
    ],
    8: [
        {"region": "top continuation of Results and Discussion", "finding": "Completes the two product-13 dr values and confirms 3eo-to-14 axial-to-central chirality transfer; evidence merged into the four related records.", "evidence_ids": [ev(8)]},
        {"region": "Conclusions and Associated Content", "finding": "Summary and SI mention only; no additional independent reaction arrow. SI was not supplied or explicitly associated.", "evidence_ids": [ev(8)]},
        {"region": "Author Information/Acknowledgments/references start", "finding": "Metadata and bibliography, not present-paper experimental records.", "evidence_ids": [ev(8)]},
    ],
    9: [{"region": "full page reference list", "finding": "Cited literature only; no present-paper experimental arrow or table.", "evidence_ids": [ev(9)]}],
    10: [{"region": "full page reference list", "finding": "Cited literature only; no present-paper experimental arrow or table.", "evidence_ids": [ev(10)]}],
    11: [
        {"region": "reference list", "finding": "Cited literature only; no present-paper experimental arrow or table.", "evidence_ids": [ev(11)]},
        {"region": "Note Added After ASAP Publication", "finding": "Correction notice for Scheme 1 structure (reposted August 9, 2021); source context, not an experiment.", "evidence_ids": [ev(11)]},
    ],
}
unresolved = {
    7: ["SI procedures and reagent quantities/times are unavailable.", "Complex structure connectivity and most complete stereochemical SMILES remain unresolved.", "Yield and selectivity measurement bases are not printed in Scheme 8."],
    8: ["The continuation prose provides no additional quantities or analytical methods; SI remains unassociated."],
    9: ["No reaction records expected: page contains references only."],
    10: ["No reaction records expected: page contains references only."],
    11: ["No reaction records expected: page contains references/correction note only; the corrected Scheme 1 lies outside this page range."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False,
            "pages": [{"page": p, "state": "attempted", "regions": regions[p],
                       "record_ids": by_page[p], "unresolved": unresolved[p]}
                      for p in range(7, 12)]}
st.write(ROOT / "runs/current/coverage/pages_0007_0011.json", coverage)

print(json.dumps({"molecules_in_incremental_bundle": len(bundle["molecules"]),
                  "new_reactions": len(bundle["reactions"]),
                  "reaction_ids": [r["reaction_id"] for r in bundle["reactions"]]},
                 ensure_ascii=False, indent=2))
