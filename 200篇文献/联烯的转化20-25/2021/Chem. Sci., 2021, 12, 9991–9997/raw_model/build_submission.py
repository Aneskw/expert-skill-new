#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for DOI 10.1039/d1sc02681j.

This is a session extraction artifact, not a chemistry tool.  It deliberately
keeps compact structures unresolved when the supplied six pages are
insufficient for a defensible concrete structure.
"""
from copy import deepcopy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_68d6021490685b69"
RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
SUBMITTED = RUN / "submitted.json"
TITLE = "Catalytic asymmetric synthesis of spirocyclobutyl oxindoles and beyond via [2+2] cycloaddition and sequential transformations"


def ev(page, kind="text"):
    return [f"{DOC}:page{page}:{kind}"]


def F(raw=None, normalized=None, page=None, method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    evidence = [] if page is None else (page if isinstance(page, list) else ev(page))
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, verification=verification, alternatives=alternatives)


def qty(number, unit, raw=None):
    return {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}


dataset = json.loads(DATASET.read_text())
dataset["documents"][0]["title"] = TITLE
DATASET.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")
bundle = st.empty(dataset["documents"])
molecules = {}


def ensure_mol(mid, label, page, *, name=None, chirality="unresolved",
               geometry="unresolved", representation="unresolved", assignment=None,
               mixture=None, smiles=None, identity_note=None):
    evidence = ev(page)
    if mid in molecules:
        m = molecules[mid]
        if not any(o["evidence_ids"] == evidence for o in m["occurrences"]):
            m["occurrences"].append({"document_id": DOC,
                "compound_label": F(label, label, page), "evidence_ids": evidence})
        return m
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC,
        "compound_label": F(label, label, page), "evidence_ids": evidence}]
    m["names"]["original"] = F(name or label, name or label, page)
    m["names"]["english"] = F(name, name, page) if name else F()
    m["names"]["abbreviations"] = F(label, [label], page)
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["stereochemistry"]["chirality_status"] = F(
        chirality, chirality, page, method="inference",
        derivation="Status preserves the source drawing/assignment; unresolved or unspecified is not treated as achiral or racemic.")
    m["stereochemistry"]["geometry_status"] = F(
        geometry, geometry, page, method="inference",
        derivation="Status preserves source E/Z information and does not infer geometry from omission.")
    m["stereochemistry"]["representation_status"] = F(
        representation, representation, page, method="inference",
        derivation="Compared with the supplied source drawing; compact structures are not forced into SMILES.")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = F(assignment, assignment, page)
    if mixture is not None:
        m["stereochemistry"]["mixture"] = F(mixture, mixture, page)
    if identity_note:
        m["identity_note"] = F(identity_note, identity_note, page)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(label, smiles, page, method="structure_derivation",
            derivation="Candidate SMILES transcribed from the fully visible source drawing; RDKit parsing/stereo diagnostics do not certify identity.")
        m["stereochemistry"]["isomeric_smiles"] = F(label, diag["canonical_isomeric_smiles"], page,
            method="structure_derivation", derivation="RDKit canonical isomeric SMILES from the evidence-linked candidate structure.")
        m["stereochemistry"]["identity_key"] = F(label, diag["stereo_graph_key"], page,
            method="structure_derivation", derivation="RDKit stereo graph key from the evidence-linked candidate structure.")
    molecules[mid] = m
    return m


# Core substrates. Only the two fully legible model substrates receive candidate SMILES.
ensure_mol("m_1a", "1a", 2, name="(E)-alkenyloxindole 1a", chirality="achiral",
    geometry="specified", representation="encoded",
    assignment={"alkene": "E (source designation)"},
    smiles="CC(C)(C)OC(=O)N1C(=O)/C(=C/C(=O)OC(C)(C)C)c2ccccc21",
    identity_note="Concrete candidate from the page-2 drawing; no CAS number is supplied.")
ensure_mol("m_2a", "2a", 2, name="N-allenamide 2a", chirality="achiral",
    geometry="not_applicable", representation="partial",
    assignment={"allene_axis": "terminal CH2 allene; no R_a/S_a assignment printed"},
    smiles="C=C=CN(c1ccccc1)S(=O)(=O)c1ccc(C)cc1",
    identity_note="Candidate connectivity from the page-2 drawing; toolkit reports an unassigned cumulene stereo element, so representation is partial.")


def product3(label, page=3):
    assignment = None
    chirality = "partially_specified"
    if label == "3aa":
        assignment = {"absolute_configuration": "(1R,2R,4S)", "method": "X-ray crystallography"}
        chirality = "specified"
    else:
        assignment = {"source_drawing": "single drawn diastereomer", "method": "assigned by CD comparison with 3aa"}
    return ensure_mol("m_" + label, label, page, name=f"spirocyclobutyl oxindole {label}",
        chirality=chirality, geometry="not_applicable", representation="not_encoded",
        assignment=assignment, identity_note="Compact source structure retained by compound label; no concrete SMILES asserted.")


def product4(label, page=3):
    chirality, geometry, assignment = "partially_specified", "unspecified", None
    if label == "4aa":
        geometry = "specified"
        assignment = {"alkene_geometry": "E", "basis": "E-4aa named in mechanistic discussion"}
    if label == "4ka":
        chirality, geometry = "specified", "specified"
        assignment = {"absolute_configuration": "R", "alkene_geometry": "E",
                      "basis": "assigned from 3ka absolute configuration and rac-4ka relative configuration"}
    if assignment is None:
        assignment = {"source_drawing": "enantioenriched product drawn with one stereochemical bond",
                      "absolute_configuration": "not explicitly named on supplied pages"}
    return ensure_mol("m_" + label, label, page, name=f"dihydropyran-fused indole {label}",
        chirality=chirality, geometry=geometry, representation="not_encoded",
        assignment=assignment, identity_note="Compact source structure retained by label; unspecified stereochemical elements are not filled by analogy.")


def product5(label, page=4):
    return ensure_mol("m_" + label, label, page, name=f"gamma-lactone substituted oxindole {label}",
        chirality="mixture", geometry="not_applicable", representation="not_encoded",
        mixture={"type": "diastereomeric mixture", "ratio": "experiment-specific; see reaction outcome"},
        identity_note="Product is reported as two diastereomers; absolute configurations are not printed on supplied pages.")


def ensure_agent(mid, label, page, name=None):
    return ensure_mol(mid, label, page, name=name or label, chirality="unresolved",
                      geometry="not_applicable", representation="not_encoded")


# Shared products and agents.
product3("3aa", 2); product4("4aa", 2); product5("5aa", 3)
ensure_mol("m_control_byproduct", "unidentified by-product", 3,
           name="unidentified by-product in 3aa to 4aa control",
           chirality="unresolved", geometry="unresolved", representation="unresolved",
           identity_note="Only described as a tiny amount; no structure or compound number is supplied on the main-text pages.")
for mid, label in [
    ("a_NiOTf2", "Ni(OTf)2"), ("a_CuOTf2", "Cu(OTf)2"),
    ("a_ZnOTf2", "Zn(OTf)2"), ("a_MgOTf2", "Mg(OTf)2"),
    ("a_CH2Cl2", "CH2Cl2"), ("a_CHCl3", "CHCl3"), ("a_DCE", "DCE"),
    ("a_EtOAc", "EtOAc"), ("a_MS4A", "4 Å molecular sieves"),
    ("a_TFA", "TFA"), ("a_H2O", "H2O"), ("a_PdC", "Pd/C"), ("a_H2", "H2")]:
    ensure_agent(mid, label, 2 if label in {"Ni(OTf)2", "Cu(OTf)2", "Zn(OTf)2", "Mg(OTf)2", "CH2Cl2", "CHCl3", "4 Å molecular sieves"} else 4)
for lig in ["L3-PiPr2", "L2-PiPr2", "L4-PiPr2", "L4-PrPr2", "L4-RaPr2", "L5-RaPr2", "L1", "L2"]:
    ensure_agent("a_" + lig.replace("-", "_"), lig, 2)


def participant(pid, mid, role, page, *, amount=None, mass=None, volume=None,
                equiv=None, molpct=None, reference=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev(page), stage=stage)
    if product_role is not None:
        p["product_role"] = F(product_role, product_role, page)
    values = {"amount": amount, "mass": mass, "volume": volume,
              "equivalents": equiv, "mol_percent": molpct}
    for key, value in values.items():
        if value is not None:
            raw = value.get("raw_text", f"{value['number']} {value['unit']}")
            p["charges"][key] = F(raw, value, page)
    if (equiv is not None or molpct is not None) and reference:
        p["equivalence_reference"] = F(reference, reference, page)
    return p


def step(step_id, page, pids, operation, temp=None, time=None, atmosphere=None,
         pressure=None, stage="reaction", temp_fact=None):
    s = st.new_step(step_id, ev(page), stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, page)
    if temp_fact is not None:
        s["temperature"] = temp_fact
    elif temp is not None:
        if isinstance(temp, str):
            s["temperature"] = F(temp, temp, page)
        else:
            s["temperature"] = F(temp.get("raw_text"), temp, page)
    if time is not None:
        s["time"] = F(time.get("raw_text"), time, page)
    if atmosphere is not None:
        s["atmosphere"] = F(atmosphere, atmosphere, page)
    if pressure is not None:
        if isinstance(pressure, str):
            s["pressure"] = F(pressure, pressure, page)
        else:
            s["pressure"] = F(pressure.get("raw_text"), pressure, page)
    return s


def outcome(oid, kind, targets, page, raw, normalized, *, basis=None,
            yield_type=None, scope=None, product_role="product"):
    o = st.new_outcome(oid, kind, targets, ev(page))
    o["product_role"] = F(product_role, product_role, page)
    o["value"] = F(raw, normalized, page)
    if yield_type is not None:
        o["yield_type"] = F(yield_type, yield_type, page)
    if basis is not None:
        o["measurement_basis"] = F(basis, basis, page)
    if scope is not None:
        o["scope"] = F(scope, scope, page)
    return o


TAX = {"C-H": 1, "C-C": 2, "C-O": 3, "N-H": 29}


def set_classification(r, page, codes, limitation):
    if codes is None:
        r["classification"]["coverage"] = F("unresolved", "unresolved", page,
            method="inference", derivation="Compact structures or missing atom correspondence prevent defensible net bond enumeration.")
        r["classification"]["labels"] = F("unresolved", [], page,
            method="structure_derivation", derivation="No bond-formation label asserted without a defensible structural comparison.")
        r["classification"]["limitations"] = [limitation]
        return
    if not codes:
        r["classification"]["coverage"] = F("not applicable", "not_applicable", page,
            method="inference", derivation="The record is context or a proposed mechanism rather than an observed net transformation.")
        r["classification"]["labels"] = F("none", [], page,
            method="structure_derivation", derivation="No observed net transformation to classify.")
        r["classification"]["limitations"] = [limitation]
        return
    changes, label_map = [], {}
    for i, code in enumerate(codes, 1):
        a, b = code.split("-") if code != "C-H" else ("C", "H")
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                  {"atom_map": None, "element": b}], "before_order": 0,
                  "after_order": 1, "change_type": "hydrogen_gain" if b == "H" else "formed"}
        if b == "H":
            change["hydrogen_count"] = 1
        changes.append(F(f"source-structure-derived {code} formation", change, page,
            method="structure_derivation", derivation="Partial net bond change from the source reaction drawing; atoms are not fully mapped."))
        label_map.setdefault(code, {"category_number": TAX[code], "code": code,
                                   "supporting_change_ids": []})["supporting_change_ids"].append(cid)
    labels = list(label_map.values())
    r["classification"]["bond_changes"] = changes
    r["classification"]["labels"] = F("partial source-structure-derived labels", labels, page,
        method="structure_derivation", derivation="Labels are derived from the explicitly drawn net transformation, without claiming complete atom mapping.")
    r["classification"]["coverage"] = F("partial", "partial", page,
        method="inference", derivation="Only unambiguous newly formed bonds are represented; hydrogen/coproduct and full atom correspondence may remain unresolved.")
    r["classification"]["limitations"] = [limitation, "No atom-mapped reaction SMILES is asserted."]


def base_reaction(suffix, page, location, record_type="specific_experiment", categories=None,
                  name=None, transformation=None):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, ev(page))
    r["record_type"] = record_type
    r["experiment_categories"] = F(str(categories or []), categories or [], page)
    r["location"] = F(location, {"page": page, "location": location}, page)
    if name:
        r["reaction_name"] = F(name, name, page)
    if transformation:
        r["transformation_type"] = F(transformation, transformation, page,
            method="inference", derivation="Normalized transformation description from the source reaction name/drawing.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_solvent(r, page, composition, concentration=None, volume=None):
    r["solvent_conditions"]["composition"] = F(composition, composition, page)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = F(concentration.get("raw_text"), concentration, page)
        r["solvent_conditions"]["concentration_basis"] = F("based on limiting substrate", "limiting substrate", page)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = F(volume.get("raw_text"), volume, page)


reactions = []

# Page 1: distinguish cited literature from the two this-work generic transformations.
for suffix, loc, name in [
    ("scheme1a_gold_prior", "Scheme 1a, upper literature route", "reported prior Au(I)-catalyzed [2+2] cycloaddition"),
    ("scheme1a_lewis_prior", "Scheme 1a, lower literature route", "reported prior Lewis-acid-catalyzed [2+2] cycloaddition")]:
    r = base_reaction(suffix, 1, loc, "literature_example", [], name, "literature context")
    set_classification(r, 1, [], "Cited prior work is not counted as a current-paper experiment.")
    r["review"]["unresolved"] = ["Generic literature structures and conditions are not fully specified on the supplied page."]
    reactions.append(r)

for suffix, loc, product, codes, name in [
    ("scheme1b_generic_2plus2", "Scheme 1b, first arrow", "3", ["C-C", "C-C"], "asymmetric [2+2] cycloaddition"),
    ("scheme1b_generic_sequential", "Scheme 1b, second arrow", "4", ["C-C", "C-O"], "one-pot sequential transformation")]:
    ensure_mol(f"m_generic_{product}", product, 1, chirality="unresolved", geometry="unresolved", representation="not_encoded")
    r = base_reaction(suffix, 1, loc, "generic_scheme", [], name, "generic reaction scheme")
    r["participants"] = [participant("p_product", f"m_generic_{product}", "product", 1, product_role="generic product")]
    set_classification(r, 1, codes, "Generic R-group scheme; classification is partial and no concrete atom mapping is possible.")
    r["review"]["unresolved"] = ["Generic R groups and full condition variants are not concrete experiments."]
    reactions.append(r)


# Table 1, 11 independent optimization experiments.
table1 = [
 (1, None, 13, "1:99", None, 0), (2, "L3-PiPr2", 58, "89:11", 91, 51),
 (3, "L3-PiPr2", 39, "1:99", None, 60), (4, "L2-PiPr2", 59, "60:40", 58, 83),
 (5, "L4-PiPr2", 55, "97:3", 95, None), (6, "L4-PrPr2", 55, "98:2", 79, None),
 (7, "L4-RaPr2", 75, "97:3", 90, None), (8, "L5-RaPr2", 74, "98:2", 95, None),
 (9, "L5-RaPr2", 94, "98:2", 95, None), (10, "L1", 29, "1:99", None, 0),
 (11, "L2", 33, "1:99", None, -11)]
for entry, lig, yld, ratio, ee3, ee4 in table1:
    r = base_reaction(f"table1_entry{entry:02d}", 2, f"Table 1 entry {entry}", categories=["optimization"],
                      name="optimization of asymmetric [2+2] cycloaddition", transformation="[2+2] cycloaddition with competing rearrangement")
    eq2, conc, temp = (1.2, 0.2, 35) if entry == 9 else (1.0, 0.1, 30)
    p = [participant("p_1a", "m_1a", "reactant", 2, amount=qty(0.05, "mmol")),
         participant("p_2a", "m_2a", "reactant", 2, equiv=qty(eq2, "equiv"), reference="m_1a"),
         participant("p_ni", "a_NiOTf2", "catalyst", 2, molpct=qty(15 if entry == 3 else 10, "mol%"), reference="m_1a"),
         participant("p_solvent", "a_CH2Cl2", "solvent", 2),
         participant("p_3aa", "m_3aa", "product", 2, product_role="[2+2] product"),
         participant("p_4aa", "m_4aa", "product", 2, product_role="rearranged product")]
    if lig:
        p.insert(3, participant("p_ligand", "a_" + lig.replace("-", "_"), "ligand", 2,
                               molpct=qty(10, "mol%"), reference="m_1a"))
    r["participants"] = p
    r["reaction_scale"] = F("1a (0.05 mmol)", {"reference_molecule_id": "m_1a", "amount": qty(0.05, "mmol")}, 2)
    set_solvent(r, 2, "CH2Cl2", qty(conc, "M"))
    r["steps"] = [step("reaction", 2, [x["participant_id"] for x in p], "stir", qty(temp, "°C"), qty(16, "h"), "N2")]
    r["outcomes"] = [outcome("yield_total", "yield", ["m_3aa", "m_4aa"], 2, f"{yld}%", qty(yld, "%"), basis="isolated", yield_type="combined isolated yield", scope="3aa + 4aa"),
                     outcome("ratio", "selectivity", ["m_3aa", "m_4aa"], 2, ratio, ratio, basis="HPLC on a chiral stationary phase", scope="3aa:4aa")]
    if ee3 is not None:
        r["outcomes"].append(outcome("ee_3aa", "selectivity", ["m_3aa"], 2, f"{ee3}% ee", qty(ee3, "%"), basis="HPLC on a chiral stationary phase", scope="ee of 3aa"))
    if ee4 is not None:
        ee4_norm = qty(ee4, "%") if ee4 >= 0 else "−11% ee (signed value; opposite major enantiomer as printed)"
        r["outcomes"].append(outcome("ee_4aa", "selectivity", ["m_4aa"], 2, f"{ee4}% ee", ee4_norm, basis="HPLC on a chiral stationary phase", scope="ee of 4aa; negative sign preserved where printed"))
    set_classification(r, 2, ["C-C", "C-C", "C-O"], "Combined-product optimization row; labels cover unambiguous formations across the two alternative products, not a single complete atom mapping.")
    r["review"]["unresolved"] = ["Combined isolated yield is not apportioned between 3aa and 4aa."]
    reactions.append(r)


# Table 2, six sequential-transformation optimization experiments.
table2 = [(1,"Ni(OTf)2",76,"82:18",92),(2,"Zn(OTf)2",54,"6:94",91),
          (3,"Mg(OTf)2",52,"34:66",91),(4,"Cu(OTf)2",49,"1:99",93),
          (5,"Cu(OTf)2",70,"1:99",93),(6,"Cu(OTf)2",70,"1:99",93)]
metal_ids = {"Ni(OTf)2":"a_NiOTf2","Zn(OTf)2":"a_ZnOTf2","Mg(OTf)2":"a_MgOTf2","Cu(OTf)2":"a_CuOTf2"}
for entry, metal, yld, ratio, ee in table2:
    r = base_reaction(f"table2_entry{entry:02d}", 2, f"Table 2 entry {entry}", categories=["optimization"],
        name="Lewis-acid optimization for sequential rearrangement", transformation="[2+2] cycloaddition followed by rearrangement")
    if entry == 6:
        ps = [participant("p_3aa", "m_3aa", "reactant", 2),
              participant("p_cu", "a_CuOTf2", "catalyst", 2, molpct=qty(5,"mol%"), reference="m_3aa"),
              participant("p_ms", "a_MS4A", "additive", 2, mass=qty(10,"mg")),
              participant("p_solvent", "a_CHCl3", "solvent", 2),
              participant("p_4aa", "m_4aa", "product", 2, product_role="major product")]
        r["reaction_scale"] = F()
        set_solvent(r,2,"CHCl3",qty(0.05,"M"))
        r["steps"] = [step("rearrangement",2,[x["participant_id"] for x in ps],"treat isolated 3aa with Cu(OTf)2 and 4 Å MS",qty(30,"°C"),qty(4,"h"),"N2")]
    else:
        ps = [participant("p_1a","m_1a","reactant",2,amount=qty(0.05,"mmol")),
              participant("p_2a","m_2a","reactant",2,equiv=qty(1.2,"equiv"),reference="m_1a"),
              participant("p_ligand","a_L5_RaPr2","ligand",2,molpct=qty(10,"mol%"),reference="m_1a"),
              participant("p_ni","a_NiOTf2","catalyst",2,molpct=qty(10,"mol%"),reference="m_1a"),
              participant("p_metal2",metal_ids[metal],"catalyst",2,molpct=qty(5,"mol%"),reference="m_1a"),
              participant("p_3aa","m_3aa","product",2,product_role="minor/residual product"),
              participant("p_4aa","m_4aa","product",2,product_role="rearranged product")]
        if entry == 5:
            ps += [participant("p_ms","a_MS4A","additive",2,mass=qty(10,"mg")), participant("p_solvent2","a_CHCl3","solvent",2)]
        else:
            ps += [participant("p_solvent1","a_CH2Cl2","solvent",2)]
        r["reaction_scale"] = F("1a (0.05 mmol)",{"reference_molecule_id":"m_1a","amount":qty(0.05,"mmol")},2)
        r["steps"] = [step("cycloaddition",2,["p_1a","p_2a","p_ligand","p_ni"],"initial asymmetric cycloaddition",qty(35,"°C"),qty(16,"h"),None),
                      step("rearrangement",2,[x["participant_id"] for x in ps],"add second metal salt" + (" and 4 Å MS" if entry==5 else ""),qty(30,"°C"),qty(3 if entry==5 else 16,"h"),"air")]
        set_solvent(r,2,"CH2Cl2 then " + ("CHCl3" if entry==5 else "CH2Cl2"),qty(0.05 if entry==5 else 0.1,"M"))
    r["participants"] = ps
    r["outcomes"] = [outcome("yield_4aa","yield",["m_4aa"],2,f"{yld}%",qty(yld,"%"),basis="isolated",yield_type="isolated yield",scope="based on 1a"),
        outcome("ratio","selectivity",["m_3aa","m_4aa"],2,ratio,ratio,basis="HPLC on a chiral stationary phase",scope="3aa:4aa"),
        outcome("ee_4aa","selectivity",["m_4aa"],2,f"{ee}% ee",qty(ee,"%"),basis="HPLC on a chiral stationary phase",scope="ee of 4aa")]
    set_classification(r,2,["C-O"],"Partial net classification emphasizes cyclobutane opening/rearrangement to the C–O ring; bond breaking and full atom correspondence remain unresolved.")
    r["review"]["unresolved"] = ["Entry-specific workup and purification are in unavailable ESI."]
    reactions.append(r)


# Table 3 scope data.
scope3 = {
"3aa":(93,95),"3ba":(85,70),"3ca":(86,78),"3da":(93,93),"3ea":(93,81),
"3fa":(91,93),"3ga":(94,94),"3ha":(77,92),"3ia":(72,88),"3ja":(81,88),"3ka":(83,90),
"3la":(84,93),"3ma":(90,93),"3na":(81,93),"3oa":(86,87),"3pa":(83,96),"3qa":(76,92),
"3ab":(94,95),"3ac":(88,95),"3ad":(85,95),"3ae":(94,94),"3af":(69,96),"3ag":(97,93)}
scope4 = {
"4aa":(70,93),"4ba":(73,70),"4ca":(83,78),"4da":(71,93),"4ea":(86,83),
"4fa":(72,91),"4ga":(79,92),"4ha":(67,92),"4ia":(55,89),"4ja":(61,88),"4ka":(57,90),
"4la":(31,91),"4ma":(62,92),"4na":(54,90),"4oa":(57,88),"4pa":(61,94),"4qa":(57,90),
"4ab":(71,92),"4ac":(63,93),"4ad":(75,91),"4ae":(83,90),"4af":(40,93),"4ag":(71,93)}


def substrate_ids(label):
    core = label[1:]
    x, y = core[0], core[1]
    m1, m2 = f"m_1{x}", f"m_2{y}"
    ensure_mol(m1, f"1{x}", 3, name=f"(E)-alkenyloxindole 1{x}", chirality="unresolved", geometry="specified", representation="not_encoded", assignment={"alkene":"E (family designation)"})
    ensure_mol(m2, f"2{y}", 3, name=f"N-allenamide 2{y}", chirality="unspecified", geometry="not_applicable", representation="not_encoded")
    return m1, m2


temp_conflict = F(["cond. 1: CH2Cl2, 35 °C, 16 h", "footnote: same as Table 1 entry 8 (30 °C)"], None, 3)
temp_conflict["data_status"] = "conflict"
temp_conflict["alternatives"] = [F("35 °C",qty(35,"°C"),3), F("Table 1 entry 8: 30 °C",qty(30,"°C"),3)]
for label,(yld,eev) in scope3.items():
    p3 = product3(label); m1,m2 = substrate_ids(label)
    r = base_reaction(f"table3_{label}",3,f"Table 3 product {label}",categories=["substrate_scope"],name="asymmetric [2+2] cycloaddition",transformation="[2+2] cycloaddition")
    r["general_scheme_ids"]=[f"{DOC}:scheme1b_generic_2plus2"]
    ps=[participant("p_1",m1,"reactant",3,amount=qty(0.05,"mmol")),participant("p_2",m2,"reactant",3,equiv=qty(1.0,"equiv"),reference=m1),
        participant("p_ligand","a_L5_RaPr2","ligand",3,molpct=qty(10,"mol%"),reference=m1),participant("p_ni","a_NiOTf2","catalyst",3,molpct=qty(10,"mol%"),reference=m1),
        participant("p_solvent","a_CH2Cl2","solvent",3),participant("p_product","m_"+label,"product",3,product_role="major [2+2] product")]
    r["participants"]=ps;r["reaction_scale"]=F("same as Table 1 entry 8: 1 (0.05 mmol)",{"reference_molecule_id":m1,"amount":qty(0.05,"mmol")},3,method="inherited_general_condition")
    set_solvent(r,3,"CH2Cl2",qty(0.1,"M"))
    r["steps"]=[step("cycloaddition",3,[x["participant_id"] for x in ps],"stir under cond. 1",time=qty(16,"h"),atmosphere="N2",temp_fact=deepcopy(temp_conflict))]
    r["outcomes"]=[outcome("yield","yield",["m_"+label],3,f"{yld}%",qty(yld,"%"),basis="isolated",yield_type="isolated yield",scope="based on substrate 1"),
        outcome("ee","selectivity",["m_"+label],3,f"{eev}% ee",qty(eev,"%"),basis="HPLC on a chiral stationary phase",scope="ee; major enantiomer configuration follows drawing/assignment where given")]
    set_classification(r,3,["C-C","C-C"],"Two new C–C connections are inferred from the [2+2] drawing; complete atom mapping and all bond-order changes are unresolved.")
    r["notes"]=[F("Source conflict: cond. 1 prints 35 °C, while footnote says conditions same as Table 1 entry 8 (30 °C).", "temperature conflict retained without choosing",3)]
    r["review"]["unresolved"]=["Temperature conflict between Table 3 cond. 1 and inherited Table 1 entry 8.","Individual workup/purification are unavailable in the supplied main text."]
    reactions.append(r)

for label,(yld,eev) in scope4.items():
    p4=product4(label);m1,m2=substrate_ids(label); code=label[1:]
    cu=10 if label=="4ea" else 7.5 if label in {"4ab","4ac","4ad","4ae","4af"} else 12.5 if label=="4ag" else 5
    r=base_reaction(f"table3_{label}",3,f"Table 3 product {label}",categories=["substrate_scope"],name="one-pot sequential transformation",transformation="[2+2] cycloaddition followed by Cu-catalyzed rearrangement")
    r["general_scheme_ids"]=[f"{DOC}:scheme1b_generic_sequential"]
    ps=[participant("p_1",m1,"reactant",3,amount=qty(0.05,"mmol")),participant("p_2",m2,"reactant",3,equiv=qty(1.2,"equiv"),reference=m1),
        participant("p_ligand","a_L5_RaPr2","ligand",3,molpct=qty(10,"mol%"),reference=m1),participant("p_ni","a_NiOTf2","catalyst",3,molpct=qty(10,"mol%"),reference=m1),
        participant("p_cu","a_CuOTf2","catalyst",3,molpct=qty(cu,"mol%"),reference=m1),participant("p_ms","a_MS4A","additive",3,mass=qty(10,"mg")),
        participant("p_solvent1","a_CH2Cl2","solvent",3),participant("p_solvent2","a_CHCl3","solvent",3),participant("p_product","m_"+label,"product",3,product_role="rearranged product")]
    side=None
    if label in {"4fa","4la"}:
        side="3"+code; product3(side); ps.append(participant("p_side","m_"+side,"product",3,product_role="residual [2+2] product"))
    r["participants"]=ps;r["reaction_scale"]=F("same as Table 2 entry 5: 1 (0.05 mmol)",{"reference_molecule_id":m1,"amount":qty(0.05,"mmol")},3,method="inherited_general_condition")
    set_solvent(r,3,"CH2Cl2 (first stage), then CHCl3",qty(0.05,"M"))
    r["steps"]=[step("cycloaddition",3,["p_1","p_2","p_ligand","p_ni","p_solvent1"],"initial asymmetric cycloaddition",qty(35,"°C"),qty(16,"h"),None),
        step("rearrangement",3,[x["participant_id"] for x in ps],"add Cu(OTf)2 and 4 Å MS in CHCl3",qty(30,"°C"),qty(3,"h"),"air")]
    r["outcomes"]=[outcome("yield","yield",["m_"+label],3,f"{yld}%",qty(yld,"%"),basis="isolated",yield_type="isolated yield",scope="based on substrate 1"),
        outcome("ee","selectivity",["m_"+label],3,f"{eev}% ee",qty(eev,"%"),basis="HPLC on a chiral stationary phase",scope="ee")]
    if side:
        ratio="19:1" if label=="4fa" else "9:1"
        r["outcomes"].append(outcome("product_ratio","selectivity",["m_"+label,"m_"+side],3,ratio,ratio,basis=None,scope=f"{label}:{side}"))
    set_classification(r,3,["C-C","C-O"],"Formal [4+2]/sequential product drawing supports partial C–C and C–O formation labels; full atom mapping is not asserted.")
    r["review"]["unresolved"]=["Individual workup/purification are unavailable in the supplied main text."]
    reactions.append(r)


# Page 3 controls reported in prose.
r=base_reaction("control_3aa_to_4aa_Ni",3,"page 3 control experiment prose",categories=["control","mechanistic_study"],name="Ni(OTf)2-promoted rearrangement control",transformation="cyclobutane rearrangement")
r["participants"]=[participant("p_3aa","m_3aa","reactant",3),participant("p_ni","a_NiOTf2","catalyst",3,molpct=qty(5,"mol%"),reference="m_3aa"),participant("p_4aa","m_4aa","product",3,product_role="rearranged product"),participant("p_byproduct","m_control_byproduct","product",3,product_role="unidentified trace by-product")]
r["steps"]=[step("reaction",3,[x["participant_id"] for x in r["participants"]],"treat isolated 3aa with Ni(OTf)2",time=qty(16,"h"),atmosphere="air")]
r["outcomes"]=[outcome("yield","yield",["m_4aa"],3,"37%",qty(37,"%"),basis=None,yield_type=None,scope="4aa; accompanied by a tiny amount of by-product"),outcome("byproduct","trace",["m_control_byproduct"],3,"tiny amount of by-product","tiny amount",basis=None,scope="unidentified by-product")]
set_classification(r,3,["C-O"],"Partial rearrangement classification; unidentified by-product and full bond changes remain unresolved.")
r["review"]["unresolved"]=["Solvent, temperature, yield measurement basis, and by-product identity require unavailable ESI Table S5/Fig. S4."]
reactions.append(r)

r=base_reaction("control_4aa_reverse_Ni",3,"page 3 reverse control prose",categories=["control","mechanistic_study"],name="reverse rearrangement control",transformation="no observed conversion to 3aa")
r["participants"]=[participant("p_4aa","m_4aa","reactant",3),participant("p_ni","a_NiOTf2","catalyst",3,molpct=qty(5,"mol%"),reference="m_4aa")]
r["steps"]=[step("reaction",3,["p_4aa","p_ni"],"subject 4aa to the same conditions as forward control",time=qty(16,"h"),atmosphere="air")]
r["outcomes"]=[outcome("no_transform","no_reaction",["m_4aa"],3,"did not transform into 3aa at all","no transformation observed",basis=None,scope="attempted 4aa to 3aa reverse transformation",product_role="reactant")]
set_classification(r,3,[],"Negative control has no observed net transformation; expected product is not fabricated.")
r["review"]["unresolved"]=["Solvent, temperature, and observation method require unavailable ESI."]
reactions.append(r)


# Scheme 2 transformations.
r=base_reaction("scheme2a_4aa_TFA_to_5aa",4,"Scheme 2a, left arrow",categories=["derivatization"],name="TFA/water conversion to gamma-lactone",transformation="deprotection, ring opening and lactonization")
r["evidence_ids"].append(f"{DOC}:page3:text")
r["participants"]=[participant("p_4aa","m_4aa","reactant",4),participant("p_tfa","a_TFA","reagent",4),participant("p_water","a_H2O","reagent",4),participant("p_5aa","m_5aa","product",4,product_role="gamma-lactone product")]
r["steps"]=[step("reaction",4,[x["participant_id"] for x in r["participants"]],"treat 4aa with TFA and H2O")]
r["outcomes"]=[outcome("yield","yield",["m_5aa"],4,"59% yield",qty(59,"%"),basis=None,scope="4aa to 5aa"),outcome("dr","selectivity",["m_5aa"],3,"65:35 dr","65:35",basis=None,scope="5aa diastereomer ratio; value supplied by accompanying prose"),outcome("ee","selectivity",["m_5aa"],3,"95% ee for each diastereomer",qty(95,"%"),basis=None,scope="each 5aa diastereomer; configuration order unassigned")]
set_classification(r,4,None,"Hydrolysis/deprotection/rearrangement is structurally complex and complete atom correspondence is unavailable.")
r["review"]["unresolved"]=["Temperature, time, reagent quantities, workup, yield basis, and selectivity methods are not stated on supplied pages."]
reactions.append(r)

scheme2b=[("3aa","5aa",70,"65:35","95/95"),("3ga","5ga",56,"69:31","94/94"),("3ia","5ia",52,"68:32","88/88"),("3na","5na",44,"67:33","93/93"),("3ab","5ab",60,"66:34","95/95")]
for sm,pm,yld,dr,eev in scheme2b:
    product3(sm);product5(pm)
    r=base_reaction(f"scheme2b_{sm}_to_{pm}",4,f"Scheme 2b {sm} to {pm}",categories=["derivatization"],name="Cu(OTf)2-promoted gamma-lactone formation",transformation="deprotection, hydrolysis and lactonization")
    if sm == "3aa":
        r["evidence_ids"].append(f"{DOC}:page3:text")
    r["participants"]=[participant("p_start","m_"+sm,"reactant",4),participant("p_cu","a_CuOTf2","catalyst",4,molpct=qty(10,"mol%"),reference="m_"+sm),participant("p_solvent","a_CH2Cl2","solvent",4),participant("p_product","m_"+pm,"product",4,product_role="gamma-lactone product")]
    set_solvent(r,4,"CH2Cl2")
    r["steps"]=[step("reaction",4,[x["participant_id"] for x in r["participants"]],"treat under air",qty(35,"°C"),atmosphere="air")]
    r["outcomes"]=[outcome("yield","yield",["m_"+pm],4,f"{yld}%",qty(yld,"%"),basis=None,scope=f"{sm} to {pm}"),outcome("dr","selectivity",["m_"+pm],4,dr,dr,basis=None,scope=f"{pm} diastereomer ratio"),outcome("ee","selectivity",["m_"+pm],4,f"{eev}% ee",eev,basis=None,scope="ee values in major/minor diastereomer order; absolute configurations unassigned")]
    set_classification(r,4,["C-O","N-H"],"Partial net lactonization/deprotection classification; hydrolysis partners, broken bonds and full H balance are unresolved.")
    r["review"]["unresolved"]=["Reaction time, workup/purification, yield basis, and selectivity methods are not stated on supplied pages."]
    reactions.append(r)

for code,R,yld,eev in [("aa","C6H5",60,97),("ab","4-MeC6H4",50,94),("ac","4-ClC6H4",48,95)]:
    ensure_mol("m_6"+code,"6"+code,4,name=f"alkoxyallene 6{code} (R={R})",chirality="unspecified",geometry="not_applicable",representation="not_encoded")
    ensure_mol("m_7"+code,"7"+code,4,name=f"dihydropyran-fused indole 7{code}",chirality="unspecified",geometry="unspecified",representation="not_encoded")
    r=base_reaction(f"scheme2c_7{code}",4,f"Scheme 2c product 7{code}",categories=["substrate_scope"],name="asymmetric formal [4+2] cycloaddition of alkoxyallene",transformation="formal [4+2] cycloaddition")
    r["participants"]=[participant("p_1a","m_1a","reactant",4),participant("p_6","m_6"+code,"reactant",4),participant("p_ligand","a_L5_RaPr2","ligand",4),participant("p_ni","a_NiOTf2","catalyst",4),participant("p_solvent","a_CH2Cl2","solvent",4),participant("p_product","m_7"+code,"product",4,product_role="formal [4+2] product")]
    set_solvent(r,4,"DCM")
    r["steps"]=[step("reaction",4,[x["participant_id"] for x in r["participants"]],"react under L5-RaPr2/Ni(OTf)2 catalysis",qty(35,"°C"))]
    r["outcomes"]=[outcome("yield","yield",["m_7"+code],4,f"{yld}% yield",qty(yld,"%"),basis=None,scope=f"product 7{code}"),outcome("ee","selectivity",["m_7"+code],4,f"{eev}% ee",qty(eev,"%"),basis=None,scope="major enantiomer configuration not assigned")]
    set_classification(r,4,["C-C","C-O"],"Partial formal [4+2] classification from the colored source drawing; atom mapping and all bond-order changes are unresolved.")
    r["review"]["unresolved"]=["Substrate amounts, catalyst loadings, time, atmosphere, workup, yield basis, and ee method require unavailable ESI Table S10."]
    reactions.append(r)


# Scheme 3 scale-ups and derivatizations.
r=base_reaction("scheme3a_scale_3aa",4,"Scheme 3a",categories=["scale_up"],name="scale-up asymmetric [2+2] cycloaddition",transformation="[2+2] cycloaddition")
r["participants"]=[participant("p_1a","m_1a","reactant",4,amount=qty(2.0,"mmol")),participant("p_2a","m_2a","reactant",4,amount=qty(2.4,"mmol")),participant("p_ligand","a_L5_RaPr2","ligand",4,molpct=qty(10,"mol%"),reference="m_1a"),participant("p_ni","a_NiOTf2","catalyst",4,molpct=qty(10,"mol%"),reference="m_1a"),participant("p_solvent","a_CH2Cl2","solvent",4),participant("p_3aa","m_3aa","product",4,product_role="scale-up product")]
r["reaction_scale"]=F("1a, 2.0 mmol",{"reference_molecule_id":"m_1a","amount":qty(2.0,"mmol")},4);set_solvent(r,4,"CH2Cl2",qty(0.2,"M"))
r["steps"]=[step("reaction",4,[x["participant_id"] for x in r["participants"]],"standard optimized conditions",qty(35,"°C"),qty(16,"h"),"N2")]
r["outcomes"]=[outcome("yield","yield",["m_3aa"],4,"89% yield",qty(89,"%"),basis=None,scope="scale-up 1a to 3aa"),outcome("amount","product_amount",["m_3aa"],4,"1.12 g",qty(1.12,"g"),basis=None,scope="isolated product amount not explicitly labelled as isolated"),outcome("ee","selectivity",["m_3aa"],4,"90% ee",qty(90,"%"),basis=None,scope="3aa")]
set_classification(r,4,["C-C","C-C"],"Partial [2+2] bond formation; complete atom mapping is not asserted.")
r["review"]["unresolved"]=["Scheme/prose does not state yield or ee measurement methods; workup/purification unavailable."]
reactions.append(r)

product4("4qa",4);ensure_mol("m_1q","1q",4,name="(E)-alkenyloxindole 1q",chirality="unresolved",geometry="specified",representation="not_encoded",assignment={"alkene":"E family substrate"})
r=base_reaction("scheme3b_scale_4qa",4,"Scheme 3b",categories=["scale_up"],name="scale-up two-step synthesis of 4qa",transformation="[2+2] cycloaddition and Cu rearrangement")
r["participants"]=[participant("p_1q","m_1q","reactant",4,amount=qty(1.60,"mmol")),participant("p_2a","m_2a","reactant",4,amount=qty(1.96,"mmol")),participant("p_ligand","a_L5_RaPr2","ligand",4,molpct=qty(10,"mol%"),reference="m_1q"),participant("p_ni","a_NiOTf2","catalyst",4,molpct=qty(10,"mol%"),reference="m_1q"),participant("p_cu","a_CuOTf2","catalyst",4,molpct=qty(7.5,"mol%"),reference="m_1q"),participant("p_ms","a_MS4A","additive",4,mass=qty(400,"mg")),participant("p_solvent1","a_CH2Cl2","solvent",4),participant("p_solvent2","a_CHCl3","solvent",4),participant("p_product","m_4qa","product",4,product_role="two-step product")]
r["reaction_scale"]=F("1q, 1.60 mmol",{"reference_molecule_id":"m_1q","amount":qty(1.60,"mmol")},4)
r["steps"]=[step("cycloaddition",4,["p_1q","p_2a","p_ligand","p_ni","p_solvent1"],"standard conditions",qty(35,"°C"),qty(16,"h"),"N2"),step("rearrangement",4,[x["participant_id"] for x in r["participants"]],"add Cu(OTf)2 and 4 Å MS in CHCl3 (0.07 M)",qty(30,"°C"),qty(3,"h"),"air")]
set_solvent(r,4,"CH2Cl2 then CHCl3",qty(0.07,"M"));r["outcomes"]=[outcome("yield","yield",["m_4qa"],4,"58% yield for two steps",qty(58,"%"),basis=None,scope="overall from 1q to 4qa over two steps"),outcome("amount","product_amount",["m_4qa"],4,"0.62 g",qty(0.62,"g"),basis=None,scope="4qa"),outcome("ee","selectivity",["m_4qa"],4,"88% ee",qty(88,"%"),basis=None,scope="4qa")]
set_classification(r,4,["C-C","C-O"],"Partial overall two-step net classification; overall yield is not assigned to the last step alone.")
r["review"]["unresolved"]=["Yield/ee methods and detailed workup are not stated; second-stage time is inherited from Table 2 entry 5."]
reactions.append(r)

product3("3ab",4);product5("5ab",4)
r=base_reaction("scheme3c_scale_5ab",4,"Scheme 3c",categories=["scale_up","derivatization"],name="scale-up gamma-lactone formation",transformation="deprotection, hydrolysis and lactonization")
r["participants"]=[participant("p_3ab","m_3ab","reactant",4,amount=qty(2.5,"mmol")),participant("p_cu","a_CuOTf2","catalyst",4,molpct=qty(10,"mol%"),reference="m_3ab"),participant("p_solvent","a_CH2Cl2","solvent",4),participant("p_product","m_5ab","product",4,product_role="gamma-lactone product")]
r["reaction_scale"]=F("3ab, 2.5 mmol",{"reference_molecule_id":"m_3ab","amount":qty(2.5,"mmol")},4);set_solvent(r,4,"CH2Cl2",qty(0.1,"M"));r["steps"]=[step("reaction",4,[x["participant_id"] for x in r["participants"]],"treat under air",qty(35,"°C"),atmosphere="air")]
r["outcomes"]=[outcome("yield","yield",["m_5ab"],4,"63% yield",qty(63,"%"),basis=None,scope="3ab to 5ab"),outcome("amount","product_amount",["m_5ab"],4,"0.72 g",qty(0.72,"g"),basis=None,scope="5ab"),outcome("dr","selectivity",["m_5ab"],4,"66:34 dr","66:34",basis=None,scope="5ab diastereomer ratio"),outcome("ee","selectivity",["m_5ab"],4,"94/94% ee","94/94%",basis=None,scope="major/minor diastereomers; configurations unassigned")]
set_classification(r,4,["C-O","N-H"],"Partial lactonization/deprotection classification; full atom and H balance unresolved.")
r["review"]["unresolved"]=["Reaction time, outcome methods, and workup/purification are not stated."]
reactions.append(r)

ensure_mol("m_8","8",4,name="hydrogenation product 8",chirality="mixture",geometry="not_applicable",representation="not_encoded",mixture={"type":"diastereomeric mixture","ratio":"85:15"})
ensure_mol("m_9","9",4,name="hydrogenation/ring-cleavage product 9",chirality="mixture",geometry="not_applicable",representation="not_encoded",mixture={"type":"diastereomeric mixture","ratio":"60:40"})
for code,start,prod,pdwt,yld,dr,eev in [("d","3aa","8",10,95,"85:15","95"),("e","4aa","9",100,75,"60:40","92/92")]:
    r=base_reaction(f"scheme3{code}_hydrogenation_{prod}",4,f"Scheme 3{code}",categories=["derivatization"],name="hydrogenation",transformation="hydrogenation" + (" with C–O bond cleavage" if code=="e" else ""))
    r["participants"]=[participant("p_start","m_"+start,"reactant",4),participant("p_pdc","a_PdC","catalyst",4),participant("p_h2","a_H2","reagent",4),participant("p_solvent","a_EtOAc","solvent",4),participant("p_product","m_"+prod,"product",4,product_role="hydrogenation product")]
    set_solvent(r,4,"EtOAc",qty(0.05,"M"));r["steps"]=[step("hydrogenation",4,[x["participant_id"] for x in r["participants"]],"hydrogenate over Pd/C",temp="r.t.",pressure="5 MPa")]
    r["outcomes"]=[outcome("yield","yield",["m_"+prod],4,f"{yld}% yield",qty(yld,"%"),basis=None,scope=f"{start} to {prod}"),outcome("dr","selectivity",["m_"+prod],4,f"{dr} dr",dr,basis=None,scope=f"product {prod} diastereomer ratio"),outcome("ee","selectivity",["m_"+prod],4,f"{eev}% ee",eev if "/" in eev else qty(int(eev),"%"),basis=None,scope="ee values in printed diastereomer order; configurations unassigned")]
    r["notes"]=[F(f"Pd/C ({pdwt}% wt)", {"agent":"Pd/C","loading":{"number":pdwt,"unit":"wt%"}}, 4)]
    set_classification(r,4,["C-H","C-H"],"Two net C–H gains across the reduced C=C are retained; reduced C=C does not create a new C–C sigma label. Additional hydrogen/proton transfers associated with C–O cleavage in Scheme 3e and full atom mapping remain unresolved.")
    r["review"]["unresolved"]=["Reaction time, substrate scale, outcome methods, and workup/purification are not stated."]
    reactions.append(r)


# Thermal control mentioned on page 4.
r=base_reaction("thermal_3aa_to_E4aa",4,"page 4 mechanistic discussion",categories=["control","mechanistic_study"],name="thermal rearrangement control",transformation="thermal cyclobutane rearrangement")
r["participants"]=[participant("p_3aa","m_3aa","reactant",4),participant("p_dce","a_DCE","solvent",4),participant("p_4aa","m_4aa","product",4,product_role="E-4aa")]
set_solvent(r,4,"DCE");r["steps"]=[step("heating",4,[x["participant_id"] for x in r["participants"]],"heat 3aa in DCE",qty(80,"°C"),qty(24,"h"))]
r["outcomes"]=[outcome("formation","other",["m_4aa"],4,"led to the generation of E-4aa","E-4aa generated",basis=None,scope="qualitative; no yield reported")]
set_classification(r,4,["C-O"],"Qualitative thermal rearrangement; no yield or complete atom mapping is supplied.")
r["review"]["unresolved"]=["Yield/conversion and observation method are not reported on supplied pages."]
reactions.append(r)


# Page 5 proposed mechanism is not an observed experiment.
r=base_reaction("figure2_proposed_mechanism",5,"Figure 2",record_type="proposed_mechanism",categories=[],name="proposed transition-state model and catalytic cycles",transformation="proposed elementary pathway")
r["participants"]=[participant("p_1a","m_1a","reactant",5),participant("p_2a","m_2a","reactant",5),participant("p_3aa","m_3aa","product",5),participant("p_4aa","m_4aa","product",5),participant("p_5aa","m_5aa","product",5),participant("p_cu","a_CuOTf2","catalyst",5)]
set_classification(r,5,[],"Arrow-pushing intermediates and DFT pathway are proposed, not independent observed net experiments.")
r["review"]["unresolved"]=["Proposed intermediates I–VI are not registered as isolated molecules; detailed computations are in unavailable ESI."]
reactions.append(r)


# Source context: non-synthetic observations and data availability, not reaction rows.
bundle["source_context"] = [
 {"context_id":f"{DOC}:abstract","document_id":DOC,"kind":"article_summary","evidence_ids":ev(1),"text":"Abstract describes asymmetric [2+2] cycloaddition and sequential Cu transformations."},
 {"context_id":f"{DOC}:figure1_dft","document_id":DOC,"kind":"computational_energy_profile","evidence_ids":ev(5),"text":"DFT energy profile for Cu(OTf)2-catalyzed 3aa to E-/Z-4aa; not a separate synthesis experiment."},
 {"context_id":f"{DOC}:data_availability","document_id":DOC,"kind":"supplement_notice","evidence_ids":ev(6),"text":"Paper states procedures, spectra, computational details and crystallographic data are in ESI, but no SI file was provided or associated in this session."},
 {"context_id":f"{DOC}:xray","document_id":DOC,"kind":"characterization","evidence_ids":ev(3),"text":"X-ray supports (1R,2R,4S)-3aa; rac-4ka and rac-5aa are cited. Characterization is not converted into a synthesis record."}
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
SUBMITTED.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

diagnostics = {
    "1a_candidate": st.stereo_inspect("CC(C)(C)OC(=O)N1C(=O)/C(=C/C(=O)OC(C)(C)C)c2ccccc21"),
    "2a_candidate": st.stereo_inspect("C=C=CN(c1ccccc1)S(=O)(=O)c1ccc(C)cc1"),
    "note": "Graph diagnostics only; compact products remain not_encoded/unresolved rather than guessed."
}
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")

inventory = {"document_id": DOC, "pages_reviewed": [1,2,3,4,5,6],
 "visual_inventory": {
  "1": ["title/abstract", "Scheme 1a two prior-work generic routes", "Scheme 1b two this-work generic transformations"],
  "2": ["Table 1 rows 1-11", "Table 2 rows 1-6", "footnotes a-e"],
  "3": ["Table 3 products 3aa-3ag (23 arrows)", "Table 3 products 4aa-4ag (23 arrows)", "forward and reverse rearrangement controls in prose"],
  "4": ["Scheme 2a two arrows", "Scheme 2b five arrows", "Scheme 2c three arrows", "Scheme 3a-e five arrows", "thermal control in prose"],
  "5": ["Figure 1 computational energy profile", "Figure 2 proposed transition-state model/mechanism", "no independent synthetic arrow counted"],
  "6": ["conclusion/data availability/authorship/references", "no current-work synthetic experiment"]},
 "reaction_records": len(reactions), "molecules": len(molecules), "semantic_accuracy_verified": False}
(RUN / "raw_model/page1_6_visual_inventory.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"reactions":len(reactions),"molecules":len(molecules)},ensure_ascii=False))
