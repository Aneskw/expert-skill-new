#!/usr/bin/env python3
"""Extend the cumulative v3.1 extraction with original PDF pages 13--18.

Pages 13--16 provide product-specific preparations under General Procedures B
and C.  Pages 17--18 provide detailed procedures for experiments already
identified on earlier pages; those identities are reused rather than duplicated.
All machine records remain needs_review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills/chemical-reaction-extraction"
sys.path.insert(0, str(SKILL / "scripts"))

from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    write,
)

DOC = "doc_c32350973ac9b683"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 22)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 22)}

dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
bundle = copy.deepcopy(dataset["standard"])
bundle["documents"] = copy.deepcopy(dataset["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def q(number, unit, qualifier=None, raw_text=None, **extra):
    out = {"number": number, "unit": unit}
    if qualifier:
        out["qualifier"] = qualifier
    if raw_text:
        out["raw_text"] = raw_text
    out.update(extra)
    return out


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None, verification=None, alternatives=()):
    return fact(raw, normalized, list(evidence), method, status, verification,
                derivation, list(alternatives))


def setf(container, key, raw, normalized, evidence, method="original_extraction",
         derivation=None, status=None, verification=None, alternatives=()):
    container[key] = F(raw, normalized, evidence, method, derivation, status,
                       verification, alternatives)


def union(values, more):
    return list(dict.fromkeys(list(values) + list(more)))


def union_fact_evidence(value, evidence):
    if isinstance(value, dict) and "evidence_ids" in value:
        value["evidence_ids"] = union(value["evidence_ids"], evidence)


def conflict_fact(alternatives, evidence, derivation):
    return F(None, None, evidence, "original_extraction", derivation,
             status="conflict", verification="unreviewed",
             alternatives=alternatives)


def add_generic_molecule(mid, label, evidence, note_text):
    ev = list(dict.fromkeys(evidence))
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    setf(m["occurrences"][0], "compound_label", label, label, ev)
    setf(m["names"], "original", label, label, ev)
    setf(m["names"], "english", label, label, ev)
    setf(m["names"], "abbreviations", None, None, ev,
         status="unreported", verification="missing")
    setf(m["names"], "aliases", None, None, ev,
         status="unreported", verification="missing")
    setf(m, "cas_numbers", None, None, ev,
         status="unreported", verification="missing")
    m["representation_type"] = "unresolved"
    setf(m, "smiles", None, None, ev, status="unreported", verification="missing")
    st = m["stereochemistry"]
    setf(st, "chirality_status", "unresolved", "unresolved", ev)
    setf(st, "geometry_status", "unresolved", "unresolved", ev)
    for key in ("atom_centers", "double_bonds", "axial", "assignment"):
        setf(st, key, None, None, ev, status="unreported", verification="missing")
    setf(st, "mixture", None, None, ev, status="not_applicable", verification="missing")
    setf(st, "representation_status", "not_encoded", "not_encoded", ev,
         method="inference", derivation="Generic variable structure is not a concrete stereochemical graph.",
         verification="inferred")
    setf(st, "isomeric_smiles", None, None, ev,
         status="unreported", verification="missing")
    setf(st, "identity_key", None, None, ev,
         status="unreported", verification="missing")
    setf(m, "identity_note", note_text, note_text, ev, method="inference",
         derivation="Evidence-linked generic identity note.", verification="inferred")
    bundle["molecules"].append(m)
    molecules[mid] = m
    return m


def touch_molecule(mid, evidence, full_name=None):
    m = molecules[mid]
    ev = list(evidence)
    occurrence = m["occurrences"][0]
    occurrence["evidence_ids"] = union(occurrence.get("evidence_ids", []), ev)
    union_fact_evidence(occurrence["compound_label"], ev)
    for value in m["names"].values():
        union_fact_evidence(value, ev)
    for key in ("cas_numbers", "smiles", "identity_note"):
        union_fact_evidence(m[key], ev)
    for value in m["stereochemistry"].values():
        union_fact_evidence(value, ev)
    if full_name:
        setf(m["names"], "original", full_name, full_name, ev)
        setf(m["names"], "english", full_name, full_name, ev)


def make_reaction(suffix, evidence, record_type="specific_experiment",
                  categories=None, location=None, reaction_name=None,
                  transformation=None, general_ids=None, unresolved=None):
    rid = f"{DOC}:{suffix}"
    ev = list(dict.fromkeys(evidence))
    r = new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    setf(r, "experiment_categories", ", ".join(categories or []), categories or [], ev)
    setf(r, "location", location, location, ev)
    setf(r, "reaction_name", reaction_name, reaction_name, ev)
    setf(r, "transformation_type", transformation, transformation, ev,
         method="inference" if transformation else "original_extraction",
         derivation="Transformation inferred from the inspected reaction drawing and product name."
         if transformation else None,
         verification="inferred" if transformation else None)
    r["general_scheme_ids"] = list(general_ids or [])
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    bundle["reactions"].append(r)
    reactions[rid] = r
    return r


def participant(r, pid, mid, role, ev, stage="reaction", *, mass=None,
                amount=None, volume=None, equiv=None, molpct=None, ref=None,
                product_role=None):
    p = new_participant(pid, mid, role, ev, stage)
    if product_role:
        setf(p, "product_role", product_role, product_role, ev)
    for key, val in (("mass", mass), ("amount", amount), ("volume", volume),
                     ("equivalents", equiv), ("mol_percent", molpct)):
        if val is not None:
            raw = val.get("raw_text") if isinstance(val, dict) else None
            setf(p["charges"], key, raw or val, val, ev)
    if ref and (equiv is not None or molpct is not None):
        setf(p, "equivalence_reference", ref, ref, ev)
    r["participants"].append(p)
    return p


def add_step(r, sid, ev, pids, operation, *, temperature=None, time=None,
             stage="reaction", atmosphere=None):
    s = new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    setf(s, "operation", operation, operation, ev)
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else temperature.get("raw_text", temperature)
        setf(s, "temperature", raw, temperature, ev)
    if time is not None:
        raw = time if isinstance(time, str) else time.get("raw_text", time)
        setf(s, "time", raw, time, ev)
    if atmosphere is not None:
        setf(s, "atmosphere", atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, ev, raw, normalized, *, basis_raw=None,
                basis=None, yield_type=None, product_role=None, scope=None):
    o = new_outcome(oid, kind, targets, ev)
    setf(o, "value", raw, normalized, ev)
    if basis_raw is not None or basis is not None:
        setf(o, "measurement_basis", basis_raw, basis, ev)
    if yield_type is not None:
        setf(o, "yield_type", yield_type, yield_type, ev)
    if product_role is not None:
        setf(o, "product_role", product_role, product_role, ev)
    if scope is not None:
        setf(o, "scope", scope, scope, ev)
    r["outcomes"].append(o)
    return o


def classify_partial(r, ev, codes, limitation):
    categories = {"C-H": 1, "C-C": 2, "C-N": 5}
    atoms = {
        "C-H": ("C", "H", 1, "hydrogen_gain"),
        "C-C": ("C", "C", 1, "formed"),
        "C-N": ("C", "N", 1, "formed"),
    }
    changes, labels = [], []
    for idx, code in enumerate(codes, 1):
        a, b, order, kind = atoms[code]
        cid = f"b{idx}"
        ch = {"change_id": cid,
              "atoms": [{"atom_map": None, "element": a},
                        {"atom_map": None, "element": b}],
              "before_order": 0, "after_order": order, "change_type": kind}
        if kind == "hydrogen_gain":
            ch["hydrogen_count"] = 1
        changes.append(F(code, ch, ev, "inference",
                         "Net formed bond inferred from the source transformation; atom mapping is unresolved.",
                         verification="inferred"))
        labels.append({"category_number": categories[code], "code": code,
                       "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = changes
    setf(r["classification"], "labels", ", ".join(codes), labels, ev,
         method="inference", derivation="Only source-supported newly formed bonds are labelled.",
         verification="inferred")
    setf(r["classification"], "coverage", "partial", "partial", ev,
         method="inference", derivation=limitation, verification="inferred")
    r["classification"]["limitations"] = [limitation, "No atom-mapped reaction SMILES was asserted."]


def note(r, text, ev):
    r["notes"].append(F(text, text, ev, "inference",
                        "Evidence-linked extraction note.", verification="inferred"))


def add_common_workup_participants(r, ev):
    present = {p["participant_id"] for p in r["participants"]}
    specs = [
        ("p_nahco3", "nahco3", "workup_agent", "workup", None, None,
         q(10, "mL", raw_text="10 mL of 5% NaHCO3")),
        ("p_dcm_extract", "dcm", "solvent", "workup", None, None,
         q(5, "mL", raw_text="2 × 5 mL", repetitions=2)),
        ("p_na2so4", "na2so4", "workup_agent", "workup", None, None, None),
        ("p_silica", "silica_gel", "purification_agent", "purification", None, None, None),
        ("p_etoac", "etoac", "purification_agent", "purification", None, None, None),
        ("p_dcm_eluent", "dcm", "purification_agent", "purification", None, None, None),
    ]
    for pid, mid, role, stage, mass, amount, volume in specs:
        if pid not in present:
            participant(r, pid, mid, role, ev, stage=stage, mass=mass,
                        amount=amount, volume=volume)


def set_charge(p, key, raw, normalized, ev):
    setf(p["charges"], key, raw, normalized, ev)


def update_product18(product, page, mass_mg, yld, eluent, heated=False,
                     ratio=None):
    rid = f"{DOC}:scheme4_scope_{product}"
    r = reactions[rid]
    ev = [E[page]]
    proc_ev = [E[12], E[page]]
    r["evidence_ids"] = union(r["evidence_ids"], proc_ev)
    gp = f"{DOC}:general_procedure_b_products_18"
    if gp not in r["general_scheme_ids"]:
        r["general_scheme_ids"].append(gp)
    by = {p["participant_id"]: p for p in r["participants"]}
    set_charge(by["p_cu"], "mass", "3.6 mg", q(3.6, "mg"), [E[12]])
    set_charge(by["p_cu"], "amount", "20 μmol", q(20, "umol"), [E[12]])
    set_charge(by["p_ligand"], "mass", "7.3 mg", q(7.3, "mg"), [E[12]])
    set_charge(by["p_ligand"], "amount", "26 μmol", q(26, "umol"), [E[12]])
    set_charge(by["p_silane"], "volume", "0.099 mL", q(0.099, "mL"), [E[12]])
    set_charge(by["p_nh4f"], "mass", "200 mg", q(200, "mg"), [E[12]])
    set_charge(by["p_meoh"], "volume", "2.5 mL", q(2.5, "mL"), [E[12]])
    by["p_tbuoh"]["charges"]["volume"] = conflict_fact(
        [F("76 μL", q(76, "uL"), [E[4]]),
         F("76.5 μL", q(76.5, "uL"), [E[12]])],
        [E[4], E[12]],
        "Scheme 4 footnote and General Procedure B print different t-BuOH volumes; neither was selected.")
    add_common_workup_participants(r, [E[12]])
    by = {p["participant_id"]: p for p in r["participants"]}
    r["steps"] = []
    add_step(r, "catalyst_mix", [E[12]],
             ["p_cu", "p_ligand", "p_solvent", "p_tbuoh"],
             "charge in an Ar-filled glovebox and stir", time=q(5, "min"),
             atmosphere="argon")
    add_step(r, "substrate_charge", [E[12]], ["p_allene", "p_imine"],
             "charge allenamide then imine, seal vial, and remove from glovebox",
             atmosphere="argon")
    add_step(r, "silane_addition", [E[12]], ["p_silane"],
             "add dimethoxymethylsilane")
    add_step(r, "reaction_hold", proc_ev,
             ["p_imine", "p_allene", "p_cu", "p_ligand", "p_tbuoh", "p_silane", "p_solvent"],
             "stir reaction mixture",
             temperature=q(65, "°C") if heated else "rt", time=q(24, "h"))
    add_step(r, "nh4f_quench", [E[12]], ["p_nh4f", "p_meoh"],
             "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"),
             stage="workup")
    add_step(r, "extract", [E[12]], ["p_nahco3", "p_dcm_extract"],
             "add 5% NaHCO3 and extract with DCM", stage="workup")
    add_step(r, "dry_concentrate", [E[12]], ["p_na2so4"],
             "dry combined organics with Na2SO4, filter, and concentrate in vacuo",
             stage="workup")
    add_step(r, "chromatography", ev,
             ["p_silica", "p_etoac", "p_dcm_eluent"],
             f"purify by silica gel chromatography ({eluent})", stage="purification")
    setf(r["procedure"], "general_procedure_reference",
         "According to General Procedure B", gp, proc_ev)
    detail = f"According to General Procedure B; {eluent} afforded {mass_mg} mg ({yld}%) of {product} as a single diastereomer"
    if ratio:
        detail += f" and as a {ratio} product mixture"
    setf(r["procedure"], "original_text", detail + ".", detail + ".", ev)
    if heated:
        setf(r["procedure"], "modifications", "stirred at 65 °C for 24 h",
             "stirred at 65 °C for 24 h", ev)
    for outcome in r["outcomes"]:
        union_fact_evidence(outcome["value"], ev)
        if outcome["outcome_id"] == "isolated_yield":
            setf(outcome, "measurement_basis",
                 f"isolated after silica gel chromatography ({eluent})",
                 f"isolated after silica gel chromatography ({eluent})", ev)
    if not any(o["outcome_id"] == "product_amount" for o in r["outcomes"]):
        add_outcome(r, "product_amount", "product_amount", [product], ev,
                    f"{mass_mg} mg", q(mass_mg, "mg"),
                    basis_raw=f"isolated after silica gel chromatography ({eluent})",
                    basis=f"isolated after silica gel chromatography ({eluent})",
                    product_role="main product")
    note(r, "Scheme 4 reports 76 μL t-BuOH, while General Procedure B reports 76.5 μL; both source values are retained as a conflict.", [E[4], E[12]])
    unresolved = "Scheme 4 and General Procedure B report t-BuOH as 76 μL and 76.5 μL, respectively."
    if unresolved not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append(unresolved)
    touch_molecule(product, ev)


product18 = {
    "18f": (13, 157, 81, "10% E.A. in DCM", False, None),
    "18g": (13, 137, 66, "10% E.A. in DCM", False, None),
    "18h": (13, 106, 53, "5% E.A. in DCM", False, None),
    "18i": (13, 88, 50, "5% E.A. in DCM", True, "91:9 18i:19i"),
    "18j": (13, 173, 82, "2.5% E.A. in DCM", True, None),
    "18k": (13, 143, 78, "40% E.A. in DCM", False, None),
    "18l": (13, 143, 73, "10% E.A. in DCM", False, None),
    "18m": (13, 127, 74, "10% E.A. in DCM", True, "95:5 18m:19m"),
    "18n": (14, 163, 81, "15% E.A. in DCM", False, None),
    "18o": (14, 151, 75, "15% E.A. in DCM", False, None),
    "18p": (14, 151, 74, "5% E.A. in DCM", False, None),
    "18q": (14, 144, 67, "5% E.A. in DCM", False, None),
    "18r": (14, 143, 76, "5% E.A. in DCM", False, None),
    "18s": (14, 154, 86, "5% E.A. in DCM", False, None),
    "18t": (14, 137, 74, "5% E.A. in DCM", False, None),
    "18u": (14, 172, 79, "3% E.A. in DCM", False, "95:5 18u:19u"),
}
for product, args in product18.items():
    update_product18(product, *args)


# General Procedure C is a source-level reusable procedure, not a scope row.
gp_c_ev = [E[14], E[15]]
add_generic_molecule("gp_c_imine", "imine used in General Procedure C", gp_c_ev,
                     "Generic imine placeholder; concrete imines are separately registered.")
add_generic_molecule("gp_c_product19", "cyclic urea product 19 from General Procedure C", gp_c_ev,
                     "Generic product placeholder; concrete urea products are separately registered.")
gp_c = make_reaction(
    "general_procedure_c_products_19", gp_c_ev, "generic_scheme", ["other"],
    "Pages 14--15, General Procedure C", "General Procedure C for the Synthesis of 19",
    "Cu-catalyzed reductive coupling followed by carbamate migration",
    [f"{DOC}:scheme2d_current_generic_reductive_coupling"],
    ["Generic procedure has no product-specific isolated outcome; concrete products are separate experiments."],
)
participant(gp_c, "p_imine", "gp_c_imine", "reactant", gp_c_ev, amount=q(400, "umol"))
participant(gp_c, "p_allene", "15a", "reactant", gp_c_ev,
            mass=q(96.6, "mg"), amount=q(480, "umol"))
participant(gp_c, "p_cu", "cu_oac2", "catalyst", gp_c_ev,
            mass=q(3.6, "mg"), amount=q(20, "umol"), molpct=q(5, "mol%"), ref="gp_c_imine")
participant(gp_c, "p_ligand", "pcy3", "ligand", gp_c_ev,
            mass=q(7.3, "mg"), amount=q(26, "umol"), molpct=q(6.5, "mol%"), ref="gp_c_imine")
participant(gp_c, "p_silane", "silane", "reagent", gp_c_ev,
            volume=q(0.099, "mL"), equiv=q(2, "equiv"), ref="gp_c_imine")
participant(gp_c, "p_toluene", "toluene", "solvent", gp_c_ev, volume=q(1.0, "mL"))
participant(gp_c, "p_product", "gp_c_product19", "product", gp_c_ev, product_role="product")
participant(gp_c, "p_nh4f", "nh4f", "workup_agent", gp_c_ev, stage="workup", mass=q(200, "mg"))
participant(gp_c, "p_meoh", "meoh", "solvent", gp_c_ev, stage="workup", volume=q(2.5, "mL"))
participant(gp_c, "p_nahco3", "nahco3", "workup_agent", gp_c_ev, stage="workup",
            volume=q(10, "mL", raw_text="10 mL of 5% NaHCO3"))
participant(gp_c, "p_dcm", "dcm", "solvent", gp_c_ev, stage="workup",
            volume=q(5, "mL", raw_text="2 × 5 mL", repetitions=2))
participant(gp_c, "p_na2so4", "na2so4", "workup_agent", gp_c_ev, stage="workup")
participant(gp_c, "p_silica", "silica_gel", "purification_agent", gp_c_ev, stage="purification")
setf(gp_c, "reaction_scale", "imine (400 μmol)",
     {"reference_molecule_id": "gp_c_imine", "amount": q(400, "umol")}, gp_c_ev)
setf(gp_c["solvent_conditions"], "composition", "toluene", ["toluene"], gp_c_ev)
setf(gp_c["solvent_conditions"], "total_volume", "1.0 mL", q(1.0, "mL"), gp_c_ev)
add_step(gp_c, "catalyst_mix", gp_c_ev, ["p_cu", "p_ligand", "p_toluene"],
         "charge in an Ar-filled glovebox and stir", time=q(5, "min"), atmosphere="argon")
add_step(gp_c, "substrate_charge", gp_c_ev, ["p_allene", "p_imine"],
         "charge allenamide then imine, seal vial, and remove from glovebox", atmosphere="argon")
add_step(gp_c, "silane_addition", gp_c_ev, ["p_silane"], "add dimethoxymethylsilane")
add_step(gp_c, "reaction_hold", gp_c_ev,
         ["p_imine", "p_allene", "p_cu", "p_ligand", "p_silane", "p_toluene"],
         "stir reaction mixture", temperature="rt", time=q(24, "h"))
add_step(gp_c, "nh4f_quench", gp_c_ev, ["p_nh4f", "p_meoh"],
         "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"), stage="workup")
add_step(gp_c, "extract", gp_c_ev, ["p_nahco3", "p_dcm"],
         "add 5% NaHCO3 and extract with DCM", stage="workup")
add_step(gp_c, "dry_concentrate", gp_c_ev, ["p_na2so4"],
         "dry combined organics with Na2SO4, filter, and concentrate in vacuo", stage="workup")
add_step(gp_c, "chromatography", gp_c_ev, ["p_silica"],
         "purify by flash chromatography on silica gel", stage="purification")
setf(gp_c["procedure"], "original_text", "General Procedure C for the Synthesis of 19",
     "General Procedure C for the Synthesis of 19", gp_c_ev)
classify_partial(gp_c, gp_c_ev, ["C-C", "C-H", "C-N"],
                 "Generic coupling/urea closure supports only net formed-bond classes; complete atom mapping is unresolved.")


product19 = {
    "19a": (15, 198, 94, "5% E.A. in DCM", False, None),
    "19b": (15, 158, 86, "10% E.A. in DCM", False, None),
    "19de": (15, 143, 71, "5% E.A. in DCM", True, None),
    "19ee": (15, 149, 74, "5% E.A. in DCM", True, None),
    "19j": (15, 131, 69, "5% E.A. in DCM", True, None),
    "19k": (15, 183, 99, "50% E.A. in DCM", False, "86:14 19k:18k"),
    "19l": (16, 184, 94, "10% E.A. in DCM", False, None),
    "19m": (16, 118, 69, "10% E.A. in DCM", True, None),
    "19n": (16, 197, 98, "10% E.A. in DCM", False, None),
    "19s": (16, 149, 83, "5% E.A. in DCM", False, None),
    "19t": (16, 184, 99, "3% E.A. in DCM", False, "92:8 19t:18t"),
    "19pe": (16, 130, 79, "5% E.A. in DCM", True, None),
}


def update_product19(product, page, mass_mg, yld, eluent, heated=False, ratio=None):
    r = reactions[f"{DOC}:table3_scope_{product}"]
    ev = [E[page]]
    proc_ev = union(gp_c_ev, ev)
    r["evidence_ids"] = union(r["evidence_ids"], proc_ev)
    gp = f"{DOC}:general_procedure_c_products_19"
    if gp not in r["general_scheme_ids"]:
        r["general_scheme_ids"].append(gp)
    by = {p["participant_id"]: p for p in r["participants"]}
    set_charge(by["p_cu"], "mass", "3.6 mg", q(3.6, "mg"), gp_c_ev)
    set_charge(by["p_cu"], "amount", "20 μmol", q(20, "umol"), gp_c_ev)
    set_charge(by["p_ligand"], "mass", "7.3 mg", q(7.3, "mg"), gp_c_ev)
    set_charge(by["p_ligand"], "amount", "26 μmol", q(26, "umol"), gp_c_ev)
    set_charge(by["p_silane"], "volume", "0.099 mL", q(0.099, "mL"), gp_c_ev)
    set_charge(by["p_nh4f"], "mass", "200 mg", q(200, "mg"), gp_c_ev)
    set_charge(by["p_meoh"], "volume", "2.5 mL", q(2.5, "mL"), gp_c_ev)
    add_common_workup_participants(r, gp_c_ev)
    r["steps"] = []
    add_step(r, "catalyst_mix", gp_c_ev, ["p_cu", "p_ligand", "p_solvent"],
             "charge in an Ar-filled glovebox and stir", time=q(5, "min"), atmosphere="argon")
    add_step(r, "substrate_charge", gp_c_ev, ["p_allene", "p_imine"],
             "charge allenamide then imine, seal vial, and remove from glovebox", atmosphere="argon")
    add_step(r, "silane_addition", gp_c_ev, ["p_silane"], "add dimethoxymethylsilane")
    add_step(r, "reaction_hold", proc_ev,
             ["p_imine", "p_allene", "p_cu", "p_ligand", "p_silane", "p_solvent"],
             "stir reaction mixture", temperature=q(65, "°C") if heated else "rt",
             time=q(24, "h"))
    add_step(r, "nh4f_quench", gp_c_ev, ["p_nh4f", "p_meoh"],
             "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"), stage="workup")
    add_step(r, "extract", gp_c_ev, ["p_nahco3", "p_dcm_extract"],
             "add 5% NaHCO3 and extract with DCM", stage="workup")
    add_step(r, "dry_concentrate", gp_c_ev, ["p_na2so4"],
             "dry combined organics with Na2SO4, filter, and concentrate in vacuo", stage="workup")
    add_step(r, "chromatography", ev, ["p_silica", "p_etoac", "p_dcm_eluent"],
             f"purify by silica gel chromatography ({eluent})", stage="purification")
    setf(r["procedure"], "general_procedure_reference",
         "According to General Procedure C", gp, proc_ev)
    detail = f"According to General Procedure C; {eluent} afforded {mass_mg} mg ({yld}%) of {product} as a single diastereomer"
    if ratio:
        detail += f" and as a {ratio} product mixture"
    setf(r["procedure"], "original_text", detail + ".", detail + ".", ev)
    if heated:
        setf(r["procedure"], "modifications", "stirred at 65 °C for 24 h",
             "stirred at 65 °C for 24 h", ev)
    for outcome in r["outcomes"]:
        union_fact_evidence(outcome["value"], ev)
        if outcome["outcome_id"] == "isolated_yield":
            setf(outcome, "measurement_basis",
                 f"isolated after silica gel chromatography ({eluent})",
                 f"isolated after silica gel chromatography ({eluent})", ev)
    if not any(o["outcome_id"] == "product_amount" for o in r["outcomes"]):
        add_outcome(r, "product_amount", "product_amount", [product], ev,
                    f"{mass_mg} mg", q(mass_mg, "mg"),
                    basis_raw=f"isolated after silica gel chromatography ({eluent})",
                    basis=f"isolated after silica gel chromatography ({eluent})",
                    product_role="main product")
    touch_molecule(product, ev)


for product, args in product19.items():
    update_product19(product, *args)


def add_specific_gp_c(product, imine, page, mass_mg, yld, full_name):
    ev = [E[page]]
    all_ev = union(gp_c_ev, ev)
    r = make_reaction(
        f"experimental_preparation_{product}", all_ev, "specific_experiment",
        ["substrate_scope"], f"Page {page}, experimental preparation of {product}",
        f"preparative synthesis of {product} according to General Procedure C",
        "Cu-catalyzed reductive coupling followed by carbamate migration",
        [f"{DOC}:scheme2d_current_generic_reductive_coupling",
         f"{DOC}:general_procedure_c_products_19"],
        ["This is a distinct preparative 400 μmol experiment, not the Table 2 NMR-yield optimization row with 106 mg of 15a and 6 mol% PCy3."],
    )
    participant(r, "p_imine", imine, "reactant", all_ev, amount=q(400, "umol"))
    participant(r, "p_allene", "15a", "reactant", all_ev,
                mass=q(96.6, "mg"), amount=q(480, "umol"))
    participant(r, "p_product", product, "product", all_ev, product_role="main product")
    participant(r, "p_cu", "cu_oac2", "catalyst", all_ev,
                mass=q(3.6, "mg"), amount=q(20, "umol"), molpct=q(5, "mol%"), ref=imine)
    participant(r, "p_ligand", "pcy3", "ligand", all_ev,
                mass=q(7.3, "mg"), amount=q(26, "umol"), molpct=q(6.5, "mol%"), ref=imine)
    participant(r, "p_silane", "silane", "reagent", all_ev,
                volume=q(0.099, "mL"), equiv=q(2, "equiv"), ref=imine)
    participant(r, "p_solvent", "toluene", "solvent", all_ev, volume=q(1.0, "mL"))
    participant(r, "p_nh4f", "nh4f", "workup_agent", all_ev, stage="workup", mass=q(200, "mg"))
    participant(r, "p_meoh", "meoh", "solvent", all_ev, stage="workup", volume=q(2.5, "mL"))
    participant(r, "p_nahco3", "nahco3", "workup_agent", all_ev, stage="workup",
                volume=q(10, "mL", raw_text="10 mL of 5% NaHCO3"))
    participant(r, "p_dcm_extract", "dcm", "solvent", all_ev, stage="workup",
                volume=q(5, "mL", raw_text="2 × 5 mL", repetitions=2))
    participant(r, "p_na2so4", "na2so4", "workup_agent", all_ev, stage="workup")
    participant(r, "p_silica", "silica_gel", "purification_agent", all_ev, stage="purification")
    participant(r, "p_etoac", "etoac", "purification_agent", all_ev, stage="purification")
    participant(r, "p_dcm_eluent", "dcm", "purification_agent", all_ev, stage="purification")
    setf(r, "reaction_scale", "imine (400 μmol)",
         {"reference_molecule_id": imine, "amount": q(400, "umol")}, all_ev)
    setf(r["solvent_conditions"], "composition", "toluene", ["toluene"], all_ev)
    setf(r["solvent_conditions"], "total_volume", "1.0 mL", q(1.0, "mL"), all_ev)
    add_step(r, "catalyst_mix", gp_c_ev, ["p_cu", "p_ligand", "p_solvent"],
             "charge in an Ar-filled glovebox and stir", time=q(5, "min"), atmosphere="argon")
    add_step(r, "substrate_charge", gp_c_ev, ["p_allene", "p_imine"],
             "charge allenamide then imine, seal vial, and remove from glovebox", atmosphere="argon")
    add_step(r, "silane_addition", gp_c_ev, ["p_silane"], "add dimethoxymethylsilane")
    add_step(r, "reaction_hold", all_ev,
             ["p_imine", "p_allene", "p_cu", "p_ligand", "p_silane", "p_solvent"],
             "stir reaction mixture", temperature=q(65, "°C"), time=q(24, "h"))
    add_step(r, "nh4f_quench", gp_c_ev, ["p_nh4f", "p_meoh"],
             "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"), stage="workup")
    add_step(r, "extract", gp_c_ev, ["p_nahco3", "p_dcm_extract"],
             "add 5% NaHCO3 and extract with DCM", stage="workup")
    add_step(r, "dry_concentrate", gp_c_ev, ["p_na2so4"],
             "dry combined organics with Na2SO4, filter, and concentrate in vacuo", stage="workup")
    add_step(r, "chromatography", ev, ["p_silica", "p_etoac", "p_dcm_eluent"],
             "purify by silica gel chromatography (5% E.A. in DCM)", stage="purification")
    setf(r["procedure"], "general_procedure_reference", "According to General Procedure C",
         f"{DOC}:general_procedure_c_products_19", all_ev)
    setf(r["procedure"], "modifications", "stirred at 65 °C for 24 h",
         "stirred at 65 °C for 24 h", ev)
    original = f"According to General Procedure C at 65 °C for 24 h; 5% E.A. in DCM afforded {mass_mg} mg ({yld}%) of {product} as a single diastereomer."
    setf(r["procedure"], "original_text", original, original, ev)
    add_outcome(r, "isolated_yield", "yield", [product], ev,
                f"{yld}%", q(yld, "%"),
                basis_raw="isolated after silica gel chromatography (5% E.A. in DCM)",
                basis="isolated after silica gel chromatography (5% E.A. in DCM)",
                yield_type="isolated yield", product_role="main product")
    add_outcome(r, "product_amount", "product_amount", [product], ev,
                f"{mass_mg} mg", q(mass_mg, "mg"),
                basis_raw="isolated after silica gel chromatography (5% E.A. in DCM)",
                basis="isolated after silica gel chromatography (5% E.A. in DCM)",
                product_role="main product")
    add_outcome(r, "diastereoselectivity", "selectivity", [product], ev,
                "single diastereomer",
                {"metric": "dr", "reported": "single diastereomer",
                 "major_configuration": "(S,S,S)", "minor_configuration": "not observed"},
                basis_raw="source reports a single diastereomer; analytical method not reported in this preparation",
                basis="single-diastereomer statement; method unreported",
                product_role="main product",
                scope={"metric": "dr", "order": "major:minor",
                       "major_configuration": "(S,S,S)", "minor_configuration": "unassigned"})
    classify_partial(r, all_ev, ["C-C", "C-H", "C-N"],
                     "Net coupling and urea closure are inferred; complete atom mapping and rearrangement cleavages are unresolved.")
    touch_molecule(product, ev, full_name)
    return r


add_specific_gp_c(
    "19pc", "9pc", 16, 125, 70,
    "(4S,5S)-1-Benzyl-3-((S)-2-hydroxy-1-phenylethyl)-5-(naphthalen-2-yl)-4-vinylimidazolidin-2-one",
)
add_specific_gp_c(
    "19pd", "9pd", 16, 132, 71,
    "(4S,5S)-1-(4-Fluorobenzyl)-3-((S)-2-hydroxy-1-phenylethyl)-5-(naphthalen-2-yl)-4-vinylimidazolidin-2-one",
)


# Pages 17--18 are contextual detail for experiment identities already present.
# Add only page-specific source context here; the records already carry these
# exact procedure evidence IDs from the previous cumulative pass.
bundle["source_context"].extend([
    {"document_id": DOC, "evidence_ids": [E[13], I[13]], "kind": "characterization_data",
     "text": "NMR/19F NMR/HRMS/Rf blocks for 18e-18m are characterization. Product-specific General Procedure B preparations for 18f-18m are linked to the existing Scheme 4 records."},
    {"document_id": DOC, "evidence_ids": [E[14], I[14]], "kind": "characterization_and_procedure",
     "text": "NMR/HRMS/Rf blocks for 18n-18u are characterization. General Procedure C begins at the bottom of the page and continues on page 15."},
    {"document_id": DOC, "evidence_ids": [E[15], I[15]], "kind": "characterization_and_procedure",
     "text": "General Procedure C and product-specific preparations 19a, 19b, 19de, 19ee, 19j, and 19k are reaction evidence; standalone spectral blocks are characterization."},
    {"document_id": DOC, "evidence_ids": [E[16], I[16]], "kind": "characterization_and_procedure",
     "text": "Preparations 19l, 19m, 19n, 19s, 19t, 19pe, 19pc, and 19pd are reaction evidence; their spectral blocks are characterization, not mass-spectrum reactions."},
    {"document_id": DOC, "evidence_ids": [E[17], I[17]], "kind": "detailed_experimental_context",
     "text": "Detailed procedures for 29, 30, gram-scale 18c, 25, 19a from 18a, and the start of 26 belong to experiment identities already established on earlier pages and are not duplicated."},
    {"document_id": DOC, "evidence_ids": [E[18], I[18]], "kind": "article_back_matter",
     "text": "The continuation of 26 and preparation of 27 are linked to existing experiments. The SI availability notice, CCDC accession, author information, acknowledgments, and references are nonreaction context; no SI attachment is established."},
])


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 13, "state": "attempted",
         "regions": [
             {"region": "top continuation of 18e characterization", "finding": "Tail of 18e NMR/HRMS characterization; no new synthesis experiment.", "evidence_ids": [E[13], I[13]]},
             {"region": "two-column product blocks 18f--18m", "finding": "Eight product-specific General Procedure B preparations, including 65 °C variants and the 18i/19i and 18m/19m mixture ratios; spectral blocks remain characterization.", "evidence_ids": [E[13], I[13]]},
         ],
         "record_ids": [f"{DOC}:scheme4_scope_18{x}" for x in "fghijklm"],
         "unresolved": ["Complete source-verified isomeric product graphs are not encoded; ratio measurement methods for 18i/19i and 18m/19m are not reported.", "Scheme 4 and General Procedure B print 76 and 76.5 μL t-BuOH, respectively; both values remain a conflict."]},
        {"page": 14, "state": "attempted",
         "regions": [
             {"region": "product blocks 18n--18u", "finding": "Eight product-specific General Procedure B preparations, including the 95:5 18u/rearranged-product mixture; NMR/HRMS/Rf blocks are characterization.", "evidence_ids": [E[14], I[14]]},
             {"region": "bottom General Procedure C opening", "finding": "Beginning of the reusable urea-product procedure: Cu(OAc)2/PCy3/toluene catalyst mixture under argon, continued on page 15.", "evidence_ids": [E[14], I[14]]},
         ],
         "record_ids": [f"{DOC}:scheme4_scope_18{x}" for x in "nopqrstu"] + [f"{DOC}:general_procedure_c_products_19"],
         "unresolved": ["Complete product graphs remain unencoded; the 18u mixture-ratio measurement method is unreported.", "The 76/76.5 μL t-BuOH conflict is retained without correction."]},
        {"page": 15, "state": "attempted",
         "regions": [
             {"region": "top General Procedure C continuation", "finding": "Ordered substrate/silane addition, room-temperature hold, NH4F/MeOH quench, extraction, drying, and chromatography procedure.", "evidence_ids": [E[15], I[15]]},
             {"region": "product blocks 19a, 19b, 19de, 19ee, 19j, 19k", "finding": "Six product-specific preparations, including three 65 °C variants and the 86:14 19k/18k mixture; spectral data are characterization.", "evidence_ids": [E[15], I[15]]},
         ],
         "record_ids": [f"{DOC}:general_procedure_c_products_19"] + [f"{DOC}:table3_scope_{x}" for x in ["19a", "19b", "19de", "19ee", "19j", "19k"]],
         "unresolved": ["Complete source-verified stereochemical graphs remain unencoded; the 19k/18k product-ratio method is unreported."]},
        {"page": 16, "state": "attempted",
         "regions": [
             {"region": "product blocks 19l, 19m, 19n, 19s, 19t, 19pe", "finding": "Six Table 3 product-specific preparations, including heated variants and the 92:8 19t/18t mixture; spectral data are characterization.", "evidence_ids": [E[16], I[16]]},
             {"region": "lower product blocks 19pc and 19pd", "finding": "Two preparative General Procedure C experiments at 65 °C for 24 h, distinct from the page-5 Table 2 NMR-yield optimization rows.", "evidence_ids": [E[16], I[16]]},
         ],
         "record_ids": [f"{DOC}:table3_scope_{x}" for x in ["19l", "19m", "19n", "19s", "19t", "19pe"]] + [f"{DOC}:experimental_preparation_19pc", f"{DOC}:experimental_preparation_19pd"],
         "unresolved": ["The detailed 19pc/19pd entries say single diastereomer but do not report the analytical method in those preparation paragraphs.", "Complete stereochemical SMILES remain unencoded; the 19t/18t product-ratio method is unreported."]},
        {"page": 17, "state": "attempted",
         "regions": [
             {"region": "top continuation of 19pd characterization", "finding": "19pd spectral continuation; not a separate synthesis.", "evidence_ids": [E[17], I[17]]},
             {"region": "detailed syntheses of 29 and 30", "finding": "Full charges, workups, isolated yields, and chiral-HPLC er bases for two already identified Scheme 7 mechanistic-study experiments.", "evidence_ids": [E[17], I[17]]},
             {"region": "gram-scale 18c, 25, 19a from 18a, and start of 26", "finding": "Detailed procedures for existing Scheme 5/6 experiment identities; the 19a-to-25 value is an overall three-stage yield.", "evidence_ids": [E[17], I[17]]},
         ],
         "record_ids": [f"{DOC}:experimental_preparation_19pd", f"{DOC}:scheme7_synthesis_29", f"{DOC}:scheme7_synthesis_30", f"{DOC}:scheme6_gram_scale_18c", f"{DOC}:scheme6_19a_to_25", f"{DOC}:scheme5_18a_to_19a", f"{DOC}:scheme6_18c_to_26"],
         "unresolved": ["Scheme 7 and detailed charges give conflicting catalyst/ligand loadings; 30 has 62% in Scheme 7 versus 60% in the detailed experiment. Existing conflict facts are retained.", "No semantic coverage beyond the stated page is inferred from cross-page links."]},
        {"page": 18, "state": "attempted",
         "regions": [
             {"region": "top continuation of 26 and synthesis of 27", "finding": "Completion of N-allylation product 26 and a distinct Hoveyda-Grubbs-II ring-closing metathesis to 27; both map to existing Scheme 6 records.", "evidence_ids": [E[18], I[18]]},
             {"region": "Associated Content, authors, acknowledgments, and references", "finding": "Back matter and cited literature, not new experiments. The paper advertises online Supporting Information, but no SI attachment is established in source mapping.", "evidence_ids": [E[18], I[18]]},
             {"region": "NMR/HRMS blocks for 26 and 27", "finding": "Characterization only; not represented as synthetic or mass-spectrum experiments.", "evidence_ids": [E[18], I[18]]},
         ],
         "record_ids": [f"{DOC}:scheme6_18c_to_26", f"{DOC}:scheme6_26_to_27"],
         "unresolved": ["No Supporting Information file is explicitly attached; none was associated.", "Complex products 26 and 27 remain without source-verified complete isomeric SMILES."]},
    ],
}


for r in bundle["reactions"]:
    r["review"]["status"] = "needs_review"
    if r["review"].get("method") in (None, "none"):
        r["review"]["method"] = "model"

write(ROOT / "runs/current/submitted.json", bundle)
write(ROOT / "runs/current/coverage/pages_0013_0018.json", coverage)
write(ROOT / "runs/current/raw_model/pages_0013_0018_inventory.json", {
    "document_id": DOC,
    "pages": [13, 14, 15, 16, 17, 18],
    "new_reaction_ids": [
        f"{DOC}:general_procedure_c_products_19",
        f"{DOC}:experimental_preparation_19pc",
        f"{DOC}:experimental_preparation_19pd",
    ],
    "updated_reaction_ids": [f"{DOC}:scheme4_scope_{x}" for x in product18]
    + [f"{DOC}:table3_scope_{x}" for x in product19],
    "existing_contextual_record_ids": [
        f"{DOC}:scheme7_synthesis_29", f"{DOC}:scheme7_synthesis_30",
        f"{DOC}:scheme6_gram_scale_18c", f"{DOC}:scheme6_19a_to_25",
        f"{DOC}:scheme5_18a_to_19a", f"{DOC}:scheme6_18c_to_26",
        f"{DOC}:scheme6_26_to_27",
    ],
    "new_reaction_count": 3,
    "updated_reaction_count": len(product18) + len(product19),
    "notes": [
        "All six original page images were viewed individually at original detail.",
        "The preparative 19pc/19pd experiments were kept separate from Table 2 NMR optimization rows.",
        "Pages 17--18 procedures were linked to existing experiment identities and not duplicated.",
        "No SI was associated or inferred.",
        "Standalone NMR/HRMS/Rf blocks remain source_context.",
    ],
})

print(json.dumps({
    "cumulative_reactions": len(bundle["reactions"]),
    "cumulative_molecules": len(bundle["molecules"]),
    "new_reactions": 3,
    "updated_reactions": len(product18) + len(product19),
}, ensure_ascii=False))
