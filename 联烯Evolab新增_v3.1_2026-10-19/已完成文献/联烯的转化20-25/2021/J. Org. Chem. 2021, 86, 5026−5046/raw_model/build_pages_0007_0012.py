#!/usr/bin/env python3
"""Extend the current v3.1 extraction with original PDF pages 7--12.

The script preserves the cumulative current-document registry, adds the
page-7--12 records, enriches same-identity records from pages 1--6, and keeps
all machine records at needs_review.  Page 17 text is used only as contextual
procedural evidence for the two Scheme 7 experiments; no page-17 coverage is
claimed.
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
    stereo_inspect,
    write,
)

DOC = "doc_c32350973ac9b683"
E = {p: f"{DOC}:page{p}:text" for p in range(1, 22)}
I = {p: f"{DOC}:page{p}:image" for p in range(1, 22)}
dataset_path = ROOT / "runs/current/dataset.json"
dataset = json.loads(dataset_path.read_text())
bundle = copy.deepcopy(dataset["standard"])
# The cumulative standard view predates page-crop registration.  Merge requires
# the immutable document/evidence payload produced by `standard init`, which is
# the current top-level dataset document payload.
bundle["documents"] = copy.deepcopy(dataset["documents"])
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}
stereo_diagnostics = {}


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


def union_evidence(obj, evidence):
    obj["evidence_ids"] = list(dict.fromkeys(obj.get("evidence_ids", []) + list(evidence)))


def conflict_fact(alternatives, evidence, derivation):
    return F(None, None, evidence, "original_extraction", derivation,
             status="conflict", verification="unreviewed", alternatives=alternatives)


def add_or_update_molecule(mid, label, evidence, name, *, smiles=None,
                           chirality="achiral", geometry="not_applicable",
                           assignment=None, mixture=None, note=None,
                           representation=None):
    ev = list(dict.fromkeys(evidence))
    if mid in molecules:
        m = molecules[mid]
        union_evidence(m["occurrences"][0], ev)
        m["occurrences"][0]["compound_label"]["evidence_ids"] = list(
            dict.fromkeys(m["occurrences"][0]["compound_label"]["evidence_ids"] + ev)
        )
    else:
        m = new_molecule(mid, DOC, ev)
        bundle["molecules"].append(m)
        molecules[mid] = m
    m["identity_status"] = "source_identified"
    setf(m["occurrences"][0], "compound_label", label, label, ev)
    setf(m["names"], "original", name, name, ev)
    setf(m["names"], "english", name, name, ev)
    if m["names"]["abbreviations"]["normalized_value"] is None:
        setf(m["names"], "abbreviations", None, None, ev,
             status="unreported", verification="missing")
    if m["names"]["aliases"]["normalized_value"] is None:
        setf(m["names"], "aliases", None, None, ev,
             status="unreported", verification="missing")
    setf(m, "cas_numbers", None, None, ev, status="unreported", verification="missing")
    setf(m["stereochemistry"], "chirality_status", chirality, chirality, ev,
         method="inference" if chirality == "achiral" else "original_extraction",
         derivation=("No stereogenic center or axis is present in the source-named candidate graph."
                     if chirality == "achiral" else None),
         verification="inferred" if chirality == "achiral" else None)
    setf(m["stereochemistry"], "geometry_status", geometry, geometry, ev)
    if smiles:
        diag = stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        setf(m, "smiles", name, diag["canonical_isomeric_smiles"], ev,
             method="structure_derivation",
             derivation="Machine candidate derived from the source name and inspected structure; RDKit parsing is not source verification.",
             verification="inferred")
        setf(m["stereochemistry"], "isomeric_smiles", name,
             diag["canonical_isomeric_smiles"], ev, method="structure_derivation",
             derivation="RDKit canonicalization of the evidence-linked machine structure candidate.",
             verification="inferred")
        setf(m["stereochemistry"], "identity_key", name,
             diag["stereo_graph_key"], ev, method="structure_derivation",
             derivation="Stereo graph key of the machine candidate only.", verification="inferred")
        setf(m["stereochemistry"], "representation_status", "encoded", "encoded", ev,
             method="inference",
             derivation="All source-reported E imine geometry is encoded in the candidate when applicable; structure remains unverified.",
             verification="inferred")
    else:
        m["representation_type"] = representation or "unresolved"
        setf(m, "smiles", None, None, ev, status="unreported", verification="missing")
        setf(m["stereochemistry"], "isomeric_smiles", None, None, ev,
             status="unreported", verification="missing")
        setf(m["stereochemistry"], "identity_key", None, None, ev,
             status="unreported", verification="missing")
        setf(m["stereochemistry"], "representation_status", "not_encoded", "not_encoded", ev,
             method="inference", derivation="A complete reliable stereochemical graph was not encoded.",
             verification="inferred")
    setf(m["stereochemistry"], "atom_centers", None, None, ev,
         status="not_applicable" if chirality == "achiral" else "unreported",
         verification="missing")
    setf(m["stereochemistry"], "double_bonds",
         assignment if geometry == "specified" else None,
         assignment if geometry == "specified" else None, ev,
         status=None if geometry == "specified" else "not_applicable",
         verification=None if geometry == "specified" else "missing")
    setf(m["stereochemistry"], "axial", None, None, ev,
         status="not_applicable" if chirality == "achiral" else "unreported",
         verification="missing")
    setf(m["stereochemistry"], "assignment", assignment, assignment, ev,
         status=None if assignment else "unreported",
         verification=None if assignment else "missing")
    setf(m["stereochemistry"], "mixture", mixture, mixture, ev,
         status=None if mixture else "not_applicable", verification=None if mixture else "missing")
    setf(m, "identity_note", note or "Identity follows the source name/drawing; no CAS number was supplied.",
         note or "Identity follows the source name/drawing; no CAS number was supplied.", ev,
         method="inference", derivation="Evidence-linked machine identity note.", verification="inferred")
    return m


def add_generic_molecule(mid, label, evidence, reason):
    return add_or_update_molecule(mid, label, evidence, label, chirality="unresolved",
                                  geometry="unresolved", representation="unresolved",
                                  note=reason)


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
         derivation="Transformation description inferred from the inspected source reaction."
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
            setf(p["charges"], key, val.get("raw_text", val), val, ev)
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
        setf(s, "temperature", temperature if isinstance(temperature, str) else temperature.get("raw_text", temperature), temperature, ev)
    if time is not None:
        setf(s, "time", time if isinstance(time, str) else time.get("raw_text", time), time, ev)
    if atmosphere is not None:
        setf(s, "atmosphere", atmosphere, atmosphere, ev)
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, ev, raw, normalized, *, basis_raw=None,
                basis=None, yield_type=None, product_role=None, scope=None,
                value_fact=None):
    o = new_outcome(oid, kind, targets, ev)
    if value_fact is not None:
        o["value"] = value_fact
    else:
        setf(o, "value", raw, normalized, ev)
    if basis is not None or basis_raw is not None:
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
    categories = {"C-H": 1, "C-C": 2, "C-N": 5, "DB-C-N": 16}
    atoms = {
        "C-H": ("C", "H", 1, "hydrogen_gain"),
        "C-C": ("C", "C", 1, "formed"),
        "C-N": ("C", "N", 1, "formed"),
        "DB-C-N": ("C", "N", 2, "formed"),
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
                         "Net change inferred from the source drawing/name; atom mapping is unresolved.",
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


def note(r, text, ev, method="inference"):
    r["notes"].append(F(text, text, ev, method,
                        "Evidence-linked extraction note." if method == "inference" else None,
                        verification="inferred" if method == "inference" else None))


# Exact reagents, workup agents, aldehydes, and amines introduced on pages 10--12.
small = [
    ("mgso4", "anhydrous magnesium sulfate", "[Mg+2].[O-]S(=O)(=O)[O-]"),
    ("nahco3", "sodium bicarbonate", "[Na+].O=C([O-])O"),
    ("na2so4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    ("etoac", "ethyl acetate", "CCOC(C)=O"),
    ("hexanes", "hexanes", "CCCCCC"),
    ("dmb_amine", "2,4-dimethoxybenzylamine", "COc1ccc(CN)c(OC)c1"),
    ("pcf3_benzylamine", "4-(trifluoromethyl)benzylamine", "NCc1ccc(C(F)(F)F)cc1"),
    ("pf_benzylamine", "4-fluorobenzylamine", "NCc1ccc(F)cc1"),
    ("ald_4cl", "4-chlorobenzaldehyde", "O=Cc1ccc(Cl)cc1"),
    ("ald_4f", "4-fluorobenzaldehyde", "O=Cc1ccc(F)cc1"),
    ("ald_4cn", "4-formylbenzonitrile", "N#Cc1ccc(C=O)cc1"),
    ("ald_4co2me", "methyl 4-formylbenzoate", "COC(=O)c1ccc(C=O)cc1"),
    ("ald_3pyr", "3-pyridinecarboxaldehyde", "O=Cc1cccnc1"),
    ("ald_piperonal", "piperonal", "O=Cc1ccc2c(c1)OCO2"),
    ("ald_dhbenzofuran", "2,3-dihydrobenzofuran-5-carbaldehyde", "O=Cc1ccc2c(c1)CCO2"),
    ("ald_2naph", "2-naphthaldehyde", "O=Cc1ccc2ccccc2c1"),
    ("ald_4biphenyl", "4-phenylbenzaldehyde", "O=Cc1ccc(-c2ccccc2)cc1"),
    ("ald_ptol", "p-tolualdehyde", "Cc1ccc(C=O)cc1"),
    ("ald_furfural", "furfural", "O=Cc1ccco1"),
    ("ald_thiophene", "thiophene-2-carboxaldehyde", "O=Cc1cccs1"),
]
for mid, name, smi in small:
    add_or_update_molecule(mid, name, [E[10], E[11], E[12]], name,
                           smiles=smi, chirality="achiral", geometry="not_applicable")

add_or_update_molecule("celite", "Celite", [E[10], E[12]], "Celite filtration aid",
                       chirality="not_applicable", geometry="not_applicable",
                       representation="unresolved",
                       note="Heterogeneous filtration aid; no unique molecular SMILES is applicable.")
add_or_update_molecule("silica_gel", "silica gel", [E[12], E[17]], "silica gel",
                       chirality="not_applicable", geometry="not_applicable",
                       representation="unresolved",
                       note="Heterogeneous purification material; no unique molecular SMILES is applicable.")
add_or_update_molecule("ph_bpe_ss", "(S,S)-Ph-BPE", [E[9], E[10], E[17]],
                       "(S,S)-Ph-BPE chiral bisphosphine ligand",
                       chirality="specified", geometry="not_applicable",
                       assignment="(S,S) ligand configuration reported in Scheme 7",
                       note="Chiral ligand is source identified; complete isomeric SMILES was not encoded.")

# Update the fifteen explicitly prepared E-imines with names and machine structure candidates.
imine_info = {
    "9d": ("(E)-1-(4-chlorophenyl)-N-(2,4-dimethoxybenzyl)methanimine", "COc1ccc(C/N=C/c2ccc(Cl)cc2)c(OC)c1"),
    "9de": ("(E)-1-(4-chlorophenyl)-N-(4-(trifluoromethyl)benzyl)methanimine", "FC(F)(F)c1ccc(C/N=C/c2ccc(Cl)cc2)cc1"),
    "9e": ("(E)-N-(2,4-dimethoxybenzyl)-1-(4-fluorophenyl)methanimine", "COc1ccc(C/N=C/c2ccc(F)cc2)c(OC)c1"),
    "9f": ("(E)-4-(2,4-dimethoxybenzyliminomethyl)benzonitrile", "COc1ccc(C/N=C/c2ccc(C#N)cc2)c(OC)c1"),
    "9g": ("methyl (E)-4-(2,4-dimethoxybenzyliminomethyl)benzoate", "COC(=O)c1ccc(/C=N/Cc2ccc(OC)cc2OC)cc1"),
    "9k": ("(E)-N-(2,4-dimethoxybenzyl)-1-(pyridin-3-yl)methanimine", "COc1ccc(C/N=C/c2cccnc2)c(OC)c1"),
    "9n": ("(E)-1-(benzo[d][1,3]dioxol-5-yl)-N-(2,4-dimethoxybenzyl)methanimine", "COc1ccc(C/N=C/c2ccc3c(c2)OCO3)c(OC)c1"),
    "9o": ("(E)-1-(2,3-dihydrobenzofuran-5-yl)-N-(2,4-dimethoxybenzyl)methanimine", "COc1ccc(C/N=C/c2ccc3c(c2)CCO3)c(OC)c1"),
    "9pa": ("(E)-N-(2,4-dimethoxybenzyl)-1-(naphthalen-2-yl)methanimine", "COc1ccc(C/N=C/c2ccc3ccccc3c2)c(OC)c1"),
    "9pd": ("(E)-N-(4-fluorobenzyl)-1-(naphthalen-2-yl)methanimine", "Fc1ccc(C/N=C/c2ccc3ccccc3c2)cc1"),
    "9pe": ("(E)-1-(naphthalen-2-yl)-N-(4-(trifluoromethyl)benzyl)methanimine", "FC(F)(F)c1ccc(C/N=C/c2ccc3ccccc3c2)cc1"),
    "9q": ("(E)-1-([1,1'-biphenyl]-4-yl)-N-(2,4-dimethoxybenzyl)methanimine", "COc1ccc(C/N=C/c2ccc(-c3ccccc3)cc2)c(OC)c1"),
    "9r": ("(E)-N-(2,4-dimethoxybenzyl)-1-(p-tolyl)methanimine", "COc1ccc(C/N=C/c2ccc(C)cc2)c(OC)c1"),
    "9s": ("(E)-N-(2,4-dimethoxybenzyl)-1-(furan-2-yl)methanimine", "COc1ccc(C/N=C/c2ccco2)c(OC)c1"),
    "9t": ("(E)-N-(2,4-dimethoxybenzyl)-1-(thiophen-2-yl)methanimine", "COc1ccc(C/N=C/c2cccs2)c(OC)c1"),
}
for mid, (name, smi) in imine_info.items():
    pages = [10] if mid in {"9d", "9de"} else ([11, 12] if mid == "9q" else [11] if mid in {"9e", "9f", "9g", "9k", "9n", "9o", "9pa", "9pd", "9pe"} else [12])
    add_or_update_molecule(mid, mid, [E[p] for p in pages], name, smiles=smi,
                           chirality="achiral", geometry="specified",
                           assignment="E imine geometry reported in the source compound name")

# Products 29 and 30: enantiomer mixtures, single relative diastereomers, not fully encoded.
add_or_update_molecule("29", "29", [E[9], E[10], E[17]],
                       "branched diamine product 29 from achiral allenamide 15b",
                       chirality="mixture", geometry="not_applicable",
                       assignment="single relative diastereomer; absolute configurations of the enantiomer pair are not assigned",
                       mixture={"metric": "er", "ratio": "57:43", "order": "first-eluting:second-eluting enantiomer not assigned"},
                       note="Scheme 7/source reports one diastereomer as a 57:43 enantiomer mixture; full isomeric SMILES is not encoded.")
add_or_update_molecule("30", "30", [E[9], E[10], E[17]],
                       "cyclic urea product 30 from achiral allenamide 15b",
                       chirality="mixture", geometry="not_applicable",
                       assignment="single relative diastereomer; absolute configurations of the enantiomer pair are not assigned",
                       mixture={"metric": "er", "ratio": "80:20", "order": "first-eluting:second-eluting enantiomer not assigned"},
                       note="Scheme 7/source reports one diastereomer as an 80:20 enantiomer mixture; full isomeric SMILES is not encoded.")


# General Procedure A as a reusable generic procedure record.
gp_a_ev = [E[10]]
add_generic_molecule("gp_a_aldehyde", "aldehyde used in General Procedure A", gp_a_ev,
                     "Generic aldehyde placeholder; concrete structures are registered in each specific experiment.")
add_generic_molecule("gp_a_amine", "benzylamine used in General Procedure A", gp_a_ev,
                     "Generic amine placeholder; concrete DMB/4-F/4-CF3 benzylamines are separately registered.")
add_generic_molecule("gp_a_imine", "E-imine product of General Procedure A", gp_a_ev,
                     "Generic product placeholder; concrete products 9d--9t are separately registered.")
gp_a = make_reaction("general_procedure_a_imines", gp_a_ev, "generic_scheme", ["other"],
                     "Page 10, General Procedure A", "General Procedure A for the Synthesis of Imines",
                     "aldehyde–amine condensation to an imine",
                     unresolved=["Generic procedure uses variable aldehydes/amines; concrete structures and yields belong to product-specific records."])
participant(gp_a, "p_aldehyde", "gp_a_aldehyde", "reactant", gp_a_ev,
            amount=q(6.0, "mmol"), equiv=q(1.0, "equiv"), ref="gp_a_aldehyde")
participant(gp_a, "p_amine", "gp_a_amine", "reactant", gp_a_ev,
            amount=q(6.0, "mmol"), equiv=q(1.0, "equiv"), ref="gp_a_aldehyde")
participant(gp_a, "p_mgso4", "mgso4", "reagent", gp_a_ev)
participant(gp_a, "p_dcm", "dcm", "solvent", gp_a_ev, volume=q(8, "mL"))
participant(gp_a, "p_product", "gp_a_imine", "product", gp_a_ev, product_role="product")
participant(gp_a, "p_celite", "celite", "workup_agent", gp_a_ev, stage="workup")
setf(gp_a, "reaction_scale", "aldehyde (6.0 mmol, 1.0 equiv)",
     {"reference_molecule_id": "gp_a_aldehyde", "amount": q(6.0, "mmol")}, gp_a_ev)
setf(gp_a["solvent_conditions"], "composition", "dichloromethane", ["dcm"], gp_a_ev)
setf(gp_a["solvent_conditions"], "total_volume", "8 mL", q(8, "mL"), gp_a_ev)
add_step(gp_a, "charge_and_add", gp_a_ev, ["p_aldehyde", "p_dcm", "p_mgso4", "p_amine"],
         "charge aldehyde and dichloromethane, add anhydrous magnesium sulfate, then add amine dropwise while stirring",
         atmosphere="nitrogen")
add_step(gp_a, "hold", gp_a_ev, ["p_aldehyde", "p_amine", "p_mgso4", "p_dcm"],
         "stir reaction mixture", temperature="room temperature", time=q(12, "h"), atmosphere="nitrogen")
add_step(gp_a, "filter", gp_a_ev, ["p_celite", "p_mgso4"],
         "filter crude reaction mixture through Celite to remove magnesium sulfate", stage="workup")
add_step(gp_a, "concentrate", gp_a_ev, [], "concentrate filtrate in vacuo", stage="workup")
setf(gp_a["procedure"], "original_text", "General Procedure A for the Synthesis of Imines", "General Procedure A for the Synthesis of Imines", gp_a_ev)
classify_partial(gp_a, gp_a_ev, ["DB-C-N"], "Generic condensation supports C=N formation but not atom-level mapping or explicit water accounting.")


imine_runs = {
    "9d": (10, "ald_4cl", 0.84, "dmb_amine", 1.0, 1.54, 89, None),
    "9de": (10, "ald_4cl", 0.84, "pcf3_benzylamine", 1.05, 1.5, 84, None),
    "9e": (11, "ald_4f", 0.742, "dmb_amine", 1.0, 1.36, 84, None),
    "9f": (11, "ald_4cn", 0.784, "dmb_amine", 1.0, 1.59, 95, None),
    "9g": (11, "ald_4co2me", 0.982, "dmb_amine", 1.0, 1.87, 100, None),
    "9k": (11, "ald_3pyr", 0.64, "dmb_amine", 1.0, 1.51, 99, None),
    "9n": (11, "ald_piperonal", 0.89, "dmb_amine", 1.0, 1.77, 99, None),
    "9o": (11, "ald_dhbenzofuran", 0.886, "dmb_amine", 1.0, 1.94, 78, 72),
    "9pa": (11, "ald_2naph", 0.937, "dmb_amine", 1.0, 1.78, 98, None),
    "9pd": (11, "ald_2naph", 0.937, "pf_benzylamine", 1.05, 1.28, 73, 90),
    "9pe": (11, "ald_2naph", 0.937, "pcf3_benzylamine", 1.05, 1.71, 91, None),
    "9q": (12, "ald_4biphenyl", 1.09, "dmb_amine", 1.0, 1.82, 86, 94),
    "9r": (12, "ald_ptol", 0.72, "dmb_amine", 1.0, 1.28, 80, None),
    "9s": (12, "ald_furfural", 0.57, "dmb_amine", 1.0, 1.6, 99, 91),
    "9t": (12, "ald_thiophene", 0.67, "dmb_amine", 1.0, 1.15, 74, None),
}
for product, (page, aldehyde, aldehyde_mass, amine, amine_mass, product_mass, yld, purity) in imine_runs.items():
    ev = [E[page]] if product != "9q" else [E[11], E[12]]
    r = make_reaction(f"imine_synthesis_{product}", ev, "specific_experiment", ["substrate_scope"],
                      f"Page {page}, experimental preparation of {product}",
                      f"synthesis of imine {product} following General Procedure A",
                      "aldehyde–amine condensation to an E-imine",
                      [f"{DOC}:general_procedure_a_imines"],
                      ["No SI is associated; source-reported characterization is retained as source context rather than reaction steps."])
    participant(r, "p_aldehyde", aldehyde, "reactant", ev,
                mass=q(aldehyde_mass, "g"), amount=q(6.0, "mmol"), equiv=q(1.0, "equiv"), ref=aldehyde)
    participant(r, "p_amine", amine, "reactant", ev,
                mass=q(amine_mass, "g"), amount=q(6.0, "mmol"), equiv=q(1.0, "equiv"), ref=aldehyde)
    participant(r, "p_mgso4", "mgso4", "reagent", ev, mass=q(2.0, "g"))
    participant(r, "p_dcm", "dcm", "solvent", ev, volume=q(8, "mL"))
    participant(r, "p_product", product, "product", ev, product_role="main product")
    participant(r, "p_celite", "celite", "workup_agent", ev, stage="workup")
    setf(r, "reaction_scale", "aldehyde (6.0 mmol, 1.0 equiv)",
         {"reference_molecule_id": aldehyde, "amount": q(6.0, "mmol")}, ev)
    setf(r["solvent_conditions"], "composition", "dichloromethane", ["dcm"], ev,
         method="inherited_general_condition")
    setf(r["solvent_conditions"], "total_volume", "8 mL", q(8, "mL"), ev)
    add_step(r, "charge_and_add", ev, ["p_aldehyde", "p_dcm", "p_mgso4", "p_amine"],
             "charge aldehyde and dichloromethane, add anhydrous magnesium sulfate, then add amine dropwise while stirring",
             atmosphere="nitrogen")
    add_step(r, "hold", ev, ["p_aldehyde", "p_amine", "p_mgso4", "p_dcm"],
             "stir reaction mixture", temperature="room temperature", time=q(12, "h"), atmosphere="nitrogen")
    add_step(r, "filter", ev, ["p_celite", "p_mgso4"],
             "filter through Celite to remove magnesium sulfate", stage="workup")
    add_step(r, "concentrate", ev, [], "concentrate filtrate in vacuo", stage="workup")
    setf(r["procedure"], "general_procedure_reference", "Following General Procedure A", f"{DOC}:general_procedure_a_imines", ev)
    setf(r["procedure"], "original_text", f"Following General Procedure A; {product} obtained as reported material.",
         f"Following General Procedure A; {product} obtained as reported material.", ev)
    basis = "material obtained after Celite filtration and concentration in vacuo"
    if purity is not None:
        basis += f"; source reports {purity}% purity and does not state whether yield is purity-corrected"
        r["review"]["unresolved"].append(f"{product}: {yld}% yield is reported with {purity}% purity; purity correction basis is unreported.")
        note(r, f"General Procedure A says the filtrate was concentrated to yield the pure imine, while this entry explicitly reports {purity}% purity; both statements are retained.", ev)
    add_outcome(r, "isolated_yield", "yield", [product], ev, f"{yld}%", q(yld, "%"),
                basis_raw=basis, basis=basis, yield_type="reported isolated yield", product_role="main product")
    add_outcome(r, "product_amount", "product_amount", [product], ev, f"{product_mass} g", q(product_mass, "g"),
                basis_raw=basis, basis=basis, product_role="main product")
    if purity is not None:
        add_outcome(r, "reported_purity", "other", [product], ev, f"{purity}% purity", q(purity, "%"),
                    basis_raw="source-reported product purity; analytical method not reported",
                    basis="source-reported purity; method unreported", product_role="main product",
                    scope={"metric": "purity"})
    classify_partial(r, ev, ["DB-C-N"], "The product name/structure supports imine C=N formation; dehydration atoms are not mapped.")


# General Procedure B: a reusable procedure record, then detailed same-ID updates for 18a--18e.
gp_b_ev = [E[12]]
add_generic_molecule("gp_b_imine", "imine used in General Procedure B", gp_b_ev,
                     "Generic imine placeholder; each concrete scope substrate is separately registered.")
add_generic_molecule("gp_b_product18", "branched product 18 from General Procedure B", gp_b_ev,
                     "Generic product placeholder; each concrete product is separately registered.")
gp_b = make_reaction("general_procedure_b_products_18", gp_b_ev, "generic_scheme", ["other"],
                     "Page 12, General Procedure B", "General Procedure B for the Synthesis of 18",
                     "Cu-catalyzed reductive coupling of an imine and chiral allenamide 15a",
                     [f"{DOC}:scheme2d_current_generic_reductive_coupling"],
                     ["Generic procedure has no product-specific isolated yield; those outcomes belong to the concrete scope experiments."])
participant(gp_b, "p_imine", "gp_b_imine", "reactant", gp_b_ev, amount=q(400, "umol"))
participant(gp_b, "p_allene", "15a", "reactant", gp_b_ev, mass=q(96.6, "mg"), amount=q(480, "umol"))
participant(gp_b, "p_cu", "cu_oac2", "catalyst", gp_b_ev, mass=q(3.6, "mg"), amount=q(20, "umol"), molpct=q(5, "mol%"), ref="gp_b_imine")
participant(gp_b, "p_pcy3", "pcy3", "ligand", gp_b_ev, mass=q(7.3, "mg"), amount=q(26, "umol"), molpct=q(6.5, "mol%"), ref="gp_b_imine")
participant(gp_b, "p_tbuoh", "tbuoh", "additive", gp_b_ev, volume=q(76.5, "uL"), equiv=q(2, "equiv"), ref="gp_b_imine")
participant(gp_b, "p_silane", "silane", "reagent", gp_b_ev, volume=q(0.099, "mL"), equiv=q(2, "equiv"), ref="gp_b_imine")
participant(gp_b, "p_toluene", "toluene", "solvent", gp_b_ev, volume=q(1.0, "mL"))
participant(gp_b, "p_product", "gp_b_product18", "product", gp_b_ev, product_role="product")
participant(gp_b, "p_nh4f", "nh4f", "workup_agent", gp_b_ev, stage="workup", mass=q(200, "mg"))
participant(gp_b, "p_meoh", "meoh", "solvent", gp_b_ev, stage="workup", volume=q(2.5, "mL"))
participant(gp_b, "p_nahco3", "nahco3", "workup_agent", gp_b_ev, stage="workup", volume=q(10, "mL", raw_text="10 mL of 5% NaHCO3"))
participant(gp_b, "p_dcm_extract", "dcm", "solvent", gp_b_ev, stage="workup", volume=q(5, "mL", raw_text="2 × 5 mL", repetitions=2))
participant(gp_b, "p_na2so4", "na2so4", "workup_agent", gp_b_ev, stage="workup")
participant(gp_b, "p_silica", "silica_gel", "purification_agent", gp_b_ev, stage="purification")
setf(gp_b, "reaction_scale", "imine (400 μmol)", {"reference_molecule_id": "gp_b_imine", "amount": q(400, "umol")}, gp_b_ev)
setf(gp_b["solvent_conditions"], "composition", "toluene", ["toluene"], gp_b_ev)
setf(gp_b["solvent_conditions"], "total_volume", "1.0 mL", q(1.0, "mL"), gp_b_ev)
add_step(gp_b, "catalyst_mix", gp_b_ev, ["p_cu", "p_pcy3", "p_toluene", "p_tbuoh"], "charge in an Ar-filled glovebox and stir", time=q(5, "min"), atmosphere="argon")
add_step(gp_b, "substrate_charge", gp_b_ev, ["p_allene", "p_imine"], "charge allenamide then imine, seal vial, and remove from glovebox", atmosphere="argon")
add_step(gp_b, "silane_addition", gp_b_ev, ["p_silane"], "add dimethoxymethylsilane")
add_step(gp_b, "reaction_hold", gp_b_ev, ["p_imine", "p_allene", "p_cu", "p_pcy3", "p_tbuoh", "p_silane", "p_toluene"], "stir reaction mixture", temperature="rt", time=q(24, "h"))
add_step(gp_b, "nh4f_quench", gp_b_ev, ["p_nh4f", "p_meoh"], "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"), stage="workup")
add_step(gp_b, "extract", gp_b_ev, ["p_nahco3", "p_dcm_extract"], "add 5% NaHCO3 and extract with DCM", stage="workup")
add_step(gp_b, "dry_concentrate", gp_b_ev, ["p_na2so4"], "dry combined organics with Na2SO4, filter, and concentrate in vacuo", stage="workup")
add_step(gp_b, "chromatography", gp_b_ev, ["p_silica"], "purify by flash chromatography on silica gel", stage="purification")
setf(gp_b["procedure"], "original_text", "General Procedure B for the Synthesis of 18", "General Procedure B for the Synthesis of 18", gp_b_ev)
classify_partial(gp_b, gp_b_ev, ["C-C", "C-H"], "Generic coupling supports net C-C formation and hydrogen gain; atom mapping is unresolved.")


product_details = {
    "18a": (180, 85, "5% E.A. in DCM"),
    "18b": (147, 80, "10% E.A. in DCM"),
    "18c": (153, 89, "10% E.A. in DCM"),
    "18d": (146, 74, "5% E.A. in DCM"),
    "18e": (147, 77, "5% E.A. in DCM"),
}
for product, (mass_mg, yld, eluent) in product_details.items():
    rid = f"{DOC}:scheme4_scope_{product}"
    r = reactions[rid]
    ev = list(dict.fromkeys(r["evidence_ids"] + [E[12]]))
    r["evidence_ids"] = ev
    if f"{DOC}:general_procedure_b_products_18" not in r["general_scheme_ids"]:
        r["general_scheme_ids"].append(f"{DOC}:general_procedure_b_products_18")
    for fact_key in ("experiment_categories", "location", "reaction_name", "transformation_type", "reaction_scale"):
        r[fact_key]["evidence_ids"] = list(dict.fromkeys(r[fact_key]["evidence_ids"] + [E[12]]))
    by_pid = {p["participant_id"]: p for p in r["participants"]}
    # Exact procedural charges from page 12.
    setf(by_pid["p_cu"]["charges"], "mass", "3.6 mg", q(3.6, "mg"), [E[12]])
    setf(by_pid["p_cu"]["charges"], "amount", "20 μmol", q(20, "umol"), [E[12]])
    setf(by_pid["p_ligand"]["charges"], "mass", "7.3 mg", q(7.3, "mg"), [E[12]])
    setf(by_pid["p_ligand"]["charges"], "amount", "26 μmol", q(26, "umol"), [E[12]])
    setf(by_pid["p_silane"]["charges"], "volume", "0.099 mL", q(0.099, "mL"), [E[12]])
    setf(by_pid["p_nh4f"]["charges"], "mass", "200 mg", q(200, "mg"), [E[12]])
    setf(by_pid["p_meoh"]["charges"], "volume", "2.5 mL", q(2.5, "mL"), [E[12]])
    # Source conflict: Scheme 4 says 76 uL; General Procedure B says 76.5 uL.
    alt_scheme = F("76 μL", q(76, "uL"), [E[4]], "original_extraction")
    alt_proc = F("76.5 μL", q(76.5, "uL"), [E[12]], "original_extraction")
    by_pid["p_tbuoh"]["charges"]["volume"] = conflict_fact(
        [alt_scheme, alt_proc], [E[4], E[12]],
        "Scheme 4 footnote and General Procedure B print different t-BuOH volumes; neither was selected.")
    for mid, pid, role, stage, volume in [
        ("nahco3", "p_nahco3", "workup_agent", "workup", q(10, "mL", raw_text="10 mL of 5% NaHCO3")),
        ("dcm", "p_dcm_extract", "solvent", "workup", q(5, "mL", raw_text="2 × 5 mL", repetitions=2)),
        ("na2so4", "p_na2so4", "workup_agent", "workup", None),
        ("silica_gel", "p_silica", "purification_agent", "purification", None),
        ("etoac", "p_etoac", "purification_agent", "purification", None),
        ("dcm", "p_dcm_eluent", "purification_agent", "purification", None),
    ]:
        if pid not in by_pid:
            by_pid[pid] = participant(r, pid, mid, role, [E[12]], stage=stage, volume=volume)
    r["steps"] = []
    add_step(r, "catalyst_mix", [E[12]], ["p_cu", "p_ligand", "p_solvent", "p_tbuoh"], "charge in an Ar-filled glovebox and stir", time=q(5, "min"), atmosphere="argon")
    add_step(r, "substrate_charge", [E[12]], ["p_allene", "p_imine"], "charge allenamide then imine, seal vial, and remove from glovebox", atmosphere="argon")
    add_step(r, "silane_addition", [E[12]], ["p_silane"], "add dimethoxymethylsilane")
    add_step(r, "reaction_hold", [E[12]], ["p_imine", "p_allene", "p_cu", "p_ligand", "p_tbuoh", "p_silane", "p_solvent"], "stir reaction mixture", temperature="rt", time=q(24, "h"))
    add_step(r, "nh4f_quench", [E[12]], ["p_nh4f", "p_meoh"], "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"), stage="workup")
    add_step(r, "extract", [E[12]], ["p_nahco3", "p_dcm_extract"], "add 5% NaHCO3 and extract with DCM", stage="workup")
    add_step(r, "dry_concentrate", [E[12]], ["p_na2so4"], "dry combined organics with Na2SO4, filter, and concentrate in vacuo", stage="workup")
    add_step(r, "chromatography", [E[12]], ["p_silica", "p_etoac", "p_dcm_eluent"], f"purify by silica gel chromatography ({eluent})", stage="purification")
    setf(r["procedure"], "general_procedure_reference", "According to General Procedure B", f"{DOC}:general_procedure_b_products_18", [E[12]])
    setf(r["procedure"], "original_text", f"According to General Procedure B; purification with {eluent} afforded {mass_mg} mg ({yld}%) of {product} as a single diastereomer.",
         f"According to General Procedure B; purification with {eluent} afforded {mass_mg} mg ({yld}%) of {product} as a single diastereomer.", [E[12]])
    # Add product amount and strengthen the existing yield basis without changing the value.
    for o in r["outcomes"]:
        if o["outcome_id"] == "isolated_yield":
            o["value"]["evidence_ids"] = list(dict.fromkeys(o["value"]["evidence_ids"] + [E[12]]))
            setf(o, "measurement_basis", f"isolated after silica gel chromatography ({eluent})",
                 f"isolated after silica gel chromatography ({eluent})", [E[12]])
    if not any(o["outcome_id"] == "product_amount" for o in r["outcomes"]):
        add_outcome(r, "product_amount", "product_amount", [product], [E[12]], f"{mass_mg} mg", q(mass_mg, "mg"),
                    basis_raw=f"isolated after silica gel chromatography ({eluent})",
                    basis=f"isolated after silica gel chromatography ({eluent})", product_role="main product")
    note(r, "Scheme 4 reports 76 μL t-BuOH, while General Procedure B reports 76.5 μL; the charge is stored as a conflict with both source values.", [E[4], E[12]])
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    if "Scheme 4 and General Procedure B report t-BuOH as 76 μL and 76.5 μL, respectively." not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append("Scheme 4 and General Procedure B report t-BuOH as 76 μL and 76.5 μL, respectively.")
    r["completeness"] = "partial"


# Scheme 7: two concrete catalyst-controlled experiments plus the separate 29 -> 30 comparison.
def scheme7_record(product, with_tbuoh, detailed_yield, scheme_yield, er_ratio, hplc):
    ev = [E[9], E[10], E[17]]
    r = make_reaction(f"scheme7_synthesis_{product}", ev, "specific_experiment",
                      ["mechanistic_study"], f"Pages 9--10, Scheme 7; contextual experimental text page 17",
                      f"synthesis of {product} from achiral allenamide 15b under (S,S)-Ph-BPE catalysis",
                      "catalyst-controlled imine–allenamide reductive coupling",
                      [f"{DOC}:scheme2d_current_generic_reductive_coupling"],
                      ["Scheme 7 prints 5 mol% Cu(OAc)2 and 6 mol% (S,S)-Ph-BPE, whereas the detailed 10 μmol/250 μmol charges correspond to 4 mol% each."])
    participant(r, "p_9a", "9a", "reactant", ev, amount=q(250, "umol"))
    participant(r, "p_15b", "15b", "reactant", ev, mass=q(37.5, "mg"), amount=q(300, "umol"))
    pcu = participant(r, "p_cu", "cu_oac2", "catalyst", ev, mass=q(1.8, "mg"), amount=q(10, "umol"), ref="9a")
    plig = participant(r, "p_ligand", "ph_bpe_ss", "ligand", ev, mass=q(5.1, "mg"), amount=q(10, "umol"), ref="9a")
    pcu["charges"]["mol_percent"] = conflict_fact(
        [F("5 mol%", q(5, "mol%"), [E[10]]), F("10 μmol relative to 250 μmol 9a", q(4, "mol%"), [E[17]], "inference", "10/250 × 100 = 4 mol%", verification="inferred")],
        [E[10], E[17]], "Scheme loading and detailed charge-derived loading conflict.")
    plig["charges"]["mol_percent"] = conflict_fact(
        [F("6 mol%", q(6, "mol%"), [E[10]]), F("10 μmol relative to 250 μmol 9a", q(4, "mol%"), [E[17]], "inference", "10/250 × 100 = 4 mol%", verification="inferred")],
        [E[10], E[17]], "Scheme loading and detailed charge-derived loading conflict.")
    participant(r, "p_silane", "silane", "reagent", ev, volume=q(0.061, "mL"), equiv=q(2, "equiv"), ref="9a")
    if with_tbuoh:
        participant(r, "p_tbuoh", "tbuoh", "additive", ev, volume=q(26.3, "uL"), amount=q(275, "umol"))
    participant(r, "p_toluene", "toluene", "solvent", ev, volume=q(2.5 if with_tbuoh else 0.5, "mL"))
    participant(r, "p_product", product, "product", ev, product_role="main product")
    participant(r, "p_nh4f", "nh4f", "workup_agent", ev, stage="workup", mass=q(100, "mg"))
    participant(r, "p_meoh", "meoh", "solvent", ev, stage="workup", volume=q(1.5, "mL"))
    participant(r, "p_nahco3", "nahco3", "workup_agent", ev, stage="workup", volume=q(5, "mL", raw_text="5 mL of 5% NaHCO3"))
    participant(r, "p_dcm", "dcm", "solvent", ev, stage="workup", volume=q(3, "mL", raw_text="2 × 3 mL", repetitions=2))
    participant(r, "p_na2so4", "na2so4", "workup_agent", ev, stage="workup")
    participant(r, "p_silica", "silica_gel", "purification_agent", ev, stage="purification")
    participant(r, "p_etoac", "etoac", "purification_agent", ev, stage="purification")
    participant(r, "p_hexanes", "hexanes", "purification_agent", ev, stage="purification")
    setf(r, "reaction_scale", "imine 9a (250 μmol)", {"reference_molecule_id": "9a", "amount": q(250, "umol")}, ev)
    setf(r["solvent_conditions"], "composition", "toluene", ["toluene"], ev)
    setf(r["solvent_conditions"], "total_volume", f"{2.5 if with_tbuoh else 0.5} mL", q(2.5 if with_tbuoh else 0.5, "mL"), [E[17]])
    mix_pids = ["p_cu", "p_ligand", "p_toluene"] + (["p_tbuoh"] if with_tbuoh else [])
    add_step(r, "catalyst_mix", [E[17]], mix_pids, "charge in an Ar-filled glovebox and stir", time=q(5, "min"), atmosphere="argon")
    add_step(r, "substrate_charge", [E[17]], ["p_15b", "p_9a"], "charge allenamide then imine, seal vial, and remove from glovebox", atmosphere="argon")
    add_step(r, "silane_addition", [E[17]], ["p_silane"], "add dimethoxymethylsilane")
    rxn_pids = ["p_9a", "p_15b", "p_cu", "p_ligand", "p_silane", "p_toluene"] + (["p_tbuoh"] if with_tbuoh else [])
    add_step(r, "reaction_hold", [E[17]], rxn_pids, "stir reaction mixture", temperature="rt", time=q(24, "h"))
    add_step(r, "nh4f_quench", [E[17]], ["p_nh4f", "p_meoh"], "add NH4F and MeOH and agitate", temperature="rt", time=q(30, "min"), stage="workup")
    add_step(r, "extract", [E[17]], ["p_nahco3", "p_dcm"], "add 5% NaHCO3 and extract with DCM", stage="workup")
    add_step(r, "dry_concentrate", [E[17]], ["p_na2so4"], "dry combined organics with Na2SO4, filter, and concentrate in vacuo", stage="workup")
    chromatography = "25% EtOAc/hexanes" if product == "29" else "50% EtOAc/hexanes"
    add_step(r, "chromatography", [E[17]], ["p_silica", "p_etoac", "p_hexanes"], f"flash chromatography on silica gel ({chromatography})", stage="purification")
    yield_value = F(f"{detailed_yield}%", q(detailed_yield, "%"), [E[10], E[17]])
    if scheme_yield != detailed_yield:
        yield_value = conflict_fact(
            [F(f"{scheme_yield}%", q(scheme_yield, "%"), [E[10]]),
             F(f"{detailed_yield}%", q(detailed_yield, "%"), [E[17]])],
            [E[10], E[17]], "Scheme 7 and the detailed experimental paragraph report different isolated yields.")
        r["review"]["unresolved"].append(f"{product}: Scheme 7 reports {scheme_yield}% but the detailed experimental paragraph reports {detailed_yield}% yield.")
    product_mass = 78 if product == "29" else 68
    add_outcome(r, "isolated_yield", "yield", [product], ev, None, None,
                value_fact=yield_value, basis_raw=f"isolated after silica gel chromatography ({chromatography})",
                basis=f"isolated after silica gel chromatography ({chromatography})",
                yield_type="isolated yield", product_role="main product")
    add_outcome(r, "product_amount", "product_amount", [product], [E[17]], f"{product_mass} mg", q(product_mass, "mg"),
                basis_raw=f"isolated after silica gel chromatography ({chromatography})",
                basis=f"isolated after silica gel chromatography ({chromatography})", product_role="main product")
    add_outcome(r, "diastereoselectivity", "selectivity", [product], [E[9], E[10], E[17]], "single diastereomer",
                {"metric": "dr", "reported": "single diastereomer", "major_configuration": "relative configuration shown in Scheme 7", "minor_configuration": "not observed"},
                basis_raw="single diastereomer; analytical method not reported for dr",
                basis="dr method unreported", product_role="main product",
                scope={"metric": "dr", "order": "major:minor", "major_configuration": "relative configuration shown", "minor_configuration": "unassigned"})
    add_outcome(r, "enantiomeric_ratio", "selectivity", [product], [E[9], E[10], E[17]], er_ratio,
                {"metric": "er", "ratio": er_ratio, "order": "first-eluting:second-eluting enantiomer; configurations unassigned"},
                basis_raw=f"chiral HPLC analysis ({hplc})", basis=f"chiral HPLC: {hplc}", product_role="main product",
                scope={"metric": "er", "order": "first-eluting:second-eluting", "major_configuration": "unassigned", "minor_configuration": "unassigned"})
    note(r, "Scheme 7 catalyst/ligand loadings conflict with the detailed charges; both are retained without choosing a corrected value.", [E[10], E[17]])
    classify_partial(r, ev, ["C-C", "C-H"] + (["C-N"] if product == "30" else []),
                     "Net coupling changes are inferred from the source drawings; complete atom mapping and rearrangement cleavages are unresolved.")
    return r


scheme7_record("29", True, 69, 69, "57:43", "Chiracel AD-3, 85:15 heptane/isopropanol, 1.50 mL/min, 254 nm")
scheme7_record("30", False, 60, 62, "80:20", "Chiracel AD-3, 90:10 heptane/isopropanol, 1.00 mL/min, 220 nm")

ev_2930 = [E[9]]
r = make_reaction("scheme7_29_to_30_nbuli", ev_2930, "specific_experiment", ["mechanistic_study"],
                  "Page 9, Scheme 7 discussion", "separate conversion of 29 to 30 using n-BuLi",
                  "base-induced carbamate migration",
                  unresolved=["The source statement gives no scale, solvent, n-BuLi dose, temperature, time, workup, product amount, or yield for this comparison experiment."])
participant(r, "p_29", "29", "reactant", ev_2930)
participant(r, "p_nbuli", "nbuli", "reagent", ev_2930)
participant(r, "p_30", "30", "product", ev_2930, product_role="comparison product")
add_step(r, "conversion", ev_2930, ["p_29", "p_nbuli"], "convert 29 to 30 using n-BuLi")
add_outcome(r, "stereochemical_correspondence", "other", ["30"], ev_2930,
            "confirmed that the same relative stereochemistry was formed in both reactions",
            "same relative stereochemistry as the independently prepared 30",
            basis_raw="comparison method not reported", basis="method unreported",
            product_role="comparison product", scope={"metric": "relative stereochemical correspondence"})
classify_partial(r, ev_2930, ["C-N"], "Urea C-N closure is inferred; detailed conditions and complete bond changes are unreported.")


# Enrich the existing proposed-mechanism record with Figures 2--5, without
# turning computational structures/energies into synthetic experiments.
mech = reactions[f"{DOC}:scheme3_proposed_catalytic_cycle"]
mech_ev = [E[7], E[8], E[9], E[10]]
mech["evidence_ids"] = list(dict.fromkeys(mech["evidence_ids"] + mech_ev))
note(mech, "Figures 2--5 report DFT hydrocupration/isomerization/C-C bond-formation pathways, distortion-interaction analysis, NCI analysis, and s-trans/s-cis conformational energetics; these are computational support for the proposed mechanism, not isolated synthetic experiments.", mech_ev, "original_extraction")
mech["review"]["status"] = "needs_review"
mech["review"]["method"] = "model"


# Non-synthetic computations and characterization are retained as context.
bundle["source_context"].extend([
    {"document_id": DOC, "evidence_ids": [E[7], I[7]], "kind": "computational_mechanistic_analysis",
     "text": "Figure 2 gives calculated hydrocupration and allylcopper isomerization pathways/transition-state structures and relative free energies; no isolated reaction experiment is created."},
    {"document_id": DOC, "evidence_ids": [E[8], I[8]], "kind": "computational_mechanistic_analysis",
     "text": "Figure 3 gives the calculated C-C bond-forming/rearrangement pathway and diastereomeric transition states; it supports the proposed mechanism rather than constituting synthetic rows."},
    {"document_id": DOC, "evidence_ids": [E[9], I[9]], "kind": "computational_mechanistic_analysis",
     "text": "Figure 4 reports distortion-interaction and noncovalent-interaction analyses. Scheme 7 discussion supplies the separate catalyst-controlled experiments recorded here."},
    {"document_id": DOC, "evidence_ids": [E[10], I[10]], "kind": "computational_and_general_context",
     "text": "Figure 5 compares s-trans/s-cis conformations. The conclusion and General analytical methods are context; Scheme 7 and General Procedure A are separately mapped to records."},
    {"document_id": DOC, "evidence_ids": [E[10], E[11], E[12]], "kind": "characterization_data",
     "text": "Melting points, Rf values, NMR chemical shifts/couplings, and HRMS calculated/found m/z blocks for 9d, 9de, 9e, 9f, 9g, 9k, 9n, 9o, 9pa, 9pd, 9pe, 9q, 9r, 9s, 9t, and 18a-18e are characterization, not independent synthesis experiments."},
    {"document_id": DOC, "evidence_ids": [E[10]], "kind": "literature_preparation_context",
     "text": "The General section states that imines 9a, 9b, 9c, 9ee, 9h, 9i, 9j, 9l, 9m, 9pb, 9pc, and 9u were synthesized as described in cited literature; these are not claimed as newly detailed experiments on pages 7-12."},
])


# Required exact six-page coverage inventory.
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted",
         "regions": [
             {"region": "Figure 2 and caption", "finding": "DFT hydrocupration/isomerization pathways and relative energies; computational mechanistic analysis, not a synthetic experiment.", "evidence_ids": [E[7], I[7]]},
             {"region": "lower two-column mechanistic text", "finding": "Author interpretation of imine coordination and reversible diastereoselective C-C bond formation; linked to the existing proposed-mechanism record.", "evidence_ids": [E[7], I[7]]},
         ], "record_ids": [f"{DOC}:scheme3_proposed_catalytic_cycle"],
         "unresolved": ["Transition-state geometries and energies are not atom-mapped synthetic records; no SI coordinates were associated."]},
        {"page": 8, "state": "attempted",
         "regions": [
             {"region": "Figure 3 and caption", "finding": "Calculated pathway, competing diastereomeric transition states, and proposed branches to 18/19; retained as support for the proposed mechanism.", "evidence_ids": [E[8], I[8]]},
             {"region": "lower discussion columns", "finding": "Distortion-interaction/NCI interpretation and product-selectivity rationale; no additional performed synthesis row.", "evidence_ids": [E[8], I[8]]},
         ], "record_ids": [f"{DOC}:scheme3_proposed_catalytic_cycle"],
         "unresolved": ["Computational conformers and proposed intermediates are not treated as isolated compounds or experiments."]},
        {"page": 9, "state": "attempted",
         "regions": [
             {"region": "Figure 4", "finding": "Distortion-interaction and noncovalent-interaction analyses; computational source context.", "evidence_ids": [E[9], I[9]]},
             {"region": "lower-left mechanistic discussion", "finding": "Calculated t-BuOH/no-t-BuOH product-selectivity branches; linked to the proposed mechanism.", "evidence_ids": [E[9], I[9]]},
             {"region": "lower-right Scheme 7 lead-in", "finding": "Two catalyst-controlled reactions (to 29 and 30) and a separate n-BuLi conversion of 29 to 30 are identified as mechanistic-study experiments.", "evidence_ids": [E[9], I[9]]},
         ], "record_ids": [f"{DOC}:scheme3_proposed_catalytic_cycle", f"{DOC}:scheme7_synthesis_29", f"{DOC}:scheme7_synthesis_30", f"{DOC}:scheme7_29_to_30_nbuli"],
         "unresolved": ["The 29→30 comparison statement gives no quantities, conditions, yield, or comparison method."]},
        {"page": 10, "state": "attempted",
         "regions": [
             {"region": "Figure 5", "finding": "Computed s-trans/s-cis enamide conformational comparison; source context, not synthesis.", "evidence_ids": [E[10], I[10]]},
             {"region": "Scheme 7", "finding": "Concrete catalyst-controlled preparations of 29 (69%, 57:43 er) and 30 (scheme: 62%, 80:20 er); mapped to independent mechanistic-study records.", "evidence_ids": [E[10], I[10]]},
             {"region": "Conclusion and General experimental methods", "finding": "Narrative/general characterization context; no extra reaction identity.", "evidence_ids": [E[10], I[10]]},
             {"region": "General Procedure A and preparations of 9d/9de", "finding": "Generic imine procedure plus two product-specific 6 mmol preparations with isolated amounts/yields.", "evidence_ids": [E[10], I[10]]},
         ], "record_ids": [f"{DOC}:scheme3_proposed_catalytic_cycle", f"{DOC}:scheme7_synthesis_29", f"{DOC}:scheme7_synthesis_30", f"{DOC}:general_procedure_a_imines", f"{DOC}:imine_synthesis_9d", f"{DOC}:imine_synthesis_9de"],
         "unresolved": ["Scheme 7 catalyst/ligand loadings conflict with detailed page-17 charges; 30 yield is 62% in Scheme 7 versus 60% in the detailed procedure.", "No SI was associated."]},
        {"page": 11, "state": "attempted",
         "regions": [
             {"region": "left/right columns, imine preparation blocks", "finding": "Product-specific General Procedure A preparations for 9e, 9f, 9g, 9k, 9n, 9o, 9pa, 9pd, and 9pe; 9q begins and continues on page 12.", "evidence_ids": [E[11], I[11]]},
             {"region": "NMR/19F NMR/HRMS/melting-point text", "finding": "Characterization retained as source context rather than separate synthesis experiments.", "evidence_ids": [E[11], I[11]]},
         ], "record_ids": [f"{DOC}:imine_synthesis_{x}" for x in ["9e", "9f", "9g", "9k", "9n", "9o", "9pa", "9pd", "9pe", "9q"]],
         "unresolved": ["For 9o and 9pd, the paper reports both yield and sub-100% purity but does not state whether yields are purity-corrected."]},
        {"page": 12, "state": "attempted",
         "regions": [
             {"region": "upper-left continuation and imine blocks", "finding": "Completion of 9q plus preparations of 9r, 9s, and 9t under General Procedure A.", "evidence_ids": [E[12], I[12]]},
             {"region": "General Procedure B", "finding": "Detailed ordered Cu-catalyzed coupling, quench, extraction, drying, and chromatography procedure for products 18.", "evidence_ids": [E[12], I[12]]},
             {"region": "product blocks 18a-18e", "finding": "Detailed product masses, yields, single-diastereomer statements, and purification eluents merged into the existing Scheme 4 scope experiment identities.", "evidence_ids": [E[12], I[12]]},
             {"region": "NMR/19F NMR/HRMS/Rf text", "finding": "Characterization retained as source context and not promoted to extra synthesis records.", "evidence_ids": [E[12], I[12]]},
         ], "record_ids": [f"{DOC}:imine_synthesis_{x}" for x in ["9q", "9r", "9s", "9t"]] + [f"{DOC}:general_procedure_b_products_18"] + [f"{DOC}:scheme4_scope_18{x}" for x in "abcde"],
         "unresolved": ["9q and 9s yields are reported with 94% and 91% purity, respectively, without a stated purity-correction basis.", "Scheme 4 reports 76 μL t-BuOH while General Procedure B reports 76.5 μL; both values are retained as a conflict."]},
    ],
}


# All machine records, including pre-existing same-document records, remain needs_review.
for rr in bundle["reactions"]:
    rr["review"]["status"] = "needs_review"
    if rr["review"]["method"] == "none":
        rr["review"]["method"] = "model"

write(ROOT / "runs/current/submitted.json", bundle)
write(ROOT / "runs/current/coverage/pages_0007_0012.json", coverage)
write(ROOT / "runs/current/raw_model/stereo_diagnostics_pages_0007_0012.json", stereo_diagnostics)
new_ids = [
    rid for rid in reactions
    if rid.startswith(f"{DOC}:imine_synthesis_")
    or rid in {f"{DOC}:general_procedure_a_imines", f"{DOC}:general_procedure_b_products_18",
               f"{DOC}:scheme7_synthesis_29", f"{DOC}:scheme7_synthesis_30",
               f"{DOC}:scheme7_29_to_30_nbuli"}
]
updated_ids = [f"{DOC}:scheme3_proposed_catalytic_cycle"] + [f"{DOC}:scheme4_scope_18{x}" for x in "abcde"]
write(ROOT / "runs/current/raw_model/pages_0007_0012_inventory.json", {
    "document_id": DOC,
    "pages": [7, 8, 9, 10, 11, 12],
    "new_reaction_ids": new_ids,
    "updated_reaction_ids": updated_ids,
    "new_reaction_count": len(new_ids),
    "updated_reaction_count": len(updated_ids),
    "notes": [
        "All six supplied original page images were viewed individually at original detail.",
        "Page 17 extracted text was used only for detailed Scheme 7 charges/workup/outcomes; page 17 is not claimed as semantically covered.",
        "No SI was associated or inferred.",
        "Independent NMR/HRMS/Rf/melting-point blocks remain source_context.",
    ],
})
print(json.dumps({"cumulative_reactions": len(bundle["reactions"]),
                  "cumulative_molecules": len(bundle["molecules"]),
                  "new_reactions": len(new_ids), "updated_reactions": len(updated_ids)}, ensure_ascii=False))
