#!/usr/bin/env python3
"""Build the v3.1 machine extraction for PDF pages 1--6.

This is an evidence-preserving construction script, not a chemistry validator.
All records remain needs_review and all structure/bond inferences are labelled.
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
    empty,
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
TITLE = "Access to Chiral Diamine Derivatives through Stereoselective Cu-Catalyzed Reductive Coupling of Imines and Allenamides"


def q(number, unit, qualifier=None, raw_text=None):
    result = {"number": number, "unit": unit}
    if qualifier:
        result["qualifier"] = qualifier
    if raw_text:
        result["raw_text"] = raw_text
    return result


def F(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, verification=None, alternatives=()):
    return fact(raw, normalized, evidence, method, status, verification, derivation, alternatives)


def setf(container, key, raw, normalized, evidence, method="original_extraction", derivation=None,
         status=None, verification=None, alternatives=()):
    container[key] = F(raw, normalized, evidence, method, derivation, status, verification, alternatives)


bundle_path = ROOT / "runs/current/submitted.json"
dataset_path = ROOT / "runs/current/dataset.json"
dataset = json.loads(dataset_path.read_text())
dataset["documents"][0]["title"] = TITLE
write(dataset_path, dataset)
bundle = empty(dataset["documents"])

molecules = {}
stereo_diagnostics = {}


def add_molecule(mid, label, evidence, name=None, chirality="unresolved", geometry="unresolved",
                 assignment=None, mixture=None, smiles=None, note=None):
    if mid in molecules:
        return mid
    ev = list(dict.fromkeys(evidence))
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    setf(m["occurrences"][0], "compound_label", label, label, ev)
    setf(m["names"], "original", name or label, name or label, ev)
    setf(m["names"], "english", name or label, name or label, ev)
    setf(m["names"], "abbreviations", None, None, ev, status="unreported", verification="missing")
    setf(m["names"], "aliases", None, None, ev, status="unreported", verification="missing")
    setf(m, "cas_numbers", None, None, ev, status="unreported", verification="missing")
    setf(m["stereochemistry"], "chirality_status", chirality, chirality, ev)
    setf(m["stereochemistry"], "geometry_status", geometry, geometry, ev)
    rep = "not_encoded"
    if smiles:
        m["representation_type"] = "smiles"
        setf(m, "smiles", label, smiles, ev, method="structure_derivation",
             derivation="Machine candidate derived from the source chemical name/label; RDKit parsing does not verify the paper drawing.",
             verification="inferred")
        diag = stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        setf(m["stereochemistry"], "isomeric_smiles", label, diag["canonical_isomeric_smiles"], ev,
             method="structure_derivation", derivation="RDKit canonicalization of the unverified machine candidate.", verification="inferred")
        setf(m["stereochemistry"], "identity_key", label, diag["stereo_graph_key"], ev,
             method="structure_derivation", derivation="SHA-256 key of the canonical isomeric graph candidate.", verification="inferred")
        rep = "partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded"
    else:
        m["representation_type"] = "unresolved"
        setf(m, "smiles", None, None, ev, status="unreported", verification="missing")
    setf(m["stereochemistry"], "representation_status", rep, rep, ev,
         method="inference", derivation="Source drawing/label was inspected, but a complete source-verified isomeric graph was not established." if not smiles else "Status follows the machine graph diagnostic; source structure remains unverified.", verification="inferred")
    if assignment is not None:
        setf(m["stereochemistry"], "assignment", assignment, assignment, ev)
        setf(m["stereochemistry"], "atom_centers", assignment, assignment, ev)
    if mixture is not None:
        setf(m["stereochemistry"], "mixture", mixture, mixture, ev)
    setf(m, "identity_note", note or "Identity follows the source label/drawing; no CAS number was supplied.",
         note or "Identity follows the source label/drawing; no CAS number was supplied.", ev,
         method="inference", derivation="Machine identity note from inspected source evidence.", verification="inferred")
    molecules[mid] = m
    return mid


# Reusable exact small-molecule candidates. Complex substrates/products remain unencoded.
simple = [
    ("cu_oac2", "Cu(OAc)2", "copper(II) acetate", "[Cu+2].CC(=O)[O-].CC(=O)[O-]"),
    ("pcy3", "PCy3", "tricyclohexylphosphine", "P(C1CCCCC1)(C1CCCCC1)C1CCCCC1"),
    ("silane", "Me(MeO)2SiH", "dimethoxymethylsilane", "CO[SiH](C)OC"),
    ("tbuoh", "t-BuOH", "tert-butanol", "CC(C)(C)O"),
    ("toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("mtbe", "MTBE", "methyl tert-butyl ether", "COC(C)(C)C"),
    ("dioxane", "dioxane", "1,4-dioxane", "C1COCCO1"),
    ("dcm", "CH2Cl2", "dichloromethane", "ClCCl"),
    ("thf", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("nh4f", "NH4F", "ammonium fluoride", "[NH4+].[F-]"),
    ("meoh", "MeOH", "methanol", "CO"),
    ("dimethyl_fumarate", "dimethyl fumarate", "dimethyl fumarate", "COC(=O)/C=C/C(=O)OC"),
    ("hcl", "HCl", "hydrogen chloride", "Cl"),
    ("ether", "Et2O", "diethyl ether", "CCOCC"),
    ("mscl", "MsCl", "methanesulfonyl chloride", "CS(=O)(=O)Cl"),
    ("et3n", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("kotbu", "KOtBu", "potassium tert-butoxide", "[K+].CC(C)(C)[O-]"),
    ("h2so4", "H2SO4", "sulfuric acid", "OS(=O)(=O)O"),
    ("allyl_br", "allyl bromide", "allyl bromide", "C=CCBr"),
    ("k2co3", "K2CO3", "potassium carbonate", "[K+].[K+].[O-]C([O-])=O"),
    ("tbai", "TBAI", "tetrabutylammonium iodide", "CCCC[N+](CCCC)(CCCC)CCCC.[I-]"),
    ("mecn", "CH3CN", "acetonitrile", "CC#N"),
]
for mid, label, name, smiles in simple:
    add_molecule(mid, label, [E[3], E[4], E[5], E[6], E[17], E[18]], name=name,
                 chirality="achiral", geometry="not_applicable", smiles=smiles)

for mid, label in [
    ("dcpe", "dcpe"), ("p_adam3", "P(adam)3"), ("xphos", "XPhos"),
    ("p_nme2_3", "P(NMe2)3"), ("p_oet3", "P(OEt)3"),
    ("ph_o2_pnme2", "(PhO)2PNMe2"), ("simes", "SIMes"),
    ("nbuli", "n-BuLi"), ("hg2", "Hoveyda–Grubbs second-generation catalyst"),
]:
    add_molecule(mid, label, [E[3], E[6], E[17], E[18]], chirality="unresolved", geometry="not_applicable",
                 note="Named reagent/catalyst is source identified; a structure candidate was not needed for this page batch.")

# Central substrates and products.
add_molecule("15a", "15a", [E[3], E[4], E[5], E[6], E[12]],
             name="chiral Evans-oxazolidinone-derived allenamide 15a", chirality="specified",
             geometry="not_applicable", assignment="(S) auxiliary stereocenter; terminal allene has no assigned axial stereochemistry",
             note="Chiral allenamide drawing inspected; complete stereochemical SMILES intentionally not encoded.")
add_molecule("15b", "15b", [E[3]], name="achiral allenamide 15b", chirality="achiral", geometry="not_applicable")

imine_labels = set(["9a", "9c", "9p_aliphatic"])
imine_labels.update(f"9{x}" for x in list("bcdefghijklmnopqrstu"))
imine_labels.update(["9pa", "9pb", "9pc", "9pd", "9pe", "9de", "9ee"])
for label in sorted(imine_labels):
    ev = [E[3], E[4], E[5], E[10], E[11], E[12]]
    geometry = "specified" if label != "9p_aliphatic" else "unspecified"
    assignment = "E imine geometry reported in the compound name/drawing" if geometry == "specified" else None
    add_molecule(label, label, ev, name=("aliphatic aldimine (Ar1 = Me)" if label == "9p_aliphatic" else f"imine {label}"),
                 chirality="achiral", geometry=geometry, assignment=assignment,
                 note="Imine identity follows source label/substituent drawing; complete SMILES was not encoded.")

for label in [f"18{x}" for x in list("abcdefghijklmnopqrstu")]:
    add_molecule(label, label, [E[4], E[5], E[12], E[13], E[14]], name=f"branched diamine synthon {label}",
                 chirality="specified", geometry="not_applicable", assignment="(S,S,S) product; source reports a single diastereomer",
                 note="Complex source drawing inspected; absolute assignment is source-reported, but no full isomeric SMILES was encoded.")
for label in ["19a", "19b", "19de", "19ee", "19j", "19k", "19l", "19m", "19n", "19s", "19t", "19pe",
              "19pa", "19pb", "19pc", "19pd", "19i", "19u"]:
    add_molecule(label, label, [E[3], E[5], E[6], E[15], E[16]], name=f"urea product {label}",
                 chirality="specified", geometry="not_applicable", assignment="(S,S,S) product; source reports a single diastereomer",
                 note="Complex source drawing inspected; full isomeric SMILES was not encoded.")
add_molecule("18c_hcl", "18c·HCl", [E[5], E[6]], name="hydrochloride salt of 18c",
             chirality="specified", geometry="not_applicable", assignment="stereochemistry assigned by single-crystal X-ray analysis")
for label in ["25", "26", "27"]:
    add_molecule(label, label, [E[6], E[17], E[18]], name=f"synthetic application product {label}",
                 chirality="partially_specified", geometry="unspecified",
                 assignment="stereochemistry retained from the single-stereoisomer precursor; the page does not print a full absolute descriptor",
                 note="Product drawing is source identified; stereochemical graph remains not encoded.")


def generic_molecule(mid, label, page):
    return add_molecule(mid, label, [E[page]], name=label, chirality="unresolved", geometry="unresolved",
                        note="Generic R-group/background structure from the scheme; no concrete molecule or SMILES is asserted.")


def make_reaction(suffix, evidence, record_type="specific_experiment", categories=None, location=None,
                  reaction_name=None, transformation=None, general_ids=None, unresolved=None):
    rid = f"{DOC}:{suffix}"
    ev = list(dict.fromkeys(evidence))
    r = new_reaction(rid, DOC, ev)
    r["record_type"] = record_type
    setf(r, "experiment_categories", ", ".join(categories or []), categories or [], ev)
    setf(r, "location", location, location, ev)
    setf(r, "reaction_name", reaction_name, reaction_name, ev)
    setf(r, "transformation_type", transformation, transformation, ev,
         method="inference" if transformation else "original_extraction",
         derivation="Transformation description inferred from the inspected reaction arrow." if transformation else None,
         verification="inferred" if transformation else None)
    r["general_scheme_ids"] = list(general_ids or [])
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    return r


def participant(r, pid, mid, role, ev, stage="reaction", mass=None, amount=None, volume=None,
                equiv=None, molpct=None, ref=None, product_role=None):
    p = new_participant(pid, mid, role, ev, stage)
    if product_role:
        setf(p, "product_role", product_role, product_role, ev)
    vals = {"mass": mass, "amount": amount, "volume": volume, "equivalents": equiv, "mol_percent": molpct}
    for key, value in vals.items():
        if value is not None:
            raw = value.get("raw_text") if isinstance(value, dict) else str(value)
            setf(p["charges"], key, raw or value, value, ev)
    if ref and (equiv is not None or molpct is not None):
        setf(p, "equivalence_reference", ref, ref, ev)
    r["participants"].append(p)
    return p


def add_step(r, sid, ev, pids, operation, temperature=None, time=None, stage="reaction", atmosphere=None):
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


def add_outcome(r, oid, kind, targets, ev, raw, normalized, basis_raw=None, basis=None,
                yield_type=None, product_role=None, scope=None):
    o = new_outcome(oid, kind, targets, ev)
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
    categories = {"C-H": 1, "C-C": 2, "C-O": 3, "C-N": 5, "DB-C-C": 13, "N-H": 29}
    changes = []
    labels = []
    atoms = {"C-H": ("C", "H", 1, "hydrogen_gain"), "C-C": ("C", "C", 1, "formed"),
             "C-O": ("C", "O", 1, "formed"), "C-N": ("C", "N", 1, "formed"), "DB-C-C": ("C", "C", 2, "formed"),
             "N-H": ("N", "H", 1, "hydrogen_gain")}
    for idx, code in enumerate(codes, 1):
        a, b, order, kind = atoms[code]
        cid = f"b{idx}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
              "before_order": 0, "after_order": order, "change_type": kind}
        if kind == "hydrogen_gain":
            ch["hydrogen_count"] = 1
        changes.append(F(code, ch, ev, method="inference",
                         derivation="Net bond-change class inferred from the source reaction drawing; atom mapping is unresolved.", verification="inferred"))
        existing = next((item for item in labels if item["code"] == code), None)
        if existing:
            existing["supporting_change_ids"].append(cid)
        else:
            labels.append({"category_number": categories[code], "code": code, "supporting_change_ids": [cid]})
    r["classification"]["bond_changes"] = changes
    setf(r["classification"], "labels", ", ".join(codes), labels, ev, method="inference",
         derivation="Labels cover only visually supported net formed/promoted bonds.", verification="inferred")
    setf(r["classification"], "coverage", "partial", "partial", ev, method="inference",
         derivation=limitation, verification="inferred")
    r["classification"]["limitations"] = [limitation, "No atom-mapped reaction SMILES was asserted."]


def note(r, text, ev, method="inference"):
    r["notes"].append(F(text, text, ev, method=method,
                        derivation="Extraction/review note tied to the cited source region." if method == "inference" else None,
                        verification="inferred" if method == "inference" else None))


def add_generic_participants(r, mids, page):
    for i, (mid, role) in enumerate(mids, 1):
        participant(r, f"p{i}", mid, role, [E[page]])


reactions = []

# Page 2, Scheme 1: conceptual prior-art arrows, retained as literature examples rather than present-work experiments.
scheme1_specs = [
    ("scheme1a_direct_diamination", "2", "1", "direct catalytic 1,2-diamination", ["C-N", "C-N"]),
    ("scheme1a_aminohydroxylation", "3", "4", "aminohydroxylation", ["C-N", "C-O"]),
    ("scheme1a_amine_substitution", "4", "1", "alcohol activation and amine substitution", ["C-N"]),
    ("scheme1a_aziridination", "3", "5", "aziridination", ["C-N", "C-N"]),
    ("scheme1a_aziridine_opening", "5", "1", "aziridine opening with an amine", ["C-N"]),
    ("scheme1b_aza_pinacol", "6", "1", "aza-pinacol coupling", ["C-C"]),
    ("scheme1b_nitro_mannich", "7 + 6", "1", "nitro-Mannich/reduction route", ["C-C"]),
    ("scheme1b_glycine_mannich", "8 + 6", "1", "glycine-Mannich/decarboxylation route", ["C-C"]),
]
for suffix, src, dst, name, _codes in scheme1_specs:
    src_mid = generic_molecule(f"generic_{suffix}_src", src, 2)
    dst_mid = generic_molecule(f"generic_{suffix}_dst", dst, 2)
    r = make_reaction(suffix, [E[2]], "literature_example", ["other"], "Page 2, Scheme 1",
                      name, name, unresolved=["Conceptual prior-art arrow; no concrete experimental conditions are supplied on this page."])
    add_generic_participants(r, [(src_mid, "reactant"), (dst_mid, "product")], 2)
    classify_partial(r, [E[2]], _codes,
                     "Generic source scheme does not support atom-level mapping or complete bond accounting.")
    reactions.append(r)

# Page 3, Scheme 2 and Scheme 3.
for suffix, rtype, name in [
    ("scheme2a_classical_imine_allylation", "literature_example", "classical allylation of imines"),
    ("scheme2b_proposed_aminoallylation", "generic_scheme", "proposed amino-substituted allyl strategy"),
    ("scheme2c_prior_carbonyl_reductive_coupling", "literature_example", "prior reductive coupling of carbonyls and allenamides"),
    ("scheme2d_current_generic_reductive_coupling", "generic_scheme", "Cu-catalyzed reductive coupling of imines and chiral allenamide 15a"),
]:
    src = generic_molecule(f"{suffix}_reactants", "reactants depicted in Scheme 2", 3)
    dst = generic_molecule(f"{suffix}_product", "product depicted in Scheme 2", 3)
    r = make_reaction(suffix, [E[3]], rtype, ["other"], "Page 3, Scheme 2", name, name,
                      unresolved=["Generic/prior-art drawing; no concrete row-level experiment is asserted."])
    add_generic_participants(r, [(src, "reactant"), (dst, "product")], 3)
    classify_partial(r, [E[3]], ["C-C"], "Only the newly connected carbon skeleton is recoverable from the generic arrow.")
    reactions.append(r)
CURRENT_SCHEME = f"{DOC}:scheme2d_current_generic_reductive_coupling"

mech = make_reaction("scheme3_proposed_catalytic_cycle", [E[3], E[4]], "proposed_mechanism", ["mechanistic_study"],
                     "Pages 3–4, Scheme 3 and accompanying discussion", "proposed LCuH catalytic cycle",
                     "proposed hydrocupration, imine addition, protonation or carbamate migration",
                     unresolved=["Author-proposed mechanism is not a discrete performed experiment; intermediates 20–24 are not isolated records."])
add_generic_participants(mech, [("15a", "reactant"), ("9a", "reactant"), ("18a", "product"), ("19a", "product")], 3)
setf(mech["classification"], "coverage", "unresolved", "unresolved", [E[3], E[4]], method="inference",
     derivation="Elementary-step atom mapping is not established from the proposed cycle.", verification="inferred")
note(mech, "The t-BuOH branch is proposed to protonate the Cu–N intermediate and give 18a; without t-BuOH, carbonyl migration gives 19a.", [E[3], E[4]], "original_extraction")
reactions.append(mech)


def add_common_coupling(r, imine, products, ev, ligand="pcy3", ligand_molpct=6.5, allene_mass=96.6,
                        allene_amount=0.48, imine_amount=0.4, tbuoh=False, solvent="toluene",
                        temperature="rt", time_h=24, silane_amount=0.8, silane_volume=99,
                        workup=True, internal_standard=False):
    participant(r, "p_imine", imine, "reactant", ev, amount=q(imine_amount, "mmol"))
    participant(r, "p_allene", "15a", "reactant", ev, mass=q(allene_mass, "mg"), amount=q(allene_amount, "mmol"))
    for i, product in enumerate(products, 1):
        participant(r, f"p_product_{i}", product, "product", ev, product_role="main product" if i == 1 else "minor/byproduct")
    participant(r, "p_cu", "cu_oac2", "catalyst", ev, molpct=q(5, "mol%"), ref=imine)
    participant(r, "p_ligand", ligand, "ligand", ev, molpct=q(ligand_molpct, "mol%"), ref=imine)
    participant(r, "p_silane", "silane", "reagent", ev,
                amount=q(silane_amount, "mmol") if silane_amount is not None else None,
                volume=q(silane_volume, "uL") if silane_volume is not None else None,
                equiv=q(2, "equiv") if silane_amount is not None else None, ref=imine)
    if tbuoh:
        participant(r, "p_tbuoh", "tbuoh", "additive", ev, amount=q(0.8, "mmol"), volume=q(76, "uL"), equiv=q(2, "equiv"), ref=imine)
    participant(r, "p_solvent", solvent, "solvent", ev, volume=q(1.0, "mL"))
    if workup:
        participant(r, "p_nh4f", "nh4f", "workup_agent", ev, stage="workup")
        participant(r, "p_meoh", "meoh", "solvent", ev, stage="workup")
    if internal_standard:
        participant(r, "p_standard", "dimethyl_fumarate", "reagent", ev, stage="characterization")
    setf(r, "reaction_scale", f"{imine_amount} mmol {imine}", {"reference_molecule_id": imine, "amount": q(imine_amount, "mmol")}, ev)
    setf(r["solvent_conditions"], "composition", molecules[solvent]["names"]["original"]["normalized_value"], [solvent], ev)
    setf(r["solvent_conditions"], "total_volume", "1.0 mL", q(1.0, "mL"), ev)
    pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    add_step(r, "reaction", ev, pids, "combine and stir", temperature, q(time_h, "h"))
    if workup:
        add_step(r, "workup", ev, ["p_nh4f", "p_meoh"], "treat with NH4F in MeOH", stage="workup")


# Page 3, Table 1 optimization (13 independent rows).
table1 = [
    (1, "dcpe", "toluene", False, "<5", 58, 20),
    (2, "pcy3", "toluene", False, "<5", 86, "<5"),
    (3, "p_adam3", "toluene", False, "<5", 28, 31),
    (4, "xphos", "toluene", False, "<5", 22, 38),
    (5, "p_nme2_3", "toluene", False, "<5", 66, "<5"),
    (6, "p_oet3", "toluene", False, "<5", 60, "<5"),
    (7, "ph_o2_pnme2", "toluene", False, "<5", 66, "<5"),
    (8, "simes", "toluene", False, "<5", 8, 67),
    (9, "pcy3", "mtbe", False, "<5", 51, "<5"),
    (10, "pcy3", "dioxane", False, "<5", 54, "<5"),
    (11, "pcy3", "dcm", False, "<5", 52, 21),
    (12, "pcy3", "thf", False, "<5", 31, 14),
    (13, "pcy3", "toluene", True, 90, "<5", "<5"),
]


def percent_value(v):
    return q(5, "%", "lt", "<5%") if v == "<5" else q(v, "%")


for row, ligand, solvent, with_tbuoh, y18, y19, rec9 in table1:
    ev = [E[3], E[4]]
    r = make_reaction(f"table1_entry_{row:02d}", ev, "specific_experiment", ["optimization"],
                      f"Page 3, Table 1, entry {row}", "ligand/solvent optimization of imine–allenamide reductive coupling",
                      "Cu-catalyzed reductive coupling", [CURRENT_SCHEME],
                      ["Silane dose is not printed in the Table 1 footnote; retained as unreported."])
    add_common_coupling(r, "9a", ["18a", "19a"], ev, ligand=ligand, ligand_molpct=6,
                        tbuoh=with_tbuoh, solvent=solvent, silane_amount=None, silane_volume=None,
                        internal_standard=True)
    add_outcome(r, "yield_18a", "yield", ["18a"], ev, str(y18), percent_value(y18),
                "1H NMR of unpurified reaction mixture using dimethyl fumarate as analytical standard",
                "1H NMR yield on unpurified mixture; dimethyl fumarate standard", "NMR yield", "product")
    add_outcome(r, "yield_19a", "yield", ["19a"], ev, str(y19), percent_value(y19),
                "1H NMR of unpurified reaction mixture using dimethyl fumarate as analytical standard",
                "1H NMR yield on unpurified mixture; dimethyl fumarate standard", "NMR yield", "product")
    add_outcome(r, "recovery_9a", "recovery", ["9a"], ev, str(rec9), percent_value(rec9),
                "1H NMR of unpurified reaction mixture using dimethyl fumarate as analytical standard",
                "unreacted imine percentage by 1H NMR; dimethyl fumarate standard", "recovered starting material")
    main = "18a" if row == 13 else "19a"
    add_outcome(r, "diastereoselectivity", "selectivity", [main], ev, "single diastereomer",
                {"metric": "dr", "reported": "single diastereomer", "major_configuration": "unassigned", "minor_configuration": "unassigned"},
                "by analysis of the unpurified reaction mixture by 1H NMR spectroscopy", "crude 1H NMR spectroscopy")
    classify_partial(r, ev, ["C-C", "C-H"] + (["C-N"] if main == "19a" else []),
                     "Product-selective net changes are partial; rearrangement cleavages and complete atom mapping remain unresolved.")
    reactions.append(r)


# Page 4 Scheme 4 scope: 21 rows.
scheme4 = {
    "18a": (85, None), "18b": (80, None), "18c": (89, None), "18d": (74, None), "18e": (77, None),
    "18f": (81, None), "18g": (66, None), "18h": (53, None), "18i": (50, "91:9 18i:19i"),
    "18j": (82, None), "18k": (78, None), "18l": (73, None), "18m": (74, "95:5 18m:19m"),
    "18n": (81, None), "18o": (75, None), "18p": (74, None), "18q": (67, None), "18r": (76, None),
    "18s": (86, None), "18t": (74, None), "18u": (79, "95:5 18u:19u"),
}
heated = {"18i", "18j", "18m"}
mixed_minor = {"18i": "19i", "18m": "19m", "18u": "19u"}
for product, (yield_no, ratio) in scheme4.items():
    letter = product[2:]
    imine = "9" + letter
    ev = [E[4], E[5]]
    prods = [product] + ([mixed_minor[product]] if product in mixed_minor else [])
    r = make_reaction(f"scheme4_scope_{product}", ev, "specific_experiment", ["substrate_scope"],
                      f"Page 4, Scheme 4, {product}", "Cu-catalyzed reductive coupling to branched diamino synthon",
                      "stereoselective imine–allenamide reductive coupling", [CURRENT_SCHEME],
                      ["Complete structures are available only as complex source drawings; SMILES remain unencoded."])
    add_common_coupling(r, imine, prods, ev, tbuoh=True, temperature=q(65, "°C") if product in heated else "rt")
    yield_targets = prods if ratio else [product]
    add_outcome(r, "isolated_yield", "yield", yield_targets, ev, f"{yield_no}%", q(yield_no, "%"),
                "isolated yield", "isolated product" if not ratio else "isolated inseparable mixture",
                "combined isolated yield" if ratio else "isolated yield", "main product")
    add_outcome(r, "dr", "selectivity", [product], ev, ">98:2 dr",
                {"metric": "dr", "ratio": ">98:2", "order": "reported major:minor", "major_configuration": "(S,S,S)", "minor_configuration": "unassigned"},
                "single diastereomer by analysis of the unpurified reaction mixture by 1H NMR spectroscopy", "crude 1H NMR spectroscopy")
    if ratio:
        add_outcome(r, "product_ratio", "selectivity", prods, ev, ratio,
                    {"metric": "product ratio", "ratio": ratio.split()[0], "order": ratio.split()[1]},
                    "isolated as an inseparable mixture; ratio measurement method not reported", "method unreported")
    classify_partial(r, ev, ["C-C", "C-H"], "Net C–C formation and hydrogen gain are inferred; complete atom mapping and imine bond-order decrease are unresolved.")
    reactions.append(r)

# Explicit qualitative failed aliphatic substrate statement on page 5.
ev = [E[4], E[5]]
r = make_reaction("scheme4_aliphatic_aldimine_no_product", ev, "specific_experiment", ["substrate_scope", "control"],
                  "Page 5, Scheme 4 discussion", "aliphatic aldimine scope limitation", "attempted reductive coupling",
                  [CURRENT_SCHEME], ["No specific product identity, yield number, or analytical detection limit is reported."])
add_common_coupling(r, "9p_aliphatic", [], ev, tbuoh=True)
add_outcome(r, "no_desired_product", "no_reaction", [], ev, "did not provide any desired products", "no desired product",
            "measurement method not reported", "method unreported")
setf(r["classification"], "coverage", "not_applicable", "not_applicable", ev, method="inference",
     derivation="No observed desired transformation; expected bonds are not classified.", verification="inferred")
r["classification"]["limitations"] = ["Negative result; no intended product or bond-formation label is invented."]
reactions.append(r)

# Page 5 Table 2, imine N-substitution study.
table2 = [("9pa", "19pa", 17), ("9pb", "19pb", 58), ("9pc", "19pc", 70), ("9pd", "19pd", 71), ("9pe", "19pe", 79)]
for idx, (imine, product, y) in enumerate(table2, 1):
    ev = [E[5]]
    r = make_reaction(f"table2_entry_{idx}_{imine}", ev, "specific_experiment", ["optimization"],
                      f"Page 5, Table 2, entry {idx}", "effect of imine N-substitution on urea-forming coupling",
                      "Cu-catalyzed reductive coupling with carbamate migration", [CURRENT_SCHEME],
                      ["Table 2 reports 106 mg for 0.480 mmol of 15a, unlike the 96.6 mg used in Schemes 4/Table 3; source value retained."])
    add_common_coupling(r, imine, [product], ev, ligand_molpct=6, allene_mass=106,
                        temperature=q(65, "°C"), internal_standard=True)
    add_outcome(r, "nmr_yield", "yield", [product], ev, f"{y}%", q(y, "%"),
                "1H NMR of unpurified reaction mixture using dimethyl fumarate as analytical standard",
                "1H NMR yield on unpurified mixture; dimethyl fumarate standard", "NMR yield", "product")
    add_outcome(r, "diastereoselectivity", "selectivity", [product], ev, "single diastereomer",
                {"metric": "dr", "reported": "single diastereomer", "major_configuration": "unassigned", "minor_configuration": "unassigned"},
                "by analysis of the unpurified reaction mixture by 1H NMR spectroscopy", "crude 1H NMR spectroscopy")
    classify_partial(r, ev, ["C-C", "C-H", "C-N"], "Net coupling and carbonyl-migration closure are partial inferred changes; complete atom mapping is unresolved.")
    reactions.append(r)

# Page 5 Table 3 urea scope.
table3 = {
    "19a": ("9a", 94, None), "19b": ("9b", 86, None), "19de": ("9de", 71, None),
    "19ee": ("9ee", 74, None), "19j": ("9j", 69, None), "19k": ("9k", 99, "86:14 19k:18k"),
    "19l": ("9l", 94, None), "19m": ("9m", 69, None), "19n": ("9n", 98, None),
    "19s": ("9s", 83, None), "19t": ("9t", 99, "92:8 19t:18t"), "19pe": ("9pe", 79, None),
}
heated3 = {"19de", "19ee", "19j", "19m", "19pe"}
minor3 = {"19k": "18k", "19t": "18t"}
for product, (imine, y, ratio) in table3.items():
    ev = [E[5]]
    prods = [product] + ([minor3[product]] if product in minor3 else [])
    r = make_reaction(f"table3_scope_{product}", ev, "specific_experiment", ["substrate_scope"],
                      f"Page 5, Table 3, {product}", "Cu-catalyzed reductive coupling to chiral urea",
                      "reductive coupling followed by carbamate migration", [CURRENT_SCHEME],
                      ["Complex structures remain unencoded; mixture ratios are retained without assigning the minor stereochemistry."])
    add_common_coupling(r, imine, prods, ev, temperature=q(65, "°C") if product in heated3 else "rt")
    targets = prods if ratio else [product]
    add_outcome(r, "isolated_yield", "yield", targets, ev, f"{y}%", q(y, "%"),
                "isolated yield", "isolated product" if not ratio else "isolated inseparable mixture",
                "combined isolated yield" if ratio else "isolated yield", "main product")
    add_outcome(r, "dr", "selectivity", [product], ev, ">98:2 dr",
                {"metric": "dr", "ratio": ">98:2", "order": "reported major:minor", "major_configuration": "(S,S,S)", "minor_configuration": "unassigned"},
                "single diastereomer by analysis of the unpurified reaction mixture by 1H NMR spectroscopy", "crude 1H NMR spectroscopy")
    if ratio:
        add_outcome(r, "product_ratio", "selectivity", prods, ev, ratio,
                    {"metric": "product ratio", "ratio": ratio.split()[0], "order": ratio.split()[1]},
                    "isolated as an inseparable mixture; ratio measurement method not reported", "method unreported")
    classify_partial(r, ev, ["C-C", "C-H", "C-N"], "Net coupling and urea closure are partial inferred changes; complete atom mapping is unresolved.")
    reactions.append(r)

# Page 6 Scheme 5 and Scheme 6 application experiments. Later experimental text is contextual evidence only.
ev = [E[6]]
r = make_reaction("scheme5_18c_hydrochloride", ev, "specific_experiment", ["derivatization"],
                  "Page 6, Scheme 5", "formation of crystalline 18c hydrochloride", "salt formation",
                  unresolved=["Reaction time, temperature, HCl amount, isolated amount, and yield are unreported."])
participant(r, "p_18c", "18c", "reactant", ev)
participant(r, "p_hcl", "hcl", "reagent", ev)
participant(r, "p_ether", "ether", "solvent", ev)
participant(r, "p_salt", "18c_hcl", "product", ev, product_role="X-ray sample")
add_step(r, "salt_formation", ev, ["p_18c", "p_hcl", "p_ether"], "treat 18c with HCl in Et2O")
add_outcome(r, "xray_assignment", "other", ["18c_hcl"], ev, "crystalline material; single-crystal X-ray analysis", "single-crystal X-ray stereochemical assignment",
            "single-crystal X-ray analysis", "X-ray crystallography")
classify_partial(r, ev, ["N-H"], "Protonation is represented as a net N–H gain; ionic association and protonation site require expert review.")
reactions.append(r)

ev = [E[6], E[17]]
r = make_reaction("scheme5_18a_to_19a", ev, "specific_experiment", ["derivatization"],
                  "Page 6, Scheme 5; contextual experimental text page 17", "conversion of 18a to urea 19a",
                  "base-induced carbamate migration", unresolved=["Complete atom mapping and all bond cleavages are unresolved."])
participant(r, "p_18a", "18a", "reactant", ev, mass=q(100, "mg"), amount=q(0.190, "mmol"))
participant(r, "p_nbuli", "nbuli", "reagent", ev, volume=q(114, "uL"), amount=q(0.285, "mmol"), equiv=q(1.5, "equiv"), ref="18a")
participant(r, "p_thf", "thf", "solvent", ev, volume=q(1.0, "mL"))
participant(r, "p_19a", "19a", "product", ev, product_role="main product")
setf(r, "reaction_scale", "190 μmol of 18a", {"reference_molecule_id": "18a", "amount": q(190, "umol")}, ev)
add_step(r, "addition", ev, ["p_18a", "p_nbuli", "p_thf"], "add n-BuLi to 18a in THF", q(-10, "°C"))
add_step(r, "warm_hold", ev, ["p_18a", "p_nbuli", "p_thf"], "warm to room temperature and stir", "rt", q(1, "h"))
add_outcome(r, "isolated_yield", "yield", ["19a"], ev, "94 mg (94%)", q(94, "%"), "isolated after flash chromatography", "isolated product", "isolated yield", "main product")
add_outcome(r, "product_amount", "product_amount", ["19a"], ev, "94 mg", q(94, "mg"), "isolated after flash chromatography", "isolated product")
classify_partial(r, ev, ["C-N"], "New urea C–N closure is inferred; carbonyl migration and broken bonds are not fully mapped.")
reactions.append(r)

# Gram-scale reaction appears on page 1 and again on page 6; one experiment, merged evidence.
ev = [E[1], E[6], E[17]]
r = make_reaction("scheme6_gram_scale_18c", ev, "specific_experiment", ["scale_up"],
                  "Pages 1 and 6, graphical abstract/Scheme 6; contextual experimental text page 17",
                  "gram-scale Cu-catalyzed reductive coupling to 18c", "stereoselective imine–allenamide reductive coupling",
                  [CURRENT_SCHEME], ["The dr basis is not explicitly stated for the scale-up scheme; it is not silently inherited from the small-scale scope footnote."])
participant(r, "p_9c", "9c", "reactant", ev, mass=q(1.0, "g"), amount=q(4.44, "mmol"))
participant(r, "p_15a", "15a", "reactant", ev, mass=q(1.07, "g"), amount=q(5.33, "mmol"), equiv=q(1.2, "equiv"), ref="9c")
participant(r, "p_cu", "cu_oac2", "catalyst", ev, mass=q(16.1, "mg"), amount=q(88.9, "umol"), molpct=q(2.0, "mol%"), ref="9c")
participant(r, "p_pcy3", "pcy3", "ligand", ev, volume=q(194, "uL"), amount=q(111, "umol"), molpct=q(2.5, "mol%"), ref="9c")
participant(r, "p_tbuoh", "tbuoh", "additive", ev, volume=q(0.85, "mL"), amount=q(8.89, "mmol"), equiv=q(2, "equiv"), ref="9c")
participant(r, "p_silane", "silane", "reagent", ev, volume=q(1.1, "mL"), amount=q(8.89, "mmol"), equiv=q(2, "equiv"), ref="9c")
participant(r, "p_toluene", "toluene", "solvent", ev, volume=q(10, "mL"))
participant(r, "p_18c", "18c", "product", ev, product_role="main product")
participant(r, "p_nh4f", "nh4f", "workup_agent", ev, stage="workup", mass=q(2, "g"))
participant(r, "p_meoh", "meoh", "solvent", ev, stage="workup", volume=q(20, "mL"))
setf(r, "reaction_scale", "1.0 g (4.44 mmol) imine 9c", {"reference_molecule_id": "9c", "amount": q(4.44, "mmol")}, ev)
setf(r["solvent_conditions"], "composition", "toluene (5 mL catalyst solution + 5 mL substrate solution)", ["toluene"], ev)
setf(r["solvent_conditions"], "total_volume", "10 mL", q(10, "mL"), ev)
add_step(r, "catalyst_prep", ev, ["p_cu", "p_pcy3", "p_toluene"], "evacuate/backfill three times; stir Cu(OAc)2 with 20 wt% PCy3 solution", "rt", q(10, "min"), atmosphere="nitrogen")
add_step(r, "reaction", ev, ["p_9c", "p_15a", "p_cu", "p_pcy3", "p_tbuoh", "p_silane", "p_toluene"], "combine catalyst and substrates, add silane, and stir", "rt", q(2, "h"), atmosphere="nitrogen")
add_step(r, "workup", ev, ["p_nh4f", "p_meoh"], "transfer into NH4F/MeOH and stir", "rt", q(30, "min"), stage="workup")
add_outcome(r, "isolated_yield", "yield", ["18c"], ev, "1.655 g (91%)", q(91, "%"), "isolated after silica gel chromatography", "isolated product", "isolated yield", "main product")
add_outcome(r, "product_amount", "product_amount", ["18c"], ev, "1.655 g", q(1.655, "g"), "isolated after silica gel chromatography", "isolated product")
add_outcome(r, "dr", "selectivity", ["18c"], [E[1], E[6]], ">98:2 dr; single isomer", {"metric": "dr", "ratio": ">98:2", "order": "reported major:minor", "major_configuration": "(S,S,S)", "minor_configuration": "unassigned"},
            "measurement method not explicitly stated for the gram-scale scheme", "method unreported")
classify_partial(r, ev, ["C-C", "C-H"], "Net C–C formation and hydrogen gain are inferred; complete atom mapping remains unresolved.")
reactions.append(r)

# Three-step overall deprotection, with stage temperatures/times paired.
ev = [E[6], E[17]]
r = make_reaction("scheme6_19a_to_25", ev, "specific_experiment", ["derivatization"],
                  "Page 6, Scheme 6; contextual experimental text page 17", "three-step conversion of 19a to urea 25",
                  "alcohol activation, elimination, and acid hydrolysis",
                  unresolved=["The 70% yield is overall from 19a to 25 and is not assigned to the final hydrolysis alone."])
participant(r, "p_19a", "19a", "reactant", ev, mass=q(207.5, "mg"), amount=q(0.393, "mmol"))
participant(r, "p_et3n", "et3n", "reagent", ev, volume=q(66, "uL"), amount=q(0.473, "mmol"))
participant(r, "p_mscl", "mscl", "reagent", ev, volume=q(30.5, "uL"), amount=q(0.393, "mmol"))
participant(r, "p_kotbu", "kotbu", "reagent", ev, volume=q(433, "uL"), amount=q(0.433, "mmol"))
participant(r, "p_h2so4", "h2so4", "reagent", ev, volume=q(788, "uL"), amount=q(3.94, "mmol"))
participant(r, "p_dcm", "dcm", "solvent", ev, volume=q(2, "mL"))
participant(r, "p_thf", "thf", "solvent", ev)
participant(r, "p_25", "25", "product", ev, product_role="main product")
setf(r, "reaction_scale", "393 μmol 19a", {"reference_molecule_id": "19a", "amount": q(393, "umol")}, ev)
add_step(r, "mesylation_addition", ev, ["p_19a", "p_et3n", "p_mscl", "p_dcm"], "add triethylamine then MsCl dropwise", q(0, "°C"))
add_step(r, "mesylation_hold", ev, ["p_19a", "p_et3n", "p_mscl", "p_dcm"], "stir", q(0, "°C"), q(30, "min"))
add_step(r, "elimination_addition", ev, ["p_kotbu", "p_thf"], "cool crude intermediate in THF and add KOtBu", q(0, "°C"))
add_step(r, "elimination_hold", ev, ["p_kotbu", "p_thf"], "warm to room temperature and stir", "rt", q(30, "min"))
add_step(r, "hydrolysis", ev, ["p_h2so4", "p_thf"], "add aqueous H2SO4, purge with argon, seal and heat", q(50, "°C"), q(3, "h"), atmosphere="argon")
add_outcome(r, "overall_yield", "yield", ["25"], ev, "112 mg (70%) from 19a", q(70, "%"), "isolated after the complete three-step sequence", "overall isolated product from 19a to 25", "three-step overall isolated yield", "main product")
add_outcome(r, "product_amount", "product_amount", ["25"], ev, "112 mg", q(112, "mg"), "isolated after flash chromatography", "isolated product")
classify_partial(r, ev, ["N-H"], "Net N–H formation in the deprotected urea is inferred; eliminated/hydrolyzed fragments and all broken bonds are unresolved.")
reactions.append(r)

ev = [E[6], E[17], E[18]]
r = make_reaction("scheme6_18c_to_26", ev, "specific_experiment", ["derivatization"],
                  "Page 6, Scheme 6; contextual experimental text pages 17–18", "N-allylation of 18c to 26", "N-allylation",
                  unresolved=["Complete atom mapping and C–Br cleavage are not encoded."])
participant(r, "p_18c", "18c", "reactant", ev, mass=q(100, "mg"), amount=q(0.233, "mmol"))
participant(r, "p_allyl_br", "allyl_br", "reagent", ev, volume=q(101, "uL"), amount=q(1.17, "mmol"), equiv=q(5, "equiv"), ref="18c")
participant(r, "p_k2co3", "k2co3", "reagent", ev, mass=q(161, "mg"), amount=q(1.17, "mmol"), equiv=q(5, "equiv"), ref="18c")
participant(r, "p_tbai", "tbai", "reagent", ev, mass=q(17.2, "mg"), amount=q(46.7, "umol"), molpct=q(20, "mol%"), ref="18c")
participant(r, "p_mecn", "mecn", "solvent", ev, volume=q(1, "mL"))
participant(r, "p_26", "26", "product", ev, product_role="main product")
setf(r, "reaction_scale", "233 μmol 18c", {"reference_molecule_id": "18c", "amount": q(233, "umol")}, ev)
add_step(r, "reaction", ev, ["p_18c", "p_allyl_br", "p_k2co3", "p_tbai", "p_mecn"], "heat reaction mixture", q(85, "°C"), q(18, "h"))
add_outcome(r, "isolated_yield", "yield", ["26"], ev, "90 mg (82%)", q(82, "%"), "isolated after flash chromatography", "isolated product", "isolated yield", "main product")
add_outcome(r, "product_amount", "product_amount", ["26"], ev, "90 mg", q(90, "mg"), "isolated after flash chromatography", "isolated product")
classify_partial(r, ev, ["C-N"], "N-allyl C–N formation is inferred; leaving-group cleavage and complete atom mapping are unresolved.")
reactions.append(r)

ev = [E[6], E[18]]
r = make_reaction("scheme6_26_to_27", ev, "specific_experiment", ["derivatization"],
                  "Page 6, Scheme 6; contextual experimental text page 18", "ring-closing metathesis of 26 to 27", "ring-closing olefin metathesis",
                  unresolved=["Complete alkene atom mapping and the volatile metathesis coproduct are not encoded."])
participant(r, "p_26", "26", "reactant", ev, mass=q(48, "mg"), amount=q(0.10, "mmol"))
participant(r, "p_hg2", "hg2", "catalyst", ev, mass=q(3.2, "mg"), amount=q(5.1, "umol"), molpct=q(5.1, "mol%"), ref="26")
participant(r, "p_toluene", "toluene", "solvent", ev, volume=q(2, "mL"))
participant(r, "p_27", "27", "product", ev, product_role="main product")
setf(r, "reaction_scale", "0.10 mmol 26", {"reference_molecule_id": "26", "amount": q(0.10, "mmol")}, ev)
add_step(r, "reaction", ev, ["p_26", "p_hg2", "p_toluene"], "seal and heat", q(90, "°C"), q(12, "h"), atmosphere="argon")
add_outcome(r, "isolated_yield", "yield", ["27"], ev, "35 mg (78%)", q(78, "%"), "isolated after flash chromatography", "isolated product", "isolated yield", "main product")
add_outcome(r, "product_amount", "product_amount", ["27"], ev, "35 mg", q(35, "mg"), "isolated after flash chromatography", "isolated product")
classify_partial(r, ev, ["DB-C-C"], "New ring-closing C=C connection is inferred; full metathesis atom mapping is unresolved.")
reactions.append(r)

# Source context preserves non-synthetic figures and exact page inventory decisions.
bundle["source_context"] = [
    {"document_id": DOC, "evidence_ids": [E[1]], "kind": "article_metadata",
     "text": TITLE + "; abstract and graphical abstract. The graphical scale-up is linked to the single Scheme 6 scale-up record."},
    {"document_id": DOC, "evidence_ids": [E[2]], "kind": "nonreaction_figure",
     "text": "Figure 1 shows selected biologically active diamine/urea compounds and is not a synthesis experiment."},
    {"document_id": DOC, "evidence_ids": [E[2]], "kind": "prior_art_inventory",
     "text": "Scheme 1 conceptual prior-art arrows are retained as literature_example records; they are not experiments of this paper."},
    {"document_id": DOC, "evidence_ids": [E[3], E[4]], "kind": "mechanistic_discussion",
     "text": "Scheme 3 is an author-proposed catalytic cycle. It is not treated as isolation of intermediates 20–24."},
    {"document_id": DOC, "evidence_ids": [E[5], E[6]], "kind": "stereochemical_context",
     "text": "Products 18/19 are reported as the (S,S,S) diastereomer; 18c·HCl was assigned by X-ray and 18a→19a established correspondence of stereoisomers."},
    {"document_id": DOC, "evidence_ids": [E[6]], "kind": "proposed_application",
     "text": "The drawn route from 27 toward aminopiperidine 28 and NK-1 antagonists is a prospective synthetic application, not a performed experiment on this page."},
]

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
write(bundle_path, bundle)
write(ROOT / "runs/current/raw_model/stereo_diagnostics_pages_0001_0006.json", stereo_diagnostics)
write(ROOT / "runs/current/raw_model/pages_0001_0006_inventory.json", {
    "document_id": DOC,
    "pages": [1, 2, 3, 4, 5, 6],
    "reaction_count": len(reactions),
    "molecule_count": len(molecules),
    "record_ids": [r["reaction_id"] for r in reactions],
    "notes": [
        "Original images for pages 1–6 were inspected individually.",
        "Later page text is used only as contextual procedural evidence for Scheme 5/6 experiments; no later-page semantic coverage is claimed.",
        "No SI was associated or inferred.",
    ],
})

ids = [r["reaction_id"] for r in reactions]
def selected(*prefixes):
    return [rid for rid in ids if any(rid.split(":", 1)[1].startswith(p) for p in prefixes)]

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 1,
            "state": "attempted",
            "regions": [
                {"region": "top title, citation, and graphical abstract", "finding": "Article metadata and the 1.0 g/91%/single-isomer scale-up graphic; linked to the same scale-up later shown in Scheme 6.", "evidence_ids": [E[1], f"{DOC}:page1:image"]},
                {"region": "abstract and introduction columns", "finding": "Motivation/background only; no additional row-level experiment.", "evidence_ids": [E[1], f"{DOC}:page1:image"]},
            ],
            "record_ids": selected("scheme6_gram_scale_18c"),
            "unresolved": ["The graphical abstract does not state the analytical method for 'single isomer'."],
        },
        {
            "page": 2,
            "state": "attempted",
            "regions": [
                {"region": "Figure 1", "finding": "Selected biologically active structures; nonreaction context, not synthesis experiments.", "evidence_ids": [E[2], f"{DOC}:page2:image"]},
                {"region": "Scheme 1A and 1B", "finding": "Eight conceptual/prior-art transformation arrows inventoried as literature_example records, not current-paper experiments.", "evidence_ids": [E[2], f"{DOC}:page2:image"]},
                {"region": "right-column discussion", "finding": "Prior asymmetric imine allylation literature context; no fully specified additional experiment.", "evidence_ids": [E[2], f"{DOC}:page2:image"]},
            ],
            "record_ids": selected("scheme1"),
            "unresolved": ["Scheme 1 uses generic R groups and supplies no row-level conditions or yields; structures remain generic/unencoded."],
        },
        {
            "page": 3,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 2A–D", "finding": "Classical/prior/proposed strategies separated from the present-work generic coupling scheme.", "evidence_ids": [E[3], f"{DOC}:page3:image"]},
                {"region": "Table 1 entries 1–13 and footnotes", "finding": "Thirteen independent ligand/solvent/additive optimization experiments extracted with product yields and unreacted-imine percentages.", "evidence_ids": [E[3], f"{DOC}:page3:image"]},
                {"region": "Scheme 3 (lower page)", "finding": "Author-proposed catalytic cycle retained as proposed_mechanism, not isolated intermediate experiments.", "evidence_ids": [E[3], f"{DOC}:page3:image"]},
            ],
            "record_ids": selected("scheme2", "table1", "scheme3"),
            "unresolved": ["Table 1 prints the silane identity but not a silane dose in its footnote; no dose was inferred.", "Structures 9a/15a/18a/19a are complex drawings and remain without source-verified SMILES."],
        },
        {
            "page": 4,
            "state": "attempted",
            "regions": [
                {"region": "upper mechanistic discussion", "finding": "Continuation of Scheme 3 rationale and Table 1 interpretation; evidence merged into existing records.", "evidence_ids": [E[4], f"{DOC}:page4:image"]},
                {"region": "Scheme 4 substrate scope and footnotes", "finding": "Twenty-one product-specific scope experiments (18a–18u) extracted separately, including heated rows and inseparable product mixtures.", "evidence_ids": [E[4], f"{DOC}:page4:image"]},
            ],
            "record_ids": selected("table1", "scheme3", "scheme4_scope_"),
            "unresolved": ["Full product graphs and minor-product ratio measurement methods are not supplied in machine-readable form; SMILES/methods remain unresolved rather than guessed."],
        },
        {
            "page": 5,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 4 discussion", "finding": "Scope interpretation and explicit failed aliphatic aldimine statement; negative experiment retained without a fabricated product.", "evidence_ids": [E[5], f"{DOC}:page5:image"]},
                {"region": "Table 2 entries 1–5 and footnote", "finding": "Five N-substituent optimization rows extracted; table-specific 106 mg 15a charge preserved.", "evidence_ids": [E[5], f"{DOC}:page5:image"]},
                {"region": "Table 3 scope and footnotes", "finding": "Twelve urea-forming scope rows extracted with isolated yields, >98:2 dr, heated cases, and inseparable product ratios.", "evidence_ids": [E[5], f"{DOC}:page5:image"]},
            ],
            "record_ids": selected("scheme4_scope_", "scheme4_aliphatic", "table2", "table3"),
            "unresolved": ["The failed aliphatic substrate statement reports no analytical method or detection limit.", "Table 2's 106 mg/0.480 mmol charge differs from the 96.6 mg/0.480 mmol charge in Scheme 4/Table 3; both source values are retained without correction."],
        },
        {
            "page": 6,
            "state": "attempted",
            "regions": [
                {"region": "Scheme 5", "finding": "18c hydrochloride/X-ray assignment and 18a→19a (94%) derivatization extracted as distinct experiments.", "evidence_ids": [E[6], f"{DOC}:page6:image"]},
                {"region": "Scheme 6", "finding": "Three-step 19a→25 overall conversion, duplicate-linked gram-scale coupling, N-allylation to 26, and metathesis to 27 extracted separately.", "evidence_ids": [E[6], f"{DOC}:page6:image"]},
                {"region": "Scheme 6 prospective right-hand route", "finding": "27→28/CP-99,994/CP-122,721 is presented as synthetic application/route context, not a performed experiment here.", "evidence_ids": [E[6], f"{DOC}:page6:image"]},
            ],
            "record_ids": selected("scheme5", "scheme6"),
            "unresolved": ["18c·HCl formation lacks quantities, time, temperature, and yield.", "The gram-scale dr method is not explicitly stated in the Scheme 6/graphical-abstract evidence.", "Later-page procedural text was used only as contextual evidence; no semantic coverage beyond page 6 is claimed."],
        },
    ],
}
write(ROOT / "runs/current/coverage/pages_0001_0006.json", coverage)
print(json.dumps({"reactions": len(reactions), "molecules": len(molecules)}, ensure_ascii=False))
