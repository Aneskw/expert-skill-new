#!/usr/bin/env python3
"""Source transcription for doc_a43f9491a3cfbff6, PDF pages 1-5.

This is the raw machine-authored construction script.  It intentionally keeps
dense source drawings unresolved and every imported record at needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st
import reaction_tools as rt

DOC = "doc_a43f9491a3cfbff6"
OUT = ROOT / "runs/current"
P1 = f"{DOC}:page1:image"
P2 = f"{DOC}:page2:image"
P3 = f"{DOC}:page3:image"
P4 = f"{DOC}:page4:image"
P5 = f"{DOC}:page5:image"
T1 = f"{DOC}:page1:text"
T2 = f"{DOC}:page2:text"
T3 = f"{DOC}:page3:text"
T4 = f"{DOC}:page4:text"
T5 = f"{DOC}:page5:text"


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status)


def derived(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def inferred(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="inference", derivation=derivation)


def set_occurrence(m, label, evidence):
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, evidence),
        "evidence_ids": list(evidence),
    }]


def set_names(m, original, english, evidence, abbreviations=None, aliases=None):
    m["names"]["original"] = f(original, original, evidence)
    m["names"]["english"] = f(english, english, evidence)
    if abbreviations:
        m["names"]["abbreviations"] = f(abbreviations, abbreviations, evidence)
    if aliases:
        m["names"]["aliases"] = f(aliases, aliases, evidence)


def set_stereo_unresolved(m, evidence, note, chirality="unresolved", geometry="unresolved",
                          representation="unresolved"):
    m["stereochemistry"]["atom_centers"] = f(
        "Source drawing reviewed; exact atom-center assignment is not fully encoded.",
        {"status": "reviewed_from_drawing", "detail": note}, evidence)
    m["stereochemistry"]["double_bonds"] = f(
        "Source drawing reviewed; E/Z status is retained only where explicitly stated.",
        {"status": geometry, "detail": note}, evidence)
    m["stereochemistry"]["axial"] = f(
        "Source drawing reviewed for allene/axial stereochemistry.",
        {"status": chirality, "detail": note}, evidence)
    m["stereochemistry"]["mixture"] = f(evidence=evidence)
    m["stereochemistry"]["assignment"] = f(note, note, evidence)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence, "Assigned conservatively from the source drawing and explicit stereochemical text.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "Assigned conservatively from the source drawing and explicit stereochemical text.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence, "No complete stereochemical graph was asserted for the dense source drawing.")


def molecule(mid, label, name, evidence, smiles=None, note=None,
             chirality=None, geometry=None, representation=None):
    evidence = list(evidence)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    set_occurrence(m, label, evidence)
    set_names(m, label, name, evidence)
    m["cas_numbers"] = f(evidence=evidence)
    m["identity_note"] = f(note, note, evidence) if note else f(evidence=evidence)
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence,
                              "Candidate connectivity transcribed from the unambiguous source label/drawing; parsing does not verify identity.")
        m["stereochemistry"]["atom_centers"] = derived(diag["atom_centers"], evidence,
                                                        "RDKit graph diagnostic on the candidate SMILES.")
        m["stereochemistry"]["double_bonds"] = derived(diag["double_bonds"], evidence,
                                                        "RDKit graph diagnostic; allene cumulated bonds do not by themselves imply E/Z.")
        m["stereochemistry"]["axial"] = derived([], evidence,
                                                 "Source structure was checked for an allene/biaryl axis; none is stereogenic for this candidate.")
        m["stereochemistry"]["mixture"] = f(evidence=evidence)
        m["stereochemistry"]["assignment"] = derived(
            "No configurational assignment is applicable to this candidate structure.", evidence,
            "No potential stereocenter or stereogenic terminal-allene axis is present in the candidate graph/source structure.")
        m["stereochemistry"]["chirality_status"] = inferred(
            chirality or "achiral", evidence,
            "Candidate graph has no potential atom stereocenter; terminal allenes have a CH2 terminus and are not axially stereogenic.")
        m["stereochemistry"]["geometry_status"] = inferred(
            geometry or "not_applicable", evidence,
            "No source-assigned E/Z element is applicable to this candidate graph.")
        m["stereochemistry"]["representation_status"] = inferred(
            representation or "not_applicable", evidence,
            "No stereochemical element requiring encoding was found for this candidate.")
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence, "Canonical isomeric SMILES from the candidate graph.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], evidence, "Stereo graph key from the candidate graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f(evidence=evidence)
        set_stereo_unresolved(m, evidence, note or "Structure retained by source label/drawing; no complete SMILES asserted.",
                              chirality or "unresolved", geometry or "unresolved",
                              representation or "not_encoded")
    return m


def generic_molecule(mid, label, name, evidence, fixed, variables, stereo_note):
    m = molecule(mid, label, name, evidence, note=stereo_note,
                 chirality="unspecified", geometry="unspecified", representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": f(f"Source-drawn attachment point for {v}", f"source attachment for {v}", evidence),
            "definition": f(desc, desc, evidence),
            "allowed_values": f(evidence=evidence),
            "constraints": f(evidence=evidence),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": f("Generic R/Ar labels are not one concrete molecule.",
                    "Preserve source variables without inventing a concrete SMILES.", evidence),
    }
    return m


def charge(p, key, raw, normalized, evidence, reference=None, method="original_extraction"):
    p["charges"][key] = f(raw, normalized, evidence, method=method)
    if key in ("equivalents", "mol_percent"):
        p["equivalence_reference"] = f(reference, reference, evidence, method=method)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    return p


def base_reaction(suffix, evidence, record_type="specific_experiment", categories=None,
                  location=None, name=None, transformation=None):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = (f(categories, categories, evidence) if categories is not None
                                      else f(evidence=evidence))
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else f(evidence=evidence)
    r["transformation_type"] = (inferred(transformation, evidence,
        "Transformation type is a machine normalization of the source-drawn reaction.") if transformation else f(evidence=evidence))
    r["general_scheme_ids"] = []
    r["participants"] = []
    r["steps"] = []
    r["outcomes"] = []
    r["notes"] = []
    r["evidence_ids"] = list(evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def step(step_id, evidence, participant_ids, operation=None, temperature=None, time=None,
         atmosphere=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(participant_ids)
    if operation is not None:
        s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


def outcome(oid, kind, targets, evidence, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence)
    if product_role is not None:
        o["product_role"] = f(product_role, product_role, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


def classification(evidence, formed=(), decreased=(), broken=(), coverage="partial", limitations=()):
    changes = []
    labels = []
    label_map = {}
    for idx, (a, b) in enumerate(formed, 1):
        cid = f"b{idx}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                            {"atom_map": None, "element": b}],
              "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(derived(ch, evidence, "Category-level bond change read from the source drawing; atom mapping is unresolved."))
        lab = st.bond_label(a, b, 1)
        if lab["code"] not in label_map:
            label_map[lab["code"]] = dict(lab, supporting_change_ids=[])
            labels.append(label_map[lab["code"]])
        label_map[lab["code"]]["supporting_change_ids"].append(cid)
    for a, b, before, after in decreased:
        cid = f"b{len(changes)+1}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                            {"atom_map": None, "element": b}],
              "before_order": before, "after_order": after, "change_type": "order_decreased"}
        changes.append(derived(ch, evidence, "Source-drawn bond-order decrease; atom mapping is unresolved."))
    for a, b in broken:
        cid = f"b{len(changes)+1}"
        ch = {"change_id": cid, "atoms": [{"atom_map": None, "element": a},
                                            {"atom_map": None, "element": b}],
              "before_order": 1, "after_order": 0, "change_type": "broken"}
        changes.append(derived(ch, evidence, "Source-drawn bond cleavage; atom mapping is unresolved."))
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(evidence=evidence),
        "bond_changes": changes,
        "labels": derived(labels, evidence, "Labels are the exact taxonomy projection of the category-level changes."),
        "coverage": inferred(coverage, evidence, "Dense structures and atom correspondence are incomplete; classification is not atom-complete."),
        "limitations": list(limitations) or [
            "No atom-mapped reaction SMILES was asserted for the dense source drawings.",
            "Hydrogen provenance and omitted coproducts are unresolved.",
        ],
        "model_votes": [],
    }


def no_classification(evidence, reason):
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(evidence=evidence),
        "bond_changes": [],
        "labels": inferred([], evidence, reason),
        "coverage": inferred("not_applicable", evidence, reason),
        "limitations": [reason],
        "model_votes": [],
    }


documents = json.loads((OUT / "dataset.json").read_text())["documents"]
bundle = st.empty(documents)
molecules = {}


def add_molecule(m):
    molecules[m["molecule_id"]] = m
    return m["molecule_id"]


# Generic source families and prior-work drawings.
add_molecule(generic_molecule("m_generic_allene", "R-substituted allene", "R-substituted allene",
                              [P1, P2], "R-CH=C=CH2 as drawn", [("R", "organic substituent")],
                              "Terminal CH2 removes axial stereogenicity for the drawn family; R remains variable."))
add_molecule(generic_molecule("m_generic_azirine", "2H-azirine 1", "2H-azirine family 1",
                              [P1, P2, P3], "three-membered 2H-azirine ring bearing Ar1/Ar2 or R1/R2",
                              [("Ar1/R1", "substituent on one azirine carbon"),
                               ("Ar2/R2", "substituent on the other azirine carbon")],
                              "Substituent-dependent center/relative stereochemistry is not specified for the generic family."))
add_molecule(generic_molecule("m_generic_borylated_product", "borylated allylaziridine 4",
                              "borylated allylaziridine product family 4", [P1, P2, P3],
                              "source-drawn allylaziridine bearing Bpin and R/Ar substituents",
                              [("R", "allene-derived substituent"), ("Ar1/R1", "azirine substituent"),
                               ("Ar2/R2", "azirine substituent")],
                              "Wedge/dash bonds specify relative stereochemistry in the source family; absolute configuration is not assigned."))
add_molecule(generic_molecule("m_generic_chlorinated_product", "chlorinated allylaziridine 5",
                              "chlorinated allylaziridine product family 5", [P1, P3],
                              "source-drawn allylaziridine bearing vinyl chloride and R/Ar substituents",
                              [("R", "allene-derived substituent"), ("Ar1/R1", "azirine substituent"),
                               ("Ar2/R2", "azirine substituent")],
                              "Wedge/dash bonds specify relative stereochemistry; absolute configuration is not assigned."))
for mid, label, name in [
    ("m_lit_hoveyda_allene", "R3-substituted allene", "R3-substituted allene (Hoveyda example)"),
    ("m_lit_hoveyda_ketimine", "unprotected ketimine hydrochloride", "unprotected ketimine hydrochloride"),
    ("m_lit_hoveyda_product", "Hoveyda borylated homoallylic amine", "borylated homoallylic amine"),
    ("m_lit_procter_allene", "R4-substituted allene", "R4-substituted allene (Procter example)"),
    ("m_lit_procter_ketiminoester", "ketiminoester", "ketiminoester"),
    ("m_lit_procter_product", "Procter borylated amino ester", "borylated amino ester"),
]:
    add_molecule(molecule(mid, label, name, [P2], note="Generic cited-work drawing with R groups; no concrete SMILES asserted.",
                          chirality="unspecified", geometry="unspecified", representation="not_encoded"))

# Common reagents, ligands, bases, and solvents.  Dense/ionic species retain labels without guessed coordination SMILES.
simple = [
    ("m_b2pin2", "B2(pin)2 (3)", "bis(pinacolato)diboron", "CC1(C)OB(B2OC(C)(C)C(C)(C)O2)OC1(C)C"),
    ("m_kotbu", "tBuOK", "potassium tert-butoxide", "CC(C)(C)[O-].[K+]"),
    ("m_naotbu", "tBuONa", "sodium tert-butoxide", "CC(C)(C)[O-].[Na+]"),
    ("m_liotbu", "tBuOLi", "lithium tert-butoxide", "CC(C)(C)[O-].[Li+]"),
    ("m_kome", "MeOK", "potassium methoxide", "C[O-].[K+]"),
    ("m_phcf3", "PhCF3", "benzotrifluoride", "FC(F)(F)c1ccccc1"),
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("m_thf", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("m_dioxane", "1,4-dioxane", "1,4-dioxane", "C1COCCO1"),
    ("m_meoh", "CH3OH", "methanol", "CO"),
    ("m_water", "H2O", "water", "O"),
    ("m_dcm", "DCM", "dichloromethane", "ClCCl"),
    ("m_benzoyl_chloride", "PhCOCl", "benzoyl chloride", "O=C(Cl)c1ccccc1"),
    ("m_triethylamine", "Et3N", "triethylamine", "CCN(CC)CC"),
    ("m_dmap", "DMAP", "4-(dimethylamino)pyridine", "CN(C)c1ccncc1"),
    ("m_iodobenzene", "PhI", "iodobenzene", "Ic1ccccc1"),
    ("m_naoh", "NaOH", "sodium hydroxide", "[Na+].[OH-]"),
    ("m_bipy", "Bipy", "2,2'-bipyridine", "n1ccccc1-c2ccccn2"),
]
for mid, label, name, smi in simple:
    add_molecule(molecule(mid, label, name, [P2, P3], smiles=smi,
                          note="Candidate SMILES derives from an unambiguous source reagent/solvent label."))
for mid, label, name in [
    ("m_cucl", "CuCl", "copper(I) chloride"), ("m_cubr", "CuBr", "copper(I) bromide"),
    ("m_cui", "CuI", "copper(I) iodide"), ("m_cu_oa2", "Cu(OAc)2", "copper(II) acetate"),
    ("m_cu_otf2", "Cu(OTf)2", "copper(II) triflate"), ("m_cubr2", "CuBr2", "copper(II) bromide"),
    ("m_cucl2", "CuCl2", "copper(II) chloride"), ("m_iprnhc", "IPrNHC", "IPr N-heterocyclic carbene ligand"),
    ("m_iprnhc_hcl", "IPrNHC·HCl", "IPr NHC hydrochloride precursor"),
    ("m_dppf", "DPPF", "1,1'-bis(diphenylphosphino)ferrocene"),
    ("m_dppbz", "DPPBz", "1,2-bis(diphenylphosphino)benzene"),
    ("m_pd_oa2", "Pd(OAc)2", "palladium(II) acetate"), ("m_sphos", "SPhos", "SPhos ligand"),
    ("m_nabo3_4h2o", "NaBO3·4H2O", "sodium perborate tetrahydrate"),
]:
    add_molecule(molecule(mid, label, name, [P2, P3], note="Source-identified reagent; coordination/ionic structure not encoded.",
                          chirality="unresolved", geometry="not_applicable", representation="not_encoded"))

# Concrete allenes 2a-2m have unambiguous terminal-allene drawings; 2n-2p are text-only labels.
allene_rows = [
    ("a", "C6H5", "C=C=Cc1ccccc1"), ("b", "4-MeC6H4", "C=C=Cc1ccc(C)cc1"),
    ("c", "4-EtC6H4", "C=C=Cc1ccc(CC)cc1"), ("d", "4-tBuC6H4", "C=C=Cc1ccc(C(C)(C)C)cc1"),
    ("e", "4-ClC6H4", "C=C=Cc1ccc(Cl)cc1"), ("f", "4-BrC6H4", "C=C=Cc1ccc(Br)cc1"),
    ("g", "4-MeOC6H4", "C=C=Cc1ccc(OC)cc1"), ("h", "4-CF3C6H4", "C=C=Cc1ccc(C(F)(F)F)cc1"),
    ("i", "3-MeC6H4", "C=C=Cc1cccc(C)c1"), ("j", "3-ClC6H4", "C=C=Cc1cccc(Cl)c1"),
    ("k", "2-ClC6H4", "C=C=Cc1ccccc1Cl"), ("l", "3,4-Me2C6H3", "C=C=Cc1ccc(C)c(C)c1"),
    ("m", "cyclohexyl", "C=C=CC1CCCCC1"),
]
for code, r3, smi in allene_rows:
    add_molecule(molecule(f"m_2{code}", f"2{code}", f"terminal allene 2{code} (R3={r3})", [P2, P3],
                          smiles=smi, note=f"Table 2 source drawing/row identifies R3={r3}; candidate connectivity preserves a CH2 terminus."))
for code, name in [("n", "alpha-substituted phenylallene"), ("o", "phenoxyallene"), ("p", "allenamide")]:
    add_molecule(molecule(f"m_2{code}", f"2{code}", name, [P3, T3],
                          note="Mentioned only by class/name in supplied main text; exact substituents/connectivity are not shown on page 3.",
                          chirality="unresolved", geometry="unresolved", representation="not_encoded"))

# Azirines 1a-1l.  The ring drawings are retained without guessed stereochemical SMILES.
azirines = {
    "a": ("C6H5", "C6H5"), "b": ("4-MeC6H4", "C6H5"), "c": ("4-nPrC6H4", "C6H5"),
    "d": ("4-ClC6H4", "C6H5"), "e": ("C6H5", "4-BrC6H4"), "f": ("C6H5", "4-MeC6H4"),
    "g": ("C6H5", "4-ClC6H4"), "h": ("C6H5", "2-ClC6H4"), "i": ("C6H5", "2-MeC6H4"),
    "j": ("C6H5", "3-MeC6H4"), "k": ("C6H5", "3-ClC6H4"), "l": ("C6H5", "H"),
}
for code, (r1, r2) in azirines.items():
    ev = [P2, P3] if code == "a" else [P3]
    add_molecule(molecule(f"m_1{code}", f"1{code}", f"2H-azirine 1{code} (R1={r1}, R2={r2})", ev,
                          note=f"Source drawing/row gives R1={r1}, R2={r2}; configuration is not explicitly assigned.",
                          chirality="unspecified", geometry="not_applicable", representation="not_encoded"))

# Products.  Product drawings carry relative wedge/dash information but not complete absolute assignments.
def product_molecule(mid, label, name, evidence, alkene=False):
    geometry = "unspecified" if alkene else "not_applicable"
    return molecule(mid, label, name, evidence,
                    note="Dense source product drawing has wedge/dash relative stereochemistry; no absolute R/S assignment is printed in the supplied main text.",
                    chirality="partially_specified", geometry=geometry, representation="not_encoded")


add_molecule(product_molecule("m_4aa", "4aa", "borylated allylaziridine 4aa", [P1, P2, P3], alkene=False))
add_molecule(product_molecule("m_5aa", "5aa", "chlorinated allylaziridine 5aa", [P1, P3], alkene=False))
t2 = [
    (1, "a", 83, 78), (2, "b", 67, 60), (3, "c", 65, 62), (4, "d", 72, 65),
    (5, "e", 75, 72), (6, "f", 80, 75), (7, "g", 83, 71), (8, "h", 63, 55),
    (9, "i", 73, 68), (10, "j", 80, 71), (11, "k", 86, 75), (12, "l", 76, 69),
    (13, "m", 36, 24),
]
for entry, code, y4, y5 in t2:
    if code != "a":
        add_molecule(product_molecule(f"m_t2e{entry}_4", f"4 (Table 2 entry {entry})",
                                      f"borylated allylaziridine intermediate, Table 2 entry {entry}", [P3]))
        add_molecule(product_molecule(f"m_5a{code}", f"5a{code}",
                                      f"chlorinated allylaziridine 5a{code}", [P3]))
t3 = [
    (1, "b", 80, 75), (2, "c", 83, 76), (3, "d", 81, 75), (4, "e", 82, 79),
    (5, "f", 80, 72), (6, "g", 78, 72), (7, "h", 79, 75), (8, "i", 83, 75),
    (9, "j", 86, 78), (10, "k", 85, 78),
]
for entry, code, y4, y5 in t3:
    add_molecule(product_molecule(f"m_t3e{entry}_4", f"4 (Table 3 entry {entry})",
                                  f"borylated allylaziridine intermediate, Table 3 entry {entry}", [P3]))
    add_molecule(product_molecule(f"m_5{code}a", f"5{code}a", f"chlorinated allylaziridine 5{code}a", [P3]))
add_molecule(product_molecule("m_4aa_bz", "N-benzoylated 4aa", "N-benzoylated derivative of 4aa", [P3]))
add_molecule(product_molecule("m_6", "6", "cross-coupled derivative 6", [P3], alkene=True))
add_molecule(product_molecule("m_7", "7", "oxidized ketone derivative 7", [P3]))

bundle["molecules"] = list(molecules.values())
reactions = []

# Cited literature examples in Scheme 1.
for suffix, mids, title in [
    ("scheme1_hoveyda_literature", ("m_lit_hoveyda_allene", "m_lit_hoveyda_ketimine", "m_lit_hoveyda_product"), "Hoveyda's work"),
    ("scheme1_procter_literature", ("m_lit_procter_allene", "m_lit_procter_ketiminoester", "m_lit_procter_product"), "Procter's work"),
]:
    r = base_reaction(suffix, [P2], "literature_example", None, "Scheme 1, page 2", title,
                      "cited allene carboboronation with an imine electrophile")
    r["participants"] = [participant("allene", mids[0], "reactant", [P2]),
                         participant("imine", mids[1], "reactant", [P2]),
                         participant("diboron", "m_b2pin2", "reagent", [P2]),
                         participant("product", mids[2], "product", [P2], product_role="cited-work product")]
    r["steps"] = [step("cited_arrow", [P2], [p["participant_id"] for p in r["participants"]],
                       operation="Conditions are shown schematically for cited work; no present-paper experiment is asserted.")]
    r["classification"] = classification([P2], formed=[("C", "B"), ("C", "C")],
                                           limitations=["Cited generic R-group drawing; classification is category-level only."])
    r["review"]["unresolved"] = ["Generic cited-work structures and charges are not recoverable from the supplied main-text scheme."]
    reactions.append(r)

# Current-work generic schemes and proposed allylcopper pathway.
r = base_reaction("general_carboboronation", [P2], "generic_scheme", None,
                  "Scheme 1(3), page 2", "Copper-catalyzed three-component carboboronation of allenes",
                  "allene carboboronation with 2H-azirines")
r["participants"] = [participant("allene", "m_generic_allene", "reactant", [P2]),
                     participant("azirine", "m_generic_azirine", "reactant", [P2]),
                     participant("diboron", "m_b2pin2", "reagent", [P2]),
                     participant("copper", "m_cu_oa2", "catalyst", [P2]),
                     participant("ligand", "m_iprnhc_hcl", "ligand", [P2]),
                     participant("base", "m_kotbu", "reagent", [P2]),
                     participant("solvent", "m_phcf3", "solvent", [P2]),
                     participant("product", "m_generic_borylated_product", "product", [P2], product_role="main product family")]
r["classification"] = classification([P2], formed=[("C", "B"), ("C", "C"), ("N", "H")],
                                       decreased=[("C", "C", 2, 1), ("C", "N", 2, 1)])
r["review"]["unresolved"] = ["Generic R groups and atom mapping remain unresolved.",
                               "Scheme 1 names IPrNHC·HCl whereas tables/prose often abbreviate the ligand as IPrNHC."]
reactions.append(r)
carbo_scheme_id = r["reaction_id"]

r = base_reaction("general_chlorination", [P3], "generic_scheme", None,
                  "Tables 2-3 scheme header, page 3", "Chlorination of borylated allylaziridines",
                  "replacement of Bpin by chlorine")
r["participants"] = [participant("substrate", "m_generic_borylated_product", "reactant", [P3]),
                     participant("cucl2", "m_cucl2", "reagent", [P3]),
                     participant("meoh", "m_meoh", "solvent", [P3]),
                     participant("water", "m_water", "solvent", [P3]),
                     participant("product", "m_generic_chlorinated_product", "product", [P3], product_role="chlorinated product family")]
r["steps"] = [step("chlorination", [P3], [p["participant_id"] for p in r["participants"]],
                   operation="Treat with CuCl2 in CH3OH/H2O", temperature=("40 °C", {"number": 40, "unit": "°C"}))]
r["classification"] = classification([P3], formed=[("C", "Cl")], broken=[("C", "B")])
r["review"]["unresolved"] = ["CuCl2 charge, solvent ratio/volume, time, atmosphere, workup, and purification are not reported in the supplied main text."]
reactions.append(r)
chlor_scheme_id = r["reaction_id"]

r = base_reaction("scheme1_proposed_pathway", [P2], "proposed_mechanism", None,
                  "Scheme 1(3), lower pathway, page 2", "In situ allylcopper addition pathway", None)
r["participants"] = [participant("allene", "m_generic_allene", "reactant", [P2]),
                     participant("azirine", "m_generic_azirine", "reactant", [P2]),
                     participant("diboron", "m_b2pin2", "reagent", [P2]),
                     participant("product", "m_generic_borylated_product", "product", [P2])]
r["notes"] = [f("Scheme depicts borylcopper addition to allene, allylcopper resonance forms, azirine addition, and product formation.",
                "Author-proposed schematic pathway; not a discrete mechanistic experiment.", [P2])]
r["classification"] = no_classification([P2], "A proposed pathway is not an observed reaction experiment for bond-formation indexing.")
r["review"]["unresolved"] = ["Intermediate structures are generic and no mechanistic validation experiment is shown on the supplied pages."]
reactions.append(r)

# Optimization Table 1: every row is independent.
opt_rows = [
    (1, "m_cucl", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 76, 0.03),
    (2, "m_cucl", "m_dppf", "m_kotbu", "m_phcf3", 25, 58, 0.03),
    (3, "m_cucl", "m_dppbz", "m_kotbu", "m_phcf3", 25, 65, 0.03),
    (4, "m_cucl", "m_bipy", "m_kotbu", "m_phcf3", 25, 40, 0.03),
    (5, "m_cucl", "m_iprnhc", "m_naotbu", "m_phcf3", 25, 74, 0.03),
    (6, "m_cucl", "m_iprnhc", "m_liotbu", "m_phcf3", 25, 69, 0.03),
    (7, "m_cucl", "m_iprnhc", "m_kome", "m_phcf3", 25, 57, 0.03),
    (8, "m_cubr", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 61, 0.03),
    (9, "m_cui", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 66, 0.03),
    (10, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 83, 0.03),
    (11, "m_cu_otf2", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 61, 0.03),
    (12, "m_cubr2", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 57, 0.03),
    (13, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_toluene", 25, 80, 0.03),
    (14, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_thf", 25, 76, 0.03),
    (15, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_dioxane", 25, 72, 0.03),
    (16, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_phcf3", 0, 49, 0.03),
    (17, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_phcf3", 60, 73, 0.03),
    (18, "m_cu_oa2", "m_iprnhc", "m_kotbu", "m_phcf3", 25, 53, 0.02),
]
for entry, copper, ligand, base, solvent, temp, yld, ligand_amount in opt_rows:
    ev = [P2]
    r = base_reaction(f"table1_entry_{entry:02d}", ev, categories=["optimization"],
                      location=f"Table 1 entry {entry}, page 2",
                      name="Optimization of copper-catalyzed carboboronation",
                      transformation="three-component allene carboboronation")
    r["general_scheme_ids"] = [carbo_scheme_id]
    parts = [participant("azirine", "m_1a", "reactant", ev), participant("allene", "m_2a", "reactant", ev),
             participant("diboron", "m_b2pin2", "reagent", ev), participant("base", base, "reagent", ev),
             participant("copper", copper, "catalyst", ev), participant("ligand", ligand, "ligand", ev),
             participant("solvent", solvent, "solvent", ev),
             participant("product", "m_4aa", "product", ev, product_role="desired product")]
    for p, raw, val in [(parts[0], "0.2 mmol", {"number": 0.2, "unit": "mmol"}),
                        (parts[1], "0.6 mmol", {"number": 0.6, "unit": "mmol"}),
                        (parts[2], "0.4 mmol", {"number": 0.4, "unit": "mmol"}),
                        (parts[3], "0.2 mmol", {"number": 0.2, "unit": "mmol"}),
                        (parts[4], "0.02 mmol", {"number": 0.02, "unit": "mmol"}),
                        (parts[5], f"{ligand_amount:g} mmol", {"number": ligand_amount, "unit": "mmol"})]:
        charge(p, "amount", raw, val, ev)
    charge(parts[4], "mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, ev, "m_1a")
    charge(parts[5], "mol_percent", "15 mol %" if entry != 18 else "10 mol %",
           {"number": 15 if entry != 18 else 10, "unit": "mol%"}, ev, "m_1a")
    charge(parts[6], "volume", "3 mL", {"number": 3, "unit": "mL"}, ev)
    r["participants"] = parts
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f(molecules[solvent]["occurrences"][0]["compound_label"]["raw_value"], [solvent], ev)
    r["solvent_conditions"]["total_volume"] = f("3 mL", {"number": 3, "unit": "mL"}, ev)
    r["steps"] = [step("reaction", ev, [p["participant_id"] for p in parts[:-1]],
                       operation="Combine listed substrates/reagents under the Table 1 conditions",
                       temperature=(f"{temp} °C", {"number": temp, "unit": "°C"}),
                       time=("12 h", {"number": 12, "unit": "h"}))]
    r["procedure"]["original_text"] = f(
        "All reactions were carried out with 1a (0.2 mmol), 2a (0.6 mmol), 3 (0.4 mmol), KOtBu (0.2 mmol), copper salt (0.02 mmol), and ligand (0.03 mmol) in 3 mL solvent at 25 °C for 12 h; entries 16-18 have footnoted modifications.",
        "Table 1 footnote plus row-specific cells", ev)
    r["outcomes"] = [
        outcome("yield_4aa", "yield", ["m_4aa"], ev, f"{yld}", {"number": yld, "unit": "%"},
                product_role="desired product", yield_type="total isolated yield", basis="total isolated yield",
                scope="Table 1 yield of 4aa"),
        outcome("dr_4aa", "selectivity", ["m_4aa"], ev, ">20:1 dr",
                {"metric": "dr", "ratio": [20, 1], "qualifier": "gt"},
                product_role="desired product", basis="1H NMR analysis",
                scope="diastereomeric ratio; major/minor configuration not named"),
    ]
    r["classification"] = classification(ev, formed=[("C", "B"), ("C", "C"), ("N", "H")],
                                           decreased=[("C", "C", 2, 1), ("C", "N", 2, 1)])
    r["review"]["unresolved"] = ["Workup/purification are unavailable because SI was not supplied.",
                                  "Dense product structure and absolute configuration are not encoded."]
    reactions.append(r)


def standard_scope_parts(ev, azirine_mid, allene_mid, product4_mid=None, product5_mid=None):
    parts = [participant("azirine", azirine_mid, "reactant", ev), participant("allene", allene_mid, "reactant", ev),
             participant("diboron", "m_b2pin2", "reagent", ev), participant("base", "m_kotbu", "reagent", ev),
             participant("copper", "m_cu_oa2", "catalyst", ev), participant("ligand", "m_iprnhc", "ligand", ev),
             participant("phcf3", "m_phcf3", "solvent", ev), participant("cucl2", "m_cucl2", "reagent", ev),
             participant("meoh", "m_meoh", "solvent", ev), participant("water", "m_water", "solvent", ev)]
    if product4_mid:
        parts.append(participant("product4", product4_mid, "product", ev, product_role="isolated borylated intermediate"))
    if product5_mid:
        parts.append(participant("product5", product5_mid, "product", ev, product_role="two-step chlorinated product"))
    for p, raw, val in [(parts[0], "0.2 mmol", {"number": 0.2, "unit": "mmol"}),
                        (parts[1], "0.6 mmol", {"number": 0.6, "unit": "mmol"}),
                        (parts[2], "0.4 mmol", {"number": 0.4, "unit": "mmol"}),
                        (parts[3], "0.2 mmol", {"number": 0.2, "unit": "mmol"})]:
        charge(p, "amount", raw, val, ev)
    charge(parts[4], "mol_percent", "10 mol %", {"number": 10, "unit": "mol%"}, ev, azirine_mid)
    parts[5]["charges"]["mol_percent"] = f(
        "15 mmol %", {"number": 15, "unit": "mol%"}, ev, method="inference",
        derivation="The Table 2/Table 3 footnote literally prints '15 mmol %'. It is retained raw and normalized as 15 mol% only because the field is a catalyst loading and page-2 optimization/prose explicitly use 15 mol %; this remains a machine inference.")
    parts[5]["equivalence_reference"] = f(azirine_mid, azirine_mid, ev)
    charge(parts[6], "volume", "3 mL", {"number": 3, "unit": "mL"}, ev)
    return parts


def fill_scope_record(r, ev, azirine_mid, product4_mid, product5_mid, y4, y5):
    parts = standard_scope_parts(ev, azirine_mid, r["participants_seed_allene"], product4_mid, product5_mid)
    r.pop("participants_seed_allene")
    r["participants"] = parts
    r["reaction_scale"] = f("azirine (0.2 mmol)", {"reference_molecule_id": azirine_mid,
                                                  "amount": {"number": 0.2, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = f("PhCF3; subsequent CH3OH/H2O", ["m_phcf3", "m_meoh", "m_water"], ev)
    r["solvent_conditions"]["total_volume"] = f("3 mL PhCF3; chlorination volumes unreported",
                                                  {"number": 3, "unit": "mL"}, ev)
    r["steps"] = [
        step("carboboronation", ev, ["azirine", "allene", "diboron", "base", "copper", "ligand", "phcf3"],
             operation="Copper-catalyzed three-component carboboronation",
             temperature=("25 °C", {"number": 25, "unit": "°C"}), time=("12 h", {"number": 12, "unit": "h"})),
        step("chlorination", ev, ["product4", "cucl2", "meoh", "water"],
             operation="Convert borylated product 4 to chlorinated product 5 with CuCl2 in CH3OH/H2O",
             temperature=("40 °C", {"number": 40, "unit": "°C"})),
    ]
    r["procedure"]["original_text"] = f(
        "Scope footnote: azirine (0.2 mmol), allene 2 (0.6 mmol), B2(pin)2 3 (0.4 mmol), tBuOK (0.2 mmol), Cu(OAc)2 (10 mol %), IPrNHC (15 mmol %) in 3 mL PhCF3 at 25 °C for 12 h; product 4 then treated with CuCl2 in CH3OH/H2O at 40 °C.",
        "Table scope footnote and scheme arrows", ev)
    r["notes"].append(f("IPrNHC (15 mmol %)",
                        "The malformed source unit is retained raw; 15 mol% is a context-based inferred normalization, not a silent correction.", ev))
    r["outcomes"] = [
        outcome("yield_product4", "yield", [product4_mid], ev, str(y4), {"number": y4, "unit": "%"},
                product_role="borylated intermediate", yield_type="isolated yield", basis="isolated yield",
                scope="first-step yield of product family 4"),
        outcome("yield_product5", "yield", [product5_mid], ev, str(y5), {"number": y5, "unit": "%"},
                product_role="chlorinated product", yield_type="total isolated yield", basis="total isolated yield",
                scope="two-step total yield from azirine/allene to product 5"),
        outcome("diastereoselectivity", "selectivity", [product4_mid], ev, "complete diastereoselectivity",
                "complete diastereoselectivity", product_role="borylated intermediate",
                scope="qualitative author statement; configurations and numeric ratio are not stated for scope entries"),
    ]
    r["classification"] = classification(ev, formed=[("C", "B"), ("C", "C"), ("N", "H"), ("C", "Cl")],
                                           decreased=[("C", "C", 2, 1), ("C", "N", 2, 1)],
                                           broken=[("C", "B")])
    r["review"]["unresolved"] = ["Chlorination reagent charge, solvent ratio/volume, time, atmosphere, workup, and purification are unreported in the supplied main text.",
                                  "Scope dr is qualitative; its measurement method and major/minor configurations are unreported.",
                                  "Dense product structures retain source drawings without complete stereochemical SMILES."]


# Table 2 scope rows.
for entry, code, y4, y5 in t2:
    ev = [P3]
    p4 = "m_4aa" if code == "a" else f"m_t2e{entry}_4"
    p5 = "m_5aa" if code == "a" else f"m_5a{code}"
    r = base_reaction(f"table2_entry_{entry:02d}", ev, categories=["substrate_scope"],
                      location=f"Table 2 entry {entry}, page 3",
                      name="Scope of allenes", transformation="carboboronation followed by chlorination")
    r["general_scheme_ids"] = [carbo_scheme_id, chlor_scheme_id]
    r["participants_seed_allene"] = f"m_2{code}"
    fill_scope_record(r, ev, "m_1a", p4, p5, y4, y5)
    reactions.append(r)

# Three explicitly tested but unsuccessful allene classes in the prose below Table 3.
for code in ("n", "o", "p"):
    ev = [P3, T3]
    r = base_reaction(f"allene_2{code}_no_reaction", ev, categories=["substrate_scope", "control"],
                      location=f"Page 3 prose, allene 2{code}", name="Unsuccessful allene scope test",
                      transformation="attempted carboboronation")
    r["general_scheme_ids"] = [carbo_scheme_id]
    parts = standard_scope_parts(ev, "m_1a", f"m_2{code}")[:7]
    for p in parts:
        for q in p["charges"].values():
            if q["normalized_value"] is not None and q["acquisition_method"] == "original_extraction":
                q["acquisition_method"] = "inherited_general_condition"
    r["participants"] = parts
    r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": "m_1a",
                                              "amount": {"number": 0.2, "unit": "mmol"}}, ev,
                                    method="inherited_general_condition")
    r["solvent_conditions"]["composition"] = f("optimal conditions / PhCF3", ["m_phcf3"], ev,
                                                   method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = f("3 mL", {"number": 3, "unit": "mL"}, ev,
                                                    method="inherited_general_condition")
    r["steps"] = [step("attempted_carboboronation", ev, [p["participant_id"] for p in parts],
                       operation="Subjected to the optimal reaction conditions",
                       temperature=("25 °C (optimal conditions)", {"number": 25, "unit": "°C"}),
                       time=("12 h (optimal conditions)", {"number": 12, "unit": "h"}))]
    # Explicitly mark inherited fields.
    for sfield in (r["steps"][0]["temperature"], r["steps"][0]["time"]):
        sfield["acquisition_method"] = "inherited_general_condition"
    r["outcomes"] = [outcome("no_desired_product", "no_reaction", [], ev,
                             "no desired products was detected", "no desired product detected",
                             scope=f"attempted reaction of allene 2{code}; detection limit/method unreported")]
    r["classification"] = no_classification(ev, "No desired product was detected; intended bonds are not classified as observed formation.")
    r["review"]["unresolved"] = ["Exact structure of the text-only allene label is unavailable on supplied pages.",
                                  "Detection method/limit, side products, workup, and purification are unreported."]
    reactions.append(r)

# Table 3 successful azirine scope rows.
for entry, code, y4, y5 in t3:
    ev = [P3]
    p4 = f"m_t3e{entry}_4"
    p5 = f"m_5{code}a"
    r = base_reaction(f"table3_entry_{entry:02d}", ev, categories=["substrate_scope"],
                      location=f"Table 3 entry {entry}, page 3", name="Scope of 2H-azirines",
                      transformation="carboboronation followed by chlorination")
    r["general_scheme_ids"] = [carbo_scheme_id, chlor_scheme_id]
    r["participants_seed_allene"] = "m_2a"
    fill_scope_record(r, ev, f"m_1{code}", p4, p5, y4, y5)
    reactions.append(r)

# Table 3 entry 11: no product and rapid substrate decomposition.
ev = [P3, T3]
r = base_reaction("table3_entry_11", ev, categories=["substrate_scope", "control"],
                  location="Table 3 entry 11 and adjacent prose, page 3", name="Scope of 2H-azirines",
                  transformation="attempted carboboronation")
r["general_scheme_ids"] = [carbo_scheme_id]
parts = standard_scope_parts(ev, "m_1l", "m_2a")[:7]
r["participants"] = parts
r["reaction_scale"] = f("1l (0.2 mmol)", {"reference_molecule_id": "m_1l", "amount": {"number": 0.2, "unit": "mmol"}}, ev)
r["solvent_conditions"]["composition"] = f("PhCF3", ["m_phcf3"], ev)
r["solvent_conditions"]["total_volume"] = f("3 mL", {"number": 3, "unit": "mL"}, ev)
r["steps"] = [step("attempted_carboboronation", ev, [p["participant_id"] for p in parts],
                   operation="Table 3 standard conditions", temperature=("25 °C", {"number": 25, "unit": "°C"}),
                   time=("12 h", {"number": 12, "unit": "h"}))]
r["outcomes"] = [
    outcome("no_desired_product", "no_reaction", [], ev, "no desired product was obtained",
            "no desired product obtained", scope="Table 3 entry 11; intended label 5la is printed but not observed"),
    outcome("substrate_decomposition", "decomposition", ["m_1l"], ev,
            "substrate decomposed rapidly under standard conditions", "rapid decomposition",
            scope="2H-azirine 1l (R1=Ph, R2=H)"),
]
r["classification"] = no_classification(ev, "No desired product was observed; decomposition products are not identified.")
r["review"]["unresolved"] = ["Decomposition products, rate, detection method/limit, and workup are unreported.",
                              "The table prints intended target label 5la, but no product molecule is asserted as observed."]
reactions.append(r)

# Scheme 2 arrow 1: scale-up carboboronation.
ev = [P3]
r = base_reaction("scheme2_scaleup_carboboronation", ev, categories=["scale_up"],
                  location="Scheme 2(1), first arrow, page 3", name="Scaled-up synthesis of 4aa",
                  transformation="three-component allene carboboronation")
r["general_scheme_ids"] = [carbo_scheme_id]
parts = [participant("azirine", "m_1a", "reactant", ev), participant("allene", "m_2a", "reactant", ev),
         participant("diboron", "m_b2pin2", "reagent", ev), participant("copper", "m_cu_oa2", "catalyst", ev),
         participant("ligand", "m_iprnhc", "ligand", ev), participant("base", "m_kotbu", "reagent", ev),
         participant("solvent", "m_phcf3", "solvent", ev),
         participant("product", "m_4aa", "product", ev, product_role="isolated scale-up product")]
charge(parts[0], "mass", "193 mg", {"number": 193, "unit": "mg"}, ev)
charge(parts[0], "amount", "1 mmol", {"number": 1, "unit": "mmol"}, ev)
charge(parts[1], "mass", "348 mg", {"number": 348, "unit": "mg"}, ev)
charge(parts[1], "amount", "3 mmol", {"number": 3, "unit": "mmol"}, ev)
r["participants"] = parts
r["reaction_scale"] = f("1a: 193 mg (1 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 1, "unit": "mmol"}}, ev)
r["solvent_conditions"]["composition"] = f("PhCF3", ["m_phcf3"], ev)
r["steps"] = [step("scaleup_carboboronation", ev, [p["participant_id"] for p in parts[:-1]],
                   operation="Scaled-up carboboronation under the depicted optimal reagent set",
                   temperature=("r.t.", "r.t."))]
r["outcomes"] = [
    outcome("yield_4aa", "yield", ["m_4aa"], ev, "61% yield", {"number": 61, "unit": "%"},
            product_role="scale-up product", yield_type="yield", scope="1a to 4aa on 1 mmol scale"),
    outcome("amount_4aa", "product_amount", ["m_4aa"], ev, "267 mg", {"number": 267, "unit": "mg"},
            product_role="scale-up product", scope="isolated/recovered product mass as printed; isolation method not stated"),
]
r["classification"] = classification(ev, formed=[("C", "B"), ("C", "C"), ("N", "H")],
                                       decreased=[("C", "C", 2, 1), ("C", "N", 2, 1)])
r["review"]["unresolved"] = ["B2(pin)2/base/catalyst/ligand charges, solvent volume, time, atmosphere, workup, purification, yield basis, and dr method are not printed for the scale-up arrow."]
reactions.append(r)

# Scheme 2 arrow 2: chlorination, retaining the printed two-step basis.
r = base_reaction("scheme2_scaleup_chlorination", ev, categories=["scale_up", "derivatization"],
                  location="Scheme 2(1), second arrow, page 3", name="Scale-up chlorination of 4aa",
                  transformation="replacement of Bpin by chlorine")
r["general_scheme_ids"] = [chlor_scheme_id]
r["participants"] = [participant("substrate", "m_4aa", "reactant", ev), participant("cucl2", "m_cucl2", "reagent", ev),
                     participant("meoh", "m_meoh", "solvent", ev), participant("water", "m_water", "solvent", ev),
                     participant("product", "m_5aa", "product", ev, product_role="chlorinated scale-up product")]
r["solvent_conditions"]["composition"] = f("CH3OH/H2O", ["m_meoh", "m_water"], ev)
r["steps"] = [step("chlorination", ev, ["substrate", "cucl2", "meoh", "water"],
                   operation="Chlorination with CuCl2 in CH3OH/H2O", temperature=("40 °C", {"number": 40, "unit": "°C"}))]
r["outcomes"] = [
    outcome("yield_5aa", "yield", ["m_5aa"], ev, "58% yield in two steps", {"number": 58, "unit": "%"},
            product_role="chlorinated product", yield_type="two-step overall yield",
            scope="overall from 1a through carboboronation and chlorination; not a standalone chlorination yield"),
    outcome("amount_5aa", "product_amount", ["m_5aa"], ev, "201 mg", {"number": 201, "unit": "mg"},
            product_role="chlorinated product", scope="printed final product mass"),
]
r["classification"] = classification(ev, formed=[("C", "Cl")], broken=[("C", "B")])
r["review"]["unresolved"] = ["4aa/CuCl2 charges, solvent ratio/volume, reaction time, atmosphere, workup, purification, and yield measurement basis are unreported.",
                              "The 58% value is explicitly a two-step overall yield and must not be treated as a standalone chlorination yield."]
reactions.append(r)

# Scheme 2(2) protection arrow.
r = base_reaction("scheme2_protection", ev, categories=["derivatization"],
                  location="Scheme 2(2), protection arrow, page 3", name="Benzoyl protection of 4aa",
                  transformation="N-benzoylation")
r["participants"] = [participant("substrate", "m_4aa", "reactant", ev),
                     participant("benzoyl_chloride", "m_benzoyl_chloride", "reagent", ev),
                     participant("triethylamine", "m_triethylamine", "reagent", ev),
                     participant("dmap", "m_dmap", "catalyst", ev), participant("dcm", "m_dcm", "solvent", ev),
                     participant("product", "m_4aa_bz", "product", ev, product_role="protected intermediate")]
r["solvent_conditions"]["composition"] = f("DCM", ["m_dcm"], ev)
r["steps"] = [step("benzoyl_protection", ev, [p["participant_id"] for p in r["participants"][:-1]],
                   operation="React 4aa with benzoyl chloride in the presence of Et3N and DMAP",
                   temperature=("rt", "rt"))]
r["outcomes"] = [outcome("protected_intermediate", "other", ["m_4aa_bz"], ev,
                         "protected product", "protected product obtained", product_role="protected intermediate",
                         scope="No standalone protection yield is printed")]
r["classification"] = classification(ev, formed=[("C", "N")], broken=[("N", "H")])
r["review"]["unresolved"] = ["The protected intermediate is also labeled 4aa in the scheme/prose; a distinct internal molecule ID is used to avoid conflation.",
                              "All charges, time, atmosphere, workup, purification, and standalone yield are unreported."]
reactions.append(r)

# Scheme 2(2a) cross-coupling arrow.
r = base_reaction("scheme2_cross_coupling", ev, categories=["derivatization"],
                  location="Scheme 2(2a), second arrow, page 3", name="Palladium-catalyzed cross-coupling to 6",
                  transformation="C-C cross-coupling of organoboron derivative")
r["participants"] = [participant("substrate", "m_4aa_bz", "reactant", ev), participant("iodobenzene", "m_iodobenzene", "reactant", ev),
                     participant("palladium", "m_pd_oa2", "catalyst", ev), participant("sphos", "m_sphos", "ligand", ev),
                     participant("naoh", "m_naoh", "reagent", ev), participant("thf", "m_thf", "solvent", ev),
                     participant("water", "m_water", "solvent", ev), participant("product", "m_6", "product", ev, product_role="cross-coupled product")]
r["solvent_conditions"]["composition"] = f("THF/H2O", ["m_thf", "m_water"], ev)
r["steps"] = [step("cross_coupling", ev, [p["participant_id"] for p in r["participants"][:-1]],
                   operation="Pd(OAc)2/SPhos/NaOH coupling with iodobenzene",
                   temperature=("60 °C", {"number": 60, "unit": "°C"}))]
r["outcomes"] = [outcome("yield_6", "yield", ["m_6"], ev, "70% yield in two steps", {"number": 70, "unit": "%"},
                         product_role="cross-coupled product", yield_type="two-step overall yield",
                         scope="overall from 4aa through protection and cross-coupling to 6")]
r["classification"] = classification(ev, formed=[("C", "C")], broken=[("C", "B")])
r["review"]["unresolved"] = ["All charges, solvent ratio/volume, time, atmosphere, workup, purification, yield basis, and alkene geometry assignment are unreported.",
                              "The 70% value is a two-step overall yield from 4aa, not a standalone coupling yield."]
reactions.append(r)

# Scheme 2(2b) oxidation arrow.
r = base_reaction("scheme2_oxidation", ev, categories=["derivatization"],
                  location="Scheme 2(2b), second arrow, page 3", name="Oxidation to ketone 7",
                  transformation="oxidation of organoboron derivative to ketone")
r["participants"] = [participant("substrate", "m_4aa_bz", "reactant", ev),
                     participant("oxidant", "m_nabo3_4h2o", "reagent", ev), participant("thf", "m_thf", "solvent", ev),
                     participant("water", "m_water", "solvent", ev), participant("product", "m_7", "product", ev, product_role="oxidized ketone product")]
r["solvent_conditions"]["composition"] = f("THF/H2O", ["m_thf", "m_water"], ev)
r["steps"] = [step("oxidation", ev, [p["participant_id"] for p in r["participants"][:-1]],
                   operation="Oxidative workup with NaBO3·4H2O in THF/H2O")]
r["outcomes"] = [outcome("yield_7", "yield", ["m_7"], ev, "74% yield in two steps", {"number": 74, "unit": "%"},
                         product_role="oxidized ketone product", yield_type="two-step overall yield",
                         scope="overall from 4aa through protection and oxidation to 7")]
r["classification"] = classification(ev, formed=[("C", "O")], broken=[("C", "B")])
# Correct the formed C-O change to the source-drawn carbonyl double bond and exact taxonomy label.
carbonyl_change = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
                   "before_order": 0, "after_order": 2, "change_type": "formed"}
r["classification"]["bond_changes"][0] = derived(carbonyl_change, ev, "Source drawing shows a newly formed carbonyl; atom mapping is unresolved.")
r["classification"]["labels"] = derived([{"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["b1"]}], ev,
                                          "Exact taxonomy projection of the source-drawn carbonyl formation.")
r["review"]["unresolved"] = ["All charges, solvent ratio/volume, temperature, time, atmosphere, workup, purification, and yield basis are unreported.",
                              "The 74% value is a two-step overall yield from 4aa, not a standalone oxidation yield."]
reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:graphical_abstract", "document_id": DOC, "kind": "duplicate_overview",
     "evidence_ids": [P1, T1], "note": "Graphical abstract duplicates the generic carboboronation/chlorination and is not a separate experiment."},
    {"context_id": f"{DOC}:nmr_characterization_limitation", "document_id": DOC, "kind": "characterization_note",
     "evidence_ids": [P2, T2], "note": "Most borylated products 4 had poor 1H splitting and missing 13C signals; chlorination was used for easier NMR analysis. No standalone spectrum is treated as a synthesis."},
    {"context_id": f"{DOC}:asymmetric_screen_summary", "document_id": DOC, "kind": "unresolved_experiment_summary",
     "evidence_ids": [P3, T3, P5, T5], "note": "The asymmetric version is reported only as poor yield and moderate enantioselectivity; exact screening experiments are deferred to unavailable SI and cannot be separated into records."},
    {"context_id": f"{DOC}:si_notice", "document_id": DOC, "kind": "supporting_information_notice",
     "evidence_ids": [P4, T4], "note": "The paper states that SI exists, but no SI file was supplied or explicitly associated in this session."},
    {"context_id": f"{DOC}:ccdc_2058388", "document_id": DOC, "kind": "crystallography_notice",
     "evidence_ids": [P2, P4, T4], "note": "CCDC 2058388 is cited for 4aa crystal data; the supplied pages do not print a complete absolute stereochemical assignment."},
    {"context_id": f"{DOC}:references_pages_4_5", "document_id": DOC, "kind": "bibliography",
     "evidence_ids": [P4, P5, T4, T5], "note": "Reference list only; cited chemistry is not counted as present-paper experiments except the two Scheme 1 arrows explicitly represented as literature examples."},
]

transcription = {
    "document_id": DOC,
    "pages_reviewed": [1, 2, 3, 4, 5],
    "table1_rows": opt_rows,
    "table2_rows": t2,
    "table3_success_rows": t3,
    "table3_entry11": {"R1": "C6H5", "R2": "H", "intended_label": "5la", "result": "no desired product; rapid decomposition"},
    "negative_allenes": ["2n alpha-substituted phenylallene", "2o phenoxyallene", "2p allenamide"],
    "scheme2_arrows": ["scale-up carboboronation", "scale-up chlorination", "benzoyl protection", "cross-coupling", "oxidation"],
    "limitations": ["No SI supplied", "dense structures not converted to guessed SMILES", "machine records remain needs_review"],
}

(OUT / "raw_model/transcription.json").write_text(json.dumps(transcription, ensure_ascii=False, indent=2) + "\n")
diagnostics = {m["molecule_id"]: st.stereo_inspect(st.value(m["smiles"]))
               for m in bundle["molecules"] if m["representation_type"] == "smiles"}
(OUT / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
(OUT / "raw_model/source_mentions.json").write_text(json.dumps(
    st.source_mentions({"documents": documents}, ["1a", "2a", "4aa", "5aa", "2n", "2o", "2p", "5la"]),
    ensure_ascii=False, indent=2) + "\n")
(OUT / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
current_dataset = json.loads((OUT / "dataset.json").read_text())
if current_dataset.get("experiments"):
    (OUT / "audit.json").write_text(json.dumps(rt.dataset_audit(current_dataset), ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
