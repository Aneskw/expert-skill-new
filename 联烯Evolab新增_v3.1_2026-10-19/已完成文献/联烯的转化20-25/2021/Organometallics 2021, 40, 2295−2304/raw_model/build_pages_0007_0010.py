#!/usr/bin/env python3
"""Build the incremental v3.1 bundle for source PDF pages 7--10.

The script merges detailed experimental prose into already-created Scheme
records by reaction_id and creates only genuinely distinct experiments from
the newly inspected pages.  All outputs remain machine records needing review.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_23c0fdb9606e2366"
DATASET = json.loads((RUN / "dataset.json").read_text())
OLD = json.loads((RUN / "standard_records.json").read_text())


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def union(items, additions):
    return list(dict.fromkeys(list(items) + list(additions)))


bundle = st.empty(DATASET["documents"])
old_molecules = {m["molecule_id"]: copy.deepcopy(m) for m in OLD["molecules"]}
old_reactions = {r["reaction_id"]: copy.deepcopy(r) for r in OLD["reactions"]}
molecules = old_molecules
reactions = []
stereo_diagnostics = {}


def touch_molecule(mid, evidence):
    m = molecules[mid]
    occurrence = next(o for o in m["occurrences"] if o["document_id"] == DOC)
    occurrence["evidence_ids"] = union(occurrence["evidence_ids"], evidence)
    occurrence["compound_label"]["evidence_ids"] = union(
        occurrence["compound_label"]["evidence_ids"], evidence)
    for key in ("original", "english", "abbreviations"):
        fact = m["names"][key]
        if fact["normalized_value"] is not None:
            fact["evidence_ids"] = union(fact["evidence_ids"], evidence)
    return m


def new_molecule(mid, name, label, evidence, smiles=None, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded",
                 identity_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(name, name, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F(label, [label], evidence) if label else F()
    m["names"]["aliases"] = F()
    m["cas_numbers"] = F()
    m["external_properties"] = []
    if smiles:
        diag = st.stereo_inspect(smiles)
        stereo_diagnostics[mid] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = F(None, smiles, evidence, method="structure_derivation",
                         derivation="Candidate connectivity encoded from the legible source name/drawing; parsing is not source verification.")
        m["stereochemistry"]["atom_centers"] = F(
            None, diag["atom_centers"], evidence, method="structure_derivation",
            derivation="RDKit diagnostic for the candidate SMILES; source configuration remains separately stated.")
        m["stereochemistry"]["double_bonds"] = F(
            None, diag["double_bonds"], evidence, method="structure_derivation",
            derivation="RDKit double-bond diagnostic; no missing source E/Z assignment is inferred.")
        m["stereochemistry"]["isomeric_smiles"] = F(
            None, diag["canonical_isomeric_smiles"], evidence, method="structure_derivation",
            derivation="Canonical isomeric graph generated without inventing unreported stereochemistry.")
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="Stereo-preserving graph key for the candidate SMILES; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        for key in ("smiles",):
            m[key] = F()
        for key in ("atom_centers", "double_bonds", "isomeric_smiles", "identity_key"):
            m["stereochemistry"][key] = F()
    m["stereochemistry"]["axial"] = F()
    m["stereochemistry"]["mixture"] = F()
    if chirality == "unspecified":
        m["stereochemistry"]["assignment"] = F(
            None, "unassigned", evidence, method="inference",
            derivation="The inspected source supplies no absolute or relative configuration; this is unassigned, not achiral or racemic.")
    else:
        m["stereochemistry"]["assignment"] = F()
    for key, value in (("chirality_status", chirality),
                       ("geometry_status", geometry),
                       ("representation_status", representation)):
        m["stereochemistry"][key] = F(
            None, value, evidence, method="inference",
            derivation="Status records the stereochemical information visible in the inspected source and candidate representation.")
    m["identity_note"] = F(
        None,
        identity_note or ("Candidate structure is machine-derived from the source name/drawing and remains needs_review."
                          if smiles else "Identity is retained from the source wording; no unique molecular representation is asserted."),
        evidence, method="inference",
        derivation="Extraction note describing identity and representation limits.")
    molecules[mid] = m
    return mid


def new_generic_molecule(mid, name, label, fixed, reason, evidence):
    new_molecule(mid, name, label, evidence, identity_note=reason)
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, evidence),
        "variables": [],
        "instance_mapping": [],
        "reason": F(None, reason, evidence, method="inference",
                    derivation="The general procedure names a chemical class rather than one unique compound.")
    }
    return mid


# Add page evidence to existing identities actually used in the new page range.
page7_mids = ["m_7a", "m_5", "m_6", "m_6p", "m_3", "m_3a", "m_3b", "m_3c", "m_3d",
              "m_nbs", "m_can", "m_k2co3", "m_etoh", "m_i2", "m_tfa", "m_ch2cl2",
              "m_chcl3", "m_fpstar_i", "m_6a", "m_7b", "m_9a", "m_9b", "m_10a", "m_10b",
              "m_dce", "m_alumina", "m_silica", "m_hexanes", "m_etoac"]
page8_mids = ["m_7b", "m_9c", "m_9d", "m_10c", "m_10d", "m_11", "m_12", "m_13", "m_14",
              "m_bf3", "m_ch2cl2", "m_hexanes", "m_etoac"]
for mid in page7_mids:
    touch_molecule(mid, [ev(7)])
for mid in page8_mids:
    touch_molecule(mid, [ev(8)])
touch_molecule("m_1a", [ev(10)])
touch_molecule("m_bf3", [ev(10)])
touch_molecule("m_dce", [ev(10)])

# The source drawing of 1a has differently substituted allene termini but gives
# no R_a/S_a or P/M descriptor.  Preserve the unassigned axis explicitly.
m1a = molecules["m_1a"]
m1a["stereochemistry"]["axial"] = F(
    "allene structure 1a shown without an axial descriptor",
    {"type": "allene_axis", "assignment": "unassigned"}, [ev(10, "image")],
    method="structure_derivation",
    derivation="Both allene termini are differently substituted in the displayed 1a drawing; no axial descriptor is printed.")
m1a["stereochemistry"]["chirality_status"] = F(
    None, "unspecified", [ev(10, "image")], method="inference",
    derivation="The displayed allene can support axial chirality, but the source does not specify an axial configuration.")
m1a["stereochemistry"]["representation_status"] = F(
    None, "not_encoded", [ev(10, "image")], method="inference",
    derivation="The organometallic allene and its unassigned axis are not encoded as a concrete SMILES.")

# DIAD contains a potentially stereogenic azo double bond; the inspected source
# gives no E/Z assignment, so it cannot be marked geometry-not-applicable.
m13 = molecules["m_13"]
m13["stereochemistry"]["geometry_status"] = F(
    None, "unspecified", [ev(8)], method="inference",
    derivation="The candidate structure has an N=N bond and the source does not assign E/Z geometry.")
m13["stereochemistry"]["representation_status"] = F(
    None, "not_encoded", [ev(8)], method="inference",
    derivation="Connectivity is encoded, but the unreported azo E/Z geometry is not.")

new_molecule("m_br2", "bromine", "Br2", [ev(7)], "BrBr", "achiral", "not_applicable", "not_applicable")
new_molecule("m_basic_alumina", "basic aluminum oxide", "basic aluminum oxide", [ev(7)])
new_molecule(
    "m_p7_alt6a_start_literal",
    "compound 6a as literally printed as the starting material in the 1.0-equivalent comparison",
    "6a (printed starting-material label)", [ev(7)],
    identity_note="The comparison sentence literally says '1.0 equiv each of 6a and I2'; context may have intended complex 6, but the machine record does not silently correct the source.")
new_generic_molecule(
    "m_generic_allenyliron_pyrazole_gp", "freshly prepared allenyliron complex (general procedure)",
    "allenyliron complex", "cyclopentadienyliron-derived eta1-allenyl complex; substituents not fixed in the general procedure",
    "The page-7 general procedure applies to several concrete complexes (7b and 11) and is not one unique molecule.", [ev(7)])
new_generic_molecule(
    "m_generic_diazonium_pyrazole_gp", "diazonium tetrafluoroborate (general procedure)",
    "diazonium tetrafluoroborate", "arenediazonium tetrafluoroborate class; aryl substituent not fixed in the general procedure",
    "The page-7 general procedure applies to several diazonium salts and is not one unique molecule.", [ev(7)])
new_molecule(
    "m_benzaldehyde_tosylimine", "benzaldehyde tosylimine", "benzaldehyde tosylimine", [ev(10)],
    "Cc1ccc(S(=O)(=O)N=Cc2ccccc2)cc1", "achiral", "unspecified", "not_encoded",
    "Candidate connectivity follows the footnote wording and displayed imine; C=N geometry is not assigned in the source.")
new_molecule(
    "m_eq1_dihydropyrrole", "N-tosyl-4,5-diphenyl-1H-2,5-dihydropyrrole",
    "N-tosyl-4,5-diphenyl-1H-2,5-dihydropyrrole", [ev(10)],
    "Cc1ccc(S(=O)(=O)N2CC=C(c3ccccc3)C2c2ccccc2)cc1",
    "unspecified", "not_applicable", "not_encoded",
    "Candidate connectivity follows the explicit product name and drawing; the potential C5 stereocenter is not assigned.")


def charge_fact(raw, normalized, evidence, method="original_extraction", derivation=None,
                status=None, alternatives=()):
    return F(raw, normalized, evidence, method=method, derivation=derivation,
             status=status, alternatives=alternatives)


def participant(pid, mid, role, evidence, stage="reaction", product_role=None,
                charges=None, reference=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = F(product_role, product_role, evidence)
    for key, spec in (charges or {}).items():
        if isinstance(spec, dict) and "raw_value" in spec:
            p["charges"][key] = spec
        else:
            raw, normalized = spec[:2]
            method = spec[2] if len(spec) > 2 else "original_extraction"
            derivation = spec[3] if len(spec) > 3 else None
            p["charges"][key] = F(raw, normalized, evidence, method=method, derivation=derivation)
    if reference:
        p["equivalence_reference"] = F(reference, reference, evidence)
    return p


def step(step_id, evidence, pids, operation, temperature=None, time=None,
         atmosphere=None, pressure=None, stage="reaction", method="original_extraction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence, method=method,
                         derivation=("Normalized operation reconstructed from the cited source passage."
                                     if method == "inference" else None))
    if temperature:
        s["temperature"] = F(temperature[0], temperature[1], evidence,
                               method=temperature[2] if len(temperature) > 2 else "original_extraction")
    if time:
        s["time"] = F(time[0], time[1], evidence,
                       method=time[2] if len(time) > 2 else "original_extraction")
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, evidence)
    if pressure:
        s["pressure"] = F(pressure[0], pressure[1], evidence)
    return s


def outcome(oid, kind, targets, evidence, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = F(raw, normalized, evidence)
    if product_role:
        o["product_role"] = F(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, evidence)
    if basis:
        braw, bnorm, bev = basis
        o["measurement_basis"] = F(braw, bnorm, bev)
    if scope:
        o["scope"] = F(scope, scope, evidence)
    return o


ISO_BASIS = ("Yields refer to isolated compounds, unless otherwise indicated.", "isolated", [ev(6)])


def classify(r, changes=None, unresolved=False, not_applicable=False):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F()
    c["model_votes"] = []
    if unresolved:
        c["bond_changes"] = []
        c["labels"] = F()
        c["coverage"] = F(None, "unresolved", r["evidence_ids"], method="inference",
                            derivation="Generic/proposed structures were not atom-mapped.")
    elif not_applicable:
        c["bond_changes"] = []
        c["labels"] = F(None, [], r["evidence_ids"], method="inference",
                          derivation="No observed transformation is classified.")
        c["coverage"] = F(None, "not_applicable", r["evidence_ids"], method="inference",
                            derivation="No observed transformation is classified.")
    else:
        facts = []
        labels = {}
        for idx, (a, b, before, after, kind) in enumerate(changes or [], 1):
            cid = f"b{idx}"
            value = {"change_id": cid,
                     "atoms": [{"atom_map": None, "element": a},
                               {"atom_map": None, "element": b}],
                     "before_order": before, "after_order": after,
                     "change_type": kind}
            facts.append(F(None, value, r["evidence_ids"], method="structure_derivation",
                           derivation="Partial atom-type change inferred from the displayed transformation; no complete atom mapping is asserted."))
            if kind in ("formed", "order_increased", "hydrogen_gain"):
                lab = st.bond_label(a, b, after)
                labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(cid)
        c["bond_changes"] = facts
        c["labels"] = F(None, list(labels.values()), r["evidence_ids"], method="structure_derivation",
                          derivation="Taxonomy labels calculated from the retained partial bond changes.")
        c["coverage"] = F(None, "partial", r["evidence_ids"], method="inference",
                            derivation="Only source-supportable net changes are retained; no complete atom/coproduct mapping is claimed.")
    c["limitations"] = [
        "No complete atom mapping was asserted for the coordination-rich or incompletely specified transformation.",
        "Hydrogen and coproduct sources remain unresolved unless explicitly reported."
    ]


def reaction(rid, evidence, record_type="specific_experiment", categories=None,
             location=None, name=None, transform=None, participants=None, steps=None,
             outcomes=None, notes=None, procedure=None, general_scheme_ids=None,
             solvents=None, reactor=None, changes=None, class_unresolved=False,
             unresolved=None):
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, evidence) if categories is not None else F()
    r["location"] = F(location, location, evidence) if location else F()
    r["reaction_name"] = F(name, name, evidence) if name else F()
    r["transformation_type"] = F(None, transform, evidence, method="inference",
                                   derivation="Concise transformation description inferred from source wording/drawing.") if transform else F()
    r["general_scheme_ids"] = list(general_scheme_ids or [])
    r["participants"] = list(participants or [])
    r["steps"] = list(steps or [])
    r["outcomes"] = list(outcomes or [])
    r["notes"] = [F(None, text, evidence, method="inference",
                    derivation="Extraction note preserving a source limitation or conflict.")
                  for text in (notes or [])]
    if solvents:
        raw, normalized, total = solvents
        r["solvent_conditions"]["composition"] = F(raw, normalized, evidence)
        if total:
            r["solvent_conditions"]["total_volume"] = F(total[0], total[1], evidence)
    if procedure:
        r["procedure"]["original_text"] = F(procedure, procedure, evidence)
    if reactor:
        r["special_conditions"]["reactor"] = F(reactor, reactor, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved or [])}
    r["completeness"] = "partial"
    classify(r, changes=changes, unresolved=class_unresolved)
    reactions.append(r)
    return r


def pids(parts):
    return [p["participant_id"] for p in parts]


# Page 7: finish the page-6 alternative 7a -> 6 experiment.
e367 = [ev(3), ev(6), ev(7)]
parts = [
    participant("p_7a", "m_7a", "reactant", e367,
                charges={"mass": ("75 mg", q(75, "mg")), "amount": ("0.20 mmol", q(0.20, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_7a"),
    participant("p_5", "m_5", "reactant", e367,
                charges={"volume": ("37 uL", q(37, "uL")), "amount": ("0.22 mmol", q(0.22, "mmol")), "equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference="m_7a"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e367, charges={"volume": ("0.6 mL", q(0.6, "mL"))}),
    participant("p_6", "m_6", "product", e367, product_role="main product"),
    participant("p_alumina", "m_alumina", "purification_agent", e367, stage="purification"),
    participant("p_hexanes", "m_hexanes", "purification_agent", e367, stage="purification"),
    participant("p_ch2cl2_pur", "m_ch2cl2", "purification_agent", e367, stage="purification")]
reaction(
    f"{DOC}:scheme4_7a_5_to_6", e367, categories=["other"],
    location={"page": 3, "scheme": "4", "arrow": "middle right", "detailed_procedure_pages": [6, 7]},
    name="alternative synthesis of iron complex 6 from 7a",
    transform="formal [3 + 2] cycloaddition with isomerization", participants=parts,
    steps=[
        step("reaction", e367, ["p_7a", "p_5", "p_ch2cl2", "p_6"], "stirred", ("40 °C", q(40, "°C")), ("16 h", q(16, "h")), "argon/glovebox"),
        step("concentration", [ev(6)], ["p_6"], "concentrated in vacuo", stage="workup"),
        step("purification", [ev(6)], ["p_6", "p_alumina", "p_hexanes", "p_ch2cl2_pur"], "flash chromatography on neutral aluminum oxide (hexanes/CH2Cl2 1/0 to 10/1 to 0/1), then recrystallization from hexanes/CH2Cl2 (10/1 v/v)", stage="purification")],
    outcomes=[
        outcome("yield", "yield", ["m_6"], e367, "75%", q(75, "%"), "main product", "reported yield", ISO_BASIS, {"product": "6"}),
        outcome("amount", "product_amount", ["m_6"], [ev(7)], "86 mg", q(86, "mg"), "main product", basis=("isolated product", "isolated", [ev(7)]), scope={"product": "6"})],
    procedure="The page-6 procedure is continued on page 7: product 6 was obtained after chromatography and recrystallization as a bright yellow solid; yield 86 mg, 75%.",
    solvents=("CH2Cl2 (0.6 mL)", ["CH2Cl2"], ("0.6 mL", q(0.6, "mL"))),
    reactor="threaded reaction tube (5 mL) with magnetic stir bar",
    changes=[("C", "C", 0, 1, "formed"), ("C", "N", 0, 1, "formed")],
    unresolved=["Exact organometallic atom mapping and absolute configuration of the product ring center are unresolved."])


# Page 7: detailed 7a + 5 -> 6' experiment.
e37 = [ev(3), ev(7)]
parts = [
    participant("p_7a", "m_7a", "reactant", e37,
                charges={"mass": ("98 mg", q(98, "mg")), "amount": ("0.26 mmol", q(0.26, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_7a"),
    participant("p_5", "m_5", "reactant", e37,
                charges={"volume": ("45 uL", q(45, "uL")), "amount": ("0.29 mmol", q(0.29, "mmol")), "equivalents": ("1.1 equiv", q(1.1, "equiv"))}, reference="m_7a"),
    participant("p_dce", "m_dce", "solvent", e37, charges={"volume": ("2 mL", q(2, "mL"))}),
    participant("p_6p", "m_6p", "product", e37, product_role="main product"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(7)], stage="purification"),
    participant("p_ch2cl2", "m_ch2cl2", "purification_agent", [ev(7)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(7)], stage="purification")]
reaction(
    f"{DOC}:scheme4_7a_5_to_6prime", e37, categories=["other"],
    location={"page": 3, "scheme": "4", "arrow": "7a + 5 to 6'", "detailed_procedure_page": 7},
    name="synthesis of iron complex 6' from 7a", transform="formal [3 + 2] cycloaddition",
    participants=parts,
    steps=[
        step("reaction", e37, ["p_7a", "p_5", "p_dce", "p_6p"], "stirred at room temperature", ("room temperature", "r.t."), ("1 h", q(1, "h")), "argon/glovebox"),
        step("purification", [ev(7)], ["p_6p", "p_hexanes", "p_ch2cl2", "p_etoac"], "loaded directly onto a column equilibrated with hexanes/CH2Cl2 (1/1 v/v) and quickly chromatographed with hexanes/CH2Cl2/EtOAc 9/1/0 to 1/1/0 to 0/49/1", stage="purification")],
    outcomes=[
        outcome("yield", "yield", ["m_6p"], e37, "43%", q(43, "%"), "main product", "reported yield", ISO_BASIS, {"product": "6'"}),
        outcome("amount", "product_amount", ["m_6p"], [ev(7)], "64 mg", q(64, "mg"), "main product", basis=("isolated amorphous solid", "isolated", [ev(7)]), scope={"product": "6'"})],
    procedure="In an argon-filled glovebox, 7a (98 mg, 0.26 mmol) in dry DCE (2 mL) was treated with tosyl isocyanate (45 μL, 0.29 mmol), capped, removed from the glovebox, stirred at room temperature for 1 h, and rapidly chromatographed to afford 6' (64 mg, 43%).",
    solvents=("dry DCE (2 mL)", ["DCE"], ("2 mL", q(2, "mL"))),
    reactor="threaded reaction tube (5 mL) with magnetic stir bar and septum cap",
    changes=[("C", "C", 0, 1, "formed"), ("C", "N", 0, 1, "formed")],
    unresolved=["Exact organometallic atom mapping and absolute configuration remain unresolved."])


def derivatization_record(suffix, product_mid, parts, steps, outcomes, procedure,
                          solvents, changes, unresolved, notes=None):
    evidence = union([ev(3)], [ev(7)])
    return reaction(
        f"{DOC}:{suffix}", evidence, categories=["derivatization"],
        location={"page": 3, "scheme": "3", "arrow": suffix, "detailed_procedure_page": 7},
        name="divergent transformation of vinyliron complex 3", transform="demetalative functionalization",
        participants=parts, steps=steps, outcomes=outcomes, notes=notes,
        procedure=procedure, solvents=solvents, changes=changes, unresolved=unresolved)


# 3 -> 3a with NBS.
e = [ev(3), ev(7)]
parts = [
    participant("p_3", "m_3", "reactant", e, charges={"mass": ("13 mg", q(13, "mg")), "amount": ("23 umol", q(23, "umol")), "equivalents": ("1.0 equiv", q(1, "equiv"))}, reference="m_3"),
    participant("p_nbs", "m_nbs", "reagent", e, charges={"mass": ("4.2 mg", q(4.2, "mg")), "amount": ("23 umol", q(23, "umol")), "equivalents": ("1.0 equiv", q(1, "equiv"))}, reference="m_3"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e, charges={"volume": ("0.5 mL", q(0.5, "mL"))}),
    participant("p_3a", "m_3a", "product", e, product_role="main product"),
    participant("p_silica", "m_silica", "purification_agent", [ev(7)], stage="purification"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(7)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(7)], stage="purification")]
derivatization_record(
    "scheme3_3_to_3a", "m_3a", parts,
    [step("reaction", e, ["p_3", "p_nbs", "p_ch2cl2", "p_3a"], "stirred under air until disappearance of starting material", ("24 °C", q(24, "°C")), ("1 h", q(1, "h")), "air"),
     step("concentration", [ev(7)], ["p_3a"], "concentrated in vacuo", stage="workup"),
     step("purification", [ev(7)], ["p_3a", "p_silica", "p_hexanes", "p_etoac"], "preparative TLC on SiO2 with hexanes/EtOAc 85/15", stage="purification")],
    [outcome("yield", "yield", ["m_3a"], e, "54%", q(54, "%"), "main product", "reported yield", ISO_BASIS, {"product": "3a"}),
     outcome("amount", "product_amount", ["m_3a"], [ev(7)], "4.6 mg", q(4.6, "mg"), "main product", basis=("isolated product", "isolated", [ev(7)]), scope={"product": "3a"})],
    "Recently chromatographed complex 3 (13 mg, 23 μmol) in dry CH2Cl2 (0.5 mL) was treated with NBS (4.2 mg, 23 μmol), stirred at 24 °C under air for 1 h, concentrated, and purified by preparative TLC to give 3a (4.6 mg, 54%).",
    ("dry CH2Cl2 (0.5 mL)", ["CH2Cl2"], ("0.5 mL", q(0.5, "mL"))),
    [("C", "Br", 0, 1, "formed")],
    ["Exact organometallic atom mapping and the unassigned stereocenter in 3a remain unresolved."])


# Distinct Br2 alternative for 3a.
e7 = [ev(7)]
parts = [participant("p_3", "m_3", "reactant", e7),
         participant("p_br2", "m_br2", "reagent", e7,
                     charges={"equivalents": ("stoichiometric Br2", None)}, reference="m_3"),
         participant("p_3a", "m_3a", "product", e7, product_role="desired product")]
reaction(
    f"{DOC}:p7_3_to_3a_br2_alternative", e7, categories=["derivatization"],
    location={"page": 7, "region": "3a experimental paragraph, Br2 alternative"},
    name="alternative bromodemetalation of complex 3 with bromine", transform="bromodemetalation",
    participants=parts,
    steps=[step("reaction", e7, pids(parts), "Br2 was used in place of NBS; other detailed conditions were not restated")],
    outcomes=[outcome("yield", "yield", ["m_3a"], e7, "15%", q(15, "%"), "desired product", "reported yield", ISO_BASIS, {"product": "3a"}),
              outcome("amount", "product_amount", ["m_3a"], e7, "1.3 mg", q(1.3, "mg"), "desired product", basis=("isolated product", "isolated", e7), scope={"product": "3a"}),
              outcome("side_products", "other", ["m_3a"], e7, "side product formation was apparent", "side products apparent", scope={"observation": "side-product formation; identities and amounts unreported"})],
    notes=["The source says only 'stoichiometric Br2' and does not state an exact equivalent or repeat the NBS procedure conditions; those conditions are not silently inherited."],
    changes=[("C", "Br", 0, 1, "formed")],
    unresolved=["Exact Br2 charge, scale, solvent, temperature, time, workup, purification, side-product identities, and measurement method are not restated."])


# 3 -> 3b.
e = [ev(3), ev(7)]
parts = [
    participant("p_3", "m_3", "reactant", e, charges={"mass": ("29.2 mg", q(29.2, "mg")), "amount": ("54.9 umol", q(54.9, "umol"))}),
    participant("p_can", "m_can", "reagent", e, charges={"mass": ("57.4 mg", q(57.4, "mg")), "amount": ("105 umol", q(105, "umol")), "equivalents": ("1.9 equiv", q(1.9, "equiv"))}, reference="m_3"),
    participant("p_k2co3", "m_k2co3", "reagent", e, charges={"mass": ("22.2 mg", q(22.2, "mg")), "amount": ("161 umol", q(161, "umol")), "equivalents": ("2.9 equiv", q(2.9, "equiv"))}, reference="m_3"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e, charges={"volume": ("0.3 mL", q(0.3, "mL"))}),
    participant("p_etoh", "m_etoh", "solvent", e, charges={"volume": ("1.0 mL", q(1.0, "mL"))}),
    participant("p_3b", "m_3b", "product", e, product_role="main product"),
    participant("p_silica", "m_silica", "workup_agent", [ev(7)], stage="workup"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(7)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(7)], stage="purification")]
derivatization_record(
    "scheme3_3_to_3b", "m_3b", parts,
    [step("reaction", e, ["p_3", "p_can", "p_k2co3", "p_ch2cl2", "p_etoh", "p_3b"], "stirred in air until disappearance of starting material", ("24 °C", q(24, "°C")), ("2 h", q(2, "h")), "air"),
     step("filtration", [ev(7)], ["p_3b", "p_silica"], "filtered through a plug of silica and concentrated in vacuo", stage="workup"),
     step("purification", [ev(7)], ["p_3b", "p_hexanes", "p_etoac"], "preparative plate chromatography on SiO2 with hexanes/EtOAc 85/15", stage="purification")],
    [outcome("yield", "yield", ["m_3b"], e, "43%", q(43, "%"), "main product", "reported yield", ISO_BASIS, {"product": "3b"}),
     outcome("amount", "product_amount", ["m_3b"], [ev(7)], "8.5 mg", q(8.5, "mg"), "main product", basis=("isolated product", "isolated", [ev(7)]), scope={"product": "3b"})],
    "Complex 3 (29.2 mg, 54.9 μmol) in dry CH2Cl2 (0.3 mL) was treated with K2CO3, EtOH, and CAN, stirred in air at 24 °C for 2 h, filtered, concentrated, and purified by preparative plate chromatography to give 3b (8.5 mg, 43%).",
    ("dry CH2Cl2 (0.3 mL) and EtOH (1.0 mL)", ["CH2Cl2", "EtOH"], None),
    [("C", "C", 0, 1, "formed")],
    ["Exact organometallic atom mapping and the unassigned stereocenter in 3b remain unresolved."])


# 3 -> 3c with staged temperature program.
e = [ev(3), ev(7)]
parts = [
    participant("p_3", "m_3", "reactant", e, charges={"mass": ("77.8 mg", q(77.8, "mg")), "amount": ("0.15 mmol", q(0.15, "mmol")), "equivalents": ("1.0 equiv", q(1, "equiv"))}, reference="m_3"),
    participant("p_i2", "m_i2", "reagent", e, charges={"volume": ("1.5 mL of 0.1 M I2 in CH2Cl2", q(1.5, "mL")), "amount": ("0.15 mmol calculated from 0.1 M × 1.5 mL", q(0.15, "mmol"), "inference", "Calculated from the explicitly reported iodine solution concentration and volume."), "equivalents": ("1.0 equiv", q(1, "equiv"))}, reference="m_3"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e, charges={"volume": ("0.5 mL reaction solvent plus 1.5 mL iodine-solution carrier", None)}),
    participant("p_3c", "m_3c", "product", e, product_role="main product"),
    participant("p_fpstar_i", "m_fpstar_i", "product", e, product_role="recovered iron fragment"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(7)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(7)], stage="purification")]
derivatization_record(
    "scheme3_3_to_3c", "m_3c", parts,
    [step("iodine_addition", e, ["p_3", "p_i2", "p_ch2cl2"], "iodine solution was added", ("0 °C", q(0, "°C"))),
     step("warm_hold", e, ["p_3", "p_i2", "p_ch2cl2", "p_3c", "p_fpstar_i"], "stirred after warming; transition duration not reported", ("24 °C", q(24, "°C")), ("2 h", q(2, "h"))),
     step("concentration", [ev(7)], ["p_3c", "p_fpstar_i"], "concentrated in vacuo", stage="workup"),
     step("purification", [ev(7)], ["p_3c", "p_fpstar_i", "p_hexanes", "p_etoac"], "flash chromatography with hexanes/EtOAc 85/15", stage="purification")],
    [outcome("yield_3c", "yield", ["m_3c"], e, "42%", q(42, "%"), "main product", "reported yield", ISO_BASIS, {"product": "3c"}),
     outcome("amount_3c", "product_amount", ["m_3c"], [ev(7)], "25 mg", q(25, "mg"), "main product", basis=("isolated product", "isolated", [ev(7)]), scope={"product": "3c"}),
     outcome("fpstar_i_recovery", "yield", ["m_fpstar_i"], [ev(7)], "63% recovered Fp*I", q(63, "%"), "recovered iron fragment", "reported recovery-like yield", basis=("source wording: recovered Fp*I", "product/released iron-fragment basis; analytical method not stated", [ev(7)]), scope={"entity": "Fp*I", "source_wording": "recovered"}),
     outcome("fpstar_i_amount", "product_amount", ["m_fpstar_i"], [ev(7)], "34 mg", q(34, "mg"), "recovered iron fragment", basis=("recovered Fp*I", "isolated/recovered material", [ev(7)]), scope={"entity": "Fp*I"})],
    "Complex 3 in dry CH2Cl2 was treated at 0 °C with 0.1 M I2 in CH2Cl2, stirred at 24 °C for 2 h, concentrated, and chromatographed to give 3c and recovered Fp*I.",
    ("dry CH2Cl2 (0.5 mL) plus CH2Cl2 carrier for iodine solution (1.5 mL)", ["CH2Cl2"], None),
    [("C", "I", 0, 1, "formed")],
    ["The duration of the 0 °C addition and the warming transition are unreported; the two charged CH2Cl2 portions are not asserted as a measured final volume.", "The basis for the reported Fp*I recovery is not stated."])


# 3 -> 3d.
e = [ev(3), ev(7)]
parts = [
    participant("p_3", "m_3", "reactant", e, charges={"mass": ("37 mg", q(37, "mg")), "amount": ("69 umol", q(69, "umol"))}),
    participant("p_tfa", "m_tfa", "reagent", e, charges={"volume": ("10.7 uL", q(10.7, "uL")), "amount": ("0.14 mmol", q(0.14, "mmol")), "equivalents": ("2.0 equiv", q(2.0, "equiv"))}, reference="m_3"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e, charges={"volume": ("0.5 mL", q(0.5, "mL"))}),
    participant("p_water", "m_water", "additive", [ev(3)]),
    participant("p_3d", "m_3d", "product", e, product_role="main product"),
    participant("p_basic_alumina", "m_basic_alumina", "workup_agent", [ev(7)], stage="workup"),
    participant("p_ch2cl2_wash", "m_ch2cl2", "workup_agent", [ev(7)], stage="workup"),
    participant("p_silica", "m_silica", "purification_agent", [ev(7)], stage="purification"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(7)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(7)], stage="purification")]
derivatization_record(
    "scheme3_3_to_3d", "m_3d", parts,
    [step("reaction", e, ["p_3", "p_tfa", "p_ch2cl2", "p_water", "p_3d"], "stirred overnight under air", ("24 °C", q(24, "°C"), "inherited_general_condition"), ("12 h", q(12, "h")), "air"),
     step("filtration", [ev(7)], ["p_3d", "p_basic_alumina", "p_ch2cl2_wash"], "filtered through a short plug of basic aluminum oxide and washed with CH2Cl2", stage="workup"),
     step("concentration", [ev(7)], ["p_3d"], "filtrate concentrated in vacuo", stage="workup"),
     step("purification", [ev(7)], ["p_3d", "p_silica", "p_hexanes", "p_etoac"], "preparative plate chromatography on SiO2 with hexanes/EtOAc 85/15", stage="purification")],
    [outcome("yield", "yield", ["m_3d"], e, "98%", q(98, "%"), "main product", "reported yield", ISO_BASIS, {"product": "3d"}),
     outcome("amount", "product_amount", ["m_3d"], [ev(7)], "21 mg", q(21, "mg"), "main product", basis=("isolated product", "isolated", [ev(7)]), scope={"product": "3d"})],
    "Complex 3 (37 mg, 69 μmol) in CH2Cl2 (0.5 mL, stored in air on the benchtop) was treated with TFA (10.7 μL, 0.14 mmol) and stirred for 12 h under air; the mixture was filtered through basic alumina, concentrated, and purified by preparative plate chromatography to give 3d (21 mg, 98%).",
    ("Scheme 3: wet CH2Cl2; page 7: CH2Cl2 (0.5 mL, stored in air on the benchtop)", ["CH2Cl2", "water (wet solvent)"], ("0.5 mL", q(0.5, "mL"))),
    [("C", "H", 0, 1, "hydrogen_gain"), ("C", "O", 1, 2, "order_increased"), ("O", "H", 0, 1, "hydrogen_gain")],
    ["The scheme states 'wet CH2Cl2', while the detailed procedure describes CH2Cl2 stored in air on the benchtop; water amount and exact proton source are unreported."],
    notes=["Scheme 3 supplies 24 °C; the page-7 detailed paragraph states overnight (12 h) under air but does not repeat the temperature."])


# 6 -> 6a main 0.80-equivalent iodine experiment.
e = [ev(3), ev(7)]
eq_conflict_6 = charge_fact(
    ["1.3 equiv (Scheme 4)", "1.0 equiv (detailed procedure)"], None, e,
    status="conflict", alternatives=[
        F("1.3 equiv", q(1.3, "equiv"), [ev(3)]),
        F("1.0 equiv", q(1.0, "equiv"), [ev(7)])])
eq_conflict_i2 = charge_fact(
    ["1.0 equiv (Scheme 4)", "0.80 equiv (detailed procedure)"], None, e,
    status="conflict", alternatives=[
        F("1.0 equiv", q(1.0, "equiv"), [ev(3)]),
        F("0.80 equiv", q(0.80, "equiv"), [ev(7)])])
parts = [
    participant("p_6", "m_6", "reactant", e, charges={"mass": ("33 mg", q(33, "mg")), "amount": ("58 umol", q(58, "umol")), "equivalents": eq_conflict_6}, reference="m_6"),
    participant("p_i2", "m_i2", "reagent", e, charges={"volume": ("460 uL of 0.1 M I2 in CHCl3", q(460, "uL")), "amount": ("46 umol calculated from 0.1 M × 460 uL", q(46, "umol"), "inference", "Calculated from the source concentration and volume."), "equivalents": eq_conflict_i2}, reference="m_6"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e, charges={"volume": ("0.5 mL", q(0.5, "mL"))}),
    participant("p_chcl3", "m_chcl3", "solvent", e, charges={"volume": ("460 uL iodine-solution carrier", q(460, "uL"))}),
    participant("p_6a", "m_6a", "product", e, product_role="main product"),
    participant("p_fpstar_i", "m_fpstar_i", "product", e, product_role="recovered iron fragment"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(7)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(7)], stage="purification")]
reaction(
    f"{DOC}:scheme4_6_to_6a", e, categories=["derivatization"],
    location={"page": 3, "scheme": "4", "arrow": "bottom", "detailed_procedure_page": 7},
    name="iododemetalation of complex 6", transform="iododemetalation", participants=parts,
    steps=[
        step("iodine_addition", e, ["p_6", "p_i2", "p_ch2cl2", "p_chcl3"], "iodine solution added dropwise", ("0 °C", q(0, "°C")), ("15 min", q(15, "min"))),
        step("low_temperature_hold", e, ["p_6", "p_i2", "p_ch2cl2", "p_chcl3", "p_6a", "p_fpstar_i"], "stirred at 0 °C after addition", ("0 °C", q(0, "°C")), ("30 min", q(30, "min"))),
        step("room_temperature_hold", e, ["p_6", "p_i2", "p_ch2cl2", "p_chcl3", "p_6a", "p_fpstar_i"], "then stirred at room temperature", ("room temperature", "r.t."), ("another 2 h", q(2, "h"))),
        step("purification", [ev(7)], ["p_6a", "p_fpstar_i", "p_hexanes", "p_etoac"], "chromatographed directly: hexanes/EtOAc 95/5 to desired product, then 85/15 to 4/1 to Fp*I", stage="purification")],
    outcomes=[
        outcome("yield_6a", "yield", ["m_6a"], e, "59%", q(59, "%"), "main product", "reported yield", ISO_BASIS, {"product": "6a"}),
        outcome("amount_6a", "product_amount", ["m_6a"], [ev(7)], "12 mg", q(12, "mg"), "main product", basis=("isolated product", "isolated", [ev(7)]), scope={"product": "6a"}),
        outcome("fpstar_i_recovery", "yield", ["m_fpstar_i"], [ev(7)], "53% recovered Fp*I", q(53, "%"), "recovered iron fragment", "reported recovery-like yield", basis=("source wording: recovered Fp*I", "product/released iron-fragment basis; analytical method not stated", [ev(7)]), scope={"entity": "Fp*I", "source_wording": "recovered"}),
        outcome("fpstar_i_amount", "product_amount", ["m_fpstar_i"], [ev(7)], "9.2 mg", q(9.2, "mg"), "recovered iron fragment", basis=("recovered Fp*I", "isolated/recovered material", [ev(7)]), scope={"entity": "Fp*I"})],
    notes=["Scheme 4 reports 6 at 1.3 equiv and I2 at 1.0 equiv, while the detailed procedure reports 6 as 1.0 equiv and I2 as 0.80 equiv; both source-specific pairs are retained as participant-level conflicts without selecting one."],
    procedure="Complex 6 (33 mg, 58 μmol) in dry CH2Cl2 (0.5 mL) at 0 °C was treated dropwise over 15 min with 0.1 M iodine in CHCl3 (460 μL, 0.80 equiv), held at 0 °C for 30 min, warmed to room temperature for 2 h, and directly chromatographed to give 6a and recovered Fp*I.",
    solvents=("dry CH2Cl2 (0.5 mL) plus CHCl3 carrier (460 uL)", ["CH2Cl2", "CHCl3"], None),
    changes=[("C", "I", 0, 1, "formed")],
    unresolved=["The basis for the reported Fp*I recovery and exact organometallic atom mapping are not stated; C5 configuration is unassigned."])


# Distinct 1.0-equivalent comparison printed with an unresolved starting label.
e7 = [ev(7)]
parts = [participant("p_start_literal", "m_p7_alt6a_start_literal", "reactant", e7,
                     charges={"equivalents": ("1.0 equiv", q(1, "equiv"))}, reference="m_i2"),
         participant("p_i2", "m_i2", "reagent", e7,
                     charges={"equivalents": ("1.0 equiv", q(1, "equiv"))}, reference="m_p7_alt6a_start_literal"),
         participant("p_6a", "m_6a", "product", e7, product_role="desired product")]
reaction(
    f"{DOC}:p7_6a_i2_1equiv_comparison", e7, categories=["derivatization"],
    location={"page": 7, "region": "6a experimental paragraph, 1.0-equivalent comparison"},
    name="1.0-equivalent iodine comparison for formation of 6a", transform="iododemetalation comparison",
    participants=parts,
    steps=[step("reaction", e7, pids(parts), "the source states that 1.0 equiv each were used instead; other conditions are not repeated")],
    outcomes=[outcome("yield", "yield", ["m_6a"], e7, "31%", q(31, "%"), "desired product", "reported yield", ISO_BASIS, {"product": "6a"}),
              outcome("amount", "product_amount", ["m_6a"], e7, "7.2 mg", q(7.2, "mg"), "desired product", basis=("isolated product", "isolated", e7), scope={"product": "6a"})],
    notes=["The sentence literally reads '1.0 equiv each of 6a and I2'; the starting-material label is therefore retained as unresolved rather than silently changed to complex 6."],
    changes=[("C", "I", 0, 1, "formed")],
    unresolved=["Starting-material identity is internally problematic as printed; solvent, temperature program, time, scale, workup, purification, and recovery result are not repeated for this comparison."])


# General procedure for N-arylpyrazoles is a generic scheme/procedure, not a
# unique concrete experiment.
gp_id = f"{DOC}:p7_general_procedure_n_arylpyrazoles"
e7 = [ev(7)]
parts = [
    participant("p_diazonium", "m_generic_diazonium_pyrazole_gp", "reactant", e7,
                charges={"amount": ("0.30 mmol", q(0.30, "mmol")), "equivalents": ("1.5 equiv", q(1.5, "equiv"))}, reference="m_generic_allenyliron_pyrazole_gp"),
    participant("p_allenyliron", "m_generic_allenyliron_pyrazole_gp", "reactant", e7,
                charges={"amount": ("0.20 mmol", q(0.20, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_generic_allenyliron_pyrazole_gp"),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", e7, charges={"volume": ("0.5 mL", q(0.5, "mL"))}),
    participant("p_silica", "m_silica", "purification_agent", e7, stage="purification")]
reaction(
    gp_id, e7, record_type="generic_scheme", location={"page": 7, "region": "General Procedure for the Synthesis of N-Arylpyrazoles"},
    name="General Procedure for the Synthesis of N-Arylpyrazoles", transform="generic formal [3 + 2] cycloaddition",
    participants=parts,
    steps=[
        step("addition", e7, ["p_diazonium", "p_allenyliron", "p_ch2cl2"], "solution of freshly prepared allenyliron complex in CH2Cl2 was added in rapid succession to diazonium tetrafluoroborate", atmosphere="argon-filled glovebox"),
        step("heated_reaction", e7, ["p_diazonium", "p_allenyliron", "p_ch2cl2"], "sealed tube placed in an oil bath preheated to 60 °C", ("60 °C", q(60, "°C")), ("typically 12 h", q(12, "h"))),
        step("cooling", e7, [], "cooled to room temperature", ("room temperature", "r.t.")),
        step("concentration", e7, [], "crude mixture concentrated in vacuo", stage="workup"),
        step("purification", e7, ["p_silica"], "purified by flash column chromatography", stage="purification")],
    procedure="In an argon-filled glovebox, diazonium tetrafluoroborate (0.30 mmol, 1.5 equiv) was charged to a sealed tube. Freshly prepared allenyliron complex (0.20 mmol, 1.0 equiv) in CH2Cl2 (0.5 mL) was added rapidly; the sealed tube was heated at 60 °C, typically 12 h, cooled, concentrated, and chromatographed.",
    solvents=("CH2Cl2 (0.5 mL)", ["CH2Cl2"], ("0.5 mL", q(0.5, "mL"))),
    reactor="sealed reaction tube with magnetic stir bar",
    class_unresolved=True,
    unresolved=["The general procedure does not define one concrete allenyliron/diazonium pair or product; concrete entries are recorded separately."])


def pyrazole_entry(suffix, page, iron_mid, iron_label, iron_mass, iron_amount,
                   diazo_mid, diazo_label, diazo_mass, diazo_amount,
                   product_mid, product_label, product_mass, yield_value,
                   chromatography, rr, rr_method, notes=None):
    evidence = [ev(4), ev(page)]
    parts = [
        participant("p_iron", iron_mid, "reactant", evidence,
                    charges={"mass": (iron_mass, q(float(iron_mass.split()[0]), "mg")), "amount": (iron_amount, q(float(iron_amount.split()[0]), "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference=iron_mid),
        participant("p_diazo", diazo_mid, "reactant", evidence,
                    charges={"mass": (diazo_mass, q(float(diazo_mass.split()[0]), "mg")), "amount": (diazo_amount, q(float(diazo_amount.split()[0]), "mmol")), "equivalents": ("1.5 equiv", q(1.5, "equiv"))}, reference=iron_mid),
        participant("p_ch2cl2", "m_ch2cl2", "solvent", evidence, charges={"volume": ("0.5 mL from general procedure", q(0.5, "mL"), "inherited_general_condition", "The entry explicitly follows the page-7 general procedure.")}),
        participant("p_product", product_mid, "product", evidence, product_role="major/main regioisomer"),
        participant("p_hexanes", "m_hexanes", "purification_agent", [ev(page)], stage="purification"),
        participant("p_etoac", "m_etoac", "purification_agent", [ev(page)], stage="purification")]
    selectivity_targets = [product_mid]
    return reaction(
        f"{DOC}:{suffix}", evidence, categories=["substrate_scope"],
        location={"page": 4, "scheme": "5", "entry": product_label, "detailed_procedure_page": page},
        name="pyrazole synthesis using cyclopentadienyl-allenyliron cycloaddition to aryldiazonium salt",
        transform="formal [3 + 2] cycloaddition", participants=parts,
        steps=[
            step("reaction", evidence, ["p_iron", "p_diazo", "p_ch2cl2", "p_product"], "heated following the general procedure", ("60 °C", q(60, "°C")), ("12 h", q(12, "h")), "argon/glovebox then sealed tube"),
            step("concentration", [ev(page)], ["p_product"], "crude residue concentrated in vacuo", stage="workup"),
            step("purification", [ev(page)], ["p_product", "p_hexanes", "p_etoac"], f"flash chromatography ({chromatography})", stage="purification")],
        outcomes=[
            outcome("yield", "yield", [product_mid], evidence, f"{yield_value}%", q(yield_value, "%"), "main isolated product", "reported yield", ISO_BASIS, {"product": product_label}),
            outcome("amount", "product_amount", [product_mid], [ev(page)], product_mass, q(float(product_mass.split()[0]), "mg"), "main isolated product", basis=("isolated product", "isolated", [ev(page)]), scope={"product": product_label}),
            outcome("rr", "selectivity", selectivity_targets, [ev(page)], rr, rr, "regioisomeric ratio", basis=(rr_method, rr_method, [ev(page)]), scope={"metric": "rr", "order": [product_label, "unobserved/unregistered regioisomer"], "sample": "crude residue"})],
        notes=notes,
        procedure=f"Following the page-7 general procedure, {iron_label} and {diazo_label} were heated at 60 °C for 12 h; the crude residue was chromatographed ({chromatography}) to give {product_label}.",
        general_scheme_ids=[gp_id], solvents=("CH2Cl2 (0.5 mL), inherited from the general procedure", ["CH2Cl2"], ("0.5 mL", q(0.5, "mL"))),
        changes=[("C", "N", 0, 1, "formed"), ("C", "N", 0, 1, "formed")],
        unresolved=["Exact organometallic atom mapping and the identity/detection limit of the unobserved minor regioisomer remain unresolved."])


pyrazole_entry("scheme5_s5_10a", 7, "m_7b", "7b", "59 mg", "0.19 mmol", "m_9a", "phenyldiazonium tetrafluoroborate", "56 mg", "0.29 mmol", "m_10a", "10a", "12 mg", 27, "hexanes/EtOAc 30/1 to 10/1", ">100:1 rr", "1H NMR of the crude residue")
pyrazole_entry("scheme5_s5_10b", 7, "m_7b", "7b", "53 mg", "0.17 mmol", "m_9b", "4-methoxybenzenediazonium tetrafluoroborate", "56 mg", "0.25 mmol", "m_10b", "10b", "16 mg", 36, "hexanes/EtOAc 30/1 to 5/1", ">100:1 rr", "1H NMR of the crude residue")
pyrazole_entry("scheme5_s5_10c", 8, "m_7b", "7b", "59 mg", "0.19 mmol", "m_9c", "4-nitrobenzenediazonium tetrafluoroborate", "70 mg", "0.29 mmol", "m_10c", "10c", "18 mg", 35, "hexanes/EtOAc 30/1 to 10/1", ">100:1 rr", "1H NMR of the crude residue")
pyrazole_entry("scheme5_s5_10d", 8, "m_7b", "7b", "53 mg", "0.17 mmol", "m_9d", "4-(methoxycarbonyl)benzenediazonium tetrafluoroborate", "64 mg", "0.25 mmol", "m_10d", "10d", "17 mg", 35, "hexanes/EtOAc 30/1 to 10/1", ">100:1 rr", "1H NMR of the crude residue")
pyrazole_entry("scheme5_11_9c_to_12", 8, "m_11", "11", "49 mg", "0.19 mmol", "m_9c", "4-nitrobenzenediazonium tetrafluoroborate", "70 mg", "0.29 mmol", "m_12", "12", "17 mg", 40, "hexanes/EtOAc 30/1 to 10/1", "9.0:1 rr", "GCMS", notes=["The connectivity of 12 was assigned by nuclear Overhauser effect spectroscopy; this is characterization context, not a separate synthesis."])


# Page 8: 7b + DIAD -> 14 with preserved 1.5 vs 1.6 equiv conflict.
e48 = [ev(4), ev(8)]
eq_conflict_13 = charge_fact(
    ["1.5 equiv (Scheme 6)", "1.6 equiv (detailed procedure)"], None, e48,
    status="conflict", alternatives=[F("1.5 equiv", q(1.5, "equiv"), [ev(4)]), F("1.6 equiv", q(1.6, "equiv"), [ev(8)])])
eq_conflict_bf3 = copy.deepcopy(eq_conflict_13)
parts = [
    participant("p_13", "m_13", "reactant", e48,
                charges={"volume": ("60 uL", q(60, "uL")), "amount": ("0.31 mmol", q(0.31, "mmol")), "equivalents": eq_conflict_13}, reference="m_7b"),
    participant("p_bf3", "m_bf3", "additive", e48,
                charges={"volume": ("38 uL", q(38, "uL")), "amount": ("0.31 mmol", q(0.31, "mmol")), "equivalents": eq_conflict_bf3}, reference="m_7b"),
    participant("p_ch2cl2_initial", "m_ch2cl2", "solvent", [ev(8)], charges={"volume": ("0.25 mL", q(0.25, "mL"))}),
    participant("p_7b", "m_7b", "reactant", e48,
                charges={"mass": ("59 mg", q(59, "mg")), "amount": ("0.19 mmol", q(0.19, "mmol")), "equivalents": ("1.0 equiv", q(1.0, "equiv"))}, reference="m_7b"),
    participant("p_ch2cl2_solution", "m_ch2cl2", "solvent", [ev(8)], charges={"volume": ("0.25 mL", q(0.25, "mL"))}),
    participant("p_14", "m_14", "product", e48, product_role="main product"),
    participant("p_hexanes", "m_hexanes", "purification_agent", [ev(8)], stage="purification"),
    participant("p_etoac", "m_etoac", "purification_agent", [ev(8)], stage="purification")]
reaction(
    f"{DOC}:scheme6_7b_13_to_14", e48, categories=["substrate_scope"],
    location={"page": 4, "scheme": "6", "detailed_procedure_page": 8},
    name="synthesis of isopropyl 5-methyl-3-phenyl-1H-pyrazole-1-carboxylate (14)",
    transform="formal [3 + 2] cycloaddition with dealkoxycarbonylation", participants=parts,
    steps=[
        step("premix", [ev(8)], ["p_13", "p_bf3", "p_ch2cl2_initial"], "DIAD and BF3·OEt2 were dissolved in CH2Cl2 in a glovebox", atmosphere="glovebox"),
        step("addition", [ev(8)], ["p_7b", "p_ch2cl2_solution", "p_13", "p_bf3", "p_ch2cl2_initial"], "solution of 7b in CH2Cl2 was added in rapid succession", atmosphere="glovebox"),
        step("reaction", e48, ["p_13", "p_bf3", "p_7b", "p_ch2cl2_initial", "p_ch2cl2_solution", "p_14"], "stirred", ("60 °C", q(60, "°C")), ("12 h", q(12, "h"))),
        step("purification", [ev(8)], ["p_14", "p_hexanes", "p_etoac"], "flash chromatography with hexanes/EtOAc 30/1 to 5/1", stage="purification")],
    outcomes=[
        outcome("yield", "yield", ["m_14"], e48, "36%", q(36, "%"), "main product", "reported yield", ISO_BASIS, {"product": "14"}),
        outcome("amount", "product_amount", ["m_14"], [ev(8)], "17 mg", q(17, "mg"), "main product", basis=("isolated product", "isolated", [ev(8)]), scope={"product": "14"}),
        outcome("rr", "selectivity", ["m_14"], [ev(8)], ">20:1 rr", ">20:1", "regioisomeric ratio", basis=("1H NMR of the crude residue", "1H NMR of crude residue", [ev(8)]), scope={"metric": "rr", "order": ["14", "unobserved/unregistered regioisomer"], "sample": "crude residue"})],
    notes=["Scheme 6 reports 13 and BF3·Et2O at 1.5 equiv each; the detailed page-8 procedure reports 1.6 equiv each. Both values are retained as a conflict without selecting one.", "The two separate 0.25 mL CH2Cl2 charges are retained as participant-level volumes; no measured final total volume is asserted."],
    procedure="In a glovebox, DIAD (60 μL, 0.31 mmol, printed 1.6 equiv) and BF3·OEt2 (38 μL, 0.31 mmol, printed 1.6 equiv) were dissolved in CH2Cl2 (0.25 mL). A solution of 7b (59 mg, 0.19 mmol) in CH2Cl2 (0.25 mL) was added rapidly; the mixture was stirred at 60 °C for 12 h and chromatographed to give 14 (17 mg, 36%).",
    solvents=("two separate CH2Cl2 portions, 0.25 mL each", ["CH2Cl2"], None),
    reactor="reaction tube with magnetic stir bar; glovebox setup",
    changes=[("C", "N", 0, 1, "formed"), ("C", "N", 0, 1, "formed")],
    unresolved=["Scheme/procedure equivalent values conflict (1.5 vs 1.6); exact organometallic atom mapping and identity/detection limit of the minor regioisomer remain unresolved."])


# Page 10 footnote 36: a current experiment, plus a distinct proposed pathway.
e10 = [ev(10), ev(10, "image")]
parts = [
    participant("p_1a", "m_1a", "reactant", e10),
    participant("p_imine", "m_benzaldehyde_tosylimine", "reactant", e10),
    participant("p_bf3", "m_bf3", "additive", e10),
    participant("p_dce", "m_dce", "solvent", e10),
    participant("p_product", "m_eq1_dihydropyrrole", "product", e10, product_role="unexpected product")]
eq1_id = f"{DOC}:p10_footnote36_eq1"
reaction(
    eq1_id, e10, categories=["mechanistic_study"],
    location={"page": 10, "footnote": 36, "equation": 1},
    name="reaction of 1a with benzaldehyde tosylimine", transform="formal [3 + 2] cycloaddition",
    participants=parts,
    steps=[step("reaction", e10, pids(parts), "heated with BF3·Et2O in DCE", ("80 °C", q(80, "°C")))],
    outcomes=[outcome("yield", "yield", ["m_eq1_dihydropyrrole"], e10, "15% yield", q(15, "%"), "unexpected product", "reported yield", ISO_BASIS, {"product": "N-tosyl-4,5-diphenyl-1H-2,5-dihydropyrrole"})],
    notes=["The experiment is reported in footnote 36 and Equation 1, not in the reference list. The detailed reagent amounts, time, workup, purification, and yield measurement method are not provided on page 10."],
    solvents=("DCE", ["DCE"], None),
    changes=[("C", "C", 0, 1, "formed"), ("C", "N", 0, 1, "formed")],
    unresolved=["Reactant/additive quantities, reaction time, atmosphere, workup, purification, and exact atom mapping are unreported; imine geometry and product C5 configuration are unassigned."])

reaction(
    f"{DOC}:p10_footnote36_proposed_propargyl_pathway", e10,
    record_type="proposed_mechanism", location={"page": 10, "footnote": 36, "equation": 1, "region": "via bracketed intermediate"},
    name="proposed propargyliron-mediated pathway for Equation 1",
    transform="proposed [3 + 2] cycloaddition through the propargyliron tautomer",
    participants=[participant("p_1a", "m_1a", "reactant", e10),
                  participant("p_imine", "m_benzaldehyde_tosylimine", "reactant", e10),
                  participant("p_product", "m_eq1_dihydropyrrole", "product", e10, product_role="depicted final product")],
    notes=["The bracketed 'via' species and statement that the result can be rationalized through propargyliron are a mechanistic proposal, not a separately observed preparation."],
    class_unresolved=True,
    unresolved=["The proposed intermediate is not registered as an observed molecule and is not treated as experimental proof."])


# Characterization and nonreaction page context.
source_context = [
    {"context_id": f"{DOC}:p7_p8_characterization", "document_id": DOC,
     "evidence_ids": [ev(7), ev(8)], "kind": "characterization",
     "description": "NMR/IR/HRMS/elemental-analysis data following the preparative paragraphs are characterization, not additional synthesis experiments."},
    {"context_id": f"{DOC}:p8_noe_12", "document_id": DOC,
     "evidence_ids": [ev(8)], "kind": "connectivity_assignment",
     "description": "The connectivity of product 12 was assigned by NOE spectroscopy; no separate reaction record is created."},
    {"context_id": f"{DOC}:p8_independent_synthesis_14", "document_id": DOC,
     "evidence_ids": [ev(8)], "kind": "literature_characterization_support",
     "description": "Product 14 connectivity was verified by an independent literature synthesis described only in unavailable SI; no experiment is invented from that statement."},
    {"context_id": f"{DOC}:p8_associated_content", "document_id": DOC,
     "evidence_ids": [ev(8)], "kind": "supplementary_information_notice",
     "description": "The article states that SI exists online, but no SI file was explicitly supplied or associated in this workspace; no SI evidence is used."},
    {"context_id": f"{DOC}:p9_references", "document_id": DOC,
     "evidence_ids": [ev(9), ev(9, "image")], "kind": "references_only",
     "description": "Page 9 contains bibliographic references only and no current experimental arrow/table."},
    {"context_id": f"{DOC}:p10_references", "document_id": DOC,
     "evidence_ids": [ev(10), ev(10, "image")], "kind": "references_with_current_footnote_experiment",
     "description": "Most of page 10 is references; footnote 36/Equation 1 is separately recorded as a current mechanistic-study experiment and proposed pathway."}
]


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "left column top: continuation of alternative synthesis of 6", "finding": "Continuation supplies 86 mg, 75% and characterization; merged into existing Scheme 4 record.", "evidence_ids": [ev(7), ev(7, "image")]},
            {"region": "left column: synthesis of 6'", "finding": "Detailed 7a/TsNCO experiment, workup, 64 mg and 43% yield merged into the existing Scheme 4 6' record.", "evidence_ids": [ev(7), ev(7, "image")]},
            {"region": "left/right columns: 3a-3d preparations", "finding": "Four detailed preparations merged into Scheme 3 records; distinct Br2 alternative for 3a recorded separately. Spectral data retained as characterization context.", "evidence_ids": [ev(7), ev(7, "image")]},
            {"region": "right column: 6a preparation", "finding": "Main 0.80-equiv I2 experiment merged into Scheme 4; the lower-yield 1.0-equiv comparison is separate with the printed starting-label ambiguity preserved.", "evidence_ids": [ev(7), ev(7, "image")]},
            {"region": "right column bottom: N-arylpyrazole general procedure, 10a and 10b", "finding": "General procedure recorded as generic; concrete 10a/10b entries updated with all charges, chromatography, product masses, yields and crude-1H-NMR rr.", "evidence_ids": [ev(7), ev(7, "image")]}
        ], "record_ids": [
            f"{DOC}:scheme4_7a_5_to_6", f"{DOC}:scheme4_7a_5_to_6prime",
            f"{DOC}:scheme3_3_to_3a", f"{DOC}:p7_3_to_3a_br2_alternative",
            f"{DOC}:scheme3_3_to_3b", f"{DOC}:scheme3_3_to_3c", f"{DOC}:scheme3_3_to_3d",
            f"{DOC}:scheme4_6_to_6a", f"{DOC}:p7_6a_i2_1equiv_comparison",
            gp_id, f"{DOC}:scheme5_s5_10a", f"{DOC}:scheme5_s5_10b"],
         "unresolved": ["Dense organometallic structures are not completely atom-mapped.", "The printed starting label in the 1.0-equiv 6a comparison is internally problematic.", "Independent spectral peak lists are not modeled as reactions."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "left column: 10c, 10d and 12", "finding": "Three concrete pyrazole experiments updated with charges, 60 °C/12 h, chromatography, product amounts/yields and rr methods; NOE for 12 retained as characterization context.", "evidence_ids": [ev(8), ev(8, "image")]},
            {"region": "left column: synthesis of 14", "finding": "Detailed DIAD/BF3·Et2O experiment merged into Scheme 6; 1.5-equiv scheme versus 1.6-equiv prose conflict preserved.", "evidence_ids": [ev(8), ev(8, "image")]},
            {"region": "Associated Content/metadata and right column", "finding": "SI availability notice, author information, acknowledgments and start of references; no additional current synthesis. SI was not associated because no SI file was supplied.", "evidence_ids": [ev(8), ev(8, "image")]}
        ], "record_ids": [f"{DOC}:scheme5_s5_10c", f"{DOC}:scheme5_s5_10d", f"{DOC}:scheme5_11_9c_to_12", f"{DOC}:scheme6_7b_13_to_14"],
         "unresolved": ["Minor regioisomer identities/detection limits are not reported.", "Product 14's independent synthesis is only referenced to unavailable SI.", "Scheme/prose equivalent values for 14 conflict."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "both columns", "finding": "Bibliographic references 3-18 only; no current reaction arrow, table, experimental paragraph or characterization dataset.", "evidence_ids": [ev(9), ev(9, "image")]}
        ], "record_ids": [], "unresolved": []},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "left column and upper right column", "finding": "Bibliographic references 18-35 only; cited prior work is not converted into current experiments.", "evidence_ids": [ev(10), ev(10, "image")]},
            {"region": "right column footnote 36 and Equation 1", "finding": "Current 1a/benzaldehyde-tosylimine experiment recorded at 80 °C in DCE with 15% product yield; bracketed propargyliron route separately classified as proposed mechanism.", "evidence_ids": [ev(10), ev(10, "image")]}
        ], "record_ids": [eq1_id, f"{DOC}:p10_footnote36_proposed_propargyl_pathway"],
         "unresolved": ["Equation 1 omits reagent amounts, reaction time, atmosphere, workup, purification and yield measurement method.", "Imine geometry and product stereocenter configuration are not assigned."]}
    ]
}

visual = {
    "document_id": DOC,
    "inspected_original_pdf_pages": [7, 8, 9, 10],
    "chemical_accuracy_verified": False,
    "page_inventory": {
        "7": {"detailed_preparations": ["6 continuation", "6'", "3a", "3b", "3c", "3d", "6a", "10a", "10b"], "distinct_comparisons": ["3a by Br2", "6a 1.0-equiv comparison"], "generic_procedures": 1},
        "8": {"detailed_preparations": ["10c", "10d", "12", "14"], "nonreaction_characterization_notes": ["NOE connectivity for 12", "independent synthesis of 14 referenced to SI"]},
        "9": {"references_only": True},
        "10": {"references": True, "current_equations": ["footnote 36, eq 1"], "proposed_pathways": 1}
    },
    "notes": ["Counts come from direct inspection of the four original page images.", "Repeated Scheme experiments were merged by reaction_id rather than duplicated."]
}

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = source_context

(RUN / "submitted_pages_0007_0010.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics_pages_0007_0010.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/visual_transcription_pages_0007_0010.json").write_text(json.dumps(visual, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0007_0010.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules_in_bundle": len(bundle["molecules"]), "incoming_reactions": len(reactions),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
                  "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in reactions),
                  "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in reactions)}, ensure_ascii=False))
