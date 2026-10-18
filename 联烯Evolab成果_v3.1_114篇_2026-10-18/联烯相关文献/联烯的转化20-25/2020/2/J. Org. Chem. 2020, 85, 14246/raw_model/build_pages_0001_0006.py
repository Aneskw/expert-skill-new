#!/usr/bin/env python3
"""Build source-bounded v3.1 records for JOC 2020, 85, 14246, PDF pages 1–6.

This script contains the machine extraction decisions for the bounded increment.
It does not certify structures or chemistry and does not alter the extraction tools.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_7ca2cfc676d4ba33"


def ev(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


def fact(raw=None, normalized=None, evidence=(), method="original_extraction",
         status=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def inferred(value, evidence, derivation, raw=None):
    return fact(raw, value, evidence, method="inference", derivation=derivation)


def derived(value, evidence, derivation, raw=None):
    return fact(raw, value, evidence, method="structure_derivation", derivation=derivation)


def quantity(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw_text:
        value["raw_text"] = raw_text
    return value


def range_quantity(minimum, maximum, unit, raw_text=None):
    value = {"minimum": minimum, "maximum": maximum, "unit": unit}
    if raw_text:
        value["raw_text"] = raw_text
    return value


def qfact(raw, number, unit, evidence, qualifier=None):
    return fact(raw, quantity(number, unit, qualifier=qualifier, raw_text=raw), evidence)


def unresolved_molecule(mid, label, name, evidence, *, chirality="unresolved",
                        geometry="unresolved", representation="unresolved",
                        assignment=None, aliases=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "unresolved"
    m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    m["names"]["original"] = fact(name, name, evidence)
    m["names"]["english"] = fact(name, name, evidence)
    if aliases:
        m["names"]["aliases"] = fact(aliases, aliases, evidence)
    m["cas_numbers"] = fact(evidence=evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = fact(evidence=evidence)
    m["generic_template"] = None
    if assignment is not None:
        m["stereochemistry"]["assignment"] = fact(assignment, assignment, evidence)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence, "Explicit status retained because no complete, source-verified structure is asserted."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "Explicit status retained because no complete, source-verified structure is asserted."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence, "No complete stereochemical representation is asserted."
    )
    m["identity_note"] = inferred(
        "Source-labelled material; structure not promoted to a concrete SMILES in this machine extraction.",
        evidence,
        "Conservative registry decision for a reagent, mixture, catalyst, or otherwise unresolved material.",
    )
    return m


def structured_molecule(mid, label, name, smiles, evidence, *, atom_centers=None,
                        axial=None, assignment=None, chirality="specified",
                        geometry="not_applicable", representation="encoded",
                        aliases=None, identity_note=None, allene_feature=None):
    diag = st.stereo_inspect(smiles)
    canonical = diag["canonical_isomeric_smiles"]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, evidence)
    m["names"]["original"] = fact(name, name, evidence)
    m["names"]["english"] = fact(name, name, evidence)
    if aliases:
        m["names"]["aliases"] = fact(aliases, aliases, evidence)
    m["cas_numbers"] = fact(evidence=evidence)
    m["representation_type"] = "smiles"
    source_raw = {"source": "compound name and/or structure drawing", "compound_label": label}
    m["smiles"] = derived(
        canonical, evidence,
        "Candidate connectivity and encoded tetrahedral stereochemistry transcribed from the source name/drawing; RDKit canonicalization is not source verification.",
        raw=source_raw,
    )
    m["generic_template"] = None
    if atom_centers is not None:
        m["stereochemistry"]["atom_centers"] = fact(atom_centers, atom_centers, evidence)
    else:
        m["stereochemistry"]["atom_centers"] = derived(
            diag["atom_centers"], evidence,
            "RDKit diagnostic on the candidate isomeric SMILES; zero-based atom indices are toolkit indices.",
            raw=source_raw,
        )
    if axial is not None:
        m["stereochemistry"]["axial"] = fact(axial, axial, evidence)
    else:
        m["stereochemistry"]["axial"] = derived(
            [], evidence, "No axial assignment was encoded in the candidate graph.", raw=source_raw
        )
    m["stereochemistry"]["double_bonds"] = derived(
        diag["double_bonds"], evidence,
        "RDKit double-bond diagnostic on the candidate graph; cyclic and cumulated-bond stereochemistry still requires source review.",
        raw=source_raw,
    )
    m["stereochemistry"]["mixture"] = fact(evidence=evidence)
    if assignment is not None:
        m["stereochemistry"]["assignment"] = fact(assignment, assignment, evidence)
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, evidence, "Status follows explicit source stereochemical naming/drawing and the retained candidate graph."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "No independent acyclic E/Z assignment is applicable or asserted for this candidate."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence,
        "Tetrahedral stereochemistry is encoded; any named allene axis is not encoded when status is partial."
    )
    m["stereochemistry"]["isomeric_smiles"] = derived(
        canonical, evidence,
        "Canonical isomeric SMILES generated from the candidate graph; source_verified=false.",
        raw=source_raw,
    )
    m["stereochemistry"]["identity_key"] = derived(
        diag["stereo_graph_key"], evidence,
        "SHA-256 stereo graph key generated by the standard stereo diagnostic; not a sample identity.",
        raw=canonical,
    )
    m["identity_note"] = inferred(
        identity_note or "Machine candidate based on the source label/drawing; parsing and CIP diagnostics do not establish transcription correctness.",
        evidence,
        "Scientific review remains required.",
    )
    if allene_feature:
        m["allene_features"] = [{
            "fragment_id": allene_feature["fragment_id"],
            "carbon_atom_maps": derived(
                allene_feature.get("carbon_atom_maps"), evidence,
                "Allene positions described from the named chain; no full atom mapping is claimed.",
                raw=allene_feature.get("raw"),
            ),
            "substitution": fact(allene_feature.get("substitution"), allene_feature.get("substitution"), evidence),
            "axial_stereochemistry": fact(allene_feature.get("axial"), allene_feature.get("axial"), evidence),
            "local_environment": fact(allene_feature.get("environment"), allene_feature.get("environment"), evidence),
        }]
    return m


def conflict_threonine():
    evidence = [ev(2, "image"), ev(4, "image"), ev(2), ev(4)]
    m = unresolved_molecule(
        "m_threonine_source_conflict", "D-/L-threonine (source conflict)",
        "threonine starting material as printed", evidence,
        chirality="unresolved", geometry="not_applicable", representation="not_encoded",
    )
    a = fact("D-threonine", "D-threonine", [ev(2, "image"), ev(2)])
    b = fact("L-threonine", "L-threonine", [ev(4, "image"), ev(4)])
    m["names"]["original"] = fact(
        ["D-threonine", "L-threonine"], None, evidence,
        status="conflict", alternatives=[a, b],
    )
    m["names"]["english"] = copy.deepcopy(m["names"]["original"])
    m["identity_note"] = inferred(
        "Scheme 2 labels D-threonine, while the page-4 experimental procedure states L-threonine; no single enantiomer identity is chosen.",
        evidence,
        "Both source statements are retained pending expert review.",
    )
    return m


def participant(pid, mid, role, evidence, *, stage="reaction", mass=None, amount=None,
                volume=None, equivalents=None, mol_percent=None, reference=None,
                product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = fact(product_role, product_role, evidence)
    values = {
        "mass": mass,
        "amount": amount,
        "volume": volume,
        "equivalents": equivalents,
        "mol_percent": mol_percent,
    }
    for key, item in values.items():
        if item is not None:
            raw, normalized = item
            p["charges"][key] = fact(raw, normalized, evidence)
    if reference:
        p["equivalence_reference"] = fact(reference, reference, evidence)
    return p


def step(sid, evidence, *, stage="reaction", pids=(), operation=None,
         temperature=None, time=None, atmosphere=None, pressure=None):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    if operation:
        s["operation"] = fact(operation, operation, evidence)
    if temperature:
        raw, normalized = temperature
        s["temperature"] = fact(raw, normalized, evidence)
    if time:
        raw, normalized = time
        s["time"] = fact(raw, normalized, evidence)
    if atmosphere:
        s["atmosphere"] = fact(atmosphere, atmosphere, evidence)
    if pressure:
        raw, normalized = pressure
        s["pressure"] = fact(raw, normalized, evidence)
    return s


def outcome(oid, kind, targets, evidence, *, raw, normalized, yield_type=None,
            measurement=None, scope=None, product_role=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = fact(raw, normalized, evidence)
    if yield_type:
        o["yield_type"] = fact(yield_type, yield_type, evidence)
    if measurement:
        o["measurement_basis"] = fact(measurement, measurement, evidence)
    if scope:
        o["scope"] = fact(scope, scope, evidence)
    if product_role:
        o["product_role"] = fact(product_role, product_role, evidence)
    return o


def product_amount(oid, target, evidence, mass_raw, mass, amount_raw, amount):
    return outcome(
        oid, "product_amount", [target], evidence,
        raw={"mass": mass_raw, "amount": amount_raw},
        normalized={"mass": mass, "amount": amount},
        measurement="isolated material after the stated purification",
        scope="isolated product amount",
        product_role="main product",
    )


def classification(changes, evidence, coverage="partial", limitations=()):
    bond_facts = []
    labels = {}
    for i, item in enumerate(changes, 1):
        change = dict(item)
        change.setdefault("change_id", f"b{i}")
        bond_facts.append(derived(
            change, evidence,
            "Partial net bond inventory inferred from source-labelled structures; no complete atom mapping is asserted.",
            raw={"source": "reactant/product drawings and names"},
        ))
        if change["change_type"] in ("formed", "order_increased", "hydrogen_gain"):
            label = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"], change["after_order"])
            labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(change["change_id"])
    if coverage == "not_applicable":
        label_value = []
    elif coverage == "unresolved" and not changes:
        label_value = None
    else:
        label_value = list(labels.values())
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": fact(evidence=evidence),
        "bond_changes": bond_facts,
        "labels": inferred(
            label_value, evidence,
            "Taxonomy labels computed from the partial net bond inventory; no model-vote ensemble or source verification is claimed.",
        ) if label_value is not None else fact(evidence=evidence),
        "coverage": inferred(
            coverage, evidence,
            "Coverage is partial or unresolved because a complete source-checked atom mapping and coproduct inventory are unavailable."
        ),
        "limitations": list(limitations) or [
            "No complete atom mapping is asserted.",
            "Implicit-H changes are net changes; elementary proton sources are not assigned.",
            "Machine bond classification requires expert review.",
        ],
        "model_votes": [],
    }


def chg(a, b, before, after, kind):
    return {"atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after, "change_type": kind}


def reaction(rid, evidence, *, record_type="specific_experiment", categories=(),
             location=None, reaction_name=None, transformation=None, participants=(),
             steps=(), outcomes=(), class_changes=(), class_coverage="partial",
             procedure=None, addition_order=None, pretreatment=None, notes=(),
             unresolved=(), solvent=None, solvent_ratio=None, total_volume=None,
             concentration=None, concentration_basis=None, general_scheme_ids=()):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = inferred(
        list(categories), evidence,
        "Experiment-category vocabulary assigned from the source presentation and record identity."
    )
    r["location"] = fact(location or rid, location or rid, evidence)
    r["general_scheme_ids"] = list(general_scheme_ids)
    if reaction_name:
        r["reaction_name"] = fact(reaction_name, reaction_name, evidence)
    if transformation:
        r["transformation_type"] = inferred(
            transformation, evidence, "Transformation description inferred from the displayed reactant/product structures."
        )
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    if solvent:
        r["solvent_conditions"]["composition"] = fact(solvent, solvent, evidence)
    if solvent_ratio:
        r["solvent_conditions"]["mixture_ratio"] = fact(solvent_ratio, solvent_ratio, evidence)
    if total_volume:
        raw, normalized = total_volume
        r["solvent_conditions"]["total_volume"] = fact(raw, normalized, evidence)
    if concentration:
        raw, normalized = concentration
        r["solvent_conditions"]["concentration"] = fact(raw, normalized, evidence)
    if concentration_basis:
        r["solvent_conditions"]["concentration_basis"] = fact(concentration_basis, concentration_basis, evidence)
    if procedure:
        r["procedure"]["original_text"] = fact(procedure, procedure, evidence)
    if addition_order:
        r["procedure"]["addition_order"] = fact(addition_order, addition_order, evidence)
    if pretreatment:
        r["procedure"]["pretreatment"] = fact(pretreatment, pretreatment, evidence)
    r["classification"] = classification(class_changes, evidence, class_coverage)
    r["notes"] = list(notes)
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved) or ["Machine extraction requires expert source and structure review."],
    }
    r["completeness"] = "partial"
    return r


def direct_note(text, evidence):
    return fact(text, text, evidence)


def review_note(text, evidence):
    return inferred(text, evidence, "Extraction/review decision, not a quotation from the paper.")


def build_molecules():
    mols = [conflict_threonine()]
    mols += [
        structured_molecule("m10", "10", "(2S,3S)-methyl 3-methyloxirane-2-carboxylate",
                            "COC(=O)[C@H]1O[C@H]1C", [ev(2), ev(4), ev(2, "image"), ev(4, "image")],
                            atom_centers=[{"position": 2, "configuration": "S"}, {"position": 3, "configuration": "S"}],
                            assignment="(2S,3S) as printed in the experimental heading; Scheme 2 drawing retained."),
        structured_molecule("m11", "11", "(2R,3R)-methyl 3-chloro-2-hydroxybutanoate",
                            "COC(=O)[C@@H](O)[C@@H](C)Cl", [ev(2), ev(4), ev(5), ev(2, "image"), ev(4, "image"), ev(5, "image")],
                            atom_centers=[{"position": 2, "configuration": "R"}, {"position": 3, "configuration": "R"}],
                            assignment="(2R,3R) as printed in the experimental heading."),
        structured_molecule("m11a", "11a", "(2R,3R)-methyl 3-chloro-2-((trimethylsilyl)oxy)butanoate",
                            "COC(=O)[C@@H](O[Si](C)(C)C)[C@@H](C)Cl", [ev(5), ev(5, "image")],
                            atom_centers=[{"position": 2, "configuration": "R"}, {"position": 3, "configuration": "R"}],
                            assignment="(2R,3R) as printed in the experimental heading."),
        structured_molecule("m7a", "7a", "(2R,3R)-3-chloro-2-((tert-butyldimethylsilyl)oxy)butanal",
                            "C[C@@H](Cl)[C@@H](C=O)O[Si](C)(C)C(C)(C)C", [ev(2), ev(2, "image")],
                            atom_centers=[{"position": 2, "configuration": "R"}, {"position": 3, "configuration": "R"}],
                            assignment="(2R,3R)-7, P=TBS, as printed in Scheme 1/2."),
        structured_molecule("m7b", "7b", "(2R,3R)-3-chloro-2-((trimethylsilyl)oxy)butanal",
                            "C[C@@H](Cl)[C@@H](C=O)O[Si](C)(C)C", [ev(1), ev(2), ev(3), ev(5), ev(1, "image"), ev(2, "image"), ev(3, "image"), ev(5, "image")],
                            atom_centers=[{"position": 2, "configuration": "R"}, {"position": 3, "configuration": "R"}],
                            assignment="(2R,3R) as printed in the experimental heading."),
        structured_molecule("m6", "6", "tert-butyl 2-methylbuta-2,3-dienoate",
                            "CC(C)(C)OC(=O)C(C)=C=C", [ev(1), ev(2), ev(3), ev(4), ev(5), ev(6)],
                            atom_centers=[], assignment="Terminal CH2 allenoate 6; no configurationally stable allene axis is assigned.",
                            chirality="achiral", representation="partial",
                            identity_note="Terminal-CH2 allenoate 6 is not assigned an allene axis; representation is marked partial because the toolkit reports an unspecified cumulated-bond stereo candidate.",
                            allene_feature={"fragment_id": "allene_6", "raw": "C2=C3=C4 in the named chain",
                                            "substitution": "terminal CH2 allene", "axial": "not stereogenic",
                                            "environment": "2-methyl allenoate"}),
        structured_molecule("m12", "12", "(3S,5R,6R,7R)-tert-butyl 7-chloro-5-hydroxy-2-methyl-6-((trimethylsilyl)oxy)octa-2,3-dienoate",
                            "CC(=C=C[C@@H](O)[C@@H](O[Si](C)(C)C)[C@@H](C)Cl)C(=O)OC(C)(C)C",
                            [ev(3), ev(5), ev(3, "image"), ev(5, "image")],
                            atom_centers=[{"position": 5, "configuration": "R"}, {"position": 6, "configuration": "R"}, {"position": 7, "configuration": "R"}],
                            axial=[{"axis": "C2=C3=C4 allene", "configuration": "3S", "encoding": "not encoded by candidate SMILES"}],
                            assignment="(3S,5R,6R,7R) as printed; product described as a single isomer.",
                            representation="partial",
                            identity_note="Tetrahedral centers are encoded; the explicitly named 3S allene axis is retained textually and is not represented by the toolkit SMILES.",
                            allene_feature={"fragment_id": "allene_12", "raw": "octa-2,3-dienoate",
                                            "substitution": "substituted allene", "axial": "3S (source name)",
                                            "environment": "allenoate gamma-addition product"}),
        structured_molecule("m3", "3", "(R)-5-((1S,2R,3R)-3-chloro-1,2-dihydroxybutyl)-3-methylfuran-2(5H)-one",
                            "CC1=C[C@H]([C@@H](O)[C@@H](O)[C@@H](C)Cl)OC1=O",
                            [ev(1), ev(2), ev(3), ev(5), ev(1, "image"), ev(2, "image"), ev(3, "image"), ev(5, "image")],
                            atom_centers=[{"lactone_C5": "R"}, {"sidechain_C1": "S"}, {"sidechain_C2": "R"}, {"sidechain_C3": "R"}],
                            assignment="(R)-5-((1S,2R,3R)-...) as printed; structure also reported as X-ray confirmed."),
        structured_molecule("m1", "1", "reported hypoxylactone stereoisomer",
                            "CC1=C[C@H]([C@H](O)[C@H]2O[C@H]2C)OC1=O",
                            [ev(1), ev(2), ev(3), ev(5), ev(1, "image"), ev(2, "image"), ev(3, "image"), ev(5, "image")],
                            atom_centers=[{"lactone_C5": "R"}, {"sidechain_C1": "S"}, {"epoxide_C2": "R"}, {"epoxide_C3": "S"}],
                            assignment="Reported hypoxylactone structure 1 as drawn; the paper states its data are inconsistent with the natural product.",
                            identity_note="Candidate absolute descriptors were derived from the drawn precursor 3 and intramolecular epoxide closure; this derivation remains needs_review."),
        structured_molecule("m_ent11", "ent-11", "(2S,3S)-methyl 3-chloro-2-hydroxybutanoate",
                            "COC(=O)[C@H](O)[C@H](C)Cl", [ev(5), ev(5, "image")],
                            atom_centers=[{"position": 2, "configuration": "S"}, {"position": 3, "configuration": "S"}],
                            assignment="(2S,3S) stated by the ent-11a heading and precursor label."),
        structured_molecule("m_ent11a", "ent-11a", "(2S,3S)-methyl 2-((tert-butyldimethylsilyl)oxy)-3-chlorobutanoate",
                            "COC(=O)[C@H](O[Si](C)(C)C(C)(C)C)[C@H](C)Cl", [ev(5), ev(5, "image")],
                            atom_centers=[{"position": 2, "configuration": "S"}, {"position": 3, "configuration": "S"}],
                            assignment="(2S,3S) as printed in the experimental heading."),
        structured_molecule("m8a", "8a", "(2S,3S)-2-((tert-butyldimethylsilyl)oxy)-3-chlorobutanal",
                            "C[C@H](Cl)[C@H](C=O)O[Si](C)(C)C(C)(C)C", [ev(1), ev(2), ev(3), ev(4), ev(5), ev(6)],
                            atom_centers=[{"position": 2, "configuration": "S"}, {"position": 3, "configuration": "S"}],
                            assignment="(2S,3S) as printed in the experimental heading."),
        structured_molecule("m13", "13", "(3S,5R,6S,7S)-tert-butyl 6-((tert-butyldimethylsilyl)oxy)-7-chloro-5-hydroxy-2-methylocta-2,3-dienoate",
                            "CC(=C=C[C@@H](O)[C@H](O[Si](C)(C)C(C)(C)C)[C@H](C)Cl)C(=O)OC(C)(C)C",
                            [ev(3), ev(4), ev(6), ev(3, "image"), ev(4, "image"), ev(6, "image")],
                            atom_centers=[{"position": 5, "configuration": "R"}, {"position": 6, "configuration": "S"}, {"position": 7, "configuration": "S"}],
                            axial=[{"axis": "C2=C3=C4 allene", "configuration": "3S", "encoding": "not encoded by candidate SMILES"}],
                            assignment="(3S,5R,6S,7S) as printed; product described as a single isomer.",
                            representation="partial",
                            identity_note="Tetrahedral centers are encoded; the explicitly named 3S allene axis is retained textually and is not represented by the toolkit SMILES.",
                            allene_feature={"fragment_id": "allene_13", "raw": "octa-2,3-dienoate",
                                            "substitution": "substituted allene", "axial": "3S (source name)",
                                            "environment": "matched allenoate gamma-addition product"}),
        structured_molecule("m14", "14", "(R)-5-((1R,2S,3S)-2-((tert-butyldimethylsilyl)oxy)-3-chloro-1-hydroxybutyl)-3-methylfuran-2(5H)-one",
                            "CC1=C[C@H]([C@@H](O)[C@H](O[Si](C)(C)C(C)(C)C)[C@H](C)Cl)OC1=O",
                            [ev(3), ev(4), ev(6), ev(3, "image"), ev(4, "image"), ev(6, "image")],
                            atom_centers=[{"lactone_C5": "R"}, {"sidechain_C1": "R"}, {"sidechain_C2": "S"}, {"sidechain_C3": "S"}],
                            assignment="(R)-5-((1R,2S,3S)-...) as printed in the experimental heading."),
        structured_molecule("m4", "4", "(+)-xylogiblactone A",
                            "CC1=C[C@H]([C@H](O)[C@H](O)[C@H](C)Cl)OC1=O",
                            [ev(1), ev(2), ev(3), ev(4), ev(6), ev(1, "image"), ev(2, "image"), ev(3, "image"), ev(4, "image"), ev(6, "image")],
                            atom_centers=[{"lactone_C5": "R"}, {"sidechain_C1": "R"}, {"sidechain_C2": "S"}, {"sidechain_C3": "S"}],
                            assignment="Product 4 follows deprotection of explicitly named 14; source drawings specify the relative stereochemistry."),
        structured_molecule("m9", "9", "corrected (+)-hypoxylactone",
                            "CC1=C[C@H]([C@@H](O)[C@@H]2O[C@@H]2C)OC1=O",
                            [ev(1), ev(2), ev(3), ev(4), ev(6), ev(1, "image"), ev(2, "image"), ev(3, "image"), ev(4, "image"), ev(6, "image")],
                            atom_centers=[{"lactone_C5": "R"}, {"sidechain_C1": "R"}, {"epoxide_C2": "S"}, {"epoxide_C3": "R"}],
                            assignment="Corrected hypoxylactone structure 9 as drawn; physical data stated identical with the natural product.",
                            identity_note="Candidate absolute descriptors were derived from the source-labelled precursor 4 and intramolecular epoxide closure; this derivation remains needs_review."),
    ]

    unresolved_defs = [
        ("m5", "5", "(R,R)-chiral bromoborane 5", "specified", "not_applicable", "not_encoded", "(R,R) source label; full boron-complex graph not encoded."),
        ("m_ligand", "(1R,2R)-bis-sulfonamide", "(+)-(1R,2R)-1,2-diphenyl-1,2-bis-p-toluenesulfonylamide", "specified", "not_applicable", "not_encoded", "(1R,2R) as printed."),
        ("m_bbr3", "BBr3", "boron tribromide", "unresolved", "not_applicable", "unresolved", None),
        ("m_ipr2net", "i-Pr2NEt", "N,N-diisopropylethylamine", "unresolved", "not_applicable", "unresolved", None),
        ("m_ch2cl2", "CH2Cl2", "dichloromethane", "unresolved", "not_applicable", "unresolved", None),
        ("m_toluene", "toluene", "toluene", "unresolved", "not_applicable", "unresolved", None),
        ("m_phcf3", "PhCF3", "benzotrifluoride", "unresolved", "not_applicable", "unresolved", None),
        ("m_ph3aucl", "Ph3PAuCl", "triphenylphosphinegold(I) chloride", "unresolved", "unresolved", "unresolved", None),
        ("m_ph3aubf4", "Ph3AuBF4", "triphenylphosphinegold(I) tetrafluoroborate", "unresolved", "unresolved", "unresolved", None),
        ("m_ph3auasf6", "Ph3AuAsF6", "triphenylphosphinegold(I) hexafluoroarsenate", "unresolved", "unresolved", "unresolved", None),
        ("m_ph3pauntf2", "Ph3PAuNTf2", "triphenylphosphinegold(I) bis(trifluoromethanesulfonyl)imide", "unresolved", "unresolved", "unresolved", None),
        ("m_agbf4", "AgBF4", "silver tetrafluoroborate", "unresolved", "unresolved", "unresolved", None),
        ("m_agasf6", "AgAsF6", "silver hexafluoroarsenate", "unresolved", "unresolved", "unresolved", None),
        ("m_ag2o", "Ag2O", "silver(I) oxide", "unresolved", "not_applicable", "unresolved", None),
        ("m_nmp", "NMP", "N-methylpyrrolidone", "unresolved", "not_applicable", "unresolved", None),
        ("m_ticl4", "TiCl4", "titanium tetrachloride", "unresolved", "not_applicable", "unresolved", None),
        ("m_kbr", "KBr", "potassium bromide", "unresolved", "not_applicable", "unresolved", None),
        ("m_hbr", "HBr", "hydrogen bromide", "unresolved", "not_applicable", "unresolved", None),
        ("m_h2so4", "H2SO4", "sulfuric acid", "unresolved", "not_applicable", "unresolved", None),
        ("m_nano2", "NaNO2", "sodium nitrite", "unresolved", "not_applicable", "unresolved", None),
        ("m_koh", "KOH", "potassium hydroxide", "unresolved", "not_applicable", "unresolved", None),
        ("m_me2so4", "Me2SO4", "dimethyl sulfate", "unresolved", "not_applicable", "unresolved", None),
        ("m_18c6", "18-C-6", "18-crown-6", "unresolved", "not_applicable", "unresolved", None),
        ("m_meoh", "MeOH", "methanol", "unresolved", "not_applicable", "unresolved", None),
        ("m_etoac", "EtOAc", "ethyl acetate", "unresolved", "not_applicable", "unresolved", None),
        ("m_brine", "brine", "brine (composition not specified)", "mixture", "not_applicable", "not_encoded", "Composition not specified."),
        ("m_na2so4", "Na2SO4", "sodium sulfate", "unresolved", "not_applicable", "unresolved", None),
        ("m_mgso4", "MgSO4", "magnesium sulfate", "unresolved", "not_applicable", "unresolved", None),
        ("m_ether", "Et2O", "diethyl ether", "unresolved", "not_applicable", "unresolved", None),
        ("m_pentane", "pentane", "pentane", "unresolved", "not_applicable", "unresolved", None),
        ("m_hexanes", "hexane(s)", "hexane/hexanes", "mixture", "not_applicable", "not_encoded", "Source alternates singular/plural solvent label; composition not expanded."),
        ("m_silica", "SiO2", "silica gel", "unresolved", "not_applicable", "unresolved", None),
        ("m_imidazole", "imidazole", "imidazole", "unresolved", "not_applicable", "unresolved", None),
        ("m_tmscl", "TMSCl", "trimethylsilyl chloride", "unresolved", "not_applicable", "unresolved", None),
        ("m_dibalh", "Dibal-H", "diisobutylaluminium hydride", "unresolved", "unresolved", "unresolved", None),
        ("m_cyclohexane", "cyclohexane", "cyclohexane", "unresolved", "not_applicable", "unresolved", None),
        ("m_rochelle", "sat. potassium sodium tartrate", "saturated aqueous potassium sodium tartrate", "mixture", "not_applicable", "not_encoded", "Concentration not reported."),
        ("m_buffer7", "pH 7 buffer", "aqueous pH 7 buffer", "mixture", "not_applicable", "not_encoded", "Buffer composition not reported."),
        ("m_tbscl", "TBSCl", "tert-butyldimethylsilyl chloride", "unresolved", "not_applicable", "unresolved", None),
        ("m_dmf", "DMF", "N,N-dimethylformamide", "unresolved", "not_applicable", "unresolved", None),
        ("m_thf", "THF", "tetrahydrofuran", "unresolved", "not_applicable", "unresolved", None),
        ("m_water", "H2O", "water", "unresolved", "not_applicable", "unresolved", None),
        ("m_celite", "Celite", "Celite filter aid", "mixture", "not_applicable", "not_encoded", "Commercial filter aid; composition not expanded."),
        ("m_hcl", "10% aq HCl", "10% aqueous hydrochloric acid", "mixture", "not_applicable", "not_encoded", "Aqueous formulation retained as reported."),
        ("m_nahco3", "NaHCO3", "sodium bicarbonate", "unresolved", "not_applicable", "unresolved", None),
    ]
    common_ev = [ev(2), ev(3), ev(4), ev(5), ev(6)]
    for mid, label, name, chir, geom, rep, assign in unresolved_defs:
        evidence = common_ev
        if mid in {"m5", "m_ligand", "m_bbr3", "m_ipr2net"}:
            evidence = [ev(2), ev(3), ev(4), ev(5), ev(6)]
        mols.append(unresolved_molecule(mid, label, name, evidence, chirality=chir,
                                        geometry=geom, representation=rep, assignment=assign))
    return mols


def table1_record(entry, t1, time1, t2, time2, solvent_mid, solvent_name, yld, ee_value):
    evidence = [ev(3), ev(3, "image")]
    if entry == 4:
        evidence += [ev(2), ev(5), ev(2, "image"), ev(5, "image")]
    ps = [
        participant(f"t1e{entry}_6", "m6", "reactant", evidence,
                    mass=("0.16 g", quantity(0.16, "g")) if entry == 4 else None,
                    amount=("1.02 mmol", quantity(1.02, "mmol")) if entry == 4 else None,
                    equivalents=("1.06 equiv", quantity(1.06, "equiv")), reference="m5"),
        participant(f"t1e{entry}_7b", "m7b", "reactant", evidence,
                    mass=("194 mg", quantity(194, "mg")) if entry == 4 else None,
                    amount=("1.0 mmol", quantity(1.0, "mmol")) if entry == 4 else None,
                    equivalents=("1.04 equiv", quantity(1.04, "equiv")), reference="m5"),
        participant(f"t1e{entry}_5", "m5", "reagent", evidence,
                    equivalents=("1.0 equiv", quantity(1.0, "equiv")), reference="m5"),
        participant(f"t1e{entry}_base", "m_ipr2net", "reagent", evidence,
                    volume=("0.35 mL", quantity(0.35, "mL")) if entry == 4 else None,
                    mass=("0.26 g", quantity(0.26, "g")) if entry == 4 else None,
                    amount=("2.01 mmol", quantity(2.01, "mmol")) if entry == 4 else None,
                    equivalents=("2.1 equiv", quantity(2.1, "equiv")), reference="m5"),
        participant(f"t1e{entry}_solv", solvent_mid, "solvent", evidence,
                    volume=("7 mL initial CH2Cl2", quantity(7, "mL")) if entry == 4 else None),
    ]
    if entry == 4:
        ps += [
            participant("t1e4_ligand", "m_ligand", "reagent", evidence,
                        mass=("0.50 g", quantity(0.50, "g")), amount=("0.96 mmol", quantity(0.96, "mmol"))),
            participant("t1e4_bbr3", "m_bbr3", "reagent", evidence,
                        volume=("1.0 mL of 1 M solution", quantity(1.0, "mL")), amount=("1.0 mmol", quantity(1.0, "mmol"))),
            participant("t1e4_ch2cl2_exchange", "m_ch2cl2", "solvent", evidence, volume=("5 mL solvent exchange", quantity(5, "mL"))),
            participant("t1e4_ch2cl2_redissolve", "m_ch2cl2", "solvent", evidence, volume=("3 mL redissolution", quantity(3, "mL"))),
            participant("t1e4_ch2cl2_6carrier", "m_ch2cl2", "solvent", evidence, volume=("1 mL carrier for 6/i-Pr2NEt", quantity(1, "mL"))),
            participant("t1e4_ch2cl2_7bcarrier", "m_ch2cl2", "solvent", evidence, volume=("1 mL carrier for 7b", quantity(1, "mL"))),
            participant("t1e4_buffer", "m_buffer7", "workup_agent", evidence, stage="workup", volume=("7 mL", quantity(7, "mL"))),
            participant("t1e4_ch2cl2_quench", "m_ch2cl2", "workup_agent", evidence, stage="workup", volume=("ca. 10 mL added after buffer", quantity(10, "mL", "approx", "ca. 10 mL"))),
            participant("t1e4_ch2cl2_extract1", "m_ch2cl2", "workup_agent", evidence, stage="workup", volume=("ca. 10 mL first extraction", quantity(10, "mL", "approx", "ca. 10 mL"))),
            participant("t1e4_ch2cl2_extract2", "m_ch2cl2", "workup_agent", evidence, stage="workup", volume=("ca. 10 mL second extraction", quantity(10, "mL", "approx", "ca. 10 mL"))),
            participant("t1e4_brine", "m_brine", "workup_agent", evidence, stage="workup"),
            participant("t1e4_na2so4", "m_na2so4", "workup_agent", evidence, stage="workup"),
            participant("t1e4_ether", "m_ether", "workup_agent", evidence, stage="workup", volume=("ca. 20 mL", quantity(20, "mL", "approx", "ca. 20 mL"))),
            participant("t1e4_silica", "m_silica", "purification_agent", evidence, stage="purification"),
        ]
    if entry != 1:
        ps.append(participant(f"t1e{entry}_12", "m12", "product", evidence, product_role="main product"))
    ss = [step(f"t1e{entry}_enolization", evidence, pids=[f"t1e{entry}_6", f"t1e{entry}_5", f"t1e{entry}_base"],
               operation="allenoate activation/enolization", temperature=(f"{t1} °C", quantity(t1, "°C")),
               time=(f"{time1} min", quantity(time1, "min")), atmosphere="dried N2")]
    if entry == 4:
        ss = [
            step("t1e4_make5_0", evidence, pids=["t1e4_ligand", "t1e4_bbr3", "t1e4_solv"],
                 operation="BBr3 added to the bis-sulfonamide in 7 mL dry CH2Cl2 and stirred",
                 temperature=("0 °C", quantity(0, "°C")), time=("30 min", quantity(30, "min")), atmosphere="N2"),
            step("t1e4_make5_25", evidence, pids=["t1e4_ligand", "t1e4_bbr3"],
                 operation="warm to 25 °C and stir", temperature=("25 °C", quantity(25, "°C")),
                 time=("1.5 h", quantity(1.5, "h")), atmosphere="N2"),
            step("t1e4_solvent_exchange", evidence, pids=["t1e4_ch2cl2_exchange"],
                 operation="concentrate at 1 mmHg; add 5 mL dry CH2Cl2 and evaporate again"),
            step("t1e4_enolization", evidence, pids=["t1e4_6", "t1e4_5", "t1e4_base", "t1e4_ch2cl2_redissolve", "t1e4_ch2cl2_6carrier"],
                 operation="redissolve 5 in 3 mL CH2Cl2; add 6/i-Pr2NEt in 1 mL CH2Cl2 dropwise and stir",
                 temperature=("−50 °C", quantity(-50, "°C")), time=("30 min", quantity(30, "min")), atmosphere="N2"),
            step("t1e4_add_7b", evidence, pids=["t1e4_7b", "t1e4_ch2cl2_7bcarrier"], operation="precooled 7b in 1 mL CH2Cl2 added by cannula along the flask wall",
                 temperature=("−50 °C", quantity(-50, "°C")), time=None, atmosphere="N2"),
            step("t1e4_warm_hold", evidence, pids=["t1e4_6", "t1e4_7b", "t1e4_5", "t1e4_base"],
                 operation="allowed to warm from −50 °C to −10 °C and proceed",
                 temperature=("−50 to −10 °C", range_quantity(-50, -10, "°C", "−50 to −10 °C")),
                 time=("4 h", quantity(4, "h")), atmosphere="N2"),
            step("t1e4_workup", evidence, stage="workup", pids=["t1e4_buffer", "t1e4_ch2cl2_quench", "t1e4_ch2cl2_extract1", "t1e4_ch2cl2_extract2", "t1e4_brine", "t1e4_na2so4", "t1e4_ether"],
                 operation="quench with pH 7 buffer; CH2Cl2 extraction; brine wash; Na2SO4 drying; evaporate; ether precipitation and 0 °C filtration",
                 temperature=("0 °C for precipitation", quantity(0, "°C")), time=("20 min", quantity(20, "min"))),
            step("t1e4_purification", evidence, stage="purification", pids=["t1e4_silica"],
                 operation="flash column chromatography, hexane/EtOAc 30:1"),
        ]
    else:
        ss.append(step(f"t1e{entry}_reaction", evidence,
                       pids=[f"t1e{entry}_6", f"t1e{entry}_7b", f"t1e{entry}_5", f"t1e{entry}_base"],
                       operation="aldehyde addition/reaction stage",
                       temperature=(f"{t2} °C", quantity(t2, "°C")),
                       time=(f"{time2} h", quantity(time2, "h")), atmosphere="dried N2"))
        if entry != 1:
            ps.append(participant(f"t1e{entry}_silica", "m_silica", "purification_agent", evidence, stage="purification"))
            ss.append(step(f"t1e{entry}_purification", evidence, stage="purification", pids=[f"t1e{entry}_silica"], operation="flash chromatography"))
    outs = []
    if entry == 1:
        outs.append(outcome("no_reaction", "no_reaction", ["m6", "m7b"], evidence,
                            raw="NR", normalized=True, measurement=None,
                            scope="Table 1 reaction outcome"))
        changes, coverage = [], "not_applicable"
    else:
        if isinstance(yld, tuple):
            raw_y, norm_y = yld
        else:
            raw_y, norm_y = f"{yld}%", quantity(yld, "%")
        outs.append(outcome(f"t1e{entry}_yield", "yield", ["m12"], evidence,
                            raw=raw_y, normalized=norm_y, yield_type="isolated",
                            measurement="isolated after flash chromatography",
                            scope="Table 1 entry product 12", product_role="main product"))
        if entry == 4:
            outs.append(product_amount("t1e4_product_amount", "m12", evidence,
                                       "115 mg", quantity(115, "mg"), "0.33 mmol", quantity(0.33, "mmol")))
        if ee_value is not None:
            raw_ee, norm_ee = ee_value
            outs.append(outcome(f"t1e{entry}_ee", "selectivity", ["m12"], evidence,
                                raw=raw_ee, normalized=norm_ee,
                                measurement="NMR analysis",
                                scope={"metric": "ee", "major_isomer": "12", "minor_isomer": "not assigned"},
                                product_role="main product"))
        changes = [chg("C", "C", 0, 1, "formed"), chg("C", "O", 2, 1, "order_decreased"), chg("O", "H", 0, 1, "hydrogen_gain")]
        coverage = "partial"
    notes = []
    unresolved = ["Candidate structures, including product-12 allene axial assignment, require expert review."]
    if entry == 4:
        notes.append(direct_note("The exact experimental paragraph and Table 1 entry 4 describe the same 34% reaction and are merged as one record.", evidence))
        unresolved.append("The 3S allene axis is retained textually but is not encoded in the candidate SMILES.")
    if entry == 7:
        notes.append(direct_note("Table 1 reports ee as NA for entry 7.", evidence))
    return reaction(
        f"table1_entry_{entry}", evidence, categories=["optimization"],
        location={"page": 3, "region": f"Table 1 entry {entry}"},
        reaction_name="allenoate γ-addition of 6 with 7b",
        transformation="asymmetric allenoate gamma-addition to an aldehyde",
        participants=ps, steps=ss, outcomes=outs, class_changes=changes,
        class_coverage=coverage,
        procedure="Table 1 conditions; entry 4 additionally uses the detailed page-5 procedure." if entry == 4 else "Table 1 row and footnotes.",
        notes=notes, unresolved=unresolved,
        solvent=solvent_name,
        solvent_ratio="CH2Cl2/PhCF3 = 4:1" if entry == 8 else None,
        general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    )


def table2_record(entry, catalyst_mid, catalyst_label, molpct, solvent_mid, solvent_name, temp, hours, yld):
    evidence = [ev(3), ev(3, "image")]
    if entry == 4:
        evidence += [ev(5), ev(5, "image")]
    ps = [
        participant(f"t2e{entry}_12", "m12", "reactant", evidence,
                    mass=("67.5 mg", quantity(67.5, "mg")) if entry == 4 else None,
                    amount=("0.19 mmol", quantity(0.19, "mmol")) if entry == 4 else None),
        participant(f"t2e{entry}_cat", catalyst_mid, "catalyst", evidence,
                    mass=("30 mg", quantity(30, "mg")) if entry == 4 else None,
                    amount=("0.019 mmol", quantity(0.019, "mmol")) if entry == 4 else None,
                    mol_percent=(f"{molpct} mol%", quantity(molpct, "mol%")), reference="m12"),
        participant(f"t2e{entry}_solv", solvent_mid, "solvent", evidence,
                    volume=("2 mL", quantity(2, "mL")) if entry == 4 else None),
        participant(f"t2e{entry}_3", "m3", "product", evidence, product_role="main product"),
        participant(f"t2e{entry}_silica", "m_silica", "purification_agent", evidence, stage="purification"),
    ]
    if entry == 2:
        ps += [participant("t2e2_ph3aucl", "m_ph3aucl", "reagent", evidence),
               participant("t2e2_agbf4", "m_agbf4", "reagent", evidence)]
    if entry == 3:
        ps += [participant("t2e3_ph3aucl", "m_ph3aucl", "reagent", evidence),
               participant("t2e3_agasf6", "m_agasf6", "reagent", evidence)]
    ss = [step(f"t2e{entry}_reaction", evidence,
               pids=[f"t2e{entry}_12", f"t2e{entry}_cat", f"t2e{entry}_solv"],
               operation="gold-catalyzed oxacyclization",
               temperature=(f"{temp} °C", quantity(temp, "°C")),
               time=(f"{hours} h", quantity(hours, "h")), atmosphere="dried N2")]
    if entry == 4:
        ss.append(step("t2e4_concentrate", evidence, stage="workup", pids=[], operation="concentrated under reduced pressure"))
    ss.append(step(f"t2e{entry}_purification", evidence, stage="purification",
                   pids=[f"t2e{entry}_silica"],
                   operation="flash chromatography" + (", hexane/EtOAc 3:1 to 1:1" if entry == 4 else "")))
    raw_y, norm_y = yld if isinstance(yld, tuple) else (f"{yld}%", quantity(yld, "%"))
    outs = [outcome(f"t2e{entry}_yield", "yield", ["m3"], evidence,
                    raw=raw_y, normalized=norm_y, yield_type="isolated",
                    measurement="isolated after flash chromatography",
                    scope=f"Table 2 entry {entry}, product 3", product_role="main product")]
    if entry == 4:
        outs.append(product_amount("t2e4_product_amount", "m3", evidence,
                                   "28 mg", quantity(28, "mg"), "0.13 mmol", quantity(0.13, "mmol")))
    notes = []
    if entry in (2, 3):
        salt = "AgBF4" if entry == 2 else "AgAsF6"
        notes.append(direct_note(f"{catalyst_label} was prepared in situ by stirring a 1:1 mixture of Ph3AuCl and {salt}.", evidence))
    if entry == 4:
        notes.append(direct_note("Table 2 entry 4 and the detailed page-5 preparation of 3 are the same 68% experiment and are merged.", evidence))
    return reaction(
        f"table2_entry_{entry}", evidence, categories=["optimization"],
        location={"page": 3, "region": f"Table 2 entry {entry}"},
        reaction_name="oxacyclization of 12 to 3 via gold catalysis",
        transformation="gold-catalyzed oxacyclization of an allenoate to a gamma-lactone",
        participants=ps, steps=ss, outcomes=outs,
        class_changes=[chg("C", "O", 0, 1, "formed"), chg("O", "H", 0, 1, "hydrogen_gain")],
        procedure="Table 2 conditions and footnotes; entry 4 additionally uses the detailed page-5 procedure." if entry == 4 else "Table 2 row and footnotes.",
        notes=notes,
        unresolved=["Net ring-forming C–O and TMS-deprotection O–H changes are partial; complete atom mapping and coproducts are unresolved."],
        solvent=solvent_name,
        total_volume=("2 mL", quantity(2, "mL")) if entry == 4 else None,
        concentration=("0.1 M", quantity(0.1, "M")),
        concentration_basis="substrate 12 (Table 2 footnote)",
        general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    )


def build_reactions():
    reactions = []
    # Conceptual strategy and proposed stereochemical models are kept separate from performed experiments.
    reactions.append(reaction(
        "scheme1_synthetic_strategy", [ev(2), ev(2, "image")], record_type="generic_scheme",
        location={"page": 2, "region": "Scheme 1"},
        reaction_name="Synthetic strategy toward reported/corrected hypoxylactone",
        transformation="conceptual allenoate gamma-addition, gold cyclization, and epoxide closure sequence",
        participants=[participant("s1_6", "m6", "reactant", [ev(2), ev(2, "image")]),
                      participant("s1_5", "m5", "reagent", [ev(2), ev(2, "image")]),
                      participant("s1_7a", "m7a", "reactant", [ev(2), ev(2, "image")]),
                      participant("s1_8a", "m8a", "reactant", [ev(2), ev(2, "image")]),
                      participant("s1_3", "m3", "product", [ev(2), ev(2, "image")], product_role="strategy intermediate"),
                      participant("s1_4", "m4", "product", [ev(2), ev(2, "image")], product_role="strategy intermediate"),
                      participant("s1_1", "m1", "product", [ev(2), ev(2, "image")], product_role="reported target"),
                      participant("s1_9", "m9", "product", [ev(2), ev(2, "image")], product_role="corrected target")],
        class_changes=(), class_coverage="unresolved",
        notes=[direct_note("Scheme 1 is a synthetic strategy, not a single performed experiment.", [ev(2), ev(2, "image")])],
        unresolved=["Scheme 1 contains conceptual arrows; concrete experiments are represented separately."],
    ))
    reactions.append(reaction(
        "scheme2_silylation_reduction_summary", [ev(2), ev(2, "image")], record_type="unresolved",
        location={"page": 2, "region": "Scheme 2 final arrow"},
        reaction_name="silylation/reduction summary from 11 to aldehydes 7a/7b",
        transformation="hydroxyl protection followed by ester-to-aldehyde reduction",
        participants=[participant("s2sum_11", "m11", "reactant", [ev(2), ev(2, "image")]),
                      participant("s2sum_7a", "m7a", "product", [ev(2), ev(2, "image")], product_role="P=TBS branch"),
                      participant("s2sum_7b", "m7b", "product", [ev(2), ev(2, "image")], product_role="P=TMS branch"),
                      participant("s2sum_imid", "m_imidazole", "reagent", [ev(2), ev(2, "image")]),
                      participant("s2sum_dibal", "m_dibalh", "reagent", [ev(2), ev(2, "image")])],
        steps=[step("s2sum_protect", [ev(2), ev(2, "image")], pids=["s2sum_11", "s2sum_imid"], operation="R3SiCl/imidazole protection"),
               step("s2sum_reduce", [ev(2), ev(2, "image")], pids=["s2sum_dibal"], operation="Dibal-H reduction",
                    temperature=("−78 °C", quantity(-78, "°C")))],
        outcomes=[outcome("s2sum_range", "other", ["m7a", "m7b"], [ev(2), ev(2, "image")],
                          raw="75–83%", normalized=range_quantity(75, 83, "%", "75–83%"),
                          scope="Scheme 2 range; not allocated to a specific branch")],
        class_changes=(), class_coverage="unresolved",
        notes=[review_note("The printed 75–83% range is not assigned to 7a versus 7b and is not substituted for the separate page-5 step yields.", [ev(2), ev(2, "image")])],
        unresolved=["Scheme 2 does not allocate 75–83% to the TBS and TMS branches; exact 7a procedure is not supplied in pages 1–6."],
    ))
    reactions.append(reaction(
        "scheme3_4_stereochemical_models", [ev(3), ev(4), ev(3, "image"), ev(4, "image")],
        record_type="proposed_mechanism", categories=[],
        location={"pages": [3, 4], "regions": ["Scheme 3 models B/C", "Scheme 4 model D"]},
        reaction_name="Felkin–Anh/anti-Felkin–Anh stereochemical models for allenoate gamma-addition",
        transformation="proposed stereochemical rationale",
        participants=[participant("model_5", "m5", "reagent", [ev(3), ev(4)]),
                      participant("model_7b", "m7b", "reactant", [ev(3)]),
                      participant("model_8a", "m8a", "reactant", [ev(4)]),
                      participant("model_12", "m12", "product", [ev(3)], product_role="mismatched-route product"),
                      participant("model_13", "m13", "product", [ev(4)], product_role="matched-route product")],
        class_changes=(), class_coverage="unresolved",
        notes=[direct_note("The authors state that exact mechanistic aspects have not been rigorously elucidated.", [ev(2), ev(3)])],
        unresolved=["Models B, C, and D are proposed stereochemical rationales, not independently verified mechanistic experiments."],
    ))

    # Real failed experiment with 7a, kept distinct from Table 1.
    evidence = [ev(2), ev(2, "image")]
    reactions.append(reaction(
        "failed_gamma_addition_6_7a", evidence, categories=["control"],
        location={"page": 2, "region": "right-column paragraph below Scheme 2"},
        reaction_name="attempted allenoate gamma-addition of 6 to 7a",
        transformation="attempted allenoate gamma-addition; no conversion to gamma-adduct",
        participants=[participant("fail_6", "m6", "reactant", evidence),
                      participant("fail_7a", "m7a", "reactant", evidence),
                      participant("fail_5", "m5", "reagent", evidence),
                      participant("fail_base", "m_ipr2net", "reagent", evidence,
                                  equivalents=("2.1 equiv", quantity(2.1, "equiv")), reference="m5"),
                      participant("fail_solv", "m_ch2cl2", "solvent", evidence)],
        steps=[step("fail_activation", evidence, pids=["fail_6", "fail_5", "fail_base"],
                    operation="reaction of 6 with in-situ prepared (R,R)-5 in CH2Cl2",
                    temperature=("−50 °C", quantity(-50, "°C")), time=("30 min", quantity(30, "min"))),
               step("fail_add_7a", evidence, pids=["fail_7a"], operation="addition of 7a",
                    temperature=("−20 °C", quantity(-20, "°C")))],
        outcomes=[outcome("fail_no_conversion", "no_reaction", ["m6", "m7a"], evidence,
                          raw="conversion to the γ-adduct completely failed", normalized=True,
                          scope="attempted gamma-addition")],
        class_changes=(), class_coverage="not_applicable",
        procedure="6, in-situ (R,R)-5 and i-Pr2NEt at −50 °C for 30 min in CH2Cl2, followed by 7a addition at −20 °C.",
        unresolved=["No product identity, assay method, or post-addition duration is reported."],
        solvent="CH2Cl2",
    ))

    # Preparation of 10: source conflicts are explicitly preserved.
    evidence = [ev(2), ev(4), ev(2, "image"), ev(4, "image")]
    t_alt_scheme = fact("−78 °C", quantity(-78, "°C"), [ev(2), ev(2, "image")])
    t_alt_exp = fact("−30 °C", quantity(-30, "°C"), [ev(4), ev(4, "image")])
    temp_conflict = fact(["−78 °C", "−30 °C"], None, evidence, status="conflict",
                         alternatives=[t_alt_scheme, t_alt_exp])
    prep10_steps = [
        step("prep10_diazotization_addition", evidence, pids=["p10_thr", "p10_kbr", "p10_h2so4", "p10_nano2"],
             operation="NaNO2 added in portions over 30 min to threonine/KBr in aqueous H2SO4",
             temperature=("0 °C", quantity(0, "°C")), time=("30 min addition", quantity(30, "min"))),
        step("prep10_diazotization_hold", evidence, pids=["p10_thr", "p10_kbr", "p10_h2so4", "p10_nano2"],
             operation="stir at the same temperature", temperature=("0 °C", quantity(0, "°C")),
             time=("2 h", quantity(2, "h"))),
        step("prep10_first_workup", evidence, stage="workup", pids=["p10_etoac", "p10_brine", "p10_na2so4"],
             operation="EtOAc extraction three times; water then brine wash; Na2SO4 drying; evaporation; crude bromoacid carried forward"),
        step("prep10_koh_addition", evidence, pids=["p10_koh", "p10_meoh"],
             operation="KOH in absolute MeOH slowly added to crude bromoacid in absolute MeOH",
             time=None),
        step("prep10_koh_temperature_conflict", evidence, pids=["p10_koh", "p10_meoh"],
             operation="cold epoxidation stage with conflicting printed temperatures"),
        step("prep10_cold_hold", [ev(4), ev(4, "image")], pids=["p10_koh", "p10_meoh"],
             operation="stir cold after KOH addition", temperature=("−30 °C", quantity(-30, "°C")),
             time=("1 h", quantity(1, "h"))),
        step("prep10_rt_hold", [ev(4), ev(4, "image")], pids=["p10_koh", "p10_meoh"],
             operation="warm to room temperature and stir", temperature=("rt", "rt"), time=("12 h", quantity(12, "h"))),
        step("prep10_methylation", [ev(4), ev(4, "image")], pids=["p10_me2so4", "p10_18c6", "p10_ch2cl2"],
             operation="Me2SO4 slowly added to crude salt/18-crown-6 in CH2Cl2",
             temperature=("rt", "rt"), time=("1 h", quantity(1, "h"))),
        step("prep10_final_workup", [ev(4), ev(4, "image")], stage="workup", pids=["p10_ether"],
             operation="Et2O precipitation; filtration; careful vacuum concentration"),
        step("prep10_purification", [ev(4), ev(4, "image")], stage="purification", pids=["p10_silica", "p10_pentane", "p10_ether_pur"],
             operation="column chromatography, pentane/Et2O 8:1"),
    ]
    prep10_steps[4]["temperature"] = temp_conflict
    prep10_ps = [
        participant("p10_thr", "m_threonine_source_conflict", "reactant", evidence,
                    mass=("2 g", quantity(2, "g")), amount=("16.8 mmol", quantity(16.8, "mmol"))),
        participant("p10_kbr", "m_kbr", "reagent", evidence, mass=("7.03 g", quantity(7.03, "g")), amount=("59.1 mmol", quantity(59.1, "mmol"))),
        participant("p10_hbr", "m_hbr", "reagent", [ev(2), ev(2, "image")]),
        participant("p10_h2so4", "m_h2so4", "reagent", evidence, volume=("35 mL", quantity(35, "mL")), amount=("43.7 mmol", quantity(43.7, "mmol"))),
        participant("p10_nano2", "m_nano2", "reagent", evidence, mass=("1.88 g", quantity(1.88, "g")), amount=("27.2 mmol", quantity(27.2, "mmol"))),
        participant("p10_meoh", "m_meoh", "solvent", evidence, volume=("28 mL + 9.3 mL", quantity(37.3, "mL"))),
        participant("p10_koh", "m_koh", "reagent", evidence, mass=("1.98 g", quantity(1.98, "g")), amount=("33.6 mmol", quantity(33.6, "mmol"))),
        participant("p10_18c6", "m_18c6", "reagent", evidence, mass=("3.95 g", quantity(3.95, "g")), amount=("14.9 mmol", quantity(14.9, "mmol"))),
        participant("p10_me2so4", "m_me2so4", "reagent", evidence, volume=("1.4 mL", quantity(1.4, "mL")), amount=("14.9 mmol", quantity(14.9, "mmol"))),
        participant("p10_ch2cl2", "m_ch2cl2", "solvent", evidence, volume=("21 mL", quantity(21, "mL"))),
        participant("p10_etoac", "m_etoac", "workup_agent", evidence, stage="workup"),
        participant("p10_brine", "m_brine", "workup_agent", evidence, stage="workup"),
        participant("p10_na2so4", "m_na2so4", "workup_agent", evidence, stage="workup"),
        participant("p10_ether", "m_ether", "workup_agent", evidence, stage="workup"),
        participant("p10_silica", "m_silica", "purification_agent", evidence, stage="purification"),
        participant("p10_pentane", "m_pentane", "purification_agent", evidence, stage="purification"),
        participant("p10_ether_pur", "m_ether", "purification_agent", evidence, stage="purification"),
        participant("p10_prod", "m10", "product", evidence, product_role="main product"),
    ]
    reactions.append(reaction(
        "prep_10", evidence, categories=["other"], location={"pages": [2, 4], "regions": ["Scheme 2 first arrow", "Preparation of 10"]},
        reaction_name="Preparation of (2S,3S)-methyl 3-methyloxirane-2-carboxylate (10)",
        transformation="three-step diazotization/bromination, epoxide formation, and methylation sequence",
        participants=prep10_ps, steps=prep10_steps,
        outcomes=[outcome("prep10_yield", "yield", ["m10"], evidence, raw="71% in 3 steps", normalized=quantity(71, "%"),
                          yield_type="three-step isolated yield", measurement="isolated after column chromatography",
                          scope={"start": "source-conflicted threonine label", "end": "10", "steps": 3}, product_role="main product"),
                  product_amount("prep10_amount", "m10", evidence, "1.38 g", quantity(1.38, "g"), "11.9 mmol", quantity(11.9, "mmol"))],
        class_changes=(), class_coverage="unresolved",
        procedure="Full page-4 experimental paragraph; Scheme 2 supplies an additional, conflicting summary.",
        notes=[review_note("Source conflict preserved: Scheme 2 says D-threonine and −78 °C; the detailed experiment says L-threonine and −30 °C.", evidence),
               review_note("HBr appears in Scheme 2 but is absent from the detailed page-4 reagent list; it is retained as a participant with unreported dose.", evidence)],
        unresolved=["Starting threonine enantiomer conflicts between Scheme 2 and the experimental paragraph.",
                    "Cold KOH-stage temperature conflicts between −78 °C (Scheme 2) and −30 °C (experiment).",
                    "Overall multistep bond mapping is unresolved; no intermediates or hidden coproducts were invented."],
    ))

    # 10 -> 11.
    evidence = [ev(2), ev(4), ev(2, "image"), ev(4, "image")]
    reactions.append(reaction(
        "prep_11", evidence, categories=["other"], location={"pages": [2, 4], "regions": ["Scheme 2 second arrow", "Preparation of 11"]},
        reaction_name="Preparation of (2R,3R)-methyl 3-chloro-2-hydroxybutanoate (11)",
        transformation="TiCl4/NMP epoxide opening to a chlorohydrin",
        participants=[participant("p11_10", "m10", "reactant", evidence, mass=("2.2 g", quantity(2.2, "g")), amount=("18.9 mmol", quantity(18.9, "mmol"))),
                      participant("p11_nmp", "m_nmp", "reagent", evidence, volume=("1.81 mL", quantity(1.81, "mL")), amount=("18.9 mmol", quantity(18.9, "mmol")), equivalents=("1.0 equiv", quantity(1.0, "equiv")), reference="m10"),
                      participant("p11_ticl4", "m_ticl4", "reagent", evidence, volume=("23.4 mL of 1.0 M solution", quantity(23.4, "mL")), amount=("23.44 mmol", quantity(23.44, "mmol")), equivalents=("1.2 equiv", quantity(1.2, "equiv")), reference="m10"),
                      participant("p11_solv", "m_ch2cl2", "solvent", evidence, volume=("14 mL", quantity(14, "mL"))),
                      participant("p11_nahco3", "m_nahco3", "workup_agent", evidence, stage="workup", volume=("50 mL of 7% aq solution", quantity(50, "mL"))),
                      participant("p11_brine", "m_brine", "workup_agent", evidence, stage="workup"),
                      participant("p11_mgso4", "m_mgso4", "workup_agent", evidence, stage="workup"),
                      participant("p11_silica", "m_silica", "purification_agent", evidence, stage="purification"),
                      participant("p11_prod", "m11", "product", evidence, product_role="main product")],
        steps=[step("p11_add", evidence, pids=["p11_10", "p11_nmp", "p11_ticl4", "p11_solv"], operation="TiCl4 added dropwise", temperature=("−78 °C", quantity(-78, "°C"))),
               step("p11_hold", evidence, pids=["p11_10", "p11_nmp", "p11_ticl4"], operation="reaction hold", temperature=("−78 °C", quantity(-78, "°C")), time=("1.5 h", quantity(1.5, "h"))),
               step("p11_quench", evidence, stage="workup", pids=["p11_nahco3"], operation="quench with 7% aqueous NaHCO3 and stir", temperature=("rt", "rt"), time=("30 min", quantity(30, "min"))),
               step("p11_workup", evidence, stage="workup", pids=["p11_brine", "p11_mgso4"], operation="CH2Cl2 extraction twice; brine wash; MgSO4 drying; vacuum concentration"),
               step("p11_purify", evidence, stage="purification", pids=["p11_silica"], operation="column chromatography, hexane/EtOAc 15:1 to 7:1")],
        outcomes=[outcome("p11_yield", "yield", ["m11"], evidence, raw="89%", normalized=quantity(89, "%"), yield_type="isolated", measurement="isolated after column chromatography", scope="10 to 11", product_role="main product"),
                  product_amount("p11_amount", "m11", evidence, "2.58 g", quantity(2.58, "g"), "16.89 mmol", quantity(16.89, "mmol"))],
        class_changes=[chg("C", "O", 1, 0, "broken"), chg("C", "Cl", 0, 1, "formed"), chg("O", "H", 0, 1, "hydrogen_gain")],
        procedure="Full page-4 preparation of 11.", solvent="CH2Cl2",
        concentration=("1.3 M", quantity(1.3, "M")), concentration_basis="compound 10",
        unresolved=["Partial epoxide-opening bond inventory; chloride/proton atom sources and coproducts are not fully mapped."],
    ))

    # 11 -> 11a and 11a -> 7b.
    e5 = [ev(5), ev(5, "image")]
    reactions.append(reaction(
        "prep_11a", e5, categories=["other"], location={"page": 5, "region": "Preparation of 11a"},
        reaction_name="Preparation of (2R,3R)-methyl 3-chloro-2-((trimethylsilyl)oxy)butanoate (11a)",
        transformation="trimethylsilyl protection of an alcohol",
        participants=[participant("p11a_11", "m11", "reactant", e5, mass=("2.58 g", quantity(2.58, "g")), amount=("16.89 mmol", quantity(16.89, "mmol"))),
                      participant("p11a_imid", "m_imidazole", "reagent", e5, mass=("1.27 g", quantity(1.27, "g")), amount=("18.58 mmol", quantity(18.58, "mmol"))),
                      participant("p11a_tmscl", "m_tmscl", "reagent", e5, volume=("2.37 mL", quantity(2.37, "mL")), amount=("18.58 mmol", quantity(18.58, "mmol"))),
                      participant("p11a_solv", "m_ch2cl2", "solvent", e5, volume=("16.9 mL", quantity(16.9, "mL"))),
                      participant("p11a_silica", "m_silica", "purification_agent", e5, stage="purification"),
                      participant("p11a_prod", "m11a", "product", e5, product_role="main product")],
        steps=[step("p11a_add", e5, pids=["p11a_11", "p11a_imid", "p11a_tmscl", "p11a_solv"], operation="TMSCl added dropwise", temperature=("0 °C", quantity(0, "°C"))),
               step("p11a_warm_hold", e5, pids=["p11a_11", "p11a_imid", "p11a_tmscl"], operation="slowly warm to room temperature and stir", temperature=("rt", "rt"), time=("30 min", quantity(30, "min"))),
               step("p11a_workup", e5, stage="workup", pids=[], operation="filter precipitate, wash with CH2Cl2, and evaporate filtrate"),
               step("p11a_purify", e5, stage="purification", pids=["p11a_silica"], operation="silica column, hexane/EtOAc 55:1 to 30:1")],
        outcomes=[outcome("p11a_yield", "yield", ["m11a"], e5, raw="83%", normalized=quantity(83, "%"), yield_type="isolated", measurement="isolated after silica chromatography", scope="11 to 11a", product_role="main product"),
                  product_amount("p11a_amount", "m11a", e5, "3.14 g", quantity(3.14, "g"), "13.97 mmol", quantity(13.97, "mmol"))],
        class_changes=[chg("O", "Si", 0, 1, "formed"), chg("O", "H", 1, 0, "hydrogen_loss")],
        procedure="Full page-5 preparation of 11a.", solvent="CH2Cl2",
        concentration=("1.0 M", quantity(1.0, "M")), concentration_basis="compound 11",
        unresolved=["Partial protection mapping; HCl/imidazolium coproducts are not mapped."],
    ))
    reactions.append(reaction(
        "prep_7b", e5, categories=["other"], location={"page": 5, "region": "Preparation of 7b"},
        reaction_name="Preparation of (2R,3R)-3-chloro-2-((trimethylsilyl)oxy)butanal (7b)",
        transformation="Dibal-H reduction of a methyl ester to an aldehyde",
        participants=[participant("p7b_11a", "m11a", "reactant", e5, mass=("1.08 g", quantity(1.08, "g")), amount=("4.79 mmol", quantity(4.79, "mmol"))),
                      participant("p7b_dibal", "m_dibalh", "reagent", e5, volume=("11 mL of 1.1 M solution", quantity(11, "mL")), amount=("10 mmol", quantity(10, "mmol"))),
                      participant("p7b_cyclo", "m_cyclohexane", "solvent", e5),
                      participant("p7b_solv", "m_ch2cl2", "solvent", e5, volume=("29 mL", quantity(29, "mL"))),
                      participant("p7b_meoh", "m_meoh", "workup_agent", e5, stage="workup"),
                      participant("p7b_rochelle", "m_rochelle", "workup_agent", e5, stage="workup"),
                      participant("p7b_brine", "m_brine", "workup_agent", e5, stage="workup"),
                      participant("p7b_na2so4", "m_na2so4", "workup_agent", e5, stage="workup"),
                      participant("p7b_silica", "m_silica", "purification_agent", e5, stage="purification"),
                      participant("p7b_prod", "m7b", "product", e5, product_role="main product")],
        steps=[step("p7b_add", e5, pids=["p7b_11a", "p7b_dibal", "p7b_solv", "p7b_cyclo"], operation="Dibal-H solution added dropwise", temperature=("−78 °C", quantity(-78, "°C"))),
               step("p7b_hold", e5, pids=["p7b_11a", "p7b_dibal"], operation="reaction hold", temperature=("−78 °C", quantity(-78, "°C")), time=("30 min", quantity(30, "min"))),
               step("p7b_quench", e5, stage="workup", pids=["p7b_meoh", "p7b_rochelle"], operation="quench with MeOH at −78 °C, warm to rt, and treat with saturated potassium sodium tartrate", temperature=("−78 °C to rt", range_quantity(-78, 25, "°C", "−78 °C to rt"))),
               step("p7b_workup", e5, stage="workup", pids=["p7b_brine", "p7b_na2so4"], operation="filter; CH2Cl2 extraction twice; brine wash; Na2SO4 drying; vacuum concentration"),
               step("p7b_purify", e5, stage="purification", pids=["p7b_silica"], operation="column chromatography, hexane/EtOAc 40:1 to 25:1")],
        outcomes=[outcome("p7b_yield", "yield", ["m7b"], e5, raw="82%", normalized=quantity(82, "%"), yield_type="isolated", measurement="isolated after column chromatography", scope="11a to 7b", product_role="main product"),
                  product_amount("p7b_amount", "m7b", e5, "770 mg", quantity(770, "mg"), "3.95 mmol", quantity(3.95, "mmol"))],
        class_changes=[chg("C", "O", 1, 0, "broken"), chg("C", "H", 0, 1, "hydrogen_gain")],
        procedure="Full page-5 preparation of 7b.", solvent="CH2Cl2; Dibal-H carrier cyclohexane retained separately",
        concentration=("0.16 M", quantity(0.16, "M")), concentration_basis="compound 11a",
        unresolved=["Partial ester-to-aldehyde mapping; methoxy/aluminium coproducts are not mapped."],
    ))

    # Table 1, one row per experiment.
    table1 = [
        (1, -78, 60, -10, 3, "m_ch2cl2", "CH2Cl2", None, None),
        (2, -50, 60, -20, 4, "m_ch2cl2", "CH2Cl2", ("<10%", quantity(10, "%", "lt", "<10%")), None),
        (3, -50, 30, -20, 4, "m_ch2cl2", "CH2Cl2", 11, (">99%", quantity(99, "%", "gt", ">99%"))),
        (4, -50, 30, -10, 4, "m_ch2cl2", "CH2Cl2", 34, (">99%", quantity(99, "%", "gt", ">99%"))),
        (5, -50, 20, -10, 3, "m_ch2cl2", "CH2Cl2", 21, (">99%", quantity(99, "%", "gt", ">99%"))),
        (6, -30, 30, 0, 3, "m_ch2cl2", "CH2Cl2", 18, (">99%", quantity(99, "%", "gt", ">99%"))),
        (7, -50, 30, -10, 4, "m_toluene", "toluene", ("<5%", quantity(5, "%", "lt", "<5%")), None),
        (8, -50, 30, -10, 3, "m_ch2cl2", "CH2Cl2/PhCF3", ("<10%", quantity(10, "%", "lt", "<10%")), (">99%", quantity(99, "%", "gt", ">99%"))),
    ]
    for args in table1:
        reactions.append(table1_record(*args))

    # Table 2, one row per experiment.
    table2 = [
        (1, "m_ph3aucl", "Ph3AuCl", 10, "m_ch2cl2", "CH2Cl2", 30, 7, ("<5%", quantity(5, "%", "lt", "<5%"))),
        (2, "m_ph3aubf4", "Ph3AuBF4", 10, "m_ch2cl2", "CH2Cl2", 20, 5, 41),
        (3, "m_ph3auasf6", "Ph3AuAsF6", 10, "m_ch2cl2", "CH2Cl2", 20, 5, 34),
        (4, "m_ph3pauntf2", "Ph3PAuNTf2", 10, "m_ch2cl2", "CH2Cl2", 0, 1, 68),
        (5, "m_ph3pauntf2", "Ph3PAuNTf2", 10, "m_phcf3", "PhCF3", 20, 5, 55),
        (6, "m_ph3pauntf2", "Ph3PAuNTf2", 5, "m_ch2cl2", "CH2Cl2", 20, 7, 47),
    ]
    for args in table2:
        reactions.append(table2_record(*args))

    # 3 -> reported 1.
    e35 = [ev(3), ev(5), ev(3, "image"), ev(5, "image")]
    reactions.append(reaction(
        "prep_reported_hypoxylactone_1", e35, categories=["other"],
        location={"pages": [3, 5], "regions": ["Scheme 3 final arrow", "(+)-Hypoxylactone (1, Reported)"]},
        reaction_name="Preparation of reported hypoxylactone stereoisomer 1",
        transformation="intramolecular chlorohydrin epoxide formation with Ag2O",
        participants=[participant("p1_3", "m3", "reactant", e35, mass=("20 mg", quantity(20, "mg")), amount=("0.091 mmol", quantity(0.091, "mmol"))),
                      participant("p1_ag2o", "m_ag2o", "reagent", e35, mass=("21 mg", quantity(21, "mg")), amount=("0.091 mmol", quantity(0.091, "mmol")), equivalents=("1.0 equiv", quantity(1.0, "equiv")), reference="m3"),
                      participant("p1_thf", "m_thf", "solvent", e35, volume=("2.3 mL", quantity(2.3, "mL"))),
                      participant("p1_water", "m_water", "solvent", e35, volume=("0.6 mL", quantity(0.6, "mL"))),
                      participant("p1_celite", "m_celite", "workup_agent", e35, stage="workup"),
                      participant("p1_silica", "m_silica", "purification_agent", e35, stage="purification"),
                      participant("p1_prod", "m1", "product", e35, product_role="main product")],
        steps=[step("p1_add", e35, pids=["p1_3", "p1_ag2o", "p1_thf", "p1_water"], operation="Ag2O added to 3 in THF/H2O", temperature=("rt", "rt")),
               step("p1_heat", e35, pids=["p1_3", "p1_ag2o"], operation="stir in an oil bath", temperature=("70 °C", quantity(70, "°C")), time=("4 h", quantity(4, "h"))),
               step("p1_filter", e35, stage="workup", pids=["p1_celite"], operation="cool to rt, filter through a short Celite pad, and concentrate filtrate"),
               step("p1_purify", e35, stage="purification", pids=["p1_silica"], operation="silica chromatography, hexane/EtOAc 3:1")],
        outcomes=[outcome("p1_yield", "yield", ["m1"], e35, raw="79%", normalized=quantity(79, "%"), yield_type="isolated", measurement="isolated after silica chromatography", scope="3 to reported structure 1", product_role="main product"),
                  product_amount("p1_amount", "m1", e35, "13 mg", quantity(13, "mg"), "0.072 mmol", quantity(0.072, "mmol"))],
        class_changes=[chg("C", "O", 0, 1, "formed"), chg("C", "Cl", 1, 0, "broken"), chg("O", "H", 1, 0, "hydrogen_loss")],
        procedure="Full page-5 preparation of reported 1.", solvent="THF/H2O",
        solvent_ratio="2.3 mL/0.6 mL (reported component volumes)", total_volume=("2.9 mL from reported component volumes", quantity(2.9, "mL")),
        notes=[direct_note("The authors state that the NMR data of synthetic 1 did not match the reported natural-product values.", [ev(3), ev(3, "image")])],
        unresolved=["Product-1 candidate absolute stereodescriptors are structure-derived and require expert verification.",
                    "Unspecified unsuccessful basic/low-temperature epoxidation trials are source context, not reconstructed experiments."],
        general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    ))

    # ent-11 preparation by referenced analogous procedure.
    reactions.append(reaction(
        "prep_ent_11_by_analogy", e5, categories=["other"], location={"page": 5, "region": "opening sentence of ent-11a preparation"},
        reaction_name="Preparation of ent-11 from L-threonine according to the procedure for 11",
        transformation="analogous preparation of the enantiomeric chlorohydrin",
        participants=[participant("pent11_thr", "m_threonine_source_conflict", "reactant", e5),
                      participant("pent11_prod", "m_ent11", "product", e5, product_role="intermediate used without reported isolated yield")],
        steps=[step("pent11_analogous", e5, pids=["pent11_thr", "pent11_prod"], operation="prepared from L-threonine according to the procedure for 11")],
        outcomes=[], class_changes=(), class_coverage="unresolved",
        procedure="Compound ent-11 was prepared from L-threonine according to the procedure for 11.",
        unresolved=["No ent-11 charge, isolated amount, yield, or explicit intermediate sequence is reported in pages 1–6; these were not borrowed from the opposite enantiomer."],
    ))

    # ent-11 -> ent-11a and ent-11a -> 8a.
    reactions.append(reaction(
        "prep_ent_11a", e5, categories=["other"], location={"page": 5, "region": "Preparation of ent-11a"},
        reaction_name="Preparation of (2S,3S)-methyl 2-((tert-butyldimethylsilyl)oxy)-3-chlorobutanoate (ent-11a)",
        transformation="tert-butyldimethylsilyl protection of an alcohol",
        participants=[participant("pe11a_ent11", "m_ent11", "reactant", e5, mass=("600 mg", quantity(600, "mg")), amount=("3.95 mmol", quantity(3.95, "mmol"))),
                      participant("pe11a_imid", "m_imidazole", "reagent", e5, mass=("538 mg", quantity(538, "mg")), amount=("7.9 mmol", quantity(7.9, "mmol"))),
                      participant("pe11a_tbscl", "m_tbscl", "reagent", e5, mass=("1.11 g", quantity(1.11, "g")), amount=("7.11 mmol", quantity(7.11, "mmol"))),
                      participant("pe11a_dmf", "m_dmf", "solvent", e5, volume=("10 mL", quantity(10, "mL"))),
                      participant("pe11a_brine", "m_brine", "workup_agent", e5, stage="workup"),
                      participant("pe11a_etoac", "m_etoac", "workup_agent", e5, stage="workup"),
                      participant("pe11a_mgso4", "m_mgso4", "workup_agent", e5, stage="workup"),
                      participant("pe11a_silica", "m_silica", "purification_agent", e5, stage="purification"),
                      participant("pe11a_prod", "m_ent11a", "product", e5, product_role="main product")],
        steps=[step("pe11a_hold", e5, pids=["pe11a_ent11", "pe11a_imid", "pe11a_tbscl", "pe11a_dmf"], operation="stir", temperature=("rt", "rt"), time=("15 h", quantity(15, "h"))),
               step("pe11a_workup", e5, stage="workup", pids=["pe11a_brine", "pe11a_etoac", "pe11a_mgso4"], operation="dilute with brine; extract with EtOAc; MgSO4 dry; evaporate"),
               step("pe11a_purify", e5, stage="purification", pids=["pe11a_silica"], operation="silica chromatography, hexane/EtOAc 40:1 to 20:1")],
        outcomes=[outcome("pe11a_yield", "yield", ["m_ent11a"], e5, raw="99%", normalized=quantity(99, "%"), yield_type="isolated", measurement="isolated after silica chromatography", scope="ent-11 to ent-11a", product_role="main product"),
                  product_amount("pe11a_amount", "m_ent11a", e5, "1.05 g", quantity(1.05, "g"), "3.93 mmol", quantity(3.93, "mmol"))],
        class_changes=[chg("O", "Si", 0, 1, "formed"), chg("O", "H", 1, 0, "hydrogen_loss")],
        procedure="Full page-5 preparation of ent-11a.", solvent="DMF",
        concentration=("0.4 M", quantity(0.4, "M")), concentration_basis="ent-11",
        unresolved=["Partial silylation mapping; coproducts are not mapped."],
    ))
    e56 = [ev(5), ev(6), ev(5, "image"), ev(6, "image")]
    reactions.append(reaction(
        "prep_8a", e56, categories=["other"], location={"pages": [5, 6], "region": "Preparation of 8a, continued"},
        reaction_name="Preparation of (2S,3S)-2-((tert-butyldimethylsilyl)oxy)-3-chlorobutanal (8a)",
        transformation="Dibal-H reduction of a methyl ester to an aldehyde",
        participants=[participant("p8a_ent11a", "m_ent11a", "reactant", e56, mass=("1.05 g", quantity(1.05, "g")), amount=("3.93 mmol", quantity(3.93, "mmol"))),
                      participant("p8a_dibal", "m_dibalh", "reagent", e56, volume=("7.4 mL of 1.1 M solution", quantity(7.4, "mL")), amount=("8.17 mmol", quantity(8.17, "mmol"))),
                      participant("p8a_cyclo", "m_cyclohexane", "solvent", e56),
                      participant("p8a_ch2cl2", "m_ch2cl2", "solvent", e56, volume=("25 mL", quantity(25, "mL"))),
                      participant("p8a_meoh", "m_meoh", "workup_agent", e56, stage="workup"),
                      participant("p8a_rochelle", "m_rochelle", "workup_agent", e56, stage="workup"),
                      participant("p8a_brine", "m_brine", "workup_agent", e56, stage="workup"),
                      participant("p8a_na2so4", "m_na2so4", "workup_agent", e56, stage="workup"),
                      participant("p8a_silica", "m_silica", "purification_agent", e56, stage="purification"),
                      participant("p8a_prod", "m8a", "product", e56, product_role="main product")],
        steps=[step("p8a_add", e56, pids=["p8a_ent11a", "p8a_dibal", "p8a_ch2cl2", "p8a_cyclo"], operation="Dibal-H solution added dropwise", temperature=("−78 °C", quantity(-78, "°C"))),
               step("p8a_hold", e56, pids=["p8a_ent11a", "p8a_dibal"], operation="reaction hold", temperature=("−78 °C", quantity(-78, "°C")), time=("30 min", quantity(30, "min"))),
               step("p8a_quench", e56, stage="workup", pids=["p8a_meoh", "p8a_rochelle"], operation="quench with MeOH at −78 °C, warm to rt, and treat with saturated potassium sodium tartrate", temperature=("−78 °C to rt", range_quantity(-78, 25, "°C", "−78 °C to rt"))),
               step("p8a_workup", e56, stage="workup", pids=["p8a_brine", "p8a_na2so4"], operation="filter; CH2Cl2 extraction twice; brine wash; Na2SO4 drying; vacuum concentration"),
               step("p8a_purify", e56, stage="purification", pids=["p8a_silica"], operation="column chromatography, hexane/EtOAc 40:1 to 25:1")],
        outcomes=[outcome("p8a_yield", "yield", ["m8a"], e56, raw="84%", normalized=quantity(84, "%"), yield_type="isolated", measurement="isolated after column chromatography", scope="ent-11a to 8a", product_role="main product"),
                  product_amount("p8a_amount", "m8a", e56, "781 mg", quantity(781, "mg"), "3.30 mmol", quantity(3.30, "mmol"))],
        class_changes=[chg("C", "O", 1, 0, "broken"), chg("C", "H", 0, 1, "hydrogen_gain")],
        procedure="Full page-5 preparation of 8a, with characterization continuing on page 6.",
        solvent="CH2Cl2; Dibal-H carrier cyclohexane retained separately",
        concentration=("0.16 M", quantity(0.16, "M")), concentration_basis="ent-11a",
        unresolved=["Partial ester-to-aldehyde mapping; methoxy/aluminium coproducts are not mapped."],
    ))

    # 6 + 8a -> 13.
    e346 = [ev(3), ev(4), ev(6), ev(3, "image"), ev(4, "image"), ev(6, "image")]
    p13 = [
        participant("p13_6", "m6", "reactant", e346, mass=("0.16 g", quantity(0.16, "g")), amount=("1.02 mmol", quantity(1.02, "mmol")), equivalents=("1.06 equiv", quantity(1.06, "equiv")), reference="m5"),
        participant("p13_8a", "m8a", "reactant", e346, mass=("237 mg", quantity(237, "mg")), amount=("1.0 mmol", quantity(1.0, "mmol")), equivalents=("1.04 equiv", quantity(1.04, "equiv")), reference="m5"),
        participant("p13_5", "m5", "reagent", e346, equivalents=("1.0 equiv", quantity(1.0, "equiv")), reference="m5"),
        participant("p13_ligand", "m_ligand", "reagent", e346, mass=("0.50 g", quantity(0.50, "g")), amount=("0.96 mmol", quantity(0.96, "mmol"))),
        participant("p13_bbr3", "m_bbr3", "reagent", e346, volume=("1.0 mL of 1 M solution", quantity(1.0, "mL")), amount=("1.0 mmol", quantity(1.0, "mmol"))),
        participant("p13_base", "m_ipr2net", "reagent", e346, volume=("0.35 mL", quantity(0.35, "mL")), mass=("0.26 g", quantity(0.26, "g")), amount=("2.01 mmol", quantity(2.01, "mmol")), equivalents=("2.1 equiv", quantity(2.1, "equiv")), reference="m5"),
        participant("p13_ch2cl2", "m_ch2cl2", "solvent", e346, volume=("7 mL initial CH2Cl2", quantity(7, "mL"))),
        participant("p13_ch2cl2_exchange", "m_ch2cl2", "solvent", e346, volume=("5 mL solvent exchange", quantity(5, "mL"))),
        participant("p13_ch2cl2_redissolve", "m_ch2cl2", "solvent", e346, volume=("3 mL redissolution", quantity(3, "mL"))),
        participant("p13_ch2cl2_6carrier", "m_ch2cl2", "solvent", e346, volume=("1 mL carrier for 6/i-Pr2NEt", quantity(1, "mL"))),
        participant("p13_ch2cl2_8acarrier", "m_ch2cl2", "solvent", e346, volume=("1 mL carrier for 8a", quantity(1, "mL"))),
        participant("p13_buffer", "m_buffer7", "workup_agent", e346, stage="workup", volume=("7 mL", quantity(7, "mL"))),
        participant("p13_ch2cl2_quench", "m_ch2cl2", "workup_agent", e346, stage="workup", volume=("ca. 10 mL added after buffer", quantity(10, "mL", "approx", "ca. 10 mL"))),
        participant("p13_ch2cl2_extract1", "m_ch2cl2", "workup_agent", e346, stage="workup", volume=("ca. 10 mL first extraction", quantity(10, "mL", "approx", "ca. 10 mL"))),
        participant("p13_ch2cl2_extract2", "m_ch2cl2", "workup_agent", e346, stage="workup", volume=("ca. 10 mL second extraction", quantity(10, "mL", "approx", "ca. 10 mL"))),
        participant("p13_brine", "m_brine", "workup_agent", e346, stage="workup"),
        participant("p13_na2so4", "m_na2so4", "workup_agent", e346, stage="workup"),
        participant("p13_ether", "m_ether", "workup_agent", e346, stage="workup", volume=("ca. 20 mL", quantity(20, "mL", "approx", "ca. 20 mL"))),
        participant("p13_silica", "m_silica", "purification_agent", e346, stage="purification"),
        participant("p13_prod", "m13", "product", e346, product_role="main product"),
    ]
    reactions.append(reaction(
        "prep_13", e346, categories=["other"], location={"pages": [3, 4, 6], "regions": ["matched-case prose", "Scheme 4 first arrow", "Preparation of 13"]},
        reaction_name="Matched allenoate gamma-addition of 6 to 8a",
        transformation="asymmetric allenoate gamma-addition to an aldehyde",
        participants=p13,
        steps=[step("p13_make5_0", e346, pids=["p13_ligand", "p13_bbr3", "p13_ch2cl2"], operation="BBr3 added to ligand in 7 mL CH2Cl2 and stirred", temperature=("0 °C", quantity(0, "°C")), time=("30 min", quantity(30, "min")), atmosphere="N2"),
               step("p13_make5_25", e346, pids=["p13_ligand", "p13_bbr3"], operation="warm and stir", temperature=("25 °C", quantity(25, "°C")), time=("1.5 h", quantity(1.5, "h")), atmosphere="N2"),
               step("p13_vacuum", e346, pids=["p13_ligand", "p13_bbr3", "p13_ch2cl2_exchange"], operation="concentrate at 1 mmHg; add 5 mL CH2Cl2 and evaporate again"),
               step("p13_enolize", e346, pids=["p13_5", "p13_6", "p13_base", "p13_ch2cl2_redissolve", "p13_ch2cl2_6carrier"], operation="redissolve 5 in 3 mL CH2Cl2; add 6/i-Pr2NEt in 1 mL CH2Cl2 dropwise and stir", temperature=("−50 °C", quantity(-50, "°C")), time=("30 min", quantity(30, "min")), atmosphere="N2"),
               step("p13_add8a", e346, pids=["p13_8a", "p13_ch2cl2_8acarrier"], operation="precooled 8a in 1 mL CH2Cl2 added by cannula while maintaining temperature", temperature=("−50 °C", quantity(-50, "°C"))),
               step("p13_warm_hold", e346, pids=["p13_5", "p13_6", "p13_8a", "p13_base"], operation="allow to warm from −50 °C to −20 °C and proceed", temperature=("−50 to −20 °C", range_quantity(-50, -20, "°C", "−50 to −20 °C")), time=("2 h", quantity(2, "h")), atmosphere="N2"),
               step("p13_workup", e346, stage="workup", pids=["p13_buffer", "p13_ch2cl2_quench", "p13_ch2cl2_extract1", "p13_ch2cl2_extract2", "p13_brine", "p13_na2so4", "p13_ether"], operation="pH 7 buffer quench; add ca. 10 mL CH2Cl2; two ca. 10 mL CH2Cl2 extractions; brine wash; Na2SO4 dry; ether precipitation; filter recovered bis-sulfonamide", temperature=("0 °C for precipitation", quantity(0, "°C")), time=("20 min", quantity(20, "min"))),
               step("p13_purify", e346, stage="purification", pids=["p13_silica"], operation="column chromatography, hexane/EtOAc 30:1")],
        outcomes=[outcome("p13_yield", "yield", ["m13"], e346, raw="75%", normalized=quantity(75, "%"), yield_type="isolated", measurement="isolated after column chromatography", scope="6 + 8a to 13", product_role="main product"),
                  product_amount("p13_amount", "m13", e346, "281 mg", quantity(281, "mg"), "0.72 mmol", quantity(0.72, "mmol")),
                  outcome("p13_single_isomer", "selectivity", ["m13"], e346, raw="as a single isomer", normalized="single regio- and stereoisomer", scope={"metric": "qualitative regio-/stereoselectivity", "major_isomer": "13", "minor_isomer": "not reported"}, product_role="main product")],
        class_changes=[chg("C", "C", 0, 1, "formed"), chg("C", "O", 2, 1, "order_decreased"), chg("O", "H", 0, 1, "hydrogen_gain")],
        procedure="Full page-6 preparation of 13; Scheme 4 and page-3 prose are duplicate summaries.",
        notes=[direct_note("The authors describe this as a matched, anti-Felkin–Anh case.", [ev(3), ev(4)])],
        unresolved=["Candidate allene axial configuration is retained textually but not encoded by the toolkit SMILES.", "Complete atom mapping and qualitative selectivity measurement method are unreported."],
        solvent="CH2Cl2", general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    ))

    # 13 -> 14.
    e346b = [ev(3), ev(4), ev(6), ev(3, "image"), ev(4, "image"), ev(6, "image")]
    reactions.append(reaction(
        "prep_14", e346b, categories=["other"], location={"pages": [3, 4, 6], "regions": ["cyclization prose", "Scheme 4 second arrow", "Preparation of 14"]},
        reaction_name="Gold-catalyzed cyclization of 13 to lactone 14",
        transformation="gold-catalyzed allenoate oxacyclization",
        participants=[participant("p14_13", "m13", "reactant", e346b, mass=("74.3 mg", quantity(74.3, "mg")), amount=("0.19 mmol", quantity(0.19, "mmol"))),
                      participant("p14_cat", "m_ph3pauntf2", "catalyst", e346b, mass=("30 mg", quantity(30, "mg")), amount=("0.019 mmol", quantity(0.019, "mmol")), mol_percent=("10 mol%", quantity(10, "mol%")), reference="m13"),
                      participant("p14_solv", "m_ch2cl2", "solvent", e346b, volume=("2 mL", quantity(2, "mL"))),
                      participant("p14_silica", "m_silica", "purification_agent", e346b, stage="purification"),
                      participant("p14_prod", "m14", "product", e346b, product_role="main product")],
        steps=[step("p14_hold", e346b, pids=["p14_13", "p14_cat", "p14_solv"], operation="cyclization", temperature=("0 °C", quantity(0, "°C")), time=("1 h", quantity(1, "h"))),
               step("p14_concentrate", e346b, stage="workup", pids=[], operation="concentrate under reduced pressure"),
               step("p14_purify", e346b, stage="purification", pids=["p14_silica"], operation="column chromatography, hexane/EtOAc 3:1 to 1:1")],
        outcomes=[outcome("p14_yield", "yield", ["m14"], e346b, raw="78%", normalized=quantity(78, "%"), yield_type="isolated", measurement="isolated after column chromatography", scope="13 to 14", product_role="main product"),
                  product_amount("p14_amount", "m14", e346b, "49.6 mg", quantity(49.6, "mg"), "0.15 mmol", quantity(0.15, "mmol"))],
        class_changes=[chg("C", "O", 0, 1, "formed")],
        procedure="Full page-6 preparation of 14; Scheme 4 and page-3 prose are duplicate summaries.", solvent="CH2Cl2",
        total_volume=("2 mL", quantity(2, "mL")), concentration=("0.095 M from reported charge and solvent volume", quantity(0.095, "M")),
        concentration_basis="inferred from 0.19 mmol 13 in 2 mL CH2Cl2",
        unresolved=["Partial ring-formation inventory; complete atom mapping, tert-butyl-derived coproducts, and all bond-order shifts are unresolved."],
        general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    ))

    # 14 -> 4.
    e46 = [ev(4), ev(6), ev(4, "image"), ev(6, "image")]
    reactions.append(reaction(
        "prep_xylogiblactone_4", e46, categories=["other"], location={"pages": [4, 6], "regions": ["Scheme 4 deprotection", "(+)-Xylogiblactone A (4)"]},
        reaction_name="Preparation of (+)-xylogiblactone A (4)",
        transformation="acidic TBS deprotection",
        participants=[participant("p4_14", "m14", "reactant", e46, mass=("30.6 mg", quantity(30.6, "mg")), amount=("0.091 mmol", quantity(0.091, "mmol"))),
                      participant("p4_hcl", "m_hcl", "reagent", e46, volume=("0.16 mL", quantity(0.16, "mL")), amount=("0.50 mmol", quantity(0.50, "mmol")), equivalents=("5.5 equiv", quantity(5.5, "equiv")), reference="m14"),
                      participant("p4_meoh", "m_meoh", "solvent", e46, volume=("1 mL", quantity(1, "mL"))),
                      participant("p4_nahco3", "m_nahco3", "workup_agent", e46, stage="workup"),
                      participant("p4_brine", "m_brine", "workup_agent", e46, stage="workup"),
                      participant("p4_na2so4", "m_na2so4", "workup_agent", e46, stage="workup"),
                      participant("p4_silica", "m_silica", "purification_agent", e46, stage="purification"),
                      participant("p4_prod", "m4", "product", e46, product_role="main product")],
        steps=[step("p4_hold", e46, pids=["p4_14", "p4_hcl", "p4_meoh"], operation="stir", temperature=("rt", "rt"), time=("overnight", "overnight")),
               step("p4_workup", e46, stage="workup", pids=["p4_nahco3", "p4_brine", "p4_na2so4"], operation="quench with saturated aqueous NaHCO3; CH2Cl2 extraction; brine wash; Na2SO4 dry; concentrate"),
               step("p4_purify", e46, stage="purification", pids=["p4_silica"], operation="flash chromatography, hexane/EtOAc 3:1")],
        outcomes=[outcome("p4_yield", "yield", ["m4"], e46, raw="90%", normalized=quantity(90, "%"), yield_type="isolated", measurement="isolated after flash chromatography", scope="14 to 4", product_role="main product"),
                  product_amount("p4_amount", "m4", e46, "18 mg", quantity(18, "mg"), "0.082 mmol", quantity(0.082, "mmol"))],
        class_changes=[chg("O", "Si", 1, 0, "broken"), chg("O", "H", 0, 1, "hydrogen_gain")],
        procedure="Full page-6 preparation of 4; page-4 prose and Scheme 4 summarize the same experiment.", solvent="MeOH",
        notes=[direct_note("The authors report that basic fluoride caused C4 epimerization in preliminary attempts; those unspecified trials are not reconstructed as concrete experiments.", [ev(3), ev(4)])],
        unresolved=["Partial deprotection inventory; silicon-containing coproducts are not mapped."],
        general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    ))

    # 4 -> corrected 9.
    e46b = [ev(4), ev(6), ev(4, "image"), ev(6, "image")]
    reactions.append(reaction(
        "prep_corrected_hypoxylactone_9", e46b, categories=["other"], location={"pages": [4, 6], "regions": ["Scheme 4 final arrow", "(+)-Hypoxylactone (9, Corrected)"]},
        reaction_name="Preparation of corrected (+)-hypoxylactone (9)",
        transformation="intramolecular chlorohydrin epoxide formation with Ag2O",
        participants=[participant("p9_4", "m4", "reactant", e46b, mass=("20 mg", quantity(20, "mg")), amount=("0.091 mmol", quantity(0.091, "mmol"))),
                      participant("p9_ag2o", "m_ag2o", "reagent", e46b, mass=("21 mg", quantity(21, "mg")), amount=("0.091 mmol", quantity(0.091, "mmol")), equivalents=("1.0 equiv", quantity(1.0, "equiv")), reference="m4"),
                      participant("p9_thf", "m_thf", "solvent", e46b, volume=("2.3 mL", quantity(2.3, "mL"))),
                      participant("p9_water", "m_water", "solvent", e46b, volume=("0.6 mL", quantity(0.6, "mL"))),
                      participant("p9_celite", "m_celite", "workup_agent", e46b, stage="workup"),
                      participant("p9_silica", "m_silica", "purification_agent", e46b, stage="purification"),
                      participant("p9_prod", "m9", "product", e46b, product_role="main product")],
        steps=[step("p9_add", e46b, pids=["p9_4", "p9_ag2o", "p9_thf", "p9_water"], operation="Ag2O added to 4 in THF/H2O", temperature=("rt", "rt")),
               step("p9_heat", e46b, pids=["p9_4", "p9_ag2o"], operation="stir in an oil bath", temperature=("70 °C", quantity(70, "°C")), time=("4 h", quantity(4, "h"))),
               step("p9_filter", e46b, stage="workup", pids=["p9_celite"], operation="cool to rt, filter through a short Celite pad, and concentrate filtrate"),
               step("p9_purify", e46b, stage="purification", pids=["p9_silica"], operation="silica chromatography, hexane/EtOAc 3:1")],
        outcomes=[outcome("p9_yield", "yield", ["m9"], e46b, raw="82%", normalized=quantity(82, "%"), yield_type="isolated", measurement="isolated after silica chromatography", scope="4 to corrected 9", product_role="main product"),
                  product_amount("p9_amount", "m9", e46b, "13.7 mg", quantity(13.7, "mg"), "0.075 mmol", quantity(0.075, "mmol"))],
        class_changes=[chg("C", "O", 0, 1, "formed"), chg("C", "Cl", 1, 0, "broken"), chg("O", "H", 1, 0, "hydrogen_loss")],
        procedure="Full page-6 preparation of corrected 9; Scheme 4 and page-4 prose summarize the same experiment.",
        solvent="THF/H2O", solvent_ratio="2.3 mL/0.6 mL (reported component volumes)", total_volume=("2.9 mL from reported component volumes", quantity(2.9, "mL")),
        notes=[direct_note("The authors state that all physical and spectroscopic data of 9 were identical with the original natural-product data.", [ev(4), ev(6)])],
        unresolved=["Product-9 candidate absolute stereodescriptors are structure-derived and require expert verification."],
        general_scheme_ids=[f"{DOC}:scheme1_synthetic_strategy"],
    ))
    return reactions


def source_context():
    return [
        {"context_id": f"{DOC}:context_page1_overview", "document_id": DOC,
         "kind": "article_overview", "evidence_ids": [ev(1), ev(1, "image")],
         "note": "Title, abstract, graphical overview, and Figure 1 natural-product relationship. Overview arrows duplicate later concrete records and are not counted as new experiments."},
        {"context_id": f"{DOC}:context_early_epoxide_opening_trials", "document_id": DOC,
         "kind": "aggregated_optimization_discussion", "evidence_ids": [ev(2), ev(2, "image")],
         "note": "Early TiCl4 or SnCl4 trials under various unspecified conditions gave 11 in 13% with impurities. Conditions cannot be separated into concrete rows, so no fabricated experiment was created."},
        {"context_id": f"{DOC}:context_unspecified_epoxidation_trials", "document_id": DOC,
         "kind": "aggregated_optimization_discussion", "evidence_ids": [ev(3), ev(3, "image")],
         "note": "Various basic and low-temperature epoxidation trials caused mixtures/epimerization or poor conversion; no unique conditions are reported."},
        {"context_id": f"{DOC}:context_xray_3", "document_id": DOC,
         "kind": "structure_characterization", "evidence_ids": [ev(3), ev(3, "image")],
         "note": "Figure 2 reports the X-ray structure of 3 (50% ellipsoid probability)."},
        {"context_id": f"{DOC}:context_characterization_pages4_6", "document_id": DOC,
         "kind": "characterization", "evidence_ids": [ev(4), ev(5), ev(6), ev(4, "image"), ev(5, "image"), ev(6, "image")],
         "note": "NMR, IR, HRMS/LC-MS/GC-MS, optical rotation, and Rf data support compound characterization and are not separate synthesis records."},
        {"context_id": f"{DOC}:context_si_availability", "document_id": DOC,
         "kind": "supplement_availability_statement", "evidence_ids": [ev(6), ev(6, "image")],
         "note": "Page 6 says Supporting Information is available online, but no SI file was supplied or explicitly associated in this workspace; no SI evidence was imported."},
    ]


def raw_inventory(reactions):
    ids = [r["reaction_id"] for r in reactions]
    return {
        "batch_id": "allene-stereo-v31-200-20261017",
        "document_id": DOC,
        "page_range": [1, 6],
        "semantic_accuracy_verified": False,
        "visual_review": {str(p): ev(p, "image") for p in range(1, 7)},
        "decisions": [
            "Same experiment repeated in prose/scheme/table is one reaction with multiple evidence IDs.",
            "Table 1 has eight independent optimization records; Table 2 has six.",
            "Table 1 entry 4 is merged with the detailed preparation of 12; Table 2 entry 4 is merged with the detailed preparation of 3.",
            "Unspecified groups of failed trials are source_context, not invented experiments.",
            "Scheme 2 threonine identity and KOH-stage temperature conflicts are preserved without choosing a side.",
            "Independent spectra and mass-spectral peak lists are characterization context, not reactions.",
            "All structures and bond classifications are machine candidates requiring review.",
        ],
        "reaction_ids": ids,
        "counts": {
            "all_records": len(reactions),
            "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
            "optimization_entries": sum("optimization" in (st.value(r["experiment_categories"]) or []) for r in reactions),
            "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in reactions),
            "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in reactions),
            "unresolved_summaries": sum(r["record_type"] == "unresolved" for r in reactions),
        },
    }


def coverage(reaction_ids):
    rid = lambda suffix: f"{DOC}:{suffix}"
    return {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title/abstract and graphical overview", "finding": "Article overview with duplicated summary arrows for mismatched/matched gamma-addition, gold cyclization, and Ag2O epoxidation; no separate experiment created.", "evidence_ids": [ev(1, "image"), ev(1)]},
                {"region": "Figure 1 and introduction", "finding": "Natural-product structures and proposed biological/chemical relationship; contextual, not a performed synthetic experiment.", "evidence_ids": [ev(1, "image"), ev(1)]}],
             "record_ids": [rid("scheme1_synthetic_strategy"), rid("table1_entry_4"), rid("table2_entry_4"), rid("prep_reported_hypoxylactone_1"), rid("prep_13"), rid("prep_14"), rid("prep_xylogiblactone_4"), rid("prep_corrected_hypoxylactone_9")],
             "unresolved": ["Graphical-abstract structures are summary occurrences; stereochemical candidates remain needs_review."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 1", "finding": "Conceptual synthetic strategy, recorded as a generic scheme rather than a performed experiment.", "evidence_ids": [ev(2, "image"), ev(2)]},
                {"region": "Scheme 2", "finding": "Preparation sequence for 10/11/7a/7b; explicit D/L-threonine and −78/−30 °C conflicts are preserved with page 4.", "evidence_ids": [ev(2, "image"), ev(2)]},
                {"region": "right-column experimental discussion", "finding": "Concrete failed 6+7a attempt and selected Table-1-entry-4 conditions; early TiCl4/SnCl4 trials remain aggregated context because unique rows are not reported.", "evidence_ids": [ev(2, "image"), ev(2)]}],
             "record_ids": [rid("scheme1_synthetic_strategy"), rid("scheme2_silylation_reduction_summary"), rid("failed_gamma_addition_6_7a"), rid("prep_10"), rid("prep_11"), rid("table1_entry_4")],
             "unresolved": ["Scheme 2's 75–83% is not allocated to 7a versus 7b.", "Scheme 2 lists HBr but the detailed procedure does not."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Scheme 3 and models B/C", "finding": "Selected gamma-addition/cyclization/epoxidation sequence plus proposed mismatched stereochemical models; duplicated concrete experiments merged with table/procedure records.", "evidence_ids": [ev(3, "image"), ev(3)]},
                {"region": "Table 1", "finding": "All eight gamma-addition optimization rows extracted as distinct experiments with yields and ee values/NA preserved.", "evidence_ids": [ev(3, "image"), ev(3)]},
                {"region": "Table 2", "finding": "All six gold-catalysis optimization rows extracted as distinct experiments; in-situ catalyst footnotes represented.", "evidence_ids": [ev(3, "image"), ev(3)]},
                {"region": "Figure 2 and lower-right prose", "finding": "X-ray characterization of 3 and unspecified epoxidation trials retained as source context; matched 8a→13 and 13→14 summaries linked to exact page-6 records.", "evidence_ids": [ev(3, "image"), ev(3)]}],
             "record_ids": [rid("scheme3_4_stereochemical_models")] + [rid(f"table1_entry_{i}") for i in range(1, 9)] + [rid(f"table2_entry_{i}") for i in range(1, 7)] + [rid("prep_reported_hypoxylactone_1"), rid("prep_13"), rid("prep_14")],
             "unresolved": ["Models B/C are proposed, not rigorously established.", "Conditions for the unsuccessful epoxidation groups are not uniquely reported."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "Scheme 4", "finding": "Four concrete arrows (6+8a→13, 13→14, 14→4, 4→9) linked to detailed records; model D is proposed mechanism.", "evidence_ids": [ev(4, "image"), ev(4)]},
                {"region": "Experimental general section", "finding": "General handling, chromatography, and analytical instrumentation retained as context, not separate reactions.", "evidence_ids": [ev(4, "image"), ev(4)]},
                {"region": "Preparations of 10 and 11", "finding": "Complete multistage operations, charges, workups, isolated amounts and yields extracted; Scheme-2 conflicts preserved.", "evidence_ids": [ev(4, "image"), ev(4)]}],
             "record_ids": [rid("scheme3_4_stereochemical_models"), rid("prep_10"), rid("prep_11"), rid("prep_13"), rid("prep_14"), rid("prep_xylogiblactone_4"), rid("prep_corrected_hypoxylactone_9")],
             "unresolved": ["General section does not establish missing conditions for any specific experiment."]},
            {"page": 5, "state": "attempted", "regions": [
                {"region": "Preparations of 11a and 7b", "finding": "Two separate experiments with ordered reaction/workup/purification steps and isolated amounts/yields.", "evidence_ids": [ev(5, "image"), ev(5)]},
                {"region": "Preparation of 12", "finding": "Exact procedure merged with Table 1 entry 4; ligand/BBr3 preformation, low-temperature addition/ramp, workup, 34% yield and >99% ee retained.", "evidence_ids": [ev(5, "image"), ev(5)]},
                {"region": "Preparation of 3 and reported 1", "finding": "Preparation of 3 merged with Table 2 entry 4; 3→1 Ag2O epoxidation recorded separately.", "evidence_ids": [ev(5, "image"), ev(5)]},
                {"region": "ent-11/ent-11a/8a section", "finding": "ent-11 analogous preparation kept partial; ent-11a experiment complete for supplied text; 8a procedure continues onto page 6.", "evidence_ids": [ev(5, "image"), ev(5)]}],
             "record_ids": [rid("prep_11a"), rid("prep_7b"), rid("table1_entry_4"), rid("table2_entry_4"), rid("prep_reported_hypoxylactone_1"), rid("prep_ent_11_by_analogy"), rid("prep_ent_11a"), rid("prep_8a")],
             "unresolved": ["ent-11 charge/yield and explicit intermediate sequence are not reported.", "Mass-spectral/NMR data are characterization, not reactions."]},
            {"page": 6, "state": "attempted", "regions": [
                {"region": "8a characterization continuation", "finding": "Continuation of characterization only; no duplicate synthesis created.", "evidence_ids": [ev(6, "image"), ev(6)]},
                {"region": "Preparation of 13", "finding": "Full matched gamma-addition operation, charges, −50→−20 °C program, workup, isolated amount, 75% yield and single-isomer statement extracted.", "evidence_ids": [ev(6, "image"), ev(6)]},
                {"region": "Preparations of 14, 4, and 9", "finding": "Three separate concrete experiments with full supplied conditions, workups, isolated amounts and yields.", "evidence_ids": [ev(6, "image"), ev(6)]},
                {"region": "Associated Content", "finding": "SI availability stated, but no SI file was supplied or explicitly associated; no SI facts imported.", "evidence_ids": [ev(6, "image"), ev(6)]}],
             "record_ids": [rid("prep_8a"), rid("prep_13"), rid("prep_14"), rid("prep_xylogiblactone_4"), rid("prep_corrected_hypoxylactone_9")],
             "unresolved": ["SI association remains unestablished.", "Candidate product stereochemistry and partial bond inventories require expert review."]},
        ],
    }


def workflow_text(bundle, raw):
    return f"""# Workflow log — allene-stereo-v31-200-20261017

## Source and bounded scope

- Document: `{DOC}`; DOI `10.1021/acs.joc.0c02194`.
- Title: *Synthesis of (+)-Hypoxylactone through Allenoate γ-Addition: Revision of Stereochemistry*.
- `runs/current/source_mapping.json` identifies the original article (`presentation_conversion=false`). SI association is `unestablished`; page 6 mentions online SI, but no SI file was supplied or imported.
- Original rendered PDF images for pages **1, 2, 3, 4, 5, and 6** were each opened at original detail. No page outside 1–6 is claimed as visually covered.
- The prior `runs/workflow.md` described a different DOI/document and was treated as stale workspace contamination, never as chemical evidence. This log replaces it.

## Page/record decisions

- Page 1: title, abstract, graphical overview, and Figure 1. The overview duplicates later arrows; natural-product/biological relations are source context.
- Page 2: Scheme 1 is a `generic_scheme`; Scheme 2's first two arrows are linked to exact preparations of 10/11, while its final 11→7a/7b range is an unresolved summary. The concrete 6+7a failure is a separate control.
- Page 3: all **8 Table 1 rows** and **6 Table 2 rows** are independent optimization records. Table 1 entry 4 is merged with the page-5 preparation of 12; Table 2 entry 4 is merged with the page-5 preparation of 3. Models B/C are proposed mechanism. Figure 2 is characterization.
- Page 4: all four Scheme 4 arrows map to separate concrete records. Model D is proposed mechanism. The full preparations of 10/11 begin here.
- Page 5: preparations of 11a, 7b, 12, 3, reported 1, ent-11a, and the start of 8a were extracted; the ent-11-by-analogy sentence is partial rather than silently borrowing the opposite-enantiomer yield.
- Page 6: the 8a characterization continuation is not a reaction; full procedures for 13, 14, 4, and corrected 9 are extracted. Associated Content does not establish a supplied SI association.
- Exact per-page regions, record IDs, and unresolved items are in `runs/current/coverage/pages_0001_0006.json`.

## v3.1 construction

- Ran `standard init` first and preserved the untouched template at `runs/current/raw_model/pages_0001_0006_template.json`.
- Built **{len(bundle['molecules'])} molecule/material records** and **{len(bundle['reactions'])} reaction records**: {raw['counts']['specific_experiments']} specific experiments, {raw['counts']['optimization_entries']} optimization rows, {raw['counts']['generic_schemes']} generic scheme, {raw['counts']['proposed_mechanisms']} proposed mechanism, and {raw['counts']['unresolved_summaries']} unresolved scheme summary.
- Every machine reaction remains `needs_review` and `partial`; `semantic_accuracy_verified=false` and `chemical_accuracy_verified=false`.
- Molecules are independently registered and reactions reference `molecule_id`. CAS numbers were not guessed. Commercial mixtures and incompletely expanded catalysts/reagents retain unresolved representations.
- Candidate isomeric SMILES preserve explicit tetrahedral stereochemistry where transcribed. Products 12/13 retain named 3S allene axial assignments textually with `representation_status=partial`; the toolkit SMILES do not encode those axes.
- Product 1 and 9 absolute descriptors are machine derivations from source-labelled chlorohydrin precursors and epoxide closure, not source certification. All stereo diagnostics are saved in `runs/current/raw_model/pages_0001_0006_standard_stereo_checks.json` with `source_verified=false`.
- Partial bond inventories label only net formed/promoted bonds. Implicit O–H/C–H formation is retained where applicable; reduced C=O bonds and broken bonds are not mislabeled as newly formed sigma bonds. No complete atom mapping or model-vote ensemble was fabricated.

## Preserved conflicts and missing information

- Scheme 2 prints **D-threonine** and **−78 °C** for the KOH stage; the page-4 experimental paragraph prints **L-threonine** and **−30 °C**. Both alternatives and evidence are retained; no single value is chosen.
- Scheme 2 lists HBr, but the detailed page-4 reagent list does not. HBr is retained with unreported dose.
- Scheme 2's **75–83%** range is not allocated to 7a versus 7b and is not substituted for the page-5 83% and 82% step yields.
- Early TiCl4/SnCl4 trials and unsuccessful epoxidation/deprotection groups are described only as “various conditions”; they remain source context rather than invented row-level experiments.
- Missing yield/ee measurement methods remain null when the paper does not report them. Table footnotes support isolated yields/NMR ee only for their stated tables.
- The exact structure/stereo transcription, axial allene assignment encoding, full atom maps, omitted coproducts, and unavailable SI remain expert-review items.

## Artifacts and verification

- Raw visual decisions: `runs/current/raw_model/pages_0001_0006_raw_extraction.json`.
- Reproducible builder: `runs/current/raw_model/build_pages_0001_0006.py`.
- Coverage: `runs/current/coverage/pages_0001_0006.json`.
- Machine submission: `runs/current/submitted.json`; formal merged output is produced at `runs/current/standard_records.json` by `standard merge`.
- Stereo diagnostics: `runs/current/raw_model/pages_0001_0006_standard_stereo_checks.json`.
- Strict validation and audit artifacts are written after merge; their results are appended below without changing the scientific rules.
"""


def main():
    dataset_path = RUN / "dataset.json"
    dataset = json.loads(dataset_path.read_text())
    doc = next(d for d in dataset["documents"] if d["document_id"] == DOC)
    doc["title"] = "Synthesis of (+)-Hypoxylactone through Allenoate γ-Addition: Revision of Stereochemistry"
    dataset_path.write_text(json.dumps(dataset, ensure_ascii=False, indent=2) + "\n")

    bundle = st.empty(dataset["documents"])
    bundle["molecules"] = build_molecules()
    bundle["reactions"] = build_reactions()
    bundle["source_context"] = source_context()
    raw = raw_inventory(bundle["reactions"])

    (RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
    (RUN / "raw_model/pages_0001_0006_raw_extraction.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
    (RUN / "coverage/pages_0001_0006.json").write_text(
        json.dumps(coverage([r["reaction_id"] for r in bundle["reactions"]]), ensure_ascii=False, indent=2) + "\n"
    )
    stereo = {m["molecule_id"]: st.stereo_inspect(st.value(m["smiles"]))
              for m in bundle["molecules"] if m["representation_type"] == "smiles"}
    (RUN / "raw_model/pages_0001_0006_standard_stereo_checks.json").write_text(json.dumps(stereo, ensure_ascii=False, indent=2) + "\n")
    (ROOT / "runs/workflow.md").write_text(workflow_text(bundle, raw))
    print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), **raw["counts"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
