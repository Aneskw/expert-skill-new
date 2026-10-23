#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1–6.

Page 9 text is used only as necessary methods context for syntheses introduced on
pages 2–3.  The output remains machine-authored and needs review.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st  # noqa: E402


DOC = "doc_ca1481599a4b0f91"
P = {n: f"{DOC}:page{n}:text" for n in range(1, 13)}
I = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}


def q(number, unit, raw_text=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    if qualifier is not None:
        value["qualifier"] = qualifier
    return value


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def not_applicable(evidence=(), reason=None):
    return st.fact(evidence=evidence, method="not_applicable",
                   status="not_applicable", derivation=reason)


def concrete_molecule(mid, label, name, smiles, evidence, abbreviations=(),
                      identity_note=None, allene_terminal=False):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, evidence),
        "evidence_ids": list(evidence),
    }]
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = (
        f(list(abbreviations), list(abbreviations), evidence)
        if abbreviations else missing(evidence)
    )
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "smiles"
    diagnostic = st.stereo_inspect(smiles)
    canonical = diagnostic["canonical_isomeric_smiles"]
    derivation = (
        "Candidate connectivity was transcribed from the source chemical name or "
        "drawing and canonicalized by RDKit; parsing does not verify source identity."
    )
    m["smiles"] = f(
        {"source_identity": name, "source_label": label}, canonical, evidence,
        method="structure_derivation", derivation=derivation,
    )
    stereo_reason = (
        "Candidate graph has no tetrahedral stereocenter."
        + (" The allene has a terminal CH2 terminus, so its axis is not stereogenic."
           if allene_terminal else "")
    )
    m["stereochemistry"]["atom_centers"] = f(
        {"basis": name}, [], evidence, method="structure_derivation",
        derivation=stereo_reason,
    )
    m["stereochemistry"]["double_bonds"] = f(
        {"basis": name},
        ([{"feature": "terminal allene", "E_Z_applicable": False}]
         if allene_terminal else []),
        evidence, method="structure_derivation",
        derivation=(
            "Carbonyl and terminal-allene double bonds in the candidate graph do "
            "not carry an E/Z assignment."
            if allene_terminal else
            "No source-assigned E/Z alkene is present in the candidate graph."
        ),
    )
    m["stereochemistry"]["axial"] = f(
        {"basis": name}, [], evidence, method="structure_derivation",
        derivation=stereo_reason,
    )
    m["stereochemistry"]["mixture"] = not_applicable(
        evidence, "The named concrete compound is not reported as a stereoisomeric mixture."
    )
    m["stereochemistry"]["assignment"] = not_applicable(
        evidence, "No stereochemical assignment is applicable to the candidate graph."
    )
    m["stereochemistry"]["chirality_status"] = f(
        {"basis": name}, "achiral", evidence, method="structure_derivation",
        derivation=stereo_reason,
    )
    m["stereochemistry"]["geometry_status"] = f(
        {"basis": name}, "not_applicable", evidence,
        method="structure_derivation",
        derivation=(
            "The terminal allene is not E/Z-assignable and the source gives no "
            "other alkene geometry."
            if allene_terminal else
            "No E/Z or ring cis/trans stereochemical element is present in the candidate graph."
        ),
    )
    m["stereochemistry"]["representation_status"] = f(
        {"basis": name}, "not_applicable", evidence,
        method="structure_derivation",
        derivation="The candidate graph has no applicable stereochemical encoding.",
    )
    m["stereochemistry"]["isomeric_smiles"] = f(
        {"candidate_smiles": smiles}, canonical, evidence,
        method="structure_derivation", derivation=derivation,
    )
    m["stereochemistry"]["identity_key"] = f(
        {"canonical_isomeric_smiles": canonical}, diagnostic["stereo_graph_key"],
        evidence, method="structure_derivation",
        derivation="SHA-256 stereo graph key generated by standard stereo diagnostics.",
    )
    m["identity_note"] = (
        f(identity_note, identity_note, evidence)
        if identity_note else
        f(
            {"source_identity": name},
            "Machine-transcribed candidate; source identity has not been independently certified.",
            evidence, method="inference",
            derivation="Named source identity and candidate structure are retained separately.",
        )
    )
    return m


def unresolved_molecule(mid, label, name, evidence, *, chirality="unresolved",
                        geometry="unresolved", representation="unresolved",
                        note, abbreviations=(), mixture=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, evidence),
        "evidence_ids": list(evidence),
    }]
    m["names"]["original"] = f(name, name, evidence)
    m["names"]["english"] = f(name, name, evidence)
    m["names"]["abbreviations"] = (
        f(list(abbreviations), list(abbreviations), evidence)
        if abbreviations else missing(evidence)
    )
    m["names"]["aliases"] = missing(evidence)
    m["cas_numbers"] = missing(evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(evidence)
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = f(
        {"source_identity": name},
        "unresolved; no complete molecular graph and configuration assignment supplied",
        evidence, method="inference", derivation=note,
    )
    m["stereochemistry"]["double_bonds"] = f(
        {"source_identity": name},
        "unresolved; no complete E/Z assignment supplied",
        evidence, method="inference", derivation=note,
    )
    m["stereochemistry"]["axial"] = f(
        {"source_identity": name},
        "unresolved; no complete axial-stereochemistry assignment supplied",
        evidence, method="inference", derivation=note,
    )
    if mixture is None:
        m["stereochemistry"]["mixture"] = missing(evidence)
    else:
        m["stereochemistry"]["mixture"] = f(
            mixture, mixture, evidence, method="inference", derivation=note
        )
    m["stereochemistry"]["assignment"] = f(
        {"source_identity": name}, note, evidence, method="inference", derivation=note
    )
    for key, value in (
        ("chirality_status", chirality),
        ("geometry_status", geometry),
        ("representation_status", representation),
    ):
        m["stereochemistry"][key] = f(
            {"source_identity": name}, value, evidence,
            method="inference", derivation=note,
        )
    m["stereochemistry"]["isomeric_smiles"] = missing(evidence)
    m["stereochemistry"]["identity_key"] = missing(evidence)
    m["identity_note"] = f(note, note, evidence)
    return m


def participant(pid, mid, role, evidence, *, stage="reaction", mass=None,
                amount=None, volume=None, equivalents=None, mol_percent=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, item in (
        ("mass", mass), ("amount", amount), ("volume", volume),
        ("equivalents", equivalents), ("mol_percent", mol_percent),
    ):
        if item is not None:
            raw, normalized = item
            p["charges"][key] = f(raw, normalized, evidence)
    if reference is not None:
        p["equivalence_reference"] = f(
            reference[0], reference[1], evidence,
        )
    return p


def step(step_id, evidence, participant_ids, operation, *, stage="reaction",
         temperature=None, time=None, atmosphere=None, pressure=None):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = f(operation, operation, evidence)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence)
    if atmosphere is not None:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    if pressure is not None:
        s["pressure"] = f(pressure[0], pressure[1], evidence)
    return s


def outcome(oid, kind, targets, evidence, value, *, product_role=None,
            yield_type=None, measurement_basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(value[0], value[1], evidence)
    if product_role is not None:
        o["product_role"] = f(product_role, product_role, evidence)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if measurement_basis is not None:
        o["measurement_basis"] = f(
            measurement_basis, measurement_basis, evidence
        )
    if scope is not None:
        o["scope"] = f(scope, scope, evidence)
    return o


def reaction_base(rid, evidence, location, name, transformation,
                  categories=("other",)):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = f(
        list(categories), list(categories), evidence
    )
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = f(
        transformation, transformation, evidence,
        method="inference",
        derivation="Transformation label summarizes the source-described reaction.",
    )
    r["classification"]["limitations"] = []
    r["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": [],
    }
    r["completeness"] = "partial"
    return r


def set_manual_classification(r, evidence, changes, limitation):
    labels = {}
    for change in changes:
        if change["change_type"] in ("formed", "hydrogen_gain", "order_increased"):
            label = st.bond_label(
                change["atoms"][0]["element"],
                change["atoms"][1]["element"],
                change["after_order"],
            )
            labels.setdefault(
                label["code"], dict(label, supporting_change_ids=[])
            )["supporting_change_ids"].append(change["change_id"])
    r["classification"]["mapped_reaction_smiles"] = missing(evidence)
    r["classification"]["bond_changes"] = [
        f(
            {"source_transformation": r["transformation_type"]["raw_value"]},
            change, evidence, method="structure_derivation",
            derivation=limitation,
        )
        for change in changes
    ]
    r["classification"]["labels"] = f(
        {"taxonomy_basis": limitation}, list(labels.values()), evidence,
        method="structure_derivation", derivation=limitation,
    )
    r["classification"]["coverage"] = f(
        {"basis": limitation}, "partial", evidence, method="inference",
        derivation=limitation,
    )
    r["classification"]["limitations"] = [limitation]
    r["classification"]["model_votes"] = []


def set_mapped_classification(r, evidence, mapped):
    calc = st.classify(mapped)
    derivation = (
        "Atom mapping is a machine derivation from the named 3-butynoic acid, "
        "APTES and BAPTES structures. Dehydration coproducts, coupling-reagent "
        "fate and proton sources are not fully mapped."
    )
    r["classification"]["mapped_reaction_smiles"] = f(
        {"mapping_basis": "source names and BAPTES structure"}, mapped,
        evidence, method="structure_derivation", derivation=derivation,
    )
    r["classification"]["bond_changes"] = [
        f({"mapping_basis": "source names and BAPTES structure"}, change,
          evidence, method="structure_derivation", derivation=derivation)
        for change in calc["bond_changes"]
    ]
    r["classification"]["labels"] = f(
        {"mapping_basis": "source names and BAPTES structure"}, calc["labels"],
        evidence, method="structure_derivation", derivation=derivation,
    )
    r["classification"]["coverage"] = f(
        {"mapping_basis": "partial net mapping"}, "partial", evidence,
        method="inference", derivation=derivation,
    )
    r["classification"]["limitations"] = calc["limitations"] + [derivation]
    r["classification"]["model_votes"] = []


def main():
    dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
    bundle = st.empty(dataset["documents"])

    m = []
    m.append(concrete_molecule(
        f"{DOC}:mol_3_butynoic_acid", "3-butynoic acid", "3-butynoic acid",
        "C#CCC(=O)O", [P[9]],
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_aptes", "APTES", "(3-aminopropyl)triethoxysilane",
        "CCO[Si](OCC)(OCC)CCCN", [P[9]], abbreviations=["APTES"],
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_cmpii", "2-chloro-1-methylpyridinium iodide",
        "2-chloro-1-methylpyridinium iodide", [P[9]],
        note="The reagent is named in the methods, but no molecular drawing was inspected in pages 1–6; a concrete SMILES was not forced.",
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_dipea", "N,N-diisopropylethylamine",
        "N,N-diisopropylethylamine", "CCN(C(C)C)C(C)C", [P[9]],
        abbreviations=["DIPEA"],
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_dcm", "dichloromethane", "dichloromethane", "ClCCl",
        [P[3], P[9]], abbreviations=["CH2Cl2"],
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_ethyl_acetate", "ethyl acetate", "ethyl acetate",
        "CCOC(C)=O", [P[9]],
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_n_hexane", "n-hexane", "n-hexane", "CCCCCC", [P[9]],
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_baptes", "BAPTES",
        "N-(3-(triethoxysilyl)propyl)buta-2,3-dienamide",
        "CCO[Si](OCC)(OCC)CCCNC(=O)C=C=C", [P[2], P[3], P[9]],
        abbreviations=["BAPTES"], allene_terminal=True,
        identity_note="The source names and depicts BAPTES. Candidate SMILES preserves the terminal allene carboxamide; it remains machine-transcribed.",
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_fe3o4_oa_np", "Fe3O4-OA", "oleic-acid-coated Fe3O4 nanoparticles",
        [P[2], P[3], P[9]], abbreviations=["Fe3O4-OA NPs"],
        note="Nanoparticle core, oleate coverage and surface attachment distribution are heterogeneous and cannot be represented as one concrete molecule.",
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_acetic_acid", "acetic acid", "acetic acid", "CC(=O)O",
        [P[9]],
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_np_baptes", "NP-BAPTES", "Fe3O4-BAPTES nanoparticles",
        [P[2], P[3], P[5], P[9]], abbreviations=["NP-BAPTES", "NP-Ctrl"],
        note="The source establishes a BAPTES silane shell but not a unique molecular composition, surface coverage or stereochemical identity.",
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_acetonitrile", "acetonitrile", "acetonitrile", "CC#N",
        [P[9]],
    ))
    m.append(concrete_molecule(
        f"{DOC}:mol_water", "nanopure water", "water", "O", [P[9]],
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_ammonium_carbonate_buffer", "1.0 M ammonium carbonate buffer pH 9.0",
        "ammonium carbonate buffer", [P[9]],
        chirality="not_applicable", geometry="not_applicable",
        note="The source reports a buffered solution rather than a uniquely speciated molecular entity; composition beyond the stated salt, concentration and pH is not expanded.",
        mixture="named aqueous buffer; exact speciation not resolved",
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_ctni_binding_peptide", "HWQIAYNEHQWQC",
        "C-terminal-cysteine cTnI-binding peptide HWQIAYNEHQWQC",
        [P[2], P[3], P[9]], abbreviations=["cTnI-binding peptide"],
        chirality="unspecified", geometry="unresolved", representation="not_encoded",
        note="The amino-acid sequence and C-terminal Cys are explicit, but the supplied source does not state residue L/D configurations or provide a fully stereospecified molecular graph.",
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_baptes_peptide_conjugate", "BAPTES-peptide conjugate",
        "BAPTES-HWQIAYNEHQWQC cysteine-thiol conjugate", [P[2]],
        chirality="unspecified", geometry="unspecified", representation="not_encoded",
        note="The source reports Cys-selective conjugation but the SI spectra/procedure are unavailable. Peptide center configurations and the residual alkene geometry are not assigned.",
    ))
    m.append(unresolved_molecule(
        f"{DOC}:mol_np_pep", "NP-Pep", "Fe3O4-BAPTES-peptide nanoparticles",
        [P[1], P[2], P[3], P[4], P[5], P[6], P[9]],
        abbreviations=["NP-Pep", "Fe3O4-BAPTES-Peptide"],
        chirality="unspecified", geometry="unspecified", representation="not_encoded",
        note="The heterogeneous nanoparticle conjugate has no unique molecular formula. The source drawing does not assign the residual alkene geometry, and peptide residue configurations are not explicit.",
    ))
    bundle["molecules"] = m

    acid = f"{DOC}:mol_3_butynoic_acid"
    aptes = f"{DOC}:mol_aptes"
    cmpii = f"{DOC}:mol_cmpii"
    dipea = f"{DOC}:mol_dipea"
    dcm = f"{DOC}:mol_dcm"
    etoac = f"{DOC}:mol_ethyl_acetate"
    hexane = f"{DOC}:mol_n_hexane"
    baptes = f"{DOC}:mol_baptes"
    feoa = f"{DOC}:mol_fe3o4_oa_np"
    hoac = f"{DOC}:mol_acetic_acid"
    npb = f"{DOC}:mol_np_baptes"
    mecn = f"{DOC}:mol_acetonitrile"
    water = f"{DOC}:mol_water"
    buffer = f"{DOC}:mol_ammonium_carbonate_buffer"
    peptide = f"{DOC}:mol_ctni_binding_peptide"
    soluble_conj = f"{DOC}:mol_baptes_peptide_conjugate"
    nppep = f"{DOC}:mol_np_pep"

    reactions = []

    ev = [P[2], I[2], P[9]]
    r = reaction_base(
        "rxn_baptes_synthesis", ev,
        {"pages": [2, 9], "section": "Results claim on page 2; Methods synthesis on page 9"},
        "Synthesis of N-(3-(triethoxysilyl)propyl)buta-2,3-dienamide",
        "amide coupling followed by propargyl-to-allene isomerization",
    )
    r["participants"] = [
        participant("acid", acid, "reactant", ev,
                    mass=("2.10 g", q(2.10, "g")), amount=("25 mmol", q(25, "mmol")),
                    equivalents=("1 eq.", q(1, "equiv")),
                    reference=("relative to 3-butynoic acid", acid)),
        participant("cmpii", cmpii, "reagent", ev,
                    mass=("6.39 g", q(6.39, "g")), amount=("25 mmol", q(25, "mmol")),
                    equivalents=("1 eq.", q(1, "equiv")),
                    reference=("relative to 3-butynoic acid", acid)),
        participant("dcm_main", dcm, "solvent", ev,
                    volume=("250 mL", q(250, "mL"))),
        participant("aptes", aptes, "reactant", ev,
                    mass=("5.53 g", q(5.53, "g")), amount=("25 mmol", q(25, "mmol")),
                    equivalents=("1 eq.", q(1, "equiv")),
                    reference=("relative to 3-butynoic acid", acid)),
        participant("dipea_initial", dipea, "reagent", ev,
                    mass=("6.46 g", q(6.46, "g")), amount=("50 mmol", q(50, "mmol")),
                    equivalents=("2 eq.", q(2, "equiv")),
                    reference=("relative to 3-butynoic acid", acid)),
        participant("dcm_aptes_solution", dcm, "solvent", ev,
                    volume=("125 mL", q(125, "mL"))),
        participant("dipea_additional", dipea, "reagent", ev,
                    volume=("6.46 mL", q(6.46, "mL")), amount=("50 mmol", q(50, "mmol")),
                    equivalents=("2 eq.", q(2, "equiv")),
                    reference=("relative to 3-butynoic acid", acid)),
        participant("ethyl_acetate", etoac, "purification_agent", ev, stage="purification"),
        participant("n_hexane", hexane, "purification_agent", ev, stage="purification"),
        participant("baptes_product", baptes, "product", ev, product_role="main product"),
    ]
    r["reaction_scale"] = f(
        "3-butynoic acid: 25 mmol", {"reference_molecule_id": acid, "amount": q(25, "mmol")}, ev
    )
    r["solvent_conditions"]["composition"] = f(
        "dichloromethane (250 mL main solution; 125 mL APTES/DIPEA solution)",
        [{"molecule_id": dcm, "role": "main reaction solvent"}], ev,
    )
    r["solvent_conditions"]["total_volume"] = missing(ev)
    r["steps"] = [
        step("baptes_s1_charge_and_reflux", ev, ["acid", "cmpii", "dcm_main"],
             "Charge 3-butynoic acid, 2-chloro-1-methylpyridinium iodide and dichloromethane to an oven-dried three-neck flask; heat to reflux.",
             temperature=("reflux", "reflux"), atmosphere="N2"),
        step("baptes_s2_prepare_and_add_aptes", ev,
             ["aptes", "dipea_initial", "dcm_aptes_solution"],
             "Dilute APTES and the initial DIPEA charge with dichloromethane in a separate flask and add by syringe to the refluxing reaction.",
             temperature=("reflux", "reflux"), atmosphere="N2"),
        step("baptes_s3_first_reflux_hold", ev,
             ["acid", "cmpii", "aptes", "dipea_initial"],
             "Allow the reaction mixture to reflux.",
             temperature=("reflux", "reflux"), time=("1 h", q(1, "h")), atmosphere="N2"),
        step("baptes_s4_additional_base", ev, ["dipea_additional"],
             "Add the additional DIPEA charge to isomerize the propargylic isomer to BAPTES.",
             temperature=("reflux", "reflux"), atmosphere="N2"),
        step("baptes_s5_overnight_reflux", ev,
             ["acid", "aptes", "dipea_additional", "baptes_product"],
             "Allow the reaction to reflux overnight.",
             temperature=("reflux", "reflux"), time=("overnight", "overnight"), atmosphere="N2"),
        step("baptes_s6_concentrate", ev, ["baptes_product"],
             "Concentrate the reaction mixture by rotary evaporation.", stage="workup"),
        step("baptes_s7_redisperse_centrifuge", ev,
             ["baptes_product", "ethyl_acetate"],
             "Redisperse the crude product in ethyl acetate, centrifuge at 5500 × g for 5 min, and concentrate the supernatant.",
             stage="workup", time=("5 min", q(5, "min"))),
        step("baptes_s8_flash_chromatography", ev,
             ["baptes_product", "ethyl_acetate", "n_hexane"],
             "Purify by flash chromatography using an ethyl acetate:n-hexane gradient from 50:50 to 100:0.",
             stage="purification"),
    ]
    r["procedure"]["original_text"] = f(
        "3-Butynoic acid (2.10 g, 1 eq., 25 mmol), 2-chloro-1-methylpyridinium iodide (6.39 g, 1 eq., 25 mmol), and dichloromethane (250 mL) ... APTES (5.53 g, 1 eq., 25 mmol) and N,N-diisopropylethylamine (6.46 g, 2 eq., 50 mmol) ... reflux for 1 h ... N,N-diisopropylethylamine (6.46 mL, 2 eq., 50 mmol) ... reflux overnight ... BAPTES ... (4.80 g, 67% yield).",
        "See ordered structured steps and charges.", ev,
    )
    r["procedure"]["addition_order"] = f(
        "Acid/coupling reagent/DCM; syringe addition of APTES/DIPEA solution; 1 h reflux; additional DIPEA; overnight reflux.",
        ["acid+coupling reagent", "APTES+initial DIPEA", "additional DIPEA"], ev,
    )
    r["procedure"]["pretreatment"] = f(
        "oven dried 500 mL three-neck flask fitted with a reflux condenser",
        "oven-dried flask with reflux condenser", ev,
    )
    r["outcomes"] = [
        outcome("baptes_mass", "product_amount", [baptes], ev,
                ("4.80 g", q(4.80, "g")), product_role="main product",
                measurement_basis="mass of purified clear orange oil"),
        outcome("baptes_yield", "yield", [baptes], ev,
                ("67% yield", q(67, "%")), product_role="main product"),
    ]
    mapped = (
        "[CH:1]#[C:2][CH2:3][C:4](=[O:5])[OH:6]."
        "[NH2:7][CH2:8][CH2:9][CH2:10][Si:11]([O:12][CH2:13][CH3:14])"
        "([O:15][CH2:16][CH3:17])[O:18][CH2:19][CH3:20]>>"
        "[CH2:1]=[C:2]=[CH:3][C:4](=[O:5])[NH:7][CH2:8][CH2:9][CH2:10]"
        "[Si:11]([O:12][CH2:13][CH3:14])([O:15][CH2:16][CH3:17])"
        "[O:18][CH2:19][CH3:20]"
    )
    set_mapped_classification(r, ev, mapped)
    r["notes"] = [
        f(
            "The source gives the initial DIPEA charge as 6.46 g (50 mmol, 2 eq.) and the later DIPEA charge as 6.46 mL (50 mmol, 2 eq.). Both are retained without correction.",
            "Both source-reported DIPEA charges retained; dimensional consistency not adjudicated.", ev,
        ),
        f(
            "Supplementary Figs. 2–6 are cited for BAPTES but no SI file is associated in this workspace.",
            "SI unavailable", [P[2]],
        ),
    ]
    r["review"]["unresolved"] = [
        "Yield type is not explicitly stated as isolated/NMR; the method only reports purified material and a 67% yield.",
        "The two DIPEA additions use different source units; no external density-based correction was made.",
        "The atom mapping omits coupling-reagent fate and coproducts and is therefore partial.",
        "Page 9 was used as methods context but is outside the declared page 1–6 coverage set.",
    ]
    reactions.append(r)

    ev = [P[2], I[2], P[3], I[3], P[9]]
    r = reaction_base(
        "rxn_np_baptes_silanization", ev,
        {"pages": [2, 3, 9], "figure": "Fig. 1a", "method": "Synthesis of Fe3O4-BAPTES NPs"},
        "Synthesis of Fe3O4-BAPTES NPs (NP-BAPTES)",
        "BAPTES silanization of oleic-acid-coated Fe3O4 nanoparticles",
        categories=("scale_up", "other"),
    )
    r["participants"] = [
        participant("fe3o4_oa", feoa, "reactant", ev,
                    mass=("6 mL from a 20 mg/mL stock (120 mg calculated)", q(120, "mg")),
                    volume=("6 mL from a 20 mg/mL stock", q(6, "mL"))),
        participant("hexane_reaction", hexane, "solvent", ev,
                    volume=("anhydrous n-hexane (300 mL)", q(300, "mL"))),
        participant("baptes_silane", baptes, "reagent", ev,
                    volume=("BAPTES (1.65 mL)", q(1.65, "mL"))),
        participant("acetic_acid", hoac, "catalyst", ev,
                    volume=("acetic acid (30 μL)", q(30, "uL"))),
        participant("hexane_wash", hexane, "workup_agent", ev, stage="workup"),
        participant("acetonitrile_wash", mecn, "workup_agent", ev, stage="workup"),
        participant("np_baptes_product", npb, "product", ev, product_role="surface-functionalized nanoparticle product"),
    ]
    r["reaction_scale"] = f(
        "Fe3O4 NPs: 6 mL of 20 mg/mL stock (120 mg by calculation)",
        {"reference_molecule_id": feoa, "amount": q(120, "mg")}, ev,
        method="inference",
        derivation="Calculated as 6 mL × 20 mg/mL; source also states 0.4 mg/mL after dilution.",
    )
    r["solvent_conditions"]["composition"] = f(
        "anhydrous n-hexane", [{"molecule_id": hexane, "anhydrous": True}], ev
    )
    r["solvent_conditions"]["total_volume"] = missing(ev)
    r["solvent_conditions"]["concentration"] = f(
        "total NP concentration of 0.4 mg/mL", "0.4 mg/mL",
        ev,
    )
    r["solvent_conditions"]["concentration_basis"] = f(
        "Fe3O4 NP mass per reaction dispersion volume", feoa, ev
    )
    r["steps"] = [
        step("npb_s1_charge_and_heat", ev, ["fe3o4_oa", "hexane_reaction"],
             "Add the Fe3O4 NP stock to anhydrous n-hexane in a 500 mL round-bottom flask and heat with stirring at 900 rpm.",
             temperature=("60 °C", q(60, "°C"))),
        step("npb_s2_dropwise_additions", ev, ["baptes_silane", "acetic_acid"],
             "Add BAPTES dropwise, followed by dropwise addition of acetic acid.",
             temperature=("60 °C", q(60, "°C"))),
        step("npb_s3_reaction_hold", ev,
             ["fe3o4_oa", "baptes_silane", "acetic_acid", "np_baptes_product"],
             "React under stirring at 900 rpm.", temperature=("60 °C", q(60, "°C")),
             time=("24 h", q(24, "h"))),
        step("npb_s4_collect_and_wash", ev,
             ["np_baptes_product", "hexane_wash", "acetonitrile_wash"],
             "Collect precipitate and wash once with n-hexane, once with n-hexane/acetonitrile (4:1 v/v), and once more with n-hexane by centrifugation at 5500 × g.",
             stage="workup", time=("10 min", q(10, "min"))),
        step("npb_s5_vacuum_dry", ev, ["np_baptes_product"],
             "Dry the nanoparticles under vacuum for later use.", stage="workup"),
    ]
    r["procedure"]["original_text"] = f(
        "Fe3O4 NPs (6 mL from a 20 mg/mL stock) ... anhydrous n-hexane (300 mL) ... heated to 60 °C with stirring (900 rpm), BAPTES (1.65 mL) ... 0.55% (v/v) ... acetic acid (30 μL) ... 0.01% (v/v) ... 24 h ... washed ... centrifugation (5500 × g, 10 min).",
        "See structured charges and ordered steps.", ev,
    )
    r["procedure"]["addition_order"] = f(
        "Fe3O4-OA dispersion; BAPTES dropwise; acetic acid dropwise",
        ["Fe3O4-OA dispersion", "BAPTES", "acetic acid"], ev,
    )
    r["outcomes"] = [
        outcome("npb_silanization_observed", "other", [npb], ev,
                ("successful silanization of the NPs with BAPTES", "NP-BAPTES obtained"),
                product_role="surface-functionalized nanoparticle product",
                measurement_basis="FTIR allene peaks and high-resolution magic-angle-spinning NMR; Fig. 1 characterization"),
    ]
    changes = [{
        "change_id": "b1",
        "atoms": [{"atom_map": None, "element": "Si"}, {"atom_map": None, "element": "O"}],
        "before_order": 0,
        "after_order": 1,
        "change_type": "formed",
    }]
    set_manual_classification(
        r, ev, changes,
        "Silanization is inferred to form Si-O surface/siloxane linkages; the heterogeneous Fe3O4 surface, hydrolysis state and exact attachment atoms are unresolved.",
    )
    r["notes"] = [
        f("BAPTES was reported as 0.55% (v/v) and acetic acid as 0.01% (v/v).",
          {"BAPTES_v_v": "0.55%", "acetic_acid_v_v": "0.01%"}, ev),
        f("The NP-BAPTES is referred to as NP-Ctrl in subsequent figures and text.",
          "NP-BAPTES alias NP-Ctrl", [P[3], I[3]]),
    ]
    r["review"]["unresolved"] = [
        "No isolated mass or chemical yield for NP-BAPTES is reported.",
        "Exact Fe-O-Si versus siloxane-network connectivity and surface coverage are unresolved.",
        "Reaction atmosphere and pressure are not reported.",
        "Page 9 was used as methods context but is outside the declared page 1–6 coverage set.",
    ]
    reactions.append(r)

    ev = [P[2], I[2]]
    r = reaction_base(
        "rxn_soluble_baptes_peptide_test", ev,
        {"page": 2, "section": "Results; BAPTES/Cys-peptide reaction analyzed by HR-MS/MS"},
        "Reaction of BAPTES with C-terminal Cys-modified cTnI-binding peptide",
        "cysteine-thiol addition to an allene carboxamide",
        categories=("mechanistic_study",),
    )
    r["participants"] = [
        participant("soluble_baptes", baptes, "reactant", ev),
        participant("soluble_peptide", peptide, "reactant", ev),
        participant("soluble_conjugate", soluble_conj, "product", ev, product_role="Cys-thiol conjugate"),
    ]
    r["steps"] = [
        step("soluble_click_s1", ev,
             ["soluble_baptes", "soluble_peptide", "soluble_conjugate"],
             "React BAPTES with the C-terminal-Cys-modified peptide; analyze the reaction by high-resolution tandem MS."
        )
    ]
    r["procedure"]["original_text"] = f(
        "we first analyzed the reaction of BAPTES with a C-terminal Cys-modified derivative of the high affinity peptide (HWQIAYNEHQWQC) using high-resolution tandem MS (MS/MS)",
        "BAPTES plus HWQIAYNEHQWQC analyzed by HR-MS/MS", ev,
    )
    r["outcomes"] = [
        outcome("soluble_click_selectivity", "selectivity", [soluble_conj], ev,
                ("occurred exclusively even in the presence of hydroxyls, amines, and carboxylates",
                 "Cys-thiol-selective conjugation; numeric selectivity not reported"),
                product_role="Cys-thiol conjugate",
                measurement_basis="high-resolution tandem mass spectrometry (MS/MS)",
                scope={"metric": "chemoselectivity", "favored_site": "C-terminal Cys thiol", "competing_nucleophiles": ["hydroxyls", "amines", "carboxylates"]}),
    ]
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    set_manual_classification(
        r, ev, changes,
        "Net thiol addition to the allene is structure-derived from the source statement and Fig. 1 chemistry; exact atom mapping and residual alkene geometry are unavailable without SI.",
    )
    r["notes"] = [
        f("Supplementary Figs. 7–8 are cited for this reaction, but no SI file is associated in the workspace.",
          "SI unavailable", ev),
    ]
    r["review"]["unresolved"] = [
        "Reactant amounts, solvent, temperature, duration, workup and product yield are not reported on pages 1–6 and the cited SI is unavailable.",
        "The residual alkene E/Z geometry and peptide residue configurations are not assigned.",
        "Chemoselectivity is qualitative; no numeric ratio or denominator is reported.",
    ]
    reactions.append(r)

    ev = [P[2], I[2], P[3], I[3], P[9]]
    r = reaction_base(
        "rxn_np_pep_coupling", ev,
        {"pages": [2, 3, 9], "figure": "Fig. 1b", "method": "Synthesis of Fe3O4-BAPTES-Peptide NPs"},
        "Synthesis of Fe3O4-BAPTES-Peptide NPs (NP-Pep)",
        "cysteine-thiol addition to surface-bound allene carboxamide",
    )
    r["participants"] = [
        participant("npb_reactant", npb, "reactant", ev,
                    mass=("10 mg of NP-BAPTES", q(10, "mg"))),
        participant("mecn_reaction", mecn, "solvent", ev,
                    volume=("2 mL of acetonitrile", q(2, "mL"))),
        participant("ctni_peptide", peptide, "reactant", ev,
                    mass=("10 mg of cTnI-binding peptide", q(10, "mg"))),
        participant("water_reaction", water, "solvent", ev,
                    volume=("8 mL of nanopure water", q(8, "mL"))),
        participant("carbonate_buffer", buffer, "additive", ev,
                    volume=("75 µL of 1.0 M ammonium carbonate buffer pH 9.0", q(75, "uL"))),
        participant("water_wash", water, "workup_agent", ev, stage="workup"),
        participant("nppep_product", nppep, "product", ev, product_role="peptide-functionalized nanoparticle product"),
    ]
    r["reaction_scale"] = f(
        "10 mg NP-BAPTES and 10 mg peptide", {"reference_molecule_id": npb, "amount": q(10, "mg")}, ev
    )
    r["solvent_conditions"]["composition"] = f(
        "NP-BAPTES in 2 mL acetonitrile; peptide in 8 mL nanopure water; 75 µL aqueous buffer added",
        [{"molecule_id": mecn, "volume": q(2, "mL")}, {"molecule_id": water, "volume": q(8, "mL")}], ev,
    )
    r["solvent_conditions"]["mixture_ratio"] = f(
        "2 mL acetonitrile and 8 mL water before 75 µL buffer addition",
        {"acetonitrile_mL": 2, "water_mL": 8, "buffer_uL": 75}, ev,
    )
    r["solvent_conditions"]["total_volume"] = missing(ev)
    r["steps"] = [
        step("npp_s1_disperse_np", ev, ["npb_reactant", "mecn_reaction"],
             "Disperse 10 mg NP-BAPTES in 2 mL acetonitrile in a 4-dram vial."),
        step("npp_s2_dissolve_peptide", ev, ["ctni_peptide", "water_reaction"],
             "Dissolve 10 mg peptide in 8 mL nanopure water in a separate 4-dram vial."),
        step("npp_s3_adjust_ph", ev, ["ctni_peptide", "carbonate_buffer"],
             "Adjust the peptide solution to pH 8.0 by adding 75 µL 1.0 M ammonium carbonate buffer pH 9.0 during water-bath sonication."),
        step("npp_s4_combine", ev,
             ["npb_reactant", "ctni_peptide", "mecn_reaction", "water_reaction"],
             "Add the pH-adjusted peptide solution to the NP dispersion under water-bath sonication."),
        step("npp_s5_sonication_reaction", ev,
             ["npb_reactant", "ctni_peptide", "nppep_product"],
             "Allow the NP reaction mixture to react under sonication.",
             time=("1 h", q(1, "h"))),
        step("npp_s6_wash_and_isolate", ev, ["nppep_product", "water_wash"],
             "Wash three times with nanopure water by centrifugation at 15,000 × g for 5 min and magnetically isolate to remove unreacted peptide.",
             stage="workup", time=("5 min per centrifugation", q(5, "min"))),
        step("npp_s7_redisperse", ev, ["nppep_product", "water_wash"],
             "Redisperse the peptide-functionalized NPs in nanopure water at 5 mg/mL.",
             stage="workup"),
    ]
    r["procedure"]["original_text"] = f(
        "10 mg of NP-BAPTES ... 2 mL of acetonitrile. 10 mg of cTnI-binding peptide (HWQIAYNEHQWQC) ... 8 mL of nanopure water ... adjusted to pH 8.0 by 75 µL of 1.0 M ammonium carbonate buffer pH 9.0 ... react under sonication for 1 h ... washed three times ... centrifugation (15,000 × g, 5 min) ... redispersed ... 5 mg/mL.",
        "See structured charges and ordered steps.", ev,
    )
    r["procedure"]["addition_order"] = f(
        "Prepare NP dispersion and peptide solution separately; adjust peptide pH; add peptide solution to NP dispersion under sonication.",
        ["NP-BAPTES/acetonitrile", "peptide/water", "buffer", "combine"], ev,
    )
    r["outcomes"] = [
        outcome("npp_coupling_observed", "other", [nppep], ev,
                ("successful consumption of the allene groups upon peptide conjugation", "NP-Pep coupling observed"),
                product_role="peptide-functionalized nanoparticle product",
                measurement_basis="decreased allene FTIR peaks; increased TGA organic content"),
        outcome("npp_surface_density", "other", [nppep], [P[3], I[3]],
                ("~0.034 µmol peptide/mg NP was inferred", "~0.034 µmol peptide/mg NP"),
                product_role="peptide-functionalized nanoparticle product",
                measurement_basis="inferred from the ~6% TGA weight-loss difference between NP-BAPTES and NP-Pep"),
    ]
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "S"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    ]
    set_manual_classification(
        r, ev, changes,
        "Net C-S and C-H formation with one allene C=C order decrease is inferred from the Fig. 1b allene-thiol product drawing. Heterogeneous surface and peptide atoms are not completely mapped.",
    )
    r["notes"] = [
        f("The source reports TGA weight loss of 26% for NP-BAPTES and 32% after peptide coupling.",
          {"NP-BAPTES": q(26, "%"), "NP-Pep": q(32, "%")}, [P[3], I[3]]),
        f("The final NP-Pep dispersion concentration was 5 mg/mL.", "5 mg/mL", [P[9]]),
    ]
    r["review"]["unresolved"] = [
        "No isolated product mass, chemical yield or reaction temperature is reported.",
        "The residual alkene E/Z geometry and peptide residue configurations are not assigned.",
        "Exact surface coverage and atom mapping are unresolved; classification is partial.",
        "Page 9 was used as methods context but is outside the declared page 1–6 coverage set.",
    ]
    reactions.append(r)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {
            "context_id": f"{DOC}:context_page1_abstract",
            "kind": "non_reaction_summary",
            "evidence_ids": [P[1], I[1]],
            "description": "Title/abstract describe the nanoproteomics platform and cTnI enrichment; no reaction arrow, table row or standalone synthetic procedure appears on page 1.",
        },
        {
            "context_id": f"{DOC}:context_fig1_characterization",
            "kind": "characterization",
            "evidence_ids": [P[2], I[2], P[3], I[3]],
            "description": "TEM, FTIR, high-resolution MAS NMR, TGA, zeta-potential and solvent-dispersibility observations characterize Fe3O4-OA, NP-BAPTES and NP-Pep. These observations support synthesis records but are not additional synthetic experiments. FTIR allene bands are reported at about 1970/1971 and 1947 cm-1; TGA weight losses are 26% and 32%.",
        },
        {
            "context_id": f"{DOC}:context_fig2_enrichment",
            "kind": "biochemical_enrichment_and_proteomics",
            "evidence_ids": [P[3], I[3], P[4], I[4], P[5], I[5]],
            "description": "Fig. 2 and surrounding prose report NP binding, magnetic isolation, washing, elution, SDS-PAGE and top-down LC/MS of cTnI proteoforms. These are affinity-enrichment/analytical workflows, not covalent synthesis records. The three NP-Pep batches are reproducibility replicates, not three newly distinct reactions in this page-limited extraction.",
        },
        {
            "context_id": f"{DOC}:context_page5_controls",
            "kind": "binding_controls",
            "evidence_ids": [P[5], I[5]],
            "description": "NP-BAPTES without peptide, acidic-pH functionalization, and negative-control-peptide NPs are described as cTnI-enrichment controls with no appreciable enrichment. The source does not report chemical conversion/yield for these controls on pages 1–6, so they are not recast as no-reaction syntheses.",
        },
        {
            "context_id": f"{DOC}:context_fig3_serum_assay",
            "kind": "bioanalytical_assay",
            "evidence_ids": [P[5], I[5], P[6], I[6]],
            "description": "Fig. 3 reports serum spike-in enrichment, SDS-PAGE, EICs, deconvoluted protein mass spectra and calibration/sensitivity data. Independent spectra, molecular masses and relative abundances are retained as analytical context and not converted into synthetic reactions.",
        },
    ]

    raw = {
        "document_id": DOC,
        "scope_pages": [1, 2, 3, 4, 5, 6],
        "page9_used_as_methods_context": True,
        "semantic_accuracy_verified": False,
        "candidate_bundle": copy.deepcopy(bundle),
        "image_review": {
            "reviewed_images": [I[n] for n in range(1, 7)],
            "detail_crop": "raw_model/crops/page3_fig1ab.png",
            "decision": "Four covalent synthesis/reaction records; protein enrichment and spectra retained in source_context.",
        },
    }
    st.write(ROOT / "runs/current/raw_model/pages_0001_0006_extraction.json", raw)
    st.write(ROOT / "runs/current/submitted.json", bundle)
    st.write(ROOT / "runs/current/raw_model/mentions_baptes_summary.json", {
        "terms": ["BAPTES", "NP-BAPTES", "NP-Pep", "HWQIAYNEHQWQC", "allene-thiol"],
        "relevant_hits": [
            {"evidence_id": P[2], "finding": "BAPTES synthesis claim, soluble BAPTES/Cys-peptide reaction, NP-BAPTES and NP-Pep preparation claims."},
            {"evidence_id": P[3], "finding": "Fig. 1 silanization/allene-thiol drawing and characterization."},
            {"evidence_id": P[5], "finding": "NP-BAPTES/NP-Pep enrichment controls; retained as non-synthetic context."},
            {"evidence_id": P[9], "finding": "Methods for BAPTES, NP-BAPTES and NP-Pep syntheses; context outside coverage pages."},
            {"evidence_id": P[10], "finding": "Enrichment methods; not used as synthesis records."},
        ],
        "limitations": "Literal text search is a locator only; page images 1–6 were inspected separately.",
    })

    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 1,
                "state": "attempted",
                "regions": [
                    {"region": "title and abstract", "finding": "Paper identity and platform summary; no reaction arrow, table, synthetic procedure, or independent characterization experiment.", "evidence_ids": [P[1], I[1]]}
                ],
                "record_ids": [],
                "unresolved": ["Abstract claims are high-level and were not expanded into reaction conditions."],
            },
            {
                "page": 2,
                "state": "attempted",
                "regions": [
                    {"region": "right column, Results—Synthesis and characterization of peptide-functionalized NPs", "finding": "BAPTES synthesis is stated; BAPTES reacts with C-terminal-Cys peptide and is analyzed by HR-MS/MS; NP-BAPTES is functionalized with peptide to form NP-Pep.", "evidence_ids": [P[2], I[2]]},
                    {"region": "right-column characterization continuation", "finding": "TEM/FTIR/MAS-NMR observations support NP-BAPTES and NP-Pep formation; spectral values are characterization, not separate synthesis records.", "evidence_ids": [P[2], I[2]]},
                ],
                "record_ids": [f"{DOC}:rxn_baptes_synthesis", f"{DOC}:rxn_np_baptes_silanization", f"{DOC}:rxn_soluble_baptes_peptide_test", f"{DOC}:rxn_np_pep_coupling"],
                "unresolved": ["Cited Supplementary Figs. 2–10 are not associated in the workspace.", "Detailed methods are on page 9, used only as context outside this page-coverage declaration."],
            },
            {
                "page": 3,
                "state": "attempted",
                "regions": [
                    {"region": "Fig. 1a silanization arrow", "finding": "Fe3O4-OA is silanized with BAPTES to NP-BAPTES; this is duplicate graphical evidence for the specific silanization record, not a second experiment.", "evidence_ids": [P[3], I[3]]},
                    {"region": "Fig. 1b allene-thiol click drawing", "finding": "Surface-bound allene carboxamide reacts with the C-terminal cysteine thiol to form the NP-Pep conjugate; residual alkene geometry is not labeled.", "evidence_ids": [P[3], I[3]]},
                    {"region": "Fig. 1c–h and surrounding text", "finding": "TEM, FTIR, TGA and solvent-dispersibility characterization; TGA supports peptide loading but does not create additional synthesis records.", "evidence_ids": [P[3], I[3]]},
                    {"region": "bottom-right start of cTnI enrichment section", "finding": "Protein enrichment workflow begins; it is biochemical/analytical context rather than a covalent synthesis.", "evidence_ids": [P[3], I[3]]},
                ],
                "record_ids": [f"{DOC}:rxn_np_baptes_silanization", f"{DOC}:rxn_np_pep_coupling"],
                "unresolved": ["Exact surface anchoring connectivity and peptide-conjugate E/Z geometry cannot be established from the heterogeneous drawing."],
            },
            {
                "page": 4,
                "state": "attempted",
                "regions": [
                    {"region": "Fig. 2a", "finding": "Affinity enrichment sequence (loading, NP binding, magnetic isolation, elution, top-down LC/MS); no covalent reaction arrow.", "evidence_ids": [P[4], I[4]]},
                    {"region": "Fig. 2b–h", "finding": "SDS-PAGE, proteoform recovery, enrichment statistics and mass spectra; analytical observations retained in source_context, not synthesis experiments.", "evidence_ids": [P[4], I[4]]},
                ],
                "record_ids": [],
                "unresolved": ["Supplementary Tables 5–6 and source-data files are cited but not supplied as associated SI."],
            },
            {
                "page": 5,
                "state": "attempted",
                "regions": [
                    {"region": "left column, NP-Pep enrichment and controls", "finding": "Protein binding/washing/elution results plus NP-BAPTES, acidic-pH and negative-peptide enrichment controls; chemical conversion and coupling yield are not reported for controls.", "evidence_ids": [P[5], I[5]]},
                    {"region": "left/right columns, LC/MS comparison and heart samples", "finding": "Proteoform distributions and NP-versus-agarose enrichment are analytical results, not new syntheses.", "evidence_ids": [P[5], I[5]]},
                    {"region": "right column, start of serum enrichment", "finding": "Serum spike-in experiment description; no covalent product is formed.", "evidence_ids": [P[5], I[5]]},
                ],
                "record_ids": [],
                "unresolved": ["Control-NP preparation details and chemical coupling conversion are not reported on pages 1–6; no no-reaction synthesis was invented."],
            },
            {
                "page": 6,
                "state": "attempted",
                "regions": [
                    {"region": "Fig. 3a–d", "finding": "SDS-PAGE, extracted-ion chromatograms and deconvoluted protein mass spectra for cTnI/HSA enrichment; characterization only.", "evidence_ids": [P[6], I[6]]},
                    {"region": "Fig. 3e–g", "finding": "Instrument-response and NP-enrichment calibration/sensitivity data; no synthetic transformation.", "evidence_ids": [P[6], I[6]]},
                    {"region": "bottom continuation text", "finding": "Discussion of serum enrichment performance; no additional reaction arrow or table row.", "evidence_ids": [P[6], I[6]]},
                ],
                "record_ids": [],
                "unresolved": ["Reported protein spectra/abundances are not converted into chemical reaction records."],
            },
        ],
    }
    st.write(ROOT / "runs/current/coverage/pages_0001_0006.json", coverage)


if __name__ == "__main__":
    main()
