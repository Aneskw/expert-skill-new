#!/usr/bin/env python3
"""Build the source-derived v3.1 machine submission for doc_041e35a0008dd7ec.

This is a reproducible machine artifact.  It does not certify chemical accuracy.
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

DOC = "doc_041e35a0008dd7ec"
P1 = f"{DOC}:page1:image"
P2 = f"{DOC}:page2:image"
P3 = f"{DOC}:page3:image"
P4 = f"{DOC}:page4:image"
T1 = f"{DOC}:page1:text"
T2 = f"{DOC}:page2:text"
T3 = f"{DOC}:page3:text"
T4 = f"{DOC}:page4:text"


def F(raw=None, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation, status=status)


def D(normalized, ev, derivation):
    return F(None, normalized, ev, method="structure_derivation", derivation=derivation)


def I(normalized, ev, derivation):
    return F(None, normalized, ev, method="inference", derivation=derivation)


def extracted(raw, normalized, ev):
    return F(raw, normalized, ev)


def q(number, unit, raw_text=None, qualifier=None):
    out = {"number": number, "unit": unit}
    if raw_text is not None:
        out["raw_text"] = raw_text
    if qualifier is not None:
        out["qualifier"] = qualifier
    return out


def molecule(mid, label, ev, *, name=None, smiles=None, chirality="unresolved",
             geometry="unresolved", representation=None, atom_centers=None,
             double_bonds=None, axial=None, mixture=None, assignment=None,
             note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = extracted(label, label, ev)
    m["names"]["original"] = extracted(name or label, name or label, ev)
    m["names"]["english"] = extracted(name, name, ev) if name else F(ev=ev)
    m["names"]["abbreviations"] = F(ev=ev)
    m["names"]["aliases"] = F(ev=ev)
    m["cas_numbers"] = F(ev=ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = D(smiles, ev, "Candidate connectivity transcribed from the source drawing/name; RDKit parsing is not source verification.")
        diag = st.stereo_inspect(smiles)
        rep = representation or ("encoded" if not any(p["specified"] == "Unspecified" for p in diag["potential_stereo"]) else "partial")
        m["stereochemistry"]["isomeric_smiles"] = D(diag["canonical_isomeric_smiles"], ev, "RDKit canonical isomeric SMILES generated from the machine-transcribed candidate.")
        m["stereochemistry"]["identity_key"] = D(diag["stereo_graph_key"], ev, "RDKit stereo graph key for the candidate SMILES; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(ev=ev)
        rep = representation or "unresolved"
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = extracted(atom_centers, atom_centers, ev) if atom_centers is not None else F(ev=ev)
    m["stereochemistry"]["double_bonds"] = extracted(double_bonds, double_bonds, ev) if double_bonds is not None else F(ev=ev)
    m["stereochemistry"]["axial"] = extracted(axial, axial, ev) if axial is not None else F(ev=ev)
    m["stereochemistry"]["mixture"] = extracted(mixture, mixture, ev) if mixture is not None else F(ev=ev)
    m["stereochemistry"]["assignment"] = extracted(assignment, assignment, ev) if assignment is not None else F(ev=ev)
    m["stereochemistry"]["chirality_status"] = I(chirality, ev, "Status assigned from the displayed structure and explicit stereochemical wording; absence of a mark is not treated as racemic.")
    m["stereochemistry"]["geometry_status"] = I(geometry, ev, "Status assigned from the displayed structure and explicit E/Z wording.")
    m["stereochemistry"]["representation_status"] = I(rep, ev, "Whether the machine candidate preserves the source-displayed stereochemistry.")
    m["identity_note"] = extracted(note, note, ev) if note else F(ev=ev)
    m["external_properties"] = []
    return m


def participant(pid, mid, role, ev, *, stage="reaction", product_role=None,
                mass=None, amount=None, volume=None, equivalents=None,
                mol_percent=None, reference=None):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    if product_role is not None:
        p["product_role"] = extracted(product_role, product_role, ev)
    for key, val in (("mass", mass), ("amount", amount), ("volume", volume),
                     ("equivalents", equivalents), ("mol_percent", mol_percent)):
        if val is not None:
            raw = val.get("raw_text") or f"{val['number']} {val['unit']}"
            p["charges"][key] = extracted(raw, val, ev)
    if reference is not None:
        p["equivalence_reference"] = extracted(reference, reference, ev)
    return p


def step(sid, ev, pids, operation, *, temperature=None, time=None, stage="reaction"):
    s = st.new_step(sid, ev, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = extracted(operation, operation, ev)
    if temperature is not None:
        raw = temperature.get("raw_text") if isinstance(temperature, dict) else str(temperature)
        s["temperature"] = extracted(raw, temperature, ev)
    if time is not None:
        raw = time.get("raw_text") if isinstance(time, dict) else str(time)
        s["time"] = extracted(raw, time, ev)
    return s


def outcome(oid, kind, targets, ev, raw, normalized, *, product_role=None,
            yield_type=None, measurement=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = extracted(raw, normalized, ev)
    if product_role is not None:
        o["product_role"] = extracted(product_role, product_role, ev)
    if yield_type is not None:
        o["yield_type"] = extracted(yield_type, yield_type, ev)
    if measurement is not None:
        o["measurement_basis"] = extracted(measurement, measurement, ev)
    if scope is not None:
        o["scope"] = extracted(scope, scope, ev)
    return o


def classification(ev, changes=None, *, coverage="partial", limitations=None):
    c = st._skeleton("classification", ev)
    c["standard_id"] = st.STANDARD
    c["mapped_reaction_smiles"] = F(ev=ev)
    c["bond_changes"] = []
    labels = {}
    for idx, spec in enumerate(changes or [], 1):
        change = {
            "change_id": f"b{idx}",
            "atoms": [{"atom_map": None, "element": spec[0]}, {"atom_map": None, "element": spec[1]}],
            "before_order": spec[2],
            "after_order": spec[3],
            "change_type": spec[4],
        }
        if len(spec) > 5:
            change["hydrogen_count"] = spec[5]
        c["bond_changes"].append(D(change, ev, "Net bond change inferred from the source-drawn reactant/product structures; atom correspondence is not a complete mapping."))
        if spec[4] in ("formed", "order_increased", "hydrogen_gain"):
            lab = st.bond_label(spec[0], spec[1], spec[3])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{idx}")
    c["labels"] = D(list(labels.values()), ev, "Taxonomy labels calculated from the explicitly listed inferred net bond changes.")
    c["coverage"] = I(coverage, ev, "Coverage is partial without a source-checked complete atom mapping." if coverage == "partial" else "No observed transformation is asserted for this conceptual record.")
    c["limitations"] = limitations or [
        "No complete atom-mapped reaction SMILES is asserted.",
        "Hydrogen/proton sources and omitted coproducts are unresolved.",
        "Bond changes are machine inferences from the displayed structures and remain needs_review.",
    ]
    c["model_votes"] = []
    return c


def reaction(rid, ev, *, record_type="specific_experiment", categories=None,
             location, name, transformation, participants=(), steps=(), outcomes=(),
             changes=None, classification_coverage="partial", notes=(), unresolved=(),
             solvent=None, concentration=None, concentration_basis=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = extracted(categories or ["other"], categories or ["other"], ev)
    r["location"] = extracted(location, {"location": location}, ev)
    r["general_scheme_ids"] = []
    r["reaction_name"] = extracted(name, name, ev)
    r["transformation_type"] = I(transformation, ev, "Machine normalization of the source-described transformation.")
    r["participants"] = list(participants)
    r["steps"] = list(steps)
    r["outcomes"] = list(outcomes)
    if solvent is not None:
        r["solvent_conditions"]["composition"] = extracted(solvent, solvent, ev)
    if concentration is not None:
        raw = concentration.get("raw_text") or f"{concentration['number']} {concentration['unit']}"
        r["solvent_conditions"]["concentration"] = extracted(raw, concentration, ev)
    if concentration_basis is not None:
        r["solvent_conditions"]["concentration_basis"] = extracted(concentration_basis, concentration_basis, ev)
    r["classification"] = classification(ev, changes, coverage=classification_coverage)
    r["notes"] = [extracted(x, x, ev) for x in notes]
    r["evidence_ids"] = list(ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


def main():
    template = json.loads((ROOT / "runs/current/submitted.json").read_text())
    mols = []
    add = mols.append
    # Source compounds.  Dense products whose exact machine transcription is unsafe remain unresolved.
    add(molecule("m12", "12", [P2], name="methyl aromatic ester 12",
                 smiles="COC(=O)c1cc(C)c(OC)cc1", chirality="achiral", geometry="not_applicable", representation="encoded"))
    add(molecule("m13", "13", [P2], name="benzylic chloride 13",
                 smiles="ClCc1cc(C)c(OC)cc1", chirality="achiral", geometry="not_applicable", representation="encoded"))
    add(molecule("m14", "14", [P2], name="ethyl alkenoate 14",
                 smiles="CCOC(=O)C=Cc1cc(C)c(OC)cc1", chirality="achiral", geometry="unspecified", representation="partial",
                 double_bonds=[{"bond": "alkene in 14", "configuration": "not explicitly labelled E/Z"}],
                 note="Candidate connectivity only; the source drawing has no textual E/Z assignment."))
    add(molecule("m15", "15", [P2], name="carboxylic acid 15",
                 smiles="O=C(O)CCc1cc(C)c(OC)cc1", chirality="achiral", geometry="not_applicable", representation="encoded"))
    add(molecule("m11", "11", [P2], name="tetralone 11",
                 smiles="COc1cc2CCC(=O)c2cc1C", chirality="achiral", geometry="not_applicable", representation="encoded"))
    add(molecule("m16", "16", [P2], name="2-aminotetralone hydrochloride 16",
                 smiles="COc1cc2CCC([NH3+])C(=O)c2cc1C.[Cl-]", chirality="unspecified", geometry="not_applicable", representation="partial",
                 atom_centers=[{"center": "C2 bearing NH3+", "configuration": "not specified"}]))
    add(molecule("m10", "10", [P1, P2, P3], name="formamide 10",
                 smiles="COc1cc2CCC(NC=O)C(=O)c2cc1C", chirality="unspecified", geometry="not_applicable", representation="partial",
                 atom_centers=[{"center": "C2 bearing formamido group", "configuration": "not specified"}]))
    add(molecule("m9", "9", [P1, P3], name="enoate 9",
                 smiles="COC(=O)C=CCC(C)C", chirality="achiral", geometry="unspecified", representation="partial",
                 double_bonds=[{"bond": "enoate C=C", "configuration": "not textually assigned"}]))
    add(molecule("m17", "17", [P2, P3], name="spirocyclic imine 17", chirality="unresolved", geometry="not_applicable",
                 representation="not_encoded", note="Dense source drawing retained without a machine SMILES; exact stereochemical encoding remains unresolved."))
    add(molecule("m18", "18", [P2, P3], name="spiro-gamma-lactam 18", chirality="specified", geometry="not_applicable",
                 representation="not_encoded", atom_centers=[{"source": "Scheme 2 wedge/dash drawing", "configuration": "source-drawn"}],
                 assignment={"statement": "absolute stereochemistry assigned through X-ray analysis", "compound": "18"},
                 note="X-ray-supported source assignment is preserved textually; no exact isomeric SMILES is asserted."))
    add(molecule("m19", "19", [P2, P3], name="spirocyclic imine 19", chirality="unresolved", geometry="not_applicable",
                 representation="not_encoded", note="Dense source drawing retained without a machine SMILES."))
    add(molecule("m20", "20", [P2, P3], name="ring-expanded lactam 20", chirality="specified", geometry="unspecified",
                 representation="not_encoded", atom_centers=[{"source": "Scheme 2 wedge/dash drawing", "configuration": "source-drawn"}],
                 assignment={"statement": "skeletal rearrangement corroborated by X-ray analysis", "compound": "20"},
                 note="X-ray-supported dense structure retained without a machine isomeric SMILES."))
    add(molecule("m21", "21", [P2, P3], name="oxazole 21", chirality="achiral", geometry="not_applicable",
                 representation="not_encoded", assignment={"statement": "identity confirmed by X-ray analysis"},
                 note="Connectivity is shown and X-ray supported, but no machine SMILES is asserted for the fused system."))
    add(molecule("m22", "22", [P2, P3], name="protonated isocyanide intermediate 22", chirality="unresolved", geometry="not_applicable",
                 representation="not_encoded", note="Proposed intermediate, not an isolated source compound."))
    add(molecule("m23", "23", [P3], name="allene 23",
                 smiles="COC(=O)C=C=CCC(C)C", chirality="unspecified", geometry="not_applicable", representation="partial",
                 axial=[{"axis": "C=C=C", "configuration": "not specified in source"}],
                 note="The unsymmetrically substituted allene is potentially axially stereogenic; no R_a/S_a or P/M assignment is reported."))
    add(molecule("m24", "24", [P3], name="extended imine 24 (E/Z mixture)", chirality="partially_specified", geometry="mixture",
                 representation="not_encoded", mixture={"type": "E/Z alkene isomers", "ratio": "unreported", "diastereomers": "one isolated diastereomer"},
                 double_bonds=[{"bond": "exocyclic alkene", "configuration": "E/Z mixture"}],
                 assignment={"statement": "one diastereomer isolated; absolute configuration not stated"},
                 note="Mixture identity is kept distinct from either pure E or pure Z isomer."))
    add(molecule("m25", "25", [P3], name="isobutyl-substituted imine 25", chirality="specified", geometry="not_applicable",
                 representation="not_encoded", atom_centers=[{"source": "Figure 4b wedge/dash drawing", "configuration": "relative source-drawn"}],
                 assignment={"statement": "source-drawn relative stereochemistry; no textual absolute assignment"}))
    add(molecule("m7", "7", [P1, P3], name="northern fragment of ansalactam A (7)", chirality="specified", geometry="not_applicable",
                 representation="not_encoded", atom_centers=[{"source": "Figure 4b wedge/dash drawing", "configuration": "source-drawn"}],
                 assignment={"statement": "single isolated diastereomer; structure displayed with X-ray model", "compound": "7"},
                 note="Dense X-ray-supported stereostructure retained without a machine isomeric SMILES."))

    # Explicit reagents and solvents.
    simple = [
        ("m_LiAlH4", "LiAlH4", "lithium aluminium hydride", "[Li+].[AlH4-]"),
        ("m_PPh3", "PPh3", "triphenylphosphine", "P(c1ccccc1)(c1ccccc1)c1ccccc1"),
        ("m_CCl4", "CCl4", "carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
        ("m_ethyl_acrylate", "ethyl acrylate", "ethyl acrylate", "CCOC(=O)C=C"),
        ("m_PdOAc2", "Pd(OAc)2", "palladium(II) acetate", "CC(=O)[O-].CC(=O)[O-].[Pd+2]"),
        ("m_H2", "H2", "hydrogen", "[H][H]"),
        ("m_NaOH", "NaOH", "sodium hydroxide", "[Na+].[OH-]"),
        ("m_oxalyl_chloride", "(COCl)2", "oxalyl chloride", "O=C(Cl)C(=O)Cl"),
        ("m_TiCl4", "TiCl4", "titanium tetrachloride", "Cl[Ti](Cl)(Cl)Cl"),
        ("m_PTAB", "PTAB", "phenyltrimethylammonium tribromide", "C[N+](C)(C)c1ccccc1.[Br-].BrBr"),
        ("m_NaN3", "NaN3", "sodium azide", "[Na+].[N-]=[N+]=[N-]"),
        ("m_HCl", "HCl", "hydrogen chloride", "Cl"),
        ("m_NaHCO3", "NaHCO3", "sodium bicarbonate", "[Na+].O=C([O-])O"),
        ("m_formic_acetic_anhydride", "formic acetic anhydride", "formic acetic anhydride", "CC(=O)OC=O"),
        ("m_Tf2O", "Tf2O", "triflic anhydride", "O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F"),
        ("m_methyl_acrylate", "methyl acrylate", "methyl acrylate", "COC(=O)C=C"),
        ("m_mCPBA", "mCPBA", "meta-chloroperoxybenzoic acid", "O=C(O)c1cccc(Cl)c1OO"),
        ("m_DCM", "CH2Cl2", "dichloromethane", "ClCCl"),
        ("m_toluene", "toluene", "toluene", "Cc1ccccc1"),
        ("m_EtOH", "EtOH", "ethanol", "CCO"),
        ("m_water", "water", "water", "O"),
    ]
    for mid, label, name, smi in simple:
        if mid == "m_LiAlH4":
            add(molecule(mid, label, [P2, P3], name=name, smiles=smi, chirality="not_applicable", geometry="not_applicable", representation="partial",
                         note="Toolkit flags tetrahedral Al as potential stereo; no molecular stereochemical claim is applicable to this ionic reagent record."))
        else:
            add(molecule(mid, label, [P2, P3], name=name, smiles=smi, chirality="achiral", geometry="not_applicable", representation="encoded"))
    add(molecule("m_PdC", "Pd/C", [P2, P3], name="palladium on carbon", chirality="not_applicable", geometry="not_applicable",
                 representation="unresolved", note="Heterogeneous supported catalyst; composition/loading not reported in the main text."))
    add(molecule("m_wet_Et3N", "wet Et3N", [P2, P3], name="wet triethylamine", chirality="not_applicable", geometry="not_applicable",
                 representation="unresolved", note="Named wet reagent mixture; water content/composition is not reported."))
    add(molecule("m_Pinnick", "Pinnick oxidation conditions", [P3], name="Pinnick oxidation reagent system", chirality="not_applicable", geometry="not_applicable",
                 representation="unresolved", note="The main article names the conditions but does not enumerate the reagent composition."))

    R = []
    # Scheme 1: each arrow/numbered sequence is one experiment; missing SI conditions remain unreported.
    R.append(reaction("scheme1_12_to_13", [P2, T2], location="PDF page 2, Scheme 1, 12 to 13", name="reduction and Appel chlorination",
        transformation="two-step ester reduction/Appel chlorination",
        participants=[participant("p12", "m12", "reactant", [P2]), participant("pLAH", "m_LiAlH4", "reagent", [P2]), participant("pPPh3", "m_PPh3", "reagent", [P2]), participant("pCCl4", "m_CCl4", "reagent", [P2]), participant("p13", "m13", "product", [P2], product_role="main product")],
        steps=[step("reduction", [P2], ["p12", "pLAH"], "1) LiAlH4"), step("appel_chlorination", [P2], ["p12", "pPPh3", "pCCl4", "p13"], "2) PPh3, CCl4")],
        outcomes=[outcome("yield13", "yield", ["m13"], [P2], "13, 99%", q(99, "%", "99%"), product_role="main product", scope={"span": "12 to 13", "steps": 2})],
        changes=[("C", "Cl", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain", 2)],
        unresolved=["Temperature, time, solvent, charges, workup, purification, and yield measurement basis are not reported in the supplied main-text pages."]))
    R.append(reaction("scheme1_13_to_14", [P2, T2], location="PDF page 2, Scheme 1, 13 to 14", name="Heck coupling",
        transformation="palladium-catalyzed coupling with ethyl acrylate",
        participants=[participant("p13", "m13", "reactant", [P2]), participant("pea", "m_ethyl_acrylate", "reactant", [P2]), participant("ppd", "m_PdOAc2", "catalyst", [P2]), participant("p14", "m14", "product", [P2], product_role="main product")],
        steps=[step("coupling", [P2], ["p13", "pea", "ppd", "p14"], "ethyl acrylate, Pd(OAc)2")],
        outcomes=[outcome("yield14", "yield", ["m14"], [P2], "14, 84%", q(84, "%", "84%"), product_role="main product")],
        changes=[("C", "C", 0, 1, "formed")], unresolved=["The source does not textually assign the E/Z geometry of 14.", "Temperature, time, solvent, quantities, and yield measurement basis are unreported."]))
    R.append(reaction("scheme1_14_to_15", [P2, T2], location="PDF page 2, Scheme 1, 14 to 15", name="hydrogenation and saponification",
        transformation="two-step alkene hydrogenation/ester hydrolysis",
        participants=[participant("p14", "m14", "reactant", [P2]), participant("ph2", "m_H2", "reagent", [P2]), participant("ppdc", "m_PdC", "catalyst", [P2]), participant("pnaoh", "m_NaOH", "reagent", [P2]), participant("ph2o", "m_water", "solvent", [P2]), participant("p15", "m15", "product", [P2], product_role="main product")],
        steps=[step("hydrogenation", [P2], ["p14", "ph2", "ppdc"], "1) H2, Pd/C"), step("saponification", [P2], ["p14", "pnaoh", "ph2o", "p15"], "2) NaOH (aq)")],
        outcomes=[outcome("yield15", "yield", ["m15"], [P2], "15, 99%", q(99, "%", "99%"), product_role="main product", scope={"span": "14 to 15", "steps": 2})],
        changes=[("C", "H", 0, 1, "hydrogen_gain", 2), ("O", "H", 0, 1, "hydrogen_gain", 1)], unresolved=["Temperature, time, pressure, quantities, and yield measurement basis are unreported."]))
    R.append(reaction("scheme1_15_to_11", [P2, T2], location="PDF page 2, Scheme 1, 15 to 11", name="Friedel-Crafts acylation",
        transformation="acid activation and intramolecular Friedel-Crafts acylation",
        participants=[participant("p15", "m15", "reactant", [P2]), participant("pox", "m_oxalyl_chloride", "reagent", [P2]), participant("pti", "m_TiCl4", "reagent", [P2]), participant("p11", "m11", "product", [P2], product_role="main product")],
        steps=[step("acid_activation", [P2], ["p15", "pox"], "(COCl)2; then"), step("cyclization", [P2], ["p15", "pti", "p11"], "TiCl4")],
        outcomes=[outcome("yield11", "yield", ["m11"], [P2], "11, 80%", q(80, "%", "80%"), product_role="main product")],
        changes=[("C", "C", 0, 1, "formed")], unresolved=["Temperature, time, solvent, quantities, and yield measurement basis are unreported."]))
    R.append(reaction("scheme1_11_to_16", [P2, T2], location="PDF page 2, Scheme 1, 11 to 16", name="bromination, azide installation, and hydrogenation",
        transformation="three-step alpha amination sequence",
        participants=[participant("p11", "m11", "reactant", [P2]), participant("pptab", "m_PTAB", "reagent", [P2]), participant("pn3", "m_NaN3", "reagent", [P2]), participant("ph2", "m_H2", "reagent", [P2]), participant("ppdc", "m_PdC", "catalyst", [P2]), participant("phcl", "m_HCl", "reagent", [P2]), participant("ph2o", "m_water", "solvent", [P2]), participant("p16", "m16", "product", [P2], product_role="main product")],
        steps=[step("bromination", [P2], ["p11", "pptab"], "1) PTAB"), step("azide_installation", [P2], ["p11", "pn3"], "2) NaN3"), step("azide_hydrogenation", [P2], ["p11", "ph2", "ppdc", "phcl", "ph2o", "p16"], "3) H2, Pd/C, HCl (aq)")],
        outcomes=[outcome("yield16", "yield", ["m16"], [P2], "16, 60%", q(60, "%", "60%"), product_role="main product", scope={"span": "11 to 16", "steps": 3})],
        changes=[("C", "N", 0, 1, "formed"), ("N", "H", 0, 1, "hydrogen_gain", 2)], unresolved=["Temperature, time, solvent, quantities, and yield measurement basis are unreported."]))
    R.append(reaction("scheme1_16_to_10", [P2, T2], location="PDF page 2, Scheme 1, 16 to 10", name="N-formylation",
        transformation="amine N-formylation",
        participants=[participant("p16", "m16", "reactant", [P2]), participant("pnahco3", "m_NaHCO3", "reagent", [P2]), participant("ph2o", "m_water", "solvent", [P2]), participant("pfaa", "m_formic_acetic_anhydride", "reagent", [P2]), participant("p10", "m10", "product", [P2], product_role="main product")],
        steps=[step("formylation", [P2], ["p16", "pnahco3", "ph2o", "pfaa", "p10"], "NaHCO3 (aq), formic acetic anhydride")],
        outcomes=[outcome("yield10", "yield", ["m10"], [P2], "10, 75%", q(75, "%", "75%"), product_role="main product")],
        changes=[("C", "N", 0, 1, "formed")], unresolved=["Temperature, time, solvent, quantities, and yield measurement basis are unreported."]))

    # Scheme 2 arrows.
    R.append(reaction("scheme2_10_to_17", [P2, P3, T2], location="PDF page 3, Scheme 2 upper cycloaddition; discussed on page 2", name="formal [3+2] cycloaddition with methyl acrylate",
        transformation="acid-promoted spirocyclic formal [3+2] cycloaddition",
        participants=[participant("p10", "m10", "reactant", [P3]), participant("pma", "m_methyl_acrylate", "reactant", [P3]), participant("ptf", "m_Tf2O", "reagent", [P3]), participant("pet", "m_wet_Et3N", "additive", [P3]), participant("p17", "m17", "product", [P3], product_role="main product")],
        steps=[step("activation_and_addition", [P3], ["p10", "pma", "ptf"], "Tf2O with methyl acrylate"), step("wet_base_addition", [P3], ["pet", "p17"], "then wet Et3N")],
        outcomes=[outcome("yield17", "yield", ["m17"], [P3], "17, 65%", q(65, "%", "65%"), product_role="main product")],
        changes=[("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed"), ("C", "N", 1, 2, "order_increased")],
        notes=["Page-2 prose says Et3N was not rigorously dried and may function as a proton shuttle/source of trace water."], unresolved=["Temperature, time, solvent, quantities, and yield measurement basis are unreported for this Scheme 2 experiment."]))
    R.append(reaction("scheme2_17_to_18", [P2, P3, T2], location="PDF page 3, Scheme 2 upper oxidation; discussed on page 2", name="mCPBA oxidation of imine 17",
        transformation="imine-to-lactam oxidation",
        participants=[participant("p17", "m17", "reactant", [P3]), participant("pmc", "m_mCPBA", "reagent", [P3]), participant("p18", "m18", "product", [P3], product_role="main product")],
        steps=[step("oxidation", [P3], ["p17", "pmc", "p18"], "mCPBA")],
        outcomes=[outcome("yield18", "yield", ["m18"], [P3], "18, 90%", q(90, "%", "90%"), product_role="main product"), outcome("dr18", "selectivity", ["m18"], [P2, P3], "single isolated diastereomer", "single isolated diastereomer", scope={"metric": "diastereoselectivity", "ratio": "not reported"})],
        changes=[("C", "O", 0, 2, "formed")], notes=["Absolute stereochemistry of 18 was assigned through X-ray analysis."], unresolved=["Reaction conditions, reagent quantity, and yield/selectivity measurement bases are unreported."]))
    R.append(reaction("scheme2_10_9_to_19", [P2, P3, T2], location="PDF page 3, Scheme 2 lower cycloaddition; discussed on page 2", name="formal [3+2] reaction with enoate 9",
        transformation="acid-promoted spirocyclic formal [3+2] cycloaddition",
        participants=[participant("p10", "m10", "reactant", [P3]), participant("p9", "m9", "reactant", [P3]), participant("ptf", "m_Tf2O", "reagent", [P3]), participant("pet", "m_wet_Et3N", "additive", [P3]), participant("p19", "m19", "product", [P3], product_role="main product")],
        steps=[step("activation_and_addition", [P3], ["p10", "p9", "ptf"], "Tf2O with enoate 9"), step("wet_base_addition", [P3], ["pet", "p19"], "then wet Et3N")],
        outcomes=[outcome("yield19", "yield", ["m19"], [P3], "19, 49%", q(49, "%", "49%"), product_role="main product")],
        changes=[("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed"), ("C", "N", 1, 2, "order_increased")], unresolved=["Temperature, time, solvent, quantities, and yield measurement basis are unreported."]))
    R.append(reaction("scheme2_19_to_20", [P2, P3, T2], location="PDF page 3, Scheme 2 lower oxidation; discussed on page 2", name="mCPBA oxidation/rearrangement of imine 19",
        transformation="oxidation accompanied by skeletal rearrangement",
        participants=[participant("p19", "m19", "reactant", [P3]), participant("pmc", "m_mCPBA", "reagent", [P3]), participant("p20", "m20", "product", [P3], product_role="main product")],
        steps=[step("oxidation_and_rearrangement", [P3], ["p19", "pmc", "p20"], "mCPBA")],
        outcomes=[outcome("yield20_step", "yield", ["m20"], [P3], "20, 57%", q(57, "%", "57%"), product_role="main product", scope={"span": "19 to 20", "step": 2}), outcome("yield20_overall", "yield", ["m20"], [T2], "20 in 28% yield over the two steps", q(28, "%", "28%"), product_role="main product", scope={"span_start_molecule_ids": ["m10", "m9"], "span_end_molecule_id": "m20", "steps": 2}), outcome("regio20", "selectivity", ["m20"], [T2], "Only one regioisomer was isolated; overall selectivity uncertain", "one regioisomer isolated; overall selectivity uncertain", scope={"metric": "regioselectivity", "ratio": "unreported"})],
        changes=[("C", "O", 0, 2, "formed"), ("C", "C", 0, 1, "formed")], notes=["The source reports a skeletal rearrangement corroborated by X-ray analysis."], unresolved=["Exact rearrangement atom mapping, conditions, quantities, and measurement bases are unresolved."]))

    # Actual mechanistic observations are distinct from the proposed pathway.
    R.append(reaction("mechanistic_isolation_21", [P2, P3, T2], location="PDF page 3, Figure 3b; discussed on page 2", name="activation of 10 without electrophile",
        transformation="mechanistic trapping/isolation of oxazole 21", categories=["mechanistic_study"],
        participants=[participant("p10", "m10", "reactant", [P3]), participant("ptf", "m_Tf2O", "reagent", [P3], equivalents=q(1.1, "equiv", "1.1 equiv"), reference="m10"), participant("pet", "m_wet_Et3N", "additive", [P3], equivalents=q(0.1, "equiv", "0.1 equiv"), reference="m10"), participant("p21", "m21", "product", [P3], product_role="isolated mechanistic product")],
        steps=[step("activation", [P3], ["p10", "ptf"], "Tf2O (1.1 equiv)"), step("wet_base_and_workup", [P3], ["pet", "p21"], "wet Et3N (0.1 equiv); isolated upon workup")],
        outcomes=[outcome("isolation21", "other", ["m21"], [P2, P3], "yielded oxazole 21; isolated upon workup", "isolated; amount/yield unreported", product_role="isolated mechanistic product", measurement="X-ray identity confirmation")],
        changes=[("C", "O", 0, 1, "formed")], notes=["The source states that activation in the absence of an electrophilic partner yielded oxazole 21, confirmed by X-ray analysis."], unresolved=["Temperature, time, solvent, isolated yield, and complete bond mapping are unreported."]))
    R.append(reaction("mechanistic_21_to_17", [P2, P3, T2], location="PDF page 3, Figure 3b; discussed on page 2", name="competence of isolated oxazole 21 in the formal [3+2] reaction",
        transformation="mechanistic competence test", categories=["mechanistic_study"],
        participants=[participant("p21", "m21", "reactant", [P3]), participant("pma", "m_methyl_acrylate", "reactant", [P3]), participant("p17", "m17", "product", [P3], product_role="observed cycloaddition product")],
        steps=[step("competence_test", [P2, P3], ["p21", "pma", "p17"], "isolated 21 reacted with methyl acrylate under the same conditions")],
        outcomes=[outcome("competent", "other", ["m17"], [T2], "competent reactant in the desired [3+2] reaction", "competent; yield unreported", product_role="observed cycloaddition product")],
        changes=[("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed")], unresolved=["The exact inherited conditions, quantities, time, temperature, yield, and analytical basis are not restated."]))
    mech = reaction("figure3_proposed_pathway", [P2, P3, T2, T3], record_type="proposed_mechanism", categories=["other"], location="PDF page 3, Figure 3a-b and associated discussion", name="proposed protonated-isocyanide/oxazole pathway to 17",
        transformation="proposed mechanism, not an independently verified experiment",
        participants=[participant("p10", "m10", "reactant", [P3]), participant("p22", "m22", "unknown", [P3]), participant("p21", "m21", "unknown", [P3]), participant("pma", "m_methyl_acrylate", "reactant", [P3]), participant("p17", "m17", "product", [P3], product_role="proposed pathway product")],
        steps=[], outcomes=[], changes=[], classification_coverage="not_applicable",
        notes=["Figure 3a is the initially envisioned reactivity; Figure 3b is the authors' proposed pathway.", "The authors explicitly state that the protonated oxazole itself cannot be ruled out as an intermediate."], unresolved=["Elementary steps and the identity of the operative intermediate are proposed rather than established."])
    mech["classification"]["labels"] = D([], [P3], "No bond-formation taxonomy is asserted for a conceptual proposed-mechanism record.")
    mech["classification"]["bond_changes"] = []
    R.append(mech)

    # Figure 4 optimization table: common activation stage plus one independent row per entry.
    rows = [
        (1, 0.1, "CH2Cl2", "m_DCM", 0.1, 16, (15, "%")),
        (2, 0.1, "toluene", "m_toluene", 0.1, 16, ("trace", None)),
        (3, 0.1, "CH2Cl2", "m_DCM", 1, 16, (26, "%")),
        (4, 0.1, "CH2Cl2", "m_DCM", 4, 16, (40, "%")),
        (5, 0.1, "CH2Cl2", "m_DCM", 4, 48, (62, "%")),
        (6, 1.0, "CH2Cl2", "m_DCM", 0.1, 16, ("trace", None)),
    ]
    for entry, eteq, solv, smid, conc, hours, y in rows:
        ev = [P3, T3]
        parts = [participant("p10", "m10", "reactant", ev), participant("ptf", "m_Tf2O", "reagent", ev, equivalents=q(1.1, "equiv", "1.1 equiv"), reference="m10"), participant("p23", "m23", "reactant", ev, equivalents=q(2, "equiv", "2 equiv"), reference="m10"), participant("pet", "m_wet_Et3N", "additive", ev, equivalents=q(eteq, "equiv", f"{eteq:g} equiv"), reference="m10"), participant("psolv", smid, "solvent", ev), participant("p24", "m24", "product", ev, product_role="E/Z product mixture")]
        pre_solvent_pid = "psolv"
        if entry == 2:
            parts.insert(-1, participant("ppre", "m_DCM", "solvent", ev))
            pre_solvent_pid = "ppre"
        steps = [step("preactivation", ev, ["p10", "ptf", pre_solvent_pid], "Tf2O (1.1 equiv), CH2Cl2, 0 °C, 45 min", temperature=q(0, "°C", "0 °C"), time=q(45, "min", "45 min")), step("allene_reaction", ev, ["p10", "p23", "pet", "psolv", "p24"], f"wet Et3N ({eteq:g} equiv), 23 (2 equiv), {solv}, rt, {hours} h", temperature="rt", time=q(hours, "h", f"{hours} h"))]
        outs = []
        if y[0] == "trace":
            outs.append(outcome(f"trace24_e{entry}", "trace", ["m24"], ev, "trace", "trace", product_role="E/Z product mixture", measurement="isolated after column chromatography", scope={"metric": "isolated amount; no numeric yield"}))
        else:
            outs.append(outcome(f"yield24_e{entry}", "yield", ["m24"], ev, f"{y[0]}%", q(y[0], "%", f"{y[0]}%"), product_role="E/Z product mixture", yield_type="isolated yield", measurement="isolated after column chromatography", scope={"mixture": "E/Z isomers"}))
        outs += [outcome(f"ez24_e{entry}", "selectivity", ["m24"], ev, "E/Z isomers isolated", "E/Z mixture; ratio unreported", measurement="column chromatography/isolation", scope={"metric": "E/Z", "isomer_order": "E/Z", "ratio": "unreported"}), outcome(f"dr24_e{entry}", "selectivity", ["m24"], [T3], "only one diastereomer isolated from these reactions", "one diastereomer isolated; ratio unreported", scope={"metric": "diastereoselectivity", "ratio": "unreported"})]
        R.append(reaction(f"figure4_opt_entry_{entry}", ev, categories=["optimization"], location=f"PDF page 3, Figure 4a table, entry {entry}", name="optimization of reaction with allene 23", transformation="acid-promoted spirocyclic formal [3+2] reaction with an allene", participants=parts, steps=steps, outcomes=outs,
            changes=[("C", "C", 0, 1, "formed"), ("C", "C", 0, 1, "formed"), ("C", "N", 1, 2, "order_increased")], solvent={"preactivation": "CH2Cl2", "post-allene conditions": solv}, concentration=q(conc, "M", f"{conc:g} M"),
            notes=["Table footnote a: wet Et3N used in all cases.", "Table footnotes b/c: E/Z isomers isolated after column chromatography; isolated yield of E/Z isomers."], unresolved=["Reaction scale, concentration basis, workup details, and exact E/Z and diastereomer ratios are unreported."] + (["Figure 4a explicitly shows CH2Cl2 for the 0 °C preactivation but lists toluene for entry 2's subsequent rt/16 h conditions; both stages are retained rather than collapsed."] if entry == 2 else [])))

    R.append(reaction("figure4b_24_to_7", [P3, T3], location="PDF page 3, Figure 4b, two-arrow sequence 24 to 25 to 7", name="hydrogenation followed by Pinnick oxidation",
        transformation="two-step exocyclic-alkene hydrogenation/imine oxidation",
        participants=[participant("p24", "m24", "reactant", [P3]), participant("ph2", "m_H2", "reagent", [P3]), participant("ppdc", "m_PdC", "catalyst", [P3]), participant("petoh", "m_EtOH", "solvent", [P3]), participant("p25", "m25", "product", [P3], product_role="intermediate product"), participant("ppin", "m_Pinnick", "reagent", [P3]), participant("p7", "m7", "product", [P3], product_role="final main product")],
        steps=[step("hydrogenation_to_25", [P3], ["p24", "ph2", "ppdc", "petoh", "p25"], "H2, Pd/C, EtOH"), step("pinnick_oxidation_to_7", [P3], ["p25", "ppin", "p7"], "Pinnick oxidation")],
        outcomes=[outcome("yield7_overall", "yield", ["m7"], [P3, T3], "7, 63% over 2 steps", q(63, "%", "63%"), product_role="final main product", scope={"span_start_molecule_id": "m24", "span_end_molecule_id": "m7", "steps": 2}), outcome("dr7", "selectivity", ["m7"], [P3, T3], "single isolated diastereomer", "single isolated diastereomer; ratio unreported", scope={"metric": "diastereoselectivity", "ratio": "unreported"})],
        changes=[("C", "H", 0, 1, "hydrogen_gain", 2), ("C", "O", 0, 2, "formed")], solvent={"step 1": "EtOH", "step 2": "unreported"}, notes=["The 63% value is retained as an overall yield from 24 to 7, not as the yield of the second arrow alone.", "The source states that mCPBA did not provide productive yields for this substrate; no discrete conditions or values are reported."], unresolved=["Times, temperatures, pressure, reagent quantities, Pinnick reagent composition, workup/purification, yield measurement basis, and numerical diastereomer ratio are unreported."]))

    template["molecules"] = mols
    template["reactions"] = R
    template["source_context"] = [
        {"context_id": f"{DOC}:page1_graphical_abstract", "evidence_ids": [P1, T1], "kind": "graphical_summary", "finding": "Graphical abstract/retrosynthetic summary linking ansalactam A, fragment 7, formamide 10 and allene 23; not counted as an additional experiment."},
        {"context_id": f"{DOC}:page1_figure1", "evidence_ids": [P1, T1], "kind": "nonreaction_structures", "finding": "Figure 1 shows ansalactams A-D and is not a synthesis experiment in this article."},
        {"context_id": f"{DOC}:page2_figure2", "evidence_ids": [P2, T2], "kind": "retrosynthesis_and_cited_hypothesis", "finding": "Figure 2a is Trauner's cited biomimetic radical-cyclization logic; Figure 2b is this paper's retrosynthetic proposal. Neither is counted as a current experimental arrow."},
        {"context_id": f"{DOC}:failed_functionalization_18", "evidence_ids": [T2], "kind": "negative_scope_statement", "finding": "All attempts to alkylate or otherwise functionalize 18 to install the isobutyl group failed; individual conditions and experiment identities are not supplied."},
        {"context_id": f"{DOC}:acid_controls", "evidence_ids": [T2], "kind": "control_summary", "finding": "Control experiments reportedly showed no reaction under neutral or basic conditions, but discrete substrates, conditions, and detection basis are not specified, so no fabricated per-condition records were created."},
        {"context_id": f"{DOC}:page4_si_notice", "evidence_ids": [P4, T4], "kind": "metadata", "finding": "Page 4 advertises an SI PDF and CCDC 2120591-2120594; no SI file was supplied or explicitly associated in this session."},
        {"context_id": f"{DOC}:page4_nonexperimental", "evidence_ids": [P4, T4], "kind": "nonexperimental_page", "finding": "Author information, acknowledgments, and references contain no additional current-paper reaction arrows or independent spectral peak tables."},
    ]
    out = ROOT / "runs/current/submitted.json"
    out.write_text(json.dumps(template, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"molecules": len(mols), "reactions": len(R), "source_context": len(template["source_context"])}, ensure_ascii=False))


if __name__ == "__main__":
    main()
