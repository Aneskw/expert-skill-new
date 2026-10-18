#!/usr/bin/env python3
"""Build the page-1--4 source-derived v3.1 submission for doc_923e441b37ab9296."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_923e441b37ab9296"
EV = {i: f"{DOC}:page{i}:image" for i in range(1, 5)}
EVT = {i: f"{DOC}:page{i}:text" for i in range(1, 5)}


def source(raw, normalized, pages, *, derivation=None):
    return st.fact(raw, normalized, [EV[p] for p in pages],
                   method="original_extraction", status="pending_review",
                   verification="unreviewed", derivation=derivation)


def derived(normalized, pages, derivation, raw=None):
    return st.fact(raw, normalized, [EV[p] for p in pages],
                   method="structure_derivation", status="pending_review",
                   verification="inferred", derivation=derivation)


def inferred(normalized, pages, derivation, raw=None):
    return st.fact(raw, normalized, [EV[p] for p in pages],
                   method="inference", status="pending_review",
                   verification="inferred", derivation=derivation)


def missing(pages, status="unreported"):
    return st.fact(None, None, [EV[p] for p in pages], status=status,
                   verification="missing")


def not_applicable(pages):
    return st.fact(None, None, [EV[p] for p in pages], method="not_applicable",
                   status="not_applicable", verification="missing")


def set_occurrence(m, label, pages):
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": source(label, label, pages),
        "evidence_ids": [EV[p] for p in pages],
    }]


def molecule(mid, label, pages, *, smiles=None, name=None, stereo="achiral",
             geometry="not_applicable", representation=None, assignment=None,
             mixture=None, axial=None):
    m = st.new_molecule(mid, DOC, [EV[p] for p in pages])
    set_occurrence(m, label, pages)
    m["identity_status"] = "source_identified"
    m["names"]["original"] = source(label, label, pages)
    m["names"]["english"] = inferred(name, pages,
        "Descriptive English registry name normalized from the source label and displayed structure; it is not an author-supplied systematic name.",
        raw={"source_label": label}) if name else missing(pages)
    m["names"]["abbreviations"] = source(label, [label], pages)
    m["names"]["aliases"] = missing(pages)
    m["cas_numbers"] = missing(pages)
    m["identity_note"] = inferred(
        f"Source label {label}; machine-created identity remains pending review.", pages,
        "Workflow identity note, not source prose.", raw={"source_label": label})
    m["external_properties"] = []
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, pages,
            "Connectivity was transcribed from the displayed structure; RDKit parsing/stereo diagnostics do not certify source identity.", raw={"source_label": label, "structure_location": "displayed molecular drawing"})
        m["generic_template"] = None
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(pages, "unreadable")
        m["generic_template"] = None
    s = m["stereochemistry"]
    s["chirality_status"] = derived(stereo, pages,
        "Assigned only to the extent supported by the displayed graph/wedges or explicit source wording.")
    s["geometry_status"] = derived(geometry, pages,
        "Assigned only to the extent supported by the displayed alkene drawing or explicit source wording.")
    rep = representation or ("not_applicable" if stereo == "achiral" and geometry == "not_applicable" else "not_encoded")
    s["representation_status"] = derived(rep, pages,
        "Toolkit representation status compared with the stereochemical information visible in the source.")
    s["atom_centers"] = missing(pages) if stereo in ("achiral", "not_applicable") else source(
        "Displayed wedges/three contiguous stereogenic centers; individual R/S labels are not printed.",
        [{"count": 3, "configuration": "not individually named"}], pages)
    s["double_bonds"] = missing(pages) if geometry == "not_applicable" else source(
        "Alkene geometry shown by the source drawing; no E/Z descriptor printed.",
        [{"configuration": "drawing-defined; descriptor not printed"}], pages)
    s["axial"] = source(axial, axial, pages) if axial else missing(pages)
    s["mixture"] = source(mixture, mixture, pages) if mixture else missing(pages)
    s["assignment"] = inferred(assignment, pages,
        "Stereochemical assignment normalized from displayed wedges, line geometry, and/or the paper's assignment statement; no additional descriptor is invented.",
        raw={"source_label": label, "structure_location": "displayed molecular drawing or assignment prose"}) if assignment else missing(pages)
    if smiles and rep == "encoded":
        diag = st.stereo_inspect(smiles)
        s["isomeric_smiles"] = derived(diag["canonical_isomeric_smiles"], pages,
            "Canonical isomeric SMILES generated by standard stereo diagnostics from the transcribed candidate SMILES.")
        s["identity_key"] = derived(diag["stereo_graph_key"], pages,
            "Stereo graph key generated by standard stereo diagnostics; it is not a sample identity.")
    else:
        s["isomeric_smiles"] = missing(pages)
        s["identity_key"] = missing(pages)
    return m


def generic_molecule(mid, label, pages, fixed, variables, *, product=False):
    m = molecule(mid, label, pages, stereo="specified" if product else "unspecified",
                 geometry="not_applicable" if product else "unspecified",
                 representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["smiles"] = missing(pages)
    m["generic_template"] = {
        "fixed_structure": source(fixed, {"description": fixed}, pages),
        "variables": [{
            "label": v,
            "attachment_points": source(f"Source drawing attachment for {v}", f"source attachment for {v}", pages),
            "definition": source(f"{v} is varied in the source drawing", f"variable substituent {v}", pages),
            "allowed_values": missing(pages),
            "constraints": missing(pages),
        } for v in variables],
        "instance_mapping": [],
        "reason": source("Generic R-group structure in the source figure.",
                         "Generic R-group structure retained without inventing a concrete molecule.", pages),
    }
    if product:
        m["stereochemistry"]["assignment"] = source(
            "Product drawing shows relative stereochemical wedges; absolute configurations are assigned by analogy to X-ray product 4s.",
            {"relative_stereochemistry": "drawing-defined", "absolute_assignment": "analogy to 4s"}, pages)
    return m


def charge(p, key, raw, value, pages):
    p["charges"][key] = source(raw, value, pages)


def participant(pid, mid, role, pages, *, amount=None, mol_percent=None,
                volume=None, reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, [EV[x] for x in pages])
    if amount:
        charge(p, "amount", amount[0], amount[1], pages)
    if mol_percent:
        charge(p, "mol_percent", mol_percent[0], mol_percent[1], pages)
    if volume:
        charge(p, "volume", volume[0], volume[1], pages)
    if reference:
        p["equivalence_reference"] = derived(reference, pages,
            "Loading denominator is the 0.10 mmol enone specified as the reaction scale in the table footnote.")
    if product_role:
        p["product_role"] = source(product_role, product_role, pages)
    return p


def classification_annulation(pages, limitation=None):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed"},
    ]
    deriv = "The [3+2] annulation drawing closes a five-membered ring through two new C-C sigma connections; atom mapping and accompanying bond-order/H shifts are unresolved."
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(pages),
        "bond_changes": [derived(c, pages, deriv) for c in changes],
        "labels": derived([{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}], pages, deriv),
        "coverage": inferred("partial", pages, deriv),
        "limitations": [limitation or "Partial visual connectivity classification; no complete atom mapping or hydrogen-source assignment."],
        "model_votes": [],
    }


def classification_isomerization(pages):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 3, "after_order": 2, "change_type": "order_decreased"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
         "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
    ]
    deriv = "Visual propargyl-to-allene isomerization: adjacent C-C single bond is promoted to double, C-C triple is reduced to double, with net one-carbon H loss and another-carbon H gain; atom identities and proton shuttle are unresolved."
    labels = [
        {"category_number": 13, "code": "DB-C-C", "supporting_change_ids": ["b1"]},
        {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b3"]},
    ]
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(pages),
            "bond_changes": [derived(c, pages, deriv) for c in changes],
            "labels": derived(labels, pages, deriv), "coverage": inferred("partial", pages, deriv),
            "limitations": ["Net visual isomerization only; atom mapping, proton donor/acceptor, and elementary mechanism are not established."],
            "model_votes": []}


def classification_unresolved(pages, why):
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": missing(pages), "bond_changes": [],
            "labels": derived([], pages, why), "coverage": inferred("unresolved", pages, why),
            "limitations": [why], "model_votes": []}


def classification_none(pages, raw):
    return {"standard_id": st.STANDARD, "mapped_reaction_smiles": not_applicable(pages), "bond_changes": [],
            "labels": derived([], pages, "No observed desired transformation; formation labels are not assigned."),
            "coverage": inferred("not_applicable", pages, "No observed desired transformation; formation labels are not assigned.", raw=raw),
            "limitations": ["Negative control: no bond-formation label inferred from the intended product."], "model_votes": []}


def outcome(oid, kind, target, pages, raw, value, *, basis_raw=None, basis=None,
            yield_type=None, role=None, scope=None):
    o = st.new_outcome(oid, kind, target, [EV[p] for p in pages])
    o["value"] = source(raw, value, pages)
    o["measurement_basis"] = source(basis_raw, basis, pages) if basis_raw is not None else missing(pages)
    o["yield_type"] = source(yield_type, yield_type, pages) if yield_type else not_applicable(pages)
    o["product_role"] = source(role, role, pages) if role else not_applicable(pages)
    o["scope"] = source(scope, scope, pages) if scope else missing(pages)
    return o


def set_common_reaction_fields(r, pages, category, location, *, record_type="specific_experiment"):
    r["record_type"] = record_type
    r["experiment_categories"] = source(category, category, pages)
    r["location"] = source(location, {"pages": pages, "location": location}, pages)
    r["reaction_name"] = source("[3 + 2] annulation", "[3+2] annulation", pages)
    r["transformation_type"] = inferred("phosphine-catalyzed [3+2] annulation", pages,
        "Normalized from the authors' reaction name and displayed transformation.")
    r["completeness"] = "partial"


def add_standard_conditions(r, pages, substrate_mid, enone_mid, product_mid, base_mid, catalyst_mid,
                            solvent_mid, *, base_label, solvent_label, product_label,
                            one_pot=False, omit_base=False, omit_catalyst=False):
    p_sub = participant("substrate", substrate_mid, "reactant", pages,
                        amount=("0.12 mmol", {"number": 0.12, "unit": "mmol"}))
    p_eno = participant("enone", enone_mid, "reactant", pages,
                        amount=("0.10 mmol", {"number": 0.10, "unit": "mmol"}))
    p_prod = participant("product", product_mid, "product", pages, product_role="major diastereomer")
    parts = [p_sub, p_eno, p_prod]
    if not omit_base:
        parts.append(participant("base", base_mid, "catalyst", pages,
                                 amount=("0.13 mmol", {"number": 0.13, "unit": "mmol"})))
    if not omit_catalyst:
        parts.append(participant("phosphine", catalyst_mid, "catalyst", pages,
                                 mol_percent=("10 mol %", {"number": 10, "unit": "mol%"}), reference=enone_mid))
    parts.append(participant("solvent", solvent_mid, "solvent", pages,
                             volume=("1 mL", {"number": 1, "unit": "mL"})))
    r["participants"] = parts
    r["reaction_scale"] = source("2 (0.10 mmol)",
        {"reference_molecule_id": enone_mid, "amount": {"number": 0.10, "unit": "mmol"}}, pages)
    r["solvent_conditions"]["composition"] = source(solvent_label, [solvent_label], pages)
    r["solvent_conditions"]["total_volume"] = source("1 mL", {"number": 1, "unit": "mL"}, pages)
    r["solvent_conditions"]["mixture_ratio"] = not_applicable(pages)
    r["solvent_conditions"]["concentration"] = missing(pages)
    r["solvent_conditions"]["concentration_basis"] = missing(pages)
    if one_pot or omit_base:
        ids = [p["participant_id"] for p in parts]
        s = st.new_step("reaction", [EV[p] for p in pages])
        s["participant_ids"] = ids
        op = "one-step, one-pot addition" if one_pot else "reaction without triethylamine pretreatment"
        s["operation"] = source(op, op, pages)
        s["temperature"] = source("rt", "room temperature", pages)
        s["time"] = source("12 h", {"number": 12, "unit": "h"}, pages)
        r["steps"] = [s]
    else:
        s1 = st.new_step("amine_isomerization", [EV[p] for p in pages])
        s1["participant_ids"] = ["substrate", "base", "solvent"]
        s1["operation"] = source(
            f"treating 1 (0.12 mmol) in {solvent_label} (1 mL) with {base_label} (0.13 mmol)",
            {"operation": "amine/base pretreatment", "addition_order": ["substrate", "solvent", "base"]}, pages)
        s1["temperature"] = source("room temperature", "room temperature", pages)
        s1["time"] = source("12 h", {"number": 12, "unit": "h"}, pages)
        s2 = st.new_step("annulation", [EV[p] for p in pages])
        s2["participant_ids"] = [x for x in ["substrate", "base", "solvent", "enone", "phosphine", "product"] if x in {p["participant_id"] for p in parts}]
        s2["operation"] = source("followed by the addition of 2 and NUSIOC-Phos (10 mol %)",
                                 {"operation": "add enone and phosphine catalyst"}, pages)
        s2["temperature"] = source("rt", "room temperature", pages)
        s2["time"] = source("12 h", {"number": 12, "unit": "h"}, pages)
        r["steps"] = [s1, s2]
    r["procedure"]["original_text"] = source(
        "Reactions were performed by treating 1 (0.12 mmol) in solvent (1 mL) with the base specified (0.13 mmol) and stirred for 12 h at room temperature, followed by the addition of 2 (0.10 mmol) and NUSIOC-Phos (10 mol %).",
        "Two-stage table-footnote procedure; entry-specific exceptions retained separately.", pages)
    r["procedure"]["addition_order"] = source(
        "base pretreatment, followed by enone and phosphine; exceptions noted",
        ["substrate/solvent/base", "12 h room temperature", "enone/phosphine", "12 h room temperature"], pages)
    r["procedure"]["pretreatment"] = source("amine/base pretreatment for 12 h at room temperature",
                                             "amine/base pretreatment for 12 h at room temperature", pages) if not one_pot and not omit_base else not_applicable(pages)
    r["procedure"]["general_procedure_reference"] = missing(pages)
    r["procedure"]["modifications"] = missing(pages)
    r["classification"] = classification_annulation(pages)
    return r


def add_results(r, pages, product_mid, product_label, dr, yld, ee, *, dr_basis, yield_basis, ee_basis,
                extra_note=None):
    outs = []
    if yld == "trace":
        outs.append(outcome("trace_product", "trace", [product_mid], pages, "trace", "trace",
                            role="desired annulation product", scope=product_label))
    elif yld is not None:
        outs.append(outcome("yield", "yield", [product_mid], pages, f"{yld}%", {"number": yld, "unit": "%"},
                            basis_raw=yield_basis, basis="isolated yield of the major diastereomer",
                            yield_type="isolated", role="major diastereomer", scope=product_label))
    if dr is not None:
        outs.append(outcome("dr", "selectivity", [product_mid], pages, dr, dr,
                            basis_raw=dr_basis, basis="crude 1H NMR analysis",
                            scope={"metric": "dr", "order": "reported component order; minor stereoisomers unassigned", "raw_ratio": dr}))
    if ee is not None:
        ev = ee if isinstance(ee, str) else {"number": ee, "unit": "%"}
        outs.append(outcome("ee", "selectivity", [product_mid], pages, f"{ee}% ee", ev,
                            basis_raw=ee_basis, basis="HPLC on a chiral stationary phase",
                            scope={"metric": "ee", "target": product_label, "absolute_configuration": "not inferred from ee"}))
    r["outcomes"] = outs
    if extra_note:
        r["notes"].append(inferred(extra_note, pages,
            "Workflow review note comparing the table cell/footnote with the surrounding source prose; not a verbatim source quotation.",
            raw={"review_basis": "table row, footnote, and surrounding prose"}))
        r["review"]["unresolved"].append(extra_note)


def main():
    dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
    b = st.empty(dataset["documents"])
    mols = []

    # Generic current-method and prior-work structures.
    mols += [
        generic_molecule(f"{DOC}:mol:1_general", "1 (general 3-butynoate)", [1, 2, 3], "R1-C≡C-CH2-CO2Et", ["R1"]),
        generic_molecule(f"{DOC}:mol:2_general", "2 (general β-trifluoromethyl enone)", [1, 2, 3], "CF3-CH=CH-C(=O)-R2", ["R2"]),
        generic_molecule(f"{DOC}:mol:4_general", "4 (general cyclopentene)", [1, 2, 3], "cyclopentene bearing R1, R2, CF3 and CO2Et", ["R1", "R2"], product=True),
        generic_molecule(f"{DOC}:mol:allenoate_general", "γ-substituted allenoate", [1, 2], "R1-CH=C=CH-CO2R2", ["R1", "R2"]),
        generic_molecule(f"{DOC}:mol:loh_enone", "Loh 2010 enone", [1], "R2-CH=CH-C(=O)-R3", ["R2", "R3"]),
        generic_molecule(f"{DOC}:mol:loh_product", "Loh 2010 cyclopentene product", [1], "cyclopentene product bearing R1, R2, R3 and CO2Et", ["R1", "R2", "R3"], product=True),
    ]

    butynoate_smiles = {
        "1a": "CCOC(=O)CC#Cc1ccccc1", "1b": "CCOC(=O)CC#Cc1ccc(Cl)cc1",
        "1c": "CCOC(=O)CC#Cc1cc(Cl)ccc1", "1d": "CCOC(=O)CC#Cc1ccc(F)cc1",
        "1e": "CCOC(=O)CC#Cc1ccccc1F", "1f": "CCOC(=O)CC#Cc1ccc(Br)cc1",
        "1g": "CCOC(=O)CC#Cc1ccc(CC)cc1", "1h": "CCOC(=O)CC#Cc1ccc(C)cc1",
        "1i": "CCOC(=O)CC#Cc1ccc(OC)cc1",
    }
    but_pages = {"1a": [1, 2, 3], **{x: [3] for x in ["1b", "1c", "1d", "1e", "1f", "1g", "1h", "1i"]}}
    for label, smi in butynoate_smiles.items():
        mols.append(molecule(f"{DOC}:mol:{label}", label, but_pages[label], smiles=smi,
                             name=f"3-butynoate {label}", representation="not_applicable"))

    enone_smiles = {
        "2a": "O=C(c1ccccc1)/C=C/C(F)(F)F", "2b": "O=C(c1ccc(Cl)cc1)/C=C/C(F)(F)F",
        "2c": "O=C(c1cc(Cl)ccc1)/C=C/C(F)(F)F", "2d": "O=C(c1ccccc1Cl)/C=C/C(F)(F)F",
        "2e": "O=C(c1ccc(Br)cc1)/C=C/C(F)(F)F", "2f": "O=C(c1ccc(F)cc1)/C=C/C(F)(F)F",
        "2g": "O=C(c1ccc(C#N)cc1)/C=C/C(F)(F)F", "2h": "O=C(c1ccc(C)cc1)/C=C/C(F)(F)F",
        "2i": "O=C(c1ccc(OC)cc1)/C=C/C(F)(F)F", "2j": "O=C(c1cc(OC)c(OC)cc1)/C=C/C(F)(F)F",
        "2k": "O=C(c1ccc2ccccc2c1)/C=C/C(F)(F)F", "2l": "O=C(c1cccs1)/C=C/C(F)(F)F",
    }
    en_pages = {"2a": [1, 2, 3], **{x: [3] for x in ["2b", "2c", "2d", "2e", "2f", "2g", "2h", "2i", "2j", "2k", "2l"]}}
    for label, smi in enone_smiles.items():
        mols.append(molecule(f"{DOC}:mol:{label}", label, en_pages[label], smiles=smi,
                             name=f"β-trifluoromethyl enone {label}", stereo="achiral",
                             geometry="specified", representation="encoded",
                             assignment="E-like connectivity transcribed from the displayed line drawing; no printed E/Z descriptor."))

    assignment = "Three contiguous stereogenic centers are displayed; absolute configurations of cyclopentenes 4 were assigned by analogy based on X-ray analysis of 4s. Individual R/S descriptors are not printed on pages 1–4."
    for idx, label in enumerate([f"4{chr(ord('a') + i)}" for i in range(20)]):
        pages = [1, 2, 3] if label == "4a" else [3]
        mols.append(molecule(f"{DOC}:mol:{label}", label, pages, name=f"cyclopentene {label}",
                             stereo="specified", geometry="not_applicable", representation="not_encoded",
                             assignment=assignment))

    mols.append(molecule(f"{DOC}:mol:1a_prime", "1a′", [2, 3], smiles="CCOC(=O)C=C=Cc1ccccc1",
                         name="γ-phenyl allenoate 1a′", stereo="racemic", geometry="not_applicable",
                         representation="partial", axial="Allene 1a′ is shown as (±) in Figure 2; axial configuration is not individually assigned.",
                         mixture="racemic (±) allene"))

    simple = [
        ("Et3N", "Et3N", [1, 2, 3, 4], "CCN(CC)CC", "triethylamine", "achiral"),
        ("DIPEA", "DIPEA", [3], "CCN(C(C)C)C(C)C", "N,N-diisopropylethylamine", "achiral"),
        ("DABCO", "DABCO", [3], "C1CN2CCN1CC2", "1,4-diazabicyclo[2.2.2]octane", "achiral"),
        ("Cs2CO3", "Cs2CO3", [3], "[Cs+].[Cs+].[O-]C([O-])=O", "cesium carbonate", "achiral"),
        ("toluene", "toluene", [2, 3, 4], "Cc1ccccc1", "toluene", "achiral"),
        ("EtOAc", "EtOAc", [3], "CCOC(C)=O", "ethyl acetate", "achiral"),
        ("CH2Cl2", "CH2Cl2", [3], "ClCCl", "dichloromethane", "achiral"),
        ("CHCl3", "CHCl3", [3], "ClC(Cl)Cl", "chloroform", "achiral"),
        ("dioxane", "dioxane", [3], "O1CCOCC1", "1,4-dioxane", "achiral"),
    ]
    for key, label, pages, smi, name, stereo in simple:
        mols.append(molecule(f"{DOC}:mol:{key}", label, pages, smiles=smi, name=name,
                             stereo=stereo, representation="not_applicable",
                             assignment=("The DABCO cage is treated as achiral by molecular symmetry; RDKit may flag bridgehead candidates without proving sample chirality."
                                         if key == "DABCO" else None)))

    catalysts = ["(R,R)-DIPAMP", "3a", "3b", "3c", "3d", "3e", "3f", "3g", "3h", "NUSIOC-Phos"]
    for cat in catalysts:
        cname = {"3a": "L-valine-derived bifunctional phosphine 3a",
                 "3b": "L-threonine-derived phosphine-amide 3b"}.get(cat, cat)
        mols.append(molecule(f"{DOC}:mol:cat:{cat}", cat, [1, 2] if cat == "(R,R)-DIPAMP" else [2],
                             name=cname, stereo="specified", geometry="not_applicable", representation="not_encoded",
                             assignment="Chiral catalyst structure/configuration is shown in Table 1; candidate SMILES was not safely encoded from the dense drawing."))

    # Failed-scope inputs and expected trace product.
    mols.append(molecule(f"{DOC}:mol:1_nPr", "n-Pr-substituted 3-butynoate", [4],
                         smiles="CCOC(=O)CC#CCCC", name="n-Pr-substituted 3-butynoate",
                         representation="not_applicable"))
    mols.append(molecule(f"{DOC}:mol:2_phenethyl", "phenethyl-substituted β-trifluoromethyl enone", [4],
                         smiles="O=C(CCc1ccccc1)/C=C/C(F)(F)F",
                         name="phenethyl-substituted β-trifluoromethyl enone", geometry="specified",
                         representation="encoded", assignment="Connectivity follows footnote 16 wording and the Table 3 enone template; exact drawn identity is absent on page 4."))
    mols.append(molecule(f"{DOC}:mol:trace_nPr_product", "trace annulation product from n-Pr substrate", [4],
                         name="unassigned trace annulation product", stereo="unresolved", geometry="unresolved",
                         representation="unresolved"))

    b["molecules"] = mols

    reactions = []
    # One current generic scheme, combining repeated graphical/table depictions.
    r = st.new_reaction(f"{DOC}:general_annulation", DOC, [EV[1], EV[2], EV[3]])
    set_common_reaction_fields(r, [1, 2, 3], ["other"], "graphical abstract and Tables 1–3 reaction schemes", record_type="generic_scheme")
    r["participants"] = [
        participant("generic_1", f"{DOC}:mol:1_general", "reactant", [1, 2, 3]),
        participant("generic_2", f"{DOC}:mol:2_general", "reactant", [1, 2, 3]),
        participant("generic_4", f"{DOC}:mol:4_general", "product", [1, 2, 3], product_role="cyclopentene product"),
        participant("generic_base", f"{DOC}:mol:Et3N", "catalyst", [1, 2, 3]),
        participant("generic_phosphine", f"{DOC}:mol:cat:NUSIOC-Phos", "catalyst", [1, 2, 3]),
    ]
    r["classification"] = classification_annulation([1, 2, 3])
    r["review"]["unresolved"] = ["Generic R-group scheme; specific structures and results are represented by separate table-row records."]
    reactions.append(r)

    # Prior literature example in Figure 1.
    r = st.new_reaction(f"{DOC}:literature_loh_2010", DOC, [EV[1]])
    set_common_reaction_fields(r, [1], ["other"], "Figure 1, boxed Loh 2010 prior-work example", record_type="literature_example")
    r["participants"] = [
        participant("prior_1", f"{DOC}:mol:1_general", "reactant", [1]),
        participant("prior_enone", f"{DOC}:mol:loh_enone", "reactant", [1]),
        participant("prior_product", f"{DOC}:mol:loh_product", "product", [1], product_role="prior-work cyclopentene"),
        participant("prior_cat", f"{DOC}:mol:cat:(R,R)-DIPAMP", "catalyst", [1]),
    ]
    r["outcomes"] = [
        outcome("prior_yield", "yield", [f"{DOC}:mol:loh_product"], [1], "66–95%", {"minimum": 66, "maximum": 95, "unit": "%", "raw_text": "66–95%"}, role="prior-work product", scope="Figure 1 literature summary"),
        outcome("prior_ee", "selectivity", [f"{DOC}:mol:loh_product"], [1], "81–99% ee", {"minimum": 81, "maximum": 99, "unit": "%", "raw_text": "81–99% ee"}, scope={"metric": "ee", "source": "Figure 1 literature summary"}),
    ]
    r["classification"] = classification_annulation([1], "Prior-work generic drawing; only the displayed ring-forming C-C connections are classified.")
    r["review"]["unresolved"] = ["This is cited Loh 2010 work, not a new experiment in the current paper; exact substrate identities and measurement bases are not supplied in Figure 1."]
    reactions.append(r)

    # Proposed working hypothesis in Figure 2.
    r = st.new_reaction(f"{DOC}:proposed_binary_catalysis", DOC, [EV[2]])
    set_common_reaction_fields(r, [2], ["mechanistic_study"], "Figure 2 working hypothesis", record_type="proposed_mechanism")
    r["participants"] = [
        participant("hyp_1", f"{DOC}:mol:1_general", "reactant", [2]),
        participant("hyp_allene", f"{DOC}:mol:allenoate_general", "reactant", [2]),
        participant("hyp_2", f"{DOC}:mol:2_general", "reactant", [2]),
        participant("hyp_4", f"{DOC}:mol:4_general", "product", [2], product_role="proposed annulation product"),
        participant("hyp_amine", f"{DOC}:mol:Et3N", "catalyst", [2]),
        participant("hyp_phosphine", f"{DOC}:mol:cat:NUSIOC-Phos", "catalyst", [2]),
    ]
    r["reaction_name"] = source("Our working hypothesis", "amine-promoted isomerization followed by phosphine-catalyzed [3+2] annulation", [2])
    r["classification"] = classification_annulation([2], "Proposed mechanism only; classification reflects the drawn net annulation, not experimentally mapped elementary steps.")
    r["review"]["unresolved"] = ["Author-proposed pathway; Figure 2 is not an independently verified mechanistic experiment."]
    reactions.append(r)

    # Table 1: 14 rows. Entry 11 is shared with Table 2 entry 1 and Table 3 entry 1.
    table1 = [
        (1, "(R,R)-DIPAMP", "3:1", 70, 70), (2, "(R,R)-DIPAMP", ">20:1", 90, 84),
        (3, "3a", "7:1", 82, 12), (4, "3b", "6:1", 74, 12), (5, "3c", ">20:1", 95, 18),
        (6, "3d", ">20:1", 91, 24), (7, "3e", "1:2", 85, 77), (8, "3f", "6:1", 76, 70),
        (9, "3g", "8:1", 81, "−65"), (10, "3h", "9:1", 84, 30),
        (11, "NUSIOC-Phos", "12:1", 87, 99), (12, "NUSIOC-Phos", None, "trace", None),
        (13, None, None, None, None), (14, "NUSIOC-Phos", "11:1:2", 75, 96),
    ]
    for entry, cat, dr, yld, ee in table1:
        rid = f"{DOC}:optimized_4a_shared" if entry == 11 else f"{DOC}:table1_entry{entry:02d}"
        pages = [2, 3] if entry == 11 else [2]
        r = st.new_reaction(rid, DOC, [EV[p] for p in pages])
        cats = ["optimization", "substrate_scope"] if entry == 11 else ["optimization"]
        loc = "Table 1 entry 11; repeated as Table 2 entry 1 and Table 3 entry 1" if entry == 11 else f"Table 1 entry {entry}"
        set_common_reaction_fields(r, pages, cats, loc)
        r["general_scheme_ids"] = [f"{DOC}:general_annulation"]
        if entry == 13:
            # Et3N-only control; observed product is allenoate, not annulation product.
            r["participants"] = [
                participant("substrate", f"{DOC}:mol:1a", "reactant", pages, amount=("0.12 mmol", {"number": 0.12, "unit": "mmol"})),
                participant("enone", f"{DOC}:mol:2a", "reactant", pages, amount=("0.10 mmol", {"number": 0.10, "unit": "mmol"})),
                participant("base", f"{DOC}:mol:Et3N", "catalyst", pages, amount=("0.13 mmol", {"number": 0.13, "unit": "mmol"})),
                participant("solvent", f"{DOC}:mol:toluene", "solvent", pages, volume=("1 mL", {"number": 1, "unit": "mL"})),
                participant("allenoate", f"{DOC}:mol:1a_prime", "product", pages, product_role="isomerization product"),
            ]
            s1 = st.new_step("amine_isomerization", [EV[2]])
            s1["participant_ids"] = ["substrate", "base", "solvent", "allenoate"]
            s1["operation"] = source("reaction was performed in the presence of only Et3N", "Et3N-only isomerization control", pages)
            s1["temperature"] = source("room temperature", "room temperature", pages)
            s1["time"] = source("12 h", {"number": 12, "unit": "h"}, pages)
            r["steps"] = [s1]
            r["reaction_scale"] = source("2a (0.10 mmol)", {"reference_molecule_id": f"{DOC}:mol:2a", "amount": {"number": 0.10, "unit": "mmol"}}, pages)
            r["solvent_conditions"]["composition"] = source("toluene", ["toluene"], pages)
            r["solvent_conditions"]["total_volume"] = source("1 mL", {"number": 1, "unit": "mL"}, pages)
            r["outcomes"] = [outcome("isomerization_only", "other", [f"{DOC}:mol:1a_prime"], pages,
                                      "only isomerization to the corresponding γ-substituted allenoate was observed",
                                      "only isomerization observed", scope="Et3N-only control")]
            r["classification"] = classification_isomerization(pages)
            r["review"]["unresolved"] = ["Table reports dashes rather than a quantified allenoate yield; whether enone 2a was added after pretreatment is not explicit in the footnote wording."]
        else:
            add_standard_conditions(r, pages, f"{DOC}:mol:1a", f"{DOC}:mol:2a", f"{DOC}:mol:4a",
                                    f"{DOC}:mol:Et3N", f"{DOC}:mol:cat:{cat}", f"{DOC}:mol:toluene",
                                    base_label="Et3N", solvent_label="toluene", product_label="4a",
                                    one_pot=(entry == 14), omit_base=(entry in (1, 12)))
            if entry == 1:
                r["procedure"]["modifications"] = source("reaction with 1a, 2a and (R,R)-DIPAMP (10 mol %) in toluene",
                                                         "direct reaction without Et3N", pages)
            if entry == 12:
                r["procedure"]["modifications"] = source("reaction carried out in the presence of NUSIOC-Phos",
                                                         "NUSIOC-Phos present; Et3N omitted", pages)
            if entry == 14:
                r["procedure"]["modifications"] = source("One-pot reaction was carried out.", "Et3N and phosphine added in one step", pages)
            note = None
            if entry == 7:
                note = "The reported dr is 1:2 while the yield footnote says isolated yield of the major diastereomer; pages 1–4 do not identify which stereoisomer is the isolated major product."
            if entry == 14:
                note = "A further unidentified diastereomer was formed; the 11:1:2 component identities are not assigned on pages 1–4."
            add_results(r, pages, f"{DOC}:mol:4a", "4a", dr, yld, ee,
                        dr_basis="determined by crude 1H NMR analysis",
                        yield_basis="isolated yield of the major diastereomer",
                        ee_basis="determined by HPLC analysis on a chiral stationary phase", extra_note=note)
            if entry == 11:
                r["notes"].append(inferred(
                    "Duplicate presentation merged into one experiment record with page-2 and page-3 evidence.", pages,
                    "Table 1 entry 11, Table 2 entry 1, and Table 3 entry 1 have the same substrates, conditions, dr, yield, and ee.",
                    raw={"Table 1 entry 11": ["12:1", "87%", "99% ee"], "Table 2 entry 1": ["12:1", "87%", "99% ee"], "Table 3 entry 1": ["12:1", "87%", "99% ee"]}))
            r["review"]["unresolved"].append("SI experimental workup/purification details were not supplied in this session.")
        reactions.append(r)

    # Table 2 entries 2-8; entry 1 is already merged above.
    table2 = [
        (2, "DIPEA", "toluene", "12:1", 85, 98), (3, "DABCO", "toluene", "5:1:17", 14, None),
        (4, "Cs2CO3", "toluene", "10:1", 71, 99), (5, "Et3N", "EtOAc", "9:1", 80, 98),
        (6, "Et3N", "CH2Cl2", "3:1", 67, 76), (7, "Et3N", "CHCl3", "4:1", 72, 98),
        (8, "Et3N", "dioxane", "11:1", 84, 97),
    ]
    for entry, base, solv, dr, yld, ee in table2:
        r = st.new_reaction(f"{DOC}:table2_entry{entry:02d}", DOC, [EV[3]])
        set_common_reaction_fields(r, [3], ["optimization"], f"Table 2 entry {entry}")
        r["general_scheme_ids"] = [f"{DOC}:general_annulation"]
        add_standard_conditions(r, [3], f"{DOC}:mol:1a", f"{DOC}:mol:2a", f"{DOC}:mol:4a",
                                f"{DOC}:mol:{base}", f"{DOC}:mol:cat:NUSIOC-Phos", f"{DOC}:mol:{solv}",
                                base_label=base, solvent_label=solv, product_label="4a")
        note = None
        if entry == 3:
            note = "Prose states another unidentified diastereomer was the major product, while Table 2 footnote calls the 14% value the isolated yield of the major diastereomer; the target identity of that isolated yield is ambiguous."
        add_results(r, [3], f"{DOC}:mol:4a", "4a", dr, yld, ee,
                    dr_basis="determined by 1H NMR analysis of the crude mixture",
                    yield_basis="isolated yield of the major diastereomer",
                    ee_basis="determined by HPLC analysis on a chiral stationary phase", extra_note=note)
        r["review"]["unresolved"].append("SI experimental workup/purification details were not supplied in this session.")
        reactions.append(r)

    # Table 3 entries 2-20; entry 1 is the shared optimized_4a record.
    scope = [
        (2, "1a", "2b", "4b", "13:1", 89, 99), (3, "1a", "2c", "4c", "13:1", 90, 99),
        (4, "1a", "2d", "4d", "10:1", 75, 91), (5, "1a", "2e", "4e", "16:1", 80, 98),
        (6, "1a", "2f", "4f", "15:1", 92, 99), (7, "1a", "2g", "4g", "15:1", 90, 95),
        (8, "1a", "2h", "4h", "14:1", 91, 99), (9, "1a", "2i", "4i", "8:1", 81, 99),
        (10, "1a", "2j", "4j", "13:1", 72, 99), (11, "1a", "2k", "4k", "17:1", 95, 98),
        (12, "1a", "2l", "4l", "16:1", 93, 98), (13, "1b", "2a", "4m", "9:1", 78, 99),
        (14, "1c", "2a", "4n", "8:1", 73, 98), (15, "1d", "2a", "4o", "10:1", 81, 99),
        (16, "1e", "2a", "4p", "8:1", 70, 99), (17, "1f", "2a", "4q", "8:1", 71, 98),
        (18, "1g", "2a", "4r", "15:1", 89, 99), (19, "1h", "2a", "4s", "16:1", 88, 99),
        (20, "1i", "2a", "4t", "12:1", 82, 99),
    ]
    for entry, sub, eno, prod, dr, yld, ee in scope:
        r = st.new_reaction(f"{DOC}:table3_entry{entry:02d}", DOC, [EV[3]])
        set_common_reaction_fields(r, [3], ["substrate_scope"], f"Table 3 entry {entry}")
        r["general_scheme_ids"] = [f"{DOC}:general_annulation"]
        add_standard_conditions(r, [3], f"{DOC}:mol:{sub}", f"{DOC}:mol:{eno}", f"{DOC}:mol:{prod}",
                                f"{DOC}:mol:Et3N", f"{DOC}:mol:cat:NUSIOC-Phos", f"{DOC}:mol:toluene",
                                base_label="Et3N", solvent_label="toluene", product_label=prod)
        add_results(r, [3], f"{DOC}:mol:{prod}", prod, dr, yld, ee,
                    dr_basis="determined by 1H NMR analysis of the crude mixture",
                    yield_basis="isolated yield of the major diastereomer",
                    ee_basis="determined by HPLC analysis on a chiral stationary phase")
        r["notes"].append(source(
            "Absolute configurations of cyclopentenes 4 were assigned by analogy based on X-ray crystallographic analysis of product 4s.",
            "Absolute stereochemical assignment by analogy to 4s; individual R/S descriptors are not printed.", [3]))
        r["review"]["unresolved"].append("Product connectivity and relative wedges are visible, but exact isomeric SMILES/R,S assignments were not safely encoded from the compact Table 3 drawing.")
        r["review"]["unresolved"].append("SI experimental workup/purification details were not supplied in this session.")
        reactions.append(r)

    # Actual NMR isomerization observation.
    r = st.new_reaction(f"{DOC}:nmr_isomerization_study", DOC, [EV[3]])
    set_common_reaction_fields(r, [3], ["mechanistic_study"], "page 3 right-column 1H NMR isomerization study")
    r["reaction_name"] = source("isomerization process of 3-butynoate 1a facilitated by triethylamine", "Et3N-facilitated propargyl-to-allene isomerization", [3])
    r["transformation_type"] = inferred("propargyl-to-allene isomerization", [3], "Normalized from the observed 1a/1a′ interconversion.")
    r["participants"] = [
        participant("starting_1a", f"{DOC}:mol:1a", "reactant", [3]),
        participant("allenoate_1a_prime", f"{DOC}:mol:1a_prime", "product", [3], product_role="equilibrium allenoate"),
        participant("Et3N", f"{DOC}:mol:Et3N", "catalyst", [3]),
    ]
    s = st.new_step("nmr_monitoring", [EV[3]])
    s["participant_ids"] = ["starting_1a", "allenoate_1a_prime", "Et3N"]
    s["operation"] = source("With the addition of triethylamine, the percentage of allenoate in the system increased rapidly within the first few hours and then reached equilibrium.",
                            "add Et3N and monitor 1a/1a′ by 1H NMR", [3])
    s["temperature"] = missing([3])
    s["time"] = source("around 12 h", {"number": 12, "unit": "h", "qualifier": "approx", "raw_text": "around 12 h"}, [3])
    r["steps"] = [s]
    r["outcomes"] = [
        outcome("initial_allenoate_fraction", "other", [f"{DOC}:mol:1a_prime"], [3], "freshly prepared 1a contains 7.6% of allenoate 1a′", {"number": 7.6, "unit": "%"}, basis_raw="studied by 1H NMR analysis", basis="1H NMR", scope="initial freshly prepared 1a mixture"),
        outcome("equilibrium_allenoate_fraction", "other", [f"{DOC}:mol:1a_prime"], [3], "reached equilibrium (83%) at around 12 h", {"number": 83, "unit": "%"}, basis_raw="studied by 1H NMR analysis", basis="1H NMR", scope="allenoate fraction at equilibrium"),
    ]
    r["classification"] = classification_isomerization([3])
    r["review"]["unresolved"] = ["Amounts, solvent, temperature, and detailed NMR protocol are in unavailable SI and are not reported on pages 1–4."]
    reactions.append(r)

    # Page-4 footnote 16: two distinct unsuccessful scope experiments.
    r = st.new_reaction(f"{DOC}:footnote16_nPr_trace", DOC, [EV[3], EV[4]])
    set_common_reaction_fields(r, [3, 4], ["substrate_scope", "control"], "Table 3 limitation discussed in text and footnote 16: n-Pr substrate")
    r["general_scheme_ids"] = [f"{DOC}:general_annulation"]
    add_standard_conditions(r, [3, 4], f"{DOC}:mol:1_nPr", f"{DOC}:mol:2a", f"{DOC}:mol:trace_nPr_product",
                            f"{DOC}:mol:Et3N", f"{DOC}:mol:cat:NUSIOC-Phos", f"{DOC}:mol:toluene",
                            base_label="Et3N", solvent_label="toluene", product_label="trace product")
    r["outcomes"] = [outcome("trace_product", "trace", [f"{DOC}:mol:trace_nPr_product"], [4],
                              "only a trace amount of product was observed", "trace",
                              scope="n-Pr-substituted 3-butynoate with enone 2a")]
    r["classification"] = classification_unresolved([4], "Trace product connectivity is not drawn or assigned in footnote 16; intended annulation labels are not asserted as observed product structure.")
    r["review"]["unresolved"] = ["Trace product identity, detection method, and SI procedure are not reported on pages 1–4."]
    reactions.append(r)

    r = st.new_reaction(f"{DOC}:footnote16_phenethyl_no_product", DOC, [EV[3], EV[4]])
    set_common_reaction_fields(r, [3, 4], ["substrate_scope", "control"], "Table 3 limitation discussed in text and footnote 16: phenethyl enone")
    r["general_scheme_ids"] = [f"{DOC}:general_annulation"]
    r["participants"] = [
        participant("substrate", f"{DOC}:mol:1a", "reactant", [3, 4], amount=("0.12 mmol", {"number": 0.12, "unit": "mmol"})),
        participant("enone", f"{DOC}:mol:2_phenethyl", "reactant", [3, 4], amount=("0.10 mmol", {"number": 0.10, "unit": "mmol"})),
        participant("base", f"{DOC}:mol:Et3N", "catalyst", [3, 4], amount=("0.13 mmol", {"number": 0.13, "unit": "mmol"})),
        participant("phosphine", f"{DOC}:mol:cat:NUSIOC-Phos", "catalyst", [3, 4], mol_percent=("10 mol %", {"number": 10, "unit": "mol%"}), reference=f"{DOC}:mol:2_phenethyl"),
        participant("solvent", f"{DOC}:mol:toluene", "solvent", [3, 4], volume=("1 mL", {"number": 1, "unit": "mL"})),
    ]
    r["reaction_scale"] = source("enone (0.10 mmol)", {"reference_molecule_id": f"{DOC}:mol:2_phenethyl", "amount": {"number": 0.10, "unit": "mmol"}}, [3, 4])
    r["solvent_conditions"]["composition"] = source("toluene", ["toluene"], [3, 4])
    r["solvent_conditions"]["total_volume"] = source("1 mL", {"number": 1, "unit": "mL"}, [3, 4])
    s1 = st.new_step("amine_isomerization", [EV[3], EV[4]])
    s1["participant_ids"] = ["substrate", "base", "solvent"]
    s1["operation"] = source("Table 3 general condition: treat 1 with Et3N in toluene", "amine pretreatment", [3, 4])
    s1["temperature"] = source("room temperature", "room temperature", [3])
    s1["time"] = source("12 h", {"number": 12, "unit": "h"}, [3])
    s2 = st.new_step("attempted_annulation", [EV[3], EV[4]])
    s2["participant_ids"] = ["substrate", "enone", "base", "phosphine", "solvent"]
    s2["operation"] = source("phenethyl-substituted β-trifluoromethyl enone was reacted with 3-butynoate 1a", "attempted annulation under Table 3 conditions", [4])
    s2["temperature"] = source("rt", "room temperature", [3])
    s2["time"] = source("12 h", {"number": 12, "unit": "h"}, [3])
    r["steps"] = [s1, s2]
    r["outcomes"] = [outcome("no_desired_product", "no_reaction", [], [4], "no desired annulation product was observed", "no desired product observed", scope="phenethyl-substituted β-trifluoromethyl enone with 1a")]
    r["classification"] = classification_none([4], "no desired annulation product was observed")
    r["review"]["unresolved"] = ["Detection method, possible byproducts, conversion, and starting-material recovery are not reported on pages 1–4."]
    reactions.append(r)

    b["reactions"] = reactions
    b["source_context"] = [
        {"context_id": f"{DOC}:context:page1_background", "document_id": DOC,
         "kind": "background_and_prior_work", "evidence_ids": [EV[1]],
         "summary": "Introduction and Figure 1 discuss literature precedents. The boxed Loh 2010 example is retained as a literature_example, not a current-paper experiment."},
        {"context_id": f"{DOC}:context:si_notice", "document_id": DOC,
         "kind": "supplementary_information_notice", "evidence_ids": [EV[3], EV[4]],
         "summary": "The article states that SI exists, but no SI file was explicitly supplied or associated in this session; no SI facts were imported."},
        {"context_id": f"{DOC}:context:page4_bibliography", "document_id": DOC,
         "kind": "bibliography", "evidence_ids": [EV[4]],
         "summary": "Page 4 is chiefly author information, acknowledgments, and references. Bibliographic citations were not converted into current-paper experiments; footnote 16's two current-paper failed reactions were extracted separately."},
        {"context_id": f"{DOC}:context:no_spectra", "document_id": DOC,
         "kind": "characterization_inventory", "evidence_ids": [EV[1], EV[2], EV[3], EV[4]],
         "summary": "No standalone mass-spectral peak table, molecular-weight table, or fragment-relative-abundance table occurs on pages 1–4."},
    ]

    # Raw page inventory and stereo diagnostics are kept separately from the formal bundle.
    raw_dir = ROOT / "runs/current/raw_model"
    raw_dir.mkdir(parents=True, exist_ok=True)
    stereo_checks = {}
    for m in mols:
        smi = st.value(m["smiles"])
        if smi:
            try:
                stereo_checks[m["molecule_id"]] = st.stereo_inspect(smi)
            except ValueError as exc:
                stereo_checks[m["molecule_id"]] = {"error": str(exc)}
    (raw_dir / "stereo_checks.json").write_text(json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")
    raw = {
        "document_id": DOC, "pages": [1, 2, 3, 4],
        "inventory": {
            "page1": ["graphical abstract/current generic scheme", "Figure 1 prior-work activation and Loh 2010 example", "introductory prose"],
            "page2": ["Figure 2 proposed working hypothesis", "Table 1 entries 1-14", "entry-specific footnotes e-h"],
            "page3": ["Table 2 entries 1-8", "Table 3 entries 1-20", "NMR isomerization study", "SI availability notice"],
            "page4": ["author/acknowledgment/reference matter", "footnote 16 two failed reactions", "footnote 18 SI-only preliminary studies"],
        },
        "deduplication": ["Table 1 entry 11 = Table 2 entry 1 = Table 3 entry 1; merged as optimized_4a_shared."],
        "formal_counts": {"molecules": len(mols), "reactions": len(reactions),
                          "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)},
        "limitations": ["No SI was explicitly supplied.", "Compact product/catalyst drawings were not forced into unsupported isomeric SMILES.",
                        "All formal records remain machine-generated and need review."],
    }
    (raw_dir / "pages_0001_0004_extraction.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2) + "\n")
    (ROOT / "runs/current/submitted.json").write_text(json.dumps(b, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(raw["formal_counts"], ensure_ascii=False))


if __name__ == "__main__":
    main()
