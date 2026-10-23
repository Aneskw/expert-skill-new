#!/usr/bin/env python3
"""Extend the source-derived v3.1 bundle with original PDF pages 7--12.

The script is paper- and page-range-specific.  It starts from the current v3
bundle only to preserve/merge prior-page identities; all facts added here are
derived from pages 7--12 (plus the explicitly referenced Conditions B on page
5) and remain machine records requiring review.
"""
from __future__ import annotations

import copy
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st
from reaction_tools import dataset_audit

DOC = "doc_dfedfa74f4eea49f"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 13)}
T = {p: f"{DOC}:page{p}:text" for p in range(1, 13)}
DIRECT7 = [E[7], T[7]]
COND_B = [E[5], T[5]]


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw:
        value["raw_text"] = raw
    return value


bundle = json.loads((RUN / "standard_records.json").read_text())
bundle = copy.deepcopy(bundle)
# The page-7--12 image evidence was staged after the previous page-range merge.
# Machine merge requires the exact currently prepared document/evidence registry.
bundle["documents"] = copy.deepcopy(
    json.loads((RUN / "dataset.json").read_text())["documents"]
)
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def add_occurrence(mid, label, evidence):
    m = molecules[mid]
    ids = list(dict.fromkeys(evidence))
    if not any(o["document_id"] == DOC and
               o["compound_label"]["normalized_value"] == label and
               o["evidence_ids"] == ids for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": f(label, label, ids),
            "evidence_ids": ids,
        })


def new_mol(mid, label, evidence, *, smiles=None, chirality="unresolved",
            geometry="unresolved", representation="not_encoded",
            assignment=None, double_bonds=None, axial=None, mixture=None,
            identity_note=None):
    if mid in molecules:
        add_occurrence(mid, label, evidence)
        return mid
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(label, label, evidence)
    if label and (label[0].isdigit() or label.startswith("(")):
        m["names"]["abbreviations"] = f(label, [label], evidence)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source_label": label, "candidate_smiles": smiles}, smiles, evidence,
            method="structure_derivation",
            derivation=("Candidate connectivity transcribed from the page-7 drawing; "
                        "RDKit parsing and stereo diagnostics do not establish source identity."),
        )
        m["stereochemistry"]["isomeric_smiles"] = f(
            smiles, diag["canonical_isomeric_smiles"], evidence,
            method="structure_derivation",
            derivation="Canonical isomeric SMILES calculated from the candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = f(
            None, diag["stereo_graph_key"], evidence,
            method="structure_derivation",
            derivation="Stereo graph key calculated from the candidate SMILES; not a sample identity.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f(evidence=evidence)
    m["stereochemistry"]["chirality_status"] = f(
        {"source_structure": label}, chirality, evidence,
        method="structure_derivation",
        derivation="Status is limited to stereochemical marks and labels visible in the cited source drawing.",
    )
    m["stereochemistry"]["geometry_status"] = f(
        {"source_structure": label}, geometry, evidence,
        method="structure_derivation",
        derivation="Status is limited to alkene/ring geometry visible or explicitly reported in the cited source.",
    )
    m["stereochemistry"]["representation_status"] = f(
        {"source_structure": label}, representation, evidence,
        method="structure_derivation",
        derivation="The machine representation does not claim more stereochemistry than is safely encoded.",
    )
    if assignment:
        m["stereochemistry"]["assignment"] = f(assignment, assignment, evidence)
    if double_bonds:
        m["stereochemistry"]["double_bonds"] = f(double_bonds, double_bonds, evidence)
    if axial:
        m["stereochemistry"]["axial"] = f(axial, axial, evidence)
    if mixture:
        m["stereochemistry"]["mixture"] = f(mixture, mixture, evidence)
    note = identity_note or (
        "Source-labelled molecule with an unreviewed candidate graph."
        if smiles else
        "Source-labelled structure; complete connectivity and source stereochemistry were not safely encoded as SMILES."
    )
    m["identity_note"] = f(note, note, evidence)
    molecules[mid] = m
    return mid


def new_rx(suffix, evidence, *, record_type="specific_experiment",
           categories=None, location=None, name=None, transformation=None,
           unresolved=None):
    rid = f"{DOC}:{suffix}"
    if rid in reactions:
        raise ValueError(f"Refuse to recreate existing reaction {rid}")
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = (
        f(categories, categories, evidence) if categories is not None
        else f(evidence=evidence)
    )
    r["location"] = f(location, location, evidence)
    if name:
        r["reaction_name"] = f(name, name, evidence)
    if transformation:
        r["transformation_type"] = f(
            None, transformation, evidence, method="inference",
            derivation="Transformation class derived from the displayed reactant/product relationship.",
        )
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None,
        "checked_fields": [], "unresolved": list(unresolved or []),
    }
    r["completeness"] = "partial"
    reactions[rid] = r
    return r


def addp(r, pid, mid, role, evidence, *, stage="reaction", amount=None,
         equivalents=None, mol_percent=None, volume=None, reference=None,
         product_role=None, method="original_extraction"):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence, method=method)
    for key, item in (("amount", amount), ("equivalents", equivalents),
                      ("mol_percent", mol_percent), ("volume", volume)):
        if item is not None:
            raw, norm = item
            p["charges"][key] = f(raw, norm, evidence, method=method)
    if any(x is not None for x in (equivalents, mol_percent)):
        p["equivalence_reference"] = f(reference, reference, evidence, method=method)
    r["participants"].append(p)
    return p


def add_step(r, sid, evidence, *, operation=None, temperature=None, time=None,
             participant_ids=None, stage="reaction", method="original_extraction"):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids or [])
    if operation is not None:
        s["operation"] = f(operation, operation, evidence, method=method)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], evidence, method=method)
    if time is not None:
        s["time"] = f(time[0], time[1], evidence, method=method)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, targets, evidence, raw, normalized, *, role=None,
            yield_type=None, basis=None, scope=None, method="original_extraction"):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence, method=method)
    if role:
        o["product_role"] = f(role, role, evidence, method=method)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence, method=method)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, evidence, method=method)
    if scope:
        o["scope"] = f(scope, scope, evidence, method=method)
    r["outcomes"].append(o)
    return o


def note(r, text, evidence, *, method="original_extraction"):
    r["notes"].append(f(text, text, evidence, method=method))


def classify(r, evidence, specs=None, coverage="partial", limitation=None):
    c = r["classification"]
    if coverage in {"unresolved", "not_applicable"}:
        c["coverage"] = f(
            None, coverage, evidence, method="inference",
            derivation="No source-observed, fully resolved bond-change experiment is asserted for this record.",
        )
        c["labels"] = f(
            None, [] if coverage == "not_applicable" else None, evidence,
            method="inference",
            derivation="No applicable observed-product formation labels are asserted.",
        )
        c["limitations"] = [limitation or "Proposed/unresolved transformation without complete atom mapping."]
        return
    changes = []
    grouped = {}
    for i, (a, b, before, after, kind) in enumerate(specs or [], 1):
        change = {
            "change_id": f"b{i}",
            "atoms": [{"atom_map": None, "element": a},
                      {"atom_map": None, "element": b}],
            "before_order": before, "after_order": after,
            "change_type": kind,
        }
        changes.append(change)
        if kind in {"formed", "order_increased", "hydrogen_gain"}:
            label = st.bond_label(a, b, after)
            grouped.setdefault(label["code"], {
                **label, "supporting_change_ids": []
            })["supporting_change_ids"].append(change["change_id"])
    derivation = "Partial bond changes derived from the page-7 displayed sequence; atom identities are not mapped."
    c["bond_changes"] = [
        f(None, change, evidence, method="structure_derivation", derivation=derivation)
        for change in changes
    ]
    c["labels"] = f(
        None, list(grouped.values()), evidence,
        method="structure_derivation", derivation=derivation,
    )
    c["coverage"] = f(
        None, coverage, evidence, method="inference",
        derivation="Only source-visible transformation classes are represented; full atom mapping is unavailable.",
    )
    c["limitations"] = [
        limitation or "No complete atom mapping; coproducts, leaving groups, and exact implicit-H balance remain unresolved."
    ]


CARBO = [
    ("C", "C", 0, 1, "formed"),
    ("C", "B", 0, 1, "formed"),
    ("C", "C", 2, 1, "order_decreased"),
]
OXIDATIVE_ISOLATION = CARBO + [
    ("C", "B", 1, 0, "broken"),
    ("C", "O", 0, 1, "formed"),
    ("O", "H", 0, 1, "formed"),
]
ALLYLBORATION = [
    ("C", "C", 0, 1, "formed"),
    ("C", "C", 1, 2, "order_increased"),
    ("C", "C", 2, 1, "order_decreased"),
    ("C", "O", 2, 1, "order_decreased"),
    ("O", "H", 0, 1, "formed"),
    ("C", "B", 1, 0, "broken"),
]


def condition_b(r, allene_mid, halide_mid):
    """Apply the explicitly cross-referenced Conditions B from Table 2."""
    ev = COND_B
    addp(r, "allene", allene_mid, "reactant", ev,
         amount=("0.5 mmol", q(0.5, "mmol")),
         equivalents=("1.0 equiv", q(1.0, "equiv")),
         reference=allene_mid, method="inherited_general_condition")
    addp(r, "halide", halide_mid, "reactant", DIRECT7,
         amount=("1.0 mmol under Conditions B", q(1.0, "mmol")),
         equivalents=("2.0 equiv", q(2.0, "equiv")),
         reference=allene_mid, method="inherited_general_condition")
    addp(r, "diboron", f"{DOC}:mol_B2pin2", "reagent", ev,
         amount=("0.6 mmol", q(0.6, "mmol")),
         equivalents=("1.2 equiv", q(1.2, "equiv")),
         reference=allene_mid, method="inherited_general_condition")
    addp(r, "catalyst", f"{DOC}:mol_CuCl", "catalyst", ev,
         amount=("0.025 mmol", q(0.025, "mmol")),
         mol_percent=("5 mol%", q(5, "mol%")),
         reference=allene_mid, method="inherited_general_condition")
    addp(r, "ligand", f"{DOC}:mol_Xantphos", "ligand", ev,
         amount=("0.025 mmol", q(0.025, "mmol")),
         mol_percent=("5 mol%", q(5, "mol%")),
         reference=allene_mid, method="inherited_general_condition")
    addp(r, "base", f"{DOC}:mol_KOtBu", "reagent", ev,
         amount=("0.6 mmol", q(0.6, "mmol")),
         equivalents=("1.2 equiv", q(1.2, "equiv")),
         reference=allene_mid, method="inherited_general_condition")
    addp(r, "solvent", f"{DOC}:mol_DMF", "solvent", ev,
         volume=("1.0 mL", q(1.0, "mL")),
         method="inherited_general_condition")
    r["reaction_scale"] = f(
        "allene (0.5 mmol) under Conditions B",
        {"reference_molecule_id": allene_mid, "amount": q(0.5, "mmol")},
        ev, method="inherited_general_condition",
    )
    r["solvent_conditions"]["composition"] = f(
        "DMF", ["DMF"], ev, method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = f(
        "1.0 mL", q(1.0, "mL"), ev, method="inherited_general_condition")
    r["solvent_conditions"]["concentration"] = f(
        "0.5 M (0.5 mmol allene in 1.0 mL DMF)", q(0.5, "M"), ev,
        method="inherited_general_condition")
    r["solvent_conditions"]["concentration_basis"] = f(
        "allene", allene_mid, ev, method="inherited_general_condition")
    add_step(
        r, "conditions_B_reaction", ev,
        operation="carboboration under Conditions B",
        temperature=("30 °C", q(30, "°C")),
        time=("24 h", q(24, "h")),
        participant_ids=[p["participant_id"] for p in r["participants"]],
        method="inherited_general_condition",
    )


# New page-7 molecule identities.
m1l = new_mol(
    f"{DOC}:mol_1l", "1l", [E[7]],
    smiles="C=C=C1CCC(c2ccccc2)CC1",
    chirality="unspecified", geometry="unspecified", representation="partial",
    axial="Terminal exocyclic allene axis has no source-assigned descriptor; toolkit axial support is limited",
    identity_note=("Candidate connectivity for the phenyl-substituted exocyclic allene 1l. "
                   "The flat drawing does not assign the potential tetrahedral/ring stereochemistry detected by graph inspection."),
)
m4la = new_mol(
    f"{DOC}:mol_4la", "4la (oxidized two-step product)", [E[7]],
    chirality="unspecified", geometry="unresolved", representation="not_encoded",
    assignment="Displayed alcohol product of Conditions B followed by NaBO3·4H2O; no E/Z descriptor is printed",
)
m5la = new_mol(
    f"{DOC}:mol_5la_regioisomer", "5la regioisomer fraction", [E[7]],
    chirality="unresolved", geometry="unresolved", representation="not_encoded",
    identity_note="Minor regioisomer implied by the reported 4:5 metric; its oxidized structure is not drawn on page 7.",
)
m1m = new_mol(
    f"{DOC}:mol_1m", "1m", [E[7]],
    chirality="specified", geometry="not_applicable", representation="not_encoded",
    assignment="Steroid substrate relative/absolute stereochemical relationships are depicted by wedge/dash bonds; no R/S text is printed",
    axial="Terminal exocyclic allene has no source-assigned axial descriptor",
)


def steroid_product(code):
    e_mid = new_mol(
        f"{DOC}:mol_E_{code}", f"(E)-{code}", [E[7]],
        chirality="specified", geometry="specified", representation="not_encoded",
        assignment=("Steroid stereochemistry is depicted; " +
                    ("the displayed structure was confirmed by single-crystal X-ray diffraction" if code == "4mb" else "no absolute R/S names are printed")),
        double_bonds="E alkene geometry explicitly labelled and drawn",
    )
    z_mid = new_mol(
        f"{DOC}:mol_Z_{code}", f"(Z)-{code} minor fraction", [E[7]],
        chirality="partially_specified", geometry="specified", representation="not_encoded",
        assignment="Minor stereoisomer is implied by the E/Z result; it is not independently drawn or assigned",
        double_bonds="Z geometry is the second component explicitly implied by the E/Z metric",
    )
    reg_mid = new_mol(
        f"{DOC}:mol_5{code[1:]}_regioisomer", f"5{code[1:]} regioisomer fraction", [E[7]],
        chirality="unresolved", geometry="unresolved", representation="not_encoded",
        identity_note="Minor regioisomer fraction implied by the 4:5 metric; its structure is not drawn on page 7.",
    )
    return e_mid, z_mid, reg_mid


e4ma, z4ma, m5ma = steroid_product("4ma")
e4mb, z4mb, m5mb = steroid_product("4mb")
m9e = new_mol(
    f"{DOC}:mol_9e_major", "9e major diastereomer", [E[7], T[7]],
    chirality="specified", geometry="not_applicable", representation="not_encoded",
    assignment=("Depicted homoallylic alcohol stereochemistry determined by single-crystal X-ray diffraction; "
                "the authors describe five contiguous stereocenters including two vicinal quaternary carbons"),
)
m9e_minor = new_mol(
    f"{DOC}:mol_9e_minor", "9e minor diastereomer", [E[7], T[7]],
    chirality="unresolved", geometry="not_applicable", representation="not_encoded",
    assignment="Minor diastereomer is implied by dr >95:5 but is not individually drawn or assigned",
)


# Figure 5A: two-step model exo-allene study.
r4la = new_rx(
    "figure5A_exo_allene_4la", DIRECT7 + COND_B,
    categories=["substrate_scope"], location="Figure 5A",
    name="Model study: borylation of an exo-allene",
    transformation="Conditions B carboboration followed by oxidative conversion of the boryl group",
    unresolved=[
        "NaBO3·4H2O charge, oxidation solvent, temperature, time, workup, purification, and yield measurement basis are not printed on page 7.",
        "Complete product/regioisomer structures and stereochemical assignments are not safely encoded.",
    ],
)
condition_b(r4la, m1l, f"{DOC}:mol_2a")
addp(r4la, "oxidant", f"{DOC}:mol_NaBO3_4H2O", "reagent", [E[7]], stage="workup")
addp(r4la, "product", m4la, "product", [E[7]], stage="workup", product_role="displayed oxidized two-step product")
addp(r4la, "regioisomer", m5la, "product", [E[7]], stage="workup", product_role="minor regioisomer fraction")
add_step(
    r4la, "oxidation", [E[7]], stage="workup",
    operation="then NaBO3·4H2O oxidation; other conditions not reported",
    participant_ids=["oxidant", "product", "regioisomer"],
)
outcome(r4la, "two_step_yield", "yield", [m4la], [E[7]], "69% (2 steps)", q(69, "%"),
        role="4la oxidized product", yield_type="reported two-step yield",
        scope="overall yield from 1l to displayed oxidized 4la; not a last-step-only yield")
outcome(r4la, "regioselectivity", "selectivity", [m4la, m5la], [E[7]], ">95:5", ">95:5",
        scope="4la:5la regioisomer ratio after the two-step sequence")
note(r4la, "Figure 5A explicitly says 'Conditions B'; page-5 Table 2 footnote c supplies the inherited first-stage charges, DMF, 30 °C, and 24 h.", [E[7], E[5]])
classify(r4la, [E[7]], OXIDATIVE_ISOLATION)


# Figure 5B: two separate steroid alkylboration experiments.
for code, halide, y, products in [
    ("4ma", f"{DOC}:mol_2a", 90, (e4ma, z4ma, m5ma)),
    ("4mb", f"{DOC}:mol_2b", 78, (e4mb, z4mb, m5mb)),
]:
    e_mid, z_mid, reg_mid = products
    r = new_rx(
        f"figure5B_steroid_{code}", DIRECT7 + COND_B,
        categories=["substrate_scope"], location=f"Figure 5B, (E)-{code}",
        name="Alkylboration of a steroid substrate",
        transformation="Conditions B alkylboration of steroid exocyclic allene 1m",
        unresolved=[
            "Workup, purification, and yield/selectivity measurement bases are not printed on page 7.",
            "The complex steroid structures and their complete absolute configurations are not safely encoded as SMILES.",
        ],
    )
    condition_b(r, m1m, halide)
    addp(r, "product_E", e_mid, "product", [E[7]], product_role="major E allylic boronate")
    addp(r, "product_Z", z_mid, "product", [E[7]], product_role="minor Z allylic-boronate fraction")
    addp(r, "regioisomer", reg_mid, "product", [E[7]], product_role="minor regioisomer fraction")
    r["steps"][0]["participant_ids"].extend(["product_E", "product_Z", "regioisomer"])
    outcome(r, "yield", "yield", [e_mid, z_mid], [E[7]], f"{y}%", q(y, "%"),
            role=f"{code} E/Z product fraction", yield_type="reported yield",
            scope=f"yield printed for (E)-{code}")
    outcome(r, "E_Z", "selectivity", [e_mid, z_mid], [E[7]], ">95:5", ">95:5",
            scope=f"E:Z ratio of {code}; target order E:Z")
    outcome(r, "4_5", "selectivity", [e_mid, z_mid, reg_mid], [E[7]], ">95:5", ">95:5",
            scope=f"{code}:5{code[1:]} regioisomer ratio")
    if code == "4mb":
        note(r, "The displayed (E)-4mb structure is paired with a single-crystal X-ray model on page 7.", [E[7]])
    note(r, "Figure 5B explicitly says 'Conditions B'; page-5 Table 2 footnote c supplies the inherited charges, DMF, 30 °C, and 24 h.", [E[7], E[5]])
    classify(r, [E[7]], CARBO)


# Figure 5C: allylboration of formaldehyde using (E)-4ma.
r9e = new_rx(
    "figure5C_formaldehyde_9e", DIRECT7,
    categories=["derivatization"], location="Figure 5C",
    name="Allylboration of a steroid substrate",
    transformation="allylboration of formaldehyde with (E)-4ma",
    unresolved=[
        "Reaction scale, workup, purification, and yield/dr measurement bases are not printed on page 7.",
        "The minor diastereomer is not individually drawn or assigned.",
    ],
)
addp(r9e, "allylic_boronate", e4ma, "reactant", [E[7]])
addp(r9e, "paraformaldehyde", f"{DOC}:mol_paraformaldehyde", "reactant", [E[7]],
     equivalents=("5 equiv", q(5, "equiv")), reference=e4ma)
addp(r9e, "solvent", f"{DOC}:mol_THF", "solvent", [E[7]])
addp(r9e, "major_product", m9e, "product", [E[7], T[7]], product_role="major homoallylic-alcohol diastereomer")
addp(r9e, "minor_product", m9e_minor, "product", [E[7], T[7]], product_role="minor diastereomer implied by dr")
r9e["solvent_conditions"]["composition"] = f("THF", ["THF"], [E[7]])
add_step(r9e, "reaction", [E[7]], operation="allylboration with (HCHO)n",
         temperature=("65 °C", q(65, "°C")), time=("24 h", q(24, "h")),
         participant_ids=[p["participant_id"] for p in r9e["participants"]])
outcome(r9e, "yield", "yield", [m9e, m9e_minor], [E[7], T[7]], "84%", q(84, "%"),
        role="9e diastereomeric product fraction", yield_type="reported yield", scope="yield of 9e")
outcome(r9e, "dr", "selectivity", [m9e, m9e_minor], [E[7], T[7]], ">95:5", ">95:5",
        scope="diastereomeric ratio of 9e; depicted/X-ray-assigned major:unassigned minor")
note(r9e, "The source states that 9e bears five contiguous stereocenters, including two vicinal quaternary carbon atoms, and that its depicted stereochemistry was determined by single-crystal X-ray diffraction.", [T[7], E[7]])
classify(r9e, [E[7]], ALLYLBORATION)


# Figure 6A: author-proposed cycle, not an observed synthesis.
for mid, label in [
    (f"{DOC}:mol_Int1", "copper(I) alkoxide Int1"),
    (f"{DOC}:mol_Int2", "boryl copper(I) intermediate Int2"),
    (f"{DOC}:mol_Int3", "alkenyl copper(I) species Int3"),
    (f"{DOC}:mol_Int4", "cuprate species Int4"),
    (f"{DOC}:mol_Int5", "proposed organocopper(III) intermediate Int5"),
]:
    new_mol(mid, label, [E[8], T[7]], chirality="unresolved", geometry="unresolved",
            representation="not_encoded",
            identity_note="Generic proposed catalytic intermediate; R groups/ligand are not expanded into a concrete molecule.")

for mid, label in [
    (f"{DOC}:mol_generic_gem_dialkylallene", "generic allene 1x"),
    (f"{DOC}:mol_generic_alkyl_halide", "generic alkyl halide 2y"),
    (f"{DOC}:mol_generic_allylic_boronate", "generic (E)-4xy product"),
    (f"{DOC}:mol_B2pin2", "diboron reagent 3"),
]:
    add_occurrence(mid, label, [E[8]])

mech = new_rx(
    "figure6A_proposed_catalytic_cycle", [E[8], T[7]],
    record_type="proposed_mechanism", categories=["mechanistic_study"],
    location="Figure 6A and adjoining page-7 discussion",
    name="Proposed reaction mechanism",
    transformation="proposed copper(I)-catalyzed alkylboration cycle",
    unresolved=[
        "The detailed alkylation mechanism is explicitly described as unclear; Int5 may be transient or a transition state of a concerted pathway.",
        "The cited stoichiometric mechanistic experiments are said to be in Supporting Information, which was not supplied or associated.",
    ],
)
for pid, mid, role in [
    ("allene", f"{DOC}:mol_generic_gem_dialkylallene", "reactant"),
    ("halide", f"{DOC}:mol_generic_alkyl_halide", "reactant"),
    ("diboron", f"{DOC}:mol_B2pin2", "reagent"),
    ("Int1", f"{DOC}:mol_Int1", "unknown"),
    ("Int2", f"{DOC}:mol_Int2", "unknown"),
    ("Int3", f"{DOC}:mol_Int3", "unknown"),
    ("Int4", f"{DOC}:mol_Int4", "unknown"),
    ("Int5", f"{DOC}:mol_Int5", "unknown"),
    ("product", f"{DOC}:mol_generic_allylic_boronate", "product"),
]:
    addp(mech, pid, mid, role, [E[8], T[7]], product_role="proposed product" if role == "product" else None)
for sid, operation, pids in [
    ("sigma_bond_metathesis", "proposed sigma-bond metathesis: Int1 to Int2", ["Int1", "diboron", "Int2"]),
    ("coordination_borylcupration", "proposed allene coordination and borylcupration: Int2 to Int3", ["Int2", "allene", "Int3"]),
    ("alkoxide_coordination", "proposed alkoxide coordination: Int3 to Int4", ["Int3", "Int4"]),
    ("oxidative_addition", "proposed SN2-type oxidative addition: Int4 to Int5", ["Int4", "halide", "Int5"]),
    ("reductive_elimination", "proposed reductive elimination to (E)-4xy and Int1", ["Int5", "product", "Int1"]),
]:
    add_step(mech, sid, [E[8], T[7]], operation=operation, participant_ids=pids)
note(mech, "The catalytic cycle is proposed from prior literature and unavailable SI mechanistic work; it is not represented as an independently observed reaction experiment on pages 7–12.", [E[8], T[7]])
classify(mech, [E[8], T[7]], coverage="unresolved",
         limitation="Author-proposed catalytic cycle with generic R groups and unresolved intermediate identities; no observed net experiment or complete atom mapping is asserted.")


# Page 9 repeats two existing optimization experiments with more granular
# three-component experimental distributions.  Merge the evidence/outcomes;
# do not create duplicate experiments.
for rid_suffix, raw, normalized, targets, scope in [
    ("table1_entry_14", "(E)-4aa/(Z)-4aa/5aa = 93:<1:7 (expt.)", "93:<1:7",
     [f"{DOC}:mol_E_4aa", f"{DOC}:mol_Z_4aa", f"{DOC}:mol_5aa_regioisomer"],
     "experimental three-component distribution in Figure 7; target order (E)-4aa:(Z)-4aa:5aa"),
    ("table1_entry_18", "(E)-4ba/(Z)-4ba/5ba = >99:<1:<1 (expt.)", ">99:<1:<1",
     [f"{DOC}:mol_E_4ba", f"{DOC}:mol_Z_4ba", f"{DOC}:mol_5ba_regioisomer"],
     "experimental three-component distribution in Figure 7; target order (E)-4ba:(Z)-4ba:5ba"),
]:
    rid = f"{DOC}:{rid_suffix}"
    r = reactions[rid]
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + [E[9], T[9]]))
    if not any(o["outcome_id"] == "figure7_three_component_distribution" for o in r["outcomes"]):
        o = outcome(r, "figure7_three_component_distribution", "selectivity", targets,
                    [E[9], E[4]], raw, normalized, scope=scope)
        o["measurement_basis"] = f(
            {"Figure 7": "experimental distribution", "Table 1 footnotes": "13C NMR for E/Z; GC for 4/5"},
            ["13C NMR analysis of roughly purified material for E/Z",
             "GC analysis of reaction mixture with internal standard for 4/5"],
            [E[9], E[4]],
        )
    note(r, "Figure 7 repeats this Table 1 experiment and reports a more granular E-major:Z-minor:regioisomer distribution; it is merged as duplicate evidence, not counted as a new experiment.", [E[9], E[4]])
    for mid in targets:
        label = next(m["occurrences"][0]["compound_label"]["normalized_value"]
                     for m in [molecules[mid]])
        add_occurrence(mid, label, [E[9]])


# Non-synthetic/computational material remains source context.
new_context = [
    {
        "context_id": f"{DOC}:figure6B_DFT_energy_profile",
        "document_id": DOC,
        "kind": "computational_mechanistic_study",
        "description": "DFT Gibbs-energy profile for coordination (TS1) and borylcupration (TS2) of simplified SIMesCuB(pin) with allene 1a; values are kcal/mol relative to EQ1.",
        "pathways": {
            "alkenyl_Z": {"TS1": 16.0, "EQ2": 9.4, "TS2": 20.1, "Int3": -26.7},
            "alkenyl_E": {"TS1": 14.7, "EQ2": 9.6, "TS2": 19.4, "Int3": -27.6},
            "allyl_E": {"TS1": 17.0, "EQ2": 12.4, "TS2": 16.4, "Int3": -25.3},
            "allyl_Z": {"TS1": 18.9, "EQ2": 13.0, "TS2": 16.5, "Int3": -28.0},
        },
        "evidence_ids": [E[8], T[7], T[8]],
    },
    {
        "context_id": f"{DOC}:figure6C_DFT_transition_state_comparison",
        "document_id": DOC,
        "kind": "computational_mechanistic_study",
        "description": "Calculated relative thermochemistry of four selectivity-determining transition states and predicted selectivities.",
        "relative_kcal_per_mol": {
            "TS1_allyl_E": {"ddG": 0.00, "ddH": 0.00, "minus_TddS": 0.00},
            "TS1_allyl_Z": {"ddG": 1.91, "ddH": 1.08, "minus_TddS": 0.83},
            "TS2_alkenyl_Z": {"ddG": 3.10, "ddH": 1.93, "minus_TddS": 1.17},
            "TS2_alkenyl_E": {"ddG": 2.44, "ddH": 1.20, "minus_TddS": 1.24},
        },
        "predicted_selectivity": {"E_Z": "97:3", "4_5": "99:1"},
        "comparison_experimental_4aa": {"E_Z": ">95:5", "4_5": "93:7"},
        "evidence_ids": [E[8], T[8], T[9]],
    },
    {
        "context_id": f"{DOC}:figure7_ligand_steric_maps",
        "document_id": DOC,
        "kind": "mechanistic_ligand_comparison",
        "description": "Steric maps compare current-work SIMes/Xantphos panels I–II with literature ligands dppf, IPr, L1, DTBM-dppbz, DTBM-SEGPHOS, and Josiphos in panels III–VIII; literature panels are not new experiments in this paper.",
        "evidence_ids": [E[9], T[9], T[10]],
    },
    {
        "context_id": f"{DOC}:page10_associated_content",
        "document_id": DOC,
        "kind": "associated_content_notice",
        "description": "The article states that experimental/characterization/computational SI and calculated structures exist and lists CCDC 2062673, 2067460, 2067461, and 2070912; no SI/CIF was supplied or associated in this workspace.",
        "evidence_ids": [E[10], T[10]],
    },
]
existing_context_ids = {c.get("context_id") for c in bundle["source_context"]}
bundle["source_context"].extend(c for c in new_context if c["context_id"] not in existing_context_ids)


bundle["molecules"] = list(molecules.values())
bundle["reactions"] = list(reactions.values())
st.write(RUN / "submitted.json", bundle)


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "Figure 5A", "finding": "One two-step model exo-allene experiment: Conditions B then NaBO3·4H2O, 4la 69% over two steps, 4:5 >95:5.", "evidence_ids": [E[7], T[7]]},
                {"region": "Figure 5B", "finding": "Two separate steroid exo-allene alkylborations recorded: (E)-4ma 90% and (E)-4mb 78%, each E/Z >95:5 and 4:5 >95:5; 4mb has an X-ray model.", "evidence_ids": [E[7], T[7]]},
                {"region": "Figure 5C and adjacent prose", "finding": "One formaldehyde allylboration of (E)-4ma recorded: 9e 84%, dr >95:5, with depicted stereochemistry assigned by X-ray.", "evidence_ids": [E[7], T[7]]},
                {"region": "right column and lower-left mechanism discussion", "finding": "Proposed catalytic cycle discussion begins; recorded with Figure 6A as a proposed mechanism, not as an observed experiment.", "evidence_ids": [E[7], T[7], E[8]]},
            ],
            "record_ids": [r4la["reaction_id"], f"{DOC}:figure5B_steroid_4ma", f"{DOC}:figure5B_steroid_4mb", r9e["reaction_id"], mech["reaction_id"]],
            "unresolved": [
                "Figure 5 yields/selectivities do not print analytical measurement bases.",
                "Complex steroid structures and full absolute configurations were not safely encoded as SMILES.",
                "Figure 5A oxidation charge/solvent/temperature/time and all workup/purification details are unavailable without SI.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "Figure 6A", "finding": "Author-proposed catalytic cycle Int1–Int5; recorded as proposed_mechanism. Int5 is explicitly uncertain/transient versus concerted transition state.", "evidence_ids": [E[8], T[7]]},
                {"region": "Figure 6B", "finding": "DFT Gibbs-energy profiles for four regio/stereochemical paths; retained as computational source_context, not a synthetic experiment.", "evidence_ids": [E[8], T[7], T[8]]},
                {"region": "Figure 6C", "finding": "Transition-state structures and relative thermochemistry; retained as computational source_context.", "evidence_ids": [E[8], T[8], T[9]]},
                {"region": "bottom continuation prose", "finding": "Interpretation of selectivity-determining steps; no independent wet-lab reaction arrow or control on this page.", "evidence_ids": [E[8], T[8]]},
            ],
            "record_ids": [mech["reaction_id"]],
            "unresolved": ["The paper refers to Xantphos/allene-1b calculations and stoichiometric mechanistic details in SI, but no SI is associated."],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "Figure 7 panels I–II", "finding": "Current-work SIMes and Xantphos steric maps repeat Table 1 entry 14/18 experiments with granular experimental distributions; merged as duplicate evidence.", "evidence_ids": [E[9], E[4], T[9]]},
                {"region": "Figure 7 panels III–VIII", "finding": "Previously reported ligand/reaction classes from Liu, Gagosz/Riant, Procter, Fujihara/Tsuji, and Hoveyda; retained as literature/mechanistic context, not new experiments.", "evidence_ids": [E[9], T[9], T[10]]},
                {"region": "lower DFT prose", "finding": "Computed relative Gibbs/enthalpy/entropy contributions and theoretical-vs-experimental selectivity comparison; retained as source_context.", "evidence_ids": [E[9], T[9]]},
            ],
            "record_ids": [f"{DOC}:table1_entry_14", f"{DOC}:table1_entry_18"],
            "unresolved": ["Figure 7 combines E/Z and regioisomer distributions whose component analytical methods are defined separately in Table 1 footnotes."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "left-column Figure 7 discussion", "finding": "Mechanistic interpretation plus cited prior-work ligand comparisons; no new experiment row or arrow.", "evidence_ids": [E[10], T[10]]},
                {"region": "Conclusions", "finding": "Summary of the reported scope and allylboration utility; no new experiment.", "evidence_ids": [E[10], T[10]]},
                {"region": "Associated Content and Accession Codes", "finding": "SI availability and four CCDC codes are stated, but no SI/CIF file was supplied or associated; retained as source_context only.", "evidence_ids": [E[10], T[10]]},
                {"region": "author information, acknowledgments, references start", "finding": "Bibliographic/nonexperimental matter.", "evidence_ids": [E[10], T[10]]},
            ],
            "record_ids": [],
            "unresolved": ["External SI and crystallographic files were not provided; their experimental details are not inferred."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "References 2–40", "finding": "Reference list only; cited publications are not current-paper experiments and no synthesis record is created.", "evidence_ids": [E[11], T[11]]},
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 12,
            "state": "attempted",
            "regions": [
                {"region": "References 40–77 and footnote 72", "finding": "Reference list plus mechanistic caveat pointing to unavailable SI; cited work is not recast as current-paper experiments.", "evidence_ids": [E[12], T[12]]},
            ],
            "record_ids": [mech["reaction_id"]],
            "unresolved": ["Footnote 72 mentions stoichiometric mechanistic studies in SI, but the SI is not associated and no separate experiment can be extracted."],
        },
    ],
}
st.write(RUN / "coverage/pages_0007_0012.json", coverage)

visual = {
    "document_id": DOC,
    "pages_inspected": [7, 8, 9, 10, 11, 12],
    "page_images": [f"assets/{DOC}/page{p}.jpg" for p in range(7, 13)],
    "figure5_arrows": [
        {"panel": "5A", "record_id": r4la["reaction_id"], "result": "4la 69% (2 steps), 4:5 >95:5"},
        {"panel": "5B", "record_id": f"{DOC}:figure5B_steroid_4ma", "result": "(E)-4ma 90%, E/Z >95:5, 4:5 >95:5"},
        {"panel": "5B", "record_id": f"{DOC}:figure5B_steroid_4mb", "result": "(E)-4mb 78%, E/Z >95:5, 4:5 >95:5"},
        {"panel": "5C", "record_id": r9e["reaction_id"], "result": "9e 84%, dr >95:5"},
    ],
    "proposed_mechanism": mech["reaction_id"],
    "duplicate_evidence": [f"{DOC}:table1_entry_14", f"{DOC}:table1_entry_18"],
    "non_experiment_pages": {"10": "discussion/conclusion/associated content", "11": "references", "12": "references and mechanistic footnote"},
}
st.write(RUN / "raw_model/visual_transcription_pages7_12.json", visual)

counts = {
    "molecules_total_in_submission": len(bundle["molecules"]),
    "reactions_total_in_submission": len(bundle["reactions"]),
    "new_reactions_pages7_12": 5,
    "new_specific_experiments_pages7_12": 4,
    "new_proposed_mechanisms_pages7_12": 1,
    "updated_existing_reactions": 2,
    "all_reaction_types": dict(Counter(r["record_type"] for r in bundle["reactions"])),
    "needs_review": sum(r["review"]["status"] == "needs_review" for r in bundle["reactions"]),
}
st.write(RUN / "raw_model/premerge_counts_pages7_12.json", counts)

marker = "## Extension: original PDF pages 7–12"
workflow_path = ROOT / "runs/workflow.md"
workflow = workflow_path.read_text() if workflow_path.exists() else "# Extraction workflow\n"
if marker not in workflow:
    workflow += f"""

{marker}

- Batch: `allene-stereo-v31-200-20261017`; document `{DOC}`; claimed visual scope is exactly original PDF pages 7–12. All six full-page JPEGs were opened at original detail. Page 7 Figure 5 and pages 8–9 Figures 6–7 were additionally cropped with lineage receipts under `runs/current/raw_model/`.
- `runs/current/source_mapping.json` identifies the main article and explicitly leaves SI association unestablished. Page 10 states that SI exists, but no SI was supplied or inferred.
- Figure 5A is one two-step experiment, not two yields: Conditions B carboboration followed by NaBO3·4H2O oxidation, with 69% overall yield from 1l to the displayed alcohol and 4:5 >95:5. The oxidation charge/solvent/time/temperature remain unreported.
- Figure 5B contains two independent experiments (4ma and 4mb). `Conditions B` is inherited only from the explicit Table 2 definition on page 5; direct page-7 values (2.0 equiv halide, yields and selectivities) retain page-7 evidence. Yield/selectivity measurement bases are not printed and remain missing.
- Figure 5C is a separate formaldehyde allylboration. The major 9e structure and X-ray assignment are distinct from the unassigned minor diastereomer implied by `dr >95:5`.
- The complex steroid structures retain wedge/dash and E/Z/X-ray source assignments as textual stereochemistry with `representation_status=not_encoded`; no unsafe full SMILES is invented. For 1l, the legible connectivity candidate was checked with `standard stereo`; a potential unassigned tetrahedral element remains, so chirality/geometry are not labelled achiral or fully encoded.
- Figure 6A is an author-proposed mechanism. The DFT energy profiles/transition states in Figure 6B/C are computational source context, not wet-lab synthesis records. Figure 7 panels III–VIII are cited prior work, not current-paper experiments.
- Figure 7 panels I–II repeat Table 1 entry 14 and entry 18. Their granular `(E):(Z):regioisomer` distributions were merged into those existing reaction IDs as new evidence/outcomes; they are not counted again.
- Pages 10–12 contain no new current-paper synthesis arrows: page 10 is mechanism discussion/conclusions/associated-content metadata; pages 11–12 are references. Footnote 72 points to unavailable SI mechanistic work and therefore cannot support a separate extracted experiment.
- All new/updated reactions remain `needs_review`; semantic and chemical accuracy flags remain false. Schema validation and audit are mechanical checks, not source verification.
"""
    workflow_path.write_text(workflow)


def finalize():
    merged = json.loads((RUN / "dataset.json").read_text())
    audit = dataset_audit(merged)
    st.write(RUN / "audit.json", audit)
    standard = merged.get("standard", {})
    final_counts = {
        "documents": len(merged.get("documents", [])),
        "molecules": len(standard.get("molecules", [])),
        "reactions": len(standard.get("reactions", [])),
        "new_reactions_pages7_12": 5,
        "updated_existing_reactions": 2,
        "reaction_types": dict(Counter(r["record_type"] for r in standard.get("reactions", []))),
        "needs_review": sum(r["review"]["status"] == "needs_review" for r in standard.get("reactions", [])),
        "warning_records": len(audit.get("warnings", [])),
        "warning_categories": dict(Counter(w.split(":", 1)[0] for row in audit.get("warnings", []) for w in row["warnings"])),
        "semantic_accuracy_verified": False,
        "chemical_accuracy_verified": False,
    }
    st.write(RUN / "raw_model/final_counts_pages7_12.json", final_counts)
    print(json.dumps(final_counts, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--finalize", action="store_true")
    args = parser.parse_args()
    if args.finalize:
        finalize()
    else:
        print(json.dumps(counts, ensure_ascii=False, indent=2))
