#!/usr/bin/env python3
"""Build the evidence-linked v3.1 submission for doc_dc20f9448ea23b61.

This is a source-specific extraction artifact, not a reusable chemistry rule.
All chemical structures are candidates derived from the supplied page drawings;
the merge step deliberately quarantines the records as needs_review.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_dc20f9448ea23b61"
RUN = ROOT / "runs/current"
E = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
T = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, verification=None, derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   verification=verification, derivation=derivation,
                   alternatives=alternatives)


def inferred(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="inference", derivation=derivation)


def derived(normalized, evidence, derivation):
    return f(None, normalized, evidence, method="structure_derivation", derivation=derivation)


def not_applicable():
    return f(method="not_applicable", status="not_applicable")


def set_status(m, key, normalized, evidence, derivation):
    m["stereochemistry"][key] = inferred(normalized, evidence, derivation)


def molecule(mid, label, evidence, *, name=None, smiles=None, stereo="unresolved",
             identity_note=None):
    ev = list(evidence)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    if name:
        m["names"]["english"] = f(name, name, ev)
    m["cas_numbers"] = f(status="unreported")
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, ev,
            "Candidate connectivity transcribed from the supplied source drawing or standard identity; RDKit parsing is only a graph check, not source certification.",
        )
        diag = st.stereo_inspect(smiles)
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "RDKit stereo_graph_key calculated from the candidate SMILES; it is not a sample identity or source verification.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f(status="unreported")

    basis = "Status assigned from the visibly inspected source drawing; absence of a named descriptor is not interpreted as achiral or racemic."
    if stereo == "achiral":
        set_status(m, "chirality_status", "achiral", ev, "Candidate graph has no potential stereogenic element; source identity remains machine-extracted.")
        set_status(m, "geometry_status", "not_applicable", ev, "No independently configurable E/Z element in the candidate graph.")
        set_status(m, "representation_status", "encoded", ev, "Candidate SMILES fully represents the graph and absence of stereogenic elements.")
        m["stereochemistry"]["isomeric_smiles"] = derived(st.stereo_inspect(smiles)["canonical_isomeric_smiles"], ev, "RDKit canonical isomeric SMILES from the candidate graph.")
    elif stereo == "alkene_E":
        set_status(m, "chirality_status", "achiral", ev, "The inspected candidate has no tetrahedral or axial stereogenic element.")
        set_status(m, "geometry_status", "specified", ev, "The source drawing shows the two non-hydrogen alkene substituents trans (E connectivity candidate).")
        set_status(m, "representation_status", "encoded", ev, "The E alkene geometry is encoded with directional SMILES bonds.")
        m["stereochemistry"]["double_bonds"] = f(
            {"location": "source alkene drawing", "configuration": "E"},
            [{"configuration": "E", "assignment_basis": "source drawing"}], ev)
        m["stereochemistry"]["isomeric_smiles"] = derived(st.stereo_inspect(smiles)["canonical_isomeric_smiles"], ev, "RDKit canonicalized candidate retaining the trans alkene geometry.")
    elif stereo == "allene_unassigned":
        set_status(m, "chirality_status", "unspecified", ev, basis)
        set_status(m, "geometry_status", "not_applicable", ev, "The stereogenic element is an allene axis rather than an E/Z alkene assignment.")
        set_status(m, "representation_status", "not_encoded", ev, "The source does not assign the potentially stereogenic allene axis; no axial descriptor was invented.")
        m["stereochemistry"]["axial"] = f(
            {"location": "allenoate C=C=C unit", "descriptor": "not reported"},
            [{"axis": "allene", "configuration": "unspecified"}], ev)
    elif stereo in {"product_assigned", "product_partial"}:
        status = "specified" if stereo == "product_assigned" else "partially_specified"
        set_status(m, "chirality_status", status, ev, basis)
        set_status(m, "geometry_status", "not_applicable", ev, "The endocyclic alkene geometry is constrained by the five-membered ring and no separate E/Z descriptor is reported.")
        set_status(m, "representation_status", "not_encoded", ev, "Wedge/dash stereochemistry in the source product drawing is retained as an assignment note but not encoded as arbitrary @/@@ tags.")
        m["stereochemistry"]["atom_centers"] = f(
            {"location": "product drawing", "feature": "wedge/dash bonds at the two newly formed pyrroline stereocenters"},
            [{"centers": "two pyrroline ring stereocenters", "configuration": "source-drawn; R/S labels not transcribed"}], ev)
        m["stereochemistry"]["assignment"] = f(
            "absolute configurations of products 3 were confirmed by analogy to X-ray analysis of 3ar",
            {"type": "absolute assignment by analogy", "reference_product": "3ar", "explicit_RS_labels": None},
            [E[3], T[3]],
        )
        if stereo == "product_partial":
            m["stereochemistry"]["mixture"] = f(
                {"location": "Figure 1", "feature": "reported dr or pre-existing stereocenters"},
                {"type": "diastereomeric or stereochemically complex sample", "composition": "experiment-specific; see selectivity outcomes"}, ev)
    elif stereo == "chiral_catalyst":
        set_status(m, "chirality_status", "specified", ev, basis)
        set_status(m, "geometry_status", "not_applicable", ev, "No source E/Z assignment is relevant to the drawn catalyst.")
        set_status(m, "representation_status", "not_encoded", ev, "Dense bicyclic catalyst stereochemistry was not converted to a speculative SMILES.")
        m["stereochemistry"]["assignment"] = f(
            {"location": "Table 1 catalyst drawing", "feature": "explicit wedge/dash stereochemistry"},
            {"type": "source-drawn configuration", "descriptor": None}, ev)
    else:
        set_status(m, "chirality_status", "unresolved", ev, basis)
        set_status(m, "geometry_status", "unresolved", ev, basis)
        set_status(m, "representation_status", "unresolved" if not smiles else "partial", ev, "The available drawing was not converted into a complete stereochemical representation.")

    if identity_note:
        m["identity_note"] = f(identity_note, identity_note, ev)
    elif not smiles:
        m["identity_note"] = f(
            "Source-labelled molecule retained without a complete SMILES because the dense or variable drawing was not safely transcribed.",
            "Source-labelled identity; complete structure unresolved", ev)
    return m


def template_molecule(mid, label, evidence, fixed, variables, mappings):
    ev = list(evidence)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, ev)
    m["names"]["original"] = f(label, label, ev)
    m["representation_type"] = "generic_template"
    m["smiles"] = f(status="unreported")
    m["generic_template"] = {
        "fixed_structure": f({"location": label}, fixed, ev),
        "variables": [],
        "instance_mapping": [],
        "reason": f("Generic R/Ar-labelled structure in the source scheme", "Variable source structure cannot be one concrete SMILES", ev),
    }
    for var in variables:
        m["generic_template"]["variables"].append({
            "label": var,
            "attachment_points": f(f"{var} attachment shown in source", {"label": var, "attachment": "source-drawn position"}, ev),
            "definition": f(f"{var} as drawn/defined in source", "source-defined substituent", ev),
            "allowed_values": f(status="unreported"),
            "constraints": f(status="unreported"),
        })
    for target, compound_label, assignments in mappings:
        m["generic_template"]["instance_mapping"].append({
            "molecule_id": target,
            "compound_label": f(compound_label, compound_label, ev),
            "assignments": f(assignments, assignments, ev),
        })
    set_status(m, "chirality_status", "unresolved", ev, "Template variables may change stereogenicity; no global assignment is possible.")
    set_status(m, "geometry_status", "unresolved", ev, "Template variables and source drawing do not define a single geometry state.")
    set_status(m, "representation_status", "not_encoded", ev, "Generic template intentionally is not represented as wildcard SMILES.")
    m["identity_note"] = f("Generic scheme template", "Not a single concrete substance", ev)
    return m


def participant(pid, mid, role, evidence, *, amount=None, volume=None, equiv=None,
                molpct=None, reference=None, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    if amount:
        p["charges"]["amount"] = f(amount[0], amount[1], evidence)
    if volume:
        p["charges"]["volume"] = f(volume[0], volume[1], evidence)
    if equiv:
        p["charges"]["equivalents"] = f(equiv[0], equiv[1], evidence)
    if molpct:
        p["charges"]["mol_percent"] = f(molpct[0], molpct[1], evidence)
    if reference:
        p["equivalence_reference"] = f(reference, reference, evidence)
    return p


def step(step_id, evidence, participant_ids, *, operation=None, temperature=None,
         time=None, stage="reaction"):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    if operation:
        s["operation"] = f(operation, operation, evidence)
    if temperature:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time:
        s["time"] = f(time[0], time[1], evidence)
    return s


def outcome(oid, kind, targets, evidence, *, raw, normalized, product_role=None,
            yield_type=None, basis=None, scope=None, status=None, alternatives=()):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = f(raw, normalized, evidence, status=status, alternatives=alternatives)
    if product_role:
        o["product_role"] = f(product_role, product_role, evidence)
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, evidence)
    if basis:
        o["measurement_basis"] = f(basis, basis, evidence)
    if scope:
        o["scope"] = f(scope, scope, evidence)
    return o


def classification_unresolved(evidence, reason):
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(status="unreported"),
        "bond_changes": [],
        "labels": f(status="unreported"),
        "coverage": inferred("unresolved", evidence, reason),
        "limitations": [reason, "No atom mapping or model votes were fabricated."],
        "model_votes": [],
    }


def classification_na(evidence, reason):
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(method="not_applicable", status="not_applicable"),
        "bond_changes": [],
        "labels": f([], [], evidence, method="inference", derivation=reason),
        "coverage": inferred("not_applicable", evidence, reason),
        "limitations": [reason],
        "model_votes": [],
    }


def classification_from_changes(evidence, changes, reason):
    labels = {}
    for ch in changes:
        if ch["change_type"] in {"formed", "hydrogen_gain", "order_increased"}:
            label = st.bond_label(ch["atoms"][0]["element"], ch["atoms"][1]["element"], ch["after_order"])
            labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(ch["change_id"])
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": f(status="unreported"),
        "bond_changes": [derived(ch, evidence, reason) for ch in changes],
        "labels": derived(list(labels.values()), evidence, reason),
        "coverage": inferred("partial", evidence, "Class-level net changes were read from the source structures; exact atom mapping, hydrogen provenance, and coproduct balance remain unresolved."),
        "limitations": [
            "No exact atom mapping is asserted.",
            "Implicit-H changes are net changes; the elementary proton-transfer sequence is not treated as experimentally proven.",
            "Unlisted coproducts and reagent-derived atoms are outside this partial classification.",
        ],
        "model_votes": [],
    }


def annulation_classification(evidence):
    changes = [
        {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "b3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        {"change_id": "b4", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        {"change_id": "b5", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
    ]
    return classification_from_changes(evidence, changes, "Partial graph comparison of the source-drawn tandem isomerization/[3+2] annulation product versus the allenoate and allylamine reactants.")


def reduction_classification(evidence):
    changes = []
    idx = 1
    for _ in range(2):
        changes += [
            {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
            {"change_id": f"b{idx+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 2},
            {"change_id": f"b{idx+2}", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ]
        idx += 3
    return classification_from_changes(evidence, changes, "Partial class-level comparison of the two ester groups in 3aa with the two alcohol groups in 4.")


def chlorination_classification(evidence):
    changes = []
    idx = 1
    for _ in range(2):
        changes += [
            {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Cl"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": f"b{idx+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
            {"change_id": f"b{idx+2}", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        ]
        idx += 3
    return classification_from_changes(evidence, changes, "Partial class-level comparison of the two alcohol groups in 4 with the two chloride substituents in 5.")


def oxidation_classification(evidence):
    changes = []
    idx = 1
    for _ in range(2):
        changes += [
            {"change_id": f"b{idx}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
            {"change_id": f"b{idx+1}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
            {"change_id": f"b{idx+2}", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
        ]
        idx += 3
    return classification_from_changes(evidence, changes, "Partial class-level comparison of the two alcohol groups in 4 with the two aldehydes in 6.")


def base_reaction(rid, evidence, record_type, categories, location, name, transformation):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence) if name else f(status="unreported")
    r["transformation_type"] = inferred(transformation, evidence, "Concise transformation description inferred from source reactants and product drawing.") if transformation else f(status="unreported")
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": ["No Supporting Information was supplied or explicitly associated; detailed workup, purification, and characterization remain unavailable."],
    }
    r["completeness"] = "partial"
    return r


def set_solvent(r, name, volume, evidence):
    r["solvent_conditions"]["composition"] = f(name, [name], evidence)
    if volume:
        r["solvent_conditions"]["total_volume"] = f(volume[0], volume[1], evidence)


def finish_procedure(r, raw, evidence, general=None):
    if raw:
        r["procedure"]["original_text"] = f(raw, raw, evidence)
    if general:
        r["procedure"]["general_procedure_reference"] = f(general, general, evidence)


def add_note(r, raw, normalized, evidence, method="original_extraction", derivation=None):
    r["notes"].append(f(raw, normalized, evidence, method=method, derivation=derivation))


def q(number, unit, **extra):
    return {"number": number, "unit": unit, **extra}


def main():
    data = json.loads((RUN / "dataset.json").read_text())
    if any(r.get("status") == "reviewed" for r in data.get("experiments", [])) or any(
        r.get("review", {}).get("status") == "source_reviewed"
        for r in data.get("standard", {}).get("reactions", [])
    ):
        raise SystemExit("Refusing to rebuild over source-reviewed records")
    doc = next(d for d in data["documents"] if d["document_id"] == DOC)
    doc["title"] = "Phosphine-Catalyzed Asymmetric Tandem Isomerization/Annulation of Allyl Amines with Allenoates: Enantioselective Annulation of a Saturated C–N Bond"
    doc["supplementary_material"] = {
        "status": "unknown",
        "note": "The article advertises Supporting Information on PDF page 4, but no SI file was supplied or explicitly associated in this session.",
    }
    data.setdefault("run_log", [])[:] = [row for row in data.get("run_log", []) if row.get("event") != "source_metadata_transcribed"]
    data["run_log"].append({
        "event": "source_metadata_transcribed",
        "document_id": DOC,
        "title_evidence": T[1],
        "note": "Title transcribed from the supplied main article; SI remains unassociated.",
    })
    (RUN / "dataset.json").write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
    bundle = st.empty(data["documents"])

    # Concrete structures safe enough to preserve as machine candidates.
    mols = []
    mols += [
        molecule(f"{DOC}:mol_1a", "1a", [E[2], E[3]], name="γ-isopropyl-substituted allenoate 1a", smiles="CCOC(=O)C=C=C(C(C)C)", stereo="allene_unassigned"),
        molecule(f"{DOC}:mol_2a", "2a", [E[2], E[3]], name="ethyl (E)-4-(tosylamino)but-2-enoate", smiles="CCOC(=O)/C=C/CNS(=O)(=O)c1ccc(C)cc1", stereo="alkene_E"),
        molecule(f"{DOC}:mol_2d", "2d", [E[3]], name="Boc-protected allylamine control 2d", smiles="CC(C)(C)OC(=O)NC/C=C/C(=O)OCC", stereo="alkene_E"),
        molecule(f"{DOC}:mol_2e", "2e", [E[3]], name="benzoyl-substituted N-tosyl allylamine control 2e", smiles="O=C(c1ccccc1)/C=C/CNS(=O)(=O)c1ccc(C)cc1", stereo="alkene_E"),
        molecule(f"{DOC}:mol_3aa", "3aa", [E[2], E[3], E[4]], name="pyrroline cycloadduct 3aa", smiles="CCOC(=O)C1=CC(C(C)C)N(S(=O)(=O)c2ccc(C)cc2)C1CCC(=O)OCC", stereo="product_assigned"),
        molecule(f"{DOC}:mol_4", "4", [E[3]], name="diol derivative 4", smiles="OCC1=CC(C(C)C)N(S(=O)(=O)c2ccc(C)cc2)C1CCCO", stereo="product_assigned"),
        molecule(f"{DOC}:mol_5", "5", [E[3]], name="dichloride derivative 5", smiles="ClCC1=CC(C(C)C)N(S(=O)(=O)c2ccc(C)cc2)C1CCCCl", stereo="product_assigned"),
        molecule(f"{DOC}:mol_6", "6", [E[3]], name="dialdehyde derivative 6", smiles="O=CC1=CC(C(C)C)N(S(=O)(=O)c2ccc(C)cc2)C1CCC=O", stereo="product_assigned"),
        molecule(f"{DOC}:mol_PPhMe2", "PPhMe2", [E[2]], name="dimethylphenylphosphine", smiles="CP(C)c1ccccc1", stereo="unresolved"),
        molecule(f"{DOC}:mol_toluene", "toluene", [E[1], E[2], E[3]], smiles="Cc1ccccc1", stereo="achiral"),
        molecule(f"{DOC}:mol_dcm", "CH2Cl2", [E[2], E[3]], name="dichloromethane", smiles="ClCCl", stereo="achiral"),
        molecule(f"{DOC}:mol_thf", "THF", [E[3]], name="tetrahydrofuran", smiles="C1CCOC1", stereo="achiral"),
        molecule(f"{DOC}:mol_lialh4", "LiAlH4", [E[3]], name="lithium aluminium hydride", smiles="[Li+].[AlH4-]", stereo="unresolved",
                 identity_note="Ionic candidate graph retained; toolkit tetrahedral diagnostics at aluminium are not interpreted as a source stereochemical assignment."),
        molecule(f"{DOC}:mol_triphosgene", "CO(OCCl3)2", [E[3]], name="triphosgene", smiles="O=C(OC(Cl)(Cl)Cl)OC(Cl)(Cl)Cl", stereo="achiral"),
        molecule(f"{DOC}:mol_et3n", "Et3N", [E[3]], name="triethylamine", smiles="CCN(CC)CC", stereo="achiral"),
        molecule(f"{DOC}:mol_d2o", "D2O", [E[3], E[4]], name="deuterium oxide", smiles="[2H]O[2H]", stereo="achiral"),
        molecule(f"{DOC}:mol_pcc", "PCC", [E[3]], name="pyridinium chlorochromate", stereo="unresolved"),
    ]

    # Chiral phosphines and source-labelled series members retained without speculative full structures.
    for label in [f"P{i}" for i in range(1, 10)]:
        mols.append(molecule(f"{DOC}:mol_{label}", label, [E[2]], name=f"source-drawn chiral phosphine {label}", stereo="chiral_catalyst"))
    for letter in "bcdefghijklmnopqrstuvwx":
        mols.append(molecule(f"{DOC}:mol_1{letter}", f"1{letter}", [E[3]], name=f"source-drawn γ-substituted allenoate 1{letter}", stereo="unresolved",
                             identity_note="Concrete allenoate identity is preserved by label and source drawing; the potentially stereogenic allene axis is not assigned and a complete SMILES was not safely transcribed."))
    for letter in "bc":
        mols.append(molecule(f"{DOC}:mol_2{letter}", f"2{letter}", [E[3]], name=f"source-drawn arylsulfonyl allylamine 2{letter}", stereo="unresolved"))
    mixture_products = {"3af", "3ag", "3ah", "3ai", "3at", "3au"}
    for suffix in [f"3a{x}" for x in "bcdefghijklmnopqrstuvwx"] + ["3ba", "3ca"]:
        mols.append(molecule(f"{DOC}:mol_{suffix}", suffix, [E[3]], name=f"source-drawn pyrroline product {suffix}",
                             stereo="product_partial" if suffix in mixture_products else "product_assigned"))
    mols.append(molecule(
        f"{DOC}:mol_3aa_D_mixture", "deuterated 3aa", [E[3], E[4]],
        name="site-partially deuterated 3aa isotopologue mixture", stereo="product_partial",
        identity_note="The source reports fractional incorporation at D1, D2 and D3; this is an isotopologue mixture, not a single fully deuterated molecule.",
    ))

    # Generic templates used only for the current general scheme.
    mols += [
        template_molecule(f"{DOC}:mol_template_1", "allenoate 1 (generic)", [E[1], E[2], E[3]],
                          {"core": "R-substituted allenoate", "fixed_feature": "C=C=C-CO2Et"}, ["R"],
                          [(f"{DOC}:mol_1a", "1a", {"R": "i-Pr"})]),
        template_molecule(f"{DOC}:mol_template_2", "allylamine 2 (generic)", [E[1], E[2], E[3]],
                          {"core": "ArSO2NH-allyl-CO2Et", "fixed_feature": "sulfonamido vinyl ester"}, ["Ar"],
                          [(f"{DOC}:mol_2a", "2a", {"ArSO2": "Ts"})]),
        template_molecule(f"{DOC}:mol_template_3", "pyrroline 3 (generic)", [E[1], E[2], E[3]],
                          {"core": "N-sulfonyl pyrroline bearing CO2Et, R, and CH2CH2CO2Et"}, ["R", "Ar"],
                          [(f"{DOC}:mol_3aa", "3aa", {"R": "i-Pr", "ArSO2": "Ts"})]),
    ]
    bundle["molecules"] = mols

    reactions = []
    general_id = f"{DOC}:general_tandem_annulation"
    r = base_reaction("general_tandem_annulation", [E[1], E[2], E[3]], "generic_scheme", ["other"],
                      "Graphical abstract; Scheme 1c; Figure 1 header", "asymmetric tandem isomerization/annulation", "phosphine-catalyzed tandem isomerization/[3+2] annulation")
    r["participants"] = [
        participant("generic_1", f"{DOC}:mol_template_1", "reactant", [E[1], E[2], E[3]]),
        participant("generic_2", f"{DOC}:mol_template_2", "reactant", [E[1], E[2], E[3]]),
        participant("generic_3", f"{DOC}:mol_template_3", "product", [E[1], E[2], E[3]], product_role="main product class"),
        participant("P6", f"{DOC}:mol_P6", "catalyst", [E[1]], molpct=("20 mol%", q(20, "mol%")), reference=f"{DOC}:mol_template_1"),
        participant("toluene", f"{DOC}:mol_toluene", "solvent", [E[1]]),
    ]
    r["steps"] = [step("general_reaction", [E[1]], [p["participant_id"] for p in r["participants"]], temperature=("25 °C", q(25, "°C")))]
    set_solvent(r, "toluene", None, [E[1]])
    add_note(r, "up to 87% yield, >99% ee", {"yield_max_percent": 87, "ee_lower_bound_percent": 99}, [E[1]])
    add_note(r, "PG = Protecting Group", "PG is a source-defined variable protecting group", [E[1]])
    r["classification"] = annulation_classification([E[1], E[2], E[3]])
    r["review"]["unresolved"].append("Generic R/Ar/PG variables and exact atom mapping remain unresolved; the graphical abstract is not counted as a concrete experiment.")
    reactions.append(r)

    literature_rows = [
        ("literature_well_developed_aromatic_imines", "Scheme 1a, left", "Prior literature: well-developed [3+2] annulation of allenes and aromatic imines", [E[2]]),
        ("literature_lu_2012_aliphatic_imines", "Scheme 1a, middle", "Prior literature: Lu 2012 asymmetric [3+2] annulation using DPP-protected aliphatic imines", [E[2], T[1]]),
        ("literature_terada_2008_relay", "Scheme 1b", "Prior literature: Terada relay-catalyzed isomerization/carbon–carbon bond formation of N-allylamines", [E[2], T[1]]),
    ]
    for rid, loc, note, ev in literature_rows:
        rr = base_reaction(rid, ev, "literature_example", ["other"], loc, None, None)
        rr["classification"] = classification_unresolved(ev, "Literature-context sketch was not re-extracted as a current-paper experiment or fully atom-mapped.")
        add_note(rr, note, note, ev)
        rr["review"]["unresolved"] = ["Literature example retained for provenance only; detailed conditions and structures were not supplied in this paper's experimental record."]
        reactions.append(rr)

    # Optimization Table 1. Entry 10 is also the 3aa scope result on page 3.
    table = [
        (1, None, "2:1", "toluene", 72, "no_reaction", None),
        (2, "PPhMe2", "2:1", "toluene", 4, 65, None),
        (3, "P1", "2:1", "CH2Cl2", 72, "trace", None),
        (4, "P2", "2:1", "CH2Cl2", 72, 40, 28),
        (5, "P3", "2:1", "CH2Cl2", 72, "trace", None),
        (6, "P4", "2:1", "CH2Cl2", 72, "trace", None),
        (7, "P5", "2:1", "CH2Cl2", 72, "trace", None),
        (8, "P6", "2:1", "CH2Cl2", 24, 43, 95),
        (9, "P6", "2:1", "toluene", 24, 46, 97),
        (10, "P6", "1:2", "toluene", 24, 75, 97),
        (11, "P6", "1:1.5", "toluene", 24, 72, 97),
        (12, "P6", "1:1.2", "toluene", 24, 69, 97),
        (13, "P7", "1:2", "toluene", 24, 49, 94),
        (14, "P8", "1:2", "toluene", 24, 47, 94),
        (15, "P9", "1:2", "toluene", 24, 58, 96),
    ]
    for entry, catalyst, ratio, solvent, hours, yld, ee in table:
        ev = [E[2], E[3]] if entry == 10 else [E[2]]
        cats = ["optimization", "substrate_scope"] if entry == 10 else ["optimization"]
        rr = base_reaction(f"table1_entry_{entry:02d}", ev, "specific_experiment", cats,
                           f"Table 1 entry {entry}" + ("; repeated as Figure 1 product 3aa" if entry == 10 else ""),
                           "asymmetric tandem isomerization/annulation", "tandem isomerization/[3+2] annulation")
        rr["general_scheme_ids"] = [general_id]
        p1 = participant("allenoate_1a", f"{DOC}:mol_1a", "reactant", ev)
        p2 = participant("allylamine_2a", f"{DOC}:mol_2a", "reactant", ev)
        if ratio == "2:1":
            p1["charges"]["amount"] = f("1a (0.20 mmol)", q(0.20, "mmol"), [E[2]])
            p2["charges"]["amount"] = f("2a (0.10 mmol)", q(0.10, "mmol"), [E[2]])
            p1["charges"]["equivalents"] = f("1a:2a ratio 2:1", q(2, "equiv"), [E[2]])
            p1["equivalence_reference"] = f(f"{DOC}:mol_2a", f"{DOC}:mol_2a", [E[2]])
            rr["reaction_scale"] = f("2a (0.10 mmol)", {"reference_molecule_id": f"{DOC}:mol_2a", "amount": q(0.10, "mmol")}, [E[2]])
            cat_ref = f"{DOC}:mol_2a"
        else:
            mult = float(ratio.split(":")[1])
            p2["charges"]["equivalents"] = f(f"1a:2a ratio {ratio}", q(mult, "equiv"), [E[2]])
            p2["equivalence_reference"] = f(f"{DOC}:mol_1a", f"{DOC}:mol_1a", [E[2]])
            cat_ref = f"{DOC}:mol_1a"
            if entry == 10:
                p1["charges"]["amount"] = f("1 (0.10 mmol)", q(0.10, "mmol"), [E[3]])
                p2["charges"]["amount"] = f("2 (0.20 mmol)", q(0.20, "mmol"), [E[3]])
                rr["reaction_scale"] = f("1 (0.10 mmol)", {"reference_molecule_id": f"{DOC}:mol_1a", "amount": q(0.10, "mmol")}, [E[3]])
        parts = [p1, p2]
        if yld != "no_reaction":
            parts.append(participant("product_3aa", f"{DOC}:mol_3aa", "product", ev, product_role="main product"))
        if catalyst:
            cp = participant("catalyst", f"{DOC}:mol_{catalyst}", "catalyst", [E[2]],
                             amount=("catalyst (0.02 mmol)", q(0.02, "mmol")),
                             molpct=("P* (20 mol%)", q(20, "mol%")), reference=cat_ref)
            parts.append(cp)
        solvent_mid = f"{DOC}:mol_toluene" if solvent == "toluene" else f"{DOC}:mol_dcm"
        parts.append(participant("solvent", solvent_mid, "solvent", [E[2]], volume=("1 mL of solvent", q(1, "mL"))))
        rr["participants"] = parts
        rr["steps"] = [step("reaction", ev, [p["participant_id"] for p in parts], operation="carried out", temperature=("25 °C", q(25, "°C")), time=(f"{hours} h", q(hours, "h")))]
        set_solvent(rr, solvent, ("1 mL of solvent", q(1, "mL")), [E[2]])
        finish_procedure(rr, "Unless indicated otherwise, reactions of 1a (0.20 mmol) and 2a (0.10 mmol) were carried out at 25 °C in the presence of a catalyst (0.02 mmol) in 1 mL of solvent.", [E[2]], "Table 1 footnote a")
        add_note(rr, f"1a:2a ratio {ratio}", {"1a_to_2a_ratio": ratio}, [E[2]])
        if yld == "no_reaction":
            rr["outcomes"].append(outcome("no_reaction", "no_reaction", [], [E[2]], raw="0% yield; text states no reaction occurred", normalized=True, scope={"observation": "no reaction"}))
            rr["classification"] = classification_na([E[2]], "No reaction was observed; no intended bond formation is classified.")
        elif yld == "trace":
            rr["outcomes"].append(outcome("trace_3aa", "trace", [f"{DOC}:mol_3aa"], [E[2]], raw="trace", normalized="trace", product_role="main product", basis="isolated yield", scope={"metric": "product yield"}))
            rr["classification"] = annulation_classification(ev)
        else:
            rr["outcomes"].append(outcome("yield_3aa", "yield", [f"{DOC}:mol_3aa"], ev, raw=f"{yld}%", normalized=q(yld, "%"), product_role="main product", yield_type="isolated", basis="isolated yield", scope={"metric": "yield"}))
            rr["classification"] = annulation_classification(ev)
        if yld not in {"no_reaction"}:
            rr["outcomes"].append(outcome("dr_3aa", "selectivity", [f"{DOC}:mol_3aa"], [E[2]], raw="dr is >20:1", normalized=">20:1", product_role="main product", basis="1H NMR analysis", scope={"metric": "dr", "order": "major:minor; configurations not named"}))
        if ee is not None:
            ee_ev = ev if entry == 10 else [E[2]]
            rr["outcomes"].append(outcome("ee_3aa", "selectivity", [f"{DOC}:mol_3aa"], ee_ev, raw=f"{ee}% ee", normalized=q(ee, "%"), product_role="main product", basis="chiral HPLC analysis", scope={"metric": "ee", "major_enantiomer": "source-drawn product; R/S label not stated in table"}))
        if entry in {11, 12, 13, 14, 15}:
            rr["review"]["unresolved"].append("The table reports a changed 1a:2a ratio but does not explicitly assign the individual mmol charges; the default footnote cannot safely resolve the changed substrate amounts.")
        reactions.append(rr)

    # Figure 1 scope excluding 3aa, which is the same experiment as Table 1 entry 10.
    scope_rows = [
        ("ab", 24, 73, 98, None, None), ("ac", 24, 83, 99, None, None),
        ("ad", 48, 31, 95, None, None), ("ae", 60, 21, 96, None, None),
        ("af", 24, 67, [">97", ">97"], "1:1", None),
        ("ag", 24, 63, [">99", ">99"], "1:1", None),
        ("ah", 24, 74, [95, 96], "1:1", None),
        ("ai", 36, 56, 99, "1.5:1", "unseparated"),
        ("aj", 72, "trace", None, None, None), ("ak", 72, "NR", None, None, None),
        ("al", 48, 35, 93, None, None), ("am", 24, 70, 96, None, None),
        ("an", 24, 65, 98, None, None), ("ao", 24, 81, 97, None, None),
        ("ap", 24, 71, 95, None, None), ("aq", 28, 79, 97, None, None),
        ("ar", 24, 87, 98, None, None), ("as", 24, 83, ">99", None, None),
        ("at", 24, 72, [96, 99], "1:1", None),
        ("au", 60, 28, None, None, ">99% de"),
        ("av", 32, 65, 79, None, None), ("aw", 72, 21, 99, None, None),
        ("ax", 48, 32, 97, None, None),
        ("ba", 48, 47, 96, None, None), ("ca", 32, 57, 98, None, None),
    ]
    for suffix, hours, yld, ee, dr, special in scope_rows:
        product_label = f"3{suffix}"
        if suffix[0] == "a":
            allene_label, amine_label = f"1{suffix[1]}", "2a"
        else:
            allene_label, amine_label = "1a", f"2{suffix[0]}"
        rr = base_reaction(f"figure1_{product_label}", [E[3]], "specific_experiment", ["substrate_scope"],
                           f"Figure 1 product {product_label}", "asymmetric tandem isomerization/annulation", "tandem isomerization/[3+2] annulation")
        rr["general_scheme_ids"] = [general_id]
        p_allene = participant("allenoate", f"{DOC}:mol_{allene_label}", "reactant", [E[3]], amount=("1 (0.10 mmol)", q(0.10, "mmol")))
        p_amine = participant("allylamine", f"{DOC}:mol_{amine_label}", "reactant", [E[3]], amount=("2 (0.20 mmol)", q(0.20, "mmol")), equiv=("1:2 substrate ratio", q(2, "equiv")), reference=f"{DOC}:mol_{allene_label}")
        parts = [p_allene, p_amine]
        if yld != "NR":
            parts.append(participant("product", f"{DOC}:mol_{product_label}", "product", [E[3]], product_role="main product"))
        parts += [
            participant("P6", f"{DOC}:mol_P6", "catalyst", [E[3]], amount=("P6 (0.02 mmol)", q(0.02, "mmol")), molpct=("P6 (20 mol%)", q(20, "mol%")), reference=f"{DOC}:mol_{allene_label}"),
            participant("toluene", f"{DOC}:mol_toluene", "solvent", [E[3]], volume=("1 mL of toluene", q(1, "mL"))),
        ]
        rr["participants"] = parts
        rr["reaction_scale"] = f("1 (0.10 mmol)", {"reference_molecule_id": f"{DOC}:mol_{allene_label}", "amount": q(0.10, "mmol")}, [E[3]])
        rr["steps"] = [step("reaction", [E[3]], [p["participant_id"] for p in parts], operation="carried out", temperature=("25 °C", q(25, "°C")), time=(f"{hours} h", q(hours, "h")))]
        set_solvent(rr, "toluene", ("1 mL of toluene", q(1, "mL")), [E[3]])
        finish_procedure(rr, "Unless indicated otherwise, reactions of 1 (0.10 mmol) and 2 (0.20 mmol) were carried out at 25 °C in the presence of P6 (0.02 mmol) in 1 mL of toluene.", [E[3]], "Figure 1 caption")
        if yld == "NR":
            rr["outcomes"].append(outcome("no_reaction", "no_reaction", [], [E[3]], raw="NR", normalized=True, scope={"observation": "no reaction"}))
            rr["classification"] = classification_na([E[3]], "No reaction was reported; the drawn target is not treated as an observed product.")
            rr["review"]["unresolved"].append("Figure 1 draws the intended 3ak target but reports NR; no product participant or bond-formation label is asserted.")
        elif yld == "trace":
            rr["outcomes"].append(outcome("trace_product", "trace", [f"{DOC}:mol_{product_label}"], [E[3]], raw="trace", normalized="trace", product_role="main product", basis="isolated yield", scope={"metric": "product yield"}))
            rr["classification"] = annulation_classification([E[3]])
        else:
            rr["outcomes"].append(outcome("yield_product", "yield", [f"{DOC}:mol_{product_label}"], [E[3]], raw=f"{yld}% yield", normalized=q(yld, "%"), product_role="main product", yield_type="isolated", basis="isolated yield", scope={"metric": "yield"}))
            rr["classification"] = annulation_classification([E[3]])
        if yld != "NR":
            dr_value = dr or ">20:1"
            rr["outcomes"].append(outcome("dr_product", "selectivity", [f"{DOC}:mol_{product_label}"], [E[3]], raw=f"dr = {dr_value}" if dr else "The dr is >20:1 unless indicated otherwise", normalized=dr_value, product_role="main product", basis="1H NMR analysis", scope={"metric": "dr", "order": "source major:minor; configurations not named"}))
        if isinstance(ee, list):
            for idx, val in enumerate(ee, 1):
                raw = f"{val}% ee" if isinstance(val, (int, float)) else f"{val}% ee"
                norm = q(float(str(val).replace(">", "")), "%", qualifier="gt") if str(val).startswith(">") else q(val, "%")
                rr["outcomes"].append(outcome(f"ee_diastereomer_{idx}", "selectivity", [f"{DOC}:mol_{product_label}"], [E[3]], raw=raw, normalized=norm, product_role="main product", basis="chiral HPLC analysis", scope={"metric": "ee", "diastereomer": idx, "configuration": "not assigned"}))
        elif ee is not None:
            raw = f"{ee}% ee"
            norm = q(float(str(ee).replace(">", "")), "%", qualifier="gt") if str(ee).startswith(">") else q(ee, "%")
            rr["outcomes"].append(outcome("ee_product", "selectivity", [f"{DOC}:mol_{product_label}"], [E[3]], raw=raw, normalized=norm, product_role="main product", basis="chiral HPLC analysis", scope={"metric": "ee", "major_enantiomer": "source-drawn product; explicit R/S label not shown"}))
        if special == "unseparated":
            add_note(rr, "99% ee, unseparated", {"ee_percent": 99, "sample": "unseparated diastereomeric mixture"}, [E[3]])
        elif special == ">99% de":
            rr["outcomes"].append(outcome("de_product", "selectivity", [f"{DOC}:mol_{product_label}"], [E[3]], raw=">99% de", normalized=q(99, "%", qualifier="gt"), product_role="main product", scope={"metric": "de", "diastereomer_order": "not specified"}))
        reactions.append(rr)

    # Scale-up reaction with the page-3 scheme/prose yield conflict retained.
    rr = base_reaction("scheme2a_scale_up_3aa", [E[3], T[3]], "specific_experiment", ["scale_up"], "Scheme 2a and following prose", "scale-up tandem isomerization/annulation", "tandem isomerization/[3+2] annulation")
    rr["general_scheme_ids"] = [general_id]
    rr["participants"] = [
        participant("1a", f"{DOC}:mol_1a", "reactant", [E[3]], amount=("1a (1 mmol)", q(1, "mmol"))),
        participant("2a", f"{DOC}:mol_2a", "reactant", [E[3]], amount=("2a (2 mmol)", q(2, "mmol")), equiv=("2a (2 mmol) relative to 1a (1 mmol)", q(2, "equiv")), reference=f"{DOC}:mol_1a"),
        participant("3aa", f"{DOC}:mol_3aa", "product", [E[3]], product_role="main product"),
        participant("P6", f"{DOC}:mol_P6", "catalyst", [E[3]], molpct=("P6 (20 mol%)", q(20, "mol%")), reference=f"{DOC}:mol_1a"),
        participant("toluene", f"{DOC}:mol_toluene", "solvent", [E[3]]),
    ]
    rr["participants"][0]["charges"]["mass"] = f("1a (0.15 g, 1 mmol)", q(0.15, "g"), [E[3]])
    rr["reaction_scale"] = f("1a (0.15 g, 1 mmol)", {"reference_molecule_id": f"{DOC}:mol_1a", "amount": q(1, "mmol"), "mass": q(0.15, "g")}, [E[3]])
    rr["steps"] = [step("reaction", [E[3]], [p["participant_id"] for p in rr["participants"]], operation="scale-up reaction", temperature=("25 °C", q(25, "°C")), time=("36 h", q(36, "h")))]
    set_solvent(rr, "toluene", None, [E[3]])
    alt_66 = f("66% yield", q(66, "%"), [E[3]])
    alt_65 = f("65% yield", q(65, "%"), [T[3]])
    rr["outcomes"] = [
        outcome("product_mass_3aa", "product_amount", [f"{DOC}:mol_3aa"], [E[3]], raw="3aa (0.31 g)", normalized=q(0.31, "g"), product_role="main product", basis="reported product mass", scope={"metric": "product amount"}),
        outcome("yield_3aa_conflict", "yield", [f"{DOC}:mol_3aa"], [E[3], T[3]], raw=["66% yield (Scheme 2a)", "65% yield (prose)"], normalized=None, product_role="main product", status="conflict", alternatives=[alt_66, alt_65], scope={"metric": "yield", "span": "1a + 2a to 3aa"}),
        outcome("ee_3aa", "selectivity", [f"{DOC}:mol_3aa"], [E[3]], raw="97% ee", normalized=q(97, "%"), product_role="main product", scope={"metric": "ee", "configuration": "source-drawn product"}),
    ]
    add_note(rr, "Scheme 2a reports 66% yield; following prose reports 65% yield.", "Unresolved source conflict retained without choosing a value.", [E[3], T[3]])
    rr["review"]["unresolved"].append("Source conflict: Scheme 2a gives 66% yield whereas the following prose gives 65%; measurement basis is not stated.")
    rr["classification"] = annulation_classification([E[3]])
    reactions.append(rr)

    # Scheme 2b derivatizations.
    deriv_specs = [
        ("scheme2b_3aa_to_4", "3aa", "4", ["LiAlH4", "THF"], 72, 98, "LiAlH4", "THF", ("0 °C", q(0, "°C")), ("1 h", q(1, "h")), reduction_classification),
        ("scheme2b_4_to_5", "4", "5", ["triphosgene", "Et3N", "CH2Cl2"], 88, 92, "triphosgene", "CH2Cl2", None, None, chlorination_classification),
        ("scheme2b_4_to_6", "4", "6", ["PCC", "CH2Cl2"], 68, 99, "PCC", "CH2Cl2", ("RT", "r.t."), ("4 h", q(4, "h")), oxidation_classification),
    ]
    reagent_map = {"LiAlH4": "lialh4", "THF": "thf", "triphosgene": "triphosgene", "Et3N": "et3n", "CH2Cl2": "dcm", "PCC": "pcc"}
    for rid, start, product, agents, yld, ee, main_reagent, solvent_name, temp, tm, class_fn in deriv_specs:
        rr = base_reaction(rid, [E[3], T[3]], "specific_experiment", ["derivatization"], f"Scheme 2b, {start} to {product}", None,
                           {"4": "reduction of two ester groups to a diol", "5": "conversion of the diol to a dichloride", "6": "oxidation of the diol to a dialdehyde"}[product])
        parts = [participant("starting_material", f"{DOC}:mol_{start}", "reactant", [E[3]]), participant("product", f"{DOC}:mol_{product}", "product", [E[3]], product_role="main product")]
        for agent in agents:
            role = "solvent" if agent in {"THF", "CH2Cl2"} else "reagent"
            parts.append(participant(agent.replace("2", "_"), f"{DOC}:mol_{reagent_map[agent]}", role, [E[3]]))
        rr["participants"] = parts
        if product == "5":
            rr["steps"] = [
                step("addition_at_0C", [E[3]], [p["participant_id"] for p in parts], operation="treat diol 4 with triphosgene and Et3N in CH2Cl2", temperature=("0 °C", q(0, "°C"))),
                step("warm_to_room_temperature", [E[3]], [p["participant_id"] for p in parts], operation="warm from 0 °C to room temperature; ramp/hold boundary not separately reported", temperature=("RT", "r.t."), time=("1 h", q(1, "h"))),
            ]
            rr["review"]["unresolved"].append("The source gives '0 °C–RT, 1 h' but does not separate addition, warming-ramp, and hold durations; 1 h is retained with the warming/RT stage.")
        else:
            rr["steps"] = [step("reaction", [E[3]], [p["participant_id"] for p in parts], operation=f"treat {start} with {main_reagent}", temperature=temp, time=tm)]
        set_solvent(rr, solvent_name, None, [E[3]])
        rr["outcomes"] = [
            outcome(f"yield_{product}", "yield", [f"{DOC}:mol_{product}"], [E[3], T[3]], raw=f"{yld}% yield", normalized=q(yld, "%"), product_role="main product", yield_type="reported yield", scope={"metric": "yield", "span": f"{start} to {product}"}),
            outcome(f"ee_{product}", "selectivity", [f"{DOC}:mol_{product}"], [E[3], T[3]], raw=f"{ee}% ee", normalized=q(ee, "%"), product_role="main product", scope={"metric": "ee", "configuration": "source-drawn product"}),
        ]
        finish_procedure(rr, f"Scheme 2b: {start} treated with {main_reagent} to give {product}.", [E[3]], "Scheme 2b")
        rr["classification"] = class_fn([E[3]])
        rr["review"]["unresolved"].append("Yield measurement basis and ee measurement method are not reported in the supplied main article.")
        reactions.append(rr)

    # Scheme 3 controls and isotope-labeling experiment.
    for suffix, substrate in [("a", "2d"), ("b", "2e")]:
        rr = base_reaction(f"scheme3{suffix}_no_reaction", [E[3], T[3]], "specific_experiment", ["control", "mechanistic_study"], f"Scheme 3{suffix}", None, "protecting-group/carbonyl control")
        rr["general_scheme_ids"] = [general_id]
        rr["participants"] = [
            participant("control_substrate", f"{DOC}:mol_{substrate}", "reactant", [E[3]]),
            participant("1a", f"{DOC}:mol_1a", "reactant", [E[3]]),
            participant("P6", f"{DOC}:mol_P6", "catalyst", [E[3]], molpct=("P6 (20 mol%)", q(20, "mol%")), reference=f"{DOC}:mol_1a"),
            participant("toluene", f"{DOC}:mol_toluene", "solvent", [E[3]]),
        ]
        rr["steps"] = [step("reaction", [E[3]], [p["participant_id"] for p in rr["participants"]], operation="control reaction", temperature=("25 °C", q(25, "°C")), time=("72 h", q(72, "h")))]
        set_solvent(rr, "toluene", None, [E[3]])
        rr["outcomes"] = [outcome("no_reaction", "no_reaction", [], [E[3], T[3]], raw="No Reaction", normalized=True, scope={"observation": "no reaction"})]
        rr["classification"] = classification_na([E[3]], "Negative control; no observed transformation and no hypothetical product is classified.")
        rr["review"]["unresolved"].append("Substrate amounts, solvent volume, observation method, and detection limit are not reported in the supplied main article.")
        reactions.append(rr)

    rr = base_reaction("scheme3c_deuteration", [E[3], E[4], T[4]], "specific_experiment", ["mechanistic_study", "control"], "Scheme 3c and page-4 discussion", "deuteration experiment", "tandem isomerization/[3+2] annulation with partial deuterium incorporation")
    rr["general_scheme_ids"] = [general_id]
    rr["participants"] = [
        participant("2a", f"{DOC}:mol_2a", "reactant", [E[3]], amount=("2a (0.2 mmol)", q(0.2, "mmol")), equiv=("2a 0.2 mmol relative to 1a 0.1 mmol", q(2, "equiv")), reference=f"{DOC}:mol_1a"),
        participant("1a", f"{DOC}:mol_1a", "reactant", [E[3]], amount=("1a (0.1 mmol)", q(0.1, "mmol"))),
        participant("D2O", f"{DOC}:mol_d2o", "additive", [E[3], E[4]], equiv=("D2O (3 equiv)", q(3, "equiv")), reference=f"{DOC}:mol_1a"),
        participant("P6", f"{DOC}:mol_P6", "catalyst", [E[3]], molpct=("P6 (20 mol%)", q(20, "mol%")), reference=f"{DOC}:mol_1a"),
        participant("toluene", f"{DOC}:mol_toluene", "solvent", [E[3]]),
        participant("deuterated_3aa", f"{DOC}:mol_3aa_D_mixture", "product", [E[3], E[4]], product_role="partially deuterated product mixture"),
    ]
    rr["reaction_scale"] = f("1a (0.1 mmol)", {"reference_molecule_id": f"{DOC}:mol_1a", "amount": q(0.1, "mmol")}, [E[3]])
    rr["steps"] = [step("reaction", [E[3]], [p["participant_id"] for p in rr["participants"]], operation="deuteration experiment under optimized reaction conditions", temperature=("25 °C", q(25, "°C")), time=("24 h", q(24, "h")))]
    set_solvent(rr, "toluene", None, [E[3]])
    rr["outcomes"] = [
        outcome("yield_deuterated_3aa", "yield", [f"{DOC}:mol_3aa_D_mixture"], [E[3]], raw="73% yield", normalized=q(73, "%"), product_role="partially deuterated product mixture", yield_type="reported yield", scope={"metric": "yield"}),
        outcome("ee_deuterated_3aa", "selectivity", [f"{DOC}:mol_3aa_D_mixture"], [E[3]], raw="97% ee", normalized=q(97, "%"), product_role="partially deuterated product mixture", scope={"metric": "ee", "major_enantiomer": "not assigned by R/S"}),
        outcome("D1_incorporation", "other", [f"{DOC}:mol_3aa_D_mixture"], [E[3]], raw="D1 (38%)", normalized=q(38, "%"), product_role="partially deuterated product mixture", scope={"metric": "deuterium incorporation", "site": "D1"}),
        outcome("D2_incorporation", "other", [f"{DOC}:mol_3aa_D_mixture"], [E[3]], raw="D2 (7%)", normalized=q(7, "%"), product_role="partially deuterated product mixture", scope={"metric": "deuterium incorporation", "site": "D2"}),
        outcome("D3_incorporation", "other", [f"{DOC}:mol_3aa_D_mixture"], [E[3]], raw="D3 (20%)", normalized=q(20, "%"), product_role="partially deuterated product mixture", scope={"metric": "deuterium incorporation", "site": "D3"}),
    ]
    add_note(rr, "successful incorporation of deuterium into the product indicated that an isomerization of substrate 2a occurred", {"author_interpretation": "supports allylamine isomerization"}, [T[4]])
    rr["classification"] = annulation_classification([E[3], E[4]])
    rr["review"]["unresolved"].append("Yield/ee/deuterium-incorporation measurement methods and the isotopologue distribution are not reported in the supplied main article.")
    reactions.append(rr)

    # Scheme 4 is an author-proposed catalytic cycle, not a set of isolated experiments.
    rr = base_reaction("scheme4_proposed_mechanism", [E[4], T[4]], "proposed_mechanism", ["other"], "Scheme 4", "plausible mechanism", "proposed tandem isomerization/annulation catalytic cycle")
    rr["participants"] = [
        participant("1a", f"{DOC}:mol_1a", "reactant", [E[4]]),
        participant("2a", f"{DOC}:mol_2a", "reactant", [E[4]]),
        participant("P6", f"{DOC}:mol_P6", "catalyst", [E[4]]),
        participant("3aa", f"{DOC}:mol_3aa", "product", [E[4]], product_role="mechanistic cycle product"),
    ]
    rr["classification"] = classification_unresolved([E[4]], "The arrows and Int-1–Int-4/2a-1–2a-3 structures are an author-proposed mechanism, not separately verified synthetic experiments.")
    add_note(rr, "1,4-H shift; Re-face attack; annulation step; 1,2-proton shift and catalyst elimination", {"sequence": ["allenoate activation", "allylamine 1,4-H shift to imine", "Re-face C–C attack", "intramolecular C–N closure", "proton shift and catalyst elimination"]}, [E[4], T[4]])
    rr["review"]["unresolved"] = ["Intermediates Int-1–Int-4 and 2a-1–2a-3 are proposed, not isolated; DFT details are in unavailable SI and no exact atom mapping is asserted."]
    reactions.append(rr)

    bundle["reactions"] = reactions
    bundle["source_context"] = [
        {"context_id": f"{DOC}:context_research_gap", "document_id": DOC, "evidence_ids": [E[2]], "kind": "discussion", "text": "Scheme 1a right is an unconquered design target, not an experiment."},
        {"context_id": f"{DOC}:context_si_notice", "document_id": DOC, "evidence_ids": [E[4]], "kind": "availability_notice", "text": "The article advertises experimental procedures, characterization, HPLC, NMR and X-ray data in SI, but no SI file was supplied or associated."},
        {"context_id": f"{DOC}:context_xray", "document_id": DOC, "evidence_ids": [E[3], E[5]], "kind": "crystallography", "text": "Product 3ar is shown with an X-ray structure; CCDC 1988551 is listed. This is not a separate synthesis experiment."},
        {"context_id": f"{DOC}:context_references", "document_id": DOC, "evidence_ids": [E[5], E[6]], "kind": "bibliography", "text": "Pages 5–6 contain author information and references; cited reactions were not counted as current-paper experiments beyond the three explicit Scheme 1 literature sketches."},
    ]

    st.write(RUN / "submitted.json", bundle)

    # Raw transcription artifact: values read directly from Table 1/Figure 1 before normalization.
    raw = {
        "document_id": DOC,
        "table1": [
            {"entry": e, "catalyst": c or "none", "ratio_1a_2a": ratio, "solvent": solv, "time_h": h, "yield": y, "ee_percent": ee}
            for e, c, ratio, solv, h, y, ee in table
        ],
        "figure1": [
            {"product": f"3{s}", "time_h": h, "yield": y, "ee": ee, "dr": dr, "special": special}
            for s, h, y, ee, dr, special in scope_rows
        ],
        "scheme2": {
            "scale_up": {"1a": "0.15 g, 1 mmol", "2a": "2 mmol", "P6": "20 mol%", "temperature": "25 °C", "time": "36 h", "product_mass": "0.31 g", "scheme_yield": "66%", "prose_yield": "65%", "ee": "97%"},
            "derivatizations": ["3aa -> 4: LiAlH4, THF, 0 °C, 1 h, 72%, 98% ee", "4 -> 5: triphosgene, Et3N, CH2Cl2, 0 °C-RT, 1 h, 88%, 92% ee", "4 -> 6: PCC, CH2Cl2, RT, 4 h, 68%, 99% ee"],
        },
        "scheme3": ["2d + 1a: no reaction", "2e + 1a: no reaction", "2a + 1a + 3 equiv D2O: 73%, 97% ee; D1 38%, D2 7%, D3 20%"],
        "limitations": "Literal transcription only; crop receipts retain image lineage. Not a source-reviewed dataset.",
    }
    st.write(RUN / "raw_model/visual_transcription.json", raw)

    stereo_checks = {
        m["molecule_id"]: {
            "input_smiles": st.value(m["smiles"]),
            "diagnostic": st.stereo_inspect(st.value(m["smiles"])),
            "source_verified": False,
        }
        for m in bundle["molecules"] if st.value(m["smiles"])
    }
    st.write(RUN / "raw_model/stereo_checks.json", {
        "command_equivalent": "reaction_tools.py standard stereo --smiles <candidate>",
        "chemical_accuracy_verified": False,
        "results": stereo_checks,
    })
    st.write(RUN / "raw_model/source_mentions.json", st.source_mentions(
        data, ["1a", "2a", "3aa", "3ar", "product 4"], context=220
    ))

    # Required page coverage: exactly pages 1–6, all attempted, no accuracy claim.
    ids_page2 = [general_id] + [f"{DOC}:{x[0]}" for x in literature_rows] + [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 16)]
    ids_scope = [f"{DOC}:table1_entry_10"] + [f"{DOC}:figure1_3{s}" for s, *_ in scope_rows]
    ids_page3 = ids_scope + [
        f"{DOC}:scheme2a_scale_up_3aa", f"{DOC}:scheme2b_3aa_to_4", f"{DOC}:scheme2b_4_to_5", f"{DOC}:scheme2b_4_to_6",
        f"{DOC}:scheme3a_no_reaction", f"{DOC}:scheme3b_no_reaction", f"{DOC}:scheme3c_deuteration",
    ]
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title, graphical abstract, and abstract", "finding": "True article title and a generic current tandem annulation scheme; graphical abstract is not a separate concrete experiment.", "evidence_ids": [E[1], T[1]]},
                {"region": "two-column introduction", "finding": "Background and cited prior work; no current-paper concrete reaction record on this page.", "evidence_ids": [E[1], T[1]]},
            ], "record_ids": [general_id], "unresolved": ["Graphical abstract uses generic R and PG variables; exact template scope is not exhaustively enumerated."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Scheme 1a–c", "finding": "Three prior-work/context panels and the current generic tandem reaction; the unconquered target is discussion, not an experiment.", "evidence_ids": [E[2], T[2]]},
                {"region": "Table 1 entries 1–15 and footnotes", "finding": "Fifteen independent optimization experiments, including entry 1 no reaction and trace outcomes in entries 3, 5, 6, and 7.", "evidence_ids": [E[2], T[2]]},
                {"region": "bottom-right scope discussion", "finding": "Textual lead-in to Figure 1; detailed per-product records are tied to page 3 image evidence.", "evidence_ids": [E[2], T[2]]},
            ], "record_ids": ids_page2, "unresolved": ["Changed ratios in Table 1 entries 11–15 do not explicitly assign individual mmol charges; these were not guessed."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Figure 1 product grid and caption", "finding": "Twenty-six scope results. Product 3aa duplicates Table 1 entry 10 and is merged there; the other 25 are separate records. 3ak is no reaction and has no observed-product participant.", "evidence_ids": [E[3], T[3]]},
                {"region": "Scheme 2a", "finding": "One scale-up experiment; scheme 66% yield conflicts with prose 65%, both retained.", "evidence_ids": [E[3], T[3]]},
                {"region": "Scheme 2b", "finding": "Three separate derivatization arrows: 3aa→4, 4→5, and 4→6.", "evidence_ids": [E[3], T[3]]},
                {"region": "Scheme 3a–c", "finding": "Two no-reaction controls and one D2O incorporation experiment with three site-specific incorporation percentages.", "evidence_ids": [E[3], T[3]]},
            ], "record_ids": ids_page3, "unresolved": ["Dense scope structures are source-labelled; many complete SMILES and explicit R/S descriptors remain unresolved without SI.", "Scheme 2/3 yield and selectivity measurement methods are generally not stated."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "Scheme 4 catalytic cycle", "finding": "Author-proposed mechanism with Int-1–Int-4 and 2a-1–2a-3; one proposed_mechanism record, not isolated experiments.", "evidence_ids": [E[4], T[4]]},
                {"region": "control-experiment continuation", "finding": "Text interpretation of D2O incorporation as evidence for substrate 2a isomerization; linked to Scheme 3c record.", "evidence_ids": [E[4], T[4]]},
                {"region": "Supporting Information notice", "finding": "SI is advertised but no SI PDF was supplied or explicitly associated.", "evidence_ids": [E[4], T[4]]},
            ], "record_ids": [f"{DOC}:scheme3c_deuteration", f"{DOC}:scheme4_proposed_mechanism"], "unresolved": ["DFT, characterization, HPLC/NMR and X-ray details cited as SI are unavailable in this session."]},
            {"page": 5, "state": "attempted", "regions": [
                {"region": "Accession Codes and author information", "finding": "CCDC 1988551 and metadata only; no synthesis experiment.", "evidence_ids": [E[5], T[5]]},
                {"region": "References 1–6", "finding": "Bibliography only; cited reactions are not current-paper experiments.", "evidence_ids": [E[5], T[5]]},
            ], "record_ids": [], "unresolved": []},
            {"page": 6, "state": "attempted", "regions": [
                {"region": "References 7–12", "finding": "Bibliography only; no current-paper experiment, table, or reaction arrow.", "evidence_ids": [E[6], T[6]]},
            ], "record_ids": [], "unresolved": []},
        ],
    }
    st.write(RUN / "coverage/pages_0001_0006.json", coverage)

    # Mark prepared page tasks attempted without claiming semantic verification.
    tasks = json.loads((RUN / "tasks.json").read_text())
    for task in tasks:
        if task["document_id"] == DOC and 1 <= int(task["evidence_id"].split(":page")[1].split(":")[0]) <= 6:
            task["state"] = "attempted"
            task["reason"] = "Original page image inspected; see coverage/pages_0001_0006.json. Semantic accuracy remains unverified."
    (RUN / "tasks.json").write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

    workflow = f"""# Workflow log — original PDF pages 1–6 of 6

- Batch: `allene-stereo-v31-200-20261017`
- Document: `{DOC}`
- Title: *Phosphine-Catalyzed Asymmetric Tandem Isomerization/Annulation of Allyl Amines with Allenoates: Enantioselective Annulation of a Saturated C–N Bond*
- Journal: *Organic Letters* 2021, 23, 9173–9178
- DOI: `10.1021/acs.orglett.1c03483`
- Source SHA-256: `dc20f9448ea23b61d83d8177c3224a0bb3938d347375f1ea6ee9f058f5cf984b`
- Source mapping: `presentation_conversion=false`; these are original journal-PDF pages, not converted slides.
- Visually inspected: original PDF pages 1, 2, 3, 4, 5, and 6, exactly the supplied range.
- SI association: `unestablished`. Page 4 advertises SI, but no SI PDF was supplied or explicitly associated.

## Isolation from stale workspace content

The pre-existing `runs/workflow.md` described a different article (`doc_19e82cae555b0ca9`, DOI `10.1021/acs.orglett.1c03482`) and was treated as cross-paper contamination. No molecule, condition, outcome, reaction ID, structure, or classification from that article was reused. The prepared current dataset contained zero candidates and zero experiments before extraction. Shared tools, schemas, scientific rules, and tests were not edited.

## Visual inspection and raw artifacts

All six supplied page images were opened at original detail. Enlarged source crops with evidence receipts are under `runs/current/raw_model/crops/`; literal Table 1, Figure 1, Scheme 2, and Scheme 3 values are in `runs/current/raw_model/visual_transcription.json`. Page decisions are in `runs/current/coverage/pages_0001_0006.json`; every page is `attempted` and `semantic_accuracy_verified=false`.

## Record boundaries

1. The page-1 graphical abstract, page-2 Scheme 1c, and page-3 Figure 1 header are one current `generic_scheme`, not three experiments.
2. Scheme 1 contributes three `literature_example` records (well-developed aromatic-imine annulation, Lu 2012 aliphatic-imine annulation, and Terada 2008 relay isomerization/C–C formation). The right-hand “unconquered” target is discussion only.
3. Table 1 contributes 15 optimization experiments, one per row. Entry 1 is a genuine no-reaction control and has no fabricated product. Trace rows retain qualitative `trace`. All reported table yields use the explicit isolated-yield footnote; dr is linked to 1H NMR and ee to chiral HPLC.
4. Figure 1 has 26 displayed scope outcomes. The 3aa result (24 h, 75%, 97% ee) is the same experiment as Table 1 entry 10 and is merged into that record with both page evidences. The other 25 outcomes are separate records. The 3ak target is drawn but marked NR, so it is not an observed product participant.
5. Scheme 2a is one scale-up record. Its graphic reports 66% yield, while the following prose says 65%; both alternatives are retained as an unresolved source conflict. Scheme 2b contains three independent arrows: 3aa→4, 4→5, and 4→6.
6. Scheme 3 contains two no-reaction controls and one D2O mechanistic experiment. The D1/D2/D3 percentages are separate site-specific outcomes on an explicitly unresolved isotopologue mixture.
7. Scheme 4 is one `proposed_mechanism`; Int-1–Int-4 and 2a-1–2a-3 are not registered as isolated products or independently verified experiments.
8. Pages 5–6 contain CCDC metadata, author information, and bibliography. No citation, molecular metadata item, or spectrum was converted into a current-paper synthesis.

## Molecules and stereochemistry

Every registry item has explicit `chirality_status`, `geometry_status`, and `representation_status`. Candidate SMILES were restricted to structures that could be transcribed with reasonable connectivity confidence (1a, 2a, 2d, 2e, 3aa, 4–6, and simple reagents/solvents) and are retained as machine-derived, not source-certified. `standard stereo` diagnostics are saved separately. The potentially stereogenic allene axis in 1a remains unassigned; no racemate was inferred. The E geometry visibly drawn for 2a/2d/2e is encoded. Products 3 and derivatives retain source wedge/dash assignments and the author's analogy to the 3ar X-ray result, but arbitrary R/S or @/@@ labels were not invented; their candidate SMILES therefore have `representation_status=not_encoded`. Products with reported low dr or stereochemically complex substituents are marked partially specified/mixture. Dense scope structures and chiral phosphines remain source-labelled with complete SMILES unresolved. No CAS number was guessed.

## Conditions, outcomes, and classification

The submission contains 52 v3.1 reactions: 47 specific experiments (15 optimization rows including the shared 3aa scope record, 25 additional scope records, 1 scale-up, 3 derivatizations, and 3 controls/mechanistic experiments), one current generic scheme, three literature examples, and one proposed mechanism. All records remain `needs_review` and `partial`.

Table and scope reagent doses, solvent volumes, temperatures, times, ratios, product-specific yields, dr/de/ee, and measurement methods are structured where explicit. Table entries 11–15 preserve their ratios while leaving individual substrate mmol charges unresolved. Scheme 2/3 measurement methods are left missing where the main article does not report them. The triphosgene arrow is split into 0 °C addition and warming/room-temperature stages; the source's combined 1 h does not distinguish ramp and hold, which is retained as unresolved.

Annulation records receive partial class-level C–C, C–N, and net C–H formation, together with N–H loss and allene bond-order decrease. Exact atom mapping, proton provenance, and coproduct balance are unresolved. The reductions receive C–H/O–H labels, dichlorination C–Cl, and PCC oxidation DB-C-O. Negative controls receive no hypothetical bond-formation label. All calculated changes are inferred and no model votes were fabricated.

## Validation and unresolved limits

After generation, `standard validate`, `standard merge`, strict dataset validation, and audit are run. Audit warnings are reviewed against the source and preserved where the main article truly omits a measurement basis. Expected unresolved limitations are: no associated SI; many dense scope/catalyst structures lack complete candidate SMILES; no exact atom maps; table entries 11–15 lack explicit individual mmol charges; Scheme 2a has a 65/66% yield conflict; Scheme 2/3 often omit yield/ee/deuterium measurement methods; 0 °C→RT timing is not split into ramp versus hold; and absolute R/S labels for the source-drawn product series are unavailable in the supplied pages. Schema validity is not chemical-accuracy certification, and attempted page coverage is not a claim that every chemical fact has been independently verified.
"""
    (ROOT / "runs/workflow.md").write_text(workflow)
    print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]), "pages": 6}, ensure_ascii=False))


if __name__ == "__main__":
    main()
