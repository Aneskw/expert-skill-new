#!/usr/bin/env python3
"""Build the source-specific v3.1 bundle for Chem. Commun. 2021, 57, 9212–9215."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

from standard_records import (  # noqa: E402
    fact,
    new_molecule,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    source_mentions,
    stereo_inspect,
)

DOC = "doc_eb3e0b8d7abf9bf7"
P = {n: f"{DOC}:page{n}:text" for n in range(1, 5)}
I = {n: f"{DOC}:page{n}:image" for n in range(1, 5)}
TITLE = "Cu-catalyzed coupling of vinylidene cyclopropanes with allyl and allenyl boronates"
GENERAL = f"{DOC}:generic_coupling"


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def F(raw=None, normalized=None, ev=(), method="original_extraction", status=None,
      verification=None, derivation=None, alternatives=()):
    return fact(raw, normalized, ev, method, status, verification, derivation, alternatives)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}
    if qualifier:
        value["qualifier"] = qualifier
    return value


def molecule(mid, label, ev, *, name=None, smiles=None, chirality="unspecified",
             geometry="unspecified", representation=None, axial=None, assignment=None,
             mixture=None, note=None):
    m = new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name or label, name or label, ev)
    m["names"]["english"] = F(name or label, name or label, ev)
    m["names"]["abbreviations"] = F([], [], ev)
    m["names"]["aliases"] = F([label], [label], ev)
    m["cas_numbers"] = F(None, None, ev, status="unreported", verification="missing")
    m["generic_template"] = None
    if smiles:
        diag = stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = F(
            f"structure drawing labelled {label}", diag["canonical_isomeric_smiles"], ev,
            method="structure_derivation", verification="inferred",
            derivation="Conservative candidate graph transcribed from the source drawing; parsing does not verify identity.",
        )
        m["stereochemistry"]["isomeric_smiles"] = F(
            f"candidate graph for {label}", diag["canonical_isomeric_smiles"], ev,
            method="structure_derivation", verification="inferred",
            derivation="RDKit canonical isomeric SMILES for the source-derived candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = F(
            None, diag["stereo_graph_key"], ev, method="structure_derivation",
            verification="inferred", derivation="Stereo graph key calculated from the candidate SMILES.",
        )
        rep = representation or "partial"
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = F(None, None, ev, status="unreported", verification="missing")
        m["stereochemistry"]["isomeric_smiles"] = F(None, None, ev, status="unreported", verification="missing")
        m["stereochemistry"]["identity_key"] = F(None, None, ev, status="unreported", verification="missing")
        rep = representation or "not_encoded"
    m["stereochemistry"]["atom_centers"] = F(
        None if chirality in ("unspecified", "unresolved") else chirality,
        None if chirality in ("unspecified", "unresolved") else chirality, ev,
        status="unreported" if chirality in ("unspecified", "unresolved") else "pending_review",
        verification="missing" if chirality in ("unspecified", "unresolved") else "unreviewed",
    )
    m["stereochemistry"]["double_bonds"] = F(
        None if geometry in ("unspecified", "unresolved", "not_applicable") else geometry,
        None if geometry in ("unspecified", "unresolved", "not_applicable") else geometry, ev,
        status="not_applicable" if geometry == "not_applicable" else "unreported" if geometry in ("unspecified", "unresolved") else "pending_review",
        verification="missing" if geometry in ("unspecified", "unresolved", "not_applicable") else "unreviewed",
    )
    m["stereochemistry"]["axial"] = F(
        axial, axial, ev,
        method="structure_derivation" if axial and "not stereogenic" in str(axial) else "original_extraction",
        status="unreported" if axial is None else None,
        verification="missing" if axial is None else "inferred" if "not stereogenic" in str(axial) else "unreviewed",
        derivation="A terminal CH2 allene end has two identical substituents, so that depicted local allene axis is not stereogenic."
        if axial and "not stereogenic" in str(axial) else None,
    )
    m["stereochemistry"]["mixture"] = F(
        mixture, mixture, ev, status="unreported" if mixture is None else None,
        verification="missing" if mixture is None else "unreviewed",
    )
    m["stereochemistry"]["assignment"] = F(
        assignment, assignment, ev, status="unreported" if assignment is None else None,
        verification="missing" if assignment is None else "unreviewed",
    )
    for key, value in (("chirality_status", chirality), ("geometry_status", geometry), ("representation_status", rep)):
        m["stereochemistry"][key] = F(
            value, value, ev, method="inference", verification="inferred",
            derivation="Status records source assignments and explicit unresolved stereochemical elements; absence of a descriptor is not treated as achiral or racemic.",
        )
    m["identity_note"] = F(
        note or "Source-labelled molecule; CAS was not reported and no external identity lookup was performed.",
        note or "Source-labelled molecule; CAS was not reported and no external identity lookup was performed.", ev,
    )
    return m


def generic_molecule(mid, label, ev, fixed, variable="R"):
    m = molecule(mid, label, ev, name=label, chirality="unresolved", geometry="unresolved", representation="not_encoded")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": F(fixed, fixed, ev),
        "variables": [{
            "label": variable,
            "attachment_points": F("shown attachment point", "shown attachment point", ev),
            "definition": F(f"source variable {variable}", f"source variable {variable}", ev),
            "allowed_values": F(None, None, ev, status="unreported", verification="missing"),
            "constraints": F(None, None, ev, status="unreported", verification="missing"),
        }],
        "instance_mapping": [],
        "reason": F("R-group family in a generic scheme", "R-group family in a generic scheme", ev),
    }
    return m


def set_charge(participant, key, raw, normalized, ev, ref=None, inherited=False):
    participant["charges"][key] = F(
        raw, normalized, ev,
        method="inherited_general_condition" if inherited else "original_extraction",
    )
    if key in ("equivalents", "mol_percent"):
        participant["equivalence_reference"] = F(
            ref, ref, ev,
            method="inherited_general_condition" if inherited else "original_extraction",
        )


def participant(pid, mid, role, ev, *, amount=None, equivalents=None, mol_percent=None,
                volume=None, ref=None, product_role=None, inherited=False):
    p = new_participant(pid, mid, role, ev)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev)
    if amount:
        set_charge(p, "amount", amount[0], amount[1], ev, inherited=inherited)
    if equivalents:
        set_charge(p, "equivalents", equivalents[0], equivalents[1], ev, ref=ref, inherited=inherited)
    if mol_percent:
        set_charge(p, "mol_percent", mol_percent[0], mol_percent[1], ev, ref=ref, inherited=inherited)
    if volume:
        set_charge(p, "volume", volume[0], volume[1], ev, inherited=inherited)
    return p


def classification(r, ev, kind="cc_h"):
    c = r["classification"]
    c["mapped_reaction_smiles"] = F(None, None, ev, status="unreported", verification="missing")
    c["model_votes"] = []
    if kind == "na":
        c["bond_changes"] = []
        c["labels"] = F([], [], ev, method="inference", verification="inferred", derivation="No observed transformation; intended bond formation is not assigned.")
        c["coverage"] = F("not_applicable", "not_applicable", ev, method="inference", verification="inferred", derivation="No observed product transformation.")
        c["limitations"] = ["No bond-formation label is inferred from an intended or failed reaction."]
        return
    if kind == "unresolved":
        c["bond_changes"] = []
        c["labels"] = F([], [], ev, method="inference", verification="inferred", derivation="Generic or mechanistic drawing was not atom-mapped.")
        c["coverage"] = F("unresolved", "unresolved", ev, method="inference", verification="inferred", derivation="No safe complete graph mapping was asserted.")
        c["limitations"] = ["Generic variables and/or proposed intermediates prevent atom-level classification."]
        return
    changes = []
    if kind == "cc_h":
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
        ]
        labels = [
            {"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1"]},
            {"category_number": 1, "code": "C-H", "supporting_change_ids": ["b2"]},
        ]
    else:  # three-component trapping
        changes = [
            {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "b2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        ]
        labels = [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}]
    c["bond_changes"] = [F(x, x, ev, method="inference", verification="inferred", derivation="Class-level change inferred from the source reaction drawing and discussion.") for x in changes]
    c["labels"] = F(labels, labels, ev, method="inference", verification="inferred", derivation="Labels cover all asserted formed bonds; exact atom identities and full mapping are unresolved.")
    c["coverage"] = F("partial", "partial", ev, method="inference", verification="inferred", derivation="Class-level only; no complete atom mapping, coproduct mapping, or hydrogen-source assignment.")
    c["limitations"] = ["Exact atom identities and complete bond-order changes are unresolved.", "Hydrogen source and coproducts are not established where a net C-H gain is asserted."]


def base_reaction(suffix, ev, record_type, categories, location, name, transformation, scheme=True):
    r = new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = F(categories, categories, ev)
    r["location"] = F(location, location, ev)
    r["reaction_name"] = F(name, name, ev)
    r["transformation_type"] = F(transformation, transformation, ev)
    r["general_scheme_ids"] = [GENERAL] if scheme and record_type == "specific_experiment" else []
    r["review"] = {
        "status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
        "unresolved": ["Supporting information was advertised but not supplied or associated; complete procedures and characterization are unavailable."]
    }
    r["completeness"] = "partial"
    return r


def reaction_step(r, ev, temp="rt", operation="combine reactants under the reported conditions", time=None):
    s = new_step("reaction", ev)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    s["operation"] = F(operation, operation, ev)
    s["temperature"] = F(temp, temp, ev) if temp is not None else F(None, None, ev, status="unreported", verification="missing")
    s["time"] = F(time, time, ev) if time is not None else F(None, None, ev, status="unreported", verification="missing")
    r["steps"] = [s]


def set_solvent(r, name, volume_raw=None, volume=None, ev=()):
    r["solvent_conditions"]["composition"] = F(name, [name] if name else None, ev, status="unreported" if name is None else None, verification="missing" if name is None else None)
    r["solvent_conditions"]["total_volume"] = F(volume_raw, volume, ev, status="unreported" if volume is None else None, verification="missing" if volume is None else None)


def yield_outcome(oid, product, percent, ev, basis=None, yield_type=None):
    o = new_outcome(oid, "yield", [product], ev)
    o["product_role"] = F("reported product", "reported product", ev)
    o["value"] = F(f"{percent}%", q(percent, "%", f"{percent}%"), ev)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    return o


def simple_outcome(oid, kind, targets, raw, normalized, ev, basis=None, scope=None):
    o = new_outcome(oid, kind, targets, ev)
    o["product_role"] = F("reported outcome", "reported outcome", ev)
    o["value"] = F(raw, normalized, ev)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev)
    if scope:
        o["scope"] = F(scope, scope, ev)
    return o


def add_common_two_component(r, sub_mid, bor_mid, product_mid, ev, *, catalyst="cucl", base="naotbu", solvent="thf", inherited=False):
    r["participants"] = [
        participant("substrate", sub_mid, "reactant", ev, amount=("0.1 mmol", q(0.1, "mmol")), equivalents=("1 equiv.", q(1, "equiv")), ref=sub_mid, inherited=inherited),
        participant("boronate", bor_mid, "reactant", ev, amount=("0.11 mmol", q(0.11, "mmol")), equivalents=("1.1 equiv.", q(1.1, "equiv")), ref=sub_mid, inherited=inherited),
        participant("catalyst", catalyst, "catalyst", ev, mol_percent=("10 mol%", q(10, "mol%")), ref=sub_mid, inherited=inherited),
        participant("base", base, "reagent", ev, amount=("0.12 mmol", q(0.12, "mmol")), equivalents=("1.2 equiv.", q(1.2, "equiv")), ref=sub_mid, inherited=inherited),
        participant("solvent", solvent, "solvent", ev, volume=("0.5 mL", q(0.5, "mL")), inherited=inherited),
        participant("product", product_mid, "product", ev, product_role="desired product"),
    ]
    r["reaction_scale"] = F("0.1 mmol substrate", {"reference_molecule_id": sub_mid, "amount": q(0.1, "mmol")}, ev,
                              method="inherited_general_condition" if inherited else "original_extraction")
    set_solvent(r, "THF", "0.5 mL", q(0.5, "mL"), ev)
    reaction_step(r, ev)
    r["procedure"]["general_procedure_reference"] = F(
        "standard conditions" if inherited else "table footnote conditions",
        "standard conditions" if inherited else "table footnote conditions", ev,
        method="inherited_general_condition" if inherited else "original_extraction",
    )


def table1_record(entry, catalyst, base, solvent, result, *, merged_scope=False):
    ev = [I[2], P[2]]
    r = base_reaction(f"table1_entry{entry}" if not merged_scope else "table1_entry2_3a", ev, "specific_experiment",
                      ["optimization"] + (["substrate_scope"] if merged_scope else []),
                      f"Table 1 entry {entry}" + (" / Table 2 product 3a" if merged_scope else ""),
                      "Cu-catalyzed coupling of vinylidene cyclopropane 1a with allylBpin 2a",
                      "allylative ring opening of a vinylidene cyclopropane")
    parts = [
        participant("substrate", "1a", "reactant", ev, amount=("0.1 mmol", q(0.1, "mmol")), equivalents=("1 equiv.", q(1, "equiv")), ref="1a"),
        participant("boronate", "2a", "reactant", ev, amount=("0.11 mmol", q(0.11, "mmol")), equivalents=("1.1 equiv.", q(1.1, "equiv")), ref="1a"),
    ]
    if catalyst:
        parts.append(participant("catalyst", catalyst, "catalyst", ev, mol_percent=("10 mol%", q(10, "mol%")), ref="1a"))
    if base:
        parts.append(participant("base", base, "reagent", ev, amount=("0.12 mmol", q(0.12, "mmol")), equivalents=("1.2 equiv.", q(1.2, "equiv")), ref="1a"))
    parts.append(participant("solvent", solvent, "solvent", ev, volume=("0.5 mL", q(0.5, "mL"))))
    if result not in ("N.R.",):
        parts.append(participant("product", "3a", "product", ev, product_role="desired product"))
    r["participants"] = parts
    r["reaction_scale"] = F("1a (0.1 mmol)", {"reference_molecule_id": "1a", "amount": q(0.1, "mmol")}, ev)
    solvent_name = {"thf": "THF", "dioxane": "dioxane", "toluene": "toluene"}[solvent]
    set_solvent(r, solvent_name, "0.5 mL", q(0.5, "mL"), ev)
    reaction_step(r, ev)
    r["procedure"]["original_text"] = F("Table 1 footnote a conditions", "Table 1 footnote a conditions", ev)
    if result == "N.R.":
        r["outcomes"] = [simple_outcome("no_reaction", "no_reaction", [], "N.R.", "N.R.", ev)]
        classification(r, ev, "na")
    elif result == "<5":
        r["outcomes"] = [simple_outcome("yield", "yield", ["3a"], "<5%", q(5, "%", "<5%", "lt"), ev,
                                                basis="isolated product yield")]
        classification(r, ev, "cc_h")
    else:
        r["outcomes"] = [yield_outcome("yield", "3a", result, ev, basis="isolated product yield", yield_type="isolated")]
        classification(r, ev, "cc_h")
    if not catalyst:
        r["notes"].append(F("No CuCl; the reaction did not proceed.", "No CuCl; the reaction did not proceed.", ev))
    if not base:
        r["notes"].append(F("No base; product 3a was not detected.", "No base; product 3a was not detected.", ev))
    if merged_scope:
        r["notes"].append(F(None, "Table 1 entry 2 and Table 2 product 3a are the same 85% standard experiment and are merged.", ev,
                              method="inference", verification="inferred", derivation="Reactants, charges, CuCl/NaOtBu/THF conditions, temperature, product, and isolated yield match."))
    return r


def scope_record(suffix, location, sub, bor, product, percent, ev, *, selectivity=None, table=2):
    r = base_reaction(suffix, ev, "specific_experiment", ["substrate_scope"], location,
                      "Cu-catalyzed coupling of a vinylidene cyclopropane with an organoboron reagent",
                      "allylative ring opening" if table == 2 else "allenylative ring opening")
    add_common_two_component(r, sub, bor, product, ev)
    r["outcomes"] = [yield_outcome("yield", product, percent, ev, basis="isolated product yield", yield_type="isolated")]
    if selectivity:
        r["outcomes"].append(simple_outcome(
            "selectivity", "selectivity", [product], selectivity,
            {"metric": "E:Z ratio", "order": ["E", "Z"], "relation": ">", "ratio": [20, 1]}, ev,
            scope={"metric": "E:Z", "order": ["E", "Z"], "major": "E", "measurement_method": None},
        ))
    classification(r, ev, "cc_h")
    return r


def scheme2_records():
    ev = [I[3], P[2], P[3]]
    a = base_reaction("scheme2a_products7_8", ev, "specific_experiment", ["substrate_scope"], "Scheme 2a",
                      "Reaction of methyl-substituted allenyl boronate 6 with vinylidene cyclopropane 1a",
                      "allenyl/propargyl coupling and vinylidene-cyclopropane ring opening")
    add_common_two_component(a, "1a", "6", "7", ev, inherited=True)
    a["participants"].append(participant("product_8", "8", "product", ev, product_role="co-product"))
    a["steps"][0]["participant_ids"].append("product_8")
    a["outcomes"] = [
        simple_outcome("combined_yield", "yield", ["7", "8"], "85% (7+8)", q(85, "%", "85% (7+8)"), ev),
        simple_outcome("ratio", "selectivity", ["7", "8"], "7:8 = 1:1", {"metric": "product ratio", "order": ["7", "8"], "ratio": [1, 1]}, ev,
                       scope={"metric": "product ratio", "order": ["7", "8"], "measurement_method": None}),
    ]
    a["review"]["unresolved"].append("The combined-yield and 7:8-ratio measurement methods are not stated on pages 1–4.")
    classification(a, ev, "cc_h")
    b = base_reaction("scheme2b_product12", ev, "specific_experiment", ["substrate_scope"], "Scheme 2b",
                      "Reaction of homoallenyl boronate 11 with vinylidene cyclopropane 1a",
                      "dienylative ring opening")
    add_common_two_component(b, "1a", "11", "12", ev, inherited=True)
    b["participants"].append(participant("unobserved_13", "13", "product", ev, product_role="unobserved alternative product"))
    b["steps"][0]["participant_ids"].append("unobserved_13")
    b["outcomes"] = [
        yield_outcome("yield", "12", 89, ev),
        simple_outcome("not_observed", "no_reaction", ["13"], "13 not observed", "not observed", ev),
    ]
    b["review"]["unresolved"].append("The 89% yield and non-observation measurement methods are not stated on pages 1–4.")
    classification(b, ev, "cc_h")
    return [a, b]


def scheme3_record(number, boronate, electrophile, percent, *, isolated=False, page4=False):
    ev = [I[3], P[3]] + ([I[4], P[4]] if page4 else [])
    rid = f"scheme3_product{number}"
    r = base_reaction(rid, ev, "specific_experiment", ["substrate_scope"], f"Scheme 3 product {number}",
                      "Three-component coupling of a boronate, vinylidene cyclopropane 1a, and an electrophile",
                      "organoboron addition, cyclopropane opening, and Cu-enolate trapping")
    r["participants"] = [
        participant("substrate", "1a", "reactant", ev),
        participant("boronate", boronate, "reactant", ev),
        participant("electrophile", electrophile, "reactant", ev),
        participant("catalyst", "cucl", "catalyst", ev),
        participant("base", "naotbu", "reagent", ev),
        participant("solvent", "thf", "solvent", ev),
        participant("product", str(number), "product", ev, product_role="reported product"),
    ]
    set_solvent(r, "THF", None, None, ev)
    reaction_step(r, ev)
    r["outcomes"] = [yield_outcome("yield", str(number), percent, ev, basis="isolated product yield" if isolated else None,
                                    yield_type="isolated" if isolated else None)]
    if not isolated:
        r["review"]["unresolved"].append(f"The {percent}% yield measurement basis is not stated on pages 1–4.")
    r["review"]["unresolved"].append("Electrophile amount, catalyst/base loading, solvent volume, and reaction time are not reported in Scheme 3 or adjacent prose.")
    classification(r, ev, "cc_cc")
    return r


def main():
    template = json.loads((RUN / "raw_model/standard_init.json").read_text())
    dataset = json.loads((RUN / "dataset.json").read_text())
    dataset["documents"][0]["title"] = TITLE
    dump(RUN / "dataset.json", dataset)
    template["documents"] = copy.deepcopy(dataset["documents"])

    molecules = []
    add = molecules.append
    # Generic reaction families and simple source-derived candidates.
    add(generic_molecule("generic_allene_A", "allene A", [I[1]], "R-substituted allene", "R"))
    add(generic_molecule("generic_product_C", "1,5-diene C", [I[1]], "generic 1,5-diene product", "R"))
    add(generic_molecule("generic_vcp_D", "vinylidene cyclopropane D", [I[1]], "vinylidene cyclopropane bearing two electron-withdrawing groups E", "R"))
    add(generic_molecule("generic_product_G", "coupled product G", [I[1]], "ring-opened coupling product", "R"))
    add(molecule("1a", "1a", [I[2], I[3]], name="dimethyl-ester vinylidene cyclopropane 1a",
                 smiles="COC(=O)C1(C(=O)OC)CC1=C=C", chirality="achiral", geometry="unspecified", representation="partial",
                 axial="terminal-CH2 allene axis is not stereogenic on the candidate graph",
                 note="Source-labelled 1a; candidate graph retains the vinylidene-cyclopropane/allene unit and two methyl esters."))
    add(molecule("2a", "2a", [I[2], I[3]], name="allylBpin 2a", smiles="C=CCB1OC(C)(C)C(C)(C)O1",
                 chirality="achiral", geometry="not_applicable", representation="encoded"))
    add(molecule("4", "4", [I[3]], name="allenyl boronate 4", smiles="C=C=CB1OC(C)(C)C(C)(C)O1",
                 chirality="achiral", geometry="not_applicable", representation="partial",
                 axial="terminal-CH2 allene axis is not stereogenic on the candidate graph"))
    add(molecule("6", "6", [I[3]], name="methyl-substituted allenyl boronate 6", smiles="C=C=C(C)B1OC(C)(C)C(C)(C)O1",
                 chirality="unspecified", geometry="not_applicable", representation="partial",
                 axial="allene axial configuration unassigned by the source"))
    add(molecule("11", "11", [I[3]], name="homoallenyl boronate 11", smiles="C=C=CCB1OC(C)(C)C(C)(C)O1",
                 chirality="achiral", geometry="not_applicable", representation="partial",
                 axial="terminal-CH2 allene axis is not stereogenic on the candidate graph"))

    simple = {
        "cucl": ("CuCl", "[Cu]Cl", "unresolved", "not_applicable"),
        "cucl2": ("CuCl2", "Cl[Cu]Cl", "unresolved", "not_applicable"),
        "liotbu": ("LiOtBu", "[Li+].[O-]C(C)(C)C", "achiral", "not_applicable"),
        "naotbu": ("NaOtBu", "[Na+].[O-]C(C)(C)C", "achiral", "not_applicable"),
        "kotbu": ("KOtBu", "[K+].[O-]C(C)(C)C", "achiral", "not_applicable"),
        "liome": ("LiOMe", "[Li+].[O-]C", "achiral", "not_applicable"),
        "naome": ("NaOMe", "[Na+].[O-]C", "achiral", "not_applicable"),
        "kome": ("KOMe", "[K+].[O-]C", "achiral", "not_applicable"),
        "thf": ("THF", "C1CCOC1", "achiral", "not_applicable"),
        "dioxane": ("dioxane", "O1CCOCC1", "achiral", "not_applicable"),
        "toluene": ("toluene", "Cc1ccccc1", "achiral", "not_applicable"),
        "mei": ("methyl iodide", "CI", "achiral", "not_applicable"),
        "benzoyl_fluoride": ("benzoyl fluoride", "O=C(F)c1ccccc1", "achiral", "not_applicable"),
        "ethyl_bromoacetate": ("ethyl bromoacetate", "CCOC(=O)CBr", "achiral", "not_applicable"),
    }
    for mid, (name, smi, chi, geo) in simple.items():
        ev = [I[2]] if mid not in ("mei", "benzoyl_fluoride", "ethyl_bromoacetate") else [I[3], I[4]]
        add(molecule(mid, name, ev, name=name, smiles=smi, chirality=chi, geometry=geo, representation="encoded" if chi == "achiral" else "partial"))

    unresolved_specs = {
        "1b": ("vinylidene cyclopropane 1b (i-Pr ester groups)", [I[2], I[3]], "unspecified", "unspecified", None),
        "vcp_3h": ("ethyl-substituted vinylidene cyclopropane for 3h", [I[2]], "unspecified", "unspecified", "allene axis unassigned"),
        "vcp_5c": ("ethyl-substituted vinylidene cyclopropane for 5c", [I[3]], "unspecified", "unspecified", "allene axis unassigned"),
        "vcp_5d": ("phenyl-substituted vinylidene cyclopropane for 5d", [I[3]], "unspecified", "unspecified", "allene axis unassigned"),
        "vcp_5e": ("n-pentyl-substituted vinylidene cyclopropane for 5e", [I[3]], "unspecified", "unspecified", "allene axis unassigned"),
        "vcp_5f": ("i-Pr-substituted vinylidene cyclopropane for 5f", [I[3]], "unspecified", "unspecified", "allene axis unassigned"),
        "2_for_3b": ("gamma-phenyl allylboronate for 3b", [I[2]], "unspecified", "specified", None),
        "2_for_3c": ("gamma-(2-fluorophenyl) allylboronate for 3c", [I[2]], "unspecified", "specified", None),
        "2_for_3d": ("gamma-(2-chlorophenyl) allylboronate for 3d", [I[2]], "unspecified", "specified", None),
        "2_for_3e": ("beta-methyl allylboronate for 3e", [I[2]], "unspecified", "unspecified", None),
        "2_for_3f": ("branched allylboronate for 3f", [I[2]], "unspecified", "unspecified", None),
        "2_for_3g": ("methyl-substituted branched allylboronate for 3g", [I[2]], "unspecified", "unspecified", None),
        "2b": ("methyl-substituted allylboronate 2b", [I[3]], "unspecified", "not_applicable", None),
    }
    for mid, (name, ev, chi, geo, ax) in unresolved_specs.items():
        add(molecule(mid, mid, ev, name=name, chirality=chi, geometry=geo, representation="not_encoded", axial=ax))

    for label in ["3a", "3e", "3f", "3g", "3h", "3i", "7", "12", "13", "16", "17", "18"]:
        ev = [I[2]] if label.startswith("3") else [I[3]]
        add(molecule(label, label, ev, name=f"product {label}", chirality="unspecified", geometry="unspecified", representation="not_encoded"))
    for label in ["3b", "3c", "3d"]:
        add(molecule(label, label, [I[2]], name=f"product {label} (E/Z mixture, E major)", chirality="unspecified", geometry="mixture", representation="not_encoded",
                     mixture="E/Z mixture; E:Z >20:1", assignment="E is the reported major alkene geometry"))
    for label in ["5a", "5b", "5c", "5d", "5e", "5f"]:
        add(molecule(label, label, [I[3]], name=f"allenyl product {label}", chirality="achiral", geometry="not_applicable", representation="not_encoded",
                     axial="terminal-CH2 allene axis is not stereogenic on the depicted local structure"))
    add(molecule("8", "8", [I[3]], name="allenic product 8", chirality="unspecified", geometry="not_applicable", representation="not_encoded",
                 axial="allene axial configuration is unassigned by the source"))
    for label in ["19", "20"]:
        add(molecule(label, label, [I[3], I[4]], name=f"three-component product {label}", chirality="unspecified", geometry="not_applicable", representation="not_encoded"))

    reactions = []
    prior = base_reaction("literature_conventional_allene_allylic_coupling", [I[1], P[1]], "literature_example", [], "Scheme 1a",
                          "Well-established coupling of an allene with allylic electrophiles",
                          "generic 1,2-functionalization of an allene", scheme=False)
    prior["participants"] = [participant("allene", "generic_allene_A", "reactant", [I[1]]), participant("product", "generic_product_C", "product", [I[1]], product_role="generic 1,5-diene")]
    classification(prior, [I[1]], "unresolved")
    reactions.append(prior)

    gen = base_reaction("generic_coupling", [I[1], I[2], I[3]], "generic_scheme", [], "Scheme 1b and Tables 2–3 headers",
                        "Cu-catalyzed coupling of vinylidene cyclopropanes with allyl/allenyl boronates",
                        "organoboron addition at the terminal allene carbon with cyclopropane ring opening", scheme=False)
    gen["participants"] = [participant("substrate", "generic_vcp_D", "reactant", [I[1]]), participant("product", "generic_product_G", "product", [I[1]], product_role="generic coupled product")]
    classification(gen, [I[1]], "unresolved")
    reactions.append(gen)

    mech = base_reaction("proposed_ring_opening_mechanism", [I[1], P[1], I[3], P[3]], "proposed_mechanism", [], "Scheme 1b and Scheme 3 generic header",
                         "Author-proposed allylcopper addition, beta-carbo elimination, and Cu-enolate interception",
                         "proposed catalytic pathway", scheme=False)
    mech["participants"] = [participant("substrate", "generic_vcp_D", "reactant", [I[1]]), participant("product", "generic_product_G", "product", [I[1]], product_role="proposed generic product")]
    mech["notes"] = [F("Proposed intermediates E, F, I, and II are mechanistic species, not isolated experiments.", "Proposed intermediates E, F, I, and II are mechanistic species, not isolated experiments.", [I[1], I[3]])]
    classification(mech, [I[1], I[3]], "unresolved")
    reactions.append(mech)

    rows = {
        1: ("cucl", "liotbu", "thf", 27), 2: ("cucl", "naotbu", "thf", 85),
        3: ("cucl", "kotbu", "thf", 54), 4: ("cucl", "liome", "thf", 71),
        5: ("cucl", "naome", "thf", 80), 6: ("cucl", "kome", "thf", 58),
        7: ("cucl2", "naotbu", "thf", 67), 8: ("cucl", "naotbu", "dioxane", 63),
        9: ("cucl", "naotbu", "toluene", "<5"), 10: ("cucl", None, "thf", "N.R."),
        11: (None, "naotbu", "thf", "N.R."),
    }
    for entry, vals in rows.items():
        reactions.append(table1_record(entry, *vals, merged_scope=(entry == 2)))

    table2 = [
        ("table2_3b", "Table 2 product 3b", "1a", "2_for_3b", "3b", 73, "E:Z >20:1"),
        ("table2_3c", "Table 2 product 3c", "1a", "2_for_3c", "3c", 85, "E:Z >20:1"),
        ("table2_3d", "Table 2 product 3d", "1a", "2_for_3d", "3d", 81, "E:Z >20:1"),
        ("table2_3e", "Table 2 product 3e", "1a", "2_for_3e", "3e", 84, None),
        ("table2_3f", "Table 2 product 3f", "1a", "2_for_3f", "3f", 72, None),
        ("table2_3g", "Table 2 product 3g", "1a", "2_for_3g", "3g", 83, None),
        ("table2_3h", "Table 2 product 3h", "vcp_3h", "2a", "3h", 83, None),
        ("table2_3i", "Table 2 product 3i", "1b", "2a", "3i", 86, None),
    ]
    for suffix, location, sub, bor, product, percent, selectivity in table2:
        reactions.append(scope_record(suffix, location, sub, bor, product, percent, [I[2], P[2]],
                                      selectivity=selectivity, table=2))

    table3 = [
        ("table3_5a", "Table 3 product 5a", "1a", "4", "5a", 76),
        ("table3_5b", "Table 3 product 5b", "1b", "4", "5b", 90),
        ("table3_5c", "Table 3 product 5c", "vcp_5c", "4", "5c", 80),
        ("table3_5d", "Table 3 product 5d", "vcp_5d", "4", "5d", 74),
        ("table3_5e", "Table 3 product 5e", "vcp_5e", "4", "5e", 86),
        ("table3_5f", "Table 3 product 5f", "vcp_5f", "4", "5f", 83),
    ]
    for row in table3:
        reactions.append(scope_record(*row, [I[3], P[3]], table=3))
    reactions.extend(scheme2_records())
    reactions.extend([
        scheme3_record(16, "2a", "mei", 84),
        scheme3_record(17, "2b", "mei", 71, isolated=True),
        scheme3_record(18, "4", "mei", 85, isolated=True),
        scheme3_record(19, "4", "benzoyl_fluoride", 67, page4=True),
        scheme3_record(20, "4", "ethyl_bromoacetate", 75, page4=True),
    ])

    template["molecules"] = molecules
    template["reactions"] = reactions
    template["source_context"] = [
        {"context_id": f"{DOC}:metadata", "document_id": DOC, "evidence_ids": [P[1], I[1]], "summary": "Title, author, DOI, abstract, and ESI availability notice. ESI is advertised but was not supplied or associated."},
        {"context_id": f"{DOC}:intro", "document_id": DOC, "evidence_ids": [P[1], I[1]], "summary": "Introduction and cited background on vinylidene cyclopropane/allene chemistry; citations were not imported as current experiments beyond the generic Scheme 1a precedent."},
        {"context_id": f"{DOC}:conclusion_references", "document_id": DOC, "evidence_ids": [P[4], I[4]], "summary": "Conclusion, acknowledgements, conflicts, and references; no new synthesis arrow occurs on page 4."},
        {"context_id": f"{DOC}:negative_characterization_inventory", "document_id": DOC, "evidence_ids": [I[1], I[2], I[3], I[4]], "summary": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block occurs on pages 1–4."},
    ]

    dump(RUN / "submitted.json", template)
    diagnostics = {}
    for m in molecules:
        smi = m["smiles"]["normalized_value"]
        if smi:
            diagnostics[m["molecule_id"]] = stereo_inspect(smi)
    dump(RUN / "raw_model/stereo_checks.json", diagnostics)
    dump(RUN / "raw_model/source_mentions.json", source_mentions(dataset, ["1a", "2a", "3a", "4", "Scheme 2", "Scheme 3"]))

    page1_ids = [f"{DOC}:literature_conventional_allene_allylic_coupling", GENERAL, f"{DOC}:proposed_ring_opening_mechanism"]
    page2_ids = [r["reaction_id"] for r in reactions if r["reaction_id"].startswith(f"{DOC}:table1_") or r["reaction_id"].startswith(f"{DOC}:table2_")]
    page3_ids = [r["reaction_id"] for r in reactions if any(x in r["reaction_id"] for x in ("table3_", "scheme2", "scheme3_"))]
    page4_ids = [f"{DOC}:scheme3_product19", f"{DOC}:scheme3_product20"]
    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 1, "state": "attempted", "regions": [
                {"region": "title/abstract and two-column introduction", "finding": "Current-paper objective, DOI, ESI notice, and cited background; no concrete current experimental row.", "evidence_ids": [I[1], P[1]]},
                {"region": "Scheme 1a", "finding": "Well-established generic allene/allylic-electrophile literature pathway retained as a literature example.", "evidence_ids": [I[1]]},
                {"region": "Scheme 1b", "finding": "Current generic coupling and author-proposed allylcopper addition/ring-opening/electrophile-trapping pathway.", "evidence_ids": [I[1], P[1]]},
            ], "record_ids": page1_ids, "unresolved": ["No SI was supplied; proposed intermediates are not isolated experiments."]},
            {"page": 2, "state": "attempted", "regions": [
                {"region": "Table 1", "finding": "All 11 optimization rows recorded; entries 10–11 are no-reaction controls and entry 9 retains <5% rather than zero.", "evidence_ids": [I[2], P[2]]},
                {"region": "Table 2", "finding": "All nine products inventoried; 3a is merged with identical Table 1 entry 2, and 3b–3i are separate scope records.", "evidence_ids": [I[2], P[2]]},
                {"region": "adjacent prose", "finding": "Regioselectivity, E-selectivity, substrate descriptions, and preview of Table 3/Scheme 2 checked against the drawings.", "evidence_ids": [I[2], P[2]]},
            ], "record_ids": page2_ids, "unresolved": ["Compact scope-substrate/product drawings were not forced into complete SMILES where transcription was not sufficiently secure.", "Reaction times are not reported for Tables 1–3."]},
            {"page": 3, "state": "attempted", "regions": [
                {"region": "Table 3", "finding": "All six allenyl-boronate scope outcomes 5a–5f recorded with isolated yields.", "evidence_ids": [I[3], P[3]]},
                {"region": "Scheme 2", "finding": "Two concrete experiments: 7/8 combined 85% at 1:1, and 12 at 89% with alternative 13 not observed; mechanistic 9/10/14/15 were not made separate syntheses.", "evidence_ids": [I[3], P[3]]},
                {"region": "Scheme 3", "finding": "Five three-component arrows to 16–20 recorded separately; page-4 prose continuation linked to 19/20.", "evidence_ids": [I[3], P[3]]},
            ], "record_ids": page3_ids, "unresolved": ["Scheme 2 yield/ratio methods and most Scheme 3 yield methods are unreported.", "Scheme 3 electrophile amounts, reaction times, and detailed loadings/volumes are absent."]},
            {"page": 4, "state": "attempted", "regions": [
                {"region": "top-left continuation", "finding": "Completes product 19 (67%) and product 20 (75%) discussion and summarizes sequential C–C bond formation.", "evidence_ids": [I[4], P[4]]},
                {"region": "summary and acknowledgements", "finding": "No additional concrete experiment; conclusion describes the established coupling and Cu-enolate interception.", "evidence_ids": [I[4], P[4]]},
                {"region": "references", "finding": "Bibliography only; cited papers were not imported as current-paper experiments.", "evidence_ids": [I[4], P[4]]},
            ], "record_ids": page4_ids, "unresolved": ["Yield measurement bases for 19 and 20 are not stated on pages 1–4."]},
        ],
    }
    dump(RUN / "coverage/pages_0001_0004.json", coverage)
    dump(RUN / "raw_model/page_extraction.json", {
        "document_id": DOC,
        "pages_attempted": [1, 2, 3, 4],
        "page_images": [I[n] for n in range(1, 5)],
        "crops": [
            "crops/page1_scheme1.jpg", "crops/page2_table1.jpg", "crops/page2_table2.jpg",
            "crops/page3_table3.jpg", "crops/page3_scheme2.jpg", "crops/page3_scheme3.jpg",
        ],
        "deduplication": ["Table 1 entry 2 and Table 2 product 3a merged as one identical 85% experiment."],
        "inventory": {"table1_rows": 11, "table2_products": 9, "table3_products": 6, "scheme2_experiments": 2, "scheme3_arrows": 5},
        "negative_inventory": "No standalone mass-spectrum peak table, molecular-weight table, or fragment-relative-abundance block on pages 1–4.",
    })
    print(json.dumps({"molecules": len(molecules), "reactions": len(reactions), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
