#!/usr/bin/env python3
"""Build the page 1–6 machine submission from visually inspected source regions.

This is a source-derived machine extraction.  It deliberately retains uncertain
structures and measurement bases instead of completing them by chemical convention.
"""

from __future__ import annotations

import copy
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_358a77039b302aac"
EV = {p: f"{DOC}:page{p}:image" for p in range(1, 7)}
EVT = {p: f"{DOC}:page{p}:text" for p in range(1, 7)}

TITLE = (
    "An approach to unsymmetrical 3,3′-diindolylmethanes through Pd-catalyzed "
    "cascade Heck cyclization of allenamides and o-ethynylanilines"
)


def fact(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
         status=None):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation, status=status)


def q(number, unit, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    return value


def mid(key):
    return f"{DOC}:mol:{key}"


bundle = st.read(RUN / "submitted.json")
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
molecules = {}


CONCRETE_SMILES = {
    "1a": "C=C=CN(c1ccccc1I)S(=O)(=O)c1ccc(C)cc1",
    "2a": "Cc1ccc(S(=O)(=O)Nc2ccccc2C#Cc2ccccc2)cc1",
    "3aa": "Cc1ccc(S(=O)(=O)n2cc(Cc3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "3al": "Cc1ccc(S(=O)(=O)n2cc(Cc3c([Si](C)(C)C)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "8": "c1ccc2[nH]cc(Cc3c[nH]c4ccccc34)c2c1",
    "9": "Cc1ccc(S(=O)(=O)n2cc(C(=O)c3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "10": "Cc1ccc(S(=O)(=O)n2cc(C(O)c3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "11": "Cc1ccc(S(=O)(=O)n2cc(C(N=[N+]=[N-])c3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "12": "Cc1ccc(S(=O)(=O)n2cc(C(C#N)c3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "13": "Cc1ccc(S(=O)(=O)n2cc(C(CC=C)c3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
    "14": "Cc1ccc(S(=O)(=O)n2cc(C(C#Cc3ccccc3)c3c(-c4ccccc4)n(S(=O)(=O)c4ccc(C)cc4)c4ccccc34)c3ccccc32)cc1",
}


def ensure_molecule(key, label, evidence, note=None, smiles=None, chirality=None,
                    geometry=None, representation=None, inferred_label=False):
    """Create or extend one paper-scoped molecule record."""
    molecule_id = mid(key)
    evidence = list(dict.fromkeys(evidence))
    if molecule_id in molecules:
        m = molecules[molecule_id]
        occ = m["occurrences"][0]
        occ["evidence_ids"] = list(dict.fromkeys(occ["evidence_ids"] + evidence))
        occ["compound_label"]["evidence_ids"] = list(dict.fromkeys(
            occ["compound_label"]["evidence_ids"] + evidence))
        return molecule_id
    m = st.new_molecule(molecule_id, DOC, evidence)
    m["identity_status"] = "source_identified"
    if inferred_label:
        label_fact = fact(
            {"inferred_from_product_code": label}, label, evidence, method="inference",
            derivation="The first/second letter in product code 3xy maps to allenamide 1x/o-alkynylaniline 2y; the standalone substrate drawing is not printed in this main-text table."
        )
    else:
        label_fact = fact(label, label, evidence)
    m["occurrences"] = [{"document_id": DOC, "compound_label": label_fact,
                         "evidence_ids": evidence}]
    m["names"]["original"] = copy.deepcopy(label_fact)
    m["names"]["english"] = fact()
    m["names"]["abbreviations"] = fact()
    m["names"]["aliases"] = fact()
    m["cas_numbers"] = fact()
    m["identity_note"] = fact(note, note, evidence) if note else fact()
    candidate = smiles or CONCRETE_SMILES.get(key)
    if candidate:
        diagnostic = st.stereo_inspect(candidate)
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            {"source": "drawn structure", "compound_label": label}, candidate, evidence,
            method="structure_derivation",
            derivation="Candidate connectivity was transcribed from the displayed structure and checked for RDKit parsing; parsing does not certify source identity."
        )
        m["generic_template"] = None
        if chirality is None:
            chirality = "unspecified" if diagnostic["atom_centers"] else "achiral"
        if geometry is None:
            geometry = "not_applicable"
        if representation is None:
            representation = "not_encoded" if chirality == "unspecified" else "not_applicable"
        m["stereochemistry"]["atom_centers"] = fact(
            {"tool_check": diagnostic["atom_centers"]}, diagnostic["atom_centers"], evidence,
            method="structure_derivation",
            derivation="standard stereo graph diagnostic on the candidate SMILES; the source drawing carries no R/S wedge assignment."
        )
        m["stereochemistry"]["double_bonds"] = fact(
            {"tool_check": diagnostic["double_bonds"]}, diagnostic["double_bonds"], evidence,
            method="structure_derivation",
            derivation="standard stereo graph diagnostic; aromatic Kekule bonds are not interpreted as source E/Z assignments."
        )
        if key == "1a":
            m["stereochemistry"]["axial"] = fact(
                "Allene shown without an axial descriptor; one terminal carbon is CH2.",
                [{"type": "allene", "assignment": None, "reason": "one terminal carbon bears two hydrogens"}],
                evidence, method="structure_derivation",
                derivation="Inspection of the displayed N-allenyl substituent; no Ra/Sa or P/M label is printed."
            )
            m["allene_features"] = [{
                "fragment_id": "allene-1",
                "carbon_atom_maps": fact("C=C=C fragment", None, evidence, method="structure_derivation",
                                         derivation="Atom maps are not asserted without an atom-mapped source structure."),
                "substitution": fact("N-substituted terminal allene with the other terminus CH2",
                                     {"terminal_1": "N-substituted", "terminal_2": "CH2"}, evidence,
                                     method="structure_derivation", derivation="Read from Table 1/2 structure."),
                "axial_stereochemistry": fact("No axial descriptor shown", None, evidence),
                "local_environment": fact("N-allenyl-2-iodo-N-tosylaniline drawing", label, evidence),
            }]
        else:
            m["stereochemistry"]["axial"] = fact()
        m["stereochemistry"]["mixture"] = fact()
        m["stereochemistry"]["assignment"] = fact()
        for field_name, value in (
            ("chirality_status", chirality),
            ("geometry_status", geometry),
            ("representation_status", representation),
        ):
            m["stereochemistry"][field_name] = fact(
                value, value, evidence, method="inference",
                derivation="Status assigned from the displayed structure plus standard stereo diagnostics; absence of a source descriptor is retained as unspecified rather than racemic."
            )
        if chirality == "achiral":
            m["stereochemistry"]["isomeric_smiles"] = fact(
                {"candidate_smiles": candidate}, diagnostic["canonical_isomeric_smiles"], evidence,
                method="structure_derivation",
                derivation="Canonical isomeric SMILES of the non-stereogenic candidate graph."
            )
        else:
            m["stereochemistry"]["isomeric_smiles"] = fact()
        m["stereochemistry"]["identity_key"] = fact(
            None, diagnostic["stereo_graph_key"], evidence, method="structure_derivation",
            derivation="SHA-256 stereo graph key generated by standard stereo; not a sample identity."
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = fact()
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = fact()
        m["stereochemistry"]["double_bonds"] = fact()
        m["stereochemistry"]["axial"] = fact()
        m["stereochemistry"]["mixture"] = fact()
        m["stereochemistry"]["assignment"] = fact()
        for field_name, value in (
            ("chirality_status", chirality or "unresolved"),
            ("geometry_status", geometry or "unresolved"),
            ("representation_status", representation or "unresolved"),
        ):
            m["stereochemistry"][field_name] = fact(
                value, value, evidence, method="inference",
                derivation="The main-text figure/label was inventoried, but a complete stereochemically safe structure was not asserted from this page."
            )
        m["stereochemistry"]["isomeric_smiles"] = fact()
        m["stereochemistry"]["identity_key"] = fact()
    m["external_properties"] = []
    molecules[molecule_id] = m
    bundle["molecules"].append(m)
    return molecule_id


def ensure_template(key, label, evidence, fixed_description, variables):
    molecule_id = mid(key)
    if molecule_id in molecules:
        return molecule_id
    m = st.new_molecule(molecule_id, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{"document_id": DOC,
                         "compound_label": fact(label, label, evidence),
                         "evidence_ids": list(evidence)}]
    m["names"]["original"] = fact(label, label, evidence)
    m["names"]["english"] = fact()
    m["names"]["abbreviations"] = fact()
    m["names"]["aliases"] = fact()
    m["cas_numbers"] = fact()
    m["representation_type"] = "generic_template"
    m["smiles"] = fact()
    var_rows = []
    for var in variables:
        var_rows.append({
            "label": var,
            "attachment_points": fact(f"Variable {var} attachment is shown graphically", None, evidence,
                                      method="structure_derivation",
                                      derivation="The source provides a graphical attachment but no atom map."),
            "definition": fact(f"{var} as drawn in Scheme 1c/Table 2", f"source-defined variable {var}", evidence),
            "allowed_values": fact(),
            "constraints": fact(),
        })
    m["generic_template"] = {
        "fixed_structure": fact(fixed_description, {"description": fixed_description}, evidence),
        "variables": var_rows,
        "instance_mapping": [],
        "reason": fact("Generic R/Ar-labelled drawing", "Retained as a template because variables are not uniquely expanded in the generic scheme.", evidence),
    }
    for field_name in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                       "isomeric_smiles", "identity_key"):
        m["stereochemistry"][field_name] = fact()
    for field_name, value in (("chirality_status", "unresolved"),
                              ("geometry_status", "unresolved"),
                              ("representation_status", "not_encoded")):
        m["stereochemistry"][field_name] = fact(
            value, value, evidence, method="inference",
            derivation="Generic variables prevent a single stereochemical identity."
        )
    m["identity_note"] = fact("Generic structure with R/Ar variables", "No concrete SMILES asserted.", evidence)
    m["external_properties"] = []
    molecules[molecule_id] = m
    bundle["molecules"].append(m)
    return molecule_id


def named(key, label, evidence, note=None):
    return ensure_molecule(key, label, evidence, note=note)


def participant(pid, molecule_id, role, evidence, *, stage="reaction", amount=None,
                mass=None, volume=None, equivalents=None, mol_percent=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, molecule_id, role, evidence, stage=stage)
    if product_role:
        p["product_role"] = fact(product_role, product_role, evidence)
    values = {"amount": amount, "mass": mass, "volume": volume,
              "equivalents": equivalents, "mol_percent": mol_percent}
    for field_name, value in values.items():
        if value is not None:
            raw = value.get("raw_text") if isinstance(value, dict) else str(value)
            raw = raw or f"{value['number']} {value['unit']}"
            p["charges"][field_name] = fact(raw, value, evidence)
    if equivalents is not None or mol_percent is not None:
        p["equivalence_reference"] = fact(
            {"table_reference": reference}, reference, evidence, method="inference",
            derivation="Loading/equivalents are interpreted relative to the listed limiting substrate for this table row."
        )
    return p


def solvent_block(solvent_ids, evidence, volume=None, ratio=None):
    return {
        "composition": fact([molecules[s]["occurrences"][0]["compound_label"]["normalized_value"] for s in solvent_ids],
                            solvent_ids, evidence),
        "mixture_ratio": fact(ratio, ratio, evidence) if ratio else fact(),
        "total_volume": fact(f"{volume['number']} {volume['unit']}", volume, evidence) if volume else fact(),
        "concentration": fact(),
        "concentration_basis": fact(),
    }


def reaction_step(step_id, participant_ids, evidence, operation, temperature=None, time=None,
                  atmosphere=None):
    s = st.new_step(step_id, evidence)
    s["participant_ids"] = participant_ids
    s["operation"] = fact(operation, operation, evidence) if operation else fact()
    if temperature is not None:
        raw = temperature if isinstance(temperature, str) else f"{temperature['number']} {temperature['unit']}"
        s["temperature"] = fact(raw, temperature, evidence)
    if time is not None:
        raw = time if isinstance(time, str) else f"{time['number']} {time['unit']}"
        s["time"] = fact(raw, time, evidence)
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, "argon" if atmosphere == "Ar" else atmosphere, evidence)
    return s


def outcome(outcome_id, kind, targets, evidence, value=None, raw=None, yield_type=None,
            measurement_basis=None, scope=None, product_role="product"):
    o = st.new_outcome(outcome_id, kind, targets, evidence)
    if product_role:
        o["product_role"] = fact(product_role, product_role, evidence)
    if value is not None or raw is not None:
        o["value"] = fact(raw if raw is not None else value, value, evidence)
    if yield_type:
        normalized_yield_type = "isolated yield" if yield_type == "Isolated Yields." else yield_type
        o["yield_type"] = fact(yield_type, normalized_yield_type, evidence)
    if measurement_basis:
        normalized_basis = "isolated" if measurement_basis == "Isolated Yields." else measurement_basis
        o["measurement_basis"] = fact(measurement_basis, normalized_basis, evidence)
    if scope:
        o["scope"] = fact(scope, scope, evidence)
    return o


def classification(evidence, labels_with_counts=None, not_applicable_raw=None, limitation=None):
    c = st._skeleton("classification", evidence)
    c["standard_id"] = st.STANDARD
    c["mapped_reaction_smiles"] = fact()
    c["bond_changes"] = []
    c["model_votes"] = []
    if labels_with_counts is None:
        c["labels"] = fact(not_applicable_raw, [], evidence, method="inference",
                           derivation="No observed product transformation is classified for an NR/0% row.")
        c["coverage"] = fact(not_applicable_raw, "not_applicable", evidence, method="inference",
                             derivation="The source reports no observed desired reaction.")
        c["limitations"] = ["Negative/zero-yield row; no bond formation inferred from the intended product."]
        return c
    label_items = {}
    change_index = 0
    for (a, b, order), count in labels_with_counts:
        label = st.bond_label(a, b, order)
        item = label_items.setdefault(label["code"], dict(label, supporting_change_ids=[]))
        for _ in range(count):
            change_index += 1
            change_id = f"b{change_index}"
            change = {"change_id": change_id,
                      "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                      "before_order": 0, "after_order": order, "change_type": "formed"}
            c["bond_changes"].append(fact(
                None, change, evidence, method="structure_derivation",
                derivation="Net bond formation interpreted from the displayed reactant/product structures; atom mapping is intentionally not asserted."
            ))
            item["supporting_change_ids"].append(change_id)
    c["labels"] = fact(
        None, list(label_items.values()), evidence, method="structure_derivation",
        derivation="Taxonomy labels generated from the listed net formed bonds, not from pre-existing product bonds."
    )
    c["coverage"] = fact(
        None, "partial", evidence, method="inference",
        derivation="Heavy-atom net formations are recorded from the figure; full atom mapping, coproducts, and implicit-H accounting are incomplete."
    )
    c["limitations"] = [
        limitation or "Partial structure-derived classification; no complete atom mapping is claimed.",
        "Implicit hydrogen changes are included only where the displayed net transformation makes them explicit.",
    ]
    return c


CASCADE_BONDS = [(('C', 'C', 1), 2), (('C', 'N', 1), 1)]


def new_reaction(rid_suffix, evidence, record_type, categories, location_raw,
                 reaction_name=None, transformation=None):
    rid = f"{DOC}:{rid_suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = fact(categories, categories, evidence)
    r["location"] = fact(location_raw, {"source": location_raw}, evidence)
    if reaction_name:
        r["reaction_name"] = fact(reaction_name, reaction_name, evidence)
    if transformation:
        r["transformation_type"] = fact(
            location_raw, transformation, evidence, method="inference",
            derivation="Transformation label assigned from the source title/caption and displayed net reaction."
        )
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def finalize_reaction(r, unresolved):
    r["review"]["unresolved"] = list(dict.fromkeys(unresolved))
    bundle["reactions"].append(r)
    return r["reaction_id"]


# Shared molecules visible or named in the paper.
m_1a = ensure_molecule("1a", "1a", [EV[2], EV[3]],
    note="N-allenyl-2-iodo-N-tosylaniline benchmark substrate; candidate structure transcribed from Tables 1–2.")
m_2a = ensure_molecule("2a", "2a", [EV[2], EV[3]],
    note="4-methyl-N-(2-(phenylethynyl)phenyl)benzenesulfonamide benchmark substrate.")
m_3aa = ensure_molecule("3aa", "3aa", [EV[2], EV[3]],
    note="Unsymmetrical 3,3′-diindolylmethane benchmark product; no stereochemical descriptor shown.")


# Scheme 1c/Table 2 generic template, kept separate from specific experiments.
g1 = ensure_template("generic_1", "allenamide 1", [EV[2], EV[3]],
    "o-Iodo-N-allenyl arylamide skeleton shown with Ar, R1 and R2 variables", ["Ar", "R1", "R2"])
g2 = ensure_template("generic_2", "o-ethynylaniline 2", [EV[2], EV[3]],
    "o-Alkynylaniline skeleton shown with Ar, R3 and R4 variables", ["Ar", "R3", "R4"])
g3 = ensure_template("generic_3", "diindolylmethane 3", [EV[2], EV[3]],
    "Unsymmetrical 3,3′-diindolylmethane skeleton shown with Ar and R1–R4 variables", ["Ar", "R1", "R2", "R3", "R4"])
generic = new_reaction("scheme1c_generic", [EV[2], EV[3]], "generic_scheme", [],
                       "Scheme 1c and Table 2 generic header",
                       "cascade Heck/cyclization", "generic cascade annulation")
generic["participants"] = [
    participant("generic_1", g1, "reactant", [EV[2], EV[3]]),
    participant("generic_2", g2, "reactant", [EV[2], EV[3]]),
    participant("generic_3", g3, "product", [EV[2], EV[3]], product_role="generic product"),
]
generic["classification"] = classification([EV[2], EV[3]], CASCADE_BONDS,
    limitation="The paper states one C–N and two C–C bonds form; generic variables prevent complete mapping.")
generic["notes"] = [fact(
    "The methodology allows formation of one C–N bond, two C–C bonds and two new rings in one-pot.",
    "Author-reported net bond-count claim for the generic method.", [EV[1], EVT[1]])]
finalize_reaction(generic, ["Generic R/Ar variables are not expanded into a concrete mapped reaction."])


# Table 1: every optimization row is an independent experiment.
table1 = [
    (1, "Pd2(dba)3", "Na2CO3", "DMF", 18, 80),
    (2, "Pd2(dba)3·CHCl3", "Na2CO3", "DMF", 30, 80),
    (3, "Pd2(dppp)2", "Na2CO3", "DMF", "NR", 80),
    (4, "Pd(PPh3)4", "Na2CO3", "DMF", 83, 80),
    (5, "Pd(CH3CN)2Cl2", "Na2CO3", "DMF", 23, 80),
    (6, "Pd(PPh3)2Cl2", "Na2CO3", "DMF", 78, 80),
    (7, "Pd(OAc)2", "Na2CO3", "DMF", 43, 80),
    (8, "Pd(PPh3)4", "K3PO4", "DMF", 81, 80),
    (9, "Pd(PPh3)4", "KHCO3", "DMF", 76, 80),
    (10, "Pd(PPh3)4", "Li2CO3", "DMF", 37, 80),
    (11, "Pd(PPh3)4", "K2CO3", "DMF", 88, 80),
    (12, "Pd(PPh3)4", "Cs2CO3", "DMF", 80, 80),
    (13, "Pd(PPh3)4", "NaOAc", "DMF", 82, 80),
    (14, "Pd(PPh3)4", "K2CO3", "MeOH", 42, 80),
    (15, "Pd(PPh3)4", "K2CO3", "DCE", 86, 80),
    (16, "Pd(PPh3)4", "K2CO3", "toluene", 37, 80),
    (17, "Pd(PPh3)4", "K2CO3", "CH3CN", 97, 80),
    (18, "Pd(PPh3)4", "K2CO3", "THF", 28, 80),
    (19, "Pd(PPh3)4", "K2CO3", "CH3CN", 81, 60),
    (20, "Pd(PPh3)4", "K2CO3", "CH3CN", 75, 100),
    (21, None, "K2CO3", "CH3CN", "NR", 80),
]
table1_ids = []
for entry, catalyst_name, base_name, solvent_name, reported, temp in table1:
    evidence = [EV[2]]
    catalyst_id = named("reagent_" + re.sub(r"\W+", "_", catalyst_name or "no_catalyst").strip("_"),
                          catalyst_name or "no palladium catalyst", evidence,
                          note="Named condition component; complete structure not asserted from the table.") if catalyst_name else None
    base_id = named("reagent_" + re.sub(r"\W+", "_", base_name).strip("_"), base_name, evidence)
    solvent_id = named("solvent_" + solvent_name.lower(), solvent_name, evidence)
    r = new_reaction(f"table1_entry_{entry:02d}", evidence, "specific_experiment", ["optimization"],
                     f"Table 1, entry {entry}", "cascade Heck/cyclization", "optimization cascade")
    parts = [
        participant("p_1a", m_1a, "reactant", evidence, amount=q(0.1, "mmol")),
        participant("p_2a", m_2a, "reactant", evidence, amount=q(0.11, "mmol")),
    ]
    if catalyst_id:
        parts.append(participant("p_catalyst", catalyst_id, "catalyst", evidence,
                                 mol_percent=q(5.0, "mol%"), reference=m_1a))
    parts.append(participant("p_base", base_id, "reagent", evidence,
                             equivalents=q(2.0, "equiv"), reference=m_1a))
    parts.append(participant("p_solvent", solvent_id, "solvent", evidence, volume=q(2.0, "mL")))
    successful = reported != "NR"
    if successful:
        parts.append(participant("p_3aa", m_3aa, "product", evidence, product_role="desired product"))
    r["participants"] = parts
    r["reaction_scale"] = fact("1a (0.1 mmol)", {"reference_molecule_id": m_1a,
                                                   "amount": q(0.1, "mmol")}, evidence)
    r["solvent_conditions"] = solvent_block([solvent_id], evidence, q(2.0, "mL"))
    r["steps"] = [reaction_step("reaction", [p["participant_id"] for p in parts if p["role"] != "product"],
                                evidence, "heat in sealed tube", q(temp, "°C"), q(24, "h"), "Ar")]
    r["special_conditions"]["reactor"] = fact("sealed tube", "sealed tube", evidence)
    r["procedure"]["original_text"] = fact(
        {"footnote_a": "Reaction conditions: 1a (0.1 mmol), 2a (0.11 mmol), catalyst (5.0 mol%), base (2.0 equiv.), solvent (2.0 mL), at 80 °C under Ar atmosphere for 24 h, sealed tube.",
         "row_temperature": f"{temp} °C",
         "temperature_footnote": "d 60 °C." if entry == 19 else "e 100 °C." if entry == 20 else None},
        "Table 1 footnote conditions with row-specific catalyst, base, solvent and temperature.", evidence)
    if successful:
        r["outcomes"] = [outcome("yield", "yield", [m_3aa], evidence, q(reported, "%"),
                                 raw=f"{reported}", yield_type="Isolated Yields.",
                                 measurement_basis="Isolated Yields.", scope="desired product 3aa")]
        r["classification"] = classification(evidence, CASCADE_BONDS)
        unresolved = ["Product structure candidate and net classification remain machine records pending expert review."]
    else:
        r["outcomes"] = [outcome("no_reaction", "no_reaction", [], evidence,
                                 "no reaction", raw="NR", product_role="intended product 3aa")]
        r["classification"] = classification(evidence, None, "NR")
        unresolved = ["The table does not report conversion, recovery, detection limit, or analysis method for NR."]
    table1_ids.append(finalize_reaction(r, unresolved))


# Table 2 scope rows. Product code 3xy maps to 1x + 2y; that mapping is retained as inference.
scope_rows = [
    ("3aa", 97, "R1 = Ts"), ("3ba", 78, "R1 = Ms"), ("3ca", 37, "R1 = Ac"), ("3da", 88, "R1 = Boc"),
    ("3ea", 80, "R = OMe"), ("3fa", 85, "R = Me"), ("3ga", 92, "R = F"), ("3ha", 80, "R = Cl"), ("3ia", 70, "R = Br"),
    ("3ja", 73, "R = Me"), ("3ka", 88, "R = F"), ("3la", 83, "R = Cl"), ("3ma", 84, "R = Br"), ("3na", 98, "R = CN"),
    ("3oa", 78, "methyl-substituted allene product drawing"), ("3pa", 92, "pyridine-containing product drawing"),
    ("3ab", 98, "R4 = 4-OMeC6H4"), ("3ac", 95, "R4 = 4-tBuC6H4"), ("3ad", 83, "R4 = 4-FC6H4"),
    ("3ae", 71, "R4 = 4-BrC6H4"), ("3af", 57, "R4 = 3-ClC6H4"), ("3ag", 78, "R4 = 4-NHTsC6H4"),
    ("3ah", 88, "R4 = 3-OTsC6H4"), ("3ai", 87, "R4 = 3,5-OTsC6H4 as printed; exact multiplicity unresolved"),
    ("3aj", 96, "R4 = 3,5-OMeC6H4 as printed; exact multiplicity unresolved"), ("3ak", 0, "R4 = nBu"),
    ("3al", 73, "R4 = TMS"), ("3am", 83, "thiophene-containing product drawing"),
    ("3an", 77, "R = Me"), ("3ao", 53, "R = OMe"), ("3ap", 78, "R = CF3"), ("3aq", 72, "R = F"), ("3ar", 91, "R = Br"),
    ("3as", 50, "methyl-substituted arylamine product drawing"), ("3at", 0, "C6-methyl o-alkynylaniline; desired product not obtained"),
    ("3au", 50, "R3 = CF3CO2 in figure; nitrogen-deprotected product obtained"),
]
scope_ids = []
pd_id = named("reagent_Pd_PPh3_4", "Pd(PPh3)4", [EV[2], EV[3]])
k2co3_id = named("reagent_K2CO3", "K2CO3", [EV[2], EV[3]])
acn_id = named("solvent_ch3cn", "CH3CN", [EV[2], EV[3]])
for product_label, reported, substitution in scope_rows:
    x, y = product_label[1], product_label[2]
    allene_id = m_1a if x == "a" else ensure_molecule(
        f"1{x}", f"1{x}", [EV[3]],
        note=f"Allenamide precursor inferred from product {product_label}; standalone exact structure is not separately printed.",
        inferred_label=True)
    alkyne_id = m_2a if y == "a" else ensure_molecule(
        f"2{y}", f"2{y}", [EV[3]],
        note=f"o-Ethynylaniline precursor inferred from product {product_label}; standalone exact structure is not separately printed.",
        inferred_label=True)
    product_id = ensure_molecule(
        product_label, product_label, [EV[3]],
        note=f"Table 2: {substitution}; yield cell {reported}%." +
             (" No product isolation is asserted for this 0% row." if reported == 0 else ""),
        chirality="unspecified" if product_label in {"3oa"} else None)
    r = new_reaction(f"table2_{product_label}", [EV[3]], "specific_experiment", ["substrate_scope"],
                     f"Table 2, {product_label}", "cascade Heck/cyclization", "substrate-scope cascade")
    parts = [
        participant("p_allenamide", allene_id, "reactant", [EV[3]], amount=q(0.1, "mmol")),
        participant("p_alkynylaniline", alkyne_id, "reactant", [EV[3]], amount=q(0.11, "mmol")),
        participant("p_catalyst", pd_id, "catalyst", [EV[3]], mol_percent=q(5.0, "mol%"), reference=allene_id),
        participant("p_base", k2co3_id, "reagent", [EV[3]], equivalents=q(2.0, "equiv"), reference=allene_id),
        participant("p_solvent", acn_id, "solvent", [EV[3]], volume=q(2.0, "mL")),
    ]
    if reported > 0:
        parts.append(participant("p_product", product_id, "product", [EV[3]], product_role="scope product"))
    r["participants"] = parts
    r["reaction_scale"] = fact("1 (0.1 mmol)", {"reference_molecule_id": allene_id,
                                                  "amount": q(0.1, "mmol")}, [EV[3]])
    r["solvent_conditions"] = solvent_block([acn_id], [EV[3]], q(2.0, "mL"))
    r["steps"] = [reaction_step("reaction", [p["participant_id"] for p in parts if p["role"] != "product"],
                                [EV[3]], "heat in sealed tube", q(80, "°C"), q(24, "h"), "Ar")]
    r["special_conditions"]["reactor"] = fact("sealed tube", "sealed tube", [EV[3]])
    r["procedure"]["original_text"] = fact(
        "1 (0.1 mmol), 2 (0.11 mmol), Pd(PPh3)4 (5.0 mol%), K2CO3 (2.0 equiv.), CH3CN (2.0 mL), at 80 °C under Ar atmosphere for 24 h, sealed tube.",
        "Table 2 footnote conditions", [EV[3]])
    if reported > 0:
        r["outcomes"] = [outcome("yield", "yield", [product_id], [EV[3]], q(reported, "%"),
                                 raw=f"{reported}% yield", yield_type="Isolated Yields.",
                                 measurement_basis="Isolated Yields.", scope=product_label)]
        r["classification"] = classification([EV[3]], CASCADE_BONDS)
        unresolved = ["Scope structures transcribed only where a stereochemically safe concrete candidate could be represented; other exact structures remain unresolved."]
    else:
        r["outcomes"] = [outcome("no_reaction", "no_reaction", [], [EV[3]],
                                 "desired product not obtained", raw=f"{product_label}, 0% yield",
                                 product_role=f"intended product {product_label}")]
        r["classification"] = classification([EV[3]], None, f"{product_label}, 0% yield")
        unresolved = ["The 0% row gives no conversion, recovery, detection limit, or analytical method."]
    if product_label == "3ab":
        r["notes"].append(fact("The structure of product 3ab was identified by X-ray crystallography; CCDC 2039429.",
                               "Author-reported X-ray identity statement.", [EV[2], EV[3]]))
    scope_ids.append(finalize_reaction(r, unresolved))


# Scheme 2: two additional scaffold-forming arrows.  Time, scale and yield basis are not printed.
scheme2_rows = [("4", "5", 90, "isoquinolone"), ("6", "7", 43, "benzofuran")]
scheme2_ids = []
for substrate_label, product_label, reported, scaffold in scheme2_rows:
    substrate_id = ensure_molecule(substrate_label, substrate_label, [EV[3]],
                                   note=f"Scheme 2 {scaffold} precursor; exact candidate SMILES not asserted from the main-text drawing.")
    product_id = ensure_molecule(product_label, product_label, [EV[3]],
                                 note=f"Scheme 2 {scaffold} product, {reported}% yield; no stereochemical descriptor shown.")
    r = new_reaction(f"scheme2_{product_label}", [EV[3]], "specific_experiment", ["substrate_scope"],
                     f"Scheme 2, product {product_label}", "cascade Heck/cyclization",
                     f"extension to {scaffold} scaffold")
    r["participants"] = [
        participant("p_allene", substrate_id, "reactant", [EV[3]]),
        participant("p_2a", m_2a, "reactant", [EV[3]]),
        participant("p_catalyst", pd_id, "catalyst", [EV[3]], mol_percent=q(5.0, "mol%"), reference=substrate_id),
        participant("p_base", k2co3_id, "reagent", [EV[3]], equivalents=q(2.0, "equiv"), reference=substrate_id),
        participant("p_solvent", acn_id, "solvent", [EV[3]]),
        participant("p_product", product_id, "product", [EV[3]], product_role=f"{scaffold} product"),
    ]
    r["solvent_conditions"] = solvent_block([acn_id], [EV[3]])
    r["steps"] = [reaction_step("reaction", [p["participant_id"] for p in r["participants"] if p["role"] != "product"],
                                [EV[3]], "heat under displayed conditions", q(80, "°C"), None, "Ar")]
    r["procedure"]["original_text"] = fact(
        "Pd(PPh3)4 (5.0 mol%), K2CO3 (2.0 equiv.), CH3CN, Ar, 80 °C", "Scheme 2 arrow conditions", [EV[3]])
    r["outcomes"] = [outcome("yield", "yield", [product_id], [EV[3]], q(reported, "%"),
                             raw=f"{reported}% yield", scope=product_label)]
    r["classification"] = classification([EV[3]], CASCADE_BONDS,
        limitation="Net C–C/C–N formations interpreted from the displayed scaffold products; complete mapping is unresolved.")
    scheme2_ids.append(finalize_reaction(r, [
        "Reaction time, substrate scale, solvent volume, workup/purification, and yield measurement basis are not reported in Scheme 2.",
        "Exact substrate/product structures remain graphical and require expert structure review."
    ]))


# Table 3 shared derived molecules and reagents.
m_3al = ensure_molecule("3al", "3al", [EV[3]], note="TMS-substituted diindolylmethane product; candidate connectivity from Table 2/3.")
for key, note, chirality in [
    ("8", "Symmetric diindolylmethane after TMS and sulfonamide removal.", "achiral"),
    ("9", "Oxidized carbonyl derivative of 3aa.", "achiral"),
    ("10", "Alcohol derivative; the newly tetrahedral carbon is shown without R/S assignment.", "unspecified"),
    ("11", "Azido derivative; tetrahedral center shown without R/S assignment.", "unspecified"),
    ("12", "Cyano derivative; tetrahedral center shown without R/S assignment.", "unspecified"),
    ("13", "Allyl derivative; tetrahedral center shown without R/S assignment.", "unspecified"),
    ("14", "Alkynyl derivative; tetrahedral center shown without R/S assignment.", "unspecified"),
]:
    ensure_molecule(key, key, [EV[3]], note=note, chirality=chirality,
                    geometry="not_applicable", representation="not_encoded" if chirality == "unspecified" else "not_applicable")


def add_table3_scaleup(tag, alkyne_label, product_label, reported, mass_g, hours):
    evidence = [EV[3]]
    alkyne_id = m_2a if alkyne_label == "2a" else ensure_molecule(
        alkyne_label, alkyne_label, evidence,
        note="TMS-substituted o-alkynylaniline precursor identified as 2l in Table 3.", inferred_label=False)
    product_id = mid(product_label)
    r = new_reaction(f"table3a_scaleup_{product_label}", evidence, "specific_experiment", ["scale_up"],
                     f"Table 3a, R={'Ph' if product_label == '3aa' else 'TMS'}",
                     "cascade Heck/cyclization", "millimole-scale cascade")
    r["participants"] = [
        participant("p_1a", m_1a, "reactant", evidence),
        participant("p_alkyne", alkyne_id, "reactant", evidence),
        participant("p_catalyst", pd_id, "catalyst", evidence, mol_percent=q(5.0, "mol%"), reference=m_1a),
        participant("p_base", k2co3_id, "reagent", evidence, equivalents=q(2.0, "equiv"), reference=m_1a),
        participant("p_solvent", acn_id, "solvent", evidence),
        participant("p_product", product_id, "product", evidence, product_role="scale-up product"),
    ]
    r["reaction_scale"] = fact("on a millimole scale", {"reference_molecule_id": m_1a,
                                                          "description": "millimole scale"}, evidence)
    r["solvent_conditions"] = solvent_block([acn_id], evidence)
    r["steps"] = [reaction_step("reaction", [p["participant_id"] for p in r["participants"] if p["role"] != "product"],
                                evidence, "heat on millimole scale", q(80, "°C"), q(hours, "h"), "Ar")]
    r["procedure"]["original_text"] = fact(
        f"R = {'Ph' if product_label == '3aa' else 'TMS'}, Pd(PPh3)4 (5.0 mol%), K2CO3 (2.0 equiv.), CH3CN, at 80 °C under Ar for {hours} h.",
        "Table 3a conditions", evidence)
    r["outcomes"] = [
        outcome("yield", "yield", [product_id], evidence, q(reported, "%"), raw=f"{reported}% yield", scope=product_label),
        outcome("product_amount", "product_amount", [product_id], evidence, q(mass_g, "g"), raw=f"{mass_g} g", scope=product_label),
    ]
    r["classification"] = classification(evidence, CASCADE_BONDS)
    return finalize_reaction(r, [
        "Exact millimole charge, solvent volume, workup/purification, and yield measurement basis are not reported in the main text.",
        "Product mass is reported, but the source does not explicitly label the yield method as isolated."
    ])


table3_ids = [
    add_table3_scaleup("ph", "2a", "3aa", 90, 1.14, 48),
    add_table3_scaleup("tms", "2l", "3al", 72, 0.45, 24),
]


def table3_deriv(rid, reactant_key, product_key, yield_value, reagent_specs, solvent_specs,
                 temp_steps, atmosphere, raw_conditions, formed, transformation):
    evidence = [EV[3]]
    reactant_id, product_id = mid(reactant_key), mid(product_key)
    panel = rid[len("table3")] if rid.startswith("table3") and len(rid) > len("table3") else "?"
    r = new_reaction(rid, evidence, "specific_experiment", ["derivatization"],
                     f"Table 3{panel}", None, transformation)
    parts = [participant("p_start", reactant_id, "reactant", evidence),
             participant("p_product", product_id, "product", evidence, product_role="derivatization product")]
    for index, spec in enumerate(reagent_specs, 1):
        key, label, role, load_kind, load = spec
        reagent_id = named(key, label, evidence)
        kwargs = {load_kind: load, "reference": reactant_id} if load_kind else {}
        parts.append(participant(f"p_reagent_{index}", reagent_id, role, evidence, **kwargs))
    solvent_ids = []
    for index, (key, label) in enumerate(solvent_specs, 1):
        solvent_id = named(key, label, evidence)
        solvent_ids.append(solvent_id)
        parts.append(participant(f"p_solvent_{index}", solvent_id, "solvent", evidence))
    r["participants"] = parts
    ratio = "3:2" if len(solvent_specs) == 2 else None
    r["solvent_conditions"] = solvent_block(solvent_ids, evidence, ratio=ratio)
    reaction_participant_ids = [p["participant_id"] for p in parts if p["role"] != "product"]
    r["steps"] = []
    for index, (operation, temperature, time_value) in enumerate(temp_steps, 1):
        r["steps"].append(reaction_step(f"reaction_{index}", reaction_participant_ids, evidence,
                                        operation, temperature, time_value,
                                        atmosphere if index == len(temp_steps) else None))
    r["procedure"]["original_text"] = fact(raw_conditions, raw_conditions, evidence)
    r["outcomes"] = [outcome("yield", "yield", [product_id], evidence, q(yield_value, "%"),
                             raw=f"{yield_value}% yield", scope=product_key)]
    r["classification"] = classification(evidence, formed)
    unresolved = [
        "Yield measurement basis, reaction scale, workup, and purification are not reported in the main-text Table 3.",
        "Full atom mapping and coproduct identities are unresolved; classification is a partial net transformation."
    ]
    if rid == "table3d_reduction":
        unresolved.append("The source states 0 °C to rt over 16 h but does not allocate the 16 h between addition, warming, and hold; separate steps retain this ambiguity.")
    return finalize_reaction(r, unresolved)


table3_ids += [
    table3_deriv("table3b_deprotection", "3al", "8", 77,
        [("reagent_NaOH", "NaOH", "reagent", "equivalents", q(6.0, "equiv"))],
        [("solvent_etoh", "EtOH")], [("heat", q(80, "°C"), q(18, "h"))], None,
        "NaOH (6.0 equiv.), EtOH, at 80 °C for 18 h.",
        [(('N', 'H', 1), 2), (('C', 'H', 1), 1)], "TMS and sulfonamide deprotection"),
    table3_deriv("table3c_oxidation", "3aa", "9", 93,
        [("reagent_SeO2", "SeO2", "reagent", "equivalents", q(20.0, "equiv"))],
        [("solvent_dioxane", "1,4-dioxane")], [("heat", q(100, "°C"), q(72, "h"))], None,
        "SeO2 (20.0 equiv.), 1,4-dioxane, at 100 °C for 72 h.",
        [(('C', 'O', 2), 1)], "oxidation of diindolylmethane to carbonyl derivative"),
    table3_deriv("table3d_reduction", "9", "10", 85,
        [("reagent_NaBH4", "NaBH4", "reagent", "equivalents", q(2.0, "equiv"))],
        [("solvent_meoh", "MeOH")],
        [("addition at 0 °C", q(0, "°C"), None), ("warm to room temperature; total reported time", "rt", q(16, "h"))], None,
        "NaBH4 (2.0 equiv.), MeOH, 0 °C to rt, 16 h.",
        [(('C', 'H', 1), 1), (('O', 'H', 1), 1)], "carbonyl reduction to secondary alcohol"),
    table3_deriv("table3e_azidation", "10", "11", 87,
        [("reagent_InBr3", "InBr3", "catalyst", "mol_percent", q(20.0, "mol%")),
         ("reagent_TMSN3", "TMSN3", "reagent", "equivalents", q(3.0, "equiv"))],
        [("solvent_dcm", "DCM")], [("stir at room temperature", "rt", q(10, "h"))], None,
        "InBr3 (20.0 mol%), TMSN3 (3.0 equiv.), DCM, rt, 10 h.",
        [(('C', 'N', 1), 1)], "direct azidation of alcohol"),
    table3_deriv("table3f_cyanation", "10", "12", 31,
        [("reagent_InBr3", "InBr3", "catalyst", "mol_percent", q(20.0, "mol%")),
         ("reagent_TMSCN", "TMSCN", "reagent", "equivalents", q(3.0, "equiv"))],
        [("solvent_dcm", "DCM")], [("stir at room temperature", "rt", q(10, "h"))], None,
        "InBr3 (20.0 mol%), TMSCN (3.0 equiv.), DCM, rt, 10 h.",
        [(('C', 'C', 1), 1)], "direct cyanation of alcohol"),
    table3_deriv("table3g_allylation", "10", "13", 60,
        [("reagent_InCl3", "InCl3", "catalyst", "mol_percent", q(20.0, "mol%")),
         ("reagent_TMSBr", "TMSBr", "reagent", "mol_percent", q(20.0, "mol%")),
         ("reagent_allylTMS", "allyltrimethylsilane (drawn)", "reagent", None, None)],
        [("solvent_hexane", "hexane"), ("solvent_dcm", "DCM")],
        [("heat", q(50, "°C"), q(2, "h"))], "Ar",
        "InCl3 (20.0 mol%), TMSBr (20.0 mol%), hexane/DCM (3:2), at 50 °C under Ar for 2 h.",
        [(('C', 'C', 1), 1)], "direct allylation of alcohol"),
    table3_deriv("table3h_alkynylation", "10", "14", 72,
        [("reagent_InCl3", "InCl3", "catalyst", "mol_percent", q(20.0, "mol%")),
         ("reagent_TMSBr", "TMSBr", "reagent", "equivalents", q(0.6, "equiv")),
         ("reagent_phenylethynylTMS", "trimethyl(phenylethynyl)silane (drawn)", "reagent", None, None)],
        [("solvent_hexane", "hexane"), ("solvent_dcm", "DCM")],
        [("heat", q(50, "°C"), q(4, "h"))], "Ar",
        "InCl3 (20.0 mol%), TMSBr (0.6 equiv.), hexane/DCM (3:2), at 50 °C under Ar for 4 h.",
        [(('C', 'C', 1), 1)], "direct alkynylation of alcohol"),
]


# Scheme 3 is an author-proposed catalytic cycle, not an experimental proof.
mechanism = new_reaction("scheme3_proposed_mechanism", [EV[3], EV[4]], "proposed_mechanism", [],
                         "Scheme 3 and continuation on page 4", "proposed reaction mechanism",
                         "proposed Pd(0)/Pd(II) cascade cycle")
mechanism["participants"] = [
    participant("p_1a", m_1a, "reactant", [EV[3], EV[4]]),
    participant("p_2a", m_2a, "reactant", [EV[3], EV[4]]),
    participant("p_3aa", m_3aa, "product", [EV[3], EV[4]], product_role="proposed-cycle product"),
    participant("p_base", k2co3_id, "reagent", [EV[3], EV[4]]),
]
mechanism["classification"] = classification([EV[3], EV[4]], CASCADE_BONDS,
    limitation="Bond formations are those drawn/labeled in the proposed cycle; this record is not an observed mechanistic experiment.")
mechanism["notes"] = [
    fact("Oxidative addition gives I; intramolecular Heck reaction gives π-allylpalladium II; alkyne coordination/ligand exchange gives III; base gives IV; C–C reductive elimination affords 3aa.",
         "Author-proposed sequence, not direct experimental proof.", [EV[3], EV[4]]),
]
mechanism_id = finalize_reaction(mechanism, [
    "Scheme 3 is explicitly proposed from literature precedents; no independent kinetic, isotope, or trapping experiment is reported on pages 1–6."
])


# Non-synthetic source context and metadata.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:title", "document_id": DOC, "kind": "bibliographic_metadata",
     "evidence_ids": [EV[1], EVT[1]], "title": TITLE, "doi": "10.1039/d0qo01539c",
     "note": "True article title transcribed from page 1."},
    {"context_id": f"{DOC}:context:figure1", "document_id": DOC, "kind": "non_experimental_figure",
     "evidence_ids": [EV[1]], "note": "Figure 1 contains representative bioactive/chromophore examples, not experiments performed in this article."},
    {"context_id": f"{DOC}:context:scheme1_prior_work", "document_id": DOC, "kind": "literature_context",
     "evidence_ids": [EV[2]], "note": "Scheme 1a–b summarizes previous approaches; it is not counted as a current-paper experiment. Scheme 1c is represented separately as a generic scheme."},
    {"context_id": f"{DOC}:context:xray_3ab", "document_id": DOC, "kind": "structure_characterization",
     "evidence_ids": [EV[2], EV[3]], "compound_label": "3ab", "ccdc": "2039429",
     "note": "Author states that 3ab was identified by X-ray crystallography; no separate synthesis record is created for the crystallographic observation."},
    {"context_id": f"{DOC}:context:esi", "document_id": DOC, "kind": "supplement_availability_statement",
     "evidence_ids": [EV[1]], "note": "Page 1 states ESI is available, but no SI file was explicitly supplied to this session; association remains unestablished."},
    {"context_id": f"{DOC}:context:references", "document_id": DOC, "kind": "bibliography",
     "evidence_ids": [EV[4], EV[5], EV[6]], "note": "Pages 4–6 contain conclusions and references. Cited literature is not imported as current-paper experimental evidence."},
]


# Raw visual inventory, separate from the normalized submission.
raw_inventory = {
    "document_id": DOC,
    "title": TITLE,
    "doi": "10.1039/d0qo01539c",
    "page_range": [1, 6],
    "source_images_visually_inspected": [f"runs/current/assets/{DOC}/page{p}.jpg" for p in range(1, 7)],
    "table1_rows": [
        {"entry": e, "catalyst": c, "base": b, "solvent": s, "yield_cell": y, "temperature_C": t}
        for e, c, b, s, y, t in table1
    ],
    "table2_rows": [
        {"product": p, "yield_percent": y, "substitution_or_visual_note": n}
        for p, y, n in scope_rows
    ],
    "scheme2_arrows": [{"product": p, "yield_percent": y, "scaffold": s}
                       for _, p, y, s in scheme2_rows],
    "table3_experiments": ["3aa scale-up", "3al scale-up", "3al→8", "3aa→9", "9→10",
                           "10→11", "10→12", "10→13", "10→14"],
    "record_counts": {"generic_scheme": 1, "optimization": 21, "substrate_scope": 36,
                      "scaffold_extension": 2, "scale_up_and_derivatization": 9,
                      "proposed_mechanism": 1, "total": 70},
    "limitations": [
        "No SI was provided or associated.",
        "Complex scope structures are retained as labels/graphical identities unless a safe candidate SMILES was checked.",
        "Schema validation does not certify chemistry or complete experimental reproducibility."
    ],
}


# Per-page coverage required by the session.
page2_records = [f"{DOC}:scheme1c_generic"] + table1_ids
page3_records = scope_ids + scheme2_ids + table3_ids + [mechanism_id]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title, abstract, introduction, Figure 1 and ESI footnote",
             "finding": "Bibliographic title/DOI and method overview; Figure 1 is non-experimental context; ESI is mentioned but not supplied.",
             "evidence_ids": [EV[1], EVT[1]]}], "record_ids": [],
         "unresolved": ["No SI file is explicitly associated; Figure 1 compound structures are context, not article experiments."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Scheme 1a–c, upper-left", "finding": "Prior-work panels a–b stored as literature context; current-work panel c stored as one generic scheme.",
             "evidence_ids": [EV[2]]},
            {"region": "Table 1, upper-right", "finding": "All 21 optimization rows recorded independently, including two NR controls and temperature footnotes.",
             "evidence_ids": [EV[2], EVT[2]]},
            {"region": "Results/discussion columns", "finding": "Narrative cross-check of benchmark identities, optimum entry 17, scope summary, Scheme 2 and Table 3 references.",
             "evidence_ids": [EV[2], EVT[2]]}], "record_ids": page2_records,
         "unresolved": ["NR rows lack conversion/recovery/detection method.", "Main text supplies no workup or purification procedures."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Table 2, left column", "finding": "All 36 scope rows inventoried and recorded; 3ak and 3at are 0% desired-product rows; isolated-yield footnote applied only to positive rows.",
             "evidence_ids": [EV[3], EVT[3]]},
            {"region": "Scheme 2, lower-left", "finding": "Two scaffold-extension arrows recorded separately (5, 90%; 7, 43%).",
             "evidence_ids": [EV[3]]},
            {"region": "Table 3, upper-right", "finding": "Two scale-up experiments and seven derivatizations recorded separately with all printed reagent loadings, temperatures, times, atmosphere and solvent ratios.",
             "evidence_ids": [EV[3], EVT[3]]},
            {"region": "Scheme 3, lower-right", "finding": "Author-proposed catalytic cycle recorded as proposed_mechanism, not an observed mechanistic experiment.",
             "evidence_ids": [EV[3], EVT[3]]}], "record_ids": page3_records,
         "unresolved": ["Scheme 2 omits time, scale and yield basis.", "Table 3 omits scale for derivatizations, workup/purification and yield basis.", "Several scope structures remain graphical/unresolved rather than assigned speculative SMILES."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "top-left mechanism continuation and conclusions", "finding": "Prose continuation of proposed Scheme 3 and general conclusion; no new experiment, control, kinetic/isotope study or reaction arrow.",
             "evidence_ids": [EV[4], EVT[4]]},
            {"region": "references 1–5", "finding": "Cited literature only; not imported as current-paper experiments.",
             "evidence_ids": [EV[4]]}], "record_ids": [mechanism_id],
         "unresolved": ["Mechanistic intermediates are proposals based on literature precedents, not directly verified observations."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 5–14 across both columns", "finding": "Bibliography continuation only; no current-paper reaction arrow, table, experimental procedure or independent spectral peak list.",
             "evidence_ids": [EV[5], EVT[5]]}], "record_ids": [],
         "unresolved": ["Cited papers' experimental details are outside the supplied source and were not extracted."]},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 14–22 across both columns", "finding": "Bibliography continuation only; no current-paper experiment or characterization table.",
             "evidence_ids": [EV[6], EVT[6]]}], "record_ids": [],
         "unresolved": ["References 18–22 motivate Table 3 transformations but do not supply this paper's missing procedures."]},
    ],
}


assert len(bundle["reactions"]) == 70, len(bundle["reactions"])
assert len(coverage["pages"]) == 6 and all(p["regions"] for p in coverage["pages"])

st.write(RUN / "raw_model/pages_0001_0006_extraction.json", raw_inventory)
st.write(RUN / "raw_model/source_mentions_1a_2a_3aa_3al_3at_10.json",
         st.source_mentions(st.read(RUN / "dataset.json"), ["1a", "2a", "3aa", "3al", "3at", "10"], 320))
st.write(RUN / "coverage/pages_0001_0006.json", coverage)
st.write(RUN / "submitted.json", bundle)
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"])}, ensure_ascii=False))
