#!/usr/bin/env python3
"""Build the machine v3.1 submission from visually inspected pages 1--4."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_782b0f765a8c7b3d"
E = {n: f"{DOC}:page{n}:text" for n in range(1, 5)}
I = {n: f"{DOC}:page{n}:image" for n in range(1, 5)}


def q(number, unit):
    return {"number": number, "unit": unit}


def direct(raw, normalized=None, ev=()):
    return st.fact(raw, raw if normalized is None else normalized, ev)


def derived(normalized, ev, why, raw=None, method="structure_derivation"):
    return st.fact(raw, normalized, ev, method=method, derivation=why)


def missing(ev=()):
    return st.fact(evidence=ev)


def set_stereo(m, ev, chirality="unresolved", geometry="unresolved",
               representation="not_encoded", assignment=None):
    why = "Stereochemistry status assigned conservatively from the supplied source drawing; no unreported configuration is inferred."
    m["stereochemistry"]["chirality_status"] = derived(chirality, ev, why, method="inference")
    m["stereochemistry"]["geometry_status"] = derived(geometry, ev, why, method="inference")
    m["stereochemistry"]["representation_status"] = derived(representation, ev, why, method="inference")
    if assignment is not None:
        m["stereochemistry"]["assignment"] = direct(assignment, assignment, ev)


def molecule(mid, label, ev, name=None, smiles=None, allene=False,
             chirality="unresolved", geometry="unresolved", assignment=None):
    m = st.new_molecule(f"{DOC}:mol:{mid}", DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, label, ev)
    m["names"]["original"] = direct(name or label, name or label, ev)
    m["names"]["english"] = direct(name or label, name or label, ev)
    m["cas_numbers"] = missing(ev)
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, ev,
            "Candidate connectivity SMILES transcribed from the supplied structure drawing; RDKit parsing is a graph diagnostic, not source identity verification.",
            raw=f"source structure drawing for {label}",
        )
        diag = st.stereo_inspect(smiles)
        representation_status = "partial" if any(x["specified"] == "Unspecified" for x in diag["potential_stereo"]) else "encoded"
        if representation_status == "partial" and geometry == "not_applicable":
            geometry = "unspecified"
        set_stereo(m, ev, chirality, geometry, representation_status)
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], ev,
            "RDKit canonical isomeric SMILES generated from the candidate source-derived graph; no additional stereochemistry was added.",
        )
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], ev,
            "RDKit stereo graph key generated from the candidate source-derived graph; it is not a sample identity.",
        )
        if allene:
            reason = ("The depicted allene has a terminal CH2 group or two identical substituents on one terminal carbon, "
                      "so its allene axis is not stereogenic for this connectivity candidate.")
            m["stereochemistry"]["axial"] = derived(
                [{"axis": "C=C=C", "configuration": "not_applicable", "reason": reason}], ev, reason
            )
            m["allene_features"] = [{
                "fragment_id": f"{mid}:allene1",
                "carbon_atom_maps": derived(None, ev, "No atom-mapped source structure is available."),
                "substitution": direct(f"source drawing for {label}", "source-depicted substitution retained in candidate SMILES", ev),
                "axial_stereochemistry": derived("not_applicable", ev, reason),
                "local_environment": direct(f"source drawing for {label}", "aryl/allenoate-substituted allene", ev),
            }]
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(ev)
        set_stereo(m, ev, chirality, geometry, "not_encoded", assignment)
    m["identity_note"] = direct(
        f"Identity retained as source label {label}; full structure is {'a machine candidate' if smiles else 'not safely encoded from the compact drawing'}.",
        f"Source label {label}; machine record needs expert structure review.", ev,
    )
    return m


def generic_molecule(mid, label, ev, fixed, variables, reason, chirality="unresolved"):
    m = molecule(mid, label, ev, name=label, chirality=chirality, geometry="unresolved")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": direct(fixed, {"description": fixed}, ev),
        "variables": [],
        "instance_mapping": [],
        "reason": direct(reason, reason, ev),
    }
    for label_, definition, allowed in variables:
        m["generic_template"]["variables"].append({
            "label": label_,
            "attachment_points": direct(f"source attachment point for {label_}", f"source-defined attachment point for {label_}", ev),
            "definition": direct(definition, definition, ev),
            "allowed_values": direct(allowed, allowed, ev),
            "constraints": missing(ev),
        })
    return m


def participant(pid, mid, role, ev, product_role=None):
    p = st.new_participant(pid, f"{DOC}:mol:{mid}", role, ev)
    if product_role:
        p["product_role"] = direct(product_role, product_role, ev)
    return p


def bond_classification(labels, ev, limitation):
    specs = {"C-C": (2, "C", "C"), "C-P": (6, "C", "P")}
    changes, out = [], []
    for i, code in enumerate(labels, 1):
        number, a, b = specs[code]
        cid = f"b{i}"
        change = {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": 0, "after_order": 1, "change_type": "formed"}
        changes.append(derived(change, ev, limitation))
        out.append({"category_number": number, "code": code, "supporting_change_ids": [cid]})
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(ev),
        "bond_changes": changes,
        "labels": derived(out, ev, limitation),
        "coverage": derived("partial" if labels else "unresolved", ev, limitation, method="inference"),
        "limitations": [limitation, "No complete atom mapping is asserted; hydrogens and omitted partners remain unresolved."],
        "model_votes": [],
    }


def no_classification(ev, reason, not_applicable=False):
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(ev),
        "bond_changes": [],
        "labels": derived([], ev, reason, method="inference"),
        "coverage": derived("not_applicable" if not_applicable else "unresolved", ev, reason, method="inference"),
        "limitations": [reason], "model_votes": [],
    }


def reaction(suffix, ev, record_type, categories, location, name, transformation,
             parts, products=(), conditions=True, labels=(), unresolved=(), notes=(), general=()):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = direct(str(categories), list(categories), ev)
    r["location"] = direct(location, {"description": location}, ev)
    r["general_scheme_ids"] = [f"{DOC}:{x}" for x in general]
    r["reaction_name"] = direct(name, name, ev)
    r["transformation_type"] = direct(transformation, transformation, ev)
    r["participants"] = [participant(pid, mid, role, ev, prod_role) for pid, mid, role, prod_role in parts]
    if conditions:
        r["solvent_conditions"]["composition"] = direct("Xylenes", [f"{DOC}:mol:xylenes"], ev)
        s = st.new_step(f"{suffix}:reaction", ev)
        s["participant_ids"] = [p[0] for p in parts if p[2] != "product"]
        s["operation"] = direct("Mo(CO)6, xylenes, reflux", "heat with Mo(CO)6 in xylenes at reflux", ev)
        s["temperature"] = direct("reflux", "reflux", ev)
        r["steps"] = [s]
    if labels:
        r["classification"] = bond_classification(labels, ev, "Net bond types inferred from the depicted cycloadduct connectivity and the authors' cycloaddition/regiochemistry discussion; exact atom identities remain unresolved.")
    else:
        r["classification"] = no_classification(ev, "The compact source drawing does not support a safe complete bond-change assignment.")
    r["notes"] = [direct(x, x, ev) for x in notes]
    base_unresolved = [
        "Reactant/reagent quantities, reaction time, atmosphere, pressure, and workup/purification details are not reported on the supplied main-text pages.",
        "No SI file was supplied; the page-4 availability statement is not an SI association.",
    ]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": base_unresolved + list(unresolved)}
    r["completeness"] = "partial"
    return r


def add_yield(r, oid, mids, value, ev, basis=None, scope=None):
    o = st.new_outcome(oid, "yield", [f"{DOC}:mol:{m}" for m in mids], ev)
    o["product_role"] = direct("product", "product", ev)
    o["value"] = direct(f"{value}%", q(value, "%"), ev)
    o["yield_type"] = direct("yield", "reported yield", ev)
    if basis:
        o["measurement_basis"] = direct(basis, basis, ev)
    if scope:
        o["scope"] = direct(scope, scope, ev)
    r["outcomes"].append(o)


def add_ratio(r, oid, mids, ratio, ev, method, order):
    o = st.new_outcome(oid, "selectivity", [f"{DOC}:mol:{m}" for m in mids], ev)
    o["product_role"] = direct("regioisomer mixture", "regioisomer mixture", ev)
    o["value"] = direct(ratio, ratio, ev)
    o["measurement_basis"] = direct(method, method, ev)
    o["scope"] = direct(f"ordered as {order}", {"metric": "regioisomer ratio", "order": order}, ev)
    r["outcomes"].append(o)


bundle = json.loads((RUN / "submitted.json").read_text())
molecules = []

# Current-paper starting complexes, allenes, products, intermediate and agents.
for mid, label, ev, ar in [
    ("1a", "1a", [I[2], E[1]], "2-thienyl"), ("1b", "1b", [I[2], I[3]], "2-naphthyl"),
    ("1c", "1c", [I[3]], "Ph"), ("1d", "1d", [I[3]], "p-MeOPh"), ("1e", "1e", [I[3]], "p-CF3Ph")]:
    molecules.append(molecule(mid, label, ev, f"1-iminylphosphirane Mo complex {label} (imino-carbon aryl = {ar})",
                              chirality="unspecified", geometry="unspecified"))

allene_data = [
    ("2a", "1,1-diphenylallene", "C=C=C(c1ccccc1)c1ccccc1", [I[2], I[3]]),
    ("2b", "1,1-bis(4-methoxyphenyl)allene", "C=C=C(c1ccc(OC)cc1)c1ccc(OC)cc1", [I[3]]),
    ("2c", "1,1-bis(4-fluorophenyl)allene", "C=C=C(c1ccc(F)cc1)c1ccc(F)cc1", [I[3]]),
    ("2d", "1-phenylallene", "C=C=Cc1ccccc1", [I[3]]),
    ("2e", "ethyl allenoate", "CCOC(=O)C=C=C", [I[3]]),
    ("2f", "1-(4-fluorophenyl)allene", "C=C=Cc1ccc(F)cc1", [I[3]]),
    ("2g", "1-(4-methoxyphenyl)allene", "C=C=Cc1ccc(OC)cc1", [I[3]]),
]
for mid, name, smi, ev in allene_data:
    molecules.append(molecule(mid, mid, ev, name, smi, allene=True, chirality="achiral", geometry="not_applicable"))

product_ids = ["3aa", "4aa", "5aa", "3ba", "6ba", "7", "3ca", "3da", "3ea", "3cb", "3cc",
               "3cd", "3cd_prime", "3dd", "3dd_prime", "4ad", "4ad_prime", "6bd", "6bd_prime",
               "6bf", "6bf_prime", "6bg", "6bg_prime"]
labels = {x: x.replace("_prime", "′") for x in product_ids}
for mid in product_ids:
    page = 2 if mid in {"3aa", "4aa", "5aa", "3ba", "6ba"} else 3
    if mid == "6bd": ev = [I[3], I[4]]
    elif mid == "7": ev = [I[3], E[1]]
    else: ev = [I[page]]
    assignment = None if mid == "7" else "Compact cycloadduct drawing may contain perspective/relative stereochemical cues, but no R/S or absolute assignment is printed; exact stereochemical interpretation remains unresolved."
    molecules.append(molecule(mid, labels[mid], ev, f"source-labelled cycloadduct {labels[mid]}",
                              chirality="unresolved", geometry="unspecified", assignment=assignment))

molecules.append(molecule("MoCO6", "Mo(CO)6", [I[2], I[3]], "molybdenum hexacarbonyl",
                          chirality="not_applicable", geometry="not_applicable"))
molecules.append(generic_molecule("xylenes", "Xylenes", [I[2], I[3]], "named solvent mixture: xylenes", [],
                                  "The source reports 'Xylenes' without specifying isomer composition; no unique concrete molecule is invented.", "not_applicable"))

# Generic/current and prior-work scheme identities.
molecules += [
    generic_molecule("generic_1c_1e", "1c-1e", [I[3]], "1-iminylphosphirane Mo complex with variable Ar",
                     [("Ar", "aryl group on imino carbon", ["Ph", "p-MeOPh", "p-CF3Ph"])], "Table 1 depicts a variable-Ar series."),
    generic_molecule("generic_2a_2e", "2a-2e", [I[3]], "allene with variable R1/R2",
                     [("R1", "allene substituent R1", ["Ph", "p-MeOPh", "p-FPh", "CO2Et"]),
                      ("R2", "allene substituent R2", ["Ph", "p-MeOPh", "p-FPh", "H"])], "Table 1 depicts a variable allene series."),
    generic_molecule("generic_3", "3", [I[3]], "[2+2] cycloadduct regioisomer 3 with variable Ar/R1/R2", [], "Generic product drawing in Table 1."),
    generic_molecule("generic_3prime", "3′", [I[3]], "alternative [2+2] cycloadduct regioisomer 3′", [], "Generic alternative product drawing in Table 1."),
    molecule("lit1a_start", "Scheme 1a starting complex", [I[2]], "1-acylphosphirane tungsten complex", chirality="unspecified", geometry="unspecified"),
    molecule("lit1a_diene", "Scheme 1a diene", [I[2]], "2,3-dimethyl-1,3-butadiene", "C=C(C)C(C)=C", chirality="achiral", geometry="not_applicable"),
    molecule("lit1a_product", "Scheme 1a product", [I[2]], "prior-work annulated phosphorus heterocycle", chirality="unresolved", geometry="unspecified"),
    molecule("lit1b_start", "Scheme 1b starting complex", [I[2]], "prior-work 1-iminylphosphirane tungsten complex", chirality="unspecified", geometry="unspecified"),
    generic_molecule("lit1b_alkene", "Scheme 1b alkene", [I[2]], "alkene with variable R", [("R", "alkene substituent", "unreported in overview")], "Scheme 1b is a generic prior-work overview."),
    molecule("lit1b_product", "Scheme 1b product", [I[2]], "prior-work 7-aza-1-phosphanorbornane complex", chirality="unresolved", geometry="unspecified"),
    molecule("lit1c_start", "Scheme 1c starting complex", [I[2]], "prior-work 1-iminylphosphirane molybdenum complex", chirality="unspecified", geometry="unspecified"),
    generic_molecule("lit1c_alkene", "Scheme 1c alkene", [I[2]], "alkene with variable R", [("R", "alkene substituent", "unreported in overview")], "Scheme 1c is a generic prior-work overview."),
    molecule("lit1c_product", "Scheme 1c product", [I[2]], "prior-work annulated heterodecalin complex", chirality="unresolved", geometry="unspecified"),
]

bundle["molecules"] = molecules
reactions = []

# Three explicitly cited prior-work transformations in Scheme 1.
for suffix, loc, parts, name, trans, note in [
    ("scheme1a_literature", "PDF page 2, Scheme 1a", [("s1a_start", "lit1a_start", "reactant", None), ("s1a_diene", "lit1a_diene", "reactant", None), ("s1a_product", "lit1a_product", "product", "reported overview product")], "Important ring expansion of phosphirane complex (prior work)", "rearrangement followed by cycloaddition", "Scheme 1a prints 130 °C; solvent, time, amounts and yield are absent from this overview."),
    ("scheme1b_literature", "PDF page 2, Scheme 1b", [("s1b_start", "lit1b_start", "reactant", None), ("s1b_alkene", "lit1b_alkene", "reactant", None), ("s1b_product", "lit1b_product", "product", "reported overview product")], "Important ring expansion of phosphirane complex (prior work)", "rearrangement followed by alkene cycloaddition", "Scheme 1b prints 130 °C; solvent, time, amounts and yield are absent from this overview."),
    ("scheme1c_literature", "PDF page 2, Scheme 1c", [("s1c_start", "lit1c_start", "reactant", None), ("s1c_alkene", "lit1c_alkene", "reactant", None), ("s1c_moco6", "MoCO6", "reagent", None), ("s1c_product", "lit1c_product", "product", "reported overview product")], "Important ring expansion of phosphirane complex (prior work)", "carbonylative ring expansion and dearomatic cycloaddition", "Scheme 1c prints 140 °C and Mo(CO)6; solvent, time, amounts and yield are absent from this overview."),
]:
    r = reaction(suffix, [I[2]], "literature_example", ["other"], loc, name, trans, parts,
                 conditions=False, notes=[note], unresolved=["Prior-work overview lacks sufficient experimental detail and safe exact product encoding."])
    temp = 140 if suffix.endswith("1c_literature") else 130
    s = st.new_step(f"{suffix}:reaction", [I[2]])
    s["participant_ids"] = [p[0] for p in parts if p[2] != "product"]
    s["operation"] = direct(f"{temp} °C" + (", Mo(CO)6" if temp == 140 else ""), "heat prior-work substrates", [I[2]])
    s["temperature"] = direct(f"{temp} °C", q(temp, "°C"), [I[2]])
    r["steps"] = [s]
    reactions.append(r)

# Scheme 2: each substrate combination is one experiment with multiple product outcomes.
r = reaction("scheme2_1a_2a", [I[2], E[1]], "specific_experiment", ["other"], "PDF page 2, Scheme 2 (upper reaction); discussed on page 1",
             "Cycloadditions of 1-iminylphosphirane complexes with allenes", "cascade carbonylative ring expansion and competing [2+2]/[4+2] cycloadditions",
             [("rx_s2a_1a", "1a", "reactant", None), ("rx_s2a_2a", "2a", "reactant", None), ("rx_s2a_moco6", "MoCO6", "reagent", None), ("rx_s2a_xylenes", "xylenes", "solvent", None),
              ("rx_s2a_3aa", "3aa", "product", "[2+2] product"), ("rx_s2a_4aa", "4aa", "product", "dearomative cycloadduct"), ("rx_s2a_5aa", "5aa", "product", "further cycloadduct/isomer identified by NMR comparison")],
             labels=["C-P", "C-C"], unresolved=["Exact all-product atom mapping and the relationship between 4aa and 5aa are not safely resolved from the main-text drawings."],
             notes=["Scheme 2 reports separate yields and 31P NMR shifts for 3aa, 4aa and 5aa; these are outcomes of one substrate combination, not three independently charged experiments."])
for mid, val in [("3aa", 31), ("4aa", 12), ("5aa", 8)]: add_yield(r, f"s2a_y_{mid}", [mid], val, [I[2], E[1]], "isolated", f"individual product {mid}")
reactions.append(r)

r = reaction("scheme2_1b_2a", [I[2], E[1]], "specific_experiment", ["other"], "PDF page 2, Scheme 2 (lower reaction); discussed on pages 1-2",
             "Cycloadditions of 1-iminylphosphirane complexes with allenes", "cascade carbonylative ring expansion and competing [2+2]/dearomatic [4+2] cycloadditions",
             [("rx_s2b_1b", "1b", "reactant", None), ("rx_s2b_2a", "2a", "reactant", None), ("rx_s2b_moco6", "MoCO6", "reagent", None), ("rx_s2b_xylenes", "xylenes", "solvent", None),
              ("rx_s2b_3ba", "3ba", "product", "[2+2] product"), ("rx_s2b_6ba", "6ba", "product", "dearomative [4+2] product")],
             labels=["C-P", "C-C"], unresolved=["The main-text Scheme 2 does not explicitly state the yield measurement basis for 3ba or 6ba."],
             notes=["The prose assigns 3ba to terminal-double-bond [2+2] cycloaddition and 6ba to dearomative [4+2] cycloaddition."])
add_yield(r, "s2b_y_3ba", ["3ba"], 21, [I[2]], None, "individual product 3ba")
add_yield(r, "s2b_y_6ba", ["6ba"], 36, [I[2]], None, "individual product 6ba")
reactions.append(r)

# Actual mechanistic experiments described in prose, kept separate from proposed Scheme 3.
r = reaction("mechanistic_1b_without_2a", [E[1], I[2], I[3]], "specific_experiment", ["mechanistic_study"], "PDF pages 1-3: reaction without 2a; Scheme 2 conditions and intermediate 7 in Scheme 3",
             "Reaction performed without 1,1-diphenylallene 2a", "carbonylative ring expansion to observed intermediate 7",
             [("mech1_1b", "1b", "reactant", None), ("mech1_moco6", "MoCO6", "reagent", None), ("mech1_xylenes", "xylenes", "solvent", None), ("mech1_7", "7", "product", "NMR-observed intermediate")],
             labels=(), unresolved=["Intermediate 7 is observed in crude solution and not reported as an isolated product; exact yield/conversion and atom mapping are unreported."],
             notes=["The crude reaction solution shows a new 31P NMR peak at +132 ppm, assigned by the authors to intermediate 7 with reference to literature [6]."])
r["solvent_conditions"]["composition"] = st.fact(
    "the reaction was performed without 2a; Scheme 2 prints Xylenes",
    [f"{DOC}:mol:xylenes"], [E[1], I[2]], method="inherited_general_condition")
r["steps"][0]["operation"] = st.fact(
    "the reaction was performed without 2a; Scheme 2 prints Mo(CO)6 / Xylenes, reflux",
    "repeat the 1b reaction while omitting 2a", [E[1], I[2]], method="inherited_general_condition")
r["steps"][0]["temperature"] = st.fact(
    "the reaction was performed without 2a; Scheme 2 prints reflux",
    "reflux", [E[1], I[2]], method="inherited_general_condition")
o = st.new_outcome("mech1_observed_7", "other", [f"{DOC}:mol:7"], [E[1], I[3]])
o["product_role"] = direct("reaction intermediate", "reaction intermediate", [E[1], I[3]])
o["value"] = direct("new peak at +132 ppm", {"observation": "31P NMR peak", "chemical_shift_ppm": 132, "sign": "+"}, [E[1], I[3]])
o["measurement_basis"] = direct("31P NMR of crude reaction solution", "31P NMR of crude reaction solution", [E[1], I[3]])
o["scope"] = direct("formation of intermediate 7", "intermediate 7", [E[1], I[3]])
r["outcomes"].append(o); reactions.append(r)

r = reaction("mechanistic_7_plus_2a", [E[1], I[3]], "specific_experiment", ["mechanistic_study"], "PDF pages 1 and 3: subsequent reaction of 7 with 2a",
             "Subsequent reaction of intermediate 7 with 2a", "[2+2] and dearomative [4+2] cycloadditions",
             [("mech2_7", "7", "reactant", None), ("mech2_2a", "2a", "reactant", None), ("mech2_3ba", "3ba", "product", "expected [2+2] product"), ("mech2_6ba", "6ba", "product", "expected [4+2] product")],
             conditions=False, labels=["C-P", "C-C"], unresolved=["Charge, solvent, temperature, time, yields and measurement basis for the subsequent reaction are not stated."],
             notes=["The prose reports formation of the expected 3ba and 6ba but gives no quantitative result."])
s = st.new_step("mechanistic_7_plus_2a:reaction", [E[1]])
s["participant_ids"] = ["mech2_7", "mech2_2a"]
s["operation"] = direct("Subsequent reaction of 7 with 2a", "react intermediate 7 with allene 2a", [E[1]])
r["steps"] = [s]
o = st.new_outcome("mech2_products_observed", "other", [f"{DOC}:mol:3ba", f"{DOC}:mol:6ba"], [E[1]])
o["product_role"] = direct("expected products", "expected products", [E[1]])
o["value"] = direct("led to the formation of the expected 3ba and 6ba", "formation reported; not quantified", [E[1]])
o["scope"] = direct("products 3ba and 6ba", [f"{DOC}:mol:3ba", f"{DOC}:mol:6ba"], [E[1]])
r["outcomes"].append(o); reactions.append(r)

# Proposed mechanism and Table 1 generic scheme.
r = reaction("scheme3_proposed_mechanism", [E[1], E[2], I[3]], "proposed_mechanism", ["other"], "PDF page 3, Scheme 3; mechanistic discussion spans pages 1-2",
             "Plausible reaction pathway", "proposed carbonylative ring expansion followed by alternative [2+2]/[4+2] pathways",
             [("s3_1b", "1b", "reactant", None), ("s3_moco6", "MoCO6", "reagent", None), ("s3_7", "7", "reactant", None), ("s3_2a", "2a", "reactant", None), ("s3_3ba", "3ba", "product", "proposed [2+2] branch"), ("s3_6ba", "6ba", "product", "proposed [4+2] branch")],
             conditions=False, labels=["C-P", "C-C"], unresolved=["Scheme 3 is explicitly labelled plausible; elementary steps and complete atom mapping are not experimentally established."],
             notes=["This record represents the authors' proposed pathway and is not counted as a specific experiment."])
reactions.append(r)

r = reaction("table1_generic_scheme", [I[3]], "generic_scheme", ["substrate_scope"], "PDF page 3, Table 1 header reaction",
             "Reaction of phosphirane complexes 1c-1e with allenes 2a-2e", "generic carbonylative ring expansion and [2+2] cycloaddition",
             [("t1g_1", "generic_1c_1e", "reactant", None), ("t1g_2", "generic_2a_2e", "reactant", None), ("t1g_moco6", "MoCO6", "reagent", None), ("t1g_xylenes", "xylenes", "solvent", None), ("t1g_3", "generic_3", "product", "regioisomer 3"), ("t1g_3p", "generic_3prime", "product", "alternative regioisomer 3′")],
             labels=["C-P", "C-C"], unresolved=["Generic R-group drawings are not complete molecules and cannot support complete atom mapping."],
             notes=["The generic Table 1 arrow is retained separately from the eight experimental rows."])
reactions.append(r)

# Eight Table 1 rows.
table_rows = [
    (1, "1c", "2a", ["3ca"], 48), (2, "1d", "2a", ["3da"], 55), (3, "1e", "2a", ["3ea"], 35),
    (4, "1c", "2b", ["3cb"], 61), (5, "1c", "2c", ["3cc"], 53),
    (6, "1c", "2d", ["3cd", "3cd_prime"], 43), (7, "1d", "2d", ["3dd", "3dd_prime"], 61),
]
for entry, one, two, prods, y in table_rows:
    parts = [(f"t1e{entry}_1", one, "reactant", None), (f"t1e{entry}_2", two, "reactant", None),
             (f"t1e{entry}_moco6", "MoCO6", "reagent", None), (f"t1e{entry}_xylenes", "xylenes", "solvent", None)]
    parts += [(f"t1e{entry}_{p}", p, "product", "isolated [2+2] cycloadduct" if len(prods)==1 else "isolated regioisomer mixture") for p in prods]
    r = reaction(f"table1_entry{entry}", [I[3], E[2], E[3]], "specific_experiment", ["substrate_scope"], f"PDF page 3, Table 1 entry {entry}",
                 "Reaction of phosphirane complexes 1c-1e with allenes 2a-2e", "regioselective [2+2] cycloaddition after carbonylative ring expansion",
                 parts, labels=["C-P", "C-C"], general=["table1_generic_scheme"],
                 unresolved=["Substrate/reagent amounts, reaction time, exact reflux temperature, atmosphere, pressure and workup are not reported in the main text."],
                 notes=["Table 1 footnote a states yields are of isolated product."])
    add_yield(r, f"t1e{entry}_yield", prods, y, [I[3]], "isolated", "combined regioisomer mixture" if len(prods)>1 else f"product {prods[0]}")
    if len(prods)>1:
        add_ratio(r, f"t1e{entry}_ratio", prods, "1:1", [I[3]], "31P NMR spectra", ":".join(labels[x] for x in prods))
    reactions.append(r)

# Table 1 entry 8: zero in the isolated-yield column, with no invented product.
r = reaction("table1_entry8", [I[3], E[2], E[3]], "specific_experiment", ["substrate_scope", "control"], "PDF page 3, Table 1 entry 8",
             "Reaction of phosphirane complex 1c with electron-poor allenoate 2e", "reported unsuccessful reaction",
             [("t1e8_1c", "1c", "reactant", None), ("t1e8_2e", "2e", "reactant", None), ("t1e8_moco6", "MoCO6", "reagent", None), ("t1e8_xylenes", "xylenes", "solvent", None)],
             labels=(), general=["table1_generic_scheme"], unresolved=["The table reports 0 in the isolated-product yield column but gives no intended product identity, conversion, recovery, decomposition or detection limit."],
             notes=["The prose states that the reaction did not proceed with electron-poor allenoate; no product is invented for this record."])
r["classification"] = no_classification([I[3]], "No observed transformation is reported for Table 1 entry 8.", not_applicable=True)
o = st.new_outcome("t1e8_no_reaction", "no_reaction", [], [I[3], E[2]])
o["value"] = direct("0 in Yield (%) column; reaction did not proceed", "no reaction / 0% isolated product yield reported", [I[3], E[2]])
o["measurement_basis"] = direct("Yields of the isolated product", "isolated-product yield column", [I[3]])
o["scope"] = direct("Table 1 entry 8", "reported unsuccessful reaction; conversion unreported", [I[3], E[2]])
r["outcomes"].append(o); reactions.append(r)

# Four Scheme 4 substrate-scope experiments; yields apply to regioisomer pairs.
scheme4 = [
    ("1a", "2d", ["4ad", "4ad_prime"], 53, "6:1"),
    ("1b", "2d", ["6bd", "6bd_prime"], 57, "8:1"),
    ("1b", "2f", ["6bf", "6bf_prime"], 40, "7:1"),
    ("1b", "2g", ["6bg", "6bg_prime"], 48, "3:1"),
]
for one, two, prods, y, ratio in scheme4:
    suffix = f"scheme4_{one}_{two}"
    ev = [I[3], E[2]] + ([I[4]] if prods[0] == "6bd" else [])
    parts = [(f"{suffix}_1", one, "reactant", None), (f"{suffix}_2", two, "reactant", None),
             (f"{suffix}_moco6", "MoCO6", "reagent", None), (f"{suffix}_xylenes", "xylenes", "solvent", None)]
    parts += [(f"{suffix}_{p}", p, "product", "regioisomer in purified [4+2] cycloadduct mixture") for p in prods]
    r = reaction(suffix, ev, "specific_experiment", ["substrate_scope"], f"PDF page 3, Scheme 4: {one} + {two}",
                 "Reactions of 1-iminylphosphirane complex 1a/1b with 1-arylallene 2", "dearomatic [4+2] cycloaddition after carbonylative ring expansion",
                 parts, labels=["C-C"], unresolved=["Scheme 4 and the prose do not explicitly identify the total-yield measurement basis as isolated, NMR or another method."],
                 notes=["The prose states that products were purified as a mixture of [4+2] cycloadducts and ratios were determined by 31P NMR."])
    add_yield(r, f"{suffix}_yield", prods, y, [I[3]], None, f"total yield of {labels[prods[0]]}+{labels[prods[1]]}")
    add_ratio(r, f"{suffix}_ratio", prods, ratio, [I[3], E[2]], "31P NMR", f"{labels[prods[0]]}:{labels[prods[1]]}")
    reactions.append(r)

bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:context:scheme2_31p", "document_id": DOC, "kind": "characterization",
     "evidence_ids": [I[2]], "text": "Scheme 2 prints 31P NMR shifts: 3aa 83.8 ppm; 4aa 14.4 ppm; 5aa 17.6 ppm; 3ba 78.3 ppm; 6ba 11.7 ppm. These are characterization observations, not separate syntheses."},
    {"context_id": f"{DOC}:context:xray_3aa", "document_id": DOC, "kind": "structure_assignment",
     "evidence_ids": [I[2], E[1]], "text": "The structure of 3aa was established by single-crystal X-ray crystallographic analysis (Fig. 1)."},
    {"context_id": f"{DOC}:context:xray_6bd", "document_id": DOC, "kind": "structure_assignment",
     "evidence_ids": [I[3], I[4], E[2], E[4]], "text": "The major product 6bd was confirmed by single-crystal X-ray crystallography (Fig. 2); this is not an additional synthesis record."},
    {"context_id": f"{DOC}:context:si", "document_id": DOC, "kind": "supplementary_material_statement",
     "evidence_ids": [I[4], E[4]], "text": "The article states that supplementary material is available online, but no SI file was supplied or explicitly associated in this run."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))
