#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for PDF pages 1--6.

This is a reproducible extraction script, not an alternate schema or validator.
Dense structures that were not safely transcribed remain unresolved by design.
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
DOC = "doc_d0f6434ef64e5c72"
E = {p: [f"{DOC}:page{p}:text", f"{DOC}:page{p}:image"] for p in range(1, 7)}


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, derivation=None):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation)


def derived(normalized, evidence, derivation, raw=None):
    return f(raw, normalized, evidence, method="structure_derivation",
             derivation=derivation)


def not_applicable():
    return f(method="not_applicable", status="not_applicable")


template = json.loads((RUN / "submitted.json").read_text())
bundle = copy.deepcopy(template)
bundle["molecules"] = []
bundle["reactions"] = []
bundle["source_context"] = []
M = {}
diagnostics = {}


def add_molecule(label, evidence, smiles=None, name=None, stereo_kind="unresolved",
                 identity_note=None):
    mid = f"{DOC}:mol:{label}"
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    m["cas_numbers"] = f()
    m["identity_note"] = f(identity_note, identity_note, evidence) if identity_note else f()
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(smiles, evidence,
            "Machine candidate transcribed from the source drawing/name; RDKit parsing is not source-structure verification.",
            raw={"source": "drawing or explicit chemical name", "label": label})
        diag = st.stereo_inspect(smiles)
        diagnostics[mid] = diag
        m["stereochemistry"]["isomeric_smiles"] = derived(
            diag["canonical_isomeric_smiles"], evidence,
            "RDKit canonical isomeric SMILES from the machine candidate; source_verified=false.")
        m["stereochemistry"]["identity_key"] = derived(
            diag["stereo_graph_key"], evidence,
            "RDKit stereo graph key from the machine candidate; not a sample identity.")
        if stereo_kind == "achiral":
            m["stereochemistry"]["chirality_status"] = derived(
                "achiral", evidence,
                "The candidate graph and source drawing contain no stereogenic center; terminal allene substrates have two identical hydrogens at one terminus.")
            m["stereochemistry"]["geometry_status"] = derived(
                "not_applicable", evidence,
                "No independently assignable E/Z element is present in the candidate graph.")
            m["stereochemistry"]["representation_status"] = derived(
                "not_applicable", evidence,
                "No stereochemical descriptor is required for this candidate graph.")
        else:
            m["stereochemistry"]["chirality_status"] = derived(
                "unresolved", evidence, "Machine graph diagnostics do not establish source stereochemistry.")
            m["stereochemistry"]["geometry_status"] = derived(
                "unresolved", evidence, "Machine graph diagnostics do not establish source geometry.")
            m["stereochemistry"]["representation_status"] = derived(
                "unresolved", evidence, "Stereo encoding has not been established from the source.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = f()
        if stereo_kind == "substrate_terminal_allene":
            m["stereochemistry"]["chirality_status"] = derived(
                "achiral", evidence,
                "The source drawing shows a terminal allene with two identical terminal hydrogens and no other stereogenic element, although the full connectivity was not safely encoded.")
            m["stereochemistry"]["geometry_status"] = derived(
                "not_applicable", evidence,
                "The terminal allene and depicted ring alkene do not carry an assigned E/Z descriptor in the source.")
            m["stereochemistry"]["representation_status"] = derived(
                "not_applicable", evidence, "No source stereodescriptor is present for the substrate.")
        elif stereo_kind == "racemic_product":
            m["stereochemistry"]["atom_centers"] = f(
                "two adjacent all-carbon quaternary stereogenic centers; relative wedges/dashes as drawn",
                [{"description": "two adjacent all-carbon quaternary stereogenic centers",
                  "configuration": "relative configuration as drawn; exact R/S unresolved"}], evidence)
            m["stereochemistry"]["mixture"] = f(
                f"rac-{label}", {"type": "racemate", "ratio": "not reported"}, evidence)
            m["stereochemistry"]["assignment"] = f(
                "racemic product; complete diastereoselectivity reported for the isolable scope series",
                {"absolute_configuration": "unassigned", "relative_configuration": "as drawn"}, evidence)
            m["stereochemistry"]["chirality_status"] = f("rac-", "racemic", evidence)
            m["stereochemistry"]["geometry_status"] = derived(
                "unresolved", evidence, "Dense polycyclic drawing was not fully encoded; alkene geometry was not separately assigned.")
            m["stereochemistry"]["representation_status"] = derived(
                "not_encoded", evidence, "Wedge/dash relative stereochemistry is retained as text because no complete SMILES was asserted.")
        else:
            for key in ("chirality_status", "geometry_status", "representation_status"):
                m["stereochemistry"][key] = derived(
                    "unresolved", evidence, "Source identity is retained, but the full structure/stereochemistry was not safely encoded.")
    bundle["molecules"].append(m)
    M[label] = mid
    return mid


def add_generic(label, evidence, fixed, variables, note):
    mid = add_molecule(label, evidence, identity_note=note)
    m = bundle["molecules"][-1]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, fixed, evidence),
        "variables": [{
            "label": v,
            "attachment_points": f(f"source drawing attachment for {v}", f"source drawing attachment for {v}", evidence),
            "definition": f(desc, desc, evidence),
            "allowed_values": f(),
            "constraints": f(),
        } for v, desc in variables],
        "instance_mapping": [],
        "reason": f(note, note, evidence),
    }
    return mid


# Generic current-work structures and the legible terminal-allene substrate candidates.
add_generic("generic_1", E[1], "benzene-tethered yne–allenone framework shown in Scheme 1c/2",
            [("R1", "substituent on the arylalkynyl/cyclohexenyl group"),
             ("R2", "substituent on the internal arene ring")],
            "R-group source drawing; not represented as wildcard SMILES.")
add_generic("generic_rac_2", E[1], "racemic hexacyclic 1-naphthol framework shown in Scheme 1c/2",
            [("R1", "two repeated substituents inherited from the dimerized substrate"),
             ("R2", "two repeated internal-arene substituents")],
            "Dense racemic R-group product drawing; relative stereochemistry retained textually.")

substrate_smiles = {
    "1a": "C=C=CC(=O)c1ccccc1C#Cc1ccccc1",
    "1b": "C=C=CC(=O)c1ccccc1C#Cc1ccc(C)cc1",
    "1c": "C=C=CC(=O)c1ccccc1C#Cc1cccc(C)c1",
    "1d": "C=C=CC(=O)c1ccccc1C#Cc1ccc(CC)cc1",
    "1e": "C=C=CC(=O)c1ccccc1C#Cc1ccc(C(C)(C)C)cc1",
    "1f": "C=C=CC(=O)c1ccccc1C#Cc1ccc(OC)cc1",
    "1g": "C=C=CC(=O)c1ccccc1C#Cc1ccc(F)cc1",
    "1h": "C=C=CC(=O)c1ccccc1C#Cc1ccc(Cl)cc1",
    "1i": "C=C=CC(=O)c1ccccc1C#Cc1ccc(Br)cc1",
    "1j": "C=C=CC(=O)c1ccccc1C#Cc1cccs1",
    "1x": "C=C=CC(=O)c1ccccc1C#CC1=CCCCC1",
}
for label in [f"1{x}" for x in "abcdefghij"] + ["1x"]:
    add_molecule(label, E[2], substrate_smiles[label], stereo_kind="achiral",
                 identity_note="Machine SMILES candidate from the source drawing; positional and connectivity review remains required.")
for label in [f"1{x}" for x in "klmnopqrstuvw"]:
    add_molecule(label, E[2], stereo_kind="substrate_terminal_allene",
                 identity_note="Internal-arene substitution position is retained by source label/drawing; full SMILES was not asserted.")

# Products 2a--2x and reductive product 3.
for label in [f"2{chr(c)}" for c in range(ord("a"), ord("x") + 1)]:
    ev = E[2] + (E[3] if label in {"2a", "2q"} else [])
    add_molecule(label, ev, stereo_kind="racemic_product",
                 identity_note="Dense polycyclic structure retained by source label and drawing; complete isomeric SMILES not asserted.")
add_molecule("3", E[3], stereo_kind="racemic_product",
             identity_note="Racemic reductive product from Scheme 3b; dense structure not fully encoded.")

# Explicit reagents/solvents. Structures are candidates only; names and doses remain source-linked.
aux = {
    "MeCN": ("acetonitrile", "CC#N", E[1] + E[2] + E[3]),
    "1,4-dioxane": ("1,4-dioxane", "O1CCOCC1", E[2]),
    "acetone": ("acetone", "CC(=O)C", E[2]),
    "THF": ("tetrahydrofuran", "C1CCOC1", E[2]),
    "DMF": ("N,N-dimethylformamide", "CN(C)C=O", E[2]),
    "toluene": ("toluene", "Cc1ccccc1", E[2]),
    "DCE": ("1,2-dichloroethane", "ClCCCl", E[2] + E[3]),
    "BHT": ("butylhydroxytoluene (BHT)", "CC(C)(C)c1cc(C)c(O)c(C)c1C(C)(C)C", E[3]),
    "triplet_quencher": ("2,5-dimethylhexa-2,4-diene", "CC(=CC=C(C)C)C", E[3]),
    "silane": ("dimethoxy(methyl)silane, Me(MeO)2SiH", "CO[SiH](C)OC", E[3]),
    "NaOMe": ("sodium methoxide", "C[O-].[Na+]", E[3]),
}
for label, (name, smiles, ev) in aux.items():
    add_molecule(label, list(dict.fromkeys(ev)), smiles, name=name, stereo_kind="achiral")
for label, name in {
    "TEMPO": "2,2,6,6-tetramethyl-1-piperidinyloxy (TEMPO)",
    "CuCN": "copper(I) cyanide (CuCN)",
    "DPPF": "1,1'-bis(diphenylphosphino)ferrocene (DPPF)",
}.items():
    add_molecule(label, E[3], name=name, identity_note="Named reagent retained; no complete source-verified SMILES asserted.")
add_generic("TEMPO_or_BHT", E[3], "one of two alternative radical scavengers, as printed 'TEMPO or BHT'",
            [("identity", "TEMPO or BHT; source does not assign the single printed 62% separately")],
            "Source-disjunctive additive identity; not represented as simultaneous addition of both scavengers.")

# Literature-example and proposed-mechanism identities.
for label in ["Yoon_substrate", "Yoon_product", "Kwon_substrate", "Kwon_product"]:
    add_molecule(label, E[1], identity_note="Prior-work structure in Scheme 1a/1b; retained as a literature example, not a current experiment.")
for label in ["A", "A_star", "B", "C", "D", "E", "F"]:
    add_molecule(label, E[4], identity_note="Author-proposed intermediate in Scheme 5; not an isolated or directly measured compound here.")


def participant(pid, mol, role, evidence, stage="reaction", amount=None,
                equiv=None, molpct=None, volume=None, ref=None, product_role=None):
    p = st.new_participant(pid, M[mol], role, evidence, stage=stage)
    if product_role:
        p["product_role"] = f(product_role, product_role, evidence)
    for key, raw, norm in [
        ("amount", None if amount is None else f"{amount[0]} {amount[1]}", None if amount is None else {"number": amount[0], "unit": amount[1]}),
        ("equivalents", None if equiv is None else f"{equiv} equiv", None if equiv is None else {"number": equiv, "unit": "equiv"}),
        ("mol_percent", None if molpct is None else f"{molpct} mol%", None if molpct is None else {"number": molpct, "unit": "mol%"}),
        ("volume", None if volume is None else f"{volume[0]} {volume[1]}", None if volume is None else {"number": volume[0], "unit": volume[1]}),
    ]:
        if norm is not None:
            p["charges"][key] = f(raw, norm, evidence)
    if ref:
        p["equivalence_reference"] = f(ref, M[ref], evidence)
    return p


def classification(kind, evidence):
    if kind == "none":
        return {
            "standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [],
            "labels": derived([], evidence, "No observed transformation; no expected-product classification was assigned."),
            "coverage": derived("not_applicable", evidence, "No reaction was observed."),
            "limitations": ["Negative control; expected-product bond formation is not inferred."], "model_votes": []}
    if kind in {"C-C", "C-H"}:
        partner = "C" if kind == "C-C" else "H"
        change_type = "formed" if kind == "C-C" else "hydrogen_gain"
        change = {"change_id": "b1", "atoms": [{"atom_map": None, "element": "C"},
                  {"atom_map": None, "element": partner}], "before_order": 0,
                  "after_order": 1, "change_type": change_type}
        label = {"category_number": 2 if kind == "C-C" else 1,
                 "code": kind, "supporting_change_ids": ["b1"]}
        return {
            "standard_id": st.STANDARD, "mapped_reaction_smiles": f(),
            "bond_changes": [derived(change, evidence,
                "At least one net bond of this element pair is evident from the depicted transformation; exact atom mapping and total bond-change inventory remain unresolved.")],
            "labels": derived([label], evidence,
                "Partial taxonomy assignment from the source transformation; not from a complete atom mapping."),
            "coverage": derived("partial", evidence,
                "Dense structures prevent a complete mapped bond-change inventory."),
            "limitations": ["Exact atom correspondence and the number of all new/broken/order-changed bonds are unresolved.",
                            "Implicit hydrogen transfers and proton source are not completely mapped."],
            "model_votes": []}
    return {
        "standard_id": st.STANDARD, "mapped_reaction_smiles": f(), "bond_changes": [],
        "labels": f(),
        "coverage": derived("unresolved", evidence, "No reliable complete structure mapping was made for this source arrow."),
        "limitations": ["Literature-example or unresolved dense structure; no bond label asserted."], "model_votes": []}


def base_reaction(suffix, evidence, record_type, categories, location, name,
                  transformation, class_kind="unresolved"):
    rid = f"{DOC}:{suffix}"
    r = st.new_reaction(rid, DOC, evidence)
    r["record_type"] = record_type
    if categories is None:
        r["experiment_categories"] = not_applicable()
    else:
        r["experiment_categories"] = f(categories, categories, evidence)
    r["location"] = f(location, location, evidence)
    r["reaction_name"] = f(name, name, evidence)
    r["transformation_type"] = derived(transformation, evidence,
        "Normalized transformation description based on the source arrow/text.")
    r["classification"] = classification(class_kind, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    return r


def set_photochemical_conditions(r, evidence, solvent="MeCN", temperature=30,
                                 time=15, atmosphere="argon", light="blue LEDs",
                                 power=30, amount=0.2, volume=2.0,
                                 inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    reactant_pid = next(p["participant_id"] for p in r["participants"] if p["role"] == "reactant")
    solvent_pid = f"{r['reaction_id']}:solvent"
    r["participants"].append(participant(solvent_pid, solvent, "solvent", evidence,
                                         volume=(volume, "mL") if volume is not None else None))
    r["reaction_scale"] = f(
        f"{amount} mmol substrate", {"reference_molecule_id": r["participants"][0]["molecule_id"],
        "amount": {"number": amount, "unit": "mmol"}}, evidence, method=method)
    r["solvent_conditions"]["composition"] = f(f"dry {solvent}", [f"dry {solvent}"], evidence, method=method)
    if volume is not None:
        r["solvent_conditions"]["total_volume"] = f(f"{volume} mL", {"number": volume, "unit": "mL"}, evidence, method=method)
    s = st.new_step(f"{r['reaction_id']}:reaction", evidence)
    s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f("irradiation under the stated conditions", "irradiation", evidence, method=method)
    s["temperature"] = f(f"{temperature} °C", {"number": temperature, "unit": "°C"}, evidence, method=method)
    s["time"] = f(f"{time} h", {"number": time, "unit": "h"}, evidence, method=method)
    s["atmosphere"] = f(atmosphere, atmosphere, evidence, method=method)
    r["steps"] = [s]
    if light is None:
        r["special_conditions"]["light_source"] = f("without blue LEDs", "no light irradiation", evidence)
    else:
        r["special_conditions"]["light_source"] = f(light, light, evidence, method=method)
    if power is not None:
        r["special_conditions"]["power"] = f(f"{power} W", {"number": power, "unit": "W"}, evidence, method=method)
    return reactant_pid


def add_yield(r, mol, value, evidence, isolated=False, outcome_id="yield"):
    o = st.new_outcome(f"{r['reaction_id']}:{outcome_id}", "yield", [M[mol]], evidence)
    o["product_role"] = f("main product", "main product", evidence)
    o["value"] = f(f"{value}%", {"number": value, "unit": "%"}, evidence)
    if isolated:
        o["yield_type"] = f("isolated yield", "isolated", evidence)
        o["measurement_basis"] = f("Isolated yield based on substrate 1a.", "isolated; based on substrate 1a", evidence)
    r["outcomes"].append(o)


def add_selectivity(r, mol, raw, evidence):
    o = st.new_outcome(f"{r['reaction_id']}:diastereoselectivity", "selectivity", [M[mol]], evidence)
    o["product_role"] = f("main product", "main product", evidence)
    o["value"] = f(raw, raw, evidence)
    o["scope"] = f("diastereoselectivity; absolute configurations not assigned", {
        "metric": "diastereoselectivity", "major_isomer": "relative configuration as drawn",
        "minor_isomer": "not reported", "absolute_configuration": "unassigned"}, evidence)
    r["outcomes"].append(o)


# Scheme 1: two prior-work arrows and one current generic transformation.
for suffix, sub, prod, light, title in [
    ("scheme1a_yoon_literature", "Yoon_substrate", "Yoon_product", "visible light", "Yoon's work"),
    ("scheme1b_kwon_literature", "Kwon_substrate", "Kwon_product", "white LED", "Kwon's work"),
]:
    r = base_reaction(suffix, E[1], "literature_example", None, f"Scheme 1 {title}", title,
                      "prior-work [2 + 2] photocycloaddition", "unresolved")
    r["participants"] = [participant(f"{r['reaction_id']}:substrate", sub, "reactant", E[1]),
                         participant(f"{r['reaction_id']}:product", prod, "product", E[1], product_role="depicted product")]
    s = st.new_step(f"{r['reaction_id']}:reaction", E[1]); s["participant_ids"] = [p["participant_id"] for p in r["participants"]]
    s["operation"] = f(light, "photochemical irradiation", E[1]); r["steps"] = [s]
    r["special_conditions"]["light_source"] = f(light, light, E[1])
    r["review"]["unresolved"] = ["Prior-work example: exact identities, quantities, time, temperature and yield are outside this paper's current experiments."]
    bundle["reactions"].append(r)

general_id = f"{DOC}:scheme1c_generic"
r = base_reaction("scheme1c_generic", E[1], "generic_scheme", None,
                  "Scheme 1c (this work)", "photoinduced double [2 + 2] cycloaddition relay",
                  "photocatalyst-free cyclodimerization of yne–allenones", "C-C")
r["participants"] = [participant(f"{r['reaction_id']}:substrate", "generic_1", "reactant", E[1]),
                     participant(f"{r['reaction_id']}:product", "generic_rac_2", "product", E[1], product_role="racemic hexacyclic product")]
set_photochemical_conditions(r, E[1], amount=0.2, volume=None)
r["reaction_scale"] = f()  # Scheme 1c does not print a scale.
r["review"]["unresolved"] = ["Generic R1/R2 structures are templates; no experiment-specific yield or full atom mapping is printed in Scheme 1c."]
bundle["reactions"].append(r)


# Table 1 optimization rows. Entry 1 is also the 2a scope point in Scheme 2.
opt = {
    1: ("blue", 65), 2: ("white", 45), 3: ("green", "trace"),
    4: ("1,4-dioxane", 50), 5: ("acetone", 55), 6: ("THF", "trace"),
    7: ("DMF", "trace"), 8: ("toluene", "trace"), 9: ("DCE", "trace"),
    10: ("40C", 60), 11: ("25C", 30), 12: ("dark", "NR"), 13: ("air", 10),
}
optimization_ids = []
for entry, (variation, result) in opt.items():
    rid_suffix = f"table1_entry_{entry:02d}"
    categories = ["optimization"] + (["substrate_scope"] if entry == 1 else [])
    loc = f"Table 1 entry {entry}" + ("; Scheme 2 product 2a is the same 65% standard-condition experiment" if entry == 1 else "")
    rx = base_reaction(rid_suffix, E[2], "specific_experiment", categories, loc,
                       "photoinduced double [2 + 2] cycloaddition relay",
                       "cyclodimerization of 1a to rac-2a", "none" if entry == 12 else "C-C")
    rx["general_scheme_ids"] = [general_id]
    rx["participants"] = [participant(f"{rx['reaction_id']}:1a", "1a", "reactant", E[2], amount=(0.2, "mmol"), equiv=1.0, ref="1a")]
    if entry != 12:
        rx["participants"].append(participant(f"{rx['reaction_id']}:2a", "2a", "product", E[2], product_role="main product"))
    solvent = variation if variation in {"1,4-dioxane", "acetone", "THF", "DMF", "toluene", "DCE"} else "MeCN"
    temp = 40 if variation == "40C" else 25 if variation == "25C" else 30
    atmosphere = "air" if variation == "air" else "argon"
    light = None if variation == "dark" else "white LEDs" if variation == "white" else "green LEDs" if variation == "green" else "blue LEDs"
    set_photochemical_conditions(rx, E[2], solvent=solvent, temperature=temp, atmosphere=atmosphere,
                                 light=light, power=None if variation == "dark" else 30, inherited=True)
    rx["procedure"]["original_text"] = f(
        "Reaction conditions: benzene-tethered yne–allenone (1a, 0.2 mmol, 1.0 equiv.), dry MeCN (2 mL), 30 °C, under argon conditions, for 15 hours; row-specific variation as printed.",
        "Table 1 footnote a plus row-specific variation", E[2])
    if isinstance(result, int):
        add_yield(rx, "2a", result, E[2], isolated=True)
        if entry == 1:
            add_selectivity(rx, "2a", "complete diastereoselectivity", E[1] + E[2])
    elif result == "trace":
        o = st.new_outcome(f"{rx['reaction_id']}:trace", "trace", [M["2a"]], E[2])
        o["product_role"] = f("expected product", "expected product", E[2]); o["value"] = f("Trace", "trace", E[2]); rx["outcomes"].append(o)
    else:
        o = st.new_outcome(f"{rx['reaction_id']}:no_reaction", "no_reaction", [], E[2])
        o["value"] = f("NR", "no reaction", E[2]); rx["outcomes"].append(o)
        rec = st.new_outcome(f"{rx['reaction_id']}:recovery", "recovery", [M["1a"]], E[2])
        rec["value"] = f("starting material 1a was recovered", "starting material recovered; amount unreported", E[2])
        rx["outcomes"].append(rec)
    rx["review"]["unresolved"] = ["Work-up and purification are not reported on covered pages; SI was not supplied."]
    if result in {"trace", "NR"}:
        rx["review"]["unresolved"].append("Detection method and detection limit for the qualitative result are unreported.")
    bundle["reactions"].append(rx); optimization_ids.append(rx["reaction_id"])


# Scheme 2 scope. 2a is already represented by Table 1 entry 1.
scope_yields = {
    "b": 78, "c": 65, "d": 80, "e": 65, "f": 55, "g": 52, "h": 45, "i": 56, "j": 70,
    "k": 45, "l": 65, "m": 60, "n": 70, "o": 72, "p": 80, "q": 61,
    "r": 75, "s": 72, "t": 55, "u": 71, "v": 82, "w": 65, "x": "trace",
}
scope_ids = []
for letter, result in scope_yields.items():
    sub, prod = f"1{letter}", f"2{letter}"
    rx = base_reaction(f"scheme2_{prod}", E[2], "specific_experiment", ["substrate_scope"],
                       f"Scheme 2, {sub} to {prod}", "photoinduced double [2 + 2] cycloaddition relay",
                       f"cyclodimerization of {sub} to rac-{prod}", "C-C")
    rx["general_scheme_ids"] = [general_id]
    rx["participants"] = [participant(f"{rx['reaction_id']}:{sub}", sub, "reactant", E[2], amount=(0.2, "mmol"), equiv=1.0, ref=sub),
                          participant(f"{rx['reaction_id']}:{prod}", prod, "product", E[2], product_role="main product")]
    set_photochemical_conditions(rx, E[2], inherited=True)
    rx["procedure"]["general_procedure_reference"] = f("standard conditions", "Table 1 entry 1 standard conditions", E[2])
    if isinstance(result, int):
        add_yield(rx, prod, result, E[2], isolated=False)
        add_selectivity(rx, prod, "single diastereoisomer" if letter == "j" else "complete diastereoselectivity", E[2] + (E[3] if letter == "q" else []))
    else:
        o = st.new_outcome(f"{rx['reaction_id']}:trace", "trace", [M[prod]], E[2])
        o["product_role"] = f("expected product", "expected product", E[2]); o["value"] = f("trace", "trace", E[2]); rx["outcomes"].append(o)
        rx["notes"].append(f("product 2x was detected in a trace amount and could not be isolated", "trace; not isolated", E[2] + E[3]))
    rx["review"]["unresolved"] = [
        "Scheme 2 does not state the yield measurement basis; isolated/NMR status is not assumed.",
        "Work-up, purification and experiment-specific procedure are unavailable because SI was not supplied.",
        "Complete atom mapping and exact stereocenter descriptors are unresolved."]
    if letter == "x": rx["review"]["unresolved"].append("Detection method and detection limit for trace 2x are unreported.")
    bundle["reactions"].append(rx); scope_ids.append(rx["reaction_id"])


# Scheme 3a scale-up.
rx = base_reaction("scheme3a_scaleup_2q", E[3], "specific_experiment", ["scale_up"],
                   "Scheme 3a", "scale-up transformation of 2q",
                   "cyclodimerization of 1q to rac-2q", "C-C")
rx["general_scheme_ids"] = [general_id]
rx["participants"] = [participant(f"{rx['reaction_id']}:1q", "1q", "reactant", E[3], amount=(2.0, "mmol"), equiv=1.0, ref="1q"),
                      participant(f"{rx['reaction_id']}:2q", "2q", "product", E[3], product_role="main product")]
set_photochemical_conditions(rx, E[3], amount=2.0, volume=None)
add_yield(rx, "2q", 54, E[3], isolated=False)
o = st.new_outcome(f"{rx['reaction_id']}:product_amount", "product_amount", [M["2q"]], E[3]); o["product_role"] = f("main product", "main product", E[3]); o["value"] = f("0.63 g", {"number": 0.63, "unit": "g"}, E[3]); rx["outcomes"].append(o)
rx["review"]["unresolved"] = ["MeCN volume/concentration, yield measurement basis, work-up and purification are unreported on covered pages.", "Full structures and atom mapping remain unresolved."]
bundle["reactions"].append(rx); scaleup_id = rx["reaction_id"]


# Scheme 3b derivatization.
rx = base_reaction("scheme3b_reduction_3", E[3], "specific_experiment", ["derivatization"],
                   "Scheme 3b", "CuH-catalyzed reduction reaction of 2q",
                   "reduction of conjugated cycloketone 2q to rac-3", "C-H")
rx["participants"] = [
    participant(f"{rx['reaction_id']}:2q", "2q", "reactant", E[3]),
    participant(f"{rx['reaction_id']}:silane", "silane", "reagent", E[3], equiv=1.5, ref="2q"),
    participant(f"{rx['reaction_id']}:CuCN", "CuCN", "catalyst", E[3], molpct=10, ref="2q"),
    participant(f"{rx['reaction_id']}:DPPF", "DPPF", "ligand", E[3], molpct=10, ref="2q"),
    participant(f"{rx['reaction_id']}:NaOMe", "NaOMe", "additive", E[3], equiv=1.5, ref="2q"),
    participant(f"{rx['reaction_id']}:DCE", "DCE", "solvent", E[3]),
    participant(f"{rx['reaction_id']}:3", "3", "product", E[3], product_role="main product"),
]
rx["solvent_conditions"]["composition"] = f("DCE", ["DCE"], E[3])
s = st.new_step(f"{rx['reaction_id']}:reaction", E[3]); s["participant_ids"] = [p["participant_id"] for p in rx["participants"]]
s["operation"] = f("CuH-catalyzed reduction", "reduction", E[3]); s["temperature"] = f("100 °C", {"number": 100, "unit": "°C"}, E[3]); s["atmosphere"] = f("Ar", "argon", E[3]); rx["steps"] = [s]
add_yield(rx, "3", 68, E[3], isolated=False)
rx["review"]["unresolved"] = ["Reaction time, substrate charge, DCE volume/concentration, yield basis, work-up and purification are unreported on covered pages.", "The exact reduced bond and full implicit-H mapping are unresolved."]
bundle["reactions"].append(rx); deriv_id = rx["reaction_id"]


# Scheme 4 mechanistic controls.
rx = base_reaction("scheme4a_scavenger_control", E[3], "specific_experiment", ["control", "mechanistic_study"],
                   "Scheme 4a", "radical-scavenger control", "1a cyclodimerization with TEMPO or BHT", "C-C")
rx["general_scheme_ids"] = [general_id]
rx["participants"] = [participant(f"{rx['reaction_id']}:1a", "1a", "reactant", E[3], amount=(0.2, "mmol"), equiv=1.0, ref="1a"),
                      participant(f"{rx['reaction_id']}:scavenger", "TEMPO_or_BHT", "additive", E[3], equiv=1.5, ref="1a"),
                      participant(f"{rx['reaction_id']}:2a", "2a", "product", E[3], product_role="main product")]
set_photochemical_conditions(rx, E[3], inherited=True)
add_yield(rx, "2a", 62, E[3], isolated=False)
rx["notes"].append(f("TEMPO or BHT", "source-disjunctive radical scavenger condition", E[3]))
rx["review"]["unresolved"] = ["The source prints one 62% outcome for 'TEMPO or BHT' and does not separately assign outcomes to the two scavengers; this record is intentionally not split.", "Yield measurement basis, work-up and purification are unreported."]
bundle["reactions"].append(rx); control_a_id = rx["reaction_id"]

rx = base_reaction("scheme4b_triplet_quencher", E[3], "specific_experiment", ["control", "mechanistic_study"],
                   "Scheme 4b", "triplet-quencher control", "1a cyclodimerization with 2,5-dimethylhexa-2,4-diene", "C-C")
rx["general_scheme_ids"] = [general_id]
rx["participants"] = [participant(f"{rx['reaction_id']}:1a", "1a", "reactant", E[3], amount=(0.2, "mmol"), equiv=1.0, ref="1a"),
                      participant(f"{rx['reaction_id']}:quencher", "triplet_quencher", "additive", E[3]),
                      participant(f"{rx['reaction_id']}:2a", "2a", "product", E[3], product_role="expected product")]
set_photochemical_conditions(rx, E[3], inherited=True)
o = st.new_outcome(f"{rx['reaction_id']}:trace", "trace", [M["2a"]], E[3]); o["product_role"] = f("expected product", "expected product", E[3]); o["value"] = f("trace", "trace", E[3]); rx["outcomes"].append(o)
rx["review"]["unresolved"] = ["Triplet-quencher loading, trace detection method/limit, work-up and purification are unreported."]
bundle["reactions"].append(rx); control_b_id = rx["reaction_id"]


# Scheme 5 author-proposed mechanism, not an experimental record.
rx = base_reaction("scheme5_proposed_mechanism", E[4], "proposed_mechanism", None,
                   "Scheme 5 and page 4 mechanism discussion", "plausible reaction pathway",
                   "proposed double [2 + 2] cycloaddition relay, ring opening, isomerization and proton transfer", "C-C")
# Build participants without treating intermediates as isolated products.
rx["participants"] = [participant(f"{rx['reaction_id']}:1a", "1a", "reactant", E[4])]
for label in ["A", "A_star", "B", "C", "D", "E", "F"]:
    rx["participants"].append(participant(f"{rx['reaction_id']}:{label}", label, "unknown", E[4]))
rx["participants"].append(participant(f"{rx['reaction_id']}:2a", "2a", "product", E[4], product_role="proposed final product"))
s = st.new_step(f"{rx['reaction_id']}:proposal", E[4]); s["participant_ids"] = [p["participant_id"] for p in rx["participants"]]
s["operation"] = f("putative pathway under visible-light/blue-LED irradiation", "author-proposed mechanism", E[4]); rx["steps"] = [s]
rx["special_conditions"]["light_source"] = f("blue LEDs / visible light irradiation", "visible-light irradiation", E[4])
rx["notes"].append(f("Although the detailed mechanism for this cyclodimerization is not clear to us yet", "authors state detailed mechanism remains unclear", E[4]))
rx["review"]["unresolved"] = ["Scheme 5 is a proposed mechanism, not direct experimental proof; intermediate identities, atom mapping and elementary hydrogen source are unresolved."]
bundle["reactions"].append(rx); mechanism_id = rx["reaction_id"]


# Non-synthetic observations and source controls.
bundle["source_context"] = [
    {"context_id": f"{DOC}:metadata", "kind": "bibliographic_metadata", "evidence_ids": E[1],
     "title": "Photoinduced double [2 + 2] cycloaddition relay of yne–allenones for highly diastereoselective synthesis of hexacyclic 1-naphthols",
     "doi": "10.1039/d0qo00917b"},
    {"context_id": f"{DOC}:esi_notice", "kind": "supplementary_information_notice", "evidence_ids": E[1],
     "observation": "The article states ESI is available, but no SI file was supplied or explicitly associated."},
    {"context_id": f"{DOC}:figure1_xray_2q", "kind": "characterization", "evidence_ids": E[1] + E[3],
     "target_molecule_id": M["2q"], "observation": "The stereo-structure of 2q was identified by single-crystal X-ray diffraction; CCDC 2015069. The article labels the product rac-2q."},
    {"context_id": f"{DOC}:figure2_emission", "kind": "photophysical_measurement", "evidence_ids": E[4],
     "target_molecule_id": M["1a"], "observation": "Emission spectrum of 1a (5 × 10^-5 M in MeCN), excitation at 427 nm; emission maximum shown near 500 nm."},
    {"context_id": f"{DOC}:figure3_absorption", "kind": "photophysical_measurement", "evidence_ids": E[4],
     "target_molecule_id": M["1a"], "observation": "Absorption spectrum of 1a (5 × 10^-5 M in MeCN); source text states absorption in the 400–450 nm range and figure marks 420 nm."},
]


# Exact page coverage. Page 2 carries all optimization/scope records; 2a is not duplicated.
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 1, "state": "attempted", "regions": [
            {"region": "title/abstract and ESI footnote", "finding": "Current article metadata and statement that ESI exists; no SI file was supplied or associated.", "evidence_ids": E[1]},
            {"region": "Scheme 1a and 1b", "finding": "Two prior-work photocycloaddition arrows recorded as literature examples, not current experiments.", "evidence_ids": E[1]},
            {"region": "Scheme 1c", "finding": "Current generic yne–allenone cyclodimerization under dry MeCN, Ar, 30 °C, 15 h and 30 W blue LEDs.", "evidence_ids": E[1]}],
         "record_ids": [f"{DOC}:scheme1a_yoon_literature", f"{DOC}:scheme1b_kwon_literature", general_id],
         "unresolved": ["ESI unavailable; generic R1/R2 structures and exact mapped bond changes remain unresolved."]},
        {"page": 2, "state": "attempted", "regions": [
            {"region": "Table 1 including footnotes", "finding": "All 13 optimization rows separately recorded; entry 12 is no reaction with recovered 1a.", "evidence_ids": E[2]},
            {"region": "Scheme 2", "finding": "Products 2a–2x inventoried; 2a is the same experiment as Table 1 entry 1, and 2b–2x give 23 additional scope records.", "evidence_ids": E[2]}],
         "record_ids": optimization_ids + scope_ids,
         "unresolved": ["Scope yield measurement basis and experiment-specific work-up/purification are not stated on this page.", "Dense product SMILES and exact stereocenter descriptors are unresolved."]},
        {"page": 3, "state": "attempted", "regions": [
            {"region": "Figure 1", "finding": "X-ray characterization of rac-2q (CCDC 2015069); not a separate synthesis.", "evidence_ids": E[3]},
            {"region": "Scheme 3a", "finding": "2 mmol scale-up of 1q to rac-2q, 0.63 g, 54%.", "evidence_ids": E[3]},
            {"region": "Scheme 3b", "finding": "CuH-type reductive derivatization of 2q to rac-3, 68%, with all printed reagent loadings retained.", "evidence_ids": E[3]},
            {"region": "Scheme 4a and 4b", "finding": "Radical-scavenger and triplet-quencher mechanistic controls recorded; TEMPO/BHT outcome not split because source gives one 62% value.", "evidence_ids": E[3]}],
         "record_ids": [scaleup_id, deriv_id, control_a_id, control_b_id],
         "unresolved": ["Scheme 3b reaction time and yield basis are unreported.", "Scheme 4a does not assign 62% separately to TEMPO and BHT.", "Triplet-quencher loading/detection basis is unreported."]},
        {"page": 4, "state": "attempted", "regions": [
            {"region": "Figures 2 and 3", "finding": "Emission and absorption spectra are non-synthetic photophysical observations retained in source_context.", "evidence_ids": E[4]},
            {"region": "Scheme 5 and mechanism prose", "finding": "Author-proposed pathway recorded as proposed_mechanism, not as an experimentally proven reaction sequence.", "evidence_ids": E[4]},
            {"region": "conclusions", "finding": "Summary repeats reaction claims but adds no independent experiment.", "evidence_ids": E[4]}],
         "record_ids": [mechanism_id],
         "unresolved": ["Authors explicitly state the detailed mechanism is unclear; intermediates and exact atom mapping remain unresolved."]},
        {"page": 5, "state": "attempted", "regions": [
            {"region": "references 1–7", "finding": "Reference list only; no current-paper experiment or independent characterization table.", "evidence_ids": E[5]}],
         "record_ids": [], "unresolved": []},
        {"page": 6, "state": "attempted", "regions": [
            {"region": "references 7–10", "finding": "Reference list only; no current-paper experiment or independent characterization table.", "evidence_ids": E[6]}],
         "record_ids": [], "unresolved": []},
    ]}

write(RUN / "raw_model/stereo_diagnostics.json", diagnostics)
write(RUN / "raw_model/extracted_bundle_pages_0001_0006.json", bundle)
write(RUN / "submitted.json", bundle)
write(RUN / "coverage/pages_0001_0006.json", coverage)

# Mark the prepared page tasks attempted without claiming semantic completion.
tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task["document_id"] == DOC and 1 <= int(task["evidence_id"].split(":page")[1].split(":")[0]) <= 6:
        task["state"] = "attempted"
        task["reason"] = "Full page image and extracted text inspected; see coverage/pages_0001_0006.json."
write(RUN / "tasks.json", tasks)

print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "source_context": len(bundle["source_context"]), "pages": len(coverage["pages"])}, ensure_ascii=False))
