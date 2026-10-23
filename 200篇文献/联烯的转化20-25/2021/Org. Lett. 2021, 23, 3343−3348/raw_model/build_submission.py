#!/usr/bin/env python3
"""Session-local transcription of Org. Lett. 2021, 23, 3343-3348.

This file contains machine transcription decisions only.  It does not assert
source review and is kept separate from the authoritative merged bundle.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DATASET = ROOT / "runs/current/dataset.json"
OUT = ROOT / "runs/current/submitted.json"
RAW = ROOT / "runs/current/raw_model"
DOC = "doc_b45c99640e5f5e97"
P = {i: f"{DOC}:page{i}:image" for i in range(1, 7)}
T = {i: f"{DOC}:page{i}:text" for i in range(1, 7)}


def f(raw=None, norm=Ellipsis, ev=(), method="original_extraction", status=None,
      derivation=None, alternatives=()):
    if norm is Ellipsis:
        norm = raw
    return st.fact(raw, norm, ev, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def q(number, unit, raw_text=None, **extra):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    value.update(extra)
    return value


def missing():
    return f()


dataset = json.loads(DATASET.read_text())
bundle = st.empty(dataset["documents"])
molecules = []
reactions = []


def add_molecule(mid, label, ev, name=None, chirality="unresolved",
                 geometry="unresolved", representation="unresolved",
                 assignment=None, mixture=None, double_bonds=None,
                 identity_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"] = [{
        "document_id": DOC,
        "compound_label": f(label, label, ev),
        "evidence_ids": list(ev),
    }]
    source_name = name or label
    m["names"]["original"] = f(source_name, source_name, ev)
    m["names"]["english"] = f(source_name, source_name, ev)
    m["names"]["abbreviations"] = f([label], [label], ev)
    m["names"]["aliases"] = f([], [], ev)
    m["cas_numbers"] = missing()
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = missing()
    m["stereochemistry"]["double_bonds"] = (
        f(double_bonds, double_bonds, ev) if double_bonds is not None else missing()
    )
    m["stereochemistry"]["axial"] = missing()
    m["stereochemistry"]["mixture"] = (
        f(mixture, mixture, ev) if mixture is not None else missing()
    )
    m["stereochemistry"]["assignment"] = (
        f(assignment, assignment, ev) if assignment is not None else missing()
    )
    m["stereochemistry"]["chirality_status"] = f(chirality, chirality, ev)
    m["stereochemistry"]["geometry_status"] = f(geometry, geometry, ev)
    m["stereochemistry"]["representation_status"] = f(representation, representation, ev)
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = missing()
    note = identity_note or "Source label retained; complete connectivity was not safely encoded from the supplied main-text drawing."
    m["identity_note"] = f(note, note, ev)
    molecules.append(m)
    return mid


def add_smiles_molecule(mid, label, name, smiles, ev):
    diag = st.stereo_inspect(smiles)
    m = add_molecule(mid, label, ev, name=name, chirality="achiral",
                     geometry="not_applicable", representation="not_applicable",
                     assignment={"basis": "common-name structure; no stereogenic element"},
                     identity_note="Candidate structure derived from the unambiguous common solvent name; parsing is not source verification.")
    m = molecules[-1]
    m["representation_type"] = "smiles"
    deriv = "Common solvent name converted to a candidate SMILES and checked with standard stereo; the source does not draw the solvent structure."
    m["smiles"] = f(name, smiles, ev, method="structure_derivation", derivation=deriv)
    m["stereochemistry"]["isomeric_smiles"] = f(
        name, diag["canonical_isomeric_smiles"], ev,
        method="structure_derivation", derivation=deriv)
    m["stereochemistry"]["identity_key"] = f(
        name, diag["stereo_graph_key"], ev,
        method="structure_derivation", derivation=deriv)
    return mid, diag


def add_template(mid, label, fixed, variables, ev, reason):
    add_molecule(mid, label, ev, name=label, chirality="unresolved",
                 geometry="unresolved", representation="not_encoded",
                 identity_note=reason)
    m = molecules[-1]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": f(fixed, {"description": fixed}, ev),
        "variables": [{
            "label": v,
            "attachment_points": f(f"Scheme attachment for {v}", f"scheme-defined attachment for {v}", ev),
            "definition": f(f"Variable {v} in the source scheme", f"source-defined {v}", ev),
            "allowed_values": missing(),
            "constraints": missing(),
        } for v in variables],
        "instance_mapping": [],
        "reason": f(reason, reason, ev),
    }
    return mid


# Generic families and substrate identities.
add_template("m_family_1", "1a-m", "allene-tethered O-propargyl oxime family shown in Scheme 2",
             ["R1", "R2"], [P[3]], "R-group family cannot be represented as one concrete SMILES.")
add_template("m_family_2", "2a-m", "trans-fused hexahydropyrrolo[3,4-b]azepin-5(4H)-one family shown in Scheme 2",
             ["R1", "R2", "exocyclic E/Z"], [P[3]], "R-group product family and E/Z mixture require a generic template.")

add_molecule("m_E1a", "E-1a", [P[2]], chirality="unspecified", geometry="specified",
             representation="not_encoded",
             double_bonds=[{"unit": "oxime C=N", "configuration": "E"}],
             assignment="Oxime configuration explicitly labelled E; propargylic stereocenter/absolute configuration not assigned.")
add_molecule("m_Z1a", "Z-1a", [P[2]], chirality="unspecified", geometry="specified",
             representation="not_encoded",
             double_bonds=[{"unit": "oxime C=N", "configuration": "Z"}],
             assignment="Oxime configuration explicitly described as Z; absolute configuration not assigned.")
add_molecule("m_ZE1a", "Z/E-1a mixture", [P[2]], chirality="unspecified", geometry="mixture",
             representation="not_encoded",
             mixture="Mixture of Z/E starting-material isomers; composition not reported.",
             assignment="Starting-material mixture explicitly described; no component ratio reported.")
for letter in "bcdefghijklm":
    add_molecule(f"m_1{letter}", f"1{letter}", [P[3]], chirality="unspecified",
                 geometry="unresolved", representation="not_encoded",
                 assignment="Scope substrate label inferred from product row; starting-material E/Z configuration is not legible as an individual structure in the main-text scope panel.")
add_molecule("m_1n", "1n", [P[4]], chirality="unspecified", geometry="unspecified",
             representation="not_encoded",
             assignment="Scheme 3 drawing does not assign absolute configuration or a unique oxime geometry.")


def add_product_set(letter, ev):
    mix = f"m_2{letter}_mix"
    z = f"m_2{letter}_Z"
    e = f"m_2{letter}_E"
    common = "trans-fused relative stereochemistry drawn; absolute configuration not assigned"
    add_molecule(mix, f"2{letter} (Z/E mixture)", ev, chirality="partially_specified",
                 geometry="mixture", representation="not_encoded",
                 mixture=f"Exocyclic Z/E mixture of 2{letter}; ratio is experiment-specific.",
                 assignment=common)
    add_molecule(z, f"trans-Z-2{letter}", ev, chirality="partially_specified",
                 geometry="specified", representation="not_encoded",
                 double_bonds=[{"unit": "exocyclic alkene", "configuration": "Z"}],
                 assignment=common + "; Z component registered separately")
    add_molecule(e, f"trans-E-2{letter}", ev, chirality="partially_specified",
                 geometry="specified", representation="not_encoded",
                 double_bonds=[{"unit": "exocyclic alkene", "configuration": "E"}],
                 assignment=common + "; E component registered separately")
    return mix, z, e


product_sets = {}
for letter in "abcdefghijklm":
    product_sets[letter] = add_product_set(letter, [P[2], P[3]] if letter == "a" else [P[3]])
add_molecule("m_complex_1n", "complex product mixture from 1n", [P[4], P[6]],
             chirality="unresolved", geometry="unresolved", representation="unresolved",
             mixture="Complex product mixture; individual components and proportions unreported.",
             assignment="The source labels the outcome as a complex mixture; the single drawn scaffold is not treated as a uniquely isolated product.")

# Catalysts and solvents.  Only common solvents receive conservative candidate SMILES.
condition_names = {
    "m_CuCl": ("CuCl", "copper(I) chloride"),
    "m_AgOTf": ("AgOTf", "silver trifluoromethanesulfonate"),
    "m_PPh3AuCl": ("PPh3AuCl", "triphenylphosphine gold(I) chloride"),
    "m_AuCl3": ("AuCl3", "gold(III) chloride"),
    "m_InOTf3": ("In(OTf)3", "indium(III) trifluoromethanesulfonate"),
    "m_CuBr": ("CuBr", "copper(I) bromide"),
    "m_CuCODCl2": ("[Cu(COD)Cl]2", "chlorobis(1,5-cyclooctadiene)copper(I) dimer"),
    "m_CuMeCN4BF4": ("Cu(CH3CN)4BF4", "tetrakis(acetonitrile)copper(I) tetrafluoroborate"),
}
for mid, (label, name) in condition_names.items():
    add_molecule(mid, label, [P[2]], name=name, chirality="unresolved",
                 geometry="not_applicable", representation="unresolved",
                 identity_note="Catalyst identity retained from the source label; coordination/ionic representation not forced into SMILES.")
stereo_diagnostics = {}
for mid, label, name, smiles in [
    ("m_toluene", "toluene", "toluene", "Cc1ccccc1"),
    ("m_MeCN", "MeCN", "acetonitrile", "CC#N"),
    ("m_THF", "THF", "tetrahydrofuran", "C1CCOC1"),
    ("m_DMF", "DMF", "N,N-dimethylformamide", "CN(C)C=O"),
]:
    _, diag = add_smiles_molecule(mid, label, name, smiles, [P[2]])
    stereo_diagnostics[mid] = diag


def participant(pid, mid, role, ev, amount=None, mol_percent=None, ref=None,
                product_role=None):
    p = st.new_participant(pid, mid, role, ev)
    if product_role is not None:
        p["product_role"] = f(product_role, product_role, ev)
    if amount is not None:
        p["charges"]["amount"] = f(amount[0], amount[1], ev)
    if mol_percent is not None:
        p["charges"]["mol_percent"] = f(mol_percent[0], mol_percent[1], ev)
    if ref is not None:
        p["equivalence_reference"] = f(ref, ref, ev)
    return p


def reaction_step(step_id, ev, pids, temperature=None, time=None, operation=None):
    s = st.new_step(step_id, ev)
    s["participant_ids"] = list(pids)
    if operation is not None:
        s["operation"] = f(operation, operation, ev)
    if temperature is not None:
        s["temperature"] = f(temperature[0], temperature[1], ev)
    if time is not None:
        s["time"] = f(time[0], time[1], ev)
    return s


def outcome(oid, kind, targets, ev, raw=None, norm=Ellipsis, basis=None,
            yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, ev)
    if raw is not None or norm is not Ellipsis:
        o["value"] = f(raw, norm, ev)
    if basis is not None:
        o["measurement_basis"] = f(basis, basis, ev)
    if yield_type is not None:
        o["yield_type"] = f(yield_type, yield_type, ev)
    if scope is not None:
        o["scope"] = f(scope, scope, ev)
    return o


def cascade_classification(ev):
    deriv = "Partial net-bond inventory derived from the authors' statement that the transformation constructs one C=O, one C=N, and two C-C bonds; no atom map is asserted."
    changes = [
        {"change_id": "c1", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "c2", "atoms": [{"element": "C"}, {"element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "c3", "atoms": [{"element": "C"}, {"element": "O"}], "before_order": 1, "after_order": 2, "change_type": "order_increased"},
        {"change_id": "c4", "atoms": [{"element": "C"}, {"element": "N"}], "before_order": 0, "after_order": 2, "change_type": "formed"},
    ]
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["c1", "c2"]},
        {"category_number": 14, "code": "DB-C-O", "supporting_change_ids": ["c3"]},
        {"category_number": 16, "code": "DB-C-N", "supporting_change_ids": ["c4"]},
    ]
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [f(None, c, ev, method="structure_derivation", derivation=deriv) for c in changes],
        "labels": f(None, labels, ev, method="structure_derivation", derivation=deriv),
        "coverage": f("Partial classification from authors' bond-count statement", "partial", ev,
                      method="inference", derivation=deriv),
        "limitations": ["No safe full atom mapping from the main-text drawings.",
                        "Bond cleavages, hydrogen transfers, and coproduct balance remain unresolved."],
        "model_votes": [],
    }


def unresolved_classification(ev, reason, not_applicable=False):
    coverage = "not_applicable" if not_applicable else "unresolved"
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": missing(),
        "bond_changes": [],
        "labels": f(reason, [], ev, method="inference", derivation=reason) if not_applicable else missing(),
        "coverage": f(reason, coverage, ev, method="inference", derivation=reason),
        "limitations": [reason],
        "model_votes": [],
    }


def base_reaction(rid, ev, record_type, categories, location, react_name,
                  transform, unresolved):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, ev)
    r["record_type"] = record_type
    r["experiment_categories"] = f(categories, categories, ev)
    r["location"] = f(location, location, ev)
    r["reaction_name"] = f(react_name, react_name, ev)
    r["transformation_type"] = f(transform, transform, ev)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


def set_common_conditions(r, ev, substrate_mid, catalyst_mid, solvent_mid,
                          temp, amount_mmol=None, concentration=None,
                          time=None, product_set=None, catalyst_loading=5):
    r["participants"] = [participant("substrate", substrate_mid, "reactant", ev)]
    if amount_mmol is not None:
        r["participants"][0]["charges"]["amount"] = f(
            f"{amount_mmol} mmol", q(amount_mmol, "mmol", f"{amount_mmol} mmol"), ev)
    r["participants"].append(participant(
        "catalyst", catalyst_mid, "catalyst", ev,
        mol_percent=(f"{catalyst_loading} mol%", q(catalyst_loading, "mol%", f"{catalyst_loading} mol%")),
        ref=substrate_mid))
    r["participants"].append(participant("solvent", solvent_mid, "solvent", ev))
    if product_set is not None:
        mix, z, e = product_set
        r["participants"] += [
            participant("product_mix", mix, "product", ev, product_role="combined isolated product mixture"),
            participant("product_Z", z, "product", ev, product_role="Z component"),
            participant("product_E", e, "product", ev, product_role="E component"),
        ]
    pids = [p["participant_id"] for p in r["participants"]]
    r["steps"] = [reaction_step(
        "reaction", ev, pids,
        temperature=(f"{temp} °C", q(temp, "°C", f"{temp} °C")),
        time=(f"{time} h", q(time, "h", f"{time} h")) if time is not None else None,
        operation="cycloisomerization under the tabulated/specified conditions")]
    solvent_label = next(m["occurrences"][0]["compound_label"]["normalized_value"] for m in molecules if m["molecule_id"] == solvent_mid)
    r["solvent_conditions"]["composition"] = f(solvent_label, [solvent_mid], ev)
    if concentration is not None:
        r["solvent_conditions"]["concentration"] = f(
            f"{concentration} M", q(concentration, "M", f"{concentration} M"), ev)
        r["solvent_conditions"]["concentration_basis"] = f(
            "1a concentration", substrate_mid, ev)
    if amount_mmol is not None:
        r["reaction_scale"] = f(
            f"{amount_mmol} mmol substrate",
            {"reference_molecule_id": substrate_mid, "amount": q(amount_mmol, "mmol", f"{amount_mmol} mmol")}, ev)


# Generic current-work diagrams.
for rid, page, location in [
    ("scheme1b_generic", 2, {"page": 2, "region": "Scheme 1b; graphical abstract on page 1"}),
    ("scheme2_generic", 3, {"page": 3, "region": "Scheme 2 generic arrow"}),
]:
    ev = [P[page]] + ([P[1]] if rid == "scheme1b_generic" else [])
    r = base_reaction(rid, ev, "generic_scheme", [], location,
                      "Copper-catalyzed cycloisomerization of unactivated allene-tethered O-propargyl oximes",
                      "[2,3]-sigmatropic rearrangement/[3+2] cycloaddition/[3,3]-sigmatropic rearrangement cascade",
                      ["Generic R-group structures are not concrete experiments.", "Exact atom mapping is unresolved."])
    r["participants"] = [
        participant("substrate_family", "m_family_1", "reactant", ev),
        participant("catalyst", "m_CuBr", "catalyst", ev,
                    mol_percent=("5 mol%", q(5, "mol%", "5 mol%")), ref="m_family_1"),
        participant("solvent", "m_toluene", "solvent", ev),
        participant("product_family", "m_family_2", "product", ev, product_role="product family"),
    ]
    r["steps"] = [reaction_step("generic_reaction", ev,
                  [p["participant_id"] for p in r["participants"]],
                  temperature=("100 °C", q(100, "°C", "100 °C")),
                  operation="generic copper-catalyzed cycloisomerization")]
    r["solvent_conditions"]["composition"] = f("toluene", ["m_toluene"], ev)
    r["classification"] = cascade_classification(ev + [P[4]])
    reactions.append(r)


# Table 1: every row is an independent optimization experiment.
rows = [
    (1, "m_CuCl", "m_toluene", 100, 44, "78:22", 12),
    (2, "m_AgOTf", "m_toluene", 100, None, None, None),
    (3, "m_PPh3AuCl", "m_toluene", 100, None, None, None),
    (4, "m_AuCl3", "m_toluene", 100, None, None, None),
    (5, "m_InOTf3", "m_toluene", 100, 12, None, None),
    (6, "m_CuBr", "m_toluene", 100, 69, "79:21", None),
    (7, "m_CuCODCl2", "m_toluene", 100, 64, "79:21", None),
    (8, "m_CuMeCN4BF4", "m_toluene", 100, 55, "78:22", None),
    (9, "m_CuBr", "m_toluene", 80, 40, "79:21", None),
    (10, "m_CuBr", "m_toluene", 120, 66, "81:19", None),
    (11, "m_CuBr", "m_MeCN", 100, 22, "86:14", None),
    (12, "m_CuBr", "m_THF", 100, 35, "80:20", None),
    (13, "m_CuBr", "m_DMF", 100, 40, "80:20", None),
    (14, None, "m_toluene", 100, None, None, None),
]
for entry, cat, solv, temp, yld, ratio, hours in rows:
    ev = [P[2], T[2]]
    r = base_reaction(f"table1_entry{entry:02d}", ev, "specific_experiment",
                      ["optimization"] + (["substrate_scope"] if entry == 6 else []),
                      {"page": 2, "region": "Table 1", "entry": entry},
                      "Optimization of the copper-catalyzed cycloisomerization",
                      "cycloisomerization", ["Workup and purification are unavailable without the SI."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1b_generic", f"{DOC}:scheme2_generic"]
    if entry == 14:
        r["participants"] = [participant("substrate", "m_E1a", "reactant", ev),
                             participant("solvent", solv, "solvent", ev)]
        r["participants"][0]["charges"]["amount"] = f("0.2 mmol", q(0.2, "mmol", "0.2 mmol"), ev)
        r["steps"] = [reaction_step("reaction", ev, ["substrate", "solvent"],
                      temperature=("100 °C", q(100, "°C", "100 °C")),
                      operation="control without catalyst")]
        r["solvent_conditions"]["composition"] = f("toluene", [solv], ev)
        r["solvent_conditions"]["concentration"] = f("0.1 M", q(0.1, "M", "0.1 M"), ev)
        r["solvent_conditions"]["concentration_basis"] = f("1a concentration", "m_E1a", ev)
        r["reaction_scale"] = f("1a (0.2 mmol)", {"reference_molecule_id": "m_E1a", "amount": q(0.2, "mmol", "0.2 mmol")}, ev)
        r["outcomes"] = [
            outcome("no_reaction", "no_reaction", [], ev, "n.r.", "n.r.", scope="Table 1 entry 14"),
            outcome("starting_material", "recovery", ["m_E1a"], ev,
                    "The starting material remained untouched", "starting material remained untouched",
                    scope="qualitative recovery; no percentage reported"),
        ]
        r["classification"] = unresolved_classification(ev, "No observable transformation; no formed-bond label is applicable.", True)
        r["review"]["unresolved"] += ["Reaction duration and quantitative recovery were not reported."]
    else:
        set_common_conditions(r, ev, "m_E1a", cat, solv, temp,
                              amount_mmol=0.2, concentration=0.1, time=hours,
                              product_set=product_sets["a"] if yld is not None else None)
        r["outcomes"] = []
        if yld is not None:
            r["outcomes"].append(outcome(
                "yield", "yield", [product_sets["a"][0]], ev,
                f"{yld}%", q(yld, "%", f"{yld}%"), basis="isolated",
                yield_type="isolated yield", scope=f"Table 1 entry {entry}"))
        if ratio is not None:
            if entry == 6:
                alt_table = f("79:21 (Z/E), Table 1 entry 6", {"order": ["Z", "E"], "values": [79, 21]}, ev)
                alt_scope = f("78:22 (Z/E), Scheme 2, 2a", {"order": ["Z", "E"], "values": [78, 22]}, [P[3]])
                value_fact = f(["79:21", "78:22"], None, [P[2], P[3]], status="conflict",
                               alternatives=[alt_table, alt_scope])
                sel = outcome("selectivity", "selectivity", [product_sets["a"][1], product_sets["a"][2]], [P[2], P[3]],
                              scope={"metric": "Z/E", "order": [product_sets["a"][1], product_sets["a"][2]]})
                sel["value"] = value_fact
                r["outcomes"].append(sel)
                r["evidence_ids"] = [P[2], T[2], P[3]]
                r["location"] = f({"table": "Table 1 entry 6", "scope": "Scheme 2 product 2a"},
                                  {"page_regions": [{"page": 2, "region": "Table 1 entry 6"}, {"page": 3, "region": "Scheme 2 product 2a"}]}, [P[2], P[3]])
                r["review"]["unresolved"] += ["Same 69% experiment reports Z/E 79:21 in Table 1 but 78:22 in Scheme 2; neither value was chosen."]
            else:
                a, b = [int(x) for x in ratio.split(":")]
                r["outcomes"].append(outcome(
                    "selectivity", "selectivity", [product_sets["a"][1], product_sets["a"][2]], ev,
                    ratio, {"order": ["Z", "E"], "values": [a, b]},
                    scope={"metric": "Z/E", "order": [product_sets["a"][1], product_sets["a"][2]]}))
        if entry in (2, 3, 4, 5):
            r["outcomes"].append(outcome(
                "decomposition", "decomposition", [], ev,
                "mixture of unknown decomposition products", "mixture of unknown decomposition products",
                scope="author discussion of Table 1 entries 2-5"))
        if entry in (2, 3, 4):
            r["classification"] = unresolved_classification(ev, "Only unknown decomposition products were reported; product identity and bond changes are unresolved.")
            r["review"]["unresolved"] += ["Table cells for yield and Z/E are blank; no product identity or measurement basis was assigned."]
        else:
            r["classification"] = cascade_classification([P[2], P[4]])
        if hours is None:
            r["review"]["unresolved"] += ["Reaction time is not stated for this table row; 12 h from entry 1 was not silently inherited."]
    row_obj = {"entry": entry, "catalyst": condition_names.get(cat, (None,))[0] if cat else None,
               "solvent": next(m["occurrences"][0]["compound_label"]["normalized_value"] for m in molecules if m["molecule_id"] == solv),
               "temperature_C": temp, "yield_percent": yld, "Z_E": ratio}
    r["procedure"]["original_text"] = f(row_obj, row_obj, [P[2]])
    reactions.append(r)


# Two additional experiments explicitly described in the page-2 prose.
for rid, substrate, raw_label, ratio_order, ratio_values in [
    ("prose_Z1a_optimized", "m_Z1a", "Z isomer of the starting material", "E/Z", [78, 22]),
    ("prose_ZE1a_optimized", "m_ZE1a", "mixture of Z/E isomers of the starting material", "Z/E", [79, 21]),
]:
    ev = [P[2], T[2]]
    r = base_reaction(rid, ev, "specific_experiment", ["optimization"],
                      {"page": 2, "region": "paragraph below Table 1"},
                      "Cycloisomerization under optimized conditions", "cycloisomerization",
                      ["Reaction scale, time, workup, purification, and yield measurement basis are unreported."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1b_generic"]
    set_common_conditions(r, ev, substrate, "m_CuBr", "m_toluene", 100,
                          product_set=product_sets["a"])
    r["outcomes"] = [
        outcome("yield", "yield", [product_sets["a"][0]], ev, "66%", q(66, "%", "66%"),
                scope=raw_label),
        outcome("selectivity", "selectivity", [product_sets["a"][1], product_sets["a"][2]], ev,
                f"{ratio_values[0]}:{ratio_values[1]} {ratio_order}",
                {"order": ratio_order.split("/"), "values": ratio_values},
                scope={"metric": ratio_order, "order": ratio_order.split("/")}),
    ]
    r["procedure"]["original_text"] = f(
        f"{raw_label}; 5 mol % CuBr in toluene at 100 °C; 66% yield; {ratio_values[0]}:{ratio_values[1]} {ratio_order}",
        {"substrate_description": raw_label, "CuBr_mol_percent": 5, "solvent": "toluene", "temperature_C": 100}, ev)
    r["classification"] = cascade_classification([P[2], P[4]])
    reactions.append(r)


# Scheme 2 substrate-scope rows 2b-2m.  Product 2a is merged with Table 1 entry 6.
scope_rows = {
    "b": (62, [76, 24]), "c": (60, [71, 29]), "d": (81, [72, 28]),
    "e": (69, [77, 23]), "f": (70, [75, 25]), "g": (51, [81, 19]),
    "h": (52, [79, 21]), "i": (64, [72, 28]), "j": (76, [74, 26]),
    "k": (62, [86, 14]), "l": (57, [85, 15]), "m": (22, [73, 27]),
}
for letter, (yld, ze) in scope_rows.items():
    ev = [P[3]]
    r = base_reaction(f"scheme2_scope_2{letter}", ev, "specific_experiment",
                      ["substrate_scope"], {"page": 3, "region": "Scheme 2", "product": f"2{letter}"},
                      "Copper-catalyzed cycloisomerization of allene-tethered O-propargyl oxime",
                      "cycloisomerization", ["Reaction scale, time, workup, purification, and yield measurement basis are not reported on the supplied pages.",
                      "Complete reactant/product connectivity and absolute stereochemistry require structure-level review or SI."])
    r["general_scheme_ids"] = [f"{DOC}:scheme2_generic"]
    set_common_conditions(r, ev, f"m_1{letter}", "m_CuBr", "m_toluene", 100,
                          product_set=product_sets[letter])
    r["outcomes"] = [
        outcome("yield", "yield", [product_sets[letter][0]], ev,
                f"{yld}%", q(yld, "%", f"{yld}%"), scope=f"Scheme 2 product 2{letter}"),
        outcome("selectivity", "selectivity", [product_sets[letter][1], product_sets[letter][2]], ev,
                f"{ze[0]}:{ze[1]}", {"order": ["Z", "E"], "values": ze},
                scope={"metric": "Z/E", "order": [product_sets[letter][1], product_sets[letter][2]]}),
    ]
    r["procedure"]["original_text"] = f(
        {"scheme": 2, "substrate": f"1{letter}", "product": f"2{letter}", "CuBr_mol_percent": 5,
         "solvent": "toluene", "temperature_C": 100, "yield_percent": yld, "Z_E": ze},
        {"scheme": 2, "substrate": f"1{letter}", "product": f"2{letter}"}, ev)
    r["classification"] = cascade_classification([P[3], P[4]])
    reactions.append(r)


# Scheme 3 is a real attempted experiment with an unresolved complex mixture.
ev = [P[4], P[6], T[6]]
r = base_reaction("scheme3_1n_complex_mixture", ev, "specific_experiment",
                  ["mechanistic_study"], {"page": 4, "region": "Scheme 3; explanatory note 15 on page 6"},
                  "Cycloisomerization of 1n under the optimized reaction condition", "cycloisomerization attempt",
                  ["Individual products, proportions, yield, reaction scale, duration, workup, and analytical basis are unreported.",
                   "The drawn possible scaffold is not treated as a uniquely identified isolated product."])
r["general_scheme_ids"] = [f"{DOC}:scheme1b_generic"]
set_common_conditions(r, ev, "m_1n", "m_CuBr", "m_toluene", 100)
r["participants"].append(participant("product_mixture", "m_complex_1n", "product", ev,
                                     product_role="unresolved complex mixture"))
r["steps"][0]["participant_ids"].append("product_mixture")
r["outcomes"] = [outcome("complex_mixture", "other", ["m_complex_1n"], ev,
                         "complex mixture", "complex product mixture",
                         scope="Scheme 3 outcome; no individual product yield")]
r["procedure"]["original_text"] = f(
    "1n, CuBr (5 mol%), toluene, 100 °C; complex mixture",
    {"substrate": "1n", "CuBr_mol_percent": 5, "solvent": "toluene", "temperature_C": 100}, ev)
r["classification"] = unresolved_classification(ev, "Complex products were not identified; net bond changes cannot be assigned.")
r["notes"] = [f("probably due to the inductive effect of the methyl group",
                "Author-proposed explanation: inductive effect of the methyl group", [P[6], T[6]])]
reactions.append(r)


# Figure 2/3 mechanism is a proposal/DFT analysis, not a preparative experiment.
ev = [P[2], P[3], P[4], T[3], T[4]]
r = base_reaction("figure2_proposed_mechanism", ev, "proposed_mechanism", ["mechanistic_study"],
                  {"pages": [2, 3, 4], "regions": ["mechanistic prose", "Figure 2", "Figure 3"]},
                  "Proposed [2,3]/[3+2]/[3,3] domino mechanism",
                  "proposed sigmatropic rearrangement/cycloaddition/rearrangement/tautomerization sequence",
                  ["Mechanistic structures and transition states are proposals/calculations, not isolated intermediates.",
                   "Exact atom mapping and hydrogen transfers remain unresolved."])
r["participants"] = [participant("starting_material", "m_E1a", "reactant", ev),
                     participant("catalyst", "m_CuBr", "catalyst", ev),
                     participant("product_mix", product_sets["a"][0], "product", ev, product_role="final product mixture")]
r["steps"] = [
    reaction_step("proposed_2_3_rearrangement", ev, ["starting_material", "catalyst"], operation="proposed [2,3]-sigmatropic rearrangement to E-3/Z-3 nitrones"),
    reaction_step("proposed_3_2_cycloaddition", ev, ["starting_material"], operation="proposed intramolecular [3+2] dipolar cycloaddition to Int.A1/Int.A2"),
    reaction_step("proposed_3_3_rearrangement", ev, ["starting_material"], operation="proposed [3,3] hetero-Claisen rearrangement to cis-E-2/cis-Z-2"),
    reaction_step("proposed_tautomerization", ev, ["product_mix"], operation="proposed tautomerization to trans-E-2/trans-Z-2"),
]
r["classification"] = cascade_classification([P[3], P[4]])
r["notes"] = [
    f("calculated Z/E based on Boltzmann populations is 71:29", {"metric": "calculated Z/E", "values": [71, 29]}, [P[3], T[3]]),
    f("1,3-DC barrier 22.2 kcal/mol; [3,3] rearrangement barrier 12.7 kcal/mol",
      {"DC_barrier_kcal_mol": 22.2, "rearrangement_barrier_kcal_mol": 12.7}, [P[3], P[4], T[3]]),
]
reactions.append(r)


bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page1_graphical_abstract", "kind": "duplicate_summary",
     "evidence_ids": [P[1]], "note": "Graphical abstract duplicates the current-work generic reaction; not a separate experiment."},
    {"context_id": f"{DOC}:page1_figure1", "kind": "background_structures",
     "evidence_ids": [P[1]], "note": "Examples of biologically active azepines; not current synthesis experiments."},
    {"context_id": f"{DOC}:scheme1a_prior_work", "kind": "literature_context",
     "evidence_ids": [P[2]], "note": "Four previous-work O-propargylic oxime transformations; not experiments from this paper."},
    {"context_id": f"{DOC}:scheme2_xray", "kind": "characterization",
     "evidence_ids": [P[3]], "note": "ORTEP/CCDC-1810670 characterization for 2a; not a separate synthesis."},
    {"context_id": f"{DOC}:figure2_transition_states", "kind": "computational_analysis",
     "evidence_ids": [P[3]], "note": "Eight calculated transition-state structures and relative Gibbs energies; no preparative records."},
    {"context_id": f"{DOC}:figure3_pes", "kind": "computational_analysis",
     "evidence_ids": [P[4]], "note": "Potential-energy surface supporting the proposed mechanism; no preparative record."},
    {"context_id": f"{DOC}:si_notice", "kind": "supplement_notice",
     "evidence_ids": [P[4]], "note": "The paper states SI exists, but no SI was supplied or explicitly associated in this workspace."},
    {"context_id": f"{DOC}:references", "kind": "bibliography",
     "evidence_ids": [P[4], P[5], P[6]], "note": "Bibliography and note 15; only note 15 supplies duplicate context for the Scheme 3 experiment."},
]

OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RAW / "stereo_diagnostics.json").write_text(json.dumps(stereo_diagnostics, ensure_ascii=False, indent=2) + "\n")
(RAW / "source_mentions_1a_2a_1n.json").write_text(json.dumps(
    st.source_mentions(dataset, ["1a", "2a", "1n"], context=320),
    ensure_ascii=False, indent=2) + "\n")
(RAW / "extraction_inventory.json").write_text(json.dumps({
    "document_id": DOC,
    "title": dataset["documents"][0]["title"],
    "pages_viewed": [1, 2, 3, 4, 5, 6],
    "molecule_records": len(molecules),
    "reaction_records": len(reactions),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in reactions),
    "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in reactions),
    "decisions": {
        "scheme2_2a": "Merged with Table 1 entry 6 because conditions and 69% yield coincide; conflicting Z/E values retained.",
        "si": "Unestablished; no supplement imported.",
        "structures": "Complex primary structures remain unresolved rather than forcing unsafe SMILES.",
    },
}, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions)}, ensure_ascii=False))
