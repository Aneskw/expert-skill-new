#!/usr/bin/env python3
"""Build the source-derived v3.1 delta for original PDF pages 7--9.

The script starts from the standard-init template and updates only records
whose source support changes on pages 7--9.  All output remains machine
generated, partial, and needs expert review.
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


DOC = "doc_03086fff1f443524"
E = {p: f"{DOC}:page{p}:image" for p in range(1, 10)}
ET = {p: f"{DOC}:page{p}:text" for p in range(1, 10)}
OLD = json.loads((ROOT / "runs/current/standard_records.json").read_text())
MOLS = {m["molecule_id"]: m for m in OLD["molecules"]}
RXNS = {r["reaction_id"]: r for r in OLD["reactions"]}
STEREO_CHECKS = {}


def source(raw, normalized=None, ev=(), *, method="original_extraction", status="pending_review"):
    if normalized is None:
        normalized = raw
    return st.fact(raw, normalized, ev, method=method, status=status)


def inferred(normalized, ev, derivation, raw=None):
    return st.fact(raw, normalized, ev, method="structure_derivation", status="pending_review",
                   verification="inferred", derivation=derivation)


def missing(status="unreported"):
    return st.fact(status=status)


def not_applicable():
    return st.fact(method="not_applicable", status="not_applicable")


def quantity(number, unit, *, raw=None, qualifier=None, extra=None):
    value = {"number": number, "unit": unit}
    if qualifier:
        value["qualifier"] = qualifier
    if raw:
        value["raw_text"] = raw
    if extra:
        value.update(extra)
    return value


def set_charge(p, key, raw, normalized, ev, *, method="original_extraction"):
    p["charges"][key] = source(raw, normalized, ev, method=method)


def participant(pid, mid, role, ev, *, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role is not None:
        p["product_role"] = source(product_role, product_role, ev)
    return p


def step(sid, ev, pids, operation, *, stage="reaction", temperature=None, time=None, method="original_extraction"):
    s = st.new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = source(operation, operation, ev, method=method)
    if temperature is not None:
        s["temperature"] = source(str(temperature), temperature, ev, method=method)
    if time is not None:
        s["time"] = source(str(time), time, ev, method=method)
    return s


def outcome(oid, kind, targets, ev, raw, normalized, *, product_role=None,
            yield_type=None, measurement_basis=None, scope=None, method="original_extraction"):
    o = st.new_outcome(oid, kind, targets, ev)
    o["value"] = source(raw, normalized, ev, method=method)
    if product_role is not None:
        o["product_role"] = source(product_role, product_role, ev, method=method)
    if yield_type is not None:
        o["yield_type"] = source(yield_type, yield_type, ev, method=method)
    if measurement_basis is not None:
        o["measurement_basis"] = source(measurement_basis, measurement_basis, ev, method=method)
    if scope is not None:
        o["scope"] = source(scope, scope, ev, method=method)
    return o


def append_evidence(items, *ids):
    for eid in ids:
        if eid not in items:
            items.append(eid)


def unresolved_molecule(mid, label, name, ev, *, chirality="unresolved", geometry="unresolved",
                        representation="not_encoded", assignment=None, mixture=None, axial=None,
                        identity_note=None):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = source(label, label, ev)
    m["names"]["original"] = source(name, name, ev)
    m["names"]["english"] = source(name, name, ev)
    m["cas_numbers"] = missing()
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = missing()
    m["stereochemistry"]["double_bonds"] = missing()
    m["stereochemistry"]["axial"] = source(axial, axial, ev) if axial else missing()
    m["stereochemistry"]["mixture"] = source(mixture, mixture, ev) if mixture else missing()
    m["stereochemistry"]["assignment"] = source(assignment, assignment, ev) if assignment else missing()
    m["stereochemistry"]["chirality_status"] = inferred(
        chirality, ev, "Status assigned conservatively after inspecting the source identity; no complete stereochemical graph was encoded.")
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, ev, "Status assigned conservatively after inspecting double-bond, ring, and axial information in the source region.")
    m["stereochemistry"]["representation_status"] = inferred(
        representation, ev, "The source identity is recognizable, but a complete source-checked stereochemical representation was not asserted.")
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = missing()
    m["identity_note"] = source(identity_note, identity_note, ev) if identity_note else missing()
    return m


def smiles_molecule(mid, label, name, smiles, ev, identity_note):
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = source(label, label, ev)
    m["names"]["original"] = source(name, name, ev)
    m["names"]["english"] = source(name, name, ev)
    m["cas_numbers"] = missing()
    diag = st.stereo_inspect(smiles)
    STEREO_CHECKS[mid] = diag
    m["representation_type"] = "smiles"
    m["smiles"] = inferred(smiles, ev, "Candidate connectivity derived from the explicit source name; parsing is not source verification.")
    m["generic_template"] = None
    m["stereochemistry"]["atom_centers"] = inferred(diag["atom_centers"], ev, "Standard stereo graph diagnostic on the candidate SMILES.")
    m["stereochemistry"]["double_bonds"] = inferred(diag["double_bonds"], ev, "Standard stereo graph diagnostic on the candidate SMILES.")
    m["stereochemistry"]["axial"] = inferred([], ev, "No axial assignment is encoded by the candidate; toolkit axial support is limited.")
    m["stereochemistry"]["mixture"] = not_applicable()
    atom_unspecified = any(x["type"].startswith("Atom_") and x["specified"] == "Unspecified" for x in diag["potential_stereo"])
    bond_unspecified = any(x["type"].startswith("Bond_") and x["specified"] == "Unspecified" for x in diag["potential_stereo"])
    m["stereochemistry"]["assignment"] = inferred(
        "Potential stereogenic elements remain unassigned in the candidate graph." if atom_unspecified or bond_unspecified
        else "No potential stereogenic element detected in the candidate graph.", ev,
        "Statement is based only on the standard stereo graph diagnostic.")
    m["stereochemistry"]["chirality_status"] = inferred(
        "unspecified" if atom_unspecified else "achiral", ev,
        "Potential atom stereocenters are unassigned." if atom_unspecified else "No potential atom stereocenter was found in the candidate graph.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "unspecified" if bond_unspecified else "not_applicable", ev,
        "Potential bond/cumulene stereochemistry is unassigned." if bond_unspecified else "No E/Z-applicable stereogenic bond was found in the candidate graph.")
    m["stereochemistry"]["representation_status"] = inferred(
        "partial" if atom_unspecified or bond_unspecified else "encoded", ev,
        "Candidate leaves potential stereochemistry unassigned." if atom_unspecified or bond_unspecified
        else "Candidate graph and canonical isomeric SMILES were inspected by the standard stereo tool.")
    m["stereochemistry"]["isomeric_smiles"] = inferred(diag["canonical_isomeric_smiles"], ev, "Canonical isomeric SMILES from standard stereo inspection.")
    m["stereochemistry"]["identity_key"] = inferred(diag["stereo_graph_key"], ev, "Stereo graph key from standard stereo inspection; not a sample identity.")
    m["identity_note"] = source(identity_note, identity_note, ev)
    return m


def generic_zw():
    ev = [E[7]]
    m = unresolved_molecule(
        "m_ZW", "ZW", "generic zwitterionic intermediate ZW", ev,
        chirality="unresolved", geometry="unresolved", representation="not_encoded",
        assignment="Generic R-group structure; configuration depends on R1/R2/R3.",
        identity_note="Scheme 9 generic zwitterion with R1/R2/R3; not a concrete compound identity.")
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": source("Scheme 9 drawing of ZW", {"source_drawing": "ZW", "description": "allylic zwitterion with N-oxide and R-group variables"}, ev),
        "variables": [
            {"label": r, "attachment_points": source(f"{r} attachment shown in Scheme 9", f"source-drawn {r} attachment", ev),
             "definition": st.fact(f"{r} not defined on page 7", None, ev, status="unreadable"),
             "allowed_values": missing(), "constraints": missing()}
            for r in ("R1", "R2", "R3")
        ],
        "instance_mapping": [],
        "reason": source("Scheme 9 uses R1/R2/R3", "Generic template retained because substituents are not instantiated.", ev),
    }
    return m


def add_occurrence(m, label, ev):
    if not any(o["evidence_ids"] == list(ev) for o in m["occurrences"]):
        m["occurrences"].append({"document_id": DOC, "compound_label": source(label, label, ev), "evidence_ids": list(ev)})


def make_reaction(suffix, ev, location, name, transformation, categories, unresolved):
    r = st.new_reaction(f"{DOC}:{suffix}", DOC, ev)
    r["record_type"] = "specific_experiment"
    r["experiment_categories"] = source(categories, categories, ev)
    r["location"] = source(location, location, ev)
    r["reaction_name"] = source(name, name, ev)
    r["transformation_type"] = inferred(transformation, ev, "Transformation class inferred from the explicit source procedure and product identity.")
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": list(unresolved)}
    r["completeness"] = "partial"
    return r


ANNULATION_CHANGES = [
    {"change_id": "b1", "stage": "nitrone addition", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b2", "stage": "aqueous workup intermediate", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
    {"change_id": "b3", "stage": "cyclization proton transfer", "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}], "before_order": 1, "after_order": 0, "change_type": "hydrogen_loss", "hydrogen_count": 1},
    {"change_id": "b4", "stage": "6-endo-trig cyclization", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    {"change_id": "b5", "stage": "cyclization proton transfer", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
]


def set_classification(r, ev, changes, limitation):
    labels = {}
    for c in changes:
        if c["change_type"] in {"formed", "hydrogen_gain", "order_increased"}:
            lab = st.bond_label(c["atoms"][0]["element"], c["atoms"][1]["element"], c["after_order"])
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(c["change_id"])
    r["classification"]["bond_changes"] = [inferred(c, ev, "Stage-level element-only change inferred from the source procedure/drawing; no complete atom mapping asserted.") for c in changes]
    r["classification"]["labels"] = inferred(list(labels.values()), ev, "Labels include only explicitly retained formed/promoted stage bonds.")
    r["classification"]["coverage"] = inferred("partial", ev, "No complete atom-mapped reaction graph was asserted.")
    r["classification"]["limitations"] = [limitation]


def add_general_charges(parts, alkoxy_id, nitrone_id, nitrone_mass, ev_specific, *, alkoxy_amount=3.50):
    ev_general = [E[7], ET[7]]
    byid = {p["participant_id"]: p for p in parts}
    set_charge(byid[alkoxy_id], "amount", f"{alkoxy_amount:.2f} mmol (general procedure)", quantity(alkoxy_amount, "mmol"), ev_general, method="inherited_general_condition")
    set_charge(byid["nBuLi"], "volume", "1.2 mL (general procedure)", quantity(1.2, "mL"), ev_general, method="inherited_general_condition")
    set_charge(byid["nBuLi"], "amount", "3.0 mmol (general procedure)", quantity(3.0, "mmol"), ev_general, method="inherited_general_condition")
    byid["nBuLi"]["charges"]["equivalents"] = inferred(quantity(3.0, "equiv"), ev_general, "3.0 mmol n-BuLi divided by the explicitly reported 1.00 mmol nitrone scale.", raw="3.0 equiv calculated from general-procedure charges")
    reference_mid = byid[nitrone_id]["molecule_id"]
    byid["nBuLi"]["equivalence_reference"] = source(reference_mid, reference_mid, ev_general, method="inherited_general_condition")
    set_charge(byid["THF_main"], "volume", "10 mL (general procedure)", quantity(10, "mL"), ev_general, method="inherited_general_condition")
    set_charge(byid[nitrone_id], "mass", nitrone_mass, quantity(float(nitrone_mass.split()[0]), "mg"), ev_specific)
    set_charge(byid[nitrone_id], "amount", "1.00 mmol", quantity(1.00, "mmol"), ev_specific)
    set_charge(byid["THF_nitrone"], "volume", "2.0 mL (general procedure)", quantity(2.0, "mL"), ev_general, method="inherited_general_condition")
    set_charge(byid["water"], "volume", "15 mL (general procedure)", quantity(15, "mL"), ev_general, method="inherited_general_condition")
    set_charge(byid["Et2O_extract"], "volume", "3 × 10 mL (general procedure)", {"number": 10, "unit": "mL", "count": 3, "raw_text": "3 × 10 mL"}, ev_general, method="inherited_general_condition")


def build_synthesis_7():
    ev = [E[7], ET[7], E[8], ET[8]]
    r = make_reaction(
        "synthesis_7", ev, "Experimental Section, page 8; general procedure begins page 7",
        "Synthesis of 4-methoxy-2-methyl-3-(p-tolyl)-3,6-dihydro-2H-1,2-oxazine (7)",
        "lithiated methoxyallene addition to nitrone followed by spontaneous 6-endo-trig cyclization",
        ["substrate_scope"],
        ["The intermediate allenyl hydroxylamine is not assigned a compound number or isolated amount.",
         "The source calls the procedure general but does not state the temperature of the 16 h cyclization period.",
         "Product 7 is labelled racemic in Scheme 5, but no numeric enantiomer ratio or measurement method is reported.",
         "Complete structures and atom mapping remain unencoded."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1_general_A_to_B"]
    r["participants"] = [
        participant("methoxyallene", "m_1", "reactant", ev),
        participant("nBuLi", "m_nBuLi", "reagent", ev),
        participant("THF_main", "m_THF", "solvent", ev),
        participant("nitrone7", "m_nitrone7", "reactant", [E[8], ET[8]]),
        participant("THF_nitrone", "m_THF", "solvent", ev),
        participant("water", "m_H2O", "workup_agent", ev, stage="workup"),
        participant("Et2O_extract", "m_Et2O", "solvent", ev, stage="workup"),
        participant("Na2SO4", "m_Na2SO4", "workup_agent", ev, stage="workup"),
        participant("silica", "m_silica", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("hexanes", "m_hexanes", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("EtOAc", "m_EtOAc", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("product7", "m_7", "product", [E[8], ET[8]], product_role="isolated target 1,2-oxazine"),
    ]
    add_general_charges(r["participants"], "methoxyallene", "nitrone7", "149 mg", [E[8], ET[8]])
    r["reaction_scale"] = source("N-methyl-C-(p-tolyl)nitrone, 1.00 mmol", {"reference_molecule_id": "m_nitrone7", "number": 1.00, "unit": "mmol"}, [E[8], ET[8]])
    r["solvent_conditions"]["composition"] = source("freshly dried THF", ["THF"], [E[7], ET[7]], method="inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = inferred(quantity(12.0, "mL"), [E[7], ET[7]], "Sum of explicitly reported 10 mL main THF and 2.0 mL nitrone-solution THF; not an author-reported measured total.", raw="10 mL + 2.0 mL THF")
    r["solvent_conditions"]["concentration"] = inferred(quantity(1.00 / 0.012, "mM"), [E[7], ET[7], E[8], ET[8]], "Nominal nitrone concentration from 1.00 mmol divided by 12.0 mL explicitly charged THF; mixing volume is not independently measured.", raw="1.00 mmol / (10 + 2.0 mL THF)")
    r["solvent_conditions"]["concentration_basis"] = inferred("nitrone 7 precursor", [E[7], ET[7], E[8], ET[8]], "Denominator is the explicitly reported nitrone scale.")
    r["steps"] = [
        step("lithiation", [E[7], ET[7]], ["methoxyallene", "nBuLi", "THF_main"], "treat methoxyallene in freshly dried THF with n-butyllithium under dry argon", temperature=quantity(-40, "°C"), time=quantity(5, "min"), method="inherited_general_condition"),
        step("cool_for_nitrone", [E[7], ET[7]], ["methoxyallene", "nBuLi", "THF_main"], "cool the resulting solution from -40 °C to -78 °C; transition time is not reported", temperature=quantity(-78, "°C"), method="inherited_general_condition"),
        step("nitrone_addition", [E[7], ET[7]], ["nitrone7", "THF_nitrone", "methoxyallene", "nBuLi"], "add the nitrone solution dropwise", temperature=quantity(-78, "°C"), method="inherited_general_condition"),
        step("low_temperature_hold", [E[8], ET[8]], ["nitrone7", "methoxyallene"], "treat N-methyl-C-(p-tolyl)nitrone with lithiated methoxyallene", temperature=quantity(-78, "°C"), time=quantity(2, "h")),
        step("aqueous_quench", [E[7], ET[7]], ["water"], "add water", stage="workup", temperature=quantity(-78, "°C"), method="inherited_general_condition"),
        step("warm_to_room_temperature", [E[7], ET[7]], ["water"], "warm the mixture from -78 °C to room temperature; ramp time is not reported", stage="workup", temperature="room temperature", method="inherited_general_condition"),
        step("extract", [E[7], ET[7]], ["Et2O_extract", "water"], "extract with diethyl ether three times", stage="workup", method="inherited_general_condition"),
        step("dry", [E[7], ET[7]], ["Na2SO4"], "dry the combined organic layers with sodium sulfate", stage="workup", method="inherited_general_condition"),
        step("cyclization", [E[8], ET[8]], ["product7"], "allow the intermediate allenyl hydroxylamine to cyclize to completion", time=quantity(16, "h")),
        step("evaporate", [E[7], ET[7]], ["product7"], "remove solvents under reduced pressure", stage="workup", method="inherited_general_condition"),
        step("column", [E[8], ET[8]], ["silica", "hexanes", "EtOAc", "product7"], "purify by silica-gel column chromatography with hexanes/EtOAc 4:1", stage="purification"),
    ]
    r["procedure"]["original_text"] = source("Following the general procedure ... treated with lithiated methoxyallene (2 h) ... cyclization ... complete (16 h) ... column chromatography (SiO2, hexanes/EtOAc, 4:1)", "Specific preparation of 7 with inherited general-procedure operations and stated 2 h/16 h durations.", [E[8], ET[8]])
    r["procedure"]["general_procedure_reference"] = source("Following the general procedure", "General procedure for synthesis of allenyl hydroxylamines and 1,2-oxazines on pages 7-8", [E[7], ET[7], E[8], ET[8]])
    r["procedure"]["addition_order"] = source("alkoxyallene + n-BuLi, then nitrone solution, then water", ["alkoxyallene", "n-BuLi", "nitrone", "water"], [E[7], ET[7]], method="inherited_general_condition")
    r["outcomes"] = [
        outcome("isolated_mass", "product_amount", ["m_7"], [E[8], ET[8]], "160 mg", quantity(160, "mg"), product_role="isolated target", measurement_basis="after silica-gel column chromatography", scope="product 7"),
        outcome("isolated_yield", "yield", ["m_7"], [E[7], ET[7], E[8], ET[8]], "73%", quantity(73, "%"), product_role="isolated target", yield_type="isolated yield; general information states analytically pure unless otherwise noted", measurement_basis="weighed after silica-gel column chromatography", scope="overall from N-methyl-C-(p-tolyl)nitrone to 7"),
    ]
    r["notes"] = [source("pale yellow solid; m.p. 29-31 °C", {"appearance": "pale yellow solid", "melting_point": {"minimum": 29, "maximum": 31, "unit": "°C"}}, [E[8], ET[8]])]
    set_classification(r, ev, ANNULATION_CHANGES, "Stage-level labels cover the observed addition/workup/cyclization sequence; complete atoms, allene bond-order changes, and proton/deuteron sources are not mapped.")
    return r


def build_synthesis_9():
    ev = [E[7], ET[7], E[8], ET[8]]
    r = make_reaction(
        "synthesis_9", ev, "Experimental Section, page 8; general procedure begins page 7",
        "Synthesis of 4-benzyloxy-2,4a-dihydropiperidino[1,2-b][1,2]-oxazine (9)",
        "lithiated benzyloxyallene addition to cyclic nitrone followed by spontaneous cyclization",
        ["substrate_scope"],
        ["The intermediate allenyl hydroxylamine is not assigned a compound number or isolated amount.",
         "The source gives ca. -100 °C for the 4 h treatment but does not report a ramp time from the general-procedure lithiation temperature.",
         "Product 9 is labelled racemic in Scheme 5, but no numeric enantiomer ratio or measurement method is reported.",
         "Complete structures and atom mapping remain unencoded."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1_general_A_to_B"]
    r["participants"] = [
        participant("benzyloxyallene", "m_benzyloxyallene", "reactant", ev),
        participant("nBuLi", "m_nBuLi", "reagent", ev),
        participant("THF_main", "m_THF", "solvent", ev),
        participant("cyclic_nitrone", "m_tetrahydropiperidine_n_oxide", "reactant", [E[8], ET[8]]),
        participant("THF_nitrone", "m_THF", "solvent", ev),
        participant("water", "m_H2O", "workup_agent", ev, stage="workup"),
        participant("Et2O_extract", "m_Et2O", "solvent", ev, stage="workup"),
        participant("Na2SO4", "m_Na2SO4", "workup_agent", ev, stage="workup"),
        participant("silica", "m_silica", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("petroleum_ether", "m_petroleum_ether", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("EtOAc", "m_EtOAc", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("product9", "m_9", "product", [E[8], ET[8]], product_role="isolated target bicyclic 1,2-oxazine"),
    ]
    add_general_charges(r["participants"], "benzyloxyallene", "cyclic_nitrone", "99.1 mg", [E[8], ET[8]])
    r["reaction_scale"] = source("2,3,4,5-tetrahydropiperidine N-oxide, 1.00 mmol", {"reference_molecule_id": "m_tetrahydropiperidine_n_oxide", "number": 1.00, "unit": "mmol"}, [E[8], ET[8]])
    r["solvent_conditions"]["composition"] = source("freshly dried THF; ca. -100 °C Et2O/CO2 bath", {"reaction_solvent": "THF", "cooling_bath": "Et2O/CO2"}, [E[7], ET[7], E[8], ET[8]])
    r["solvent_conditions"]["total_volume"] = inferred(quantity(12.0, "mL"), [E[7], ET[7]], "Sum of 10 mL main THF and 2.0 mL nitrone-solution THF inherited from the general procedure; the cooling-bath ether is not included.", raw="10 mL + 2.0 mL THF")
    r["solvent_conditions"]["concentration"] = inferred(quantity(1.00 / 0.012, "mM"), [E[7], ET[7], E[8], ET[8]], "Nominal nitrone concentration from 1.00 mmol divided by 12.0 mL charged THF; not an author-reported concentration.", raw="1.00 mmol / (10 + 2.0 mL THF)")
    r["solvent_conditions"]["concentration_basis"] = inferred("2,3,4,5-tetrahydropiperidine N-oxide", [E[7], ET[7], E[8], ET[8]], "Denominator is the explicitly reported cyclic nitrone scale.")
    r["special_conditions"]["reactor"] = source("Et2O/CO2 bath", "diethyl ether/carbon dioxide cooling bath", [E[8], ET[8]])
    r["steps"] = [
        step("lithiation", [E[7], ET[7]], ["benzyloxyallene", "nBuLi", "THF_main"], "treat benzyloxyallene in freshly dried THF with n-butyllithium under dry argon", temperature=quantity(-40, "°C"), time=quantity(5, "min"), method="inherited_general_condition"),
        step("cool_for_nitrone", [E[7], ET[7], E[8], ET[8]], ["benzyloxyallene", "nBuLi", "THF_main"], "cool the lithiated benzyloxyallene from the general-procedure lithiation temperature to ca. -100 °C; ramp time is not reported", temperature=quantity(-100, "°C", raw="ca. -100 °C", qualifier="approx")),
        step("nitrone_addition", [E[7], ET[7], E[8], ET[8]], ["cyclic_nitrone", "THF_nitrone", "benzyloxyallene", "nBuLi"], "add the cyclic nitrone solution to lithiated benzyloxyallene", temperature=quantity(-100, "°C", raw="ca. -100 °C", qualifier="approx")),
        step("low_temperature_hold", [E[8], ET[8]], ["cyclic_nitrone", "benzyloxyallene"], "treat 2,3,4,5-tetrahydropiperidine N-oxide with lithiated benzyloxyallene in an Et2O/CO2 bath", temperature=quantity(-100, "°C", raw="ca. -100 °C", qualifier="approx"), time=quantity(4, "h")),
        step("aqueous_quench", [E[7], ET[7]], ["water"], "add water at -78 °C as stated by the referenced general procedure; the specific procedure does not restate the quench temperature", stage="workup", temperature=quantity(-78, "°C"), method="inherited_general_condition"),
        step("warm_to_room_temperature", [E[7], ET[7]], ["water"], "warm the mixture from -78 °C to room temperature; ramp time is not reported", stage="workup", temperature="room temperature", method="inherited_general_condition"),
        step("extract", [E[7], ET[7]], ["Et2O_extract", "water"], "extract with diethyl ether three times", stage="workup", method="inherited_general_condition"),
        step("dry", [E[7], ET[7]], ["Na2SO4"], "dry the combined organic layers with sodium sulfate", stage="workup", method="inherited_general_condition"),
        step("cyclization", [E[8], ET[8]], ["product9"], "allow the resulting allenyl hydroxylamine to cyclize completely", time=quantity(2, "d")),
        step("evaporate", [E[7], ET[7]], ["product9"], "remove solvents under reduced pressure", stage="workup", method="inherited_general_condition"),
        step("column", [E[8], ET[8]], ["silica", "petroleum_ether", "EtOAc", "product9"], "purify by silica-gel column chromatography with petroleum ether/EtOAc 5:1", stage="purification"),
    ]
    r["procedure"]["original_text"] = source("Following the general procedure ... lithiated benzyloxyallene at ca. -100 °C (Et2O/CO2 bath) for 4 h ... required 2 days ... column chromatography (SiO2, petroleum ether/EtOAc, 5:1)", "Specific preparation of 9 with low-temperature and purification modifications.", [E[8], ET[8]])
    r["procedure"]["general_procedure_reference"] = source("Following the general procedure", "General procedure for synthesis of allenyl hydroxylamines and 1,2-oxazines on pages 7-8", [E[7], ET[7], E[8], ET[8]])
    r["procedure"]["modifications"] = source("ca. -100 °C (Et2O/CO2 bath) for 4 h; 2 days for complete cyclization; petroleum ether/EtOAc 5:1", {"addition_hold_temperature": "ca. -100 °C", "addition_hold_time": "4 h", "cyclization_time": "2 d", "column_eluent": "petroleum ether/EtOAc 5:1"}, [E[8], ET[8]])
    r["procedure"]["addition_order"] = source("benzyloxyallene + n-BuLi, then cyclic nitrone, then water", ["benzyloxyallene", "n-BuLi", "cyclic nitrone", "water"], [E[7], ET[7], E[8], ET[8]])
    r["outcomes"] = [
        outcome("isolated_mass", "product_amount", ["m_9"], [E[8], ET[8]], "106 mg", quantity(106, "mg"), product_role="isolated target", measurement_basis="after silica-gel column chromatography", scope="product 9"),
        outcome("isolated_yield", "yield", ["m_9"], [E[7], ET[7], E[8], ET[8]], "44%", quantity(44, "%"), product_role="isolated target", yield_type="isolated yield; general information states analytically pure unless otherwise noted", measurement_basis="weighed after silica-gel column chromatography", scope="overall from cyclic nitrone to 9"),
    ]
    r["notes"] = [source("colorless solid; m.p. 100-103 °C", {"appearance": "colorless solid", "melting_point": {"minimum": 100, "maximum": 103, "unit": "°C"}}, [E[8], ET[8]])]
    set_classification(r, ev, ANNULATION_CHANGES, "Stage-level labels cover addition/workup/cyclization; complete atom mapping, ring stereochemistry, allene bond-order changes, and exact proton source are unresolved.")
    return r


def build_4d1_update():
    ev = [E[2], ET[2], E[7], ET[7], E[8], ET[8]]
    r = make_reaction(
        "scheme3_prepare_3d1", ev, "Scheme 3 and Experimental Section, synthesis of 3-d1 and 4-d1",
        "Synthesis of O-deuterated allenyl hydroxylamine 3-d1 and C-5-deuterated 1,2-oxazine 4-d1",
        "lithiated methoxyallene addition to nitrone, D2O quench, and spontaneous 6-endo-trig cyclization",
        ["mechanistic_study"],
        ["The specific procedure reports 4.00 mmol of 'lithiated methoxyallene' but does not separately state its n-BuLi charge; the general procedure uses 3.50 mmol alkoxyallene and 3.0 mmol n-BuLi and is not silently substituted for the specific value.",
         "Intermediate 3-d1 is generated in situ and not assigned an isolated mass or yield.",
         "The source reports >95% O-deuteration of 3-d1 via subsequent product evidence and >95% deuteration at C-5 of 4-d1; these locations are retained separately.",
         "Complete stereochemical structures and atom mapping remain unencoded."])
    r["general_scheme_ids"] = [f"{DOC}:scheme1_general_A_to_B"]
    r["participants"] = [
        participant("lithiated_methoxyallene", "m_lithiated_methoxyallene", "reactant", [E[8], ET[8]]),
        participant("nitrone2", "m_2", "reactant", ev),
        participant("THF", "m_THF", "solvent", [E[8], ET[8]]),
        participant("D2O_quench", "m_D2O", "workup_agent", [E[8], ET[8]], stage="workup"),
        participant("D2O_wash", "m_D2O", "workup_agent", [E[8], ET[8]], stage="workup"),
        participant("Na2SO4", "m_Na2SO4", "workup_agent", [E[8], ET[8]], stage="workup"),
        participant("silica_TLC", "m_silica", "purification_agent", [E[8], ET[8]], stage="characterization"),
        participant("p_anisaldehyde", "m_p_anisaldehyde", "reagent", [E[8], ET[8]], stage="characterization"),
        participant("silica_column", "m_silica", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("petroleum_ether", "m_petroleum_ether", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("EtOAc", "m_EtOAc", "purification_agent", [E[8], ET[8]], stage="purification"),
        participant("intermediate_3d1", "m_3d1", "product", ev, product_role="in situ O-deuterated allenyl hydroxylamine intermediate"),
        participant("product_4d1", "m_4d1", "product", ev, product_role="isolated C-5-deuterated 1,2-oxazine target"),
    ]
    byid = {p["participant_id"]: p for p in r["participants"]}
    set_charge(byid["lithiated_methoxyallene"], "amount", "4.00 mmol", quantity(4.00, "mmol"), [E[8], ET[8]])
    set_charge(byid["nitrone2"], "mass", "235 mg", quantity(235, "mg"), [E[8], ET[8]])
    set_charge(byid["nitrone2"], "amount", "1.00 mmol", quantity(1.00, "mmol"), [E[8], ET[8]])
    set_charge(byid["THF"], "volume", "12 mL", quantity(12, "mL"), [E[8], ET[8]])
    set_charge(byid["D2O_quench"], "volume", "3.0 mL", quantity(3.0, "mL"), [E[8], ET[8]])
    set_charge(byid["D2O_wash"], "volume", "2 × 2 mL", {"number": 2, "unit": "mL", "count": 2, "raw_text": "2 × 2 mL"}, [E[8], ET[8]])
    set_charge(byid["Na2SO4"], "mass", "1.5 g", quantity(1.5, "g"), [E[8], ET[8]])
    r["reaction_scale"] = source("D-glyceraldehyde-derived nitrone, 1.00 mmol", {"reference_molecule_id": "m_2", "number": 1.00, "unit": "mmol"}, [E[8], ET[8]])
    r["solvent_conditions"]["composition"] = source("freshly dried THF", ["THF"], [E[8], ET[8]])
    r["solvent_conditions"]["total_volume"] = source("12 mL", quantity(12, "mL"), [E[8], ET[8]])
    r["solvent_conditions"]["concentration"] = inferred(quantity(1.00 / 0.012, "mM"), [E[8], ET[8]], "Nominal nitrone concentration from 1.00 mmol divided by 12 mL THF; not an author-reported concentration.", raw="1.00 mmol / 12 mL THF")
    r["solvent_conditions"]["concentration_basis"] = inferred("D-glyceraldehyde-derived nitrone 2", [E[8], ET[8]], "Denominator is the explicitly reported nitrone scale.")
    r["steps"] = [
        step("low_temperature_reaction", [E[8], ET[8]], ["lithiated_methoxyallene", "nitrone2", "THF"], "treat lithiated methoxyallene with D-glyceraldehyde-derived nitrone in freshly dried THF", temperature=quantity(-78, "°C"), time=quantity(2, "h")),
        step("D2O_quench", [E[8], ET[8]], ["D2O_quench", "intermediate_3d1"], "quench the resulting mixture with D2O", stage="workup", temperature=quantity(-78, "°C")),
        step("warm_to_room_temperature", [E[8], ET[8]], ["D2O_quench", "intermediate_3d1"], "allow the mixture to warm from -78 °C to room temperature; ramp time is not reported", stage="workup", temperature="room temperature"),
        step("separate_and_D2O_wash", [E[8], ET[8]], ["D2O_wash", "intermediate_3d1"], "separate the layers and wash the organic layer twice with D2O", stage="workup"),
        step("dry_and_cyclize", [E[8], ET[8]], ["Na2SO4", "intermediate_3d1", "product_4d1", "silica_TLC", "p_anisaldehyde"], "add solid Na2SO4 and monitor cyclization by TLC (petroleum ether/EtOAc 7:1; p-anisaldehyde stain)", temperature="room temperature", time=quantity(24, "h")),
        step("filter", [E[8], ET[8]], ["Na2SO4", "product_4d1"], "filter the mixture", stage="workup"),
        step("evaporate", [E[8], ET[8]], ["product_4d1"], "remove solvent in vacuo", stage="workup"),
        step("column", [E[8], ET[8]], ["silica_column", "petroleum_ether", "EtOAc", "product_4d1"], "purify by silica-gel column chromatography with petroleum ether/EtOAc 6:1", stage="purification"),
    ]
    r["procedure"]["original_text"] = source("lithiated methoxyallene (4.00 mmol) ... nitrone (235 mg, 1.00 mmol) in THF (12 mL) at -78 °C for 2 h ... D2O (3.0 mL) ... D2O (2 × 2 mL) ... Na2SO4 (1.5 g) ... after 24 h ... SiO2, petroleum ether/EtOAc 6:1", "Specific preparation of in situ 3-d1 and isolated 4-d1.", [E[8], ET[8]])
    r["procedure"]["general_procedure_reference"] = source("Following the general procedure", "General procedure on pages 7-8, with all printed specific quantities and modifications retained", [E[7], ET[7], E[8], ET[8]])
    r["procedure"]["modifications"] = source("D2O quench and two D2O washes; 24 h cyclization; petroleum ether/EtOAc 6:1", {"quench": "D2O", "wash": "D2O, 2 × 2 mL", "cyclization_time": "24 h", "column_eluent": "petroleum ether/EtOAc 6:1"}, [E[8], ET[8]])
    r["procedure"]["addition_order"] = source("lithiated methoxyallene + nitrone, D2O quench, D2O washes, Na2SO4, filtration, evaporation, column", ["lithiated methoxyallene", "nitrone", "D2O quench", "D2O wash", "Na2SO4", "filtration", "evaporation", "column"], [E[8], ET[8]])
    r["outcomes"] = [
        outcome("intermediate_O_deuteration", "selectivity", ["m_3d1"], [E[2], ET[2], E[8], ET[8]], ">95% O-deuteration", quantity(95, "%", raw=">95%", qualifier="gt"), product_role="in situ intermediate", measurement_basis="inferred by the authors from subsequent cyclization to 4-d1 and its 13C NMR C-D coupling", scope="O-deuteration of intermediate 3-d1"),
        outcome("isolated_mass", "product_amount", ["m_4d1"], [E[8], ET[8]], "191 mg", quantity(191, "mg"), product_role="isolated target", measurement_basis="after silica-gel column chromatography", scope="product 4-d1"),
        outcome("isolated_yield", "yield", ["m_4d1"], [E[7], ET[7], E[8], ET[8]], "63%", quantity(63, "%"), product_role="isolated target", yield_type="isolated yield; general information states analytically pure unless otherwise noted", measurement_basis="weighed after silica-gel column chromatography", scope="overall from nitrone 2 to 4-d1"),
        outcome("C5_deuteration", "selectivity", ["m_4d1"], [E[8], ET[8]], "degree of deuteration at C-5 >95%", quantity(95, "%", raw=">95%", qualifier="gt"), product_role="isolated target", measurement_basis="13C NMR C-5 signal at 92.5 ppm with JC-D = 24.1 Hz", scope="C-5 deuterium incorporation in 4-d1"),
    ]
    r["notes"] = [
        source("colorless oil; [alpha]D20 = +33.1 (c = 0.26, CHCl3)", {"appearance": "colorless oil", "optical_rotation": {"value": 33.1, "temperature": 20, "concentration": 0.26, "solvent": "CHCl3"}}, [E[8], ET[8]]),
        st.fact(
            {
                "general_procedure": "alkoxyallene 3.50 mmol and n-BuLi 3.0 mmol",
                "specific_4d1_procedure": "lithiated methoxyallene (4.00 mmol)",
            },
            None,
            [E[7], ET[7], E[8], ET[8]],
            method="original_extraction",
            status="conflict",
            alternatives=[
                source("alkoxyallene 3.50 mmol and n-BuLi 3.0 mmol", {"alkoxyallene": quantity(3.50, "mmol"), "nBuLi": quantity(3.0, "mmol")}, [E[7], ET[7]]),
                source("lithiated methoxyallene (4.00 mmol)", {"lithiated_methoxyallene": quantity(4.00, "mmol")}, [E[8], ET[8]]),
            ],
        ),
    ]
    set_classification(r, ev, ANNULATION_CHANGES, "Stage-level C-C/O-H/C-O/C-H labels include the isolated hydroxylamine stage and cyclization; isotope-specific and complete atom mapping are unresolved.")
    return r


def build_scheme9():
    ev = [E[7], ET[7]]
    r = st.new_reaction(f"{DOC}:scheme9_zwitterionic_summary", DOC, ev)
    r["record_type"] = "proposed_mechanism"
    r["experiment_categories"] = not_applicable()
    r["location"] = source("Scheme 9", "Scheme 9", ev)
    r["reaction_name"] = source("Formation of 1,2-oxazines B from A via zwitterion ZW", "Formation of 1,2-oxazines B from A via zwitterion ZW", ev)
    r["transformation_type"] = inferred("proposed proton-transfer/cyclization and fragmentation-recombination network", ev, "Mechanistic class interpreted from Scheme 9 and the adjacent conclusions.")
    r["general_scheme_ids"] = [f"{DOC}:scheme1_general_A_to_B"]
    r["participants"] = [
        participant("A", "m_A", "reactant", ev),
        participant("ZW", "m_ZW", "unknown", ev),
        participant("B", "m_B", "product", ev, product_role="proposed 1,2-oxazine product"),
        participant("G", "m_G", "product", ev, product_role="proposed diene fragmentation product"),
        participant("H", "m_H", "product", ev, product_role="proposed nitroso fragmentation product"),
    ]
    r["notes"] = [
        source("A -> ZW -> B; ZW -> G + H; G + H ?-> B", {"branches": ["A to ZW", "ZW to B", "ZW to G + H", "G + H uncertain recombination to B"]}, ev),
        source("Experiments should be designed to prove or disprove this alternative pathway to 1,2-oxazine B.", "The G/H recombination branch is proposed and explicitly unproven.", ev),
    ]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Scheme 9 is a generic proposed mechanism, not a laboratory execution.",
                                  "R1/R2/R3 are not instantiated and the question-marked G/H recombination branch is explicitly unproven.",
                                  "Complete atom mapping and fragmentation-bond accounting are unresolved."]}
    r["completeness"] = "partial"
    set_classification(r, ev, ANNULATION_CHANGES[3:5], "Only the proposed ZW-to-B cyclization branch is element-level classified; proton-transfer, fragmentation, and question-marked recombination branches are not completely mapped.")
    return r


def update_scheme8():
    r = copy.deepcopy(RXNS[f"{DOC}:scheme8_protonated_network"])
    append_evidence(r["evidence_ids"], E[7], ET[7])
    r["location"] = source("Scheme 8 and discussion continuing on page 7", "Scheme 8 and pages 6-7 discussion", [E[6], E[7], ET[7]])
    if "Page 7 continues the calculated protonated-species pathways and remains computational, not experimental." not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append("Page 7 continues the calculated protonated-species pathways and remains computational, not experimental.")
    return r


def update_kinetic_record(rid):
    r = copy.deepcopy(RXNS[rid])
    ev = [E[8], ET[8]]
    append_evidence(r["evidence_ids"], *ev)
    substrate = next(p for p in r["participants"] if p["participant_id"] == "substrate")
    is_d = substrate["molecule_id"] == "m_3d1"
    reserved = {"sample_quench", "cold_Et2O", "sample_wash", "drying_agent"}
    r["participants"] = [p for p in r["participants"] if p["participant_id"] not in reserved]
    r["participants"].append(participant("sample_quench", "m_D2O" if is_d else "m_H2O", "workup_agent", ev, stage="workup"))
    if not is_d:
        r["participants"].append(participant("cold_Et2O", "m_Et2O", "solvent", ev, stage="workup"))
    r["participants"].append(participant("sample_wash", "m_D2O" if is_d else "m_H2O", "workup_agent", ev, stage="workup"))
    drying = participant("drying_agent", "m_Na2SO4", "workup_agent", ev, stage="workup")
    set_charge(drying, "mass", "1.5 g", quantity(1.5, "g"), ev)
    r["participants"].append(drying)
    solvent = next(p for p in r["participants"] if p["participant_id"] == "solvent")
    set_charge(solvent, "volume", "0.8 mL", quantity(0.8, "mL"), ev)
    r["solvent_conditions"]["total_volume"] = source("0.8 mL corresponding deuterated solvent", quantity(0.8, "mL"), ev)
    remove_ids = {"sample_quench", "sample_warm_to_zero", "sample_separate_wash_dry", "sample_aliquot", "sample_evaporate", "sample_redissolve"}
    old_steps = [s for s in r["steps"] if s["step_id"] not in remove_ids]
    prep_steps = [
        step("sample_quench", ev, ["substrate", "sample_quench"], f"quench the hydroxylamine preparation with {'D2O' if is_d else 'water'}", stage="workup", temperature=quantity(-78, "°C")),
        step("sample_warm_to_zero", ev, ["substrate", "sample_quench"], "allow the quenched mixture to warm from -78 °C to 0 °C; ramp time is not reported", stage="workup", temperature=quantity(0, "°C")),
        step("sample_separate_wash_dry", ev, ["substrate", "sample_wash", "drying_agent"] + ([] if is_d else ["cold_Et2O"]), "separate the organic layer" + ("" if is_d else ", adding cold diethyl ether to achieve phase separation") + f"; wash rapidly with {'D2O' if is_d else 'water'} and dry with Na2SO4", stage="workup", temperature=quantity(4, "°C"), time=quantity(15, "min")),
        step("sample_aliquot", ev, ["substrate"], f"transfer {'0.2-0.3 mL' if is_d else '0.6-1.0 mL'} of the dried organic solution by syringe into a flame-dried flask", stage="workup"),
        step("sample_evaporate", ev, ["substrate"], "remove solvent in vacuo to give a thick oil", stage="workup", temperature="room temperature"),
        step("sample_redissolve", ev, ["substrate", "solvent"], "dissolve the residue in 0.8 mL of the corresponding deuterated solvent, optionally with the tabulated additive", stage="reaction"),
    ]
    r["steps"] = prep_steps + old_steps
    r["procedure"]["general_procedure_reference"] = source("Preparation of allenyl hydroxylamine samples for kinetic measurements", "Page 8 kinetic-sample preparation protocol", ev)
    r["procedure"]["pretreatment"] = source("quench at -78 °C; warm to 0 °C; separate/wash; Na2SO4 1.5 g at 4 °C for ca. 15 min; aliquot; evaporate at room temperature; dissolve in 0.8 mL deuterated solvent", "Ordered sample preparation before NMR monitoring", ev)
    r["notes"] = [n for n in r["notes"] if not (n.get("raw_value") or "").startswith("Kinetic sample aliquot")]
    r["notes"].append(source(
        f"Kinetic sample aliquot: {'0.2-0.3 mL for 3-d1' if is_d else '0.6-1.0 mL for 3'}; remaining solution stored with drying agent in a dry-ice Dewar until required",
        {"sample_aliquot": {"minimum": 0.2 if is_d else 0.6, "maximum": 0.3 if is_d else 1.0, "unit": "mL"}, "remaining_solution_storage": "with drying agent in dry-ice Dewar"}, ev))
    missing_note = "The page-8 sample-preparation protocol does not report exact water/D2O wash volumes, cold-Et2O volume, substrate amount in each aliquot, or the dry-ice storage duration."
    if missing_note not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append(missing_note)
    return r


def main():
    template = json.loads((ROOT / "runs/current/raw_model/init_pages_0007_0009.json").read_text())
    molecules = []

    # New generic/mechanistic and experimental identities.
    molecules.append(generic_zw())
    molecules.append(smiles_molecule("m_nitrone7", "N-methyl-C-(p-tolyl)nitrone", "N-methyl-C-(p-tolyl)nitrone", "Cc1ccc(C=[N+]([O-])C)cc1", [E[8], ET[8]], "Candidate connectivity from the explicit source name; nitrone C=N geometry is not assigned."))
    molecules.append(smiles_molecule("m_benzyloxyallene", "benzyloxyallene", "benzyloxyallene", "C=C=COCc1ccccc1", [E[8], ET[8]], "Candidate connectivity from the explicit source name; terminal-allene axial stereogenicity is not asserted."))
    molecules.append(smiles_molecule("m_tetrahydropiperidine_n_oxide", "2,3,4,5-tetrahydropiperidine N-oxide", "2,3,4,5-tetrahydropiperidine N-oxide", "[O-][N+]1=CCCCC1", [E[8], ET[8]], "Candidate cyclic nitrone connectivity from the explicit source name; ring/double-bond geometry is not asserted."))
    molecules.append(unresolved_molecule("m_lithiated_methoxyallene", "lithiated methoxyallene", "lithiated methoxyallene", [E[8], ET[8]], chirality="unresolved", geometry="unresolved", representation="not_encoded", axial="Allene-containing organolithium reagent; axial configuration is not reported.", identity_note="Specific procedure reports this charged reagent without a complete ionic/counterion structure."))
    molecules.append(smiles_molecule("m_Na2SO4", "Na2SO4", "sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]", [E[7], ET[7], E[8], ET[8]], "Drying agent; hydration state is not stated."))
    molecules.append(smiles_molecule("m_EtOAc", "EtOAc", "ethyl acetate", "CCOC(C)=O", [E[8], ET[8]], "Chromatography eluent component."))
    molecules.append(smiles_molecule("m_p_anisaldehyde", "p-anisaldehyde", "p-anisaldehyde", "COc1ccc(C=O)cc1", [E[8], ET[8]], "TLC visualization stain, not a reaction-stage reagent."))
    molecules.append(unresolved_molecule("m_silica", "SiO2", "silica gel", [E[7], ET[7], E[8], ET[8]], chirality="not_applicable", geometry="not_applicable", representation="not_applicable", identity_note="Stationary phase; no discrete molecular SMILES asserted."))
    molecules.append(unresolved_molecule("m_hexanes", "hexanes", "hexanes", [E[8], ET[8]], chirality="mixture", geometry="unresolved", representation="not_encoded", mixture="Commercial hexanes mixture; composition is not specified.", identity_note="Named chromatography solvent mixture; no unique molecular identity inferred."))
    molecules.append(unresolved_molecule("m_petroleum_ether", "petroleum ether", "petroleum ether", [E[8], ET[8]], chirality="mixture", geometry="unresolved", representation="not_encoded", mixture="Petroleum ether mixture; boiling range/composition is not specified.", identity_note="Named chromatography solvent mixture; no unique molecular identity inferred."))

    # Existing molecule identities gain page-8 occurrences without collapsing stereo states.
    for mid, label in (("m_3d1", "3-d1"), ("m_4d1", "4-d1"), ("m_7", "7"), ("m_9", "9"),
                       ("m_1", "methoxyallene"), ("m_2", "D-glyceraldehyde-derived nitrone"),
                       ("m_nBuLi", "n-butyllithium"), ("m_THF", "THF"), ("m_H2O", "H2O"),
                       ("m_D2O", "D2O"), ("m_Et2O", "Et2O")):
        m = copy.deepcopy(MOLS[mid])
        add_occurrence(m, label, [E[8], ET[8]])
        molecules.append(m)

    # Replace detailed names with explicit experimental-section names while retaining earlier occurrences/stereo evidence.
    bymid = {m["molecule_id"]: m for m in molecules}
    bymid["m_7"]["names"]["original"] = source("4-methoxy-2-methyl-3-(p-tolyl)-3,6-dihydro-2H-1,2-oxazine (7)", "4-methoxy-2-methyl-3-(p-tolyl)-3,6-dihydro-2H-1,2-oxazine", [E[8], ET[8]])
    bymid["m_7"]["names"]["english"] = copy.deepcopy(bymid["m_7"]["names"]["original"])
    bymid["m_9"]["names"]["original"] = source("4-benzyloxy-2,4a-dihydropiperidino[1,2-b][1,2]-oxazine (9)", "4-benzyloxy-2,4a-dihydropiperidino[1,2-b][1,2]-oxazine", [E[8], ET[8]])
    bymid["m_9"]["names"]["english"] = copy.deepcopy(bymid["m_9"]["names"]["original"])

    # A submitted bundle must be internally closed. Retain the existing v3.1
    # registry entries needed by referenced participants; merge will update only
    # records whose IDs are present with changed content.
    present_mids = {m["molecule_id"] for m in molecules}
    molecules.extend(copy.deepcopy(m) for mid, m in MOLS.items() if mid not in present_mids)

    reactions = [build_4d1_update(), build_synthesis_7(), build_synthesis_9(), build_scheme9(), update_scheme8()]
    kinetic_ids = [f"{DOC}:table1_entry_{i:02d}" for i in range(1, 13)]
    kinetic_ids += [f"{DOC}:table2_entry_{i:02d}" for i in range(2, 8)]
    kinetic_ids += [f"{DOC}:pyridine_3d1_questionable"]
    reactions.extend(update_kinetic_record(rid) for rid in kinetic_ids)
    # General-scheme references are part of the closed submission dependency set.
    incoming_rids = {r["reaction_id"] for r in reactions}
    for dependency_id in (f"{DOC}:scheme1_general_A_to_B", f"{DOC}:scheme4_general_3_to_4_5"):
        if dependency_id not in incoming_rids:
            reactions.append(copy.deepcopy(RXNS[dependency_id]))

    template["molecules"] = molecules
    template["reactions"] = reactions
    template["source_context"] = [
        {"context_id": f"{DOC}:page7_conclusions", "document_id": DOC, "evidence_ids": [E[7], ET[7]], "kind": "conclusions_and_mechanistic_interpretation", "note": "Conclusion text summarizes isotope/additive results and proposes the zwitterionic pathway; it is not a new execution."},
        {"context_id": f"{DOC}:pages7_8_general_information", "document_id": DOC, "evidence_ids": [E[7], ET[7], E[8], ET[8]], "kind": "general_experimental_information", "note": "Reagent handling, THF drying, argon/flame-dried glassware, chromatography, analytical-purity yield convention, and the general hydroxylamine/oxazine procedure used as explicitly inherited evidence."},
        {"context_id": f"{DOC}:page8_characterization_4d1_7_9", "document_id": DOC, "evidence_ids": [E[8], ET[8]], "kind": "characterization", "note": "Optical rotation, NMR, IR, MS/HRMS, elemental analyses and melting points for 4-d1, 7 and 9; no separate synthesis records were created for spectral peak lists."},
        {"context_id": f"{DOC}:page8_kinetic_sample_preparation", "document_id": DOC, "evidence_ids": [E[8], ET[8]], "kind": "shared_experimental_procedure", "note": "Shared quench/wash/dry/aliquot/evaporate/redissolve protocol was linked to each existing kinetic experiment rather than counted as a separate reaction."},
        {"context_id": f"{DOC}:page8_quantum_calculation_method", "document_id": DOC, "evidence_ids": [E[8], ET[8]], "kind": "computational_method", "note": "TPSSTPSS/def2-TZVP+GD3BJ, Gaussian 09/16, PCM solvent, reaction-path/QST2/IRC description; not a laboratory reaction."},
        {"context_id": f"{DOC}:page9_back_matter", "document_id": DOC, "evidence_ids": [E[9], ET[9]], "kind": "back_matter_and_references", "note": "Acknowledgment continuation, keywords, references 1-28 and received date; cited reactions were not reclassified as experiments in this paper."},
    ]

    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {"page": 7, "state": "attempted", "regions": [
                {"region": "left column: Scheme 8 discussion and Conclusions", "finding": "computational/proposed mechanism discussion; continuation linked to the existing Scheme 8 record, not a new laboratory experiment", "evidence_ids": [E[7], ET[7]]},
                {"region": "right column: Scheme 9 and adjacent Conclusions", "finding": "generic proposed A-to-ZW-to-B cyclization plus ZW fragmentation to G/H and question-marked recombination; new proposed-mechanism record", "evidence_ids": [E[7], ET[7]]},
                {"region": "right column lower: Experimental Section general information and general procedure", "finding": "shared experimental conditions beginning the cross-page general procedure; inherited by the 4-d1, 7 and 9 preparations, not counted as an independent execution", "evidence_ids": [E[7], ET[7]]}],
             "record_ids": [f"{DOC}:scheme8_protonated_network", f"{DOC}:scheme9_zwitterionic_summary", f"{DOC}:scheme3_prepare_3d1", f"{DOC}:synthesis_7", f"{DOC}:synthesis_9"],
             "unresolved": ["Scheme 9 R1/R2/R3 identities and complete mapped graphs are not given.", "The question-marked G/H recombination to B is explicitly unproven.", "The general procedure continues onto page 8; its ramp time from -78 °C to room temperature and Na2SO4 amount are not reported."]},
            {"page": 8, "state": "attempted", "regions": [
                {"region": "top-left continuation of general procedure", "finding": "workup, drying, cyclization monitoring, solvent removal and column purification shared by specific preparations", "evidence_ids": [E[8], ET[8]]},
                {"region": "left column: synthesis of 3-d1 and 4-d1", "finding": "specific isotope-labeling preparation merged into the existing Scheme 3 record; 191 mg, 63%, C-5 deuteration >95%", "evidence_ids": [E[8], ET[8]]},
                {"region": "left column: synthesis of 7", "finding": "specific substrate-scope preparation with 2 h low-temperature treatment, 16 h cyclization, 160 mg and 73% isolated yield", "evidence_ids": [E[8], ET[8]]},
                {"region": "lower-left/right top: synthesis of 9", "finding": "specific substrate-scope preparation at ca. -100 °C for 4 h, 2 d cyclization, 106 mg and 44% isolated yield", "evidence_ids": [E[8], ET[8]]},
                {"region": "right column: preparation of kinetic samples", "finding": "shared sample-preparation procedure linked to the 12 Table 1 rows, six nonduplicate Table 2 additive rows, and prose-only pyridine isotope experiment; not a separate reaction", "evidence_ids": [E[8], ET[8]]},
                {"region": "analytical paragraphs for 4-d1, 7 and 9", "finding": "NMR/IR/MS/elemental analysis/melting-point characterization retained as source context, not separate syntheses", "evidence_ids": [E[8], ET[8]]},
                {"region": "right column lower: quantum chemical calculations", "finding": "computational method only; no new experimental reaction", "evidence_ids": [E[8], ET[8]]}],
             "record_ids": [f"{DOC}:scheme3_prepare_3d1", f"{DOC}:synthesis_7", f"{DOC}:synthesis_9"] + kinetic_ids,
             "unresolved": ["Specific 4-d1 text does not separately report the n-BuLi charge used to generate the stated 4.00 mmol lithiated methoxyallene.", "Exact temperature of the 16 h and 2 d post-workup cyclizations is not printed.", "Kinetic sample protocol omits exact water/D2O wash volume, cold-Et2O volume, substrate amount per aliquot and storage duration.", "Complete source-checked SMILES for 3-d1, 4-d1, 7 and 9 remain unencoded; stereo assignments are retained textually."]},
            {"page": 9, "state": "attempted", "regions": [
                {"region": "top-left acknowledgment continuation and keywords", "finding": "back matter; no reaction or result record", "evidence_ids": [E[9], ET[9]]},
                {"region": "two-column references 1-28", "finding": "cited prior literature only; no cited procedure was promoted to an experiment of this paper", "evidence_ids": [E[9], ET[9]]},
                {"region": "bottom-right received date", "finding": "bibliographic metadata; no reaction content", "evidence_ids": [E[9], ET[9]]}],
             "record_ids": [],
             "unresolved": ["No experimental arrows, tables or procedures occur on page 9; cited papers were not inspected or associated as source documents."]},
        ],
    }

    inventory = {
        "document_id": DOC,
        "page_range": [7, 9],
        "visual_review": "all three registered original page images viewed at original detail",
        "observations": {
            "proposed_mechanism": ["Scheme 9"],
            "specific_preparations": ["merged Scheme 3 / synthesis of 3-d1 and 4-d1", "synthesis of 7", "synthesis of 9"],
            "shared_procedures": ["general procedure pages 7-8", "kinetic sample preparation page 8"],
            "characterization_only": ["4-d1", "7", "9"],
            "nonreaction_page": [9],
            "new_reaction_records": 3,
            "updated_existing_reaction_records": 21,
        },
        "chemical_accuracy_verified": False,
    }

    st.write(ROOT / "runs/current/submitted.json", template)
    st.write(ROOT / "runs/current/raw_model/pages_0007_0009_extraction.json", inventory)
    st.write(ROOT / "runs/current/raw_model/stereo_checks_pages_0007_0009.json", STEREO_CHECKS)
    st.write(ROOT / "runs/current/coverage/pages_0007_0009.json", coverage)
    print(json.dumps({"bundle_molecules": len(molecules), "bundle_reactions": len(reactions),
                      "source_changed_reactions": 24,
                      "new_reactions": 3, "updated_reactions": 21,
                      "stereo_checks": len(STEREO_CHECKS)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
