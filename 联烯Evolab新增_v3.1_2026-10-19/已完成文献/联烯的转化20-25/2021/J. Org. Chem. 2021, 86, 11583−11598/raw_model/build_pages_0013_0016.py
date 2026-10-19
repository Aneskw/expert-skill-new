#!/usr/bin/env python3
"""Build the v3.1 incremental submission for source PDF pages 13--16.

Pages 13--16 add no new experiment identities. Page 13 completes four
Scheme 5 product procedures, the 7ha procedure, and three existing scale-up
records. Pages 14--16 contain crystallographic/bibliographic context only.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as sr  # noqa: E402


DOC = "doc_ee191f89470cf89e"
BATCH = "allene-stereo-v31-200-20261017"


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def evidence(page: int) -> list[str]:
    return [ev(page, "image"), ev(page, "text")]


def q(number, unit, qualifier=None, raw_text=None):
    value = {"number": number, "unit": unit}
    if qualifier is not None:
        value["qualifier"] = qualifier
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


def qraw(value):
    if not isinstance(value, dict):
        return str(value)
    if value.get("raw_text"):
        return value["raw_text"]
    qualifier = {"lt": "<", "le": "≤", "gt": ">", "ge": "≥", "approx": "~"}.get(value.get("qualifier"), "")
    return f"{qualifier}{value['number']:g} {value['unit']}"


def sf(raw=None, normalized=None, evs=(), method="original_extraction",
       status=None, derivation=None, alternatives=()):
    return sr.fact(raw=raw, normalized=normalized, evidence=evs, method=method,
                   status=status, derivation=derivation, alternatives=alternatives)


def unique(values):
    return list(dict.fromkeys(values))


def merge_evidence(fact, evs):
    if not isinstance(fact, dict) or "evidence_ids" not in fact:
        return
    fact["evidence_ids"] = unique(fact["evidence_ids"] + list(evs))


def mid(label):
    return f"{DOC}:mol:{label}"


template = json.loads((RUN / "raw_model/pages_0013_0016_init.json").read_text())
current = json.loads((RUN / "standard_records.json").read_text())
bundle = copy.deepcopy(template)
bundle["molecules"] = copy.deepcopy(current["molecules"])
bundle["reactions"] = []
bundle["source_context"] = []

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: copy.deepcopy(r) for r in current["reactions"]}


def get_reaction(suffix):
    return reactions[f"{DOC}:{suffix}"]


def include_reaction(reaction):
    for index, item in enumerate(bundle["reactions"]):
        if item["reaction_id"] == reaction["reaction_id"]:
            bundle["reactions"][index] = reaction
            return
    bundle["reactions"].append(reaction)


def participant(pid, molecule_id, role, evs, stage="reaction", mass=None,
                amount=None, volume=None, equiv=None, molpct=None, ref=None,
                product_role=None):
    p = sr.new_participant(pid, molecule_id, role, evs, stage=stage)
    if mass is not None:
        p["charges"]["mass"] = sf(qraw(mass), mass, evs)
    if amount is not None:
        p["charges"]["amount"] = sf(qraw(amount), amount, evs)
    if volume is not None:
        p["charges"]["volume"] = sf(qraw(volume), volume, evs)
    if equiv is not None:
        p["charges"]["equivalents"] = sf(qraw(equiv), equiv, evs)
    if molpct is not None:
        p["charges"]["mol_percent"] = sf(qraw(molpct), molpct, evs)
    if ref is not None:
        p["equivalence_reference"] = sf(ref, ref, evs)
    if product_role is not None:
        p["product_role"] = sf(product_role, product_role, evs)
    return p


def add_step(reaction, sid, evs, pids, operation, temperature=None, time=None,
             stage="reaction"):
    step = sr.new_step(sid, evs, stage=stage)
    step["participant_ids"] = list(pids)
    step["operation"] = sf(operation, operation, evs)
    if temperature is not None:
        raw = qraw(temperature) if isinstance(temperature, dict) else str(temperature)
        step["temperature"] = sf(raw, temperature, evs)
    if time is not None:
        step["time"] = sf(qraw(time), time, evs)
    reaction["steps"].append(step)


def add_outcome(reaction, oid, kind, target, evs, value, raw, basis=None,
                yield_type=None, product_role="desired product", scope=None):
    outcome = sr.new_outcome(oid, kind, [target], evs)
    outcome["value"] = sf(raw, value, evs)
    if basis is not None:
        outcome["measurement_basis"] = sf(basis, basis, evs)
    if yield_type is not None:
        outcome["yield_type"] = sf(yield_type, yield_type, evs)
    if product_role is not None:
        outcome["product_role"] = sf(product_role, product_role, evs)
    if scope is not None:
        outcome["scope"] = sf(scope, scope, evs)
    reaction["outcomes"].append(outcome)


def add_note(reaction, text, evs, method="original_extraction", derivation=None):
    for note in reaction["notes"]:
        if note.get("normalized_value") == text:
            merge_evidence(note, evs)
            return
    reaction["notes"].append(sf(
        text if method == "original_extraction" else None,
        text,
        evs,
        method=method,
        derivation=derivation,
    ))


def conflict_fact(alternatives):
    return sf(None, None, (), status="conflict", alternatives=[
        sf(raw, normalized, evs) for raw, normalized, evs in alternatives
    ])


def find_participant(reaction, pid):
    return next(p for p in reaction["participants"] if p["participant_id"] == pid)


def set_charge(reaction, pid, field, value, evs, raw=None):
    find_participant(reaction, pid)["charges"][field] = sf(raw or qraw(value), value, evs)


def replace_product_amount(reaction, target, mass, evs):
    reaction["outcomes"] = [o for o in reaction["outcomes"] if o["outcome_id"] != "product_amount"]
    add_outcome(reaction, "product_amount", "product_amount", target, evs,
                mass, qraw(mass), basis="reported isolated mass")


def enrich_yield(reaction, evs):
    for outcome in reaction["outcomes"]:
        if outcome["kind"] == "yield":
            merge_evidence(outcome["value"], evs)
            outcome["measurement_basis"] = sf("isolated yield", "isolated yield", evs)
            outcome["yield_type"] = sf("isolated", "isolated", evs)


def workup_participants(extract_volume, evs, ratio):
    return [
        participant("p_water", mid("H2O"), "workup_agent", evs, stage="workup", volume=q(10, "mL")),
        participant("p_extract_EtOAc", mid("EtOAc"), "workup_agent", evs, stage="workup",
                    volume=q(extract_volume, "mL", raw_text=f"3 × {extract_volume} mL")),
        participant("p_brine", mid("brine solution"), "workup_agent", evs, stage="workup", volume=q(20, "mL")),
        participant("p_drying", mid("Na2SO4"), "workup_agent", evs, stage="workup"),
        participant("p_column_EtOAc", mid("EtOAc"), "purification_agent", evs, stage="purification"),
        participant("p_column_hexane", mid("hexane"), "purification_agent", evs, stage="purification"),
    ]


def append_workup_steps(reaction, evs, extract_volume, ratio):
    add_step(reaction, "quench", evs, ["p_water"],
             "quench by addition of water (10 mL)", stage="workup")
    add_step(reaction, "extraction", evs, ["p_extract_EtOAc"],
             f"extract aqueous layer with ethyl acetate (3 × {extract_volume} mL)", stage="workup")
    add_step(reaction, "wash_dry", evs, ["p_brine", "p_drying"],
             "wash combined organic layer with brine (20 mL), dry over anhydrous Na2SO4, and concentrate",
             stage="workup")
    add_step(reaction, "purification", evs,
             ["p_column_EtOAc", "p_column_hexane", "p_product"],
             f"silica gel chromatography using ethyl acetate/hexane {ratio}",
             stage="purification")


def set_z_geometry(label, evs):
    molecule = molecules[mid(label)]
    stereo = molecule["stereochemistry"]
    stereo["double_bonds"] = sf(
        f"Experimental heading explicitly identifies (Z)-{label}.",
        [{"location": "2'-tosylimino C=N", "assignment": "Z"}],
        evs,
    )
    stereo["geometry_status"] = sf(
        None,
        "specified",
        evs,
        method="inference",
        derivation="Normalized from the explicit (Z) descriptor in the experimental heading.",
    )
    stereo["representation_status"] = sf(
        None,
        "not_encoded",
        evs,
        method="inference",
        derivation="The explicit Z assignment is retained textually; no complete source-verified isomeric SMILES was assigned.",
    )
    merge_evidence(stereo["chirality_status"], evs)


def set_7ha_geometry(evs):
    molecule = molecules[mid("7ha")]
    stereo = molecule["stereochemistry"]
    stereo["double_bonds"] = sf(
        "Experimental heading explicitly identifies (2Z,4E)-7ha.",
        [{"position": 2, "assignment": "Z"}, {"position": 4, "assignment": "E"}],
        evs,
    )
    stereo["geometry_status"] = sf(
        None,
        "specified",
        evs,
        method="inference",
        derivation="Normalized from the explicit (2Z,4E) descriptor in the experimental heading.",
    )
    stereo["representation_status"] = sf(
        None,
        "not_encoded",
        evs,
        method="inference",
        derivation="The explicit alkene assignments are retained textually; no complete source-verified isomeric SMILES was assigned.",
    )
    merge_evidence(stereo["chirality_status"], evs)


E11 = evidence(11)
E13 = evidence(13)
E14 = evidence(14)
E15 = evidence(15)
E16 = evidence(16)


# Complete the four Scheme 5 product procedures that begin on page 13.
scheme5_details = [
    ("6da", "1d", "5a", 75.3, 77.4, 96.0, 75, (75, 77), None),
    ("6bf", "1b", "5f", 60.0, 79.0, 83.0, 73, (66, 68), 0.04),
    ("6ga", "1g", "5a", 69.0, 77.4, 86.1, 71, (96, 98), None),
    ("6hb", "1h", "5b", 77.2, 84.5, 114.0, 84, (63, 65), 0.04),
]
scheme5_ids = []
for product, sub1, sub2, sub1_mass, sub2_mass, product_mass, yld, mp, printed_amount in scheme5_details:
    reaction = get_reaction(f"scheme5_{product}")
    reaction["evidence_ids"] = unique(reaction["evidence_ids"] + E11 + E13)
    reaction["location"] = sf(
        f"Scheme 5 and experimental section vi, {product} (page 13)",
        f"Scheme 5 and experimental section vi, {product}",
        E11 + E13,
    )
    set_z_geometry(product, E13)
    set_charge(reaction, "p_sub1", "mass", q(sub1_mass, "mg"), E13)
    set_charge(reaction, "p_sub2", "mass", q(sub2_mass, "mg"), E13)
    set_charge(reaction, "p_cat", "mass", q(12.2, "mg"), E13)
    if printed_amount is not None:
        amount_conflict = conflict_fact([
            ("Scheme 5/general procedure vi: sulfonamidoindole (0.20 mmol)", q(0.20, "mmol"), evidence(5) + E11),
            (f"Specific {product} paragraph: {sub1} (0.04 mmol)", q(printed_amount, "mmol"), E13),
        ])
        find_participant(reaction, "p_sub1")["charges"]["amount"] = amount_conflict
        reaction["reaction_scale"] = conflict_fact([
            ("Scheme 5/general procedure vi: sulfonamidoindole (0.20 mmol)",
             {"reference_molecule_id": mid(sub1), "amount": q(0.20, "mmol")}, evidence(5) + E11),
            (f"Specific {product} paragraph: {sub1} (0.04 mmol)",
             {"reference_molecule_id": mid(sub1), "amount": q(printed_amount, "mmol")}, E13),
        ])
    else:
        set_charge(reaction, "p_sub1", "amount", q(0.20, "mmol"), E13)
        reaction["reaction_scale"] = sf(
            f"{sub1} (0.20 mmol)",
            {"reference_molecule_id": mid(sub1), "amount": q(0.20, "mmol")},
            E13,
        )
    set_charge(reaction, "p_sub2", "amount", q(0.24, "mmol"), E13)
    set_charge(reaction, "p_cat", "amount", q(0.04, "mmol"), E13)
    set_charge(reaction, "p_cat", "mol_percent", q(20, "mol%"), E11)
    if any(p["participant_id"] == "p_solvent" for p in reaction["participants"]):
        find_participant(reaction, "p_solvent")["participant_id"] = "p_solvent_charge"
    set_charge(reaction, "p_solvent_charge", "volume", q(1.0, "mL"), E11)
    reaction["participants"] = [p for p in reaction["participants"] if p["participant_id"] not in {
        "p_solvent_feed", "p_water", "p_extract_EtOAc", "p_brine", "p_drying",
        "p_column_EtOAc", "p_column_hexane",
    }]
    reaction["participants"].append(participant(
        "p_solvent_feed", mid("PhMe"), "solvent", E11, volume=q(1.0, "mL")))
    set_charge(reaction, "p_product", "mass", q(product_mass, "mg"), E13)
    reaction["participants"].extend(workup_participants(15, E11, "1:4"))
    reaction["solvent_conditions"]["composition"] = sf("toluene", [mid("PhMe")], E11)
    reaction["solvent_conditions"]["total_volume"] = sf(
        "2.0 mL total (1.0 mL charge + 1.0 mL feed)", q(2.0, "mL"), E11)
    reaction["procedure"]["original_text"] = sf(
        "Charge sulfonamidoindole, tri(p-tolyl)phosphine, and 1.0 mL toluene; add beta-prime-acetoxyallenoate in 1.0 mL toluene over 30 min at 80 °C; hold 36 h; quench with water; extract with EtOAc; wash with brine; dry; concentrate; chromatograph (EtOAc/hexane 1:4).",
        "Staged 80 °C allenoate addition, 36 h hold, aqueous workup, and silica purification.",
        E11 + E13,
    )
    reaction["procedure"]["general_procedure_reference"] = sf(
        "vi. General Procedure for the Synthesis of Compounds 6aa-6ag, 6ba-6da, 6bf, 6ga, and 6hb",
        "experimental section vi",
        E11 + E13,
    )
    reaction["procedure"]["addition_order"] = sf(
        "sulfonamidoindole + tri(p-tolyl)phosphine + toluene; then beta-prime-acetoxyallenoate/toluene feed",
        ["sulfonamidoindole", "tri(p-tolyl)phosphine", "toluene", "beta-prime-acetoxyallenoate solution"],
        E11,
    )
    reaction["steps"] = []
    add_step(reaction, "initial_charge", E11,
             ["p_sub1", "p_cat", "p_solvent_charge"],
             "charge sulfonamidoindole, tri(p-tolyl)phosphine, and toluene")
    add_step(reaction, "allenoate_addition", E11,
             ["p_sub2", "p_solvent_feed"],
             "add beta-prime-acetoxyallenoate solution gradually",
             q(80, "°C"), q(30, "min"))
    add_step(reaction, "reaction_hold", E11,
             ["p_sub1", "p_sub2", "p_cat", "p_solvent_charge", "p_solvent_feed"],
             "hold after addition and monitor by TLC",
             q(80, "°C"), q(36, "h"))
    append_workup_steps(reaction, E11, 15, "1:4")
    replace_product_amount(reaction, mid(product), q(product_mass, "mg"), E13)
    enrich_yield(reaction, E13)
    add_note(reaction, f"Specific paragraph reports {product} as a " +
             ("white" if product in {"6da", "6hb"} else "brown") +
             f" solid, mp {mp[0]}-{mp[1]} °C.", E13)
    if printed_amount is not None:
        conflict_text = (
            f"Source conflict: Scheme 5 and general procedure vi state 0.20 mmol {sub1}, "
            f"whereas the specific {product} paragraph prints 0.04 mmol; both values are retained."
        )
        add_note(reaction, conflict_text, E13)
        reaction["review"]["unresolved"].append(conflict_text)
    reaction["review"]["unresolved"] = unique(reaction["review"]["unresolved"] + [
        "Absolute configuration remains unassigned; the explicit Z imine geometry and relative drawing are not converted into an absolute stereochemical claim.",
        "The diastereomeric-ratio measurement method is not reported in Scheme 5.",
        "Amounts of drying agent and chromatography solvents are unreported.",
    ])
    reaction["review"]["status"] = "needs_review"
    include_reaction(reaction)
    scheme5_ids.append(reaction["reaction_id"])


# Complete the specific 7ha procedure.
reaction = get_reaction("scheme6_7ha")
reaction["evidence_ids"] = unique(reaction["evidence_ids"] + E13)
reaction["location"] = sf(
    "Scheme 6 and experimental section vii, page 13",
    "Scheme 6 and experimental section vii, 7ha",
    evidence(5) + E13,
)
set_7ha_geometry(E13)
set_charge(reaction, "p_1h", "mass", q(77.2, "mg"), E13)
set_charge(reaction, "p_1h", "amount", q(0.20, "mmol"), E13)
set_charge(reaction, "p_2a", "mass", q(62.4, "mg"), E13)
set_charge(reaction, "p_2a", "amount", q(0.24, "mmol"), E13)
set_charge(reaction, "p_cat", "mass", q(10.5, "mg"), E13)
set_charge(reaction, "p_cat", "amount", q(0.04, "mmol"), E13)
set_charge(reaction, "p_cat", "mol_percent", q(20, "mol%"), E13)
set_charge(reaction, "p_solvent", "volume", q(2.0, "mL"), E13)
set_charge(reaction, "p_product", "mass", q(91.4, "mg"), E13)
reaction["participants"] = [p for p in reaction["participants"] if p["participant_id"] not in {
    "p_water", "p_extract_EtOAc", "p_brine", "p_drying", "p_column_EtOAc", "p_column_hexane",
}]
reaction["participants"].extend(workup_participants(10, E13, "2:8"))
reaction["reaction_scale"] = sf(
    "1h (0.20 mmol)",
    {"reference_molecule_id": mid("1h"), "amount": q(0.20, "mmol")},
    E13,
)
reaction["solvent_conditions"]["composition"] = sf("toluene", [mid("PhMe")], E13)
reaction["solvent_conditions"]["total_volume"] = sf("2.0 mL", q(2.0, "mL"), E13)
reaction["procedure"]["original_text"] = sf(
    "Charge 1h, delta-acetoxy allenoate 2a, and toluene; add PPh3 at room temperature (25 °C); stir 36 h; quench with water; extract with EtOAc; wash with brine; dry; concentrate; chromatograph (EtOAc/hexane 2:8).",
    "Room-temperature alpha-umpolung addition, aqueous workup, and silica purification.",
    E13,
)
reaction["procedure"]["addition_order"] = sf(
    "1h + 2a + toluene; then PPh3",
    ["1h", "2a", "toluene", "PPh3"],
    E13,
)
reaction["steps"] = []
add_step(reaction, "charge", E13, ["p_1h", "p_2a", "p_solvent"],
         "charge substrates and toluene", q(25, "°C"))
add_step(reaction, "catalyst_addition", E13, ["p_cat"],
         "add PPh3", q(25, "°C"))
add_step(reaction, "reaction_hold", E13,
         ["p_1h", "p_2a", "p_cat", "p_solvent"],
         "stir and monitor by TLC", q(25, "°C"), q(36, "h"))
append_workup_steps(reaction, E13, 10, "2:8")
replace_product_amount(reaction, mid("7ha"), q(91.4, "mg"), E13)
enrich_yield(reaction, E13)
add_note(reaction, "Specific paragraph reports 7ha as a white solid, mp 103-105 °C.", E13)
reaction["review"]["unresolved"] = unique(reaction["review"]["unresolved"] + [
    "Complete source-verified structure and atom mapping remain unresolved; explicit 2Z,4E geometry is retained textually.",
    "Amounts of drying agent and chromatography solvents are unreported.",
])
reaction["review"]["status"] = "needs_review"
include_reaction(reaction)
scheme6_id = reaction["reaction_id"]


def rebuild_scaleup(product, sub2, catalyst, additive, temp, yld, product_mass,
                    sub2_mass, catalyst_mass, ratio):
    reaction = get_reaction(f"scaleup_{product}")
    reaction["evidence_ids"] = unique(reaction["evidence_ids"] + E13)
    reaction["location"] = sf(
        f"Scheme 9 and experimental section viii, page 13, {product}",
        f"Scheme 9 and scale-up procedure, {product}",
        evidence(7) + E13,
    )
    parts = [
        participant("p_1a", mid("1a"), "reactant", E13,
                    mass=q(314.0, "mg"), amount=q(1.0, "mmol")),
        participant("p_sub2", mid(sub2), "reactant", E13,
                    mass=q(sub2_mass, "mg"), amount=q(1.20, "mmol")),
        participant("p_cat", mid(catalyst), "catalyst", E13,
                    mass=q(catalyst_mass, "mg"), amount=q(0.20, "mmol"),
                    molpct=q(20, "mol%"), ref=mid("1a")),
    ]
    if additive:
        parts.append(participant(
            "p_base", mid(additive), "additive", E13,
            mass=q(651.6, "mg"), amount=q(2.0, "mmol"),
            equiv=q(2.0, "equiv"), ref=mid("1a")))
    if product == "3aa":
        parts.append(participant("p_solvent", mid("PhMe"), "solvent", E13, volume=q(10.0, "mL")))
    else:
        parts.extend([
            participant("p_solvent_charge", mid("PhMe"), "solvent", E13, volume=q(5.0, "mL")),
            participant("p_solvent_feed", mid("PhMe"), "solvent", E13, volume=q(5.0, "mL")),
        ])
    parts.extend(workup_participants(10, E13, ratio))
    parts.append(participant(
        "p_product", mid(product), "product", E13,
        mass=q(product_mass, "mg"), product_role="desired product"))
    reaction["participants"] = parts
    reaction["reaction_scale"] = sf(
        "1a (1.0 mmol)",
        {"reference_molecule_id": mid("1a"), "amount": q(1.0, "mmol")},
        E13,
    )
    reaction["solvent_conditions"]["composition"] = sf("toluene", [mid("PhMe")], E13)
    reaction["solvent_conditions"]["total_volume"] = sf("10.0 mL", q(10.0, "mL"), E13)
    reaction["special_conditions"]["reactor"] = sf(
        "25 mL round-bottomed flask", "25 mL round-bottomed flask", E13)
    reaction["procedure"]["original_text"] = sf(
        f"One-millimole scale synthesis of {product} in a 25 mL round-bottomed flask, followed by water quench, ethyl acetate extraction, brine wash, Na2SO4 drying, concentration, and silica chromatography (EtOAc/hexane {ratio}).",
        f"One-millimole scale synthesis of {product} with staged reaction and workup.",
        E13,
    )
    reaction["steps"] = []
    if product == "3aa":
        add_step(reaction, "initial_charge", E13, ["p_1a", "p_sub2", "p_solvent"],
                 "charge substrates and toluene", "room temperature")
        add_step(reaction, "catalyst_addition", E13, ["p_cat"],
                 "add PPh3", "room temperature")
        add_step(reaction, "reaction_hold", E13,
                 ["p_1a", "p_sub2", "p_cat", "p_solvent"],
                 "stir", "room temperature", q(36, "h"))
    else:
        initial_ids = ["p_1a", "p_cat", "p_solvent_charge"]
        if additive:
            initial_ids.append("p_base")
        add_step(reaction, "initial_charge", E13, initial_ids,
                 "charge initial components and toluene")
        add_step(reaction, "allenoate_addition", E13,
                 ["p_sub2", "p_solvent_feed"],
                 "add allenoate solution gradually in an oil bath",
                 q(temp, "°C"), q(45, "min"))
        add_step(reaction, "reaction_hold", E13,
                 initial_ids + ["p_sub2", "p_solvent_feed"],
                 "hold after addition", q(temp, "°C"), q(36, "h"))
    append_workup_steps(reaction, E13, 10, ratio)
    replace_product_amount(reaction, mid(product), q(product_mass, "mg"), E13)
    enrich_yield(reaction, E13)
    reaction["review"]["unresolved"] = unique(reaction["review"]["unresolved"] + [
        "Complete source-verified structures and atom mapping remain unresolved.",
        "Amounts of drying agent and chromatography solvents are unreported.",
    ])
    reaction["review"]["status"] = "needs_review"
    include_reaction(reaction)
    return reaction["reaction_id"]


scaleup_ids = [
    rebuild_scaleup("3aa", "2a", "PPh3", None, 25, 79, 406.2, 312.0, 52.5, "1:9"),
    rebuild_scaleup("4aa", "2a", "PPh3", "Cs2CO3", 80, 72, 369.0, 312.0, 52.5, "1:9"),
    rebuild_scaleup("6aa", "5a", "pTol3P", None, 80, 75, 432.0, 386.5, 61.0, "1:4"),
]


# General records referenced by the enriched experiments are included unchanged.
supporting_generic = []
for suffix in ("generic_3_series", "generic_4_series", "generic_6_series"):
    generic = get_reaction(suffix)
    include_reaction(generic)
    supporting_generic.append(generic["reaction_id"])


bundle["source_context"].extend([
    {
        "context_id": f"{DOC}:characterization_page13",
        "kind": "characterization_context",
        "description": "Page 13 contains NMR, IR, HRMS, physical-form, and melting-point characterization for 6ca, 6da, 6bf, 6ga, 6hb, and 7ha. Peak lists, calculated masses, and found masses are not synthetic experiments.",
        "evidence_ids": E13,
    },
    {
        "context_id": f"{DOC}:associated_content_notice",
        "kind": "supporting_information_notice",
        "description": "The article states that Supporting Information is available and lists its contents. No SI file was explicitly provided or associated in this workspace, so no SI evidence or experiments are imported.",
        "evidence_ids": E13,
    },
    {
        "context_id": f"{DOC}:ccdc_and_author_information",
        "kind": "crystallographic_and_bibliographic_metadata",
        "description": "Page 14 reports CCDC 2064415-2064419, author/contact information, acknowledgments, and the beginning of the references; it contains no current-paper synthetic experiment.",
        "evidence_ids": E14,
    },
    {
        "context_id": f"{DOC}:references_pages14_16",
        "kind": "bibliography",
        "description": "Pages 14-16 are literature references and cited prior work. They are not counted as experiments performed in this paper.",
        "evidence_ids": unique(E14 + E15 + E16),
    },
])


enriched_ids = scheme5_ids + [scheme6_id] + scaleup_ids
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {
                    "region": "top-left continuation from page 12",
                    "finding": "Analytical characterization continuation for 6ca; no new synthesis record was created.",
                    "evidence_ids": E13,
                },
                {
                    "region": "left column, compound sections 6da, 6bf, 6ga, and 6hb",
                    "finding": "Four existing Scheme 5 experiments were enriched with compound-specific masses, explicit Z headings, isolated product masses/yields, physical forms, and melting points. The printed 0.04 mmol values for 1b in 6bf and 1h in 6hb conflict with the 0.20 mmol Scheme 5/general procedure and are retained as conflicts.",
                    "evidence_ids": E13,
                },
                {
                    "region": "right column, experimental section vii",
                    "finding": "The existing 7ha experiment was enriched with ordered charge, catalyst-addition, 25 °C/36 h hold, workup, purification, isolated mass/yield, and explicit 2Z,4E geometry.",
                    "evidence_ids": E13,
                },
                {
                    "region": "right column, experimental section viii",
                    "finding": "Three existing 1.0 mmol scale-up experiments for 3aa, 4aa, and 6aa were matched to their full procedures; 45 min additions and subsequent 36 h holds are separate steps for 4aa and 6aa.",
                    "evidence_ids": E13,
                },
                {
                    "region": "bottom-right Associated Content",
                    "finding": "Supporting Information availability is mentioned, but no SI file was supplied or associated; the notice is retained only as source context.",
                    "evidence_ids": E13,
                },
            ],
            "record_ids": enriched_ids,
            "unresolved": [
                "6bf and 6hb contain source amount conflicts: specific paragraphs print 0.04 mmol while Scheme 5/general procedure vi state 0.20 mmol.",
                "Complete source-verified structures/isomeric SMILES and atom mappings are unavailable for the complex products; explicit geometry is retained textually.",
                "Scheme 5 does not report the method used to determine dr.",
                "The Supporting Information document was not supplied and is not associated.",
            ],
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page: Accession Codes, Author Information, Acknowledgments, and references 1-6",
                    "finding": "Crystallographic/bibliographic metadata and cited prior work only; no current-paper synthetic experiment or reaction arrow is present.",
                    "evidence_ids": E14,
                },
            ],
            "record_ids": [],
            "unresolved": ["CCDC entries are external crystallographic records and were not imported as reaction experiments."],
        },
        {
            "page": 15,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page: references 6-14",
                    "finding": "Literature citations and cited prior work only; no experiment performed in the current paper is reported on this page.",
                    "evidence_ids": E15,
                },
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 16,
            "state": "attempted",
            "regions": [
                {
                    "region": "full page: references 14-18",
                    "finding": "Literature citations and cited prior work only; no experiment performed in the current paper is reported on this page.",
                    "evidence_ids": E16,
                },
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}


inventory = {
    "document_id": DOC,
    "batch_id": BATCH,
    "visually_inspected_pages": [13, 14, 15, 16],
    "new_reactions": [],
    "enriched_reactions": enriched_ids,
    "supporting_generic_records_in_submission": supporting_generic,
    "submitted_reaction_count": len(bundle["reactions"]),
    "registered_molecule_count_in_bundle": len(bundle["molecules"]),
    "source_context_added": [item["context_id"] for item in bundle["source_context"]],
    "chemistry_accuracy_verified": False,
    "notes": [
        "All eight page-13 experiments were matched to existing experiment identities; no duplicate experiment was created.",
        "All machine records remain needs_review.",
        "Characterization peak lists and bibliography entries were not converted into reactions.",
        "No Supporting Information document was associated.",
        "Printed 6bf/6hb amount conflicts were retained rather than corrected.",
    ],
}


for reaction in bundle["reactions"]:
    reaction["review"]["status"] = "needs_review"
    for fact in sr.facts(reaction):
        fact["evidence_ids"] = unique(fact["evidence_ids"])


(RUN / "coverage").mkdir(parents=True, exist_ok=True)
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0013_0016.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages_0013_0016_extraction.json").write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "new_reactions": 0,
    "enriched_reactions": len(enriched_ids),
    "submitted_reactions": len(bundle["reactions"]),
    "molecules": len(bundle["molecules"]),
}, ensure_ascii=False))
