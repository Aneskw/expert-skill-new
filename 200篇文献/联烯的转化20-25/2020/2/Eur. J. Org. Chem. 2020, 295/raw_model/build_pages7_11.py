#!/usr/bin/env python3
"""Build the page-7--11 v3.1 incoming bundle for doc_88f54768a8114d5d.

The script starts from the standard-init template, copies only the records that
are continued on pages 7--11, enriches them from the visible experimental text,
and adds the explicitly reported 7d -> 8a execution.  It does not elevate any
machine record beyond needs_review.
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
import standard_records as st
import reaction_tools as rt

DOC = "doc_88f54768a8114d5d"
old = json.loads((RUN / "standard_records.json").read_text())
bundle = json.loads((RUN / "raw_model/submitted_pages7_11_init.json").read_text())
old_mols = {m["molecule_id"]: m for m in old["molecules"]}
old_rxns = {r["reaction_id"]: r for r in old["reactions"]}
molmap: dict[str, dict] = {}
reactions: list[dict] = []


def ev(*pages: int) -> list[str]:
    out = []
    for page in pages:
        out.extend([f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"])
    return out


def q(number, unit, raw_text=None, **extra):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    value.update(extra)
    return value


def fact(raw, normalized, evidence, method="original_extraction", status=None,
         derivation=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, status=status,
                   derivation=derivation, alternatives=alternatives)


def direct(value, evidence):
    return fact(value, value, evidence)


def inferred(value, evidence, derivation, raw=None):
    return fact(raw, value, evidence, method="inference", derivation=derivation)


def inherited(value, evidence, derivation="Inherited only where the product paragraph explicitly cites the named procedure."):
    return fact(value, value, evidence, method="inherited_general_condition", derivation=derivation)


def missing(status="unreported"):
    return fact(None, None, (), method="not_applicable" if status == "not_applicable" else "original_extraction", status=status)


def conflict(old_fact, new_raw, new_norm, evidence, label):
    new_fact = direct(new_norm, evidence)
    new_fact["raw_value"] = new_raw
    return fact({"conflict": label}, None,
                list(dict.fromkeys(old_fact.get("evidence_ids", []) + list(evidence))),
                status="conflict", alternatives=[copy.deepcopy(old_fact), new_fact])


def take_mol(mid: str) -> dict:
    if mid not in molmap:
        molmap[mid] = copy.deepcopy(old_mols[mid])
    return molmap[mid]


def ensure_unresolved_molecule(mid: str, label: str, name: str, pages: tuple[int, ...],
                               chirality="unresolved", assignment=None,
                               geometry="not_applicable", mixture=None) -> dict:
    evidence = ev(*pages)
    if mid in old_mols:
        m = take_mol(mid)
    else:
        m = st.new_molecule(mid, DOC, evidence)
        molmap[mid] = m
    m["identity_status"] = "source_identified"
    if not m.get("occurrences"):
        m["occurrences"] = [{"document_id": DOC, "compound_label": direct(label, evidence), "evidence_ids": evidence}]
    else:
        occ = m["occurrences"][0]
        occ["compound_label"] = direct(label, evidence)
        occ["evidence_ids"] = list(dict.fromkeys(occ.get("evidence_ids", []) + evidence))
    m["names"]["original"] = direct(name, evidence)
    m["names"]["english"] = direct(name, evidence)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = direct(chirality, evidence)
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence,
        "The named products contain no assigned E/Z alkene; terminal vinyl groups do not define E/Z."
        if geometry == "not_applicable" else "Geometry status retained conservatively from the source description.")
    m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", evidence,
        "The source heading gives stereochemical descriptors, but no stereo-safe candidate SMILES was asserted from these text-only experimental pages.")
    m["stereochemistry"]["atom_centers"] = direct(
        assignment if assignment is not None else "stereogenic centers unresolved from this page", evidence)
    m["stereochemistry"]["double_bonds"] = inferred(
        [], evidence, "No E/Z descriptor is applicable to the terminal vinyl group; other unsaturation is retained only in the source name.")
    m["stereochemistry"]["axial"] = inferred(
        "not_applicable", evidence, "No stereogenic allene or biaryl axis is present in the named product.")
    m["stereochemistry"]["assignment"] = direct(assignment, evidence) if assignment else missing()
    m["stereochemistry"]["mixture"] = direct(mixture, evidence) if mixture else missing("not_applicable")
    m["stereochemistry"]["isomeric_smiles"] = missing()
    m["stereochemistry"]["identity_key"] = missing()
    m["identity_note"] = inferred(
        "Source-identified by compound heading; full connectivity and stereochemical encoding remain unverified.",
        evidence, "Identity is retained from the printed compound heading without treating name parsing as structure verification.")
    m["external_properties"] = []
    return m


def ensure_agent(mid: str, label: str, pages=(7, 8, 9, 10)) -> dict:
    if mid in old_mols:
        return take_mol(mid)
    evidence = ev(*pages)
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = direct(label, evidence)
    m["names"]["original"] = direct(label, evidence)
    m["names"]["english"] = direct(label, evidence)
    m["names"]["abbreviations"] = missing()
    m["names"]["aliases"] = missing()
    m["cas_numbers"] = missing()
    m["representation_type"] = "unresolved"
    m["smiles"] = missing()
    m["generic_template"] = None
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment",
                "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing()
    m["stereochemistry"]["chirality_status"] = inferred(
        "unresolved", evidence, "No structure-based stereochemical claim was made for this reagent/work-up material.")
    m["stereochemistry"]["geometry_status"] = inferred(
        "unresolved", evidence, "No structure-based geometry claim was made for this reagent/work-up material.")
    m["stereochemistry"]["representation_status"] = inferred(
        "unresolved", evidence, "The agent is retained by source label; no candidate graph was needed for the reaction record.")
    m["identity_note"] = direct("Source label retained without CAS or structure expansion.", evidence)
    m["external_properties"] = []
    molmap[mid] = m
    return m


AGENTS = {
    "agent_diethyl_ether": "diethyl ether", "agent_celite": "Celite",
    "agent_methanol": "methanol", "agent_hcl": "HCl",
    "agent_petroleum_ether": "petroleum ether", "agent_benzene": "benzene",
}
for mid, label in AGENTS.items():
    ensure_agent(mid, label)


def take_rxn(suffix: str, *pages: int) -> dict:
    rid = f"{DOC}:{suffix}"
    r = copy.deepcopy(old_rxns[rid])
    evidence = ev(*pages)
    r["evidence_ids"] = list(dict.fromkeys(r.get("evidence_ids", []) + evidence))
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": []}
    r["completeness"] = "partial"
    reactions.append(r)
    return r


def pget(r, pid):
    return next(p for p in r["participants"] if p["participant_id"] == pid)


def set_charge(p, key, value, evidence, *, raw=None, as_conflict=False, conflict_label="source values differ"):
    if as_conflict and p["charges"][key]["normalized_value"] is not None:
        p["charges"][key] = conflict(p["charges"][key], raw if raw is not None else value, value, evidence, conflict_label)
    else:
        p["charges"][key] = fact(raw if raw is not None else value, value, evidence)


def add_participant(r, pid, mid, role, evidence, stage="reaction", product_role=None):
    take_mol(mid)
    try:
        p = pget(r, pid)
    except StopIteration:
        p = st.new_participant(pid, mid, role, evidence, stage)
        r["participants"].append(p)
    p["molecule_id"] = mid
    p["role"] = role
    p["stage"] = stage
    p["product_role"] = direct(product_role, evidence) if product_role else missing("not_applicable")
    return p


def add_agent_participant(r, pid, mid, role, evidence, stage="reaction"):
    ensure_agent(mid, AGENTS.get(mid, old_mols.get(mid, {}).get("names", {}).get("original", {}).get("normalized_value", mid)))
    return add_participant(r, pid, mid, role, evidence, stage)


def new_step(step_id, stage, participant_ids, evidence, operation,
             temperature=None, time=None, atmosphere=None, pressure=None,
             inherit_values=False):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = list(participant_ids)
    maker = inherited if inherit_values else direct
    s["operation"] = maker(operation, evidence)
    s["temperature"] = maker(temperature, evidence) if temperature is not None else missing()
    s["time"] = maker(time, evidence) if time is not None else missing()
    s["atmosphere"] = maker(atmosphere, evidence) if atmosphere is not None else missing()
    s["pressure"] = maker(pressure, evidence) if pressure is not None else missing()
    return s


def outcome(oid, kind, targets, value, evidence, product_role=None, yield_type=None,
            basis=None, basis_inferred=False, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    o["value"] = direct(value, evidence)
    o["product_role"] = direct(product_role, evidence) if product_role else missing("not_applicable")
    o["yield_type"] = direct(yield_type, evidence) if yield_type else missing()
    if basis is None:
        o["measurement_basis"] = missing()
    elif basis_inferred:
        o["measurement_basis"] = inferred(basis, evidence, "The material was obtained after the stated work-up/purification; the source does not use a separate analytical-yield label.")
    else:
        o["measurement_basis"] = direct(basis, evidence)
    o["scope"] = direct(scope, evidence) if scope is not None else missing()
    return o


def set_procedure(r, text, evidence, gp=None, addition=None, pretreatment=None, modifications=None):
    r["procedure"]["original_text"] = direct(text, evidence)
    r["procedure"]["general_procedure_reference"] = direct(gp, evidence) if gp else missing()
    r["procedure"]["addition_order"] = direct(addition, evidence) if addition else missing()
    r["procedure"]["pretreatment"] = direct(pretreatment, evidence) if pretreatment else missing()
    r["procedure"]["modifications"] = direct(modifications, evidence) if modifications else missing()


def set_solvent(r, composition, evidence, volume=None, ratio=None, conflict_volume=False):
    r["solvent_conditions"]["composition"] = direct(composition, evidence)
    r["solvent_conditions"]["mixture_ratio"] = direct(ratio, evidence) if ratio else missing()
    if volume is None:
        r["solvent_conditions"]["total_volume"] = missing()
    elif conflict_volume and r["solvent_conditions"]["total_volume"]["normalized_value"] is not None:
        r["solvent_conditions"]["total_volume"] = conflict(
            r["solvent_conditions"]["total_volume"], volume, volume, evidence,
            "Table footnote and product-specific experimental paragraph report different solvent volumes")
    else:
        r["solvent_conditions"]["total_volume"] = direct(volume, evidence)
    r["solvent_conditions"]["concentration"] = missing()
    r["solvent_conditions"]["concentration_basis"] = missing()


def set_scale(r, mid, amount, evidence, as_conflict=False):
    normalized = {"reference_molecule_id": mid, "amount": amount}
    if as_conflict and r["reaction_scale"]["normalized_value"] is not None:
        r["reaction_scale"] = conflict(r["reaction_scale"], amount, normalized, evidence,
                                       "Table footnote and product-specific paragraph report different scales")
    else:
        r["reaction_scale"] = fact(amount, normalized, evidence)


def add_note(r, text, evidence, method="original_extraction", derivation=None):
    r["notes"].append(fact(text if method == "original_extraction" else None, text, evidence,
                           method=method, derivation=derivation))


def classify(r, changes, evidence, limitation):
    labels = {}
    facts = []
    for i, (a, b, before, after, change_type) in enumerate(changes, 1):
        change = {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
                  "before_order": before, "after_order": after, "change_type": change_type}
        facts.append(fact("net transformation in source name/drawing", change, evidence,
                          method="structure_derivation",
                          derivation="Partial abstract bond change inferred from the named/displayed net transformation; no atom mapping was asserted."))
        if change_type in ("formed", "hydrogen_gain", "order_increased"):
            label = st.bond_label(a, b, after)
            labels.setdefault(label["code"], {**label, "supporting_change_ids": []})["supporting_change_ids"].append(f"b{i}")
    r["classification"]["mapped_reaction_smiles"] = missing()
    r["classification"]["bond_changes"] = facts
    r["classification"]["labels"] = fact(
        "net transformation in source name/drawing", list(labels.values()), evidence,
        method="structure_derivation", derivation="Formation labels generated from the partial abstract bond-change list.")
    r["classification"]["coverage"] = inferred(
        "partial" if changes else "unresolved", evidence,
        "No complete atom mapping was asserted; only source-supported net bond classes were retained.")
    r["classification"]["limitations"] = [limitation]
    r["classification"]["model_votes"] = []


def finalize(r, unresolved):
    r["review"]["unresolved"] = list(unresolved)
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["completeness"] = "partial"


# Stereo-aware identities named on pages 7--10.  No text-only name is promoted
# to a full candidate graph.
stereo_names = {
    "m_7a": ("7a", "(4R,8aS)-2-Benzyl-4-vinylhexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (7,), "specified", "(4R,8aS) stated in heading"),
    "m_7b": ("7b", "(4S,8aS)-2-Benzyl-4-vinylhexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (7,), "specified", "(4S,8aS) stated in heading; major product also described as trans by NOESY on page 2"),
    "m_7d": ("7d", "(4R,8aS)-2-(4-Methoxybenzyl)-4-vinylhexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (7,), "specified", "(4R,8aS) stated in heading"),
    "m_7c": ("7c", "(4S,8aS)-2-(4-Methoxybenzyl)-4-vinylhexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (7,), "specified", "(4S,8aS) stated in heading"),
    "m_8": ("8", "(S)-2-(4-Methoxybenzyl)-4-vinyl-3,4-dihydropyrrolo[1,2-a]pyrazin-1(2H)-one", (7,), "specified", "(S) stated in heading"),
    "m_9a": ("9a", "(S)-4-(2-Hydroxyethyl)-2-(4-methoxybenzyl)-3,4-dihydropyrrolo[1,2-a]pyrazin-1(2H)-one", (7, 8), "specified", "(S) stated in heading"),
    "m_9": ("9", "(S)-2-(1-Oxo-1,2,3,4-tetrahydropyrrolo[1,2-a]pyrazin-4-yl)ethyl 2,2,2-trifluoroacetate", (7, 8), "specified", "(S) stated in heading"),
    "m_10": ("10", "(S)-2-(6,7-Dibromo-1-oxo-1,2,3,4-tetrahydropyrrolo[1,2-a]pyrazin-4-yl)ethyl 2,2,2-trifluoroacetate", (8,), "specified", "(S) stated in heading"),
    "m_11": ("11", "(S)-6,7-Dibromo-4-(2-hydroxyethyl)-3,4-dihydropyrrolo[1,2-a]pyrazin-1(2H)-one", (8,), "specified", "(S) stated in heading"),
    "m_longamide_b": ("Longamide B", "Longamide B", (8,), "specified", "Single enantiomeric product heading and optical rotation reported; descriptor retained from the route, not encoded here"),
    "m_stylisine_d": ("Stylisine D", "Stylisine D", (8,), "specified", "Single enantiomeric product heading and optical rotation reported; descriptor retained from the route, not encoded here"),
    "m_12aa": ("12aa", "(4R,8aS)-4-(1-(3-Chloro-4-methylphenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (8,), "specified", "(4R,8aS) stated in heading"),
    "m_12ba": ("12ba", "(4S,8aS)-4-(1-(3-Chloro-4-methylphenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (8, 9), "specified", "(4S,8aS) stated in heading"),
    "m_12ab": ("12ab", "(4R,8aS)-4-(1-(3,4-Dichlorophenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9,), "specified", "(4R,8aS) stated in heading"),
    "m_12bb": ("12bb", "(4S,8aS)-4-(1-(3,4-Dichlorophenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9,), "specified", "(4S,8aS) stated in heading"),
    "m_12ac": ("12ac", "(4R,8aS)-4-(1-(3,5-Dimethoxyphenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9,), "specified", "(4R,8aS) stated in heading"),
    "m_12bc": ("12bc", "(4S,8aS)-4-(1-(3,5-Dimethoxyphenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9,), "specified", "(4S,8aS) stated in heading"),
    "m_12ad": ("12ad", "(4R,8aS)-4-(1-(4-Fluorophenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9,), "specified", "(4R,8aS) stated in heading"),
    "m_12bd": ("12bd", "(4S,8aS)-4-(1-(4-Fluorophenyl)vinyl)-2-(4-methoxybenzyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9,), "specified", "(4S,8aS) stated in heading"),
    "m_12ae": ("12ae", "(4R,8aS)-2-Benzyl-4-(1-(naphthalen-2-yl)vinyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9, 10), "specified", "(4R,8aS) stated in heading"),
    "m_12be": ("12be", "(4S,8aS)-2-Benzyl-4-(1-(naphthalen-2-yl)vinyl)hexahydropyrrolo[1,2-a]pyrazin-1(2H)-one", (9, 10), "specified", "(4S,8aS) stated in heading"),
    "m_13": ("13", "(S)-4-(1-(3-Chloro-4-methylphenyl)vinyl)-2-(4-methoxybenzyl)-3,4-dihydropyrrolo[1,2-a]pyrazin-1(2H)-one", (10,), "specified", "(S) stated in heading"),
    "m_14a": ("14a", "(S)-4-((S)-1-(3-Chloro-4-methylphenyl)-2-hydroxyethyl)-2-(4-methoxybenzyl)-3,4-dihydropyrrolo[1,2-a]pyrazin-1(2H)-one", (10,), "specified", "two (S) descriptors stated in heading"),
    "m_14": ("14", "(S)-2-(3-Chloro-4-methylphenyl)-2-((S)-1-oxo-1,2,3,4-tetrahydropyrrolo[1,2-a]pyrazin-4-yl)ethyl 2,2,2-trifluoroacetate", (10,), "specified", "two (S) descriptors stated in heading"),
    "m_15": ("15", "(S)-2-(3-Chloro-4-methylphenyl)-2-((S)-6,7-dibromo-1-oxo-1,2,3,4-tetrahydropyrrolo[1,2-a]pyrazin-4-yl)ethyl 2,2,2-trifluoroacetate", (10,), "specified", "two (S) descriptors stated in heading"),
    "m_16": ("16", "(S)-6,7-Dibromo-4-((S)-1-(3-chloro-4-methylphenyl)-2-hydroxyethyl)-3,4-dihydropyrrolo[1,2-a]pyrazin-1(2H)-one", (10,), "specified", "two (S) descriptors stated in heading"),
}
for mid, (label, name, pages, chirality, assignment) in stereo_names.items():
    ensure_unresolved_molecule(mid, label, name, pages, chirality, assignment)

ensure_unresolved_molecule("m_8a", "8a", "other enantiomer of compound 8 (8a)", (7,),
                           "partially_specified", "Author calls 8a the other enantiomer of 8; no absolute descriptor is printed for 8a on page 7")
ensure_unresolved_molecule(
    "m_17", "17",
    "2-(3-Chloro-4-methylphenyl)-2-((S)-6,7-dibromo-1-oxo-1,2,3,4-tetrahydropyrrolo[1,2-a]pyrazin-4-yl)acetic acid",
    (10, 11), "mixture", "Ring (S) descriptor stated; newly sensitive alpha stereocenter is an unassigned diastereomeric mixture",
    mixture={"type": "diastereomeric mixture", "ratio": [44, 56], "order": ["unassigned diastereomer 1", "unassigned diastereomer 2"], "inseparable": True})

# Existing reactants and source-labelled agents used by the incoming reactions.
for mid in ["m_5b", "m_6b", "m_6a"]:
    take_mol(mid)


# 5b -> 6b continuation on page 7.
r = take_rxn("scheme2_5b_to_6b", 7)
e = ev(6, 7)
set_charge(pget(r, "substrate"), "mass", q(400, "mg"), ev(7))
set_charge(pget(r, "substrate"), "amount", q(1.12, "mmol"), ev(7))
set_charge(pget(r, "reagent_2"), "volume", q(0.86, "mL"), ev(7))
set_charge(pget(r, "reagent_2"), "amount", q(11.2, "mmol"), ev(7))
set_charge(pget(r, "product"), "mass", q(264, "mg"), ev(7))
set_charge(pget(r, "product"), "amount", q(1.03, "mmol"), ev(7))
set_scale(r, "m_5b", q(1.12, "mmol"), ev(7))
r["steps"] = [
    new_step("reaction", "reaction", ["substrate", "reagent_1", "reagent_2"], e,
             "prepared following the procedure for 6a", "room temperature", q(12, "h"), inherit_values=True),
    new_step("workup", "workup", ["substrate", "reagent_1", "reagent_2"], e,
             "evaporate; redissolve in ethyl acetate; wash with saturated NaHCO3 and brine; dry over Na2SO4; evaporate", inherit_values=True),
]
r["outcomes"] = [outcome("yield", "yield", ["m_6b"], q(92, "%"), ev(7), "main product", "reported product yield", None,
                                 scope={"start": "m_5b", "end": "m_6b"})]
set_solvent(r, ["CH2Cl2"], e)
set_procedure(r, "The title compound was prepared following the procedure for 6a starting from 5b (400 mg, 1.12 mmol) and TFA (0.86 mL, 11.2 mmol) in 92 % yield (264 mg, 1.03 mmol).", ev(7), gp="procedure for 6a")
classify(r, [("N", "H", 0, 1, "hydrogen_gain")], e, "Boc deprotection is represented only as net N-H gain; leaving fragments are not atom mapped.")
finalize(r, ["The product-specific CH2Cl2 volume and inherited work-up volumes are not printed.", "Yield measurement basis is not explicitly named."])


# Detailed Ag cyclisation of 6b: same 83%/17:83 execution as Table 1 entry r.
r = take_rxn("table1_r", 7)
e = ev(3, 7)
for pid, key, val in [("substrate", "mass", q(50, "mg")), ("substrate", "amount", q(0.195, "mmol")),
                      ("additive_1", "mass", q(15.6, "mg")), ("additive_1", "amount", q(0.156, "mmol")),
                      ("metal", "mass", q(26.6, "mg")), ("metal", "amount", q(0.156, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(7))
set_charge(pget(r, "metal"), "mol_percent", q(80, "mol%"), e)
set_charge(pget(r, "additive_1"), "mol_percent", q(80, "mol%"), ev(7))
pget(r, "additive_1")["charges"]["mol_percent"] = inferred(
    q(80, "mol%"), ev(7),
    "Calculated from the explicitly reported 0.156 mmol CaCO3 and 0.195 mmol 6b doses.",
    raw={"CaCO3": {"number": 0.156, "unit": "mmol"}, "6b": {"number": 0.195, "unit": "mmol"}},
)
pget(r, "additive_1")["equivalence_reference"] = direct("m_6b", ev(7))
for pid, mass, amount in [("product_7a", q(7, "mg"), q(0.027, "mmol")), ("product_7b", q(34.4, "mg"), q(0.134, "mmol"))]:
    set_charge(pget(r, pid), "mass", mass, ev(7)); set_charge(pget(r, pid), "amount", amount, ev(7))
for pid, mid, role, stage in [("ether_workup", "agent_diethyl_ether", "workup_agent", "workup"),
                              ("drying_agent", "agent_anhydrous_na2so4", "workup_agent", "workup"),
                              ("silica", "agent_silica_gel", "purification_agent", "purification")]:
    add_agent_participant(r, pid, mid, role, ev(7), stage)
set_charge(pget(r, "ether_workup"), "volume", {"number": 15, "unit": "mL", "raw_text": "3 × 5 mL"}, ev(7), raw="3 × 5 mL")
set_scale(r, "m_6b", q(0.195, "mmol"), ev(7))
set_solvent(r, ["acetone", "water"], ev(7), q(5, "mL"), "1:1")
r["steps"] = [
    new_step("charge_and_hold", "reaction", ["substrate", "metal", "additive_1", "solvent_1", "solvent_2"], ev(7),
             "dissolve 6b in acetone/water in a foil-wrapped flask; add CaCO3 and AgNO3; stir", "room temperature", "overnight"),
    new_step("filter_and_remove_acetone", "workup", ["substrate", "metal", "additive_1", "solvent_1", "solvent_2"], ev(7),
             "filter and remove acetone under reduced pressure; do not heat over 30 °C", {"number": 30, "unit": "°C", "qualifier": "lt", "raw_text": "not ... heated over 30 °C"}),
    new_step("extract", "workup", ["ether_workup", "drying_agent"], ev(7),
             "extract remaining aqueous solution with diethyl ether (3 × 5 mL), dry with Na2SO4, concentrate"),
    new_step("purify", "purification", ["silica", "product_7a", "product_7b"], ev(7), "silica-gel chromatography"),
]
r["outcomes"] = [
    outcome("yield_combined", "yield", ["m_7a", "m_7b"], q(83, "%"), e, "combined diastereomers", "combined yield", "isolated", scope={"start": "m_6b", "end": ["m_7a", "m_7b"]}),
    outcome("yield_7a", "yield", ["m_7a"], q(14, "%"), ev(7), "diastereomer 7a", "product-specific yield", "isolated after chromatography", True, {"start": "m_6b", "end": "m_7a"}),
    outcome("yield_7b", "yield", ["m_7b"], q(69, "%"), ev(7), "diastereomer 7b", "product-specific yield", "isolated after chromatography", True, {"start": "m_6b", "end": "m_7b"}),
    outcome("ratio", "selectivity", ["m_7a", "m_7b"], {"metric": "diastereomer_ratio", "order": ["7a", "7b"], "ratio": [17, 83]}, e, "diastereomer ratio", basis=None, scope={"metric": "dr", "order": ["m_7a", "m_7b"]}),
]
set_procedure(r, "6b (50 mg, 0.195 mmol) was dissolved in acetone and water (5 mL, 1:1); CaCO3 (15.6 mg, 0.156 mmol) and AgNO3 (26.6 mg, 0.156 mmol) were added and the mixture was stirred overnight at room temperature.", ev(7), pretreatment="foil-wrapped round-bottom flask")
add_note(r, "The reaction mixture should not at any point be heated over 30 °C due to formation of a purple colored side product.", ev(7))
classify(r, [("C", "N", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "C", 2, 1, "order_decreased")], e,
         "Cyclisation bond classes are partial; no complete atom mapping or proton source assignment was asserted.")
finalize(r, ["The method used to establish the 7a/7b ratio is not explicitly stated; the Table 1 NMR heading applies to yield.", "No pressure is reported."])


# 6a -> 7d + 7c, inherited from the detailed 7a/7b procedure.
r = take_rxn("scheme3_6a_to_7c", 7)
e = ev(3, 7)
pget(r, "agent_1")["role"] = "catalyst"; pget(r, "agent_2")["role"] = "additive"
for pid, key, val in [("substrate", "mass", q(1.3, "g")), ("substrate", "amount", q(4.54, "mmol")),
                      ("agent_2", "mass", q(363, "mg")), ("agent_2", "amount", q(3.63, "mmol")),
                      ("agent_1", "mass", q(618, "mg")), ("agent_1", "amount", q(3.63, "mmol")),
                      ("product", "mass", q(896, "mg")), ("product", "amount", q(3.13, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(7))
minor = add_participant(r, "product_7d", "m_7d", "product", ev(7), product_role="minor diastereomer 7d")
set_charge(minor, "mass", q(184, "mg"), ev(7)); set_charge(minor, "amount", q(0.64, "mmol"), ev(7))
set_scale(r, "m_6a", q(4.54, "mmol"), ev(7))
set_solvent(r, ["acetone", "water"], ev(3, 7), ratio="1:1")
r["steps"] = [
    new_step("reaction", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4"], ev(7),
             "prepared following the procedure for 7a/7b", "room temperature", "overnight", inherit_values=True),
    new_step("workup_and_purification", "purification", ["product_7d", "product"], ev(7),
             "work-up and silica-gel chromatography following the procedure for 7a/7b", inherit_values=True),
]
r["outcomes"] = [
    outcome("yield_7d", "yield", ["m_7d"], q(14, "%"), ev(7), "minor diastereomer 7d", "product-specific yield", "isolated after chromatography", True, {"start": "m_6a", "end": "m_7d"}),
    outcome("yield_7c", "yield", ["m_7c"], q(69, "%"), ev(7), "major diastereomer 7c", "product-specific yield", "isolated after chromatography", True, {"start": "m_6a", "end": "m_7c"}),
]
set_procedure(r, "The title compounds were prepared following the procedure for 7a/7b starting from 6a (1.3 g, 4.54 mmol), CaCO3 (363 mg, 3.63 mmol) and AgNO3 (618 mg, 3.63 mmol).", ev(7), gp="procedure for 7a/7b")
add_note(r, "7c quickly turns purple upon heat and air exposure.", ev(7))
classify(r, [("C", "N", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("C", "C", 2, 1, "order_decreased")], e,
         "Cyclisation bond classes are partial; exact product atom mapping was not asserted.")
finalize(r, ["The product paragraph does not repeat the acetone/water volume; only the 1:1 procedure relationship is inherited.", "No dr measurement method is reported for this scale."])


# 7c -> 8 detailed oxidation.
r = take_rxn("scheme3_7c_to_8", 7)
e = ev(3, 7)
for pid, key, val in [("substrate", "mass", q(800, "mg")), ("substrate", "amount", q(2.80, "mmol")),
                      ("agent_1", "mass", q(6.16, "g")), ("agent_1", "amount", q(70.80, "mmol")),
                      ("agent_2", "volume", q(20, "mL")), ("product", "mass", q(450, "mg")), ("product", "amount", q(1.60, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(7))
add_agent_participant(r, "celite", "agent_celite", "workup_agent", ev(7), "workup")
add_agent_participant(r, "silica", "agent_silica_gel", "purification_agent", ev(7), "purification")
set_scale(r, "m_7c", q(2.80, "mmol"), ev(7)); set_solvent(r, ["THF"], ev(7), q(20, "mL"))
r["steps"] = [
    new_step("charge", "reaction", ["substrate", "agent_1", "agent_2"], ev(7), "charge 7c, freshly distilled THF and MnO2 to a flame-dried flask", atmosphere=None),
    new_step("heat", "reaction", ["substrate", "agent_1", "agent_2"], ev(7), "heat and stir", q(85, "°C"), "overnight"),
    new_step("first_filtration", "workup", ["celite"], ev(7), "filter through Celite and wash the Celite with hot THF"),
    new_step("celite_reflux", "workup", ["celite", "agent_2"], ev(7), "place Celite in a round-bottom flask and reflux in THF", "reflux", q(30, "min")),
    new_step("second_filtration", "workup", ["celite"], ev(7), "filter through a new pad of Celite and concentrate combined filtrate under reduced pressure"),
    new_step("purify", "purification", ["silica", "product"], ev(7), "silica-gel chromatography"),
]
r["outcomes"] = [outcome("yield", "yield", ["m_8"], q(57, "%"), ev(7), "main product", "product yield", "isolated after chromatography", True, {"start": "m_7c", "end": "m_8"})]
set_procedure(r, "7c (800 mg, 2.80 mmol), freshly distilled THF (20 mL) and MnO2 (6.16 g, 70.80 mmol) were heated to 85 °C and stirred overnight.", ev(7), pretreatment="flame-dried round-bottom flask; freshly distilled THF")
classify(r, [("C", "C", 1, 2, "order_increased")], e, "Aromatizing oxidation is represented by a partial C=C bond-order increase; the full dehydrogenation mapping is unresolved.")
finalize(r, ["The number and exact positions of all bond-order changes and hydrogen losses require structure-level expert mapping."])


# Explicit low-stability 7d -> 8a execution, reported without scale or yield.
r = copy.deepcopy(r)
r["reaction_id"] = f"{DOC}:page7_7d_to_8a"
r["record_type"] = "specific_experiment"
r["experiment_categories"] = direct(["derivatization"], ev(7))
r["location"] = direct({"page": 7, "location": "footnote below compounds 7d/7c and characterization of 8a"}, ev(7)); r["location"]["raw_value"] = "Page 7 footnote: Due to its low stability ... 8a"
r["participants"] = []
add_participant(r, "substrate", "m_7d", "reactant", ev(7))
add_participant(r, "oxidant", "agent_mno2", "reagent", ev(7))
add_participant(r, "solvent", "agent_thf", "solvent", ev(7))
add_participant(r, "product", "m_8a", "product", ev(7), product_role="other enantiomer of 8")
r["reaction_scale"] = missing(); set_solvent(r, ["THF"], ev(7))
r["steps"] = [new_step("reaction", "reaction", ["substrate", "oxidant", "solvent"], ev(7), "oxidized as described below for compound 8", q(85, "°C"), "overnight", inherit_values=True),
              new_step("workup_and_purification", "purification", ["product"], ev(7), "work-up and purification as described for compound 8", inherit_values=True)]
r["outcomes"] = []
set_procedure(r, "Due to its low stability the compound was oxidized as described below (compound 8) and stored and further characterized as other enantiomer of 8 named 8a.", ev(7), gp="procedure for compound 8")
r["notes"] = [direct("Due to its low stability the compound was oxidized as described below (compound 8) and stored and further characterized as other enantiomer of 8 named 8a.", ev(7))]
classify(r, [("C", "C", 1, 2, "order_increased")], ev(7), "Only the inherited aromatizing oxidation class is retained; scale, yield and complete mapping are not reported.")
finalize(r, ["The amounts of 7d, MnO2 and THF and the yield of 8a are not reported.", "8a is called the other enantiomer of 8, but no explicit absolute descriptor is printed on this page."])
reactions.append(r)


# 8 -> 9a -> 9 two-stage overall experiment across pages 7--8.
r = take_rxn("scheme3_8_to_9", 7, 8)
e = ev(3, 7, 8)
for pid, key, val in [("substrate", "mass", q(400, "mg")), ("substrate", "amount", q(1.42, "mmol")),
                      ("agent_2", "volume", q(15, "mL")), ("agent_1", "amount", q(4.26, "mmol")), ("agent_1", "volume", q(8.51, "mL")),
                      ("agent_3", "volume", q(12.76, "mL")), ("agent_4", "volume", q(1.54, "mL")),
                      ("agent_6", "volume", q(5, "mL")), ("agent_5", "volume", q(0.26, "mL")), ("agent_5", "amount", q(2.37, "mmol")),
                      ("product", "mass", q(216, "mg")), ("product", "amount", q(0.78, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(7, 8))
set_charge(pget(r, "agent_1"), "amount", q(4.26, "mmol"), ev(7), raw="9-BBN in THF (0.5 M, 4.26 mmol, 8.51 mL)")
pget(r, "agent_3")["charges"]["amount"] = inferred(
    q(12.76, "mmol"), ev(7),
    "Calculated from the explicitly reported 1.0 M concentration and 12.76 mL volume.",
    raw="NaOH (1.0 M, 12.76 mL)",
)
inter = add_participant(r, "intermediate_9a", "m_9a", "product", ev(7, 8), product_role="crude intermediate consumed in next stage")
set_charge(inter, "mass", q(352, "mg"), ev(7))
for pid, mid, role, stage in [("sat_bicarbonate", "agent_saturated_nahco3", "workup_agent", "workup"),
                              ("ethyl_acetate", "agent_ethyl_acetate", "workup_agent", "workup"),
                              ("drying_agent", "agent_anhydrous_na2so4", "workup_agent", "workup"),
                              ("silica", "agent_silica_gel", "purification_agent", "purification")]:
    add_agent_participant(r, pid, mid, role, ev(7, 8), stage)
set_scale(r, "m_8", q(1.42, "mmol"), ev(7)); set_solvent(r, ["THF", "TFA"], ev(7, 8))
r["steps"] = [
    new_step("hydroboration", "reaction", ["substrate", "agent_1", "agent_2"], ev(7), "add 9-BBN solution and monitor by TLC until alkene consumption", None, q(3, "h")),
    new_step("oxidant_addition", "reaction", ["agent_3", "agent_4"], ev(7), "cool, then add NaOH and H2O2", q(0, "°C")),
    new_step("oxidation_hold", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4"], ev(7), "stir", "room temperature", "overnight"),
    new_step("stage1_workup", "workup", ["sat_bicarbonate", "ethyl_acetate", "drying_agent"], ev(7), "transfer with saturated NaHCO3 (10 mL), extract with ethyl acetate (3 × 10 mL), dry, filter, concentrate"),
    new_step("stage1_purification", "purification", ["silica", "intermediate_9a"], ev(7), "silica-gel chromatography; obtain crude 9a slightly contaminated with 9-BBN residues"),
    new_step("deprotection_esterification", "reaction", ["intermediate_9a", "agent_6", "agent_5"], ev(7, 8), "under nitrogen dissolve crude oil in TFA, add anisole and stir", q(60, "°C"), "overnight", "nitrogen"),
    new_step("stage2_workup", "workup", ["ethyl_acetate", "sat_bicarbonate", "drying_agent"], ev(8), "evaporate; dilute with ethyl acetate (15 mL); wash with NaHCO3 (3 × 15 mL); dry, filter and concentrate"),
    new_step("stage2_purification", "purification", ["silica", "product"], ev(8), "silica-gel chromatography"),
]
r["outcomes"] = [
    outcome("crude_9a_mass", "product_amount", ["m_9a"], q(352, "mg"), ev(7), "crude intermediate consumed in next stage", basis="crude mass; slightly contaminated with 9-BBN residues", scope={"stage": "hydroboration/oxidation", "start": "m_8", "end": "m_9a"}),
    outcome("overall_yield", "yield", ["m_9"], q(55, "%"), ev(8), "main product", "overall two-step yield", "isolated after chromatography", True, {"start": "m_8", "end": "m_9", "stages": 2}),
]
set_procedure(r, "To a solution of 8 (400 mg, 1.42 mmol) in THF (15 mL), 9-BBN in THF (0.5 M, 4.26 mmol, 8.51 mL) was added; after 3 h the mixture was cooled to 0 °C and NaOH/H2O2 were added, followed by overnight stirring. Crude 9a was then treated with anisole in TFA at 60 °C overnight to afford 9 in 55 % yield over 2 steps.", ev(7, 8), addition="9-BBN; then at 0 °C NaOH and H2O2; after work-up TFA and anisole")
add_note(r, "The oil was used without further purification in the following step.", ev(7))
classify(r, [("C", "O", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("N", "H", 0, 1, "hydrogen_gain"), ("C", "O", 0, 1, "formed")], e,
         "Hydroboration/oxidation, PMB cleavage and trifluoroacetate formation are represented only by partial net bond classes.")
finalize(r, ["No isolated yield is reported for 9a; 352 mg is a contaminated crude mass.", "The 55% value is explicitly an overall two-step yield and is not assigned to the second step alone."])


# 9 -> 10.
r = take_rxn("scheme3_9_to_10", 8); e = ev(3, 8)
for pid, key, val in [("substrate", "mass", q(120, "mg")), ("substrate", "amount", q(0.44, "mmol")),
                      ("agent_2", "volume", q(5, "mL")), ("agent_1", "mass", q(162, "mg")), ("agent_1", "amount", q(0.913, "mmol")),
                      ("product", "mass", q(160, "mg")), ("product", "amount", q(0.37, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(8))
add_agent_participant(r, "thf_carrier", "agent_thf", "solvent", ev(8)); set_charge(pget(r, "thf_carrier"), "volume", q(1, "mL"), ev(8))
add_agent_participant(r, "silica", "agent_silica_gel", "purification_agent", ev(8), "purification")
set_scale(r, "m_9", q(0.44, "mmol"), ev(8)); set_solvent(r, ["THF"], ev(8), {"number": 6, "unit": "mL", "raw_text": "5 mL reaction THF + 1 mL NBS-solution THF"})
r["solvent_conditions"]["total_volume"] = inferred(
    q(6, "mL", "5 mL reaction THF + 1 mL NBS-solution THF"), ev(8),
    "Summed from the separately reported 5 mL reaction portion and 1 mL NBS-solution portion.",
    raw="5 mL + 1 mL THF",
)
r["steps"] = [new_step("addition_and_hold", "reaction", ["substrate", "agent_1", "agent_2", "thf_carrier"], ev(8), "add a solution of NBS in THF and stir", "room temperature", "overnight"),
              new_step("concentrate", "workup", ["substrate", "agent_1", "agent_2", "thf_carrier"], ev(8), "evaporate under reduced pressure"),
              new_step("purify", "purification", ["silica", "product"], ev(8), "silica-gel chromatography")]
r["outcomes"] = [outcome("yield", "yield", ["m_10"], q(85, "%"), ev(8), "main product", "product yield", "isolated after chromatography", True, {"start": "m_9", "end": "m_10"})]
set_procedure(r, "To 9 (120 mg, 0.44 mmol) in THF (5 mL) was added NBS (162 mg, 0.913 mmol) in THF (1 mL) at room temperature; the reaction was stirred overnight.", ev(8), addition="NBS solution in THF added to solution of 9")
classify(r, [("C", "Br", 0, 1, "formed"), ("C", "Br", 0, 1, "formed")], e, "Two brominations are inferred from 9 -> dibrominated product 10; atom positions are not mapped.")
finalize(r, ["No atmosphere or pressure is reported."])


# 10 -> 11.
r = take_rxn("scheme3_10_to_11", 8); e = ev(3, 8)
for pid, key, val in [("substrate", "mass", q(50, "mg")), ("substrate", "amount", q(0.115, "mmol")),
                      ("agent_1", "volume", q(0.46, "mL")), ("product", "mass", q(38.5, "mg")), ("product", "amount", q(0.114, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(8))
add_agent_participant(r, "methanol", "agent_methanol", "solvent", ev(8));
add_agent_participant(r, "ethyl_acetate", "agent_ethyl_acetate", "workup_agent", ev(8), "workup"); set_charge(pget(r, "ethyl_acetate"), "volume", {"number": 15, "unit": "mL", "raw_text": "3 × 5 mL"}, ev(8), raw="3 × 5 mL")
add_agent_participant(r, "drying_agent", "agent_anhydrous_na2so4", "workup_agent", ev(8), "workup")
set_scale(r, "m_10", q(0.115, "mmol"), ev(8)); set_solvent(r, ["methanol", "THF", "water"], ev(8), q(5, "mL"), "1:2:1")
r["steps"] = [new_step("add_base", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "methanol"], ev(8), "add NaOH (1 M) dropwise", "room temperature", q(3, "h")),
              new_step("remove_organic_solvents", "workup", ["agent_2", "methanol"], ev(8), "evaporate methanol and THF under reduced pressure"),
              new_step("extract", "workup", ["ethyl_acetate", "drying_agent"], ev(8), "extract remaining aqueous solution with ethyl acetate (3 × 5 mL), dry, filter and concentrate")]
r["outcomes"] = [outcome("yield", "yield", ["m_11"], q(99, "%"), ev(8), "main product", "product yield", "isolated after work-up", True, {"start": "m_10", "end": "m_11"})]
set_procedure(r, "To 10 (50 mg, 0.115 mmol) in methanol/THF/water (5 mL, 1:2:1), NaOH (0.46 mL of 1 M solution) was added dropwise. After 3 h methanol and THF were evaporated.", ev(8), addition="1 M NaOH added dropwise")
classify(r, [("O", "H", 0, 1, "hydrogen_gain")], e, "Trifluoroacetate hydrolysis is represented as net O-H gain; leaving-group bonds are not mapped.")
finalize(r, ["The paragraph does not independently state the reaction temperature; room temperature is supplied by Scheme 3 on page 3."])


# 11 -> Longamide B.
r = take_rxn("scheme3_11_to_longamide_b", 8); e = ev(3, 8)
for pid, key, val in [("substrate", "mass", q(15, "mg")), ("substrate", "amount", q(0.045, "mmol")),
                      ("agent_2", "mass", q(36, "mg")), ("agent_2", "amount", q(0.112, "mmol")),
                      ("agent_3", "mass", q(15, "mg")), ("agent_3", "amount", q(0.176, "mmol")),
                      ("agent_1", "mass", q(1.6, "mg")), ("agent_1", "amount", q(0.011, "mmol")),
                      ("product", "mass", q(12.5, "mg")), ("product", "amount", q(0.037, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(8))
for pid, mid, role, stage in [("sat_bicarbonate", "agent_saturated_nahco3", "workup_agent", "workup"),
                              ("ethyl_acetate", "agent_ethyl_acetate", "workup_agent", "workup"),
                              ("hcl", "agent_hcl", "workup_agent", "workup"),
                              ("drying_agent", "agent_anhydrous_na2so4", "workup_agent", "workup"),
                              ("dcm", "agent_ch2cl2", "purification_agent", "purification"),
                              ("petroleum_ether", "agent_petroleum_ether", "purification_agent", "purification"),
                              ("methanol", "agent_methanol", "purification_agent", "purification")]:
    add_agent_participant(r, pid, mid, role, ev(8), stage)
set_charge(pget(r, "sat_bicarbonate"), "volume", q(5, "mL"), ev(8))
set_scale(r, "m_11", q(0.045, "mmol"), ev(8)); set_solvent(r, ["MeCN", "water"], ev(8), q(3, "mL"), "1:1")
r["steps"] = [new_step("oxidation", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4", "agent_5"], ev(8), "stir", "room temperature", q(18, "h")),
              new_step("basic_workup", "workup", ["sat_bicarbonate", "ethyl_acetate"], ev(8), "add saturated NaHCO3 (5 mL) and wash aqueous phase with ethyl acetate (3 × 5 mL)"),
              new_step("acidify_extract", "workup", ["hcl", "ethyl_acetate", "drying_agent"], ev(8), "acidify remaining water solution with 2 M HCl to pH 2; extract with ethyl acetate (3 × 5 mL), dry, filter and concentrate"),
              new_step("wash_crude", "purification", ["dcm", "petroleum_ether"], ev(8), "wash crude solid with dichloromethane/petroleum ether (2:8)"),
              new_step("recrystallize", "purification", ["dcm", "methanol", "product"], ev(8), "recrystallize from dichloromethane/methanol (98:2)")]
r["outcomes"] = [outcome("yield", "yield", ["m_longamide_b"], q(83, "%"), ev(8), "main product", "product yield", "isolated after recrystallization", True, {"start": "m_11", "end": "m_longamide_b"})]
set_procedure(r, "11 (15 mg, 0.045 mmol), BAIB (36 mg, 0.112 mmol), NaHCO3 (15 mg, 0.176 mmol) and TEMPO (1.6 mg, 0.011 mmol) in acetonitrile/water (3 mL, 1:1) were stirred for 18 h at room temperature.", ev(8))
classify(r, [("C", "O", 1, 2, "order_increased"), ("C", "O", 0, 1, "formed"), ("O", "H", 0, 1, "hydrogen_gain")], e,
         "Primary-alcohol oxidation to carboxylic acid is represented by partial C=O, C-O and O-H changes without full atom mapping.")
finalize(r, ["No pressure or atmosphere is reported."])


# Longamide B -> Stylisine D.
r = take_rxn("scheme3_longamide_b_to_stylisine_d", 8); e = ev(3, 8)
for pid, key, val in [("substrate", "mass", q(10, "mg")), ("substrate", "amount", q(0.028, "mmol")),
                      ("agent_2", "mass", q(8.2, "mg")), ("agent_2", "amount", q(0.043, "mmol")),
                      ("agent_3", "mass", q(5.8, "mg")), ("agent_3", "amount", q(0.043, "mmol")),
                      ("agent_4", "volume", q(8.4, "uL")), ("agent_4", "amount", q(0.048, "mmol")),
                      ("agent_1", "mass", q(6.8, "mg")), ("agent_1", "amount", q(0.048, "mmol")),
                      ("agent_5", "volume", q(5, "mL")), ("product", "mass", q(8.1, "mg")), ("product", "amount", q(0.023, "mmol"))]:
    set_charge(pget(r, pid), key, val, ev(8))
set_charge(pget(r, "agent_1"), "mass", q(6.8, "mg"), ev(8), raw="25 % ammonium hydroxide (6.8 mg, 0.048 mmol)")
for pid, mid, role, stage in [("ethyl_acetate", "agent_ethyl_acetate", "workup_agent", "workup"),
                              ("hcl", "agent_hcl", "workup_agent", "workup"),
                              ("bicarbonate", "agent_nahco3", "workup_agent", "workup"),
                              ("brine", "agent_brine", "workup_agent", "workup"),
                              ("drying_agent", "agent_anhydrous_na2so4", "workup_agent", "workup"),
                              ("silica", "agent_silica_gel", "purification_agent", "purification")]:
    add_agent_participant(r, pid, mid, role, ev(8), stage)
set_scale(r, "m_longamide_b", q(0.028, "mmol"), ev(8)); set_solvent(r, ["THF"], ev(8), q(5, "mL"))
r["steps"] = [new_step("activation", "reaction", ["substrate", "agent_2", "agent_3", "agent_4", "agent_5"], ev(8), "add EDCI, HOBt and DIPEA at 0 °C and stir", q(0, "°C"), q(5, "min")),
              new_step("amine_addition_hold", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4", "agent_5"], ev(8), "add 25% ammonium hydroxide solution and stir", "room temperature", q(18, "h")),
              new_step("workup", "workup", ["ethyl_acetate", "hcl", "bicarbonate", "brine", "drying_agent"], ev(8), "evaporate; redissolve in ethyl acetate (5 mL); wash with 0.5 M HCl (3 × 5 mL), NaHCO3 (3 × 5 mL), then brine (2 × 5 mL); dry and concentrate"),
              new_step("purify", "purification", ["silica", "product"], ev(8), "silica-gel chromatography")]
r["outcomes"] = [outcome("yield", "yield", ["m_stylisine_d"], q(81, "%"), ev(8), "main product", "product yield", "isolated after chromatography", True, {"start": "m_longamide_b", "end": "m_stylisine_d"})]
set_procedure(r, "Longamide B (10 mg, 0.028 mmol), EDCI (8.2 mg, 0.043 mmol), HOBt (5.8 mg, 0.043 mmol) and DIPEA (8.4 uL, 0.048 mmol) were combined in THF at 0 °C; after 5 min 25% ammonium hydroxide was added and the solution was stirred for 18 h at room temperature.", ev(8), addition="EDCI/HOBt/DIPEA at 0 °C; after 5 min add 25% NH4OH")
classify(r, [("C", "N", 0, 1, "formed"), ("C", "O", 1, 0, "broken")], e, "Amide formation is represented by net C-N formation and loss of the acid C-O connection; activated intermediates are not mapped.")
finalize(r, ["The source does not separately name a yield measurement method beyond isolation after chromatography."])


# Pd scope records a-e: product-specific procedures enrich the Table 3 rows.
scope = {
    "a": {"pages": (8, 9), "allene": "m_6a", "allene_mass": q(1114, "mg"), "allene_amount": q(3.90, "mmol"),
          "aryl_mass": q(982.8, "mg"), "aryl_amount": q(3.90, "mmol"), "base_mass": q(1075, "mg"), "base_amount": q(7.79, "mmol"),
          "pd_mass": q(87.3, "mg"), "pd_amount": q(0.39, "mmol"), "lig_mass": q(225.5, "mg"), "lig_amount": q(0.39, "mmol"),
          "solvent_volume": q(10, "mL"), "pa": "m_12aa", "pb": "m_12ba", "ya": 32, "yb": 26,
          "ma": q(512.9, "mg"), "aa": q(1.25, "mmol"), "mb": q(416.8, "mg"), "ab": q(1.01, "mmol")},
    "b": {"pages": (9,), "allene": "m_6a", "allene_mass": q(30.0, "mg"), "allene_amount": q(0.105, "mmol"), "aryl_mass": q(28.7, "mg"), "aryl_amount": q(0.105, "mmol"), "base_amount": q(0.210, "mmol"), "pd_mass": q(2.4, "mg"), "pd_amount": q(0.0105, "mmol"), "lig_mass": q(6.1, "mg"), "lig_amount": q(0.0105, "mmol"), "pa": "m_12ab", "pb": "m_12bb", "ya": 33, "yb": 26, "ma": q(15.1, "mg"), "aa": q(0.035, "mmol"), "mb": q(11.8, "mg"), "ab": q(0.028, "mmol")},
    "c": {"pages": (9,), "allene": "m_6a", "allene_mass": q(30.0, "mg"), "allene_amount": q(0.105, "mmol"), "aryl_mass": q(27.7, "mg"), "aryl_amount": q(0.105, "mmol"), "base_amount": q(0.210, "mmol"), "pd_mass": q(2.4, "mg"), "pd_amount": q(0.0105, "mmol"), "lig_mass": q(6.1, "mg"), "lig_amount": q(0.0105, "mmol"), "pa": "m_12ac", "pb": "m_12bc", "ya": 29, "yb": 24, "ma": q(12.8, "mg"), "aa": q(0.030, "mmol"), "mb": q(10.6, "mg"), "ab": q(0.025, "mmol")},
    "d": {"pages": (9,), "allene": "m_6a", "allene_mass": q(30.0, "mg"), "allene_amount": q(0.105, "mmol"), "aryl_mass": q(23.3, "mg"), "aryl_amount": q(0.105, "mmol"), "base_amount": q(0.210, "mmol"), "pd_mass": q(2.4, "mg"), "pd_amount": q(0.0105, "mmol"), "lig_mass": q(6.1, "mg"), "lig_amount": q(0.0105, "mmol"), "pa": "m_12ad", "pb": "m_12bd", "ya": 39, "yb": 29, "ma": q(15.4, "mg"), "aa": q(0.041, "mmol"), "mb": q(11.7, "mg"), "ab": q(0.031, "mmol")},
    "e": {"pages": (9, 10), "allene": "m_6b", "allene_mass": q(26.9, "mg"), "allene_amount": q(0.105, "mmol"), "aryl_mass": q(26.7, "mg"), "aryl_amount": q(0.105, "mmol"), "base_amount": q(0.210, "mmol"), "pd_mass": q(2.4, "mg"), "pd_amount": q(0.0105, "mmol"), "lig_mass": q(6.1, "mg"), "lig_amount": q(0.0105, "mmol"), "pa": "m_12ae", "pb": "m_12be", "ya": 43, "yb": 32, "ma": q(17.3, "mg"), "aa": q(0.045, "mmol"), "mb": q(12.7, "mg"), "ab": q(0.033, "mmol")},
}
for row, d in scope.items():
    r = take_rxn(f"table3_{row}", *d["pages"]); evidence = ev(4, *d["pages"])
    product_evidence = ev(*d["pages"])
    # Paragraph values are preserved against the shared 0.1 mmol table footnote.
    set_charge(pget(r, "allene"), "mass", d["allene_mass"], product_evidence)
    set_charge(pget(r, "allene"), "amount", d["allene_amount"], product_evidence, as_conflict=True,
               conflict_label="shared Table 3 footnote says 0.1 mmol; product paragraph gives a product-specific amount")
    set_charge(pget(r, "aryl_iodide"), "mass", d["aryl_mass"], product_evidence)
    set_charge(pget(r, "aryl_iodide"), "amount", d["aryl_amount"], product_evidence, as_conflict=True,
               conflict_label="shared Table 3 footnote says 0.1 mmol; product paragraph gives a product-specific amount")
    if d.get("base_mass"):
        set_charge(pget(r, "base"), "mass", d["base_mass"], product_evidence)
    set_charge(pget(r, "base"), "amount", d["base_amount"], product_evidence, as_conflict=True,
               conflict_label="shared Table 3 footnote says 0.23 mmol; product paragraph gives a product-specific amount")
    for pid, mass_key, amt_key in [("pd", "pd_mass", "pd_amount"), ("ligand", "lig_mass", "lig_amount")]:
        set_charge(pget(r, pid), "mass", d[mass_key], product_evidence)
        set_charge(pget(r, pid), "amount", d[amt_key], product_evidence, as_conflict=True,
                   conflict_label="shared Table 3 footnote and product paragraph report different scales")
    set_charge(pget(r, "product_a"), "mass", d["ma"], product_evidence); set_charge(pget(r, "product_a"), "amount", d["aa"], product_evidence)
    set_charge(pget(r, "product_b"), "mass", d["mb"], product_evidence); set_charge(pget(r, "product_b"), "amount", d["ab"], product_evidence)
    set_scale(r, d["allene"], d["allene_amount"], product_evidence, as_conflict=True)
    if row == "a":
        set_charge(pget(r, "solvent"), "volume", d["solvent_volume"], product_evidence, as_conflict=True,
                   conflict_label="Table 3 footnote gives MeCN 3 mL; product-specific paragraph gives 10 mL")
        set_solvent(r, ["MeCN"], product_evidence, d["solvent_volume"], conflict_volume=True)
        steps_evidence = product_evidence
        r["steps"] = [new_step("charge", "reaction", ["allene", "aryl_iodide", "pd", "ligand", "base", "solvent"], steps_evidence,
                               "under nitrogen add K2CO3, Pd(OAc)2, XantPhos and aryl iodide to the allene solution", atmosphere="nitrogen"),
                      new_step("reflux", "reaction", ["allene", "aryl_iodide", "pd", "ligand", "base", "solvent"], steps_evidence, "reflux", "reflux", q(18, "h"), "nitrogen"),
                      new_step("concentrate", "workup", ["allene", "aryl_iodide", "pd", "ligand", "base", "solvent"], steps_evidence, "concentrate under reduced pressure"),
                      new_step("purify", "purification", ["product_a", "product_b"], steps_evidence, "silica-gel chromatography")]
        set_procedure(r, "6a (1114 mg, 3.90 mmol) in MeCN (10 mL), K2CO3 (1075 mg, 7.79 mmol), Pd(OAc)2 (87.3 mg, 0.39 mmol), XantPhos (225.5 mg, 0.39 mmol) and aryl iodide (982.8 mg, 3.90 mmol) were refluxed for 18 h under nitrogen.", product_evidence, addition="K2CO3, Pd(OAc)2, XantPhos and aryl iodide added to 6a solution", pretreatment="nitrogen atmosphere")
    else:
        # Exact solvent volume and operations come from the shared Table 3 footnote / cited procedure, not from the product paragraph.
        set_solvent(r, ["MeCN"], ev(4), q(3, "mL"))
        r["steps"] = [new_step("reaction", "reaction", ["allene", "aryl_iodide", "pd", "ligand", "base", "solvent"], evidence,
                               "prepared following the procedure for 12aa/12ba; reflux", "reflux", q(18, "h"), "nitrogen", inherit_values=True),
                      new_step("workup_and_purification", "purification", ["product_a", "product_b"], evidence,
                               "concentrate and purify by silica-gel chromatography following the procedure for 12aa/12ba", inherit_values=True)]
        set_procedure(r, f"The title compounds were prepared following the procedure for 12aa/12ba starting from the printed product-specific charges for Table 3 entry {row}.", product_evidence, gp="procedure for 12aa/12ba")
    # Keep the explicitly isolated combined table yield and add the two isolated product yields.
    existing_combined = next(o for o in r["outcomes"] if o["outcome_id"] == "yield")
    existing_ratio = next(o for o in r["outcomes"] if o["outcome_id"] == "ratio")
    existing_combined["measurement_basis"] = direct("isolated combined yield", ev(4))
    existing_ratio["measurement_basis"] = direct("1H NMR", ev(4))
    r["outcomes"] = [existing_combined, existing_ratio,
                     outcome("yield_product_a", "yield", [d["pa"]], q(d["ya"], "%"), product_evidence, f"diastereomer {d['pa'][2:]}", "product-specific yield", "isolated after chromatography", True, {"start": d["allene"], "end": d["pa"]}),
                     outcome("yield_product_b", "yield", [d["pb"]], q(d["yb"], "%"), product_evidence, f"diastereomer {d['pb'][2:]}", "product-specific yield", "isolated after chromatography", True, {"start": d["allene"], "end": d["pb"]})]
    if row != "a":
        add_note(r, "K2CO3 is printed as '(29.0, 0.210 mmol)' without a mass unit; the 0.210 mmol amount is retained and no mass unit is invented.", product_evidence)
    add_note(r, "Shared Table 3 footnote and product-specific experimental paragraph differ in scale; both evidence-linked values are retained as conflicts.", evidence, method="inference", derivation="Comparison of page-4 footnote with the product paragraph.")
    classify(r, [("C", "C", 0, 1, "formed"), ("C", "N", 0, 1, "formed")], evidence,
             "Pd-catalysed arylative cyclisation has partial C-C and C-N formation classes; products have no asserted atom mapping.")
    unresolved = ["Table 3 shared-footnote and product-paragraph scales differ; neither value was silently preferred.",
                  "Product structures retain heading stereodescriptors but no stereo-safe SMILES was asserted."]
    if row == "a": unresolved.append("MeCN volume is 3 mL in the Table 3 footnote and 10 mL in the product-specific paragraph.")
    else: unresolved.append("The K2CO3 mass is printed without a unit; only 0.210 mmol is structured.")
    finalize(r, unresolved)


# Downstream Scheme 5 records, enriched from page 10.
r = take_rxn("scheme5_12ba_to_13", 10); e = ev(5, 10)
for pid, key, val in [("substrate", "mass", q(320, "mg")), ("substrate", "amount", q(0.78, "mmol")),
                      ("agent_1", "mass", q(1716, "mg")), ("agent_1", "amount", q(19.5, "mmol")),
                      ("product", "mass", q(178, "mg")), ("product", "amount", q(0.44, "mmol"))]: set_charge(pget(r, pid), key, val, ev(10))
set_scale(r, "m_12ba", q(0.78, "mmol"), ev(10)); set_solvent(r, ["THF"], ev(5, 10))
r["steps"] = [new_step("reaction", "reaction", ["substrate", "agent_1", "agent_2"], ev(5, 7, 10), "prepared following the procedure for 8; heat and stir", q(85, "°C"), "overnight", inherit_values=True),
              new_step("workup_and_purification", "purification", ["product"], ev(7, 10), "Celite filtration, concentration and silica-gel chromatography following procedure for 8", inherit_values=True)]
r["outcomes"] = [outcome("yield", "yield", ["m_13"], q(56, "%"), ev(10), "main product", "product yield", "isolated after chromatography", True, {"start": "m_12ba", "end": "m_13"})]
set_procedure(r, "The title compound was prepared following the procedure for 8 starting from 12ba (320 mg, 0.78 mmol) and MnO2 (1716 mg, 19.5 mmol) to afford 13 in a 56 % yield (178 mg, 0.44 mmol).", ev(10), gp="procedure for 8")
classify(r, [("C", "C", 1, 2, "order_increased")], e, "Aromatizing oxidation is represented by a partial C=C bond-order increase.")
finalize(r, ["The page-10 paragraph does not restate the THF volume; exact product-specific solvent volume is unreported.", "Complete dehydrogenation mapping is unresolved."])

r = take_rxn("scheme5_13_to_14", 10); e = ev(5, 10)
for pid, key, val in [("substrate", "mass", q(170, "mg")), ("substrate", "amount", q(0.42, "mmol")),
                      ("agent_1", "volume", q(2.52, "mL")), ("agent_1", "amount", q(1.26, "mmol")),
                      ("agent_3", "volume", q(3.77, "mL")), ("agent_4", "volume", q(0.46, "mL")),
                      ("agent_5", "mass", q(90.1, "mg")), ("agent_5", "amount", q(0.84, "mmol")),
                      ("agent_6", "volume", q(5, "mL")), ("product", "mass", q(53.9, "mg")), ("product", "amount", q(0.134, "mmol"))]: set_charge(pget(r, pid), key, val, ev(10))
inter = add_participant(r, "intermediate_14a", "m_14a", "product", ev(10), product_role="crude intermediate consumed in next stage")
set_scale(r, "m_13", q(0.42, "mmol"), ev(10)); set_solvent(r, ["THF", "TFA"], ev(10))
r["steps"] = [new_step("hydroboration", "reaction", ["substrate", "agent_1", "agent_2"], ev(7, 10), "prepared following procedure for 9; add 9-BBN and hold until alkene consumption", None, q(3, "h"), inherit_values=True),
              new_step("oxidant_addition", "reaction", ["agent_3", "agent_4"], ev(7, 10), "cool then add NaOH and H2O2", q(0, "°C"), inherit_values=True),
              new_step("oxidation_hold", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4"], ev(7, 10), "stir", "room temperature", "overnight", inherit_values=True),
              new_step("isolate_crude_14a", "workup", ["intermediate_14a"], ev(10), "work-up following procedure for 9 to afford crude 14a", inherit_values=True),
              new_step("deprotection_esterification", "reaction", ["intermediate_14a", "agent_5", "agent_6"], ev(8, 10), "treat crude 14a with anisole in TFA", q(60, "°C"), "overnight", inherit_values=True),
              new_step("workup_purification", "purification", ["product"], ev(8, 10), "work-up and silica-gel chromatography following procedure for 9", inherit_values=True)]
r["outcomes"] = [outcome("overall_yield", "yield", ["m_14"], q(32, "%"), ev(10), "main product", "overall two-step yield", "isolated after chromatography", True, {"start": "m_13", "end": "m_14", "stages": 2})]
set_procedure(r, "13 (170 mg, 0.42 mmol) was treated with 9-BBN (0.5 M in THF, 2.52 mL, 1.26 mmol), followed by NaOH (1 M, 3.77 mL) and H2O2 (35 wt.-%/vol, 0.46 mL) to afford crude 14a, then anisole (90.1 mg, 0.84 mmol) in TFA (5 mL) afforded 14 in an overall 32 % yield.", ev(10), gp="procedure for 9")
classify(r, [("C", "O", 0, 1, "formed"), ("C", "H", 0, 1, "hydrogen_gain"), ("N", "H", 0, 1, "hydrogen_gain"), ("C", "O", 0, 1, "formed")], e,
         "Two-stage hydroboration/oxidation and deprotection/esterification are represented only by partial net classes.")
finalize(r, ["No amount or yield is reported for crude 14a.", "The 32% value is explicitly overall and is not a second-stage-only yield.", "The initial THF volume is not restated for this substrate."])

r = take_rxn("scheme5_14_to_15", 10); e = ev(5, 10)
for pid, key, val in [("substrate", "mass", q(25, "mg")), ("substrate", "amount", q(0.062, "mmol")),
                      ("agent_1", "mass", q(23, "mg")), ("agent_1", "amount", q(0.124, "mmol")),
                      ("product", "mass", q(26, "mg")), ("product", "amount", q(0.046, "mmol"))]: set_charge(pget(r, pid), key, val, ev(10))
set_scale(r, "m_14", q(0.062, "mmol"), ev(10)); set_solvent(r, ["THF"], ev(10))
r["steps"] = [new_step("reaction", "reaction", ["substrate", "agent_1", "agent_2"], ev(8, 10), "prepared following procedure for 10; add NBS and stir", "room temperature", "overnight", inherit_values=True),
              new_step("workup_purification", "purification", ["product"], ev(8, 10), "evaporate and purify by silica-gel chromatography following procedure for 10", inherit_values=True)]
r["outcomes"] = [outcome("yield", "yield", ["m_15"], q(75, "%"), ev(10), "main product", "product yield", "isolated after chromatography", True, {"start": "m_14", "end": "m_15"})]
set_procedure(r, "The title compound was prepared following the procedure for 10 starting from 14 (25 mg, 0.062 mmol) and NBS (23 mg, 0.124 mmol) to afford 15 in a 75 % yield (26 mg, 0.046 mmol).", ev(10), gp="procedure for 10")
classify(r, [("C", "Br", 0, 1, "formed"), ("C", "Br", 0, 1, "formed")], e, "Two brominations are retained as partial bond classes without atom mapping.")
finalize(r, ["THF volumes are not restated for this substrate-specific experiment."])

r = take_rxn("scheme5_15_to_16", 10); e = ev(5, 10)
for pid, key, val in [("substrate", "mass", q(15, "mg")), ("substrate", "amount", q(0.029, "mmol")),
                      ("agent_1", "volume", q(0.1, "mL")), ("product", "mass", q(13.2, "mg")), ("product", "amount", q(0.029, "mmol"))]: set_charge(pget(r, pid), key, val, ev(10))
set_scale(r, "m_15", q(0.029, "mmol"), ev(10)); set_solvent(r, ["methanol", "THF", "water"], ev(8, 10), ratio="1:2:1")
r["steps"] = [new_step("reaction", "reaction", ["substrate", "agent_1", "agent_2", "agent_3"], ev(8, 10), "prepared following procedure for 11; add 1 M NaOH dropwise and hold", "room temperature", q(3, "h"), inherit_values=True),
              new_step("workup", "workup", ["product"], ev(8, 10), "remove methanol/THF, extract, dry, filter and concentrate following procedure for 11", inherit_values=True)]
r["outcomes"] = [outcome("yield", "yield", ["m_16"], q(99, "%"), ev(10), "main product", "product yield", "isolated after work-up", True, {"start": "m_15", "end": "m_16"})]
set_procedure(r, "The title compound was prepared following the procedure for 11 starting from 15 (15 mg, 0.029 mmol) and NaOH (0.1 mL, 1 M) to afford 16 in a 99 % yield (13.2 mg, 0.029 mmol).", ev(10), gp="procedure for 11")
classify(r, [("O", "H", 0, 1, "hydrogen_gain")], e, "Trifluoroacetate hydrolysis is represented by net O-H gain without leaving-group mapping.")
finalize(r, ["The total methanol/THF/water volume is not restated for this substrate-specific experiment."])

r = take_rxn("scheme5_16_to_17", 10, 11); e = ev(5, 8, 10, 11)
for pid, key, val in [("substrate", "mass", q(4.00, "mg")), ("substrate", "amount", q(0.009, "mmol")),
                      ("agent_2", "mass", q(7.24, "mg")), ("agent_2", "amount", q(0.022, "mmol")),
                      ("agent_3", "mass", q(3.00, "mg")), ("agent_3", "amount", q(0.036, "mmol")),
                      ("agent_1", "mass", q(0.15, "mg")), ("agent_1", "amount", q(0.001, "mmol")),
                      ("product", "mass", q(3.56, "mg")), ("product", "amount", q(0.007, "mmol"))]: set_charge(pget(r, pid), key, val, ev(10))
add_agent_participant(r, "silica", "agent_silica_gel", "purification_agent", ev(10), "purification")
set_scale(r, "m_16", q(0.009, "mmol"), ev(10)); set_solvent(r, ["MeCN", "water"], ev(8, 10), ratio="1:1")
r["steps"] = [new_step("oxidation", "reaction", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4", "agent_5"], ev(8, 10), "prepared following the procedure for Longamide B; stir", "room temperature", q(18, "h"), inherit_values=True),
              new_step("workup", "workup", ["substrate", "agent_1", "agent_2", "agent_3", "agent_4", "agent_5"], ev(8, 10), "bicarbonate work-up, acidification and extraction following Longamide B procedure", inherit_values=True),
              new_step("purification", "purification", ["silica", "product"], ev(10), "column chromatography")]
r["outcomes"] = [outcome("yield", "yield", ["m_17"], q(83, "%"), ev(10), "inseparable diastereomeric product mixture", "product yield", "isolated after column chromatography", True, {"start": "m_16", "end": "m_17"}),
                 outcome("diastereomer_ratio", "selectivity", ["m_17"], {"metric": "diastereomer_ratio", "order": ["unassigned diastereomer 1", "unassigned diastereomer 2"], "ratio": [44, 56]}, ev(10), "inseparable diastereomeric mixture", basis="1H NMR mixture integration", scope={"metric": "dr", "order": ["unassigned diastereomer 1", "unassigned diastereomer 2"], "configurations": "unassigned"})]
set_procedure(r, "The title compound was prepared following the procedure for Longamide B starting from 16 (4.00 mg, 0.009 mmol), BAIB (7.24 mg, 0.022 mmol), NaHCO3 (3.00 mg, 0.036 mmol) and TEMPO (0.15 mg, 0.001 mmol); column chromatography afforded 17 as an inseparable mixture of diastereoisomers (44%:56%) in 83% yield.", ev(10), gp="procedure for Longamide B")
classify(r, [("C", "O", 1, 2, "order_increased"), ("C", "O", 0, 1, "formed"), ("O", "H", 0, 1, "hydrogen_gain")], e,
         "Alcohol oxidation is represented by partial carboxylic-acid bond changes; epimerization at the alpha center is not encoded as a bond change.")
finalize(r, ["The two product diastereomers are not assigned absolute configurations and remain one explicitly inseparable mixture identity.", "The product-specific MeCN/water volume is not restated."])


# Add every molecule referenced by the incoming records, preserving the v3
# independent registry contract.
core_reaction_ids = [r["reaction_id"] for r in reactions]

# Partial-bundle validation requires referenced generic schemes to travel with
# the concrete records.  They are dependencies, not additional page-7--11
# experimental candidates.
dependency_ids = sorted({sid for r in reactions for sid in r.get("general_scheme_ids", [])})
for rid in dependency_ids:
    dependency = copy.deepcopy(old_rxns[rid])
    dependency["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                            "checked_fields": [], "unresolved": dependency.get("review", {}).get("unresolved", [])}
    reactions.append(dependency)

for r in reactions:
    for p in r["participants"]:
        if p["molecule_id"] not in molmap:
            take_mol(p["molecule_id"])

bundle["molecules"] = list(molmap.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {"context_id": f"{DOC}:page7_characterization", "document_id": DOC, "evidence_ids": ev(7), "kind": "characterization",
     "summary": "NMR/IR/HRMS/optical-rotation/TLC data for 6b, 7a, 7b, 7c, 7d, 8 and 8a accompany preparations but are not additional synthetic experiments."},
    {"context_id": f"{DOC}:page8_characterization", "document_id": DOC, "evidence_ids": ev(8), "kind": "characterization",
     "summary": "Spectral/TLC/optical-rotation blocks for 9a, 9, 10, 11, Longamide B, Stylisine D and 12aa are retained as characterization, not separate reactions."},
    {"context_id": f"{DOC}:page9_characterization", "document_id": DOC, "evidence_ids": ev(9), "kind": "characterization",
     "summary": "Spectral/TLC/optical-rotation blocks for 12ba and Table 3 products accompany the scope preparations; no mass-spectrum peak table is treated as synthesis."},
    {"context_id": f"{DOC}:page10_characterization", "document_id": DOC, "evidence_ids": ev(10), "kind": "characterization",
     "summary": "Spectral/TLC/optical-rotation blocks for 12ae, 12be and 13--17 accompany their preparations and are not additional synthesis records."},
    {"context_id": f"{DOC}:page11_back_matter", "document_id": DOC, "evidence_ids": ev(11), "kind": "characterization_and_back_matter",
     "summary": "Page 11 begins with the IR/HRMS tail for product 17, followed by acknowledgments, keywords, references and received date; no new reaction execution occurs."},
]

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/pages7_11_candidates.json").write_text(json.dumps({
    "document_id": DOC,
    "page_range": [7, 11],
    "candidate_reaction_ids": core_reaction_ids,
    "new_reaction_ids": [f"{DOC}:page7_7d_to_8a"],
    "note": "Machine candidates generated from visual/text review; all remain needs_review."
}, ensure_ascii=False, indent=2) + "\n")

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 7, "state": "attempted", "regions": [
            {"region": "left column top and compound 6b continuation", "finding": "Completion of the 5b->6b deprotection paragraph with substrate/TFA/product quantities and 92% reported yield; spectral block is characterization.", "evidence_ids": ev(7)},
            {"region": "left column middle, 7a/7b preparation", "finding": "Concrete AgNO3/CaCO3 cyclisation of 6b with separate isolated 7a and 7b yields, work-up, heat-sensitivity warning and Table 1 entry-r correspondence.", "evidence_ids": ev(7)},
            {"region": "left column bottom to right column top, 7d/7c preparation and footnote", "finding": "Concrete inherited cyclisation of 6a gives 7d and 7c separately; footnote additionally reports actual oxidation of unstable 7d to the other enantiomer 8a without scale/yield.", "evidence_ids": ev(7)},
            {"region": "right column middle, compound 8 preparation", "finding": "Concrete MnO2 oxidation of 7c with ordered heating, two Celite filtrations/reflux and chromatography.", "evidence_ids": ev(7)},
            {"region": "right column bottom, start of 8->9a->9", "finding": "First stage hydroboration/oxidation through crude 9a and start of TFA/anisole stage; paragraph continues on page 8.", "evidence_ids": ev(7)}],
         "record_ids": [f"{DOC}:scheme2_5b_to_6b", f"{DOC}:table1_r", f"{DOC}:scheme3_6a_to_7c", f"{DOC}:scheme3_7c_to_8", f"{DOC}:page7_7d_to_8a", f"{DOC}:scheme3_8_to_9"],
         "unresolved": ["7d->8a scale and yield are unreported.", "Several inherited procedure volumes are not product-specifically restated.", "Stereo-safe SMILES for fused products were not asserted from experimental-page text alone."]},
        {"page": 8, "state": "attempted", "regions": [
            {"region": "left column top, continuation of compound 9", "finding": "Completion of TFA/anisole stage, work-up and explicit 55% overall two-step yield from 8 to 9; 9a is only a contaminated crude mass.", "evidence_ids": ev(8)},
            {"region": "left column middle and bottom", "finding": "Concrete preparations 9->10, 10->11 and start of 11->Longamide B with complete stated charges and work-up sequence.", "evidence_ids": ev(8)},
            {"region": "right column top and middle", "finding": "Completion of Longamide B isolation/recrystallization and concrete Longamide B->Stylisine D amidation.", "evidence_ids": ev(8)},
            {"region": "right column bottom", "finding": "Product-specific scale-up/scope execution 6a->12aa+12ba with nitrogen, 18 h reflux and separate isolated yields; linked to Table 3a, not duplicated.", "evidence_ids": ev(8)},
            {"region": "spectral and TLC blocks throughout", "finding": "Characterization accompanying preparations; no independent synthesis records created from NMR/IR/HRMS values.", "evidence_ids": ev(8)}],
         "record_ids": [f"{DOC}:scheme3_8_to_9", f"{DOC}:scheme3_9_to_10", f"{DOC}:scheme3_10_to_11", f"{DOC}:scheme3_11_to_longamide_b", f"{DOC}:scheme3_longamide_b_to_stylisine_d", f"{DOC}:table3_a"],
         "unresolved": ["The 9a intermediate has no isolated yield and is explicitly contaminated.", "Table 3a footnote and product paragraph disagree on scale and MeCN volume; both are retained as conflicts."]},
        {"page": 9, "state": "attempted", "regions": [
            {"region": "left column top", "finding": "12ba characterization continues from page 8; linked to Table 3a but not counted as a new experiment.", "evidence_ids": ev(9)},
            {"region": "left and right columns, product headings 12ab/12bb through 12ad/12bd", "finding": "Three concrete Table 3 scope preparations b-d with product-specific charges and separate isolated diastereomer yields; inherited procedure retained.", "evidence_ids": ev(9)},
            {"region": "right column bottom", "finding": "Start of Table 3e preparation from 6b and 2-iodonaphthalene; product sentence and characterization continue on page 10.", "evidence_ids": ev(9)},
            {"region": "NMR/IR/HRMS/TLC/optical rotation blocks", "finding": "Characterization only; not separate reactions.", "evidence_ids": ev(9)}],
         "record_ids": [f"{DOC}:table3_a", f"{DOC}:table3_b", f"{DOC}:table3_c", f"{DOC}:table3_d", f"{DOC}:table3_e"],
         "unresolved": ["K2CO3 is printed as '(29.0, 0.210 mmol)' without a mass unit in scope paragraphs; no mass unit was invented.", "Shared Table 3 footnote uses 0.1 mmol/0.23 mmol whereas product paragraphs use 0.105 mmol/0.210 mmol; retained as conflicts."]},
        {"page": 10, "state": "attempted", "regions": [
            {"region": "left column top", "finding": "Completion and characterization of Table 3e products 12ae/12be; same experiment as page 9.", "evidence_ids": ev(10)},
            {"region": "left column middle and bottom", "finding": "Concrete 12ba->13 oxidation and two-stage 13->14a->14 hydroboration/oxidation/deprotection sequence; 32% explicitly overall.", "evidence_ids": ev(10)},
            {"region": "right column middle", "finding": "Concrete inherited preparations 14->15 and 15->16 with product-specific charges and yields.", "evidence_ids": ev(10)},
            {"region": "right column bottom", "finding": "Concrete 16->17 oxidation with separate 83% yield and inseparable 44:56 unassigned diastereomer ratio; characterization continues on page 11.", "evidence_ids": ev(10)},
            {"region": "spectral/TLC/optical rotation blocks", "finding": "Characterization accompanying preparations; not additional synthesis experiments.", "evidence_ids": ev(10)}],
         "record_ids": [f"{DOC}:table3_e", f"{DOC}:scheme5_12ba_to_13", f"{DOC}:scheme5_13_to_14", f"{DOC}:scheme5_14_to_15", f"{DOC}:scheme5_15_to_16", f"{DOC}:scheme5_16_to_17"],
         "unresolved": ["Crude 14a has no amount or yield.", "The 17 diastereomers have no absolute configuration assignments.", "Inherited solvent volumes are not restated for several product-specific procedures."]},
        {"page": 11, "state": "attempted", "regions": [
            {"region": "top left first two lines", "finding": "Only the IR/HRMS tail of product 17 characterization; linked to the page-10 reaction, not a new synthesis.", "evidence_ids": ev(11)},
            {"region": "remaining left column and right column", "finding": "Acknowledgments, keywords, references 1-57 and received date; cited literature is not reclassified as current-paper experiments.", "evidence_ids": ev(11)}],
         "record_ids": [f"{DOC}:scheme5_16_to_17"],
         "unresolved": ["No new reaction execution occurs on this page; bibliography was not imported as experiments."]},
    ]
}
(RUN / "coverage/pages_0007_0011.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/page_visual_inventory_pages7_11.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")

merged_dataset_path = RUN / "dataset.json"
if merged_dataset_path.exists():
    merged_dataset = json.loads(merged_dataset_path.read_text())
    if any(r.get("record_id") == f"{DOC}:page7_7d_to_8a" for r in merged_dataset.get("experiments", [])):
        audit_result = rt.dataset_audit(merged_dataset)
        (RUN / "raw_model/audit_pages7_11.json").write_text(
            json.dumps(audit_result, ensure_ascii=False, indent=2) + "\n"
        )

print(json.dumps({"molecules_in_bundle": len(bundle["molecules"]), "reactions_in_bundle": len(reactions),
                  "new_reactions": 1, "updated_core_reactions": len(core_reaction_ids)-1,
                  "generic_dependencies": len(dependency_ids)}, ensure_ascii=False))
