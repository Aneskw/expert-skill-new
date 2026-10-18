#!/usr/bin/env python3
"""Build the page 7--11 v3.1 extraction increment for doc_2992316b1e64fc22.

The script deliberately updates the reaction IDs already established from the
schemes.  It adds the experimental-section quantities and ordered operations
without creating duplicate reactions.  Biological assays and characterization
blocks are retained as source_context rather than synthetic experiments.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DOC = "doc_2992316b1e64fc22"
RAW = ROOT / "runs/current/raw_model"

spec = importlib.util.spec_from_file_location(
    "standard_records", ROOT / "skills/chemical-reaction-extraction/scripts/standard_records.py"
)
st = importlib.util.module_from_spec(spec)
assert spec.loader
spec.loader.exec_module(st)

init = json.loads((RAW / "pages_0007_0011_init.json").read_text())
old = json.loads((ROOT / "runs/current/standard_records.json").read_text())
old_reactions = {r["reaction_id"]: r for r in old["reactions"]}
old_molecules = {m["molecule_id"]: m for m in old["molecules"]}


def eid(page: int, kind: str = "text") -> str:
    return f"{DOC}:page{page}:{kind}"


P9 = [eid(9, "image"), eid(9)]
P10 = [eid(10, "image"), eid(10)]


def q(raw, number, unit, ev):
    return st.fact(raw=raw, normalized={"number": number, "unit": unit}, evidence=ev)


def range_q(raw, lower, upper, unit, ev):
    return st.fact(raw=raw, normalized={"minimum": lower, "maximum": upper, "unit": unit}, evidence=ev)


def text_fact(raw, ev, normalized=None):
    return st.fact(raw=raw, normalized=raw if normalized is None else normalized, evidence=ev)


def missing(ev, status="unreported"):
    return st.fact(evidence=ev, status=status, verification="missing")


def conflict(raw, alternatives, ev):
    return st.fact(raw=raw, normalized=None, evidence=ev, status="conflict", alternatives=alternatives)


def add_images_to_facts(obj):
    """Pair page 9/10 extracted text evidence with the images actually inspected."""
    if isinstance(obj, dict):
        if "normalized_value" in obj and "evidence_ids" in obj:
            for page in (9, 10):
                tx, im = eid(page), eid(page, "image")
                if tx in obj["evidence_ids"] and im not in obj["evidence_ids"]:
                    obj["evidence_ids"].append(im)
        for value in obj.values():
            add_images_to_facts(value)
    elif isinstance(obj, list):
        for value in obj:
            add_images_to_facts(value)


def molid(suffix):
    return f"{DOC}:{suffix}"


def named_molecule(suffix, name, ev=P9, identity_status="source_identified"):
    mid = molid(suffix)
    m = st.new_molecule(mid, DOC, ev)
    m["identity_status"] = identity_status
    m["occurrences"][0]["compound_label"] = text_fact(name, ev)
    m["names"]["original"] = text_fact(name, ev)
    m["names"]["english"] = text_fact(name, ev)
    m["cas_numbers"] = missing(ev)
    m["representation_type"] = "unresolved"
    m["smiles"] = missing(ev)
    for key in ("atom_centers", "double_bonds", "axial", "mixture", "assignment", "isomeric_smiles", "identity_key"):
        m["stereochemistry"][key] = missing(ev)
    derivation = "The source names this procedural component, but no complete structure or stereochemical representation was asserted."
    for key in ("chirality_status", "geometry_status"):
        m["stereochemistry"][key] = st.fact(
            raw="unresolved", normalized="unresolved", evidence=ev,
            method="inference", verification="inferred", derivation=derivation,
        )
    m["stereochemistry"]["representation_status"] = st.fact(
        raw="not encoded", normalized="not_encoded", evidence=ev,
        method="inference", verification="inferred", derivation=derivation,
    )
    m["identity_note"] = text_fact(
        "Named source component; CAS and structure were not supplied or externally completed.", ev
    )
    return m


extra_molecules = {
    m["molecule_id"]: m
    for m in [
        named_molecule("reagent_nabh4", "sodium borohydride"),
        named_molecule("reagent_eda", "ethylenediamine"),
        named_molecule("reagent_etoac", "ethyl acetate"),
        named_molecule("agent_celite", "Celite"),
        named_molecule("reagent_hexane", "n-hexane"),
        named_molecule("reagent_triethylamine", "triethylamine"),
        named_molecule("reagent_hcl", "hydrochloric acid"),
        named_molecule("reagent_water", "water"),
        named_molecule("reagent_na2so4", "anhydrous sodium sulfate"),
        named_molecule("reagent_nh4cl", "ammonium chloride"),
        named_molecule("reagent_brine", "brine"),
        named_molecule("reagent_mgso4", "anhydrous magnesium sulfate"),
        named_molecule("reagent_nahco3", "sodium bicarbonate"),
        named_molecule("reagent_argon", "argon"),
        named_molecule("agent_silica_gel", "silica gel"),
    ]
}


def amine_conflict_molecule():
    ev = [eid(2, "image"), eid(3), eid(3, "image"), *P9]
    m = named_molecule("reagent_amine_9_to_10_conflict", "amine used for 9→10", ev, "unresolved")
    alt1 = st.fact(
        raw="HN(C6H11)2 / dicyclohexylamine", normalized="dicyclohexylamine",
        evidence=[eid(2, "image"), *P9],
    )
    alt2 = st.fact(raw="dihexylamine", normalized="dihexylamine", evidence=[eid(3), eid(3, "image")])
    cf = conflict(
        ["HN(C6H11)2 / dicyclohexylamine (scheme and experimental section)", "dihexylamine (discussion)"],
        [alt1, alt2], ev,
    )
    m["occurrences"][0]["compound_label"] = copy.deepcopy(cf)
    m["names"]["original"] = copy.deepcopy(cf)
    m["names"]["english"] = copy.deepcopy(cf)
    m["identity_note"] = text_fact(
        "The source is internally inconsistent: the drawn formula and Experimental Section indicate dicyclohexylamine, while the discussion says dihexylamine; no identity was selected.", ev
    )
    return m


amine_mol = amine_conflict_molecule()
extra_molecules[amine_mol["molecule_id"]] = amine_mol


def base_reaction(suffix, page):
    rid = f"{DOC}:{suffix}"
    r = copy.deepcopy(old_reactions[rid])
    pe = [eid(page, "image"), eid(page)]
    for x in pe:
        if x not in r["evidence_ids"]:
            r["evidence_ids"].append(x)
    add_images_to_facts(r)
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["review"]["reviewer"] = None
    r["review"]["checked_fields"] = []
    r["completeness"] = "partial"
    # Product-amount outcomes are regenerated below so rerunning this builder
    # after a merge is idempotent.
    r["outcomes"] = [o for o in r["outcomes"] if o["kind"] != "product_amount"]
    return r


def participant(pid, mid, role, ev, stage="reaction", product_role=None):
    p = st.new_participant(pid, mid, role, ev, stage)
    if product_role is not None:
        p["product_role"] = text_fact(product_role, ev)
    return p


def set_charge(p, kind, raw, number, unit, ev):
    p["charges"][kind] = q(raw, number, unit, ev)


def step(step_id, stage, pids, operation, ev, temperature=None, time=None, atmosphere=None, pressure=None):
    s = st.new_step(step_id, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = text_fact(operation, ev)
    if temperature is not None:
        s["temperature"] = temperature if isinstance(temperature, dict) else text_fact(temperature, ev)
    if time is not None:
        s["time"] = time if isinstance(time, dict) else text_fact(time, ev)
    if atmosphere is not None:
        s["atmosphere"] = text_fact(atmosphere, ev)
    if pressure is not None:
        s["pressure"] = text_fact(pressure, ev)
    return s


def outcome(outcome_id, kind, target, raw, normalized, ev, role="main product", basis=None, yield_type=None):
    o = st.new_outcome(outcome_id, kind, [target], ev)
    o["product_role"] = text_fact(role, ev)
    o["value"] = st.fact(raw=raw, normalized=normalized, evidence=ev)
    if basis is not None:
        o["measurement_basis"] = text_fact(basis, ev)
    if yield_type is not None:
        o["yield_type"] = text_fact(yield_type, ev)
    return o


def product_amount(oid, target, raw, number, unit, ev, basis):
    return outcome(oid, "product_amount", target, raw, {"number": number, "unit": unit}, ev, basis=basis)


def reaction_scale(mid, raw, number, unit, ev):
    value = {"reference_molecule_id": mid, "amount": {"number": number, "unit": unit}}
    return st.fact(raw={"reference_molecule_id": mid, "amount": raw}, normalized=value, evidence=ev)


def add_unresolved(r, *messages):
    for message in messages:
        if message not in r["review"]["unresolved"]:
            r["review"]["unresolved"].append(message)


reactions = []

# 4/3 -> 5: page 9 reports the one-pot experiment beginning from compound 3.
r = base_reaction("scheme2_4_to_5", 9)
p3 = participant("s3_experimental", molid("mol_3"), "reactant", P9)
set_charge(p3, "mass", "20 g", 20, "g", P9); set_charge(p3, "amount", "80 mmol", 80, "mmol", P9)
etmg = participant("r_etmgbr", molid("reagent_etmgbr"), "reagent", P9)
set_charge(etmg, "volume", "32.0 mL", 32.0, "mL", P9); set_charge(etmg, "amount", "96.0 mmol", 96.0, "mmol", P9)
thf1 = participant("sol_thf_initial", molid("reagent_thf"), "solvent", P9); set_charge(thf1, "volume", "50 mL", 50, "mL", P9)
p4 = participant("in_situ_4", molid("mol_4"), "unknown", P9, product_role="in-situ organomagnesium intermediate")
eo = participant("r_ethylene_oxide", molid("reagent_ethylene_oxide"), "reagent", P9)
set_charge(eo, "mass", "5.85 g", 5.85, "g", P9); set_charge(eo, "amount", "130 mmol", 130, "mmol", P9)
thf2 = participant("sol_thf_ethylene_oxide", molid("reagent_thf"), "solvent", P9); set_charge(thf2, "volume", "30 mL", 30, "mL", P9)
hcl = participant("workup_hcl", molid("reagent_hcl"), "workup_agent", P9, "workup"); set_charge(hcl, "volume", "50 mL of 5% HCl", 50, "mL", P9)
etoac = participant("extract_etoac", molid("reagent_etoac"), "workup_agent", P9, "workup"); set_charge(etoac, "volume", "2 × 100 mL", 200, "mL", P9)
mgso4 = participant("dry_mgso4", molid("reagent_mgso4"), "workup_agent", P9, "workup")
p5 = participant("p5", molid("mol_5"), "product", P9, product_role="main product")
r["participants"] = [p3, etmg, thf1, p4, eo, thf2, hcl, etoac, mgso4, p5]
r["reaction_scale"] = reaction_scale(molid("mol_3"), "80 mmol compound 3", 80, "mmol", P9)
r["steps"] = [
    step("add_etmgbr", "reaction", ["s3_experimental", "r_etmgbr", "sol_thf_initial"], "slowly add EtMgBr solution to octadec-1-yne in THF", P9),
    step("reflux_grignard", "reaction", ["s3_experimental", "r_etmgbr", "sol_thf_initial"], "heat under reflux", P9, "reflux", q("4 h", 4, "h", P9)),
    step("cool_over_ice", "reaction", ["in_situ_4"], "cool over ice", P9, "ice cooling"),
    step("add_ethylene_oxide", "reaction", ["in_situ_4", "r_ethylene_oxide", "sol_thf_ethylene_oxide"], "slowly add ice-cold ethylene oxide solution; slight exotherm reported", P9, "ice cooled"),
    step("warm_and_stir_overnight", "reaction", ["in_situ_4", "r_ethylene_oxide"], "allow to warm and stir", P9, "room temperature", "overnight"),
    step("acid_hydrolysis", "workup", ["workup_hcl"], "cool over ice and carefully add 5% HCl", P9, "ice cooled"),
    step("extract_and_dry", "workup", ["extract_etoac", "dry_mgso4"], "extract twice with EtOAc, dry over MgSO4 and evaporate", P9),
]
existing_yield = copy.deepcopy(r["outcomes"][0])
existing_yield["measurement_basis"] = text_fact(
    "Scheme 2 gives 78% without a stated basis; the Experimental Section gives 16.0 g and 68% after extraction, drying and evaporation.",
    [eid(2, "image"), *P9],
)
r["outcomes"] = [existing_yield, product_amount("amount_5", molid("mol_5"), "16.0 g", 16.0, "g", P9, "material after extraction, drying and evaporation")]
r["procedure"]["original_text"] = text_fact("Octadec-1-yne/EtMgBr in THF was refluxed 4 h; ethylene oxide was added cold, the mixture warmed and stirred overnight, then hydrolyzed with 5% HCl and extracted.", P9)
r["procedure"]["addition_order"] = text_fact("EtMgBr to 3; after reflux and cooling, ethylene oxide solution; after overnight reaction and cooling, 5% HCl.", P9)
add_unresolved(r, "The page-9 experiment is one-pot from compound 3, whereas the scheme separates transient 4→5.", "The Scheme 2 yield (78%) and Experimental Section yield (68%, 16.0 g) conflict and remain unselected.")
reactions.append(r)

# 5 -> 6 hydrogenation.
r = base_reaction("scheme2_5_to_6", 9)
p5 = participant("s5", molid("mol_5"), "reactant", P9); set_charge(p5, "mass", "12.6 g", 12.6, "g", P9); set_charge(p5, "amount", "43 mmol", 43, "mmol", P9)
ni = participant("cat_ni", molid("reagent_ni_oac2"), "catalyst", P9); set_charge(ni, "mass", "10.7 g", 10.7, "g", P9); set_charge(ni, "amount", "42.7 mmol", 42.7, "mmol", P9)
nabh4 = participant("r_nabh4", molid("reagent_nabh4"), "reagent", P9); set_charge(nabh4, "mass", "1.2 g", 1.2, "g", P9); set_charge(nabh4, "amount", "31 mmol", 31, "mmol", P9)
etoh1 = participant("sol_etoh_ni", molid("reagent_etoh"), "solvent", P9); set_charge(etoh1, "volume", "86 mL", 86, "mL", P9)
etoh2 = participant("sol_etoh_nabh4", molid("reagent_etoh"), "solvent", P9); set_charge(etoh2, "volume", "43 mL", 43, "mL", P9)
h2 = participant("r_h2", molid("reagent_h2"), "reagent", P9)
eda = participant("r_eda", molid("reagent_eda"), "additive", P9); set_charge(eda, "mass", "10.32 g", 10.32, "g", P9); set_charge(eda, "amount", "171.5 mmol", 171.5, "mmol", P9)
etoac = participant("dilute_etoac", molid("reagent_etoac"), "workup_agent", P9, "workup")
celite = participant("filter_celite", molid("agent_celite"), "purification_agent", P9, "purification")
silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P9, "purification")
hexane = participant("eluent_hexane", molid("reagent_hexane"), "solvent", P9, "purification")
etoac_el = participant("eluent_etoac", molid("reagent_etoac"), "solvent", P9, "purification")
p6 = participant("p6", molid("mol_6"), "product", P9, product_role="main product")
r["participants"] = [p5, ni, nabh4, etoh1, etoh2, h2, eda, etoac, celite, silica, hexane, etoac_el, p6]
r["reaction_scale"] = reaction_scale(molid("mol_5"), "43 mmol compound 5", 43, "mmol", P9)
hold_time = conflict(
    ["24 h (Scheme 2)", "12 h (Experimental Section)"],
    [q("24 h", 24, "h", [eid(2, "image")]), q("12 h", 12, "h", P9)],
    [eid(2, "image"), *P9],
)
r["steps"] = [
    step("prepare_nickel_boron", "reaction", ["cat_ni", "r_nabh4", "sol_etoh_ni", "sol_etoh_nabh4", "r_h2"], "add ethanolic NaBH4 to Ni(OAc)2·4H2O in EtOH and stir under hydrogen", P9, "room temperature", q("60 min", 60, "min", P9), "hydrogen"),
    step("add_eda_then_5", "reaction", ["r_eda", "s5"], "add ethylenediamine and alkynol 5 sequentially", P9),
    step("hydrogenation_hold", "reaction", ["s5", "cat_ni", "r_h2", "r_eda"], "stir at the same temperature under hydrogen", [eid(2, "image"), *P9], "room temperature", hold_time, "hydrogen"),
    step("concentrate_filter", "workup", ["dilute_etoac", "filter_celite"], "concentrate, dilute with EtOAc, filter through Celite and reconcentrate", P9),
    step("silica_purification", "purification", ["purify_silica", "eluent_hexane", "eluent_etoac"], "silica gel chromatography with n-hexane/EtOAc 5:1", P9),
]
r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P9)
r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P9)
r["outcomes"].append(product_amount("amount_6", molid("mol_6"), "12.24 g", 12.24, "g", P9, "product isolated after silica gel chromatography"))
r["procedure"]["original_text"] = text_fact("Ni(OAc)2·4H2O/NaBH4 was stirred 60 min under H2; ethylenediamine and 5 were added sequentially and the mixture stirred 12 h, then filtered and chromatographed.", P9)
add_unresolved(r, "Hydrogenation duration conflicts between Scheme 2 (24 h) and Experimental Section (12 h).")
reactions.append(r)

# 6 -> 7 mesylation.
r = base_reaction("scheme2_6_to_7", 9)
p6 = participant("s6", molid("mol_6"), "reactant", P9); set_charge(p6, "mass", "10.0 g", 10.0, "g", P9); set_charge(p6, "amount", "34.0 mmol", 34.0, "mmol", P9)
mscl = participant("r_mscl", molid("reagent_mscl"), "reagent", P9); set_charge(mscl, "mass", "4.93 g", 4.93, "g", P9); set_charge(mscl, "amount", "43.0 mmol", 43.0, "mmol", P9)
tea = participant("r_triethylamine", molid("reagent_triethylamine"), "reagent", P9); set_charge(tea, "volume", "8.22 mL", 8.22, "mL", P9)
dcm = participant("sol_dcm", molid("reagent_dcm"), "solvent", P9); set_charge(dcm, "volume", "52.2 mL", 52.2, "mL", P9)
dcmwash = participant("wash_dcm", molid("reagent_dcm"), "workup_agent", P9, "workup"); set_charge(dcmwash, "volume", "40 mL", 40, "mL", P9)
hcl = participant("wash_hcl", molid("reagent_hcl"), "workup_agent", P9, "workup"); set_charge(hcl, "volume", "30 mL diluted HCl", 30, "mL", P9)
water = participant("wash_water", molid("reagent_water"), "workup_agent", P9, "workup"); set_charge(water, "volume", "2 × 25 mL", 50, "mL", P9)
na2so4 = participant("dry_na2so4", molid("reagent_na2so4"), "workup_agent", P9, "workup")
p7 = participant("p7", molid("mol_7"), "product", P9, product_role="main product")
r["participants"] = [p6, mscl, tea, dcm, dcmwash, hcl, water, na2so4, p7]
r["reaction_scale"] = reaction_scale(molid("mol_6"), "34.0 mmol compound 6", 34.0, "mmol", P9)
temp = conflict(["room temperature (Scheme 2)", "0 °C (Experimental Section)"], [text_fact("room temperature", [eid(2, "image")]), q("0 °C", 0, "°C", P9)], [eid(2, "image"), *P9])
tm = conflict(["4 h (discussion)", "1 h (Experimental Section)"], [q("4 h", 4, "h", [eid(3), eid(3, "image")]), q("1 h", 1, "h", P9)], [eid(3), eid(3, "image"), *P9])
r["steps"] = [
    step("mesylation", "reaction", ["s6", "r_mscl", "r_triethylamine", "sol_dcm"], "stir the mixture", [eid(2, "image"), eid(3), eid(3, "image"), *P9], temp, tm),
    step("filter_and_wash", "workup", ["wash_dcm"], "filter and wash the residue with CH2Cl2", P9),
    step("aqueous_washes", "workup", ["wash_hcl", "wash_water", "dry_na2so4"], "wash combined CH2Cl2 layer with diluted HCl and twice with water; dry over Na2SO4", P9),
    step("remove_solvent", "purification", ["p7"], "remove solvent under reduced pressure", P9),
]
r["outcomes"][0]["measurement_basis"] = text_fact("mass of residue after workup and solvent removal", P9)
r["outcomes"].append(product_amount("amount_7", molid("mol_7"), "11.88 g", 11.88, "g", P9, "residue after workup and solvent removal"))
add_unresolved(r, "Mesylation temperature and duration conflict across the scheme/discussion and the Experimental Section; alternatives remain unselected.")
reactions.append(r)

# 7 -> 8 bromide substitution.
r = base_reaction("scheme2_7_to_8", 9)
p7 = participant("s7", molid("mol_7"), "reactant", P9); set_charge(p7, "mass", "11.88 g", 11.88, "g", P9)
libr = participant("r_libr", molid("reagent_libr"), "reagent", P9); set_charge(libr, "mass", "6.0 g", 6.0, "g", P9); set_charge(libr, "amount", "69.2 mmol", 69.2, "mmol", P9)
acetone = participant("sol_acetone", molid("reagent_acetone"), "solvent", P9); set_charge(acetone, "volume", "104 mL", 104, "mL", P9)
water = participant("quench_water", molid("reagent_water"), "workup_agent", P9, "workup"); set_charge(water, "volume", "60 mL ice-cold water", 60, "mL", P9)
hexane = participant("extract_hexane", molid("reagent_hexane"), "workup_agent", P9, "workup"); set_charge(hexane, "volume", "2 × 50 mL", 100, "mL", P9)
waterwash = participant("wash_water", molid("reagent_water"), "workup_agent", P9, "workup"); set_charge(waterwash, "volume", "20 mL", 20, "mL", P9)
na2so4 = participant("dry_na2so4", molid("reagent_na2so4"), "workup_agent", P9, "workup")
silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P9, "purification")
ether = participant("eluent_ether", molid("reagent_ether"), "solvent", P9, "purification")
p8 = participant("p8", molid("mol_8"), "product", P9, product_role="main product")
r["participants"] = [p7, libr, acetone, water, hexane, waterwash, na2so4, silica, ether, p8]
r["steps"] = [
    step("reflux_with_libr", "reaction", ["s7", "r_libr", "sol_acetone"], "add LiBr to the acetone solution and reflux", P9, "reflux", q("3 h", 3, "h", P9)),
    step("add_ice_water", "workup", ["quench_water"], "add ice-cold water and separate the oil", P9),
    step("extract_wash_dry", "workup", ["extract_hexane", "wash_water", "dry_na2so4"], "extract aqueous layer twice with n-hexane; combine with oil, wash with water and dry", P9),
    step("chromatography", "purification", ["purify_silica", "extract_hexane", "eluent_ether"], "evaporate and chromatograph on silica gel with n-hexane/Et2O 95:5", P9),
]
r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P9)
r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P9)
r["outcomes"].append(product_amount("amount_8", molid("mol_8"), "9.8 g", 9.8, "g", P9, "product isolated after silica gel chromatography"))
reactions.append(r)

# 8 -> 9 acetylide substitution.
r = base_reaction("scheme2_8_to_9", 9)
p8 = participant("s8", molid("mol_8"), "reactant", P9); set_charge(p8, "mass", "8.0 g", 8.0, "g", P9); set_charge(p8, "amount", "22.2 mmol", 22.2, "mmol", P9)
liace = participant("r_liace", molid("reagent_li_acetylide_eda"), "reagent", P9); set_charge(liace, "mass", "3.1 g", 3.1, "g", P9); set_charge(liace, "amount", "36.9 mmol", 36.9, "mmol", P9)
dmso1 = participant("sol_dmso_suspension", molid("reagent_dmso"), "solvent", P9); set_charge(dmso1, "volume", "13.6 mL", 13.6, "mL", P9)
dmso2 = participant("sol_dmso_substrate", molid("reagent_dmso"), "solvent", P9); set_charge(dmso2, "volume", "11.2 mL", 11.2, "mL", P9)
nh4cl = participant("quench_nh4cl", molid("reagent_nh4cl"), "workup_agent", P9, "workup"); set_charge(nh4cl, "volume", "10 mL saturated NH4Cl", 10, "mL", P9)
ether = participant("extract_ether", molid("reagent_ether"), "workup_agent", P9, "workup")
brine = participant("wash_brine", molid("reagent_brine"), "workup_agent", P9, "workup")
mgso4 = participant("dry_mgso4", molid("reagent_mgso4"), "workup_agent", P9, "workup")
silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P9, "purification")
hexane = participant("eluent_hexane", molid("reagent_hexane"), "solvent", P9, "purification")
p9 = participant("p9", molid("mol_9"), "product", P9, product_role="main product")
r["participants"] = [p8, liace, dmso1, dmso2, nh4cl, ether, brine, mgso4, silica, hexane, p9]
r["reaction_scale"] = reaction_scale(molid("mol_8"), "22.2 mmol compound 8", 22.2, "mmol", P9)
r["steps"] = [
    step("add_8_to_acetylide", "reaction", ["s8", "r_liace", "sol_dmso_suspension", "sol_dmso_substrate"], "add the DMSO solution of 8 to the lithium acetylide/ethylenediamine suspension", P9, q("0 °C", 0, "°C", P9)),
    step("reaction_hold", "reaction", ["s8", "r_liace"], "stir", P9, q("23 °C", 23, "°C", P9), q("12 h", 12, "h", P9)),
    step("quench_extract", "workup", ["quench_nh4cl", "extract_ether", "wash_brine", "dry_mgso4"], "add saturated NH4Cl, extract with Et2O, wash with brine, dry and concentrate", P9),
    step("chromatography", "purification", ["purify_silica", "eluent_hexane"], "chromatograph over silica gel with n-hexane", P9),
]
r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P9)
r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P9)
r["outcomes"].append(product_amount("amount_9", molid("mol_9"), "4.5 g", 4.5, "g", P9, "product isolated after silica gel chromatography"))
reactions.append(r)

# 9 -> 10 allene rearrangement.  The amine identity is kept explicitly conflicted.
r = base_reaction("schemes2_3_9_to_10", 9)
p9 = participant("s9", molid("mol_9"), "reactant", P9); set_charge(p9, "mass", "5.0 g", 5.0, "g", P9); set_charge(p9, "amount", "16.4 mmol", 16.4, "mmol", P9)
para = participant("r_paraform", molid("reagent_paraformaldehyde"), "reagent", P9); set_charge(para, "mass", "1.0 g", 1.0, "g", P9)
cui = participant("cat_cui", molid("reagent_cui"), "catalyst", P9); set_charge(cui, "mass", "1.26 g", 1.26, "g", P9); set_charge(cui, "amount", "6.8 mmol", 6.8, "mmol", P9)
amine = participant("r_amine_conflict", amine_mol["molecule_id"], "reagent", [eid(2, "image"), eid(3), eid(3, "image"), *P9]); set_charge(amine, "volume", "75.7 mL", 75.7, "mL", P9)
dioxane = participant("sol_dioxane", molid("reagent_dioxane"), "solvent", P9); set_charge(dioxane, "volume", "75.7 mL", 75.7, "mL", P9)
hcl = participant("workup_hcl", molid("reagent_hcl"), "workup_agent", P9, "workup"); set_charge(hcl, "volume", "50 mL of 2 M HCl", 50, "mL", P9)
ether = participant("extract_ether", molid("reagent_ether"), "workup_agent", P9, "workup")
nahco3 = participant("wash_nahco3", molid("reagent_nahco3"), "workup_agent", P9, "workup")
water = participant("wash_water", molid("reagent_water"), "workup_agent", P9, "workup")
brine = participant("wash_brine", molid("reagent_brine"), "workup_agent", P9, "workup")
mgso4 = participant("dry_mgso4", molid("reagent_mgso4"), "workup_agent", P9, "workup")
silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P9, "purification")
hexane = participant("eluent_hexane", molid("reagent_hexane"), "solvent", P9, "purification")
etoac = participant("eluent_etoac", molid("reagent_etoac"), "solvent", P9, "purification")
p10 = participant("p10", molid("mol_10"), "product", P9, product_role="main product")
r["participants"] = [p9, para, cui, amine, dioxane, hcl, ether, nahco3, water, brine, mgso4, silica, hexane, etoac, p10]
r["reaction_scale"] = reaction_scale(molid("mol_9"), "16.4 mmol compound 9", 16.4, "mmol", P9)
r["steps"] = [
    step("sequential_addition", "reaction", ["s9", "r_paraform", "cat_cui", "r_amine_conflict", "sol_dioxane"], "add paraformaldehyde, CuI and the amine sequentially to the dioxane solution of 9", P9),
    step("reflux_hold", "reaction", ["s9", "r_paraform", "cat_cui", "r_amine_conflict", "sol_dioxane"], "heat under reflux (Scheme 2 additionally labels 100 °C)", [eid(2, "image"), *P9], conflict(["100 °C (Scheme 2)", "reflux (Experimental Section)"], [q("100 °C", 100, "°C", [eid(2, "image")]), text_fact("reflux", P9)], [eid(2, "image"), *P9]), q("8 h", 8, "h", [eid(2, "image"), *P9])),
    step("acid_workup", "workup", ["workup_hcl", "extract_ether"], "evaporate solvent, add 2 M HCl and extract with Et2O", P9),
    step("wash_dry", "workup", ["wash_nahco3", "wash_water", "wash_brine", "dry_mgso4"], "wash with NaHCO3, water and brine; dry over MgSO4", P9),
    step("chromatography", "purification", ["purify_silica", "eluent_hexane", "eluent_etoac"], "silica gel chromatography with n-hexane/EtOAc 30:1", P9),
]
r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P9)
r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P9)
r["outcomes"].append(product_amount("amount_10", molid("mol_10"), "3.96 g", 3.96, "g", P9, "product isolated after silica gel chromatography"))
amine_note = text_fact("The amine identity is internally inconsistent: HN(C6H11)2/dicyclohexylamine in the scheme and Experimental Section versus dihexylamine in the discussion.", [eid(2, "image"), eid(3), eid(3, "image"), *P9])
if amine_note not in r["notes"]:
    r["notes"].append(amine_note)
add_unresolved(r, "The 9→10 amine identity conflicts within the source and remains unresolved.")
reactions.append(r)


def build_cross(suffix, substrate_suffix, product_suffix, product_mass, product_yield):
    r = base_reaction(suffix, 9)
    # Preserve previously encoded scheme/prose/experimental conflicts for catalyst/base loadings.
    original = {p["participant_id"]: copy.deepcopy(p) for p in r["participants"]}
    s10 = original["s10"]; set_charge(s10, "mass", "1.0 g", 1.0, "g", P9); set_charge(s10, "amount", "4.8 mmol", 4.8, "mmol", P9)
    s15 = original["s15"]; set_charge(s15, "amount", "4.0 mmol", 4.0, "mmol", P9)
    cat = original["cat_cp2ticl2"]; set_charge(cat, "mass", "0.1 g", 0.1, "g", P9); set_charge(cat, "amount", "0.4 mmol", 0.4, "mmol", P9)
    etmg = original["r_etmgbr"]; set_charge(etmg, "amount", "40.0 mmol", 40.0, "mmol", P9)
    mg = original["r_mg"]; set_charge(mg, "mass", "0.7 g", 0.7, "g", P9); set_charge(mg, "amount", "30.0 mmol", 30.0, "mmol", P9)
    ether = original["sol_ether"]; set_charge(ether, "volume", "30 mL", 30, "mL", P9)
    argon = participant("atm_argon", molid("reagent_argon"), "additive", P9)
    nh4cl = participant("workup_nh4cl", molid("reagent_nh4cl"), "workup_agent", P9, "workup"); set_charge(nh4cl, "volume", "20 mL of 5% NH4Cl", 20, "mL", P9)
    extether = participant("extract_ether", molid("reagent_ether"), "workup_agent", P9, "workup"); set_charge(extether, "volume", "2 × 100 mL", 200, "mL", P9)
    mgso4 = participant("dry_mgso4", molid("reagent_mgso4"), "workup_agent", P9, "workup")
    silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P9, "purification")
    hexane = participant("eluent_hexane", molid("reagent_hexane"), "solvent", P9, "purification")
    etoac = participant("eluent_etoac", molid("reagent_etoac"), "solvent", P9, "purification")
    product = original["p_main"]
    r["participants"] = [s10, s15, cat, etmg, mg, ether, argon, nh4cl, extether, mgso4, silica, hexane, etoac, product]
    r["reaction_scale"] = reaction_scale(molid(substrate_suffix), f"4.0 mmol compound {substrate_suffix[-3:]}", 4.0, "mmol", P9)
    hold_time = conflict(["10 h (discussion)", "6 h (Experimental Section)"], [q("10 h", 10, "h", [eid(3), eid(3, "image")]), q("6 h", 6, "h", P9)], [eid(3), eid(3, "image"), *P9])
    r["steps"] = [
        step("charge_at_zero", "reaction", ["s10", "s15", "cat_cp2ticl2", "r_etmgbr", "r_mg", "sol_ether", "atm_argon"], "place all reaction components in a glass reactor with stirring under argon", P9, st.fact(raw="approximately 0 °C", normalized={"number": 0, "unit": "°C", "qualifier": "approx", "raw_text": "approximately 0 °C"}, evidence=P9), atmosphere="argon"),
        step("warm_to_room_temperature", "reaction", ["s10", "s15", "cat_cp2ticl2", "r_etmgbr", "r_mg", "sol_ether"], "warm to room temperature", P9, range_q("20–22 °C", 20, 22, "°C", P9)),
        step("cross_cyclomagnesiation_hold", "reaction", ["s10", "s15", "cat_cp2ticl2", "r_etmgbr", "r_mg", "sol_ether"], "stir", [eid(3), eid(3, "image"), *P9], range_q("20–22 °C", 20, 22, "°C", P9), hold_time, "argon"),
        step("nh4cl_workup", "workup", ["workup_nh4cl", "extract_ether"], "treat with 5% aqueous NH4Cl and extract twice with Et2O", P9),
        step("dry_filter_evaporate", "workup", ["dry_mgso4"], "dry combined organic phases over MgSO4, filter and remove solvents", P9),
        step("chromatography", "purification", ["purify_silica", "eluent_hexane", "eluent_etoac"], "silica gel chromatography with n-hexane/EtOAc 35:1", P9),
    ]
    r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P9)
    r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P9)
    r["outcomes"].append(product_amount(f"amount_{product_suffix}", molid(product_suffix), f"{product_mass} g", product_mass, "g", P9, "product isolated after silica gel chromatography"))
    r["procedure"]["original_text"] = text_fact("Components were charged at approximately 0 °C under argon, warmed to 20–22 °C and stirred 6 h; the mixture was treated with 5% NH4Cl, extracted and chromatographed.", P9)
    add_unresolved(r, "Cp2TiCl2, EtMgBr and Mg loadings conflict across Scheme 4, discussion ratio and Experimental Section.", "Reaction duration conflicts (10 h discussion versus 6 h Experimental Section).", "Warm-up duration was not reported.")
    return r


reactions.append(build_cross("scheme4_10_15a_to_17", "mol_15a", "mol_17", 1.8, 81))
reactions.append(build_cross("scheme4_10_15b_to_18", "mol_15b", "mol_18", 2.14, 89))

# Page-10 literature-procedure references: no unreported quantities are invented.
r = base_reaction("scheme4_18_to_19", 10)
r["procedure"]["general_procedure_reference"] = text_fact("known method 13", P10)
r["procedure"]["original_text"] = text_fact("THP deprotection of ether 18 with p-TsOH in CHCl3/MeOH using known method 13 afforded alcohol 19.", P10)
add_unresolved(r, "Page 10 supplies no reagent quantities or product yield for 18→19; prior scheme conditions remain machine-extracted and need expert review.")
reactions.append(r)

r = base_reaction("scheme4_19_to_17", 10)
r["procedure"]["general_procedure_reference"] = text_fact("known procedure 13", P10)
r["procedure"]["original_text"] = text_fact("Oxidation of alcohol 19 with Dess–Martin periodinane was carried out according to known procedure 13.", P10)
add_unresolved(r, "Page 10 supplies no experimental quantities for 19→17; the reported approximately 64% is a two-step 18→19→17 yield and is not assigned as a standalone last-step yield.")
reactions.append(r)

# 17 -> 20 staged low-temperature acetylide addition.
r = base_reaction("scheme5_17_to_20", 10)
tms = participant("r_tms_acetylene", molid("reagent_tms_acetylene"), "reagent", P10); set_charge(tms, "mass", "0.183 g", 0.183, "g", P10); set_charge(tms, "amount", "1.87 mmol", 1.87, "mmol", P10)
thf1 = participant("sol_thf_tms", molid("reagent_thf"), "solvent", P10); set_charge(thf1, "volume", "2.13 mL", 2.13, "mL", P10)
nbuli = participant("r_nbuli", molid("reagent_nbuli"), "reagent", P10); set_charge(nbuli, "volume", "0.67 mL of 2.5 M solution in n-hexane", 0.67, "mL", P10)
s17 = participant("s17", molid("mol_17"), "reactant", P10); set_charge(s17, "mass", "0.529 g", 0.529, "g", P10); set_charge(s17, "amount", "1.00 mmol", 1.00, "mmol", P10)
thf2 = participant("sol_thf_aldehyde", molid("reagent_thf"), "solvent", P10); set_charge(thf2, "volume", "2.13 mL", 2.13, "mL", P10)
li_tms = participant("r_li_tms_in_situ", molid("reagent_li_tms_acetylide"), "unknown", P10, product_role="in-situ lithium trimethylsilylacetylide")
nh4cl = participant("workup_nh4cl", molid("reagent_nh4cl"), "workup_agent", P10, "workup"); set_charge(nh4cl, "volume", "20 mL of 5% NH4Cl", 20, "mL", P10)
ether = participant("extract_ether", molid("reagent_ether"), "workup_agent", P10, "workup"); set_charge(ether, "volume", "30 mL", 30, "mL", P10)
mgso4 = participant("dry_mgso4", molid("reagent_mgso4"), "workup_agent", P10, "workup")
silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P10, "purification")
p20 = participant("p20", molid("mol_20"), "product", P10, product_role="main product")
r["participants"] = [tms, thf1, nbuli, li_tms, s17, thf2, nh4cl, ether, mgso4, silica, p20]
r["reaction_scale"] = reaction_scale(molid("mol_17"), "1.00 mmol compound 17", 1.00, "mmol", P10)
r["steps"] = [
    step("add_nbuli_minus78", "reaction", ["r_tms_acetylene", "r_nbuli", "sol_thf_tms"], "add n-BuLi solution dropwise to trimethylsilylacetylene in THF", P10, q("−78 °C", -78, "°C", P10)),
    step("adjust_minus20", "reaction", ["r_li_tms_in_situ"], "adjust temperature", P10, q("−20 °C", -20, "°C", P10)),
    step("add_aldehyde_minus20", "reaction", ["r_li_tms_in_situ", "s17", "sol_thf_aldehyde"], "add the THF solution of aldehyde 17", P10, q("−20 °C", -20, "°C", P10)),
    step("warm_to_room_temperature", "reaction", ["r_li_tms_in_situ", "s17"], "warm to room temperature", P10, range_q("20–22 °C", 20, 22, "°C", P10)),
    step("room_temperature_hold", "reaction", ["r_li_tms_in_situ", "s17"], "stir", P10, range_q("20–22 °C", 20, 22, "°C", P10), q("3 days", 3, "day", P10)),
    step("quench_extract", "workup", ["workup_nh4cl", "extract_ether"], "treat with 5% NH4Cl and extract with Et2O", P10),
    step("dry_and_isolate", "purification", ["dry_mgso4", "purify_silica"], "dry, filter, remove solvent and chromatograph on silica gel", P10),
]
r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P10)
r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P10)
r["outcomes"].append(product_amount("amount_20", molid("mol_20"), "0.55 g", 0.55, "g", P10, "product isolated after silica gel chromatography"))
r["procedure"]["addition_order"] = text_fact("n-BuLi was added to TMS-acetylene at −78 °C; after adjustment to −20 °C, the aldehyde 17 solution was added; the mixture was then warmed and held 3 days.", P10)
add_unresolved(r, "The duration of each low-temperature addition and of the warm-up was not reported.")
reactions.append(r)

# 20 -> 21 desilylation.
r = base_reaction("scheme5_20_to_21", 10)
s20 = participant("s20", molid("mol_20"), "reactant", P10); set_charge(s20, "mass", "0.55 g", 0.55, "g", P10); set_charge(s20, "amount", "0.9 mmol", 0.9, "mmol", P10)
thf = participant("sol_thf", molid("reagent_thf"), "solvent", P10); set_charge(thf, "volume", "3.74 mL", 3.74, "mL", P10)
tbaf = participant("r_tbaf", molid("reagent_tbaf"), "reagent", P10); set_charge(tbaf, "equivalents", "1.2 equiv (1 M in THF)", 1.2, "equiv", P10); tbaf["equivalence_reference"] = text_fact(molid("mol_20"), P10)
nh4cl = participant("workup_nh4cl", molid("reagent_nh4cl"), "workup_agent", P10, "workup")
ether = participant("extract_ether", molid("reagent_ether"), "workup_agent", P10, "workup"); set_charge(ether, "volume", "2 × 10 mL", 20, "mL", P10)
mgso4 = participant("dry_mgso4", molid("reagent_mgso4"), "workup_agent", P10, "workup")
silica = participant("purify_silica", molid("agent_silica_gel"), "purification_agent", P10, "purification")
p21 = participant("p21", molid("mol_21"), "product", P10, product_role="main product")
r["participants"] = [s20, thf, tbaf, nh4cl, ether, mgso4, silica, p21]
r["reaction_scale"] = reaction_scale(molid("mol_20"), "0.9 mmol compound 20", 0.9, "mmol", P10)
hold_time = conflict(["4 h (discussion)", "2 h (Experimental Section)"], [q("4 h", 4, "h", [eid(3), eid(3, "image")]), q("2 h", 2, "h", P10)], [eid(3), eid(3, "image"), *P10])
r["steps"] = [
    step("desilylation_hold", "reaction", ["s20", "r_tbaf", "sol_thf"], "add TBAF and stir", [eid(3), eid(3, "image"), *P10], "room temperature", hold_time),
    step("quench_extract", "workup", ["workup_nh4cl", "extract_ether"], "treat with saturated aqueous NH4Cl and extract twice with Et2O", P10),
    step("dry_and_isolate", "purification", ["dry_mgso4", "purify_silica"], "dry, filter, remove solvent and chromatograph on silica gel", P10),
]
r["outcomes"][0]["measurement_basis"] = text_fact("product isolated after silica gel chromatography", P10)
r["outcomes"][0]["yield_type"] = text_fact("isolated yield", P10)
r["outcomes"].append(product_amount("amount_21", molid("mol_21"), "0.48 g", 0.48, "g", P10, "product isolated after silica gel chromatography"))
add_unresolved(r, "Desilylation duration conflicts between discussion (4 h) and Experimental Section (2 h).")
reactions.append(r)

r = base_reaction("scheme5_21_to_22", 10)
r["procedure"]["general_procedure_reference"] = text_fact("published procedure 18", P10)
r["procedure"]["original_text"] = text_fact("Oxidation of alcohol 21 with Dess–Martin periodinane was carried out according to published procedure 18.", P10)
add_unresolved(r, "Page 10 gives no quantities or isolation details for 21→22; the 87% yield comes from Scheme 5/discussion evidence and its measurement basis remains unreported.")
reactions.append(r)

r = base_reaction("scheme5_22_to_1", 10)
r["procedure"]["general_procedure_reference"] = text_fact("published procedure 14", P10)
r["procedure"]["original_text"] = text_fact("Stereoselective reduction of ketone 22 with (R)-2-methyl-CBS-oxazaborolidine followed published procedure 14; lembehyne C was reported in 82% yield and 96% ee.", P10)
for o in r["outcomes"]:
    if eid(10, "image") not in o["value"]["evidence_ids"]:
        o["value"]["evidence_ids"].extend(P10)
    if o["kind"] in ("yield", "selectivity"):
        o["measurement_basis"] = missing(P10)
add_unresolved(r, "The source reports 96% ee for lembehyne C but does not report the ee measurement method or absolute R/S descriptor on pages 7–11.")
reactions.append(r)


# Collect every molecule referenced by the increment.  Existing stereochemical
# registrations are reused unchanged; new procedural agents remain unresolved.
referenced = {p["molecule_id"] for r in reactions for p in r["participants"]}
molecules = []
for mid in sorted(referenced):
    if mid in extra_molecules:
        molecules.append(extra_molecules[mid])
    else:
        molecules.append(copy.deepcopy(old_molecules[mid]))


source_context = [
    {
        "context_id": f"{DOC}:page7_figure5_mitochondrial_potential",
        "document_id": DOC,
        "category": "biological_mechanistic_assay",
        "description": "Figure 5 and surrounding prose report mitochondrial membrane-potential/apoptosis measurements in Jurkat cells; this is not a synthetic reaction.",
        "evidence_ids": [eid(7, "image"), eid(7)],
    },
    {
        "context_id": f"{DOC}:page8_figure6_ros_apoptosis",
        "document_id": DOC,
        "category": "biological_mechanistic_assay",
        "description": "Figure 6 reports ROS and apoptotic-cell-death flow-cytometry panels for lembehyne C and controls; no chemical product is defined.",
        "evidence_ids": [eid(8, "image"), eid(8)],
    },
    {
        "context_id": f"{DOC}:page8_figure7_signaling",
        "document_id": DOC,
        "category": "biological_mechanistic_assay",
        "description": "Figure 7 reports total/phosphorylated signaling-protein histograms; it is retained as biological source context, not synthesis.",
        "evidence_ids": [eid(8, "image"), eid(8)],
    },
    {
        "context_id": f"{DOC}:page9_general_analytical_methods",
        "document_id": DOC,
        "category": "characterization_methods",
        "description": "General IR, NMR, MS, elemental analysis, TLC and chromatography methods precede the experimental syntheses.",
        "evidence_ids": [eid(9, "image"), eid(9)],
    },
    {
        "context_id": f"{DOC}:pages9_10_product_characterization",
        "document_id": DOC,
        "category": "compound_characterization",
        "description": "NMR/IR/HRMS/elemental-analysis blocks for compounds 5, 6, 8, 9, 10, 17–22 and lembehyne C are characterization data, not independent reactions or yields.",
        "evidence_ids": [eid(9, "image"), eid(9), eid(10, "image"), eid(10)],
    },
    {
        "context_id": f"{DOC}:pages10_11_biology_methods",
        "document_id": DOC,
        "category": "biological_methods",
        "description": "Cell culture, cytotoxicity, apoptosis, cell-cycle, H2A.X, mitochondrial-damage and multiplex assay methods are retained separately from chemical synthesis records.",
        "evidence_ids": [eid(10, "image"), eid(10), eid(11, "image"), eid(11)],
    },
    {
        "context_id": f"{DOC}:page11_supporting_information_statement",
        "document_id": DOC,
        "category": "source_availability",
        "description": "The main article states that Supporting Information exists, but no SI file was explicitly provided or associated in this workspace; no SI linkage is inferred.",
        "evidence_ids": [eid(11, "image"), eid(11)],
    },
]

bundle = copy.deepcopy(init)
bundle["molecules"] = molecules
bundle["reactions"] = reactions
bundle["source_context"] = source_context

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "upper two-column mitochondrial-apoptosis discussion", "finding": "Biological mechanism discussion and quantitative observations; no synthetic arrow or chemical product isolation.", "evidence_ids": [eid(7, "image"), eid(7)]},
                {"region": "Figure 5 panels A–E and caption", "finding": "Jurkat mitochondrial membrane-potential/apoptosis flow-cytometry assay; retained in source_context.", "evidence_ids": [eid(7, "image"), eid(7)]},
                {"region": "lower ROS/signaling-pathway discussion", "finding": "Biological mechanistic interpretation continuing to page 8; not a synthesis experiment.", "evidence_ids": [eid(7, "image"), eid(7)]},
            ],
            "record_ids": [],
            "unresolved": ["Figure-panel numeric values were not converted into chemical reaction outcomes."],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "Figure 6 panels A–E and caption", "finding": "ROS/apoptotic cell-death assay with controls and incubation times; biological source context only.", "evidence_ids": [eid(8, "image"), eid(8)]},
                {"region": "Figure 7 histograms and caption", "finding": "Signaling-protein assay; biological mechanistic experiment, not chemical synthesis.", "evidence_ids": [eid(8, "image"), eid(8)]},
                {"region": "conclusion text", "finding": "Whole-route and anticancer conclusions; no new reaction arrow or independently specified run.", "evidence_ids": [eid(8, "image"), eid(8)]},
            ],
            "record_ids": [],
            "unresolved": ["Plotted biological values are not represented as chemistry outcomes in the current schema increment."],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "General Experimental Procedures", "finding": "Analytical and chromatography methods; retained as source_context rather than a reaction.", "evidence_ids": [eid(9, "image"), eid(9)]},
                {"region": "Synthesis of 5 and synthesis of (Z)-6", "finding": "Two detailed synthetic procedures; quantities, ordered temperature/time stages, workups and product amounts merged into existing 4→5 and 5→6 records.", "evidence_ids": [eid(9, "image"), eid(9)]},
                {"region": "Synthesis of mesylate 7 and bromide 8", "finding": "A two-stage printed procedure; 6→7 and 7→8 remain independent records and conflicting scheme/text conditions are retained.", "evidence_ids": [eid(9, "image"), eid(9)]},
                {"region": "Synthesis of 9 and allene 10", "finding": "Two detailed procedures; internal dicyclohexylamine/dihexylamine identity conflict retained without choosing.", "evidence_ids": [eid(9, "image"), eid(9)]},
                {"region": "Cross-cyclomagnesiation general procedure", "finding": "The 15a→17 and 15b→18 branches share a procedure but remain two reaction records; loading/time conflicts retained.", "evidence_ids": [eid(9, "image"), eid(9)]},
            ],
            "record_ids": [
                f"{DOC}:scheme2_4_to_5", f"{DOC}:scheme2_5_to_6", f"{DOC}:scheme2_6_to_7", f"{DOC}:scheme2_7_to_8",
                f"{DOC}:scheme2_8_to_9", f"{DOC}:schemes2_3_9_to_10", f"{DOC}:scheme4_10_15a_to_17", f"{DOC}:scheme4_10_15b_to_18",
            ],
            "unresolved": ["Several source conflicts are explicit in records; no conflicting condition or yield was selected.", "The Experimental Section starts 5 from compound 3 via in-situ 4, while the scheme gives a separate 4→5 arrow."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "characterization of 17, 18 and 19", "finding": "Spectral/elemental data are retained as characterization source_context, not synthetic experiments.", "evidence_ids": [eid(10, "image"), eid(10)]},
                {"region": "18→19 THP deprotection and 19→17 oxidation statements", "finding": "Known-method citations only; no missing quantities were invented and the approximate 64% remains a two-step result.", "evidence_ids": [eid(10, "image"), eid(10)]},
                {"region": "Preparation and characterization of 20", "finding": "Detailed staged −78 °C addition, −20 °C addition, warm-up and 3-day hold merged into existing 17→20 record.", "evidence_ids": [eid(10, "image"), eid(10)]},
                {"region": "Preparation of rac-21", "finding": "Detailed TBAF procedure merged into existing 20→21 record; 2 h/4 h conflict retained.", "evidence_ids": [eid(10, "image"), eid(10)]},
                {"region": "21→22 and 22→1 statements plus characterization", "finding": "Published-procedure references linked to existing reactions; 82% yield and 96% ee target lembehyne C, with measurement method left missing.", "evidence_ids": [eid(10, "image"), eid(10)]},
                {"region": "Biology methods beginning at lower right", "finding": "Cell-culture and assay methods; source_context only.", "evidence_ids": [eid(10, "image"), eid(10)]},
            ],
            "record_ids": [f"{DOC}:scheme4_18_to_19", f"{DOC}:scheme4_19_to_17", f"{DOC}:scheme5_17_to_20", f"{DOC}:scheme5_20_to_21", f"{DOC}:scheme5_21_to_22", f"{DOC}:scheme5_22_to_1"],
            "unresolved": ["References 13, 14 and 18 are cited for procedures, but their full experimental details are not present on this page and were not backfilled.", "Absolute R/S assignment and ee measurement method for lembehyne C are not given on pages 7–11."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "H2A.X, mitochondrial-damage and multiplex assay methods", "finding": "Biological experimental methods; no synthetic reaction arrow.", "evidence_ids": [eid(11, "image"), eid(11)]},
                {"region": "Associated Content / Supporting Information", "finding": "The article states SI availability, but no SI document is explicitly supplied in the workspace; association remains unestablished.", "evidence_ids": [eid(11, "image"), eid(11)]},
                {"region": "author information, acknowledgments and references", "finding": "Bibliographic and funding material; cited literature is not treated as this paper's experiments.", "evidence_ids": [eid(11, "image"), eid(11)]},
            ],
            "record_ids": [],
            "unresolved": ["No Supporting Information file was linked or extracted in this batch."],
        },
    ],
}

(RAW / "pages_0007_0011_raw_extraction.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(ROOT / "runs/current/submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(ROOT / "runs/current/coverage/pages_0007_0011.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
dataset_path = ROOT / "runs/current/dataset.json"
if dataset_path.exists():
    sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
    from reaction_tools import dataset_audit

    audit_result = dataset_audit(json.loads(dataset_path.read_text()))
    (RAW / "pages_0007_0011_audit.json").write_text(json.dumps(audit_result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"reactions": len(reactions), "molecules": len(molecules), "source_context": len(source_context), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
