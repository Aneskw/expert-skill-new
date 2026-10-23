#!/usr/bin/env python3
"""Build the bounded page 7--9 v3.1 extraction/update bundle."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st  # noqa: E402


DOC = "doc_60acc8016b4c0614"
BASE = ROOT / "runs/current"
OLD = st.read(BASE / "standard_records.json")
TEMPLATE = st.read(BASE / "raw_model/pages_0007_0009_template.json")
OUT = BASE / "submitted.json"
RAW = BASE / "raw_model/pages_0007_0009_extraction.json"
COVERAGE = BASE / "coverage/pages_0007_0009.json"


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def fact(raw=None, normalized=None, pages=(), method="original_extraction", *,
         status=None, verification=None, derivation=None, alternatives=()):
    return st.fact(
        raw,
        normalized,
        [eid for page in pages for eid in ev(page)],
        method=method,
        status=status,
        verification=verification,
        derivation=derivation,
        alternatives=alternatives,
    )


def missing(pages=(), *, status=None):
    return fact(pages=pages, status=status)


bundle = copy.deepcopy(TEMPLATE)
bundle["molecules"] = copy.deepcopy(OLD["molecules"])
old_molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
old_reactions = {r["reaction_id"]: r for r in OLD["reactions"]}


def add_named_molecule(mid: str, label: str, name: str, pages: tuple[int, ...], smiles=None):
    evidence = [eid for page in pages for eid in ev(page)]
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, pages)
    m["names"]["original"] = fact(name, name, pages)
    m["names"]["english"] = fact(name, name, pages)
    m["names"]["abbreviations"] = missing(pages)
    m["names"]["aliases"] = fact(label, [label], pages)
    m["cas_numbers"] = missing(pages)
    m["external_properties"] = []
    m["identity_note"] = fact(
        "The page identifies this participant by name; no CAS number is reported.",
        "Source-named participant; CAS unreported.",
        pages,
    )
    if smiles:
        diag = st.stereo_inspect(smiles)
        iso = diag["canonical_isomeric_smiles"]
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            {"source_name": name},
            iso,
            pages,
            method="structure_derivation",
            derivation="Candidate connectivity was derived from the explicit source chemical name and canonicalized with the bundled RDKit helper; parsing does not certify source identity.",
        )
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = fact(
            {"candidate_smiles": iso}, diag["atom_centers"], pages,
            method="structure_derivation",
            derivation="Bundled standard stereo diagnostics on the evidence-linked candidate graph.",
        )
        m["stereochemistry"]["double_bonds"] = fact(
            {"candidate_smiles": iso}, diag["double_bonds"], pages,
            method="structure_derivation",
            derivation="Bundled standard stereo diagnostics; carbonyl and aromatic bonds are not E/Z assignments.",
        )
        m["stereochemistry"]["axial"] = fact(
            {"candidate_smiles": iso}, [], pages,
            method="structure_derivation",
            derivation="No allene or biaryl stereogenic axis is present in the candidate graph.",
        )
        m["stereochemistry"]["mixture"] = missing(pages)
        potential = diag.get("potential_stereo", [])
        has_unspecified_atom_stereo = any(
            item.get("type", "").startswith("Atom_") and item.get("specified") == "Unspecified"
            for item in potential
        )
        assignment = (
            "toolkit flags potential tetrahedral nitrogen centers; source gives no stereochemical assignment"
            if has_unspecified_atom_stereo
            else "achiral candidate graph"
        )
        m["stereochemistry"]["assignment"] = fact(
            {"candidate_smiles": iso, "potential_stereo": potential}, assignment, pages,
            method="structure_derivation",
            derivation=(
                "The bundled diagnostic flags potential tetrahedral nitrogen centers; no source assignment is available, so the graph is not asserted to encode complete stereochemistry."
                if has_unspecified_atom_stereo
                else "No potential tetrahedral, axial, or E/Z stereogenic element is present in the candidate graph; for 2-methyl-2-butene, one alkene terminus has identical methyl substituents."
            ),
        )
        m["stereochemistry"]["chirality_status"] = fact(
            {"candidate_smiles": iso, "potential_stereo": potential},
            "unspecified" if has_unspecified_atom_stereo else "achiral",
            pages,
            method="structure_derivation",
            derivation=(
                "Potential tetrahedral nitrogen stereochemistry is left unspecified because the source supplies no assignment."
                if has_unspecified_atom_stereo
                else "The candidate graph has no potential atom or axial stereogenic element."
            ),
        )
        m["stereochemistry"]["geometry_status"] = fact(
            {"candidate_smiles": iso}, "not_applicable", pages,
            method="structure_derivation",
            derivation="The candidate graph has no double bond for which E/Z is applicable.",
        )
        m["stereochemistry"]["representation_status"] = fact(
            {"candidate_smiles": iso, "potential_stereo": potential},
            "partial" if has_unspecified_atom_stereo else "encoded",
            pages,
            method="structure_derivation",
            derivation=(
                "Connectivity is encoded, but toolkit-flagged potential tetrahedral nitrogen stereochemistry is unspecified."
                if has_unspecified_atom_stereo
                else "All stereochemical content of this achiral candidate graph is represented in the isomeric SMILES."
            ),
        )
        m["stereochemistry"]["isomeric_smiles"] = fact(
            {"candidate_smiles": iso}, iso, pages,
            method="structure_derivation",
            derivation="RDKit canonical isomeric SMILES for the candidate graph.",
        )
        m["stereochemistry"]["identity_key"] = fact(
            None, diag["stereo_graph_key"], pages,
            method="structure_derivation",
            derivation="SHA-256 stereo graph key calculated by the bundled standard stereo helper.",
        )
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(pages)
        m["generic_template"] = None
        m["stereochemistry"]["atom_centers"] = missing(pages)
        m["stereochemistry"]["double_bonds"] = missing(pages)
        m["stereochemistry"]["axial"] = missing(pages)
        m["stereochemistry"]["mixture"] = missing(pages)
        m["stereochemistry"]["assignment"] = missing(pages)
        m["stereochemistry"]["chirality_status"] = fact(
            {"source_name": name}, "unresolved", pages, method="inference",
            derivation="The bounded page names this reagent/workup material but does not provide a reviewed structure from which stereogenic elements can be audited.",
        )
        m["stereochemistry"]["geometry_status"] = fact(
            {"source_name": name}, "unresolved", pages, method="inference",
            derivation="No source drawing was used to establish E/Z applicability for this named reagent/workup material.",
        )
        m["stereochemistry"]["representation_status"] = fact(
            {"source_name": name}, "not_encoded", pages, method="inference",
            derivation="The participant is retained by source name without a concrete reviewed SMILES.",
        )
        m["stereochemistry"]["isomeric_smiles"] = missing(pages)
        m["stereochemistry"]["identity_key"] = missing(pages)
    bundle["molecules"].append(m)
    old_molecules[mid] = m


new_molecules = [
    ("m24_start", "2-chloro-1,4-dimethoxybenzene", "2-chloro-1,4-dimethoxybenzene", (7,), "COc1c(Cl)cc(OC)cc1"),
    ("m24_acid", "4-chloro-2,5-dimethoxybenzoic acid", "4-chloro-2,5-dimethoxybenzoic acid", (7,), "O=C(O)c1c(OC)cc(Cl)c(OC)c1"),
    ("hmta", "hexamethylenetetramine", "hexamethylenetetramine", (7,), "C1N2CN3CN1CN(C2)C3"),
    ("tbuoh", "t-BuOH", "tert-butanol", (7,), "CC(C)(C)O"),
    ("2m2b", "2-methyl-2-butene", "2-methyl-2-butene", (7,), "CC=C(C)C"),
    ("nah2po4_2h2o", "NaH2PO4·2H2O", "sodium dihydrogen phosphate dihydrate", (7,), None),
    ("naclo2", "NaClO2", "sodium chlorite", (7,), None),
    ("naoh", "NaOH aqueous solution", "aqueous sodium hydroxide", (7,), None),
    ("nahco3", "NaHCO3", "sodium bicarbonate", (7,), None),
    ("hcl_aq", "1 M aqueous HCl", "1 M aqueous hydrochloric acid", (7, 8), None),
]
for args in new_molecules:
    if args[0] not in old_molecules:
        add_named_molecule(*args)


def reaction_copy(suffix: str):
    return copy.deepcopy(old_reactions[f"{DOC}:{suffix}"])


def append_evidence(r, *pages):
    for page in pages:
        for eid in ev(page):
            if eid not in r["evidence_ids"]:
                r["evidence_ids"].append(eid)


def participant(r, pid, mid, role, pages, *, stage="reaction", mass=None,
                amount=None, volume=None, equiv=None, molpct=None,
                reference=None, product_role=None):
    p = st.new_participant(pid, mid, role, [eid for page in pages for eid in ev(page)], stage=stage)
    if product_role is not None:
        p["product_role"] = fact(product_role, product_role, pages)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume),
                       ("equivalents", equiv), ("mol_percent", molpct)):
        if value is not None:
            p["charges"][key] = fact(value[0], value[1], pages)
    if reference is not None:
        p["equivalence_reference"] = fact(reference, reference, pages)
    r["participants"].append(p)
    return p


def pmap(r):
    return {p["participant_id"]: p for p in r["participants"]}


def set_charge(p, key, raw, normalized, pages):
    p["charges"][key] = fact(raw, normalized, pages)


def add_step(r, sid, pages, operation, *, stage="reaction", participants=(),
             temperature=None, time=None, pressure=None, atmosphere=None):
    s = st.new_step(sid, [eid for page in pages for eid in ev(page)], stage=stage)
    s["participant_ids"] = list(participants)
    s["operation"] = fact(operation, operation, pages)
    if temperature is not None:
        s["temperature"] = fact(temperature[0], temperature[1], pages)
    if time is not None:
        s["time"] = fact(time[0], time[1], pages)
    if pressure is not None:
        s["pressure"] = fact(pressure[0], pressure[1], pages)
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, atmosphere, pages)
    r["steps"].append(s)
    return s


def add_outcome(r, oid, kind, targets, pages, raw, normalized, *,
                role=None, basis=None, yield_type=None, scope=None):
    o = st.new_outcome(oid, kind, targets, [eid for page in pages for eid in ev(page)])
    o["value"] = fact(raw, normalized, pages)
    if role is not None:
        o["product_role"] = fact(role, role, pages)
    if basis is not None:
        o["measurement_basis"] = fact(basis, basis, pages)
    if yield_type is not None:
        o["yield_type"] = fact(yield_type, yield_type, pages)
    if scope is not None:
        o["scope"] = fact(scope, scope, pages)
    r["outcomes"].append(o)
    return o


def set_partial_classification(r, pages, changes, limitations):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing(pages)
    c["bond_changes"] = []
    labels = {}
    for i, ch in enumerate(changes, 1):
        item = {"change_id": f"b{i}", **ch}
        c["bond_changes"].append(fact(
            {"source": "named/drawn transformation"}, item, pages,
            method="structure_derivation",
            derivation="Transformation-level bond change inferred from the source-named substrates/products; atom identities are not mapped.",
        ))
        if item["change_type"] in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(item["atoms"][0]["element"], item["atoms"][1]["element"], item["after_order"])
            labels.setdefault(lab["code"], {**lab, "supporting_change_ids": []})["supporting_change_ids"].append(item["change_id"])
    c["labels"] = fact(
        {"source": "partial transformation-level bond changes"}, list(labels.values()), pages,
        method="structure_derivation",
        derivation="Taxonomy labels were calculated from all and only the partial formed/promoted bonds listed above.",
    )
    c["coverage"] = fact(
        None, "partial", pages, method="inference",
        derivation="Connectivity class is clear, but no source-checked atom mapping or complete coproduct mapping is asserted.",
    )
    c["limitations"] = list(limitations)
    c["model_votes"] = []


CC_FORMED = {
    "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
    "before_order": 0, "after_order": 1, "change_type": "formed",
}
CH_GAIN = {
    "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}],
    "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1,
}
CO_FORMED = {
    "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
    "before_order": 0, "after_order": 1, "change_type": "formed",
}
OH_GAIN = {
    "atoms": [{"atom_map": None, "element": "O"}, {"atom_map": None, "element": "H"}],
    "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1,
}


# Page 7: preparation of aldehyde 24, previously missing as an experiment record.
r24 = st.new_reaction(f"{DOC}:preparation_24_from_2_chloro_1_4_dimethoxybenzene", DOC, ev(7))
r24["record_type"] = "specific_experiment"
r24["experiment_categories"] = fact(["scale_up"], ["scale_up"], (7,))
r24["location"] = fact(
    {"page": 7, "section": "4-Chloro-2,5-dimethoxybenzaldehyde (24)"},
    {"page": 7, "section": "4-Chloro-2,5-dimethoxybenzaldehyde (24)"}, (7,))
r24["reaction_name"] = fact(
    "4-Chloro-2,5-dimethoxybenzaldehyde (24)",
    "preparation of 4-chloro-2,5-dimethoxybenzaldehyde (24)", (7,))
r24["transformation_type"] = fact(
    None, "aromatic formylation", (7,), method="inference",
    derivation="Normalized from conversion of the named dimethoxybenzene to the named benzaldehyde; the paper does not label the reaction type in this paragraph.")
participant(r24, "s24_start", "m24_start", "reactant", (7,),
            volume=("4.2 mL", {"number": 4.2, "unit": "mL"}),
            amount=("30 mmol", {"number": 30, "unit": "mmol"}))
participant(r24, "hmta", "hmta", "reagent", (7,),
            mass=("4.2 g", {"number": 4.2, "unit": "g"}),
            amount=("30 mmol", {"number": 30, "unit": "mmol"}))
participant(r24, "tfa", "tfa", "solvent", (7,),
            volume=("45 mL", {"number": 45, "unit": "mL"}))
participant(r24, "p24", "m24", "product", (7,), product_role="benzaldehyde product")
participant(r24, "ice_water", "h2o", "workup_agent", (7,), stage="workup")
participant(r24, "nahco3", "nahco3", "workup_agent", (7,), stage="workup",
            mass=("70 g", {"number": 70, "unit": "g"}))
participant(r24, "hexane", "nhexane", "workup_agent", (7,), stage="workup")
participant(r24, "etoac", "etoac", "workup_agent", (7,), stage="workup")
participant(r24, "mgso4", "mgso4", "workup_agent", (7,), stage="workup")
r24["reaction_scale"] = fact(
    "2-chloro-1,4-dimethoxybenzene: 30 mmol",
    {"reference_molecule_id": "m24_start", "amount": {"number": 30, "unit": "mmol"}}, (7,))
r24["solvent_conditions"]["composition"] = fact("TFA", "TFA", (7,))
r24["solvent_conditions"]["total_volume"] = fact("45 mL TFA", {"number": 45, "unit": "mL"}, (7,))
add_step(r24, "tfa_addition", (7,), "TFA was added to the substrate/HMTA mixture at room temperature",
         participants=("s24_start", "hmta", "tfa"), temperature=("room temperature", "r.t."))
add_step(r24, "reflux_hold", (7,), "mixture was heated to reflux",
         participants=("s24_start", "hmta", "tfa", "p24"),
         temperature=("reflux", "reflux"), time=("14 h", {"number": 14, "unit": "h"}))
add_step(r24, "quench_filter", (7,), "hot solution was poured onto crushed ice; solid NaHCO3 was added cautiously; the yellow solid was filtered",
         stage="workup", participants=("ice_water", "nahco3"))
add_step(r24, "wash_extract", (7,), "solid was washed with water and hexane, dissolved in EtOAc, washed with water, dried over MgSO4, filtered, and concentrated in vacuo",
         stage="workup", participants=("ice_water", "hexane", "etoac", "mgso4", "p24"))
add_outcome(r24, "isolated_yield_24", "yield", ["m24"], (7,), "87%", {"number": 87, "unit": "%"},
            role="4-chloro-2,5-dimethoxybenzaldehyde", basis="isolated after filtration/extraction workup", yield_type="isolated")
add_outcome(r24, "product_mass_24", "product_amount", ["m24"], (7,), "5.2 g", {"number": 5.2, "unit": "g"},
            role="4-chloro-2,5-dimethoxybenzaldehyde", basis="isolated material after workup")
r24["procedure"]["original_text"] = fact(
    "To a mixture of 2-chloro-1,4-dimethoxybenzene (4.2 mL, 30 mmol) and hexamethylenetetramine (4.2 g, 30 mmol) was added TFA (45 mL) at room temperature, and the mixture was heated to reflux for 14 h.",
    "Page-7 experimental procedure for aldehyde 24", (7,))
set_partial_classification(r24, (7,), [CC_FORMED, CH_GAIN], [
    "No atom-mapped reaction SMILES; formyl carbon provenance and coproducts are not mapped.",
    "C–C/C–H formation is recorded at transformation level and remains needs_review.",
])
r24["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                 "checked_fields": [], "unresolved": [
                     "The named reaction participants are source-derived, but the complete atom mapping and mechanistic coproducts are unresolved.",
                     "The source reports TFA as a 45 mL addition; its reagent-versus-solvent role is represented as solvent and remains reviewable.",
                 ]}
r24["completeness"] = "partial"
bundle["reactions"].append(r24)


# Page 7: enrich the existing two-step 24 -> 25 record with the complete experimental procedure.
r = reaction_copy("scheme2_24_to_25")
append_evidence(r, 7)
r["location"] = fact(
    {"scheme": "Scheme 2", "experimental_page": 7, "sections": ["4-Chloro-2,5-dimethoxybenzoic Acid", "Methyl 4-Chloro-2,5-dimethoxybenzoate (25)"]},
    {"scheme": "Scheme 2", "experimental_page": 7}, (4, 7))
r["participants"] = []
participant(r, "s24", "m24", "reactant", (7,),
            mass=("2.01 g", {"number": 2.01, "unit": "g"}), amount=("10.0 mmol", {"number": 10.0, "unit": "mmol"}))
participant(r, "thf", "thf", "solvent", (7,), volume=("20 mL", {"number": 20, "unit": "mL"}))
participant(r, "tbuoh", "tbuoh", "solvent", (7,), volume=("20 mL", {"number": 20, "unit": "mL"}))
participant(r, "water_rxn", "h2o", "solvent", (7,), volume=("10 mL", {"number": 10, "unit": "mL"}))
participant(r, "2m2b", "2m2b", "reagent", (7,),
            volume=("10 mL", {"number": 10, "unit": "mL"}), amount=("100 mmol", {"number": 100, "unit": "mmol"}),
            equiv=("10 equiv", {"number": 10, "unit": "equiv"}), reference="m24")
participant(r, "phosphate", "nah2po4_2h2o", "reagent", (7,),
            mass=("9.4 g", {"number": 9.4, "unit": "g"}), amount=("60 mmol", {"number": 60, "unit": "mmol"}),
            equiv=("6.0 equiv", {"number": 6.0, "unit": "equiv"}), reference="m24")
participant(r, "chlorite", "naclo2", "reagent", (7,),
            mass=("2.7 g", {"number": 2.7, "unit": "g"}), amount=("30 mmol", {"number": 30, "unit": "mmol"}),
            equiv=("3.0 equiv", {"number": 3.0, "unit": "equiv"}), reference="m24")
participant(r, "acid_intermediate", "m24_acid", "product", (7,), product_role="crude stage-1 intermediate consumed in esterification")
participant(r, "naoh", "naoh", "workup_agent", (7,), stage="workup")
participant(r, "water_workup", "h2o", "workup_agent", (7,), stage="workup")
participant(r, "etoac_workup", "etoac", "workup_agent", (7,), stage="workup")
participant(r, "et2o", "et2o", "workup_agent", (7,), stage="workup")
participant(r, "mgso4", "mgso4", "workup_agent", (7,), stage="workup")
participant(r, "meoh", "meoh", "solvent", (7,), volume=("30 mL", {"number": 30, "unit": "mL"}))
participant(r, "h2so4", "h2so4", "reagent", (7,), volume=("1.5 mL", {"number": 1.5, "unit": "mL"}))
participant(r, "nahco3", "nahco3", "workup_agent", (7,), stage="workup")
participant(r, "p25", "m25", "product", (4, 7), product_role="two-step isolated ester product")
r["reaction_scale"] = fact(
    "24: 10.0 mmol", {"reference_molecule_id": "m24", "amount": {"number": 10.0, "unit": "mmol"}}, (7,))
r["solvent_conditions"]["composition"] = fact("THF/t-BuOH/H2O (20/20/10 mL), then MeOH (30 mL)", "THF/t-BuOH/H2O, then MeOH", (7,))
r["solvent_conditions"]["total_volume"] = fact(
    "stage 1: THF 20 mL + t-BuOH 20 mL + H2O 10 mL; stage 2: MeOH 30 mL",
    {"stage_1_components": [{"solvent": "THF", "volume": {"number": 20, "unit": "mL"}}, {"solvent": "t-BuOH", "volume": {"number": 20, "unit": "mL"}}, {"solvent": "H2O", "volume": {"number": 10, "unit": "mL"}}], "stage_2": {"solvent": "MeOH", "volume": {"number": 30, "unit": "mL"}}}, (7,))
r["steps"] = []
add_step(r, "pinnick_addition", (7,), "NaH2PO4·2H2O and NaClO2 were added at 0 °C",
         participants=("s24", "thf", "tbuoh", "water_rxn", "2m2b", "phosphate", "chlorite"), temperature=("0 °C", {"number": 0, "unit": "°C"}))
add_step(r, "pinnick_hold", (7,), "mixture was stirred at 0 °C",
         participants=("s24", "thf", "tbuoh", "water_rxn", "2m2b", "phosphate", "chlorite", "acid_intermediate"),
         temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("1 h", {"number": 1, "unit": "h"}))
add_step(r, "acid_workup", (7,), "water was added and the mixture extracted with EtOAc; organics were concentrated; the residue was dissolved in aqueous NaOH, washed twice with Et2O, acidified to pH 1, extracted with EtOAc, dried over MgSO4, filtered, and concentrated",
         stage="workup", participants=("water_workup", "etoac_workup", "naoh", "et2o", "mgso4", "acid_intermediate"))
add_step(r, "h2so4_addition", (7,), "concentrated H2SO4 was added dropwise to the crude acid in MeOH",
         participants=("acid_intermediate", "meoh", "h2so4"))
add_step(r, "esterification_reflux", (7,), "solution was heated to reflux",
         participants=("acid_intermediate", "meoh", "h2so4", "p25"), temperature=("reflux", "reflux"), time=("1.5 h", {"number": 1.5, "unit": "h"}))
add_step(r, "ester_workup", (7,), "cooled to room temperature, diluted with water, concentrated, extracted with EtOAc, washed with aqueous NaHCO3, dried over MgSO4, filtered, and concentrated",
         stage="workup", participants=("water_workup", "etoac_workup", "nahco3", "mgso4", "p25"))
r["outcomes"] = []
add_outcome(r, "yield_25", "yield", ["m25"], (4, 7), "81% (2 steps)", {"number": 81, "unit": "%"},
            role="methyl 4-chloro-2,5-dimethoxybenzoate", basis="isolated after the two-step sequence", yield_type="isolated, two steps",
            scope={"start_molecule_id": "m24", "intermediate_molecule_id": "m24_acid", "end_molecule_id": "m25"})
add_outcome(r, "product_mass_25", "product_amount", ["m25"], (7,), "1.86 g", {"number": 1.86, "unit": "g"},
            role="methyl 4-chloro-2,5-dimethoxybenzoate", basis="isolated material after two-step workup",
            scope={"start_molecule_id": "m24", "end_molecule_id": "m25"})
r["procedure"]["original_text"] = fact(
    "24 (2.01 g, 10.0 mmol) was oxidized with NaH2PO4·2H2O/NaClO2 at 0 °C for 1 h; the crude acid was esterified in MeOH with concentrated H2SO4 at reflux for 1.5 h to afford 25 (1.86 g, 81%, 2 steps).",
    "Page-7 two-step experimental sequence from 24 to 25", (7,))
set_partial_classification(r, (4, 7), [CO_FORMED], [
    "No atom-mapped reaction SMILES; oxidant coproducts and the crude acid intermediate are not completely mapped.",
    "The C–O label represents the net aldehyde-to-methyl-ester sequence and is not an atom-specific assignment.",
])
r["review"]["unresolved"] = [
    "The crude 4-chloro-2,5-dimethoxybenzoic acid was carried forward without isolation; no stage-1 yield is reported.",
    "The 81% value is explicitly a two-step isolated yield from 24 to 25 and is not assigned to esterification alone.",
    "Complete atom mapping and oxidant/coproduct identities remain unresolved.",
]
r["completeness"] = "partial"
bundle["reactions"].append(r)


# Page 7: enrich 25 -> 26 with the preparative procedure.
r = reaction_copy("scheme2_preparative_25_to_26")
append_evidence(r, 7)
r["location"] = fact(
    {"scheme": "Scheme 2", "experimental_page": 7, "section": "Methyl 4-Chloro-3-hydroxy-2,5-dimethoxybenzoate (26)"},
    {"scheme": "Scheme 2", "experimental_page": 7}, (4, 7))
r["participants"] = []
participant(r, "s25", "m25", "reactant", (7,), mass=("1.04 g", {"number": 1.04, "unit": "g"}), amount=("4.5 mmol", {"number": 4.5, "unit": "mmol"}))
participant(r, "tfa", "tfa", "solvent", (7,), volume=("10 mL", {"number": 10, "unit": "mL"}))
participant(r, "acoh", "acoh", "solvent", (7,), volume=("1.0 mL", {"number": 1.0, "unit": "mL"}))
participant(r, "pida", "pida", "reagent", (7,), mass=("1.74 g", {"number": 1.74, "unit": "g"}), amount=("5.4 mmol", {"number": 5.4, "unit": "mmol"}), equiv=("1.2 equiv", {"number": 1.2, "unit": "equiv"}), reference="m25")
participant(r, "k2co3", "k2co3", "reagent", (7,))
participant(r, "meoh", "meoh", "solvent", (7,), volume=("10 mL", {"number": 10, "unit": "mL"}))
participant(r, "p26", "m26", "product", (4, 7), product_role="phenol product")
participant(r, "etoac", "etoac", "workup_agent", (7,), stage="workup")
participant(r, "nahco3", "nahco3", "workup_agent", (7,), stage="workup")
participant(r, "mgso4", "mgso4", "workup_agent", (7,), stage="workup")
participant(r, "hexane", "nhexane", "purification_agent", (7,), stage="purification")
participant(r, "silica", "silica", "purification_agent", (7,), stage="purification")
r["reaction_scale"] = fact("25: 4.5 mmol", {"reference_molecule_id": "m25", "amount": {"number": 4.5, "unit": "mmol"}}, (7,))
r["solvent_conditions"]["composition"] = fact("TFA (10 mL)/AcOH (1.0 mL), then MeOH (10 mL)", "TFA/AcOH (10/1), then MeOH", (7,))
r["solvent_conditions"]["mixture_ratio"] = fact("TFA/AcOH = 10/1", "10/1", (7,))
r["solvent_conditions"]["total_volume"] = fact(
    "TFA 10 mL + AcOH 1.0 mL; then MeOH 10 mL",
    {"oxidation": {"TFA": {"number": 10, "unit": "mL"}, "AcOH": {"number": 1.0, "unit": "mL"}}, "hydrolysis": {"MeOH": {"number": 10, "unit": "mL"}}}, (7,))
r["steps"] = []
add_step(r, "pida_addition", (7,), "PIDA was added to 25 in TFA/AcOH at 0 °C",
         participants=("s25", "tfa", "acoh", "pida"), temperature=("0 °C", {"number": 0, "unit": "°C"}))
add_step(r, "oxidation_hold", (7,), "mixture was stirred at 0 °C",
         participants=("s25", "tfa", "acoh", "pida"), temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("1 h", {"number": 1, "unit": "h"}))
add_step(r, "basification", (7,), "reaction mixture was basified to pH 9 with aqueous K2CO3 and MeOH was added",
         participants=("k2co3", "meoh"))
add_step(r, "hydrolysis_hold", (7,), "mixture was stirred at 70 °C",
         participants=("k2co3", "meoh", "p26"), temperature=("70 °C", {"number": 70, "unit": "°C"}), time=("1 h", {"number": 1, "unit": "h"}))
add_step(r, "workup", (7,), "mixture was neutralized, concentrated, extracted with EtOAc, washed with aqueous NaHCO3, dried over MgSO4, filtered, and concentrated",
         stage="workup", participants=("etoac", "nahco3", "mgso4", "p26"))
add_step(r, "purification", (7,), "flash silica-gel chromatography (n-hexane/EtOAc = 10/1 to 3/1)",
         stage="purification", participants=("hexane", "etoac", "silica", "p26"))
r["outcomes"] = []
add_outcome(r, "yield_26", "yield", ["m26"], (4, 7), "68%", {"number": 68, "unit": "%"},
            role="methyl 4-chloro-3-hydroxy-2,5-dimethoxybenzoate", basis="isolated after flash silica-gel chromatography", yield_type="isolated")
add_outcome(r, "product_mass_26", "product_amount", ["m26"], (7,), "752 mg", {"number": 752, "unit": "mg"},
            role="methyl 4-chloro-3-hydroxy-2,5-dimethoxybenzoate", basis="isolated after flash silica-gel chromatography")
r["procedure"]["original_text"] = fact(
    "25 (1.04 g, 4.5 mmol), PIDA (1.74 g, 5.4 mmol, 1.2 equiv), TFA (10 mL), and AcOH (1.0 mL) were reacted at 0 °C for 1 h; after K2CO3/MeOH treatment at 70 °C for 1 h and chromatography, 26 (752 mg, 68%) was obtained.",
    "Page-7 preparative one-pot oxidation/hydrolysis procedure", (7,))
set_partial_classification(r, (4, 7), [CO_FORMED, OH_GAIN], [
    "No atom-mapped reaction SMILES; oxygen and proton sources are not atom-mapped.",
    "C–O and O–H formation are retained at transformation level for the phenolic OH product.",
])
r["review"]["unresolved"] = [
    "The quantity of aqueous K2CO3 used to reach pH 9 and neutralization details are not reported.",
    "The oxygen/proton sources and complete atom mapping remain unresolved.",
]
r["completeness"] = "partial"
bundle["reactions"].append(r)


def enrich_table3(entry: int):
    suffix = f"table3_entry_{entry:02d}"
    r = reaction_copy(suffix)
    append_evidence(r, 7)
    p = pmap(r)
    if entry in (1, 2, 3, 4):
        set_charge(p["s4"], "mass", "11.5 mg", {"number": 11.5, "unit": "mg"}, (7,))
        set_charge(p["s4"], "amount", "0.030 mmol", {"number": 0.030, "unit": "mmol"}, (7,))
        set_charge(p["s2"], "mass", "10.4 mg", {"number": 10.4, "unit": "mg"}, (7,))
        set_charge(p["s2"], "amount", "0.040 mmol", {"number": 0.040, "unit": "mmol"}, (7,))
        set_charge(p["solv"], "volume", "1.0 mL", {"number": 1.0, "unit": "mL"}, (7,))
        r["reaction_scale"] = fact("(−)-4: 0.030 mmol", {"reference_molecule_id": "m4", "amount": {"number": 0.030, "unit": "mmol"}}, (7,))
        r["solvent_conditions"]["total_volume"] = fact("CH2Cl2 (1.0 mL)", {"number": 1.0, "unit": "mL"}, (7,))
        r["special_conditions"]["reactor"] = fact("1.0 mL PTFE vial in high-pressure apparatus", "1.0 mL PTFE vial; high-pressure apparatus", (7,))
        # Preserve table temperature and experimental prose together; these are compatible descriptions, not silently collapsed evidence.
        r["steps"][0]["temperature"] = fact(
            {"Table 3": "25 °C", "general procedure": "room temperature"},
            {"number": 25, "unit": "°C"}, (5, 7))
        r["steps"].append(st.new_step("depressurize_remove_solvent", ev(7), stage="workup"))
        r["steps"][-1]["participant_ids"] = ["solv"]
        r["steps"][-1]["operation"] = fact("After depressurizing, the solvent was removed in vacuo.", "depressurize and remove solvent in vacuo", (7,))
        participant(r, "hexane_pur", "nhexane", "purification_agent", (7,), stage="purification")
        participant(r, "etoac_pur", "etoac", "purification_agent", (7,), stage="purification")
        participant(r, "silica_pur", "silica", "purification_agent", (7,), stage="purification")
        add_step(r, "flash_chromatography", (7,), "flash silica-gel chromatography (hexane/AcOEt = 3/1 to 2/1)",
                 stage="purification", participants=("p6", "p7", "p11", "hexane_pur", "etoac_pur", "silica_pur"))
        for o in r["outcomes"]:
            o["measurement_basis"] = fact("The yield was calculated from 1H NMR.", "1H NMR", (7,))
        r["procedure"]["general_procedure_reference"] = fact(
            "General Procedure for High-Pressure Conditions (Table 3, Entries 1−4)",
            "General Procedure for High-Pressure Conditions (Table 3, Entries 1−4)", (7,))
        r["review"]["unresolved"] = [
            "Complex reactant/product structures and complete atom mapping remain unresolved.",
            "The table states 25 °C while the general procedure says room temperature; both source expressions are retained in the temperature fact.",
            "Individual isolated product masses are not reported because the general procedure gives only an NMR-calculated mixture yield.",
        ]
    elif entry == 5:
        set_charge(p["s4"], "mass", "76.4 mg", {"number": 76.4, "unit": "mg"}, (7,))
        set_charge(p["s4"], "amount", "0.20 mmol", {"number": 0.20, "unit": "mmol"}, (7,))
        set_charge(p["s2"], "mass", "66.3 mg", {"number": 66.3, "unit": "mg"}, (7,))
        set_charge(p["s2"], "amount", "0.25 mmol", {"number": 0.25, "unit": "mmol"}, (7,))
        set_charge(p["solv"], "volume", "1.0 mL", {"number": 1.0, "unit": "mL"}, (7,))
        r["reaction_scale"] = fact("(−)-4: 0.20 mmol", {"reference_molecule_id": "m4", "amount": {"number": 0.20, "unit": "mmol"}}, (7,))
        r["solvent_conditions"]["total_volume"] = fact("CH2Cl2 (1.0 mL)", {"number": 1.0, "unit": "mL"}, (7,))
        r["special_conditions"]["reactor"] = fact("1.0 mL PTFE vial in high-pressure apparatus", "1.0 mL PTFE vial; high-pressure apparatus", (7,))
        r["steps"][0]["temperature"] = fact({"Table 3": "25 °C", "procedure": "room temperature"}, {"number": 25, "unit": "°C"}, (5, 7))
        r["steps"].insert(1, st.new_step("depressurize_concentrate", ev(7), stage="workup"))
        r["steps"][1]["participant_ids"] = ["solv"]
        r["steps"][1]["operation"] = fact("After depressurizing, the solvent was removed in vacuo.", "depressurize and remove CH2Cl2 in vacuo", (7,))
        r["steps"][2]["operation"] = fact("The mixture was stirred for 9 h at 60 °C.", "neat carbonyl-ene conversion after solvent removal", (7,))
        r["steps"].append(st.new_step("consumption_check", ev(7), stage="characterization"))
        r["steps"][-1]["participant_ids"] = ["p7"]
        r["steps"][-1]["operation"] = fact("After the full consumption of 7 was checked by 1H NMR", "full consumption of 7 checked by 1H NMR", (7,))
        participant(r, "hexane_pur", "nhexane", "purification_agent", (7,), stage="purification")
        participant(r, "etoac_pur", "etoac", "purification_agent", (7,), stage="purification")
        participant(r, "silica_pur", "silica", "purification_agent", (7,), stage="purification")
        add_step(r, "flash_chromatography", (7,), "flash silica-gel chromatography (n-hexane/EtOAc = 2/1)",
                 stage="purification", participants=("p6", "p11", "hexane_pur", "etoac_pur", "silica_pur"))
        for o in r["outcomes"]:
            o["measurement_basis"] = fact("isolated after flash silica-gel chromatography", "isolated after flash silica-gel chromatography", (7,))
        add_outcome(r, "product_mass_11", "product_amount", ["m11"], (7,), "88.0 mg", {"number": 88.0, "unit": "mg"}, role="chloropupukeanolide D", basis="isolated after flash silica-gel chromatography")
        add_outcome(r, "product_mass_6", "product_amount", ["m6"], (7,), "26.8 mg", {"number": 26.8, "unit": "mg"}, role="chloropestolide B", basis="isolated after flash silica-gel chromatography")
        r["procedure"]["original_text"] = fact(
            "(−)-4 (76.4 mg, 0.20 mmol) and (+)-2 (66.3 mg, 0.25 mmol) in CH2Cl2 (1.0 mL) were held at 1.0 GPa at room temperature for 64 h, concentrated, and heated at 60 °C for 9 h; chromatography gave 11 (88.0 mg, 69%) and 6 (26.8 mg, 21%).",
            "Page-7 Table 3 entry 5 procedure", (7,))
        r["notes"].append(fact(
            {"Table 3 ratio 2/4": "1.3/1.0", "procedure charges": "0.25/0.20 mmol"},
            "The table reports a rounded 2/4 ratio of 1.3/1.0; the exact listed charges correspond to 1.25/1.0. Both are retained without treating either as an author error.",
            (5, 7), method="inference",
            derivation="Direct comparison of the table ratio with the explicit experimental charges."))
        r["review"]["unresolved"] = [
            "Complex reactant/product structures and complete atom mapping remain unresolved.",
            "Table 3 gives 2/4 = 1.3/1.0 whereas 0.25/0.20 mmol equals 1.25/1.0; the rounded table value and exact charges are retained separately.",
            "The table says 25 °C and the procedure says room temperature; both expressions are retained.",
        ]
    elif entry == 7:
        set_charge(p["s4"], "mass", "37.9 mg", {"number": 37.9, "unit": "mg"}, (7,))
        set_charge(p["s4"], "amount", "0.10 mmol", {"number": 0.10, "unit": "mmol"}, (7,))
        set_charge(p["s2"], "mass", "26.2 mg", {"number": 26.2, "unit": "mg"}, (7,))
        set_charge(p["s2"], "amount", "0.10 mmol", {"number": 0.10, "unit": "mmol"}, (7,))
        participant(r, "ch2cl2_pretreat", "ch2cl2", "solvent", (7,))
        r["reaction_scale"] = fact("(−)-4 and (+)-2: 0.10 mmol each", {"reference_molecule_id": "m4", "amount": {"number": 0.10, "unit": "mmol"}}, (7,))
        r["special_conditions"]["reactor"] = fact("sealed tube", "sealed tube", (7,))
        r["procedure"]["pretreatment"] = fact(
            "To make a homogeneous mixture, the mixture was dissolved in a small amount of CH2Cl2, and CH2Cl2 was removed in vacuo using a rotary evaporator.",
            "temporary CH2Cl2 homogenization followed by complete solvent removal", (7,))
        r["steps"].insert(0, st.new_step("homogenize_remove_solvent", ev(7), stage="reaction"))
        r["steps"][0]["participant_ids"] = ["s2", "s4", "ch2cl2_pretreat"]
        r["steps"][0]["operation"] = fact("dissolved in a small amount of CH2Cl2, then CH2Cl2 was removed in vacuo", "temporary CH2Cl2 homogenization and solvent removal", (7,))
        r["steps"][1]["temperature"] = fact({"Table 3": "25 °C", "procedure": "room temperature"}, {"number": 25, "unit": "°C"}, (5, 7))
        r["steps"][1]["operation"] = fact("The mixture was placed for 120 h at room temperature.", "neat reaction in sealed tube", (7,))
        r["steps"][2]["operation"] = fact("the mixture was stirred for 68 h at 60 °C", "neat carbonyl-ene conversion", (7,))
        r["steps"].append(st.new_step("starting_material_check", ev(7), stage="characterization"))
        r["steps"][-1]["participant_ids"] = ["s2", "s4"]
        r["steps"][-1]["operation"] = fact("full consumption of 2 and 4 was checked by 1H NMR", "full consumption of 2 and 4 checked by 1H NMR", (7,))
        r["steps"].append(st.new_step("intermediate_check", ev(7), stage="characterization"))
        r["steps"][-1]["participant_ids"] = ["p7"]
        r["steps"][-1]["operation"] = fact("After the full consumption of 7 was checked", "full consumption of 7 checked; method unreported", (7,))
        participant(r, "hexane_pur", "nhexane", "purification_agent", (7,), stage="purification")
        participant(r, "etoac_pur", "etoac", "purification_agent", (7,), stage="purification")
        participant(r, "silica_pur", "silica", "purification_agent", (7,), stage="purification")
        add_step(r, "flash_chromatography", (7,), "flash silica-gel chromatography (n-hexane/EtOAc = 2/1)",
                 stage="purification", participants=("p6", "p11", "hexane_pur", "etoac_pur", "silica_pur"))
        for o in r["outcomes"]:
            o["measurement_basis"] = fact("isolated after flash silica-gel chromatography", "isolated after flash silica-gel chromatography", (7,))
        add_outcome(r, "product_mass_11", "product_amount", ["m11"], (7,), "36.5 mg", {"number": 36.5, "unit": "mg"}, role="chloropupukeanolide D", basis="isolated after flash silica-gel chromatography")
        add_outcome(r, "product_mass_6", "product_amount", ["m6"], (7,), "16.0 mg", {"number": 16.0, "unit": "mg"}, role="chloropestolide B", basis="isolated after flash silica-gel chromatography")
        r["procedure"]["original_text"] = fact(
            "(−)-4 (37.9 mg, 0.10 mmol) and (+)-2 (26.2 mg, 0.10 mmol) were homogenized with a small amount of CH2Cl2, made neat, held 120 h at room temperature, then 68 h at 60 °C; chromatography gave 11 (36.5 mg, 57%) and 6 (16.0 mg, 25%).",
            "Page-7 Table 3 entry 7 procedure", (7,))
        r["review"]["unresolved"] = [
            "Complex reactant/product structures and complete atom mapping remain unresolved.",
            "The amount of temporary CH2Cl2 used for homogenization is unreported.",
            "The method used for the second check of full consumption of intermediate 7 is not stated.",
            "The table states 25 °C and the procedure says room temperature; both expressions are retained.",
        ]
    r["completeness"] = "partial"
    return r


for entry in (1, 2, 3, 4, 5, 7):
    bundle["reactions"].append(enrich_table3(entry))


# Page 7: complete procedure for the Scheme 3a migration.
r = reaction_copy("scheme3a_migration_11_to_1")
append_evidence(r, 7)
p = pmap(r)
set_charge(p["s11"], "mass", "24.0 mg", {"number": 24.0, "unit": "mg"}, (7,))
set_charge(p["s11"], "amount", "0.037 mmol", {"number": 0.037, "unit": "mmol"}, (7,))
set_charge(p["dmf"], "volume", "1.0 mL", {"number": 1.0, "unit": "mL"}, (7,))
set_charge(p["base"], "volume", "0.055 mL of 1.0 M solution", {"number": 0.055, "unit": "mL"}, (7,))
set_charge(p["base"], "amount", "0.055 mmol", {"number": 0.055, "unit": "mmol"}, (7,))
set_charge(p["base"], "equivalents", "1.5 equiv", {"number": 1.5, "unit": "equiv"}, (7,))
p["base"]["equivalence_reference"] = fact("m11", "m11", (7,))
participant(r, "thf_carrier", "thf", "solvent", (7,))
participant(r, "hcl", "hcl_aq", "workup_agent", (7,), stage="workup")
participant(r, "etoac", "etoac", "workup_agent", (7,), stage="workup")
participant(r, "mgso4", "mgso4", "workup_agent", (7,), stage="workup")
participant(r, "hexane_pur", "nhexane", "purification_agent", (7,), stage="purification")
participant(r, "silica_pur", "silica", "purification_agent", (7,), stage="purification")
r["reaction_scale"] = fact("11: 0.037 mmol", {"reference_molecule_id": "m11", "amount": {"number": 0.037, "unit": "mmol"}}, (7,))
r["solvent_conditions"]["composition"] = fact("DMF (1.0 mL); t-BuOK delivered as 1.0 M solution in THF", "DMF with THF reagent carrier", (7,))
r["solvent_conditions"]["total_volume"] = fact(
    "DMF 1.0 mL plus 0.055 mL of t-BuOK/THF solution",
    {"DMF": {"number": 1.0, "unit": "mL"}, "tBuOK_solution_in_THF": {"number": 0.055, "unit": "mL"}}, (7,))
r["steps"] = []
add_step(r, "base_addition", (7,), "t-BuOK solution in THF was added to 11 in DMF at 0 °C",
         participants=("s11", "dmf", "base", "thf_carrier"), temperature=("0 °C", {"number": 0, "unit": "°C"}))
add_step(r, "migration_hold", (7,), "mixture was stirred at 0 °C",
         participants=("s11", "dmf", "base", "thf_carrier", "p1"), temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("30 min", {"number": 30, "unit": "min"}))
add_step(r, "acid_quench_extract", (7,), "quenched with 1 M aqueous HCl, extracted with EtOAc, dried over MgSO4, filtered, and concentrated in vacuo",
         stage="workup", participants=("hcl", "etoac", "mgso4", "p1"))
add_step(r, "flash_chromatography", (7,), "flash silica-gel chromatography (n-hexane/EtOAc = 2/1)",
         stage="purification", participants=("hexane_pur", "etoac", "silica_pur", "p1"))
r["outcomes"] = []
add_outcome(r, "yield_1", "yield", ["m1"], (5, 6, 7), "83%", {"number": 83, "unit": "%"}, role="chloropupukeananin", basis="isolated after flash silica-gel chromatography", yield_type="isolated")
add_outcome(r, "product_mass_1", "product_amount", ["m1"], (7,), "20.0 mg", {"number": 20.0, "unit": "mg"}, role="chloropupukeananin", basis="isolated after flash silica-gel chromatography")
r["procedure"]["original_text"] = fact(
    "11 (24.0 mg, 0.037 mmol) in DMF (1.0 mL) was treated with t-BuOK in THF (1.0 M, 0.055 mL, 0.055 mmol, 1.5 equiv) at 0 °C for 30 min; workup and chromatography gave 1 (20.0 mg, 83%).",
    "Page-7 migration procedure", (7,))
r["review"]["unresolved"] = [
    "Complex rearrangement lacks atom mapping; classification remains partial.",
    "The volume of aqueous HCl and extraction/drying solvent quantities are unreported.",
]
r["completeness"] = "partial"
bundle["reactions"].append(r)


# Page 8: complete one-pot procedure; same experiment as Scheme 3b/abstract, so update rather than duplicate.
r = reaction_copy("scheme3b_one_pot_2_4_to_1_6")
append_evidence(r, 8)
p = pmap(r)
set_charge(p["s4"], "mass", "47.0 mg", {"number": 47.0, "unit": "mg"}, (8,))
set_charge(p["s4"], "amount", "0.124 mmol", {"number": 0.124, "unit": "mmol"}, (8,))
set_charge(p["s2"], "mass", "32.5 mg", {"number": 32.5, "unit": "mg"}, (8,))
set_charge(p["s2"], "amount", "0.124 mmol", {"number": 0.124, "unit": "mmol"}, (8,))
set_charge(p["dmf"], "volume", "3.0 mL", {"number": 3.0, "unit": "mL"}, (8,))
set_charge(p["base"], "volume", "0.18 mL of 1.0 M solution", {"number": 0.18, "unit": "mL"}, (8,))
set_charge(p["base"], "amount", "0.18 mmol", {"number": 0.18, "unit": "mmol"}, (8,))
set_charge(p["base"], "equivalents", "1.5 equiv", {"number": 1.5, "unit": "equiv"}, (8,))
p["base"]["equivalence_reference"] = fact("m4", "m4", (8,))
participant(r, "ch2cl2_pretreat", "ch2cl2", "solvent", (8,))
participant(r, "thf_carrier", "thf", "solvent", (8,))
participant(r, "hcl", "hcl_aq", "workup_agent", (8,), stage="workup")
participant(r, "etoac", "etoac", "workup_agent", (8,), stage="workup")
participant(r, "mgso4", "mgso4", "workup_agent", (8,), stage="workup")
participant(r, "hexane_pur", "nhexane", "purification_agent", (8,), stage="purification")
participant(r, "silica_pur", "silica", "purification_agent", (8,), stage="purification")
r["reaction_scale"] = fact("(−)-4 and (+)-2: 0.124 mmol each", {"reference_molecule_id": "m4", "amount": {"number": 0.124, "unit": "mmol"}}, (8,))
r["solvent_conditions"]["composition"] = fact("neat after temporary CH2Cl2 homogenization; then DMF (3.0 mL) with THF reagent carrier", "neat, then DMF/THF carrier", (8,))
r["special_conditions"]["reactor"] = fact("sealed tube; refrigerator for 5 °C stage", "sealed tube; refrigerator", (8,))
r["procedure"]["pretreatment"] = fact(
    "To make a homogeneous mixture, the mixture was dissolved in a small amount of CH2Cl2, and CH2Cl2 was removed in vacuo using a rotary evaporator.",
    "temporary CH2Cl2 homogenization followed by complete solvent removal", (8,))
r["steps"] = []
add_step(r, "homogenize_remove_solvent", (8,), "dissolved in a small amount of CH2Cl2 and CH2Cl2 was removed in vacuo using a rotary evaporator",
         participants=("s2", "s4", "ch2cl2_pretreat"))
add_step(r, "neat_cycloaddition", (8,), "sealed neat mixture was placed in a refrigerator",
         participants=("s2", "s4"), temperature=("5 °C", {"number": 5, "unit": "°C"}), time=("168 h", {"number": 168, "unit": "h"}), pressure=("atmospheric pressure", "0.1 MPa"))
add_step(r, "starting_material_check", (8,), "full consumption of 2 and 4 was checked by 1H NMR",
         stage="characterization", participants=("s2", "s4"))
add_step(r, "carbonyl_ene_heating", (8,), "mixture was heated",
         participants=("s2", "s4"), temperature=("60 °C", {"number": 60, "unit": "°C"}), time=("9 h", {"number": 9, "unit": "h"}), pressure=("atmospheric pressure", "0.1 MPa"))
add_step(r, "base_addition", (8,), "DMF and t-BuOK solution in THF were added at 0 °C",
         participants=("dmf", "base", "thf_carrier"), temperature=("0 °C", {"number": 0, "unit": "°C"}))
add_step(r, "migration_hold", (8,), "mixture was stirred at 0 °C",
         participants=("dmf", "base", "thf_carrier", "p1", "p6"), temperature=("0 °C", {"number": 0, "unit": "°C"}), time=("30 min", {"number": 30, "unit": "min"}))
add_step(r, "acid_quench_extract", (8,), "quenched with 1 M aqueous HCl, extracted with EtOAc, dried over MgSO4, filtered, and concentrated in vacuo",
         stage="workup", participants=("hcl", "etoac", "mgso4", "p1", "p6"))
add_step(r, "flash_chromatography", (8,), "flash silica-gel chromatography (n-hexane/EtOAc = 2/1)",
         stage="purification", participants=("hexane_pur", "etoac", "silica_pur", "p1", "p6"))
r["outcomes"] = []
add_outcome(r, "yield_1", "yield", ["m1"], (1, 5, 6, 8), "64%", {"number": 64, "unit": "%"}, role="chloropupukeananin", basis="isolated after flash silica-gel chromatography", yield_type="isolated")
add_outcome(r, "yield_6", "yield", ["m6"], (5, 6, 8), "20%", {"number": 20, "unit": "%"}, role="chloropestolide B", basis="isolated after flash silica-gel chromatography", yield_type="isolated")
add_outcome(r, "product_mass_1", "product_amount", ["m1"], (8,), "51.0 mg", {"number": 51.0, "unit": "mg"}, role="chloropupukeananin", basis="isolated after flash silica-gel chromatography")
add_outcome(r, "product_mass_6", "product_amount", ["m6"], (8,), "15.9 mg", {"number": 15.9, "unit": "mg"}, role="chloropestolide B", basis="isolated after flash silica-gel chromatography")
r["procedure"]["original_text"] = fact(
    "(−)-4 (47.0 mg, 0.124 mmol) and (+)-2 (32.5 mg, 0.124 mmol) were made neat, held at 5 °C for 168 h and 60 °C for 9 h, then treated with DMF (3.0 mL) and t-BuOK/THF (1.0 M, 0.18 mL, 0.18 mmol, 1.5 equiv) at 0 °C for 30 min; chromatography gave 1 (51.0 mg, 64%) and 6 (15.9 mg, 20%).",
    "Page-8 one-pot experimental procedure", (8,))
r["review"]["unresolved"] = [
    "The amount of temporary CH2Cl2 used for homogenization is unreported.",
    "Complex cascade/migration atom mapping remains unresolved; classification records transformation-level bond classes only.",
    "The volume of aqueous HCl and extraction/drying solvent quantities are unreported.",
]
r["completeness"] = "partial"
bundle["reactions"].append(r)


bundle["source_context"] = [
    {
        "context_id": f"{DOC}:page7_characterization_24_acid_25_26_11_6",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": ev(7),
        "text": "Page 7 reports Rf/mp/NMR/IR/HRMS and optical-rotation characterization for 24, the crude acid, 25, 26, 11, and 6. These peak/property listings are not separate synthesis experiments.",
    },
    {
        "context_id": f"{DOC}:page8_characterization_1",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": ev(8),
        "text": "Page 8 reports Rf/mp/NMR/IR/HRMS/optical rotation for (+)-chloropupukeananin (1); this is characterization of the one-pot product, not another synthesis record.",
    },
    {
        "context_id": f"{DOC}:page8_si_availability_statement",
        "document_id": DOC,
        "kind": "supporting_information_notice",
        "evidence_ids": ev(8),
        "text": "The article states that Supporting Information is available and lists its contents, but no SI PDF was supplied in this workspace; no SI document association was made.",
    },
    {
        "context_id": f"{DOC}:page8_author_acknowledgment_references",
        "document_id": DOC,
        "kind": "article_metadata_and_references",
        "evidence_ids": ev(8),
        "text": "Author information, acknowledgments, and the opening reference list contain no additional current-paper experiment.",
    },
    {
        "context_id": f"{DOC}:page9_references",
        "document_id": DOC,
        "kind": "references",
        "evidence_ids": ev(9),
        "text": "Page 9 is the continuation of the cited reference list; cited paper titles are literature context, not experiments of this article.",
    },
    {
        "context_id": f"{DOC}:page9_optical_rotation_note_24",
        "document_id": DOC,
        "kind": "author_interpretation_of_characterization",
        "evidence_ids": ev(9),
        "text": "Reference note 24 discusses the sign difference between synthetic 11 ([alpha]D24 -3.8) and a reported natural sample (+13.2) that contained (+)-1. This author interpretation is retained as context and is not used to invent a new reaction or overwrite unresolved atomwise stereochemistry.",
    },
]


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "left column, continuation and 4-chloro-2,5-dimethoxybenzaldehyde (24)", "finding": "independent preparative formylation of 2-chloro-1,4-dimethoxybenzene to 24 extracted; characterization kept as source context", "evidence_ids": ev(7)},
                {"region": "left column, 4-chloro-2,5-dimethoxybenzoic acid and methyl ester 25", "finding": "two-step 24-to-25 experiment merged into the existing Scheme 2 record; crude acid stage has no separate yield", "evidence_ids": ev(7)},
                {"region": "left column, methyl 4-chloro-3-hydroxy-2,5-dimethoxybenzoate (26)", "finding": "preparative 25-to-26 oxidation/hydrolysis merged into the existing Scheme 2 record with full staged conditions", "evidence_ids": ev(7)},
                {"region": "column break/right column, Table 3 procedures", "finding": "general procedure for entries 1-4 and specific procedures for entries 5 and 7 merged into their existing table-row records; entry 6 has no detailed procedure here", "evidence_ids": ev(7)},
                {"region": "right column, characterization of 11 and 6", "finding": "spectral/property data classified as characterization, not synthetic experiments", "evidence_ids": ev(7)},
                {"region": "right column bottom, (+)-chloropupukeananin (1)", "finding": "complete 11-to-1 migration procedure merged into existing Scheme 3a record", "evidence_ids": ev(7)},
            ],
            "record_ids": [
                f"{DOC}:preparation_24_from_2_chloro_1_4_dimethoxybenzene",
                f"{DOC}:scheme2_24_to_25",
                f"{DOC}:scheme2_preparative_25_to_26",
                *[f"{DOC}:table3_entry_{i:02d}" for i in (1, 2, 3, 4, 5, 7)],
                f"{DOC}:scheme3a_migration_11_to_1",
            ],
            "unresolved": [
                "Complex natural-product structures remain without source-checked atom-mapped/isomeric SMILES.",
                "Table 3 uses 25 °C while procedures say room temperature; both expressions are retained.",
                "Entry 5 table ratio 2/4 = 1.3/1.0 differs from the exact 0.25/0.20 mmol charge ratio (1.25/1.0); both are retained.",
                "No SI file was supplied or associated.",
            ],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "left column top, One-Pot Procedure from (+)-2 and (−)-4", "finding": "complete one-pot experiment merged with existing Scheme 3b record; no duplicate created", "evidence_ids": ev(8)},
                {"region": "left column, (+)-chloropupukeananin (1) characterization", "finding": "spectral/property data retained as source context, not a new synthesis", "evidence_ids": ev(8)},
                {"region": "left column lower half, Associated Content and author information", "finding": "SI availability statement and metadata only; no SI PDF associated", "evidence_ids": ev(8)},
                {"region": "right column, acknowledgments and references 1-7", "finding": "no current-paper experiment; citations are literature context", "evidence_ids": ev(8)},
            ],
            "record_ids": [f"{DOC}:scheme3b_one_pot_2_4_to_1_6"],
            "unresolved": [
                "Complex structures and atom mapping for the cascade remain unresolved.",
                "Temporary CH2Cl2 volume and several workup solvent quantities are unreported.",
                "The article advertises SI, but no SI document was provided in the workspace.",
            ],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "both columns, references 7-24", "finding": "reference list only; cited experiments were not reclassified as experiments of this paper", "evidence_ids": ev(9)},
                {"region": "right column, reference note 24", "finding": "author discussion of optical rotation for synthetic 11 versus an inseparable natural 11/(+)-1 mixture; retained as characterization context", "evidence_ids": ev(9)},
            ],
            "record_ids": [],
            "unresolved": [
                "The optical-rotation discussion does not supply atomwise stereochemical assignments and was not used to invent them.",
                "No synthetic experiment occurs on this page.",
            ],
        },
    ],
}


raw_inventory = {
    "document_id": DOC,
    "bounded_pages": [7, 8, 9],
    "visual_inspection": "page7.jpg, page8.jpg, and page9.jpg were each inspected at original detail",
    "new_reaction_ids": [f"{DOC}:preparation_24_from_2_chloro_1_4_dimethoxybenzene"],
    "updated_reaction_ids": [
        f"{DOC}:scheme2_24_to_25",
        f"{DOC}:scheme2_preparative_25_to_26",
        *[f"{DOC}:table3_entry_{i:02d}" for i in (1, 2, 3, 4, 5, 7)],
        f"{DOC}:scheme3a_migration_11_to_1",
        f"{DOC}:scheme3b_one_pot_2_4_to_1_6",
    ],
    "new_molecule_ids": [row[0] for row in new_molecules if row[0] in old_molecules],
    "source_context_ids": [c["context_id"] for c in bundle["source_context"]],
    "decisions": [
        "Cross-page procedures were merged into existing reaction IDs; distinct table rows remain distinct records.",
        "Independent characterization and mass-spectral/NMR/property data were retained as source_context.",
        "The SI availability notice was recorded, but no SI document was associated because no SI PDF was supplied.",
        "All machine records remain needs_review; schema validation is not chemistry verification.",
    ],
}


st.write(OUT, bundle)
st.write(RAW, raw_inventory)
st.write(COVERAGE, coverage)
print(json.dumps({
    "bundle_molecules": len(bundle["molecules"]),
    "new_or_updated_reactions": len(bundle["reactions"]),
    "source_context": len(bundle["source_context"]),
    "coverage_pages": [p["page"] for p in coverage["pages"]],
}, ensure_ascii=False))
