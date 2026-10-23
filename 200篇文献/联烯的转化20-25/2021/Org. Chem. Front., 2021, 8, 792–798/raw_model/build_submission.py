#!/usr/bin/env python3
"""Build the source-derived v3.1 bundle for doc_45ed223b4e2508bb.

This is a session-local serialization aid.  It does not alter extraction tools,
schemas, taxonomies, or validation rules.  All chemistry remains needs_review.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_45ed223b4e2508bb"
RUN = ROOT / "runs/current"
P = {n: f"{DOC}:page{n}:image" for n in range(1, 7)}
PT = {n: f"{DOC}:page{n}:text" for n in range(1, 7)}


def ext(raw, normalized=None, evidence=(), method="original_extraction"):
    if normalized is None:
        normalized = raw
    return st.fact(raw, normalized, evidence, method=method)


def derived(normalized, evidence, reason):
    return st.fact(None, normalized, evidence, method="structure_derivation", derivation=reason)


def inferred(normalized, evidence, reason):
    return st.fact(None, normalized, evidence, method="inference", derivation=reason)


def missing():
    return st.fact()


bundle = st.empty(json.loads((RUN / "dataset.json").read_text())["documents"])
molecules = {}


def molecule(mid, label, name=None, evidence=(P[3],), smiles=None,
             chirality="unresolved", geometry="unresolved", rep_status="unresolved",
             stereo_note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = ext(label, label, evidence)
    m["names"]["original"] = ext(name or label, name or label, evidence)
    m["identity_note"] = ext(
        {"source_label": label, "note": stereo_note or "Identity is retained at the resolution printed in the source."},
        {"source_label": label, "note": stereo_note or "Identity is retained at the resolution printed in the source."},
        evidence,
    )
    if smiles:
        m["representation_type"] = "smiles"
        m["smiles"] = derived(
            smiles, evidence,
            "Candidate connectivity was transcribed from the displayed structure/name; RDKit parsing does not verify source identity.",
        )
    else:
        m["representation_type"] = "unresolved"
    for key, value in (
        ("chirality_status", chirality),
        ("geometry_status", geometry),
        ("representation_status", rep_status),
    ):
        m["stereochemistry"][key] = inferred(
            value, evidence,
            "Status assigned from the visible source drawing and labels; absence of a label is not treated as racemic or as proof of achirality.",
        )
    if stereo_note:
        m["stereochemistry"]["assignment"] = ext(
            {"source_region": label, "observation": stereo_note},
            {"source_region": label, "observation": stereo_note}, evidence,
        )
    molecules[mid] = m
    return mid


def generic_molecule(mid, label, fixed, variables, evidence):
    molecule(mid, label, label, evidence=evidence, chirality="unresolved",
             geometry="unresolved", rep_status="not_encoded")
    m = molecules[mid]
    m["representation_type"] = "generic_template"
    m["generic_template"] = {
        "fixed_structure": ext({"drawing": fixed}, {"drawing": fixed}, evidence),
        "variables": [
            {
                "label": var,
                "attachment_points": ext({"source": "generic drawing", "variable": var}, {"source": "generic drawing", "variable": var}, evidence),
                "definition": ext(definition, definition, evidence),
                "allowed_values": missing(),
                "constraints": missing(),
            }
            for var, definition in variables
        ],
        "instance_mapping": [],
        "reason": ext("Generic R-group structure in Scheme 1/8", "Generic R-group structure in Scheme 1/8", evidence),
    }
    return mid


# Core structures for which a candidate graph can be transcribed safely enough to aid review.
molecule("m_1a", "1a", "N-(2-pyridyl)sulfonyl-2-iodoaniline", (P[1], P[2], P[3], P[6]),
         "O=S(=O)(Nc1ccccc1I)c1ccccn1", "achiral", "not_applicable", "not_applicable",
         "No stereogenic element is present in the candidate graph.")
molecule("m_2a", "2a", "phenylallene", (P[3], P[4], P[5], P[6]),
         "C=C=Cc1ccccc1", "achiral", "not_applicable", "not_applicable",
         "The displayed allene has a terminal CH2 and therefore does not define an allene stereogenic axis.")
molecule("m_2k", "2k", "3-methoxyphenylallene", (P[2], P[3]),
         "C=C=Cc1cccc(OC)c1", "achiral", "not_applicable", "not_applicable",
         "The displayed allene has a terminal CH2 and therefore does not define an allene stereogenic axis.")
molecule("m_3aa", "3aa", "cyclocarbonylation product 3aa", (P[2], P[3], P[4], P[5]),
         "O=C1C(=Cc2ccccc2)CN(S(=O)(=O)c2ccccn2)c2ccccc21",
         "achiral", "unspecified", "partial",
         "The source confirms the connectivity by X-ray (CCDC 2023302), but no E/Z descriptor is printed for the exocyclic alkene; geometry is not invented.")

# Common catalysts, ligands, solvents, reagents, and purification agents.
common = {
    "m_tfben": ("TFBen", "benzene-1,3,5-triyl triformate"),
    "m_et3n": ("Et3N", "triethylamine"),
    "m_pd_oac2": ("Pd(OAc)2", "palladium(II) acetate"),
    "m_pd_cl2": ("PdCl2", "palladium(II) chloride"),
    "m_pd_tfa2": ("Pd(TFA)2", "palladium(II) trifluoroacetate"),
    "m_pd_acac2": ("Pd(acac)2", "palladium(II) acetylacetonate"),
    "m_pd_pph3cl2": ("Pd(PPh3)2Cl2", "bis(triphenylphosphine)palladium(II) dichloride"),
    "m_pph3": ("PPh3", "triphenylphosphine"),
    "m_pcy3": ("PCy3", "tricyclohexylphosphine"),
    "m_xantphos": ("Xantphos", "Xantphos"),
    "m_dppp": ("DPPP", "1,3-bis(diphenylphosphino)propane"),
    "m_dppb": ("DPPB", "1,4-bis(diphenylphosphino)butane"),
    "m_dioxane": ("Dioxane", "dioxane"),
    "m_mecn": ("MeCN", "acetonitrile"),
    "m_toluene": ("Toluene", "toluene"),
    "m_dmso": ("DMSO", "dimethyl sulfoxide"),
    "m_dmf": ("DMF", "N,N-dimethylformamide"),
    "m_pe": ("petroleum ether", "petroleum ether"),
    "m_etoac": ("ethyl acetate", "ethyl acetate"),
    "m_thf": ("THF", "tetrahydrofuran"),
    "m_etoh": ("EtOH", "ethanol"),
    "m_nbs": ("NBS", "N-bromosuccinimide"),
    "m_phco2h": ("PhCOOH", "benzoic acid"),
    "m_bnh2": ("BnNH2", "benzylamine"),
    "m_i2": ("I2", "iodine"),
    "m_tbhp": ("TBHP", "tert-butyl hydroperoxide"),
    "m_alpha_keto_ester": ("alpha-keto ester (drawn)", "alpha-keto ester reagent shown in Scheme 4"),
    "m_thiazolium": ("thiazolium catalyst", "thiazolium catalyst shown in Scheme 4"),
    "m_zn": ("Zn", "zinc powder"),
    "m_nh4cl": ("NH4Cl(aq)", "aqueous ammonium chloride"),
}
for mid, (label, name) in common.items():
    ev = (P[2], P[3], P[6]) if mid in {"m_tfben", "m_et3n", "m_pd_acac2", "m_dppp", "m_toluene", "m_pe", "m_etoac"} else (P[2], P[3])
    molecule(mid, label, name, ev, chirality="unresolved", geometry="not_applicable", rep_status="not_encoded")

# Allene scope inputs and products.  Source drawings are retained without forcing definitive SMILES.
scope2 = [
    ("a", 87, 24), ("b", 83, 24), ("c", 74, 24), ("d", 69, 24),
    ("e", 87, 24), ("f", 75, 24), ("g", 95, 24), ("h", 84, 24),
    ("i", 90, 24), ("j", 86, 24), ("k", 85, 24), ("l", 87, 24),
    ("m", 91, 24), ("n", 90, 24), ("o", 86, 24), ("p", 80, 24),
    ("q", 85, 24), ("r", 88, 36), ("s", 93, 24), ("t", 93, 36),
    ("u", 85, 36), ("v", 73, 36), ("w", 82, 36), ("x", 74, 24),
]
for letter, _, _ in scope2:
    if letter not in {"a", "k"}:
        molecule(f"m_2{letter}", f"2{letter}", f"allene 2{letter}", (P[3],),
                 chirality="unspecified", geometry="not_applicable", rep_status="not_encoded",
                 stereo_note="The allene drawing is visible but no axial configuration is assigned; substitution-dependent axial chirality is left unspecified.")
    if letter != "a":
        molecule(f"m_3a{letter}", f"3a{letter}", f"product 3a{letter}", (P[2], P[3]) if letter == "k" else (P[3],),
                 chirality="unresolved", geometry="unspecified", rep_status="not_encoded",
                 stereo_note="The product drawing contains an exocyclic alkene without a printed E/Z assignment.")
molecule("m_3ak_prime", "3ak'", "regioisomer 3ak'", (P[2],), chirality="unresolved", geometry="unspecified", rep_status="not_encoded")

# Aniline scope inputs/products.
scope3 = [("b", 85), ("c", 80), ("d", 90), ("e", 91), ("f", 93), ("g", 82),
          ("h", 87), ("i", 88), ("j", 93), ("k", 82), ("l", 82)]
for letter, _ in scope3:
    molecule(f"m_1{letter}", f"1{letter}", f"(N-SO2Py)-2-iodoaniline 1{letter}", (P[3],),
             chirality="unresolved", geometry="not_applicable", rep_status="not_encoded")
    molecule(f"m_3{letter}a", f"3{letter}a", f"product 3{letter}a", (P[3],),
             chirality="unresolved", geometry="unspecified", rep_status="not_encoded",
             stereo_note="The product drawing contains an exocyclic alkene without a printed E/Z assignment.")

# Derivatization and control compounds.
for label in ["4", "5", "6", "7"]:
    molecule(f"m_{label}", label, f"compound {label}", (P[3],),
             chirality="unresolved", geometry="unresolved", rep_status="not_encoded",
             stereo_note="The source drawing is retained, but exact stereochemical identity is not assigned in the main-text pages.")
for label, name in [
    ("8", "2-iodoaniline"), ("9", "N-methyl-2-iodoaniline"), ("10", "N-Boc-2-iodoaniline"),
    ("11", "N-acetyl-2-iodoaniline"), ("12", "aldehyde 12"), ("13", "target product 13"),
    ("14", "N-benzenesulfonyl-2-iodoaniline"), ("15", "cyclocarbonylation regioisomer 15"),
    ("15_prime", "cyclocarbonylation regioisomer 15'"), ("16", "N-picolinyl-2-iodoaniline"),
    ("17", "possible reaction intermediate 17"), ("18", "intramolecular Heck cyclization product 18"),
]:
    source_label = "15'" if label == "15_prime" else label
    molecule(f"m_{label}", source_label, name, (P[4],),
             chirality="unresolved", geometry="unresolved", rep_status="not_encoded")

generic_molecule("m_generic_1", "generic 1", "N-SO2Py-2-iodoaniline bearing R1", [("R1", "aryl-ring substituent")], (P[1], P[3], P[4]))
generic_molecule("m_generic_2", "generic 2", "substituted allene bearing R1/R2/R3", [("R1", "allene substituent"), ("R2", "allene substituent"), ("R3", "allene substituent")], (P[1], P[3], P[4]))
generic_molecule("m_generic_3", "generic 3", "cyclocarbonylation product bearing R1/R2/R3", [("R1", "product substituent"), ("R2", "product substituent"), ("R3", "product substituent")], (P[1], P[3], P[4]))


def charge_fact(raw, number, unit, evidence, method="original_extraction"):
    return st.fact(raw, {"number": number, "unit": unit, "raw_text": raw}, evidence, method=method)


def participant(pid, mid, role, evidence, stage="reaction", amount=None, equiv=None,
                molpct=None, volume=None, ref=None, product_role=None, inherited=False):
    p = st.new_participant(pid, mid, role, evidence, stage=stage)
    method = "inherited_general_condition" if inherited else "original_extraction"
    if amount is not None:
        p["charges"]["amount"] = charge_fact(f"{amount} mmol", amount, "mmol", evidence, method)
    if equiv is not None:
        p["charges"]["equivalents"] = charge_fact(f"{equiv} equiv.", equiv, "equiv", evidence, method)
    if molpct is not None:
        p["charges"]["mol_percent"] = charge_fact(f"{molpct} mol%", molpct, "mol%", evidence, method)
    if volume is not None:
        p["charges"]["volume"] = charge_fact(f"{volume} mL", volume, "mL", evidence, method)
    if ref is not None:
        p["equivalence_reference"] = ext(ref, ref, evidence, method)
    if product_role is not None:
        p["product_role"] = ext(product_role, product_role, evidence)
    return p


def base_reaction(rid, evidence, record_type="specific_experiment", categories=("other",), location=None):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = ext(list(categories), list(categories), evidence)
    r["location"] = ext(location or {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence})},
                                location or {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence})}, evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": ["Machine-extracted record; chemistry and structure assignments require expert review."]}
    return r


def add_step(r, sid, participant_ids, operation=None, temperature=None, time=None,
             atmosphere=None, evidence=(), stage="reaction", inherited=False):
    s = st.new_step(sid, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    method = "inherited_general_condition" if inherited else "original_extraction"
    if operation is not None:
        s["operation"] = ext(operation, operation, evidence, method)
    if temperature is not None:
        raw, norm = temperature
        s["temperature"] = st.fact(raw, norm, evidence, method=method)
    if time is not None:
        raw, norm = time
        s["time"] = st.fact(raw, norm, evidence, method=method)
    if atmosphere is not None:
        s["atmosphere"] = ext(atmosphere, atmosphere, evidence, method)
    r["steps"].append(s)


def add_outcome(r, oid, kind, targets, raw=None, norm=None, evidence=(), yield_type=None,
                basis=None, scope=None):
    o = st.new_outcome(oid, kind, targets, evidence)
    if raw is not None or norm is not None:
        o["value"] = ext(raw, norm if norm is not None else raw, evidence)
    if targets:
        o["product_role"] = ext("target product" if kind in {"yield", "product_amount"} else "measured entity",
                                "target product" if kind in {"yield", "product_amount"} else "measured entity", evidence)
    if yield_type is not None:
        o["yield_type"] = ext(yield_type, yield_type, evidence)
    if basis is not None:
        o["measurement_basis"] = ext(basis, basis, evidence)
    if scope is not None:
        o["scope"] = ext(scope, scope, evidence)
    r["outcomes"].append(o)


def classify(r, changes, evidence, coverage="partial", limitation=None):
    c = r["classification"]
    c["mapped_reaction_smiles"] = missing()
    labels = {}
    c["bond_changes"] = []
    for change in changes:
        c["bond_changes"].append(derived(change, evidence, "Net bond change inferred from the displayed reactant/product drawings; exact atom mapping is unresolved."))
        if change["change_type"] in {"formed", "order_increased", "hydrogen_gain"}:
            lab = st.bond_label(change["atoms"][0]["element"], change["atoms"][1]["element"], change["after_order"])
            labels.setdefault(lab["code"], dict(lab, supporting_change_ids=[]))["supporting_change_ids"].append(change["change_id"])
    c["labels"] = derived(list(labels.values()), evidence, "Taxonomy labels were calculated from the explicitly recorded partial net bond changes.")
    c["coverage"] = inferred(coverage, evidence, "Exact atom mapping and hydrogen accounting are unresolved; classification is intentionally partial.")
    c["limitations"] = [limitation or "Partial net classification from the source drawing; no complete atom mapping is asserted."]


def classify_main(r, evidence):
    classify(r, [
        {"change_id": "cc1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "cc2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
        {"change_id": "cn1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
    ], evidence, limitation="Cyclocarbonylation is assigned two net C-C formations and one C-N ring closure from the displayed product; exact atom origins and implicit-H changes remain unresolved.")


def classify_na(r, evidence, reason):
    c = r["classification"]
    c["labels"] = inferred([], evidence, reason)
    c["coverage"] = inferred("not_applicable", evidence, reason)
    c["limitations"] = [reason]


def classify_unresolved(r, evidence, reason):
    c = r["classification"]
    c["labels"] = inferred([], evidence, reason)
    c["coverage"] = inferred("unresolved", evidence, reason)
    c["limitations"] = [reason]


def add_scope_participants(r, aniline_mid, allene_mid, product_mid, evidence, hours=24, include_full_procedure=True):
    inherited = True
    r["participants"] = [
        participant("p_aniline", aniline_mid, "reactant", evidence, amount=0.5, equiv=1.0, ref=aniline_mid, inherited=inherited),
        participant("p_allene", allene_mid, "reactant", evidence, amount=1.25, equiv=2.5, ref=aniline_mid, inherited=inherited),
        participant("p_tfben", "m_tfben", "reagent", evidence, amount=0.5, equiv=1.0, ref=aniline_mid, inherited=inherited),
        participant("p_pd", "m_pd_acac2", "catalyst", evidence, amount=0.05, molpct=10, ref=aniline_mid, inherited=inherited),
        participant("p_dppp", "m_dppp", "ligand", evidence, amount=0.05, molpct=10, ref=aniline_mid, inherited=inherited),
        participant("p_et3n", "m_et3n", "reagent", evidence, amount=1.0, equiv=2.0, ref=aniline_mid, inherited=inherited),
        participant("p_toluene", "m_toluene", "solvent", evidence, volume=2.5, inherited=inherited),
        participant("p_product", product_mid, "product", evidence, product_role="target product"),
    ]
    if include_full_procedure:
        r["participants"].extend([
            participant("p_pe", "m_pe", "purification_agent", (P[6],), stage="purification", inherited=True),
            participant("p_etoac", "m_etoac", "purification_agent", (P[6],), stage="purification", inherited=True),
        ])
    r["reaction_scale"] = ext("1 or 1a (0.5 mmol)", {"reference_molecule_id": aniline_mid, "amount": {"number": 0.5, "unit": "mmol"}}, evidence, "inherited_general_condition")
    r["solvent_conditions"]["composition"] = ext("toluene", ["m_toluene"], evidence, "inherited_general_condition")
    r["solvent_conditions"]["total_volume"] = charge_fact("2.5 mL", 2.5, "mL", evidence, "inherited_general_condition")
    reaction_pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    add_step(r, "reaction", reaction_pids,
             operation="Oven-dried 15 mL tube; evacuate/refill with nitrogen three times; add allene, Et3N and toluene by syringe; seal and stir.",
             temperature=("90 °C", {"number": 90, "unit": "°C", "raw_text": "90 °C"}),
             time=(f"{hours} h", {"number": hours, "unit": "h", "raw_text": f"{hours} h"}),
             atmosphere="nitrogen", evidence=evidence, inherited=True)
    if include_full_procedure:
        add_step(r, "workup", [], operation="After completion, concentrate the resulting mixture under vacuum.", evidence=(P[6],), stage="workup", inherited=True)
        add_step(r, "purification", ["p_pe", "p_etoac"], operation="Silica gel column chromatography (petroleum ether/ethyl acetate = 5:1).", evidence=(P[6],), stage="purification", inherited=True)
        r["procedure"]["general_procedure_reference"] = ext("Experimental, General procedure", {"page": 6, "section": "General procedure"}, (P[6],), "inherited_general_condition")
        r["procedure"]["addition_order"] = ext("1/Pd(acac)2/DPPP/TFBen, N2 cycling, then allene/Et3N/toluene by syringe", ["1", "Pd(acac)2", "DPPP", "TFBen", "N2 cycling", "allene", "Et3N", "toluene"], (P[6],), "inherited_general_condition")
        r["procedure"]["pretreatment"] = ext("oven-dried tube; vacuum/refill nitrogen three times", "oven-dried tube; vacuum/refill nitrogen three times", (P[6],), "inherited_general_condition")


reactions = []

# Current-work generic scheme and proposed mechanism are kept distinct from experiments.
r = base_reaction("scheme1_current_generic", (P[1],), "generic_scheme", ("other",), {"page": 1, "scheme": "Scheme 1b"})
r["reaction_name"] = ext("Carbonylative cyclization of 2-iodoanilines with allenes", "Carbonylative cyclization of 2-iodoanilines with allenes", (P[1],))
r["transformation_type"] = inferred("palladium-catalyzed cyclocarbonylation", (P[1],), "Normalized from the title and Scheme 1 caption.")
r["participants"] = [participant("p1", "m_generic_1", "reactant", (P[1],)), participant("p2", "m_generic_2", "reactant", (P[1],)), participant("p3", "m_generic_3", "product", (P[1],), product_role="generic target")]
classify_main(r, (P[1],))
reactions.append(r)

r = base_reaction("scheme8_proposed_mechanism", (P[4], P[5]), "proposed_mechanism", ("mechanistic_study",), {"pages": [4, 5], "scheme": "Scheme 8 and Figure 1"})
r["reaction_name"] = ext("Proposed mechanism", "Proposed mechanism", (P[4],))
r["transformation_type"] = ext("palladium-catalyzed cyclocarbonylation", "palladium-catalyzed cyclocarbonylation", (P[4],))
r["participants"] = [participant("p1", "m_generic_1", "reactant", (P[4],)), participant("p2", "m_generic_2", "reactant", (P[4],)), participant("p3", "m_generic_3", "product", (P[4],), product_role="mechanistic product")]
r["notes"] = [ext("Path I is calculated to be favored; allene migratory insertion is the rate-determining step (29.3 kcal mol−1).", {"favored_path": "I", "rate_determining_step": "allene migratory insertion", "barrier_kcal_mol": 29.3}, (P[4], P[5]))]
classify_main(r, (P[4], P[5]))
reactions.append(r)

# Table 1 optimization entries. Entry 16 is the same experiment/results as Scheme 2 product 3ak.
table = [
    (1, "m_pd_oac2", "m_pph3", "m_dioxane", "nr", None, None),
    (2, "m_pd_oac2", "m_pcy3", "m_dioxane", "nr", None, None),
    (3, "m_pd_oac2", "m_xantphos", "m_dioxane", "trace", None, None),
    (4, "m_pd_oac2", "m_dppp", "m_dioxane", "yield", 35, "1:1"),
    (5, "m_pd_oac2", "m_dppb", "m_dioxane", "trace", None, None),
    (6, "m_pd_oac2", "m_dppp", "m_mecn", "yield", 41, "5:1"),
    (7, "m_pd_oac2", "m_dppp", "m_toluene", "yield", 53, "3ak only"),
    (8, "m_pd_oac2", "m_dppp", "m_dmso", "trace", None, None),
    (9, "m_pd_oac2", "m_dppp", "m_dmf", "trace", None, None),
    (10, "m_pd_cl2", "m_dppp", "m_toluene", "yield", 40, "3ak only"),
    (11, "m_pd_tfa2", "m_dppp", "m_toluene", "yield", 55, "3ak only"),
    (12, "m_pd_acac2", "m_dppp", "m_toluene", "yield", 61, "3ak only"),
    (13, "m_pd_pph3cl2", "m_dppp", "m_toluene", "yield", 34, "3ak only"),
    (14, "m_pd_acac2", "m_dppp", "m_toluene", "yield", 87, "3ak only"),
    (15, "m_pd_acac2", "m_dppp", "m_toluene", "yield", 70, "3ak only"),
    (16, "m_pd_acac2", "m_dppp", "m_toluene", "yield", 85, "3ak only"),
    (17, "m_pd_acac2", "m_dppp", "m_toluene", "yield", 47, "3ak only"),
]
for entry, pd, lig, solv, result_kind, value, ratio in table:
    ev = (P[2], P[3], P[6]) if entry == 16 else (P[2],)
    rid = "table1_entry16_scope_3ak" if entry == 16 else f"table1_entry{entry:02d}"
    cats = ("optimization", "substrate_scope") if entry == 16 else ("optimization",)
    r = base_reaction(rid, ev, "specific_experiment", cats, {"page": 2, "table": 1, "entry": entry, "duplicate_scope_location": "Scheme 2, 3ak" if entry == 16 else None})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["reaction_name"] = ext("Screening of the reaction conditions", "palladium-catalyzed cyclocarbonylation", ev)
    r["transformation_type"] = inferred("cyclocarbonylation/annulation", ev, "Normalized from the table caption and product drawing.")
    ligand_molpct = 20 if lig in {"m_pph3", "m_pcy3"} else 10
    pd_molpct = 1 if entry == 15 else 10
    tf_equiv = 1.0 if entry == 16 else 0.5 if entry == 17 else 4.0
    solvent_volume = 2.5 if entry >= 14 else 1.0
    products = [] if result_kind == "nr" else [participant("p_3ak", "m_3ak", "product", ev, product_role="target regioisomer")]
    if ratio in {"1:1", "5:1"}:
        products.append(participant("p_3ak_prime", "m_3ak_prime", "product", ev, product_role="regioisomer"))
    r["participants"] = [
        participant("p_1a", "m_1a", "reactant", ev, amount=0.5, equiv=1.0, ref="m_1a"),
        participant("p_2k", "m_2k", "reactant", ev, amount=1.25 if entry == 16 else None, equiv=2.5, ref="m_1a"),
        participant("p_tfben", "m_tfben", "reagent", ev, amount=0.5 if entry == 16 else None, equiv=tf_equiv, ref="m_1a"),
        participant("p_pd", pd, "catalyst", ev, amount=0.05 if entry == 16 else None, molpct=pd_molpct, ref="m_1a"),
        participant("p_lig", lig, "ligand", ev, amount=0.05 if entry == 16 else None, molpct=ligand_molpct, ref="m_1a"),
        participant("p_et3n", "m_et3n", "reagent", ev, amount=1.0 if entry == 16 else None, equiv=2.0, ref="m_1a"),
        participant("p_solvent", solv, "solvent", ev, volume=solvent_volume),
    ] + products
    r["reaction_scale"] = ext("1a (0.5 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 0.5, "unit": "mmol"}}, ev)
    r["solvent_conditions"]["composition"] = ext(molecules[solv]["names"]["original"]["normalized_value"], [solv], ev)
    r["solvent_conditions"]["total_volume"] = charge_fact(f"{solvent_volume} mL", solvent_volume, "mL", ev)
    pids = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
    add_step(r, "reaction", pids, operation="Screening experiment under Table 1 conditions",
             temperature=("90 °C", {"number": 90, "unit": "°C", "raw_text": "90 °C"}),
             time=("24 h", {"number": 24, "unit": "h", "raw_text": "24 h"}), evidence=ev)
    if entry == 16:
        r["participants"].extend([
            participant("p_pe", "m_pe", "purification_agent", (P[6],), stage="purification", inherited=True),
            participant("p_etoac", "m_etoac", "purification_agent", (P[6],), stage="purification", inherited=True),
        ])
        r["steps"][0]["atmosphere"] = ext("nitrogen", "nitrogen", (P[6],), "inherited_general_condition")
        add_step(r, "workup", [], operation="After completion, concentrate the resulting mixture under vacuum.", evidence=(P[6],), stage="workup", inherited=True)
        add_step(r, "purification", ["p_pe", "p_etoac"], operation="Silica gel column chromatography (petroleum ether/ethyl acetate = 5:1).", evidence=(P[6],), stage="purification", inherited=True)
        r["procedure"]["general_procedure_reference"] = ext("Experimental, General procedure", {"page": 6, "section": "General procedure"}, (P[6],), "inherited_general_condition")
        r["procedure"]["addition_order"] = ext("1a/Pd(acac)2/DPPP/TFBen, N2 cycling, then 2k/Et3N/toluene by syringe", ["1a", "Pd(acac)2", "DPPP", "TFBen", "N2 cycling", "2k", "Et3N", "toluene"], (P[6],), "inherited_general_condition")
        r["procedure"]["pretreatment"] = ext("oven-dried tube; vacuum/refill nitrogen three times", "oven-dried tube; vacuum/refill nitrogen three times", (P[6],), "inherited_general_condition")
        r["notes"].append(ext("Table 1 entry 16 and Scheme 2 product 3ak report the same substrates, conditions, 85% isolated yield, and single regioisomer; represented once with both locations.", "duplicate source locations merged", (P[2], P[3])))
    if result_kind == "nr":
        add_outcome(r, "o_nr", "no_reaction", [], "n.r.", "no reaction", ev, scope="Table 1 entry")
        classify_na(r, ev, "No observed reaction; no intended bond-formation labels are asserted.")
    elif result_kind == "trace":
        add_outcome(r, "o_trace", "trace", ["m_3ak"], "Trace", "trace", ev, scope="Table 1 yield/ratio cell")
        classify_main(r, ev)
    else:
        targets = ["m_3ak", "m_3ak_prime"] if ratio in {"1:1", "5:1"} else ["m_3ak"]
        add_outcome(r, "o_yield", "yield", targets, f"{value}%", {"number": value, "unit": "%", "raw_text": f"{value}%"}, ev,
                    yield_type="isolated yield", basis="isolated", scope="Table 1 entry")
        if ratio:
            add_outcome(r, "o_ratio", "selectivity", targets, ratio, ratio, ev, scope="3ak:3ak' regioisomer ratio")
        classify_main(r, ev)
    reactions.append(r)

# Scheme 2 allene scope; 3ak is already represented by Table 1 entry 16.
for letter, yld, hours in scope2:
    if letter == "k":
        continue
    allene_mid = f"m_2{letter}"
    product_mid = "m_3aa" if letter == "a" else f"m_3a{letter}"
    rid = f"scheme2_3a{letter}"
    r = base_reaction(rid, (P[3], P[6]), "specific_experiment", ("substrate_scope",), {"page": 3, "scheme": 2, "product": f"3a{letter}"})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["reaction_name"] = ext("Scope of allenes", "palladium-catalyzed cyclocarbonylation", (P[3],))
    r["transformation_type"] = inferred("regioselective cyclocarbonylation/annulation", (P[3],), "Normalized from Scheme 2.")
    add_scope_participants(r, "m_1a", allene_mid, product_mid, (P[3], P[6]), hours)
    add_outcome(r, "o_yield", "yield", [product_mid], f"{yld}%", {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, (P[3],), yield_type="isolated yield", basis="isolated", scope=f"Scheme 2 product 3a{letter}")
    classify_main(r, (P[3],))
    reactions.append(r)

# Scheme 3 aniline scope.
for letter, yld in scope3:
    aniline_mid, product_mid = f"m_1{letter}", f"m_3{letter}a"
    r = base_reaction(f"scheme3_3{letter}a", (P[3], P[6]), "specific_experiment", ("substrate_scope",), {"page": 3, "scheme": 3, "product": f"3{letter}a"})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["reaction_name"] = ext("Scope of (N-SO2Py)-2-iodoanilines", "palladium-catalyzed cyclocarbonylation", (P[3],))
    r["transformation_type"] = inferred("regioselective cyclocarbonylation/annulation", (P[3],), "Normalized from Scheme 3.")
    add_scope_participants(r, aniline_mid, "m_2a", product_mid, (P[3], P[6]), 24)
    add_outcome(r, "o_yield", "yield", [product_mid], f"{yld}%", {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, (P[3],), yield_type="isolated yield", basis="isolated", scope=f"Scheme 3 product 3{letter}a")
    classify_main(r, (P[3],))
    reactions.append(r)

# Scheme 4 and Scheme 5 derivatizations.
derivs = [
    ("scheme4_eq1_3aa_to_4", "m_3aa", "m_4", [("p_nbs", "m_nbs", "reagent", 3.0), ("p_acid", "m_phco2h", "reagent", 0.4)], "m_thf", 50, 24, 90, "dibromination"),
    ("scheme4_eq2_3aa_to_5", "m_3aa", "m_5", [("p_bnh2", "m_bnh2", "reagent", 6.0), ("p_i2", "m_i2", "reagent", 0.2), ("p_tbhp", "m_tbhp", "reagent", 4.0)], "m_etoac", "r.t.", 12, 75, "aziridination"),
    ("scheme4_eq3_3aa_to_6", "m_3aa", "m_6", [("p_pyruvate", "m_alpha_keto_ester", "reactant", 2.0), ("p_cat", "m_thiazolium", "catalyst", 0.4), ("p_et3n", "m_et3n", "reagent", 1.0)], "m_etoh", 50, 24, 87, "Stetter reaction"),
    ("scheme5_6_to_7", "m_6", "m_7", [("p_zn", "m_zn", "reagent", None), ("p_nh4cl", "m_nh4cl", "reagent", None)], "m_thf", 60, 24, 78, "removal of the N-SO2Py directing group"),
]
for rid, start, prod, reagents, solvent, temp, hours, yld, rxname in derivs:
    scheme = 5 if rid.startswith("scheme5") else 4
    r = base_reaction(rid, (P[3],), "specific_experiment", ("derivatization",), {"page": 3, "scheme": scheme, "arrow": rid})
    r["reaction_name"] = ext(rxname, rxname, (P[3],))
    r["transformation_type"] = inferred(rxname, (P[3],), "Normalized from the Scheme 4/5 drawing and accompanying prose.")
    r["participants"] = [participant("p_start", start, "reactant", (P[3],), equiv=1.0, ref=start), participant("p_product", prod, "product", (P[3],), product_role="target product")]
    for pid, mid, role, eq in reagents:
        p = participant(pid, mid, role, (P[3],), equiv=eq, ref=start) if eq is not None else participant(pid, mid, role, (P[3],))
        if mid == "m_zn":
            p["charges"]["equivalents"] = ext("excess", "excess", (P[3],))
            p["equivalence_reference"] = ext(start, start, (P[3],))
        r["participants"].append(p)
    r["participants"].append(participant("p_solvent", solvent, "solvent", (P[3],)))
    if rid.startswith("scheme5"):
        r["solvent_conditions"]["composition"] = ext("THF/NH4Cl(aq)", ["m_thf", "m_nh4cl"], (P[3],))
        r["solvent_conditions"]["mixture_ratio"] = ext("1:1", "1:1", (P[3],))
    else:
        r["solvent_conditions"]["composition"] = ext(molecules[solvent]["names"]["original"]["normalized_value"], [solvent], (P[3],))
    tnorm = temp if isinstance(temp, str) else {"number": temp, "unit": "°C", "raw_text": f"{temp} °C"}
    traw = temp if isinstance(temp, str) else f"{temp} °C"
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], operation=rxname,
             temperature=(traw, tnorm), time=(f"{hours} h", {"number": hours, "unit": "h", "raw_text": f"{hours} h"}), evidence=(P[3],))
    add_outcome(r, "o_yield", "yield", [prod], f"{yld}%", {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, (P[3],), scope=f"Scheme {scheme} arrow")
    if rid.endswith("to_4"):
        classify(r, [
            {"change_id": "cbr1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "cbr2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "Br"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "cc_down", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ], (P[3],))
    elif rid.endswith("to_5"):
        classify(r, [
            {"change_id": "cn1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "cn2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "N"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "cc_down", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ], (P[3],))
    elif rid.endswith("to_6"):
        classify(r, [
            {"change_id": "cc1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"},
            {"change_id": "ch1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "cc_down", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 2, "after_order": 1, "change_type": "order_decreased"},
        ], (P[3],))
    else:
        classify(r, [
            {"change_id": "nh1", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "H"}], "before_order": 0, "after_order": 1, "change_type": "hydrogen_gain", "hydrogen_count": 1},
            {"change_id": "ns_break", "atoms": [{"atom_map": None, "element": "N"}, {"atom_map": None, "element": "S"}], "before_order": 1, "after_order": 0, "change_type": "broken"},
        ], (P[3],))
    r["review"]["unresolved"].append("The main-text scheme reports the yield but not its measurement/isolation basis; measurement_basis remains unreported.")
    reactions.append(r)

# Scheme 6 directing-group controls under the inherited standard conditions.
controls = [
    (1, "m_8", "complicated mixtures"), (2, "m_9", "complicated mixtures"), (3, "m_10", "complicated mixtures"),
    (4, "m_11", "aldehyde"), (5, "m_14", "regioisomers"), (6, "m_16", "nr"),
]
for eq, start, result in controls:
    r = base_reaction(f"scheme6_eq{eq}", (P[4], P[6]), "specific_experiment", ("control", "mechanistic_study"), {"page": 4, "scheme": 6, "equation": eq})
    r["general_scheme_ids"] = [f"{DOC}:scheme1_current_generic"]
    r["reaction_name"] = ext("Effect of directing groups", "directing-group control for cyclocarbonylation", (P[4],))
    r["transformation_type"] = inferred("directing-group control", (P[4],), "Normalized from Scheme 6 caption.")
    # Use standard quantities, with the altered aniline as the 0.5 mmol reference.
    product_mid = "m_12" if eq == 4 else "m_15" if eq == 5 else "m_3aa"
    add_scope_participants(r, start, "m_2a", product_mid, (P[4], P[6]), 24, include_full_procedure=False)
    if result == "complicated mixtures":
        # Remove the intended product participant: no product identity was assigned to the mixture.
        r["participants"] = [p for p in r["participants"] if p["role"] != "product"]
        r["steps"][0]["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
        add_outcome(r, "o_mixture", "other", [], "complicated mixtures", "complicated mixtures", (P[4],), scope=f"Scheme 6 eq. ({eq})")
        classify_unresolved(r, (P[4],), "Complicated mixtures were reported without assigned product structures; net bond changes cannot be classified.")
    elif result == "nr":
        r["participants"] = [p for p in r["participants"] if p["role"] != "product"]
        r["steps"][0]["participant_ids"] = [p["participant_id"] for p in r["participants"] if p["stage"] == "reaction"]
        add_outcome(r, "o_nr", "no_reaction", [], "N.R.", "no reaction", (P[4],), scope="Scheme 6 eq. (6)")
        classify_na(r, (P[4],), "No reaction was observed; no intended bond-formation labels are asserted.")
    elif result == "aldehyde":
        r["participants"].append(participant("p_13", "m_13", "product", (P[4],), product_role="intended but not observed cyclocarbonylation product"))
        add_outcome(r, "o_yield12", "yield", ["m_12"], "80%", {"number": 80, "unit": "%", "raw_text": "80%"}, (P[4],), scope="aldehyde 12")
        add_outcome(r, "o_13_not_observed", "no_reaction", ["m_13"], "not observed", "not observed", (P[4],), scope="target product 13")
        classify_unresolved(r, (P[4],), "Aldehyde 12 is an alternative observed product, but exact atom origins/net changes are not safely resolved from the main-text raster.")
        r["review"]["unresolved"].append("Yield measurement/isolation basis for aldehyde 12 is not reported on pages 1–6.")
    else:
        r["participants"].append(participant("p_15prime", "m_15_prime", "product", (P[4],), product_role="regioisomer"))
        add_outcome(r, "o_yield", "yield", ["m_15", "m_15_prime"], "31%", {"number": 31, "unit": "%", "raw_text": "31%"}, (P[4],), scope="combined regioisomers 15 and 15'")
        add_outcome(r, "o_ratio", "selectivity", ["m_15", "m_15_prime"], "6:1", "6:1", (P[4],), scope="15:15' regioisomer ratio")
        classify_main(r, (P[4],))
        r["review"]["unresolved"].append("Yield and regioisomer-ratio measurement bases are not reported on pages 1–6.")
    reactions.append(r)

# Scheme 7: with and without phenylallene are separate control experiments.
for with_allene, yld in [(False, 30), (True, 80)]:
    suffix = "with_2a" if with_allene else "without_2a"
    r = base_reaction(f"scheme7_{suffix}", (P[4], P[6]), "specific_experiment", ("control", "mechanistic_study"), {"page": 4, "scheme": 7, "variant": suffix})
    r["reaction_name"] = ext("Control experiments", "intramolecular Heck cyclization control", (P[4],))
    r["transformation_type"] = inferred("intramolecular Heck cyclization", (P[4],), "Author identifies product 18 as the intramolecular Heck cyclization product.")
    r["participants"] = [participant("p_17", "m_17", "reactant", (P[4],)), participant("p_18", "m_18", "product", (P[4],), product_role="observed product"), participant("p_3aa", "m_3aa", "product", (P[4],), product_role="not observed cyclocarbonylation product")]
    if with_allene:
        r["participants"].append(participant("p_2a", "m_2a", "reactant", (P[4],), equiv=2.5, ref="m_17"))
    # Standard-condition agents are retained, but no unsupported absolute charge for 17 is invented.
    for pid, mid, role, eq, molpct in [("p_tfben", "m_tfben", "reagent", 1.0, None), ("p_pd", "m_pd_acac2", "catalyst", None, 10), ("p_dppp", "m_dppp", "ligand", None, 10), ("p_et3n", "m_et3n", "reagent", 2.0, None)]:
        p = participant(pid, mid, role, (P[4], P[6]), equiv=eq, molpct=molpct, ref="m_17", inherited=True)
        r["participants"].append(p)
    r["participants"].append(participant("p_toluene", "m_toluene", "solvent", (P[4], P[6]), volume=2.5, inherited=True))
    add_step(r, "reaction", [p["participant_id"] for p in r["participants"]], operation="Apply the standard reaction conditions to intermediate 17" + (" with phenylallene 2a" if with_allene else " without phenylallene 2a"),
             temperature=("90 °C", {"number": 90, "unit": "°C", "raw_text": "90 °C"}), time=("24 h", {"number": 24, "unit": "h", "raw_text": "24 h"}), atmosphere="nitrogen", evidence=(P[4], P[6]), inherited=True)
    add_outcome(r, "o_yield18", "yield", ["m_18"], f"{yld}%", {"number": yld, "unit": "%", "raw_text": f"{yld}%"}, (P[4],), scope=f"product 18, {suffix}")
    add_outcome(r, "o_3aa_not_observed", "no_reaction", ["m_3aa"], "not observed", "not observed", (P[4],), scope="cyclocarbonylated product 3aa")
    classify(r, [{"change_id": "cc1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"}], (P[4],), limitation="Only the intramolecular C-C ring closure is classified; alkene migration/H changes and exact atom mapping remain unresolved.")
    r["review"]["unresolved"].append("Yield measurement/isolation basis for product 18 is not reported on pages 1–6.")
    reactions.append(r)

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:bibliography",
        "kind": "bibliographic_metadata",
        "title": "Palladium-catalyzed directing group assisted and regioselectivity reversed cyclocarbonylation of arylallenes with 2-iodoanilines",
        "authors": ["Jian-Shu Wang", "Lingyun Yao", "Jun Ying", "Xiaoling Luo", "Xiao-Feng Wu"],
        "journal": "Organic Chemistry Frontiers",
        "year": 2021,
        "volume": 8,
        "pages": "792–798",
        "doi": "10.1039/d0qo01404d",
        "evidence_ids": [P[1]],
    },
    {
        "context_id": f"{DOC}:si_notice",
        "kind": "supplementary_information_notice",
        "observation": "Page 1 states that ESI is available, but no SI file was supplied or explicitly associated in this workspace.",
        "evidence_ids": [P[1]],
    },
    {
        "context_id": f"{DOC}:scheme1_prior_work",
        "kind": "cited_prior_work",
        "observation": "Scheme 1a reproduces Alper and Ye's prior ionic-liquid-promoted cyclocarbonylation and example yields (80%, 29%, 55%); these are not experiments performed in the current paper.",
        "evidence_ids": [P[1]],
    },
    {
        "context_id": f"{DOC}:xray_3aa",
        "kind": "characterization",
        "observation": "The connectivity of 3aa is stated to be confirmed by X-ray crystallography, CCDC 2023302.",
        "evidence_ids": [P[2], P[3]],
    },
    {
        "context_id": f"{DOC}:dft_figures",
        "kind": "computational_mechanistic_context",
        "observation": "Figures 1–2 and Table 2 report DFT energy profiles, regioselectivity analysis, transition-state geometries, and distortion/interaction energies; they are not synthetic experiments.",
        "evidence_ids": [P[5], P[6]],
    },
    {
        "context_id": f"{DOC}:page6_references",
        "kind": "bibliography_only_region",
        "observation": "The right column of page 6 begins acknowledgements and cited references; cited reaction titles were not imported as this paper's experiments.",
        "evidence_ids": [P[6]],
    },
]

# Graph diagnostics are retained separately; they do not promote source verification.
diagnostics = {}
for mid in ("m_1a", "m_2a", "m_2k", "m_3aa"):
    diagnostics[mid] = st.stereo_inspect(st.value(molecules[mid]["smiles"]))

(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/stereo_diagnostics.json").write_text(json.dumps(diagnostics, ensure_ascii=False, indent=2) + "\n")
(RUN / "raw_model/record_inventory.json").write_text(json.dumps({
    "document_id": DOC,
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions),
    "formal_records": len(reactions),
    "molecules": len(molecules),
    "record_ids": [r["reaction_id"] for r in reactions],
    "notes": [
        "Table 1 entry 16 and Scheme 2 product 3ak are represented once with both source locations.",
        "Cited prior work, SI notice, X-ray statement, and DFT-only regions are source_context, not current-paper synthesis experiments.",
    ],
}, ensure_ascii=False, indent=2) + "\n")

page_regions = {
    1: [
        {"region": "title/abstract/bibliographic header", "finding": "Current article metadata and abstract; ESI availability is stated, but no SI file is supplied.", "evidence_ids": [P[1]]},
        {"region": "Scheme 1a", "finding": "Cited Alper/Ye prior work with example outcomes; retained as cited-prior-work context, not current-paper experiments.", "evidence_ids": [P[1]]},
        {"region": "Scheme 1b", "finding": "Generic current-work cyclocarbonylation scheme; recorded separately as a generic scheme.", "evidence_ids": [P[1]]},
    ],
    2: [
        {"region": "Table 1 including footnotes a-e", "finding": "All 17 optimization rows inventoried independently; n.r., trace, isolated yields, regioisomer ratios, solvent-volume and catalyst/TFBen exceptions retained.", "evidence_ids": [P[2]]},
        {"region": "right-column scope discussion", "finding": "Narrative cross-check of Scheme 2/3 yields and regioselectivity; exact product-specific values are taken from the visible Scheme drawings on page 3.", "evidence_ids": [P[2]]},
        {"region": "right-column Scheme 4 discussion", "finding": "Text cross-check for three product-3aa derivatizations and the page-3 arrows; no additional independent experiment was created.", "evidence_ids": [P[2]]},
    ],
    3: [
        {"region": "Scheme 2", "finding": "Twenty-four allene-scope product entries (3aa-3ax) inventoried. Product 3ak duplicates Table 1 entry 16 and is represented once with both sources; footnote-b 36 h is preserved for 3ar and 3at-3aw.", "evidence_ids": [P[3]]},
        {"region": "Scheme 3", "finding": "Eleven (N-SO2Py)-2-iodoaniline scope entries (3ba-3la) inventoried with isolated yields.", "evidence_ids": [P[3]]},
        {"region": "Scheme 4", "finding": "Three distinct derivatization arrows (dibromination, aziridination, Stetter reaction) recorded with every printed reagent loading, solvent, temperature, time and yield.", "evidence_ids": [P[3]]},
        {"region": "Scheme 5", "finding": "One directing-group-removal experiment (6 to 7) recorded; exact structures and yield measurement basis remain unresolved from the supplied main-text pages.", "evidence_ids": [P[3]]},
    ],
    4: [
        {"region": "Scheme 6 equations 1-6", "finding": "Six directing-group controls recorded separately: three complicated-mixture outcomes, aldehyde 12 with target 13 not observed, regioisomers 15/15' (31%, 6:1), and one no-reaction result.", "evidence_ids": [P[4]]},
        {"region": "Scheme 7", "finding": "Intermediate 17 tested with and without 2a as two distinct experiments; product 18 yields 80% and 30%, respectively, while 3aa is not observed.", "evidence_ids": [P[4]]},
        {"region": "Scheme 8 and mechanism prose", "finding": "Author-proposed pathways retained as one proposed-mechanism record, not as experimentally proven reaction steps.", "evidence_ids": [P[4]]},
    ],
    5: [
        {"region": "Figure 1", "finding": "DFT free-energy profile and computed intermediates/transition states; computational mechanism context, not a synthetic experiment.", "evidence_ids": [P[5]]},
        {"region": "Figure 2 and surrounding prose", "finding": "Computed regioselectivity comparison; 33-ts is stated to be 2.5 kcal mol-1 above 27-ts. Retained as computational context.", "evidence_ids": [P[5]]},
    ],
    6: [
        {"region": "Table 2 and top-left discussion", "finding": "DFT distortion/interaction energies and allene angles; computational observations, not synthesis records.", "evidence_ids": [P[6]]},
        {"region": "Conclusions", "finding": "Summary only; no independent experiment added.", "evidence_ids": [P[6]]},
        {"region": "Experimental—General procedure", "finding": "Full printed charge, atmosphere, addition-order, vessel, temperature/time, concentration/workup and chromatography details inherited only by matching optimized/scope/standard-condition records.", "evidence_ids": [P[6]]},
        {"region": "Acknowledgements and references", "finding": "Non-experimental and bibliography-only content; cited titles were not imported as current-paper experiments.", "evidence_ids": [P[6]]},
    ],
}
page_unresolved = {
    1: ["ESI is mentioned but no supplementary file is explicitly associated; no SI content was inferred."],
    2: ["Structures in the dense scope discussion are resolved from page-3 drawings only; exact systematic names/CAS are not supplied."],
    3: ["Most complex scope structures remain drawing/label identities without definitive SMILES; product alkene E/Z assignments are not printed.", "Scheme 4/5 yield measurement bases are not stated on pages 1-6."],
    4: ["Measurement bases for Scheme 6/7 yields and regioisomer ratio are not stated.", "Complicated mixtures have no assigned product identities; exact net bond changes are unresolved."],
    5: ["Computed structures/energies do not establish experimental intermediates; SI-only path-II details are unavailable."],
    6: ["The general procedure gives no per-product characterization and cannot resolve complex structures or missing measurement bases without associated SI."],
}
coverage_pages = []
for page in range(1, 7):
    eid = P[page]
    coverage_pages.append({
        "page": page,
        "state": "attempted",
        "regions": page_regions[page],
        "record_ids": [r["reaction_id"] for r in reactions if eid in r["evidence_ids"]],
        "unresolved": page_unresolved[page],
    })
(RUN / "coverage/pages_0001_0006.json").write_text(json.dumps({
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": coverage_pages,
}, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(molecules), "reactions": len(reactions), "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in reactions)}, ensure_ascii=False))
