#!/usr/bin/env python3
"""Build the page-1--6 v3.1 machine submission from the visually transcribed facts."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

DOC = "doc_484c4de81c12e684"
RUN = ROOT / "runs/current"
RAW = json.loads((RUN / "raw_model/page1_6_extraction.json").read_text())
BUNDLE = json.loads((RUN / "submitted.json").read_text())


def ev(*pages: int) -> list[str]:
    ids: list[str] = []
    for page in pages:
        ids += [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]
    return ids


def src(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def derived(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="structure_derivation", derivation=derivation)


def inferred(normalized, evidence, derivation, raw=None):
    return st.fact(raw, normalized, evidence, method="inference", derivation=derivation)


def not_applicable(evidence=(), derivation="Field is not applicable to the represented source entity."):
    return st.fact(evidence=evidence, method="not_applicable", status="not_applicable", derivation=derivation)


def make_molecule(mid, label, evidence, name=None, stereo_kind="unresolved", geometry="unresolved",
                  representation="unresolved", note=None):
    m = st.new_molecule(mid, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = src(label, label, evidence)
    m["names"]["original"] = src(name, name, evidence) if name else st.fact(evidence=evidence)
    m["names"]["english"] = src(name, name, evidence) if name else st.fact(evidence=evidence)
    m["names"]["abbreviations"] = src(label, [label], evidence)
    m["names"]["aliases"] = st.fact(evidence=evidence)
    m["cas_numbers"] = st.fact(evidence=evidence)
    m["representation_type"] = "unresolved"
    m["smiles"] = st.fact(evidence=evidence)
    m["generic_template"] = None
    m["stereochemistry"]["chirality_status"] = inferred(
        stereo_kind, evidence, "Stereochemical status was assessed from the supplied source drawing/name; no absent label was converted to racemic."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        geometry, evidence, "Geometry status was assessed from the supplied source drawing and explicit E/Z wording."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        representation, evidence, "No concrete stereochemical SMILES is asserted for this unresolved source structure."
    )
    m["identity_note"] = src(note, note, evidence) if note else st.fact(evidence=evidence)
    m["external_properties"] = []
    return m


def set_product_stereo(m, evidence):
    m["stereochemistry"]["atom_centers"] = src(
        "A stereogenic carbon is marked with an asterisk in the product drawing; no R/S assignment is given.",
        [{"site": "asterisked ring carbon", "configuration": "unassigned"}], evidence
    )
    m["stereochemistry"]["double_bonds"] = src(
        "products 3 were assigned to be E; Figure 1 labels 3 as E-isomer only",
        [{"site": "exocyclic ylidene double bond", "configuration": "E"}], evidence
    )
    m["stereochemistry"]["assignment"] = src(
        "The E configuration was assigned by conversion to known 5 and NMR comparison; absolute configuration by X-ray was unsuccessful.",
        {"geometry": "E", "absolute_configuration": "unassigned"}, evidence
    )
    m["stereochemistry"]["chirality_status"] = inferred(
        "partially_specified", evidence,
        "The source marks a stereogenic center and reports ee, but does not assign its absolute configuration."
    )
    m["stereochemistry"]["geometry_status"] = src("E-isomer", "specified", evidence)
    m["stereochemistry"]["representation_status"] = inferred(
        "not_encoded", evidence, "The source-defined E geometry and unassigned stereocenter are retained textually; no candidate product SMILES is asserted."
    )


def make_allenoate(mid, label, name, smiles, evidence, key):
    m = make_molecule(mid, label, evidence, name, "achiral", "not_applicable", "not_applicable",
                      "Concrete connectivity candidate derived from the source drawing; parsing does not verify drawing recognition.")
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, evidence, "SMILES candidate transcribed from the source allenoate drawing.", raw=f"source drawing of {label}")
    m["stereochemistry"]["axial"] = derived(
        [{"axis": "C=C=C", "configuration": "not stereogenic because one terminus is CH2"}], evidence,
        "The terminal allene carbon bears two hydrogens in the source structure."
    )
    m["stereochemistry"]["assignment"] = inferred(
        {"allene_axis": "not stereogenic; terminal CH2"}, evidence,
        "Axial chirality requires two different substituents at both termini; the source shows CH2 at one end."
    )
    m["stereochemistry"]["chirality_status"] = inferred(
        "achiral", evidence, "The allene has a terminal CH2 group and no other stereogenic element in the depicted molecule."
    )
    m["stereochemistry"]["geometry_status"] = inferred(
        "not_applicable", evidence, "Ordinary E/Z geometry is not assigned to this terminal cumulene in the source."
    )
    m["stereochemistry"]["representation_status"] = inferred(
        "not_applicable", evidence, "No applicable source stereochemistry needs encoding for this achiral terminal allene candidate."
    )
    m["stereochemistry"]["identity_key"] = derived(
        key, evidence, "Generated by standard stereo from the candidate SMILES; graph diagnostic only."
    )
    return m


def make_template(mid, label, evidence, product=False, previous=False):
    m = make_molecule(mid, label, evidence, label, "partially_specified" if product else "achiral",
                      "specified" if product else "not_applicable", "not_encoded")
    m["representation_type"] = "generic_template"
    m["smiles"] = st.fact(evidence=evidence)
    vars_ = [
        {"label": "R1", "attachment_points": src("coumarin benzene ring", "coumarin benzene ring", evidence),
         "definition": src("R1 as drawn", "substituent(s) on coumarin benzene ring", evidence),
         "allowed_values": st.fact(evidence=evidence), "constraints": st.fact(evidence=evidence)},
        {"label": "R2" if previous else "Ar", "attachment_points": src("3-aroyl carbonyl substituent", "3-aroyl carbonyl substituent", evidence),
         "definition": src("R2" if previous else "Ar", "aryl/acyl substituent as drawn", evidence),
         "allowed_values": st.fact(evidence=evidence), "constraints": st.fact(evidence=evidence)}
    ]
    m["generic_template"] = {
        "fixed_structure": src(
            "Figure 1 generic product scaffold" if product else "Figure 1 generic 3-aroylcoumarine scaffold",
            {"scaffold": "dihydrocoumarin-fused dihydropyranone" if product else "3-aroylcoumarine"}, evidence
        ),
        "variables": vars_, "instance_mapping": [],
        "reason": src("R-group drawing", "The source uses R/Ar variables, so no concrete SMILES is asserted.", evidence)
    }
    if product:
        set_product_stereo(m, evidence)
    return m


def set_concrete_simple(m, smiles, evidence, derivation_text):
    diag = st.stereo_inspect(smiles)
    m["representation_type"] = "smiles"
    m["smiles"] = derived(smiles, evidence, derivation_text, raw=m["names"]["original"]["raw_value"] or m["occurrences"][0]["compound_label"]["raw_value"])
    m["stereochemistry"]["chirality_status"] = inferred("achiral", evidence, "Candidate graph has no stereogenic atom or axis relevant to the named material.")
    m["stereochemistry"]["geometry_status"] = inferred("not_applicable", evidence, "No source E/Z or ring cis/trans element applies to this candidate graph.")
    m["stereochemistry"]["representation_status"] = inferred("not_applicable", evidence, "No applicable stereochemistry needs encoding.")
    m["stereochemistry"]["identity_key"] = derived(diag["stereo_graph_key"], evidence, "Generated by standard stereo from the candidate SMILES.")


MOLECULES = {}


def add(m):
    MOLECULES[m["molecule_id"]] = m
    return m


# Core templates and allenoates.
add(make_template("m_generic_1", "generic 3-aroylcoumarine 1", ev(2, 5, 6), product=False))
add(make_template("m_generic_3", "generic product 3", ev(2, 5, 6), product=True))
add(make_template("m_prev_1", "previous-work generic substrate 1", ev(1, 2), product=False, previous=True))
add(make_template("m_prev_3", "previous-work generic product 3", ev(1, 2), product=True, previous=True))
add(make_allenoate("m_2a", "2a", "ethyl 2,3-butadienoate", "C=C=CC(=O)OCC", ev(1, 2, 5),
                   "b296108b7dc16b6a6a6899022a9a40aaffa779c671ceff63703b00d10db18b30"))
add(make_allenoate("m_2b", "2b", "benzyl buta-2,3-dienoate", "C=C=CC(=O)OCc1ccccc1", ev(1, 2, 3, 5, 6),
                   "ee1b789569b20b530310b933abe72ec86000e320a33e7bc27b014124d8de8289"))

# Source-specific substrate/product identities. No product SMILES is invented.
product_notes = {
    "a": "unsubstituted coumarin ring; 4-phenyl product",
    "b": "methoxy-substituted coumarin ring (product name on page 6: 9-methoxy)",
    "c": "methoxy-substituted coumarin ring at a different position",
    "d": "methoxy-substituted coumarin ring at a different position",
    "e": "methyl-substituted coumarin ring",
    "f": "di-tert-butyl-substituted coumarin ring",
    "g": "fluoro-substituted coumarin ring",
    "h": "chloro-substituted coumarin ring",
    "i": "bromo-substituted coumarin ring",
    "j": "bromo-substituted coumarin ring at a different position",
    "k": "p-tolyl aroyl-derived substituent",
    "l": "m-tolyl aroyl-derived substituent",
    "m": "p-chlorophenyl aroyl-derived substituent",
    "n": "m-chlorophenyl aroyl-derived substituent",
    "o": "[1,1'-biphenyl]-4-yl aroyl-derived substituent",
    "p": "thiophene-derived aroyl substituent; exact source naming should be rechecked"
}
for letter in "abcdefghijklmnop":
    substrate_ev = ev(2, 3, 5, 6) if letter == "a" else ev(4, 5, 6)
    sname = "3-benzoylcoumarine" if letter == "a" else None
    sm = make_molecule(f"m_1{letter}", f"1{letter}", substrate_ev, sname, "achiral", "not_applicable", "not_encoded",
                       "Substrate identity follows the source label and scope drawing; no unverified concrete SMILES is asserted. " + product_notes[letter])
    add(sm)
    pname = None
    if letter == "a":
        pname = "Benzyl (E)-2-(5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate"
    elif letter == "b":
        pname = "Benzyl (E)-2-(9-methoxy-5-oxo-4-phenyl-1,10b-dihydro-2H,5H-pyrano[3,4-c]chromen-2-ylidene)acetate"
    pm = make_molecule(f"m_3{letter}", f"3{letter}", ev(2, 4, 5, 6), pname,
                       "partially_specified", "specified", "not_encoded", product_notes[letter])
    set_product_stereo(pm, ev(2, 4, 5))
    add(pm)

# Product 5 shares the assigned E geometry but has no assigned absolute configuration.
m5 = make_molecule("m_5", "5", ev(4, 5), None, "partially_specified", "specified", "not_encoded",
                   "Ethyl-ester product used for E/Z assignment; absolute configuration remains unassigned.")
set_product_stereo(m5, ev(4, 5))
add(m5)

# Catalysts and intermediates.
catalyst_names = {
    "4a": "Quinine", "4b": "Cinchonine", "4c": "(DHQD)2PHAL", "4d": "(DHQD)2AQN",
    "4e": "(DHQ)2PHAL", "4f": "(DHQ)2AQN", "4g": "beta-isocupreidine (R = OH)",
    "4h": "6'-phenyl-beta-iso-cinchonine", "4i": "6'-(2-naphthyl)-beta-iso-cinchonine",
    "4j": "6'-(3,5-di-tert-butylphenyl)-beta-iso-cinchonine",
    "4k": "6'-(3,5-bis(trifluoromethyl)phenyl)-beta-iso-cinchonine",
    "4l": "6'-(4-biphenyl)-beta-iso-cinchonine",
    "4m": "6'-(terphenyl)-beta-iso-cinchonine"
}
for label, name in catalyst_names.items():
    cm = make_molecule(f"m_{label}", label, ev(2, 3, 5, 6), name, "specified", "not_applicable", "not_encoded",
                       "Chiral catalyst structure/name is source-identified, but a concrete stereochemical SMILES was not transcribed.")
    cm["stereochemistry"]["atom_centers"] = src(
        "source catalyst drawing contains wedge/dash stereochemistry",
        [{"configuration": "source-depicted; R/S not transcribed"}], ev(3)
    )
    cm["stereochemistry"]["assignment"] = src(
        "chiral cinchona-derived catalyst as drawn", {"absolute_configuration": "source drawing retained; not encoded"}, ev(3)
    )
    add(cm)

# Other named materials. Simple solvents/reagents remain source-identified without invented CAS data.
named = {
    "m_dabco": ("DABCO", "DABCO"), "m_dce": ("DCE", "1,2-dichloroethane"),
    "m_cyclohexane": ("n-C6H12", "cyclohexane"), "m_toluene": ("toluene", "toluene"),
    "m_thf": ("THF", "tetrahydrofuran"), "m_acetone": ("Acetone", "acetone"),
    "m_etoac": ("EtOAc", "ethyl acetate"), "m_acetonitrile": ("CH3CN", "acetonitrile"),
    "m_dmf": ("DMF", "N,N-dimethylformamide"), "m_meoh": ("MeOH", "methanol"),
    "m_ch2cl2": ("CH2Cl2", "dichloromethane"), "m_pdc": ("Pd/C", "palladium on carbon"),
    "m_h2": ("H2", "hydrogen"), "m_dcc": ("DCC", "dicyclohexylcarbodiimide"),
    "m_dmap": ("DMAP", "4-dimethylaminopyridine"), "m_etoh": ("EtOH", "ethanol"),
    "m_pe": ("PE", "petroleum ether"), "m_intermediate_I": ("I", "proposed zwitterionic intermediate I"),
    "m_intermediate_II": ("II", "proposed intermediate II")
}
for mid, (label, name) in named.items():
    pages = ev(1, 2) if mid in {"m_dabco", "m_intermediate_I", "m_intermediate_II"} else ev(3, 5, 6)
    add(make_molecule(mid, label, pages, name, "unresolved", "unresolved", "unresolved",
                      "Named source material; structure was not independently transcribed for this page-bounded batch."))


def participant(pid, mid, role, evidence, stage="reaction"):
    return st.new_participant(pid, mid, role, evidence, stage)


def put_charge(p, key, raw, normalized, evidence, method="original_extraction", derivation=None):
    p["charges"][key] = st.fact(raw, normalized, evidence, method=method, derivation=derivation)


def put_reference(p, mid, raw, evidence, method="original_extraction", derivation=None):
    p["equivalence_reference"] = st.fact(raw, mid, evidence, method=method, derivation=derivation)


def classification(evidence, extra_esterification=False):
    changes = [
        {"change_id": "c1", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}],
         "before_order": 0, "after_order": 1, "change_type": "formed", "source_basis": "annulation ring closure"},
        {"change_id": "c2", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
         "before_order": 0, "after_order": 1, "change_type": "formed", "source_basis": "pyran ring closure"}
    ]
    if extra_esterification:
        changes.append({"change_id": "c3", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "O"}],
                        "before_order": 0, "after_order": 1, "change_type": "formed", "source_basis": "ethyl ester formation"})
    labels = [
        {"category_number": 2, "code": "C-C", "supporting_change_ids": ["c1"]},
        {"category_number": 3, "code": "C-O", "supporting_change_ids": ["c2"] + (["c3"] if extra_esterification else [])}
    ]
    return {
        "standard_id": st.STANDARD,
        "mapped_reaction_smiles": st.fact(evidence=evidence),
        "bond_changes": [derived(c, evidence, "Manual source-drawing comparison; atom mapping is unresolved.") for c in changes],
        "labels": derived(labels, evidence, "Taxonomy labels follow the manually identified newly formed bonds."),
        "coverage": inferred("partial", evidence, "Source drawings support net C-C/C-O formation, but complete atom mapping and implicit-H accounting were not established."),
        "limitations": [
            "No complete atom-mapped reaction SMILES was asserted.",
            "Bond-order shifts, proton transfers, coproducts, and implicit-hydrogen changes remain unresolved.",
            "Classification is a machine derivation from the supplied drawings, not source verification."
        ],
        "model_votes": []
    }


def base_reaction(rid, evidence, record_type="specific_experiment", categories=None, location=None,
                  name="asymmetric [4 + 2] cyclization", transformation="[4 + 2] annulation"):
    r = st.new_reaction(f"{DOC}:{rid}", DOC, evidence)
    r["record_type"] = record_type
    r["experiment_categories"] = src(categories or ["other"], categories or ["other"], evidence)
    r["location"] = src(location or rid, {"pages": sorted({int(x.split(":page")[1].split(":")[0]) for x in evidence})}, evidence)
    r["reaction_name"] = src(name, name, evidence)
    r["transformation_type"] = inferred(transformation, evidence, "Transformation type normalized from the source reaction title/caption.")
    r["classification"] = classification(evidence)
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
                   "unresolved": ["Machine-extracted page-bounded record; chemical identity and full atom mapping require expert review."]}
    r["completeness"] = "partial"
    return r


def step(step_id, stage, pids, operation, evidence, temperature=None, time=None,
         operation_method="original_extraction", condition_method="original_extraction"):
    s = st.new_step(step_id, evidence, stage)
    s["participant_ids"] = pids
    s["operation"] = st.fact(operation, operation, evidence, method=operation_method,
                              derivation="Inherited from the page-6 general procedure." if operation_method == "inherited_general_condition" else None)
    if temperature is not None:
        raw_temp, norm_temp = temperature
        s["temperature"] = st.fact(raw_temp, norm_temp, evidence, method=condition_method,
                                    derivation="Inherited from the stated standard conditions/general procedure." if condition_method == "inherited_general_condition" else None)
    if time is not None:
        raw_time, norm_time = time
        s["time"] = st.fact(raw_time, norm_time, evidence, method=condition_method,
                            derivation="Inherited from the stated standard conditions/general procedure." if condition_method == "inherited_general_condition" else None)
    return s


def outcome(oid, kind, target, raw_value, normalized, evidence, basis, scope=None, yield_type=None):
    o = st.new_outcome(oid, kind, [target] if target else [], evidence)
    o["product_role"] = src("desired product", "desired product", evidence) if target else st.fact(evidence=evidence)
    o["value"] = src(raw_value, normalized, evidence)
    o["yield_type"] = src(yield_type, yield_type, evidence) if yield_type else st.fact(evidence=evidence)
    o["measurement_basis"] = src(basis, basis, evidence) if basis else st.fact(evidence=evidence)
    o["scope"] = src(scope, scope, evidence) if scope is not None else st.fact(evidence=evidence)
    return o


def conflict_outcome(oid, target, evidence, alternatives, basis, scope):
    o = st.new_outcome(oid, "selectivity", [target], evidence)
    o["product_role"] = src("desired product", "desired product", evidence)
    alt_facts = [src(raw, {"number": value, "unit": "%"}, alt_ev) for raw, value, alt_ev in alternatives]
    o["value"] = st.fact([a[0] for a in alternatives], None, evidence, status="conflict", alternatives=alt_facts)
    o["yield_type"] = st.fact(evidence=evidence)
    o["measurement_basis"] = src(basis, basis, evidence)
    o["scope"] = src(scope, scope, evidence)
    return o


def add_manual_workup(r, evidence, inherited=False):
    method = "inherited_general_condition" if inherited else "original_extraction"
    r["participants"] += [
        participant("p_pe_purification", "m_pe", "purification_agent", evidence, "purification"),
        participant("p_etoac_purification", "m_etoac", "purification_agent", evidence, "purification")
    ]
    r["steps"].append(step("workup_evaporation", "workup", [], "solvent removed by rotary evaporation", evidence,
                           operation_method=method))
    r["steps"].append(step("purification_column", "purification", ["p_pe_purification", "p_etoac_purification"],
                           "direct silica-gel column chromatography; PE/EtOAc 20:1–15:1", evidence,
                           operation_method=method))


REACTIONS = []

# Current generic scheme and the page-6 procedure.
gen_ev = ev(2, 5, 6)
g = base_reaction("scheme_current_general", gen_ev, "generic_scheme", ["substrate_scope", "other"],
                  "Figure 1b, Figure 2 header, and section 3.2 general procedure")
g["participants"] = [
    participant("p_generic_1", "m_generic_1", "reactant", gen_ev),
    participant("p_2b", "m_2b", "reactant", gen_ev),
    participant("p_4l", "m_4l", "catalyst", gen_ev),
    participant("p_etoac", "m_etoac", "solvent", gen_ev),
    participant("p_generic_3", "m_generic_3", "product", gen_ev)
]
put_charge(g["participants"][0], "amount", "1 (0.2 mmol)", {"number": 0.2, "unit": "mmol"}, ev(5, 6))
put_charge(g["participants"][1], "amount", "2b (0.3 mmol)", {"number": 0.3, "unit": "mmol"}, ev(6))
put_charge(g["participants"][1], "equivalents", "1.5 equiv.", {"number": 1.5, "unit": "equiv"}, ev(5, 6))
put_reference(g["participants"][1], "m_generic_1", "relative to 1", ev(5, 6))
put_charge(g["participants"][2], "amount", "4l (0.02 mmol)", {"number": 0.02, "unit": "mmol"}, ev(6))
put_charge(g["participants"][2], "mol_percent", "10 mol%", {"number": 10, "unit": "mol%"}, ev(2, 5, 6))
put_reference(g["participants"][2], "m_generic_1", "relative to 1", ev(2, 5, 6))
put_charge(g["participants"][3], "volume", "EtOAc (1.0 mL)", {"number": 1.0, "unit": "mL"}, ev(6))
g["reaction_scale"] = src("3-aroylcoumarines 1 (0.2 mmol)", {"reference_molecule_id": "m_generic_1", "amount": {"number": 0.2, "unit": "mmol"}}, ev(5, 6))
g["solvent_conditions"]["composition"] = src("EtOAc", ["ethyl acetate"], gen_ev)
g["solvent_conditions"]["total_volume"] = src("1.0 mL", {"number": 1.0, "unit": "mL"}, ev(6))
g["solvent_conditions"]["concentration"] = src("EtOAc (0.2 M)", {"number": 0.2, "unit": "M"}, ev(5))
g["solvent_conditions"]["concentration_basis"] = inferred("relative to substrate 1", ev(5, 6), "0.2 mmol substrate in 1.0 mL EtOAc matches the figure's 0.2 M label.")
g["steps"] = [
    step("pre_stir", "reaction", ["p_generic_1", "p_4l", "p_etoac"],
         "sequentially add 1, 4l, and EtOAc; stir before charging 2b", ev(6),
         ("0 °C", {"number": 0, "unit": "°C"}), ("15 min", {"number": 15, "unit": "min"})),
    step("annulation", "reaction", ["p_generic_1", "p_2b", "p_4l", "p_etoac", "p_generic_3"],
         "charge 2b and monitor by TLC until completion", ev(6),
         ("0 °C", {"number": 0, "unit": "°C"}), None)
]
g["procedure"]["original_text"] = src(
    "To a 4 mL vial was sequentially added 3-aroylcoumarines 1 (0.2 mmol), catalyst 4l (0.02 mmol, 10 mol%), and EtOAc (1.0 mL); the mixture was stirred at 0 °C for 15 min before benzyl buta-2,3-dienoate 2b (0.3 mmol, 1.5 equiv.) was charged.",
    "page-6 general procedure", ev(6))
g["procedure"]["addition_order"] = src("1, then 4l, then EtOAc; pre-stir; then charge 2b", ["1", "4l", "EtOAc", "pre-stir", "2b"], ev(6))
g["procedure"]["pretreatment"] = src("stirred at 0 °C for 15 min before 2b", "0 °C pre-stir for 15 min", ev(6))
add_manual_workup(g, ev(6))
g["review"]["unresolved"] += ["Reaction completion time is substrate-dependent and not specified in the general procedure.", "Atmosphere and pressure are unreported."]
REACTIONS.append(g)

# Previous literature example and its explicitly proposed mechanism.
prev_ev = ev(1, 2)
pr = base_reaction("literature_shi_2012_scheme", prev_ev, "literature_example", ["other"],
                   "Introduction and Figure 1a", "DABCO-catalyzed [4 + 2] cyclization", "literature [4 + 2] annulation")
pr["participants"] = [participant("p_prev_1", "m_prev_1", "reactant", prev_ev), participant("p_2a", "m_2a", "reactant", prev_ev),
                      participant("p_dabco", "m_dabco", "catalyst", prev_ev), participant("p_dce", "m_dce", "solvent", prev_ev),
                      participant("p_prev_3", "m_prev_3", "product", prev_ev)]
put_charge(pr["participants"][2], "mol_percent", "20 mol%", {"number": 20, "unit": "mol%"}, ev(2))
put_reference(pr["participants"][2], "m_prev_1", "relative to 1", ev(2))
pr["solvent_conditions"]["composition"] = src("CH2Cl2", ["dichloromethane"], ev(2))
pr["steps"] = [step("literature_annulation", "reaction", [p["participant_id"] for p in pr["participants"]],
                    "DABCO-catalyzed annulation", ev(2), ("rt", "rt"), None)]
pr["outcomes"] = [outcome("yield_range", "yield", "m_prev_3", "79–95%", {"minimum": 79, "maximum": 95, "unit": "%", "raw_text": "79–95%"},
                                  ev(1, 2), None, {"reported_examples": 12}, "reported literature yield range")]
pr["review"]["unresolved"] += ["Cited prior work is retained as literature context, not counted as a current-paper experiment.", "Individual 12 examples are not enumerated in this paper."]
REACTIONS.append(pr)

mech = base_reaction("proposed_mechanism_figure1a", prev_ev, "proposed_mechanism", ["other"],
                     "Introduction prose and Figure 1a possible mechanism", "possible mechanism for DABCO-catalyzed [4 + 2] cyclization",
                     "proposed zwitterionic addition/cyclization mechanism")
mech["participants"] = [participant("p_2a", "m_2a", "reactant", prev_ev), participant("p_prev_1", "m_prev_1", "reactant", prev_ev),
                        participant("p_dabco", "m_dabco", "catalyst", prev_ev), participant("p_I", "m_intermediate_I", "unknown", prev_ev),
                        participant("p_II", "m_intermediate_II", "unknown", prev_ev), participant("p_prev_3", "m_prev_3", "product", prev_ev)]
mech["steps"] = [step("nucleophilic_addition", "reaction", ["p_2a", "p_dabco", "p_I"], "tertiary amine attacks 2a to form zwitterion I", prev_ev),
                 step("conjugate_addition", "reaction", ["p_I", "p_prev_1", "p_II"], "gamma-carbanion attacks the enone beta-carbon to give Z intermediate II", prev_ev),
                 step("ring_closure", "reaction", ["p_II", "p_prev_3", "p_dabco"], "intramolecular nucleophilic substitution gives 3 and regenerates tertiary amine", prev_ev)]
mech["notes"] = [src("Possible mechanism", "Author-proposed mechanism; not an independently verified mechanistic experiment.", ev(2))]
mech["review"]["unresolved"] += ["Arrow-pushing proposal is not treated as experimental proof."]
REACTIONS.append(mech)

# Optimization table: every row is a distinct experiment. Entry 22 is also the repeated 3a scope example.
solvent_mid = {"DCE": "m_dce", "n-C6H12": "m_cyclohexane", "toluene": "m_toluene", "THF": "m_thf",
               "Acetone": "m_acetone", "EtOAc": "m_etoac", "CH3CN": "m_acetonitrile", "DMF": "m_dmf", "MeOH": "m_meoh"}
scope_general_id = f"{DOC}:scheme_current_general"
for row in RAW["table1"]:
    n = row["entry"]
    evidence = (ev(2) if n <= 13 else []) + ev(3) + (ev(5, 6) if n == 22 else [])
    cats = ["optimization"] + (["substrate_scope"] if n == 22 else [])
    r = base_reaction(f"table1_entry{n:02d}", evidence, "specific_experiment", cats, f"Table 1 entry {n}")
    r["general_scheme_ids"] = [scope_general_id]
    r["participants"] = [participant("p_1a", "m_1a", "reactant", evidence), participant("p_2b", "m_2b", "reactant", evidence),
                         participant("p_cat", f"m_{row['catalyst']}", "catalyst", evidence),
                         participant("p_solvent", solvent_mid[row["solvent"]], "solvent", evidence),
                         participant("p_3a", "m_3a", "product", evidence)]
    put_charge(r["participants"][2], "mol_percent", "Cat. (10 mol%)", {"number": 10, "unit": "mol%"}, ev(3))
    put_reference(r["participants"][2], "m_1a", "relative to 1a", ev(3))
    if n == 22:
        put_charge(r["participants"][0], "amount", "1 (0.2 mmol)", {"number": 0.2, "unit": "mmol"}, ev(5, 6), "inherited_general_condition",
                   "Entry 22 is the same 3a standard-condition example repeated in Figure 2/general procedure.")
        put_charge(r["participants"][1], "amount", "2b (0.3 mmol)", {"number": 0.3, "unit": "mmol"}, ev(6), "inherited_general_condition",
                   "Inherited from the page-6 general procedure for the repeated standard-condition example.")
        put_charge(r["participants"][1], "equivalents", "2b (1.5 equiv.)", {"number": 1.5, "unit": "equiv"}, ev(5, 6), "inherited_general_condition",
                   "Inherited from Figure 2/general procedure.")
        put_reference(r["participants"][1], "m_1a", "relative to 1", ev(5, 6), "inherited_general_condition", "Figure 2 reference substrate.")
        put_charge(r["participants"][2], "amount", "4l (0.02 mmol)", {"number": 0.02, "unit": "mmol"}, ev(6), "inherited_general_condition",
                   "Inherited from the page-6 general procedure.")
        put_charge(r["participants"][3], "volume", "EtOAc (1.0 mL)", {"number": 1.0, "unit": "mL"}, ev(6), "inherited_general_condition",
                   "Inherited from the page-6 general procedure.")
        r["reaction_scale"] = st.fact("1 (0.2 mmol)", {"reference_molecule_id": "m_1a", "amount": {"number": 0.2, "unit": "mmol"}}, ev(5, 6),
                                      method="inherited_general_condition", derivation="Repeated 3a scope/standard-condition experiment.")
    r["solvent_conditions"]["composition"] = src(row["solvent"], [row["solvent"]], ev(3))
    if n == 22:
        r["solvent_conditions"]["total_volume"] = st.fact("1.0 mL", {"number": 1.0, "unit": "mL"}, ev(6), method="inherited_general_condition",
                                                          derivation="Inherited from the general procedure.")
        r["solvent_conditions"]["concentration"] = st.fact("EtOAc (0.2 M)", {"number": 0.2, "unit": "M"}, ev(5), method="inherited_general_condition",
                                                            derivation="Figure 2 standard-condition concentration.")
        r["solvent_conditions"]["concentration_basis"] = inferred("relative to 1a", ev(5, 6), "0.2 mmol in 1.0 mL.")
        r["steps"] = [
            step("pre_stir", "reaction", ["p_1a", "p_cat", "p_solvent"], "pre-stir before charging 2b", ev(6),
                 ("0 °C", {"number": 0, "unit": "°C"}), ("15 min", {"number": 15, "unit": "min"}),
                 "inherited_general_condition", "inherited_general_condition"),
            step("annulation", "reaction", [p["participant_id"] for p in r["participants"]], "[4 + 2] cyclization", evidence,
                 ("0 °C", {"number": 0, "unit": "°C"}), ("1 h", {"number": 1, "unit": "h"}))
        ]
        add_manual_workup(r, ev(6), inherited=True)
    else:
        r["steps"] = [step("annulation", "reaction", [p["participant_id"] for p in r["participants"]], "Table 1 model reaction", ev(3),
                           (f"{row['temperature_c']} °C", {"number": row["temperature_c"], "unit": "°C"}),
                           (f"{row['time_h']} h", {"number": row["time_h"], "unit": "h"}))]
    if row["yield"] == "trace":
        r["outcomes"] = [outcome("trace_product", "trace", "m_3a", "trace", "trace", ev(3), "pure isolated substances (global statement on page 6)",
                                  {"metric": "product formation"})]
    else:
        r["outcomes"] = [outcome("isolated_yield", "yield", "m_3a", f"{row['yield']}%", {"number": row["yield"], "unit": "%"}, ev(3),
                                  "pure isolated substances (global statement on page 6)", {"product": "3a"}, "isolated")]
    if n == 18:
        r["outcomes"].append(conflict_outcome("ee", "m_3a", ev(3),
            [("Table 1 entry 18: 79% ee", 79, ev(3)), ("prose below Table 1: 80% ee", 80, ev(3))],
            "chiral HPLC", {"metric": "ee", "major_enantiomer": "unassigned"}))
        r["review"]["unresolved"].append("Entry 18 has 79% ee in the table but 80% ee in the adjacent prose; neither value was selected.")
    elif n == 22:
        r["outcomes"].append(conflict_outcome("ee", "m_3a", ev(3, 5),
            [("Table 1 entry 22: 79% ee", 79, ev(3)), ("Figure 2 product 3a: 80% ee", 80, ev(5))],
            "chiral HPLC; 3a method on page 6", {"metric": "ee", "major_enantiomer": "unassigned"}))
        r["review"]["unresolved"].append("The repeated standard-condition 3a experiment has 79% ee in Table 1 and 80% ee in Figure 2; neither value was selected.")
    elif row["ee"] is not None:
        r["outcomes"].append(outcome("ee", "selectivity", "m_3a", str(row["ee"]) if isinstance(row["ee"], str) else f"{row['ee']}% ee",
                                      row["ee"] if isinstance(row["ee"], str) else {"number": row["ee"], "unit": "%"}, ev(3),
                                      "chiral HPLC; 3a method on page 6", {"metric": "ee", "major_enantiomer": "unassigned"}))
    r["procedure"]["general_procedure_reference"] = src("Table 1 model reaction / section 3.2", "section 3.2" if n == 22 else "Table 1", evidence)
    r["review"]["unresolved"] += ["Optimization-table substrate/allenoate charges and solvent volume are not printed for this row." if n != 22 else "Atmosphere and pressure are unreported."]
    REACTIONS.append(r)

# Scope examples 3b--3p. 3a is already represented by Table 1 entry 22.
for row in RAW["scope"][1:]:
    letter = row["product"][-1]
    evidence = ev(4, 5, 6)
    r = base_reaction(f"scope_{row['product']}", evidence, "specific_experiment", ["substrate_scope"], f"Figure 2 product {row['product']}")
    r["general_scheme_ids"] = [scope_general_id]
    r["participants"] = [participant(f"p_1{letter}", f"m_1{letter}", "reactant", evidence), participant("p_2b", "m_2b", "reactant", evidence),
                         participant("p_4l", "m_4l", "catalyst", evidence), participant("p_etoac", "m_etoac", "solvent", evidence),
                         participant(f"p_3{letter}", f"m_3{letter}", "product", evidence)]
    put_charge(r["participants"][0], "amount", "1 (0.2 mmol)", {"number": 0.2, "unit": "mmol"}, ev(5, 6), "inherited_general_condition",
               "Figure 2 header and page-6 general procedure apply to the scope series.")
    put_charge(r["participants"][1], "amount", "2b (0.3 mmol)", {"number": 0.3, "unit": "mmol"}, ev(6), "inherited_general_condition",
               "Page-6 general procedure applies to the scope series.")
    put_charge(r["participants"][1], "equivalents", "2b (1.5 equiv.)", {"number": 1.5, "unit": "equiv"}, ev(5, 6), "inherited_general_condition",
               "Figure 2 header and page-6 general procedure.")
    put_reference(r["participants"][1], f"m_1{letter}", "relative to substrate 1", ev(5, 6), "inherited_general_condition", "Scope-series denominator.")
    put_charge(r["participants"][2], "amount", "4l (0.02 mmol)", {"number": 0.02, "unit": "mmol"}, ev(6), "inherited_general_condition",
               "Page-6 general procedure.")
    put_charge(r["participants"][2], "mol_percent", "4l (10 mol%)", {"number": 10, "unit": "mol%"}, ev(5, 6), "inherited_general_condition",
               "Figure 2 header and page-6 general procedure.")
    put_reference(r["participants"][2], f"m_1{letter}", "relative to substrate 1", ev(5, 6), "inherited_general_condition", "Scope-series denominator.")
    put_charge(r["participants"][3], "volume", "EtOAc (1.0 mL)", {"number": 1.0, "unit": "mL"}, ev(6), "inherited_general_condition",
               "Page-6 general procedure.")
    r["reaction_scale"] = st.fact("1 (0.2 mmol)", {"reference_molecule_id": f"m_1{letter}", "amount": {"number": 0.2, "unit": "mmol"}}, ev(5, 6),
                                  method="inherited_general_condition", derivation="Figure 2/general procedure scope scale.")
    r["solvent_conditions"]["composition"] = st.fact("EtOAc", ["ethyl acetate"], ev(5, 6), method="inherited_general_condition", derivation="Scope standard conditions.")
    r["solvent_conditions"]["total_volume"] = st.fact("1.0 mL", {"number": 1.0, "unit": "mL"}, ev(6), method="inherited_general_condition", derivation="General procedure.")
    r["solvent_conditions"]["concentration"] = st.fact("EtOAc (0.2 M)", {"number": 0.2, "unit": "M"}, ev(5), method="inherited_general_condition", derivation="Figure 2 header.")
    r["solvent_conditions"]["concentration_basis"] = inferred(f"relative to 1{letter}", ev(5, 6), "0.2 mmol substrate in 1.0 mL EtOAc.")
    r["steps"] = [
        step("pre_stir", "reaction", [f"p_1{letter}", "p_4l", "p_etoac"], "pre-stir 1, 4l, and EtOAc before charging 2b", ev(6),
             ("0 °C", {"number": 0, "unit": "°C"}), ("15 min", {"number": 15, "unit": "min"}),
             "inherited_general_condition", "inherited_general_condition"),
        step("annulation", "reaction", [p["participant_id"] for p in r["participants"]], "charge 2b; monitor by TLC until completion", evidence,
             ("0 °C", {"number": 0, "unit": "°C"}), (f"{row['time_h']} h", {"number": row["time_h"], "unit": "h"}),
             "inherited_general_condition", "original_extraction")
    ]
    add_manual_workup(r, ev(6), inherited=True)
    r["procedure"]["general_procedure_reference"] = src("section 3.2 general procedure", "section 3.2", ev(6))
    r["outcomes"] = [
        outcome("isolated_yield", "yield", f"m_3{letter}", f"{row['yield']}%", {"number": row["yield"], "unit": "%"}, ev(5, 6),
                "pure isolated substances", {"product": row["product"]}, "isolated"),
        outcome("ee", "selectivity", f"m_3{letter}", f"{row['ee']}% ee", {"number": row["ee"], "unit": "%"}, ev(5, 6),
                "chiral HPLC", {"metric": "ee", "major_enantiomer": "unassigned"})
    ]
    if row["product"] in {"3b", "3c"}:
        r["notes"].append(src("page-4 prose says 3b–e have 80–88% ee, while Figure 2 gives 3b 76% ee and 3c 79% ee",
                              "aggregate prose range conflicts with the product-specific Figure 2 value", ev(4, 5)))
        r["review"]["unresolved"].append("The aggregate page-4 80–88% ee statement does not include the product-specific 76%/79% values; Figure 2 value is retained and the discrepancy is not silently corrected.")
    r["review"]["unresolved"].append("Absolute configuration is unassigned; ee does not determine R/S.")
    REACTIONS.append(r)

# Figure 3 direct ethyl-allenoate route.
fig3_ev = ev(4, 5)
direct = base_reaction("figure3_direct_to_5", fig3_ev, "specific_experiment", ["mechanistic_study"], "Figure 3 upper route",
                       "standard-condition [4 + 2] cyclization for E/Z assignment")
direct["general_scheme_ids"] = [scope_general_id]
direct["participants"] = [participant("p_1m", "m_1m", "reactant", fig3_ev), participant("p_2a", "m_2a", "reactant", fig3_ev),
                          participant("p_4l", "m_4l", "catalyst", fig3_ev), participant("p_etoac", "m_etoac", "solvent", fig3_ev),
                          participant("p_5", "m_5", "product", fig3_ev)]
put_charge(direct["participants"][1], "equivalents", "2a (1.5 equiv.)", {"number": 1.5, "unit": "equiv"}, ev(5))
put_reference(direct["participants"][1], "m_1m", "relative to 1m", ev(5))
put_charge(direct["participants"][2], "mol_percent", "standard conditions: 4l (10 mol%)", {"number": 10, "unit": "mol%"}, ev(5),
           "inherited_general_condition", "Standard conditions are defined in Figure 2.")
put_reference(direct["participants"][2], "m_1m", "relative to 1m", ev(5), "inherited_general_condition", "Standard-condition denominator.")
direct["solvent_conditions"]["composition"] = src("EtOAc", ["ethyl acetate"], ev(5))
direct["steps"] = [step("direct_annulation", "reaction", [p["participant_id"] for p in direct["participants"]], "standard conditions", fig3_ev,
                        ("0 °C", {"number": 0, "unit": "°C"}), None)]
direct["outcomes"] = [outcome("isolated_yield", "yield", "m_5", "86%", {"number": 86, "unit": "%"}, ev(5, 6), "pure isolated substances", {"product": "5"}, "isolated"),
                      outcome("ee", "selectivity", "m_5", "45% ee", {"number": 45, "unit": "%"}, ev(4, 5, 6), "chiral HPLC", {"metric": "ee", "major_enantiomer": "unassigned"})]
direct["review"]["unresolved"] += ["Reaction time and exact scale are not printed in Figure 3.", "Absolute configuration remains unassigned."]
REACTIONS.append(direct)

# Figure 3 three-step route. The total yield is kept across the full 1m-to-5 sequence.
route = base_reaction("figure3_three_step_to_5", fig3_ev, "specific_experiment", ["mechanistic_study", "derivatization"], "Figure 3 lower route",
                      "three-step conversion used for E/Z assignment", "annulation, hydrogenolysis, and esterification")
route["general_scheme_ids"] = [scope_general_id]
route["classification"] = classification(fig3_ev, extra_esterification=True)
route["participants"] = [
    participant("p_1m", "m_1m", "reactant", fig3_ev), participant("p_2b", "m_2b", "reactant", fig3_ev),
    participant("p_4l", "m_4l", "catalyst", fig3_ev), participant("p_etoac", "m_etoac", "solvent", fig3_ev),
    participant("p_3m", "m_3m", "product", fig3_ev), participant("p_pdc", "m_pdc", "catalyst", fig3_ev),
    participant("p_h2", "m_h2", "reagent", fig3_ev), participant("p_dcc", "m_dcc", "reagent", fig3_ev),
    participant("p_dmap", "m_dmap", "catalyst", fig3_ev), participant("p_etoh", "m_etoh", "reagent", fig3_ev),
    participant("p_ch2cl2", "m_ch2cl2", "solvent", fig3_ev), participant("p_5", "m_5", "product", fig3_ev)
]
put_charge(route["participants"][1], "equivalents", "standard conditions use 2b (1.5 equiv.)", {"number": 1.5, "unit": "equiv"}, ev(5),
           "inherited_general_condition", "Figure 2 standard conditions.")
put_reference(route["participants"][1], "m_1m", "relative to 1m", ev(5), "inherited_general_condition", "Standard-condition denominator.")
put_charge(route["participants"][2], "mol_percent", "standard conditions: 4l (10 mol%)", {"number": 10, "unit": "mol%"}, ev(5),
           "inherited_general_condition", "Figure 2 standard conditions.")
put_reference(route["participants"][2], "m_1m", "relative to 1m", ev(5), "inherited_general_condition", "Standard-condition denominator.")
route["steps"] = [
    step("annulation", "reaction", ["p_1m", "p_2b", "p_4l", "p_etoac", "p_3m"], "standard conditions form 3m", fig3_ev,
         ("0 °C", {"number": 0, "unit": "°C"}), None),
    step("hydrogenolysis", "reaction", ["p_3m", "p_pdc", "p_h2"], "then Pd/C, H2", fig3_ev),
    step("esterification", "reaction", ["p_dcc", "p_dmap", "p_etoh", "p_ch2cl2", "p_5"], "DCC, DMAP; EtOH/CH2Cl2", fig3_ev)
]
route["outcomes"] = [
    outcome("total_yield", "yield", "m_5", "48% total yield", {"number": 48, "unit": "%"}, ev(5, 6), "pure isolated substances",
            {"start": "1m", "end": "5", "span": "three-step sequence"}, "three-step total isolated yield"),
    outcome("ee", "selectivity", "m_5", "79% ee", {"number": 79, "unit": "%"}, ev(4, 5, 6), "chiral HPLC",
            {"metric": "ee", "major_enantiomer": "unassigned"})
]
route["notes"] = [src("page-4 prose says a 3-step sequence from 3m and 2b; Figure 3 lower route starts from a structure labeled 1m under standard conditions",
                      "source label/prose ambiguity retained; sequence modeled as 1m + 2b -> 3m -> 5", fig3_ev)]
route["review"]["unresolved"] += ["The page-4 phrase 'from 3m and 2b' conflicts with the lower scheme's starting label 1m; no authorial correction is asserted.",
                                   "Step-specific quantities, temperatures, and times after annulation are not printed on pages 1–6.",
                                   "The 48% value is a three-step total yield, not a last-step yield."]
REACTIONS.append(route)

# Figure 4 gram-scale experiment.
scale_ev = ev(4, 6)
scale = base_reaction("figure4_gram_scale_3f", scale_ev, "specific_experiment", ["scale_up"], "Figure 4 gram-scale synthesis")
scale["general_scheme_ids"] = [scope_general_id]
scale["participants"] = [participant("p_1f", "m_1f", "reactant", scale_ev), participant("p_2b", "m_2b", "reactant", scale_ev),
                         participant("p_4l", "m_4l", "catalyst", scale_ev), participant("p_etoac", "m_etoac", "solvent", scale_ev),
                         participant("p_3f", "m_3f", "product", scale_ev)]
put_charge(scale["participants"][0], "amount", "1f (2.5 mmol)", {"number": 2.5, "unit": "mmol"}, ev(4, 6))
put_charge(scale["participants"][1], "amount", "2b (3.75 mmol)", {"number": 3.75, "unit": "mmol"}, ev(6))
put_charge(scale["participants"][1], "equivalents", "1.5 equivalents of 2b", {"number": 1.5, "unit": "equiv"}, ev(4, 6))
put_reference(scale["participants"][1], "m_1f", "relative to 1f", ev(4, 6))
put_charge(scale["participants"][2], "mol_percent", "4l (2.5 mol%)", {"number": 2.5, "unit": "mol%"}, ev(4, 6))
put_reference(scale["participants"][2], "m_1f", "relative to 1f", ev(4, 6))
scale["reaction_scale"] = src("1f (2.5 mmol)", {"reference_molecule_id": "m_1f", "amount": {"number": 2.5, "unit": "mmol"}}, scale_ev)
scale["solvent_conditions"]["composition"] = src("EtOAc", ["ethyl acetate"], ev(6))
scale["solvent_conditions"]["concentration"] = src("EtOAc (0.4 M)", {"number": 0.4, "unit": "M"}, ev(6))
scale["solvent_conditions"]["concentration_basis"] = inferred("relative to 1f", ev(6), "Figure 4 labels the reaction EtOAc (0.4 M).")
scale["steps"] = [step("gram_scale_annulation", "reaction", [p["participant_id"] for p in scale["participants"]], "gram-scale [4 + 2] cyclization", scale_ev,
                       ("0 °C", {"number": 0, "unit": "°C"}), ("10 h", {"number": 10, "unit": "h"}))]
scale["outcomes"] = [
    outcome("isolated_yield", "yield", "m_3f", "91%", {"number": 91, "unit": "%"}, scale_ev, "pure isolated substances", {"product": "3f"}, "isolated"),
    outcome("product_mass", "product_amount", "m_3f", "1.215 g", {"number": 1.215, "unit": "g"}, ev(4, 6), "weighed isolated product", {"product": "3f"}),
    outcome("ee", "selectivity", "m_3f", "93% ee", {"number": 93, "unit": "%"}, scale_ev, "chiral HPLC", {"metric": "ee", "major_enantiomer": "unassigned"})
]
scale["review"]["unresolved"] += ["Exact EtOAc volume, catalyst amount, workup, and purification details for the scaled experiment are not printed on pages 1–6.",
                                   "Absolute configuration remains unassigned."]
REACTIONS.append(scale)

# Tighten field-level provenance: result values cite their result figure/table,
# while measurement method and isolated-yield type cite the page-6 global method.
for reaction in REACTIONS:
    rid = reaction["reaction_id"]
    for result in reaction["outcomes"]:
        if result["measurement_basis"]["normalized_value"] is not None:
            result["measurement_basis"]["evidence_ids"] = ev(6)
        if result["yield_type"]["normalized_value"] is not None:
            result["yield_type"]["evidence_ids"] = ev(6)
        if ":table1_entry" in rid:
            result["value"]["evidence_ids"] = ev(3, 5) if rid.endswith("entry22") else ev(3)
        elif ":scope_3" in rid:
            result["value"]["evidence_ids"] = ev(5)
        elif ":figure3_" in rid:
            result["value"]["evidence_ids"] = ev(4, 5)

BUNDLE["molecules"] = list(MOLECULES.values())
BUNDLE["reactions"] = REACTIONS
BUNDLE["source_context"] = [
    {"context_id": f"{DOC}:title", "kind": "bibliographic", "evidence_ids": ev(1), "title": RAW["title"], "doi": "10.3390/molecules26020489"},
    {"context_id": f"{DOC}:page1_prior_work", "kind": "literature_context", "evidence_ids": ev(1),
     "note": "Introduction surveys cited prior [4+2] methods; citations other than the explicitly diagrammed Shi example were not imported as current-paper experiments."},
    {"context_id": f"{DOC}:configuration_assignment", "kind": "stereochemistry", "evidence_ids": ev(4, 5),
     "note": "Products 3 were assigned E by conversion to known compound 5 and NMR comparison. Absolute configuration by recrystallization/X-ray was unsuccessful."},
    {"context_id": f"{DOC}:global_methods", "kind": "characterization_and_method", "evidence_ids": ev(6),
     "note": "Chemical yields refer to pure isolated substances; chiral HPLC was used for enantiomeric purity. Independent NMR/HRMS peak lists are characterization, not separate synthesis experiments."},
    {"context_id": f"{DOC}:page4_scope_discussion", "kind": "source_discrepancy", "evidence_ids": ev(4, 5),
     "note": "Page-4 aggregate 80–88% ee statement for 3b–e conflicts with Figure 2 values for 3b (76%) and 3c (79%); page-4 para/meta product labels also appear inconsistent with Figure 2 structures. No silent correction was made."},
    {"context_id": f"{DOC}:si_status", "kind": "source_limitation", "evidence_ids": ev(6),
     "note": "The main paper states supplementary spectra/HPLC are available online, but no SI was explicitly supplied or associated in this workspace."}
]

(RUN / "submitted.json").write_text(json.dumps(BUNDLE, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(BUNDLE["molecules"]), "reactions": len(BUNDLE["reactions"]),
                  "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in BUNDLE["reactions"])}, ensure_ascii=False))
