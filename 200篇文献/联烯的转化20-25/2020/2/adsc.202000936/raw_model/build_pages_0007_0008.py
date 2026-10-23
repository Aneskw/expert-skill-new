#!/usr/bin/env python3
"""Augment adsc.202000936 records with visually checked PDF pages 7-8.

Page 7 supplies the general procedure for the Table 4 scope experiments.
It is merged into those existing experiments rather than emitted as a duplicate
reaction. Page 8 is references only and contributes source context/coverage.
All formal records remain machine records needing expert review.
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

DOC = "doc_2abd3d7449bbf0b6"
RUN = ROOT / "runs/current"
P5 = [f"{DOC}:page5:image", f"{DOC}:page5:text"]
P7 = [f"{DOC}:page7:image", f"{DOC}:page7:text"]
P8 = [f"{DOC}:page8:image", f"{DOC}:page8:text"]


def fact(raw=None, normalized=None, ev=(), method="original_extraction",
         derivation=None, status=None, alternatives=()):
    return st.fact(raw, normalized, ev, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def missing(ev=()):
    return st.fact(evidence=ev)


def not_applicable(ev=()):
    return st.fact(evidence=ev, method="not_applicable", status="not_applicable")


# Start from the separately preserved pages 1-6 raw candidate so rerunning this
# batch cannot stack duplicate page-7 participants or steps.
bundle = json.loads((RUN / "raw_model/extracted_bundle_pages_0001_0006.json").read_text())
dataset = json.loads((RUN / "dataset.json").read_text())
bundle["documents"] = dataset["documents"]
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}


def add_material(mid: str, label: str, name: str, *, mixture=None,
                 smiles=None, chirality="not_applicable",
                 geometry="not_applicable", representation="not_applicable"):
    if mid in molecules:
        return molecules[mid]
    m = st.new_molecule(mid, DOC, P7)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = fact(label, label, P7)
    m["names"]["original"] = fact(name, name, P7)
    m["names"]["english"] = fact(name, name, P7)
    m["names"]["abbreviations"] = missing(P7)
    m["names"]["aliases"] = missing(P7)
    m["cas_numbers"] = missing(P7)
    m["external_properties"] = []
    m["generic_template"] = None
    m["stereochemistry"]["axial"] = not_applicable(P7)
    m["stereochemistry"]["assignment"] = not_applicable(P7)
    m["stereochemistry"]["chirality_status"] = fact(
        {"source_label": label}, chirality, P7, method="inference",
        derivation="The named workup or purification material has no applicable molecular stereochemical assignment; this is not a claim about every component of a formulation.")
    m["stereochemistry"]["geometry_status"] = fact(
        {"source_label": label}, geometry, P7, method="inference",
        derivation="No E/Z or ring-relative geometry applies to the source-named workup or purification material as registered here.")
    m["stereochemistry"]["representation_status"] = fact(
        {"source_label": label}, representation, P7, method="inference",
        derivation="The source names a formulation/material rather than a unique stereochemical molecular graph." if not smiles else
                   "The candidate molecular graph contains no stereochemical element requiring encoding.")
    if mixture is None:
        m["stereochemistry"]["mixture"] = not_applicable(P7)
    else:
        m["stereochemistry"]["mixture"] = fact(label, mixture, P7)
    if smiles:
        diag = st.stereo_inspect(smiles)
        m["representation_type"] = "smiles"
        m["smiles"] = fact(
            {"source_formula": label}, smiles, P7, method="structure_derivation",
            derivation="Candidate formula-level SMILES derived from the printed Na2SO4 label; parsing and stereo diagnostics do not establish hydration state or source verification.")
        m["stereochemistry"]["atom_centers"] = fact(
            {"candidate_smiles": smiles}, diag["atom_centers"], P7,
            method="structure_derivation", derivation="Standard stereo diagnostics for the candidate drying-agent graph.")
        m["stereochemistry"]["double_bonds"] = fact(
            {"candidate_smiles": smiles}, diag["double_bonds"], P7,
            method="structure_derivation", derivation="Standard stereo diagnostics for the candidate drying-agent graph.")
        m["stereochemistry"]["isomeric_smiles"] = fact(
            {"candidate_smiles": smiles}, diag["canonical_isomeric_smiles"], P7,
            method="structure_derivation", derivation="Canonical isomeric SMILES from standard stereo diagnostics.")
        m["stereochemistry"]["identity_key"] = fact(
            {"candidate_smiles": smiles}, diag["stereo_graph_key"], P7,
            method="structure_derivation", derivation="Stereo graph key for comparing the candidate graph; not a sample identity.")
    else:
        m["representation_type"] = "unresolved"
        m["smiles"] = missing(P7)
        m["stereochemistry"]["atom_centers"] = not_applicable(P7)
        m["stereochemistry"]["double_bonds"] = not_applicable(P7)
        m["stereochemistry"]["isomeric_smiles"] = missing(P7)
        m["stereochemistry"]["identity_key"] = missing(P7)
    m["identity_note"] = fact(
        {"source_label": label},
        "Source-named workup/purification material from page 7; exact formulation or physical state is preserved as unresolved where the source does not define it.",
        P7)
    molecules[mid] = m
    bundle["molecules"].append(m)
    return m


add_material(
    "m_saturated_nh4cl_solution", "saturated NH4Cl solution",
    "saturated ammonium chloride solution",
    mixture={"kind": "aqueous formulation", "components": ["ammonium chloride", "water"], "concentration": "saturated; numeric concentration unreported"})
add_material(
    "m_brine", "brine", "brine",
    mixture={"kind": "aqueous formulation", "components": ["sodium chloride", "water"], "composition": "unreported"})
add_material(
    "m_sodium_sulfate", "Na2SO4", "sodium sulfate drying agent",
    smiles="[Na+].[Na+].[O-]S(=O)(=O)[O-]", chirality="achiral",
    representation="encoded")
add_material(
    "m_silica_gel", "silica gel", "silica gel",
    mixture={"kind": "heterogeneous chromatographic stationary phase", "composition": "not molecularly specified"})
add_material(
    "m_petroleum_ether", "petroleum ether", "petroleum ether",
    mixture={"kind": "hydrocarbon mixture", "composition": "unreported"})


def participant_by_id(r, participant_id):
    return next(p for p in r["participants"] if p["participant_id"] == participant_id)


def add_participant(r, pid, mid, role, stage, *, volume=None, ev=P7,
                    quantity_method="inherited_general_condition"):
    p = st.new_participant(pid, mid, role, ev, stage=stage)
    p["product_role"] = not_applicable(ev)
    if volume is not None:
        raw, normalized = volume
        p["charges"]["volume"] = fact(raw, normalized, ev, method=quantity_method,
            derivation="Inherited from the printed general procedure for the ordinary Table 4 scope experiment." if quantity_method == "inherited_general_condition" else None)
    r["participants"].append(p)
    return p


def volume(raw, number):
    return raw, {"number": number, "unit": "mL"}


def set_amount(p, raw, number, *, ev=P7, method="inherited_general_condition", derivation=None):
    p["charges"]["amount"] = fact(raw, {"number": number, "unit": "mmol"}, ev,
        method=method, derivation=derivation or "Inherited from the printed general procedure for the ordinary Table 4 scope experiment.")


def append_evidence(f, ev):
    f["evidence_ids"] = list(dict.fromkeys(f["evidence_ids"] + list(ev)))


def conflict_scale(r):
    substrate = participant_by_id(r, "substrate")
    table_amount = fact("0.5 mmol scale", {"number": 0.5, "unit": "mmol"}, P5)
    prose_amount = fact("aldehyde (0.47 mmol)", {"number": 0.47, "unit": "mmol"}, P7)
    substrate["charges"]["amount"] = fact(
        {"Table 4 footnote": "0.5 mmol scale", "General Procedure": "aldehyde (0.47 mmol)"},
        None, P5 + P7, status="conflict", alternatives=[table_amount, prose_amount])
    r["reaction_scale"] = fact(
        {"Table 4 footnote": "0.5 mmol scale", "General Procedure": "aldehyde (0.47 mmol)"},
        None, P5 + P7, status="conflict",
        alternatives=[
            fact("0.5 mmol scale", {"reference_molecule_id": substrate["molecule_id"], "amount": {"number": 0.5, "unit": "mmol"}}, P5),
            fact("aldehyde (0.47 mmol)", {"reference_molecule_id": substrate["molecule_id"], "amount": {"number": 0.47, "unit": "mmol"}}, P7),
        ])
    r["notes"].append(fact(
        {"Table 4 footnote": "0.5 mmol scale", "General Procedure": "aldehyde (0.47 mmol)"},
        "The table-level nominal scale and the detailed general-procedure aldehyde charge disagree; both source values are retained without selecting one.",
        P5 + P7, method="inference",
        derivation="Direct comparison of the visually inspected Table 4 footnote and General Procedure."))
    msg = "Source conflict retained: Table 4 states a 0.5 mmol scale, whereas the General Procedure charges 0.47 mmol aldehyde."
    if msg not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append(msg)


GENERAL_TEXT = (
    "Allenyltrichlorosilane (0.71 mmol, 1.5 equiv) in 0.5 mL DCM was added to catalyst "
    "(0.0047 mmol, 1 mol%), DIPEA (0.71 mmol, 1.5 equiv), and aldehyde (0.47 mmol) in "
    "1.5 mL DCM under argon at -60 °C; held 8 h; quenched with saturated NH4Cl (2 mL); "
    "extracted with Et2O (3 × 20 mL); washed with brine; dried over Na2SO4; concentrated; "
    "and flash chromatographed on silica gel with petroleum ether/ethyl acetate 9:1."
)


def step(step_id, stage, operation_raw, operation_norm, participant_ids, *,
         temperature=None, time=None, atmosphere=None, ev=P7,
         method="inherited_general_condition", derivation=None):
    s = st.new_step(step_id, ev, stage=stage)
    s["participant_ids"] = participant_ids
    s["operation"] = fact(operation_raw, operation_norm, ev, method=method,
        derivation=derivation or ("Inherited from the printed general procedure." if method == "inherited_general_condition" else None))
    if temperature is not None:
        s["temperature"] = fact(temperature[0], temperature[1], ev, method=method,
            derivation=derivation or ("Inherited from the printed general procedure." if method == "inherited_general_condition" else None))
    if time is not None:
        s["time"] = fact(time[0], time[1], ev, method=method,
            derivation=derivation or ("Inherited from the printed general procedure." if method == "inherited_general_condition" else None))
    if atmosphere is not None:
        s["atmosphere"] = fact(atmosphere, atmosphere, ev, method=method,
            derivation=derivation or ("Inherited from the printed general procedure." if method == "inherited_general_condition" else None))
    return s


scope_ids = [rid for rid in reactions if rid.startswith(f"{DOC}:table4_")]
if len(scope_ids) != 29:
    raise RuntimeError(f"Expected 29 existing Table 4 experiments, found {len(scope_ids)}")

for rid in scope_ids:
    r = reactions[rid]
    is_scale = rid.endswith("table4_3a_scale")
    is_modified = rid.endswith("table4_3aa")
    r["evidence_ids"] = list(dict.fromkeys(r["evidence_ids"] + P7))
    if is_scale:
        procedure_raw = {
            "page_7_general_procedure": GENERAL_TEXT,
            "Table_4_scale_exception": "1.0 g (9.4 mmol) scale",
        }
        procedure_norm = (
            "General-procedure operation sequence with the printed 9.4 mmol scale exception; "
            "scale-specific reagent, solvent, and workup quantities are unreported."
        )
    elif is_modified:
        procedure_raw = {
            "page_7_general_procedure": GENERAL_TEXT,
            "Table_4_condition_exception": "catalyst 7g (10 mol%) in EtCN at -60 °C for 18 h",
        }
        procedure_norm = (
            "General-procedure operation/workup sequence with the printed 3aa exceptions: "
            "catalyst 7g (10 mol%), EtCN, -60 °C, and 18 h; EtCN volume is unreported."
        )
    else:
        procedure_raw = GENERAL_TEXT
        procedure_norm = GENERAL_TEXT
    r["procedure"]["original_text"] = fact(
        procedure_raw, procedure_norm, P5 + P7 if (is_scale or is_modified) else P7,
        method="inherited_general_condition",
        derivation="The printed General Procedure is associated with the Table 4 asymmetric-propargylation series; explicit Table 4 exceptions remain controlling.")
    r["procedure"]["general_procedure_reference"] = fact(
        "General Procedure for Asymmetric Propargylation, page 7",
        {"document_id": DOC, "page": 7, "section": "Experimental Section"},
        P7, method="inherited_general_condition",
        derivation="Cross-page continuation from the Table 4 scope to the paper's named general procedure.")
    r["procedure"]["addition_order"] = fact(
        "The solution of allenyltrichlorosilane ... was added to a solution of catalyst, diisopropylethylamine and aldehyde",
        ["prepare allenyltrichlorosilane solution", "prepare catalyst/DIPEA/aldehyde solution", "add allenyltrichlorosilane solution to the latter at -60 °C under argon"],
        P7, method="inherited_general_condition",
        derivation="Inherited addition order; Table 4 exceptions modify catalyst/solvent/scale but do not print a replacement order of addition.")

    substrate = participant_by_id(r, "substrate")
    allenylsilane = participant_by_id(r, "allenylsilane")
    base = participant_by_id(r, "base")
    catalyst = participant_by_id(r, "catalyst")
    solvent = participant_by_id(r, "solvent")

    # Page 7 independently confirms the 1.5 equiv reagent doses and ordinary
    # 1 mol% catalyst loading. Exact small-scale charges are not imposed on the
    # 9.4 mmol exception, and the ordinary catalyst charge is not imposed on 3aa.
    append_evidence(allenylsilane["charges"]["equivalents"], P7)
    append_evidence(base["charges"]["equivalents"], P7)
    if not is_modified:
        append_evidence(catalyst["charges"]["mol_percent"], P7)

    exact_small_scale = not is_scale
    if exact_small_scale:
        conflict_scale(r)
        set_amount(allenylsilane, "0.71 mmol", 0.71)
        set_amount(base, "0.71 mmol", 0.71)
        if not is_modified:
            set_amount(catalyst, "0.0047 mmol", 0.0047)

    # Exact DCM portions only belong to the ordinary DCM procedure. The EtCN
    # exception and scale-up keep their unreported experiment-specific volumes.
    reaction_solvent_ids = ["solvent"]
    if not is_scale and not is_modified:
        solvent["charges"]["volume"] = fact(
            "DCM (1.5 mL)", {"number": 1.5, "unit": "mL"}, P7,
            method="inherited_general_condition",
            derivation="Main reaction-solution DCM portion in the printed general procedure.")
        carrier = add_participant(r, "allenylsilane_carrier_dcm", "m_dcm", "solvent", "reaction",
                                  volume=volume("0.5 mL DCM", 0.5))
        reaction_solvent_ids.append(carrier["participant_id"])
        r["solvent_conditions"]["composition"] = fact(
            "allenyltrichlorosilane in 0.5 mL DCM; catalyst/DIPEA/aldehyde in 1.5 mL DCM",
            [{"molecule_id": "m_dcm", "participant_id": "allenylsilane_carrier_dcm", "volume": {"number": 0.5, "unit": "mL"}},
             {"molecule_id": "m_dcm", "participant_id": "solvent", "volume": {"number": 1.5, "unit": "mL"}}],
            P7, method="inherited_general_condition",
            derivation="The two source-reported solvent portions are kept distinct; no measured total volume is asserted.")
        r["solvent_conditions"]["total_volume"] = missing(P7)
    elif is_scale:
        msg = "The 9.4 mmol scale-up gives no scale-specific solvent or workup volumes; page-7 small-scale volumes were not multiplied or copied."
        if msg not in r["review"]["unresolved"]:
            r["review"]["unresolved"].append(msg)
    else:
        msg = "For 3aa, Table 4 explicitly substitutes EtCN for DCM; page-7 DCM portion volumes were not copied to the modified experiment."
        if msg not in r["review"]["unresolved"]:
            r["review"]["unresolved"].append(msg)

    quantity_method = "inherited_general_condition"
    quench = add_participant(r, "quench_nh4cl", "m_saturated_nh4cl_solution", "workup_agent", "workup",
                             volume=None if is_scale else volume("2 mL", 2.0), quantity_method=quantity_method)
    extractant = add_participant(r, "extractant_ether", "m_ether", "solvent", "workup",
                                volume=None if is_scale else ("3 × 20 mL", {"number": 20, "unit": "mL", "count": 3, "total": {"number": 60, "unit": "mL"}}),
                                quantity_method=quantity_method)
    brine = add_participant(r, "brine_wash", "m_brine", "workup_agent", "workup")
    drying = add_participant(r, "drying_agent", "m_sodium_sulfate", "workup_agent", "workup")
    silica = add_participant(r, "silica_gel", "m_silica_gel", "purification_agent", "purification")
    pet_ether = add_participant(r, "petroleum_ether", "m_petroleum_ether", "purification_agent", "purification")
    etoac = add_participant(r, "ethyl_acetate_eluent", "m_etoac", "purification_agent", "purification")

    reaction_ids = ["substrate", "allenylsilane", "base", "catalyst"] + reaction_solvent_ids
    modified_derivation = None
    modified_ev = P7
    if is_modified:
        modified_ev = P5 + P7
        modified_derivation = "Page 7 supplies the general addition sequence; Table 4 footnote c explicitly replaces catalyst/loading, solvent, and hold time with 7g (10 mol%), EtCN, and 18 h."
    if is_scale:
        modified_ev = P5 + P7
        modified_derivation = "Page 7 supplies the operation sequence; Table 4 footnote b explicitly replaces the small-scale aldehyde charge with 1.0 g (9.4 mmol)."
    hold_hours = 18 if is_modified else 8
    r["steps"] = [
        step("low_temperature_addition", "reaction",
             "allenyltrichlorosilane solution was added to catalyst/DIPEA/aldehyde solution under argon at -60 °C",
             "low-temperature addition in the printed order", reaction_ids,
             temperature=("-60 °C", {"number": -60, "unit": "°C"}), atmosphere="argon",
             ev=modified_ev, derivation=modified_derivation),
        step("constant_temperature_hold", "reaction",
             "The mixture was stirred at the same temperature for 8 hours" if not is_modified else
             "Table 4 footnote c: -60 °C for 18 h (modification of the general procedure)",
             "stir/hold at -60 °C", reaction_ids,
             temperature=("-60 °C", {"number": -60, "unit": "°C"}),
             time=(f"{hold_hours} h", {"number": hold_hours, "unit": "h"}), atmosphere="argon",
             ev=modified_ev, derivation=modified_derivation),
        step("quench", "workup",
             "general-procedure quench with saturated NH4Cl solution; scale-specific volume unreported" if is_scale else
             "quenched with saturated NH4Cl solution (2 mL)",
             "quench with saturated ammonium chloride solution", [quench["participant_id"]]),
        step("extraction", "workup",
             "general-procedure extraction with Et2O; scale-specific extraction volume unreported" if is_scale else
             "aqueous layer extracted with Et2O (3 × 20 mL)",
             "extract aqueous layer three times with diethyl ether", [extractant["participant_id"]]),
        step("brine_wash", "workup", "combined organic extracts were washed with brine",
             "wash combined organic extracts with brine", [brine["participant_id"]]),
        step("drying", "workup", "dried over Na2SO4", "dry organic extracts over sodium sulfate",
             [drying["participant_id"]]),
        step("solvent_removal", "workup", "The solvent was removed in vacuo",
             "remove solvent under reduced pressure", []),
        step("flash_chromatography", "purification",
             "purified by flash chromatography on silica gel with petroleum ether-ethyl acetate 9:1",
             "flash chromatography on silica gel; petroleum ether/ethyl acetate 9:1",
             [silica["participant_id"], pet_ether["participant_id"], etoac["participant_id"]]),
        step("stereochemical_characterization", "characterization",
             "enantiopurity determined by GC or HPLC with chiral sorbents; absolute configurations assigned by comparing optical rotation values to the literature",
             "determine enantiopurity by chiral GC or HPLC and assign absolute configuration by literature optical-rotation comparison",
             ["product"]),
    ]

    # Eluent ratio is an explicit purification condition, not a reaction-solvent ratio.
    r["notes"].append(fact(
        "petroleum ether-ethyl acetate mixture 9:1",
        {"stage": "purification", "components": ["petroleum ether", "ethyl acetate"], "ratio": "9:1"},
        P7))

    # Page 7 resolves the method class and the assignment basis, but the source
    # still leaves the product-specific choice of GC versus HPLC in unavailable SI.
    for outcome in r["outcomes"]:
        if outcome["kind"] == "selectivity":
            append_evidence(outcome["measurement_basis"], P7)
            outcome["measurement_basis"]["raw_value"] = (
                "Table 4: chiral GC or chiral HPLC, see SI; General Procedure: GC or HPLC with chiral sorbents; absolute configurations assigned by optical rotation comparison to literature"
            )
            outcome["measurement_basis"]["normalized_value"] = {
                "enantiopurity_method": "chiral GC or chiral HPLC; product-specific method unavailable without SI",
                "absolute_configuration_assignment": "comparison of optical rotation values to literature",
            }

    # The general procedure adds workup/purification that was previously absent;
    # retain only still-valid limitations.
    r["review"]["unresolved"] = [
        x for x in r["review"]["unresolved"]
        if x != "SI was advertised but not supplied or explicitly associated; compound-specific operations and the choice of chiral GC versus HPLC remain unavailable."
    ]
    si_msg = "SI remains unassociated; product-specific characterization details and the choice between chiral GC and HPLC remain unavailable."
    if si_msg not in r["review"]["unresolved"]:
        r["review"]["unresolved"].insert(0, si_msg)
    r["review"]["status"] = "needs_review"
    r["completeness"] = "partial"


contexts = {c.get("context_id") for c in bundle["source_context"]}
new_contexts = [
    {
        "context_id": f"{DOC}:page7_nonreaction_material",
        "document_id": DOC,
        "kind": "nonreaction_context",
        "evidence_ids": P7,
        "note": "Page 7 also contains the end of the conclusion, acknowledgements, and references [1]-[16]; these do not create additional experiments. The General Procedure is linked to the existing Table 4 experiments rather than duplicated as a new reaction.",
    },
    {
        "context_id": f"{DOC}:page8_references",
        "document_id": DOC,
        "kind": "bibliography",
        "evidence_ids": P8,
        "note": "Page 8 contains only references [17]-[21], including software citations; no reaction arrow, table, experimental procedure, or characterization peak list is present.",
    },
]
bundle["source_context"].extend(c for c in new_contexts if c["context_id"] not in contexts)

# Preserve a raw pre-merge candidate separately from the formal submission.
raw_path = RUN / "raw_model/extracted_bundle_pages_0007_0008.json"
raw_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
(RUN / "submitted.json").write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

diag = st.stereo_inspect("[Na+].[Na+].[O-]S(=O)(=O)[O-]")
(RUN / "raw_model/stereo_diagnostics_pages_0007_0008.json").write_text(
    json.dumps({"m_sodium_sulfate": diag}, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "updated_scope_reactions": len(scope_ids),
    "new_molecules": 5,
    "total_reactions": len(bundle["reactions"]),
    "total_molecules": len(bundle["molecules"]),
    "source_context": len(bundle["source_context"]),
}, ensure_ascii=False))
