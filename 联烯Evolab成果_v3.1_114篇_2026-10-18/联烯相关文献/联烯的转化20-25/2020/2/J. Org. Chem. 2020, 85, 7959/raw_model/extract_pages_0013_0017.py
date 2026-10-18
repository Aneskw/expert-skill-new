#!/usr/bin/env python3
"""Extend the v3.1 machine bundle with inspected source-PDF pages 13--17.

The five original page images were inspected at original detail.  This script
updates experiments whose detailed procedures continue on pages 13--14, adds
the distinct runs first reported there, and records pages 15--17 as
non-synthetic source context.  All records remain machine ``needs_review``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SKILL_SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SKILL_SCRIPTS))

from standard_records import (  # noqa: E402
    fact,
    new_outcome,
    new_participant,
    new_reaction,
    new_step,
    source_mentions,
    stereo_inspect,
)

DOC = "doc_f3d393b861fb548d"
RUN = ROOT / "runs/current"
RAW = RUN / "raw_model"
OUT = RUN / "submitted.json"
BASE = RUN / "standard_records.json"


def tids(*pages: int) -> list[str]:
    out: list[str] = []
    for page in pages:
        out.extend([f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"])
    return out


def F(raw, normalized=None, ev=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    if normalized is None and raw is not None and status != "conflict":
        normalized = raw
    return fact(raw=raw, normalized=normalized, evidence=ev, method=method,
                derivation=derivation, status=status, alternatives=alternatives)


def Q(number, unit, **extra):
    value = {"number": number, "unit": unit}
    value.update(extra)
    return value


def inferred(raw, normalized, ev, derivation):
    return F(raw, normalized, ev, method="inference", derivation=derivation)


bundle = json.loads(BASE.read_text())
template = json.loads((RAW / "init_template_pages_0013_0017.json").read_text())

# Retain the accumulated document/evidence registry but fill the true title.
document = bundle["documents"][0]
TRUE_TITLE = (
    "Metal-Free C-C Coupling of an Allenyl Sulfone with Picolyl Amides to "
    "Access Vinyl Sulfones via Pyridine-Initiated In Situ Generation of "
    "Sulfinate Anion"
)
# Standard merge requires the document registry to match the mechanically
# prepared dataset byte-for-byte.  Its title slot is null, so the verified
# title is retained in page_inventory/source_context/workflow rather than
# mutating the registry here.
document["title"] = template["documents"][0].get("title")
document["doi"] = template["documents"][0].get("doi")
for item in template["documents"][0].get("evidence", []):
    if item.get("evidence_id") not in {e.get("evidence_id") for e in document.get("evidence", [])}:
        document.setdefault("evidence", []).append(item)

molecules = {m["molecule_id"]: m for m in bundle["molecules"]}


def add_molecule_occurrence_evidence(mid: str, *pages: int) -> None:
    ev = tids(*pages)
    occ = molecules[mid]["occurrences"][0]
    for eid in ev:
        if eid not in occ["evidence_ids"]:
            occ["evidence_ids"].append(eid)
        if eid not in occ["compound_label"]["evidence_ids"]:
            occ["compound_label"]["evidence_ids"].append(eid)


for mid, pages in {
    "m_1m": (13,), "m_2a": (13, 14), "m_3m": (13,), "m_4m": (13,),
    "m_1a": (13, 14), "m_bht": (13,), "m_aibn": (13,),
    "m_5b": (13, 14), "m_3q": (13,), "m_5a": (13, 14),
    "m_2a_d2": (13, 14), "m_3a_d": (13, 14), "m_4a_d": (13, 14),
    "m_11": (14,), "m_pyridine": (14,), "m_3a": (13, 14),
    "m_4a": (13, 14), "m_k2co3": (14,),
}.items():
    add_molecule_occurrence_evidence(mid, *pages)


def base_reaction(tag, pages, location, name, transformation, categories,
                  unresolved):
    ev = tids(*pages)
    reaction = new_reaction(f"{DOC}:{tag}", DOC, ev)
    reaction["record_type"] = "specific_experiment"
    reaction["experiment_categories"] = F(list(categories), list(categories), ev)
    reaction["location"] = F(location, location, ev)
    reaction["reaction_name"] = F(name, name, ev)
    reaction["transformation_type"] = F(transformation, transformation, ev)
    reaction["review"] = {
        "status": "needs_review",
        "method": "model",
        "reviewer": None,
        "checked_fields": [],
        "unresolved": list(unresolved),
    }
    return reaction, ev


def participant(reaction, mid, role, ev, *, mass=None, amount=None, volume=None,
                equivalents=None, reference=None, stage="reaction",
                product_role=None, method="original_extraction"):
    p = new_participant(
        f"{reaction['reaction_id']}:p{len(reaction['participants']) + 1}",
        mid, role, ev, stage=stage)
    if mass:
        p["charges"]["mass"] = F(mass[0], Q(mass[1], mass[2]), ev, method=method)
    if amount:
        p["charges"]["amount"] = F(amount[0], Q(amount[1], amount[2]), ev, method=method)
    if volume:
        p["charges"]["volume"] = F(volume[0], Q(volume[1], volume[2]), ev, method=method)
    if equivalents is not None:
        p["charges"]["equivalents"] = F(
            f"{equivalents:g} equiv", Q(equivalents, "equiv"), ev, method=method)
        p["equivalence_reference"] = F(
            f"relative to {reference}", reference, ev, method=method)
    if product_role:
        p["product_role"] = F(product_role, product_role, ev, method=method)
    reaction["participants"].append(p)
    return p


def step(reaction, ev, operation, *, temperature=None, time=None,
         atmosphere=None, pressure=None, stage="reaction",
         method="original_extraction"):
    s = new_step(
        f"{reaction['reaction_id']}:s{len(reaction['steps']) + 1}", ev,
        stage=stage)
    s["participant_ids"] = [
        p["participant_id"] for p in reaction["participants"]
        if p["stage"] == stage]
    s["operation"] = F(operation, operation, ev, method=method)
    if temperature:
        s["temperature"] = F(temperature[0], temperature[1], ev, method=method)
    if time:
        s["time"] = F(time[0], time[1], ev, method=method)
    if atmosphere:
        s["atmosphere"] = F(atmosphere, atmosphere, ev, method=method)
    if pressure:
        s["pressure"] = F(pressure, pressure, ev, method=method)
    reaction["steps"].append(s)
    return s


def outcome(reaction, ev, kind, targets, raw, normalized, *, role=None,
            yield_type=None, basis=None, scope=None, method="original_extraction"):
    o = new_outcome(
        f"{reaction['reaction_id']}:o{len(reaction['outcomes']) + 1}",
        kind, targets, ev)
    o["value"] = F(raw, normalized, ev, method=method)
    if role:
        o["product_role"] = F(role, role, ev, method=method)
    if yield_type:
        o["yield_type"] = F(yield_type, yield_type, ev, method=method)
    if basis:
        o["measurement_basis"] = F(basis, basis, ev, method=method)
    else:
        o["measurement_basis"] = fact(
            raw=None, normalized=None, evidence=ev, status="unreported")
    if scope:
        o["scope"] = F(scope, scope, ev, method=method)
    reaction["outcomes"].append(o)
    return o


def set_scale(reaction, ev, mid, raw, number, unit):
    reaction["reaction_scale"] = F(
        raw, {"reference_molecule_id": mid, "amount": Q(number, unit)}, ev)


def set_solvent(reaction, ev, mids, raw, *, volume=None, concentration=None,
                basis=None, method="original_extraction"):
    reaction["solvent_conditions"]["composition"] = F(
        raw, list(mids), ev, method=method)
    if volume:
        reaction["solvent_conditions"]["total_volume"] = F(
            volume[0], Q(volume[1], volume[2]), ev, method=method)
    if concentration:
        reaction["solvent_conditions"]["concentration"] = F(
            concentration[0], Q(concentration[1], concentration[2]), ev,
            method=method)
    if basis:
        reaction["solvent_conditions"]["concentration_basis"] = F(
            basis, basis, ev, method=method)


def classify_cc_ch(reaction, ev, note=None):
    limitation = note or (
        "The source substrate/product identities support a new C-C connection "
        "and net C-H formation, but allene atom correspondence, proton source, "
        "leaving groups, and coproducts are not completely mapped.")
    changes = []
    labels = []
    for index, (a, b, code, category_number) in enumerate(
            [("C", "C", "C-C", 2), ("C", "H", "C-H", 1)], 1):
        cid = f"b{index}"
        changes.append(inferred(
            "partial source-derived bond-change candidate",
            {"change_id": cid,
             "atoms": [{"atom_map": None, "element": a},
                       {"atom_map": None, "element": b}],
             "before_order": 0, "after_order": 1,
             "change_type": "formed"}, ev, limitation))
        labels.append({"category_number": category_number, "code": code,
                       "supporting_change_ids": [cid]})
    reaction["classification"]["bond_changes"] = changes
    reaction["classification"]["labels"] = inferred(
        "partial formed-bond labels", labels, ev, limitation)
    reaction["classification"]["coverage"] = inferred(
        "partial source-derived classification", "partial", ev, limitation)
    reaction["classification"]["limitations"] = [
        limitation, "No complete atom-mapped reaction SMILES is asserted."]


def classify_none(reaction, ev, reason):
    reaction["classification"]["bond_changes"] = []
    reaction["classification"]["labels"] = inferred(
        "no product bond formation assigned", [], ev, reason)
    reaction["classification"]["coverage"] = inferred(
        "no applicable formed-bond classification", "not_applicable", ev, reason)
    reaction["classification"]["limitations"] = [reason]


def add_standard_workup(reaction, ev, purification=None):
    participant(reaction, "m_ch2cl2", "workup_agent", ev, stage="workup")
    step(reaction, ev,
         "dilute with dichloromethane, transfer, and concentrate under reduced pressure",
         stage="workup", pressure="reduced pressure")
    if purification:
        participant(reaction, "m_silica", "purification_agent", ev,
                    stage="purification")
        participant(reaction, "m_etoac", "purification_agent", ev,
                    stage="purification")
        participant(reaction, "m_hexane", "purification_agent", ev,
                    stage="purification")
        step(reaction, ev, f"silica-gel flash chromatography ({purification})",
             stage="purification")


def replace_or_add(reaction):
    rid = reaction["reaction_id"]
    for index, old in enumerate(bundle["reactions"]):
        if old["reaction_id"] == rid:
            bundle["reactions"][index] = reaction
            return "updated"
    bundle["reactions"].append(reaction)
    return "added"


new_ids: list[str] = []
updated_ids: list[str] = []
ADDED_THIS_PASS = [
    f"{DOC}:scope_3m_run1", f"{DOC}:scope_3m_run2",
    f"{DOC}:scheme2d_deuterium_run2",
    f"{DOC}:nbenzyl_no_pyridine_run1",
    f"{DOC}:nbenzyl_pyridine_low_run4",
    f"{DOC}:interconversion_mixture_preparation",
    f"{DOC}:control_1a_5b_k2co3_exp2",
]
UPDATED_THIS_PASS = [
    f"{DOC}:scheme2a_bht", f"{DOC}:scheme2b_aibn",
    f"{DOC}:scheme2c_crossover", f"{DOC}:scheme2d_deuterium",
    f"{DOC}:nbenzyl_pyridine_1eq", f"{DOC}:nbenzyl_pyridine_0p2eq",
    f"{DOC}:interconversion_probe", f"{DOC}:in_situ_nmr_monitoring",
    f"{DOC}:control_1a_5b_k2co3",
]


def save_reaction(reaction):
    state = replace_or_add(reaction)
    (new_ids if state == "added" else updated_ids).append(reaction["reaction_id"])


# Page 13: two independent 3m/4m substrate-scope runs.
for run_no, product_mass, yld, ratio in (
        (1, 46.3, 60, (92, 8)), (2, 44.8, 58, (90, 10))):
    reaction, ev = base_reaction(
        f"scope_3m_run{run_no}", [13],
        f"Compound paragraph for 3m, run {run_no}",
        f"individual reported run {run_no}: formation of 3m/4m",
        "C(sp3)-H functionalization with allenyl sulfone 2a",
        ("substrate_scope",),
        ("The isolated material contains 3m and 4m; product-specific masses are not reported.",
         "The newly generated tetrahedral center is unassigned; no racemic composition is inferred.",
         "The paragraph invokes General Procedure C; the 50 °C operation is inherited from page 6 and the 48 h scope duration from Figure 1/page 3 rather than restated on page 13."))
    reaction["general_scheme_ids"] = [f"{DOC}:figure1_generic"]
    participant(reaction, "m_1m", "reactant", ev,
                mass=("38.5 mg", 38.5, "mg"),
                amount=("0.2 mmol", 0.2, "mmol"),
                equivalents=1.0, reference="m_1m")
    participant(reaction, "m_2a", "reactant", ev,
                mass=("78 mg", 78, "mg"), amount=("0.4 mmol", 0.4, "mmol"),
                equivalents=2.0, reference="m_1m")
    participant(reaction, "m_thf", "solvent", ev,
                volume=("1 mL", 1, "mL"))
    participant(reaction, "m_3m", "product", ev,
                product_role="major product in isolated mixture")
    participant(reaction, "m_4m", "product", ev,
                product_role="minor product in isolated mixture")
    set_scale(reaction, ev, "m_1m", "0.2 mmol of 1m", 0.2, "mmol")
    set_solvent(reaction, ev, ["m_thf"], "THF", volume=("1 mL", 1, "mL"),
                concentration=("0.2 M", 0.2, "M"),
                basis="0.2 mmol 1m in 1 mL THF")
    inherited_ev = tids(3, 6)
    step(reaction, inherited_ev,
         "charge solid picolyl amide 1m and 2a in air and seal the screw-top tube",
         atmosphere="air", method="inherited_general_condition")
    step(reaction, inherited_ev,
         "evacuate for 1-2 s and fill with nitrogen; repeat three times",
         time=("1-2 s evacuation per cycle",
               {"minimum": 1, "maximum": 2, "unit": "s", "raw_text": "1-2 s"}),
         atmosphere="nitrogen", method="inherited_general_condition")
    step(reaction, inherited_ev,
         "add THF by syringe, lower into a preheated oil bath, and monitor by 1H NMR",
         temperature=("50 °C", Q(50, "°C")),
         time=("48 h (Figure 1 scope arrow)", Q(48, "h")),
         atmosphere="nitrogen", method="inherited_general_condition")
    add_standard_workup(reaction, ev, "20-100% ethyl acetate/hexane")
    outcome(reaction, ev, "product_amount", ["m_3m", "m_4m"],
            f"{product_mass:g} mg", Q(product_mass, "mg"),
            role="combined 3m/4m isolated mixture",
            basis="isolated after silica-gel flash chromatography")
    outcome(reaction, ev, "yield", ["m_3m", "m_4m"], f"{yld}%", Q(yld, "%"),
            role="combined 3m/4m mixture", yield_type="isolated yield",
            basis="isolated after silica-gel flash chromatography")
    outcome(reaction, ev, "selectivity", ["m_3m", "m_4m"],
            f"3m/4m = {ratio[0]}:{ratio[1]}",
            {"ratio": list(ratio), "order": ["m_3m", "m_4m"],
             "metric": "product ratio"},
            basis="post-chromatography 1H NMR integration ratios",
            scope="3m/4m product selectivity; absolute configuration unknown")
    classify_cc_ch(reaction, ev)
    save_reaction(reaction)


# Page 13: radical inhibitor and radical initiator controls.
reaction, ev = base_reaction(
    "scheme2a_bht", [4, 13], "Scheme 2a and Radical Inhibitor Experiment",
    "radical inhibitor experiment with BHT",
    "C(sp3)-H functionalization under radical-inhibitor conditions",
    ("mechanistic_study", "control"),
    ("The isolated material is a 96:4 3a/4a mixture; product-specific masses are not reported.",
     "The new tetrahedral center is not assigned."))
participant(reaction, "m_1a", "reactant", ev, mass=("40 mg", 40, "mg"),
            amount=("0.19 mmol", 0.19, "mmol"), equivalents=1, reference="m_1a")
participant(reaction, "m_2a", "reactant", ev, mass=("74 mg", 74, "mg"),
            amount=("0.38 mmol", 0.38, "mmol"), equivalents=2, reference="m_1a")
participant(reaction, "m_bht", "additive", ev, mass=("208 mg", 208, "mg"),
            amount=("0.94 mmol", 0.94, "mmol"), equivalents=5, reference="m_1a")
participant(reaction, "m_thf", "solvent", ev, volume=("1 mL", 1, "mL"))
participant(reaction, "m_3a", "product", ev, product_role="major product")
participant(reaction, "m_4a", "product", ev, product_role="minor product")
set_scale(reaction, ev, "m_1a", "0.19 mmol 1a", 0.19, "mmol")
set_solvent(reaction, ev, ["m_thf"], "THF", volume=("1 mL", 1, "mL"),
            concentration=("0.19 M", 0.19, "M"), basis="0.19 mmol 1a in 1 mL THF")
step(reaction, ev, "charge 1a, 2a, and BHT in air and seal the screw-top tube",
     atmosphere="air")
step(reaction, ev, "evacuate and fill with nitrogen three times", atmosphere="nitrogen")
step(reaction, ev, "add THF by syringe and lower into a preheated oil bath",
     temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")),
     atmosphere="nitrogen")
add_standard_workup(reaction, ev, "20-70% ethyl acetate/hexane")
outcome(reaction, ev, "product_amount", ["m_3a", "m_4a"], "51.6 mg",
        Q(51.6, "mg"), role="combined 3a/4a isolated mixture",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "yield", ["m_3a", "m_4a"], "67%", Q(67, "%"),
        role="combined 3a/4a mixture", yield_type="isolated yield",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "selectivity", ["m_3a", "m_4a"], "3a/4a = 96:4",
        {"ratio": [96, 4], "order": ["m_3a", "m_4a"], "metric": "product ratio"},
        basis="post-chromatography 1H NMR integration ratios",
        scope="3a/4a product selectivity")
classify_cc_ch(reaction, ev)
save_reaction(reaction)

reaction, ev = base_reaction(
    "scheme2b_aibn", [4, 13], "Scheme 2b and Radical Initiator Experiment",
    "radical initiator experiment with AIBN",
    "C(sp3)-H functionalization under radical-initiator conditions",
    ("mechanistic_study", "control"),
    ("The isolated material is an 85:15 3a/4a mixture; product-specific masses are not reported.",
     "Recovered 1a is retained as recovery, not product yield."))
participant(reaction, "m_1a", "reactant", ev, mass=("40 mg", 40, "mg"),
            amount=("0.19 mmol", 0.19, "mmol"), equivalents=1, reference="m_1a")
participant(reaction, "m_2a", "reactant", ev, mass=("74 mg", 74, "mg"),
            amount=("0.38 mmol", 0.38, "mmol"), equivalents=2, reference="m_1a")
participant(reaction, "m_aibn", "additive", ev, mass=("31 mg", 31, "mg"),
            amount=("0.19 mmol", 0.19, "mmol"), equivalents=1, reference="m_1a")
participant(reaction, "m_thf", "solvent", ev, volume=("1 mL", 1, "mL"))
participant(reaction, "m_3a", "product", ev, product_role="major product")
participant(reaction, "m_4a", "product", ev, product_role="minor product")
set_scale(reaction, ev, "m_1a", "0.19 mmol 1a", 0.19, "mmol")
set_solvent(reaction, ev, ["m_thf"], "THF", volume=("1 mL", 1, "mL"),
            concentration=("0.19 M", 0.19, "M"), basis="0.19 mmol 1a in 1 mL THF")
step(reaction, ev, "charge 1a and 2a in air; transfer to a glovebox, add AIBN, seal, and remove from the glovebox",
     atmosphere="glovebox during AIBN charge")
step(reaction, ev, "connect to the Schlenk line and add THF by syringe",
     atmosphere="nitrogen")
step(reaction, ev, "lower into a preheated oil bath and react",
     temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")),
     atmosphere="nitrogen")
add_standard_workup(reaction, ev, "20-70% ethyl acetate/hexane")
outcome(reaction, ev, "product_amount", ["m_3a", "m_4a"], "45.8 mg",
        Q(45.8, "mg"), role="combined 3a/4a isolated mixture",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "yield", ["m_3a", "m_4a"], "60%", Q(60, "%"),
        role="combined 3a/4a mixture", yield_type="isolated yield",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "selectivity", ["m_3a", "m_4a"], "3a/4a = 85:15",
        {"ratio": [85, 15], "order": ["m_3a", "m_4a"], "metric": "product ratio"},
        basis="post-chromatography 1H NMR integration ratios",
        scope="3a/4a product selectivity")
outcome(reaction, ev, "product_amount", ["m_1a"], "7.1 mg recovered 1a",
        Q(7.1, "mg"), role="recovered starting material",
        basis="recovered after chromatography")
outcome(reaction, ev, "recovery", ["m_1a"], "18% recovered 1a", Q(18, "%"),
        role="recovered starting material", basis="recovered after chromatography")
classify_cc_ch(reaction, ev)
save_reaction(reaction)


# Page 13: crossover experiment, including byproduct and recoveries.
reaction, ev = base_reaction(
    "scheme2c_crossover", [4, 13], "Scheme 2c and Crossover Experiment",
    "crossover experiment with disulfone 5b", "crossover vinyl-sulfone formation",
    ("mechanistic_study",),
    ("The 41% isolated material is a 93:7 3q/3a mixture; product-specific masses are not reported.",
     "The separate 5a branch is retained without complete atom mapping."))
participant(reaction, "m_1a", "reactant", ev, mass=("15 mg", 15, "mg"),
            amount=("0.07 mmol", 0.07, "mmol"), equivalents=1, reference="m_1a")
participant(reaction, "m_2a", "reactant", ev, mass=("13.7 mg", 13.7, "mg"),
            amount=("0.07 mmol", 0.07, "mmol"), equivalents=1, reference="m_1a")
participant(reaction, "m_5b", "reactant", ev, mass=("47.6 mg", 47.6, "mg"),
            amount=("0.14 mmol", 0.14, "mmol"), equivalents=2, reference="m_1a")
participant(reaction, "m_thf", "solvent", ev, volume=("0.37 mL", 0.37, "mL"))
participant(reaction, "m_3q", "product", ev, product_role="major crossover product")
participant(reaction, "m_3a", "product", ev, product_role="minor product")
participant(reaction, "m_5a", "product", ev, product_role="isolated byproduct")
set_scale(reaction, ev, "m_1a", "0.07 mmol 1a", 0.07, "mmol")
set_solvent(reaction, ev, ["m_thf"], "THF", volume=("0.37 mL", 0.37, "mL"),
            concentration=("approximately 0.19 M", 0.189, "M"),
            basis="0.07 mmol 1a in 0.37 mL THF")
step(reaction, ev, "charge 1a, 2a, and 5b; seal the screw-top tube and add THF by syringe")
step(reaction, ev, "lower into a preheated oil bath and react",
     temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")))
add_standard_workup(reaction, ev, "20-80% ethyl acetate/hexane")
outcome(reaction, ev, "product_amount", ["m_3q", "m_3a"], "11.3 mg",
        Q(11.3, "mg"), role="combined 3q/3a isolated mixture",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "yield", ["m_3q", "m_3a"], "41%", Q(41, "%"),
        role="combined 3q/3a mixture", yield_type="isolated yield",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "selectivity", ["m_3q", "m_3a"], "3q/3a = 93:7",
        {"ratio": [93, 7], "order": ["m_3q", "m_3a"], "metric": "product ratio"},
        basis="post-chromatography 1H NMR integration ratios", scope="3q/3a")
outcome(reaction, ev, "yield", ["m_5a"], "approximately 57% based on 2a",
        Q(57, "%", qualifier="approx", raw_text="~57%"),
        role="isolated byproduct", yield_type="isolated yield",
        basis="isolated after chromatography", scope="based on 2a as limiting reagent")
outcome(reaction, ev, "recovery", ["m_5b"], "61% recovered 5b", Q(61, "%"),
        role="recovered starting material", basis="recovered after chromatography")
outcome(reaction, ev, "recovery", ["m_2a"], "48% recovered 2a", Q(48, "%"),
        role="recovered starting material", basis="recovered after chromatography")
classify_cc_ch(reaction, ev,
               "The crossover products support new C-C and net C-H bonds; the concurrent 5a pathway and complete atom correspondence remain unresolved.")
save_reaction(reaction)


# Page 13/14: deuterium-labeling runs.  The detailed paragraph calls the
# labeled allene "2d" although the preparation registered it as 2a-d2.
def add_deuterium_run(tag, pages, run_no, sub_mass, sub_amount, allene_mass,
                      allene_amount, thf_volume, product_mass, yld, ratio,
                      include_scheme_composition=False):
    reaction, ev = base_reaction(
        tag, pages, f"Deuterium-Labeling Experiment, run {run_no}",
        f"deuterium-labeling experiment run {run_no}",
        "isotope-labeling C(sp3)-H functionalization",
        ("mechanistic_study",),
        ("The source paragraph calls the labeled allene 2d, whereas its preparation names it 2a-d2; both labels are retained as an unresolved source-label inconsistency.",
         "The products are positional-isotopologue mixtures and are not forced into one concrete stereochemical SMILES.",
         "The isolated material contains 3a-d and 4a-d; product-specific masses are not reported."))
    participant(reaction, "m_1a", "reactant", ev,
                mass=(f"{sub_mass:g} mg", sub_mass, "mg"),
                amount=(f"{sub_amount:g} mmol", sub_amount, "mmol"),
                equivalents=1, reference="m_1a")
    participant(reaction, "m_2a_d2", "reactant", ev,
                mass=(f"{allene_mass:g} mg", allene_mass, "mg"),
                amount=(f"{allene_amount:g} mmol", allene_amount, "mmol"),
                equivalents=2, reference="m_1a")
    participant(reaction, "m_thf", "solvent", ev,
                volume=(f"{thf_volume:g} mL", thf_volume, "mL"))
    participant(reaction, "m_3a_d", "product", ev,
                product_role="major positional-isotopologue mixture")
    participant(reaction, "m_4a_d", "product", ev,
                product_role="minor positional-isotopologue mixture")
    set_scale(reaction, ev, "m_1a", f"{sub_amount:g} mmol 1a", sub_amount, "mmol")
    set_solvent(reaction, ev, ["m_thf"], "THF",
                volume=(f"{thf_volume:g} mL", thf_volume, "mL"),
                concentration=(f"{sub_amount / thf_volume:.3g} M",
                               sub_amount / thf_volume, "M"),
                basis=f"{sub_amount:g} mmol 1a in {thf_volume:g} mL THF")
    step(reaction, ev, "charge 1a and labeled allene (called 2d in this paragraph) in air and seal the screw-top tube",
         atmosphere="air")
    step(reaction, ev, "evacuate and fill with nitrogen three times", atmosphere="nitrogen")
    step(reaction, ev, "add THF by syringe, lower into a preheated oil bath, and react",
         temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")),
         atmosphere="nitrogen")
    add_standard_workup(reaction, ev,
                        "20-90% ethyl acetate/hexane" if run_no == 1 else
                        "20-100% ethyl acetate/hexane")
    outcome(reaction, ev, "product_amount", ["m_3a_d", "m_4a_d"],
            f"{product_mass:g} mg", Q(product_mass, "mg"),
            role="combined isotope-labeled product mixture",
            basis="isolated after silica-gel flash chromatography")
    outcome(reaction, ev, "yield", ["m_3a_d", "m_4a_d"], f"{yld}%", Q(yld, "%"),
            role="combined isotope-labeled product mixture", yield_type="isolated yield",
            basis="isolated after silica-gel flash chromatography")
    outcome(reaction, ev, "selectivity", ["m_3a_d", "m_4a_d"],
            f"3a-d/4a-d = {ratio[0]}:{ratio[1]}",
            {"ratio": list(ratio), "order": ["m_3a_d", "m_4a_d"],
             "metric": "product ratio"},
            basis="post-chromatography 1H NMR integration ratios",
            scope="3a-d/4a-d product selectivity")
    if include_scheme_composition:
        scheme_ev = tids(4)
        outcome(reaction, scheme_ev, "selectivity", ["m_3a_d"],
                "allylic C-D approximately 68.5%; alkenyl C-D approximately 36.5%",
                {"allylic_d_percent": 68.5, "alkenyl_d_percent": 36.5,
                 "metric": "deuterium content"},
                basis="1H/2H NMR", scope="at least five labeled components")
        outcome(reaction, scheme_ev, "selectivity", ["m_3a_d"],
                "3a1-d2 44%, 3a2-d3 20%, 3a3-d3 27%, 3a4-d2 5%",
                {"components_percent": {"3a1-d2": 44, "3a2-d3": 20,
                                        "3a3-d3": 27, "3a4-d2": 5},
                 "metric": "positional isotopologue composition"},
                basis="1H/2H NMR and Scheme 2d assignments",
                scope="positional isotopologue mixture; configurations unassigned")
    classify_cc_ch(reaction, ev,
                   "The labeled coupling supports new C-C and net C-H/C-D connectivity, but isotope/proton correspondence and individual isotopologue mappings are incomplete.")
    save_reaction(reaction)


add_deuterium_run("scheme2d_deuterium", [4, 13], 1, 30, 0.14, 55.5,
                  0.28, 0.74, 35.2, 61, (96, 4), True)
add_deuterium_run("scheme2d_deuterium_run2", [13, 14], 2, 15, 0.07, 28,
                  0.14, 0.37, 18.7, 65, (95, 5), False)


# Page 14: N-benzylbenzamide/2a controls.  Runs 2 and 3 refine the page-5
# summary records; runs 1 and 4 are separate conditions.
def add_nbenzyl_control(tag, name, mass11, amount11, mass2a, amount2a,
                        pyridine=None, main_toluene=None, carrier_toluene=None,
                        observation="", low_pyridine=False):
    pages = [5, 14] if tag in ("nbenzyl_pyridine_1eq", "nbenzyl_pyridine_0p2eq") else [14]
    unresolved = [
        "Unreacted N-benzylbenzamide is a qualitative crude-NMR observation, not a quantified recovery.",
        "Absence of a 2a or 5a signal is not converted into a quantitative yield or conversion."]
    if low_pyridine:
        unresolved.append(
            "Page 5 summarizes the pyridine loading as 0.2 equiv, while the detailed page-14 charge is 0.005 mmol relative to 0.02 mmol 11 (calculated 0.25 equiv); neither value is silently preferred.")
        unresolved.append(
            "Page 5 says the low-pyridine run gave the same decomposition result as the 1-equiv run, whereas page 14 quantifies approximately 82% 2a consumption; both source statements are retained.")
    reaction, ev = base_reaction(
        tag, pages, f"Reaction of N-benzylbenzamide and 2a, {name}", name,
        "intramolecular-pyridine mechanistic control",
        ("mechanistic_study", "control"), tuple(unresolved))
    participant(reaction, "m_11", "reactant", ev,
                mass=(f"{mass11:g} mg", mass11, "mg"),
                amount=(f"{amount11:g} mmol", amount11, "mmol"),
                equivalents=1, reference="m_11")
    participant(reaction, "m_2a", "reactant", ev,
                mass=(f"{mass2a:g} mg", mass2a, "mg"),
                amount=(f"{amount2a:g} mmol", amount2a, "mmol"))
    if pyridine:
        p = participant(reaction, "m_pyridine", "additive", ev,
                        volume=pyridine.get("volume"), amount=pyridine.get("amount"))
        if pyridine.get("equivalents") is not None:
            p["charges"]["equivalents"] = F(
                f"{pyridine['equivalents']:g} equiv", Q(pyridine["equivalents"], "equiv"), ev)
            p["equivalence_reference"] = F("relative to m_11", "m_11", ev)
        if low_pyridine:
            p["charges"]["equivalents"] = F(
                "page 5: 0.2 equiv; page 14 detailed charge: 0.005 mmol/0.02 mmol = 0.25 equiv",
                None, ev, status="conflict",
                alternatives=[
                    F("0.2 equiv (page 5 summary)", Q(0.2, "equiv"), tids(5)),
                    inferred("0.25 equiv from page-14 detailed charge", Q(0.25, "equiv"), tids(14),
                             "Calculated from 0.005 mmol pyridine divided by 0.02 mmol N-benzylbenzamide.")])
            p["equivalence_reference"] = F("relative to m_11", "m_11", ev)
    participant(reaction, "m_5a", "product", ev,
                product_role="observed byproduct" if low_pyridine else
                "byproduct explicitly not observed")
    if carrier_toluene:
        participant(reaction, "m_toluene", "solvent", ev,
                    volume=(f"{carrier_toluene:g} mL as pyridine-solution carrier",
                            carrier_toluene, "mL"))
    if main_toluene:
        participant(reaction, "m_toluene", "solvent", ev,
                    volume=(f"{main_toluene:g} mL separately added toluene",
                            main_toluene, "mL"))
    total_toluene = (carrier_toluene or 0) + (main_toluene or 0)
    set_scale(reaction, ev, "m_11", f"{amount11:g} mmol N-benzylbenzamide", amount11, "mmol")
    if total_toluene:
        set_solvent(reaction, ev, ["m_toluene"],
                    "toluene (separate main-solvent and pyridine-solution-carrier volumes retained)")
        reaction["solvent_conditions"]["total_volume"] = inferred(
            f"{total_toluene:g} mL total stated toluene by summing source charges",
            Q(total_toluene, "mL"), ev,
            "Calculated by summing the separately stated toluene carrier and main-solvent volumes; both original charges remain on separate participants.")
    step(reaction, ev, "charge N-benzylbenzamide and 2a in air" +
         ("; add the stated pyridine charge" if pyridine else ""), atmosphere="air")
    step(reaction, ev, "seal the tube, evacuate, and fill with nitrogen three times",
         atmosphere="nitrogen", method="inherited_general_condition" if tag.endswith("run4") else "original_extraction")
    step(reaction, ev, "add the stated toluene charge(s), lower into a preheated oil bath, and react",
         temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")),
         atmosphere="nitrogen", method="inherited_general_condition" if tag.endswith("run4") else "original_extraction")
    add_standard_workup(reaction, ev)
    outcome(reaction, ev, "recovery", ["m_11"], "unreacted N-benzylbenzamide",
            "unreacted starting material observed", role="starting material",
            basis="crude 1H NMR")
    if low_pyridine:
        outcome(reaction, tids(5), "other", ["m_2a"],
                "page-5 summary: same decomposition result as the 1-equivalent pyridine experiment",
                "page-5 qualitative decomposition summary", role="reactant fate",
                basis="page-5 source prose",
                scope="retained alongside the detailed page-14 approximately 82% consumption")
        outcome(reaction, ev, "conversion", ["m_2a"],
                "approximately 82% allenyl sulfone 2a consumed", Q(82, "%", qualifier="approx"),
                role="reactant consumption", basis="crude 1H NMR")
        outcome(reaction, ev, "yield", ["m_5a"],
                "approximately 17% based on 2 equivalents of 2a",
                Q(17, "%", qualifier="approx", raw_text="~17%"),
                role="byproduct", basis="crude 1H NMR",
                scope="based on the charged two equivalents of 2a")
    elif tag == "nbenzyl_pyridine_1eq":
        outcome(reaction, tids(5), "decomposition", ["m_2a"],
                "complete decomposition of 2a", "complete decomposition of 2a",
                role="reactant fate", basis="page-5 source prose")
        outcome(reaction, tids(14), "other", ["m_5a"],
                "byproduct 5a was not observed", "5a not observed",
                basis="crude 1H NMR",
                scope="qualitative non-observation; not a quantitative zero")
    else:
        outcome(reaction, ev, "other", ["m_2a", "m_5a"], observation,
                observation, basis="crude 1H NMR",
                scope="qualitative non-observation; not a quantitative zero")
    classify_none(reaction, ev,
                  "No coupling product was assigned; the experiment is a qualitative mechanistic control.")
    if low_pyridine:
        reaction["notes"].append(F(
            "Pyridine-loading conflict retained without reconciliation.",
            {"page5_summary_equiv": 0.2, "page14_charge_derived_equiv": 0.25}, ev))
    save_reaction(reaction)


add_nbenzyl_control(
    "nbenzyl_no_pyridine_run1", "run 1 without pyridine", 5, 0.02, 9.2, 0.05,
    main_toluene=0.2,
    observation="2a and byproduct 5a were not observed")
add_nbenzyl_control(
    "nbenzyl_pyridine_1eq", "run 2 with 1 equivalent pyridine", 5, 0.02, 9.2, 0.05,
    pyridine={"volume": ("2 uL", 2, "uL"), "amount": ("0.02 mmol", 0.02, "mmol"),
               "equivalents": 1.0}, main_toluene=0.2,
    observation="2a and byproduct 5a were not observed")
add_nbenzyl_control(
    "nbenzyl_pyridine_0p2eq", "run 3 with low pyridine loading", 5, 0.02, 9.2, 0.05,
    pyridine={"volume": ("0.15 mL of 0.031 M solution in toluene", 0.15, "mL"),
               "amount": ("0.005 mmol", 0.005, "mmol")},
    carrier_toluene=0.15, main_toluene=0.05, low_pyridine=True)
add_nbenzyl_control(
    "nbenzyl_pyridine_low_run4", "run 4 with low pyridine loading", 20, 0.09, 37, 0.19,
    pyridine={"volume": ("0.38 mL of 0.012 M solution in toluene", 0.38, "mL"),
               "amount": ("0.005 mmol", 0.005, "mmol")},
    carrier_toluene=0.38, main_toluene=0.41,
    observation="2a and byproduct 5a were not observed")


# Page 14: product-interconversion preparation, probe, and in situ NMR.
reaction, ev = base_reaction(
    "interconversion_mixture_preparation", [14],
    "Probing for Conversion of Allyl Sulfone 4a to Vinyl Sulfone 3a, Experiment 1a",
    "preparation of the 59:41 3a/4a mixture used in the interconversion probe",
    "C(sp3)-H functionalization in ethanol",
    ("mechanistic_study",),
    ("The paragraph says it follows General Procedure C but explicitly substitutes EtOH for THF; EtOH is retained as the experiment-specific solvent.",
     "General Procedure C supplies 50 °C and monitoring to disappearance but no fixed duration; duration remains unreported."))
participant(reaction, "m_1a", "reactant", ev, mass=("80 mg", 80, "mg"),
            amount=("0.38 mmol", 0.38, "mmol"), equivalents=1, reference="m_1a")
participant(reaction, "m_2a", "reactant", ev, mass=("146 mg", 146, "mg"),
            amount=("0.75 mmol", 0.75, "mmol"))
participant(reaction, "m_ethanol", "solvent", ev, volume=("2 mL", 2, "mL"))
participant(reaction, "m_3a", "product", ev, product_role="major product")
participant(reaction, "m_4a", "product", ev, product_role="minor product")
set_scale(reaction, ev, "m_1a", "0.38 mmol 1a", 0.38, "mmol")
set_solvent(reaction, ev, ["m_ethanol"], "EtOH", volume=("2 mL", 2, "mL"),
            concentration=("0.19 M", 0.19, "M"), basis="0.38 mmol 1a in 2 mL EtOH")
step(reaction, tids(6, 14),
     "perform General Procedure C with experiment-specific EtOH replacing THF: charge in air, seal, evacuate/fill with nitrogen three times, and monitor at 50 °C to disappearance",
     temperature=("50 °C", Q(50, "°C")), atmosphere="nitrogen",
     method="inherited_general_condition")
add_standard_workup(reaction, ev, "30-70% ethyl acetate/hexane")
outcome(reaction, ev, "product_amount", ["m_3a", "m_4a"], "106.5 mg", Q(106.5, "mg"),
        role="combined 3a/4a isolated mixture",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "yield", ["m_3a", "m_4a"], "70%", Q(70, "%"),
        role="combined 3a/4a mixture", yield_type="isolated yield",
        basis="isolated after silica-gel flash chromatography")
outcome(reaction, ev, "selectivity", ["m_3a", "m_4a"], "3a/4a = 59:41",
        {"ratio": [59, 41], "order": ["m_3a", "m_4a"], "metric": "product ratio"},
        basis="post-chromatography 1H NMR integration ratios", scope="3a/4a")
classify_cc_ch(reaction, ev)
save_reaction(reaction)

reaction, ev = base_reaction(
    "interconversion_probe", [4, 14],
    "Probing for Conversion of Allyl Sulfone 4a to Vinyl Sulfone 3a, Experiment 1b",
    "probe for conversion of allyl sulfone 4a to vinyl sulfone 3a",
    "product-interconversion control with concurrent coupling",
    ("mechanistic_study", "control"),
    ("The final ratio contains both newly formed material and the charged 59:41 3a/4a mixture; it is not a direct product-only interconversion yield.",))
participant(reaction, "m_1a", "reactant", ev, mass=("20 mg", 20, "mg"),
            amount=("0.09 mmol", 0.09, "mmol"))
participant(reaction, "m_2a", "reactant", ev, mass=("37 mg", 37, "mg"),
            amount=("0.19 mmol", 0.19, "mmol"))
participant(reaction, "m_3a", "product", ev, mass=("part of 38.3 mg charged 3a/4a mixture", 38.3, "mg"),
            product_role="charged and observed product")
participant(reaction, "m_4a", "product", ev,
            product_role="charged and observed product")
participant(reaction, "m_thf", "solvent", ev, volume=("0.5 mL", 0.5, "mL"))
set_scale(reaction, ev, "m_1a", "0.09 mmol 1a", 0.09, "mmol")
set_solvent(reaction, ev, ["m_thf"], "THF", volume=("0.5 mL", 0.5, "mL"))
step(reaction, ev, "charge 1a, 2a, and 38.3 mg (0.09 mmol) of a 59:41 3a/4a mixture in air and seal",
     atmosphere="air")
step(reaction, ev, "evacuate and fill with nitrogen three times", atmosphere="nitrogen")
step(reaction, ev, "add THF, lower into a preheated oil bath, and react",
     temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")),
     atmosphere="nitrogen")
add_standard_workup(reaction, ev)
outcome(reaction, ev, "selectivity", ["m_3a", "m_4a"], "charged 3a/4a = 59:41",
        {"ratio": [59, 41], "order": ["m_3a", "m_4a"], "metric": "charged mixture ratio"},
        basis="post-chromatography 1H NMR ratio of the separately prepared charged mixture",
        scope="initial charged mixture")
outcome(reaction, ev, "selectivity", ["m_3a", "m_4a"], "crude after 24 h: 3a/4a = 78:22",
        {"ratio": [78, 22], "order": ["m_3a", "m_4a"], "metric": "product ratio"},
        basis="crude 1H NMR", scope="whole mixture after 24 h")
classify_cc_ch(reaction, ev)
save_reaction(reaction)

reaction, ev = base_reaction(
    "in_situ_nmr_monitoring", [4, 14],
    "Probing for Conversion, Experiment 2", "in situ 1H NMR monitoring of 1a with 2a",
    "reaction time-course monitoring", ("mechanistic_study",),
    ("No quantitative yield is reported.",
     "4a was not observed at the experiment's sensitivity; this is not converted into a zero yield."))
participant(reaction, "m_1a", "reactant", ev, mass=("5 mg", 5, "mg"),
            amount=("0.02 mmol", 0.02, "mmol"))
participant(reaction, "m_2a", "reactant", ev, mass=("9.2 mg", 9.2, "mg"),
            amount=("0.05 mmol", 0.05, "mmol"))
participant(reaction, "m_toluene_d8", "solvent", ev, volume=("0.5 mL", 0.5, "mL"))
participant(reaction, "m_3a", "product", ev, product_role="observed product")
participant(reaction, "m_4a", "product", ev, product_role="not observed at low abundance")
set_scale(reaction, ev, "m_1a", "0.02 mmol 1a", 0.02, "mmol")
set_solvent(reaction, ev, ["m_toluene_d8"], "toluene-d8", volume=("0.5 mL", 0.5, "mL"))
step(reaction, ev, "charge 1a and 2a in air, add toluene-d8, transfer to an NMR tube, and seal",
     atmosphere="air")
step(reaction, ev, "evacuate and fill with nitrogen three times", atmosphere="nitrogen")
step(reaction, ev, "insert into a 600 MHz NMR probe preheated to 50 °C and collect 24 spectra",
     temperature=("50 °C", Q(50, "°C")), time=("25.5 h", Q(25.5, "h")),
     atmosphere="nitrogen")
outcome(reaction, ev, "other", ["m_1a", "m_3a", "m_4a"],
        "3a increased as 1a decreased; 4a signal was not observed at any time because of low abundance",
        "qualitative time-course observation", basis="24 in situ 1H NMR spectra",
        scope="non-observation of 4a is not a quantified zero")
classify_cc_ch(reaction, ev)
save_reaction(reaction)


# Page 14: two K2CO3/5b negative controls.
def add_5b_control(tag, experiment_no, k_mass, k_amount, observation):
    pages = [3, 14] if experiment_no == 1 else [14]
    unresolved = [
        "No coupling product is reported; no product is invented.",
        "5b loss/decomposition is qualitative and is not assigned a concrete product structure."]
    if experiment_no == 1:
        unresolved.append(
            "Page 3 summarizes recovered starting materials, while the detailed page-14 paragraph reports unreacted 1a and decomposition of 5b; both statements are retained as a source conflict.")
    reaction, ev = base_reaction(
        tag, pages, f"Reaction of Picolyl Amide 1a with Disulfone 5b, Experiment {experiment_no}",
        f"1a/5b/K2CO3 negative control, experiment {experiment_no}",
        "negative mechanistic control", ("mechanistic_study", "control"), tuple(unresolved))
    participant(reaction, "m_1a", "reactant", ev, mass=("10 mg", 10, "mg"),
                amount=("0.05 mmol", 0.05, "mmol"))
    participant(reaction, "m_5b", "reactant", ev, mass=("16 mg", 16, "mg"),
                amount=("0.05 mmol", 0.05, "mmol"))
    participant(reaction, "m_k2co3", "reagent", ev,
                mass=(f"{k_mass:g} mg", k_mass, "mg"),
                amount=(f"{k_amount:g} mmol", k_amount, "mmol"))
    participant(reaction, "m_acetonitrile", "solvent", ev,
                volume=("0.5 mL", 0.5, "mL"))
    set_scale(reaction, ev, "m_1a", "0.05 mmol 1a", 0.05, "mmol")
    set_solvent(reaction, ev, ["m_acetonitrile"], "acetonitrile",
                volume=("0.5 mL", 0.5, "mL"), concentration=("0.1 M", 0.1, "M"),
                basis="0.05 mmol 1a in 0.5 mL acetonitrile")
    step(reaction, ev, "charge 1a, 5b, and potassium carbonate in air and seal",
         atmosphere="air")
    step(reaction, ev, "evacuate and fill with nitrogen three times", atmosphere="nitrogen")
    step(reaction, ev, "add acetonitrile, lower into a preheated oil bath, and react",
         temperature=("50 °C", Q(50, "°C")), time=("24 h", Q(24, "h")),
         atmosphere="nitrogen")
    add_standard_workup(reaction, ev)
    outcome(reaction, ev, "recovery", ["m_1a"], "unreacted 1a",
            "unreacted starting material observed", role="starting material",
            basis="crude 1H NMR")
    outcome(reaction, ev, "decomposition" if experiment_no == 1 else "other",
            ["m_5b"], observation, observation, role="reactant fate",
            basis="crude 1H NMR",
            scope="qualitative observation; no identified product")
    classify_none(reaction, ev, "No coupling product formed; the fate of 5b is not structurally assigned.")
    if experiment_no == 1:
        reaction["notes"].append(F(
            "Source conflict retained: page 3 says recovered starting materials; page 14 says unreacted 1a and decomposition of 5b.",
            {"page3": "recovered starting materials",
             "page14": "unreacted 1a and decomposition of 5b"}, ev))
    save_reaction(reaction)


add_5b_control("control_1a_5b_k2co3", 1, 13, 0.09, "decomposition of 5b")
add_5b_control("control_1a_5b_k2co3_exp2", 2, 6.5, 0.05,
               "almost complete disappearance of 5b")


# Characterization and reference material are retained as source context, not
# synthetic reactions.
context_ids = {
    f"{DOC}:context_pages13_14_characterization": {
        "context_id": f"{DOC}:context_pages13_14_characterization",
        "document_id": DOC,
        "evidence_ids": tids(13, 14),
        "kind": "characterization_or_analysis",
        "description": "NMR, 2H NMR, 13C NMR, HRMS, TLC, and distinguishable minor-product resonances on pages 13-14 characterize 3l, 3m/4m, and the labeled 3a-d/4a-d mixtures; spectral peaks are not standalone syntheses.",
    },
    f"{DOC}:context_page14_associated_content": {
        "context_id": f"{DOC}:context_page14_associated_content",
        "document_id": DOC,
        "evidence_ids": tids(14),
        "kind": "supplementary_material_notice",
        "description": "The article advertises external Supporting Information and crystallographic files, but no SI document is explicitly supplied or associated in the workspace mapping; no SI facts were imported.",
    },
    f"{DOC}:context_pages15_17_back_matter": {
        "context_id": f"{DOC}:context_pages15_17_back_matter",
        "document_id": DOC,
        "evidence_ids": tids(15, 16, 17),
        "kind": "back_matter_and_cited_literature",
        "description": "Pages 15-17 contain acknowledgments, additional spectral-note markers, and references to prior work. Cited literature transformations are not current-paper experiments and no new reaction record is created.",
    },
}
bundle["source_context"] = [
    c for c in bundle["source_context"] if c["context_id"] not in context_ids]
bundle["source_context"].extend(context_ids.values())


coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {"region": "top-left continuation of 3l and 3m/4m scope paragraph",
                 "finding": "The 3l continuation is characterization only. Two independent 3m/4m runs were recorded with their charges, isolated mixture masses, yields, and NMR product ratios.",
                 "evidence_ids": tids(13)},
                {"region": "middle Radical Inhibitor, Radical Initiator, and Crossover experiment paragraphs",
                 "finding": "Existing Scheme 2a-c records were extended with concrete charges, ordered operations, isolated masses/yields, selectivity methods, recovered 1a, crossover byproduct 5a, and recovered 5b/2a.",
                 "evidence_ids": tids(13)},
                {"region": "lower Deuterium-Labeling Experiment paragraphs",
                 "finding": "The detailed first labeling run was merged into the existing Scheme 2d record and the second run was created separately. The source-label inconsistency 2d versus 2a-d2 was preserved.",
                 "evidence_ids": tids(13)},
            ],
            "record_ids": [
                f"{DOC}:scope_3m_run1", f"{DOC}:scope_3m_run2",
                f"{DOC}:scheme2a_bht", f"{DOC}:scheme2b_aibn",
                f"{DOC}:scheme2c_crossover", f"{DOC}:scheme2d_deuterium",
                f"{DOC}:scheme2d_deuterium_run2"],
            "unresolved": [
                "Product-specific masses and absolute configurations for the 3/4 mixtures are unreported.",
                "The labeled-allene paragraph calls the reagent 2d, whereas its preparation names it 2a-d2; no silent relabeling was made."],
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of deuterium run 2",
                 "finding": "Only isotopologue characterization continues here; it supports the run-2 record and is not a separate synthesis.",
                 "evidence_ids": tids(14)},
                {"region": "N-benzylbenzamide/2a runs 1-4",
                 "finding": "Four control conditions were kept distinct. Runs 2 and 3 refine existing page-5 summaries; runs 1 and 4 are new records. The 0.2-equivalent summary versus 0.25-equivalent charge-derived pyridine conflict is explicit.",
                 "evidence_ids": tids(14)},
                {"region": "product-interconversion Experiments 1a, 1b, and in situ NMR Experiment 2",
                 "finding": "The EtOH preparation of the 59:41 mixture, the charged-mixture probe ending at 78:22, and the 25.5 h/24-spectrum NMR time course are three distinct records.",
                 "evidence_ids": tids(14)},
                {"region": "1a/5b/K2CO3 Experiments 1-2 and Associated Content",
                 "finding": "Both negative controls were recorded without invented products. The page-3/page-14 starting-material conflict is retained. The SI notice is source context because no SI file is associated.",
                 "evidence_ids": tids(14)},
            ],
            "record_ids": [
                f"{DOC}:scheme2d_deuterium_run2",
                f"{DOC}:nbenzyl_no_pyridine_run1",
                f"{DOC}:nbenzyl_pyridine_1eq",
                f"{DOC}:nbenzyl_pyridine_0p2eq",
                f"{DOC}:nbenzyl_pyridine_low_run4",
                f"{DOC}:interconversion_mixture_preparation",
                f"{DOC}:interconversion_probe",
                f"{DOC}:in_situ_nmr_monitoring",
                f"{DOC}:control_1a_5b_k2co3",
                f"{DOC}:control_1a_5b_k2co3_exp2"],
            "unresolved": [
                "Qualitative non-observation of 2a, 4a, or 5a was not converted into a quantitative zero.",
                "The SI is advertised but not supplied or associated; no SI-only detail was inferred."],
        },
        {
            "page": 15,
            "state": "attempted",
            "regions": [
                {"region": "acknowledgments, Additional Notes, and references 1-9",
                 "finding": "No current-paper synthesis is reported. Superscripts a/b explain distinguishable isomer/diastereomer signals, and the remaining material is acknowledgments or cited prior literature.",
                 "evidence_ids": tids(15)}],
            "record_ids": [],
            "unresolved": ["Cited prior-work transformations were not re-extracted as experiments in this paper."],
        },
        {
            "page": 16,
            "state": "attempted",
            "regions": [
                {"region": "references 9-20",
                 "finding": "The full page is bibliography. It contains cited prior work but no current-paper experimental arrow, table row, or procedure.",
                 "evidence_ids": tids(16)}],
            "record_ids": [],
            "unresolved": ["No current-paper reaction record is warranted for bibliographic entries."],
        },
        {
            "page": 17,
            "state": "attempted",
            "regions": [
                {"region": "references 21-48 and reference footnotes",
                 "finding": "The page is bibliography and explanatory citation text. The K2CO3 citation note is rationale for a page-14 control, not a new experiment.",
                 "evidence_ids": tids(17)}],
            "record_ids": [],
            "unresolved": ["No current-paper reaction record is warranted for cited literature or citation rationale."],
        },
    ],
}


page_inventory = {
    "document_id": DOC,
    "pages": [13, 14, 15, 16, 17],
    "source_title": TRUE_TITLE,
    "presentation_conversion": False,
    "new_reaction_ids": ADDED_THIS_PASS,
    "updated_reaction_ids": UPDATED_THIS_PASS,
    "notes": [
        "All five original source-PDF page images were inspected at original detail.",
        "Pages 15-17 contain no current-paper experiments.",
        "No Supporting Information file was associated or used."],
}

stereo_checks = {}
for mid in ("m_1m", "m_2a", "m_3m", "m_4m", "m_1a", "m_3a", "m_4a",
            "m_3q", "m_5a", "m_5b", "m_2a_d2", "m_11"):
    smiles = molecules[mid]["smiles"]["normalized_value"]
    if smiles:
        stereo_checks[mid] = stereo_inspect(smiles)

dataset = json.loads((RUN / "dataset.json").read_text())
mentions = source_mentions(
    dataset,
    ["3m", "4m", "Radical Inhibitor", "Radical Initiator",
     "Crossover Experiment", "Deuterium-Labeling Experiment",
     "N-Benzylbenzamide", "Probing for Conversion",
     "Reaction of Picolyl Amide 1a with Disulfone 5b"], context=260)

(RAW / "page_inventory_0013_0017.json").write_text(
    json.dumps(page_inventory, ensure_ascii=False, indent=2) + "\n")
(RAW / "source_mentions_pages_0013_0017.json").write_text(
    json.dumps(mentions, ensure_ascii=False, indent=2) + "\n")
(RAW / "stereo_checks_pages_0013_0017.json").write_text(
    json.dumps(stereo_checks, ensure_ascii=False, indent=2) + "\n")
(RUN / "coverage/pages_0013_0017.json").write_text(
    json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
OUT.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "molecules": len(bundle["molecules"]),
    "reactions": len(bundle["reactions"]),
    "new_reactions": len(ADDED_THIS_PASS),
    "updated_reactions": len(UPDATED_THIS_PASS),
    "coverage_pages": [p["page"] for p in coverage["pages"]],
}, ensure_ascii=False))
