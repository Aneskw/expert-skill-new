#!/usr/bin/env python3
"""Build evidence-linked v3.1 records for original PDF pages 7--12."""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as sr

RUN = ROOT / "runs/current"
DOC = "doc_0341d857f4f70e0f"
DATASET = RUN / "dataset.json"
OLD_STANDARD = RUN / "standard_records.json"
SUBMITTED = RUN / "submitted.json"
RAW = RUN / "raw_model/pages_0007_0012_extraction.json"


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads(DATASET.read_text())
old = json.loads(OLD_STANDARD.read_text())
bundle = sr.empty(dataset["documents"])
molecules = {m["molecule_id"]: copy.deepcopy(m) for m in old["molecules"]}
old_reactions = {r["reaction_id"]: r for r in old["reactions"]}
diagnostics = {}


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def f(raw=None, normalized=None, evidence=(), method="original_extraction",
      derivation=None, status=None):
    if normalized is None and raw is not None:
        normalized = raw
    return sr.fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                   derivation=derivation, status=status)


def mid(key: str) -> str:
    return f"{DOC}:mol:{key}"


def add_molecule(key: str, label: str, page: int, name: str | None = None,
                 smiles: str | None = None, chirality="unresolved",
                 geometry="unresolved", representation="unresolved",
                 identity_note: str | None = None):
    molecule_id = mid(key)
    if molecule_id in molecules:
        if smiles:
            diagnostics[molecule_id] = sr.stereo_inspect(smiles)
        return molecule_id
    evidence = ev(page)
    m = sr.new_molecule(molecule_id, DOC, evidence)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = f(label, label, evidence)
    m["names"]["original"] = f(name or label, name or label, evidence)
    if identity_note:
        m["identity_note"] = f(identity_note, identity_note, evidence)
    if smiles:
        diag = sr.stereo_inspect(smiles)
        diagnostics[molecule_id] = diag
        m["representation_type"] = "smiles"
        m["smiles"] = f(
            {"source_label": label, "machine_candidate": smiles},
            diag["canonical_isomeric_smiles"], evidence, "structure_derivation",
            "Canonicalized from a machine structure candidate; parsing does not establish source correctness."
        )
        m["stereochemistry"]["isomeric_smiles"] = f(
            diag["canonical_isomeric_smiles"], diag["canonical_isomeric_smiles"], evidence,
            "structure_derivation", "RDKit canonical isomeric SMILES from the machine candidate."
        )
        m["stereochemistry"]["identity_key"] = f(
            diag["stereo_graph_key"], diag["stereo_graph_key"], evidence,
            "structure_derivation", "Stereo graph key calculated from the machine candidate."
        )
    reason = (
        "Assigned after inspection of the page label/name and any shown stereochemical notation; "
        "unreported configuration was not converted to achiral or racemic."
    )
    for field, value in (("chirality_status", chirality), ("geometry_status", geometry),
                         ("representation_status", representation)):
        m["stereochemistry"][field] = f(None, value, evidence, "inference", reason)
    molecules[molecule_id] = m
    return molecule_id


SIMPLE = {
    "acetone": ("acetone", "CC(C)=O"),
    "propargyl_bromide": ("3-bromopropyne", "C#CCBr"),
    "ethyl_acetate": ("ethyl acetate", "CCOC(C)=O"),
    "sodium_sulfate": ("anhydrous sodium sulfate", "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    "methyl_iodide": ("methyl iodide", "CI"),
    "benzyl_bromide": ("benzyl bromide", "BrCc1ccccc1"),
    "allyl_bromide": ("allyl bromide", "C=CCBr"),
    "DCM": ("dichloromethane", "ClCCl"),
    "NaHCO3": ("sodium bicarbonate", "[Na+].O=C([O-])O"),
}
for key, (name, smiles) in SIMPLE.items():
    add_molecule(key, name, 7 if key in {"acetone", "propargyl_bromide", "ethyl_acetate", "sodium_sulfate"} else 8,
                 name, smiles, "achiral", "not_applicable", "not_applicable")
for key, name, page in [
    ("petroleum_ether", "petroleum ether (60–80 °C fraction)", 7),
    ("brine", "saturated salt water / brine", 7),
    ("silica_gel", "silica gel (100–200 mesh)", 7),
    ("unknown_butylating_reagent", "butylating reagent not listed in the general-procedure alternatives", 8),
]:
    add_molecule(key, name, page, name, None, "not_applicable", "not_applicable", "unresolved",
                 "Mixture/material or exact reagent identity is not fully specified on inspected pages.")

# Generic source labels and product-specific unresolved precursors.
for key, label, page in [
    ("S1_general", "S1", 7), ("S2_general", "S2", 7),
    ("S3_general", "S3", 8), ("6_general", "6", 7),
    ("7_general", "7 / ArI", 7), ("1_general", "1", 7),
    ("2_general", "2", 12), ("9_general", "9", 12),
]:
    add_molecule(key, label, page, f"generic source-labelled compound {label}", None,
                 "unresolved", "unresolved", "unresolved",
                 "Generic R/Ar-labelled structure; no concrete SMILES assigned.")
for suffix in "bcde":
    ch = "unspecified" if suffix in "cde" else "achiral"
    add_molecule(f"S1{suffix}", f"S1{suffix} (corresponding precursor)", 7,
                 f"S1 precursor corresponding to S2{suffix}", None, ch,
                 "not_applicable", "not_encoded" if ch == "unspecified" else "not_applicable",
                 "The individual S1 label/structure is not shown on pages 7–12; identity is retained only as the corresponding unresolved precursor.")
for suffix, page, ch in [("a",7,"achiral"),("b",7,"achiral"),("c",7,"unspecified"),("d",7,"unspecified"),("e",8,"unspecified")]:
    add_molecule(f"S2{suffix}", f"S2{suffix}", page, f"acetylene S2{suffix}", None, ch,
                 "not_applicable", "not_encoded" if ch == "unspecified" else "not_applicable")
for suffix, page, ch in [
    ("a",8,"achiral"),("b",8,"achiral"),("c",8,"achiral"),("d",8,"achiral"),("e",8,"unspecified"),
    ("f",8,"unspecified"),("g",8,"unspecified"),("i",8,"achiral")
]:
    add_molecule(f"S3{suffix}", f"S3{suffix}", page, f"indole-N-substituted acetylene S3{suffix}", None,
                 ch, "not_applicable", "not_encoded" if ch == "unspecified" else "not_applicable")


def participant(r, molecule_id, role, pid, evidence, charges=None, reference=None,
                stage="reaction", method="original_extraction"):
    existing = next((p for p in r["participants"] if p["participant_id"] == pid), None)
    if existing is not None:
        if existing["molecule_id"] != molecule_id:
            raise ValueError(f"participant {pid} already refers to {existing['molecule_id']}, not {molecule_id}")
        p = existing
        p["role"] = role
        p["stage"] = stage
    else:
        p = sr.new_participant(pid, molecule_id, role, evidence, stage)
    for field, (raw, norm) in (charges or {}).items():
        p["charges"][field] = f(raw, norm, evidence, method)
    if reference:
        p["equivalence_reference"] = f(reference, reference, evidence, method)
    if existing is None:
        r["participants"].append(p)
    return p


def step(r, step_id, page, stage, pids, operation, temperature=None, time=None,
         atmosphere=None, method="original_extraction"):
    s = sr.new_step(step_id, ev(page), stage)
    s["participant_ids"] = list(pids)
    s["operation"] = f(operation, operation, ev(page), method)
    if temperature:
        s["temperature"] = f(temperature[0], temperature[1], ev(page), method)
    if time:
        s["time"] = f(time[0], time[1], ev(page), method)
    if atmosphere:
        s["atmosphere"] = f(atmosphere, atmosphere, ev(page), method)
    r["steps"].append(s)
    return s


def outcome(r, oid, kind, target, raw, norm, page, basis=None, yield_type=None, scope=None):
    o = sr.new_outcome(oid, kind, [target], ev(page))
    o["product_role"] = f("product", "product", ev(page)) if kind in {"yield", "product_amount"} else sr.fact(evidence=ev(page))
    o["value"] = f(raw, norm, ev(page))
    if yield_type:
        o["yield_type"] = f(yield_type, yield_type, ev(page))
    if basis:
        o["measurement_basis"] = f(basis, basis, ev(page))
    if scope:
        o["scope"] = f(scope, scope, ev(page))
    existing = next((i for i, item in enumerate(r["outcomes"]) if item["outcome_id"] == oid), None)
    if existing is None:
        r["outcomes"].append(o)
    else:
        r["outcomes"][existing] = o
    return o


def reaction(suffix, page, record_type, name, categories=None, extra_pages=()):
    evidence = ev(page)
    for p in extra_pages:
        evidence += ev(p)
    r = sr.new_reaction(f"{DOC}:{suffix}", DOC, evidence)
    r["record_type"] = record_type
    r["reaction_name"] = f(name, name, evidence)
    r["location"] = f(
        "PDF page " + ", ".join(map(str, (page,) + tuple(extra_pages))),
        {"pdf_pages": [page, *extra_pages]}, evidence
    )
    if categories is not None:
        r["experiment_categories"] = f(categories, categories, evidence)
    r["review"]["unresolved"] = [
        "Machine extraction remains needs_review; schema validity does not establish chemical correctness."
    ]
    return r


def classify(r, page, mode=None, note=None):
    evidence = ev(page)
    if mode is None:
        r["classification"]["labels"] = f([], [], evidence, "inference",
            note or "No single net experimental transformation is asserted for this generic/mechanistic record.")
        r["classification"]["coverage"] = f(None, "unresolved", evidence, "inference",
            note or "No atom mapping; classification unresolved.")
        r["classification"]["limitations"] = [note or "No atom mapping; classification unresolved."]
        return
    spec = {
        "CN": ("C", "N", 0, 1, "formed", 5, "C-N"),
        "CC": ("C", "C", 0, 1, "formed", 2, "C-C"),
        "DBCC": ("C", "C", 1, 2, "order_increased", 13, "DB-C-C"),
    }
    a, b, before, after, change_type, number, code = spec[mode]
    change = {"change_id":"b1", "atoms":[{"atom_map":None,"element":a},{"atom_map":None,"element":b}],
              "before_order":before, "after_order":after, "change_type":change_type}
    r["classification"]["bond_changes"] = [f(change, change, evidence, "structure_derivation",
        "Partial net bond change inferred from the named/drawn source transformation without atom mapping.")]
    labels = [{"category_number":number,"code":code,"supporting_change_ids":["b1"]}]
    r["classification"]["labels"] = f(labels, labels, evidence, "structure_derivation",
        "Taxonomy label inferred from the net source transformation; mapping remains incomplete.")
    r["classification"]["coverage"] = f(None, "partial", evidence, "inference",
        "Only the clearly indicated formed/promoted bond is classified; no atom-mapped reaction SMILES was supplied.")
    r["classification"]["limitations"] = [
        "No atom mapping; leaving-group, implicit-H, and other bond-order changes may be incomplete."
    ]


def set_procedure(r, page, text, addition=None, pretreatment=None, method="original_extraction"):
    r["procedure"]["original_text"] = f(text, text, ev(page), method)
    if addition:
        r["procedure"]["addition_order"] = f(addition, addition, ev(page), method)
    if pretreatment:
        r["procedure"]["pretreatment"] = f(pretreatment, pretreatment, ev(page), method)


def workup_participants(r, page, prefix="", extraction_volume=None, purification="25–35% ethyl acetate–petroleum ether (v/v)"):
    evidence = ev(page)
    water = mid("water")
    et = mid("ethyl_acetate")
    brine = mid("brine")
    drying = mid("sodium_sulfate")
    silica = mid("silica_gel")
    pet = mid("petroleum_ether")
    if not any(p["participant_id"] == prefix+"p_water" for p in r["participants"]):
        participant(r, water, "workup_agent", prefix+"p_water", evidence, stage="workup")
        charges = {"volume": extraction_volume} if extraction_volume else None
        participant(r, et, "workup_agent", prefix+"p_EtOAc", evidence, charges=charges, stage="workup")
        participant(r, brine, "workup_agent", prefix+"p_brine", evidence, stage="workup")
        participant(r, drying, "workup_agent", prefix+"p_Na2SO4", evidence, stage="workup")
        participant(r, silica, "purification_agent", prefix+"p_silica", evidence, stage="purification")
        participant(r, pet, "purification_agent", prefix+"p_petroleum_ether", evidence, stage="purification")
    wp = [p["participant_id"] for p in r["participants"] if p["stage"] == "workup"]
    pp = [p["participant_id"] for p in r["participants"] if p["stage"] == "purification"]
    step(r, prefix+"workup", page, "workup", wp, "quenched/extracted, washed, dried, and concentrated as reported")
    step(r, prefix+"purification", page, "purification", pp,
         f"silica-gel column chromatography; eluent {purification}")


# Scheme 6 is a proposed mechanism, not an experimental arrow.
mech = reaction("scheme6_proposed_mechanism", 7, "proposed_mechanism",
                "Scheme 6 proposed Pd-catalyzed cyclocondensation/migration mechanism")
for key, role, pid in [("6_general","reactant","p_6"),("7_general","reactant","p_7"),("1_general","product","p_1")]:
    participant(mech, mid(key), role, pid, ev(7))
mech["procedure"]["original_text"] = f(
    "Proposed sequence A–G: oxidative addition, allene addition, intramolecular indole attack, reductive elimination, allylic migration, and deprotonation.",
    "Proposed mechanism only", ev(7))
mech["notes"].append(f("Intermediates A–G and elementary arrows are hypotheses, not separate reported experiments.",
                        "Mechanistic proposal, not experimental records", ev(7)))
classify(mech, 7, None, "Scheme 6 is a proposed multistep mechanism rather than one atom-mapped net experiment.")
bundle["reactions"].append(mech)

# Explicit Ag2CO3 mechanistic control.
control = reaction("control_9a_ag2co3_to_1ha", 7, "specific_experiment",
                   "Ag2CO3-promoted conversion of spiroindolenine 9a to THBC 1ha",
                   ["control", "mechanistic_study"])
participant(control, mid("9a"), "reactant", "p_9a", ev(7), {"amount":("1 mmol",{"number":1,"unit":"mmol"})})
participant(control, mid("Ag2CO3"), "reagent", "p_Ag2CO3", ev(7), {"amount":("1.1 mmol",{"number":1.1,"unit":"mmol"})})
participant(control, mid("DMF"), "solvent", "p_DMF", ev(7), {"volume":("2.5 mL",{"number":2.5,"unit":"mL"})})
participant(control, mid("1ha"), "product", "p_1ha", ev(7))
control["reaction_scale"] = f("9a (1 mmol)", {"reference_molecule_id":mid("9a"),"amount":{"number":1,"unit":"mmol"}}, ev(7))
step(control, "heating", 7, "reaction", ["p_9a","p_Ag2CO3","p_DMF","p_1ha"],
     "heated in DMF with Ag2CO3", None, ("24 h",{"number":24,"unit":"h"}))
set_procedure(control, 7, "Spiroindolenine 9a (1 mmol) was heated in DMF (2.5 mL) with Ag2CO3 (1.1 mmol); 1ha was isolated in 40% yield after 24 h.")
outcome(control, "yield_1ha", "yield", mid("1ha"), "40%", {"number":40,"unit":"%"}, 7,
        "isolated product", "isolated")
classify(control, 7, "CC")
control["review"]["unresolved"].append("Heating temperature and workup/purification are not reported in the page-7 control description.")
bundle["reactions"].append(control)


def add_s2_participants(r, s1, s2, method="original_extraction"):
    evidence = ev(7)
    participant(r, s1, "reactant", "p_S1", evidence, {"amount":("10 mmol",{"number":10,"unit":"mmol"})}, method=method)
    participant(r, mid("K2CO3"), "reagent", "p_K2CO3", evidence,
                {"mass":("2.8 g",{"number":2.8,"unit":"g"}),"amount":("20 mmol",{"number":20,"unit":"mmol"})}, method=method)
    participant(r, mid("propargyl_bromide"), "reagent", "p_propargyl_bromide", evidence,
                {"volume":("1.3 mL",{"number":1.3,"unit":"mL"}),"amount":("15 mmol",{"number":15,"unit":"mmol"})}, method=method)
    participant(r, mid("acetone"), "solvent", "p_acetone", evidence,
                {"volume":("25 mL",{"number":25,"unit":"mL"})}, method=method)
    participant(r, s2, "product", "p_S2", evidence, method=method)


def add_s2_steps(r, method="original_extraction"):
    step(r, "addition", 7, "reaction", ["p_S1","p_K2CO3","p_propargyl_bromide","p_acetone"],
         "K2CO3 was added to S1 in acetone under nitrogen, followed by 3-bromopropyne", atmosphere="nitrogen", method=method)
    step(r, "reflux", 7, "reaction", ["p_S1","p_K2CO3","p_propargyl_bromide","p_acetone","p_S2"],
         "heated to reflux until TLC completion", ("reflux", "reflux"),
         ("24–48 h",{"minimum":24,"maximum":48,"unit":"h"}), "nitrogen", method)
    workup_participants(r, 7)


s2gen = reaction("general_synthesis_S2", 7, "generic_scheme", "General procedure for acetylenes S2")
add_s2_participants(s2gen, mid("S1_general"), mid("S2_general"))
add_s2_steps(s2gen)
set_procedure(s2gen, 7, "S1 (10 mmol) in acetone (25 mL), K2CO3 (2.8 g, 20 mmol), then 3-bromopropyne (1.3 mL, 15 mmol); reflux 24–48 h; aqueous workup and silica chromatography.",
              "K2CO3, then 3-bromopropyne")
s2gen["notes"].append(f("General isolated-yield range for S2a–e: 50–98%; not assigned as an entry-specific yield.",
                         "General range only", ev(7)))
classify(s2gen, 7, "CN")
bundle["reactions"].append(s2gen)

s2_data = {
    "b": (7, "3.6 g", 95), "c": (7, "2.6 g", 60),
    "d": (7, "3.1 g", 68), "e": (8, "2.3 g", 50),
}
for suffix, (page, mass, yld) in s2_data.items():
    r = reaction(f"synthesis_S2{suffix}", page, "specific_experiment", f"Synthesis of acetylene S2{suffix}", ["other"], extra_pages=(7,) if page != 7 else ())
    r["general_scheme_ids"] = [s2gen["reaction_id"]]
    add_s2_participants(r, mid(f"S1{suffix}"), mid(f"S2{suffix}"), "inherited_general_condition")
    add_s2_steps(r, "inherited_general_condition")
    r["reaction_scale"] = f("S1 (10 mmol)", {"reference_molecule_id":mid(f"S1{suffix}"),"amount":{"number":10,"unit":"mmol"}}, ev(7), "inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f("General Procedure for the Synthesis of Acetylenes S2 from S1", s2gen["reaction_id"], ev(7), "inherited_general_condition")
    outcome(r, f"amount_S2{suffix}", "product_amount", mid(f"S2{suffix}"), mass,
            {"number":float(mass.split()[0]),"unit":mass.split()[1]}, page, "isolated product")
    outcome(r, f"yield_S2{suffix}", "yield", mid(f"S2{suffix}"), f"{yld}%",
            {"number":yld,"unit":"%"}, page, "isolated product", "isolated")
    classify(r, 7, "CN")
    r["review"]["unresolved"].append("The corresponding S1 precursor structure/individual label is not shown on pages 7–12; no concrete precursor SMILES was assigned.")
    bundle["reactions"].append(r)


def s3_start_and_reagent(suffix):
    if suffix == "b": return mid("S2b"), mid("methyl_iodide")
    if suffix == "c": return mid("S2a"), mid("benzyl_bromide")
    if suffix == "d": return mid("S2a"), mid("allyl_bromide")
    if suffix == "e": return mid("S2c"), mid("methyl_iodide")
    if suffix == "f": return mid("S2d"), mid("methyl_iodide")
    if suffix == "g": return mid("S2e"), mid("methyl_iodide")
    return mid("S2a"), mid("unknown_butylating_reagent")


def add_s3_participants(r, start, product, reagent, method="original_extraction"):
    evidence = ev(8)
    participant(r, start, "reactant", "p_S2", evidence, {"amount":("1 mmol",{"number":1,"unit":"mmol"})}, method=method)
    participant(r, mid("DMF"), "solvent", "p_DMF", evidence, {"volume":("10.0 mL",{"number":10.0,"unit":"mL"})}, method=method)
    participant(r, mid("NaH"), "reagent", "p_NaH", evidence, {"amount":("1.5 mmol",{"number":1.5,"unit":"mmol"})}, method=method)
    participant(r, reagent, "reagent", "p_alkylating_agent", evidence,
                {"amount":("1.3 mmol",{"number":1.3,"unit":"mmol"})} if reagent != mid("unknown_butylating_reagent") else None, method=method)
    participant(r, product, "product", "p_S3", evidence, method=method)


def add_s3_steps(r, method="original_extraction"):
    step(r, "cool_and_add", 8, "reaction", ["p_S2","p_DMF","p_NaH","p_alkylating_agent"],
         "cooled to 0 °C; NaH and alkylating agent were added successively",
         ("0 °C",{"number":0,"unit":"°C"}), method=method)
    step(r, "cold_hold", 8, "reaction", ["p_S2","p_DMF","p_NaH","p_alkylating_agent"],
         "held after additions before warming", ("0 °C",{"number":0,"unit":"°C"}),
         ("5 min",{"number":5,"unit":"min"}), method=method)
    step(r, "warm_to_rt", 8, "reaction", ["p_S2","p_DMF","p_NaH","p_alkylating_agent"],
         "allowed to reach room temperature gradually", ("room temperature","room temperature"), method=method)
    step(r, "rt_hold", 8, "reaction", ["p_S2","p_DMF","p_NaH","p_alkylating_agent","p_S3"],
         "stirred at room temperature until TLC completion", ("room temperature","room temperature"), method=method)
    # Page 8 gives exact quench/extraction volumes.
    participant(r, mid("water"), "workup_agent", "p_water", ev(8), {"volume":("2.0 mL",{"number":2.0,"unit":"mL"})}, stage="workup", method=method)
    participant(r, mid("ethyl_acetate"), "workup_agent", "p_EtOAc", ev(8),
                {"volume":("3 × 20 mL",{"number":20,"unit":"mL","count":3})}, stage="workup", method=method)
    participant(r, mid("brine"), "workup_agent", "p_brine", ev(8), {"volume":("10 mL",{"number":10,"unit":"mL"})}, stage="workup", method=method)
    participant(r, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", ev(8), stage="workup", method=method)
    participant(r, mid("silica_gel"), "purification_agent", "p_silica", ev(8), stage="purification", method=method)
    participant(r, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(8), stage="purification", method=method)
    step(r, "workup", 8, "workup", ["p_water","p_EtOAc","p_brine","p_Na2SO4"],
         "quenched with water, extracted with ethyl acetate, washed with brine, dried, and evaporated", method=method)
    step(r, "purification", 8, "purification", ["p_silica","p_petroleum_ether"],
         "silica-gel chromatography with petroleum ether:ethyl acetate = 4:1 (v/v)", method=method)


s3gen = reaction("general_synthesis_S3", 8, "generic_scheme", "General procedure for indole-N-substituted acetylenes S3")
add_s3_participants(s3gen, mid("S2_general"), mid("S3_general"), mid("methyl_iodide"))
add_s3_steps(s3gen)
set_procedure(s3gen, 8, "S2 (1 mmol) in DMF (10.0 mL) was cooled to 0 °C; NaH (60% oil suspension, 1.5 mmol) and MeI/benzyl bromide/allyl bromide (1.3 mmol) were added successively; after 5 min, warm gradually to room temperature and stir until TLC completion.",
              "NaH, then MeI/benzyl bromide/allyl bromide at 0 °C")
classify(s3gen, 8, "CN")
s3gen["review"]["unresolved"].append("The page lists three alternative alkylating agents but does not map them entry-by-entry; S3i requires an unlisted butylating reagent.")
bundle["reactions"].append(s3gen)

s3_data = {
    "b":("301.7 mg",76), "c":("433.2 mg",98), "d":("345.0 mg",88),
    "e":("331.5 mg",75), "f":("146.1 mg",30), "g":("236 mg",50), "i":("367.2 mg",90),
}
for suffix, (mass, yld) in s3_data.items():
    start, reagent = s3_start_and_reagent(suffix)
    r = reaction(f"synthesis_S3{suffix}", 8, "specific_experiment", f"Synthesis of acetylene S3{suffix}", ["other"])
    r["general_scheme_ids"] = [s3gen["reaction_id"]]
    add_s3_participants(r, start, mid(f"S3{suffix}"), reagent, "inherited_general_condition")
    add_s3_steps(r, "inherited_general_condition")
    r["reaction_scale"] = f("S2 (1 mmol)", {"reference_molecule_id":start,"amount":{"number":1,"unit":"mmol"}}, ev(8), "inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f("General Procedure for the Synthesis of Indole-N-substituted Acetylenes S3", s3gen["reaction_id"], ev(8), "inherited_general_condition")
    outcome(r, f"amount_S3{suffix}", "product_amount", mid(f"S3{suffix}"), mass,
            {"number":float(mass.split()[0]),"unit":"mg"}, 8, "isolated product")
    outcome(r, f"yield_S3{suffix}", "yield", mid(f"S3{suffix}"), f"{yld}%", {"number":yld,"unit":"%"}, 8,
            "isolated product", "isolated")
    classify(r, 8, "CN")
    if suffix == "i":
        r["notes"].append(f("S3i is N-butyl, but the listed general-procedure alternatives are MeI, benzyl bromide, and allyl bromide; no butylating reagent is named.",
                            "Source omission/conflict retained; no reagent guessed", ev(8)))
        r["review"]["unresolved"].append("Exact alkylating reagent and its amount for S3i are unreported in the listed alternatives; the 1.3 mmol alternative dose was not assigned to the unknown reagent.")
    bundle["reactions"].append(r)


def add_allene_participants(r, start, product, method="original_extraction"):
    evidence = ev(8) + ev(9)
    participant(r, start, "reactant", "p_S3", evidence, {"amount":("0.25 mmol",{"number":0.25,"unit":"mmol"})}, method=method)
    participant(r, mid("THF"), "solvent", "p_THF", evidence, {"volume":("3 mL",{"number":3,"unit":"mL"})}, method=method)
    participant(r, mid("KOtBu"), "reagent", "p_KOtBu", evidence,
                {"mass":("42 mg",{"number":42,"unit":"mg"}),"amount":("0.38 mmol",{"number":0.38,"unit":"mmol"})}, method=method)
    participant(r, product, "product", "p_6", evidence, method=method)


def add_allene_steps(r, method="original_extraction"):
    step(r, "cold_addition", 8, "reaction", ["p_S3","p_THF","p_KOtBu"],
         "KOtBu was added to S3 in dry THF under argon at ice-cold condition",
         ("ice-cold","ice-cold"), atmosphere="argon", method=method)
    step(r, "rt_hold", 8, "reaction", ["p_S3","p_THF","p_KOtBu","p_6"],
         "stirred at room temperature until TLC completion",
         ("room temperature","room temperature"),
         ("5 min to 2 h",{"minimum":5,"maximum":120,"unit":"min","raw_text":"5 min to 2 h"}),
         "argon", method)
    workup_participants(r, 8, purification="15–20% ethyl acetate/petroleum ether (v/v)")


allgen = reaction("general_synthesis_allenamides_6", 8, "generic_scheme",
                  "General procedure for N-allenamides 6", extra_pages=(9,))
add_allene_participants(allgen, mid("S3_general"), mid("6_general"))
add_allene_steps(allgen)
set_procedure(allgen, 8, "S3 (0.25 mmol) in dry THF (3 mL) under argon; add KOtBu (42 mg, 0.38 mmol) ice-cold, then stir at room temperature for 5 min to 2 h depending on substrate; aqueous workup and chromatography.",
              "KOtBu added at ice-cold condition")
classify(allgen, 8, "DBCC")
allgen["review"]["unresolved"].append("The 5 min–2 h range is substrate-dependent; entry-specific reaction times are not reported.")
bundle["reactions"].append(allgen)

allene_data = {
    "a":("72.3 mg",79,"achiral"), "b":("51.4 mg",52,"achiral"),
    "c":("88.4 mg",80,"achiral"), "d":("88.2 mg",90,"achiral"),
    "e":("72.9 mg",66,"unspecified"), "i":("76.5 mg",75,"achiral"),
}
for suffix, (mass, yld, _) in allene_data.items():
    r = reaction(f"synthesis_6{suffix}", 9, "specific_experiment", f"Synthesis of N-allenamide 6{suffix}", ["other"], extra_pages=(8,))
    r["general_scheme_ids"] = [allgen["reaction_id"]]
    add_allene_participants(r, mid(f"S3{suffix}"), mid(f"6{suffix}"), "inherited_general_condition")
    add_allene_steps(r, "inherited_general_condition")
    r["reaction_scale"] = f("S3 (0.25 mmol)", {"reference_molecule_id":mid(f"S3{suffix}"),"amount":{"number":0.25,"unit":"mmol"}}, ev(8), "inherited_general_condition")
    r["procedure"]["general_procedure_reference"] = f("General Procedure for the Synthesis of N-Allenamides 6", allgen["reaction_id"], ev(8), "inherited_general_condition")
    outcome(r, f"amount_6{suffix}", "product_amount", mid(f"6{suffix}"), mass,
            {"number":float(mass.split()[0]),"unit":"mg"}, 9, "isolated product")
    outcome(r, f"yield_6{suffix}", "yield", mid(f"6{suffix}"), f"{yld}%", {"number":yld,"unit":"%"}, 9,
            "isolated product", "isolated")
    classify(r, 9, "DBCC")
    r["review"]["unresolved"].append("Entry-specific time within the reported 5 min–2 h general range is unreported.")
    bundle["reactions"].append(r)


# Present-work THBC general procedure and the distinct gram-scale run.
thbcgen = reaction("general_synthesis_THBC_1", 9, "generic_scheme", "General procedure for 1-vinyltetrahydro-β-carbolines 1")
participant(thbcgen, mid("6_general"), "reactant", "p_6", ev(9), {"amount":("0.1 mmol",{"number":0.1,"unit":"mmol"})})
participant(thbcgen, mid("7_general"), "reactant", "p_7", ev(9), {"amount":("0.13 mmol",{"number":0.13,"unit":"mmol"})})
participant(thbcgen, mid("Pd2(dba)3"), "catalyst", "p_Pd2dba3", ev(9),
            {"mass":("4.6 mg",{"number":4.6,"unit":"mg"}),"mol_percent":("5 mol%",{"number":5,"unit":"mol%"})}, mid("6_general"))
participant(thbcgen, mid("t-Bu-XantPhos"), "ligand", "p_ligand", ev(9),
            {"mass":("6.0 mg",{"number":6.0,"unit":"mg"}),"mol_percent":("12 mol%",{"number":12,"unit":"mol%"})}, mid("6_general"))
participant(thbcgen, mid("Ag2CO3"), "reagent", "p_Ag2CO3", ev(9),
            {"mass":("55.2 mg",{"number":55.2,"unit":"mg"}),"amount":("0.2 mmol",{"number":0.2,"unit":"mmol"})})
participant(thbcgen, mid("DMF"), "solvent", "p_DMF", ev(9), {"volume":("2 mL",{"number":2,"unit":"mL"})})
participant(thbcgen, mid("1_general"), "product", "p_1", ev(9))
step(thbcgen, "catalyst_premix", 9, "reaction", ["p_Pd2dba3","p_ligand","p_DMF"],
     "stirred catalyst and ligand in dry DMF", ("room temperature","room temperature"), ("30 min",{"number":30,"unit":"min"}), "argon")
step(thbcgen, "aryl_iodide_premix", 9, "reaction", ["p_Pd2dba3","p_ligand","p_DMF","p_7"],
     "added aryl iodide and continued stirring", ("room temperature","room temperature"), ("30 min",{"number":30,"unit":"min"}), "argon")
step(thbcgen, "reaction", 9, "reaction", ["p_6","p_7","p_Pd2dba3","p_ligand","p_Ag2CO3","p_DMF","p_1"],
     "Ag2CO3 and allenamide were added successively, then the mixture was stirred in an oil bath until TLC completion",
     ("80 °C",{"number":80,"unit":"°C"}), ("1 h",{"number":1,"unit":"h"}), "argon")
participant(thbcgen, mid("DCM"), "workup_agent", "p_DCM", ev(9), {"volume":("3 × 20 mL",{"number":20,"unit":"mL","count":3})}, stage="workup")
participant(thbcgen, mid("water"), "workup_agent", "p_water", ev(9), {"volume":("10 mL",{"number":10,"unit":"mL"})}, stage="workup")
participant(thbcgen, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", ev(9), stage="workup")
participant(thbcgen, mid("silica_gel"), "purification_agent", "p_silica", ev(9), stage="purification")
participant(thbcgen, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(9), stage="purification")
step(thbcgen, "workup", 9, "workup", ["p_DCM","p_water","p_Na2SO4"], "extracted with DCM, washed with water, dried, and concentrated")
step(thbcgen, "purification", 9, "purification", ["p_silica","p_petroleum_ether"], "silica-gel chromatography using 5–15% ethyl acetate–petroleum ether (v/v)")
set_procedure(thbcgen, 9, "Pd2(dba)3 (4.6 mg, 5 mol%) and t-Bu-XantPhos (6.0 mg, 12 mol%) in dry DMF (2 mL), 30 min at rt under argon; add aryl iodide (0.13 mmol), 30 min; add Ag2CO3 (55.2 mg, 0.2 mmol) and allenamide 6 (0.1 mmol) successively; 80 °C, 1 h.",
              "Catalyst/ligand premix; aryl iodide; Ag2CO3; allenamide")
classify(thbcgen, 9, "CC")
bundle["reactions"].append(thbcgen)

scale = reaction("gram_scale_1aa", 9, "specific_experiment", "Gram-scale synthesis of THBC 1aa", ["scale_up"])
for mol, role, pid, charges, ref in [
    ("6a","reactant","p_6a",{"mass":("1.0 g",{"number":1.0,"unit":"g"}),"amount":("2.73 mmol",{"number":2.73,"unit":"mmol"})},None),
    ("7a","reactant","p_7a",{"mass":("0.724 g",{"number":0.724,"unit":"g"}),"volume":("0.40 mL",{"number":0.40,"unit":"mL"}),"amount":("3.55 mmol",{"number":3.55,"unit":"mmol"})},None),
    ("Pd2(dba)3","catalyst","p_Pd2dba3",{"mass":("125.14 mg",{"number":125.14,"unit":"mg"}),"mol_percent":("5 mol%",{"number":5,"unit":"mol%"})},mid("6a")),
    ("t-Bu-XantPhos","ligand","p_ligand",{"mass":("163.36 mg",{"number":163.36,"unit":"mg"}),"mol_percent":("12 mol%",{"number":12,"unit":"mol%"})},mid("6a")),
    ("Ag2CO3","reagent","p_Ag2CO3",{"mass":("1.5 g",{"number":1.5,"unit":"g"}),"amount":("5.46 mmol",{"number":5.46,"unit":"mmol"})},None),
    ("DMF","solvent","p_DMF",{"volume":("55 mL",{"number":55,"unit":"mL"})},None),
    ("1aa","product","p_1aa",{},None),
]: participant(scale, mid(mol), role, pid, ev(9), charges, ref)
scale["reaction_scale"] = f("allenamide 6a (1.0 g, 2.73 mmol)", {"reference_molecule_id":mid("6a"),"amount":{"number":2.73,"unit":"mmol"}}, ev(9))
step(scale, "catalyst_premix", 9, "reaction", ["p_Pd2dba3","p_ligand","p_DMF"], "stirred catalyst and ligand in dry DMF", ("room temperature","room temperature"), ("30 min",{"number":30,"unit":"min"}), "argon")
step(scale, "aryl_iodide_premix", 9, "reaction", ["p_Pd2dba3","p_ligand","p_DMF","p_7a"], "added 7a and continued stirring", ("room temperature","room temperature"), ("30 min",{"number":30,"unit":"min"}), "argon")
step(scale, "reaction", 9, "reaction", ["p_6a","p_7a","p_Pd2dba3","p_ligand","p_Ag2CO3","p_DMF","p_1aa"], "Ag2CO3 and 6a were added successively; stirred in an oil bath until TLC completion", ("80 °C",{"number":80,"unit":"°C"}), ("1 h",{"number":1,"unit":"h"}), "argon")
participant(scale, mid("DCM"), "workup_agent", "p_DCM", ev(9), {"volume":("3 × 30 mL",{"number":30,"unit":"mL","count":3})}, stage="workup")
participant(scale, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", ev(9), stage="workup")
participant(scale, mid("silica_gel"), "purification_agent", "p_silica", ev(9), stage="purification")
participant(scale, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(9), stage="purification")
step(scale, "workup", 9, "workup", ["p_DCM","p_Na2SO4"], "extracted with dichloromethane, dried, filtered, and concentrated")
step(scale, "purification", 9, "purification", ["p_silica","p_petroleum_ether"], "silica-gel chromatography using 5% ethyl acetate–petroleum ether (v/v)")
set_procedure(scale, 9, "Gram-scale run: catalyst/ligand premix in DMF, add 7a, then Ag2CO3 and 6a; 80 °C for 1 h; DCM extraction and chromatography.", "Catalyst/ligand premix; 7a; Ag2CO3; 6a")
outcome(scale, "amount_1aa", "product_amount", mid("1aa"), "0.787 g", {"number":0.787,"unit":"g"}, 9, "isolated product")
outcome(scale, "yield_1aa", "yield", mid("1aa"), "65%", {"number":65,"unit":"%"}, 9, "isolated product", "isolated")
classify(scale, 9, "CC")
bundle["reactions"].append(scale)


# Generic procedures on page 12.
betagen = reaction("general_synthesis_beta_carbolines_2", 12, "generic_scheme", "General procedure for β-carbolines 2")
participant(betagen, mid("1_general"), "reactant", "p_1", ev(12), {"amount":("0.05 mmol",{"number":0.05,"unit":"mmol"})})
participant(betagen, mid("DDQ"), "reagent", "p_DDQ", ev(12), {"mass":("34.1 mg",{"number":34.1,"unit":"mg"}),"amount":("0.15 mmol",{"number":0.15,"unit":"mmol"})})
participant(betagen, mid("toluene"), "solvent", "p_toluene", ev(12), {"volume":("1.5 mL",{"number":1.5,"unit":"mL"})})
participant(betagen, mid("2_general"), "product", "p_2", ev(12))
step(betagen, "reaction", 12, "reaction", ["p_1","p_DDQ","p_toluene","p_2"], "heated in an oil bath until TLC completion", ("100 °C",{"number":100,"unit":"°C"}), atmosphere="argon")
participant(betagen, mid("NaHCO3"), "workup_agent", "p_NaHCO3", ev(12), stage="workup")
participant(betagen, mid("ethyl_acetate"), "workup_agent", "p_EtOAc", ev(12), {"volume":("3 × 20 mL",{"number":20,"unit":"mL","count":3})}, stage="workup")
participant(betagen, mid("brine"), "workup_agent", "p_brine", ev(12), stage="workup")
participant(betagen, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", ev(12), stage="workup")
participant(betagen, mid("silica_gel"), "purification_agent", "p_silica", ev(12), stage="purification")
participant(betagen, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(12), stage="purification")
step(betagen, "workup", 12, "workup", ["p_NaHCO3","p_EtOAc","p_brine","p_Na2SO4"], "evaporated solvent, quenched with saturated sodium bicarbonate, extracted, washed with brine, dried, and concentrated")
step(betagen, "purification", 12, "purification", ["p_silica","p_petroleum_ether"], "silica-gel chromatography using 15–25% ethyl acetate–petroleum ether (v/v)")
set_procedure(betagen, 12, "1 (0.05 mmol) in dry toluene (1.5 mL), DDQ (34.1 mg, 0.15 mmol), argon, 100 °C until TLC completion; bicarbonate quench, extraction, and chromatography.", "DDQ added to 1 in dry toluene")
# Aromatization has at least C-C and C-N bond-order promotions; preserve both labels.
classify(betagen, 12, "DBCC")
change2 = {"change_id":"b2","atoms":[{"atom_map":None,"element":"C"},{"atom_map":None,"element":"N"}],"before_order":1,"after_order":2,"change_type":"order_increased"}
betagen["classification"]["bond_changes"].append(f(change2, change2, ev(12), "structure_derivation", "Second promoted bond inferred from aromatization without atom mapping."))
labs = [{"category_number":13,"code":"DB-C-C","supporting_change_ids":["b1"]},{"category_number":16,"code":"DB-C-N","supporting_change_ids":["b2"]}]
betagen["classification"]["labels"] = f(labs, labs, ev(12), "structure_derivation", "Partial aromatization labels inferred from the source transformation.")
bundle["reactions"].append(betagen)

ninegen = reaction("general_synthesis_spiroindolenines_9", 12, "generic_scheme", "General procedure for aza-spiroindolenines 9")
for mol, role, pid, charges, ref in [
    ("6h","reactant","p_6h",{"amount":("0.1 mmol",{"number":0.1,"unit":"mmol"})},None),
    ("7_general","reactant","p_7",{"amount":("0.11 mmol",{"number":0.11,"unit":"mmol"})},None),
    ("Pd(PPh3)4","catalyst","p_Pd",{"mass":("5.8 mg",{"number":5.8,"unit":"mg"}),"mol_percent":("5 mol%",{"number":5,"unit":"mol%"})},mid("6h")),
    ("PCy3","ligand","p_PCy3",{"mass":("5.9 mg",{"number":5.9,"unit":"mg"}),"mol_percent":("12 mol%",{"number":12,"unit":"mol%"})},mid("6h")),
    ("Cs2CO3","reagent","p_Cs2CO3",{"mass":("130.4 mg",{"number":130.4,"unit":"mg"}),"amount":("0.4 mmol",{"number":0.4,"unit":"mmol"})},None),
    ("DMF","solvent","p_DMF",{"volume":("2 mL",{"number":2,"unit":"mL"})},None),
    ("9_general","product","p_9",{},None),
]: participant(ninegen, mid(mol), role, pid, ev(12), charges, ref)
step(ninegen, "premix", 12, "reaction", ["p_Pd","p_PCy3","p_DMF","p_7"], "Pd(PPh3)4/PCy3 solution was treated with aryl iodide and stirred", ("room temperature","room temperature"), ("30 min",{"number":30,"unit":"min"}))
step(ninegen, "reaction", 12, "reaction", ["p_6h","p_7","p_Pd","p_PCy3","p_Cs2CO3","p_DMF","p_9"], "Cs2CO3 and 6h were added successively, then heated for the requisite time", ("100 °C",{"number":100,"unit":"°C"}))
participant(ninegen, mid("water"), "workup_agent", "p_water", ev(12), {"volume":("10 mL",{"number":10,"unit":"mL"})}, stage="workup")
participant(ninegen, mid("ethyl_acetate"), "workup_agent", "p_EtOAc", ev(12), {"volume":("3 × 10 mL",{"number":10,"unit":"mL","count":3})}, stage="workup")
participant(ninegen, mid("brine"), "workup_agent", "p_brine", ev(12), {"volume":("10 mL",{"number":10,"unit":"mL"})}, stage="workup")
participant(ninegen, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", ev(12), stage="workup")
participant(ninegen, mid("silica_gel"), "purification_agent", "p_silica", ev(12), stage="purification")
participant(ninegen, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(12), stage="purification")
step(ninegen, "workup", 12, "workup", ["p_water","p_EtOAc","p_brine","p_Na2SO4"], "evaporated DMF, added water, extracted with ethyl acetate, washed with brine, dried, and concentrated")
step(ninegen, "purification", 12, "purification", ["p_silica","p_petroleum_ether"], "silica-gel chromatography using 30–50% ethyl acetate–petroleum ether (v/v)")
set_procedure(ninegen, 12, "Pd(PPh3)4 (5.8 mg, 5 mol%) and PCy3 (5.9 mg, 12 mol%) in DMF (2 mL), add aryl iodide (0.11 mmol), stir 30 min at rt, then Cs2CO3 (130.4 mg, 0.4 mmol) and 6h (0.1 mmol), 100 °C for the requisite time.", "Aryl iodide; after 30 min Cs2CO3, then 6h")
classify(ninegen, 12, "CC")
ninegen["review"]["unresolved"].append("Entry-specific reaction times and the continuation of the diastereomer discussion occur outside page 12; existing page-6 scope records retain their reported times/dr.")
bundle["reactions"].append(ninegen)


# Update existing records with same-experiment experimental/characterization evidence.
def append_evidence(r, page):
    for eid in ev(page):
        if eid not in r["evidence_ids"]:
            r["evidence_ids"].append(eid)


table2_mass = {
    "1aa":(9,"33.2 mg"), "1ab":(9,"47.5 mg"), "1ac":(10,"29.5 mg"),
    "1ad":(10,"24.6 mg"), "1ae":(10,"39.8 mg"), "1af":(10,"29.9 mg"),
    "1ag":(10,"22.8 mg"), "1ah":(10,"22.7 mg"), "1ba":(10,"27.0 mg"),
    "1bb":(10,"32.9 mg"), "1ca":(10,"31.1 mg"), "1ch":(10,"25.8 mg"),
    "1db":(11,"23.7 mg"), "1de":(11,"19.8 mg"), "1eb":(11,"36.9 mg"),
    "1fb":(11,"43.5 mg"), "1ga":(11,"41.1 mg"), "1gb":(11,"37.0 mg"),
    "1ha":(11,"17.1 mg"), "1hb":(11,"21.9 mg"), "1hh":(11,"16.0 mg"),
    "1ia":(11,"36.8 mg"),
}
for label, (page, mass) in table2_mass.items():
    rid = f"{DOC}:table2_{label}"
    r = copy.deepcopy(old_reactions[rid])
    append_evidence(r, page)
    r["general_scheme_ids"] = list(dict.fromkeys(r["general_scheme_ids"] + [thbcgen["reaction_id"]]))
    r["procedure"]["general_procedure_reference"] = f("General Experimental Procedure for the Synthesis of 1-Vinyltetrahydro-β-carbolines 1", thbcgen["reaction_id"], ev(9), "inherited_general_condition")
    outcome(r, f"amount_{label}", "product_amount", mid(label), mass, {"number":float(mass.split()[0]),"unit":"mg"}, page, "isolated product")
    bundle["reactions"].append(r)

# Scheme 2: retain the page-4 identity and add page-11 quantities/workup.
r = copy.deepcopy(old_reactions[f"{DOC}:scheme2_6a_to_1a"])
append_evidence(r, 11)
by_pid = {p["participant_id"]:p for p in r["participants"]}
by_pid["p_6a"]["charges"]["amount"] = f("0.1 mmol",{"number":0.1,"unit":"mmol"},ev(11))
by_pid["p_pTsOH"]["charges"]["mass"] = f("1.9 mg",{"number":1.9,"unit":"mg"},ev(11))
by_pid["p_pTsOH"]["charges"]["amount"] = f("0.01 mmol",{"number":0.01,"unit":"mmol"},ev(11))
by_pid["p_THF"]["charges"]["volume"] = f("2 mL",{"number":2,"unit":"mL"},ev(11))
r["reaction_scale"] = f("6a (0.1 mmol)",{"reference_molecule_id":mid("6a"),"amount":{"number":0.1,"unit":"mmol"}},ev(11))
r["steps"] = []
step(r, "reaction", 11, "reaction", ["p_6a","p_pTsOH","p_THF","p_1a"], "stirred at room temperature until TLC completion", ("room temperature","room temperature"), atmosphere="argon")
r["steps"][-1]["time"] = f("overnight; detailed procedure: until completion (TLC)", "overnight", ev(4) + ev(11))
r["notes"].append(f("The discussion/Scheme 2 states overnight; the detailed page-11 procedure states until TLC completion. These complementary source expressions are both retained.",
                    "overnight / until TLC completion", ev(4) + ev(11)))
participant(r, mid("silica_gel"), "purification_agent", "p_silica", ev(11), stage="purification")
participant(r, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(11), stage="purification")
step(r, "solvent_removal", 11, "workup", [], "solvents removed under reduced pressure")
step(r, "purification", 11, "purification", ["p_silica","p_petroleum_ether"], "silica-gel chromatography using 10% ethyl acetate–petroleum ether (v/v)")
set_procedure(r, 11, "6a (0.1 mmol) in anhydrous THF (2 mL), p-TsOH (1.9 mg, 0.01 mmol), room temperature until TLC completion; evaporate and chromatograph.", "p-TsOH added to 6a in THF")
outcome(r, "amount_1a", "product_amount", mid("1a"), "15.7 mg", {"number":15.7,"unit":"mg"}, 11, "isolated product")
bundle["reactions"].append(r)

# Table 4: add the page-12 general procedure and isolated product masses to the same scope records.
table4_mass = {"2a":"11.4 mg","2b":"11.8 mg","2c":"12.1 mg","2d":"14.5 mg","2e":"13.3 mg","2f":"12.3 mg","2g":"16.0 mg","2h":"18.1 mg"}
for label, mass in table4_mass.items():
    r = copy.deepcopy(old_reactions[f"{DOC}:table4_{label}"])
    append_evidence(r, 12)
    r["general_scheme_ids"] = list(dict.fromkeys(r["general_scheme_ids"] + [betagen["reaction_id"]]))
    r["procedure"]["general_procedure_reference"] = f("General Experimental Procedure for the Synthesis of β-Carbolines 2", betagen["reaction_id"], ev(12), "inherited_general_condition")
    outcome(r, f"amount_{label}", "product_amount", mid(label), mass, {"number":float(mass.split()[0]),"unit":"mg"}, 12, "isolated product")
    bundle["reactions"].append(r)

# Table 3 entry 3 is the same product-8 experiment described in detail on page 12.
r = copy.deepcopy(old_reactions[f"{DOC}:table3_entry_03"])
append_evidence(r, 12)
by_pid = {p["participant_id"]:p for p in r["participants"]}
by_pid["p_1aa"]["charges"]["mass"] = f("44 mg",{"number":44,"unit":"mg"},ev(12))
by_pid["p_1aa"]["charges"]["amount"] = f("0.1 mmol",{"number":0.1,"unit":"mmol"},ev(12))
by_pid["p_KOtBu"]["charges"]["mass"] = f("22.4 mg",{"number":22.4,"unit":"mg"},ev(12))
by_pid["p_KOtBu"]["charges"]["amount"] = f("0.2 mmol",{"number":0.2,"unit":"mmol"},ev(12))
by_pid["p_dioxane"]["charges"]["volume"] = f("2 mL",{"number":2,"unit":"mL"},ev(12))
r["reaction_scale"] = f("1aa (44 mg, 0.1 mmol)",{"reference_molecule_id":mid("1aa"),"amount":{"number":0.1,"unit":"mmol"}},ev(12))
participant(r, mid("ethyl_acetate"), "workup_agent", "p_EtOAc", ev(12), {"volume":("3 × 20 mL",{"number":20,"unit":"mL","count":3})}, stage="workup")
participant(r, mid("brine"), "workup_agent", "p_brine", ev(12), stage="workup")
participant(r, mid("sodium_sulfate"), "workup_agent", "p_Na2SO4", ev(12), stage="workup")
participant(r, mid("silica_gel"), "purification_agent", "p_silica", ev(12), stage="purification")
participant(r, mid("petroleum_ether"), "purification_agent", "p_petroleum_ether", ev(12), stage="purification")
r["steps"] = []
step(r, "reaction", 12, "reaction", ["p_1aa","p_KOtBu","p_dioxane","p_2a","p_8"], "heated at reflux in an oil bath until TLC completion", ("110 °C",{"number":110,"unit":"°C"}), ("5 h",{"number":5,"unit":"h"}), "argon")
r["steps"][-1]["time"] = f("5 h (Table 3); detailed procedure: until completion (TLC)", {"number":5,"unit":"h"}, ev(5) + ev(12))
step(r, "workup", 12, "workup", ["p_EtOAc","p_brine","p_Na2SO4"], "evaporated solvent, extracted with ethyl acetate, washed with brine, dried, and concentrated")
step(r, "purification", 12, "purification", ["p_silica","p_petroleum_ether"], "column chromatography using 15% ethyl acetate–petroleum ether (v/v)")
set_procedure(r, 12, "1aa (44 mg, 0.1 mmol) in dry 1,4-dioxane (2 mL), KOtBu (22.4 mg, 0.2 mmol), argon, reflux at 110 °C until TLC completion; extraction and chromatography.", "KOtBu added to 1aa in dioxane")
outcome(r, "amount_8", "product_amount", mid("8"), "38.0 mg", {"number":38.0,"unit":"mg"}, 12, "isolated product")
bundle["reactions"].append(r)

# Scheme 3 scope entries receive the detailed page-12 general-procedure link.
for label in ("9a","9b","9i","9g"):
    r = copy.deepcopy(old_reactions[f"{DOC}:scheme3_{label}"])
    append_evidence(r, 12)
    r["general_scheme_ids"] = list(dict.fromkeys(r["general_scheme_ids"] + [ninegen["reaction_id"]]))
    r["procedure"]["general_procedure_reference"] = f("Experimental Procedure for the Synthesis of Products 9", ninegen["reaction_id"], ev(12), "inherited_general_condition")
    bundle["reactions"].append(r)


bundle["molecules"] = list(molecules.values())
bundle["source_context"] = [
    {"context_id":"page7_general_information_and_spectra","document_id":DOC,"evidence_ids":ev(7),"type":"characterization_context","text":"General analytical methods and NMR/HRMS characterization blocks for S2b, S2c, and the start of S2d are characterization, not separate synthesis experiments."},
    {"context_id":"page8_spectral_blocks","document_id":DOC,"evidence_ids":ev(8),"type":"characterization_context","text":"Continuation of S2d plus S2e and S3b–S3i characterization (NMR/HRMS, melting points) is retained as source context; isolated masses/yields are attached to the corresponding synthesis records."},
    {"context_id":"page9_spectral_blocks","document_id":DOC,"evidence_ids":ev(9),"type":"characterization_context","text":"NMR/HRMS and melting-point data for 6a–6i and THBC 1aa/1ab/1ac are characterization; isolated masses/yields support the same synthesis identities and are not new reactions."},
    {"context_id":"page10_thbc_spectra","document_id":DOC,"evidence_ids":ev(10),"type":"characterization_context","text":"Page 10 consists of THBC product characterization continuations; these repeat Table 2 experiment identities rather than defining new experiments."},
    {"context_id":"page11_thbc_spectra","document_id":DOC,"evidence_ids":ev(11),"type":"characterization_context","text":"Page 11 continues THBC characterization. The 6a→1a experimental paragraph is merged into the existing Scheme 2 record."},
    {"context_id":"page12_product_spectra","document_id":DOC,"evidence_ids":ev(12),"type":"characterization_context","text":"NMR/HRMS blocks for 1a and β-carbolines 2a–2h are characterization; product masses/yields are linked to the existing Scheme 2/Table 4 reactions. Product 8 and product-9 general procedures are merged with existing experiment identities."},
]

new_batch_ids = {
    f"{DOC}:scheme6_proposed_mechanism", f"{DOC}:control_9a_ag2co3_to_1ha",
    f"{DOC}:general_synthesis_S2", f"{DOC}:general_synthesis_S3",
    f"{DOC}:general_synthesis_allenamides_6", f"{DOC}:general_synthesis_THBC_1",
    f"{DOC}:gram_scale_1aa", f"{DOC}:general_synthesis_beta_carbolines_2",
    f"{DOC}:general_synthesis_spiroindolenines_9",
    *{f"{DOC}:synthesis_S2{x}" for x in "bcde"},
    *{f"{DOC}:synthesis_S3{x}" for x in "bcdefgi"},
    *{f"{DOC}:synthesis_6{x}" for x in "abcdei"},
}
raw_inventory = {
    "document_id": DOC,
    "page_range": [7, 12],
    "semantic_accuracy_verified": False,
    "new_reaction_ids": [r["reaction_id"] for r in bundle["reactions"] if r["reaction_id"] in new_batch_ids],
    "updated_reaction_ids": [r["reaction_id"] for r in bundle["reactions"] if r["reaction_id"] not in new_batch_ids],
    "notes": [
        "All six page images were visually inspected; extracted text was used only alongside image evidence.",
        "Spectral blocks are source context, not synthetic experiments.",
        "Complex structures without a sufficiently reliable encoding retain null SMILES and explicit stereochemistry/representation states.",
        "No SI association was inferred."
    ]
}
dump(RAW, raw_inventory)
dump(RUN / "raw_model/stereo_diagnostics_pages_0007_0012.json", diagnostics)
dump(SUBMITTED, bundle)
print(json.dumps({"submitted_reactions":len(bundle["reactions"]),
                  "new_reactions":len(raw_inventory["new_reaction_ids"]),
                  "updated_reactions":len(raw_inventory["updated_reaction_ids"]),
                  "bundle_molecules":len(bundle["molecules"])}, ensure_ascii=False))
