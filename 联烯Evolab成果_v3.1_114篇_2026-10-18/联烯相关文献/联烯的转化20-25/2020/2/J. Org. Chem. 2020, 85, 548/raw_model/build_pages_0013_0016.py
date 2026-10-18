#!/usr/bin/env python3
"""Merge visually reviewed PDF pages 13-16 into the v3.1 records.

Pages 13-14 provide detailed procedures for existing Table 3 experiments.
Pages 15-16 contain references only.  No new reaction is created here.
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from standard_records import fact, new_molecule, new_participant, new_step, read, stereo_inspect, write  # noqa: E402

DOC = "doc_66a75af3fbad9e13"
OUT = ROOT / "runs" / "current"
BASE = read(OUT / "standard_records.json")
SUBMITTED = read(OUT / "submitted_pages_0013_0016.json")
P = {i: f"{DOC}:page{i}:image" for i in range(13, 17)}
P7 = f"{DOC}:page7:image"


def F(raw=None, norm=None, ev=(), method="original_extraction", status=None,
      verification=None, derivation=None, alternatives=()):
    return fact(raw, norm, ev, method, status, verification, derivation, alternatives)


def NA(ev=()):
    return fact(evidence=ev, method="not_applicable", status="not_applicable", verification="missing")


def q(number, unit, raw=None):
    return {"number": number, "unit": unit, "raw_text": raw or f"{number} {unit}"}


def qr(minimum, maximum, unit, raw):
    return {"minimum": minimum, "maximum": maximum, "unit": unit, "raw_text": raw}


old_molecules = {m["molecule_id"]: m for m in BASE["molecules"]}
old_reactions = {r["reaction_id"]: r for r in BASE["reactions"]}
new_molecules = {}
updated_reactions = []
stereo_checks = {}


def mid(label):
    return f"{DOC}:mol_{label}"


def add_simple_molecule(label, name, smiles, ev):
    molecule_id = mid(label)
    diag = stereo_inspect(smiles)
    stereo_checks[molecule_id] = diag
    m = new_molecule(molecule_id, DOC, ev)
    m["identity_status"] = "source_identified"
    m["occurrences"][0]["compound_label"] = F(label, label, ev)
    m["names"]["original"] = F(name, name, ev)
    m["names"]["english"] = F(name, name, ev)
    m["names"]["abbreviations"] = F(label, [label], ev)
    m["names"]["aliases"] = fact(evidence=ev)
    m["cas_numbers"] = fact(evidence=ev)
    m["external_properties"] = []
    m["generic_template"] = None
    m["representation_type"] = "smiles"
    derivation = (
        "Candidate connectivity was assigned from the explicit reagent name and checked "
        "with the bundled stereo inspector; parsing is not source verification."
    )
    m["smiles"] = F(name, diag["canonical_isomeric_smiles"], ev, "structure_derivation", derivation=derivation)
    st = m["stereochemistry"]
    st["isomeric_smiles"] = F(name, diag["canonical_isomeric_smiles"], ev, "structure_derivation", derivation=derivation)
    st["identity_key"] = F(name, diag["stereo_graph_key"], ev, "structure_derivation", derivation=derivation)
    st["atom_centers"] = F(name, diag["atom_centers"], ev, "structure_derivation", derivation=derivation)
    st["double_bonds"] = F(name, diag["double_bonds"], ev, "structure_derivation", derivation=derivation)
    st["axial"] = NA(ev)
    st["mixture"] = fact(evidence=ev)
    st["assignment"] = F(name, {"tool_diagnostic": diag["potential_stereo"]}, ev, "structure_derivation", derivation=derivation)
    st["chirality_status"] = F(name, "achiral", ev, "structure_derivation", derivation="No stereogenic center or axis is present in the candidate inorganic reagent graph.")
    st["geometry_status"] = F(name, "not_applicable", ev, "structure_derivation", derivation="No E/Z or ring cis/trans element is applicable.")
    st["representation_status"] = F(name, "encoded", ev, "structure_derivation", derivation="The named reagent connectivity is encoded; chemical correctness remains pending review.")
    m["identity_note"] = F(
        "Named workup reagent from the detailed procedure; no CAS number was inferred.",
        "Named workup reagent from the detailed procedure; no CAS number was inferred.", ev, "inference",
        derivation="Machine extraction decision.")
    new_molecules[molecule_id] = m
    return molecule_id


add_simple_molecule("nah2po4", "sodium dihydrogen phosphate", "[Na+].OP(=O)(O)[O-]", [P[13], P[14]])


def append_unique(seq, value):
    if value not in seq:
        seq.append(value)


def add_evidence(f, ev):
    for eid in ev:
        append_unique(f["evidence_ids"], eid)


def get_reaction(entry, pages):
    rid = f"{DOC}:table3_entry_{entry}"
    r = copy.deepcopy(old_reactions[rid])
    ev = [P[p] for p in pages]
    for eid in ev:
        append_unique(r["evidence_ids"], eid)
    add_evidence(r["experiment_categories"], ev)
    add_evidence(r["reaction_name"], ev)
    add_evidence(r["transformation_type"], ev)
    r["location"] = F(
        f"PDF page 7, Table 3, entry {entry}; detailed procedure on PDF page(s) {', '.join(map(str, pages))}",
        f"PDF page 7, Table 3, entry {entry}; detailed procedure on PDF page(s) {', '.join(map(str, pages))}",
        [P7, *ev])
    r["review"]["status"] = "needs_review"
    r["review"]["method"] = "model"
    r["completeness"] = "partial"
    r["notes"].append(F(
        "The detailed experimental paragraph was merged into the existing Table 3 row; it is not a second experiment.",
        "Detailed-procedure duplicate merged", ev, "inference",
        derivation="Identity matched by substrate, product, nucleophile, conditions, and source heading."))
    updated_reactions.append(r)
    return r


def charge(r, participant_id, ev, *, mass=None, amount=None, volume=None, equivalents=None):
    p = next(p for p in r["participants"] if p["participant_id"] == participant_id)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume), ("equivalents", equivalents)):
        if value is not None:
            raw, norm = value
            p["charges"][key] = F(raw, norm, ev)
    return p


def participant(r, pid, molecule, role, ev, stage, *, mass=None, amount=None, volume=None):
    molecule_id = mid(molecule)
    if molecule_id not in old_molecules and molecule_id not in new_molecules:
        raise KeyError(molecule_id)
    p = new_participant(pid, molecule_id, role, ev, stage)
    for key, value in (("mass", mass), ("amount", amount), ("volume", volume)):
        if value is not None:
            raw, norm = value
            p["charges"][key] = F(raw, norm, ev)
    r["participants"].append(p)
    return p


def step(r, sid, ev, pids, operation, *, temp=None, time=None, stage="reaction"):
    s = new_step(sid, ev, stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, ev)
    if temp is not None:
        s["temperature"] = F(temp[0], temp[1], ev)
    if time is not None:
        s["time"] = F(time[0], time[1], ev)
    r["steps"].append(s)
    return s


def exact_yield(r, outcome_id, product_mass, yield_value, ev, basis):
    o = next(o for o in r["outcomes"] if o["outcome_id"] == outcome_id)
    o["value"] = F(f"{yield_value}%", q(yield_value, "%"), ev)
    o["yield_type"] = F("isolated yield", "isolated", ev)
    o["measurement_basis"] = F(basis, "isolated", ev)
    product = next(p for p in r["participants"] if p["role"] == "product")
    product["charges"]["mass"] = F(product_mass[0], product_mass[1], ev)


def conflict_fact(raw_a, norm_a, ev_a, raw_b, norm_b, ev_b):
    return F(None, None, [*ev_a, *ev_b], status="conflict", verification="missing", alternatives=[
        F(raw_a, norm_a, ev_a),
        F(raw_b, norm_b, ev_b),
    ])


def common_workup(r, prefix, ev, specs, operation):
    ids = []
    for label, role, volume in specs:
        pid = f"{prefix}_{label}"
        participant(r, pid, label, role, ev, "workup" if role == "workup_agent" else "purification",
                    volume=volume)
        ids.append(pid)
    step(r, f"{prefix}_workup", ev, [pid for pid in ids if "purification" not in next(p["stage"] for p in r["participants"] if p["participant_id"] == pid)], operation, stage="workup")
    return ids


def reset_detailed_steps(r):
    # Replace the earlier table-summary steps with the ordered detailed procedure.
    r["steps"] = []


# Entry 2: product 46.  Table and paragraph disagree on time and yield.
r = get_reaction(2, [13])
ev = [P[13]]
charge(r, "t3e2_substrate", ev, mass=("0.150 g", q(0.150, "g")), amount=("0.43 mmol", q(0.43, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e2_reagent", ev, volume=("2.04 mL of 0.62 M solution in THF", q(2.04, "mL")), amount=("1.28 mmol", q(1.28, "mmol")), equivalents=("3.0 equiv", q(3.0, "equiv")))
charge(r, "t3e2_solvent", ev, volume=("5.5 mL", q(5.5, "mL")))
reset_detailed_steps(r)
step(r, "add_propylmagnesium_bromide", ev, ["t3e2_substrate", "t3e2_reagent", "t3e2_solvent"], "add propylmagnesium bromide dropwise", time=("30 s", q(30, "s")))
s = step(r, "heat_at_reflux", [P7, P[13]], [], "heat at reflux")
s["temperature"] = F("reflux", "reflux", [P7, P[13]])
s["time"] = conflict_fact("13 h (Table 3)", q(13, "h"), [P7], "12 h (experimental paragraph)", q(12, "h"), ev)
participant(r, "t3e2_nh4cl", "nh4cl", "workup_agent", ev, "workup", volume=("2.5 mL saturated aqueous solution", q(2.5, "mL")))
participant(r, "t3e2_water", "water", "workup_agent", ev, "workup", volume=("3 mL", q(3, "mL")))
participant(r, "t3e2_etoac", "etoac", "workup_agent", ev, "workup", volume=("three 10 mL portions", q(30, "mL")))
participant(r, "t3e2_nacl", "nacl", "workup_agent", ev, "workup", volume=("30 mL saturated aqueous solution", q(30, "mL")))
participant(r, "t3e2_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e2_silica", "silica", "purification_agent", ev, "purification", mass=("15 g", q(15, "g")))
participant(r, "t3e2_ether", "ether", "purification_agent", ev, "purification")
participant(r, "t3e2_hexanes", "hexanes", "purification_agent", ev, "purification")
step(r, "workup", ev, ["t3e2_nh4cl", "t3e2_water", "t3e2_etoac", "t3e2_nacl", "t3e2_mgso4"], "cool, quench, extract, wash, dry, filter, and concentrate", stage="workup")
step(r, "purification", ev, ["t3e2_silica", "t3e2_ether", "t3e2_hexanes"], "silica-gel column chromatography with 2% ether/hexanes", stage="purification")
o = next(o for o in r["outcomes"] if o["outcome_id"] == "yield_46")
o["value"] = conflict_fact("77-82% (Table 3)", qr(77, 82, "%", "77-82%"), [P7], "80% (experimental paragraph)", q(80, "%"), ev)
o["measurement_basis"] = F("isolated after silica-gel column chromatography", "isolated", [P7, P[13]])
charge(r, "t3e2_product", ev, mass=("0.082 g", q(0.082, "g")))
r["review"]["unresolved"].extend([
    "Table 3 reports 13 h, while the detailed paragraph reports 12 h; both are retained without selection.",
    "Table 3 reports 77-82%, while the detailed paragraph reports 80%; both are retained without selection.",
])
r["notes"].append(F("Table 3: 13 h and 77-82%; experimental paragraph: 12 h and 80%.", "Source conflict retained", [P7, P[13]], "inference", derivation="The table and detailed paragraph were transcribed independently; the actual conflicting fields retain null selected values and both alternatives."))
r["procedure"]["original_text"] = F("Pyridine 16 and propylmagnesium bromide in THF; reflux; aqueous workup; silica-gel chromatography.", "detailed experimental procedure", ev)


# Entry 3: product 47.
r = get_reaction(3, [13])
ev = [P[13]]
charge(r, "t3e3_substrate", ev, mass=("0.070 g", q(0.070, "g")), amount=("0.199 mmol", q(0.199, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e3_reagent", ev, volume=("0.60 mL of 2.0 M solution in THF", q(0.60, "mL")), amount=("1.20 mmol", q(1.20, "mmol")), equivalents=("6.0 equiv", q(6.0, "equiv")))
charge(r, "t3e3_solvent", ev, volume=("1.4 mL", q(1.4, "mL")))
reset_detailed_steps(r)
step(r, "add_phenylmagnesium_chloride", ev, ["t3e3_substrate", "t3e3_reagent", "t3e3_solvent"], "add phenylmagnesium chloride")
step(r, "heat_at_reflux", ev, [], "heat at reflux", temp=("reflux", "reflux"), time=("16 h", q(16, "h")))
participant(r, "t3e3_nah2po4", "nah2po4", "workup_agent", ev, "workup", volume=("1 mL half-saturated aqueous solution", q(1, "mL")))
participant(r, "t3e3_etoac", "etoac", "workup_agent", ev, "workup", volume=("25 mL dilution plus three 20 mL portions", q(85, "mL")))
participant(r, "t3e3_water", "water", "workup_agent", ev, "workup", volume=("25 mL", q(25, "mL")))
participant(r, "t3e3_nacl", "nacl", "workup_agent", ev, "workup")
participant(r, "t3e3_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e3_silica", "silica", "purification_agent", ev, "purification", mass=("6 g", q(6, "g")))
participant(r, "t3e3_etoac_eluent", "etoac", "purification_agent", ev, "purification")
participant(r, "t3e3_hexanes", "hexanes", "purification_agent", ev, "purification")
step(r, "workup", ev, ["t3e3_nah2po4", "t3e3_etoac", "t3e3_water", "t3e3_nacl", "t3e3_mgso4"], "quench, dilute, extract, wash with brine, dry, filter, and concentrate", stage="workup")
step(r, "purification", ev, ["t3e3_silica", "t3e3_etoac_eluent", "t3e3_hexanes"], "silica-gel column chromatography, 0-10% EtOAc/hexanes", stage="purification")
exact_yield(r, "yield_47", ("0.047 g", q(0.047, "g")), 87, ev, "isolated after silica-gel column chromatography")
r["procedure"]["original_text"] = F("Pyridine 16 and PhMgCl in THF; reflux 16 h; aqueous workup; silica-gel chromatography.", "detailed experimental procedure", ev)


# Entry 4: product 48.
r = get_reaction(4, [13])
ev = [P[13]]
charge(r, "t3e4_substrate", ev, mass=("0.170 g", q(0.170, "g")), amount=("0.484 mmol", q(0.484, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e4_reagent", ev, mass=("0.118 g", q(0.118, "g")), amount=("2.418 mmol", q(2.418, "mmol")), equivalents=("5.0 equiv", q(5.0, "equiv")))
charge(r, "t3e4_extra", ev, mass=("0.640 g", q(0.640, "g")), amount=("2.418 mmol", q(2.418, "mmol")), equivalents=("5.0 equiv", q(5.0, "equiv")))
charge(r, "t3e4_solvent", ev, volume=("5 mL", q(5, "mL")))
reset_detailed_steps(r)
step(r, "charge_reagents", ev, ["t3e4_substrate", "t3e4_reagent", "t3e4_extra", "t3e4_solvent"], "charge pyridine, sodium cyanide, 18-crown-6, and DMF")
step(r, "heat_at_reflux", ev, [], "heat at reflux", temp=("reflux", "reflux"), time=("48 h", q(48, "h")))
participant(r, "t3e4_water", "water", "workup_agent", ev, "workup", volume=("100 mL dilution plus five 50 mL washes", q(350, "mL")))
participant(r, "t3e4_etoac", "etoac", "workup_agent", ev, "workup", volume=("100 mL", q(100, "mL")))
participant(r, "t3e4_nacl", "nacl", "workup_agent", ev, "workup", volume=("50 mL saturated aqueous solution", q(50, "mL")))
participant(r, "t3e4_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e4_ch2cl2", "ch2cl2", "purification_agent", ev, "purification", volume=("10 mL", q(10, "mL")))
participant(r, "t3e4_silica_load", "silica", "purification_agent", ev, "purification", mass=("0.60 g for dry loading", q(0.60, "g")))
participant(r, "t3e4_silica_column", "silica", "purification_agent", ev, "purification", mass=("4 g column", q(4, "g")))
participant(r, "t3e4_etoac_eluent", "etoac", "purification_agent", ev, "purification")
participant(r, "t3e4_hexanes", "hexanes", "purification_agent", ev, "purification", volume=("three 1 mL portions then three 0.25 mL portions for trituration", q(3.75, "mL")))
step(r, "workup", ev, ["t3e4_water", "t3e4_etoac", "t3e4_nacl", "t3e4_mgso4"], "cool, dilute, separate, wash, dry, filter, and concentrate", stage="workup")
step(r, "column", ev, ["t3e4_ch2cl2", "t3e4_silica_load", "t3e4_silica_column", "t3e4_etoac_eluent", "t3e4_hexanes"], "dry-load and chromatograph on silica with 0-10% EtOAc/hexanes", stage="purification")
step(r, "trituration", ev, ["t3e4_hexanes"], "triturate the first crop and concentrated mother liquor with hexanes", stage="purification")
exact_yield(r, "yield_48", ("0.093 g total", q(0.093, "g")), 86, ev, "isolated after silica-gel chromatography and hexanes trituration")
r["procedure"]["original_text"] = F("Pyridine 16, NaCN, and 18-crown-6 in DMF; reflux 48 h; aqueous workup; chromatography and trituration.", "detailed experimental procedure", ev)


# Entry 5: product 49.
r = get_reaction(5, [13])
ev = [P[13]]
charge(r, "t3e5_substrate", ev, mass=("0.189 g", q(0.189, "g")), amount=("0.54 mmol", q(0.54, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e5_reagent", ev, mass=("0.145 g", q(0.145, "g")), amount=("1.292 mmol", q(1.292, "mmol")), equivalents=("2.4 equiv", q(2.4, "equiv")))
charge(r, "t3e5_extra", ev, volume=("75 uL", q(75, "uL")), mass=("0.060 g", q(0.060, "g")), amount=("1.29 mmol", q(1.29, "mmol")), equivalents=("2.4 equiv", q(2.4, "equiv")))
charge(r, "t3e5_solvent", ev, volume=("3.8 mL", q(3.8, "mL")))
reset_detailed_steps(r)
step(r, "add_ethanol", ev, ["t3e5_reagent", "t3e5_extra", "t3e5_solvent"], "add ethanol dropwise to potassium tert-butoxide in THF", temp=("rt", "rt"))
step(r, "pre_stir", ev, [], "stir the ethoxide-forming mixture", temp=("rt", "rt"), time=("10 min", q(10, "min")))
step(r, "add_substrate", ev, ["t3e5_substrate"], "add sulfonylpyridine 16 via an open neck", temp=("rt", "rt"))
step(r, "reaction_hold", ev, [], "stir the resulting mixture", temp=("rt", "rt"), time=("40 min", q(40, "min")))
participant(r, "t3e5_nh4cl", "nh4cl", "workup_agent", ev, "workup", volume=("2 mL saturated aqueous solution", q(2, "mL")))
participant(r, "t3e5_water", "water", "workup_agent", ev, "workup", volume=("3 mL then 10 mL", q(13, "mL")))
participant(r, "t3e5_etoac", "etoac", "workup_agent", ev, "workup", volume=("10 mL dilution plus three 10 mL portions", q(40, "mL")))
participant(r, "t3e5_nacl", "nacl", "workup_agent", ev, "workup", volume=("30 mL saturated aqueous solution", q(30, "mL")))
participant(r, "t3e5_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e5_silica", "silica", "purification_agent", ev, "purification", mass=("2.5 g", q(2.5, "g")))
participant(r, "t3e5_etoac_eluent", "etoac", "purification_agent", ev, "purification")
participant(r, "t3e5_hexanes", "hexanes", "purification_agent", ev, "purification")
step(r, "workup", ev, ["t3e5_nh4cl", "t3e5_water", "t3e5_etoac", "t3e5_nacl", "t3e5_mgso4"], "quench, dilute, extract, wash, dry, filter, and concentrate", stage="workup")
step(r, "purification", ev, ["t3e5_silica", "t3e5_etoac_eluent", "t3e5_hexanes"], "silica-gel column chromatography, 0-1% EtOAc/hexanes", stage="purification")
exact_yield(r, "yield_49", ("0.118 g", q(0.118, "g")), 91, ev, "isolated after silica-gel column chromatography")
r["procedure"]["original_text"] = F("Generate potassium ethoxide from KOtBu and ethanol in THF; add 16; stir 40 min at rt; work up and chromatograph.", "detailed experimental procedure", ev)


# Entry 6: product 50.  Entry 7 remains a distinct table-only experiment.
r = get_reaction(6, [13])
ev = [P[13]]
charge(r, "t3e6_substrate", ev, mass=("0.065 g", q(0.065, "g")), amount=("0.20 mmol", q(0.20, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e6_reagent", ev, volume=("0.40 mL of 3.0 M solution in ether", q(0.40, "mL")), amount=("1.18 mmol", q(1.18, "mmol")), equivalents=("6.0 equiv", q(6.0, "equiv")))
charge(r, "t3e6_solvent", ev, volume=("2 mL", q(2, "mL")))
reset_detailed_steps(r)
step(r, "add_methylmagnesium_bromide", ev, ["t3e6_substrate", "t3e6_reagent", "t3e6_solvent"], "add methylmagnesium bromide dropwise", temp=("rt", "rt"), time=("1 min", q(1, "min")))
step(r, "heat", ev, [], "heat the sealed reaction solution", temp=("70 °C", q(70, "°C")), time=("1.5 h", q(1.5, "h")))
participant(r, "t3e6_nh4cl", "nh4cl", "workup_agent", ev, "workup", volume=("2 mL saturated aqueous solution", q(2, "mL")))
participant(r, "t3e6_water", "water", "workup_agent", ev, "workup", volume=("2 mL", q(2, "mL")))
participant(r, "t3e6_etoac", "etoac", "workup_agent", ev, "workup", volume=("three 10 mL portions", q(30, "mL")))
participant(r, "t3e6_nacl", "nacl", "workup_agent", ev, "workup", volume=("20 mL saturated aqueous solution", q(20, "mL")))
participant(r, "t3e6_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e6_silica", "silica", "purification_agent", ev, "purification")
participant(r, "t3e6_etoac_eluent", "etoac", "purification_agent", ev, "purification")
step(r, "workup", ev, ["t3e6_nh4cl", "t3e6_water", "t3e6_etoac", "t3e6_nacl", "t3e6_mgso4"], "cool, quench, extract, wash, dry, filter, and concentrate", stage="workup")
step(r, "purification", ev, ["t3e6_silica", "t3e6_etoac_eluent"], "preparative thin-layer chromatography with EtOAc", stage="purification")
exact_yield(r, "yield_50", ("0.020 g", q(0.020, "g")), 53, ev, "isolated after preparative thin-layer chromatography")
r["procedure"]["original_text"] = F("Sulfonylpyridine 7 and MeMgBr in THF; 70 °C for 1.5 h; workup and preparative TLC.", "detailed experimental procedure", ev)


# Entry 8: product 51.  Table and paragraph disagree on reaction time.
r = get_reaction(8, [13])
ev = [P[13]]
charge(r, "t3e8_substrate", ev, mass=("0.078 g", q(0.078, "g")), amount=("0.226 mmol", q(0.226, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e8_reagent", ev, volume=("0.9 mL of 3.0 M solution in THF", q(0.9, "mL")), amount=("2.7 mmol", q(2.7, "mmol")), equivalents=("12 equiv", q(12, "equiv")))
charge(r, "t3e8_solvent", ev, volume=("1 mL", q(1, "mL")))
reset_detailed_steps(r)
step(r, "add_methylmagnesium_chloride", ev, ["t3e8_substrate", "t3e8_reagent", "t3e8_solvent"], "add methylmagnesium chloride dropwise", temp=("rt", "rt"), time=("0.5 min", q(0.5, "min")))
s = step(r, "reaction_hold", [P7, P[13]], [], "stir the resulting solution", temp=("rt", "rt"))
s["time"] = conflict_fact("16 h (Table 3)", q(16, "h"), [P7], "18 h (experimental paragraph)", q(18, "h"), ev)
participant(r, "t3e8_nah2po4", "nah2po4", "workup_agent", ev, "workup", volume=("2 mL half-saturated aqueous solution", q(2, "mL")))
participant(r, "t3e8_water", "water", "workup_agent", ev, "workup", volume=("50 mL", q(50, "mL")))
participant(r, "t3e8_etoac", "etoac", "workup_agent", ev, "workup", volume=("50 mL dilution plus three 15 mL portions", q(95, "mL")))
participant(r, "t3e8_nacl", "nacl", "workup_agent", ev, "workup", volume=("20 mL saturated aqueous solution", q(20, "mL")))
participant(r, "t3e8_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e8_silica", "silica", "purification_agent", ev, "purification")
participant(r, "t3e8_etoac_eluent", "etoac", "purification_agent", ev, "purification")
participant(r, "t3e8_hexanes", "hexanes", "purification_agent", ev, "purification")
step(r, "workup", ev, ["t3e8_nah2po4", "t3e8_water", "t3e8_etoac", "t3e8_nacl", "t3e8_mgso4"], "quench, dilute, extract, wash, dry, filter, and concentrate", stage="workup")
step(r, "purification", ev, ["t3e8_silica", "t3e8_etoac_eluent", "t3e8_hexanes"], "preparative thin-layer chromatography with 10% EtOAc/hexanes", stage="purification")
exact_yield(r, "yield_51", ("0.036 g", q(0.036, "g")), 78, ev, "isolated after preparative thin-layer chromatography")
r["review"]["unresolved"].append("Table 3 reports 16 h, while the detailed paragraph reports 18 h; both are retained without selection.")
r["notes"].append(F("Table 3: 16 h; experimental paragraph: 18 h.", "Source conflict retained", [P7, P[13]], "inference", derivation="The table and detailed paragraph were transcribed independently; the actual conflicting field retains a null selected value and both alternatives."))
r["procedure"]["original_text"] = F("Sulfonylpyridine 26 and MeMgCl in THF at rt; workup and preparative TLC.", "detailed experimental procedure", ev)


# Entry 9: product 52, spanning pages 13-14.
r = get_reaction(9, [13, 14])
ev = [P[13], P[14]]
charge(r, "t3e9_substrate", ev, mass=("0.390 g", q(0.390, "g")), amount=("0.930 mmol", q(0.930, "mmol")), equivalents=("1.0 equiv", q(1.0, "equiv")))
charge(r, "t3e9_reagent", ev, volume=("5.6 mL of 1.0 M solution in THF", q(5.6, "mL")), amount=("5.6 mmol", q(5.6, "mmol")), equivalents=("6.0 equiv", q(6.0, "equiv")))
charge(r, "t3e9_extra", ev, volume=("33 uL", q(33, "uL")), mass=("0.033 g", q(0.033, "g")), amount=("1.860 mmol", q(1.860, "mmol")), equivalents=("2.0 equiv", q(2.0, "equiv")))
charge(r, "t3e9_solvent", ev, volume=("4 mL", q(4, "mL")))
reset_detailed_steps(r)
step(r, "add_potassium_tert_butoxide", [P[13]], ["t3e9_substrate", "t3e9_reagent", "t3e9_extra", "t3e9_solvent"], "add potassium tert-butoxide solution to pyridine 13 and water in THF")
step(r, "heat_at_reflux", [P[13]], [], "heat at reflux", temp=("reflux", "reflux"), time=("16 h", q(16, "h")))
participant(r, "t3e9_nah2po4", "nah2po4", "workup_agent", ev, "workup", volume=("3 mL half-saturated aqueous solution", q(3, "mL")))
participant(r, "t3e9_water", "water", "workup_agent", ev, "workup", volume=("50 mL", q(50, "mL")))
participant(r, "t3e9_etoac", "etoac", "workup_agent", ev, "workup", volume=("50 mL dilution plus three 25 mL portions", q(125, "mL")))
participant(r, "t3e9_nacl", "nacl", "workup_agent", ev, "workup")
participant(r, "t3e9_mgso4", "mgso4", "workup_agent", ev, "workup")
participant(r, "t3e9_hexanes", "hexanes", "purification_agent", [P[14]], "purification", volume=("three 1 mL portions", q(3, "mL")))
step(r, "workup", ev, ["t3e9_nah2po4", "t3e9_water", "t3e9_etoac", "t3e9_nacl", "t3e9_mgso4"], "cool, quench, dilute, extract, wash, dry, filter, and concentrate", stage="workup")
step(r, "purification", [P[14]], ["t3e9_hexanes"], "triturate with hexanes and concentrate", stage="purification")
exact_yield(r, "yield_52", ("0.196 g", q(0.196, "g")), 75, ev, "isolated after hexanes trituration and concentration")
r["procedure"]["original_text"] = F("Pyridine 13, water, and KOtBu in THF; reflux 16 h; aqueous workup; hexanes trituration.", "detailed experimental procedure", ev)


# Re-run the already-saved graph checks for products 46-52 to document that the
# current batch did not collapse stereoisomers or reinterpret unspecified stereo.
for label in ("46", "47", "48", "49", "50", "51", "52"):
    m = old_molecules[mid(label)]
    smiles = m["stereochemistry"]["isomeric_smiles"]["normalized_value"]
    if smiles:
        stereo_checks[mid(label)] = stereo_inspect(smiles)


source_context = [
    {"context_id": f"{DOC}:pages13_14_characterization", "document_id": DOC,
     "kind": "characterization_data", "evidence_ids": [P[13], P[14]],
     "note": "IR, NMR, HRMS, melting-point, crude-mass, and calculated/found mass data accompany compounds 46-52; they are characterization context, not independent synthesis experiments."},
    {"context_id": f"{DOC}:page14_si_notice", "document_id": DOC,
     "kind": "supporting_information_notice", "evidence_ids": [P[14]],
     "note": "The paper advertises Supporting Information, but source_mapping.json does not establish a supplied SI association; no SI was linked."},
    {"context_id": f"{DOC}:pages14_16_back_matter", "document_id": DOC,
     "kind": "bibliographic_context", "evidence_ids": [P[14], P[15], P[16]],
     "note": "Author information, acknowledgments, and references are back matter or cited prior work, not current-paper experiments."},
]


referenced = {p["molecule_id"] for r in updated_reactions for p in r["participants"]}
all_molecules = []
for molecule_id in sorted(referenced):
    all_molecules.append(copy.deepcopy(new_molecules.get(molecule_id, old_molecules.get(molecule_id))))

SUBMITTED["molecules"] = all_molecules
SUBMITTED["reactions"] = updated_reactions
SUBMITTED["source_context"] = source_context
write(OUT / "submitted_pages_0013_0016.json", SUBMITTED)
write(OUT / "submitted.json", SUBMITTED)
write(OUT / "raw_model" / "stereo_checks_pages_0013_0016.json", stereo_checks)

coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {"page": 13, "state": "attempted", "regions": [
            {"region": "experimental section, compounds 46-51", "finding": "Detailed procedures for Table 3 entries 2-6 and 8 were matched to and merged into the existing row-level experiments. Entry 7 remains a distinct table-only condition and was not conflated with the detailed compound-50 preparation.", "evidence_ids": [P[13]]},
            {"region": "lower-right start of compound 52 procedure", "finding": "The preparation corresponding to Table 3 entry 9 begins here and continues on page 14; both pages are linked to one experiment.", "evidence_ids": [P[13]]},
            {"region": "analytical paragraphs for compounds 46-51", "finding": "IR/NMR/HRMS peak lists and calculated/found masses were classified as characterization context, not synthesis experiments.", "evidence_ids": [P[13]]},
        ], "record_ids": [f"{DOC}:table3_entry_{i}" for i in (2, 3, 4, 5, 6, 8, 9)],
         "unresolved": ["Table 3 versus detailed-text conflicts remain unresolved for entry 2 time/yield (13 h and 77-82% versus 12 h and 80%) and entry 8 time (16 h versus 18 h)."]},
        {"page": 14, "state": "attempted", "regions": [
            {"region": "upper-left continuation of compound 52", "finding": "Workup, hexanes trituration, 0.196 g product, and 75% isolated yield complete the existing Table 3 entry 9 record.", "evidence_ids": [P[14]]},
            {"region": "compound 52 characterization", "finding": "Melting point, IR, NMR, and HRMS data were retained as characterization context rather than a reaction.", "evidence_ids": [P[14]]},
            {"region": "Associated Content, author information, acknowledgments, and references", "finding": "Back matter and cited prior work are not current-paper experiments. The SI notice was recorded, but no SI was associated because source_mapping.json leaves that relationship unestablished.", "evidence_ids": [P[14]]},
        ], "record_ids": [f"{DOC}:table3_entry_9"],
         "unresolved": ["The advertised Supporting Information was not supplied or explicitly associated, so it was not used."]},
        {"page": 15, "state": "attempted", "regions": [
            {"region": "full page, references 7-24", "finding": "Bibliographic citations and descriptions of prior work only; no current-paper experimental arrow, table row, or procedure was found.", "evidence_ids": [P[15]]},
        ], "record_ids": [], "unresolved": []},
        {"page": 16, "state": "attempted", "regions": [
            {"region": "full page, references 25-41", "finding": "Bibliographic citations and descriptions of prior work only; no current-paper experimental arrow, table row, or procedure was found.", "evidence_ids": [P[16]]},
        ], "record_ids": [], "unresolved": []},
    ],
}
write(OUT / "coverage" / "pages_0013_0016.json", coverage)

raw_summary = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "paper_title": "Synthesis of Highly Substituted Pyridines via [4 + 2] Cycloadditions of Vinylallenes and Sulfonyl Cyanides",
    "doi": "10.1021/acs.joc.9b02628",
    "page_range": [13, 16],
    "visual_pages_inspected": [13, 14, 15, 16],
    "semantic_accuracy_verified": False,
    "new_reaction_count": 0,
    "updated_existing_reaction_count": len(updated_reactions),
    "new_molecule_count": len(new_molecules),
    "decisions": {
        "deduplication": "Detailed procedures for compounds 46-52 were merged into existing Table 3 rows.",
        "table3_entry7": "Kept separate and unchanged; page 13 compound 50 is the entry-6 condition, not entry 7.",
        "characterization": "Independent spectral and mass data were retained only in source_context.",
        "pages15_16": "References only; no experiments invented.",
        "si": "Not associated; the source mapping is unestablished.",
    },
    "unresolved": [
        "Entry 2 table/text conflicts: 13 h and 77-82% versus 12 h and 80%.",
        "Entry 8 table/text conflict: 16 h versus 18 h.",
        "All machine records remain needs_review; validation is not chemical verification.",
    ],
}
write(OUT / "raw_model" / "raw_extraction_pages_0013_0016.json", raw_summary)

print(json.dumps({
    "new_reactions": 0,
    "updated_reactions": len(updated_reactions),
    "bundle_reactions": len(updated_reactions),
    "bundle_molecules": len(all_molecules),
    "new_molecules": len(new_molecules),
}, ensure_ascii=False))
