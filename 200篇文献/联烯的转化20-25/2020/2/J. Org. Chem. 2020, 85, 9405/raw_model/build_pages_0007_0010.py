#!/usr/bin/env python3
"""Enrich v3.1 records from inspected PDF pages 7--10."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

DOC = "doc_106b2f72bb5b02d2"
BASE = ROOT / "runs/current/standard_records.json"
DATASET = ROOT / "runs/current/dataset.json"
OUT_RAW = ROOT / "runs/current/raw_model/extracted_bundle_pages_0007_0010.json"
OUT_SUBMITTED = ROOT / "runs/current/submitted.json"
OUT_COVERAGE = ROOT / "runs/current/coverage/pages_0007_0010.json"


def ev(page):
    return [f"{DOC}:page{page}:text", f"{DOC}:page{page}:image"]


def evs(*pages):
    return [eid for page in pages for eid in ev(page)]


def uniq(items):
    seen = set()
    result = []
    for item in items:
        key = json.dumps(item, ensure_ascii=False, sort_keys=True) if isinstance(item, (dict, list)) else item
        if key not in seen:
            seen.add(key)
            result.append(item)
    return result


def f(raw=None, normalized=None, evidence=(), method="original_extraction", derivation=None,
      status=None, alternatives=()):
    return st.fact(raw, normalized, evidence, method=method, derivation=derivation,
                   status=status, alternatives=alternatives)


def q(number, unit, raw=None, qualifier=None):
    value = {"number": number, "unit": unit}
    if raw is not None:
        value["raw_text"] = raw
    if qualifier:
        value["qualifier"] = qualifier
    return value


dataset = json.loads(DATASET.read_text())
base = json.loads(BASE.read_text())
all_molecules = {m["molecule_id"]: m for m in base["molecules"]}
all_reactions = {r["reaction_id"]: r for r in base["reactions"]}
bundle = st.empty(dataset["documents"])
selected_reactions = {}
selected_molecules = {}


def mid(key):
    return f"{DOC}:mol_{key}"


def rid(key):
    return f"{DOC}:{key}"


def reaction(key):
    full = rid(key)
    if full not in selected_reactions:
        selected_reactions[full] = copy.deepcopy(all_reactions[full])
    return selected_reactions[full]


def molecule(key):
    full = mid(key)
    if full not in selected_molecules:
        selected_molecules[full] = copy.deepcopy(all_molecules[full])
    return selected_molecules[full]


def extend_evidence(obj, evidence):
    obj["evidence_ids"] = uniq(obj.get("evidence_ids", []) + list(evidence))


def touch_molecule(key, evidence):
    m = molecule(key)
    extend_evidence(m["occurrences"][0], evidence)
    return m


def charge(participant, field, raw, normalized, evidence):
    participant["charges"][field] = f(raw, normalized, evidence)


def part(r, participant_id):
    return next(p for p in r["participants"] if p["participant_id"] == participant_id)


def step(step_id, evidence, operation, participant_ids=(), stage="reaction",
         temperature=None, time=None, atmosphere=None):
    s = st.new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(participant_ids)
    s["operation"] = f(operation, operation, evidence)
    if temperature:
        s["temperature"] = f(temperature[0], temperature[1], evidence)
    if time:
        s["time"] = f(time[0], time[1], evidence)
    if atmosphere:
        s["atmosphere"] = f(atmosphere, atmosphere, evidence)
    return s


def amount_outcome(outcome_id, target, raw, value, evidence, scope="isolated product"):
    o = st.new_outcome(outcome_id, "product_amount", [target], evidence)
    o["value"] = f(raw, value, evidence)
    o["scope"] = f(scope, scope, evidence)
    return o


def set_isolated_yield(outcome, raw, number, evidence, scope=None):
    outcome["value"] = f(raw, q(number, "%", raw), evidence)
    outcome["yield_type"] = f("isolated yield", "isolated yield", evidence)
    outcome["measurement_basis"] = f("isolated after silica-gel chromatography",
                                      "isolated", evidence)
    if scope:
        outcome["scope"] = f(scope, scope, evidence)


def enrich(r, page, location, procedure, steps, notes=(), unresolved=()):
    evidence = ev(page)
    extend_evidence(r, evidence)
    old_pages = r["location"].get("normalized_value", {}).get("pages", [])
    r["location"] = f(location, {"pages": sorted(set(old_pages + [page]))}, evidence)
    r["procedure"]["original_text"] = f(procedure, procedure, evidence)
    r["steps"] = steps
    note_facts = [f(note, note, evidence) if isinstance(note, str) else note for note in notes]
    r["notes"] = uniq(r.get("notes", []) + note_facts)
    old_unresolved = [x for x in r["review"].get("unresolved", [])
                      if "measurement basis" not in x.lower() and "exact charges" not in x.lower()]
    r["review"] = {"status": "needs_review", "method": "model", "reviewer": None,
                   "checked_fields": [], "unresolved": uniq(old_unresolved + list(unresolved))}
    r["completeness"] = "partial"


def add_scope_details(letter, page, mass, yield_pct, eluent):
    r = reaction(f"scheme3_scope_2{letter}_to_3{letter}")
    evidence = ev(page)
    target = mid(f"3{letter}")
    y = next(o for o in r["outcomes"] if o["kind"] == "yield")
    set_isolated_yield(y, f"{yield_pct}%", yield_pct, evidence)
    r["outcomes"].append(amount_outcome("isolated_mass", target, f"{mass} mg",
                                         q(mass, "mg", f"{mass} mg"), evidence))
    procedure = (f"Tf2NH-catalyzed aza-Ferrier reaction furnished 3{letter} ({mass} mg, "
                 f"{yield_pct}%); silica-gel chromatography used n-hexane/EtOAc {eluent}.")
    steps = [
        step("reaction", evidence, f"stir 2{letter} with 10 mol% Tf2NH in CH2Cl2",
             ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")),
             time=("1 h", q(1, "h", "1 h"))),
        step("purification", evidence, f"silica-gel chromatography, n-hexane/EtOAc {eluent}",
             stage="purification"),
    ]
    enrich(r, page, f"pages 3 and {page}, Scheme 3 and compound 3{letter} paragraph",
           procedure, steps,
           unresolved=["Run-specific starting-material and catalyst charges are not restated for this scope example."])
    touch_molecule(f"3{letter}", evidence)


# Page 7: overall preparation of 1l by the explicitly cross-referenced three-step 1k route.
e67 = evs(6, 7)
m2pip = st.new_molecule(mid("piperidinone"), DOC, ev(7))
m2pip["identity_status"] = "source_identified"
m2pip["occurrences"][0]["compound_label"] = f("2-piperidone", "2-piperidone", ev(7))
m2pip["names"]["original"] = f("2-piperidone", "2-piperidone", ev(7))
m2pip["names"]["english"] = f("2-piperidone", "2-piperidone", ev(7))
m2pip["representation_type"] = "smiles"
m2pip["smiles"] = f("2-piperidone", "O=C1CCCCN1", ev(7), method="structure_derivation",
                     derivation="Candidate SMILES derived from the explicit compound name; parsing is diagnostic only.")
diag = st.stereo_inspect("O=C1CCCCN1")
m2pip["stereochemistry"]["atom_centers"] = f("standard stereo diagnostic", diag["atom_centers"], ev(7),
    method="structure_derivation", derivation="RDKit zero-based stereocenter diagnostic.")
m2pip["stereochemistry"]["double_bonds"] = f("standard stereo diagnostic", [], ev(7),
    method="structure_derivation", derivation="The carbonyl has no E/Z assignment.")
m2pip["stereochemistry"]["axial"] = f("source drawing and graph", [], ev(7), method="structure_derivation",
    derivation="No stereogenic axis is present in the candidate graph.")
for key, value in (("chirality_status", "achiral"), ("geometry_status", "not_applicable"),
                   ("representation_status", "not_applicable")):
    m2pip["stereochemistry"][key] = f("source name plus standard stereo diagnostic", value, ev(7),
        method="structure_derivation", derivation="Candidate graph contains no unassigned stereogenic element.")
m2pip["stereochemistry"]["identity_key"] = f("standard stereo diagnostic", diag["stereo_graph_key"], ev(7),
    method="structure_derivation", derivation="Hash of the canonical candidate graph; not a sample identity.")
selected_molecules[m2pip["molecule_id"]] = m2pip

prep1l = st.new_reaction(rid("prep_1l_overall_from_2_piperidone"), DOC, e67)
prep1l["record_type"] = "specific_experiment"
prep1l["experiment_categories"] = f(["other"], ["other"], e67)
prep1l["location"] = f("page 7 compound 1l paragraph, referring to page 6 procedure for 1k",
                       {"pages": [6, 7]}, e67)
prep1l["reaction_name"] = f("Preparation of 1l", "Preparation of 1l", e67)
prep1l["transformation_type"] = f("three-step preparation of N,O-propargyl acetal 1l",
    "three-step preparation of N,O-propargyl acetal", e67, method="inference",
    derivation="The page-7 paragraph explicitly invokes the three-step procedure for 1k with 2-piperidone substituted in Step 1.")
participant_specs = [
    ("lactam", "piperidinone", "reactant"), ("nbuli", "nbuli", "reagent"),
    ("boc2o", "boc2o", "reagent"), ("redal", "redal", "reagent"),
    ("propargyl", "propargyl_alcohol", "reactant"), ("sc", "sc_otf3", "catalyst"),
    ("sieves", "molecular_sieves", "additive"), ("thf", "thf", "solvent"),
    ("ch2cl2", "ch2cl2", "solvent"), ("product", "1l", "product"),
]
prep1l["participants"] = [st.new_participant(pid, mid(mkey), role, e67) for pid, mkey, role in participant_specs]
prep1l["steps"] = [
    step("lactam_deprotonation", e67, "apply 1k Step 1: treat 2-piperidone in THF with nBuLi", ["lactam", "nbuli", "thf"], temperature=("0 °C", q(0, "°C", "0 °C")), time=("30 min", q(30, "min", "30 min")), atmosphere="Ar"),
    step("boc_addition", e67, "add (Boc)2O", ["boc2o"], temperature=("0 °C", q(0, "°C", "0 °C"))),
    step("boc_room_temperature_hold", e67, "continue the Boc-protection step", ["lactam", "nbuli", "boc2o"], temperature=("room temperature", "room temperature"), time=("2 h", q(2, "h", "2 h"))),
    step("lactam_reduction", e67, "apply 1k Step 2: add Red-Al to the Boc lactam in CH2Cl2 and hold", ["redal", "ch2cl2"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("2 h", q(2, "h", "2 h")), atmosphere="Ar"),
    step("acetalization", e67, "apply 1k Step 3: Sc(OTf)3/propargyl alcohol/3 Å molecular-sieve acetalization", ["propargyl", "sc", "sieves", "ch2cl2"], temperature=("room temperature", "room temperature"), time=("18 h", q(18, "h", "18 h")), atmosphere="Ar"),
    step("purification", ev(7), "silica-gel chromatography, n-hexane/EtOAc 10:1 to 7:1", stage="purification"),
]
prep1l["procedure"]["original_text"] = f(
    "1l was prepared by the same procedure as 1k, using 2-piperidone instead of 2-pyrrolidinone in Step 1; chromatography gave 309 mg (57% overall yield).",
    "Three-step procedure inherited only where explicitly cross-referenced; 309 mg, 57% overall from 2-piperidone.", e67)
prep1l["procedure"]["general_procedure_reference"] = f("same procedure with 1k", rid("prep_11k_boc_protection"), e67)
prep1l["procedure"]["modifications"] = f("2-piperidone instead of 2-pyrrolidinone in Step 1",
    "replace 2-pyrrolidinone with 2-piperidone in Step 1", ev(7))
prep1l["outcomes"] = [
    amount_outcome("isolated_mass", mid("1l"), "309 mg", q(309, "mg", "309 mg"), ev(7), "isolated colorless crystals"),
    st.new_outcome("overall_yield", "yield", [mid("1l")], ev(7)),
]
set_isolated_yield(prep1l["outcomes"][1], "57% overall yield", 57, ev(7), "2-piperidone to 1l over three steps")
prep1l["classification"]["standard_id"] = "bond-formation-20261016"
formed = [("b1", "C", "N", "C-N"), ("b2", "C", "H", "C-H"),
          ("b3", "C", "O", "C-O")]
prep1l["classification"]["bond_changes"] = [f("partial three-step structure comparison",
    {"change_id": cid, "atoms": [{"atom_map": None, "element": a}, {"atom_map": None, "element": b}],
     "before_order": 0, "after_order": 1, "change_type": "formed" if b != "H" else "hydrogen_gain",
     **({"hydrogen_count": 1} if b == "H" else {})}, e67, method="structure_derivation",
    derivation="Partial overall-route bond assignment from the explicit route reference; no full atom mapping asserted.")
    for cid, a, b, _ in formed]
tax = {"C-N": 5, "C-H": 1, "C-O": 3}
prep1l["classification"]["labels"] = f("partial three-step source comparison",
    [{"category_number": tax[code], "code": code, "supporting_change_ids": [cid]} for cid, _, _, code in formed],
    e67, method="structure_derivation", derivation="Multilabel classification preserves the distinct bonds formed across the reported three-step span.")
prep1l["classification"]["coverage"] = f("partial overall-route comparison", "partial", e67,
    method="inference", derivation="Intermediates and complete atom mapping are not supplied for the 1l run.")
prep1l["classification"]["limitations"] = ["Overall three-step comparison; individual intermediate structures and omitted coproducts are not fully mapped."]
prep1l["review"] = {"status": "needs_review", "method": "model", "reviewer": None, "checked_fields": [],
    "unresolved": ["Run-specific charges for the 1l sequence are not restated; only the explicit substitution and inherited 1k operations are retained.",
                   "The 57% value is an overall three-step yield from 2-piperidone, not the yield of Step 3 alone."]}
prep1l["completeness"] = "partial"
selected_reactions[prep1l["reaction_id"]] = prep1l
touch_molecule("1l", ev(7))

# Page 7: 1l -> 2l exact isolated amount/yield.
r = reaction("scheme2_isomerization_1l_to_2l")
set_isolated_yield(next(o for o in r["outcomes"] if o["kind"] == "yield"), "55% yield", 55, ev(7))
r["outcomes"].append(amount_outcome("isolated_mass", mid("2l"), "166 mg", q(166, "mg", "166 mg"), ev(7)))
enrich(r, 7, "pages 2 and 7, Scheme 2 and compound 2l paragraph",
       "2l was prepared by the same procedure as 2a; chromatography (n-hexane/EtOAc 15:1) gave 166 mg, 55%.",
       [step("isomerization", ev(7), "apply the procedure for 2a to 1l", ["substrate", "triton", "dmso"],
             temperature=("room temperature", "room temperature"), time=("16 h", q(16, "h", "16 h")), atmosphere="Ar"),
        step("purification", ev(7), "silica-gel chromatography, n-hexane/EtOAc 15:1", stage="purification")],
       unresolved=["Run-specific 1l and Triton B charges are not restated in the 2l paragraph."])
touch_molecule("2l", ev(7))

# Page 7 representative entry 6, with the Table 1/procedure scale disagreement retained.
r = reaction("table1_entry_6")
charge(part(r, "substrate"), "mass", "160 mg", q(160, "mg", "160 mg"), ev(7))
charge(part(r, "substrate"), "amount", "0.529 mmol", q(0.529, "mmol", "0.529 mmol"), ev(7))
charge(part(r, "catalyst"), "amount", "0.053 mmol Tf2NH in the charged 0.53 mL aliquot", q(0.053, "mmol", "0.053 mmol"), ev(7))
charge(part(r, "catalyst"), "volume", "0.53 mL aliquot", q(0.53, "mL", "0.53 mL"), ev(7))
charge(part(r, "catalyst"), "mol_percent", "10 mol%", q(10, "mol%", "10 mol%"), evs(2, 7))
charge(part(r, "solvent"), "volume", "5.9 mL reaction solution (2.7 mL dilution plus 3.2 mL substrate solution)",
       q(5.9, "mL", "2.7 mL + 3.2 mL"), ev(7))
r["reaction_scale"] = f("Table 1: approximately 0.3 mmol; experimental paragraph: 0.529 mmol", None, evs(2, 7),
    status="conflict", alternatives=[f("approximately 0.3 mmol", {"reference_molecule_id": mid("2a"), "amount": q(0.3, "mmol", "approximately 0.3 mmol", "approximately")}, ev(2)),
                                     f("0.529 mmol", {"reference_molecule_id": mid("2a"), "amount": q(0.529, "mmol", "0.529 mmol")}, ev(7))])
r["solvent_conditions"]["total_volume"] = f("2.7 mL dilution plus 3.2 mL substrate solution after the 0.53 mL aliquot was charged",
    q(5.9, "mL", "2.7 mL + 3.2 mL"), ev(7), method="inference",
    derivation="Sum of the explicitly stated reaction-solution volumes; the 0.9 mL stock-preparation volume is not added again.")
set_isolated_yield(next(o for o in r["outcomes"] if o["kind"] == "yield"), "93% yield", 93, ev(7))
r["outcomes"].append(amount_outcome("isolated_mass", mid("3a"), "149 mg", q(149, "mg", "149 mg"), ev(7)))
entry6_steps = [
    step("catalyst_stock", ev(7), "dissolve Tf2NH (26 mg, 0.09 mmol) in CH2Cl2 (0.9 mL) at room temperature under Ar", ["catalyst", "solvent"], temperature=("room temperature", "room temperature"), atmosphere="Ar"),
    step("aliquot_dilution", ev(7), "charge 0.53 mL of stock solution (0.053 mmol Tf2NH) and dilute with CH2Cl2 (2.7 mL)", ["catalyst", "solvent"]),
    step("low_temperature_addition", ev(7), "add catalyst solution to 2a (160 mg, 0.529 mmol) in CH2Cl2 (3.2 mL)", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), atmosphere="Ar"),
    step("low_temperature_hold", ev(7), "stir at the same temperature", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("1 h", q(1, "h", "1 h")), atmosphere="Ar"),
    step("workup", ev(7), "pour into saturated aqueous NaHCO3, extract with CH2Cl2, wash with brine, dry over Na2SO4 and evaporate", stage="workup"),
    step("purification", ev(7), "silica-gel chromatography, n-hexane/EtOAc 3:1", stage="purification"),
]
enrich(r, 7, "pages 2 and 7, Table 1 entry 6 and representative procedure",
       "Tf2NH stock solution was aliquoted and diluted, added to 2a in CH2Cl2 at −78 °C under Ar, held 1 h, worked up and chromatographed to give 3a (149 mg, 93%).",
       entry6_steps,
       notes=["Source conflict retained: Table 1 describes approximately 0.3 mmol of 2a, whereas the representative procedure uses 0.529 mmol."],
       unresolved=["The Table 1/procedure scale discrepancy is unresolved and neither value was silently preferred."])
touch_molecule("2a", ev(7)); touch_molecule("3a", ev(7)); touch_molecule("tf2nh", ev(7)); touch_molecule("ch2cl2", ev(7))

# Page 7 scale-up entry 7; preserve 1.25 vs 1.26 mmol.
r = reaction("table1_entry_7")
charge(part(r, "substrate"), "amount", "Table 1: 1.25 mmol; experimental paragraph: 1.26 mmol", None, evs(2, 7))
part(r, "substrate")["charges"]["amount"]["data_status"] = "conflict"
part(r, "substrate")["charges"]["amount"]["alternatives"] = [f("1.25 mmol", q(1.25, "mmol", "1.25 mmol"), ev(2)), f("1.26 mmol", q(1.26, "mmol", "1.26 mmol"), ev(7))]
r["reaction_scale"] = f("Table 1: 1.25 mmol; experimental paragraph: 1.26 mmol", None, evs(2, 7), status="conflict",
    alternatives=[f("1.25 mmol", {"reference_molecule_id": mid("2a"), "amount": q(1.25, "mmol", "1.25 mmol")}, ev(2)),
                  f("1.26 mmol", {"reference_molecule_id": mid("2a"), "amount": q(1.26, "mmol", "1.26 mmol")}, ev(7))])
set_isolated_yield(next(o for o in r["outcomes"] if o["kind"] == "yield"), "93% yield", 93, ev(7))
r["outcomes"].append(amount_outcome("isolated_mass", mid("3a"), "354 mg", q(354, "mg", "354 mg"), ev(7)))
enrich(r, 7, "pages 2 and 7, Table 1 entry 7 and scale-up note",
       "The over-1 mmol experiment gave 3a (354 mg, 93%); Table 1 states 1.25 mmol of 2a and the experimental paragraph states 1.26 mmol.",
       [step("reaction", evs(2, 7), "apply the optimized 10 mol% Tf2NH aza-Ferrier conditions to the scale-up", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("1 h", q(1, "h", "1 h"))),
        step("purification", ev(7), "isolate 3a after chromatography", stage="purification")],
       notes=["Source conflict retained: 1.25 mmol in Table 1 versus 1.26 mmol in the experimental paragraph."],
       unresolved=["Detailed run-specific solvent volume and catalyst mass are not restated for entry 7.", "The 1.25/1.26 mmol discrepancy remains unresolved."])

# Scope products characterized on pages 7 and 8.
for args in [("b", 7, 132, 91, "5:1"), ("c", 7, 60.0, 86, "5:1"),
             ("d", 7, 73.8, 91, "5:1 to 3:1"), ("e", 7, 48.9, 89, "4:1 to 2:1"),
             ("f", 7, 60.1, 75, "3:1"), ("g", 7, 86.1, 69, "3:1"),
             ("h", 7, 128, 89, "5:1 to 3:1"), ("i", 8, 81.4, 81, "7:1 to 5:1"),
             ("j", 8, 3.0, 2, "5:1 to 3:1"), ("k", 8, 178, 94, "5:1 to 3:1"),
             ("l", 8, 154, 93, "5:1 to 3:1")]:
    add_scope_details(*args)

# Separate scale experiments for 2e and 2i.
for letter, page, mmol, mass, ypct in [("e", 7, 1.49, 279, 88), ("i", 8, 1.14, 196, 86)]:
    r = reaction(f"scheme3_scaleup_2{letter}_to_3{letter}")
    charge(part(r, "substrate"), "amount", f"{mmol} mmol", q(mmol, "mmol", f"{mmol} mmol"), ev(page))
    r["reaction_scale"] = f(f"{mmol} mmol 2{letter}", {"reference_molecule_id": mid(f"2{letter}"), "amount": q(mmol, "mmol", f"{mmol} mmol")}, ev(page))
    set_isolated_yield(next(o for o in r["outcomes"] if o["kind"] == "yield"), f"{ypct}% yield", ypct, ev(page))
    r["outcomes"].append(amount_outcome("isolated_mass", mid(f"3{letter}"), f"{mass} mg", q(mass, "mg", f"{mass} mg"), ev(page)))
    enrich(r, page, f"pages 3 and {page}, Scheme 3 footnote and compound 3{letter} paragraph",
           f"Over-1 mmol run from {mmol} mmol of 2{letter} furnished 3{letter} ({mass} mg, {ypct}%).",
           [step("reaction", evs(3, page), f"apply optimized aza-Ferrier conditions to {mmol} mmol of 2{letter}", ["substrate", "catalyst", "solvent"], temperature=("−78 °C", q(-78, "°C", "−78 °C")), time=("1 h", q(1, "h", "1 h"))),
            step("isolation", ev(page), f"isolate 3{letter} after chromatography", stage="purification")],
           unresolved=["Run-specific catalyst mass and solvent volume are not reported for the scale experiment."])

# Page 8 exact follow-up transformations.
def exact_reaction(key, page, charges, steps, procedure, product_key, mass, ypct, eluent, unresolved=()):
    r = reaction(key); evidence = ev(page)
    for pid, fields in charges.items():
        for field, (raw, norm) in fields.items(): charge(part(r, pid), field, raw, norm, evidence)
    set_isolated_yield(next(o for o in r["outcomes"] if o["kind"] == "yield"), f"{ypct}% yield", ypct, evidence)
    r["outcomes"].append(amount_outcome("isolated_mass", mid(product_key), f"{mass} mg", q(mass, "mg", f"{mass} mg"), evidence))
    steps = steps + [step("purification", evidence, f"silica-gel chromatography, n-hexane/EtOAc {eluent}", stage="purification")]
    enrich(r, page, f"pages 3 and {page}, Scheme 4 and experimental paragraph", procedure, steps, unresolved=unresolved)
    for p in r["participants"]: touch_molecule(p["molecule_id"].split(":mol_")[-1], evidence)
    return r

r = exact_reaction("scheme4a_wittig_3e_to_4", 8, {
    "substrate": {"mass": ("106 mg", q(106, "mg", "106 mg")), "amount": ("0.497 mmol", q(0.497, "mmol", "0.497 mmol"))},
    "phosphonium": {"mass": ("404 mg", q(404, "mg", "404 mg")), "amount": ("0.992 mmol", q(0.992, "mmol", "0.992 mmol"))},
    "base": {"mass": ("112 mg", q(112, "mg", "112 mg")), "amount": ("1.00 mmol", q(1.00, "mmol", "1.00 mmol"))},
    "solvent": {"volume": ("7.5 mL total (5 mL plus 2.5 mL)", q(7.5, "mL", "5 mL + 2.5 mL"))}}, [
    step("ylide_formation", ev(8), "treat methyltriphenylphosphonium iodide in Et2O (5 mL) with tBuOK", ["phosphonium", "base", "solvent"], temperature=("room temperature", "room temperature"), time=("15 min", q(15, "min", "15 min")), atmosphere="Ar"),
    step("cool_and_add", ev(8), "cool to 0 °C and add 3e in Et2O (2.5 mL)", ["substrate", "solvent"], temperature=("0 °C", q(0, "°C", "0 °C"))),
    step("room_temperature_hold", ev(8), "stir after substrate addition", ["substrate", "phosphonium", "base"], temperature=("room temperature", "room temperature"), time=("18 h", q(18, "h", "18 h"))),
    step("workup", ev(8), "quench with saturated NH4Cl, extract with EtOAc, wash brine, dry Na2SO4 and evaporate", stage="workup")],
    "Methyltriphenylphosphonium iodide and tBuOK were stirred in Et2O 15 min at room temperature under Ar; 3e was added at 0 °C, the mixture was stirred 18 h at room temperature, worked up and chromatographed to give 4 (82.3 mg, 78%).",
    "4", 82.3, 78, "10:1")
r["reaction_scale"] = f("0.497 mmol 3e", {"reference_molecule_id": mid("3e"), "amount": q(0.497, "mmol", "0.497 mmol")}, ev(8))

r = exact_reaction("scheme4a_diels_alder_4_to_5", 8, {
    "diene": {"mass": ("61.7 mg", q(61.7, "mg", "61.7 mg")), "amount": ("0.292 mmol", q(0.292, "mmol", "0.292 mmol"))},
    "dienophile": {"mass": ("46 mg", q(46, "mg", "46 mg")), "amount": ("0.29 mmol", q(0.29, "mmol", "0.29 mmol"))},
    "solvent": {"volume": ("3 mL", q(3, "mL", "3 mL"))}}, [
    step("reflux", ev(8), "reflux 4 and 1,4-naphthoquinone in toluene", ["diene", "dienophile", "solvent"], temperature=("reflux", "reflux"), time=("1 day", q(1, "day", "1 day")), atmosphere="Ar"),
    step("cool_and_concentrate", ev(8), "cool to room temperature and concentrate by evaporation", stage="workup", temperature=("room temperature", "room temperature"))],
    "4 and 1,4-naphthoquinone were refluxed in toluene for 1 day under Ar; cooling, evaporation and chromatography gave 5 (91.0 mg, 85%), reported as a 1:1 diastereomer mixture.",
    "5", 91.0, 85, "3:1 to 2:1", unresolved=["Diastereomer identities/order and selectivity measurement method are not reported."])
r["reaction_scale"] = f("0.292 mmol 4", {"reference_molecule_id": mid("4"), "amount": q(0.292, "mmol", "0.292 mmol")}, ev(8))
dr = next(o for o in r["outcomes"] if o["kind"] == "selectivity")
dr["value"] = f("Scheme 4: 5/5 dr; experimental paragraph: 1:1 mixture of diastereomers", {"metric": "dr", "ratio": [1, 1], "order": "unassigned", "configurations": ["unknown", "unknown"]}, evs(3, 8))
dr["measurement_basis"] = f(None, None, evs(3, 8), status="unreported")
dr["scope"] = f("diastereomers of compound 5", "diastereomers of compound 5", ev(8))
# Diels-Alder creates two C-C sigma bonds; retain both changes under one label.
r["classification"]["bond_changes"] = [f("source Scheme 4 Diels-Alder drawing", {"change_id": f"b{i}", "atoms": [{"atom_map": None, "element": "C"}, {"atom_map": None, "element": "C"}], "before_order": 0, "after_order": 1, "change_type": "formed"}, evs(3, 8), method="structure_derivation", derivation="Two distinct cycloaddition C-C bonds are present in the source product drawing; no full atom mapping asserted.") for i in (1, 2)]
r["classification"]["labels"] = f("two C-C bonds formed in the [4+2] cycloaddition", [{"category_number": 2, "code": "C-C", "supporting_change_ids": ["b1", "b2"]}], evs(3, 8), method="structure_derivation", derivation="Both newly formed sigma bonds support the C-C category.")
m5 = touch_molecule("5", ev(8)); m5["stereochemistry"]["mixture"] = f("1:1 mixture of diastereomers", "1:1 mixture; configurations/order unassigned", ev(8)); m5["stereochemistry"]["chirality_status"] = f("1:1 mixture of diastereomers", "mixture", ev(8)); m5["stereochemistry"]["geometry_status"] = f("relative stereochemistry not assigned", "unresolved", ev(8)); m5["stereochemistry"]["representation_status"] = f("diastereomers not individually encoded", "not_encoded", ev(8), method="inference", derivation="The source reports a mixture without individual configuration assignments.")

r = exact_reaction("scheme4b_reduction_3l_to_6", 8, {
    "substrate": {"mass": ("45.5 mg", q(45.5, "mg", "45.5 mg")), "amount": ("0.190 mmol", q(0.190, "mmol", "0.190 mmol"))},
    "cecl3": {"mass": ("85 mg", q(85, "mg", "85 mg")), "amount": ("0.23 mmol", q(0.23, "mmol", "0.23 mmol"))},
    "nabh4": {"mass": ("12 mg", q(12, "mg", "12 mg")), "amount": ("0.32 mmol", q(0.32, "mmol", "0.32 mmol"))},
    "solvent": {"volume": ("0.95 mL", q(0.95, "mL", "0.95 mL"))}}, [
    step("low_temperature_addition", ev(8), "treat 3l and CeCl3·7H2O in MeOH with NaBH4", ["substrate", "cecl3", "nabh4", "solvent"], temperature=("0 °C", q(0, "°C", "0 °C")), atmosphere="Ar"),
    step("low_temperature_hold", ev(8), "stir at 0 °C", ["substrate", "cecl3", "nabh4"], temperature=("0 °C", q(0, "°C", "0 °C")), time=("30 min", q(30, "min", "30 min")), atmosphere="Ar"),
    step("workup", ev(8), "quench with saturated NH4Cl, extract EtOAc, wash brine, dry Na2SO4 and evaporate", stage="workup")],
    "3l and CeCl3·7H2O in MeOH were treated with NaBH4 at 0 °C and stirred 30 min under Ar; workup and chromatography gave 6 (40.4 mg, 88%).",
    "6", 40.4, 88, "2:1")
r["reaction_scale"] = f("0.190 mmol 3l", {"reference_molecule_id": mid("3l"), "amount": q(0.190, "mmol", "0.190 mmol")}, ev(8))

r = exact_reaction("scheme4b_cyclization_6_to_7", 8, {
    "substrate": {"mass": ("40.2 mg", q(40.2, "mg", "40.2 mg")), "amount": ("0.167 mmol", q(0.167, "mmol", "0.167 mmol"))},
    "nah": {"mass": ("8 mg of 60 wt% dispersion in oil", q(8, "mg", "8 mg")), "amount": ("0.2 mmol", q(0.2, "mmol", "0.2 mmol"))},
    "solvent": {"volume": ("0.85 mL", q(0.85, "mL", "0.85 mL"))}}, [
    step("low_temperature_addition", ev(8), "treat 6 in THF with NaH (60 wt% in oil)", ["substrate", "nah", "solvent"], temperature=("0 °C", q(0, "°C", "0 °C"))),
    step("low_temperature_hold", ev(8), "stir at 0 °C", ["substrate", "nah"], temperature=("0 °C", q(0, "°C", "0 °C")), time=("1 h", q(1, "h", "1 h"))),
    step("room_temperature_hold", ev(8), "continue stirring at room temperature", ["substrate", "nah"], temperature=("room temperature", "room temperature"), time=("2 h", q(2, "h", "2 h"))),
    step("workup", ev(8), "quench with saturated NH4Cl at 0 °C, extract EtOAc, wash brine, dry Na2SO4 and evaporate", stage="workup", temperature=("0 °C", q(0, "°C", "0 °C")))],
    "6 in THF was treated with NaH at 0 °C, held 1 h at 0 °C and 2 h at room temperature; workup and chromatography gave 7 (20.6 mg, 74%).",
    "7", 20.6, 74, "1:2", unresolved=["Compound 7 connectivity/stereochemistry remains unencoded pending expert review of the source drawing."])
r["reaction_scale"] = f("0.167 mmol 6", {"reference_molecule_id": mid("6"), "amount": q(0.167, "mmol", "0.167 mmol")}, ev(8))

# Page 9 Michael addition.
r = exact_reaction("scheme4c_michael_3i_to_8", 9, {
    "substrate": {"mass": ("61.6 mg", q(61.6, "mg", "61.6 mg")), "amount": ("0.309 mmol", q(0.309, "mmol", "0.309 mmol"))},
    "malonate": {"volume": ("47 μL", q(47, "μL", "47 μL")), "amount": ("0.31 mmol", q(0.31, "mmol", "0.31 mmol"))},
    "catalyst": {"mass": ("2 mg", q(2, "mg", "2 mg")), "amount": ("0.03 mmol", q(0.03, "mmol", "0.03 mmol"))},
    "solvent": {"volume": ("3.0 mL total (1.5 mL plus 1.5 mL)", q(3.0, "mL", "1.5 mL + 1.5 mL"))}}, [
    step("malonate_activation", ev(9), "treat diethyl malonate in EtOH (1.5 mL) with EtONa", ["malonate", "catalyst", "solvent"], temperature=("room temperature", "room temperature"), time=("5 min", q(5, "min", "5 min")), atmosphere="Ar"),
    step("substrate_addition", ev(9), "add a solution of 3i in EtOH (1.5 mL)", ["substrate", "solvent"], temperature=("room temperature", "room temperature")),
    step("room_temperature_hold", ev(9), "stir the combined mixture", ["substrate", "malonate", "catalyst"], temperature=("room temperature", "room temperature"), time=("6 h", q(6, "h", "6 h"))),
    step("workup", ev(9), "quench with saturated NH4Cl, extract EtOAc, wash brine, dry Na2SO4 and evaporate", stage="workup")],
    "Diethyl malonate and EtONa were stirred in EtOH for 5 min under Ar; 3i in EtOH was added and stirred 6 h at room temperature, then worked up and chromatographed to give 8 (87.6 mg, 79%).",
    "8", 87.6, 79, "4:1 to 2:1")
r["reaction_scale"] = f("0.309 mmol 3i", {"reference_molecule_id": mid("3i"), "amount": q(0.309, "mmol", "0.309 mmol")}, ev(9))

# Characterization and reference material are source context, not experiments.
bundle["source_context"] = [
    {"context_id": f"{DOC}:context_characterization_pages7_9", "kind": "characterization",
     "evidence_ids": evs(7, 8, 9),
     "text": "Physical descriptions, melting points, NMR, IR and HRMS for 1l, 2l, 3a-l and 4-8; mass-spectral peaks and analytical data were not promoted to synthetic experiments. 13C NMR for 3j was explicitly not performed because of the small amount."},
    {"context_id": f"{DOC}:context_associated_content_page9", "kind": "associated_content_notice",
     "evidence_ids": ev(9),
     "text": "The paper states that Supporting Information is available online, but no SI file was supplied in this workspace; no SI association was inferred."},
    {"context_id": f"{DOC}:context_references_pages9_10", "kind": "references",
     "evidence_ids": evs(9, 10),
     "text": "Reference-list material only; cited literature and publication metadata were not treated as current-paper experiments."},
]

# Pull every registry molecule referenced by the selected records into the atomic merge bundle.
# The standalone submission must also contain the generic scheme targeted by
# general_scheme_ids; it is unchanged and is not a newly claimed page-7--10 result.
selected_reactions[rid("scheme1d_current_generic")] = copy.deepcopy(all_reactions[rid("scheme1d_current_generic")])
needed = set(selected_molecules)
for r in selected_reactions.values():
    needed.update(p["molecule_id"] for p in r["participants"])
    for o in r["outcomes"]:
        needed.update(o.get("target_molecule_ids", []))
for full in sorted(needed):
    if full not in selected_molecules:
        selected_molecules[full] = copy.deepcopy(all_molecules[full])

bundle["molecules"] = list(selected_molecules.values())
bundle["reactions"] = list(selected_reactions.values())
OUT_RAW.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
OUT_SUBMITTED.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

regions = {
    7: [
        {"region": "top-left continuation of 2k characterization", "finding": "Characterization only; retained in source_context, not a synthesis record.", "evidence_ids": ev(7)},
        {"region": "left column compounds 1l and 2l", "finding": "Three-step overall preparation of 1l (309 mg, 57% overall) and 1l→2l isomerization (166 mg, 55%); distinct records.", "evidence_ids": ev(7)},
        {"region": "left-to-right representative 2a→3a procedure", "finding": "Exact Tf2NH stock/aliquot, −78 °C addition and 1 h hold, workup, purification, 149 mg/93%; Table 1 scale conflict retained.", "evidence_ids": ev(7)},
        {"region": "product paragraphs 3b–3h", "finding": "Isolated masses, yields and chromatography conditions merged into the corresponding Scheme 3 records; 3e scale-up remains separate.", "evidence_ids": ev(7)},
    ],
    8: [
        {"region": "top continuation of 3h and paragraphs 3i–3l", "finding": "Characterization plus isolated masses/yields/purification for 3i–3l; 3i scale-up remains separate and 3j small-amount/mixture limitation retained.", "evidence_ids": ev(8)},
        {"region": "Wittig reaction 3e→4", "finding": "Exact reagent doses, 15 min ylide formation, 0 °C addition, 18 h room-temperature hold, workup, chromatography, 82.3 mg/78%.", "evidence_ids": ev(8)},
        {"region": "Diels–Alder reaction 4→5", "finding": "Exact doses, toluene reflux for 1 day, workup, 91.0 mg/85%, and source-reported 1:1 unassigned diastereomer mixture.", "evidence_ids": ev(8)},
        {"region": "reduction 3l→6 and cyclization 6→7", "finding": "Exact doses and ordered temperature/time stages; 6→7 retains separate 0 °C 1 h and room-temperature 2 h holds.", "evidence_ids": ev(8)},
    ],
    9: [
        {"region": "top-left tail of compound 7 characterization", "finding": "Analytical characterization only; not a separate synthesis experiment.", "evidence_ids": ev(9)},
        {"region": "Michael addition 3i→8", "finding": "Exact diethyl malonate/EtONa activation, substrate addition, 6 h room-temperature hold, workup, chromatography, 87.6 mg/79%.", "evidence_ids": ev(9)},
        {"region": "Associated Content and references", "finding": "SI availability notice retained without associating absent SI; reference-list entries are not current experiments.", "evidence_ids": ev(9)},
    ],
    10: [
        {"region": "entire page, continuation of references 10–13", "finding": "Reference-list material only; no current-paper experimental arrow, table row, or procedure was present.", "evidence_ids": ev(10)},
    ],
}
page_records = {
    7: [rid("prep_1l_overall_from_2_piperidone"), rid("scheme2_isomerization_1l_to_2l"), rid("table1_entry_6"), rid("table1_entry_7")] +
       [rid(f"scheme3_scope_2{x}_to_3{x}") for x in "bcdefgh"] + [rid("scheme3_scaleup_2e_to_3e")],
    8: [rid(f"scheme3_scope_2{x}_to_3{x}") for x in "hijkl"] + [rid("scheme3_scaleup_2i_to_3i"),
       rid("scheme4a_wittig_3e_to_4"), rid("scheme4a_diels_alder_4_to_5"), rid("scheme4b_reduction_3l_to_6"), rid("scheme4b_cyclization_6_to_7")],
    9: [rid("scheme4c_michael_3i_to_8")],
    10: [],
}
unresolved = {
    7: ["Table 1 versus experimental-procedure scale conflicts (approximately 0.3 vs 0.529 mmol; 1.25 vs 1.26 mmol) remain unresolved.", "Run-specific charges for same-procedure scope preparations are not restated."],
    8: ["Compound 5 diastereomer identities/order and measurement method are unreported; compound 7 structure remains unencoded.", "Scope-record starting-material and catalyst charges are not restated in the product paragraphs."],
    9: ["No SI file was supplied; the availability notice does not establish an SI association.", "Compound 8 stereocenter configuration is not assigned by the source and remains unspecified."],
    10: ["No current-paper experiment is present on this reference-only page."],
}
coverage = {"document_id": DOC, "semantic_accuracy_verified": False,
            "pages": [{"page": p, "state": "attempted", "regions": regions[p],
                       "record_ids": uniq(page_records[p]), "unresolved": unresolved[p]}
                      for p in range(7, 11)]}
OUT_COVERAGE.parent.mkdir(parents=True, exist_ok=True)
OUT_COVERAGE.write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"molecules": len(bundle["molecules"]), "reactions": len(bundle["reactions"]),
                  "raw": str(OUT_RAW), "submitted": str(OUT_SUBMITTED),
                  "coverage": str(OUT_COVERAGE)}, ensure_ascii=False))
