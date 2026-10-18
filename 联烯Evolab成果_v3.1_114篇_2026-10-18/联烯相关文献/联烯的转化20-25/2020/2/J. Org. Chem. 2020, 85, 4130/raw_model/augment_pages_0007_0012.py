#!/usr/bin/env python3
"""Merge visually inspected PDF pages 7--12 into the existing v3.1 bundle.

This script only adds source-derived evidence and resolves no chemistry by fiat.
Where the paper is internally inconsistent, both readings are retained as a
conflict and the reaction remains needs_review.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/current"
DOC = "doc_5909d9a1e9253598"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

from standard_records import fact, new_outcome, new_participant, new_step  # noqa: E402


def ev(page: int) -> list[str]:
    return [f"{DOC}:page{page}:image"]


def q(number, unit, raw_text=None):
    value = {"number": number, "unit": unit}
    if raw_text is not None:
        value["raw_text"] = raw_text
    return value


def F(raw=None, normalized=None, evidence=(), method="original_extraction",
      status=None, alternatives=(), derivation=None):
    return fact(raw=raw, normalized=normalized, evidence=evidence, method=method,
                status=status, alternatives=alternatives, derivation=derivation)


bundle_path = OUT / "submitted.json"
bundle = json.loads(bundle_path.read_text())
# The previous page-1--6 submission predated the staged page-7--12 image
# evidence.  Formal submissions must use the dataset evidence registry exactly.
bundle["documents"] = json.loads((OUT / "dataset.json").read_text())["documents"]
if any(c.get("context_id") == f"{DOC}:context:characterization_pages7_12"
       for c in bundle.get("source_context", [])):
    raise SystemExit("pages 7--12 augmentation is already present; refusing to duplicate participants or steps")
reactions = {r["reaction_id"]: r for r in bundle["reactions"]}
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}


def mid(label: str) -> str:
    return f"{DOC}:mol:{label}"


def rid(suffix: str) -> str:
    return f"{DOC}:{suffix}"


def uniq_extend(target: list, values):
    for value in values:
        if value not in target:
            target.append(value)


def add_page_evidence(r, page):
    uniq_extend(r["evidence_ids"], ev(page))


def add_occurrence(label: str, page: int):
    m = molecules[mid(label)]
    eid = ev(page)
    if not any(o["document_id"] == DOC and o["compound_label"]["normalized_value"] == label
               and o["evidence_ids"] == eid for o in m["occurrences"]):
        m["occurrences"].append({
            "document_id": DOC,
            "compound_label": F(label, label, eid),
            "evidence_ids": eid,
        })


def participant(r, label: str):
    return next(p for p in r["participants"] if p["molecule_id"] == mid(label))


def set_mass(p, number, evidence):
    p["charges"]["mass"] = F(f"{number} mg", q(number, "mg"), evidence)


def set_amount(p, number, evidence, method="original_extraction"):
    p["charges"]["amount"] = F(f"{number} mmol", q(number, "mmol"), evidence, method=method)


def add_product_mass(r, product, number, evidence):
    oid = f"o_product_amount_{product}"
    if any(o["outcome_id"] == oid for o in r["outcomes"]):
        return
    o = new_outcome(oid, "product_amount", [mid(product)], evidence)
    o["product_role"] = F("isolated product", "product", evidence)
    o["value"] = F(f"{number} mg", q(number, "mg"), evidence)
    o["measurement_basis"] = F("isolated product mass", "isolated product mass", evidence)
    r["outcomes"].append(o)


def enrich_yield_evidence(r, evidence):
    for o in r["outcomes"]:
        if o["kind"] == "yield":
            for key in ("value", "yield_type", "measurement_basis", "product_role"):
                uniq_extend(o[key]["evidence_ids"], evidence)


def add_unresolved(r, text):
    if text not in r["review"]["unresolved"]:
        r["review"]["unresolved"].append(text)


def add_note(r, raw, normalized, evidence, *, status=None, alternatives=()):
    candidate = F(raw, normalized, evidence, status=status, alternatives=alternatives)
    if not any(n["raw_value"] == raw for n in r["notes"]):
        r["notes"].append(candidate)


def add_participant(r, label, role, evidence, stage, *, volume_raw=None,
                    volume_number=None):
    base = f"p_{stage}_{label}"
    pid = base
    counter = 2
    ids = {p["participant_id"] for p in r["participants"]}
    while pid in ids:
        pid = f"{base}_{counter}"
        counter += 1
    p = new_participant(pid, mid(label), role, evidence, stage=stage)
    if volume_raw is not None:
        p["charges"]["volume"] = F(volume_raw, q(volume_number, "mL", volume_raw), evidence,
                                           method="inherited_general_condition")
    r["participants"].append(p)
    return p


def add_step(r, step_id, operation, evidence, stage, pids, *, temperature=None,
             time=None, method="inherited_general_condition", prepend=False):
    if any(s["step_id"] == step_id for s in r["steps"]):
        return
    s = new_step(step_id, evidence, stage=stage)
    s["participant_ids"] = list(pids)
    s["operation"] = F(operation, operation, evidence, method=method)
    if temperature is not None:
        raw = "rt (25 °C)" if temperature == 25 else f"{temperature} °C"
        s["temperature"] = F(raw, q(temperature, "°C"), evidence, method=method)
    if time is not None:
        number, unit = time
        s["time"] = F(f"{number} {unit}", q(number, unit), evidence, method=method)
    if prepend:
        r["steps"].insert(0, s)
    else:
        r["steps"].append(s)


def add_scope_workup(r, evidence, eluent, *, water=10, extraction="3 × 5 mL",
                     extraction_total=15, brine_raw=None):
    water_p = add_participant(r, "water", "workup_agent", evidence, "workup",
                              volume_raw=f"{water} mL", volume_number=water)
    etoac_w = add_participant(r, "EtOAc", "workup_agent", evidence, "workup",
                              volume_raw=extraction, volume_number=extraction_total)
    brine_p = add_participant(r, "brine", "workup_agent", evidence, "workup",
                              volume_raw=brine_raw, volume_number=20 if brine_raw else None)
    sulfate_p = add_participant(r, "Na2SO4", "workup_agent", evidence, "workup")
    silica_p = add_participant(r, "silica_gel", "purification_agent", evidence, "purification")
    etoac_p = add_participant(r, "EtOAc", "purification_agent", evidence, "purification")
    hexane_p = add_participant(r, "hexane", "purification_agent", evidence, "purification")
    add_step(r, "workup_quench_extract", f"quench with water ({water} mL); extract the aqueous layer with ethyl acetate ({extraction})",
             evidence, "workup", [water_p["participant_id"], etoac_w["participant_id"]])
    wash = "wash the combined organic layer with brine"
    if brine_raw:
        wash += f" ({brine_raw})"
    add_step(r, "workup_wash_dry_concentrate",
             wash + "; dry over anhydrous Na2SO4; concentrate under reduced pressure",
             evidence, "workup", [brine_p["participant_id"], sulfate_p["participant_id"]])
    add_step(r, "purification_column",
             f"silica gel column chromatography using ethyl acetate/hexane ({eluent})",
             evidence, "purification",
             [silica_p["participant_id"], etoac_p["participant_id"], hexane_p["participant_id"]])


# Individually reported Table 2 products on pages 7--8.  Page 7 begins with
# characterization of 3ae; it is retained as context, not a new experiment.
scope3_pages = {
    "3af": 7, "3ba": 7, "3bb": 7, "3bc": 7, "3bd": 7,
    "3ca": 7, "3cb": 7, "3cc": 7, "3da": 7,
    "3db": 8, "3dc": 8, "3ea": 8,
}
scope3_masses = {
    "3af": (39.4, 63.8, 15.2, 31.8, 46.7),
    "3ba": (42.2, 62.4, 15.2, 31.8, 56.2),
    "3bb": (42.2, 69.6, 15.2, 31.8, 56.5),
    "3bc": (42.2, 70.6, 15.2, 31.8, 63.2),
    "3bd": (42.2, 73.2, 15.2, 31.8, 56.4),
    "3ca": (55.0, 62.4, 15.2, 31.8, 61.0),
    "3cb": (55.0, 69.6, 15.2, 31.8, 59.1),
    "3cc": (55.0, 70.6, 15.2, 31.8, 62.3),
    "3da": (42.2, 62.4, 30.4, 31.8, 52.0),
    "3db": (42.2, 69.6, 30.4, 31.8, 54.3),
    "3dc": (42.2, 70.6, 30.4, 31.8, 58.0),
    "3ea": (49.4, 62.4, 15.2, 31.8, 53.0),
}
for product, page in scope3_pages.items():
    r = reactions[rid(f"scope_table2_{product}")]
    evidence = ev(page) + (ev(8) if product == "3da" else [])
    add_page_evidence(r, page)
    if product == "3da":
        add_page_evidence(r, 8)
    add_occurrence(product, page)
    imine_p, allene_p = r["participants"][0], r["participants"][1]
    dbu_p = participant(r, "DBU")
    carbonate_p = participant(r, "Na2CO3")
    imine_mass, allene_mass, dbu_mass, carbonate_mass, product_mass = scope3_masses[product]
    set_mass(imine_p, imine_mass, evidence)
    set_mass(allene_p, allene_mass, evidence)
    set_mass(dbu_p, dbu_mass, evidence)
    set_mass(carbonate_p, carbonate_mass, evidence)
    add_product_mass(r, product, product_mass, evidence)
    enrich_yield_evidence(r, evidence)
    # The inspected individual entries and General Procedure ii say only
    # "stipulated time"; no numeric duration is attached to these entries.
    r["steps"][0]["time"] = F(evidence=evidence)
    add_unresolved(r, "The inspected individual entry/general procedure says 'stipulated time' but gives no numeric reaction duration.")

# Explicit 3bb table/general-procedure versus individual-prose conflict.
r = reactions[rid("scope_table2_3bb")]
p = participant(r, "2b")
p["charges"]["amount"] = F(
    {"reported_values": ["0.24 mmol", "0.20 mmol"]}, None, ev(3) + ev(6) + ev(7),
    status="conflict",
    alternatives=[
        F("0.24 mmol", q(0.24, "mmol"), ev(3) + ev(6), method="inherited_general_condition"),
        F("0.20 mmol", q(0.20, "mmol"), ev(7)),
    ],
)
add_note(r, "Table 2/General Procedure ii give allenoate 2b as 0.24 mmol, whereas the individual 3bb paragraph prints 0.20 mmol.",
         None, ev(3) + ev(6) + ev(7), status="conflict",
         alternatives=[F("0.24 mmol", q(0.24, "mmol"), ev(3) + ev(6)),
                       F("0.20 mmol", q(0.20, "mmol"), ev(7))])
add_unresolved(r, "Source conflict: 2b amount is 0.24 mmol in the table/general procedure but 0.20 mmol in the individual 3bb paragraph.")


# General Procedures iii and iv and their individual scope paragraphs.
scope4 = {"4aa": (8, 39.4, 62.4, 57.9), "4ba": (8, 42.2, 62.4, 62.5),
          "4ca": (8, 55.0, 62.4, 67.4), "4af": (8, 39.4, 63.8, 52.4)}
scope5 = {
    "5aa": (9, 39.4, 62.4, 51.5), "5ab": (9, 39.4, 69.6, 52.2),
    "5ac": (9, 39.4, 70.6, 58.4), "5ad": (9, 39.4, 73.2, 62.0),
    "5ae": (9, 39.4, 95.8, 59.4), "5af": (9, 39.4, 63.8, 44.0),
    "5ba": (9, 42.2, 62.4, 49.1), "5bb": (9, 42.2, 69.6, 50.7),
    "5bc": (9, 42.2, 70.6, 55.0), "5bd": (10, 42.2, 73.2, 58.0),
    "5ca": (10, 55.0, 62.4, 58.0), "5cb": (10, 55.0, 69.6, 59.0),
    "5cc": (10, 55.0, 70.6, 63.0), "5cd": (10, 55.0, 73.2, 67.0),
    "5da": (10, 42.2, 62.4, 48.0),
}
for product, (page, imine_mass, allene_mass, product_mass) in {**scope4, **scope5}.items():
    r = reactions[rid(f"scope_table3_{product}")]
    gp_page = 8
    evidence = ev(gp_page) + (ev(9) if product == "5bd" else []) + ev(page)
    add_page_evidence(r, page)
    if product == "5bd":
        add_page_evidence(r, 9)
        add_occurrence(product, 9)
    add_page_evidence(r, gp_page)
    add_occurrence(product, page)
    set_mass(r["participants"][0], imine_mass, evidence)
    set_mass(r["participants"][1], allene_mass, evidence)
    set_mass(participant(r, "PPh3"), 10.5, evidence)
    add_product_mass(r, product, product_mass, ev(page))
    enrich_yield_evidence(r, ev(page))
    family4 = product.startswith("4")
    gp_label = "General Procedure iii for compounds 4aa-4ca" if family4 else "General Procedure iv for compounds 5aa-5da"
    r["procedure"]["general_procedure_reference"] = F(gp_label, gp_label.split(" for")[0], ev(8))
    r["procedure"]["addition_order"] = F(
        "charge cyclic sulfamidate imine, PPh3 and 1.0 mL toluene; add allenoate in 1.0 mL toluene gradually over 30 min",
        ["cyclic sulfamidate imine", "PPh3", "1.0 mL toluene", "allenoate in 1.0 mL toluene over 30 min"],
        ev(8), method="inherited_general_condition")
    addition_temp = 25 if family4 else 80
    add_step(r, "allenoate_addition",
             "add δ-acetoxy allenoate solution gradually over 30 min",
             ev(8), "reaction", [r["participants"][1]["participant_id"], participant(r, "toluene")["participant_id"]],
             temperature=addition_temp, time=(30, "min"), prepend=True)
    add_scope_workup(r, ev(8), "30:70" if family4 else "10:90")


# Scale-up experiments: page 10 supplies the complete operations omitted from
# the earlier Scheme 3-only record.
for product in ("3ab", "5ab"):
    r = reactions[rid(f"scaleup_{product}")]
    add_page_evidence(r, 10)
    evidence = ev(3) + ev(10)
    for p in r["participants"]:
        uniq_extend(p["charges"]["amount"]["evidence_ids"], ev(10))
        uniq_extend(p["charges"]["mass"]["evidence_ids"], ev(10))
    enrich_yield_evidence(r, ev(10))
    for o in r["outcomes"]:
        if o["kind"] == "product_amount":
            for key in ("value", "measurement_basis", "product_role"):
                uniq_extend(o[key]["evidence_ids"], ev(10))
    if product == "3ab":
        set_mass(participant(r, "Na2CO3"), 159.0, ev(10))
        participant(r, "DBU")["charges"]["volume"] = F("75 μL", q(75, "μL"), ev(10))
        r["procedure"]["addition_order"] = F(
            "charge 1a, 2b, Na2CO3 and 8.0 mL toluene; then add DBU in 2.0 mL toluene",
            ["1a", "2b", "Na2CO3", "8.0 mL toluene", "DBU in 2.0 mL toluene"], ev(10))
        add_scope_workup(r, ev(10), "10:90", water=30, extraction="3 × 15 mL",
                         extraction_total=45, brine_raw="20 mL")
    else:
        set_mass(participant(r, "PPh3"), 52.5, ev(10))
        r["procedure"]["addition_order"] = F(
            "charge 1a, PPh3 and 5.0 mL toluene; add 2b in 5.0 mL toluene gradually over 45 min at 80 °C",
            ["1a", "PPh3", "5.0 mL toluene", "2b in 5.0 mL toluene over 45 min at 80 °C"], ev(10))
        add_step(r, "allenoate_addition", "add 2b solution gradually over 45 min",
                 ev(10), "reaction", [participant(r, "2b")["participant_id"], participant(r, "toluene")["participant_id"]],
                 temperature=80, time=(45, "min"), method="original_extraction", prepend=True)
        add_scope_workup(r, ev(10), "10:90", water=30, extraction="3 × 5 mL",
                         extraction_total=15, brine_raw="20 mL")
    r["review"]["unresolved"] = [u for u in r["review"]["unresolved"] if "outside the page-1--6" not in u]


# General Procedure vi and Table 4 individual paragraphs.  The general
# procedure says DBU 0.20 mmol, whereas each printed individual example says
# 30.4 mg (0.10 mmol); both values are kept.
scope7 = {
    "7aa": (10, 36.2, 62.4, 45.0), "7ab": (11, 36.2, 69.6, 48.0),
    "7ac": (11, 36.2, 70.6, 51.2), "7ba": (11, 39.0, 62.4, 45.7),
    "7bb": (11, 39.0, 69.6, 46.2), "7bc": (11, 39.0, 70.6, 51.1),
    "7ae": (11, 36.2, 95.8, 52.0), "7af": (11, 36.2, 63.8, 40.0),
}
for product, (page, imine_mass, allene_mass, product_mass) in scope7.items():
    r = reactions[rid(f"scope_table4_{product}")]
    evidence = ev(4) + ev(10) + ev(page)
    add_page_evidence(r, 10)
    add_page_evidence(r, page)
    add_occurrence(product, page)
    set_mass(r["participants"][0], imine_mass, ev(page))
    set_mass(r["participants"][1], allene_mass, ev(page))
    set_mass(participant(r, "DBU"), 30.4, ev(page))
    set_mass(participant(r, "Na2CO3"), 31.8, ev(page))
    participant(r, "DBU")["charges"]["amount"] = F(
        {"reported_values": ["0.20 mmol", "0.10 mmol"]}, None, evidence,
        status="conflict",
        alternatives=[
            F("0.20 mmol", q(0.20, "mmol"), ev(4) + ev(10), method="inherited_general_condition"),
            F("0.10 mmol", q(0.10, "mmol"), ev(page)),
        ])
    add_product_mass(r, product, product_mass, ev(page))
    enrich_yield_evidence(r, ev(page))
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure vi for compounds 7aa-7bc", "General Procedure vi", ev(10))
    r["procedure"]["addition_order"] = F(
        "charge cyclic sulfonyl imine, allenoate, Na2CO3 and 1.5 mL toluene; then add DBU in 0.5 mL toluene",
        ["cyclic sulfonyl imine", "allenoate", "Na2CO3", "1.5 mL toluene", "DBU in 0.5 mL toluene"],
        ev(10), method="inherited_general_condition")
    add_note(r,
             "General Procedure vi gives DBU as 0.20 mmol, while the individual paragraph prints DBU 30.4 mg (0.10 mmol).",
             None, ev(10) + ev(page), status="conflict",
             alternatives=[F("0.20 mmol", q(0.20, "mmol"), ev(10)),
                           F("0.10 mmol", q(0.10, "mmol"), ev(page))])
    add_unresolved(r, "Source conflict: DBU is 0.20 mmol in General Procedure vi but 0.10 mmol in the individual paragraph; the printed 30.4 mg is retained separately.")
    if product != "7aa":
        r["steps"][0]["time"] = F(evidence=ev(4) + ev(10) + ev(page))
        add_unresolved(r, "The scope/general procedure says 'stipulated time' but gives no numeric duration for this individual experiment.")
    add_scope_workup(r, ev(10), "5:95")

r = reactions[rid("scope_table4_7bc")]
add_note(r,
         "The section heading and Table 4 identify product 7bc, but the individual paragraph says 'afforded 3ad'.",
         None, ev(4) + ev(11), status="conflict",
         alternatives=[F("7bc", mid("7bc"), ev(4) + ev(11)),
                       F("3ad", mid("3ad"), ev(11))])
add_unresolved(r, "Source label conflict: heading/Table 4 give 7bc, while the individual prose prints 3ad.")


# Isolated mechanistic intermediate 8bb: page 11 provides the full preparation.
r = reactions[rid("mechanistic_8bb_isolation")]
add_page_evidence(r, 11)
add_occurrence("8bb", 11)
set_mass(participant(r, "6b"), 39.0, ev(11))
set_mass(participant(r, "2b"), 69.6, ev(11))
set_mass(participant(r, "DBU"), 46.0, ev(11))
add_product_mass(r, "8bb", 63.0, ev(11))
enrich_yield_evidence(r, ev(11))
r["procedure"]["addition_order"] = F(
    "charge 6b, 2b and 1.5 mL toluene; then add DBU in 0.5 mL toluene",
    ["6b", "2b", "1.5 mL toluene", "DBU in 0.5 mL toluene"], ev(11))
add_scope_workup(r, ev(11), "10:90")

# The separate 8bb -> 7bb arrow gives no numeric time in the inspected source.
r = reactions[rid("mechanistic_8bb_to_7bb")]
r["steps"][0]["time"] = F(evidence=ev(3) + ev(4))
add_unresolved(r, "No numeric duration is reported for conversion of 8bb to 7bb; the earlier machine value was removed.")


# General Procedure viii and Table 5 individual paragraphs.
scope9 = {
    "9aa": (12, 36.2, 62.4, 58.0), "9ab": (12, 36.2, 69.6, 58.4),
    "9ac": (12, 36.2, 70.6, 65.6), "9ad": (12, 36.2, 73.2, 70.0),
    "9af": (12, 36.2, 63.8, 54.2), "9ba": (12, 39.0, 62.4, 54.5),
}
for product, (page, imine_mass, allene_mass, product_mass) in scope9.items():
    r = reactions[rid(f"scope_table5_{product}")]
    evidence = ev(4) + ev(11) + ev(12)
    add_page_evidence(r, 11)
    add_page_evidence(r, 12)
    add_occurrence(product, page)
    set_mass(r["participants"][0], imine_mass, ev(12))
    set_mass(r["participants"][1], allene_mass, ev(12))
    set_mass(participant(r, "PPh3"), 10.5, ev(12))
    add_product_mass(r, product, product_mass, ev(12))
    enrich_yield_evidence(r, ev(12))
    r["procedure"]["general_procedure_reference"] = F(
        "General Procedure viii for compounds 9aa-9ba", "General Procedure viii", ev(11))
    r["procedure"]["addition_order"] = F(
        "charge cyclic sulfonyl imine, PPh3 and 1.0 mL toluene; add allenoate in 1.0 mL toluene gradually over 30 min at 80 °C",
        ["cyclic sulfonyl imine", "PPh3", "1.0 mL toluene", "allenoate in 1.0 mL toluene over 30 min at 80 °C"],
        ev(11), method="inherited_general_condition")
    add_step(r, "allenoate_addition", "add δ-acetoxy allenoate solution gradually over 30 min",
             ev(11), "reaction", [r["participants"][1]["participant_id"], participant(r, "toluene")["participant_id"]],
             temperature=80, time=(30, "min"), prepend=True)
    # "stipulated time" has no numeric value in General Procedure viii or the
    # individual paragraphs.
    reaction_step = next(s for s in r["steps"] if s["step_id"] == "reaction")
    reaction_step["time"] = F(evidence=evidence)
    add_unresolved(r, "General Procedure viii says 'stipulated time'; no numeric hold time is reported in the inspected pages.")
    add_scope_workup(r, ev(11), "30:70")

r = reactions[rid("scope_table5_9ba")]
add_note(r,
         "Table 5 and the 9ba label map the substrate to 6b, while the individual 9ba paragraph prints 'cyclic sulfonyl imine 6a (39.0 mg, 0.20 mmol)'.",
         None, ev(4) + ev(12), status="conflict",
         alternatives=[F("6b", mid("6b"), ev(4)), F("6a", mid("6a"), ev(12))])
add_unresolved(r, "Source reactant-label conflict for 9ba: Table 5 indicates 6b, but the individual paragraph prints 6a; no silent correction was made.")


# Source-context records: characterization blocks and general procedures are
# not counted as separate synthetic experiments.
new_context = [
    {
        "context_id": f"{DOC}:context:characterization_pages7_12",
        "document_id": DOC,
        "pages": [7, 8, 9, 10, 11, 12],
        "kind": "characterization",
        "evidence_ids": [e for page in range(7, 13) for e in ev(page)],
        "note": "NMR, IR, HRMS and melting-point blocks characterize named compounds; spectral peak lists, molecular-ion values and fragments are not synthesis experiments.",
    },
    {
        "context_id": f"{DOC}:context:general_procedures_iii_iv_vi_viii",
        "document_id": DOC,
        "pages": [8, 10, 11],
        "kind": "general_procedure",
        "evidence_ids": ev(8) + ev(10) + ev(11),
        "note": "General Procedures iii, iv, vi and viii were inherited only by matching scope records; they are not additional experiments.",
    },
    {
        "context_id": f"{DOC}:context:associated_content_page12",
        "document_id": DOC,
        "pages": [12],
        "kind": "associated_content_notice",
        "evidence_ids": ev(12),
        "note": "The paper states that Supporting Information is available online, but no SI file was supplied in this batch; no SI facts or association were imported.",
    },
]
known_context = {c.get("context_id") for c in bundle["source_context"]}
bundle["source_context"].extend(c for c in new_context if c["context_id"] not in known_context)


# Raw visual transcription kept separately from the formal bundle.
raw = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "visual_pages": list(range(7, 13)),
    "individual_product_masses_mg": {
        **{k: v[-1] for k, v in scope3_masses.items()},
        **{k: v[-1] for k, v in scope4.items()},
        **{k: v[-1] for k, v in scope5.items()},
        **{k: v[-1] for k, v in scope7.items()},
        "8bb": 63.0,
        **{k: v[-1] for k, v in scope9.items()},
    },
    "source_conflicts": [
        {"record": rid("scope_table2_3bb"), "field": "2b amount", "values": ["0.24 mmol", "0.20 mmol"], "pages": [3, 6, 7]},
        {"records": [rid(f"scope_table4_{p}") for p in scope7], "field": "DBU amount", "values": ["0.20 mmol", "0.10 mmol"], "pages": [10, 11]},
        {"record": rid("scope_table4_7bc"), "field": "product label", "values": ["7bc", "3ad"], "pages": [4, 11]},
        {"record": rid("scope_table5_9ba"), "field": "reactant label", "values": ["6b", "6a"], "pages": [4, 12]},
    ],
    "limitations": [
        "Visual transcription is not a source-reviewed chemistry record.",
        "Characterization peak lists were not converted into synthesis experiments.",
        "The Supporting Information was mentioned by the article but not supplied in this batch.",
    ],
}
(OUT / "raw_model/visual_transcription_pages_0007_0012.json").write_text(
    json.dumps(raw, ensure_ascii=False, indent=2) + "\n")


# Coverage maps every inspected page to the already existing experiment IDs;
# general procedures and characterization have explicit non-experiment regions.
R = lambda suffix: rid(suffix)
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 7,
            "state": "attempted",
            "regions": [
                {"region": "top continuation and compound 3af", "finding": "3ae characterization continuation followed by the individual 3af synthesis; characterization is source_context, not a new reaction", "evidence_ids": ev(7)},
                {"region": "middle and lower compound paragraphs", "finding": "individual syntheses 3ba-3cc with isolated masses/yields and analytical characterization", "evidence_ids": ev(7)},
                {"region": "bottom continuation", "finding": "3da synthesis begins and continues on page 8", "evidence_ids": ev(7)},
            ],
            "record_ids": [R(f"scope_table2_{p}") for p in ["3af", "3ba", "3bb", "3bc", "3bd", "3ca", "3cb", "3cc", "3da"]],
            "unresolved": ["3bb reports 2b as 0.20 mmol, conflicting with the 0.24 mmol Table 2/general-procedure value.", "Numeric reaction durations are not printed in these individual entries."],
        },
        {
            "page": 8,
            "state": "attempted",
            "regions": [
                {"region": "top compound paragraphs", "finding": "3da continuation plus individual syntheses 3db, 3dc and 3ea", "evidence_ids": ev(8)},
                {"region": "General Procedure iii and compounds 4aa-4af", "finding": "30 min gradual allenoate addition at rt, 12 h hold inherited from Table 3, workup/purification, and four individual product examples", "evidence_ids": ev(8)},
                {"region": "General Procedure iv and 5aa start", "finding": "30 min gradual allenoate addition at 80 °C, workup/purification; 5aa paragraph continues on page 9", "evidence_ids": ev(8)},
            ],
            "record_ids": [R(f"scope_table2_{p}") for p in ["3da", "3db", "3dc", "3ea"]] + [R(f"scope_table3_{p}") for p in ["4aa", "4ba", "4ca", "4af", "5aa"]],
            "unresolved": ["Full structures and stereochemical representations remain unverified source-drawing candidates."],
        },
        {
            "page": 9,
            "state": "attempted",
            "regions": [
                {"region": "5aa continuation through 5af", "finding": "individual product paragraphs with isolated masses/yields and characterization", "evidence_ids": ev(9)},
                {"region": "5ba through 5bd", "finding": "individual product paragraphs; 5bd characterization continues on page 10", "evidence_ids": ev(9)},
            ],
            "record_ids": [R(f"scope_table3_{p}") for p in ["5aa", "5ab", "5ac", "5ad", "5ae", "5af", "5ba", "5bb", "5bc", "5bd"]],
            "unresolved": ["Spectral assignments were retained as characterization context, not independently structure-verified."],
        },
        {
            "page": 10,
            "state": "attempted",
            "regions": [
                {"region": "5bd continuation and compounds 5ca-5da", "finding": "individual scope syntheses with isolated masses/yields and characterization", "evidence_ids": ev(10)},
                {"region": "Scale-Up Experiments", "finding": "complete 1.0 mmol preparations of 3ab and 5ab, including 45 min heated addition for 5ab, workup, chromatography and isolated mass/yield", "evidence_ids": ev(10)},
                {"region": "General Procedure vi and 7aa", "finding": "100 °C one-pot procedure and individual 7aa example; DBU amount conflicts between general and individual text", "evidence_ids": ev(10)},
            ],
            "record_ids": [R(f"scope_table3_{p}") for p in ["5bd", "5ca", "5cb", "5cc", "5cd", "5da"]] + [R("scaleup_3ab"), R("scaleup_5ab"), R("scope_table4_7aa")],
            "unresolved": ["General Procedure vi gives DBU 0.20 mmol, while individual examples print 30.4 mg (0.10 mmol)."],
        },
        {
            "page": 11,
            "state": "attempted",
            "regions": [
                {"region": "compounds 7ab-7af", "finding": "individual Table 4 syntheses and characterization; 7bc prose product label conflicts with its heading/Table 4", "evidence_ids": ev(11)},
                {"region": "compound 8bb preparation", "finding": "mechanistic-study isolation of 8bb at 60 °C for 12 h, with full workup, 63.0 mg and 74% isolated yield", "evidence_ids": ev(11)},
                {"region": "General Procedure viii and 9aa start", "finding": "30 min allenoate addition at 80 °C and workup/purification; hold time is only called stipulated; 9aa continues on page 12", "evidence_ids": ev(11)},
            ],
            "record_ids": [R(f"scope_table4_{p}") for p in ["7ab", "7ac", "7ba", "7bb", "7bc", "7ae", "7af"]] + [R("mechanistic_8bb_isolation"), R("scope_table5_9aa")],
            "unresolved": ["7bc paragraph says 'afforded 3ad' although its heading and Table 4 identify 7bc.", "General Procedure viii does not give a numeric hold time."],
        },
        {
            "page": 12,
            "state": "attempted",
            "regions": [
                {"region": "compounds 9aa-9ba", "finding": "individual Table 5 syntheses with isolated masses/yields and characterization", "evidence_ids": ev(12)},
                {"region": "9ba paragraph", "finding": "reactant is printed as 6a, conflicting with Table 5/9ba mapping to 6b", "evidence_ids": ev(12)},
                {"region": "Associated Content and article metadata", "finding": "SI availability notice, authors and acknowledgments; no supplied SI file and no additional synthesis record", "evidence_ids": ev(12)},
            ],
            "record_ids": [R(f"scope_table5_{p}") for p in ["9aa", "9ab", "9ac", "9ad", "9af", "9ba"]],
            "unresolved": ["9ba reactant label is internally inconsistent (6b in Table 5 versus 6a in the individual paragraph).", "SI is mentioned but not supplied; association remains unestablished.", "Numeric hold times for the Table 5 examples are unreported in the inspected pages."],
        },
    ],
}
(OUT / "coverage/pages_0007_0012.json").parent.mkdir(parents=True, exist_ok=True)
(OUT / "coverage/pages_0007_0012.json").write_text(json.dumps(coverage, ensure_ascii=False, indent=2) + "\n")


# The compatibility schema also requires evidence lists to be set-like.  Some
# cross-page records intentionally cite the same page as both general and
# individual evidence; collapse only exact duplicate IDs, preserving order.
def dedupe_evidence(node):
    if isinstance(node, dict):
        if isinstance(node.get("evidence_ids"), list):
            node["evidence_ids"] = list(dict.fromkeys(node["evidence_ids"]))
        for value in node.values():
            dedupe_evidence(value)
    elif isinstance(node, list):
        for value in node:
            dedupe_evidence(value)


dedupe_evidence(bundle)
bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "updated_reactions": 48,
    "total_reactions": len(bundle["reactions"]),
    "total_molecules": len(bundle["molecules"]),
    "pages": [7, 8, 9, 10, 11, 12],
}, ensure_ascii=False))
