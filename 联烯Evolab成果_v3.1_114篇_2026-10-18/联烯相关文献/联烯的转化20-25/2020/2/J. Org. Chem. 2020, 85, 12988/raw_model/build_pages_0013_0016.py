#!/usr/bin/env python3
"""Build the page-13--16 v3.1 extraction increment after visual inspection."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import standard_records as st

RUN = ROOT / "runs/current"
DOC = "doc_2a371ced54a752e8"
BATCH = "allene-stereo-v31-200-20261017"
TITLE = (
    "Bi(cyclopentyl)diol-Derived Boronates in Highly Enantioselective "
    "Chiral Phosphoric Acid-Catalyzed Allylation, Propargylation, and Crotylation of Aldehydes"
)


def ev(page: int, kind: str = "image") -> str:
    return f"{DOC}:page{page}:{kind}"


def evid(*pages: int) -> list[str]:
    out: list[str] = []
    for page in pages:
        out.extend([ev(page), ev(page, "text")])
    return out


def src(raw, normalized=None, evidence=()):
    return st.fact(raw, raw if normalized is None else normalized, evidence)


def append_unique(seq: list, values) -> None:
    for value in values:
        if value not in seq:
            seq.append(value)


def add_fact_evidence(fact: dict, evidence) -> None:
    append_unique(fact.setdefault("evidence_ids", []), evidence)


def add_note(record: dict, raw: str, normalized: str, evidence) -> None:
    if not any(n.get("raw_value") == raw for n in record["notes"]):
        record["notes"].append(src(raw, normalized, evidence))


def add_unresolved(record: dict, item: str) -> None:
    if item not in record["review"]["unresolved"]:
        record["review"]["unresolved"].append(item)


def step(step_id: str, evidence, participant_ids, operation: str, *, stage="reaction", temperature=None, time=None):
    row = st.new_step(step_id, evidence, stage=stage)
    row["participant_ids"] = list(participant_ids)
    row["operation"] = src(operation, operation, evidence)
    if temperature is not None:
        row["temperature"] = src(temperature[0], temperature[1], evidence)
    if time is not None:
        row["time"] = src(time[0], time[1], evidence)
    return row


standard = st.read(RUN / "standard_records.json")
dataset = st.read(RUN / "dataset.json")
existing_reactions = {r["reaction_id"]: r for r in standard["reactions"]}
existing_molecules = {m["molecule_id"]: m for m in standard["molecules"]}
bundle = st.empty(dataset["documents"])
reactions: list[dict] = []
molecules: dict[str, dict] = {}
updated_ids: list[str] = []
stereo_diagnostics: dict[str, dict] = {}


def copied_reaction(short_id: str, new_evidence) -> dict:
    full = f"{DOC}:{short_id}"
    row = copy.deepcopy(existing_reactions[full])
    append_unique(row["evidence_ids"], new_evidence)
    add_fact_evidence(row["location"], new_evidence)
    pages = set(row["location"]["normalized_value"].get("pages", []))
    pages.update(int(e.split(":page")[1].split(":")[0]) for e in new_evidence)
    row["location"]["normalized_value"]["pages"] = sorted(pages)
    row["review"].update(status="needs_review", method="model", reviewer=None)
    updated_ids.append(full)
    return row


def touch_molecule(mid: str, label: str, evidence) -> None:
    molecule = molecules.setdefault(mid, copy.deepcopy(existing_molecules[mid]))
    matching = [o for o in molecule["occurrences"] if o["compound_label"].get("normalized_value") == label]
    if matching:
        append_unique(matching[0]["evidence_ids"], evidence)
        add_fact_evidence(matching[0]["compound_label"], evidence)
    else:
        molecule["occurrences"].append({
            "document_id": DOC,
            "compound_label": src(label, label, evidence),
            "evidence_ids": list(evidence),
        })
    for fact in (molecule["names"]["original"], molecule["names"]["english"]):
        add_fact_evidence(fact, evidence)
    stereo = molecule["stereochemistry"]
    for key in ("chirality_status", "geometry_status", "representation_status", "assignment", "isomeric_smiles"):
        add_fact_evidence(stereo[key], evidence)
    smiles = molecule.get("smiles", {}).get("normalized_value")
    if smiles:
        stereo_diagnostics[mid] = st.stereo_inspect(smiles)


def outcome(record: dict, oid: str) -> dict:
    return next(o for o in record["outcomes"] if o["outcome_id"] == oid)


# Page 13: product 35 characterization supplements its existing Scheme 3 experiment.
e13 = evid(13)
r = copied_reaction("scheme3_product_35_method_b", e13)
touch_molecule("m_product_35", "35", e13)
for oid in ("yield", "ee"):
    o = outcome(r, oid)
    add_fact_evidence(o["value"], e13)
    add_fact_evidence(o["measurement_basis"], e13)
outcome(r, "ee")["measurement_basis"] = src(
    "Chiralpak AD-H; hexane/iPrOH 90/10; 1.0 mL/min; tmajor 6.273 min, tminor 6.518 min",
    "Chiralpak AD-H; hexane/iPrOH 90/10; 1.0 mL/min; tmajor 6.273 min, tminor 6.518 min",
    e13,
)
reactions.append(r)


# Page 13: full Section 4.6 procedure and individual HPLC bases for products 36--45.
hplc = {
    36: "Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 10.40 min, tminor 12.37 min",
    37: "Chiralcel OD-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 20.57 min, tminor 23.23 min",
    38: "Chiralcel OJ-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 29.91 min, tmajor 33.41 min",
    39: "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tmajor 8.72 min, tminor 10.75 min",
    40: "Chiralcel OD-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 36.16 min, tminor 41.49 min",
    41: "Chiralpak AD-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 16.07 min, tminor 20.40 min",
    42: "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tmajor 11.88 min, tminor 16.97 min",
    43: "Chiralcel OD-H; hexane/iPrOH 95/5; 1.0 mL/min; tmajor 15.89 min, tminor 22.40 min",
    44: "Chiralcel OJ-H; hexane/iPrOH 90/10; 1.0 mL/min; tminor 16.64 min, tmajor 18.13 min",
    45: "Chiralcel OD-H; hexane/iPrOH 90/10; 1.0 mL/min; tmajor 9.51 min, tminor 13.17 min",
}
general_propargylation = (
    "load a screw-cap reaction tube with a stir bar and 4 Å molecular sieves (50 mg); evacuate, "
    "flame-dry and back-fill with argon; add (R)-TRIP-PA PA1 (5 mol%), freshly distilled aldehyde "
    "(0.1 mmol) and dry cyclohexane (1.0 mL); add B10 (0.12 mmol) by syringe dropwise over 30 s; "
    "stir 8 h at room temperature; load directly onto silica and flash chromatograph with ethyl "
    "acetate/hexane 1:9"
)
for number, basis in hplc.items():
    r = copied_reaction(f"scheme4_product_{number}", e13)
    touch_molecule(f"m_product_{number}", str(number), e13)
    for oid in ("yield", "ee"):
        add_fact_evidence(outcome(r, oid)["value"], e13)
    outcome(r, "ee")["measurement_basis"] = src(basis, basis, e13)
    r["procedure"]["original_text"] = src(general_propargylation, general_propargylation, e13)
    r["procedure"]["general_procedure_reference"] = src(
        "Section 4.6, General Procedure for the Propargylation of Aldehydes",
        "Section 4.6, General Procedure for the Propargylation of Aldehydes",
        e13,
    )
    r["procedure"]["addition_order"] = src(
        "4 Å MS; PA1, aldehyde and cyclohexane; then B10 dropwise over 30 s",
        ["4 Å MS", "PA1", "aldehyde", "cyclohexane", "B10 dropwise over 30 s"],
        e13,
    )
    active = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    r["steps"] = [
        step("prepare_under_argon", e13, active, "evacuate, flame-dry and back-fill the screw-cap tube with argon"),
        step("charge_catalyst_aldehyde_solvent", e13, active, "add PA1, freshly distilled aldehyde and dry cyclohexane"),
        step("add_b10", e13, active, "add B10 by syringe dropwise over 30 s", time=("30 s", {"number": 30, "unit": "s"})),
        step("room_temperature_hold", e13, active, "stir the reaction mixture", temperature=("room temperature", "room temperature"), time=("8 h", {"number": 8, "unit": "h"})),
        step("silica_purification", e13, active, "load directly onto silica and flash chromatograph with ethyl acetate/hexane 1:9", stage="purification"),
    ]
    r["steps"][0]["atmosphere"] = src("argon", "argon", e13)
    r["review"]["unresolved"] = [
        item for item in r["review"]["unresolved"]
        if "procedure details" not in item.lower()
    ]
    add_unresolved(r, "The source does not state a separate temperature for the 30 s B10 addition; it remains unreported for that step.")
    reactions.append(r)


# Pages 13--14: application compounds supplement the three existing Scheme 5 experiments.
application_specs = {
    "page6_19_to_47": {
        "product": "m_product_47", "label": "47", "yield_basis": "isolated after flash chromatography (hexane/ethyl acetate 8:1)",
        "ee_basis": "Chiralpak IF; hexane/iPrOH 80/20; 1.0 mL/min; tminor 7.233 min, tmajor 8.502 min",
        "reference": "literature procedure 15b", "pages": evid(13, 14),
    },
    "page6_19_to_49": {
        "product": "m_product_49", "label": "49", "yield_basis": "isolated after flash chromatography (hexane/ethyl acetate 10:1)",
        "ee_basis": "Chiralcel OJ-H; hexane/iPrOH 98/2; 1.0 mL/min; tmajor 5.293 min, tminor 8.161 min",
        "reference": "literature procedure 27", "pages": evid(14),
    },
    "page6_49_to_50": {
        "product": "m_product_50", "label": "50", "yield_basis": "isolated after flash chromatography (hexane/ethyl acetate 2:1)",
        "ee_basis": "Chiralpak AD-H; hexane/iPrOH 75/25; 1.0 mL/min; tmajor 6.826 min, tminor 7.334 min",
        "reference": "literature procedure 27", "pages": evid(14),
    },
}
for short_id, spec in application_specs.items():
    r = copied_reaction(short_id, spec["pages"])
    touch_molecule(spec["product"], spec["label"], spec["pages"])
    outcome(r, "yield")["measurement_basis"] = src(spec["yield_basis"], spec["yield_basis"], evid(14))
    outcome(r, "ee")["measurement_basis"] = src(spec["ee_basis"], spec["ee_basis"], evid(14))
    add_fact_evidence(outcome(r, "yield")["value"], evid(14))
    add_fact_evidence(outcome(r, "ee")["value"], evid(14))
    r["procedure"]["general_procedure_reference"] = src(spec["reference"], spec["reference"], spec["pages"])
    r["review"]["unresolved"] = [
        item for item in r["review"]["unresolved"]
        if "all outcome measurement bases" not in item.lower()
    ]
    add_unresolved(r, "Reagent charges and reaction duration are not reported in the inspected application characterization section.")
    if short_id in ("page6_19_to_49", "page6_49_to_50"):
        statement = "The page-14 preparation preface states that compounds 49 and 50 used alcohol 19 with 89% ee and dr >20:1."
        add_note(r, statement, statement, evid(14))
        add_unresolved(r, "Starting-alcohol statement conflicts with the 99% ee characterization of alcohol 19 elsewhere; the application-specific 89% ee statement is retained without reconciliation.")
    reactions.append(r)


# Page 14, Section 4.8: add the reported 2 h boronate addition, isolated masses and purification.
gram_specs = {
    "page6_gram_allylation": ("m_product_2", "m_B2", "1.43 g"),
    "page6_gram_crotylation": ("m_product_22", "m_B3", "2.13 g"),
    "page6_gram_propargylation": ("m_product_36", "m_B10", "1.41 g"),
}
for short_id, (product_mid, boronate_mid, mass_raw) in gram_specs.items():
    r = copied_reaction(short_id, evid(14))
    active = [p["participant_id"] for p in r["participants"] if p["role"] != "product"]
    mass_value = {"number": float(mass_raw.split()[0]), "unit": "g"}
    amount = st.new_outcome("isolated_product_amount", "product_amount", [product_mid], evid(14))
    amount["value"] = src(mass_raw, mass_value, evid(14))
    r["outcomes"] = [o for o in r["outcomes"] if o["outcome_id"] != "isolated_product_amount"] + [amount]
    outcome(r, "yield")["measurement_basis"] = src(
        "isolated after silica-gel flash chromatography (ethyl acetate/hexane 1:9)",
        "isolated after silica-gel flash chromatography (ethyl acetate/hexane 1:9)",
        evid(14),
    )
    add_fact_evidence(outcome(r, "yield")["value"], evid(14))
    add_fact_evidence(outcome(r, "ee")["value"], evid(14))
    existing_hold = copy.deepcopy(next(s for s in r["steps"] if s["step_id"] == "scale_up_hold"))
    r["steps"] = [
        step("boronate_addition", evid(14), active, f"add {boronate_mid.removeprefix('m_')} to the reaction system over 2 h", time=("2 h", {"number": 2, "unit": "h"})),
        existing_hold,
        step("silica_purification", evid(14), active, "purify by silica-gel flash chromatography with ethyl acetate/hexane 1:9", stage="purification"),
    ]
    r["procedure"]["original_text"] = src(
        "conduct on a 10.0 mmol scale of benzaldehyde according to the above general allylboration method; add B2, B3 or B10 over 2 h; purify by silica-gel flash chromatography with ethyl acetate/hexane 1:9",
        "conduct on a 10.0 mmol scale of benzaldehyde according to the above general allylboration method; add the row-specific boronate over 2 h; purify by silica-gel flash chromatography with ethyl acetate/hexane 1:9",
        evid(14),
    )
    r["procedure"]["addition_order"] = src(
        f"{boronate_mid.removeprefix('m_')} added to the reaction system over 2 h",
        [f"{boronate_mid.removeprefix('m_')} added over 2 h"],
        evid(14),
    )
    r["review"]["unresolved"] = [
        item for item in r["review"]["unresolved"]
        if "outcome measurement bases" not in item.lower()
    ]
    add_unresolved(r, "The gram-scale section does not report the ee/dr measurement method; those measurement bases remain unreported.")
    add_unresolved(r, "The gram-scale section does not explicitly pair a temperature with the 2 h boronate-addition stage.")
    statement = "Section 4.8 collectively describes a 10.0 mmol benzaldehyde scale and use of the above general allylboration method."
    add_note(r, statement, statement, evid(14))
    if short_id == "page6_gram_crotylation":
        add_unresolved(r, "Source conflict retained: Section 4.8 says benzaldehyde, while the Scheme 5c crotylation row depicts methyl 4-formylbenzoate/product 52.")
    reactions.append(r)


# Every molecule referenced by the increment must be present in the incoming bundle.
for r in reactions:
    for participant in r["participants"]:
        mid = participant["molecule_id"]
        molecules.setdefault(mid, copy.deepcopy(existing_molecules[mid]))
    for o in r["outcomes"]:
        for mid in o["target_molecule_ids"]:
            molecules.setdefault(mid, copy.deepcopy(existing_molecules[mid]))

bundle["molecules"] = list(molecules.values())
bundle["reactions"] = reactions
bundle["source_context"] = [
    {
        "context_id": f"{DOC}:page13_spectral_characterization",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": evid(13),
        "finding": "The opening continuation for compound 34 and the NMR/HRMS/optical-rotation blocks for 35--45 are characterization context, not additional synthetic experiments.",
    },
    {
        "context_id": f"{DOC}:page14_spectral_characterization",
        "document_id": DOC,
        "kind": "characterization",
        "evidence_ids": evid(14),
        "finding": "Spectral, melting-point and optical-rotation data for 47, 49 and 50 are characterization context; yield, purification and chiral-HPLC facts supplement existing application records.",
    },
    {
        "context_id": f"{DOC}:page14_si_availability",
        "document_id": DOC,
        "kind": "associated_content_notice",
        "evidence_ids": evid(14),
        "finding": "The article states that Supporting Information is available online, but no SI file was supplied or associated in this workspace.",
    },
    {
        "context_id": f"{DOC}:pages15_16_references",
        "document_id": DOC,
        "kind": "references",
        "evidence_ids": evid(15, 16),
        "finding": "Pages 15--16 contain bibliographic references only; cited experiments were not promoted to current-paper experiments.",
    },
]

st.write(RUN / "submitted.json", bundle)
st.write(RUN / "raw_model/pages_0013_0016_standard_stereo_checks.json", stereo_diagnostics)

raw = {
    "batch_id": BATCH,
    "document_id": DOC,
    "title": TITLE,
    "visual_pages_inspected": [13, 14, 15, 16],
    "new_record_ids": [],
    "updated_record_ids": updated_ids,
    "new_record_count": 0,
    "updated_record_count": len(updated_ids),
    "decisions": [
        "Characterization paragraphs for 35--45, 47, 49 and 50 supplement existing experiment records and are not duplicate experiments.",
        "Section 4.6 was structured into vessel preparation, charge, 30 s B10 addition, 8 h room-temperature hold and purification stages.",
        "Section 4.8 retains three distinct gram-scale rows and adds the row-specific 2 h boronate addition and isolated product mass.",
        "The application-specific 89% ee statement for starting alcohol 19 and the gram-crotylation aldehyde conflict are retained unresolved.",
        "Standalone NMR/HRMS/optical-rotation data and pages 15--16 references are source context, not synthesis records.",
        "The article's SI availability notice is not treated as a supplied or linked SI document.",
    ],
    "records": [
        {
            "reaction_id": r["reaction_id"],
            "record_type": r["record_type"],
            "location": r["location"]["raw_value"],
            "evidence_ids": r["evidence_ids"],
        }
        for r in reactions
    ],
    "molecule_stereo_review": [
        {
            "molecule_id": mid,
            "chirality_status": molecules[mid]["stereochemistry"]["chirality_status"]["normalized_value"],
            "geometry_status": molecules[mid]["stereochemistry"]["geometry_status"]["normalized_value"],
            "representation_status": molecules[mid]["stereochemistry"]["representation_status"]["normalized_value"],
            "isomeric_smiles": molecules[mid]["stereochemistry"]["isomeric_smiles"]["normalized_value"],
        }
        for mid in sorted(stereo_diagnostics)
    ],
}
st.write(RUN / "raw_model/pages_0013_0016_raw_extraction.json", raw)

page13_records = [f"{DOC}:scheme3_product_34_method_b", f"{DOC}:scheme3_product_35_method_b"]
page13_records.extend(f"{DOC}:scheme4_product_{n}" for n in range(36, 46))
page13_records.append(f"{DOC}:page6_19_to_47")
page14_records = [
    f"{DOC}:page6_19_to_47", f"{DOC}:page6_19_to_49", f"{DOC}:page6_49_to_50",
    f"{DOC}:page6_gram_allylation", f"{DOC}:page6_gram_crotylation", f"{DOC}:page6_gram_propargylation",
]
coverage = {
    "document_id": DOC,
    "semantic_accuracy_verified": False,
    "pages": [
        {
            "page": 13,
            "state": "attempted",
            "regions": [
                {"region": "top continuation of Sections 4.5.18--4.5.19", "finding": "Compound 34 spectral continuation and compound 35 characterization; only 35 supplies new yield/HPLC evidence for an existing reaction record.", "evidence_ids": evid(13)},
                {"region": "Section 4.6 general procedure", "finding": "Propargylation setup, quantities, 30 s B10 addition, 8 h room-temperature hold and silica purification were structured as ordered stages.", "evidence_ids": evid(13)},
                {"region": "Sections 4.6.1--4.6.10", "finding": "Products 36--45 yields and compound-specific chiral-HPLC methods supplement existing Scheme 4 experiments; spectral peak lists remain characterization context.", "evidence_ids": evid(13)},
                {"region": "bottom start of Section 4.7", "finding": "Compound 47 is stated to follow literature procedure 15b; no unprinted literature conditions were imported.", "evidence_ids": evid(13)},
            ],
            "record_ids": page13_records,
            "unresolved": ["Separate temperatures for the PA1/aldehyde charge and 30 s B10-addition stages are not stated.", "The opening compound-34 spectral continuation was retained as context and not mistaken for a new synthesis."],
        },
        {
            "page": 14,
            "state": "attempted",
            "regions": [
                {"region": "Sections 4.7.1--4.7.3", "finding": "Application products 47, 49 and 50 provide isolated yields, purification ratios and chiral-HPLC bases for existing experiments.", "evidence_ids": evid(14)},
                {"region": "Section 4.8 gram-scale procedure", "finding": "Three existing gram-scale experiments were supplemented with 10 mmol context, row-specific 2 h boronate addition, isolated masses and silica chromatography details.", "evidence_ids": evid(14)},
                {"region": "Associated Content", "finding": "SI availability is mentioned, but no SI document was supplied; no SI association was inferred.", "evidence_ids": evid(14)},
                {"region": "spectral and physical-property blocks", "finding": "NMR/HRMS/optical-rotation/melting-point data were retained as characterization context, not independent reactions.", "evidence_ids": evid(14)},
            ],
            "record_ids": page14_records,
            "unresolved": ["Compounds 49/50 are introduced as using alcohol 19 with 89% ee, conflicting with the 99% ee characterization elsewhere.", "Section 4.8 says benzaldehyde collectively, whereas the Scheme 5c crotylation row depicts methyl 4-formylbenzoate/product 52.", "The gram-scale section does not give ee/dr measurement methods or a temperature specifically paired with the 2 h boronate addition."],
        },
        {
            "page": 15,
            "state": "attempted",
            "regions": [
                {"region": "full page, References 3--10", "finding": "Bibliographic references and cited prior work only; no present-paper experiment, reaction arrow, table row or footnote requiring a record.", "evidence_ids": evid(15)},
            ],
            "record_ids": [],
            "unresolved": [],
        },
        {
            "page": 16,
            "state": "attempted",
            "regions": [
                {"region": "full page, References 10--27", "finding": "Bibliographic references and cited prior work only; no present-paper experiment, reaction arrow, table row or footnote requiring a record.", "evidence_ids": evid(16)},
            ],
            "record_ids": [],
            "unresolved": [],
        },
    ],
}
st.write(RUN / "coverage/pages_0013_0016.json", coverage)
print(json.dumps({"new_records": 0, "updated_records": len(updated_ids), "molecules_in_increment": len(molecules), "coverage_pages": len(coverage["pages"])}, ensure_ascii=False))
