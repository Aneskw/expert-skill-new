#!/usr/bin/env python3
"""Incremental source extraction for doc_518d3edb75482abc, PDF pages 7-8.

The script does not alter tools, schemas, or quality rules.  It first registers
the supplied page-7/page-8 images as prepared-document evidence, then builds an
incremental v3.1 submission from a fresh ``standard init`` template.  Existing
reaction IDs are reused only to append cross-page evidence and source-derived
product details; no prior record is counted as a new page-7/page-8 experiment.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import standard_records as st


DOC = "doc_518d3edb75482abc"
P = {i: f"{DOC}:page{i}:text" for i in range(1, 9)}
IMG = {i: f"{DOC}:page{i}:image" for i in range(1, 9)}
SELECTED_RIDS = [f"{DOC}:scope{i}" for i in range(28, 35)]


def load(path: Path):
    return json.loads(path.read_text())


def dump(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_ids(sequence, *values):
    for value in values:
        if value not in sequence:
            sequence.append(value)


def augment_evidence():
    dataset_path = RUN / "dataset.json"
    dataset = load(dataset_path)
    document = next(d for d in dataset["documents"] if d["document_id"] == DOC)
    source = RUN / document["source_pdf"]
    if sha256(source) != document["sha256"]:
        raise ValueError("Prepared source PDF hash no longer matches the document registry")

    known = {e["evidence_id"] for e in document["evidence"]}
    for page in (7, 8):
        asset = RUN / f"assets/{DOC}/page{page}.jpg"
        if not asset.is_file() or asset.stat().st_size == 0:
            raise ValueError(f"Missing supplied page image: {asset}")
        evidence_id = IMG[page]
        if evidence_id not in known:
            document["evidence"].append({
                "evidence_id": evidence_id,
                "document_id": DOC,
                "page": page,
                "kind": "page_image",
                "text": None,
                "asset": f"assets/{DOC}/page{page}.jpg",
                "bbox": None,
                "bbox_units": None,
                "provenance": "rendered_source_pdf",
            })
            known.add(evidence_id)

    tasks_path = RUN / "tasks.json"
    tasks = load(tasks_path)
    for task in tasks:
        if task.get("evidence_id") in (P[7], P[8]):
            page = 7 if task["evidence_id"] == P[7] else 8
            task["visual_evidence_id"] = IMG[page]
            task["state"] = "attempted"
            task["reason"] = (
                "Original PDF page image inspected at original detail; findings, existing-record links, "
                "and unresolved items are recorded in coverage/pages_0007_0008.json."
            )

    dump(dataset_path, dataset)
    dump(tasks_path, tasks)
    print(json.dumps({"augmented_image_evidence": [IMG[7], IMG[8]]}, ensure_ascii=False))


def inherited(raw, normalized, evidence, derivation):
    return st.fact(
        raw,
        normalized,
        evidence,
        method="inherited_general_condition",
        derivation=derivation,
    )


def direct(raw, normalized, evidence):
    return st.fact(raw, normalized, evidence)


def missing(evidence=()):
    return st.fact(evidence=evidence)


def add_fact_evidence(fact, *evidence_ids):
    add_ids(fact["evidence_ids"], *evidence_ids)


def ensure_participant(reaction, participant_id, molecule_id, role, stage, evidence):
    existing = next((p for p in reaction["participants"] if p["participant_id"] == participant_id), None)
    if existing:
        return existing
    participant = st.new_participant(participant_id, molecule_id, role, evidence, stage=stage)
    reaction["participants"].append(participant)
    return participant


def set_charge(participant, key, raw, normalized, evidence, method="inherited_general_condition"):
    participant["charges"][key] = st.fact(
        raw,
        normalized,
        evidence,
        method=method,
        derivation=(
            "The product paragraph explicitly invokes General Procedure A; this charge is inherited from that procedure."
            if method == "inherited_general_condition" else None
        ),
    )


def make_step(step_id, stage, participant_ids, operation, evidence, temperature=None,
              time=None, atmosphere=None, method="inherited_general_condition"):
    step = st.new_step(step_id, evidence, stage=stage)
    step["participant_ids"] = participant_ids
    derivation = (
        "The product paragraph explicitly invokes General Procedure A; this operation is inherited from that procedure."
        if method == "inherited_general_condition" else None
    )
    step["operation"] = st.fact(operation, operation, evidence, method=method, derivation=derivation)
    if temperature is not None:
        step["temperature"] = st.fact(
            temperature[0], temperature[1], evidence, method=method, derivation=derivation
        )
    if time is not None:
        step["time"] = st.fact(time[0], time[1], evidence, method=method, derivation=derivation)
    if atmosphere is not None:
        step["atmosphere"] = st.fact(atmosphere, atmosphere, evidence, method=method, derivation=derivation)
    return step


def product_amount(mid, mg, evidence):
    outcome = st.new_outcome(f"amount_{mid}", "product_amount", [mid], evidence)
    outcome["product_role"] = direct("isolated product", "isolated product", evidence)
    outcome["value"] = direct(f"{mg} mg", {"number": mg, "unit": "mg"}, evidence)
    outcome["measurement_basis"] = direct(
        "after flash chromatography", "after flash chromatography", evidence
    )
    outcome["scope"] = direct("product preparation paragraph", "product preparation paragraph", evidence)
    return outcome


def add_full_general_procedure(reaction, product_page_evidence, purification_text,
                               purification_molecules):
    general_ev = [P[5], IMG[5]]
    apparatus_ev = [P[4], IMG[4], P[5], IMG[5]]
    substrate = next(p for p in reaction["participants"] if p["participant_id"] == "p_substrate")
    set_charge(substrate, "amount", "0.215 mmol", {"number": 0.215, "unit": "mmol"}, general_ev)
    reaction["reaction_scale"] = inherited(
        "N-Tosyl alpha-AAP (0.215 mmol)",
        {"reference_molecule_id": substrate["molecule_id"], "amount": {"number": 0.215, "unit": "mmol"}},
        general_ev,
        "The product paragraph explicitly invokes General Procedure A, whose reference substrate charge is 0.215 mmol.",
    )
    malonate = next(p for p in reaction["participants"] if p["participant_id"] == "p_malonate")
    set_charge(malonate, "amount", "0.966 mmol", {"number": 0.966, "unit": "mmol"}, general_ev)
    reaction["solvent_conditions"]["total_volume"] = inherited(
        "25 mL of electrolyte solution",
        {"number": 25, "unit": "mL"},
        general_ev,
        "The product paragraph explicitly invokes General Procedure A, which reports 25 mL of electrolyte solution.",
    )

    ensure_participant(reaction, "p_argon", "a_argon", "reagent", "reaction", general_ev)
    ensure_participant(reaction, "p_workup_etoac", "a_etoac", "workup_agent", "workup", general_ev)
    ensure_participant(reaction, "p_celite", "a_celite", "workup_agent", "workup", general_ev)
    ensure_participant(reaction, "p_silica", "a_silica", "purification_agent", "purification", [P[4], IMG[4]])
    purification_ids = ["p_silica"]
    for molecule_id in purification_molecules:
        participant_id = "p_purification_" + molecule_id.removeprefix("a_")
        ensure_participant(
            reaction, participant_id, molecule_id, "purification_agent", "purification", product_page_evidence
        )
        purification_ids.append(participant_id)

    reaction["steps"] = [
        make_step(
            "s1_setup", "reaction", ["p_substrate", "p_malonate", "p_electrolyte"],
            "Add the N-Tosyl alpha-AAP and diethyl malonate to the electrolyte solution in a divided cell.",
            general_ev,
        ),
        make_step(
            "s2_degas", "reaction", ["p_argon"],
            "Degas the solution under argon before electrolysis.", general_ev,
            time=("15 min", {"number": 15, "unit": "min"}), atmosphere="argon",
        ),
        make_step(
            "s3_electrolysis", "reaction", ["p_substrate", "p_malonate", "p_electrolyte", "p_argon"],
            "Apply constant potential while maintaining argon; stop when current is less than 5% of its starting value.",
            general_ev,
            temperature=("room temperature", "room temperature"),
            time=("usually 5 h", {"number": 5, "unit": "h"}),
            atmosphere="argon",
        ),
        make_step("s4_remove_solvent", "workup", [], "Remove solvent under vacuum.", general_ev),
        make_step(
            "s5_wash_filter", "workup", ["p_workup_etoac", "p_celite"],
            "Wash the residue with ethyl acetate and filter the suspension through a pad of Celite.", general_ev,
        ),
        make_step("s6_evaporate", "workup", [], "Evaporate solvent under vacuum.", general_ev),
        make_step(
            "s7_product_flash_chromatography", "purification", purification_ids,
            f"Flash chromatography: {purification_text}.", product_page_evidence,
            method="original_extraction",
        ),
    ]
    reaction["procedure"]["general_procedure_reference"] = direct(
        "Prepared according to the general procedure A.", "general procedure A", product_page_evidence
    )
    reaction["procedure"]["addition_order"] = inherited(
        "N-Tosyl alpha-AAP and diethyl malonate were added to the electrolyte solution; the solution was degassed before electrolysis.",
        ["substrate", "diethyl malonate", "electrolyte solution", "argon degassing", "electrolysis"],
        general_ev,
        "The product paragraph explicitly invokes General Procedure A; the ordered operations are inherited from it.",
    )
    reaction["procedure"]["pretreatment"] = inherited(
        "degassed under argon for 15 min", "argon degassing for 15 min", general_ev,
        "The product paragraph explicitly invokes General Procedure A, which explicitly reports this pretreatment.",
    )
    reaction["procedure"]["modifications"] = direct(
        f"Flash Chromatography: {purification_text}",
        {"purification": "flash chromatography", "eluent": purification_text},
        product_page_evidence,
    )
    reaction["special_conditions"]["electrodes"] = inherited(
        "GC working electrode ca. 120 x 10 x 2 mm3; SCE reference; Pt-gauze counter electrode; divided cell",
        {"working": "glassy carbon stick ca. 120 x 10 x 2 mm", "reference": "SCE", "counter": "Pt gauze"},
        apparatus_ev,
        "The product paragraph invokes General Procedure A; the apparatus is inherited from the article's general electrochemical setup.",
    )
    scheme_charge = copy.deepcopy(reaction["special_conditions"]["current"])
    general_charge = inherited(
        "usually after 4-5 F/mol of starting material",
        {"minimum": 4, "maximum": 5, "unit": "F/mol", "raw_text": "4-5 F/mol"},
        general_ev,
        "The product paragraph invokes General Procedure A, which reports a usual 4-5 F/mol endpoint.",
    )
    reaction["special_conditions"]["current"] = st.fact(
        [scheme_charge["raw_value"], general_charge["raw_value"]],
        None,
        list(dict.fromkeys(scheme_charge["evidence_ids"] + general_charge["evidence_ids"])),
        status="conflict",
        alternatives=[scheme_charge, general_charge],
    )
    unresolved = "Scheme 6 reports 4 F/mol, while invoked General Procedure A reports a usual 4-5 F/mol endpoint; both are retained without selecting one."
    if unresolved not in reaction["review"]["unresolved"]:
        reaction["review"]["unresolved"].append(unresolved)


def build_submission(init_path: Path):
    dataset = load(RUN / "dataset.json")
    old = dataset.get("standard")
    if not old:
        raise ValueError("Existing v3.1 records are required for evidence-preserving incremental update")
    bundle = load(init_path)
    if bundle["documents"] != dataset["documents"]:
        raise ValueError("Fresh standard-init template does not match the augmented prepared dataset")

    molecules = {m["molecule_id"]: copy.deepcopy(m) for m in old["molecules"]}
    reactions = {r["reaction_id"]: copy.deepcopy(r) for r in old["reactions"]}
    selected = [reactions[rid] for rid in SELECTED_RIDS]

    product_details = {
        "m28": {"mg": 49, "evidence": [P[6], IMG[6], P[7], IMG[7]], "purification": "EtOAc 100", "agents": ["a_etoac"]},
        "m29": {"mg": 36, "evidence": [P[7], IMG[7]], "purification": "EtOAc 100", "agents": ["a_etoac"]},
        "m30": {"mg": 65, "evidence": [P[7], IMG[7]], "purification": "EtOAc 100", "agents": ["a_etoac"]},
        "m31": {"mg": 61, "evidence": [P[7], IMG[7]], "purification": "EtOAc (5% v/v NH3 aq.) 100", "agents": ["a_etoac", "a_nh3"]},
        "m32": {"mg": 39, "evidence": [P[7], IMG[7]], "purification": "EtOAc (5 wt.-% NH3)/PE 60:40", "agents": ["a_etoac", "a_nh3", "a_pe"]},
        "m33": {"mg": 38, "evidence": [P[7], IMG[7]], "purification": "EtOAc/PE = 50:50", "agents": ["a_etoac", "a_pe"]},
        "m34": {"mg": 56, "evidence": [P[7], IMG[7]], "purification": "EtOAc/PE = 50:50", "agents": ["a_etoac", "a_pe"]},
    }

    for index, reaction in enumerate(selected, start=28):
        mid = f"m{index}"
        details = product_details[mid]
        add_ids(reaction["evidence_ids"], P[7], IMG[7])
        add_fact_evidence(reaction["location"], P[7], IMG[7])
        add_fact_evidence(reaction["reaction_name"], P[7], IMG[7])
        yield_outcome = next(o for o in reaction["outcomes"] if o["kind"] == "yield" and mid in o["target_molecule_ids"])
        add_fact_evidence(yield_outcome["value"], *details["evidence"])
        add_fact_evidence(yield_outcome["yield_type"], *details["evidence"])
        if not any(o["kind"] == "product_amount" and mid in o["target_molecule_ids"] for o in reaction["outcomes"]):
            reaction["outcomes"].append(product_amount(mid, details["mg"], details["evidence"]))
        for outcome in reaction["outcomes"]:
            if outcome["kind"] == "conversion" and index in {28, 29, 31, 32}:
                add_fact_evidence(outcome["value"], P[7], IMG[7])
        add_full_general_procedure(
            reaction,
            details["evidence"],
            details["purification"],
            details["agents"],
        )
        old_unresolved = "The covered pages give scheme conditions/outcomes but not product-specific scale, time, work-up, purification, or product amount."
        reaction["review"]["unresolved"] = [u for u in reaction["review"]["unresolved"] if u != old_unresolved]

        molecule = molecules[mid]
        add_ids(molecule["occurrences"][0]["evidence_ids"], P[7], IMG[7])
        add_fact_evidence(molecule["occurrences"][0]["compound_label"], P[7], IMG[7])
        add_fact_evidence(molecule["names"]["original"], P[7], IMG[7])
        add_fact_evidence(molecule["names"]["english"], P[7], IMG[7])

    # The product-32 heading conflicts with Scheme 6 and its own F-containing
    # analytical formula/19F NMR.  Preserve all source expressions without
    # silently correcting the printed title.
    m32 = molecules["m32"]
    printed_name = direct(
        "Diethyl 1-(2-Methoxybenzylamino)-3-methylbut-2-enylphosphonate (32)",
        "Diethyl 1-(2-methoxybenzylamino)-3-methylbut-2-enylphosphonate",
        [P[7], IMG[7]],
    )
    structure_name = st.fact(
        "Scheme 6 draws a para-CF3 benzyl substituent; the product paragraph reports 19F NMR and HRMS for C17H26F3NO4P.",
        "Diethyl 1-[(4-trifluoromethylbenzyl)amino]-3-methylbut-2-enylphosphonate",
        [P[4], IMG[4], P[7], IMG[7]],
        method="structure_derivation",
        derivation="Candidate English identity follows the Scheme 6 drawing and the F3 analytical formula/19F NMR, which conflict with the printed page-7 heading.",
    )
    m32["names"]["original"] = printed_name
    m32["names"]["english"] = st.fact(
        [printed_name["raw_value"], structure_name["raw_value"]],
        None,
        [P[4], IMG[4], P[7], IMG[7]],
        status="conflict",
        alternatives=[printed_name, structure_name],
    )
    m32["identity_note"] = st.fact(
        "The printed product-32 heading says 2-methoxybenzyl, whereas Scheme 6 shows para-CF3 and the same paragraph reports 19F NMR/HRMS for an F3 formula.",
        None,
        [P[4], IMG[4], P[7], IMG[7]],
        status="conflict",
        alternatives=[printed_name, structure_name],
    )
    scope32 = reactions[f"{DOC}:scope32"]
    conflict_text = "Product 32 has a source identity conflict: page-7 heading says 2-methoxybenzyl, while Scheme 6 and the page-7 19F/HRMS data indicate para-CF3; no source value is silently corrected."
    if conflict_text not in scope32["review"]["unresolved"]:
        scope32["review"]["unresolved"].append(conflict_text)
    scope32["notes"].append(st.fact(
        [printed_name["raw_value"], structure_name["raw_value"]],
        None,
        [P[4], IMG[4], P[7], IMG[7]],
        status="conflict",
        alternatives=[printed_name, structure_name],
    ))

    required_mids = set()
    for reaction in selected:
        required_mids.update(p["molecule_id"] for p in reaction["participants"])
    bundle["molecules"] = [molecules[mid] for mid in sorted(required_mids)]
    bundle["reactions"] = selected
    bundle["source_context"] = [
        {
            "context_id": f"{DOC}:context_characterization_page7",
            "document_id": DOC,
            "evidence_ids": [P[7], IMG[7]],
            "kind": "compound_preparation_and_characterization",
            "text": "Product paragraphs for 28-34 provide masses, yields, conversions where stated, purification eluents, Rf, NMR and HRMS. Preparation facts are linked to existing scope records; spectral peak lists remain characterization context, not new reactions.",
        },
        {
            "context_id": f"{DOC}:context_references_pages7_8",
            "document_id": DOC,
            "evidence_ids": [P[7], IMG[7], P[8], IMG[8]],
            "kind": "bibliography_and_received_date",
            "text": "Acknowledgments, keywords, references 1-29, and Received: April 19, 2020. These regions contain no new current-paper synthesis execution.",
        },
    ]
    dump(RUN / "submitted.json", bundle)

    stereo = {}
    for mid in [f"m{i}" for i in range(28, 35)]:
        smiles = st.value(molecules[mid]["smiles"])
        stereo[mid] = st.stereo_inspect(smiles) if smiles else {"error": "SMILES missing"}
    dump(RUN / "raw_model/stereo_checks_pages7_8.json", stereo)

    coverage = {
        "document_id": DOC,
        "semantic_accuracy_verified": False,
        "pages": [
            {
                "page": 7,
                "state": "attempted",
                "regions": [
                    {
                        "region": "top-left continuation of product 28 paragraph",
                        "finding": "Completes the 81% conversion statement and gives characterization for the already-recorded scope-28 preparation; no separate experiment is created.",
                        "evidence_ids": [P[7], IMG[7]],
                    },
                    {
                        "region": "product preparation paragraphs 29-34",
                        "finding": "Prepared according to General Procedure A; product masses, isolated yields, conversions where printed, and flash-chromatography eluents are appended to existing Scheme-6 scope records 29-34.",
                        "evidence_ids": [P[7], IMG[7]],
                    },
                    {
                        "region": "product 32 heading and analytical block",
                        "finding": "Source conflict preserved: the heading names 2-methoxybenzyl, while Scheme 6 and the F3 HRMS/19F NMR indicate a para-CF3 benzyl structure.",
                        "evidence_ids": [P[7], IMG[7], P[4], IMG[4]],
                    },
                    {
                        "region": "Rf, NMR and HRMS blocks for products 28-34",
                        "finding": "Characterization evidence only; spectra and mass peaks are retained in source_context and are not counted as synthesis experiments.",
                        "evidence_ids": [P[7], IMG[7]],
                    },
                    {
                        "region": "acknowledgments, keywords, and references 1-14",
                        "finding": "Non-experimental article matter and cited literature; no new current-paper experiment.",
                        "evidence_ids": [P[7], IMG[7]],
                    },
                ],
                "record_ids": SELECTED_RIDS,
                "unresolved": [
                    "Product 32 identity conflict between the printed heading and Scheme 6/analytical data requires expert review.",
                    "No R/S assignment is printed for the product alpha stereocenters; unspecified is not interpreted as achiral or racemic.",
                    "No SI file was explicitly supplied or associated.",
                ],
            },
            {
                "page": 8,
                "state": "attempted",
                "regions": [
                    {
                        "region": "references 14-29 and received-date line",
                        "finding": "Bibliographic continuation only. It contains cited literature and the received date, not a current-paper reaction arrow, table row, procedure, control, or mechanism experiment.",
                        "evidence_ids": [P[8], IMG[8]],
                    }
                ],
                "record_ids": [],
                "unresolved": [
                    "References are not expanded into reaction records because page 8 reports no executable experiments from this paper."
                ],
            },
        ],
    }
    dump(RUN / "coverage/pages_0007_0008.json", coverage)
    dump(RUN / "raw_model/page_visual_inventory_pages7_8.json", {
        "document_id": DOC,
        "visually_inspected_pages": [7, 8],
        "image_detail": "original",
        "pages": [
            {"page": 7, "regions": ["product-28 continuation", "product paragraphs 29-34", "characterization blocks", "acknowledgments/keywords/references 1-14"]},
            {"page": 8, "regions": ["references 14-29", "received-date line", "blank remainder"]},
        ],
    })
    print(json.dumps({
        "incremental_molecules": len(bundle["molecules"]),
        "updated_existing_reactions": len(bundle["reactions"]),
        "new_reactions": 0,
        "coverage_pages": [7, 8],
    }, ensure_ascii=False))


def write_audit_artifacts():
    from reaction_tools import dataset_audit

    dataset = load(RUN / "dataset.json")
    audit = dataset_audit(dataset)
    mentions = st.source_mentions(
        dataset,
        [
            "Prepared according to the general procedure A",
            "C17H26F3NO4P",
            "19F NMR",
            "Received: April 19, 2020",
        ],
        context=220,
    )
    dump(RUN / "quality_checks.json", audit)
    dump(RUN / "raw_model/audit_pages7_8.json", audit)
    dump(RUN / "raw_model/source_mentions_pages7_8.json", mentions)
    warning_types = sorted({
        warning.split(":", 1)[0]
        for item in audit["warnings"]
        for warning in item["warnings"]
    })
    print(json.dumps({
        "audit_warning_records": len(audit["warnings"]),
        "audit_warning_types": warning_types,
        "source_mentions": len(mentions["matches"]),
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("augment-evidence")
    sub.add_parser("audit-artifacts")
    build = sub.add_parser("build")
    build.add_argument("--init", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "augment-evidence":
        augment_evidence()
    elif args.command == "audit-artifacts":
        write_audit_artifacts()
    else:
        build_submission(args.init)


if __name__ == "__main__":
    main()
