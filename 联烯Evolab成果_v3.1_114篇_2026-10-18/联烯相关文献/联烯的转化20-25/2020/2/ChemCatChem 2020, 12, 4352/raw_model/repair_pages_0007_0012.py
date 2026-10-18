#!/usr/bin/env python3
"""One-time normalization after extending the already-built page 1--6 bundle."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from standard_records import fact, new_molecule, new_participant, stereo_inspect  # noqa: E402


def F(raw=None, normalized=None, evidence=(), method="original_extraction", *, derivation=None):
    if normalized is None and raw is not None:
        normalized = raw
    return fact(raw, normalized, evidence, method, derivation=derivation)


def q(number, unit):
    return {"number": number, "unit": unit}


bundle_path = ROOT / "runs/current/submitted.json"
bundle = json.loads(bundle_path.read_text())
dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
bundle["documents"] = dataset["documents"]

# Deduplicate page-7/10 identities against labels already registered on pages
# 1--6.  Reactions keep one molecule registry identity and gain the new source
# occurrence/structure evidence.
molecules = {m["molecule_id"]: m for m in bundle["molecules"]}
dedup = {
    "mol_benzoyl_acetonitrile": "mol_19",
    "mol_s9": "mol_pmosnco",
    "mol_s10": "mol_pnosnco",
}
for duplicate_id, canonical_id in dedup.items():
    if duplicate_id not in molecules:
        continue
    duplicate = molecules[duplicate_id]
    canonical = molecules[canonical_id]
    canonical["occurrences"].extend(
        occurrence for occurrence in duplicate["occurrences"]
        if not any(existing["compound_label"]["normalized_value"] == occurrence["compound_label"]["normalized_value"]
                   and set(existing["evidence_ids"]) == set(occurrence["evidence_ids"])
                   for existing in canonical["occurrences"])
    )
    canonical["representation_type"] = duplicate["representation_type"]
    canonical["smiles"] = duplicate["smiles"]
    canonical["stereochemistry"] = duplicate["stereochemistry"]
    for reaction in bundle["reactions"]:
        for participant in reaction["participants"]:
            if participant["molecule_id"] == duplicate_id:
                participant["molecule_id"] = canonical_id
                participant["participant_id"] = participant["participant_id"].replace(duplicate_id, canonical_id)
        for step in reaction["steps"]:
            step["participant_ids"] = [pid.replace(duplicate_id, canonical_id) for pid in step["participant_ids"]]
        for outcome in reaction["outcomes"]:
            outcome["target_molecule_ids"] = [canonical_id if mid == duplicate_id else mid for mid in outcome["target_molecule_ids"]]
    bundle["molecules"] = [m for m in bundle["molecules"] if m["molecule_id"] != duplicate_id]
    del molecules[duplicate_id]


def ensure_simple(mid, label, name, smiles, evidence):
    if mid in molecules:
        return
    diagnostic = stereo_inspect(smiles)
    canonical = diagnostic["canonical_isomeric_smiles"]
    m = new_molecule(mid, "doc_2401da5b5fbe4542", evidence)
    m["identity_status"] = "source_identified"
    m["representation_type"] = "smiles"
    m["occurrences"][0]["compound_label"] = F(label, label, evidence)
    m["names"]["original"] = F(label, label, evidence)
    m["names"]["english"] = F(name, name, evidence)
    m["names"]["abbreviations"] = F([label], [label], evidence)
    m["names"]["aliases"] = fact(evidence=evidence, derivation="No additional alias reported.")
    m["cas_numbers"] = fact(evidence=evidence, derivation="CAS number not reported; none guessed.")
    m["generic_template"] = None
    m["smiles"] = F(smiles, canonical, evidence, "structure_derivation", derivation="Conservative structure derivation from an unambiguous source name; parsing is not source verification.")
    st = m["stereochemistry"]
    st["chirality_status"] = F("no stereogenic element in the named molecule", "achiral", evidence, "structure_derivation", derivation="Source name and graph inspected; machine record remains needs_review.")
    st["geometry_status"] = F("not applicable", "not_applicable", evidence, "structure_derivation", derivation="No E/Z, axial or ring relationship applies to this graph.")
    st["representation_status"] = F("complete graph representation", "encoded", evidence, "inference", derivation="The encoded graph contains no unresolved stereogenic element.")
    st["atom_centers"] = F(diagnostic["atom_centers"], diagnostic["atom_centers"], evidence, "structure_derivation", derivation="Standard stereo graph diagnostic.")
    st["double_bonds"] = F(diagnostic["double_bonds"], diagnostic["double_bonds"], evidence, "structure_derivation", derivation="Standard stereo graph diagnostic.")
    st["axial"] = fact(evidence=evidence)
    st["mixture"] = fact(evidence=evidence)
    st["assignment"] = fact(evidence=evidence)
    st["isomeric_smiles"] = F(canonical, canonical, evidence, "structure_derivation", derivation="Canonical isomeric graph; source identity remains needs_review.")
    st["identity_key"] = F(diagnostic["stereo_graph_key"], diagnostic["stereo_graph_key"], evidence, "structure_derivation", derivation="Toolkit graph key, not sample identity verification.")
    m["identity_note"] = F(name, name, evidence)
    bundle["molecules"].append(m)
    molecules[mid] = m


page_evidence = lambda *pages: [item for page in pages for item in (f"doc_2401da5b5fbe4542:page{page}:image", f"doc_2401da5b5fbe4542:page{page}:text")]
ensure_simple("mol_water", "H2O", "water", "O", page_evidence(9, 10, 11))
ensure_simple("mol_hexanes", "hexanes", "hexanes", "CCCCCC", page_evidence(9, 10, 11, 12))


def add_workup_participant(reaction_id, molecule_id, evidence, raw_volume=None, number=None, *, role="workup_agent", stage="workup"):
    reaction = next(r for r in bundle["reactions"] if r["reaction_id"] == reaction_id)
    existing = next((p for p in reaction["participants"] if p["molecule_id"] == molecule_id and p["stage"] == stage), None)
    if existing:
        return existing["participant_id"]
    pid = f"{reaction_id}:{molecule_id}:{len(reaction['participants']) + 1}"
    participant = new_participant(pid, molecule_id, role, evidence, stage=stage)
    if raw_volume is not None:
        participant["charges"]["volume"] = F(raw_volume, q(number, "mL"), evidence)
    reaction["participants"].append(participant)
    return pid


def link_to_step(reaction_id, phrase, participant_ids):
    reaction = next(r for r in bundle["reactions"] if r["reaction_id"] == reaction_id)
    target = next(step for step in reaction["steps"] if phrase in (step["operation"]["normalized_value"] or ""))
    target["participant_ids"] = list(dict.fromkeys(target["participant_ids"] + participant_ids))


# Explicit workup/extraction volumes are participants, not only prose.
rid = "doc_2401da5b5fbe4542:rxn_mbh_7a"
ids = [
    add_workup_participant(rid, "mol_nh4cl", page_evidence(10), "50 mL saturated aqueous NH4Cl", 50),
    add_workup_participant(rid, "mol_etoac", page_evidence(10), "3 × 50 mL", 150),
    add_workup_participant(rid, "mol_nacl", page_evidence(10), "50 mL brine", 50),
    add_workup_participant(rid, "mol_na2so4", page_evidence(10)),
]
link_to_step(rid, "quench with saturated NH4Cl", ids)

rid = "doc_2401da5b5fbe4542:rxn_mbh_7b"
ids = [
    add_workup_participant(rid, "mol_nh4cl", page_evidence(11), "10 mL saturated aqueous NH4Cl", 10),
    add_workup_participant(rid, "mol_etoac", page_evidence(11), "3 × 10 mL", 30),
    add_workup_participant(rid, "mol_water", page_evidence(11), "10 mL", 10),
    add_workup_participant(rid, "mol_nacl", page_evidence(11), "10 mL brine", 10),
    add_workup_participant(rid, "mol_na2so4", page_evidence(11)),
]
link_to_step(rid, "quench with saturated NH4Cl", ids)

rid = "doc_2401da5b5fbe4542:rxn_s7_from_s6"
water = add_workup_participant(rid, "mol_water", page_evidence(9, 10), "200 mL + 300 mL", 500)
hexanes = add_workup_participant(rid, "mol_hexanes", page_evidence(10), "200 mL + 2 × 200 mL", 600)
link_to_step(rid, "add water, separate", [water])
link_to_step(rid, "cool, filter through Celite", [water, hexanes])

rid = "doc_2401da5b5fbe4542:rxn_s8"
ids = [
    add_workup_participant(rid, "mol_water", page_evidence(10), "20 mL", 20),
    add_workup_participant(rid, "mol_etoac", page_evidence(10), "40 mL + 3 × 40 mL", 160),
    add_workup_participant(rid, "mol_nh4cl", page_evidence(10), "3 × 40 mL saturated aqueous NH4Cl", 120),
    add_workup_participant(rid, "mol_nacl", page_evidence(10), "20 mL brine", 20),
    add_workup_participant(rid, "mol_na2so4", page_evidence(10)),
]
link_to_step(rid, "add water and EtOAc", ids)

rid = "doc_2401da5b5fbe4542:rxn_allenylation_7c"
ids = [
    add_workup_participant(rid, "mol_et2o", page_evidence(11), "15 mL dilution + 4 × 15 mL extraction", 75),
    add_workup_participant(rid, "mol_nh4cl", page_evidence(11), "15 mL saturated aqueous NH4Cl", 15),
    add_workup_participant(rid, "mol_nacl", page_evidence(11), "15 mL brine", 15),
    add_workup_participant(rid, "mol_na2so4", page_evidence(11)),
]
link_to_step(rid, "cool to 0 °C, dilute with ether", ids)

# The Scheme 11 first arrow explicitly gives a solvent mixture concentration.
scheme11 = next(r for r in bundle["reactions"] if r["reaction_id"] == "doc_2401da5b5fbe4542:rxn_scheme11_3m_s12")
scheme11["solvent_conditions"]["composition"] = F("15% TFA in CH2Cl2", ["TFA", "CH2Cl2"], page_evidence(8))
scheme11["solvent_conditions"]["mixture_ratio"] = F("15% TFA", q(15, "%"), page_evidence(8))

category_map = {
    "proposed_mechanism": "other",
    "cited_literature": "other",
    "substrate_preparation": "other",
    "application_synthesis": "derivatization",
    "scope": "substrate_scope",
}
for reaction in bundle["reactions"]:
    fact = reaction["experiment_categories"]
    if any(e.startswith("doc_2401da5b5fbe4542:page") and
           int(e.split(":page", 1)[1].split(":", 1)[0]) in range(7, 13)
           for e in reaction["evidence_ids"]):
        raw = fact["raw_value"] or []
        fact["normalized_value"] = [category_map.get(category, category) for category in raw]

# Multiple apparatus temperatures belong in the raw/normalized descriptive
# value rather than a malformed single-temperature quantity.
s8 = next(r for r in bundle["reactions"] if r["reaction_id"] == "doc_2401da5b5fbe4542:rxn_s8")
distill = next(s for s in s8["steps"] if "distill under reduced pressure" in (s["operation"]["normalized_value"] or ""))
distill["temperature"]["raw_value"] = "head temperature 100-110 °C"
distill["temperature"]["normalized_value"] = {"minimum": 100, "maximum": 110, "unit": "°C"}

s9 = next(r for r in bundle["reactions"] if r["reaction_id"] == "doc_2401da5b5fbe4542:rxn_s9")
distill_steps = [s for s in s9["steps"] if "distill" in (s["operation"]["normalized_value"] or "")]
distill_steps[0]["temperature"]["raw_value"] = "head temperature 120 °C"
distill_steps[0]["temperature"]["normalized_value"] = {"number": 120, "unit": "°C"}
distill_steps[1]["temperature"]["raw_value"] = "head temperature 160-165 °C"
distill_steps[1]["temperature"]["normalized_value"] = {"minimum": 160, "maximum": 165, "unit": "°C"}

bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")

tasks_path = ROOT / "runs/current/tasks.json"
tasks = json.loads(tasks_path.read_text())
target_evidence = {
    f"doc_2401da5b5fbe4542:page{page}:text" for page in range(7, 13)
}
for task in tasks:
    if task.get("document_id") == "doc_2401da5b5fbe4542" and task.get("evidence_id") in target_evidence:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected at original detail; records and contexts mapped in coverage/pages_0007_0012.json. Machine records remain needs_review."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")
