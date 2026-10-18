#!/usr/bin/env python3
"""Apply the reproducible v3.1 outcome-container cleanup after page extension."""
import json
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(root / "skills/chemical-reaction-extraction/scripts"))
from standard_records import fact, new_participant
path = root / "runs/current/submitted.json"
bundle = json.loads(path.read_text())
removed = 0
for reaction in bundle["reactions"]:
    for outcome in reaction["outcomes"]:
        if "evidence_ids" in outcome:
            outcome.pop("evidence_ids")
            removed += 1

doc = "doc_2401da5b5fbe4542"
evidence = [f"{doc}:page15:image", f"{doc}:page15:text", f"{doc}:page16:image", f"{doc}:page16:text"]
new_ids = {"mol_13aa", "mol_13ab", "mol_13ba", "mol_13bb", "mol_13ca", "mol_13cb",
           "mol_cf3ph_snbu3", "mol_chlorosulfonyl_isocyanate", "mol_chlorobenzene",
           "mol_26lutidine", "mol_acetyl_chloride", "mol_acetic_acid"}
for molecule in bundle["molecules"]:
    if molecule["molecule_id"] not in new_ids:
        continue
    for stereo_fact in molecule["stereochemistry"].values():
        if (isinstance(stereo_fact, dict) and stereo_fact.get("normalized_value") is not None
                and stereo_fact.get("acquisition_method") in {"inference", "structure_derivation"}
                and not stereo_fact.get("derivation")):
            stereo_fact["derivation"] = "Standard stereo graph/source-drawing inspection; derived status is not source identity verification."

for reaction in bundle["reactions"]:
    # A recovered mass is an amount, not a percentage recovery outcome.
    if reaction["reaction_id"] == f"{doc}:rxn_carb_1bb":
        reaction["participants"] = [p for p in reaction["participants"]
                                    if not (p["molecule_id"] == "mol_7b" and p["role"] == "product")]
        for item in reaction["outcomes"]:
            if item["target_molecule_ids"] == ["mol_7b"] and item["value"].get("normalized_value", {}).get("unit") == "mg":
                item["kind"] = "product_amount"
                item["yield_type"]["raw_value"] = "recovered starting material amount"
                item["yield_type"]["normalized_value"] = "recovered starting material amount"
    # Selectivity targets must also be reaction participants.
    for item in reaction["outcomes"]:
        if item["kind"] != "selectivity":
            continue
        for target in item["target_molecule_ids"]:
            if target not in new_ids or any(p["molecule_id"] == target for p in reaction["participants"]):
                continue
            participant_id = f"{reaction['reaction_id']}:{target}:{len(reaction['participants']) + 1}"
            p = new_participant(participant_id, target, "product", evidence)
            p["product_role"] = fact("component identity of isolated E/Z mixture", "component identity of isolated E/Z mixture", evidence)
            reaction["participants"].append(p)
    if reaction["reaction_id"] == f"{doc}:rxn_scheme7_9_19":
        for p in reaction["participants"]:
            if p["molecule_id"] == "mol_benzoyl_acetonitrile":
                old = p["participant_id"]
                p["molecule_id"] = "mol_19"
                p["participant_id"] = old.replace("mol_benzoyl_acetonitrile", "mol_19")
                for s in reaction["steps"]:
                    s["participant_ids"] = [p["participant_id"] if x == old else x for x in s["participant_ids"]]
    if reaction["reaction_id"] in {f"{doc}:rxn_scheme9_39_pph3_catalytic", f"{doc}:rxn_scheme9_39_pph3_stoichiometric"}:
        e17 = [f"{doc}:page17:image", f"{doc}:page17:text"]
        stoichiometric = reaction["reaction_id"].endswith("stoichiometric")
        for p in reaction["participants"]:
            if p["molecule_id"] == "mol_39":
                p["charges"]["amount"] = fact("0.350 mmol", {"number": 0.350, "unit": "mmol"}, e17)
                p["charges"]["equivalents"] = fact("1.0 equiv", {"number": 1.0, "unit": "equiv"}, e17)
                p["equivalence_reference"] = fact("mol_39", "mol_39", e17)
            elif p["molecule_id"] == "mol_pph3":
                eq = 1.0 if stoichiometric else 0.2
                p["charges"]["equivalents"] = fact(f"{eq:g} equiv", {"number": eq, "unit": "equiv"}, e17)
                p["equivalence_reference"] = fact("mol_39", "mol_39", e17)
            elif p["molecule_id"] == "mol_mecn":
                p["charges"]["volume"] = fact("10 mL", {"number": 10, "unit": "mL"}, e17)
        for note in reaction["notes"]:
            if isinstance(note.get("raw_value"), str) and note["raw_value"].startswith("Detailed control:"):
                text = "The detailed page-17 control quantities and the catalytic/stoichiometric PPh3 distinction are structured in participants; TLC, 1H NMR and 31P NMR support the author-reported no-reaction outcome."
                note["raw_value"] = text
                note["normalized_value"] = text
path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"removed_invalid_outcome_container_fields": removed}))
