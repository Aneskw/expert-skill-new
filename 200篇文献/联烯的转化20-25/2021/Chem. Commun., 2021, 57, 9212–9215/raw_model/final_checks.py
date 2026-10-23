#!/usr/bin/env python3
"""Persist advisory audit and reproducible source-specific integrity checks."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))

from reaction_tools import dataset_audit  # noqa: E402


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


data = json.loads((RUN / "dataset.json").read_text())
standard = data["standard"]
audit = dataset_audit(data)
dump(RUN / "audit.json", audit)

checks = {
    "table1_entry2_3a": (
        RUN / "records/doc_eb3e0b8d7abf9bf7__table1_entry2_3a.md",
        ["0.1 mmol", "0.11 mmol", "1.1 equiv.", "10 mol%", "0.12 mmol", "1.2 equiv.", "0.5 mL", "85%", "same 85% standard experiment and are merged"],
    ),
    "scheme2a_products7_8": (
        RUN / "records/doc_eb3e0b8d7abf9bf7__scheme2a_products7_8.md",
        ["0.1 mmol", "0.11 mmol", "0.12 mmol", "0.5 mL", "85% (7+8)", "7:8 = 1:1", "standard conditions"],
    ),
    "scheme3_product19": (
        RUN / "records/doc_eb3e0b8d7abf9bf7__scheme3_product19.md",
        ["benzoyl fluoride", "67%", "yield measurement basis is not stated", "Electrophile amount"],
    ),
    "scheme3_product20": (
        RUN / "records/doc_eb3e0b8d7abf9bf7__scheme3_product20.md",
        ["ethyl bromoacetate", "75%", "yield measurement basis is not stated", "reaction time are not reported"],
    ),
}
verification = {"records": {}, "global": {}}
for key, (path, tokens) in checks.items():
    text = path.read_text()
    verification["records"][key] = {
        "markdown": str(path.relative_to(ROOT)),
        "required_tokens": tokens,
        "missing_tokens": [token for token in tokens if token not in text],
    }

coverage = json.loads((RUN / "coverage/pages_0001_0004.json").read_text())
reaction_ids = [r["reaction_id"] for r in standard["reactions"]]
molecule_ids = [m["molecule_id"] for m in standard["molecules"]]
verification["global"] = {
    "reaction_count": len(reaction_ids),
    "specific_experiment_count": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
    "molecule_count": len(molecule_ids),
    "unique_reaction_ids": len(reaction_ids) == len(set(reaction_ids)),
    "unique_molecule_ids": len(molecule_ids) == len(set(molecule_ids)),
    "all_machine_records_need_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "all_molecule_stereo_statuses_present": all(
        all(m["stereochemistry"][k]["normalized_value"] is not None for k in ("chirality_status", "geometry_status", "representation_status"))
        for m in standard["molecules"]
    ),
    "coverage_pages_exact": [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4],
    "coverage_regions_nonempty": all(p["regions"] for p in coverage["pages"]),
    "audit_warning_records": len(audit["warnings"]),
    "audit_warning_categories": sorted({w.split(":", 1)[0] for row in audit["warnings"] for w in row["warnings"]}),
}
dump(RUN / "raw_model/verification_summary.json", verification)

tasks = json.loads((RUN / "tasks.json").read_text())
for task in tasks:
    if task.get("document_id") == "doc_eb3e0b8d7abf9bf7" and task.get("evidence_id") in {
        "doc_eb3e0b8d7abf9bf7:page1:text", "doc_eb3e0b8d7abf9bf7:page2:text",
        "doc_eb3e0b8d7abf9bf7:page3:text", "doc_eb3e0b8d7abf9bf7:page4:text",
    }:
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; page inventory and record mapping saved in coverage/pages_0001_0004.json."
dump(RUN / "tasks.json", tasks)
print(json.dumps(verification["global"], ensure_ascii=False, indent=2))
