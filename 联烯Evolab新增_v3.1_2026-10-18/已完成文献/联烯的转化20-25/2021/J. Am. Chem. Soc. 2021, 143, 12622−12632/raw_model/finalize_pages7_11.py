#!/usr/bin/env python3
"""Persist the post-merge audit/counts and close page tasks 7--11."""
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt
import standard_records as st

data = st.read(ROOT / "runs/current/dataset.json")
audit = rt.dataset_audit(data)
st.write(ROOT / "runs/current/audit.json", audit)

tasks_path = ROOT / "runs/current/tasks.json"
tasks = st.read(tasks_path)
for task in tasks:
    match = re.search(r":page(\d+):", task.get("evidence_id", ""))
    if match and 7 <= int(match.group(1)) <= 11:
        task["state"] = "completed"
        task["reason"] = "Original page image visually inspected; decisions and record links are in coverage/pages_0007_0011.json."
st.write(tasks_path, tasks)

standard = data["standard"]
warning_types = Counter()
new_warning_types = Counter()
new_warning_records = []
for row in audit["warnings"]:
    is_new = ":scheme8_" in row["record_id"]
    if is_new:
        new_warning_records.append(row["record_id"])
    for warning in row["warnings"]:
        kind = warning.split(":", 1)[0]
        warning_types[kind] += 1
        if is_new:
            new_warning_types[kind] += 1

counts = {
    "documents": len(data["documents"]),
    "molecules": len(standard["molecules"]),
    "reactions": len(standard["reactions"]),
    "new_pages7_11_reactions": sum(":scheme8_" in r["reaction_id"] for r in standard["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in standard["reactions"]),
    "literature_examples": sum(r["record_type"] == "literature_example" for r in standard["reactions"]),
    "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in standard["reactions"]),
    "needs_review": sum(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "warning_records": len(audit["warnings"]),
    "warning_types": dict(warning_types),
    "new_pages7_11_warning_records": len(new_warning_records),
    "new_pages7_11_warning_types": dict(new_warning_types),
    "reaction_types": dict(Counter(r["record_type"] for r in standard["reactions"])),
    "experiment_categories": dict(Counter(
        c for r in standard["reactions"]
        for c in (st.value(r["experiment_categories"]) or []))),
}
st.write(ROOT / "runs/current/raw_model/final_counts.json", counts)
print(counts)
