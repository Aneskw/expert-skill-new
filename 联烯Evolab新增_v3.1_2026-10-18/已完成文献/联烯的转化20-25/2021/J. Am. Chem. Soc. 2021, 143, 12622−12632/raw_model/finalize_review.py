#!/usr/bin/env python3
"""Write deterministic audit artifacts and extraction counts after standard merge."""
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

data = st.read(ROOT / "runs/current/dataset.json")
audit = rt.dataset_audit(data)
st.write(ROOT / "runs/current/audit.json", audit)

tasks_path = ROOT / "runs/current/tasks.json"
tasks = st.read(tasks_path)
for task in tasks:
    match = __import__("re").search(r":page(\d+):", task.get("evidence_id", ""))
    if match and 1 <= int(match.group(1)) <= 6:
        task["state"] = "completed"
        task["reason"] = "Original page image visually inspected; decisions and record links are in coverage/pages_0001_0006.json."
st.write(tasks_path, tasks)

standard = data["standard"]
warning_types = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        warning_types[warning.split(":", 1)[0]] += 1

counts = {
    "documents": len(data["documents"]),
    "molecules": len(standard["molecules"]),
    "reactions": len(standard["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
    "generic_schemes": sum(r["record_type"] == "generic_scheme" for r in standard["reactions"]),
    "literature_examples": sum(r["record_type"] == "literature_example" for r in standard["reactions"]),
    "proposed_mechanisms": sum(r["record_type"] == "proposed_mechanism" for r in standard["reactions"]),
    "needs_review": sum(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "warning_records": len(audit["warnings"]),
    "warning_types": dict(warning_types),
    "reaction_types": dict(Counter(r["record_type"] for r in standard["reactions"])),
    "experiment_categories": dict(Counter(
        c for r in standard["reactions"]
        for c in (st.value(r["experiment_categories"]) or []))),
}
st.write(ROOT / "runs/current/raw_model/final_counts.json", counts)
print(counts)
