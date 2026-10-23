"""Persist deterministic audit and mark this run's five prepared page tasks attempted."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit

run = Path(__file__).resolve().parents[1]
dataset = json.loads((run / "dataset.json").read_text())
audit = dataset_audit(dataset)
(run / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

tasks_path = run / "tasks.json"
tasks = json.loads(tasks_path.read_text())
for task in tasks:
    if task["document_id"] == "doc_1cf04303b61f633d" and task["task_id"].endswith(":text"):
        task["state"] = "attempted"
        task["reason"] = "Original page image inspected; extraction decisions are documented in coverage/pages_0001_0005.json and runs/workflow.md."
tasks_path.write_text(json.dumps(tasks, ensure_ascii=False, indent=2) + "\n")

issues = {}
for item in audit["warnings"]:
    for warning in item["warnings"]:
        key = warning.split(":", 1)[0]
        issues[key] = issues.get(key, 0) + 1
print(json.dumps({"warning_records": len(audit["warnings"]), "warning_instances": issues}, ensure_ascii=False))
