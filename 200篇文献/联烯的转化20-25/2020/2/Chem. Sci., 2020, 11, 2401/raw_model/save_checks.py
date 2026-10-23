#!/usr/bin/env python3
"""Persist deterministic validation/audit results without altering records."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records

dataset_path = ROOT / "runs/current/dataset.json"
standard_path = ROOT / "runs/current/standard_records.json"
dataset = json.loads(dataset_path.read_text())
standard = json.loads(standard_path.read_text())

strict_errors = reaction_tools.validate(
    dataset,
    dataset_path.parent,
    strict_record_ids={r["record_id"] for r in dataset["experiments"]},
)
standard_errors = standard_records.validate(standard, standard_path.parent)
audit = reaction_tools.dataset_audit(dataset)

standard_records.write(ROOT / "runs/current/raw_model/validation.json", {
    "standard_valid": not standard_errors,
    "standard_errors": standard_errors,
    "strict_dataset_valid": not strict_errors,
    "strict_dataset_errors": strict_errors,
    "chemical_accuracy_verified": False,
})
standard_records.write(ROOT / "runs/current/raw_model/audit.json", audit)
print(json.dumps({
    "standard_errors": len(standard_errors),
    "strict_errors": len(strict_errors),
    "audit_warning_records": len(audit["warnings"]),
}, ensure_ascii=False))
