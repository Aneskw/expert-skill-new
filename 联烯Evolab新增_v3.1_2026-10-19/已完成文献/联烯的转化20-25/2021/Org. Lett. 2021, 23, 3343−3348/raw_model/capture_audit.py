#!/usr/bin/env python3
"""Persist advisory audit and validator results without modifying records."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records

dataset_path = ROOT / "runs/current/dataset.json"
bundle_path = ROOT / "runs/current/standard_records.json"
dataset = json.loads(dataset_path.read_text())
bundle = json.loads(bundle_path.read_text())
standard_errors = standard_records.validate(bundle, bundle_path.parent)
strict_errors = reaction_tools.validate(
    dataset, dataset_path.parent,
    strict_record_ids={r["record_id"] for r in dataset["experiments"]},
)
audit = reaction_tools.dataset_audit(dataset)
(ROOT / "runs/current/raw_model/audit.json").write_text(
    json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
(ROOT / "runs/current/raw_model/validation.json").write_text(json.dumps({
    "standard_valid": not standard_errors,
    "standard_errors": standard_errors,
    "strict_dataset_valid": not strict_errors,
    "strict_dataset_errors": strict_errors,
    "chemical_accuracy_verified": False,
}, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({
    "standard_errors": len(standard_errors),
    "strict_errors": len(strict_errors),
    "audit_flagged_records": len(audit["warnings"]),
}, ensure_ascii=False))
