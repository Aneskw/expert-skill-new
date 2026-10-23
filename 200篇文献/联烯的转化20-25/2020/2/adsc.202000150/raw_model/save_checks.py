#!/usr/bin/env python3
"""Save deterministic validation and audit snapshots for this extraction run."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records

run = ROOT / "runs/current"
dataset = reaction_tools.read(run / "dataset.json")
standard = standard_records.read(run / "standard_records.json")
standard_errors = standard_records.validate(standard, run)
strict_errors = reaction_tools.validate(dataset, run, strict_record_ids={r["record_id"] for r in dataset["experiments"]})
audit = reaction_tools.dataset_audit(dataset)

snapshots = {
    "standard_validation_pages_0001_0006.json": {
        "valid": not standard_errors,
        "errors": standard_errors,
        "chemical_accuracy_verified": False,
    },
    "strict_validation_pages_0001_0006.json": {
        "valid": not strict_errors,
        "errors": strict_errors,
        "standard_v3_present": "standard" in dataset,
        "chemical_accuracy_verified": False,
    },
    "audit_pages_0001_0006.json": audit,
}
for name, payload in snapshots.items():
    (run / "raw_model" / name).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")

print(json.dumps({
    "standard_valid": not standard_errors,
    "strict_valid": not strict_errors,
    "warning_records": len(audit["warnings"]),
    "measurement_basis_warnings": sum(
        any(w.startswith("measurement_basis:") for w in row["warnings"])
        for row in audit["warnings"]
    ),
    "structure_gap_warnings": sum(
        any(w.startswith("structure_gap:") for w in row["warnings"])
        for row in audit["warnings"]
    ),
}, ensure_ascii=False))
