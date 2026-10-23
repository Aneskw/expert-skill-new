#!/usr/bin/env python3
"""Save deterministic validation/audit output beside raw session artifacts."""

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
strict_ids = {r["record_id"] for r in dataset["experiments"]}
legacy_errors = reaction_tools.validate(dataset, run, strict_record_ids=strict_ids)
standard_errors = standard_records.validate(dataset["standard"], run)
validation = {
    "valid": not legacy_errors and not standard_errors,
    "legacy_strict_errors": legacy_errors,
    "standard_v3_errors": standard_errors,
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
audit = reaction_tools.dataset_audit(dataset)
(run / "raw_model/validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2) + "\n")
(run / "raw_model/audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"valid": validation["valid"], "audit_warning_records": len(audit["warnings"])}, ensure_ascii=False))
