#!/usr/bin/env python3
"""Persist validation/audit reports because the harness commands print only."""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "runs/current"
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

from reaction_tools import dataset_audit, read as read_dataset, validate as legacy_validate  # noqa: E402
from standard_records import validate as standard_validate, write  # noqa: E402

dataset = read_dataset(OUT / "dataset.json")
standard = json.loads((OUT / "standard_records.json").read_text())

standard_errors = standard_validate(standard, OUT)
strict_errors = legacy_validate(
    dataset,
    OUT,
    strict_record_ids={r["record_id"] for r in dataset["experiments"]},
)

write(OUT / "standard_validate.json", {
    "valid": not standard_errors,
    "errors": standard_errors,
    "chemical_accuracy_verified": False,
})
write(OUT / "strict_validate.json", {
    "valid": not strict_errors,
    "errors": strict_errors,
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
})
write(OUT / "audit.json", dataset_audit(dataset))

print(json.dumps({
    "standard_valid": not standard_errors,
    "strict_valid": not strict_errors,
    "audit_warning_records": len(dataset_audit(dataset)["warnings"]),
}, ensure_ascii=False))
