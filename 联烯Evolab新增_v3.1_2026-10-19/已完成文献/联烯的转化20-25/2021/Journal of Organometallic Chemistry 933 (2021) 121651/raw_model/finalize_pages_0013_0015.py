#!/usr/bin/env python3
"""Persist strict-validation and audit outputs for the pages 13--15 batch."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(TOOLS))

from reaction_tools import dataset_audit, validate, write  # noqa: E402
from standard_records import validate as validate_standard  # noqa: E402


RUN = ROOT / "runs" / "current"
dataset = json.loads((RUN / "dataset.json").read_text())
strict_ids = {record["record_id"] for record in dataset["experiments"]}
errors = validate(dataset, RUN, strict_record_ids=strict_ids)
errors.extend(validate_standard(dataset["standard"], RUN))

strict_result = {
    "valid": not errors,
    "errors": errors,
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
audit_result = dataset_audit(dataset)

write(RUN / "raw_model" / "validate_strict_pages_0013_0015.json", strict_result)
write(RUN / "raw_model" / "audit_pages_0013_0015.json", audit_result)

print(
    json.dumps(
        {
            "strict_valid": strict_result["valid"],
            "audit_warning_records": len(audit_result["warnings"]),
        },
        ensure_ascii=False,
    )
)

if errors:
    raise SystemExit(1)
