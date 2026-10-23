#!/usr/bin/env python3
"""Persist strict-validation and audit artifacts for pages 25--28."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from reaction_tools import dataset_audit, validate  # noqa: E402
from standard_records import validate as standard_validate  # noqa: E402


dataset = json.loads((RUN / "dataset.json").read_text())
strict_ids = {record["record_id"] for record in dataset["experiments"]}
errors = validate(dataset, RUN, strict_record_ids=strict_ids)
errors += standard_validate(dataset["standard"], RUN)
strict_result = {
    "valid": not errors,
    "errors": sorted(set(errors)),
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
(RUN / "raw_model" / "pages_0025_0028_strict_validate.json").write_text(
    json.dumps(strict_result, ensure_ascii=False, indent=2) + "\n"
)

audit = dataset_audit(dataset)
encoded = json.dumps(audit, ensure_ascii=False, indent=2) + "\n"
(RUN / "quality_checks.json").write_text(encoded)
(RUN / "raw_model" / "pages_0025_0028_audit.json").write_text(encoded)

raw = json.loads((RUN / "raw_model" / "pages_0025_0028_raw_extraction.json").read_text())
batch_ids = {
    record_id
    for page in raw["pages"]
    for record_id in page.get("linked_record_ids", [])
}
batch_warnings = [row for row in audit["warnings"] if row["record_id"] in batch_ids]
summary = {
    "valid": strict_result["valid"],
    "record_kind_counts": audit["record_kind_counts"],
    "full_warning_rows": len(audit["warnings"]),
    "batch_record_count": len(batch_ids),
    "batch_warning_rows": len(batch_warnings),
    "batch_warning_types": sorted(
        {
            warning.split(":", 1)[0]
            for row in batch_warnings
            for warning in row["warnings"]
        }
    ),
    "batch_reaction_identity_warnings": sum(
        any(warning.startswith("reaction_identity:") for warning in row["warnings"])
        for row in batch_warnings
    ),
    "batch_outcome_target_warnings": sum(
        any(warning.startswith("outcome_target:") for warning in row["warnings"])
        for row in batch_warnings
    ),
    "batch_measurement_basis_warnings": sum(
        any(warning.startswith("measurement_basis:") for warning in row["warnings"])
        for row in batch_warnings
    ),
    "batch_structure_gap_warnings": sum(
        any(warning.startswith("structure_gap:") for warning in row["warnings"])
        for row in batch_warnings
    ),
    "chemical_accuracy_verified": audit["chemical_accuracy_verified"],
}
(RUN / "raw_model" / "pages_0025_0028_audit_summary.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2) + "\n"
)
print(json.dumps(summary, ensure_ascii=False))
