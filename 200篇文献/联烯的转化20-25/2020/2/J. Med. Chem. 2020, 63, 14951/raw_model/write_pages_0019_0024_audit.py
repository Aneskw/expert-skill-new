#!/usr/bin/env python3
"""Regenerate final strict-validation and audit artifacts for pages 19--24."""

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
strict_ids = {r["record_id"] for r in dataset["experiments"]}
errors = validate(dataset, RUN, strict_record_ids=strict_ids)
errors += standard_validate(dataset["standard"], RUN)
strict_result = {
    "valid": not errors,
    "errors": sorted(set(errors)),
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
(RUN / "raw_model" / "pages_0019_0024_strict_validate.json").write_text(
    json.dumps(strict_result, ensure_ascii=False, indent=2) + "\n"
)

audit = dataset_audit(dataset)
encoded = json.dumps(audit, ensure_ascii=False, indent=2) + "\n"
(RUN / "quality_checks.json").write_text(encoded)
(RUN / "raw_model" / "pages_0019_0024_audit.json").write_text(encoded)

inventory = json.loads((RUN / "raw_model" / "pages_0019_0024_inventory.json").read_text())
batch_ids = set(inventory["new_reaction_ids"] + inventory["updated_existing_reaction_ids"])
batch_warnings = [w for w in audit["warnings"] if w["record_id"] in batch_ids]
print(json.dumps({
    "valid": strict_result["valid"],
    "record_kind_counts": audit["record_kind_counts"],
    "full_warning_rows": len(audit["warnings"]),
    "batch_warning_rows": len(batch_warnings),
    "batch_warning_types": sorted({warning.split(":", 1)[0] for row in batch_warnings for warning in row["warnings"]}),
    "batch_reaction_identity_warnings": sum(any(w.startswith("reaction_identity:") for w in row["warnings"]) for row in batch_warnings),
    "batch_outcome_target_warnings": sum(any(w.startswith("outcome_target:") for w in row["warnings"]) for row in batch_warnings),
    "batch_measurement_basis_warnings": sum(any(w.startswith("measurement_basis:") for w in row["warnings"]) for row in batch_warnings),
    "batch_structure_gap_warnings": sum(any(w.startswith("structure_gap:") for w in row["warnings"]) for row in batch_warnings),
    "chemical_accuracy_verified": audit["chemical_accuracy_verified"],
}, ensure_ascii=False))
