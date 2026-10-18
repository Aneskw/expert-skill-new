#!/usr/bin/env python3
"""Save deterministic validation/audit outputs without modifying source records."""

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt  # noqa: E402
import standard_records as st  # noqa: E402

run = ROOT / "runs/current"
raw = run / "raw_model"
data = json.loads((run / "dataset.json").read_text())
standard = json.loads((run / "standard_records.json").read_text())

standard_errors = st.validate(standard, run)
strict_errors = rt.validate(data, run, strict_record_ids={r["record_id"] for r in data["experiments"]})
audit = rt.dataset_audit(data)

outputs = {
    "validation_standard.json": {
        "valid": not standard_errors,
        "errors": standard_errors,
        "chemical_accuracy_verified": False,
    },
    "validation_strict.json": {
        "valid": not strict_errors,
        "errors": strict_errors,
        "standard_v3_present": "standard" in data,
        "chemical_accuracy_verified": False,
    },
    "audit.json": audit,
}
for name, value in outputs.items():
    (raw / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")

warning_counts = {}
for row in audit["warnings"]:
    for warning in row["warnings"]:
        key = warning.split(":", 1)[0]
        warning_counts[key] = warning_counts.get(key, 0) + 1

print(json.dumps({
    "standard_valid": not standard_errors,
    "strict_valid": not strict_errors,
    "records": len(data["experiments"]),
    "source_reviewed_records": audit["source_reviewed_records"],
    "warning_record_count": len(audit["warnings"]),
    "warning_counts": warning_counts,
}, ensure_ascii=False, indent=2))
