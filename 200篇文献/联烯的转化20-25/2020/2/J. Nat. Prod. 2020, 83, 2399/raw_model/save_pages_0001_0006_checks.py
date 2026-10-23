#!/usr/bin/env python3
"""Persist read-only validation/audit results for the pages 1--6 increment."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt
import standard_records as st

dataset_path = ROOT / "runs/current/dataset.json"
dataset = json.loads(dataset_path.read_text())
errors = rt.validate(
    dataset,
    dataset_path.parent,
    strict_record_ids={record["record_id"] for record in dataset["experiments"]},
)
strict = {
    "valid": not errors,
    "errors": errors,
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
audit = rt.dataset_audit(dataset)
summary = {
    "warning_records": len(audit["warnings"]),
    "warning_kinds": sorted({warning.split(":", 1)[0]
                             for row in audit["warnings"]
                             for warning in row["warnings"]}),
    "source_reviewed_records": audit["source_reviewed_records"],
    "chemical_accuracy_verified": False,
}
st.write(ROOT / "runs/current/raw_model/pages_0001_0006_strict_validate.json", strict)
st.write(ROOT / "runs/current/raw_model/pages_0001_0006_audit.json", audit)
st.write(ROOT / "runs/current/raw_model/pages_0001_0006_audit_summary.json", summary)
st.write(ROOT / "runs/current/quality_checks.json", audit)
print(json.dumps({"strict": strict, "audit_summary": summary}, ensure_ascii=False))
