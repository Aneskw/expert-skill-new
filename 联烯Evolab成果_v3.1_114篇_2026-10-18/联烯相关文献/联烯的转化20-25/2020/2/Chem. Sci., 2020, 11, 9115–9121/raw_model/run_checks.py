#!/usr/bin/env python3
"""Persist deterministic validation/audit results beside the extraction provenance."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

run = ROOT / "runs/current"
bundle = st.read(run / "standard_records.json")
data = rt.read(run / "dataset.json")
standard_errors = st.validate(bundle, run)
strict_errors = rt.validate(data, run, strict_record_ids={r["record_id"] for r in data["experiments"]})
if "standard" in data:
    strict_errors += st.validate(data["standard"], run)

standard_result = {"valid": not standard_errors, "errors": standard_errors, "chemical_accuracy_verified": False}
strict_result = {"valid": not strict_errors, "errors": sorted(set(strict_errors)), "standard_v3_present": "standard" in data, "chemical_accuracy_verified": False}
audit_result = rt.dataset_audit(data)
st.write(run / "raw_model/validation_standard.json", standard_result)
st.write(run / "raw_model/validation_strict.json", strict_result)
st.write(run / "raw_model/audit.json", audit_result)
print(json.dumps({
    "standard_valid": standard_result["valid"],
    "strict_valid": strict_result["valid"],
    "audit_warning_records": len(audit_result["warnings"]),
}, ensure_ascii=False))
