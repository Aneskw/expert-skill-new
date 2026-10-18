#!/usr/bin/env python3
"""Persist strict validation and audit results after the page-13--16 merge."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

run = ROOT / "runs/current"
data = rt.read(run / "dataset.json")
strict_ids = {record["record_id"] for record in data["experiments"]}
errors = rt.validate(data, run, strict_record_ids=strict_ids)
errors.extend(st.validate(data["standard"], run))
validation = {
    "valid": not errors,
    "errors": errors,
    "standard_v3_present": "standard" in data,
    "chemical_accuracy_verified": False,
}
audit = rt.dataset_audit(data)
st.write(run / "raw_model/pages_0013_0016_strict_validate.json", validation)
st.write(run / "raw_model/pages_0013_0016_audit.json", audit)
st.write(run / "quality_checks.json", audit)
print({"valid": validation["valid"], "audit_warning_records": len(audit["warnings"])})
