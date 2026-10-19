#!/usr/bin/env python3
"""Persist post-merge validation and advisory audit without changing records."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

base = ROOT / "runs/current"
data = json.loads((base / "dataset.json").read_text())
strict_ids = {r["record_id"] for r in data["experiments"]}
legacy_errors = rt.validate(data, base, strict_record_ids=strict_ids)
standard_errors = st.validate(data["standard"], base)
audit = rt.dataset_audit(data)
result = {
    "strict_dataset_validation": {
        "valid": not (legacy_errors or standard_errors),
        "legacy_errors": legacy_errors,
        "standard_errors": standard_errors,
        "standard_v3_present": "standard" in data,
        "chemical_accuracy_verified": False,
    },
    "audit": audit,
}
st.write(base / "quality_checks.json", result)
st.write(base / "raw_model/pages_0001_0006_audit.json", audit)
print(json.dumps({
    "valid": result["strict_dataset_validation"]["valid"],
    "warnings": sum(len(x["warnings"]) for x in audit["warnings"]),
    "flagged_records": len(audit["warnings"]),
}, ensure_ascii=False))
