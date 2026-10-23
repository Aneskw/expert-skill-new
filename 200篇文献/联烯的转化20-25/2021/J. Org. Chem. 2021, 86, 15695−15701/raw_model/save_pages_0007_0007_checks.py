#!/usr/bin/env python3
"""Persist strict validation and advisory audit after the page-7 delta."""
from __future__ import annotations

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
strict_ids = {record["record_id"] for record in data["experiments"]}
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
    "page_0007_0007_delta": {
        "page_range": [7, 7],
        "new_reactions": 0,
        "new_molecules": 0,
        "new_source_context_items": 1,
        "semantic_accuracy_verified": False,
    },
    "audit": audit,
}

st.write(base / "quality_checks.json", result)
st.write(base / "raw_model/pages_0007_0007_audit.json", audit)
print(
    json.dumps(
        {
            "valid": result["strict_dataset_validation"]["valid"],
            "warnings": sum(len(item["warnings"]) for item in audit["warnings"]),
            "flagged_records": len(audit["warnings"]),
        },
        ensure_ascii=False,
    )
)
