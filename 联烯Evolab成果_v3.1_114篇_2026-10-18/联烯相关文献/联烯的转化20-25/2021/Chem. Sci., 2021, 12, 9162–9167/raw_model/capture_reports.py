#!/usr/bin/env python3
"""Capture final deterministic validation and advisory audit reports."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt
import standard_records as st

run = ROOT / "runs/current"
dataset = rt.read(run / "dataset.json")
errors = rt.validate(dataset, run, strict_record_ids={r["record_id"] for r in dataset["experiments"]})
if "standard" in dataset:
    errors += st.validate(dataset["standard"], run)
validation = {
    "valid": not errors,
    "errors": sorted(set(errors)),
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
audit = rt.dataset_audit(dataset)
counts = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        counts[warning.split(":", 1)[0]] += 1
audit["warning_summary"] = {
    "records_with_warnings": len(audit["warnings"]),
    "categories": dict(counts),
    "disposition": {
        "structure_gap": "Retained: compact main-text structures were not converted into guessed complete SMILES; labels, yields, conditions, and stereochemical uncertainty are preserved.",
        "measurement_basis": "Retained: Fig. 5 and Fig. 6a do not state isolated/NMR/other yield basis on the supplied pages.",
    },
}
st.write(run / "validate_strict.json", validation)
st.write(run / "audit.json", audit)
print(json.dumps({"validation": validation, "warning_summary": audit["warning_summary"]}, ensure_ascii=False, indent=2))
