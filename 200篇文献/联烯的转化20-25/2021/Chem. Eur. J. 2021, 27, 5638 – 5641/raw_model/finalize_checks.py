#!/usr/bin/env python3
"""Persist advisory audit and representative Markdown completeness checks."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

RUN = ROOT / "runs/current"
data = rt.read(RUN / "dataset.json")
audit = rt.dataset_audit(data)
rt.write(RUN / "audit.json", audit)
rt.write(RUN / "raw_model/source_mentions.json",
         st.source_mentions(data, ["1d", "2d", "1h", "product 4", "1p", "product 5"], context=220))
standard_errors = st.validate(rt.read(RUN / "standard_records.json"), RUN)
strict_errors = rt.validate(data, RUN, strict_record_ids={r["record_id"] for r in data["experiments"]})
rt.write(RUN / "standard_validation.json", {
    "valid": not standard_errors,
    "errors": standard_errors,
    "chemical_accuracy_verified": False,
})
rt.write(RUN / "strict_validation.json", {
    "valid": not strict_errors,
    "errors": strict_errors,
    "standard_v3_present": "standard" in data,
    "chemical_accuracy_verified": False,
})

checks = {
    "table1_entry9_nmr": ["0.1", "mmol", "2", "mol%", "1 atm", "DCM", "0.1", "M", "RT", "12", "h", "76", "1H-NMR", "1,3,5-trimethoxybenzene"],
    "table2_entry9": ["0.3", "mmol", "5", "mol%", "1 atm", "DCM", "0.1", "M", "0", "°C", "12", "h", "RT", "24", "h", "56", "isolated yield", "average of two runs"],
    "scheme3_product4": ["1", "equiv", "2", "equiv", "5", "mol%", "1 atm", "DCE", "80", "°C", "15", "h", "23", "syringe pump"],
    "table2_entry1": ["70", "isolated yield", "Cross-table duplicate merged", "0.1 mmol NMR run is separate"],
}
results = {}
for suffix, tokens in checks.items():
    path = RUN / "records" / f"doc_18eac4838543527a__{suffix}.md"
    text = path.read_text()
    results[suffix] = {
        "path": str(path.relative_to(ROOT)),
        "missing_tokens": [token for token in tokens if token not in text],
    }

warning_categories = {}
for row in audit["warnings"]:
    for warning in row["warnings"]:
        key = warning.split(":", 1)[0]
        warning_categories[key] = warning_categories.get(key, 0) + 1

summary = {
    "markdown_checks": results,
    "all_markdown_tokens_present": all(not row["missing_tokens"] for row in results.values()),
    "audit_warning_categories": warning_categories,
    "reaction_identity_warnings": sum(1 for row in audit["warnings"] for w in row["warnings"] if w.startswith("reaction_identity:")),
    "outcome_target_warnings": sum(1 for row in audit["warnings"] for w in row["warnings"] if w.startswith("outcome_target:")),
    "source_reviewed_records": audit["source_reviewed_records"],
    "chemical_accuracy_verified": False,
}
rt.write(RUN / "raw_model/verification_summary.json", summary)
print(json.dumps(summary, ensure_ascii=False, indent=2))
