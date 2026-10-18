#!/usr/bin/env python3
"""Persist final validation/audit evidence without altering extraction rules."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

run = ROOT / "runs/current"
dataset = st.read(run / "dataset.json")
standard = dataset["standard"]
standard_errors = st.validate(standard, run)
legacy_errors = rt.validate(dataset, run, strict_record_ids={r["record_id"] for r in dataset["experiments"]})
audit = rt.dataset_audit(dataset)
st.write(run / "audit.json", audit)

status_counts = Counter(r["review"]["status"] for r in standard["reactions"])
completeness_counts = Counter(r["completeness"] for r in standard["reactions"])
warning_prefixes = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        warning_prefixes[warning.split(":", 1)[0]] += 1

checks = {
    "standard_validate": {"valid": not standard_errors, "errors": standard_errors},
    "strict_dataset_validate": {"valid": not legacy_errors, "errors": legacy_errors},
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
    "review_status_counts": dict(status_counts),
    "completeness_counts": dict(completeness_counts),
    "audit_warning_records": len(audit["warnings"]),
    "audit_warning_prefix_counts": dict(warning_prefixes),
    "checks": [],
}

targets = {
    "records/doc_08515448f9baa759__scheme2_acetone_rt_10min_2equiv.md": ["2.0 equiv", "acetone, rt, 10 min"],
    "records/doc_08515448f9baa759__scheme7_2f_to_7fa.md": ["5 mol%", "1,4-dioxane:H2O (2:1)", "60 °C", "2 h"],
    "records/doc_08515448f9baa759__scheme7_2h_to_11h.md": ["H2 (1 atm)", "72 h"],
    "records/doc_08515448f9baa759__scheme7_15b_to_16b.md": ["IPrAuCl", "AgOTf", "5 mol%", "80 °C", "24 h"],
    "records/doc_08515448f9baa759__scheme5_2a_to_7a.md": ["allocation between rt annulation and 0 °C base treatment is not stated"],
}
for rel, tokens in targets.items():
    text = (run / rel).read_text()
    missing = [token for token in tokens if token not in text]
    checks["checks"].append({"file": rel, "required_tokens": tokens, "missing": missing, "passed": not missing})

coverage = st.read(run / "coverage/pages_0001_0005.json")
checks["coverage"] = {
    "pages": [p["page"] for p in coverage["pages"]],
    "all_attempted": all(p["state"] == "attempted" for p in coverage["pages"]),
    "all_regions_nonempty": all(bool(p["regions"]) for p in coverage["pages"]),
}
checks["all_checks_passed"] = (
    not standard_errors and not legacy_errors and
    status_counts == {"needs_review": len(standard["reactions"])} and
    completeness_counts == {"partial": len(standard["reactions"])} and
    all(row["passed"] for row in checks["checks"]) and
    checks["coverage"] == {"pages": [1, 2, 3, 4, 5], "all_attempted": True, "all_regions_nonempty": True}
)
st.write(run / "raw_model/verification_summary.json", checks)
print(json.dumps(checks, ensure_ascii=False, indent=2))
if not checks["all_checks_passed"]:
    raise SystemExit(1)
