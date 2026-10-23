#!/usr/bin/env python3
"""Persist page-7 coverage, strict-validation, and audit checks."""
from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools as rt  # noqa: E402
import standard_records as st  # noqa: E402


base = ROOT / "runs/current"
data = st.read(base / "dataset.json")
standard = st.read(base / "standard_records.json")
coverage = st.read(base / "coverage/pages_0007_0007.json")

standard_errors = st.validate(standard, base)
strict_ids = {record["record_id"] for record in data["experiments"]}
legacy_errors = rt.validate(data, base, strict_record_ids=strict_ids)

coverage_errors = []
if coverage.get("document_id") != "doc_3f18b691f7dce01d":
    coverage_errors.append("coverage document_id mismatch")
pages = coverage.get("pages", [])
if [page.get("page") for page in pages] != [7]:
    coverage_errors.append("coverage must contain exactly page 7")
for page in pages:
    if page.get("state") != "attempted":
        coverage_errors.append("page 7 state must be attempted")
    if not page.get("regions"):
        coverage_errors.append("page 7 regions must be nonempty")
    if page.get("record_ids"):
        coverage_errors.append("page 7 bibliography must not claim reaction records")
    for region in page.get("regions", []):
        if not region.get("evidence_ids"):
            coverage_errors.append("page 7 region lacks evidence_ids")

context_id = "doc_3f18b691f7dce01d:context:page7_references"
contexts = {item.get("context_id"): item for item in standard.get("source_context", [])}
if context_id not in contexts:
    coverage_errors.append("page 7 references source_context missing after merge")
elif set(contexts[context_id].get("evidence_ids", [])) != {
    "doc_3f18b691f7dce01d:page7:image",
    "doc_3f18b691f7dce01d:page7:text",
}:
    coverage_errors.append("page 7 references source_context evidence mismatch")

if len(standard.get("reactions", [])) != 63:
    coverage_errors.append("cumulative standard reaction count changed unexpectedly")
if len(standard.get("molecules", [])) != 147:
    coverage_errors.append("cumulative molecule count changed unexpectedly")
if any(record.get("review", {}).get("status") != "needs_review" for record in standard.get("reactions", [])):
    coverage_errors.append("one or more machine reactions are not needs_review")

audit = rt.dataset_audit(data)
validation = {
    "page_range": [7, 7],
    "standard": {"valid": not standard_errors, "errors": standard_errors},
    "legacy_strict": {"valid": not legacy_errors, "errors": legacy_errors},
    "coverage": {"valid": not coverage_errors, "errors": coverage_errors},
    "page7_new_molecules": 0,
    "page7_new_reactions": 0,
    "page7_new_source_context": 1,
    "cumulative_molecules": len(standard.get("molecules", [])),
    "cumulative_reactions": len(standard.get("reactions", [])),
    "audit_warning_records": len(audit.get("warnings", [])),
    "audit_warning_count": sum(len(item.get("warnings", [])) for item in audit.get("warnings", [])),
    "chemical_accuracy_verified": False,
    "semantic_accuracy_verified": False,
}

st.write(base / "raw_model/strict_validation_pages_0007_0007.json", validation)
st.write(base / "raw_model/audit.json", audit)
st.write(base / "quality_checks.json", audit)
print(json.dumps(validation, ensure_ascii=False, indent=2))

if standard_errors or legacy_errors or coverage_errors:
    raise SystemExit(1)
