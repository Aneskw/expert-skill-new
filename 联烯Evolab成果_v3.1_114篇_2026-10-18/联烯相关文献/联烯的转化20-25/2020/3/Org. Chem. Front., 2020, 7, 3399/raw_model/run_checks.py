#!/usr/bin/env python3
"""Persist validation and advisory audit snapshots for this extraction run."""
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
standard_errors = st.validate(standard, base)
strict_ids = {r["record_id"] for r in data["experiments"]}
legacy_errors = rt.validate(data, base, strict_record_ids=strict_ids)
validation = {
    "standard": {"valid": not standard_errors, "errors": standard_errors},
    "legacy_strict": {"valid": not legacy_errors, "errors": legacy_errors},
    "chemical_accuracy_verified": False,
}
audit = rt.dataset_audit(data)
coverage = st.read(base / "coverage/pages_0001_0006.json")
coverage_errors = []
if coverage.get("document_id") != "doc_3f18b691f7dce01d":
    coverage_errors.append("coverage document_id mismatch")
if [p.get("page") for p in coverage.get("pages", [])] != [1, 2, 3, 4, 5, 6]:
    coverage_errors.append("coverage must contain exactly pages 1..6")
known_ids = {r["reaction_id"] for r in standard["reactions"]}
for page in coverage.get("pages", []):
    if page.get("state") != "attempted" or not page.get("regions"):
        coverage_errors.append(f"page {page.get('page')} lacks attempted/nonempty regions")
    missing = set(page.get("record_ids", [])) - known_ids
    if missing:
        coverage_errors.append(f"page {page.get('page')} has unknown record_ids: {sorted(missing)}")
preview = (base / "规范记录预览.md").read_text() + "\n" + "\n".join(
    p.read_text() for p in sorted((base / "records").glob("*.md"))
)
required_markdown = [
    "0.6 mmol of 1a was used in two portions", "CuBr·SMe2", "HMPA", "TMSCl",
    "DIBAL-H", "Phthalimide", "DEAD", "2.2 equiv", "50 mol%",
    "57% yield over two steps", "inconsistent with the separately depicted 48% alpha and 32% gamma",
]
missing_markdown = [term for term in required_markdown if term not in preview]
if coverage_errors or missing_markdown:
    validation["coverage_errors"] = coverage_errors
    validation["missing_markdown_terms"] = missing_markdown
st.write(base / "raw_model/strict_validation.json", validation)
st.write(base / "raw_model/audit.json", audit)
print(json.dumps({
    "validation": validation,
    "audit_warning_records": len(audit["warnings"]),
    "audit_warning_count": sum(len(x["warnings"]) for x in audit["warnings"]),
}, ensure_ascii=False, indent=2))
if standard_errors or legacy_errors or coverage_errors or missing_markdown:
    raise SystemExit(1)
