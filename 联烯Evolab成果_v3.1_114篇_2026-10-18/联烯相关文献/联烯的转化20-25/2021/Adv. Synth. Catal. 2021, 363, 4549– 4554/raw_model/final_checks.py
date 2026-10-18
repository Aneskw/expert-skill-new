#!/usr/bin/env python3
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
run = ROOT / "runs/current"
data = json.loads((run / "dataset.json").read_text())
standard = json.loads((run / "standard_records.json").read_text())
coverage = json.loads((run / "coverage/pages_0001_0006.json").read_text())
audit = json.loads((run / "audit.json").read_text())

assert standard == data["standard"]
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5, 6]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert len({r["reaction_id"] for r in standard["reactions"]}) == len(standard["reactions"])
assert len({m["molecule_id"] for m in standard["molecules"]}) == len(standard["molecules"])
assert all(r["review"]["status"] == "needs_review" for r in standard["reactions"])
assert len({c["context_id"] for c in standard["source_context"]}) == len(standard["source_context"])

checks = {
    "table2_4a": {
        "path": "records/doc_c690e758ccda5a7c__table2_4a.md",
        "tokens": ["0.2 mmol", "1.0 equiv", "10 mol%", "2.0 mL", "12 h", "100 mol%",
                   "NaHCO3", "EtOAc", "hexanes", "DCM", "silica gel", "TLC"],
    },
    "scheme4_4a_to_5a": {
        "path": "records/doc_c690e758ccda5a7c__scheme4_a_hydrogenation.md",
        "tokens": ["50 mol%", "H2 (balloon)", "IPA", "95%"],
    },
    "scheme2_3b_recovery": {
        "path": "records/doc_c690e758ccda5a7c__scheme2_4b_cyclization.md",
        "tokens": ["almost quantitative", "complex mixture", "no desired product or desulfonation product"],
    },
}
for item in checks.values():
    text = (run / item["path"]).read_text()
    item["missing_tokens"] = [t for t in item["tokens"] if t not in text]
    assert not item["missing_tokens"], (item["path"], item["missing_tokens"])

warning_categories = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        warning_categories[warning.split(":", 1)[0]] += 1

pdf = run / data["documents"][0]["source_pdf"]
assert hashlib.sha256(pdf.read_bytes()).hexdigest() == data["documents"][0]["sha256"]

summary = {
    "document_id": data["documents"][0]["document_id"],
    "pages": [1, 2, 3, 4, 5, 6],
    "reactions": len(standard["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
    "molecules": len(standard["molecules"]),
    "source_context": len(standard["source_context"]),
    "reaction_types": dict(Counter(r["record_type"] for r in standard["reactions"])),
    "all_needs_review": True,
    "markdown_checks": checks,
    "audit_warning_records": len(audit["warnings"]),
    "audit_warning_categories": dict(warning_categories),
    "source_hash_verified": True,
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
(run / "raw_model/verification_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(summary, ensure_ascii=False, indent=2))
