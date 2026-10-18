#!/usr/bin/env python3
"""Post-merge assertions and audit capture for the current extraction only."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

RUN = ROOT / "runs/current"
data = json.loads((RUN / "dataset.json").read_text())
standard = data["standard"]
coverage = json.loads((RUN / "coverage/pages_0001_0004.json").read_text())
reaction_ids = {r["reaction_id"] for r in standard["reactions"]}
covered = {rid for page in coverage["pages"] for rid in page["record_ids"]}

assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert covered <= reaction_ids
assert reaction_ids <= covered
assert all(r["review"]["status"] == "needs_review" for r in standard["reactions"])
assert all(r["completeness"] == "partial" for r in standard["reactions"])
assert all(all(m["stereochemistry"][k]["normalized_value"] is not None
               for k in ("chirality_status", "geometry_status", "representation_status"))
           for m in standard["molecules"])
assert len(list((RUN / "records").glob("*.json"))) == len(standard["reactions"])

audit = rt.dataset_audit(data)
st.write(RUN / "audit.json", audit)
warning_classes = sorted({w.split(":", 1)[0]
                          for row in audit["warnings"] for w in row["warnings"]})
assert warning_classes == ["measurement_basis"]
assert not any(any(w.startswith(("reaction_identity:", "outcome_target:", "structure_gap:", "missing_steps:"))
                   for w in row["warnings"]) for row in audit["warnings"])

checks = {
    "scope_3aa": {
        "path": "records/doc_739c1c0b69cb742a__scheme2_3aa.md",
        "tokens": ["0.2 mmol", "0.24 mmol", "5 mol%", "1 equiv", "0.2 M", "100 °C", "4 h", "73%", "isolated yield", "E only"],
    },
    "mechanistic_3aa_D2O": {
        "path": "records/doc_739c1c0b69cb742a__scheme3a_3aa_D2O_annulation.md",
        "tokens": ["5 mol%", "1 equiv", "DCE/D2O", "100 °C", "0.5 h", "58%", "30", "81"],
    },
    "derivatization_11": {
        "path": "records/doc_739c1c0b69cb742a__eq2_3aa_to_11.md",
        "tokens": ["LiAlH4 (3 equiv)", "MeOH (4 equiv)", "THF", "rt", "3 h", "61%"],
    },
    "derivatization_12": {
        "path": "records/doc_739c1c0b69cb742a__eq2_3aa_to_12.md",
        "tokens": ["RMgBr (2 equiv)", "THF", "rt", "12 h", "under N2", "51%"],
    },
}
for item in checks.values():
    text = (RUN / item["path"]).read_text()
    item["missing_tokens"] = [token for token in item["tokens"] if token not in text]
    assert not item["missing_tokens"], (item["path"], item["missing_tokens"])

result = {
    "document_id": data["documents"][0]["document_id"],
    "source_sha256": data["documents"][0]["sha256"],
    "reactions": len(standard["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
    "molecules": len(standard["molecules"]),
    "source_context": len(standard["source_context"]),
    "coverage_pages": [p["page"] for p in coverage["pages"]],
    "warning_record_count": len(audit["warnings"]),
    "warning_classes": warning_classes,
    "source_reviewed_records": audit["source_reviewed_records"],
    "markdown_checks": checks,
    "chemical_accuracy_verified": False,
    "semantic_accuracy_verified": False,
}
st.write(RUN / "raw_model/verification_summary.json", result)
print(json.dumps(result, ensure_ascii=False, indent=2))
