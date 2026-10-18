#!/usr/bin/env python3
"""Final source-specific integrity checks and audit serialization."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools

base = ROOT / "runs/current"
data = json.loads((base / "dataset.json").read_text())
audit = reaction_tools.dataset_audit(data)
(base / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

standard = json.loads((base / "standard_records.json").read_text())
coverage = json.loads((base / "coverage/pages_0001_0006.json").read_text())
assert data["standard"] == standard
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5, 6]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert len(standard["reactions"]) == 92 and len(standard["molecules"]) == 150
assert all(r["review"]["status"] == "needs_review" for r in standard["reactions"])
assert data["documents"][0]["title"].startswith("Diastereodivergent Synthesis")
assert data["documents"][0]["supplementary_material"]["status"] == "unknown"
assert all(st["normalized_value"] is not None for m in standard["molecules"] for key, st in m["stereochemistry"].items()
           if key in ("chirality_status", "geometry_status", "representation_status"))

checks = {
    "table1_entry10": {
        "path": base / "records/doc_1ba437add780cf5e__table1_entry10.md",
        "tokens": ["0.1", "0.2", "4", "6", "5", "200", "0.5", "30", "48", "14:1", "70%", "71%", "conflict"],
    },
    "table2_3ja": {
        "path": base / "records/doc_1ba437add780cf5e__table2_3ja.md",
        "tokens": ["83%", "56%", "conflict", "isolated", "4:1", ">99"],
    },
    "scheme5_5d": {
        "path": base / "records/doc_1ba437add780cf5e__scheme5_5d.md",
        "tokens": ["10", "1 atm balloon", "5", "78", ">20:1", "Pd/C", "H2"],
    },
}
verified = {}
for name, spec in checks.items():
    text = spec["path"].read_text()
    missing = [t for t in spec["tokens"] if t not in text]
    assert not missing, (name, missing)
    verified[name] = {"path": str(spec["path"].relative_to(ROOT)), "missing_tokens": missing}

issue_counts = Counter()
for row in audit["warnings"]:
    issue_counts.update({w.split(":", 1)[0] for w in row["warnings"]})

summary = {
    "documents": len(standard["documents"]),
    "pages_attempted": [1, 2, 3, 4, 5, 6],
    "reactions": len(standard["reactions"]),
    "molecules": len(standard["molecules"]),
    "source_context": len(standard["source_context"]),
    "source_reviewed_records": sum(r["review"]["status"] == "source_reviewed" for r in standard["reactions"]),
    "warning_records": len(audit["warnings"]),
    "warning_categories": dict(issue_counts),
    "reaction_identity_warnings": issue_counts.get("reaction_identity", 0),
    "outcome_target_warnings": issue_counts.get("outcome_target", 0),
    "markdown_checks": verified,
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
(base / "raw_model/verification_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(summary, ensure_ascii=False, indent=2))
