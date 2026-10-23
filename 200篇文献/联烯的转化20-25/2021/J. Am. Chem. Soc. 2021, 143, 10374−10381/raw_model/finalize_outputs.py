#!/usr/bin/env python3
"""Persist read-only audit and literal source-mention review artifacts."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools
import standard_records

run = ROOT / "runs/current"
data = json.loads((run / "dataset.json").read_text())
audit = reaction_tools.dataset_audit(data)
(run / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

terms = ["2a", "3a", "4a", "10u", "14i", "carfilzomib", "ACP (65"]
mentions = standard_records.source_mentions(data, terms, context=160)
(run / "raw_model/source_mentions.json").write_text(
    json.dumps(mentions, ensure_ascii=False, indent=2) + "\n"
)

issues = Counter(w.split(":", 1)[0] for row in audit["warnings"] for w in row["warnings"])
summary = json.loads((run / "summary.json").read_text())
summary["v3_molecules"] = len(data["standard"]["molecules"])
summary["v3_reactions"] = len(data["standard"]["reactions"])
summary["audit_warning_records"] = len(audit["warnings"])
summary["audit_issue_counts"] = dict(issues)
summary["semantic_coverage_complete"] = False
(run / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"mentions": len(mentions["matches"]), "warning_records": len(audit["warnings"]),
                  "issues": dict(issues)}, ensure_ascii=False))
