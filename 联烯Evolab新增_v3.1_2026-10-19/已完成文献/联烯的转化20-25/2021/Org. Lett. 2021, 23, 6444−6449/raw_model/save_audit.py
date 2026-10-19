#!/usr/bin/env python3
"""Save the final advisory audit without modifying source-derived facts."""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit

run = ROOT / "runs/current"
data = json.loads((run / "dataset.json").read_text())
audit = dataset_audit(data)
issues = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        issues[warning.split(":", 1)[0]] += 1
audit["warning_summary"] = {
    "records": len(audit["warnings"]),
    "items": sum(issues.values()),
    "by_issue": dict(issues),
}
(run / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(audit["warning_summary"], ensure_ascii=False))
