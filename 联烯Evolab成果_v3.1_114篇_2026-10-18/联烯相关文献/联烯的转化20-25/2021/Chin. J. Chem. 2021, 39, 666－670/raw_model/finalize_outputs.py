#!/usr/bin/env python3
"""Persist deterministic audit output without modifying extraction tools."""
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
(run / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
issues = Counter()
for row in audit["warnings"]:
    issues.update(w.split(":", 1)[0] for w in row["warnings"])
summary = {
    "warning_records": len(audit["warnings"]),
    "warning_incidences": dict(issues),
    "reaction_identity": issues.get("reaction_identity", 0),
    "outcome_target": issues.get("outcome_target", 0),
    "measurement_basis": issues.get("measurement_basis", 0),
    "structure_gap": issues.get("structure_gap", 0),
    "chemical_accuracy_verified": False,
}
(run / "raw_model/audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(summary, ensure_ascii=False))
