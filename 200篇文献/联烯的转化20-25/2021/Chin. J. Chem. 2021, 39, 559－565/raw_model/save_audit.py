#!/usr/bin/env python3
"""Persist the deterministic post-merge audit without modifying records."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools

run = ROOT / "runs/current"
data = json.loads((run / "dataset.json").read_text())
audit = reaction_tools.dataset_audit(data)
(run / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
(run / "raw_model/audit_stdout.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"warning_records": len(audit["warnings"]), "source_reviewed_records": audit["source_reviewed_records"]}, ensure_ascii=False))
