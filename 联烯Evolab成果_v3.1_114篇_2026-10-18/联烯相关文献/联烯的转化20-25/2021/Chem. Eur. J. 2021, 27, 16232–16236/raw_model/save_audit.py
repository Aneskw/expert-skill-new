#!/usr/bin/env python3
"""Persist the read-only advisory audit beside the merged dataset."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit

dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
result = dataset_audit(dataset)
(ROOT / "runs/current/audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"warning_records": len(result["warnings"]), "source_reviewed_records": result["source_reviewed_records"]}, ensure_ascii=False))
