#!/usr/bin/env python3
"""Persist the read-only advisory audit beside other raw machine artifacts."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools

dataset = reaction_tools.read(ROOT / "runs/current/dataset.json")
result = reaction_tools.dataset_audit(dataset)
(ROOT / "runs/current/raw_model/audit.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + "\n"
)
print(json.dumps({"warning_records": len(result["warnings"]),
                  "source_reviewed_records": result["source_reviewed_records"]},
                 ensure_ascii=False))
