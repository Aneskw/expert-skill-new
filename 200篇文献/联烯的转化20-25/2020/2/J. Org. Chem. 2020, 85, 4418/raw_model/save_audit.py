#!/usr/bin/env python3
"""Save the deterministic advisory audit without modifying source records."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit

dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
result = dataset_audit(dataset)
(ROOT / "runs/current/raw_model/audit.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + "\n"
)
print(json.dumps({"warnings": len(result["warnings"]),
                  "source_reviewed_records": result["source_reviewed_records"]}, ensure_ascii=False))
