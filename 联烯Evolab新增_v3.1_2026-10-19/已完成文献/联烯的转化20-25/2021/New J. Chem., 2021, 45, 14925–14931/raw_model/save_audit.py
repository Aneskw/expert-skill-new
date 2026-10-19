#!/usr/bin/env python3
import json
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(root / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit

dataset = json.loads((root / "runs/current/dataset.json").read_text())
result = dataset_audit(dataset)
target = root / "runs/current/raw_model/audit.json"
target.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
