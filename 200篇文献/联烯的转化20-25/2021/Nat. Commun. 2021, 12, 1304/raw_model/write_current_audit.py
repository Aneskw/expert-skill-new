#!/usr/bin/env python3
"""Persist the current advisory audit without changing audit rules."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
from reaction_tools import dataset_audit

dataset_path = ROOT / "runs/current/dataset.json"
output_path = ROOT / "runs/current/quality_checks.json"
dataset = json.loads(dataset_path.read_text(encoding="utf-8"))
result = dataset_audit(dataset)
output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
