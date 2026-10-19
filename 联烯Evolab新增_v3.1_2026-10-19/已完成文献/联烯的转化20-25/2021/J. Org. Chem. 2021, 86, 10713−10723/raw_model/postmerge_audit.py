#!/usr/bin/env python3
"""Persist the unmodified advisory audit for the current merged dataset."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

from reaction_tools import dataset_audit
from standard_records import write

run = ROOT / "runs/current"
dataset = json.loads((run / "dataset.json").read_text())
result = dataset_audit(dataset)
write(run / "quality_checks.json", result)
write(run / "raw_model/audit.json", result)
print(json.dumps({
    "records": sum(result["record_kind_counts"].values()),
    "warning_records": len(result["warnings"]),
    "source_reviewed_records": result["source_reviewed_records"],
}, ensure_ascii=False))
