#!/usr/bin/env python3
"""Persist the read-only advisory audit for this bounded extraction."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools  # noqa: E402

data = reaction_tools.read(ROOT / "runs/current/dataset.json")
audit = reaction_tools.dataset_audit(data)
reaction_tools.write(ROOT / "runs/current/quality_checks.json", audit)
reaction_tools.write(ROOT / "runs/current/raw_model/pages_0001_0006_audit.json", audit)
print({"warning_records": len(audit["warnings"]), "source_reviewed_records": audit["source_reviewed_records"]})
