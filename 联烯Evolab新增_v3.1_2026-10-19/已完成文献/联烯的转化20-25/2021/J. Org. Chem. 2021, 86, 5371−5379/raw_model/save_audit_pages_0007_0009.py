#!/usr/bin/env python3
"""Persist the advisory audit output after the pages 7--9 merge."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records as st

dataset = json.loads((ROOT / "runs/current/dataset.json").read_text())
result = reaction_tools.dataset_audit(dataset)
st.write(ROOT / "runs/current/audit_pages_0007_0009/audit.json", result)
print(json.dumps({
    "audit": str(ROOT / "runs/current/audit_pages_0007_0009/audit.json"),
    "warnings": len(result["warnings"]),
    "chemical_accuracy_verified": result["chemical_accuracy_verified"],
}, ensure_ascii=False, indent=2))
