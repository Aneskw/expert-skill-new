#!/usr/bin/env python3
"""Reset only this session's machine import before a clean final merge.

The prepared dataset began with zero experiments and no v3 standard bundle.
This removes iterative imports created by build_submission.py; documents,
evidence, candidates, source assets, and prior run-log entries are preserved.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools

path = ROOT / "runs/current/dataset.json"
data = reaction_tools.read(path)
data["experiments"] = []
data.pop("standard", None)
data["run_log"] = [row for row in data.get("run_log", []) if row.get("event") != "records_merged"]
reaction_tools.write(path, data)
print("Cleared iterative session imports; source/evidence retained.")
