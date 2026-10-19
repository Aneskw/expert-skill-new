#!/usr/bin/env python3
"""Save deterministic audit output and mark the bounded page task queue attempted."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt

run = ROOT / "runs/current"
dataset = rt.read(run / "dataset.json")
rt.write(run / "raw_model/audit.json", rt.dataset_audit(dataset))
tasks = rt.read(run / "tasks.json")
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original full-page image and text inspected; exact region-to-record decisions are in coverage/pages_0001_0004.json."
rt.write(run / "tasks.json", tasks)
