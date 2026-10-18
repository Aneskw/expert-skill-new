#!/usr/bin/env python3
"""Persist deterministic audit and mark the five prepared page tasks attempted."""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))

from reaction_tools import dataset_audit, read, write


data = read(RUN / "dataset.json")
write(RUN / "audit.json", dataset_audit(data))
tasks = read(RUN / "tasks.json")
for task in tasks:
    task["state"] = "attempted"
    task["reason"] = "Original page image and extracted text were inspected; detailed region findings and record links are in coverage/pages_0001_0005.json."
write(RUN / "tasks.json", tasks)
print(json.dumps({"audit_warnings": len(read(RUN / "audit.json")["warnings"]), "tasks_attempted": len(tasks)}, ensure_ascii=False))
