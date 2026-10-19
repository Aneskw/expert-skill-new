#!/usr/bin/env python3
"""Persist the post-merge audit for the bounded page 7--9 extraction."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
from reaction_tools import dataset_audit  # noqa: E402


dataset = json.loads((ROOT / "runs/current/dataset.json").read_text(encoding="utf-8"))
audit = dataset_audit(dataset)
payload = json.dumps(audit, ensure_ascii=False, indent=2) + "\n"
(ROOT / "runs/current/quality_checks.json").write_text(payload, encoding="utf-8")
(ROOT / "runs/current/raw_model/pages_0007_0009_audit.json").write_text(
    payload, encoding="utf-8"
)

print(
    json.dumps(
        {
            "warning_records": len(audit.get("warnings", [])),
            "source_reviewed_records": audit.get("source_reviewed_records"),
            "chemical_accuracy_verified": audit.get("chemical_accuracy_verified"),
        },
        ensure_ascii=False,
    )
)
