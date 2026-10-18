#!/usr/bin/env python3
"""Regenerate the final audit artifacts from the merged v3.1 dataset."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
SKILL_SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SKILL_SCRIPTS))

from reaction_tools import dataset_audit  # noqa: E402


dataset = json.loads((RUN / "dataset.json").read_text())
audit = dataset_audit(dataset)
encoded = json.dumps(audit, ensure_ascii=False, indent=2) + "\n"

(RUN / "quality_checks.json").write_text(encoded)
(RUN / "raw_model" / "pages_0007_0012_audit.json").write_text(encoded)

print(
    json.dumps(
        {
            "record_kind_counts": audit["record_kind_counts"],
            "warning_rows": len(audit["warnings"]),
            "source_reviewed_records": audit["source_reviewed_records"],
            "chemical_accuracy_verified": audit["chemical_accuracy_verified"],
        },
        ensure_ascii=False,
    )
)
