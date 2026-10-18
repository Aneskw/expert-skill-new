#!/usr/bin/env python3
"""Regenerate final audit artifacts after merging pages 13-18."""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from reaction_tools import dataset_audit  # noqa: E402


dataset = json.loads((RUN / "dataset.json").read_text())
audit = dataset_audit(dataset)
encoded = json.dumps(audit, ensure_ascii=False, indent=2) + "\n"
(RUN / "quality_checks.json").write_text(encoded)
(RUN / "raw_model" / "pages_0013_0018_audit.json").write_text(encoded)

batch_prefixes = tuple(f"doc_832eb8cea3728513:rxn:rxn_p{page}_" for page in range(13, 19))
batch_warnings = [w for w in audit["warnings"] if w["record_id"].startswith(batch_prefixes)]
print(json.dumps({
    "record_kind_counts": audit["record_kind_counts"],
    "full_warning_rows": len(audit["warnings"]),
    "batch_warning_rows": len(batch_warnings),
    "batch_warning_types": sorted({warning.split(":", 1)[0] for row in batch_warnings for warning in row["warnings"]}),
    "chemical_accuracy_verified": audit["chemical_accuracy_verified"],
}, ensure_ascii=False))
