#!/usr/bin/env python3
"""Persist the deterministic post-merge audit without changing its rules."""

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

run = ROOT / "runs/current"
audit = rt.dataset_audit(json.loads((run / "dataset.json").read_text()))
st.write(run / "raw_model/pages_0007_0012_audit.json", audit)
st.write(run / "quality_checks.json", audit)
print(json.dumps({"warnings": len(audit["warnings"]), "record_kind_counts": audit["record_kind_counts"]}, ensure_ascii=False))
