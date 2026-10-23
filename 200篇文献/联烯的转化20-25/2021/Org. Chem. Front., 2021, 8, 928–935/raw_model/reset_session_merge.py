#!/usr/bin/env python3
"""Remove only this session's first machine merge before corrected deduplication."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st  # noqa: E402

path = ROOT / "runs/current/dataset.json"
data = json.loads(path.read_text())
doc = "doc_755164270d1e5816"
expected_duplicate = f"{doc}:scheme2_3a"
records = data.get("experiments", [])
standard = data.get("standard")
if len(records) != 67 or not standard or len(standard.get("reactions", [])) != 67:
    raise SystemExit("Refuse reset: dataset is not the 67-record machine merge created by this session")
if not all(r.get("record_id", "").startswith(doc + ":") and r.get("status") == "needs_review" for r in records):
    raise SystemExit("Refuse reset: unexpected/non-machine experiment present")
if expected_duplicate not in {r["record_id"] for r in records}:
    raise SystemExit("Refuse reset: expected duplicate is absent")
data["experiments"] = []
data.pop("standard")
data.setdefault("run_log", []).append({"event": "session_machine_merge_reset_for_exact_duplicate_removal",
                                         "record_id": expected_duplicate})
st.write(path, data)
print("Reset only the 67 machine records from this session; prepared documents and candidates preserved.")
