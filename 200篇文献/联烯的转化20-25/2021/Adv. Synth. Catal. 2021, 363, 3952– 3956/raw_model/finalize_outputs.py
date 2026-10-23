#!/usr/bin/env python3
"""Persist audit, source-mention hints, and generated-Markdown checks."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records as st

run = ROOT / "runs/current"
data = json.loads((run / "dataset.json").read_text())
(run / "audit.json").write_text(
    json.dumps(reaction_tools.dataset_audit(data), ensure_ascii=False, indent=2) + "\n"
)
(run / "raw_model/source_mentions.json").write_text(
    json.dumps(st.source_mentions(data, ["1a", "2a", "3a", "1b", "3b", "4a"]), ensure_ascii=False, indent=2) + "\n"
)

checks = {
    "doc_e6c9ddf245b60941__scheme_scope_3b.md": [
        "34.8 mg", "0.18 mmol", "59.1 mg", "0.36 mmol", "3.3 mg",
        "0.018 mmol", "61.2 mg", "2 mL", "3 × 2 mL", "7:1",
        "37 mg", "62%", "until starting material disappeared by TLC",
    ],
    "doc_e6c9ddf245b60941__table1_entry2.md": [
        "0.15 mmol", "0.30 mmol", "200 mol%", "10 mol%", "1 mL",
        "68%", "cross-page duplicate merged",
    ],
}
result = {}
for name, tokens in checks.items():
    text = (run / "records" / name).read_text()
    result[name] = {"required_tokens": tokens, "missing_tokens": [t for t in tokens if t not in text]}
    if result[name]["missing_tokens"]:
        raise SystemExit(f"Generated Markdown missing tokens: {result[name]}")
(run / "raw_model/verification_summary.json").write_text(
    json.dumps(result, ensure_ascii=False, indent=2) + "\n"
)

coverage = json.loads((run / "coverage/pages_0001_0005.json").read_text())
assert coverage["document_id"] == "doc_e6c9ddf245b60941"
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
