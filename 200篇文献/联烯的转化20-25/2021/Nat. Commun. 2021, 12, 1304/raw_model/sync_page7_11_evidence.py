#!/usr/bin/env python3
"""Copy prepared page 7-11 image evidence into the evolving submission document."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[3]
run = ROOT / "runs/current"
submission_path = run / "submitted.json"
template_path = run / "raw_model/blank_template_pages7_11.json"

submission = json.loads(submission_path.read_text(encoding="utf-8"))
template = json.loads(template_path.read_text(encoding="utf-8"))
source_doc = {d["document_id"]: d for d in template["documents"]}
for document in submission["documents"]:
    prepared = source_doc[document["document_id"]]
    existing = {item["evidence_id"] for item in document.get("evidence", [])}
    for item in prepared.get("evidence", []):
        if item["evidence_id"] not in existing:
            document.setdefault("evidence", []).append(item)
            existing.add(item["evidence_id"])

submission_path.write_text(json.dumps(submission, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
