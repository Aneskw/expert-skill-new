#!/usr/bin/env python3
"""Remove the duplicate source-context object created by the first page-7 merge."""
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[3]
dataset_path = ROOT / "runs/current/dataset.json"
data = json.loads(dataset_path.read_text())

xray_id = "doc_cb03cf42d5a5f24b:context_xray_3ba"
contexts = data["standard"]["source_context"]
original = [
    context for context in contexts
    if context.get("context_id") == xray_id
    and context.get("evidence_ids") == [
        "doc_cb03cf42d5a5f24b:page2:image",
        "doc_cb03cf42d5a5f24b:page4:image",
    ]
]
assert len(original) == 1

cleaned = []
for context in contexts:
    if context.get("context_id") == xray_id and context is not original[0]:
        continue
    cleaned.append(context)

assert len(contexts) == 6
assert len(cleaned) == 5
data["standard"]["source_context"] = cleaned
dataset_path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
print({"source_context_before": len(contexts), "source_context_after": len(cleaned)})
