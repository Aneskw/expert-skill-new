#!/usr/bin/env python3
"""Remove accidental duplicate source-context IDs before the final merge.

The standard merge intentionally preserves older source context.  An earlier
intermediate run changed the Figure-1 context in place, causing the merge to
retain both versions under one context_id.  This mechanical repair keeps the
first pre-existing item for each ID; the page-9 X-ray occurrence is imported
under its own new ID by build_pages_0007_0009.py.
"""
import json
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "dataset.json"
data = json.loads(path.read_text())
contexts = data.get("standard", {}).get("source_context", [])
seen = set()
unique = []
for item in contexts:
    key = item.get("context_id")
    if key in seen:
        continue
    seen.add(key)
    unique.append(item)
data["standard"]["source_context"] = unique
tmp = path.with_name(path.name + ".tmp")
tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
tmp.replace(path)
print(json.dumps({"before": len(contexts), "after": len(unique)}, ensure_ascii=False))
