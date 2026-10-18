#!/usr/bin/env python3
"""Remove one erroneous, unused pre-merge candidate before the final standard merge."""

import json
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "dataset.json"
data = json.loads(path.read_text())
target = "doc_fdf0107b6851842e:mol:na2so4"
standard = data.get("standard")
if standard:
    referenced = {p["molecule_id"] for r in standard["reactions"] for p in r["participants"]}
    if target in referenced:
        raise RuntimeError("refuse to remove a referenced molecule")
    before = len(standard["molecules"])
    standard["molecules"] = [m for m in standard["molecules"] if m["molecule_id"] != target]
    if len(standard["molecules"]) not in (before, before - 1):
        raise RuntimeError("unexpected removal count")
tmp = path.with_name(path.name + ".tmp")
tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False))
tmp.replace(path)
