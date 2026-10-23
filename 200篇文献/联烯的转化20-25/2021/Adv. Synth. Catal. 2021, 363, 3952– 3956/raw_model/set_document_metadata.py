#!/usr/bin/env python3
"""Set source-derived bibliographic metadata before standard init."""
import json
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "dataset.json"
data = json.loads(path.read_text())
doc = data["documents"][0]
assert doc["document_id"] == "doc_e6c9ddf245b60941"
doc["title"] = (
    "Metal-Catalyzed Reactivity Reversal in the Sulfonylation Reaction of "
    "α-Allenols: Controlled Synthesis of 4-(Arylsulfonyl)-2,5-Dihydrofurans"
)
tmp = path.with_suffix(".json.tmp")
tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
tmp.replace(path)
