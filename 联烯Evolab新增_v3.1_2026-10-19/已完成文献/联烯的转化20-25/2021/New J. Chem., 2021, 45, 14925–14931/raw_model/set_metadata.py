#!/usr/bin/env python3
import json
from pathlib import Path

path = Path("runs/current/dataset.json")
data = json.loads(path.read_text())
doc = next(d for d in data["documents"] if d["document_id"] == "doc_0414a2e90dea4adf")
doc["title"] = "Copper(II)-catalyzed protoboration of allenes in aqueous media and open air"
doc["doi"] = "10.1039/d0nj02010a"
tmp = path.with_name(path.name + ".tmp")
tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
tmp.replace(path)
