"""Apply source-visible article metadata before standard initialization."""
import json
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "dataset.json"
data = json.loads(path.read_text())
doc = next(d for d in data["documents"] if d["document_id"] == "doc_1cf04303b61f633d")
doc["title"] = "Copper Catalyzed Protosilylation/Protoborylation of gem-Difluoroallenes"
doc["supplementary_material"] = {
    "status": "mentioned_not_supplied",
    "note": "Page 1 explicitly states that Supporting Information is available, but no SI file was supplied or explicitly associated in this session.",
}
tmp = path.with_suffix(".json.tmp")
tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")
tmp.replace(path)
