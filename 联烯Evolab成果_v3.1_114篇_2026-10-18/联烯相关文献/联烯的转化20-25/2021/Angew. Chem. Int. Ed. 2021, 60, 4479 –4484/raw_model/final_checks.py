#!/usr/bin/env python3
import json
from pathlib import Path

root = Path(__file__).resolve().parents[3]
run = root / "runs/current"
dataset = json.loads((run / "dataset.json").read_text())
standard = json.loads((run / "standard_records.json").read_text())
coverage = json.loads((run / "coverage/pages_0001_0006.json").read_text())

assert dataset["standard"] == standard
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5, 6]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert len(standard["reactions"]) == 91
assert len(standard["molecules"]) == 158
assert sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]) == 87
assert all(r["review"]["status"] == "needs_review" for r in standard["reactions"])
assert all(all(m["stereochemistry"][k]["normalized_value"] is not None
               for k in ("chirality_status", "geometry_status", "representation_status"))
           for m in standard["molecules"])

checks = {
    "table1_entry4": {
        "file": f"{dataset['documents'][0]['document_id']}__table1_entry4.md",
        "tokens": ["2.0 equiv", "5 mol% catalyst", "25 °C", "3 min", "10 min", "68%", "<2%", "Table cell for 3a' is a dash"],
    },
    "scheme2_hydroboration_5c": {
        "file": f"{dataset['documents'][0]['document_id']}__scheme2_hydroboration_5c.md",
        "tokens": ["BH3", "H2O2", "THF", "60 °C", "3 h", "72%", "dr >20:1"],
    },
    "table4_entry7": {
        "file": f"{dataset['documents'][0]['document_id']}__table4_entry7.md",
        "tokens": ["2.2 mol%", "2.0 mol%", "chiral phosphoric acid 7d", "4.0 mol%", "0 °C", "3 h", "1 h", "73%", "93:7"],
    },
    "table5_entry17": {
        "file": f"{dataset['documents'][0]['document_id']}__table5_entry17_6q.md",
        "tokens": ["0.125 M", "2.2 mol%", "2.0 mol%", "4.0 mol%", "1.2 equiv", "70%", "92:8", "HPLC analysis on a chiral stationary phase"],
    },
    "eq2_conflict": {
        "file": f"{dataset['documents'][0]['document_id']}__eq2_water_competition.md",
        "tokens": ["1.1 equiv", "5 equiv", "25 °C", "3 min", "3a (8%)", "2k-O (62%)", "identity not silently merged"],
    },
}

for key, spec in checks.items():
    text = (run / "records" / spec["file"]).read_text()
    missing = [token for token in spec["tokens"] if token not in text]
    spec["missing_tokens"] = missing
    assert not missing, (key, missing)

result = {
    "document_id": dataset["documents"][0]["document_id"],
    "title": dataset["documents"][0]["title"],
    "counts": {"reactions": 91, "specific_experiments": 87, "molecules": 158},
    "coverage_pages": [1, 2, 3, 4, 5, 6],
    "all_machine_records_needs_review": True,
    "standard_matches_dataset": True,
    "markdown_checks": checks,
    "semantic_accuracy_verified": False,
}
(run / "raw_model/verification_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(result, ensure_ascii=False, indent=2))
