#!/usr/bin/env python3
"""Save deterministic validation/audit snapshots for this extraction run."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records

run = ROOT / "runs/current"
dataset = reaction_tools.read(run / "dataset.json")
standard = standard_records.read(run / "standard_records.json")
standard_errors = standard_records.validate(standard, run)
legacy_errors = reaction_tools.validate(
    dataset, run, strict_record_ids={r["record_id"] for r in dataset["experiments"]})
validation = {
    "standard_valid": not standard_errors,
    "strict_legacy_valid": not legacy_errors,
    "standard_errors": standard_errors,
    "strict_legacy_errors": legacy_errors,
    "standard_v3_present": "standard" in dataset,
    "chemical_accuracy_verified": False,
}
audit = reaction_tools.dataset_audit(dataset)
mentions = standard_records.source_mentions(dataset, ["1a", "2a", "2q", "TEMPO", "triplet"], 100)

def write(name, value):
    (run / "raw_model" / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")

write("strict_validation.json", validation)
write("audit.json", audit)
write("source_mentions_1a_2a_2q_controls.json", mentions)
print(json.dumps({"validation": validation,
                  "warning_records": len(audit["warnings"]),
                  "mention_hits": len(mentions["matches"])}, ensure_ascii=False))
