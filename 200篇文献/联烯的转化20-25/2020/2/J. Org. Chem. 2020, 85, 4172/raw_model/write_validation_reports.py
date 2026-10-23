#!/usr/bin/env python3
"""Write reproducible validation and audit reports for this extraction run."""

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools as legacy
import standard_records as standard


def write(name, value):
    (RUN / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


dataset = json.loads((RUN / "dataset.json").read_text())
bundle = json.loads((RUN / "standard_records.json").read_text())

standard_errors = standard.validate(bundle, RUN)
write(
    "standard_validate.json",
    {
        "valid": not standard_errors,
        "errors": standard_errors,
        "chemical_accuracy_verified": False,
    },
)

strict_ids = {record["record_id"] for record in dataset["experiments"]}
strict_errors = legacy.validate(dataset, RUN, strict_record_ids=strict_ids)
write(
    "strict_validate.json",
    {
        "valid": not strict_errors,
        "errors": strict_errors,
        "standard_v3_present": "standard" in dataset,
        "chemical_accuracy_verified": False,
    },
)

write("audit.json", legacy.dataset_audit(dataset))
