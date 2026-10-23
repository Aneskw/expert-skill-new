#!/usr/bin/env python3
"""Persist validation, audit, coverage, and Markdown checks for pages 7--9."""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as legacy
import standard_records as standard


RUN = ROOT / "runs/current"
DATASET = RUN / "dataset.json"
DOC = "doc_64f4a881f331038c"


def dump(path: Path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


data = json.loads(DATASET.read_text())
standard_errors = standard.validate(data["standard"], RUN)
strict_errors = legacy.validate(
    data, RUN, strict_record_ids={r["record_id"] for r in data["experiments"]})
audit = legacy.dataset_audit(data)
coverage = json.loads((RUN / "coverage/pages_0007_0009.json").read_text())

markdown_specs = {
    "scale_up_all_charges_and_results": {
        "suffix": "scale_up_5_mmol_1a_2a_to_4a",
        "snippets": [
            "1.01 g", "5.0 mmol", "865 mg", "6.0 mmol", "25 mL",
            "50 mL", "1.0% (v/v) water", "12.5 g", "3 × 25 mL",
            "petroleum ether/ethyl acetate = 10:1", "1.27 g", "87%",
            "No scale-specific stereoselectivity measurement reported",
        ],
        "minimum_12h_mentions": 2,
    },
    "table4_entry5_mass_and_hplc": {
        "suffix": "table4_allenone_scope_entry_05",
        "snippets": ["40.8 mg", "Chiralpak IC", "97:3", "7.12", "8.45"],
    },
    "table4_entry9_mass_and_hplc": {
        "suffix": "table4_allenone_scope_entry_09",
        "snippets": ["50.2 mg", "Chiralpak OJ-H", "95:5", "9.28", "9.94"],
    },
}
markdown_checks = {}
for name, spec in markdown_specs.items():
    text = (RUN / "records" / f"{DOC}__{spec['suffix']}.md").read_text()
    passed = all(snippet in text for snippet in spec["snippets"])
    if "minimum_12h_mentions" in spec:
        passed = passed and text.count("12 h") >= spec["minimum_12h_mentions"]
    markdown_checks[name] = {
        "passed": passed,
        "required_snippets": spec["snippets"],
        "minimum_12h_mentions": spec.get("minimum_12h_mentions"),
    }

warning_prefixes = Counter(
    warning.split(":", 1)[0]
    for record in audit["warnings"]
    for warning in record["warnings"]
)
bounded_ids = {
    f"{DOC}:scale_up_5_mmol_1a_2a_to_4a",
    *(f"{DOC}:table4_allenone_scope_entry_{entry:02d}" for entry in range(4, 10)),
}
bounded_warning_rows = [
    row for row in audit["warnings"] if row.get("record_id") in bounded_ids
]

quality = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": DOC,
    "bounded_pages": [7, 8, 9],
    "standard_v3_present": "standard" in data,
    "standard_valid": not standard_errors,
    "standard_errors": standard_errors,
    "strict_compatibility_valid": not strict_errors,
    "strict_compatibility_errors": strict_errors,
    "chemical_accuracy_verified": False,
    "semantic_accuracy_verified": False,
    "counts": {
        "molecules": len(data["standard"]["molecules"]),
        "reactions_total": len(data["standard"]["reactions"]),
        "new_reactions_pages_7_9": 1,
        "updated_existing_reactions_pages_7_9": 6,
        "source_context_total": len(data["standard"]["source_context"]),
        "source_context_added_pages_7_9": 4,
        "needs_review": sum(
            r["review"]["status"] == "needs_review"
            for r in data["standard"]["reactions"]
        ),
    },
    "coverage_contract": {
        "exact_pages_7_through_9": [p["page"] for p in coverage["pages"]] == [7, 8, 9],
        "all_states_attempted": all(p["state"] == "attempted" for p in coverage["pages"]),
        "all_pages_have_regions": all(p["regions"] for p in coverage["pages"]),
    },
    "full_record_markdown_checks": markdown_checks,
    "source_conflicts_detected_on_pages_7_9": 0,
    "conflict_note_required": False,
    "audit_warning_prefix_counts": dict(warning_prefixes),
    "focus_warning_counts": {
        key: warning_prefixes.get(key, 0)
        for key in ("reaction_identity", "outcome_target", "measurement_basis", "structure_gap")
    },
    "bounded_record_audit_warnings": bounded_warning_rows,
    "documents_without_linked_si": audit["documents_without_linked_si"],
    "limitations": [
        "Schema and audit success do not establish chemistry correctness.",
        "No supporting-information file was supplied or explicitly associated.",
        "Scale-up stereoselectivity and several operational details are unreported; the record remains partial and needs_review.",
    ],
}

dump(RUN / "quality_checks.json", quality)
dump(RUN / "raw_model/pages_0007_0009_audit.json", audit)
print(json.dumps(quality, ensure_ascii=False, indent=2))
