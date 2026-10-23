#!/usr/bin/env python3
"""Persist deterministic validation and audit results for this bounded extraction."""
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


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


data = json.loads(DATASET.read_text())
standard_errors = standard.validate(data["standard"], RUN)
legacy_errors = legacy.validate(
    data, RUN, strict_record_ids={r["record_id"] for r in data["experiments"]}
)
audit = legacy.dataset_audit(data)
coverage = json.loads((RUN / "coverage/pages_0001_0006.json").read_text())

markdown_checks = {}
for name, snippets in {
    "preparation_1a_common_reagent_doses": [
        "{'number': 238, 'unit': 'mg'}", "{'number': 1.1, 'unit': 'mmol'}",
        "{'number': 185, 'unit': 'mg'}", "{'number': 2.2, 'unit': 'mmol'}",
        "{'number': 115, 'unit': '°C'}", "{'number': 23, 'unit': 'h'}",
    ],
    "preparation_1j_multistep_doses": [
        "{'number': 159, 'unit': 'mg'}", "{'number': 4.2, 'unit': 'mmol'}",
        "{'number': 684, 'unit': 'mg'}", "four-step isolated yield from S4",
    ],
    "preparation_1l_values_and_conflict_note": [
        "{'number': 115, 'unit': 'mg'}", "{'number': 64, 'unit': '%'}",
        "95:5", "about 32%", "no correction selected",
    ],
    "table2_wet_silica_general_conditions": [
        "{'number': 500, 'unit': 'mg'}", "1.0% (v/v) water",
        "3 × 10 mL", "petroleum ether/ethyl acetate = 10:1",
    ],
}.items():
    suffix = {
        "preparation_1a_common_reagent_doses": "preparation_1a",
        "preparation_1j_multistep_doses": "preparation_1j_s4_to_1j",
        "preparation_1l_values_and_conflict_note": "preparation_1l_s6_to_1l",
        "table2_wet_silica_general_conditions": "table2_n_substituent_entry_01",
    }[name]
    text = (RUN / "records" / f"doc_64f4a881f331038c__{suffix}.md").read_text()
    markdown_checks[name] = {
        "passed": all(snippet in text for snippet in snippets),
        "required_snippets": snippets,
    }

warning_prefixes = Counter(
    warning.split(":", 1)[0]
    for record in audit["warnings"]
    for warning in record["warnings"]
)
quality = {
    "batch_id": "allene-stereo-v31-200-20261017",
    "document_id": "doc_64f4a881f331038c",
    "bounded_pages": [1, 2, 3, 4, 5, 6],
    "standard_v3_present": "standard" in data,
    "standard_valid": not standard_errors,
    "standard_errors": standard_errors,
    "strict_compatibility_valid": not legacy_errors,
    "strict_compatibility_errors": legacy_errors,
    "chemical_accuracy_verified": False,
    "semantic_accuracy_verified": False,
    "counts": {
        "molecules": len(data["standard"]["molecules"]),
        "reactions": len(data["standard"]["reactions"]),
        "source_context": len(data["standard"]["source_context"]),
        "needs_review": sum(
            r["review"]["status"] == "needs_review"
            for r in data["standard"]["reactions"]
        ),
    },
    "coverage_contract": {
        "exact_pages_1_through_6": [p["page"] for p in coverage["pages"]]
        == [1, 2, 3, 4, 5, 6],
        "all_states_attempted": all(p["state"] == "attempted" for p in coverage["pages"]),
        "all_pages_have_regions": all(p["regions"] for p in coverage["pages"]),
    },
    "full_record_markdown_checks": markdown_checks,
    "audit_warning_prefix_counts": dict(warning_prefixes),
    "focus_warning_counts": {
        key: warning_prefixes.get(key, 0)
        for key in ("reaction_identity", "outcome_target", "measurement_basis", "structure_gap")
    },
    "documents_without_linked_si": audit["documents_without_linked_si"],
}

write(RUN / "quality_checks.json", quality)
write(RUN / "raw_model/pages_0001_0006_audit.json", audit)
print(json.dumps(quality, ensure_ascii=False, indent=2))
