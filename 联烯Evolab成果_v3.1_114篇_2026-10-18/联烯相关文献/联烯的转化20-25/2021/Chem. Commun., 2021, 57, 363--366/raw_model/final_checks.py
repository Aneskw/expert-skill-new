#!/usr/bin/env python3
"""Persist the live advisory audit and source-specific integrity checks."""

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs" / "current"
sys.path.insert(0, str(ROOT / "skills" / "chemical-reaction-extraction" / "scripts"))

import reaction_tools


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


data = json.loads((RUN / "dataset.json").read_text())
standard = json.loads((RUN / "standard_records.json").read_text())
coverage = json.loads((RUN / "coverage" / "pages_0001_0004.json").read_text())
audit = reaction_tools.dataset_audit(data)
write(RUN / "audit.json", audit)

checks = {
    "dataset_standard_matches_standard_records": data.get("standard") == standard,
    "coverage_pages_exact": [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4],
    "coverage_all_attempted": all(p["state"] == "attempted" for p in coverage["pages"]),
    "coverage_regions_nonempty": all(p["regions"] for p in coverage["pages"]),
    "unique_reaction_ids": len({r["reaction_id"] for r in standard["reactions"]}) == len(standard["reactions"]),
    "all_reactions_need_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "all_molecule_stereo_statuses_present": all(
        all(m["stereochemistry"][key]["normalized_value"] is not None
            for key in ("chirality_status", "geometry_status", "representation_status"))
        for m in standard["molecules"]
    ),
    "document_title_present": bool(standard["documents"][0].get("title")),
    "document_id_correct": standard["documents"][0]["document_id"] == "doc_1fb7cef55cdd0f91",
}

token_sets = {
    "table1_entry7": ["0.1 mmol", "2.5", "10", "1 mL", "0.5", "66", "18", "99", "chiral HPLC"],
    "scheme3_gram_scale_4aa": ["2.5", "6.3", "0.25", "0.97", "66", "97"],
    "table2_entry1": ["DABCO", "24", "1:3", "14", "31", "99", "1H NMR", "chiral HPLC"],
    "scheme5_parallel_resolution_af": ["0.1 mmol", "1 mL", "12", "45", "43", "99", "91", "chiral HPLC"],
}
markdown_checks = {}
for suffix, tokens in token_sets.items():
    path = RUN / "records" / f"doc_1fb7cef55cdd0f91__{suffix}.md"
    text = path.read_text()
    markdown_checks[suffix] = {
        "path": str(path.relative_to(ROOT)),
        "missing_tokens": [token for token in tokens if token not in text],
    }

warning_types = Counter()
for row in audit["warnings"]:
    for warning in row["warnings"]:
        warning_types[warning.split(":", 1)[0]] += 1

result = {
    "counts": {
        "molecules": len(standard["molecules"]),
        "reactions": len(standard["reactions"]),
        "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
        "source_context": len(standard["source_context"]),
        "audit_warning_records": len(audit["warnings"]),
        "audit_warning_messages_by_type": dict(warning_types),
    },
    "checks": checks,
    "markdown_checks": markdown_checks,
    "all_checks_pass": all(checks.values()) and all(not row["missing_tokens"] for row in markdown_checks.values()),
}
write(RUN / "raw_model" / "verification_summary.json", result)
print(json.dumps(result, ensure_ascii=False, indent=2))
if not result["all_checks_pass"]:
    raise SystemExit(1)
