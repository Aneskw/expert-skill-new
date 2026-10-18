#!/usr/bin/env python3
"""Persist the live audit and source-specific post-merge checks."""

from collections import Counter
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt

RUN = ROOT / "runs/current"
data = json.loads((RUN / "dataset.json").read_text())
standard = json.loads((RUN / "standard_records.json").read_text())
audit = rt.dataset_audit(data)
rt.write(RUN / "audit.json", audit)

checks = {
    "doc_923e87cf958fa06a__table1_entry11.md": [
        "0.2", "0.24", "0.3", "0.4", "5.0", "2", "20", "95", "11.3:1", "cross-page duplicate merged"
    ],
    "doc_923e87cf958fa06a__scheme3_5u.md": [
        "0.2", "0.24", "0.3", "5.0", "12", "99", "8.5:1", "X-ray"
    ],
    "doc_923e87cf958fa06a__scheme4_methylation_8.md": [
        "methyl iodide", "silver(I) oxide", "dichloromethane", "r.t.", "48", "85", ">20:1"
    ],
}
markdown_checks = {}
for filename, tokens in checks.items():
    text = (RUN / "records" / filename).read_text()
    markdown_checks[filename] = {
        "checked_tokens": tokens,
        "missing_tokens": [token for token in tokens if token not in text],
    }

coverage = json.loads((RUN / "coverage/pages_0001_0004.json").read_text())
warning_categories = Counter()
for row in audit["warnings"]:
    warning_categories.update(message.split(":", 1)[0] for message in row["warnings"])

summary = {
    "standard_matches_dataset": standard == data.get("standard"),
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
    "source_context_count": len(standard["source_context"]),
    "all_reactions_need_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "source_reviewed_records": audit["source_reviewed_records"],
    "coverage_pages": [p["page"] for p in coverage["pages"]],
    "all_coverage_regions_nonempty": all(p["regions"] for p in coverage["pages"]),
    "warning_bearing_records": len(audit["warnings"]),
    "warning_categories": dict(warning_categories),
    "markdown_checks": markdown_checks,
}
assert summary["standard_matches_dataset"]
assert summary["reaction_count"] == 64
assert summary["molecule_count"] == 115
assert summary["all_reactions_need_review"]
assert summary["source_reviewed_records"] == 0
assert summary["coverage_pages"] == [1, 2, 3, 4]
assert summary["all_coverage_regions_nonempty"]
assert all(not row["missing_tokens"] for row in markdown_checks.values())
rt.write(RUN / "raw_model/verification_summary.json", summary)
print(json.dumps(summary, ensure_ascii=False, indent=2))
