#!/usr/bin/env python3
import collections
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_a996ed1ac1a2b8f0"

dataset = json.loads((RUN / "dataset.json").read_text())
standard = json.loads((RUN / "standard_records.json").read_text())
audit = json.loads((RUN / "audit.json").read_text())
coverage = json.loads((RUN / "coverage/pages_0001_0004.json").read_text())

type_counts = collections.Counter(r["record_type"] for r in standard["reactions"])
category_counts = collections.Counter(c for r in standard["reactions"] for c in (r["experiment_categories"]["normalized_value"] or []))
issue_counts = collections.Counter(w.split(":", 1)[0] for row in audit["warnings"] for w in row["warnings"])

checks = {
    "standard_matches_dataset": dataset["standard"] == standard,
    "document_id": standard["documents"][0]["document_id"] == DOC,
    "source_hash": hashlib.sha256((RUN / standard["documents"][0]["source_pdf"]).read_bytes()).hexdigest() == standard["documents"][0]["sha256"],
    "coverage_exact_pages": [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4],
    "coverage_regions_nonempty": all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"]),
    "all_needs_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "stereo_status_complete": all(all(m["stereochemistry"][k]["normalized_value"] is not None for k in ("chirality_status", "geometry_status", "representation_status")) for m in standard["molecules"]),
    "model_votes_empty": all(not r["classification"]["model_votes"] for r in standard["reactions"]),
    "si_not_associated": standard["documents"][0]["supplementary_material"]["status"] == "unknown",
}

markdown_expectations = {
    "table1_entry1": ["0.2 mmol", "0.4 mmol", "2.5 mL", "60 °C", "1.5 h", "81 %", "isolated yield based on 1a", "were merged"],
    "table1_entry6": ["Co(acac)3", "Co(acac)2", "0.4 mmol", "59 %", "source conflict retained"],
    "table2_entry1": ["1.5 equiv", "10 mol%", "2.5 mL", "50 °C", "3 h", "90 %", "isolated yield based on 1a", "were merged"],
    "scheme2_scale_up_3a": ["4.0 mmol", "0.69 g", "70 %", "EtOH volume is not reported"],
    "scheme4_tempo_cu": ["TEMPO", "2.0 equiv", "79 %", "measurement basis is unreported"],
}
markdown_missing = {}
for suffix, tokens in markdown_expectations.items():
    text = (RUN / "records" / f"{DOC}__{suffix}.md").read_text()
    markdown_missing[suffix] = [t for t in tokens if t not in text]
checks["representative_markdown_tokens"] = all(not v for v in markdown_missing.values())

result = {
    "checks": checks,
    "counts": {
        "molecules": len(standard["molecules"]),
        "reactions": len(standard["reactions"]),
        "record_types": dict(type_counts),
        "experiment_categories": dict(category_counts),
        "source_context": len(standard["source_context"]),
        "audit_warning_records": len(audit["warnings"]),
        "audit_issue_messages": dict(issue_counts),
    },
    "markdown_missing_tokens": markdown_missing,
    "audit_warning_record_ids": [row["record_id"] for row in audit["warnings"]],
}
assert all(checks.values()), checks
(RUN / "raw_model/verification_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(result["counts"], ensure_ascii=False))
