#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_dbfa3548c0baadcd"


def record_md(suffix):
    return (RUN / "records" / f"{DOC}__{suffix}.md").read_text()


checks = {
    "table1_optimized_3a": (
        record_md("table1_entry4_panisyl"),
        ["0.24 mmol", "1.2 equiv", "0.2 mmol", "0.04 mmol", "20 mol%", "3 mL", "rt", "12 h", "71%", "isolated yield", "second entry 4", "entry 5"],
    ),
    "scheme5_scaleup_3g": (
        record_md("scheme5_scaleup_3g"),
        ["1.2 equiv", "0.25 g", "1.0 mmol", "20 mol%", "1.2 equiv", "rt", "12 h", "0.28 g", "57%"],
    ),
    "scheme5_heck_to6": (
        record_md("scheme5_heck_to6"),
        ["10 mol%", "40 mol%", "2.0 equiv", "MeCN", "85 °C", "6 h", "73%"],
    ),
    "no_phosphine_control": (
        record_md("no_phosphine_product5_control"),
        ["absence of phosphine catalyst", "no product 5 could be detected", "detection method"],
    ),
}

summary = {"markdown": {}}
for name, (text, tokens) in checks.items():
    summary["markdown"][name] = {
        "checked_tokens": tokens,
        "missing_tokens": [token for token in tokens if token not in text],
    }

dataset = json.loads((RUN / "dataset.json").read_text())
standard = json.loads((RUN / "standard_records.json").read_text())
coverage = json.loads((RUN / "coverage/pages_0001_0004.json").read_text())
audit = json.loads((RUN / "audit.json").read_text())
summary["bundle"] = {
    "reactions": len(standard["reactions"]),
    "specific_experiments": sum(r["record_type"] == "specific_experiment" for r in standard["reactions"]),
    "molecules": len(standard["molecules"]),
    "source_context": len(standard["source_context"]),
    "needs_review": sum(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "dataset_standard_matches_file": dataset["standard"] == standard,
    "coverage_pages": [p["page"] for p in coverage["pages"]],
    "coverage_all_attempted": all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"]),
    "warning_records": len(audit["warnings"]),
}
summary["source_hash"] = {
    "mapped": dataset["documents"][0]["sha256"],
    "actual": hashlib.sha256((RUN / dataset["documents"][0]["source_pdf"]).read_bytes()).hexdigest(),
}
summary["all_checks_pass"] = (
    all(not row["missing_tokens"] for row in summary["markdown"].values())
    and summary["bundle"]["reactions"] == 38
    and summary["bundle"]["specific_experiments"] == 34
    and summary["bundle"]["molecules"] == 64
    and summary["bundle"]["needs_review"] == 38
    and summary["bundle"]["dataset_standard_matches_file"]
    and summary["bundle"]["coverage_pages"] == [1, 2, 3, 4]
    and summary["bundle"]["coverage_all_attempted"]
    and summary["source_hash"]["mapped"] == summary["source_hash"]["actual"]
)
(RUN / "raw_model/verification_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(summary, ensure_ascii=False, indent=2))
if not summary["all_checks_pass"]:
    raise SystemExit(1)
