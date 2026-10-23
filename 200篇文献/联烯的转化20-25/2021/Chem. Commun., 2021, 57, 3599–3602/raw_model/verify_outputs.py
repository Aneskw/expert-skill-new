#!/usr/bin/env python3
"""Source-specific post-merge checks for generated records and Markdown."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_eaeb2de1d7d95d0f"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_record(suffix, tokens):
    path = RUN / "records" / f"{DOC}__{suffix}.md"
    text = path.read_text()
    return {"record_id": f"{DOC}:{suffix}", "markdown": str(path.relative_to(ROOT)),
            "checked_tokens": tokens, "missing_tokens": [x for x in tokens if x not in text]}


dataset = json.loads((RUN / "dataset.json").read_text())
standard = json.loads((RUN / "standard_records.json").read_text())
coverage = json.loads((RUN / "coverage/pages_0001_0004.json").read_text())

checks = [
    check_record("table1_entry1", ["0.2", "mmol", "100", "mg", "5", "mol%", "2.0", "mL", "10", "min", "94", "90", "1,3,5-trimethoxybenzene", "Duplicate sources merged"]),
    check_record("scheme2_3ba", ["0.2", "0.22", "100", "25", "5", "60", "92", "isolated yield"]),
    check_record("scheme3_3aj", ["0.2", "0.22", "100", "25", "5", "60", "70", "3:1", "rotamer"]),
    check_record("scheme4_scaleup_3ab", ["3.0", "0.78", "3.15", "0.917", "2.36", "1.299", "79", "absolute solvent and molecular-sieve amounts at scale are unreported"]),
    check_record("scheme4_derivatization_4", ["lithium aluminium hydride", "zinc chloride", "diethyl ether", "0", "50", "yield measurement basis"]),
    check_record("scheme5_D_to_3aa_without_AgSbF6", ["without AgSbF6", "60", "yield measurement basis are unreported"]),
]

result = {
    "all_markdown_checks_pass": all(not row["missing_tokens"] for row in checks),
    "checks": checks,
    "coverage_pages": [p["page"] for p in coverage["pages"]],
    "coverage_states": [p["state"] for p in coverage["pages"]],
    "coverage_nonempty_regions": [bool(p["regions"]) for p in coverage["pages"]],
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
    "all_machine_needs_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "all_stereo_statuses_present": all(
        all(m["stereochemistry"][k]["normalized_value"] is not None
            for k in ("chirality_status", "geometry_status", "representation_status"))
        for m in standard["molecules"]),
    "dataset_standard_equals_file": dataset["standard"] == standard,
    "source_sha256": sha(RUN / dataset["documents"][0]["source_pdf"]),
    "mapped_source_sha256": dataset["documents"][0]["sha256"],
}
(RUN / "raw_model/verification_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(result, ensure_ascii=False, indent=2))
