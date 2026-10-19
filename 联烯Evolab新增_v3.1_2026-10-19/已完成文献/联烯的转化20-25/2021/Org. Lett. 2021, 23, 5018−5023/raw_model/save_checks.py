#!/usr/bin/env python3
"""Persist deterministic validation/audit results beside the raw transcription."""
import json
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools
import standard_records

base = ROOT / "runs/current"
dataset = json.loads((base / "dataset.json").read_text())
standard = json.loads((base / "standard_records.json").read_text())

audit = reaction_tools.dataset_audit(dataset)
standard_errors = standard_records.validate(standard, base)
legacy_errors = reaction_tools.validate(
    dataset, base, strict_record_ids={r["record_id"] for r in dataset["experiments"]}
)
standard_records.write(base / "raw_model/audit.json", audit)
standard_records.write(base / "raw_model/validation.json", {
    "standard_v3_valid": not standard_errors,
    "standard_v3_errors": standard_errors,
    "strict_dataset_valid": not legacy_errors,
    "strict_dataset_errors": legacy_errors,
    "chemical_accuracy_verified": False,
})

checks = {
    "all_reactions_need_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "reaction_json_count": len(list((base / "records").glob("*.json"))),
    "reaction_markdown_count": len(list((base / "records").glob("*.md"))),
    "coverage_pages": [p["page"] for p in json.loads((base / "coverage/pages_0001_0006.json").read_text())["pages"]],
    "coverage_regions_nonempty": all(p["regions"] for p in json.loads((base / "coverage/pages_0001_0006.json").read_text())["pages"]),
    "markdown_contains_table1_loadings": all(
        token in (base / f"records/{dataset['documents'][0]['document_id']}__table1_entry01.md").read_text()
        for token in ["20 mol %", "40 mol %", "0.2 mmol", "0.6 mmol", "2 mL"]
    ),
    "markdown_retains_99_vs_quantitative_note": "quantitative" in (
        base / f"records/{dataset['documents'][0]['document_id']}__scheme2_hydrogenation_to_5.md"
    ).read_text(),
    "source_pdf_hash_matches": hashlib.sha256((base / dataset["documents"][0]["source_pdf"]).read_bytes()).hexdigest()
        == dataset["documents"][0]["sha256"],
    "resource_hashes_match": all(
        hashlib.sha256((ROOT / "skills/chemical-reaction-extraction" / rel).read_bytes()).hexdigest() == expected
        for rel, expected in dataset["extraction_profile"]["resource_sha256"].items()
    ),
}
standard_records.write(base / "raw_model/integrity_checks.json", checks)
print(json.dumps(checks, ensure_ascii=False, indent=2))
