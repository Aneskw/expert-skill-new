#!/usr/bin/env python3
"""Save advisory audit and reproducible final consistency checks."""

import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
from reaction_tools import dataset_audit


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


dataset = read(RUN / "dataset.json")
standard = read(RUN / "standard_records.json")
coverage = read(RUN / "coverage/pages_0001_0005.json")
mapping = read(RUN / "source_mapping.json")
audit = dataset_audit(dataset)
(RUN / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")

token_checks = {
    "scope_3a": {
        "file": "records/doc_27ee161674b4e4ee__scope_3a.md",
        "tokens": ["0.2 mmol", "1.5 equiv", "1 mol%", "2 equiv", "0.5 mL", "1.5 mL", "N2", "rt", "5 W", "47 h", "72%"],
    },
    "scope_4c_footnotes": {
        "file": "records/doc_27ee161674b4e4ee__scope_4c.md",
        "tokens": ["3 equiv", "2 mol%", "1.5 mL", "0.3 mL", "H2O", "12 h", "60%"],
    },
    "figure3_5d": {
        "file": "records/doc_27ee161674b4e4ee__scope_5d.md",
        "tokens": ["0.2 mmol", "1.5 equiv", "1 mol%", "2 equiv", "36 h", "71%", "dr = 1:1"],
    },
    "figure4_6j": {
        "file": "records/doc_27ee161674b4e4ee__scope_6j.md",
        "tokens": ["0.2 mmol", "3 equiv", "2 mol%", "1 equiv", "DMF", "2 mL", "N2", "rt", "30 W", "20 h", "50%", "2.7:1"],
    },
    "lindlar": {
        "file": "records/doc_27ee161674b4e4ee__scheme2_lindlar_8.md",
        "tokens": ["5 wt %", "20 mol%", "0.5 equiv", "1 atm H2", "dry EtOAc", "rt", "15 min", "68%"],
    },
    "zirconium_deoxygenation": {
        "file": "records/doc_27ee161674b4e4ee__scheme2_deoxygenation_9.md",
        "tokens": ["1.6 equiv", "n-BuLi (2.5 M)", "3.2 equiv", "THF", "-78 °C", "1 h", "rt", "3 h", "60%", "warming duration not stated"],
    },
    "qualitative_optimization": {
        "file": "records/doc_27ee161674b4e4ee__optimization_meoh_only.md",
        "tokens": ["Methanol as sole solvent entirely impeded the transformation", "unavailable SI"],
    },
}
for item in token_checks.values():
    text = (RUN / item["file"]).read_text()
    item["missing_tokens"] = [token for token in item["tokens"] if token not in text]

evidence_ids = {e["evidence_id"] for d in dataset["documents"] for e in d["evidence"]}
coverage_evidence = {eid for page in coverage["pages"] for region in page["regions"] for eid in region["evidence_ids"]}
warnings = audit["warnings"]
warning_categories = Counter(w.split(":", 1)[0] for row in warnings for w in row["warnings"])
verification = {
    "document_id": "doc_27ee161674b4e4ee",
    "source_sha256_matches_mapping": digest(ROOT / "papers/current/source.pdf") == mapping["sha256"],
    "prepared_source_sha256_matches_dataset": digest(RUN / dataset["documents"][0]["source_pdf"]) == dataset["documents"][0]["sha256"],
    "standard_records_matches_dataset_standard": standard == dataset["standard"],
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
    "source_context_count": len(standard["source_context"]),
    "record_type_counts": dict(Counter(r["record_type"] for r in standard["reactions"])),
    "all_machine_records_need_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "all_stereo_statuses_present": all(
        m["stereochemistry"][key]["normalized_value"] is not None
        for m in standard["molecules"]
        for key in ("chirality_status", "geometry_status", "representation_status")
    ),
    "coverage_pages_exact": [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5],
    "coverage_all_attempted": all(p["state"] == "attempted" for p in coverage["pages"]),
    "coverage_region_counts": [len(p["regions"]) for p in coverage["pages"]],
    "coverage_evidence_valid": coverage_evidence <= evidence_ids,
    "coverage_record_ids_valid": all(
        rid in {r["reaction_id"] for r in standard["reactions"]}
        for p in coverage["pages"] for rid in p["record_ids"]
    ),
    "audit_warning_records": len(warnings),
    "audit_warning_categories": dict(warning_categories),
    "audit_non_measurement_warnings": [
        {"record_id": row["record_id"], "warnings": [w for w in row["warnings"] if not w.startswith("measurement_basis:")]}
        for row in warnings if any(not w.startswith("measurement_basis:") for w in row["warnings"])
    ],
    "markdown_checks": token_checks,
    "markdown_missing_tokens": [
        {"check": key, "missing": item["missing_tokens"]}
        for key, item in token_checks.items() if item["missing_tokens"]
    ],
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
(RAW := RUN / "raw_model").mkdir(exist_ok=True)
(RAW / "verification_summary.json").write_text(json.dumps(verification, ensure_ascii=False, indent=2) + "\n")

failed = [k for k, v in verification.items() if isinstance(v, bool) and not v and k not in ("semantic_accuracy_verified", "chemical_accuracy_verified")]
if verification["markdown_missing_tokens"]:
    failed.append("markdown_missing_tokens")
if verification["audit_non_measurement_warnings"]:
    failed.append("audit_non_measurement_warnings")
print(json.dumps({"failed_checks": failed, **{k: verification[k] for k in (
    "reaction_count", "molecule_count", "record_type_counts", "audit_warning_records",
    "audit_warning_categories", "coverage_region_counts")}}, ensure_ascii=False, indent=2))
raise SystemExit(1 if failed else 0)
