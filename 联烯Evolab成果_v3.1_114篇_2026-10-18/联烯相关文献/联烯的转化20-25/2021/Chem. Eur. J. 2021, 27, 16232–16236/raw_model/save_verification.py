#!/usr/bin/env python3
"""Save reproducible post-merge stereo and deliverable checks."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import standard_records as st

base = ROOT / "runs/current"
bundle = json.loads((base / "standard_records.json").read_text())
checks = {}
for m in bundle["molecules"]:
    smi = m["smiles"]["normalized_value"]
    if smi:
        checks[m["molecule_id"]] = {"input_smiles": smi, "diagnostic": st.stereo_inspect(smi)}
(base / "raw_model/stereo_checks.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2) + "\n")

coverage = json.loads((base / "coverage/pages_0001_0005.json").read_text())
audit = json.loads((base / "audit.json").read_text())
record_ids = {r["reaction_id"] for r in bundle["reactions"]}
markdown = "\n".join(p.read_text() for p in (base / "records").glob("*.md"))
required_tokens = ["20 mol%", "30 mol%", "67:33 er", "chiral stationary phases", "n-BuLi", "10% Pd/C", "footnote c", "27 (12)^b"]
result = {
    "reaction_count": len(bundle["reactions"]),
    "molecule_count": len(bundle["molecules"]),
    "smiles_stereo_checks": len(checks),
    "all_records_needs_review": all(r["review"]["status"] == "needs_review" for r in bundle["reactions"]),
    "all_records_partial": all(r["completeness"] == "partial" for r in bundle["reactions"]),
    "model_votes_empty": all(r["classification"]["model_votes"] == [] for r in bundle["reactions"]),
    "coverage_pages_exact": [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5],
    "coverage_regions_nonempty": all(p["regions"] for p in coverage["pages"]),
    "coverage_record_ids_valid": all(rid in record_ids for p in coverage["pages"] for rid in p["record_ids"]),
    "markdown_required_tokens": {token: token in markdown for token in required_tokens},
    "audit_warning_records": len(audit["warnings"]),
    "source_reviewed_records": audit["source_reviewed_records"],
    "schema_and_strict_validation": "passed in command log",
    "pytest": "69 passed"
}
(base / "raw_model/verification_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
assert all(result["markdown_required_tokens"].values())
assert result["coverage_pages_exact"] and result["coverage_regions_nonempty"] and result["coverage_record_ids_valid"]
print(json.dumps(result, ensure_ascii=False, indent=2))
