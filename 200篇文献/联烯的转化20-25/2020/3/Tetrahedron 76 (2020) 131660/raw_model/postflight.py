#!/usr/bin/env python3
"""Save deterministic postflight checks for the pages 1-6 extraction."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "skills/chemical-reaction-extraction/scripts"))
import reaction_tools as rt
import standard_records as st

RUN = ROOT / "runs/current"
data = rt.read(RUN / "dataset.json")
strict_ids = {r["record_id"] for r in data["experiments"]}
errors = rt.validate(data, RUN, strict_record_ids=strict_ids)
if "standard" in data:
    errors += st.validate(data["standard"], RUN)
validation = {"valid": not errors, "errors": sorted(set(errors)), "standard_v3_present": "standard" in data,
              "chemical_accuracy_verified": False}
audit = rt.dataset_audit(data)
coverage = rt.read(RUN / "coverage/pages_0001_0006.json")
standard = data["standard"]
conflict = next(r for r in standard["reactions"] if r["reaction_id"].endswith(":scheme10_14_to_35"))["outcomes"][0]["value"]
checks = {
    "validation_passed": validation["valid"],
    "coverage_pages_exact": [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5, 6],
    "coverage_regions_nonempty": all(p["regions"] for p in coverage["pages"]),
    "all_standard_reactions_need_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "all_legacy_records_need_review": all(r["status"] == "needs_review" for r in data["experiments"]),
    "source_reviewed_records": audit["source_reviewed_records"],
    "title_present": data["documents"][0]["title"] == "A stereotetrad-centered approach toward pironetin: Dead ends, Detour, and evolution of the synthetic strategy",
    "compound_35_conflict_preserved": conflict["data_status"] == "conflict" and len(conflict["alternatives"]) == 2,
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
}
st.write(RUN / "raw_model/strict_validation.json", validation)
st.write(RUN / "raw_model/audit.json", audit)
st.write(RUN / "raw_model/postflight.json", checks)
if errors or not all(v is True or k in {"source_reviewed_records", "reaction_count", "molecule_count"}
                     for k, v in checks.items()):
    raise SystemExit("postflight failed")
print(checks)
