#!/usr/bin/env python3
"""Save advisory audit and deterministic post-merge checks."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools
import standard_records as st

run = ROOT / "runs/current"
data = st.read(run / "dataset.json")
audit = reaction_tools.dataset_audit(data)
st.write(run / "audit.json", audit)

coverage = st.read(run / "coverage/pages_0001_0006.json")
assert [p["page"] for p in coverage["pages"]] == [1, 2, 3, 4, 5, 6]
assert all(p["state"] == "attempted" and p["regions"] for p in coverage["pages"])
assert len(data["standard"]["reactions"]) == 64
assert len(data["standard"]["molecules"]) == 99
assert all(r["review"]["status"] == "needs_review" for r in data["standard"]["reactions"])

record_md = (run / "records/doc_cb03cf42d5a5f24b__scheme2_3aa.md").read_text()
tokens = ["11.9", "124", "0.9", "3 mL", "110", "18 h", "DCM", "H2O",
          "Na2SO4", "petroleum ether", "ethyl acetate", "dichloromethane",
          "Source conflict", "0.2/0.4", "0.3/0.6"]
summary = st.read(run / "raw_model/verification_summary.json")
summary["markdown_checked_record"] = "doc_cb03cf42d5a5f24b:scheme2_3aa"
summary["markdown_missing_tokens"] = [token for token in tokens if token not in record_md]
summary["coverage_pages"] = [p["page"] for p in coverage["pages"]]
summary["coverage_regions_per_page"] = [len(p["regions"]) for p in coverage["pages"]]
summary["needs_review_reactions"] = sum(r["review"]["status"] == "needs_review" for r in data["standard"]["reactions"])
summary["audit_warning_records"] = len(audit["warnings"])
summary["audit_warnings"] = audit["warnings"]
st.write(run / "raw_model/verification_summary.json", summary)
print(summary)
