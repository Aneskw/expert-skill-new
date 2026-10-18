#!/usr/bin/env python3
"""Persist advisory audit and extraction-specific verification checks."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))
import reaction_tools as rt
import standard_records as st

RUN = ROOT / "runs/current"
data = st.read(RUN / "dataset.json")
audit = rt.dataset_audit(data)
st.write(RUN / "audit.json", audit)

record_md = RUN / "records/doc_b7f7182105e32edc__scheme2_scope_3a.md"
text = record_md.read_text()
required_markdown_tokens = [
    "48.8 mg", "0.2 mmol", "1 equiv", "107.36 mg", "0.44 mmol",
    "2.2 equiv", "10 mol%", "2.82 mg", "2.5 mL", "60 °C", "12 h",
    "petroleum ether/ethyl acetate = 5:1",
]

# Participant quantities were initially inherited for scope rows. Page 5 directly
# supplies the concrete 3a masses, so confirm the normalized record includes them.
standard = data["standard"]
scope3a = next(r for r in standard["reactions"] if r["reaction_id"].endswith(":scheme2_scope_3a"))
by_pid = {p["participant_id"]: p for p in scope3a["participants"]}
direct_page5_quantities = {
    "scope_1": {"mass": {"number": 48.8, "unit": "mg"}},
    "scope_2": {"mass": {"number": 107.36, "unit": "mg"}},
    "scope_bf3": {"mass": {"number": 2.82, "unit": "mg"}},
}

# This script audits only; it does not alter the merged record. The source-specific
# builder/refinement is responsible for record edits.
summary = {
    "reaction_count": len(standard["reactions"]),
    "molecule_count": len(standard["molecules"]),
    "source_context_count": len(standard["source_context"]),
    "all_reactions_needs_review": all(r["review"]["status"] == "needs_review" for r in standard["reactions"]),
    "audit_warning_record_count": len(audit["warnings"]),
    "audit_structure_gap_records": [w["record_id"] for w in audit["warnings"] if any(x.startswith("structure_gap:") for x in w["warnings"])],
    "audit_reaction_identity_records": [w["record_id"] for w in audit["warnings"] if any(x.startswith("reaction_identity:") for x in w["warnings"])],
    "audit_outcome_target_records": [w["record_id"] for w in audit["warnings"] if any(x.startswith("outcome_target:") for x in w["warnings"])],
    "audit_measurement_basis_records": [w["record_id"] for w in audit["warnings"] if any(x.startswith("measurement_basis:") for x in w["warnings"])],
    "markdown_file": str(record_md.relative_to(ROOT)),
    "markdown_missing_tokens": [token for token in required_markdown_tokens if token not in text],
    "page5_direct_quantities_expected": direct_page5_quantities,
    "page5_direct_quantities_present": {
        pid: {name: st.value(by_pid[pid]["charges"][name]) for name in fields}
        for pid, fields in direct_page5_quantities.items()
    },
    "coverage_pages": [p["page"] for p in st.read(RUN / "coverage/pages_0001_0006.json")["pages"]],
}
st.write(RUN / "raw_model/verification_summary.json", summary)
print(json.dumps(summary, ensure_ascii=False, indent=2))
