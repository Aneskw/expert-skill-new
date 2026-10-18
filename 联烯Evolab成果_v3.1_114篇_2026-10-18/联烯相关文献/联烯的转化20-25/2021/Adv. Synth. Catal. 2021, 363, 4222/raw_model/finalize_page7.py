#!/usr/bin/env python3
"""Run deterministic post-merge checks for the page-7 extraction increment."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills/chemical-reaction-extraction/scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools
import standard_records as st


run = ROOT / "runs/current"
data = st.read(run / "dataset.json")
standard = data["standard"]
coverage = st.read(run / "coverage/pages_0007_0007.json")

assert coverage["document_id"] == "doc_cb03cf42d5a5f24b"
assert coverage["semantic_accuracy_verified"] is False
assert [page["page"] for page in coverage["pages"]] == [7]
assert all(page["state"] == "attempted" and page["regions"] for page in coverage["pages"])
assert coverage["pages"][0]["record_ids"] == []

assert len(standard["reactions"]) == 64
assert len(standard["molecules"]) == 99
assert len(standard["source_context"]) == 5
assert all(reaction["review"]["status"] == "needs_review" for reaction in standard["reactions"])

status_fields = ("chirality_status", "geometry_status", "representation_status")
for molecule in standard["molecules"]:
    for field in status_fields:
        fact = molecule["stereochemistry"][field]
        assert fact["normalized_value"] is not None

page7_contexts = [
    context for context in standard["source_context"]
    if "doc_cb03cf42d5a5f24b:page7:image" in context.get("evidence_ids", [])
]
assert {context["context_id"] for context in page7_contexts} == {
    "doc_cb03cf42d5a5f24b:context_page7_references",
}

audit = reaction_tools.dataset_audit(data)
st.write(run / "audit.json", audit)

warning_text = [
    warning
    for item in audit["warnings"]
    for warning in item.get("warnings", [])
]
summary = {
    "document_id": "doc_cb03cf42d5a5f24b",
    "pages_checked": [7],
    "page7_regions": len(coverage["pages"][0]["regions"]),
    "page7_reaction_ids": coverage["pages"][0]["record_ids"],
    "new_reactions": 0,
    "new_molecules": 0,
    "total_reactions": len(standard["reactions"]),
    "total_molecules": len(standard["molecules"]),
    "source_context_records": len(standard["source_context"]),
    "page7_context_ids": [context["context_id"] for context in page7_contexts],
    "needs_review_reactions": sum(
        reaction["review"]["status"] == "needs_review"
        for reaction in standard["reactions"]
    ),
    "molecules_with_required_stereo_statuses": len(standard["molecules"]),
    "audit_warning_records": len(audit["warnings"]),
    "reaction_identity_warnings": [text for text in warning_text if "reaction_identity" in text],
    "outcome_target_warnings": [text for text in warning_text if "outcome_target" in text],
    "measurement_basis_warnings": [text for text in warning_text if "measurement_basis" in text],
    "structure_gap_warnings": [text for text in warning_text if "structure_gap" in text],
    "semantic_accuracy_verified": False,
    "chemical_accuracy_verified": False,
}
st.write(run / "raw_model/verification_summary_pages_0007_0007.json", summary)
print(summary)
