#!/usr/bin/env python3
"""Synchronize source-derived title metadata and persist the final advisory audit."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "skills" / "chemical-reaction-extraction" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import reaction_tools as rt
import standard_records as st

TITLE = "TfOH-Catalyzed Phosphinylation of 2,3-Allenols into γ-Ketophosphine Oxides"
DOC = "doc_4e4afb29c8afe53f"


def set_title(documents):
    for document in documents:
        if document.get("document_id") == DOC:
            document["title"] = TITLE


current = ROOT / "runs" / "current"
dataset_path = current / "dataset.json"
dataset = st.read(dataset_path)
set_title(dataset["documents"])
set_title(dataset["standard"]["documents"])
st.write(dataset_path, dataset)

for path in [
    current / "submitted.json",
    current / "standard_records.json",
    current / "raw_model" / "extracted_bundle_pages_0001_0006.json",
]:
    bundle = st.read(path)
    set_title(bundle["documents"])
    st.write(path, bundle)

# Re-render only from the synchronized standard bundle. These are generated
# views; source records and evidence remain untouched.
standard = st.read(current / "standard_records.json")
st.render(standard, current)
rt.report(current)

audit = rt.dataset_audit(st.read(dataset_path))
st.write(current / "raw_model" / "audit_pages_0001_0006.json", audit)

print({
    "title": TITLE,
    "warning_records": len(audit["warnings"]),
    "source_reviewed_records": audit["source_reviewed_records"],
})
