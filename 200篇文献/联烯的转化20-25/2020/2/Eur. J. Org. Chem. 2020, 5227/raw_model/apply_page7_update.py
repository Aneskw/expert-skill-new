#!/usr/bin/env python3
"""Mechanically add the page-7 evidence/context to the existing v3.1 bundle."""

from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[3]
RUN = ROOT / "runs/current"
DOC = "doc_3bed13a021d538c6"
BUNDLE = RUN / "submitted.json"
DATASET = RUN / "dataset.json"

bundle = json.loads(BUNDLE.read_text())
dataset = json.loads(DATASET.read_text())

dataset_doc = next(d for d in dataset["documents"] if d["document_id"] == DOC)
bundle_doc = next(d for d in bundle["documents"] if d["document_id"] == DOC)
bundle_evidence_ids = {e["evidence_id"] for e in bundle_doc["evidence"]}
for evidence in dataset_doc["evidence"]:
    if evidence.get("page") == 7 and evidence["evidence_id"] not in bundle_evidence_ids:
        bundle_doc["evidence"].append(evidence)

context_id = f"{DOC}:context:page7_references_and_ccdc"
context = {
    "context_id": context_id,
    "kind": "bibliography_and_crystallographic_availability",
    "evidence_ids": [f"{DOC}:page7:text", f"{DOC}:page7:image"],
    "summary": (
        "PDF page 7 contains only references [6]-[19], the CCDC deposition "
        "1981291 availability statement for compound 2o, and the received date. "
        "It contains no reaction arrow, table row, synthetic procedure, or new "
        "experiment; cited literature is not treated as an execution in this paper."
    ),
}
contexts = [c for c in bundle["source_context"] if c.get("context_id") != context_id]
contexts.append(context)
bundle["source_context"] = contexts

BUNDLE.write_text(json.dumps(bundle, ensure_ascii=False, indent=2) + "\n")
